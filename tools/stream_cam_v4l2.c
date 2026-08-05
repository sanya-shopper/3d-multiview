/* Linux (V4L2) camera streamer for the multiview live protocol -- the
 * Ubuntu twin of tools/stream_cam.swift: webcam -> grayscale -> TCP to
 * livehub.  Only the capture shim differs from the Swift tool; the wire
 * format and pacing behaviour are identical.
 *
 * PROTOCOL (little-endian, one message per frame, back to back):
 *   "MVFR" | u32 camid | u32 w | u32 h | u64 seq | f64 t_mono | w*h gray
 * t_mono is CLOCK_MONOTONIC seconds of THIS machine as an IEEE-754
 * double (each camera is its own clock domain; the displayed counter is
 * the shared clock, per the paper's temporal model).  seq counts frames
 * actually sent; a frame is either shipped complete or not at all --
 * the byte stream is always parseable at 32+w*h boundaries.
 *
 * UBUNTU RUNBOOK:
 *   packages : build-essential only (V4L2 headers ship with libc-dev);
 *              v4l-utils is optional but handy:
 *                  v4l2-ctl --list-devices
 *                  v4l2-ctl -d /dev/video0 --list-formats-ext
 *   perms    : the user must be in the `video` group
 *              (sudo usermod -aG video $USER; re-login) or run via sudo.
 *   build    : make stream_cam_v4l2
 *   run      : ./stream_cam_v4l2 <hub-host> <port> <camid> [fps] [device]
 *              fps defaults to 5, device to /dev/video0.
 *   The tool requests YUYV 1280x720 and accepts whatever size the
 *   driver falls back to, but NOT another pixel format: if the device
 *   only does MJPEG at the chosen size it exits with a clear message
 *   (JPEG decoding is out of scope) -- pick a YUYV-capable mode.
 *   No auto-reconnect: on a lost hub connection it exits and the
 *   runbook restarts it.
 *
 * OWNERSHIP (parallel build): this file belongs to the v4l2 work item.
 */
#ifdef __linux__

#define _POSIX_C_SOURCE 200809L

#include <errno.h>
#include <fcntl.h>
#include <netdb.h>
#include <stdarg.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <poll.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/ioctl.h>
#include <sys/mman.h>
#include <sys/select.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#include <linux/videodev2.h>

/* t_mono is memcpy'd through an unsigned long long; both must be 8
 * bytes for the wire format to hold (file-scope so no unused warning) */
typedef char mv_assert_double8[sizeof(double) == 8 ? 1 : -1];
typedef char mv_assert_ull8[sizeof(unsigned long long) == 8 ? 1 : -1];

#define NBUFS 4              /* mmap streaming buffers requested */
#define SEND_BUDGET_MS 500   /* blocking window for one frame's send */
#define SEND_STALL_MS 10000  /* mid-frame hard cap before giving up */
#define IDLE_TIMEOUT_S 2     /* select() slice on the capture fd */
#define IDLE_MAX 5           /* consecutive empty slices before exit */

static volatile sig_atomic_t g_stop = 0;

static void on_sig(int sig)
{
    (void)sig;
    g_stop = 1;
}

/* the classic xioctl: retry the ioctl for as long as it is interrupted
 * by a signal (unless that signal was our own shutdown request) */
static int xioctl(int fd, unsigned long req, void *arg)
{
    int r;
    do {
        r = ioctl(fd, req, arg);
    } while (r == -1 && errno == EINTR && !g_stop);
    return r;
}

static double now_mono(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

/* Spoken feedback, matching the Swift twin: during the live session
 * this laptop's screen faces the volume, so the state changes are
 * announced aloud (speech-dispatcher, falling back to espeak; silent
 * if neither is installed). Text is program-authored only. MV_MUTE=1
 * disables. */
static void say(const char *fmt, ...)
{
    char msg[200], cmd[480];
    int rc;
    va_list ap;
    if (getenv("MV_MUTE"))
        return;
    va_start(ap, fmt);
    vsnprintf(msg, sizeof(msg), fmt, ap);
    va_end(ap);
    snprintf(cmd, sizeof(cmd),
             "{ spd-say \"%s\" 2>/dev/null || espeak \"%s\"; } "
             ">/dev/null 2>&1 &", msg, msg);
    rc = system(cmd);
    (void)rc;
}

static void put32(unsigned char *p, unsigned v)
{
    p[0] = (unsigned char)(v & 255);
    p[1] = (unsigned char)((v >> 8) & 255);
    p[2] = (unsigned char)((v >> 16) & 255);
    p[3] = (unsigned char)((v >> 24) & 255);
}

static void put64(unsigned char *p, unsigned long long v)
{
    int i;
    for (i = 0; i < 8; i++)
        p[i] = (unsigned char)((v >> (8 * i)) & 255);
}

static void fourcc_str(unsigned v, char out[5])
{
    int i;
    for (i = 0; i < 4; i++) {
        int c = (int)((v >> (8 * i)) & 255);
        out[i] = (c >= 32 && c < 127) ? (char)c : '?';
    }
    out[4] = 0;
}

/* ---- capture device state (globals so signal-path teardown is easy) */
static int g_vfd = -1;
static void *g_map[NBUFS];
static size_t g_maplen[NBUFS];
static int g_nmapped = 0;
static int g_streaming = 0;

static void teardown(void)
{
    int i;
    if (g_vfd >= 0) {
        if (g_streaming) {
            enum v4l2_buf_type type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
            xioctl(g_vfd, VIDIOC_STREAMOFF, &type);
            g_streaming = 0;
        }
        for (i = 0; i < g_nmapped; i++)
            munmap(g_map[i], g_maplen[i]);
        g_nmapped = 0;
        close(g_vfd);
        g_vfd = -1;
    }
}

/* ---- TCP client -------------------------------------------------- */
static int connect_hub(const char *host, const char *port)
{
    struct addrinfo hints, *res;
    int s, one = 1, fl;

    memset(&hints, 0, sizeof(hints));
    hints.ai_family = AF_INET;
    hints.ai_socktype = SOCK_STREAM;
    if (getaddrinfo(host, port, &hints, &res) != 0) {
        fprintf(stderr, "cannot resolve %s\n", host);
        return -1;
    }
    s = socket(res->ai_family, res->ai_socktype, res->ai_protocol);
    if (s < 0 || connect(s, res->ai_addr, res->ai_addrlen) != 0) {
        fprintf(stderr, "cannot connect to %s:%s\n", host, port);
        if (s >= 0)
            close(s);
        freeaddrinfo(res);
        return -1;
    }
    freeaddrinfo(res);
    setsockopt(s, IPPROTO_TCP, TCP_NODELAY, &one, sizeof(one));
    /* non-blocking AFTER connect: sender-side liveness is handled by
     * the poll() logic in send_frame, never by blocking in send() */
    fl = fcntl(s, F_GETFL, 0);
    if (fl < 0 || fcntl(s, F_SETFL, fl | O_NONBLOCK) < 0) {
        fprintf(stderr, "cannot set socket non-blocking\n");
        close(s);
        return -1;
    }
    return s;
}

/* Ship one complete frame or none at all -- the TORN-FRAME INVARIANT:
 * a receiver must never see a header without its full payload.
 *   - If a previous frame left the socket backlogged, probe with a
 *     zero-timeout poll(): still unwritable -> drop this frame whole
 *     (return 1) before a single byte goes out.
 *   - Otherwise send with a SEND_BUDGET_MS blocking window (100 ms
 *     poll slices).  If the budget expires with off == 0, nothing was
 *     sent: drop the frame, mark backlogged.  Once off > 0 we are
 *     committed and keep polling up to SEND_STALL_MS; only a dead or
 *     fully stalled hub ends the process (return -1), which closes the
 *     socket -- still a parseable stream end, never a torn frame.
 * Returns 0 sent, 1 dropped (nothing sent), -1 fatal. */
static int send_frame(int s, const unsigned char *buf, size_t len,
                      int *backlogged)
{
    size_t off = 0;
    double t0 = now_mono();

    if (*backlogged) {
        struct pollfd p;
        int r;
        p.fd = s;
        p.events = POLLOUT;
        p.revents = 0;
        r = poll(&p, 1, 0);
        if (r < 0 && errno != EINTR)
            return -1;
        if (r > 0 && (p.revents & (POLLERR | POLLHUP | POLLNVAL)))
            return -1; /* dead socket must not drop-loop forever */
        if (r <= 0 || !(p.revents & POLLOUT))
            return 1; /* socket not drained yet: drop this frame */
        *backlogged = 0;
    }
    while (off < len) {
        ssize_t k = send(s, buf + off, len - off, MSG_NOSIGNAL);
        if (k > 0) {
            off += (size_t)k;
            continue;
        }
        if (k < 0 && errno == EINTR)
            continue;
        if (k < 0 && (errno == EAGAIN || errno == EWOULDBLOCK)) {
            struct pollfd p;
            int r;
            double waited_ms = (now_mono() - t0) * 1000.0;
            if (off == 0 && waited_ms >= SEND_BUDGET_MS) {
                *backlogged = 1;
                return 1; /* whole frame dropped, stream unharmed */
            }
            if (waited_ms >= SEND_STALL_MS)
                return -1; /* hub wedged mid-frame: give up cleanly */
            p.fd = s;
            p.events = POLLOUT;
            p.revents = 0;
            r = poll(&p, 1, 100);
            if (r < 0 && errno != EINTR)
                return -1;
            if (r > 0 && (p.revents & (POLLERR | POLLHUP | POLLNVAL)))
                return -1;
            continue;
        }
        return -1; /* EPIPE, ECONNRESET, k == 0, ... */
    }
    if ((now_mono() - t0) * 1000.0 >= SEND_BUDGET_MS)
        *backlogged = 1; /* frame completed late: drain before more */
    return 0;
}

/* ---- hub -> camera channel (clock sync) --------------------------
 * The hub piggybacks 12-byte messages on this same TCP connection:
 *   "MVAK" | u32 rx | u32 decoded    liveness ack (ignored here)
 *   "MVPB" | f64 t_hub_send          clock-sync probe (~1 Hz)
 * A probe is answered immediately with
 *   "MVTS" | u32 camid | f64 t_hub_echo | f64 t_cam
 * where t_hub_echo is the probe's payload copied BIT-EXACTLY (raw
 * bytes, never through a double) and t_cam is CLOCK_MONOTONIC at
 * reply time -- the same clock the MVFR headers carry, which is what
 * lets the hub map capture timestamps onto its own clock.
 *
 * Incoming bytes can arrive split at ANY boundary (the mirror image
 * of the torn-frame invariant above), so a small carry buffer keeps
 * the partial tail (< 12 bytes) between drain calls; unrecognized
 * bytes are skipped one at a time until the next magic. */
static unsigned char g_rxbuf[256];
static size_t g_rxlen;

/* the reply goes out whole or not at all on this non-blocking socket:
 * if the FIRST byte cannot be written the probe is skipped (nothing
 * torn, the next probe is a second away); once committed, short polls
 * finish the message -- only a hub dead for >2 s aborts mid-message,
 * and that connection is done for anyway */
static void send_mvts(int s, const unsigned char *p, size_t len)
{
    size_t off = 0;
    int polls = 0;
    while (off < len) {
        ssize_t k = send(s, p + off, len - off, MSG_NOSIGNAL);
        if (k > 0) {
            off += (size_t)k;
            continue;
        }
        if (k < 0 && errno == EINTR)
            continue;
        if (k < 0 && (errno == EAGAIN || errno == EWOULDBLOCK)) {
            struct pollfd pf;
            if (off == 0)
                return; /* backlogged: skip this probe entirely */
            if (++polls > 20)
                return; /* hub wedged mid-reply: give up */
            pf.fd = s;
            pf.events = POLLOUT;
            pf.revents = 0;
            (void)poll(&pf, 1, 100);
            continue;
        }
        return; /* dead socket: the frame path will notice */
    }
}

static void drain_hub(int s, unsigned camid)
{
    size_t i = 0;

    if (g_rxlen < sizeof(g_rxbuf)) {
        ssize_t k = recv(s, g_rxbuf + g_rxlen,
                         sizeof(g_rxbuf) - g_rxlen, MSG_DONTWAIT);
        if (k > 0)
            g_rxlen += (size_t)k;
    }
    while (i + 12 <= g_rxlen) {
        if (memcmp(g_rxbuf + i, "MVAK", 4) == 0) {
            i += 12;
            continue;
        }
        if (memcmp(g_rxbuf + i, "MVPB", 4) == 0) {
            unsigned char rep[24];
            unsigned long long bits;
            double t = now_mono();
            rep[0] = 'M'; rep[1] = 'V'; rep[2] = 'T'; rep[3] = 'S';
            put32(rep + 4, camid);
            memcpy(rep + 8, g_rxbuf + i + 4, 8); /* bit-exact echo */
            memcpy(&bits, &t, 8);
            put64(rep + 16, bits);
            send_mvts(s, rep, sizeof(rep));
            i += 12;
            continue;
        }
        i++; /* resync byte by byte */
    }
    memmove(g_rxbuf, g_rxbuf + i, g_rxlen - i);
    g_rxlen -= i;
}

int main(int argc, char **argv)
{
    const char *device = "/dev/video0";
    double fps = 5.0, min_gap, last_sent = 0.0;
    unsigned camid, w, h, bpl;
    int sock, i, idle = 0, dqerr = 0, backlogged = 0;
    unsigned long long seq = 0, dropped = 0, skipped_short = 0;
    unsigned char *msg = NULL;
    size_t msglen;
    struct sigaction sa;

    if (argc < 4 || argc > 6) {
        fprintf(stderr,
                "usage: %s <hub-host> <port> <camid> [fps] [device]\n"
                "       fps defaults to 5, device to /dev/video0\n",
                argv[0]);
        return 1;
    }
    camid = (unsigned)atoi(argv[3]);
    if (argc >= 5) {
        fps = atof(argv[4]);
        if (fps <= 0.0)
            fps = 5.0;
    }
    if (argc >= 6)
        device = argv[5];
    min_gap = 1.0 / fps;

    memset(&sa, 0, sizeof(sa));
    sa.sa_handler = on_sig;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = 0; /* no SA_RESTART: blocking calls must see EINTR */
    sigaction(SIGINT, &sa, NULL);
    sigaction(SIGTERM, &sa, NULL);
    signal(SIGPIPE, SIG_IGN);

    sock = connect_hub(argv[1], argv[2]);
    if (sock < 0) {
        say("camera %u cannot reach the hub", camid);
        return 1;
    }
    printf("connected to %s:%s as camera %u, %.1f fps\n", argv[1],
           argv[2], camid, fps);

    /* ---- open + configure the capture device (blocking fd; the read
     * loop paces itself with select(), so DQBUF never spins) -------- */
    g_vfd = open(device, O_RDWR);
    if (g_vfd < 0) {
        fprintf(stderr, "cannot open %s: %s (in the `video` group? "
                        "try v4l2-ctl --list-devices)\n",
                device, strerror(errno));
        return 1;
    }
    {
        struct v4l2_capability cap;
        memset(&cap, 0, sizeof(cap));
        if (xioctl(g_vfd, VIDIOC_QUERYCAP, &cap) < 0) {
            fprintf(stderr, "%s: VIDIOC_QUERYCAP failed: %s (not a "
                            "V4L2 device?)\n", device, strerror(errno));
            teardown();
            return 1;
        }
        if (!(cap.capabilities & V4L2_CAP_VIDEO_CAPTURE)
            || !(cap.capabilities & V4L2_CAP_STREAMING)) {
            fprintf(stderr, "%s (%s) lacks capture+streaming caps "
                            "(0x%08x)\n", device, (const char *)cap.card,
                    (unsigned)cap.capabilities);
            teardown();
            return 1;
        }
        printf("device %s: %s\n", device, (const char *)cap.card);
    }
    {
        struct v4l2_format fmt;
        char fcc[5];
        memset(&fmt, 0, sizeof(fmt));
        fmt.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
        fmt.fmt.pix.width = 1280;
        fmt.fmt.pix.height = 720;
        fmt.fmt.pix.pixelformat = V4L2_PIX_FMT_YUYV;
        fmt.fmt.pix.field = V4L2_FIELD_NONE;
        if (xioctl(g_vfd, VIDIOC_S_FMT, &fmt) < 0) {
            fprintf(stderr, "%s: VIDIOC_S_FMT failed: %s\n", device,
                    strerror(errno));
            teardown();
            return 1;
        }
        /* the driver writes the format it actually granted back into
         * fmt: size fallback is fine, a foreign pixel format is not */
        if (fmt.fmt.pix.pixelformat != V4L2_PIX_FMT_YUYV) {
            fourcc_str(fmt.fmt.pix.pixelformat, fcc);
            fprintf(stderr, "%s delivers fourcc '%s' (%ux%u), not "
                            "YUYV; JPEG/compressed decoding is out of "
                            "scope for this streamer.  Pick a "
                            "YUYV-capable mode (v4l2-ctl -d %s "
                            "--list-formats-ext)\n", device, fcc,
                    (unsigned)fmt.fmt.pix.width,
                    (unsigned)fmt.fmt.pix.height, device);
            teardown();
            return 1;
        }
        w = fmt.fmt.pix.width;
        h = fmt.fmt.pix.height;
        bpl = fmt.fmt.pix.bytesperline;
        if (bpl == 0) /* rare driver shortcut: rows are packed tight */
            bpl = 2 * w;
        if (w < 64 || h < 64 || w > 4096 || h > 4096 || bpl < 2 * w) {
            fprintf(stderr, "%s: unusable geometry %ux%u stride %u\n",
                    device, w, h, bpl);
            teardown();
            return 1;
        }
        if (w != 1280 || h != 720)
            printf("driver fell back to %ux%u (asked 1280x720)\n", w,
                   h);
    }
    {
        struct v4l2_requestbuffers req;
        memset(&req, 0, sizeof(req));
        req.count = NBUFS;
        req.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
        req.memory = V4L2_MEMORY_MMAP;
        if (xioctl(g_vfd, VIDIOC_REQBUFS, &req) < 0) {
            fprintf(stderr, "%s: VIDIOC_REQBUFS(mmap) failed: %s\n",
                    device, strerror(errno));
            teardown();
            return 1;
        }
        if (req.count < 2) {
            fprintf(stderr, "%s: driver granted only %u buffers\n",
                    device, (unsigned)req.count);
            teardown();
            return 1;
        }
        for (i = 0; i < (int)req.count && i < NBUFS; i++) {
            struct v4l2_buffer buf;
            memset(&buf, 0, sizeof(buf));
            buf.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
            buf.memory = V4L2_MEMORY_MMAP;
            buf.index = (unsigned)i;
            if (xioctl(g_vfd, VIDIOC_QUERYBUF, &buf) < 0) {
                fprintf(stderr, "%s: VIDIOC_QUERYBUF %d failed: %s\n",
                        device, i, strerror(errno));
                teardown();
                return 1;
            }
            g_map[i] = mmap(NULL, buf.length, PROT_READ | PROT_WRITE,
                            MAP_SHARED, g_vfd, (off_t)buf.m.offset);
            if (g_map[i] == MAP_FAILED) {
                fprintf(stderr, "%s: mmap buffer %d failed: %s\n",
                        device, i, strerror(errno));
                teardown();
                return 1;
            }
            g_maplen[i] = buf.length;
            g_nmapped = i + 1;
        }
        for (i = 0; i < g_nmapped; i++) {
            struct v4l2_buffer buf;
            memset(&buf, 0, sizeof(buf));
            buf.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
            buf.memory = V4L2_MEMORY_MMAP;
            buf.index = (unsigned)i;
            if (xioctl(g_vfd, VIDIOC_QBUF, &buf) < 0) {
                fprintf(stderr, "%s: initial VIDIOC_QBUF %d failed: "
                                "%s\n", device, i, strerror(errno));
                teardown();
                return 1;
            }
        }
    }
    {
        enum v4l2_buf_type type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
        if (xioctl(g_vfd, VIDIOC_STREAMON, &type) < 0) {
            fprintf(stderr, "%s: VIDIOC_STREAMON failed: %s\n", device,
                    strerror(errno));
            teardown();
            return 1;
        }
        g_streaming = 1;
    }

    msglen = 32 + (size_t)w * h;
    msg = malloc(msglen);
    if (!msg) {
        fprintf(stderr, "out of memory for %zu-byte frame\n", msglen);
        teardown();
        return 1;
    }
    printf("streaming from %s (%ux%u YUYV, stride %u); Ctrl-C to "
           "stop\n", device, w, h, bpl);
    say("camera %u streaming to the hub", camid);

    /* ---- capture / pace / convert / send loop --------------------- */
    while (!g_stop) {
        fd_set fds;
        struct timeval tv;
        struct v4l2_buffer buf;
        int r;
        double now;

        /* answer any pending clock-sync probes once per iteration;
         * non-blocking either way, so the capture cadence is kept */
        drain_hub(sock, camid);

        FD_ZERO(&fds);
        FD_SET(g_vfd, &fds);
        tv.tv_sec = IDLE_TIMEOUT_S;
        tv.tv_usec = 0;
        r = select(g_vfd + 1, &fds, NULL, NULL, &tv);
        if (r < 0) {
            if (errno == EINTR)
                continue; /* g_stop check at loop top */
            fprintf(stderr, "select on %s failed: %s\n", device,
                    strerror(errno));
            break;
        }
        if (r == 0) {
            idle++;
            if (idle >= IDLE_MAX) {
                fprintf(stderr, "no frames from %s for %d s; giving "
                                "up\n", device, idle * IDLE_TIMEOUT_S);
                say("camera %u gets no frames from its device",
                    camid);
                break;
            }
            fprintf(stderr, "waiting for frames from %s...\n", device);
            continue;
        }
        idle = 0;

        memset(&buf, 0, sizeof(buf));
        buf.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
        buf.memory = V4L2_MEMORY_MMAP;
        if (xioctl(g_vfd, VIDIOC_DQBUF, &buf) < 0) {
            int err = errno; /* fprintf below may clobber errno */
            if (g_stop)
                break;
            if (err == EAGAIN)
                continue; /* raced another wakeup; select again */
            dqerr++;
            fprintf(stderr, "VIDIOC_DQBUF: %s\n", strerror(err));
            if (err == EIO && dqerr < 10)
                continue; /* transient (e.g. USB glitch): retry */
            break;
        }
        dqerr = 0;
        if (buf.index >= (unsigned)g_nmapped) {
            fprintf(stderr, "driver returned bogus buffer index %u\n",
                    (unsigned)buf.index);
            break;
        }

        /* pace by dropping, exactly like the Swift twin's minGap: a
         * too-early frame goes straight back to the driver.  A short
         * read (driver delivered fewer bytes than the last row needs,
         * i.e. under (h-1)*bpl + 2*w) is skipped the same way. */
        now = now_mono();
        if (now - last_sent < min_gap
            || buf.bytesused < (size_t)(h - 1) * bpl + 2 * w) {
            struct v4l2_buffer qb;
            if (now - last_sent >= min_gap)
                skipped_short++;
            memset(&qb, 0, sizeof(qb));
            qb.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
            qb.memory = V4L2_MEMORY_MMAP;
            qb.index = buf.index;
            if (xioctl(g_vfd, VIDIOC_QBUF, &qb) < 0) {
                fprintf(stderr, "VIDIOC_QBUF: %s\n", strerror(errno));
                break;
            }
            continue;
        }
        last_sent = now;

        /* YUYV 4:2:2 packed: each pixel pair is [Y0 U Y1 V], so every
         * pixel owns 2 bytes and its luma is the FIRST of them: the Y
         * of pixel x in row r lives at r*bytesperline + 2*x.
         * Worked stride example: w=1280 means 2*1280=2560 data bytes
         * per row; a driver padding rows to bytesperline=2592 puts the
         * Y of pixel (r=1, x=3) at 1*2592 + 6 = 2598, not 1*2560 + 6 =
         * 2566 -- hence rows step by bpl while the output steps by w. */
        {
            const unsigned char *base = (const unsigned char *)
                                        g_map[buf.index];
            unsigned char *dst = msg + 32;
            unsigned rr, x;
            for (rr = 0; rr < h; rr++) {
                const unsigned char *src = base + (size_t)rr * bpl;
                unsigned char *drow = dst + (size_t)rr * w;
                for (x = 0; x < w; x++)
                    drow[x] = src[2 * x];
            }
        }
        /* frame consumed from the mmap'd buffer: hand it back to the
         * driver BEFORE the (possibly slow) send, so capture never
         * starves while we wait on the network */
        {
            struct v4l2_buffer qb;
            memset(&qb, 0, sizeof(qb));
            qb.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
            qb.memory = V4L2_MEMORY_MMAP;
            qb.index = buf.index;
            if (xioctl(g_vfd, VIDIOC_QBUF, &qb) < 0) {
                fprintf(stderr, "VIDIOC_QBUF: %s\n", strerror(errno));
                break;
            }
        }

        msg[0] = 'M'; msg[1] = 'V'; msg[2] = 'F'; msg[3] = 'R';
        put32(msg + 4, camid);
        put32(msg + 8, w);
        put32(msg + 12, h);
        put64(msg + 16, seq);
        {
            unsigned long long bits;
            memcpy(&bits, &now, 8);
            put64(msg + 24, bits);
        }
        r = send_frame(sock, msg, msglen, &backlogged);
        if (r < 0) {
            fprintf(stderr, "hub connection lost\n");
            say("camera %u lost the hub", camid);
            teardown();
            free(msg);
            close(sock);
            return 1;
        }
        if (r == 1) {
            dropped++;
            continue;
        }
        seq++;
        if (seq % 50 == 0)
            printf("sent %llu frames (%ux%u), dropped %llu\n", seq, w,
                   h, dropped);
    }

    teardown();
    close(sock);
    free(msg);
    printf("stopped: sent %llu frames, dropped %llu, short-reads "
           "%llu\n", seq, dropped, skipped_short);
    return 0;
}

#else /* !__linux__ */
#include <stdio.h>
int main(void)
{
    fprintf(stderr, "stream_cam_v4l2 is Linux-only; use stream_cam "
                    "(Swift) on macOS\n");
    return 1;
}
#endif
