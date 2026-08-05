/* Protocol torture tool for livehub: a deterministic malformed-traffic
 * client that reproduces, on demand, the wire damage a real WiFi
 * deployment produces -- torn headers, byte-at-a-time delivery, noise
 * bursts, absurd dimensions, interleaved streams, floods and stalls --
 * so parser and robustness bugs surface here rather than in the field.
 *
 * OWNERSHIP (parallel build): this file belongs to the nettest item.
 *
 * Frame protocol under test (little-endian throughout):
 *   "MVFR" | u32 camid | u32 w | u32 h | u64 seq | f64 t_mono
 *          | w*h bytes of 8-bit gray
 * The 32-byte header is followed immediately by the payload.  This tool
 * is written against the PROTOCOL, not against any one implementation of
 * it, so it stays valid as the hub is edited.
 *
 * Usage:
 *   nettest <host> <port> [seed]   run the battery against a live hub
 *   nettest --selftest             offline check of the frame encoder
 *                                  against a built-in parser (no hub)
 *   nettest --fast-drip <n> ...    cap scenario 2's byte-at-a-time drip
 *                                  at n payload bytes (default: whole
 *                                  frame, which takes about 77 s)
 *
 * Everything is deterministic given the seed (default 42); the noise
 * generator is the repo's standard LCG.  Each scenario runs on its own
 * fresh TCP connection and reports what it sent and whether the hub kept
 * reading or dropped the connection.  The bottom line is whether the hub
 * is still alive and parsing afterwards.
 *
 * Portable POSIX C: macOS and Ubuntu alike.
 * Build: make nettest */

#define _POSIX_C_SOURCE 200809L
#ifdef __APPLE__
#define _DARWIN_C_SOURCE 1
#endif

#include <errno.h>
#include <netdb.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/select.h>
#include <sys/socket.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>

#define HDR 32
#define DRIP_W 320
#define DRIP_H 240
#define FLOOD_W 640
#define FLOOD_H 480
#define FLOOD_N 200
#define LORIS_S 10
#define MAXSCEN 32

/* verdicts */
#define V_ALIVE 1
#define V_DROP 0
#define V_CHATTY 2
#define V_ERR (-1)

/* ------------------------------------------------------------------ */
/* deterministic noise: the repo's standard LCG                        */

static unsigned long long rng_next(unsigned long long *s)
{
    *s = *s * 6364136223846793005ULL + 1442695040888963407ULL;
    return *s;
}

static unsigned char rng_byte(unsigned long long *s)
{
    return (unsigned char)((rng_next(s) >> 33) & 255u);
}

/* ------------------------------------------------------------------ */
/* frame encoding                                                      */

static void put32(unsigned char *p, unsigned v)
{
    p[0] = (unsigned char)(v & 255u);
    p[1] = (unsigned char)((v >> 8) & 255u);
    p[2] = (unsigned char)((v >> 16) & 255u);
    p[3] = (unsigned char)((v >> 24) & 255u);
}

static unsigned get32(const unsigned char *p)
{
    return (unsigned)p[0] | ((unsigned)p[1] << 8)
         | ((unsigned)p[2] << 16) | ((unsigned)p[3] << 24);
}

static void put64(unsigned char *p, unsigned long long v)
{
    int i;
    for (i = 0; i < 8; i++)
        p[i] = (unsigned char)((v >> (8 * i)) & 255u);
}

static unsigned long long get64(const unsigned char *p)
{
    unsigned long long v = 0;
    int i;
    for (i = 7; i >= 0; i--)
        v = (v << 8) | (unsigned long long)p[i];
    return v;
}

/* the wire carries the raw IEEE-754 bits of t_mono, little-endian */
static void putf64(unsigned char *p, double t)
{
    unsigned long long bits;
    memcpy(&bits, &t, 8);
    put64(p, bits);
}

static double getf64(const unsigned char *p)
{
    unsigned long long bits = get64(p);
    double t;
    memcpy(&t, &bits, 8);
    return t;
}

static void hdr_make(unsigned char *hdr, const char *magic, unsigned camid,
                     unsigned w, unsigned h, unsigned long long seq,
                     double t)
{
    memcpy(hdr, magic, 4);
    put32(hdr + 4, camid);
    put32(hdr + 8, w);
    put32(hdr + 12, h);
    put64(hdr + 16, seq);
    putf64(hdr + 24, t);
}

/* synthetic gradient: deterministic in (x, y, seq), never all-zero */
static void fill_gradient(unsigned char *p, unsigned w, unsigned h,
                          unsigned long long seq)
{
    unsigned x, y, s = (unsigned)(seq * 7u);
    for (y = 0; y < h; y++)
        for (x = 0; x < w; x++)
            p[(size_t)y * w + x] = (unsigned char)((x * 3u + y * 5u + s)
                                                   & 255u);
}

/* a complete frame (header + gradient payload) in one malloc'd block */
static unsigned char *frame_make(unsigned camid, unsigned w, unsigned h,
                                 unsigned long long seq, double t,
                                 size_t *len)
{
    size_t n = HDR + (size_t)w * h;
    unsigned char *f = malloc(n);
    if (!f)
        return NULL;
    hdr_make(f, "MVFR", camid, w, h, seq, t);
    fill_gradient(f + HDR, w, h, seq);
    *len = n;
    return f;
}

static double now_mono(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

static void sleep_ms(double ms)
{
    struct timespec d;
    d.tv_sec = (time_t)(ms / 1000.0);
    d.tv_nsec = (long)((ms - 1000.0 * (double)d.tv_sec) * 1e6);
    if (d.tv_nsec < 0)
        d.tv_nsec = 0;
    nanosleep(&d, NULL);
}

/* ------------------------------------------------------------------ */
/* tiny reference parser: an independent implementation of the frame    */
/* protocol, used by --selftest to validate the encoder above without   */
/* any hub running.  Deliberately a separate copy -- if it and the      */
/* encoder disagree, one of them is wrong.                              */

typedef void (*tp_cb)(void *ud, unsigned camid, unsigned w, unsigned h,
                      unsigned long long seq, double t,
                      const unsigned char *pix);

typedef struct {
    unsigned char *buf;
    size_t len, cap;
    unsigned long frames, resyncs;
    tp_cb on_frame;
    void *ud;
} tinyparse;

static void tp_init(tinyparse *tp, tp_cb cb, void *ud)
{
    memset(tp, 0, sizeof(*tp));
    tp->on_frame = cb;
    tp->ud = ud;
}

static void tp_free(tinyparse *tp)
{
    free(tp->buf);
    tp->buf = NULL;
    tp->len = tp->cap = 0;
}

/* feed an arbitrary chunk; returns -1 only on allocation failure */
static int tp_feed(tinyparse *tp, const unsigned char *data, size_t n)
{
    if (tp->len + n > tp->cap) {
        size_t nc = tp->cap ? tp->cap : 4096;
        unsigned char *nb;
        while (nc < tp->len + n)
            nc *= 2;
        nb = realloc(tp->buf, nc);
        if (!nb)
            return -1;
        tp->buf = nb;
        tp->cap = nc;
    }
    memcpy(tp->buf + tp->len, data, n);
    tp->len += n;

    for (;;) {
        unsigned camid, w, h;
        size_t need;
        if (tp->len < HDR)
            break;
        if (memcmp(tp->buf, "MVFR", 4) != 0) {
            memmove(tp->buf, tp->buf + 1, --tp->len);
            tp->resyncs++;
            continue;
        }
        camid = get32(tp->buf + 4);
        w = get32(tp->buf + 8);
        h = get32(tp->buf + 12);
        if (w < 64 || h < 64 || w > 4096 || h > 4096) {
            memmove(tp->buf, tp->buf + 1, --tp->len);
            tp->resyncs++;
            continue;
        }
        need = HDR + (size_t)w * h;
        if (tp->len < need)
            break;
        if (tp->on_frame)
            tp->on_frame(tp->ud, camid, w, h, get64(tp->buf + 16),
                         getf64(tp->buf + 24), tp->buf + HDR);
        tp->frames++;
        memmove(tp->buf, tp->buf + need, tp->len - need);
        tp->len -= need;
    }
    return 0;
}

/* ------------------------------------------------------------------ */
/* socket helpers                                                      */

static char lasterr[256];

static int connect_hub(const char *host, const char *port)
{
    struct addrinfo hints, *res, *ai;
    struct timeval tv;
    int s = -1, rc, one = 1;

    lasterr[0] = '\0';
    memset(&hints, 0, sizeof(hints));
    hints.ai_family = AF_INET;
    hints.ai_socktype = SOCK_STREAM;
    rc = getaddrinfo(host, port, &hints, &res);
    if (rc != 0) {
        snprintf(lasterr, sizeof(lasterr), "cannot resolve %s:%s (%s)",
                 host, port, gai_strerror(rc));
        return -1;
    }
    for (ai = res; ai; ai = ai->ai_next) {
        s = socket(ai->ai_family, ai->ai_socktype, ai->ai_protocol);
        if (s < 0) {
            snprintf(lasterr, sizeof(lasterr), "socket: %s",
                     strerror(errno));
            continue;
        }
        if (connect(s, ai->ai_addr, ai->ai_addrlen) == 0)
            break;
        snprintf(lasterr, sizeof(lasterr), "cannot connect to %s:%s (%s)",
                 host, port, strerror(errno));
        close(s);
        s = -1;
    }
    freeaddrinfo(res);
    if (s < 0)
        return -1;
    /* byte-at-a-time scenarios must really hit the wire one byte at a
     * time, so Nagle has to go; the send timeout keeps a wedged hub from
     * hanging this tool forever (it is itself a finding). */
    setsockopt(s, IPPROTO_TCP, TCP_NODELAY, &one, sizeof(one));
    tv.tv_sec = 20;
    tv.tv_usec = 0;
    setsockopt(s, SOL_SOCKET, SO_SNDTIMEO, &tv, sizeof(tv));
    return s;
}

/* returns 0 on success, -1 on error (lasterr set) */
static int send_all(int s, const void *buf, size_t n)
{
    const char *p = (const char *)buf;
    while (n > 0) {
        ssize_t k = send(s, p, n, 0);
        if (k <= 0) {
            snprintf(lasterr, sizeof(lasterr), "send: %s",
                     k == 0 ? "zero-length write" : strerror(errno));
            return -1;
        }
        p += k;
        n -= (size_t)k;
    }
    return 0;
}

/* Liveness probe.  The hub never speaks on this socket, so "alive" means
 * the socket is neither at EOF nor reset: a peer that closed makes the
 * socket readable with a zero-length read, a peer that reset makes it
 * readable with ECONNRESET.  With poke set we also write one filler byte
 * (which the hub discards on its resync path) to draw out a reset the
 * FIN did not carry; scenarios whose byte count is part of the test run
 * with poke clear so the wire stays exactly as specified.  Note this
 * only proves the hub kept reading -- a hub that reads and silently
 * discards looks identical from out here. */
static int conn_alive(int s, int poke)
{
    unsigned char probe = 0x00, peek;
    int i;
    for (i = 0; i < 2; i++) {
        fd_set rd;
        struct timeval tv;
        int n;
        FD_ZERO(&rd);
        FD_SET(s, &rd);
        tv.tv_sec = 0;
        tv.tv_usec = 150000;
        n = select(s + 1, &rd, NULL, NULL, &tv);
        if (n > 0) {
            ssize_t k = recv(s, &peek, 1, MSG_PEEK);
            if (k <= 0)
                return V_DROP;
            return V_CHATTY;
        }
        if (n < 0 && errno != EINTR)
            return V_DROP;
        if (poke && send(s, &probe, 1, 0) < 0)
            return V_DROP;
    }
    return V_ALIVE;
}

/* ------------------------------------------------------------------ */
/* reporting                                                           */

static struct {
    const char *name;
    int verdict;
    char note[96];
} results[MAXSCEN];
static int nresults;

static const char *verdict_str(int v)
{
    switch (v) {
    case V_ALIVE:
        return "SURVIVED";
    case V_DROP:
        return "DROPPED";
    case V_CHATTY:
        return "SURVIVED (hub sent unexpected bytes)";
    default:
        return "ERROR";
    }
}

static void scen_head(const char *tag, const char *name, const char *what)
{
    printf("\n--- %s %s\n    sends: %s\n", tag, name, what);
    fflush(stdout);
}

static void scen_done(const char *name, int verdict, const char *note)
{
    printf("    connection: %s%s%s\n", verdict_str(verdict),
           note && note[0] ? " -- " : "", note ? note : "");
    fflush(stdout);
    if (nresults < MAXSCEN) {
        results[nresults].name = name;
        results[nresults].verdict = verdict;
        snprintf(results[nresults].note, sizeof(results[nresults].note),
                 "%s", note ? note : "");
        nresults++;
    }
}

/* close after probing so the hub sees a normal peer departure */
static void scen_close(const char *name, int s, const char *note)
{
    int v = conn_alive(s, 1);
    close(s);
    scen_done(name, v, note);
}

/* same, but leaves the wire byte-exact: no filler byte is written */
static void scen_close_exact(const char *name, int s, const char *note)
{
    int v = conn_alive(s, 0);
    close(s);
    scen_done(name, v, note);
}

static int scen_open(const char *name, const char *host, const char *port)
{
    int s = connect_hub(host, port);
    if (s < 0)
        scen_done(name, V_ERR, lasterr);
    return s;
}

/* ------------------------------------------------------------------ */
/* the battery                                                         */

/* 1. baseline: one well-formed frame at the smallest legal size */
static void sc_valid_small(const char *host, const char *port)
{
    const char *nm = "valid-small";
    unsigned char *f;
    size_t n;
    int s;
    scen_head("1.", nm, "one well-formed 64x64 gradient frame, camid 1");
    s = scen_open(nm, host, port);
    if (s < 0)
        return;
    f = frame_make(1, 64, 64, 0, now_mono(), &n);
    if (!f) {
        close(s);
        scen_done(nm, V_ERR, "out of memory");
        return;
    }
    if (send_all(s, f, n) != 0) {
        free(f);
        close(s);
        scen_done(nm, V_DROP, lasterr);
        return;
    }
    free(f);
    printf("    wrote %lu bytes (32 header + 4096 payload)\n",
           (unsigned long)n);
    scen_close(nm, s, "baseline: a healthy hub must accept this");
}

/* 2. one byte at a time: partial-header and partial-payload reassembly */
static void sc_drip_feed(const char *host, const char *port, long cap)
{
    const char *nm = "drip-feed";
    unsigned char *f;
    size_t n, i, ndrip;
    double t0;
    int s;
    char note[96];
    scen_head("2.", nm, "a valid 320x240 frame, 1 byte per write, 1 ms "
                        "apart");
    s = scen_open(nm, host, port);
    if (s < 0)
        return;
    f = frame_make(1, DRIP_W, DRIP_H, 1, now_mono(), &n);
    if (!f) {
        close(s);
        scen_done(nm, V_ERR, "out of memory");
        return;
    }
    ndrip = n;
    if (cap > 0 && (size_t)cap + HDR < n)
        ndrip = (size_t)cap + HDR;
    printf("    dripping %lu of %lu bytes (about %.0f s), then the rest "
           "in bulk\n", (unsigned long)ndrip, (unsigned long)n,
           (double)ndrip / 1000.0);
    fflush(stdout);
    t0 = now_mono();
    for (i = 0; i < ndrip; i++) {
        if (send_all(s, f + i, 1) != 0) {
            char msg[96];
            snprintf(msg, sizeof(msg), "dropped after %lu of %lu bytes",
                     (unsigned long)i, (unsigned long)n);
            free(f);
            close(s);
            scen_done(nm, V_DROP, msg);
            return;
        }
        sleep_ms(1.0);
    }
    if (ndrip < n && send_all(s, f + ndrip, n - ndrip) != 0) {
        free(f);
        close(s);
        scen_done(nm, V_DROP, lasterr);
        return;
    }
    free(f);
    snprintf(note, sizeof(note), "%lu single-byte writes in %.1f s",
             (unsigned long)ndrip, now_mono() - t0);
    scen_close(nm, s, note);
}

/* 3. a header that stops mid-field, then the peer vanishes */
static void sc_torn_header(const char *host, const char *port)
{
    const char *nm = "torn-header";
    unsigned char hdr[HDR];
    int s;
    scen_head("3.", nm, "17 bytes of a valid header, then close");
    s = scen_open(nm, host, port);
    if (s < 0)
        return;
    hdr_make(hdr, "MVFR", 1, 640, 480, 2, now_mono());
    if (send_all(s, hdr, 17) != 0) {
        close(s);
        scen_done(nm, V_DROP, lasterr);
        return;
    }
    printf("    wrote 17 bytes (header truncated inside the seq field)\n");
    scen_close_exact(nm, s,
                     "hub must not block waiting for the missing 15 B");
}

/* 4. noise burst then a good frame: the resync path */
static void sc_garbage_then_valid(const char *host, const char *port,
                                  unsigned long long *rng)
{
    const char *nm = "garbage-then-valid";
    unsigned char noise[4096], *f;
    size_t n, i;
    int s;
    scen_head("4.", nm, "4 KB of LCG noise, then a valid 128x128 frame");
    s = scen_open(nm, host, port);
    if (s < 0)
        return;
    for (i = 0; i < sizeof(noise); i++)
        noise[i] = rng_byte(rng);
    f = frame_make(1, 128, 128, 3, now_mono(), &n);
    if (!f) {
        close(s);
        scen_done(nm, V_ERR, "out of memory");
        return;
    }
    if (send_all(s, noise, sizeof(noise)) != 0
        || send_all(s, f, n) != 0) {
        free(f);
        close(s);
        scen_done(nm, V_DROP, lasterr);
        return;
    }
    free(f);
    printf("    wrote %lu bytes of noise + %lu bytes of frame\n",
           (unsigned long)sizeof(noise), (unsigned long)n);
    scen_close(nm, s, "hub must resync and land on the good frame");
}

/* 5. right shape, wrong magic, over and over */
static void sc_bad_magic_flood(const char *host, const char *port)
{
    const char *nm = "bad-magic-flood";
    unsigned char hdr[HDR];
    int s, i;
    scen_head("5.", nm, "100 well-formed headers whose magic is "
                        "\"MVFQ\"");
    s = scen_open(nm, host, port);
    if (s < 0)
        return;
    for (i = 0; i < 100; i++) {
        hdr_make(hdr, "MVFQ", 1, 320, 240, (unsigned long long)i,
                 now_mono());
        if (send_all(s, hdr, HDR) != 0) {
            close(s);
            scen_done(nm, V_DROP, lasterr);
            return;
        }
    }
    printf("    wrote %d bytes, none of it a frame\n", 100 * HDR);
    scen_close(nm, s, "hub must discard, not accumulate forever");
}

/* 6. dimensions no camera could produce */
static void sc_absurd_dims(const char *host, const char *port)
{
    static const struct {
        const char *tag, *name, *what;
        unsigned w, h;
        size_t pay;
        int shortclose;
    } cases[] = {
        { "6a.", "absurd-dims/zero", "valid magic, w=0 h=0",
          0, 0, 0, 0 },
        { "6b.", "absurd-dims/huge", "valid magic, w=100000 h=100000",
          100000, 100000, 0, 0 },
        { "6c.", "absurd-dims/starved",
          "valid magic, w=4096 h=4096, only 10 payload bytes, then close",
          4096, 4096, 10, 1 }
    };
    unsigned char hdr[HDR], pay[16];
    size_t i, j;
    for (i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
        int s;
        scen_head(cases[i].tag, cases[i].name, cases[i].what);
        s = scen_open(cases[i].name, host, port);
        if (s < 0)
            continue;
        hdr_make(hdr, "MVFR", 1, cases[i].w, cases[i].h, 4 + i,
                 now_mono());
        for (j = 0; j < sizeof(pay); j++)
            pay[j] = (unsigned char)(j * 17u);
        if (send_all(s, hdr, HDR) != 0
            || (cases[i].pay
                && send_all(s, pay, cases[i].pay) != 0)) {
            close(s);
            scen_done(cases[i].name, V_DROP, lasterr);
            continue;
        }
        printf("    wrote %lu bytes\n",
               (unsigned long)(HDR + cases[i].pay));
        scen_close_exact(cases[i].name, s,
                         cases[i].shortclose
                         ? "16 MB promised, 10 B delivered, then EOF"
                         : "header must be rejected, not trusted");
    }
}

/* 7. camera identity abuse: extreme id, and more ids than slots */
static void sc_giant_camid(const char *host, const char *port)
{
    const char *nm1 = "giant-camid";
    const char *nm2 = "camid-exhaustion";
    unsigned char *f;
    size_t n;
    int s, i;
    static const unsigned ids[6] = { 900, 901, 902, 903, 904, 905 };

    scen_head("7a.", nm1, "a valid 64x64 frame with camid 0xFFFFFFFF");
    s = scen_open(nm1, host, port);
    if (s >= 0) {
        f = frame_make(0xFFFFFFFFu, 64, 64, 10, now_mono(), &n);
        if (!f) {
            close(s);
            scen_done(nm1, V_ERR, "out of memory");
        } else if (send_all(s, f, n) != 0) {
            free(f);
            close(s);
            scen_done(nm1, V_DROP, lasterr);
        } else {
            free(f);
            printf("    wrote %lu bytes\n", (unsigned long)n);
            scen_close(nm1, s, "camid is opaque; no id may be special");
        }
    }

    scen_head("7b.", nm2, "six valid frames with six distinct camids "
                          "(hub has 4 camera slots)");
    s = scen_open(nm2, host, port);
    if (s < 0)
        return;
    for (i = 0; i < 6; i++) {
        f = frame_make(ids[i], 64, 64, (unsigned long long)(11 + i),
                       now_mono(), &n);
        if (!f) {
            close(s);
            scen_done(nm2, V_ERR, "out of memory");
            return;
        }
        if (send_all(s, f, n) != 0) {
            free(f);
            close(s);
            scen_done(nm2, V_DROP, lasterr);
            return;
        }
        free(f);
    }
    printf("    wrote 6 frames, camids %u..%u\n", ids[0], ids[5]);
    scen_close(nm2, s, "surplus cameras must be ignored, not fatal");
}

/* 8. two frames' bytes shuffled together on one connection: an
 * out-of-spec sender the hub can only treat as corruption */
static void sc_interleave(const char *host, const char *port)
{
    const char *nm = "interleave-abuse";
    unsigned char *a, *b;
    size_t na, nb, ia = 0, ib = 0, chunk = 64, total = 0;
    int s;
    scen_head("8.", nm, "one 128x128 and one 160x160 frame, alternating "
                        "64-byte chunks");
    s = scen_open(nm, host, port);
    if (s < 0)
        return;
    a = frame_make(1, 128, 128, 20, now_mono(), &na);
    b = frame_make(1, 160, 160, 21, now_mono(), &nb);
    if (!a || !b) {
        free(a);
        free(b);
        close(s);
        scen_done(nm, V_ERR, "out of memory");
        return;
    }
    while (ia < na || ib < nb) {
        size_t k;
        if (ia < na) {
            k = na - ia < chunk ? na - ia : chunk;
            if (send_all(s, a + ia, k) != 0)
                break;
            ia += k;
            total += k;
        }
        if (ib < nb) {
            k = nb - ib < chunk ? nb - ib : chunk;
            if (send_all(s, b + ib, k) != 0)
                break;
            ib += k;
            total += k;
        }
    }
    free(a);
    free(b);
    if (ia < na || ib < nb) {
        char msg[96];
        snprintf(msg, sizeof(msg), "dropped after %lu bytes",
                 (unsigned long)total);
        close(s);
        scen_done(nm, V_DROP, msg);
        return;
    }
    printf("    wrote %lu bytes as %s\n", (unsigned long)total,
           "interleaved 64-byte chunks");
    scen_close(nm, s, "corrupt payload: resync or discard, never crash");
}

/* 9. as fast as the socket will take it: backpressure behaviour */
static void sc_flood(const char *host, const char *port)
{
    const char *nm = "flood";
    unsigned char *f;
    size_t n;
    double t0, dt;
    int s, i, sent = 0;
    char note[96];
    scen_head("9.", nm, "200 valid 640x480 frames back to back");
    s = scen_open(nm, host, port);
    if (s < 0)
        return;
    f = frame_make(1, FLOOD_W, FLOOD_H, 100, now_mono(), &n);
    if (!f) {
        close(s);
        scen_done(nm, V_ERR, "out of memory");
        return;
    }
    t0 = now_mono();
    for (i = 0; i < FLOOD_N; i++) {
        put64(f + 16, (unsigned long long)(100 + i));
        putf64(f + 24, now_mono());
        if (send_all(s, f, n) != 0)
            break;
        sent++;
    }
    dt = now_mono() - t0;
    free(f);
    if (dt <= 0.0)
        dt = 1e-9;
    printf("    wrote %d of %d frames, %.1f MB in %.2f s = %.1f MB/s, "
           "%.1f fps\n", sent, FLOOD_N,
           (double)sent * (double)n / 1048576.0, dt,
           (double)sent * (double)n / 1048576.0 / dt, sent / dt);
    if (sent < FLOOD_N) {
        snprintf(note, sizeof(note), "stalled at frame %d: %s", sent,
                 lasterr);
        close(s);
        scen_done(nm, V_DROP, note);
        return;
    }
    snprintf(note, sizeof(note), "%.1f MB/s offered, hub kept draining",
             (double)sent * (double)n / 1048576.0 / dt);
    scen_close(nm, s, note);
}

/* 10. a connection that says almost nothing for ten seconds */
static void sc_slow_loris(const char *host, const char *port)
{
    const char *nm = "slow-loris";
    unsigned char *f;
    size_t n;
    int s;
    scen_head("10.", nm, "\"MVF\", then silence for 10 s, then the rest "
                         "of a valid 128x128 frame");
    s = scen_open(nm, host, port);
    if (s < 0)
        return;
    f = frame_make(1, 128, 128, 30, now_mono(), &n);
    if (!f) {
        close(s);
        scen_done(nm, V_ERR, "out of memory");
        return;
    }
    if (send_all(s, f, 3) != 0) {
        free(f);
        close(s);
        scen_done(nm, V_DROP, lasterr);
        return;
    }
    printf("    wrote 3 bytes, holding the connection open for %d s\n",
           LORIS_S);
    fflush(stdout);
    sleep_ms(1000.0 * LORIS_S);
    if (send_all(s, f + 3, n - 3) != 0) {
        free(f);
        close(s);
        scen_done(nm, V_DROP, lasterr);
        return;
    }
    free(f);
    printf("    completed the frame (%lu bytes total)\n",
           (unsigned long)n);
    /* LIMITATION: from outside we can only prove this connection stayed
     * open.  Whether the hub kept serving its other connections while
     * this one stalled is not observable here -- run a replaycam against
     * the same hub concurrently and watch its frame counter to check. */
    scen_close(nm, s, "only proves this socket lived; see LIMITATION");
}

/* final: is the hub still there, and still parsing? */
static int sc_alive_check(const char *host, const char *port)
{
    const char *nm = "final-alive-check";
    unsigned char *f;
    size_t n;
    int s, v;
    scen_head("11.", nm, "one well-formed 64x64 frame on a fresh "
                         "connection");
    s = scen_open(nm, host, port);
    if (s < 0)
        return 0;
    f = frame_make(1, 64, 64, 999, now_mono(), &n);
    if (!f) {
        close(s);
        scen_done(nm, V_ERR, "out of memory");
        return 0;
    }
    if (send_all(s, f, n) != 0) {
        free(f);
        close(s);
        scen_done(nm, V_DROP, lasterr);
        return 0;
    }
    free(f);
    printf("    wrote %lu bytes\n", (unsigned long)n);
    v = conn_alive(s, 1);
    close(s);
    scen_done(nm, v, "hub accepted a fresh connection and kept reading");
    return v == V_ALIVE || v == V_CHATTY;
}

/* ------------------------------------------------------------------ */
/* self-test: the encoder against the built-in parser, no hub needed    */

typedef struct {
    int n;
    unsigned camid, w, h;
    unsigned long long seq;
    double t;
    unsigned long long sum;
    int payload_ok;
} capture;

static unsigned long long fnv(const unsigned char *p, size_t n)
{
    unsigned long long h = 14695981039346656037ULL;
    size_t i;
    for (i = 0; i < n; i++) {
        h ^= (unsigned long long)p[i];
        h *= 1099511628211ULL;
    }
    return h;
}

static void cap_frame(void *ud, unsigned camid, unsigned w, unsigned h,
                      unsigned long long seq, double t,
                      const unsigned char *pix)
{
    capture *c = (capture *)ud;
    unsigned char *ref = malloc((size_t)w * h);
    c->n++;
    c->camid = camid;
    c->w = w;
    c->h = h;
    c->seq = seq;
    c->t = t;
    c->sum = fnv(pix, (size_t)w * h);
    c->payload_ok = 0;
    if (ref) {
        fill_gradient(ref, w, h, seq);
        c->payload_ok = memcmp(ref, pix, (size_t)w * h) == 0;
        free(ref);
    }
}

static int st_fail;

static void st_check(int ok, const char *what)
{
    printf("  %-4s %s\n", ok ? "ok" : "FAIL", what);
    if (!ok)
        st_fail++;
}

/* push a buffer through the parser in fixed-size chunks */
static void st_push(tinyparse *tp, const unsigned char *p, size_t n,
                    size_t chunk)
{
    size_t i;
    for (i = 0; i < n; i += chunk) {
        size_t k = n - i < chunk ? n - i : chunk;
        if (tp_feed(tp, p + i, k) != 0) {
            printf("  FAIL parser allocation failed\n");
            st_fail++;
            return;
        }
    }
}

static int selftest(unsigned long long seed)
{
    static const unsigned dims[3][2] = { { 64, 64 }, { 320, 240 },
                                         { 640, 480 } };
    static const size_t chunks[7] = { 1, 2, 3, 5, 7, 13, 65536 };
    unsigned long long rng = seed, s2 = seed;
    unsigned char noise[4096];
    unsigned char hdr[HDR];
    unsigned char *f, *g;
    size_t n, m, i, j;
    capture cap;
    tinyparse tp;
    double t = 1234.56789;

    printf("nettest --selftest: frame encoder vs. built-in parser "
           "(seed %llu)\n", seed);

    /* 1. round-trip every header field and the payload, at each size */
    for (i = 0; i < 3; i++) {
        char msg[96];
        f = frame_make(7u, dims[i][0], dims[i][1], 42, t, &n);
        if (!f) {
            st_check(0, "allocate frame");
            return 1;
        }
        memset(&cap, 0, sizeof(cap));
        tp_init(&tp, cap_frame, &cap);
        st_push(&tp, f, n, 65536);
        snprintf(msg, sizeof(msg), "round-trip %ux%u: one frame, all "
                 "fields exact", dims[i][0], dims[i][1]);
        st_check(cap.n == 1 && cap.camid == 7u && cap.w == dims[i][0]
                 && cap.h == dims[i][1] && cap.seq == 42
                 && cap.t == t && cap.payload_ok
                 && n == HDR + (size_t)dims[i][0] * dims[i][1]
                 && tp.len == 0, msg);
        tp_free(&tp);
        free(f);
    }

    /* 2. reassembly: the same frame split at every awkward boundary */
    f = frame_make(3u, 320, 240, 5, t, &n);
    if (!f) {
        st_check(0, "allocate frame");
        return 1;
    }
    for (j = 0; j < sizeof(chunks) / sizeof(chunks[0]); j++) {
        char msg[96];
        memset(&cap, 0, sizeof(cap));
        tp_init(&tp, cap_frame, &cap);
        st_push(&tp, f, n, chunks[j]);
        snprintf(msg, sizeof(msg),
                 "reassembly in %lu-byte chunks: exactly one frame",
                 (unsigned long)chunks[j]);
        st_check(cap.n == 1 && cap.payload_ok && tp.resyncs == 0
                 && tp.len == 0, msg);
        tp_free(&tp);
    }

    /* 3. noise burst then a good frame: the resync path finds it */
    for (i = 0; i < sizeof(noise); i++)
        noise[i] = rng_byte(&rng);
    memset(&cap, 0, sizeof(cap));
    tp_init(&tp, cap_frame, &cap);
    st_push(&tp, noise, sizeof(noise), 512);
    st_check(cap.n == 0, "4 KB of LCG noise yields no frame");
    st_push(&tp, f, n, 512);
    st_check(cap.n == 1 && cap.payload_ok && cap.seq == 5,
             "noise then a valid frame: resync recovers the frame");
    st_check(tp.resyncs >= sizeof(noise) - HDR,
             "resync consumed the noise a byte at a time");
    tp_free(&tp);

    /* 4. wrong magic is never a frame */
    memset(&cap, 0, sizeof(cap));
    tp_init(&tp, cap_frame, &cap);
    for (i = 0; i < 100; i++) {
        hdr_make(hdr, "MVFQ", 1, 320, 240, i, t);
        st_push(&tp, hdr, HDR, HDR);
    }
    st_check(cap.n == 0, "100 \"MVFQ\" headers yield no frame");
    tp_free(&tp);

    /* 5. absurd dimensions are rejected without trusting the length */
    {
        static const unsigned bad[3][2] = { { 0, 0 }, { 100000, 100000 },
                                            { 4096, 4096 } };
        for (i = 0; i < 3; i++) {
            char msg[96];
            unsigned char pay[10];
            memset(&cap, 0, sizeof(cap));
            memset(pay, 0xAB, sizeof(pay));
            tp_init(&tp, cap_frame, &cap);
            hdr_make(hdr, "MVFR", 1, bad[i][0], bad[i][1], 1, t);
            st_push(&tp, hdr, HDR, HDR);
            st_push(&tp, pay, sizeof(pay), sizeof(pay));
            snprintf(msg, sizeof(msg),
                     "dims %ux%u with 10 payload bytes: no frame, no "
                     "over-read", bad[i][0], bad[i][1]);
            st_check(cap.n == 0 && tp.len <= HDR + sizeof(pay), msg);
            tp_free(&tp);
        }
    }

    /* 6. interleaved frames: corruption, but bounded and terminating */
    g = frame_make(1u, 160, 160, 6, t, &m);
    if (!g) {
        st_check(0, "allocate frame");
        free(f);
        return 1;
    }
    {
        size_t ia = 0, ib = 0;
        memset(&cap, 0, sizeof(cap));
        tp_init(&tp, cap_frame, &cap);
        while (ia < n || ib < m) {
            size_t k;
            if (ia < n) {
                k = n - ia < 64 ? n - ia : 64;
                st_push(&tp, f + ia, k, 64);
                ia += k;
            }
            if (ib < m) {
                k = m - ib < 64 ? m - ib : 64;
                st_push(&tp, g + ib, k, 64);
                ib += k;
            }
        }
        st_check(cap.n <= 2 && tp.len <= n + m,
                 "interleaved frames: parser terminates, buffer bounded");
        if (cap.n > 0)
            st_check(cap.w >= 64 && cap.w <= 4096 && cap.h >= 64
                     && cap.h <= 4096,
                     "any frame recovered from the interleave has sane "
                     "dimensions");
        tp_free(&tp);
    }
    free(g);
    free(f);

    /* 7. little-endian field placement is what the spec says */
    hdr_make(hdr, "MVFR", 0xFFFFFFFFu, 0x01020304u, 64,
             0x0102030405060708ULL, t);
    st_check(hdr[0] == 'M' && hdr[1] == 'V' && hdr[2] == 'F'
             && hdr[3] == 'R', "magic occupies bytes 0..3");
    st_check(hdr[4] == 0xFF && hdr[5] == 0xFF && hdr[6] == 0xFF
             && hdr[7] == 0xFF, "camid 0xFFFFFFFF little-endian at 4");
    st_check(hdr[8] == 0x04 && hdr[9] == 0x03 && hdr[10] == 0x02
             && hdr[11] == 0x01, "w little-endian at 8");
    st_check(hdr[16] == 0x08 && hdr[23] == 0x01,
             "seq u64 little-endian at 16");
    st_check(getf64(hdr + 24) == t, "t_mono f64 round-trips bit-exactly");

    /* 8. the noise generator is the repo LCG and is seed-deterministic */
    {
        unsigned long long a, b;
        for (i = 0; i < sizeof(noise); i++)
            noise[i] = rng_byte(&s2);
        a = fnv(noise, sizeof(noise));
        s2 = seed;
        for (i = 0; i < sizeof(noise); i++)
            noise[i] = rng_byte(&s2);
        b = fnv(noise, sizeof(noise));
        st_check(a == b, "noise is reproducible from the seed");
        s2 = 42;
        st_check(rng_next(&s2) == 42ULL * 6364136223846793005ULL
                                  + 1442695040888963407ULL,
                 "LCG uses the repo constants 6364136223846793005 / "
                 "1442695040888963407");
        s2 = 42;
        for (i = 0; i < sizeof(noise); i++)
            noise[i] = rng_byte(&s2);
        st_check(fnv(noise, sizeof(noise)) == 9582443110044643055ULL,
                 "seed 42 produces the expected 4 KB noise block");
    }

    printf("%s: %d check%s failed\n", st_fail ? "SELFTEST FAILED"
                                              : "SELFTEST PASSED",
           st_fail, st_fail == 1 ? "" : "s");
    return st_fail ? 1 : 0;
}

/* ------------------------------------------------------------------ */

static void usage(const char *argv0)
{
    fprintf(stderr,
            "usage: %s <host> <port> [seed]     torture a running hub\n"
            "       %s --selftest [seed]        offline encoder check\n"
            "options: --fast-drip <n>  cap scenario 2's byte-at-a-time\n"
            "                          drip at n payload bytes\n",
            argv0, argv0);
}

int main(int argc, char **argv)
{
    const char *host = NULL, *port = NULL;
    unsigned long long seed = 42, rng;
    long drip_cap = 0;
    int i, self = 0, alive, npos = 0;

    for (i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--selftest") == 0) {
            self = 1;
        } else if (strcmp(argv[i], "--fast-drip") == 0 && i + 1 < argc) {
            drip_cap = atol(argv[++i]);
        } else if (strcmp(argv[i], "-h") == 0
                   || strcmp(argv[i], "--help") == 0) {
            usage(argv[0]);
            return 0;
        } else if (argv[i][0] == '-' && argv[i][1] != '\0') {
            fprintf(stderr, "unknown option %s\n", argv[i]);
            usage(argv[0]);
            return 2;
        } else {
            if (npos == 0 && !self)
                host = argv[i];
            else if (npos == 1 && !self)
                port = argv[i];
            else
                seed = strtoull(argv[i], NULL, 10);
            npos++;
        }
    }
    if (self)
        return selftest(seed);
    if (!host || !port) {
        usage(argv[0]);
        return 2;
    }
    signal(SIGPIPE, SIG_IGN);
    rng = seed;

    printf("nettest: protocol torture against %s:%s (seed %llu)\n", host,
           port, seed);
    printf("protocol: \"MVFR\" | u32 camid | u32 w | u32 h | u64 seq | "
           "f64 t_mono | w*h gray\n");
    printf("each scenario runs on its own fresh connection; \"SURVIVED\""
           " means the hub\nkept reading, \"DROPPED\" means it closed or "
           "reset on us.\n");
    printf("expect minutes, not seconds: the drip is %s,\nthe flood runs"
           " until the hub drains %.0f MB, the loris waits %d s.\n",
           drip_cap > 0 ? "capped" : "77 s of 1 ms writes",
           (double)FLOOD_N * (HDR + (double)FLOOD_W * FLOOD_H)
           / 1048576.0, LORIS_S);

    /* preflight: fail cleanly if nothing is listening */
    {
        int s = connect_hub(host, port);
        if (s < 0) {
            fflush(stdout);
            fprintf(stderr, "\nnettest: %s\n", lasterr);
            fprintf(stderr, "nettest: no hub answered -- start one with "
                            "./livehub <port> <pitch_mm>\n");
            return 1;
        }
        close(s);
        printf("preflight: hub is accepting connections\n");
    }

    sc_valid_small(host, port);
    sc_drip_feed(host, port, drip_cap);
    sc_torn_header(host, port);
    sc_garbage_then_valid(host, port, &rng);
    sc_bad_magic_flood(host, port);
    sc_absurd_dims(host, port);
    sc_giant_camid(host, port);
    sc_interleave(host, port);
    sc_flood(host, port);
    sc_slow_loris(host, port);
    alive = sc_alive_check(host, port);

    printf("\nsummary\n");
    for (i = 0; i < nresults; i++) {
        printf("  %-22s %s\n", results[i].name,
               verdict_str(results[i].verdict));
        /* repeat the detail only where something went wrong, so a bad
         * run is diagnosable from the summary alone */
        if (results[i].verdict != V_ALIVE && results[i].note[0])
            printf("  %-22s   %s\n", "", results[i].note);
    }
    printf("\nBOTTOM LINE: hub %s\n", alive ? "ALIVE" : "DEAD");
    if (alive)
        printf("             (it accepted a fresh connection and a valid"
               " frame after the battery)\n");
    return alive ? 0 : 1;
}
