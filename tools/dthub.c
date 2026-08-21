/* dthub -- the ALIGNED_ASSESSMENT.md section-8 measurement harness.
 *
 * A stripped-down hub that answers ONE question: what is the
 * frame-arrival delta-t distribution between two cameras at full
 * capture rate, BEFORE any decode?  It speaks the same wire protocol
 * as livehub (MVFR frame headers + payload, MVLG log lines), stamps
 * CLOCK_MONOTONIC the moment each 32-byte MVFR header is complete,
 * discards every pixel unread, and prints the distribution at the end
 * (tools/dtstats.c).  It sends nothing back: no MVAK acks, no MVPB
 * clock probes -- both stream_cam and replaycam poll their hub channel
 * with MSG_DONTWAIT and are content with silence.
 *
 * Usage:  dthub <port> [duration-s] [arrivals.csv]
 *   duration-s   stop and report after this many seconds (0 or absent:
 *                run until Ctrl-C, then report)
 *   arrivals.csv stream "camid,seq,t_arrival_s" rows as frames land,
 *                for offline re-analysis
 *
 * Run it exactly like livehub in the deployment runbook, at the
 * cameras' full rate:
 *   hub:     dthub 9900 120 dt.csv
 *   cam N:   ./stream_cam <hub-ip> 9900 <camid> 0
 * Loopback rehearsal without cameras: genframes + two replaycams at a
 * short interval (see RETURN.md, "the next experiment").
 *
 * Deliberately single-threaded (one select loop): every timestamp is
 * taken by the same thread on the same clock, so the measurement can
 * not be skewed by scheduler interleaving between per-camera threads.
 * Payload bytes are drained in bounded chunks per loop turn so one
 * camera's bulk transfer cannot starve the other's header stamp. */

#define _POSIX_C_SOURCE 200809L
#define _DARWIN_C_SOURCE 1

#include <arpa/inet.h>
#include <errno.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/select.h>
#include <sys/socket.h>
#include <time.h>
#include <unistd.h>

#include "dtstats.h"

#define MAXCONN   8
#define HDRMAX    32
#define DRAINCHUNK 65536
#define MAXDIM    16384

typedef struct {
    int fd;                     /* -1 = free slot */
    unsigned char hdr[HDRMAX];  /* header bytes collected so far */
    size_t hlen;                /* how many are valid */
    size_t need;                /* total header bytes wanted (4 until
                                   the magic is known) */
    size_t skip;                /* payload bytes left to discard */
} conn;

static volatile sig_atomic_t stop;

static void on_sigint(int sig)
{
    (void)sig;
    stop = 1;
}

static double now_mono(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

static unsigned get32(const unsigned char *p)
{
    return (unsigned)p[0] | ((unsigned)p[1] << 8)
         | ((unsigned)p[2] << 16) | ((unsigned)p[3] << 24);
}

static unsigned long long get64(const unsigned char *p)
{
    unsigned long long v = 0;
    int i;
    for (i = 7; i >= 0; i--)
        v = (v << 8) | p[i];
    return v;
}

/* header length for a known magic, 0 if unknown */
static size_t hdr_len(const unsigned char *m)
{
    if (memcmp(m, "MVFR", 4) == 0)
        return 32;
    if (memcmp(m, "MVLG", 4) == 0)
        return 12;
    if (memcmp(m, "MVTS", 4) == 0)
        return 24;   /* never solicited, but tolerate it */
    return 0;
}

int main(int argc, char **argv)
{
    int lfd, port, i;
    double duration = 0.0, t0, tlast_status;
    FILE *csv = NULL;
    conn cn[MAXCONN];
    dtstats *st;
    unsigned long long nframes = 0;
    struct sockaddr_in sa;
    int one = 1;

    if (argc < 2 || argc > 4) {
        fprintf(stderr,
                "usage: %s <port> [duration-s] [arrivals.csv]\n",
                argv[0]);
        return 1;
    }
    port = atoi(argv[1]);
    if (argc >= 3)
        duration = atof(argv[2]);
    if (argc >= 4) {
        csv = fopen(argv[3], "w");
        if (!csv) {
            fprintf(stderr, "cannot write %s\n", argv[3]);
            return 1;
        }
        fprintf(csv, "camid,seq,t_arrival_s\n");
    }
    signal(SIGINT, on_sigint);
    signal(SIGTERM, on_sigint);
    signal(SIGPIPE, SIG_IGN);

    st = dt_new();
    if (!st) {
        fprintf(stderr, "out of memory\n");
        return 1;
    }
    for (i = 0; i < MAXCONN; i++) {
        cn[i].fd = -1;
        cn[i].hlen = 0;
        cn[i].need = 4;
        cn[i].skip = 0;
    }

    lfd = socket(AF_INET, SOCK_STREAM, 0);
    if (lfd < 0) {
        fprintf(stderr, "socket: %s\n", strerror(errno));
        return 1;
    }
    setsockopt(lfd, SOL_SOCKET, SO_REUSEADDR, &one, sizeof(one));
    memset(&sa, 0, sizeof(sa));
    sa.sin_family = AF_INET;
    sa.sin_addr.s_addr = htonl(INADDR_ANY);
    sa.sin_port = htons((unsigned short)port);
    if (bind(lfd, (struct sockaddr *)&sa, sizeof(sa)) != 0
        || listen(lfd, 8) != 0) {
        fprintf(stderr, "bind/listen on %d: %s\n", port,
                strerror(errno));
        return 1;
    }
    printf("dthub: listening on %d%s; measuring arrival dt, "
           "discarding pixels\n",
           port, duration > 0.0 ? " (timed run)" : " (Ctrl-C to stop)");
    fflush(stdout);

    t0 = now_mono();
    tlast_status = t0;
    while (!stop) {
        fd_set rd;
        struct timeval tv;
        int maxfd = lfd, rv;
        double tnow = now_mono();

        if (duration > 0.0 && tnow - t0 >= duration)
            break;
        FD_ZERO(&rd);
        FD_SET(lfd, &rd);
        for (i = 0; i < MAXCONN; i++)
            if (cn[i].fd >= 0) {
                FD_SET(cn[i].fd, &rd);
                if (cn[i].fd > maxfd)
                    maxfd = cn[i].fd;
            }
        tv.tv_sec = 0;
        tv.tv_usec = 250000;
        rv = select(maxfd + 1, &rd, NULL, NULL, &tv);
        if (rv < 0) {
            if (errno == EINTR)
                continue;
            fprintf(stderr, "select: %s\n", strerror(errno));
            break;
        }
        if (FD_ISSET(lfd, &rd)) {
            int fd = accept(lfd, NULL, NULL);
            if (fd >= 0) {
                for (i = 0; i < MAXCONN && cn[i].fd >= 0; i++)
                    ;
                if (i == MAXCONN) {
                    close(fd);
                } else {
                    setsockopt(fd, IPPROTO_TCP, TCP_NODELAY, &one,
                               sizeof(one));
                    cn[i].fd = fd;
                    cn[i].hlen = 0;
                    cn[i].need = 4;
                    cn[i].skip = 0;
                    printf("dthub: connection %d up\n", i);
                    fflush(stdout);
                }
            }
        }
        for (i = 0; i < MAXCONN; i++) {
            conn *c = &cn[i];
            int dead = 0;
            if (c->fd < 0 || !FD_ISSET(c->fd, &rd))
                continue;
            if (c->skip > 0) {
                /* discard payload, bounded per turn (fairness) */
                static unsigned char bin[DRAINCHUNK];
                size_t want = c->skip < sizeof(bin) ? c->skip
                                                    : sizeof(bin);
                ssize_t k = recv(c->fd, bin, want, 0);
                if (k <= 0)
                    dead = 1;
                else
                    c->skip -= (size_t)k;
            } else {
                ssize_t k = recv(c->fd, c->hdr + c->hlen,
                                 c->need - c->hlen, 0);
                if (k <= 0) {
                    dead = 1;
                } else {
                    c->hlen += (size_t)k;
                    if (c->hlen == 4 && c->need == 4) {
                        size_t hl = hdr_len(c->hdr);
                        if (hl == 0) {
                            /* resync byte by byte */
                            memmove(c->hdr, c->hdr + 1, 3);
                            c->hlen = 3;
                        } else {
                            c->need = hl;
                        }
                    }
                    if (c->hlen == c->need && c->need > 4) {
                        if (memcmp(c->hdr, "MVFR", 4) == 0) {
                            double ta = now_mono() - t0;
                            unsigned camid = get32(c->hdr + 4);
                            unsigned w = get32(c->hdr + 8);
                            unsigned h = get32(c->hdr + 12);
                            unsigned long long seq = get64(c->hdr + 16);
                            if (w == 0 || h == 0 || w > MAXDIM
                                || h > MAXDIM) {
                                fprintf(stderr, "dthub: insane frame "
                                        "%ux%u, dropping conn %d\n",
                                        w, h, i);
                                dead = 1;
                            } else {
                                c->skip = (size_t)w * (size_t)h;
                                if (dt_add(st, camid, seq, ta) != 0)
                                    fprintf(stderr, "dthub: recorder "
                                            "full, frame dropped\n");
                                if (csv)
                                    fprintf(csv, "%u,%llu,%.6f\n",
                                            camid, seq, ta);
                                nframes++;
                            }
                        } else if (memcmp(c->hdr, "MVLG", 4) == 0) {
                            c->skip = get32(c->hdr + 8);
                        }
                        /* MVTS: fixed size, no payload */
                        c->hlen = 0;
                        c->need = 4;
                    }
                }
            }
            if (dead) {
                printf("dthub: connection %d closed\n", i);
                fflush(stdout);
                close(c->fd);
                c->fd = -1;
            }
        }
        if (tnow - tlast_status >= 5.0) {
            printf("dthub: %llu frames so far, t=%.0f s\n", nframes,
                   tnow - t0);
            fflush(stdout);
            tlast_status = tnow;
        }
    }

    if (csv)
        fclose(csv);
    {
        dt_report r;
        if (dt_analyze(st, &r) == 0) {
            dt_print(&r, stdout);
        } else {
            printf("dthub: not enough data to analyze (need >= 2 "
                   "cameras, >= 8 frames each; got %llu frames)\n",
                   nframes);
        }
    }
    dt_free(st);
    close(lfd);
    for (i = 0; i < MAXCONN; i++)
        if (cn[i].fd >= 0)
            close(cn[i].fd);
    return 0;
}
