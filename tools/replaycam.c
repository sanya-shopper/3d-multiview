/* Replay recorded PGM frames over the live-streaming protocol -- the
 * loopback tester and debugging stand-in for stream_cam: everything
 * downstream (livehub) sees exactly what a live camera would send.
 *
 * Usage: replaycam <hub-host> <port> <camid> <interval-s> <frames.pgm...>
 * Portable POSIX C (macOS and Ubuntu alike).
 *
 * Build: make replaycam */

#include <math.h>
#include <netdb.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/socket.h>
#include <time.h>
#include <unistd.h>

#include "mv/mv.h"

static int send_all(int s, const void *buf, size_t n)
{
    const char *p = (const char *)buf;
    while (n > 0) {
        ssize_t k = send(s, p, n, 0);
        if (k <= 0)
            return -1;
        p += k;
        n -= (size_t)k;
    }
    return 0;
}

static void put32(unsigned char *p, unsigned v)
{
    p[0] = (unsigned char)(v & 255);
    p[1] = (unsigned char)((v >> 8) & 255);
    p[2] = (unsigned char)((v >> 16) & 255);
    p[3] = (unsigned char)((v >> 24) & 255);
}

int main(int argc, char **argv)
{
    struct addrinfo hints, *res;
    unsigned camid;
    double interval;
    int sock, f;
    unsigned long long seq = 0;

    if (argc < 6) {
        fprintf(stderr, "usage: %s <hub-host> <port> <camid> "
                        "<interval-s> <frames.pgm...>\n", argv[0]);
        return 1;
    }
    camid = (unsigned)atoi(argv[3]);
    interval = atof(argv[4]);
    signal(SIGPIPE, SIG_IGN);

    memset(&hints, 0, sizeof(hints));
    hints.ai_family = AF_INET;
    hints.ai_socktype = SOCK_STREAM;
    if (getaddrinfo(argv[1], argv[2], &hints, &res) != 0) {
        fprintf(stderr, "cannot resolve %s\n", argv[1]);
        return 1;
    }
    sock = socket(res->ai_family, res->ai_socktype, res->ai_protocol);
    if (sock < 0 || connect(sock, res->ai_addr, res->ai_addrlen) != 0) {
        fprintf(stderr, "cannot connect to %s:%s\n", argv[1], argv[2]);
        return 1;
    }
    freeaddrinfo(res);
    printf("replaycam %u: %d frames every %.2f s -> %s:%s\n", camid,
           argc - 5, interval, argv[1], argv[2]);

    for (f = 5; f < argc; f++) {
        unsigned char *img, hdr[32];
        int w, h;
        struct timespec ts;
        double t;
        unsigned long long bits;
        int i;
        if (mv_pgm_read(argv[f], &img, &w, &h) != MV_OK) {
            fprintf(stderr, "skip unreadable %s\n", argv[f]);
            continue;
        }
        clock_gettime(CLOCK_MONOTONIC, &ts);
        t = ts.tv_sec + 1e-9 * ts.tv_nsec;
        hdr[0] = 'M'; hdr[1] = 'V'; hdr[2] = 'F'; hdr[3] = 'R';
        put32(hdr + 4, camid);
        put32(hdr + 8, (unsigned)w);
        put32(hdr + 12, (unsigned)h);
        for (i = 0; i < 8; i++)
            hdr[16 + i] = (unsigned char)((seq >> (8 * i)) & 255);
        memcpy(&bits, &t, 8);
        for (i = 0; i < 8; i++)
            hdr[24 + i] = (unsigned char)((bits >> (8 * i)) & 255);
        if (send_all(sock, hdr, 32) != 0
            || send_all(sock, img, (size_t)w * h) != 0) {
            fprintf(stderr, "hub connection lost\n");
            free(img);
            return 1;
        }
        free(img);
        seq++;
        if (seq % 25 == 0)
            printf("  cam %u: sent %llu frames\n", camid, seq);
        {
            struct timespec d;
            d.tv_sec = (time_t)interval;
            d.tv_nsec = (long)((interval - (double)d.tv_sec) * 1e9);
            nanosleep(&d, NULL);
        }
    }
    printf("cam %u: replay done (%llu frames)\n", camid, seq);
    close(sock);
    return 0;
}
