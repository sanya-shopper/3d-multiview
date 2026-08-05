/* Live processing hub: receives camera streams over TCP, blind-decodes
 * the calibration pattern continuously, calibrates each camera as views
 * accumulate (robust Zhang + LM refinement, live), and estimates
 * camera-to-camera extrinsics from counter-matched pattern poses.
 * Portable POSIX C: runs on macOS and Ubuntu unchanged.
 *
 * DEPLOYMENT RUNBOOK (live test, three MacBooks):
 *  1. Pattern MacBook: copy tools/pattern.html over, open in a browser,
 *     press f (fullscreen) then 6 (mux: serves near AND far cameras).
 *     Inhibit sleep (macOS: caffeinate -dimsu).
 *  2. Hub (this machine): make livehub && ./livehub 9900 0.1133 rec/
 *     (port, display pixel pitch in mm, optional record directory).
 *  3. Each camera MacBook: copy the stream_cam binary over
 *     (swiftc -O tools/stream_cam.swift -o stream_cam; binaries are
 *     architecture-specific -- rebuild on Intel Macs) and run
 *         ./stream_cam <hub-ip> 9900 <camid> 5
 *     with distinct camids (1, 2, ...). Grant camera permission.
 *  4. Carry the pattern laptop through the volume with dwells near
 *     each camera (the round-five protocol); watch the hub converge.
 *  5. Ctrl-C the hub when done; if recording, rerun offline tools on
 *     the recorded frames for the full-precision pass.
 *  Loopback test without hardware:  ./livehub 9900 0.1133  plus two
 *  replaycam instances streaming recorded PGM frames.
 *
 * Ubuntu delta: only the capture shim differs (tools/stream_cam_v4l2.c);
 * this hub and replaycam compile unchanged. */

/* glibc hides POSIX declarations (clock_gettime, getaddrinfo,
 * struct timespec) under -std=c99 -pedantic without this; Darwin
 * exposes them unconditionally -- found by deploy/ static audit */
#define _POSIX_C_SOURCE 200809L
#define _DARWIN_C_SOURCE 1

#include <arpa/inet.h>
#include <errno.h>
#include <ifaddrs.h>
#include <math.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <sys/select.h>
#include <sys/socket.h>
#include <sys/stat.h>
#include <stdarg.h>
#include <time.h>
#include <unistd.h>

#include "mv/mv.h"

#define MAXCONN 8
#define MAXCAMS 4
#define MAXV 220
#define MAXANCH 512
#define MAXPAIR 256

typedef struct {
    int used;
    unsigned camid;
    int conn;               /* owning connection index, -1 = none:
                             * binds a camid to one socket so two
                             * cameras sharing an id cannot merge */
    int w, h;               /* dims of `latest` (mutex-guarded) */
    unsigned char *latest;  /* written by network thread (mutex) */
    unsigned char *work;    /* DECODER-OWNED: allocated and read only by
                             * the decoder thread; the network thread
                             * must never touch it */
    int workcap;            /* capacity of work in bytes */
    int busy;               /* decoder mid-decode: slot must not be
                             * reclaimed (mutex-guarded flag) */
    double t_latest;
    int fresh;
    long frames_rx, decodes, reads_ok, valid_ctr;
    /* calibration state */
    double (*objs)[2 * MV_READ_MAXC];
    double (*imgs)[2 * MV_READ_MAXC];
    mv_calib_view views[MAXV];
    int nviews, last_calib_n;
    int calibrated;
    double K[9], kr[2], rms;
    /* extrinsic anchors: display pose at decoded counter. Ring buffer:
     * late, wide-baseline dwells must not be lost once the array fills
     * (the informative poses arrive over minutes). */
    struct anchor {
        double k;      /* counter (full for fine, mod 256 for coarse) */
        int tier;
        int dwell;     /* display static vs previous anchor */
        double t_host; /* hub receive time */
        mv_camera pose;
    } anch[MAXANCH];
    int nanch;         /* valid anchors, capped at MAXANCH */
    int anhead;        /* ring write index */
    int last_tier;
    unsigned last_counter;
    double last_dist;
} cam_t;

typedef struct {
    int sock;
    unsigned char *buf;
    size_t len, cap;
} conn_t;

static cam_t cams[MAXCAMS];
static conn_t conns[MAXCONN];
static pthread_mutex_t mbx = PTHREAD_MUTEX_INITIALIZER;
static double pitch;
static const char *recdir = NULL;   /* the per-session log directory */
static FILE *reclog = NULL;
static FILE *hublog = NULL;          /* session copy of the hub console */
static volatile sig_atomic_t g_stop = 0;

/* print to the console AND the session hub.log, so the log is
 * self-contained regardless of how the hub was launched */
static void hlog(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    vprintf(fmt, ap);
    va_end(ap);
    if (hublog) {
        va_start(ap, fmt);
        vfprintf(hublog, fmt, ap);
        va_end(ap);
        fflush(hublog);
    }
}

static void on_signal(int s)
{
    (void)s;
    g_stop = 1;
}

/* print this Mac's LAN IPv4 address(es) in a banner, so the operator
 * knows exactly what to type on the other camera Mac */
static void print_hub_ips(int port)
{
    struct ifaddrs *ifap, *p;
    printf("\n");
    printf("  +--------------------------------------------------+\n");
    printf("  |  HUB READY -- on the OTHER camera Mac, enter:     |\n");
    if (getifaddrs(&ifap) == 0) {
        for (p = ifap; p; p = p->ifa_next) {
            char ip[INET_ADDRSTRLEN];
            struct sockaddr_in *sa;
            if (!p->ifa_addr || p->ifa_addr->sa_family != AF_INET)
                continue;
            sa = (struct sockaddr_in *)p->ifa_addr;
            inet_ntop(AF_INET, &sa->sin_addr, ip, sizeof(ip));
            if (strcmp(ip, "127.0.0.1") == 0)
                continue;
            printf("  |    hub IP:  %-15s   port %-5d          |\n",
                   ip, port);
        }
        freeifaddrs(ifap);
    }
    printf("  +--------------------------------------------------+\n\n");
}

static double now_mono(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

static unsigned get32(const unsigned char *p)
{
    return (unsigned)p[0] | ((unsigned)p[1] << 8)
         | ((unsigned)p[2] << 16) | ((unsigned)p[3] << 24);
}

static void put32_le(unsigned char *p, unsigned v)
{
    p[0] = (unsigned char)(v & 255);
    p[1] = (unsigned char)((v >> 8) & 255);
    p[2] = (unsigned char)((v >> 16) & 255);
    p[3] = (unsigned char)((v >> 24) & 255);
}

/* append a camera's pushed log line to rec/camlog_<id>.txt (called from
 * the network thread only). This is how remote cameras' diagnostics --
 * fps, drops, errors -- reach the hub without any SSH or file transfer:
 * they ride the same TCP connection as the frames (MVLG messages). */
static void camlog(unsigned camid, const unsigned char *text, int len)
{
    char path[512];
    FILE *f;
    if (!recdir)
        return;
    snprintf(path, sizeof(path), "%s/camlog_%u.txt", recdir, camid);
    f = fopen(path, "a");
    if (!f)
        return;
    if (ftell(f) == 0) { /* new file: identify it */
        time_t nw = time(NULL);
        char st[64];
        strftime(st, sizeof(st), "%Y-%m-%d %H:%M:%S", localtime(&nw));
        fprintf(f, "# multiview camera %u log (pushed to hub) | opened "
                "%s\n", camid, st);
    }
    fprintf(f, "[%.3f] ", now_mono());
    fwrite(text, 1, (size_t)len, f);
    fputc('\n', f);
    fclose(f);
}

/* called from the network thread with mbx HELD. ci = the connection
 * index the frame arrived on; a camid is bound to one connection so a
 * second camera claiming an in-use id is refused, not silently merged. */
static cam_t *cam_slot(unsigned camid, int ci)
{
    int i;
    for (i = 0; i < MAXCAMS; i++)
        if (cams[i].used && cams[i].camid == camid) {
            if (cams[i].conn == ci)
                return &cams[i];
            if (cams[i].conn < 0) { /* prior owner disconnected */
                cams[i].conn = ci;
                return &cams[i];
            }
            hlog("[hub] camid %u already owned by another live "
                 "connection; frame refused\n", camid);
            return NULL;
        }
    for (i = 0; i < MAXCAMS; i++)
        if (!cams[i].used) {
            double (*o)[2 * MV_READ_MAXC] =
                malloc((size_t)MAXV * sizeof(*o));
            double (*m)[2 * MV_READ_MAXC] =
                malloc((size_t)MAXV * sizeof(*m));
            if (!o || !m) { /* do not burn the slot on failure */
                free(o);
                free(m);
                return NULL;
            }
            memset(&cams[i], 0, sizeof(cams[i]));
            cams[i].used = 1;
            cams[i].camid = camid;
            cams[i].conn = ci;
            cams[i].objs = o;
            cams[i].imgs = m;
            hlog("[hub] camera %u joined\n", camid);
            return &cams[i];
        }
    /* table full: reclaim a slot that went quiet (a camera restarted
     * under a new camid must not lock the rig out -- found by nettest
     * slot-exhaustion scenario), else log the refusal visibly */
    {
        int stale = -1;
        double tnow = now_mono();
        for (i = 0; i < MAXCAMS; i++)
            if (!cams[i].busy && tnow - cams[i].t_latest > 60.0
                && (stale < 0
                    || cams[i].t_latest < cams[stale].t_latest))
                stale = i;
        if (stale >= 0) {
            hlog("[hub] camera %u replaces quiet camera %u\n",
                 camid, cams[stale].camid);
            free(cams[stale].latest);
            free(cams[stale].work);
            cams[stale].latest = NULL;
            cams[stale].work = NULL;
            cams[stale].workcap = 0;
            cams[stale].camid = camid;
            cams[stale].conn = ci;
            cams[stale].w = cams[stale].h = 0;
            cams[stale].fresh = 0;
            cams[stale].frames_rx = cams[stale].decodes = 0;
            cams[stale].reads_ok = cams[stale].valid_ctr = 0;
            cams[stale].nviews = cams[stale].last_calib_n = 0;
            cams[stale].calibrated = 0;
            cams[stale].nanch = cams[stale].anhead = 0;
            cams[stale].last_tier = 0;
            return &cams[stale];
        }
    }
    hlog("[hub] camera %u refused: %d slots busy\n", camid, MAXCAMS);
    return NULL;
}

/* undistort corner pixels under the camera's current K, kr */
static void undo_distort(double *out, const double *in, int n,
                         const double K[9], const double kr[2])
{
    double Kinv[9];
    int i, it;
    mv_mat_inv3(Kinv, K);
    for (i = 0; i < n; i++) {
        double u = in[2 * i], v = in[2 * i + 1];
        double xn = Kinv[0] * u + Kinv[1] * v + Kinv[2];
        double yn = Kinv[3] * u + Kinv[4] * v + Kinv[5];
        double x = xn, y = yn;
        for (it = 0; it < 20; it++) {
            double r2 = x * x + y * y;
            double rad = 1.0 + r2 * (kr[0] + r2 * kr[1]);
            if (fabs(rad) < 1e-6)
                break;
            x = xn / rad;
            y = yn / rad;
        }
        out[2 * i] = K[0] * x + K[1] * y + K[2];
        out[2 * i + 1] = K[4] * y + K[5];
    }
}

/* live recalibration once enough new views have arrived */
static void try_calibrate(cam_t *c)
{
    static mv_camera est[MAXV];
    unsigned char inl[MAXV];
    double K[9], kr[2];
    if (c->nviews < 12 || c->nviews - c->last_calib_n < 10)
        return;
    c->last_calib_n = c->nviews;
    if (mv_calib_planar_robust(K, est, kr, c->views, c->nviews, 1, inl)
        != MV_OK)
        return;
    /* LM-polish on the accepted views (the robust fit already set the
     * inlier mask; refine over all views would re-admit bad reads) */
    {
        static mv_calib_view iv[MAXV];
        static mv_camera ie[MAXV];
        int i, m = 0;
        for (i = 0; i < c->nviews; i++)
            if (inl[i]) {
                iv[m] = c->views[i];
                ie[m] = est[i];
                m++;
            }
        if (m >= 6
            && mv_calib_refine(K, ie, kr, iv, m) == MV_OK) {
            double rms = mv_calib_reproj_rms(ie, iv, m);
            pthread_mutex_lock(&mbx);
            memcpy(c->K, K, sizeof(K));
            memcpy(c->kr, kr, sizeof(kr));
            c->rms = rms;
            c->calibrated = 1;
            pthread_mutex_unlock(&mbx);
            hlog("[cal] cam %u: fx %.1f fy %.1f cx %.1f cy %.1f "
                   "k1 %+.3f k2 %+.3f | RMS %.3f px over %d views\n",
                   c->camid, K[0], K[4], K[2], K[5], kr[0], kr[1],
                   rms, m);
        }
    }
}

/* accumulate pairwise relative poses; report chordal mean */
static void try_extrinsics(void)
{
    int a, b;
    for (a = 0; a < MAXCAMS; a++)
        for (b = a + 1; b < MAXCAMS; b++) {
            cam_t *A = &cams[a], *B = &cams[b];
            double Rsum[9] = { 0 }, tacc[3][MAXPAIR];
            int npair = 0, i, j, k;
            if (!A->used || !B->used || !A->calibrated
                || !B->calibrated)
                continue;
            for (i = 0; i < A->nanch && npair < MAXPAIR; i++) {
                for (j = 0; j < B->nanch && npair < MAXPAIR; j++) {
                    struct anchor *pa = &A->anch[i], *pb = &B->anch[j];
                    double dk;
                    int relaxed;
                    if (fabs(pa->t_host - pb->t_host) > 4.0)
                        continue;
                    if (pa->tier == 1 && pb->tier == 1)
                        dk = fabs(pa->k - pb->k);
                    else {
                        double d8 = fmod(fmod(pa->k, 256.0)
                                         - fmod(pb->k, 256.0) + 384.0,
                                         256.0) - 128.0;
                        dk = fabs(d8);
                    }
                    relaxed = pa->dwell && pb->dwell;
                    if (dk > (relaxed ? 120.0 : 2.0))
                        continue;
                    {
                        double Rbt[9], R[9], t[3];
                        mv_mat_transpose(Rbt, pb->pose.R, 3, 3);
                        mv_mat_mul(R, pa->pose.R, Rbt, 3, 3, 3);
                        mv_mat_mul(t, R, pb->pose.t, 3, 3, 1);
                        for (k = 0; k < 3; k++)
                            t[k] = pa->pose.t[k] - t[k];
                        for (k = 0; k < 9; k++)
                            Rsum[k] += R[k];
                        for (k = 0; k < 3; k++)
                            tacc[k][npair] = t[k];
                        npair++;
                    }
                }
            }
            if (npair >= 3) {
                double U[9], S[3], V[9], Vt[9], Rm[9], tm[3];
                double dev = 0.0;
                memcpy(U, Rsum, sizeof(Rsum));
                if (mv_svd(U, S, V, 3, 3) != MV_OK)
                    continue;
                mv_mat_transpose(Vt, V, 3, 3);
                mv_mat_mul(Rm, U, Vt, 3, 3, 3);
                for (k = 0; k < 3; k++) {
                    double col[MAXPAIR];
                    int q, r;
                    memcpy(col, tacc[k],
                           (size_t)npair * sizeof(double));
                    for (q = 1; q < npair; q++)
                        for (r = q; r > 0 && col[r] < col[r - 1]; r--) {
                            double tmp = col[r];
                            col[r] = col[r - 1];
                            col[r - 1] = tmp;
                        }
                    tm[k] = col[npair / 2];
                }
                for (i = 0; i < npair; i++) {
                    double d = 0.0;
                    for (k = 0; k < 3; k++)
                        d += (tacc[k][i] - tm[k])
                             * (tacc[k][i] - tm[k]);
                    dev += sqrt(d);
                }
                hlog("[ext] cam %u <-> cam %u: baseline %.3f m  "
                       "t=(%+.3f %+.3f %+.3f)  %d pairs, mean |dt| "
                       "%.1f mm\n", A->camid, B->camid,
                       sqrt(tm[0] * tm[0] + tm[1] * tm[1]
                            + tm[2] * tm[2]), tm[0], tm[1], tm[2],
                       npair, 1000.0 * dev / npair);
            }
        }
}

/* decode a camera's frame. img/w/h are a decoder-private snapshot (the
 * shared c->work/c->w/c->h must never be read here -- the network
 * thread can resize them concurrently); results are published into c
 * under mbx. */
static void process_cam(cam_t *c, const unsigned char *img, int w, int h)
{
    mv_read_result rr;
    int tier = 1, i;

    int read_ok;
    read_ok = mv_read_pattern(&rr, img, w, h) == MV_OK;
    if (!read_ok && mv_read_coarse(&rr, img, w, h) == MV_OK) {
        read_ok = 1;
        tier = 2;
    }
    pthread_mutex_lock(&mbx);
    c->decodes++;
    if (!read_ok) {
        c->last_tier = 0;
        pthread_mutex_unlock(&mbx);
        return;
    }
    c->reads_ok++;
    c->last_tier = tier;
    c->last_counter = rr.counter_valid ? rr.counter : 0;
    if (rr.counter_valid)
        c->valid_ctr++;
    pthread_mutex_unlock(&mbx);

    if (recdir && reclog) {
        char path[512], camname[32];
        snprintf(path, sizeof(path), "%s/cam%u_%06ld.pgm", recdir,
                 c->camid, c->reads_ok);
        snprintf(camname, sizeof(camname), "%u", c->camid);
        mv_pgm_write(path, img, w, h);
        mv_session_frm(reclog, camname, c->t_latest, (int)c->reads_ok);
        mv_session_read(reclog, camname, c->t_latest, &rr, NULL);
    }

    /* calibration view: skip near-duplicates (dwell frames) once we
     * have a base set, so 220 slots span poses rather than time */
    if ((tier == 1 ? rr.n >= 50 : rr.n >= 8) && c->nviews < MAXV) {
        int take = 1;
        if (c->nviews > 10) {
            const mv_calib_view *lv = &c->views[c->nviews - 1];
            if (lv->n == rr.n) {
                double dm = 0.0;
                for (i = 0; i < rr.n; i++)
                    dm += fabs(lv->img[2 * i] - rr.uv[2 * i])
                        + fabs(lv->img[2 * i + 1] - rr.uv[2 * i + 1]);
                if (dm / rr.n < 3.0)
                    take = 0;
            }
        }
        if (take) {
            for (i = 0; i < rr.n; i++) {
                double xy[2];
                if (tier == 1)
                    mv_pattern_corner_px(rr.id[i] % MV_PAT_CORNER_COLS,
                                         rr.id[i]
                                         / MV_PAT_CORNER_COLS, xy);
                else
                    mv_pattern2_corner_px(rr.id[i]
                                          % MV_PAT2_CORNER_COLS,
                                          rr.id[i]
                                          / MV_PAT2_CORNER_COLS, xy);
                c->objs[c->nviews][2 * i] = xy[0] * pitch;
                c->objs[c->nviews][2 * i + 1] = xy[1] * pitch;
                c->imgs[c->nviews][2 * i] = rr.uv[2 * i];
                c->imgs[c->nviews][2 * i + 1] = rr.uv[2 * i + 1];
            }
            c->views[c->nviews].obj = c->objs[c->nviews];
            c->views[c->nviews].img = c->imgs[c->nviews];
            c->views[c->nviews].n = rr.n;
            pthread_mutex_lock(&mbx);
            c->nviews++;
            pthread_mutex_unlock(&mbx);
            try_calibrate(c);
        }
    }

    /* extrinsic anchor */
    if (c->calibrated && rr.counter_valid) {
        static double obj[2 * MV_READ_MAXC], und[2 * MV_READ_MAXC];
        mv_camera pose;
        for (i = 0; i < rr.n; i++) {
            double xy[2];
            if (tier == 1)
                mv_pattern_corner_px(rr.id[i] % MV_PAT_CORNER_COLS,
                                     rr.id[i] / MV_PAT_CORNER_COLS,
                                     xy);
            else
                mv_pattern2_corner_px(rr.id[i] % MV_PAT2_CORNER_COLS,
                                      rr.id[i] / MV_PAT2_CORNER_COLS,
                                      xy);
            obj[2 * i] = xy[0] * pitch;
            obj[2 * i + 1] = xy[1] * pitch;
        }
        undo_distort(und, rr.uv, rr.n, c->K, c->kr);
        if (mv_calib_plane_pose(&pose, c->K, obj, und, rr.n) == MV_OK) {
            double cx = MV_PAT_W / 2.0 * pitch;
            double cy = MV_PAT_H / 2.0 * pitch;
            struct anchor an;
            an.k = (double)rr.counter;
            an.tier = tier;
            an.t_host = now_mono();
            an.pose = pose;
            an.dwell = 0;
            /* append to the ring under the lock; try_extrinsics reads
             * every camera's anchors/state, so it runs locked too */
            pthread_mutex_lock(&mbx);
            c->last_dist = pose.R[6] * cx + pose.R[7] * cy + pose.t[2];
            if (c->nanch > 0) {
                int prev = (c->anhead - 1 + MAXANCH) % MAXANCH;
                struct anchor *pv = &c->anch[prev];
                double d = 0.0;
                for (i = 0; i < 3; i++)
                    d += (pv->pose.t[i] - pose.t[i])
                         * (pv->pose.t[i] - pose.t[i]);
                an.dwell = sqrt(d) < 0.01
                           && an.t_host - pv->t_host < 5.0;
            }
            c->anch[c->anhead] = an;
            c->anhead = (c->anhead + 1) % MAXANCH;
            if (c->nanch < MAXANCH)
                c->nanch++;
            try_extrinsics();
            pthread_mutex_unlock(&mbx);
        }
    }
}

/* close connection ci and release any camid it owned, so a camera that
 * reconnects (same or new id) is accepted rather than locked out */
static void release_conn(int ci)
{
    int j;
    close(conns[ci].sock);
    conns[ci].sock = -1;
    conns[ci].len = 0;
    pthread_mutex_lock(&mbx);
    for (j = 0; j < MAXCAMS; j++)
        if (cams[j].used && cams[j].conn == ci)
            cams[j].conn = -1;
    pthread_mutex_unlock(&mbx);
}

/* decoder thread: the reader is non-reentrant (static buffers), so one
 * thread serializes all decodes; it always works on each camera's
 * FRESHEST frame, and the network thread never blocks on it */
static void *decoder_main(void *arg)
{
    int rr_next = 0;
    (void)arg;
    while (!g_stop) {
        cam_t *pick = NULL;
        int i, dw = 0, dh = 0;
        pthread_mutex_lock(&mbx);
        for (i = 0; i < MAXCAMS; i++) {
            cam_t *c = &cams[(rr_next + i) % MAXCAMS];
            if (c->used && c->fresh && c->latest) {
                size_t sz = (size_t)c->w * c->h;
                /* work is decoder-owned: grow it here (under the lock,
                 * where c->w/c->h are stable) and never from the
                 * network thread */
                if (c->workcap < (int)sz) {
                    unsigned char *nb = realloc(c->work, sz);
                    if (nb) {
                        c->work = nb;
                        c->workcap = (int)sz;
                    }
                }
                if (c->workcap >= (int)sz) {
                    memcpy(c->work, c->latest, sz);
                    c->fresh = 0;
                    c->busy = 1;
                    dw = c->w;
                    dh = c->h;
                    pick = c;
                    rr_next = ((rr_next + i) % MAXCAMS) + 1;
                }
                break;
            }
        }
        pthread_mutex_unlock(&mbx);
        if (pick) {
            /* pick->work stays valid: busy blocks reclaim, and the
             * network thread never frees/resizes work */
            process_cam(pick, pick->work, dw, dh);
            pthread_mutex_lock(&mbx);
            pick->busy = 0;
            pthread_mutex_unlock(&mbx);
        }
        else {
            struct timespec d = { 0, 50000000L };
            nanosleep(&d, NULL);
        }
    }
    return NULL;
}

static void status(void)
{
    int i;
    for (i = 0; i < MAXCAMS; i++) {
        cam_t *c = &cams[i];
        if (!c->used)
            continue;
        double age = now_mono() - c->t_latest;
        hlog("[cam %u] rx %ld | decoded %ld/%ld | ctr %ld | views %d"
             " | %s", c->camid, c->frames_rx, c->reads_ok,
             c->decodes, c->valid_ctr, c->nviews,
             c->calibrated ? "CALIBRATED" : "collecting");
        if (c->calibrated)
            hlog(" fx %.0f rms %.2f", c->K[0], c->rms);
        if (c->last_tier)
            hlog(" | last tier %d ctr %u", c->last_tier,
                 c->last_counter);
        /* loud flag when a camera that WAS flowing has gone quiet */
        if (c->frames_rx > 0 && age > 2.0)
            hlog(" | *** STALLED: no frames for %.0f s ***", age);
        else
            hlog(" | frame age %.1f s", age);
        hlog("\n");
    }
}

int main(int argc, char **argv)
{
    int port, lsock, i;
    struct sockaddr_in addr;
    double last_status = 0.0;
    pthread_t decoder_th;

    static char sess[600];
    const char *base;
    char host[128], stamp[32], started[64];
    time_t now = time(NULL);
    struct tm *lt = localtime(&now);

    setvbuf(stdout, NULL, _IOLBF, 0); /* live logs must survive kill */
    /* every argument is optional, with the deployment defaults baked in:
     *   hubengine [port] [pitch_mm] [logs_base_dir]   (9900 0.1133 logs) */
    port = argc > 1 ? atoi(argv[1]) : 9900;
    pitch = (argc > 2 ? atof(argv[2]) : 0.1133) * 1e-3;
    base = argc > 3 ? argv[3] : "logs";

    /* one clean, self-identifying directory per session:
     * <base>/<host>-<YYYYMMDD-HHMMSS>/ -- keeps runs from mixing and
     * makes it obvious which log is which, from where, and when */
    if (gethostname(host, sizeof(host)) != 0)
        strcpy(host, "hub");
    { char *d = strchr(host, '.'); if (d) *d = 0; }
    strftime(stamp, sizeof(stamp), "%Y%m%d-%H%M%S", lt);
    strftime(started, sizeof(started), "%Y-%m-%d %H:%M:%S", lt);
    mkdir(base, 0755);
    snprintf(sess, sizeof(sess), "%s/%s-%s", base, host, stamp);
    if (mkdir(sess, 0755) == 0 || 1) {
        char path[700];
        recdir = sess;
        snprintf(path, sizeof(path), "%s/hub.log", sess);
        hublog = fopen(path, "w");
        if (hublog) {
            fprintf(hublog, "# multiview HUB log | host %s | started %s "
                    "| port %d\n", host, started, port);
            fflush(hublog);
        }
        snprintf(path, sizeof(path), "%s/records.txt", sess);
        reclog = fopen(path, "w");
        if (reclog) {
            fprintf(reclog, "# multiview records | host %s | started %s\n",
                    host, started);
            mv_session_ver(reclog, MV_PAT_SPEC_VERSION, "hub");
        }
    }
    signal(SIGPIPE, SIG_IGN);
    {   /* clean shutdown on Ctrl-C / SIGTERM: no SA_RESTART so select
         * returns EINTR and the loop notices g_stop */
        struct sigaction sa;
        memset(&sa, 0, sizeof(sa));
        sa.sa_handler = on_signal;
        sigaction(SIGINT, &sa, NULL);
        sigaction(SIGTERM, &sa, NULL);
    }

    lsock = socket(AF_INET, SOCK_STREAM, 0);
    {
        int one = 1;
        setsockopt(lsock, SOL_SOCKET, SO_REUSEADDR, &one, sizeof(one));
    }
    memset(&addr, 0, sizeof(addr));
    addr.sin_family = AF_INET;
    addr.sin_addr.s_addr = htonl(INADDR_ANY);
    addr.sin_port = htons((unsigned short)port);
    if (bind(lsock, (struct sockaddr *)&addr, sizeof(addr)) != 0
        || listen(lsock, 8) != 0) {
        fprintf(stderr, "cannot listen on port %d\n", port);
        return 1;
    }
    for (i = 0; i < MAXCONN; i++)
        conns[i].sock = -1;
    {
        if (pthread_create(&decoder_th, NULL, decoder_main, NULL) != 0) {
            fprintf(stderr, "cannot start decoder thread\n");
            return 1;
        }
    }
    hlog("[hub] host %s | started %s | listening on %d "
         "(pitch %.4f mm)\n", host, started, port, pitch * 1000.0);
    hlog("[hub] logs -> %s/\n", recdir);
    print_hub_ips(port);
    fflush(stdout);

    while (!g_stop) {
        fd_set rd;
        struct timeval tv;
        int maxfd = lsock, n;
        static double last_ack = 0.0;
        FD_ZERO(&rd);
        FD_SET(lsock, &rd);
        for (i = 0; i < MAXCONN; i++)
            if (conns[i].sock >= 0) {
                FD_SET(conns[i].sock, &rd);
                if (conns[i].sock > maxfd)
                    maxfd = conns[i].sock;
            }
        tv.tv_sec = 0;
        tv.tv_usec = 100000;
        n = select(maxfd + 1, &rd, NULL, NULL, &tv);
        if (n < 0) {
            /* EINTR is routine; other errors (EBADF from a raced fd,
             * transient ENOMEM) must not tear the whole rig down --
             * log and retry rather than exit mid-deployment */
            if (errno != EINTR)
                fprintf(stderr, "[hub] select: %s (continuing)\n",
                        strerror(errno));
            continue;
        }

        if (n > 0 && FD_ISSET(lsock, &rd)) {
            int s = accept(lsock, NULL, NULL);
            if (s >= 0) {
                for (i = 0; i < MAXCONN; i++)
                    if (conns[i].sock < 0) {
                        if (!conns[i].buf) {
                            conns[i].cap = 1 << 22;
                            conns[i].buf = malloc(conns[i].cap);
                            if (!conns[i].buf) { /* refuse cleanly */
                                printf("[hub] no memory for connection "
                                       "buffer; refusing\n");
                                close(s);
                                break;
                            }
                        }
                        conns[i].sock = s;
                        conns[i].len = 0;
                        printf("[hub] connection accepted\n");
                        break;
                    }
                if (i == MAXCONN) {
                    printf("[hub] all %d connection slots busy; "
                           "refusing\n", MAXCONN);
                    close(s);
                }
            }
        }

        for (i = 0; i < MAXCONN; i++) {
            conn_t *cn = &conns[i];
            ssize_t k;
            if (cn->sock < 0 || !FD_ISSET(cn->sock, &rd))
                continue;
            if (cn->len + 65536 > cn->cap) {
                unsigned char *nb = realloc(cn->buf, cn->cap * 2);
                if (!nb) { /* drop the connection, not the hub */
                    printf("[hub] out of buffer memory; dropping "
                           "connection\n");
                    release_conn(i);
                    continue;
                }
                cn->cap *= 2;
                cn->buf = nb;
            }
            k = recv(cn->sock, cn->buf + cn->len, 65536, 0);
            if (k <= 0) {
                printf("[hub] connection closed\n");
                release_conn(i);
                continue;
            }
            cn->len += (size_t)k;
            /* extract complete frames */
            for (;;) {
                unsigned camid, w, h;
                size_t need;
                cam_t *c;
                if (cn->len < 12)
                    break;
                /* MVLG: a pushed camera log line -- camid, len, text */
                if (memcmp(cn->buf, "MVLG", 4) == 0) {
                    unsigned lc = get32(cn->buf + 4);
                    unsigned ln = get32(cn->buf + 8);
                    if (ln > 4096) { /* bogus: resync */
                        unsigned char *m = memchr(cn->buf + 1, 'M',
                                                  cn->len - 1);
                        size_t sk = m ? (size_t)(m - cn->buf) : cn->len;
                        memmove(cn->buf, cn->buf + sk, cn->len - sk);
                        cn->len -= sk;
                        continue;
                    }
                    if (cn->len < 12 + (size_t)ln)
                        break;
                    camlog(lc, cn->buf + 12, (int)ln);
                    memmove(cn->buf, cn->buf + 12 + ln,
                            cn->len - 12 - ln);
                    cn->len -= 12 + ln;
                    continue;
                }
                if (cn->len < 32)
                    break;
                if (memcmp(cn->buf, "MVFR", 4) != 0) {
                    /* resync: jump to the next candidate 'M' in ONE
                     * move -- per-byte moves are quadratic and WiFi
                     * corruption plus a deep backlog turns that into
                     * a CPU DoS (found by nettest) */
                    unsigned char *m = memchr(cn->buf + 1, 'M',
                                              cn->len - 1);
                    size_t skip = m ? (size_t)(m - cn->buf) : cn->len;
                    memmove(cn->buf, cn->buf + skip, cn->len - skip);
                    cn->len -= skip;
                    continue;
                }
                camid = get32(cn->buf + 4);
                w = get32(cn->buf + 8);
                h = get32(cn->buf + 12);
                if (w < 64 || h < 64 || w > 4096 || h > 4096
                    || (double)w * h > 4.5e6) {
                    /* dims out of range (the pixel cap also bounds
                     * attacker-chosen decode and buffer cost); resync
                     * one byte in, then jump to the next 'M' */
                    unsigned char *m = memchr(cn->buf + 1, 'M',
                                              cn->len - 1);
                    size_t skip = m ? (size_t)(m - cn->buf) : cn->len;
                    memmove(cn->buf, cn->buf + skip, cn->len - skip);
                    cn->len -= skip;
                    continue;
                }
                need = 32 + (size_t)w * h;
                if (cn->len < need)
                    break;
                pthread_mutex_lock(&mbx);
                c = cam_slot(camid, i);
                if (c) {
                    /* manage ONLY latest here; work is decoder-owned and
                     * must never be freed/resized from this thread */
                    if (!c->latest || c->w != (int)w
                        || c->h != (int)h) {
                        unsigned char *nb = malloc((size_t)w * h);
                        if (nb) {
                            free(c->latest);
                            c->latest = nb;
                            c->w = (int)w;
                            c->h = (int)h;
                        }
                    }
                    if (c->latest && c->w == (int)w && c->h == (int)h) {
                        memcpy(c->latest, cn->buf + 32,
                               (size_t)w * h);
                        c->t_latest = now_mono();
                        c->fresh = 1;
                        c->frames_rx++;
                        if (c->frames_rx == 1)
                            hlog("[hub] >>> camera %u CONNECTED, "
                                 "frames flowing (%ux%u) <<<\n",
                                 camid, w, h);
                    }
                }
                pthread_mutex_unlock(&mbx);
                memmove(cn->buf, cn->buf + need, cn->len - need);
                cn->len -= need;
            }
        }

        /* ACK each camera ~1/s: send back what we have received and
         * decoded, so the sender knows the hub is really getting its
         * frames (not just that TCP accepted them). MVAK | rx | reads. */
        if (now_mono() - last_ack > 1.0) {
            last_ack = now_mono();
            pthread_mutex_lock(&mbx);
            for (i = 0; i < MAXCAMS; i++) {
                cam_t *c = &cams[i];
                unsigned char ak[12];
                if (!c->used || c->conn < 0
                    || conns[c->conn].sock < 0)
                    continue;
                ak[0] = 'M'; ak[1] = 'V'; ak[2] = 'A'; ak[3] = 'K';
                put32_le(ak + 4, (unsigned)c->frames_rx);
                put32_le(ak + 8, (unsigned)c->reads_ok);
                send(conns[c->conn].sock, ak, 12, 0);
            }
            pthread_mutex_unlock(&mbx);
        }

        if (now_mono() - last_status > 5.0) {
            last_status = now_mono();
            pthread_mutex_lock(&mbx);
            status();
            pthread_mutex_unlock(&mbx);
            if (reclog)
                fflush(reclog);
        }
    }

    /* clean shutdown: stop the decoder, close everything, free all
     * heap so Valgrind can distinguish real leaks from live state */
    hlog("\n[hub] shutting down\n");
    pthread_join(decoder_th, NULL);

    /* write a post-mortem summary into the recording dir: per-camera
     * throughput, decode rates, calibration, and the extrinsics -- the
     * one file to read to see what happened */
    if (recdir) {
        char path[512];
        FILE *f;
        snprintf(path, sizeof(path), "%s/summary.txt", recdir);
        f = fopen(path, "w");
        if (f) {
            int a, b;
            fprintf(f, "multiview live session summary\n");
            fprintf(f, "pitch %.4f mm\n\n", pitch * 1000.0);
            for (i = 0; i < MAXCAMS; i++) {
                cam_t *c = &cams[i];
                if (!c->used)
                    continue;
                fprintf(f, "camera %u: rx %ld frames, decoded %ld/%ld"
                        " (%.0f%%), valid-counter %ld, views %d\n",
                        c->camid, c->frames_rx, c->reads_ok, c->decodes,
                        c->decodes ? 100.0 * c->reads_ok / c->decodes : 0,
                        c->valid_ctr, c->nviews);
                if (c->calibrated)
                    fprintf(f, "  CALIBRATED fx %.1f fy %.1f cx %.1f "
                            "cy %.1f k1 %+.3f k2 %+.3f | RMS %.3f px\n",
                            c->K[0], c->K[4], c->K[2], c->K[5],
                            c->kr[0], c->kr[1], c->rms);
                else
                    fprintf(f, "  NOT calibrated (need ~12+ views with "
                            "pose diversity)\n");
            }
            fprintf(f, "\nextrinsics (camera pairs):\n");
            for (a = 0; a < MAXCAMS; a++)
                for (b = a + 1; b < MAXCAMS; b++)
                    if (cams[a].used && cams[b].used
                        && cams[a].calibrated && cams[b].calibrated)
                        fprintf(f, "  cam %u<->%u: %d anchors each; see"
                                " hub.log [ext] lines for the pose\n",
                                cams[a].camid, cams[b].camid,
                                cams[a].nanch < cams[b].nanch
                                    ? cams[a].nanch : cams[b].nanch);
            fclose(f);
            printf("[hub] wrote %s -- collection is in %s/ (frames, "
                   "records.txt, camlog_*.txt, hub.log, summary.txt)\n",
                   path, recdir);
        }
    }
    for (i = 0; i < MAXCONN; i++) {
        if (conns[i].sock >= 0)
            close(conns[i].sock);
        free(conns[i].buf);
    }
    for (i = 0; i < MAXCAMS; i++) {
        free(cams[i].latest);
        free(cams[i].work);
        free(cams[i].objs);
        free(cams[i].imgs);
    }
    close(lsock);
    if (reclog)
        fclose(reclog);
    if (hublog)
        fclose(hublog);
    return 0;
}
