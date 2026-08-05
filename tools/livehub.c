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

#include <arpa/inet.h>
#include <errno.h>
#include <math.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/select.h>
#include <sys/socket.h>
#include <sys/stat.h>
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
    int w, h;
    unsigned char *latest;
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
    /* extrinsic anchors: display pose at decoded counter */
    struct anchor {
        double k;      /* counter (full for fine, mod 256 for coarse) */
        int tier;
        int dwell;     /* display static vs previous anchor */
        double t_host; /* hub receive time */
        mv_camera pose;
    } anch[MAXANCH];
    int nanch;
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
static double pitch;
static const char *recdir = NULL;
static FILE *reclog = NULL;

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

static cam_t *cam_slot(unsigned camid)
{
    int i;
    for (i = 0; i < MAXCAMS; i++)
        if (cams[i].used && cams[i].camid == camid)
            return &cams[i];
    for (i = 0; i < MAXCAMS; i++)
        if (!cams[i].used) {
            memset(&cams[i], 0, sizeof(cams[i]));
            cams[i].used = 1;
            cams[i].camid = camid;
            cams[i].objs = malloc((size_t)MAXV
                                  * sizeof(*cams[i].objs));
            cams[i].imgs = malloc((size_t)MAXV
                                  * sizeof(*cams[i].imgs));
            printf("[hub] camera %u joined\n", camid);
            return (cams[i].objs && cams[i].imgs) ? &cams[i] : NULL;
        }
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
            memcpy(c->K, K, sizeof(K));
            memcpy(c->kr, kr, sizeof(kr));
            c->rms = mv_calib_reproj_rms(ie, iv, m);
            c->calibrated = 1;
            printf("[cal] cam %u: fx %.1f fy %.1f cx %.1f cy %.1f "
                   "k1 %+.3f k2 %+.3f | RMS %.3f px over %d views\n",
                   c->camid, K[0], K[4], K[2], K[5], kr[0], kr[1],
                   c->rms, m);
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
                printf("[ext] cam %u <-> cam %u: baseline %.3f m  "
                       "t=(%+.3f %+.3f %+.3f)  %d pairs, mean |dt| "
                       "%.1f mm\n", A->camid, B->camid,
                       sqrt(tm[0] * tm[0] + tm[1] * tm[1]
                            + tm[2] * tm[2]), tm[0], tm[1], tm[2],
                       npair, 1000.0 * dev / npair);
            }
        }
}

/* decode the freshest frame of one camera; all downstream state flows
 * from here */
static void process_cam(cam_t *c)
{
    mv_read_result rr;
    int tier = 1, i;

    c->fresh = 0;
    c->decodes++;
    if (mv_read_pattern(&rr, c->latest, c->w, c->h) != MV_OK) {
        if (mv_read_coarse(&rr, c->latest, c->w, c->h) != MV_OK) {
            c->last_tier = 0;
            return;
        }
        tier = 2;
    }
    c->reads_ok++;
    c->last_tier = tier;
    c->last_counter = rr.counter_valid ? rr.counter : 0;
    if (rr.counter_valid)
        c->valid_ctr++;

    if (recdir && reclog) {
        char path[512], camname[32];
        snprintf(path, sizeof(path), "%s/cam%u_%06ld.pgm", recdir,
                 c->camid, c->reads_ok);
        snprintf(camname, sizeof(camname), "%u", c->camid);
        mv_pgm_write(path, c->latest, c->w, c->h);
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
            c->nviews++;
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
            c->last_dist = pose.R[6] * cx + pose.R[7] * cy + pose.t[2];
            if (c->nanch < MAXANCH) {
                struct anchor *an = &c->anch[c->nanch];
                an->k = (double)rr.counter;
                an->tier = tier;
                an->t_host = now_mono();
                an->pose = pose;
                an->dwell = 0;
                if (c->nanch > 0) {
                    struct anchor *pv = &c->anch[c->nanch - 1];
                    double d = 0.0;
                    for (i = 0; i < 3; i++)
                        d += (pv->pose.t[i] - pose.t[i])
                             * (pv->pose.t[i] - pose.t[i]);
                    an->dwell = sqrt(d) < 0.01
                                && an->t_host - pv->t_host < 5.0;
                }
                c->nanch++;
                try_extrinsics();
            }
        }
    }
}

static void status(void)
{
    int i;
    for (i = 0; i < MAXCAMS; i++) {
        cam_t *c = &cams[i];
        if (!c->used)
            continue;
        printf("[cam %u] rx %ld | decoded %ld/%ld | ctr %ld | views %d"
               " | %s", c->camid, c->frames_rx, c->reads_ok,
               c->decodes, c->valid_ctr, c->nviews,
               c->calibrated ? "CALIBRATED" : "collecting");
        if (c->calibrated)
            printf(" fx %.0f rms %.2f", c->K[0], c->rms);
        if (c->last_tier)
            printf(" | last: tier %d ctr %u dist %.2f m", c->last_tier,
                   c->last_counter, c->last_dist);
        printf("\n");
    }
}

int main(int argc, char **argv)
{
    int port, lsock, i;
    struct sockaddr_in addr;
    double last_status = 0.0;
    int rr_next = 0;

    if (argc < 3 || argc > 4) {
        fprintf(stderr, "usage: %s <port> <pitch_mm> [record_dir]\n",
                argv[0]);
        return 1;
    }
    port = atoi(argv[1]);
    pitch = atof(argv[2]) * 1e-3;
    if (argc == 4) {
        char path[512];
        recdir = argv[3];
        mkdir(recdir, 0755);
        snprintf(path, sizeof(path), "%s/records.txt", recdir);
        reclog = fopen(path, "w");
        if (reclog)
            mv_session_ver(reclog, MV_PAT_SPEC_VERSION, "hub");
    }
    signal(SIGPIPE, SIG_IGN);

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
    printf("[hub] listening on %d (pitch %.4f mm)%s%s\n", port,
           pitch * 1000.0, recdir ? ", recording to " : "",
           recdir ? recdir : "");

    for (;;) {
        fd_set rd;
        struct timeval tv;
        int maxfd = lsock, n;
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
        if (n < 0 && errno != EINTR)
            break;

        if (n > 0 && FD_ISSET(lsock, &rd)) {
            int s = accept(lsock, NULL, NULL);
            if (s >= 0) {
                for (i = 0; i < MAXCONN; i++)
                    if (conns[i].sock < 0) {
                        conns[i].sock = s;
                        conns[i].len = 0;
                        if (!conns[i].buf) {
                            conns[i].cap = 1 << 22;
                            conns[i].buf = malloc(conns[i].cap);
                        }
                        printf("[hub] connection accepted\n");
                        break;
                    }
                if (i == MAXCONN)
                    close(s);
            }
        }

        for (i = 0; i < MAXCONN; i++) {
            conn_t *cn = &conns[i];
            ssize_t k;
            if (cn->sock < 0 || !FD_ISSET(cn->sock, &rd))
                continue;
            if (cn->len + 65536 > cn->cap) {
                cn->cap *= 2;
                cn->buf = realloc(cn->buf, cn->cap);
            }
            k = recv(cn->sock, cn->buf + cn->len, 65536, 0);
            if (k <= 0) {
                printf("[hub] connection closed\n");
                close(cn->sock);
                cn->sock = -1;
                continue;
            }
            cn->len += (size_t)k;
            /* extract complete frames */
            for (;;) {
                unsigned camid, w, h;
                size_t need;
                cam_t *c;
                if (cn->len < 32)
                    break;
                if (memcmp(cn->buf, "MVFR", 4) != 0) {
                    /* resync: drop one byte */
                    memmove(cn->buf, cn->buf + 1, --cn->len);
                    continue;
                }
                camid = get32(cn->buf + 4);
                w = get32(cn->buf + 8);
                h = get32(cn->buf + 12);
                if (w < 64 || h < 64 || w > 4096 || h > 4096) {
                    memmove(cn->buf, cn->buf + 1, --cn->len);
                    continue;
                }
                need = 32 + (size_t)w * h;
                if (cn->len < need)
                    break;
                c = cam_slot(camid);
                if (c) {
                    if (!c->latest || c->w != (int)w
                        || c->h != (int)h) {
                        free(c->latest);
                        c->latest = malloc((size_t)w * h);
                        c->w = (int)w;
                        c->h = (int)h;
                    }
                    if (c->latest) {
                        memcpy(c->latest, cn->buf + 32,
                               (size_t)w * h);
                        c->t_latest = now_mono();
                        c->fresh = 1;
                        c->frames_rx++;
                    }
                }
                memmove(cn->buf, cn->buf + need, cn->len - need);
                cn->len -= need;
            }
        }

        /* decode ONE fresh frame per loop, round-robin: decode is much
         * slower than frame arrival, so always work on the newest */
        for (i = 0; i < MAXCAMS; i++) {
            cam_t *c = &cams[(rr_next + i) % MAXCAMS];
            if (c->used && c->fresh && c->latest) {
                process_cam(c);
                rr_next = ((rr_next + i) % MAXCAMS) + 1;
                break;
            }
        }

        if (now_mono() - last_status > 5.0) {
            last_status = now_mono();
            status();
            if (reclog)
                fflush(reclog);
        }
    }
    return 0;
}
