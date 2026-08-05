/* Live processing hub: receives camera streams over TCP, blind-decodes
 * the calibration pattern continuously, calibrates each camera as views
 * accumulate (robust Zhang + LM refinement, live), and estimates
 * camera-to-camera extrinsics from counter-matched pattern poses.
 * Portable POSIX C: runs on macOS and Ubuntu unchanged.
 *
 * DEPLOYMENT RUNBOOK (live test, three MacBooks):
 *  1. Pattern MacBook: copy tools/pattern.html over, open in a browser,
 *     click once -- it goes fullscreen in mux mode (near + far) by
 *     itself; no keys needed.  Inhibit sleep (macOS: caffeinate -dimsu).
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
 *  AUDIO: the load-bearing events (hub ready, camera connected /
 *  calibrated / stalled, pair solved, shutdown) are also spoken aloud
 *  (macOS `say`, Linux `spd-say`/`espeak`) -- during the session the
 *  laptop screens face the volume, not the operator.  MV_MUTE=1
 *  silences all speech.
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

#include "hub_clock.h"
#include "hub_solve.h"

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
    double t_cam_latest; /* sender-clock capture time of `latest` */
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
        double t_cam;  /* sender capture time (its own clock domain;
                        * meaningful across cameras only through
                        * hub_clock_map once clock sync is live) */
        mv_camera pose;
    } anch[MAXANCH];
    int nanch;         /* valid anchors, capped at MAXANCH */
    int anhead;        /* ring write index */
    int last_tier;
    unsigned last_counter;
    double last_dist;
    double t_ctr;      /* when the last VALID counter decoded (0 =
                        * never): drives the spoken guidance */
    int stall_said;    /* spoken-stall latch: announce each stall and
                        * each recovery once, not every status tick */
    int badcal_said;   /* implausible-calibration voice latch */
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

/* Spoken feedback: during a live session the laptops face each other
 * (screens are the pattern and the cameras' subjects), so the operator
 * cannot see any display -- announce the load-bearing events out loud.
 * Rare one-shot events only (connect / calibrate / pair / stall), so
 * the fork cost is irrelevant; the speech itself runs in background.
 * Text is program-authored (no network input reaches the shell).
 * MV_MUTE=1 disables. */
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
#ifdef __APPLE__
    snprintf(cmd, sizeof(cmd), "say \"%s\" >/dev/null 2>&1 &", msg);
#else
    snprintf(cmd, sizeof(cmd),
             "{ spd-say \"%s\" 2>/dev/null || espeak \"%s\"; } "
             ">/dev/null 2>&1 &", msg, msg);
#endif
    rc = system(cmd);
    (void)rc;
}

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
            cams[stale].stall_said = 0;
            cams[stale].t_ctr = 0.0;
            cams[stale].badcal_said = 0;
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
            int first, cw, ch;
            pthread_mutex_lock(&mbx);
            cw = c->w;
            ch = c->h;
            pthread_mutex_unlock(&mbx);
            /* sanity-gate the fit before it becomes the camera's K:
             * two live sessions produced degenerate 6-inlier refits
             * (fx/fy split 10%+, principal point off the sensor) that
             * the flag accepted and that would poison every anchor
             * pose.  Reject, keep the previous solution, say so. */
            if (!(K[0] > 0.0 && K[4] > 0.0
                  && K[0] / K[4] < 1.05 && K[4] / K[0] < 1.05
                  && K[2] > 0.2 * cw && K[2] < 0.8 * cw
                  && K[5] > 0.2 * ch && K[5] < 0.8 * ch)) {
                hlog("[cal] cam %u: REJECTED implausible fit fx %.1f "
                     "fy %.1f cx %.1f cy %.1f | RMS %.3f px over %d "
                     "views\n", c->camid, K[0], K[4], K[2], K[5],
                     rms, m);
                if (!c->badcal_said) {
                    c->badcal_said = 1;
                    say("camera %u calibration failed a sanity check."
                        " show it the pattern up close at more "
                        "angles", c->camid);
                }
                return;
            }
            c->badcal_said = 0;
            pthread_mutex_lock(&mbx);
            first = !c->calibrated;
            memcpy(c->K, K, sizeof(K));
            memcpy(c->kr, kr, sizeof(kr));
            c->rms = rms;
            c->calibrated = 1;
            pthread_mutex_unlock(&mbx);
            hlog("[cal] cam %u: fx %.1f fy %.1f cx %.1f cy %.1f "
                   "k1 %+.3f k2 %+.3f | RMS %.3f px over %d views\n",
                   c->camid, K[0], K[4], K[2], K[5], kr[0], kr[1],
                   rms, m);
            if (first)
                say("camera %u calibrated, r m s %.2f pixels",
                    c->camid, rms);
        }
    }
}

/* accumulate pairwise relative poses; report chordal mean */
/* pairing progress state, file-scope so the spoken guidance can see
 * which pairs are still unsolved (all guarded by mbx like the cams) */
static unsigned char pair_said[MAXCAMS][MAXCAMS];
static unsigned char pair_done[MAXCAMS][MAXCAMS];
static int npair_best[MAXCAMS][MAXCAMS];
static int npair_ann[MAXCAMS][MAXCAMS];
static double last_pose_say;

/* SCENE MODE: once the rig is solved the hub does not stop -- it
 * flushes the calibration artifacts and switches the decoder thread
 * to live stereo on the solved pair: detect features in a
 * near-simultaneous frame pair, match, epipolar-gate with the rig's
 * F, triangulate, and grow a metric point cloud written periodically
 * to <recdir>/scene.ply.  All rig and scene fields are mbx-guarded
 * except the cloud, which only the decoder thread touches. */
#define SCENE_MAXFEAT 800
#define SCENE_PAIR_DT 0.35   /* s: max capture-time skew of a pair */
#define SCENE_EPI_PX 3.0     /* px: symmetric epipolar gate */
#define SCENE_RMS_PX 2.0     /* px: triangulation reprojection gate */
#define SCENE_ZMIN 0.2       /* m, in camera-A frame */
#define SCENE_ZMAX 6.0
#define SCENE_MAXPTS 300000
#define SCENE_SETTLE_S 30.0  /* extra refinement time after solve */
static double rig_R9[9], rig_t3[3]; /* x_a = R x_b + t (solved pair) */
static int rig_a = -1, rig_b = -1;  /* cam table indices of the pair */
static int rig_have;
static double t_solved;
static int scene_mode;
static long scene_pairs, scene_npts;

static void try_extrinsics(void)
{
    int a, b;
    for (a = 0; a < MAXCAMS; a++)
        for (b = a + 1; b < MAXCAMS; b++) {
            cam_t *A = &cams[a], *B = &cams[b];
            static hub_obs obs[MAXPAIR];
            int npair = 0, i, j, k;
            int synced;
            if (!A->used || !B->used || !A->calibrated
                || !B->calibrated)
                continue;
            /* when the clock-sync module has both cameras mapped, gate
             * on mapped capture times; identity fallback = hub receive
             * times, i.e. the behavior before clock sync existed */
            synced = hub_clock_err(A->camid) >= 0.0
                     && hub_clock_err(B->camid) >= 0.0;
            for (i = 0; i < A->nanch && npair < MAXPAIR; i++) {
                for (j = 0; j < B->nanch && npair < MAXPAIR; j++) {
                    struct anchor *pa = &A->anch[i], *pb = &B->anch[j];
                    double dk, taT, tbT;
                    int relaxed = pa->dwell && pb->dwell;
                    if (synced) {
                        taT = hub_clock_map(A->camid, pa->t_cam);
                        tbT = hub_clock_map(B->camid, pb->t_cam);
                    }
                    else {
                        taT = pa->t_host;
                        tbT = pb->t_host;
                    }
                    /* both dwelling = display verified still: allow
                     * the slow decoder's 7-15 s anchor spacing to
                     * still produce cross-camera pairs */
                    if (fabs(taT - tbT) > (relaxed ? 15.0 : 4.0))
                        continue;
                    if (pa->tier == 1 && pb->tier == 1)
                        dk = fabs(pa->k - pb->k);
                    else {
                        double d8 = fmod(fmod(pa->k, 256.0)
                                         - fmod(pb->k, 256.0) + 384.0,
                                         256.0) - 128.0;
                        dk = fabs(d8);
                    }
                    if (dk > (relaxed ? 120.0 : 2.0))
                        continue;
                    memcpy(obs[npair].Ra, pa->pose.R,
                           sizeof(obs[0].Ra));
                    memcpy(obs[npair].ta, pa->pose.t,
                           sizeof(obs[0].ta));
                    memcpy(obs[npair].Rb, pb->pose.R,
                           sizeof(obs[0].Rb));
                    memcpy(obs[npair].tb, pb->pose.t,
                           sizeof(obs[0].tb));
                    npair++;
                }
            }
            (void)k;
            /* audible progress below the npair>=3 threshold: the
             * operator holding a dwell cannot see that matches are
             * (or are not) accumulating */
            if (npair > npair_best[a][b]) {
                npair_best[a][b] = npair;
                if (npair < 3 && !pair_said[a][b])
                    say("cameras %u and %u matched, keep holding",
                        A->camid, B->camid);
            }
            if (npair >= 3) {
                double Rm[9], tm[3], devmm;
                if (hub_solve_pair(Rm, tm, &devmm, obs, npair) != 0)
                    continue;
                hlog("[ext] cam %u <-> cam %u: baseline %.3f m  "
                       "t=(%+.3f %+.3f %+.3f)  %d pairs, mean |dt| "
                       "%.1f mm\n", A->camid, B->camid,
                       sqrt(tm[0] * tm[0] + tm[1] * tm[1]
                            + tm[2] * tm[2]), tm[0], tm[1], tm[2],
                       npair, devmm);
                /* latest solution becomes THE rig for scene mode */
                memcpy(rig_R9, Rm, sizeof(rig_R9));
                memcpy(rig_t3, tm, sizeof(rig_t3));
                rig_a = a;
                rig_b = b;
                rig_have = 1;
                /* the session goal, announced once per pair: after
                 * this the operator knows the rig is solved */
                if (!pair_said[a][b]) {
                    pair_said[a][b] = 1;
                    t_solved = now_mono();
                    npair_ann[a][b] = npair;
                    last_pose_say = now_mono();
                    say("cameras %u and %u paired, baseline %.2f "
                        "meters. change the angle and hold again",
                        A->camid, B->camid,
                        sqrt(tm[0] * tm[0] + tm[1] * tm[1]
                             + tm[2] * tm[2]));
                }
                /* dwell throughput feedback: a held pose typically
                 * banks a burst of matched pairs; once this pose has
                 * contributed, tell the operator to move on rather
                 * than leave them holding a solved angle */
                else if (!pair_done[a][b]
                         && npair - npair_ann[a][b] >= 4
                         && now_mono() - last_pose_say > 15.0) {
                    npair_ann[a][b] = npair;
                    last_pose_say = now_mono();
                    say("got this angle, next one please");
                }
                if (!pair_done[a][b] && npair >= MAXPAIR) {
                    pair_done[a][b] = 1;
                    say("cameras %u and %u have enough pairs, "
                        "pairing done", A->camid, B->camid);
                }
            }
        }
}

/* decode a camera's frame. img/w/h are a decoder-private snapshot (the
 * shared c->work/c->w/c->h must never be read here -- the network
 * thread can resize them concurrently); results are published into c
 * under mbx. */
static void process_cam(cam_t *c, const unsigned char *img, int w, int h,
                        double t_cam)
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
    if (rr.counter_valid) {
        c->valid_ctr++;
        c->t_ctr = now_mono();
    }
    {
        /* counter lock is the prerequisite for pairing (anchors need
         * the shared clock), and its absence is invisible without a
         * screen: announce the first one (live-session debug: cam 2
         * decoded corners for 40 min while an occluded counter strip
         * silently starved the extrinsics of anchors) */
        int first_ctr = rr.counter_valid && c->valid_ctr == 1;
        pthread_mutex_unlock(&mbx);
        if (first_ctr)
            say("camera %u reading the clock", c->camid);
    }

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
            an.t_cam = t_cam;
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
                /* 30 s window, not 5: at full-HD the decoder takes
                 * ~4 s per frame, so consecutive counter-valid
                 * anchors arrive 7-15 s apart and a 5 s rule can
                 * never fire (measured live 2026-08-05; the <1 cm
                 * pose gate still guarantees actual stillness) */
                an.dwell = sqrt(d) < 0.01
                           && an.t_host - pv->t_host < 30.0;
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

/* post-mortem / milestone summary into the recording dir: per-camera
 * throughput, decode rates, calibration, extrinsics, scene stats --
 * the one file to read to see what happened.  Called at scene-mode
 * entry (auto-flush: the calibration result is on disk from that
 * moment) and again at shutdown. */
static void write_summary(void)
{
    char path[512];
    FILE *f;
    int i, a, b;
    if (!recdir)
        return;
    snprintf(path, sizeof(path), "%s/summary.txt", recdir);
    f = fopen(path, "w");
    if (!f)
        return;
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
    if (rig_have) {
        fprintf(f, "\nsolved rig (x_a = R x_b + t, a=cam %u b=cam %u)"
                ":\n  R = [%+.6f %+.6f %+.6f; %+.6f %+.6f %+.6f; "
                "%+.6f %+.6f %+.6f]\n  t = (%+.4f %+.4f %+.4f) m  "
                "baseline %.4f m\n",
                cams[rig_a].camid, cams[rig_b].camid,
                rig_R9[0], rig_R9[1], rig_R9[2], rig_R9[3], rig_R9[4],
                rig_R9[5], rig_R9[6], rig_R9[7], rig_R9[8],
                rig_t3[0], rig_t3[1], rig_t3[2],
                sqrt(rig_t3[0] * rig_t3[0] + rig_t3[1] * rig_t3[1]
                     + rig_t3[2] * rig_t3[2]));
    }
    if (scene_mode)
        fprintf(f, "\nscene mode: %ld stereo pairs, %ld points -> "
                "scene.ply\n", scene_pairs, scene_npts);
    fclose(f);
    printf("[hub] wrote %s\n", path);
}

/* one live-stereo step on the solved pair: snapshot a near-
 * simultaneous frame pair, feature-match, epipolar-gate under the
 * solved rig, triangulate, accumulate the cloud.  Runs on the decoder
 * thread; the pattern reader is no longer called in scene mode. */
static void scene_step(void)
{
    static mv_feature fA[SCENE_MAXFEAT], fB[SCENE_MAXFEAT];
    static int idx2[SCENE_MAXFEAT];
    static mv_cloud cloud;
    static int cloud_init, first_said;
    static double last_ply;
    cam_t *A, *B;
    int wA, hA, wB, hB, nA, nB, nm, i, ok = 0;
    mv_camera camA, camB;
    double F[9];

    pthread_mutex_lock(&mbx);
    A = &cams[rig_a];
    B = &cams[rig_b];
    if (A->used && B->used && A->latest && B->latest
        && (A->fresh || B->fresh)
        && fabs(A->t_latest - B->t_latest) < SCENE_PAIR_DT) {
        size_t sA = (size_t)A->w * A->h, sB = (size_t)B->w * B->h;
        if (A->workcap < (int)sA) {
            unsigned char *nb = realloc(A->work, sA);
            if (nb) {
                A->work = nb;
                A->workcap = (int)sA;
            }
        }
        if (B->workcap < (int)sB) {
            unsigned char *nb = realloc(B->work, sB);
            if (nb) {
                B->work = nb;
                B->workcap = (int)sB;
            }
        }
        if (A->workcap >= (int)sA && B->workcap >= (int)sB) {
            memcpy(A->work, A->latest, sA);
            memcpy(B->work, B->latest, sB);
            A->fresh = B->fresh = 0;
            A->busy = B->busy = 1;
            wA = A->w;
            hA = A->h;
            wB = B->w;
            hB = B->h;
            ok = 1;
        }
    }
    pthread_mutex_unlock(&mbx);
    if (!ok) {
        struct timespec d = { 0, 50000000L };
        nanosleep(&d, NULL);
        return;
    }

    nA = mv_feat_detect(fA, SCENE_MAXFEAT, A->work, wA, hA);
    nB = mv_feat_detect(fB, SCENE_MAXFEAT, B->work, wB, hB);
    nm = 0;
    if (nA > 0 && nB > 0)
        nm = mv_feat_match(idx2, fA, nA, fB, nB, 0.56);
    if (nm > 0) {
        /* world = camera A frame: A at identity, B from the rig
         * (x_a = R x_b + t  =>  x_b = R^T x_a - R^T t) */
        double Rt[9];
        int k;
        memset(&camA, 0, sizeof(camA));
        memset(&camB, 0, sizeof(camB));
        memcpy(camA.K, A->K, sizeof(camA.K));
        mv_cam_set_identity_pose(&camA);
        memcpy(camB.K, B->K, sizeof(camB.K));
        mv_mat_transpose(Rt, rig_R9, 3, 3);
        memcpy(camB.R, Rt, sizeof(camB.R));
        mv_mat_mul(camB.t, Rt, rig_t3, 3, 3, 1);
        for (k = 0; k < 3; k++)
            camB.t[k] = -camB.t[k];
        mv_fundamental_from_cams(F, &camA, &camB);

        if (!cloud_init) {
            mv_cloud_init(&cloud, 1);
            cloud_init = 1;
        }
        for (i = 0; i < nA; i++) {
            double ua[2], ub[2], da[2], db[2], X[3], uv4[4], rms;
            const mv_camera *cc[2] = { &camA, &camB };
            unsigned char rgb[3];
            if (idx2[i] < 0 || cloud.n >= SCENE_MAXPTS)
                continue;
            da[0] = fA[i].u;
            da[1] = fA[i].v;
            db[0] = fB[idx2[i]].u;
            db[1] = fB[idx2[i]].v;
            undo_distort(ua, da, 1, A->K, A->kr);
            undo_distort(ub, db, 1, B->K, B->kr);
            if (mv_sym_epipolar_dist(F, ua, ub) > SCENE_EPI_PX)
                continue;
            uv4[0] = ua[0];
            uv4[1] = ua[1];
            uv4[2] = ub[0];
            uv4[3] = ub[1];
            if (mv_triangulate(X, cc, uv4, 2) != MV_OK)
                continue;
            if (!(X[2] > SCENE_ZMIN && X[2] < SCENE_ZMAX))
                continue;
            rms = mv_reproj_rms(cc, uv4, 2, X);
            if (!(rms < SCENE_RMS_PX))
                continue;
            rgb[0] = rgb[1] = rgb[2] =
                A->work[(int)(da[1] + 0.5) * wA + (int)(da[0] + 0.5)];
            mv_cloud_push(&cloud, X, rgb);
        }
    }

    pthread_mutex_lock(&mbx);
    scene_pairs++;
    scene_npts = cloud_init ? cloud.n : 0;
    A->busy = B->busy = 0;
    pthread_mutex_unlock(&mbx);

    if (!first_said && cloud_init && cloud.n > 300) {
        first_said = 1;
        say("stereo is live. %d points and counting", cloud.n);
    }
    if (recdir && cloud_init && cloud.n > 0
        && now_mono() - last_ply > 30.0) {
        char path[512];
        last_ply = now_mono();
        snprintf(path, sizeof(path), "%s/scene.ply", recdir);
        mv_cloud_write_ply(path, &cloud);
    }
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
        double dtc = 0.0;
        /* rig solved and settled -> auto-flush and switch this thread
         * to live stereo; the pattern has done its job */
        if (!scene_mode) {
            int enter;
            pthread_mutex_lock(&mbx);
            enter = rig_have && pair_said[rig_a][rig_b]
                    && now_mono() - t_solved > SCENE_SETTLE_S;
            pthread_mutex_unlock(&mbx);
            if (enter) {
                write_summary();
                if (reclog)
                    fflush(reclog);
                hlog("[hub] >>> RIG SOLVED -- entering SCENE MODE "
                     "(live stereo on cams %u and %u) <<<\n",
                     cams[rig_a].camid, cams[rig_b].camid);
                say("rig saved. scene mode. remove the pattern, "
                    "stereo is live");
                scene_mode = 1;
            }
        }
        if (scene_mode) {
            scene_step();
            continue;
        }
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
                    dtc = c->t_cam_latest;
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
            process_cam(pick, pick->work, dw, dh, dtc);
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

/* Spoken guidance: name the state the rig is WAITING in and the
 * operator action that unblocks it, every ~30 s while it persists --
 * the operator is eyes-free and otherwise just stands there guessing.
 * Called from status() with mbx held. */
static void guidance(void)
{
    static double last_say;
    static long last_pts;
    double now = now_mono();
    int i, a, b;
    int used = 0, ncal = 0, uncal = -1, noclk = -1, only = -1;
    int unpaired = 0;

    if (now - last_say < 30.0)
        return;
    if (scene_mode) {
        /* live stereo progress, only while it is actually growing */
        if (scene_npts > last_pts) {
            say("%ld points", scene_npts);
            last_pts = scene_npts;
            last_say = now;
        }
        return;
    }
    for (i = 0; i < MAXCAMS; i++) {
        cam_t *c = &cams[i];
        if (!c->used || c->frames_rx == 0)
            continue;
        used++;
        only = i;
        if (!c->calibrated) {
            if (uncal < 0)
                uncal = i;
        }
        else {
            ncal++;
            /* "reading the clock" = a valid counter in the last 20 s */
            if ((c->t_ctr == 0.0 || now - c->t_ctr > 20.0)
                && noclk < 0)
                noclk = i;
        }
    }
    for (a = 0; a < MAXCAMS; a++)
        for (b = a + 1; b < MAXCAMS; b++)
            if (cams[a].used && cams[b].used && cams[a].calibrated
                && cams[b].calibrated) {
                if (!pair_said[a][b])
                    unpaired = 1;
            }

    if (used == 0)
        return;
    if (uncal >= 0)
        say("camera %u still calibrating. show it the pattern at "
            "different tilts", cams[uncal].camid);
    else if (used == 1 && ncal == 1)
        say("camera %u calibrated. waiting for a second camera",
            cams[only].camid);
    else if (ncal >= 2 && noclk >= 0)
        say("calibration done. waiting for time alignment. camera %u "
            "is not reading the clock. face the pattern toward it, "
            "closer", cams[noclk].camid);
    else if (ncal >= 2 && unpaired)
        say("calibration done, clocks reading. waiting for time "
            "alignment. hold the pattern still where both cameras "
            "see it");
    else if (ncal >= 2) {
        /* every pair solved: make success unmissable (live lesson:
         * the one-shot 'paired' announcement was missed and the
         * operator kept standing by, thinking it was stuck) */
        static int solved_said;
        if (solved_said)
            return;
        solved_said = 1;
        say("the rig is solved. more angles only sharpen it. "
            "you can stop when you like");
    }
    else
        return; /* nothing pending: stay quiet, keep the timer */
    last_say = now;
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
        if (scene_mode && (i == rig_a || i == rig_b))
            hlog(" | SCENE pairs %ld pts %ld", scene_pairs,
                 scene_npts);
        /* loud flag when a camera that WAS flowing has gone quiet */
        if (c->frames_rx > 0 && age > 2.0) {
            hlog(" | *** STALLED: no frames for %.0f s ***", age);
            if (!c->stall_said) {
                c->stall_said = 1;
                say("camera %u stalled", c->camid);
            }
        }
        else {
            hlog(" | frame age %.1f s", age);
            if (c->stall_said) {
                c->stall_said = 0;
                say("camera %u back", c->camid);
            }
        }
        hlog("\n");
    }
    guidance();
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
            /* the ACTUAL pitch in use, not the reference constant: a
             * replay reconstructs metric scale from this field */
            mv_session_ver(reclog, MV_PAT_SPEC_VERSION, pitch * 1000.0,
                           "hub");
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
    say("hub ready on port %d", port);

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
                /* MVTS: clock-sync probe reply (see hub_clock.h) */
                if (memcmp(cn->buf, "MVTS", 4) == 0) {
                    if (cn->len < 24)
                        break;
                    hub_clock_on_msg(cn->buf + 4, now_mono());
                    memmove(cn->buf, cn->buf + 24, cn->len - 24);
                    cn->len -= 24;
                    continue;
                }
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
                        unsigned long long tb = 0;
                        double tcam;
                        int bi;
                        for (bi = 7; bi >= 0; bi--)
                            tb = (tb << 8) | cn->buf[24 + bi];
                        memcpy(&tcam, &tb, 8);
                        memcpy(c->latest, cn->buf + 32,
                               (size_t)w * h);
                        c->t_cam_latest = tcam;
                        c->t_latest = now_mono();
                        c->fresh = 1;
                        c->frames_rx++;
                        if (c->frames_rx == 1) {
                            hlog("[hub] >>> camera %u CONNECTED, "
                                 "frames flowing (%ux%u) <<<\n",
                                 camid, w, h);
                            say("camera %u connected", camid);
                        }
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
                hub_clock_probe(conns[c->conn].sock);
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
    {
        int nc = 0, ncal = 0;
        for (i = 0; i < MAXCAMS; i++)
            if (cams[i].used) {
                nc++;
                ncal += cams[i].calibrated != 0;
            }
        say("hub stopped. %d cameras, %d calibrated", nc, ncal);
    }

    write_summary();
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
