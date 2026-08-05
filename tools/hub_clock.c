/* Camera clock synchronization over the existing TCP connection.
 *
 * NTP-style: the hub sends "MVPB"|f64 t_hub_send probes (~1 Hz); a
 * camera replies "MVTS"|u32 camid|f64 t_hub_echo|f64 t_cam.  Each
 * reply yields RTT = t_rx - t_hub_echo and one correspondence
 *   t_cam  <->  (t_hub_echo + t_rx) / 2
 * (the camera stamped t_cam somewhere inside the round trip; the
 * midpoint is unbiased when the two legs are symmetric).  Per camera
 * we keep a ring of recent correspondences and robust-fit
 *   hub_time = phi + rate * t_cam
 * so capture timestamps can be mapped onto the hub clock to ~1-2 ms
 * instead of the ~100 ms jitter of hub receive times.
 *
 * Robustness gates (each one exists because live Wi-Fi produces the
 * misbehavior it guards against):
 *   - RTT > 25 ms or > 3x the ring's median RTT: sample excluded from
 *     the fit (asymmetric queueing delay poisons the midpoint).
 *   - fewer than 5 accepted samples: no fit, identity fallback.
 *   - fitted |rate - 1| > 1 %: garbage data, identity fallback (real
 *     crystal skew is tens of ppm).
 *   - t_cam moving backwards: the camera (machine) restarted, its old
 *     epoch is meaningless -> that camera's history is reset.
 * Fit = Theil-Sen slope + median intercept, then a MAD-gated least
 * squares refinement on the inliers.  Reported err = median accepted
 * RTT / 2 (the midpoint's worst-case asymmetry bound).
 *
 * CONCURRENCY: livehub calls hub_clock_on_msg() from the network
 * thread and hub_clock_map()/hub_clock_err() from the decoder thread;
 * all state here is guarded by the module mutex g_mx.  hub_clock_probe
 * touches no shared state and sends with MSG_DONTWAIT so the network
 * thread can never block on a camera's full socket buffer (a skipped
 * probe costs nothing at 1 Hz).
 *
 * DEPENDENCIES: deliberately none beyond libm/libpthread -- the
 * test_clock_sync binary links exactly this object and no libmv, so
 * the small robust fit lives here instead of calling mv_clock_fit,
 * and cmp_double below intentionally duplicates mv_cmp_double
 * (src/mat.c) for the same reason.
 *
 * OWNERSHIP (parallel build): clock-sync work item ONLY. */

/* glibc hides clock_gettime under -std=c99 -pedantic without this;
 * Darwin exposes it unconditionally */
#define _POSIX_C_SOURCE 200809L
#define _DARWIN_C_SOURCE 1

#include <fcntl.h>
#include <math.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>
#include <sys/socket.h>
#include <time.h>

#include "hub_clock.h"

/* the wire carries doubles as raw IEEE-754 bits through a u64 */
typedef char hc_assert_double8[sizeof(double) == 8 ? 1 : -1];
typedef char hc_assert_ull8[sizeof(unsigned long long) == 8 ? 1 : -1];

#define HC_CAMS 8            /* distinct camids tracked */
#define HC_RING 64           /* samples kept per camera */
#define HC_MINSAMP 5         /* accepted samples before a fit exists */
#define HC_RTT_ABS 0.025     /* s: absolute RTT admission gate */
#define HC_RTT_MEDX 3.0      /* fit-time gate: RTT <= this x median */
#define HC_RATE_TOL 0.01     /* sanity: fitted |rate-1| must be < 1 % */
#define HC_MAD_K 4.4478      /* 3 * 1.4826: 3-sigma MAD gate */
#define HC_MAD_FLOOR 5e-4    /* s: never gate tighter than 0.5 ms */

struct hc_sample {
    double t_cam;            /* camera CLOCK_MONOTONIC / systemUptime */
    double t_hub;            /* midpoint estimate on the hub clock */
    double rtt;              /* measured round trip, seconds */
};

struct hc_cam {
    int used;
    unsigned camid;
    struct hc_sample ring[HC_RING];
    int n;                   /* samples stored (<= HC_RING) */
    int head;                /* next write slot */
    double t_cam_last;       /* restart (backwards-jump) detection */
    int dirty;               /* ring changed since last fit */
    int have;                /* phi/rate/err below are valid */
    double phi, rate, err;
};

static struct hc_cam g_cams[HC_CAMS];
static pthread_mutex_t g_mx = PTHREAD_MUTEX_INITIALIZER;

static double now_mono(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

static unsigned get32_le(const unsigned char *p)
{
    return (unsigned)p[0] | ((unsigned)p[1] << 8)
           | ((unsigned)p[2] << 16) | ((unsigned)p[3] << 24);
}

static double getf64_le(const unsigned char *p)
{
    unsigned long long b = 0;
    double d;
    int i;
    for (i = 7; i >= 0; i--)
        b = (b << 8) | p[i];
    memcpy(&d, &b, 8);
    return d;
}

/* local copy of mv_cmp_double (src/mat.c): this object must link
 * without libmv (see header comment) */
static int cmp_double(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

/* median of v[0..n-1]; sorts v in place; 0.0 on an empty set (all
 * callers gate on HC_MINSAMP first, the guard is defensive) */
static double median_inplace(double *v, int n)
{
    if (n < 1)
        return 0.0;
    qsort(v, (size_t)n, sizeof(double), cmp_double);
    return (n & 1) ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

/* refit c from its ring; caller holds g_mx */
static void hc_refit(struct hc_cam *c)
{
    /* scratch guarded by g_mx (too big for comfort on a thread stack) */
    static double slopes[HC_RING * (HC_RING - 1) / 2];
    double rtts[HC_RING], tc[HC_RING], th[HC_RING], keep_rtt[HC_RING];
    double off[HC_RING], res[HC_RING];
    double med_rtt, rate, phi, gate;
    int i, j, m = 0, ns = 0;

    c->have = 0;
    if (c->n < HC_MINSAMP)
        return;

    for (i = 0; i < c->n; i++)
        rtts[i] = c->ring[i].rtt;
    med_rtt = median_inplace(rtts, c->n);

    /* fit-time gate: relative to the ring's own median RTT */
    for (i = 0; i < c->n; i++) {
        if (c->ring[i].rtt <= HC_RTT_MEDX * med_rtt) {
            tc[m] = c->ring[i].t_cam;
            th[m] = c->ring[i].t_hub;
            keep_rtt[m] = c->ring[i].rtt;
            m++;
        }
    }
    if (m < HC_MINSAMP)
        return;

    /* Theil-Sen: median of pairwise slopes */
    for (i = 0; i < m; i++)
        for (j = i + 1; j < m; j++) {
            double dt = tc[j] - tc[i];
            if (fabs(dt) > 1e-6)
                slopes[ns++] = (th[j] - th[i]) / dt;
        }
    if (ns < 1)
        return;
    rate = median_inplace(slopes, ns);
    if (!(fabs(rate - 1.0) < HC_RATE_TOL)) /* positive form: NaN fails */
        return;

    for (i = 0; i < m; i++) {
        off[i] = th[i] - rate * tc[i];
        res[i] = off[i];
    }
    phi = median_inplace(res, m);

    /* MAD gate around the robust line, then LS refine on inliers */
    for (i = 0; i < m; i++)
        res[i] = fabs(off[i] - phi);
    gate = HC_MAD_K * median_inplace(res, m);
    if (gate < HC_MAD_FLOOR)
        gate = HC_MAD_FLOOR;
    {
        double xm = 0.0, ym = 0.0;
        int k = 0;
        for (i = 0; i < m; i++)
            if (fabs(off[i] - phi) <= gate) {
                tc[k] = tc[i];
                th[k] = th[i];
                keep_rtt[k] = keep_rtt[i];
                k++;
            }
        if (k >= HC_MINSAMP) {
            double sxx = 0.0, sxy = 0.0;
            for (i = 0; i < k; i++) {
                xm += tc[i];
                ym += th[i];
            }
            xm /= k;
            ym /= k;
            for (i = 0; i < k; i++) {
                sxx += (tc[i] - xm) * (tc[i] - xm);
                sxy += (tc[i] - xm) * (th[i] - ym);
            }
            if (sxx > 0.0) {
                double r2 = sxy / sxx, p2 = ym - (sxy / sxx) * xm;
                if (isfinite(r2) && isfinite(p2)
                    && fabs(r2 - 1.0) < HC_RATE_TOL) {
                    rate = r2;
                    phi = p2;
                }
            }
            m = k; /* err reflects the samples actually used */
        }
    }

    c->rate = rate;
    c->phi = phi;
    c->err = 0.5 * median_inplace(keep_rtt, m);
    c->have = 1;
}

/* find (or allocate) the slot for camid; caller holds g_mx; NULL if
 * the table is full of other cameras (excess camids are ignored) */
static struct hc_cam *hc_slot(unsigned camid, int alloc)
{
    int i;
    for (i = 0; i < HC_CAMS; i++)
        if (g_cams[i].used && g_cams[i].camid == camid)
            return &g_cams[i];
    if (!alloc)
        return NULL;
    for (i = 0; i < HC_CAMS; i++)
        if (!g_cams[i].used) {
            memset(&g_cams[i], 0, sizeof(g_cams[i]));
            g_cams[i].used = 1;
            g_cams[i].camid = camid;
            return &g_cams[i];
        }
    return NULL;
}

void hub_clock_probe(int sock)
{
    unsigned char msg[12];
    unsigned long long bits;
    double t;
    int i, fl;
    if (sock < 0)
        return;
    t = now_mono();
    msg[0] = 'M'; msg[1] = 'V'; msg[2] = 'P'; msg[3] = 'B';
    memcpy(&bits, &t, 8);
    for (i = 0; i < 8; i++)
        msg[4 + i] = (unsigned char)((bits >> (8 * i)) & 255);
    /* non-blocking: a full socket buffer skips this probe rather than
     * stalling the hub's network thread.  Darwin does NOT honor
     * MSG_DONTWAIT on send() for every socket type (verified: AF_UNIX
     * blocks), so O_NONBLOCK is toggled around the send as well; this
     * is race-free because only livehub's network thread -- the
     * caller -- touches camera sockets.  A partial send can tear the
     * probe; camera-side scanners resync on the next magic, so the
     * cost is one lost probe. */
    fl = fcntl(sock, F_GETFL, 0);
    if (fl >= 0 && !(fl & O_NONBLOCK))
        (void)fcntl(sock, F_SETFL, fl | O_NONBLOCK);
    (void)send(sock, msg, sizeof(msg), MSG_DONTWAIT);
    if (fl >= 0 && !(fl & O_NONBLOCK))
        (void)fcntl(sock, F_SETFL, fl);
}

void hub_clock_on_msg(const unsigned char *body20, double t_rx)
{
    unsigned camid;
    double t_echo, t_cam, rtt;
    struct hc_cam *c;

    if (!body20)
        return;
    camid = get32_le(body20);
    t_echo = getf64_le(body20 + 4);
    t_cam = getf64_le(body20 + 12);
    if (!isfinite(t_echo) || !isfinite(t_cam) || !isfinite(t_rx))
        return;
    rtt = t_rx - t_echo;
    /* positive form so a NaN rtt is rejected, not admitted */
    if (!(rtt > 0.0) || !(rtt <= HC_RTT_ABS))
        return;

    pthread_mutex_lock(&g_mx);
    c = hc_slot(camid, 1);
    if (c) {
        if (c->n > 0 && t_cam < c->t_cam_last) {
            /* monotonic time went backwards: the camera restarted;
             * its previous epoch would poison the fit */
            c->n = 0;
            c->head = 0;
            c->have = 0;
        }
        c->t_cam_last = t_cam;
        c->ring[c->head].t_cam = t_cam;
        c->ring[c->head].t_hub = 0.5 * (t_echo + t_rx);
        c->ring[c->head].rtt = rtt;
        c->head = (c->head + 1) % HC_RING;
        if (c->n < HC_RING)
            c->n++;
        c->dirty = 1;
    }
    pthread_mutex_unlock(&g_mx);
}

double hub_clock_map(unsigned camid, double t)
{
    struct hc_cam *c;
    double out = t;
    pthread_mutex_lock(&g_mx);
    c = hc_slot(camid, 0);
    if (c) {
        if (c->dirty) {
            hc_refit(c);
            c->dirty = 0;
        }
        if (c->have)
            out = c->phi + c->rate * t;
    }
    pthread_mutex_unlock(&g_mx);
    return out;
}

double hub_clock_err(unsigned camid)
{
    struct hc_cam *c;
    double e = -1.0;
    pthread_mutex_lock(&g_mx);
    c = hc_slot(camid, 0);
    if (c) {
        if (c->dirty) {
            hc_refit(c);
            c->dirty = 0;
        }
        if (c->have)
            e = c->err;
    }
    pthread_mutex_unlock(&g_mx);
    return e;
}
