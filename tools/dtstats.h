#ifndef DTSTATS_H
#define DTSTATS_H

/* Frame-arrival delta-t statistics: the measurement behind
 * ALIGNED_ASSESSMENT.md section 8.
 *
 * The section-6 phase argument says two free-running ~30 fps cameras
 * give nearly every frame a partner within half a frame period
 * (nearest-partner |dt| uniform on [0, T/2], mean T/4 ~ 8 ms).  This
 * module turns a recorded arrival log -- (camid, seq, t_arrival) per
 * frame, no pixels, no decode -- into the distribution that confirms
 * or kills that argument.  tools/dthub.c feeds it live; the unit test
 * feeds it synthetic streams.
 *
 * Pairing model: for each frame of the SPARSER of the two busiest
 * streams, the nearest-in-arrival-time frame of the other stream
 * (monotone scan, partners may repeat).  That is deliberately the
 * "does a simultaneous partner exist" question pair-then-decode asks,
 * not the one-to-one matching hub_pair.c later enforces.
 *
 * Clock: all timestamps must come from ONE clock (dthub uses the
 * hub host's CLOCK_MONOTONIC at header arrival).  Sender-side capture
 * stamps from two machines are NOT comparable here; that comparison
 * needs hub_clock's sync and is a different measurement. */

#include <stdio.h>

#define DT_MAXCAMS  4
#define DT_HBUCKETS 24   /* histogram buckets across [-T/2, +T/2) */

typedef struct dtstats dtstats;   /* opaque recorder */

typedef struct {
    unsigned camid;
    unsigned long long nframes;
    double t_first, t_last;   /* seconds, recorder clock */
    double fps;               /* (n-1)/(t_last-t_first) */
    double period_med;        /* median inter-arrival, s */
    double jitter_p90;        /* p90 of |inter-arrival - period_med| */
} dt_camsum;

typedef struct {
    int ncams;                    /* distinct camids seen */
    dt_camsum cam[DT_MAXCAMS];    /* per camera, ascending camid */

    /* pairing between the two busiest cameras */
    unsigned cam_a, cam_b;        /* cam_a is the lower camid */
    unsigned long long npairs;    /* one per frame of the sparser cam */
    double period;                /* min of the two median periods, s */
    double dt_mean, dt_median;    /* signed t(cam_a) - t(cam_b), s */
    double adt_p50, adt_p90, adt_p99, adt_max;   /* |dt| percentiles */
    double frac_half;             /* fraction with |dt| <  period/2 */
    double frac_quarter;          /* fraction with |dt| <= period/4 */

    /* signed-dt histogram over [-period/2, +period/2), plus tails */
    unsigned long hist[DT_HBUCKETS];
    unsigned long hist_lo, hist_hi;   /* dt < -T/2, dt >= +T/2 */
} dt_report;

dtstats *dt_new(void);
void dt_free(dtstats *s);

/* Record one frame arrival.  seq is kept only for CSV passthrough by
 * callers; analysis orders by t.  Frames must be added in
 * nondecreasing t per camera (a single-threaded receiver guarantees
 * this).  Returns 0, or -1 on allocation failure / too many cameras. */
int dt_add(dtstats *s, unsigned camid, unsigned long long seq, double t);

/* Analyze.  Needs >= 2 cameras with >= 8 frames each; returns 0 and
 * fills *r, or -1 if there is not enough data. */
int dt_analyze(const dtstats *s, dt_report *r);

/* Human-readable report: per-camera rates, |dt| percentiles against
 * the section-6 prediction, and an ASCII histogram. */
void dt_print(const dt_report *r, FILE *f);

#endif
