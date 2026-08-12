#ifndef HUB_PAIR_H
#define HUB_PAIR_H

/* Temporal frame-pair collector: the pair-then-decode contract of
 * ALIGNED_ASSESSMENT.md sections 6-7.
 *
 * OWNERSHIP (parallel build): this header, tools/hub_pair.c, and
 * tests/test_hub_pair.c belong to the PAIRING work item, together with
 * the livehub.c integration points (frame ring at ingest, pair-first
 * decode scheduling, direct obs emission).
 *
 * The problem this retires: the hub used to decode frames one camera
 * at a time and then match the resulting sparse anchors ACROSS cameras
 * by display counter, with a dwell-relaxed gate admitting up to 15 s of
 * skew -- observations that were never simultaneous fed hub_solve_pair
 * as if they were one rigid measurement.  This module inverts the
 * order: frames are matched on CAPTURE TIMESTAMPS at arrival rate,
 * cheaply, before any decode; only frames that already have a
 * simultaneous partner are worth decoding, and the display counter
 * returns to being an observed quantity (and a sanity check) rather
 * than the correspondence signal.
 *
 * Pure timestamp bookkeeping: frames are opaque integer handles owned
 * by the caller (the hub maps them to raw-frame ring slots).  Two
 * free-running 30 fps sensors have a nearest-partner |dt| that is
 * uniform over a half frame period -- mean T/4 ~ 8.3 ms -- so with a
 * gate of half a period essentially every frame has a true partner and
 * pair rate equals capture rate (assessment, "pair-then-decode").
 * Streams whose clocks never land within the gate produce ZERO pairs
 * by design: replay of unrelated sessions must report nothing rather
 * than manufacture a rig estimate. */

#define HP_MAXCAMS 4
#define HP_RING    16

typedef struct {
    int handle;     /* caller's opaque frame id, -1 = empty */
    double t;       /* capture time on the COMMON clock */
    int consumed;   /* already used in a pair */
} hp_slot;

typedef struct {
    double gate;    /* max |dt| for a pair, seconds */
    hp_slot ring[HP_MAXCAMS][HP_RING];
    int head[HP_MAXCAMS];        /* next write index per camera */
    long npush[HP_MAXCAMS];
    long npairs;                 /* pairs formed (all camera pairs) */
    double sum_dt;               /* sum of |dt| over formed pairs */
} hub_pair;

/* gate_s: maximum |dt| between the two frames of a pair.  For 30 fps
 * free-running cameras half a frame period (1/60 s) admits nearly every
 * frame while bounding simultaneity error at ~17 ms. */
void hub_pair_init(hub_pair *hp, double gate_s);

/* Register an arrived frame (cam in [0,HP_MAXCAMS), t on the common
 * clock).  Returns the handle EVICTED to make room (oldest slot,
 * consumed or not), or -1 if a free slot was used -- the caller
 * recycles the evicted frame buffer. */
int hub_pair_push(hub_pair *hp, int cam, int handle, double t);

/* Take the freshest unconsumed cross-camera pair between cama and camb
 * with |dt| <= gate: scans cama newest-first, matches each frame to its
 * nearest unconsumed camb frame.  On success writes the two handles and
 * the signed dt (tb - ta), marks both frames consumed, and returns 1;
 * returns 0 when no admissible pair exists. */
int hub_pair_take(hub_pair *hp, int cama, int camb,
                  int *ha, int *hb, double *dt);

/* pairing statistics: pairs formed and their mean |dt| (0 if none) --
 * the measured number the assessment's open experiment asks for */
long hub_pair_formed(const hub_pair *hp);
double hub_pair_mean_dt(const hub_pair *hp);

#endif /* HUB_PAIR_H */
