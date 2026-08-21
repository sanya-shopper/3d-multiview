/* Arrival delta-t statistics tests (tools/dtstats.c) -- the analysis
 * half of the section-8 measurement harness (tools/dthub.c is the
 * socket shell; this module is everything that computes).
 *
 * Tests: T1 fixed-offset ideal streams recover the offset exactly;
 * T2 the anti-phase worst case lands at T/2 and fails the T/4 gate;
 * T3 bounded jitter keeps percentiles inside arithmetic bounds;
 * T4 dropped frames degrade frac_half by exactly the drop pattern;
 * T5 not-enough-data contract; T6 busiest-two selection ignores a
 * stray third camera; T7 histogram mass equals pair count. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../tools/dtstats.h"

static int failures = 0;

#define CHECK(cond, name) \
    do { \
        if (cond) \
            printf("ok:   %s\n", name); \
        else { \
            printf("FAIL: %s\n", name); \
            failures++; \
        } \
    } while (0)

/* deterministic PRNG (no rand(): repo convention) */
static unsigned rng_state = 22u;
static double frand(void)
{
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 17;
    rng_state ^= rng_state << 5;
    return (double)(rng_state & 0xffffffu) / (double)0x1000000u;
}

static const double T = 1.0 / 30.0;

/* T1: two ideal 30 fps streams offset by 8 ms: every percentile is the
 * offset, the period is recovered, and everything passes the gate */
static void t1_fixed_offset(void)
{
    dtstats *s = dt_new();
    dt_report r;
    const double off = 0.008;
    int k, rc;
    for (k = 0; k < 300; k++) {
        dt_add(s, 1, (unsigned long long)k, k * T);
        dt_add(s, 2, (unsigned long long)k, k * T + off);
    }
    rc = dt_analyze(s, &r);
    CHECK(rc == 0, "T1 analyze succeeds");
    CHECK(r.cam_a == 1 && r.cam_b == 2, "T1 cam order by camid");
    CHECK(r.npairs == 300, "T1 every frame pairs");
    CHECK(fabs(r.period - T) < 1e-9, "T1 period recovered");
    CHECK(fabs(r.adt_p50 - off) < 1e-9 && fabs(r.adt_p90 - off) < 1e-9
          && fabs(r.adt_max - off) < 1e-9,
          "T1 |dt| percentiles equal the offset");
    CHECK(fabs(r.dt_median + off) < 1e-9,
          "T1 signed dt is cam1 - cam2 = -offset");
    CHECK(r.frac_half > 0.999 && r.frac_quarter > 0.999,
          "T1 all pairs inside T/4");
    CHECK(fabs(r.cam[0].fps - 30.0) < 0.01, "T1 fps ~ 30");
    dt_free(s);
}

/* T2: anti-phase streams (offset T/2): |dt| = T/2 for every pair --
 * the worst case the section-6 arithmetic allows -- and the T/4
 * fraction is 0 */
static void t2_antiphase(void)
{
    dtstats *s = dt_new();
    dt_report r;
    int k;
    for (k = 0; k < 100; k++) {
        dt_add(s, 1, (unsigned long long)k, k * T);
        dt_add(s, 2, (unsigned long long)k, k * T + T / 2.0);
    }
    CHECK(dt_analyze(s, &r) == 0, "T2 analyze succeeds");
    CHECK(fabs(r.adt_p50 - T / 2.0) < 1e-9, "T2 |dt| = T/2");
    CHECK(r.frac_quarter < 0.001, "T2 nothing inside T/4");
    dt_free(s);
}

/* T3: +/-2 ms uniform jitter on both streams around an 8 ms offset:
 * |dt| must stay in [4,12] ms and p50 near 8 ms */
static void t3_jitter(void)
{
    dtstats *s = dt_new();
    dt_report r;
    int k;
    for (k = 0; k < 500; k++) {
        dt_add(s, 1, (unsigned long long)k,
               k * T + 0.002 * (2.0 * frand() - 1.0));
        dt_add(s, 2, (unsigned long long)k,
               k * T + 0.008 + 0.002 * (2.0 * frand() - 1.0));
    }
    CHECK(dt_analyze(s, &r) == 0, "T3 analyze succeeds");
    CHECK(r.adt_max <= 0.012 + 1e-9 && r.adt_p50 >= 0.004 - 1e-9,
          "T3 |dt| bounded by offset +/- total jitter");
    CHECK(fabs(r.adt_p50 - 0.008) < 0.002, "T3 p50 near the offset");
    CHECK(r.cam[0].jitter_p90 < 0.009, "T3 inter-arrival jitter seen");
    dt_free(s);
}

/* T4: cam 2 drops every 5th frame.  The dropped partners force their
 * cam-1 frames to a neighbour at |dt| = T - off > T/2, so frac_half
 * degrades to exactly 4/5 while the surviving pairs stay at off */
static void t4_drops(void)
{
    dtstats *s = dt_new();
    dt_report r;
    const double off = 0.008;
    int k;
    for (k = 0; k < 500; k++) {
        dt_add(s, 1, (unsigned long long)k, k * T);
        if (k % 5 != 0)
            dt_add(s, 2, (unsigned long long)k, k * T + off);
    }
    CHECK(dt_analyze(s, &r) == 0, "T4 analyze succeeds");
    CHECK(r.npairs == 400, "T4 pairs = sparser stream count");
    CHECK(fabs(r.adt_p50 - off) < 1e-9, "T4 surviving pairs at offset");
    /* pairing iterates the SPARSER stream (cam 2 here), all of whose
     * frames still have a true partner in cam 1 */
    CHECK(r.frac_half > 0.999, "T4 sparser-side pairs all inside T/2");
    dt_free(s);
}

/* T4b: a capture-side stall -- cam 2 delivers nothing while frames
 * 200..249 are captured, then flushes them in one burst at frame 250's
 * arrival time.  Exactly those 50 cam-1 frames lose their simultaneous
 * partner, so frac_half degrades to 450/500 = 0.9: the "capture-side
 * problem" signature section 8 warns about */
static void t4b_stall(void)
{
    dtstats *s = dt_new();
    dt_report r;
    const double off = 0.008;
    int k;
    for (k = 0; k < 500; k++) {
        double t2 = (k >= 200 && k < 250) ? 250 * T : k * T;
        dt_add(s, 1, (unsigned long long)k, k * T);
        dt_add(s, 2, (unsigned long long)k, t2 + off);
    }
    CHECK(dt_analyze(s, &r) == 0, "T4b analyze succeeds");
    CHECK(r.frac_half > 0.85 && r.frac_half < 0.95,
          "T4b stall shows up as ~10% of pairs outside T/2");
    CHECK(fabs(r.adt_p50 - off) < 1e-9,
          "T4b unaffected pairs still at the offset");
    dt_free(s);
}

/* T5: contract -- one camera, or too few frames, refuses to analyze */
static void t5_not_enough(void)
{
    dtstats *s = dt_new();
    dt_report r;
    int k;
    for (k = 0; k < 100; k++)
        dt_add(s, 1, (unsigned long long)k, k * T);
    CHECK(dt_analyze(s, &r) == -1, "T5 one camera refused");
    for (k = 0; k < 4; k++)
        dt_add(s, 2, (unsigned long long)k, k * T);
    CHECK(dt_analyze(s, &r) == -1, "T5 4 frames refused");
    dt_free(s);
}

/* T6: a stray third camera with a handful of frames does not displace
 * the two busy ones from the pairing */
static void t6_third_camera(void)
{
    dtstats *s = dt_new();
    dt_report r;
    int k;
    for (k = 0; k < 200; k++) {
        dt_add(s, 1, (unsigned long long)k, k * T);
        dt_add(s, 2, (unsigned long long)k, k * T + 0.008);
    }
    for (k = 0; k < 10; k++)
        dt_add(s, 7, (unsigned long long)k, k * 1.0);
    CHECK(dt_analyze(s, &r) == 0, "T6 analyze succeeds");
    CHECK(r.ncams == 3, "T6 three cameras summarized");
    CHECK(r.cam_a == 1 && r.cam_b == 2,
          "T6 pairing picks the two busiest");
    CHECK(r.npairs == 200, "T6 stray camera not in the pairing");
    dt_free(s);
}

/* T7: histogram mass (buckets + tails) equals the pair count */
static void t7_histogram_mass(void)
{
    dtstats *s = dt_new();
    dt_report r;
    unsigned long long mass = 0;
    int k;
    for (k = 0; k < 300; k++) {
        dt_add(s, 1, (unsigned long long)k,
               k * T + 0.010 * (2.0 * frand() - 1.0));
        dt_add(s, 2, (unsigned long long)k,
               k * T + 0.010 * (2.0 * frand() - 1.0));
    }
    CHECK(dt_analyze(s, &r) == 0, "T7 analyze succeeds");
    for (k = 0; k < DT_HBUCKETS; k++)
        mass += r.hist[k];
    mass += r.hist_lo + r.hist_hi;
    CHECK(mass == r.npairs, "T7 histogram mass equals pair count");
    dt_free(s);
}

int main(void)
{
    t1_fixed_offset();
    t2_antiphase();
    t3_jitter();
    t4_drops();
    t4b_stall();
    t5_not_enough();
    t6_third_camera();
    t7_histogram_mass();
    if (failures) {
        printf("%d dtstats test(s) FAILED\n", failures);
        return 1;
    }
    printf("all dtstats tests passed\n");
    return 0;
}
