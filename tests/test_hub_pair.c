/* Temporal frame-pair collector tests (hub_pair.c).
 * OWNERSHIP (parallel build): pairing work item ONLY.
 *
 * The collector pairs frames on capture timestamps before any decode
 * (ALIGNED_ASSESSMENT.md).  Tests: T1 exact offset streams; T2 the
 * mean nearest-partner |dt| ~ T/4 claim over uniform phases; T3 the
 * zero-pairs contract for unrelated clocks; T4 jitter + drops never
 * exceed the gate and never reuse a frame; T5 ring eviction hands
 * every displaced handle back exactly once; T6 determinism. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../tools/hub_pair.h"

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

/* T1: two 30 fps streams, fixed 8.3 ms phase offset, generous gate:
 * every frame pairs and every |dt| equals the offset */
static void t1_fixed_offset(void)
{
    hub_pair hp;
    const double T = 1.0 / 30.0, off = 0.0083;
    int k, ha, hb, pairs = 0, ok_dt = 1;
    double dt;
    hub_pair_init(&hp, T / 2.0);
    for (k = 0; k < 200; k++) {
        hub_pair_push(&hp, 0, 2 * k, k * T);
        hub_pair_push(&hp, 1, 2 * k + 1, k * T + off);
        while (hub_pair_take(&hp, 0, 1, &ha, &hb, &dt)) {
            pairs++;
            if (fabs(fabs(dt) - off) > 1e-12)
                ok_dt = 0;
        }
    }
    CHECK(pairs == 200, "T1: every frame of a fixed-offset stream pairs");
    CHECK(ok_dt, "T1: every pair's |dt| equals the true phase offset");
    CHECK(fabs(hub_pair_mean_dt(&hp) - off) < 1e-12,
          "T1: mean |dt| statistic equals the offset");
}

/* T2: uniform relative phase -> mean nearest-partner |dt| = T/4
 * (the assessment's ~8.3 ms at 30 fps) */
static void t2_uniform_phase(void)
{
    const double T = 1.0 / 30.0;
    double sum = 0.0;
    int nph = 400, p, k, ha, hb;
    long n = 0;
    for (p = 0; p < nph; p++) {
        hub_pair hp;
        double phase = (p + 0.5) / nph * T, dt;
        hub_pair_init(&hp, T / 2.0);
        for (k = 0; k < 40; k++) {
            hub_pair_push(&hp, 0, 2 * k, k * T);
            hub_pair_push(&hp, 1, 2 * k + 1, k * T + phase);
            while (hub_pair_take(&hp, 0, 1, &ha, &hb, &dt)) {
                sum += fabs(dt);
                n++;
            }
        }
    }
    /* nearest-partner |dt| for phase f is min(f, T-f): mean T/4 */
    CHECK(fabs(sum / (double)n - T / 4.0) < 0.02 * T,
          "T2: mean nearest-partner |dt| over uniform phase is T/4");
}

/* T3: unrelated clocks (replay of two different sessions) -> ZERO
 * pairs, not a misleading rig estimate */
static void t3_zero_pairs(void)
{
    hub_pair hp;
    int k, ha, hb;
    double dt;
    long pairs = 0;
    hub_pair_init(&hp, 0.008);
    for (k = 0; k < 100; k++) {
        hub_pair_push(&hp, 0, 2 * k, k / 30.0);
        hub_pair_push(&hp, 1, 2 * k + 1, 5000.0 + k / 30.0);
        while (hub_pair_take(&hp, 0, 1, &ha, &hb, &dt))
            pairs++;
    }
    CHECK(pairs == 0, "T3: clocks outside the gate form zero pairs");
    CHECK(hub_pair_formed(&hp) == 0, "T3: formed-pairs statistic is zero");
}

/* T4: jitter and 10% drops -- pairs never exceed the gate, no handle
 * is consumed twice, pair rate tracks the surviving frame rate */
static void t4_jitter_drops(void)
{
    hub_pair hp;
    const double T = 1.0 / 30.0;
    int k, ha, hb, used[4000];
    long pairs = 0, pushedb = 0;
    int reuse = 0, over = 0;
    double dt;
    memset(used, 0, sizeof(used));
    rng_state = 22u;
    hub_pair_init(&hp, T / 2.0);
    for (k = 0; k < 1000; k++) {
        double ja = (frand() - 0.5) * 0.006;
        double jb = (frand() - 0.5) * 0.006;
        hub_pair_push(&hp, 0, 2 * k, k * T + ja);
        if (frand() > 0.10) {
            hub_pair_push(&hp, 1, 2 * k + 1, k * T + 0.004 + jb);
            pushedb++;
        }
        while (hub_pair_take(&hp, 0, 1, &ha, &hb, &dt)) {
            pairs++;
            if (fabs(dt) > T / 2.0 + 1e-12)
                over = 1;
            if (used[ha] || used[hb])
                reuse = 1;
            used[ha] = used[hb] = 1;
        }
    }
    CHECK(!over, "T4: no pair exceeds the gate under jitter");
    CHECK(!reuse, "T4: no frame handle is consumed twice");
    CHECK(pairs >= pushedb - 20,
          "T4: pair rate tracks the surviving frame rate under drops");
}

/* T5: eviction accounting -- every pushed handle is exactly one of
 * (still stored | consumed in a pair | evicted back to the caller) */
static void t5_eviction(void)
{
    hub_pair hp;
    int k, ha, hb, ev;
    double dt;
    int state[300]; /* 0 = outstanding, 1 = consumed, 2 = evicted */
    int stored = 0, bad = 0, i;
    memset(state, 0, sizeof(state));
    hub_pair_init(&hp, 0.001); /* tight gate: nothing pairs */
    for (k = 0; k < 150; k++) {
        ev = hub_pair_push(&hp, 0, k, k * 1.0);
        if (ev >= 0) {
            if (state[ev] != 0)
                bad = 1; /* evicted something not outstanding */
            state[ev] = 2;
        }
    }
    (void)hub_pair_take(&hp, 0, 1, &ha, &hb, &dt);
    for (i = 0; i < 150; i++)
        if (state[i] == 0)
            stored++;
    CHECK(!bad, "T5: eviction never returns a handle twice");
    CHECK(stored == HP_RING,
          "T5: exactly one ring of handles remains outstanding");
}

/* T6: determinism -- identical inputs give identical pairings */
static void t6_determinism(void)
{
    long p1, p2;
    double m1, m2;
    int rep;
    for (rep = 0; rep < 2; rep++) {
        hub_pair hp;
        int k, ha, hb;
        double dt;
        rng_state = 77u;
        hub_pair_init(&hp, 1.0 / 60.0);
        for (k = 0; k < 500; k++) {
            hub_pair_push(&hp, 0, 2 * k, k / 30.0 + (frand() - 0.5) * 0.004);
            hub_pair_push(&hp, 1, 2 * k + 1, k / 30.0 + (frand() - 0.5) * 0.004);
            while (hub_pair_take(&hp, 0, 1, &ha, &hb, &dt))
                ;
        }
        if (rep == 0) {
            p1 = hub_pair_formed(&hp);
            m1 = hub_pair_mean_dt(&hp);
        }
        else {
            p2 = hub_pair_formed(&hp);
            m2 = hub_pair_mean_dt(&hp);
        }
    }
    CHECK(p1 == p2 && m1 == m2, "T6: identical inputs, identical pairing");
}

int main(void)
{
    t1_fixed_offset();
    t2_uniform_phase();
    t3_zero_pairs();
    t4_jitter_drops();
    t5_eviction();
    t6_determinism();
    if (failures) {
        printf("%d hub-pair test(s) FAILED\n", failures);
        return 1;
    }
    printf("all hub-pair tests passed\n");
    return 0;
}
