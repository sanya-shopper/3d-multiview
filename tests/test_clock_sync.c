/* Clock-sync estimator tests (tools/hub_clock.c).
 *
 * Links ONLY tools/hub_clock.o (no libmv), so everything here is
 * self-contained.  Synthetic MVTS wire bytes are constructed exactly
 * as a camera would build them (little-endian, raw IEEE-754 bits) and
 * fed to hub_clock_on_msg; ground truth is the known simulated
 * mapping hub_time = PHI + RATE * t_cam, so every accuracy assertion
 * is against independent truth and CAN fail (a broken fit, a missing
 * gate, or a missing restart-reset each trips a specific check).
 *
 * Determinism: fixed LCG, no rand()/time seeds.
 * OWNERSHIP (parallel build): clock-sync work item ONLY. */

#define _POSIX_C_SOURCE 200809L
#define _DARWIN_C_SOURCE 1

#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/socket.h>
#include <unistd.h>

#include "../tools/hub_clock.h"

static int g_fail;

#define CHECK(cond, msg) \
    do { \
        if (cond) { \
            printf("ok:   %s\n", msg); \
        } else { \
            printf("FAIL: %s\n", msg); \
            g_fail = 1; \
        } \
    } while (0)

/* ---- deterministic PRNG (fixed LCG; no wall clock anywhere) ------ */
static unsigned long g_lcg = 0x2545F491UL;

static double frand(void) /* uniform [0,1) */
{
    g_lcg = (g_lcg * 1664525UL + 1013904223UL) & 0xFFFFFFFFUL;
    return (double)g_lcg / 4294967296.0;
}

static double grand(void) /* ~N(0,1), Irwin-Hall of 12 uniforms */
{
    double s = 0.0;
    int i;
    for (i = 0; i < 12; i++)
        s += frand();
    return s - 6.0;
}

/* ---- wire construction, byte-for-byte as a camera builds it ------ */
static void put32_le(unsigned char *p, unsigned v)
{
    p[0] = (unsigned char)(v & 255);
    p[1] = (unsigned char)((v >> 8) & 255);
    p[2] = (unsigned char)((v >> 16) & 255);
    p[3] = (unsigned char)((v >> 24) & 255);
}

static void putf64_le(unsigned char *p, double d)
{
    unsigned long long b;
    int i;
    memcpy(&b, &d, 8);
    for (i = 0; i < 8; i++)
        p[i] = (unsigned char)((b >> (8 * i)) & 255);
}

/* deliver one MVTS as the hub's network loop would: body after magic */
static void deliver(unsigned camid, double t_echo, double t_cam,
                    double t_rx)
{
    unsigned char body[20];
    put32_le(body, camid);
    putf64_le(body + 4, t_echo);
    putf64_le(body + 12, t_cam);
    hub_clock_on_msg(body, t_rx);
}

/* ---- probe-exchange simulator ------------------------------------
 * Truth: hub_time = phi + rate * t_cam.  The hub probes at hub time
 * t_send; the forward leg takes frac*rtt, the camera replies with no
 * processing delay, the return leg takes (1-frac)*rtt.  The estimator
 * therefore sees a midpoint error of (0.5 - frac) * rtt. */
static void probe(unsigned camid, double phi, double rate,
                  double t_send, double rtt, double frac)
{
    double t_cam = (t_send + frac * rtt - phi) / rate;
    deliver(camid, t_send, t_cam, t_send + rtt);
}

/* jittered near-symmetric round trip: base 2 ms + |N(0,1)| * 2 ms,
 * clipped to the 2-8 ms band; split 30/70..70/30 */
static double good_rtt(void)
{
    double r = 0.002 + fabs(grand()) * 0.002;
    if (r > 0.008)
        r = 0.008;
    return r;
}

static double good_frac(void)
{
    return 0.3 + 0.4 * frand();
}

/* worst mapping error against truth over [t0,t1] in camera time */
static double span_err(unsigned camid, double phi, double rate,
                       double t0, double t1)
{
    double worst = 0.0;
    int i;
    for (i = 0; i <= 100; i++) {
        double t = t0 + (t1 - t0) * i / 100.0;
        double e = fabs(hub_clock_map(camid, t) - (phi + rate * t));
        if (e > worst)
            worst = e;
    }
    return worst;
}

/* ---- scenarios --------------------------------------------------- */

/* offset-only camera with jitter and gross+moderate outliers */
static void test_offset_and_outliers(void)
{
    /* camera clock = hub - 12345.678 s (hub = phi + 1.0 * t_cam);
     * probing starts at hub time 20000 so t_cam stays positive like a
     * real uptime */
    const double phi = 12345.678, rate = 1.0;
    const unsigned id = 2;
    double t;
    int i;

    CHECK(hub_clock_err(id) < 0.0, "no sync before any samples");

    for (i = 0, t = 20000.0; i < 60; i++, t += 1.0) {
        double u = frand();
        if (u < 0.10) {
            /* gross outlier: RTT 50-200 ms, one-sided (Wi-Fi burst);
             * must fall to the 25 ms absolute gate */
            probe(id, phi, rate, t, 0.050 + 0.150 * frand(), 0.9);
        } else if (u < 0.20) {
            /* moderate outlier: RTT 15-24 ms, heavily asymmetric;
             * inside the absolute gate, must fall to the 3x-median
             * gate (would bias the midpoint by ~+6..10 ms) */
            probe(id, phi, rate, t, 0.015 + 0.009 * frand(), 0.9);
        } else {
            probe(id, phi, rate, t, good_rtt(), good_frac());
        }
    }
    CHECK(hub_clock_err(id) >= 0.0, "sync established after samples");
    CHECK(hub_clock_err(id) < 0.006,
          "err estimate ~ median good RTT/2, not inflated by outliers");
    CHECK(span_err(id, phi, rate, 20000.0 - phi, 20060.0 - phi)
              < 0.0015,
          "known offset recovered to < 1.5 ms despite 20% outliers");
    printf("      measured: offset case worst error %.3f ms, "
           "err estimate %.3f ms\n",
           1e3 * span_err(id, phi, rate, 20000.0 - phi, 20060.0 - phi),
           1e3 * hub_clock_err(id));
}

/* the RTT gate is load-bearing: 70% gross one-sided outliers exceed
 * the breakdown point of ANY median-style fit, so only the RTT
 * admission gate can save the estimate.  The fixture first proves it
 * can fail (the ungated median offset is off by tens of ms), then
 * checks the gated result is still < 1.5 ms. */
static void test_gating_matters(void)
{
    const double phi = -777.25, rate = 1.0;
    const unsigned id = 3;
    double t, biased_med, resid[60];
    int i, n = 0;

    for (i = 0, t = 500.0; i < 60; i++, t += 1.0) {
        double rtt, frac;
        if (frand() < 0.7) { /* 70% congested, one-sided delays */
            rtt = 0.050 + 0.150 * frand();
            frac = 0.9;
        } else {
            rtt = good_rtt();
            frac = good_frac();
        }
        probe(id, phi, rate, t, rtt, frac);
        /* what an ungated estimator would see: midpoint offsets */
        {
            double t_cam = (t + frac * rtt - phi) / rate;
            resid[n++] = (t + 0.5 * rtt) - t_cam; /* offset sample */
        }
    }
    /* median of ALL midpoint offsets (the ungated fit's intercept);
     * insertion sort keeps this test free of library dependencies */
    for (i = 1; i < n; i++) {
        int j = i;
        double v = resid[i];
        while (j > 0 && resid[j - 1] > v) {
            resid[j] = resid[j - 1];
            j--;
        }
        resid[j] = v;
    }
    biased_med = resid[n / 2] - phi;
    CHECK(fabs(biased_med) > 0.005,
          "fixture can fail: ungated median is off by > 5 ms");
    CHECK(hub_clock_err(id) >= 0.0,
          "sync established from the fast minority");
    CHECK(span_err(id, phi, rate, 500.0 - phi, 560.0 - phi) < 0.0015,
          "RTT gate keeps < 1.5 ms under 70% outliers");
}

/* +80 ppm skew over a simulated 10-minute span */
static void test_skew_10min(void)
{
    const double phi = -3600.5, rate = 1.0 + 80e-6;
    const unsigned id = 4;
    double t;
    int i;

    /* one probe every ~9.4 s -> 64 samples cover the full 600 s and
     * all fit in the ring, so the fit sees the whole span */
    for (i = 0, t = 1000.0; i < 64; i++, t += 600.0 / 63.0) {
        double u = frand();
        if (u < 0.15)
            probe(id, phi, rate, t, 0.050 + 0.150 * frand(), 0.9);
        else
            probe(id, phi, rate, t, good_rtt(), good_frac());
    }
    CHECK(hub_clock_err(id) >= 0.0, "sync established with skew");
    CHECK(span_err(id, phi, rate, (1000.0 - phi) / rate,
                   (1600.0 - phi) / rate) < 0.0015,
          "+80 ppm skew: mapping < 1.5 ms across the full 10 min");
    printf("      measured: skew case worst error over 600 s: "
           "%.3f ms\n",
           1e3 * span_err(id, phi, rate, (1000.0 - phi) / rate,
                          (1600.0 - phi) / rate));
    /* the two ends specifically (48 ms of accumulated skew if the
     * rate term were ignored) */
    {
        double ta = (1000.0 - phi) / rate, tb = (1600.0 - phi) / rate;
        CHECK(fabs(hub_clock_map(id, ta) - (phi + rate * ta)) < 0.0015
              && fabs(hub_clock_map(id, tb) - (phi + rate * tb))
                 < 0.0015,
              "skew: both span endpoints < 1.5 ms");
    }
}

/* two cameras with wildly different uptime epochs, interleaved */
static void test_two_epochs(void)
{
    const double phiA = 4995.0, phiB = -1999999.5;
    const unsigned idA = 5, idB = 6;
    double t;
    int i;

    /* camera A booted ~seconds ago (t_cam ~ small), camera B has
     * ~23 days of uptime (t_cam ~ 2e6): same hub, same instant */
    for (i = 0, t = 5000.0; i < 30; i++, t += 1.0) {
        probe(idA, phiA, 1.0, t, good_rtt(), good_frac());
        probe(idB, phiB, 1.0, t, good_rtt(), good_frac());
    }
    CHECK(hub_clock_err(idA) >= 0.0 && hub_clock_err(idB) >= 0.0,
          "both epoch-divergent cameras synced");
    CHECK(span_err(idA, phiA, 1.0, 5000.0 - phiA, 5030.0 - phiA)
              < 0.0015
          && span_err(idB, phiB, 1.0, 5000.0 - phiB, 5030.0 - phiB)
              < 0.0015,
          "small-uptime and 23-day-uptime cameras both < 1.5 ms");
    /* cross-check: the same hub instant maps consistently */
    {
        double hub_a = hub_clock_map(idA, 5015.0 - phiA);
        double hub_b = hub_clock_map(idB, 5015.0 - phiB);
        CHECK(fabs(hub_a - hub_b) < 0.003,
              "cross-camera pairing agrees to < 3 ms");
    }
}

/* camera app restart: t_cam jumps backwards -> history must reset */
static void test_restart(void)
{
    const double phi1 = 4242.0, phi2 = 7019.5; /* new epoch: the
        camera rebooted just before hub time 7025, so its uptime
        restarts near zero (t_cam ~ 5.5 s, far below the old ~2780) */
    const unsigned id = 7;
    double t;
    int i;

    for (i = 0, t = 7000.0; i < 20; i++, t += 1.0)
        probe(id, phi1, 1.0, t, good_rtt(), good_frac());
    CHECK(hub_clock_err(id) >= 0.0, "synced before restart");

    /* restart: same camid reconnects, uptime back near zero, so the
     * new correspondence is hub = phi2 + t_cam with phi2 != phi1 */
    for (i = 0, t = 7025.0; i < 3; i++, t += 1.0)
        probe(id, phi2, 1.0, t, good_rtt(), good_frac());
    CHECK(hub_clock_err(id) < 0.0,
          "backwards t_cam wipes history (only 3 new samples)");
    for (i = 0, t = 7028.0; i < 17; i++, t += 1.0)
        probe(id, phi2, 1.0, t, good_rtt(), good_frac());
    CHECK(hub_clock_err(id) >= 0.0, "re-synced after restart");
    CHECK(span_err(id, phi2, 1.0, 7025.0 - phi2, 7045.0 - phi2)
              < 0.0015,
          "post-restart mapping uses the NEW epoch (< 1.5 ms)");
}

/* unknown camid: identity map, no error estimate */
static void test_unknown_cam(void)
{
    CHECK(hub_clock_err(42) < 0.0, "unknown camid: err < 0");
    CHECK(hub_clock_map(42, 123.456) == 123.456,
          "unknown camid: identity map");
}

/* malformed and hostile inputs must be ignored without disturbing an
 * established sync */
static void test_malformed(void)
{
    const double phi = 55.5;
    const unsigned id = 1;
    double t, before;
    int i;
    unsigned char junk[20];

    for (i = 0, t = 9000.0; i < 20; i++, t += 1.0)
        probe(id, phi, 1.0, t, good_rtt(), good_frac());
    CHECK(hub_clock_err(id) >= 0.0, "baseline sync for cam 1");
    before = hub_clock_map(id, 9010.0 - phi);

    deliver(id, NAN, 9020.0 - phi, 9020.0);        /* NaN echo */
    deliver(id, 9021.0, NAN, 9021.001);            /* NaN t_cam */
    deliver(id, INFINITY, 9022.0 - phi, 9022.0);   /* Inf echo */
    deliver(id, 9023.5, 9023.0 - phi, 9023.0);     /* negative RTT */
    deliver(id, 9024.0, 9024.0 - phi, 9024.0);     /* zero RTT */
    deliver(id, 8000.0, 9025.0 - phi, 9025.0);     /* RTT >> gate */
    hub_clock_on_msg(NULL, 9026.0);                /* NULL body */
    for (i = 0; i < 20; i++) /* garbage bytes as f64s */
        junk[i] = (unsigned char)(0xA5 ^ (i * 37));
    put32_le(junk, id);
    hub_clock_on_msg(junk, 9027.0);

    CHECK(hub_clock_err(id) >= 0.0
              && fabs(hub_clock_map(id, 9010.0 - phi) - before)
                 < 1e-9,
          "malformed inputs ignored; established fit undisturbed");
}

/* more camids than table slots: extras degrade to identity, tracked
 * cameras keep working */
static void test_table_full(void)
{
    double t;
    int i;
    /* ids 1..7 are in use by earlier scenarios; id 8 takes the last
     * slot, ids 9 and 10 must be turned away gracefully */
    for (i = 0, t = 11000.0; i < 10; i++, t += 1.0) {
        probe(8, 1.0, 1.0, t, good_rtt(), good_frac());
        probe(9, 2.0, 1.0, t, good_rtt(), good_frac());
        probe(10, 3.0, 1.0, t, good_rtt(), good_frac());
    }
    CHECK(hub_clock_err(8) >= 0.0, "8th camera fits in the table");
    CHECK(hub_clock_err(9) < 0.0 && hub_clock_err(10) < 0.0
              && hub_clock_map(9, 77.0) == 77.0,
          "camids beyond the table: identity fallback, no crash");
    CHECK(hub_clock_err(1) >= 0.0, "earlier cameras unaffected");
}

/* hub_clock_probe: writes exactly "MVPB"|f64 and never blocks */
static void test_probe_wire(void)
{
    int sv[2];
    unsigned char rx[64];
    ssize_t n;

    if (socketpair(AF_UNIX, SOCK_STREAM, 0, sv) != 0) {
        printf("FAIL: socketpair unavailable\n");
        g_fail = 1;
        return;
    }
    hub_clock_probe(sv[0]);
    n = recv(sv[1], rx, sizeof(rx), 0);
    CHECK(n == 12 && memcmp(rx, "MVPB", 4) == 0,
          "probe is exactly 12 bytes of MVPB");
    {
        unsigned long long b = 0;
        double d;
        int i;
        for (i = 7; i >= 0; i--)
            b = (b << 8) | rx[4 + i];
        memcpy(&d, &b, 8);
        CHECK(isfinite(d) && d > 0.0,
              "probe timestamp is a finite positive double");
    }
    /* fill the pipe: probe must return (skip), not block the caller;
     * reaching the CHECK at all is the assertion.  The fill loop sets
     * O_NONBLOCK itself because Darwin ignores MSG_DONTWAIT on
     * AF_UNIX send(); the probe is then called on the RESTORED
     * blocking socket, exactly the hub's situation. */
    {
        unsigned char blob[4096];
        int spins = 0, fl = fcntl(sv[0], F_GETFL, 0);
        memset(blob, 0, sizeof(blob));
        (void)fcntl(sv[0], F_SETFL, fl | O_NONBLOCK);
        while (send(sv[0], blob, sizeof(blob), 0) > 0
               && ++spins < 10000)
            ;
        (void)fcntl(sv[0], F_SETFL, fl);
        hub_clock_probe(sv[0]);
        CHECK(1, "probe on a full socket buffer returns immediately");
    }
    close(sv[0]);
    close(sv[1]);
    hub_clock_probe(-1); /* must not crash */
}

int main(void)
{
    test_probe_wire();
    test_unknown_cam();
    test_offset_and_outliers();
    test_gating_matters();
    test_skew_10min();
    test_two_epochs();
    test_restart();
    test_malformed();
    test_table_full();
    if (g_fail) {
        printf("CLOCK-SYNC TESTS FAILED\n");
        return 1;
    }
    printf("all clock-sync tests passed\n");
    return 0;
}
