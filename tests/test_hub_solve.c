/* Rig pair-solver tests (hub_solve.c).
 * OWNERSHIP (parallel build): joint-solve work item ONLY.
 *
 * The solver estimates x_a = R x_b + t from matched display-pose
 * pairs.  Tests: T1 exactness on zero-noise synthetic rigs; T2 noise
 * accuracy vs the v0 chordal-mean + median reference; T3 gross random
 * outliers; T4 clustered outliers (the briefly-degenerate-calibration
 * failure mode); T5 degenerate inputs; T6 determinism. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mv/mv.h"

#include "../tools/hub_solve.h"

#define PI 3.14159265358979323846

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

/* ---- deterministic noise (project-standard LCG) ------------------ */

static unsigned long long rng_state;

static double urand(void)
{
    rng_state = rng_state * 6364136223846793005ULL
              + 1442695040888963407ULL;
    return (double)(rng_state >> 33) / 2147483648.0;
}

/* gaussian-ish: sum of 12 uniforms - 6 (unit variance) */
static double grand12(void)
{
    double g = -6.0;
    int i;
    for (i = 0; i < 12; i++)
        g += urand();
    return g;
}

/* ---- rotation helpers -------------------------------------------- */

static void rot_axis_deg(double R[9], double ax, double ay, double az,
                         double deg)
{
    double r[3];
    double nn, th;
    r[0] = ax;
    r[1] = ay;
    r[2] = az;
    nn = mv_norm(r, 3);
    th = deg * PI / 180.0;
    r[0] *= th / nn;
    r[1] *= th / nn;
    r[2] *= th / nn;
    mv_rot_exp(R, r);
}

/* random axis (gaussian direction), angle uniform in [mindeg, maxdeg] */
static void rand_rot(double R[9], double mindeg, double maxdeg)
{
    double r[3];
    double th;
    int k;
    do {
        for (k = 0; k < 3; k++)
            r[k] = grand12();
    } while (mv_norm(r, 3) < 1e-6);
    th = (mindeg + (maxdeg - mindeg) * urand()) * PI / 180.0;
    mv_normalize(r, 3);
    for (k = 0; k < 3; k++)
        r[k] *= th;
    mv_rot_exp(R, r);
}

/* left-perturb R by a gaussian axis-angle step, sigma per axis (rad) */
static void perturb_rot(double R[9], double sig_rad)
{
    double e[3], E[9], Rn[9];
    int k;
    for (k = 0; k < 3; k++)
        e[k] = sig_rad * grand12();
    mv_rot_exp(E, e);
    mv_mat_mul(Rn, E, R, 3, 3, 3);
    memcpy(R, Rn, sizeof(Rn));
}

static double rot_err_deg(const double A[9], const double B[9])
{
    double At[9], D[9], r[3];
    mv_mat_transpose(At, A, 3, 3);
    mv_mat_mul(D, At, B, 3, 3, 3);
    mv_rot_log(r, D);
    return mv_norm(r, 3) * 180.0 / PI;
}

static double t_err_mm(const double a[3], const double b[3])
{
    double d[3];
    int k;
    for (k = 0; k < 3; k++)
        d[k] = a[k] - b[k];
    return 1000.0 * mv_norm(d, 3);
}

/* ---- synthetic observations -------------------------------------- */

#define MAXOBS_T 300

static hub_obs g_obs[MAXOBS_T];

/* one matched pair: display pose in camera b sampled like the live
 * sessions (tilts to 40 deg, distances 0.3-0.7 m), camera a pose
 * composed exactly through the ground truth, then both sides
 * perturbed by sig_rot_deg (per axis) and sig_t_mm (per axis). */
static void make_obs(hub_obs *o, const double Rstar[9],
                     const double tstar[3], double sig_rot_deg,
                     double sig_t_mm)
{
    double tb[3], tt[3];
    int k;
    rand_rot(o->Rb, 0.0, 40.0);
    tb[0] = -0.25 + 0.5 * urand();
    tb[1] = -0.25 + 0.5 * urand();
    tb[2] = 0.3 + 0.4 * urand();
    memcpy(o->tb, tb, sizeof(tb));
    mv_mat_mul(o->Ra, Rstar, o->Rb, 3, 3, 3);
    mv_mat_mul(tt, Rstar, tb, 3, 3, 1);
    for (k = 0; k < 3; k++)
        o->ta[k] = tt[k] + tstar[k];
    if (sig_rot_deg > 0.0) {
        perturb_rot(o->Ra, sig_rot_deg * PI / 180.0);
        perturb_rot(o->Rb, sig_rot_deg * PI / 180.0);
    }
    if (sig_t_mm > 0.0)
        for (k = 0; k < 3; k++) {
            o->ta[k] += 0.001 * sig_t_mm * grand12();
            o->tb[k] += 0.001 * sig_t_mm * grand12();
        }
}

static void make_obs_set(hub_obs *obs, int n, const double Rstar[9],
                         const double tstar[3], double sig_rot_deg,
                         double sig_t_mm, unsigned long long seed)
{
    int i;
    rng_state = seed;
    for (i = 0; i < n; i++)
        make_obs(&obs[i], Rstar, tstar, sig_rot_deg, sig_t_mm);
}

/* gross random outlier: camera-a pose knocked 30-180 deg off with a
 * 0.2-1 m translation offset, direction random */
static void make_gross_outlier(hub_obs *o)
{
    double Rbad[9], Rn[9], dir[3];
    double mag;
    int k;
    rand_rot(Rbad, 30.0, 180.0);
    mv_mat_mul(Rn, Rbad, o->Ra, 3, 3, 3);
    memcpy(o->Ra, Rn, sizeof(Rn));
    do {
        for (k = 0; k < 3; k++)
            dir[k] = grand12();
    } while (mv_norm(dir, 3) < 1e-6);
    mv_normalize(dir, 3);
    mag = 0.2 + 0.8 * urand();
    for (k = 0; k < 3; k++)
        o->ta[k] += mag * dir[k];
}

/* clustered outlier: the briefly-degenerate-calibration mode -- every
 * outlier shares the same wrong transform (dR ~40 deg, dt ~0.3 m off
 * the truth) plus small jitter, forming a coherent second mode. */
static void make_clustered_outlier(hub_obs *o, const double Rstar[9],
                                   const double tstar[3])
{
    double dR[9], Rbad[9], tbad[3], tt[3];
    int k;
    rot_axis_deg(dR, 0.2, -1.0, 0.5, 40.0);
    mv_mat_mul(Rbad, dR, Rstar, 3, 3, 3);
    tbad[0] = tstar[0] + 0.25;
    tbad[1] = tstar[1] - 0.10;
    tbad[2] = tstar[2] + 0.15;
    mv_mat_mul(o->Ra, Rbad, o->Rb, 3, 3, 3);
    mv_mat_mul(tt, Rbad, o->tb, 3, 3, 1);
    for (k = 0; k < 3; k++)
        o->ta[k] = tt[k] + tbad[k];
    perturb_rot(o->Ra, 0.2 * PI / 180.0);
    for (k = 0; k < 3; k++)
        o->ta[k] += 0.002 * grand12();
}

/* ---- v0 reference: chordal mean + per-axis median, verbatim copy of
 * the pre-robust solver, kept here so the tests can demonstrate what
 * the joint estimate buys ------------------------------------------ */

#define REF_MAXOBS 256

static int ref_v0(double R[9], double t[3], double *dev_mm,
                  const hub_obs *obs, int n)
{
    double Rsum[9] = { 0 }, tacc[3][REF_MAXOBS];
    double U[9], S[3], V[9], Vt[9];
    double dev = 0.0;
    int i, k;

    if (n < 3)
        return -1;
    if (n > REF_MAXOBS)
        n = REF_MAXOBS;
    for (i = 0; i < n; i++) {
        double Rbt[9], Ri[9], ti[3];
        mv_mat_transpose(Rbt, obs[i].Rb, 3, 3);
        mv_mat_mul(Ri, obs[i].Ra, Rbt, 3, 3, 3);
        mv_mat_mul(ti, Ri, obs[i].tb, 3, 3, 1);
        for (k = 0; k < 3; k++)
            ti[k] = obs[i].ta[k] - ti[k];
        for (k = 0; k < 9; k++)
            Rsum[k] += Ri[k];
        for (k = 0; k < 3; k++)
            tacc[k][i] = ti[k];
    }
    memcpy(U, Rsum, sizeof(Rsum));
    if (mv_svd(U, S, V, 3, 3) != MV_OK)
        return -1;
    mv_mat_transpose(Vt, V, 3, 3);
    mv_mat_mul(R, U, Vt, 3, 3, 3);
    for (k = 0; k < 3; k++) {
        double col[REF_MAXOBS];
        int q, r;
        memcpy(col, tacc[k], (size_t)n * sizeof(double));
        for (q = 1; q < n; q++)
            for (r = q; r > 0 && col[r] < col[r - 1]; r--) {
                double tmp = col[r];
                col[r] = col[r - 1];
                col[r - 1] = tmp;
            }
        t[k] = col[n / 2];
    }
    for (i = 0; i < n; i++) {
        double d = 0.0;
        for (k = 0; k < 3; k++)
            d += (tacc[k][i] - t[k]) * (tacc[k][i] - t[k]);
        dev += sqrt(d);
    }
    if (dev_mm)
        *dev_mm = 1000.0 * dev / n;
    return 0;
}

/* ---- ground-truth rigs ------------------------------------------- */

enum { NGT = 4 };

static void gt_case(int c, double Rstar[9], double tstar[3])
{
    switch (c) {
    case 0: /* identity rotation */
        rot_axis_deg(Rstar, 1, 0, 0, 0.0);
        tstar[0] = 0.10;
        tstar[1] = -0.05;
        tstar[2] = 0.02;
        break;
    case 1: /* the real rig: 56 deg about a skew axis */
        rot_axis_deg(Rstar, 1.0, 2.0, -1.0, 56.0);
        tstar[0] = 0.45;
        tstar[1] = 0.03;
        tstar[2] = -0.12;
        break;
    case 2: /* near-pi, exercises the theta ~ pi log branch */
        rot_axis_deg(Rstar, 0.3, -1.0, 0.2, 179.0);
        tstar[0] = 0.20;
        tstar[1] = 0.20;
        tstar[2] = 0.20;
        break;
    default: /* small angle */
        rot_axis_deg(Rstar, 0.0, 0.0, 1.0, 3.0);
        tstar[0] = -0.02;
        tstar[1] = 0.50;
        tstar[2] = 0.01;
        break;
    }
}

/* ---- T1: exactness ----------------------------------------------- */

static void test_exact(void)
{
    static const int ns[] = { 3, 4, 5, 10, 37, 70 };
    double worst_r = 0.0, worst_t = 0.0;
    int c, u;

    for (c = 0; c < NGT; c++)
        for (u = 0; u < (int)(sizeof(ns) / sizeof(ns[0])); u++) {
            double Rstar[9], tstar[3], R[9], t[3], dev;
            double er, et;
            gt_case(c, Rstar, tstar);
            make_obs_set(g_obs, ns[u], Rstar, tstar, 0.0, 0.0,
                         0xC0FFEEULL + (unsigned long long)(c * 100 + u));
            if (hub_solve_pair(R, t, &dev, g_obs, ns[u]) != 0) {
                CHECK(0, "exact: solver returned success");
                return;
            }
            er = rot_err_deg(R, Rstar) * PI / 180.0;
            et = t_err_mm(t, tstar) / 1000.0;
            if (er > worst_r)
                worst_r = er;
            if (et > worst_t)
                worst_t = et;
        }
    printf("      exact: worst rot err %.3e rad, worst t err %.3e m\n",
           worst_r, worst_t);
    CHECK(worst_r < 1e-9, "exact: rotation recovered to 1e-9 rad");
    CHECK(worst_t < 1e-9, "exact: translation recovered to 1e-9 m");
}

/* ---- T2: noise only ---------------------------------------------- */

static void test_noise(void)
{
    /* {n, sigma_rot_deg, sigma_t_mm} spanning the live sessions,
     * including the small-n regime where robustness must degrade
     * gracefully toward v0 */
    static const double cfg[][3] = {
        { 3, 0.2, 1.0 },
        { 5, 0.3, 3.0 },
        { 20, 0.5, 5.0 },
        { 70, 1.0, 10.0 }
    };
    enum { NSEED = 8 };
    int c;

    for (c = 0; c < (int)(sizeof(cfg) / sizeof(cfg[0])); c++) {
        int n = (int)cfg[c][0];
        double sr = cfg[c][1], st = cfg[c][2];
        double rob_r = 0.0, rob_t = 0.0, v0_r = 0.0, v0_t = 0.0;
        double Rstar[9], tstar[3];
        char name[128];
        int s;

        gt_case(1, Rstar, tstar);
        for (s = 0; s < NSEED; s++) {
            double R[9], t[3], Rr[9], tr[3], dev;
            make_obs_set(g_obs, n, Rstar, tstar, sr, st,
                         0xB00B5ULL + (unsigned long long)(c * 977 + s));
            if (hub_solve_pair(R, t, &dev, g_obs, n) != 0
                || ref_v0(Rr, tr, &dev, g_obs, n) != 0) {
                CHECK(0, "noise: both solvers succeed");
                return;
            }
            rob_r += rot_err_deg(R, Rstar);
            rob_t += t_err_mm(t, tstar);
            v0_r += rot_err_deg(Rr, Rstar);
            v0_t += t_err_mm(tr, tstar);
        }
        rob_r /= NSEED;
        rob_t /= NSEED;
        v0_r /= NSEED;
        v0_t /= NSEED;
        printf("      noise n=%2d sig=(%.1f deg,%4.1f mm): "
               "robust (%.4f deg, %.3f mm)  v0 (%.4f deg, %.3f mm)\n",
               n, sr, st, rob_r, rob_t, v0_r, v0_t);
        snprintf(name, sizeof(name), "noise n=%d: robust at least as good as v0", n);
        CHECK(rob_r <= v0_r * 1.02 + 1e-12 && rob_t <= v0_t * 1.02 + 1e-9,
              name);
        snprintf(name, sizeof(name), "noise n=%d: absolute accuracy sensible", n);
        /* mean-of-n scaling: per-axis sigma sr deg -> angle error of
         * the mean ~ sr*sqrt(3/n); allow 3x for finite-sample spread */
        CHECK(rob_r < 3.0 * sr * sqrt(3.0 / n)
              && rob_t < 3.0 * st * sqrt(3.0 / n) + 0.5, name);
    }
}

/* ---- T3: gross random outliers ----------------------------------- */

static void test_outliers(void)
{
    static const int cfg[][2] = { { 10, 3 }, { 50, 13 } }; /* 30%, 26% */
    int c;

    for (c = 0; c < 2; c++) {
        int n = cfg[c][0], nout = cfg[c][1];
        double Rstar[9], tstar[3], R[9], t[3], Rr[9], tr[3], dev;
        static hub_solve_diag diag[MAXOBS_T];
        double win = 0.0, wout = 0.0;
        char name[128];
        int i;

        gt_case(1, Rstar, tstar);
        make_obs_set(g_obs, n, Rstar, tstar, 0.2, 2.0,
                     0xDEADULL + (unsigned long long)c);
        /* corrupt the first nout pairs (position must not matter) */
        for (i = 0; i < nout; i++)
            make_gross_outlier(&g_obs[i]);

        CHECK(hub_solve_pair_diag(R, t, &dev, diag, g_obs, n) == 0,
              "outliers: robust solve succeeds");
        CHECK(ref_v0(Rr, tr, &dev, g_obs, n) == 0,
              "outliers: v0 solve succeeds");
        printf("      gross outliers n=%2d (%d bad): "
               "robust (%.4f deg, %.3f mm)  v0 (%.3f deg, %.2f mm)\n",
               n, nout, rot_err_deg(R, Rstar), t_err_mm(t, tstar),
               rot_err_deg(Rr, Rstar), t_err_mm(tr, tstar));
        snprintf(name, sizeof(name), "outliers n=%d: robust R within 0.3 deg", n);
        CHECK(rot_err_deg(R, Rstar) < 0.3, name);
        snprintf(name, sizeof(name), "outliers n=%d: robust t within 5 mm", n);
        CHECK(t_err_mm(t, tstar) < 5.0, name);
        snprintf(name, sizeof(name), "outliers n=%d: v0 visibly degraded", n);
        CHECK(rot_err_deg(Rr, Rstar) > 5.0 * rot_err_deg(R, Rstar)
              && rot_err_deg(Rr, Rstar) > 1.0, name);
        for (i = 0; i < n; i++) {
            if (i < nout)
                wout += diag[i].weight;
            else
                win += diag[i].weight;
        }
        win /= n - nout;
        wout /= nout;
        printf("      diag weights: inlier mean %.3f, outlier mean %.4f\n",
               win, wout);
        snprintf(name, sizeof(name), "outliers n=%d: diag downweights the bad pairs", n);
        CHECK(wout < 0.1 * win && win > 0.5, name);
    }
}

/* ---- T4: clustered outliers -------------------------------------- */

/* The briefly-degenerate-calibration mode: 20-30% of pairs share ONE
 * consistent wrong transform, so the residual population is bimodal
 * and a badly-seeded Huber fit can converge onto the wrong mode.  The
 * solver's documented behavior: score both the v0 mean and the medoid
 * candidate by trimmed cost and seed from the winner, which is always
 * the majority mode; IRLS then stays there.  These checks pin that
 * behavior down at n=5/10/50. */
static void test_clustered(void)
{
    static const int cfg[][2] = { { 5, 1 }, { 10, 3 }, { 50, 14 } };
    int c;

    for (c = 0; c < 3; c++) {
        int n = cfg[c][0], nout = cfg[c][1];
        double Rstar[9], tstar[3], R[9], t[3], Rr[9], tr[3], dev;
        char name[128];
        int i;

        gt_case(1, Rstar, tstar);
        make_obs_set(g_obs, n, Rstar, tstar, 0.2, 2.0,
                     0xFACEULL + (unsigned long long)c);
        /* corrupt a spread of indices, not just a prefix */
        for (i = 0; i < nout; i++)
            make_clustered_outlier(&g_obs[(i * n) / nout + n / (2 * nout)],
                                   Rstar, tstar);

        CHECK(hub_solve_pair(R, t, &dev, g_obs, n) == 0,
              "clustered: robust solve succeeds");
        CHECK(ref_v0(Rr, tr, &dev, g_obs, n) == 0,
              "clustered: v0 solve succeeds");
        printf("      clustered n=%2d (%d bad): "
               "robust (%.4f deg, %.3f mm)  v0 (%.3f deg, %.2f mm)\n",
               n, nout, rot_err_deg(R, Rstar), t_err_mm(t, tstar),
               rot_err_deg(Rr, Rstar), t_err_mm(tr, tstar));
        snprintf(name, sizeof(name), "clustered n=%d: robust R within 0.3 deg", n);
        CHECK(rot_err_deg(R, Rstar) < 0.3, name);
        snprintf(name, sizeof(name), "clustered n=%d: robust t within 5 mm", n);
        CHECK(t_err_mm(t, tstar) < 5.0, name);
        snprintf(name, sizeof(name), "clustered n=%d: robust beats v0 clearly", n);
        CHECK(rot_err_deg(Rr, Rstar) > 5.0 * rot_err_deg(R, Rstar), name);
    }
}

/* ---- T5: degenerate inputs --------------------------------------- */

static void test_degenerate(void)
{
    double Rstar[9], tstar[3], R[9], t[3], R2[9], t2[3], dev, dev2;
    int i;

    gt_case(1, Rstar, tstar);

    /* n < 3 rejected */
    make_obs_set(g_obs, 2, Rstar, tstar, 0.0, 0.0, 7ULL);
    CHECK(hub_solve_pair(R, t, &dev, g_obs, 2) == -1,
          "degenerate: n=2 returns -1");
    CHECK(hub_solve_pair(R, t, &dev, g_obs, 0) == -1,
          "degenerate: n=0 returns -1");

    /* all-identical observations: exact recovery, zero scatter */
    make_obs_set(g_obs, 1, Rstar, tstar, 0.0, 0.0, 11ULL);
    for (i = 1; i < 8; i++)
        g_obs[i] = g_obs[0];
    CHECK(hub_solve_pair(R, t, &dev, g_obs, 8) == 0,
          "degenerate: identical obs solve");
    CHECK(rot_err_deg(R, Rstar) * PI / 180.0 < 1e-9
          && t_err_mm(t, tstar) < 1e-6 && dev < 1e-6,
          "degenerate: identical obs recover exactly");

    /* n past the internal cap (256): entries beyond the cap must never
     * be read -- poison them with NaN so a read would show up */
    make_obs_set(g_obs, 256, Rstar, tstar, 0.5, 5.0, 13ULL);
    for (i = 256; i < 300; i++) {
        double qnan = nan("");
        int k;
        for (k = 0; k < 9; k++) {
            g_obs[i].Ra[k] = qnan;
            g_obs[i].Rb[k] = qnan;
        }
        for (k = 0; k < 3; k++) {
            g_obs[i].ta[k] = qnan;
            g_obs[i].tb[k] = qnan;
        }
    }
    CHECK(hub_solve_pair(R, t, &dev, g_obs, 300) == 0,
          "degenerate: n=300 solves (capped)");
    CHECK(hub_solve_pair(R2, t2, &dev2, g_obs, 256) == 0,
          "degenerate: n=256 solves");
    CHECK(memcmp(R, R2, sizeof(R)) == 0 && memcmp(t, t2, sizeof(t)) == 0
          && dev == dev2,
          "degenerate: n=300 identical to n=256 (cap, no NaN read)");
}

/* ---- T6: determinism --------------------------------------------- */

static void test_determinism(void)
{
    double Rstar[9], tstar[3];
    double Ra[9], ta[3], deva, Rb[9], tb[3], devb;
    static hub_solve_diag da[MAXOBS_T], db[MAXOBS_T];
    int i;

    gt_case(2, Rstar, tstar);
    make_obs_set(g_obs, 50, Rstar, tstar, 0.5, 5.0, 0xABCDULL);
    for (i = 0; i < 12; i++)
        make_gross_outlier(&g_obs[i * 4]);

    CHECK(hub_solve_pair_diag(Ra, ta, &deva, da, g_obs, 50) == 0
          && hub_solve_pair_diag(Rb, tb, &devb, db, g_obs, 50) == 0,
          "determinism: both runs succeed");
    CHECK(memcmp(Ra, Rb, sizeof(Ra)) == 0
          && memcmp(ta, tb, sizeof(ta)) == 0
          && deva == devb
          && memcmp(da, db, 50 * sizeof(da[0])) == 0,
          "determinism: bit-identical outputs and diagnostics");
}

int main(void)
{
    test_exact();
    test_noise();
    test_outliers();
    test_clustered();
    test_degenerate();
    test_determinism();
    if (failures) {
        printf("%d hub-solve test(s) FAILED\n", failures);
        return 1;
    }
    printf("all hub-solve tests passed\n");
    return 0;
}
