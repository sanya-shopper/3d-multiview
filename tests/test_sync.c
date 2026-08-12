/* Software synchronization correction (mv/sync.h) against the numbers
 * of doc/multiview.tex "Software synchronization": the reference rig
 * (f = 800 px, B = 0.5 m, Z ~ 4.7 m) watching a mover at 1 m/s with
 * cameras half a frame (16.7 ms) out of step.  Naive nearest-in-time
 * pairing shows the first-order phantom-depth error f v dt / Z through
 * the depth law; interpolating the tracked detections onto the common
 * clock leaves only the second-order curvature term a dt^2 / 8.
 * OWNERSHIP (parallel build): SYNC-CORRECTION work item ONLY. */

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mv/core.h"
#include "mv/cam.h"
#include "mv/triangulate.h"
#include "mv/sync.h"

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

/* ------------------------------------------------------------------ */
/* fixture: the reference rig and a laterally moving point            */

#define FPS 30.0
#define PERIOD (1.0 / FPS)
#define DT_HALF (1.0 / 60.0) /* half a frame period, the doc's 16.7 ms */
#define NSAMP 25

#define FOCAL 800.0
#define BASELINE 0.5
#define DEPTH 4.7

/* mover: x(t) = X0 + v t + a t^2 / 2, constant y and z */
#define X0 (-0.4)
#define Y0 0.1

static void make_rig(mv_camera *ca, mv_camera *cb)
{
    double C[3] = { BASELINE, 0.0, 0.0 };
    mv_cam_set_K(ca, FOCAL, FOCAL, 320.0, 240.0);
    mv_cam_set_identity_pose(ca);
    memset(ca->k, 0, sizeof(ca->k));
    *cb = *ca;
    mv_cam_set_pose_yaw(cb, 0.0, C); /* pure translation along x */
}

static void point_at(double X[3], double t, double v, double a)
{
    X[0] = X0 + v * t + 0.5 * a * t * t;
    X[1] = Y0;
    X[2] = DEPTH;
}

/* sample the mover through cam at t0 + i * PERIOD; returns MV_OK */
static int make_track(double *t, double *uv, int n, const mv_camera *cam,
                      double t0, double v, double a)
{
    double X[3];
    int i;
    for (i = 0; i < n; i++) {
        t[i] = t0 + i * PERIOD;
        point_at(X, t[i], v, a);
        if (mv_cam_project(uv + 2 * i, cam, X) != MV_OK)
            return MV_ERR;
    }
    return MV_OK;
}

/* ------------------------------------------------------------------ */
/* baseline: naive nearest-in-time pairing + triangulation            */

static double naive_pair_err(double *mean_err, const mv_camera *ca,
                             const mv_camera *cb, const double *ta,
                             const double *uva, int na, const double *tb,
                             const double *uvb, int nb, double v,
                             double acc)
{
    const mv_camera *cams[2];
    double obs[4], Xtri[3], Xtrue[3], maxe = 0.0, sum = 0.0;
    int i, j;
    cams[0] = ca;
    cams[1] = cb;
    for (i = 0; i < na; i++) {
        double e, dx, dy, dz;
        int jbest = 0;
        for (j = 1; j < nb; j++)
            if (fabs(tb[j] - ta[i]) < fabs(tb[jbest] - ta[i]))
                jbest = j;
        obs[0] = uva[2 * i + 0];
        obs[1] = uva[2 * i + 1];
        obs[2] = uvb[2 * jbest + 0];
        obs[3] = uvb[2 * jbest + 1];
        if (mv_triangulate(Xtri, cams, obs, 2) != MV_OK)
            return -1.0;
        point_at(Xtrue, ta[i], v, acc);
        dx = Xtri[0] - Xtrue[0];
        dy = Xtri[1] - Xtrue[1];
        dz = Xtri[2] - Xtrue[2];
        e = sqrt(dx * dx + dy * dy + dz * dz);
        sum += e;
        if (e > maxe)
            maxe = e;
    }
    if (mean_err)
        *mean_err = sum / na;
    return maxe;
}

/* corrected pipeline error: mv_sync_pair_triangulate vs ground truth
 * at each pair's own timestamp; returns max 3-D error, np via pointer */
static double corrected_err(int *np_out, const mv_camera *ca,
                            const mv_camera *cb, const mv_sync_track *a,
                            const mv_sync_track *b, double v, double acc)
{
    double X[3 * NSAMP], tp[NSAMP], Xtrue[3], maxe = 0.0;
    int np, i;
    np = mv_sync_pair_triangulate(X, tp, ca, cb, a, b, 0.5);
    *np_out = np;
    if (np <= 0)
        return -1.0;
    for (i = 0; i < np; i++) {
        double dx, dy, dz, e;
        point_at(Xtrue, tp[i], v, acc);
        dx = X[3 * i + 0] - Xtrue[0];
        dy = X[3 * i + 1] - Xtrue[1];
        dz = X[3 * i + 2] - Xtrue[2];
        e = sqrt(dx * dx + dy * dy + dz * dz);
        if (e > maxe)
            maxe = e;
    }
    return maxe;
}

int main(void)
{
    mv_camera ca, cb;
    double ta[NSAMP], uva[2 * NSAMP], tb[NSAMP], uvb[2 * NSAMP];
    mv_sync_track tra, trb;
    double naive_lin, naive_lin_mean, corr_lin;
    double naive_q, naive_q_mean, corr_q;
    int np;

    make_rig(&ca, &cb);
    tra.t = ta;
    tra.uv = uva;
    tra.n = NSAMP;
    trb.t = tb;
    trb.uv = uvb;
    trb.n = NSAMP;

    /* ---- (a)+(b) linear mover at 1 m/s, B lags by half a frame ---- */
    CHECK(make_track(ta, uva, NSAMP, &ca, 0.0, 1.0, 0.0) == MV_OK
              && make_track(tb, uvb, NSAMP, &cb, DT_HALF, 1.0, 0.0) == MV_OK,
          "sync: fixture tracks project cleanly");

    naive_lin = naive_pair_err(&naive_lin_mean, &ca, &cb, ta, uva, NSAMP,
                               tb, uvb, NSAMP, 1.0, 0.0);
    printf("      (naive nearest-in-time, dt=16.7ms, v=1m/s: max 3-D"
           " error %.1f mm, mean %.1f mm; doc quotes ~150 mm)\n",
           naive_lin * 1e3, naive_lin_mean * 1e3);
    CHECK(naive_lin > 0.05 && naive_lin < 0.5,
          "sync: naive pairing shows phantom depth error of order"
          " f*v*dt/Z (tens to hundreds of mm)");

    corr_lin = corrected_err(&np, &ca, &cb, &tra, &trb, 1.0, 0.0);
    printf("      (corrected: %d pairs, max 3-D error %.2e mm)\n", np,
           corr_lin * 1e3);
    CHECK(np >= NSAMP - 1, "sync: corrected pairing covers the track");
    CHECK(corr_lin >= 0.0 && corr_lin < 1e-6,
          "sync: linear motion corrected to sub-micrometre (interpolation"
          " exact)");
    CHECK(corr_lin > 0.0 && naive_lin / corr_lin > 100.0,
          "sync: correction wins >= two orders of magnitude");

    /* interpolation itself reproduces the exact detection mid-frame */
    {
        double uv[2], uvt[2], Xq[3], tq = 10 * PERIOD; /* between B samples */
        point_at(Xq, tq, 1.0, 0.0);
        CHECK(mv_cam_project(uvt, &cb, Xq) == MV_OK
                  && mv_sync_interp(uv, &trb, tq, 0.5) == MV_OK
                  && fabs(uv[0] - uvt[0]) < 1e-9
                  && fabs(uv[1] - uvt[1]) < 1e-9,
              "sync: mv_sync_interp matches exact projection on linear"
              " track");
    }

    /* ---- (c) quadratic mover, a = 2 m/s^2: residual is 2nd order ---- */
    CHECK(make_track(ta, uva, NSAMP, &ca, 0.0, 1.0, 2.0) == MV_OK
              && make_track(tb, uvb, NSAMP, &cb, DT_HALF, 1.0, 2.0) == MV_OK,
          "sync: quadratic fixture tracks project cleanly");

    naive_q = naive_pair_err(&naive_q_mean, &ca, &cb, ta, uva, NSAMP, tb,
                             uvb, NSAMP, 1.0, 2.0);
    corr_q = corrected_err(&np, &ca, &cb, &tra, &trb, 1.0, 2.0);
    printf("      (quadratic a=2m/s^2: naive max %.1f mm, corrected max"
           " %.3f mm over %d pairs)\n",
           naive_q * 1e3, corr_q * 1e3, np);
    /* curvature term: every query lands mid-frame, so the apparent
     * lateral offset is a*T^2/8 = 278 um (doc's a*dt^2/8 = 70 um for
     * its dt = 16.7 ms convention -- same second order); the depth law
     * amplifies laterally-coded depth by Z/B ~ 9.4, predicting ~2.6 mm
     * for interior pairs.  The max lands on the first pair, which
     * extrapolates B half a frame back: a/2*dt*(dt+T) ~ 3x the chord
     * term, ~7.8 mm -- still second-order, still sub-cm. */
    CHECK(corr_q > 5e-5 && corr_q < 1e-2,
          "sync: quadratic residual is second-order small (sub-cm, order"
          " a*dt^2/8 through the depth law)");
    CHECK(corr_q > 0.0 && naive_q / corr_q > 20.0,
          "sync: correction still wins big on curved motion");
    {
        /* detection-space curvature residual, converted to metres at the
         * mover: expect ~a*T^2/8 = 278 um, doc order of magnitude */
        double uv[2], uvt[2], Xq[3], maxdx = 0.0;
        int i, allok = 1;
        for (i = 1; i < NSAMP; i++) { /* ta[0] extrapolates; skip */
            double dx;
            point_at(Xq, ta[i], 1.0, 2.0);
            if (mv_cam_project(uvt, &cb, Xq) != MV_OK
                || mv_sync_interp(uv, &trb, ta[i], 0.5) != MV_OK) {
                allok = 0;
                break;
            }
            dx = fabs(uv[0] - uvt[0]) * DEPTH / FOCAL;
            if (dx > maxdx)
                maxdx = dx;
        }
        CHECK(allok, "sync: quadratic interp queries succeed");
        printf("      (detection-space curvature residual: max %.0f um;"
               " a*T^2/8 = %.0f um)\n",
               maxdx * 1e6, 2.0 * PERIOD * PERIOD / 8.0 * 1e6);
        CHECK(maxdx > 3e-5 && maxdx < 1e-3,
              "sync: curvature residual matches a*dt^2/8 order (tens to"
              " hundreds of um)");
    }

    /* ---- (d) dt = 0 degenerates to the static (simultaneous) case ---- */
    CHECK(make_track(ta, uva, NSAMP, &ca, 0.0, 1.0, 0.0) == MV_OK
              && make_track(tb, uvb, NSAMP, &cb, 0.0, 1.0, 0.0) == MV_OK,
          "sync: dt=0 fixture tracks project cleanly");
    corr_lin = corrected_err(&np, &ca, &cb, &tra, &trb, 1.0, 0.0);
    naive_lin = naive_pair_err(NULL, &ca, &cb, ta, uva, NSAMP, tb, uvb,
                               NSAMP, 1.0, 0.0);
    printf("      (dt=0: %d pairs, corrected max %.2e mm, naive max"
           " %.2e mm)\n", np, corr_lin * 1e3, naive_lin * 1e3);
    CHECK(np == NSAMP, "sync: dt=0 pairs every sample");
    CHECK(corr_lin >= 0.0 && corr_lin < 1e-9
              && fabs(corr_lin - naive_lin) < 1e-9,
          "sync: dt=0 reproduces the static result exactly");

    /* ---- (e) MV_ERR paths ---- */
    {
        static const double tbad[3] = { 0.0, 1.0, 1.0 }; /* not increasing */
        static const double uvbad[6] = { 0, 0, 1, 1, 2, 2 };
        static const double one_t[1] = { 5.0 };
        static const double one_uv[2] = { 7.0, 9.0 };
        double uv[2];
        mv_sync_track bad, single, empty;
        double po[2 * NSAMP], pt[2 * NSAMP], tt[NSAMP];
        bad.t = tbad;
        bad.uv = uvbad;
        bad.n = 3;
        single.t = one_t;
        single.uv = one_uv;
        single.n = 1;
        empty.t = tbad;
        empty.uv = uvbad;
        empty.n = 0;

        CHECK(mv_sync_interp(NULL, &trb, 0.1, 0.5) == MV_ERR
                  && mv_sync_interp(uv, NULL, 0.1, 0.5) == MV_ERR
                  && mv_sync_interp(uv, &empty, 0.1, 0.5) == MV_ERR
                  && mv_sync_interp(uv, &bad, 0.5, 0.5) == MV_ERR,
              "sync: interp rejects null/empty/non-increasing tracks");
        CHECK(mv_sync_interp(uv, &trb, -1.0, 0.5) == MV_ERR
                  && mv_sync_interp(uv, &trb,
                                    tb[NSAMP - 1] + 0.6 * PERIOD,
                                    0.5) == MV_ERR
                  && mv_sync_interp(uv, &trb,
                                    tb[NSAMP - 1] + 0.4 * PERIOD,
                                    0.5) == MV_OK,
              "sync: extrapolation limit enforced at both ends");
        CHECK(mv_sync_interp(uv, &trb, 0.1, -0.1) == MV_ERR
                  && mv_sync_interp(uv, &trb, nan(""), 0.5) == MV_ERR,
              "sync: negative extrap and NaN query rejected");
        CHECK(mv_sync_interp(uv, &single, 5.0, 0.5) == MV_OK
                  && uv[0] == 7.0 && uv[1] == 9.0
                  && mv_sync_interp(uv, &single, 5.0001, 0.5) == MV_ERR,
              "sync: single-sample track answers only its own instant");
        CHECK(mv_sync_pair(NULL, NULL, po, &tra, &trb, 0.5) == MV_ERR
                  && mv_sync_pair(NULL, po, NULL, &tra, &trb, 0.5) == MV_ERR
                  && mv_sync_pair(tt, po, pt, &bad, &trb, 0.5) == MV_ERR,
              "sync: pair rejects bad arguments");
        CHECK(mv_sync_pair_triangulate(NULL, NULL, &ca, &cb, &tra, &trb,
                                       0.5) == MV_ERR,
              "sync: pair_triangulate rejects bad arguments");
        /* disjoint coverage: zero pairs is an answer, not an error */
        CHECK(mv_sync_pair(tt, po, pt, &tra, &single, 0.5) == 0,
              "sync: disjoint tracks pair to zero, not MV_ERR");
    }

    if (failures) {
        printf("\n%d FAILURE(S)\n", failures);
        return 1;
    }
    printf("\nall sync tests passed\n");
    return 0;
}
