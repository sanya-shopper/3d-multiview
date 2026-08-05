/* Rig pair solver -- robust joint estimate of the rigid transform
 * x_a = R x_b + t from matched display-pose pairs.
 *
 * Pipeline: (1) seed from the better of the v0 chordal-mean+median
 * estimate and the medoid candidate; (2) IRLS with Huber weights over
 * a combined SE(3) residual to reach the majority mode robustly;
 * (3) a final detect-and-refit -- hard-trim residuals beyond
 * max(4 x median, floor) and take the plain mean over the surviving
 * inliers.  Step 3 exists because Huber's soft rejection leaves a
 * systematic bias when outliers are coherent (clustered anchors from
 * a briefly-degenerate calibration all pull the same way); the trim
 * removes them outright, and on clean data it keeps everything, so
 * the result gracefully reduces to the chordal mean + mean
 * translation.  Everything is deterministic: no RNG, and every
 * ordering used is a total order on the input values.
 *
 * OWNERSHIP (parallel build): joint-solve work item ONLY. */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "mv/mv.h"

#include "hub_solve.h"

#define SOLVE_MAXOBS 256

/* Translation lever arm (m): converts translation residuals into
 * radian-equivalents so rotation and translation share one robust
 * scale.  0.5 m matches the display distances of the live sessions
 * (0.3-0.7 m), where a 1 mrad pose error moves the plane ~0.5 mm. */
#define SOLVE_LEVER_M 0.5

/* Huber threshold delta = SOLVE_HUBER_K * median combined residual,
 * floored so exact or near-exact data cannot produce delta = 0 (which
 * would collapse every weight to 0/0).  1e-6 rad-equivalent is 0.2 mdeg
 * or 0.5 um -- far below any real signal, so the floor never shields an
 * actual outlier. */
#define SOLVE_HUBER_K 1.5
#define SOLVE_DELTA_FLOOR 1e-6

#define SOLVE_MAX_ITER 30
#define SOLVE_CONV_TOL 1e-13

/* Final hard trim: reject residuals beyond SOLVE_TRIM_K x the median,
 * floored at SOLVE_TRIM_FLOOR (0.05 rad-equivalent = 2.9 deg / 25 mm)
 * so plausible inlier noise (<= ~1 deg / 10 mm) is never trimmed even
 * when the median residual is tiny, while gross outliers (>= 30 deg /
 * 0.2 m, the live failure mode) always are. */
#define SOLVE_TRIM_K 4.0
#define SOLVE_TRIM_FLOOR 0.05

/* Combined SE(3) distance between two transforms, in radian
 * equivalents: geodesic rotation angle (via the canonical mv_rot_log,
 * see CODING_STANDARDS "one kernel") plus lever-scaled translation
 * distance, combined in quadrature. */
static double pair_dist(const double Ra[9], const double ta[3],
                        const double Rb[9], const double tb[3])
{
    double At[9], D[9], r[3], d[3];
    double ang, dt;
    int k;

    mv_mat_transpose(At, Ra, 3, 3);
    mv_mat_mul(D, At, Rb, 3, 3, 3);
    mv_rot_log(r, D);
    ang = mv_norm(r, 3);
    for (k = 0; k < 3; k++)
        d[k] = tb[k] - ta[k];
    dt = mv_norm(d, 3) / SOLVE_LEVER_M;
    return sqrt(ang * ang + dt * dt);
}

/* Sum of the h smallest combined residuals of the n candidates against
 * (R, t).  Trimming makes the score blind to any minority of outliers,
 * so it ranks competing seeds by how well they fit the majority. */
static double trimmed_cost(const double R[9], const double t[3],
                           const double (*Rc)[9], const double (*tc)[3],
                           int n, int h)
{
    double col[SOLVE_MAXOBS] = { 0 };
    double s = 0.0;
    int i;

    for (i = 0; i < n; i++)
        col[i] = pair_dist(R, t, Rc[i], tc[i]);
    qsort(col, (size_t)n, sizeof(double), mv_cmp_double);
    for (i = 0; i < h; i++)
        s += col[i];
    return s;
}

int hub_solve_pair_diag(double R[9], double t[3], double *dev_mm,
                        hub_solve_diag *diag, const hub_obs *obs, int n)
{
    double Rc[SOLVE_MAXOBS][9], tc[SOLVE_MAXOBS][3];
    double res[SOLVE_MAXOBS] = { 0 }, w[SOLVE_MAXOBS] = { 0 };
    double col[SOLVE_MAXOBS] = { 0 };
    double Rv0[9], tv0[3], Rmed[9], tmed[3];
    double dev = 0.0;
    int i, k, h, it;

    if (n < 3)
        return -1;
    if (n > SOLVE_MAXOBS)
        n = SOLVE_MAXOBS;

    /* Per-pair candidate transforms: R_i = Ra_i Rb_i^T maps camera b to
     * camera a; t_i = ta_i - R_i tb_i. */
    for (i = 0; i < n; i++) {
        double Rbt[9], tt[3];
        mv_mat_transpose(Rbt, obs[i].Rb, 3, 3);
        mv_mat_mul(Rc[i], obs[i].Ra, Rbt, 3, 3, 3);
        mv_mat_mul(tt, Rc[i], obs[i].tb, 3, 3, 1);
        for (k = 0; k < 3; k++)
            tc[i][k] = obs[i].ta[k] - tt[k];
    }

    /* Trim count: keep the smallest 60% (>= 2) residuals when scoring
     * seeds, matching the worst outlier fraction we must survive
     * (~30%) with margin. */
    h = (3 * n) / 5;
    if (h < 2)
        h = 2;

    /* Seed A: v0 chordal mean + per-axis median (statistically
     * efficient when the data are clean). */
    {
        double Rsum[9] = { 0 };
        for (i = 0; i < n; i++)
            for (k = 0; k < 9; k++)
                Rsum[k] += Rc[i][k];
        if (mv_rot_project(Rv0, Rsum) != MV_OK)
            return -1;
        for (k = 0; k < 3; k++) {
            for (i = 0; i < n; i++)
                col[i] = tc[i][k];
            qsort(col, (size_t)n, sizeof(double), mv_cmp_double);
            tv0[k] = col[n / 2];
        }
    }

    /* Seed B: the medoid candidate -- the observation minimizing the
     * trimmed sum of distances to all others.  For any outlier
     * fraction < 50% the medoid is an inlier. */
    {
        double best = HUGE_VAL;
        int m = 0, j;
        for (i = 0; i < n; i++) {
            double s = 0.0;
            for (j = 0; j < n; j++)
                col[j] = pair_dist(Rc[i], tc[i], Rc[j], tc[j]);
            qsort(col, (size_t)n, sizeof(double), mv_cmp_double);
            for (j = 0; j < h; j++)
                s += col[j];
            if (s < best) {
                best = s;
                m = i;
            }
        }
        memcpy(Rmed, Rc[m], sizeof(Rmed));
        memcpy(tmed, tc[m], sizeof(tmed));
    }

    /* Clustered outliers (the live failure mode: a subset of anchors
     * from a briefly-degenerate calibration shares one systematically
     * wrong transform) make the population bimodal, and Huber IRLS
     * converges to whichever mode dominates its seed.  The v0 chordal
     * mean lands between the modes, so seeding from it alone is
     * unsafe.  Both seeds are therefore scored by the trimmed cost,
     * which the majority mode always wins, and IRLS then stays in that
     * mode: once the majority carries weight ~1 the minority cannot
     * dominate the weighted sums.  With small clean n (3-5) the seeds
     * nearly coincide and the refinement reduces to v0-like weighted
     * averaging rather than rejecting anything. */
    if (trimmed_cost(Rmed, tmed,
                     (const double (*)[9])Rc, (const double (*)[3])tc,
                     n, h) <
        trimmed_cost(Rv0, tv0,
                     (const double (*)[9])Rc, (const double (*)[3])tc,
                     n, h)) {
        memcpy(R, Rmed, sizeof(Rmed));
        memcpy(t, tmed, sizeof(tmed));
    } else {
        memcpy(R, Rv0, sizeof(Rv0));
        memcpy(t, tv0, sizeof(tv0));
    }

    /* IRLS refinement: Huber weights on the combined residual, then a
     * weighted chordal mean (SVD projection) and weighted mean
     * translation.  Huber weights are strictly positive, so no
     * observation is ever fully discarded and the weighted sums stay
     * well-conditioned. */
    for (it = 0; it < SOLVE_MAX_ITER; it++) {
        double Rsum[9] = { 0 }, tnew[3] = { 0 }, Rnew[9];
        double wsum = 0.0, step, delta;

        for (i = 0; i < n; i++)
            res[i] = pair_dist(R, t, Rc[i], tc[i]);
        memcpy(col, res, (size_t)n * sizeof(double));
        qsort(col, (size_t)n, sizeof(double), mv_cmp_double);
        delta = SOLVE_HUBER_K * col[n / 2];
        if (delta < SOLVE_DELTA_FLOOR)
            delta = SOLVE_DELTA_FLOOR;
        for (i = 0; i < n; i++) {
            w[i] = (res[i] <= delta) ? 1.0 : delta / res[i];
            wsum += w[i];
            for (k = 0; k < 9; k++)
                Rsum[k] += w[i] * Rc[i][k];
            for (k = 0; k < 3; k++)
                tnew[k] += w[i] * tc[i][k];
        }
        for (k = 0; k < 3; k++)
            tnew[k] /= wsum;
        if (mv_rot_project(Rnew, Rsum) != MV_OK)
            return -1;
        step = pair_dist(R, t, Rnew, tnew);
        memcpy(R, Rnew, sizeof(Rnew));
        memcpy(t, tnew, sizeof(tnew));
        if (step < SOLVE_CONV_TOL)
            break;
    }

    /* Detect-and-refit: with the majority mode established, hard-trim
     * anything beyond the threshold and take the plain (unweighted)
     * mean over the surviving inliers -- removing the residual bias
     * Huber's strictly-positive weights leave when the outliers are
     * coherent.  Two passes so the inlier set is re-derived from a
     * post-trim estimate.  If fewer than 3 observations survive, the
     * trim is skipped and the IRLS estimate stands (small-n graceful
     * degradation: never reject down to a degenerate set). */
    for (it = 0; it < 2; it++) {
        double Rsum[9] = { 0 }, tnew[3] = { 0 }, Rnew[9];
        double thr;
        int nin = 0;

        for (i = 0; i < n; i++)
            res[i] = pair_dist(R, t, Rc[i], tc[i]);
        memcpy(col, res, (size_t)n * sizeof(double));
        qsort(col, (size_t)n, sizeof(double), mv_cmp_double);
        thr = SOLVE_TRIM_K * col[n / 2];
        if (thr < SOLVE_TRIM_FLOOR)
            thr = SOLVE_TRIM_FLOOR;
        for (i = 0; i < n; i++) {
            w[i] = (res[i] <= thr) ? 1.0 : 0.0;
            if (w[i] > 0.0)
                nin++;
        }
        if (nin < 3) {
            for (i = 0; i < n; i++)
                w[i] = 1.0;
            break;
        }
        for (i = 0; i < n; i++) {
            if (w[i] <= 0.0)
                continue;
            for (k = 0; k < 9; k++)
                Rsum[k] += Rc[i][k];
            for (k = 0; k < 3; k++)
                tnew[k] += tc[i][k];
        }
        for (k = 0; k < 3; k++)
            tnew[k] /= nin;
        if (mv_rot_project(Rnew, Rsum) != MV_OK)
            return -1;
        memcpy(R, Rnew, sizeof(Rnew));
        memcpy(t, tnew, sizeof(tnew));
    }

    /* Final residuals against the final estimate, for diagnostics and
     * the scatter figure.  dev_mm keeps the v0 contract: unweighted
     * mean |t_i - t| in mm over all n observations.  diag weights are
     * the trim weights actually used by the final refit (1 = kept,
     * 0 = rejected as outlier). */
    for (i = 0; i < n; i++) {
        double At[9], D[9], r[3], d[3], dn;
        mv_mat_transpose(At, R, 3, 3);
        mv_mat_mul(D, At, Rc[i], 3, 3, 3);
        mv_rot_log(r, D);
        for (k = 0; k < 3; k++)
            d[k] = tc[i][k] - t[k];
        dn = mv_norm(d, 3);
        dev += dn;
        if (diag) {
            diag[i].rot_err = mv_norm(r, 3);
            diag[i].t_err = dn;
            diag[i].weight = w[i];
        }
    }
    if (dev_mm)
        *dev_mm = 1000.0 * dev / n;
    return 0;
}

int hub_solve_pair(double R[9], double t[3], double *dev_mm,
                   const hub_obs *obs, int n)
{
    return hub_solve_pair_diag(R, t, dev_mm, NULL, obs, n);
}
