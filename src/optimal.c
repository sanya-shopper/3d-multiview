#include <math.h>
#include <string.h>

#include "mv/core.h"
#include "mv/mat.h"
#include "mv/epipolar.h"
#include "mv/triangulate.h"
#include "mv/optimal.h"

#define OPT_PI 3.14159265358979323846

/* --- Hartley-Sturm optimal two-view correction --------------------------- */

/* out (na+nb-1 coeffs, out[i] multiplies t^i) = a * b; no aliasing. */
static void poly_mul(double *out, const double *a, int na,
                     const double *b, int nb)
{
    int i, j;
    for (i = 0; i < na + nb - 1; i++)
        out[i] = 0.0;
    for (i = 0; i < na; i++)
        for (j = 0; j < nb; j++)
            out[i + j] += a[i] * b[j];
}

static double poly_eval(const double *p, int n, double t)
{
    double v = 0.0;
    int i;
    for (i = n - 1; i >= 0; i--)
        v = v * t + p[i];
    return v;
}

/* Cost s(t): squared distances from both (origin-translated) measured
 * points to the epipolar-line pair indexed by t (HS eq. for d1^2+d2^2). */
static double hs_cost(double t, double f1, double f2,
                      double pa, double pb, double pc, double pd)
{
    double p = pa * t + pb, q = pc * t + pd;
    double den2 = p * p + f2 * f2 * q * q;
    double cost = t * t / (1.0 + f1 * f1 * t * t);
    if (den2 > 0.0)
        cost += q * q / den2;
    else if (q != 0.0)
        return HUGE_VAL;
    return cost;
}

/* Bisect g (degree ng-1 polynomial) on [lo,hi] with sign(g(lo)) = sign(glo)
 * opposite to sign(g(hi)). */
static double bisect_root(const double *g, int ng, double lo, double hi,
                          double glo)
{
    int i;
    for (i = 0; i < 120; i++) {
        double mid = 0.5 * (lo + hi);
        double gm = poly_eval(g, ng, mid);
        if (gm == 0.0)
            return mid;
        if ((gm > 0.0) == (glo > 0.0)) {
            lo = mid;
            glo = gm;
        } else {
            hi = mid;
        }
    }
    return 0.5 * (lo + hi);
}

int mv_epipolar_correct(double out1[2], double out2[2], const double F[9],
                        const double uv1[2], const double uv2[2])
{
    double T1i[9], T2it[9], tmp[9], Ft[9], Ftt[9];
    double e1[3], e2[3], s1, s2, sc;
    double R1[9], R2[9], R1t[9], R2t[9], R2F[9], Fr[9];
    double f1, f2, pa, pb, pc, pd, k;
    double P[2], Q[2], PP[3], QQ[3], S[3], S2[5];
    double Rp[3], RR[5], PQ[3], RRPQ[7], g[7];
    double best_t, best_cost, prev_t = 0.0, prev_g = 0.0;
    double t, l1[3], l2[3], x1h[3], x2h[3], y1[3], y2[3];
    int i, have_prev;
    const int nsamp = 2048;

    if (!out1 || !out2 || !F || !uv1 || !uv2)
        return MV_ERR;

    /* translate both measured points to their image origins:
     * Ft = T2^-T F T1^-1 */
    memset(T1i, 0, sizeof(T1i));
    T1i[0] = T1i[4] = T1i[8] = 1.0;
    T1i[2] = uv1[0];
    T1i[5] = uv1[1];
    memset(T2it, 0, sizeof(T2it));
    T2it[0] = T2it[4] = T2it[8] = 1.0;
    T2it[6] = uv2[0];
    T2it[7] = uv2[1];
    mv_mat_mul(tmp, F, T1i, 3, 3, 3);
    mv_mat_mul(Ft, T2it, tmp, 3, 3, 3);

    /* epipoles of the translated F: Ft e1 = 0, Ft^T e2 = 0 */
    if (mv_nullvec(e1, Ft, 3, 3) != MV_OK)
        return MV_ERR;
    mv_mat_transpose(Ftt, Ft, 3, 3);
    if (mv_nullvec(e2, Ftt, 3, 3) != MV_OK)
        return MV_ERR;
    s1 = e1[0] * e1[0] + e1[1] * e1[1];
    s2 = e2[0] * e2[0] + e2[1] * e2[1];
    /* epipole at (or numerically indistinguishable from) a measured
     * point: the 1-D reduction is undefined */
    if (s1 < 1e-16 || s2 < 1e-16)
        return MV_ERR;
    sc = 1.0 / sqrt(s1);
    for (i = 0; i < 3; i++)
        e1[i] *= sc;
    sc = 1.0 / sqrt(s2);
    for (i = 0; i < 3; i++)
        e2[i] *= sc;

    /* rotate both epipoles onto the x-axes: e -> (1, 0, f) */
    memset(R1, 0, sizeof(R1));
    R1[0] = e1[0];
    R1[1] = e1[1];
    R1[3] = -e1[1];
    R1[4] = e1[0];
    R1[8] = 1.0;
    memset(R2, 0, sizeof(R2));
    R2[0] = e2[0];
    R2[1] = e2[1];
    R2[3] = -e2[1];
    R2[4] = e2[0];
    R2[8] = 1.0;
    mv_mat_transpose(R1t, R1, 3, 3);
    mv_mat_transpose(R2t, R2, 3, 3);
    mv_mat_mul(R2F, R2, Ft, 3, 3, 3);
    mv_mat_mul(Fr, R2F, R1t, 3, 3, 3);

    f1 = e1[2];
    f2 = e2[2];
    pa = Fr[4];
    pb = Fr[5];
    pc = Fr[7];
    pd = Fr[8];

    /* g(t) = t (p^2 + f2^2 q^2)^2 - (ad-bc)(1 + f1^2 t^2)^2 p q,
     * p = a t + b, q = c t + d: numerator of s'(t), degree 6 */
    P[0] = pb;
    P[1] = pa;
    Q[0] = pd;
    Q[1] = pc;
    poly_mul(PP, P, 2, P, 2);
    poly_mul(QQ, Q, 2, Q, 2);
    for (i = 0; i < 3; i++)
        S[i] = PP[i] + f2 * f2 * QQ[i];
    poly_mul(S2, S, 3, S, 3);
    Rp[0] = 1.0;
    Rp[1] = 0.0;
    Rp[2] = f1 * f1;
    poly_mul(RR, Rp, 3, Rp, 3);
    poly_mul(PQ, P, 2, Q, 2);
    poly_mul(RRPQ, RR, 5, PQ, 3);
    k = pa * pd - pb * pc;
    for (i = 0; i < 7; i++)
        g[i] = -k * RRPQ[i];
    for (i = 0; i < 5; i++)
        g[i + 1] += S2[i];

    /* every strict local minimum of s is a sign change of g: sample the
     * whole real line through t = tan(theta), bisect each bracket, and
     * keep the cheapest candidate (plus t = 0 and t -> infinity) */
    best_t = 0.0;
    best_cost = hs_cost(0.0, f1, f2, pa, pb, pc, pd);
    if (f1 != 0.0) {
        double den = pa * pa + f2 * f2 * pc * pc;
        double cinf = 1.0 / (f1 * f1);
        if (den > 0.0)
            cinf += pc * pc / den;
        if (cinf < best_cost) {
            best_cost = cinf;
            best_t = 1e12;
        }
    }
    have_prev = 0;
    for (i = 0; i < nsamp; i++) {
        double th = -0.5 * OPT_PI + OPT_PI * (i + 0.5) / nsamp;
        double ts = tan(th);
        double gv = poly_eval(g, 7, ts);
        double cv = hs_cost(ts, f1, f2, pa, pb, pc, pd);
        if (cv < best_cost) {
            best_cost = cv;
            best_t = ts;
        }
        if (have_prev && gv != 0.0 && prev_g != 0.0 &&
            (gv > 0.0) != (prev_g > 0.0)) {
            double r = bisect_root(g, 7, prev_t, ts, prev_g);
            double cr = hs_cost(r, f1, f2, pa, pb, pc, pd);
            if (cr < best_cost) {
                best_cost = cr;
                best_t = r;
            }
        }
        prev_t = ts;
        prev_g = gv;
        have_prev = 1;
    }
    if (!(best_cost < HUGE_VAL))
        return MV_ERR;

    /* closest points to the two origins on the optimal line pair */
    t = best_t;
    l1[0] = t * f1;
    l1[1] = 1.0;
    l1[2] = -t;
    l2[0] = -f2 * (pc * t + pd);
    l2[1] = pa * t + pb;
    l2[2] = pc * t + pd;
    x1h[0] = -l1[0] * l1[2];
    x1h[1] = -l1[1] * l1[2];
    x1h[2] = l1[0] * l1[0] + l1[1] * l1[1];
    x2h[0] = -l2[0] * l2[2];
    x2h[1] = -l2[1] * l2[2];
    x2h[2] = l2[0] * l2[0] + l2[1] * l2[1];

    /* undo the rotations and translations */
    mv_mat_mul(y1, R1t, x1h, 3, 3, 1);
    mv_mat_mul(y2, R2t, x2h, 3, 3, 1);
    if (fabs(y1[2]) <= 1e-14 * (fabs(y1[0]) + fabs(y1[1])) ||
        fabs(y2[2]) <= 1e-14 * (fabs(y2[0]) + fabs(y2[1])))
        return MV_ERR;
    out1[0] = y1[0] / y1[2] + uv1[0];
    out1[1] = y1[1] / y1[2] + uv1[1];
    out2[0] = y2[0] / y2[2] + uv2[0];
    out2[1] = y2[1] / y2[2] + uv2[1];
    if (!isfinite(out1[0]) || !isfinite(out1[1]) ||
        !isfinite(out2[0]) || !isfinite(out2[1]))
        return MV_ERR;
    return MV_OK;
}

int mv_triangulate_optimal(double X[3], const mv_camera *c1,
                           const mv_camera *c2, const double uv1[2],
                           const double uv2[2])
{
    double F[9], uv[4];
    const mv_camera *cams[2];

    if (!X || !c1 || !c2 || !uv1 || !uv2)
        return MV_ERR;
    cams[0] = c1;
    cams[1] = c2;
    mv_fundamental_from_cams(F, c1, c2);
    if (mv_epipolar_correct(uv, uv + 2, F, uv1, uv2) != MV_OK) {
        /* degenerate epipole-at-point geometry: plain triangulation */
        uv[0] = uv1[0];
        uv[1] = uv1[1];
        uv[2] = uv2[0];
        uv[3] = uv2[1];
    }
    return mv_triangulate(X, cams, uv, 2);
}

/* --- Mahalanobis-weighted known-radius circle fit ------------------------ */

/* Unweighted Gauss-Newton fit of a known-radius circle (the estimator
 * used by demo/track_robot.c); serves as the initializer. */
static void circle_fit_r(double c[2], const double *pts, int n, double R)
{
    int i, it;
    c[0] = c[1] = 0.0;
    for (i = 0; i < n; i++) {
        c[0] += pts[2 * i];
        c[1] += pts[2 * i + 1];
    }
    c[0] /= n;
    c[1] /= n;
    for (it = 0; it < 50; it++) {
        double a = 0.0, b = 0.0, d2 = 0.0, gx = 0.0, gz = 0.0, det, sx, sz;
        for (i = 0; i < n; i++) {
            double dx = pts[2 * i] - c[0], dz = pts[2 * i + 1] - c[1];
            double d = sqrt(dx * dx + dz * dz), rr, jx, jz;
            if (d < 1e-9)
                continue;
            rr = d - R;
            jx = -dx / d;
            jz = -dz / d;
            a += jx * jx;
            b += jx * jz;
            d2 += jz * jz;
            gx += jx * rr;
            gz += jz * rr;
        }
        det = a * d2 - b * b;
        if (det <= 1e-12 * (a + d2) * (a + d2))
            break;
        sx = -(d2 * gx - b * gz) / det;
        sz = -(a * gz - b * gx) / det;
        c[0] += sx;
        c[1] += sz;
        if (fabs(sx) + fabs(sz) < 1e-14)
            break;
    }
}

/* Squared Mahalanobis distance from the (center-relative) point (px,pz)
 * to the circle point at angle th, for inverse covariance
 * [S00 S01; S01 S11]. */
static double maha_at(double px, double pz, double R, double th,
                      double S00, double S01, double S11)
{
    double dx = px - R * cos(th), dz = pz - R * sin(th);
    return dx * dx * S00 + 2.0 * dx * dz * S01 + dz * dz * S11;
}

int mv_circle_fit_weighted(double center[2], const double *pts,
                           const double *dirs, const double *sig,
                           int n, double radius)
{
    int i, it, j, q;
    const int ncoarse = 128;
    const double gr = 0.6180339887498949; /* golden ratio - 1 */

    if (!center || !pts || !dirs || !sig || n < 2 || radius <= 0.0)
        return MV_ERR;
    circle_fit_r(center, pts, n, radius);
    for (it = 0; it < 100; it++) {
        double A = 0.0, B = 0.0, D = 0.0, g0 = 0.0, g1 = 0.0;
        double det, sx, sz;
        for (i = 0; i < n; i++) {
            double ax = dirs[2 * i], az = dirs[2 * i + 1];
            double va = sig[2 * i] * sig[2 * i];
            double vc = sig[2 * i + 1] * sig[2 * i + 1];
            double ia, ic, S00, S01, S11;
            double px, pz, bt, bm, lo, hi, x1, x2, f1, f2, dx, dz;
            if (va < 1e-18)
                va = 1e-18;
            if (vc < 1e-18)
                vc = 1e-18;
            ia = 1.0 / va;
            ic = 1.0 / vc;
            /* Sigma^-1 = a a^T / va + b b^T / vc, b = perp(a) */
            S00 = ax * ax * ia + az * az * ic;
            S01 = ax * az * (ia - ic);
            S11 = az * az * ia + ax * ax * ic;
            px = pts[2 * i] - center[0];
            pz = pts[2 * i + 1] - center[1];
            /* Mahalanobis-closest circle point: coarse scan over the
             * whole circle, then golden-section polish on the bracket */
            bt = 0.0;
            bm = HUGE_VAL;
            for (j = 0; j < ncoarse; j++) {
                double th = 2.0 * OPT_PI * j / ncoarse;
                double m = maha_at(px, pz, radius, th, S00, S01, S11);
                if (m < bm) {
                    bm = m;
                    bt = th;
                }
            }
            lo = bt - 2.0 * OPT_PI / ncoarse;
            hi = bt + 2.0 * OPT_PI / ncoarse;
            x1 = hi - gr * (hi - lo);
            x2 = lo + gr * (hi - lo);
            f1 = maha_at(px, pz, radius, x1, S00, S01, S11);
            f2 = maha_at(px, pz, radius, x2, S00, S01, S11);
            for (q = 0; q < 60; q++) {
                if (f1 < f2) {
                    hi = x2;
                    x2 = x1;
                    f2 = f1;
                    x1 = hi - gr * (hi - lo);
                    f1 = maha_at(px, pz, radius, x1, S00, S01, S11);
                } else {
                    lo = x1;
                    x1 = x2;
                    f1 = f2;
                    x2 = lo + gr * (hi - lo);
                    f2 = maha_at(px, pz, radius, x2, S00, S01, S11);
                }
            }
            bt = 0.5 * (lo + hi);
            dx = px - radius * cos(bt);
            dz = pz - radius * sin(bt);
            /* envelope theorem: d(cost_i)/d(center) = -2 Sigma^-1 dv,
             * GN Hessian contribution = Sigma^-1 */
            A += S00;
            B += S01;
            D += S11;
            g0 += S00 * dx + S01 * dz;
            g1 += S01 * dx + S11 * dz;
        }
        det = A * D - B * B;
        if (det <= 1e-12 * (A + D) * (A + D))
            break;
        sx = (D * g0 - B * g1) / det;
        sz = (A * g1 - B * g0) / det;
        center[0] += sx;
        center[1] += sz;
        if (fabs(sx) + fabs(sz) < 1e-14)
            break;
    }
    if (!isfinite(center[0]) || !isfinite(center[1]))
        return MV_ERR;
    return MV_OK;
}
