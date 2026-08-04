#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "mv/core.h"
#include "mv/mat.h"
#include "mv/calib.h"

/* similarity normalization: centroid 0, mean radius sqrt(2) */
static int norm_transform(double T[9], const double *pts, int n)
{
    double cx = 0.0, cy = 0.0, s = 0.0;
    int i;
    for (i = 0; i < n; i++) {
        cx += pts[2 * i];
        cy += pts[2 * i + 1];
    }
    cx /= n;
    cy /= n;
    for (i = 0; i < n; i++) {
        double dx = pts[2 * i] - cx, dy = pts[2 * i + 1] - cy;
        s += sqrt(dx * dx + dy * dy);
    }
    if (s < 1e-12)
        return MV_ERR;
    s = sqrt(2.0) * n / s;
    memset(T, 0, 9 * sizeof(double));
    T[0] = s; T[2] = -s * cx;
    T[4] = s; T[5] = -s * cy;
    T[8] = 1.0;
    return MV_OK;
}

int mv_homography_dlt(double H[9], const double *xy_plane, const double *uv,
                      int n)
{
    double Tp[9], Ti[9], Tiinv[9], h[9], Hn[9], tmp[9];
    double *A;
    int i, ret;

    if (n < 4)
        return MV_ERR;
    if (norm_transform(Tp, xy_plane, n) != MV_OK
        || norm_transform(Ti, uv, n) != MV_OK)
        return MV_ERR;

    A = (double *)malloc((size_t)(2 * n) * 9 * sizeof(double));
    if (!A)
        return MV_ERR;
    for (i = 0; i < n; i++) {
        double X = Tp[0] * xy_plane[2 * i] + Tp[2];
        double Y = Tp[4] * xy_plane[2 * i + 1] + Tp[5];
        double u = Ti[0] * uv[2 * i] + Ti[2];
        double v = Ti[4] * uv[2 * i + 1] + Ti[5];
        double *r0 = A + (2 * i) * 9;
        double *r1 = A + (2 * i + 1) * 9;
        r0[0] = X;   r0[1] = Y;   r0[2] = 1.0;
        r0[3] = 0.0; r0[4] = 0.0; r0[5] = 0.0;
        r0[6] = -u * X; r0[7] = -u * Y; r0[8] = -u;
        r1[0] = 0.0; r1[1] = 0.0; r1[2] = 0.0;
        r1[3] = X;   r1[4] = Y;   r1[5] = 1.0;
        r1[6] = -v * X; r1[7] = -v * Y; r1[8] = -v;
    }
    ret = mv_nullvec(h, A, 2 * n, 9);
    free(A);
    if (ret != MV_OK)
        return MV_ERR;
    memcpy(Hn, h, sizeof(Hn));

    /* denormalize: H = Ti^-1 Hn Tp */
    memset(Tiinv, 0, sizeof(Tiinv));
    Tiinv[0] = 1.0 / Ti[0]; Tiinv[2] = -Ti[2] / Ti[0];
    Tiinv[4] = 1.0 / Ti[4]; Tiinv[5] = -Ti[5] / Ti[4];
    Tiinv[8] = 1.0;
    mv_mat_mul(tmp, Tiinv, Hn, 3, 3, 3);
    mv_mat_mul(H, tmp, Tp, 3, 3, 3);

    if (fabs(H[8]) > 1e-12)
        for (i = 0; i < 9; i++)
            H[i] /= H[8];
    return MV_OK;
}

/* row of the Zhang constraint system: v_ij from columns i,j of H */
static void vij(double v[6], const double H[9], int i, int j)
{
    double h1i = H[i], h2i = H[3 + i], h3i = H[6 + i];
    double h1j = H[j], h2j = H[3 + j], h3j = H[6 + j];
    v[0] = h1i * h1j;
    v[1] = h1i * h2j + h2i * h1j;
    v[2] = h2i * h2j;
    v[3] = h3i * h1j + h1i * h3j;
    v[4] = h3i * h2j + h2i * h3j;
    v[5] = h3i * h3j;
}

/* intrinsics from B = K^-T K^-1 (up to scale), b = (B11,B12,B22,B13,B23,B33) */
static int intrinsics_from_b(double K[9], double b[6])
{
    double den, v0, lam, alpha, beta, gamma, u0;
    int i;
    if (b[0] < 0.0) /* nullvec sign is arbitrary; B11 must be positive */
        for (i = 0; i < 6; i++)
            b[i] = -b[i];
    den = b[0] * b[2] - b[1] * b[1];
    if (fabs(b[0]) < 1e-300 || fabs(den) < 1e-300)
        return MV_ERR;
    v0 = (b[1] * b[3] - b[0] * b[4]) / den;
    lam = b[5] - (b[3] * b[3] + v0 * (b[1] * b[3] - b[0] * b[4])) / b[0];
    if (lam / b[0] <= 0.0 || lam * b[0] / den <= 0.0)
        return MV_ERR;
    alpha = sqrt(lam / b[0]);
    beta = sqrt(lam * b[0] / den);
    gamma = -b[1] * alpha * alpha * beta / lam;
    u0 = gamma * v0 / beta - b[3] * alpha * alpha / lam;
    memset(K, 0, 9 * sizeof(double));
    K[0] = alpha; K[1] = gamma; K[2] = u0;
    K[4] = beta;  K[5] = v0;
    K[8] = 1.0;
    return MV_OK;
}

/* pose of the target plane from H and K: [r1 r2 t] = K^-1 H up to scale */
static int extrinsics_from_H(mv_camera *cam, const double K[9],
                             const double H[9])
{
    double Kinv[9], h[3], r1[3], r2[3], r3[3], t[3];
    double Q[9], U[9], S[3], V[9], Vt[9], R[9];
    double lam, det;
    int i;

    if (mv_mat_inv3(Kinv, K) != MV_OK)
        return MV_ERR;
    for (i = 0; i < 3; i++)
        h[i] = H[i * 3 + 0];
    mv_mat_mul(r1, Kinv, h, 3, 3, 1);
    for (i = 0; i < 3; i++)
        h[i] = H[i * 3 + 1];
    mv_mat_mul(r2, Kinv, h, 3, 3, 1);
    for (i = 0; i < 3; i++)
        h[i] = H[i * 3 + 2];
    mv_mat_mul(t, Kinv, h, 3, 3, 1);

    lam = mv_norm(r1, 3);
    if (lam < 1e-300)
        return MV_ERR;
    lam = 2.0 / (lam + mv_norm(r2, 3)); /* average the two column scales */
    for (i = 0; i < 3; i++) {
        r1[i] *= lam;
        r2[i] *= lam;
        t[i] *= lam;
    }
    if (t[2] < 0.0) { /* target must be in front of the camera */
        for (i = 0; i < 3; i++) {
            r1[i] = -r1[i];
            r2[i] = -r2[i];
            t[i] = -t[i];
        }
    }
    mv_cross3(r3, r1, r2);

    /* nearest rotation to [r1 r2 r3] via SVD: R = U V^T */
    for (i = 0; i < 3; i++) {
        Q[i * 3 + 0] = r1[i];
        Q[i * 3 + 1] = r2[i];
        Q[i * 3 + 2] = r3[i];
    }
    memcpy(U, Q, sizeof(Q));
    if (mv_svd(U, S, V, 3, 3) != MV_OK)
        return MV_ERR;
    mv_mat_transpose(Vt, V, 3, 3);
    mv_mat_mul(R, U, Vt, 3, 3, 3);
    det = R[0] * (R[4] * R[8] - R[5] * R[7])
        - R[1] * (R[3] * R[8] - R[5] * R[6])
        + R[2] * (R[3] * R[7] - R[4] * R[6]);
    if (det < 0.0)
        return MV_ERR;

    memcpy(cam->K, K, 9 * sizeof(double));
    memcpy(cam->R, R, sizeof(R));
    memcpy(cam->t, t, sizeof(t));
    memset(cam->k, 0, sizeof(cam->k));
    return MV_OK;
}

int mv_calib_plane_pose(mv_camera *cam, const double K[9],
                        const double *obj, const double *img, int n)
{
    double H[9];
    if (n < 4)
        return MV_ERR;
    if (mv_homography_dlt(H, obj, img, n) != MV_OK)
        return MV_ERR;
    return extrinsics_from_H(cam, K, H);
}

/* linear least-squares k1,k2 given intrinsics and per-view poses:
 * u_d - u = (u - u0)(k1 r^2 + k2 r^4), likewise for v (Zhang sect. 3.3) */
static int radial_linear(double k[2], const double K[9],
                         const mv_camera *cams, const mv_calib_view *views,
                         int nviews)
{
    double AtA[4] = { 0, 0, 0, 0 }, Atb[2] = { 0, 0 };
    double u0 = K[2], v0 = K[5];
    double det;
    int vw, i, j;

    for (vw = 0; vw < nviews; vw++) {
        const mv_camera *c = &cams[vw];
        for (i = 0; i < views[vw].n; i++) {
            double X[3], Xc[3], x, y, r2, u, v;
            double row[2], rhs;
            X[0] = views[vw].obj[2 * i];
            X[1] = views[vw].obj[2 * i + 1];
            X[2] = 0.0;
            for (j = 0; j < 3; j++)
                Xc[j] = c->R[j * 3 + 0] * X[0] + c->R[j * 3 + 1] * X[1]
                      + c->R[j * 3 + 2] * X[2] + c->t[j];
            if (Xc[2] <= 1e-12)
                continue;
            x = Xc[0] / Xc[2];
            y = Xc[1] / Xc[2];
            r2 = x * x + y * y;
            u = K[0] * x + K[1] * y + u0;
            v = K[4] * y + v0;
            /* u equation */
            row[0] = (u - u0) * r2;
            row[1] = (u - u0) * r2 * r2;
            rhs = views[vw].img[2 * i] - u;
            AtA[0] += row[0] * row[0];
            AtA[1] += row[0] * row[1];
            AtA[3] += row[1] * row[1];
            Atb[0] += row[0] * rhs;
            Atb[1] += row[1] * rhs;
            /* v equation */
            row[0] = (v - v0) * r2;
            row[1] = (v - v0) * r2 * r2;
            rhs = views[vw].img[2 * i + 1] - v;
            AtA[0] += row[0] * row[0];
            AtA[1] += row[0] * row[1];
            AtA[3] += row[1] * row[1];
            Atb[0] += row[0] * rhs;
            Atb[1] += row[1] * rhs;
        }
    }
    AtA[2] = AtA[1];
    det = AtA[0] * AtA[3] - AtA[1] * AtA[2];
    if (fabs(det) < 1e-300)
        return MV_ERR;
    k[0] = (AtA[3] * Atb[0] - AtA[1] * Atb[1]) / det;
    k[1] = (AtA[0] * Atb[1] - AtA[2] * Atb[0]) / det;
    return MV_OK;
}

/* one pass of homographies + Zhang closed form + extrinsics */
static int calib_once(double K[9], mv_camera *cams,
                      const mv_calib_view *views, int nviews, int zero_skew)
{
    double *Hs, *V;
    double b[6];
    int i, rows, ret = MV_ERR;

    Hs = (double *)malloc((size_t)nviews * 9 * sizeof(double));
    rows = 2 * nviews + (zero_skew ? 1 : 0);
    V = (double *)calloc((size_t)rows * 6, sizeof(double));
    if (!Hs || !V)
        goto done;

    for (i = 0; i < nviews; i++)
        if (mv_homography_dlt(Hs + 9 * i, views[i].obj, views[i].img,
                              views[i].n) != MV_OK)
            goto done;

    /* Zhang constraints: v12^T b = 0 and (v11 - v22)^T b = 0 per view */
    for (i = 0; i < nviews; i++) {
        double v12[6], v11[6], v22[6];
        int j;
        vij(v12, Hs + 9 * i, 0, 1);
        vij(v11, Hs + 9 * i, 0, 0);
        vij(v22, Hs + 9 * i, 1, 1);
        for (j = 0; j < 6; j++) {
            V[(2 * i) * 6 + j] = v12[j];
            V[(2 * i + 1) * 6 + j] = v11[j] - v22[j];
        }
    }
    if (zero_skew)
        V[(2 * nviews) * 6 + 1] = 1.0; /* B12 = 0 <=> zero skew */

    if (mv_nullvec(b, V, rows, 6) != MV_OK)
        goto done;
    if (intrinsics_from_b(K, b) != MV_OK)
        goto done;
    if (zero_skew)
        K[1] = 0.0;

    for (i = 0; i < nviews; i++)
        if (extrinsics_from_H(&cams[i], K, Hs + 9 * i) != MV_OK)
            goto done;
    ret = MV_OK;
done:
    free(Hs);
    free(V);
    return ret;
}

/* remove radial distortion from observed pixels given current K, k1, k2 */
static int undistort_points(double *out, const double *in, int n,
                            const double K[9], const double k[2])
{
    double Kinv[9];
    int i, it;
    if (mv_mat_inv3(Kinv, K) != MV_OK)
        return MV_ERR;
    for (i = 0; i < n; i++) {
        double u = in[2 * i], v = in[2 * i + 1];
        double xn = Kinv[0] * u + Kinv[1] * v + Kinv[2];
        double yn = Kinv[3] * u + Kinv[4] * v + Kinv[5];
        double x = xn, y = yn;
        for (it = 0; it < 20; it++) {
            double r2 = x * x + y * y;
            double rad = 1.0 + r2 * (k[0] + r2 * k[1]);
            if (fabs(rad) < 1e-6)
                break;
            x = xn / rad;
            y = yn / rad;
        }
        out[2 * i] = K[0] * x + K[1] * y + K[2];
        out[2 * i + 1] = K[4] * y + K[5];
    }
    return MV_OK;
}

int mv_calib_planar(double K[9], mv_camera *cams, double k_radial[2],
                    const mv_calib_view *views, int nviews, int zero_skew)
{
    int i;

    if (nviews < 2 || (nviews == 2 && !zero_skew) || !cams)
        return MV_ERR;
    for (i = 0; i < nviews; i++)
        if (views[i].n < 4)
            return MV_ERR;

    if (!k_radial)
        return calib_once(K, cams, views, nviews, zero_skew);

    /* Zhang alternation: fit geometry, estimate k1 k2, undistort the
     * observations, refit; a few rounds converge for moderate distortion */
    {
        mv_calib_view *work;
        double *buf;
        double k[2] = { 0.0, 0.0 };
        int total = 0, off = 0, round, ret = MV_ERR;

        for (i = 0; i < nviews; i++)
            total += views[i].n;
        work = (mv_calib_view *)malloc((size_t)nviews * sizeof(*work));
        buf = (double *)malloc((size_t)total * 2 * sizeof(double));
        if (!work || !buf) {
            free(work);
            free(buf);
            return MV_ERR;
        }
        for (i = 0; i < nviews; i++) {
            work[i] = views[i];
            memcpy(buf + off, views[i].img,
                   (size_t)views[i].n * 2 * sizeof(double));
            work[i].img = buf + off;
            off += 2 * views[i].n;
        }
        for (round = 0; round < 5; round++) {
            if (calib_once(K, cams, work, nviews, zero_skew) != MV_OK)
                goto out;
            /* k against the ORIGINAL (distorted) observations */
            if (radial_linear(k, K, cams, views, nviews) != MV_OK)
                goto out;
            off = 0;
            for (i = 0; i < nviews; i++) {
                if (undistort_points(buf + off, views[i].img, views[i].n,
                                     K, k) != MV_OK)
                    goto out;
                off += 2 * views[i].n;
            }
        }
        k_radial[0] = k[0];
        k_radial[1] = k[1];
        for (i = 0; i < nviews; i++) {
            cams[i].k[0] = k[0];
            cams[i].k[1] = k[1];
        }
        ret = MV_OK;
out:
        free(work);
        free(buf);
        return ret;
    }
}

/* ---- robust view-subset fit (see calib.h) ------------------------ */

static unsigned long long rob_next(unsigned long long *s)
{
    *s = *s * 6364136223846793005ULL + 1442695040888963407ULL;
    return *s >> 33;
}

/* reprojection RMS of one view under candidate K, pose from its own H */
static double view_rms(const double K[9], const double H[9],
                       const mv_calib_view *vw)
{
    mv_camera cam;
    if (extrinsics_from_H(&cam, K, H) != MV_OK)
        return HUGE_VAL;
    return mv_calib_reproj_rms(&cam, vw, 1);
}

static double median_of(double *tmp, const double *x, int n)
{
    int i, j;
    memcpy(tmp, x, (size_t)n * sizeof(double));
    for (i = 1; i < n; i++) { /* insertion sort: n is a view count */
        double v = tmp[i];
        for (j = i; j > 0 && tmp[j - 1] > v; j--)
            tmp[j] = tmp[j - 1];
        tmp[j] = v;
    }
    return (n & 1) ? tmp[n / 2] : 0.5 * (tmp[n / 2 - 1] + tmp[n / 2]);
}

static double score_K(const double K[9], const double *Hs, const int *ok,
                      const mv_calib_view *views, int nviews, double *rms,
                      double *tmp)
{
    int i, nu = 0;
    for (i = 0; i < nviews; i++) {
        rms[i] = ok[i] ? view_rms(K, Hs + 9 * i, &views[i]) : HUGE_VAL;
        if (ok[i])
            tmp[nu++] = rms[i];
    }
    return nu ? median_of(tmp, tmp, nu) : HUGE_VAL;
}

int mv_calib_planar_robust(double K[9], mv_camera *cams, double k_radial[2],
                           const mv_calib_view *views, int nviews,
                           int zero_skew, unsigned char *inlier)
{
    enum { TRIALS = 64, SUBSET = 6 };
    unsigned long long seed = 20260804ULL;
    double *Hs = NULL, *rms = NULL, *tmp = NULL;
    int *ok = NULL, *fitidx = NULL;
    mv_calib_view *fitv = NULL;
    mv_camera *fitcams = NULL;
    double bestK[9], bestmed = HUGE_VAL, thr, kk[2] = { 0.0, 0.0 };
    int bestpick[SUBSET];
    int i, j, t, nusable = 0, nfit = 0, ret = MV_ERR;

    if (nviews < 2 || !cams)
        return MV_ERR;
    if (nviews <= SUBSET) { /* too few views to subsample: plain fit */
        if (mv_calib_planar(K, cams, k_radial, views, nviews,
                            zero_skew) != MV_OK)
            return MV_ERR;
        if (inlier)
            memset(inlier, 1, (size_t)nviews);
        return MV_OK;
    }

    Hs = (double *)malloc((size_t)nviews * 9 * sizeof(double));
    rms = (double *)malloc((size_t)nviews * sizeof(double));
    tmp = (double *)malloc((size_t)nviews * sizeof(double));
    ok = (int *)calloc((size_t)nviews, sizeof(int));
    fitidx = (int *)malloc((size_t)nviews * sizeof(int));
    fitv = (mv_calib_view *)malloc((size_t)nviews * sizeof(*fitv));
    fitcams = (mv_camera *)malloc((size_t)nviews * sizeof(*fitcams));
    if (!Hs || !rms || !tmp || !ok || !fitidx || !fitv || !fitcams)
        goto done;

    for (i = 0; i < nviews; i++) {
        ok[i] = views[i].n >= 4
             && mv_homography_dlt(Hs + 9 * i, views[i].obj, views[i].img,
                                  views[i].n) == MV_OK;
        nusable += ok[i];
    }
    if (nusable <= SUBSET) { /* not enough to subsample either */
        if (mv_calib_planar(K, cams, k_radial, views, nviews,
                            zero_skew) != MV_OK)
            goto done;
        if (inlier)
            memset(inlier, 1, (size_t)nviews);
        ret = MV_OK;
        goto done;
    }

    for (t = 0; t < TRIALS; t++) {
        int pick[SUBSET], np = 0;
        mv_calib_view sub[SUBSET];
        mv_camera subcams[SUBSET];
        double Kt[9], med;
        while (np < SUBSET) {
            int c = (int)(rob_next(&seed) % (unsigned long long)nviews);
            int dup = !ok[c];
            for (j = 0; j < np; j++)
                if (pick[j] == c)
                    dup = 1;
            if (!dup)
                pick[np++] = c;
        }
        for (j = 0; j < SUBSET; j++)
            sub[j] = views[pick[j]];
        if (calib_once(Kt, subcams, sub, SUBSET, zero_skew) != MV_OK)
            continue;
        med = score_K(Kt, Hs, ok, views, nviews, rms, tmp);
        if (med < bestmed) {
            bestmed = med;
            memcpy(bestK, Kt, sizeof(bestK));
            memcpy(bestpick, pick, sizeof(bestpick));
        }
    }
    if (bestmed == HUGE_VAL)
        goto done;

    /* inliers under the best subset's K, then refit on them */
    score_K(bestK, Hs, ok, views, nviews, rms, tmp);
    thr = 3.0 * bestmed + 1e-9;
    for (i = 0; i < nviews; i++)
        if (ok[i] && rms[i] <= thr)
            fitidx[nfit++] = i;
    memcpy(K, bestK, sizeof(bestK));
    if (nfit >= 3) {
        double K2[9], med2;
        for (i = 0; i < nfit; i++)
            fitv[i] = views[fitidx[i]];
        if (mv_calib_planar(K2, fitcams, k_radial ? kk : NULL, fitv, nfit,
                            zero_skew) == MV_OK) {
            med2 = score_K(K2, Hs, ok, views, nviews, rms, tmp);
            if (med2 <= 1.1 * bestmed + 1e-9)
                memcpy(K, K2, sizeof(K2));
            else
                nfit = 0; /* refit worsened things: keep bestK */
        } else {
            nfit = 0;
        }
    } else {
        nfit = 0;
    }
    if (!nfit) { /* fall back to refitting on the winning subset itself */
        for (i = 0; i < SUBSET; i++)
            fitv[i] = views[bestpick[i]];
        if (mv_calib_planar(K, fitcams, k_radial ? kk : NULL, fitv, SUBSET,
                            zero_skew) != MV_OK)
            goto done;
        memcpy(fitidx, bestpick, sizeof(bestpick));
        nfit = SUBSET;
    }

    /* poses for every view: fitted ones from the accepted fit, the rest
     * by decomposing their homography under the final K */
    if (inlier)
        memset(inlier, 0, (size_t)nviews);
    for (i = 0; i < nviews; i++)
        if (!ok[i] || extrinsics_from_H(&cams[i], K, Hs + 9 * i) != MV_OK)
            memset(&cams[i], 0, sizeof(cams[i]));
    for (i = 0; i < nfit; i++) {
        cams[fitidx[i]] = fitcams[i];
        if (inlier)
            inlier[fitidx[i]] = 1;
    }
    for (i = 0; i < nviews; i++) {
        cams[i].k[0] = kk[0];
        cams[i].k[1] = kk[1];
    }
    if (k_radial) {
        k_radial[0] = kk[0];
        k_radial[1] = kk[1];
    }
    ret = MV_OK;
done:
    free(Hs);
    free(rms);
    free(tmp);
    free(ok);
    free(fitidx);
    free(fitv);
    free(fitcams);
    return ret;
}

double mv_calib_reproj_rms(const mv_camera *cams, const mv_calib_view *views,
                           int nviews)
{
    double se = 0.0;
    int vw, i, n = 0;
    for (vw = 0; vw < nviews; vw++) {
        for (i = 0; i < views[vw].n; i++) {
            double X[3], p[2], du, dv;
            X[0] = views[vw].obj[2 * i];
            X[1] = views[vw].obj[2 * i + 1];
            X[2] = 0.0;
            if (mv_cam_project(p, &cams[vw], X) != MV_OK)
                continue;
            du = p[0] - views[vw].img[2 * i];
            dv = p[1] - views[vw].img[2 * i + 1];
            se += du * du + dv * dv;
            n++;
        }
    }
    if (!n)
        return HUGE_VAL;
    return sqrt(se / n);
}
