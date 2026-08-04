#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "mv/core.h"
#include "mv/mat.h"
#include "mv/epipolar.h"
#include "mv/triangulate.h"

/* Unit Frobenius norm, canonical sign: largest-|.| entry positive. */
static void canonicalize(double F[9])
{
    int i, mx = 0;
    double nm = mv_norm(F, 9);
    if (nm < 1e-300)
        return;
    for (i = 0; i < 9; i++)
        F[i] /= nm;
    for (i = 1; i < 9; i++)
        if (fabs(F[i]) > fabs(F[mx]))
            mx = i;
    if (F[mx] < 0.0)
        for (i = 0; i < 9; i++)
            F[i] = -F[i];
}

void mv_fundamental_from_cams(double F[9], const mv_camera *c1,
                              const mv_camera *c2)
{
    double P1[12], P2[12], P1t[12], PPt[9], PPtinv[9], P1plus[12];
    double C1[3], C1h[4], e2[3], E2x[9], M[9];

    mv_cam_projection(P1, c1);
    mv_cam_projection(P2, c2);

    /* P1^+ = P1^T (P1 P1^T)^-1 (4x3) */
    mv_mat_transpose(P1t, P1, 3, 4);
    mv_mat_mul(PPt, P1, P1t, 3, 4, 3);
    mv_mat_inv3(PPtinv, PPt);
    mv_mat_mul(P1plus, P1t, PPtinv, 4, 3, 3);

    /* epipole in image 2: e2 = P2 [C1; 1] */
    mv_cam_center(C1, c1);
    C1h[0] = C1[0]; C1h[1] = C1[1]; C1h[2] = C1[2]; C1h[3] = 1.0;
    mv_mat_mul(e2, P2, C1h, 3, 4, 1);

    /* F = [e2]_x P2 P1^+ */
    mv_skew3(E2x, e2);
    mv_mat_mul(M, P2, P1plus, 3, 4, 3);
    mv_mat_mul(F, E2x, M, 3, 3, 3);
    canonicalize(F);
}

int mv_fundamental_8point(double F[9], const double *uv1, const double *uv2,
                          int n)
{
    double c1x = 0.0, c1y = 0.0, c2x = 0.0, c2y = 0.0, s1 = 0.0, s2 = 0.0;
    double T1[9], T2[9], T2t[9], f[9], Fn[9], U[9], S[3], V[9], Vt[9];
    double US[9], F2[9], tmp[9];
    double *A;
    int i;

    if (n < 8)
        return MV_ERR;

    for (i = 0; i < n; i++) {
        c1x += uv1[2 * i]; c1y += uv1[2 * i + 1];
        c2x += uv2[2 * i]; c2y += uv2[2 * i + 1];
    }
    c1x /= n; c1y /= n; c2x /= n; c2y /= n;
    for (i = 0; i < n; i++) {
        double dx1 = uv1[2 * i] - c1x, dy1 = uv1[2 * i + 1] - c1y;
        double dx2 = uv2[2 * i] - c2x, dy2 = uv2[2 * i + 1] - c2y;
        s1 += sqrt(dx1 * dx1 + dy1 * dy1);
        s2 += sqrt(dx2 * dx2 + dy2 * dy2);
    }
    if (s1 < 1e-12 || s2 < 1e-12)
        return MV_ERR;
    s1 = sqrt(2.0) * n / s1; /* mean distance -> sqrt(2) */
    s2 = sqrt(2.0) * n / s2;

    memset(T1, 0, sizeof(T1));
    T1[0] = s1; T1[2] = -s1 * c1x; T1[4] = s1; T1[5] = -s1 * c1y; T1[8] = 1.0;
    memset(T2, 0, sizeof(T2));
    T2[0] = s2; T2[2] = -s2 * c2x; T2[4] = s2; T2[5] = -s2 * c2y; T2[8] = 1.0;

    A = (double *)malloc((size_t)n * 9 * sizeof(double));
    if (!A)
        return MV_ERR;
    for (i = 0; i < n; i++) {
        double x1 = s1 * (uv1[2 * i] - c1x);
        double y1 = s1 * (uv1[2 * i + 1] - c1y);
        double x2 = s2 * (uv2[2 * i] - c2x);
        double y2 = s2 * (uv2[2 * i + 1] - c2y);
        double *row = A + 9 * i;
        row[0] = x2 * x1; row[1] = x2 * y1; row[2] = x2;
        row[3] = y2 * x1; row[4] = y2 * y1; row[5] = y2;
        row[6] = x1;      row[7] = y1;      row[8] = 1.0;
    }
    i = mv_nullvec(f, A, n, 9);
    free(A);
    if (i != MV_OK)
        return MV_ERR;
    memcpy(Fn, f, sizeof(Fn));

    /* enforce rank 2: zero the smallest singular value */
    memcpy(U, Fn, sizeof(U));
    if (mv_svd(U, S, V, 3, 3) != MV_OK)
        return MV_ERR;
    S[2] = 0.0;
    for (i = 0; i < 3; i++) {
        US[i * 3 + 0] = U[i * 3 + 0] * S[0];
        US[i * 3 + 1] = U[i * 3 + 1] * S[1];
        US[i * 3 + 2] = 0.0;
    }
    mv_mat_transpose(Vt, V, 3, 3);
    mv_mat_mul(F2, US, Vt, 3, 3, 3);

    /* denormalize: F = T2^T F2 T1 */
    mv_mat_transpose(T2t, T2, 3, 3);
    mv_mat_mul(tmp, T2t, F2, 3, 3, 3);
    mv_mat_mul(F, tmp, T1, 3, 3, 3);
    canonicalize(F);
    return MV_OK;
}

void mv_essential_from_fundamental(double E[9], const double F[9],
                                   const double K1[9], const double K2[9])
{
    double K2t[9], tmp[9];
    mv_mat_transpose(K2t, K2, 3, 3);
    mv_mat_mul(tmp, K2t, F, 3, 3, 3);
    mv_mat_mul(E, tmp, K1, 3, 3, 3);
}

double mv_sym_epipolar_dist(const double F[9], const double uv1[2],
                            const double uv2[2])
{
    double x1[3], x2[3], l2[3], l1[3], v, n1, n2;
    int i;
    x1[0] = uv1[0]; x1[1] = uv1[1]; x1[2] = 1.0;
    x2[0] = uv2[0]; x2[1] = uv2[1]; x2[2] = 1.0;
    for (i = 0; i < 3; i++) {
        l2[i] = F[i * 3 + 0] * x1[0] + F[i * 3 + 1] * x1[1] + F[i * 3 + 2];
        l1[i] = F[0 * 3 + i] * x2[0] + F[1 * 3 + i] * x2[1] + F[2 * 3 + i];
    }
    v = x2[0] * l2[0] + x2[1] * l2[1] + l2[2];
    n2 = sqrt(l2[0] * l2[0] + l2[1] * l2[1]);
    n1 = sqrt(l1[0] * l1[0] + l1[1] * l1[1]);
    if (n1 < 1e-300 || n2 < 1e-300)
        return HUGE_VAL;
    return 0.5 * (fabs(v) / n1 + fabs(v) / n2);
}

int mv_essential_pose(double R[9], double t[3], const double E[9],
                      const double *x1, const double *x2, int n)
{
    static const double Wm[9] = { 0, -1, 0, 1, 0, 0, 0, 0, 1 };
    double U[9], S[3], V[9], Vt[9], Wt[9], tmp[9];
    double Rc[2][9], tc[2][3];
    int best = -1, bestvotes = -1, ri, ti, i, j;

    memcpy(U, E, sizeof(U));
    if (mv_svd(U, S, V, 3, 3) != MV_OK)
        return MV_ERR;
    /* force proper rotations out of the factors */
    {
        double dU = U[0] * (U[4] * U[8] - U[5] * U[7])
                  - U[1] * (U[3] * U[8] - U[5] * U[6])
                  + U[2] * (U[3] * U[7] - U[4] * U[6]);
        double dV = V[0] * (V[4] * V[8] - V[5] * V[7])
                  - V[1] * (V[3] * V[8] - V[5] * V[6])
                  + V[2] * (V[3] * V[7] - V[4] * V[6]);
        if (dU < 0.0)
            for (i = 0; i < 3; i++)
                U[i * 3 + 2] = -U[i * 3 + 2];
        if (dV < 0.0)
            for (i = 0; i < 3; i++)
                V[i * 3 + 2] = -V[i * 3 + 2];
    }
    mv_mat_transpose(Vt, V, 3, 3);
    mv_mat_transpose(Wt, Wm, 3, 3);
    mv_mat_mul(tmp, U, Wm, 3, 3, 3);
    mv_mat_mul(Rc[0], tmp, Vt, 3, 3, 3);
    mv_mat_mul(tmp, U, Wt, 3, 3, 3);
    mv_mat_mul(Rc[1], tmp, Vt, 3, 3, 3);
    for (i = 0; i < 3; i++) {
        tc[0][i] = U[i * 3 + 2];
        tc[1][i] = -U[i * 3 + 2];
    }

    for (ri = 0; ri < 2; ri++)
        for (ti = 0; ti < 2; ti++) {
            mv_camera c1, c2;
            const mv_camera *cams[2];
            int votes = 0;
            mv_cam_set_K(&c1, 1.0, 1.0, 0.0, 0.0);
            mv_cam_set_identity_pose(&c1);
            memset(c1.k, 0, sizeof(c1.k));
            c2 = c1;
            memcpy(c2.R, Rc[ri], sizeof(c2.R));
            for (j = 0; j < 3; j++)
                c2.t[j] = tc[ti][j];
            cams[0] = &c1;
            cams[1] = &c2;
            for (i = 0; i < n; i++) {
                double X[3], uvs[4], z2;
                uvs[0] = x1[2 * i];
                uvs[1] = x1[2 * i + 1];
                uvs[2] = x2[2 * i];
                uvs[3] = x2[2 * i + 1];
                if (mv_triangulate(X, cams, uvs, 2) != MV_OK)
                    continue;
                z2 = c2.R[6] * X[0] + c2.R[7] * X[1] + c2.R[8] * X[2]
                     + c2.t[2];
                if (X[2] > 0.0 && z2 > 0.0)
                    votes++;
            }
            if (votes > bestvotes) {
                bestvotes = votes;
                best = 2 * ri + ti;
            }
        }
    if (best < 0 || bestvotes < n / 2 || bestvotes < 2)
        return MV_ERR;
    memcpy(R, Rc[best / 2], 9 * sizeof(double));
    memcpy(t, tc[best % 2], 3 * sizeof(double));
    return MV_OK;
}

static int fcmp_dbl(const void *a, const void *b)
{
    double d = *(const double *)a - *(const double *)b;
    return d < 0 ? -1 : d > 0 ? 1 : 0;
}

int mv_fundamental_8point_robust(double F[9], double *uv1, double *uv2,
                                 int *n)
{
    int round, i, m;
    if (mv_fundamental_8point(F, uv1, uv2, *n) != MV_OK)
        return MV_ERR;
    for (round = 0; round < 2; round++) {
        double *resid, medr;
        resid = (double *)malloc((size_t)*n * 2 * sizeof(double));
        if (!resid)
            return MV_ERR;
        for (i = 0; i < *n; i++) {
            resid[i] = mv_sym_epipolar_dist(F, uv1 + 2 * i, uv2 + 2 * i);
            resid[*n + i] = resid[i];
        }
        qsort(resid + *n, (size_t)*n, sizeof(double), fcmp_dbl);
        medr = resid[*n + *n / 2];
        m = 0;
        for (i = 0; i < *n; i++)
            if (resid[i] <= 3.0 * medr + 1e-9) {
                uv1[2 * m] = uv1[2 * i];
                uv1[2 * m + 1] = uv1[2 * i + 1];
                uv2[2 * m] = uv2[2 * i];
                uv2[2 * m + 1] = uv2[2 * i + 1];
                m++;
            }
        free(resid);
        if (m < 8)
            return MV_ERR;
        *n = m;
        if (mv_fundamental_8point(F, uv1, uv2, *n) != MV_OK)
            return MV_ERR;
    }
    return MV_OK;
}
