#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "mv/core.h"
#include "mv/mat.h"
#include "mv/cam.h"
#include "mv/feat.h"
#include "mv/triangulate.h"
#include "mv/photo.h"

int mv_resect_dlt(double P[12], const double *X3, const double *uv, int n)
{
    double c3[3] = { 0, 0, 0 }, c2[2] = { 0, 0 }, s3 = 0.0, s2 = 0.0;
    double *A, p[12];
    int i, j;

    if (n < 6)
        return MV_ERR;
    for (i = 0; i < n; i++) {
        for (j = 0; j < 3; j++)
            c3[j] += X3[3 * i + j];
        c2[0] += uv[2 * i];
        c2[1] += uv[2 * i + 1];
    }
    for (j = 0; j < 3; j++)
        c3[j] /= n;
    c2[0] /= n;
    c2[1] /= n;
    for (i = 0; i < n; i++) {
        double d3 = 0.0, d2 = 0.0;
        for (j = 0; j < 3; j++)
            d3 += (X3[3 * i + j] - c3[j]) * (X3[3 * i + j] - c3[j]);
        d2 = (uv[2 * i] - c2[0]) * (uv[2 * i] - c2[0])
           + (uv[2 * i + 1] - c2[1]) * (uv[2 * i + 1] - c2[1]);
        s3 += sqrt(d3);
        s2 += sqrt(d2);
    }
    s3 = (s3 > 1e-12) ? sqrt(3.0) * n / s3 : 1.0;
    s2 = (s2 > 1e-12) ? sqrt(2.0) * n / s2 : 1.0;

    A = (double *)calloc((size_t)2 * n * 12, sizeof(double));
    if (!A)
        return MV_ERR;
    for (i = 0; i < n; i++) {
        double X = (X3[3 * i] - c3[0]) * s3;
        double Y = (X3[3 * i + 1] - c3[1]) * s3;
        double Z = (X3[3 * i + 2] - c3[2]) * s3;
        double u = (uv[2 * i] - c2[0]) * s2;
        double v = (uv[2 * i + 1] - c2[1]) * s2;
        double *r0 = A + (size_t)(2 * i) * 12;
        double *r1 = A + (size_t)(2 * i + 1) * 12;
        r0[0] = X; r0[1] = Y; r0[2] = Z; r0[3] = 1.0;
        r0[8] = -u * X; r0[9] = -u * Y; r0[10] = -u * Z; r0[11] = -u;
        r1[4] = X; r1[5] = Y; r1[6] = Z; r1[7] = 1.0;
        r1[8] = -v * X; r1[9] = -v * Y; r1[10] = -v * Z; r1[11] = -v;
    }
    if (mv_nullvec(p, A, 2 * n, 12) != MV_OK) {
        free(A);
        return MV_ERR;
    }
    free(A);

    /* denormalize: P = T2^{-1} Pn T3 with
     * T3 = [s3 I | -s3 c3; 0 1], T2^{-1} = [1/s2 I | c2; 0 1] */
    {
        double Pn[12], T3r[16], tmp[12];
        int r, c, k;
        memcpy(Pn, p, sizeof(Pn));
        memset(T3r, 0, sizeof(T3r));
        T3r[0] = T3r[5] = T3r[10] = s3;
        T3r[3] = -s3 * c3[0];
        T3r[7] = -s3 * c3[1];
        T3r[11] = -s3 * c3[2];
        T3r[15] = 1.0;
        for (r = 0; r < 3; r++)
            for (c = 0; c < 4; c++) {
                double a = 0.0;
                for (k = 0; k < 4; k++)
                    a += Pn[r * 4 + k] * T3r[k * 4 + c];
                tmp[r * 4 + c] = a;
            }
        for (c = 0; c < 4; c++) {
            P[0 * 4 + c] = tmp[0 * 4 + c] / s2 + c2[0] * tmp[2 * 4 + c];
            P[1 * 4 + c] = tmp[1 * 4 + c] / s2 + c2[1] * tmp[2 * 4 + c];
            P[2 * 4 + c] = tmp[2 * 4 + c];
        }
    }
    return MV_OK;
}

/* RQ decomposition of the 3x3 left block by Givens rotations from the
 * right: zero m21 (using columns y,z), then m20 (x,z), then m10 (x,y).
 * Accumulates Q so that M = K Q with K upper triangular. */
static int rq3(double K[9], double R[9], const double M[9])
{
    double A[9], Q[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    int step;
    memcpy(A, M, sizeof(A));
    for (step = 0; step < 3; step++) {
        /* (row, colA, colB): zero A[row*3+colA] rotating colA into colB */
        static const int tgt[3][3] = {
            { 2, 1, 2 }, { 2, 0, 2 }, { 1, 0, 1 }
        };
        int r = tgt[step][0], ca = tgt[step][1], cb = tgt[step][2];
        double a = A[r * 3 + cb], b = A[r * 3 + ca], n, c, s;
        int i;
        n = sqrt(a * a + b * b);
        if (n < 1e-15)
            continue;
        c = a / n;
        s = b / n;
        for (i = 0; i < 3; i++) {
            double va = A[i * 3 + ca], vb = A[i * 3 + cb];
            A[i * 3 + ca] = c * va - s * vb;
            A[i * 3 + cb] = s * va + c * vb;
            va = Q[ca * 3 + i];
            vb = Q[cb * 3 + i];
            Q[ca * 3 + i] = c * va - s * vb;
            Q[cb * 3 + i] = s * va + c * vb;
        }
    }
    /* force positive diagonal on K, pushing signs into Q */
    {
        int j, i;
        for (j = 0; j < 3; j++)
            if (A[j * 3 + j] < 0.0)
                for (i = 0; i < 3; i++) {
                    A[i * 3 + j] = -A[i * 3 + j];
                    Q[j * 3 + i] = -Q[j * 3 + i];
                }
    }
    memcpy(K, A, sizeof(A));
    memcpy(R, Q, 9 * sizeof(double));
    return (fabs(K[0]) > 1e-12 && fabs(K[4]) > 1e-12
            && fabs(K[8]) > 1e-12) ? MV_OK : MV_ERR;
}

int mv_cam_from_projection(mv_camera *cam, const double P[12])
{
    double M[9], Ps[12], K[9], R[9], det;
    int i;

    for (i = 0; i < 3; i++) {
        M[3 * i] = P[4 * i];
        M[3 * i + 1] = P[4 * i + 1];
        M[3 * i + 2] = P[4 * i + 2];
    }
    det = M[0] * (M[4] * M[8] - M[5] * M[7])
        - M[1] * (M[3] * M[8] - M[5] * M[6])
        + M[2] * (M[3] * M[7] - M[4] * M[6]);
    if (fabs(det) < 1e-15)
        return MV_ERR;
    /* a projective nullvector has arbitrary sign: a camera with det(KR)
     * < 0 is the mirror solution, so flip P */
    for (i = 0; i < 12; i++)
        Ps[i] = (det < 0.0) ? -P[i] : P[i];
    for (i = 0; i < 3; i++) {
        M[3 * i] = Ps[4 * i];
        M[3 * i + 1] = Ps[4 * i + 1];
        M[3 * i + 2] = Ps[4 * i + 2];
    }
    if (rq3(K, R, M) != MV_OK)
        return MV_ERR;
    {
        double kinv[9], p4[3], t[3];
        int j;
        for (i = 0; i < 9; i++)
            K[i] /= K[8];
        if (mv_mat_inv3(kinv, K) != MV_OK)
            return MV_ERR;
        for (i = 0; i < 3; i++)
            p4[i] = Ps[4 * i + 3];
        /* K was rescaled; rescale the translation solve consistently:
         * P ~ K[R|t] holds up to the same global scale as M = K R, so
         * recover the scale from M's decomposition */
        {
            double KR[9], sc = 0.0;
            mv_mat_mul(KR, K, R, 3, 3, 3);
            for (i = 0; i < 9; i++)
                if (fabs(KR[i]) > 1e-9) {
                    sc = M[i] / KR[i];
                    break;
                }
            if (fabs(sc) < 1e-15)
                return MV_ERR;
            for (i = 0; i < 3; i++)
                p4[i] /= sc;
        }
        mv_mat_mul(t, kinv, p4, 3, 3, 1);
        memcpy(cam->K, K, sizeof(K));
        memcpy(cam->R, R, 9 * sizeof(double));
        for (j = 0; j < 3; j++)
            cam->t[j] = t[j];
        memset(cam->k, 0, sizeof(cam->k));
    }
    return MV_OK;
}

static int rcmp_dbl(const void *a, const void *b)
{
    double d = *(const double *)a - *(const double *)b;
    return d < 0 ? -1 : d > 0 ? 1 : 0;
}

/* geometric polish: Gauss-Newton with LM damping over {fx, fy, cx,
 * cy, Rodrigues r, t} minimizing true reprojection error -- the DLT
 * minimizes an algebraic surrogate whose bias shows up mostly as a
 * correlated focal/depth shift (the gold-standard two-step of
 * Hartley-Zisserman: linear init, geometric refine) */
static void rodrigues_R(double R[9], const double r[3])
{
    double th = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
    double kx, ky, kz, c, sn, v;
    if (th < 1e-12) {
        memset(R, 0, 9 * sizeof(double));
        R[0] = R[4] = R[8] = 1.0;
        R[1] = -r[2]; R[2] = r[1];
        R[3] = r[2];  R[5] = -r[0];
        R[6] = -r[1]; R[7] = r[0];
        return;
    }
    kx = r[0] / th; ky = r[1] / th; kz = r[2] / th;
    c = cos(th); sn = sin(th); v = 1.0 - c;
    R[0] = c + kx * kx * v;      R[1] = kx * ky * v - kz * sn;
    R[2] = kx * kz * v + ky * sn;
    R[3] = ky * kx * v + kz * sn; R[4] = c + ky * ky * v;
    R[5] = ky * kz * v - kx * sn;
    R[6] = kz * kx * v - ky * sn; R[7] = kz * ky * v + kx * sn;
    R[8] = c + kz * kz * v;
}

static void R_rodrigues(double r[3], const double R[9])
{
    double c = (R[0] + R[4] + R[8] - 1.0) / 2.0, th, s;
    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;
    th = acos(c);
    s = sin(th);
    if (fabs(s) < 1e-9) {
        /* theta ~ 0: no rotation. theta ~ pi: sin -> 0 but the rotation
         * is a real half-turn -- recover the axis from the diagonal
         * (R = 2 a a^T - I there), sign from the off-diagonal, rather
         * than collapsing it to zero. */
        if (c > 0.0) {
            r[0] = r[1] = r[2] = 0.0;
        } else {
            double ax = sqrt((R[0] + 1.0) > 0.0 ? (R[0] + 1.0) / 2.0 : 0),
                   ay = sqrt((R[4] + 1.0) > 0.0 ? (R[4] + 1.0) / 2.0 : 0),
                   az = sqrt((R[8] + 1.0) > 0.0 ? (R[8] + 1.0) / 2.0 : 0);
            double PI = 3.14159265358979323846;
            if (R[1] + R[3] < 0.0) ay = -ay;
            if (R[2] + R[6] < 0.0) az = -az;
            /* pin the largest-magnitude axis positive for a stable sign */
            if (ax >= ay && ax >= az) {
                if (R[1] + R[3] < 0.0 && ay > 0.0) ay = -ay;
            }
            r[0] = PI * ax;
            r[1] = PI * ay;
            r[2] = PI * az;
        }
        return;
    }
    r[0] = th * (R[7] - R[5]) / (2.0 * s);
    r[1] = th * (R[2] - R[6]) / (2.0 * s);
    r[2] = th * (R[3] - R[1]) / (2.0 * s);
}

static void resect_residuals(const double q[10], const double *X3,
                             const double *uv, int n, double *r)
{
    mv_camera c;
    int i;
    mv_cam_set_K(&c, q[0], q[1], q[2], q[3]);
    rodrigues_R(c.R, q + 4);
    c.t[0] = q[7]; c.t[1] = q[8]; c.t[2] = q[9];
    memset(c.k, 0, sizeof(c.k));
    for (i = 0; i < n; i++) {
        double p[2];
        if (mv_cam_project(p, &c, X3 + 3 * i) != MV_OK) {
            r[2 * i] = 1e6;
            r[2 * i + 1] = 1e6;
            continue;
        }
        r[2 * i] = p[0] - uv[2 * i];
        r[2 * i + 1] = p[1] - uv[2 * i + 1];
    }
}

/* solve A x = b for a small symmetric system by Gaussian elimination
 * with partial pivoting; A is m x m row-major, destroyed */
static int solve_small(double *A, double *b, int m)
{
    int i, j, k;
    for (i = 0; i < m; i++) {
        int piv = i;
        double t;
        for (j = i + 1; j < m; j++)
            if (fabs(A[j * m + i]) > fabs(A[piv * m + i]))
                piv = j;
        if (fabs(A[piv * m + i]) < 1e-15)
            return MV_ERR;
        if (piv != i) {
            for (k = 0; k < m; k++) {
                t = A[i * m + k];
                A[i * m + k] = A[piv * m + k];
                A[piv * m + k] = t;
            }
            t = b[i];
            b[i] = b[piv];
            b[piv] = t;
        }
        for (j = i + 1; j < m; j++) {
            double f = A[j * m + i] / A[i * m + i];
            for (k = i; k < m; k++)
                A[j * m + k] -= f * A[i * m + k];
            b[j] -= f * b[i];
        }
    }
    for (i = m - 1; i >= 0; i--) {
        double a = b[i];
        for (k = i + 1; k < m; k++)
            a -= A[i * m + k] * b[k];
        b[i] = a / A[i * m + i];
    }
    return MV_OK;
}

/* full LM over the free subset of {fx,fy,cx,cy,r,t}: numeric Jacobian,
 * normal equations, multiplicative damping -- the couplings between
 * focal, depth, and rotation demand the full system */
static void resect_polish_m(mv_camera *cam, const double *X3,
                            const double *uv, int n,
                            const unsigned char freem[10])
{
    double q[10], lam = 1e-3, cost = 0.0;
    double *r0, *r1, *J;
    int map[10], nf = 0, it, i, a;

    for (a = 0; a < 10; a++)
        if (freem[a])
            map[nf++] = a;
    /* nf is in [1,10] by construction; state it so the -O2 fortify
     * range analysis can bound the fixed A[100]/g[10] copies below */
    if (nf < 1 || nf > 10 || n < 6)
        return;
    r0 = (double *)malloc((size_t)2 * n * sizeof(double));
    r1 = (double *)malloc((size_t)2 * n * sizeof(double));
    J = (double *)malloc((size_t)2 * n * nf * sizeof(double));
    if (!r0 || !r1 || !J) {
        free(r0);
        free(r1);
        free(J);
        return;
    }
    q[0] = cam->K[0]; q[1] = cam->K[4];
    q[2] = cam->K[2]; q[3] = cam->K[5];
    R_rodrigues(q + 4, cam->R);
    q[7] = cam->t[0]; q[8] = cam->t[1]; q[9] = cam->t[2];
    resect_residuals(q, X3, uv, n, r0);
    for (i = 0; i < 2 * n; i++)
        cost += r0[i] * r0[i];
    for (it = 0; it < 100; it++) {
        double A[100] = { 0 }, g[10] = { 0 }, qn[10], cn;
        int b, c;
        for (a = 0; a < nf; a++) {
            double step = 1e-6 * (1.0 + fabs(q[map[a]]));
            memcpy(qn, q, sizeof(qn));
            qn[map[a]] = q[map[a]] + step;
            resect_residuals(qn, X3, uv, n, r1);
            for (i = 0; i < 2 * n; i++)
                J[(size_t)i * nf + a] = (r1[i] - r0[i]) / step;
        }
        for (b = 0; b < nf; b++) {
            g[b] = 0.0;
            for (i = 0; i < 2 * n; i++)
                g[b] += J[(size_t)i * nf + b] * r0[i];
            for (c = 0; c < nf; c++) {
                double h = 0.0;
                for (i = 0; i < 2 * n; i++)
                    h += J[(size_t)i * nf + b] * J[(size_t)i * nf + c];
                A[b * nf + c] = h;
            }
        }
        {
            double Ad[100], gd[10];
            memcpy(Ad, A, (size_t)nf * nf * sizeof(double));
            memcpy(gd, g, (size_t)nf * sizeof(double));
            for (b = 0; b < nf; b++)
                Ad[b * nf + b] *= 1.0 + lam;
            if (solve_small(Ad, gd, nf) != MV_OK)
                break;
            memcpy(qn, q, sizeof(qn));
            for (b = 0; b < nf; b++)
                qn[map[b]] -= gd[b];
        }
        resect_residuals(qn, X3, uv, n, r1);
        cn = 0.0;
        for (i = 0; i < 2 * n; i++)
            cn += r1[i] * r1[i];
        if (cn < cost) {
            int done = (cost - cn) < 1e-12 * (1.0 + cost);
            memcpy(q, qn, sizeof(q));
            memcpy(r0, r1, (size_t)2 * n * sizeof(double));
            cost = cn;
            lam *= 0.3;
            if (lam < 1e-12)
                lam = 1e-12;
            if (done)
                break;
        } else {
            lam *= 5.0;
            if (lam > 1e10)
                break;
        }
    }
    free(r0);
    free(r1);
    free(J);
    mv_cam_set_K(cam, q[0], q[1], q[2], q[3]);
    rodrigues_R(cam->R, q + 4);
    cam->t[0] = q[7]; cam->t[1] = q[8]; cam->t[2] = q[9];
    memset(cam->k, 0, sizeof(cam->k));
}

static void resect_polish(mv_camera *cam, const double *X3,
                          const double *uv, int n)
{
    static const unsigned char all_free[10] = { 1, 1, 1, 1, 1, 1, 1, 1,
                                                1, 1 };
    resect_polish_m(cam, X3, uv, n, all_free);
}

int mv_resect_robust(mv_camera *cam, double *X3, double *uv, int *n)
{
    double P[12];
    int round, i, m;

    /* RANSAC seeding [FischlerBolles1981]: descriptor matching hands us
     * correspondence sets with genuine gross outliers, and a DLT fit on
     * the full set can be wrong enough that a median gate keeps the
     * garbage.  Minimal 6-point hypotheses, deterministic sampling,
     * 2.5 px inlier gate; the winner's inlier set then goes through the
     * usual fit-trim-refit polish below. */
    if (*n >= 12) {
        unsigned long long seed = 20260806ULL;
        int best = -1, trial;
        static unsigned char bestin[8192];
        unsigned char *in;
        if (*n > 8192)
            *n = 8192;
        in = (unsigned char *)malloc((size_t)*n);
        if (!in)
            return MV_ERR;
        for (trial = 0; trial < 256; trial++) {
            double Xs[18], us[12], Pt[12];
            mv_camera ct;
            int pick[6], np = 0, k, votes = 0;
            while (np < 6) {
                int c, dup = 0;
                seed = seed * 6364136223846793005ULL
                       + 1442695040888963407ULL;
                c = (int)((seed >> 33) % (unsigned long long)*n);
                for (k = 0; k < np; k++)
                    if (pick[k] == c)
                        dup = 1;
                if (!dup)
                    pick[np++] = c;
            }
            for (k = 0; k < 6; k++) {
                memcpy(Xs + 3 * k, X3 + 3 * pick[k],
                       3 * sizeof(double));
                memcpy(us + 2 * k, uv + 2 * pick[k],
                       2 * sizeof(double));
            }
            if (mv_resect_dlt(Pt, Xs, us, 6) != MV_OK
                || mv_cam_from_projection(&ct, Pt) != MV_OK)
                continue;
            for (i = 0; i < *n; i++) {
                double p[2], du, dv;
                in[i] = 0;
                if (mv_cam_project(p, &ct, X3 + 3 * i) != MV_OK)
                    continue;
                du = p[0] - uv[2 * i];
                dv = p[1] - uv[2 * i + 1];
                if (du * du + dv * dv < 6.25) {
                    in[i] = 1;
                    votes++;
                }
            }
            if (votes > best) {
                best = votes;
                memcpy(bestin, in, (size_t)*n);
            }
        }
        free(in);
        if (best < 6)
            return MV_ERR;
        m = 0;
        for (i = 0; i < *n; i++)
            if (bestin[i]) {
                memmove(X3 + 3 * m, X3 + 3 * i, 3 * sizeof(double));
                memmove(uv + 2 * m, uv + 2 * i, 2 * sizeof(double));
                m++;
            }
        *n = m;
    }

    if (mv_resect_dlt(P, X3, uv, *n) != MV_OK
        || mv_cam_from_projection(cam, P) != MV_OK)
        return MV_ERR;
    for (round = 0; round < 2; round++) {
        double *resid, medr;
        resid = (double *)malloc((size_t)*n * 2 * sizeof(double));
        if (!resid)
            return MV_ERR;
        for (i = 0; i < *n; i++) {
            double p[2], du, dv;
            if (mv_cam_project(p, cam, X3 + 3 * i) != MV_OK) {
                resid[i] = 1e9;
            } else {
                du = p[0] - uv[2 * i];
                dv = p[1] - uv[2 * i + 1];
                resid[i] = sqrt(du * du + dv * dv);
            }
            resid[*n + i] = resid[i];
        }
        qsort(resid + *n, (size_t)*n, sizeof(double), rcmp_dbl);
        medr = resid[*n + *n / 2];
        m = 0;
        for (i = 0; i < *n; i++)
            if (resid[i] <= 3.0 * medr + 1e-9) {
                memmove(X3 + 3 * m, X3 + 3 * i, 3 * sizeof(double));
                memmove(uv + 2 * m, uv + 2 * i, 2 * sizeof(double));
                m++;
            }
        free(resid);
        if (m < 6)
            return MV_ERR;
        *n = m;
        if (mv_resect_dlt(P, X3, uv, *n) != MV_OK
            || mv_cam_from_projection(cam, P) != MV_OK)
            return MV_ERR;
    }
    resect_polish(cam, X3, uv, *n);
    return MV_OK;
}

#define PH_MAXF 500

int mv_photo_register(mv_camera *photo, const mv_camera *c1,
                      const mv_camera *c2, const unsigned char *img1,
                      const unsigned char *img2, const unsigned char *imgp,
                      int w, int h, int wp, int hp, int *nanchor,
                      int *ninl)
{
    mv_feature *f1, *f2, *fp, *fa;
    double *ax; /* anchor 3D points, parallel to fa */
    int *idx;
    int n1, n2, np, na = 0, i, ret = MV_ERR;

    f1 = (mv_feature *)malloc(PH_MAXF * sizeof(mv_feature));
    f2 = (mv_feature *)malloc(PH_MAXF * sizeof(mv_feature));
    fp = (mv_feature *)malloc(PH_MAXF * sizeof(mv_feature));
    fa = (mv_feature *)malloc(PH_MAXF * sizeof(mv_feature));
    ax = (double *)malloc((size_t)PH_MAXF * 3 * sizeof(double));
    idx = (int *)malloc(PH_MAXF * sizeof(int));
    if (!f1 || !f2 || !fp || !fa || !ax || !idx)
        goto done;

    n1 = mv_feat_detect(f1, PH_MAXF, img1, w, h);
    n2 = mv_feat_detect(f2, PH_MAXF, img2, w, h);
    np = mv_feat_detect(fp, PH_MAXF, imgp, wp, hp);
    if (n1 < 20 || n2 < 20 || np < 20)
        goto done;

    /* rig pair -> triangulated anchor bank; each anchor keeps camera
     * 1's descriptor.  Gate on positive depth in both views and a
     * 2 px reprojection sanity bound. */
    if (mv_feat_match(idx, f1, n1, f2, n2, 0.6) < 8)
        goto done;
    for (i = 0; i < n1 && na < PH_MAXF; i++) {
        const mv_camera *cams[2];
        double uvs[4], X[3], p1[2], p2[2], e1, e2, z1, z2;
        if (idx[i] < 0)
            continue;
        cams[0] = c1;
        cams[1] = c2;
        uvs[0] = f1[i].u;
        uvs[1] = f1[i].v;
        uvs[2] = f2[idx[i]].u;
        uvs[3] = f2[idx[i]].v;
        if (mv_triangulate(X, cams, uvs, 2) != MV_OK)
            continue;
        z1 = c1->R[6] * X[0] + c1->R[7] * X[1] + c1->R[8] * X[2]
             + c1->t[2];
        z2 = c2->R[6] * X[0] + c2->R[7] * X[1] + c2->R[8] * X[2]
             + c2->t[2];
        if (z1 <= 0.0 || z2 <= 0.0)
            continue;
        if (mv_cam_project(p1, c1, X) != MV_OK
            || mv_cam_project(p2, c2, X) != MV_OK)
            continue;
        e1 = fabs(p1[0] - uvs[0]) + fabs(p1[1] - uvs[1]);
        e2 = fabs(p2[0] - uvs[2]) + fabs(p2[1] - uvs[3]);
        if (e1 > 2.0 || e2 > 2.0)
            continue;
        fa[na] = f1[i];
        ax[3 * na] = X[0];
        ax[3 * na + 1] = X[1];
        ax[3 * na + 2] = X[2];
        na++;
    }
    if (nanchor)
        *nanchor = na;
    if (na < 10)
        goto done;

    /* photo -> anchor bank -> robust resection */
    {
        double *X3, *uv;
        int nc = 0;
        if (mv_feat_match(idx, fp, np, fa, na, 0.6) < 6)
            goto done;
        X3 = (double *)malloc((size_t)np * 3 * sizeof(double));
        uv = (double *)malloc((size_t)np * 2 * sizeof(double));
        if (!X3 || !uv) {
            free(X3);
            free(uv);
            goto done;
        }
        for (i = 0; i < np; i++) {
            if (idx[i] < 0)
                continue;
            memcpy(X3 + 3 * nc, ax + 3 * idx[i], 3 * sizeof(double));
            uv[2 * nc] = fp[i].u;
            uv[2 * nc + 1] = fp[i].v;
            nc++;
        }
        if (nc >= 6 && mv_resect_robust(photo, X3, uv, &nc) == MV_OK) {
            /* casual-photo conditioning: the principal point is nearly
             * unobservable from one resection over a shallow scene (it
             * trades against lateral pose), so pin it at the image
             * center -- the standard SfM practice [Snavely2006] -- and
             * re-polish focal + pose under that constraint */
            static const unsigned char pinned[10] = { 1, 1, 0, 0, 1, 1,
                                                      1, 1, 1, 1 };
            photo->K[2] = wp / 2.0;
            photo->K[5] = hp / 2.0;
            resect_polish_m(photo, X3, uv, nc, pinned);
            if (ninl)
                *ninl = nc;
            ret = MV_OK;
        }
        free(X3);
        free(uv);
    }
done:
    free(f1);
    free(f2);
    free(fp);
    free(fa);
    free(ax);
    free(idx);
    return ret;
}
