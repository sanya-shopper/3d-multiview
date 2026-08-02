#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "mv/core.h"
#include "mv/mat.h"

void mv_mat_mul(double *C, const double *A, const double *B,
                int m, int n, int p)
{
    int i, j, k;
    for (i = 0; i < m; i++) {
        for (j = 0; j < p; j++) {
            double s = 0.0;
            for (k = 0; k < n; k++)
                s += A[i * n + k] * B[k * p + j];
            C[i * p + j] = s;
        }
    }
}

void mv_mat_transpose(double *T, const double *A, int m, int n)
{
    int i, j;
    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            T[j * m + i] = A[i * n + j];
}

int mv_mat_inv3(double inv[9], const double A[9])
{
    double det = A[0] * (A[4] * A[8] - A[5] * A[7])
               - A[1] * (A[3] * A[8] - A[5] * A[6])
               + A[2] * (A[3] * A[7] - A[4] * A[6]);
    double id;
    if (fabs(det) < 1e-300)
        return MV_ERR;
    id = 1.0 / det;
    inv[0] = (A[4] * A[8] - A[5] * A[7]) * id;
    inv[1] = (A[2] * A[7] - A[1] * A[8]) * id;
    inv[2] = (A[1] * A[5] - A[2] * A[4]) * id;
    inv[3] = (A[5] * A[6] - A[3] * A[8]) * id;
    inv[4] = (A[0] * A[8] - A[2] * A[6]) * id;
    inv[5] = (A[2] * A[3] - A[0] * A[5]) * id;
    inv[6] = (A[3] * A[7] - A[4] * A[6]) * id;
    inv[7] = (A[1] * A[6] - A[0] * A[7]) * id;
    inv[8] = (A[0] * A[4] - A[1] * A[3]) * id;
    return MV_OK;
}

void mv_skew3(double M[9], const double v[3])
{
    M[0] = 0.0;   M[1] = -v[2]; M[2] = v[1];
    M[3] = v[2];  M[4] = 0.0;   M[5] = -v[0];
    M[6] = -v[1]; M[7] = v[0];  M[8] = 0.0;
}

void mv_cross3(double c[3], const double a[3], const double b[3])
{
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

double mv_dot(const double *a, const double *b, int n)
{
    double s = 0.0;
    int i;
    for (i = 0; i < n; i++)
        s += a[i] * b[i];
    return s;
}

double mv_norm(const double *a, int n)
{
    return sqrt(mv_dot(a, a, n));
}

int mv_normalize(double *a, int n)
{
    double nm = mv_norm(a, n);
    int i;
    if (nm < 1e-300)
        return MV_ERR;
    for (i = 0; i < n; i++)
        a[i] /= nm;
    return MV_OK;
}

int mv_svd(double *A, double *S, double *V, int m, int n)
{
    const double eps = 1e-15;
    int i, j, p, q, k, sweep;

    if (m < n)
        return MV_ERR;

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            V[i * n + j] = (i == j) ? 1.0 : 0.0;

    for (sweep = 0; sweep < 100; sweep++) {
        int rotated = 0;
        for (p = 0; p < n - 1; p++) {
            for (q = p + 1; q < n; q++) {
                double alpha = 0.0, beta = 0.0, gamma = 0.0;
                double zeta, tt, c, s;
                for (k = 0; k < m; k++) {
                    alpha += A[k * n + p] * A[k * n + p];
                    beta  += A[k * n + q] * A[k * n + q];
                    gamma += A[k * n + p] * A[k * n + q];
                }
                if (gamma == 0.0 || fabs(gamma) <= eps * sqrt(alpha * beta))
                    continue;
                rotated = 1;
                zeta = (beta - alpha) / (2.0 * gamma);
                tt = (zeta >= 0.0 ? 1.0 : -1.0)
                     / (fabs(zeta) + sqrt(1.0 + zeta * zeta));
                c = 1.0 / sqrt(1.0 + tt * tt);
                s = c * tt;
                for (k = 0; k < m; k++) {
                    double ap = A[k * n + p], aq = A[k * n + q];
                    A[k * n + p] = c * ap - s * aq;
                    A[k * n + q] = s * ap + c * aq;
                }
                for (k = 0; k < n; k++) {
                    double vp = V[k * n + p], vq = V[k * n + q];
                    V[k * n + p] = c * vp - s * vq;
                    V[k * n + q] = s * vp + c * vq;
                }
            }
        }
        if (!rotated)
            break;
    }

    for (j = 0; j < n; j++) {
        double nm = 0.0;
        for (k = 0; k < m; k++)
            nm += A[k * n + j] * A[k * n + j];
        nm = sqrt(nm);
        S[j] = nm;
        if (nm > 0.0)
            for (k = 0; k < m; k++)
                A[k * n + j] /= nm;
    }

    /* order singular values descending, permuting U and V alongside */
    for (i = 0; i < n - 1; i++) {
        int mx = i;
        double tmp;
        for (j = i + 1; j < n; j++)
            if (S[j] > S[mx])
                mx = j;
        if (mx == i)
            continue;
        tmp = S[i]; S[i] = S[mx]; S[mx] = tmp;
        for (k = 0; k < m; k++) {
            tmp = A[k * n + i];
            A[k * n + i] = A[k * n + mx];
            A[k * n + mx] = tmp;
        }
        for (k = 0; k < n; k++) {
            tmp = V[k * n + i];
            V[k * n + i] = V[k * n + mx];
            V[k * n + mx] = tmp;
        }
    }
    return MV_OK;
}

int mv_nullvec(double *x, const double *A, int m, int n)
{
    int rows = (m > n) ? m : n; /* zero-pad so the Jacobi SVD sees rows >= n */
    double *W = (double *)calloc((size_t)rows * n, sizeof(double));
    double *S = (double *)malloc((size_t)n * sizeof(double));
    double *V = (double *)malloc((size_t)n * n * sizeof(double));
    int k, ret = MV_ERR;

    if (W && S && V) {
        memcpy(W, A, (size_t)m * n * sizeof(double));
        if (mv_svd(W, S, V, rows, n) == MV_OK) {
            for (k = 0; k < n; k++)
                x[k] = V[k * n + (n - 1)];
            ret = MV_OK;
        }
    }
    free(W);
    free(S);
    free(V);
    return ret;
}
