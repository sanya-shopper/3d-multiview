#ifndef MV_MAT_H
#define MV_MAT_H

/* Small dense linear algebra, sized for multiview geometry problems
 * (design matrices are at most a few hundred rows by <= 9 columns).
 * No external dependencies; all routines are deterministic. */

/* C(mxp) = A(mxn) * B(nxp). C must not alias A or B. */
void mv_mat_mul(double *C, const double *A, const double *B,
                int m, int n, int p);

/* T(nxm) = A(mxn)^T. T must not alias A. */
void mv_mat_transpose(double *T, const double *A, int m, int n);

/* inv = A^-1 for 3x3 A; returns MV_ERR if singular. */
int mv_mat_inv3(double inv[9], const double A[9]);

/* M = [v]_x, the 3x3 skew-symmetric cross-product matrix. */
void mv_skew3(double M[9], const double v[3]);

void mv_cross3(double c[3], const double a[3], const double b[3]);
double mv_dot(const double *a, const double *b, int n);
double mv_norm(const double *a, int n);
/* Scale a to unit norm; MV_ERR if the norm is ~0. */
int mv_normalize(double *a, int n);

/* Thin SVD by one-sided (Hestenes) Jacobi, for m >= n.
 * On entry A holds the m x n matrix; on exit A holds U (m x n, orthonormal
 * columns), S the n singular values in descending order, V the n x n right
 * singular vectors (columns), so that A_in = U * diag(S) * V^T. */
int mv_svd(double *A, double *S, double *V, int m, int n);

/* x (length n) = right singular vector of A (m x n) with smallest singular
 * value, i.e. the least-squares solution of A x = 0 with |x| = 1.
 * Works for m < n as well (rows are zero-padded internally). A is untouched. */
int mv_nullvec(double *x, const double *A, int m, int n);

#endif /* MV_MAT_H */
