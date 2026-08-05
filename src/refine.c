#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "mv/core.h"
#include "mv/mat.h"
#include "mv/refine.h"

/* Joint Levenberg-Marquardt refinement of plane-based calibration.
 *
 * Parameter vector layout (np = NS + NB*nviews):
 *   p[0..5]           shared {fx, fy, cx, cy, k1, k2} (skew fixed at 0)
 *   p[NS + NB*i ...]  view i {rodrigues r[3], t[3]}
 *
 * The residual vector stacks, view by view, the x and y reprojection
 * errors of the DISTORTED projection (mv_cam_project) against the raw
 * observations, so the minimum is the maximum-likelihood estimate under
 * i.i.d. isotropic pixel noise.  The normal matrix is block-arrow:
 * a shared 6x6 block U, per-view 6x6 blocks V_i, couplings W_i; the
 * Schur complement over the V_i reduces each solve to nviews damped
 * 6x6 inversions plus one 6x6 system, linear in nviews. */

enum { NS = 6, NB = 6 }; /* shared / per-view parameter block sizes */

#define REFINE_MAX_ITER   100    /* outer (Jacobian) iteration cap */
#define REFINE_REL_TOL    1e-12  /* relative cost-change convergence */
#define REFINE_FD_REL     1e-6   /* fwd-diff relative step ~ sqrt(eps):
                                  * balances truncation vs roundoff */
#define REFINE_LAMBDA0    1e-3   /* customary Marquardt start */
#define REFINE_LAMBDA_UP  10.0   /* on reject */
#define REFINE_LAMBDA_DN  0.1    /* on accept */
#define REFINE_LAMBDA_MIN 1e-12  /* floor: keep the damping meaningful */
#define REFINE_LAMBDA_MAX 1e12   /* past this the step is ~0: either a
                                  * cost plateau (converged) or a system
                                  * damping cannot fix (structural) */

/* ---- Rodrigues <-> rotation matrix ------------------------------- */

/* R = exp([r]x), row-major.  Below THETA_SERIES the closed form loses
 * digits to cancellation; the truncated series is exact to ~1e-16
 * there (next term is theta^4/120 < 1e-17). */
#define THETA_SERIES 1e-4

static void rod_to_mat(double R[9], const double r[3])
{
    double th2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
    double th = sqrt(th2);
    double a, b; /* a = sin(th)/th, b = (1-cos(th))/th^2 */
    if (th < THETA_SERIES) {
        a = 1.0 - th2 / 6.0;
        b = 0.5 - th2 / 24.0;
    } else {
        a = sin(th) / th;
        b = (1.0 - cos(th)) / th2;
    }
    R[0] = 1.0 - b * (r[1] * r[1] + r[2] * r[2]);
    R[1] = b * r[0] * r[1] - a * r[2];
    R[2] = b * r[0] * r[2] + a * r[1];
    R[3] = b * r[0] * r[1] + a * r[2];
    R[4] = 1.0 - b * (r[0] * r[0] + r[2] * r[2]);
    R[5] = b * r[1] * r[2] - a * r[0];
    R[6] = b * r[0] * r[2] - a * r[1];
    R[7] = b * r[1] * r[2] + a * r[0];
    R[8] = 1.0 - b * (r[0] * r[0] + r[1] * r[1]);
}

/* r = log(R) for orthonormal R.  Three regimes: generic (axis from the
 * antisymmetric part), theta ~ 0 (r = w exact to O(theta^3)), theta ~ pi
 * (antisymmetric part vanishes; axis from the diagonal of R, sign from
 * the residual antisymmetric part; at exactly pi the sign is a genuine
 * ambiguity and either choice is the same rotation). */
static void mat_to_rod(double r[3], const double R[9])
{
    double w[3]; /* = sin(theta) * axis */
    double c = 0.5 * (R[0] + R[4] + R[8] - 1.0);
    double s, th;
    w[0] = 0.5 * (R[7] - R[5]);
    w[1] = 0.5 * (R[2] - R[6]);
    w[2] = 0.5 * (R[3] - R[1]);
    if (c > 1.0)
        c = 1.0;
    if (c < -1.0)
        c = -1.0;
    s = mv_norm(w, 3);
    th = atan2(s, c);
    if (s > 1e-6) { /* generic: sin(theta) well away from 0 */
        double f = th / s;
        r[0] = f * w[0];
        r[1] = f * w[1];
        r[2] = f * w[2];
    } else if (c > 0.0) { /* theta ~ 0 */
        r[0] = w[0];
        r[1] = w[1];
        r[2] = w[2];
    } else { /* theta ~ pi: R ~ 2 a a^T - I, so a_i^2 = (R_ii + 1)/2 */
        double a[3], na;
        int imax = 0, j;
        if (R[4] > R[0])
            imax = 1;
        if (R[8] > R[imax * 3 + imax])
            imax = 2;
        a[imax] = sqrt(0.5 * (R[imax * 3 + imax] + 1.0));
        for (j = 0; j < 3; j++)
            if (j != imax) /* symmetrized off-diagonal: R_ij = 2 a_i a_j */
                a[j] = 0.25 * (R[imax * 3 + j] + R[j * 3 + imax]) / a[imax];
        na = mv_norm(a, 3);
        if (na > 0.0) {
            for (j = 0; j < 3; j++)
                a[j] /= na;
        } else { /* unreachable for orthonormal input; keep r finite */
            a[0] = 1.0;
            a[1] = a[2] = 0.0;
        }
        if (a[0] * w[0] + a[1] * w[1] + a[2] * w[2] < 0.0)
            for (j = 0; j < 3; j++)
                a[j] = -a[j];
        for (j = 0; j < 3; j++)
            r[j] = th * a[j];
    }
}

/* ---- small Cholesky ---------------------------------------------- */

/* In-place lower Cholesky A = L L^T, n <= 6.  The damped normal blocks
 * are positive definite unless some parameter has no effect on any
 * residual, which multiplicative damping cannot repair. */
static int chol_factor(double *A, int n)
{
    int i, j, k;
    for (j = 0; j < n; j++) {
        double d = A[j * n + j];
        for (k = 0; k < j; k++)
            d -= A[j * n + k] * A[j * n + k];
        if (d <= 0.0)
            return MV_ERR;
        d = sqrt(d);
        A[j * n + j] = d;
        for (i = j + 1; i < n; i++) {
            double v = A[i * n + j];
            for (k = 0; k < j; k++)
                v -= A[i * n + k] * A[j * n + k];
            A[i * n + j] = v / d;
        }
    }
    return MV_OK;
}

static void chol_solve(const double *L, double *x, const double *b, int n)
{
    int i, k;
    for (i = 0; i < n; i++) { /* L y = b */
        double v = b[i];
        for (k = 0; k < i; k++)
            v -= L[i * n + k] * x[k];
        x[i] = v / L[i * n + i];
    }
    for (i = n - 1; i >= 0; i--) { /* L^T x = y */
        double v = x[i];
        for (k = i + 1; k < n; k++)
            v -= L[k * n + i] * x[k];
        x[i] = v / L[i * n + i];
    }
}

static int inv6_chol(double inv[36], const double A[36])
{
    double L[36], bcol[6], xcol[6];
    int i, j;
    memcpy(L, A, sizeof(L));
    if (chol_factor(L, 6) != MV_OK)
        return MV_ERR;
    for (j = 0; j < 6; j++) {
        for (i = 0; i < 6; i++)
            bcol[i] = (i == j) ? 1.0 : 0.0;
        chol_solve(L, xcol, bcol, 6);
        for (i = 0; i < 6; i++)
            inv[i * 6 + j] = xcol[i];
    }
    return MV_OK;
}

/* ---- residual evaluation ----------------------------------------- */

static void cam_from_params(mv_camera *c, const double *s, const double *vb,
                            const double kfix[3])
{
    memset(c->K, 0, sizeof(c->K));
    c->K[0] = s[0];
    c->K[4] = s[1];
    c->K[2] = s[2];
    c->K[5] = s[3];
    c->K[8] = 1.0; /* skew K[1] fixed at 0 */
    rod_to_mat(c->R, vb);
    c->t[0] = vb[3];
    c->t[1] = vb[4];
    c->t[2] = vb[5];
    c->k[0] = s[4];
    c->k[1] = s[5];
    c->k[2] = kfix[0]; /* p1, p2, k3 held at their input values */
    c->k[3] = kfix[1];
    c->k[4] = kfix[2];
}

/* res[2n] = distorted projection minus raw observation, view vw */
static int view_residuals(double *res, const double *s, const double *vb,
                          const double kfix[3], const mv_calib_view *vw)
{
    mv_camera c;
    int j;
    cam_from_params(&c, s, vb, kfix);
    for (j = 0; j < vw->n; j++) {
        double X[3], uv[2];
        X[0] = vw->obj[2 * j];
        X[1] = vw->obj[2 * j + 1];
        X[2] = 0.0;
        if (mv_cam_project(uv, &c, X) != MV_OK)
            return MV_ERR;
        res[2 * j] = uv[0] - vw->img[2 * j];
        res[2 * j + 1] = uv[1] - vw->img[2 * j + 1];
    }
    return MV_OK;
}

/* full residual vector; cost (sum of squares) optional */
static int eval_all(double *cost, double *res, const double *p,
                    const double kfix[3], const mv_calib_view *views,
                    int nviews, const int *off)
{
    int i, r;
    for (i = 0; i < nviews; i++)
        if (view_residuals(res + off[i], p, p + NS + NB * i, kfix,
                           &views[i]) != MV_OK)
            return MV_ERR;
    if (cost) {
        double c = 0.0;
        for (i = 0; i < nviews; i++)
            for (r = 0; r < 2 * views[i].n; r++)
                c += res[off[i] + r] * res[off[i] + r];
        *cost = c;
    }
    return MV_OK;
}

/* ---- damped Schur solve ------------------------------------------ */

/* Solve (J^T J + lambda diag) delta = -J^T e using the block-arrow
 * structure.  U, V, W, ga, gb are the undamped accumulations; Vinv
 * (36*nviews) receives the damped V_i inverses for the back-substitution.
 * MV_ERR: a Cholesky failed at this damping level. */
static int schur_step(double *delta, const double *U, const double *V,
                      const double *W, const double *ga, const double *gb,
                      double *Vinv, int nviews, double lambda)
{
    double S[36], L[36], rhs[6], da[6];
    int i, j;

    memcpy(S, U, sizeof(S));
    for (j = 0; j < 6; j++)
        S[j * 6 + j] *= 1.0 + lambda;
    memcpy(rhs, ga, sizeof(rhs));
    for (i = 0; i < nviews; i++) {
        const double *Wi = W + 36 * i;
        double Vd[36], Wt[36], M[36], MWt[36], t6[6];
        memcpy(Vd, V + 36 * i, sizeof(Vd));
        for (j = 0; j < 6; j++)
            Vd[j * 6 + j] *= 1.0 + lambda;
        if (inv6_chol(Vinv + 36 * i, Vd) != MV_OK)
            return MV_ERR;
        mv_mat_mul(M, Wi, Vinv + 36 * i, 6, 6, 6);
        mv_mat_transpose(Wt, Wi, 6, 6);
        mv_mat_mul(MWt, M, Wt, 6, 6, 6);
        for (j = 0; j < 36; j++)
            S[j] -= MWt[j]; /* S = U* - sum W V*^-1 W^T */
        mv_mat_mul(t6, M, gb + 6 * i, 6, 6, 1);
        for (j = 0; j < 6; j++)
            rhs[j] -= t6[j];
    }
    memcpy(L, S, sizeof(L));
    if (chol_factor(L, 6) != MV_OK)
        return MV_ERR;
    chol_solve(L, da, rhs, 6);
    memcpy(delta, da, sizeof(da));
    for (i = 0; i < nviews; i++) { /* delta_b = V*^-1 (gb - W^T da) */
        double Wt[36], t6[6], db[6];
        mv_mat_transpose(Wt, W + 36 * i, 6, 6);
        mv_mat_mul(t6, Wt, da, 6, 6, 1);
        for (j = 0; j < 6; j++)
            t6[j] = gb[6 * i + j] - t6[j];
        mv_mat_mul(db, Vinv + 36 * i, t6, 6, 6, 1);
        memcpy(delta + NS + NB * i, db, sizeof(db));
    }
    return MV_OK;
}

/* ---- driver ------------------------------------------------------ */

int mv_calib_refine(double K[9], mv_camera *cams, double k_radial[2],
                    const mv_calib_view *views, int nviews)
{
    double *buf = NULL, *p, *ptry, *delta, *e, *etry, *epert, *A, *B;
    double *V, *W, *Vinv, *gb;
    double U[36], ga[6], kfix[3];
    double cost, lambda = REFINE_LAMBDA0;
    int *off = NULL;
    int np, ntot = 0, m, i, j, c, it, done = 0, ret = MV_ERR;

    if (!K || !cams || !views || nviews < 1)
        return MV_ERR;
    for (i = 0; i < nviews; i++) {
        if (!views[i].obj || !views[i].img || views[i].n < 4)
            return MV_ERR;
        ntot += views[i].n;
    }
    np = NS + NB * nviews;
    m = 2 * ntot;

    /* p, ptry, delta | e, etry, epert | A (m x NS), B (m x NB) |
     * V, W, Vinv (36 each), gb (6) per view */
    buf = (double *)malloc((size_t)(3 * np + 3 * m + m * (NS + NB)
                                    + 114 * nviews) * sizeof(double));
    off = (int *)malloc((size_t)nviews * sizeof(int));
    if (!buf || !off)
        goto out;
    p = buf;
    ptry = p + np;
    delta = ptry + np;
    e = delta + np;
    etry = e + m;
    epert = etry + m;
    A = epert + m;
    B = A + m * NS;
    V = B + m * NB;
    W = V + 36 * nviews;
    Vinv = W + 36 * nviews;
    gb = Vinv + 36 * nviews;

    off[0] = 0;
    for (i = 1; i < nviews; i++)
        off[i] = off[i - 1] + 2 * views[i - 1].n;

    kfix[0] = cams[0].k[2];
    kfix[1] = cams[0].k[3];
    kfix[2] = cams[0].k[4];
    p[0] = K[0];
    p[1] = K[4];
    p[2] = K[2];
    p[3] = K[5];
    p[4] = k_radial ? k_radial[0] : cams[0].k[0];
    p[5] = k_radial ? k_radial[1] : cams[0].k[1];
    for (i = 0; i < nviews; i++) {
        mat_to_rod(p + NS + NB * i, cams[i].R);
        p[NS + NB * i + 3] = cams[i].t[0];
        p[NS + NB * i + 4] = cams[i].t[1];
        p[NS + NB * i + 5] = cams[i].t[2];
    }

    if (eval_all(&cost, e, p, kfix, views, nviews, off) != MV_OK)
        goto out; /* initial values do not project: bad inputs */

    for (it = 0; it < REFINE_MAX_ITER && !done; it++) {
        /* forward-difference Jacobian about p (baseline e reused) */
        for (c = 0; c < NS; c++) {
            double save = p[c], h;
            p[c] = save + REFINE_FD_REL * (1.0 + fabs(save));
            h = p[c] - save; /* representable step actually taken */
            if (eval_all(NULL, epert, p, kfix, views, nviews, off)
                != MV_OK) {
                p[c] = save;
                goto out; /* domain boundary at infinitesimal step */
            }
            p[c] = save;
            for (j = 0; j < m; j++)
                A[j * NS + c] = (epert[j] - e[j]) / h;
        }
        for (i = 0; i < nviews; i++) {
            for (c = 0; c < NB; c++) {
                int pc = NS + NB * i + c;
                double save = p[pc], h;
                p[pc] = save + REFINE_FD_REL * (1.0 + fabs(save));
                h = p[pc] - save;
                if (view_residuals(epert + off[i], p, p + NS + NB * i,
                                   kfix, &views[i]) != MV_OK) {
                    p[pc] = save;
                    goto out;
                }
                p[pc] = save;
                for (j = 0; j < 2 * views[i].n; j++)
                    B[(off[i] + j) * NB + c] =
                        (epert[off[i] + j] - e[off[i] + j]) / h;
            }
        }

        /* accumulate normal blocks and gradient */
        memset(U, 0, sizeof(U));
        memset(ga, 0, sizeof(ga));
        memset(V, 0, (size_t)(36 * nviews) * sizeof(double));
        memset(W, 0, (size_t)(36 * nviews) * sizeof(double));
        memset(gb, 0, (size_t)(6 * nviews) * sizeof(double));
        for (i = 0; i < nviews; i++) {
            int r, a, b;
            for (r = 0; r < 2 * views[i].n; r++) {
                const double *Ar = A + (off[i] + r) * NS;
                const double *Br = B + (off[i] + r) * NB;
                double er = e[off[i] + r];
                for (a = 0; a < NS; a++) {
                    for (b = 0; b < NS; b++)
                        U[a * NS + b] += Ar[a] * Ar[b];
                    for (b = 0; b < NB; b++)
                        W[36 * i + a * NB + b] += Ar[a] * Br[b];
                    ga[a] -= Ar[a] * er;
                }
                for (a = 0; a < NB; a++) {
                    for (b = 0; b < NB; b++)
                        V[36 * i + a * NB + b] += Br[a] * Br[b];
                    gb[6 * i + a] -= Br[a] * er;
                }
            }
        }

        /* damping search on this Jacobian */
        for (;;) {
            int solved = schur_step(delta, U, V, W, ga, gb, Vinv,
                                    nviews, lambda);
            if (solved == MV_OK) {
                double ctry;
                for (j = 0; j < np; j++)
                    ptry[j] = p[j] + delta[j];
                if (eval_all(&ctry, etry, ptry, kfix, views, nviews, off)
                        == MV_OK
                    && ctry < cost) {
                    done = (cost - ctry) <= REFINE_REL_TOL * cost;
                    memcpy(p, ptry, (size_t)np * sizeof(double));
                    memcpy(e, etry, (size_t)m * sizeof(double));
                    cost = ctry;
                    lambda *= REFINE_LAMBDA_DN;
                    if (lambda < REFINE_LAMBDA_MIN)
                        lambda = REFINE_LAMBDA_MIN;
                    break;
                }
            }
            lambda *= REFINE_LAMBDA_UP;
            if (lambda > REFINE_LAMBDA_MAX) {
                if (solved != MV_OK)
                    goto out; /* singular at max damping: structural */
                done = 1; /* solvable but no descent: plateau */
                break;
            }
        }
    }

    K[0] = p[0];
    K[1] = 0.0;
    K[2] = p[2];
    K[3] = 0.0;
    K[4] = p[1];
    K[5] = p[3];
    K[6] = 0.0;
    K[7] = 0.0;
    K[8] = 1.0;
    if (k_radial) {
        k_radial[0] = p[4];
        k_radial[1] = p[5];
    }
    for (i = 0; i < nviews; i++)
        cam_from_params(&cams[i], p, p + NS + NB * i, kfix);
    ret = MV_OK;
out:
    free(buf);
    free(off);
    return ret;
}
