#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "mv/core.h"
#include "mv/mat.h"
#include "mv/rot.h"
#include "mv/bundle.h"

/* Joint Levenberg-Marquardt bundle adjustment [Triggs2000].
 *
 * Parameter vector layout (npar = NC*ncam + 3*npts):
 *   p[NC*i ...]           camera i {rodrigues r[3], t[3], fx, fy, cx, cy}
 *   p[NC*ncam + 3*j ...]  point j {x, y, z}
 * Slots masked by MV_BUNDLE_FIX_* are carried in p but never free.
 *
 * The residual vector stacks, observation by observation, the x and y
 * errors of mv_cam_project against the raw pixels.  An observation
 * touches one camera block and one point block only, so the normal
 * matrix has the classic BA structure: camera-diagonal U, per-point
 * 3x3 blocks V_j, couplings W.  Points vastly outnumber cameras, so
 * each solve takes the Schur complement over the damped V_j, leaving a
 * dense reduced camera system of at most NC*ncam (ncam is small). */

enum { NC = 10 }; /* camera slots: r[3], t[3], fx, fy, cx, cy */

#define BUNDLE_MAX_ITER   100    /* outer (Jacobian) iteration cap */
#define BUNDLE_REL_TOL    1e-12  /* relative cost-change convergence */
#define BUNDLE_FD_REL     1e-6   /* fwd-diff relative step ~ sqrt(eps) */
#define BUNDLE_LAMBDA0    1e-3   /* customary Marquardt start */
#define BUNDLE_LAMBDA_UP  10.0   /* on reject */
#define BUNDLE_LAMBDA_DN  0.1    /* on accept */
#define BUNDLE_LAMBDA_MIN 1e-12  /* floor: keep the damping meaningful */
#define BUNDLE_LAMBDA_MAX 1e12   /* past this the step is ~0: either a
                                  * cost plateau (converged) or a system
                                  * damping cannot fix (structural) */

/* ---- Rodrigues <-> rotation matrix (as src/refine.c) -------------- */

/* R = exp([r]x), row-major.  Below THETA_SERIES the closed form loses
 * digits to cancellation; the truncated series is exact to ~1e-16
 * there (next term is theta^4/120 < 1e-17). */
#define THETA_SERIES 1e-4


/* r = log(R) for orthonormal R.  Three regimes: generic (axis from the
 * antisymmetric part), theta ~ 0 (r = w exact to O(theta^3)), theta ~ pi
 * (antisymmetric part vanishes; axis from the diagonal of R, sign from
 * the residual antisymmetric part; at exactly pi the sign is a genuine
 * ambiguity and either choice is the same rotation). */

/* ---- small Cholesky (as src/refine.c, size run-time n) ------------ */

/* In-place lower Cholesky A = L L^T.  The damped normal blocks are
 * positive definite unless some parameter has no effect on any
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

static int inv3_chol(double inv[9], const double A[9])
{
    double L[9], bcol[3], xcol[3];
    int i, j;
    memcpy(L, A, sizeof(L));
    if (chol_factor(L, 3) != MV_OK)
        return MV_ERR;
    for (j = 0; j < 3; j++) {
        for (i = 0; i < 3; i++)
            bcol[i] = (i == j) ? 1.0 : 0.0;
        chol_solve(L, xcol, bcol, 3);
        for (i = 0; i < 3; i++)
            inv[i * 3 + j] = xcol[i];
    }
    return MV_OK;
}

/* ---- residual evaluation ----------------------------------------- */

/* Working camera i from base (fixed blocks, k1/k2, any fixed skew) and
 * its parameter slots pc (free blocks). */
static void cam_apply(mv_camera *w, const mv_camera *base, const double *pc,
                      unsigned flags)
{
    *w = *base;
    if (!(flags & MV_BUNDLE_FIX_POSE)) {
        mv_rot_exp(w->R, pc);
        w->t[0] = pc[3];
        w->t[1] = pc[4];
        w->t[2] = pc[5];
    }
    if (!(flags & MV_BUNDLE_FIX_K)) {
        memset(w->K, 0, sizeof(w->K));
        w->K[0] = pc[6];
        w->K[4] = pc[7];
        w->K[2] = pc[8];
        w->K[5] = pc[9];
        w->K[8] = 1.0; /* skew fixed at 0 for free-K cameras */
    }
}

static int obs_residual(double r[2], const mv_camera *w, const double *Xp,
                        const mv_bundle_obs *o)
{
    double uv[2];
    if (mv_cam_project(uv, w, Xp) != MV_OK)
        return MV_ERR;
    r[0] = uv[0] - o->u;
    r[1] = uv[1] - o->v;
    return MV_OK;
}

/* full residual vector under parameters p; wk receives the ncam working
 * cameras; cost (sum of squares) optional */
static int eval_all(double *cost, double *e, mv_camera *wk,
                    const mv_camera *base, const unsigned *camflags,
                    int ncam, const double *p, const mv_bundle_obs *obs,
                    int nobs)
{
    const double *po = p + NC * ncam;
    int i, k;
    for (i = 0; i < ncam; i++)
        cam_apply(&wk[i], &base[i], p + NC * i, camflags[i]);
    for (k = 0; k < nobs; k++)
        if (obs_residual(e + 2 * k, &wk[obs[k].cam],
                         po + 3 * obs[k].pt, &obs[k]) != MV_OK)
            return MV_ERR;
    if (cost) {
        double c = 0.0;
        for (k = 0; k < 2 * nobs; k++)
            c += e[k] * e[k];
        *cost = c;
    }
    return MV_OK;
}

/* ---- driver ------------------------------------------------------ */

int mv_bundle_adjust(mv_camera *cams, const unsigned *camflags, int ncam,
                     double *X, int npts, const mv_bundle_obs *obs,
                     int nobs, double *rms_out)
{
    double *p = NULL, *ptry = NULL, *e = NULL, *etry = NULL;
    double *Jc = NULL, *Jp = NULL, *Wobs = NULL;
    double *U = NULL, *S = NULL, *gcam = NULL, *rhs = NULL, *dc = NULL;
    double *Vp = NULL, *Vinv = NULL, *gp = NULL, *Wc = NULL;
    mv_camera *work = NULL, *wtry = NULL;
    int *cmap = NULL, *coff = NULL, *cidx = NULL;
    int *poff = NULL, *pidx = NULL, *mark = NULL, *clist = NULL;
    double cost, lambda = BUNDLE_LAMBDA0;
    int npar, ncf = 0, i, j, k, s, it, done = 0, ret = MV_ERR;

    if (!cams || !X || !obs || ncam < 1 || npts < 1 || nobs < 1)
        return MV_ERR;
    if (!camflags)
        return MV_ERR; /* no flags: no gauge anchor possible */
    for (i = 0, j = 0; i < ncam; i++)
        if (camflags[i] & MV_BUNDLE_FIX_POSE)
            j = 1;
    if (!j)
        return MV_ERR; /* gauge would be loose without a fixed pose */
    for (k = 0; k < nobs; k++)
        if (obs[k].cam < 0 || obs[k].cam >= ncam || obs[k].pt < 0
            || obs[k].pt >= npts)
            return MV_ERR;

    npar = NC * ncam + 3 * npts;
    cmap = (int *)malloc((size_t)NC * ncam * sizeof(int));
    coff = (int *)malloc((size_t)(ncam + 1) * sizeof(int));
    cidx = (int *)malloc((size_t)nobs * sizeof(int));
    poff = (int *)malloc((size_t)(npts + 1) * sizeof(int));
    pidx = (int *)malloc((size_t)nobs * sizeof(int));
    mark = (int *)malloc((size_t)ncam * sizeof(int));
    clist = (int *)malloc((size_t)ncam * sizeof(int));
    if (!cmap || !coff || !cidx || !poff || !pidx || !mark || !clist)
        goto out;

    /* free-slot map: global reduced-system index per camera slot */
    for (i = 0; i < ncam; i++)
        for (s = 0; s < NC; s++) {
            int fixed = (s < 6) ? (camflags[i] & MV_BUNDLE_FIX_POSE)
                                : (camflags[i] & MV_BUNDLE_FIX_K);
            cmap[NC * i + s] = fixed ? -1 : ncf++;
        }

    /* observation buckets by camera and by point (counting sort) */
    memset(coff, 0, (size_t)(ncam + 1) * sizeof(int));
    memset(poff, 0, (size_t)(npts + 1) * sizeof(int));
    for (k = 0; k < nobs; k++) {
        coff[obs[k].cam + 1]++;
        poff[obs[k].pt + 1]++;
    }
    for (i = 0; i < ncam; i++)
        coff[i + 1] += coff[i];
    for (j = 0; j < npts; j++)
        poff[j + 1] += poff[j];
    {
        int *cfill = (int *)malloc((size_t)ncam * sizeof(int));
        int *pfill = (int *)malloc((size_t)npts * sizeof(int));
        if (!cfill || !pfill) {
            free(cfill);
            free(pfill);
            goto out;
        }
        memcpy(cfill, coff, (size_t)ncam * sizeof(int));
        memcpy(pfill, poff, (size_t)npts * sizeof(int));
        for (k = 0; k < nobs; k++) {
            cidx[cfill[obs[k].cam]++] = k;
            pidx[pfill[obs[k].pt]++] = k;
        }
        free(cfill);
        free(pfill);
    }

    p = (double *)malloc((size_t)npar * sizeof(double));
    ptry = (double *)malloc((size_t)npar * sizeof(double));
    e = (double *)malloc((size_t)(2 * nobs) * sizeof(double));
    etry = (double *)malloc((size_t)(2 * nobs) * sizeof(double));
    Jc = (double *)calloc((size_t)(2 * NC * nobs), sizeof(double));
    Jp = (double *)malloc((size_t)(6 * nobs) * sizeof(double));
    Wobs = (double *)malloc((size_t)(3 * NC * nobs) * sizeof(double));
    U = (double *)malloc((size_t)(ncf * ncf + 1) * sizeof(double));
    S = (double *)malloc((size_t)(ncf * ncf + 1) * sizeof(double));
    gcam = (double *)malloc((size_t)(ncf + 1) * sizeof(double));
    rhs = (double *)malloc((size_t)(ncf + 1) * sizeof(double));
    dc = (double *)malloc((size_t)(ncf + 1) * sizeof(double));
    Vp = (double *)malloc((size_t)(9 * npts) * sizeof(double));
    Vinv = (double *)malloc((size_t)(9 * npts) * sizeof(double));
    gp = (double *)malloc((size_t)(3 * npts) * sizeof(double));
    Wc = (double *)malloc((size_t)(3 * NC * ncam) * sizeof(double));
    work = (mv_camera *)malloc((size_t)ncam * sizeof(mv_camera));
    wtry = (mv_camera *)malloc((size_t)ncam * sizeof(mv_camera));
    if (!p || !ptry || !e || !etry || !Jc || !Jp || !Wobs || !U || !S
        || !gcam || !rhs || !dc || !Vp || !Vinv || !gp || !Wc || !work
        || !wtry)
        goto out;

    /* initial parameter vector from the cameras and points as given */
    for (i = 0; i < ncam; i++) {
        double *pc = p + NC * i;
        mv_rot_log(pc, cams[i].R);
        pc[3] = cams[i].t[0];
        pc[4] = cams[i].t[1];
        pc[5] = cams[i].t[2];
        pc[6] = cams[i].K[0];
        pc[7] = cams[i].K[4];
        pc[8] = cams[i].K[2];
        pc[9] = cams[i].K[5];
    }
    memcpy(p + NC * ncam, X, (size_t)(3 * npts) * sizeof(double));
    for (i = 0; i < ncam; i++)
        mark[i] = -1;

    if (eval_all(&cost, e, work, cams, camflags, ncam, p, obs, nobs)
        != MV_OK)
        goto out; /* initial values do not project: bad inputs */

    for (it = 0; it < BUNDLE_MAX_ITER && !done; it++) {
        double *po = p + NC * ncam;

        /* forward-difference Jacobian about p, bucket by bucket:
         * a camera slot touches only that camera's observations, a
         * point coordinate only that point's (baseline e reused) */
        for (i = 0; i < ncam; i++)
            for (s = 0; s < NC; s++) {
                int pc = NC * i + s, kk;
                double save, h, pr[2];
                mv_camera tmp;
                if (cmap[pc] < 0)
                    continue;
                save = p[pc];
                p[pc] = save + BUNDLE_FD_REL * (1.0 + fabs(save));
                h = p[pc] - save; /* representable step actually taken */
                cam_apply(&tmp, &cams[i], p + NC * i, camflags[i]);
                for (kk = coff[i]; kk < coff[i + 1]; kk++) {
                    k = cidx[kk];
                    if (obs_residual(pr, &tmp, po + 3 * obs[k].pt,
                                     &obs[k]) != MV_OK) {
                        p[pc] = save;
                        goto out; /* domain edge at infinitesimal step */
                    }
                    Jc[2 * NC * k + s] = (pr[0] - e[2 * k]) / h;
                    Jc[2 * NC * k + NC + s] = (pr[1] - e[2 * k + 1]) / h;
                }
                p[pc] = save;
            }
        for (j = 0; j < npts; j++)
            for (s = 0; s < 3; s++) {
                int pc = NC * ncam + 3 * j + s, kk;
                double save = p[pc], h, pr[2];
                p[pc] = save + BUNDLE_FD_REL * (1.0 + fabs(save));
                h = p[pc] - save;
                for (kk = poff[j]; kk < poff[j + 1]; kk++) {
                    k = pidx[kk];
                    if (obs_residual(pr, &work[obs[k].cam], po + 3 * j,
                                     &obs[k]) != MV_OK) {
                        p[pc] = save;
                        goto out;
                    }
                    Jp[6 * k + s] = (pr[0] - e[2 * k]) / h;
                    Jp[6 * k + 3 + s] = (pr[1] - e[2 * k + 1]) / h;
                }
                p[pc] = save;
            }

        /* accumulate normal blocks, gradient, and per-obs couplings */
        memset(U, 0, (size_t)(ncf * ncf) * sizeof(double));
        memset(gcam, 0, (size_t)ncf * sizeof(double));
        memset(Vp, 0, (size_t)(9 * npts) * sizeof(double));
        memset(gp, 0, (size_t)(3 * npts) * sizeof(double));
        for (k = 0; k < nobs; k++) {
            const double *J0 = Jc + 2 * NC * k, *J1 = J0 + NC;
            const double *P0 = Jp + 6 * k, *P1 = P0 + 3;
            const int *cm = cmap + NC * obs[k].cam;
            double e0 = e[2 * k], e1 = e[2 * k + 1];
            int a, b;
            for (a = 0; a < NC; a++) {
                int g1 = cm[a];
                if (g1 < 0)
                    continue;
                gcam[g1] -= J0[a] * e0 + J1[a] * e1;
                for (b = 0; b < NC; b++) {
                    int g2 = cm[b];
                    if (g2 >= 0)
                        U[g1 * ncf + g2] += J0[a] * J0[b] + J1[a] * J1[b];
                }
            }
            j = obs[k].pt;
            for (a = 0; a < 3; a++) {
                gp[3 * j + a] -= P0[a] * e0 + P1[a] * e1;
                for (b = 0; b < 3; b++)
                    Vp[9 * j + 3 * a + b] += P0[a] * P0[b] + P1[a] * P1[b];
            }
            for (a = 0; a < NC; a++) /* W row a is 0 for fixed slots */
                for (b = 0; b < 3; b++)
                    Wobs[3 * NC * k + 3 * a + b] =
                        J0[a] * P0[b] + J1[a] * P1[b];
        }

        /* damping search on this Jacobian */
        for (;;) {
            int solved = MV_OK, kk, a, b, c;

            /* reduced camera system S = U* - sum_j W_j V_j*^-1 W_j^T,
             * built point by point (only the cameras observing a point
             * couple through it) */
            memcpy(S, U, (size_t)(ncf * ncf) * sizeof(double));
            for (a = 0; a < ncf; a++)
                S[a * ncf + a] *= 1.0 + lambda;
            memcpy(rhs, gcam, (size_t)ncf * sizeof(double));
            for (j = 0; j < npts && solved == MV_OK; j++) {
                double Vd[9], Y[3 * NC];
                int nl = 0, l1, l2;
                if (poff[j] == poff[j + 1])
                    continue; /* unobserved point: held fixed */
                memcpy(Vd, Vp + 9 * j, sizeof(Vd));
                for (a = 0; a < 3; a++)
                    Vd[3 * a + a] *= 1.0 + lambda;
                if (inv3_chol(Vinv + 9 * j, Vd) != MV_OK) {
                    solved = MV_ERR;
                    break;
                }
                for (kk = poff[j]; kk < poff[j + 1]; kk++) {
                    k = pidx[kk];
                    c = obs[k].cam;
                    if (mark[c] < 0) {
                        mark[c] = nl;
                        clist[nl] = c;
                        memset(Wc + 3 * NC * nl, 0,
                               (size_t)(3 * NC) * sizeof(double));
                        nl++;
                    }
                    for (a = 0; a < 3 * NC; a++)
                        Wc[3 * NC * mark[c] + a] += Wobs[3 * NC * k + a];
                }
                for (l1 = 0; l1 < nl; l1++) {
                    const double *W1 = Wc + 3 * NC * l1;
                    const int *cm1 = cmap + NC * clist[l1];
                    for (a = 0; a < NC; a++)
                        for (b = 0; b < 3; b++)
                            Y[3 * a + b] =
                                W1[3 * a + 0] * Vinv[9 * j + 0 + b]
                                + W1[3 * a + 1] * Vinv[9 * j + 3 + b]
                                + W1[3 * a + 2] * Vinv[9 * j + 6 + b];
                    for (a = 0; a < NC; a++) {
                        int g1 = cm1[a];
                        if (g1 < 0)
                            continue;
                        rhs[g1] -= Y[3 * a + 0] * gp[3 * j]
                                 + Y[3 * a + 1] * gp[3 * j + 1]
                                 + Y[3 * a + 2] * gp[3 * j + 2];
                        for (l2 = 0; l2 < nl; l2++) {
                            const double *W2 = Wc + 3 * NC * l2;
                            const int *cm2 = cmap + NC * clist[l2];
                            for (b = 0; b < NC; b++) {
                                int g2 = cm2[b];
                                if (g2 >= 0)
                                    S[g1 * ncf + g2] -=
                                        Y[3 * a + 0] * W2[3 * b + 0]
                                        + Y[3 * a + 1] * W2[3 * b + 1]
                                        + Y[3 * a + 2] * W2[3 * b + 2];
                            }
                        }
                    }
                }
                for (l1 = 0; l1 < nl; l1++)
                    mark[clist[l1]] = -1;
            }
            if (solved == MV_OK && ncf > 0
                && chol_factor(S, ncf) != MV_OK)
                solved = MV_ERR;
            if (solved == MV_OK) {
                double ctry;
                if (ncf > 0)
                    chol_solve(S, dc, rhs, ncf);
                memcpy(ptry, p, (size_t)npar * sizeof(double));
                for (a = 0; a < NC * ncam; a++)
                    if (cmap[a] >= 0)
                        ptry[a] += dc[cmap[a]];
                /* back-substitute: dp_j = V_j*^-1 (g_j - W_j^T dcam) */
                for (j = 0; j < npts; j++) {
                    double t3[3];
                    if (poff[j] == poff[j + 1])
                        continue;
                    t3[0] = gp[3 * j];
                    t3[1] = gp[3 * j + 1];
                    t3[2] = gp[3 * j + 2];
                    for (kk = poff[j]; kk < poff[j + 1]; kk++) {
                        const double *Wk;
                        const int *cm;
                        k = pidx[kk];
                        Wk = Wobs + 3 * NC * k;
                        cm = cmap + NC * obs[k].cam;
                        for (a = 0; a < NC; a++) {
                            if (cm[a] < 0)
                                continue;
                            for (b = 0; b < 3; b++)
                                t3[b] -= Wk[3 * a + b] * dc[cm[a]];
                        }
                    }
                    for (a = 0; a < 3; a++)
                        ptry[NC * ncam + 3 * j + a] +=
                            Vinv[9 * j + 3 * a] * t3[0]
                            + Vinv[9 * j + 3 * a + 1] * t3[1]
                            + Vinv[9 * j + 3 * a + 2] * t3[2];
                }
                if (eval_all(&ctry, etry, wtry, cams, camflags, ncam,
                             ptry, obs, nobs) == MV_OK
                    && ctry < cost) {
                    done = (cost - ctry) <= BUNDLE_REL_TOL * cost;
                    memcpy(p, ptry, (size_t)npar * sizeof(double));
                    memcpy(e, etry, (size_t)(2 * nobs) * sizeof(double));
                    memcpy(work, wtry, (size_t)ncam * sizeof(mv_camera));
                    cost = ctry;
                    lambda *= BUNDLE_LAMBDA_DN;
                    if (lambda < BUNDLE_LAMBDA_MIN)
                        lambda = BUNDLE_LAMBDA_MIN;
                    break;
                }
            }
            lambda *= BUNDLE_LAMBDA_UP;
            if (lambda > BUNDLE_LAMBDA_MAX) {
                if (solved != MV_OK)
                    goto out; /* singular at max damping: structural */
                done = 1; /* solvable but no descent: plateau */
                break;
            }
        }
    }

    /* write back the refined free blocks (work already holds them) */
    for (i = 0; i < ncam; i++)
        cams[i] = work[i];
    memcpy(X, p + NC * ncam, (size_t)(3 * npts) * sizeof(double));
    if (rms_out)
        *rms_out = sqrt(cost / (double)nobs);
    ret = MV_OK;
out:
    free(cmap);
    free(coff);
    free(cidx);
    free(poff);
    free(pidx);
    free(mark);
    free(clist);
    free(p);
    free(ptry);
    free(e);
    free(etry);
    free(Jc);
    free(Jp);
    free(Wobs);
    free(U);
    free(S);
    free(gcam);
    free(rhs);
    free(dc);
    free(Vp);
    free(Vinv);
    free(gp);
    free(Wc);
    free(work);
    free(wtry);
    return ret;
}
