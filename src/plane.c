#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "mv/core.h"
#include "mv/mat.h"
#include "mv/plane.h"

/* Paper: doc/multiview.tex, section "The representation ladder", rung 4.
 * OWNERSHIP (parallel build): this file, include/mv/plane.h, and
 * tests/test_plane.c belong to the TSDF-completion work item ONLY. */

/* house LCG (same constants as photo.c's RANSAC sampler) */
static unsigned long long pl_next(unsigned long long *s)
{
    *s = *s * 6364136223846793005ULL + 1442695040888963407ULL;
    return *s >> 33;
}

/* flip so d >= 0 (ties broken toward a positive leading normal
 * component), making (n,d) and (-n,-d) compare equal */
static void pl_canon(double plane[4])
{
    int a;
    int flip = plane[3] < 0.0;
    if (plane[3] == 0.0) {
        for (a = 0; a < 3; a++) {
            if (plane[a] != 0.0) {
                flip = plane[a] < 0.0;
                break;
            }
        }
    }
    if (flip)
        for (a = 0; a < 4; a++)
            plane[a] = -plane[a];
}

/* total least squares on the masked points: centroid + smallest
 * eigenvector of the 3x3 scatter matrix (mv_nullvec = smallest
 * singular vector, which for a symmetric PSD matrix is that
 * eigenvector).  mask == NULL means all points. */
static int pl_fit_ls(const double *pts, int n, const unsigned char *mask,
                     double plane[4])
{
    double c[3] = { 0.0, 0.0, 0.0 }, S[9], nvec[3];
    int i, a, b, m = 0;
    for (i = 0; i < n; i++) {
        if (mask && !mask[i])
            continue;
        for (a = 0; a < 3; a++)
            c[a] += pts[3 * i + a];
        m++;
    }
    if (m < 3)
        return MV_ERR;
    for (a = 0; a < 3; a++)
        c[a] /= m;
    memset(S, 0, sizeof(S));
    for (i = 0; i < n; i++) {
        double e[3];
        if (mask && !mask[i])
            continue;
        for (a = 0; a < 3; a++)
            e[a] = pts[3 * i + a] - c[a];
        for (a = 0; a < 3; a++)
            for (b = 0; b < 3; b++)
                S[a * 3 + b] += e[a] * e[b];
    }
    if (mv_nullvec(nvec, S, 3, 3) != MV_OK)
        return MV_ERR;
    for (a = 0; a < 3; a++)
        plane[a] = nvec[a];
    plane[3] = -mv_dot(nvec, c, 3);
    pl_canon(plane);
    return MV_OK;
}

static int pl_count(const double *pts, int n, const double plane[4],
                    double tol, unsigned char *mask)
{
    int i, votes = 0;
    for (i = 0; i < n; i++) {
        double r = fabs(mv_dot(plane, pts + 3 * i, 3) + plane[3]);
        /* admit on the comparison, so a NaN residual is rejected */
        int in = r <= tol;
        if (mask)
            mask[i] = (unsigned char)in;
        votes += in;
    }
    return votes;
}

int mv_plane_ransac(const double *pts, int n, int iters, double tol,
                    unsigned seed, double plane[4],
                    unsigned char *inliers)
{
    unsigned long long s;
    double best_pl[4];
    int best = -1, trial;
    unsigned char *mask = inliers;

    if (!pts || n < 3 || iters < 1 || !(tol > 0.0))
        return MV_ERR;
    if (!mask) {
        mask = (unsigned char *)malloc((size_t)n);
        if (!mask)
            return MV_ERR;
    }
    s = (unsigned long long)seed * 6364136223846793005ULL
        + 1442695040888963407ULL;
    for (trial = 0; trial < iters; trial++) {
        double e1[3], e2[3], nvec[3], cand[4];
        int pick[3], np = 0, a, votes;
        while (np < 3) {
            int c = (int)(pl_next(&s) % (unsigned long long)n);
            int k, dup = 0;
            for (k = 0; k < np; k++)
                if (pick[k] == c)
                    dup = 1;
            if (!dup)
                pick[np++] = c;
        }
        for (a = 0; a < 3; a++) {
            e1[a] = pts[3 * pick[1] + a] - pts[3 * pick[0] + a];
            e2[a] = pts[3 * pick[2] + a] - pts[3 * pick[0] + a];
        }
        mv_cross3(nvec, e1, e2);
        /* near-collinear sample: relative area gate, skip the trial */
        if (!(mv_norm(nvec, 3)
              > 1e-12 * mv_norm(e1, 3) * mv_norm(e2, 3)))
            continue;
        if (mv_normalize(nvec, 3) != MV_OK)
            continue;
        for (a = 0; a < 3; a++)
            cand[a] = nvec[a];
        cand[3] = -mv_dot(nvec, pts + 3 * pick[0], 3);
        votes = pl_count(pts, n, cand, tol, NULL);
        if (votes > best) {
            best = votes;
            memcpy(best_pl, cand, sizeof(best_pl));
        }
    }
    if (best < 3) {
        if (!inliers)
            free(mask);
        return MV_ERR;
    }
    /* refit on the consensus set; the output mask is recomputed against
     * the refit plane so plane[] and inliers[] are consistent */
    pl_count(pts, n, best_pl, tol, mask);
    pl_canon(best_pl);
    if (pl_fit_ls(pts, n, mask, best_pl) == MV_OK)
        pl_count(pts, n, best_pl, tol, mask);
    memcpy(plane, best_pl, sizeof(best_pl));
    if (!inliers)
        free(mask);
    return MV_OK;
}

int mv_planes_extract(const double *pts, int n, int maxplanes,
                      int min_inliers, int iters, double tol,
                      unsigned seed, double *planes, int *nplanes,
                      int *labels)
{
    double *sub;
    int *idx;
    unsigned char *mask;
    int i, remaining, pl;

    *nplanes = 0;
    if (!pts || n < 1 || maxplanes < 1 || iters < 1 || !(tol > 0.0))
        return MV_ERR;
    if (min_inliers < 3)
        min_inliers = 3;
    for (i = 0; i < n; i++)
        labels[i] = -1;
    sub = (double *)malloc((size_t)n * 3 * sizeof(double));
    idx = (int *)malloc((size_t)n * sizeof(int));
    mask = (unsigned char *)malloc((size_t)n);
    if (!sub || !idx || !mask) {
        free(sub);
        free(idx);
        free(mask);
        return MV_ERR;
    }
    for (i = 0; i < n; i++)
        idx[i] = i;
    remaining = n;
    for (pl = 0; pl < maxplanes; pl++) {
        double plane[4];
        unsigned sd = seed + 0x9E3779B9u * (unsigned)(pl + 1);
        int votes = 0, kept;
        if (remaining < min_inliers)
            break;
        for (i = 0; i < remaining; i++)
            memcpy(sub + 3 * i, pts + 3 * idx[i], 3 * sizeof(double));
        if (mv_plane_ransac(sub, remaining, iters, tol, sd, plane,
                            mask) != MV_OK)
            break;
        for (i = 0; i < remaining; i++)
            votes += mask[i];
        if (votes < min_inliers)
            break;
        for (i = 0; i < remaining; i++)
            if (mask[i])
                labels[idx[i]] = *nplanes;
        memcpy(planes + 4 * *nplanes, plane, 4 * sizeof(double));
        (*nplanes)++;
        kept = 0;
        for (i = 0; i < remaining; i++)
            if (!mask[i])
                idx[kept++] = idx[i];
        remaining = kept;
    }
    free(sub);
    free(idx);
    free(mask);
    return MV_OK;
}
