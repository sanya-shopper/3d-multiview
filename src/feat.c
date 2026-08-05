#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mv/core.h"
#include "mv/mat.h"
#include "mv/feat.h"

/* Tuning constants; margins below depend on them. */
#define FEAT_WIN_R 2    /* structure-matrix window radius (5x5) */
#define FEAT_NMS_R 5    /* non-maximum suppression radius */
#define FEAT_ORI_R 7    /* orientation window radius */
#define FEAT_ORI_SIG 3.5
#define FEAT_STEP 2.0   /* descriptor grid spacing in pixels */
/* descriptor reach: 3.5 * FEAT_STEP * sqrt(2) ~ 9.9, +0.5 sub-pixel,
 * +1 bilinear -> 12 px border margin */
#define FEAT_MARGIN 12

static double feat_bilinear(const unsigned char *img, int w, int h,
                            double x, double y)
{
    int ix, iy;
    double fx, fy;
    const unsigned char *p;
    if (x < 0.0) x = 0.0;
    if (y < 0.0) y = 0.0;
    if (x > w - 1.001) x = w - 1.001;
    if (y > h - 1.001) y = h - 1.001;
    ix = (int)x;
    iy = (int)y;
    fx = x - ix;
    fy = y - iy;
    p = img + iy * w + ix;
    return (1 - fx) * (1 - fy) * p[0] + fx * (1 - fy) * p[1]
         + (1 - fx) * fy * p[w] + fx * fy * p[w + 1];
}

/* 1-D quadratic peak offset from samples (m, c, p) at -1, 0, +1;
 * clamped to [-0.5, 0.5]. */
static double feat_subpix(double m, double c, double p)
{
    double den = m - 2.0 * c + p;
    double d;
    if (fabs(den) < 1e-12)
        return 0.0;
    d = 0.5 * (m - p) / den;
    if (d < -0.5) d = -0.5;
    if (d > 0.5) d = 0.5;
    return d;
}

/* Gaussian-weighted mean gradient orientation at (x, y). */
static double feat_orientation(const float *gx, const float *gy, int w,
                               int x, int y)
{
    double sx = 0.0, sy = 0.0;
    int dx, dy;
    for (dy = -FEAT_ORI_R; dy <= FEAT_ORI_R; dy++)
        for (dx = -FEAT_ORI_R; dx <= FEAT_ORI_R; dx++) {
            double wgt = exp(-(double)(dx * dx + dy * dy)
                             / (2.0 * FEAT_ORI_SIG * FEAT_ORI_SIG));
            int o = (y + dy) * w + (x + dx);
            sx += wgt * gx[o];
            sy += wgt * gy[o];
        }
    return atan2(sy, sx);
}

/* Oriented, illumination-normalized 8x8 patch descriptor. */
static void feat_descriptor(double desc[64], const unsigned char *img,
                            int w, int h, double u, double v, double theta)
{
    double c = cos(theta), s = sin(theta);
    double mean = 0.0, var = 0.0, sd;
    int i, j, k = 0;
    for (i = 0; i < 8; i++)
        for (j = 0; j < 8; j++) {
            double lx = (j - 3.5) * FEAT_STEP, ly = (i - 3.5) * FEAT_STEP;
            double sx = u + c * lx - s * ly;
            double sy = v + s * lx + c * ly;
            desc[k++] = feat_bilinear(img, w, h, sx, sy);
        }
    for (k = 0; k < 64; k++)
        mean += desc[k];
    mean /= 64.0;
    for (k = 0; k < 64; k++) {
        desc[k] -= mean;
        var += desc[k] * desc[k];
    }
    sd = sqrt(var / 64.0);
    if (sd < 1e-9) { /* flat patch: descriptor carries no signal */
        for (k = 0; k < 64; k++)
            desc[k] = 0.0;
        return;
    }
    for (k = 0; k < 64; k++)
        desc[k] /= sd;
}

/* Candidate corner prior to refinement. */
typedef struct {
    int x, y;
    float resp;
} feat_cand;

static int feat_cand_cmp(const void *a, const void *b)
{
    const feat_cand *ca = (const feat_cand *)a;
    const feat_cand *cb = (const feat_cand *)b;
    if (ca->resp > cb->resp) return -1;
    if (ca->resp < cb->resp) return 1;
    /* deterministic tie-break by scan order */
    if (ca->y != cb->y) return ca->y - cb->y;
    return ca->x - cb->x;
}

int mv_feat_detect(mv_feature *out, int maxf, const unsigned char *img,
                   int w, int h)
{
    float *gx, *gy, *resp;
    feat_cand *cand;
    int ncand = 0, nout = 0;
    int x, y, i;
    float maxresp = 0.0f, thr;

    if (!out || !img || maxf <= 0 || w < 2 * FEAT_MARGIN + 2
        || h < 2 * FEAT_MARGIN + 2)
        return MV_ERR;

    gx = (float *)calloc((size_t)w * h, sizeof(float));
    gy = (float *)calloc((size_t)w * h, sizeof(float));
    resp = (float *)calloc((size_t)w * h, sizeof(float));
    if (!gx || !gy || !resp) {
        free(gx);
        free(gy);
        free(resp);
        return MV_ERR;
    }

    /* Sobel gradients (interior only; borders stay zero) */
    for (y = 1; y < h - 1; y++)
        for (x = 1; x < w - 1; x++) {
            const unsigned char *p = img + y * w + x;
            gx[y * w + x] = (float)((p[-w + 1] + 2 * p[1] + p[w + 1])
                                    - (p[-w - 1] + 2 * p[-1] + p[w - 1]));
            gy[y * w + x] = (float)((p[w - 1] + 2 * p[w] + p[w + 1])
                                    - (p[-w - 1] + 2 * p[-w] + p[-w + 1]));
        }

    /* Shi-Tomasi response: min eigenvalue of the 5x5 gradient-moment
     * matrix.  Only where the descriptor footprint fits. */
    for (y = FEAT_MARGIN; y < h - FEAT_MARGIN; y++)
        for (x = FEAT_MARGIN; x < w - FEAT_MARGIN; x++) {
            double sxx = 0.0, sxy = 0.0, syy = 0.0, tr, det2;
            int dx, dy;
            float r;
            for (dy = -FEAT_WIN_R; dy <= FEAT_WIN_R; dy++)
                for (dx = -FEAT_WIN_R; dx <= FEAT_WIN_R; dx++) {
                    double a = gx[(y + dy) * w + (x + dx)];
                    double b = gy[(y + dy) * w + (x + dx)];
                    sxx += a * a;
                    sxy += a * b;
                    syy += b * b;
                }
            tr = 0.5 * (sxx + syy);
            det2 = sqrt(0.25 * (sxx - syy) * (sxx - syy) + sxy * sxy);
            r = (float)(tr - det2);
            resp[y * w + x] = r;
            if (r > maxresp)
                maxresp = r;
        }
    thr = 0.01f * maxresp;

    /* Non-maximum suppression: strict local maxima over FEAT_NMS_R.
     * Two passes: count, then fill. */
    cand = NULL;
    {
        int pass;
        for (pass = 0; pass < 2; pass++) {
            int n = 0;
            for (y = FEAT_MARGIN; y < h - FEAT_MARGIN; y++)
                for (x = FEAT_MARGIN; x < w - FEAT_MARGIN; x++) {
                    float r = resp[y * w + x];
                    int dx, dy, ismax = 1;
                    if (r <= thr)
                        continue;
                    for (dy = -FEAT_NMS_R; dy <= FEAT_NMS_R && ismax; dy++)
                        for (dx = -FEAT_NMS_R; dx <= FEAT_NMS_R; dx++) {
                            int nx = x + dx, ny = y + dy;
                            float q;
                            if ((dx == 0 && dy == 0) || nx < 0 || ny < 0
                                || nx >= w || ny >= h)
                                continue;
                            q = resp[ny * w + nx];
                            /* ties resolved by scan order */
                            if (q > r || (q == r && (dy < 0
                                          || (dy == 0 && dx < 0)))) {
                                ismax = 0;
                                break;
                            }
                        }
                    if (!ismax)
                        continue;
                    if (pass == 1) {
                        cand[n].x = x;
                        cand[n].y = y;
                        cand[n].resp = r;
                    }
                    n++;
                }
            if (pass == 0) {
                ncand = n;
                if (ncand == 0)
                    break;
                cand = (feat_cand *)malloc((size_t)ncand
                                           * sizeof(feat_cand));
                if (!cand) {
                    free(gx);
                    free(gy);
                    free(resp);
                    return MV_ERR;
                }
            }
        }
    }

    if (ncand > 0)
        qsort(cand, (size_t)ncand, sizeof(feat_cand), feat_cand_cmp);

    for (i = 0; i < ncand && nout < maxf; i++) {
        int cx = cand[i].x, cy = cand[i].y, o = cy * w + cx;
        double du = feat_subpix(resp[o - 1], resp[o], resp[o + 1]);
        double dv = feat_subpix(resp[o - w], resp[o], resp[o + w]);
        double theta = feat_orientation(gx, gy, w, cx, cy);
        mv_feature *f = &out[nout++];
        f->u = cx + du;
        f->v = cy + dv;
        f->resp = cand[i].resp;
        feat_descriptor(f->desc, img, w, h, f->u, f->v, theta);
    }

    free(cand);
    free(gx);
    free(gy);
    free(resp);
    return nout;
}

static double feat_ssd(const double *a, const double *b)
{
    double s = 0.0;
    int k;
    for (k = 0; k < 64; k++) {
        double d = a[k] - b[k];
        s += d * d;
    }
    return s;
}

int mv_feat_match(int *idx2, const mv_feature *f1, int n1,
                  const mv_feature *f2, int n2, double ratio)
{
    int *best2; /* best2[j] = nearest f1 index for f2[j] (mutuality) */
    int i, j, count = 0;

    if (!idx2 || !f1 || !f2 || n1 <= 0 || n2 <= 0 || ratio <= 0.0)
        return MV_ERR;

    for (i = 0; i < n1; i++)
        idx2[i] = -1;
    if (n2 < 2)
        return 0; /* no second-nearest for the ratio test */

    best2 = (int *)malloc((size_t)n2 * sizeof(int));
    if (!best2)
        return MV_ERR;
    for (j = 0; j < n2; j++) {
        double bd = HUGE_VAL;
        int bi = -1;
        for (i = 0; i < n1; i++) {
            double d = feat_ssd(f1[i].desc, f2[j].desc);
            if (d < bd) {
                bd = d;
                bi = i;
            }
        }
        best2[j] = bi;
    }

    for (i = 0; i < n1; i++) {
        double d1 = HUGE_VAL, d2 = HUGE_VAL;
        int bj = -1;
        for (j = 0; j < n2; j++) {
            double d = feat_ssd(f1[i].desc, f2[j].desc);
            if (d < d1) {
                d2 = d1;
                d1 = d;
                bj = j;
            } else if (d < d2) {
                d2 = d;
            }
        }
        if (bj >= 0 && d1 < ratio * d2 && best2[bj] == i) {
            idx2[i] = bj;
            count++;
        }
    }

    free(best2);
    return count;
}
