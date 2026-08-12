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
        f->scale = 1.0; /* single level: full resolution */
        feat_descriptor(f->desc, img, w, h, f->u, f->v, theta);
    }

    free(cand);
    free(gx);
    free(gy);
    free(resp);
    return nout;
}

/* ---- multi-scale detection ---------------------------------------- */

/* 3x3 binomial blur ([1 2 1]/4 separable, edges clamped): the
 * prefilter before each sqrt(2) downsample.  tmp holds the horizontal
 * pass; dst, tmp, src are all w x h and dst may not alias src. */
static void feat_blur3(unsigned char *dst, unsigned char *tmp,
                       const unsigned char *src, int w, int h)
{
    int x, y;
    for (y = 0; y < h; y++)
        for (x = 0; x < w; x++) {
            int xm = (x > 0) ? x - 1 : 0;
            int xp = (x < w - 1) ? x + 1 : w - 1;
            const unsigned char *r = src + (size_t)y * w;
            tmp[(size_t)y * w + x] =
                (unsigned char)((r[xm] + 2 * r[x] + r[xp] + 2) >> 2);
        }
    for (y = 0; y < h; y++)
        for (x = 0; x < w; x++) {
            int ym = (y > 0) ? y - 1 : 0;
            int yp = (y < h - 1) ? y + 1 : h - 1;
            dst[(size_t)y * w + x] =
                (unsigned char)((tmp[(size_t)ym * w + x]
                                 + 2 * tmp[(size_t)y * w + x]
                                 + tmp[(size_t)yp * w + x] + 2) >> 2);
        }
}

/* Bilinear resample: dst pixel (x, y) samples src at (x*r, y*r), so a
 * dst coordinate maps to level-0 through the product of the per-step
 * ratios -- exactly the absolute scale bookkeeping the caller keeps. */
static void feat_resample(unsigned char *dst, int dw, int dh,
                          const unsigned char *src, int sw, int sh,
                          double r)
{
    int x, y;
    for (y = 0; y < dh; y++)
        for (x = 0; x < dw; x++) {
            double g = feat_bilinear(src, sw, sh, x * r, y * r);
            int gi = (int)(g + 0.5);
            if (gi < 0) gi = 0;
            if (gi > 255) gi = 255;
            dst[(size_t)y * dw + x] = (unsigned char)gi;
        }
}

/* Merged-pyramid ordering: response descending; ties broken by scale
 * (finer first), then position, so the selection is deterministic. */
static int feat_ms_cmp(const void *a, const void *b)
{
    const mv_feature *fa = (const mv_feature *)a;
    const mv_feature *fb = (const mv_feature *)b;
    if (fa->resp > fb->resp) return -1;
    if (fa->resp < fb->resp) return 1;
    if (fa->scale < fb->scale) return -1;
    if (fa->scale > fb->scale) return 1;
    if (fa->v < fb->v) return -1;
    if (fa->v > fb->v) return 1;
    if (fa->u < fb->u) return -1;
    if (fa->u > fb->u) return 1;
    return 0;
}

int mv_feat_detect_ms(mv_feature *out, int maxf, const unsigned char *img,
                      int w, int h, int nlevels)
{
    mv_feature *all;
    unsigned char *cur, *buf, *tmp;
    double s = 1.0;
    int cw = w, ch = h;
    int nall = 0, nout, L;

    if (!out || !img || maxf <= 0 || nlevels < 1 || nlevels > 12
        || w < 2 * FEAT_MARGIN + 2 || h < 2 * FEAT_MARGIN + 2)
        return MV_ERR;
    if (nlevels == 1)
        return mv_feat_detect(out, maxf, img, w, h);

    all = (mv_feature *)malloc((size_t)nlevels * maxf
                               * sizeof(mv_feature));
    cur = (unsigned char *)malloc((size_t)w * h);
    buf = (unsigned char *)malloc((size_t)w * h);
    tmp = (unsigned char *)malloc((size_t)w * h);
    if (!all || !cur || !buf || !tmp) {
        free(all);
        free(cur);
        free(buf);
        free(tmp);
        return MV_ERR;
    }
    memcpy(cur, img, (size_t)w * h);

    for (L = 0; L < nlevels; L++) {
        int n, i;
        if (cw < 2 * FEAT_MARGIN + 2 || ch < 2 * FEAT_MARGIN + 2)
            break; /* image outgrown by the descriptor footprint */
        n = mv_feat_detect(all + nall, maxf, cur, cw, ch);
        if (n < 0) {
            free(all);
            free(cur);
            free(buf);
            free(tmp);
            return MV_ERR;
        }
        /* positions back to level-0 pixels; descriptor stays as
         * sampled on this level (that IS the detection scale) */
        for (i = 0; i < n; i++) {
            all[nall + i].u *= s;
            all[nall + i].v *= s;
            all[nall + i].scale = s;
        }
        nall += n;
        if (L + 1 < nlevels) {
            /* absolute scales from the original dims, so rounding
             * never accumulates across levels */
            double sn = pow(sqrt(2.0), L + 1);
            int nw = (int)(w / sn), nh = (int)(h / sn);
            feat_blur3(buf, tmp, cur, cw, ch);
            feat_resample(cur, nw, nh, buf, cw, ch, sn / s);
            cw = nw;
            ch = nh;
            s = sn;
        }
    }

    if (nall > 0)
        qsort(all, (size_t)nall, sizeof(mv_feature), feat_ms_cmp);
    nout = (nall < maxf) ? nall : maxf;
    memcpy(out, all, (size_t)nout * sizeof(mv_feature));

    free(all);
    free(cur);
    free(buf);
    free(tmp);
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
