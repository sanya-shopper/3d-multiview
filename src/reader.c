#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mv/core.h"
#include "mv/mat.h"
#include "mv/calib.h"
#include "mv/pattern.h"
#include "mv/reader.h"

#define MAXCAND 6144
#define RDBG(...) do { if (getenv("MV_READER_DEBUG")) \
        fprintf(stderr, __VA_ARGS__); } while (0)
#define GRIDN 48          /* lattice bookkeeping extent */
#define GOFF 24

/* ---------------- saddle detection ---------------------------------- */

#define SADR 5

struct offs {
    int n;
    int dx[128], dy[128];
    double t[4][128]; /* template value per phase */
};

/* Frames at or above this pixel count take the downscale-first ladder
 * (decode on a half-resolution working image, then re-localize corners
 * at full resolution): detection cost scales with pixels, and full-HD
 * frames put the validated ~20 px cell regime at half resolution. */
#define BIGFRAME_PX 1500000L

static void build_offsets(struct offs *o)
{
    const double ph[4] = { 0.0, 0.3927, 0.7854, 1.1781 };
    int dx, dy, p;
    o->n = 0;
    for (dy = -SADR; dy <= SADR; dy++)
        for (dx = -SADR; dx <= SADR; dx++) {
            if (dx * dx + dy * dy > SADR * SADR || (dx == 0 && dy == 0))
                continue;
            o->dx[o->n] = dx;
            o->dy[o->n] = dy;
            for (p = 0; p < 4; p++) {
                double c = cos(ph[p]), s = sin(ph[p]);
                double a = c * dx + s * dy, b = -s * dx + c * dy;
                o->t[p][o->n] = (a * b > 0.0) ? 1.0 : (a * b < 0.0 ? -1.0
                                                                   : 0.0);
            }
            o->n++;
        }
}

struct cand {
    double u, v, resp;
    int used;
};

static const struct offs *offsets(void)
{
    static struct offs o;
    static int o_init = 0;
    if (!o_init) {
        build_offsets(&o);
        o_init = 1;
    }
    return &o;
}

static int detect_saddles(struct cand *cd, const unsigned char *img,
                          int w, int h, double relthr)
{
    const struct offs *o = offsets();
    float *R;
    short *acc; /* four phase-accumulator rows */
    double rmax = 0.0;
    int x, y, k, p, n = 0;

    if (w < 2 * SADR + 4 || h < 2 * SADR + 4)
        return 0;
    R = (float *)calloc((size_t)w * h, sizeof(float));
    acc = (short *)malloc(4 * (size_t)w * sizeof(short));
    if (!R || !acc) {
        free(R);
        free(acc);
        return 0;
    }
    /* Response pass, restructured for throughput: accumulate each
     * template tap across a whole output row (contiguous loads the
     * compiler vectorizes) instead of gathering all taps per pixel.
     * Every sum is an exact integer with |s| <= o->n * 255 < 32767,
     * so int16 accumulation reproduces the per-pixel double formula
     * bit for bit -- detection behavior is unchanged. */
    for (y = SADR + 1; y < h - SADR - 1; y++) {
        const int x0 = SADR + 1, x1 = w - SADR - 1;
        memset(acc, 0, 4 * (size_t)w * sizeof(short));
        for (p = 0; p < 4; p++) {
            short *ap = acc + (size_t)p * w;
            for (k = 0; k < o->n; k++) {
                const unsigned char *src =
                    img + (size_t)(y + o->dy[k]) * w + o->dx[k];
                if (o->t[p][k] > 0.0) {
                    for (x = x0; x < x1; x++)
                        ap[x] = (short)(ap[x] + src[x]);
                } else if (o->t[p][k] < 0.0) {
                    for (x = x0; x < x1; x++)
                        ap[x] = (short)(ap[x] - src[x]);
                }
            }
        }
        for (x = SADR + 1; x < w - SADR - 1; x++) {
            int best = 0;
            double b;
            for (p = 0; p < 4; p++) {
                int s = acc[(size_t)p * w + x];
                if (s < 0)
                    s = -s;
                if (s > best)
                    best = s;
            }
            b = (double)best / o->n;
            R[(size_t)y * w + x] = (float)b;
            if (b > rmax)
                rmax = b;
        }
    }

    /* NMS + threshold + sub-pixel; if the pool fills, keep strongest */
    for (y = SADR + 2; y < h - SADR - 2; y++)
        for (x = SADR + 2; x < w - SADR - 2; x++) {
            float r0 = R[y * w + x];
            int i, j, peak = 1;
            double dxs, dys, den;
            if (r0 < relthr * rmax || r0 < 6.0)
                continue;
            for (j = -2; j <= 2 && peak; j++)
                for (i = -2; i <= 2 && peak; i++) {
                    float rn;
                    if (!i && !j)
                        continue;
                    rn = R[(y + j) * w + (x + i)];
                    /* ties (exact plateaus on clean synthetic input)
                     * break to the scan-earlier pixel, else a plateau
                     * emits one candidate per pixel */
                    if (rn > r0 || (rn == r0 && (j < 0 || (j == 0
                                                           && i < 0))))
                        peak = 0;
                }
            if (!peak)
                continue;
            den = R[y * w + x - 1] - 2.0 * r0 + R[y * w + x + 1];
            dxs = (den < 0.0)
                  ? 0.5 * (R[y * w + x - 1] - R[y * w + x + 1]) / den : 0.0;
            den = R[(y - 1) * w + x] - 2.0 * r0 + R[(y + 1) * w + x];
            dys = (den < 0.0)
                  ? 0.5 * (R[(y - 1) * w + x] - R[(y + 1) * w + x]) / den
                  : 0.0;
            if (fabs(dxs) > 1.0)
                dxs = 0.0;
            if (fabs(dys) > 1.0)
                dys = 0.0;
            if (n < MAXCAND) {
                cd[n].u = x + dxs;
                cd[n].v = y + dys;
                cd[n].resp = r0;
                cd[n].used = 0;
                n++;
            } else {
                int wk = 0, q;
                for (q = 1; q < n; q++)
                    if (cd[q].resp < cd[wk].resp)
                        wk = q;
                if (r0 > cd[wk].resp) {
                    cd[wk].u = x + dxs;
                    cd[wk].v = y + dys;
                    cd[wk].resp = r0;
                    cd[wk].used = 0;
                }
            }
        }
    free(acc);
    free(R);
    return n;
}

/* Gradient-orthogonality sub-pixel corner solve (Foerstner / the
 * cornerSubPix normal equations): at a checkerboard corner q the image
 * gradient g(p) at every window pixel p lies on an edge through q or
 * is ~zero, so g(p) . (p - q) = 0; q solves the weighted least-squares
 * system (sum w g g^T) q = sum w g g^T p.  Iterates with window
 * recentering; returns 0 (leaving q untouched) when the window leaves
 * the image, the normal matrix degenerates, or the solution runs away
 * (> 2 px from the start -- then the initial estimate was not on a
 * corner and must be kept). */
#define FSR 5  /* half-window, px; < half a cell so the fine tier's
                * center dots stay outside the window */
#define FPAD 3 /* patch margin: gradient (1) + two blur passes (2) */
#define FPW (2 * (FSR + FPAD) + 1)

static int foerstner_step(const unsigned char *img, int w, int h,
                          double *u, double *v)
{
    static double wtab[2 * FSR + 1][2 * FSR + 1];
    static int wtab_init = 0;
    double qu = *u, qv = *v;
    int it, i, j;

    if (!wtab_init) {
        const double sig = 0.5 * FSR;
        for (j = -FSR; j <= FSR; j++)
            for (i = -FSR; i <= FSR; i++)
                wtab[j + FSR][i + FSR] =
                    exp(-(i * i + j * j) / (2.0 * sig * sig));
        wtab_init = 1;
    }
    for (it = 0; it < 8; it++) {
        int cx = (int)floor(qu + 0.5), cy = (int)floor(qv + 0.5);
        double P[FPW][FPW], Q[FPW][FPW];
        double A0 = 0.0, A1 = 0.0, A3 = 0.0, b0 = 0.0, b1 = 0.0;
        double det, nu, nv, du, dv;
        if (cx - FSR - FPAD < 0 || cy - FSR - FPAD < 0
            || cx + FSR + FPAD >= w || cy + FSR + FPAD >= h)
            return 0;
        /* Smooth the patch (two separable [1 2 1]/4 passes) before
         * taking gradients: sampled/rendered edges are nearly hard
         * steps whose sub-pixel position lives in a single transition
         * pixel, and smoothing spreads it into gradients the normal
         * equations can use.  The kernel is symmetric, so the solve
         * stays unbiased. */
        for (j = 0; j < FPW; j++)
            for (i = 0; i < FPW; i++)
                P[j][i] = img[(size_t)(cy - FSR - FPAD + j) * w
                              + cx - FSR - FPAD + i];
        for (j = 0; j < 2; j++) {
            int a, b;
            for (a = 0; a < FPW; a++)
                for (b = 1; b < FPW - 1; b++)
                    Q[a][b] = 0.25 * P[a][b - 1] + 0.5 * P[a][b]
                              + 0.25 * P[a][b + 1];
            for (b = 1; b < FPW - 1; b++)
                for (a = 1; a < FPW - 1; a++)
                    P[a][b] = 0.25 * Q[a - 1][b] + 0.5 * Q[a][b]
                              + 0.25 * Q[a + 1][b];
        }
        for (j = -FSR; j <= FSR; j++)
            for (i = -FSR; i <= FSR; i++) {
                int pi = i + FSR + FPAD, pj = j + FSR + FPAD;
                double gx = 0.5 * (P[pj][pi + 1] - P[pj][pi - 1]);
                double gy = 0.5 * (P[pj + 1][pi] - P[pj - 1][pi]);
                double wt = wtab[j + FSR][i + FSR];
                double gxx = wt * gx * gx, gxy = wt * gx * gy;
                double gyy = wt * gy * gy;
                A0 += gxx;
                A1 += gxy;
                A3 += gyy;
                b0 += gxx * (cx + i) + gxy * (cy + j);
                b1 += gxy * (cx + i) + gyy * (cy + j);
            }
        det = A0 * A3 - A1 * A1;
        if (!(det > 1e-6))
            return 0;
        nu = (A3 * b0 - A1 * b1) / det;
        nv = (A0 * b1 - A1 * b0) / det;
        du = nu - qu;
        dv = nv - qv;
        qu = nu;
        qv = nv;
        if (du * du + dv * dv < 1e-6)
            break;
    }
    if ((qu - *u) * (qu - *u) + (qv - *v) * (qv - *v) > 4.0)
        return 0;
    *u = qu;
    *v = qv;
    return 1;
}

/* Re-localize one corner at full resolution after a decode on a
 * downscaled working image: evaluate the same saddle response on a
 * small pixel grid around the scaled-up estimate, snap to the peak
 * (with the detector's quadratic sub-pixel step), then polish with the
 * gradient-orthogonality solve above.  The estimate error after a
 * 2x/4x downscale is ~1 px, well inside the search radius and far from
 * the >= 15 px saddle spacing, so the peak is the same physical
 * corner.  Leaves the estimate untouched when the window would cross
 * the image border. */
#define RFR 3 /* search radius, px */

static void refine_saddle(const unsigned char *img, int w, int h,
                          double *u, double *v)
{
    const struct offs *o = offsets();
    double Rw[2 * RFR + 3][2 * RFR + 3];
    int cx = (int)floor(*u + 0.5), cy = (int)floor(*v + 0.5);
    int i, j, k, p, bi = 0, bj = 0;
    double bbest = -1.0, dxs, dys, den, r0;

    if (cx - RFR - 1 - SADR < 0 || cy - RFR - 1 - SADR < 0
        || cx + RFR + 1 + SADR >= w || cy + RFR + 1 + SADR >= h)
        return;
    for (j = -RFR - 1; j <= RFR + 1; j++)
        for (i = -RFR - 1; i <= RFR + 1; i++) {
            const unsigned char *c0 = img + (size_t)(cy + j) * w + cx + i;
            int best = 0;
            for (p = 0; p < 4; p++) {
                int s = 0;
                for (k = 0; k < o->n; k++) {
                    int g = c0[o->dy[k] * w + o->dx[k]];
                    if (o->t[p][k] > 0.0)
                        s += g;
                    else if (o->t[p][k] < 0.0)
                        s -= g;
                }
                if (s < 0)
                    s = -s;
                if (s > best)
                    best = s;
            }
            Rw[j + RFR + 1][i + RFR + 1] = (double)best / o->n;
        }
    for (j = -RFR; j <= RFR; j++)
        for (i = -RFR; i <= RFR; i++)
            if (Rw[j + RFR + 1][i + RFR + 1] > bbest) {
                bbest = Rw[j + RFR + 1][i + RFR + 1];
                bi = i;
                bj = j;
            }
    r0 = Rw[bj + RFR + 1][bi + RFR + 1];
    den = Rw[bj + RFR + 1][bi + RFR] - 2.0 * r0
          + Rw[bj + RFR + 1][bi + RFR + 2];
    dxs = (den < 0.0)
          ? 0.5 * (Rw[bj + RFR + 1][bi + RFR]
                   - Rw[bj + RFR + 1][bi + RFR + 2]) / den : 0.0;
    den = Rw[bj + RFR][bi + RFR + 1] - 2.0 * r0
          + Rw[bj + RFR + 2][bi + RFR + 1];
    dys = (den < 0.0)
          ? 0.5 * (Rw[bj + RFR][bi + RFR + 1]
                   - Rw[bj + RFR + 2][bi + RFR + 1]) / den : 0.0;
    if (fabs(dxs) > 1.0)
        dxs = 0.0;
    if (fabs(dys) > 1.0)
        dys = 0.0;
    *u = cx + bi + dxs;
    *v = cy + bj + dys;
    foerstner_step(img, w, h, u, v);
}

/* Image-pixel size of one pattern cell at the pattern center under H
 * (the smaller of the two axis steps): the gate for the sub-pixel
 * polish, whose FSR window must stay inside a half cell. */
static double cell_size_px(const double H[9], double cx, double cy,
                           double cell)
{
    double W0 = H[6] * cx + H[7] * cy + H[8];
    double Wx = H[6] * (cx + cell) + H[7] * cy + H[8];
    double Wy = H[6] * cx + H[7] * (cy + cell) + H[8];
    double u0, v0, ux, vx, uy, vy, dx, dy;
    if (fabs(W0) < 1e-12 || fabs(Wx) < 1e-12 || fabs(Wy) < 1e-12)
        return 0.0;
    u0 = (H[0] * cx + H[1] * cy + H[2]) / W0;
    v0 = (H[3] * cx + H[4] * cy + H[5]) / W0;
    ux = (H[0] * (cx + cell) + H[1] * cy + H[2]) / Wx;
    vx = (H[3] * (cx + cell) + H[4] * cy + H[5]) / Wx;
    uy = (H[0] * cx + H[1] * (cy + cell) + H[2]) / Wy;
    vy = (H[3] * cx + H[4] * (cy + cell) + H[5]) / Wy;
    dx = sqrt((ux - u0) * (ux - u0) + (vx - v0) * (vx - v0));
    dy = sqrt((uy - u0) * (uy - u0) + (vy - v0) * (vy - v0));
    return dx < dy ? dx : dy;
}

/* ---------------- grid growth ---------------------------------------- */

struct lattice {
    int have[GRIDN][GRIDN];
    double p[GRIDN][GRIDN][2];
};

static int grow_from_seed(struct lattice *L, struct cand *cd, int n,
                          int seed)
{
    int qi[GRIDN * GRIDN][2], qh = 0, qt = 0;
    double d1[2], d2[2];
    int k, count = 0;

    memset(L->have, 0, sizeof(L->have));
    for (k = 0; k < n; k++)
        cd[k].used = 0;

    /* d1: nearest neighbor; d2: nearest at 55-125 degrees to d1 */
    {
        double su = cd[seed].u, sv = cd[seed].v;
        int b1 = -1, b2 = -1;
        double bd1 = 1e18, bd2 = 1e18;
        for (k = 0; k < n; k++) {
            double du, dv, d;
            if (k == seed)
                continue;
            du = cd[k].u - su;
            dv = cd[k].v - sv;
            d = du * du + dv * dv;
            if (d < 64.0 || d > 6400.0)
                continue;
            if (d < bd1) {
                bd1 = d;
                b1 = k;
            }
        }
        if (b1 < 0)
            return 0;
        d1[0] = cd[b1].u - su;
        d1[1] = cd[b1].v - sv;
        for (k = 0; k < n; k++) {
            double du, dv, d, dot, cross, ang;
            if (k == seed || k == b1)
                continue;
            du = cd[k].u - su;
            dv = cd[k].v - sv;
            d = du * du + dv * dv;
            if (d < 64.0 || d > 4.0 * bd1)
                continue;
            dot = du * d1[0] + dv * d1[1];
            cross = d1[0] * dv - d1[1] * du;
            ang = fabs(atan2(cross, dot));
            if (ang < 0.96 || ang > 2.18) /* outside 55..125 deg */
                continue;
            if (d < bd2) {
                bd2 = d;
                b2 = k;
            }
        }
        if (b2 < 0)
            return 0;
        d2[0] = cd[b2].u - su;
        d2[1] = cd[b2].v - sv;
        /* fix chirality: force cross(d1,d2) > 0 */
        if (d1[0] * d2[1] - d1[1] * d2[0] < 0.0) {
            d2[0] = -d2[0];
            d2[1] = -d2[1];
        }
        L->have[GOFF][GOFF] = 1;
        L->p[GOFF][GOFF][0] = su;
        L->p[GOFF][GOFF][1] = sv;
        cd[seed].used = 1;
        qi[qt][0] = GOFF;
        qi[qt][1] = GOFF;
        qt++;
        count = 1;
    }

    while (qh < qt) {
        int ci = qi[qh][0], cj = qi[qh][1];
        static const int dir[4][2] = { {1,0},{-1,0},{0,1},{0,-1} };
        int d;
        qh++;
        for (d = 0; d < 4; d++) {
            int ni = ci + dir[d][0], nj = cj + dir[d][1];
            int oi = ci - dir[d][0], oj = cj - dir[d][1];
            double pu, pv, step;
            if (ni < 1 || nj < 1 || ni >= GRIDN - 1 || nj >= GRIDN - 1)
                continue;
            if (L->have[ni][nj])
                continue;
            if (L->have[oi][oj]) {
                pu = 2.0 * L->p[ci][cj][0] - L->p[oi][oj][0];
                pv = 2.0 * L->p[ci][cj][1] - L->p[oi][oj][1];
            } else {
                pu = L->p[ci][cj][0] + dir[d][0] * d1[0] + dir[d][1] * d2[0];
                pv = L->p[ci][cj][1] + dir[d][0] * d1[1] + dir[d][1] * d2[1];
            }
            step = sqrt((pu - L->p[ci][cj][0]) * (pu - L->p[ci][cj][0])
                        + (pv - L->p[ci][cj][1]) * (pv - L->p[ci][cj][1]));
            {
                int best = -1, k2;
                double gate = 0.35 * step, bd = gate * gate;
                for (k2 = 0; k2 < n; k2++) {
                    double du, dv, dd;
                    if (cd[k2].used)
                        continue;
                    du = cd[k2].u - pu;
                    dv = cd[k2].v - pv;
                    dd = du * du + dv * dv;
                    if (dd < bd) {
                        bd = dd;
                        best = k2;
                    }
                }
                if (best >= 0) {
                    L->have[ni][nj] = 1;
                    L->p[ni][nj][0] = cd[best].u;
                    L->p[ni][nj][1] = cd[best].v;
                    cd[best].used = 1;
                    qi[qt][0] = ni;
                    qi[qt][1] = nj;
                    qt++;
                    count++;
                }
            }
        }
    }
    return count;
}

static int grow_grid(struct lattice *L, struct cand *cd, int n)
{
    static char tried[MAXCAND];
    int attempt;
    memset(tried, 0, sizeof(tried));
    for (attempt = 0; attempt < 24; attempt++) {
        int seed = -1, k, count;
        for (k = 0; k < n; k++)
            if (!tried[k] && (seed < 0 || cd[k].resp > cd[seed].resp))
                seed = k;
        if (seed < 0)
            break;
        tried[seed] = 1;
        count = grow_from_seed(L, cd, n, seed);
        if (count >= 30)
            return count;
    }
    return 0;
}

/* ---------------- sampling helpers ------------------------------------ */

static double bilin(const unsigned char *img, int w, int h, double u,
                    double v)
{
    int x = (int)u, y = (int)v;
    double fx = u - x, fy = v - y;
    const unsigned char *p;
    if (x < 0 || y < 0 || x >= w - 1 || y >= h - 1)
        return -1.0;
    p = img + y * w + x;
    return (1 - fx) * (1 - fy) * p[0] + fx * (1 - fy) * p[1]
         + (1 - fx) * fy * p[w] + fx * fy * p[w + 1];
}

static double sample3(const unsigned char *img, int w, int h, double u,
                      double v)
{
    double s = 0.0;
    int i, j, n = 0;
    for (j = -1; j <= 1; j++)
        for (i = -1; i <= 1; i++) {
            double g = bilin(img, w, h, u + i, v + j);
            if (g >= 0.0) {
                s += g;
                n++;
            }
        }
    return n ? s / n : -1.0;
}

static double med4(double a, double b, double c, double d)
{
    double v[4] = { a, b, c, d }, t;
    int i, j;
    for (i = 1; i < 4; i++)
        for (j = i; j > 0 && v[j] < v[j - 1]; j--) {
            t = v[j];
            v[j] = v[j - 1];
            v[j - 1] = t;
        }
    return 0.5 * (v[1] + v[2]);
}

/* ---------------- edge-intersection corner polish --------------------- */

static void proj_H(const double H[9], double X, double Y, double *u,
                   double *v); /* defined with the display-outline code */

/* Collect sub-pixel edge crossings along the checker edge through
 * pattern point (X0,Y0) in pattern direction (ex,ey), both ways, at
 * arc offsets lo..hi pattern px from the corner: at each spot, sample
 * the image profile perpendicular to the projected edge and take the
 * mid-level crossing nearest the prediction by linear interpolation.
 * Writes up to maxp (u,v) points; returns the count. */
static int edge_crossings(const unsigned char *img, int w, int h,
                          const double H[9], double X0, double Y0,
                          double ex, double ey, double lo, double hi,
                          double *pts, int maxp)
{
    int side, k, np = 0;

    for (side = 0; side < 2; side++) {
        double sgn = side ? -1.0 : 1.0;
        for (k = 0; k < 8 && np < maxp; k++) {
            double s = sgn * (lo + (hi - lo) * k / 7.0);
            double uc, vc, u1, v1, tx, ty, tn, px, py;
            double g[13], gmin = 255.0, gmax = 0.0, mid;
            double bestd = 1e18, rstar = 0.0;
            int m, bad = 0, found = 0;
            proj_H(H, X0 + s * ex, Y0 + s * ey, &uc, &vc);
            proj_H(H, X0 + (s + 1.0) * ex, Y0 + (s + 1.0) * ey,
                   &u1, &v1);
            tx = u1 - uc;
            ty = v1 - vc;
            tn = sqrt(tx * tx + ty * ty);
            if (tn < 1e-9)
                continue;
            px = -ty / tn;
            py = tx / tn;
            for (m = 0; m < 13 && !bad; m++) {
                double r = -2.4 + 0.4 * m;
                g[m] = bilin(img, w, h, uc + r * px, vc + r * py);
                if (g[m] < 0.0)
                    bad = 1;
                else {
                    if (g[m] < gmin)
                        gmin = g[m];
                    if (g[m] > gmax)
                        gmax = g[m];
                }
            }
            if (bad || gmax - gmin < 25.0)
                continue;
            mid = 0.5 * (gmin + gmax);
            for (m = 0; m < 12; m++) {
                double a = g[m] - mid, b = g[m + 1] - mid;
                double f, r, d;
                if (a == 0.0 && b == 0.0)
                    continue;
                if ((a <= 0.0 && b >= 0.0) || (a >= 0.0 && b <= 0.0)) {
                    f = a / (a - b);
                    r = -2.4 + 0.4 * (m + f);
                    d = fabs(r);
                    if (d < bestd) {
                        bestd = d;
                        rstar = r;
                        found = 1;
                    }
                }
            }
            if (!found)
                continue;
            pts[2 * np] = uc + rstar * px;
            pts[2 * np + 1] = vc + rstar * py;
            np++;
        }
    }
    return np;
}

/* Total-least-squares line through a point scatter: centroid plus the
 * principal direction. */
static int fit_line(const double *pts, int n, double c[2], double d[2])
{
    double mx = 0.0, my = 0.0, sxx = 0.0, sxy = 0.0, syy = 0.0, ang;
    int k;
    if (n < 8)
        return MV_ERR;
    for (k = 0; k < n; k++) {
        mx += pts[2 * k];
        my += pts[2 * k + 1];
    }
    mx /= n;
    my /= n;
    for (k = 0; k < n; k++) {
        double dx = pts[2 * k] - mx, dy = pts[2 * k + 1] - my;
        sxx += dx * dx;
        sxy += dx * dy;
        syy += dy * dy;
    }
    ang = 0.5 * atan2(2.0 * sxy, sxx - syy);
    c[0] = mx;
    c[1] = my;
    d[0] = cos(ang);
    d[1] = sin(ang);
    return MV_OK;
}

/* Refine one corner as the intersection of its two fitted edge lines.
 * Checker edges are exactly straight under a homography, so ~30
 * crossings per edge average the per-pixel sampling quantization far
 * below what any local window can reach.  Rejects (keeping *u,*v)
 * when either edge collects too few crossings, the lines are near
 * parallel, or the intersection runs away from the start (> 2.5 px:
 * then the corner was not where H says and the local estimate is the
 * safer answer). */
static int refine_corner_edges(const unsigned char *img, int w, int h,
                               const double H[9], double X0, double Y0,
                               double cell, double *u, double *v)
{
    double p1[2 * 16], p2[2 * 16];
    double c1[2], d1[2], c2[2], d2[2];
    double lo = 0.12 * cell, hi = 0.78 * cell;
    double det, a, iu, iv;
    int n1, n2;

    n1 = edge_crossings(img, w, h, H, X0, Y0, 1.0, 0.0, lo, hi, p1, 16);
    n2 = edge_crossings(img, w, h, H, X0, Y0, 0.0, 1.0, lo, hi, p2, 16);
    if (fit_line(p1, n1, c1, d1) != MV_OK
        || fit_line(p2, n2, c2, d2) != MV_OK)
        return MV_ERR;
    det = d1[0] * (-d2[1]) - (-d2[0]) * d1[1];
    if (fabs(det) < 0.34) /* < ~20 degrees apart */
        return MV_ERR;
    a = ((c2[0] - c1[0]) * (-d2[1]) - (-d2[0]) * (c2[1] - c1[1])) / det;
    iu = c1[0] + a * d1[0];
    iv = c1[1] + a * d1[1];
    if ((iu - *u) * (iu - *u) + (iv - *v) * (iv - *v) > 6.25)
        return MV_ERR;
    *u = iu;
    *v = iv;
    return MV_OK;
}

/* ---------------- main entry ------------------------------------------ */

static int read_pattern_once(mv_read_result *res, const unsigned char *img,
                             int w, int h, int allow_down2,
                             int *used_down2)
{
    static struct cand cd[MAXCAND];
    static struct lattice L;
    static signed char bit[GRIDN][GRIDN];   /* -1 unknown */
    static double base[GRIDN][GRIDN];
    int ncand, ngrid, i, j;
    double blackref = 255.0, whiteref = 0.0;

    /* escalation ladder: threshold x scale; small-in-frame patterns
     * (cells ~10 px) need the 2x upsample to re-enter the validated
     * ~20 px regime; big frames (full HD) get two half-resolution
     * rungs FIRST -- the common in-range decode then costs a quarter
     * of the detection work, and corners are re-localized at full
     * resolution afterwards so no precision is lost */
    enum { WM_FULL, WM_UP2, WM_DOWN2 };
    double relthr[6];
    int wmode[6], npass = 0;
    unsigned char *up = NULL, *dn = NULL;
    const unsigned char *im;
    int iw, ih, wm = WM_FULL;
    int pass, decoded = 0;
    int best_m = -1, best_u = 0, best_v = 0;

    memset(res, 0, sizeof(*res));
    if (used_down2)
        *used_down2 = 0;

    if (allow_down2 && (long)w * h >= BIGFRAME_PX) {
        wmode[npass] = WM_DOWN2; relthr[npass++] = 0.55;
        wmode[npass] = WM_DOWN2; relthr[npass++] = 0.30;
    }
    wmode[npass] = WM_FULL; relthr[npass++] = 0.55;
    wmode[npass] = WM_FULL; relthr[npass++] = 0.30;
    wmode[npass] = WM_UP2;  relthr[npass++] = 0.55;
    wmode[npass] = WM_UP2;  relthr[npass++] = 0.30;

    /* Detection-threshold escalation: the strict threshold suppresses
     * dot L-corners when the pattern is large in frame (the sim regime);
     * when the pattern is small, blur can drop alternate saddle
     * polarities below it, leaving only the 45-degree sublattice -- so
     * on decode failure retry with the permissive threshold (tiny dots
     * cannot produce significant L-corner clutter). */
    for (pass = 0; pass < npass && !decoded; pass++) {
    memset(bit, -1, sizeof(bit));

    wm = wmode[pass];
    if (wm == WM_FULL) {
        im = img;
        iw = w;
        ih = h;
    } else if (wm == WM_DOWN2) {
        int dw = w / 2, dh = h / 2;
        if (!dn) {
            int x2, y2;
            dn = (unsigned char *)malloc((size_t)dw * dh);
            if (!dn)
                break;
            /* 2x2 box mean, truncating, matching the coarse tier's
             * downsample convention (full-res center of downsampled
             * pixel x is 2x + 0.5) */
            for (y2 = 0; y2 < dh; y2++)
                for (x2 = 0; x2 < dw; x2++) {
                    const unsigned char *q = img + (size_t)2 * y2 * w
                                             + 2 * x2;
                    dn[(size_t)y2 * dw + x2] = (unsigned char)
                        ((q[0] + q[1] + q[w] + q[w + 1]) / 4);
                }
        }
        im = dn;
        iw = dw;
        ih = dh;
    } else {
        if (!up) {
            int x2, y2;
            up = (unsigned char *)malloc((size_t)w * h * 4);
            if (!up)
                break;
            for (y2 = 0; y2 < 2 * h; y2++)
                for (x2 = 0; x2 < 2 * w; x2++) {
                    double sx = x2 * 0.5, sy = y2 * 0.5;
                    int ix = (int)sx, iy2 = (int)sy;
                    double fx = sx - ix, fy = sy - iy2;
                    const unsigned char *q;
                    if (ix >= w - 1) { ix = w - 2; fx = 1.0; }
                    if (iy2 >= h - 1) { iy2 = h - 2; fy = 1.0; }
                    q = img + iy2 * w + ix;
                    up[y2 * 2 * w + x2] = (unsigned char)
                        ((1 - fx) * (1 - fy) * q[0] + fx * (1 - fy) * q[1]
                         + (1 - fx) * fy * q[w] + fx * fy * q[w + 1]
                         + 0.5);
                }
        }
        im = up;
        iw = 2 * w;
        ih = 2 * h;
    }
    ncand = detect_saddles(cd, im, iw, ih, relthr[pass]);
    RDBG("reader: pass %d (mode %d, %dx%d) candidates %d\n", pass, wm,
         iw, ih, ncand);
    if (ncand < 12)
        continue;
    ngrid = grow_grid(&L, cd, ncand);
    RDBG("reader: lattice corners %d\n", ngrid);
    if (ngrid < 12)
        continue;

    /* --- decode cell bits (lattice square (i,j): corners i..i+1,j..j+1)
     * Two passes: measure all cells first, then set bits with a threshold
     * adaptive to the MEASURED contrast -- real cameras attenuate the
     * small dots (blur, compression) far below the ideal 255 swing. */
    {
        static double bases[GRIDN * GRIDN], dots[GRIDN * GRIDN];
        static int ci_[GRIDN * GRIDN], cj_[GRIDN * GRIDN];
        int nb = 0, k;
        double bmin = 255.0, bmax = 0.0, thr;
        for (j = 0; j < GRIDN - 1; j++)
            for (i = 0; i < GRIDN - 1; i++) {
                double cu, cv, bs[4], b, dot;
                if (!L.have[i][j] || !L.have[i + 1][j] || !L.have[i][j + 1]
                    || !L.have[i + 1][j + 1])
                    continue;
                /* the cell's PROJECTIVE center is the intersection of
                 * the quad's diagonals (the corner centroid is biased
                 * under strong perspective, enough to miss the dot) */
                {
                    const double *a = L.p[i][j], *b2 = L.p[i + 1][j + 1];
                    const double *c2 = L.p[i + 1][j], *e = L.p[i][j + 1];
                    double dx1 = b2[0] - a[0], dy1 = b2[1] - a[1];
                    double dx2 = e[0] - c2[0], dy2 = e[1] - c2[1];
                    double den = dx1 * dy2 - dy1 * dx2, t2;
                    if (fabs(den) < 1e-9)
                        continue;
                    t2 = ((c2[0] - a[0]) * dy2 - (c2[1] - a[1]) * dx2)
                         / den;
                    cu = a[0] + t2 * dx1;
                    cv = a[1] + t2 * dy1;
                }
                for (k = 0; k < 4; k++) {
                    const double *q =
                        L.p[i + (k & 1)][j + (k >> 1)];
                    bs[k] = bilin(im, iw, ih, cu + 0.32 * (q[0] - cu),
                                  cv + 0.32 * (q[1] - cv));
                }
                b = med4(bs[0], bs[1], bs[2], bs[3]);
                dot = sample3(im, iw, ih, cu, cv);
                if (b < 0.0 || dot < 0.0)
                    continue;
                base[i][j] = b;
                bases[nb] = b;
                dots[nb] = fabs(dot - b);
                ci_[nb] = i;
                cj_[nb] = j;
                nb++;
                if (b < bmin)
                    bmin = b;
                if (b > bmax)
                    bmax = b;
            }
        RDBG("reader: cells %d, base range %.0f..%.0f\n", nb, bmin, bmax);
        if (nb < 16 || bmax - bmin < 40.0)
            continue;
        /* per-cell normalization by LOCAL contrast (neighbor cells are
         * the opposite color, so |base - neighbor base| measures the
         * local black-white range under shading/vignetting), then a
         * 1-D 2-means split classifies dotted vs plain cells */
        {
            static double nrm[GRIDN * GRIDN];
            double m0, m1;
            int it;
            for (k = 0; k < nb; k++) {
                int i2 = ci_[k], j2 = cj_[k], nn = 0;
                double lc = 0.0;
                static const int d4[4][2] = { {1,0},{-1,0},{0,1},{0,-1} };
                int q;
                for (q = 0; q < 4; q++) {
                    int a2 = i2 + d4[q][0], b2 = j2 + d4[q][1];
                    if (a2 < 0 || b2 < 0 || a2 >= GRIDN || b2 >= GRIDN)
                        continue;
                    if (bit[a2][b2] >= -1 && base[a2][b2] > 0.0
                        && (L.have[a2][b2] || 1)) {
                        double df = fabs(base[i2][j2] - base[a2][b2]);
                        if (df > 5.0) {
                            lc += df;
                            nn++;
                        }
                    }
                }
                lc = nn ? lc / nn : 0.5 * (bmax - bmin);
                if (lc < 12.0)
                    lc = 12.0;
                nrm[k] = dots[k] / lc;
            }
            /* 2-means on nrm */
            m0 = 0.05;
            m1 = 0.6;
            for (it = 0; it < 12; it++) {
                double s0 = 0, s1 = 0;
                int n0 = 0, n1 = 0;
                for (k = 0; k < nb; k++) {
                    if (fabs(nrm[k] - m0) < fabs(nrm[k] - m1)) {
                        s0 += nrm[k];
                        n0++;
                    } else {
                        s1 += nrm[k];
                        n1++;
                    }
                }
                if (n0)
                    m0 = s0 / n0;
                if (n1)
                    m1 = s1 / n1;
            }
            RDBG("reader: dotness clusters %.2f / %.2f\n", m0, m1);
            if (m1 - m0 < 0.12)
                continue; /* modes not separable */
            thr = 0.5 * (m0 + m1);
            for (k = 0; k < nb; k++)
                bit[ci_[k]][cj_[k]] = (nrm[k] > thr) ? 1 : 0;
        }
        /* --- checkerboard-coherence pruning: keep only the largest
         * connected region whose cells alternate black/white with real
         * contrast; grown-in clutter (fabric, background) fails the
         * alternation test and is dropped before voting --- */
        {
            static signed char cls[GRIDN][GRIDN];
            static int comp[GRIDN][GRIDN];
            static int qx[GRIDN * GRIDN], qy[GRIDN * GRIDN];
            static const int d4[4][2] = { {1,0},{-1,0},{0,1},{0,-1} };
            int ncomp = 0, bestc = -1, bestn = 0, kept = 0;
            int i2, j2, q;
            memset(cls, -1, sizeof(cls));
            memset(comp, 0, sizeof(comp));
            for (k = 0; k < nb; k++) {
                double nm = 0.0;
                int nn = 0;
                i2 = ci_[k];
                j2 = cj_[k];
                for (q = 0; q < 4; q++) {
                    int a2 = i2 + d4[q][0], b2 = j2 + d4[q][1];
                    if (a2 >= 0 && b2 >= 0 && a2 < GRIDN && b2 < GRIDN
                        && bit[a2][b2] >= 0) {
                        nm += base[a2][b2];
                        nn++;
                    }
                }
                if (nn >= 2 && fabs(base[i2][j2] - nm / nn) > 15.0)
                    cls[i2][j2] = base[i2][j2] > nm / nn;
            }
            for (j2 = 0; j2 < GRIDN; j2++)
                for (i2 = 0; i2 < GRIDN; i2++) {
                    int qh2 = 0, qt2 = 0, sz = 0;
                    if (cls[i2][j2] < 0 || comp[i2][j2])
                        continue;
                    ncomp++;
                    comp[i2][j2] = ncomp;
                    qx[qt2] = i2;
                    qy[qt2] = j2;
                    qt2++;
                    while (qh2 < qt2) {
                        int cx2 = qx[qh2], cy2 = qy[qh2];
                        qh2++;
                        sz++;
                        for (q = 0; q < 4; q++) {
                            int a2 = cx2 + d4[q][0], b2 = cy2 + d4[q][1];
                            if (a2 < 0 || b2 < 0 || a2 >= GRIDN
                                || b2 >= GRIDN)
                                continue;
                            if (cls[a2][b2] < 0 || comp[a2][b2])
                                continue;
                            if (cls[a2][b2] == cls[cx2][cy2])
                                continue; /* must alternate */
                            comp[a2][b2] = ncomp;
                            qx[qt2] = a2;
                            qy[qt2] = b2;
                            qt2++;
                        }
                    }
                    if (sz > bestn) {
                        bestn = sz;
                        bestc = ncomp;
                    }
                }
            for (j2 = 0; j2 < GRIDN; j2++)
                for (i2 = 0; i2 < GRIDN; i2++)
                    if (bit[i2][j2] >= 0
                        && (cls[i2][j2] < 0 || comp[i2][j2] != bestc))
                        bit[i2][j2] = -1;
                    else if (bit[i2][j2] >= 0)
                        kept++;
            RDBG("reader: coherence kept %d cells (largest of %d comps)\n",
                 kept, ncomp);
            if (kept < 20)
                continue;
            if (getenv("MV_READER_DEBUG")) {
                int i2, j2, jmin = GRIDN, jmax = 0, imin = GRIDN, imax = 0;
                for (k = 0; k < nb; k++) {
                    if (ci_[k] < imin) imin = ci_[k];
                    if (ci_[k] > imax) imax = ci_[k];
                    if (cj_[k] < jmin) jmin = cj_[k];
                    if (cj_[k] > jmax) jmax = cj_[k];
                }
                fprintf(stderr, "reader: bit matrix (lattice %d..%d x %d..%d):\n",
                        imin, imax, jmin, jmax);
                for (j2 = jmin; j2 <= jmax; j2++) {
                    fputs("  ", stderr);
                    for (i2 = imin; i2 <= imax; i2++)
                        fputc(bit[i2][j2] < 0 ? '.' : '0' + bit[i2][j2],
                              stderr);
                    fputc('\n', stderr);
                }
            }
        }
        if (getenv("MV_READER_DEBUG")) {
            double s[GRIDN * GRIDN];
            int a2, b2;
            memcpy(s, dots, (size_t)nb * sizeof(double));
            for (a2 = 1; a2 < nb; a2++)
                for (b2 = a2; b2 > 0 && s[b2] < s[b2 - 1]; b2--) {
                    double t2 = s[b2];
                    s[b2] = s[b2 - 1];
                    s[b2 - 1] = t2;
                }
            RDBG("reader: |dot-base| quartiles %.0f/%.0f/%.0f/%.0f/%.0f,"
                 " thr %.0f\n", s[0], s[nb / 4], s[nb / 2],
                 s[3 * nb / 4], s[nb - 1], thr);
        }
        /* refs for counter decoding */
        {
            double mid = 0.5 * (bmin + bmax), sb = 0, sw = 0;
            int cb = 0, cw = 0, k;
            for (k = 0; k < nb; k++) {
                if (bases[k] < mid) {
                    sb += bases[k];
                    cb++;
                } else {
                    sw += bases[k];
                    cw++;
                }
            }
            if (cb)
                blackref = sb / cb;
            if (cw)
                whiteref = sw / cw;
        }
    }

    /* --- window votes for (rot m, offsets u,v) --- */
    {
        int votes_m[4] = { 0, 0, 0, 0 };
        int best_n = 0;
        int vu[256], vv[256], vm[256], vc[256], nv = 0, k;
        best_m = -1;
        for (j = 0; j + MV_PAT_WIN - 1 < GRIDN - 1; j++)
            for (i = 0; i + MV_PAT_WIN - 1 < GRIDN - 1; i++) {
                unsigned code = 0;
                int r, c, ok = 1, pc, pr, m, u = 0, v = 0;
                for (r = 0; r < MV_PAT_WIN && ok; r++)
                    for (c = 0; c < MV_PAT_WIN && ok; c++) {
                        if (bit[i + c][j + r] < 0)
                            ok = 0;
                        else
                            code = (code << 1)
                                 | (unsigned)bit[i + c][j + r];
                    }
                if (!ok)
                    continue;
                if (mv_marray_lookup(code, &pc, &pr, &m) != MV_OK)
                    continue;
                switch (m) {
                case 0: u = pc - i; v = pr - j; break;
                case 1: u = pc + j + 3; v = pr - i; break;
                case 2: u = pc + i + 3; v = pr + j + 3; break;
                default: u = pc - j; v = pr + i + 3; break;
                }
                votes_m[m]++;
                for (k = 0; k < nv; k++)
                    if (vm[k] == m && vu[k] == u && vv[k] == v) {
                        vc[k]++;
                        break;
                    }
                if (k == nv && nv < 256) {
                    vm[nv] = m;
                    vu[nv] = u;
                    vv[nv] = v;
                    vc[nv] = 1;
                    nv++;
                }
            }
        for (k = 0; k < nv; k++)
            if (vc[k] > best_n) {
                best_n = vc[k];
                best_m = vm[k];
                best_u = vu[k];
                best_v = vv[k];
            }
        RDBG("reader: votes m=[%d %d %d %d], best (m=%d,u=%d,v=%d) x%d of %d\n",
             votes_m[0], votes_m[1], votes_m[2], votes_m[3],
             best_m, best_u, best_v, best_n, nv);
        if (best_m < 0 || best_n < 3)
            continue;
        res->rot = best_m;
        decoded = 1;
    }
    } /* escalation loop */
    free(up);
    free(dn);
    up = NULL;
    dn = NULL;
    if (!decoded)
        return MV_ERR;
    if (used_down2 && wm == WM_DOWN2)
        *used_down2 = 1;
    {

        /* --- corners: map each lattice corner via its 4 squares --- */
        for (j = 1; j < GRIDN - 1; j++)
            for (i = 1; i < GRIDN - 1; i++) {
                int a, b, pi_min = 999, pj_min = 999, ok = 1;
                if (!L.have[i][j])
                    continue;
                for (b = -1; b <= 0 && ok; b++)
                    for (a = -1; a <= 0 && ok; a++) {
                        int si = i + a, sj = j + b, pi, pj;
                        switch (best_m) {
                        case 0: pi = si + best_u; pj = sj + best_v; break;
                        case 1: pi = best_u - sj; pj = best_v + si; break;
                        case 2: pi = best_u - si; pj = best_v - sj; break;
                        default: pi = best_u + sj; pj = best_v - si; break;
                        }
                        if (pi < 0 || pi >= MV_PAT_GRID_COLS || pj < 0
                            || pj >= MV_PAT_GRID_ROWS)
                            ok = 0;
                        if (pi < pi_min)
                            pi_min = pi;
                        if (pj < pj_min)
                            pj_min = pj;
                    }
                if (!ok)
                    continue;
                if (pi_min < 0 || pi_min >= MV_PAT_CORNER_COLS
                    || pj_min < 0 || pj_min >= MV_PAT_CORNER_ROWS)
                    continue;
                if (res->n < MV_READ_MAXC) {
                    double uu = L.p[i][j][0], vv = L.p[i][j][1];
                    if (wm == WM_UP2) {
                        uu *= 0.5;
                        vv *= 0.5;
                    } else if (wm == WM_DOWN2) {
                        uu = 2.0 * uu + 0.5;
                        vv = 2.0 * vv + 0.5;
                        refine_saddle(img, w, h, &uu, &vv);
                    }
                    res->id[res->n] = pj_min * MV_PAT_CORNER_COLS + pi_min;
                    res->uv[2 * res->n] = uu;
                    res->uv[2 * res->n + 1] = vv;
                    res->n++;
                }
            }
    }
    RDBG("reader: identified corners %d (rot %d)\n", res->n, res->rot);
    if (res->n < 8)
        return MV_ERR;

    /* --- homography pattern px -> image --- */
    {
        static double obj[2 * MV_READ_MAXC];
        for (i = 0; i < res->n; i++) {
            double xy[2];
            mv_pattern_corner_px(res->id[i] % MV_PAT_CORNER_COLS,
                                 res->id[i] / MV_PAT_CORNER_COLS, xy);
            obj[2 * i] = xy[0];
            obj[2 * i + 1] = xy[1];
        }
        if (mv_homography_dlt(res->H, obj, res->uv, res->n) != MV_OK)
            return MV_ERR;
        /* residual gate: lattice growth can capture auxiliary screen
         * elements (version block, nested reference grid) whose
         * lattice positions map into valid corner ids; a corner far
         * from the consensus homography is such a capture error, not
         * noise.  Two rounds: a wide gate first so catastrophic
         * outliers stop skewing the fit, then a tight one. */
        {
            static const double GATE[2] = { 20.0, 3.0 };
            int round, k;
            for (round = 0; round < 2; round++) {
                int nk = 0;
                for (k = 0; k < res->n; k++) {
                    double pu, pv, du, dv;
                    proj_H(res->H, obj[2 * k], obj[2 * k + 1], &pu,
                           &pv);
                    du = res->uv[2 * k] - pu;
                    dv = res->uv[2 * k + 1] - pv;
                    if (du * du + dv * dv
                        <= GATE[round] * GATE[round]) {
                        res->id[nk] = res->id[k];
                        res->uv[2 * nk] = res->uv[2 * k];
                        res->uv[2 * nk + 1] = res->uv[2 * k + 1];
                        obj[2 * nk] = obj[2 * k];
                        obj[2 * nk + 1] = obj[2 * k + 1];
                        nk++;
                    }
                }
                if (nk == res->n)
                    continue;
                RDBG("reader: residual gate %.0f px dropped %d of %d\n",
                     GATE[round], res->n - nk, res->n);
                res->n = nk;
                if (res->n < 8)
                    return MV_ERR;
                if (mv_homography_dlt(res->H, obj, res->uv, res->n)
                    != MV_OK)
                    return MV_ERR;
            }
        }
        /* sub-pixel polish at full resolution: when the cells are big
         * enough that a corner's polish window sees only its own two
         * edges (>= 14 px keeps the M-array dots clear), refine every
         * corner by edge-line intersection (gradient-orthogonality
         * fallback) and refit.  Cuts corner error from the detector's
         * ~0.35 px quadratic-peak floor to well under 0.1 px;
         * tiny-cell (upsampled) decodes keep the detector estimate. */
        if (cell_size_px(res->H,
                         MV_PAT_GRID_X0
                         + 0.5 * MV_PAT_GRID_COLS * MV_PAT_CELL,
                         MV_PAT_GRID_Y0
                         + 0.5 * MV_PAT_GRID_ROWS * MV_PAT_CELL,
                         MV_PAT_CELL) >= 14.0) {
            int moved = 0;
            for (i = 0; i < res->n; i++) {
                if (refine_corner_edges(img, w, h, res->H, obj[2 * i],
                                        obj[2 * i + 1], MV_PAT_CELL,
                                        &res->uv[2 * i],
                                        &res->uv[2 * i + 1]) == MV_OK)
                    moved++;
                else
                    moved += foerstner_step(img, w, h, &res->uv[2 * i],
                                            &res->uv[2 * i + 1]);
            }
            if (moved
                && mv_homography_dlt(res->H, obj, res->uv, res->n)
                   != MV_OK)
                return MV_ERR;
            /* corner recovery: ids the lattice never reached (blur or
             * noise broke a saddle at working resolution) now have a
             * predicted position under H; measure each candidate with
             * the edge-line intersection ALONE (no gradient fallback:
             * recovery demands the strong evidence) and accept only if
             * it lands within 1 px of the prediction.  Occluded and
             * out-of-frame corners fail the crossing checks and stay
             * absent. */
            {
                static char have[MV_PAT_CORNER_COLS
                                 * MV_PAT_CORNER_ROWS];
                int k, added = 0;
                memset(have, 0, sizeof(have));
                for (k = 0; k < res->n; k++)
                    have[res->id[k]] = 1;
                for (k = 0;
                     k < MV_PAT_CORNER_COLS * MV_PAT_CORNER_ROWS
                     && res->n < MV_READ_MAXC; k++) {
                    double xy[2], pu, pv, ru, rv, du, dv;
                    if (have[k])
                        continue;
                    mv_pattern_corner_px(k % MV_PAT_CORNER_COLS,
                                         k / MV_PAT_CORNER_COLS, xy);
                    proj_H(res->H, xy[0], xy[1], &pu, &pv);
                    if (pu < 8.0 || pv < 8.0 || pu >= w - 8.0
                        || pv >= h - 8.0)
                        continue;
                    ru = pu;
                    rv = pv;
                    if (refine_corner_edges(img, w, h, res->H, xy[0],
                                            xy[1], MV_PAT_CELL, &ru,
                                            &rv) != MV_OK)
                        continue;
                    du = ru - pu;
                    dv = rv - pv;
                    if (du * du + dv * dv > 1.0)
                        continue;
                    res->id[res->n] = k;
                    res->uv[2 * res->n] = ru;
                    res->uv[2 * res->n + 1] = rv;
                    obj[2 * res->n] = xy[0];
                    obj[2 * res->n + 1] = xy[1];
                    res->n++;
                    added++;
                }
                if (added) {
                    RDBG("reader: recovered %d corners via H\n", added);
                    if (mv_homography_dlt(res->H, obj, res->uv, res->n)
                        != MV_OK)
                        return MV_ERR;
                }
            }
        }
    }

    /* --- counter strip --- */
    {
        unsigned bits[MV_PAT_CTR_CELLS], g = 0, parity = 0;
        double mid = 0.5 * (blackref + whiteref);
        double halfrange = 0.5 * (whiteref - blackref);
        double conf = 1.0;
        int ok = 1, c;
        for (c = 0; c < MV_PAT_CTR_CELLS && ok; c++) {
            double xy[2], W, u, v, s, m;
            mv_pattern_ctr_cell_px(c, xy);
            W = res->H[6] * xy[0] + res->H[7] * xy[1] + res->H[8];
            if (fabs(W) < 1e-12) {
                ok = 0;
                break;
            }
            u = (res->H[0] * xy[0] + res->H[1] * xy[1] + res->H[2]) / W;
            v = (res->H[3] * xy[0] + res->H[4] * xy[1] + res->H[5]) / W;
            s = sample3(img, w, h, u, v);
            if (s < 0.0) {
                ok = 0;
                break;
            }
            bits[c] = (s < mid) ? 1u : 0u; /* black = bit 1 */
            m = fabs(s - mid) / (halfrange > 1.0 ? halfrange : 1.0);
            if (m > 1.0)
                m = 1.0;
            if (m < conf)
                conf = m;
        }
        if (ok && bits[0] == 1u && bits[1] == 0u) {
            for (c = 0; c < MV_PAT_CTR_BITS; c++) {
                g = (g << 1) | bits[2 + c];
                parity ^= bits[2 + c];
            }
            if (bits[22] == parity && bits[23] == (parity ^ 1u)) {
                unsigned b = g;
                b ^= b >> 1;
                b ^= b >> 2;
                b ^= b >> 4;
                b ^= b >> 8;
                b ^= b >> 16;
                res->counter = b;
                res->counter_valid = 1;
                res->counter_conf = conf;
            }
        }
    }
    return MV_OK;
}

int mv_read_pattern(mv_read_result *res, const unsigned char *img,
                    int w, int h)
{
    int used_down2 = 0;
    int r = read_pattern_once(res, img, w, h, 1, &used_down2);
    /* A half-resolution decode that identified well under the corners
     * its own homography says are in view has lost real corners to
     * the downscale (strong tilt foreshortens the far cells below the
     * detector's regime); the pattern is demonstrably present, so one
     * native-resolution retry is cheap relative to a lost anchor, and
     * the richer result wins.  Partial visibility does NOT trigger
     * this: out-of-frame corners are not counted as expected. */
    if (r == MV_OK && used_down2) {
        int vis = 0, i, j;
        for (j = 0; j < MV_PAT_CORNER_ROWS; j++)
            for (i = 0; i < MV_PAT_CORNER_COLS; i++) {
                double xy[2], u, v;
                mv_pattern_corner_px(i, j, xy);
                proj_H(res->H, xy[0], xy[1], &u, &v);
                if (u >= 8.0 && v >= 8.0 && u < w - 8.0 && v < h - 8.0)
                    vis++;
            }
        if (10 * res->n < 9 * vis) {
            mv_read_result full;
            RDBG("reader: down2 kept %d of %d visible, native retry\n",
                 res->n, vis);
            if (read_pattern_once(&full, img, w, h, 0, NULL) == MV_OK
                && full.n > res->n)
                *res = full;
        }
    }
    return r;
}

/* ================= coarse tier (spec v2) ============================= */

/* Map an observed lattice cell (a,b) within a di x dj bounding box to
 * pattern corner (i 0..4, j 0..1) under rot clockwise quarter-turns.
 * Landscape observations (di=5) pair with rot 0/2, portrait (di=2)
 * with rot 1/3. */
static void coarse_map(int rot, int a, int b, int *i, int *j)
{
    switch (rot) {
    case 0:  *i = a;     *j = b;     break;
    case 1:  *i = b;     *j = 1 - a; break;
    case 2:  *i = 4 - a; *j = 1 - b; break;
    default: *i = 4 - b; *j = a;     break;
    }
}

/* Coarse counter strip: [1,0] sync + 8 Gray bits + parity + ~parity,
 * sampled through H on the given image against the checker-derived
 * mid level.  Sets counter/counter_valid/counter_conf on a clean
 * decode; leaves res untouched otherwise.  Single implementation for
 * the in-validation decode and the full-resolution retry. */
static void coarse_strip_decode(mv_read_result *res,
                                const unsigned char *img, int w, int h,
                                const double H[9], double mid,
                                double halfrange)
{
    unsigned bits[MV_PAT2_CTR_CELLS], g2 = 0, parity = 0;
    double conf = 1.0;
    int okc = 1, c;

    for (c = 0; c < MV_PAT2_CTR_CELLS && okc; c++) {
        double X = MV_PAT2_CTR_X0 + (c + 0.5) * MV_PAT2_CTR_CELL;
        double Y = MV_PAT2_CTR_Y0 + 0.5 * MV_PAT2_CTR_H;
        double W = H[6] * X + H[7] * Y + H[8], u, v, s, m;
        if (fabs(W) < 1e-12) {
            okc = 0;
            break;
        }
        u = (H[0] * X + H[1] * Y + H[2]) / W;
        v = (H[3] * X + H[4] * Y + H[5]) / W;
        s = sample3(img, w, h, u, v);
        if (s < 0.0) {
            okc = 0;
            break;
        }
        bits[c] = (s < mid) ? 1u : 0u;
        m = fabs(s - mid) / (halfrange > 1.0 ? halfrange : 1.0);
        if (m > 1.0)
            m = 1.0;
        if (m < conf)
            conf = m;
    }
    if (okc)
        RDBG("coarse strip: bits %u%u %u%u%u%u%u%u%u%u p %u%u "
             "(mid %.0f)\n", bits[0], bits[1], bits[2], bits[3],
             bits[4], bits[5], bits[6], bits[7], bits[8], bits[9],
             bits[10], bits[11], mid);
    if (okc && bits[0] == 1u && bits[1] == 0u) {
        for (c = 0; c < MV_PAT2_CTR_BITS; c++) {
            g2 = (g2 << 1) | bits[2 + c];
            parity ^= bits[2 + c];
        }
        if (bits[10] == parity && bits[11] == (parity ^ 1u)) {
            unsigned bv = g2;
            bv ^= bv >> 1;
            bv ^= bv >> 2;
            bv ^= bv >> 4;
            res->counter = bv;
            res->counter_valid = 1;
            res->counter_conf = conf;
        }
    }
}

/* Validate one (lattice, rot) hypothesis on the working image: fit the
 * homography, require the 6x3 checker phase to match exactly, require
 * the observed orientation-mark set to equal the spec's asymmetric L,
 * then decode the counter strip.  Fills res (uv in working-image
 * coordinates) on success. */
static int coarse_validate(mv_read_result *res, const unsigned char *img,
                           int w, int h, struct lattice *L,
                           int imin, int jmin, int di, int dj, int rot)
{
    double obj[2 * 10], uvs[2 * 10], H[9];
    int ids[10];
    int n = 0, a, b, k;

    for (b = 0; b < dj; b++)
        for (a = 0; a < di; a++) {
            int pi, pj;
            double xy[2];
            if (!L->have[imin + a][jmin + b])
                continue;
            coarse_map(rot, a, b, &pi, &pj);
            if (pi < 0 || pi > 4 || pj < 0 || pj > 1)
                return MV_ERR;
            mv_pattern2_corner_px(pi, pj, xy);
            obj[2 * n] = xy[0];
            obj[2 * n + 1] = xy[1];
            uvs[2 * n] = L->p[imin + a][jmin + b][0];
            uvs[2 * n + 1] = L->p[imin + a][jmin + b][1];
            ids[n] = pj * MV_PAT2_CORNER_COLS + pi;
            n++;
        }
    if (n < 6)
        return MV_ERR;
    if (mv_homography_dlt(H, obj, uvs, n) != MV_OK)
        return MV_ERR;

    {
        double base[MV_PAT2_GRID_ROWS][MV_PAT2_GRID_COLS];
        double cen[MV_PAT2_GRID_ROWS][MV_PAT2_GRID_COLS];
        double mb = 0.0, mw = 0.0, mid;
        int nb = 0, nw = 0, r, c;
        static const double off[4][2] = {
            { -0.30, -0.30 }, { 0.30, -0.30 },
            { -0.30, 0.30 }, { 0.30, 0.30 }
        };
        for (r = 0; r < MV_PAT2_GRID_ROWS; r++)
            for (c = 0; c < MV_PAT2_GRID_COLS; c++) {
                double cx = MV_PAT2_GRID_X0 + (c + 0.5) * MV_PAT2_CELL;
                double cy = MV_PAT2_GRID_Y0 + (r + 0.5) * MV_PAT2_CELL;
                double s = 0.0;
                for (k = 0; k < 4; k++) {
                    double X = cx + off[k][0] * MV_PAT2_CELL;
                    double Y = cy + off[k][1] * MV_PAT2_CELL;
                    double W = H[6] * X + H[7] * Y + H[8], u, v, g;
                    if (fabs(W) < 1e-12)
                        return MV_ERR;
                    u = (H[0] * X + H[1] * Y + H[2]) / W;
                    v = (H[3] * X + H[4] * Y + H[5]) / W;
                    g = sample3(img, w, h, u, v);
                    if (g < 0.0)
                        return MV_ERR;
                    s += g;
                }
                base[r][c] = s / 4.0;
                {
                    double W = H[6] * cx + H[7] * cy + H[8], u, v, g;
                    if (fabs(W) < 1e-12)
                        return MV_ERR;
                    u = (H[0] * cx + H[1] * cy + H[2]) / W;
                    v = (H[3] * cx + H[4] * cy + H[5]) / W;
                    g = sample3(img, w, h, u, v);
                    if (g < 0.0)
                        return MV_ERR;
                    cen[r][c] = g;
                }
                if ((r + c) % 2 == 0) {
                    mb += base[r][c];
                    nb++;
                } else {
                    mw += base[r][c];
                    nw++;
                }
            }
        mb /= nb;
        mw /= nw;
        if (mw - mb < 25.0)
            return MV_ERR;
        mid = 0.5 * (mb + mw);
        for (r = 0; r < MV_PAT2_GRID_ROWS; r++)
            for (c = 0; c < MV_PAT2_GRID_COLS; c++) {
                int expblack = ((r + c) % 2 == 0);
                int marked;
                if ((base[r][c] < mid) != expblack)
                    return MV_ERR;
                marked = fabs(cen[r][c] - base[r][c]) > 0.5 * (mw - mb);
                if (marked != mv_pattern2_mark(c, r))
                    return MV_ERR;
            }

        /* counter strip: [1,0] + 8 Gray + parity + ~parity */
        res->counter_valid = 0;
        coarse_strip_decode(res, img, w, h, H, mid, 0.5 * (mw - mb));
    }

    res->n = n;
    for (k = 0; k < n; k++) {
        res->id[k] = ids[k];
        res->uv[2 * k] = uvs[2 * k];
        res->uv[2 * k + 1] = uvs[2 * k + 1];
    }
    memcpy(res->H, H, sizeof(H));
    res->rot = rot;
    return MV_OK;
}

static int coarse_try(mv_read_result *res, const unsigned char *img,
                      int w, int h, double relthr)
{
    static struct cand cd[MAXCAND];
    static char tried[MAXCAND];
    static struct lattice L;
    int n, attempt;

    n = detect_saddles(cd, img, w, h, relthr);
    RDBG("coarse: %d saddle candidates (thr %.2f)\n", n, relthr);
    if (n < 8)
        return MV_ERR;
    memset(tried, 0, sizeof(tried));
    for (attempt = 0; attempt < 24; attempt++) {
        int seed = -1, k, count, i, j, di, dj, rot;
        int imin = GRIDN, imax = -1, jmin = GRIDN, jmax = -1;
        for (k = 0; k < n; k++)
            if (!tried[k] && (seed < 0 || cd[k].resp > cd[seed].resp))
                seed = k;
        if (seed < 0)
            break;
        tried[seed] = 1;
        count = grow_from_seed(&L, cd, n, seed);
        if (count < 8 || count > 40)
            continue;
        for (j = 0; j < GRIDN; j++)
            for (i = 0; i < GRIDN; i++)
                if (L.have[i][j]) {
                    if (i < imin)
                        imin = i;
                    if (i > imax)
                        imax = i;
                    if (j < jmin)
                        jmin = j;
                    if (j > jmax)
                        jmax = j;
                }
        di = imax - imin + 1;
        dj = jmax - jmin + 1;
        RDBG("coarse: seed %d grew %d corners, bbox %dx%d\n", seed,
             count, di, dj);
        if (di > 14 || dj > 14)
            continue;
        /* the grown lattice may include stray captures (strip gutters,
         * mark corners, foreshortened neighbors); do not demand a clean
         * bounding box -- enumerate candidate 5x2 / 2x5 windows inside
         * it and let the exact validator (checker phase + mark L +
         * counter) decide */
        {
            int wi, wj, ww, wh, shape;
            for (shape = 0; shape < 2; shape++) {
                ww = shape ? 2 : 5;
                wh = shape ? 5 : 2;
                if (di < ww || dj < wh)
                    continue;
                for (wj = jmin; wj + wh - 1 <= jmax; wj++)
                    for (wi = imin; wi + ww - 1 <= imax; wi++) {
                        int a, b, present = 0;
                        for (b = 0; b < wh; b++)
                            for (a = 0; a < ww; a++)
                                present += L.have[wi + a][wj + b];
                        if (present < 8)
                            continue;
                        for (rot = 0; rot < 4; rot++) {
                            if ((rot % 2 == 0) != (ww == 5))
                                continue;
                            if (coarse_validate(res, img, w, h, &L, wi,
                                                wj, ww, wh, rot)
                                == MV_OK) {
                                RDBG("coarse: validated rot %d, %d "
                                     "corners, counter %svalid\n", rot,
                                     res->n,
                                     res->counter_valid ? "" : "in");
                                return MV_OK;
                            }
                        }
                    }
            }
        }
    }
    return MV_ERR;
}

/* Shared full-resolution polish for a successful coarse read: when the
 * cells are big enough for the crossing profiles to stay inside their
 * own cells, refine each corner by edge-line intersection (gradient
 * fallback) and refit the homography. */
static void coarse_polish(mv_read_result *res, const unsigned char *img,
                          int w, int h)
{
    static double obj[2 * 10];
    int i, moved = 0;

    if (cell_size_px(res->H,
                     MV_PAT2_GRID_X0
                     + 0.5 * MV_PAT2_GRID_COLS * MV_PAT2_CELL,
                     MV_PAT2_GRID_Y0
                     + 0.5 * MV_PAT2_GRID_ROWS * MV_PAT2_CELL,
                     MV_PAT2_CELL) < 14.0)
        return;
    for (i = 0; i < res->n; i++) {
        double xy[2];
        mv_pattern2_corner_px(res->id[i] % MV_PAT2_CORNER_COLS,
                              res->id[i] / MV_PAT2_CORNER_COLS, xy);
        obj[2 * i] = xy[0];
        obj[2 * i + 1] = xy[1];
        if (refine_corner_edges(img, w, h, res->H, xy[0], xy[1],
                                MV_PAT2_CELL, &res->uv[2 * i],
                                &res->uv[2 * i + 1]) == MV_OK)
            moved++;
        else
            moved += foerstner_step(img, w, h, &res->uv[2 * i],
                                    &res->uv[2 * i + 1]);
    }
    if (moved)
        mv_homography_dlt(res->H, obj, res->uv, res->n);
    /* corner recovery, as in the fine tier: measure missing ids at
     * their H-predicted positions with the edge-line intersection
     * alone, accept within 1 px */
    {
        char have[MV_PAT2_CORNER_COLS * MV_PAT2_CORNER_ROWS];
        int k, added = 0;
        memset(have, 0, sizeof(have));
        for (k = 0; k < res->n; k++)
            have[res->id[k]] = 1;
        for (k = 0; k < MV_PAT2_CORNER_COLS * MV_PAT2_CORNER_ROWS
                    && res->n < MV_READ_MAXC; k++) {
            double xy[2], pu, pv, ru, rv, du, dv;
            if (have[k])
                continue;
            mv_pattern2_corner_px(k % MV_PAT2_CORNER_COLS,
                                  k / MV_PAT2_CORNER_COLS, xy);
            proj_H(res->H, xy[0], xy[1], &pu, &pv);
            if (pu < 8.0 || pv < 8.0 || pu >= w - 8.0 || pv >= h - 8.0)
                continue;
            ru = pu;
            rv = pv;
            if (refine_corner_edges(img, w, h, res->H, xy[0], xy[1],
                                    MV_PAT2_CELL, &ru, &rv) != MV_OK)
                continue;
            du = ru - pu;
            dv = rv - pv;
            if (du * du + dv * dv > 1.0)
                continue;
            res->id[res->n] = k;
            res->uv[2 * res->n] = ru;
            res->uv[2 * res->n + 1] = rv;
            obj[2 * res->n] = xy[0];
            obj[2 * res->n + 1] = xy[1];
            res->n++;
            added++;
        }
        if (added) {
            RDBG("coarse: recovered %d corners via H\n", added);
            mv_homography_dlt(res->H, obj, res->uv, res->n);
        }
    }
}

/* Full-resolution counter retry for a decode won on a downsampled
 * working image: the 60 px-high strip can blur below the margin there
 * while remaining crisp at native resolution.  Re-derives the checker
 * brightness references through the final H and reruns the SAME strip
 * decode; only ever turns an invalid counter into a valid one. */
static void coarse_counter_fullres(mv_read_result *res,
                                   const unsigned char *img, int w,
                                   int h)
{
    static const double off[4][2] = {
        { -0.30, -0.30 }, { 0.30, -0.30 },
        { -0.30, 0.30 }, { 0.30, 0.30 }
    };
    double mb = 0.0, mw = 0.0;
    int nb = 0, nw = 0, r, c, k, ok = 1;

    for (r = 0; r < MV_PAT2_GRID_ROWS && ok; r++)
        for (c = 0; c < MV_PAT2_GRID_COLS && ok; c++) {
            double cx = MV_PAT2_GRID_X0 + (c + 0.5) * MV_PAT2_CELL;
            double cy = MV_PAT2_GRID_Y0 + (r + 0.5) * MV_PAT2_CELL;
            double s = 0.0;
            for (k = 0; k < 4 && ok; k++) {
                double u, v, g;
                proj_H(res->H, cx + off[k][0] * MV_PAT2_CELL,
                       cy + off[k][1] * MV_PAT2_CELL, &u, &v);
                g = sample3(img, w, h, u, v);
                if (g < 0.0)
                    ok = 0;
                else
                    s += g;
            }
            if (!ok)
                break;
            if ((r + c) % 2 == 0) {
                mb += s / 4.0;
                nb++;
            } else {
                mw += s / 4.0;
                nw++;
            }
        }
    if (!ok || !nb || !nw)
        return;
    mb /= nb;
    mw /= nw;
    if (mw - mb < 25.0)
        return;
    coarse_strip_decode(res, img, w, h, res->H, 0.5 * (mb + mw),
                        0.5 * (mw - mb));
}

int mv_read_coarse(mv_read_result *res, const unsigned char *img,
                   int w, int h)
{
    /* Big frames try the cheap downsampled rungs first (detection cost
     * scales with pixels; corners are re-localized at full resolution
     * on success, so precision does not follow the working scale).
     * Small frames keep the native-resolution-first order: their s=1
     * pass is cheap and the far-range regime lives there. */
    static const int SC_NEAR[3] = { 1, 2, 4 };
    static const int SC_BIG[3] = { 2, 4, 1 };
    static const double THR[2] = { 0.55, 0.30 };
    const int *SC = ((long)w * h >= BIGFRAME_PX) ? SC_BIG : SC_NEAR;
    mv_read_result keep;
    int have_keep = 0;
    int si, ti, i;

    memset(res, 0, sizeof(*res));
    for (si = 0; si < 3; si++) {
        int s = SC[si], ws = w / s, hs = h / s;
        int got = 0;
        unsigned char *ds;
        if (ws < 32 || hs < 32)
            continue;
        if (s == 1)
            ds = (unsigned char *)img;
        else {
            int x, y, a, b;
            ds = (unsigned char *)malloc((size_t)ws * hs);
            if (!ds)
                continue;
            for (y = 0; y < hs; y++)
                for (x = 0; x < ws; x++) {
                    int acc = 0;
                    for (b = 0; b < s; b++)
                        for (a = 0; a < s; a++)
                            acc += img[(y * s + b) * w + x * s + a];
                    ds[y * ws + x] = (unsigned char)(acc / (s * s));
                }
        }
        for (ti = 0; ti < 2 && !got; ti++)
            got = (coarse_try(res, ds, ws, hs, THR[ti]) == MV_OK);
        if (got && s != 1) {
            static double obj[2 * 10];
            for (i = 0; i < res->n; i++) {
                double xy[2];
                res->uv[2 * i] = res->uv[2 * i] * s + (s - 1) * 0.5;
                res->uv[2 * i + 1] = res->uv[2 * i + 1] * s
                                     + (s - 1) * 0.5;
                refine_saddle(img, w, h, &res->uv[2 * i],
                              &res->uv[2 * i + 1]);
                mv_pattern2_corner_px(res->id[i] % MV_PAT2_CORNER_COLS,
                                      res->id[i] / MV_PAT2_CORNER_COLS,
                                      xy);
                obj[2 * i] = xy[0];
                obj[2 * i + 1] = xy[1];
            }
            mv_homography_dlt(res->H, obj, res->uv, res->n);
        }
        if (s != 1)
            free(ds);
        if (!got)
            continue;
        coarse_polish(res, img, w, h);
        if (s != 1 && !res->counter_valid)
            coarse_counter_fullres(res, img, w, h);
        if (res->counter_valid)
            return MV_OK;
        /* Decode is good but the counter strip would not read at this
         * working scale (blur can erase the 60 px-high strip before it
         * touches the checker): keep the first such result and let the
         * remaining scales try -- a later rung may read the strip, and
         * an anchor without a counter cannot time-align. */
        if (!have_keep) {
            keep = *res;
            have_keep = 1;
        }
        RDBG("coarse: scale %d decoded but counter invalid, "
             "escalating\n", s);
    }
    if (have_keep) {
        *res = keep;
        return MV_OK;
    }
    return MV_ERR;
}

/* ================= display outline (see reader.h) ==================== */

static void proj_H(const double H[9], double X, double Y, double *u,
                   double *v)
{
    double W = H[6] * X + H[7] * Y + H[8];
    if (fabs(W) < 1e-12)
        W = 1e-12;
    *u = (H[0] * X + H[1] * Y + H[2]) / W;
    *v = (H[3] * X + H[4] * Y + H[5]) / W;
}

int mv_display_outline(double quad[8], const mv_read_result *rr,
                       const unsigned char *img, int w, int h)
{
    double whiteref = 0.0, blackref = 0.0, thr;
    int nw = 0, nb = 0;
    int fine = rr->n > 10;
    int gx0 = fine ? MV_PAT_GRID_X0 : MV_PAT2_GRID_X0;
    int gy0 = fine ? MV_PAT_GRID_Y0 : MV_PAT2_GRID_Y0;
    int cell = fine ? MV_PAT_CELL : MV_PAT2_CELL;
    int cols = fine ? MV_PAT_GRID_COLS : MV_PAT2_GRID_COLS;
    int rows = fine ? MV_PAT_GRID_ROWS : MV_PAT2_GRID_ROWS;
    int r, c, i, side;
    double lim[4]; /* panel bounds in pattern coords: xl, xr, yt, yb */

    /* brightness references from the pattern's own cells, sampled off
     * center to dodge marks and dots */
    for (r = 0; r < rows; r++)
        for (c = 0; c < cols; c++) {
            double u, v, g;
            proj_H(rr->H, gx0 + (c + 0.30) * cell, gy0 + (r + 0.30) * cell,
                   &u, &v);
            g = sample3(img, w, h, u, v);
            if (g < 0.0)
                continue;
            if ((r + c) % 2 == 0) {
                blackref += g;
                nb++;
            } else {
                whiteref += g;
                nw++;
            }
        }
    if (nw < 4 || nb < 4)
        return MV_ERR;
    whiteref /= nw;
    blackref /= nb;
    if (whiteref - blackref < 20.0)
        return MV_ERR;
    thr = blackref + 0.35 * (whiteref - blackref);

    /* march rays OUTWARD in pattern-plane coordinates, where the lit
     * panel's edges are axis-aligned straight lines by construction.
     * Each side takes the MEDIAN stop coordinate over its rays, which
     * is immune to leaks through any single bright bridge (a lit arm,
     * a bright wall touching the panel in the image) -- the failure
     * mode of region flooding. */
    for (side = 0; side < 4; side++) {
        double stops[64];
        int ns = 0;
        for (i = 0; i < 48; i++) {
            /* ray anchor spread along the pattern's extent */
            double along = (i + 0.5) / 48.0;
            double px, py, dx = 0.0, dy = 0.0, pos, g, u, v;
            int dark = 0, steps;
            switch (side) {
            case 0: /* left: start in left margin, march -x */
                px = gx0 * 0.5;
                py = MV_PAT_H * along;
                dx = -8.0;
                break;
            case 1: /* right */
                px = (gx0 + cols * cell + MV_PAT_W) * 0.5;
                py = MV_PAT_H * along;
                dx = 8.0;
                break;
            case 2: /* top */
                px = MV_PAT_W * along;
                py = gy0 * 0.5;
                dy = -8.0;
                break;
            default: /* bottom */
                px = MV_PAT_W * along;
                py = (gy0 + rows * cell + MV_PAT_H) * 0.5;
                dy = 8.0;
                break;
            }
            /* the anchor itself must be lit (it sits in the pattern's
             * white margin); otherwise skip this ray */
            proj_H(rr->H, px, py, &u, &v);
            g = sample3(img, w, h, u, v);
            if (g < thr)
                continue;
            pos = (side < 2) ? px : py;
            for (steps = 0; steps < 400 && dark < 3; steps++) {
                px += dx;
                py += dy;
                proj_H(rr->H, px, py, &u, &v);
                if (u < 1.0 || v < 1.0 || u > w - 2.0 || v > h - 2.0)
                    break; /* image border: panel extends past view */
                g = sample3(img, w, h, u, v);
                if (g < thr) {
                    dark++;
                } else {
                    dark = 0;
                    pos = (side < 2) ? px : py;
                }
            }
            if (ns < 64)
                stops[ns++] = pos;
        }
        if (ns < 8)
            return MV_ERR;
        {
            int j, k;
            for (j = 1; j < ns; j++)
                for (k = j; k > 0 && stops[k] < stops[k - 1]; k--) {
                    double t = stops[k];
                    stops[k] = stops[k - 1];
                    stops[k - 1] = t;
                }
            lim[side] = stops[ns / 2];
        }
    }
    if (lim[0] >= lim[1] || lim[2] >= lim[3])
        return MV_ERR;

    /* rectangle in pattern coords -> image quad through H */
    {
        static const int cx[4] = { 0, 1, 1, 0 }, cy[4] = { 2, 2, 3, 3 };
        for (i = 0; i < 4; i++)
            proj_H(rr->H, lim[cx[i]], lim[cy[i]], &quad[2 * i],
                   &quad[2 * i + 1]);
    }
    return MV_OK;
}
