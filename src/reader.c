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

static int detect_saddles(struct cand *cd, const unsigned char *img,
                          int w, int h, double relthr)
{
    static struct offs o;
    static int o_init = 0;
    float *R = (float *)calloc((size_t)w * h, sizeof(float));
    double rmax = 0.0;
    int x, y, k, p, n = 0;

    if (!R)
        return 0;
    if (!o_init) {
        build_offsets(&o);
        o_init = 1;
    }
    for (y = SADR + 1; y < h - SADR - 1; y++)
        for (x = SADR + 1; x < w - SADR - 1; x++) {
            double best = 0.0;
            for (p = 0; p < 4; p++) {
                double s = 0.0;
                for (k = 0; k < o.n; k++)
                    s += o.t[p][k]
                       * img[(y + o.dy[k]) * w + (x + o.dx[k])];
                if (fabs(s) > best)
                    best = fabs(s);
            }
            best /= o.n;
            R[y * w + x] = (float)best;
            if (best > rmax)
                rmax = best;
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
    free(R);
    return n;
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

/* ---------------- main entry ------------------------------------------ */

int mv_read_pattern(mv_read_result *res, const unsigned char *img,
                    int w, int h)
{
    static struct cand cd[MAXCAND];
    static struct lattice L;
    static signed char bit[GRIDN][GRIDN];   /* -1 unknown */
    static double base[GRIDN][GRIDN];
    int ncand, ngrid, i, j;
    double blackref = 255.0, whiteref = 0.0;

    /* escalation ladder: threshold x scale; small-in-frame patterns
     * (cells ~10 px) need the 2x upsample to re-enter the validated
     * ~20 px regime */
    static const double RELTHR[4] = { 0.55, 0.30, 0.55, 0.30 };
    static const int SCALE[4] = { 1, 1, 2, 2 };
    unsigned char *up = NULL;
    const unsigned char *im;
    int iw, ih, scale = 1;
    int pass, decoded = 0;
    int best_m = -1, best_u = 0, best_v = 0;

    memset(res, 0, sizeof(*res));

    /* Detection-threshold escalation: the strict threshold suppresses
     * dot L-corners when the pattern is large in frame (the sim regime);
     * when the pattern is small, blur can drop alternate saddle
     * polarities below it, leaving only the 45-degree sublattice -- so
     * on decode failure retry with the permissive threshold (tiny dots
     * cannot produce significant L-corner clutter). */
    for (pass = 0; pass < 4 && !decoded; pass++) {
    memset(bit, -1, sizeof(bit));

    scale = SCALE[pass];
    if (scale == 1) {
        im = img;
        iw = w;
        ih = h;
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
    ncand = detect_saddles(cd, im, iw, ih, RELTHR[pass]);
    RDBG("reader: pass %d (scale %d) candidates %d\n", pass, scale,
         ncand);
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
    up = NULL;
    if (!decoded)
        return MV_ERR;
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
                    res->id[res->n] = pj_min * MV_PAT_CORNER_COLS + pi_min;
                    res->uv[2 * res->n] = L.p[i][j][0] / scale;
                    res->uv[2 * res->n + 1] = L.p[i][j][1] / scale;
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
        {
            unsigned bits[MV_PAT2_CTR_CELLS], g2 = 0, parity = 0;
            double conf = 1.0, halfrange = 0.5 * (mw - mb);
            int okc = 1;
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

int mv_read_coarse(mv_read_result *res, const unsigned char *img,
                   int w, int h)
{
    static const int SC[3] = { 1, 2, 4 };
    static const double THR[2] = { 0.55, 0.30 };
    int si, ti, i;

    memset(res, 0, sizeof(*res));
    for (si = 0; si < 3; si++) {
        int s = SC[si], ws = w / s, hs = h / s;
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
        for (ti = 0; ti < 2; ti++) {
            if (coarse_try(res, ds, ws, hs, THR[ti]) == MV_OK) {
                if (s != 1) {
                    static double obj[2 * 10];
                    for (i = 0; i < res->n; i++) {
                        double xy[2];
                        res->uv[2 * i] = res->uv[2 * i] * s
                                         + (s - 1) * 0.5;
                        res->uv[2 * i + 1] = res->uv[2 * i + 1] * s
                                             + (s - 1) * 0.5;
                        mv_pattern2_corner_px(res->id[i]
                                              % MV_PAT2_CORNER_COLS,
                                              res->id[i]
                                              / MV_PAT2_CORNER_COLS, xy);
                        obj[2 * i] = xy[0];
                        obj[2 * i + 1] = xy[1];
                    }
                    mv_homography_dlt(res->H, obj, res->uv, res->n);
                    free(ds);
                }
                return MV_OK;
            }
        }
        if (s != 1)
            free(ds);
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
