#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mv/core.h"
#include "mv/mat.h"
#include "mv/calib.h"
#include "mv/pattern.h"
#include "mv/reader.h"

#define MAXCAND 512
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
                          int w, int h)
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

    /* NMS + threshold + sub-pixel */
    for (y = SADR + 2; y < h - SADR - 2 && n < MAXCAND; y++)
        for (x = SADR + 2; x < w - SADR - 2 && n < MAXCAND; x++) {
            float r0 = R[y * w + x];
            int i, j, peak = 1;
            double dxs, dys, den;
            if (r0 < 0.55 * rmax || r0 < 6.0)
                continue;
            for (j = -2; j <= 2 && peak; j++)
                for (i = -2; i <= 2 && peak; i++)
                    if ((i || j) && R[(y + j) * w + (x + i)] > r0)
                        peak = 0;
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
            cd[n].u = x + dxs;
            cd[n].v = y + dys;
            cd[n].resp = r0;
            cd[n].used = 0;
            n++;
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
    char tried[MAXCAND] = { 0 };
    int attempt;
    for (attempt = 0; attempt < 8; attempt++) {
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

    memset(res, 0, sizeof(*res));
    memset(bit, -1, sizeof(bit));

    ncand = detect_saddles(cd, img, w, h);
    RDBG("reader: candidates %d\n", ncand);
    if (ncand < 12)
        return MV_ERR;
    ngrid = grow_grid(&L, cd, ncand);
    RDBG("reader: lattice corners %d\n", ngrid);
    if (ngrid < 12)
        return MV_ERR;

    /* --- decode cell bits (lattice square (i,j): corners i..i+1,j..j+1) */
    {
        double bases[GRIDN * GRIDN];
        int nb = 0;
        double bmin = 255.0, bmax = 0.0;
        for (j = 0; j < GRIDN - 1; j++)
            for (i = 0; i < GRIDN - 1; i++) {
                double cu, cv, bs[4], b, dot;
                int k;
                if (!L.have[i][j] || !L.have[i + 1][j] || !L.have[i][j + 1]
                    || !L.have[i + 1][j + 1])
                    continue;
                cu = 0.25 * (L.p[i][j][0] + L.p[i + 1][j][0]
                             + L.p[i][j + 1][0] + L.p[i + 1][j + 1][0]);
                cv = 0.25 * (L.p[i][j][1] + L.p[i + 1][j][1]
                             + L.p[i][j + 1][1] + L.p[i + 1][j + 1][1]);
                for (k = 0; k < 4; k++) {
                    const double *q =
                        L.p[i + (k & 1)][j + (k >> 1)];
                    bs[k] = bilin(img, w, h, cu + 0.32 * (q[0] - cu),
                                  cv + 0.32 * (q[1] - cv));
                }
                b = med4(bs[0], bs[1], bs[2], bs[3]);
                dot = sample3(img, w, h, cu, cv);
                if (b < 0.0 || dot < 0.0)
                    continue;
                base[i][j] = b;
                bases[nb++] = b;
                if (b < bmin)
                    bmin = b;
                if (b > bmax)
                    bmax = b;
                bit[i][j] = (fabs(dot - b) > 0.30 * 200.0) ? 1 : 0;
            }
        RDBG("reader: cells %d, base range %.0f..%.0f\n", nb, bmin, bmax);
        if (nb < 16 || bmax - bmin < 60.0)
            return MV_ERR;
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
        int best_m = -1, best_u = 0, best_v = 0, best_n = 0;
        int vu[256], vv[256], vm[256], vc[256], nv = 0, k;
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
            return MV_ERR;
        res->rot = best_m;

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
                    res->uv[2 * res->n] = L.p[i][j][0];
                    res->uv[2 * res->n + 1] = L.p[i][j][1];
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
