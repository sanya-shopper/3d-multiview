#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mv/core.h"
#include "mv/pattern.h"
#include "mv/session.h"

/* ---- record emission ---------------------------------------------- */

void mv_session_ver(FILE *f, int spec_version, double pitch_mm,
                    const char *cam_id)
{
    fprintf(f, "VER spec=%d w=%d h=%d pitch_mm=%g Td_ms=%g",
            spec_version, MV_PAT_W, MV_PAT_H,
            pitch_mm, MV_SESSION_TD_MS);
    if (cam_id)
        fprintf(f, " cam=%s", cam_id);
    fputc('\n', f);
}

void mv_session_frm(FILE *f, const char *cam_id, double t_mono,
                    int frame_idx)
{
    fprintf(f, "FRM cam=%s k=%d t=%.6f\n", cam_id, frame_idx, t_mono);
}

void mv_session_read(FILE *f, const char *cam_id, double t_mono,
                     const mv_read_result *rr, const double phs[5])
{
    int i;
    if (rr->counter_valid)
        fprintf(f, "CTR cam=%s t=%.6f count=%u conf=%.4f\n",
                cam_id, t_mono, rr->counter, rr->counter_conf);
    for (i = 0; i < rr->n; i++)
        fprintf(f, "CRN cam=%s t=%.6f id=%d u=%.3f v=%.3f\n",
                cam_id, t_mono, rr->id[i],
                rr->uv[2 * i], rr->uv[2 * i + 1]);
    if (phs)
        for (i = 0; i < 5; i++)
            fprintf(f, "PHS cam=%s t=%.6f patch=%d g=%.6f\n",
                    cam_id, t_mono, i, phs[i]);
}

/* ---- record parsing ------------------------------------------------ */

int mv_session_parse(mv_session_record *rec, const char *line)
{
    const char *p = line;
    int i;

    rec->type[0] = '\0';
    rec->nf = 0;
    while (*p == ' ' || *p == '\t')
        p++;
    if (*p == '\0' || *p == '\n' || *p == '#')
        return MV_OK; /* blank/comment: type "" */

    for (i = 0; i < 3 && isalnum((unsigned char)p[i]); i++)
        rec->type[i] = p[i];
    rec->type[i] = '\0';
    if (i == 0 || (p[i] != '\0' && !isspace((unsigned char)p[i])))
        return MV_ERR;
    p += i;

    for (;;) {
        const char *eq, *end;
        size_t kl, vl;
        while (*p == ' ' || *p == '\t')
            p++;
        if (*p == '\0' || *p == '\n' || *p == '\r')
            return MV_OK;
        for (end = p; *end && !isspace((unsigned char)*end); end++)
            ;
        eq = memchr(p, '=', (size_t)(end - p));
        if (!eq || eq == p)
            return MV_ERR;
        kl = (size_t)(eq - p);
        vl = (size_t)(end - eq - 1);
        if (rec->nf >= MV_SESSION_MAXF || kl >= MV_SESSION_KEYSZ
            || vl >= MV_SESSION_VALSZ)
            return MV_ERR;
        memcpy(rec->key[rec->nf], p, kl);
        rec->key[rec->nf][kl] = '\0';
        memcpy(rec->val[rec->nf], eq + 1, vl);
        rec->val[rec->nf][vl] = '\0';
        rec->nf++;
        p = end;
    }
}

const char *mv_session_field(const mv_session_record *rec,
                             const char *key)
{
    int i;
    for (i = 0; i < rec->nf; i++)
        if (strcmp(rec->key[i], key) == 0)
            return rec->val[i];
    return NULL;
}

int mv_session_num(const mv_session_record *rec, const char *key,
                   double *out)
{
    const char *v = mv_session_field(rec, key);
    char *end;
    if (!v)
        return MV_ERR;
    *out = strtod(v, &end);
    return (end != v && *end == '\0') ? MV_OK : MV_ERR;
}

/* ---- reader-side samplers ----------------------------------------- */

static double bilin(const unsigned char *img, int w, int h, double u,
                    double v)
{
    int x = (int)u, y = (int)v;
    double fx = u - x, fy = v - y;
    const unsigned char *p;
    if (u < 0.0 || v < 0.0 || x >= w - 1 || y >= h - 1)
        return -1.0;
    p = img + y * w + x;
    return (1 - fx) * (1 - fy) * p[0] + fx * (1 - fy) * p[1]
         + (1 - fx) * fy * p[w] + fx * fy * p[w + 1];
}

/* 3x3 mean of bilinear samples; -1 if every tap is out of bounds */
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

static void proj_H(const double H[9], double X, double Y, double *u,
                   double *v)
{
    double W = H[6] * X + H[7] * Y + H[8];
    if (fabs(W) < 1e-12)
        W = 1e-12;
    *u = (H[0] * X + H[1] * Y + H[2]) / W;
    *v = (H[3] * X + H[4] * Y + H[5]) / W;
}

/* white/black references from the checker cells of the FINE pattern,
 * sampled off center (0.30 of a cell from the corner) to dodge the
 * M-array dots. MV_ERR if too few cells land in the image or the
 * contrast is below 20 gray levels. */
static int pattern_refs(double *blackref, double *whiteref,
                        const mv_read_result *rr,
                        const unsigned char *img, int w, int h)
{
    double sb = 0.0, sw = 0.0;
    int nb = 0, nw = 0, r, c;
    for (r = 0; r < MV_PAT_GRID_ROWS; r++)
        for (c = 0; c < MV_PAT_GRID_COLS; c++) {
            double u, v, g;
            proj_H(rr->H, MV_PAT_GRID_X0 + (c + 0.30) * MV_PAT_CELL,
                   MV_PAT_GRID_Y0 + (r + 0.30) * MV_PAT_CELL, &u, &v);
            g = sample3(img, w, h, u, v);
            if (g < 0.0)
                continue;
            if ((r + c) % 2 == 0) {
                sb += g;
                nb++;
            } else {
                sw += g;
                nw++;
            }
        }
    if (nb < 4 || nw < 4)
        return MV_ERR;
    *blackref = sb / nb;
    *whiteref = sw / nw;
    return (*whiteref - *blackref >= 20.0) ? MV_OK : MV_ERR;
}

int mv_read_version(const mv_read_result *rr, const unsigned char *img,
                    int w, int h)
{
    double blackref, whiteref, mid;
    unsigned bits[MV_PAT_VER_CELLS], g = 0, b;
    int c;

    if (pattern_refs(&blackref, &whiteref, rr, img, w, h) != MV_OK)
        return -1;
    mid = 0.5 * (blackref + whiteref);
    for (c = 0; c < MV_PAT_VER_CELLS; c++) {
        double xy[2], u, v, s;
        mv_pattern_version_cell_px(c, xy);
        proj_H(rr->H, xy[0], xy[1], &u, &v);
        s = sample3(img, w, h, u, v);
        if (s < 0.0)
            return -1;
        bits[c] = (s < mid) ? 1u : 0u; /* black = bit 1 */
    }
    if (bits[0] != 1u || bits[1] != 0u)
        return -1;
    for (c = 2; c < MV_PAT_VER_CELLS; c++)
        g = (g << 1) | bits[c]; /* MSB first */
    b = g;
    b ^= b >> 1;
    b ^= b >> 2;
    b ^= b >> 4;
    return (int)b;
}

int mv_read_phases(double phs[5], const mv_read_result *rr,
                   const unsigned char *img, int w, int h)
{
    double blackref, whiteref;
    int i;

    if (pattern_refs(&blackref, &whiteref, rr, img, w, h) != MV_OK)
        return MV_ERR;
    for (i = 0; i < 5; i++) {
        double xy[2], u, v, s, g;
        mv_pattern_patch_px(i, xy);
        proj_H(rr->H, xy[0], xy[1], &u, &v);
        s = sample3(img, w, h, u, v);
        if (s < 0.0)
            return MV_ERR;
        g = (s - blackref) / (whiteref - blackref);
        if (g < 0.0)
            g = 0.0;
        if (g > 1.0)
            g = 1.0;
        phs[i] = g;
    }
    return MV_OK;
}

/* ---- robust clock fit --------------------------------------------- */

static int cmp_double(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

static double median_inplace(double *v, int n)
{
    qsort(v, (size_t)n, sizeof(double), cmp_double);
    return (n % 2) ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

#define CLOCK_SUB 200 /* Theil-Sen subsample cap */

int mv_clock_fit(double *phi, double *T, const double *t,
                 const double *k, int n)
{
    static double slopes[CLOCK_SUB * (CLOCK_SUB - 1) / 2];
    double *ws;
    double T0, phi0, medr, mad, thr;
    int i, j, m, ns = 0, nin = 0;
    int idx[CLOCK_SUB];

    if (n < 2)
        return MV_ERR;

    /* pairwise slopes on <= CLOCK_SUB evenly spaced samples
     * (deterministic; preserves the full time span) */
    m = (n > CLOCK_SUB) ? CLOCK_SUB : n;
    for (i = 0; i < m; i++)
        idx[i] = (int)(((long long)i * (n - 1)) / (m - 1));
    for (i = 0; i < m; i++)
        for (j = i + 1; j < m; j++) {
            double dk = k[idx[j]] - k[idx[i]];
            if (fabs(dk) > 0.0)
                slopes[ns++] = (t[idx[j]] - t[idx[i]]) / dk;
        }
    if (ns == 0)
        return MV_ERR;
    T0 = median_inplace(slopes, ns);

    ws = (double *)malloc((size_t)n * sizeof(double));
    if (!ws)
        return MV_ERR;

    /* intercept: median over ALL samples */
    for (i = 0; i < n; i++)
        ws[i] = t[i] - T0 * k[i];
    phi0 = median_inplace(ws, n);

    /* 3-MAD gate on residuals (floor guards the noise-free case) */
    for (i = 0; i < n; i++)
        ws[i] = t[i] - phi0 - T0 * k[i];
    {
        double *tmp = (double *)malloc((size_t)n * sizeof(double));
        if (!tmp) {
            free(ws);
            return MV_ERR;
        }
        memcpy(tmp, ws, (size_t)n * sizeof(double));
        medr = median_inplace(tmp, n);
        for (i = 0; i < n; i++)
            tmp[i] = fabs(ws[i] - medr);
        mad = median_inplace(tmp, n);
        free(tmp);
    }
    thr = 3.0 * mad;
    if (thr < 1e-9)
        thr = 1e-9;

    /* least squares on inliers, centered for conditioning */
    {
        double sk = 0.0, st = 0.0, skk = 0.0, skt = 0.0, kb, tb;
        for (i = 0; i < n; i++)
            if (fabs(ws[i] - medr) <= thr) {
                sk += k[i];
                st += t[i];
                nin++;
            }
        if (nin < 2) {
            free(ws);
            return MV_ERR;
        }
        kb = sk / nin;
        tb = st / nin;
        for (i = 0; i < n; i++)
            if (fabs(ws[i] - medr) <= thr) {
                skk += (k[i] - kb) * (k[i] - kb);
                skt += (k[i] - kb) * (t[i] - tb);
            }
        free(ws);
        if (skk <= 0.0)
            return MV_ERR;
        *T = skt / skk;
        *phi = tb - *T * kb;
    }
    return MV_OK;
}
