#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mv/core.h"
#include "mv/mat.h"
#include "mv/tsdf.h"

int mv_tsdf_init(mv_tsdf *t, double x0, double y0, double z0,
                 double sx, double sy, double sz, double voxel, double tau)
{
    size_t n;
    int i;
    if (voxel <= 0.0 || tau <= 0.0 || sx <= 0.0 || sy <= 0.0 || sz <= 0.0)
        return MV_ERR;
    t->nx = (int)ceil(sx / voxel) + 1;
    t->ny = (int)ceil(sy / voxel) + 1;
    t->nz = (int)ceil(sz / voxel) + 1;
    t->x0 = x0;
    t->y0 = y0;
    t->z0 = z0;
    t->voxel = voxel;
    t->tau = tau;
    n = (size_t)t->nx * t->ny * t->nz;
    t->d = (float *)malloc(n * sizeof(float));
    t->w = (float *)calloc(n, sizeof(float));
    if (!t->d || !t->w) {
        free(t->d);
        free(t->w);
        return MV_ERR;
    }
    for (i = 0; (size_t)i < n; i++)
        t->d[i] = (float)tau;
    return MV_OK;
}

void mv_tsdf_free(mv_tsdf *t)
{
    free(t->d);
    free(t->w);
    memset(t, 0, sizeof(*t));
}

static long vidx(const mv_tsdf *t, int i, int j, int k)
{
    return ((long)k * t->ny + j) * t->nx + i;
}

void mv_tsdf_fuse(mv_tsdf *t, const double p[3], const double o[3],
                  double wgt)
{
    double dir[3], L, s, step;
    int a;
    for (a = 0; a < 3; a++)
        dir[a] = p[a] - o[a];
    L = mv_norm(dir, 3);
    if (L < 1e-9 || wgt <= 0.0)
        return;
    for (a = 0; a < 3; a++)
        dir[a] /= L;
    step = 0.7 * t->voxel;
    for (s = L - t->tau; s <= L + t->tau; s += step) {
        double q[3], eta = L - s;
        int i, j, k;
        long id;
        for (a = 0; a < 3; a++)
            q[a] = o[a] + s * dir[a];
        i = (int)floor((q[0] - t->x0) / t->voxel + 0.5);
        j = (int)floor((q[1] - t->y0) / t->voxel + 0.5);
        k = (int)floor((q[2] - t->z0) / t->voxel + 0.5);
        if (i < 0 || j < 0 || k < 0 || i >= t->nx || j >= t->ny
            || k >= t->nz)
            continue;
        id = vidx(t, i, j, k);
        t->d[id] = (float)((t->w[id] * t->d[id] + wgt * eta)
                           / (t->w[id] + wgt));
        t->w[id] += (float)wgt;
    }
}

/* ---- voxel-major fast path ----------------------------------------- */

int mv_tsdf_table_build(const mv_tsdf *t, const mv_camera *cam,
                        int w, int h, mv_tsdf_table *tab)
{
    long n, id;
    int i, j, k;
    if (w <= 0 || h <= 0 || !t->d)
        return MV_ERR;
    n = (long)t->nx * t->ny * t->nz;
    tab->nvox = n;
    tab->w = w;
    tab->h = h;
    tab->pix = (int *)malloc((size_t)n * sizeof(int));
    tab->lv = (float *)malloc((size_t)n * sizeof(float));
    tab->scale = (float *)malloc((size_t)n * sizeof(float));
    if (!tab->pix || !tab->lv || !tab->scale) {
        mv_tsdf_table_free(tab);
        return MV_ERR;
    }
    id = 0;
    for (k = 0; k < t->nz; k++)
        for (j = 0; j < t->ny; j++)
            for (i = 0; i < t->nx; i++, id++) {
                double q[3], Xc[3], uv[2], L;
                int u, v, a;
                q[0] = t->x0 + i * t->voxel;
                q[1] = t->y0 + j * t->voxel;
                q[2] = t->z0 + k * t->voxel;
                tab->pix[id] = -1;
                tab->lv[id] = 0.0f;
                tab->scale[id] = 1.0f;
                for (a = 0; a < 3; a++)
                    Xc[a] = cam->R[a * 3 + 0] * q[0]
                          + cam->R[a * 3 + 1] * q[1]
                          + cam->R[a * 3 + 2] * q[2] + cam->t[a];
                if (mv_cam_project(uv, cam, q) != MV_OK)
                    continue;
                u = (int)floor(uv[0] + 0.5);
                v = (int)floor(uv[1] + 0.5);
                if (u < 0 || v < 0 || u >= w || v >= h)
                    continue;
                L = mv_norm(Xc, 3);
                tab->pix[id] = v * w + u;
                tab->lv[id] = (float)L;
                tab->scale[id] = (float)(L / Xc[2]);
            }
    return MV_OK;
}

void mv_tsdf_table_free(mv_tsdf_table *tab)
{
    free(tab->pix);
    free(tab->lv);
    free(tab->scale);
    memset(tab, 0, sizeof(*tab));
}

int mv_tsdf_fuse_depthmap(mv_tsdf *t, const mv_tsdf_table *tab,
                          const float *depth, int w, int h, double sigma0)
{
    long n = (long)t->nx * t->ny * t->nz;
    long id;
    if (w != tab->w || h != tab->h || n != tab->nvox || !tab->pix)
        return MV_ERR;
    for (id = 0; id < n; id++) {
        double z, eta, wgt;
        int p = tab->pix[id];
        if (p < 0)
            continue;
        z = depth[p];
        /* positive-form gate: NaN fails z > 0 and is rejected */
        if (!(z > 0.0) || !isfinite(z))
            continue;
        eta = z * tab->scale[id] - tab->lv[id];
        if (eta < -t->tau)
            continue; /* occluded: the ray says nothing about it */
        if (eta > t->tau)
            eta = t->tau; /* free space: carve with the clamped value */
        if (sigma0 > 0.0) {
            double sz = sigma0 * z * z;
            wgt = 1.0 / (sz * sz);
        } else {
            wgt = 1.0;
        }
        t->d[id] = (float)((t->w[id] * t->d[id] + wgt * eta)
                           / (t->w[id] + wgt));
        t->w[id] += (float)wgt;
    }
    return MV_OK;
}

double mv_tsdf_query(const mv_tsdf *t, const double q[3])
{
    double fx = (q[0] - t->x0) / t->voxel;
    double fy = (q[1] - t->y0) / t->voxel;
    double fz = (q[2] - t->z0) / t->voxel;
    int i = (int)floor(fx), j = (int)floor(fy), k = (int)floor(fz);
    double v = 0.0;
    int a, b, c;
    if (i < 0 || j < 0 || k < 0 || i >= t->nx - 1 || j >= t->ny - 1
        || k >= t->nz - 1)
        return HUGE_VAL;
    fx -= i;
    fy -= j;
    fz -= k;
    for (c = 0; c < 2; c++)
        for (b = 0; b < 2; b++)
            for (a = 0; a < 2; a++) {
                long id = vidx(t, i + a, j + b, k + c);
                if (t->w[id] <= 0.0f)
                    return HUGE_VAL;
                v += t->d[id] * (a ? fx : 1 - fx) * (b ? fy : 1 - fy)
                   * (c ? fz : 1 - fz);
            }
    return v;
}

/* ---- extraction: marching tetrahedra ------------------------------- */

/* cube corner offsets, binary order (x,y,z bits) */
static const int CO[8][3] = {
    { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 1, 1, 0 },
    { 0, 0, 1 }, { 1, 0, 1 }, { 0, 1, 1 }, { 1, 1, 1 }
};

/* six tetrahedra sharing the main diagonal 0-7 */
static const int TETS[6][4] = {
    { 0, 7, 1, 3 }, { 0, 7, 3, 2 }, { 0, 7, 2, 6 },
    { 0, 7, 6, 4 }, { 0, 7, 4, 5 }, { 0, 7, 5, 1 }
};

struct meshbuf {
    double *tris;
    int n, cap;
};

static int push_tri(struct meshbuf *m, const double a[3], const double b[3],
                    const double c[3], const double out_ref[3])
{
    double u[3], v[3], nrm[3], r[3], cen[3];
    int i;
    if (m->n == m->cap) {
        int ncap = m->cap ? m->cap * 2 : 1024;
        double *nt = (double *)realloc(m->tris,
                                       (size_t)ncap * 9 * sizeof(double));
        if (!nt)
            return MV_ERR;
        m->tris = nt;
        m->cap = ncap;
    }
    /* orient outward: normal must point toward out_ref */
    for (i = 0; i < 3; i++) {
        u[i] = b[i] - a[i];
        v[i] = c[i] - a[i];
        cen[i] = (a[i] + b[i] + c[i]) / 3.0;
        r[i] = out_ref[i] - cen[i];
    }
    mv_cross3(nrm, u, v);
    if (mv_dot(nrm, r, 3) < 0.0) {
        const double *tmp = b;
        b = c;
        c = tmp;
    }
    memcpy(m->tris + 9 * m->n + 0, a, 3 * sizeof(double));
    memcpy(m->tris + 9 * m->n + 3, b, 3 * sizeof(double));
    memcpy(m->tris + 9 * m->n + 6, c, 3 * sizeof(double));
    m->n++;
    return MV_OK;
}

static void edge_cross(double out[3], const double pa[3], double da,
                       const double pb[3], double db)
{
    double f = da / (da - db);
    int i;
    for (i = 0; i < 3; i++)
        out[i] = pa[i] + f * (pb[i] - pa[i]);
}

int mv_tsdf_mesh(const mv_tsdf *t, double **tris, int *ntri)
{
    struct meshbuf m = { NULL, 0, 0 };
    int i, j, k, e, c;

    for (k = 0; k < t->nz - 1; k++)
        for (j = 0; j < t->ny - 1; j++)
            for (i = 0; i < t->nx - 1; i++) {
                double P[8][3], D[8];
                int ok = 1;
                for (c = 0; c < 8 && ok; c++) {
                    long id = vidx(t, i + CO[c][0], j + CO[c][1],
                                   k + CO[c][2]);
                    if (t->w[id] <= 0.0f)
                        ok = 0;
                    P[c][0] = t->x0 + (i + CO[c][0]) * t->voxel;
                    P[c][1] = t->y0 + (j + CO[c][1]) * t->voxel;
                    P[c][2] = t->z0 + (k + CO[c][2]) * t->voxel;
                    D[c] = t->d[id];
                }
                if (!ok)
                    continue;
                for (e = 0; e < 6; e++) {
                    const int *tv = TETS[e];
                    int inside[4], nin = 0, v;
                    double x[4][3], cen_in[3] = { 0, 0, 0 };
                    double cen_out[3] = { 0, 0, 0 };
                    int nin_c = 0, nout_c = 0;
                    for (v = 0; v < 4; v++) {
                        inside[v] = D[tv[v]] < 0.0;
                        nin += inside[v];
                    }
                    if (nin == 0 || nin == 4)
                        continue;
                    for (v = 0; v < 4; v++) {
                        int a2;
                        for (a2 = 0; a2 < 3; a2++) {
                            if (inside[v])
                                cen_in[a2] += P[tv[v]][a2];
                            else
                                cen_out[a2] += P[tv[v]][a2];
                        }
                        if (inside[v])
                            nin_c++;
                        else
                            nout_c++;
                    }
                    for (v = 0; v < 3; v++) {
                        cen_in[v] /= nin_c;
                        cen_out[v] /= nout_c;
                    }
                    if (nin == 1 || nin == 3) {
                        int lone = -1, others[3], no = 0;
                        for (v = 0; v < 4; v++) {
                            if (inside[v] == (nin == 1))
                                lone = v;
                            else
                                others[no++] = v;
                        }
                        for (v = 0; v < 3; v++)
                            edge_cross(x[v], P[tv[lone]], D[tv[lone]],
                                       P[tv[others[v]]], D[tv[others[v]]]);
                        if (push_tri(&m, x[0], x[1], x[2], cen_out)
                            != MV_OK)
                            goto fail;
                    } else { /* 2 in, 2 out: quad -> two triangles */
                        int in2[2], out2[2], ni = 0, no2 = 0;
                        for (v = 0; v < 4; v++) {
                            if (inside[v])
                                in2[ni++] = v;
                            else
                                out2[no2++] = v;
                        }
                        edge_cross(x[0], P[tv[in2[0]]], D[tv[in2[0]]],
                                   P[tv[out2[0]]], D[tv[out2[0]]]);
                        edge_cross(x[1], P[tv[in2[0]]], D[tv[in2[0]]],
                                   P[tv[out2[1]]], D[tv[out2[1]]]);
                        edge_cross(x[2], P[tv[in2[1]]], D[tv[in2[1]]],
                                   P[tv[out2[1]]], D[tv[out2[1]]]);
                        edge_cross(x[3], P[tv[in2[1]]], D[tv[in2[1]]],
                                   P[tv[out2[0]]], D[tv[out2[0]]]);
                        if (push_tri(&m, x[0], x[1], x[2], cen_out)
                            != MV_OK)
                            goto fail;
                        if (push_tri(&m, x[0], x[2], x[3], cen_out)
                            != MV_OK)
                            goto fail;
                    }
                }
            }
    *tris = m.tris;
    *ntri = m.n;
    return MV_OK;
fail:
    free(m.tris);
    return MV_ERR;
}

int mv_tsdf_write_ply(const mv_tsdf *t, const char *path)
{
    double *tris;
    int ntri, i;
    FILE *f;
    if (mv_tsdf_mesh(t, &tris, &ntri) != MV_OK)
        return MV_ERR;
    f = fopen(path, "w");
    if (!f) {
        free(tris);
        return MV_ERR;
    }
    fprintf(f, "ply\nformat ascii 1.0\nelement vertex %d\n", 3 * ntri);
    fprintf(f, "property float x\nproperty float y\nproperty float z\n");
    fprintf(f, "element face %d\n", ntri);
    fprintf(f, "property list uchar int vertex_indices\nend_header\n");
    for (i = 0; i < 3 * ntri; i++)
        fprintf(f, "%g %g %g\n", tris[3 * i], tris[3 * i + 1],
                tris[3 * i + 2]);
    for (i = 0; i < ntri; i++)
        fprintf(f, "3 %d %d %d\n", 3 * i, 3 * i + 1, 3 * i + 2);
    fclose(f);
    free(tris);
    return MV_OK;
}
