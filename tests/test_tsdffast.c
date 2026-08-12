/* Voxel-major TSDF fast path: table sweep vs analytic truth and vs the
 * sample fold, plus free-space carving and a sweep/fold timing ratio.
 * OWNERSHIP (parallel build): TSDF-completion work item ONLY. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "mv/core.h"
#include "mv/mat.h"
#include "mv/cam.h"
#include "mv/tsdf.h"

static int failures = 0;

#define CHECK(cond, name) \
    do { \
        if (cond) \
            printf("ok:   %s\n", name); \
        else { \
            printf("FAIL: %s\n", name); \
            failures++; \
        } \
    } while (0)

#define IMW 640
#define IMH 480

static void make_cam(mv_camera *c)
{
    mv_cam_set_K(c, 600.0, 600.0, 319.5, 239.5);
    mv_cam_set_identity_pose(c);
    memset(c->k, 0, sizeof(c->k));
}

/* depth map of the plane z = z0 (fronto-parallel: constant depth) */
static void depth_plane(float *depth, double z0)
{
    int i;
    for (i = 0; i < IMW * IMH; i++)
        depth[i] = (float)z0;
}

/* depth map of a sphere at ctr with radius r; 0 = miss */
static void depth_sphere(float *depth, const mv_camera *cam,
                         const double ctr[3], double r)
{
    int x, y;
    for (y = 0; y < IMH; y++)
        for (x = 0; x < IMW; x++) {
            double uv[2] = { (double)x, (double)y };
            double orig[3], dir[3], oc[3], b, c, disc, s;
            int a;
            depth[y * IMW + x] = 0.0f;
            if (mv_cam_ray(orig, dir, cam, uv) != MV_OK)
                continue;
            for (a = 0; a < 3; a++)
                oc[a] = orig[a] - ctr[a];
            b = mv_dot(oc, dir, 3);
            c = mv_dot(oc, oc, 3) - r * r;
            disc = b * b - c;
            if (disc < 0.0)
                continue;
            s = -b - sqrt(disc);
            if (s <= 0.0)
                continue;
            depth[y * IMW + x] = (float)(s * dir[2]);
        }
}

/* fold the same depth map through the enriched-sample interface */
static void fold_depthmap(mv_tsdf *t, const mv_camera *cam,
                          const float *depth, double wgt)
{
    int x, y;
    for (y = 0; y < IMH; y++)
        for (x = 0; x < IMW; x++) {
            double uv[2] = { (double)x, (double)y };
            double orig[3], dir[3], p[3], z, s;
            int a;
            z = depth[y * IMW + x];
            if (!(z > 0.0))
                continue;
            if (mv_cam_ray(orig, dir, cam, uv) != MV_OK)
                continue;
            s = z / dir[2]; /* depth -> along-ray distance */
            for (a = 0; a < 3; a++)
                p[a] = orig[a] + s * dir[a];
            mv_tsdf_fuse(t, p, orig, wgt);
        }
}

static long vid(const mv_tsdf *t, int i, int j, int k)
{
    return ((long)k * t->ny + j) * t->nx + i;
}

/* ---- (a),(b): known plane, analytic signed distance + carving ------- */

static void test_plane_scene(void)
{
    mv_camera cam;
    mv_tsdf t;
    mv_tsdf_table tab;
    static float depth[IMW * IMH];
    const double z0 = 2.0, voxel = 0.01, tau = 0.03;
    double band_max = 0.0;
    int i, j, k, nband = 0, carve_bad = 0, ncarve = 0;
    int behind_bad = 0, nbehind = 0;

    make_cam(&cam);
    CHECK(mv_tsdf_init(&t, -0.4, -0.4, 1.7, 0.8, 0.8, 0.6, voxel, tau)
          == MV_OK, "plane: tsdf init");
    CHECK(mv_tsdf_table_build(&t, &cam, IMW, IMH, &tab) == MV_OK,
          "plane: table build");
    depth_plane(depth, z0);
    CHECK(mv_tsdf_fuse_depthmap(&t, &tab, depth, IMW, IMH, 0.0) == MV_OK,
          "plane: fuse depth map");

    for (k = 0; k < t.nz; k++)
        for (j = 0; j < t.ny; j++)
            for (i = 0; i < t.nx; i++) {
                long id = vid(&t, i, j, k);
                double zv = t.z0 + k * voxel;
                double sd = z0 - zv; /* Euclidean signed distance */
                if (tab.pix[id] < 0)
                    continue;
                /* margins keep grid points that land exactly on the
                 * truncation boundary out of every bucket */
                if (fabs(sd) < tau - 1e-6) {
                    double e = fabs(t.d[id] - sd);
                    if (t.w[id] > 0.0f) {
                        nband++;
                        if (e > band_max)
                            band_max = e;
                    }
                } else if (sd > tau + 1e-6) { /* in front: carved */
                    ncarve++;
                    if (!(t.w[id] > 0.0f)
                        || fabs(t.d[id] - tau) > 1e-7)
                        carve_bad++;
                } else if (sd < -tau - 1e-6) { /* behind: untouched */
                    nbehind++;
                    if (t.w[id] != 0.0f)
                        behind_bad++;
                }
            }
    printf("plane: %d band voxels, max |D - sd| = %.2f mm "
           "(half voxel = %.2f mm)\n", nband, band_max * 1e3,
           0.5 * voxel * 1e3);
    CHECK(nband > 10000 && band_max < 0.5 * voxel,
          "plane: band D matches analytic sd to half a voxel");
    CHECK(ncarve > 10000 && carve_bad == 0,
          "plane: observed free space carved to +tau");
    CHECK(nbehind > 10000 && behind_bad == 0,
          "plane: voxels beyond tau behind the surface untouched");

    /* zero level set: signed queries and the axis zero crossing */
    {
        double qf[3] = { 0.1, -0.1, z0 - 0.01 };
        double qb[3] = { 0.1, -0.1, z0 + 0.01 };
        double vf = mv_tsdf_query(&t, qf), vb = mv_tsdf_query(&t, qb);
        CHECK(vf < HUGE_VAL && vf > 0.0 && vb < HUGE_VAL && vb < 0.0,
              "plane: query signs straddle the surface");
    }
    {
        int ic = t.nx / 2, jc = t.ny / 2, found = 0;
        for (k = 0; k + 1 < t.nz && !found; k++) {
            float d0 = t.d[vid(&t, ic, jc, k)];
            float d1 = t.d[vid(&t, ic, jc, k + 1)];
            if (t.w[vid(&t, ic, jc, k)] > 0.0f
                && t.w[vid(&t, ic, jc, k + 1)] > 0.0f
                && d0 >= 0.0f && d1 < 0.0f) {
                double zc = t.z0 + k * voxel
                    + voxel * d0 / (d0 - d1);
                CHECK(fabs(zc - z0) < 0.5 * voxel,
                      "plane: axis zero crossing at z0");
                found = 1;
            }
        }
        CHECK(found, "plane: zero crossing exists on the axis");
    }
    mv_tsdf_table_free(&tab);
    mv_tsdf_free(&t);
}

/* ---- (a),(b),(c): sphere, sweep vs fold agreement + carving --------- */

static void test_sphere_scene(void)
{
    mv_camera cam;
    mv_tsdf ts, tf;
    mv_tsdf_table tab;
    static float depth[IMW * IMH];
    const double ctr[3] = { 0.0, 0.0, 2.0 };
    const double r = 0.25, voxel = 0.01, tau = 0.03;
    double agree_max = 0.0, agree_se = 0.0, graze_max = 0.0;
    double cap_max = 0.0;
    int i, j, k, nagree = 0, nsteep = 0, ncap = 0;

    make_cam(&cam);
    CHECK(mv_tsdf_init(&ts, -0.35, -0.35, 1.65, 0.7, 0.7, 0.7, voxel,
                       tau) == MV_OK, "sphere: sweep tsdf init");
    CHECK(mv_tsdf_init(&tf, -0.35, -0.35, 1.65, 0.7, 0.7, 0.7, voxel,
                       tau) == MV_OK, "sphere: fold tsdf init");
    CHECK(mv_tsdf_table_build(&ts, &cam, IMW, IMH, &tab) == MV_OK,
          "sphere: table build");
    depth_sphere(depth, &cam, ctr, r);
    CHECK(mv_tsdf_fuse_depthmap(&ts, &tab, depth, IMW, IMH, 0.0)
          == MV_OK, "sphere: sweep fuse");
    fold_depthmap(&tf, &cam, depth, 1.0);

    for (k = 0; k < ts.nz; k++)
        for (j = 0; j < ts.ny; j++)
            for (i = 0; i < ts.nx; i++) {
                long id = vid(&ts, i, j, k);
                double q[3], rho, sd;
                q[0] = ts.x0 + i * voxel;
                q[1] = ts.y0 + j * voxel;
                q[2] = ts.z0 + k * voxel;
                rho = sqrt((q[0] - ctr[0]) * (q[0] - ctr[0])
                           + (q[1] - ctr[1]) * (q[1] - ctr[1])
                           + (q[2] - ctr[2]) * (q[2] - ctr[2]));
                sd = rho - r; /* + outside = in front, near the cap */
                /* (c) sweep vs fold on common band voxels.  Near the
                 * silhouette the two traversals read almost-parallel
                 * grazing rays whose depths differ by many mm per
                 * pixel, so the pointwise bound is stated away from
                 * grazing incidence (< 60 deg) and the rms over all. */
                if (ts.w[id] > 0.0f && tf.w[id] > 0.0f
                    && fabs(ts.d[id]) < 0.9 * tau
                    && fabs(tf.d[id]) < 0.9 * tau) {
                    double e = fabs(ts.d[id] - tf.d[id]);
                    double L = mv_norm(q, 3), ci = 0.0;
                    nagree++;
                    agree_se += e * e;
                    if (rho > 1e-9 && L > 1e-9)
                        ci = fabs(((q[0] - ctr[0]) * q[0]
                                   + (q[1] - ctr[1]) * q[1]
                                   + (q[2] - ctr[2]) * q[2])
                                  / (rho * L));
                    if (ci > 0.5) {
                        nsteep++;
                        if (e > agree_max)
                            agree_max = e;
                    } else if (e > graze_max) {
                        graze_max = e;
                    }
                }
                /* (a) analytic check on the near-normal front cap */
                if (ts.w[id] > 0.0f && fabs(sd) < 0.9 * tau
                    && q[2] < ctr[2]
                    && q[0] * q[0] + q[1] * q[1] < 0.08 * 0.08) {
                    double e = fabs(ts.d[id] - sd);
                    ncap++;
                    if (e > cap_max)
                        cap_max = e;
                }
            }
    printf("sphere: %d common band voxels (%d non-grazing), sweep-fold "
           "rms %.2f mm, max %.2f mm non-grazing / %.2f mm grazing; "
           "%d cap voxels, max |D - sd| = %.2f mm\n", nagree, nsteep,
           nagree ? sqrt(agree_se / nagree) * 1e3 : 0.0,
           agree_max * 1e3, graze_max * 1e3, ncap, cap_max * 1e3);
    CHECK(nagree > 1000 && nsteep > 1000
          && agree_max < 1.5 * voxel
          && sqrt(agree_se / nagree) < 0.5 * voxel,
          "sphere: sweep and fold agree on band voxels");
    CHECK(ncap > 100 && cap_max < 0.5 * voxel,
          "sphere: front-cap D matches analytic sd to half a voxel");

    /* (b) carving: on-axis voxel well in front of the sphere is carved
     * by the sweep, untouched by the band-only fold; the sphere's
     * center (far behind the surface) is untouched by both */
    {
        long id_front = vid(&ts, ts.nx / 2, ts.ny / 2, 1);
        long id_ctr = vid(&ts, ts.nx / 2, ts.ny / 2,
                          (int)floor((ctr[2] - ts.z0) / voxel + 0.5));
        CHECK(ts.w[id_front] > 0.0f
              && fabs(ts.d[id_front] - tau) < 1e-7,
              "sphere: sweep carves observed free space to +tau");
        CHECK(tf.w[id_front] == 0.0f,
              "sphere: band-only fold leaves free space untouched");
        CHECK(ts.w[id_ctr] == 0.0f && tf.w[id_ctr] == 0.0f,
              "sphere: interior beyond tau stays unobserved");
    }
    /* zero level set on the axis: front surface at ctr.z - r */
    {
        double qf[3] = { 0.0, 0.0, ctr[2] - r - 0.015 };
        double qb[3] = { 0.0, 0.0, ctr[2] - r + 0.015 };
        double vf = mv_tsdf_query(&ts, qf), vb = mv_tsdf_query(&ts, qb);
        CHECK(vf < HUGE_VAL && vf > 0.0 && vb < HUGE_VAL && vb < 0.0,
              "sphere: sweep zero level set at the front surface");
    }
    mv_tsdf_table_free(&tab);
    mv_tsdf_free(&ts);
    mv_tsdf_free(&tf);
}

/* ---- weights: sigma0 depth law ------------------------------------- */

static void test_weights(void)
{
    mv_camera cam;
    mv_tsdf t;
    mv_tsdf_table tab;
    static float depth[IMW * IMH];
    const double z0 = 2.0, sigma0 = 0.005;
    double expect;
    long id;

    make_cam(&cam);
    CHECK(mv_tsdf_init(&t, -0.1, -0.1, 1.9, 0.2, 0.2, 0.2, 0.01, 0.03)
          == MV_OK, "weights: init");
    CHECK(mv_tsdf_table_build(&t, &cam, IMW, IMH, &tab) == MV_OK,
          "weights: table build");
    depth_plane(depth, z0);
    CHECK(mv_tsdf_fuse_depthmap(&t, &tab, depth, IMW, IMH, sigma0)
          == MV_OK, "weights: fuse");
    /* every observed voxel saw one sample of weight 1/(sigma0*z0^2)^2 */
    expect = 1.0 / (sigma0 * z0 * z0 * sigma0 * z0 * z0);
    id = vid(&t, t.nx / 2, t.ny / 2, t.nz / 2);
    CHECK(t.w[id] > 0.0f
          && fabs(t.w[id] - expect) < 1e-4 * expect,
          "weights: W = 1/sigma_Z^2 with sigma_Z = sigma0 Z^2");
    /* invalid pixels are no measurement */
    {
        mv_tsdf t2;
        long nvox, i, nobs = 0;
        CHECK(mv_tsdf_init(&t2, -0.1, -0.1, 1.9, 0.2, 0.2, 0.2, 0.01,
                           0.03) == MV_OK, "weights: init 2");
        for (i = 0; i < IMW * IMH; i++)
            depth[i] = (i % 2) ? -1.0f : (float)NAN;
        CHECK(mv_tsdf_fuse_depthmap(&t2, &tab, depth, IMW, IMH, 0.0)
              == MV_OK, "weights: fuse invalid map");
        nvox = (long)t2.nx * t2.ny * t2.nz;
        for (i = 0; i < nvox; i++)
            nobs += t2.w[i] > 0.0f;
        CHECK(nobs == 0, "weights: depth <= 0 and NaN update nothing");
        mv_tsdf_free(&t2);
    }
    CHECK(mv_tsdf_fuse_depthmap(&t, &tab, depth, IMW / 2, IMH, 0.0)
          == MV_ERR, "weights: size mismatch rejected");
    mv_tsdf_table_free(&tab);
    mv_tsdf_free(&t);
}

/* ---- (d): timing, sweep vs fold on a 129^3 grid --------------------- */

static void test_timing(void)
{
    mv_camera cam;
    mv_tsdf ts, tf;
    mv_tsdf_table tab;
    static float depth[IMW * IMH];
    const double ctr[3] = { 0.0, 0.0, 2.0 };
    const double r = 0.4, voxel = 0.01, tau = 0.03;
    double t_sweep, t_fold, ratio;
    long nvox;
    clock_t c0;
    int rep;

    make_cam(&cam);
    CHECK(mv_tsdf_init(&ts, -0.64, -0.64, 1.36, 1.28, 1.28, 1.28,
                       voxel, tau) == MV_OK, "timing: sweep init");
    CHECK(mv_tsdf_init(&tf, -0.64, -0.64, 1.36, 1.28, 1.28, 1.28,
                       voxel, tau) == MV_OK, "timing: fold init");
    CHECK(mv_tsdf_table_build(&ts, &cam, IMW, IMH, &tab) == MV_OK,
          "timing: table build");
    nvox = (long)ts.nx * ts.ny * ts.nz;
    depth_sphere(depth, &cam, ctr, r);

    c0 = clock();
    for (rep = 0; rep < 8; rep++)
        mv_tsdf_fuse_depthmap(&ts, &tab, depth, IMW, IMH, 0.0);
    t_sweep = (double)(clock() - c0) / CLOCKS_PER_SEC / 8.0;

    c0 = clock();
    for (rep = 0; rep < 2; rep++)
        fold_depthmap(&tf, &cam, depth, 1.0);
    t_fold = (double)(clock() - c0) / CLOCKS_PER_SEC / 2.0;

    ratio = t_sweep > 0.0 ? t_fold / t_sweep : 0.0;
    printf("timing: %ld voxels; sweep %.1f ms/frame (%.0fM voxel "
           "visits/s), fold %.1f ms/frame; fold/sweep = %.1fx\n",
           nvox, t_sweep * 1e3, nvox / t_sweep / 1e6, t_fold * 1e3,
           ratio);
    /* correctness is asserted elsewhere; here only that both ran */
    CHECK(t_sweep > 0.0 && t_fold > 0.0 && ratio > 0.0,
          "timing: both paths measured");
    mv_tsdf_table_free(&tab);
    mv_tsdf_free(&ts);
    mv_tsdf_free(&tf);
}

int main(void)
{
    test_plane_scene();
    test_sphere_scene();
    test_weights();
    test_timing();
    if (failures) {
        printf("\n%d FAILURE(S)\n", failures);
        return 1;
    }
    printf("\nall tests passed\n");
    return 0;
}
