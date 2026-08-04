/* Room-model and change-detection simulation: numbers quoted in
 * doc/multiview.tex ("Basic results").
 *
 * A textured synthetic room (floor, back wall, side wall, and a box) is
 * imaged by the standard rig through the REAL dense pipeline: rendered
 * stereo pairs -> rectification warp -> SAD block matching -> depth ->
 * per-pixel temporal median background -> enriched samples -> TSDF ->
 * marching-tetrahedra mesh. Then the box moves and a single new frame is
 * compared against the background depth for change detection. Scored
 * against the renderer's ground-truth depth and the true plane set.
 *
 * Run: make demo_room && ./demo_room     (~20 s)
 * Out: out_room_left.pgm   raw left view          (Preview opens PGM)
 *      out_room_rect.pgm   rectified left view
 *      out_room_disp.pgm   disparity map
 *      out_room_truth.pgm  shaded view of the TRUE generating model
 *      out_room_mesh.pgm   shaded view of the reconstructed mesh
 *      out_room_change.pgm change mask after the box moves
 *      out_room.ply        the room mesh itself */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mv/mv.h"

#define MV_PI 3.14159265358979323846
#define W 640
#define H 480
#define NBG 9              /* background frames */
#define SIGMA_I 2.0        /* image noise, gray levels */
#define DMIN 58
#define DMAX 130
#define HWIN 3
#define VOXEL 0.015
#define TAU 0.06
#define FLOOR_Y 0.8
#define WALL_Z 6.2
#define WALL_X (-1.4)

static unsigned long long rng = 20260804;

static double urand(void)
{
    rng = rng * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)((rng >> 11) & ((1ULL << 53) - 1)) / 9007199254740992.0;
}

/* speckle texture with strong local contrast for the matcher */
static void make_tex(unsigned char *t, int w, int h)
{
    int i, x, y;
    unsigned char *tmp = (unsigned char *)malloc((size_t)w * h);
    for (i = 0; i < w * h; i++)
        t[i] = (unsigned char)(40 + 175 * (urand() > 0.5));
    /* one box blur to soften aliasing */
    memcpy(tmp, t, (size_t)w * h);
    for (y = 1; y < h - 1; y++)
        for (x = 1; x < w - 1; x++) {
            int s = 0, a, b;
            for (b = -1; b <= 1; b++)
                for (a = -1; a <= 1; a++)
                    s += tmp[(y + b) * w + (x + a)];
            t[y * w + x] = (unsigned char)(s / 9);
        }
    free(tmp);
}

static void set_plane(mv_plane *p, const unsigned char *tex, int tw, int th,
                      double pitch, const double r1[3], const double r2[3],
                      const double t[3])
{
    double r3[3];
    int i;
    mv_cross3(r3, r1, r2);
    for (i = 0; i < 3; i++) {
        p->R[i * 3 + 0] = r1[i];
        p->R[i * 3 + 1] = r2[i];
        p->R[i * 3 + 2] = r3[i];
        p->t[i] = t[i];
    }
    p->tex = tex;
    p->tw = tw;
    p->th = th;
    p->pitch = pitch;
}

/* box: 3 visible faces (front, top, camera-side), 0.5 x 0.45 x 0.4 m */
#define BOXW 0.5
#define BOXH 0.45
#define BOXD 0.4
static void add_box(mv_plane *pls, int *n, const unsigned char *tex,
                    double bx, double bz)
{
    static const double ex[3] = { 1, 0, 0 }, ey[3] = { 0, 1, 0 },
                        ez[3] = { 0, 0, 1 };
    double t[3];
    t[0] = bx - BOXW / 2; t[1] = FLOOR_Y - BOXH; t[2] = bz;
    set_plane(&pls[(*n)++], tex, (int)(BOXW / 0.02), (int)(BOXH / 0.02),
              0.02, ex, ey, t); /* front */
    set_plane(&pls[(*n)++], tex, (int)(BOXW / 0.02), (int)(BOXD / 0.02),
              0.02, ex, ez, t); /* top */
    t[0] = bx - BOXW / 2;
    set_plane(&pls[(*n)++], tex, (int)(BOXD / 0.02), (int)(BOXH / 0.02),
              0.02, ez, ey, t); /* left side */
}

/* distance of a point to the true phase-1 scene (planes with extents) */
static double scene_dist(const double p[3], double bx, double bz)
{
    double d = HUGE_VAL, c;
    /* floor, back wall, side wall (extent checks loose) */
    c = fabs(p[1] - FLOOR_Y);
    if (c < d) d = c;
    c = fabs(p[2] - WALL_Z);
    if (c < d) d = c;
    c = fabs(p[0] - WALL_X);
    if (c < d) d = c;
    /* box faces */
    if (p[2] > bz - 0.1 && p[2] < bz + BOXD + 0.1
        && p[0] > bx - BOXW / 2 - 0.1 && p[0] < bx + BOXW / 2 + 0.1) {
        c = fabs(p[1] - (FLOOR_Y - BOXH));
        if (c < d) d = c;
        c = fabs(p[2] - bz);
        if (c < d) d = c;
        c = fabs(p[0] - (bx - BOXW / 2));
        if (c < d) d = c;
    }
    return d;
}

/* invalidate disparities where the left window has too little texture
 * (the doc's "texture desert" failure: SAD returns garbage there) */
static void texture_gate(float *disp, const unsigned char *img)
{
    int x, y, a, b;
    for (y = HWIN; y < H - HWIN; y++)
        for (x = HWIN; x < W - HWIN; x++) {
            double s = 0.0, s2 = 0.0;
            int n = 0;
            for (b = -HWIN; b <= HWIN; b++)
                for (a = -HWIN; a <= HWIN; a++) {
                    double v = img[(y + b) * W + (x + a)];
                    s += v;
                    s2 += v * v;
                    n++;
                }
            if (s2 / n - (s / n) * (s / n) < 36.0)
                disp[y * W + x] = -1.0f;
        }
}

static void write_norm_pgm(const char *path, const float *v, int w, int h,
                           float lo, float hi)
{
    unsigned char *img = (unsigned char *)malloc((size_t)w * h);
    int i;
    for (i = 0; i < w * h; i++) {
        float x = v[i];
        if (x >= HUGE_VALF || x < 0.0f)
            img[i] = 0;
        else {
            float g = (x - lo) / (hi - lo) * 255.0f;
            img[i] = (unsigned char)(g < 0 ? 0 : g > 255 ? 255 : g);
        }
    }
    mv_pgm_write(path, img, w, h);
    free(img);
}

/* ---- shaded software preview of a triangle soup ---------------------- */
static void mesh_preview(const char *path, const double *tris, int ntri,
                         const double *alb)
{
    mv_camera pc;
    static unsigned char img[W * H];
    static float zb[W * H];
    double Cpos[3] = { 1.5, -1.0, 2.0 }; /* above-right, looking in */
    double light[3] = { -0.45, 0.65, 0.6 };
    int i, k;

    mv_normalize(light, 3);
    mv_cam_set_K(&pc, 700.0, 700.0, 320.0, 240.0);
    mv_cam_set_pose_yaw(&pc, -30.0 * MV_PI / 180.0, Cpos);
    /* add a downward pitch by composing a small rotation about x */
    {
        double a = 18.0 * MV_PI / 180.0, Rx[9], Rn[9];
        Rx[0] = 1; Rx[1] = 0; Rx[2] = 0;
        Rx[3] = 0; Rx[4] = cos(a); Rx[5] = -sin(a);
        Rx[6] = 0; Rx[7] = sin(a); Rx[8] = cos(a);
        mv_mat_mul(Rn, Rx, pc.R, 3, 3, 3);
        memcpy(pc.R, Rn, sizeof(Rn));
        for (i = 0; i < 3; i++)
            pc.t[i] = -(Rn[i * 3 + 0] * Cpos[0] + Rn[i * 3 + 1] * Cpos[1]
                        + Rn[i * 3 + 2] * Cpos[2]);
    }
    memset(img, 15, sizeof(img));
    for (i = 0; i < W * H; i++)
        zb[i] = HUGE_VALF;

    for (k = 0; k < ntri; k++) {
        const double *a = tris + 9 * k, *b = a + 3, *c = a + 6;
        double ua[2], ub[2], uc[2], n[3], e1[3], e2[3], shade;
        double za, zbv, zc2;
        int xmin, xmax, ymin, ymax, x, y, j;
        if (mv_cam_project(ua, &pc, a) != MV_OK
            || mv_cam_project(ub, &pc, b) != MV_OK
            || mv_cam_project(uc, &pc, c) != MV_OK)
            continue;
        for (j = 0; j < 3; j++) {
            e1[j] = b[j] - a[j];
            e2[j] = c[j] - a[j];
        }
        mv_cross3(n, e1, e2);
        if (mv_normalize(n, 3) != MV_OK)
            continue;
        shade = (45.0 + 190.0 * fabs(mv_dot(n, light, 3)))
              * (alb ? alb[k] : 1.0);
        za = pc.R[6] * a[0] + pc.R[7] * a[1] + pc.R[8] * a[2] + pc.t[2];
        zbv = pc.R[6] * b[0] + pc.R[7] * b[1] + pc.R[8] * b[2] + pc.t[2];
        zc2 = pc.R[6] * c[0] + pc.R[7] * c[1] + pc.R[8] * c[2] + pc.t[2];
        xmin = (int)fmin(ua[0], fmin(ub[0], uc[0]));
        xmax = (int)ceil(fmax(ua[0], fmax(ub[0], uc[0])));
        ymin = (int)fmin(ua[1], fmin(ub[1], uc[1]));
        ymax = (int)ceil(fmax(ua[1], fmax(ub[1], uc[1])));
        if (xmin < 0) xmin = 0;
        if (ymin < 0) ymin = 0;
        if (xmax >= W) xmax = W - 1;
        if (ymax >= H) ymax = H - 1;
        {
            double d = (ub[0] - ua[0]) * (uc[1] - ua[1])
                     - (uc[0] - ua[0]) * (ub[1] - ua[1]);
            if (fabs(d) < 1e-9)
                continue;
            for (y = ymin; y <= ymax; y++)
                for (x = xmin; x <= xmax; x++) {
                    double l2 = ((x - ua[0]) * (uc[1] - ua[1])
                                 - (uc[0] - ua[0]) * (y - ua[1])) / d;
                    double l3 = ((ub[0] - ua[0]) * (y - ua[1])
                                 - (x - ua[0]) * (ub[1] - ua[1])) / d;
                    double l1 = 1.0 - l2 - l3, z;
                    if (l1 < 0 || l2 < 0 || l3 < 0)
                        continue;
                    z = l1 * za + l2 * zbv + l3 * zc2;
                    if (z > 0 && z < zb[y * W + x]) {
                        double att = 3.0 / z; /* mild depth cue */
                        double v2;
                        if (att < 0.55) att = 0.55;
                        if (att > 1.15) att = 1.15;
                        v2 = shade * att;
                        zb[y * W + x] = (float)z;
                        img[y * W + x] =
                            (unsigned char)(v2 > 255 ? 255 : v2);
                    }
                }
        }
    }
    mv_pgm_write(path, img, W, H);
}

int main(void)
{
    static unsigned char texA[140 * 160], texB[160 * 80], texC[160 * 80],
        texD[40 * 40];
    static unsigned char imgL[W * H], imgR[W * H], rectL[W * H],
        rectR[W * H];
    static float disp[W * H], depth[W * H], bg[W * H], zTrue[W * H],
        zTrue2[W * H];
    static float stack[NBG][W * H];
    mv_camera c1, c2, r1, r2;
    double C2pos[3] = { 0.5, 0.0, 0.0 }, H1[9], H2[9], Kinv[9], C1[3];
    mv_plane pls[8];
    mv_tsdf t;
    double bx1 = 0.35, bz1 = 4.1, bx2 = -0.55, bz2 = 4.7;
    double fnew, Bnew;
    int npl = 0, f, i, x, y;

    printf("multiview room-model experiment (seed 20260804)\n");
    printf("-----------------------------------------------\n");

    make_tex(texA, 140, 160);
    make_tex(texB, 160, 80);
    make_tex(texC, 160, 80);
    make_tex(texD, 40, 40);
    {
        static const double ex[3] = { 1, 0, 0 }, ey[3] = { 0, 1, 0 },
                            ez[3] = { 0, 0, 1 };
        double t0[3];
        t0[0] = -1.4; t0[1] = FLOOR_Y; t0[2] = 3.0;
        set_plane(&pls[npl++], texA, 140, 160, 0.02, ex, ez, t0); /* floor */
        t0[0] = -1.6; t0[1] = -0.8; t0[2] = WALL_Z;
        set_plane(&pls[npl++], texB, 160, 80, 0.02, ex, ey, t0); /* back */
        t0[0] = WALL_X; t0[1] = -0.8; t0[2] = 3.0;
        set_plane(&pls[npl++], texC, 160, 80, 0.02, ez, ey, t0); /* side */
    }
    add_box(pls, &npl, texD, bx1, bz1);

    mv_cam_set_K(&c1, 800.0, 800.0, 320.0, 240.0);
    mv_cam_set_identity_pose(&c1);
    memset(c1.k, 0, sizeof(c1.k));
    c2 = c1;
    mv_cam_set_pose_yaw(&c2, -6.0 * MV_PI / 180.0, C2pos);

    if (mv_rectify_pair(&c1, &c2, &r1, &r2, H1, H2) != MV_OK)
        return 1;
    fnew = r1.K[0];
    Bnew = mv_baseline(&r1, &r2);
    mv_mat_inv3(Kinv, r1.K);
    mv_cam_center(C1, &r1);

    /* ---- phase 1: background frames through the dense pipeline ----- */
    for (f = 0; f < NBG; f++) {
        mv_render_scene(imgL, f == 0 ? zTrue : NULL, W, H, &c1, pls, npl,
                        10, SIGMA_I, &rng);
        mv_render_scene(imgR, NULL, W, H, &c2, pls, npl, 10, SIGMA_I,
                        &rng);
        if (f == 0)
            mv_pgm_write("out_room_left.pgm", imgL, W, H);
        mv_warp_homography(rectL, W, H, imgL, W, H, H1, 0);
        mv_warp_homography(rectR, W, H, imgR, W, H, H2, 0);
        if (f == 0)
            mv_pgm_write("out_room_rect.pgm", rectL, W, H);
        mv_stereo_sad(disp, rectL, rectR, W, H, HWIN, DMIN, DMAX);
        texture_gate(disp, rectL);
        for (i = 0; i < W * H; i++)
            stack[f][i] = disp[i] > 0.0f
                          ? (float)(fnew * Bnew / disp[i]) : -1.0f;
        if (f == 0)
            write_norm_pgm("out_room_disp.pgm", disp, W, H, DMIN, DMAX);
    }
    /* per-pixel median of valid depths */
    {
        int nvalid_px = 0;
        for (i = 0; i < W * H; i++) {
            float v[NBG], tmp;
            int n = 0, a, b;
            for (f = 0; f < NBG; f++)
                if (stack[f][i] > 0.0f)
                    v[n++] = stack[f][i];
            if (n < 5) {
                bg[i] = -1.0f;
                continue;
            }
            for (a = 1; a < n; a++)
                for (b = a; b > 0 && v[b] < v[b - 1]; b--) {
                    tmp = v[b];
                    v[b] = v[b - 1];
                    v[b - 1] = tmp;
                }
            bg[i] = v[n / 2];
            nvalid_px++;
        }
        printf("background coverage : %.1f%% of pixels\n",
               100.0 * nvalid_px / (W * H));
    }

    /* ---- shaded view of the TRUE generating model ------------------- */
    {
        double gt[8 * 2 * 9]; /* two triangles per plane */
        double alb[16];
        /* distinct albedos so a human can tell the surfaces apart:
         * floor light, back wall mid, side wall dark, box brightest */
        static const double PALB[6] = { 1.0, 0.72, 0.5, 1.45, 1.3, 1.15 };
        int ngt = 0, p, j;
        for (p = 0; p < npl; p++) {
            const mv_plane *pl = &pls[p];
            double c00[3], c10[3], c01[3], c11[3];
            for (j = 0; j < 3; j++) {
                c00[j] = pl->t[j];
                c10[j] = pl->t[j] + pl->R[j * 3 + 0] * pl->tw * pl->pitch;
                c01[j] = pl->t[j] + pl->R[j * 3 + 1] * pl->th * pl->pitch;
                c11[j] = c10[j] + pl->R[j * 3 + 1] * pl->th * pl->pitch;
            }
            memcpy(gt + 9 * ngt + 0, c00, 24);
            memcpy(gt + 9 * ngt + 3, c10, 24);
            memcpy(gt + 9 * ngt + 6, c11, 24);
            alb[ngt] = PALB[p];
            ngt++;
            memcpy(gt + 9 * ngt + 0, c00, 24);
            memcpy(gt + 9 * ngt + 3, c11, 24);
            memcpy(gt + 9 * ngt + 6, c01, 24);
            alb[ngt] = PALB[p];
            ngt++;
        }
        mesh_preview("out_room_truth.pgm", gt, ngt, alb);
    }

    /* ---- fuse background into TSDF, extract, preview, score -------- */
    if (mv_tsdf_init(&t, -1.5, -0.7, 3.1, 3.0, 1.6, 3.2, VOXEL, TAU)
        != MV_OK)
        return 1;
    for (y = 0; y < H; y++)
        for (x = 0; x < W; x++) {
            float Z = bg[y * W + x];
            double pc[3], pw[3], sigZ, wgt;
            int j;
            if (Z <= 0.0f)
                continue;
            pc[0] = (Kinv[0] * x + Kinv[1] * y + Kinv[2]) * Z;
            pc[1] = (Kinv[4] * y + Kinv[5]) * Z;
            pc[2] = Z;
            for (j = 0; j < 3; j++)
                pw[j] = r1.R[0 * 3 + j] * pc[0] + r1.R[1 * 3 + j] * pc[1]
                      + r1.R[2 * 3 + j] * pc[2] + C1[j];
            sigZ = (double)Z * Z * 0.35 / (fnew * Bnew);
            wgt = 1.0 / (sigZ * sigZ);
            mv_tsdf_fuse(&t, pw, C1, wgt);
        }
    {
        double *tris;
        int ntri;
        double se = 0.0;
        long ns = 0;
        if (mv_tsdf_mesh(&t, &tris, &ntri) != MV_OK)
            return 1;
        for (i = 0; i < 3 * ntri; i++) {
            double e = scene_dist(tris + 3 * i, bx1, bz1);
            if (e < 0.5) { /* exclude fringe artifacts from the RMS */
                se += e * e;
                ns++;
            }
        }
        printf("room mesh           : %d triangles, surface RMS %.1f mm "
               "(%.1f%% of vertices scored)\n", ntri,
               1000.0 * sqrt(se / ns), 100.0 * ns / (3.0 * ntri));
        mv_tsdf_write_ply(&t, "out_room.ply");
        mesh_preview("out_room_mesh.pgm", tris, ntri, NULL);
        free(tris);
    }
    mv_tsdf_free(&t);

    /* ---- phase 2: move the box, detect change ----------------------- */
    npl = 3;
    add_box(pls, &npl, texD, bx2, bz2);
    mv_render_scene(imgL, zTrue2, W, H, &c1, pls, npl, 10, SIGMA_I, &rng);
    mv_render_scene(imgR, NULL, W, H, &c2, pls, npl, 10, SIGMA_I, &rng);
    mv_warp_homography(rectL, W, H, imgL, W, H, H1, 0);
    mv_warp_homography(rectR, W, H, imgR, W, H, H2, 0);
    mv_stereo_sad(disp, rectL, rectR, W, H, HWIN, DMIN, DMAX);
    texture_gate(disp, rectL);
    for (i = 0; i < W * H; i++)
        depth[i] = disp[i] > 0.0f ? (float)(fnew * Bnew / disp[i]) : -1.0f;
    {
        static unsigned char mask[W * H], truemask[W * H];
        long tp = 0, fp = 0, fn = 0, truechg = 0;
        for (y = 0; y < H; y++)
            for (x = 0; x < W; x++) {
                double Zb = bg[y * W + x], Zn = depth[y * W + x];
                double sigZ, thr;
                int est = 0, tru;
                i = y * W + x;
                /* ground-truth change in the ORIGINAL left image, mapped
                 * through H1 for comparison in rectified coords */
                {
                    double Hi[9], X, Y, Wd;
                    mv_mat_inv3(Hi, H1);
                    X = Hi[0] * x + Hi[1] * y + Hi[2];
                    Y = Hi[3] * x + Hi[4] * y + Hi[5];
                    Wd = Hi[6] * x + Hi[7] * y + Hi[8];
                    tru = 0;
                    if (fabs(Wd) > 1e-12) {
                        int ox = (int)(X / Wd + 0.5),
                            oy = (int)(Y / Wd + 0.5);
                        if (ox >= 0 && oy >= 0 && ox < W && oy < H
                            && zTrue[oy * W + ox] < HUGE_VALF
                            && zTrue2[oy * W + ox] < HUGE_VALF)
                            tru = fabs(zTrue[oy * W + ox]
                                       - zTrue2[oy * W + ox]) > 0.05;
                    }
                }
                if (Zb > 0.0 && Zn > 0.0) {
                    sigZ = Zn * Zn * 0.35 / (fnew * Bnew);
                    thr = 4.0 * sigZ > 0.08 ? 4.0 * sigZ : 0.08;
                    est = fabs(Zn - Zb) > thr;
                }
                mask[i] = est ? 255 : 0;
                truemask[i] = (unsigned char)tru;
                if (tru) {
                    truechg++;
                    if (est) tp++; else fn++;
                } else if (est)
                    fp++;
            }
        /* 3x3 majority clean: isolated flickers are not changes */
        {
            static unsigned char m2[W * H];
            memcpy(m2, mask, sizeof(m2));
            for (y = 1; y < H - 1; y++)
                for (x = 1; x < W - 1; x++) {
                    int a, b, n = 0;
                    for (b = -1; b <= 1; b++)
                        for (a = -1; a <= 1; a++)
                            n += m2[(y + b) * W + (x + a)] > 0;
                    mask[y * W + x] = (n >= 5) ? 255 : 0;
                }
            tp = fp = fn = 0;
            for (i = 0; i < W * H; i++) {
                int est = mask[i] > 0, tru = truemask[i];
                if (tru) {
                    if (est) tp++; else fn++;
                } else if (est)
                    fp++;
            }
        }
        printf("change detection    : %ld true changed px; recall %.1f%%, "
               "false alarms %.2f%% of image\n", truechg,
               100.0 * tp / (tp + fn), 100.0 * fp / (W * H));
        mv_pgm_write("out_room_change.pgm", mask, W, H);
    }
    printf("\nwrote out_room_left/disp/mesh/change.pgm and out_room.ply\n");
    return 0;
}
