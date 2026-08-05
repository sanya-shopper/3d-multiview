/* Scene point cloud from a single handheld video of the calibration
 * pattern: the pattern gives every decoded frame a metric pose (162
 * known screen points -> homography -> R, t under the calibrated K), so
 * frame pairs with a few centimetres of accidental baseline become
 * calibrated stereo pairs.  One phone waved in front of a laptop is a
 * many-camera rig on a time budget.
 *
 * Usage: scenecloud <pitch_mm> <fx> <fy> <cx> <cy> <k1> <k2>
 *                   <out_prefix> <frames.pgm...>
 *   intrinsics from calibreal on the same video.  Writes
 *   <out_prefix>.ply (gray-textured points, metric, screen plane = z=0)
 *   and <out_prefix>_view{0,1,2}.pgm splat renders.
 *
 * Build: make scenecloud */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mv/mv.h"

#define MV_PI 3.14159265358979323846
#define MAXF 220
#define MAXPAIR 10
#define ZMIN 0.25 /* m: accept triangulated depths in this range */
#define ZMAX 2.5

typedef struct {
    mv_camera cam;      /* pinhole at HALF resolution, no distortion */
    unsigned char *img; /* undistorted half-res gray */
    double C[3];        /* camera center */
    int w, h;
} frame_t;

/* pose of the screen plane from decoded corners: Zhang extrinsics from
 * the homography, with the rotation snapped to the nearest orthonormal
 * matrix via SVD (R = U V^T) */
static int pose_from_read(mv_camera *cam, const double K[9],
                          const mv_read_result *rr, double pitch)
{
    double obj[2 * MV_READ_MAXC], H[9], Kinv[9];
    double r1[3], r2[3], r3[3], t[3], l1 = 0, l2 = 0, lam;
    double A[9], S[3], V[9], U[9], Vt[9];
    int i, j;

    for (i = 0; i < rr->n; i++) {
        double xy[2];
        mv_pattern_corner_px(rr->id[i] % MV_PAT_CORNER_COLS,
                             rr->id[i] / MV_PAT_CORNER_COLS, xy);
        obj[2 * i] = xy[0] * pitch;
        obj[2 * i + 1] = xy[1] * pitch;
    }
    if (mv_homography_dlt(H, obj, rr->uv, rr->n) != MV_OK)
        return MV_ERR;
    if (mv_mat_inv3(Kinv, K) != MV_OK)
        return MV_ERR;
    for (i = 0; i < 3; i++) {
        r1[i] = Kinv[3 * i] * H[0] + Kinv[3 * i + 1] * H[3]
              + Kinv[3 * i + 2] * H[6];
        r2[i] = Kinv[3 * i] * H[1] + Kinv[3 * i + 1] * H[4]
              + Kinv[3 * i + 2] * H[7];
        t[i] = Kinv[3 * i] * H[2] + Kinv[3 * i + 1] * H[5]
             + Kinv[3 * i + 2] * H[8];
        l1 += r1[i] * r1[i];
        l2 += r2[i] * r2[i];
    }
    lam = 2.0 / (sqrt(l1) + sqrt(l2));
    if (t[2] * lam < 0.0)
        lam = -lam; /* plane must be in front of the camera */
    for (i = 0; i < 3; i++) {
        r1[i] *= lam;
        r2[i] *= lam;
        t[i] *= lam;
    }
    r3[0] = r1[1] * r2[2] - r1[2] * r2[1];
    r3[1] = r1[2] * r2[0] - r1[0] * r2[2];
    r3[2] = r1[0] * r2[1] - r1[1] * r2[0];
    for (i = 0; i < 3; i++) {
        A[3 * i + 0] = r1[i];
        A[3 * i + 1] = r2[i];
        A[3 * i + 2] = r3[i];
    }
    memcpy(U, A, sizeof(U));
    if (mv_svd(U, S, V, 3, 3) != MV_OK)
        return MV_ERR;
    mv_mat_transpose(Vt, V, 3, 3);
    mv_mat_mul(cam->R, U, Vt, 3, 3, 3);
    memcpy(cam->K, K, 9 * sizeof(double));
    for (j = 0; j < 3; j++)
        cam->t[j] = t[j];
    memset(cam->k, 0, sizeof(cam->k));
    return MV_OK;
}

/* 2x2 box downsample, then remove radial distortion by inverse mapping
 * (output is an ideal pinhole image under Khalf) */
static unsigned char *half_undistort(const unsigned char *img, int w,
                                     int h, const double K[9],
                                     const double kr[2], double Khalf[9])
{
    int hw = w / 2, hh = h / 2, x, y;
    unsigned char *tmp, *out;
    double Kinv[9];

    tmp = (unsigned char *)malloc((size_t)hw * hh);
    out = (unsigned char *)malloc((size_t)hw * hh);
    if (!tmp || !out) {
        free(tmp);
        free(out);
        return NULL;
    }
    for (y = 0; y < hh; y++)
        for (x = 0; x < hw; x++)
            tmp[y * hw + x] = (unsigned char)((img[2 * y * w + 2 * x]
                + img[2 * y * w + 2 * x + 1]
                + img[(2 * y + 1) * w + 2 * x]
                + img[(2 * y + 1) * w + 2 * x + 1] + 2) / 4);
    memset(Khalf, 0, 9 * sizeof(double));
    Khalf[0] = K[0] / 2.0;
    Khalf[2] = (K[2] + 0.5) / 2.0 - 0.5;
    Khalf[4] = K[4] / 2.0;
    Khalf[5] = (K[5] + 0.5) / 2.0 - 0.5;
    Khalf[8] = 1.0;
    if (mv_mat_inv3(Kinv, Khalf) != MV_OK) {
        free(tmp);
        free(out);
        return NULL;
    }
    for (y = 0; y < hh; y++)
        for (x = 0; x < hw; x++) {
            double xn = Kinv[0] * x + Kinv[1] * y + Kinv[2];
            double yn = Kinv[4] * y + Kinv[5];
            double r2 = xn * xn + yn * yn;
            double rad = 1.0 + r2 * (kr[0] + r2 * kr[1]);
            double u = Khalf[0] * xn * rad + Khalf[2];
            double v = Khalf[4] * yn * rad + Khalf[5];
            int u0 = (int)floor(u), v0 = (int)floor(v);
            double fu = u - u0, fv = v - v0;
            unsigned char g = 0;
            if (u0 >= 0 && v0 >= 0 && u0 + 1 < hw && v0 + 1 < hh) {
                double g00 = tmp[v0 * hw + u0], g10 = tmp[v0 * hw + u0 + 1];
                double g01 = tmp[(v0 + 1) * hw + u0],
                       g11 = tmp[(v0 + 1) * hw + u0 + 1];
                g = (unsigned char)((1 - fv) * ((1 - fu) * g00 + fu * g10)
                                    + fv * ((1 - fu) * g01 + fu * g11)
                                    + 0.5);
            }
            out[y * hw + x] = g;
        }
    free(tmp);
    return out;
}

/* reject disparities that disagree with the median of their valid 3x3
 * neighbours by more than 2 px (SAD speckle) */
static void speckle_gate(float *disp, int w, int h)
{
    float *cp;
    int x, y, a, b;
    cp = (float *)malloc((size_t)w * h * sizeof(float));
    if (!cp)
        return;
    memcpy(cp, disp, (size_t)w * h * sizeof(float));
    for (y = 1; y < h - 1; y++)
        for (x = 1; x < w - 1; x++) {
            float v[9], t;
            int n = 0, i, j2;
            if (cp[y * w + x] <= 0.0f)
                continue;
            for (b = -1; b <= 1; b++)
                for (a = -1; a <= 1; a++) {
                    t = cp[(y + b) * w + x + a];
                    if (t > 0.0f)
                        v[n++] = t;
                }
            if (n < 5) {
                disp[y * w + x] = -1.0f;
                continue;
            }
            for (i = 1; i < n; i++)
                for (j2 = i; j2 > 0 && v[j2] < v[j2 - 1]; j2--) {
                    t = v[j2];
                    v[j2] = v[j2 - 1];
                    v[j2 - 1] = t;
                }
            if (fabsf(cp[y * w + x] - v[n / 2]) > 2.0f)
                disp[y * w + x] = -1.0f;
        }
    free(cp);
}

static int cmp_dbl(const void *a, const void *b)
{
    double d = *(const double *)a - *(const double *)b;
    return d < 0 ? -1 : d > 0 ? 1 : 0;
}

/* invalidate disparities where the matching window is textureless */
static void texture_gate(float *disp, const unsigned char *img, int w,
                         int h)
{
    int x, y, a, b;
    for (y = 4; y < h - 4; y++)
        for (x = 4; x < w - 4; x++) {
            double s = 0.0, ss = 0.0, v;
            for (b = -4; b <= 4; b++)
                for (a = -4; a <= 4; a++) {
                    v = img[(y + b) * w + x + a];
                    s += v;
                    ss += v * v;
                }
            v = ss / 81.0 - (s / 81.0) * (s / 81.0);
            if (v < 40.0)
                disp[y * w + x] = -1.0f;
        }
}

/* splat render of the cloud from a look-at camera; gray from the cloud,
 * simple z-buffer, 2x2 splats */
static void splat_view(const char *path, const mv_cloud *cl,
                       const double eye[3], const double at[3], double f,
                       int w, int h)
{
    mv_camera cam;
    unsigned char *img;
    float *zb;
    double fw[3], rt[3], up0[3] = { 0.0, 1.0, 0.0 }, up[3], n;
    int i, j, x, y;

    img = (unsigned char *)malloc((size_t)w * h);
    zb = (float *)malloc((size_t)w * h * sizeof(float));
    if (!img || !zb) {
        free(img);
        free(zb);
        return;
    }
    memset(img, 16, (size_t)w * h);
    for (i = 0; i < w * h; i++)
        zb[i] = HUGE_VALF;
    for (j = 0; j < 3; j++)
        fw[j] = at[j] - eye[j];
    n = sqrt(fw[0] * fw[0] + fw[1] * fw[1] + fw[2] * fw[2]);
    if (!(n > 1e-9)) { /* degenerate framing (eye == center): skip */
        free(img);
        free(zb);
        return;
    }
    for (j = 0; j < 3; j++)
        fw[j] /= n;
    rt[0] = up0[1] * fw[2] - up0[2] * fw[1];
    rt[1] = up0[2] * fw[0] - up0[0] * fw[2];
    rt[2] = up0[0] * fw[1] - up0[1] * fw[0];
    n = sqrt(rt[0] * rt[0] + rt[1] * rt[1] + rt[2] * rt[2]);
    for (j = 0; j < 3; j++)
        rt[j] /= n;
    up[0] = fw[1] * rt[2] - fw[2] * rt[1];
    up[1] = fw[2] * rt[0] - fw[0] * rt[2];
    up[2] = fw[0] * rt[1] - fw[1] * rt[0];
    mv_cam_set_K(&cam, f, f, w / 2.0, h / 2.0);
    for (j = 0; j < 3; j++) {
        cam.R[0 * 3 + j] = rt[j];
        cam.R[1 * 3 + j] = up[j];
        cam.R[2 * 3 + j] = fw[j];
    }
    for (j = 0; j < 3; j++)
        cam.t[j] = -(cam.R[j * 3 + 0] * eye[0] + cam.R[j * 3 + 1] * eye[1]
                     + cam.R[j * 3 + 2] * eye[2]);
    memset(cam.k, 0, sizeof(cam.k));
    for (i = 0; i < cl->n; i++) {
        double X[3], uv[2], Z;
        X[0] = cl->xyz[3 * i];
        X[1] = cl->xyz[3 * i + 1];
        X[2] = cl->xyz[3 * i + 2];
        if (mv_cam_project(uv, &cam, X) != MV_OK)
            continue;
        Z = cam.R[6] * X[0] + cam.R[7] * X[1] + cam.R[8] * X[2] + cam.t[2];
        for (y = (int)uv[1]; y <= (int)uv[1] + 1; y++)
            for (x = (int)uv[0]; x <= (int)uv[0] + 1; x++) {
                if (x < 0 || y < 0 || x >= w || y >= h)
                    continue;
                if (Z < zb[y * w + x]) {
                    zb[y * w + x] = (float)Z;
                    img[y * w + x] = cl->rgb ? cl->rgb[3 * i] : 200;
                }
            }
    }
    mv_pgm_write(path, img, w, h);
    free(img);
    free(zb);
}

int main(int argc, char **argv)
{
    static frame_t fr[MAXF];
    double K[9], kr[2], pitch;
    const char *prefix;
    mv_cloud cloud;
    unsigned char *tag = NULL;
    int ntag = 0, captag = 0;
    int nf = 0, f, i, npairs = 0, lastused = -1;

    if (argc < 10) {
        fprintf(stderr, "usage: %s <pitch_mm> <fx> <fy> <cx> <cy> <k1> "
                        "<k2> <out_prefix> <frames.pgm...>\n", argv[0]);
        return 1;
    }
    pitch = atof(argv[1]) * 1e-3;
    memset(K, 0, sizeof(K));
    K[0] = atof(argv[2]);
    K[4] = atof(argv[3]);
    K[2] = atof(argv[4]);
    K[5] = atof(argv[5]);
    K[8] = 1.0;
    kr[0] = atof(argv[6]);
    kr[1] = atof(argv[7]);
    prefix = argv[8];

    for (f = 9; f < argc && nf < MAXF; f++) {
        unsigned char *img;
        int w, h;
        mv_read_result rr;
        mv_camera cam;
        double Khalf[9];
        if (mv_pgm_read(argv[f], &img, &w, &h) != MV_OK)
            continue;
        if (mv_read_pattern(&rr, img, w, h) != MV_OK || rr.n < 100) {
            free(img);
            continue;
        }
        /* pose under the FULL-res K (corners are full-res pixels) */
        if (pose_from_read(&cam, K, &rr, pitch) != MV_OK) {
            free(img);
            continue;
        }
        fr[nf].img = half_undistort(img, w, h, K, kr, Khalf);
        free(img);
        if (!fr[nf].img)
            continue;
        fr[nf].cam = cam;
        memcpy(fr[nf].cam.K, Khalf, sizeof(Khalf));
        fr[nf].w = w / 2;
        fr[nf].h = h / 2;
        mv_cam_center(fr[nf].C, &fr[nf].cam);
        printf("%-26s: pose ok, C = (%+.3f %+.3f %+.3f)\n", argv[f],
               fr[nf].C[0], fr[nf].C[1], fr[nf].C[2]);
        nf++;
    }
    printf("\nframes with pose: %d\n", nf);
    if (nf < 2) {
        fprintf(stderr, "need at least 2 posed frames\n");
        return 1;
    }

    mv_cloud_init(&cloud, 1);
    for (i = 1; i < nf && npairs < MAXPAIR; i++) {
        int j, best = -1;
        double bestscore = HUGE_VAL;
        for (j = lastused + 1; j < i; j++) {
            double B, tr = 0.0, ang, score;
            double Rt[9], Rd[9];
            B = mv_baseline(&fr[j].cam, &fr[i].cam);
            mv_mat_transpose(Rt, fr[j].cam.R, 3, 3);
            mv_mat_mul(Rd, fr[i].cam.R, Rt, 3, 3, 3);
            tr = Rd[0] + Rd[4] + Rd[8];
            ang = acos(fmax(-1.0, fmin(1.0, (tr - 1.0) / 2.0)));
            if (B < 0.02 || B > 0.12 || ang > 10.0 * MV_PI / 180.0)
                continue;
            score = fabs(B - 0.06) + 0.3 * ang;
            if (score < bestscore) {
                bestscore = score;
                best = j;
            }
        }
        if (best < 0)
            continue;
        /* stereo on the pair (best, i) */
        {
            frame_t *L = &fr[best], *R = &fr[i];
            int w = L->w, h = L->h, x, y;
            mv_camera r1, r2;
            double H1[9], H2[9], Kinv[9], C1[3], fnew, Bnew;
            unsigned char *rectL, *rectR;
            float *disp;
            int np0 = cloud.n;
            if (mv_rectify_pair(&L->cam, &R->cam, &r1, &r2, H1, H2)
                != MV_OK)
                continue;
            fnew = r1.K[0];
            Bnew = mv_baseline(&r1, &r2);
            mv_mat_inv3(Kinv, r1.K);
            mv_cam_center(C1, &r1);
            rectL = (unsigned char *)malloc((size_t)w * h);
            rectR = (unsigned char *)malloc((size_t)w * h);
            disp = (float *)malloc((size_t)w * h * sizeof(float));
            if (!rectL || !rectR || !disp) {
                free(rectL);
                free(rectR);
                free(disp);
                break;
            }
            mv_warp_homography(rectL, w, h, L->img, w, h, H1, 0);
            mv_warp_homography(rectR, w, h, R->img, w, h, H2, 0);
            mv_stereo_sad(disp, rectL, rectR, w, h, 4, 2, 240);
            texture_gate(disp, rectL, w, h);
            speckle_gate(disp, w, h);
            for (y = 0; y < h; y += 2)
                for (x = 0; x < w; x += 2) {
                    double Z, pc[3], pw[3];
                    unsigned char rgb[3];
                    int j2;
                    float d = disp[y * w + x];
                    if (d <= 0.0f)
                        continue;
                    Z = mv_disp_to_depth(fnew, Bnew, d);
                    if (!(Z >= ZMIN && Z <= ZMAX)) /* also drops NaN */
                        continue;
                    pc[0] = (Kinv[0] * x + Kinv[1] * y + Kinv[2]) * Z;
                    pc[1] = (Kinv[4] * y + Kinv[5]) * Z;
                    pc[2] = Z;
                    for (j2 = 0; j2 < 3; j2++)
                        pw[j2] = r1.R[0 * 3 + j2] * pc[0]
                               + r1.R[1 * 3 + j2] * pc[1]
                               + r1.R[2 * 3 + j2] * pc[2] + C1[j2];
                    rgb[0] = rgb[1] = rgb[2] = rectL[y * w + x];
                    if (mv_cloud_push(&cloud, pw, rgb) == MV_OK) {
                        if (ntag == captag) {
                            captag = captag ? 2 * captag : 65536;
                            tag = (unsigned char *)realloc(tag,
                                                           (size_t)captag);
                        }
                        if (tag)
                            tag[ntag++] = (unsigned char)npairs;
                    }
                }
            printf("pair (%d,%d): B = %.3f m -> %d points\n", best, i,
                   Bnew, cloud.n - np0);
            free(rectL);
            free(rectR);
            free(disp);
            npairs++;
            lastused = i;
        }
    }
    printf("\ntotal points: %d from %d pairs\n", cloud.n, npairs);
    if (!cloud.n)
        return 1;

    /* cross-pair confirmation: a stereo false match is random along its
     * ray, so it will not be reproduced by an independent pair; real
     * structure will.  Keep a point only if its 10 mm voxel neighborhood
     * contains points from at least two distinct pairs. */
    if (npairs >= 2 && tag && ntag == cloud.n) {
        enum { HBITS = 21 };
        const int HSIZE = 1 << HBITS;
        const double VOX = 0.010;
        long long *keys;
        unsigned short *bits;
        int p, kept = 0;
        keys = (long long *)calloc((size_t)HSIZE, sizeof(long long));
        bits = (unsigned short *)calloc((size_t)HSIZE,
                                        sizeof(unsigned short));
        if (keys && bits) {
            for (p = 0; p < cloud.n; p++) {
                long long ix = (long long)floor(cloud.xyz[3 * p] / VOX)
                               + 4096;
                long long iy = (long long)floor(cloud.xyz[3 * p + 1] / VOX)
                               + 4096;
                long long iz = (long long)floor(cloud.xyz[3 * p + 2] / VOX)
                               + 4096;
                long long k = ((ix << 26) | (iy << 13) | iz) + 1;
                unsigned int h = (unsigned int)((unsigned long long)k
                                 * 2654435761ULL >> (32 - HBITS))
                                 & (unsigned int)(HSIZE - 1);
                while (keys[h] && keys[h] != k)
                    h = (h + 1) & (unsigned int)(HSIZE - 1);
                keys[h] = k;
                bits[h] |= (unsigned short)(1u << tag[p]);
            }
            for (p = 0; p < cloud.n; p++) {
                long long ix = (long long)floor(cloud.xyz[3 * p] / VOX)
                               + 4096;
                long long iy = (long long)floor(cloud.xyz[3 * p + 1] / VOX)
                               + 4096;
                long long iz = (long long)floor(cloud.xyz[3 * p + 2] / VOX)
                               + 4096;
                unsigned int m = 0;
                int dx, dy, dz, nb = 0;
                for (dz = -1; dz <= 1; dz++)
                    for (dy = -1; dy <= 1; dy++)
                        for (dx = -1; dx <= 1; dx++) {
                            long long k = (((ix + dx) << 26)
                                           | ((iy + dy) << 13)
                                           | (iz + dz)) + 1;
                            unsigned int h = (unsigned int)
                                ((unsigned long long)k * 2654435761ULL
                                 >> (32 - HBITS))
                                & (unsigned int)(HSIZE - 1);
                            while (keys[h] && keys[h] != k)
                                h = (h + 1) & (unsigned int)(HSIZE - 1);
                            if (keys[h] == k)
                                m |= bits[h];
                        }
                while (m) {
                    nb += (int)(m & 1u);
                    m >>= 1;
                }
                if (nb >= 2) {
                    cloud.xyz[3 * kept] = cloud.xyz[3 * p];
                    cloud.xyz[3 * kept + 1] = cloud.xyz[3 * p + 1];
                    cloud.xyz[3 * kept + 2] = cloud.xyz[3 * p + 2];
                    cloud.rgb[3 * kept] = cloud.rgb[3 * p];
                    cloud.rgb[3 * kept + 1] = cloud.rgb[3 * p + 1];
                    cloud.rgb[3 * kept + 2] = cloud.rgb[3 * p + 2];
                    kept++;
                }
            }
            printf("cross-pair confirmation: %d of %d points kept\n",
                   kept, cloud.n);
            cloud.n = kept;
        }
        free(keys);
        free(bits);
    }
    free(tag);
    if (!cloud.n)
        return 1;
    {
        char path[512];
        double lo[3], hi[3], ctr[3], eye[3], span;
        double *ax;
        int p;
        int j;
        snprintf(path, sizeof(path), "%s.ply", prefix);
        mv_cloud_write_ply(path, &cloud);
        printf("wrote %s\n", path);
        /* frame on the 3rd..97th percentile so stragglers cannot zoom
         * the camera out to nothing */
        ax = (double *)malloc((size_t)cloud.n * sizeof(double));
        for (j = 0; j < 3; j++) {
            if (!ax)
                break;
            for (p = 0; p < cloud.n; p++)
                ax[p] = cloud.xyz[3 * p + j];
            qsort(ax, (size_t)cloud.n, sizeof(double), cmp_dbl);
            lo[j] = ax[(int)(0.03 * (cloud.n - 1))];
            hi[j] = ax[(int)(0.97 * (cloud.n - 1))];
        }
        free(ax);
        for (j = 0; j < 3; j++)
            ctr[j] = 0.5 * (lo[j] + hi[j]);
        span = fmax(hi[0] - lo[0], fmax(hi[1] - lo[1], hi[2] - lo[2]));
        /* three orbit views around the cloud center */
        for (i = 0; i < 3; i++) {
            double az = (i - 1) * 0.7; /* -40, 0, +40 degrees-ish */
            eye[0] = ctr[0] + 1.6 * span * sin(az);
            eye[1] = ctr[1] - 0.35 * span;
            eye[2] = ctr[2] - 1.6 * span * cos(az);
            snprintf(path, sizeof(path), "%s_view%d.pgm", prefix, i);
            splat_view(path, &cloud, eye, ctr, 900.0, 960, 720);
            printf("wrote %s\n", path);
        }
    }
    mv_cloud_free(&cloud);
    for (f = 0; f < nf; f++)
        free(fr[f].img);
    return 0;
}
