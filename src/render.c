#include <math.h>

#include "mv/core.h"
#include "mv/mat.h"
#include "mv/render.h"

#define MV_PI_L 3.14159265358979323846

static double urand(unsigned long long *s)
{
    *s = *s * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)((*s >> 11) & ((1ULL << 53) - 1)) / 9007199254740992.0;
}

static double grand_(unsigned long long *s)
{
    double u1 = urand(s), u2 = urand(s);
    if (u1 < 1e-300)
        u1 = 1e-300;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * MV_PI_L * u2);
}

int mv_render_plane(unsigned char *out, int w, int h, const mv_camera *cam,
                    const unsigned char *src, int sw, int sh, double pitch,
                    unsigned char background, double sigma,
                    unsigned long long *seed)
{
    /* homography image <- plane(metres): H = K [r1 r2 t]; invert it and
     * map each camera pixel back into the source */
    double H[9], Hinv[9];
    int x, y, i;

    for (i = 0; i < 3; i++) {
        double col1 = cam->R[i * 3 + 0], col2 = cam->R[i * 3 + 1];
        H[i * 3 + 0] = col1;
        H[i * 3 + 1] = col2;
        H[i * 3 + 2] = cam->t[i];
    }
    {
        double KH[9];
        mv_mat_mul(KH, cam->K, H, 3, 3, 3);
        if (mv_mat_inv3(Hinv, KH) != MV_OK)
            return MV_ERR;
    }

    for (y = 0; y < h; y++) {
        for (x = 0; x < w; x++) {
            double X = Hinv[0] * x + Hinv[1] * y + Hinv[2];
            double Y = Hinv[3] * x + Hinv[4] * y + Hinv[5];
            double W = Hinv[6] * x + Hinv[7] * y + Hinv[8];
            double v = background;
            if (fabs(W) > 1e-12) {
                double sx = (X / W) / pitch, sy = (Y / W) / pitch;
                if (sx >= 0.0 && sx < sw - 1.0 && sy >= 0.0
                    && sy < sh - 1.0) {
                    int ix = (int)sx, iy = (int)sy;
                    double fx = sx - ix, fy = sy - iy;
                    const unsigned char *p = src + iy * sw + ix;
                    v = (1 - fx) * (1 - fy) * p[0] + fx * (1 - fy) * p[1]
                      + (1 - fx) * fy * p[sw] + fx * fy * p[sw + 1];
                }
            }
            v += sigma * grand_(seed);
            out[y * w + x] = (unsigned char)(v < 0.0 ? 0.0
                                             : v > 255.0 ? 255.0 : v);
        }
    }
    return MV_OK;
}

int mv_render_scene(unsigned char *img, float *zc, int w, int h,
                    const mv_camera *cam, const mv_plane *pls, int n,
                    unsigned char background, double sigma,
                    unsigned long long *seed)
{
    double C[3];
    int x, y, i, p;

    mv_cam_center(C, cam);
    for (y = 0; y < h; y++)
        for (x = 0; x < w; x++) {
            double uv[2] = { (double)x, (double)y };
            double o[3], d[3];
            double best = HUGE_VAL, v = background;
            if (mv_cam_ray(o, d, cam, uv) != MV_OK)
                continue;
            for (p = 0; p < n; p++) {
                const mv_plane *pl = &pls[p];
                double op[3], dp[3], s, px, py, sx, sy;
                int ix, iy;
                for (i = 0; i < 3; i++) {
                    op[i] = pl->R[0 + i] * (C[0] - pl->t[0])
                          + pl->R[3 + i] * (C[1] - pl->t[1])
                          + pl->R[6 + i] * (C[2] - pl->t[2]);
                    dp[i] = pl->R[0 + i] * d[0] + pl->R[3 + i] * d[1]
                          + pl->R[6 + i] * d[2];
                }
                if (fabs(dp[2]) < 1e-9)
                    continue;
                s = -op[2] / dp[2];
                if (s <= 1e-6 || s >= best)
                    continue;
                px = op[0] + s * dp[0];
                py = op[1] + s * dp[1];
                sx = px / pl->pitch;
                sy = py / pl->pitch;
                if (sx < 0.0 || sy < 0.0 || sx >= pl->tw - 1.0
                    || sy >= pl->th - 1.0)
                    continue;
                ix = (int)sx;
                iy = (int)sy;
                {
                    double fx = sx - ix, fy = sy - iy;
                    const unsigned char *q = pl->tex + iy * pl->tw + ix;
                    v = (1 - fx) * (1 - fy) * q[0] + fx * (1 - fy) * q[1]
                      + (1 - fx) * fy * q[pl->tw] + fx * fy * q[pl->tw + 1];
                    best = s;
                }
            }
            if (zc) {
                if (best < HUGE_VAL) {
                    double Z = cam->R[6] * d[0] + cam->R[7] * d[1]
                             + cam->R[8] * d[2];
                    zc[y * w + x] = (float)(best * Z);
                } else
                    zc[y * w + x] = HUGE_VALF;
            }
            v += sigma * grand_(seed);
            img[y * w + x] = (unsigned char)(v < 0 ? 0 : v > 255 ? 255 : v);
        }
    return MV_OK;
}

int mv_warp_homography(unsigned char *dst, int dw, int dh,
                       const unsigned char *src, int sw, int sh,
                       const double H[9], unsigned char fill)
{
    double Hinv[9];
    int x, y;
    if (mv_mat_inv3(Hinv, H) != MV_OK)
        return MV_ERR;
    for (y = 0; y < dh; y++)
        for (x = 0; x < dw; x++) {
            double X = Hinv[0] * x + Hinv[1] * y + Hinv[2];
            double Y = Hinv[3] * x + Hinv[4] * y + Hinv[5];
            double W = Hinv[6] * x + Hinv[7] * y + Hinv[8];
            double v = fill;
            if (fabs(W) > 1e-12) {
                double sx = X / W, sy = Y / W;
                if (sx >= 0.0 && sy >= 0.0 && sx < sw - 1.0
                    && sy < sh - 1.0) {
                    int ix = (int)sx, iy = (int)sy;
                    double fx = sx - ix, fy = sy - iy;
                    const unsigned char *q = src + iy * sw + ix;
                    v = (1 - fx) * (1 - fy) * q[0] + fx * (1 - fy) * q[1]
                      + (1 - fx) * fy * q[sw] + fx * fy * q[sw + 1];
                }
            }
            dst[y * dw + x] = (unsigned char)(v + 0.5);
        }
    return MV_OK;
}
