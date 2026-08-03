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
