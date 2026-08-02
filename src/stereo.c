#include <limits.h>

#include "mv/core.h"
#include "mv/stereo.h"

static long sad_cost(const unsigned char *left, const unsigned char *right,
                     int w, int x, int y, int d, int halfwin)
{
    long cost = 0;
    int i, j;
    for (j = -halfwin; j <= halfwin; j++) {
        const unsigned char *lrow = left + (y + j) * w;
        const unsigned char *rrow = right + (y + j) * w;
        for (i = -halfwin; i <= halfwin; i++) {
            int diff = (int)lrow[x + i] - (int)rrow[x - d + i];
            cost += (diff < 0) ? -diff : diff;
        }
    }
    return cost;
}

int mv_stereo_sad(float *disp, const unsigned char *left,
                  const unsigned char *right, int w, int h,
                  int halfwin, int dmin, int dmax)
{
    int x, y, d, i;

    if (w <= 0 || h <= 0 || halfwin < 0 || dmin < 0 || dmax < dmin)
        return MV_ERR;

    for (i = 0; i < w * h; i++)
        disp[i] = -1.0f;

    for (y = halfwin; y < h - halfwin; y++) {
        for (x = halfwin; x < w - halfwin; x++) {
            long best_cost = LONG_MAX;
            int best_d = -1;
            for (d = dmin; d <= dmax; d++) {
                long cost;
                if (x - d - halfwin < 0)
                    break;
                cost = sad_cost(left, right, w, x, y, d, halfwin);
                if (cost < best_cost) {
                    best_cost = cost;
                    best_d = d;
                }
            }
            if (best_d < 0)
                continue;
            disp[y * w + x] = (float)best_d;
            /* sub-pixel refinement: parabola through the costs at
             * best_d - 1, best_d, best_d + 1 */
            if (best_d > dmin && best_d < dmax
                && x - (best_d + 1) - halfwin >= 0) {
                long cp = sad_cost(left, right, w, x, y, best_d - 1, halfwin);
                long cn = sad_cost(left, right, w, x, y, best_d + 1, halfwin);
                long denom = cp - 2 * best_cost + cn;
                if (denom > 0) {
                    double frac = 0.5 * (double)(cp - cn) / (double)denom;
                    if (frac > -1.0 && frac < 1.0)
                        disp[y * w + x] = (float)(best_d + frac);
                }
            }
        }
    }
    return MV_OK;
}
