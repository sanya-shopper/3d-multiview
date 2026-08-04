#include <stdlib.h>
#include "mv/core.h"
#include "mv/graycode.h"

int mv_graycode_bits(int extent)
{
    int b = 0;
    while ((1 << b) < extent && b < 31)
        b++;
    return b;
}

unsigned mv_gray_encode(unsigned v)
{
    return v ^ (v >> 1);
}

unsigned mv_gray_decode(unsigned g)
{
    unsigned b = g;
    b ^= b >> 1;
    b ^= b >> 2;
    b ^= b >> 4;
    b ^= b >> 8;
    b ^= b >> 16;
    return b;
}

void mv_graycode_frame(unsigned char *img, int w, int h, int axis, int bit,
                       int nbits, int inverse)
{
    int x, y;
    int shift = nbits - 1 - bit;
    for (y = 0; y < h; y++) {
        for (x = 0; x < w; x++) {
            unsigned c = (unsigned)(axis ? y : x);
            unsigned on = (mv_gray_encode(c) >> shift) & 1u;
            if (inverse)
                on = !on;
            img[y * w + x] = on ? 255 : 0;
        }
    }
}

int mv_graycode_decode(int *coord, unsigned char *valid,
                       const unsigned char *const *frames,
                       const unsigned char *const *inv_frames,
                       int nbits, int npix, int min_contrast)
{
    int p, b, nvalid = 0;
    for (p = 0; p < npix; p++) {
        unsigned g = 0;
        int ok = 1;
        for (b = 0; b < nbits; b++) {
            int a = frames[b][p], i = inv_frames[b][p];
            int diff = a - i;
            if (diff < 0)
                diff = -diff;
            if (diff < min_contrast) {
                ok = 0;
                break;
            }
            g = (g << 1) | (unsigned)(frames[b][p] > inv_frames[b][p]);
        }
        if (ok) {
            coord[p] = (int)mv_gray_decode(g);
            valid[p] = 1;
            nvalid++;
        } else {
            coord[p] = -1;
            valid[p] = 0;
        }
    }
    return nvalid;
}

int mv_graycode_match(double *uv1, double *uv2, int maxm,
                      const int *colc1, const int *rowc1,
                      const unsigned char *valid1, int w1, int h1,
                      const int *colc2, const int *rowc2,
                      const unsigned char *valid2, int w2, int h2,
                      int cw, int ch)
{
    /* codes are binned into BLOCK x BLOCK projector-pixel patches; a
     * camera pixel spans a few projector pixels anyway, and the two
     * cameras' centroids of one patch correspond to sub-pixel level */
    enum { BLOCK = 2 };
    int bw = (cw + BLOCK - 1) / BLOCK, bh = (ch + BLOCK - 1) / BLOCK;
    double *acc; /* per bin: su1 sv1 n1 su2 sv2 n2 */
    int i, x, y, nm = 0;

    acc = (double *)calloc((size_t)bw * bh * 6, sizeof(double));
    if (!acc)
        return 0;
    for (y = 0; y < h1; y++)
        for (x = 0; x < w1; x++) {
            int p = y * w1 + x, b;
            if (!valid1[p])
                continue;
            if (colc1[p] < 0 || colc1[p] >= cw || rowc1[p] < 0
                || rowc1[p] >= ch)
                continue;
            b = (rowc1[p] / BLOCK) * bw + colc1[p] / BLOCK;
            acc[6 * b + 0] += x;
            acc[6 * b + 1] += y;
            acc[6 * b + 2] += 1.0;
        }
    for (y = 0; y < h2; y++)
        for (x = 0; x < w2; x++) {
            int p = y * w2 + x, b;
            if (!valid2[p])
                continue;
            if (colc2[p] < 0 || colc2[p] >= cw || rowc2[p] < 0
                || rowc2[p] >= ch)
                continue;
            b = (rowc2[p] / BLOCK) * bw + colc2[p] / BLOCK;
            acc[6 * b + 3] += x;
            acc[6 * b + 4] += y;
            acc[6 * b + 5] += 1.0;
        }
    for (i = 0; i < bw * bh && nm < maxm; i++) {
        if (acc[6 * i + 2] < 1.0 || acc[6 * i + 5] < 1.0)
            continue;
        uv1[2 * nm] = acc[6 * i + 0] / acc[6 * i + 2];
        uv1[2 * nm + 1] = acc[6 * i + 1] / acc[6 * i + 2];
        uv2[2 * nm] = acc[6 * i + 3] / acc[6 * i + 5];
        uv2[2 * nm + 1] = acc[6 * i + 4] / acc[6 * i + 5];
        nm++;
    }
    free(acc);
    return nm;
}
