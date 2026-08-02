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
