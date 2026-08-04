/* Annotate a camera frame with everything the system detected:
 * green crosses  = identified pattern corners (sub-pixel)
 * cyan lines     = the decoded corner lattice
 * yellow quad    = the pattern's extent, projected through the fitted
 *                  homography
 * red quad       = the physical display's lit-panel outline, DETECTED
 *                  by mv_display_outline (flood of the panel glow +
 *                  boundary line fits) -- found in the image, not
 *                  assumed from the spec.
 *
 * Usage: annotate <frame.pgm> <out.ppm>
 * Build: make annotate */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "mv/mv.h"

static unsigned char *rgb;
static int W, H;

static void put(int x, int y, int r, int g, int b)
{
    int i;
    if (x < 0 || y < 0 || x >= W || y >= H)
        return;
    i = 3 * (y * W + x);
    rgb[i] = (unsigned char)r;
    rgb[i + 1] = (unsigned char)g;
    rgb[i + 2] = (unsigned char)b;
}

static void thickline(double x0, double y0, double x1, double y1,
                      int r, int g, int b, int rad)
{
    double dx = x1 - x0, dy = y1 - y0;
    double len = sqrt(dx * dx + dy * dy);
    int n = (int)(len + 1.0), i, a, c;
    for (i = 0; i <= n; i++) {
        int cx = (int)(x0 + dx * i / n + 0.5);
        int cy = (int)(y0 + dy * i / n + 0.5);
        for (a = -rad; a <= rad; a++)
            for (c = -rad; c <= rad; c++)
                if (a * a + c * c <= rad * rad)
                    put(cx + c, cy + a, r, g, b);
    }
}

static void projH(const double Hm[9], double X, double Y, double *u,
                  double *v)
{
    double Wd = Hm[6] * X + Hm[7] * Y + Hm[8];
    if (fabs(Wd) < 1e-12)
        Wd = 1e-12;
    *u = (Hm[0] * X + Hm[1] * Y + Hm[2]) / Wd;
    *v = (Hm[3] * X + Hm[4] * Y + Hm[5]) / Wd;
}

int main(int argc, char **argv)
{
    unsigned char *img;
    mv_read_result rr;
    int tier = 1, i, x, y;

    if (argc != 3) {
        fprintf(stderr, "usage: %s <frame.pgm> <out.ppm>\n", argv[0]);
        return 1;
    }
    if (mv_pgm_read(argv[1], &img, &W, &H) != MV_OK) {
        fprintf(stderr, "cannot read %s\n", argv[1]);
        return 1;
    }
    if (mv_read_pattern(&rr, img, W, H) != MV_OK) {
        if (mv_read_coarse(&rr, img, W, H) != MV_OK) {
            fprintf(stderr, "no pattern found in %s\n", argv[1]);
            return 1;
        }
        tier = 2;
    }
    printf("%s: %d corners (%s tier), rot %d, counter %s%u\n", argv[1],
           rr.n, tier == 1 ? "fine" : "coarse", rr.rot,
           rr.counter_valid ? "" : "INVALID ", rr.counter);

    rgb = (unsigned char *)malloc((size_t)W * H * 3);
    if (!rgb)
        return 1;
    for (y = 0; y < H; y++)
        for (x = 0; x < W; x++) {
            unsigned char g = img[y * W + x];
            rgb[3 * (y * W + x)] = g;
            rgb[3 * (y * W + x) + 1] = g;
            rgb[3 * (y * W + x) + 2] = g;
        }

    /* red: DETECTED lit-panel outline (drawn first, under the rest) */
    {
        double quad[8];
        if (mv_display_outline(quad, &rr, img, W, H) == MV_OK) {
            for (i = 0; i < 4; i++) {
                int j = (i + 1) % 4;
                thickline(quad[2 * i], quad[2 * i + 1], quad[2 * j],
                          quad[2 * j + 1], 255, 40, 40, 2);
            }
            printf("display outline: (%.0f,%.0f) (%.0f,%.0f) "
                   "(%.0f,%.0f) (%.0f,%.0f)\n", quad[0], quad[1],
                   quad[2], quad[3], quad[4], quad[5], quad[6], quad[7]);
        } else {
            printf("display outline: not found\n");
        }
    }

    /* yellow: pattern extent through the homography */
    {
        static const double ext[4][2] = {
            { 0, 0 }, { MV_PAT_W, 0 }, { MV_PAT_W, MV_PAT_H },
            { 0, MV_PAT_H }
        };
        for (i = 0; i < 4; i++) {
            int j = (i + 1) % 4;
            double u0, v0, u1, v1;
            projH(rr.H, ext[i][0], ext[i][1], &u0, &v0);
            projH(rr.H, ext[j][0], ext[j][1], &u1, &v1);
            thickline(u0, v0, u1, v1, 255, 220, 40, 1);
        }
    }

    /* cyan: the identified lattice */
    {
        int cols = tier == 1 ? MV_PAT_CORNER_COLS : MV_PAT2_CORNER_COLS;
        int rows = tier == 1 ? MV_PAT_CORNER_ROWS : MV_PAT2_CORNER_ROWS;
        static double px[256][2];
        static int have[256];
        int a, b;
        for (i = 0; i < 256; i++)
            have[i] = 0;
        for (i = 0; i < rr.n; i++) {
            px[rr.id[i]][0] = rr.uv[2 * i];
            px[rr.id[i]][1] = rr.uv[2 * i + 1];
            have[rr.id[i]] = 1;
        }
        for (b = 0; b < rows; b++)
            for (a = 0; a < cols; a++) {
                int id = b * cols + a;
                if (!have[id])
                    continue;
                if (a + 1 < cols && have[id + 1])
                    thickline(px[id][0], px[id][1], px[id + 1][0],
                              px[id + 1][1], 40, 200, 255, 0);
                if (b + 1 < rows && have[id + cols])
                    thickline(px[id][0], px[id][1], px[id + cols][0],
                              px[id + cols][1], 40, 200, 255, 0);
            }
    }

    /* green: sub-pixel corners on top */
    for (i = 0; i < rr.n; i++) {
        double u = rr.uv[2 * i], v = rr.uv[2 * i + 1];
        int t;
        for (t = -3; t <= 3; t++) {
            put((int)(u + t + 0.5), (int)(v + 0.5), 40, 255, 60);
            put((int)(u + 0.5), (int)(v + t + 0.5), 40, 255, 60);
        }
    }

    {
        FILE *fp = fopen(argv[2], "wb");
        if (!fp)
            return 1;
        fprintf(fp, "P6\n%d %d\n255\n", W, H);
        fwrite(rgb, 3, (size_t)W * H, fp);
        fclose(fp);
    }
    printf("wrote %s\n", argv[2]);
    free(rgb);
    free(img);
    return 0;
}
