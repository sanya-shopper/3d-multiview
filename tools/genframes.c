/* Synthetic camera-view generator: renders the fine-tier pattern from a
 * sweep of poses (the test_calib tilt/distance idiom) into PGM frames.
 * Purpose: CI and offline stress runs need camera footage without any
 * camera -- replaycam streams these into livehub and the hub must
 * calibrate. Deterministic.
 *
 * Usage: genframes <outdir> <count>
 * Build: make genframes */

#define _POSIX_C_SOURCE 200809L
#define _DARWIN_C_SOURCE 1

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "mv/mv.h"

int main(int argc, char **argv)
{
    static unsigned char pat[MV_PAT_W * MV_PAT_H];
    static unsigned char img[640 * 480];
    const double pitch = 0.0002745;
    unsigned long long seed = 20260807ULL;
    double center[3] = { MV_PAT_W / 2.0 * 0.0002745,
                         MV_PAT_H / 2.0 * 0.0002745, 0.0 };
    int count, f, j;

    if (argc != 3) {
        fprintf(stderr, "usage: %s <outdir> <count>\n", argv[0]);
        return 1;
    }
    count = atoi(argv[2]);
    mkdir(argv[1], 0755);

    for (f = 0; f < count; f++) {
        mv_camera cam;
        /* pose sweep: tilts to ~20 degrees on both axes, distance
         * stepping 0.75..0.95 m -- comfortably diverse for Zhang */
        double ax = 0.35 * sin(0.7 * f), ay = 0.30 * cos(1.1 * f);
        double cax = cos(ax), sax = sin(ax), cay = cos(ay),
               say = sin(ay);
        double Rx[9], Ry[9];
        char path[512];
        Rx[0] = 1; Rx[1] = 0;   Rx[2] = 0;
        Rx[3] = 0; Rx[4] = cax; Rx[5] = -sax;
        Rx[6] = 0; Rx[7] = sax; Rx[8] = cax;
        Ry[0] = cay; Ry[1] = 0; Ry[2] = -say;
        Ry[3] = 0;   Ry[4] = 1; Ry[5] = 0;
        Ry[6] = say; Ry[7] = 0; Ry[8] = cay;
        mv_cam_set_K(&cam, 800.0, 800.0, 320.0, 240.0);
        memset(cam.k, 0, sizeof(cam.k));
        mv_mat_mul(cam.R, Rx, Ry, 3, 3, 3);
        for (j = 0; j < 3; j++)
            cam.t[j] = -(cam.R[j * 3 + 0] * center[0]
                         + cam.R[j * 3 + 1] * center[1]
                         + cam.R[j * 3 + 2] * center[2]);
        cam.t[2] += 0.75 + 0.02 * (f % 11);
        mv_pattern_render(pat, (unsigned)(100000 + 30 * f));
        if (mv_render_plane(img, 640, 480, &cam, pat, MV_PAT_W,
                            MV_PAT_H, pitch, 40, 1.0, &seed) != MV_OK) {
            fprintf(stderr, "render failed at frame %d\n", f);
            return 1;
        }
        snprintf(path, sizeof(path), "%s/f%04d.pgm", argv[1], f);
        if (mv_pgm_write(path, img, 640, 480) != MV_OK) {
            fprintf(stderr, "cannot write %s\n", path);
            return 1;
        }
    }
    printf("wrote %d frames to %s\n", count, argv[1]);
    return 0;
}
