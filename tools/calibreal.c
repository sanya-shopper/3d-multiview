/* Real-data calibration: PGM frames of the spec-v1 pattern (shown by
 * tools/pattern.html on any display) -> blind reads -> Zhang.
 *
 * Usage: calibreal <pitch_mm> <frame.pgm> [more.pgm ...]
 *   pitch_mm: the display's physical pixel pitch in mm (measure one
 *   80-px grid square on screen and divide by 80; intrinsics K do not
 *   depend on it, only the reported metric distances do).
 *
 * Prints per-frame read results (corners, orientation, decoded display
 * counter) and, from all frames with a solid read, the camera
 * intrinsics, distortion, reprojection RMS, and the camera-to-screen
 * distance per view (verifiable with a tape measure).
 *
 * Build: make calibreal */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "mv/mv.h"

#define MAXV 220

int main(int argc, char **argv)
{
    static double objs[MAXV][2 * MV_READ_MAXC];
    static double imgs[MAXV][2 * MV_READ_MAXC];
    mv_calib_view views[MAXV];
    mv_camera est[MAXV];
    double K[9], k_radial[2], pitch;
    int nviews = 0, i, f;

    if (argc < 3) {
        fprintf(stderr, "usage: %s <pitch_mm> <frames.pgm...>\n", argv[0]);
        return 1;
    }
    pitch = atof(argv[1]) * 1e-3;

    for (f = 2; f < argc; f++) {
        unsigned char *img;
        int w, h;
        mv_read_result rr;
        if (mv_pgm_read(argv[f], &img, &w, &h) != MV_OK) {
            printf("%-28s : unreadable\n", argv[f]);
            continue;
        }
        if (mv_read_pattern(&rr, img, w, h) != MV_OK) {
            printf("%-28s : no pattern found\n", argv[f]);
            free(img);
            continue;
        }
        printf("%-28s : %3d corners, rot %d, counter %s%u (conf %.2f)\n",
               argv[f], rr.n, rr.rot, rr.counter_valid ? "" : "INVALID ",
               rr.counter, rr.counter_conf);
        if (getenv("MV_DUMP_CORNERS")) {
            for (i = 0; i < rr.n; i++)
                fprintf(stderr, "CRN %d %d %.3f %.3f\n", f - 2,
                        rr.id[i], rr.uv[2 * i], rr.uv[2 * i + 1]);
        }
        if (rr.n >= 50 && nviews < MAXV) {
            for (i = 0; i < rr.n; i++) {
                double xy[2];
                mv_pattern_corner_px(rr.id[i] % MV_PAT_CORNER_COLS,
                                     rr.id[i] / MV_PAT_CORNER_COLS, xy);
                objs[nviews][2 * i] = xy[0] * pitch;
                objs[nviews][2 * i + 1] = xy[1] * pitch;
                imgs[nviews][2 * i] = rr.uv[2 * i];
                imgs[nviews][2 * i + 1] = rr.uv[2 * i + 1];
            }
            views[nviews].obj = objs[nviews];
            views[nviews].img = imgs[nviews];
            views[nviews].n = rr.n;
            nviews++;
        }
        free(img);
    }

    printf("\nviews accepted: %d\n", nviews);
    if (nviews < 3) {
        printf("need >= 3 good views for calibration\n");
        return 1;
    }
    if (mv_calib_planar(K, est, k_radial, views, nviews, 1) != MV_OK) {
        printf("calibration FAILED\n");
        return 1;
    }
    printf("\nintrinsics (pixels)\n");
    printf("  fx = %.1f   fy = %.1f\n", K[0], K[4]);
    printf("  cx = %.1f   cy = %.1f\n", K[2], K[5]);
    printf("  k1 = %+.4f  k2 = %+.4f\n", k_radial[0], k_radial[1]);
    printf("  reprojection RMS = %.3f px over %d views\n",
           mv_calib_reproj_rms(est, views, nviews), nviews);
    printf("\nper-view camera-to-screen-center distance "
           "(at pitch %.4f mm):\n  ", pitch * 1000.0);
    for (i = 0; i < nviews; i++) {
        double c[3] = { 0, 0, 0 }, Xc2;
        c[0] = MV_PAT_W / 2.0 * pitch;
        c[1] = MV_PAT_H / 2.0 * pitch;
        Xc2 = est[i].R[6] * c[0] + est[i].R[7] * c[1] + est[i].t[2];
        printf("%.2f ", Xc2);
    }
    printf("m\n");
    return 0;
}
