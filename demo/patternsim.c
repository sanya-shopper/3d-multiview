/* End-to-end calibration-subsystem simulation: numbers quoted in
 * doc/multiview.tex ("Basic results").
 *
 * The spec-v1 pattern is rendered (mv_pattern_render), shown at metric
 * scale on a virtual display, imaged by the standard camera under six
 * poses (mv_render_plane, bilinear + sensor noise), and decoded blind by
 * the reader (mv_read_pattern). Checks: corner identification vs ground
 * truth, sub-pixel accuracy, counter decode, orientation recovery under a
 * 180-degree camera roll, and finally full image-tier Zhang calibration
 * through the rendered pipeline. Deterministic. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mv/mv.h"

#define MV_PI 3.14159265358979323846
#define CAMW 640
#define CAMH 480
#define PITCH 0.0002745      /* m per display pixel */
#define COUNTER 48213u
#define SIGMA 2.0            /* sensor noise, gray levels */

static void make_pose(mv_camera *cam, double pitch_deg, double yaw_deg,
                      double dx, double dy, double dist, int roll180)
{
    double cx = cos(pitch_deg * MV_PI / 180.0);
    double sx = sin(pitch_deg * MV_PI / 180.0);
    double cy = cos(yaw_deg * MV_PI / 180.0);
    double sy = sin(yaw_deg * MV_PI / 180.0);
    double Rx[9], Ry[9], R[9];
    double center[3] = { MV_PAT_W / 2.0 * PITCH, MV_PAT_H / 2.0 * PITCH,
                         0.0 };
    int i;
    Rx[0] = 1; Rx[1] = 0;  Rx[2] = 0;
    Rx[3] = 0; Rx[4] = cx; Rx[5] = -sx;
    Rx[6] = 0; Rx[7] = sx; Rx[8] = cx;
    Ry[0] = cy; Ry[1] = 0; Ry[2] = -sy;
    Ry[3] = 0;  Ry[4] = 1; Ry[5] = 0;
    Ry[6] = sy; Ry[7] = 0; Ry[8] = cy;
    mv_cam_set_K(cam, 800.0, 800.0, 320.0, 240.0);
    mv_mat_mul(R, Rx, Ry, 3, 3, 3);
    if (roll180)
        for (i = 0; i < 6; i++)
            R[i] = -R[i]; /* Rz(pi) * R negates first two rows */
    memcpy(cam->R, R, sizeof(R));
    for (i = 0; i < 3; i++)
        cam->t[i] = -(R[i * 3 + 0] * center[0] + R[i * 3 + 1] * center[1]);
    cam->t[0] += dx;
    cam->t[1] += dy;
    cam->t[2] += dist;
    memset(cam->k, 0, sizeof(cam->k));
}

/* ground-truth image position of pattern corner id under cam */
static int gt_corner(double uv[2], const mv_camera *cam, int id)
{
    double xy[2], X[3];
    mv_pattern_corner_px(id % MV_PAT_CORNER_COLS,
                         id / MV_PAT_CORNER_COLS, xy);
    X[0] = xy[0] * PITCH;
    X[1] = xy[1] * PITCH;
    X[2] = 0.0;
    return mv_cam_project(uv, cam, X);
}

int main(void)
{
    static unsigned char pat[MV_PAT_W * MV_PAT_H];
    static unsigned char img[CAMW * CAMH];
    /* pitch_deg, yaw_deg, dx, dy, dist, roll180 */
    const double poses[7][6] = {
        { 0, 0, 0.00, 0.00, 0.85, 0 },
        { 25, 0, 0.02, -0.01, 0.80, 0 },
        { -25, 10, -0.02, 0.01, 0.80, 0 },
        { 15, -30, 0.01, 0.02, 0.78, 0 },
        { -15, -20, -0.01, -0.02, 0.80, 0 },
        { 30, 20, 0.02, 0.01, 0.85, 0 },
        { 0, 0, 0.00, 0.00, 0.85, 1 } /* rolled 180 deg */
    };
    static double objs[7][2 * MV_READ_MAXC], imgs[7][2 * MV_READ_MAXC];
    mv_calib_view views[6];
    unsigned long long seed = 99;
    int p, i, nviews = 0;
    int total_corners = 0, total_id_err = 0, ctr_ok = 0;
    double se = 0.0;
    long nse = 0;

    printf("multiview pattern-pipeline experiment (seed 99)\n");
    printf("-----------------------------------------------\n");
    printf("pattern             : spec v1, counter %u, pitch %.4f mm\n",
           COUNTER, PITCH * 1000.0);
    printf("camera              : %dx%d, f=800, noise %.1f gray\n\n",
           CAMW, CAMH, SIGMA);

    if (mv_pattern_selftest() != MV_OK) {
        fprintf(stderr, "M-array selftest FAILED\n");
        return 1;
    }
    mv_pattern_render(pat, COUNTER);

    for (p = 0; p < 7; p++) {
        mv_camera cam;
        mv_read_result rr;
        int id_bad = 0, matched = 0;
        double rms;

        make_pose(&cam, poses[p][0], poses[p][1], poses[p][2], poses[p][3],
                  poses[p][4], (int)poses[p][5]);
        if (mv_render_plane(img, CAMW, CAMH, &cam, pat, MV_PAT_W, MV_PAT_H,
                            PITCH, 128, SIGMA, &seed) != MV_OK) {
            fprintf(stderr, "render failed\n");
            return 1;
        }
        if (p == 0)
            mv_pgm_write("out_patternsim.pgm", img, CAMW, CAMH);

        if (mv_read_pattern(&rr, img, CAMW, CAMH) != MV_OK) {
            printf("pose %d: READ FAILED\n", p);
            continue;
        }
        se = 0.0;
        nse = 0;
        for (i = 0; i < rr.n; i++) {
            double guv[2], du, dv, e2;
            if (gt_corner(guv, &cam, rr.id[i]) != MV_OK)
                continue;
            du = rr.uv[2 * i] - guv[0];
            dv = rr.uv[2 * i + 1] - guv[1];
            e2 = du * du + dv * dv;
            if (e2 > 9.0) {
                id_bad++;
                continue;
            }
            se += e2;
            nse++;
            matched++;
        }
        rms = nse ? sqrt(se / nse) : -1.0;
        printf("pose %d (%+3.0f,%+3.0f deg%s): corners %3d  "
               "id errors %d  subpix RMS %.3f px  rot %d  counter %s\n",
               p, poses[p][0], poses[p][1],
               (int)poses[p][5] ? ", rolled" : "", rr.n, id_bad, rms,
               rr.rot,
               (rr.counter_valid && rr.counter == COUNTER) ? "ok"
                                                           : "BAD");
        total_corners += matched;
        total_id_err += id_bad;
        if (rr.counter_valid && rr.counter == COUNTER)
            ctr_ok++;

        if (p < 6) {
            for (i = 0; i < rr.n; i++) {
                double xy[2];
                mv_pattern_corner_px(rr.id[i] % MV_PAT_CORNER_COLS,
                                     rr.id[i] / MV_PAT_CORNER_COLS, xy);
                objs[p][2 * i] = xy[0] * PITCH;
                objs[p][2 * i + 1] = xy[1] * PITCH;
                imgs[p][2 * i] = rr.uv[2 * i];
                imgs[p][2 * i + 1] = rr.uv[2 * i + 1];
            }
            views[nviews].obj = objs[p];
            views[nviews].img = imgs[p];
            views[nviews].n = rr.n;
            nviews++;
        }
    }

    printf("\ntotals: %d corners matched, %d misidentified, "
           "counter %d/7 frames\n", total_corners, total_id_err, ctr_ok);

    /* --- end-to-end Zhang through the image pipeline --- */
    if (nviews >= 3) {
        mv_camera est[6];
        double K[9], k_radial[2];
        if (mv_calib_planar(K, est, k_radial, views, nviews, 1) != MV_OK) {
            fprintf(stderr, "calibration failed\n");
            return 1;
        }
        printf("\nimage-tier Zhang calibration (%d rendered views)\n",
               nviews);
        printf("  fx : %8.2f  (800.0)  error %+.2f px\n", K[0],
               K[0] - 800.0);
        printf("  fy : %8.2f  (800.0)  error %+.2f px\n", K[4],
               K[4] - 800.0);
        printf("  cx : %8.2f  (320.0)  error %+.2f px\n", K[2],
               K[2] - 320.0);
        printf("  cy : %8.2f  (240.0)  error %+.2f px\n", K[5],
               K[5] - 240.0);
        printf("  k1, k2 (true 0)     : %+.4f, %+.4f\n", k_radial[0],
               k_radial[1]);
        printf("  reprojection RMS    : %.3f px\n",
               mv_calib_reproj_rms(est, views, nviews));
    }
    return 0;
}
