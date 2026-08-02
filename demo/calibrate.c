/* Synthetic calibration experiment: the numbers this prints are quoted in
 * doc/multiview.tex (section "Basic results").
 *
 * A camera with known ground-truth intrinsics photographs the letter-page
 * checkerboard target (7x9 inner corners, 24 mm squares) in six poses.
 * Corner projections are corrupted with Gaussian pixel noise; Zhang
 * calibration then recovers the intrinsics and per-view poses, which are
 * compared against ground truth. Also renders the printable target to
 * target_letter.pgm. Fully deterministic (fixed-seed LCG). */

#include <math.h>
#include <stdio.h>

#include "mv/mv.h"

#define NVIEWS 6
#define NOISE_SIGMA 0.3
#define MV_PI 3.14159265358979323846

static unsigned long long rng_state = 7;

static double urand(void)
{
    rng_state = rng_state * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)((rng_state >> 11) & ((1ULL << 53) - 1)) / 9007199254740992.0;
}

static double grand(void)
{
    double u1 = urand(), u2 = urand();
    if (u1 < 1e-300)
        u1 = 1e-300;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * MV_PI * u2);
}

/* R = Rx(pitch) * Ry(yaw) */
static void rot_pitch_yaw(double R[9], double pitch, double yaw)
{
    double cx = cos(pitch), sx = sin(pitch);
    double cy = cos(yaw), sy = sin(yaw);
    double Rx[9], Ry[9];
    Rx[0] = 1; Rx[1] = 0;  Rx[2] = 0;
    Rx[3] = 0; Rx[4] = cx; Rx[5] = -sx;
    Rx[6] = 0; Rx[7] = sx; Rx[8] = cx;
    Ry[0] = cy; Ry[1] = 0; Ry[2] = -sy;
    Ry[3] = 0;  Ry[4] = 1; Ry[5] = 0;
    Ry[6] = sy; Ry[7] = 0; Ry[8] = cy;
    mv_mat_mul(R, Rx, Ry, 3, 3, 3);
}

static double rot_angle_deg(const double Ra[9], const double Rb[9])
{
    double Rbt[9], D[9], tr, c;
    mv_mat_transpose(Rbt, Rb, 3, 3);
    mv_mat_mul(D, Ra, Rbt, 3, 3, 3);
    tr = D[0] + D[4] + D[8];
    c = 0.5 * (tr - 1.0);
    if (c > 1.0)
        c = 1.0;
    if (c < -1.0)
        c = -1.0;
    return acos(c) * 180.0 / MV_PI;
}

int main(void)
{
    /* ground-truth intrinsics (zero skew, slightly anisotropic) */
    const double fx = 800.0, fy = 805.0, cx = 322.0, cy = 238.0;
    const double pose_deg[NVIEWS][2] = {
        { 0.0, 0.0 }, { 30.0, 0.0 }, { -30.0, 15.0 },
        { 20.0, -35.0 }, { -20.0, -25.0 }, { 35.0, 25.0 }
    };
    double obj[2 * MV_TARGET_LETTER_COLS * MV_TARGET_LETTER_ROWS];
    double img[NVIEWS][2 * MV_TARGET_LETTER_COLS * MV_TARGET_LETTER_ROWS];
    mv_calib_view views[NVIEWS];
    mv_camera gt[NVIEWS], est[NVIEWS];
    double K[9], k_radial[2];
    double center[3], rot_err = 0.0, t_err = 0.0;
    int ncorners, v, i;

    ncorners = mv_target_checkerboard(obj, MV_TARGET_LETTER_COLS,
                                      MV_TARGET_LETTER_ROWS,
                                      MV_TARGET_LETTER_SQUARE);
    center[0] = 0.5 * (MV_TARGET_LETTER_COLS - 1) * MV_TARGET_LETTER_SQUARE;
    center[1] = 0.5 * (MV_TARGET_LETTER_ROWS - 1) * MV_TARGET_LETTER_SQUARE;
    center[2] = 0.0;

    printf("multiview calibration experiment (seed 7)\n");
    printf("-----------------------------------------\n");
    printf("target              : letter page, %dx%d inner corners, "
           "%.0f mm squares\n", MV_TARGET_LETTER_COLS, MV_TARGET_LETTER_ROWS,
           MV_TARGET_LETTER_SQUARE * 1000.0);
    printf("views / noise sigma : %d / %.2f px\n", NVIEWS, NOISE_SIGMA);

    for (v = 0; v < NVIEWS; v++) {
        double R[9];
        int j;
        rot_pitch_yaw(R, pose_deg[v][0] * MV_PI / 180.0,
                      pose_deg[v][1] * MV_PI / 180.0);
        mv_cam_set_K(&gt[v], fx, fy, cx, cy);
        for (j = 0; j < 9; j++)
            gt[v].R[j] = R[j];
        /* place the board center in front of the camera at 0.55-0.8 m */
        for (j = 0; j < 3; j++)
            gt[v].t[j] = -(R[j * 3 + 0] * center[0]
                           + R[j * 3 + 1] * center[1]);
        gt[v].t[0] += -0.03 + 0.06 * urand();
        gt[v].t[1] += -0.03 + 0.06 * urand();
        gt[v].t[2] += 0.45 + 0.20 * urand();
        for (j = 0; j < 5; j++)
            gt[v].k[j] = 0.0;

        for (i = 0; i < ncorners; i++) {
            double X[3], p[2];
            X[0] = obj[2 * i];
            X[1] = obj[2 * i + 1];
            X[2] = 0.0;
            if (mv_cam_project(p, &gt[v], X) != MV_OK) {
                fprintf(stderr, "corner behind camera\n");
                return 1;
            }
            img[v][2 * i] = p[0] + NOISE_SIGMA * grand();
            img[v][2 * i + 1] = p[1] + NOISE_SIGMA * grand();
        }
        views[v].obj = obj;
        views[v].img = img[v];
        views[v].n = ncorners;
    }

    if (mv_calib_planar(K, est, k_radial, views, NVIEWS, 1) != MV_OK) {
        fprintf(stderr, "calibration failed\n");
        return 1;
    }

    printf("\nintrinsics (recovered vs true)\n");
    printf("  fx : %8.2f  (%.1f)   error %+.2f px\n", K[0], fx, K[0] - fx);
    printf("  fy : %8.2f  (%.1f)   error %+.2f px\n", K[4], fy, K[4] - fy);
    printf("  cx : %8.2f  (%.1f)   error %+.2f px\n", K[2], cx, K[2] - cx);
    printf("  cy : %8.2f  (%.1f)   error %+.2f px\n", K[5], cy, K[5] - cy);
    printf("  k1, k2 (true 0)     : %+.4f, %+.4f\n", k_radial[0], k_radial[1]);

    for (v = 0; v < NVIEWS; v++) {
        double dt[3];
        int j;
        rot_err += rot_angle_deg(gt[v].R, est[v].R);
        for (j = 0; j < 3; j++)
            dt[j] = gt[v].t[j] - est[v].t[j];
        t_err += mv_norm(dt, 3);
    }
    printf("\nextrinsics (mean over %d views)\n", NVIEWS);
    printf("  rotation error      : %.3f deg\n", rot_err / NVIEWS);
    printf("  translation error   : %.2f mm\n", 1000.0 * t_err / NVIEWS);
    printf("\nreprojection RMS      : %.3f px\n",
           mv_calib_reproj_rms(est, views, NVIEWS));

    if (mv_target_render_pgm("target_letter.pgm", MV_TARGET_LETTER_COLS,
                             MV_TARGET_LETTER_ROWS, 240, 60) == MV_OK)
        printf("\nwrote target_letter.pgm (print at 254 dpi = 10 px/mm "
               "for 24 mm squares)\n");
    return 0;
}
