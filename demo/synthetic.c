/* Synthetic two-camera experiment: the numbers this prints are quoted in
 * doc/multiview.tex (section "Basic results").
 *
 * Setup: two calibrated cameras with a 0.5 m baseline observe points in an
 * overlapping volume at 3.5-6 m depth. Projections are corrupted with
 * Gaussian pixel noise; we then estimate the fundamental matrix with the
 * normalized 8-point algorithm, triangulate every correspondence, rectify
 * the pair, and report accuracy against ground truth. Fully deterministic
 * (fixed-seed LCG). */

#include <math.h>
#include <stdio.h>

#include "mv/mv.h"

#define NPOINTS 250
#define NOISE_SIGMA 0.3 /* pixels, per coordinate */
#define IMG_W 640
#define IMG_H 480
#define MV_PI 3.14159265358979323846

static unsigned long long rng_state = 42;

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

static int in_bounds(const double uv[2])
{
    return uv[0] >= 0.0 && uv[0] < IMG_W && uv[1] >= 0.0 && uv[1] < IMG_H;
}

int main(void)
{
    mv_camera c1, c2, r1, r2;
    const mv_camera *cams[2];
    double C2pos[3] = { 0.5, 0.0, 0.0 };
    double Xtrue[NPOINTS][3];
    double uv1c[NPOINTS][2], uv2c[NPOINTS][2]; /* clean projections */
    double uv1n[2 * NPOINTS], uv2n[2 * NPOINTS]; /* noisy, kept points */
    double Fa[9], F8[9], H1[9], H2[9];
    double sum_dv_before = 0.0, sum_dv_after = 0.0;
    double se3d = 0.0, max3d = 0.0, se_reproj = 0.0;
    double frob = 0.0, mean_epi_a = 0.0, mean_epi_8 = 0.0;
    mv_cloud cloud;
    int i, kept = 0, tri_ok = 0;

    /* camera 1 at the origin looking down +z; camera 2 half a metre to the
     * right, yawed 6 degrees inward so the frusta overlap at 3.5-6 m */
    mv_cam_set_K(&c1, 800.0, 800.0, 320.0, 240.0);
    mv_cam_set_identity_pose(&c1);
    for (i = 0; i < 5; i++)
        c1.k[i] = 0.0;
    c2 = c1;
    mv_cam_set_pose_yaw(&c2, -6.0 * MV_PI / 180.0, C2pos);
    cams[0] = &c1;
    cams[1] = &c2;

    /* sample scene points visible in both views */
    for (i = 0; i < NPOINTS; i++) {
        double X[3], p1[2], p2[2];
        X[0] = -1.0 + 2.0 * urand();
        X[1] = -0.75 + 1.5 * urand();
        X[2] = 3.5 + 2.5 * urand();
        if (mv_cam_project(p1, &c1, X) != MV_OK || !in_bounds(p1))
            continue;
        if (mv_cam_project(p2, &c2, X) != MV_OK || !in_bounds(p2))
            continue;
        Xtrue[kept][0] = X[0]; Xtrue[kept][1] = X[1]; Xtrue[kept][2] = X[2];
        uv1c[kept][0] = p1[0]; uv1c[kept][1] = p1[1];
        uv2c[kept][0] = p2[0]; uv2c[kept][1] = p2[1];
        uv1n[2 * kept]     = p1[0] + NOISE_SIGMA * grand();
        uv1n[2 * kept + 1] = p1[1] + NOISE_SIGMA * grand();
        uv2n[2 * kept]     = p2[0] + NOISE_SIGMA * grand();
        uv2n[2 * kept + 1] = p2[1] + NOISE_SIGMA * grand();
        kept++;
    }

    printf("multiview synthetic experiment (seed 42)\n");
    printf("----------------------------------------\n");
    printf("baseline            : %.3f m\n", mv_baseline(&c1, &c2));
    printf("points sampled/kept : %d / %d\n", NPOINTS, kept);
    printf("pixel noise sigma   : %.2f px\n", NOISE_SIGMA);

    /* --- fundamental matrix: analytic vs estimated ------------------- */
    mv_fundamental_from_cams(Fa, &c1, &c2);
    if (mv_fundamental_8point(F8, uv1n, uv2n, kept) != MV_OK) {
        fprintf(stderr, "8-point estimation failed\n");
        return 1;
    }
    for (i = 0; i < 9; i++)
        frob += (Fa[i] - F8[i]) * (Fa[i] - F8[i]);
    frob = sqrt(frob);
    for (i = 0; i < kept; i++) {
        mean_epi_a += mv_sym_epipolar_dist(Fa, uv1n + 2 * i, uv2n + 2 * i);
        mean_epi_8 += mv_sym_epipolar_dist(F8, uv1n + 2 * i, uv2n + 2 * i);
    }
    mean_epi_a /= kept;
    mean_epi_8 /= kept;
    printf("\nfundamental matrix\n");
    printf("  ||F_analytic - F_8point||_F : %.3e (unit-norm F)\n", frob);
    printf("  mean sym epi dist, analytic : %.3f px\n", mean_epi_a);
    printf("  mean sym epi dist, 8-point  : %.3f px\n", mean_epi_8);

    /* --- triangulation ------------------------------------------------ */
    mv_cloud_init(&cloud, 0);
    for (i = 0; i < kept; i++) {
        double uv[4], X[3], d[3];
        double e;
        uv[0] = uv1n[2 * i]; uv[1] = uv1n[2 * i + 1];
        uv[2] = uv2n[2 * i]; uv[3] = uv2n[2 * i + 1];
        if (mv_triangulate(X, cams, uv, 2) != MV_OK)
            continue;
        d[0] = X[0] - Xtrue[i][0];
        d[1] = X[1] - Xtrue[i][1];
        d[2] = X[2] - Xtrue[i][2];
        e = mv_norm(d, 3);
        se3d += e * e;
        if (e > max3d)
            max3d = e;
        e = mv_reproj_rms(cams, uv, 2, X);
        se_reproj += e * e;
        mv_cloud_push(&cloud, X, NULL);
        tri_ok++;
    }
    printf("\ntriangulation (%d/%d points)\n", tri_ok, kept);
    printf("  RMS 3D error        : %.4f m\n", sqrt(se3d / tri_ok));
    printf("  max 3D error        : %.4f m\n", max3d);
    printf("  RMS reprojection    : %.3f px\n", sqrt(se_reproj / tri_ok));
    mv_cloud_write_ply("out_cloud.ply", &cloud);
    printf("  wrote out_cloud.ply : %d points\n", cloud.n);
    mv_cloud_free(&cloud);

    /* --- rectification ------------------------------------------------ */
    if (mv_rectify_pair(&c1, &c2, &r1, &r2, H1, H2) != MV_OK) {
        fprintf(stderr, "rectification failed\n");
        return 1;
    }
    for (i = 0; i < kept; i++) {
        double p1[2], p2[2];
        sum_dv_before += fabs(uv1c[i][1] - uv2c[i][1]);
        mv_cam_project(p1, &r1, Xtrue[i]);
        mv_cam_project(p2, &r2, Xtrue[i]);
        sum_dv_after += fabs(p1[1] - p2[1]);
    }
    printf("\nrectification (Fusiello et al.)\n");
    printf("  mean |v1 - v2| before : %.3f px\n", sum_dv_before / kept);
    printf("  mean |v1 - v2| after  : %.2e px\n", sum_dv_after / kept);
    return 0;
}
