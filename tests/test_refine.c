#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mv/mv.h"

static int failures = 0;

#define CHECK(cond, name) \
    do { \
        if (cond) \
            printf("ok:   %s\n", name); \
        else { \
            printf("FAIL: %s\n", name); \
            failures++; \
        } \
    } while (0)

/* ---- deterministic noise (project-standard LCG) ------------------ */

static unsigned long long rng_state;

static double urand(void)
{
    rng_state = rng_state * 6364136223846793005ULL
              + 1442695040888963407ULL;
    return (double)(rng_state >> 33) / 2147483648.0;
}

/* gaussian-ish: sum of 12 uniforms - 6 (unit variance) */
static double grand12(void)
{
    double g = -6.0;
    int i;
    for (i = 0; i < 12; i++)
        g += urand();
    return g;
}

/* ---- shared synthetic geometry ----------------------------------- */

/* 6 views of the letter checkerboard, tilts up to +-0.3 rad, distances
 * ~0.6 m; ground truth fx=800 fy=805 cx=322 cy=238, k1=-0.10 k2=+0.02
 * applied in projection.  View 0 is fronto-parallel (identity rotation),
 * so the Rodrigues theta ~ 0 branch is exercised by every run. */
enum { NV = 6, NPTS = MV_TARGET_LETTER_COLS * MV_TARGET_LETTER_ROWS };
#define GT_FX 800.0
#define GT_FY 805.0
#define GT_CX 322.0
#define GT_CY 238.0
#define GT_K1 -0.10
#define GT_K2 0.02

static double g_obj[2 * NPTS];
static double g_img[NV][2 * NPTS];

static void make_views(mv_calib_view *views, mv_camera *gt, double sigma,
                       unsigned long long seed)
{
    static const double ang[NV][2] = {
        { 0.0, 0.0 },   { 0.3, 0.0 },   { -0.3, 0.15 },
        { 0.2, -0.3 },  { -0.25, -0.2 }, { 0.3, 0.25 }
    };
    int n, v, i, j;

    rng_state = seed;
    n = mv_target_checkerboard(g_obj, MV_TARGET_LETTER_COLS,
                               MV_TARGET_LETTER_ROWS,
                               MV_TARGET_LETTER_SQUARE);
    for (v = 0; v < NV; v++) {
        double cax = cos(ang[v][0]), sax = sin(ang[v][0]);
        double cay = cos(ang[v][1]), say = sin(ang[v][1]);
        double Rx[9], Ry[9];
        Rx[0] = 1; Rx[1] = 0;   Rx[2] = 0;
        Rx[3] = 0; Rx[4] = cax; Rx[5] = -sax;
        Rx[6] = 0; Rx[7] = sax; Rx[8] = cax;
        Ry[0] = cay; Ry[1] = 0; Ry[2] = -say;
        Ry[3] = 0;   Ry[4] = 1; Ry[5] = 0;
        Ry[6] = say; Ry[7] = 0; Ry[8] = cay;
        mv_cam_set_K(&gt[v], GT_FX, GT_FY, GT_CX, GT_CY);
        mv_mat_mul(gt[v].R, Rx, Ry, 3, 3, 3);
        /* aim at the target centre (0.072, 0.096) as in test_calib */
        for (j = 0; j < 3; j++)
            gt[v].t[j] = -(gt[v].R[j * 3 + 0] * 0.072
                           + gt[v].R[j * 3 + 1] * 0.096);
        gt[v].t[2] += 0.6 + 0.05 * v;
        memset(gt[v].k, 0, sizeof(gt[v].k));
        gt[v].k[0] = GT_K1;
        gt[v].k[1] = GT_K2;
        for (i = 0; i < n; i++) {
            double X[3] = { g_obj[2 * i], g_obj[2 * i + 1], 0.0 };
            mv_cam_project(g_img[v] + 2 * i, &gt[v], X);
            g_img[v][2 * i] += sigma * grand12();
            g_img[v][2 * i + 1] += sigma * grand12();
        }
        views[v].obj = g_obj;
        views[v].img = g_img[v];
        views[v].n = n;
    }
}

/* ---- 1. sigma = 0 exactness -------------------------------------- */

static void test_exact(void)
{
    mv_calib_view views[NV];
    mv_camera gt[NV], est[NV];
    double K[9], k_radial[2], kerr, perr = 0.0;
    int v, j;

    make_views(views, gt, 0.0, 1ULL); /* seed unused at sigma = 0 */
    CHECK(mv_calib_planar(K, est, k_radial, views, NV, 1) == MV_OK,
          "exact: alternation initializer succeeds");
    CHECK(mv_calib_refine(K, est, k_radial, views, NV) == MV_OK,
          "exact: refinement succeeds");
    CHECK(mv_calib_reproj_rms(est, views, NV) < 1e-8,
          "exact: refined reprojection rms < 1e-8");
    kerr = fabs(K[0] - GT_FX) + fabs(K[4] - GT_FY) + fabs(K[2] - GT_CX)
         + fabs(K[5] - GT_CY);
    CHECK(kerr < 1e-5, "exact: intrinsics within 1e-5 of truth");
    CHECK(fabs(k_radial[0] - GT_K1) + fabs(k_radial[1] - GT_K2) < 1e-6,
          "exact: k1 k2 within 1e-6 of truth");
    /* pose recovery doubles as the Rodrigues roundtrip exactness test
     * for the generic and theta ~ 0 branches (view 0 is identity) */
    for (v = 0; v < NV; v++) {
        for (j = 0; j < 9; j++)
            perr += fabs(gt[v].R[j] - est[v].R[j]);
        for (j = 0; j < 3; j++)
            perr += fabs(gt[v].t[j] - est[v].t[j]);
    }
    CHECK(perr < 1e-5, "exact: poses within 1e-5 of truth");
}

/* ---- 2. Rodrigues near-pi branch --------------------------------- */

/* A single view whose world-to-camera rotation is within 1e-7 of a
 * half-turn: refinement is started AT the exact solution, so any error
 * in the mat->rod->mat roundtrip through the theta ~ pi branch would
 * show up directly as a nonzero residual and a moved rotation. */
static void test_rodrigues_pi(void)
{
    const double th = 3.14159255358979; /* pi - 1e-7 */
    mv_calib_view vw;
    mv_camera cam;
    double img[2 * NPTS], K[9], k_radial[2], rerr = 0.0;
    int n, i, j;

    n = mv_target_checkerboard(g_obj, MV_TARGET_LETTER_COLS,
                               MV_TARGET_LETTER_ROWS,
                               MV_TARGET_LETTER_SQUARE);
    mv_cam_set_K(&cam, GT_FX, GT_FY, GT_CX, GT_CY);
    cam.R[0] = 1; cam.R[1] = 0;        cam.R[2] = 0;
    cam.R[3] = 0; cam.R[4] = cos(th);  cam.R[5] = -sin(th);
    cam.R[6] = 0; cam.R[7] = sin(th);  cam.R[8] = cos(th);
    cam.t[0] = -0.072;
    cam.t[1] = 0.096;
    cam.t[2] = 0.6;
    memset(cam.k, 0, sizeof(cam.k));
    cam.k[0] = GT_K1;
    cam.k[1] = GT_K2;
    for (i = 0; i < n; i++) {
        double X[3] = { g_obj[2 * i], g_obj[2 * i + 1], 0.0 };
        mv_cam_project(img + 2 * i, &cam, X);
    }
    vw.obj = g_obj;
    vw.img = img;
    vw.n = n;
    memcpy(K, cam.K, sizeof(K));
    k_radial[0] = GT_K1;
    k_radial[1] = GT_K2;
    {
        mv_camera est = cam;
        CHECK(mv_calib_refine(K, &est, k_radial, &vw, 1) == MV_OK,
              "rodrigues pi: refinement succeeds");
        CHECK(mv_calib_reproj_rms(&est, &vw, 1) < 1e-9,
              "rodrigues pi: residual stays at zero");
        for (j = 0; j < 9; j++)
            rerr += fabs(est.R[j] - cam.R[j]);
        CHECK(rerr < 1e-7, "rodrigues pi: rotation roundtrip exact");
    }
}

/* ---- 3. T5 acceptance + never-worsens ---------------------------- */

static void test_t5(void)
{
    /* Pre-registered acceptance (doc/multiview.tex, tripwire T5): the
     * alternation estimator leaves RMS/predicted ~ 1.17; joint ML must
     * bring the ensemble mean to 1.00 +- 0.03.  predicted per-point
     * distance RMS = sqrt(2) sigma sqrt(1 - p/n) with p = 6 + 6*NV
     * fitted parameters and n = 2*NV*NPTS scalar residuals. */
    enum { NSEEDS = 10 };
    const double sigma = 0.3;
    const int pfit = 6 + 6 * NV;
    const int nres = 2 * NV * NPTS;
    const double predicted =
        sqrt(2.0) * sigma * sqrt(1.0 - (double)pfit / (double)nres);
    mv_calib_view views[NV];
    mv_camera gt[NV], est[NV];
    double K[9], k_radial[2];
    double sum_alt = 0.0, sum_ref = 0.0;
    int s, ok_all = 1, never_worse = 1;

    for (s = 0; s < NSEEDS; s++) {
        double rms_alt, rms_ref;
        make_views(views, gt, sigma, 20260805ULL + 1000ULL * (unsigned)s);
        if (mv_calib_planar(K, est, k_radial, views, NV, 1) != MV_OK
            || (rms_alt = mv_calib_reproj_rms(est, views, NV),
                mv_calib_refine(K, est, k_radial, views, NV) != MV_OK)) {
            ok_all = 0;
            break;
        }
        rms_ref = mv_calib_reproj_rms(est, views, NV);
        if (rms_ref > rms_alt)
            never_worse = 0;
        sum_alt += rms_alt / predicted;
        sum_ref += rms_ref / predicted;
    }
    CHECK(ok_all, "t5: all ensemble members calibrate and refine");
    if (ok_all) {
        double mean_alt = sum_alt / NSEEDS, mean_ref = sum_ref / NSEEDS;
        printf("t5: mean ratio alternation %.4f  refined %.4f "
               "(predicted rms %.4f px)\n", mean_alt, mean_ref, predicted);
        CHECK(mean_ref >= 0.97 && mean_ref <= 1.03,
              "t5: refined ratio within 1.00 +- 0.03");
        CHECK(mean_ref < mean_alt,
              "t5: refinement beats alternation on the ensemble mean");
        CHECK(never_worse,
              "t5: refinement never worsens rms on any member");
    }
}

int main(void)
{
    test_exact();
    test_rodrigues_pi();
    test_t5();
    if (failures) {
        printf("\n%d FAILURE(S)\n", failures);
        return 1;
    }
    printf("\nall refine tests passed\n");
    return 0;
}
