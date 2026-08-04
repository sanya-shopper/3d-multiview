#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mv/mv.h"

#define MV_PI 3.14159265358979323846

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

static void make_pair(mv_camera *c1, mv_camera *c2)
{
    double C2pos[3] = { 0.5, 0.0, 0.0 };
    mv_cam_set_K(c1, 800.0, 800.0, 320.0, 240.0);
    mv_cam_set_identity_pose(c1);
    memset(c1->k, 0, sizeof(c1->k));
    *c2 = *c1;
    mv_cam_set_pose_yaw(c2, -6.0 * MV_PI / 180.0, C2pos);
}

/* deterministic scene points spread over the shared volume */
static void scene_point(double X[3], int i)
{
    X[0] = -0.8 + 0.17 * (i % 10);
    X[1] = -0.6 + 0.23 * (i % 6);
    X[2] = 3.6 + 0.19 * (i % 12);
}

static void test_inv3(void)
{
    double A[9] = { 4, 7, 2, 3, 6, 1, 2, 5, 3 };
    double Ai[9], I[9];
    double err = 0.0;
    int i;
    CHECK(mv_mat_inv3(Ai, A) == MV_OK, "inv3 succeeds");
    mv_mat_mul(I, A, Ai, 3, 3, 3);
    for (i = 0; i < 9; i++)
        err += fabs(I[i] - ((i % 4 == 0) ? 1.0 : 0.0));
    CHECK(err < 1e-12, "inv3: A * A^-1 = I");
}

static void test_svd(void)
{
    double A[15] = { 1, 2, 3, 4, 5, 6, 7, 8, 10, 2, 1, 0, 3, 3, 3 };
    double U[15], S[3], V[9], Vt[9], US[15], R[15];
    double err = 0.0, orth = 0.0;
    int i, j;
    memcpy(U, A, sizeof(A));
    CHECK(mv_svd(U, S, V, 5, 3) == MV_OK, "svd succeeds");
    CHECK(S[0] >= S[1] && S[1] >= S[2], "svd: singular values descending");
    for (i = 0; i < 5; i++)
        for (j = 0; j < 3; j++)
            US[i * 3 + j] = U[i * 3 + j] * S[j];
    mv_mat_transpose(Vt, V, 3, 3);
    mv_mat_mul(R, US, Vt, 5, 3, 3);
    for (i = 0; i < 15; i++)
        err += fabs(R[i] - A[i]);
    CHECK(err < 1e-10, "svd: U S V^T reconstructs A");
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++) {
            double d = 0.0;
            int k;
            for (k = 0; k < 5; k++)
                d += U[k * 3 + i] * U[k * 3 + j];
            orth += fabs(d - (i == j ? 1.0 : 0.0));
        }
    CHECK(orth < 1e-10, "svd: U columns orthonormal");
}

static void test_camera(void)
{
    mv_camera c1, c2;
    double C[3], orig[3], dir[3], uv[2], X[3] = { 0.3, -0.2, 4.0 };
    double px[3], d[3];
    double s;
    int i;
    make_pair(&c1, &c2);
    mv_cam_center(C, &c2);
    CHECK(fabs(C[0] - 0.5) + fabs(C[1]) + fabs(C[2]) < 1e-12,
          "camera center recovered from R,t");
    CHECK(mv_cam_project(uv, &c2, X) == MV_OK, "projection in front");
    CHECK(mv_cam_ray(orig, dir, &c2, uv) == MV_OK, "ray cast");
    /* X must lie on the ray: X = orig + s * dir */
    for (i = 0; i < 3; i++)
        d[i] = X[i] - orig[i];
    s = mv_dot(d, dir, 3);
    for (i = 0; i < 3; i++)
        px[i] = orig[i] + s * dir[i] - X[i];
    CHECK(mv_norm(px, 3) < 1e-9, "back-projected ray passes through X");
}

static void test_epipolar(void)
{
    mv_camera c1, c2;
    double Fa[9], F8[9];
    double uv1[2 * 24], uv2[2 * 24];
    double maxd = 0.0, frob = 0.0;
    int i, n = 24;
    make_pair(&c1, &c2);
    for (i = 0; i < n; i++) {
        double X[3], p1[2], p2[2];
        scene_point(X, i);
        mv_cam_project(p1, &c1, X);
        mv_cam_project(p2, &c2, X);
        uv1[2 * i] = p1[0]; uv1[2 * i + 1] = p1[1];
        uv2[2 * i] = p2[0]; uv2[2 * i + 1] = p2[1];
    }
    mv_fundamental_from_cams(Fa, &c1, &c2);
    for (i = 0; i < n; i++) {
        double d = mv_sym_epipolar_dist(Fa, uv1 + 2 * i, uv2 + 2 * i);
        if (d > maxd)
            maxd = d;
    }
    CHECK(maxd < 1e-8, "analytic F: exact points on epipolar lines");
    CHECK(mv_fundamental_8point(F8, uv1, uv2, n) == MV_OK,
          "8-point succeeds");
    for (i = 0; i < 9; i++)
        frob += (Fa[i] - F8[i]) * (Fa[i] - F8[i]);
    CHECK(sqrt(frob) < 1e-6, "8-point matches analytic F on exact data");
}

static void test_triangulate(void)
{
    mv_camera c1, c2;
    const mv_camera *cams[2];
    double X[3] = { -0.1, 0.25, 4.7 }, Xr[3], uv[4], d[3];
    make_pair(&c1, &c2);
    cams[0] = &c1;
    cams[1] = &c2;
    mv_cam_project(uv, &c1, X);
    mv_cam_project(uv + 2, &c2, X);
    CHECK(mv_triangulate(Xr, cams, uv, 2) == MV_OK, "triangulation succeeds");
    d[0] = Xr[0] - X[0]; d[1] = Xr[1] - X[1]; d[2] = Xr[2] - X[2];
    CHECK(mv_norm(d, 3) < 1e-8, "triangulation exact on noiseless data");
    CHECK(mv_reproj_rms(cams, uv, 2, Xr) < 1e-8, "reprojection ~ 0");
}

static void test_rectify(void)
{
    mv_camera c1, c2, r1, r2;
    double H1[9], H2[9], X[3], p1[2], p2[2];
    double maxdv = 0.0;
    int i;
    make_pair(&c1, &c2);
    CHECK(mv_rectify_pair(&c1, &c2, &r1, &r2, H1, H2) == MV_OK,
          "rectification succeeds");
    for (i = 0; i < 12; i++) {
        scene_point(X, i);
        mv_cam_project(p1, &r1, X);
        mv_cam_project(p2, &r2, X);
        if (fabs(p1[1] - p2[1]) > maxdv)
            maxdv = fabs(p1[1] - p2[1]);
    }
    CHECK(maxdv < 1e-8, "rectified rows aligned (v1 == v2)");
    CHECK(fabs(mv_baseline(&r1, &r2) - 0.5) < 1e-12,
          "rectification preserves baseline");
}

static void test_homography(void)
{
    /* exact homography recovery: plane at a known pose, no noise */
    mv_camera c;
    double C[3] = { 0.05, -0.02, 0.0 };
    double obj[2 * 12], uv[2 * 12], H[9];
    double maxerr = 0.0;
    int i;
    mv_cam_set_K(&c, 800.0, 805.0, 320.0, 240.0);
    mv_cam_set_pose_yaw(&c, 0.2, C);
    c.t[2] += 0.6;
    memset(c.k, 0, sizeof(c.k));
    mv_target_checkerboard(obj, 4, 3, 0.03);
    for (i = 0; i < 12; i++) {
        double X[3] = { obj[2 * i], obj[2 * i + 1], 0.0 };
        mv_cam_project(uv + 2 * i, &c, X);
    }
    CHECK(mv_homography_dlt(H, obj, uv, 12) == MV_OK, "homography succeeds");
    for (i = 0; i < 12; i++) {
        double w = H[6] * obj[2 * i] + H[7] * obj[2 * i + 1] + H[8];
        double u = (H[0] * obj[2 * i] + H[1] * obj[2 * i + 1] + H[2]) / w;
        double v = (H[3] * obj[2 * i] + H[4] * obj[2 * i + 1] + H[5]) / w;
        double e = fabs(u - uv[2 * i]) + fabs(v - uv[2 * i + 1]);
        if (e > maxerr)
            maxerr = e;
    }
    CHECK(maxerr < 1e-8, "homography exact on noiseless data");
}

static void test_calib(void)
{
    /* Zhang calibration recovers K and poses exactly from noiseless views */
    const double fx = 800.0, fy = 805.0, cxx = 322.0, cyy = 238.0;
    double obj[2 * MV_TARGET_LETTER_COLS * MV_TARGET_LETTER_ROWS];
    double img[5][2 * MV_TARGET_LETTER_COLS * MV_TARGET_LETTER_ROWS];
    const double angles[5][2] = {
        { 0.0, 0.0 }, { 0.25, 0.0 }, { -0.2, 0.15 },
        { 0.15, -0.3 }, { -0.15, -0.2 }
    };
    mv_calib_view views[5];
    mv_camera gt[5], est[5];
    double K[9], k_radial[2];
    double kerr, perr = 0.0;
    int n, v, i, j;

    n = mv_target_checkerboard(obj, MV_TARGET_LETTER_COLS,
                               MV_TARGET_LETTER_ROWS,
                               MV_TARGET_LETTER_SQUARE);
    for (v = 0; v < 5; v++) {
        double cax = cos(angles[v][0]), sax = sin(angles[v][0]);
        double cay = cos(angles[v][1]), say = sin(angles[v][1]);
        double Rx[9], Ry[9];
        Rx[0] = 1; Rx[1] = 0;   Rx[2] = 0;
        Rx[3] = 0; Rx[4] = cax; Rx[5] = -sax;
        Rx[6] = 0; Rx[7] = sax; Rx[8] = cax;
        Ry[0] = cay; Ry[1] = 0; Ry[2] = -say;
        Ry[3] = 0;   Ry[4] = 1; Ry[5] = 0;
        Ry[6] = say; Ry[7] = 0; Ry[8] = cay;
        mv_cam_set_K(&gt[v], fx, fy, cxx, cyy);
        mv_mat_mul(gt[v].R, Rx, Ry, 3, 3, 3);
        for (j = 0; j < 3; j++)
            gt[v].t[j] = -(gt[v].R[j * 3 + 0] * 0.072
                           + gt[v].R[j * 3 + 1] * 0.096);
        gt[v].t[2] += 0.6 + 0.05 * v;
        memset(gt[v].k, 0, sizeof(gt[v].k));
        for (i = 0; i < n; i++) {
            double X[3] = { obj[2 * i], obj[2 * i + 1], 0.0 };
            mv_cam_project(img[v] + 2 * i, &gt[v], X);
        }
        views[v].obj = obj;
        views[v].img = img[v];
        views[v].n = n;
    }
    CHECK(mv_calib_planar(K, est, k_radial, views, 5, 1) == MV_OK,
          "zhang calibration succeeds");
    kerr = fabs(K[0] - fx) + fabs(K[4] - fy) + fabs(K[2] - cxx)
         + fabs(K[5] - cyy);
    CHECK(kerr < 1e-5, "zhang: intrinsics exact on noiseless data");
    CHECK(fabs(k_radial[0]) + fabs(k_radial[1]) < 1e-8,
          "zhang: radial ~ 0 on undistorted data");
    for (v = 0; v < 5; v++) {
        for (j = 0; j < 9; j++)
            perr += fabs(gt[v].R[j] - est[v].R[j]);
        for (j = 0; j < 3; j++)
            perr += fabs(gt[v].t[j] - est[v].t[j]);
    }
    CHECK(perr < 1e-5, "zhang: poses exact on noiseless data");
    CHECK(mv_calib_reproj_rms(est, views, 5) < 1e-6,
          "zhang: reprojection ~ 0");
}

static void test_radial(void)
{
    /* distorted synthetic views: linear step recovers k1, k2 */
    double obj[2 * MV_TARGET_LETTER_COLS * MV_TARGET_LETTER_ROWS];
    double img[3][2 * MV_TARGET_LETTER_COLS * MV_TARGET_LETTER_ROWS];
    mv_calib_view views[3];
    mv_camera gt[3];
    double K[9], k_radial[2];
    int n, v, i, j;

    n = mv_target_checkerboard(obj, MV_TARGET_LETTER_COLS,
                               MV_TARGET_LETTER_ROWS,
                               MV_TARGET_LETTER_SQUARE);
    for (v = 0; v < 3; v++) {
        double C[3] = { 0.0, 0.0, 0.0 };
        mv_cam_set_K(&gt[v], 800.0, 800.0, 320.0, 240.0);
        mv_cam_set_pose_yaw(&gt[v], 0.15 * (v - 1), C);
        for (j = 0; j < 3; j++)
            gt[v].t[j] = -(gt[v].R[j * 3 + 0] * 0.072
                           + gt[v].R[j * 3 + 1] * 0.096);
        gt[v].t[2] += 0.55 + 0.1 * v;
        memset(gt[v].k, 0, sizeof(gt[v].k));
        gt[v].k[0] = -0.10; /* k1 */
        gt[v].k[1] = 0.02;  /* k2 */
        for (i = 0; i < n; i++) {
            double X[3] = { obj[2 * i], obj[2 * i + 1], 0.0 };
            mv_cam_project(img[v] + 2 * i, &gt[v], X);
        }
        views[v].obj = obj;
        views[v].img = img[v];
        views[v].n = n;
    }
    /* k1 and k2 are individually near-unidentifiable over the small radius
     * range a letter-page target covers; what must hold is that the fitted
     * distortion CURVE matches the data (small reprojection residual), and
     * that it genuinely models distortion (dropping it makes things much
     * worse) */
    {
        mv_camera est[3];
        double rms_with, rms_without;
        CHECK(mv_calib_planar(K, est, k_radial, views, 3, 1) == MV_OK,
              "calibration on distorted data succeeds");
        rms_with = mv_calib_reproj_rms(est, views, 3);
        CHECK(rms_with < 0.1, "radial: distortion curve fits (rms < 0.1 px)");
        for (v = 0; v < 3; v++)
            memset(est[v].k, 0, sizeof(est[v].k));
        rms_without = mv_calib_reproj_rms(est, views, 3);
        CHECK(rms_without > 2.0 * rms_with,
              "radial: modeled distortion is significant");
    }
}

static void test_calib_robust(void)
{
    /* A video-like view set: 36 near-fronto-parallel frames, 4 tilted
     * ones, one corrupted read, noise sigma = 0.3 px.  Plain Zhang is
     * swamped by the redundant frames (Sturm-Maybank degeneracy) and
     * poisoned by the corrupt view; the robust fit must recover K from
     * the informative minority and reject the corrupt view. */
    enum { NV = 40, TILT0 = 36, BAD = 20 };
    const double fx = 1600.0, fy = 1600.0, cxx = 540.0, cyy = 960.0;
    static double obj[2 * MV_TARGET_LETTER_COLS * MV_TARGET_LETTER_ROWS];
    static double img[NV][2 * MV_TARGET_LETTER_COLS * MV_TARGET_LETTER_ROWS];
    mv_calib_view views[NV];
    static mv_camera gt[NV], est[NV];
    unsigned char inl[NV];
    double K[9], k_radial[2];
    unsigned long long s = 20260804ULL;
    int n, v, i, j;

#define TCR_U() ((double)((s = s * 6364136223846793005ULL \
                               + 1442695040888963407ULL) >> 33) / 2147483648.0)
    n = mv_target_checkerboard(obj, MV_TARGET_LETTER_COLS,
                               MV_TARGET_LETTER_ROWS,
                               MV_TARGET_LETTER_SQUARE);
    for (v = 0; v < NV; v++) {
        double ax, ay;
        double cax, sax, cay, say, Rx[9], Ry[9];
        if (v < TILT0) { /* fronto-parallel to within ~1.5 degrees */
            ax = 0.025 * sin(0.7 * v);
            ay = 0.025 * cos(1.1 * v);
        } else { /* the informative minority: ~25 degree tilts */
            ax = (v & 1) ? 0.45 : -0.45;
            ay = (v & 2) ? 0.40 : -0.40;
        }
        cax = cos(ax); sax = sin(ax); cay = cos(ay); say = sin(ay);
        Rx[0] = 1; Rx[1] = 0;   Rx[2] = 0;
        Rx[3] = 0; Rx[4] = cax; Rx[5] = -sax;
        Rx[6] = 0; Rx[7] = sax; Rx[8] = cax;
        Ry[0] = cay; Ry[1] = 0; Ry[2] = -say;
        Ry[3] = 0;   Ry[4] = 1; Ry[5] = 0;
        Ry[6] = say; Ry[7] = 0; Ry[8] = cay;
        mv_cam_set_K(&gt[v], fx, fy, cxx, cyy);
        mv_mat_mul(gt[v].R, Rx, Ry, 3, 3, 3);
        for (j = 0; j < 3; j++)
            gt[v].t[j] = -(gt[v].R[j * 3 + 0] * 0.072
                           + gt[v].R[j * 3 + 1] * 0.096);
        gt[v].t[2] += 0.9 + 0.003 * v;
        memset(gt[v].k, 0, sizeof(gt[v].k));
        for (i = 0; i < n; i++) {
            double X[3] = { obj[2 * i], obj[2 * i + 1], 0.0 };
            double g0, g1;
            mv_cam_project(img[v] + 2 * i, &gt[v], X);
            /* gaussian-ish noise: sum of 12 uniforms - 6 */
            g0 = g1 = -6.0;
            for (j = 0; j < 12; j++) {
                g0 += TCR_U();
                g1 += TCR_U();
            }
            img[v][2 * i] += 0.3 * g0;
            img[v][2 * i + 1] += 0.3 * g1;
        }
        views[v].obj = obj;
        views[v].img = img[v];
        views[v].n = n;
    }
#undef TCR_U
    /* corrupt one fronto view by swapping adjacent point pairs.  NOTE:
     * reversing the whole point order would NOT corrupt it - a checker
     * grid is centro-symmetric, so reversal is a 180-degree rotation of
     * the plane, still a perfectly valid homography (the same trap as
     * the reader's rotation ambiguity).  Adjacent swaps are not a
     * projective map, so this is a genuinely bad decode. */
    for (i = 0; i + 1 < n; i += 2) {
        double t0 = img[BAD][2 * i], t1 = img[BAD][2 * i + 1];
        img[BAD][2 * i] = img[BAD][2 * (i + 1)];
        img[BAD][2 * i + 1] = img[BAD][2 * (i + 1) + 1];
        img[BAD][2 * (i + 1)] = t0;
        img[BAD][2 * (i + 1) + 1] = t1;
    }
    CHECK(mv_calib_planar_robust(K, est, k_radial, views, NV, 1, inl)
              == MV_OK, "robust calibration succeeds");
    CHECK(fabs(K[0] - fx) / fx < 0.02 && fabs(K[4] - fy) / fy < 0.02,
          "robust: focal within 2% on degenerate-heavy set");
    CHECK(fabs(K[2] - cxx) < 25.0 && fabs(K[5] - cyy) < 25.0,
          "robust: principal point within 25 px");
    CHECK(!inl[BAD], "robust: corrupted view rejected");
    {
        int nin = 0;
        for (v = 0; v < NV; v++)
            nin += inl[v];
        CHECK(nin >= 10, "robust: healthy views retained");
    }
    /* the discriminating half: plain Zhang on the same data must do
     * WORSE than robust, or this test guards nothing */
    {
        static mv_camera pest[NV];
        double Kp[9], kp[2], errp, errr;
        if (mv_calib_planar(Kp, pest, kp, views, NV, 1) == MV_OK) {
            errp = fabs(Kp[0] - fx) + fabs(Kp[4] - fy);
            errr = fabs(K[0] - fx) + fabs(K[4] - fy);
            CHECK(errp > errr,
                  "robust: beats plain Zhang on the degenerate-heavy set");
        } else {
            CHECK(1, "robust: plain Zhang fails outright here");
        }
    }
}

static void test_graycode(void)
{
    /* encode/decode roundtrip for a 1280-wide display axis */
    enum { W = 1280 };
    int nbits = mv_graycode_bits(W);
    unsigned char *frames[11], *invs[11];
    const unsigned char *cf[11], *ci[11];
    int coord[W];
    unsigned char valid[W];
    int b, x, ok = 1, nvalid;

    CHECK(nbits == 11, "graycode: 11 bits for 1280 columns");
    for (b = 0; b < nbits; b++) {
        frames[b] = (unsigned char *)malloc(W);
        invs[b] = (unsigned char *)malloc(W);
        mv_graycode_frame(frames[b], W, 1, 0, b, nbits, 0);
        mv_graycode_frame(invs[b], W, 1, 0, b, nbits, 1);
        cf[b] = frames[b];
        ci[b] = invs[b];
    }
    nvalid = mv_graycode_decode(coord, valid, cf, ci, nbits, W, 10);
    CHECK(nvalid == W, "graycode: all pixels decode as valid");
    for (x = 0; x < W; x++)
        if (coord[x] != x)
            ok = 0;
    CHECK(ok, "graycode: decoded coordinate equals pixel coordinate");
    for (b = 0; b < nbits; b++) {
        free(frames[b]);
        free(invs[b]);
    }
}

static unsigned test_rot90(unsigned code)
{
    unsigned out = 0;
    int r, c;
    for (r = 0; r < 4; r++)
        for (c = 0; c < 4; c++)
            out = (out << 1)
                | ((code >> (15 - ((3 - c) * 4 + r))) & 1u);
    return out;
}

static void test_pattern(void)
{
    static unsigned char pat[MV_PAT_W * MV_PAT_H];
    static unsigned char img[640 * 480];
    unsigned code = 0;
    int r, c, pc, pr, rot;

    CHECK(mv_pattern_selftest() == MV_OK, "pattern: M-array selftest");

    for (r = 0; r < 4; r++)
        for (c = 0; c < 4; c++)
            code = (code << 1) | (unsigned)mv_marray_bit(5 + c, 2 + r);
    CHECK(mv_marray_lookup(code, &pc, &pr, &rot) == MV_OK
          && pc == 5 && pr == 2 && rot == 0,
          "pattern: window lookup, unrotated");
    CHECK(mv_marray_lookup(test_rot90(code), &pc, &pr, &rot) == MV_OK
          && pc == 5 && pr == 2 && rot != 0,
          "pattern: window lookup under rotation");

    /* end-to-end: render frontal view, read blind, verify exactly */
    {
        mv_camera cam;
        mv_read_result rr;
        unsigned long long seed = 1;
        double center[3] = { MV_PAT_W / 2.0 * 0.0002745,
                             MV_PAT_H / 2.0 * 0.0002745, 0.0 };
        int i, id_bad = 0;
        mv_cam_set_K(&cam, 800.0, 800.0, 320.0, 240.0);
        mv_cam_set_identity_pose(&cam);
        memset(cam.k, 0, sizeof(cam.k));
        for (i = 0; i < 3; i++)
            cam.t[i] = -center[i];
        cam.t[2] += 0.85;
        mv_pattern_render(pat, 777777u);
        CHECK(mv_render_plane(img, 640, 480, &cam, pat, MV_PAT_W,
                              MV_PAT_H, 0.0002745, 128, 0.0, &seed)
              == MV_OK, "pattern: render succeeds");
        CHECK(mv_read_pattern(&rr, img, 640, 480) == MV_OK,
              "pattern: blind read succeeds");
        CHECK(rr.n == 162, "pattern: all 162 corners identified");
        {
            double se = 0.0;
            for (i = 0; i < rr.n; i++) {
                double xy[2], X[3], uv[2], du, dv, e2;
                mv_pattern_corner_px(rr.id[i] % MV_PAT_CORNER_COLS,
                                     rr.id[i] / MV_PAT_CORNER_COLS, xy);
                X[0] = xy[0] * 0.0002745;
                X[1] = xy[1] * 0.0002745;
                X[2] = 0.0;
                mv_cam_project(uv, &cam, X);
                du = rr.uv[2 * i] - uv[0];
                dv = rr.uv[2 * i + 1] - uv[1];
                e2 = du * du + dv * dv;
                se += e2;
                if (e2 > 9.0)
                    id_bad++;
            }
            CHECK(id_bad == 0, "pattern: every corner id correct");
            CHECK(sqrt(se / rr.n) < 0.5,
                  "pattern: localization RMS < 0.5 px");
        }
        CHECK(rr.counter_valid && rr.counter == 777777u,
              "pattern: counter decodes exactly");
    }
}

static void test_tsdf(void)
{
    mv_camera c1, c2;
    const mv_camera *cams[2];
    double C2pos[3] = { 0.5, 0.0, 0.0 };
    double ctr[3] = { 0.0, 0.1, 4.7 };
    mv_tsdf t;
    double *tris;
    int ntri, ci, x, y, i;

    mv_cam_set_K(&c1, 800.0, 800.0, 320.0, 240.0);
    mv_cam_set_identity_pose(&c1);
    memset(c1.k, 0, sizeof(c1.k));
    c2 = c1;
    mv_cam_set_pose_yaw(&c2, -6.0 * MV_PI / 180.0, C2pos);
    cams[0] = &c1;
    cams[1] = &c2;
    CHECK(mv_tsdf_init(&t, ctr[0] - 0.45, ctr[1] - 0.45, ctr[2] - 0.45,
                       0.9, 0.9, 0.9, 0.01, 0.05) == MV_OK,
          "tsdf: init");
    /* fuse noiseless depth of a r=0.3 sphere at ctr */
    for (ci = 0; ci < 2; ci++) {
        double O[3];
        mv_cam_center(O, cams[ci]);
        for (y = 0; y < 480; y++)
            for (x = 0; x < 640; x++) {
                double uv[2] = { (double)x, (double)y };
                double orig[3], dir[3], oc[3], b, c, disc, s, p[3];
                if (mv_cam_ray(orig, dir, cams[ci], uv) != MV_OK)
                    continue;
                for (i = 0; i < 3; i++)
                    oc[i] = orig[i] - ctr[i];
                b = mv_dot(oc, dir, 3);
                c = mv_dot(oc, oc, 3) - 0.09;
                disc = b * b - c;
                if (disc < 0.0)
                    continue;
                s = -b - sqrt(disc);
                for (i = 0; i < 3; i++)
                    p[i] = orig[i] + s * dir[i];
                mv_tsdf_fuse(&t, p, O, 1.0);
            }
    }
    /* signed queries just outside / inside the front surface */
    {
        double qo[3] = { 0.0, 0.1, 4.7 - 0.3 - 0.03 };
        double qi[3] = { 0.0, 0.1, 4.7 - 0.3 + 0.03 };
        double vo = mv_tsdf_query(&t, qo), vi = mv_tsdf_query(&t, qi);
        CHECK(vo < HUGE_VAL && vo > 0.01 && vo < 0.05,
              "tsdf: positive outside the surface");
        CHECK(vi < HUGE_VAL && vi < -0.01 && vi > -0.05,
              "tsdf: negative inside the surface");
    }
    {
        double far_q[3] = { 0.4, 0.4, 4.7 }; /* off the sphere band */
        CHECK(mv_tsdf_query(&t, far_q) == HUGE_VAL,
              "tsdf: unobserved voxels report unknown");
    }
    CHECK(mv_tsdf_mesh(&t, &tris, &ntri) == MV_OK && ntri > 1000,
          "tsdf: mesh extraction");
    {
        double se = 0.0;
        for (i = 0; i < 3 * ntri; i++) {
            double d[3], e;
            int a;
            for (a = 0; a < 3; a++)
                d[a] = tris[3 * i + a] - ctr[a];
            e = mv_norm(d, 3) - 0.3;
            se += e * e;
        }
        CHECK(sqrt(se / (3 * ntri)) < 0.003,
              "tsdf: noiseless surface RMS < 3 mm");
    }
    free(tris);
    mv_tsdf_free(&t);
}

static void test_warp_scene(void)
{
    /* warp: identity homography reproduces the image exactly */
    static unsigned char src[64 * 48], dst[64 * 48];
    double Hid[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    int i, ok = 1;
    for (i = 0; i < 64 * 48; i++)
        src[i] = (unsigned char)((i * 37) & 0xff);
    CHECK(mv_warp_homography(dst, 64, 48, src, 64, 48, Hid, 0) == MV_OK,
          "warp: runs");
    for (i = 0; i < 64 * 48; i++)
        if (dst[i] != src[i] && i % 64 < 63 && i / 64 < 47)
            ok = 0;
    CHECK(ok, "warp: identity homography is exact");

    /* scene renderer: frontal plane at Z=2 gives depth 2 at center */
    {
        static unsigned char tex[32 * 32], img[64 * 48];
        static float zc[64 * 48];
        mv_camera cam;
        mv_plane pl;
        static const double ex[3] = { 1, 0, 0 }, ey[3] = { 0, 1, 0 };
        double r3[3], t0[3] = { -0.32, -0.24, 2.0 };
        int j;
        for (i = 0; i < 32 * 32; i++)
            tex[i] = (unsigned char)(i & 0xff);
        mv_cam_set_K(&cam, 100.0, 100.0, 32.0, 24.0);
        mv_cam_set_identity_pose(&cam);
        memset(cam.k, 0, sizeof(cam.k));
        for (j = 0; j < 3; j++) {
            pl.R[j * 3 + 0] = ex[j];
            pl.R[j * 3 + 1] = ey[j];
            pl.t[j] = t0[j];
        }
        mv_cross3(r3, ex, ey);
        for (j = 0; j < 3; j++)
            pl.R[j * 3 + 2] = r3[j];
        pl.tex = tex;
        pl.tw = 32;
        pl.th = 32;
        pl.pitch = 0.02;
        {
            unsigned long long seed = 5;
            CHECK(mv_render_scene(img, zc, 64, 48, &cam, &pl, 1, 7, 0.0,
                                  &seed) == MV_OK, "scene render: runs");
        }
        CHECK(fabs(zc[24 * 64 + 32] - 2.0) < 1e-9,
              "scene render: center depth exact");
        CHECK(zc[0] == HUGE_VALF || zc[0] == 2.0f,
              "scene render: miss marked or hit");
    }
}

int main(void)
{
    test_inv3();
    test_svd();
    test_camera();
    test_epipolar();
    test_triangulate();
    test_rectify();
    test_homography();
    test_calib();
    test_radial();
    test_calib_robust();
    test_graycode();
    test_pattern();
    test_tsdf();
    test_warp_scene();
    if (failures) {
        printf("\n%d FAILURE(S)\n", failures);
        return 1;
    }
    printf("\nall tests passed\n");
    return 0;
}
