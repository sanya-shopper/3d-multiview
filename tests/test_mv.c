#include <math.h>
#include <stdio.h>
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

int main(void)
{
    test_inv3();
    test_svd();
    test_camera();
    test_epipolar();
    test_triangulate();
    test_rectify();
    if (failures) {
        printf("\n%d FAILURE(S)\n", failures);
        return 1;
    }
    printf("\nall tests passed\n");
    return 0;
}
