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

/* deterministic LCG; gaussian-ish noise = sum of 12 uniforms - 6 */
static unsigned long long rs = 1ULL;

static double urand(void)
{
    rs = rs * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)(rs >> 33) / 2147483648.0;
}

static double grand_(void)
{
    double g = -6.0;
    int i;
    for (i = 0; i < 12; i++)
        g += urand();
    return g;
}

/* same rig as tests/test_mv.c: fx 800, 0.5 m baseline, 6 deg vergence */
static void make_pair(mv_camera *c1, mv_camera *c2)
{
    double C2pos[3] = { 0.5, 0.0, 0.0 };
    mv_cam_set_K(c1, 800.0, 800.0, 320.0, 240.0);
    mv_cam_set_identity_pose(c1);
    memset(c1->k, 0, sizeof(c1->k));
    *c2 = *c1;
    mv_cam_set_pose_yaw(c2, -6.0 * MV_PI / 180.0, C2pos);
}

/* deterministic scene points in the shared volume, depth 3.5-4.5 m */
static void scene_point(double X[3], int i)
{
    X[0] = -0.8 + 0.17 * (i % 10);
    X[1] = -0.5 + 0.21 * (i % 6);
    X[2] = 3.5 + 0.09 * (i % 12);
}

static double dist3(const double a[3], const double b[3])
{
    double d[3];
    d[0] = a[0] - b[0];
    d[1] = a[1] - b[1];
    d[2] = a[2] - b[2];
    return mv_norm(d, 3);
}

static double dist2sq(const double a[2], const double b[2])
{
    double du = a[0] - b[0], dv = a[1] - b[1];
    return du * du + dv * dv;
}

/* --- (A) optimal triangulation ------------------------------------------- */

static void test_optimal_exact(void)
{
    mv_camera c1, c2;
    double F[9];
    double maxe = 0.0, maxcorr = 0.0;
    int i, ok = 0;
    make_pair(&c1, &c2);
    mv_fundamental_from_cams(F, &c1, &c2);
    for (i = 0; i < 60; i++) {
        double X[3], p1[2], p2[2], Xr[3], o1[2], o2[2], e, m;
        scene_point(X, i);
        if (mv_cam_project(p1, &c1, X) != MV_OK ||
            mv_cam_project(p2, &c2, X) != MV_OK)
            continue;
        if (mv_triangulate_optimal(Xr, &c1, &c2, p1, p2) != MV_OK)
            continue;
        e = dist3(Xr, X);
        if (e > maxe)
            maxe = e;
        if (mv_epipolar_correct(o1, o2, F, p1, p2) != MV_OK)
            continue;
        m = sqrt(dist2sq(o1, p1)) + sqrt(dist2sq(o2, p2));
        if (m > maxcorr)
            maxcorr = m;
        ok++;
    }
    CHECK(ok == 60, "sigma=0: all points triangulated and corrected");
    CHECK(maxe < 1e-8, "sigma=0: optimal triangulation exact (< 1e-8)");
    CHECK(maxcorr < 1e-7, "sigma=0: exact pairs are left unchanged");
}

static void test_optimal_noisy(void)
{
    mv_camera c1, c2;
    const mv_camera *cams[2];
    double F[9];
    double sum_opt = 0.0, sum_plain = 0.0, maxep = 0.0;
    double mean_opt, mean_plain;
    int i, n = 0, cost_violations = 0, n_improved = 0;
    const double sigma = 1.0;
    char line[128];
    make_pair(&c1, &c2);
    cams[0] = &c1;
    cams[1] = &c2;
    mv_fundamental_from_cams(F, &c1, &c2);
    rs = 20260805ULL;
    for (i = 0; i < 300; i++) {
        double X[3], p1[2], p2[2], m[4], o1[2], o2[2];
        double Xo[3], Xp[3], q1[2], q2[2], d, cost_opt, cost_dlt;
        scene_point(X, i);
        if (mv_cam_project(p1, &c1, X) != MV_OK ||
            mv_cam_project(p2, &c2, X) != MV_OK)
            continue;
        m[0] = p1[0] + sigma * grand_();
        m[1] = p1[1] + sigma * grand_();
        m[2] = p2[0] + sigma * grand_();
        m[3] = p2[1] + sigma * grand_();
        if (mv_epipolar_correct(o1, o2, F, m, m + 2) != MV_OK)
            continue;
        d = mv_sym_epipolar_dist(F, o1, o2);
        if (d > maxep)
            maxep = d;
        if (mv_triangulate_optimal(Xo, &c1, &c2, m, m + 2) != MV_OK ||
            mv_triangulate(Xp, cams, m, 2) != MV_OK)
            continue;
        sum_opt += dist3(Xo, X);
        sum_plain += dist3(Xp, X);
        if (dist3(Xo, X) <= dist3(Xp, X))
            n_improved++;
        /* optimality: reprojecting the plain-DLT point gives another
         * exactly consistent pair; the HS pair must cost no more */
        cost_opt = dist2sq(o1, m) + dist2sq(o2, m + 2);
        if (mv_cam_project(q1, &c1, Xp) == MV_OK &&
            mv_cam_project(q2, &c2, Xp) == MV_OK) {
            cost_dlt = dist2sq(q1, m) + dist2sq(q2, m + 2);
            if (cost_opt > cost_dlt + 1e-9)
                cost_violations++;
        }
        n++;
    }
    mean_opt = sum_opt / n;
    mean_plain = sum_plain / n;
    printf("  noisy n=%d sigma=%.1f px: mean 3-D err optimal %.3f mm, "
           "plain %.3f mm (ratio %.5f), optimal <= plain for %d/%d, "
           "max epi dist %.2e px\n",
           n, sigma, 1e3 * mean_opt, 1e3 * mean_plain,
           mean_opt / mean_plain, n_improved, n, maxep);
    sprintf(line, "noisy: ensemble size %d >= 200", n);
    CHECK(n >= 200, line);
    CHECK(maxep < 1e-9, "noisy: corrected pairs on epipolar constraint (< 1e-9)");
    CHECK(mean_opt <= mean_plain * 1.005,
          "noisy: optimal mean 3-D error <= plain * 1.005");
    CHECK(cost_violations == 0,
          "noisy: HS correction cost <= any consistent competitor (optimality)");
}

static void test_optimal_degenerate(void)
{
    mv_camera c1, c2;
    const mv_camera *cams[2];
    double F[9], Ft[9], e2h[3], uv[4], Xf[3], Xp[3];
    int rc_c, rc_o, rc_p;
    make_pair(&c1, &c2);
    cams[0] = &c1;
    cams[1] = &c2;
    mv_fundamental_from_cams(F, &c1, &c2);
    /* place the image-2 measurement exactly at the image-2 epipole */
    mv_mat_transpose(Ft, F, 3, 3);
    CHECK(mv_nullvec(e2h, Ft, 3, 3) == MV_OK, "degenerate: epipole found");
    CHECK(fabs(e2h[2]) > 1e-12, "degenerate: epipole is finite for this rig");
    uv[0] = 320.0;
    uv[1] = 240.0;
    uv[2] = e2h[0] / e2h[2];
    uv[3] = e2h[1] / e2h[2];
    rc_c = mv_epipolar_correct(Xf, Xf + 1, F, uv, uv + 2);
    CHECK(rc_c == MV_ERR, "degenerate: epipole-at-point rejected by correction");
    rc_o = mv_triangulate_optimal(Xf, &c1, &c2, uv, uv + 2);
    rc_p = mv_triangulate(Xp, cams, uv, 2);
    CHECK(rc_o == rc_p &&
          (rc_o != MV_OK || dist3(Xf, Xp) < 1e-9),
          "degenerate: falls back to plain triangulation");
}

/* --- (B) weighted circle fit --------------------------------------------- */

/* unweighted known-radius Gauss-Newton fit (mirror of the estimator in
 * demo/diagnose.c / demo/track_robot.c), for independent comparison */
static void circle_fit_unw(double c[2], const double *pts, int n, double R)
{
    int i, it;
    c[0] = c[1] = 0.0;
    for (i = 0; i < n; i++) {
        c[0] += pts[2 * i];
        c[1] += pts[2 * i + 1];
    }
    c[0] /= n;
    c[1] /= n;
    for (it = 0; it < 50; it++) {
        double a = 0.0, b = 0.0, d2 = 0.0, gx = 0.0, gz = 0.0, det, sx, sz;
        for (i = 0; i < n; i++) {
            double dx = pts[2 * i] - c[0], dz = pts[2 * i + 1] - c[1];
            double d = sqrt(dx * dx + dz * dz), rr, jx, jz;
            if (d < 1e-9)
                continue;
            rr = d - R;
            jx = -dx / d;
            jz = -dz / d;
            a += jx * jx;
            b += jx * jz;
            d2 += jz * jz;
            gx += jx * rr;
            gz += jz * rr;
        }
        det = a * d2 - b * b;
        if (fabs(det) < 1e-12)
            break;
        sx = -(d2 * gx - b * gz) / det;
        sz = -(a * gz - b * gx) / det;
        c[0] += sx;
        c[1] += sz;
        if (fabs(sx) + fabs(sz) < 1e-14)
            break;
    }
}

/* 20 points on the 180-degree near half of the circle as seen from the
 * camera at the origin (the visible silhouette arc of a disk robot);
 * the dominant noise axis is the per-point viewing (depth) direction */
static void make_arc(double *pts, double *dirs, double *sig, int n,
                     const double ctr[2], double R,
                     double s_along, double s_across, int noisy)
{
    double base = atan2(-ctr[1], -ctr[0]);
    int i;
    for (i = 0; i < n; i++) {
        double phi = base + ((double)i / (n - 1) - 0.5)
                          * (180.0 * MV_PI / 180.0);
        double px = ctr[0] + R * cos(phi);
        double pz = ctr[1] + R * sin(phi);
        double dn = sqrt(px * px + pz * pz);
        double ax = px / dn, az = pz / dn;
        double g1 = noisy ? grand_() : 0.0;
        double g2 = noisy ? grand_() : 0.0;
        pts[2 * i] = px + s_along * g1 * ax + s_across * g2 * (-az);
        pts[2 * i + 1] = pz + s_along * g1 * az + s_across * g2 * ax;
        dirs[2 * i] = ax;
        dirs[2 * i + 1] = az;
        sig[2 * i] = s_along;
        sig[2 * i + 1] = s_across;
    }
}

static void test_circle_noisy(void)
{
    /* radius and anisotropy match the documented sigma=0.6 px regime:
     * ~40 mm depth-elongated triangulation noise, 6:1 along the view
     * ray (doc: depth RMS 40 mm vs ~15 mm under the isotropy
     * assumption) */
    const double R = 0.17;
    const double ctrue[2] = { 0.30, 4.20 };
    const double s_along = 0.040, s_across = 0.040 / 6.0;
    double sum_u = 0.0, sum_w = 0.0, sum_f = 0.0, factor;
    int seed, ok = 0;
    for (seed = 0; seed < 10; seed++) {
        double pts[40], dirs[40], sig[40], cu[2], cw[2], eu, ew;
        rs = 977ULL * (unsigned long long)(seed + 1);
        make_arc(pts, dirs, sig, 20, ctrue, R, s_along, s_across, 1);
        circle_fit_unw(cu, pts, 20, R);
        if (mv_circle_fit_weighted(cw, pts, dirs, sig, 20, R) == MV_OK)
            ok++;
        eu = sqrt((cu[0] - ctrue[0]) * (cu[0] - ctrue[0])
                  + (cu[1] - ctrue[1]) * (cu[1] - ctrue[1]));
        ew = sqrt((cw[0] - ctrue[0]) * (cw[0] - ctrue[0])
                  + (cw[1] - ctrue[1]) * (cw[1] - ctrue[1]));
        sum_u += eu;
        sum_w += ew;
        sum_f += eu / ew;
    }
    factor = sum_f / 10.0;
    printf("  circle fit, 10 seeds: mean center err unweighted %.1f mm, "
           "weighted %.1f mm (ratio of means %.2fx); "
           "per-seed factor avg %.2fx\n",
           1e3 * sum_u / 10.0, 1e3 * sum_w / 10.0, sum_u / sum_w, factor);
    CHECK(ok == 10, "circle: weighted fit succeeds on all seeds");
    CHECK(sum_w < sum_u, "circle: weighted center error < unweighted");
    CHECK(factor >= 1.3, "circle: improvement factor >= 1.3 over 10 seeds");
}

static void test_circle_exact(void)
{
    const double R = 0.17;
    const double ctrue[2] = { 0.30, 4.20 };
    double pts[40], dirs[40], sig[40], cw[2], e;
    make_arc(pts, dirs, sig, 20, ctrue, R, 0.040, 0.040 / 6.0, 0);
    CHECK(mv_circle_fit_weighted(cw, pts, dirs, sig, 20, R) == MV_OK,
          "circle sigma=0: weighted fit succeeds");
    e = sqrt((cw[0] - ctrue[0]) * (cw[0] - ctrue[0])
             + (cw[1] - ctrue[1]) * (cw[1] - ctrue[1]));
    CHECK(e < 1e-9, "circle sigma=0: weighted fit exact");
}

int main(void)
{
    test_optimal_exact();
    test_optimal_noisy();
    test_optimal_degenerate();
    test_circle_noisy();
    test_circle_exact();
    if (failures) {
        printf("\n%d FAILURE(S)\n", failures);
        return 1;
    }
    printf("\nall optimal tests passed\n");
    return 0;
}
