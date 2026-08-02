/* Error-budget diagnostics: the numbers this prints are quoted in
 * doc/multiview.tex ("Error budget and consistency").
 *
 * Three checks that separate noise-driven error from bias (a bias that
 * survives sigma -> 0, or error that fails to scale linearly in sigma,
 * indicates a defect):
 *   A. triangulation noise sweep, measured vs per-point analytic
 *      prediction sqrt(E[Z^4]) * sqrt(2) sigma / (f B);
 *   B. robot-tracking noise sweep, with the circle-fit "arc factor"
 *      (measured error over naive sqrt(n) average) tracked across sigma;
 *   C. calibration seed sweep: the distribution of focal-length errors
 *      across noise realizations, to judge whether any single run is
 *      typical or anomalous.
 * Deterministic given the seeds below. */

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mv/mv.h"

#define MV_PI 3.14159265358979323846

struct rng { unsigned long long s; };

static double urand(struct rng *r)
{
    r->s = r->s * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)((r->s >> 11) & ((1ULL << 53) - 1)) / 9007199254740992.0;
}

static double grand(struct rng *r)
{
    double u1 = urand(r), u2 = urand(r);
    if (u1 < 1e-300)
        u1 = 1e-300;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * MV_PI * u2);
}

static void make_rig(mv_camera *c1, mv_camera *c2)
{
    double C2pos[3] = { 0.5, 0.0, 0.0 };
    mv_cam_set_K(c1, 800.0, 800.0, 320.0, 240.0);
    mv_cam_set_identity_pose(c1);
    memset(c1->k, 0, sizeof(c1->k));
    *c2 = *c1;
    mv_cam_set_pose_yaw(c2, -6.0 * MV_PI / 180.0, C2pos);
}

static int in_img(const double uv[2])
{
    return uv[0] >= 0.0 && uv[0] < 640.0 && uv[1] >= 0.0 && uv[1] < 480.0;
}

/* --- A: triangulation sweep --------------------------------------------- */

static void tri_stats(double sigma, unsigned long long seed,
                      double *rms_z, double *rms_lat, double *pred_z)
{
    mv_camera c1, c2;
    const mv_camera *cams[2];
    struct rng r = { seed };
    double sez = 0.0, sel = 0.0, sp = 0.0;
    int i, n = 0;

    make_rig(&c1, &c2);
    cams[0] = &c1;
    cams[1] = &c2;
    for (i = 0; i < 250; i++) {
        double X[3], p1[2], p2[2], uv[4], Xr[3];
        X[0] = -1.0 + 2.0 * urand(&r);
        X[1] = -0.75 + 1.5 * urand(&r);
        X[2] = 3.5 + 2.5 * urand(&r);
        if (mv_cam_project(p1, &c1, X) != MV_OK || !in_img(p1))
            continue;
        if (mv_cam_project(p2, &c2, X) != MV_OK || !in_img(p2))
            continue;
        uv[0] = p1[0] + sigma * grand(&r);
        uv[1] = p1[1] + sigma * grand(&r);
        uv[2] = p2[0] + sigma * grand(&r);
        uv[3] = p2[1] + sigma * grand(&r);
        if (mv_triangulate(Xr, cams, uv, 2) != MV_OK)
            continue;
        sez += (Xr[2] - X[2]) * (Xr[2] - X[2]);
        sel += (Xr[0] - X[0]) * (Xr[0] - X[0])
             + (Xr[1] - X[1]) * (Xr[1] - X[1]);
        /* per-point depth-law prediction, sigma_d = sqrt(2) sigma */
        sp += pow(X[2] * X[2] * sqrt(2.0) * sigma / (800.0 * 0.5), 2.0);
        n++;
    }
    *rms_z = sqrt(sez / n);
    *rms_lat = sqrt(sel / n);
    *pred_z = sqrt(sp / n);
}

/* --- B: tracking sweep --------------------------------------------------- */

static void fit_circle_known_r(double c[2], const double *xz, int n, double R)
{
    int i, it;
    c[0] = c[1] = 0.0;
    for (i = 0; i < n; i++) {
        c[0] += xz[2 * i];
        c[1] += xz[2 * i + 1];
    }
    c[0] /= n;
    c[1] /= n;
    for (it = 0; it < 20; it++) {
        double a = 0.0, b = 0.0, d2 = 0.0, gx = 0.0, gz = 0.0, det, sx, sz;
        for (i = 0; i < n; i++) {
            double dx = xz[2 * i] - c[0], dz = xz[2 * i + 1] - c[1];
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
        if (fabs(sx) + fabs(sz) < 1e-12)
            break;
    }
}

static void track_stats(double sigma, unsigned long long seed,
                        double *rms_x, double *rms_z, double *arc_factor)
{
    mv_camera c1, c2;
    const mv_camera *cams[2];
    struct rng r = { seed };
    double cx = 0.0, cz = 4.7, dirx, dirz, turn_at = 0.0;
    double sex = 0.0, sez = 0.0, snaive = 0.0;
    int f, k, nf = 0;

    make_rig(&c1, &c2);
    cams[0] = &c1;
    cams[1] = &c2;
    {
        double psi = 2.0 * MV_PI * urand(&r);
        dirx = cos(psi);
        dirz = sin(psi);
    }
    for (f = 0; f < 600; f++) {
        double xz[2 * 36], est[2];
        double dt = 1.0 / 30.0, tnow = f * dt;
        double szc = 0.0;
        int nv = 0;
        if (tnow >= turn_at) {
            double psi = 2.0 * MV_PI * urand(&r);
            dirx = cos(psi);
            dirz = sin(psi);
            turn_at = tnow + 2.0 + 2.0 * urand(&r);
        }
        cx += 0.3 * dt * dirx;
        cz += 0.3 * dt * dirz;
        if (cx < -1.0 + 0.17 || cx > 1.0 - 0.17) {
            dirx = -dirx;
            cx += 2.0 * 0.3 * dt * dirx;
        }
        if (cz < 3.7 + 0.17 || cz > 5.8 - 0.17) {
            dirz = -dirz;
            cz += 2.0 * 0.3 * dt * dirz;
        }
        for (k = 0; k < 36; k++) {
            double th = 2.0 * MV_PI * (k % 18) / 18.0;
            double y = (k < 18) ? 0.71 : 0.78;
            double P[3], n[3], C[3], d[3], p1[2], p2[2], uv[4], X[3];
            int j;
            P[0] = cx + 0.17 * cos(th);
            P[1] = y;
            P[2] = cz + 0.17 * sin(th);
            n[0] = cos(th);
            n[1] = 0.0;
            n[2] = sin(th);
            mv_cam_center(C, &c1);
            for (j = 0; j < 3; j++)
                d[j] = C[j] - P[j];
            if (mv_dot(n, d, 3) <= 0.0)
                continue;
            mv_cam_center(C, &c2);
            for (j = 0; j < 3; j++)
                d[j] = C[j] - P[j];
            if (mv_dot(n, d, 3) <= 0.0)
                continue;
            if (mv_cam_project(p1, &c1, P) != MV_OK || !in_img(p1))
                continue;
            if (mv_cam_project(p2, &c2, P) != MV_OK || !in_img(p2))
                continue;
            uv[0] = p1[0] + sigma * grand(&r);
            uv[1] = p1[1] + sigma * grand(&r);
            uv[2] = p2[0] + sigma * grand(&r);
            uv[3] = p2[1] + sigma * grand(&r);
            if (mv_triangulate(X, cams, uv, 2) != MV_OK)
                continue;
            xz[2 * nv] = X[0];
            xz[2 * nv + 1] = X[2];
            szc += pow(P[2] * P[2] * sqrt(2.0) * sigma / 400.0, 2.0);
            nv++;
        }
        if (nv < 4)
            continue;
        fit_circle_known_r(est, xz, nv, 0.17);
        sex += (est[0] - cx) * (est[0] - cx);
        sez += (est[1] - cz) * (est[1] - cz);
        /* naive sqrt(n) average of the per-point depth prediction */
        snaive += szc / (double)(nv * nv);
        nf++;
    }
    *rms_x = sqrt(sex / nf);
    *rms_z = sqrt(sez / nf);
    *arc_factor = (sigma > 0.0) ? sqrt(sez / nf) / sqrt(snaive / nf) : 0.0;
}

/* --- C: calibration seed sweep ------------------------------------------ */

static double calib_fx_err(double sigma, unsigned long long seed,
                           double *reproj)
{
    const double pose_deg[6][2] = {
        { 0.0, 0.0 }, { 30.0, 0.0 }, { -30.0, 15.0 },
        { 20.0, -35.0 }, { -20.0, -25.0 }, { 35.0, 25.0 }
    };
    struct rng r = { seed };
    double obj[2 * MV_TARGET_LETTER_COLS * MV_TARGET_LETTER_ROWS];
    double img[6][2 * MV_TARGET_LETTER_COLS * MV_TARGET_LETTER_ROWS];
    mv_calib_view views[6];
    mv_camera gt[6], est[6];
    double K[9], k_radial[2], center[2];
    int n, v, i, j;

    n = mv_target_checkerboard(obj, MV_TARGET_LETTER_COLS,
                               MV_TARGET_LETTER_ROWS,
                               MV_TARGET_LETTER_SQUARE);
    center[0] = 0.5 * (MV_TARGET_LETTER_COLS - 1) * MV_TARGET_LETTER_SQUARE;
    center[1] = 0.5 * (MV_TARGET_LETTER_ROWS - 1) * MV_TARGET_LETTER_SQUARE;
    for (v = 0; v < 6; v++) {
        double cxa = cos(pose_deg[v][0] * MV_PI / 180.0);
        double sxa = sin(pose_deg[v][0] * MV_PI / 180.0);
        double cya = cos(pose_deg[v][1] * MV_PI / 180.0);
        double sya = sin(pose_deg[v][1] * MV_PI / 180.0);
        double Rx[9], Ry[9];
        Rx[0] = 1; Rx[1] = 0;   Rx[2] = 0;
        Rx[3] = 0; Rx[4] = cxa; Rx[5] = -sxa;
        Rx[6] = 0; Rx[7] = sxa; Rx[8] = cxa;
        Ry[0] = cya; Ry[1] = 0; Ry[2] = -sya;
        Ry[3] = 0;   Ry[4] = 1; Ry[5] = 0;
        Ry[6] = sya; Ry[7] = 0; Ry[8] = cya;
        mv_cam_set_K(&gt[v], 800.0, 805.0, 322.0, 238.0);
        mv_mat_mul(gt[v].R, Rx, Ry, 3, 3, 3);
        for (j = 0; j < 3; j++)
            gt[v].t[j] = -(gt[v].R[j * 3 + 0] * center[0]
                           + gt[v].R[j * 3 + 1] * center[1]);
        gt[v].t[0] += -0.03 + 0.06 * urand(&r);
        gt[v].t[1] += -0.03 + 0.06 * urand(&r);
        gt[v].t[2] += 0.45 + 0.20 * urand(&r);
        memset(gt[v].k, 0, sizeof(gt[v].k));
        for (i = 0; i < n; i++) {
            double X[3] = { obj[2 * i], obj[2 * i + 1], 0.0 };
            double p[2];
            mv_cam_project(p, &gt[v], X);
            img[v][2 * i] = p[0] + sigma * grand(&r);
            img[v][2 * i + 1] = p[1] + sigma * grand(&r);
        }
        views[v].obj = obj;
        views[v].img = img[v];
        views[v].n = n;
    }
    if (mv_calib_planar(K, est, k_radial, views, 6, 1) != MV_OK)
        return 1e9;
    if (reproj)
        *reproj = mv_calib_reproj_rms(est, views, 6);
    return K[0] - 800.0;
}

int main(void)
{
    const double sigmas[4] = { 0.0, 0.1, 0.3, 0.6 };
    int i;

    printf("multiview error-budget diagnostics\n");
    printf("==================================\n");

    printf("\nA. triangulation: measured vs depth-law prediction\n");
    printf("   sigma   rms_z (mm)   pred_z (mm)   ratio   rms_lat (mm)\n");
    for (i = 0; i < 4; i++) {
        double rz, rl, pz;
        tri_stats(sigmas[i], 42, &rz, &rl, &pz);
        printf("   %.2f  %10.3f  %11.3f  %6s  %11.3f\n",
               sigmas[i], 1000.0 * rz, 1000.0 * pz,
               sigmas[i] > 0.0 ? "" : "-", 1000.0 * rl);
        if (sigmas[i] > 0.0)
            printf("%56.2f\n", rz / pz);
    }

    printf("\nB. tracking: circle-fit error across noise levels\n");
    printf("   sigma   rms_x (mm)   rms_z (mm)   arc factor (z)\n");
    for (i = 0; i < 4; i++) {
        double rx, rz, af;
        track_stats(sigmas[i], 12345, &rx, &rz, &af);
        printf("   %.2f  %10.4f  %10.4f  %10.2f\n",
               sigmas[i], 1000.0 * rx, 1000.0 * rz, af);
    }

    printf("\nC. calibration: fx error across noise realizations "
           "(sigma = 0.3)\n   seeds 1..10: ");
    {
        double sum = 0.0, sum2 = 0.0, rsum = 0.0, rp;
        for (i = 1; i <= 10; i++) {
            double e = calib_fx_err(0.3, (unsigned long long)i, &rp);
            printf("%+.1f ", e);
            sum += e;
            sum2 += e * e;
            rsum += rp;
        }
        printf("\n   mean %+.2f px, sample std %.2f px\n", sum / 10.0,
               sqrt((sum2 - sum * sum / 10.0) / 9.0));
        printf("   (demo seed 7 gave -6.95 px)\n");
        /* residual-consistency (chi-square-style) check: fitted residual
         * RMS vs sqrt(2) sigma sqrt(1 - p/n): p = 42 params (4 K +
         * 2 radial + 6x6 pose), n = 756 residual coordinates */
        printf("   mean reprojection RMS: %.3f px; "
               "noise-model prediction: %.3f px; ratio %.3f\n",
               rsum / 10.0,
               sqrt(2.0) * 0.3 * sqrt(1.0 - 42.0 / 756.0),
               (rsum / 10.0)
               / (sqrt(2.0) * 0.3 * sqrt(1.0 - 42.0 / 756.0)));
    }
    printf("\nC'. calibration linearity: mean |fx err| at sigma = 0.15: ");
    {
        double s = 0.0;
        for (i = 1; i <= 10; i++)
            s += fabs(calib_fx_err(0.15, (unsigned long long)i, NULL));
        printf("%.2f px (expect ~half of sigma=0.3 value)\n", s / 10.0);
    }
    return 0;
}
