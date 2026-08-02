/* Simulated tracking of a disk-like vacuum robot: the numbers this prints
 * are quoted in doc/multiview.tex (section "Basic results").
 *
 * A disk robot (radius 0.17 m, height 0.09 m) wanders the monitored volume
 * at constant speed with random turns, observed by the standard two-camera
 * rig. Feature points sit on its cylindrical bumper band; the far side of
 * the band is self-occluded per camera (normal test), so each frame sees
 * only a partial arc — the realistic hard case. Per frame: project visible
 * points with pixel noise, gate correspondences with the analytic F,
 * triangulate, then estimate the robot center two ways:
 *   (a) naive centroid of the triangulated cloud (partial-view biased),
 *   (b) known-radius circle fit in the ground plane (Gauss-Newton).
 * Both are scored against the known simulated trajectory. Deterministic
 * (fixed-seed LCG). Writes out_track.ply (ground truth green, estimate
 * red) for visual inspection. */

#include <math.h>
#include <stdio.h>

#include "mv/mv.h"

#define MV_PI 3.14159265358979323846

#define FPS 30.0
#define NFRAMES 600            /* 20 seconds */
#define SPEED 0.3              /* m/s */
#define NOISE_SIGMA 0.3        /* px */
#define ROBOT_R 0.17           /* m */
#define BAND_TOP 0.71          /* y of upper band ring (y is down, floor 0.8) */
#define BAND_BOT 0.78          /* y of lower band ring */
#define NRING 18               /* points per band ring */
#define IMG_W 640
#define IMG_H 480
#define XMIN (-1.0)
#define XMAX 1.0
#define ZMIN 3.7
#define ZMAX 5.8
#define GATE_PX 1.5            /* epipolar gate */

static unsigned long long rng_state = 12345;

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

/* rim point visible to a camera if its outward normal faces the camera */
static int rim_visible(const double p[3], const double n[3],
                       const mv_camera *cam)
{
    double C[3], d[3];
    int i;
    mv_cam_center(C, cam);
    for (i = 0; i < 3; i++)
        d[i] = C[i] - p[i];
    return mv_dot(n, d, 3) > 0.0;
}

/* Gauss-Newton fit of a circle with KNOWN radius R to ground-plane points:
 * minimize sum (|p_i - c| - R)^2 over the 2D center c = (x, z) */
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
        double a = 0.0, b = 0.0, d2 = 0.0, gx = 0.0, gz = 0.0, det;
        double dx0, dz0;
        for (i = 0; i < n; i++) {
            double dx = xz[2 * i] - c[0], dz = xz[2 * i + 1] - c[1];
            double d = sqrt(dx * dx + dz * dz);
            double r, jx, jz;
            if (d < 1e-9)
                continue;
            r = d - R;
            jx = -dx / d;
            jz = -dz / d;
            a += jx * jx;
            b += jx * jz;
            d2 += jz * jz;
            gx += jx * r;
            gz += jz * r;
        }
        det = a * d2 - b * b;
        if (fabs(det) < 1e-12)
            break;
        dx0 = -(d2 * gx - b * gz) / det;
        dz0 = -(a * gz - b * gx) / det;
        c[0] += dx0;
        c[1] += dz0;
        if (fabs(dx0) + fabs(dz0) < 1e-10)
            break;
    }
}

int main(void)
{
    mv_camera c1, c2;
    const mv_camera *cams[2];
    double C2pos[3] = { 0.5, 0.0, 0.0 };
    double F[9];
    double gt[NFRAMES][2], est[NFRAMES][2], cent[NFRAMES][2];
    int nvis[NFRAMES];
    double cx = 0.0, cz = 4.7, dirx, dirz, turn_at = 0.0, tnow = 0.0;
    double se_x = 0.0, se_z = 0.0, max_e = 0.0, se_cent = 0.0;
    double se_speed = 0.0;
    double sum_vis = 0.0, max_epi = 0.0;
    int gated_out = 0, nspeed = 0;
    mv_cloud cloud;
    int f, i, k;

    mv_cam_set_K(&c1, 800.0, 800.0, 320.0, 240.0);
    mv_cam_set_identity_pose(&c1);
    for (i = 0; i < 5; i++)
        c1.k[i] = 0.0;
    c2 = c1;
    mv_cam_set_pose_yaw(&c2, -6.0 * MV_PI / 180.0, C2pos);
    cams[0] = &c1;
    cams[1] = &c2;
    mv_fundamental_from_cams(F, &c1, &c2);

    {
        double psi = 2.0 * MV_PI * urand();
        dirx = cos(psi);
        dirz = sin(psi);
    }

    for (f = 0; f < NFRAMES; f++) {
        double xz[2 * 2 * NRING];
        double dt = 1.0 / FPS;
        int nv = 0;

        /* --- simulate motion: constant speed, random turns, wall bounce */
        tnow = f * dt;
        if (tnow >= turn_at) {
            double psi = 2.0 * MV_PI * urand();
            dirx = cos(psi);
            dirz = sin(psi);
            turn_at = tnow + 2.0 + 2.0 * urand();
        }
        cx += SPEED * dt * dirx;
        cz += SPEED * dt * dirz;
        if (cx < XMIN + ROBOT_R || cx > XMAX - ROBOT_R) {
            dirx = -dirx;
            cx += 2.0 * SPEED * dt * dirx;
        }
        if (cz < ZMIN + ROBOT_R || cz > ZMAX - ROBOT_R) {
            dirz = -dirz;
            cz += 2.0 * SPEED * dt * dirz;
        }
        gt[f][0] = cx;
        gt[f][1] = cz;

        /* --- observe: bumper-band points, self-occlusion, noise, gating */
        for (k = 0; k < 2 * NRING; k++) {
            double th = 2.0 * MV_PI * (k % NRING) / NRING;
            double y = (k < NRING) ? BAND_TOP : BAND_BOT;
            double P[3], n[3], uv[4], p1[2], p2[2], X[3], d;
            P[0] = cx + ROBOT_R * cos(th);
            P[1] = y;
            P[2] = cz + ROBOT_R * sin(th);
            n[0] = cos(th);
            n[1] = 0.0;
            n[2] = sin(th);
            if (!rim_visible(P, n, &c1) || !rim_visible(P, n, &c2))
                continue;
            if (mv_cam_project(p1, &c1, P) != MV_OK || !in_bounds(p1))
                continue;
            if (mv_cam_project(p2, &c2, P) != MV_OK || !in_bounds(p2))
                continue;
            p1[0] += NOISE_SIGMA * grand();
            p1[1] += NOISE_SIGMA * grand();
            p2[0] += NOISE_SIGMA * grand();
            p2[1] += NOISE_SIGMA * grand();
            d = mv_sym_epipolar_dist(F, p1, p2);
            if (d > max_epi)
                max_epi = d;
            if (d > GATE_PX) {
                gated_out++;
                continue;
            }
            uv[0] = p1[0]; uv[1] = p1[1];
            uv[2] = p2[0]; uv[3] = p2[1];
            if (mv_triangulate(X, cams, uv, 2) != MV_OK)
                continue;
            xz[2 * nv] = X[0];
            xz[2 * nv + 1] = X[2];
            nv++;
        }
        nvis[f] = nv;
        sum_vis += nv;

        /* --- estimate center: naive centroid vs known-radius circle fit */
        {
            double cc[2] = { 0.0, 0.0 };
            for (i = 0; i < nv; i++) {
                cc[0] += xz[2 * i];
                cc[1] += xz[2 * i + 1];
            }
            cent[f][0] = cc[0] / nv;
            cent[f][1] = cc[1] / nv;
        }
        fit_circle_known_r(est[f], xz, nv, ROBOT_R);
    }

    /* --- score against ground truth ------------------------------------ */
    for (f = 0; f < NFRAMES; f++) {
        double ex = est[f][0] - gt[f][0];
        double ez = est[f][1] - gt[f][1];
        double e = sqrt(ex * ex + ez * ez);
        double ecx = cent[f][0] - gt[f][0];
        double ecz = cent[f][1] - gt[f][1];
        se_x += ex * ex;
        se_z += ez * ez;
        se_cent += ecx * ecx + ecz * ecz;
        if (e > max_e)
            max_e = e;
    }
    /* speed by central difference over +-2 frames, vs same on ground truth */
    for (f = 2; f < NFRAMES - 2; f++) {
        double dt4 = 4.0 / FPS;
        double vex = (est[f + 2][0] - est[f - 2][0]) / dt4;
        double vez = (est[f + 2][1] - est[f - 2][1]) / dt4;
        double vgx = (gt[f + 2][0] - gt[f - 2][0]) / dt4;
        double vgz = (gt[f + 2][1] - gt[f - 2][1]) / dt4;
        double dv = sqrt((vex - vgx) * (vex - vgx) + (vez - vgz) * (vez - vgz));
        se_speed += dv * dv;
        nspeed++;
    }

    printf("multiview robot-tracking experiment (seed 12345)\n");
    printf("------------------------------------------------\n");
    printf("robot               : disk r=%.2f m, band features, "
           "%.1f m/s, random turns\n", ROBOT_R, SPEED);
    printf("frames              : %d at %.0f fps (%.0f s)\n",
           NFRAMES, FPS, NFRAMES / FPS);
    printf("visible band points : %.1f mean per frame (of %d)\n",
           sum_vis / NFRAMES, 2 * NRING);
    printf("epipolar gate       : max sym dist %.2f px, %d rejected\n",
           max_epi, gated_out);
    printf("\ncenter estimate, known-radius circle fit\n");
    printf("  RMS lateral (x)   : %.1f mm\n",
           1000.0 * sqrt(se_x / NFRAMES));
    printf("  RMS depth (z)     : %.1f mm\n",
           1000.0 * sqrt(se_z / NFRAMES));
    printf("  RMS horizontal    : %.1f mm\n",
           1000.0 * sqrt((se_x + se_z) / NFRAMES));
    printf("  max error         : %.1f mm\n", 1000.0 * max_e);
    printf("\ncenter estimate, naive centroid (for contrast)\n");
    printf("  RMS horizontal    : %.1f mm\n",
           1000.0 * sqrt(se_cent / NFRAMES));
    printf("\nspeed estimate (central difference, +-2 frames)\n");
    printf("  RMS speed error   : %.1f mm/s (true speed %.0f mm/s)\n",
           1000.0 * sqrt(se_speed / nspeed), 1000.0 * SPEED);

    /* --- trajectory cloud: ground truth green, estimate red ------------ */
    mv_cloud_init(&cloud, 1);
    for (f = 0; f < NFRAMES; f++) {
        unsigned char green[3] = { 40, 200, 40 };
        unsigned char red[3] = { 220, 50, 50 };
        double p[3];
        p[0] = gt[f][0]; p[1] = 0.745; p[2] = gt[f][1];
        mv_cloud_push(&cloud, p, green);
        p[0] = est[f][0]; p[1] = 0.745; p[2] = est[f][1];
        mv_cloud_push(&cloud, p, red);
    }
    mv_cloud_write_ply("out_track.ply", &cloud);
    printf("\nwrote out_track.ply (%d points: truth green, estimate red)\n",
           cloud.n);
    mv_cloud_free(&cloud);
    return 0;
}
