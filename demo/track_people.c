/* Human occupancy: who is present, when, and along which trajectories.
 * Numbers quoted in doc/multiview.tex ("Scenario studies").
 *
 * A small known cohort (three residents of distinct heights) makes
 * scheduled visits to the monitored room: enter at the door, walk to a
 * target (desk or window), linger, walk back, exit. Each person is a
 * cylinder with torso-band features (self-occluded far side) plus a head
 * point. Per frame: triangulate, known-radius circle fit for the ground
 * position, head height for identity. Scored: episode detection, identity
 * accuracy (nearest known height), main-trajectory classification,
 * position RMS, entry/exit timing. Deterministic. */

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mv/mv.h"

#define MV_PI 3.14159265358979323846
#define FPS 10.0
#define BODY_R 0.18
#define FLOOR_Y 0.8
#define NEP 13

static const double HEIGHTS[3] = { 1.62, 1.75, 1.88 };
static const char *NAMES[3] = { "A", "B", "C" };
static const double DOOR[2] = { 0.0, 5.7 };
static const double TARGET[2][2] = { { -0.7, 4.1 }, { 0.75, 4.4 } };
static const char *TNAME[2] = { "desk", "window" };

/* schedule: person, nominal hour, target */
static const struct { int person; double hour; int target; } SCHED[NEP] = {
    { 0, 7.1, 0 }, { 0, 8.2, 1 }, { 0, 19.3, 0 }, { 0, 21.5, 0 },
    { 1, 9.7, 1 }, { 1, 11.4, 1 }, { 1, 13.2, 0 }, { 1, 15.8, 1 },
    { 1, 16.9, 1 },
    { 2, 18.4, 0 }, { 2, 20.1, 1 }, { 2, 22.3, 0 }, { 2, 23.2, 0 }
};

static unsigned long long rng_state = 31337;

static double urand(void)
{
    rng_state = rng_state * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)((rng_state >> 11) & ((1ULL << 53) - 1)) / 9007199254740992.0;
}

static double grand_(void)
{
    double u1 = urand(), u2 = urand();
    if (u1 < 1e-300)
        u1 = 1e-300;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * MV_PI * u2);
}

static int in_img(const double uv[2])
{
    return uv[0] >= 0.0 && uv[0] < 640.0 && uv[1] >= 0.0 && uv[1] < 480.0;
}

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
        double a = 0, b = 0, d2 = 0, gx = 0, gz = 0, det, sx, sz;
        for (i = 0; i < n; i++) {
            double dx = xz[2 * i] - c[0], dz = xz[2 * i + 1] - c[1];
            double d = sqrt(dx * dx + dz * dz), rr, jx, jz;
            if (d < 1e-9)
                continue;
            rr = d - R;
            jx = -dx / d;
            jz = -dz / d;
            a += jx * jx; b += jx * jz; d2 += jz * jz;
            gx += jx * rr; gz += jz * rr;
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

int main(void)
{
    mv_camera c1, c2;
    const mv_camera *cams[2];
    double C2pos[3] = { 0.5, 0.0, 0.0 };
    int detected = 0, id_ok = 0, path_ok = 0;
    double se_h = 0.0, se_pos = 0.0, lat_in = 0.0, lat_out = 0.0;
    long n_pos = 0;
    double minutes_gt[3] = { 0, 0, 0 }, minutes_est[3] = { 0, 0, 0 };
    int e, i, k;

    mv_cam_set_K(&c1, 800.0, 800.0, 320.0, 240.0);
    mv_cam_set_identity_pose(&c1);
    memset(c1.k, 0, sizeof(c1.k));
    c2 = c1;
    mv_cam_set_pose_yaw(&c2, -6.0 * MV_PI / 180.0, C2pos);
    cams[0] = &c1;
    cams[1] = &c2;

    printf("multiview occupancy experiment (seed 31337)\n");
    printf("-------------------------------------------\n");
    printf("cohort heights      : %.2f %.2f %.2f m\n",
           HEIGHTS[0], HEIGHTS[1], HEIGHTS[2]);
    printf("episodes scheduled  : %d over one simulated day\n\n", NEP);

    for (e = 0; e < NEP; e++) {
        int person = SCHED[e].person, target = SCHED[e].target;
        double h = HEIGHTS[person];
        double tx = TARGET[target][0], tz = TARGET[target][1];
        double len = sqrt((tx - DOOR[0]) * (tx - DOOR[0])
                          + (tz - DOOR[1]) * (tz - DOOR[1]));
        double linger = 5.0 + 10.0 * urand();
        double t_total = 2.0 * len / 1.0 + linger;
        int nf = (int)(t_total * FPS);
        double dt = 1.0 / FPS;
        double sum_hgt = 0.0, far_d = -1.0, far_p[2] = { 0, 0 };
        int nh = 0, first_conf = -1, last_seen = -1, run = 0;
        int id;
        double ep_min;

        for (i = 0; i < nf; i++) {
            double t = i * dt, u, cx, cz;
            double xz[2 * 24], hx[2], hp1[2], hp2[2], huv[4], H3[3];
            int nv = 0, j;
            /* position along door->target->door with linger */
            if (t < len)
                u = t / len;
            else if (t < len + linger)
                u = 1.0;
            else
                u = 1.0 - (t - len - linger) / len;
            cx = DOOR[0] + u * (tx - DOOR[0]) + 0.02 * grand_();
            cz = DOOR[1] + u * (tz - DOOR[1]) + 0.02 * grand_();

            /* torso band */
            for (k = 0; k < 24; k++) {
                double th = 2.0 * MV_PI * (k % 12) / 12.0;
                double y = (k < 12) ? FLOOR_Y - 1.2 : FLOOR_Y - 1.05;
                double P[3], n[3], C[3], d[3], p1[2], p2[2], uv[4], X[3];
                P[0] = cx + BODY_R * cos(th);
                P[1] = y;
                P[2] = cz + BODY_R * sin(th);
                n[0] = cos(th); n[1] = 0.0; n[2] = sin(th);
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
                uv[0] = p1[0] + 0.3 * grand_();
                uv[1] = p1[1] + 0.3 * grand_();
                uv[2] = p2[0] + 0.3 * grand_();
                uv[3] = p2[1] + 0.3 * grand_();
                if (mv_triangulate(X, cams, uv, 2) != MV_OK)
                    continue;
                xz[2 * nv] = X[0];
                xz[2 * nv + 1] = X[2];
                nv++;
            }
            if (nv < 4) {
                run = 0;
                continue;
            }
            run++;
            last_seen = i;
            if (run >= 3 && first_conf < 0)
                first_conf = i;

            fit_circle_known_r(hx, xz, nv, BODY_R);
            {
                double d0 = hx[0] - DOOR[0], d1 = hx[1] - DOOR[1];
                double dd = sqrt(d0 * d0 + d1 * d1);
                if (dd > far_d) {
                    far_d = dd;
                    far_p[0] = hx[0];
                    far_p[1] = hx[1];
                }
            }
            se_pos += (hx[0] - cx) * (hx[0] - cx)
                    + (hx[1] - cz) * (hx[1] - cz);
            n_pos++;

            /* head point -> height */
            {
                double Ph[3];
                Ph[0] = cx; Ph[1] = FLOOR_Y - h; Ph[2] = cz;
                if (mv_cam_project(hp1, &c1, Ph) == MV_OK && in_img(hp1)
                    && mv_cam_project(hp2, &c2, Ph) == MV_OK && in_img(hp2)) {
                    huv[0] = hp1[0] + 0.3 * grand_();
                    huv[1] = hp1[1] + 0.3 * grand_();
                    huv[2] = hp2[0] + 0.3 * grand_();
                    huv[3] = hp2[1] + 0.3 * grand_();
                    if (mv_triangulate(H3, cams, huv, 2) == MV_OK) {
                        sum_hgt += FLOOR_Y - H3[1];
                        nh++;
                    }
                }
            }
        }

        if (first_conf < 0 || nh == 0)
            continue;
        detected++;
        {
            double hest = sum_hgt / nh, bd = 1e9;
            id = 0;
            for (k = 0; k < 3; k++)
                if (fabs(hest - HEIGHTS[k]) < bd) {
                    bd = fabs(hest - HEIGHTS[k]);
                    id = k;
                }
            se_h += (hest - h) * (hest - h);
            if (id == person)
                id_ok++;
        }
        {
            double d0 = far_p[0] - TARGET[0][0], d1 = far_p[1] - TARGET[0][1];
            double e0 = sqrt(d0 * d0 + d1 * d1);
            d0 = far_p[0] - TARGET[1][0];
            d1 = far_p[1] - TARGET[1][1];
            if ((e0 < sqrt(d0 * d0 + d1 * d1) ? 0 : 1) == target)
                path_ok++;
        }
        lat_in += first_conf * dt;
        lat_out += (nf - 1 - last_seen) * dt;
        ep_min = t_total / 60.0;
        minutes_gt[person] += ep_min;
        minutes_est[id] += (last_seen - first_conf) * dt / 60.0;
        printf("  ep %2d  %02.0f:%02.0f  gt %s->%s  est %s->%s  %s\n",
               e, floor(SCHED[e].hour),
               60.0 * (SCHED[e].hour - floor(SCHED[e].hour)),
               NAMES[person], TNAME[target], NAMES[id],
               TNAME[(far_d > 0
                      && (fabs(far_p[0] - TARGET[0][0])
                          + fabs(far_p[1] - TARGET[0][1]))
                       < (fabs(far_p[0] - TARGET[1][0])
                          + fabs(far_p[1] - TARGET[1][1]))) ? 0 : 1],
               (id == person) ? "ok" : "ID-ERR");
    }

    printf("\nepisodes detected   : %d/%d\n", detected, NEP);
    printf("identity accuracy   : %d/%d\n", id_ok, detected);
    printf("path accuracy       : %d/%d\n", path_ok, detected);
    printf("height RMS error    : %.1f mm\n",
           1000.0 * sqrt(se_h / detected));
    printf("ground position RMS : %.1f mm\n",
           1000.0 * sqrt(se_pos / n_pos));
    printf("mean entry latency  : %.2f s   mean exit latency: %.2f s\n",
           lat_in / detected, lat_out / detected);
    printf("\noccupancy minutes (gt/est): A %.1f/%.1f  B %.1f/%.1f  "
           "C %.1f/%.1f\n", minutes_gt[0], minutes_est[0],
           minutes_gt[1], minutes_est[1], minutes_gt[2], minutes_est[2]);
    return 0;
}
