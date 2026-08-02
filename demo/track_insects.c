/* Insect counting and entry/exit statistics: numbers quoted in
 * doc/multiview.tex ("Scenario studies").
 *
 * Point-like targets (insects) enter the monitored box through a random
 * face, fly an Ornstein-Uhlenbeck random walk, and leave. Each camera
 * also sees uniform clutter detections. Stereo association pairs the two
 * cameras' detections by the epipolar gate (the analytic F), pairs are
 * triangulated, and a nearest-neighbor tracker with M-of-N confirmation
 * and miss-based termination maintains tracks. Confirmed tracks are
 * counted and their birth/death positions mapped to box faces; both
 * are scored against ground truth. Deterministic. */

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mv/mv.h"

#define MV_PI 3.14159265358979323846
#define FPS 30.0
#define DURATION 300.0          /* seconds */
#define SPAWN_RATE 0.1          /* insects per second */
#define P_DET 0.9               /* joint detection probability per frame */
#define CLUTTER_PER_CAM 2.0     /* mean false detections per camera frame */
#define GATE_PX 1.5
#define ASSOC_GATE 0.12         /* m */
#define CONFIRM_HITS 3
#define KILL_MISSES 5
#define MAXB 32                 /* max simultaneous ground-truth insects */
#define MAXT 128                /* max simultaneous tracks */
#define MAXD 64                 /* max detections per frame */

/* monitored box */
#define BX0 (-1.0)
#define BX1 1.0
#define BY0 (-0.5)
#define BY1 0.79
#define BZ0 3.7
#define BZ1 5.8

static unsigned long long rng_state = 777;

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

/* face ids: 0 x-, 1 x+, 2 y-, 3 y+, 4 z-, 5 z+ */
static const char *face_name[6] = { "x-", "x+", "y-(top)", "y+(bot)",
                                    "z-(near)", "z+(far)" };

static int nearest_face(const double p[3])
{
    double d[6];
    int i, best = 0;
    d[0] = p[0] - BX0; d[1] = BX1 - p[0];
    d[2] = p[1] - BY0; d[3] = BY1 - p[1];
    d[4] = p[2] - BZ0; d[5] = BZ1 - p[2];
    for (i = 1; i < 6; i++)
        if (d[i] < d[best])
            best = i;
    return best;
}

struct insect {
    double p[3], v[3];
    int alive, entry_face;
};

struct track {
    double p[3], born[3], last[3];
    int hits, misses, confirmed, alive;
};

static void spawn_insect(struct insect *b)
{
    int f = (int)(urand() * 6.0) % 6;
    double s = 0.3 + 0.5 * urand();
    double n[3] = { 0, 0, 0 };
    b->p[0] = BX0 + (BX1 - BX0) * urand();
    b->p[1] = BY0 + (BY1 - BY0) * urand();
    b->p[2] = BZ0 + (BZ1 - BZ0) * urand();
    switch (f) {
    case 0: b->p[0] = BX0; n[0] = 1; break;
    case 1: b->p[0] = BX1; n[0] = -1; break;
    case 2: b->p[1] = BY0; n[1] = 1; break;
    case 3: b->p[1] = BY1; n[1] = -1; break;
    case 4: b->p[2] = BZ0; n[2] = 1; break;
    default: b->p[2] = BZ1; n[2] = -1; break;
    }
    b->v[0] = s * n[0] + 0.2 * grand_();
    b->v[1] = s * n[1] + 0.2 * grand_();
    b->v[2] = s * n[2] + 0.2 * grand_();
    b->alive = 1;
    b->entry_face = f;
}

int main(void)
{
    mv_camera c1, c2;
    const mv_camera *cams[2];
    double C2pos[3] = { 0.5, 0.0, 0.0 };
    double F[9];
    struct insect bugs[MAXB];
    struct track trk[MAXT];
    int gt_entry[6] = { 0 }, gt_exit[6] = { 0 };
    int est_entry[6] = { 0 }, est_exit[6] = { 0 };
    int nframes = (int)(DURATION * FPS);
    int spawned = 0, counted = 0, clutter_passed = 0;
    double se_pos = 0.0;
    long n_pos = 0;
    int f, i, j, k;

    memset(bugs, 0, sizeof(bugs));
    memset(trk, 0, sizeof(trk));

    mv_cam_set_K(&c1, 800.0, 800.0, 320.0, 240.0);
    mv_cam_set_identity_pose(&c1);
    memset(c1.k, 0, sizeof(c1.k));
    c2 = c1;
    mv_cam_set_pose_yaw(&c2, -6.0 * MV_PI / 180.0, C2pos);
    cams[0] = &c1;
    cams[1] = &c2;
    mv_fundamental_from_cams(F, &c1, &c2);

    for (f = 0; f < nframes; f++) {
        double dt = 1.0 / FPS;
        double det1[MAXD][2], det2[MAXD][2];
        int is_clutter1[MAXD], is_clutter2[MAXD];
        double X[MAXD][3];
        int used2[MAXD] = { 0 }, assigned[MAXD];
        int n1 = 0, n2 = 0, nx = 0;

        /* --- ground truth: spawn, move, exit --- */
        if (urand() < SPAWN_RATE * dt) {
            for (i = 0; i < MAXB; i++)
                if (!bugs[i].alive) {
                    spawn_insect(&bugs[i]);
                    gt_entry[bugs[i].entry_face]++;
                    spawned++;
                    break;
                }
        }
        for (i = 0; i < MAXB; i++) {
            struct insect *b = &bugs[i];
            double sp;
            if (!b->alive)
                continue;
            for (k = 0; k < 3; k++) {
                b->v[k] += (-0.8 * b->v[k]) * dt + 0.6 * sqrt(dt) * grand_();
                b->p[k] += b->v[k] * dt;
            }
            sp = mv_norm(b->v, 3);
            if (sp > 1.0)
                for (k = 0; k < 3; k++)
                    b->v[k] *= 1.0 / sp;
            if (b->p[0] < BX0 || b->p[0] > BX1 || b->p[1] < BY0
                || b->p[1] > BY1 || b->p[2] < BZ0 || b->p[2] > BZ1) {
                gt_exit[nearest_face(b->p)]++;
                b->alive = 0;
            }
        }

        /* --- detections: true insects + per-camera clutter --- */
        for (i = 0; i < MAXB && n1 < MAXD - 4 && n2 < MAXD - 4; i++) {
            double p1[2], p2[2];
            if (!bugs[i].alive || urand() > P_DET)
                continue;
            if (mv_cam_project(p1, &c1, bugs[i].p) != MV_OK || !in_img(p1))
                continue;
            if (mv_cam_project(p2, &c2, bugs[i].p) != MV_OK || !in_img(p2))
                continue;
            det1[n1][0] = p1[0] + 0.3 * grand_();
            det1[n1][1] = p1[1] + 0.3 * grand_();
            is_clutter1[n1] = 0;
            n1++;
            det2[n2][0] = p2[0] + 0.3 * grand_();
            det2[n2][1] = p2[1] + 0.3 * grand_();
            is_clutter2[n2] = 0;
            n2++;
        }
        for (k = 0; k < 4 && n1 < MAXD; k++)
            if (urand() < CLUTTER_PER_CAM / 4.0) {
                det1[n1][0] = 640.0 * urand();
                det1[n1][1] = 480.0 * urand();
                is_clutter1[n1] = 1;
                n1++;
            }
        for (k = 0; k < 4 && n2 < MAXD; k++)
            if (urand() < CLUTTER_PER_CAM / 4.0) {
                det2[n2][0] = 640.0 * urand();
                det2[n2][1] = 480.0 * urand();
                is_clutter2[n2] = 1;
                n2++;
            }

        /* --- stereo association by epipolar gate (greedy best pair) --- */
        for (i = 0; i < n1 && nx < MAXD; i++) {
            double best = 1e9;
            int bj = -1;
            for (j = 0; j < n2; j++) {
                double d;
                if (used2[j])
                    continue;
                d = mv_sym_epipolar_dist(F, det1[i], det2[j]);
                if (d < best) {
                    best = d;
                    bj = j;
                }
            }
            if (bj >= 0 && best < GATE_PX) {
                double uv[4];
                uv[0] = det1[i][0]; uv[1] = det1[i][1];
                uv[2] = det2[bj][0]; uv[3] = det2[bj][1];
                if (mv_triangulate(X[nx], cams, uv, 2) == MV_OK) {
                    used2[bj] = 1;
                    if (is_clutter1[i] || is_clutter2[bj])
                        clutter_passed++;
                    nx++;
                }
            }
        }

        /* --- track update: greedy NN association --- */
        for (i = 0; i < nx; i++)
            assigned[i] = 0;
        for (j = 0; j < MAXT; j++) {
            double best = 1e9;
            int bi = -1;
            if (!trk[j].alive)
                continue;
            for (i = 0; i < nx; i++) {
                double d[3], dist;
                if (assigned[i])
                    continue;
                for (k = 0; k < 3; k++)
                    d[k] = X[i][k] - trk[j].p[k];
                dist = mv_norm(d, 3);
                if (dist < best) {
                    best = dist;
                    bi = i;
                }
            }
            if (bi >= 0 && best < ASSOC_GATE) {
                assigned[bi] = 1;
                memcpy(trk[j].p, X[bi], sizeof(trk[j].p));
                memcpy(trk[j].last, X[bi], sizeof(trk[j].p));
                trk[j].hits++;
                trk[j].misses = 0;
                if (!trk[j].confirmed && trk[j].hits >= CONFIRM_HITS) {
                    trk[j].confirmed = 1;
                    counted++;
                    est_entry[nearest_face(trk[j].born)]++;
                }
            } else {
                trk[j].misses++;
                if (trk[j].misses >= KILL_MISSES) {
                    if (trk[j].confirmed)
                        est_exit[nearest_face(trk[j].last)]++;
                    trk[j].alive = 0;
                }
            }
        }
        for (i = 0; i < nx; i++) {
            if (assigned[i])
                continue;
            for (j = 0; j < MAXT; j++)
                if (!trk[j].alive) {
                    memset(&trk[j], 0, sizeof(trk[j]));
                    trk[j].alive = 1;
                    trk[j].hits = 1;
                    memcpy(trk[j].p, X[i], sizeof(trk[j].p));
                    memcpy(trk[j].born, X[i], sizeof(trk[j].p));
                    memcpy(trk[j].last, X[i], sizeof(trk[j].p));
                    break;
                }
        }

        /* --- accuracy of confirmed-track positions vs nearest truth --- */
        for (j = 0; j < MAXT; j++) {
            double best = 1e9;
            if (!trk[j].alive || !trk[j].confirmed || trk[j].misses)
                continue;
            for (i = 0; i < MAXB; i++) {
                double d[3];
                if (!bugs[i].alive)
                    continue;
                for (k = 0; k < 3; k++)
                    d[k] = trk[j].p[k] - bugs[i].p[k];
                if (mv_norm(d, 3) < best)
                    best = mv_norm(d, 3);
            }
            if (best < 1e8) {
                se_pos += best * best;
                n_pos++;
            }
        }
    }
    /* flush surviving confirmed tracks as exits at last position */
    for (j = 0; j < MAXT; j++)
        if (trk[j].alive && trk[j].confirmed)
            est_exit[nearest_face(trk[j].last)]++;

    printf("multiview insect-tracking experiment (seed 777)\n");
    printf("-----------------------------------------------\n");
    printf("duration            : %.0f s at %.0f fps\n", DURATION, FPS);
    printf("insects entered (gt): %d\n", spawned);
    printf("tracks confirmed    : %d\n", counted);
    printf("clutter pairs passing epipolar gate over %d frames: %d\n",
           nframes, clutter_passed);
    printf("track position RMS  : %.1f mm (%ld samples)\n",
           1000.0 * sqrt(se_pos / n_pos), n_pos);
    printf("\nface       entered(gt/est)   exited(gt/est)\n");
    for (i = 0; i < 6; i++)
        printf("%-9s      %2d / %-2d          %2d / %-2d\n", face_name[i],
               gt_entry[i], est_entry[i], gt_exit[i], est_exit[i]);
    return 0;
}
