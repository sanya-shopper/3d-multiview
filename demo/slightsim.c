/* Structured-light integration sim: fabricates the exact captures the
 * slreal tool expects -- two cameras filming a projector-lit scene,
 * frames extracted at 0.1 s, slots held 0.5 s, two full M5 cycles with
 * white/black marker slots -- and writes them as PGM files plus spec
 * files.  This validates the ENTIRE runbook path (timing recovery,
 * decode, match, pose, cloud) with known ground truth before any real
 * projector exists.
 *
 * Setup & run:
 *     make demo_slight && ./demo_slight
 *     ./slreal 0.6 /tmp/slsim/camA.spec /tmp/slsim/camB.spec /tmp/slsim/out
 * Expected: rotation ~6.9 deg, baseline direction ~(-0.98 0 +0.2),
 * tens of thousands of points in /tmp/slsim/out.ply shaped like two
 * walls.  Ground truth is printed below for comparison. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "mv/mv.h"

enum { CW = 1920, CH = 1080, IW = 480, IH = 360, TW = 500, TH = 340 };
enum { NBITS = 11, NSLOTS = 46, FPS = 10, SLOTF = 5 };

static unsigned char code[CW * CH];
static unsigned char texA[TW * TH], texB[TW * TH];
static unsigned char img[IW * IH];

static void set_plane(mv_plane *p, unsigned char *tex, double pitch,
                      const double r1[3], const double r2[3],
                      const double t[3])
{
    double r3[3];
    int i;
    mv_cross3(r3, r1, r2);
    for (i = 0; i < 3; i++) {
        p->R[i * 3 + 0] = r1[i];
        p->R[i * 3 + 1] = r2[i];
        p->R[i * 3 + 2] = r3[i];
        p->t[i] = t[i];
    }
    p->tex = tex;
    p->tw = TW;
    p->th = TH;
    p->pitch = pitch;
}

static void project_onto(unsigned char *tex, const mv_plane *pl,
                         const mv_camera *proj)
{
    int i, j;
    for (j = 0; j < pl->th; j++)
        for (i = 0; i < pl->tw; i++) {
            double acc = 0.0;
            int k, a, b, ns = 0;
            for (b = 0; b < 3; b++)
                for (a = 0; a < 3; a++) {
                    double X[3], uv[2];
                    int u, v;
                    for (k = 0; k < 3; k++)
                        X[k] = pl->R[k * 3 + 0]
                                   * (i + (a + 0.5) / 3.0) * pl->pitch
                             + pl->R[k * 3 + 1]
                                   * (j + (b + 0.5) / 3.0) * pl->pitch
                             + pl->t[k];
                    if (mv_cam_project(uv, proj, X) != MV_OK) {
                        acc += 8.0;
                        ns++;
                        continue;
                    }
                    u = (int)(uv[0] + 0.5);
                    v = (int)(uv[1] + 0.5);
                    acc += (u >= 0 && v >= 0 && u < CW && v < CH)
                           ? code[v * CW + u] : 8.0;
                    ns++;
                }
            tex[j * pl->tw + i] = (unsigned char)(acc / ns + 0.5);
        }
}

/* the M5 slot content, mirroring pattern.html exactly */
static void slot_code(int s)
{
    if (s == 0)
        memset(code, 255, sizeof(code));
    else if (s == 1)
        memset(code, 0, sizeof(code));
    else {
        int i = s - 2, axis = i < 22 ? 0 : 1, j = axis ? i - 22 : i;
        mv_graycode_frame(code, CW, CH, axis, j >> 1, NBITS, j & 1);
    }
}

int main(void)
{
    mv_camera proj, c1, c2;
    mv_plane pls[2];
    unsigned long long seed = 20260805ULL;
    FILE *fa, *fb;
    int slot, f, idx = 0, cyc;

    printf("structured-light integration sim (seed 20260805)\n");
    mkdir("/tmp/slsim", 0755);

    mv_cam_set_K(&proj, 1800.0, 1800.0, CW / 2.0, CH / 2.0);
    mv_cam_set_identity_pose(&proj);
    memset(proj.k, 0, sizeof(proj.k));
    mv_cam_set_K(&c1, 600.0, 600.0, IW / 2.0, IH / 2.0);
    mv_cam_set_identity_pose(&c1);
    memset(c1.k, 0, sizeof(c1.k));
    c1.t[0] = 0.2;
    c2 = c1;
    {
        double C2pos[3] = { 0.4, 0.0, 0.0 };
        mv_cam_set_pose_yaw(&c2, -0.12, C2pos);
    }
    {
        static const double ex[3] = { 1, 0, 0 }, ey[3] = { 0, 1, 0 };
        static const double exs[3] = { 0.9701, 0, -0.2425 };
        double t0[3];
        t0[0] = -1.7; t0[1] = -1.2; t0[2] = 3.4;
        set_plane(&pls[0], texA, 0.007, ex, ey, t0);
        t0[0] = 0.3; t0[1] = -1.2; t0[2] = 2.7;
        set_plane(&pls[1], texB, 0.007, exs, ey, t0);
    }
    {
        double R1t[9], Rrel[9], trel[3], tn;
        int j;
        mv_mat_transpose(R1t, c1.R, 3, 3);
        mv_mat_mul(Rrel, c2.R, R1t, 3, 3, 3);
        mv_mat_mul(trel, Rrel, c1.t, 3, 3, 1);
        for (j = 0; j < 3; j++)
            trel[j] = c2.t[j] - trel[j];
        tn = mv_norm(trel, 3);
        printf("ground truth: rotation %.2f deg, baseline %.3f m, "
               "direction (%+.3f %+.3f %+.3f)\n",
               0.12 * 180.0 / 3.14159265358979324, tn, trel[0] / tn,
               trel[1] / tn, trel[2] / tn);
    }

    fa = fopen("/tmp/slsim/camA.spec", "w");
    fb = fopen("/tmp/slsim/camB.spec", "w");
    if (!fa || !fb)
        return 1;
    fprintf(fa, "600 600 %d %d 0 0\ninterval 0.1\n", IW / 2, IH / 2);
    fprintf(fb, "600 600 %d %d 0 0\ninterval 0.1\n", IW / 2, IH / 2);

    /* a short pre-roll of code-slot light so the first cycle marker is
     * an interior transition, then two full cycles */
    for (cyc = 0; cyc < 2; cyc++)
        for (slot = 0; slot < NSLOTS; slot++) {
            slot_code(slot);
            project_onto(texA, &pls[0], &proj);
            project_onto(texB, &pls[1], &proj);
            for (f = 0; f < SLOTF; f++) {
                char p[128];
                mv_render_scene(img, NULL, IW, IH, &c1, pls, 2, 4, 1.0,
                                &seed);
                snprintf(p, sizeof(p), "/tmp/slsim/a%04d.pgm", idx);
                mv_pgm_write(p, img, IW, IH);
                fprintf(fa, "%s\n", p);
                mv_render_scene(img, NULL, IW, IH, &c2, pls, 2, 4, 1.0,
                                &seed);
                snprintf(p, sizeof(p), "/tmp/slsim/b%04d.pgm", idx);
                mv_pgm_write(p, img, IW, IH);
                fprintf(fb, "%s\n", p);
                idx++;
            }
        }
    fclose(fa);
    fclose(fb);
    printf("wrote %d frames per camera + specs under /tmp/slsim\n", idx);
    printf("now run:\n  ./slreal 0.6 /tmp/slsim/camA.spec "
           "/tmp/slsim/camB.spec /tmp/slsim/out\n");
    return 0;
}
