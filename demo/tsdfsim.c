/* TSDF fusion validation: numbers quoted in doc/multiview.tex
 * ("Basic results").
 *
 * A sphere of known radius sits in the working volume. Both rig cameras
 * produce per-pixel depth of it with the rig's depth-law noise
 * sigma_Z = Z^2 sqrt(2) sigma_px / (f B); the noisy depth samples are
 * back-projected into enriched samples (point, origin, weight) and fused
 * into a TSDF; the zero level set is extracted by marching tetrahedra
 * and every mesh vertex is scored against the true sphere. Run at two
 * frame counts to verify the sqrt(n) averaging law (a T3-style scaling
 * tripwire). Deterministic.
 *
 * Run: make demo_tsdf && ./demo_tsdf
 * Out: out_tsdf.ply (mesh; open in MeshLab/CloudCompare/Blender) */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mv/mv.h"

#define MV_PI 3.14159265358979323846
#define CAMW 640
#define CAMH 480
#define SIGMA_PX 0.3
#define VOXEL 0.01
#define TAU 0.05

static const double CTR[3] = { 0.0, 0.1, 4.7 };
static const double RAD = 0.3;

static unsigned long long rng_state = 4242;

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

/* first ray-sphere intersection distance, or -1 */
static double hit_sphere(const double o[3], const double d[3])
{
    double oc[3], b, c, disc;
    int i;
    for (i = 0; i < 3; i++)
        oc[i] = o[i] - CTR[i];
    b = mv_dot(oc, d, 3);
    c = mv_dot(oc, oc, 3) - RAD * RAD;
    disc = b * b - c;
    if (disc < 0.0)
        return -1.0;
    return -b - sqrt(disc);
}

static double run_fusion(int nframes, int write_ply, int *ntri_out)
{
    mv_camera c1, c2;
    const mv_camera *cams[2];
    double C2pos[3] = { 0.5, 0.0, 0.0 };
    mv_tsdf t;
    double *tris;
    int ntri, f, ci, x, y, i;
    double se = 0.0;
    long ns = 0;

    mv_cam_set_K(&c1, 800.0, 800.0, 320.0, 240.0);
    mv_cam_set_identity_pose(&c1);
    memset(c1.k, 0, sizeof(c1.k));
    c2 = c1;
    mv_cam_set_pose_yaw(&c2, -6.0 * MV_PI / 180.0, C2pos);
    cams[0] = &c1;
    cams[1] = &c2;

    if (mv_tsdf_init(&t, CTR[0] - 0.45, CTR[1] - 0.45, CTR[2] - 0.45,
                     0.9, 0.9, 0.9, VOXEL, TAU) != MV_OK)
        return -1.0;

    for (f = 0; f < nframes; f++) {
        for (ci = 0; ci < 2; ci++) {
            double O[3];
            mv_cam_center(O, cams[ci]);
            /* sample every 2nd pixel for speed; dense enough */
            for (y = 0; y < CAMH; y += 2)
                for (x = 0; x < CAMW; x += 2) {
                    double uv[2] = { (double)x, (double)y };
                    double orig[3], dir[3], s, p[3], Zc, sigZ, w;
                    if (mv_cam_ray(orig, dir, cams[ci], uv) != MV_OK)
                        continue;
                    s = hit_sphere(orig, dir);
                    if (s <= 0.0)
                        continue;
                    /* depth-law noise applied along the ray */
                    for (i = 0; i < 3; i++)
                        p[i] = orig[i] + s * dir[i];
                    Zc = cams[ci]->R[6] * (p[0] - O[0])
                       + cams[ci]->R[7] * (p[1] - O[1])
                       + cams[ci]->R[8] * (p[2] - O[2]);
                    sigZ = Zc * Zc * sqrt(2.0) * SIGMA_PX / (800.0 * 0.5);
                    s += sigZ * grand_();
                    for (i = 0; i < 3; i++)
                        p[i] = orig[i] + s * dir[i];
                    w = 1.0 / (sigZ * sigZ);
                    mv_tsdf_fuse(&t, p, O, w);
                }
        }
    }

    if (mv_tsdf_mesh(&t, &tris, &ntri) != MV_OK) {
        mv_tsdf_free(&t);
        return -1.0;
    }
    for (i = 0; i < 3 * ntri; i++) {
        double d[3], e;
        int a;
        for (a = 0; a < 3; a++)
            d[a] = tris[3 * i + a] - CTR[a];
        e = mv_norm(d, 3) - RAD;
        se += e * e;
        ns++;
    }
    if (write_ply)
        mv_tsdf_write_ply(&t, "out_tsdf.ply");
    if (ntri_out)
        *ntri_out = ntri;
    free(tris);
    mv_tsdf_free(&t);
    return ns ? sqrt(se / ns) : -1.0;
}

int main(void)
{
    double r10, r160;
    int nt10, nt160;

    printf("multiview TSDF-fusion experiment (seed 4242)\n");
    printf("--------------------------------------------\n");
    printf("scene              : sphere r=%.2f m at depth %.1f m\n",
           RAD, CTR[2]);
    printf("depth noise        : depth law, sigma_px %.2f "
           "(sigma_Z ~ %.0f mm at sphere)\n", SIGMA_PX,
           1000.0 * 4.4 * 4.4 * sqrt(2.0) * SIGMA_PX / 400.0);
    printf("grid               : %.0f mm voxels, tau %.0f mm\n\n",
           VOXEL * 1000.0, TAU * 1000.0);

    r10 = run_fusion(10, 0, &nt10);
    printf("10 frames  : surface RMS %.2f mm  (%d triangles)\n",
           1000.0 * r10, nt10);
    r160 = run_fusion(160, 1, &nt160);
    printf("160 frames : surface RMS %.2f mm  (%d triangles)\n",
           1000.0 * r160, nt160);
    printf("\nscaling check      : ratio %.2f, sqrt(160/10) = 4.00\n",
           r10 / r160);
    printf("wrote out_tsdf.ply\n");
    return 0;
}
