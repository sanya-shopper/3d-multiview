/* Structured light on real captures: two cameras film a scene lit by a
 * projector running pattern.html mode M5 (Gray-code stripes, 46 slots
 * of ~0.5 s each, cycle ~23.5 s).  Output: the cameras' relative pose
 * from the dense code correspondences, and a triangulated point cloud.
 *
 * RUNBOOK (solo operation, no assistance needed):
 *  1. Connect the projector, open tools/pattern.html on it fullscreen,
 *     press 5.  Darken the room reasonably.
 *  2. Fix both cameras so their views overlap on the lit surfaces.
 *     Record >= 60 seconds of video on each (>= 2 full cycles).
 *  3. Calibrate each camera's intrinsics first if not already done
 *     (film the pattern in mode 1 or 4 separately; run calibreal).
 *  4. Extract frames at 0.1 s:
 *       /tmp/extract_frames camA.mov /tmp/slA 0.1
 *       /tmp/extract_frames camB.mov /tmp/slB 0.1
 *  5. Write a spec file per camera:
 *       line 1: fx fy cx cy k1 k2        (from calibreal)
 *       line 2: interval 0.1             (the extraction interval)
 *       rest:   the frame paths in order (ls /tmp/slA/f*.pgm)
 *  6. Run:  ./slreal <baseline_m> camA.spec camB.spec out
 *     where baseline_m is the camera-to-camera distance from a tape
 *     measure (or 0 for a unit-baseline result; scale is the only
 *     thing the codes cannot provide).
 *  7. Outputs: relative pose printout (compare the rotation angle and
 *     baseline direction against the physical setup), out.ply (open in
 *     MeshLab/CloudCompare/Blender), and per-stage diagnostics.  If
 *     "cycle detection" fails, the video was too short or the room too
 *     bright -- the white/black marker slots must be clearly the
 *     brightest/darkest frames.
 *
 * Build: make slreal */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mv/mv.h"

#define MAXFR 4096
#define NSLOTS 46
#define NBITS 11
#define MAXM 200000

typedef struct {
    double K[9], kr[2], interval;
    char paths[MAXFR][300];
    double mean[MAXFR];
    int nfr;
    /* per-slot representative frames (loaded on demand) */
    unsigned char *rep[NSLOTS];
    int w, h;
} slcam_t;

static int load_spec(slcam_t *c, const char *spec)
{
    FILE *fp = fopen(spec, "r");
    char line[512];
    if (!fp) {
        fprintf(stderr, "cannot open %s\n", spec);
        return MV_ERR;
    }
    memset(c, 0, sizeof(*c));
    if (!fgets(line, sizeof(line), fp)
        || sscanf(line, "%lf %lf %lf %lf %lf %lf", &c->K[0], &c->K[4],
                  &c->K[2], &c->K[5], &c->kr[0], &c->kr[1]) != 6) {
        fprintf(stderr, "%s: bad intrinsics line\n", spec);
        fclose(fp);
        return MV_ERR;
    }
    c->K[8] = 1.0;
    if (!fgets(line, sizeof(line), fp)
        || sscanf(line, "interval %lf", &c->interval) != 1) {
        fprintf(stderr, "%s: line 2 must be 'interval <s>'\n", spec);
        fclose(fp);
        return MV_ERR;
    }
    while (fgets(line, sizeof(line), fp) && c->nfr < MAXFR) {
        char *nl = strchr(line, '\n');
        if (nl)
            *nl = 0;
        if (line[0])
            snprintf(c->paths[c->nfr++], 300, "%s", line);
    }
    fclose(fp);
    return c->nfr > 100 ? MV_OK : MV_ERR;
}

/* pass 1: mean brightness per frame (subsampled) */
static int scan_means(slcam_t *c)
{
    int i;
    for (i = 0; i < c->nfr; i++) {
        unsigned char *img;
        int w, h, p, n = 0;
        long s = 0;
        if (mv_pgm_read(c->paths[i], &img, &w, &h) != MV_OK)
            return MV_ERR;
        c->w = w;
        c->h = h;
        for (p = 0; p < w * h; p += 7) {
            s += img[p];
            n++;
        }
        c->mean[i] = (double)s / n;
        free(img);
    }
    return MV_OK;
}

/* find cycle starts: a run of bright frames followed by a run of dark
 * frames (the WHITE and BLACK marker slots).  Returns slot duration in
 * frames and writes the first cycle-start frame index. */
static int find_cycle(const slcam_t *c, double *t0f, double *slotf)
{
    double hi = -1e9, lo = 1e9, thrh, thrl;
    double starts[32];
    int i, ns = 0;
    for (i = 0; i < c->nfr; i++) {
        if (c->mean[i] > hi)
            hi = c->mean[i];
        if (c->mean[i] < lo)
            lo = c->mean[i];
    }
    thrh = lo + 0.75 * (hi - lo);
    thrl = lo + 0.30 * (hi - lo);
    for (i = 1; i < c->nfr - 1 && ns < 32; i++) {
        /* bright frame whose successor region turns dark = W->B edge */
        if (c->mean[i] >= thrh && c->mean[i + 1] < thrh) {
            int j = i + 1, dark = 0;
            while (j < c->nfr && j < i + 10) {
                if (c->mean[j] <= thrl)
                    dark = 1;
                j++;
            }
            if (dark)
                starts[ns++] = i + 0.5; /* boundary between W and B */
        }
    }
    if (ns < 1)
        return MV_ERR;
    if (ns >= 2) {
        /* slot duration from the cycle period (NSLOTS slots/cycle) */
        double per = (starts[ns - 1] - starts[0]) / (ns - 1);
        *slotf = per / NSLOTS;
    } else {
        *slotf = 0.5 / c->interval; /* nominal 0.5 s slots */
    }
    /* starts[] marks the W|B boundary = end of slot 0 */
    *t0f = starts[0] - *slotf;
    printf("  cycle: %d marker(s), slot = %.2f frames (%.3f s)\n", ns,
           *slotf, *slotf * c->interval);
    return MV_OK;
}

static int load_reps(slcam_t *c)
{
    double t0f, slotf;
    int s;
    if (find_cycle(c, &t0f, &slotf) != MV_OK) {
        fprintf(stderr, "cycle detection failed\n");
        return MV_ERR;
    }
    for (s = 0; s < NSLOTS; s++) {
        int idx = (int)(t0f + (s + 0.5) * slotf + 0.5);
        int w, h;
        if (idx < 0 || idx >= c->nfr)
            return MV_ERR;
        if (mv_pgm_read(c->paths[idx], &c->rep[s], &w, &h) != MV_OK)
            return MV_ERR;
    }
    return MV_OK;
}

static void undistort_uv(double *out, const double *in, int n,
                         const double K[9], const double kr[2])
{
    double Kinv[9];
    int i, it;
    mv_mat_inv3(Kinv, K);
    for (i = 0; i < n; i++) {
        double u = in[2 * i], v = in[2 * i + 1];
        double xn = Kinv[0] * u + Kinv[1] * v + Kinv[2];
        double yn = Kinv[3] * u + Kinv[4] * v + Kinv[5];
        double x = xn, y = yn;
        for (it = 0; it < 20; it++) {
            double r2 = x * x + y * y;
            double rad = 1.0 + r2 * (kr[0] + r2 * kr[1]);
            if (fabs(rad) < 1e-6)
                break;
            x = xn / rad;
            y = yn / rad;
        }
        out[2 * i] = K[0] * x + K[1] * y + K[2];
        out[2 * i + 1] = K[4] * y + K[5];
    }
}

static int decode_cam(slcam_t *c, int **colc, int **rowc,
                      unsigned char **valid)
{
    const unsigned char *fr[NBITS], *fi[NBITS];
    unsigned char *vb;
    int npix = c->w * c->h, b, i;
    *colc = (int *)malloc((size_t)npix * sizeof(int));
    *rowc = (int *)malloc((size_t)npix * sizeof(int));
    *valid = (unsigned char *)malloc((size_t)npix);
    vb = (unsigned char *)malloc((size_t)npix);
    if (!*colc || !*rowc || !*valid || !vb)
        return MV_ERR;
    for (b = 0; b < NBITS; b++) {
        fr[b] = c->rep[2 + 2 * b];
        fi[b] = c->rep[2 + 2 * b + 1];
    }
    mv_graycode_decode(*colc, *valid, fr, fi, NBITS, npix, 12);
    for (b = 0; b < NBITS; b++) {
        fr[b] = c->rep[24 + 2 * b];
        fi[b] = c->rep[24 + 2 * b + 1];
    }
    mv_graycode_decode(*rowc, vb, fr, fi, NBITS, npix, 12);
    for (i = 0; i < npix; i++)
        (*valid)[i] = (*valid)[i] && vb[i];
    free(vb);
    return MV_OK;
}

int main(int argc, char **argv)
{
    static slcam_t A, B;
    double baseline;
    double *uv1, *uv2, *x1, *x2;
    int *colcA, *rowcA, *colcB, *rowcB;
    unsigned char *validA, *validB;
    double F[9], E[9], R[9], t[3], K1i[9], K2i[9];
    int nm, i;

    if (argc != 5) {
        fprintf(stderr, "usage: %s <baseline_m|0> <camA.spec> "
                        "<camB.spec> <out_prefix>\n", argv[0]);
        return 1;
    }
    baseline = atof(argv[1]);

    if (load_spec(&A, argv[2]) != MV_OK
        || load_spec(&B, argv[3]) != MV_OK)
        return 1;
    printf("camera A: %d frames; camera B: %d frames\n", A.nfr, B.nfr);
    if (scan_means(&A) != MV_OK || scan_means(&B) != MV_OK)
        return 1;
    printf("camera A:\n");
    if (load_reps(&A) != MV_OK)
        return 1;
    printf("camera B:\n");
    if (load_reps(&B) != MV_OK)
        return 1;

    if (decode_cam(&A, &colcA, &rowcA, &validA) != MV_OK
        || decode_cam(&B, &colcB, &rowcB, &validB) != MV_OK)
        return 1;
    {
        int nva = 0, nvb = 0, npixA = A.w * A.h, npixB = B.w * B.h;
        for (i = 0; i < npixA; i++)
            nva += validA[i];
        for (i = 0; i < npixB; i++)
            nvb += validB[i];
        printf("decoded pixels: A %d (%.1f%%), B %d (%.1f%%)\n", nva,
               100.0 * nva / npixA, nvb, 100.0 * nvb / npixB);
    }

    uv1 = (double *)malloc(2 * MAXM * sizeof(double));
    uv2 = (double *)malloc(2 * MAXM * sizeof(double));
    x1 = (double *)malloc(2 * MAXM * sizeof(double));
    x2 = (double *)malloc(2 * MAXM * sizeof(double));
    if (!uv1 || !uv2 || !x1 || !x2)
        return 1;
    nm = mv_graycode_match(uv1, uv2, MAXM, colcA, rowcA, validA, A.w,
                           A.h, colcB, rowcB, validB, B.w, B.h, 1920,
                           1080);
    printf("code correspondences: %d\n", nm);
    if (nm < 100) {
        fprintf(stderr, "too few matches -- views must overlap on lit "
                        "surfaces\n");
        return 1;
    }
    undistort_uv(uv1, uv1, nm, A.K, A.kr);
    undistort_uv(uv2, uv2, nm, B.K, B.kr);

    if (mv_fundamental_8point_robust(F, uv1, uv2, &nm) != MV_OK) {
        fprintf(stderr, "F estimation failed\n");
        return 1;
    }
    printf("after robust trim: %d correspondences\n", nm);
    mv_essential_from_fundamental(E, F, A.K, B.K);
    mv_mat_inv3(K1i, A.K);
    mv_mat_inv3(K2i, B.K);
    for (i = 0; i < nm; i++) {
        x1[2 * i] = K1i[0] * uv1[2 * i] + K1i[2];
        x1[2 * i + 1] = K1i[4] * uv1[2 * i + 1] + K1i[5];
        x2[2 * i] = K2i[0] * uv2[2 * i] + K2i[2];
        x2[2 * i + 1] = K2i[4] * uv2[2 * i + 1] + K2i[5];
    }
    if (mv_essential_pose(R, t, E, x1, x2, nm) != MV_OK) {
        fprintf(stderr, "essential pose failed\n");
        return 1;
    }
    {
        double tr = (R[0] + R[4] + R[8] - 1.0) / 2.0, ang;
        if (tr > 1.0)
            tr = 1.0;
        if (tr < -1.0)
            tr = -1.0;
        ang = acos(tr) * 180.0 / 3.14159265358979324;
        printf("\nrelative pose (camera B w.r.t. camera A):\n");
        printf("  rotation %.2f deg\n", ang);
        printf("  baseline direction t = (%+.4f %+.4f %+.4f)\n", t[0],
               t[1], t[2]);
        printf("  R = [%+.4f %+.4f %+.4f; %+.4f %+.4f %+.4f; "
               "%+.4f %+.4f %+.4f]\n", R[0], R[1], R[2], R[3], R[4],
               R[5], R[6], R[7], R[8]);
        if (baseline > 0.0)
            printf("  metric baseline applied: %.3f m\n", baseline);
        else {
            baseline = 1.0;
            printf("  NO baseline given: cloud is in unit-baseline "
                   "scale\n");
        }
    }

    /* triangulate all surviving correspondences */
    {
        mv_camera c1, c2;
        const mv_camera *cams[2];
        mv_cloud cloud;
        char path[512];
        int kept = 0;
        mv_cam_set_K(&c1, A.K[0], A.K[4], A.K[2], A.K[5]);
        mv_cam_set_identity_pose(&c1);
        memset(c1.k, 0, sizeof(c1.k));
        mv_cam_set_K(&c2, B.K[0], B.K[4], B.K[2], B.K[5]);
        memcpy(c2.R, R, sizeof(c2.R));
        for (i = 0; i < 3; i++)
            c2.t[i] = t[i] * baseline;
        memset(c2.k, 0, sizeof(c2.k));
        cams[0] = &c1;
        cams[1] = &c2;
        mv_cloud_init(&cloud, 1);
        for (i = 0; i < nm; i++) {
            double X[3], uvs[4];
            unsigned char rgb[3];
            double z2;
            uvs[0] = uv1[2 * i];
            uvs[1] = uv1[2 * i + 1];
            uvs[2] = uv2[2 * i];
            uvs[3] = uv2[2 * i + 1];
            if (mv_triangulate(X, cams, uvs, 2) != MV_OK)
                continue;
            z2 = c2.R[6] * X[0] + c2.R[7] * X[1] + c2.R[8] * X[2]
                 + c2.t[2];
            if (X[2] <= 0.0 || z2 <= 0.0)
                continue;
            {
                int px = (int)(uv1[2 * i] + 0.5),
                    py = (int)(uv1[2 * i + 1] + 0.5);
                unsigned char g = 200;
                if (px >= 0 && py >= 0 && px < A.w && py < A.h)
                    g = A.rep[0][py * A.w + px]; /* white frame = albedo */
                rgb[0] = rgb[1] = rgb[2] = g;
            }
            mv_cloud_push(&cloud, X, rgb);
            kept++;
        }
        snprintf(path, sizeof(path), "%s.ply", argv[4]);
        mv_cloud_write_ply(path, &cloud);
        printf("\ntriangulated %d points -> %s\n", kept, path);
        mv_cloud_free(&cloud);
    }
    return 0;
}
