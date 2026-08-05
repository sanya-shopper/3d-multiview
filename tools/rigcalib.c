/* Multi-camera rig extrinsics from simultaneous views of the pattern.
 *
 * Each camera films the same (moving) display; every decoded frame
 * yields the display-plane pose in that camera AND the display's own
 * frame counter -- a shared clock.  Frames from two cameras whose
 * counters differ by <= MAXDK display refreshes (~50 ms) are treated
 * as simultaneous; each such pair gives one estimate of the relative
 * transform  x_a = R_ab x_b + t_ab  between the two cameras, and the
 * estimates are averaged (chordal mean rotation, median-trimmed mean
 * translation) with per-pair deviations reported as the quality gate.
 *
 * Usage: rigcalib <pitch_mm> <cam.spec> <cam.spec> [more.spec...]
 *   spec file: first line "fx fy cx cy k1 k2" (from calibreal),
 *   remaining lines: paths to that camera's frame PGMs.
 *   Camera 0 (first spec) anchors the rig frame.
 *
 * Build: make rigcalib */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mv/mv.h"

#define MAXCAM 8
#define MAXOBS 400
#define MAXDK 2 /* display refreshes: pairing tolerance (~33 ms) */

typedef struct {
    double K[9], kr[2];
    unsigned counter[MAXOBS];
    mv_camera pose[MAXOBS]; /* display-plane -> camera */
    int tier[MAXOBS];       /* 1 = fine pattern, 2 = coarse */
    int seq[MAXOBS];        /* frame order, for coarse wrap gating */
    int nobs;
} cam_t;

/* counter distance between observation i of camera a and j of camera b.
 * Fine counters compare directly; when either side is coarse (8-bit,
 * wraps every 256 refreshes ~ 4.3 s) compare mod 256 and gate by frame
 * order so a wrap cannot alias two distant moments into one. */
static long obs_dk(const cam_t *a, int i, const cam_t *b, int j)
{
    long dk;
    if (a->tier[i] == 1 && b->tier[j] == 1) {
        dk = (long)a->counter[i] - (long)b->counter[j];
    } else {
        long sd = a->seq[i] - b->seq[j];
        if (sd < -8 || sd > 8)
            return 1000;
        dk = ((long)(a->counter[i] & 255u) - (long)(b->counter[j] & 255u)
              + 384) % 256 - 128;
    }
    return dk < 0 ? -dk : dk;
}

static int undistort_uv(double *out, const double *in, int n,
                        const double K[9], const double kr[2])
{
    double Kinv[9];
    int i, it;
    if (mv_mat_inv3(Kinv, K) != MV_OK)
        return MV_ERR;
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
    return MV_OK;
}

static int load_camera(cam_t *c, const char *spec, double pitch)
{
    FILE *fp = fopen(spec, "r");
    char line[1024];
    int nread = 0;

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
    while (fgets(line, sizeof(line), fp) && c->nobs < MAXOBS) {
        unsigned char *img;
        int w, h, i;
        mv_read_result rr;
        double obj[2 * MV_READ_MAXC], und[2 * MV_READ_MAXC];
        char *nl = strchr(line, '\n');
        if (nl)
            *nl = 0;
        if (!line[0])
            continue;
        int tier = 1;
        nread++;
        if (mv_pgm_read(line, &img, &w, &h) != MV_OK)
            continue;
        if (mv_read_pattern(&rr, img, w, h) != MV_OK || rr.n < 100
            || !rr.counter_valid) {
            if (mv_read_coarse(&rr, img, w, h) != MV_OK || rr.n < 8
                || !rr.counter_valid) {
                free(img);
                continue;
            }
            tier = 2;
        }
        free(img);
        for (i = 0; i < rr.n; i++) {
            double xy[2];
            if (tier == 1)
                mv_pattern_corner_px(rr.id[i] % MV_PAT_CORNER_COLS,
                                     rr.id[i] / MV_PAT_CORNER_COLS, xy);
            else
                mv_pattern2_corner_px(rr.id[i] % MV_PAT2_CORNER_COLS,
                                      rr.id[i] / MV_PAT2_CORNER_COLS,
                                      xy);
            obj[2 * i] = xy[0] * pitch;
            obj[2 * i + 1] = xy[1] * pitch;
        }
        if (undistort_uv(und, rr.uv, rr.n, c->K, c->kr) != MV_OK)
            continue;
        if (mv_calib_plane_pose(&c->pose[c->nobs], c->K, obj, und, rr.n)
            != MV_OK)
            continue;
        c->counter[c->nobs] = rr.counter;
        c->tier[c->nobs] = tier;
        c->seq[c->nobs] = nread;
        c->nobs++;
    }
    fclose(fp);
    printf("%s: %d frames read, %d posed with valid counter\n", spec,
           nread, c->nobs);
    return c->nobs >= 3 ? MV_OK : MV_ERR;
}

static int cmp_dbl(const void *a, const void *b)
{
    double d = *(const double *)a - *(const double *)b;
    return d < 0 ? -1 : d > 0 ? 1 : 0;
}

/* relative transform of pair (a_obs, b_obs): x_a = R x_b + t */
static void rel_pose(double R[9], double t[3], const mv_camera *a,
                     const mv_camera *b)
{
    double Rbt[9];
    int i;
    mv_mat_transpose(Rbt, b->R, 3, 3);
    mv_mat_mul(R, a->R, Rbt, 3, 3, 3);
    mv_mat_mul(t, R, b->t, 3, 3, 1);
    for (i = 0; i < 3; i++)
        t[i] = a->t[i] - t[i];
}

int main(int argc, char **argv)
{
    static cam_t cams[MAXCAM];
    double pitch;
    int ncam, c;

    if (argc < 4) {
        fprintf(stderr, "usage: %s <pitch_mm> <cam.spec> <cam.spec> "
                        "[more...]\n", argv[0]);
        return 1;
    }
    pitch = atof(argv[1]) * 1e-3;
    ncam = argc - 2;
    if (ncam > MAXCAM)
        ncam = MAXCAM;
    for (c = 0; c < ncam; c++)
        if (load_camera(&cams[c], argv[c + 2], pitch) != MV_OK)
            return 1;

    /* counter-rate sanity: median counter step between successive
     * observations, per camera */
    for (c = 0; c < ncam; c++) {
        double steps[MAXOBS];
        int i, ns = 0;
        for (i = 1; i < cams[c].nobs; i++)
            if (cams[c].counter[i] > cams[c].counter[i - 1])
                steps[ns++] = (double)(cams[c].counter[i]
                                       - cams[c].counter[i - 1]);
        if (ns) {
            qsort(steps, (size_t)ns, sizeof(double), cmp_dbl);
            printf("camera %d: median counter step %.0f refreshes "
                   "between samples\n", c, steps[ns / 2]);
        }
    }

    for (c = 1; c < ncam; c++) {
        double Rsum[9] = { 0 }, Racc[9], U[9], S[3], V[9], Vt[9];
        double R[9], t[3], tacc[3][256], tmean[3] = { 0, 0, 0 };
        double angdev[256], tdev[256];
        int npair = 0, i, j, k, n2 = 0;

        for (i = 0; i < cams[0].nobs && npair < 256; i++) {
            int bestj = -1;
            long bestdk = MAXDK + 1;
            for (j = 0; j < cams[c].nobs; j++) {
                long dk = obs_dk(&cams[0], i, &cams[c], j);
                if (dk < bestdk) {
                    bestdk = dk;
                    bestj = j;
                }
            }
            if (bestj < 0 || bestdk > MAXDK)
                continue;
            rel_pose(R, t, &cams[0].pose[i], &cams[c].pose[bestj]);
            for (k = 0; k < 9; k++)
                Rsum[k] += R[k];
            for (k = 0; k < 3; k++)
                tacc[k][npair] = t[k];
            npair++;
        }
        printf("\npair (0,%d): %d simultaneous observations "
               "(|dk| <= %d)\n", c, npair, MAXDK);
        if (npair < 3) {
            printf("  too few simultaneous views -- cannot solve\n");
            continue;
        }
        /* chordal mean rotation: project the summed matrix to SO(3).
         * U*Vt is the nearest ORTHOGONAL matrix -- a rotation only when
         * its determinant is +1; if det(Rsum) < 0 (inconsistent pairs)
         * the bare product is a reflection, so flip the last column of
         * U, exactly as the pose solver in src/calib.c guards. */
        memcpy(U, Rsum, sizeof(Rsum));
        if (mv_svd(U, S, V, 3, 3) != MV_OK)
            continue;
        mv_mat_transpose(Vt, V, 3, 3);
        mv_mat_mul(Racc, U, Vt, 3, 3, 3);
        {
            double det = Racc[0] * (Racc[4] * Racc[8] - Racc[5] * Racc[7])
                       - Racc[1] * (Racc[3] * Racc[8] - Racc[5] * Racc[6])
                       + Racc[2] * (Racc[3] * Racc[7] - Racc[4] * Racc[6]);
            if (det < 0.0) {
                int r;
                for (r = 0; r < 3; r++)
                    U[r * 3 + 2] = -U[r * 3 + 2];
                mv_mat_mul(Racc, U, Vt, 3, 3, 3);
            }
        }
        /* median-trimmed translation: per-axis median */
        for (k = 0; k < 3; k++) {
            double col[256];
            memcpy(col, tacc[k], (size_t)npair * sizeof(double));
            qsort(col, (size_t)npair, sizeof(double), cmp_dbl);
            tmean[k] = col[npair / 2];
        }
        /* per-pair deviations from the mean */
        for (i = 0; i < npair; i++) {
            double dt = 0.0;
            for (k = 0; k < 3; k++)
                dt += (tacc[k][i] - tmean[k]) * (tacc[k][i] - tmean[k]);
            tdev[i] = sqrt(dt);
        }
        /* recompute rotations for deviation stats */
        {
            for (i = 0; i < cams[0].nobs && n2 < npair; i++) {
                int bestj = -1;
                long bestdk = MAXDK + 1;
                for (j = 0; j < cams[c].nobs; j++) {
                    long dk = obs_dk(&cams[0], i, &cams[c], j);
                    if (dk < bestdk) {
                        bestdk = dk;
                        bestj = j;
                    }
                }
                if (bestj < 0 || bestdk > MAXDK)
                    continue;
                rel_pose(R, t, &cams[0].pose[i], &cams[c].pose[bestj]);
                {
                    double Rt[9], D[9], tr;
                    mv_mat_transpose(Rt, Racc, 3, 3);
                    mv_mat_mul(D, R, Rt, 3, 3, 3);
                    tr = D[0] + D[4] + D[8];
                    tr = (tr - 1.0) / 2.0;
                    if (tr > 1.0)
                        tr = 1.0;
                    if (tr < -1.0)
                        tr = -1.0;
                    angdev[n2] = acos(tr) * 180.0 / 3.14159265358979324;
                }
                n2++;
            }
        }
        /* n2 == npair by construction (identical matching); index by
         * n2 so the analyzer sees only written entries reach qsort */
        qsort(angdev, (size_t)n2, sizeof(double), cmp_dbl);
        qsort(tdev, (size_t)npair, sizeof(double), cmp_dbl);
        {
            double ax, ay, az, ang, tr;
            tr = (Racc[0] + Racc[4] + Racc[8] - 1.0) / 2.0;
            if (tr > 1.0)
                tr = 1.0;
            if (tr < -1.0)
                tr = -1.0;
            ang = acos(tr) * 180.0 / 3.14159265358979324;
            ax = Racc[7] - Racc[5];
            ay = Racc[2] - Racc[6];
            az = Racc[3] - Racc[1];
            printf("  baseline |t| = %.3f m   t = (%+.3f %+.3f %+.3f)\n",
                   sqrt(tmean[0] * tmean[0] + tmean[1] * tmean[1]
                        + tmean[2] * tmean[2]),
                   tmean[0], tmean[1], tmean[2]);
            printf("  relative rotation %.1f deg about axis "
                   "(%+.2f %+.2f %+.2f)\n", ang, ax, ay, az);
            printf("  consistency: median |dR| %.2f deg, median |dt| "
                   "%.1f mm over %d pairs\n", angdev[(n2 > 0 ? n2 : 1) / 2],
                   1000.0 * tdev[npair / 2], n2 > 0 ? n2 : npair);
            printf("  R = [%+.4f %+.4f %+.4f; %+.4f %+.4f %+.4f; "
                   "%+.4f %+.4f %+.4f]\n", Racc[0], Racc[1], Racc[2],
                   Racc[3], Racc[4], Racc[5], Racc[6], Racc[7], Racc[8]);
        }
    }
    return 0;
}
