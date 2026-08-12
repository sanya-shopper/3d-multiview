/* densereal: the dense pipeline at REAL tier -- rectification, SAD
 * block matching, texture gating, TSDF fusion and mesh extraction over
 * real recorded frames (roadmap: "first real-image room model" step).
 *
 * The recorded sessions are single-camera, but the 2026-08-06 session
 * sweeps the camera/display distance through ~8 cm, so frame PAIRS
 * carry a genuine accidental baseline.  The calibration display in
 * view anchors everything: each frame's decoded pattern gives that
 * frame's metric camera pose in the display frame (mv_calib_plane_pose
 * at the known pixel pitch), so the two-frame rig needs no essential
 * matrix and no scale guess -- the display is ruler and world origin.
 *
 * Ground truth for free: the display surface IS the plane z = 0 of the
 * reconstruction frame, so every reconstructed point that projects
 * inside the decoded pattern quad scores its own |z| as a flatness
 * residual.  That number is the real-tier answer to the sim-tier
 * plane-RMS results of demo_room.
 *
 *   usage: densereal <pitch_mm> <fx> <fy> <cx> <cy> <k1> <k2>
 *                    <out_prefix> <frames.pgm...>
 *
 * Calibrate first (tools/calibreal on the same frames) and pass its
 * K/k here.  The tool decodes every frame, picks the pair with the
 * largest LATERAL accidental baseline (depth-only baselines are
 * useless for stereo), and runs the full dense chain on it. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mv/mv.h"
#include "mv/plane.h"

#define MV_PI 3.14159265358979323846
#define MAXFRAMES 256
#define HWIN 4              /* SAD window (2*4+1)^2 */
#define TEXTURE_MIN 6.0     /* mean |grad| gate, as in demo_room */
#define VOXEL 0.004         /* m */
#define TAU   0.012         /* m */

typedef struct {
    const char *path;
    mv_camera cam;          /* pose in the display frame (undistorted) */
    double C[3];            /* camera center */
    double uvmin[2], uvmax[2]; /* decoded-corner bbox, raw pixels */
    int ncorners;
} frameinfo;

/* iterative undistortion of pixel points (same as the hub's) */
static void undo_distort(double *out, const double *in, int n,
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

/* full-resolution undistorted copy: output pixel (u,v) samples the
 * source at the forward-distorted position (no iteration needed) */
static unsigned char *undistort_img(const unsigned char *img, int w,
                                    int h, const double K[9],
                                    const double kr[2])
{
    unsigned char *out = malloc((size_t)w * h);
    double Kinv[9];
    int x, y;
    if (!out)
        return NULL;
    mv_mat_inv3(Kinv, K);
    for (y = 0; y < h; y++)
        for (x = 0; x < w; x++) {
            double xn = Kinv[0] * x + Kinv[1] * y + Kinv[2];
            double yn = Kinv[4] * y + Kinv[5];
            double r2 = xn * xn + yn * yn;
            double rad = 1.0 + r2 * (kr[0] + r2 * kr[1]);
            double u = K[0] * xn * rad + K[2];
            double v = K[4] * yn * rad + K[5];
            int u0 = (int)floor(u), v0 = (int)floor(v);
            double fu = u - u0, fv = v - v0;
            unsigned char g = 0;
            if (u0 >= 0 && v0 >= 0 && u0 + 1 < w && v0 + 1 < h) {
                double g00 = img[v0 * w + u0];
                double g10 = img[v0 * w + u0 + 1];
                double g01 = img[(v0 + 1) * w + u0];
                double g11 = img[(v0 + 1) * w + u0 + 1];
                g = (unsigned char)((1 - fv) * ((1 - fu) * g00 + fu * g10)
                                    + fv * ((1 - fu) * g01 + fu * g11)
                                    + 0.5);
            }
            out[y * w + x] = g;
        }
    return out;
}

/* plane-guided block matching: search each left pixel only within
 * +/- band of its EXPECTED display-plane disparity (from the z = 0
 * homography between the rectified views).  The pattern repeats every
 * ~26 px; a global range spanning the tilted display's ~90 px of
 * disparity lets SAD lock whole periods off (measured: 23.5 mm RMS
 * flatness); a band under half a period makes that impossible. */
static void guided_sad(float *disp, const unsigned char *L,
                       const unsigned char *R, int w, int h,
                       const double Hlr[9], int band)
{
    int x, y, a, b, d;
    for (y = 0; y < h; y++)
        for (x = 0; x < w; x++)
            disp[y * w + x] = -1.0f;
    for (y = HWIN; y < h - HWIN; y++)
        for (x = HWIN; x < w - HWIN; x++) {
            double dn = Hlr[6] * x + Hlr[7] * y + Hlr[8];
            double ux, dexp, best = 1e30, c0 = 1e30, cm = 1e30,
                   cp = 1e30;
            int dlo, dhi, dbest = -1;
            if (fabs(dn) < 1e-12)
                continue;
            ux = (Hlr[0] * x + Hlr[1] * y + Hlr[2]) / dn;
            dexp = x - ux;
            dlo = (int)floor(dexp) - band;
            dhi = (int)ceil(dexp) + band;
            for (d = dlo; d <= dhi; d++) {
                double c = 0.0;
                if (x - d - HWIN < 0 || x - d + HWIN >= w)
                    continue;
                for (b = -HWIN; b <= HWIN; b++)
                    for (a = -HWIN; a <= HWIN; a++)
                        c += fabs((double)L[(y + b) * w + x + a]
                                  - R[(y + b) * w + x - d + a]);
                if (c < best) {
                    best = c;
                    dbest = d;
                }
            }
            if (dbest <= dlo || dbest >= dhi)
                continue;
            /* parabolic sub-pixel refinement around the minimum */
            {
                double c;
                int dd;
                for (dd = dbest - 1; dd <= dbest + 1; dd++) {
                    if (x - dd - HWIN < 0 || x - dd + HWIN >= w)
                        continue;
                    c = 0.0;
                    for (b = -HWIN; b <= HWIN; b++)
                        for (a = -HWIN; a <= HWIN; a++)
                            c += fabs((double)L[(y + b) * w + x + a]
                                      - R[(y + b) * w + x - dd + a]);
                    if (dd < dbest)
                        cm = c;
                    else if (dd > dbest)
                        cp = c;
                    else
                        c0 = c;
                }
                if (cm < 1e30 && cp < 1e30
                    && cm + cp - 2.0 * c0 > 1e-9)
                    disp[y * w + x] = (float)(dbest
                        + 0.5 * (cm - cp) / (cm + cp - 2.0 * c0));
                else
                    disp[y * w + x] = (float)dbest;
            }
        }
}

/* invalidate disparities where the left window has too little texture
 * (demo_room's texture-desert gate) */
static void texture_gate(float *disp, const unsigned char *img, int w,
                         int h)
{
    int x, y, a, b;
    for (y = HWIN; y < h - HWIN; y++)
        for (x = HWIN; x < w - HWIN; x++) {
            double g = 0.0;
            for (b = -HWIN; b <= HWIN; b++)
                for (a = -HWIN; a < HWIN; a++)
                    g += fabs((double)img[(y + b) * w + x + a + 1]
                              - img[(y + b) * w + x + a]);
            g /= (2 * HWIN + 1) * (2 * HWIN);
            if (g < TEXTURE_MIN)
                disp[y * w + x] = -1.0f;
        }
}

static void write_norm_pgm(const char *path, const float *d, int w,
                           int h, double lo, double hi)
{
    unsigned char *img = malloc((size_t)w * h);
    int i;
    if (!img)
        return;
    for (i = 0; i < w * h; i++) {
        double v = d[i];
        img[i] = v < 0.0 ? 0
                 : (unsigned char)(20.0 + 235.0
                                   * ((v < lo ? lo : v > hi ? hi : v) - lo)
                                   / (hi - lo));
    }
    mv_pgm_write(path, img, w, h);
    free(img);
}

/* decode one frame and recover its metric pose in the display frame */
static int frame_pose(frameinfo *fi, const unsigned char *img, int w,
                      int h, const double K[9], const double kr[2],
                      double pitch)
{
    mv_read_result rr;
    static double obj[2 * MV_READ_MAXC], und[2 * MV_READ_MAXC];
    int i;
    if (mv_read_pattern(&rr, img, w, h) != MV_OK || rr.n < 40)
        return MV_ERR;
    for (i = 0; i < rr.n; i++) {
        double xy[2];
        mv_pattern_corner_px(rr.id[i] % MV_PAT_CORNER_COLS,
                             rr.id[i] / MV_PAT_CORNER_COLS, xy);
        obj[2 * i] = xy[0] * pitch;
        obj[2 * i + 1] = xy[1] * pitch;
    }
    undo_distort(und, rr.uv, rr.n, K, kr);
    if (mv_calib_plane_pose(&fi->cam, K, obj, und, rr.n) != MV_OK)
        return MV_ERR;
    memset(fi->cam.k, 0, sizeof(fi->cam.k));
    mv_cam_center(fi->C, &fi->cam);
    fi->uvmin[0] = fi->uvmin[1] = 1e9;
    fi->uvmax[0] = fi->uvmax[1] = -1e9;
    for (i = 0; i < rr.n; i++) {
        double u = rr.uv[2 * i], v = rr.uv[2 * i + 1];
        if (u < fi->uvmin[0]) fi->uvmin[0] = u;
        if (v < fi->uvmin[1]) fi->uvmin[1] = v;
        if (u > fi->uvmax[0]) fi->uvmax[0] = u;
        if (v > fi->uvmax[1]) fi->uvmax[1] = v;
    }
    fi->ncorners = rr.n;
    return MV_OK;
}

int main(int argc, char **argv)
{
    double pitch, K[9], kr[2];
    const char *prefix;
    static frameinfo fr[MAXFRAMES];
    int nfr = 0, i, j, x, y, w = 0, h = 0;
    int besti = -1, bestj = -1;
    double bestlat = 0.0;

    if (argc < 11) {
        fprintf(stderr, "usage: %s <pitch_mm> <fx> <fy> <cx> <cy> "
                        "<k1> <k2> <out_prefix> <frames.pgm...>\n",
                argv[0]);
        return 1;
    }
    pitch = atof(argv[1]) * 1e-3;
    memset(K, 0, sizeof(K));
    K[0] = atof(argv[2]);
    K[4] = atof(argv[3]);
    K[2] = atof(argv[4]);
    K[5] = atof(argv[5]);
    K[8] = 1.0;
    kr[0] = atof(argv[6]);
    kr[1] = atof(argv[7]);
    prefix = argv[8];

    /* ---- decode every frame; keep the ones with a solid pose ------- */
    for (i = 9; i < argc && nfr < MAXFRAMES; i++) {
        unsigned char *img;
        int fw, fh;
        if (mv_pgm_read(argv[i], &img, &fw, &fh) != MV_OK)
            continue;
        if (w == 0) {
            w = fw;
            h = fh;
        }
        if (fw == w && fh == h) {
            fr[nfr].path = argv[i];
            if (frame_pose(&fr[nfr], img, w, h, K, kr, pitch) == MV_OK)
                nfr++;
        }
        free(img);
    }
    fprintf(stderr, "[densereal] %d/%d frames decoded to poses\n",
            nfr, argc - 9);
    if (nfr < 2)
        return 1;

    /* ---- pick the best USABLE pair: the largest lateral baseline
     * whose expected display disparity stays in SAD range and whose
     * inter-frame rotation is small enough for appearance matching
     * (a 30 cm baseline at 35 cm range is geometrically wonderful and
     * photometrically hopeless) ---------------------------------- */
    for (i = 0; i < nfr; i++)
        for (j = i + 1; j < nfr; j++) {
            double bvec[3], zax[3], dotp = 0.0, lat2 = 0.0;
            double zdist, dexp, tr, ang;
            int k;
            for (k = 0; k < 3; k++)
                bvec[k] = fr[j].C[k] - fr[i].C[k];
            for (k = 0; k < 3; k++)
                zax[k] = fr[i].cam.R[6 + k]; /* optical axis, world */
            for (k = 0; k < 3; k++)
                dotp += bvec[k] * zax[k];
            for (k = 0; k < 3; k++) {
                double l = bvec[k] - dotp * zax[k];
                lat2 += l * l;
            }
            zdist = sqrt(fr[i].C[0] * fr[i].C[0]
                         + fr[i].C[1] * fr[i].C[1]
                         + fr[i].C[2] * fr[i].C[2]);
            dexp = K[0] * sqrt(lat2) / (zdist > 0.05 ? zdist : 0.05);
            if (dexp < 15.0 || dexp > 110.0)
                continue;
            /* relative rotation angle between the two frames */
            tr = 0.0;
            for (k = 0; k < 3; k++) {
                int m;
                for (m = 0; m < 3; m++)
                    tr += fr[i].cam.R[3 * k + m] * fr[j].cam.R[3 * k + m];
            }
            ang = acos(((tr - 1.0) / 2.0) > 1.0 ? 1.0
                       : ((tr - 1.0) / 2.0) < -1.0 ? -1.0
                       : (tr - 1.0) / 2.0);
            if (ang > 10.0 * MV_PI / 180.0)
                continue;
            if (lat2 > bestlat) {
                bestlat = lat2;
                besti = i;
                bestj = j;
            }
        }
    if (besti < 0) {
        fprintf(stderr, "[densereal] no usable pair (lateral baseline "
                        "vs range vs rotation)\n");
        return 1;
    }
    {
        double B;
        frameinfo *A = &fr[besti], *Bf = &fr[bestj];
        unsigned char *rawL, *rawR, *undL, *undR, *rectL, *rectR;
        int fw, fh;
        mv_camera r1, r2;
        double H1[9], H2[9], Kinv[9], C1[3], fnew, Bnew;
        float *disp, *depth;
        int dmin, dmax, nvalid = 0;
        mv_cloud cloud;
        mv_tsdf t;
        double zsum = 0.0, z2sum = 0.0;
        long nz = 0;
        double lo[3] = { 1e9, 1e9, 1e9 }, hi[3] = { -1e9, -1e9, -1e9 };

        B = sqrt((Bf->C[0] - A->C[0]) * (Bf->C[0] - A->C[0])
                 + (Bf->C[1] - A->C[1]) * (Bf->C[1] - A->C[1])
                 + (Bf->C[2] - A->C[2]) * (Bf->C[2] - A->C[2]));
        fprintf(stderr, "[densereal] pair: %s + %s\n", A->path, Bf->path);
        fprintf(stderr, "[densereal] accidental baseline %.1f mm "
                        "(lateral %.1f mm)\n",
                B * 1e3, sqrt(bestlat) * 1e3);

        if (mv_pgm_read(A->path, &rawL, &fw, &fh) != MV_OK
            || mv_pgm_read(Bf->path, &rawR, &fw, &fh) != MV_OK)
            return 1;
        undL = undistort_img(rawL, w, h, K, kr);
        undR = undistort_img(rawR, w, h, K, kr);
        if (!undL || !undR)
            return 1;

        if (mv_rectify_pair(&A->cam, &Bf->cam, &r1, &r2, H1, H2)
            != MV_OK) {
            fprintf(stderr, "[densereal] rectification failed\n");
            return 1;
        }
        /* re-center: the rectifying rotation can swing the content
         * off the canvas entirely.  Shift the (shared) rectified
         * principal point so the display center lands mid-canvas --
         * the SAME shift for both cameras, so rows stay aligned and
         * disparity is preserved. */
        {
            double Xc[3], ua[2], ub[2], dx, dy, T[9], Hn[9];
            Xc[0] = MV_PAT_W / 2.0 * pitch;
            Xc[1] = MV_PAT_H / 2.0 * pitch;
            Xc[2] = 0.0;
            if (mv_cam_project(ua, &r1, Xc) == MV_OK
                && mv_cam_project(ub, &r2, Xc) == MV_OK) {
                dx = w / 2.0 - (ua[0] + ub[0]) / 2.0;
                dy = h / 2.0 - (ua[1] + ub[1]) / 2.0;
                r1.K[2] += dx;
                r1.K[5] += dy;
                r2.K[2] += dx;
                r2.K[5] += dy;
                memset(T, 0, sizeof(T));
                T[0] = T[4] = T[8] = 1.0;
                T[2] = dx;
                T[5] = dy;
                mv_mat_mul(Hn, T, H1, 3, 3, 3);
                memcpy(H1, Hn, sizeof(Hn));
                mv_mat_mul(Hn, T, H2, 3, 3, 3);
                memcpy(H2, Hn, sizeof(Hn));
            }
        }
        rectL = malloc((size_t)w * h);
        rectR = malloc((size_t)w * h);
        disp = malloc((size_t)w * h * sizeof(float));
        depth = malloc((size_t)w * h * sizeof(float));
        if (!rectL || !rectR || !disp || !depth)
            return 1;
        mv_warp_homography(rectL, w, h, undL, w, h, H1, 0);
        mv_warp_homography(rectR, w, h, undR, w, h, H2, 0);
        fnew = r1.K[0];
        Bnew = mv_baseline(&r1, &r2);

        /* disparity search range measured, not guessed: project the
         * display quad's corners (z = 0 plane, known metric extent)
         * through both rectified cameras and bracket their true
         * disparities with headroom for scene in front and behind */
        {
            double dlo = 1e9, dhi = -1e9;
            double qx[4] = { 0.0, MV_PAT_W * pitch, MV_PAT_W * pitch,
                             0.0 };
            double qy[4] = { 0.0, 0.0, MV_PAT_H * pitch,
                             MV_PAT_H * pitch };
            int k;
            for (k = 0; k < 4; k++) {
                double Xw[3], ua[2], ub[2];
                Xw[0] = qx[k];
                Xw[1] = qy[k];
                Xw[2] = 0.0;
                if (mv_cam_project(ua, &r1, Xw) != MV_OK
                    || mv_cam_project(ub, &r2, Xw) != MV_OK)
                    continue;
                if (ua[0] - ub[0] < dlo)
                    dlo = ua[0] - ub[0];
                if (ua[0] - ub[0] > dhi)
                    dhi = ua[0] - ub[0];
            }
            /* TIGHT band around the display's measured disparities:
             * the pattern is periodic (~26 px), and a wide search lets
             * SAD lock one period off -- the doc's repetitive-texture
             * failure mode, worth a ~45 mm depth error here.  A band
             * narrower than half a period makes aliasing impossible
             * where it matters; scene far outside the band is simply
             * not reconstructed in this pass. */
            dmin = (int)floor(dlo) - 15;
            dmax = (int)ceil(dhi) + 15;
            if (dmin < 0)
                dmin = 0;
            if (dmax > w / 2)
                dmax = w / 2;
            if (dmax <= dmin) {
                fprintf(stderr, "[densereal] degenerate disparity "
                                "range\n");
                return 1;
            }
        }
        fprintf(stderr, "[densereal] rectified f %.1f px, B %.1f mm, "
                        "disparity range [%d, %d]\n",
                fnew, Bnew * 1e3, dmin, dmax);

        {
            char path[512];
            snprintf(path, sizeof(path), "%s_rectL.pgm", prefix);
            mv_pgm_write(path, rectL, w, h);
            snprintf(path, sizeof(path), "%s_rectR.pgm", prefix);
            mv_pgm_write(path, rectR, w, h);
        }
        /* homography of the display plane (z = 0) between the two
         * rectified views: H_wl = K [r1 r2 t], Hlr = H_wr * H_wl^-1 */
        {
            double Hwl[9], Hwr[9], Hwli[9], Hlr[9];
            int k;
            memset(Hwl, 0, sizeof(Hwl));
            for (k = 0; k < 3; k++) {
                int m;
                for (m = 0; m < 3; m++) {
                    double col0 = r1.R[3 * m], col1 = r1.R[3 * m + 1];
                    double tcm = r1.t[m];
                    Hwl[3 * k] += r1.K[3 * k + m] * col0;
                    Hwl[3 * k + 1] += r1.K[3 * k + m] * col1;
                    Hwl[3 * k + 2] += r1.K[3 * k + m] * tcm;
                }
            }
            memset(Hwr, 0, sizeof(Hwr));
            for (k = 0; k < 3; k++) {
                int m;
                for (m = 0; m < 3; m++) {
                    Hwr[3 * k] += r2.K[3 * k + m] * r2.R[3 * m];
                    Hwr[3 * k + 1] += r2.K[3 * k + m] * r2.R[3 * m + 1];
                    Hwr[3 * k + 2] += r2.K[3 * k + m] * r2.t[m];
                }
            }
            if (mv_mat_inv3(Hwli, Hwl) != MV_OK)
                return 1;
            mv_mat_mul(Hlr, Hwr, Hwli, 3, 3, 3);
            /* guided matching wants left->right: invert the naming --
             * Hlr here maps LEFT rectified pixels onto RIGHT ones */
            guided_sad(disp, rectL, rectR, w, h, Hlr, 12);
        }
        texture_gate(disp, rectL, w, h);

        mv_mat_inv3(Kinv, r1.K);
        mv_cam_center(C1, &r1);
        mv_cloud_init(&cloud, 1);

        for (i = 0; i < w * h; i++)
            depth[i] = disp[i] > 0.0f
                       ? (float)(fnew * Bnew / disp[i]) : -1.0f;

        /* back-project to the display frame; collect extent */
        for (y = 0; y < h; y += 1)
            for (x = 0; x < w; x += 1) {
                double Z = depth[y * w + x], pc[3], pw[3];
                int k;
                if (Z <= 0.0)
                    continue;
                pc[0] = (Kinv[0] * x + Kinv[1] * y + Kinv[2]) * Z;
                pc[1] = (Kinv[4] * y + Kinv[5]) * Z;
                pc[2] = Z;
                for (k = 0; k < 3; k++)
                    pw[k] = r1.R[0 + k] * pc[0] + r1.R[3 + k] * pc[1]
                          + r1.R[6 + k] * pc[2] + C1[k];
                for (k = 0; k < 3; k++) {
                    if (pw[k] < lo[k]) lo[k] = pw[k];
                    if (pw[k] > hi[k]) hi[k] = pw[k];
                }
                nvalid++;
            }
        fprintf(stderr, "[densereal] %d valid depth px (%.1f%%)\n",
                nvalid, 100.0 * nvalid / (w * h));
        if (nvalid < 1000)
            return 1;

        if (mv_tsdf_init(&t, lo[0] - 0.02, lo[1] - 0.02, lo[2] - 0.02,
                         hi[0] - lo[0] + 0.04, hi[1] - lo[1] + 0.04,
                         hi[2] - lo[2] + 0.04, VOXEL, TAU) != MV_OK)
            return 1;

        for (y = 0; y < h; y++)
            for (x = 0; x < w; x++) {
                double Z = depth[y * w + x], pc[3], pw[3], uv[2], px[2];
                double sigZ, wgt;
                unsigned char rgb[3];
                int k;
                if (Z <= 0.0)
                    continue;
                pc[0] = (Kinv[0] * x + Kinv[1] * y + Kinv[2]) * Z;
                pc[1] = (Kinv[4] * y + Kinv[5]) * Z;
                pc[2] = Z;
                for (k = 0; k < 3; k++)
                    pw[k] = r1.R[0 + k] * pc[0] + r1.R[3 + k] * pc[1]
                          + r1.R[6 + k] * pc[2] + C1[k];
                sigZ = Z * Z * 0.35 / (fnew * Bnew);
                wgt = 1.0 / (sigZ * sigZ);
                mv_tsdf_fuse(&t, pw, C1, wgt);
                if ((x % 2) == 0 && (y % 2) == 0) {
                    rgb[0] = rgb[1] = rgb[2] = rectL[y * w + x];
                    mv_cloud_push(&cloud, pw, rgb);
                }
                /* display-flatness score: project into the ORIGINAL
                 * left camera; inside the decoded pattern bbox, the
                 * truth is the plane z = 0 of the display frame */
                if (mv_cam_project(uv, &A->cam, pw) == MV_OK) {
                    px[0] = uv[0];
                    px[1] = uv[1];
                    if (px[0] > A->uvmin[0] && px[0] < A->uvmax[0]
                        && px[1] > A->uvmin[1] && px[1] < A->uvmax[1]) {
                        zsum += fabs(pw[2]);
                        z2sum += pw[2] * pw[2];
                        nz++;
                    }
                }
            }

        {
            char path[512];
            double *tris;
            int ntri;
            snprintf(path, sizeof(path), "%s_disp.pgm", prefix);
            write_norm_pgm(path, disp, w, h, dmin, dmax);
            snprintf(path, sizeof(path), "%s_rect.pgm", prefix);
            mv_pgm_write(path, rectL, w, h);
            snprintf(path, sizeof(path), "%s_cloud.ply", prefix);
            mv_cloud_write_ply(path, &cloud);
            snprintf(path, sizeof(path), "%s_mesh.ply", prefix);
            mv_tsdf_write_ply(&t, path);
            if (mv_tsdf_mesh(&t, &tris, &ntri) == MV_OK) {
                printf("mesh triangles      : %d\n", ntri);
                free(tris);
            }
        }
        /* RANSAC plane on the reconstructed cloud (mv/plane.h): does
         * the dominant plane rediscover the display, unsupervised? */
        {
            double *pts = malloc((size_t)cloud.n * 3 * sizeof(double));
            if (pts) {
                double plane[4];
                unsigned char *inl = malloc((size_t)cloud.n);
                int k, ninl = 0;
                for (k = 0; k < cloud.n; k++) {
                    pts[3 * k] = cloud.xyz[3 * k];
                    pts[3 * k + 1] = cloud.xyz[3 * k + 1];
                    pts[3 * k + 2] = cloud.xyz[3 * k + 2];
                }
                if (inl && mv_plane_ransac(pts, cloud.n, 500, 0.004,
                                           20260812u, plane, inl)
                    == MV_OK) {
                    double ang = acos(fabs(plane[2])) * 180.0 / MV_PI;
                    for (k = 0; k < cloud.n; k++)
                        ninl += inl[k];
                    printf("RANSAC dominant plane: normal tilt vs "
                           "display %.2f deg, offset %.2f mm, "
                           "%d/%d inliers\n",
                           ang, fabs(plane[3]) * 1e3, ninl, cloud.n);
                }
                free(pts);
                free(inl);
            }
        }
        printf("frames decoded      : %d\n", nfr);
        printf("accidental baseline : %.1f mm (lateral %.1f mm)\n",
               B * 1e3, sqrt(bestlat) * 1e3);
        printf("valid depth pixels  : %d (%.1f%% of image)\n",
               nvalid, 100.0 * nvalid / (w * h));
        printf("cloud points        : %d\n", cloud.n);
        if (nz > 0)
            printf("display flatness    : %ld px on the display quad, "
                   "median-free RMS %.2f mm, mean |z| %.2f mm "
                   "(truth: the display IS z = 0)\n",
                   nz, sqrt(z2sum / nz) * 1e3, zsum / nz * 1e3);
        free(rawL);
        free(rawR);
        free(undL);
        free(undR);
        free(rectL);
        free(rectR);
        free(disp);
        free(depth);
        mv_tsdf_free(&t);
        mv_cloud_free(&cloud);
    }
    return 0;
}
