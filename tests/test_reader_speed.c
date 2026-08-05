/* Full-HD reader throughput and decode-parity tests (speedup work
 * item).  The live hub decodes 1920x1080 frames on one thread for two
 * cameras, so per-frame cost is a shipped quantity: these tests pin it
 * as a ratio to a 640x480 decode (portable across build
 * configurations -- see decode_ref), and pin decode QUALITY at full HD --
 * corner ids, counter, and sub-0.1 px corner localization against
 * rendered ground truth -- so the fast paths (half-resolution decode
 * plus full-resolution corner polish) cannot silently trade precision
 * for speed.  Covers both tiers, strong tilt, near-180-degree roll
 * (historical bug regime), and partial pattern visibility.
 * Deterministic (fixed seeds). */

#define _POSIX_C_SOURCE 199309L
#define _DARWIN_C_SOURCE 1

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "mv/mv.h"

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

#define PITCH 0.0002745 /* display m/px, as in tests/test_mv.c */
#define IW 1920
#define IH 1080
#define K_FINE 654321u /* fits the 20-bit fine counter */
#define K_COARSE 173u

static unsigned char pat[MV_PAT_W * MV_PAT_H];
static unsigned char img[IW * IH];

static double now_s(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

/* Camera at `dist` metres from the pattern center, tilted about x by
 * tiltx and rolled about the axis by roll (radians), then shifted
 * laterally by dx metres (for partial-visibility views). */
static void pose_cam(mv_camera *cam, double dist, double f,
                     double tiltx, double roll, double dx)
{
    double center[3] = { MV_PAT_W / 2.0 * PITCH,
                         MV_PAT_H / 2.0 * PITCH, 0.0 };
    double ca = cos(tiltx), sa = sin(tiltx);
    double cr = cos(roll), sr = sin(roll);
    double Rz[9], Rx[9];
    int i;
    mv_cam_set_K(cam, f, f, IW / 2.0, IH / 2.0);
    mv_cam_set_identity_pose(cam);
    memset(cam->k, 0, sizeof(cam->k));
    Rz[0] = cr;  Rz[1] = -sr; Rz[2] = 0;
    Rz[3] = sr;  Rz[4] = cr;  Rz[5] = 0;
    Rz[6] = 0;   Rz[7] = 0;   Rz[8] = 1;
    Rx[0] = 1; Rx[1] = 0;  Rx[2] = 0;
    Rx[3] = 0; Rx[4] = ca; Rx[5] = -sa;
    Rx[6] = 0; Rx[7] = sa; Rx[8] = ca;
    mv_mat_mul(cam->R, Rz, Rx, 3, 3, 3);
    for (i = 0; i < 3; i++)
        cam->t[i] = -(cam->R[i * 3 + 0] * center[0]
                      + cam->R[i * 3 + 1] * center[1]
                      + cam->R[i * 3 + 2] * center[2]);
    cam->t[0] += dx;
    cam->t[2] += dist;
}

/* Ground-truth image position of a pattern corner: mv_render_plane
 * interpolates source samples as POINTS at integer continuum
 * coordinates, so the photometric cell boundary the reader measures
 * sits at display coordinate k - 0.5 (the existing 0.5 px-tolerance
 * tests absorb this; at 0.1 px it must be modeled). */
static void gt_corner(double uv[2], const mv_camera *cam,
                      const double xy[2])
{
    double X[3];
    X[0] = (xy[0] - 0.5) * PITCH;
    X[1] = (xy[1] - 0.5) * PITCH;
    X[2] = 0.0;
    mv_cam_project(uv, cam, X);
}

/* Score a decode against ground truth: every id must project within
 * 0.5 px (id correctness) and the localization RMS must be under
 * 0.1 px.  min_n additionally gates how much of the pattern must have
 * been identified. */
static void score_view(const char *name, const mv_read_result *rr,
                       const mv_camera *cam, int fine, int min_n)
{
    double se = 0.0, emax = 0.0;
    int i, id_bad = 0;
    char buf[128];

    for (i = 0; i < rr->n; i++) {
        double xy[2], uv[2], du, dv, e2;
        if (fine)
            mv_pattern_corner_px(rr->id[i] % MV_PAT_CORNER_COLS,
                                 rr->id[i] / MV_PAT_CORNER_COLS, xy);
        else
            mv_pattern2_corner_px(rr->id[i] % MV_PAT2_CORNER_COLS,
                                  rr->id[i] / MV_PAT2_CORNER_COLS, xy);
        gt_corner(uv, cam, xy);
        du = rr->uv[2 * i] - uv[0];
        dv = rr->uv[2 * i + 1] - uv[1];
        e2 = du * du + dv * dv;
        se += e2;
        if (e2 > emax)
            emax = e2;
        if (e2 > 0.25)
            id_bad++;
    }
    snprintf(buf, sizeof(buf), "%s: %d corners identified (>= %d)",
             name, rr->n, min_n);
    CHECK(rr->n >= min_n, buf);
    snprintf(buf, sizeof(buf), "%s: every corner id correct", name);
    CHECK(rr->n > 0 && id_bad == 0, buf);
    snprintf(buf, sizeof(buf),
             "%s: localization RMS %.3f px < 0.1 (max %.3f)", name,
             rr->n ? sqrt(se / rr->n) : 999.0, sqrt(emax));
    CHECK(rr->n > 0 && sqrt(se / rr->n) < 0.1, buf);
}

static void test_fine_fullhd_parity(void)
{
    /* dist, f, tiltx, roll, sigma: frontal-ish, strong tilt, and
     * >180-degree roll (the historical rotation-bug regime) */
    static const double V[3][5] = {
        { 0.70, 1400.0, 0.00, 0.08, 0.0 },
        { 0.60, 1400.0, 0.30, -0.12, 1.0 },
        { 0.80, 1600.0, -0.25, 3.24, 0.0 },
    };
    int v;

    mv_pattern_render(pat, K_FINE);
    for (v = 0; v < 3; v++) {
        mv_camera cam;
        mv_read_result rr;
        unsigned long long seed = 100 + v;
        char buf[128];
        pose_cam(&cam, V[v][0], V[v][1], V[v][2], V[v][3], 0.0);
        CHECK(mv_render_plane(img, IW, IH, &cam, pat, MV_PAT_W,
                              MV_PAT_H, PITCH, 128, V[v][4], &seed)
              == MV_OK, "fine: full-HD render succeeds");
        snprintf(buf, sizeof(buf), "fine view %d: blind read succeeds",
                 v);
        CHECK(mv_read_pattern(&rr, img, IW, IH) == MV_OK, buf);
        snprintf(buf, sizeof(buf), "fine view %d", v);
        score_view(buf, &rr, &cam, 1, 162);
        snprintf(buf, sizeof(buf),
                 "fine view %d: counter decodes exactly", v);
        CHECK(rr.counter_valid && rr.counter == K_FINE, buf);
    }
}

static void test_fine_partial_visibility(void)
{
    /* lateral shift pushes part of the lattice out of frame: the
     * reader must decode the visible subset with correct ids and
     * unchanged precision */
    mv_camera cam;
    mv_read_result rr;
    unsigned long long seed = 55;

    mv_pattern_render(pat, K_FINE);
    pose_cam(&cam, 0.70, 1400.0, 0.05, 0.10, 0.30);
    CHECK(mv_render_plane(img, IW, IH, &cam, pat, MV_PAT_W, MV_PAT_H,
                          PITCH, 128, 0.0, &seed) == MV_OK,
          "fine partial: render succeeds");
    CHECK(mv_read_pattern(&rr, img, IW, IH) == MV_OK,
          "fine partial: blind read succeeds");
    CHECK(rr.n < 162, "fine partial: pattern is genuinely clipped");
    score_view("fine partial", &rr, &cam, 1, 60);
}

static void test_coarse_fullhd_parity(void)
{
    /* room-range coarse views, incl. a >180-degree roll with noise */
    static const double V[3][5] = {
        { 2.00, 1400.0, 0.00, 0.10, 0.0 },
        { 1.80, 1400.0, 0.35, 3.30, 1.0 },
        { 2.40, 1600.0, -0.20, 0.50, 0.0 },
    };
    int v;

    mv_pattern2_render(pat, K_COARSE);
    for (v = 0; v < 3; v++) {
        mv_camera cam;
        mv_read_result rr;
        unsigned long long seed = 200 + v;
        char buf[128];
        pose_cam(&cam, V[v][0], V[v][1], V[v][2], V[v][3], 0.0);
        CHECK(mv_render_plane(img, IW, IH, &cam, pat, MV_PAT_W,
                              MV_PAT_H, PITCH, 128, V[v][4], &seed)
              == MV_OK, "coarse: full-HD render succeeds");
        snprintf(buf, sizeof(buf),
                 "coarse view %d: blind read succeeds", v);
        CHECK(mv_read_coarse(&rr, img, IW, IH) == MV_OK, buf);
        snprintf(buf, sizeof(buf), "coarse view %d", v);
        score_view(buf, &rr, &cam, 0, 10);
        snprintf(buf, sizeof(buf),
                 "coarse view %d: counter decodes exactly", v);
        CHECK(rr.counter_valid && rr.counter == (K_COARSE & 255u), buf);
    }
}

/* Machine-and-build speed reference: blind-decode one 640x480 fine
 * frame.  Absolute seconds are NOT portable here -- at -O2 the hub
 * worst case runs in 0.34 s, but the same code takes 13.2 s at
 * -O0 --coverage (40x) and more under ASan, so an absolute ceiling
 * passes on the dev machine and fails the coverage, sanitizer and
 * slow-runner CI jobs.  That is exactly how this test first went red:
 * a flaky gate that says nothing about the reader.
 *
 * The reference is deliberately the SAME reader on a SMALLER frame
 * (the sub-1.5 Mpx path, unchanged by the speedup work), not a
 * different kind of loop: instrumentation and optimization level scale
 * both sides together, so the ratio measures algorithmic cost rather
 * than the build.  A render-based reference was tried first and
 * rejected -- being a shallow per-pixel loop, it is far less sensitive
 * to instrumentation than the reader's nested search, and the ratio
 * still moved 5x between builds.
 *
 * What the ratio pins: full-HD decode must cost only a small multiple
 * of a 640x480 decode.  Before the downscale-first work it was ~40x;
 * after, ~3x. */
static unsigned char ref_img[640 * 480];

static double decode_ref(void)
{
    mv_camera cam;
    mv_read_result rr;
    unsigned long long seed = 3;
    double t0, best = 0.0;
    int i;
    mv_pattern_render(pat, K_FINE);
    mv_cam_set_K(&cam, 460.0, 460.0, 320.0, 240.0);
    mv_cam_set_identity_pose(&cam);
    memset(cam.k, 0, sizeof(cam.k));
    {   /* frame the pattern at 640x480, mildly tilted */
        double center[3] = { MV_PAT_W / 2.0 * PITCH,
                             MV_PAT_H / 2.0 * PITCH, 0.0 };
        int k;
        for (k = 0; k < 3; k++)
            cam.t[k] = -(cam.R[k * 3 + 0] * center[0]
                         + cam.R[k * 3 + 1] * center[1]
                         + cam.R[k * 3 + 2] * center[2]);
        cam.t[2] += 0.42;
    }
    mv_render_plane(ref_img, 640, 480, &cam, pat, MV_PAT_W, MV_PAT_H,
                    PITCH, 128, 0.0, &seed);
    for (i = 0; i < 3; i++) {
        double dt;
        t0 = now_s();
        mv_read_pattern(&rr, ref_img, 640, 480);
        dt = now_s() - t0;
        if (i == 0 || dt < best)
            best = dt;
    }
    return best > 1e-9 ? best : 1e-9;
}

static void test_speed_fullhd(void)
{
    /* Budgets in units of one 640x480 decode (see decode_ref).
     * Measured, and stable across build configurations because both
     * sides are the same code:
     *            fine   coarse   hub worst case
     *   -O2       2.1x     2.0x       73x
     *   -O0 cov   1.8x     1.8x       75x
     * The pre-speedup reader scored ~106x / ~106x / ~1160x (its
     * 640x480 path, and hence the reference, was already what it is
     * today), so these limits leave 4-8x headroom above the slowest
     * observed build while still catching any return to that regime
     * by an order of magnitude. */
    const double LIM_FINE = 8.0, LIM_COARSE = 12.0, LIM_HUB = 150.0;
    mv_camera cam;
    mv_read_result rr;
    unsigned long long seed = 7;
    double t0, t_fine, t_coarse, t_fallback, t_ref;
    char buf[200];
    int r;

    t_ref = decode_ref();
    printf("      (speed reference: %.3f s per 640x480 fine decode)\n",
           t_ref);

    mv_pattern_render(pat, K_FINE);
    pose_cam(&cam, 0.70, 1400.0, 0.10, 0.08, 0.0);
    mv_render_plane(img, IW, IH, &cam, pat, MV_PAT_W, MV_PAT_H, PITCH,
                    128, 1.0, &seed);
    t0 = now_s();
    r = mv_read_pattern(&rr, img, IW, IH);
    t_fine = now_s() - t0;
    CHECK(r == MV_OK, "speed: fine full-HD decode succeeds");
    snprintf(buf, sizeof(buf),
             "speed: fine full-HD decode %.3f s = %.1fx ref < %.0f",
             t_fine, t_fine / t_ref, LIM_FINE);
    CHECK(t_fine / t_ref < LIM_FINE, buf);

    mv_pattern2_render(pat, K_COARSE);
    pose_cam(&cam, 2.00, 1400.0, 0.10, 0.08, 0.0);
    seed = 11;
    mv_render_plane(img, IW, IH, &cam, pat, MV_PAT_W, MV_PAT_H, PITCH,
                    128, 1.0, &seed);
    t0 = now_s();
    r = mv_read_pattern(&rr, img, IW, IH); /* hub order: fine first */
    t_fallback = now_s() - t0;
    CHECK(r != MV_OK, "speed: fine reader cleanly rejects coarse frame");
    t0 = now_s();
    r = mv_read_coarse(&rr, img, IW, IH);
    t_coarse = now_s() - t0;
    CHECK(r == MV_OK, "speed: coarse full-HD decode succeeds");
    snprintf(buf, sizeof(buf),
             "speed: coarse full-HD decode %.3f s = %.1fx ref "
             "< %.0f", t_coarse, t_coarse / t_ref, LIM_COARSE);
    CHECK(t_coarse / t_ref < LIM_COARSE, buf);
    snprintf(buf, sizeof(buf),
             "speed: hub worst case (failed fine %.3f s + coarse) "
             "%.3f s = %.1fx ref < %.0f", t_fallback,
             t_fallback + t_coarse, (t_fallback + t_coarse) / t_ref,
             LIM_HUB);
    CHECK((t_fallback + t_coarse) / t_ref < LIM_HUB, buf);
}

int main(void)
{
    test_fine_fullhd_parity();
    test_fine_partial_visibility();
    test_coarse_fullhd_parity();
    test_speed_fullhd();
    if (failures) {
        printf("\n%d FAILURE(S)\n", failures);
        return 1;
    }
    printf("\nall reader-speed tests passed\n");
    return 0;
}
