/* Temporal multiplexing of the two pattern tiers (mv/pattern.h):
 * schedule, per-tier round trips through synthetic cameras at each
 * tier's working range, cross-tier counter consistency on the shared
 * timeline, and clean wrong-tier failure. Deterministic (fixed seeds). */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

/* K_FINE / MV_PAT_MUX_BLOCK = 15626 (even -> fine block); K_COARSE is
 * the next block (odd -> coarse). Both fit the 20-bit fine counter. */
#define K_FINE 1000064u
#define K_COARSE (K_FINE + MV_PAT_MUX_BLOCK)

#define PITCH 0.0002745 /* display m/px, as in tests/test_mv.c */

static unsigned char pat[MV_PAT_W * MV_PAT_H];
static unsigned char img_near[640 * 480];
static unsigned char img_far[1280 * 720];

static unsigned dec_fine;   /* decoded by test_mux_near */
static unsigned dec_coarse; /* decoded by test_mux_far */

/* close-range camera of test_pattern in tests/test_mv.c */
static void near_cam(mv_camera *cam)
{
    double center[3] = { MV_PAT_W / 2.0 * PITCH,
                         MV_PAT_H / 2.0 * PITCH, 0.0 };
    int i;
    mv_cam_set_K(cam, 800.0, 800.0, 320.0, 240.0);
    mv_cam_set_identity_pose(cam);
    memset(cam->k, 0, sizeof(cam->k));
    for (i = 0; i < 3; i++)
        cam->t[i] = -center[i];
    cam->t[2] += 0.85;
}

/* room-range camera of test_pattern_coarse in tests/test_mv.c */
static void far_cam(mv_camera *cam)
{
    double center[3] = { MV_PAT_W / 2.0 * PITCH,
                         MV_PAT_H / 2.0 * PITCH, 0.0 };
    int i;
    mv_cam_set_K(cam, 1100.0, 1100.0, 640.0, 360.0);
    mv_cam_set_identity_pose(cam);
    memset(cam->k, 0, sizeof(cam->k));
    for (i = 0; i < 3; i++)
        cam->t[i] = -center[i];
    cam->t[2] += 3.6;
}

static void test_mux_schedule(void)
{
    unsigned k;
    int nfine = 0, ncoarse = 0;
    CHECK(mv_pattern_mux_tier(0) == 1, "mux: tier(0) is fine");
    CHECK(mv_pattern_mux_tier(63) == 1, "mux: tier(63) is fine");
    CHECK(mv_pattern_mux_tier(64) == 2, "mux: tier(64) is coarse");
    CHECK(mv_pattern_mux_tier(127) == 2, "mux: tier(127) is coarse");
    CHECK(mv_pattern_mux_tier(128) == 1, "mux: tier(128) is fine");
    for (k = 0; k < 1024; k++) {
        if (mv_pattern_mux_tier(k) == 1)
            nfine++;
        else
            ncoarse++;
    }
    CHECK(nfine == 512 && ncoarse == 512,
          "mux: each tier gets exactly half of k=0..1023");
}

static void test_mux_near(void)
{
    mv_camera cam;
    mv_read_result rr;
    unsigned long long seed = 1;
    CHECK(mv_pattern_mux_tier(K_FINE) == 1, "mux near: K_FINE is a fine block");
    near_cam(&cam);
    mv_pattern_mux_render(pat, K_FINE);
    CHECK(mv_render_plane(img_near, 640, 480, &cam, pat, MV_PAT_W,
                          MV_PAT_H, PITCH, 128, 0.0, &seed) == MV_OK,
          "mux near: render succeeds");
    CHECK(mv_read_pattern(&rr, img_near, 640, 480) == MV_OK,
          "mux near: fine blind read succeeds");
    CHECK(rr.n == 162, "mux near: all 162 corners identified");
    CHECK(rr.counter_valid && rr.counter == K_FINE,
          "mux near: counter == k exactly");
    dec_fine = rr.counter;
}

static void test_mux_far(void)
{
    mv_camera cam;
    mv_read_result rr;
    unsigned long long seed = 9;
    CHECK(mv_pattern_mux_tier(K_COARSE) == 2,
          "mux far: K_COARSE is a coarse block");
    far_cam(&cam);
    mv_pattern_mux_render(pat, K_COARSE);
    CHECK(mv_render_plane(img_far, 1280, 720, &cam, pat, MV_PAT_W,
                          MV_PAT_H, PITCH, 128, 0.0, &seed) == MV_OK,
          "mux far: render succeeds");
    CHECK(mv_read_coarse(&rr, img_far, 1280, 720) == MV_OK,
          "mux far: coarse blind read succeeds at room range");
    CHECK(rr.n == 10, "mux far: all 10 corners found");
    CHECK(rr.counter_valid && rr.counter == (K_COARSE & 255u),
          "mux far: counter == k mod 256");
    dec_coarse = rr.counter;
}

static void test_mux_consistency(void)
{
    /* Both decodes derive from the same display counter: the fine read
     * gives k in full, the coarse read k mod 256. The wrap-aware
     * mod-256 distance of tools/rigcalib.c (obs_dk: +384 keeps C's %
     * nonnegative, equivalent to ((d + 128) mod 256) - 128) must
     * recover the true block offset K_FINE - K_COARSE. */
    long dk = ((long)(dec_fine & 255u) - (long)(dec_coarse & 255u)
               + 384) % 256 - 128;
    CHECK(dk == -(long)MV_PAT_MUX_BLOCK,
          "mux: decoded counters land on one shared timeline");
}

static void test_mux_wrong_tier(void)
{
    /* Fine reader on a coarse-block frame: failing outright is the
     * expected clean outcome; a decode must never surface a valid
     * counter that aliases a fine counter. */
    mv_camera cam;
    mv_read_result rr;
    unsigned long long seed = 1;
    int r;
    near_cam(&cam);
    mv_pattern_mux_render(pat, K_COARSE);
    CHECK(mv_render_plane(img_near, 640, 480, &cam, pat, MV_PAT_W,
                          MV_PAT_H, PITCH, 128, 0.0, &seed) == MV_OK,
          "mux wrong-tier: render succeeds");
    rr.counter_valid = 1; /* poison: only a real decode may set it */
    r = mv_read_pattern(&rr, img_near, 640, 480);
    CHECK(r != MV_OK || rr.counter_valid == 0,
          "mux wrong-tier: fine reader yields no aliased counter");
}

int main(void)
{
    test_mux_schedule();
    test_mux_near();
    test_mux_far();
    test_mux_consistency();
    test_mux_wrong_tier();
    if (failures) {
        printf("\n%d FAILURE(S)\n", failures);
        return 1;
    }
    printf("\nall mux tests passed\n");
    return 0;
}
