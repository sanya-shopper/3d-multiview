#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mv/mv.h"

#define MV_PI 3.14159265358979323846

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

/* ------------------------------------------------------------------ */
/* deterministic synthetic scene: three speckle-textured planes seen
 * by a stereo-ish pair (same idiom as tests/test_mv.c make_pair and
 * demo/roomsim.c make_tex) */

#define W 640
#define H 480
#define MAXF 500

static unsigned long long rng = 12345;

static double urand(void)
{
    rng = rng * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)((rng >> 11) & ((1ULL << 53) - 1)) / 9007199254740992.0;
}

/* speckle texture with strong local contrast, one box blur */
static void make_tex(unsigned char *t, int w, int h)
{
    static unsigned char tmp[160 * 160];
    int i, x, y;
    for (i = 0; i < w * h; i++)
        t[i] = (unsigned char)(40 + 175 * (urand() > 0.5));
    memcpy(tmp, t, (size_t)w * h);
    for (y = 1; y < h - 1; y++)
        for (x = 1; x < w - 1; x++) {
            int s = 0, a, b;
            for (b = -1; b <= 1; b++)
                for (a = -1; a <= 1; a++)
                    s += tmp[(y + b) * w + (x + a)];
            t[y * w + x] = (unsigned char)(s / 9);
        }
}

static void set_plane(mv_plane *p, const unsigned char *tex, int tw, int th,
                      double pitch, const double r1[3], const double r2[3],
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
    p->tw = tw;
    p->th = th;
    p->pitch = pitch;
}

static unsigned char texA[160 * 80];  /* back wall */
static unsigned char texB[140 * 160]; /* floor */
static unsigned char texC[160 * 80];  /* side wall */
static mv_plane pls[3];

static void make_scene(void)
{
    static const double ex[3] = { 1, 0, 0 }, ey[3] = { 0, 1, 0 },
                        ez[3] = { 0, 0, 1 };
    double t0[3];
    rng = 12345;
    make_tex(texA, 160, 80);
    make_tex(texB, 140, 160);
    make_tex(texC, 160, 80);
    t0[0] = -1.6; t0[1] = -0.8; t0[2] = 4.2;
    set_plane(&pls[0], texA, 160, 80, 0.02, ex, ey, t0);  /* back wall */
    t0[0] = -1.4; t0[1] = 0.9; t0[2] = 2.2;
    set_plane(&pls[1], texB, 140, 160, 0.02, ex, ez, t0); /* floor */
    t0[0] = 1.7; t0[1] = -0.8; t0[2] = 2.2;
    set_plane(&pls[2], texC, 160, 80, 0.02, ez, ey, t0);  /* side wall */
}

static void make_pair(mv_camera *c1, mv_camera *c2)
{
    double C2pos[3] = { 0.5, 0.0, 0.0 };
    mv_cam_set_K(c1, 800.0, 800.0, 320.0, 240.0);
    mv_cam_set_identity_pose(c1);
    memset(c1->k, 0, sizeof(c1->k));
    *c2 = *c1;
    mv_cam_set_pose_yaw(c2, -6.0 * MV_PI / 180.0, C2pos);
}

/* compose a roll about the camera z axis: Xc' = Rz Xc, camera center
 * unchanged */
static void roll_camera(mv_camera *c, double angle)
{
    double ca = cos(angle), sa = sin(angle);
    double Rz[9], Rn[9], tn[3];
    int i;
    Rz[0] = ca; Rz[1] = -sa; Rz[2] = 0.0;
    Rz[3] = sa; Rz[4] = ca;  Rz[5] = 0.0;
    Rz[6] = 0.0; Rz[7] = 0.0; Rz[8] = 1.0;
    mv_mat_mul(Rn, Rz, c->R, 3, 3, 3);
    mv_mat_mul(tn, Rz, c->t, 3, 3, 1);
    for (i = 0; i < 9; i++)
        c->R[i] = Rn[i];
    for (i = 0; i < 3; i++)
        c->t[i] = tn[i];
}

static int cmp_double(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

static double median(double *v, int n)
{
    qsort(v, (size_t)n, sizeof(double), cmp_double);
    return (n % 2) ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

/* detect + match both views, collect matched pixel pairs; returns the
 * match count */
static int match_views(double *uv1, double *uv2,
                       const unsigned char *im1, const unsigned char *im2,
                       int *n1out, int *n2out)
{
    static mv_feature f1[MAXF], f2[MAXF];
    static int idx2[MAXF];
    int n1, n2, nm = 0, i;
    n1 = mv_feat_detect(f1, MAXF, im1, W, H);
    n2 = mv_feat_detect(f2, MAXF, im2, W, H);
    if (n1out) *n1out = n1;
    if (n2out) *n2out = n2;
    if (n1 <= 0 || n2 <= 0)
        return 0;
    /* 0.55 SSD ratio ~ 0.74 in L2, near Lowe's working point */
    if (mv_feat_match(idx2, f1, n1, f2, n2, 0.55) < 0)
        return 0;
    for (i = 0; i < n1; i++)
        if (idx2[i] >= 0) {
            uv1[2 * nm + 0] = f1[i].u;
            uv1[2 * nm + 1] = f1[i].v;
            uv2[2 * nm + 0] = f2[idx2[i]].u;
            uv2[2 * nm + 1] = f2[idx2[i]].v;
            nm++;
        }
    return nm;
}

/* robust F fit on the matches (compacts uv arrays in place), then the
 * median symmetric epipolar distance of the survivors under the ANALYTIC
 * F of the camera pair; returns survivor count, median via *med */
static int verify(double *med, double *uv1, double *uv2, int nm,
                  const mv_camera *c1, const mv_camera *c2)
{
    double F[9], Fa[9];
    static double d[MAXF];
    int n = nm, i;
    if (nm < 8)
        return 0;
    if (mv_fundamental_8point_robust(F, uv1, uv2, &n) != MV_OK)
        return 0;
    mv_fundamental_from_cams(Fa, c1, c2);
    for (i = 0; i < n; i++)
        d[i] = mv_sym_epipolar_dist(Fa, uv1 + 2 * i, uv2 + 2 * i);
    *med = median(d, n);
    return n;
}

static unsigned char im1[W * H], im2[W * H];
static double uv1[2 * MAXF], uv2[2 * MAXF];

int main(void)
{
    mv_camera c1, c2;
    unsigned long long seed;
    int n1, n2, nm, nsurv;
    double med;

    make_scene();
    make_pair(&c1, &c2);
    seed = 777;
    mv_render_scene(im1, NULL, W, H, &c1, pls, 3, 128, 1.0, &seed);
    mv_render_scene(im2, NULL, W, H, &c2, pls, 3, 128, 1.0, &seed);

    /* 1. detection determinism and count */
    {
        static mv_feature fa[MAXF], fb[MAXF];
        int na = mv_feat_detect(fa, MAXF, im1, W, H);
        int nb = mv_feat_detect(fb, MAXF, im1, W, H);
        CHECK(na == nb && na > 0
              && memcmp(fa, fb, (size_t)na * sizeof(mv_feature)) == 0,
              "feat: detection is deterministic");
        printf("      (view1: %d features)\n", na);
    }

    nm = match_views(uv1, uv2, im1, im2, &n1, &n2);
    printf("      (view1 %d, view2 %d features, %d matches)\n", n1, n2, nm);
    CHECK(n1 >= 100 && n2 >= 100, "feat: >= 100 features per 640x480 view");

    /* 2. mutual ratio-test matches */
    CHECK(nm >= 40, "feat: >= 40 mutual ratio-test matches");

    /* 3. geometric verification against the analytic F */
    nsurv = verify(&med, uv1, uv2, nm, &c1, &c2);
    printf("      (robust F: %d survivors, median analytic residual"
           " %.3f px)\n", nsurv, med);
    CHECK(nsurv >= 30, "feat: >= 30 robust-F survivors");
    CHECK(med < 1.0, "feat: median analytic epipolar residual < 1 px");

    /* 4. rotation robustness: 20-degree camera roll on view 2 */
    {
        mv_camera c2r = c2;
        roll_camera(&c2r, 20.0 * MV_PI / 180.0);
        seed = 777;
        mv_render_scene(im2, NULL, W, H, &c2r, pls, 3, 128, 1.0, &seed);
        nm = match_views(uv1, uv2, im1, im2, &n1, &n2);
        nsurv = verify(&med, uv1, uv2, nm, &c1, &c2r);
        printf("      (roll 20deg: %d matches, %d survivors, median"
               " %.3f px)\n", nm, nsurv, med);
        CHECK(nsurv >= 25, "feat: >= 25 verified matches under 20deg roll");
    }

    /* 5. noise-free localization: sub-pixel detection accuracy */
    {
        seed = 777;
        mv_render_scene(im1, NULL, W, H, &c1, pls, 3, 128, 0.0, &seed);
        mv_render_scene(im2, NULL, W, H, &c2, pls, 3, 128, 0.0, &seed);
        nm = match_views(uv1, uv2, im1, im2, &n1, &n2);
        nsurv = verify(&med, uv1, uv2, nm, &c1, &c2);
        printf("      (sigma=0: %d matches, %d survivors, median"
               " %.3f px)\n", nm, nsurv, med);
        CHECK(nsurv >= 30 && med < 0.5,
              "feat: sigma=0 median analytic residual < 0.5 px");
    }

    if (failures) {
        printf("\n%d FAILURE(S)\n", failures);
        return 1;
    }
    printf("\nall feat tests passed\n");
    return 0;
}
