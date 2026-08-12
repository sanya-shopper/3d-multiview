#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mv/mv.h"
#include "mv/track.h"

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
/* 1. track building: hand-built match graphs with known components   */

static void test_tracks_build(void)
{
    /* 3 images x 3 features.  Global ids: img0 = 0..2, img1 = 3..5,
     * img2 = 6..8.  Chains: {0,4,6} and {2,5,8}; 1, 3, 7 singletons. */
    static const int nfeat[3] = { 3, 3, 3 };
    static const int m01[3] = { 1, -1, 2 };  /* 0:0->1:1, 0:2->1:2 */
    static const int m12[3] = { -1, 0, 2 };  /* 1:1->2:0, 1:2->2:2 */
    static const int m02[3] = { -1, 0, -1 }; /* 0:1->2:0: the conflict */
    mv_match_pair pairs[3];
    mv_track *tr = NULL;
    int ntr = 0;

    pairs[0].a = 0; pairs[0].b = 1; pairs[0].idx2 = m01;
    pairs[1].a = 1; pairs[1].b = 2; pairs[1].idx2 = m12;

    CHECK(mv_tracks_build(&tr, &ntr, nfeat, 3, pairs, 2) == MV_OK
          && ntr == 2, "tracks: two chains -> two tracks");
    if (ntr == 2) {
        CHECK(tr[0].len == 3
              && tr[0].img[0] == 0 && tr[0].feat[0] == 0
              && tr[0].img[1] == 1 && tr[0].feat[1] == 1
              && tr[0].img[2] == 2 && tr[0].feat[2] == 0,
              "tracks: first chain observations exact");
        CHECK(tr[1].len == 3
              && tr[1].img[0] == 0 && tr[1].feat[0] == 2
              && tr[1].img[1] == 1 && tr[1].feat[1] == 2
              && tr[1].img[2] == 2 && tr[1].feat[2] == 2,
              "tracks: second chain observations exact");
    }
    mv_tracks_free(tr, ntr);

    /* adding 0:1->2:0 merges feature 1 of image 0 into the first
     * chain, which already holds feature 0 of image 0: the whole
     * component must be dropped, the clean chain must survive */
    pairs[2].a = 0; pairs[2].b = 2; pairs[2].idx2 = m02;
    tr = NULL;
    CHECK(mv_tracks_build(&tr, &ntr, nfeat, 3, pairs, 3) == MV_OK
          && ntr == 1 && tr[0].len == 3 && tr[0].feat[0] == 2,
          "tracks: conflicted component dropped whole");
    mv_tracks_free(tr, ntr);

    /* determinism: identical input -> identical output */
    {
        mv_track *ta = NULL, *tb = NULL;
        int na = 0, nb = 0, same = 1, i, k;
        mv_tracks_build(&ta, &na, nfeat, 3, pairs, 2);
        mv_tracks_build(&tb, &nb, nfeat, 3, pairs, 2);
        if (na != nb)
            same = 0;
        for (i = 0; same && i < na; i++) {
            if (ta[i].len != tb[i].len)
                same = 0;
            for (k = 0; same && k < ta[i].len; k++)
                if (ta[i].img[k] != tb[i].img[k]
                    || ta[i].feat[k] != tb[i].feat[k])
                    same = 0;
        }
        CHECK(same, "tracks: build is deterministic");
        mv_tracks_free(ta, na);
        mv_tracks_free(tb, nb);
    }

    /* corrupt input: match index out of range must be MV_ERR */
    {
        static const int bad[3] = { 5, -1, -1 };
        mv_match_pair pb;
        pb.a = 0; pb.b = 1; pb.idx2 = bad;
        tr = NULL;
        CHECK(mv_tracks_build(&tr, &ntr, nfeat, 3, &pb, 1) == MV_ERR,
              "tracks: out-of-range match rejected");
    }
}

/* ------------------------------------------------------------------ */
/* 2. incremental driver on exact projections (no images): the third
 * view must be recovered to numerical precision                      */

static void cloud_point(double X[3], int i)
{
    X[0] = -0.9 + 0.19 * (i % 10);
    X[1] = -0.7 + 0.21 * ((i / 2) % 7);
    X[2] = 3.2 + 0.23 * (i % 9);
}

static void test_driver_exact(void)
{
    enum { NP = 40 };
    static mv_feature fv[3][NP];
    static const mv_feature *feats[3];
    static int ident[NP];
    mv_camera gt[3], cams[3];
    mv_match_pair pairs[2];
    mv_track *tr = NULL;
    int nfeat[3] = { NP, NP, NP };
    unsigned char reg[3];
    unsigned char *xvalid;
    double *X, rms[3];
    int i, j, ntr = 0, nreg;

    mv_cam_set_K(&gt[0], 800.0, 800.0, 320.0, 240.0);
    mv_cam_set_identity_pose(&gt[0]);
    memset(gt[0].k, 0, sizeof(gt[0].k));
    gt[1] = gt[0];
    gt[2] = gt[0];
    {
        double C1[3] = { 0.4, 0.0, 0.0 }, C2[3] = { -0.35, 0.05, 0.25 };
        mv_cam_set_pose_yaw(&gt[1], -5.0 * MV_PI / 180.0, C1);
        mv_cam_set_pose_yaw(&gt[2], 6.0 * MV_PI / 180.0, C2);
    }
    for (i = 0; i < 3; i++) {
        feats[i] = fv[i];
        for (j = 0; j < NP; j++) {
            double P[3], uv[2];
            cloud_point(P, j);
            memset(&fv[i][j], 0, sizeof(mv_feature));
            mv_cam_project(uv, &gt[i], P);
            fv[i][j].u = uv[0];
            fv[i][j].v = uv[1];
            fv[i][j].scale = 1.0;
        }
    }
    /* identity matches through the chain 0-1-2: 40 three-view tracks */
    for (j = 0; j < NP; j++)
        ident[j] = j;
    pairs[0].a = 0; pairs[0].b = 1; pairs[0].idx2 = ident;
    pairs[1].a = 1; pairs[1].b = 2; pairs[1].idx2 = ident;
    CHECK(mv_tracks_build(&tr, &ntr, nfeat, 3, pairs, 2) == MV_OK
          && ntr == NP, "sfm-exact: 40 three-view tracks built");

    X = (double *)malloc((size_t)ntr * 3 * sizeof(double));
    xvalid = (unsigned char *)malloc((size_t)ntr);
    if (!X || !xvalid) {
        printf("FAIL: sfm-exact: alloc\n");
        failures++;
        free(X);
        free(xvalid);
        mv_tracks_free(tr, ntr);
        return;
    }
    cams[0] = gt[0];
    cams[1] = gt[1];             /* the calibrated seed pair */
    cams[2] = gt[0];             /* intrinsics prior, pose unknown */
    mv_cam_set_identity_pose(&cams[2]);
    nreg = mv_sfm_register_incremental(cams, reg, X, xvalid, rms,
                                       feats, nfeat, 3, tr, ntr, 0, 1);
    CHECK(nreg == 3 && reg[0] && reg[1] && reg[2],
          "sfm-exact: all three views registered");
    {
        double C[3], Cg[3], d, emax = 0.0;
        mv_cam_center(C, &cams[2]);
        mv_cam_center(Cg, &gt[2]);
        d = sqrt((C[0] - Cg[0]) * (C[0] - Cg[0])
                 + (C[1] - Cg[1]) * (C[1] - Cg[1])
                 + (C[2] - Cg[2]) * (C[2] - Cg[2]));
        for (i = 0; i < ntr; i++) {
            double P[3], e;
            cloud_point(P, i);
            e = fabs(X[3 * i] - P[0]) + fabs(X[3 * i + 1] - P[1])
              + fabs(X[3 * i + 2] - P[2]);
            if (!xvalid[i])
                e = 1.0;
            if (e > emax)
                emax = e;
        }
        CHECK(d < 1e-6, "sfm-exact: third view center exact (< 1e-6 m)");
        CHECK(emax < 1e-6, "sfm-exact: all track points exact");
        CHECK(rms[0] < 1e-6 && rms[1] < 1e-6 && rms[2] < 1e-6,
              "sfm-exact: reprojection RMS ~ 0");
    }
    free(X);
    free(xvalid);
    mv_tracks_free(tr, ntr);
}

/* ------------------------------------------------------------------ */
/* 3. end-to-end on rendered imagery: 6 views over a ~1.5x            *
 * magnification spread, multi-scale features, incremental            *
 * registration off the rig seed pair                                 */

#define W 640
#define H 480
#define NVIEW 6
#define NPAIR (NVIEW * (NVIEW - 1) / 2)
#define MAXF 400
#define NLEV 4

static unsigned long long rng = 12345;

static double urand(void)
{
    rng = rng * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)((rng >> 11) & ((1ULL << 53) - 1)) / 9007199254740992.0;
}

/* speckle texture with strong local contrast, one box blur (same
 * fixture idiom as tests/test_feat.c) */
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

/* ground-truth cameras: v0/v1 = the rig pair; v2..v5 walk toward the
 * scene, back wall at 4.2 m -> distances 3.5, 3.0, 2.7, 3.85 m, i.e.
 * a magnification spread of ~1.55x against the seed pair */
static void make_views(mv_camera gt[NVIEW])
{
    static const double Cs[4][3] = {
        { 0.2, -0.1, 0.7 }, { -0.3, 0.1, 1.2 },
        { 0.6, -0.1, 1.5 }, { -0.5, 0.0, 0.35 }
    };
    static const double yawdeg[4] = { -4.0, 5.0, -9.0, 7.0 };
    int i;
    mv_cam_set_K(&gt[0], 800.0, 800.0, 320.0, 240.0);
    mv_cam_set_identity_pose(&gt[0]);
    memset(gt[0].k, 0, sizeof(gt[0].k));
    for (i = 1; i < NVIEW; i++)
        gt[i] = gt[0];
    {
        double C2pos[3] = { 0.5, 0.0, 0.0 };
        mv_cam_set_pose_yaw(&gt[1], -6.0 * MV_PI / 180.0, C2pos);
    }
    for (i = 0; i < 4; i++)
        mv_cam_set_pose_yaw(&gt[2 + i], yawdeg[i] * MV_PI / 180.0,
                            Cs[i]);
}

static unsigned char imgs[NVIEW][W * H];
static mv_feature fv[NVIEW][MAXF];
static int midx[NPAIR][MAXF];

static void test_e2e(void)
{
    static const mv_feature *feats[NVIEW];
    mv_camera gt[NVIEW], cams[NVIEW], cams2[NVIEW];
    mv_match_pair pairs[NPAIR];
    mv_track *tr = NULL;
    int nfeat[NVIEW];
    unsigned char reg[NVIEW], reg2[NVIEW];
    unsigned char *xvalid, *xvalid2;
    double *X, *X2, rms[NVIEW], rms2[NVIEW];
    unsigned long long seed = 777;
    int i, j, np = 0, ntr = 0, nreg;

    make_scene();
    make_views(gt);
    for (i = 0; i < NVIEW; i++)
        mv_render_scene(imgs[i], NULL, W, H, &gt[i], pls, 3, 128, 1.0,
                        &seed);

    for (i = 0; i < NVIEW; i++) {
        feats[i] = fv[i];
        nfeat[i] = mv_feat_detect_ms(fv[i], MAXF, imgs[i], W, H, NLEV);
        if (nfeat[i] < 0)
            nfeat[i] = 0;
    }
    CHECK(nfeat[0] >= 200 && nfeat[NVIEW - 1] >= 200,
          "e2e: multi-scale features on every view");

    /* all 15 unordered pairs, house ratio 0.55 */
    for (i = 0; i < NVIEW; i++)
        for (j = i + 1; j < NVIEW; j++) {
            pairs[np].a = i;
            pairs[np].b = j;
            pairs[np].idx2 = midx[np];
            mv_feat_match(midx[np], fv[i], nfeat[i], fv[j], nfeat[j],
                          0.55);
            np++;
        }
    CHECK(mv_tracks_build(&tr, &ntr, nfeat, NVIEW, pairs, NPAIR) == MV_OK
          && ntr >= 100, "e2e: track set built (>= 100 tracks)");
    {
        int nlong = 0;
        for (i = 0; i < ntr; i++)
            if (tr[i].len >= 3)
                nlong++;
        printf("      (%d tracks, %d spanning >= 3 views)\n", ntr, nlong);
        CHECK(nlong >= 30, "e2e: >= 30 tracks span 3+ views");
    }

    X = (double *)malloc((size_t)ntr * 3 * sizeof(double));
    X2 = (double *)malloc((size_t)ntr * 3 * sizeof(double));
    xvalid = (unsigned char *)malloc((size_t)ntr);
    xvalid2 = (unsigned char *)malloc((size_t)ntr);
    if (!X || !X2 || !xvalid || !xvalid2) {
        printf("FAIL: e2e: alloc\n");
        failures++;
        free(X);
        free(X2);
        free(xvalid);
        free(xvalid2);
        mv_tracks_free(tr, ntr);
        return;
    }

    /* seed pair = the rig cameras with known poses; the rest carry
     * only the shared intrinsics prior */
    cams[0] = gt[0];
    cams[1] = gt[1];
    for (i = 2; i < NVIEW; i++) {
        cams[i] = gt[0];
        mv_cam_set_identity_pose(&cams[i]);
    }
    memcpy(cams2, cams, sizeof(cams2));

    nreg = mv_sfm_register_incremental(cams, reg, X, xvalid, rms,
                                       feats, nfeat, NVIEW, tr, ntr,
                                       0, 1);
    CHECK(nreg == NVIEW, "e2e: all 6 views registered");

    /* per-view pose error against ground truth: the single-scale
     * benchmark cost ~0.11 m of position error at this magnification
     * mismatch (tests/test_photo.c); multi-scale must land far below */
    {
        double dmax = 0.0, amax = 0.0, rmax = 0.0;
        for (i = 2; i < NVIEW; i++) {
            double C[3], Cg[3], d, Rt[9], D[9], trc, ang;
            if (!reg[i]) {
                dmax = 1e9;
                continue;
            }
            mv_cam_center(C, &cams[i]);
            mv_cam_center(Cg, &gt[i]);
            d = sqrt((C[0] - Cg[0]) * (C[0] - Cg[0])
                     + (C[1] - Cg[1]) * (C[1] - Cg[1])
                     + (C[2] - Cg[2]) * (C[2] - Cg[2]));
            mv_mat_transpose(Rt, gt[i].R, 3, 3);
            mv_mat_mul(D, cams[i].R, Rt, 3, 3, 3);
            trc = (D[0] + D[4] + D[8] - 1.0) / 2.0;
            if (trc > 1.0)
                trc = 1.0;
            if (trc < -1.0)
                trc = -1.0;
            ang = acos(trc) * 180.0 / MV_PI;
            printf("      view %d: pos err %.1f mm, rot err %.3f deg,"
                   " rms %.3f px\n", i, 1000.0 * d, ang, rms[i]);
            if (d > dmax)
                dmax = d;
            if (ang > amax)
                amax = ang;
            if (rms[i] > rmax)
                rmax = rms[i];
        }
        CHECK(dmax < 0.03,
              "e2e: every position error < 0.03 m (vs 0.11 m benchmark)");
        CHECK(amax < 0.5, "e2e: every orientation error < 0.5 deg");
        CHECK(rmax < 1.5, "e2e: every per-view reprojection RMS < 1.5 px");
    }

    /* the whole driver must be bit-deterministic */
    nreg = mv_sfm_register_incremental(cams2, reg2, X2, xvalid2, rms2,
                                       feats, nfeat, NVIEW, tr, ntr,
                                       0, 1);
    CHECK(nreg == NVIEW
          && memcmp(cams, cams2, sizeof(cams)) == 0
          && memcmp(reg, reg2, sizeof(reg)) == 0
          && memcmp(X, X2, (size_t)ntr * 3 * sizeof(double)) == 0
          && memcmp(xvalid, xvalid2, (size_t)ntr) == 0
          && memcmp(rms, rms2, sizeof(rms)) == 0,
          "e2e: incremental registration is deterministic");

    free(X);
    free(X2);
    free(xvalid);
    free(xvalid2);
    mv_tracks_free(tr, ntr);
}

int main(void)
{
    test_tracks_build();
    test_driver_exact();
    test_e2e();
    if (failures) {
        printf("\n%d FAILURE(S)\n", failures);
        return 1;
    }
    printf("\nall track tests passed\n");
    return 0;
}
