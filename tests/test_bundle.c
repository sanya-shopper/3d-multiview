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

/* ---- deterministic noise (project-standard LCG) ------------------ */

static unsigned long long rng;

static double urand(void)
{
    rng = rng * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)(rng >> 33) / 2147483648.0;
}

/* gaussian-ish: sum of 12 uniforms - 6 (unit variance) */
static double grand12(void)
{
    double g = -6.0;
    int i;
    for (i = 0; i < 12; i++)
        g += urand();
    return g;
}

/* ---- shared synthetic geometry ----------------------------------- */

/* Rig pair as in tests/test_photo.c plus one free "photo" camera with
 * a distinct K; the two rig cameras are the gauge anchor. */
enum { NPTS = 60, NCAM = 3, NOBS = NCAM * NPTS };

static const unsigned g_flags[NCAM] = {
    MV_BUNDLE_FIX_POSE | MV_BUNDLE_FIX_K,
    MV_BUNDLE_FIX_POSE | MV_BUNDLE_FIX_K,
    0
};

/* deterministic 3-D cloud in front of the cameras (test_photo's) */
static void cloud_point(double X[3], int i)
{
    X[0] = -0.9 + 0.19 * (i % 10);
    X[1] = -0.7 + 0.21 * ((i / 2) % 7);
    X[2] = 3.2 + 0.23 * (i % 9);
}

/* R = axis-angle rotation exp([r]x), th > 0 in all uses here */
static void aa_rot(double R[9], const double r[3])
{
    double th = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
    double kx = r[0] / th, ky = r[1] / th, kz = r[2] / th;
    double c = cos(th), s = sin(th), v = 1.0 - c;
    R[0] = c + kx * kx * v;
    R[1] = kx * ky * v - kz * s;
    R[2] = kx * kz * v + ky * s;
    R[3] = ky * kx * v + kz * s;
    R[4] = c + ky * ky * v;
    R[5] = ky * kz * v - kx * s;
    R[6] = kz * kx * v - ky * s;
    R[7] = kz * ky * v + kx * s;
    R[8] = c + kz * kz * v;
}

static void make_photo_cam(mv_camera *c)
{
    /* a "different phone": distinct focal, rolled, translated (as
     * tests/test_photo.c) */
    double Rz[9], Ry[9], ca, sa;
    mv_cam_set_K(c, 720.0, 728.0, 322.0, 242.0);
    mv_cam_set_identity_pose(c);
    memset(c->k, 0, sizeof(c->k));
    ca = cos(0.26);
    sa = sin(0.26);
    Rz[0] = ca; Rz[1] = -sa; Rz[2] = 0;
    Rz[3] = sa; Rz[4] = ca;  Rz[5] = 0;
    Rz[6] = 0;  Rz[7] = 0;   Rz[8] = 1;
    ca = cos(-0.15);
    sa = sin(-0.15);
    Ry[0] = ca;  Ry[1] = 0; Ry[2] = -sa;
    Ry[3] = 0;   Ry[4] = 1; Ry[5] = 0;
    Ry[6] = sa;  Ry[7] = 0; Ry[8] = ca;
    mv_mat_mul(c->R, Rz, Ry, 3, 3, 3);
    {
        double C[3] = { 0.9, -0.25, 0.3 };
        int j;
        for (j = 0; j < 3; j++)
            c->t[j] = -(c->R[j * 3 + 0] * C[0] + c->R[j * 3 + 1] * C[1]
                        + c->R[j * 3 + 2] * C[2]);
    }
}

static void make_cams(mv_camera cams[NCAM])
{
    double C2pos[3] = { 0.5, 0.0, 0.0 };
    mv_cam_set_K(&cams[0], 800.0, 800.0, 320.0, 240.0);
    mv_cam_set_identity_pose(&cams[0]);
    memset(cams[0].k, 0, sizeof(cams[0].k));
    cams[1] = cams[0];
    mv_cam_set_pose_yaw(&cams[1], -6.0 * 3.14159265358979324 / 180.0,
                        C2pos);
    make_photo_cam(&cams[2]);
}

/* one observation per camera-point pair, pixel noise sigma from rng */
static void make_obs(mv_bundle_obs *obs, const mv_camera cams[NCAM],
                     const double *X, double sigma)
{
    int c, i, k = 0;
    for (c = 0; c < NCAM; c++)
        for (i = 0; i < NPTS; i++) {
            double uv[2];
            mv_cam_project(uv, &cams[c], X + 3 * i);
            obs[k].cam = c;
            obs[k].pt = i;
            obs[k].u = uv[0] + sigma * grand12();
            obs[k].v = uv[1] + sigma * grand12();
            k++;
        }
}

/* perturbed initialization: free camera pose off by ~2 deg / ~5 cm,
 * focal 3%, principal point a few px; every point off by ~2 cm */
static void perturb_init(mv_camera *cam, double *X)
{
    static const double dr[3] = { 0.025, -0.020, 0.015 }; /* 2.0 deg */
    double dR[9], Rn[9];
    int i;
    aa_rot(dR, dr);
    mv_mat_mul(Rn, dR, cam->R, 3, 3, 3);
    memcpy(cam->R, Rn, sizeof(Rn));
    cam->t[0] += 0.03;
    cam->t[1] -= 0.04;
    cam->t[2] += 0.02;
    cam->K[0] *= 1.03;
    cam->K[4] *= 1.03;
    cam->K[2] += 3.0;
    cam->K[5] -= 2.0;
    for (i = 0; i < 3 * NPTS; i++)
        X[i] += 0.02 * (2.0 * urand() - 1.0);
}

static double scene_rms(const mv_camera *cams, const double *X,
                        const mv_bundle_obs *obs, int nobs)
{
    double ss = 0.0;
    int k;
    for (k = 0; k < nobs; k++) {
        double uv[2], du, dv;
        if (mv_cam_project(uv, &cams[obs[k].cam], X + 3 * obs[k].pt)
            != MV_OK)
            return 1e9;
        du = uv[0] - obs[k].u;
        dv = uv[1] - obs[k].v;
        ss += du * du + dv * dv;
    }
    return sqrt(ss / nobs);
}

/* ---- 1. sigma = 0 exactness -------------------------------------- */

static void test_exact(void)
{
    mv_camera gt[NCAM], cams[NCAM];
    static double Xgt[3 * NPTS], Xw[3 * NPTS];
    static mv_bundle_obs obs[NOBS];
    double rms = 1e9;
    int i;

    rng = 1ULL;
    make_cams(gt);
    for (i = 0; i < NPTS; i++)
        cloud_point(Xgt + 3 * i, i);
    make_obs(obs, gt, Xgt, 0.0);
    memcpy(cams, gt, sizeof(cams));
    memcpy(Xw, Xgt, sizeof(Xw));
    perturb_init(&cams[2], Xw);

    CHECK(mv_bundle_adjust(cams, g_flags, NCAM, Xw, NPTS, obs, NOBS,
                           &rms) == MV_OK,
          "exact: bundle adjustment succeeds");
    CHECK(rms < 1e-8, "exact: reprojection rms < 1e-8");
    {
        double kerr = fabs(cams[2].K[0] - 720.0)
                    + fabs(cams[2].K[4] - 728.0)
                    + fabs(cams[2].K[2] - 322.0)
                    + fabs(cams[2].K[5] - 242.0);
        CHECK(kerr < 1e-4, "exact: free camera K within 1e-4");
    }
    {
        double perr = 0.0;
        for (i = 0; i < 9; i++)
            perr += fabs(cams[2].R[i] - gt[2].R[i]);
        for (i = 0; i < 3; i++)
            perr += fabs(cams[2].t[i] - gt[2].t[i]);
        CHECK(perr < 1e-6, "exact: free camera pose within 1e-6");
    }
    {
        double pmax = 0.0;
        for (i = 0; i < 3 * NPTS; i++)
            if (fabs(Xw[i] - Xgt[i]) > pmax)
                pmax = fabs(Xw[i] - Xgt[i]);
        CHECK(pmax < 1e-6, "exact: points within 1e-6");
    }
    {
        mv_camera fixed = gt[0];
        double ferr = 0.0;
        for (i = 0; i < 9; i++)
            ferr += fabs(cams[0].R[i] - fixed.R[i])
                  + fabs(cams[0].K[i] - fixed.K[i]);
        CHECK(ferr == 0.0, "exact: fixed camera untouched");
    }
}

/* ---- 2. noise + statistical acceptance (T5 discipline) ----------- */

static void test_t5(void)
{
    /* Pre-registered acceptance: after joint ML the ensemble mean of
     * RMS/predicted must sit in [0.90, 1.06], predicted per-point
     * distance RMS = sqrt(2) sigma sqrt(1 - p/n) with p the FREE
     * scalar parameters (10 photo-camera + 3*NPTS) and n = 2*NOBS. */
    enum { NSEEDS = 10 };
    const double sigma = 0.4;
    const int pfit = 10 + 3 * NPTS;
    const int nres = 2 * NOBS;
    const double predicted =
        sqrt(2.0) * sigma * sqrt(1.0 - (double)pfit / (double)nres);
    mv_camera gt[NCAM], cams[NCAM];
    static double Xgt[3 * NPTS], Xw[3 * NPTS];
    static mv_bundle_obs obs[NOBS];
    double sum_ratio = 0.0;
    int s, i, ok_all = 1, always_reduces = 1;

    for (s = 0; s < NSEEDS; s++) {
        double rms_init, rms_ba = 1e9;
        rng = 20260805ULL + 1000ULL * (unsigned)s;
        make_cams(gt);
        for (i = 0; i < NPTS; i++)
            cloud_point(Xgt + 3 * i, i);
        make_obs(obs, gt, Xgt, sigma);
        memcpy(cams, gt, sizeof(cams));
        memcpy(Xw, Xgt, sizeof(Xw));
        perturb_init(&cams[2], Xw);
        rms_init = scene_rms(cams, Xw, obs, NOBS);
        if (mv_bundle_adjust(cams, g_flags, NCAM, Xw, NPTS, obs, NOBS,
                             &rms_ba) != MV_OK) {
            ok_all = 0;
            break;
        }
        if (rms_ba >= rms_init)
            always_reduces = 0;
        sum_ratio += rms_ba / predicted;
    }
    CHECK(ok_all, "t5: all ensemble members adjust");
    if (ok_all) {
        double mean_ratio = sum_ratio / NSEEDS;
        printf("t5: mean ratio %.4f (predicted rms %.4f px)\n",
               mean_ratio, predicted);
        CHECK(mean_ratio >= 0.90 && mean_ratio <= 1.06,
              "t5: rms/predicted ratio within [0.90, 1.06]");
        CHECK(always_reduces,
              "t5: cost strictly reduced on every seed");
    }
}

/* ---- 3. loose gauge is rejected ---------------------------------- */

static void test_gauge(void)
{
    static const unsigned loose[NCAM] = { 0, 0, 0 };
    static const unsigned konly[NCAM] = {
        MV_BUNDLE_FIX_K, MV_BUNDLE_FIX_K, 0
    };
    mv_camera cams[NCAM];
    static double Xw[3 * NPTS];
    static mv_bundle_obs obs[NOBS];
    int i;

    make_cams(cams);
    for (i = 0; i < NPTS; i++)
        cloud_point(Xw + 3 * i, i);
    make_obs(obs, cams, Xw, 0.0);
    CHECK(mv_bundle_adjust(cams, loose, NCAM, Xw, NPTS, obs, NOBS, NULL)
              == MV_ERR,
          "gauge: no fixed pose rejected");
    CHECK(mv_bundle_adjust(cams, konly, NCAM, Xw, NPTS, obs, NOBS, NULL)
              == MV_ERR,
          "gauge: fixed K alone is not a gauge anchor");
    CHECK(mv_bundle_adjust(cams, NULL, NCAM, Xw, NPTS, obs, NOBS, NULL)
              == MV_ERR,
          "gauge: NULL camflags rejected");
}

/* ---- 4. photo-registration integration --------------------------- */

/* The pipeline of tests/test_photo.c: register the photo, rebuild the
 * anchor/match observations exactly as src/photo.c does internally,
 * bundle-adjust with the rig fixed, and demand the photo camera's
 * position error improve.
 *
 * Plane placement deviates from test_photo deliberately.  Measured on
 * that scene, the rig pair's anchor bank is 100% coplanar (all 314
 * anchors on the fronto-parallel back plane; the grazing floor/side
 * planes yield no rig matches even at a 2000-feature budget), and a
 * free-K resection of coplanar points is structurally under-determined
 * (10 parameters against the 8 observable homography DOF -- the same
 * degeneracy src/photo.c documents by pinning the principal point), so
 * ML bundle adjustment correctly walks the flat focal/depth family to
 * a spurious solution metres away at lower residual.  Here the back
 * wall is folded into a shallow V (two walls tilted +-0.35 rad about
 * vertical), which gives the anchor bank ~0.8 m of depth relief with
 * opposite gradients: the free camera becomes observable and BA can
 * tighten the registration, principal point included. */

static unsigned char texA[160 * 80], texB[140 * 160], texC[160 * 80];
static mv_plane pls[3];

static void make_tex(unsigned char *t, int w, int h)
{
    int i;
    for (i = 0; i < w * h; i++)
        t[i] = (unsigned char)(30 + 200 * urand());
}

static void set_plane(mv_plane *p, const unsigned char *tex, int tw,
                      int th, double pitch, const double r1[3],
                      const double r2[3], const double t[3])
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

#define BA_MAXF 500

static void test_photo_ba(void)
{
    static unsigned char img1[640 * 480], img2[640 * 480],
        imgp[640 * 480];
    static mv_feature f1[BA_MAXF], f2[BA_MAXF], fp[BA_MAXF],
        fa[BA_MAXF];
    static double ax[3 * BA_MAXF], uv1a[2 * BA_MAXF], uv2a[2 * BA_MAXF];
    static double Xba[3 * BA_MAXF];
    static mv_bundle_obs obs[3 * BA_MAXF];
    static int idx[BA_MAXF];
    mv_camera c1, c2, gt, est, bc[3];
    unsigned long long seed = 3;
    double pre, post, rms = 1e9;
    int n1, n2, np, na = 0, nb = 0, nobs = 0, i;

    /* seed 7: a draw whose registration error (~77 mm) is at the
     * documented scale of the pipeline (doc/multiview.tex: 0.11 m) */
    rng = 7ULL;
    {
        static const double ex[3] = { 1, 0, 0 }, ey[3] = { 0, 1, 0 },
                            ez[3] = { 0, 0, 1 };
        const double ca = cos(0.35), sa = sin(0.35);
        double t0[3], r1t[3];
        make_tex(texA, 160, 80);
        make_tex(texB, 140, 160);
        make_tex(texC, 160, 80);
        r1t[0] = ca; r1t[1] = 0.0; r1t[2] = sa; /* left V wall */
        t0[0] = -1.9; t0[1] = -0.8; t0[2] = 4.6 - 1.6 * sa;
        set_plane(&pls[0], texA, 160, 80, 0.02, r1t, ey, t0);
        r1t[0] = ca; r1t[1] = 0.0; r1t[2] = -sa; /* right V wall */
        t0[0] = -2.5 + 3.2 * ca; t0[1] = -0.8; t0[2] = 4.4 + 1.6 * sa;
        set_plane(&pls[2], texC, 160, 80, 0.02, r1t, ey, t0);
        t0[0] = -1.4; t0[1] = 0.8; t0[2] = 2.6;
        set_plane(&pls[1], texB, 140, 160, 0.02, ex, ez, t0); /* floor */
    }
    mv_cam_set_K(&c1, 800.0, 800.0, 320.0, 240.0);
    mv_cam_set_identity_pose(&c1);
    memset(c1.k, 0, sizeof(c1.k));
    c2 = c1;
    {
        double C2pos[3] = { 0.5, 0.0, 0.0 };
        mv_cam_set_pose_yaw(&c2, -6.0 * 3.14159265358979324 / 180.0,
                            C2pos);
    }
    make_photo_cam(&gt);

    mv_render_scene(img1, NULL, 640, 480, &c1, pls, 3, 10, 1.0, &seed);
    mv_render_scene(img2, NULL, 640, 480, &c2, pls, 3, 10, 1.0, &seed);
    mv_render_scene(imgp, NULL, 640, 480, &gt, pls, 3, 10, 1.0, &seed);

    CHECK(mv_photo_register(&est, &c1, &c2, img1, img2, imgp, 640, 480,
                            640, 480, NULL, NULL) == MV_OK,
          "photo: registration succeeds");
    {
        double C[3], Cg[3];
        mv_cam_center(C, &est);
        mv_cam_center(Cg, &gt);
        pre = sqrt((C[0] - Cg[0]) * (C[0] - Cg[0])
                   + (C[1] - Cg[1]) * (C[1] - Cg[1])
                   + (C[2] - Cg[2]) * (C[2] - Cg[2]));
    }

    /* anchor bank exactly as src/photo.c: rig match, triangulate, gate
     * on positive depth in both views and a 2 px reprojection bound;
     * keep camera 1's descriptor plus BOTH rig pixels */
    n1 = mv_feat_detect(f1, BA_MAXF, img1, 640, 480);
    n2 = mv_feat_detect(f2, BA_MAXF, img2, 640, 480);
    np = mv_feat_detect(fp, BA_MAXF, imgp, 640, 480);
    CHECK(n1 >= 20 && n2 >= 20 && np >= 20, "photo: features detected");
    CHECK(mv_feat_match(idx, f1, n1, f2, n2, 0.6) >= 8,
          "photo: rig pair matches");
    for (i = 0; i < n1 && na < BA_MAXF; i++) {
        const mv_camera *cams2[2];
        double uvs[4], X[3], p1[2], p2[2], e1, e2, z1, z2;
        if (idx[i] < 0)
            continue;
        cams2[0] = &c1;
        cams2[1] = &c2;
        uvs[0] = f1[i].u;
        uvs[1] = f1[i].v;
        uvs[2] = f2[idx[i]].u;
        uvs[3] = f2[idx[i]].v;
        if (mv_triangulate(X, cams2, uvs, 2) != MV_OK)
            continue;
        z1 = c1.R[6] * X[0] + c1.R[7] * X[1] + c1.R[8] * X[2] + c1.t[2];
        z2 = c2.R[6] * X[0] + c2.R[7] * X[1] + c2.R[8] * X[2] + c2.t[2];
        if (z1 <= 0.0 || z2 <= 0.0)
            continue;
        if (mv_cam_project(p1, &c1, X) != MV_OK
            || mv_cam_project(p2, &c2, X) != MV_OK)
            continue;
        e1 = fabs(p1[0] - uvs[0]) + fabs(p1[1] - uvs[1]);
        e2 = fabs(p2[0] - uvs[2]) + fabs(p2[1] - uvs[3]);
        if (e1 > 2.0 || e2 > 2.0)
            continue;
        fa[na] = f1[i];
        ax[3 * na] = X[0];
        ax[3 * na + 1] = X[1];
        ax[3 * na + 2] = X[2];
        uv1a[2 * na] = uvs[0];
        uv1a[2 * na + 1] = uvs[1];
        uv2a[2 * na] = uvs[2];
        uv2a[2 * na + 1] = uvs[3];
        na++;
    }
    CHECK(na >= 30, "photo: anchor bank populated");

    /* photo -> anchor matches; gate at the pipeline's 2.5 px RANSAC
     * inlier radius against the registered camera (the equivalent of
     * mv_resect_robust's inlier selection, which photo.c discards) */
    CHECK(mv_feat_match(idx, fp, np, fa, na, 0.6) >= 6,
          "photo: photo-to-anchor matches");
    for (i = 0; i < np; i++) {
        double uv[2], du, dv;
        int j = idx[i];
        if (j < 0)
            continue;
        if (mv_cam_project(uv, &est, ax + 3 * j) != MV_OK)
            continue;
        du = uv[0] - fp[i].u;
        dv = uv[1] - fp[i].v;
        if (du * du + dv * dv >= 6.25)
            continue;
        memcpy(Xba + 3 * nb, ax + 3 * j, 3 * sizeof(double));
        obs[nobs].cam = 0;
        obs[nobs].pt = nb;
        obs[nobs].u = uv1a[2 * j];
        obs[nobs].v = uv1a[2 * j + 1];
        nobs++;
        obs[nobs].cam = 1;
        obs[nobs].pt = nb;
        obs[nobs].u = uv2a[2 * j];
        obs[nobs].v = uv2a[2 * j + 1];
        nobs++;
        obs[nobs].cam = 2;
        obs[nobs].pt = nb;
        obs[nobs].u = fp[i].u;
        obs[nobs].v = fp[i].v;
        nobs++;
        nb++;
    }
    printf("photo: joint set %d points, %d observations\n", nb, nobs);
    CHECK(nb >= 15, "photo: joint point set populated");

    bc[0] = c1;
    bc[1] = c2;
    bc[2] = est;
    CHECK(mv_bundle_adjust(bc, g_flags, 3, Xba, nb, obs, nobs, &rms)
              == MV_OK,
          "photo: bundle adjustment succeeds");
    {
        double C[3], Cg[3];
        mv_cam_center(C, &bc[2]);
        mv_cam_center(Cg, &gt);
        post = sqrt((C[0] - Cg[0]) * (C[0] - Cg[0])
                    + (C[1] - Cg[1]) * (C[1] - Cg[1])
                    + (C[2] - Cg[2]) * (C[2] - Cg[2]));
    }
    printf("photo: position error pre-BA %.1f mm, post-BA %.1f mm "
           "(rms %.3f px)\n", 1000.0 * pre, 1000.0 * post, rms);
    CHECK(post <= pre * 0.8, "photo: BA improves position error >= 20%");
}

int main(void)
{
    test_exact();
    test_t5();
    test_gauge();
    test_photo_ba();
    if (failures) {
        printf("\n%d FAILURE(S)\n", failures);
        return 1;
    }
    printf("\nall bundle tests passed\n");
    return 0;
}
