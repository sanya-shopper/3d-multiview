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

static unsigned long long rng = 20260806ULL;
static double urand(void)
{
    rng = rng * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)(rng >> 33) / 2147483648.0;
}
static double grand(void)
{
    double s = -6.0;
    int i;
    for (i = 0; i < 12; i++)
        s += urand();
    return s;
}

/* deterministic 3-D cloud in front of the cameras */
static void cloud_point(double X[3], int i)
{
    X[0] = -0.9 + 0.19 * (i % 10);
    X[1] = -0.7 + 0.21 * ((i / 2) % 7);
    X[2] = 3.2 + 0.23 * (i % 9);
}

static void make_photo_cam(mv_camera *c)
{
    /* a "different phone": distinct focal (~10% from the rig's 800,
     * inside the single-scale envelope mv/feat.h documents), principal
     * point near center (true of real phone cameras -- and the
     * registration pipeline pins pp at center, the standard casual-photo
     * convention), rolled and translated viewpoint */
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

static void test_resection(void)
{
    mv_camera gt, est;
    double P[12];
    static double X3[3 * 40], uv[2 * 40];
    int i, n = 40;

    make_photo_cam(&gt);
    for (i = 0; i < n; i++) {
        cloud_point(X3 + 3 * i, i);
        mv_cam_project(uv + 2 * i, &gt, X3 + 3 * i);
    }
    CHECK(mv_resect_dlt(P, X3, uv, n) == MV_OK, "resect: DLT succeeds");
    CHECK(mv_cam_from_projection(&est, P) == MV_OK,
          "resect: RQ decomposition succeeds");
    {
        double e = fabs(est.K[0] - 720.0) + fabs(est.K[4] - 728.0)
                 + fabs(est.K[2] - 322.0) + fabs(est.K[5] - 242.0)
                 + fabs(est.K[1]);
        CHECK(e < 1e-6, "resect: K exact on noiseless data");
    }
    {
        double e = 0.0;
        for (i = 0; i < 9; i++)
            e += fabs(est.R[i] - gt.R[i]);
        for (i = 0; i < 3; i++)
            e += fabs(est.t[i] - gt.t[i]);
        CHECK(e < 1e-6, "resect: pose exact on noiseless data");
    }

    /* robust: 25% gross outliers must be rejected */
    {
        static double X3b[3 * 40], uvb[2 * 40];
        int m = n;
        memcpy(X3b, X3, sizeof(X3b));
        memcpy(uvb, uv, sizeof(uvb));
        for (i = 0; i < 10; i++) {
            uvb[2 * (4 * i)] += 40.0 + 90.0 * urand();
            uvb[2 * (4 * i) + 1] -= 35.0 + 80.0 * urand();
        }
        CHECK(mv_resect_robust(&est, X3b, uvb, &m) == MV_OK,
              "resect: robust fit succeeds");
        CHECK(m >= 28 && m <= 32,
              "resect: outliers trimmed, inliers kept");
        {
            double e = fabs(est.K[0] - 720.0) / 720.0;
            double C[3], Cg[3], d;
            mv_cam_center(C, &est);
            mv_cam_center(Cg, &gt);
            d = sqrt((C[0] - Cg[0]) * (C[0] - Cg[0])
                     + (C[1] - Cg[1]) * (C[1] - Cg[1])
                     + (C[2] - Cg[2]) * (C[2] - Cg[2]));
            CHECK(e < 1e-3 && d < 1e-3,
                  "resect: robust recovery exact after trim");
        }
    }

    /* noisy resection: sigma = 0.5 px, focal within 2%, center 3 cm */
    {
        static double X3n[3 * 40], uvn[2 * 40];
        int m = n;
        memcpy(X3n, X3, sizeof(X3n));
        for (i = 0; i < n; i++) {
            uvn[2 * i] = uv[2 * i] + 0.5 * grand();
            uvn[2 * i + 1] = uv[2 * i + 1] + 0.5 * grand();
        }
        CHECK(mv_resect_robust(&est, X3n, uvn, &m) == MV_OK,
              "resect: noisy fit succeeds");
        {
            double C[3], Cg[3], d;
            mv_cam_center(C, &est);
            mv_cam_center(Cg, &gt);
            d = sqrt((C[0] - Cg[0]) * (C[0] - Cg[0])
                     + (C[1] - Cg[1]) * (C[1] - Cg[1])
                     + (C[2] - Cg[2]) * (C[2] - Cg[2]));
            printf("      noisy: fx err %.2f%%, center err %.1f mm\n",
                   100.0 * fabs(est.K[0] - 720.0) / 720.0, 1000.0 * d);
            CHECK(fabs(est.K[0] - 720.0) / 720.0 < 0.02 && d < 0.03,
                  "resect: noisy recovery within 2% / 3 cm");
        }
    }
}

/* ---- full pipeline on rendered imagery ------------------------------ */

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

static void test_pipeline(void)
{
    static unsigned char img1[640 * 480], img2[640 * 480],
        imgp[640 * 480];
    mv_camera c1, c2, gt, est;
    unsigned long long seed = 3;
    int nanchor = 0, ninl = 0;

    {
        static const double ex[3] = { 1, 0, 0 }, ey[3] = { 0, 1, 0 },
                            ez[3] = { 0, 0, 1 };
        double t0[3];
        make_tex(texA, 160, 80);
        make_tex(texB, 140, 160);
        make_tex(texC, 160, 80);
        t0[0] = -1.6; t0[1] = -0.8; t0[2] = 4.4;
        set_plane(&pls[0], texA, 160, 80, 0.02, ex, ey, t0); /* back */
        t0[0] = -1.4; t0[1] = 0.8; t0[2] = 2.6;
        set_plane(&pls[1], texB, 140, 160, 0.02, ex, ez, t0); /* floor */
        t0[0] = 1.4; t0[1] = -0.8; t0[2] = 2.6;
        set_plane(&pls[2], texC, 160, 80, 0.02, ez, ey, t0); /* side */
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
    /* the "photo": different K, rolled, displaced -- but similar
     * magnification (the v1 single-scale feature limitation) */
    make_photo_cam(&gt);

    mv_render_scene(img1, NULL, 640, 480, &c1, pls, 3, 10, 1.0, &seed);
    mv_render_scene(img2, NULL, 640, 480, &c2, pls, 3, 10, 1.0, &seed);
    mv_render_scene(imgp, NULL, 640, 480, &gt, pls, 3, 10, 1.0, &seed);

    CHECK(mv_photo_register(&est, &c1, &c2, img1, img2, imgp, 640, 480,
                            640, 480, &nanchor, &ninl) == MV_OK,
          "photo: registration succeeds");
    printf("      anchors %d, resection inliers %d\n", nanchor, ninl);
    CHECK(nanchor >= 30, "photo: anchor bank populated");
    CHECK(ninl >= 15, "photo: resection inliers");
    {
        double C[3], Cg[3], d, tr, ang, Rt[9], D[9];
        mv_cam_center(C, &est);
        mv_cam_center(Cg, &gt);
        d = sqrt((C[0] - Cg[0]) * (C[0] - Cg[0])
                 + (C[1] - Cg[1]) * (C[1] - Cg[1])
                 + (C[2] - Cg[2]) * (C[2] - Cg[2]));
        mv_mat_transpose(Rt, gt.R, 3, 3);
        mv_mat_mul(D, est.R, Rt, 3, 3, 3);
        tr = (D[0] + D[4] + D[8] - 1.0) / 2.0;
        if (tr > 1.0)
            tr = 1.0;
        if (tr < -1.0)
            tr = -1.0;
        ang = acos(tr) * 180.0 / 3.14159265358979324;
        printf("      recovered fx %.1f (true 720), center err %.1f mm, "
               "rot err %.2f deg\n", est.K[0], 1000.0 * d, ang);
        CHECK(fabs(est.K[0] - 720.0) / 720.0 < 0.05,
              "photo: unknown focal recovered within 5%");
        /* the ~0.1 m position error is the measured cost of descriptor
         * matching under a 10% magnification mismatch (single-scale v1
         * features); bundle adjustment over the registered set is the
         * documented next stage for tightening it */
        CHECK(d < 0.15, "photo: camera position within 15 cm");
        CHECK(ang < 1.5, "photo: orientation within 1.5 deg");
    }
}

static void test_halfturn(void)
{
    /* a photo camera rotated ~180 deg about a tilted axis from the rig
     * frame: exercises R_rodrigues's theta~pi branch inside the resect
     * polish (which previously collapsed the rotation to zero) */
    mv_camera gt, est;
    static double X3[3 * 40], uv[2 * 40];
    int i, n = 40, m;
    double ca = cos(3.0), sa = sin(3.0); /* ~172 deg */
    double Rz[9] = { ca, -sa, 0, sa, ca, 0, 0, 0, 1 }, Rx[9], cb, sb;
    cb = cos(0.2); sb = sin(0.2);
    Rx[0]=1;Rx[1]=0;Rx[2]=0;Rx[3]=0;Rx[4]=cb;Rx[5]=-sb;Rx[6]=0;Rx[7]=sb;Rx[8]=cb;
    mv_cam_set_K(&gt, 700.0, 705.0, 320.0, 240.0);
    memset(gt.k, 0, sizeof(gt.k));
    mv_mat_mul(gt.R, Rz, Rx, 3, 3, 3);
    {
        double C[3] = { 0.1, 0.05, -0.6 };
        int j;
        for (j = 0; j < 3; j++)
            gt.t[j] = -(gt.R[j*3+0]*C[0] + gt.R[j*3+1]*C[1] + gt.R[j*3+2]*C[2]);
    }
    for (i = 0; i < n; i++) {
        cloud_point(X3 + 3 * i, i);
        mv_cam_project(uv + 2 * i, &gt, X3 + 3 * i);
    }
    m = n;
    CHECK(mv_resect_robust(&est, X3, uv, &m) == MV_OK,
          "halfturn: resect succeeds near 180 deg");
    {
        double Rt[9], D[9], tr, ang;
        mv_mat_transpose(Rt, gt.R, 3, 3);
        mv_mat_mul(D, est.R, Rt, 3, 3, 3);
        tr = (D[0]+D[4]+D[8]-1.0)/2.0;
        if (tr>1.0) tr=1.0; if (tr<-1.0) tr=-1.0;
        ang = acos(tr) * 180.0 / 3.14159265358979324;
        CHECK(ang < 0.5, "halfturn: rotation recovered (pi branch works)");
    }
}

int main(void)
{
    test_resection();
    test_pipeline();
    test_halfturn();
    if (failures) {
        printf("\n%d FAILURE(S)\n", failures);
        return 1;
    }
    printf("\nall photo tests passed\n");
    return 0;
}
