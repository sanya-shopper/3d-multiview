/* RANSAC plane extraction: synthetic room (floor + two walls) with
 * Gaussian noise and gross uniform outliers; exact recovery at
 * sigma = 0, accuracy bands at sigma = 2 mm, determinism.
 * OWNERSHIP (parallel build): TSDF-completion work item ONLY. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mv/core.h"
#include "mv/mat.h"
#include "mv/plane.h"

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

static unsigned long long rng = 20260811ULL;
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

/* room: floor y = 0, wall x = -2, wall z = 5; canonical (d >= 0) truth */
static const double TRUE_PL[3][4] = {
    { 0.0, 1.0, 0.0, 0.0 },  /* floor  */
    { 1.0, 0.0, 0.0, 2.0 },  /* wall x = -2 */
    { 0.0, 0.0, -1.0, 5.0 }  /* wall z = 5 */
};

#define NPP 400  /* points per plane */
#define NOUT 500 /* gross outliers (~30%) */
#define NTOT (3 * NPP + NOUT)

/* pts gets NTOT xyz points: NPP per plane plus NOUT uniform in the
 * room box; truth[i] = plane index or -1 for outliers */
static void make_room(double *pts, int *truth, double sigma)
{
    int i, a;
    for (i = 0; i < NPP; i++) { /* floor */
        pts[3 * i + 0] = -2.0 + 4.0 * urand();
        pts[3 * i + 1] = 0.0;
        pts[3 * i + 2] = 1.0 + 4.0 * urand();
        truth[i] = 0;
    }
    for (i = NPP; i < 2 * NPP; i++) { /* wall x = -2 */
        pts[3 * i + 0] = -2.0;
        pts[3 * i + 1] = 2.5 * urand();
        pts[3 * i + 2] = 1.0 + 4.0 * urand();
        truth[i] = 1;
    }
    for (i = 2 * NPP; i < 3 * NPP; i++) { /* wall z = 5 */
        pts[3 * i + 0] = -2.0 + 4.0 * urand();
        pts[3 * i + 1] = 2.5 * urand();
        pts[3 * i + 2] = 5.0;
        truth[i] = 2;
    }
    if (sigma > 0.0)
        for (i = 0; i < 3 * NPP; i++)
            for (a = 0; a < 3; a++)
                pts[3 * i + a] += sigma * grand();
    for (i = 3 * NPP; i < NTOT; i++) { /* gross outliers in the box */
        pts[3 * i + 0] = -2.0 + 4.0 * urand();
        pts[3 * i + 1] = 2.5 * urand();
        pts[3 * i + 2] = 1.0 + 4.0 * urand();
        truth[i] = -1;
    }
}

/* angle between unit normals, degrees, sign-insensitive */
static double normal_deg(const double a[4], const double b[4])
{
    double c = fabs(mv_dot(a, b, 3));
    if (c > 1.0)
        c = 1.0;
    return acos(c) * 180.0 / MV_PI;
}

/* |d| error after aligning the sign of b to a */
static double offset_err(const double a[4], const double b[4])
{
    double s = mv_dot(a, b, 3) < 0.0 ? -1.0 : 1.0;
    return fabs(a[3] - s * b[3]);
}

static void test_single_exact(void)
{
    static double pts[3 * NTOT];
    static unsigned char in[NTOT];
    static int truth[NTOT];
    double pl[4];
    int i, nin = 0;

    make_room(pts, truth, 0.0);
    /* floor points only, no noise, no outliers: exact recovery */
    CHECK(mv_plane_ransac(pts, NPP, 100, 0.001, 7u, pl, in) == MV_OK,
          "single: ransac succeeds");
    CHECK(normal_deg(pl, TRUE_PL[0]) < 1e-6
          && offset_err(pl, TRUE_PL[0]) < 1e-9,
          "single: exact plane at sigma = 0");
    for (i = 0; i < NPP; i++)
        nin += in[i];
    CHECK(nin == NPP, "single: every point is an inlier");
    CHECK(fabs(mv_norm(pl, 3) - 1.0) < 1e-12 && pl[3] >= 0.0,
          "single: unit normal, canonical d >= 0");
    /* degenerate input */
    CHECK(mv_plane_ransac(pts, 2, 100, 0.001, 7u, pl, in) == MV_ERR,
          "single: n < 3 rejected");
}

static void test_room(double sigma, double ang_tol_deg, double off_tol,
                      const char *tag)
{
    static double pts[3 * NTOT];
    static int truth[NTOT], labels[NTOT];
    double planes[4 * 5];
    int map[3]; /* true plane -> found plane index */
    int nplanes, i, p;
    char name[128];

    rng = 20260811ULL; /* same scene for every tolerance run */
    make_room(pts, truth, sigma);
    CHECK(mv_planes_extract(pts, NTOT, 5, 100, 300, 0.006, 42u, planes,
                            &nplanes, labels) == MV_OK
          && nplanes == 3, "room: exactly three planes found");
    if (nplanes != 3)
        return;
    for (p = 0; p < 3; p++) {
        double best = 1e9;
        int q, bq = -1;
        for (q = 0; q < nplanes; q++)
            if (normal_deg(TRUE_PL[p], planes + 4 * q) < best) {
                best = normal_deg(TRUE_PL[p], planes + 4 * q);
                bq = q;
            }
        map[p] = bq;
        snprintf(name, sizeof(name),
                 "room %s: plane %d normal %.4f deg (< %g), "
                 "offset %.2f mm (< %g mm)", tag, p, best, ang_tol_deg,
                 offset_err(TRUE_PL[p], planes + 4 * bq) * 1e3,
                 off_tol * 1e3);
        CHECK(best < ang_tol_deg
              && offset_err(TRUE_PL[p], planes + 4 * bq) < off_tol,
              name);
    }
    CHECK(map[0] != map[1] && map[1] != map[2] && map[0] != map[2],
          "room: the three matches are distinct");
    /* labelling quality: most true-plane points claimed by the matching
     * plane, gross outliers mostly unassigned */
    {
        int good = 0, out_free = 0;
        for (i = 0; i < NTOT; i++) {
            if (truth[i] >= 0)
                good += labels[i] == map[truth[i]];
            else
                out_free += labels[i] == -1;
        }
        snprintf(name, sizeof(name),
                 "room %s: %d/%d plane points labelled, %d/%d "
                 "outliers unassigned", tag, good, 3 * NPP, out_free,
                 NOUT);
        CHECK(good > 3 * NPP * 90 / 100 && out_free > NOUT * 95 / 100,
              name);
    }
}

static void test_determinism(void)
{
    static double pts[3 * NTOT];
    static int truth[NTOT], la[NTOT], lb[NTOT];
    double pa[4 * 5], pb[4 * 5];
    int na, nb;

    rng = 20260811ULL;
    make_room(pts, truth, 0.002);
    CHECK(mv_planes_extract(pts, NTOT, 5, 100, 300, 0.006, 99u, pa, &na,
                            la) == MV_OK, "determinism: run 1");
    CHECK(mv_planes_extract(pts, NTOT, 5, 100, 300, 0.006, 99u, pb, &nb,
                            lb) == MV_OK, "determinism: run 2");
    CHECK(na == nb
          && memcmp(pa, pb, (size_t)(4 * na) * sizeof(double)) == 0
          && memcmp(la, lb, sizeof(la)) == 0,
          "determinism: same seed, identical planes and labels");
    /* a different seed still finds the same three planes (values may
     * differ in the last bits; the count must not) */
    CHECK(mv_planes_extract(pts, NTOT, 5, 100, 300, 0.006, 7u, pb, &nb,
                            lb) == MV_OK && nb == na,
          "determinism: other seed, same plane count");
}

int main(void)
{
    test_single_exact();
    /* sigma = 0 with outliers: recovery to numerical precision */
    test_room(0.0, 0.01, 1e-4, "exact");
    /* sigma = 2 mm, ~30% outliers: 0.5 deg / 2 mm bands */
    test_room(0.002, 0.5, 0.002, "noisy");
    test_determinism();
    if (failures) {
        printf("\n%d FAILURE(S)\n", failures);
        return 1;
    }
    printf("\nall tests passed\n");
    return 0;
}
