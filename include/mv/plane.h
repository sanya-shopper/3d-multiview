#ifndef MV_PLANE_H
#define MV_PLANE_H

/* Paper: doc/multiview.tex, section "The representation ladder", rung 4
 * (geometric primitives): RANSAC plane extraction [FischlerBolles1981].
 * OWNERSHIP (parallel build): this header, src/plane.c, and
 * tests/test_plane.c belong to the TSDF-completion work item ONLY. */

/* Robust plane fitting for 3-D point clouds (floors, walls, table
 * tops).  A plane is double[4] = (nx, ny, nz, d) with unit normal and
 * the convention n.x + d = 0; the stored form is canonicalized to
 * d >= 0, so (n, d) and (-n, -d) compare equal.  All routines are
 * deterministic for a given seed (local LCG, no rand()). */

/* Fit one plane to pts (n points, xyz interleaved) by RANSAC: iters
 * minimal 3-point hypotheses, inlier gate |n.p + d| <= tol, then a
 * total-least-squares refit (centroid + smallest eigenvector of the
 * scatter matrix) on the winning consensus set.  On success plane[]
 * holds the refit plane and inliers[i] (if inliers is non-NULL) is
 * 1/0 against the *refit* plane.  MV_ERR if n < 3, iters < 1,
 * tol <= 0, every sample was degenerate, or out of memory. */
int mv_plane_ransac(const double *pts, int n, int iters, double tol,
                    unsigned seed, double plane[4],
                    unsigned char *inliers);

/* Greedy multi-plane extraction: repeatedly RANSAC the not-yet-claimed
 * points, stop when a winner has fewer than min_inliers inliers (or
 * maxplanes planes were found).  planes holds 4*maxplanes doubles,
 * *nplanes the number filled in; labels (length n) gets the plane
 * index per point, -1 = unassigned.  MV_ERR on bad arguments or OOM;
 * finding zero planes is a success with *nplanes = 0. */
int mv_planes_extract(const double *pts, int n, int maxplanes,
                      int min_inliers, int iters, double tol,
                      unsigned seed, double *planes, int *nplanes,
                      int *labels);

#endif /* MV_PLANE_H */
