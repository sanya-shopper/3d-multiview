#ifndef MV_EPIPOLAR_H
#define MV_EPIPOLAR_H

/* Paper: doc/multiview.tex, section "Two-view epipolar geometry". */

#include "mv/cam.h"

/* Two-view epipolar geometry. All fundamental matrices returned by this
 * module are normalized to unit Frobenius norm with a canonical sign
 * (largest-magnitude entry positive) so results are directly comparable. */

/* Analytic F from two calibrated cameras: F = [e2]_x P2 P1^+.
 * Satisfies x2^T F x1 = 0 for ideal (distortion-free) projections. */
void mv_fundamental_from_cams(double F[9], const mv_camera *c1,
                              const mv_camera *c2);

/* Normalized 8-point algorithm (Hartley) with rank-2 enforcement.
 * uv1, uv2: n corresponding pixels as flat arrays u0,v0,u1,v1,...; n >= 8.
 * Points are assumed undistorted. */
int mv_fundamental_8point(double F[9], const double *uv1, const double *uv2,
                          int n);

/* E = K2^T F K1. */
void mv_essential_from_fundamental(double E[9], const double F[9],
                                   const double K1[9], const double K2[9]);

/* Symmetric point-to-epipolar-line distance in pixels: the mean of the
 * distance from x2 to F x1 and from x1 to F^T x2. */
double mv_sym_epipolar_dist(const double F[9], const double uv1[2],
                            const double uv2[2]);

/* 8-point with two rounds of median-residual trimming: fit, drop
 * correspondences beyond 3x the median symmetric epipolar distance,
 * refit.  uv1/uv2 are compacted in place; *n becomes the surviving
 * count.  The tool for correspondence sets with a fat outlier tail
 * (dense code matching, feature matching). */
int mv_fundamental_8point_robust(double F[9], double *uv1, double *uv2,
                                 int *n);

/* Decompose an essential matrix into the relative pose of camera 2
 * with respect to camera 1: x2cam = R x1cam + t, |t| = 1 (the scale is
 * unobservable from E; multiply t by a known baseline to go metric).
 * Of the four (R, t) candidates the one placing the most
 * correspondences in front of BOTH cameras wins (cheirality).
 * x1, x2: n corresponding points in NORMALIZED image coordinates
 * (pixels through K^{-1}, distortion removed), flat x0,y0,x1,y1,...
 * Returns MV_ERR if no candidate wins clearly. */
int mv_essential_pose(double R[9], double t[3], const double E[9],
                      const double *x1, const double *x2, int n);

#endif /* MV_EPIPOLAR_H */
