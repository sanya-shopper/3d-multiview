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

#endif /* MV_EPIPOLAR_H */
