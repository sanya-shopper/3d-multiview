#ifndef MV_OPTIMAL_H
#define MV_OPTIMAL_H

/* Paper: doc/multiview.tex -- Optimal and weighted estimators: Hartley-Sturm triangulation, weighted circle fit.
 * OWNERSHIP (parallel build): this header, src/optimal.c, and
 * tests/test_optimal.c belong to the optimal work item ONLY. */

#include "mv/cam.h"

/* Hartley-Sturm optimal correction (refs/HartleySturm1997.pdf): replace
 * the measured pair (uv1, uv2) by the pair (out1, out2) EXACTLY
 * satisfying out2^T F out1 = 0 that minimizes the sum of squared image
 * distances to the measurements.  F convention: x2^T F x1 = 0 (as
 * produced by mv_fundamental_from_cams).  Points are assumed
 * undistorted.  Returns MV_ERR when an epipole coincides with a
 * measured point or the geometry is otherwise degenerate; the inputs
 * are then left uncorrected. */
int mv_epipolar_correct(double out1[2], double out2[2], const double F[9],
                        const double uv1[2], const double uv2[2]);

/* Optimal two-view triangulation: F is derived analytically from the
 * two cameras, the observations are corrected with
 * mv_epipolar_correct, and the corrected (exactly consistent) pair is
 * triangulated with mv_triangulate.  Degenerate epipole-at-point cases
 * fall back to plain triangulation of the measured pair. */
int mv_triangulate_optimal(double X[3], const mv_camera *c1,
                           const mv_camera *c2, const double uv1[2],
                           const double uv2[2]);

/* Mahalanobis-weighted known-radius circle fit for anisotropic
 * (depth-elongated) point noise.
 * pts:  n 2-D points x0,z0,x1,z1,... (ground-plane projections)
 * dirs: per-point UNIT 2-D direction of the dominant (depth) noise axis
 * sig:  per-point pair (sigma_along, sigma_across), same units as pts
 * Minimizes over the center c the sum of exact squared Mahalanobis
 * distances min_q (p_i - q)^T Sigma_i^-1 (p_i - q) over circle points
 * q, with Sigma_i = sa_i^2 a_i a_i^T + sc_i^2 b_i b_i^T (b = perp(a));
 * this also removes the curvature bias that the unweighted radial fit
 * incurs when the noise is comparable to the radius.  Gauss-Newton
 * (envelope theorem for the inner closest-point), initialized from the
 * unweighted fit.  Returns MV_ERR on invalid input or divergence. */
int mv_circle_fit_weighted(double center[2], const double *pts,
                           const double *dirs, const double *sig,
                           int n, double radius);

#endif /* MV_OPTIMAL_H */
