#ifndef MV_ROT_H
#define MV_ROT_H

/* Paper: doc/multiview.tex, section "Rotations".
 *
 * Canonical rotation kernels. These were previously copied across
 * calib, refine, bundle, photo, rigcalib, and scenecloud -- and the
 * copies diverged into bugs: a Rodrigues log that collapsed to zero at
 * theta ~ pi, and a nearest-rotation projection missing its
 * determinant flip (so it could return a reflection). One tested
 * implementation lives here; see CODING_STANDARDS.md "one kernel". */

/* Rodrigues exponential: axis-angle vector r (direction = axis,
 * magnitude = angle in radians) -> rotation matrix R (row-major).
 * Exact through theta = 0. */
void mv_rot_exp(double R[9], const double r[3]);

/* Rodrigues logarithm: rotation matrix R -> axis-angle vector r.
 * Handles both degenerate regimes correctly: theta ~ 0 (from the
 * antisymmetric part) and theta ~ pi (axis from the diagonal, where
 * sin(theta) -> 0 and the naive formula divides by zero). */
void mv_rot_log(double r[3], const double R[9]);

/* Nearest rotation to an arbitrary 3x3 matrix M (row-major), by SVD:
 * R = U V^T with a determinant flip so the result is a proper rotation
 * (det +1), never a reflection. Returns MV_ERR if the SVD fails. */
int mv_rot_project(double R[9], const double M[9]);

#endif /* MV_ROT_H */
