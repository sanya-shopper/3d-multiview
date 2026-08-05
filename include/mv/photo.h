#ifndef MV_PHOTO_H
#define MV_PHOTO_H

/* Paper: doc/multiview.tex, section "Auxiliary photo collections". */

/* Registration of a casual photograph -- arbitrary viewpoint, UNKNOWN
 * camera -- into the rig's metric world frame, following the recipe of
 * that section: the calibrated pair triangulates a bank of 3D anchor
 * points carrying feature descriptors; the photo's features match into
 * the bank; the 2D-3D correspondences are resected by normalized DLT
 * into a full projection matrix; and its RQ decomposition yields the
 * photo's own K, R, t.  The rig anchors the metric scale, so the
 * classic SfM scale ambiguity never arises.  Scale caveat inherited
 * from mv/feat.h: the v1 features are single-scale, so the photo must
 * image the scene at roughly the rig cameras' magnification. */

#include "mv/cam.h"

/* Normalized DLT resection: projection matrix P (3x4, row-major) from
 * n >= 6 world points X3 (x0,y0,z0,x1,...) and pixels uv (u0,v0,...).
 * Points and pixels are Hartley-normalized internally.  P is scaled so
 * the decomposed K has K[8] = 1. */
int mv_resect_dlt(double P[12], const double *X3, const double *uv, int n);

/* Decompose P = K [R | t] by RQ: K upper-triangular with positive
 * diagonal and K[8] = 1, R a proper rotation, cam->k zeroed.
 * MV_ERR if the left 3x3 block is singular. */
int mv_cam_from_projection(mv_camera *cam, const double P[12]);

/* Robust resection: DLT fit, drop correspondences beyond 3x the median
 * reprojection error, refit (two rounds -- the same trimming discipline
 * as mv_fundamental_8point_robust).  X3 and uv are compacted in place;
 * *n becomes the surviving count.  cam receives the decomposed camera. */
int mv_resect_robust(mv_camera *cam, double *X3, double *uv, int *n);

/* The full pipeline: register photograph imgp (wp x hp) against the
 * calibrated pair c1/c2 viewing img1/img2 (both w x h).  On success
 * photo holds the recovered camera (K, R, t; distortion zero -- casual
 * photos are treated as pinhole in v1).  nanchor (may be NULL) receives
 * the size of the triangulated anchor bank, ninl the resection inliers.
 * MV_ERR if any stage starves (too few features, matches, or inliers:
 * the photo must overlap the pair's common field of view). */
int mv_photo_register(mv_camera *photo, const mv_camera *c1,
                      const mv_camera *c2, const unsigned char *img1,
                      const unsigned char *img2, const unsigned char *imgp,
                      int w, int h, int wp, int hp, int *nanchor,
                      int *ninl);

#endif /* MV_PHOTO_H */
