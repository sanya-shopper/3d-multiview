#ifndef MV_CALIB_H
#define MV_CALIB_H

/* Paper: doc/multiview.tex, section "Calibration from planar targets": homography, Zhang
 * constraints, radial alternation. */

#include "mv/cam.h"

/* Plane-based camera calibration (Zhang's method).
 *
 * Input: several views of a planar target (printed checkerboard, or a
 * display decoded via mv/graycode.h), each view being correspondences
 * between metric plane coordinates and observed pixels. Output: the shared
 * intrinsic matrix K, the per-view pose (extrinsics R, t of the target
 * plane in each shot), and optionally a linear estimate of radial
 * distortion k1, k2. */

typedef struct {
    const double *obj; /* n plane points X0,Y0,X1,Y1,... (metres, Z=0) */
    const double *img; /* n observed pixels u0,v0,... */
    int n;             /* n >= 4 */
} mv_calib_view;

/* Homography H (plane -> image) by normalized DLT, n >= 4 points. */
int mv_homography_dlt(double H[9], const double *xy_plane, const double *uv,
                      int n);

/* Zhang calibration. nviews >= 3, or >= 2 with zero_skew nonzero.
 * On success K holds the intrinsics and cams[i] is a fully usable camera
 * for view i (K, R, t set; distortion zeroed). If k_radial is non-NULL, a
 * linear least-squares estimate of k1, k2 is computed afterwards and also
 * written into every cams[i].k. */
int mv_calib_planar(double K[9], mv_camera *cams, double k_radial[2],
                    const mv_calib_view *views, int nviews, int zero_skew);

/* RMS reprojection error over all views/points, in pixels. */
double mv_calib_reproj_rms(const mv_camera *cams, const mv_calib_view *views,
                           int nviews);

#endif /* MV_CALIB_H */
