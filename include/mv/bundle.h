#ifndef MV_BUNDLE_H
#define MV_BUNDLE_H

/* Paper: doc/multiview.tex, section "Auxiliary photo collections" --
 * bundle adjustment [Triggs2000] over registered views, the designated
 * tightening stage after mv_photo_register.
 * OWNERSHIP (parallel build): this header, src/bundle.c, and
 * tests/test_bundle.c belong to the bundle work item ONLY. */

#include "mv/cam.h"

/* One image observation: point pt seen by camera cam at raw pixel
 * (u, v).  Residuals are formed against mv_cam_project, so a camera's
 * k1/k2 (held at their input values) are applied; the casual-photo
 * pipeline hands in pinhole cameras with k zeroed. */
typedef struct {
    int cam;      /* camera index */
    int pt;       /* point index */
    double u, v;  /* observed pixel (raw) */
} mv_bundle_obs;

#define MV_BUNDLE_FIX_POSE   1u  /* camera pose held fixed (gauge anchor) */
#define MV_BUNDLE_FIX_K      2u  /* camera intrinsics held fixed */

/* Joint Levenberg-Marquardt refinement of ncam cameras and npts 3-D
 * points (X = x0,y0,z0,x1,...) from their initial values, minimizing
 * total squared reprojection error over the nobs observations.  Free
 * parameters per camera, unless masked by camflags[i]: pose as
 * Rodrigues r + t (6) and intrinsics fx, fy, cx, cy (4; skew fixed at
 * 0, distortion k1/k2 held at input values).  Each point contributes 3.
 * The gauge must be fixed by the caller flagging at least one camera
 * MV_BUNDLE_FIX_POSE (typically the two rig cameras, pose and K);
 * MV_ERR if none is.  The solve exploits the classic BA sparsity:
 * Schur complement over the per-point 3x3 blocks leaves a dense reduced
 * camera system of at most 10*ncam.  rms_out (may be NULL) receives
 * the final reprojection RMS in pixels. */
int mv_bundle_adjust(mv_camera *cams, const unsigned *camflags, int ncam,
                     double *X, int npts, const mv_bundle_obs *obs,
                     int nobs, double *rms_out);

#endif /* MV_BUNDLE_H */
