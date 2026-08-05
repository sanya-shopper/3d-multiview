#ifndef MV_REFINE_H
#define MV_REFINE_H

/* Paper: doc/multiview.tex -- Joint Levenberg-Marquardt refinement of intrinsics, distortion, and per-view poses.
 * OWNERSHIP (parallel build): this header, src/refine.c, and
 * tests/test_refine.c belong to the refine work item ONLY. */

#include "mv/calib.h"
#include "mv/cam.h"

/* Joint maximum-likelihood polish of a plane-based calibration.
 *
 * mv_calib_planar alternates Zhang's closed form with a linear radial
 * step; it minimizes algebraic surrogates, not the reprojection error,
 * and leaves ~17% systematic residual above the noise floor (the T5
 * tripwire of doc/multiview.tex).  This routine takes that output as
 * the starting point and runs Levenberg-Marquardt on the true cost:
 * the sum of squared pixel residuals between the DISTORTED projection
 * (mv_cam_project conventions) and the raw observations views[i].img.
 *
 * Parameters refined jointly, p = 6 + 6*nviews:
 *   shared    fx, fy, cx, cy (skew fixed at 0), k1, k2
 *   per view  rotation (Rodrigues 3-vector) and translation t (3)
 * Higher distortion terms p1, p2, k3 are held fixed at cams[0].k[2..4].
 *
 * K, cams (nviews of them), and k_radial are refined IN PLACE from the
 * initial values they already hold; k_radial may be NULL, in which case
 * the initial k1, k2 are taken from cams[0].k.  On success every
 * cams[i] is a fully usable camera (K copied, R, t, k set).
 *
 * The normal equations are solved by the Schur complement over the
 * per-view 6x6 blocks, so the cost per iteration is linear in nviews.
 * Returns MV_ERR only on structural failure (bad inputs, or a normal
 * system still singular at maximum damping); running out of iterations
 * at a cost plateau is convergence, not failure. */
int mv_calib_refine(double K[9], mv_camera *cams, double k_radial[2],
                    const mv_calib_view *views, int nviews);

#endif /* MV_REFINE_H */
