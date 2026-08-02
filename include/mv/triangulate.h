#ifndef MV_TRIANGULATE_H
#define MV_TRIANGULATE_H

#include "mv/cam.h"

/* N-view linear (DLT) triangulation. Observations are assumed undistorted
 * (apply mv_cam_ray / undistortion upstream if the cameras have distortion).
 *
 * cams:   array of nviews camera pointers
 * uv:     flat array u0,v0,u1,v1,... of the observation in each view
 * X:      output Euclidean 3D point
 * Returns MV_ERR on degenerate input (nviews < 2 or point at infinity). */
int mv_triangulate(double X[3], const mv_camera *const *cams,
                   const double *uv, int nviews);

/* RMS reprojection error of X against the observations, in pixels. */
double mv_reproj_rms(const mv_camera *const *cams, const double *uv,
                     int nviews, const double X[3]);

#endif /* MV_TRIANGULATE_H */
