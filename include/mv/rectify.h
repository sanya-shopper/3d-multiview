#ifndef MV_RECTIFY_H
#define MV_RECTIFY_H

#include "mv/cam.h"

/* Planar stereo rectification (compact algorithm of Fusiello, Trucco &
 * Verri): rotate both cameras about their centers onto a common orientation
 * whose x axis is the baseline, so epipolar lines become horizontal scan
 * lines with equal v in both images.
 *
 * Outputs the two rectified cameras (shared K, shared R, original centers)
 * and the 3x3 pixel homographies H1, H2 mapping original image coordinates
 * to rectified image coordinates. Distortion must be removed beforehand.
 * Returns MV_ERR if the camera centers coincide. */
int mv_rectify_pair(const mv_camera *c1, const mv_camera *c2,
                    mv_camera *r1, mv_camera *r2,
                    double H1[9], double H2[9]);

/* Baseline length |C2 - C1|. */
double mv_baseline(const mv_camera *c1, const mv_camera *c2);

/* Depth from disparity for a rectified pair: Z = f * B / d. */
double mv_disp_to_depth(double focal_px, double baseline, double disp);

#endif /* MV_RECTIFY_H */
