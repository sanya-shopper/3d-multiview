#ifndef MV_STEREO_H
#define MV_STEREO_H

/* Dense correspondence on a *rectified* pair by block matching.
 * Local SAD cost, winner-take-all with sub-pixel parabolic refinement.
 *
 * left, right: 8-bit grayscale images, w x h, row-major
 * disp:        output, w*h floats; matched disparity d >= 0 such that
 *              right pixel (u - d, v) corresponds to left pixel (u, v);
 *              -1 where no match was evaluated (borders, search cut short)
 * halfwin:     matching window is (2*halfwin+1)^2
 * dmin, dmax:  inclusive disparity search range, 0 <= dmin <= dmax */
int mv_stereo_sad(float *disp, const unsigned char *left,
                  const unsigned char *right, int w, int h,
                  int halfwin, int dmin, int dmax);

#endif /* MV_STEREO_H */
