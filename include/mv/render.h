#ifndef MV_RENDER_H
#define MV_RENDER_H

/* Paper: doc/multiview.tex, section "The calibration subsystem" (validation tier: synthetic views
 * for developing the reader against exact ground truth). */

#include "mv/cam.h"

/* Synthetic view rendering: project a planar source image (the calibration
 * display) into a camera. The plane's own frame has X right, Y down, Z=0,
 * with metric scale = pitch metres per source pixel; cam holds the
 * plane-to-camera pose (Xc = R Xplane + t). Camera pixels that do not hit
 * the source rectangle get `background`; Gaussian intensity noise of
 * sigma gray levels is added everywhere (deterministic given seed). */
int mv_render_plane(unsigned char *out, int w, int h, const mv_camera *cam,
                    const unsigned char *src, int sw, int sh, double pitch,
                    unsigned char background, double sigma,
                    unsigned long long *seed);

#endif /* MV_RENDER_H */
