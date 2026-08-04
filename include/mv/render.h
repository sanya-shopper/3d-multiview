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

/* A textured finite plane placed in the world: plane frame has X right,
 * Y down, Z = 0; X_world = R * X_plane + t; metric size =
 * (tw * pitch) x (th * pitch). */
typedef struct {
    const unsigned char *tex;
    int tw, th;
    double pitch;
    double R[9], t[3];
} mv_plane;

/* Render a scene of n textured planes with a z-buffer. If zc is non-NULL
 * it receives the true camera-frame depth per pixel (HUGE_VALF where no
 * plane is hit) -- the ground truth the sim experiments score against. */
int mv_render_scene(unsigned char *img, float *zc, int w, int h,
                    const mv_camera *cam, const mv_plane *pls, int n,
                    unsigned char background, double sigma,
                    unsigned long long *seed);

/* Warp src into dst through the pixel homography H (src coords -> dst
 * coords, e.g. the rectifying H of mv_rectify_pair); inverse-mapped,
 * bilinear; dst pixels with no source get `fill`. */
int mv_warp_homography(unsigned char *dst, int dw, int dh,
                       const unsigned char *src, int sw, int sh,
                       const double H[9], unsigned char fill);

#endif /* MV_RENDER_H */
