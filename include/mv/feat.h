#ifndef MV_FEAT_H
#define MV_FEAT_H

/* Paper: doc/multiview.tex -- Natural-feature detection and robust matching.
 * OWNERSHIP (parallel build): this header, src/feat.c, tests/test_feat.c,
 * and the track module (include/mv/track.h, src/track.c,
 * tests/test_track.c) belong to the SfM connective-tissue work item
 * ONLY. */

#include "mv/cam.h"

/* Natural-texture interest points for two-view correspondence.  This is
 * the front end that feeds mv_fundamental_8point_robust (mv/epipolar.h):
 * detect corners in each view, match descriptors, hand the matched pixel
 * pairs to the robust epipolar fit.
 *
 * Invariance, honestly stated:
 *  - affine lighting: YES (descriptor normalized to zero mean / unit
 *    variance over the patch)
 *  - in-plane rotation: YES (sampling grid rotated to the local dominant
 *    gradient direction before extraction)
 *  - scale: YES via mv_feat_detect_ms (half-octave image pyramid;
 *    magnification differences up to ~sqrt(2)^(nlevels-1) between
 *    views).  mv_feat_detect alone remains single-scale: views must
 *    then image the scene at roughly the same magnification (within
 *    ~20-30%)
 *  - viewpoint / perspective foreshortening: only mildly, via the small
 *    patch footprint; strongly slanted surfaces will lose matches
 * It is a correspondence engine for narrow-baseline photo registration,
 * not a SIFT replacement [Lowe2004]. */

typedef struct {
    double u, v;     /* pixel position, sub-pixel refined, always in
                      * full-resolution (level-0) coordinates */
    double resp;     /* Shi-Tomasi corner response (min eigenvalue) */
    double desc[64]; /* normalized oriented patch descriptor, sampled at
                      * the detection scale */
    double scale;    /* detection scale: 1.0 = full resolution,
                      * sqrt(2)^L for pyramid level L.  Appended field:
                      * source-compatible with pre-scale callers. */
} mv_feature;

/* Detect up to maxf corners in an 8-bit grayscale image (row-major,
 * stride == w).  Response is the smaller eigenvalue of the 2x2
 * gradient-moment (structure) matrix accumulated over a 5x5 window from
 * Sobel gradients (Shi-Tomasi "good features to track"); candidates are
 * local maxima over a 5 px radius, above 1% of the strongest response,
 * far enough from the border for the descriptor footprint, refined to
 * sub-pixel by a 1-D quadratic fit per axis.  The strongest maxf survive,
 * sorted by descending response.
 *
 * Descriptor: 8x8 intensity samples on a grid of ~2 px spacing centered
 * on the corner, the grid rotated by the dominant local gradient
 * orientation (Gaussian-weighted mean gradient), then normalized to zero
 * mean and unit variance.
 *
 * Returns the number of features written (>= 0), or MV_ERR on bad
 * arguments.  Every feature gets scale = 1.0. */
int mv_feat_detect(mv_feature *out, int maxf, const unsigned char *img,
                   int w, int h);

/* Multi-scale detection over a half-octave image pyramid: level L is
 * the input downsampled by sqrt(2)^L (3x3 binomial prefilter before
 * each sqrt(2) resample, so coarse levels are anti-aliased), and
 * mv_feat_detect runs on each level.  Positions are mapped back to
 * level-0 pixel coordinates; each descriptor is sampled on its own
 * level, i.e. at the detection scale, which makes descriptors
 * comparable between views that image the scene at different
 * magnifications (up to ~sqrt(2)^(nlevels-1) apart).
 *
 * The same physical corner may be reported at more than one scale --
 * deliberately, as in SIFT [Lowe2004]: the adjacent-scale descriptors
 * see genuinely different footprints, and the mutual ratio-test match
 * selects whichever scale pairing corresponds across views.  Features
 * from all levels compete for the maxf slots by descending response
 * (ties broken by scale, then position -- deterministic).
 *
 * nlevels == 1 is exactly mv_feat_detect.  Levels whose downsampled
 * image becomes too small for the descriptor footprint are skipped.
 * Returns the number of features written (>= 0), or MV_ERR on bad
 * arguments. */
int mv_feat_detect_ms(mv_feature *out, int maxf, const unsigned char *img,
                      int w, int h, int nlevels);

/* Match f1 against f2 by descriptor SSD.  A feature i in f1 is matched
 * to its nearest neighbour j in f2 iff
 *  - Lowe's ratio test [Lowe2004]: SSD(i,j) < ratio * SSD(i, second
 *    nearest).  `ratio` compares SSD values, so an L2-distance ratio r
 *    corresponds to r*r here (e.g. 0.75 L2 -> 0.56 SSD); and
 *  - the match is mutual: i is also the nearest neighbour of j in f1.
 * idx2[i] = matched index in f2, or -1 if unmatched (n1 entries always
 * written).  Returns the number of accepted matches (>= 0), or MV_ERR
 * on bad arguments.  Expect a residual outlier fraction; feed the pairs
 * to mv_fundamental_8point_robust for geometric verification. */
int mv_feat_match(int *idx2, const mv_feature *f1, int n1,
                  const mv_feature *f2, int n2, double ratio);

#endif /* MV_FEAT_H */
