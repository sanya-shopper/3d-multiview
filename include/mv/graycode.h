#ifndef MV_GRAYCODE_H
#define MV_GRAYCODE_H

/* Active calibration target: a computer display we control shows a
 * synchronized sequence of binary Gray-code stripe patterns. Each camera
 * pixel that sees the display decodes, over time, the display coordinate it
 * observes — dense display<->camera correspondences with no corner
 * detector. Each bit is shown together with its inverse so decoding is a
 * per-pixel comparison, robust to illumination and display gamma.
 *
 * The display plane with its known pixel pitch (metres/pixel) then acts as
 * a metric planar target for mv_calib_planar. */

/* Number of bits needed to encode coordinates 0..extent-1. */
int mv_graycode_bits(int extent);

unsigned mv_gray_encode(unsigned v);
unsigned mv_gray_decode(unsigned g);

/* Fill img (w x h, 0/255) with stripe pattern for the given bit.
 * axis: 0 encodes the x coordinate, 1 the y coordinate.
 * bit: 0 is the most significant of mv_graycode_bits(extent) bits.
 * inverse: nonzero renders the complement frame. */
void mv_graycode_frame(unsigned char *img, int w, int h, int axis, int bit,
                       int nbits, int inverse);

/* Decode one axis from captured frames. frames[b] and inv_frames[b] are the
 * camera images (npix pixels each) of bit b and its inverse, b = 0 (MSB)
 * .. nbits-1. For each pixel: coord[p] = decoded display coordinate,
 * valid[p] = 1 if every bit had at least min_contrast gray-level difference
 * between frame and inverse (0 otherwise, coord[p] = -1). Returns the
 * number of valid pixels. */
int mv_graycode_decode(int *coord, unsigned char *valid,
                       const unsigned char *const *frames,
                       const unsigned char *const *inv_frames,
                       int nbits, int npix, int min_contrast);

#endif /* MV_GRAYCODE_H */
