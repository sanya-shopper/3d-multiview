#ifndef MV_IMG_H
#define MV_IMG_H

/* Paper: utility I/O (no dedicated section); PLY/PGM appear in
 * "From depth to models" and throughout the demos. */

/* Minimal binary PGM (P5, 8-bit) image I/O — enough to feed the stereo
 * module from files produced by any common tool. */

/* Allocates *data with malloc; caller frees. Returns MV_ERR on parse
 * failure or maxval > 255. */
int mv_pgm_read(const char *path, unsigned char **data, int *w, int *h);

int mv_pgm_write(const char *path, const unsigned char *data, int w, int h);

#endif /* MV_IMG_H */
