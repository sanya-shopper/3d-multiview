#ifndef MV_READER_H
#define MV_READER_H

/* Paper: doc/multiview.tex, section "The calibration subsystem" (reader algorithms, record
 * consumers). */

/* Calibration pattern reader (spec v1, doc section "The calibration
 * subsystem"): consumes one camera frame, returns identified corners in
 * pattern coordinates, the decoded frame counter, and the fitted
 * pattern->image homography. Blind: no pose prior needed.
 * Not reentrant: uses static working buffers; one reader call at a time. */

#define MV_READ_MAXC 256

typedef struct {
    int n;                       /* identified corners */
    int id[MV_READ_MAXC];        /* pattern corner id = 18*j + i */
    double uv[2 * MV_READ_MAXC]; /* raw pixel coordinates (sub-pixel) */
    int counter_valid;
    unsigned counter;            /* decoded display frame counter */
    double counter_conf;         /* min sampling margin, 0..1 */
    double H[9];                 /* pattern px -> image (from all corners) */
    int rot;                     /* lattice rotations applied (0..3) */
} mv_read_result;

/* Returns MV_OK if at least 8 corners were identified (enough for a
 * homography with slack); MV_ERR otherwise. */
int mv_read_pattern(mv_read_result *res, const unsigned char *img,
                    int w, int h);

#endif /* MV_READER_H */
