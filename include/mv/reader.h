#ifndef MV_READER_H
#define MV_READER_H

/* Paper: doc/multiview.tex, section "The calibration subsystem" (reader algorithms, record
 * consumers). */

/* Calibration pattern reader (spec v1, doc section "The calibration
 * subsystem"): consumes one camera frame, returns identified corners in
 * pattern coordinates, the decoded frame counter, and the fitted
 * pattern->image homography. Blind: no pose prior needed.
 * Not reentrant: uses static working buffers; one reader call at a time.
 *
 * Setup & run: needs only an 8-bit mono frame (any source; the demo
 * renders one synthetically). Build and exercise end to end with
 *     make demo_patternsim && ./demo_patternsim
 * (renders the spec-v1 screen under seven poses, reads it blind, scores
 * corners/counter/orientation against ground truth, then runs Zhang
 * through the full pipeline; writes out_patternsim.pgm for inspection).
 * Set MV_READER_DEBUG=1 for per-stage tracing on stderr. Unit coverage:
 * make check (blind-read round trip among the 47 assertions). */

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

/* Coarse-tier reader (pattern spec v2, mv/pattern.h): same result
 * struct, different conventions -- id = 5*j + i over the 5x2 coarse
 * corner lattice, counter is 8-bit and wraps every 256 refreshes
 * (~4.3 s), H maps coarse pattern px -> image.  Orientation comes
 * from the checker phase plus the asymmetric mark L (validated
 * exactly: all 18 cell colors and all 18 mark states must match).
 * Internally escalates over 2x and 4x downsampling for close views.
 * Callers that accept both tiers should try mv_read_pattern first
 * (162 corners beat 10) and fall back to this. */
int mv_read_coarse(mv_read_result *res, const unsigned char *img,
                   int w, int h);

#endif /* MV_READER_H */
