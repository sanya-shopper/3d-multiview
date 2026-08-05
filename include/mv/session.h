#ifndef MV_SESSION_H
#define MV_SESSION_H

/* Paper: doc/multiview.tex -- Session records (the recordable-oracle grammar) and robust clock fitting.
 * OWNERSHIP (parallel build): this header, src/session.c, and
 * tests/test_session.c belong to the session work item ONLY. */

#include <stdio.h>

#include "mv/reader.h"

/* Calibration session records (doc, table "Calibration record schema"):
 * ASCII, one record per line, a 3-letter type tag followed by
 * whitespace-separated key=value fields; times are seconds as decimal
 * doubles on the host CLOCK_MONOTONIC; image coordinates are raw
 * distorted pixels. Emitted grammar (choices where the doc is silent
 * are noted):
 *
 *   VER spec=1 w=1920 h=1080 pitch_mm=0.2745 Td_ms=16.667 [cam=ID]
 *   FRM cam=ID k=K t=T
 *   CTR cam=ID t=T count=C conf=F
 *   CRN cam=ID t=T id=I u=U v=V
 *   PHS cam=ID t=T patch=P g=G
 *
 * Deviations/choices:
 *  - The doc keys per-frame records by camera frame index k; the read
 *    wrapper (mv_session_read) carries no driver frame index, so CTR/
 *    CRN/PHS are keyed by capture timestamp t instead, and FRM (which
 *    carries both) provides the k <-> t binding for consumers.
 *  - VER: w/h/Td_ms are the spec-v1 reference display constants
 *    below; pitch_mm is the ACTUAL pixel pitch of the display used
 *    (callers pass MV_SESSION_PITCH_MM when it is the reference
 *    display) -- a replay reconstructs metric scale from this field,
 *    so recording the constant while calibrating with a different
 *    physical pitch would scale every replayed length by their ratio
 *    (shipped bug: the live hub ran a 0.1133 mm laptop panel but
 *    stamped 0.2745, a 2.4x metric error for any consumer);
 *    cam= is appended only when a camera id is given (the doc's VER
 *    is rig-global).
 *  - CTR validity is encoded by presence, per the doc's "emitted when
 *    strip decodes": mv_session_read emits no CTR line when the
 *    counter failed to decode.
 *  - Field precision: t %.6f (microseconds), u/v %.3f (the doc's
 *    example precision; round-trips to 5e-4 px), g %.6f, conf %.4f. */

#define MV_SESSION_PITCH_MM 0.2745 /* spec-v1 reference display */
#define MV_SESSION_TD_MS 16.667

/* Session-start VER record. pitch_mm = the physical pixel pitch of
 * the display actually showing the pattern (MV_SESSION_PITCH_MM for
 * the reference display); cam_id may be NULL for a rig-global one. */
void mv_session_ver(FILE *f, int spec_version, double pitch_mm,
                    const char *cam_id);

/* Frame-arrival FRM record (t_mono = CLOCK_MONOTONIC seconds). */
void mv_session_frm(FILE *f, const char *cam_id, double t_mono,
                    int frame_idx);

/* Emit one frame's read result: CTR (only if the counter decoded),
 * one CRN per identified corner, and 5 PHS lines if phs is non-NULL
 * (phs[i] = normalized gray of patch i, mv_read_phases). */
void mv_session_read(FILE *f, const char *cam_id, double t_mono,
                     const mv_read_result *rr, const double phs[5]);

/* ---- parsing (round-trip consumer side) --------------------------- */

#define MV_SESSION_MAXF 12  /* max key=value fields per record */
#define MV_SESSION_KEYSZ 16
#define MV_SESSION_VALSZ 40

typedef struct {
    char type[4];                                /* "VER", "CRN", ... */
    int nf;
    char key[MV_SESSION_MAXF][MV_SESSION_KEYSZ];
    char val[MV_SESSION_MAXF][MV_SESSION_VALSZ];
} mv_session_record;

/* Parse one line (trailing newline ok). MV_ERR on malformed lines
 * (no 1-3 char type tag, a field without '=', or overflow). Blank
 * lines and lines starting with '#' parse to type "" with nf 0. */
int mv_session_parse(mv_session_record *rec, const char *line);

/* Value string for key, or NULL if the record has no such field. */
const char *mv_session_field(const mv_session_record *rec,
                             const char *key);

/* Numeric field: MV_OK and *out on success, MV_ERR if absent or not
 * a number. */
int mv_session_num(const mv_session_record *rec, const char *key,
                   double *out);

/* ---- reader-side samplers ----------------------------------------- */

/* Decode the spec version block through rr->H (fine tier only).
 * Brightness references come from the pattern's own checker cells;
 * black cell = bit 1, [1,0] marker required. Returns the version
 * (0..255), or -1 if any cell falls outside the image, contrast is
 * insufficient, or the marker check fails. */
int mv_read_version(const mv_read_result *rr, const unsigned char *img,
                    int w, int h);

/* Sample the 5 phase-patch centers through rr->H. phs[i] is the
 * patch's gray level normalized by the pattern's own white/black
 * references and clamped to [0,1]: 0 = black (patch bit set),
 * 1 = white. MV_ERR if any patch falls outside the image or the
 * references lack contrast. */
int mv_read_phases(double phs[5], const mv_read_result *rr,
                   const unsigned char *img, int w, int h);

/* ---- robust clock fit --------------------------------------------- */

/* Fit the doc's clock model (eq. "clock", exposure/row terms aside):
 * capture time t = phi + k*T for display counter value k, i.e. phi is
 * the camera clock offset in seconds and T the display frame period
 * in seconds per count. Theil-Sen median of pairwise slopes
 * (deterministically subsampled to 200 evenly spaced points when
 * n > 200), then a 3-MAD residual gate, then least squares on the
 * inliers. Needs n >= 2 with distinct k. MV_OK/MV_ERR. */
int mv_clock_fit(double *phi, double *T, const double *t,
                 const double *k, int n);

#endif /* MV_SESSION_H */
