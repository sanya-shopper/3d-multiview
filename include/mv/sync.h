#ifndef MV_SYNC_H
#define MV_SYNC_H

/* Paper: doc/multiview.tex, subsection "Software synchronization:
 * measuring and correcting the clocks" -- the "Correct" layer:
 * interpolate tracked detections onto a common clock before pairing and
 * triangulating, turning the first-order phantom-depth error f v dt / Z
 * into the second-order curvature residual a dt^2 / 8.
 * OWNERSHIP (parallel build): this header, src/sync.c, and
 * tests/test_sync.c belong to the SYNC-CORRECTION work item ONLY. */

#include "mv/cam.h"

/* One camera's timestamped 2-D detection track.  Timestamps must already
 * be expressed on the COMMON clock: the producer is the per-camera clock
 * model of tools/hub_clock.h -- apply hub_clock_map(camid, t_cam) to
 * each capture time before building the track.  t strictly increasing. */
typedef struct {
    const double *t;  /* n timestamps, strictly increasing */
    const double *uv; /* n detections, flat u0,v0,u1,v1,... */
    int n;
} mv_sync_track;

/* Evaluate the track at query time tq by linear interpolation between
 * the bracketing samples.  Small extrapolation is allowed: tq may lie
 * outside [t[0], t[n-1]] by up to extrap times the adjacent sample
 * spacing (extrap = 0.5 allows half a frame period past either end);
 * the boundary segment's line is extended.  A single-sample track only
 * answers tq == t[0].  Returns MV_ERR on bad arguments, non-increasing
 * t, extrap < 0 or non-finite, or tq outside the allowed window (a NaN
 * tq is outside every window). */
int mv_sync_interp(double uv[2], const mv_sync_track *tr, double tq,
                   double extrap);

/* Software-synchronized pairing: for each sample time of track a that
 * lies within track b's coverage (widened by the extrap margin above),
 * pair a's detection as observed with b interpolated at that instant.
 * Output arrays are caller-allocated for up to a->n entries, indexed by
 * pair:
 *   t_out  common-clock time of the pair (may be NULL)
 *   uva    a's detection, flat u,v
 *   uvb    b's interpolated detection, flat u,v
 * Returns the number of pairs (>= 0; 0 when the tracks do not overlap
 * in time) or MV_ERR on bad arguments. */
int mv_sync_pair(double *t_out, double *uva, double *uvb,
                 const mv_sync_track *a, const mv_sync_track *b,
                 double extrap);

/* mv_sync_pair followed by two-view DLT triangulation of every pair
 * (mv/triangulate.h; observations assumed undistorted).  X receives up
 * to a->n Euclidean 3-D points x,y,z flat; t_out (may be NULL) their
 * common-clock times.  Pairs whose triangulation is degenerate are
 * dropped.  Returns the number of points (>= 0) or MV_ERR on bad
 * arguments. */
int mv_sync_pair_triangulate(double *X, double *t_out,
                             const mv_camera *ca, const mv_camera *cb,
                             const mv_sync_track *a,
                             const mv_sync_track *b, double extrap);

#endif /* MV_SYNC_H */
