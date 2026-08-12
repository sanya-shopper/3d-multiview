#ifndef MV_TRACK_H
#define MV_TRACK_H

/* Paper: doc/multiview.tex -- SfM ledger items (ii) feature tracks and
 * (iii) incremental registration.
 * OWNERSHIP (parallel build): this header, src/track.c, and
 * tests/test_track.c belong to the SfM connective-tissue work item
 * ONLY (see mv/feat.h). */

#include "mv/cam.h"
#include "mv/feat.h"

/* The connective tissue between pairwise feature matching (mv/feat.h)
 * and multi-view reconstruction: pairwise matches are chained into
 * feature tracks (one track = one candidate 3D point seen in several
 * images), and the tracks drive incremental camera registration on top
 * of the calibrated rig pair via the robust resection of mv/photo.h. */

/* One feature track: len observations, observation k is feature
 * feat[k] of image img[k].  img[] is strictly ascending (a valid track
 * never contains two features of the same image -- see
 * mv_tracks_build).  img and feat are owned by the track set and freed
 * by mv_tracks_free. */
typedef struct {
    int len;
    int *img;
    int *feat;
} mv_track;

/* Pairwise match input for mv_tracks_build, same idx2 convention as
 * mv_feat_match: idx2 has nfeat[a] entries; idx2[i] = matched feature
 * index in image b, or -1 if feature i of image a is unmatched. */
typedef struct {
    int a, b;
    const int *idx2;
} mv_match_pair;

/* Chain pairwise matches into tracks by union-find over the features
 * of nimg images (nfeat[i] features in image i).  A connected
 * component of the match graph becomes a track iff it spans >= 2
 * images and is conflict-free.
 *
 * CONFLICT REJECTION: a component containing two different features of
 * the same image is DROPPED whole, not split.  A conflict proves at
 * least one constituent match is wrong, and the match graph carries no
 * evidence of which; any split would keep a guess.  Dropping loses a
 * few true observations but never launders a wrong association into
 * the reconstruction (the same discipline as the robust fits: prefer
 * starving to lying).
 *
 * *tracks receives a malloc'd array of *ntracks tracks, ordered by
 * their smallest global feature id (deterministic); each track's
 * observations are ordered by image index.  Free with mv_tracks_free.
 * Returns MV_OK, or MV_ERR on bad arguments, out-of-range match
 * indices, or allocation failure. */
int mv_tracks_build(mv_track **tracks, int *ntracks,
                    const int *nfeat, int nimg,
                    const mv_match_pair *pairs, int npairs);

void mv_tracks_free(mv_track *tracks, int ntracks);

/* Incremental structure-from-motion registration on top of a
 * calibrated seed pair.
 *
 * Inputs: per-image feature arrays feats[i] (nfeat[i] entries, only
 * u,v are read) for nimg images; the track set from mv_tracks_build;
 * seed1/seed2 = the two seed images (the rig pair), whose cams entries
 * carry known K, R, t on entry.  The remaining cams entries carry
 * their known intrinsics (per-image or all the same K); their pose is
 * ignored on entry.  Intrinsics are AUTHORITATIVE: cams[i].K and .k
 * are never changed, only poses are estimated.
 *
 * Loop: triangulate every track with >= 2 registered observations
 * (mv_triangulate N-view, gated on cheirality and a 2 px reprojection
 * RMS); repeatedly pick the unregistered image observing the most
 * triangulated tracks (ties -> lowest index: deterministic); resect it
 * with mv_resect_robust (mv/photo.h), keeping only the pose -- the
 * free-K DLT resection trades focal against depth, so its pose is an
 * initial value; then refine the registered set with a fixed-K bundle
 * adjustment (mv/bundle.h; seed pair fully fixed as the gauge anchor),
 * retriangulate with the enlarged set, and iterate until no image sees
 * >= 6 triangulated tracks or all are registered.  An image whose
 * resection fails is skipped permanently.  A final fixed-K bundle
 * consolidates the full registered set.
 *
 * Outputs: cams[i] pose for every registered image (K untouched);
 * registered[i] = 1/0 (nimg entries); X = 3*ntracks coords of the
 * triangulated track points, valid where xvalid[k] = 1; rms[i] = the
 * reprojection RMS in pixels of image i over the valid tracks it
 * observes, or -1.0 if unregistered or observing none.  The whole
 * procedure is deterministic.
 *
 * Returns the number of registered images (>= 2: the seeds), or
 * MV_ERR on bad arguments or allocation failure. */
int mv_sfm_register_incremental(mv_camera *cams, unsigned char *registered,
                                double *X, unsigned char *xvalid,
                                double *rms,
                                const mv_feature *const *feats,
                                const int *nfeat, int nimg,
                                const mv_track *tracks, int ntracks,
                                int seed1, int seed2);

#endif /* MV_TRACK_H */
