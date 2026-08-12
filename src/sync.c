/* Software synchronization correction: interpolate tracked detections
 * onto a common clock before pairing/triangulation.  See mv/sync.h and
 * doc/multiview.tex "Software synchronization" for the error model.
 * OWNERSHIP (parallel build): SYNC-CORRECTION work item ONLY. */

#include <math.h>
#include <stddef.h>

#include "mv/core.h"
#include "mv/cam.h"
#include "mv/triangulate.h"
#include "mv/sync.h"

/* Preconditions shared by every entry point.  The O(n) monotonicity
 * sweep is cheap at track scale and its positive-form comparison also
 * rejects NaN timestamps (NaN > x is false). */
static int track_ok(const mv_sync_track *tr)
{
    int i;
    if (!tr || !tr->t || !tr->uv || tr->n < 1)
        return 0;
    for (i = 1; i < tr->n; i++)
        if (!(tr->t[i] > tr->t[i - 1]))
            return 0;
    return 1;
}

/* Allowed query window [lo, hi]: the track's span widened at each end
 * by extrap times the adjacent sample spacing.  A single-sample track
 * has no spacing to scale, so its window is the point t[0]. */
static void query_window(double *lo, double *hi, const mv_sync_track *tr,
                         double extrap)
{
    if (tr->n == 1) {
        *lo = *hi = tr->t[0];
        return;
    }
    *lo = tr->t[0] - extrap * (tr->t[1] - tr->t[0]);
    *hi = tr->t[tr->n - 1] + extrap * (tr->t[tr->n - 1] - tr->t[tr->n - 2]);
}

/* Linear interpolation at tq; caller guarantees track_ok and that tq is
 * inside the query window, so at most the boundary segment's line is
 * extended (never beyond the extrap margin). */
static void interp_at(double uv[2], const mv_sync_track *tr, double tq)
{
    int lo = 0, hi = tr->n - 1;
    double s;
    if (tr->n == 1) {
        uv[0] = tr->uv[0];
        uv[1] = tr->uv[1];
        return;
    }
    while (hi - lo > 1) {
        int mid = lo + (hi - lo) / 2;
        if (tr->t[mid] <= tq)
            lo = mid;
        else
            hi = mid;
    }
    s = (tq - tr->t[lo]) / (tr->t[lo + 1] - tr->t[lo]);
    uv[0] = tr->uv[2 * lo + 0] + s * (tr->uv[2 * lo + 2] - tr->uv[2 * lo + 0]);
    uv[1] = tr->uv[2 * lo + 1] + s * (tr->uv[2 * lo + 3] - tr->uv[2 * lo + 1]);
}

int mv_sync_interp(double uv[2], const mv_sync_track *tr, double tq,
                   double extrap)
{
    double lo, hi;
    if (!uv || !track_ok(tr) || !isfinite(extrap) || extrap < 0.0)
        return MV_ERR;
    query_window(&lo, &hi, tr, extrap);
    if (!(tq >= lo && tq <= hi))
        return MV_ERR;
    interp_at(uv, tr, tq);
    return MV_OK;
}

int mv_sync_pair(double *t_out, double *uva, double *uvb,
                 const mv_sync_track *a, const mv_sync_track *b,
                 double extrap)
{
    double lo, hi;
    int i, np = 0;
    if (!uva || !uvb || !track_ok(a) || !track_ok(b)
        || !isfinite(extrap) || extrap < 0.0)
        return MV_ERR;
    query_window(&lo, &hi, b, extrap);
    for (i = 0; i < a->n; i++) {
        double tq = a->t[i];
        if (!(tq >= lo && tq <= hi))
            continue;
        uva[2 * np + 0] = a->uv[2 * i + 0];
        uva[2 * np + 1] = a->uv[2 * i + 1];
        interp_at(uvb + 2 * np, b, tq);
        if (t_out)
            t_out[np] = tq;
        np++;
    }
    return np;
}

int mv_sync_pair_triangulate(double *X, double *t_out,
                             const mv_camera *ca, const mv_camera *cb,
                             const mv_sync_track *a, const mv_sync_track *b,
                             double extrap)
{
    const mv_camera *cams[2];
    double lo, hi, obs[4];
    int i, np = 0;
    if (!X || !ca || !cb || !track_ok(a) || !track_ok(b)
        || !isfinite(extrap) || extrap < 0.0)
        return MV_ERR;
    cams[0] = ca;
    cams[1] = cb;
    query_window(&lo, &hi, b, extrap);
    for (i = 0; i < a->n; i++) {
        double tq = a->t[i];
        if (!(tq >= lo && tq <= hi))
            continue;
        obs[0] = a->uv[2 * i + 0];
        obs[1] = a->uv[2 * i + 1];
        interp_at(obs + 2, b, tq);
        /* degenerate instant (e.g. point at infinity): drop the pair,
         * keep the rest of the track */
        if (mv_triangulate(X + 3 * np, cams, obs, 2) != MV_OK)
            continue;
        if (t_out)
            t_out[np] = tq;
        np++;
    }
    return np;
}
