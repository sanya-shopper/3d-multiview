#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "mv/core.h"
#include "mv/cam.h"
#include "mv/feat.h"
#include "mv/triangulate.h"
#include "mv/photo.h"
#include "mv/bundle.h"
#include "mv/track.h"

/* ---- track building ------------------------------------------------ */

/* Union-find with path halving. */
static int uf_find(int *parent, int x)
{
    while (parent[x] != x) {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    return x;
}

/* Union by smaller root id: the canonical root of a component is then
 * its minimum global feature id, which is what makes the output track
 * order deterministic. */
static void uf_union(int *parent, int a, int b)
{
    int ra = uf_find(parent, a), rb = uf_find(parent, b);
    if (ra == rb)
        return;
    if (ra < rb)
        parent[rb] = ra;
    else
        parent[ra] = rb;
}

int mv_tracks_build(mv_track **tracks, int *ntracks,
                    const int *nfeat, int nimg,
                    const mv_match_pair *pairs, int npairs)
{
    int *off, *parent, *imgof, *cnt, *start, *order;
    mv_track *out = NULL;
    long long total = 0;
    int N, i, g, p, nout = 0, ncomp = 0;

    if (!tracks || !ntracks || !nfeat || nimg < 2 || npairs < 0
        || (npairs > 0 && !pairs))
        return MV_ERR;
    *tracks = NULL;
    *ntracks = 0;
    for (i = 0; i < nimg; i++) {
        if (nfeat[i] < 0)
            return MV_ERR;
        total += nfeat[i];
    }
    if (total > INT_MAX / 2)
        return MV_ERR;
    N = (int)total;
    if (N == 0)
        return MV_OK;

    off = (int *)malloc((size_t)(nimg + 1) * sizeof(int));
    parent = (int *)malloc((size_t)N * sizeof(int));
    imgof = (int *)malloc((size_t)N * sizeof(int));
    cnt = (int *)calloc((size_t)N, sizeof(int));
    start = (int *)malloc((size_t)(N + 1) * sizeof(int));
    order = (int *)malloc((size_t)N * sizeof(int));
    if (!off || !parent || !imgof || !cnt || !start || !order)
        goto fail;
    off[0] = 0;
    for (i = 0; i < nimg; i++)
        off[i + 1] = off[i] + nfeat[i];
    for (i = 0; i < nimg; i++)
        for (g = off[i]; g < off[i + 1]; g++)
            imgof[g] = i;
    for (g = 0; g < N; g++)
        parent[g] = g;

    for (p = 0; p < npairs; p++) {
        int a = pairs[p].a, b = pairs[p].b;
        const int *idx2 = pairs[p].idx2;
        if (a < 0 || a >= nimg || b < 0 || b >= nimg || a == b || !idx2)
            goto fail;
        for (i = 0; i < nfeat[a]; i++) {
            int j = idx2[i];
            if (j < 0)
                continue;
            if (j >= nfeat[b])
                goto fail; /* corrupt match array */
            uf_union(parent, off[a] + i, off[b] + j);
        }
    }

    /* group members by root: counting sort keyed on the root id, so
     * components come out ordered by their minimum feature id and
     * members ascending within each */
    for (g = 0; g < N; g++)
        cnt[uf_find(parent, g)]++;
    start[0] = 0;
    for (g = 0; g < N; g++)
        start[g + 1] = start[g] + cnt[g];
    {
        int *fill = (int *)malloc((size_t)N * sizeof(int));
        if (!fill)
            goto fail;
        memcpy(fill, start, (size_t)N * sizeof(int));
        for (g = 0; g < N; g++)
            order[fill[uf_find(parent, g)]++] = g;
        free(fill);
    }
    for (g = 0; g < N; g++)
        if (cnt[g] >= 2)
            ncomp++;
    if (ncomp == 0)
        goto done;
    out = (mv_track *)calloc((size_t)ncomp, sizeof(mv_track));
    if (!out)
        goto fail;

    for (g = 0; g < N; g++) {
        int len = cnt[g], k, conflict = 0;
        mv_track *t;
        if (len < 2)
            continue;
        /* members are ascending, so images are non-decreasing: a
         * repeat is adjacent.  Conflict rejection: DROP the whole
         * component -- see mv/track.h for why not split. */
        for (k = 1; k < len; k++)
            if (imgof[order[start[g] + k]]
                == imgof[order[start[g] + k - 1]])
                conflict = 1;
        if (conflict)
            continue;
        t = &out[nout];
        t->img = (int *)malloc((size_t)len * sizeof(int));
        t->feat = (int *)malloc((size_t)len * sizeof(int));
        if (!t->img || !t->feat) {
            free(t->img);
            free(t->feat);
            t->img = NULL;
            t->feat = NULL;
            goto fail;
        }
        t->len = len;
        for (k = 0; k < len; k++) {
            int m = order[start[g] + k];
            t->img[k] = imgof[m];
            t->feat[k] = m - off[imgof[m]];
        }
        nout++;
    }
    *tracks = out;
    *ntracks = nout;
    out = NULL;
done:
    free(off);
    free(parent);
    free(imgof);
    free(cnt);
    free(start);
    free(order);
    return MV_OK;
fail:
    if (out)
        mv_tracks_free(out, nout);
    free(off);
    free(parent);
    free(imgof);
    free(cnt);
    free(start);
    free(order);
    return MV_ERR;
}

void mv_tracks_free(mv_track *tracks, int ntracks)
{
    int i;
    if (!tracks)
        return;
    for (i = 0; i < ntracks; i++) {
        free(tracks[i].img);
        free(tracks[i].feat);
    }
    free(tracks);
}

/* ---- incremental registration -------------------------------------- */

/* (Re)triangulate every track from the currently registered views.
 * A point is accepted only if the N-view DLT succeeds, it is finite,
 * it lies in front of every registered view observing it, and its
 * reprojection RMS beats 2 px (gates written in positive form so NaN
 * is rejected, not admitted).  tcams/tuv are caller scratch of nimg
 * entries. */
static void sfm_triangulate_all(double *X, unsigned char *xvalid,
                                const mv_camera *cams,
                                const unsigned char *registered,
                                const mv_feature *const *feats,
                                const mv_track *tracks, int ntracks,
                                const mv_camera **tcams, double *tuv)
{
    int k, o, j;
    for (k = 0; k < ntracks; k++) {
        const mv_track *t = &tracks[k];
        double P[3], rms;
        int nv = 0, good = 1;
        xvalid[k] = 0;
        for (o = 0; o < t->len; o++) {
            const mv_feature *f;
            if (!registered[t->img[o]])
                continue;
            f = &feats[t->img[o]][t->feat[o]];
            tcams[nv] = &cams[t->img[o]];
            tuv[2 * nv] = f->u;
            tuv[2 * nv + 1] = f->v;
            nv++;
        }
        if (nv < 2)
            continue;
        if (mv_triangulate(P, tcams, tuv, nv) != MV_OK)
            continue;
        for (j = 0; j < 3; j++)
            if (!isfinite(P[j]))
                good = 0;
        for (o = 0; o < nv && good; o++) {
            const mv_camera *c = tcams[o];
            double z = c->R[6] * P[0] + c->R[7] * P[1] + c->R[8] * P[2]
                       + c->t[2];
            if (!(z > 0.0))
                good = 0;
        }
        if (!good)
            continue;
        rms = mv_reproj_rms(tcams, tuv, nv, P);
        if (!(rms < 2.0))
            continue;
        X[3 * k] = P[0];
        X[3 * k + 1] = P[1];
        X[3 * k + 2] = P[2];
        xvalid[k] = 1;
    }
}

/* The observation of image i in track t, or -1. */
static int sfm_obs_of(const mv_track *t, int i)
{
    int o;
    for (o = 0; o < t->len; o++)
        if (t->img[o] == i)
            return o;
    return -1;
}

/* Fixed-intrinsics bundle refinement of the registered set: the seed
 * pair is the gauge anchor (pose and K fixed), every other registered
 * camera keeps its known K and frees only its pose, the valid track
 * points are free.  This is what makes the caller's intrinsics
 * authoritative: mv_resect_robust estimates a full P whose focal
 * trades against depth along the optical axis, so its pose alone --
 * re-anchored to the known K here -- is only an initial value.
 * bc/flags/camid are nimg-sized scratch, ptid ntracks-sized, Xs
 * 3*ntracks, obs sized for every observation of every track.  On any
 * failure the unrefined state is simply kept. */
static void sfm_bundle(mv_camera *cams, const unsigned char *registered,
                       double *X, const unsigned char *xvalid,
                       const mv_feature *const *feats,
                       const mv_track *tracks, int ntracks,
                       int seed1, int seed2, int nimg,
                       mv_camera *bc, unsigned *flags, int *camid,
                       int *ptid, double *Xs, mv_bundle_obs *obs)
{
    int i, k, nsub = 0, npts = 0, nobs = 0;

    for (i = 0; i < nimg; i++) {
        camid[i] = -1;
        if (!registered[i])
            continue;
        bc[nsub] = cams[i];
        flags[nsub] = (i == seed1 || i == seed2)
                      ? (MV_BUNDLE_FIX_POSE | MV_BUNDLE_FIX_K)
                      : MV_BUNDLE_FIX_K;
        camid[i] = nsub++;
    }
    for (k = 0; k < ntracks; k++) {
        const mv_track *t = &tracks[k];
        int o;
        ptid[k] = -1;
        if (!xvalid[k])
            continue;
        ptid[k] = npts;
        memcpy(Xs + 3 * npts, X + 3 * k, 3 * sizeof(double));
        for (o = 0; o < t->len; o++) {
            const mv_feature *f;
            if (!registered[t->img[o]])
                continue;
            f = &feats[t->img[o]][t->feat[o]];
            obs[nobs].cam = camid[t->img[o]];
            obs[nobs].pt = npts;
            obs[nobs].u = f->u;
            obs[nobs].v = f->v;
            nobs++;
        }
        npts++;
    }
    if (npts < 1)
        return;
    if (mv_bundle_adjust(bc, flags, nsub, Xs, npts, obs, nobs, NULL)
        != MV_OK)
        return;
    for (i = 0; i < nimg; i++)
        if (camid[i] >= 0)
            cams[i] = bc[camid[i]];
    for (k = 0; k < ntracks; k++)
        if (ptid[k] >= 0)
            memcpy(X + 3 * k, Xs + 3 * ptid[k], 3 * sizeof(double));
}

int mv_sfm_register_incremental(mv_camera *cams, unsigned char *registered,
                                double *X, unsigned char *xvalid,
                                double *rms,
                                const mv_feature *const *feats,
                                const int *nfeat, int nimg,
                                const mv_track *tracks, int ntracks,
                                int seed1, int seed2)
{
    const mv_camera **tcams;
    double *tuv, *X3, *uv, *Xs;
    unsigned char *failed;
    mv_camera *bc;
    mv_bundle_obs *obs;
    unsigned *flags;
    int *camid, *ptid;
    long long totobs = 0;
    int i, k, nreg = 2;

    if (!cams || !registered || !X || !xvalid || !rms || !feats || !nfeat
        || nimg < 2 || ntracks < 0 || (ntracks > 0 && !tracks)
        || seed1 < 0 || seed1 >= nimg || seed2 < 0 || seed2 >= nimg
        || seed1 == seed2)
        return MV_ERR;
    for (i = 0; i < nimg; i++)
        if (!feats[i] || nfeat[i] < 0)
            return MV_ERR;
    /* the tracks must actually index into these images and features */
    for (k = 0; k < ntracks; k++) {
        int o;
        if (tracks[k].len < 2 || !tracks[k].img || !tracks[k].feat)
            return MV_ERR;
        for (o = 0; o < tracks[k].len; o++) {
            int im = tracks[k].img[o], f = tracks[k].feat[o];
            if (im < 0 || im >= nimg || f < 0 || f >= nfeat[im])
                return MV_ERR;
        }
        totobs += tracks[k].len;
    }

    tcams = (const mv_camera **)malloc((size_t)nimg
                                       * sizeof(const mv_camera *));
    tuv = (double *)malloc((size_t)nimg * 2 * sizeof(double));
    X3 = (double *)malloc((size_t)(ntracks > 0 ? ntracks : 1) * 3
                          * sizeof(double));
    uv = (double *)malloc((size_t)(ntracks > 0 ? ntracks : 1) * 2
                          * sizeof(double));
    Xs = (double *)malloc((size_t)(ntracks > 0 ? ntracks : 1) * 3
                          * sizeof(double));
    failed = (unsigned char *)calloc((size_t)nimg, 1);
    bc = (mv_camera *)malloc((size_t)nimg * sizeof(mv_camera));
    flags = (unsigned *)malloc((size_t)nimg * sizeof(unsigned));
    camid = (int *)malloc((size_t)nimg * sizeof(int));
    ptid = (int *)malloc((size_t)(ntracks > 0 ? ntracks : 1)
                         * sizeof(int));
    obs = (mv_bundle_obs *)malloc((size_t)(totobs > 0 ? totobs : 1)
                                  * sizeof(mv_bundle_obs));
    if (!tcams || !tuv || !X3 || !uv || !Xs || !failed || !bc || !flags
        || !camid || !ptid || !obs) {
        free(tcams);
        free(tuv);
        free(X3);
        free(uv);
        free(Xs);
        free(failed);
        free(bc);
        free(flags);
        free(camid);
        free(ptid);
        free(obs);
        return MV_ERR;
    }
    memset(registered, 0, (size_t)nimg);
    registered[seed1] = 1;
    registered[seed2] = 1;
    memset(X, 0, (size_t)(ntracks > 0 ? ntracks : 1) * 3
           * sizeof(double));
    if (ntracks > 0)
        memset(xvalid, 0, (size_t)ntracks);
    for (i = 0; i < nimg; i++)
        rms[i] = -1.0;

    sfm_triangulate_all(X, xvalid, cams, registered, feats, tracks,
                        ntracks, tcams, tuv);

    /* next-best-view loop: most visible triangulated tracks wins,
     * lowest index on ties (deterministic) */
    for (;;) {
        int best = -1, bestcnt = 0, n;
        mv_camera cam;
        for (i = 0; i < nimg; i++) {
            int c = 0;
            if (registered[i] || failed[i])
                continue;
            for (k = 0; k < ntracks; k++)
                if (xvalid[k] && sfm_obs_of(&tracks[k], i) >= 0)
                    c++;
            if (c > bestcnt) {
                bestcnt = c;
                best = i;
            }
        }
        if (best < 0 || bestcnt < 6)
            break;
        n = 0;
        for (k = 0; k < ntracks; k++) {
            int o;
            if (!xvalid[k])
                continue;
            o = sfm_obs_of(&tracks[k], best);
            if (o < 0)
                continue;
            memcpy(X3 + 3 * n, X + 3 * k, 3 * sizeof(double));
            uv[2 * n] = feats[best][tracks[k].feat[o]].u;
            uv[2 * n + 1] = feats[best][tracks[k].feat[o]].v;
            n++;
        }
        if (mv_resect_robust(&cam, X3, uv, &n) == MV_OK && n >= 6) {
            /* the resected pose is the initial value only; the known
             * intrinsics stay authoritative, and the fixed-K bundle
             * below repairs the focal-vs-depth slide the free-K DLT
             * pose carries */
            memcpy(cams[best].R, cam.R, sizeof(cam.R));
            memcpy(cams[best].t, cam.t, sizeof(cam.t));
            registered[best] = 1;
            nreg++;
            sfm_bundle(cams, registered, X, xvalid, feats, tracks,
                       ntracks, seed1, seed2, nimg, bc, flags, camid,
                       ptid, Xs, obs);
            sfm_triangulate_all(X, xvalid, cams, registered, feats,
                                tracks, ntracks, tcams, tuv);
        } else {
            failed[best] = 1; /* never retried: keeps the loop finite */
        }
    }
    /* final consolidation over the full registered set */
    if (nreg > 2)
        sfm_bundle(cams, registered, X, xvalid, feats, tracks, ntracks,
                   seed1, seed2, nimg, bc, flags, camid, ptid, Xs, obs);

    /* per-image reprojection RMS over the valid tracks it observes */
    for (i = 0; i < nimg; i++) {
        double ss = 0.0;
        int c = 0;
        if (!registered[i])
            continue;
        for (k = 0; k < ntracks; k++) {
            int o;
            double p[2], du, dv;
            if (!xvalid[k])
                continue;
            o = sfm_obs_of(&tracks[k], i);
            if (o < 0)
                continue;
            if (mv_cam_project(p, &cams[i], X + 3 * k) != MV_OK) {
                du = 1e6; /* behind the camera: keep it visible in rms */
                dv = 1e6;
            } else {
                du = p[0] - feats[i][tracks[k].feat[o]].u;
                dv = p[1] - feats[i][tracks[k].feat[o]].v;
            }
            ss += du * du + dv * dv;
            c++;
        }
        if (c > 0)
            rms[i] = sqrt(ss / c);
    }

    free(tcams);
    free(tuv);
    free(X3);
    free(uv);
    free(Xs);
    free(failed);
    free(bc);
    free(flags);
    free(camid);
    free(ptid);
    free(obs);
    return nreg;
}
