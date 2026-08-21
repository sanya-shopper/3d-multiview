/* Frame-arrival delta-t statistics -- see tools/dtstats.h for the
 * contract and ALIGNED_ASSESSMENT.md section 8 for why this exists.
 * Pure bookkeeping and arithmetic: no sockets, no threads, no libmv. */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "dtstats.h"

typedef struct {
    unsigned camid;
    unsigned long long n, cap;
    double *t;
} dt_cam;

struct dtstats {
    int ncams;
    dt_cam cam[DT_MAXCAMS];
};

dtstats *dt_new(void)
{
    return (dtstats *)calloc(1, sizeof(dtstats));
}

void dt_free(dtstats *s)
{
    int i;
    if (!s)
        return;
    for (i = 0; i < s->ncams; i++)
        free(s->cam[i].t);
    free(s);
}

int dt_add(dtstats *s, unsigned camid, unsigned long long seq, double t)
{
    dt_cam *c = NULL;
    int i;
    (void)seq;
    for (i = 0; i < s->ncams; i++)
        if (s->cam[i].camid == camid)
            c = &s->cam[i];
    if (!c) {
        if (s->ncams >= DT_MAXCAMS)
            return -1;
        c = &s->cam[s->ncams++];
        c->camid = camid;
    }
    if (c->n == c->cap) {
        unsigned long long ncap = c->cap ? 2 * c->cap : 1024;
        double *nt = (double *)realloc(c->t, (size_t)ncap * sizeof(double));
        if (!nt)
            return -1;
        c->t = nt;
        c->cap = ncap;
    }
    c->t[c->n++] = t;
    return 0;
}

/* percentile by nearest-rank on a sorted array (q in [0,1]) */
static double pctl(const double *v, unsigned long long n, double q)
{
    unsigned long long k;
    if (n == 0)
        return 0.0;
    k = (unsigned long long)ceil(q * (double)n);
    if (k > 0)
        k--;
    if (k >= n)
        k = n - 1;
    return v[k];
}

static int dblcmp(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

static void camsum(const dt_cam *c, dt_camsum *o)
{
    unsigned long long i, m = c->n - 1;
    double *gap = (double *)malloc((size_t)(m ? m : 1) * sizeof(double));
    o->camid = c->camid;
    o->nframes = c->n;
    o->t_first = c->t[0];
    o->t_last = c->t[c->n - 1];
    o->fps = (o->t_last > o->t_first)
                 ? (double)m / (o->t_last - o->t_first) : 0.0;
    o->period_med = 0.0;
    o->jitter_p90 = 0.0;
    if (gap && m > 0) {
        for (i = 0; i < m; i++)
            gap[i] = c->t[i + 1] - c->t[i];
        qsort(gap, (size_t)m, sizeof(double), dblcmp);
        o->period_med = pctl(gap, m, 0.5);
        for (i = 0; i < m; i++)
            gap[i] = fabs(gap[i] - o->period_med);
        qsort(gap, (size_t)m, sizeof(double), dblcmp);
        o->jitter_p90 = pctl(gap, m, 0.9);
    }
    free(gap);
}

int dt_analyze(const dtstats *s, dt_report *r)
{
    const dt_cam *a, *b, *sp, *dn;
    double *dt, *adt, sum = 0.0;
    unsigned long long i, j, n;
    int i0 = -1, i1 = -1, k, sign;

    memset(r, 0, sizeof(*r));
    r->ncams = s->ncams;
    if (s->ncams < 2)
        return -1;

    /* per-camera summaries in ascending camid order */
    {
        int order[DT_MAXCAMS], oi, oj;
        for (oi = 0; oi < s->ncams; oi++)
            order[oi] = oi;
        for (oi = 0; oi < s->ncams; oi++)
            for (oj = oi + 1; oj < s->ncams; oj++)
                if (s->cam[order[oj]].camid < s->cam[order[oi]].camid) {
                    int t = order[oi];
                    order[oi] = order[oj];
                    order[oj] = t;
                }
        for (oi = 0; oi < s->ncams; oi++) {
            if (s->cam[order[oi]].n < 2)
                return -1;
            camsum(&s->cam[order[oi]], &r->cam[oi]);
        }
    }

    /* the two busiest cameras; ties broken toward lower camid */
    for (k = 0; k < s->ncams; k++)
        if (i0 < 0 || s->cam[k].n > s->cam[i0].n)
            i0 = k;
    for (k = 0; k < s->ncams; k++)
        if (k != i0 && (i1 < 0 || s->cam[k].n > s->cam[i1].n))
            i1 = k;
    a = &s->cam[i0];
    b = &s->cam[i1];
    if (a->camid > b->camid) {
        const dt_cam *t = a;
        a = b;
        b = t;
    }
    if (a->n < 8 || b->n < 8)
        return -1;
    r->cam_a = a->camid;
    r->cam_b = b->camid;

    /* iterate the sparser stream, monotone nearest scan in the denser;
     * sign convention is always t(cam_a) - t(cam_b) */
    if (a->n <= b->n) {
        sp = a;
        dn = b;
        sign = 1;
    } else {
        sp = b;
        dn = a;
        sign = -1;
    }
    n = sp->n;
    dt = (double *)malloc((size_t)n * sizeof(double));
    adt = (double *)malloc((size_t)n * sizeof(double));
    if (!dt || !adt) {
        free(dt);
        free(adt);
        return -1;
    }
    j = 0;
    for (i = 0; i < n; i++) {
        double d;
        while (j + 1 < dn->n
               && fabs(dn->t[j + 1] - sp->t[i])
                  <= fabs(dn->t[j] - sp->t[i]))
            j++;
        d = (double)sign * (sp->t[i] - dn->t[j]);
        dt[i] = d;
        adt[i] = fabs(d);
        sum += d;
    }
    r->npairs = n;
    r->dt_mean = sum / (double)n;

    {
        double pa = 0.0, pb = 0.0;
        for (k = 0; k < s->ncams; k++) {
            if (r->cam[k].camid == r->cam_a)
                pa = r->cam[k].period_med;
            if (r->cam[k].camid == r->cam_b)
                pb = r->cam[k].period_med;
        }
        r->period = (pa < pb) ? pa : pb;
    }

    /* histogram over signed dt before sorting */
    if (r->period > 0.0) {
        double half = r->period / 2.0;
        for (i = 0; i < n; i++) {
            if (dt[i] < -half)
                r->hist_lo++;
            else if (dt[i] >= half)
                r->hist_hi++;
            else {
                int bk = (int)((dt[i] + half) / r->period
                               * (double)DT_HBUCKETS);
                if (bk < 0)
                    bk = 0;
                if (bk >= DT_HBUCKETS)
                    bk = DT_HBUCKETS - 1;
                r->hist[bk]++;
            }
        }
        {
            unsigned long long nh = 0, nq = 0;
            for (i = 0; i < n; i++) {
                if (adt[i] < half)
                    nh++;
                if (adt[i] <= r->period / 4.0)
                    nq++;
            }
            r->frac_half = (double)nh / (double)n;
            r->frac_quarter = (double)nq / (double)n;
        }
    }

    qsort(dt, (size_t)n, sizeof(double), dblcmp);
    qsort(adt, (size_t)n, sizeof(double), dblcmp);
    r->dt_median = pctl(dt, n, 0.5);
    r->adt_p50 = pctl(adt, n, 0.5);
    r->adt_p90 = pctl(adt, n, 0.9);
    r->adt_p99 = pctl(adt, n, 0.99);
    r->adt_max = adt[n - 1];
    free(dt);
    free(adt);
    return 0;
}

void dt_print(const dt_report *r, FILE *f)
{
    int k;
    unsigned long hmax = 1;
    double tq = r->period / 4.0, th = r->period / 2.0;

    fprintf(f, "== frame-arrival dt report (assessment section 8) ==\n");
    for (k = 0; k < r->ncams; k++) {
        const dt_camsum *c = &r->cam[k];
        fprintf(f, "cam %u: %llu frames in %.1f s -> %.2f fps, "
                   "median period %.2f ms, inter-arrival jitter p90 "
                   "%.2f ms\n",
                c->camid, c->nframes, c->t_last - c->t_first, c->fps,
                1e3 * c->period_med, 1e3 * c->jitter_p90);
    }
    fprintf(f, "pairing cam %u vs cam %u: %llu pairs, frame period T = "
               "%.2f ms\n",
            r->cam_a, r->cam_b, r->npairs, 1e3 * r->period);
    fprintf(f, "  signed dt (cam%u - cam%u): mean %+.2f ms, median "
               "%+.2f ms\n",
            r->cam_a, r->cam_b, 1e3 * r->dt_mean, 1e3 * r->dt_median);
    fprintf(f, "  |dt|: p50 %.2f ms  p90 %.2f ms  p99 %.2f ms  max "
               "%.2f ms\n",
            1e3 * r->adt_p50, 1e3 * r->adt_p90, 1e3 * r->adt_p99,
            1e3 * r->adt_max);
    fprintf(f, "  within T/4 (%.2f ms): %.1f%%   within T/2 (%.2f ms): "
               "%.1f%%\n",
            1e3 * tq, 100.0 * r->frac_quarter, 1e3 * th,
            100.0 * r->frac_half);
    fprintf(f, "  section-6 prediction for free-running cameras: |dt| "
               "uniform on [0,T/2],\n"
               "  mean T/4 = %.2f ms; DECISION: p90 well under T/2 and "
               ">=95%% within T/2\n"
               "  -> phase argument holds, pair-then-decode refactor "
               "justified.\n",
            1e3 * tq);
    for (k = 0; k < DT_HBUCKETS; k++)
        if (r->hist[k] > hmax)
            hmax = r->hist[k];
    fprintf(f, "  signed-dt histogram [-T/2,+T/2), %d buckets "
               "(tails: %lu low, %lu high):\n",
            DT_HBUCKETS, r->hist_lo, r->hist_hi);
    for (k = 0; k < DT_HBUCKETS; k++) {
        double lo = -r->period / 2.0
                    + r->period * (double)k / (double)DT_HBUCKETS;
        int bar = (int)(50.0 * (double)r->hist[k] / (double)hmax), c2;
        fprintf(f, "  %+7.2f ms |", 1e3 * lo);
        for (c2 = 0; c2 < bar; c2++)
            fputc('#', f);
        fprintf(f, " %lu\n", r->hist[k]);
    }
}
