/* Temporal frame-pair collector -- see hub_pair.h for the contract and
 * ALIGNED_ASSESSMENT.md for the design record. */

#include "hub_pair.h"

#include <math.h>
#include <string.h>

void hub_pair_init(hub_pair *hp, double gate_s)
{
    int c, i;
    memset(hp, 0, sizeof(*hp));
    hp->gate = gate_s;
    for (c = 0; c < HP_MAXCAMS; c++)
        for (i = 0; i < HP_RING; i++)
            hp->ring[c][i].handle = -1;
}

int hub_pair_push(hub_pair *hp, int cam, int handle, double t)
{
    hp_slot *s;
    int evicted;
    if (cam < 0 || cam >= HP_MAXCAMS || handle < 0)
        return -1;
    s = &hp->ring[cam][hp->head[cam]];
    evicted = s->handle;
    s->handle = handle;
    s->t = t;
    s->consumed = 0;
    hp->head[cam] = (hp->head[cam] + 1) % HP_RING;
    hp->npush[cam]++;
    return evicted;
}

/* newest-first walk of one camera's ring */
static hp_slot *slot_at_age(hub_pair *hp, int cam, int age)
{
    int idx = (hp->head[cam] - 1 - age + 2 * HP_RING) % HP_RING;
    return &hp->ring[cam][idx];
}

int hub_pair_take(hub_pair *hp, int cama, int camb,
                  int *ha, int *hb, double *dt)
{
    int agea, ageb;
    if (cama < 0 || cama >= HP_MAXCAMS || camb < 0 || camb >= HP_MAXCAMS
        || cama == camb)
        return 0;
    for (agea = 0; agea < HP_RING; agea++) {
        hp_slot *sa = slot_at_age(hp, cama, agea);
        hp_slot *best = NULL;
        double bestd = 0.0;
        if (sa->handle < 0 || sa->consumed)
            continue;
        for (ageb = 0; ageb < HP_RING; ageb++) {
            hp_slot *sb = slot_at_age(hp, camb, ageb);
            double d;
            if (sb->handle < 0 || sb->consumed)
                continue;
            d = fabs(sb->t - sa->t);
            if (!best || d < bestd) {
                best = sb;
                bestd = d;
            }
        }
        if (best && bestd <= hp->gate) {
            *ha = sa->handle;
            *hb = best->handle;
            *dt = best->t - sa->t;
            sa->consumed = 1;
            best->consumed = 1;
            hp->npairs++;
            hp->sum_dt += bestd;
            return 1;
        }
    }
    return 0;
}

long hub_pair_formed(const hub_pair *hp)
{
    return hp->npairs;
}

double hub_pair_mean_dt(const hub_pair *hp)
{
    return hp->npairs ? hp->sum_dt / (double)hp->npairs : 0.0;
}
