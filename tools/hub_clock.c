/* Camera clock synchronization over the existing TCP connection --
 * v0 STUBS: identity mapping, no probes sent, replies ignored, so the
 * hub behaves exactly as before this module is implemented.
 * OWNERSHIP (parallel build): clock-sync work item ONLY. */

#include "hub_clock.h"

void hub_clock_probe(int sock)
{
    (void)sock;
}

void hub_clock_on_msg(const unsigned char *body20, double t_rx)
{
    (void)body20;
    (void)t_rx;
}

double hub_clock_map(unsigned camid, double t)
{
    (void)camid;
    return t;
}

double hub_clock_err(unsigned camid)
{
    (void)camid;
    return -1.0;
}
