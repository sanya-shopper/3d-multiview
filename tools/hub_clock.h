#ifndef HUB_CLOCK_H
#define HUB_CLOCK_H

/* OWNERSHIP (parallel build): this header, tools/hub_clock.c, and the
 * camera-side halves in tools/stream_cam.swift, tools/mvhub.swift, and
 * tools/stream_cam_v4l2.c belong to the CLOCK-SYNC work item ONLY.
 * livehub.c is FROZEN during the parallel phase; the hooks below are
 * already wired there (probe on the 1 Hz ACK tick, MVTS parse in the
 * network loop, hub_clock_map() at anchor-gating time) and v0 stubs
 * keep behavior identical until this module is implemented.
 *
 * WIRE FORMAT (fixed by this scaffold; little-endian):
 *   hub -> camera, ~1 Hz, piggybacked on the ACK tick:
 *     "MVPB" | f64 t_hub_send          (12 bytes)
 *   camera -> hub, immediate reply:
 *     "MVTS" | u32 camid | f64 t_hub_echo | f64 t_cam  (24 bytes)
 *   Cameras that do not implement MVPB simply never reply; the hub
 *   falls back to identity mapping (hub receive time), i.e. today's
 *   behavior.  Unknown 'M...' messages are skipped by the existing
 *   resync logic on both sides. */

/* hub -> camera probe; called with the camera's socket ~1/s */
void hub_clock_probe(int sock);

/* network thread hands over a complete 24-byte MVTS body (after the
 * "MVTS" magic); t_rx = hub CLOCK_MONOTONIC at receive */
void hub_clock_on_msg(const unsigned char *body20, double t_rx);

/* map a camera-reported capture timestamp onto the hub clock; camid
 * with no sync data returns t unchanged (identity fallback) */
double hub_clock_map(unsigned camid, double t);

/* offset quality estimate in seconds for status lines (<0 = no sync) */
double hub_clock_err(unsigned camid);

#endif /* HUB_CLOCK_H */
