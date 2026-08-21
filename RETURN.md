# RETURN.md — solo re-entry guide

Written 2026-08-21 for a return several weeks later, assuming only this
document, this repo, and the hardware.  Everything in section 3 was
re-verified green on 2026-08-21 on the dev Mac.

## 1. What this project is

`multiview` is a portable, dependency-free C99 library (plus tools, a
LaTeX paper in `doc/`, and two GitHub Pages web companions) for
multiview-geometry scene analysis with **two fixed cameras** watching an
overlapping volume.  The rig calibrates itself live: each camera MacBook
runs `stream_cam` and streams grayscale frames over TCP to a hub laptop
running `hubengine` (tools/livehub.c), which blind-decodes a Gray-code
pattern shown on a third, carried laptop, calibrates each camera (Zhang +
LM), and solves camera-to-camera extrinsics.  A live session has already
solved a real rig at a 0.319 m baseline.  The strategic document is
`ALIGNED_ASSESSMENT.md`: it decides (§4–§7) to refactor the hub from
decode-then-pair to **pair-then-decode**, and leaves ONE open measurement
(§8) that justifies or kills that refactor.  That measurement is the next
experiment (section 4 below).

## 2. Exact state as of 2026-08-21

What works, verified today:

- `make check` green locally: 17 C test suites, target-consistency,
  bibliography and web checks (519 assertions), plus the paper builds.
- `sh deploy/build-linux.sh` green locally (also under ASan+UBSan).
- The build was recently migrated to write ALL objects/binaries/pdf into
  `../_buildoutput/3d-multiview/` (CLAUDE.md T2).  The migration had
  broken every bare `make <tool>` name that CI and the deploy scripts
  use; fixed today (bare aliases in the Makefile, scripts locate
  binaries in the build tree, `tests/check_targets.py` re-tightened so
  this cannot silently regress again).
- New today: `dthub` (tools/dthub.c + tools/dtstats.c), the §8
  measurement harness, with its own suite `test_dtstats` in
  `make check`.  Verified end-to-end over loopback with `replaycam`.
- `hub_pair.c` (the pair collector) exists and is unit-tested, but is
  NOT yet integrated into livehub — that integration is exactly what §8
  decides.

Unverified / not touched for weeks:

- Anything requiring cameras: `stream_cam` on the camera MacBooks,
  live calibration, the recorded-session replay flows.  The
  `stream_cam` binaries previously copied to the camera MacBooks are
  from an older commit — rebuild and re-copy (see traps).
- GitHub CI: pushed green fixes today; if the badge is red when you
  return, start at the `build-consistency` job log.

## 3. Rebuild and test (verified 2026-08-21)

Everything builds into `../_buildoutput/3d-multiview/` (created on
demand; safe to delete wholesale).  `BUILD_TARGET_PREFIX` comes from
`~/.zshenv` and defaults to the repo's parent directory if unset.

    cd ~/Claude/Projects/3d-multiview
    make check          # builds lib + 17 suites + doc, runs everything
    sh deploy/build-linux.sh   # the CI script, runnable locally too

Individual binaries by bare name (they land in the build tree, not here):

    make hubengine replaycam genframes dthub
    OUT=../_buildoutput/3d-multiview     # where to find them

The paper alone: `make doc` (needs pdflatex+bibtex; pdf is copied back
to `doc/multiview.pdf`, which IS versioned — commit it when it changes).

## 4. THE next experiment (ALIGNED_ASSESSMENT.md §8)

**Question:** at full capture rate, before any decode, what is the
frame-arrival Δt distribution between the two cameras at the hub?
§6 predicts nearest-partner |Δt| uniform on [0, T/2] (mean T/4 ≈ 8 ms at
30 fps).  If measurement agrees, the pair-then-decode refactor
(integrating `hub_pair.c` into livehub) is justified; if not, the real
problem is capture-side and must be fixed first.

The harness is `dthub`: a hub that speaks the normal camera protocol,
stamps each frame's arrival on one monotonic clock, discards all pixels,
and prints the distribution plus the decision line.

### 4a. Rehearsal, no cameras needed (do this first)

    make dthub replaycam genframes
    OUT=../_buildoutput/3d-multiview
    "$OUT/genframes" /tmp/dtframes 40
    "$OUT/dthub" 9944 12 /tmp/dt.csv &
    "$OUT/replaycam" 127.0.0.1 9944 1 0.1 /tmp/dtframes/f*.pgm &
    "$OUT/replaycam" 127.0.0.1 9944 2 0.1 /tmp/dtframes/f*.pgm

Expected (reproduced today): both cams ~9.9 fps, period ~101 ms, |Δt|
p90 ≈ 2 ms, "within T/2: 100%", a histogram concentrated in the central
buckets, and `/tmp/dt.csv` with one `camid,seq,t_arrival_s` row per
frame.  If this works, the tooling is healthy.

### 4b. The real measurement (two camera MacBooks + hub)

1. On EACH camera MacBook, disable macOS camera Reactions/video effects
   (Control Center → Video Effects) — they throttle capture.
2. Rebuild and re-copy `stream_cam` (architecture-specific, and the
   camera side must match the hub's commit):
   `swiftc -O tools/stream_cam.swift -o stream_cam`
3. Hub laptop: `make dthub`, then
   `"$OUT/dthub" 9900 120 dt-$(date +%Y%m%d).csv`
   (120 s timed run; or omit the duration and Ctrl-C to stop).
4. Camera MacBooks (distinct camids, full rate — the fps argument is a
   throttle, so ask for 30):
   `./stream_cam <hub-ip> 9900 1 30`   and   `./stream_cam <hub-ip> 9900 2 30`
5. Point both cameras at anything; content is irrelevant (nothing is
   decoded).  Let the timed run finish.

### 4c. Reading the result and the decision rule

`dthub` prints per-camera achieved fps / median period T / inter-arrival
jitter, then pairing stats.  Read it in this order:

1. **Per-camera fps first.**  If either camera is far below its target
   rate (e.g. 12 fps when asked for 30), stop: the problem is
   capture-side (camera, effects overlay, or network bandwidth — see
   traps) and the Δt distribution is not yet meaningful.
2. **|Δt| percentiles.**  Decision rule, also printed by the tool:
   p90(|Δt|) well under T/2 AND ≥95% of pairs within T/2 →
   **§6 holds; proceed with the pair-then-decode refactor** (wire
   `hub_pair.c` into livehub: frame ring at ingest, pair-first decode
   scheduling — the integration points are listed in hub_pair.h).
3. **Failure signature.**  A flat or bimodal histogram, heavy tails, or
   `frac_half` well below 1 with healthy per-camera fps → arrival times
   at the hub do not preserve capture simultaneity; investigate
   transport (per-connection buffering, WiFi) before any refactor.

**Record:** append a dated addendum to `ALIGNED_ASSESSMENT.md` §8 with
the printed report (per-cam fps, T, p50/p90/p99, frac T/2, verdict), and
commit the CSV under `rec/` only if the run is the decisive one.
Regardless of outcome, §8's question is then closed — update it.

## 5. Known traps

- **Binaries live in `../_buildoutput/3d-multiview/`, not in the repo.**
  `./hubengine` in old notes means `$OUT/hubengine` now.  Deleting
  `_buildoutput` at any time is safe; the next make recreates it.
- **Stale camera-side binaries.**  `stream_cam` on the camera MacBooks
  does not update itself; rebuild with swiftc ON that machine (Intel vs
  Apple Silicon) after every hub-side protocol change.
- **macOS Reactions/video effects** silently cap or perturb camera
  output; disable in Control Center on every camera Mac before quoting
  any rate or Δt number.
- **Bandwidth at full rate:** 1280×720 gray at 30 fps is ~27.6 MB/s per
  camera, ~55 MB/s total — beyond most WiFi/hotspot links (the old
  deployment used an iPhone hotspot, 172.20.10.x).  Use wired Ethernet
  or a strong 5 GHz link; otherwise the measurement is of the network,
  not the cameras.  If bandwidth-limited, drop both cameras to the same
  lower fps and say so next to the recorded numbers.
- **`make check` needs TeX** (it rebuilds the paper's aux for the bib
  check) — on a machine without TeX use `sh deploy/build-linux.sh`,
  which skips doc/bib on purpose.
- **The hub talks; the harness doesn't.**  `hubengine` speaks status
  aloud (`MV_MUTE=1` silences); `dthub` is silent by design and sends
  the cameras nothing — both camera programs tolerate that.
- **Replay feed-rate caveat** (deploy/replay-ab.sh header): replayed
  sessions compress the pattern-counter timeline; never quote replay
  pairing numbers as if they were live.
