# Ubuntu validation kit

The owner's ultimate deployment target is Ubuntu; development happens on
macOS. This kit exists to surface Linux portability problems *now*,
before any Ubuntu machine exists, by building and testing the whole C
side of the repo inside an `ubuntu:24.04` container.

## Contents

| File | Role |
|---|---|
| `Dockerfile` | `ubuntu:24.04` + `build-essential` (+ `python3`), copies the repo in, hands off to the build script |
| `Dockerfile.dockerignore` | keeps `.git`, `refs/`, PDFs, and host build artifacts out of the image (BuildKit) |
| `build-linux.sh` | in-container (or bare-Ubuntu) staged build+test; per-stage PASS/FAIL; exits nonzero on any failure |
| `validate.sh` | host-side driver: finds docker/podman, builds the image, runs the container, prints the summary |
| `LINUX-ISSUES.md` | every Linux issue found, with root cause and suggested one-line fix |

## How to run (one command)

    sh deploy/validate.sh

Exit 0 = Ubuntu validation passed; exit 1 = it ran and failed (see
`deploy/LINUX-ISSUES.md`); exit 2 = no container runtime on this host
(the script then prints exact instructions for running
`deploy/build-linux.sh` on a real Ubuntu box instead).

## What it validates

Inside Ubuntu 24.04 (gcc 13, glibc 2.39), `build-linux.sh` runs:

1. `make libmv.a` — the core library;
2. every Linux-relevant C tool, one target at a time: `calibreal`,
   `scenecloud`, `rigcalib`, `annotate`, `slreal`, `replaycam`,
   `livehub`, `nettest`;
3. `stream_cam_v4l2` — a declared stub that `#error`s on `__linux__`;
   reported as XFAIL, never fails the run;
4. builds the eight C test binaries `test_mv test_refine test_optimal
   test_feat test_session test_photo test_bundle test_mux`;
5. runs each test binary directly (deliberately *not* `make check`: the
   python bib checker validates the LaTeX doc, which is out of scope
   in-container).

## Current measured status (2026-08-05, macOS host)

- **Container run: NOT YET EXECUTED — no runtime available.** This host
  has neither `docker` nor `podman` (nor colima/lima, OrbStack, or
  Docker Desktop). `validate.sh` was executed and correctly took the
  no-runtime path (exit 2) with fallback instructions. No Linux results
  are being claimed.
- **Script logic smoke-tested on macOS**: `sh deploy/build-linux.sh` was
  run on the host, where all stages pass (26 passed, 0 failed;
  `stream_cam_v4l2` builds its non-Linux stub there). This validates the
  staging/summary/exit-code machinery, not Linux portability.
- **Static audit predicts 2 hard compile failures on Ubuntu**, both the
  same root cause (no `_POSIX_C_SOURCE` under `-std=c99 -pedantic`, so
  glibc hides `struct timespec`/`clock_gettime`/`CLOCK_MONOTONIC` and
  `struct addrinfo`/`getaddrinfo`): `tools/livehub.c` and
  `tools/replaycam.c`. Library, tests, and the other six tools are pure
  ISO C99 and expected clean. Details and one-line fixes:
  [`LINUX-ISSUES.md`](LINUX-ISSUES.md).

Once a runtime is installed (`brew install colima docker` or
`brew install podman`, or Docker Desktop), re-run `sh deploy/validate.sh`
and update this status block plus `LINUX-ISSUES.md` with measured
output.
