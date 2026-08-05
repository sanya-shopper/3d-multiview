# Linux (Ubuntu 24.04) portability issues

Status of this file: **static-analysis findings, not yet container-measured.**
No container runtime exists on this host (checked: `docker`, `podman`,
colima/lima, OrbStack, Docker.app — none installed), so
`deploy/validate.sh` could not execute the Ubuntu build here. The issues
below come from a line-by-line audit of every source the Linux build
compiles, against glibc/gcc-13 semantics (Ubuntu 24.04's
`build-essential`). When the container run first happens, replace the
"predicted" wording with the measured compiler output.

The build uses `CFLAGS = -std=c99 -pedantic -Wall -Wextra -O2` with **no
feature-test macros anywhere in the tree** (verified:
`grep -rn '_POSIX_C_SOURCE\|_DEFAULT_SOURCE\|_GNU_SOURCE' src tools tests include`
finds nothing). On macOS this is harmless because Darwin headers expose
POSIX/BSD APIs regardless of `-std=`. On glibc, `-std=c99` defines
`__STRICT_ANSI__`, which stops `features.h` from auto-enabling
`_DEFAULT_SOURCE`, so `<time.h>` and `<netdb.h>` hide everything newer
than ISO C. That is the classic Mac-to-Linux trap, and it bites twice
here.

## Issue 1: `tools/livehub.c` — predicted hard compile FAIL

- `tools/livehub.c:94` — `struct timespec ts;` — type is not defined by
  glibc `<time.h>` under strict C99 (needs `_POSIX_C_SOURCE >= 199309L`).
  Predicted gcc error: `storage size of 'ts' isn't known`.
- `tools/livehub.c:95` — `clock_gettime(CLOCK_MONOTONIC, &ts);` —
  `CLOCK_MONOTONIC` undeclared (hard error) and `clock_gettime`
  implicitly declared (warning on gcc 13, error on gcc >= 14).

Root cause: missing feature-test macro, as above. Everything else in the
file survives strict mode (`select`/`fd_set`/`struct timeval` are
unconditional in `<sys/select.h>`, `mkdir` at line 426 is unconditional
in `<sys/stat.h>`, `SO_REUSEADDR`/`setsockopt`/`signal(SIGPIPE, ...)`
are fine, and `<netinet/tcp.h>` provides `TCP_NODELAY` on Linux).

Suggested one-line fix: add as line 1 of `tools/livehub.c`, before any
`#include`:

    #define _POSIX_C_SOURCE 200809L

## Issue 2: `tools/replaycam.c` — predicted hard compile FAIL

- `tools/replaycam.c:45` — `struct addrinfo hints, *res;` — in glibc
  `<netdb.h>`, `struct addrinfo` (and `getaddrinfo`/`freeaddrinfo`,
  lines 63 and 72) are guarded by `__USE_XOPEN2K`, which strict C99
  leaves unset. Predicted errors: `storage size of 'hints' isn't known`
  plus invalid use of undefined type on `res->ai_family` etc. (line 66).
- `tools/replaycam.c:79,109` — `struct timespec` — same hiding as
  Issue 1.
- `tools/replaycam.c:87` — `clock_gettime(CLOCK_MONOTONIC, ...)` —
  `CLOCK_MONOTONIC` undeclared (hard error).
- `tools/replaycam.c:112` — `nanosleep(&d, NULL);` — implicit
  declaration (warning on gcc 13, error on gcc >= 14).

Root cause: identical to Issue 1. `AF_INET`/`SOCK_STREAM`/`socket`/
`connect` come from `<sys/socket.h>` unconditionally and are fine.

Suggested one-line fix: add as line 1 of `tools/replaycam.c`, before any
`#include`:

    #define _POSIX_C_SOURCE 200809L

## Preferred single fix for both (main-thread decision)

Instead of two per-file defines, one Makefile line fixes both and
future-proofs new tools (`Makefile:3` area):

    CPPFLAGS += -D_POSIX_C_SOURCE=200809L

This is safe for the pure-ISO-C files (it only *adds* declarations) and
matches the `-std=c99 -pedantic` policy. Either form works; pick one.

## Issue 3 (informational): `tools/stream_cam_v4l2.c` — intentional stub

`tools/stream_cam_v4l2.c:3-4` is `#ifdef __linux__` /
`#error "work item under construction"`: the Linux build of this target
fails **by design** until the v4l2 work item lands.
`deploy/build-linux.sh` treats it as XFAIL (reported, does not fail the
run). No action needed from the portability side; noted so a future
container FAIL on this target is not misread as a regression.

## Audited and expected clean

- `src/*.c`, `tests/*.c`, `include/` — pure ISO C99; zero POSIX/BSD
  headers or calls (verified by grep for `unistd|sys/|netdb|arpa|
  netinet|signal|time.h` includes and for `strdup|fdopen|popen|u_int|
  qsort_r|usleep|drand48` uses). `libmv.a` and all eight test binaries
  should build and link cleanly on Ubuntu.
- Link order — every link line ends `libmv.a $(LDLIBS)` with
  `LDLIBS = -lm`; correct for GNU ld's left-to-right resolution.
- Byte order — wire helpers (`get32`/`put32` in the live tools) are
  explicit little-endian byte shuffles, and both platforms are LE; the
  `memcpy` of a double into 8 LE-serialized bytes in
  `tools/replaycam.c:94-96` is IEEE-754-safe on both.
- No `u_int*` types, no `qsort_r`, no `MSG_NOSIGNAL`/`SO_NOSIGPIPE`
  divergence (SIGPIPE handled via `signal(SIGPIPE, SIG_IGN)`), no
  `sys/errno.h`-style Darwinisms.
