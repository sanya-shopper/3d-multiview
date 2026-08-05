# multiview coding standards

These rules are distilled from bugs this codebase actually shipped and
then had to fix. Each rule carries the defect it would have prevented,
because "bound allocations from input" is forgettable but "or you get
the 4.4 GB PGM bomb again" is not. Follow them in new code; when you
touch old code that violates one, fix it in passing.

## The gate is the toughest toolchain, and it is automatic

The development platform (macOS clang) is permissive; nearly every bug
below was invisible there and caught the instant a stricter gate ran.
CI is therefore a standing requirement, not optional:

- Build and test on the **deployment platform** (Ubuntu), not just dev.
  → the `clock_gettime` / `getaddrinfo` / `mkstemp` compile breaks that
  glibc hides under `-std=c99 -pedantic` and Darwin does not.
- `-Werror` across **gcc and clang at -O0/-O2/-O3**. → dead functions,
  misleading indentation, and `-Werror=stringop-overflow` at -O2 that
  only one compiler at one opt level emits.
- **ASan+UBSan, ThreadSanitizer, Valgrind memcheck** on the suites and
  on the live tools. → the hub use-after-free, the data races, leaks.
- **libFuzzer + cppcheck + scan-build** on the input surface. → the
  PGM decompression-bomb DoS.
- Run `cppcheck` locally before pushing; it is fast and installed.

## Input and memory safety

- **Never allocate a size derived from untrusted input without bounding
  it first.** → `mv_pgm_read` malloc'd `w*h` from the header; a crafted
  `P5 66350 66350` asked for 4.4 GB. Cap before the malloc.
- **Every `malloc`/`realloc` return is checked, and `realloc` goes
  through a temporary** so the old block survives failure. →
  `tag = realloc(tag, ...)` leaked the old buffer on OOM.
- **Long-running services use bounded/ring buffers, never
  append-and-stop.** → extrinsics silently froze once the anchor array
  filled a minute into a run.

## Numerical robustness

- **Gate values that flow to output with `isfinite`, and write range
  checks in positive form** (`!(x > 0)`, not `x <= 0`) so NaN is
  rejected rather than silently admitted (`NaN <= 0` is false). → NaN
  coordinates reaching a PLY.
- **Each numeric kernel exists once, in the library, and is tested at
  its degenerate cases.** → Rodrigues log collapsed to zero at θ≈π in
  one of three copies; the nearest-SO(3) projection was missing its
  determinant flip in one of four copies (`mv/rot.h` now holds the
  single tested implementation of both). If you are about to re-derive
  SVD-to-rotation, Rodrigues, pixel undistortion, or a qsort
  comparator, call the existing one.

## Concurrency

- **Every shared mutable field has one documented owner and a stated
  lock discipline. Slow work (decode, solve) runs on a private snapshot
  taken under the lock — never on a shared buffer.** → the hub's
  network thread freed the decoder's image buffer mid-decode (confirmed
  SEGV). Write the ownership rule as a comment on the shared struct.

## Output honesty and test fixtures

- **Any user-facing measured quantity gets a test that ties it to
  independent ground truth — and the synthetic fixture must be able to
  fail.** → `calibreal` printed optical-axis depth while labelling it
  the tape-measured distance; the sim masked it because the frame
  generator happened to put the target on the optical axis. A green
  test that cannot fail is worse than no test.

## House style

- Strict C99 + libm only; row-major doubles; `MV_OK`/`MV_ERR` returns.
- Declarations at block top; comments state a constraint the code
  cannot show, not a narration of the next line.
- Determinism is a feature: fixed LCG seeds, no wall-clock or RNG in
  library code.
- Tests are source: update them with the code, keep `make check` green
  before every commit.

## Tracked consolidation debt (do when touching these files)

Consolidated this round into `mv/rot.h` (Rodrigues, SO(3) projection —
both had shipped bugs) and `mv_cmp_double` (six copies). Still
duplicated, lower-risk, migrate opportunistically:

- **Iterative pixel undistortion** (`1 + r2*(k1 + r2*k2)` inverse) is
  copied in calibreal, rigcalib, livehub, slreal, and scenecloud's
  half-res variant. A public `mv_cam_undistort_points` should absorb
  them; `scenecloud`'s is fused with downsampling and stays bespoke.
- **Feature-test macros** (`_POSIX_C_SOURCE` / `_DARWIN_C_SOURCE`) are
  repeated atop six POSIX-using tools. Correct as-is (the per-file
  combo is the robust form), but a shared `tools/compat.h` would
  centralize the policy.
