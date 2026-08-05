#!/bin/sh
# build-linux.sh -- staged Linux build + test for the multiview repo.
#
# Runs inside the Ubuntu container built from deploy/Dockerfile, or
# directly on any Ubuntu/POSIX box from anywhere:
#
#     sh deploy/build-linux.sh
#
# Stages:
#   1. make libmv.a                      (core library)
#   2. each C tool target, individually  (calibreal scenecloud rigcalib
#                                         annotate slreal replaycam
#                                         livehub nettest)
#   3. stream_cam_v4l2                   (real V4L2 code now: must
#                                         '#error work item under
#                                         construction' on __linux__;
#                                         reported as XFAIL, not FAIL)
#   4. build the C test binaries         (test_mv .. test_mux)
#   5. run each test binary
#
# The python bib check (tests/check_bib.py) is deliberately skipped:
# it validates the LaTeX doc, which is not built in this image, and is
# not a Linux-portability concern.  Test binaries are run directly
# rather than via 'make check' for the same reason.
#
# Exit status: 0 = all required stages passed; 1 = at least one
# required stage failed (details above the summary, analysis in
# deploy/LINUX-ISSUES.md).

set -u

# Run from the repo root regardless of invocation directory.
cd "$(dirname "$0")/.." || exit 2

LOG="${TMPDIR:-/tmp}/mv_stage_$$.log"
trap 'rm -f "$LOG"' EXIT

pass=0
fail=0
xfail=0
failed_names=""

# stage <name> <command...> : run a required stage.
stage() {
    name="$1"; shift
    printf '=== %-28s ' "$name"
    if "$@" >"$LOG" 2>&1; then
        echo "PASS"
        pass=$((pass + 1))
    else
        echo "FAIL"
        fail=$((fail + 1))
        failed_names="$failed_names $name"
        sed 's/^/    | /' "$LOG"
    fi
}

# stage_xfail <name> <command...> : a stage that is expected to fail on
# Linux (known stub).  Failure is reported but does not fail the run;
# an unexpected pass is also fine and reported as PASS.
stage_xfail() {
    name="$1"; shift
    printf '=== %-28s ' "$name"
    if "$@" >"$LOG" 2>&1; then
        echo "PASS (stub built; unexpected on Linux)"
        pass=$((pass + 1))
    else
        echo "XFAIL (known stub: '#error work item under construction')"
        xfail=$((xfail + 1))
    fi
}

echo "--- multiview Linux validation ($(uname -s) $(uname -m), $(cc --version 2>/dev/null | head -n 1)) ---"

# Stage 1: core library.
stage "libmv.a" make libmv.a

# Stage 2: C tools, one target at a time so each failure is attributed.
for t in calibreal scenecloud rigcalib annotate slreal replaycam \
         livehub nettest; do
    stage "tool $t" make "$t"
done

# Stage 3: the V4L2 streamer is real code now -- a required build.
stage "tool stream_cam_v4l2" make stream_cam_v4l2

# Stage 4: build test binaries.
for t in test_mv test_refine test_optimal test_feat test_session \
         test_photo test_bundle test_mux; do
    stage "build $t" make "$t"
done

# Stage 5: run test binaries (only those that were built).
for t in test_mv test_refine test_optimal test_feat test_session \
         test_photo test_bundle test_mux; do
    if [ -x "./$t" ]; then
        stage "run $t" "./$t"
    else
        printf '=== %-28s SKIP (binary missing: build stage failed)\n' \
               "run $t"
    fi
done

echo "--- summary: $pass passed, $fail failed, $xfail expected-fail (stub) ---"
if [ "$fail" -ne 0 ]; then
    echo "FAILED stages:$failed_names"
    echo "See deploy/LINUX-ISSUES.md for root-cause analysis and fixes."
    exit 1
fi
echo "RESULT: PASS (Linux build and tests clean; bib/doc checks skipped by design)"
exit 0
