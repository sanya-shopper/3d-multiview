#!/bin/sh
# validate.sh -- host-side driver for the Ubuntu validation kit.
#
# One command, from anywhere in the repo:
#
#     sh deploy/validate.sh
#
# Finds docker or podman, builds the Ubuntu 24.04 image from
# deploy/Dockerfile, runs deploy/build-linux.sh inside it, and prints a
# compact PASS/FAIL summary per stage.
#
# Exit status: 0 = Linux validation passed
#              1 = validation ran and FAILED (see deploy/LINUX-ISSUES.md)
#              2 = could not run (no container runtime on this host)

set -u

cd "$(dirname "$0")/.." || exit 2
IMAGE=multiview-linux-validate

# --- Stage 0: locate a container runtime -------------------------------
RUNTIME=""
if command -v docker >/dev/null 2>&1; then
    RUNTIME=docker
elif command -v podman >/dev/null 2>&1; then
    RUNTIME=podman
fi

if [ -z "$RUNTIME" ]; then
    echo "=== runtime check              FAIL (neither docker nor podman found)"
    echo ""
    echo "No container runtime is installed on this host, so the Ubuntu"
    echo "build cannot be exercised here.  To validate on a real Ubuntu"
    echo "box instead:"
    echo ""
    echo "  1. Copy the repo over, e.g.:"
    echo "       rsync -a --exclude .git --exclude refs \\"
    echo "           \"$(pwd)/\" ubuntu-host:multiview/"
    echo "  2. On the Ubuntu box:"
    echo "       sudo apt-get install -y build-essential"
    echo "       cd multiview && sh deploy/build-linux.sh"
    echo ""
    echo "  build-linux.sh prints per-stage PASS/FAIL and exits nonzero"
    echo "  on any failure.  Known/predicted Linux issues are analyzed"
    echo "  in deploy/LINUX-ISSUES.md."
    echo ""
    echo "Alternatively, install a runtime here (e.g. Docker Desktop,"
    echo "'brew install colima docker' or 'brew install podman') and"
    echo "re-run: sh deploy/validate.sh"
    exit 2
fi
echo "=== runtime check              PASS ($RUNTIME)"

# --- Stage 1: build the Ubuntu image -----------------------------------
BLOG="${TMPDIR:-/tmp}/mv_imgbuild_$$.log"
trap 'rm -f "$BLOG"' EXIT
if "$RUNTIME" build -f deploy/Dockerfile -t "$IMAGE" . >"$BLOG" 2>&1; then
    echo "=== image build (ubuntu:24.04) PASS"
else
    echo "=== image build (ubuntu:24.04) FAIL"
    tail -n 40 "$BLOG" | sed 's/^/    | /'
    echo "--- summary: image build failed before any compile stage ran ---"
    exit 1
fi

# --- Stage 2: run the in-container build + tests -----------------------
echo "--- container stages (deploy/build-linux.sh on Ubuntu 24.04) ---"
if "$RUNTIME" run --rm "$IMAGE"; then
    echo "--- summary: Ubuntu validation PASS ---"
    exit 0
else
    echo "--- summary: Ubuntu validation FAIL ---"
    echo "Per-stage detail above; root causes and suggested one-line"
    echo "fixes are in deploy/LINUX-ISSUES.md."
    exit 1
fi
