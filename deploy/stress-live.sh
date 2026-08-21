#!/bin/sh
# Live-stack stress: fuzz the hub, then stream synthetic pattern frames
# through it and (optionally) demand a live calibration. Reusable
# locally and in CI, plain or under sanitizers.
#   sh deploy/stress-live.sh [require_cal 0|1]  (default 1)
# Expects hubengine replaycam nettest genframes already built into the
# _buildoutput tree (CLAUDE.md T2) -- located below the same way the
# Makefile computes $(OUT).
# The replay uses camid 1 DELIBERATELY: nettest's baseline scenario
# occupies that slot, so this also proves a real camera can reuse its
# slot immediately after the fuzz battery.
cd "$(dirname "$0")/.." || exit 2
BUILD_TARGET_PREFIX="${BUILD_TARGET_PREFIX:-$(cd .. && pwd)}"
OUT="$BUILD_TARGET_PREFIX/_buildoutput/3d-multiview"
REQUIRE_CAL=${1:-1}
PORT=9911
FDIR=/tmp/mv_stress_frames
HLOG=/tmp/mv_stress_hub.log
NLOG=/tmp/mv_stress_nettest.log

"$OUT/genframes" "$FDIR" 60 || exit 1
"$OUT/hubengine" $PORT 0.1133 > "$HLOG" 2>&1 &
HUB=$!
sleep 2

echo "--- fuzz battery ---"
"$OUT/nettest" 127.0.0.1 $PORT --fast-drip 300 > "$NLOG" 2>&1
NS=$?
grep "BOTTOM LINE" "$NLOG" || true
if [ $NS -ne 0 ]; then
    echo "STRESS FAIL: hub dead after fuzz battery"
    cat "$NLOG"
    echo "--- hub log (sanitizer reports live here) ---"
    cat "$HLOG"
    kill -TERM $HUB 2>/dev/null
    exit 1
fi

echo "--- live calibration stream (3 passes over 60 frames) ---"
FR=$(ls "$FDIR"/f*.pgm)
"$OUT/replaycam" 127.0.0.1 $PORT 1 0.5 $FR $FR $FR > /dev/null 2>&1
sleep 5
kill -TERM $HUB 2>/dev/null
wait $HUB 2>/dev/null

echo "--- hub log tail ---"
tail -5 "$HLOG"
if [ "$REQUIRE_CAL" = "1" ]; then
    if grep -q "\[cal\] cam 1" "$HLOG"; then
        echo "STRESS PASS: hub survived fuzz and calibrated live"
    else
        echo "STRESS FAIL: no live calibration after fuzz"
        cat "$HLOG"
        exit 1
    fi
else
    if grep -q "decoded" "$HLOG"; then
        echo "STRESS PASS: hub survived fuzz and kept decoding"
    else
        echo "STRESS FAIL: hub produced no status output"
        cat "$HLOG"
        exit 1
    fi
fi
