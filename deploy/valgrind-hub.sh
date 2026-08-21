#!/bin/sh
# Memcheck the live hub on the streaming path: synthetic camera frames
# in -> decode -> clean SIGTERM shutdown, all under Valgrind. The GOAL
# is memory correctness of the network+decode+shutdown paths; a full
# live calibration is NOT required here (memcheck's ~30x slowdown makes
# decoding enough views impractical in CI -- the plain stress job and
# the loopback already prove calibration). Success = memcheck clean
# (--error-exitcode=1) AND the hub actually received/decoded frames AND
# it shut down cleanly.
cd "$(dirname "$0")/.." || exit 2
BUILD_TARGET_PREFIX="${BUILD_TARGET_PREFIX:-$(cd .. && pwd)}"
OUT="$BUILD_TARGET_PREFIX/_buildoutput/3d-multiview"
PORT=9931
FDIR=/tmp/mv_vg_frames
VLOG=/tmp/mv_vg_hub.log
"$OUT/genframes" "$FDIR" 24
valgrind --leak-check=full --show-leak-kinds=definite,indirect \
         --errors-for-leak-kinds=definite,indirect --error-exitcode=1 \
         "$OUT/hubengine" $PORT 0.1133 > "$VLOG" 2>&1 &
VG=$!
sleep 12   # valgrind startup
"$OUT/replaycam" 127.0.0.1 $PORT 1 0.4 $(ls "$FDIR"/f*.pgm) >/dev/null 2>&1
sleep 10   # let the hub drain and decode some frames under memcheck
kill -TERM $VG
wait $VG
RC=$?
echo "--- valgrind summary ---"
grep -E "ERROR SUMMARY|definitely lost|indirectly lost" "$VLOG" || true
if ! grep -qE "decoded [1-9]|rx [1-9]" "$VLOG"; then
    echo "VALGRIND-HUB FAIL: hub processed no frames"; cat "$VLOG"; exit 1
fi
if ! grep -q "shutting down" "$VLOG"; then
    echo "VALGRIND-HUB FAIL: no clean shutdown"; cat "$VLOG"; exit 1
fi
[ $RC -eq 0 ] && echo "VALGRIND-HUB PASS (memcheck clean, frames processed, clean exit)" \
  || { echo "VALGRIND-HUB FAIL: memcheck errors (rc=$RC)"; cat "$VLOG"; exit 1; }
