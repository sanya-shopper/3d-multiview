#!/bin/sh
# Run livehub under Valgrind memcheck on the normal live path (synthetic
# camera stream -> live calibration -> clean SIGTERM shutdown). The fuzz
# battery is skipped here: memcheck's ~30x slowdown makes the flood
# scenario impractical, and ASan already covers the fuzz path.
set -e
PORT=9931
FDIR=/tmp/mv_vg_frames
VLOG=/tmp/mv_vg_hub.log
./genframes "$FDIR" 40
valgrind --leak-check=full --show-leak-kinds=definite,indirect \
         --errors-for-leak-kinds=definite,indirect --error-exitcode=1 \
         --track-fds=yes \
         ./livehub $PORT 0.1133 > "$VLOG" 2>&1 &
VG=$!
sleep 12   # valgrind startup is slow
./replaycam 127.0.0.1 $PORT 1 0.5 $(ls "$FDIR"/f*.pgm) >/dev/null 2>&1
sleep 8
kill -TERM $VG
wait $VG
RC=$?
echo "--- valgrind summary ---"
grep -E "ERROR SUMMARY|definitely lost|indirectly lost|FILE DESCRIPTORS|Open file" "$VLOG" || true
grep -q "\[cal\] cam 1" "$VLOG" && echo "hub calibrated live under valgrind" \
  || { echo "VALGRIND-HUB FAIL: no live calibration"; cat "$VLOG"; exit 1; }
[ $RC -eq 0 ] && echo "VALGRIND-HUB PASS (memcheck clean)" \
  || { echo "VALGRIND-HUB FAIL: memcheck errors (rc=$RC)"; cat "$VLOG"; exit 1; }
