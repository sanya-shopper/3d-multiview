#!/bin/bash
# replay-ab.sh -- measure TIME TO CALIBRATE and TIME TO RIG SOLVE by
# replaying a recorded two-camera session into a built hub.
#
# The hub's own logs carry no wall-clock stamps, so this script
# timestamps hub stdout itself and reports, relative to hub start:
#   (a) first [cal] line per camera   -- time-to-calibrate
#   (b) first [ext] line              -- time-to-rig-solved (headline)
#   (c) decodes/min and counter-valid rate at that moment (throughput)
#   (d) fx / RMS per camera and the solved baseline |t| (quality)
# A run that does not solve inside the cap is reported as CENSORED --
# a "never solved" outcome is a result, not a reason to wait longer.
#
# Usage:
#   deploy/replay-ab.sh <stack-dir> <frames-dir> <feed-s> <cap-s> [outdir]
#     stack-dir   directory holding a built `hubengine` and `replaycam`
#                 (use the same commit for both: the camera side answers
#                 clock probes and old cameras do not)
#     frames-dir  recorded session with cam1_*.pgm and cam2_*.pgm
#     feed-s      inter-frame interval fed to each replaycam
#     cap-s       give up after this many seconds
#
# A/B: build two stacks (e.g. `git worktree add ... <baseline>` plus
# `make hubengine replaycam` in each) and run this against both at the
# SAME feed rate, at least 3 times each, and compare medians.
#
# CAVEAT, read before quoting any number from this script: replay is not
# live.  A recorded frame set holds only the frames the recording hub
# managed to decode, so feeding them faster than they were captured
# compresses the pattern's displayed counter timeline by the same
# factor.  Cross-camera anchor pairing is gated on that counter (see
# try_extrinsics in tools/livehub.c), so the feed rate changes pairing
# behaviour directly.  Feeding at the recorded capture cadence keeps the
# counter timeline honest but starves a fast decoder; feeding faster
# feeds the decoder but distorts pairing.  State which regime a number
# came from.

set -u
if [ $# -lt 4 ]; then
    sed -n '2,36p' "$0"
    exit 1
fi
STACK=$1
FRAMES=$2
FEED=$3
CAP=$4
OUT=${5:-replay-ab-$(date +%Y%m%d-%H%M%S)}
PORT=${MV_AB_PORT:-9931}

for b in hubengine replaycam; do
    if [ ! -x "$STACK/$b" ]; then
        echo "missing $STACK/$b (make hubengine replaycam)" >&2
        exit 1
    fi
done

mkdir -p "$OUT/logs"
F1=$(ls "$FRAMES"/cam1_*.pgm)
F2=$(ls "$FRAMES"/cam2_*.pgm)
if [ -z "$F1" ] || [ -z "$F2" ]; then
    echo "no cam1_*.pgm / cam2_*.pgm in $FRAMES" >&2
    exit 1
fi
# loop the frame list so the stream never runs dry inside the cap
NF=$(echo "$F1" | wc -l)
LOOPS=$(python3 -c "import math;print(max(1,math.ceil($CAP/($FEED*$NF))+1))")
L1=""; L2=""; i=0
while [ $i -lt "$LOOPS" ]; do L1="$L1 $F1"; L2="$L2 $F2"; i=$((i+1)); done

T0=$(python3 -c 'import time;print("%.3f"%time.time())')
# MV_MUTE=1: the hub speaks aloud otherwise, and speech forks cost time
MV_MUTE=1 "$STACK/hubengine" "$PORT" 0.1133 "$OUT/logs" 2>&1 \
  | python3 -u -c '
import sys, time
for l in sys.stdin:
    sys.stdout.write("%.3f %s" % (time.time(), l)); sys.stdout.flush()' \
  > "$OUT/hub.ts" &
sleep 1
HUB=$(pgrep -f "hubengine $PORT " | head -1)
"$STACK/replaycam" 127.0.0.1 "$PORT" 1 "$FEED" $L1 > "$OUT/cam1.log" 2>&1 &
R1=$!
"$STACK/replaycam" 127.0.0.1 "$PORT" 2 "$FEED" $L2 > "$OUT/cam2.log" 2>&1 &
R2=$!

END=$(python3 -c "import time;print('%.3f'%(time.time()+$CAP))")
OUTCOME=CENSORED
while :; do
    if grep -q '\[ext\]' "$OUT/hub.ts" 2>/dev/null; then
        OUTCOME=SOLVED
        break
    fi
    if [ "$(python3 -c "import time;print(1 if time.time()>$END else 0)")" \
         = "1" ]; then
        break
    fi
    sleep 0.25
done
kill "$R1" "$R2" 2>/dev/null
[ -n "$HUB" ] && kill -INT "$HUB" 2>/dev/null
sleep 2
pkill -f "hubengine $PORT " 2>/dev/null
sleep 0.5
cp "$OUT"/logs/*/summary.txt "$OUT/summary.txt" 2>/dev/null
cp "$OUT"/logs/*/hub.log "$OUT/hub.log" 2>/dev/null
# the hub records every decoded frame as a full-res PGM; keep the
# metrics, drop the multi-GB dump
rm -rf "$OUT/logs"

T0=$T0 OUTCOME=$OUTCOME FEED=$FEED python3 - "$OUT" <<'PY'
import os, re, sys
out = sys.argv[1]
t0 = float(os.environ["T0"])
cal, ext, status = {}, [], []
for line in open(os.path.join(out, "hub.ts")):
    try:
        ts, rest = line.split(" ", 1)
        ts = float(ts) - t0
    except ValueError:
        continue
    m = re.match(r"\[cal\] cam (\d+): fx ([\d.]+) .*RMS ([\d.]+) px", rest)
    if m:
        cal.setdefault(int(m.group(1)),
                       (ts, float(m.group(2)), float(m.group(3))))
    m = re.match(r"\[ext\] cam (\d+) <-> cam (\d+): baseline ([\d.]+) m"
                 r".*?(\d+) pairs, mean \|dt\| ([\d.]+) mm", rest)
    if m:
        ext.append((ts, float(m.group(3)), int(m.group(4)),
                    float(m.group(5))))
    m = re.match(r"\[cam (\d+)\] rx (\d+) \| decoded (\d+)/(\d+) \| "
                 r"ctr (\d+) \| views (\d+)", rest)
    if m:
        status.append((ts,) + tuple(int(g) for g in m.groups()))
cut = ext[0][0] if ext else (status[-1][0] if status else 0.0)
per = {}
for row in status:
    if row[0] <= cut + 0.01:
        per[row[1]] = row
print("feed %s s | outcome %s" % (os.environ["FEED"], os.environ["OUTCOME"]))
for cid in sorted(cal):
    print("  (a) cam %d first [cal] at %7.1f s  fx %.1f  RMS %.3f px"
          % (cid, cal[cid][0], cal[cid][1], cal[cid][2]))
for cid in sorted(set(cal) | set(per)):
    if cid not in cal:
        print("  (a) cam %d NOT calibrated" % cid)
if ext:
    print("  (b) first [ext] at %7.1f s  baseline %.3f m  %d pairs  "
          "mean |dt| %.1f mm" % (ext[0][0], ext[0][1], ext[0][2],
                                 ext[0][3]))
    print("      last  [ext] at %7.1f s  baseline %.3f m  %d pairs"
          % (ext[-1][0], ext[-1][1], ext[-1][2]))
else:
    print("  (b) NO [ext]: rig NOT solved within the cap")
for cid in sorted(per):
    t, _, rx, ok, dec, ctr, views = per[cid]
    print("  (c) cam %d at %.1f s: rx %d, decodes %d (%.1f/min), "
          "counter-valid %d/%d (%.2f), views %d"
          % (cid, t, rx, dec, 60.0 * dec / t if t else 0, ctr, dec,
             (ctr / dec) if dec else 0.0, views))
sp = os.path.join(out, "summary.txt")
if os.path.exists(sp):
    txt = open(sp).read()
    for m in re.finditer(r"CALIBRATED fx ([\d.]+) .*?RMS ([\d.]+) px", txt):
        print("  (d) fx %s RMS %s px" % (m.group(1), m.group(2)))
    m = re.search(r"baseline ([\d.]+) m", txt)
    if m:
        print("  (d) solved rig baseline |t| = %s m" % m.group(1))
PY
echo "artifacts in $OUT/ (hub.ts, hub.log, summary.txt)"
