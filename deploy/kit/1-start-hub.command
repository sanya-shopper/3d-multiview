#!/bin/sh
# Double-click me on the HUB MacBook (also acts as camera 1).
cd "$(dirname "$0")"
xattr -dr com.apple.quarantine . 2>/dev/null
IP=$(ipconfig getifaddr en0 2>/dev/null || ipconfig getifaddr en1 2>/dev/null || echo "?")
printf '\n==================================================\n'
printf '  On the OTHER MacBook, 2-start-camera.command\n'
printf '  will ask for this HUB IP:   %s\n' "$IP"
printf '==================================================\n\n'
caffeinate -dimsu >/dev/null 2>&1 &        # keep this Mac awake
CAFF=$!
./stream_cam 127.0.0.1 9900 1 5 >/tmp/mvcam1.log 2>&1 &
CAM=$!
trap 'kill $CAM $CAFF 2>/dev/null' INT TERM EXIT
./livehub 9900 0.1133 rec                  # foreground; Ctrl-C to finish
