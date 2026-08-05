#!/bin/sh
# Double-click me on the SECOND camera MacBook.
cd "$(dirname "$0")"
xattr -dr com.apple.quarantine . 2>/dev/null
printf 'Hub IP (shown on the hub Mac window): '
read HUBIP
[ -z "$HUBIP" ] && { echo "no IP entered"; exit 1; }
caffeinate -dimsu >/dev/null 2>&1 &
trap 'kill %1 2>/dev/null' INT TERM EXIT
exec ./stream_cam "$HUBIP" 9900 2 5
