#!/bin/sh
# Double-click on the HUB MacBook. One binary: hub + camera 1 (preview).
cd "$(dirname "$0")"
xattr -dr com.apple.quarantine . 2>/dev/null
exec ./mvhub
