#!/bin/sh
# Double-click me on the PATTERN MacBook, then click the page once.
cd "$(dirname "$0")"
xattr -dr com.apple.quarantine . 2>/dev/null
caffeinate -dimsu >/dev/null 2>&1 &
open pattern.html
