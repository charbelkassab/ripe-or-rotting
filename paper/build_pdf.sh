#!/usr/bin/env bash
# Render paper/print.html to paper/ripe_or_rotting.pdf with headless Chrome or Chromium.
set -euo pipefail
cd "$(dirname "$0")"
for c in "${CHROME:-}" "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome" google-chrome chromium chromium-browser; do
  if [ -n "$c" ] && command -v "$c" >/dev/null 2>&1 || [ -x "$c" ]; then CHROME_BIN="$c"; break; fi
done
: "${CHROME_BIN:?Chrome or Chromium not found; set CHROME=/path/to/chrome}"
"$CHROME_BIN" --headless=new --disable-gpu --no-pdf-header-footer --virtual-time-budget=10000 \
  --print-to-pdf="$PWD/ripe_or_rotting.pdf" "file://$PWD/print.html"
echo "wrote paper/ripe_or_rotting.pdf"
