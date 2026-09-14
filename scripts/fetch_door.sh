#!/usr/bin/env bash
# Fetch the DoOR 2.0 odorant-response tables (Muench & Galizia 2016) used by simulation/receptors.py.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p ref/door
BASE=https://raw.githubusercontent.com/ropensci/DoOR.data/master/data
for f in door_response_matrix.csv odor.csv; do
  if [ -s "ref/door/$f" ]; then echo "have $f"; continue; fi
  echo "downloading $f"
  curl -fL --progress-bar -o "ref/door/$f" "$BASE/$f"
done
