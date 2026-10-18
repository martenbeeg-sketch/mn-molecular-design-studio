#!/usr/bin/env bash
set -euo pipefail
kit_dir="${1:?Usage: build_images.sh /path/to/uplifting-biomolecular-modeling}"
recipe_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
docker build -f "$recipe_dir/ESM-IF1.Dockerfile" -t mn-esm-if1-opt:2.0.1-2b369911 "$kit_dir"
docker build -f "$recipe_dir/Caliby.Dockerfile" -t mn-caliby-opt:0.1-41d31560 "$kit_dir"
