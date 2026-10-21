#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
image="${1:-mn-rfd4-proteina:0.0.3-cu130}"

docker build \
  --build-arg RFD4_COMMIT=9d35ee4f6328fd9b17f5b84756ee5db99d6d89ac \
  --tag "$image" \
  --file "$repo_root/scripts/rfd4-proteina/Dockerfile" \
  "$repo_root"

printf 'Built %s (CUDA 13; shared by RTX 4090 and RTX 5090 hosts).\n' "$image"
