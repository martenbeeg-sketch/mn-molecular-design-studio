#!/usr/bin/env bash
set -euo pipefail

studio_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
patch_directory="$studio_root/scripts/bindcraft2-pae"
base_image="${BINDCRAFT2_BASE_IMAGE:-mn-bindcraft2:c0ee638-cu13}"
image_tag="${BINDCRAFT2_PAE_IMAGE:-mn-bindcraft2:c0ee638-cu13-pae-v1}"

docker image inspect "$base_image" >/dev/null
docker build \
  --build-arg "BASE_IMAGE=$base_image" \
  --tag "$image_tag" \
  --file "$patch_directory/Dockerfile" \
  "$patch_directory"

echo "Built $image_tag from $base_image"
