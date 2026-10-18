#!/usr/bin/env bash
set -euo pipefail

studio_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
upstream_root="$studio_root/tools/rfantibody"
image_tag="${RFANTIBODY_IMAGE:-mn-rfantibody:3e5a5bb-cu118-v1}"

if [[ ! -f "$upstream_root/pyproject.toml" || ! -f "$upstream_root/uv.lock" ]]; then
  echo "Pinned RFantibody source is missing from $upstream_root" >&2
  exit 1
fi

docker build \
  --tag "$image_tag" \
  --file "$upstream_root/Dockerfile.studio" \
  "$upstream_root"

docker image inspect "$image_tag" --format 'Built {{.RepoTags}} · {{.Size}} bytes'
