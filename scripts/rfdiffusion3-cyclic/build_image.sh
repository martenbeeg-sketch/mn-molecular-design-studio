#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
commit=8bffbcb9bec7eb8c3be0813f87d33c98faa0f5d2
image="${1:-mn-foundry:rfd3-cyclic-8bffbcb-cu128}"
build_dir="$(mktemp -d)"
trap 'rm -rf "$build_dir"' EXIT
source_dir="${FOUNDRY_SOURCE:-$build_dir/source}"
if [[ -z "${FOUNDRY_SOURCE:-}" ]]; then
    git init --quiet "$source_dir"
    git -C "$source_dir" fetch --quiet --depth 1 https://github.com/RosettaCommons/foundry.git "$commit"
fi
mkdir "$build_dir/context"
git -C "$source_dir" archive "$commit" src/foundry models/rfd3/src/rfd3 models/rfd3/configs | tar -x -C "$build_dir/context"
docker build --build-arg "FOUNDRY_COMMIT=$commit" --tag "$image" \
    --file "$script_dir/Dockerfile" "$build_dir/context"
