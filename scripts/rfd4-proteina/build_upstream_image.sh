#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
commit=6306a210dd17dcae8e4e0f381b5b74a719d9b30e
image="${1:-mn-rfd4-proteina:0.0.3-6306a21-cu130}"
build_dir="$(mktemp -d)"
trap 'rm -rf "$build_dir"' EXIT
source_dir="${RFD4_SOURCE:-$build_dir/source}"
if [[ -z "${RFD4_SOURCE:-}" ]]; then
    git init --quiet "$source_dir"
    git -C "$source_dir" fetch --quiet --depth 1 https://github.com/RosettaCommons/RFD4-Proteina.git "$commit"
fi
mkdir -p "$build_dir/context/upstream"
git -C "$source_dir" archive "$commit" pyproject.toml README.md LICENSE .project-root src configs data/benchmarks/monomer data/cath_label_mapping.pkl |
    tar -x -C "$build_dir/context/upstream"
docker build --build-arg "RFD4_COMMIT=$commit" --tag "$image" \
    --file "$script_dir/Dockerfile.upstream" "$build_dir/context"
