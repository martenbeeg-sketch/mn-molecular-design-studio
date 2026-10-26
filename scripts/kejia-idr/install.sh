#!/usr/bin/env bash
set -euo pipefail
studio_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
reference_dir="${1:-$(cd -- "$studio_dir" && python -c 'from mn_molecular_design_studio.runtime import reference_root; print(reference_root())')}"
source_dir="$reference_dir/kejia-idr/source"
revision=7f65c18cc57b3e72d290c885c540276c7d9c6ab8
if [[ ! -e "$source_dir" ]]; then
    mkdir -p -- "$reference_dir/kejia-idr"
    git clone --filter=blob:none --sparse https://github.com/drhicks/Kejia_peptide_binders.git "$source_dir"
fi
[[ "$(git -C "$source_dir" remote get-url origin)" == https://github.com/drhicks/Kejia_peptide_binders.git ]] || { echo "Unexpected source repository at $source_dir" >&2; exit 1; }
[[ -z "$(git -C "$source_dir" status --porcelain)" ]] || { echo "Preserving modified source checkout at $source_dir" >&2; exit 1; }
git -C "$source_dir" sparse-checkout set templates threading
git -C "$source_dir" checkout --detach "$revision"
docker build -f "$studio_dir/scripts/kejia-idr/Dockerfile" -t mn-kejia-idr:7f65c18 "$studio_dir"
