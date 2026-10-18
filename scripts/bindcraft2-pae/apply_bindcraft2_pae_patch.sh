#!/usr/bin/env bash
set -euo pipefail

bindcraft_root="${1:-/opt/bindcraft}"
patch_file="${2:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/bindcraft2-pae.patch}"

if [[ ! -f "$bindcraft_root/bindcraft/af2.py" || ! -f "$bindcraft_root/bindcraft/MPNN_stage.py" || ! -f "$bindcraft_root/bindcraft/campaign.py" ]]; then
  echo "Expected BindCraft2 source files under: $bindcraft_root" >&2
  exit 2
fi

if grep -Fq "metrics['pae_raw'] = pae" "$bindcraft_root/bindcraft/af2.py"; then
  if grep -Fq 'def write_pae_matrix(' "$bindcraft_root/bindcraft/MPNN_stage.py" \
    && grep -Fq 'write_pae_matrix' "$bindcraft_root/bindcraft/campaign.py"; then
    echo "BindCraft2 PAE patch is already applied."
    exit 0
  fi
  echo "BindCraft2 source has a partial PAE patch; refusing to continue." >&2
  exit 3
fi

cd "$bindcraft_root"
patch --batch --forward --dry-run -p1 -i "$patch_file"
patch --batch --forward -p1 -i "$patch_file"
