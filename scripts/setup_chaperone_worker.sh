#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
studio_root="$(cd -- "$script_dir/.." && pwd)"
lock_file="$studio_root/requirements/chaperone-worker-py312-linux.lock"
expected_revision="$(sed -n 's/^# chaperone-source-revision: //p' "$lock_file")"
chaperone_repo="${CHAPERONE_REPO:-$(cd -- "$studio_root/../tools-to-implement/chaperone" 2>/dev/null && pwd || true)}"
python_bin="${CHAPERONE_SETUP_PYTHON:-python3.12}"

if [[ -z "$chaperone_repo" || ! -d "$chaperone_repo/.git" ]]; then
  echo "Set CHAPERONE_REPO to the Chaperone checkout used by this Studio." >&2
  exit 2
fi
actual_revision="$(git -C "$chaperone_repo" rev-parse HEAD)"
if [[ "$actual_revision" != "$expected_revision" ]]; then
  echo "Chaperone revision mismatch: expected $expected_revision, found $actual_revision." >&2
  echo "Check out the locked source revision in a separate Chaperone worktree, then retry." >&2
  exit 2
fi
"$python_bin" -c 'import sys; assert sys.version_info[:2] == (3, 12), sys.version' || {
  echo "The Chaperone worker requires Python 3.12." >&2
  exit 2
}

app_home="$(PYTHONPATH="$studio_root" "$python_bin" -c 'from mn_molecular_design_studio.runtime import app_home; print(app_home())')"
worker_dir="${MDS_CHAPERONE_WORKER_DIR:-$app_home/chaperone-worker-${expected_revision:0:8}}"
worker_python="$worker_dir/bin/python"
mkdir -p "$(dirname -- "$worker_dir")"
if [[ ! -x "$worker_python" ]]; then
  "$python_bin" -m venv "$worker_dir"
fi
"$worker_python" -m pip install -r "$lock_file"
"$worker_python" -c 'from chaperone.data import PipelineDataHandler; from blocks.util_blocks.filter import FilterBlock; print("Chaperone handler and FilterBlock import successfully")' \
  2>/dev/null || PYTHONPATH="$chaperone_repo/src:$chaperone_repo/src/blocks" "$worker_python" -c \
  'from chaperone.data import PipelineDataHandler; from blocks.util_blocks.filter import FilterBlock; print("Chaperone handler and FilterBlock import successfully")'

config_path="$app_home/chaperone-worker.json"
mkdir -p "$app_home"
CONFIG_PATH="$config_path" WORKER_PYTHON="$worker_python" CHAPERONE_REPO="$chaperone_repo" \
CHAPERONE_REVISION="$actual_revision" "$worker_python" - <<'PY'
import json
import os
from pathlib import Path

repo = Path(os.environ["CHAPERONE_REPO"]).resolve()
target = Path(os.environ["CONFIG_PATH"])
payload = {
    # Keep the virtualenv path. Resolving its python symlink would point at
    # the base interpreter and bypass the locked worker site-packages.
    "python": str(Path(os.environ["WORKER_PYTHON"]).absolute()),
    "source": str(repo / "src"),
    "blocks": str(repo / "src" / "blocks"),
    "source_revision": os.environ["CHAPERONE_REVISION"],
    "python_version": "3.12",
}
temporary = target.with_suffix(target.suffix + ".tmp")
temporary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
temporary.replace(target)
print(f"Configured Chaperone worker at {target}")
PY
