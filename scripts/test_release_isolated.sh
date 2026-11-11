#!/usr/bin/env bash
set -euo pipefail
repository=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
test_root=${MDS_TEST_ROOT:-$(mktemp -d /tmp/mds-release-tests.XXXXXX)}
python=${MDS_TEST_PYTHON:-python}
mkdir -p "$test_root/studio" "$test_root/runs" "$test_root/reference" "$test_root/cluster"
extra_mounts=()
if [[ -d /data ]]; then
  extra_mounts+=(--tmpfs /data)
fi
exec bwrap --ro-bind / / --bind /tmp /tmp --bind "$repository" "$repository" \
  --dev /dev --unshare-net --tmpfs /run --tmpfs /mnt "${extra_mounts[@]}" --die-with-parent -- \
  env MN_MOLECULAR_DESIGN_STUDIO_HOME="$test_root/studio" \
  MN_MOLECULAR_DESIGN_STUDIO_APP_HOME="$test_root/studio" \
  MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR="$test_root/runs" \
  MN_MOLECULAR_DESIGN_STUDIO_REFERENCE_DIR="$test_root/reference" \
  MDS_CLUSTER_HOME="$test_root/cluster" \
  MN_MOLECULAR_DESIGN_STUDIO_AUTH_MODE=production \
  MN_MOLECULAR_DESIGN_STUDIO_SCHEDULER_MODE=local \
  MDS_GPU_TESTS=0 MDS_QC_DOCKER_TESTS=0 MDS_RFD4_DOCKER_TESTS=0 \
  PYTHONDONTWRITEBYTECODE=1 PYTHONPATH="$repository:$repository/tests" \
  "$python" -m pytest "$@"
