"""Opt-in native-container regression tests for RFD4 input preparation.

Set MDS_RFD4_DOCKER_TESTS=1 to run tests requiring the pinned container image.
These checks use --prepare-only and never launch generation or inference.
"""
from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest

from mn_molecular_design_studio.engine.workflow_engine.core.docker_runner import build_docker_command
from mn_molecular_design_studio.engine.workflow_engine.core.manifests import load_manifest
from mn_molecular_design_studio import rfd4_proteina_input as inputs


REPO = Path(__file__).resolve().parents[1]
PDL1 = REPO / "mn_molecular_design_studio/assets/proteina_complexa/PD-L1.pdb"


@pytest.mark.skipif(os.environ.get("MDS_RFD4_DOCKER_TESTS") != "1", reason="requires opt-in RFD4 Docker image")
def test_nanobody_prepare_preserves_target_segment_defaults_and_framework_loops(tmp_path):
    parameters = {
        "task": "nanobody",
        "input_mode": "guided",
        "target_structure": str(PDL1),
        "target_chains": ["A"],
        "num_designs": 1,
        "batch_size": 1,
        "nsteps": 60,
        "length_min": 60,
        "length_max": 120,
        "options": {"nanobody_framework": "Vobarilizumab"},
    }
    prepared = inputs.prepare_request(parameters, run_dir=tmp_path)
    runtime_dir = prepared["root"] / "runtime"
    runtime_dir.mkdir(parents=True, exist_ok=True)
    runtime_source = Path(__file__).parents[1] / "mn_molecular_design_studio/engine/workflow_engine/workflows/rfd4_runtime.py"
    shutil.copy2(runtime_source, runtime_dir / runtime_source.name)

    check_script = runtime_dir / "assert_nanobody_prepare.py"
    check_script.write_text(
        '''import json, sys
from pathlib import Path
import numpy as np
import rfd4_runtime as runtime
from atomworks.io.utils.standard_annotations import S_SEGMIN, S_SEGMAX

request = json.loads(Path(sys.argv[1]).read_text())
runtime.run(request, prepare_only=True)
inventory = json.loads((Path(request["root"]) / "prepared_inputs.json").read_text())
assert len(inventory) == 1, inventory
item = inventory[0]
assert item["target_chains"] == ["A"], item
assert item["designed_chains"] == ["B"], item
aa = runtime.load(item["path"], annotated=True)
target = aa.chain_id == "A"
segmin = S_SEGMIN.annotation(aa, default="generate")
segmax = S_SEGMAX.annotation(aa, default="generate")
assert target.any(), "prepared structure lost target chain A"
assert not S_SEGMIN.mask(aa, default="generate")[target].any(), "target has active segment-min conditions"
assert not S_SEGMAX.mask(aa, default="generate")[target].any(), "target has active segment-max conditions"
assert np.all(segmin[target] == -1), "target S_SEGMIN did not retain native -1 sentinel"
assert np.all(segmax[target] == -1), "target S_SEGMAX did not retain native -1 sentinel"
framework = aa.chain_id == "B"
expsegmin = aa.get_annotation("annotation_expsegmin_1_residue")
assert np.any(expsegmin[framework] > 0), "framework CDR expandable segments were lost"
runtime._require_output_target_chains("nanobody", {"target_chains": ["A"]}, {"A": "ACD", "B": "EFG"}, "result.cif")
try:
    runtime._require_output_target_chains("nanobody", {"target_chains": ["A"]}, {"B": "EFG"}, "broken.cif")
except ValueError as exc:
    assert "missing expected target chain(s): A" in str(exc), str(exc)
else:
    raise AssertionError("PPI output integrity check accepted a missing target chain")
# Unindexed motif tasks can legitimately merge source chains, so they are exempt.
runtime._require_output_target_chains("motif_scaffolding", {"target_chains": ["A"]}, {"B": "EFG"}, "motif.cif")
'''
    )

    manifest = {**load_manifest("rfd4_proteina"), "entrypoint": "python", "gpu": False}
    args = ["/work/" + check_script.relative_to(tmp_path).as_posix(),
            "/work/" + prepared["request_path"].relative_to(tmp_path).as_posix()]
    command = build_docker_command(manifest, tmp_path, args, environment={"PYTHONUNBUFFERED": "1"})
    result = subprocess.run(command, capture_output=True, text=True, check=False, timeout=300)
    assert result.returncode == 0, (result.stdout + "\n" + result.stderr)[-6000:]


@pytest.mark.skipif(os.environ.get("MDS_RFD4_DOCKER_TESTS") != "1", reason="requires opt-in RFD4 Docker image")
def test_condition_spec_internal_group_uses_resolved_backbone_atoms(tmp_path):
    """A source CIF's unresolved hydrogens cannot belong to finite C_CRDGRP targets."""
    parameters = {
        "task": "custom_conditions",
        "input_mode": "condition_spec",
        "input_files": [str(PDL1)],
        "num_designs": 1,
        "batch_size": 1,
        "nsteps": 60,
        "options": {"filter_conditions": False},
        "condition_spec": """global_args:
  input: PD-L1.pdb
binder:
  contig: "30-30,/0,A1-20"
  conditions:
    C_SEQ:
      true:
        - select: "sel('A')"
    C_CRDGRP:
      1:
        - select: "sel('A/*/*/[N,CA,C,O]')"
""",
    }
    prepared = inputs.prepare_request(parameters, run_dir=tmp_path)
    runtime_dir = prepared["root"] / "runtime"
    runtime_dir.mkdir(parents=True, exist_ok=True)
    runtime_source = Path(__file__).parents[1] / "mn_molecular_design_studio/engine/workflow_engine/workflows/rfd4_runtime.py"
    shutil.copy2(runtime_source, runtime_dir / runtime_source.name)

    check_script = runtime_dir / "assert_condition_spec_prepare.py"
    check_script.write_text(
        '''import json, sys
from pathlib import Path
import numpy as np
import rfd4_runtime as runtime
from rfproteina.conditions import C_CRDGRP

request = json.loads(Path(sys.argv[1]).read_text())
runtime.run(request, prepare_only=True)
inventory = json.loads((Path(request["root"]) / "prepared_inputs.json").read_text())
assert len(inventory) == 1, inventory
item = inventory[0]
aa = runtime.load(item["path"], annotated=True)
groups = C_CRDGRP.annotation(aa, default="generate")
target = np.isin(aa.chain_id, item["target_chains"])
assert target.any(), "prepared condition spec lost its target chain"
assert np.any(groups[target] == 1), "resolved target backbone was not assigned coordinate group 1"
assert np.all(np.isfinite(aa.coord[groups == 1])), "coordinate group 1 contains unresolved atoms"
assert np.any((groups == 1) & np.isin(aa.atom_name, ["N", "CA", "C", "O"]))
assert not np.any((groups == 1) & np.char.startswith(aa.atom_name.astype(str), "H"))
'''
    )

    manifest = {**load_manifest("rfd4_proteina"), "entrypoint": "python", "gpu": False}
    args = ["/work/" + check_script.relative_to(tmp_path).as_posix(),
            "/work/" + prepared["request_path"].relative_to(tmp_path).as_posix()]
    command = build_docker_command(manifest, tmp_path, args, environment={"PYTHONUNBUFFERED": "1"})
    result = subprocess.run(command, capture_output=True, text=True, check=False, timeout=300)
    assert result.returncode == 0, (result.stdout + "\n" + result.stderr)[-6000:]
