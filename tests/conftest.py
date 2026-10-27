"""Small, isolated design artifacts for validation and API contract tests."""
import json
from pathlib import Path

import numpy as np
import pytest
from Bio.PDB import Atom, Chain, MMCIFIO, Model, Residue, Structure
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.SeqUtils import seq3


@pytest.fixture
def design_source_factory(tmp_path):
    def create(run_dir=None, *, status="completed", success=True, passed=True,
               project_id="test-project", missing_occupancy=True):
        run_dir = Path(run_dir or tmp_path / "design/source")
        run_dir.mkdir(parents=True)
        structure = Structure.Structure("design")
        model = Model.Model(0); structure.add(model)
        coordinates = {
            "A": ("ACDEFG", [(0, 0, 0), (2, 0, 0), (2, 3, 0), (3, 4, 1), (5, 4, 2), (5, 6, 2)]),
            "B": ("KLMNP", [(0, 8, 0), (2, 8, 0), (2, 9, 1), (4, 10, 1), (5, 10, 2)]),
        }
        for chain_id, (sequence, coords) in coordinates.items():
            chain = Chain.Chain(chain_id); model.add(chain)
            for i, (letter, coordinate) in enumerate(zip(sequence, coords), 1):
                residue = Residue.Residue((" ", i, " "), seq3(letter).upper(), " ")
                residue.add(Atom.Atom("CA", np.array(coordinate, dtype=float), 90., 1., " ", " CA ", i, element="C"))
                chain.add(residue)
        writer = MMCIFIO(); writer.set_structure(structure)
        path = run_dir / "designed.cif"; writer.save(str(path))
        if missing_occupancy:
            payload = MMCIF2Dict(str(path)); del payload["_atom_site.occupancy"]
            writer.set_dict(payload); writer.save(str(path))
        candidate = {
            "candidate_id": "passing-design", "source_tool": "bindcraft2", "stage": "complex_refolding",
            "complex_pdb": "designed.cif", "binder_sequence": "KLMNP",
            "target_chains": ["A"], "binder_chains": ["B"],
            "metrics": {"native_pass_filters": passed}, "raw_metadata": {},
        }
        normalized = run_dir / "artifacts/normalized_candidates"
        normalized.mkdir(parents=True)
        (normalized / "candidates.jsonl").write_text(json.dumps(candidate) + "\n")
        for name, payload in {
            "metadata.json": {"run_id": run_dir.name, "task_group": "design", "status": status,
                              "tool": "bindcraft2", "job_type": "bindcraft2_design"},
            "input.json": {"tool": "bindcraft2", "params": {}},
            "result.json": {"success": success},
            "molecular_design_studio.json": {"project_id": project_id, "studio_owner_id": "test-user"},
        }.items():
            (run_dir / name).write_text(json.dumps(payload))
        return run_dir, candidate
    return create


def pytest_collection_modifyitems(config, items):
    """GPU inference is opt-in even when the whole test suite is collected."""
    import os
    if os.environ.get("MDS_GPU_TESTS") == "1":
        return
    skip = pytest.mark.skip(reason="real GPU inference requires MDS_GPU_TESTS=1")
    for item in items:
        if item.get_closest_marker("gpu") or item.get_closest_marker("gpu_e2e"):
            item.add_marker(skip)


@pytest.fixture
def gpu_runtime(request, monkeypatch):
    """Use production workers/containers and retain real inference artifacts."""
    import datetime
    import os
    import subprocess
    from mn_molecular_design_studio.runtime import reference_root

    if os.environ.get("MDS_GPU_TESTS") != "1":
        pytest.skip("real GPU inference requires MDS_GPU_TESTS=1")
    device = os.environ.get("MDS_TEST_GPU", "0")
    try:
        inventory = subprocess.run(
            ["nvidia-smi", "--query-gpu=index,name,uuid,driver_version,memory.total", "--format=csv"],
            capture_output=True, text=True, check=True, timeout=30,
        ).stdout
        subprocess.run(["docker", "info"], capture_output=True, check=True, timeout=30)
    except (OSError, subprocess.SubprocessError) as exc:
        pytest.skip(f"GPU/container runtime unavailable: {exc}")
    reference = reference_root()
    root = Path(os.environ.get("MDS_GPU_ARTIFACTS", str(request.config.rootpath / ".tmp/gpu-tests")))
    root = root.resolve() / datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%S%fZ") / request.node.name
    root.mkdir(parents=True)
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_REFERENCE_DIR", str(reference))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR", str(root / "runs"))
    # Installed checkpoints only. A missing checkpoint is a skip, never a download.
    monkeypatch.setenv("HF_HUB_OFFLINE", "1")
    (root / "runtime.json").write_text(json.dumps({
        "gpu_inventory": inventory, "selected_device": device, "reference_root": str(reference),
        "test": request.node.nodeid,
    }, indent=2))

    def require_image(image):
        result = subprocess.run(["docker", "image", "inspect", image], capture_output=True, text=True, timeout=30)
        if result.returncode:
            pytest.skip(f"required local image unavailable: {image}")
        (root / (image.replace("/", "_").replace(":", "_") + ".json")).write_text(result.stdout)

    return {"root": root, "reference": reference, "device": device, "require_image": require_image}
