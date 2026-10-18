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
