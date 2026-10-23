from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
from Bio.PDB import Atom, Chain, MMCIFIO, Model, Residue, Structure

from mn_molecular_design_studio.engine.workflow_engine.workflows import design


def _write_cif(path: Path, binder_bad: bool, include_binder: bool = True) -> None:
    structure = Structure.Structure("test")
    model = Model.Model(0)
    structure.add(model)
    chain_rows = [("A", [1])] + ([ ("B", [100, 102]) ] if include_binder else [])
    for chain_id, residue_ids in chain_rows:
        chain = Chain.Chain(chain_id)
        model.add(chain)
        for index, residue_id in enumerate(residue_ids):
            residue = Residue.Residue((" ", residue_id, " "), "GLY", " ")
            chain.add(residue)
            if chain_id == "A":
                c_coord = np.asarray([0.0, 5.0, 0.0])
                n_coord = np.asarray([1.329, 5.0, 0.0])
            elif index == 0:
                c_coord = np.asarray([0.0, 0.0, 0.0])
                n_coord = np.asarray([-1.0, 0.0, 0.0])
            else:
                n_coord = np.asarray([10.0 if binder_bad else 1.329, 0.0, 0.0])
                c_coord = np.asarray([11.0, 0.0, 0.0])
            residue.add(Atom.Atom("N", n_coord, 20.0, 1.0, " ", " N  ", residue_id * 2, element="N"))
            residue.add(Atom.Atom("C", c_coord, 20.0, 1.0, " ", " C  ", residue_id * 2 + 1, element="C"))
    writer = MMCIFIO()
    writer.set_structure(structure)
    writer.save(str(path))


def test_integrity_filter_keeps_source_chain_and_preserves_rejected_backbone(tmp_path):
    generated = tmp_path / "intermediate_designs"
    generated.mkdir()
    target = tmp_path / "target.cif"
    _write_cif(target, binder_bad=False, include_binder=False)
    # Source target chain A has a bad geometry and is intentionally ignored.
    _write_cif(generated / "design.cif", binder_bad=True)
    # BoltzGen remaps target chain letters on output; the token design mask
    # identifies the binder from residue order without comparing chain IDs.
    np.savez(generated / "design.npz", design_mask=np.asarray([0, 1, 1], dtype=np.float32))
    report_path = tmp_path / "mimosa_integrity_results.json"
    script = Path(__file__).parents[1] / "scripts/check_mimosa_backbone.py"
    subprocess.run([sys.executable, str(script), str(generated), str(target), str(report_path)], check=True)

    report = json.loads(report_path.read_text())
    assert report["candidate_count"] == 1
    assert report["rejected_count"] == 1
    assert report["candidates"][0]["binder_chains"] == ["B"]
    assert not (generated / "design.cif").exists()
    assert (tmp_path / "rejected_backbones/design.cif").exists()
    assert (tmp_path / "rejected_backbones/design.npz").exists()


def test_mimosa_generation_metrics_survive_candidate_normalization_and_stay_stage_scoped(tmp_path):
    from mn_molecular_design_studio.output_vocabulary import candidate_output

    raw = tmp_path / "artifacts/raw/boltzgen/run-vanilla/intermediate_designs"
    raw.mkdir(parents=True)
    target = tmp_path / "artifacts/inputs/target.cif"
    target.parent.mkdir(parents=True)
    _write_cif(target, binder_bad=False, include_binder=False)
    generated_cif = raw / "design.cif"
    _write_cif(generated_cif, binder_bad=False)
    diagnostics = {
        "algorithm_revision": "mimosa-si-v1",
        "requested_particles": 4,
        "effective_forward_batch_size": 4,
        "events": [{"step": 0, "energies": [2.5, 3.0, 4.0, 5.0], "weights": [0.6, 0.2, 0.1, 0.1]}],
    }
    np.savez(
        raw / "design.npz",
        mimosa_generated_metrics=json.dumps({
            "generated_motif_rmsd_angstrom": 1.25,
            "potential_half_sum_angstrom_squared": 3.125,
        }),
        mimosa_diagnostics=json.dumps(diagnostics),
        mimosa_fixed_aa=np.asarray([0, -1, 3]),
    )

    generated = design._normalize_boltzgen_generation_candidates(
        tmp_path,
        {"mimosa_enabled": True, "target_chains": ["A"], "hotspots": ""},
        target,
        write_output=False,
    )
    assert len(generated) == 1
    assert generated[0]["metrics"]["generated_motif_rmsd_angstrom"] == 1.25
    assert generated[0]["metrics"]["generated_motif_potential_half_sum_angstrom2"] == 3.125
    assert generated[0]["metrics"]["fixed_motif_residue_count"] == 2
    assert generated[0]["raw_metadata"]["mimosa_diagnostics"]["effective_forward_batch_size"] == 4
    normalized_generation = candidate_output(generated[0])["metrics"]
    assert normalized_generation["mimosa.generated_motif_rmsd"] == 1.25
    assert normalized_generation["mimosa.generated_potential"] == 3.125

    refolded = {
        **generated[0],
        "candidate_id": "refolded-1",
        "stage": "complex_refolding",
        "metrics": {
            "generated_motif_rmsd_angstrom": 99.0,
            "refolded_motif_rmsd_angstrom": 2.5,
        },
    }
    normalized_refold = candidate_output(refolded)["metrics"]
    assert "mimosa.generated_motif_rmsd" not in normalized_refold
    assert normalized_refold["mimosa.refolded_motif_rmsd"] == 2.5
