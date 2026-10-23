from mn_molecular_design_studio import output_vocabulary
from mn_molecular_design_studio.engine.workflow_engine.workflows import design
import json
import numpy as np


def test_boltzgen_final_candidate_uses_full_binder_sequence_from_complex(tmp_path, monkeypatch):
    complex_path = (
        tmp_path
        / "artifacts/raw/boltzgen/run-vanilla/final_ranked_designs/final_1_designs/rank1_target_binder.cif"
    )
    complex_path.parent.mkdir(parents=True)
    complex_path.write_text("data_test\n")
    target_path = tmp_path / "artifacts/raw/boltzgen/input/target.pdb"
    target_path.parent.mkdir(parents=True)
    target_path.write_text("END\n")

    row = {
        "id": "target_binder",
        "file_name": "target_binder.cif",
        "designed_sequence": "CDR_ONLY",
        "designed_chain_sequence": "FULL_VHH_SEQUENCE",
        "design_to_target_iptm": "0.7",
    }
    monkeypatch.setattr(design, "_boltzgen_metric_rows", lambda _run, _budget: [row])
    monkeypatch.setattr(
        design, "_boltzgen_final_design_path", lambda _run, _budget, _row: complex_path
    )
    monkeypatch.setattr(design, "_boltzgen_design_pae_index", lambda _run: {})
    monkeypatch.setattr(design, "_boltzgen_design_pae_path", lambda *_args: None)
    monkeypatch.setattr(
        design,
        "_infer_generated_chain_roles",
        lambda *_args: (["B"], ["A"], "sequence_match"),
    )
    monkeypatch.setattr(
        design,
        "_sequences_by_chain",
        lambda _path: {"A": "TARGET_SEQUENCE", "B": "FULL_VHH_SEQUENCE"},
    )
    monkeypatch.setattr(design, "_normalize_boltzgen_generation_candidates", lambda *_a, **_k: [])
    monkeypatch.setattr(
        design,
        "_merge_candidate_pools",
        lambda **kwargs: kwargs["native_candidates"],
    )
    monkeypatch.setattr(design, "write_candidates", lambda _run, _tool, candidates: candidates)

    candidates = design._normalize_boltzgen_vanilla_candidates(
        tmp_path,
        {"budget": 1, "target_chains": ["X"], "hotspots": "X334"},
        target_path,
    )

    assert len(candidates) == 1
    assert candidates[0]["binder_sequence"] == "FULL_VHH_SEQUENCE"
    normalized = output_vocabulary.candidate_output(candidates[0])
    assert normalized["designed_sequence"] == "FULL_VHH_SEQUENCE"
    assert normalized["capabilities"]["sequence"] is True


def test_mimosa_sidecar_metrics_normalize_with_separate_scopes(tmp_path):
    sidecar = tmp_path / "generated.npz"
    np.savez(
        sidecar,
        mimosa_generated_metrics=np.asarray(json.dumps({
            "potential_half_sum_angstrom_squared": 5.0,
            "generated_motif_rmsd_angstrom": 1.25,
        })),
        mimosa_fixed_aa=np.asarray([0, -1, 3]),
        mimosa_diagnostics=np.asarray(json.dumps({"num_particles_per_output": 4})),
    )
    metrics, diagnostics = design._mimosa_generation_metadata(sidecar)
    assert metrics["fixed_motif_residue_count"] == 2
    assert diagnostics["num_particles_per_output"] == 4
    candidate = {
        "candidate_id": "mimosa_1", "source_tool": "boltzgen", "stage": "generation.backbone_sequence",
        "binder_sequence": "AAAA", "metrics": metrics,
    }
    normalized = output_vocabulary.candidate_output(candidate)
    assert normalized["metrics"]["mimosa.generated_motif_rmsd"] == 1.25
    assert normalized["metrics"]["mimosa.generated_potential"] == 5.0
    assert normalized["metrics"]["mimosa.fixed_motif_residues"] == 2
    candidate["metrics"]["refolded_motif_rmsd_angstrom"] = 0.5
    candidate["stage"] = "complex_refolding"
    normalized_refold = output_vocabulary.candidate_output(candidate)
    assert normalized_refold["metrics"]["mimosa.refolded_motif_rmsd"] == 0.5
    assert "mimosa.generated_motif_rmsd" not in normalized_refold["metrics"]
