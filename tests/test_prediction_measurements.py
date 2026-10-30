import json

import numpy as np
import pytest
from Bio.PDB import Atom, Chain, Model, PDBIO, Residue, Structure

from mn_molecular_design_studio.cofolding_results import enrich_candidate
from mn_molecular_design_studio.engine.workflow_engine.workflows.analysis import _analysis_source_candidates
from mn_molecular_design_studio.engine.workflow_engine.workflows.optimized_prediction import ENGINE_PROFILES
from mn_molecular_design_studio.output_vocabulary import PREFIX, candidate_output, vocabulary
from mn_molecular_design_studio.prediction_measurements import (
    PAE_MATRIX_ENGINES,
    interface_pae_from_matrix,
    prediction_engine_catalog,
)
from mn_molecular_design_studio.validation_acceptance import METRICS, apply_acceptance, available_metrics


def _two_chain_structure(path):
    structure = Structure.Structure("fixture")
    model = Model.Model(0)
    structure.add(model)
    serial = 1
    for chain_id, x in (("A", 0.0), ("B", 8.0)):
        chain = Chain.Chain(chain_id)
        model.add(chain)
        residue = Residue.Residue((" ", 1, " "), "ALA", " ")
        chain.add(residue)
        residue.add(Atom.Atom("CA", np.array([x, 0.0, 0.0]), 90.0, 1.0, " ", " CA ", serial, element="C"))
        serial += 1
    writer = PDBIO()
    writer.set_structure(structure)
    writer.save(str(path))


def test_engine_catalog_matches_installed_profiles_and_marks_conditional_outputs():
    catalog = prediction_engine_catalog()
    by_engine = {row["engine"]: row for row in catalog["engines"]}
    assert set(ENGINE_PROFILES) <= set(by_engine)
    assert {"bindcraft2", "rfd4_proteina", "af2_initial_guess"} <= set(by_engine)
    assert set(ENGINE_PROFILES) - {"esmfold"} <= PAE_MATRIX_ENGINES
    assert "esmfold" not in PAE_MATRIX_ENGINES
    assert catalog["normalization"]["normalized_pae"] == "raw PAE / 31 Å; a shared numeric scale, not a model calibration"
    assert by_engine["mol1"]["measurements"][3]["availability"] == "conditional saved per-model sidecar"
    esmfold_pae = next(row for row in by_engine["esmfold"]["measurements"] if row["metric"] == "PAE matrix")
    assert esmfold_pae["available"] is False
    assert "No PAE value is inferred" in esmfold_pae["note"]
    bc2_ipae = next(row for row in by_engine["bindcraft2"]["measurements"] if row["metric"] == "iPAE")
    assert bc2_ipae["unit"] == "fraction of 31 Å PAE scale"
    rfd4_ipae = next(row for row in by_engine["rfd4_proteina"]["measurements"] if row["metric"] == "iPAE")
    assert rfd4_ipae["unit"] == "Å"
    assert "chain of interest" in rfd4_ipae["note"]
    assert vocabulary()["prediction_measurements"]["version"] == catalog["version"]


def test_pae_directional_means_keep_raw_angstroms_and_31_scale_separate():
    result = interface_pae_from_matrix([[0, 3], [5, 0]], binder_indices=[1], target_indices=[0])
    assert result["prediction_ipae_binder_to_target"] == 5.0
    assert result["prediction_ipae_target_to_binder"] == 3.0
    assert result["prediction_ipae_bidirectional"] == 4.0
    assert result["prediction_ipae_binder_to_target_norm31"] == pytest.approx(5 / 31, abs=1e-6)
    assert result["prediction_ipae_bidirectional_norm31"] == pytest.approx(4 / 31, abs=1e-6)
    assert interface_pae_from_matrix([[0, float("nan")], [1, 0]], [1], [0]) is None


def test_saved_prediction_pae_is_available_to_acceptance_without_reprediction(tmp_path):
    run_dir = tmp_path / "refolding-validation" / "boltz-run"
    run_dir.mkdir(parents=True)
    structure_path = run_dir / "complex.pdb"
    _two_chain_structure(structure_path)
    (run_dir / "pae.json").write_text(json.dumps({"predicted_aligned_error": [[0, 3], [5, 0]]}))
    candidate = {
        "candidate_id": "saved-boltz-complex",
        "source_tool": "optimized_prediction_boltz2",
        "stage": "complex_refolding",
        "complex_pdb": "complex.pdb",
        "binder_sequence": "A",
        "binder_chains": ["B"],
        "target_chains": ["A"],
        "metrics": {"iptm": 0.76},
        "raw_metadata": {"pae_path": "pae.json"},
    }
    candidate_path = run_dir / "artifacts/normalized_candidates/candidates.jsonl"
    candidate_path.parent.mkdir(parents=True)
    candidate_path.write_text(json.dumps(candidate) + "\n")

    enriched = _analysis_source_candidates(run_dir, candidate_path)[0]
    output = candidate_output(enriched)
    assert output["metrics"]["confidence.validation.interface.pae_binder_to_target"] == 5.0
    assert output["metrics"]["confidence.validation.interface.pae_target_to_binder"] == 3.0
    assert output["metrics"]["confidence.validation.interface.pae_bidirectional"] == 4.0
    assert output["metrics"]["confidence.validation.interface.pae_bidirectional_norm31"] == pytest.approx(4 / 31, abs=1e-6)
    provenance = enriched["raw_metadata"]["prediction_pae_measurement"]
    assert provenance["engine"] == "boltz2"
    assert provenance["unit"] == "Å before normalization"
    assert provenance["normalization"].endswith("not model calibration")

    choices = dict(available_metrics("miniprotein", "structured", "complex", engine="boltz2"))
    assert "predicted_pae_binder_to_target" in choices
    assert "predicted_pae_bidirectional_norm31" in choices
    field = METRICS["predicted_pae_bidirectional_norm31"]["field"]
    accepted, summary = apply_acceptance([enriched], {
        "binder_format": "miniprotein", "target_type": "structured",
        "mandatory_filters": [{"field": field, "operator": "<=", "value": 0.14}],
    })
    decision = accepted[0]["raw_metadata"]["validation_acceptance"]
    assert summary["accepted"] == 1
    assert decision["measurements"][field]["value"] == pytest.approx(4 / 31, abs=1e-6)
    assert decision["measurements"][field]["status"] == "measured"
    assert PREFIX + "confidence.validation.interface.pae_bidirectional_norm31" == field
