from __future__ import annotations

import json
from pathlib import Path

from mn_molecular_design_studio.target_identity import protein_sequence_reference, prediction_input_systems


def example_chains():
    path = Path(__file__).parents[1] / "mn_molecular_design_studio/components/cofolding_editor/examples.json"
    return [card["value"] for card in json.loads(path.read_text())["8JYR"]["cards"]]


def test_reference_match_preserves_chain_copies_and_does_not_guess_roles():
    sequences = example_chains()
    reference = protein_sequence_reference(list(reversed(sequences)))
    assert reference["pdb_id"] == "8JYR"
    assert reference["matched_chain_count"] == reference["reference_chain_count"] == 3
    assert protein_sequence_reference(sequences[:2])["matched_chain_count"] == 2
    assert not protein_sequence_reference([*sequences, sequences[0]])
    assert not protein_sequence_reference([sequences[0][1:], *sequences[1:]])


def test_historical_sequence_inputs_recover_reference_without_claiming_a_template():
    chains = [{"id": chr(65 + i), "sequence": sequence, "role": "unspecified"}
              for i, sequence in enumerate(example_chains())]
    data = {"inputs": {"tasks": [{"protein_chains": chains, "has_target_structure": False}]}}
    result = prediction_input_systems(data)
    assert len(result) == 1
    assert result[0]["label"] == "8JYR system"
    assert result[0]["reference_role"] == "sequence_source"
    assert result[0]["identity_kind"] == "system"
    assert result[0]["chain_count"] == 3
    assert "sequence" not in result[0]
    chains[0]["role"] = "target"
    assert prediction_input_systems(data) == []


def test_unmatched_protein_system_is_sequence_defined_and_mixed_system_is_not_mislabeled():
    task = {"protein_chains": [{"id": "A", "sequence": "AAAA", "role": "unspecified"}]}
    result = prediction_input_systems({"inputs": {"tasks": [task]}})[0]
    assert result["reference_type"] == "sequence_only"
    assert result["pdb_id"] == ""
    task["ligand_smiles"] = ["CCO"]
    assert prediction_input_systems({"inputs": {"tasks": [task]}}) == []


def test_editor_payload_retains_reference_only_for_unchanged_matching_sequences():
    from mn_molecular_design_studio.cofolding_input import model_to_payload

    model = {"cards": [{"id": i + 1, "kind": "protein", "value": sequence, "copies": 1}
                       for i, sequence in enumerate(example_chains())]}
    assert model_to_payload(model)["sequence_source_reference"]["pdb_id"] == "8JYR"
    model["cards"][0]["value"] = "AAA"
    assert "sequence_source_reference" not in model_to_payload(model)


def test_api_prediction_normalization_retains_reference_with_unspecified_chain_roles():
    from mn_molecular_design_studio.engine.workflow_engine.workflows.optimized_prediction import _prepare_prediction

    prepared = _prepare_prediction(engine="localfold", tasks=[{"protein_sequences": ":".join(example_chains())}],
        samples=1, seed=42, engine_options={"model": "mol1", "use_msa": False})
    task = prepared.normalized_tasks[0]
    assert task["sequence_source_reference"]["pdb_id"] == "8JYR"
    assert task["target_protein_ids"] == []
    assert task["binder_protein_ids"] == []
