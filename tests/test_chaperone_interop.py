from __future__ import annotations

import copy
import hashlib
from pathlib import Path
import shutil
import sys

import pytest

pytest.importorskip("chaperone")

from mn_molecular_design_studio.chaperone_adapter import (
    apply_filter_block,
    filter_candidates,
    roundtrip_candidates,
    to_pipeline_data_handler,
)
from mn_molecular_design_studio import candidate_sets
from mn_molecular_design_studio.engine.workflow_engine.workflows import graph_pipeline


def candidates(tmp_path):
    """Two outputs from one parent, kept distinct by model and artifact path."""
    source_structure = Path(__file__).parents[1] / "mn_molecular_design_studio/assets/vina/1iep.pdb"
    source_digest = hashlib.sha256(source_structure.read_bytes()).hexdigest()
    structure_paths = []
    for model in range(2):
        path = tmp_path / f"prediction-model-{model}.pdb"
        shutil.copyfile(source_structure, path)
        structure_paths.append(path)
    assert hashlib.sha256(source_structure.read_bytes()).hexdigest() == source_digest
    return [
        {
            "candidate_id": "variant-7-boltz-model-0",
            "source_tool": "boltz2",
            "stage": "complex_refolding",
            "complex_pdb": str(structure_paths[0]),
            "binder_sequence": "ACDEFGHIK",
            "target_chains": ["A"],
            "binder_chains": ["B"],
            "parents": ["variant-7"],
            "metrics": {"iptm": 0.82, "native_pass_filters": True},
            "raw_metadata": {
                "model_index": 0,
                "prediction_id": "boltz-run-model-0",
                "target_sequences": {"A": "GGHHII"},
                "custom_source_payload": {"kept": [1, 2]},
            },
        },
        {
            "candidate_id": "variant-7-boltz-model-1",
            "source_tool": "boltz2",
            "stage": "complex_refolding",
            "complex_pdb": str(structure_paths[1]),
            "binder_sequence": "ACDEFGHIK",
            "target_chains": ["A"],
            "binder_chains": ["B"],
            "parents": ["variant-7"],
            "metrics": {"native_pass_filters": True},
            "raw_metadata": {"model_index": 1, "prediction_id": "boltz-run-model-1",
                             "target_sequences": {"A": "GGHHII"}},
        },
    ]


def test_real_chaperone_handler_roundtrip_preserves_studio_payload_and_models(tmp_path):
    source = candidates(tmp_path)
    handler, sidecar = to_pipeline_data_handler(source)
    assert len(handler) == 2
    assert handler._df["design_id"].to_list() == [row["candidate_id"] for row in source]
    assert handler._df["parent_design_id"].to_list() == ["variant-7", "variant-7"]
    assert handler._df["sequences"].to_list() == [
        '{"A": "GGHHII", "B": "ACDEFGHIK"}',
        '{"A": "GGHHII", "B": "ACDEFGHIK"}',
    ]

    restored = roundtrip_candidates(source)
    assert [{k: v for k, v in row.items() if k != "raw_metadata"} for row in restored] == [
        {k: v for k, v in row.items() if k != "raw_metadata"} for row in source
    ]
    assert [row["raw_metadata"]["prediction_id"] for row in restored] == [
        "boltz-run-model-0", "boltz-run-model-1"
    ]
    assert [row["raw_metadata"]["chaperone_interop"]["model_id"] for row in restored] == [0, 1]
    assert [row["complex_pdb"] for row in restored] == [row["complex_pdb"] for row in source]
    returned_paths = [row["raw_metadata"]["chaperone_interop"]["chaperone_record"]["structures"][0]["path"]
                      for row in restored]
    assert returned_paths == [row["complex_pdb"] for row in source]
    assert len(set(returned_paths)) == 2 and all(Path(path).is_file() for path in returned_paths)
    assert restored[0]["raw_metadata"]["custom_source_payload"] == {"kept": [1, 2]}
    assert sidecar["variant-7-boltz-model-0"] == source[0]


def test_real_chaperone_filter_marks_missing_required_metric_not_passed_and_keeps_review(tmp_path):
    source = candidates(tmp_path)
    original = copy.deepcopy(source)
    passing, decisions = filter_candidates(source, "studio.iptm >= 0.7")
    assert [row["candidate_id"] for row in passing] == ["variant-7-boltz-model-0"]
    assert decisions == [
        {"candidate_id": "variant-7-boltz-model-0", "status": "passed", "filter_expression": "studio.iptm >= 0.7"},
        {"candidate_id": "variant-7-boltz-model-1", "status": "not_passed", "filter_expression": "studio.iptm >= 0.7"},
    ]
    retained = roundtrip_candidates(source, filter_expression="studio.iptm >= 0.7")
    assert len(retained) == 2
    assert [row["raw_metadata"]["chaperone_interop"]["filter"]["status"] for row in retained] == ["passed", "not_passed"]
    assert source == original
    assert all(row["raw_metadata"]["chaperone_interop"]["filter"]["authority"].startswith("auxiliary")
               for row in retained)


def test_registered_canvas_filter_can_use_native_chaperone_filterblock(tmp_path):
    source = candidates(tmp_path)
    receipt = apply_filter_block(source, {"metric": "iptm", "operator": ">=", "threshold": .7})
    assert receipt["status"] == "applied"
    assert receipt["implementation"] == "blocks.util_blocks.filter.FilterBlock"
    assert receipt["passed_candidate_ids"] == ["variant-7-boltz-model-0"]
    assert receipt["candidate_count"] == 2


def test_chaperone_handler_is_used_at_generic_prediction_tool_handoff(tmp_path, monkeypatch):
    from mn_molecular_design_studio import chaperone_adapter

    monkeypatch.setattr(chaperone_adapter, "worker_configuration", lambda: {
        "python": sys.executable, "source": "", "blocks": "", "source_revision": "test-worker",
    })
    original = candidates(tmp_path)
    target_path = original[0]["complex_pdb"]
    from mn_molecular_design_studio.design_validation import protein_chains, read_structure, sequence
    actual_binder_sequence = sequence(protein_chains(read_structure(Path(target_path)))["B"])
    for row in original:
        row["binder_sequence"] = actual_binder_sequence
    target = {"target_pdb": target_path, "target_chains": ["A"]}
    source = candidate_sets.write_collection(tmp_path, "prediction-source", original, target)

    handed_off = candidate_sets.chaperone_handoff(tmp_path, "boltz2", source)
    records = list(candidate_sets.iter_candidates(handed_off))
    assert [row["candidate_id"] for row in records] == [row["candidate_id"] for row in original]
    assert [row["parents"] for row in records] == [row["parents"] for row in original]
    assert [row["binder_sequence"] for row in records] == [row["binder_sequence"] for row in original]
    assert [row["complex_pdb"] for row in records] == [row["complex_pdb"] for row in original]
    assert [row["metrics"] for row in records] == [row["metrics"] for row in original]
    assert all(row["raw_metadata"]["chaperone_interop"]["artifact_kind"] == "structure" for row in records)
    assert all(row["raw_metadata"]["chaperone_interop"]["model_id"] == row["raw_metadata"]["model_index"]
               for row in records)

    request = graph_pipeline.node_request("boltz2", {
        "samples": 1, "num_recycles": 1, "seed": 17, "optimization_mode": "off",
        "target_use_msa": False,
    }, handed_off, "Chaperone handoff fixture")
    tasks = request["parameters"]["tasks"]
    assert [task["source_candidate_id"] for task in tasks] == [row["candidate_id"] for row in original]
    assert all(task["binder_sequence"] == row["binder_sequence"] for task, row in zip(tasks, original))
