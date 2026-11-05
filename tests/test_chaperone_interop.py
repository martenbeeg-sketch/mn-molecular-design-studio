from __future__ import annotations

import copy
import hashlib
from pathlib import Path
import shutil

import pytest

pytest.importorskip("chaperone")

from mn_molecular_design_studio.chaperone_adapter import (
    filter_candidates,
    roundtrip_candidates,
    to_pipeline_data_handler,
)


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
