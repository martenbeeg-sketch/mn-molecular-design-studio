"""Controller-level smoke tests for active-learning stop conditions.

The docking and model steps are replaced with deterministic fixtures so these
tests exercise the real campaign loop and persisted stop decision without
launching UniDock, GPUs, or remote workers.
"""
from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from mn_molecular_design_studio import pipeline_graph
from mn_molecular_design_studio.engine.workflow_engine.workflows import active_learning_controller as controller


def _compound(key: int, library_id: str = "screening-library") -> dict:
    return {
        "library_id": library_id,
        "compound_key": key,
        "compound_id": f"compound-{key}",
        "compound_name": f"Compound {key}",
        "source_smiles": "CCO",
        "fingerprint": bytes(128),
        "selection_reason": "centroid_seed",
        "predicted_probability": None,
    }


def test_top_n_stability_uses_ranked_identity_and_configured_score_tolerance():
    previous = [
        {"library_id": "screen", "compound_key": 1, "compound_id": "one", "score": -9.0},
        {"library_id": "screen", "compound_key": 2, "compound_id": "two", "score": -8.0},
    ]
    ranked = [
        {"library_id": "screen", "compound_key": 1, "compound_id": "one", "best_score": -9.04},
        {"library_id": "screen", "compound_key": 2, "compound_id": "two", "best_score": -8.06},
    ]

    current, stable, maximum_change, same_order = controller._top_vector(ranked, 2, previous, 0.05)

    assert stable is False
    assert same_order is True
    assert maximum_change == pytest.approx(0.06)
    assert len(current) == 2

    ranked[1]["best_score"] = -8.04
    _, stable, maximum_change, same_order = controller._top_vector(ranked, 2, previous, 0.05)

    assert stable is True
    assert same_order is True
    assert maximum_change == pytest.approx(0.04)


@pytest.mark.parametrize(
    ("scenario", "expected_stop", "expected_docking_rounds", "expected_model_steps"),
    [
        ("stable", "stable_top_scores", 4, 3),
        ("iteration_limit", "maximum_iterations", 2, 1),
    ],
)
def test_campaign_persists_stop_reason_and_does_not_queue_another_round(
    monkeypatch, tmp_path, scenario, expected_stop, expected_docking_rounds, expected_model_steps,
):
    campaign = tmp_path / "campaign"
    campaign.mkdir()
    (campaign / "metadata.json").write_text(json.dumps({"status": "running"}))
    target_pdb = tmp_path / "receptor.pdb"
    target_pdb.write_text("fixture")
    feature_manifest = tmp_path / "features.json"
    feature_manifest.write_text("{}")
    centroid_manifest = tmp_path / "centroids.json"
    centroid_manifest.write_text("{}")

    study = pipeline_graph.new_graph("unidock_active_learning")
    nodes = {node["data"]["kind"]: node for node in study["nodes"]}
    nodes["active_library"]["data"]["params"].update(
        library_id="screening-library", feature_manifest=str(feature_manifest),
    )
    nodes["active_seed"]["data"]["params"].update(
        library_id="centroid-library", centroid_manifest=str(centroid_manifest), sample_size=4,
    )
    nodes["active_target"]["data"]["params"].update(target_pdb=str(target_pdb))
    nodes["active_decision"]["data"]["params"].update(
        top_n=2, stable_rounds=2 if scenario == "stable" else 100,
        score_tolerance=0.05, max_iterations=0 if scenario == "stable" else 2,
    )

    docking_rounds = []
    model_steps = []
    completed = []

    def fake_dock_round(*, campaign, round_number, rows, **_kwargs):
        docking_rounds.append(round_number)
        if scenario == "stable":
            round_scores = {
                1: [-9.0, -8.0, -6.0, -5.0],
                2: [-10.0, -7.0],
                3: [-4.0, -3.0],
                4: [-2.0, -1.0],
            }[round_number]
        else:
            round_scores = {
                1: [-9.0, -8.0, -7.0, -6.0],
                2: [-12.0, -11.0],
            }[round_number]
        scored = []
        for row, score in zip(rows, round_scores, strict=True):
            scored.append({
                **row,
                "score": score,
                "pose_path": "",
                "worker_node_id": "local-worker",
                "worker_node_name": "Local fixture worker",
                "worker_host": "",
                "worker_run_dir": str(campaign / f"worker-{round_number}"),
                "worker_job_code": f"W{round_number:04d}",
                "worker_cpu_cores": 2,
                "worker_gpu_device": "0",
            })
        return scored, []

    def fake_model_step(*, campaign, round_number, **_kwargs):
        model_steps.append(round_number)
        model_run = campaign / f"model-{round_number}"
        model_run.mkdir()
        (model_run / "input.json").write_text(json.dumps({"params": {"cpu_cores": 2, "gpu_device": "0"}}))
        (model_run / "metadata.json").write_text(json.dumps({"job_code": f"M{round_number:04d}"}))
        next_start = 5 + 2 * (round_number - 1)
        next_rows = [_compound(next_start), _compound(next_start + 1)]
        selection = model_run / "selected.parquet"
        write_selection(selection, next_rows)
        return model_run, selection

    def json_default(value):
        if isinstance(value, bytes):
            return {"__bytes__": value.hex()}
        raise TypeError(f"Cannot serialize {type(value).__name__}")

    def json_hook(value):
        if set(value) == {"__bytes__"}:
            return bytes.fromhex(value["__bytes__"])
        return value

    def read_selection(path):
        return json.loads(Path(path).read_text(), object_hook=json_hook)

    def write_selection(path, rows, *_args, **_kwargs):
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        Path(path).write_text(json.dumps(rows, default=json_default))

    def collect_poses(campaign, _scores, _count, **_kwargs):
        manifest = campaign / "artifacts" / "active_learning" / "retained_poses.json"
        manifest.parent.mkdir(parents=True, exist_ok=True)
        manifest.write_text(json.dumps({"retained_count": 0, "poses": []}))
        return []

    monkeypatch.setattr(controller, "validate_graph", lambda graph, **_kwargs: {"graph": graph})
    monkeypatch.setattr(controller, "create_job", lambda *_args, **_kwargs: SimpleNamespace(run_dir=campaign))
    monkeypatch.setattr(controller, "update_status", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(controller, "finish_job", lambda *args, **kwargs: completed.append((args, kwargs)))
    monkeypatch.setattr(controller, "_library_sources", lambda *_args: {
        "feature_manifest": {"row_count": 100},
        "feature_manifest_path": feature_manifest,
        "centroid_manifest": centroid_manifest,
        "seed_sample": [_compound(key, "centroid-library") for key in range(1, 5)],
    })
    monkeypatch.setattr(controller, "_selected_worker_specs", lambda _node: [{
        "node_id": "local-worker", "node_name": "Local fixture worker", "local": True,
        "cpu_cores": 2, "gpu_device": "0",
    }])
    monkeypatch.setattr(controller, "_dock_round", fake_dock_round)
    monkeypatch.setattr(controller, "_queue_model_step", fake_model_step)
    monkeypatch.setattr(controller, "_read_selection_table", read_selection)
    monkeypatch.setattr(controller, "_write_selection_table", write_selection)
    monkeypatch.setattr(controller, "_collect_best_poses", collect_poses)
    monkeypatch.setattr(controller, "_cleanup_docking_artifacts", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(controller, "_save_score_table", lambda path, _rows: Path(path).write_text("fixture"))
    monkeypatch.setattr(controller, "_save_failed_compound_table", lambda path, _rows: Path(path).write_text("fixture"))
    from mn_molecular_design_studio import bulk_screening_results
    monkeypatch.setattr(bulk_screening_results, "publish_final_selection",
                        lambda *args, **kwargs: {"imported_count": 0, "requested_limit": 0, "group_id": None})

    from mn_molecular_design_studio import workspace
    from mn_molecular_design_studio.engine.workflow_engine.workflows import ligand_tools

    monkeypatch.setattr(workspace, "mark_studio_job", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(ligand_tools, "validate_target_pdb", lambda _path: None)

    result = controller.run_active_learning_pipeline(
        study, "Stop smoke", "test-owner", "test-project",
    )

    assert result == campaign
    assert docking_rounds == list(range(1, expected_docking_rounds + 1))
    assert model_steps == list(range(1, expected_model_steps + 1))
    state = json.loads((campaign / "artifacts/active_learning/state.json").read_text())
    history = json.loads((campaign / "artifacts/active_learning/round_history.json").read_text())
    graph_state = json.loads((campaign / "graph_execution.json").read_text())
    decision_id = nodes["active_decision"]["id"]
    assert state["stop_reason"] == history["stop_reason"] == expected_stop
    assert history["rounds"][-1]["stop_reason"] == expected_stop
    assert graph_state["nodes"][decision_id]["status"] == "completed"
    assert graph_state["nodes"][decision_id]["stop_reason"] == expected_stop
    assert completed
    assert completed[-1][0][2]["metrics"]["stop_reason"] == expected_stop

    if scenario == "stable":
        assert [row["top_n_stable"] for row in history["rounds"]] == [False, False, True, True]
        assert [row["stable_round_count"] for row in history["rounds"]] == [0, 0, 1, 2]
        assert history["rounds"][1]["same_top_order"] is False
