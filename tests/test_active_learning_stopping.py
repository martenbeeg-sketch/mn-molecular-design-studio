"""Controller-level smoke tests for active-learning stop conditions.

The docking and model steps are replaced with deterministic fixtures so these
tests exercise the real campaign loop and persisted stop decision without
launching UniDock, GPUs, or remote workers.
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path
from types import SimpleNamespace

import pytest
import pyarrow as pa
import pyarrow.parquet as pq

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


@pytest.mark.parametrize("saved_import_limit", [0, 1000])
@pytest.mark.parametrize(
    ("scenario", "expected_stop", "expected_docking_rounds", "expected_model_steps"),
    [
        ("stable", "stable_top_scores", 4, 3),
        ("iteration_limit", "maximum_iterations", 2, 1),
    ],
)
def test_campaign_persists_stop_reason_and_does_not_queue_another_round(
    monkeypatch, tmp_path, saved_import_limit, scenario, expected_stop, expected_docking_rounds, expected_model_steps,
):
    campaign = tmp_path / "campaign"
    campaign.mkdir()
    (campaign / "molecular_design_studio.json").write_text(json.dumps({
        "studio_owner_id": "test-owner", "project_id": "test-project",
    }))
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
    nodes["active_results"]["data"]["params"]["final_database_import_limit"] = saved_import_limit

    docking_rounds = []
    model_steps = []
    completed = []
    final_import_limits = []
    cleaned_worker_runs = []

    def fake_dock_round(*, campaign, round_number, rows, **kwargs):
        docking_rounds.append(round_number)
        assert kwargs["docking"]["batch_retained_pose_count"] == 0
        worker_run = campaign / f"worker-{round_number}"
        for artifact in ("poses", "ligands", "prepared"):
            artifact_dir = worker_run / "artifacts" / artifact
            artifact_dir.mkdir(parents=True, exist_ok=True)
            (artifact_dir / "fixture.dat").write_text(f"round-{round_number}")
        staging = kwargs["work_directory"] / "worker-00" / "batch-0001"
        staging.mkdir(parents=True, exist_ok=True)
        (staging / "staged-input.json").write_text(f"round-{round_number}")
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
                "worker_run_dir": str(worker_run),
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

    def cleanup_worker_runs(worker_runs, **_kwargs):
        for item in worker_runs:
            for run_dir in item.get("run_dirs") or [item.get("run_dir")]:
                if not run_dir:
                    continue
                run_path = Path(run_dir)
                cleaned_worker_runs.append(run_path.name)
                for artifact in ("poses", "ligands", "prepared"):
                    shutil.rmtree(run_path / "artifacts" / artifact, ignore_errors=True)

    monkeypatch.setattr(controller, "_cleanup_docking_artifacts", cleanup_worker_runs)
    def write_scores(path, rows):
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        pq.write_table(pa.Table.from_pylist([
            {**row, "smiles": row.get("smiles") or row.get("source_smiles") or "CCO"}
            for row in rows
        ]), path)

    monkeypatch.setattr(controller, "_save_score_table", write_scores)
    monkeypatch.setattr(controller, "_save_failed_compound_table", lambda path, _rows: Path(path).write_text("fixture"))
    from mn_molecular_design_studio import bulk_screening_results
    from mn_molecular_design_studio.molecular_database import MolecularDatabase
    database = MolecularDatabase(tmp_path / "database")
    real_publisher = bulk_screening_results.publish_final_selection

    def publish_results_only(run, user, limit, **kwargs):
        final_import_limits.append(limit)
        return real_publisher(run, user, limit, database=database, **kwargs)

    monkeypatch.setattr(bulk_screening_results, "publish_final_selection", publish_results_only)

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
    assert cleaned_worker_runs == [f"worker-{round_number}" for round_number in range(1, expected_docking_rounds)]
    final_worker = campaign / f"worker-{expected_docking_rounds}" / "artifacts"
    assert (final_worker / "poses" / "fixture.dat").is_file()
    assert (final_worker / "ligands" / "fixture.dat").is_file()
    assert (final_worker / "prepared" / "fixture.dat").is_file()
    rounds_dir = campaign / "artifacts" / "active_learning" / "rounds"
    final_round = rounds_dir / f"round-{expected_docking_rounds:04d}"
    assert (final_round / "worker-00" / "batch-0001" / "staged-input.json").is_file()
    assert (rounds_dir / f"round-{expected_docking_rounds:04d}-selection.parquet").is_file()
    for round_number in range(1, expected_docking_rounds):
        assert not (rounds_dir / f"round-{round_number:04d}" / "worker-00").exists()
        assert not (rounds_dir / f"round-{round_number:04d}-selection.parquet").exists()
    state = json.loads((campaign / "artifacts/active_learning/state.json").read_text())
    history = json.loads((campaign / "artifacts/active_learning/round_history.json").read_text())
    graph_state = json.loads((campaign / "graph_execution.json").read_text())
    decision_id = nodes["active_decision"]["id"]
    assert state["stop_reason"] == history["stop_reason"] == expected_stop
    assert history["rounds"][-1]["stop_reason"] == expected_stop
    for row in history["rounds"]:
        summaries = {group["top_n"]: group for group in row["score_summaries"]}
        assert {100, 1000, 10000}.issubset(summaries)
        assert summaries[2]["stable_round_count"] == row["stable_round_count"]
        assert summaries[2]["maximum_rank_score_change"] == row["maximum_rank_score_change"]
    assert graph_state["nodes"][decision_id]["status"] == "completed"
    assert graph_state["nodes"][decision_id]["stop_reason"] == expected_stop
    assert completed
    assert completed[-1][0][2]["metrics"]["stop_reason"] == expected_stop
    assert final_import_limits == [0]
    saved_graph = json.loads((campaign / "graph_snapshot.json").read_text())
    saved_limit = next(node["data"]["params"]["final_database_import_limit"]
                       for node in saved_graph["nodes"] if node["data"]["kind"] == "active_results")
    assert saved_limit == saved_import_limit
    external_manifest = json.loads((campaign / "artifacts/active_learning/external_results.json").read_text())
    assert external_manifest["database_selection"]["saved_import_limit_preserved_but_not_applied"] == saved_import_limit
    assert external_manifest["compound_count"] == pq.ParquetFile(
        campaign / "artifacts/active_learning/docking_scores.parquet",
    ).metadata.num_rows
    with database.connection() as connection:
        for table in ("molecules", "molecule_members", "molecule_sources", "molecule_observations"):
            assert connection.execute(f"SELECT count(*) FROM {table}").fetchone()[0] == 0

    if scenario == "stable":
        assert [row["top_n_stable"] for row in history["rounds"]] == [False, False, True, True]
        assert [row["stable_round_count"] for row in history["rounds"]] == [0, 0, 1, 2]
        assert history["rounds"][1]["same_top_order"] is False
