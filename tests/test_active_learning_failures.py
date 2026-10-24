import csv

from mn_molecular_design_studio.engine.workflow_engine.workflows.active_learning_controller import (
    _count_library_compounds,
    _persist_round_failures,
    _persist_round_scores,
    _read_selection_table,
    _read_worker_scores,
    _score_database,
    _summarize_round_workers,
    _write_selection_table,
)


def _compound(key, engine_id):
    return {
        "library_id": "screening-library",
        "compound_key": key,
        "compound_id": f"compound-{key}",
        "compound_name": f"Compound {key}",
        "source_smiles": "CCO",
        "fingerprint": bytes(128),
        "selection_reason": "random_exploration",
        "predicted_probability": None,
        "engine_id": engine_id,
    }


def test_worker_scores_separate_failures_and_missing_rows(tmp_path):
    score_path = tmp_path / "scores.csv"
    with score_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["compound_id", "status", "best_score_kcal_mol", "error"])
        writer.writeheader()
        writer.writerow({"compound_id": "L0000001", "status": "completed", "best_score_kcal_mol": "-8.1"})
        writer.writerow({"compound_id": "L0000002", "status": "failed", "error": "Meeko failed"})

    rows = {"L0000001": _compound(1, "L0000001"), "L0000002": _compound(2, "L0000002"),
            "L0000003": _compound(3, "L0000003")}
    scored, failures = _read_worker_scores({
        "score_path": str(score_path), "node_id": "worker-a", "node_name": "Worker A",
        "worker_run_dir": "/runs/job-a", "job_code": "job-a", "cpu_cores": 8,
    }, rows)

    assert [(row["compound_key"], row["score"]) for row in scored] == [(1, -8.1)]
    assert [row["compound_key"] for row in failures] == [2, 3]
    assert failures[0]["error"] == "Meeko failed"
    assert failures[1]["error"] == "Docking worker did not return a score-table row."


def test_failure_rows_are_counted_but_never_used_as_docking_labels(tmp_path):
    db = _score_database(tmp_path / "scores.sqlite3")
    success = {**_compound(1, "L0000001"), "score": -8.1,
               "worker_node_id": "worker-a", "worker_node_name": "Worker A"}
    failed = {**_compound(2, "L0000002"), "error": "Ligand preparation failed",
              "worker_node_id": "worker-b", "worker_node_name": "Worker B"}

    _persist_round_scores(db, [success], round_number=1, score_direction="Lower is better")
    _persist_round_failures(db, [failed], round_number=1)

    assert _count_library_compounds(db, "screening-library") == 2
    assert db.execute("SELECT COUNT(*) FROM compounds").fetchone()[0] == 1
    assert db.execute("SELECT COUNT(*) FROM failed_compounds").fetchone()[0] == 1
    summary = _summarize_round_workers([success], [failed])
    assert {row["node_id"]: (row["scored_count"], row["failed_count"]) for row in summary} == {
        "worker-a": (1, 0), "worker-b": (0, 1),
    }

    recovered = {**success, "compound_key": 2, "compound_id": "compound-2", "score": -7.4}
    _persist_round_scores(db, [recovered], round_number=2, score_direction="Lower is better")
    assert db.execute("SELECT COUNT(*) FROM failed_compounds").fetchone()[0] == 0
    assert db.execute("SELECT COUNT(*) FROM compounds").fetchone()[0] == 2
    db.close()


def test_selection_table_preserves_predicted_and_random_reasons(tmp_path):
    rows = [
        {**_compound(1, "predicted-1"), "selection_reason": "predicted_hit",
         "predicted_probability": 0.99},
        {**_compound(2, "random-1"), "selection_reason": "random_exploration",
         "predicted_probability": None},
    ]
    selection_path = tmp_path / "selected.parquet"

    _write_selection_table(selection_path, rows, "screening-library", reason="model_selected")

    saved_rows = _read_selection_table(selection_path)
    assert [row["selection_reason"] for row in saved_rows] == ["predicted_hit", "random_exploration"]
    assert [row["predicted_probability"] for row in saved_rows] == [0.99, None]
