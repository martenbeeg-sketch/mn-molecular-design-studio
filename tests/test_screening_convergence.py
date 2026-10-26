import json

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from mn_molecular_design_studio.screening_convergence import summarize_rankings, rebuild_completed_convergence


def vector(scores, offset=0):
    return [{"library_id": "library", "compound_key": index + offset, "score": score}
            for index, score in enumerate(scores)]


def test_summary_compares_energy_and_order_for_each_cutoff():
    before = vector([-10, -8, -6, -4])
    now = vector([-10, -8, -7, -4])
    rows = summarize_rankings(now, before, sizes=(1, 2, 4), tolerance=0.05,
        previous_summaries=[{"top_n": 2, "stable_round_count": 2}])
    assert [row["mean_score"] for row in rows] == [-10, -9, -7.25]
    assert [row["maximum_rank_score_change"] for row in rows] == [0, 0, 1]
    assert [row["stable_round_count"] for row in rows] == [1, 3, 0]
    changed_order = summarize_rankings(vector([-10, -8], offset=1), before, sizes=(2,))[0]
    assert changed_order["maximum_rank_score_change"] == 0
    assert changed_order["same_top_order"] is False and changed_order["stable_round_count"] == 0


def test_partial_top_set_does_not_claim_convergence():
    row = summarize_rankings(vector([-10]), vector([-10]), sizes=(100,))[0]
    assert row["mean_score"] == -10 and row["count"] == 1
    assert row["maximum_rank_score_change"] is None
    assert row["stable_round_count"] == 0


def campaign(tmp_path, *, updated=False):
    data = tmp_path / "artifacts/active_learning"
    data.mkdir(parents=True)
    (tmp_path / "metadata.json").write_text(json.dumps({"status": "completed"}))
    (data / "round_history.json").write_text(json.dumps({"score_direction": "Lower is better", "rounds": [
        {"round": 1, "unique_tested_count": 2, "top_scores": vector([-9, -7])},
        {"round": 2, "unique_tested_count": 3, "top_scores": [vector([-11], offset=2)[0], *vector([-9, -7])]},
    ]}))
    (tmp_path / "graph_snapshot.json").write_text(json.dumps({"nodes": [{"data": {"kind": "active_decision",
        "params": {"top_n": 2, "score_tolerance": 0.05}}}]}))
    pq.write_table(pa.Table.from_pylist([
        {**vector([-11], offset=2)[0], "first_round": 1 if updated else 2, "best_round": 2},
        {**vector([-9])[0], "first_round": 1, "best_round": 1},
        {**vector([-7], offset=1)[0], "first_round": 1, "best_round": 1},
    ]), data / "docking_scores.parquet")
    return data


def test_reconstruction_checks_original_counts_and_vectors(tmp_path):
    data = campaign(tmp_path)
    result = rebuild_completed_convergence(tmp_path)
    first, second = [row["score_summaries"][0] for row in result["rounds"]]
    assert first["top_n"] == 2 and first["mean_score"] == -8
    assert second["mean_score"] == -10 and second["maximum_rank_score_change"] == 2
    assert (data / "score_convergence.json").is_file()
    assert len(result["rounds"][0]["score_summaries"]) == 4


def test_reconstruction_refuses_to_rewrite_history_after_rescoring(tmp_path):
    data = campaign(tmp_path, updated=True)
    with pytest.raises(ValueError, match="Scores changed"):
        rebuild_completed_convergence(tmp_path)
    assert not (data / "score_convergence.json").exists()
