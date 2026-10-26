"""Compact cumulative ranking summaries, kept outside the molecular Database."""
from __future__ import annotations

import math
from pathlib import Path

DEFAULT_TOP_COUNTS = (100, 1000, 10000)


def summarize_rankings(current, previous, *, sizes=DEFAULT_TOP_COUNTS, tolerance=0.0, previous_summaries=()):
    """Compare ordered docking score vectors at each requested cutoff."""
    prior_counts = {item["top_n"]: item.get("stable_round_count", 0) for item in previous_summaries}
    summaries = []
    for size in sorted(set(sizes)):
        selected = current[:size]
        old = previous[:size] if previous is not None else []
        complete = len(selected) == size and len(old) == size
        same_order = all((a.get("library_id"), a.get("compound_key")) ==
                         (b.get("library_id"), b.get("compound_key")) for a, b in zip(selected, old)) if complete else None
        change = max(abs(float(a["score"]) - float(b["score"])) for a, b in zip(selected, old)) if complete else None
        stable = same_order is True and change is not None and change <= tolerance
        summaries.append({"top_n": size, "count": len(selected),
            "mean_score": math.fsum(float(row["score"]) for row in selected) / len(selected) if selected else None,
            "boundary_score": float(selected[-1]["score"]) if selected else None,
            "maximum_rank_score_change": change, "same_top_order": same_order,
            "stable_round_count": int(prior_counts.get(size, 0)) + 1 if stable else 0})
    return summaries


def rebuild_completed_convergence(run_dir):
    """One-time reconstruction for campaigns that did not save these summaries.

    This is exact only when scores were never improved after their first round.
    Check every row and the saved top vectors; refuse to invent missing history.
    The UI subsequently reads a small sidecar instead of scanning this table.
    """
    import pyarrow.parquet as pq
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json, write_json

    run_dir = Path(run_dir)
    if str(read_json(run_dir / "metadata.json").get("status")).casefold() not in {"completed", "succeeded"}:
        raise ValueError("Only completed campaigns can have their convergence history reconstructed.")
    data = run_dir / "artifacts/active_learning"
    history = read_json(data / "round_history.json")
    rounds = sorted(history.get("rounds") or [], key=lambda row: int(row["round"]))
    if not rounds:
        raise ValueError("No iteration history is available.")
    graph = read_json(run_dir / "graph_snapshot.json")
    decision = next((node.get("data", {}).get("params", {}) for node in graph.get("nodes", [])
                     if node.get("data", {}).get("kind") == "active_decision"), {})
    sizes = sorted({*DEFAULT_TOP_COUNTS, int(decision.get("top_n") or 100)})
    vectors = {int(row["round"]): [] for row in rounds}
    counts = dict.fromkeys(vectors, 0)
    ascending = history.get("score_direction", "Lower is better") == "Lower is better"
    last_score = None
    table = pq.ParquetFile(data / "docking_scores.parquet")
    for batch in table.iter_batches(columns=["library_id", "compound_key", "score", "first_round", "best_round"]):
        for row in batch.to_pylist():
            first = int(row["first_round"])
            if first != int(row["best_round"]):
                raise ValueError("Scores changed after their first round; earlier rankings cannot be reconstructed from the final table.")
            score = float(row["score"])
            if not math.isfinite(score) or first not in vectors:
                raise ValueError("The final score table contains invalid scores or iteration numbers.")
            if last_score is not None and ((ascending and score < last_score) or (not ascending and score > last_score)):
                raise ValueError("The final score table is not in the declared ranking order.")
            last_score = score
            for number, vector in vectors.items():
                if first <= number:
                    counts[number] += 1
                    if len(vector) < max(sizes):
                        vector.append({key: row[key] for key in ("library_id", "compound_key", "score")})
    previous, previous_summaries, output = None, [], []
    for saved in rounds:
        number = int(saved["round"])
        vector = vectors[number]
        if saved.get("unique_tested_count") is not None and counts[number] != int(saved["unique_tested_count"]):
            raise ValueError(f"Reconstructed compound count differs from saved iteration {number}.")
        for old, new in zip(saved.get("top_scores") or [], vector):
            if any(old.get(key) != new.get(key) for key in ("library_id", "compound_key", "score")):
                raise ValueError(f"Reconstructed ranking differs from saved iteration {number}.")
        summaries = summarize_rankings(vector, previous, sizes=sizes,
            tolerance=float(decision.get("score_tolerance") or 0), previous_summaries=previous_summaries)
        output.append({"round": number, "score_summaries": summaries})
        previous, previous_summaries = vector, summaries
    payload = {"schema_version": "studio.screening-convergence.v1", "source": str(data / "docking_scores.parquet"),
               "reconstruction": "Validated cumulative rankings; no score updated after first docking", "rounds": output}
    write_json(data / "score_convergence.json", payload)
    return payload
