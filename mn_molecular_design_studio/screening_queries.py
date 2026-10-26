"""Bounded, streaming reads of score-ordered external screening tables."""
from __future__ import annotations

import math
from pathlib import Path


FILES = {"scores": "docking_scores.parquet", "failures": "failed_compounds.parquet",
         "rounds": "round_history.json", "convergence": "score_convergence.json",
         "manifest": "external_results.json", "model": "model.pt"}


def result_file(run_dir, user, kind):
    from mn_molecular_design_studio.molecular_database_results import _authorized_run
    from mn_molecular_design_studio.bulk_screening_results import RESULT_DIRECTORY
    selected = _authorized_run(run_dir, user)[0]
    if kind not in FILES:
        raise FileNotFoundError("Screening result file not found.")
    root = selected / RESULT_DIRECTORY
    path = root / FILES[kind]
    if not path.is_file() or not path.resolve().is_relative_to(root.resolve()):
        raise FileNotFoundError("Screening result file not found.")
    return path


def query_results(run_dir, user, *, offset=0, limit=100, query="", score_min=None,
                  score_max=None, compound_ids=None):
    """Preserve cumulative docking rank; never load the entire result table.

    Unfiltered pages skip whole Parquet row groups using metadata. Filtered
    pages scan only until the requested page plus one matching row is found.
    Their total is deliberately unknown rather than forcing a complete scan.
    """
    import pyarrow.parquet as pq
    if type(offset) is not int or offset < 0 or type(limit) is not int or not 1 <= limit <= 1000:
        raise ValueError("Use a non-negative offset and page size between 1 and 1,000.")
    for value in (score_min, score_max):
        if value is not None and not math.isfinite(value):
            raise ValueError("Score bounds must be finite.")
    if score_min is not None and score_max is not None and score_min > score_max:
        raise ValueError("score_min cannot exceed score_max.")
    ids = set(compound_ids or [])
    text = str(query).strip().casefold()
    filtered = bool(text or ids or score_min is not None or score_max is not None)
    path = result_file(run_dir, user, "scores")
    parquet = pq.ParquetFile(path)
    total = parquet.metadata.num_rows
    rows, rank_base, matched = [], 0, 0
    for group in range(parquet.metadata.num_row_groups):
        count = parquet.metadata.row_group(group).num_rows
        if not filtered and rank_base + count <= offset:
            rank_base += count
            continue
        for batch in parquet.iter_batches(batch_size=2048, row_groups=[group]):
            for index, row in enumerate(batch.to_pylist()):
                rank = rank_base + index + 1
                if not filtered:
                    if rank <= offset:
                        continue
                else:
                    if ids and str(row.get("compound_id")) not in ids:
                        continue
                    if text and not any(text in str(row.get(key) or "").casefold()
                                        for key in ("compound_id", "compound_name", "name", "smiles")):
                        continue
                    score = row.get("score", row.get("best_score"))
                    if score_min is not None or score_max is not None:
                        if score is None or not math.isfinite(float(score)):
                            continue
                        if score_min is not None and float(score) < score_min:
                            continue
                        if score_max is not None and float(score) > score_max:
                            continue
                    matched += 1
                    if matched <= offset:
                        continue
                rows.append({**row, "rank": rank})
                if len(rows) > limit:
                    return _page(rows, offset, limit, total, filtered, parquet.schema_arrow.names)
            rank_base += batch.num_rows
    return _page(rows, offset, limit, total, filtered, parquet.schema_arrow.names)


def _page(rows, offset, limit, total, filtered, columns):
    has_more = len(rows) > limit
    return {"items": rows[:limit], "offset": offset, "limit": limit, "has_more": has_more,
            "next_offset": offset + limit if has_more else None,
            "total": None if filtered else total, "source_total": total,
            "ranking": "cumulative docking score order", "columns": ["rank", *columns]}
