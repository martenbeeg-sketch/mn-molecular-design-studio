"""Curated discovery dates for tools shown in the Studio catalog.

Add an ``introduced_at`` date when a tool first becomes available in the catalog.
Update ``modified_at`` when its user-facing capability changes materially. Dates
use ISO format and are shared with the application code across workstations.
"""

from __future__ import annotations

from datetime import date, datetime, timedelta, timezone


TOOL_CATALOG_HISTORY = {
    "design.rfd4_proteina": {"introduced_at": "2026-10-21"},
    # Protein and ligand QC tiles added to the Studio catalog.
    "analysis.sequence_qc": {"introduced_at": "2026-10-21"},
    "analysis.proteinsol": {"introduced_at": "2026-10-21"},
    "analysis.biotite_interface": {"introduced_at": "2026-10-21"},
    "analysis.sequence_clustering": {"introduced_at": "2026-10-21"},
    "analysis.interface_clustering": {"introduced_at": "2026-10-21"},
    "analysis.foldseek": {"introduced_at": "2026-10-21"},
    "analysis.pep_patch": {"introduced_at": "2026-10-21"},
    "analysis.esm1v_score": {"introduced_at": "2026-10-21"},
    "analysis.esmif_score": {"introduced_at": "2026-10-21"},
    "analysis.plip": {"introduced_at": "2026-10-21"},
    "analysis.pandamap": {"introduced_at": "2026-10-21"},
    "analysis.posebusters": {"introduced_at": "2026-10-21"},
}


def tool_discovery_category(
    row: dict,
    highlight_days: int = 30,
    *,
    today: date | None = None,
) -> str | None:
    """Return ``New`` or ``Modified`` while a catalog event is recent."""
    try:
        days = int(highlight_days)
    except (TypeError, ValueError):
        days = 30
    if isinstance(highlight_days, bool) or not 1 <= days <= 365:
        days = 30

    current_date = today or datetime.now(timezone.utc).date()
    cutoff = current_date - timedelta(days=days - 1)
    workflow_ids = row.get("workflow_ids") or ()
    keys = [str(workflow_id) for workflow_id in workflow_ids]
    if row.get("key"):
        keys.insert(0, str(row["key"]))

    dates: dict[str, list[date]] = {"introduced_at": [], "modified_at": []}
    for key in keys:
        metadata = TOOL_CATALOG_HISTORY.get(key, {})
        if not isinstance(metadata, dict):
            continue
        for event in dates:
            raw_date = metadata.get(event)
            if not raw_date:
                continue
            try:
                dates[event].append(date.fromisoformat(str(raw_date)))
            except ValueError:
                continue

    if any(cutoff <= event_date <= current_date for event_date in dates["introduced_at"]):
        return "New"
    if any(cutoff <= event_date <= current_date for event_date in dates["modified_at"]):
        return "Modified"
    return None
