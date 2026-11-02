"""Resolve a candidate's source-native accept/reject decision without guessing.

Different engines publish their native filter result under different names.
This adapter preserves generic Studio fields and recognizes RFD4's native PPI
gate only when RFD4 reports that all required PPI metrics were available.
"""
from __future__ import annotations

from collections.abc import Iterator


def _boolean(value):
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)) and value in (0, 1):
        return bool(value)
    if isinstance(value, str):
        normalized = value.strip().casefold()
        if normalized in {"true", "1", "pass", "passed", "accepted", "yes"}:
            return True
        if normalized in {"false", "0", "fail", "failed", "rejected", "no"}:
            return False
    return None


def _candidate_sources(candidate: dict) -> Iterator[dict]:
    """Yield the normalized candidate, then its preserved source candidate."""
    yield candidate
    raw = candidate.get("raw_metadata")
    source = raw.get("source_candidate") if isinstance(raw, dict) else None
    if isinstance(source, dict):
        yield source


def native_filter_evidence(candidate: dict) -> dict:
    """Return a tri-state native filter result and its evidence source.

    ``value`` is ``True`` or ``False`` only when the source explicitly reports
    a complete decision. Missing or incomplete evidence remains unavailable.
    """
    if not isinstance(candidate, dict):
        return {"value": None, "status": "unavailable", "source": None,
                "reason": "candidate evidence is missing"}

    sources = list(_candidate_sources(candidate))
    for source in sources:
        metrics = source.get("metrics")
        if not isinstance(metrics, dict):
            continue
        for field in ("native_pass_filters", "pass_filters", "passes_filters", "passed_filters"):
            if field not in metrics or metrics[field] is None:
                continue
            value = _boolean(metrics[field])
            if value is not None:
                return {"value": value, "status": "passed" if value else "failed",
                        "source": field, "reason": ""}
            return {"value": None, "status": "unavailable", "source": field,
                    "reason": "native filter value is not boolean"}

    for source in sources:
        metrics = source.get("metrics")
        if not isinstance(metrics, dict) or "ppi.pass" not in metrics:
            continue
        complete = _boolean(metrics.get("ppi.all_metrics_present"))
        if complete is not True:
            return {"value": None, "status": "unavailable", "source": "ppi.pass",
                    "reason": "native PPI metrics incomplete"}
        value = _boolean(metrics.get("ppi.pass"))
        if value is None:
            return {"value": None, "status": "unavailable", "source": "ppi.pass",
                    "reason": "native PPI pass value is invalid"}
        return {"value": value, "status": "passed" if value else "failed",
                "source": "ppi.pass", "reason": ""}

    for source in sources:
        raw = source.get("raw_metadata")
        if not isinstance(raw, dict):
            continue
        kind = str(raw.get("bindcraft_output_kind") or "").strip().casefold()
        if kind in {"accepted", "rejected", "filtered"}:
            value = kind == "accepted"
            return {"value": value, "status": "passed" if value else "failed",
                    "source": "bindcraft_output_kind", "reason": ""}
        campaign = raw.get("distributed_campaign")
        if isinstance(campaign, dict):
            value = campaign.get("native_acceptance")
            if value in {"accepted", "rejected"}:
                passed = value == "accepted"
                return {"value": passed, "status": "passed" if passed else "failed",
                        "source": "distributed_campaign.native_acceptance", "reason": ""}
            if value not in (None, "unavailable"):
                return {"value": None, "status": "unresolved",
                        "source": "distributed_campaign.native_acceptance",
                        "reason": "native campaign acceptance value is unsupported"}

    return {"value": None, "status": "unavailable", "source": None,
            "reason": "native filter evidence is not present"}
