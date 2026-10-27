"""Scientific validity is independent of process completion and acceptance."""
from __future__ import annotations

import math


def finite_number(value, *, allow_bool=False):
    if value is None or (isinstance(value, bool) and not allow_bool):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError):
        return None
    return number if math.isfinite(number) else None


def measurement_outcome(metrics, required, error="", *, boolean_fields=()):
    fields = {key: "missing" if metrics.get(key) is None else
              "valid" if finite_number(metrics[key], allow_bool=key in boolean_fields) is not None else "invalid"
              for key in required}
    invalid = [f"{key}: {state}" for key, state in fields.items() if state != "valid"]
    return {"status": "error" if error else "invalid" if "invalid" in fields.values() else
            "missing" if invalid else "measured", "fields": fields,
            "error": str(error or ("Required measurements unavailable: " + "; ".join(invalid) if invalid else ""))}


def measurement_summary(outcomes):
    valid = sum(item["status"] == "measured" for item in outcomes)
    return {"candidate_count": len(outcomes), "scored_count": valid,
            "error_count": sum(item["status"] == "error" for item in outcomes),
            "missing_count": sum(item["status"] == "missing" for item in outcomes),
            "invalid_count": sum(item["status"] == "invalid" for item in outcomes),
            "scientific_status": "failed" if not valid else "valid" if valid == len(outcomes) else "partial"}
