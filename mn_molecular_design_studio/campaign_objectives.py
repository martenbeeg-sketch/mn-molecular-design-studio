"""Stopping semantics shared by tool campaigns and graph campaign scopes.

The existing controllers own dispatch, receipts and recovery. These functions
only normalize objectives and evaluate count/time limits against durable state.
"""
from __future__ import annotations

import math
from datetime import datetime, timezone


OBJECTIVE_LABELS = {
    "generated": "Generate a fixed count",
    "accepted": "Reach a passing-design count",
    "time": "Run until time is up",
}


def timestamp(value: str) -> datetime:
    try:
        parsed = datetime.fromisoformat(str(value).replace("Z", "+00:00"))
    except (TypeError, ValueError) as exc:
        raise ValueError("deadline_at must be an ISO date/time with a timezone offset.") from exc
    if parsed.tzinfo is None or parsed.utcoffset() is None:
        raise ValueError("deadline_at must include a timezone offset, for example +01:00.")
    return parsed.astimezone(timezone.utc)


def normalize_objective(config: dict) -> dict:
    result = dict(config)
    kind = str(result.get("objective_type", "accepted")).strip().lower()
    if kind not in OBJECTIVE_LABELS:
        raise ValueError("objective_type must be generated, accepted or time.")
    result["objective_type"] = kind
    duration = result.get("max_runtime_seconds")
    if duration is not None and (
        type(duration) not in {int, float} or not math.isfinite(duration) or duration <= 0
    ):
        raise ValueError("max_runtime_seconds must be finite and positive when provided.")
    if result.get("deadline_at") is not None:
        result["deadline_at"] = timestamp(result["deadline_at"]).isoformat()
    if kind == "time":
        if duration is None and not result.get("deadline_at"):
            raise ValueError("A time objective requires max_runtime_seconds or deadline_at.")
        if result.get("accepted_target") is not None or result.get("generated_target") is not None:
            raise ValueError("A time objective has no count target; remove accepted_target/generated_target.")
    return result


def count_target(config: dict) -> int | None:
    if config.get("objective_type", "accepted") == "time":
        return None
    field = "generated_target" if config.get("objective_type") == "generated" else "accepted_target"
    return int(config.get(field, config.get("accepted_target", 1)))


def count_reached(config: dict, count: int) -> bool:
    target = count_target(config)
    return target is not None and count >= target


def time_progress(config: dict, created_at: str, *, now: datetime | None = None) -> dict:
    now = now or datetime.now(timezone.utc)
    if now.tzinfo is None:
        raise ValueError("Campaign time evaluation requires a timezone-aware clock.")
    elapsed = max(0.0, (now - timestamp(created_at)).total_seconds())
    remaining = []
    if config.get("max_runtime_seconds") is not None:
        remaining.append(float(config["max_runtime_seconds"]) - elapsed)
    if config.get("deadline_at"):
        remaining.append((timestamp(config["deadline_at"]) - now).total_seconds())
    seconds_left = max(0.0, min(remaining)) if remaining else None
    return {
        "elapsed_seconds": elapsed,
        "remaining_seconds": seconds_left,
        "time_expired": bool(remaining) and min(remaining) <= 0,
        "deadline_at": config.get("deadline_at"),
        "max_runtime_seconds": config.get("max_runtime_seconds"),
    }


def objective_progress(config: dict, count: int, created_at: str, *, now: datetime | None = None) -> dict:
    clock = time_progress(config, created_at, now=now)
    count_met = count_reached(config, count)
    return {
        **clock, "objective_type": config.get("objective_type", "accepted"),
        "objective_target": count_target(config), "objective_count": count,
        "objective_met": clock["time_expired"] if config.get("objective_type") == "time" else count_met,
        "count_target_met": count_met,
    }


def configuration_schema() -> dict:
    """Public stopping fields; engine-specific resources/stages remain available."""
    return {
        "type": ["object", "null"],
        "description": "Optional campaign scope. Count/time objectives need no acceptance filter.",
        "properties": {
            "enabled": {"type": "boolean"},
            "objective_type": {"enum": list(OBJECTIVE_LABELS)},
            "accepted_target": {"type": ["integer", "null"], "minimum": 1,
                                "description": "Generated or passing count for legacy native campaign APIs; null for time."},
            "max_attempts": {"type": "integer", "minimum": 1,
                             "description": "Maximum requested design outputs across the campaign, independent of job and engine microbatch sizes."},
            "native_designs_per_job": {"type": "integer", "minimum": 1,
                                        "description": "Native designs requested in one scheduler job; currently configurable for RFD4."},
            "max_runtime_seconds": {"type": "number", "exclusiveMinimum": 0},
            "deadline_at": {"type": "string", "format": "date-time"},
            "validation_stages": {"type": "array", "items": {"type": "object"}},
        },
    }
