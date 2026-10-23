"""Live progress readers for scientific workflows that persist campaign state."""
from __future__ import annotations

import json
from pathlib import Path


_BINDCRAFT2_STATE = Path("artifacts/raw/bindcraft2/output/.campaign_state.json")
_BINDCRAFT2_SETTINGS = (
    Path("artifacts/raw/bindcraft2/output/workers/campaign_settings.json"),
    Path("artifacts/raw/bindcraft2/settings.json"),
    Path("input.json"),
)


def bindcraft2_progress_mtime_ns(run_dir: Path) -> int:
    """Return the state-file timestamp used to trigger a fresh workstation heartbeat."""
    try:
        return (run_dir / _BINDCRAFT2_STATE).stat().st_mtime_ns
    except OSError:
        return 0


def _read_object(path: Path) -> dict:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError, UnicodeDecodeError):
        return {}
    return value if isinstance(value, dict) else {}


def _setting(settings: dict, key: str) -> int | None:
    value = settings.get(key)
    if value is None:
        for container in ("params", "parameters"):
            nested = settings.get(container)
            if isinstance(nested, dict) and nested.get(key) is not None:
                value = nested[key]
                break
    try:
        return max(0, int(value)) if value is not None else None
    except (TypeError, ValueError, OverflowError):
        return None


def bindcraft2_progress(run_dir: Path) -> dict | None:
    """Read BindCraft2's atomic state file and the requested campaign limits.

    BindCraft increments ``trajectories`` when a worker claims an attempt. The
    per-stage ``terminated`` counters tell us which of those attempts have
    finished, so active trajectories can be shown separately.
    """
    run_dir = Path(run_dir)
    state_path = run_dir / _BINDCRAFT2_STATE
    state = _read_object(state_path)
    if not state:
        return None

    settings = {}
    for relative_path in _BINDCRAFT2_SETTINGS:
        settings = _read_object(run_dir / relative_path)
        if settings:
            break

    attempted = _setting(state, "trajectories")
    if attempted is None:
        attempts = state.get("attempted")
        attempted = len(attempts) if isinstance(attempts, list) else 0

    rejections = state.get("rejections")
    terminated = rejections.get("terminated") if isinstance(rejections, dict) else {}
    completed = sum(
        max(0, int(count))
        for count in terminated.values()
        if isinstance(count, (int, float)) and not isinstance(count, bool)
    ) if isinstance(terminated, dict) else 0

    accepted = _setting(state, "accepted") or 0
    metadata = _read_object(run_dir / "metadata.json")
    terminal = str(metadata.get("status") or "").casefold() in {
        "completed", "succeeded", "failed", "cancelled", "canceled",
    }
    unfinished = max(0, attempted - completed)
    try:
        updated_at = state_path.stat().st_mtime
    except OSError:
        updated_at = None

    return {
        "engine": "bindcraft2",
        "trajectories_started": attempted,
        "trajectories_completed": completed,
        "trajectories_in_progress": 0 if terminal else unfinished,
        "trajectories_interrupted": unfinished if terminal else 0,
        "trajectory_limit": _setting(settings, "max_trajectories"),
        "designs_accepted": accepted,
        "design_target": _setting(settings, "number_of_final_designs"),
        "updated_at": updated_at,
    }
