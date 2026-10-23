"""Opt-in collection boundaries for bulk campaigns and their child jobs.

No policy preserves existing automatic collection. A restricted run can only be
promoted through a separately created, explicitly selected collection snapshot.
Context scopes are process-local: a scheduler controller must persist/inherit the
policy before dispatching a child; context variables do not cross processes.
"""
from __future__ import annotations

import json
from contextlib import contextmanager
from contextvars import ContextVar
from pathlib import Path

_MODES = {"automatic": 0, "manual": 1, "none": 2}
_CURRENT = ContextVar("studio_database_collection_policy", default=None)
POLICY_KEY = "database_collection_policy"
PARENT_KEY = "database_collection_parent"


def validate_policy(mode):
    if not isinstance(mode, str) or mode not in _MODES:
        raise ValueError("Database collection policy must be automatic, manual, or none.")
    return mode


@contextmanager
def collection_scope(mode="manual"):
    mode = validate_policy(mode)
    outer = _CURRENT.get()
    mode = max((mode, outer or "automatic"), key=_MODES.__getitem__)
    token = _CURRENT.set(mode)
    try:
        yield
    finally:
        _CURRENT.reset(token)


def tag_created_metadata(metadata):
    mode = _CURRENT.get()
    if mode is not None:
        existing = validate_policy(metadata.get(POLICY_KEY, "automatic"))
        metadata[POLICY_KEY] = max((mode, existing), key=_MODES.__getitem__)


def _object(path):
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        return {}
    if not isinstance(data, dict):
        raise ValueError(f"Invalid collection policy metadata: {path}")
    return data


def policy_for_run(run_dir, _seen=None):
    run_dir = Path(run_dir).resolve()
    seen = set() if _seen is None else _seen
    if run_dir in seen or len(seen) >= 32:
        raise ValueError("Cycle or excessive depth in collection policy inheritance.")
    seen.add(run_dir)
    metadata = _object(run_dir / "metadata.json")
    marker = _object(run_dir / "molecular_design_studio.json")
    modes = [validate_policy(data[POLICY_KEY]) for data in (metadata, marker) if POLICY_KEY in data]
    if parent := metadata.get(PARENT_KEY):
        parent = Path(parent)
        if not parent.is_absolute() or not (parent / "metadata.json").is_file():
            raise ValueError("The collection-policy parent is unavailable.")
        modes.append(policy_for_run(parent, seen))
    return max(modes or ["automatic"], key=_MODES.__getitem__)


def automatic_collection_allowed(run_dir):
    return policy_for_run(run_dir) == "automatic"


def set_run_policy(run_dir, mode="manual", *, parent_run_dir=None):
    """Persist before execution/dispatch. Restriction cannot be loosened here."""
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json, write_json
    run_dir = Path(run_dir).resolve()
    metadata_path = run_dir / "metadata.json"
    if not metadata_path.is_file():
        raise ValueError("Create the job before assigning a collection policy.")
    metadata = read_json(metadata_path)
    modes = [validate_policy(mode), policy_for_run(run_dir)]
    if parent_run_dir is not None:
        parent = Path(parent_run_dir).resolve()
        if parent == run_dir or not (parent / "metadata.json").is_file():
            raise ValueError("Choose an existing, different parent job.")
        # Also checks for cycles through the prospective parent.
        modes.append(policy_for_run(parent, {run_dir}))
        metadata[PARENT_KEY] = str(parent)
    metadata[POLICY_KEY] = max(modes, key=_MODES.__getitem__)
    write_json(metadata_path, metadata)
    return metadata[POLICY_KEY]
