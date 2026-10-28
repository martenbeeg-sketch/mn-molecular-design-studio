"""Durable acceptance accounting for distributed BindCraft2 campaigns.

The ledger lives beside its parent workflow run and uses an OS file lock plus
atomic replacement. It is deliberately not a second scheduler or database.
"""
from __future__ import annotations

import fcntl
import hashlib
import json
import os
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterator
from uuid import uuid4

from mn_molecular_design_studio.engine.workflow_engine.core.jobs import utc_now


SCHEMA_VERSION = 1
ATTEMPT_STATES = {
    "authorized", "queued", "running", "completed", "failed", "interrupted",
    "cancelled", "not_started",
}
TERMINATION_REASONS = {
    "target_reached", "budget_exhausted", "runtime_exhausted", "cancelled", "failed",
}


def _atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    temporary = Path(name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


@contextmanager
def _locked(path: Path) -> Iterator[None]:
    lock_path = path.with_suffix(path.suffix + ".lock")
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("a+b") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def _read(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise RuntimeError(f"Campaign ledger is missing or unreadable: {path}") from exc
    if not isinstance(value, dict) or value.get("schema_version") != SCHEMA_VERSION:
        raise RuntimeError(f"Unsupported campaign ledger schema: {path}")
    return value


def validate_campaign_config(config: dict[str, Any]) -> dict[str, Any]:
    """Validate global counters and execution bounds without inventing science thresholds."""
    if not isinstance(config, dict):
        raise ValueError("distributed_campaign must be an object.")
    acceptance_level = str(config.get("acceptance_level", "native")).strip().lower()
    if acceptance_level not in {"native", "validated", "diversity_qualified"}:
        raise ValueError("acceptance_level must be native, validated, or diversity_qualified.")

    def positive_int(key: str, default: int) -> int:
        value = config.get(key, default)
        if type(value) is not int or value < 1:
            raise ValueError(f"{key} must be a positive integer.")
        return value

    maximum = positive_int("max_attempts", 1)
    target = positive_int("accepted_target", 1)
    batch = positive_int("batch_size", 1)
    concurrency = positive_int("max_concurrent_batches", 1)
    if target > maximum:
        raise ValueError("accepted_target cannot exceed max_attempts.")
    seed = config.get("seed", 0)
    if type(seed) is not int or seed < 0:
        raise ValueError("seed must be a nonnegative integer.")
    runtime = config.get("max_runtime_seconds")
    if runtime is not None and (type(runtime) not in {int, float} or runtime <= 0):
        raise ValueError("max_runtime_seconds must be positive when provided.")
    requirements = config.get("validation_stages", [])
    if not isinstance(requirements, list) or any(not isinstance(item, dict) for item in requirements):
        raise ValueError("validation_stages must be a list of workflow objects.")
    diversity = config.get("diversity", {})
    if not isinstance(diversity, dict):
        raise ValueError("diversity must be an object.")
    identity = diversity.get("sequence_identity_max")
    if identity is not None and (
        type(identity) not in {float, int} or not 0 < float(identity) <= 1
    ):
        raise ValueError("diversity.sequence_identity_max must be between 0 and 1.")
    cluster_cap = diversity.get("max_per_cluster")
    if cluster_cap is not None and (type(cluster_cap) is not int or cluster_cap < 1):
        raise ValueError("diversity.max_per_cluster must be a positive integer.")
    return {
        **config,
        "acceptance_level": acceptance_level,
        "max_attempts": maximum,
        "accepted_target": target,
        "batch_size": batch,
        "max_concurrent_batches": concurrency,
        "seed": seed,
        "validation_stages": requirements,
        "diversity": diversity,
    }


class CampaignLedger:
    """Atomic attempt reservations, idempotent result receipts, and global counts."""

    def __init__(self, path: Path):
        self.path = Path(path)

    @classmethod
    def create(cls, path: Path, config: dict[str, Any], *, campaign_id: str | None = None) -> "CampaignLedger":
        ledger = cls(path)
        validated = validate_campaign_config(config)
        with _locked(ledger.path):
            if ledger.path.exists():
                existing = _read(ledger.path)
                if existing.get("config") != validated:
                    raise ValueError("A campaign ledger already exists with different scientific settings.")
                return ledger
            now = utc_now()
            _atomic_json(ledger.path, {
                "schema_version": SCHEMA_VERSION,
                "campaign_id": campaign_id or str(uuid4()),
                "config": validated,
                "status": "running",
                "termination_reason": None,
                "created_at": now,
                "updated_at": now,
                "next_ordinal": 1,
                "attempts": {},
                "receipts": {},
                "candidates": {},
                "events": [{"type": "campaign_created", "at": now}],
            })
        return ledger

    def snapshot(self) -> dict[str, Any]:
        with _locked(self.path):
            return _read(self.path)

    def counters(self) -> dict[str, int]:
        state = self.snapshot()
        attempts = list(state["attempts"].values())
        candidates = list(state["candidates"].values())
        return {
            "attempts_authorized": len(attempts),
            "attempts_started": sum(bool(item.get("started_at")) for item in attempts),
            "attempts_completed": sum(item.get("status") == "completed" for item in attempts),
            "native_accepted_candidates": sum(item.get("native_acceptance") == "accepted" for item in candidates),
            "validated_accepted_candidates": sum(item.get("validated_acceptance") == "accepted" for item in candidates),
            "diversity_qualified_accepted_candidates": sum(item.get("diversity_acceptance") == "accepted" for item in candidates),
            "rejected_candidates": sum(
                item.get("native_acceptance") == "rejected"
                or item.get("validated_acceptance") == "rejected"
                or item.get("diversity_acceptance") == "rejected"
                for item in candidates
            ),
            "failed_interrupted_attempts": sum(item.get("status") in {"failed", "interrupted"} for item in attempts),
        }

    @staticmethod
    def _seed(campaign_id: str, configured_seed: int, ordinal: int) -> int:
        digest = hashlib.sha256(f"{campaign_id}:{configured_seed}:{ordinal}".encode()).digest()
        return int.from_bytes(digest[:4], "big") % 2_147_483_646 + 1

    def authorize(self, requested: int, *, batch_id: str) -> list[dict[str, Any]]:
        if type(requested) is not int or requested < 1 or not batch_id:
            raise ValueError("Attempt authorization needs a positive count and batch identity.")
        with _locked(self.path):
            state = _read(self.path)
            if state.get("status") != "running":
                return []
            remaining = state["config"]["max_attempts"] - len(state["attempts"])
            count = min(requested, max(0, remaining))
            created = []
            for _ in range(count):
                ordinal = int(state["next_ordinal"])
                state["next_ordinal"] = ordinal + 1
                attempt_id = f"{state['campaign_id']}:a{ordinal:08d}"
                attempt = {
                    "attempt_id": attempt_id,
                    "ordinal": ordinal,
                    "batch_id": batch_id,
                    "seed": self._seed(state["campaign_id"], state["config"]["seed"], ordinal),
                    "status": "authorized",
                    "authorized_at": utc_now(),
                    "started_at": None,
                    "completed_at": None,
                    "child_run_dir": None,
                    "trajectory_index": None,
                    "retry_count": 0,
                    "receipts": [],
                }
                state["attempts"][attempt_id] = attempt
                created.append(dict(attempt))
            state["updated_at"] = utc_now()
            state["events"].append({
                "type": "attempts_authorized", "batch_id": batch_id,
                "attempt_ids": [item["attempt_id"] for item in created], "at": state["updated_at"],
            })
            _atomic_json(self.path, state)
            return created

    def transition(self, attempt_id: str, status: str, **details: Any) -> dict[str, Any]:
        if status not in ATTEMPT_STATES:
            raise ValueError(f"Unknown attempt status: {status}")
        with _locked(self.path):
            state = _read(self.path)
            attempt = state["attempts"].get(attempt_id)
            if attempt is None:
                raise KeyError(f"Unknown campaign attempt: {attempt_id}")
            previous = attempt["status"]
            terminal = {"completed", "failed", "interrupted", "cancelled", "not_started"}
            if previous in terminal and status != previous:
                raise ValueError(f"Terminal attempt {attempt_id} cannot transition from {previous} to {status}.")
            if status == "running" and not attempt.get("started_at"):
                attempt["started_at"] = utc_now()
            if status in terminal:
                attempt["completed_at"] = utc_now()
            if status == "authorized" and previous != "authorized":
                raise ValueError("An existing attempt cannot be reset to authorized.")
            attempt.update(details)
            attempt["status"] = status
            state["updated_at"] = utc_now()
            state["events"].append({
                "type": "attempt_transition", "attempt_id": attempt_id,
                "from": previous, "to": status, "at": state["updated_at"],
            })
            _atomic_json(self.path, state)
            return dict(attempt)

    def retry(self, attempt_id: str, *, child_run_dir: str | None = None) -> dict[str, Any]:
        with _locked(self.path):
            state = _read(self.path)
            attempt = state["attempts"].get(attempt_id)
            if attempt is None:
                raise KeyError(f"Unknown campaign attempt: {attempt_id}")
            if attempt["status"] not in {"failed", "interrupted", "queued"}:
                raise ValueError("Only failed, interrupted, or queued attempts can be retried.")
            attempt["retry_count"] = int(attempt.get("retry_count", 0)) + 1
            attempt.setdefault("retry_history", []).append({
                "at": utc_now(), "child_run_dir": child_run_dir or attempt.get("child_run_dir"),
                "seed": attempt["seed"],
            })
            attempt["child_run_dir"] = child_run_dir or attempt.get("child_run_dir")
            attempt["status"] = "queued"
            attempt["started_at"] = None
            attempt["completed_at"] = None
            state["updated_at"] = utc_now()
            state["events"].append({"type": "attempt_retry", "attempt_id": attempt_id, "at": state["updated_at"]})
            _atomic_json(self.path, state)
            return dict(attempt)

    def record_receipt(
        self,
        attempt_id: str,
        receipt_id: str,
        *,
        candidates: list[dict[str, Any]],
        native_state: dict[str, Any],
    ) -> bool:
        """Record one verified child result. Duplicate receipts/results are no-ops."""
        if not receipt_id or not isinstance(native_state, dict):
            raise ValueError("A result receipt needs an identity and verified native state.")
        with _locked(self.path):
            state = _read(self.path)
            if receipt_id in state["receipts"]:
                return False
            attempt = state["attempts"].get(attempt_id)
            if attempt is None:
                raise KeyError(f"Unknown campaign attempt: {attempt_id}")
            if attempt.get("status") == "completed":
                state["receipts"][receipt_id] = {
                    "attempt_id": attempt_id, "duplicate_attempt_result": True, "at": utc_now(),
                }
                _atomic_json(self.path, state)
                return False
            attempt["native_state"] = dict(native_state)
            attempt["candidate_ids"] = []
            for candidate in candidates:
                candidate_id = str(candidate.get("candidate_id") or "")
                if not candidate_id:
                    raise ValueError("Verified campaign candidates need nonempty stable IDs.")
                if candidate_id in state["candidates"]:
                    continue
                native_acceptance = (
                    "accepted" if candidate.get("metrics", {}).get("native_pass_filters") is True else "rejected"
                )
                state["candidates"][candidate_id] = {
                    "candidate_id": candidate_id,
                    "attempt_id": attempt_id,
                    "native_acceptance": native_acceptance,
                    "validated_acceptance": "pending",
                    "diversity_acceptance": "pending",
                }
                attempt["candidate_ids"].append(candidate_id)
            attempt["receipts"].append(receipt_id)
            attempt["status"] = "completed"
            attempt["completed_at"] = utc_now()
            state["receipts"][receipt_id] = {"attempt_id": attempt_id, "at": utc_now()}
            state["updated_at"] = utc_now()
            state["events"].append({
                "type": "result_verified", "attempt_id": attempt_id, "receipt_id": receipt_id,
                "candidate_count": len(attempt["candidate_ids"]), "at": state["updated_at"],
            })
            _atomic_json(self.path, state)
            return True

    def set_candidate_acceptance(
        self, candidate_id: str, *, validated: str | None = None,
        diversity: str | None = None, reason: str = "",
    ) -> None:
        allowed = {None, "pending", "accepted", "rejected", "failed"}
        if validated not in allowed or diversity not in allowed:
            raise ValueError("Candidate acceptance must be pending, accepted, rejected, or failed.")
        with _locked(self.path):
            state = _read(self.path)
            candidate = state["candidates"].get(candidate_id)
            if candidate is None:
                raise KeyError(f"Unknown campaign candidate: {candidate_id}")
            if validated is not None:
                candidate["validated_acceptance"] = validated
            if diversity is not None:
                candidate["diversity_acceptance"] = diversity
            if reason:
                candidate["decision_reason"] = reason
            state["updated_at"] = utc_now()
            state["events"].append({
                "type": "candidate_acceptance", "candidate_id": candidate_id,
                "validated": validated, "diversity": diversity, "reason": reason,
                "at": state["updated_at"],
            })
            _atomic_json(self.path, state)

    def terminate(self, reason: str, *, detail: str = "") -> dict[str, Any]:
        if reason not in TERMINATION_REASONS:
            raise ValueError(f"Unknown campaign termination reason: {reason}")
        with _locked(self.path):
            state = _read(self.path)
            existing = state.get("termination_reason")
            if existing and existing != reason:
                return state
            state["status"] = "completed" if reason in {
                "target_reached", "budget_exhausted", "runtime_exhausted",
            } else reason
            state["termination_reason"] = reason
            state["termination_detail"] = detail
            state["terminated_at"] = utc_now()
            state["updated_at"] = state["terminated_at"]
            state["events"].append({
                "type": "campaign_terminated", "reason": reason,
                "detail": detail, "at": state["updated_at"],
            })
            _atomic_json(self.path, state)
            return state
