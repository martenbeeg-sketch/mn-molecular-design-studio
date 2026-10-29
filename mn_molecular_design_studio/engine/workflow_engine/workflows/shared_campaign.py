"""Durable coordination shared by graph-level generator campaigns.

Native generation and validation remain owned by each generator's ordinary
campaign controller. This file only reserves shared budgets and aggregates
candidate identity and downstream acceptance across those controllers.
"""
from __future__ import annotations

import copy
import fcntl
import hashlib
import json
import os
import tempfile
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ACTIVE_ATTEMPT_STATES = {"reserved", "authorized", "submitting", "queued", "preparing", "running"}
TERMINAL_ATTEMPT_STATES = {"completed", "failed", "interrupted", "cancelled", "not_started"}


def sequence_identity(candidate: dict[str, Any]) -> str:
    sequence = "".join(str(candidate.get("binder_sequence") or "").split()).upper()
    if not sequence or any(letter not in "ACDEFGHIKLMNPQRSTVWYBXZJUO" for letter in sequence):
        raise ValueError("Shared protein-design objectives require a valid binder_sequence for every candidate.")
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def _read(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
        return value if isinstance(value, dict) else {}
    except (OSError, ValueError):
        return {}


def _atomic_write(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(value, handle, indent=2, sort_keys=True, default=str)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass


@contextmanager
def _locked(path: Path):
    lock_path = path.with_suffix(path.suffix + ".lock")
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("a+b") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


class SharedCampaignLedger:
    """Atomic graph-level target, generation budgets, and sequence registry."""

    def __init__(self, path: Path):
        self.path = Path(path).expanduser().resolve()

    @classmethod
    def create(cls, path: Path, config: dict[str, Any]) -> "SharedCampaignLedger":
        ledger = cls(path)
        with _locked(ledger.path):
            state = _read(ledger.path)
            if state:
                if state.get("config") != config:
                    raise ValueError("A shared campaign ledger already exists with different graph settings.")
                return ledger
            now = _now()
            _atomic_write(ledger.path, {
                "schema_version": 1,
                "config": copy.deepcopy(config),
                "status": "running",
                "termination_reason": None,
                "created_at": now,
                "updated_at": now,
                "attempts": {},
                "campaigns": {},
                "sequences": {},
                "participants": {},
                "events": [{"type": "shared_campaign_created", "at": now}],
            })
        return ledger

    def snapshot(self) -> dict[str, Any]:
        with _locked(self.path):
            return _read(self.path)

    def reserve_attempt(
        self, campaign_id: str, attempt_id: str, generator_node: str, generator_id: str,
    ) -> bool:
        key = f"{campaign_id}:{attempt_id}"
        with _locked(self.path):
            state = _read(self.path)
            if state.get("status") != "running":
                return False
            attempts = state["attempts"]
            if key in attempts:
                return attempts[key].get("status") != "released"
            config = state["config"]
            budget = int(config["max_attempts"])
            reserved = [a for a in attempts.values() if a.get("status") != "released"]
            active = [a for a in attempts.values() if a.get("status") in ACTIVE_ATTEMPT_STATES]
            validations = [v for v in state.get("validations", {}).values() if v.get("status") in ACTIVE_ATTEMPT_STATES]
            limits = config.get("generator_limits", {}).get(generator_node, {})
            generator_attempts = [a for a in reserved if a.get("generator_node") == generator_node]
            generator_active = [a for a in active if a.get("generator_node") == generator_node]
            max_generator_attempts = int(limits.get("max_attempts", budget))
            max_generator_active = int(limits.get("max_concurrent_generation_jobs", config["max_concurrent_generation_jobs"]))
            compute_limit = config.get("max_compute_seconds")
            reported_compute = sum(float(row.get("compute_seconds", 0) or 0) for row in state["campaigns"].values())
            if (len(reserved) >= budget
                    or len(active) + len(validations) >= int(config["max_concurrent_jobs"])
                    or len(generator_attempts) >= max_generator_attempts
                    or len(generator_active) >= max_generator_active
                    or compute_limit is not None and reported_compute >= float(compute_limit)
                    or self._accepted_count(state) >= int(config["accepted_target"])
            ):
                self._refresh_terminal(state)
                state["updated_at"] = _now()
                _atomic_write(self.path, state)
                return False
            now = _now()
            attempts[key] = {
                "campaign_id": campaign_id, "attempt_id": attempt_id,
                "generator_node": generator_node, "generator_id": generator_id,
                "status": "reserved", "reserved_at": now,
            }
            state["events"].append({"type": "shared_attempt_reserved", "attempt_key": key,
                                    "generator_node": generator_node, "at": now})
            state["updated_at"] = now
            _atomic_write(self.path, state)
            return True

    def release_attempt(self, campaign_id: str, attempt_id: str, reason: str) -> None:
        key = f"{campaign_id}:{attempt_id}"
        with _locked(self.path):
            state = _read(self.path)
            attempt = state.get("attempts", {}).get(key)
            if attempt and attempt.get("status") == "reserved":
                attempt.update(status="released", release_reason=reason, completed_at=_now())
                state["events"].append({"type": "shared_attempt_released", "attempt_key": key, "at": _now()})
                state["updated_at"] = _now()
                self._refresh_terminal(state)
                _atomic_write(self.path, state)

    def sync_campaign(
        self, campaign_id: str, generator_node: str, local_state: dict[str, Any],
        compute_seconds: float = 0.0,
    ) -> dict[str, Any]:
        with _locked(self.path):
            state = _read(self.path)
            campaign = state.setdefault("campaigns", {}).setdefault(campaign_id, {
                "generator_node": generator_node, "generator_id": local_state.get("config", {}).get("generator_id"),
                "compute_seconds": 0.0,
            })
            campaign["compute_seconds"] = max(0.0, float(compute_seconds))
            local_attempts = local_state.get("attempts", {})
            for attempt_key, attempt in list(state.get("attempts", {}).items()):
                if attempt.get("campaign_id") != campaign_id:
                    continue
                local = local_attempts.get(attempt.get("attempt_id"))
                if local is None:
                    # The controller lock is held by the caller, so a reservation
                    # without a local authorization is an interrupted pre-launch step.
                    if attempt.get("status") == "reserved":
                        attempt.update(status="released", release_reason="No local attempt authorization was recovered.", completed_at=_now())
                    continue
                attempt["status"] = str(local.get("status") or "authorized")
                attempt["child_run_dir"] = local.get("child_run_dir")
                attempt["updated_at"] = _now()
            local_candidates = local_state.get("candidates", {})
            active_validations = 0
            pending_validations = 0
            local_validation_keys = set()
            for candidate_id, candidate in local_candidates.items():
                stages = candidate.get("validation_stages") or {}
                for stage_index, stage in stages.items():
                    key = f"{campaign_id}:{candidate_id}:{stage_index}"
                    local_validation_keys.add(key)
                    status = str(stage.get("status") or "pending")
                    reservation = state.setdefault("validations", {}).get(key)
                    if reservation:
                        reservation["status"] = status
                        reservation["run_dir"] = stage.get("run_dir")
                        reservation["updated_at"] = _now()
                    if status in ACTIVE_ATTEMPT_STATES:
                        active_validations += 1
                shared = candidate.get("shared_objective") or {}
                if (candidate.get("native_acceptance") in {"accepted", "unavailable"}
                        and candidate.get("validated_acceptance") not in {"accepted", "rejected", "failed"}
                        and shared.get("role") not in {"duplicate", "native_rejected", "not_evaluated"}):
                    pending_validations += 1
            for key, validation in state.get("validations", {}).items():
                if validation.get("campaign_id") != campaign_id:
                    continue
                if key not in local_validation_keys and validation.get("status") == "reserved":
                    validation.update(status="released", updated_at=_now())
            campaign["active_validations"] = active_validations
            campaign["pending_validations"] = pending_validations
            state["updated_at"] = _now()
            self._refresh_terminal(state)
            _atomic_write(self.path, state)
            return copy.deepcopy(state)

    def claim_candidate(
        self, campaign_id: str, generator_node: str, candidate: dict[str, Any],
        *, native_acceptance: str,
    ) -> dict[str, Any]:
        candidate_id = str(candidate.get("candidate_id") or "")
        if not candidate_id:
            raise ValueError("Shared-campaign candidates need stable candidate IDs.")
        identity = sequence_identity(candidate)
        participant_id = f"{campaign_id}:{candidate_id}"
        with _locked(self.path):
            state = _read(self.path)
            participants = state.setdefault("participants", {})
            existing_participant = participants.get(participant_id)
            sequence = state.setdefault("sequences", {}).get(identity)
            if existing_participant:
                if (state.get("status") == "running"
                        and existing_participant.get("role") == "duplicate" and sequence
                        and sequence.get("status") == "failed" and native_acceptance != "rejected"):
                    previous_owner = sequence.get("owner")
                    sequence.update(status="pending", owner=participant_id, campaign_id=campaign_id,
                                    candidate_id=candidate_id, generator_node=generator_node,
                                    previous_owner=previous_owner, retry_claimed_at=_now())
                    existing_participant.update(
                        role="owner", status="pending", decision="pending", duplicate_of=None,
                        owner=participant_id, retry_claimed_at=_now(),
                    )
                    state["events"].append({"type": "shared_candidate_retry_claimed",
                                            "participant_id": participant_id, "at": _now()})
                    self._atomic_refresh_write(state)
                    return copy.deepcopy(existing_participant)
                return copy.deepcopy(existing_participant)
            now = _now()
            if native_acceptance != "rejected" and state.get("status") != "running":
                row = {
                    "participant_id": participant_id, "sequence_identity": identity,
                    "campaign_id": campaign_id, "generator_node": generator_node,
                    "candidate_id": candidate_id, "role": "not_evaluated",
                    "status": "not_evaluated", "decision": state.get("termination_reason"),
                    "created_at": now,
                }
                participants[participant_id] = row
                state["updated_at"] = now
                self._atomic_refresh_write(state)
                return copy.deepcopy(row)
            if native_acceptance == "rejected":
                if sequence is None:
                    sequence = {"sequence_identity": identity, "status": "native_rejected", "created_at": now}
                    state["sequences"][identity] = sequence
                role = "native_rejected"
            elif sequence is None and state.get("status") != "running":
                row = {
                    "participant_id": participant_id, "sequence_identity": identity,
                    "campaign_id": campaign_id, "generator_node": generator_node,
                    "candidate_id": candidate_id, "role": "not_evaluated",
                    "status": "not_evaluated", "decision": state.get("termination_reason"),
                    "created_at": now,
                }
                participants[participant_id] = row
                state["updated_at"] = now
                self._atomic_refresh_write(state)
                return copy.deepcopy(row)
            elif sequence is None or sequence.get("status") == "native_rejected":
                sequence = {
                    "sequence_identity": identity, "status": "pending",
                    "owner": participant_id, "created_at": now,
                    "generator_node": generator_node, "campaign_id": campaign_id,
                    "candidate_id": candidate_id,
                }
                state["sequences"][identity] = sequence
                role = "owner"
            elif sequence.get("status") == "failed" and native_acceptance != "rejected":
                previous_owner = sequence.get("owner")
                sequence.update(status="pending", owner=participant_id, campaign_id=campaign_id,
                                candidate_id=candidate_id, generator_node=generator_node,
                                previous_owner=previous_owner, retry_claimed_at=now)
                role = "owner"
            else:
                role = "duplicate"
            row = {
                "participant_id": participant_id, "sequence_identity": identity,
                "campaign_id": campaign_id, "generator_node": generator_node,
                "candidate_id": candidate_id, "role": role,
                "status": "native_rejected" if role == "native_rejected" else "pending" if role == "owner" else "duplicate",
                "owner": sequence.get("owner"), "decision": sequence.get("status"),
                "duplicate_of": sequence.get("candidate_id") if role == "duplicate" else None,
                "created_at": now,
            }
            participants[participant_id] = row
            state["events"].append({"type": "shared_candidate_registered", "participant_id": participant_id,
                                    "sequence_identity": identity, "role": role, "at": now})
            state["updated_at"] = now
            self._refresh_terminal(state)
            _atomic_write(self.path, state)
            return copy.deepcopy(row)

    def reserve_validation(
        self, campaign_id: str, generator_node: str, candidate_id: str, stage_index: int,
    ) -> bool:
        key = f"{campaign_id}:{candidate_id}:{stage_index}"
        with _locked(self.path):
            state = _read(self.path)
            if state.get("status") != "running":
                return False
            validations = state.setdefault("validations", {})
            if key in validations:
                existing = validations[key]
                if existing.get("status") not in {"released", "reserved", "submitting"}:
                    return False
                existing.update(
                    status="reserved", reserved_at=_now(),
                    retry_count=int(existing.get("retry_count", 0)) + 1,
                )
                state["events"].append({"type": "shared_validation_rereserved", "validation_key": key, "at": _now()})
                state["updated_at"] = _now()
                _atomic_write(self.path, state)
                return True
            active = [row for row in validations.values() if row.get("status") in ACTIVE_ATTEMPT_STATES]
            all_active = sum(row.get("status") in ACTIVE_ATTEMPT_STATES for row in state.get("attempts", {}).values())
            compute_limit = state["config"].get("max_compute_seconds")
            reported_compute = sum(float(row.get("compute_seconds", 0) or 0) for row in state.get("campaigns", {}).values())
            if (all_active + len(active) >= int(state["config"]["max_concurrent_jobs"])
                    or len(active) >= int(state["config"].get("max_concurrent_validations", 1))):
                return False
            if compute_limit is not None and reported_compute >= float(compute_limit):
                self._refresh_terminal(state)
                state["updated_at"] = _now()
                _atomic_write(self.path, state)
                return False
            now = _now()
            validations[key] = {"campaign_id": campaign_id, "generator_node": generator_node,
                                "candidate_id": candidate_id, "stage_index": stage_index,
                                "status": "reserved", "reserved_at": now}
            state["events"].append({"type": "shared_validation_reserved", "validation_key": key, "at": now})
            state["updated_at"] = now
            _atomic_write(self.path, state)
            return True

    def release_validation(self, campaign_id: str, candidate_id: str, stage_index: int) -> None:
        key = f"{campaign_id}:{candidate_id}:{stage_index}"
        with _locked(self.path):
            state = _read(self.path)
            validation = state.get("validations", {}).get(key)
            if validation and validation.get("status") == "reserved":
                validation.update(status="released", updated_at=_now())
                state["updated_at"] = _now()
                _atomic_write(self.path, state)

    def _atomic_refresh_write(self, state: dict[str, Any]) -> None:
        self._refresh_terminal(state)
        state["updated_at"] = _now()
        _atomic_write(self.path, state)

    def finish_candidate(self, campaign_id: str, candidate_id: str, status: str, reason: str = "") -> None:
        if status not in {"accepted", "rejected", "failed"}:
            raise ValueError("Shared downstream decision must be accepted, rejected, or failed.")
        participant_id = f"{campaign_id}:{candidate_id}"
        with _locked(self.path):
            state = _read(self.path)
            participant = state.get("participants", {}).get(participant_id)
            if not participant or participant.get("role") != "owner":
                return
            sequence = state.get("sequences", {}).get(participant["sequence_identity"])
            if not sequence or sequence.get("owner") != participant_id:
                return
            sequence.update(status=status, decision_reason=str(reason or ""), completed_at=_now())
            participant.update(status=status, decision=status, decision_reason=str(reason or ""), completed_at=_now())
            state["events"].append({"type": "shared_candidate_decided", "participant_id": participant_id,
                                    "decision": status, "at": _now()})
            state["updated_at"] = _now()
            self._refresh_terminal(state)
            _atomic_write(self.path, state)

    def participant(self, campaign_id: str, candidate_id: str) -> dict[str, Any] | None:
        state = self.snapshot()
        participant = state.get("participants", {}).get(f"{campaign_id}:{candidate_id}")
        if participant is None:
            return None
        participant = copy.deepcopy(participant)
        sequence = state.get("sequences", {}).get(participant.get("sequence_identity"), {})
        participant["decision"] = sequence.get("status", participant.get("decision"))
        participant["decision_reason"] = sequence.get("decision_reason", participant.get("decision_reason", ""))
        return participant

    def summary(self) -> dict[str, Any]:
        state = self.snapshot()
        attempts = [row for row in state.get("attempts", {}).values() if row.get("status") != "released"]
        per_generator: dict[str, dict[str, int]] = {}
        for row in attempts:
            counters = per_generator.setdefault(row["generator_node"], {"attempts": 0, "active_attempts": 0})
            counters["attempts"] += 1
            counters["active_attempts"] += row.get("status") in ACTIVE_ATTEMPT_STATES
        for row in state.get("participants", {}).values():
            counters = per_generator.setdefault(row["generator_node"], {"attempts": 0, "active_attempts": 0})
            counters[f"candidates_{row.get('status', 'pending')}"] = counters.get(f"candidates_{row.get('status', 'pending')}", 0) + 1
        sequences = list(state.get("sequences", {}).values())
        return {
            "status": state.get("status", "running"),
            "termination_reason": state.get("termination_reason"),
            "accepted_target": state["config"]["accepted_target"],
            "unique_accepted": sum(row.get("status") == "accepted" for row in sequences),
            "unique_rejected": sum(row.get("status") == "rejected" for row in sequences),
            "unique_validation_failed": sum(row.get("status") == "failed" for row in sequences),
            "native_rejected_unique": sum(row.get("status") == "native_rejected" for row in sequences),
            "duplicate_candidates": sum(row.get("role") == "duplicate" for row in state.get("participants", {}).values()),
            "attempts_authorized": len(attempts),
            "attempts_active": sum(row.get("status") in ACTIVE_ATTEMPT_STATES for row in attempts),
            "validations_active": sum(row.get("status") in ACTIVE_ATTEMPT_STATES for row in state.get("validations", {}).values()),
            "max_attempts": state["config"]["max_attempts"],
            "max_concurrent_generation_jobs": state["config"]["max_concurrent_generation_jobs"],
            "compute_seconds": sum(float(row.get("compute_seconds", 0) or 0) for row in state.get("campaigns", {}).values()),
            "per_generator": per_generator,
            "updated_at": state.get("updated_at"),
        }

    def branch_budget_exhausted(self, generator_node: str) -> bool:
        state = self.snapshot()
        config = state.get("config", {})
        maximum = int(config.get("generator_limits", {}).get(generator_node, {}).get(
            "max_attempts", config.get("max_attempts", 1),
        ))
        rows = [row for row in state.get("attempts", {}).values()
                if row.get("generator_node") == generator_node and row.get("status") != "released"]
        return (len(rows) >= maximum
                and not any(row.get("status") in ACTIVE_ATTEMPT_STATES for row in rows))

    def _accepted_count(self, state: dict[str, Any]) -> int:
        return sum(row.get("status") == "accepted" for row in state.get("sequences", {}).values())

    def _refresh_terminal(self, state: dict[str, Any]) -> None:
        if state.get("status") != "running":
            return
        if self._accepted_count(state) >= int(state["config"]["accepted_target"]):
            state.update(status="target_reached", termination_reason="target_reached", completed_at=_now())
            state["events"].append({"type": "shared_target_reached", "at": _now()})
            return
        attempts = [row for row in state.get("attempts", {}).values() if row.get("status") != "released"]
        active = [row for row in attempts if row.get("status") in ACTIVE_ATTEMPT_STATES]
        pending = (any(row.get("status") == "pending" for row in state.get("sequences", {}).values())
                   or any(int(row.get("pending_validations", 0)) for row in state.get("campaigns", {}).values()))
        config = state["config"]
        compute = sum(float(row.get("compute_seconds", 0) or 0) for row in state.get("campaigns", {}).values())
        compute_exhausted = config.get("max_compute_seconds") is not None and compute >= float(config["max_compute_seconds"])
        attempt_exhausted = len(attempts) >= int(config["max_attempts"])
        limits = config.get("generator_limits", {})
        branch_exhausted = all(
            sum(a.get("generator_node") == node and a.get("status") != "released" for a in attempts)
            >= int(limits.get(node, {}).get("max_attempts", config["max_attempts"]))
            for node in config["generator_nodes"]
        )
        active_validations = any(int(row.get("active_validations", 0)) for row in state.get("campaigns", {}).values())
        if compute_exhausted and not active and not active_validations:
            for row in state.get("sequences", {}).values():
                if row.get("status") == "pending":
                    row.update(status="unresolved_compute_budget", completed_at=_now())
            reason = "compute_exhausted"
            state.update(status=reason, termination_reason=reason, completed_at=_now())
            state["events"].append({"type": "shared_budget_exhausted", "reason": reason, "at": _now()})
            return
        if (attempt_exhausted or branch_exhausted) and not active and not active_validations and not pending:
            reason = "compute_exhausted" if compute_exhausted else "budget_exhausted"
            state.update(status=reason, termination_reason=reason, completed_at=_now())
            state["events"].append({"type": "shared_budget_exhausted", "reason": reason, "at": _now()})


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()
