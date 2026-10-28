"""Durable acceptance accounting for distributed BindCraft2 campaigns.

The ledger lives beside its parent workflow run and uses an OS file lock plus
atomic replacement. It is deliberately not a second scheduler or database.
"""
from __future__ import annotations

import fcntl
import hashlib
import json
import math
import os
import tempfile
import time
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
    "target_reached", "budget_exhausted", "runtime_exhausted", "compute_exhausted", "cancelled", "failed",
}


class CampaignControllerAlreadyRunning(RuntimeError):
    """A duplicate controller worker should exit without changing parent state."""


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
        directory_fd = os.open(path.parent, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
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
    """Validate engine capability, global counters, and execution bounds."""
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
    if runtime is not None and (
        type(runtime) not in {int, float} or not math.isfinite(float(runtime)) or runtime <= 0
    ):
        raise ValueError("max_runtime_seconds must be positive when provided.")
    compute = config.get("max_compute_seconds")
    if compute is not None and (
        type(compute) not in {int, float} or not math.isfinite(float(compute)) or compute <= 0
    ):
        raise ValueError("max_compute_seconds must be positive when provided.")
    max_retries = config.get("max_retries", 1)
    if type(max_retries) is not int or max_retries < 0 or max_retries > 5:
        raise ValueError("max_retries must be an integer from 0 to 5.")
    poll_interval = config.get("poll_interval_seconds", 2.0)
    if type(poll_interval) not in {int, float} or not 0.25 <= float(poll_interval) <= 30:
        raise ValueError("poll_interval_seconds must be between 0.25 and 30.")
    placement = str(config.get("workstation_selection", "local")).strip().lower()
    if placement not in {"local", "auto", "eligible"}:
        raise ValueError("workstation_selection must be local, auto, or eligible.")
    gpu_selection = str(config.get("gpu_selection", "auto")).strip().lower()
    if gpu_selection not in {"auto", "eligible", "all"}:
        gpu_parts = [part.strip() for part in gpu_selection.removeprefix("device=").split(",")]
        if not gpu_parts or any(not part.isdigit() for part in gpu_parts):
            raise ValueError("gpu_selection must be auto, eligible, all, or a comma-separated list of GPU IDs.")
    resources = config.get("resources", {})
    if not isinstance(resources, dict):
        raise ValueError("resources must be an object.")
    unknown_resources = set(resources) - {"cpu_cores", "min_vram_gb", "ram_gb", "scratch_gb"}
    if unknown_resources:
        raise ValueError(f"Unknown resource field(s): {', '.join(sorted(unknown_resources))}.")
    if "cpu_cores" in resources and (type(resources["cpu_cores"]) is not int or resources["cpu_cores"] < 1):
        raise ValueError("resources.cpu_cores must be a positive integer.")
    for key in ("min_vram_gb", "ram_gb", "scratch_gb"):
        value = resources.get(key, 0)
        if type(value) not in {int, float} or not math.isfinite(float(value)) or value < 0:
            raise ValueError(f"resources.{key} must be finite and nonnegative.")
    requirements = config.get("validation_stages", [])
    if not isinstance(requirements, list) or any(not isinstance(item, dict) for item in requirements):
        raise ValueError("validation_stages must be a list of workflow objects.")
    allowed_validation_workflows = {"refolding.optimized_prediction", "analysis.ipsae"}
    has_validation_criteria = False
    for index, stage in enumerate(requirements, start=1):
        workflow = str(stage.get("workflow") or "")
        if workflow not in allowed_validation_workflows:
            raise ValueError(
                f"validation_stages[{index - 1}].workflow must be one of "
                f"{', '.join(sorted(allowed_validation_workflows))}."
            )
        if not isinstance(stage.get("parameters", {}), dict):
            raise ValueError(f"validation_stages[{index - 1}].parameters must be an object.")
        required_metrics = stage.get("required_metrics", [])
        if not isinstance(required_metrics, list) or any(not isinstance(item, str) or not item.strip() for item in required_metrics):
            raise ValueError(f"validation_stages[{index - 1}].required_metrics must be a list of metric names.")
        rules = stage.get("rules", [])
        if not isinstance(rules, list):
            raise ValueError(f"validation_stages[{index - 1}].rules must be a list.")
        for rule in rules:
            if (
                not isinstance(rule, dict) or not str(rule.get("field") or "").strip()
                or rule.get("operator") not in {">=", "<=", ">", "<"}
                or type(rule.get("value")) not in {int, float}
                or not math.isfinite(float(rule["value"]))
            ):
                raise ValueError(f"validation_stages[{index - 1}] has an invalid scientific rule.")
        has_validation_criteria = has_validation_criteria or bool(required_metrics or rules)
        stage_resources = stage.get("resources", {})
        if not isinstance(stage_resources, dict) or set(stage_resources) - {"cpu_cores", "gpu_device", "min_vram_gb"}:
            raise ValueError(f"validation_stages[{index - 1}].resources has unsupported fields.")
        if stage_resources.get("cpu_cores") is not None and (
            type(stage_resources["cpu_cores"]) is not int or stage_resources["cpu_cores"] < 1
        ):
            raise ValueError(f"validation_stages[{index - 1}].resources.cpu_cores must be positive.")
        if stage_resources.get("min_vram_gb") is not None and (
            type(stage_resources["min_vram_gb"]) not in {int, float}
            or not math.isfinite(float(stage_resources["min_vram_gb"])) or stage_resources["min_vram_gb"] < 0
        ):
            raise ValueError(f"validation_stages[{index - 1}].resources.min_vram_gb must be finite and nonnegative.")
        if stage["workflow"] == "analysis.ipsae" and (
            index < 2 or requirements[index - 2].get("workflow") != "refolding.optimized_prediction"
        ):
            raise ValueError("analysis.ipsae must follow a refolding.optimized_prediction stage.")
    generator_id = str(config.get("generator_id") or "bindcraft2").strip().lower()
    from mn_molecular_design_studio.engine.workflow_engine.workflows.generator_capabilities import (
        generator_capabilities,
        validate_stopping_policy,
    )

    capability = generator_capabilities(generator_id)
    validate_stopping_policy(
        generator_id, acceptance_level, requirements,
        has_validation_criteria=has_validation_criteria,
    )
    if acceptance_level == "validated" and (not requirements or not has_validation_criteria):
        raise ValueError("Validated acceptance requires downstream stages with explicit measurements or scientific rules.")
    if capability.minimum_vram_gb is not None:
        requested_minimum = float(resources.get("min_vram_gb", 0))
        if requested_minimum and requested_minimum < capability.minimum_vram_gb:
            raise ValueError(
                f"{capability.label} requires at least {capability.minimum_vram_gb:g} GiB of GPU memory for campaign attempts."
            )
        resources = {**resources, "min_vram_gb": max(requested_minimum, capability.minimum_vram_gb)}
    validation_concurrency = config.get("max_concurrent_validations", 2)
    if type(validation_concurrency) is not int or validation_concurrency < 1:
        raise ValueError("max_concurrent_validations must be a positive integer.")
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
    structural_distance = diversity.get("binder_ca_rmsd_min_angstrom")
    if structural_distance is not None and (
        type(structural_distance) not in {int, float} or not math.isfinite(float(structural_distance))
        or structural_distance <= 0
    ):
        raise ValueError("diversity.binder_ca_rmsd_min_angstrom must be positive when provided.")
    interface_field = diversity.get("interface_cluster_field")
    if interface_field is not None and (not isinstance(interface_field, str) or not interface_field.strip()):
        raise ValueError("diversity.interface_cluster_field must be a nonempty measurement selector.")
    if interface_field is not None and not interface_field.startswith("evidence/"):
        raise ValueError("diversity.interface_cluster_field must use a stage-qualified evidence/<stage>/<tool>/<metric> selector.")
    return {
        **config,
        "generator_id": generator_id,
        "acceptance_level": acceptance_level,
        "max_attempts": maximum,
        "accepted_target": target,
        "batch_size": batch,
        "max_concurrent_batches": concurrency,
        "seed": seed,
        "validation_stages": requirements,
        "diversity": diversity,
        "max_retries": max_retries,
        "poll_interval_seconds": float(poll_interval),
        "workstation_selection": placement,
        "gpu_selection": gpu_selection,
        "resources": resources,
        "max_concurrent_validations": validation_concurrency,
    }


def _state_accepted_count(state: dict[str, Any]) -> int:
    acceptance_field = {
        "native": "native_acceptance",
        "validated": "validated_acceptance",
        "diversity_qualified": "diversity_acceptance",
    }[state["config"]["acceptance_level"]]
    return sum(
        candidate.get(acceptance_field) == "accepted"
        for candidate in state.get("candidates", {}).values()
    )


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
                "generation_budget_exhausted": False,
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

    @contextmanager
    def controller_lock(self) -> Iterator[None]:
        lock_path = self.path.with_suffix(self.path.suffix + ".controller.lock")
        lock_path.parent.mkdir(parents=True, exist_ok=True)
        with lock_path.open("a+b") as handle:
            try:
                fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as exc:
                raise CampaignControllerAlreadyRunning("Another controller process already owns this campaign.") from exc
            try:
                yield
            finally:
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


    def counters(self) -> dict[str, int]:
        state = self.snapshot()
        attempts = list(state["attempts"].values())
        candidates = list(state["candidates"].values())
        return {
            "attempts_authorized": len(attempts),
            "attempts_started": sum(bool(item.get("first_started_at")) for item in attempts),
            "attempts_completed": sum(item.get("status") == "completed" for item in attempts),
            "native_rejected_trajectories": sum(
                max(
                    0,
                    int((item.get("native_state") or {}).get("trajectories", 0))
                    - int((item.get("native_state") or {}).get("accepted", 0)),
                )
                for item in attempts
                if item.get("status") == "completed"
                and isinstance(item.get("native_state"), dict)
                and (item.get("native_state") or {}).get("authoritative") is not False
            ),
            "native_acceptance_unavailable_attempts": sum(
                item.get("status") == "completed"
                and isinstance(item.get("native_state"), dict)
                and (item.get("native_state") or {}).get("authoritative") is False
                for item in attempts
            ),
            "native_accepted_candidates": sum(item.get("native_acceptance") == "accepted" for item in candidates),
            "native_acceptance_unavailable_candidates": sum(item.get("native_acceptance") == "unavailable" for item in candidates),
            "validated_accepted_candidates": sum(item.get("validated_acceptance") == "accepted" for item in candidates),
            "diversity_qualified_accepted_candidates": sum(item.get("diversity_acceptance") == "accepted" for item in candidates),
            "rejected_candidates": sum(
                item.get("native_acceptance") == "rejected"
                or item.get("validated_acceptance") == "rejected"
                or item.get("diversity_acceptance") == "rejected"
                for item in candidates
            ),
            "failed_interrupted_attempts": sum(
                len(item.get("retry_history") or []) + int(item.get("status") in {"failed", "interrupted"})
                for item in attempts
            ),
        }

    def update_candidate(self, candidate_id: str, **details: Any) -> None:
        with _locked(self.path):
            state = _read(self.path)
            candidate = state["candidates"].get(candidate_id)
            if candidate is None:
                raise KeyError(f"Unknown campaign candidate: {candidate_id}")
            if all(candidate.get(key) == value for key, value in details.items()):
                return
            candidate.update(details)
            state["updated_at"] = utc_now()
            state["events"].append({
                "type": "candidate_updated", "candidate_id": candidate_id,
                "fields": sorted(details), "at": state["updated_at"],
            })
            _atomic_json(self.path, state)

    def record_compute(self, attempt_id: str, run_id: str, gpu_seconds: float) -> bool:
        if not math.isfinite(float(gpu_seconds)) or gpu_seconds < 0:
            raise ValueError("GPU compute duration must be finite and nonnegative.")
        with _locked(self.path):
            state = _read(self.path)
            attempt = state["attempts"].get(attempt_id)
            if attempt is None:
                raise KeyError(f"Unknown campaign attempt: {attempt_id}")
            runs = attempt.setdefault("compute_runs", {})
            if run_id in runs:
                return False
            runs[run_id] = float(gpu_seconds)
            state["updated_at"] = utc_now()
            state["events"].append({
                "type": "compute_recorded", "attempt_id": attempt_id,
                "run_id": run_id, "gpu_seconds": float(gpu_seconds), "at": state["updated_at"],
            })
            _atomic_json(self.path, state)
            return True

    def gpu_compute_seconds(self) -> float:
        return sum(
            float(seconds)
            for attempt in self.snapshot()["attempts"].values()
            for seconds in (attempt.get("compute_runs") or {}).values()
        )

    def mark_generation_budget_exhausted(self) -> dict[str, Any]:
        """Persist that all generation attempts have been authorized, even if validation is pending."""
        with _locked(self.path):
            state = _read(self.path)
            if len(state["attempts"]) < int(state["config"]["max_attempts"]):
                return state
            if state.get("generation_budget_exhausted"):
                return state
            now = utc_now()
            state["generation_budget_exhausted"] = True
            state["generation_budget_exhausted_at"] = now
            state["updated_at"] = now
            state["events"].append({
                "type": "generation_budget_exhausted", "at": now,
            })
            _atomic_json(self.path, state)
            return state

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
                attempt.setdefault("first_started_at", attempt["started_at"])
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
        gpu_seconds: float = 0.0,
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
            attempt.setdefault("compute_runs", {})
            if attempt.get("child_run_dir"):
                attempt["compute_runs"][Path(str(attempt["child_run_dir"])).name] = max(0.0, float(gpu_seconds))
            attempt["candidate_ids"] = []
            from mn_molecular_design_studio.engine.workflow_engine.workflows.generator_capabilities import (
                generator_capabilities,
            )

            capability = generator_capabilities(str(state["config"].get("generator_id") or "bindcraft2"))
            for candidate in candidates:
                candidate_id = str(candidate.get("candidate_id") or "")
                if not candidate_id:
                    raise ValueError("Verified campaign candidates need nonempty stable IDs.")
                if candidate_id in state["candidates"]:
                    continue
                if capability.native_acceptance_available:
                    native_acceptance = (
                        "accepted" if candidate.get("metrics", {}).get("native_pass_filters") is True else "rejected"
                    )
                else:
                    native_acceptance = "unavailable"
                state["candidates"][candidate_id] = {
                    "candidate_id": candidate_id,
                    "attempt_id": attempt_id,
                    "native_acceptance": native_acceptance,
                    "validated_acceptance": "pending",
                    "diversity_acceptance": "pending",
                    "candidate": dict(candidate),
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
            if validated == "accepted" and candidate.get("native_acceptance") not in {"accepted", "unavailable"}:
                raise ValueError("Validated acceptance requires a native-accepted candidate or an engine without native acceptance.")
            if diversity == "accepted" and candidate.get("native_acceptance") not in {"accepted", "unavailable"}:
                raise ValueError("Diversity acceptance requires a native-accepted candidate or an engine without native acceptance.")
            if (
                diversity == "accepted"
                and state["config"].get("validation_stages")
                and candidate.get("validated_acceptance") != "accepted"
            ):
                raise ValueError("Diversity acceptance requires configured validation stages to pass.")
            if (
                (validated is None or candidate.get("validated_acceptance") == validated)
                and (diversity is None or candidate.get("diversity_acceptance") == diversity)
                and (not reason or candidate.get("decision_reason") == reason)
            ):
                return
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
            promotes_budget_outcome = (
                existing == "budget_exhausted"
                and reason == "target_reached"
                and _state_accepted_count(state) >= int(state["config"]["accepted_target"])
            )
            if existing and existing != reason and not promotes_budget_outcome:
                return state
            if len(state["attempts"]) >= int(state["config"]["max_attempts"]):
                state["generation_budget_exhausted"] = True
            state["status"] = "completed" if reason in {
                "target_reached", "budget_exhausted", "runtime_exhausted", "compute_exhausted",
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


def _campaign_ledger_path(parent_run_dir: Path) -> Path:
    return Path(parent_run_dir) / "artifacts" / "distributed_campaign" / "ledger.json"


def _gpu_for_attempt(config: dict[str, Any], ordinal: int) -> str:
    selection = str(config.get("gpu_selection") or "auto").strip().lower().removeprefix("device=")
    if selection in {"all", "auto", "eligible"}:
        return "auto"
    ids = [part.strip() for part in selection.split(",") if part.strip()]
    return ids[(ordinal - 1) % len(ids)]


def _peer_campaign_mode(config: dict[str, Any]) -> bool:
    """Whether this campaign should submit eligible attempts through peer dispatch."""
    placement = str(config.get("workstation_selection") or "local").strip().lower()
    if placement == "local":
        return False
    from mn_molecular_design_studio.cluster import config as cluster_config

    mode = cluster_config.load().get("mode")
    if mode in {"peer", "coordinator"}:
        return True
    if placement == "eligible":
        raise ValueError("Eligible workstation execution requires an enrolled peer or coordinator.")
    return False


def _campaign_actor(parent_run_dir: Path) -> tuple[str, str, bool]:
    """Resolve the campaign owner/project from durable metadata or its API marker."""
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json

    metadata = read_json(Path(parent_run_dir) / "metadata.json")
    marker = read_json(Path(parent_run_dir) / "molecular_design_studio.json")
    owner_id = str(metadata.get("studio_owner_id") or marker.get("studio_owner_id") or "")
    project_id = str(metadata.get("studio_project_id") or metadata.get("project_id") or marker.get("project_id") or "")
    if not owner_id or not project_id:
        raise ValueError("Remote campaign attempts require a project-owned API submission.")
    development_mode = bool(metadata.get("studio_development_mode") or marker.get("studio_development_mode"))
    return owner_id, project_id, development_mode


def _submit_rfd4_peer_attempt(
    parent_run_dir: Path,
    ledger: CampaignLedger,
    attempt: dict[str, Any],
    child_run_dir: Path,
) -> dict[str, Any]:
    """Submit or recover one idempotent remote RFD4 attempt for its local proxy run."""
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json, write_json
    from mn_molecular_design_studio.cluster import candidate_batches

    parent_run_dir, child_run_dir = Path(parent_run_dir), Path(child_run_dir)
    state = ledger.snapshot()
    campaign_id = state["campaign_id"]
    peer = dict(attempt.get("cluster_submission") or {})
    peer.setdefault("manifest", str(child_run_dir / "artifacts" / "distributed_campaign" / "peer_submission.json"))
    peer.setdefault("transfers", {})
    peer.setdefault("submission_key", hashlib.sha256(
        f"{campaign_id}:{attempt['attempt_id']}:rfd4".encode()
    ).hexdigest())
    owner_id, project_id, development_mode = _campaign_actor(parent_run_dir)
    child_payload = read_json(child_run_dir / "input.json")
    parameters = dict(child_payload.get("params") or {})
    parameters["gpu_device"] = "auto"
    # Queue routing metadata is stored alongside workflow parameters in the
    # local child input, but it is not part of the public RFD4 request schema.
    parameters.pop("task_group", None)
    parameters.pop("distributed_campaign", None)
    metadata = read_json(child_run_dir / "metadata.json")
    resources = dict(metadata.get("resource_request") or {})
    # The scheduler's boolean GPU reservation marker is internal metadata;
    # the remote workflow API infers GPU use from the workflow and device.
    resources.pop("gpu", None)
    resources.update(gpu_device="auto")
    request = {
        "schema_version": 1,
        "workflow": "design.rfd4_proteina",
        "parameters": parameters,
        "resources": resources,
    }

    def persist() -> None:
        ledger.transition(
            attempt["attempt_id"], "queued", child_run_dir=str(child_run_dir),
            cluster_submission=peer, remote_submission_key=peer["submission_key"],
        )

    result = candidate_batches.submit(
        request, peer["manifest"], peer["submission_key"], {"placement": "auto"},
        owner_id=owner_id, project_id=project_id, state=peer, persist=persist,
        allow_development=development_mode,
    )
    peer["cluster_job_id"] = result["cluster_job_id"]
    peer["origin_node_id"] = result.get("origin_node_id")
    metadata.update(
        cluster_remote=True,
        cluster_job_id=peer["cluster_job_id"],
        remote_submission_key=peer["submission_key"],
        requested_workstation_selection=state["config"].get("workstation_selection"),
    )
    write_json(child_run_dir / "metadata.json", metadata)
    ledger.transition(
        attempt["attempt_id"], "queued", child_run_dir=str(child_run_dir),
        cluster_job_id=peer["cluster_job_id"], cluster_submission=peer,
        remote_submission_key=peer["submission_key"],
    )
    return peer


def _materialize_rfd4_peer_output(
    parent_run_dir: Path,
    ledger: CampaignLedger,
    attempt: dict[str, Any],
    child_run_dir: Path,
    peer: dict[str, Any],
    output: dict[str, Any],
    remote_job: dict[str, Any],
) -> bool:
    """Copy SHA-verified remote artifacts into the local attempt proxy."""
    import shutil

    from mn_molecular_design_studio.cluster import candidate_batches
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json, write_json

    owner_id, project_id, development_mode = _campaign_actor(parent_run_dir)
    collected = Path(peer["manifest"]).with_suffix(".collected")
    job_id = str(peer.get("cluster_job_id") or "")
    # These files are needed by the campaign's existing target/parameter
    # identity checks. Candidate structures and normalized rows are fetched by
    # candidate_batches.poll() from their own artifact references.
    request_relative = "artifacts/rfd4_proteina/request.json"
    if candidate_batches._artifact(
        peer, job_id, request_relative,
        candidate_batches.actor(owner_id, project_id, allow_development=development_mode), collected,
    ) is None:
        return False
    remote_request = read_json(collected / request_relative)
    target_value = remote_request.get("target_structure")
    if not target_value and isinstance(remote_request.get("parameters"), dict):
        target_value = remote_request["parameters"].get("target_structure")
    target_name = Path(str(target_value or "")).name
    if not target_name or target_name in {".", ".."}:
        raise RuntimeError("Remote RFD4 request is missing its staged target identity.")
    for relative in (
        f"artifacts/rfd4_proteina/input/{target_name}",
        "artifacts/rfd4_proteina/native_command.json",
        "artifacts/rfd4_proteina/effective_config.yaml",
        "artifacts/rfd4_proteina/prepared_inputs.json",
    ):
        if candidate_batches._artifact(
            peer, job_id, relative,
            candidate_batches.actor(owner_id, project_id, allow_development=development_mode), collected,
        ) is None:
            return False

    source_artifacts = collected / "artifacts"
    if not source_artifacts.is_dir():
        raise RuntimeError("Remote RFD4 completed without a collected artifact directory.")
    shutil.copytree(source_artifacts, child_run_dir / "artifacts", dirs_exist_ok=True)
    rows_path = Path(output["candidates_jsonl"])
    if not rows_path.is_file():
        raise RuntimeError("Remote RFD4 completed without normalized candidates.")

    def rebase(value: Any) -> Any:
        if isinstance(value, dict):
            return {key: rebase(item) for key, item in value.items()}
        if isinstance(value, list):
            return [rebase(item) for item in value]
        if isinstance(value, str) and (value == str(collected) or value.startswith(str(collected) + os.sep)):
            return str(child_run_dir) + value[len(str(collected)):]
        return value

    rows = [rebase(json.loads(line)) for line in rows_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    normalized = child_run_dir / "artifacts" / "normalized_candidates" / "candidates.jsonl"
    normalized.parent.mkdir(parents=True, exist_ok=True)
    normalized.write_text("".join(json.dumps(row, sort_keys=True) + "\n" for row in rows), encoding="utf-8")
    remote_metadata = read_json(collected / "metadata.json")
    remote_result = read_json(collected / "result.json")
    if remote_metadata.get("status") != "completed" or remote_result.get("success") is not True:
        raise RuntimeError("Remote completion disagrees with the transferred RFD4 result.")
    write_json(child_run_dir / "result.json", remote_result)
    from mn_molecular_design_studio.cache_provenance import file_digest
    native_command = read_json(child_run_dir / "artifacts" / "rfd4_proteina" / "native_command.json")
    effective_config_path = child_run_dir / "artifacts" / "rfd4_proteina" / "effective_config.yaml"
    allocation = remote_metadata.get("resource_allocation") if isinstance(remote_metadata.get("resource_allocation"), dict) else {}
    write_json(child_run_dir / "artifacts" / "distributed_campaign" / "remote_execution.json", {
        "cluster_job_id": job_id,
        "worker_node_id": remote_job.get("node_id"),
        "worker_name": remote_job.get("node_name") or remote_job.get("name"),
        "status": remote_metadata.get("status"),
        "resource_allocation": remote_metadata.get("resource_allocation") or {},
        "software_environment": remote_metadata.get("software_environment") or {},
        "checkpoint_revision": remote_result.get("metrics", {}).get("checkpoint_revision"),
        "image": remote_result.get("metrics", {}).get("image"),
        "started_at": remote_metadata.get("worker_started_at") or remote_metadata.get("started_at"),
        "created_at": remote_metadata.get("created_at"),
        "queue_wait_seconds": remote_metadata.get("queue_wait_seconds"),
        "allocated_at": allocation.get("allocated_at"),
        "completed_at": remote_metadata.get("completed_at"),
        "native_invocation": native_command,
        "effective_config_sha256": file_digest(effective_config_path) if effective_config_path.is_file() else None,
        "input_transfer_bytes": peer.get("input_transfer_bytes", 0),
        "output_transfer_bytes": peer.get("output_transfer_bytes", 0),
        "transfers": peer.get("transfers", {}),
    })
    local_metadata = read_json(child_run_dir / "metadata.json")
    local_metadata.update(
        status="completed",
        cluster_job_id=job_id,
        execution_node_id=remote_job.get("node_id"),
        scheduler_worker_id=remote_metadata.get("scheduler_worker_id") or remote_metadata.get("worker_id") or remote_job.get("node_id"),
        resource_allocation=allocation,
        worker_started_at=remote_metadata.get("worker_started_at") or allocation.get("allocated_at"),
        completed_at=remote_metadata.get("completed_at"),
        remote_execution_node_id=remote_job.get("node_id"),
        remote_status="completed",
    )
    write_json(child_run_dir / "metadata.json", local_metadata)
    attempt_state = ledger.snapshot()["attempts"].get(attempt["attempt_id"], {})
    peer_state = dict(attempt_state.get("cluster_submission") or peer)
    peer_state["collected_output"] = output
    peer_state["worker_node_id"] = remote_job.get("node_id")
    peer_state["resource_allocation"] = allocation
    ledger.transition(
        attempt["attempt_id"], "running", cluster_submission=peer_state,
        worker_node_id=remote_job.get("node_id"), resource_allocation=allocation,
        remote_completed_at=remote_metadata.get("completed_at"),
    )
    return True


def _attempt_child_inputs(parent_inputs: dict[str, Any], params: dict[str, Any]) -> dict[str, Any]:
    """Rebase staged target references into each isolated native child run."""
    target_type = str(params.get("target_type") or parent_inputs.get("target_type") or "structure")
    target_file = Path(str(params.get("target_file") or "target.pdb")).name
    if target_file not in {"target.pdb", "target.cif", "target.fasta"}:
        raise ValueError(f"Unsupported staged BindCraft2 target filename: {target_file!r}.")
    target_ref = f"artifacts/raw/bindcraft2/{target_file}"
    child_inputs = dict(parent_inputs)
    child_inputs.update(
        target_type=target_type,
        target_chains=[str(value) for value in params.get("target_chains") or parent_inputs.get("target_chains") or []],
        target_structure=target_ref if target_type == "structure" else None,
        target_pdb=target_ref if target_type == "structure" and target_file == "target.pdb" else None,
        target_fasta=target_ref if target_type == "sequence" else None,
        second_target_pdb=(
            "artifacts/raw/bindcraft2/second_target.pdb"
            if parent_inputs.get("second_target_pdb") else None
        ),
    )
    return child_inputs


def _bounded_child_settings(parent_raw: Path, child_params: dict[str, Any], attempt: dict[str, Any]) -> dict[str, Any]:
    """Read raw settings so paths relative to BindCraft2's working folder stay relative."""
    from mn_molecular_design_studio.engine.workflow_engine.workflows.bindcraft2 import bounded_bindcraft2_settings

    settings = json.loads((parent_raw / "settings.json").read_text(encoding="utf-8"))
    settings = bounded_bindcraft2_settings(settings, int(attempt["seed"]))
    settings["campaign_name"] = child_params["campaign_name"]
    settings["workers_per_gpu"] = 1
    settings["project_folder"] = "output"
    for target in settings.get("targets", []):
        target["target_path"] = Path(str(target.get("target_path") or "target.pdb")).name
    return settings


def _sequence_is_subsequence(candidate_sequence: str, target_sequence: str) -> bool:
    candidate_sequence = "".join(str(candidate_sequence or "").split()).upper()
    target_sequence = "".join(str(target_sequence or "").split()).upper()
    if not candidate_sequence or not target_sequence:
        return False
    position = 0
    for residue in target_sequence:
        if position < len(candidate_sequence) and candidate_sequence[position] == residue:
            position += 1
    return position == len(candidate_sequence)


def _native_candidate_rows(
    child_run_dir: Path,
    campaign_id: str,
    attempt: dict[str, Any],
    native_state: dict[str, Any],
) -> list[dict[str, Any]]:
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import (
        _resolve_candidate_structure_path,
        read_candidates,
    )
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json
    from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import portable_path

    candidates_path = child_run_dir / "artifacts" / "normalized_candidates" / "candidates.jsonl"
    candidates = read_candidates(candidates_path)
    passing = [
        row for row in candidates
        if row.get("metrics", {}).get("native_pass_filters") is True
    ]
    if len(passing) != int(native_state.get("accepted", -1)):
        raise RuntimeError(
            "BindCraft2 native acceptance does not reconcile with normalized candidates "
            f"for {child_run_dir.name}: state={native_state.get('accepted')} candidates={len(passing)}."
        )
    metadata = read_json(child_run_dir / "metadata.json")
    allocation = metadata.get("resource_allocation") if isinstance(metadata.get("resource_allocation"), dict) else {}
    params = read_json(child_run_dir / "input.json").get("params") or {}
    expected_target_chains = [str(value) for value in params.get("target_chains") or []]
    result = []
    seen_native_ids: set[str] = set()
    for row in candidates:
        native_id = str(row.get("candidate_id") or "")
        if not native_id or native_id in seen_native_ids:
            raise RuntimeError(f"BindCraft2 returned a missing or duplicate candidate ID in {child_run_dir.name}.")
        seen_native_ids.add(native_id)
        if row.get("metrics", {}).get("native_pass_filters") is True:
            if not row.get("binder_sequence") or not row.get("binder_chains"):
                raise RuntimeError(f"Native-accepted candidate {native_id} is missing its sequence or binder chain.")
            structure = _resolve_candidate_structure_path(child_run_dir, row)
            if structure is None or not structure.is_file():
                raise RuntimeError(f"Native-accepted candidate {native_id} has no verified complex structure.")
            if expected_target_chains and sorted(row.get("target_chains") or []) != sorted(expected_target_chains):
                raise RuntimeError(f"Native-accepted candidate {native_id} has mismatched target-chain identity.")
            from mn_molecular_design_studio.design_validation import protein_chains, read_structure, sequence

            structure_chains = protein_chains(read_structure(structure))
            binder_ids = [str(value) for value in row.get("binder_chains") or []]
            target_ids = [str(value) for value in row.get("target_chains") or []]
            if any(value not in structure_chains for value in binder_ids + target_ids):
                raise RuntimeError(f"Native-accepted candidate {native_id} references a missing protein chain.")
            actual_binder = ":".join(sequence(structure_chains[value]) for value in binder_ids)
            expected_binder = str(row.get("binder_sequence") or "").replace(" ", "").upper()
            if actual_binder.replace(":", "") != expected_binder.replace(":", ""):
                raise RuntimeError(f"Native-accepted candidate {native_id} binder sequence does not match its complex structure.")
            raw_target_root = child_run_dir / "artifacts" / "raw" / "bindcraft2"
            target_input = next((raw_target_root / name for name in ("target.pdb", "target.cif") if (raw_target_root / name).is_file()), None)
            requested_target_ids = [str(value) for value in params.get("target_chains") or []]
            if target_input is not None and requested_target_ids:
                expected_structure_chains = protein_chains(read_structure(target_input))
                if any(value not in expected_structure_chains for value in requested_target_ids):
                    raise RuntimeError("The staged BindCraft2 target structure no longer contains the requested chain identity.")
                expected_targets = [sequence(expected_structure_chains[value]) for value in requested_target_ids]
                actual_targets = [sequence(structure_chains[value]) for value in target_ids]
                if len(expected_targets) != len(actual_targets) or any(
                    not _sequence_is_subsequence(actual, expected)
                    for actual, expected in zip(actual_targets, expected_targets)
                ):
                    raise RuntimeError(f"Native-accepted candidate {native_id} target sequence is not an exact crop of the configured target.")
        import copy

        candidate = copy.deepcopy(row)
        candidate["candidate_id"] = f"dc_{int(attempt['ordinal']):08d}_{hashlib.sha256(f'{campaign_id}:{native_id}'.encode()).hexdigest()[:16]}"
        raw = candidate.setdefault("raw_metadata", {})
        raw["source_run_dir"] = portable_path(child_run_dir)
        raw["native_candidate_id"] = native_id
        raw["distributed_campaign"] = {
            "campaign_id": campaign_id,
            "attempt_id": attempt["attempt_id"],
            "batch_id": attempt["batch_id"],
            "attempt_seed": attempt["seed"],
            "trajectory_index": 1,
            "worker_run_id": child_run_dir.name,
            "worker": metadata.get("scheduler_worker_id") or metadata.get("worker_id") or "",
            "gpu_device": allocation.get("gpu_device") or params.get("gpu_device") or "",
            "gpu_ids": allocation.get("gpu_ids") or [],
            "native_state": native_state,
            "native_acceptance": "accepted" if row.get("metrics", {}).get("native_pass_filters") is True else "rejected",
        }
        result.append(candidate)
    return result


def _rfd4_candidate_rows(
    child_run_dir: Path,
    campaign_id: str,
    attempt: dict[str, Any],
    native_state: dict[str, Any],
) -> list[dict[str, Any]]:
    """Verify one RFD4 co-designed sequence/structure output without inventing native acceptance."""
    import copy

    from mn_molecular_design_studio.design_validation import protein_chains, read_structure, sequence
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import (
        _resolve_candidate_structure_path,
        read_candidates,
    )
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json
    from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import portable_path, resolve_stored_path
    from mn_molecular_design_studio.rfd4_proteina_input import _cif_label_mapping

    child_run_dir = Path(child_run_dir)
    candidates_path = child_run_dir / "artifacts" / "normalized_candidates" / "candidates.jsonl"
    candidates = read_candidates(candidates_path)
    if len(candidates) > 1:
        raise RuntimeError(
            f"RFD4 attempt {attempt['attempt_id']} requested one output but normalized {len(candidates)} candidates."
        )
    payload = read_json(child_run_dir / "input.json")
    params = payload.get("params") if isinstance(payload.get("params"), dict) else {}
    target_path = resolve_stored_path(params.get("target_structure"), must_exist=True)
    if target_path is None or not target_path.is_file():
        raise RuntimeError(f"RFD4 attempt {attempt['attempt_id']} has no readable staged target structure.")
    target_chains = protein_chains(read_structure(target_path))
    selected_source_ids = [str(value) for value in params.get("target_chains") or []]
    if not selected_source_ids or any(value not in target_chains for value in selected_source_ids):
        raise RuntimeError(f"RFD4 attempt {attempt['attempt_id']} has invalid source target-chain identities.")
    source_to_label, _residue_mapping = _cif_label_mapping(target_path)
    selected_target_sequences: dict[str, tuple[str, str]] = {}
    for source_id in selected_source_ids:
        label_id = source_to_label.get(source_id)
        if not label_id or label_id in selected_target_sequences:
            raise RuntimeError(f"RFD4 attempt {attempt['attempt_id']} has an ambiguous source-to-label target-chain mapping.")
        selected_target_sequences[label_id] = (source_id, sequence(target_chains[source_id]))
    metadata = read_json(child_run_dir / "metadata.json")
    allocation = metadata.get("resource_allocation") if isinstance(metadata.get("resource_allocation"), dict) else {}
    result: list[dict[str, Any]] = []
    seen_native_ids: set[str] = set()
    for row in candidates:
        native_id = str(row.get("candidate_id") or "")
        if not native_id or native_id in seen_native_ids:
            raise RuntimeError(f"RFD4 returned a missing or duplicate candidate ID in {child_run_dir.name}.")
        seen_native_ids.add(native_id)
        binder_ids = [str(value) for value in row.get("binder_chains") or []]
        target_ids = [str(value) for value in row.get("target_chains") or []]
        if not binder_ids or not target_ids:
            raise RuntimeError(f"RFD4 candidate {native_id} is missing explicit target/binder roles.")
        structure = _resolve_candidate_structure_path(child_run_dir, row)
        if structure is None or not structure.is_file():
            raise RuntimeError(f"RFD4 candidate {native_id} has no readable co-designed complex structure.")
        structure_chains = protein_chains(read_structure(structure))
        if any(chain_id not in structure_chains for chain_id in binder_ids + target_ids):
            raise RuntimeError(f"RFD4 candidate {native_id} references a missing protein chain.")
        actual_binder = ":".join(sequence(structure_chains[chain_id]) for chain_id in binder_ids)
        expected_binder = str(row.get("binder_sequence") or "").replace(" ", "").upper()
        if not expected_binder or actual_binder.replace(":", "") != expected_binder.replace(":", ""):
            raise RuntimeError(f"RFD4 candidate {native_id} sequence differs from its co-designed complex structure.")
        if len(target_ids) != len(selected_source_ids):
            raise RuntimeError(f"RFD4 candidate {native_id} changed the number of selected target chains.")
        target_mapping = []
        unmatched = dict(selected_target_sequences)
        for prediction_chain in target_ids:
            predicted_sequence = sequence(structure_chains[prediction_chain])
            matches = [
                label_chain for label_chain, (_source_chain, source_sequence) in unmatched.items()
                if label_chain == prediction_chain and _sequence_is_subsequence(predicted_sequence, source_sequence)
            ]
            if len(matches) != 1:
                raise RuntimeError(
                    f"RFD4 candidate {native_id} target chain {prediction_chain} does not map uniquely to a selected source chain."
                )
            label_chain = matches[0]
            source_chain, source_sequence = unmatched.pop(label_chain)
            target_mapping.append({
                "source_chain": source_chain,
                "prediction_chain": prediction_chain,
                "source_sequence": source_sequence,
                "prediction_sequence": predicted_sequence,
            })
        minimum, maximum = int(params.get("length_min", 60)), int(params.get("length_max", 120))
        if not minimum <= len(actual_binder.replace(":", "")) <= maximum:
            raise RuntimeError(f"RFD4 candidate {native_id} is outside its requested binder-length range.")

        candidate = copy.deepcopy(row)
        candidate["candidate_id"] = (
            f"dc_{int(attempt['ordinal']):08d}_"
            f"{hashlib.sha256(f'{campaign_id}:{native_id}'.encode()).hexdigest()[:16]}"
        )
        raw = candidate.setdefault("raw_metadata", {})
        raw["source_run_dir"] = portable_path(child_run_dir)
        raw["native_candidate_id"] = native_id
        raw["source_to_prediction_components"] = {
            "target": target_mapping,
            "designed_binder": [{
                "source_chain": None,
                "prediction_chain": chain_id,
                "prediction_sequence": sequence(structure_chains[chain_id]),
            } for chain_id in binder_ids],
        }
        raw["distributed_campaign"] = {
            "campaign_id": campaign_id,
            "generator_id": "rfd4_proteina",
            "attempt_id": attempt["attempt_id"],
            "batch_id": attempt["batch_id"],
            "attempt_seed": attempt["seed"],
            "attempt_unit": "one native design request",
            "worker_run_id": child_run_dir.name,
            "worker": metadata.get("scheduler_worker_id") or metadata.get("worker_id") or "",
            "gpu_device": allocation.get("gpu_device") or params.get("gpu_device") or "",
            "gpu_ids": allocation.get("gpu_ids") or [],
            "native_state": native_state,
            "native_acceptance": "unavailable",
            "sequence_structure_mode": "co_design",
        }
        result.append(candidate)
    return result


def _create_rfd4_attempt_child(
    parent_run_dir: Path,
    ledger: CampaignLedger,
    attempt: dict[str, Any],
    config: dict[str, Any],
) -> Path:
    """Queue one independently seeded, one-output RFD4 co-design request."""
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import (
        mark_internal_job,
        read_json,
        submission_tags,
        write_json,
    )
    from mn_molecular_design_studio.engine.workflow_engine.core.local_worker import spawn_worker_for_run
    from mn_molecular_design_studio.engine.workflow_engine.core.workflow_queue import (
        enqueue_workflow_call,
        execution_resources,
    )
    from mn_molecular_design_studio.engine.workflow_engine.workflows import rfd4_proteina

    parent_run_dir = Path(parent_run_dir)
    parent_payload = read_json(parent_run_dir / "input.json")
    parent_params = parent_payload.get("params") if isinstance(parent_payload.get("params"), dict) else {}
    if parent_params.get("task") != "protein_binder" or parent_params.get("input_mode") not in {"guided", "guided_binder"}:
        raise ValueError("RFD4 campaign attempts support the guided protein-binder workflow only.")
    if parent_params.get("input_files") or parent_params.get("avoid_structures") or parent_params.get("also_bind_structures"):
        raise ValueError("RFD4 campaign attempts currently require one target structure without secondary specificity files.")
    request_path = parent_run_dir / "artifacts" / "rfd4_proteina" / "request.json"
    request = read_json(request_path)
    staged_target = rfd4_proteina.local_path(parent_run_dir, request.get("target_structure"))
    if staged_target is None or not staged_target.is_file():
        raise ValueError("RFD4 campaign parent has no staged target structure to share with attempts.")

    peer_execution = _peer_campaign_mode(config)
    child_kwargs = dict(parent_params)
    child_kwargs.update(
        target_structure=str(staged_target),
        input_files=[],
        avoid_structures=[],
        also_bind_structures=[],
        num_designs=1,
        batch_size=1,
        seed=int(attempt["seed"]),
        gpu_device="auto" if peer_execution else _gpu_for_attempt(config, int(attempt["ordinal"])),
        cpu_cores=int((config.get("resources") or {}).get("cpu_cores") or parent_params.get("cpu_cores") or 8),
        job_name=f"{parent_params.get('job_name') or parent_run_dir.name}-{attempt['ordinal']:08d}",
        distributed_campaign=None,
        task_group="design",
    )
    resource_options = config.get("resources") if isinstance(config.get("resources"), dict) else {}
    resource_request = {
        "gpu": True,
        "gpu_device": child_kwargs["gpu_device"],
        "cpu_cores": child_kwargs["cpu_cores"],
        "min_vram_gb": float(resource_options.get("min_vram_gb") or 16.0),
        "ram_gb": float(resource_options.get("ram_gb") or 0),
        "scratch_gb": float(resource_options.get("scratch_gb") or 0),
    }
    campaign_id = ledger.snapshot()["campaign_id"]
    with submission_tags(
        distributed_campaign_id=campaign_id,
        distributed_attempt_id=attempt["attempt_id"],
        distributed_batch_id=attempt["batch_id"],
        parent_campaign_run_dir=str(parent_run_dir),
    ), execution_resources(resource_request):
        child_run_dir = enqueue_workflow_call(
            rfd4_proteina.run_rfd4_proteina_design,
            task_group="design",
            tool="rfd4_proteina",
            job_type="rfd4_proteina_design",
            kwargs=child_kwargs,
            cpu_cores=int(child_kwargs["cpu_cores"]),
            gpu_device=child_kwargs["gpu_device"],
            dispatch=False,
        )
    mark_internal_job(
        child_run_dir, parent_run_dir=parent_run_dir,
        parent_task_group=parent_run_dir.parent.name,
        role="distributed_rfd4_attempt", engine="rfd4_proteina",
    )
    metadata = read_json(child_run_dir / "metadata.json")
    metadata.update(
        distributed_campaign_id=campaign_id,
        distributed_attempt_id=attempt["attempt_id"],
        distributed_batch_id=attempt["batch_id"],
        benchmark_engine="rfd4_proteina",
        resource_request=resource_request,
        cluster_remote=peer_execution,
    )
    write_json(child_run_dir / "metadata.json", metadata)
    ledger.transition(
        attempt["attempt_id"], "queued", child_run_dir=str(child_run_dir),
        native_campaign_seed=int(attempt["seed"]), attempt_unit="one native design request",
        requested_outputs=1, cluster_remote=peer_execution,
    )
    if peer_execution:
        current_attempt = ledger.snapshot()["attempts"][attempt["attempt_id"]]
        _submit_rfd4_peer_attempt(parent_run_dir, ledger, current_attempt, child_run_dir)
    else:
        spawn_worker_for_run(child_run_dir)
    return child_run_dir


def _repair_rfd4_attempt_child(
    parent_run_dir: Path,
    child_run_dir: Path,
    attempt: dict[str, Any],
    config: dict[str, Any],
) -> None:
    """Repair an interrupted RFD4 queue write before dispatching the child."""
    import inspect

    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import (
        mark_internal_job,
        read_json,
        write_json,
    )
    from mn_molecular_design_studio.engine.workflow_engine.core.local_worker import spawn_worker_for_run
    from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import portable_path
    from mn_molecular_design_studio.engine.workflow_engine.runtime import runs_root
    from mn_molecular_design_studio.engine.workflow_engine.workflows.rfd4_proteina import run_rfd4_proteina_design

    parent_run_dir, child_run_dir = Path(parent_run_dir), Path(child_run_dir)
    child_payload = read_json(child_run_dir / "input.json")
    params = child_payload.get("params") if isinstance(child_payload.get("params"), dict) else {}
    peer_execution = _peer_campaign_mode(config)
    params.update(
        num_designs=1, batch_size=1, seed=int(attempt["seed"]),
        gpu_device="auto" if peer_execution else _gpu_for_attempt(config, int(attempt["ordinal"])),
        distributed_campaign=None,
    )
    child_payload["params"] = params
    write_json(child_run_dir / "input.json", child_payload)
    metadata = read_json(child_run_dir / "metadata.json")
    resource_request = metadata.get("resource_request") if isinstance(metadata.get("resource_request"), dict) else {}
    campaign_id = CampaignLedger(_campaign_ledger_path(parent_run_dir)).snapshot()["campaign_id"]
    mark_internal_job(
        child_run_dir, parent_run_dir=parent_run_dir,
        parent_task_group=parent_run_dir.parent.name,
        role="distributed_rfd4_attempt", engine="rfd4_proteina",
    )
    metadata = read_json(child_run_dir / "metadata.json")
    metadata.update(
        distributed_campaign_id=campaign_id,
        distributed_attempt_id=attempt["attempt_id"],
        distributed_batch_id=attempt["batch_id"],
        benchmark_engine="rfd4_proteina",
        resource_request=resource_request,
        cluster_remote=peer_execution,
    )
    write_json(child_run_dir / "metadata.json", metadata)
    request_path = child_run_dir / "worker_request.json"
    if not request_path.is_file():
        kwargs = {
            key: value for key, value in params.items()
            if key in inspect.signature(run_rfd4_proteina_design).parameters
        }
        write_json(child_run_dir / "worker_request.json", {
            "kind": "workflow_call",
            "module": "mn_molecular_design_studio.engine.workflow_engine.workflows.rfd4_proteina",
            "function": "run_rfd4_proteina_design",
            "args": [], "kwargs": kwargs,
            "resources": {**resource_request, "gpu_device": params.get("gpu_device"), "cpu_cores": params.get("cpu_cores")},
        })
        write_json(child_run_dir / "command.json", {
            "mode": "local_worker",
            "command": [
                "python", "-m", "mn_molecular_design_studio.engine.workflow_engine.core.local_worker",
                "--run-dir", portable_path(child_run_dir, run_dir=runs_root()),
            ],
        })
    if metadata.get("status") == "queued" and not metadata.get("scheduler_dispatch_requested"):
        if peer_execution:
            _submit_rfd4_peer_attempt(parent_run_dir, CampaignLedger(_campaign_ledger_path(parent_run_dir)), attempt, child_run_dir)
        else:
            spawn_worker_for_run(child_run_dir)


def _create_attempt_child(
    parent_run_dir: Path,
    ledger: CampaignLedger,
    attempt: dict[str, Any],
    config: dict[str, Any],
) -> Path:
    if str(config.get("generator_id") or "bindcraft2") == "rfd4_proteina":
        return _create_rfd4_attempt_child(parent_run_dir, ledger, attempt, config)
    import copy
    import shutil

    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import (
        JobPaths,
        mark_internal_job,
        read_json,
        submission_tags,
        write_json,
    )
    from mn_molecular_design_studio.engine.workflow_engine.core.local_worker import spawn_worker_for_run
    from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import portable_path
    from mn_molecular_design_studio.engine.workflow_engine.runtime import runs_root
    from mn_molecular_design_studio.engine.workflow_engine.workflows.bindcraft2 import (
        BINDCRAFT2_TOOL,
        DESIGN_GROUP,
        bounded_bindcraft2_settings,
    )
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import create_job

    parent_payload = read_json(parent_run_dir / "input.json")
    parent_params = parent_payload.get("params") if isinstance(parent_payload.get("params"), dict) else {}
    child_params = copy.deepcopy(parent_params)
    child_params.pop("distributed_campaign", None)
    child_params["target_file"] = Path(str(child_params.get("target_file") or "target.pdb")).name
    child_params["campaign_name"] = f"{parent_params.get('campaign_name') or parent_run_dir.name}-{attempt['ordinal']:08d}"
    child_params["number_of_final_designs"] = 1
    child_params["max_trajectories"] = 1
    child_params["campaign_seed"] = int(attempt["seed"])
    child_params["gpu_device"] = _gpu_for_attempt(config, int(attempt["ordinal"]))
    child_params["workers_per_gpu"] = 1
    resource_config = config.get("resources") if isinstance(config.get("resources"), dict) else {}
    resource_request = {
        "gpu": True,
        "gpu_device": child_params["gpu_device"],
        "cpu_cores": int(resource_config.get("cpu_cores") or parent_params.get("cpu_cores") or 4),
        "min_vram_gb": float(resource_config.get("min_vram_gb") or 0),
        "ram_gb": float(resource_config.get("ram_gb") or 0),
        "scratch_gb": float(resource_config.get("scratch_gb") or 0),
    }
    child_params["cpu_cores"] = resource_request["cpu_cores"]
    child_params["resource_request"] = resource_request
    campaign_id = ledger.snapshot()["campaign_id"]
    with submission_tags(
        distributed_campaign_id=campaign_id,
        distributed_attempt_id=attempt["attempt_id"],
        distributed_batch_id=attempt["batch_id"],
        parent_campaign_run_dir=str(parent_run_dir),
    ):
        child_job = create_job(
            DESIGN_GROUP,
            "bindcraft2_design",
            BINDCRAFT2_TOOL,
            inputs=_attempt_child_inputs(parent_payload.get("inputs") or {}, child_params),
            params=child_params,
        )
    child_run_dir = child_job.run_dir
    mark_internal_job(
        child_run_dir,
        parent_run_dir=parent_run_dir,
        parent_task_group=parent_run_dir.parent.name,
        role="distributed_bindcraft2_attempt",
        engine="bindcraft2",
    )
    parent_raw = parent_run_dir / "artifacts" / "raw" / "bindcraft2"
    child_raw = child_run_dir / "artifacts" / "raw" / "bindcraft2"
    child_raw.mkdir(parents=True, exist_ok=True)
    for filename in ("target.pdb", "target.cif", "target.fasta", "second_target.pdb"):
        source = parent_raw / filename
        if source.is_file():
            shutil.copyfile(source, child_raw / filename)
    settings = _bounded_child_settings(parent_raw, child_params, attempt)
    write_json(child_raw / "settings.json", settings)
    metadata = read_json(child_run_dir / "metadata.json")
    metadata.update(
        distributed_campaign_id=campaign_id,
        distributed_attempt_id=attempt["attempt_id"],
        distributed_batch_id=attempt["batch_id"],
        resource_request=resource_request,
        queue_resource=("gpu" if child_params["gpu_device"] in {"auto", "all", "eligible"} else f"gpu:{child_params['gpu_device']}"),
    )
    write_json(child_run_dir / "metadata.json", metadata)
    write_json(
        child_run_dir / "worker_request.json",
        {"kind": "bindcraft2_design", "kwargs": {}, "path_kwargs": [], "resources": resource_request},
    )
    write_json(
        child_run_dir / "command.json",
        {
            "mode": "local_worker",
            "command": [
                "python", "-m", "mn_molecular_design_studio.engine.workflow_engine.core.local_worker",
                "--run-dir", portable_path(child_run_dir, run_dir=runs_root()),
            ],
        },
    )
    ledger.transition(
        attempt["attempt_id"], "queued", child_run_dir=str(child_run_dir),
        native_campaign_seed=int(attempt["seed"]), trajectory_index=1,
    )
    spawn_worker_for_run(child_run_dir)
    return child_run_dir


def _repair_attempt_child(
    parent_run_dir: Path,
    child_run_dir: Path,
    attempt: dict[str, Any],
    config: dict[str, Any],
) -> None:
    """Finish durable staging if the controller stopped between job creation and dispatch."""
    if str(config.get("generator_id") or "bindcraft2") == "rfd4_proteina":
        return _repair_rfd4_attempt_child(parent_run_dir, child_run_dir, attempt, config)
    import shutil

    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import (
        mark_internal_job,
        read_json,
        write_json,
    )
    from mn_molecular_design_studio.engine.workflow_engine.core.local_worker import spawn_worker_for_run
    from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import portable_path
    from mn_molecular_design_studio.engine.workflow_engine.runtime import runs_root
    from mn_molecular_design_studio.engine.workflow_engine.workflows.bindcraft2 import bounded_bindcraft2_settings

    child_payload = read_json(child_run_dir / "input.json")
    params = child_payload.get("params") if isinstance(child_payload.get("params"), dict) else {}
    parent_payload = read_json(parent_run_dir / "input.json")
    parent_params = parent_payload.get("params") if isinstance(parent_payload.get("params"), dict) else {}
    params["target_file"] = Path(str(params.get("target_file") or parent_params.get("target_file") or "target.pdb")).name
    child_payload["params"] = params
    child_payload["inputs"] = _attempt_child_inputs(parent_payload.get("inputs") or {}, params)
    write_json(child_run_dir / "input.json", child_payload)
    campaign_id = CampaignLedger(_campaign_ledger_path(parent_run_dir)).snapshot()["campaign_id"]
    child_raw = child_run_dir / "artifacts" / "raw" / "bindcraft2"
    parent_raw = parent_run_dir / "artifacts" / "raw" / "bindcraft2"
    child_raw.mkdir(parents=True, exist_ok=True)
    for filename in ("target.pdb", "target.cif", "target.fasta", "second_target.pdb"):
        source = parent_raw / filename
        destination = child_raw / filename
        if source.is_file() and not destination.is_file():
            shutil.copyfile(source, destination)
    settings = _bounded_child_settings(
        parent_raw,
        {**params, "campaign_name": params.get("campaign_name") or child_run_dir.name},
        attempt,
    )
    write_json(child_raw / "settings.json", settings)
    mark_internal_job(
        child_run_dir, parent_run_dir=parent_run_dir,
        parent_task_group=parent_run_dir.parent.name,
        role="distributed_bindcraft2_attempt", engine="bindcraft2",
    )
    metadata = read_json(child_run_dir / "metadata.json")
    metadata.update(
        distributed_campaign_id=campaign_id,
        distributed_attempt_id=attempt["attempt_id"],
        distributed_batch_id=attempt["batch_id"],
        resource_request=params.get("resource_request") or {},
    )
    write_json(child_run_dir / "metadata.json", metadata)
    if not (child_run_dir / "worker_request.json").is_file():
        write_json(
            child_run_dir / "worker_request.json",
            {
                "kind": "bindcraft2_design", "kwargs": {}, "path_kwargs": [],
                "resources": params.get("resource_request") or {},
            },
        )
        write_json(
            child_run_dir / "command.json",
            {
                "mode": "local_worker",
                "command": [
                    "python", "-m", "mn_molecular_design_studio.engine.workflow_engine.core.local_worker",
                    "--run-dir", portable_path(child_run_dir, run_dir=runs_root()),
                ],
            },
        )
    if metadata.get("status") == "queued" and not metadata.get("scheduler_dispatch_requested"):
        spawn_worker_for_run(child_run_dir)


def _find_unlinked_children(ledger: CampaignLedger) -> dict[str, Path]:
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json
    from mn_molecular_design_studio.engine.workflow_engine.runtime import runs_root

    state = ledger.snapshot()
    unlinked = {
        attempt_id for attempt_id, item in state["attempts"].items()
        if not item.get("child_run_dir") and item.get("status") not in {"completed", "failed", "interrupted", "cancelled", "not_started"}
    }
    if not unlinked:
        return {}
    found: dict[str, Path] = {}
    design_root = runs_root() / "design"
    if not design_root.is_dir():
        return found
    for child in design_root.iterdir():
        if not child.is_dir():
            continue
        metadata = read_json(child / "metadata.json")
        if metadata.get("distributed_campaign_id") != state["campaign_id"]:
            continue
        attempt_id = str(metadata.get("distributed_attempt_id") or "")
        if attempt_id in unlinked:
            if attempt_id in found:
                raise RuntimeError(f"Multiple generator child runs claim attempt {attempt_id}.")
            found[attempt_id] = child
    return found


def _attempt_child_paths(ledger: CampaignLedger) -> dict[str, Path]:
    from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path

    state = ledger.snapshot()
    found = _find_unlinked_children(ledger)
    children: dict[str, Path] = {}
    for attempt_id, attempt in state["attempts"].items():
        child_text = str(attempt.get("child_run_dir") or "")
        child = resolve_stored_path(child_text, must_exist=True) if child_text else found.get(attempt_id)
        if child is not None and child.is_dir():
            children[attempt_id] = child
    return children


def _accepted_count(ledger: CampaignLedger, config: dict[str, Any]) -> int:
    counters = ledger.counters()
    key = {
        "native": "native_accepted_candidates",
        "validated": "validated_accepted_candidates",
        "diversity_qualified": "diversity_qualified_accepted_candidates",
    }[config["acceptance_level"]]
    return counters[key]


def _rfd4_retry_safe(child_run_dir: Path) -> bool:
    """Only retry an RFD4 child if its native Docker process was never launched."""
    command_path = Path(child_run_dir) / "command.json"
    if not command_path.is_file():
        return True
    try:
        command = json.loads(command_path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return False
    return command.get("mode") != "docker"


def _reconcile_attempts(
    ledger: CampaignLedger,
    parent_run_dir: Path,
    config: dict[str, Any],
) -> tuple[dict[str, Path], list[dict[str, Any]], list[dict[str, Any]]]:
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json, resume_job
    generator_id = str(config.get("generator_id") or "bindcraft2")
    if generator_id == "bindcraft2":
        from mn_molecular_design_studio.engine.workflow_engine.workflows.bindcraft2 import read_bindcraft2_campaign_state
    else:
        read_bindcraft2_campaign_state = None

    state = ledger.snapshot()
    children = _attempt_child_paths(ledger)
    for attempt_id, attempt in list(ledger.snapshot()["attempts"].items()):
        if attempt.get("status") != "authorized":
            continue
        child = children.get(attempt_id)
        if child is None:
            try:
                child = _create_attempt_child(parent_run_dir, ledger, attempt, config)
            except Exception as exc:
                ledger.transition(attempt_id, "failed", failure=f"{type(exc).__name__}: {exc}")
                continue
            children[attempt_id] = child
        elif not (child / "worker_request.json").is_file():
            _repair_attempt_child(parent_run_dir, child, attempt, config)
            ledger.transition(attempt_id, "queued", child_run_dir=str(child))
        elif str(read_json(child / "metadata.json").get("status") or "") == "queued":
            dispatch_meta = read_json(child / "metadata.json")
            if not dispatch_meta.get("scheduler_dispatch_requested") and not dispatch_meta.get("scheduler_job_id"):
                if dispatch_meta.get("cluster_remote"):
                    _submit_rfd4_peer_attempt(parent_run_dir, ledger, attempt, child)
                else:
                    from mn_molecular_design_studio.engine.workflow_engine.core.local_worker import spawn_worker_for_run

                    spawn_worker_for_run(child)
                ledger.transition(attempt_id, "queued", child_run_dir=str(child))
    active: list[dict[str, Any]] = []
    candidate_rows: list[dict[str, Any]] = []
    for attempt_id, attempt in state["attempts"].items():
        child = children.get(attempt_id)
        if child is None:
            if attempt["status"] == "authorized":
                continue
            if attempt["status"] in {"queued", "running"}:
                retry_allowed = (
                    generator_id != "rfd4_proteina"
                    and int(attempt.get("retry_count", 0)) < int(config["max_retries"])
                )
                if retry_allowed:
                    if attempt["status"] != "interrupted":
                        ledger.transition(attempt_id, "interrupted", failure="Origin child run is unavailable for reconciliation.")
                    retried = ledger.retry(attempt_id)
                    try:
                        child = _create_attempt_child(parent_run_dir, ledger, retried, config)
                        children[attempt_id] = child
                    except Exception as exc:
                        ledger.transition(attempt_id, "failed", failure=f"{type(exc).__name__}: {exc}")
                else:
                    ledger.transition(attempt_id, "interrupted", failure="Origin child run is unavailable for reconciliation.")
            latest_attempt = ledger.snapshot()["attempts"].get(attempt_id, {})
            if latest_attempt.get("status") in {"authorized", "queued", "running"}:
                active.append(latest_attempt)
            continue
        attempt = ledger.snapshot()["attempts"].get(attempt_id, attempt)
        child_metadata = read_json(child / "metadata.json")
        if (
            generator_id == "rfd4_proteina"
            and attempt.get("status") not in {"completed", "failed", "interrupted", "cancelled", "not_started"}
            and (attempt.get("cluster_submission") or child_metadata.get("cluster_remote"))
        ):
            from mn_molecular_design_studio.cluster import candidate_batches

            peer = dict(attempt.get("cluster_submission") or {})
            try:
                if not peer.get("cluster_job_id"):
                    peer = _submit_rfd4_peer_attempt(parent_run_dir, ledger, attempt, child)
                    attempt = ledger.snapshot()["attempts"][attempt_id]
                    peer = dict(attempt.get("cluster_submission") or peer)
                owner_id, project_id, development_mode = _campaign_actor(parent_run_dir)
                peer_status, peer_output, peer_job = candidate_batches.poll(
                    peer, owner_id=owner_id, project_id=project_id,
                    allow_development=development_mode,
                )
            except Exception as exc:
                failure = f"{type(exc).__name__}: {exc}"
                child_metadata.update(status="failed", worker_error=failure, completed_at=utc_now())
                from mn_molecular_design_studio.engine.workflow_engine.core.jobs import write_json
                write_json(child / "metadata.json", child_metadata)
                ledger.transition(attempt_id, "failed", child_run_dir=str(child), failure=failure)
                continue
            peer["cluster_job_id"] = peer.get("cluster_job_id") or peer_job.get("cluster_job_id")
            peer["worker_node_id"] = peer_job.get("node_id")
            peer["wait_reason"] = peer_job.get("wait_reason") or peer.get("wait_reason") or ""
            remote_metadata = peer_job.get("metadata") if isinstance(peer_job.get("metadata"), dict) else {}
            ledger.transition(
                attempt_id,
                "running" if peer_status == "running" else "queued",
                cluster_submission=peer,
                cluster_job_id=peer.get("cluster_job_id"),
                worker_node_id=peer.get("worker_node_id"),
                wait_reason=peer.get("wait_reason"),
            )
            child_metadata.update(
                cluster_remote=True,
                cluster_job_id=peer.get("cluster_job_id"),
                execution_node_id=peer.get("worker_node_id"),
                remote_status=peer_status,
                remote_wait_reason=peer.get("wait_reason"),
                resource_allocation=remote_metadata.get("resource_allocation") or child_metadata.get("resource_allocation") or {},
                worker_started_at=remote_metadata.get("worker_started_at") or child_metadata.get("worker_started_at"),
                scheduler_worker_id=remote_metadata.get("scheduler_worker_id") or peer.get("worker_node_id"),
            )
            if peer_status == "completed" and peer_output is not None:
                if not _materialize_rfd4_peer_output(
                    parent_run_dir, ledger, attempt, child, peer, peer_output, peer_job,
                ):
                    latest = ledger.snapshot()["attempts"][attempt_id]
                    active.append({**latest, "worker_status": "collecting", "worker_node_id": peer.get("worker_node_id")})
                    continue
                attempt = ledger.snapshot()["attempts"][attempt_id]
            elif peer_status in {"failed", "cancelled", "lost"}:
                terminal = "cancelled" if peer_status == "cancelled" else "failed"
                child_metadata.update(status=terminal, worker_error=peer_job.get("wait_reason") or peer_status, completed_at=utc_now())
                from mn_molecular_design_studio.engine.workflow_engine.core.jobs import write_json
                write_json(child / "metadata.json", child_metadata)
                ledger.transition(
                    attempt_id, terminal, child_run_dir=str(child),
                    failure=peer_job.get("wait_reason") or peer_status,
                )
                continue
            else:
                child_metadata["status"] = "running" if peer_status == "running" else "queued"
                from mn_molecular_design_studio.engine.workflow_engine.core.jobs import write_json
                write_json(child / "metadata.json", child_metadata)
                latest = ledger.snapshot()["attempts"][attempt_id]
                active.append({
                    **latest,
                    "worker_status": peer_status,
                    "worker_node_id": peer.get("worker_node_id"),
                    "wait_reason": peer.get("wait_reason"),
                })
                continue
        metadata = read_json(child / "metadata.json")
        status = str(metadata.get("status") or "unknown")
        result = read_json(child / "result.json")
        if status in {"completed", "failed", "cancelled"}:
            ledger.record_compute(attempt_id, child.name, _worker_gpu_seconds(metadata))
        if status in {"queued", "preparing", "running"}:
            if attempt["status"] != status and attempt["status"] in {"authorized", "queued", "running"}:
                ledger.transition(attempt_id, "running" if status == "running" else "queued", child_run_dir=str(child))
            if read_bindcraft2_campaign_state is not None:
                try:
                    progress = read_bindcraft2_campaign_state(child)
                except RuntimeError:
                    progress = None  # Native writes this small progress file while the worker is running.
            else:
                progress = {
                    "authoritative": False,
                    "phase": metadata.get("current_phase") or metadata.get("progress_label") or status,
                    "requested_outputs": 1,
                }
            active.append({**attempt, "worker_status": status, "native_progress": progress or {}})
            continue
        if attempt["status"] == "completed":
            continue
        if status == "completed" and result.get("success") is True:
            if generator_id == "bindcraft2":
                native_state = read_bindcraft2_campaign_state(child)
                if native_state is None:
                    raise RuntimeError(f"Completed BindCraft2 attempt {attempt_id} has no native .campaign_state.json.")
                if native_state["trajectories"] != 1:
                    raise RuntimeError(
                        f"Bounded BindCraft2 attempt {attempt_id} reports {native_state['trajectories']} native trajectories; expected 1."
                    )
                rows = _native_candidate_rows(child, state["campaign_id"], attempt, native_state)
            elif generator_id == "rfd4_proteina":
                from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates

                candidate_path = child / "artifacts" / "normalized_candidates" / "candidates.jsonl"
                output_count = len(read_candidates(candidate_path))
                native_state = {
                    "engine": "rfd4_proteina",
                    "authoritative": False,
                    "attempt_unit": "one native design request",
                    "requested_designs": 1,
                    "output_candidates": output_count,
                    "seed": int(attempt["seed"]),
                }
                rows = _rfd4_candidate_rows(child, state["campaign_id"], attempt, native_state)
            else:
                raise ValueError(f"No campaign reconciliation adapter is registered for {generator_id}.")
            receipt_id = hashlib.sha256(
                f"{child.name}:{metadata.get('completed_at') or result.get('completed_at') or result}".encode()
            ).hexdigest()
            ledger.record_receipt(
                attempt_id, receipt_id, candidates=rows, native_state=native_state,
                gpu_seconds=_worker_gpu_seconds(metadata),
            )
            continue
        if status in {"failed", "cancelled"}:
            retry_allowed = (
                not attempt.get("cluster_submission")
                and int(attempt.get("retry_count", 0)) < int(config["max_retries"])
            )
            if generator_id == "rfd4_proteina" and not _rfd4_retry_safe(child):
                retry_allowed = False
            if status == "failed" and retry_allowed:
                if attempt["status"] != "queued":
                    ledger.retry(attempt_id, child_run_dir=str(child))
                metadata = read_json(child / "metadata.json")
                if metadata.get("status") in {"failed", "cancelled", "paused", "holding"}:
                    resume_job(child, gpu_device=_gpu_for_attempt(config, int(attempt["ordinal"])))
                active.append({**attempt, "worker_status": "queued", "retrying": True})
            else:
                terminal_status = "interrupted" if metadata.get("worker_error") or status == "cancelled" else "failed"
                if attempt["status"] not in {"failed", "interrupted", "cancelled", "not_started"}:
                    ledger.transition(
                        attempt_id, terminal_status, child_run_dir=str(child),
                        failure=str(metadata.get("worker_error") or result.get("metrics", {}).get("worker_error") or status),
                    )
        elif attempt["status"] not in {"failed", "interrupted", "cancelled", "not_started"}:
            ledger.transition(attempt_id, "interrupted", child_run_dir=str(child), failure=f"Unexpected child status {status}")
    children = _attempt_child_paths(ledger)
    state = ledger.snapshot()
    attempt_order = {
        attempt_id: attempt["ordinal"]
        for attempt_id, attempt in state["attempts"].items()
    }
    candidate_rows = [
        dict(record["candidate"])
        for _candidate_id, record in sorted(
            state["candidates"].items(),
            key=lambda item: (attempt_order.get(item[1].get("attempt_id"), 0), item[0]),
        )
        if isinstance(record.get("candidate"), dict)
    ]
    missing_attempts = {
        str(record.get("attempt_id") or "")
        for record in state["candidates"].values()
        if not isinstance(record.get("candidate"), dict)
    }
    known_ids = {str(row.get("candidate_id") or "") for row in candidate_rows}
    for attempt_id in sorted(missing_attempts, key=lambda key: attempt_order.get(key, 0)):
        attempt = state["attempts"].get(attempt_id, {})
        child = children.get(attempt_id)
        if attempt.get("status") != "completed" or child is None:
            continue
        if generator_id == "bindcraft2":
            recovered = _native_candidate_rows(child, state["campaign_id"], attempt, attempt.get("native_state") or {})
        elif generator_id == "rfd4_proteina":
            recovered = _rfd4_candidate_rows(child, state["campaign_id"], attempt, attempt.get("native_state") or {})
        else:
            raise ValueError(f"No campaign recovery adapter is registered for {generator_id}.")
        candidate_rows.extend(row for row in recovered if str(row.get("candidate_id") or "") not in known_ids)
    return children, active, candidate_rows


def _sync_parent_results(
    parent_run_dir: Path,
    ledger: CampaignLedger,
    children: dict[str, Path],
    candidates: list[dict[str, Any]],
) -> None:
    import copy

    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import write_candidates
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json, write_json

    state = ledger.snapshot()
    metadata = read_json(parent_run_dir / "metadata.json")
    old_ids = metadata.get("distributed_candidate_ids") or []
    candidate_ids = [str(row.get("candidate_id") or "") for row in candidates]
    enriched = []
    for candidate in candidates:
        row = copy.deepcopy(candidate)
        record = state["candidates"].get(str(row.get("candidate_id") or ""), {})
        if record:
            raw = row.setdefault("raw_metadata", {})
            raw["evidence_candidates"] = copy.deepcopy(record.get("validation_evidence") or [])
            raw.setdefault("distributed_campaign", {}).update({
                "native_acceptance": record.get("native_acceptance"),
                "validated_acceptance": record.get("validated_acceptance"),
                "diversity_acceptance": record.get("diversity_acceptance"),
                "decision_reason": record.get("decision_reason", ""),
                "validation_decision_reason": record.get("validation_decision_reason", ""),
                "diversity_decision_reason": record.get("diversity_decision_reason", ""),
                "attempt_id": record.get("attempt_id"),
            })
        enriched.append(row)
    candidate_digest = hashlib.sha256(json.dumps(enriched, sort_keys=True, default=str).encode()).hexdigest()
    if (
        old_ids != candidate_ids
        or metadata.get("distributed_candidate_digest") != candidate_digest
        or not (parent_run_dir / "artifacts" / "normalized_candidates" / "campaign_result.json").is_file()
    ):
        write_candidates(parent_run_dir, str(state["config"].get("generator_id") or "bindcraft2"), enriched)
        metadata["distributed_candidate_ids"] = candidate_ids
        metadata["distributed_candidate_digest"] = candidate_digest
    attempts = []
    for attempt_id, attempt in sorted(state["attempts"].items(), key=lambda item: item[1]["ordinal"]):
        child = children.get(attempt_id)
        child_metadata = read_json(child / "metadata.json") if child else {}
        attempts.append({
            "attempt_id": attempt_id,
            "batch_id": attempt["batch_id"],
            "seed": attempt["seed"],
            "status": attempt["status"],
            "retry_count": attempt.get("retry_count", 0),
            "child_run_dir": str(child) if child else attempt.get("child_run_dir"),
            "worker": child_metadata.get("scheduler_worker_id") or child_metadata.get("worker_id") or "",
            "resource_allocation": child_metadata.get("resource_allocation") or {},
        })
    summary = {
        "schema_version": 1,
        "campaign_id": state["campaign_id"],
        "generator_id": state["config"].get("generator_id") or "bindcraft2",
        "acceptance_level": state["config"]["acceptance_level"],
        "accepted_target": state["config"]["accepted_target"],
        "termination_reason": state.get("termination_reason"),
        "generation_budget_exhausted": bool(
            state.get("generation_budget_exhausted")
            or len(state["attempts"]) >= int(state["config"]["max_attempts"])
        ),
        "counters": ledger.counters(),
        "gpu_compute_seconds": ledger.gpu_compute_seconds(),
        "attempts": attempts,
        "candidate_count": len(enriched),
        "candidates_jsonl": "artifacts/normalized_candidates/candidates.jsonl",
    }
    write_json(parent_run_dir / "artifacts" / "distributed_campaign" / "summary.json", summary)
    metadata["distributed_campaign_summary"] = {
        "campaign_id": state["campaign_id"],
        "termination_reason": state.get("termination_reason"),
        "generation_budget_exhausted": summary["generation_budget_exhausted"],
        "counters": summary["counters"],
        "candidate_count": len(candidates),
    }
    write_json(parent_run_dir / "metadata.json", metadata)


def _parent_progress(
    parent_run_dir: Path,
    ledger: CampaignLedger,
    config: dict[str, Any],
    active: list[dict[str, Any]],
) -> None:
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import update_status
    from mn_molecular_design_studio.engine.workflow_engine.workflows.generator_capabilities import generator_capabilities

    counters = ledger.counters()
    capability = generator_capabilities(str(config.get("generator_id") or "bindcraft2"))
    selected = _accepted_count(ledger, config)
    state = ledger.snapshot()
    attempts_remaining = max(0, int(config["max_attempts"]) - counters["attempts_authorized"])
    validation_pending = _validation_pending_candidates(state, config)
    active_ids = [str(item.get("attempt_id") or "") for item in active]
    label = (
        f"{capability.label} · {selected}/{config['accepted_target']} accepted · "
        f"{counters['attempts_started']}/{config['max_attempts']} attempts started · "
        f"{len(active_ids)} active"
    )
    if state.get("generation_budget_exhausted") and validation_pending:
        label += f" · generation budget consumed · {validation_pending} awaiting validation"
    update_status(
        parent_run_dir,
        "running",
        current_phase=label,
        current_engine=f"{capability.label} acceptance-driven campaign",
        progress_label=label,
        campaign_progress={
            **counters,
            "accepted_target": config["accepted_target"],
            "acceptance_level": config["acceptance_level"],
            "attempts_remaining": attempts_remaining,
            "generation_budget_exhausted": bool(state.get("generation_budget_exhausted")),
            "validation_pending_candidates": validation_pending,
            "active_attempt_ids": active_ids,
            "gpu_compute_seconds": _campaign_compute_seconds(
                ledger, _attempt_child_paths(ledger), active,
            ),
        },
    )


def _dispatch_attempts(
    parent_run_dir: Path,
    ledger: CampaignLedger,
    config: dict[str, Any],
    *,
    active_count: int,
    accepted_count: int,
) -> None:
    state = ledger.snapshot()
    max_concurrent = int(config["batch_size"]) * int(config["max_concurrent_batches"])
    slots = max(0, max_concurrent - active_count)
    target_remaining = max(0, int(config["accepted_target"]) - accepted_count - active_count)
    authorized_remaining = max(0, int(config["max_attempts"]) - len(state["attempts"]))
    requested = min(slots, target_remaining, authorized_remaining)
    if requested <= 0:
        return
    batch_base = (int(state["next_ordinal"]) - 1) // int(config["batch_size"]) + 1
    for offset in range(requested):
        batch_number = batch_base + offset // int(config["batch_size"])
        batch_id = f"b{batch_number:06d}"
        attempts = ledger.authorize(1, batch_id=batch_id)
        if not attempts:
            return
        attempt = attempts[0]
        try:
            _create_attempt_child(parent_run_dir, ledger, attempt, config)
        except Exception as exc:
            failure = f"{type(exc).__name__}: {exc}"
            ledger.transition(
                attempt["attempt_id"], "failed",
                failure=failure,
            )
            child_text = ledger.snapshot()["attempts"][attempt["attempt_id"]].get("child_run_dir")
            if child_text:
                from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json, write_json

                child = Path(child_text)
                metadata = read_json(child / "metadata.json")
                if metadata.get("status") in {"queued", "preparing"}:
                    metadata.update(status="failed", worker_error=failure, completed_at=utc_now())
                    write_json(child / "metadata.json", metadata)


def _campaign_elapsed_seconds(state: dict[str, Any]) -> float:
    from datetime import datetime, timezone

    created = str(state.get("created_at") or "")
    try:
        stamp = datetime.fromisoformat(created.replace("Z", "+00:00"))
        if stamp.tzinfo is None:
            stamp = stamp.replace(tzinfo=timezone.utc)
        return max(0.0, time.time() - stamp.timestamp())
    except ValueError:
        return 0.0


def _timestamp_seconds(value: Any) -> float | None:
    from datetime import datetime

    try:
        text = str(value or "")
        if not text:
            return None
        stamp = datetime.fromisoformat(text.replace("Z", "+00:00"))
        return stamp.timestamp()
    except (TypeError, ValueError, OverflowError):
        return None


def _worker_gpu_seconds(metadata: dict[str, Any], *, include_active: bool = False) -> float:
    allocation = metadata.get("resource_allocation") if isinstance(metadata.get("resource_allocation"), dict) else {}
    start = _timestamp_seconds(allocation.get("allocated_at") or metadata.get("worker_started_at"))
    if start is None:
        return 0.0
    end_text = metadata.get("completed_at") or metadata.get("updated_at")
    end = _timestamp_seconds(end_text)
    if include_active or end is None:
        end = time.time()
    device_text = str(allocation.get("gpu_device") or "")
    devices = allocation.get("gpu_ids") or allocation.get("gpu_devices")
    gpu_count = len(devices) if isinstance(devices, list) and devices else len([
        part for part in device_text.replace("device=", "").split(",") if part.strip() and part.strip().lower() not in {"cpu", "none"}
    ])
    if gpu_count == 0:
        return 0.0
    return max(0.0, end - start) * gpu_count


def _campaign_compute_seconds(
    ledger: CampaignLedger, children: dict[str, Path], active: list[dict[str, Any]],
) -> float:
    total = ledger.gpu_compute_seconds()
    for attempt in active:
        child = children.get(str(attempt.get("attempt_id") or ""))
        if child is None:
            continue
        from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json

        metadata = read_json(child / "metadata.json")
        if metadata.get("status") in {"queued", "preparing", "running"}:
            total += _worker_gpu_seconds(metadata, include_active=True)
    state = ledger.snapshot()
    for candidate in state["candidates"].values():
        attempt_id = str(candidate.get("attempt_id") or "")
        for stage_index, stage in (candidate.get("validation_stages") or {}).items():
            run_text = str(stage.get("run_dir") or "")
            if not run_text:
                continue
            run_dir = Path(run_text).expanduser()
            if not run_dir.is_dir():
                continue
            from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json

            metadata = read_json(run_dir / "metadata.json")
            if metadata.get("status") in {"queued", "preparing", "running"}:
                total += _worker_gpu_seconds(metadata, include_active=True)
    return total


def _cancel_queued_attempts(ledger: CampaignLedger, children: dict[str, Path], reason: str) -> None:
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json, stop_job

    state = ledger.snapshot()
    for attempt_id, attempt in state["attempts"].items():
        if attempt["status"] not in {"authorized", "queued"}:
            continue
        child = children.get(attempt_id)
        if child is not None:
            metadata = read_json(child / "metadata.json")
            if attempt.get("cluster_submission"):
                _cancel_campaign_attempt(ledger, attempt, child, reason)
            elif metadata.get("status") in {"queued", "preparing"}:
                stop_job(child, reason)
        ledger.transition(attempt_id, "not_started", stop_reason=reason)


def _cancel_campaign_attempt(
    ledger: CampaignLedger,
    attempt: dict[str, Any],
    child: Path,
    reason: str,
) -> None:
    """Cancel the originating peer job when an attempt has been dispatched remotely."""
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json, write_json, stop_job

    peer = attempt.get("cluster_submission") if isinstance(attempt.get("cluster_submission"), dict) else {}
    if peer.get("cluster_job_id"):
        from mn_molecular_design_studio.cluster import candidate_batches

        parent = Path(child).parent
        # Attempt runs live below the same results root as their campaign
        # parent, recorded durably in their own metadata.
        metadata = read_json(child / "metadata.json")
        parent_text = str(metadata.get("parent_run_dir") or "")
        if parent_text:
            parent = Path(parent_text)
        owner_id, project_id, development_mode = _campaign_actor(parent)
        candidate_batches.cancel(
            peer, owner_id=owner_id, project_id=project_id,
            allow_development=development_mode,
        )
        current = read_json(child / "metadata.json")
        current.update(status="cancelled", cancel_reason=reason, completed_at=utc_now())
        write_json(child / "metadata.json", current)
        return
    metadata = read_json(child / "metadata.json")
    if metadata.get("status") in {"queued", "preparing", "running"}:
        stop_job(child, reason)


def _cancel_validation_runs(ledger: CampaignLedger, reason: str) -> None:
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json, stop_job

    state = ledger.snapshot()
    for candidate_id, candidate in state["candidates"].items():
        stages = dict(candidate.get("validation_stages") or {})
        changed = False
        for stage in stages.values():
            run_text = str(stage.get("run_dir") or "")
            if not run_text:
                continue
            run_dir = Path(run_text).expanduser()
            if not run_dir.is_dir():
                continue
            status = str(read_json(run_dir / "metadata.json").get("status") or "")
            if status in {"queued", "preparing", "running"}:
                stop_job(run_dir, reason)
            if stage.get("status") in {"submitting", "queued", "preparing", "running"}:
                stage["status"] = "cancelled"
                stage["cancel_reason"] = reason
                changed = True
        if changed:
            ledger.update_candidate(candidate_id, validation_stages=stages)


def _finish_campaign(parent_run_dir: Path, ledger: CampaignLedger, candidates: list[dict[str, Any]]) -> Path:
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import finish_job

    state = ledger.snapshot()
    counters = ledger.counters()
    accepted = _accepted_count(ledger, state["config"])
    target_met = accepted >= int(state["config"]["accepted_target"])
    if target_met and state.get("termination_reason") == "budget_exhausted":
        state = ledger.terminate(
            "target_reached",
            detail=f"Reached {accepted} {state['config']['acceptance_level']} accepted candidates.",
        )
    _sync_parent_results(parent_run_dir, ledger, _attempt_child_paths(ledger), candidates)
    state = ledger.snapshot()
    reason = str(state.get("termination_reason") or "failed")
    generation_budget_exhausted = bool(
        state.get("generation_budget_exhausted")
        or len(state["attempts"]) >= int(state["config"]["max_attempts"])
    )
    success = reason in {"target_reached", "budget_exhausted", "runtime_exhausted", "compute_exhausted"}
    result = {
        "outputs": {
            "candidates": candidates,
            "candidates_jsonl": "artifacts/normalized_candidates/candidates.jsonl",
            "campaign_result": "artifacts/normalized_candidates/campaign_result.json",
            "distributed_campaign_summary": "artifacts/distributed_campaign/summary.json",
        },
        "metrics": {
            **counters,
            "generator_id": state["config"].get("generator_id") or "bindcraft2",
            "accepted_target": state["config"]["accepted_target"],
            "acceptance_level": state["config"]["acceptance_level"],
            "campaign_success": target_met,
            "generation_budget_exhausted": generation_budget_exhausted,
            "attempts_remaining": max(0, int(state["config"]["max_attempts"]) - counters["attempts_authorized"]),
            "termination_reason": reason,
            "target_met": target_met,
            "gpu_compute_seconds": ledger.gpu_compute_seconds(),
        },
        "downstream_artifacts": {
            "campaign_ledger": "artifacts/distributed_campaign/ledger.json",
            "attempt_summary": "artifacts/distributed_campaign/summary.json",
            "native_worker_runs": [item.get("child_run_dir") for item in state["attempts"].values() if item.get("child_run_dir")],
        },
    }
    finish_job(parent_run_dir, success, result)
    return parent_run_dir


def _validation_tagged_run(campaign_id: str, candidate_id: str, stage_index: int) -> Path | None:
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json
    from mn_molecular_design_studio.engine.workflow_engine.runtime import runs_root

    root = runs_root()
    for metadata_path in root.glob("*/*/metadata.json"):
        metadata = read_json(metadata_path)
        if (
            metadata.get("distributed_campaign_id") == campaign_id
            and metadata.get("distributed_candidate_id") == candidate_id
            and metadata.get("distributed_validation_stage") == stage_index
        ):
            return metadata_path.parent
    return None


def _validation_request(
    parent_run_dir: Path,
    candidate: dict[str, Any],
    stage: dict[str, Any],
    stage_index: int,
    campaign_id: str,
    previous_run: Path | None,
) -> tuple[dict[str, Any], dict[str, Any]]:
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import (
        _resolve_candidate_structure_path, candidates_jsonl_path,
    )

    workflow = str(stage["workflow"])
    parameters = dict(stage.get("parameters") or {})
    stage_resources = dict(stage.get("resources") or {})
    cpu_cores = int(stage_resources.get("cpu_cores") or 4)
    gpu_device = str(stage_resources.get("gpu_device") or "auto")
    if workflow == "refolding.optimized_prediction":
        from Bio.PDB import MMCIFIO, Select
        from mn_molecular_design_studio.engine.workflow_engine.workflows.optimized_prediction import ENGINE_PROFILES

        engine = str(parameters.pop("engine", "boltz2"))
        if engine not in ENGINE_PROFILES:
            raise ValueError(f"Unsupported distributed validation prediction engine: {engine}.")
        if "tasks" in parameters or "validation_sources" in parameters:
            raise ValueError("Distributed optimized-prediction stages construct their candidate inputs automatically.")
        structure = _resolve_candidate_structure_path(parent_run_dir, candidate)
        if structure is None:
            raise RuntimeError(f"Candidate {candidate.get('candidate_id')} has no readable complex structure.")
        from mn_molecular_design_studio.design_validation import protein_chains, read_structure, sequence

        chains = protein_chains(read_structure(structure))
        target_ids = [str(value) for value in candidate.get("target_chains") or []]
        binder_ids = [str(value) for value in candidate.get("binder_chains") or []]
        if not target_ids or not binder_ids or any(value not in chains for value in target_ids + binder_ids):
            raise RuntimeError("Candidate validation requires explicit target and binder protein chains.")
        binder_sequence = ":".join(sequence(chains[value]) for value in binder_ids)
        expected_binder = ":".join(str(candidate.get("binder_sequence") or "").split(":"))
        if binder_sequence.replace(":", "") != expected_binder.replace(":", "").replace(" ", ""):
            raise RuntimeError("Candidate binder sequence differs from its normalized structure.")
        target_sequences = ":".join(sequence(chains[value]) for value in target_ids)
        chain_ids = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
        if len(target_ids) + len(binder_ids) > len(chain_ids):
            raise ValueError("Distributed validation supports at most 52 protein chains.")
        target_use_msa = bool(parameters.pop("target_use_msa", True))
        use_target_template = bool(parameters.pop("use_target_template", False))
        templates = []
        if use_target_template and ENGINE_PROFILES[engine].get("supports_templates"):
            structure_obj = read_structure(structure)
            for index, target_id in enumerate(target_ids):
                class TargetChain(Select):
                    def accept_chain(self, chain):
                        return chain.id == target_id

                stream = __import__("io").StringIO()
                writer = MMCIFIO()
                writer.set_structure(structure_obj)
                writer.save(stream, TargetChain())
                templates.append({
                    "name": f"target_{chain_ids[index]}.cif", "text": stream.getvalue(),
                    "query_chain_id": chain_ids[index], "template_chain_id": target_id,
                })
        options = parameters.pop("engine_options", {})
        if not isinstance(options, dict):
            raise ValueError("Validation engine_options must be an object.")
        optimization_mode = str(parameters.pop("optimization_mode", "off"))
        seed = int(parameters.pop("seed", CampaignLedger._seed(
            campaign_id,
            int((candidate.get("raw_metadata") or {}).get("distributed_campaign", {}).get("attempt_seed") or 1),
            stage_index + 1,
        )))
        samples = int(parameters.pop("samples", 1))
        num_recycles = int(parameters.pop("num_recycles", 10))
        accept_terms = bool(parameters.pop("accept_alphafold3_terms", False))
        if parameters:
            raise ValueError(f"Unsupported optimized-prediction validation settings: {', '.join(sorted(parameters))}.")
        parameters = {
            "engine": engine,
            "optimization_mode": optimization_mode,
            "job_name": f"campaign-{campaign_id[:8]}-{candidate['candidate_id'][-12:]}-stage-{stage_index + 1}",
            "seed": seed,
            "samples": samples,
            "num_recycles": num_recycles,
            "cpu_cores": cpu_cores,
            "accept_alphafold3_terms": accept_terms,
            "use_target_template": use_target_template,
            "target_use_msa": target_use_msa,
            "engine_options": options,
            "tasks": [{
                "name": str(candidate["candidate_id"]),
                "target_protein_sequences": target_sequences,
                "binder_sequence": binder_sequence,
                "source_candidate_id": str(candidate["candidate_id"]),
                "source_run_dir": str(parent_run_dir),
                "structure_templates": templates,
                **({"molecule_chain_options": [
                    {"id": chain_ids[index], "use_msa": target_use_msa}
                    for index, _ in enumerate(target_sequences.split(":"))
                ] + [{
                    "id": chain_ids[len(target_sequences.split(":")) + index],
                    "use_msa": False,
                } for index, _ in enumerate(binder_sequence.split(":"))]} if engine == "boltz2" else {}),
            }],
        }
        request_resources = {"cpu_cores": cpu_cores, "gpu_device": gpu_device}
        if "min_vram_gb" in stage_resources:
            request_resources["min_vram_gb"] = float(stage_resources["min_vram_gb"])
        return {"schema_version": 1, "workflow": workflow, "parameters": parameters, "resources": request_resources}, stage_resources

    if workflow == "analysis.ipsae":
        if previous_run is None:
            raise ValueError("ipSAE validation requires a preceding complex-refolding prediction stage.")
        parameters.pop("cpu_cores", None)
        parameters.pop("gpu_device", None)
        parameters["source_run_dir"] = str(previous_run)
        parameters["candidates_jsonl"] = str(candidates_jsonl_path(previous_run))
        request_resources = {"cpu_cores": cpu_cores, "gpu_device": "none"}
        return {"schema_version": 1, "workflow": workflow, "parameters": parameters, "resources": request_resources}, stage_resources

    raise ValueError(f"Unsupported validation workflow: {workflow}.")


def _submit_validation_stage(
    parent_run_dir: Path,
    ledger: CampaignLedger,
    candidate: dict[str, Any],
    stage: dict[str, Any],
    stage_index: int,
    previous_run: Path | None,
) -> Path:
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import (
        read_json, submission_tags,
    )
    from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import submit_workflow_request

    state = ledger.snapshot()
    campaign_id = state["campaign_id"]
    candidate_id = str(candidate["candidate_id"])
    existing = _validation_tagged_run(campaign_id, candidate_id, stage_index)
    if existing is not None:
        return existing
    payload, _ = _validation_request(
        parent_run_dir, candidate, stage, stage_index, campaign_id, previous_run,
    )
    requests_dir = parent_run_dir / "artifacts" / "distributed_campaign" / "requests"
    requests_dir.mkdir(parents=True, exist_ok=True)
    request_path = requests_dir / f"{hashlib.sha256(f'{candidate_id}:{stage_index}'.encode()).hexdigest()[:24]}.json"
    request_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    with submission_tags(
        distributed_campaign_id=campaign_id,
        distributed_candidate_id=candidate_id,
        distributed_validation_stage=stage_index,
        parent_campaign_run_dir=str(parent_run_dir),
    ):
        result = submit_workflow_request(request_path)
    run_dir = Path(str(result["job"]["run_dir"])).expanduser().resolve()
    return run_dir


def _validation_evidence_row(run_dir: Path, source_candidate_id: str, sequence: str) -> dict[str, Any]:
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import candidates_jsonl_path, read_candidates
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json

    from mn_molecular_design_studio.output_vocabulary import candidate_evidence

    candidates = read_candidates(candidates_jsonl_path(run_dir))
    matching = [row for row in candidates if (
        any(record["candidate_id"] == source_candidate_id and record["same_sequence"] for record in candidate_evidence(row))
        and str(row.get("binder_sequence") or "").replace(":", "") == sequence.replace(":", "")
    )]
    if len(matching) != 1:
        raise RuntimeError(
            f"Validation run {run_dir.name} returned {len(matching)} same-sequence descendants for "
            f"candidate {source_candidate_id}; expected exactly one."
        )
    candidate = matching[0]
    result = read_json(run_dir / "result.json")
    if result.get("success") is not True:
        raise RuntimeError(f"Validation run {run_dir.name} did not complete successfully.")
    if str(candidate.get("stage") or "") not in {"complex_refolding", "analysis"}:
        raise RuntimeError(f"Validation run {run_dir.name} returned an incompatible candidate stage.")
    records = candidate_evidence(candidate)
    if not any(
        row["stage"] == "complex_refolding"
        and row["same_sequence"]
        and (row["candidate"].get("raw_metadata") or {}).get("identity_check", {}).get("sequence_chain_identity") == "pass"
        for row in records
    ):
        raise RuntimeError("Validation lineage has no sequence-verified complex prediction for the same binder.")
    return candidate


def _stage_selectors(candidate: dict[str, Any], stage_output: dict[str, Any], stage_config: dict[str, Any]) -> dict[str, Any]:
    from mn_molecular_design_studio.candidate_selection import evaluate_selection

    stage_name = str(stage_output.get("stage") or "")
    tool = str(stage_output.get("source_tool") or stage_output.get("tool") or "")
    output_id = str(stage_output.get("candidate_id") or "")
    fields = [str(value) for value in stage_config.get("required_metrics") or []]
    rules = []
    for rule in stage_config.get("rules") or []:
        rules.append({**rule, "field": str(rule["field"])})
    selectors = [
        field if field.startswith("evidence/") else f"evidence/{stage_name}/{tool}/{output_id}/{field}"
        for field in fields
    ]
    for rule in rules:
        field = str(rule["field"])
        rule["field"] = field if field.startswith("evidence/") else f"evidence/{stage_name}/{tool}/{output_id}/{field}"
    if not selectors and not rules:
        return {"status": "pass", "selected": False, "reasons": [], "measurements": {}, "policy": {}}
    selection = evaluate_selection(
        [stage_output],
        {"required_metrics": selectors, "rules": rules},
    )[0]["raw_metadata"]["selection"]
    return selection


def _validation_active_count(state: dict[str, Any]) -> int:
    count = 0
    for candidate in state["candidates"].values():
        for stage in (candidate.get("validation_stages") or {}).values():
            if stage.get("status") in {"queued", "preparing", "running"}:
                count += 1
    return count


def _validation_pending_candidates(state: dict[str, Any], config: dict[str, Any]) -> int:
    if not config.get("validation_stages"):
        return 0
    pending = 0
    for candidate in state["candidates"].values():
        if candidate.get("native_acceptance") not in {"accepted", "unavailable"}:
            continue
        if candidate.get("validated_acceptance") in {"accepted", "rejected", "failed"}:
            continue
        pending += 1
    return pending


def _mark_validation_pending(ledger: CampaignLedger, candidate_id: str) -> None:
    """Keep requeued validation visible as pending until its evidence is reconciled."""
    ledger.set_candidate_acceptance(candidate_id, validated="pending")
    ledger.update_candidate(candidate_id, validation_decision_reason="", decision_reason="")


def _drive_validation(
    parent_run_dir: Path,
    ledger: CampaignLedger,
    candidates: list[dict[str, Any]],
    config: dict[str, Any],
    *,
    dispatch: bool = True,
) -> list[dict[str, Any]]:
    if not config.get("validation_stages"):
        if config["acceptance_level"] == "validated":
            raise ValueError("Validated acceptance has no configured downstream validation stage.")
        return candidates
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json

    candidate_rows = {str(row["candidate_id"]): row for row in candidates}
    state = ledger.snapshot()
    stage_specs = config["validation_stages"]
    active_count = _validation_active_count(state)
    for candidate_id, ledger_candidate in sorted(
        state["candidates"].items(),
        key=lambda item: (state["attempts"][item[1]["attempt_id"]]["ordinal"], item[0]),
    ):
        native = candidate_rows.get(candidate_id)
        if native is None or ledger_candidate.get("native_acceptance") not in {"accepted", "unavailable"}:
            continue
        stage_results = dict(ledger_candidate.get("validation_stages") or {})
        previous_run: Path | None = None
        failed_stage = False
        for stage_index, stage_config in enumerate(stage_specs):
            key = str(stage_index)
            stage_result = dict(stage_results.get(key) or {})
            run_dir = Path(str(stage_result.get("run_dir") or "")).expanduser() if stage_result.get("run_dir") else None
            if run_dir is None or not run_dir.exists():
                run_dir = _validation_tagged_run(state["campaign_id"], candidate_id, stage_index)
            if run_dir is not None:
                metadata = read_json(run_dir / "metadata.json")
                worker_status = str(metadata.get("status") or "unknown")
                if worker_status in {"queued", "preparing", "running"}:
                    _mark_validation_pending(ledger, candidate_id)
                    stage_result.update(status=worker_status, run_dir=str(run_dir))
                    stage_results[key] = stage_result
                    ledger.update_candidate(candidate_id, validation_stages=stage_results)
                    previous_run = run_dir
                    continue
                if worker_status in {"completed", "succeeded"}:
                    try:
                        ledger.record_compute(
                            str(ledger_candidate["attempt_id"]),
                            f"validation:{stage_index}:{run_dir.name}",
                            _worker_gpu_seconds(metadata),
                        )
                        evidence = _validation_evidence_row(
                            run_dir, candidate_id, str(native.get("binder_sequence") or ""),
                        )
                        selection = _stage_selectors(native, evidence, stage_config)
                        status = "accepted" if selection.get("status") == "pass" else "rejected"
                        stage_result.update(
                            status="completed", run_dir=str(run_dir), output_candidate_id=evidence["candidate_id"],
                            evidence_candidate=evidence, selection=selection,
                        )
                        recovered_error = stage_result.pop("error", None)
                        if recovered_error:
                            stage_result.setdefault("recovered_submission_errors", []).append(recovered_error)
                        ledger_candidate.setdefault("validation_evidence", [])
                        evidence_rows = [
                            row for row in ledger_candidate["validation_evidence"]
                            if row.get("candidate_id") != evidence["candidate_id"]
                        ]
                        evidence_rows.append(evidence)
                        ledger_candidate["validation_evidence"] = evidence_rows
                        stage_results[key] = stage_result
                        ledger.update_candidate(
                            candidate_id,
                            validation_stages=stage_results,
                            validation_evidence=evidence_rows,
                        )
                        if status == "rejected":
                            failed_stage = True
                            ledger.set_candidate_acceptance(
                                candidate_id, validated="rejected",
                                reason="; ".join(selection.get("reasons") or ["Validation criteria did not pass."]),
                            )
                            ledger.update_candidate(
                                candidate_id,
                                validation_decision_reason="; ".join(selection.get("reasons") or ["Validation criteria did not pass."]),
                            )
                            break
                        previous_run = run_dir
                        continue
                    except Exception as exc:
                        stage_result.update(status="failed", run_dir=str(run_dir), error=f"{type(exc).__name__}: {exc}")
                        stage_results[key] = stage_result
                        ledger.update_candidate(candidate_id, validation_stages=stage_results)
                        ledger.set_candidate_acceptance(
                            candidate_id, validated="failed",
                            reason=f"Validation stage failed: {type(exc).__name__}: {exc}",
                        )
                        ledger.update_candidate(candidate_id, validation_decision_reason=f"Validation stage failed: {type(exc).__name__}: {exc}")
                        failed_stage = True
                        break
                worker_error = str(metadata.get("worker_error") or worker_status)
                retry_count = int(stage_result.get("retry_count", 0))
                if worker_status == "failed" and retry_count < int(config["max_retries"]):
                    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import resume_job

                    stage_result.update(
                        status="queued",
                        run_dir=str(run_dir),
                        retry_count=retry_count + 1,
                        last_retry_error=worker_error,
                    )
                    stage_results[key] = stage_result
                    ledger.update_candidate(candidate_id, validation_stages=stage_results)
                    _mark_validation_pending(ledger, candidate_id)
                    try:
                        retry_gpu = str((stage_config.get("resources") or {}).get("gpu_device") or "auto")
                        resume_job(run_dir, gpu_device=retry_gpu)
                    except Exception as exc:
                        stage_result.update(status="failed", error=f"{type(exc).__name__}: {exc}")
                        stage_results[key] = stage_result
                        ledger.update_candidate(candidate_id, validation_stages=stage_results)
                        ledger.set_candidate_acceptance(
                            candidate_id, validated="failed",
                            reason=f"Validation retry failed: {type(exc).__name__}: {exc}",
                        )
                        ledger.update_candidate(
                            candidate_id,
                            validation_decision_reason=f"Validation retry failed: {type(exc).__name__}: {exc}",
                        )
                        failed_stage = True
                        break
                    active_count += 1
                    previous_run = run_dir
                    continue

                stage_result.update(status="failed", run_dir=str(run_dir), error=worker_error)
                ledger.record_compute(
                    str(ledger_candidate["attempt_id"]),
                    f"validation:{stage_index}:{run_dir.name}",
                    _worker_gpu_seconds(metadata),
                )
                stage_results[key] = stage_result
                ledger.update_candidate(candidate_id, validation_stages=stage_results)
                ledger.set_candidate_acceptance(
                    candidate_id, validated="failed", reason=f"Validation stage ended with status {worker_status}.",
                )
                ledger.update_candidate(candidate_id, validation_decision_reason=f"Validation stage ended with status {worker_status}.")
                failed_stage = True
                break

            if stage_result.get("status") == "completed":
                recovered_error = stage_result.pop("error", None)
                if recovered_error:
                    stage_result.setdefault("recovered_submission_errors", []).append(recovered_error)
                    stage_results[key] = stage_result
                    ledger.update_candidate(candidate_id, validation_stages=stage_results)
                previous_run = Path(str(stage_result["run_dir"]))
                continue
            if failed_stage:
                break
            if stage_index and (previous_run is None or (stage_results.get(str(stage_index - 1)) or {}).get("status") != "completed"):
                break
            if active_count >= int(config["max_concurrent_validations"]):
                break
            if not dispatch:
                break
            stage_result.update(status="submitting", workflow=stage_config["workflow"])
            stage_results[key] = stage_result
            ledger.update_candidate(candidate_id, validation_stages=stage_results)
            try:
                run_dir = _submit_validation_stage(
                    parent_run_dir, ledger, native, stage_config, stage_index, previous_run,
                )
                _mark_validation_pending(ledger, candidate_id)
                stage_result.update(status="queued", run_dir=str(run_dir))
                stage_results[key] = stage_result
                ledger.update_candidate(candidate_id, validation_stages=stage_results)
                active_count += 1
                previous_run = run_dir
            except Exception as exc:
                stage_result.update(status="failed", error=f"{type(exc).__name__}: {exc}")
                stage_results[key] = stage_result
                ledger.update_candidate(candidate_id, validation_stages=stage_results)
                ledger.set_candidate_acceptance(
                    candidate_id, validated="failed",
                    reason=f"Validation submission failed: {type(exc).__name__}: {exc}",
                )
                ledger.update_candidate(candidate_id, validation_decision_reason=f"Validation submission failed: {type(exc).__name__}: {exc}")
                failed_stage = True
                break

        latest = ledger.snapshot()["candidates"].get(candidate_id, {})
        latest_stages = latest.get("validation_stages") or {}
        if not failed_stage and len(latest_stages) == len(stage_specs) and all(
            (latest_stages.get(str(index)) or {}).get("status") == "completed"
            for index in range(len(stage_specs))
        ):
            ledger.set_candidate_acceptance(candidate_id, validated="accepted", reason="All configured validation stages and criteria passed.")
            ledger.update_candidate(candidate_id, validation_decision_reason="All configured validation stages and criteria passed.")
    updated = ledger.snapshot()
    by_id = {row.get("candidate_id"): row for row in candidates}
    for candidate_id, row in by_id.items():
        ledger_candidate = updated["candidates"].get(str(candidate_id), {})
        if ledger_candidate.get("native_acceptance") in {"accepted", "unavailable"}:
            raw = row.setdefault("raw_metadata", {})
            raw["evidence_candidates"] = ledger_candidate.get("validation_evidence") or []
            raw.setdefault("distributed_campaign", {})["validated_acceptance"] = ledger_candidate.get("validated_acceptance", "pending")
            if ledger_candidate.get("decision_reason"):
                raw["distributed_campaign"]["validation_decision_reason"] = ledger_candidate["decision_reason"]
    return candidates


def _sequence_identity(left: str, right: str) -> float:
    from Bio.Align import PairwiseAligner

    left = str(left or "").replace(":", "").upper()
    right = str(right or "").replace(":", "").upper()
    if not left or not right:
        return 0.0
    aligner = PairwiseAligner(mode="global", match_score=2, mismatch_score=-1, open_gap_score=-8, extend_gap_score=-0.5)
    alignment = aligner.align(left, right)[0]
    matches = sum(
        left[left_start + offset] == right[right_start + offset]
        for (left_start, left_end), (right_start, right_end) in zip(*alignment.aligned)
        for offset in range(min(left_end - left_start, right_end - right_start))
    )
    return matches / max(len(left), len(right))


def _binder_target_aligned_rmsd(parent_run_dir: Path, left: dict[str, Any], right: dict[str, Any]) -> float | None:
    from Bio.Align import PairwiseAligner
    from Bio.PDB import Superimposer
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import _resolve_candidate_structure_path
    from mn_molecular_design_studio.design_validation import protein_chains, read_structure

    left_path = _resolve_candidate_structure_path(parent_run_dir, left)
    right_path = _resolve_candidate_structure_path(parent_run_dir, right)
    if left_path is None or right_path is None:
        return None
    left_chains = protein_chains(read_structure(left_path))
    right_chains = protein_chains(read_structure(right_path))
    left_targets, right_targets = left.get("target_chains") or [], right.get("target_chains") or []
    left_binders, right_binders = left.get("binder_chains") or [], right.get("binder_chains") or []
    if len(left_targets) != len(right_targets) or len(left_binders) != len(right_binders):
        return None
    fixed, moving = [], []
    for left_id, right_id in zip(left_targets, right_targets):
        left_res, right_res = left_chains.get(str(left_id), []), right_chains.get(str(right_id), [])
        if len(left_res) != len(right_res) or not left_res:
            return None
        fixed.extend(res["CA"] for res in left_res if "CA" in res)
        moving.extend(res["CA"] for res in right_res if "CA" in res)
    if len(fixed) < 3 or len(fixed) != len(moving):
        return None
    superimposer = Superimposer()
    superimposer.set_atoms(fixed, moving)
    rotation, translation = superimposer.rotran
    left_binder_atoms, right_binder_atoms = [], []
    left_binder_sequence, right_binder_sequence = [], []
    for left_id, right_id in zip(left_binders, right_binders):
        left_res, right_res = left_chains.get(str(left_id), []), right_chains.get(str(right_id), [])
        if not left_res or not right_res:
            return None
        left_binder_atoms.extend(res["CA"] for res in left_res if "CA" in res)
        right_binder_atoms.extend(res["CA"] for res in right_res if "CA" in res)
        left_binder_sequence.extend(sequence([res]) for res in left_res if "CA" in res)
        right_binder_sequence.extend(sequence([res]) for res in right_res if "CA" in res)
    if len(left_binder_atoms) < 3 or len(right_binder_atoms) < 3:
        return None
    right_xyz = [atom.coord @ rotation + translation for atom in right_binder_atoms]
    aligner = PairwiseAligner(mode="global", match_score=2, mismatch_score=-1, open_gap_score=-8, extend_gap_score=-0.5)
    alignment = aligner.align("".join(left_binder_sequence), "".join(right_binder_sequence))[0]
    aligned_pairs = [
        (left_start + offset, right_start + offset)
        for (left_start, left_end), (right_start, right_end) in zip(*alignment.aligned)
        for offset in range(min(left_end - left_start, right_end - right_start))
    ]
    if len(aligned_pairs) < 3:
        return None
    return math.sqrt(sum(
        float(((left_binder_atoms[left_index].coord - right_xyz[right_index]) ** 2).sum())
        for left_index, right_index in aligned_pairs
    ) / len(aligned_pairs))


def _interface_cluster_key(candidate: dict[str, Any], field: str) -> str | None:
    from mn_molecular_design_studio.candidate_selection import evaluate_selection, metric_value
    from mn_molecular_design_studio.output_vocabulary import candidate_evidence

    decision = evaluate_selection([candidate], {"required_metrics": [field]})[0]["raw_metadata"]["selection"]
    if decision.get("status") != "pass":
        return None
    _, stage, tool, *rest = field.split("/")
    output_id = rest[-2] if len(rest) == 2 else ""
    records = [row for row in candidate_evidence(candidate) if row["same_sequence"]
               and row["stage"] == stage and row["source_tool"] == tool
               and (not output_id or row["candidate_id"] == output_id)]
    if len(records) != 1:
        return None
    item = records[0]["candidate"]
    raw = item.get("raw_metadata") or {}
    value = metric_value(candidate, field)
    return hashlib.sha256(json.dumps({
        "stage": stage, "tool": tool, "protocol": raw.get("analysis_parameters") or raw.get("prediction_protocol") or {},
        "cluster_id": raw.get("clustering_run_id"), "value": value,
    }, sort_keys=True, default=str).encode()).hexdigest()


def _apply_diversity(
    parent_run_dir: Path,
    ledger: CampaignLedger,
    candidates: list[dict[str, Any]],
    config: dict[str, Any],
    *,
    settled: bool = True,
) -> list[dict[str, Any]]:
    diversity = config.get("diversity") or {}
    if config["acceptance_level"] != "diversity_qualified" and not diversity:
        return candidates
    if not settled:
        return candidates
    state = ledger.snapshot()
    candidate_by_id = {str(row.get("candidate_id")): row for row in candidates}
    ordered = sorted(
        state["candidates"].items(),
        key=lambda item: (state["attempts"][item[1]["attempt_id"]]["ordinal"], item[0]),
    )
    use_validation = bool(config.get("validation_stages"))
    qualified: list[tuple[str, dict[str, Any], int]] = []
    clusters: list[list[tuple[str, dict[str, Any], int]]] = []
    for candidate_id, ledger_candidate in ordered:
        row = candidate_by_id.get(candidate_id)
        if row is None or ledger_candidate.get("native_acceptance") not in {"accepted", "unavailable"}:
            continue
        if use_validation and ledger_candidate.get("validated_acceptance") != "accepted":
            if ledger_candidate.get("validated_acceptance") in {"rejected", "failed"}:
                ledger.set_candidate_acceptance(candidate_id, diversity="rejected", reason=ledger_candidate.get("decision_reason", "Validation did not pass."))
                ledger.update_candidate(candidate_id, diversity_decision_reason=ledger_candidate.get("decision_reason", "Validation did not pass."))
            continue
        similar_cluster = None
        similarity_reason = None
        interface_key = None
        if diversity.get("interface_cluster_field"):
            interface_key = _interface_cluster_key(row, str(diversity["interface_cluster_field"]))
            if interface_key is None:
                ledger.set_candidate_acceptance(candidate_id, diversity="rejected", reason="Configured interface diversity evidence is missing or invalid.")
                ledger.update_candidate(candidate_id, diversity_decision_reason="Configured interface diversity evidence is missing or invalid.")
                continue
        for cluster_index, cluster in enumerate(clusters):
            for previous_id, previous_row, _ in cluster:
                if diversity.get("deduplicate_exact", True) and str(row.get("binder_sequence") or "").replace(":", "").replace(" ", "").upper() == str(previous_row.get("binder_sequence") or "").replace(":", "").replace(" ", "").upper():
                    similar_cluster, similarity_reason = cluster_index, "Exact sequence duplicate"
                    break
                threshold = diversity.get("sequence_identity_max")
                if threshold is not None and _sequence_identity(row.get("binder_sequence"), previous_row.get("binder_sequence")) > float(threshold):
                    similar_cluster, similarity_reason = cluster_index, f"Sequence identity exceeds {float(threshold):g}"
                    break
                structural_min = diversity.get("binder_ca_rmsd_min_angstrom")
                if structural_min is not None:
                    rmsd = _binder_target_aligned_rmsd(parent_run_dir, row, previous_row)
                    if rmsd is None:
                        similar_cluster, similarity_reason = cluster_index, "Configured structural diversity could not be evaluated"
                        break
                    if rmsd < float(structural_min):
                        similar_cluster, similarity_reason = cluster_index, f"Binder Cα RMSD {rmsd:.3f} Å is below {float(structural_min):g} Å"
                        break
                if interface_key is not None and interface_key == _interface_cluster_key(previous_row, str(diversity["interface_cluster_field"])):
                    similar_cluster, similarity_reason = cluster_index, "Interface/binding-mode cluster matches an earlier candidate"
                    break
            if similar_cluster is not None:
                break
        if similar_cluster is None:
            clusters.append([])
            similar_cluster = len(clusters) - 1
        cluster = clusters[similar_cluster]
        max_per_cluster = int(diversity.get("max_per_cluster") or (1 if diversity.get("sequence_identity_max") is not None or diversity.get("binder_ca_rmsd_min_angstrom") is not None or diversity.get("interface_cluster_field") or diversity.get("deduplicate_exact", True) else 1_000_000))
        if len(cluster) >= max_per_cluster:
            ledger.set_candidate_acceptance(
                candidate_id, diversity="rejected",
                reason=similarity_reason or f"Diversity cluster {similar_cluster + 1} reached its {max_per_cluster}-candidate quota.",
            )
            ledger.update_candidate(candidate_id, diversity_decision_reason=similarity_reason or f"Diversity cluster {similar_cluster + 1} reached its {max_per_cluster}-candidate quota.")
            continue
        cluster.append((candidate_id, row, similar_cluster))
        qualified.append((candidate_id, row, similar_cluster))
        ledger.set_candidate_acceptance(candidate_id, diversity="accepted", reason=f"Diversity cluster {similar_cluster + 1}.")
        ledger.update_candidate(candidate_id, diversity_decision_reason=f"Diversity cluster {similar_cluster + 1}.")
    state = ledger.snapshot()
    for candidate_id, row in candidate_by_id.items():
        ledger_candidate = state["candidates"].get(candidate_id, {})
        if ledger_candidate:
            row.setdefault("raw_metadata", {}).setdefault("distributed_campaign", {})["diversity_acceptance"] = ledger_candidate.get("diversity_acceptance", "pending")
            row["raw_metadata"]["distributed_campaign"]["diversity_decision_reason"] = ledger_candidate.get("decision_reason", "")
            row["raw_metadata"]["distributed_campaign"]["diversity_decision_reason"] = ledger_candidate.get("diversity_decision_reason", ledger_candidate.get("decision_reason", ""))
            row["raw_metadata"]["distributed_campaign"]["diversity_cluster"] = next(
                (index + 1 for index, cluster in enumerate(clusters) if any(item[0] == candidate_id for item in cluster)),
                None,
            )
    return candidates


def _run_distributed_campaign(parent_run_dir: Path, expected_generator_id: str) -> Path:
    """Run one durable campaign controller with engine-specific attempt adapters."""
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json, update_status
    from mn_molecular_design_studio.engine.workflow_engine.workflows.generator_capabilities import generator_capabilities

    parent_run_dir = Path(parent_run_dir).expanduser().resolve()
    parent_payload = read_json(parent_run_dir / "input.json")
    parent_params = parent_payload.get("params") if isinstance(parent_payload.get("params"), dict) else {}
    raw_config = parent_params.get("distributed_campaign")
    if not isinstance(raw_config, dict) or not raw_config.get("enabled"):
        raise ValueError("Distributed campaign settings are missing or disabled.")
    raw_config = {**raw_config, "generator_id": str(raw_config.get("generator_id") or expected_generator_id)}
    if raw_config["generator_id"] != expected_generator_id:
        raise ValueError(
            f"Campaign parent is for {expected_generator_id}, but settings select {raw_config['generator_id']}."
        )
    config = validate_campaign_config(raw_config)
    capability = generator_capabilities(expected_generator_id)
    ledger = CampaignLedger.create(_campaign_ledger_path(parent_run_dir), config)
    with ledger.controller_lock():
        state = ledger.snapshot()
        update_status(
            parent_run_dir, "running",
            distributed_campaign_id=state["campaign_id"],
            distributed_campaign_acceptance_level=config["acceptance_level"],
            distributed_campaign_generator=config["generator_id"],
        )
        while True:
            state = ledger.snapshot()
            if len(state["attempts"]) >= int(config["max_attempts"]):
                state = ledger.mark_generation_budget_exhausted()
            parent_metadata = read_json(parent_run_dir / "metadata.json")
            if state.get("termination_reason") == "failed":
                raise RuntimeError(state.get("termination_detail") or "Campaign controller previously failed.")
            if parent_metadata.get("status") == "cancelled" or state.get("termination_reason") == "cancelled":
                ledger.terminate("cancelled", detail="Campaign cancelled by user.")
                children = _attempt_child_paths(ledger)
                for attempt_id, attempt in ledger.snapshot()["attempts"].items():
                    child = children.get(attempt_id)
                    if child is None or attempt.get("status") not in {"authorized", "queued", "running"}:
                        continue
                    _cancel_campaign_attempt(
                        ledger, attempt, child, f"Parent {capability.label} campaign cancelled",
                    )
                    if attempt.get("status") in {"authorized", "queued", "running"}:
                        ledger.transition(attempt_id, "cancelled", stop_reason="Parent campaign cancelled")
                _cancel_validation_runs(ledger, f"Parent {capability.label} campaign cancelled")
                _, _, candidates = _reconcile_attempts(ledger, parent_run_dir, config)
                _sync_parent_results(parent_run_dir, ledger, children, candidates)
                return parent_run_dir

            children, active, candidates = _reconcile_attempts(ledger, parent_run_dir, config)
            # Downstream prediction workers resolve source_candidate_id from the
            # parent's normalized collection, so publish native receipts first.
            _sync_parent_results(parent_run_dir, ledger, children, candidates)
            current_state = ledger.snapshot()
            runtime_limit = config.get("max_runtime_seconds")
            compute_limit = config.get("max_compute_seconds")
            validation_dispatch_allowed = not (
                runtime_limit is not None and _campaign_elapsed_seconds(current_state) >= float(runtime_limit)
            ) and not (
                compute_limit is not None
                and _campaign_compute_seconds(ledger, children, active) >= float(compute_limit)
            )
            candidates = _drive_validation(
                parent_run_dir, ledger, candidates, config,
                dispatch=validation_dispatch_allowed,
            )
            validation_pending = _validation_pending_candidates(ledger.snapshot(), config)
            validation_active = _validation_active_count(ledger.snapshot())
            diversity_settled = not active and not validation_pending and not validation_active
            candidates = _apply_diversity(
                parent_run_dir, ledger, candidates, config, settled=diversity_settled,
            )
            _sync_parent_results(parent_run_dir, ledger, children, candidates)
            counters = ledger.counters()
            selected = _accepted_count(ledger, config)
            target_reached = selected >= int(config["accepted_target"])
            runtime_budget = config.get("max_runtime_seconds")
            runtime_exhausted = runtime_budget is not None and _campaign_elapsed_seconds(state) >= float(runtime_budget)
            compute_budget = config.get("max_compute_seconds")
            compute_spent = _campaign_compute_seconds(ledger, children, active)
            compute_exhausted = compute_budget is not None and compute_spent >= float(compute_budget)
            if target_reached:
                _cancel_queued_attempts(ledger, children, "Global accepted-design target reached")
                children, active, candidates = _reconcile_attempts(ledger, parent_run_dir, config)
                _sync_parent_results(parent_run_dir, ledger, children, candidates)
                candidates = _drive_validation(parent_run_dir, ledger, candidates, config)
                diversity_settled = not active and not _validation_pending_candidates(ledger.snapshot(), config) and not _validation_active_count(ledger.snapshot())
                candidates = _apply_diversity(
                    parent_run_dir, ledger, candidates, config, settled=diversity_settled,
                )
                _sync_parent_results(parent_run_dir, ledger, children, candidates)
                _parent_progress(parent_run_dir, ledger, config, active)
                pending_validation = _validation_pending_candidates(ledger.snapshot(), config)
                selected = _accepted_count(ledger, config)
                if not active and not pending_validation and selected >= int(config["accepted_target"]):
                    ledger.terminate("target_reached", detail=f"Reached {selected} {config['acceptance_level']} accepted candidates.")
                    return _finish_campaign(parent_run_dir, ledger, candidates)
            elif runtime_exhausted or compute_exhausted:
                _cancel_queued_attempts(
                    ledger, children,
                    "Campaign runtime budget reached" if runtime_exhausted else "Campaign compute budget reached",
                )
                children, active, candidates = _reconcile_attempts(ledger, parent_run_dir, config)
                _sync_parent_results(parent_run_dir, ledger, children, candidates)
                candidates = _drive_validation(
                    parent_run_dir, ledger, candidates, config, dispatch=False,
                )
                diversity_settled = not active and not _validation_active_count(ledger.snapshot())
                candidates = _apply_diversity(
                    parent_run_dir, ledger, candidates, config, settled=diversity_settled,
                )
                _sync_parent_results(parent_run_dir, ledger, children, candidates)
                _parent_progress(parent_run_dir, ledger, config, active)
                if not active and not _validation_active_count(ledger.snapshot()):
                    reason = "runtime_exhausted" if runtime_exhausted else "compute_exhausted"
                    detail = "Configured runtime budget expired." if runtime_exhausted else "Configured GPU compute budget expired."
                    ledger.terminate(reason, detail=detail)
                    return _finish_campaign(parent_run_dir, ledger, candidates)
            else:
                _dispatch_attempts(
                    parent_run_dir, ledger, config,
                    active_count=len(active) + _validation_pending_candidates(ledger.snapshot(), config),
                    accepted_count=selected,
                )
                # Newly created child jobs are scheduler-owned immediately,
                # even before the scheduler records their first status. Reconcile
                # them before deciding that the attempt budget is exhausted.
                children, active, candidates = _reconcile_attempts(ledger, parent_run_dir, config)
                _sync_parent_results(parent_run_dir, ledger, children, candidates)
                candidates = _drive_validation(parent_run_dir, ledger, candidates, config)
                _sync_parent_results(parent_run_dir, ledger, children, candidates)
                state = ledger.snapshot()
                if len(state["attempts"]) >= int(config["max_attempts"]):
                    state = ledger.mark_generation_budget_exhausted()
                if (
                    len(state["attempts"]) >= int(config["max_attempts"])
                    and not active
                ):
                    pending_validation = _validation_pending_candidates(state, config)
                    selected = _accepted_count(ledger, config)
                    if not pending_validation and selected < int(config["accepted_target"]):
                        ledger.terminate("budget_exhausted", detail=f"Attempt budget exhausted at {selected}/{config['accepted_target']} accepted.")
                        _, _, candidates = _reconcile_attempts(ledger, parent_run_dir, config)
                        _sync_parent_results(parent_run_dir, ledger, _attempt_child_paths(ledger), candidates)
                        return _finish_campaign(parent_run_dir, ledger, candidates)
                _parent_progress(parent_run_dir, ledger, config, active)
            time.sleep(float(config["poll_interval_seconds"]))


def _run_distributed_campaign_with_failure_state(parent_run_dir: Path, expected_generator_id: str) -> Path:
    parent_run_dir = Path(parent_run_dir).expanduser().resolve()
    try:
        return _run_distributed_campaign(parent_run_dir, expected_generator_id)
    except CampaignControllerAlreadyRunning:
        return parent_run_dir
    except Exception as exc:
        path = _campaign_ledger_path(parent_run_dir)
        if path.is_file():
            try:
                CampaignLedger(path).terminate("failed", detail=f"{type(exc).__name__}: {exc}")
            except Exception:
                pass
        raise


def run_bindcraft2_distributed_campaign(parent_run_dir: Path) -> Path:
    return _run_distributed_campaign_with_failure_state(parent_run_dir, "bindcraft2")


def run_rfd4_proteina_distributed_campaign(parent_run_dir: Path) -> Path:
    return _run_distributed_campaign_with_failure_state(parent_run_dir, "rfd4_proteina")


def cancel_distributed_campaign(parent_run_dir: Path, reason: str = "Campaign cancelled by user") -> None:
    """Cancel only children owned by this campaign; the controller remains origin-bound."""
    path = _campaign_ledger_path(Path(parent_run_dir))
    if not path.is_file():
        return
    ledger = CampaignLedger(path)
    ledger.terminate("cancelled", detail=reason)
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json, stop_job

    children = _attempt_child_paths(ledger)
    for attempt_id, child in children.items():
        attempt = ledger.snapshot()["attempts"].get(attempt_id, {})
        status = str(read_json(child / "metadata.json").get("status") or "")
        if status in {"queued", "preparing", "running"}:
            _cancel_campaign_attempt(ledger, attempt, child, reason)
        if attempt.get("status") in {"authorized", "queued", "running"}:
            ledger.transition(attempt_id, "cancelled", stop_reason=reason)
    _cancel_validation_runs(ledger, reason)


def cancel_bindcraft2_distributed_campaign(parent_run_dir: Path, reason: str = "Campaign cancelled by user") -> None:
    """Backward-compatible BindCraft2 campaign cancellation entry point."""
    cancel_distributed_campaign(parent_run_dir, reason)
