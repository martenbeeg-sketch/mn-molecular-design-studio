"""Schema-v2 execution accounting. No scheduler, native process, or inference calls.

Explicit multi-member reservations are accounting-only. The production adapter
reserves single-member jobs through the legacy authorization API.
"""
from __future__ import annotations

from copy import deepcopy
import hashlib
import json
import math

from mn_molecular_design_studio.engine.workflow_engine.core.jobs import utc_now


class LedgerConflict(ValueError):
    """An identity was replayed with different immutable content."""


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def attempt_seed(campaign_id, configured_seed, ordinal):
    raw = hashlib.sha256(f"{campaign_id}:{configured_seed}:{ordinal}".encode()).digest()
    return int.from_bytes(raw[:4], "big") % 2_147_483_646 + 1


def manifest_entries(campaign_id, configured_seed, first_ordinal, source_inputs):
    """Pure ordered membership, shared by durable reservations and plan previews."""
    return [{"attempt_id": f"{campaign_id}:a{ordinal:08d}", "ordinal": ordinal,
             "job_position": position, "seed": attempt_seed(campaign_id, configured_seed, ordinal),
             "source_input": deepcopy(source),
             "native_rng_mapping": {"kind": "independent_requested_seed", "requested_seed": attempt_seed(campaign_id, configured_seed, ordinal),
                                    "native_mapping_verified": False}}
            for position, (ordinal, source) in enumerate(zip(range(first_ordinal, first_ordinal + len(source_inputs)), source_inputs))]


def _check_contract(state):
    """Detect changed configuration/membership and orphaned attempts on read."""
    if state.get("schema_version") != 2:
        return
    members = set()
    for job_id, job in state["execution_jobs"].items():
        if job.get("record_version") != 1 or job.get("mode") not in {"accounting_only", "single_attempt_adapter"}:
            raise LedgerConflict("Unsupported execution-job record version or adapter mode.")
        if not job["manifest"] or (job["mode"] == "single_attempt_adapter" and len(job["manifest"]) != 1):
            raise LedgerConflict("A production adapter job requires exactly one authorized attempt.")
        if job_id != job["execution_job_id"] or digest(job["configuration"]) != job["configuration_digest"]:
            raise LedgerConflict("Execution-job configuration digest mismatch.")
        if job["configuration"]["campaign"] != state["config"]:
            raise LedgerConflict("Campaign configuration changed after authorization.")
        if digest(job["manifest"]) != job["manifest_digest"]:
            raise LedgerConflict("Execution-job manifest digest mismatch.")
        for position, member in enumerate(job["manifest"]):
            attempt = state["attempts"][member["attempt_id"]]
            if member["attempt_id"] in members or member["job_position"] != position:
                raise LedgerConflict("Attempt membership must be unique and ordered.")
            members.add(member["attempt_id"])
            if attempt.get("execution_job_id") != job_id or any(attempt.get(k) != v for k, v in member.items()):
                raise LedgerConflict("Immutable attempt membership or provenance changed.")
    if members != set(state["attempts"]):
        raise LedgerConflict("Schema-v2 attempts require execution-job membership.")


def check_contract(state):
    try:
        _check_contract(state)
    except LedgerConflict:
        raise
    except (KeyError, TypeError, IndexError, ValueError) as exc:
        raise LedgerConflict("Malformed execution-job ledger contract.") from exc


def _execution(job, execution_id=None):
    key = execution_id or job["current_execution_id"]
    for execution in job["executions"]:
        if execution["execution_id"] == key:
            return execution
    raise LedgerConflict("Unknown execution identity for this job.")


def _new_execution(job, retry_key):
    identity = f"{job['execution_job_id']}:e{len(job['executions']) + 1:04d}"
    execution = {"execution_id": identity, "retry_key": retry_key,
                 "submission_key": digest({"execution_id": identity, "manifest_digest": job["manifest_digest"]}),
                 "submission_state": "reserved", "launch_state": "not_exposed", "contact_state": "unknown",
                 "scheduler_identifiers": {}, "worker_identity": {}, "resource_allocation": {},
                 "lease_receipts": {}, "created_at": utc_now()}
    job["executions"].append(execution)
    job["current_execution_id"] = identity
    return execution


def reserve(state, request_key, sources, native_configuration, resources, *, mode, batch_id=None):
    if state["schema_version"] != 2:
        raise ValueError("Execution-job operations require schema v2; historical ledgers are not migrated.")
    if not isinstance(request_key, str) or not request_key.strip() or not isinstance(sources, list) or not sources or any(not isinstance(s, dict) or not s for s in sources):
        raise ValueError("A reservation requires a request key and nonempty source-input objects.")
    if not isinstance(native_configuration, dict) or not isinstance(resources, dict):
        raise ValueError("Native configuration and resources must be objects.")
    declaration = {"sources": sources, "native_configuration": native_configuration, "resources": resources, "mode": mode}
    reservation_digest = digest(declaration)
    for job in state["execution_jobs"].values():
        if job["request_key"] == request_key:
            if job["reservation_digest"] != reservation_digest:
                raise LedgerConflict("Reservation request key conflicts with its durable content.")
            return deepcopy(job)
    if state["status"] != "running":
        raise ValueError("A stopped campaign cannot authorize a new job.")
    if len(state["attempts"]) + len(sources) > state["config"]["max_attempts"]:
        raise ValueError("The entire manifest must fit the remaining scientific attempt budget.")
    first = state["next_ordinal"]
    job_id = f"{state['campaign_id']}:j{first:08d}"
    manifest = manifest_entries(state["campaign_id"], state["config"]["seed"], first, sources)
    configuration = {"campaign": deepcopy(state["config"]), "native": deepcopy(native_configuration), "resource_requirements": deepcopy(resources)}
    job = {"record_version": 1, "execution_job_id": job_id, "request_key": request_key,
           "reservation_digest": reservation_digest, "mode": mode, "configuration": configuration,
           "configuration_digest": digest(configuration), "manifest": manifest, "manifest_digest": digest(manifest),
           "authorized_at": utc_now(), "executions": [], "output_receipts": {}}
    _new_execution(job, "initial")
    state["execution_jobs"][job_id] = job
    for entry in manifest:
        state["attempts"][entry["attempt_id"]] = {**deepcopy(entry), "execution_job_id": job_id,
            "batch_id": batch_id, "status": "authorized", "native_status": "authorized", "native_started": None,
            "authorized_at": job["authorized_at"], "started_at": None, "completed_at": None,
            "child_run_dir": None, "trajectory_index": None, "retry_count": 0, "receipts": []}
    state["next_ordinal"] += len(sources)
    state["events"].append({"type": "execution_job_authorized", "execution_job_id": job_id,
                            "manifest_digest": job["manifest_digest"], "at": utc_now()})
    return deepcopy(job)


def _identity_update(execution, field, values):
    for key, value in values.items():
        if value is None or value == "":
            continue
        old = execution[field].get(key)
        if old is not None and old != value:
            raise LedgerConflict(f"Conflicting {field}: {key}.")
        execution[field][key] = deepcopy(value)


def _prepare(job):
    execution = _execution(job)
    if execution["submission_state"] == "reserved":
        execution.update(submission_state="submission_pending", launch_state="may_have_started", submission_intent_at=utc_now())
    return deepcopy(execution)


def _retry(state, job, retry_key):
    for execution in job["executions"]:
        if execution["retry_key"] == retry_key:
            return deepcopy(execution)
    previous = _execution(job)
    if previous["launch_state"] not in {"not_exposed", "not_started_confirmed"} or job["output_receipts"]:
        raise ValueError("Ambiguous or post-launch stochastic work cannot be retried.")
    if previous["submission_state"] not in {"failed", "interrupted", "cancelled"}:
        raise ValueError("Execution retry requires a settled unsuccessful execution.")
    for entry in job["manifest"]:
        attempt = state["attempts"][entry["attempt_id"]]
        if attempt["native_started"] is True:
            raise ValueError("Started scientific attempts cannot be replayed as an execution retry.")
        attempt.update(native_status="authorized", native_started=None, status="authorized", started_at=None, completed_at=None)
    return deepcopy(_new_execution(job, retry_key))


def _result(state, job, execution_id, receipt_id, attempt_id, manifest_digest, outcome,
            native_state, native_output_ids, candidates, evidence):
    if not isinstance(receipt_id, str) or not receipt_id.strip() or not evidence or not isinstance(evidence, dict):
        raise ValueError("A verified receipt needs an identity and evidence object.")
    if manifest_digest != job["manifest_digest"]:
        raise LedgerConflict("Receipt manifest digest does not match authorization.")
    execution = _execution(job, execution_id)
    attempt = state["attempts"][attempt_id]
    if attempt["execution_job_id"] != job["execution_job_id"]:
        raise LedgerConflict("Receipt attempt is not a member of this execution job.")
    if execution.get("terminal_receipt", {}).get("not_started"):
        raise LedgerConflict("Native result conflicts with evidence that the execution never started.")
    if outcome not in {"completed", "failed"} or not isinstance(native_state, dict):
        raise ValueError("Native receipts describe verified completed or failed attempts.")
    if not isinstance(native_output_ids, list) or any(not isinstance(v, str) or not v for v in native_output_ids) or len(set(native_output_ids)) != len(native_output_ids):
        raise ValueError("Native output identifiers must be unique nonempty strings.")
    content = {"execution_job_id": job["execution_job_id"], "execution_id": execution_id, "attempt_id": attempt_id,
               "manifest_digest": manifest_digest, "outcome": outcome, "native_state": native_state,
               "native_output_ids": native_output_ids, "candidates": candidates, "evidence": evidence}
    content_digest = digest(content)
    if receipt_id in state["receipts"]:
        if state["receipts"][receipt_id].get("content_digest") != content_digest:
            raise LedgerConflict("Receipt identity replayed with conflicting content.")
        return False
    if attempt.get("result_digest"):
        if attempt["result_digest"] != content_digest:
            raise LedgerConflict("Attempt already has a different verified result.")
        return False
    if execution_id != job["current_execution_id"]:
        raise LedgerConflict("Result from a superseded execution requires manual reconciliation.")
    if attempt["native_status"] in {"cancelled", "not_started"}:
        raise LedgerConflict("Native result conflicts with evidence that the attempt never started.")
    if not isinstance(candidates, list) or (outcome == "failed" and candidates):
        raise ValueError("Only completed attempt receipts may contain candidates.")
    ids = [c.get("candidate_id") for c in candidates]
    if any(not isinstance(v, str) or not v for v in ids) or len(set(ids)) != len(ids) or any(v in state["candidates"] for v in ids):
        raise LedgerConflict("Candidate identities must be unique across scientific attempts.")
    from .generator_capabilities import generator_capabilities
    capability = generator_capabilities(state["config"]["generator_id"])
    for candidate in candidates:
        flag = candidate.get("metrics", {}).get("native_pass_filters")
        acceptance = "unavailable"
        if capability.native_acceptance_available and native_state.get("authoritative") is not False:
            acceptance = "accepted" if flag is True else "rejected" if flag is False else "unavailable"
        state["candidates"][candidate["candidate_id"]] = {
            "candidate_id": candidate["candidate_id"], "attempt_id": attempt_id,
            "execution_job_id": job["execution_job_id"], "execution_id": execution_id,
            "native_output_ids": deepcopy(native_output_ids), "receipt_id": receipt_id,
            "native_acceptance": acceptance, "validated_acceptance": "pending", "diversity_acceptance": "pending",
            "selection_state": "pending", "candidate": deepcopy(candidate)}
    now = utc_now()
    attempt.update(status=outcome, native_status=outcome, native_started=True, native_state=deepcopy(native_state),
                   native_output_ids=deepcopy(native_output_ids), candidate_ids=ids, result_digest=content_digest,
                   result_execution_id=execution_id, completed_at=now, result_evidence=deepcopy(evidence))
    attempt.setdefault("native_start_evidence", deepcopy(evidence))
    execution["launch_state"] = "may_have_started"
    if execution["submission_state"] == "reserved":
        execution["submission_state"] = "running"  # Native evidence can arrive before a scheduler acknowledgement.
    attempt["receipts"].append(receipt_id)
    receipt = {**deepcopy(content), "content_digest": content_digest, "at": now}
    state["receipts"][receipt_id] = receipt
    job["output_receipts"][receipt_id] = receipt
    state["events"].append({"type": "attempt_result_verified", "attempt_id": attempt_id, "receipt_id": receipt_id, "at": now})
    return True


def observe_single_attempt(state, attempt, status, details):
    """Bridge existing worker states without claiming native launch from a PID."""
    job = state["execution_jobs"][attempt["execution_job_id"]]
    if job["mode"] != "single_attempt_adapter" or len(job["manifest"]) != 1:
        raise ValueError("Accounting-only execution jobs cannot enter campaign dispatch.")
    execution = _execution(job)
    if status in {"queued", "running"}:
        _prepare(job)
        execution["submission_state"] = status
    if details.get("child_run_dir"):
        _identity_update(execution, "scheduler_identifiers", {"child_run_dir": details["child_run_dir"]})
    _identity_update(execution, "scheduler_identifiers", {k: details[k] for k in ("cluster_job_id", "remote_submission_key") if k in details})
    _identity_update(execution, "worker_identity", {"worker_node_id": details.get("worker_node_id")})
    if status in {"failed", "interrupted", "cancelled", "not_started"}:
        execution["submission_state"] = "cancelled" if status == "not_started" else status
        # A legacy worker cancellation/failure is not proof that native work never ran.
        if execution["launch_state"] == "not_exposed":
            attempt.update(native_status="cancelled" if status == "cancelled" else "not_started", native_started=False)
        elif attempt.get("native_status") not in {"completed", "failed"}:
            attempt["native_status"] = "unresolved"  # Keep any verified native-start evidence.


def record_lease(job, execution_id, lease_id, gpu_seconds):
    if not isinstance(lease_id, str) or not lease_id.strip() or type(gpu_seconds) not in {int, float} or not math.isfinite(gpu_seconds) or gpu_seconds < 0:
        raise ValueError("A job lease needs an identity and finite nonnegative GPU seconds.")
    receipts = _execution(job, execution_id)["lease_receipts"]
    if lease_id in receipts:
        if receipts[lease_id] != float(gpu_seconds):
            raise LedgerConflict("Lease receipt conflicts with its recorded duration.")
        return False
    receipts[lease_id] = float(gpu_seconds)
    return True


class ExecutionJobLedgerMixin:
    """Mutations are supplied the parent ledger's locked atomic transaction."""

    def reserve_execution_job(self, request_key, source_inputs, *, native_configuration=None, resource_requirements=None):
        return self._job_transaction(lambda state: reserve(state, request_key, source_inputs, {} if native_configuration is None else native_configuration,
                                                         {} if resource_requirements is None else resource_requirements, mode="accounting_only"))

    def prepare_execution_submission(self, job_id):
        """Persist intent; replay returns the same key, never submits to a scheduler."""
        return self._job_transaction(lambda state: _prepare(state["execution_jobs"][job_id]))

    def observe_execution(self, job_id, execution_id, *, submission_key, status, scheduler_identifiers=None,
                          worker_identity=None, resource_allocation=None, contact_state="connected"):
        if status not in {"queued", "running"} or contact_state not in {"connected", "disconnected", "unknown"}:
            raise ValueError("Invalid submission observation.")
        def change(state):
            job = state["execution_jobs"][job_id]
            execution = _execution(job, execution_id)
            if execution_id != job["current_execution_id"] or execution["submission_key"] != submission_key:
                raise LedgerConflict("Submission identity does not match the current execution.")
            if execution["submission_state"] not in {"submission_pending", "queued", "running"}:
                raise ValueError("Observation requires a durable submission intent and active execution.")
            if execution["submission_state"] == "running" and status == "queued":
                raise ValueError("An execution cannot move backwards from running to queued.")
            _identity_update(execution, "scheduler_identifiers", scheduler_identifiers or {})
            _identity_update(execution, "worker_identity", worker_identity or {})
            _identity_update(execution, "resource_allocation", resource_allocation or {})
            execution.update(submission_state=status, contact_state=contact_state, launch_state="may_have_started")
            return deepcopy(execution)
        return self._job_transaction(change)

    def record_native_start(self, job_id, execution_id, attempt_id, *, manifest_digest, evidence):
        def change(state):
            job = state["execution_jobs"][job_id]
            execution = _execution(job, execution_id)
            attempt = state["attempts"][attempt_id]
            if (execution_id != job["current_execution_id"] or manifest_digest != job["manifest_digest"]
                    or attempt["execution_job_id"] != job_id):
                raise LedgerConflict("Native start is not in the current authorized manifest.")
            if not isinstance(evidence, dict) or not evidence:
                raise ValueError("Native start requires explicit native evidence.")
            if attempt["native_status"] not in {"authorized", "running"}:
                raise ValueError("Native start conflicts with the settled attempt state.")
            if attempt.get("native_start_evidence") and attempt["native_start_evidence"] != evidence:
                raise LedgerConflict("Conflicting native start evidence.")
            execution["launch_state"] = "may_have_started"
            if execution["submission_state"] == "reserved":
                execution["submission_state"] = "running"
            attempt.update(native_status="running", native_started=True, native_start_evidence=deepcopy(evidence))
        return self._job_transaction(change)

    def record_attempt_result(self, job_id, execution_id, receipt_id, attempt_id, *, manifest_digest,
                              native_state, native_output_ids, candidates, evidence, outcome="completed"):
        return self._job_transaction(lambda state: _result(state, state["execution_jobs"][job_id], execution_id, receipt_id,
            attempt_id, manifest_digest, outcome, native_state, native_output_ids, candidates, evidence))

    def finish_execution(self, job_id, execution_id, *, status, evidence, not_started=False):
        if status not in {"completed", "failed", "interrupted", "cancelled"} or not isinstance(evidence, dict) or not evidence or type(not_started) is not bool:
            raise ValueError("Execution completion requires a terminal state and evidence.")
        def change(state):
            job = state["execution_jobs"][job_id]
            execution = _execution(job, execution_id)
            terminal = {"status": status, "evidence": evidence, "not_started": not_started}
            if execution.get("terminal_receipt"):
                if execution["terminal_receipt"] != terminal:
                    raise LedgerConflict("Execution completion conflicts with its terminal receipt.")
                return False
            if execution_id != job["current_execution_id"]:
                raise LedgerConflict("Cannot finalize a superseded execution.")
            attempts = [state["attempts"][entry["attempt_id"]] for entry in job["manifest"]]
            if execution["launch_state"] == "not_started_confirmed" and not not_started:
                raise LedgerConflict("Completion cannot erase confirmed no-launch evidence.")
            if not_started and any(a["native_started"] is True for a in attempts):
                raise LedgerConflict("No-launch evidence conflicts with verified native work.")
            execution.update(submission_state=status, terminal_receipt=deepcopy(terminal), completed_at=utc_now())
            if not_started:
                execution["launch_state"] = "not_started_confirmed"
            for attempt in attempts:
                if attempt["native_status"] in {"completed", "failed"}:
                    continue
                if not_started:
                    attempt.update(native_status="cancelled" if status == "cancelled" else "not_started", native_started=False,
                                   status="cancelled" if status == "cancelled" else "not_started")
                else:
                    attempt.update(native_status="unresolved", status="interrupted")
            return True
        return self._job_transaction(change)

    def cancel_execution_job(self, job_id):
        def change(state):
            job = state["execution_jobs"][job_id]
            execution = _execution(job)
            if execution["submission_state"] == "cancelled":
                return deepcopy(execution)
            if execution["submission_state"] == "reserved" and execution["launch_state"] == "not_exposed":
                execution.update(submission_state="cancelled", launch_state="not_started_confirmed")
                for member in job["manifest"]:
                    state["attempts"][member["attempt_id"]].update(status="cancelled", native_status="cancelled", native_started=False)
            else:
                execution["cancel_requested"] = True
            return deepcopy(execution)
        return self._job_transaction(change)

    def retry_execution_job(self, job_id, *, retry_key):
        if not isinstance(retry_key, str) or not retry_key.strip() or retry_key == "initial":
            raise ValueError("Execution retries require a new durable retry key.")
        return self._job_transaction(lambda state: _retry(state, state["execution_jobs"][job_id], retry_key))

    def record_job_compute(self, job_id, execution_id, lease_id, gpu_seconds):
        def change(state):
            for other_job in state["execution_jobs"].values():
                for execution in other_job["executions"]:
                    if lease_id in execution["lease_receipts"] and (other_job["execution_job_id"], execution["execution_id"]) != (job_id, execution_id):
                        raise LedgerConflict("A scheduler lease cannot belong to multiple executions.")
            return record_lease(state["execution_jobs"][job_id], execution_id, lease_id, gpu_seconds)
        return self._job_transaction(change)

    def execution_job_summary(self):
        state = self.snapshot()
        if state["schema_version"] == 1:
            return {"schema_version": 1, "execution_jobs": None, "interpretation": "Historical batch_id is dispatch grouping; native job membership is unavailable."}
        rows = []
        for job in state["execution_jobs"].values():
            counts = {status: sum(state["attempts"][m["attempt_id"]]["native_status"] == status for m in job["manifest"])
                      for status in ("authorized", "running", "completed", "failed", "not_started", "cancelled", "unresolved")}
            rows.append({"execution_job_id": job["execution_job_id"], "manifest_digest": job["manifest_digest"],
                         "attempt_counts": counts, "execution_retries": len(job["executions"]) - 1,
                         "submission_state": _execution(job)["submission_state"]})
        return {"schema_version": 2, "execution_jobs": rows}
