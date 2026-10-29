"""Read-only campaign evidence projections for the P6 results interface."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any


_VALIDATION_METRICS = ("iptm", "ptm", "plddt_mean")


def _read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    return value if isinstance(value, dict) else {}


def _resolved_path(value: object, *, base: Path) -> Path | None:
    if not value:
        return None
    from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path

    try:
        resolved = resolve_stored_path(str(value), run_dir=base, must_exist=True)
    except (OSError, ValueError):
        return None
    return resolved.resolve() if resolved is not None else None


def _candidate_index(ledger_root: Path) -> tuple[dict[str, int], bool]:
    path = ledger_root / "artifacts" / "normalized_candidates" / "candidates.jsonl"
    if not path.is_file():
        return {}, False
    positions: dict[str, int] = {}
    try:
        with path.open(encoding="utf-8") as handle:
            for index, line in enumerate(handle, start=1):
                try:
                    item = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if isinstance(item, dict) and item.get("candidate_id"):
                    positions.setdefault(str(item["candidate_id"]), index)
    except OSError:
        return {}, False
    return positions, True


def _number(value: object, metric: str) -> float | None:
    if isinstance(value, bool):
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(result):
        return None
    if metric in {"iptm", "ptm"} and not 0 <= result <= 1:
        return None
    if metric == "plddt_mean" and not 0 <= result <= 100:
        return None
    return result


def _metric_projection(stage: dict[str, Any], required: list[str]) -> dict[str, Any]:
    evidence = stage.get("evidence_candidate")
    metrics = evidence.get("metrics") if isinstance(evidence, dict) else None
    metrics = metrics if isinstance(metrics, dict) else {}
    required = required or list(_VALIDATION_METRICS)
    values = {name: _number(metrics.get(name), name) for name in required}
    valid = bool(values) and all(value is not None for value in values.values())
    return {
        "measurement_status": "Valid measurements" if valid else "Measurements unavailable or incomplete",
        "metrics": {name: value for name, value in values.items() if value is not None},
        "required_metrics": required,
    }


def _resource_allocation(run_dir: object, base: Path) -> dict[str, Any]:
    path = _resolved_path(run_dir, base=base)
    if path is None:
        return {}
    metadata = _read_json(path / "metadata.json")
    allocation = metadata.get("resource_allocation")
    return allocation if isinstance(allocation, dict) else {}


def _gpu_devices(allocation: dict[str, Any]) -> str:
    values = allocation.get("gpu_devices") or allocation.get("gpu_ids") or []
    if not values:
        device = str(allocation.get("gpu_device") or "").strip()
        return device if device and device.casefold() not in {"cpu", "none", "auto"} else "—"
    return ", ".join(str(value) for value in values)


def _acceptance_label(status: object, reason: str) -> tuple[str, bool]:
    normalized = str(status or "unknown").strip().casefold()
    text = reason.casefold()
    lineage_problem = any(
        marker in text
        for marker in ("same-sequence descendant", "normalized candidate", "candidate index", "lineage", "index metadata")
    )
    if lineage_problem:
        return "Unresolved · candidate index/lineage metadata", True
    if normalized == "accepted":
        return "Accepted", False
    if normalized == "rejected":
        return "Rejected by campaign decision", False
    if normalized in {"failed", "error"}:
        return "Validation failed · no scientific rejection recorded", False
    if normalized == "unavailable":
        return "Unavailable at this acceptance level", False
    if normalized in {"pending", "running", "queued"}:
        return normalized.title(), False
    return normalized.replace("_", " ").title(), False


def _native_batch_size(
    run_dir: Path,
    input_payload: dict[str, Any],
    candidate_records: dict[str, dict[str, Any]],
    ledger_root: Path,
) -> tuple[object, str]:
    params = input_payload.get("params") if isinstance(input_payload.get("params"), dict) else {}
    value = params.get("batch_size")
    if value is not None:
        return value, run_dir.name

    # Validation-only records live separately from their retained generation inputs.
    for record in candidate_records.values():
        candidate = record.get("candidate") if isinstance(record.get("candidate"), dict) else {}
        raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
        source_run = _resolved_path(raw.get("source_run_dir"), base=ledger_root)
        source_input = _read_json(source_run / "input.json") if source_run else {}
        source_params = source_input.get("params") if isinstance(source_input.get("params"), dict) else {}
        value = source_params.get("batch_size")
        if value is not None:
            return value, source_run.name

    manifest = _read_json(ledger_root / "manifest.json")
    source_candidates = _resolved_path(manifest.get("source_candidates_path"), base=ledger_root)
    if source_candidates is not None:
        source_root = source_candidates.parent.parent.parent
        source_input = _read_json(source_root / "input.json")
        source_params = source_input.get("params") if isinstance(source_input.get("params"), dict) else {}
        value = source_params.get("batch_size")
        if value is not None:
            return value, source_root.name
    return None, ""


def campaign_context_for_job(
    run_dir: str | Path,
    metadata: dict[str, Any] | None = None,
    input_payload: dict[str, Any] | None = None,
) -> dict[str, Any] | None:
    """Return a compact campaign projection for a campaign or its validation child."""
    current_dir = Path(run_dir).expanduser().resolve()
    metadata = metadata if isinstance(metadata, dict) else {}
    input_payload = input_payload if isinstance(input_payload, dict) else _read_json(current_dir / "input.json")
    ledger_root = current_dir
    ledger_path = ledger_root / "artifacts" / "distributed_campaign" / "ledger.json"

    if not ledger_path.is_file():
        parent = _resolved_path(metadata.get("parent_run_dir"), base=current_dir)
        if parent is None:
            return None
        parent_ledger = parent / "artifacts" / "distributed_campaign" / "ledger.json"
        if not parent_ledger.is_file():
            return None
        expected_id = str(metadata.get("distributed_campaign_id") or "")
        parent_ledger_value = _read_json(parent_ledger)
        if expected_id and str(parent_ledger_value.get("campaign_id") or "") != expected_id:
            return None
        ledger_root, ledger_path = parent, parent_ledger

    ledger = _read_json(ledger_path)
    campaign_id = str(ledger.get("campaign_id") or "")
    if not campaign_id:
        return None
    config = ledger.get("config") if isinstance(ledger.get("config"), dict) else {}
    attempts_raw = ledger.get("attempts") if isinstance(ledger.get("attempts"), dict) else {}
    candidates_raw = ledger.get("candidates") if isinstance(ledger.get("candidates"), dict) else {}
    attempt_statuses = [str(value.get("status") or "unknown").casefold() for value in attempts_raw.values() if isinstance(value, dict)]
    attempts_completed = sum(status == "completed" for status in attempt_statuses)
    active_attempts = sum(status in {"authorized", "queued", "running", "retrying", "submitted"} for status in attempt_statuses)
    accepted_target = config.get("accepted_target")
    attempt_budget = config.get("max_attempts")
    index_by_candidate, index_file_present = _candidate_index(ledger_root)
    acceptance_level = str(config.get("acceptance_level") or "validated")
    required_by_stage = {
        str(index): list(stage.get("required_metrics") or _VALIDATION_METRICS)
        for index, stage in enumerate(config.get("validation_stages") or [])
        if isinstance(stage, dict)
    }

    # Candidate records are compactly projected; raw structures and full metrics stay in Results.
    source_campaign_ids: set[str] = set()
    candidates: list[dict[str, Any]] = []
    stage_rows: list[dict[str, Any]] = []
    accepted_count = rejected_count = unresolved_count = 0
    valid_measurement_count = 0

    for candidate_id, raw_record in candidates_raw.items():
        if not isinstance(raw_record, dict):
            continue
        record_candidate = raw_record.get("candidate") if isinstance(raw_record.get("candidate"), dict) else {}
        raw_metadata = record_candidate.get("raw_metadata") if isinstance(record_candidate.get("raw_metadata"), dict) else {}
        distributed = raw_metadata.get("distributed_campaign") if isinstance(raw_metadata.get("distributed_campaign"), dict) else {}
        source_id = str(distributed.get("campaign_id") or "")
        if source_id and source_id != campaign_id:
            source_campaign_ids.add(source_id)

        reason = str(raw_record.get("validation_decision_reason") or raw_record.get("decision_reason") or "")
        acceptance = raw_record.get(f"{acceptance_level}_acceptance")
        acceptance_label, lineage_error = _acceptance_label(acceptance, reason)
        if acceptance_label == "Accepted":
            accepted_count += 1
        elif acceptance_label == "Rejected by campaign decision":
            rejected_count += 1

        stages_raw = raw_record.get("validation_stages") if isinstance(raw_record.get("validation_stages"), dict) else {}
        projected_stages: list[dict[str, Any]] = []
        for stage_index, stage in sorted(stages_raw.items(), key=lambda pair: str(pair[0])):
            if not isinstance(stage, dict):
                continue
            required = required_by_stage.get(str(stage_index), list(_VALIDATION_METRICS))
            measurement = _metric_projection(stage, required)
            if measurement["measurement_status"] == "Valid measurements":
                valid_measurement_count += 1
            stage_run = _resolved_path(stage.get("run_dir"), base=ledger_root)
            allocation = _resource_allocation(stage.get("run_dir"), ledger_root)
            stage_row = {
                "Candidate ID": str(candidate_id),
                "Stage": str(stage_index),
                "Status": str(stage.get("status") or "unknown").replace("_", " ").title(),
                "Prediction measurements": measurement["measurement_status"],
                **{name: value for name, value in measurement["metrics"].items()},
                "Validation run": stage_run.name if stage_run else str(stage.get("run_dir") or "—"),
                "GPU device(s)": _gpu_devices(allocation),
                "Detail": str(stage.get("error") or ""),
            }
            stage_rows.append(stage_row)
            projected_stages.append({
                "stage": str(stage_index),
                "status": str(stage.get("status") or "unknown"),
                "measurement_status": measurement["measurement_status"],
                "metrics": measurement["metrics"],
                "run_id": stage_run.name if stage_run else "",
                "gpu_devices": _gpu_devices(allocation),
                "error": str(stage.get("error") or ""),
            })

        attempt = attempts_raw.get(str(raw_record.get("attempt_id") or ""), {})
        attempt = attempt if isinstance(attempt, dict) else {}
        native_state = attempt.get("native_state") if isinstance(attempt.get("native_state"), dict) else {}
        allocation = attempt.get("resource_allocation")
        allocation = allocation if isinstance(allocation, dict) else _resource_allocation(attempt.get("child_run_dir"), ledger_root)
        worker = str(attempt.get("worker") or native_state.get("source_worker") or distributed.get("worker") or "—")
        worker_run_id = str(
            attempt.get("worker_run_id") or native_state.get("source_worker_run_id") or distributed.get("worker_run_id") or "—"
        )
        position = index_by_candidate.get(str(candidate_id))
        if position is None:
            index_label = "Missing index file" if not index_file_present else "Not found in index"
        else:
            index_label = f"Row {position}"
        index_missing = not index_file_present or position is None
        metadata_issue = lineage_error or index_missing
        if metadata_issue:
            unresolved_count += 1

        candidates.append({
            "candidate_id": str(candidate_id),
            "candidate_index": position,
            "candidate_index_status": index_label,
            "acceptance_level": acceptance_level,
            "acceptance_status": str(acceptance or "unknown"),
            "acceptance_label": acceptance_label,
            "lineage_error": lineage_error,
            "index_missing": index_missing,
            "metadata_issue": metadata_issue,
            "decision_reason": reason,
            "source_campaign_id": source_id,
            "source_candidate_id": str(candidate_id) if source_id else "",
            "worker": worker,
            "worker_run_id": worker_run_id,
            "attempt_id": str(raw_record.get("attempt_id") or ""),
            "gpu_devices": _gpu_devices(allocation),
            "validation_stages": projected_stages,
        })

    attempts: list[dict[str, Any]] = []
    for attempt_id, attempt in attempts_raw.items():
        if not isinstance(attempt, dict):
            continue
        native_state = attempt.get("native_state") if isinstance(attempt.get("native_state"), dict) else {}
        allocation = attempt.get("resource_allocation")
        allocation = allocation if isinstance(allocation, dict) else _resource_allocation(attempt.get("child_run_dir"), ledger_root)
        matching = next((item for item in candidates if item["attempt_id"] == str(attempt_id)), {})
        attempts.append({
            "Attempt": str(attempt.get("ordinal") or attempt_id),
            "Batch": str(attempt.get("batch_id") or "—"),
            "Status": str(attempt.get("status") or "unknown").replace("_", " ").title(),
            "Candidate(s)": ", ".join(str(value) for value in attempt.get("candidate_ids", [])) or "—",
            "Worker": str(attempt.get("worker") or native_state.get("source_worker") or matching.get("worker") or "—"),
            "Worker run": str(attempt.get("worker_run_id") or native_state.get("source_worker_run_id") or matching.get("worker_run_id") or "—"),
            "GPU device(s)": _gpu_devices(allocation),
            "Scheduler allocation": str(allocation.get("scheduler_job_id") or "—"),
            "Attempt unit": str(attempt.get("attempt_unit") or native_state.get("attempt_unit") or "—"),
        })

    validation_stage_configs = config.get("validation_stages") if isinstance(config.get("validation_stages"), list) else []
    for attempt in attempts_raw.values():
        native_state = attempt.get("native_state") if isinstance(attempt, dict) and isinstance(attempt.get("native_state"), dict) else {}
        source_id = str(native_state.get("source_campaign_id") or "")
        if source_id and source_id != campaign_id:
            source_campaign_ids.add(source_id)
    validation_max = config.get("max_concurrent_validations")
    dispatch_batch = config.get("batch_size")
    concurrent_batches = config.get("max_concurrent_batches")
    dispatch_capacity = None
    if isinstance(dispatch_batch, int) and isinstance(concurrent_batches, int):
        dispatch_capacity = dispatch_batch * concurrent_batches

    native_batch_size, native_source_run_id = _native_batch_size(current_dir, input_payload, candidates_raw, ledger_root)
    objective_type = str(config.get("objective_type") or "accepted")
    if objective_type == "generated":
        objective_count = len(candidates_raw)
    elif config.get("unique_accepted"):
        field = f"{acceptance_level}_acceptance"
        unique_sequences = {
            "".join(str((record.get("candidate") or {}).get("binder_sequence") or "").split()).upper()
            for record in candidates_raw.values() if isinstance(record, dict) and record.get(field) == "accepted"
        }
        unique_sequences.discard("")
        objective_count = len(unique_sequences)
    else:
        objective_count = accepted_count
    validation_only = bool(source_campaign_ids) or any(
        "retained" in str(
            (attempt.get("native_state") if isinstance(attempt.get("native_state"), dict) else {}).get("attempt_unit") or ""
        ).casefold()
        for attempt in attempts_raw.values() if isinstance(attempt, dict)
    )

    return {
        "campaign_id": campaign_id,
        "campaign_kind": "Validation-only follow-up" if validation_only else "Original generation campaign",
        "status": str(ledger.get("status") or "unknown"),
        "termination_reason": str(ledger.get("termination_reason") or ""),
        "accepted_target": accepted_target,
        "objective_type": objective_type,
        "objective_count": objective_count,
        "attempt_budget": attempt_budget,
        "attempts_authorized": len(attempts_raw),
        "attempts_completed": attempts_completed,
        "active_attempts": active_attempts,
        "source_campaign_ids": sorted(source_campaign_ids),
        "acceptance_level": acceptance_level,
        "index_file_present": index_file_present,
        "index_candidate_count": len(index_by_candidate),
        "controller_batch_size": dispatch_batch,
        "max_concurrent_batches": concurrent_batches,
        "controller_dispatch_capacity": dispatch_capacity,
        "max_concurrent_validations": validation_max,
        "native_inference_batch_size": native_batch_size,
        "native_inference_source_run_id": native_source_run_id,
        "validation_stage_count": len(validation_stage_configs),
        "candidate_count": len(candidates),
        "accepted_count": accepted_count,
        "rejected_count": rejected_count,
        "unresolved_count": unresolved_count,
        "valid_measurement_count": valid_measurement_count,
        "candidates": candidates,
        "attempts": attempts,
        "validation_stages": stage_rows,
    }


def render_campaign_context(st: Any, context: dict[str, Any] | None) -> None:
    """Render the campaign projection in local or remotely inspected Results."""
    if not context:
        return
    with st.expander("Distributed campaign evidence", expanded=True):
        kind = str(context.get("campaign_kind") or "Campaign")
        st.markdown(f"**{kind}** · Campaign `{context.get('campaign_id') or 'unknown'}` · Status: {str(context.get('status') or 'unknown').replace('_', ' ').title()}")
        if kind == "Validation-only follow-up":
            sources = ", ".join(f"`{value}`" for value in context.get("source_campaign_ids", [])) or "the retained source campaign"
            st.info(f"This run validated retained candidates from {sources}; it did not generate new candidates. Its measurements are separate from the original campaign evidence, which remains unchanged.")

        dispatch = context.get("controller_dispatch_capacity")
        batch = context.get("controller_batch_size")
        batches = context.get("max_concurrent_batches")
        dispatch_text = f"{dispatch} max in flight ({batch} × {batches})" if dispatch is not None else "Not recorded"
        generation_text = "No generation in this run" if kind == "Validation-only follow-up" else dispatch_text
        native_batch = context.get("native_inference_batch_size")
        native_text = str(native_batch) if native_batch is not None else "Not recorded"
        concurrency_columns = st.columns(3)
        concurrency_columns[0].metric("Controller dispatch batch size", str(batch) if batch is not None else "Not recorded")
        concurrency_columns[1].metric("Generation concurrency", generation_text)
        concurrency_columns[2].metric("Downstream validation concurrency", str(context.get("max_concurrent_validations") or "Not recorded"))
        st.caption(f"Controller attempt dispatch capacity: {dispatch_text}. This ceiling is distinct from downstream validation concurrency.")
        st.caption(f"Native generator inference batch_size: {native_text} (source run {context.get('native_inference_source_run_id') or 'not recorded'}). Controller dispatch size and native inference batch size are separate settings.")

        accepted_target = context.get("accepted_target")
        attempt_budget = context.get("attempt_budget")
        if isinstance(accepted_target, int) and accepted_target > 0:
            accepted = int(context.get("objective_count", context.get("accepted_count")) or 0)
            attempt_text = ""
            if isinstance(attempt_budget, int) and attempt_budget > 0:
                attempt_text = (
                    f" · {context.get('attempts_completed', 0)}/{attempt_budget} attempts completed"
                    f" · {context.get('active_attempts', 0)} active"
                )
            objective_label = "generation" if context.get("objective_type") == "generated" else "acceptance"
            label = f"Global {objective_label} progress: {accepted}/{accepted_target}{attempt_text}"
            if kind == "Validation-only follow-up":
                label = "Validation-only · retained candidates · " + label
            st.caption(label)
            st.progress(min(1.0, accepted / accepted_target))

        summary_columns = st.columns(4)
        summary_columns[0].metric("Candidates", str(context.get("candidate_count", 0)))
        summary_columns[1].metric("Accepted at campaign level", str(context.get("accepted_count", 0)))
        summary_columns[2].metric("Valid prediction measurements", str(context.get("valid_measurement_count", 0)))
        summary_columns[3].metric("Unresolved index/lineage", str(context.get("unresolved_count", 0)))

        candidate_rows = []
        for candidate in context.get("candidates", []):
            stage_values = candidate.get("validation_stages") or []
            metrics = []
            measurement_status = "No current campaign measurements"
            for stage in stage_values:
                metric_text = ", ".join(f"{key}={value:g}" for key, value in stage.get("metrics", {}).items())
                if metric_text:
                    metrics.append(metric_text)
                if stage.get("measurement_status") == "Valid measurements":
                    measurement_status = "Valid measurements"
            interpretation = candidate.get("acceptance_label") or "Unknown"
            if candidate.get("lineage_error"):
                interpretation = "Unresolved · candidate-index/lineage metadata (not a scientific rejection)"
            source = candidate.get("source_campaign_id")
            if source:
                source = f"{source} · {candidate.get('source_candidate_id') or candidate.get('candidate_id')}"
            candidate_rows.append({
                "Candidate ID": candidate.get("candidate_id"),
                "Normalized candidate index": candidate.get("candidate_index_status"),
                "Prediction measurements": measurement_status + (" · " + "; ".join(metrics) if metrics else ""),
                f"Recorded global {context.get('acceptance_level', 'validated')} acceptance": str(candidate.get("acceptance_status") or "unknown").replace("_", " ").title(),
                "Acceptance interpretation": interpretation,
                "Source campaign · candidate": source or "Original campaign candidate",
                "Worker · run": f"{candidate.get('worker') or '—'} · {candidate.get('worker_run_id') or '—'}",
                "Generation GPU": candidate.get("gpu_devices") or "—",
        })
        if candidate_rows:
            st.dataframe(candidate_rows, hide_index=True, width="stretch")
        if not context.get("index_file_present"):
            st.warning("The normalized candidate index file is missing. This is a candidate identity/lineage metadata gap, not a scientific candidate rejection.")

        if context.get("validation_stages"):
            st.markdown("**Validation jobs and measurements**")
            st.dataframe(context["validation_stages"], hide_index=True, width="stretch")
        if context.get("attempts"):
            with st.expander("Generation attempts and worker provenance"):
                st.dataframe(context["attempts"], hide_index=True, width="stretch")
