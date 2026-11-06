"""Run-first landing page for the existing Complex Analysis workspace.

This module is a presentation and inventory layer over Studio's current Results,
Database, workflow registry and Evaluation records. It never imports a Results
run or submits work during page rendering.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import time
from functools import lru_cache
from pathlib import Path

import pandas as pd
import streamlit as st


_CATALOG_TTL_SECONDS = 20
_METHOD_LABELS = {
    "sequence_qc": "Sequence Composition",
    "biotite_interface": "Biotite interface geometry",
    "pep_patch": "PEP-Patch surface chemistry",
    "plip": "PLIP ligand interactions",
    "pandamap": "PandaMap ligand interactions",
    "posebusters": "PoseBusters pose checks",
    "sequence_liability": "Sequence and structure liabilities",
}
_UNCONFIRMED_JOB = "Completed job; Database Evaluation not confirmed"
_CANDIDATE_SUMMARY_CACHE: dict[str, dict] = {}


def _needs_execution(status: str) -> bool:
    return status in {"Not run", "Failed", "Stale or nonmatching"}


def _run_action_label(row: dict) -> str:
    states = {str(value).casefold() for value in row.get("analysis_job_statuses") or []}
    resumable = {"queued", "pending", "submitted", "scheduled", "starting", "preparing",
                 "running", "staging", "waiting", "failed"}
    if states & resumable:
        return "Continue analysis"
    return "Open dashboard" if int(row.get("analysis_job_count") or 0) else "Analyze"


def _run_catalog(user: dict, *, refresh: bool = False) -> list[dict]:
    key = f"complex-analysis-run-catalog-{user.get('id') or ''}"
    cached = st.session_state.get(key)
    if (not refresh and isinstance(cached, dict)
            and time.monotonic() - float(cached.get("loaded_at") or 0) < _CATALOG_TTL_SECONDS):
        return list(cached.get("runs") or [])
    from mn_molecular_design_studio.logical_runs import catalog

    runs = catalog(user, token=str(st.session_state.get("studio-session-token") or "development-mode"))
    st.session_state[key] = {"loaded_at": time.monotonic(), "runs": runs}
    return runs


def _result_sources(runs: list[dict]) -> list[dict]:
    from mn_molecular_design_studio.cluster import config

    local_node = config.identity()
    output_rows = []
    for logical_run in runs:
        activities = [job for job in logical_run.get("jobs") or []
                      if str(job.get("task_group") or "") == "analysis"]
        for output in logical_run.get("outputs") or []:
            if output.get("id") != "results":
                continue
            sources = []
            seen_sources = set()
            for source in output.get("sources") or []:
                if source.get("task_group") == "analysis":
                    continue
                token = (str(source.get("node_id") or ""), str(source.get("task_group") or ""),
                         str(source.get("run_id") or ""))
                if not token[1] or not token[2] or token in seen_sources:
                    continue
                seen_sources.add(token)
                source_row = dict(source)
                source_dir = Path(str(source_row.get("run_dir") or "")).expanduser()
                source_row["is_local"] = (
                    (not source_row.get("node_id") or str(source_row.get("node_id")) == local_node)
                    and source_dir.is_dir()
                )
                source_row["source_summary"] = _source_summary(source_row)
                linked_analyses = [analysis for analysis in activities
                    if str(analysis.get("upstream_task_group") or "") == token[1]
                    and str(analysis.get("upstream_run_id") or "") == token[2]
                    and str(analysis.get("node_id") or "") == token[0]]
                source_row["analysis_job_count"] = len(linked_analyses)
                source_row["analysis_job_statuses"] = [str(analysis.get("status") or "unknown").casefold()
                                                        for analysis in linked_analyses]
                source_row["last_analysis_activity"] = max(
                    (str(analysis.get("updated_at") or analysis.get("completed_at") or analysis.get("created_at") or "")
                     for analysis in linked_analyses), default="")
                sources.append(source_row)
            if not sources:
                continue

            summaries = [source.get("source_summary") or {} for source in sources]
            confirmed_structures = sum(int(item.get("structure_count") or 0) for item in summaries
                                       if item.get("structure_count") is not None)
            unknown_structure_candidates = sum(int(item.get("candidate_count") or 0) for item in summaries
                                               if item.get("structure_count") is None)
            # Complex Analysis lists only logical runs with at least one
            # confirmed candidate/pose structure. Sequence-only and empty jobs
            # belong in Results/Database, not this analysis queue.
            if confirmed_structures <= 0:
                continue

            target_systems = {}
            for summary in summaries:
                for target in summary.get("target_systems") or []:
                    target_systems[str(target.get("identity") or target.get("label") or "")] = target
            target_values = [value for key, value in target_systems.items() if key]
            if len(target_values) == 1:
                target_summary = target_values[0]
            elif len(target_values) > 1:
                target_summary = {"identity": "multiple", "label": f"Multiple targets ({len(target_values)})",
                                  "sequence_count": len(target_values)}
            else:
                target_summary = {"identity": "", "label": "Target not recorded"}

            analysis_jobs = sum(int(source.get("analysis_job_count") or 0) for source in sources)
            analysis_statuses = [status for source in sources for status in source.get("analysis_job_statuses", [])]
            engine_names = sorted({str(source.get("display_tool") or source.get("tool") or source.get("job_type") or "Unknown")
                                   for source in sources}, key=str.casefold)
            workstations = sorted({str(source.get("node_name") or "This workstation") for source in sources}, key=str.casefold)
            candidate_counts = [item.get("candidate_count") for item in summaries if item.get("candidate_count") is not None]
            acceptance_observed = sum(int(item.get("acceptance_observed") or 0) for item in summaries)
            accepted_count = sum(int(item.get("accepted_count") or 0) for item in summaries)
            rejected_count = sum(int(item.get("rejected_count") or 0) for item in summaries)
            unresolved_count = sum(int(item.get("unresolved_count") or 0) for item in summaries)
            not_evaluated_count = sum(int(item.get("not_evaluated_count") or 0) for item in summaries)
            root = dict(logical_run.get("root") or {})
            row = {
                **root,
                "logical_run_id": str(logical_run.get("id") or ""),
                "logical_run_name": str(logical_run.get("name") or "Run"),
                "logical_run_status": str(logical_run.get("status") or "unknown"),
                "logical_run_job_count": int(logical_run.get("job_count") or len(sources)),
                "output_id": str(output.get("id") or "results"),
                "output_name": str(output.get("name") or "Results"),
                "run_sources": sources,
                "source_job_count": len(sources),
                "structure_count": confirmed_structures,
                "unknown_structure_candidates": unknown_structure_candidates,
                "candidate_count": sum(map(int, candidate_counts)) if candidate_counts else None,
                "target_summary": target_summary,
                "engine_names": engine_names,
                "workstation_names": workstations,
                "analysis_job_count": analysis_jobs,
                "analysis_job_statuses": analysis_statuses,
                "last_analysis_activity": max((source.get("last_analysis_activity") or "" for source in sources), default=""),
                "acceptance_observed": acceptance_observed,
                "accepted_count": accepted_count,
                "rejected_count": rejected_count,
                "unresolved_count": unresolved_count,
                "not_evaluated_count": not_evaluated_count,
                "is_local": any(source.get("is_local") for source in sources),
                "updated_at": max((str(source.get("updated_at") or source.get("created_at") or "") for source in sources), default=""),
            }
            # Preserve the source tool for simple filters while retaining the
            # complete engine set for the grouped run row.
            if len(engine_names) == 1:
                row["tool"] = engine_names[0]
                row["display_tool"] = engine_names[0]
            output_rows.append(row)
    return output_rows


def _job_ref(row: dict) -> str:
    if row.get("logical_run_id"):
        return f"{row['logical_run_id']}~{row.get('output_id') or 'results'}"
    return f"{row.get('task_group') or 'pending'}/{row.get('run_id') or row.get('cluster_job_id') or ''}"


def _source_job_ref(row: dict) -> str:
    return f"{row.get('task_group') or 'pending'}/{row.get('run_id') or row.get('cluster_job_id') or ''}"


def _target_system_label(target: dict) -> str:
    label = str(target.get("label") or "")
    if target.get("reference_type") == "structure_backed" and label.startswith("Target sequence · "):
        return label.replace("Target sequence · ", "Target structure · ", 1)
    return label


def _target_label(row: dict) -> str:
    target = row.get("target_summary") if isinstance(row.get("target_summary"), dict) else {}
    if target.get("label"):
        return _target_system_label(target)
    for key in ("target_name", "target_id", "system_name", "molecule_name"):
        value = str(row.get(key) or "").strip()
        if value:
            return value
    summary = row.get("candidate_summary") if isinstance(row.get("candidate_summary"), dict) else {}
    for key in ("target_name", "target_id", "system"):
        value = str(summary.get(key) or "").strip()
        if value:
            return value
    return "Not recorded"


def _read_json(path: Path) -> dict:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
        return value if isinstance(value, dict) else {}
    except (OSError, UnicodeError, ValueError):
        return {}


def _target_from_saved_inputs(run_dir: Path, input_data: dict, metadata: dict) -> dict:
    """Resolve an explicitly supplied sequence target without inventing a pose."""
    from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path
    from mn_molecular_design_studio.molecule_imports import parse_sequences

    mappings = [metadata, input_data, input_data.get("inputs") or {}, input_data.get("params") or {}]
    sequences: list[str] = []
    label = ""
    for mapping in mappings:
        if not isinstance(mapping, dict):
            continue
        label = label or next((str(mapping.get(key) or "").strip() for key in
            ("target_name", "target_label", "target_id", "pdb_id", "receptor") if mapping.get(key)), "")
        value = next((mapping.get(key) for key in
                      ("target_sequence", "target_protein_sequences", "target_fasta") if mapping.get(key)), None)
        if value is None:
            continue
        if isinstance(value, list):
            candidates = [str(item.get("sequence") or "") if isinstance(item, dict) else str(item) for item in value]
            sequences.extend(item.strip().upper() for item in candidates if item.strip())
            continue
        if isinstance(value, dict):
            value = value.get("sequence") or value.get("fasta") or value.get("path") or ""
        text = str(value or "").strip()
        if not text:
            continue
        path = resolve_stored_path(text, run_dir=run_dir, must_exist=True)
        if path and path.is_file():
            try:
                text = path.read_text(encoding="utf-8-sig")
            except (OSError, UnicodeError):
                continue
        try:
            parsed = parse_sequences(text)
            sequences.extend(str(component.get("sequence") or "").upper()
                             for row in parsed for component in row.get("components") or [])
        except (KeyError, TypeError, ValueError):
            compact = "".join(text.split()).upper()
            if compact and all(char in "ACDEFGHIKLMNPQRSTVWYBXZJUO:*-." for char in compact):
                sequences.extend(part for part in compact.split(":") if part)
        if sequences:
            break
    sequences = ["".join(sequence.split()).upper() for sequence in sequences if sequence]
    if not sequences:
        return {}
    joined = ":".join(sequences)
    digest = hashlib.sha256(":".join(sorted(sequences)).encode()).hexdigest()
    return {"identity": "sequence:" + digest, "sequence": joined, "label": label or
            f"Target sequence · {sum(map(len, sequences))} aa · {digest[:8]}"}


def _local_candidate_summary(run_dir: str) -> dict:
    path = Path(run_dir).expanduser().resolve() / "artifacts" / "normalized_candidates" / "candidates.jsonl"
    try:
        stat = path.stat()
    except OSError:
        return {}
    cache_key = f"complex-analysis-candidate-summary-v4-{hashlib.sha256(str(path).encode()).hexdigest()[:20]}"
    cached = _CANDIDATE_SUMMARY_CACHE.get(cache_key)
    if isinstance(cached, dict) and cached.get("mtime_ns") == stat.st_mtime_ns and cached.get("size") == stat.st_size:
        return dict(cached.get("summary") or {})
    cached = st.session_state.get(cache_key)
    if isinstance(cached, dict) and cached.get("mtime_ns") == stat.st_mtime_ns and cached.get("size") == stat.st_size:
        return dict(cached.get("summary") or {})
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates

    candidates = read_candidates(path)
    from mn_molecular_design_studio.target_identity import target_reference_type

    run_path = Path(run_dir).expanduser().resolve()
    input_data = _read_json(run_path / "input.json")
    metadata_data = _read_json(run_path / "metadata.json")
    target_reference = target_reference_type(input_data, metadata_data, run_dir=run_path)
    observed_acceptance = 0
    accepted = rejected = unresolved = not_evaluated = 0
    structures = complexes = 0
    target_systems = {}
    from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path

    def candidate_structure(candidate, run_path: Path) -> tuple[bool, bool]:
        raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
        complex_exists = False
        structure_exists = False
        for key in ("complex_pdb", "binder_pdb", "structure"):
            value = candidate.get(key) or raw.get("structure_file")
            path_value = resolve_stored_path(value, run_dir=run_path, must_exist=True) if value else None
            if path_value is None:
                continue
            try:
                resolved = path_value.resolve(strict=True)
                resolved.relative_to(run_path.resolve())
            except (OSError, ValueError):
                continue
            if resolved.is_file():
                structure_exists = True
                complex_exists = complex_exists or key == "complex_pdb"
        return structure_exists, complex_exists

    def add_target(identity_value: str, sequence: str, name: str = "") -> None:
        identity_value = str(identity_value or "").strip()
        sequence = "".join(str(sequence or "").split()).upper()
        if not identity_value and sequence:
            identity_value = "sequence:" + hashlib.sha256(sequence.encode()).hexdigest()
        if not identity_value:
            return
        label = str(name or "").strip() or (
            f"Target sequence · {len(sequence)} aa · {identity_value.rsplit(':', 1)[-1][:8]}" if sequence
            else "Target system")
        target_systems[identity_value] = {"identity": identity_value, "label": label,
            "sequence_length": len(sequence) if sequence else None,
            "identity_kind": "sequence" if sequence else "named target",
            "reference_type": target_reference}

    for candidate in candidates:
        if not isinstance(candidate, dict):
            continue
        raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
        metrics = candidate.get("metrics") if isinstance(candidate.get("metrics"), dict) else {}
        has_structure, has_complex = candidate_structure(candidate, Path(run_dir).expanduser().resolve())
        if has_complex:
            complexes += 1
        if has_structure:
            structures += 1
        # Use only sequence explicitly assigned to target chains. Never guess
        # roles from chain order or treat the designed binder as its target.
        target_chains = [str(value) for value in candidate.get("target_chains") or raw.get("target_chains") or []]
        sequences_by_chain = raw.get("sequences_by_chain") if isinstance(raw.get("sequences_by_chain"), dict) else {}
        target_sequences = [str(sequences_by_chain[chain]) for chain in target_chains if sequences_by_chain.get(chain)]
        if not target_sequences:
            target_sequences = [str(component.get("sequence") or "") for component in raw.get("components") or []
                                if isinstance(component, dict) and component.get("sequence")
                                and str(component.get("role") or "").casefold() == "target"]
        if target_sequences:
            normalized_sequences = ["".join(sequence.split()).upper() for sequence in target_sequences]
            normalized = ":".join(normalized_sequences)
            digest = hashlib.sha256(":".join(sorted(normalized_sequences)).encode()).hexdigest()
            name = next((str(raw.get(key) or "").strip() for key in
                         ("target_name", "target_label", "target_id", "pdb_id") if raw.get(key)), "")
            add_target("sequence:" + digest, normalized, name)
        native = raw.get("validation_acceptance") if isinstance(raw.get("validation_acceptance"), dict) else {}
        status = str(native.get("status") or "").casefold()
        if status in {"accepted", "rejected", "unresolved", "not_evaluated"}:
            observed_acceptance += 1
            accepted += status == "accepted"
            rejected += status == "rejected"
            unresolved += status == "unresolved"
            not_evaluated += status == "not_evaluated"
        elif metrics.get("native_pass_filters") in (True, False):
            observed_acceptance += 1
            accepted += metrics["native_pass_filters"] is True
            rejected += metrics["native_pass_filters"] is False
    # Older generators may preserve an unstructured target only in the job's
    # input record. Inspect explicit target sequence/FASTA fields as a fallback.
    if not target_systems:
        target = _target_from_saved_inputs(run_path, input_data, metadata_data)
        if target:
            add_target(target["identity"], target.get("sequence") or "", target.get("label") or "")
    if not target_systems:
        from mn_molecular_design_studio.target_identity import prediction_input_systems

        target_systems.update({system["identity"]: system for system in prediction_input_systems(input_data)})
    summary = {
        "candidate_count": len(candidates), "structure_count": structures, "complex_count": complexes,
        "acceptance_observed": observed_acceptance, "accepted_count": accepted,
        "rejected_count": rejected, "unresolved_count": unresolved,
        "not_evaluated_count": not_evaluated, "target_systems": list(target_systems.values()),
    }
    cache_entry = {"mtime_ns": stat.st_mtime_ns, "size": stat.st_size, "summary": summary}
    _CANDIDATE_SUMMARY_CACHE[cache_key] = cache_entry
    st.session_state[cache_key] = cache_entry
    return summary


def _source_summary(row: dict) -> dict:
    from mn_molecular_design_studio.logical_runs import _collection_properties

    if row.get("is_local") and row.get("run_dir"):
        summary = _local_candidate_summary(str(row["run_dir"]))
        if summary:
            return summary
        path = Path(str(row["run_dir"])) / "artifacts" / "normalized_candidates" / "candidates.jsonl"
        try:
            stat = path.stat()
            counts, _ = _collection_properties(str(path), stat.st_mtime_ns, stat.st_size)
            return {"candidate_count": counts["total"], "structure_count": counts["structure"],
                    "complex_count": counts["complex"], "acceptance_observed": 0}
        except OSError:
            pass
    candidate_summary = row.get("candidate_summary") if isinstance(row.get("candidate_summary"), dict) else {}
    acceptance = candidate_summary.get("acceptance_counts") if isinstance(candidate_summary.get("acceptance_counts"), dict) else {}
    accepted = candidate_summary.get("accepted_count", acceptance.get("accepted"))
    rejected = candidate_summary.get("rejected_count", acceptance.get("rejected"))
    unresolved = candidate_summary.get("unresolved_count", acceptance.get("unresolved"))
    not_evaluated = candidate_summary.get("not_evaluated_count", acceptance.get("not_evaluated"))
    observed = candidate_summary.get("acceptance_observed")
    if observed is None:
        observed = sum(int(value or 0) for value in (accepted, rejected, unresolved, not_evaluated))
    return {
        "candidate_count": candidate_summary.get("candidate_count", row.get("candidate_count")),
        "structure_count": candidate_summary.get("structure_count", row.get("structure_count")),
        "complex_count": candidate_summary.get("complex_count", row.get("complex_count")),
        "acceptance_observed": int(observed or 0), "accepted_count": accepted,
        "rejected_count": rejected, "unresolved_count": unresolved,
        "not_evaluated_count": not_evaluated,
        "target_systems": candidate_summary.get("target_systems") or row.get("target_systems") or [],
    }


def _format_recorded_acceptance(row: dict) -> str:
    if not int(row.get("acceptance_observed") or 0):
        return "Not recorded"
    parts = []
    for key, label in (("accepted_count", "accepted"), ("rejected_count", "rejected"),
                       ("unresolved_count", "unresolved"), ("not_evaluated_count", "not evaluated")):
        count = int(row.get(key) or 0)
        if count:
            parts.append(f"{count:,} {label}")
    return " · ".join(parts) if parts else "Recorded; status unavailable"


def _table_rows(rows: list[dict], selected_ref: str = "") -> list[dict]:
    result = []
    for row in rows:
        stats = row
        recorded_acceptance = _format_recorded_acceptance(stats)
        structure_count = stats.get("structure_count")
        if structure_count is not None and int(stats.get("unknown_structure_candidates") or 0):
            structure_count = f"{int(structure_count)} confirmed · {int(stats['unknown_structure_candidates'])} remote unverified"
        activity_count = int(row.get("analysis_job_count") or 0)
        coverage = ("Run review open" if _job_ref(row) == selected_ref else
                    f"{activity_count} linked analysis job(s)" if activity_count
                    else "No linked analyses")
        engines = row.get("engine_names") or []
        run_name = row.get("logical_run_name") or row.get("campaign_name") or row.get("display_tool") or row.get("tool") or "Results run"
        run_id = _job_ref(row)
        short_id = run_id.split("~", 1)[0].rsplit(":", 1)[-1].rsplit("/", 1)[-1][-8:]
        updated = row.get("updated_at") or row.get("created_at") or ""
        display_status = ("Analysis running" if set(row.get("analysis_job_statuses") or []) & {
            "queued", "pending", "submitted", "scheduled", "starting", "preparing", "running", "staging", "waiting"
        } else "Review analyses" if activity_count else str(row.get("logical_run_status") or row.get("status") or "Unknown"))
        result.append({
            "Run (ID)": f"{run_name} · {short_id} · {int(row.get('logical_run_job_count') or row.get('source_job_count') or 1)} jobs",
            "Target / system": _target_label(row),
            "Engine": ", ".join(engines) if engines else (row.get("display_tool") or row.get("tool") or row.get("job_type") or "Unknown"),
            "Structures": structure_count if structure_count is not None else "—",
            "Recorded acceptance": recorded_acceptance,
            "Analysis coverage": coverage,
            "Status": display_status,
            "Updated": updated[:10] if updated else "—",
            "Last analysis activity": row.get("last_analysis_activity") or "—",
            "_ref": _job_ref(row),
        })
    return result


def _target_rows(rows: list[dict]) -> tuple[list[dict], int]:
    grouped: dict[str, dict] = {}
    unrecorded_runs = 0
    for row in rows:
        targets = {}
        for source in row.get("run_sources") or []:
            for target in (source.get("source_summary") or {}).get("target_systems") or []:
                identity = str(target.get("identity") or "").strip()
                if identity:
                    targets[identity] = target
        if not targets:
            target = row.get("target_summary") if isinstance(row.get("target_summary"), dict) else {}
            identity = str(target.get("identity") or "").strip()
            if identity and identity != "multiple":
                targets[identity] = target
        if not targets:
            unrecorded_runs += 1
            continue
        for identity, target in targets.items():
            entry = grouped.setdefault(identity, {
                "Target / system": _target_system_label(target) or "Recorded target",
                "Identity": ("sequence:" + identity.partition(":")[2][:12]
                             if identity.startswith("sequence:") else identity),
                "Identity type": target.get("identity_kind") or ("sequence" if identity.startswith("sequence:") else "named target"),
                "Input target": set(), "Reference roles": set(), "Runs": set(), "Result structures": 0, "Engines": set(),
            })
            reference_type = str(target.get("reference_type") or "unknown")
            entry["Input target"].add(reference_type if reference_type in {
                "sequence_only", "structure_backed", "unknown"} else "unknown")
            entry["Reference roles"].add(str(target.get("reference_role") or ""))
            entry["Runs"].add(_job_ref(row))
            entry["Result structures"] += int(row.get("structure_count") or 0)
            entry["Engines"].update(row.get("engine_names") or [])
    result = []
    for entry in grouped.values():
        reference_types = entry["Input target"]
        labels = {
            "sequence_only": "Sequence only",
            "structure_backed": "Structure-backed",
        }
        known = [labels[value] for value in ("sequence_only", "structure_backed") if value in reference_types]
        if len(known) == 2:
            input_target = "Sequence-only and structure-backed runs"
            for prefix in ("Target sequence · ", "Target structure · "):
                if str(entry["Target / system"]).startswith(prefix):
                    entry["Target / system"] = str(entry["Target / system"]).replace(prefix, "Target · ", 1)
        elif known and "unknown" in reference_types:
            input_target = known[0] + " · some runs unrecorded"
        elif known:
            input_target = known[0]
        else:
            input_target = "Not recorded"
        if entry["Reference roles"] == {"sequence_source"}:
            input_target = "PDB sequence source"
        result.append({**entry, "Input target": input_target,
                       "Runs": len(entry["Runs"]),
                       "Engines": ", ".join(sorted(entry["Engines"], key=str.casefold))})
    return sorted(result, key=lambda item: str(item["Target / system"]).casefold()), unrecorded_runs


def _filter_runs(rows: list[dict]) -> list[dict]:
    options = st.columns([2.2, 1.2, 1.1, 1.2, 1.2, 1.3])
    query = options[0].text_input("Search runs", placeholder="Run, exact ID, target, or engine", key="complex-analysis-run-search").casefold().strip()
    targets = sorted({_target_label(row) for row in rows}, key=str.casefold)
    projects = sorted({str(row.get("project_name") or row.get("project_id") or "No project") for row in rows}, key=str.casefold)
    engines = sorted({engine for row in rows for engine in
                      (row.get("engine_names") or [str(row.get("display_tool") or row.get("tool") or row.get("job_type") or "Unknown")])},
                     key=str.casefold)
    workstations = sorted({str(row.get("node_name") or "This workstation") for row in rows}, key=str.casefold)
    target = options[1].selectbox("Target", ["All", *targets], key="complex-analysis-run-target")
    project = options[2].selectbox("Project", ["All", *projects], key="complex-analysis-run-project")
    engine = options[3].selectbox("Engine", ["All", *engines], key="complex-analysis-run-engine")
    workstation = options[4].selectbox("Workstation", ["All", *workstations], key="complex-analysis-run-workstation")
    date_window = options[5].text_input("Created since", placeholder="YYYY-MM-DD", key="complex-analysis-run-date")
    filtered = []
    for row in rows:
        project_value = str(row.get("project_name") or row.get("project_id") or "No project")
        engine_value = str(row.get("display_tool") or row.get("tool") or row.get("job_type") or "Unknown")
        workstation_value = str(row.get("node_name") or "This workstation")
        haystack = " ".join(str(value or "") for value in (
            row.get("logical_run_name"), _job_ref(row), _target_label(row), engine_value,
            " ".join(row.get("engine_names") or []), " ".join(row.get("workstation_names") or []),
        )).casefold()
        created = str(row.get("created_at") or "")[:10]
        if query and query not in haystack:
            continue
        if target != "All" and target != _target_label(row):
            continue
        if project != "All" and project != project_value:
            continue
        if engine != "All" and engine != engine_value and engine not in (row.get("engine_names") or []):
            continue
        if workstation != "All" and workstation != workstation_value:
            continue
        if date_window and (len(date_window) != 10 or created < date_window):
            continue
        filtered.append(row)
    return filtered


def _read_database_links(user: dict, job: dict) -> list[dict]:
    """Read already indexed run associations without triggering Results ingestion."""
    from mn_molecular_design_studio.molecular_database import MolecularDatabase

    database = MolecularDatabase()
    task_group, run_id = str(job.get("task_group") or ""), str(job.get("run_id") or "")
    if not task_group or not run_id:
        return []
    with database.connection() as connection:
        rows = connection.execute(
            "SELECT DISTINCT molecule_id FROM molecule_sources "
            "WHERE json_extract(provenance,'$.task_group')=? AND json_extract(provenance,'$.run_id')=?",
            (task_group, run_id),
        ).fetchall()
    molecule_ids = [str(row["molecule_id"]) for row in rows if row["molecule_id"]]
    if not molecule_ids:
        return []
    return list(database.list_molecules(user, molecule_ids=molecule_ids, collection="all",
                                        limit=min(100_000, len(molecule_ids))).get("items") or [])


def _evaluation_matches_candidate(observation: dict, descriptor: dict, candidate: dict,
                                  structure_digest: str, structure_id: str = "") -> tuple[bool, bool]:
    context = observation.get("context") if isinstance(observation.get("context"), dict) else {}
    saved_candidate = context.get("candidate") if isinstance(context.get("candidate"), dict) else {}
    raw = saved_candidate.get("raw_metadata") if isinstance(saved_candidate.get("raw_metadata"), dict) else {}
    source = raw.get("source_candidate") if isinstance(raw.get("source_candidate"), dict) else {}
    candidate_ids = {
        str(descriptor.get("candidate_id") or ""), str(saved_candidate.get("candidate_id") or ""),
        str(source.get("candidate_id") or ""),
    }
    candidate_ids.update(str(item.get("candidate_id") or "") for item in saved_candidate.get("parents", []) if isinstance(item, dict))
    source_ids = {str(candidate.get("candidate_id") or "")}
    source_ids.update(str(item.get("candidate_id") or "") for item in candidate.get("parents", []) if isinstance(item, dict))
    id_match = bool((candidate_ids - {""}) & (source_ids - {""}))
    if not id_match:
        id_match = any(observed.startswith(requested + "_")
                       for observed in candidate_ids - {""} for requested in source_ids - {""})
    saved_digest = str(context.get("source_structure_digest") or context.get("structure_digest") or
                       descriptor.get("representation") or raw.get("structure_digest") or "")
    saved_structure = str(context.get("source_structure_id") or context.get("representation_id") or
                          descriptor.get("representation_id") or "")
    expected_raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
    expected_model = str(candidate.get("model_index") or expected_raw.get("model_index") or "")
    selection = raw.get("complex_analysis_selection") if isinstance(raw.get("complex_analysis_selection"), dict) else {}
    saved_model = str(selection.get("model_id") or context.get("model_id") or "")
    exact = (bool(structure_id and saved_structure == structure_id) or
             bool(structure_digest and saved_digest == structure_digest and (id_match or not candidate_ids - {""})))
    if expected_model and (not saved_model or saved_model != expected_model):
        exact = False
    same_candidate = id_match or bool(structure_digest and saved_digest == structure_digest)
    return exact, same_candidate


def _database_evaluations(records: list[dict], structures: list[dict]) -> dict[tuple[str, str], list[dict]]:
    from mn_molecular_design_studio.evaluations import evaluation_descriptor

    result: dict[tuple[str, str], list[dict]] = {}
    for item in structures:
        candidate = item.get("candidate") if isinstance(item.get("candidate"), dict) else {}
        digest = str(item.get("digest") or "")
        structure_id = str(item.get("database_structure_id") or "")
        for record in records:
            for observation in record.get("observations") or []:
                try:
                    descriptor = evaluation_descriptor(observation)
                except (KeyError, TypeError, ValueError):
                    continue
                exact, same_candidate = _evaluation_matches_candidate(observation, descriptor, candidate, digest, structure_id)
                if str(descriptor.get("evaluator") or "") == "sequence_qc":
                    saved_candidate = ((observation.get("context") or {}).get("candidate") or {})
                    source_sequence = str(candidate.get("binder_sequence") or "").strip().upper()
                    saved_sequence = str(saved_candidate.get("binder_sequence") or "").strip().upper()
                    exact = bool(source_sequence and saved_sequence == source_sequence)
                    same_candidate = same_candidate or bool(saved_sequence and saved_sequence == source_sequence)
                if exact:
                    result.setdefault((str(item.get("candidate_id") or ""), str(descriptor.get("evaluator") or "")), []).append(
                        {"descriptor": descriptor, "observation": observation, "match": "Exact"})
                elif same_candidate and descriptor.get("evaluator"):
                    result.setdefault((str(item.get("candidate_id") or ""), str(descriptor.get("evaluator") or "")), []).append(
                        {"descriptor": descriptor, "observation": observation, "match": "Stale or nonmatching structure"})
    return result


def _discover_structures(user: dict, job: dict) -> list[dict]:
    from mn_molecular_design_studio.complex_analysis_ui import _job_structure_options, _structure_digest, component_inventory

    options = _job_structure_options(user, job)
    structures = []
    seen = set()
    for option in options:
        path = Path(option["path"])
        candidate = option.get("candidate") if isinstance(option.get("candidate"), dict) else {}
        candidate_id = str(candidate.get("candidate_id") or option.get("candidate_id") or "")
        # job_artifacts also indexes targets and references under input folders.
        # Run-first coverage is about produced candidate poses, so never turn
        # an unassociated input PDB into a candidate or invite its analysis.
        if not candidate_id:
            continue
        digest = _structure_digest(path)
        identity = (candidate_id, digest, str(option.get("relative_path") or ""))
        if identity in seen:
            continue
        seen.add(identity)
        raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
        try:
            components = component_inventory(path.read_bytes(), path.name, raw.get("components"))
        except (OSError, ValueError, TypeError):
            components = []
        structures.append({
            **option, "candidate": candidate, "candidate_id": candidate_id,
            "digest": digest, "model_id": str(candidate.get("model_index") or raw.get("model_index") or path.name),
            "components": components,
        })
    return structures


def _candidate_acceptance_status(candidate: dict) -> str:
    """Return only an explicitly saved candidate acceptance decision."""
    raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
    decision = raw.get("validation_acceptance") if isinstance(raw.get("validation_acceptance"), dict) else {}
    status = str(decision.get("status") or "").casefold()
    if status in {"accepted", "rejected", "unresolved", "not_evaluated"}:
        return status
    metrics = candidate.get("metrics") if isinstance(candidate.get("metrics"), dict) else {}
    native_pass = metrics.get("native_pass_filters")
    if native_pass is True:
        return "accepted"
    if native_pass is False:
        return "rejected"
    return "unknown"


def _filter_structures_by_acceptance(structures: list[dict], selection: str) -> list[dict]:
    """Filter poses by saved acceptance metadata, preserving unknowns as a separate choice."""
    status_groups = {
        "Accepted only": {"accepted"},
        "Rejected only": {"rejected"},
        "Unresolved / not evaluated": {"unresolved", "not_evaluated"},
        "Acceptance unknown": {"unknown"},
    }
    statuses = status_groups.get(selection)
    if statuses is None:
        return list(structures)
    return [item for item in structures
            if _candidate_acceptance_status(item.get("candidate") or {}) in statuses]


def _component_profile(item: dict) -> dict:
    candidate = item.get("candidate") if isinstance(item.get("candidate"), dict) else {}
    raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
    protein_chains = sorted({str(component.get("source_chain") or component.get("chain") or "")
                             for component in item.get("components") or []
                             if component.get("type") == "protein" and (component.get("source_chain") or component.get("chain"))})
    binder = [str(value) for value in candidate.get("binder_chains") or raw.get("binder_chains") or []]
    target = [str(value) for value in candidate.get("target_chains") or raw.get("target_chains") or []]
    roles_valid = bool(binder and target and not set(binder) & set(target)
                       and set(binder + target).issubset(set(protein_chains)))
    ligand_components = [component for component in item.get("components") or [] if component.get("type") == "small_molecule"]
    ligand_identity = None
    if len(ligand_components) == 1:
        component = ligand_components[0]
        smiles_value = component.get("smiles") or raw.get("ligand_smiles") or raw.get("smiles") or ""
        if isinstance(smiles_value, list):
            smiles_value = smiles_value[0] if len(smiles_value) == 1 else ""
        if isinstance(smiles_value, dict):
            smiles_value = smiles_value.get("smiles") or smiles_value.get("value") or ""
        smiles = str(smiles_value).strip()
        ligand_identity = {
            "chain": str(component.get("source_chain") or component.get("chain") or ""),
            "residue_key": str(component.get("residue_key") or ""),
            "ccd": str(component.get("ccd") or ""), "smiles": smiles,
        }
        ligand_identity = {key: value for key, value in ligand_identity.items() if value}
        if not (ligand_identity.get("chain") or ligand_identity.get("residue_key")) or not ligand_identity.get("smiles"):
            ligand_identity = None
    kind = "protein-ligand" if ligand_components and protein_chains else (
        "protein-protein" if len(protein_chains) >= 2 else "standalone protein" if protein_chains else "other")
    return {"kind": kind, "protein_chains": protein_chains, "binder_chains": binder,
            "target_chains": target, "roles_valid": roles_valid,
            "ligand_selection": ligand_identity if len(ligand_components) == 1 else None,
            "ligand_count": len(ligand_components)}


def _analysis_job_index(user: dict, job: dict) -> dict:
    """Index local, source-linked batch analysis jobs once for a selected run."""
    if job.get("run_sources"):
        combined: dict[tuple, list[dict]] = {}
        for source in job.get("run_sources") or []:
            for key, values in _analysis_job_index(user, source).items():
                combined.setdefault(key, []).extend(values)
        return combined
    from mn_molecular_design_studio.auth import accessible_project_ids, can_access_job
    from mn_molecular_design_studio.job_catalog import get_job_access_record, search_jobs

    cache_key = f"complex-analysis-batch-index-{user.get('id')}-{job.get('task_group')}-{job.get('run_id')}"
    cached = st.session_state.get(cache_key)
    if isinstance(cached, dict) and time.monotonic() - float(cached.get("loaded_at") or 0) < 20:
        return dict(cached.get("index") or {})
    tools = ["sequence_qc", "biotite_interface", "pep_patch", "plip", "pandamap", "posebusters"]
    rows = search_jobs(tools=tools, task_group="analysis", limit=1000, sort_by="active_desc",
                       access_owner_id=str(user.get("id") or ""),
                       accessible_project_ids=accessible_project_ids(user),
                       access_admin=user.get("role") == "admin").get("items") or []
    index: dict[tuple, list[dict]] = {}
    for row in rows:
        task_group, run_id = str(row.get("task_group") or ""), str(row.get("run_id") or "")
        if not task_group or not run_id or not can_access_job(user, get_job_access_record(task_group, run_id)):
            continue
        run_dir = Path(str(row.get("run_dir") or ""))
        if not run_dir.is_dir():
            continue
        try:
            metadata = json.loads((run_dir / "metadata.json").read_text())
            if (str(metadata.get("upstream_task_group") or "") != str(job.get("task_group") or "")
                    or str(metadata.get("upstream_run_id") or "") != str(job.get("run_id") or "")):
                continue
            manifest = json.loads((run_dir / "artifacts" / "analysis" / "batch_manifest.json").read_text())
        except (OSError, ValueError, TypeError):
            continue
        tool = str(manifest.get("tool") or "")
        if tool not in tools:
            continue
        selection = manifest.get("selection") if isinstance(manifest.get("selection"), dict) else {}
        status_data = {}
        try:
            status_data = json.loads((run_dir / "artifacts" / "analysis" / "batch_status.json").read_text())
        except (OSError, ValueError, TypeError):
            pass
        status_by_item = {str(item.get("id") or ""): item for item in status_data.get("items", []) if isinstance(item, dict)}
        job_status = str(row.get("status") or metadata.get("status") or "unknown").casefold()
        for item in manifest.get("items") or []:
            if not isinstance(item, dict):
                continue
            outcome = status_by_item.get(str(item.get("id") or ""), {})
            item_status = str(outcome.get("status") or job_status).casefold()
            if item_status in {"completed", "succeeded"} and row.get("status") not in {"completed", "succeeded"}:
                item_status = job_status
            model_id = str(item.get("model_id") or selection.get("model_id") or "")
            item_digest = str(item.get("source_structure_digest") or item.get("structure_digest") or "")
            binder = tuple(sorted(map(str, item.get("binder_chains") or selection.get("binder_chains") or [])))
            target = tuple(sorted(map(str, item.get("target_chains") or selection.get("target_chains") or [])))
            ligand = item.get("ligand_selection") or selection.get("ligand_selection") or {}
            if tool in {"pep_patch", "plip", "pandamap", "posebusters", "sequence_qc"}:
                target = ()
            if tool in {"plip", "pandamap", "posebusters", "sequence_qc"}:
                binder = ()
            if tool == "sequence_qc":
                sequence = str(item.get("sequence") or "").strip().upper()
                item_digest = hashlib.sha256(sequence.encode()).hexdigest() if sequence else ""
                model_id = ""
            settings = manifest.get("settings") if isinstance(manifest.get("settings"), dict) else {}
            relevant_setting = (settings.get("contact_distance") if tool == "biotite_interface" else
                                settings.get("ph") if tool in {"pep_patch", "sequence_qc"} else None)
            key = (tool, str(item.get("candidate_id") or ""),
                   item_digest,
                   model_id, binder, target, json.dumps(ligand, sort_keys=True), relevant_setting)
            index.setdefault(key, []).append({"status": item_status, "run_id": run_id,
                "error": str(outcome.get("error") or metadata.get("error") or ""),
                "created_at": str(row.get("created_at") or "")})
    st.session_state[cache_key] = {"loaded_at": time.monotonic(), "index": index}
    return index


@lru_cache(maxsize=16)
def _qc_runtime_identity(method: str) -> tuple[str, str]:
    """Return the current QC protocol and locally installed container digests."""
    try:
        from mn_molecular_design_studio.engine.workflow_engine.workflows.quality import protocol_revision
        from mn_molecular_design_studio.quality_tools import TOOLS
        image = str(TOOLS[method]["image"])
        inspected = subprocess.run(["docker", "image", "inspect", image, "--format", "{{.Id}}"],
                                   capture_output=True, text=True, check=False, timeout=5)
        image_digest = inspected.stdout.strip() if inspected.returncode == 0 else ""
        return protocol_revision(), image_digest
    except (KeyError, OSError, RuntimeError, subprocess.SubprocessError, TypeError, ValueError):
        return "", ""


def _method_status(database_index: dict, job_index: dict, item: dict, profile: dict,
                   method: str, *, contact_distance: float = 5.0, ph: float = 7.4) -> tuple[str, str]:
    candidate_id = str(item.get("candidate_id") or "")
    digest = str(item.get("digest") or "")
    model_id = "" if method == "sequence_qc" else str(item.get("model_id") or "")
    if method == "sequence_qc":
        candidate = item.get("candidate") if isinstance(item.get("candidate"), dict) else {}
        sequence = str(candidate.get("binder_sequence") or "").strip().upper()
        digest = hashlib.sha256(sequence.encode()).hexdigest() if sequence else ""
    parameters = {
        "sequence_qc": ((), (), ph),
        "biotite_interface": (tuple(sorted(profile.get("binder_chains") or [])), tuple(sorted(profile.get("target_chains") or [])), contact_distance),
        "pep_patch": (tuple(sorted(profile.get("binder_chains") or [])), (), ph),
        "sequence_liability": (
            tuple(sorted(profile.get("binder_chains") or profile.get("protein_chains") or [])),
            tuple(sorted(profile.get("target_chains") or [])), None),
        "plip": ((), (), json.dumps(profile.get("ligand_selection") or {}, sort_keys=True), None),
        "pandamap": ((), (), json.dumps(profile.get("ligand_selection") or {}, sort_keys=True), None),
        "posebusters": ((), (), json.dumps(profile.get("ligand_selection") or {}, sort_keys=True), None),
    }
    expected = parameters.get(method)
    if expected:
        binder, target, setting = expected[:3]
        exact_matches = []
        for value in database_index.get((candidate_id, method), []):
            if value.get("match") != "Exact":
                continue
            context = value["observation"].get("context") or {}
            saved_candidate = context.get("candidate") if isinstance(context.get("candidate"), dict) else {}
            raw = saved_candidate.get("raw_metadata") if isinstance(saved_candidate.get("raw_metadata"), dict) else {}
            selection = raw.get("complex_analysis_selection") if isinstance(raw.get("complex_analysis_selection"), dict) else {}
            settings = context.get("settings") if isinstance(context.get("settings"), dict) else {}
            saved_binder = settings.get("binder_chains") or saved_candidate.get("binder_chains") or []
            saved_target = settings.get("target_chains") or saved_candidate.get("target_chains") or []
            if method != "sequence_qc" and tuple(sorted(map(str, saved_binder))) != binder:
                continue
            if method in {"biotite_interface", "sequence_liability"} and tuple(sorted(map(str, saved_target))) != target:
                continue
            expected_model = "" if method == "sequence_qc" else str(item.get("model_id") or "")
            saved_model = str(selection.get("model_id") or selection.get("model_index")
                              or saved_candidate.get("model_id") or saved_candidate.get("model_index")
                              or context.get("model_id") or "")
            if expected_model and saved_model != expected_model:
                continue
            try:
                if method == "biotite_interface" and (
                        "contact_distance" not in settings or float(settings["contact_distance"]) != contact_distance):
                    continue
                if method in {"pep_patch", "sequence_qc"} and (
                        "ph" not in settings or float(settings["ph"]) != ph):
                    continue
            except (TypeError, ValueError):
                continue
            if method in {"plip", "pandamap", "posebusters"} and json.dumps(
                    selection.get("ligand_selection") or {}, sort_keys=True) != setting:
                continue
            if method in {"sequence_qc", "biotite_interface", "pep_patch", "plip", "pandamap", "posebusters"}:
                expected_protocol, expected_image = _qc_runtime_identity(method)
                saved_protocol = str(settings.get("protocol_revision") or raw.get("analysis_parameters", {}).get("protocol_revision")
                                     if isinstance(raw.get("analysis_parameters"), dict) else settings.get("protocol_revision") or "")
                saved_image = str(raw.get("image_digest") or context.get("image_digest")
                                  or value["descriptor"].get("image_digest") or "")
                if not expected_protocol or saved_protocol != expected_protocol:
                    continue
                if not expected_image or saved_image != expected_image:
                    continue
            exact_matches.append(value)
        completed = [value for value in exact_matches if str(value["descriptor"].get("status") or "").casefold() in {"completed", "succeeded"}
                     and not str((value["observation"].get("context") or {}).get("error") or "")]
        if completed:
            return "Completed", "Compatible Evaluation is saved in the Molecular Database."
        if exact_matches:
            context = exact_matches[0]["observation"].get("context") or {}
            status = str(context.get("status") or exact_matches[0]["descriptor"].get("status") or "Not run").replace("_", " ").title()
            return status, str(context.get("error") or "A matching Evaluation record exists but has not completed.")
    if database_index.get((candidate_id, method)):
        return "Stale or nonmatching", "A same-candidate Evaluation exists, but its structure, model, component selection, or settings do not match this analysis plan."
    if method in {"biotite_interface", "pep_patch", "sequence_qc"}:
        setting = contact_distance if method == "biotite_interface" else ph
        binder = tuple(sorted(profile.get("binder_chains") or [])) if method != "sequence_qc" else ()
        target = tuple(sorted(profile.get("target_chains") or [])) if method == "biotite_interface" else ()
        job_key = (method, candidate_id, digest, model_id, binder, target, "{}", setting)
        matches = job_index.get(job_key) or []
        if matches:
            preferred = next((row for row in matches if row["status"] in {"completed", "succeeded"}), matches[0])
            state = preferred["status"]
            if state in {"completed", "succeeded"}:
                return "Completed job; Database Evaluation not confirmed", "Exact Results analysis job exists; verify its ordinary Evaluation link before counting coverage."
            if state in {"queued", "pending", "submitted", "scheduled", "starting", "preparing", "running", "staging", "waiting"}:
                return "Queued" if state in {"queued", "pending", "submitted", "scheduled", "waiting"} else "Running", "Exact matching analysis job is active."
            if state in {"failed", "cancelled", "canceled"}:
                return "Failed" if state == "failed" else "Not run", preferred.get("error") or "Matching analysis attempt did not complete."
    elif method in {"plip", "pandamap", "posebusters"}:
        job_key = (method, candidate_id, digest, model_id, (), (),
                   json.dumps(profile.get("ligand_selection") or {}, sort_keys=True), None)
        matches = job_index.get(job_key) or []
        if matches:
            preferred = next((row for row in matches if row["status"] in {"completed", "succeeded"}), matches[0])
            state = preferred["status"]
            if state in {"completed", "succeeded"}:
                return "Completed job; Database Evaluation not confirmed", "Exact Results analysis job exists; verify its ordinary Evaluation link before counting coverage."
            if state in {"queued", "pending", "submitted", "scheduled", "starting", "preparing", "running", "staging", "waiting"}:
                return "Queued" if state in {"queued", "pending", "submitted", "scheduled", "waiting"} else "Running", "Exact matching analysis job is active."
            if state in {"failed", "cancelled", "canceled"}:
                return "Failed" if state == "failed" else "Not run", preferred.get("error") or "Matching analysis attempt did not complete."
    return "Not run", "No compatible saved Evaluation or exact active analysis job was found."


def _status_frame(structures: list[dict], database_index: dict, job_index: dict) -> pd.DataFrame:
    methods = ["sequence_qc", "biotite_interface", "pep_patch", "plip", "pandamap", "posebusters", "sequence_liability"]
    rows = []
    for item in structures:
        profile = _component_profile(item)
        structure_label = _structure_label(item)
        candidate = item.get("candidate") if isinstance(item.get("candidate"), dict) else {}
        has_sequence = bool(str(candidate.get("binder_sequence") or "").strip())
        for method in methods:
            status = "Not applicable"
            detail = "This method does not apply to the detected molecular content."
            if method == "sequence_qc" and has_sequence:
                status, detail = _method_status(database_index, job_index, item, profile, method)
            elif method in {"biotite_interface", "pep_patch"} and profile["kind"] == "protein-protein":
                if not profile["roles_valid"]:
                    status, detail = "Unavailable", "Explicit binder/target chain roles are missing; Studio will not infer them from chain order."
                else:
                    status, detail = _method_status(database_index, job_index, item, profile, method)
            elif method in {"plip", "pandamap", "posebusters"} and profile["kind"] == "protein-ligand":
                if not profile.get("ligand_selection"):
                    status, detail = "Unavailable", "A unique ligand residue and chemical identity are required for this pose."
                else:
                    status, detail = _method_status(database_index, job_index, item, profile, method)
            elif method == "sequence_liability" and profile["kind"] in {"protein-protein", "standalone protein"}:
                liability_chains = profile.get("binder_chains") if profile["kind"] == "protein-protein" else profile.get("protein_chains")
                if not liability_chains or (profile["kind"] == "protein-protein" and not profile["roles_valid"]):
                    status, detail = "Unavailable", "A selected protein component and its exact chain role are required for liability screening."
                elif not item.get("database_structure_id"):
                    status, detail = "Unavailable", "The registered liability workflow requires this exact pose to be linked to a Database structure."
                else:
                    status, detail = _method_status(database_index, job_index, item, profile, method)
                    if status in {"Not run", "Failed"}:
                        status, detail = "Unavailable", (
                            "This registered screen requires an explicit component and construct-format choice in the Database Inspector; "
                            "the current run-first bulk action does not submit it automatically."
                        )
                    elif status == "Stale or nonmatching":
                        detail += " A compatible rerun requires explicit Database Inspector settings."
            rows.append({"Candidate / pose": structure_label,
                         "Molecular content": profile["kind"], "Method": _METHOD_LABELS[method],
                         "Status": status, "Evidence / requirement": detail})
        if profile["kind"] in {"protein-protein", "standalone protein"} and not has_sequence:
            rows.append({"Candidate / pose": structure_label, "Molecular content": profile["kind"],
                         "Method": "Sequence Composition", "Status": "Unavailable",
                         "Evidence / requirement": "The Results candidate has no saved protein sequence for this calculation."})
    return pd.DataFrame(rows, columns=["Candidate / pose", "Molecular content", "Method", "Status", "Evidence / requirement"])


def _structure_label(item: dict) -> str:
    source = item.get("source_job") if isinstance(item.get("source_job"), dict) else {}
    source_suffix = str(source.get("run_id") or "")[-8:]
    label = f"{item.get('candidate_id') or 'Structure'} · {item.get('relative_path') or item.get('label') or ''}"
    return f"{label} · {source_suffix}" if source_suffix else label


def _render_job_detail(user: dict, job: dict, *, workflow_schema=None, project_selector=None, submit=None) -> None:
    from mn_molecular_design_studio.complex_analysis_ui import open_result_complex
    from mn_molecular_design_studio.molecular_database import MolecularDatabase

    st.markdown(f"### {job.get('logical_run_name') or job.get('campaign_name') or job.get('tool') or 'Results run'}")
    source_jobs = list(job.get("run_sources") or [job])
    local_sources = [source for source in source_jobs if source.get("is_local")]
    remote_sources = [source for source in source_jobs if not source.get("is_local")]
    engine_labels = ", ".join(job.get("engine_names") or [job.get("display_tool") or job.get("tool") or "Unknown engine"])
    workstation_labels = ", ".join(job.get("workstation_names") or [job.get("node_name") or "This workstation"])
    st.caption(f"Run `{_job_ref(job)}` · {int(job.get('logical_run_job_count') or len(source_jobs))} related jobs · {engine_labels} · {workstation_labels}")

    structures = []
    seen_structures = set()
    for source in local_sources:
        for item in _discover_structures(user, source):
            item["source_job"] = source
            identity = (str(item.get("candidate_id") or ""), str(item.get("digest") or ""))
            if identity in seen_structures:
                continue
            seen_structures.add(identity)
            structures.append(item)
    remote_structure_count = sum(int((source.get("source_summary") or {}).get("structure_count") or 0)
                                 for source in remote_sources)
    remote_unverified_count = sum(int((source.get("source_summary") or {}).get("candidate_count") or 0)
                                  for source in remote_sources
                                  if (source.get("source_summary") or {}).get("structure_count") is None)
    if not structures and remote_sources:
        remote_names = ", ".join(sorted({str(source.get("node_name") or source.get("node_id") or "another workstation")
                                         for source in remote_sources}))
        st.info(f"This Run has source jobs on {remote_names}, where its structure files remain. Switch to one of those workstations to inspect or analyze its poses. Synchronized metadata reports {remote_structure_count} structures and {remote_unverified_count} candidates whose structure count has not been reported.")
        if st.button("Back to run list", key="complex-analysis-detail-back-remote"):
            st.session_state.pop("complex-analysis-selected-run", None)
            st.rerun()
        return
    if not structures:
        st.warning("The selected Run has no normalized candidate or pose structures available for analysis. Sequence-only and empty jobs are excluded from this list; input/reference PDB files are never counted as generated poses.")
        if st.button("Back to run list", key="complex-analysis-detail-back-remote"):
            st.session_state.pop("complex-analysis-selected-run", None)
            st.session_state["complex-analysis-detail-visible"] = False
            st.rerun()
        return
    all_structures = structures
    candidate_statuses = {}
    for item in all_structures:
        candidate = item.get("candidate") if isinstance(item.get("candidate"), dict) else {}
        candidate_id = str(item.get("candidate_id") or candidate.get("candidate_id") or item.get("relative_path") or "")
        if candidate_id:
            source = item.get("source_job") if isinstance(item.get("source_job"), dict) else {}
            candidate_statuses[f"{_source_job_ref(source)}::{candidate_id}"] = _candidate_acceptance_status(candidate)
    local_acceptance_counts = {status: list(candidate_statuses.values()).count(status)
                                for status in ("accepted", "rejected", "unresolved", "not_evaluated", "unknown")}
    filter_options = ["All structures", "Accepted only", "Rejected only",
                      "Unresolved / not evaluated", "Acceptance unknown"]
    has_saved_acceptance = any(local_acceptance_counts[status]
                               for status in ("accepted", "rejected", "unresolved", "not_evaluated"))
    default_filter = "Accepted only" if local_acceptance_counts["accepted"] else "All structures"
    filter_key = "complex-analysis-acceptance-filter-" + hashlib.sha256(_job_ref(job).encode()).hexdigest()[:12]
    if filter_key not in st.session_state:
        st.session_state[filter_key] = default_filter
    selected_filter = st.selectbox(
        "Candidate acceptance filter", filter_options,
        key=filter_key,
        help="The run list keeps every Results run. This filter narrows the poses used in this review and analysis plan. Saved passing decisions default to accepted candidates; other statuses remain available here.",
        disabled=not has_saved_acceptance,
    )
    structures = _filter_structures_by_acceptance(all_structures, selected_filter)
    if has_saved_acceptance:
        st.caption(
            f"Showing {len(structures):,} of {len(all_structures):,} local pose structures · "
            f"{local_acceptance_counts['accepted']:,} accepted · {local_acceptance_counts['rejected']:,} rejected · "
            f"{local_acceptance_counts['unresolved'] + local_acceptance_counts['not_evaluated']:,} unresolved/not evaluated · "
            f"{local_acceptance_counts['unknown']:,} without a saved decision."
        )
    else:
        st.caption("No per-candidate acceptance decisions are saved with these local pose files; all structures remain available.")
    if not structures:
        st.info("This acceptance filter has no matching pose structures. Choose All structures or another acceptance category to inspect the run.")
        if st.button("Back to run list", key="complex-analysis-detail-back-empty-filter"):
            st.session_state["complex-analysis-detail-visible"] = False
            st.rerun()
        return
    options = structures
    unique_candidates = {}
    for option in options:
        cid = str(option.get("candidate_id") or option.get("relative_path") or "")
        unique_candidates.setdefault(cid, option)
    st.markdown("#### Run inventory")
    a, b, c, d = st.columns(4)
    a.metric("Structures in this view", len(options))
    a.metric("Distinct candidates / poses", len(unique_candidates))
    b.metric("Run status", str(job.get("logical_run_status") or job.get("status") or "Unknown").replace("_", " ").title())
    acceptance_count = int(job.get("acceptance_observed") or 0)
    c.metric("Recorded acceptance", _format_recorded_acceptance(job) if acceptance_count else "Not recorded")
    d.metric("Project", job.get("project_name") or job.get("project_id") or "No project")
    if remote_sources:
        st.info(f"This grouped Run also has {remote_structure_count} reported structures and {remote_unverified_count} candidates with unconfirmed structure counts on other workstations. Those remote files are not included in the local Evaluation inventory.")
    if remote_unverified_count:
        st.caption("Remote structure counts are shown as unverified until each source workstation reports its structure inventory; candidate totals are not treated as structure totals.")

    db = MolecularDatabase()
    database_records_by_id = {}
    for source in source_jobs:
        for record in _read_database_links(user, source):
            database_records_by_id[str(record.get("id") or "")] = record
    database_records = list(database_records_by_id.values())
    source_has_targets = any((source.get("source_summary") or {}).get("target_systems") for source in source_jobs)
    database_target_linked = any(
        isinstance(record.get("metadata"), dict) and record["metadata"].get("design_target_molecule_id")
        for record in database_records
    )
    target_link_missing = bool(source_has_targets and local_sources and not database_target_linked)
    structure_by_id = {}
    for record in database_records:
        for representation in record.get("representations") or record.get("structures") or []:
            context = representation.get("context") if isinstance(representation.get("context"), dict) else {}
            source = item_source = next((candidate.get("source_job") for candidate in structures
                if str(candidate.get("candidate_id") or "") == str(context.get("candidate_id") or "")
                and str(candidate.get("source_job", {}).get("run_id") or "") == str(context.get("run_id") or "")
                and str(candidate.get("source_job", {}).get("task_group") or "") == str(context.get("task_group") or "")), None)
            if source:
                structure_by_id[(str(record.get("id") or ""), str(representation.get("id") or ""))] = {
                    "record": record, "representation": representation, "source_job": source,
                }
    for item in structures:
        context_matches = []
        cid = str(item.get("candidate_id") or "")
        for (molecule_id, structure_id), record_and_rep in structure_by_id.items():
            rep_context = record_and_rep["representation"].get("context") or {}
            source = item.get("source_job") or {}
            if (cid and str(rep_context.get("candidate_id") or "") == cid
                    and str(rep_context.get("run_id") or "") == str(source.get("run_id") or "")
                    and str(rep_context.get("task_group") or "") == str(source.get("task_group") or "")):
                context_matches.append((molecule_id, structure_id))
        item["database_structure_id"] = context_matches[0][1] if len(context_matches) == 1 else ""
        item["database_molecule_id"] = context_matches[0][0] if len(context_matches) == 1 else ""
    database_index = _database_evaluations(database_records, structures)
    job_index = _analysis_job_index(user, job)
    matrix = _status_frame(structures, database_index, job_index)
    applicable = matrix[matrix["Status"] != "Not applicable"] if not matrix.empty else matrix
    coverage_cache = dict(st.session_state.get("complex-analysis-coverage-cache") or {})
    coverage_cache[_job_ref(job)] = {
        "signature": _coverage_signature(job),
        "coverage": {"completed": int((applicable["Status"] == "Completed").sum()) if not applicable.empty else 0,
                     "total": len(applicable),
                     "local_only": bool(remote_structure_count or remote_unverified_count)},
    }
    st.session_state["complex-analysis-coverage-cache"] = coverage_cache
    linked = sum(bool(item.get("database_structure_id")) for item in structures)
    by_structure = matrix.groupby("Candidate / pose")["Status"].agg(list) if not matrix.empty else pd.Series(dtype=object)
    applicable_by_structure = [[value for value in statuses if value != "Not applicable"] for statuses in by_structure]
    complete = int(sum(bool(statuses) and all(value == "Completed" for value in statuses)
                       for statuses in applicable_by_structure))
    analyzed = int(sum(any(value == "Completed" for value in statuses) for statuses in by_structure))
    partial = int(sum(any(value == "Completed" for value in statuses)
                      and any(value not in {"Completed", "Not applicable"} for value in statuses)
                      for statuses in by_structure))
    pending_failed = int(sum(any(value in {"Queued", "Running", "Failed", _UNCONFIRMED_JOB}
                                 for value in statuses) for statuses in by_structure))
    unavailable = int(sum(any(value == "Unavailable" for value in statuses) for statuses in by_structure))
    stale = int(sum(any(value == "Stale or nonmatching" for value in statuses) for statuses in by_structure))
    no_applicable = int(sum(not statuses for statuses in applicable_by_structure))
    not_analyzed = max(0, len(structures) - analyzed - pending_failed - unavailable - stale - no_applicable)
    st.markdown("#### Saved analysis inventory")
    cols = st.columns(8)
    cols[0].metric("Eligible structures", len(structures))
    cols[1].metric("Linked to Database", linked)
    cols[2].metric("Already analyzed", analyzed)
    cols[3].metric("Complete minimum", complete)
    cols[4].metric("Partial / pending / failed", partial + pending_failed)
    cols[5].metric("Not analyzed", not_analyzed)
    cols[6].metric("Stale / nonmatching", stale)
    cols[7].metric("Unavailable / N/A", unavailable + no_applicable)
    st.caption(f"Exact coverage summary: {analyzed} with a compatible saved Evaluation · {complete} complete minimum · {partial} partial · {not_analyzed} not analyzed · {pending_failed} queued/running/failed/unconfirmed · {stale} requiring a matching rerun · {unavailable} requiring unavailable inputs · {no_applicable} with no applicable registered minimum.")
    if matrix.empty:
        st.info("The run has structures, but no minimum analysis profile is currently applicable to the detected molecular content.")
    else:
        st.dataframe(matrix, hide_index=True, use_container_width=True, key="complex-analysis-evaluation-inventory")

    if linked < len(structures) or target_link_missing:
        reason = ("At least one exact Results structure lacks a unique Database representation link"
                  if linked < len(structures) else "The Results run has an explicit target sequence that is not linked to a Database design record")
        st.info(reason + ". Opening this page did not import anything. The explicit link action uses Studio's idempotent importer and stores unstructured targets as sequence-only molecule identities.")
        if st.button("Link missing run structures to the Molecular Database", type="secondary", key="complex-analysis-link-run"):
            try:
                from mn_molecular_design_studio.molecular_database_ingest import ingest_run
                imported = 0
                imported_jobs = 0
                for source in local_sources:
                    source_dir = Path(str(source.get("run_dir") or ""))
                    if source_dir.is_dir():
                        imported += ingest_run(source_dir, database=db, force=bool(database_records))
                        imported_jobs += 1
                st.session_state.pop(f"complex-analysis-run-catalog-{user.get('id') or ''}", None)
                st.success(f"Indexed {imported_jobs} source job(s) in this Run. {imported} Database record(s) were added or updated; matching target sequences are stored as sequence identities when no target structure exists.")
                st.rerun()
            except (OSError, PermissionError, ValueError, RuntimeError) as exc:
                st.error(f"Could not link this run to the Molecular Database: {exc}")

    st.markdown("#### Analyze missing minimum")
    st.caption("Studio proposes only methods supported by the saved molecular inputs. No analysis starts until you submit the plan.")
    profiles = [_component_profile(item) for item in structures]
    ppi_ready = [item for item, profile in zip(structures, profiles)
                 if profile["kind"] == "protein-protein" and profile["roles_valid"] and item.get("candidate_id")]
    ppi_unresolved = [item for item, profile in zip(structures, profiles)
                      if profile["kind"] == "protein-protein" and not profile["roles_valid"]]
    ligand_unresolved = [item for item, profile in zip(structures, profiles)
                         if profile["kind"] == "protein-ligand" and not profile.get("ligand_selection")]
    sequence_ready = [item for item, profile in zip(structures, profiles)
                      if profile["kind"] in {"protein-protein", "protein-ligand", "standalone protein"}
                      and item.get("candidate_id")
                      and str((item.get("candidate") or {}).get("binder_sequence") or "").strip()]
    ligand = [item for item, profile in zip(structures, profiles)
              if profile["kind"] == "protein-ligand" and profile.get("ligand_selection")]
    liability_ready = [item for item, profile in zip(structures, profiles)
                       if profile["kind"] == "standalone protein"
                       or (profile["kind"] == "protein-protein" and profile["roles_valid"])]
    def plan_row(method_label: str, items: list[dict], input_rule: str, rationale: str) -> dict:
        labels = {_structure_label(item) for item in items}
        method_rows = matrix[(matrix["Method"] == method_label)
                             & matrix["Candidate / pose"].isin(labels)] if not matrix.empty else pd.DataFrame()
        states = list(method_rows["Status"]) if not method_rows.empty else []
        return {
            "Method": method_label, "Eligible": len(items),
            "Reuse completed": sum(value == "Completed" for value in states),
            "Submit missing": sum(_needs_execution(value) for value in states),
            "Queued / running": sum(value in {"Queued", "Running", _UNCONFIRMED_JOB} for value in states),
            "Unavailable": sum(value == "Unavailable" for value in states),
            "Input rule": input_rule, "Why baseline": rationale,
        }
    plan_rows = []
    if sequence_ready:
        plan_rows.append(plan_row("Sequence Composition", sequence_ready,
            "Saved protein sequence; no structure or interface assumption.",
            "Records sequence length, charge/pI, GRAVY and related composition measures."))
    if ppi_ready:
        plan_rows.extend([
            plan_row(_METHOD_LABELS["biotite_interface"], ppi_ready,
                "Explicit, non-overlapping binder/target chains; 5 Å contacts.",
                "Quantifies pose contact geometry and residue maps; it does not establish experimental binding."),
            plan_row(_METHOD_LABELS["pep_patch"], ppi_ready,
                "Explicit binder chain group; default pH 7.4.",
                "Characterizes modeled hydrophobic/electrostatic surface patches, not binding affinity."),
        ])
    if ppi_unresolved:
        plan_rows.extend([
            {"Method": _METHOD_LABELS["biotite_interface"], "Eligible": len(ppi_unresolved),
             "Reuse completed": 0, "Submit missing": 0, "Queued / running": 0,
             "Unavailable": len(ppi_unresolved), "Input rule": "Explicit component roles are missing.",
             "Why baseline": "Skipped safely; no chain-order heuristic is used."},
            {"Method": _METHOD_LABELS["pep_patch"], "Eligible": len(ppi_unresolved),
             "Reuse completed": 0, "Submit missing": 0, "Queued / running": 0,
             "Unavailable": len(ppi_unresolved), "Input rule": "Explicit component roles are missing.",
             "Why baseline": "Skipped safely; no chain-order heuristic is used."},
        ])
    if liability_ready:
        plan_rows.append(plan_row(_METHOD_LABELS["sequence_liability"], liability_ready,
            "Exact Database representation plus an explicit component and construct-format selection.",
            "Screens sequence motifs and structure-supported liabilities; the current registered workflow is a single-structure Database action."))
    for method in ("plip", "pandamap", "posebusters"):
        if ligand:
            rationale = ("Finds named protein–ligand contacts." if method in {"plip", "pandamap"}
                         else "Checks ligand pose geometry and chemistry.")
            plan_rows.append(plan_row(_METHOD_LABELS[method], ligand,
                "Exact ligand residue and chemical identity for every pose.", rationale))
        if ligand_unresolved:
            plan_rows.append({"Method": _METHOD_LABELS[method], "Eligible": len(ligand_unresolved),
                "Reuse completed": 0, "Submit missing": 0, "Queued / running": 0,
                "Unavailable": len(ligand_unresolved),
                "Input rule": "A unique ligand residue and chemical identity are missing.",
                "Why baseline": "Reported unavailable until the exact input identity is established."})
    if plan_rows:
        st.dataframe(pd.DataFrame(plan_rows), hide_index=True, use_container_width=True)
    with st.expander("Optional and currently input-dependent methods", expanded=False):
        st.markdown(
            "Rosetta interface energy, ipSAE, reference placement, environmental context and tag/fusion assessment remain optional or unavailable unless their exact required inputs are present. "
            "ipSAE requires a compatible pairwise PAE matrix; reference and environment workflows require validated reference context. "
            "Sequence liability screening is Database-only, single-structure, and requires an explicit exact component and construct-format selection, so it is not submitted by this bulk run-first action. "
            "No missing PAE, ligand identity, receptor or topology is approximated."
        )
    if not job.get("project_id"):
        st.warning("New Results analyses require a source project assignment. This run has none, so the plan is view-only until Studio can preserve its access boundary.")
    elif workflow_schema is None or submit is None:
        st.info("The existing workflow registry or scheduler submission callback is unavailable in this page context.")
    elif ppi_ready or sequence_ready or ligand:
        groups: dict[tuple, list[dict]] = {}
        for item in ppi_ready:
            profile = _component_profile(item)
            source_key = _source_job_ref(item.get("source_job") or {})
            group_key = (source_key, tuple(sorted(profile["binder_chains"])), tuple(sorted(profile["target_chains"])))
            groups.setdefault(group_key, []).append(item)
        st.caption(f"Prepared from the whole selected Run: {len(groups)} exact protein chain-group batch(es), {len(sequence_ready)} sequence candidate(s), and {len(ligand)} exact ligand pose(s). Batches stay attached to their original source job, and compatible Evaluations are reused.")
        submit_minimum = st.button("Submit missing minimum", type="primary", key="complex-analysis-submit-minimum",
                                   disabled=not job.get("project_id"))
        if submit_minimum:
            errors = []
            submitted = 0
            project_id = str(job.get("project_id") or "")

            def submit_batch(source: dict, workflow_id: str, label: str, parameters: dict) -> None:
                nonlocal submitted
                schema = workflow_schema(workflow_id)
                submit(schema, workflow_id, {"source_job": _source_job_ref(source), **parameters}, project_id, {}, label)
                submitted += 1

            # Keep candidates from separate source jobs in separate scheduler
            # requests because identical candidate IDs may occur in different
            # batches and each exact Results manifest lives with its source job.
            sequence_groups: dict[str, list[dict]] = {}
            for item in sequence_ready:
                sequence_groups.setdefault(_source_job_ref(item.get("source_job") or {}), []).append(item)
            source_lookup = {_source_job_ref(source): source for source in local_sources}
            for source_ref, values in sequence_groups.items():
                unique_values = list({str(item["candidate_id"]): item for item in values}.values())
                missing = [item for item in unique_values if _needs_execution(
                    _method_status(database_index, job_index, item, {}, "sequence_qc")[0])]
                if missing:
                    try:
                        submit_batch(source_lookup[source_ref], "analysis.sequence_qc", "Sequence Composition",
                                     {"candidate_ids": [str(item["candidate_id"]) for item in missing], "ph": 7.4})
                    except (OSError, RuntimeError, TypeError, ValueError) as exc:
                        errors.append(f"Sequence Composition ({source_ref}): {exc}")

            for (source_ref, binder, target), values in groups.items():
                unique_values = list({str(item["candidate_id"]): item for item in values}.values())
                model_ids = {str(item["candidate_id"]): str(item.get("model_id") or "") for item in unique_values}
                for workflow_id, label, method in (
                    ("analysis.biotite_interface", "Biotite Interface", "biotite_interface"),
                    ("analysis.pep_patch", "PEP-Patch", "pep_patch"),
                ):
                    needed = [item for item in unique_values if _needs_execution(_method_status(
                        database_index, job_index, item,
                        {"binder_chains": list(binder), "target_chains": list(target)}, method,
                    )[0])]
                    if not needed:
                        continue
                    needed_ids = [str(item["candidate_id"]) for item in needed]
                    parameters = {"candidate_ids": needed_ids,
                                  "selected_model_ids": {key: model_ids[key] for key in needed_ids}}
                    if method == "biotite_interface":
                        parameters.update(binder_chains=list(binder), target_chains=list(target), contact_distance=5.0)
                    else:
                        parameters.update(binder_chains=list(binder), ph=7.4)
                    try:
                        submit_batch(source_lookup[source_ref], workflow_id, label, parameters)
                    except (OSError, RuntimeError, TypeError, ValueError) as exc:
                        errors.append(f"{label} ({source_ref}): {exc}")

            ligand_groups: dict[str, list[dict]] = {}
            for item in ligand:
                ligand_groups.setdefault(_source_job_ref(item.get("source_job") or {}), []).append(item)
            for source_ref, values in ligand_groups.items():
                ligand_values = list({str(item["candidate_id"]): item for item in values}.values())
                for workflow_id, label, method in (
                    ("analysis.plip", "PLIP", "plip"),
                    ("analysis.pandamap", "PandaMap", "pandamap"),
                    ("analysis.posebusters", "PoseBusters", "posebusters"),
                ):
                    needed = [item for item in ligand_values if _needs_execution(_method_status(
                        database_index, job_index, item, _component_profile(item), method,
                    )[0])]
                    if not needed:
                        continue
                    needed_ids = [str(item["candidate_id"]) for item in needed]
                    parameters = {
                        "candidate_ids": needed_ids,
                        "selected_model_ids": {str(item["candidate_id"]): str(item.get("model_id") or "") for item in needed},
                        "ligand_selections": {str(item["candidate_id"]): _component_profile(item)["ligand_selection"] for item in needed},
                    }
                    try:
                        submit_batch(source_lookup[source_ref], workflow_id, label, parameters)
                    except (OSError, RuntimeError, TypeError, ValueError) as exc:
                        errors.append(f"{label} ({source_ref}): {exc}")
            if errors:
                st.error("Some analysis batches could not be submitted: " + " · ".join(errors))
            if submitted:
                st.success(f"Submitted {submitted} bounded analysis batch(es) through Studio's registered workflows and scheduler.")
                st.rerun()
            elif not errors:
                st.info("No missing minimum batches remain. Exact pending or completed attempts were preserved.")
    else:
        st.info("No sequence, protein-interface or chemically identified ligand inputs are ready for a baseline batch.")

    selected_option_refs = [f"{_source_job_ref(item.get('source_job') or {})}::{item.get('relative_path') or ''}" for item in structures]
    selected_ref_key = "complex-analysis-detail-structure"
    selected_path_ref = st.selectbox("Inspect a structure or pose", selected_option_refs,
        format_func=lambda value: next((f"{item.get('candidate_id') or 'Structure'} · {item.get('relative_path') or item.get('label') or ''} · {_source_job_ref(item.get('source_job') or {})}" for item, ref in zip(structures, selected_option_refs) if ref == value), value),
        key=selected_ref_key) if selected_option_refs else ""
    selected_item = next((item for item, ref in zip(structures, selected_option_refs) if ref == selected_path_ref), None)
    actions = st.columns([1, 1, 2])
    if selected_item and actions[0].button("Open detailed Inspector", type="secondary", key="complex-analysis-detail-open-inspector"):
        open_result_complex(Path(str(selected_item["source_job"]["run_dir"])), Path(selected_item["path"]),
                            selected_item.get("candidate"), user=user)
    if actions[1].button("Back to run list", key="complex-analysis-detail-back"):
        st.session_state.pop("complex-analysis-selected-run", None)
        st.rerun()
    if selected_item:
        actions[2].caption(f"Exact artifact: {selected_item.get('relative_path')} · SHA-256 {str(selected_item.get('digest') or '')[:16]}…")


def _coverage_signature(row: dict) -> str:
    return hashlib.sha256(json.dumps([
        row.get("structure_count"), row.get("updated_at"), row.get("last_analysis_activity"),
        row.get("analysis_job_count"), row.get("analysis_job_statuses"),
    ], sort_keys=True).encode()).hexdigest()


def _landing_rows(rows: list[dict], coverage_cache: dict | None = None) -> list[dict]:
    """Send presentation data only; never infer Evaluation coverage from job counts."""
    engine_labels = {"bindcraft2": "BindCraft2", "rfd4_proteina": "RFD4", "rfd4": "RFD4",
                     "optimized_prediction_localfold": "LocalFold", "gnina": "GNINA",
                     "optimized_prediction_opendde_native": "OpenDDE", "openfold3": "OpenFold3"}
    active = {"queued", "pending", "submitted", "scheduled", "starting", "preparing",
              "running", "staging", "waiting"}
    result = []
    for row, table in zip(rows, _table_rows(rows)):
        target = row.get("target_summary") or {}
        reference = str(target.get("reference_type") or "unknown")
        target_label = _target_label(row)
        label_parts = target_label.split(" · ")
        generic_target = target_label.startswith(("Target sequence · ", "Target structure · "))
        reference_label = {"sequence_only": "Sequence only", "structure_backed": "Structure supplied"}.get(reference, "Input type not recorded")
        if target.get("reference_role") == "sequence_source":
            reference_label = "PDB sequence source"
        if generic_target:
            target_name = "Structured target" if reference == "structure_backed" else "Sequence target"
        else:
            target_name = "Not recorded" if target_label == "Target not recorded" else label_parts[0]
        target_detail = " · ".join(label_parts[1:2])
        if target.get("identity_kind") == "system":
            target_detail = f"{int(target.get('chain_count') or 1)} protein chain(s)"
        if reference != "unknown":
            target_detail = " · ".join(filter(None, [target_detail, reference_label]))
        statuses = {str(value).casefold() for value in row.get("analysis_job_statuses") or []}
        activity_count = int(row.get("analysis_job_count") or 0)
        cached = (coverage_cache or {}).get(_job_ref(row)) or {}
        coverage = cached.get("coverage") if cached.get("signature") == _coverage_signature(row) else None
        if statuses & active:
            status, tone = "In progress", "blue"
        elif statuses & {"failed", "cancelled", "canceled"}:
            status, tone = "Needs attention", "amber"
        elif coverage and coverage["total"] and coverage["completed"] == coverage["total"] and not coverage["local_only"]:
            status, tone = "Analyzed", "green"
        elif coverage and coverage["total"] and not coverage["completed"]:
            status, tone = "Missing analysis", "blue"
        elif coverage or activity_count:
            status, tone = "Review results", "blue"
        else:
            status, tone = "Not checked", ""
        sources = row.get("run_sources") or []
        exact_id = str(sources[0].get("run_id") or "") if len(sources) == 1 else _job_ref(row).split("~", 1)[0].rsplit(":", 1)[-1][-8:]
        result.append({
            "id": _job_ref(row), "name": str(row.get("logical_run_name") or row.get("campaign_name") or row.get("tool") or "Results run"),
            "run_id": exact_id, "jobs": int(row.get("logical_run_job_count") or len(sources) or 1),
            "target": target_label, "target_name": target_name, "target_detail": target_detail,
            "target_kind": reference, "target_reference": reference_label,
            "engines": [engine_labels.get(str(value).casefold(), str(value)) for value in
                        (row.get("engine_names") or [row.get("tool") or "Unknown"])],
            "structures": int(row.get("structure_count") or 0),
            "unverified": int(row.get("unknown_structure_candidates") or 0),
            "acceptance": table["Recorded acceptance"], "analysis_jobs": activity_count,
            "coverage": coverage, "status": status, "status_tone": tone,
            "project": str(row.get("project_name") or row.get("project_id") or "No project"),
            "workstations": list(row.get("workstation_names") or [row.get("node_name") or "This workstation"]),
            "created": str(row.get("created_at") or ""), "updated": str(row.get("updated_at") or ""),
            "last_analysis": str(row.get("last_analysis_activity") or ""),
        })
    return result


def _apply_landing_event(event: dict | None, rows: list[dict], state) -> bool:
    """Validate client navigation against the current authorized run catalog."""
    if not isinstance(event, dict) or not isinstance(event.get("token"), str):
        return False
    token = event["token"][:128]
    if not token or token == state.get("complex-analysis-landing-event"):
        return False
    known = {_job_ref(row) for row in rows}
    ref = str(event.get("ref") or "")
    selected = event.get("selected_refs")
    selected = list(dict.fromkeys(value for value in selected if isinstance(value, str) and value in known)) if isinstance(selected, list) else []
    action = event.get("action")
    if action == "open_run":
        if ref not in known:
            return False
        state["complex-analysis-selected-run"] = ref
        state["complex-analysis-detail-visible"] = True
        state["complex-analysis-detail-intent"] = str(event.get("intent") or "review")
        state["complex-analysis-secondary-view"] = "analyze"
        state["complex-analysis-main-tab"] = "Runs"
    elif action == "compare_runs":
        if not 2 <= len(selected) <= 5:
            return False
        state["complex-analysis-comparison-runs"] = selected
        state["complex-analysis-secondary-view"] = "compare"
        state["complex-analysis-detail-visible"] = False
    elif action == "view":
        view = event.get("view")
        if view not in {"analyze", "compare", "collections", "analyses"}:
            return False
        state["complex-analysis-secondary-view"] = view
        state["complex-analysis-detail-visible"] = False
        state.pop("complex-analysis-comparison-runs", None)
        if view == "analyze":
            state["complex-analysis-main-tab"] = "Runs"
    elif action == "tab":
        tab = event.get("tab")
        if tab not in {"Runs", "Targets", "Molecules", "Collections"}:
            return False
        state["complex-analysis-main-tab"] = tab
        state["complex-analysis-secondary-view"] = "analyze"
        state["complex-analysis-detail-visible"] = False
    elif action == "refresh":
        state["complex-analysis-refresh-requested"] = True
    else:
        return False
    saved = event.get("view_state")
    saved = saved if isinstance(saved, dict) else {}
    view_state = {key: str(saved.get(key) or "")[:500] for key in
                  ("query", "project", "engine", "target", "workstation", "status", "since")}
    view_state.update({"selected_refs": selected, "focus_ref": ref if ref in known else "",
                       "sort": saved.get("sort") if saved.get("sort") in {
                           "name", "target", "engine", "structures", "coverage", "status", "updated"} else "updated",
                       "direction": "asc" if saved.get("direction") == "asc" else "desc",
                       "page": min(max(1, saved["page"]), 10000) if isinstance(saved.get("page"), int) else 1,
                       "more": bool(saved.get("more"))})
    state["complex-analysis-landing-state"] = view_state
    state["complex-analysis-landing-event"] = token
    return True


def _render_secondary(user: dict, mode: str, *, workflow_schema=None, project_selector=None, submit=None,
                      run_rows: list[dict] | None = None) -> None:
    if mode == "collections":
        from mn_molecular_design_studio.complex_analysis_collection_ui import render_collection_workspace
        render_collection_workspace(user)
    elif mode == "compare":
        st.markdown("#### Compare saved structures across runs")
        st.caption("The existing comparison view resolves structures by their Molecular Database identity and target/molecule provenance.")
        from mn_molecular_design_studio.complex_analysis_comparison_ui import render_comparison
        selected_refs = st.session_state.get("complex-analysis-comparison-runs") or []
        by_ref = {_job_ref(row): row for row in run_rows or []}
        groups = [by_ref[ref] for ref in selected_refs if ref in by_ref]
        if groups:
            st.caption("Selected runs: " + " · ".join(str(group.get("logical_run_name") or "Run") for group in groups))
        render_comparison(user, workflow_schema=workflow_schema, project_selector=project_selector,
                          submit=submit, run_groups=groups or None)
    elif mode == "analyses":
        st.markdown("#### Existing analysis activity")
        from mn_molecular_design_studio.logical_runs import catalog
        runs = catalog(user, token=str(st.session_state.get("studio-session-token") or "development-mode"))
        jobs = [job for run in runs for job in run.get("jobs", []) if str(job.get("task_group") or "") == "analysis"]
        if not jobs:
            st.info("No accessible registered analysis jobs were found.")
        else:
            frame = pd.DataFrame([{"Analysis": job.get("campaign_name") or job.get("tool") or "Analysis",
                                   "Exact job": _job_ref(job), "Method": job.get("tool") or job.get("job_type"),
                                   "Status": job.get("status"), "Source run": f"{job.get('upstream_task_group') or ''}/{job.get('upstream_run_id') or ''}",
                                   "Updated": job.get("updated_at") or job.get("completed_at") or job.get("created_at"),
                                   "Workstation": job.get("node_name") or "This workstation"} for job in jobs])
            st.dataframe(frame, hide_index=True, use_container_width=True)


def render_complex_analysis_dashboard(user: dict, *, workflow_schema=None, project_selector=None, submit=None) -> None:
    from mn_molecular_design_studio.complex_analysis_landing import render_landing

    st.html("""<style>
    [data-testid="stMainBlockContainer"]{max-width:1800px;padding:1.35rem 1.6rem 3rem}
    [data-testid="stHeader"]{height:0;background:transparent}
    [data-testid="stToolbar"]{display:none}
    [data-testid="stSidebar"]{min-width:240px!important;max-width:240px!important;background:linear-gradient(120deg,#f7f9fc,#fafbfc);border-right:1px solid #edf0f5}
    [data-testid="stSidebarHeader"]{height:12px;min-height:12px}
    [data-testid="stSidebarUserContent"]{padding-top:.25rem}
    [data-testid="stSidebarCollapseButton"]{position:absolute;right:3px;top:6px;opacity:.4}
    [data-testid="stSidebarCollapseButton"]:hover{opacity:1}
    [data-testid="stSidebar"] h3{font-size:1rem;letter-spacing:-.3px}
    [data-testid="stSidebar"] [data-testid="stVerticalBlock"]{gap:.8rem}
    [data-testid="stSidebar"] hr{margin:1.1rem 0}
    iframe[title*="studio_analysis_landing"]{display:block;border:0;color-scheme:light}
    @media(max-width:800px){[data-testid="stMainBlockContainer"]{padding:1rem}}
    </style>""")
    next_tab = st.session_state.pop("complex-analysis-next-tab", None)
    if next_tab:
        st.session_state["complex-analysis-main-tab"] = next_tab
    refresh = bool(st.session_state.pop("complex-analysis-refresh-requested", False))
    try:
        runs = _run_catalog(user, refresh=refresh)
        rows = _result_sources(runs)
    except (OSError, RuntimeError, ValueError, PermissionError) as exc:
        st.error(f"Could not load accessible Results runs: {exc}")
        return
    selected_ref = str(st.session_state.get("complex-analysis-selected-run") or "")
    selected = next((row for row in rows if _job_ref(row) == selected_ref), None)
    if st.session_state.get("complex-analysis-detail-visible") and selected:
        st.title("Complex Analysis")
        if st.button("← Back to runs", key="complex-analysis-run-detail-back-top"):
            st.session_state["complex-analysis-detail-visible"] = False
            st.rerun()
        if st.session_state.get("complex-analysis-detail-intent") == "additional":
            st.info("Review the run inventory below. To configure additional analyses for a specific pose, choose a structure and open its detailed Inspector.")
        _render_job_detail(user, selected, workflow_schema=workflow_schema,
                           project_selector=project_selector, submit=submit)
        return
    st.session_state["complex-analysis-detail-visible"] = False
    selected_view = str(st.session_state.get("complex-analysis-secondary-view") or "analyze")
    nav = str(st.session_state.get("complex-analysis-main-tab") or "Runs")
    event = render_landing(
        rows=_landing_rows(rows, st.session_state.get("complex-analysis-coverage-cache")),
        view=nav, selected_view=selected_view,
        view_state=st.session_state.get("complex-analysis-landing-state") or {},
        key="complex-analysis-landing",
    )
    if _apply_landing_event(event, rows, st.session_state):
        st.rerun()
    if selected_view in {"compare", "collections", "analyses"}:
        _render_secondary(user, selected_view, workflow_schema=workflow_schema,
                          project_selector=project_selector, submit=submit, run_rows=rows)
        return
    if nav == "Targets":
        st.subheader("Targets")
        target_rows, unrecorded_count = _target_rows(rows)
        st.caption("Targets are grouped by exact recorded sequence identity. Input target distinguishes a supplied structure from sequence-only design; generated IDR sequence carriers are not treated as target structures.")
        if target_rows:
            target_frame = pd.DataFrame(target_rows)
            target_columns = [column for column in
                              ("Target / system", "Input target", "Identity", "Runs", "Result structures", "Engines")
                              if column in target_frame.columns]
            st.dataframe(target_frame[target_columns], hide_index=True, use_container_width=True)
        else:
            st.info("The accessible structure-bearing runs do not contain a recorded target identity yet.")
        if unrecorded_count:
            st.caption(f"{unrecorded_count} structure-bearing run(s) do not record a separate target identity; some are standalone prediction runs. Chain roles are not inferred from order.")
    elif nav == "Molecules":
        _render_secondary(user, "compare", workflow_schema=workflow_schema,
                          project_selector=project_selector, submit=submit, run_rows=rows)
    elif nav == "Collections":
        _render_secondary(user, "collections", workflow_schema=workflow_schema,
                          project_selector=project_selector, submit=submit, run_rows=rows)
    elif not rows:
        st.info("No accessible Results runs with verified pose structures are available. Empty and sequence-only jobs remain available in Results and the Database.")
