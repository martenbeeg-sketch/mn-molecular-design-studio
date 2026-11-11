"""Traceable, method-specific export of reviewed design cohorts.

The export projects already selected authorized Explorer samples. It does not
retrieve structures, change acceptance, run analyses, or promote Results into
the Molecular Database.
"""
from __future__ import annotations

import csv
import io
import json
import math
from datetime import datetime, timezone
from typing import Any


SHORTLIST_SCHEMA = "studio.experimental-shortlist.v1"
SHORTLISTED = "Shortlisted for experimental testing"
NOT_SHORTLISTED = "Not shortlisted for experimental testing"
REVIEW_STATUSES = (
    "Unreviewed", "Needs review", "Reviewed", SHORTLISTED, NOT_SHORTLISTED, "Not applicable",
)

CSV_FIELDS = (
    "candidate_id", "native_candidate_id", "molecule_id", "structure_id", "model_id", "structure_sha256",
    "generator", "prediction_engine", "target_identity", "target_molecule_id", "target_chains",
    "target_sequences_by_chain", "binder_chains", "binder_sequences_by_chain", "binder_sequence",
    "source_kind", "project_id", "source_reference", "original_artifact_reference", "source_job", "source_workstation", "structure_relative_path",
    "native_acceptance", "studio_acceptance", "final_computational_acceptance", "recorded_acceptance_json",
    "source_lineage_json",
    "review_status", "review_rationale", "evaluation_count", "evaluation_methods", "candidate_metrics_json",
    "evaluations_json", "assessment_coverage_json", "missing_evidence_json", "limitations_json",
)


def _safe(value: Any):
    if isinstance(value, dict):
        return {str(key): _safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_safe(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


def _sequence_by_chain(sample: dict, chain_ids: list[str], *, side: str) -> dict[str, str]:
    tracks = sample.get("tracks") if isinstance(sample.get("tracks"), list) else []
    by_chain = {str(row.get("chain_id") or ""): str(row.get("sequence") or "")
                for row in tracks if isinstance(row, dict)}
    if not by_chain:
        numbering = sample.get(f"{side}_numbering")
        if isinstance(numbering, list):
            by_chain = {str(row.get("chain_id") or ""): str(row.get("sequence") or "")
                        for row in numbering if isinstance(row, dict)}
    return {chain: by_chain[chain] for chain in chain_ids if chain in by_chain}


def _evaluation_record(row: dict, *, model_id: str = "", structure_digest: str = "") -> dict:
    provenance = row.get("provenance") if isinstance(row.get("provenance"), dict) else {}
    parameters = row.get("parameters") or provenance.get("settings") or {}
    source_candidate = provenance.get("candidate") if isinstance(provenance.get("candidate"), dict) else {}
    source_raw = source_candidate.get("raw_metadata") if isinstance(source_candidate.get("raw_metadata"), dict) else {}
    interface_result = source_raw.get("binder_interface_chemistry")
    interface_result = interface_result if isinstance(interface_result, dict) else {}
    residue_evidence = None
    if interface_result:
        residue_evidence = {
            "selected_chains": interface_result.get("selected_chains") or {},
            "contacts": interface_result.get("contacts") or [],
            "hydrogen_bonds": interface_result.get("hydrogen_bonds") or [],
            "residue_annotations": interface_result.get("residue_annotations") or [],
            "methodological_assumptions": interface_result.get("methodological_assumptions") or [],
            "attribution": interface_result.get("attribution") or {},
        }
    source_keys = (
        "source_candidate_id", "source_model_id", "source_structure_digest", "source_run_id",
        "source_job", "source_workstation", "workstation", "source_node", "task_group", "run_id",
        "source_kind", "source_molecule_id", "source_structure_id", "result_reference", "target_id",
        "target_molecule_id", "target_chains", "binder_chains",
    )
    source = {key: provenance[key] for key in source_keys if key in provenance}
    limitations = []
    for container in (row, provenance):
        for key in ("methodological_limitations", "interpretation_limitations", "limitations", "missing_data"):
            value = container.get(key)
            if value:
                limitations.append({"field": key, "value": value})
    if interface_result.get("methodological_assumptions"):
        limitations.append({"field": "methodological_assumptions",
                            "value": interface_result["methodological_assumptions"]})
    return _safe({
        "evaluation_id": row.get("evaluation_id") or row.get("id") or "",
        "observation_id": row.get("observation_id") or "",
        "method": row.get("method") or row.get("evaluator") or "Unknown",
        "method_version": row.get("method_version") or row.get("version")
                         or (parameters.get("method_version") if isinstance(parameters, dict) else None)
                         or "Not recorded",
        "source_engine": row.get("source_engine") or "Not recorded",
        "status": row.get("status") or "completed",
        "run_id": row.get("run_id") or "",
        "candidate_id": row.get("candidate_id") or "",
        "model_id": row.get("source_model_id") or row.get("model_id") or model_id,
        "structure_sha256": row.get("structure_digest") or row.get("source_structure_digest") or structure_digest,
        "parameters": parameters,
        "measurements": row.get("metrics") or {},
        "metric_definitions": row.get("metric_definitions") or {},
        "residue_level_evidence": residue_evidence,
        "evidence_summary": ({"contact_pairs": len(residue_evidence["contacts"]),
                              "hydrogen_bonds": len(residue_evidence["hydrogen_bonds"]),
                              "residue_annotations": len(residue_evidence["residue_annotations"])}
                             if residue_evidence is not None else {}),
        "source_provenance": source,
        "limitations_and_missing_data": limitations,
    })


def _portable_artifact_reference(value: Any) -> str:
    """Keep registered transport URIs; never put a workstation file path in exports."""
    reference = str(value or "").strip()
    if not reference or reference.startswith(("/", "\\", "~")) or (len(reference) > 2 and reference[1:3] == ":\\"):
        return ""
    return reference


def _source_lineage(record: dict, candidate: dict, native_candidate_id: str) -> list[dict]:
    candidate_id = str(candidate.get("candidate_id") or "")
    wanted = {value for value in (candidate_id, native_candidate_id) if value and value != "Not recorded"}
    parents = {str(value) for value in candidate.get("parents") or [] if value}
    wanted.update(parents)
    entries = []
    for source in record.get("sources") or []:
        if not isinstance(source, dict):
            continue
        source_candidate_id = str(source.get("candidate_id") or "")
        source_parents = [str(value) for value in source.get("parents") or [] if value]
        if source_candidate_id not in wanted and not wanted.intersection(source_parents):
            continue
        entries.append({key: source.get(key) for key in (
            "source_key", "candidate_id", "native_candidate_id", "source_tool", "tool", "task_group",
            "run_id", "workstation", "source_node", "parents", "source_generating_run_id",
            "source_candidate_id", "source_link_kind", "job_reference",
        ) if source.get(key) not in (None, "", [])})
    return entries


def build_shortlist_record(sample: dict) -> dict:
    """Create one export row without inventing molecular identities or decisions."""
    candidate = sample.get("candidate") if isinstance(sample.get("candidate"), dict) else {}
    context = sample.get("context") if isinstance(sample.get("context"), dict) else {}
    raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else context
    record = sample.get("record") if isinstance(sample.get("record"), dict) else {}
    structure = sample.get("structure") if isinstance(sample.get("structure"), dict) else {}
    job = sample.get("job") if isinstance(sample.get("job"), dict) else {}
    source_kind = str(sample.get("source_kind") or "")
    if not source_kind and str(sample.get("kind") or "").casefold() in {"database", "db", "results"}:
        source_kind = str(sample.get("kind"))
    source_kind = source_kind or "Unknown"
    is_database = source_kind.casefold() in {"database", "db"}
    source_kind = "Database" if is_database else "Results" if source_kind.casefold() == "results" else source_kind

    candidate_id = str(sample.get("candidate_id") or candidate.get("candidate_id") or context.get("candidate_id") or "")
    source_candidate = raw.get("source_candidate") if isinstance(raw.get("source_candidate"), dict) else {}
    native_candidate_id = str(raw.get("native_candidate_id") or context.get("native_candidate_id")
                              or raw.get("source_candidate_id") or context.get("source_candidate_id")
                              or source_candidate.get("candidate_id") or candidate.get("native_candidate_id") or "Not recorded")
    if is_database and native_candidate_id == candidate_id and not (
        raw.get("native_candidate_id") or context.get("native_candidate_id")
        or raw.get("source_candidate_id") or context.get("source_candidate_id")
        or source_candidate.get("candidate_id")
    ):
        # The Database Explorer's transient CandidateSet view can copy its
        # normalized ID into native_candidate_id as a convenience fallback.
        # Do not export that fallback as an asserted original native ID.
        native_candidate_id = "Not recorded"
    molecule_id = str(sample.get("molecule_id") or record.get("id") or "")
    structure_id = str(sample.get("structure_id") or structure.get("id") or "")
    model_id = str(sample.get("model_id") or context.get("model_id") or context.get("model_index") or "")
    digest = str(sample.get("digest") or sample.get("structure_digest") or structure.get("sha256") or "")
    binder_chains = [str(value) for value in (sample.get("binder_chains") or candidate.get("binder_chains")
                     or raw.get("binder_chains") or context.get("binder_chains") or [])]
    target_chains = [str(value) for value in (sample.get("target_chains") or candidate.get("target_chains")
                     or raw.get("target_chains") or context.get("target_chains") or [])]
    binder_sequences = _sequence_by_chain(sample, binder_chains, side="binder")
    target_sequences = _sequence_by_chain(sample, target_chains, side="target")

    from mn_molecular_design_studio.complex_analysis_collection import acceptance_decisions
    from mn_molecular_design_studio.native_filter_metrics import native_filter_evidence

    decisions = acceptance_decisions(candidate)
    final = next((row for row in decisions if row.get("stage_id") == "final"), None)
    native = native_filter_evidence(candidate) if candidate else {"status": "not recorded"}
    studio = raw.get("validation_acceptance") if isinstance(raw.get("validation_acceptance"), dict) else {}
    evaluations = [_evaluation_record(row, model_id=model_id, structure_digest=digest)
                   for row in sample.get("evaluation_entries", []) if isinstance(row, dict)]
    review = sample.get("review_annotation") if isinstance(sample.get("review_annotation"), dict) else {}
    if not review:
        review = raw.get("human_review") if isinstance(raw.get("human_review"), dict) else {}

    methods_present = sorted({str(row.get("method") or row.get("evaluator") or "Unknown")
                              for row in sample.get("evaluation_entries", []) if isinstance(row, dict)})
    missing_evidence = sample.get("missing_evidence")
    if not isinstance(missing_evidence, list):
        from mn_molecular_design_studio.complex_analysis_batch import ASSESSMENT_METHODS, COMPREHENSIVE_PROFILE
        observed = set(methods_present)
        missing_evidence = [
            {"method_id": method_id, "method": ASSESSMENT_METHODS[method_id]["label"],
             "status": "No matching saved Evaluation listed; applicability/settings are not inferred by export."}
            for method_id in COMPREHENSIVE_PROFILE
            if method_id in ASSESSMENT_METHODS
            and not (observed & set(ASSESSMENT_METHODS[method_id].get("evaluators") or set()))
        ]
    coverage = sample.get("assessment_coverage") if isinstance(sample.get("assessment_coverage"), (dict, list)) else {
        "profile": "Saved Evaluation inventory",
        "recorded_methods": methods_present,
        "note": "Absent methods are not inferred to be applicable, compatible, failed, or measured; use the assessment coverage panel for parameter-aware reuse decisions.",
    }
    limitations = []
    for entry in evaluations:
        limitations.extend(entry.get("limitations_and_missing_data") or [])
    for key in ("interpretation_limitations", "methodological_assumptions", "unsupported_measurements"):
        value = raw.get(key) or context.get(key)
        if value:
            limitations.append({"field": key, "value": value})
    diversity = [entry for entry in evaluations if any(
        token in str(entry.get("method") or "").casefold() for token in ("interface_clustering", "foldseek")
    )]
    evaluation_target_ids = sorted({str((row.get("provenance") or {}).get("target_id")
                                      or (row.get("provenance") or {}).get("target_molecule_id") or "")
                                    for row in sample.get("evaluation_entries", []) if isinstance(row, dict)
                                    and isinstance(row.get("provenance"), dict)} - {""})
    target_identity = (raw.get("target_identity") or raw.get("target_id") or context.get("target_identity")
                       or context.get("target_id") or context.get("design_target_name")
                       or sample.get("target_identity") or (evaluation_target_ids[0] if len(evaluation_target_ids) == 1 else "Not recorded"))
    target_molecule_id = str(raw.get("target_molecule_id") or context.get("target_molecule_id")
                             or context.get("design_target_molecule_id")
                             or (evaluation_target_ids[0] if len(evaluation_target_ids) == 1 else ""))
    artifact_reference = _portable_artifact_reference(
        sample.get("artifact_reference") or context.get("artifact_reference")
        or context.get("structure_uri") or context.get("job_reference")
    )
    source_reference = str(sample.get("source_reference") or sample.get("source_ref") or "")
    if not source_reference and is_database and molecule_id:
        source_reference = f"database:{molecule_id}"
    source_job = str(sample.get("source_job") or "")
    if not source_job and job.get("task_group") and job.get("run_id"):
        source_job = f"{job['task_group']}/{job['run_id']}"
    if not source_job and context.get("task_group") and context.get("run_id"):
        source_job = f"{context['task_group']}/{context['run_id']}"
    project_id = str(sample.get("project_id") or record.get("project_id") or context.get("project_id") or "")
    workstation = str(sample.get("source_workstation") or job.get("node_name") or context.get("workstation")
                      or context.get("source_workstation") or ("Database" if is_database else "Not recorded"))
    generator = str(sample.get("generator_engine") or candidate.get("source_tool") or raw.get("generator_engine")
                    or raw.get("design_engine") or source_candidate.get("source_tool") or "Not recorded")
    prediction = str(sample.get("prediction_engine") or raw.get("prediction_engine") or context.get("prediction_engine") or "Not recorded")
    relative_path = str(sample.get("relative_path") or "")
    structure_reference = {
        "kind": source_kind,
        "molecule_id": molecule_id,
        "structure_id": structure_id,
        "model_id": model_id,
        "sha256": digest,
        "source_reference": source_reference,
        "job_relative_path": relative_path,
    }
    return _safe({
        "candidate_id": candidate_id,
        "native_candidate_id": native_candidate_id,
        "molecule_id": molecule_id,
        "structure_id": structure_id,
        "model_id": model_id,
        "structure_sha256": digest,
        "structure_reference": structure_reference,
        "original_artifact_reference": artifact_reference,
        "generator": generator,
        "prediction_engine": prediction,
        "target": {"identity": target_identity, "molecule_id": target_molecule_id,
                   "chains": target_chains, "sequences_by_chain": target_sequences},
        "binder": {"chains": binder_chains, "sequences_by_chain": binder_sequences,
                   "sequence": candidate.get("binder_sequence") or ""},
        "source": {"kind": source_kind, "project_id": project_id, "reference": source_reference,
                   "job": source_job, "workstation": workstation},
        "source_lineage": _source_lineage(record, candidate, native_candidate_id),
        "acceptance": {"native": native.get("status") or "not recorded",
                       "studio": studio.get("status") or "not evaluated",
                       "final_computational": final.get("status") if final else "not recorded",
                       "recorded_decisions": decisions},
        "human_review": {"status": review.get("status") or "Unreviewed",
                         "rationale": review.get("note") or ""},
        "candidate_metrics": candidate.get("metrics") or {},
        "evaluations": evaluations,
        "assessment_coverage": coverage,
        "missing_evidence": missing_evidence,
        "saved_diversity_assignments": diversity,
        "limitations_and_missing_evidence": limitations,
    })


def build_shortlist_manifest(samples: list[dict], *, exported_at: str | None = None) -> dict:
    """Build a machine-readable record set for exactly the selected source rows."""
    records = [build_shortlist_record(sample) for sample in samples]
    identities = [(row.get("source", {}).get("kind"), row.get("molecule_id"), row.get("structure_id"),
                   row.get("model_id"), row.get("structure_sha256")) for row in records]
    if len(identities) != len(set(identities)):
        raise ValueError("The selected cohort contains duplicate source/model/structure identities.")
    timestamp = exported_at or datetime.now(timezone.utc).isoformat(timespec="seconds")
    return _safe({
        "schema": SHORTLIST_SCHEMA,
        "exported_at": timestamp,
        "record_count": len(records),
        "shortlisted_count": sum(row["human_review"]["status"] == SHORTLISTED for row in records),
        "records": records,
        "interpretation": "Evidence and decisions are method-specific. No experimental affinity or universal quality score is inferred.",
    })


def manifest_json(manifest: dict) -> bytes:
    return (json.dumps(_safe(manifest), ensure_ascii=False, sort_keys=True, indent=2, allow_nan=False) + "\n").encode("utf-8")


def manifest_csv(manifest: dict, *, shortlisted_only: bool = False) -> bytes:
    rows = list(manifest.get("records") or [])
    if shortlisted_only:
        rows = [row for row in rows if (row.get("human_review") or {}).get("status") == SHORTLISTED]
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=CSV_FIELDS, extrasaction="ignore")
    writer.writeheader()
    for row in rows:
        target, binder, source = row.get("target") or {}, row.get("binder") or {}, row.get("source") or {}
        acceptance, review = row.get("acceptance") or {}, row.get("human_review") or {}
        values = {
            "candidate_id": row.get("candidate_id"), "native_candidate_id": row.get("native_candidate_id"),
            "molecule_id": row.get("molecule_id"), "structure_id": row.get("structure_id"),
            "model_id": row.get("model_id"), "structure_sha256": row.get("structure_sha256"),
            "generator": row.get("generator"), "prediction_engine": row.get("prediction_engine"),
            "target_identity": target.get("identity"), "target_molecule_id": target.get("molecule_id"),
            "target_chains": target.get("chains"), "target_sequences_by_chain": target.get("sequences_by_chain"),
            "binder_chains": binder.get("chains"), "binder_sequences_by_chain": binder.get("sequences_by_chain"),
            "binder_sequence": binder.get("sequence"), "source_kind": source.get("kind"),
            "project_id": source.get("project_id"), "source_reference": source.get("reference"),
            "original_artifact_reference": row.get("original_artifact_reference"),
            "source_job": source.get("job"), "source_workstation": source.get("workstation"),
            "structure_relative_path": (row.get("structure_reference") or {}).get("job_relative_path"),
            "native_acceptance": acceptance.get("native"), "studio_acceptance": acceptance.get("studio"),
            "final_computational_acceptance": acceptance.get("final_computational"),
            "recorded_acceptance_json": acceptance.get("recorded_decisions"),
            "source_lineage_json": row.get("source_lineage"),
            "review_status": review.get("status"), "review_rationale": review.get("rationale"),
            "evaluation_count": len(row.get("evaluations") or []),
            "evaluation_methods": [entry.get("method") for entry in row.get("evaluations") or []],
            "candidate_metrics_json": row.get("candidate_metrics"),
            "evaluations_json": [{key: value for key, value in entry.items() if key != "residue_level_evidence"}
                                 for entry in row.get("evaluations") or []],
            "assessment_coverage_json": row.get("assessment_coverage"),
            "missing_evidence_json": row.get("missing_evidence"),
            "limitations_json": row.get("limitations_and_missing_evidence"),
        }
        writer.writerow({key: json.dumps(value, ensure_ascii=False, sort_keys=True, default=str)
                         if isinstance(value, (dict, list)) else value for key, value in values.items()})
    return output.getvalue().encode("utf-8")


def _database_identity(sample: dict) -> tuple[str, str] | None:
    source_kind = str(sample.get("source_kind") or "").casefold()
    if not source_kind and str(sample.get("kind") or "").casefold() in {"database", "db"}:
        source_kind = str(sample.get("kind")).casefold()
    if source_kind not in {"database", "db"}:
        return None
    record = sample.get("record") if isinstance(sample.get("record"), dict) else {}
    structure = sample.get("structure") if isinstance(sample.get("structure"), dict) else {}
    molecule_id = str(sample.get("molecule_id") or record.get("id") or "")
    structure_id = str(sample.get("structure_id") or structure.get("id") or "")
    return (molecule_id, structure_id) if molecule_id and structure_id else None


def render_shortlist_review_export(user: dict, samples: list[dict], *, key_prefix: str = "shortlist") -> None:
    """Add per-design persistent review decisions and two export formats to comparison."""
    import streamlit as st

    st.markdown("#### Experimental shortlist and export")
    st.caption("Review decisions are account-scoped human annotations, separate from computational acceptance and Evaluations. Export preserves source and model provenance; it creates no score or Database promotion.")
    database_samples = [(sample, _database_identity(sample)) for sample in samples]
    database_samples = [(sample, identity) for sample, identity in database_samples if identity]
    from mn_molecular_design_studio.complex_analysis_review_ui import render_review_editor, mark_review_saved
    database, editors, unavailable = None, [], False
    if database_samples:
        from mn_molecular_design_studio.molecular_database import MolecularDatabase
        database = MolecularDatabase()
        for sample, (molecule_id, structure_id) in database_samples:
            structure = sample.get("structure") if isinstance(sample.get("structure"), dict) else {}
            group_id = str(structure.get("group_id") or "")
            identity = {"source_kind": "database", "subject_id": molecule_id, "structure_ref": structure_id}
            try:
                if group_id:
                    sample["project_id"] = str(database.group(group_id, user).get("project_id") or "")
                current = database.get_complex_review(user, **identity)
            except (OSError, PermissionError, ValueError, TypeError) as exc:
                st.warning(f"Cannot load saved review for {sample.get('label') or molecule_id}: {exc}")
                unavailable = True
                continue
            # Export the current persisted value, never a cached or unsaved draft.
            sample["review_annotation"] = current
            with st.expander(f"Review decision · {sample.get('label') or molecule_id}", expanded=False):
                decision, base_key, dirty = render_review_editor(user, identity, current,
                    key_prefix=key_prefix, status_label="Experimental review decision", note_label="Rationale or exclusion reason",
                    placeholder="Record the evidence and practical reason for advancing or excluding this design.")
                if dirty:
                    editors.append((sample, decision, base_key))
    if st.button("Save selected-design review decisions", key=f"{key_prefix}-save-review",
                  disabled=not bool(database_samples) or unavailable):
        errors = [f"{sample.get('label') or decision['subject_id']}: enter a rationale for this decision."
                  for sample, decision, _ in editors
                  if decision["status"] in {SHORTLISTED, NOT_SHORTLISTED} and not decision["note"].strip()]
        if errors:
            for error in errors:
                st.error(error)
        elif editors:
            try:
                saved = database.save_complex_reviews(user, [decision for _, decision, _ in editors])
            except (OSError, PermissionError, ValueError, TypeError) as exc:
                st.error(str(exc))
            else:
                for (sample, _, base_key), review in zip(editors, saved):
                    sample["review_annotation"] = review
                    mark_review_saved(base_key, review)
                st.success("Saved the review decisions against the exact selected Database structures.")
        else:
            st.info("No review changes to save.")
    if unavailable:
        st.info("Review export is unavailable until access to the current saved reviews is restored.")
        return
    results_count = sum(1 for sample in samples if str(sample.get("source_kind") or sample.get("kind") or "").casefold() == "results")
    if results_count:
        st.info("Results-source structures remain in the export with their original job references. This panel persists shortlist annotations only for existing Database structures; it does not promote Results into the Database.")

    manifest = build_shortlist_manifest(samples)
    csv_all = manifest_csv(manifest)
    csv_shortlisted = manifest_csv(manifest, shortlisted_only=True)
    left, middle, right = st.columns(3)
    left.download_button("Download review manifest (JSON)", manifest_json(manifest),
                         "studio-design-review.json", "application/json", key=f"{key_prefix}-json")
    middle.download_button("Download review table (CSV)", csv_all,
                           "studio-design-review.csv", "text/csv", key=f"{key_prefix}-csv")
    right.download_button("Download shortlisted designs (CSV)", csv_shortlisted,
                           "studio-experimental-shortlist.csv", "text/csv", key=f"{key_prefix}-shortlisted-csv",
                           disabled=manifest["shortlisted_count"] == 0)
    st.caption(f"Exported {manifest['record_count']} selected source/model row(s); {manifest['shortlisted_count']} currently have a persisted experimental-shortlist decision.")
