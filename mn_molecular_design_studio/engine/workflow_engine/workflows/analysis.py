from __future__ import annotations

import csv
import copy
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any

import pandas as pd

from mn_molecular_design_studio.engine.workflow_engine.core.candidates import (
    STAGE_ANALYSIS,
    read_candidates,
    write_candidates,
)
from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path
from mn_molecular_design_studio.engine.workflow_engine.core.hotspot_metrics import calculate_hotspot_metrics
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import create_job, finish_job, read_json, update_status, write_json
from mn_molecular_design_studio.engine.workflow_engine.core.scheduler import apply_docker_cpu_limits_to_steps
from mn_molecular_design_studio.measurement_status import finite_number, measurement_outcome, measurement_summary, RANGES


ANALYSIS_GROUP = "analysis"
IPSAE_RUNNER = Path("/opt/mn-ipsae/run_ipsae.py")
IPSAE_FULL_RUNNER = Path("/opt/de_novo_binder_scoring/scripts/ipsae_w_ipae.py")
IPSAE_BATCH_RUNNER = Path(__file__).with_name("ipsae_batch_runner.py")
IPSAE_IMAGE = "mn-protein-scoring-ipsae:latest"
PYROSETTA_IMAGE = "mn-protein-scoring-pyrosetta:latest"
PYROSETTA_WORKDIR = Path("/opt/de_novo_binder_scoring")
PYROSETTA_SCRIPT = PYROSETTA_WORKDIR / "scripts" / "compute_rosetta_metrics.py"


# Acceptance is task-specific; an empty policy is not an affirmative pass.
DEFAULT_THRESHOLDS = {}

METRIC_ALIASES = {
    "binder_plddt": [
        "binder_plddt",
        "af2_binder_plddt",
    ],
    "confidence": [
        "confidence",
        "confidence_score",
        "boltz2_confidence_score",
        "boltz2_confidence",
    ],
    "iptm": [
        "iptm",
        "i_ptm",
        "boltz2_iptm",
        "boltz2_pair_chains_iptm",
    ],
    "ipae": [
        "ipae",
        "i_pae",
        "binder_pae",
        "boltz2_ipae",
    ],
    "ipde": [
        "ipde",
        "boltz2_ipde",
        "boltz2_complex_ipde",
    ],
    "binder_rmsd": [
        "target_aligned_binder_rmsd",
        "binder_rmsd",
        "boltz2_target_aligned_binder_rmsd",
        "design_backbone_rmsd",
    ],
    "monomer_rmsd": [
        "monomer_refolding_rmsd",
        "monomer_rmsd",
        "binder_monomer_rmsd",
    ],
    "ipsae": [
        "ipsae",
        "ipSAE",
        "ipSAE_min",
        "ipsae_min",
    ],
    "ipsae_min": [
        "ipsae_min",
        "ipSAE_min",
    ],
    "ipsae_max": [
        "ipsae_max",
        "ipSAE_max",
    ],
    "ipsae_avg": [
        "ipsae_avg",
        "ipSAE_avg",
    ],
    "lis": [
        "lis",
        "LIS",
    ],
    "ipsae_min_in_calculation": [
        "ipsae_min_in_calculation",
        "ipSAE_min_in_calculation",
    ],
    "ipsae_d0chn": [
        "ipsae_d0chn",
        "ipSAE_d0chn",
    ],
    "ipsae_d0dom": [
        "ipsae_d0dom",
        "ipSAE_d0dom",
    ],
    "pdockq_min": [
        "pdockq_min",
        "pDockQ_min",
    ],
    "pdockq_max": [
        "pdockq_max",
        "pDockQ_max",
    ],
    "pdockq2_min": [
        "pdockq2_min",
        "pDockQ2_min",
    ],
    "pdockq2_max": [
        "pdockq2_max",
        "pDockQ2_max",
    ],
    "hotspot_contact_fraction": ["hotspot_contact_fraction"],
    "hotspots_contacted": ["hotspots_contacted"],
    "hotspot_contact_pairs": ["hotspot_contact_pairs"],
    "min_binder_to_hotspot_distance": ["min_binder_to_hotspot_distance"],
    "mean_binder_to_hotspot_distance": ["mean_binder_to_hotspot_distance"],
    "binder_interface_contacts": ["binder_interface_contacts"],
    "hotspot_interface_contact_fraction": ["hotspot_interface_contact_fraction"],
    "binder_radius_of_gyration": ["binder_radius_of_gyration"],
    "binder_end_to_end_distance": ["binder_end_to_end_distance"],
}

BINDCRAFT_METRICS = [
    "Protocol",
    "Length",
    "Seed",
    "Helicity",
    "Target_Hotspot",
    "InterfaceResidues",
    "pAE",
    "i_pLDDT",
    "ss_pLDDT",
    "Unrelaxed_Clashes",
    "Relaxed_Clashes",
    "Binder_Energy_Score",
    "Surface_Hydrophobicity",
    "ShapeComplementarity",
    "PackStat",
    "dG",
    "dSASA",
    "dG/dSASA",
    "Interface_SASA_%",
    "Interface_Hydrophobicity",
    "n_InterfaceResidues",
    "n_InterfaceHbonds",
    "InterfaceHbondsPercentage",
    "n_InterfaceUnsatHbonds",
    "InterfaceUnsatHbondsPercentage",
    "Interface_Helix%",
    "Interface_BetaSheet%",
    "Interface_Loop%",
    "Binder_Helix%",
    "Binder_BetaSheet%",
    "Binder_Loop%",
    "InterfaceAAs",
    "Target_RMSD",
    "TrajectoryTime",
    "Notes",
]


def _candidate_artifacts(run_dir: Path) -> list[dict[str, str]]:
    artifacts = []
    for path in sorted((run_dir / "artifacts" / "analysis").glob("*")):
        if path.is_file():
            artifacts.append(
                {
                    "name": path.stem,
                    "path": str(path.relative_to(run_dir)),
                    "type": "analysis_table" if path.suffix == ".csv" else "analysis_summary",
                }
            )
    for path in sorted((run_dir / "artifacts" / "normalized_candidates").glob("*")):
        if path.is_file():
            artifacts.append(
                {
                    "name": path.stem,
                    "path": str(path.relative_to(run_dir)),
                    "type": "normalized_candidates" if path.suffix == ".jsonl" else "campaign_result",
                }
            )
    for path in sorted((run_dir / "artifacts" / "validation").glob("*")):
        if path.is_file():
            artifacts.append({
                "name": path.stem,
                "path": str(path.relative_to(run_dir)),
                "type": "validation_report" if path.suffix == ".csv" else "validation_summary",
            })
    return artifacts


def _analysis_source_candidates(source_run_dir: Path, candidates_jsonl: Path) -> list[dict[str, Any]]:
    """Retain inputs and enrich them from their exact saved prediction sidecars."""
    candidates = read_candidates(candidates_jsonl)
    if not candidates:
        raise ValueError("The selected results contain no candidates.")
    for candidate in candidates:
        for key in ("complex_pdb", "binder_pdb", "target_pdb"):
            if candidate.get(key):
                candidate[key] = str(_resolve_path(source_run_dir, candidate[key]))
    # Validation & Acceptance is saved-evidence-only. Recover score summaries
    # and inter-chain PAE from the selected model's own files; never launch a
    # predictor or borrow confidence data from a neighbouring sample.
    from mn_molecular_design_studio.cofolding_results import enrich_candidate

    return [enrich_candidate(source_run_dir, candidate) for candidate in candidates]


def _select_complex_analysis_candidates(
    source_run_dir: Path,
    candidates: list[dict[str, Any]],
    *,
    candidate_ids: list[str] | None = None,
    selected_structure_path: Path | str | None = None,
    selected_model_id: str | None = None,
    binder_chains: list[str] | None = None,
    target_chains: list[str] | None = None,
) -> list[dict[str, Any]]:
    """Apply an exact, authorized Results selection to an existing analysis input."""
    if candidate_ids is not None:
        requested = [str(value).strip() for value in candidate_ids if str(value).strip()]
        if not requested or len(requested) != len(set(requested)):
            raise ValueError("candidate_ids must contain one or more unique candidate IDs.")
        by_id = {str(row.get("candidate_id") or ""): row for row in candidates}
        missing = [identifier for identifier in requested if identifier not in by_id]
        if missing:
            raise ValueError("Selected candidate IDs are unavailable in the source collection: " + ", ".join(missing))
        candidates = [by_id[identifier] for identifier in requested]
    if selected_structure_path is not None and len(candidates) != 1:
        raise ValueError("An exact selected structure can only be applied to one candidate per analysis request.")

    normalized_groups: dict[str, list[str] | None] = {}
    for name, chains in (("binder_chains", binder_chains), ("target_chains", target_chains)):
        if chains is None:
            normalized_groups[name] = None
            continue
        normalized = [str(chain).strip() for chain in chains if str(chain).strip()]
        if not normalized or len(normalized) != len(set(normalized)):
            raise ValueError(f"{name} must contain one or more unique chain identifiers.")
        normalized_groups[name] = normalized
    binder = normalized_groups["binder_chains"]
    target = normalized_groups["target_chains"]
    if binder is not None and target is not None and set(binder) & set(target):
        raise ValueError("The two interface chain groups must not overlap.")

    selected_structure = None
    if selected_structure_path is not None:
        source_root = Path(source_run_dir).expanduser().resolve()
        selected_structure = Path(selected_structure_path).expanduser().resolve(strict=True)
        try:
            selected_structure.relative_to(source_root)
        except ValueError as exc:
            raise ValueError("The selected structure must belong to the authorized source Results job.") from exc
        if not selected_structure.is_file():
            raise ValueError("The selected structure is no longer available.")

    for candidate in candidates:
        if selected_structure is None and binder is None and target is None and not selected_model_id:
            continue
        raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
        if selected_structure is not None:
            candidate["complex_pdb"] = str(selected_structure)
        if binder is not None:
            candidate["binder_chains"] = list(binder)
        if target is not None:
            candidate["target_chains"] = list(target)
        try:
            with Path(candidate.get("complex_pdb") or selected_structure).open("rb") as stream:
                digest = hashlib.file_digest(stream, "sha256").hexdigest()
        except (OSError, TypeError, ValueError):
            digest = ""
        raw["complex_analysis_selection"] = {
            "structure_path": str(selected_structure or candidate.get("complex_pdb") or ""),
            "structure_digest": digest,
            "model_id": str(selected_model_id or raw.get("model_index") or ""),
            "candidate_id": str(candidate.get("candidate_id") or ""),
            "binder_chains": list(candidate.get("binder_chains") or []),
            "target_chains": list(candidate.get("target_chains") or []),
        }
        candidate["raw_metadata"] = raw
    return candidates


def _write_complex_analysis_selection(
    job: Any,
    source_run_dir: Path,
    candidate: dict[str, Any],
    *,
    workflow: str,
    parameters: dict[str, Any],
) -> None:
    """Persist exact pose and chain provenance for reuse and audit."""
    raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
    selection = raw.get("complex_analysis_selection") if isinstance(raw.get("complex_analysis_selection"), dict) else None
    if not selection:
        return
    source_metadata = read_json(Path(source_run_dir) / "metadata.json")
    source_kind = str(source_metadata.get("source_kind") or "results")
    source_manifest = read_json(Path(source_run_dir) / "complex_analysis_source.json") if source_kind == "database" else {}
    write_json(job.run_dir / "artifacts" / "analysis" / "complex_analysis_selection.json", {
        "schema": "studio.complex_analysis_selection.v1",
        "workflow": workflow,
        "source_kind": source_kind,
        "source_task_group": str(source_metadata.get("task_group") or "") if source_kind != "database" else "",
        "source_run_id": str(source_metadata.get("run_id") or "") if source_kind != "database" else "",
        "database_source": source_manifest,
        "candidate_id": selection.get("candidate_id"),
        "structure_path": selection.get("structure_path"),
        "structure_digest": selection.get("structure_digest"),
        "model_id": selection.get("model_id"),
        "binder_chains": selection.get("binder_chains") or [],
        "target_chains": selection.get("target_chains") or [],
        "parameters": parameters,
    })


def _analysis_job(
    source_run_dir: Path,
    candidates_jsonl: Path,
    *,
    tool: str,
    params: dict[str, Any],
    note: str,
    mode: str,
) -> Any:
    source_metadata = read_json(source_run_dir / "metadata.json")
    source_identity = {**read_json(source_run_dir / "molecular_design_studio.json"), **source_metadata}
    job = create_job(
        ANALYSIS_GROUP,
        job_type="analysis",
        tool=tool,
        inputs={"source_run_dir": str(source_run_dir), "candidates_jsonl": str(candidates_jsonl)},
        params={**params, "backend": mode},
    )
    inherited_metadata = {
        "campaign_name": params.get("job_name") or source_metadata.get("campaign_name") or "",
        "campaign_id": source_metadata.get("campaign_id") or source_metadata.get("run_id") or "",
        "upstream_task_group": source_metadata.get("task_group") or source_run_dir.parent.name,
        "upstream_run_id": source_metadata.get("run_id") or source_run_dir.name,
        "upstream_job_code": source_metadata.get("job_code") or "",
    }
    metadata = read_json(job.run_dir / "metadata.json")
    # Direct internal analysis calls retain source ownership as API calls do.
    # Explicit submission ownership and graph/Run context take precedence.
    for field in ("studio_owner_id", "project_id", "studio_project_id", "studio_development_mode"):
        if not metadata.get(field) and source_identity.get(field):
            metadata[field] = source_identity[field]
    source_kind = str(source_metadata.get("source_kind") or "results")
    if source_kind == "database":
        source_manifest = read_json(source_run_dir / "complex_analysis_source.json")
        inherited_metadata.update({
            "upstream_source_kind": "database",
            "upstream_task_group": "",
            "upstream_run_id": "",
            "upstream_job_code": "",
            "database_source": source_manifest,
        })
    elif source_kind == "complex_analysis_collection":
        collection_manifest = read_json(source_run_dir / "complex_analysis_collection.json")
        inherited_metadata.update({
            "upstream_source_kind": "complex_analysis_collection",
            "upstream_task_group": "",
            "upstream_run_id": "",
            "upstream_job_code": "",
            "upstream_collection_id": str(collection_manifest.get("collection_id") or ""),
        })
    else:
        inherited_metadata["upstream_source_kind"] = "results"
    write_json(job.run_dir / "metadata.json", {**metadata, **inherited_metadata})
    write_json(job.run_dir / "command.json", {"mode": mode, "command": [], "note": note})
    update_status(job.run_dir, "running")
    return job


def _analysis_output_candidates(
    candidates: list[dict[str, Any]],
    *,
    tool: str,
) -> list[dict[str, Any]]:
    """Write each tool's scores as a child result without changing upstream ranking data."""
    outputs: list[dict[str, Any]] = []
    for index, source in enumerate(candidates, start=1):
        source_id = str(source.get("candidate_id") or f"candidate_{index}")
        metrics = dict(source.get("metrics") or {})
        metrics.update({"analysis_source_candidate_id": source_id, "analysis_tool": tool})
        outputs.append(
            {
                **source,
                "candidate_id": f"{source_id}_{tool}_{index:03d}",
                "stage": STAGE_ANALYSIS,
                "source_tool": tool,
                "tool": tool,
                "metrics": metrics,
                "parents": [source_id],
                "raw_metadata": {
                    **dict(source.get("raw_metadata") or {}),
                    "result_kind": "app_reanalysis",
                    "analysis_tool": tool,
                    "measurement_provenance": {
                        **dict((source.get("raw_metadata") or {}).get("measurement_provenance") or {}),
                        "tool": tool,
                        "analysis_parameters": (source.get("raw_metadata") or {}).get("analysis_parameters", {}),
                        "source_candidate_id": source_id,
                        "structure": source.get("complex_pdb"),
                        "binder_chains": list(source.get("binder_chains") or []),
                        "target_chains": list(source.get("target_chains") or []),
                        "status": metrics.get("analysis_status", "error" if metrics.get("ipsae_error" if tool == "ipsae" else "rosetta_error") else "unverified"),
                    },
                    "source_candidate": source.get("_measurement_source", source),
                },
            }
        )
    return outputs


def _analysis_metrics_table(candidates: list[dict[str, Any]]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for candidate in candidates:
        row = {
            "candidate_id": candidate.get("candidate_id"),
            "source_tool": candidate.get("source_tool"),
            "stage": candidate.get("stage"),
            "binder_sequence": candidate.get("binder_sequence"),
            "complex_pdb": candidate.get("complex_pdb"),
        }
        row.update(_flatten_metrics(candidate.get("metrics") or {}))
        rows.append(row)
    return pd.DataFrame(rows)


def run_validation_acceptance_analysis(
    source_run_dir: Path,
    candidates_jsonl: Path,
    acceptance_policy: dict[str, Any],
    candidate_ids: list[str] | None = None,
    candidate_limit: int | None = None,
    job_name: str = "Validation & Acceptance",
) -> Path:
    """Filter saved design evidence without invoking any prediction engine.

    Source candidate IDs and generator provenance are preserved in Results.
    Structure-derived interface geometry is measured directly from each saved
    complex, while model scores are read from the candidate's existing native
    metrics and retain their original producing workflow.
    """
    from mn_molecular_design_studio.design_validation import interface_filter_measurements
    from mn_molecular_design_studio.validation_acceptance import apply_acceptance, normalize_policy, write_report

    source_run_dir = Path(source_run_dir)
    candidates_jsonl = Path(candidates_jsonl)
    policy = normalize_policy(acceptance_policy)
    if candidate_limit is not None and (type(candidate_limit) is not int or candidate_limit < 1):
        raise ValueError("candidate_limit must be a positive integer when provided.")
    candidates = _analysis_source_candidates(source_run_dir, candidates_jsonl)
    if candidate_ids is not None:
        requested = [str(value).strip() for value in candidate_ids if str(value).strip()]
        if not requested or len(requested) != len(set(requested)):
            raise ValueError("candidate_ids must contain one or more unique candidate IDs.")
        by_id = {str(row.get("candidate_id") or ""): row for row in candidates}
        missing = [identifier for identifier in requested if identifier not in by_id]
        if missing:
            raise ValueError("Selected candidate ID(s) are absent from this source run: " + ", ".join(missing[:10]))
        candidates = [by_id[identifier] for identifier in requested]
    if candidate_limit is not None:
        candidates = candidates[:candidate_limit]
    if not candidates:
        raise ValueError("The selected source contains no candidates to filter.")

    source_metadata = read_json(source_run_dir / "metadata.json")
    job = _analysis_job(
        source_run_dir, candidates_jsonl, tool="validation_acceptance",
        params={"job_name": job_name,
                "acceptance_policy": policy,
                "candidate_ids": [str(row.get("candidate_id") or "") for row in candidates],
                "candidate_limit": candidate_limit,
                "prediction_performed": False},
        note="Evaluate the configured acceptance rules using existing native measurements and saved structures only. No prediction engine is run.",
        mode="saved-evidence-only",
    )
    prepared = []
    geometry_errors = {}
    for source in candidates:
        row = copy.deepcopy(source)
        metrics = dict(row.get("metrics") or {})
        raw = row.setdefault("raw_metadata", {})
        source_id = str(row.get("candidate_id") or "")
        source_tool = str(row.get("source_tool") or row.get("tool") or "")
        raw["source_run_dir"] = str(source_run_dir.resolve())
        raw["validation_acceptance_provenance"] = {
            "source_candidate_id": source_id,
            "source_tool": source_tool,
            "source_run_id": source_metadata.get("run_id") or source_run_dir.name,
            "method": "saved native metrics + Studio geometry from saved candidate structure",
            "prediction_performed": False,
            "structure": row.get("complex_pdb"),
        }
        structure_path = row.get("complex_pdb") or row.get("binder_pdb")
        if structure_path and row.get("binder_chains") and row.get("target_chains"):
            try:
                measurements = interface_filter_measurements(
                    Path(structure_path), list(row.get("binder_chains") or []),
                    list(row.get("target_chains") or []),
                )
                metrics["validation_interface_residue_count_4a"] = measurements["validation_interface_residue_count_4a"]
                metrics["validation_backbone_ca_clashes_2_5a"] = measurements["validation_backbone_ca_clashes_2_5a"]
                if source_tool == "rfd4_proteina":
                    metrics["rfd4_native_interface_residue_count_4a"] = measurements["validation_interface_residue_count_4a"]
                    metrics["rfd4_native_backbone_ca_clashes_2_5a"] = measurements["validation_backbone_ca_clashes_2_5a"]
            except Exception as exc:
                geometry_errors[source_id] = f"{type(exc).__name__}: {exc}"
                metrics["validation_acceptance_geometry_error"] = geometry_errors[source_id]
        else:
            geometry_errors[source_id] = "Saved complex structure or binder/target chain roles are unavailable."
            metrics["validation_acceptance_geometry_error"] = geometry_errors[source_id]
        row["metrics"] = metrics
        prepared.append(row)

    accepted_rows, summary = apply_acceptance(prepared, policy)
    for row in accepted_rows:
        decision = row.setdefault("raw_metadata", {}).get("validation_acceptance") or {}
        provenance = row["raw_metadata"].setdefault("validation_acceptance_provenance", {})
        provenance.update({
            "criteria_source": "configured Validation & Acceptance policy",
            "measurements": {
                field: {key: value for key, value in observation.items()
                        if key in {"value", "unit", "status", "source_tool", "stage", "method", "passed", "threshold", "operator"}}
                for field, observation in (decision.get("measurements") or {}).items()
            },
        })
    summary.update({
        "tool": "validation_acceptance",
        "source_run_id": source_metadata.get("run_id") or source_run_dir.name,
        "source_tool": source_metadata.get("tool") or "",
        "prediction_performed": False,
        "geometry_measurement_errors": geometry_errors,
        "execution_status": "completed",
        "scientific_status": "partial" if summary["unresolved"] else "valid",
        "source_candidate_ids_preserved": True,
    })
    analysis_dir = job.run_dir / "artifacts" / "analysis"
    analysis_dir.mkdir(parents=True, exist_ok=True)
    acceptance_report = write_report(job.run_dir, accepted_rows, summary)
    table = _analysis_metrics_table(accepted_rows)
    for row in accepted_rows:
        decision = (row.get("raw_metadata") or {}).get("validation_acceptance") or {}
        table.loc[table["candidate_id"] == row.get("candidate_id"), "acceptance_status"] = decision.get("status")
        table.loc[table["candidate_id"] == row.get("candidate_id"), "acceptance_reasons"] = "; ".join(decision.get("reasons") or [])
        for field, observation in (decision.get("measurements") or {}).items():
            table.loc[table["candidate_id"] == row.get("candidate_id"), f"{field} · measured"] = observation.get("value")
            table.loc[table["candidate_id"] == row.get("candidate_id"), f"{field} · passed"] = observation.get("passed")
    table_path = analysis_dir / "acceptance_measurements.csv"
    table.to_csv(table_path, index=False)
    summary_path = job.run_dir / "artifacts" / "validation" / "acceptance_summary.json"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    write_json(summary_path, summary)
    normalized = write_candidates(job.run_dir, "validation_acceptance", accepted_rows)
    artifacts = _candidate_artifacts(job.run_dir)
    finish_job(job.run_dir, True, {
        "outputs": {
            "artifacts": artifacts,
            "candidates": normalized,
            "acceptance_report": str(acceptance_report.relative_to(job.run_dir)),
            "acceptance_measurements": str(table_path.relative_to(job.run_dir)),
        },
        "metrics": {"validation_acceptance": summary},
        "execution_status": "completed",
        "scientific_status": summary["scientific_status"],
        "downstream_artifacts": {
            "candidates_jsonl": "artifacts/normalized_candidates/candidates.jsonl",
            "campaign_result": "artifacts/normalized_candidates/campaign_result.json",
            "acceptance_report": str(acceptance_report.relative_to(job.run_dir)),
            "acceptance_summary": str(summary_path.relative_to(job.run_dir)),
        },
        "summary": summary,
    })
    return job.run_dir


def _finish_metric_analysis(
    job: Any,
    candidates: list[dict[str, Any]],
    *,
    tool: str,
    table_name: str,
    metric_names: list[str],
    required_metrics: list[str] | None = None,
    error_key: str | None = None,
) -> Path:
    analysis_dir = job.run_dir / "artifacts" / "analysis"
    analysis_dir.mkdir(parents=True, exist_ok=True)
    table_path = analysis_dir / table_name
    error_key = error_key or ("ipsae_error" if tool == "ipsae" else "rosetta_error")
    outcomes = []
    required = required_metrics or (["ipsae"] if tool == "ipsae" else ["rosetta_interface_dG"])
    for candidate in candidates:
        metrics = candidate.setdefault("metrics", {})
        outcome = measurement_outcome(metrics, required, metrics.get(error_key))
        structure_measured = tool == "ipsae" and metrics.get("structure_metrics_status") == "measured"
        if tool == "ipsae" and outcome["status"] != "measured" and structure_measured:
            outcome = {**outcome, "status": "partial",
                       "error": str(metrics.get(error_key) or "ipSAE/PAE-derived measurements were unavailable; structure analysis completed")}
        outcomes.append(outcome)
        if outcome["error"]:
            metrics[error_key] = outcome["error"]
        metrics["analysis_status"] = outcome["status"] if outcome["status"] in {"measured", "partial"} else "error"
        candidate.setdefault("raw_metadata", {})["measurement_status"] = outcome
    table = _analysis_metrics_table(candidates)
    table.to_csv(table_path, index=False)
    outputs = _analysis_output_candidates(candidates, tool=tool)
    normalized = write_candidates(job.run_dir, tool, outputs)
    counts = measurement_summary(outcomes)
    if tool == "ipsae":
        partial_count = sum(item["status"] == "partial" for item in outcomes)
        structure_count = sum((candidate.get("metrics") or {}).get("structure_metrics_status") == "measured" for candidate in candidates)
        counts["partial_count"] = partial_count
        counts["structure_measured_count"] = structure_count
        counts["ipsae_measured_count"] = sum(item["status"] == "measured" for item in outcomes)
        counts["scored_count"] += partial_count
        if not counts["scored_count"]:
            counts["scientific_status"] = "failed"
        elif partial_count or counts["error_count"] or counts["missing_count"] or counts["invalid_count"]:
            counts["scientific_status"] = "partial"
        else:
            counts["scientific_status"] = "valid"
    summary = {
        "tool": tool,
        **counts,
        "execution_status": "completed",
        "metrics": metric_names,
        "source_ranking_preserved": True,
    }
    summary_path = analysis_dir / f"{tool}_summary.json"
    write_json(summary_path, summary)
    artifacts = _candidate_artifacts(job.run_dir)
    finish_job(
        job.run_dir,
        counts["scored_count"] > 0,
        {
            "outputs": {"artifacts": artifacts, "candidates": normalized},
            "metrics": counts,
            "execution_status": "completed",
            "scientific_status": counts["scientific_status"],
            "downstream_artifacts": {
                "candidates_jsonl": "artifacts/normalized_candidates/candidates.jsonl",
                "campaign_result": "artifacts/normalized_candidates/campaign_result.json",
            },
            "summary": summary,
        },
    )
    return job.run_dir


def run_reference_placement_contract(
    source_run_dir: Path,
    candidates_jsonl: Path,
    candidate_ids: list[str] | None = None,
    selected_structure_path: Path | None = None,
    selected_model_id: str | None = None,
    source_component_chain: str | None = None,
    reference_context_id: str | None = None,
    reference_molecule_id: str | None = None,
    reference_structure_id: str | None = None,
    reference_structure_digest: str | None = None,
    reference_model_id: str | None = None,
    reference_component_chain: str | None = None,
    method: str | None = None,
    method_version: str | None = None,
) -> Path:
    """Fit one Database complex into its explicitly associated reference frame.

    The Database-source adapter is responsible for authorizing both saved
    structures and staging the existing digest-bound reference association.
    This registered evaluator rechecks that provenance before fitting.
    """
    from mn_molecular_design_studio.complex_analysis_placement import (
        PLACEMENT_METHOD,
        PLACEMENT_VERSION,
        calculate_reference_placement,
    )
    from mn_molecular_design_studio.evaluations import file_digest

    source_run_dir = Path(source_run_dir).expanduser().resolve(strict=True)
    candidates_jsonl = Path(candidates_jsonl).expanduser().resolve(strict=True)
    candidates = _analysis_source_candidates(source_run_dir, candidates_jsonl)
    candidates = _select_complex_analysis_candidates(
        source_run_dir,
        candidates,
        candidate_ids=candidate_ids,
        selected_structure_path=selected_structure_path,
        selected_model_id=selected_model_id,
    )
    if len(candidates) != 1:
        raise ValueError("Reference placement requires exactly one selected Database structure.")
    source_metadata = read_json(source_run_dir / "metadata.json")
    if str(source_metadata.get("source_kind") or "") != "database":
        raise ValueError("Reference placement is available only from an authorized Database source.")
    manifest = read_json(source_run_dir / "complex_analysis_source.json")
    if manifest.get("schema") != "studio.complex_analysis_database_source.v1":
        raise ValueError("The Database source manifest is missing or has an unsupported schema.")
    reference_context = manifest.get("reference_context")
    if not isinstance(reference_context, dict) or reference_context.get("status") != "validated":
        raise ValueError("Reference placement requires an explicitly saved, digest-current, validated association.")
    source_chain = str(source_component_chain or "").strip()
    if not source_chain or source_chain != str(reference_context.get("source_component_chain") or ""):
        raise ValueError("The selected source component does not match the saved reference association.")
    if str(reference_context.get("source_structure_id") or "") != str(manifest.get("structure_id") or ""):
        raise ValueError("The saved reference association points to a different source Database structure.")
    if str(reference_context.get("source_structure_digest") or "") != str(manifest.get("structure_digest") or ""):
        raise ValueError("The saved reference association is stale for the selected source structure.")
    if str(reference_context.get("source_model_id") or "") != str(manifest.get("model_id") or ""):
        raise ValueError("The saved reference association is bound to a different source model.")
    if selected_model_id and str(selected_model_id) != str(manifest.get("model_id") or ""):
        raise ValueError("The selected source model does not match the Database representation.")

    selected_path = Path(str(candidates[0].get("complex_pdb") or "")).resolve(strict=True)
    try:
        selected_path.relative_to(source_run_dir)
    except ValueError as exc:
        raise ValueError("The selected complex must be staged inside the authorized Database source bundle.") from exc
    if file_digest(selected_path) != str(manifest.get("structure_digest") or ""):
        raise ValueError("The selected source structure changed after Database authorization.")
    reference_path = Path(str(reference_context.get("structure_path") or "")).resolve(strict=True)
    try:
        reference_path.relative_to(source_run_dir)
    except ValueError as exc:
        raise ValueError("The associated reference must be staged inside the authorized Database source bundle.") from exc
    if file_digest(reference_path) != str(reference_context.get("reference_structure_digest") or ""):
        raise ValueError("The associated reference structure changed after Database authorization.")

    mapping = reference_context.get("mapping") if isinstance(reference_context.get("mapping"), dict) else {}
    source_identity = {
        "molecule_id": str(manifest.get("molecule_id") or ""),
        "structure_id": str(manifest.get("structure_id") or ""),
        "group_id": str(manifest.get("group_id") or ""),
        "structure_digest": str(manifest.get("structure_digest") or ""),
        "model_id": str(manifest.get("model_id") or ""),
        "chain_id": source_chain,
    }
    reference_identity = {
        "molecule_id": str(reference_context.get("reference_molecule_id") or ""),
        "structure_id": str(reference_context.get("reference_structure_id") or ""),
        "group_id": str(reference_context.get("reference_group_id") or ""),
        "structure_digest": str(reference_context.get("reference_structure_digest") or ""),
        "model_id": str(reference_context.get("reference_model_id") or ""),
        "chain_id": str(reference_context.get("reference_component_chain") or ""),
    }
    expected_public_inputs = {
        "reference_context_id": str(reference_context.get("id") or ""),
        "reference_molecule_id": reference_identity["molecule_id"],
        "reference_structure_id": reference_identity["structure_id"],
        "reference_structure_digest": reference_identity["structure_digest"],
        "reference_model_id": reference_identity["model_id"],
        "reference_component_chain": reference_identity["chain_id"],
        "method": PLACEMENT_METHOD,
        "method_version": PLACEMENT_VERSION,
    }
    requested_public_inputs = {
        "reference_context_id": reference_context_id,
        "reference_molecule_id": reference_molecule_id,
        "reference_structure_id": reference_structure_id,
        "reference_structure_digest": reference_structure_digest,
        "reference_model_id": reference_model_id,
        "reference_component_chain": reference_component_chain,
        "method": method,
        "method_version": method_version,
    }
    if any(value is not None and str(value) != expected_public_inputs[key]
           for key, value in requested_public_inputs.items()):
        raise ValueError("The requested placement identity or method does not match the authorized saved association.")
    if any(not value for identity in (source_identity, reference_identity) for value in identity.values()):
        raise ValueError("The saved reference association is missing a required structure identity field.")
    result = calculate_reference_placement(
        selected_path, reference_path,
        source_chain=source_chain,
        reference_chain=reference_identity["chain_id"],
        source_model_id=source_identity["model_id"],
        reference_model_id=reference_identity["model_id"],
        source_digest=source_identity["structure_digest"],
        reference_digest=reference_identity["structure_digest"],
        mapping=mapping,
        source_identity=source_identity,
        reference_identity=reference_identity,
    )
    result.update({
        "reference_context_id": str(reference_context.get("id") or ""),
        "association_source": str(mapping.get("association_source") or "user_selected"),
    })

    params = {
        "candidate_ids": [str(candidates[0].get("candidate_id") or "")],
        "selected_structure_path": str(selected_path),
        "selected_model_id": source_identity["model_id"],
        "source_component_chain": source_chain,
        "reference_context_id": result["reference_context_id"],
        "reference_molecule_id": reference_identity["molecule_id"],
        "reference_structure_id": reference_identity["structure_id"],
        "reference_structure_digest": reference_identity["structure_digest"],
        "reference_model_id": reference_identity["model_id"],
        "reference_component_chain": reference_identity["chain_id"],
        "method": PLACEMENT_METHOD,
        "method_version": PLACEMENT_VERSION,
    }
    job = _analysis_job(
        source_run_dir, candidates_jsonl,
        tool="reference_placement", params=params,
        note="Fit the explicitly associated source component to its unique exact reference-chain correspondence using mapped C-alpha atoms. This is a coordinate-frame transform only; no biological accessibility or context score is calculated.",
        mode="saved-structure-coordinate-fit",
    )
    _write_complex_analysis_selection(
        job, source_run_dir, candidates[0], workflow="analysis.reference_placement",
        parameters={key: value for key, value in params.items() if key not in {"candidate_ids", "selected_structure_path"}},
    )
    placement_dir = job.run_dir / "artifacts" / "analysis" / "reference_placement"
    placement_dir.mkdir(parents=True, exist_ok=True)
    write_json(placement_dir / "reference_placement.json", result)
    candidate = copy.deepcopy(candidates[0])
    metrics = dict(candidate.get("metrics") or {})
    metrics.update({
        "reference_placement_aligned_ca_pair_count": result["aligned_pair_count"],
        "reference_placement_source_coverage": result["source_coverage"],
        "reference_placement_reference_coverage": result["reference_coverage"],
        "reference_placement_ca_rmsd_angstrom": result["rmsd_angstrom"],
        "reference_placement_mapping_method": result["mapping_method"],
        "reference_placement_status": "measured",
    })
    candidate["metrics"] = metrics
    raw = candidate.setdefault("raw_metadata", {})
    raw["reference_placement"] = result
    raw["measurement_provenance"] = {
        **dict(raw.get("measurement_provenance") or {}),
        "reference_placement": {
            "method": result["method"],
            "method_version": result["method_version"],
            "parameters": result["parameters"],
            "reference_context_id": result["reference_context_id"],
            "source": result["source"],
            "reference": result["reference"],
            "result_artifact": "artifacts/analysis/reference_placement/reference_placement.json",
            "biological_accessibility_claimed": False,
        },
    }
    raw["analysis_parameters"] = params
    return _finish_metric_analysis(
        job, [candidate], tool="reference_placement",
        table_name="reference_placement_metrics.csv",
        metric_names=[
            "reference_placement_aligned_ca_pair_count",
            "reference_placement_source_coverage",
            "reference_placement_reference_coverage",
            "reference_placement_ca_rmsd_angstrom",
            "reference_placement_mapping_method",
        ],
        required_metrics=["reference_placement_ca_rmsd_angstrom"],
        error_key="reference_placement_error",
    )


def run_reference_environment_contract(
    source_run_dir: Path, candidates_jsonl: Path,
    candidate_ids: list[str] | None = None, selected_structure_path: Path | None = None,
    selected_model_id: str | None = None, source_component_chain: str | None = None,
    selected_reference_chains: list[str] | None = None, assembly_id: str | None = None,
) -> Path:
    """Measure static atom proximity to omitted and explicitly selected reference context."""
    from mn_molecular_design_studio.complex_analysis_environment import (
        METHOD, METHOD_VERSION, calculate_reference_environment,
    )
    from mn_molecular_design_studio.evaluations import file_digest

    source_run_dir = Path(source_run_dir).expanduser().resolve(strict=True)
    candidates_jsonl = Path(candidates_jsonl).expanduser().resolve(strict=True)
    candidates = _analysis_source_candidates(source_run_dir, candidates_jsonl)
    candidates = _select_complex_analysis_candidates(
        source_run_dir, candidates, candidate_ids=candidate_ids,
        selected_structure_path=selected_structure_path, selected_model_id=selected_model_id,
    )
    if len(candidates) != 1:
        raise ValueError("Reference-environment evaluation requires exactly one selected Database complex.")
    metadata = read_json(source_run_dir / "metadata.json")
    manifest = read_json(source_run_dir / "complex_analysis_source.json")
    context = manifest.get("reference_environment_context")
    if metadata.get("source_kind") != "database" or manifest.get("schema") != "studio.complex_analysis_database_source.v1":
        raise ValueError("Reference-environment evaluation requires an authorized Database source bundle.")
    if not isinstance(context, dict) or context.get("status") != "validated":
        raise ValueError("An explicit, current, validated Database reference association is required.")
    source_chain = str(source_component_chain or "")
    if source_chain != str(context.get("source_component_chain") or ""):
        raise ValueError("The selected associated source component does not match the saved reference context.")
    source_path = Path(str(candidates[0].get("complex_pdb") or "")).resolve(strict=True)
    reference_path = Path(str(context.get("structure_path") or "")).resolve(strict=True)
    for path in (source_path, reference_path):
        try:
            path.relative_to(source_run_dir)
        except ValueError as exc:
            raise ValueError("Source and reference coordinates must be staged inside the authorized Database bundle.") from exc
    source_identity = {
        "molecule_id": str(manifest.get("molecule_id") or ""),
        "structure_id": str(manifest.get("structure_id") or ""),
        "group_id": str(manifest.get("group_id") or ""),
        "structure_digest": str(manifest.get("structure_digest") or ""),
        "model_id": str(manifest.get("model_id") or ""), "chain_id": source_chain,
    }
    reference_identity = {
        "molecule_id": str(context.get("reference_molecule_id") or ""),
        "structure_id": str(context.get("reference_structure_id") or ""),
        "group_id": str(context.get("reference_group_id") or ""),
        "structure_digest": str(context.get("reference_structure_digest") or ""),
        "model_id": str(context.get("reference_model_id") or ""),
        "chain_id": str(context.get("reference_component_chain") or ""),
    }
    if not all(source_identity.values()) or not all(reference_identity.values()):
        raise ValueError("The saved association is missing exact source or reference identity fields.")
    if file_digest(source_path) != source_identity["structure_digest"] or file_digest(reference_path) != reference_identity["structure_digest"]:
        raise ValueError("The source or associated reference digest changed after authorization.")
    placement = context.get("placement_result") if isinstance(context.get("placement_result"), dict) else {}
    if (str(placement.get("reference_context_id") or "") != str(context.get("id") or "")
            or str(placement.get("evaluation_id") or "") != str(context.get("placement_evaluation_id") or "")
            or str(context.get("placement_run_id") or "") == ""):
        raise ValueError("A completed and ingested Reference Placement Evaluation matching this exact association is required.")
    result = calculate_reference_environment(
        source_path, reference_path, placement=placement,
        source_identity=source_identity, reference_identity=reference_identity,
        reference_context_id=str(context.get("id") or ""),
        selected_reference_chains=list(selected_reference_chains or []),
        assembly_id=str(assembly_id or "").strip() or None,
    )
    result["placement"]["evaluation_id"] = str(context.get("placement_evaluation_id") or "")
    result["placement"]["run_id"] = str(context.get("placement_run_id") or "")
    params = {
        "candidate_ids": [str(candidates[0].get("candidate_id") or "")],
        "selected_structure_path": str(source_path), "selected_model_id": source_identity["model_id"],
        "source_component_chain": source_chain, "reference_context_id": str(context.get("id") or ""),
        "reference_molecule_id": reference_identity["molecule_id"],
        "reference_structure_id": reference_identity["structure_id"],
        "reference_structure_digest": reference_identity["structure_digest"],
        "reference_model_id": reference_identity["model_id"],
        "reference_component_chain": reference_identity["chain_id"],
        "placement_evaluation_id": result["placement"]["evaluation_id"],
        "placement_run_id": result["placement"]["run_id"],
        "placement_mapping_method": result["placement"]["mapping_method"],
        "selected_reference_chains": list(selected_reference_chains or []),
        "assembly_id": str(assembly_id or "").strip() or None,
        "method": METHOD, "method_version": METHOD_VERSION,
    }
    job = _analysis_job(
        source_run_dir, candidates_jsonl, tool="reference_environment", params=params,
        note="Measure static atom-pair proximity and van-der-Waals envelope overlap against resolved omitted reference regions and explicitly selected reference components. No biological accessibility score is calculated.",
        mode="saved-structure-reference-context-geometry",
    )
    _write_complex_analysis_selection(
        job, source_run_dir, candidates[0], workflow="analysis.reference_environment",
        parameters={key: value for key, value in params.items() if key not in {"candidate_ids", "selected_structure_path"}},
    )
    output_dir = job.run_dir / "artifacts" / "analysis" / "reference_environment"
    output_dir.mkdir(parents=True, exist_ok=True)
    result_path = output_dir / "reference_environment.json"
    write_json(result_path, result)
    candidate = copy.deepcopy(candidates[0])
    metrics = dict(candidate.get("metrics") or {})
    measurements = result["measurements"]
    metrics.update({
        "reference_environment_context_residue_count": measurements["context_residue_count"],
        "reference_environment_residues_with_pairs_within_8A": measurements["residue_count_with_atom_pairs_within_proximity_cutoff"],
        "reference_environment_residues_with_vdw_overlap": measurements["residue_count_with_vdw_overlap_above_epsilon"],
        "reference_environment_closest_distance_angstrom": measurements["closest_distance_angstrom"],
        "reference_environment_source_component_closest_distance_angstrom": measurements["associated_source_component_closest_distance_angstrom"],
        "reference_environment_other_component_closest_distance_angstrom": measurements["other_source_complex_closest_distance_angstrom"],
        "reference_environment_atom_pairs_within_8A": measurements["atom_pairs_within_proximity_cutoff"],
        "reference_environment_atom_pairs_with_vdw_overlap": measurements["atom_pairs_with_vdw_overlap_above_epsilon"],
        "reference_environment_other_component_vdw_overlap_pairs": measurements["other_source_complex_atom_pairs_with_vdw_overlap_above_epsilon"],
        "reference_environment_status": "measured",
    })
    candidate["metrics"] = metrics
    raw = candidate.setdefault("raw_metadata", {})
    raw["reference_environment"] = result
    raw["analysis_parameters"] = params
    raw["measurement_provenance"] = {
        **dict(raw.get("measurement_provenance") or {}),
        "reference_environment": {
            "method": METHOD, "method_version": METHOD_VERSION,
            "parameters": result["parameters"], "source": source_identity,
            "reference": reference_identity, "reference_context_id": result["reference_context_id"],
            "placement": result["placement"],
            "result_artifact": "artifacts/analysis/reference_environment/reference_environment.json",
            "biological_accessibility_claimed": False,
        },
    }
    return _finish_metric_analysis(
        job, [candidate], tool="reference_environment",
        table_name="reference_environment_metrics.csv",
        metric_names=[
            "reference_environment_context_residue_count",
            "reference_environment_residues_with_pairs_within_8A",
            "reference_environment_residues_with_vdw_overlap",
            "reference_environment_closest_distance_angstrom",
            "reference_environment_source_component_closest_distance_angstrom",
            "reference_environment_other_component_closest_distance_angstrom",
            "reference_environment_atom_pairs_within_8A",
            "reference_environment_atom_pairs_with_vdw_overlap",
            "reference_environment_other_component_vdw_overlap_pairs",
        ],
        required_metrics=["reference_environment_closest_distance_angstrom"],
        error_key="reference_environment_error",
    )


def _run_tag_fusion_contract(
    source_run_dir: Path,
    candidates_jsonl: Path,
    *,
    contextual: bool,
    candidate_ids: list[str] | None = None,
    selected_structure_path: Path | None = None,
    selected_model_id: str | None = None,
    binder_chains: list[str] | None = None,
    target_chains: list[str] | None = None,
    source_component_chain: str | None = None,
    interface_contact_cutoff: float = 4.5,
    termini_confirmed: bool = False,
    component_format: str = "single_chain_linear_protein",
    b_factor_interpretation: str = "uninterpreted",
) -> Path:
    from mn_molecular_design_studio.complex_analysis_tag_fusion import calculate_tag_fusion_suitability
    from mn_molecular_design_studio.evaluations import file_digest

    source_run_dir = Path(source_run_dir).expanduser().resolve(strict=True)
    candidates_jsonl = Path(candidates_jsonl).expanduser().resolve(strict=True)
    source_candidates = _analysis_source_candidates(source_run_dir, candidates_jsonl)
    selected = _select_complex_analysis_candidates(
        source_run_dir, source_candidates, candidate_ids=candidate_ids,
        selected_structure_path=selected_structure_path, selected_model_id=selected_model_id,
        binder_chains=binder_chains, target_chains=target_chains,
    )
    if len(selected) != 1:
        raise ValueError("Tag and Fusion Suitability requires one selected Database structure.")
    source_metadata = read_json(source_run_dir / "metadata.json")
    if str(source_metadata.get("source_kind") or "") != "database":
        raise ValueError("Tag and Fusion Suitability is currently available only for authorized Database structures.")
    manifest = read_json(source_run_dir / "complex_analysis_source.json")
    if manifest.get("schema") != "studio.complex_analysis_database_source.v1":
        raise ValueError("The authorized Database source manifest is missing or unsupported.")
    selected_path = Path(str(selected[0].get("complex_pdb") or "")).resolve(strict=True)
    try:
        selected_path.relative_to(source_run_dir)
    except ValueError as exc:
        raise ValueError("The selected complex must remain inside the authorized Database source bundle.") from exc
    if file_digest(selected_path) != str(manifest.get("structure_digest") or ""):
        raise ValueError("The selected source structure changed after Database authorization.")
    if str(selected_model_id or "") != str(manifest.get("model_id") or ""):
        raise ValueError("The selected model does not match the authorized Database representation.")
    component_group = list(binder_chains or [])
    partner_group = list(target_chains or [])
    if len(component_group) != 1 or not partner_group:
        raise ValueError("Choose one protein component and one or more disjoint interface-partner chains.")
    if set(component_group) & set(partner_group):
        raise ValueError("The assessed protein component and interface-partner chains must be disjoint.")

    reference_context = None
    contextual_source = manifest.get("tag_fusion_context") if contextual else None
    if contextual:
        if not isinstance(contextual_source, dict) or contextual_source.get("status") != "validated":
            raise ValueError("Contextual Tag and Fusion analysis requires a staged, explicitly validated reference association.")
        if (not str(source_component_chain or "")
                or str(contextual_source.get("source_component_chain") or "") != str(source_component_chain)):
            raise ValueError("The saved placement Evaluation must belong to the explicitly selected reference-associated source chain.")
        ref_path = Path(str(contextual_source.get("structure_path") or "")).resolve(strict=True)
        try:
            ref_path.relative_to(source_run_dir)
        except ValueError as exc:
            raise ValueError("The associated reference must remain inside the authorized Database source bundle.") from exc
        if file_digest(ref_path) != str(contextual_source.get("reference_structure_digest") or ""):
            raise ValueError("The associated reference structure changed after Database authorization.")
        reference_context = {**contextual_source, "structure_path": str(ref_path)}
    elif source_component_chain:
        raise ValueError("Reference identity is accepted only by the separately registered contextual evaluation.")

    source_identity = {
        key: str(manifest.get(key) or "") for key in (
            "molecule_id", "structure_id", "group_id", "structure_digest", "model_id", "evaluation_id",
        )
    }
    if any(not source_identity[key] for key in ("molecule_id", "structure_id", "group_id", "structure_digest", "model_id")):
        raise ValueError("The Database source manifest is missing exact structure identity provenance.")
    source_identity["structure_digest"] = file_digest(selected_path)
    result = calculate_tag_fusion_suitability(
        selected_path, component_chain=component_group[0], interface_chains=partner_group,
        model_id=str(selected_model_id), source_identity=source_identity,
        interface_contact_cutoff=float(interface_contact_cutoff), termini_confirmed=bool(termini_confirmed),
        component_format=str(component_format), b_factor_interpretation=str(b_factor_interpretation),
        saved_confidence_context=(manifest.get("selected_evaluation_confidence_context")
                                  if isinstance(manifest.get("selected_evaluation_confidence_context"), dict) else None),
        reference_context=reference_context,
    )
    workflow = "analysis.tag_fusion_context" if contextual else "analysis.tag_fusion_suitability"
    tool = "tag_fusion_context" if contextual else "tag_fusion_suitability"
    params = {
        "candidate_ids": [str(selected[0].get("candidate_id") or "")],
        "selected_structure_path": str(selected_path), "selected_model_id": str(selected_model_id),
        "binder_chains": component_group, "target_chains": partner_group,
        "interface_contact_cutoff": float(interface_contact_cutoff),
        "termini_confirmed": bool(termini_confirmed), "component_format": str(component_format),
        "b_factor_interpretation": str(b_factor_interpretation),
    }
    if contextual:
        params["source_component_chain"] = str(source_component_chain or "")
    job = _analysis_job(
        source_run_dir, candidates_jsonl, tool=tool, params=params,
        note="Evaluate saved N/C terminus geometry, accessibility evidence, and advisory tag/fusion clearance. This calculation does not simulate a tag or establish experimental accessibility.",
        mode="saved-database-structure-geometry",
    )
    _write_complex_analysis_selection(
        job, source_run_dir, selected[0], workflow=workflow,
        parameters={key: value for key, value in params.items() if key not in {"candidate_ids", "selected_structure_path"}},
    )
    analysis_dir = job.run_dir / "artifacts" / "analysis"
    analysis_dir.mkdir(parents=True, exist_ok=True)
    result_path = analysis_dir / "tag_fusion_suitability.json"
    write_json(result_path, result)
    candidate = copy.deepcopy(selected[0])
    metrics = dict(candidate.get("metrics") or {})
    metrics["tag_fusion_assessment_status"] = str(result.get("assessment_status") or "unknown")
    metrics["tag_fusion_assessment_scope"] = str(result.get("assessment_scope") or "intrinsic_complex")
    metrics["tag_fusion_component_chain"] = component_group[0]
    metrics["tag_fusion_partner_chains"] = partner_group
    if isinstance(result.get("interface"), dict):
        metrics["tag_fusion_interface_contact_residue_count"] = result["interface"].get("component_residue_count_contacting_partner")
    for end, prefix in (("N", "n_terminus"), ("C", "c_terminus")):
        terminal = (result.get("terminal_geometry") or {}).get(end) or {}
        metrics[f"tag_fusion_{prefix}_status"] = str(terminal.get("assessment_status") or "not_assessed")
        for field, suffix in (("distance_to_selected_interface_angstrom", "distance_to_partner_angstrom"),
                              ("terminal_residue_sasa_angstrom2", "residue_sasa_angstrom2"),
                              ("terminal_window_mean_plddt", "window_mean_plddt")):
            value = terminal.get(field)
            if value is not None:
                metrics[f"tag_fusion_{prefix}_{suffix}"] = value
    analysis_status = str(result.get("assessment_status") or "not_assessed")
    metrics["analysis_status"] = "measured" if analysis_status == "measured" else "partial"
    candidate["metrics"] = metrics
    raw = candidate.setdefault("raw_metadata", {})
    raw["tag_fusion_suitability"] = result
    raw["analysis_parameters"] = params
    raw["measurement_status"] = {
        "status": "measured" if analysis_status == "measured" else "partial",
        "assessment_status": analysis_status,
        "error": "",
    }
    raw["measurement_provenance"] = {
        **dict(raw.get("measurement_provenance") or {}),
        tool: {
            "method": result.get("method"), "method_version": result.get("method_version"),
            "parameters": result.get("parameters"), "source": result.get("source"),
            "result_artifact": str(result_path.relative_to(job.run_dir)),
            "attribution": result.get("adaptation_attribution"),
            "assessment_scope": result.get("assessment_scope"),
            "reference_context": result.get("reference_context") or {},
            "experimental_claims": False,
        },
    }
    table_path = analysis_dir / "tag_fusion_suitability_metrics.csv"
    _analysis_metrics_table([candidate]).to_csv(table_path, index=False)
    outputs = _analysis_output_candidates([candidate], tool=tool)
    normalized = write_candidates(job.run_dir, tool, outputs)
    scientific_status = "valid" if analysis_status == "measured" else "partial"
    summary = {
        "tool": tool, "evaluation_count": 1, "execution_status": "completed",
        "scientific_status": scientific_status, "assessment_status": analysis_status,
        "assessment_scope": result.get("assessment_scope"), "result_artifact": str(result_path.relative_to(job.run_dir)),
        "source_ranking_preserved": True,
    }
    write_json(analysis_dir / "tag_fusion_suitability_summary.json", summary)
    artifacts = _candidate_artifacts(job.run_dir)
    finish_job(job.run_dir, True, {
        "outputs": {"artifacts": artifacts, "candidates": normalized},
        "metrics": {"evaluation_count": 1, "scientific_status": scientific_status,
                    "assessment_status": analysis_status},
        "execution_status": "completed", "scientific_status": scientific_status,
        "downstream_artifacts": {
            "candidates_jsonl": "artifacts/normalized_candidates/candidates.jsonl",
            "campaign_result": "artifacts/normalized_candidates/campaign_result.json",
            "tag_fusion_suitability": str(result_path.relative_to(job.run_dir)),
        },
        "summary": summary,
    })
    return job.run_dir


def run_tag_fusion_suitability_contract(
    source_run_dir: Path, candidates_jsonl: Path,
    candidate_ids: list[str] | None = None, selected_structure_path: Path | None = None,
    selected_model_id: str | None = None, binder_chains: list[str] | None = None,
    target_chains: list[str] | None = None, interface_contact_cutoff: float = 4.5,
    termini_confirmed: bool = False, component_format: str = "single_chain_linear_protein",
    b_factor_interpretation: str = "uninterpreted",
) -> Path:
    return _run_tag_fusion_contract(
        source_run_dir, candidates_jsonl, contextual=False, candidate_ids=candidate_ids,
        selected_structure_path=selected_structure_path, selected_model_id=selected_model_id,
        binder_chains=binder_chains, target_chains=target_chains,
        interface_contact_cutoff=interface_contact_cutoff, termini_confirmed=termini_confirmed,
        component_format=component_format, b_factor_interpretation=b_factor_interpretation,
    )


def run_tag_fusion_context_contract(
    source_run_dir: Path, candidates_jsonl: Path,
    candidate_ids: list[str] | None = None, selected_structure_path: Path | None = None,
    selected_model_id: str | None = None, binder_chains: list[str] | None = None,
    target_chains: list[str] | None = None, source_component_chain: str | None = None,
    interface_contact_cutoff: float = 4.5, termini_confirmed: bool = False,
    component_format: str = "single_chain_linear_protein",
    b_factor_interpretation: str = "uninterpreted",
) -> Path:
    return _run_tag_fusion_contract(
        source_run_dir, candidates_jsonl, contextual=True, candidate_ids=candidate_ids,
        selected_structure_path=selected_structure_path, selected_model_id=selected_model_id,
        binder_chains=binder_chains, target_chains=target_chains,
        source_component_chain=source_component_chain,
        interface_contact_cutoff=interface_contact_cutoff, termini_confirmed=termini_confirmed,
        component_format=component_format, b_factor_interpretation=b_factor_interpretation,
    )


def _safe_float(value: Any) -> float | None:
    return finite_number(value)


def _allocated_cpu_cores(requested: int | None = None) -> int:
    """Use the CPU slots granted by the shared scheduler for batch parallelism."""
    if requested not in (None, ""):
        cores = int(requested)
    else:
        try:
            cores = int(os.getenv("MN_MOLECULAR_DESIGN_STUDIO_JOB_CPU_CORES", "1"))
        except (TypeError, ValueError):
            cores = 1
    if cores < 1:
        raise ValueError("CPU allocation must be at least one core.")
    granted = os.getenv("MN_MOLECULAR_DESIGN_STUDIO_JOB_CPU_CORES")
    if granted:
        try:
            cores = min(cores, max(1, int(granted)))
        except ValueError:
            pass
    return cores


def _resolve_path(base_dir: Path, path_text: str | None) -> Path | None:
    if not path_text:
        return None
    return resolve_stored_path(path_text, run_dir=base_dir)


def _safe_id(candidate: dict[str, Any]) -> str:
    return str(candidate.get("candidate_id") or "candidate").replace("/", "_").replace(" ", "_")


def _candidate_design_reference(source_run_dir: Path, candidate: dict[str, Any]) -> tuple[Path | None, list[str]]:
    """Find a separate designed binder structure without confusing it with the cofold."""
    raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
    pairing = raw.get("analysis_design_pairing") if isinstance(raw.get("analysis_design_pairing"), dict) else {}
    if pairing.get("status") in {"unmatched", "ambiguous"}:
        return None, []
    predicted = _resolve_path(source_run_dir, candidate.get("complex_pdb"))
    pending: list[tuple[dict[str, Any], Path, bool]] = [(candidate, source_run_dir, False)]
    seen: set[int] = set()
    while pending:
        item, root, inherited = pending.pop(0)
        if id(item) in seen:
            continue
        seen.add(id(item))
        raw = item.get("raw_metadata") if isinstance(item.get("raw_metadata"), dict) else {}
        nested_root_text = raw.get("upstream_source_run_dir") or raw.get("source_run_dir") or raw.get("design_run_dir")
        nested_root = Path(nested_root_text) if nested_root_text else root
        keys = ["design_structure_path", "reference_structure_path", "designed_structure_path", "input_structure_path"]
        if inherited or str(item.get("stage") or "").startswith(("generation.", "sequence_design")):
            keys.extend(["binder_pdb", "complex_pdb"])
        for key in keys:
            value = item.get(key) or raw.get(key)
            path = _resolve_path(nested_root, value)
            if not path or not path.is_file():
                continue
            if predicted and path.resolve() == predicted.resolve():
                continue
            binder_chains = [str(value) for value in item.get("binder_chains") or [] if str(value)]
            if not binder_chains:
                binder_chains = [str(value) for value in (raw.get("binder_chains") or []) if str(value)]
            return path, binder_chains
        nested = raw.get("source_candidate")
        if isinstance(nested, dict):
            pending.append((nested, nested_root, True))
        evidence = raw.get("evidence_candidates")
        if isinstance(evidence, list):
            pending.extend((entry, nested_root, True) for entry in evidence if isinstance(entry, dict))
    return None, []


def _candidate_design_complex(source_run_dir: Path, candidate: dict[str, Any]) -> tuple[Path | None, list[str], list[str]]:
    """Locate the upstream generated complex for contact-pattern comparison."""
    raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
    pairing = raw.get("analysis_design_pairing") if isinstance(raw.get("analysis_design_pairing"), dict) else {}
    if pairing.get("status") in {"unmatched", "ambiguous"}:
        return None, [], []
    predicted = _resolve_path(source_run_dir, candidate.get("complex_pdb"))
    pending: list[tuple[dict[str, Any], Path, bool]] = [(candidate, source_run_dir, False)]
    seen: set[int] = set()
    while pending:
        item, root, inherited = pending.pop(0)
        if id(item) in seen:
            continue
        seen.add(id(item))
        raw = item.get("raw_metadata") if isinstance(item.get("raw_metadata"), dict) else {}
        nested_root_text = raw.get("upstream_source_run_dir") or raw.get("source_run_dir") or raw.get("design_run_dir")
        nested_root = Path(nested_root_text) if nested_root_text else root
        keys = ["design_complex_path", "reference_complex_path"]
        if inherited or str(item.get("stage") or "").startswith(("generation.", "sequence_design")):
            keys.append("complex_pdb")
        for key in keys:
            path = _resolve_path(nested_root, item.get(key) or raw.get(key))
            if not path or not path.is_file() or (predicted and path.resolve() == predicted.resolve()):
                continue
            binder_chains = [str(value) for value in item.get("binder_chains") or raw.get("binder_chains") or [] if str(value)]
            target_chains = [str(value) for value in item.get("target_chains") or raw.get("target_chains") or [] if str(value)]
            return path, binder_chains, target_chains
        nested = raw.get("source_candidate")
        if isinstance(nested, dict):
            pending.append((nested, nested_root, True))
        evidence = raw.get("evidence_candidates")
        if isinstance(evidence, list):
            pending.extend((entry, nested_root, True) for entry in evidence if isinstance(entry, dict))
    return None, [], []


def _candidate_hotspots(candidate: dict[str, Any]) -> str:
    pending = [candidate]
    seen: set[int] = set()
    while pending:
        item = pending.pop(0)
        if id(item) in seen:
            continue
        seen.add(id(item))
        value = item.get("hotspots")
        if value:
            return ",".join(map(str, value)) if isinstance(value, (list, tuple, set)) else str(value)
        raw = item.get("raw_metadata") if isinstance(item.get("raw_metadata"), dict) else {}
        for key in ("hotspots", "hotspot_residues", "target_hotspots"):
            value = raw.get(key)
            if value:
                return ",".join(map(str, value)) if isinstance(value, (list, tuple, set)) else str(value)
        nested = raw.get("source_candidate")
        if isinstance(nested, dict):
            pending.append(nested)
    return ""


def _flatten_metrics(metrics: dict[str, Any], prefix: str = "") -> dict[str, Any]:
    flat: dict[str, Any] = {}
    for key, value in metrics.items():
        flat_key = f"{prefix}{key}" if not prefix else f"{prefix}.{key}"
        if isinstance(value, dict):
            flat.update(_flatten_metrics(value, flat_key))
        elif isinstance(value, (str, int, float, bool)) or value is None:
            flat[flat_key] = value
    return flat


def _find_pae_path(source_run_dir: Path, candidate: dict[str, Any], structure_path: Path | None) -> Path | None:
    search_dirs: list[Path] = []
    raw = candidate.get("raw_metadata") or {}
    explicit_pae = raw.get("pae_path") or candidate.get("pae_path")
    if explicit_pae:
        resolved = _resolve_path(source_run_dir, str(explicit_pae))
        if resolved and resolved.exists():
            return resolved
    if structure_path:
        direct_matches = [
            structure_path.with_name(f"{structure_path.stem}_pae.json"),
            structure_path.with_name(f"{structure_path.stem}.pae.json"),
            structure_path.with_name(f"pae_{structure_path.stem}.json"),
            structure_path.with_suffix(".npz"),
        ]
        for match in direct_matches:
            if match.exists():
                return match
    prediction_dir = raw.get("prediction_dir")
    if prediction_dir:
        resolved = _resolve_path(source_run_dir, str(prediction_dir))
        if resolved:
            search_dirs.append(resolved)
    if structure_path:
        search_dirs.append(structure_path.parent)
    search_dirs.append(source_run_dir)

    patterns = [
        "pae*.npz",
        "*pae*.npz",
        "pae*.json",
        "*pae*.json",
        "scores*.json",
        "*scores*.json",
    ]
    seen: set[Path] = set()
    for directory in search_dirs:
        if not directory or not directory.exists() or directory in seen:
            continue
        seen.add(directory)
        if structure_path:
            stem = structure_path.stem
            for pattern in [f"{stem}*pae*.json", f"{stem}*pae*.npz", f"*{stem}*pae*.json", f"*{stem}*pae*.npz"]:
                matches = sorted(directory.glob(pattern))
                if matches:
                    return matches[0]
        for pattern in patterns:
            matches = sorted(directory.glob(pattern))
            if matches:
                return matches[0]
    return None


def _ipsae_engine(candidate: dict[str, Any], pae_path: Path | None) -> tuple[str, str]:
    source_tool = str(candidate.get("source_tool") or candidate.get("tool") or "").lower()
    if pae_path and pae_path.suffix.lower() == ".npz":
        return "boltz", "boltz_npz"
    if "boltz" in source_tool:
        return "boltz", "boltz_npz"
    return "af2", "af2_json"


def _parse_ipsae_row(path: Path) -> dict[str, Any]:
    header = [
        "ID",
        "ipsae",
        "ipsae_d0chn",
        "ipsae_d0dom",
        "iptm_af",
        "iptm_d0chn",
        "binder_chain",
        "target_chain",
        "type",
        "engine",
        "pae_format",
        "error",
        "ipae_binder_to_target",
        "ipae_target_to_binder",
        "ipae_bidirectional",
        "ipae_binder_to_target_norm31",
        "ipae_target_to_binder_norm31",
        "ipae_bidirectional_norm31",
        "ipae_contact_pair_count",
        "ipae_contact_cutoff",
        "bindcraft2_ipdae_recomputed",
        "bindcraft2_ipdae_binder_anchor_max",
        "bindcraft2_ipdae_target_anchor_max",
        "bindcraft2_ipdae_contact_pair_count",
        "bindcraft2_ipdae_contact_cutoff",
        "pae_measurement_error",
        "structure_interface_atom_contact_pair_count",
        "structure_interface_residue_contact_pair_count",
        "structure_binder_interface_residue_count",
        "structure_target_interface_residue_count",
        "structure_ca_clash_pair_count_2_5a",
        "structure_binder_ca_bfactor_mean",
        "structure_target_ca_bfactor_mean",
        "structure_complex_ca_bfactor_mean",
        "structure_contact_cutoff",
        "structure_hotspot_count",
        "structure_hotspot_found_count",
        "structure_hotspot_contacted_count",
        "structure_hotspot_contact_fraction_found",
        "structure_hotspot_contact_fraction_requested",
        "structure_hotspot_min_distance",
        "structure_hotspot_mean_min_distance",
        "structure_hotspot_contact_cutoff",
        "structure_hotspots_contacted",
        "structure_hotspots_missing",
        "structure_metrics_status",
        "structure_metrics_error",
        "design_vs_cofold_rmsd_status",
        "design_vs_cofold_rmsd_error",
        "structure_binder_design_to_cofold_fitted_ca_rmsd",
        "structure_binder_rmsd_ca_count",
        "structure_binder_rmsd_alignment_coverage",
        "target_aligned_pose_rmsd_status",
        "target_aligned_pose_rmsd_error",
        "structure_target_aligned_binder_pose_ca_rmsd",
        "structure_target_alignment_ca_rmsd",
        "structure_target_alignment_ca_count",
        "structure_target_aligned_binder_ca_count",
        "structure_target_aligned_binder_sequence_coverage",
        "hotspot_analysis_status",
        "hotspot_analysis_error",
        "interaction_preservation_status",
        "interaction_preservation_error",
        "structure_design_interface_residue_contact_pairs",
        "structure_cofold_interface_residue_contact_pairs",
        "structure_interface_contact_pairs_preserved",
        "structure_design_interaction_preservation_fraction",
        "structure_interface_contact_jaccard",
        "structure_binder_interface_residue_preservation_fraction",
        "structure_target_interface_residue_preservation_fraction",
        "structure_interaction_sequence_mapping_coverage",
    ]
    if not path.exists():
        return {"ipsae_error": "IPSAE output row missing"}
    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        rows = list(csv.reader(handle))
    if not rows:
        return {"ipsae_error": "IPSAE output row empty"}
    row = dict(zip(header, rows[0]))
    parsed: dict[str, Any] = {}
    for key, value in row.items():
        if key == "ID":
            continue
        if key == "error":
            if value:
                parsed["ipsae_error"] = value
            continue
        if key == "pae_measurement_error":
            if value:
                parsed[key] = value
            continue
        if key in {"binder_chain", "target_chain", "type", "engine", "pae_format", "structure_metrics_status",
                   "structure_metrics_error", "design_vs_cofold_rmsd_status", "design_vs_cofold_rmsd_error",
                   "target_aligned_pose_rmsd_status", "target_aligned_pose_rmsd_error",
                   "structure_hotspots_contacted", "structure_hotspots_missing", "hotspot_analysis_status", "hotspot_analysis_error",
                   "interaction_preservation_status", "interaction_preservation_error"}:
            parsed[key] = value
            continue
        number = _safe_float(value)
        parsed[key] = number if number is not None else value
    return parsed


def _find_ipsae_summary(structure: Path, pae_cutoff: float, dist_cutoff: float) -> Path | None:
    pae_text = f"{int(pae_cutoff):02d}"
    dist_text = f"{int(dist_cutoff):02d}"
    candidates = [
        structure.parent / f"{structure.stem}_{pae_text}_{dist_text}.txt",
        structure.parent / f"{structure.with_suffix('').name}_{pae_text}_{dist_text}.txt",
    ]
    for path in candidates:
        if path.exists():
            return path
    matches = sorted(
        [
            path
            for path in structure.parent.glob(f"{structure.stem}_*_*.txt")
            if not path.name.endswith("_byres.txt") and not path.name.endswith("_done.txt")
        ],
        key=lambda path: path.stat().st_mtime,
        reverse=True,
    )
    return matches[0] if matches else None


def _float_field(row: dict[str, str], key: str) -> float | None:
    return _safe_float(row.get(key))


def _parse_ipsae_summary_scores(summary_path: Path, focus_chain: str = "A") -> dict[str, Any]:
    with summary_path.open("r", encoding="utf-8", errors="replace") as handle:
        lines = [line.strip() for line in handle if line.strip() and not line.startswith("#")]
    if not lines:
        return {"ipsae_error": "IPSAE summary file empty"}

    header_index = None
    headers: list[str] = []
    for index, line in enumerate(lines):
        if line.startswith("Chn1") and "Chn2" in line and "ipSAE" in line:
            header_index = index
            headers = line.split()
            break
    if header_index is None:
        return {"ipsae_error": f"could not parse IPSAE summary: {summary_path.name}"}

    rows: list[dict[str, str]] = []
    for line in lines[header_index + 1 :]:
        parts = line.split()
        if len(parts) < len(headers):
            continue
        row = dict(zip(headers, parts))
        if focus_chain in {row.get("Chn1"), row.get("Chn2")}:
            rows.append(row)
    if not rows:
        return {"ipsae_error": f"no IPSAE chain pairs found for chain {focus_chain}"}

    data_max: dict[str, dict[str, list[float]]] = {}
    data_asym: dict[str, list[float]] = {}
    ipae_by_partner: dict[str, list[float]] = {}
    pdockq_by_partner: dict[str, dict[str, list[float]]] = {}
    qualified_contact_rows = 0

    for row in rows:
        chain1 = row.get("Chn1")
        chain2 = row.get("Chn2")
        if not chain1 or not chain2 or focus_chain not in {chain1, chain2}:
            continue
        partner = chain2 if chain1 == focus_chain else chain1
        row_type = (row.get("Type") or "").lower()

        if row_type == "max":
            ipae = _float_field(row, "ipae")
            if ipae is not None:
                ipae_by_partner.setdefault(partner, []).append(ipae)

        dist1 = _float_field(row, "dist1")
        dist2 = _float_field(row, "dist2")
        has_qualified_contacts = dist1 not in {None, 0.0} and dist2 not in {None, 0.0}
        if has_qualified_contacts:
            qualified_contact_rows += 1

        p_bucket = pdockq_by_partner.setdefault(partner, {"pDockQ": [], "pDockQ2": []})
        for source_key in p_bucket:
            value = _float_field(row, source_key)
            if value is not None:
                p_bucket[source_key].append(value)

        if row_type == "max":
            bucket = data_max.setdefault(
                partner,
                {
                    "ipSAE": [],
                    "ipSAE_avg": [],
                    "LIS": [],
                    "ipSAE_min_in_calculation": [],
                    "ipSAE_d0chn": [],
                    "ipSAE_d0dom": [],
                },
            )
            for source_key in bucket:
                value = _float_field(row, source_key)
                if value is not None:
                    bucket[source_key].append(value)
        elif row_type == "asym":
            value = _float_field(row, "ipSAE")
            if value is not None:
                data_asym.setdefault(partner, []).append(value)

    def avg(values: list[float]) -> float | None:
        return sum(values) / len(values) if values else None

    ipsae_min = avg([min(values) for values in data_asym.values() if values])
    ipsae_max = avg([max(values["ipSAE"]) for values in data_max.values() if values["ipSAE"]])
    ipsae_avg = avg([min(values["ipSAE_avg"]) for values in data_max.values() if values["ipSAE_avg"]])
    lis = avg([min(values["LIS"]) for values in data_max.values() if values["LIS"]])
    ipsae_min_calc = avg(
        [
            min(values["ipSAE_min_in_calculation"])
            for values in data_max.values()
            if values["ipSAE_min_in_calculation"]
        ]
    )
    ipsae_d0chn = avg([min(values["ipSAE_d0chn"]) for values in data_max.values() if values["ipSAE_d0chn"]])
    ipsae_d0dom = avg([min(values["ipSAE_d0dom"]) for values in data_max.values() if values["ipSAE_d0dom"]])
    ipae = avg([min(values) for values in ipae_by_partner.values() if values])

    pdockq_min = avg([min(values["pDockQ"]) for values in pdockq_by_partner.values() if values["pDockQ"]])
    pdockq_max = avg([max(values["pDockQ"]) for values in pdockq_by_partner.values() if values["pDockQ"]])
    pdockq2_min = avg([min(values["pDockQ2"]) for values in pdockq_by_partner.values() if values["pDockQ2"]])
    pdockq2_max = avg([max(values["pDockQ2"]) for values in pdockq_by_partner.values() if values["pDockQ2"]])

    result = {
        "ipsae": ipsae_min,
        "ipsae_min": ipsae_min,
        "ipsae_max": ipsae_max,
        "ipsae_avg": ipsae_avg,
        "lis": lis,
        "ipsae_min_in_calculation": ipsae_min_calc,
        "ipsae_d0chn": ipsae_d0chn,
        "ipsae_d0dom": ipsae_d0dom,
        "ipae": ipae,
        "pdockq_min": pdockq_min,
        "pdockq_max": pdockq_max,
        "pdockq2_min": pdockq2_min,
        "pdockq2_max": pdockq2_max,
        "ipsae_summary": str(summary_path),
    }
    if ipsae_min is None and ipsae_max is None:
        result["ipsae_error"] = "IPSAE summary has no parsable chain-pair scores"
    elif qualified_contact_rows == 0:
        result["ipsae_error"] = "no PAE-qualified interface contacts at current IPSAE cutoffs"
    return result


def _pae_matrix_from_json(path: Path) -> list[list[float]] | None:
    try:
        payload = json.loads(path.read_text())
    except Exception:
        return None
    matrix = payload.get("predicted_aligned_error") or payload.get("pae") or payload.get("token_pair_pae")
    if not isinstance(matrix, list) or not matrix:
        return None
    rows: list[list[float]] = []
    for row in matrix:
        if not isinstance(row, list):
            return None
        values = [_safe_float(value) for value in row]
        if any(value is None for value in values):
            return None
        rows.append([float(value) for value in values if value is not None])
    if any(len(row) != len(rows) for row in rows):
        return None
    return rows


def _chain_ca_coordinates(structure_path: Path) -> dict[str, list[tuple[float, float, float]]]:
    chains: dict[str, list[tuple[float, float, float]]] = {}
    seen: set[tuple[str, str, str]] = set()
    with structure_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")) or line[12:16].strip() != "CA":
                continue
            chain = line[21].strip()
            residue_key = (chain, line[22:26].strip(), line[26].strip())
            if residue_key in seen:
                continue
            seen.add(residue_key)
            try:
                coord = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                continue
            chains.setdefault(chain, []).append(coord)
    return chains


def _chain_roles(candidate: dict[str, Any]) -> tuple[str, str]:
    binder_chains, target_chains = _chain_role_groups(candidate)
    if not binder_chains or not target_chains:
        raise ValueError("Explicit binder and target chain assignments are required.")
    return (binder_chains[0], target_chains[0])


def _chain_role_groups(candidate: dict[str, Any]) -> tuple[list[str], list[str]]:
    binder_chains = [str(chain) for chain in candidate.get("binder_chains") or [] if str(chain)]
    target_chains = [str(chain) for chain in candidate.get("target_chains") or [] if str(chain)]
    raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
    role_payload = candidate.get("chain_roles")
    if not isinstance(role_payload, dict):
        role_payload = raw.get("chain_roles")
    role_payload = role_payload if isinstance(role_payload, dict) else {}
    roles = role_payload.get("roles")
    roles = roles if isinstance(roles, dict) else {}
    if not binder_chains:
        binder_chains = [
            str(chain) for chain in role_payload.get("binder_chains") or [] if str(chain)
        ] or [str(chain) for chain, role in roles.items() if role == "binder"]
    if not target_chains:
        target_chains = [
            str(chain) for chain in role_payload.get("target_chains") or [] if str(chain)
        ] or [str(chain) for chain, role in roles.items() if role == "target"]
    entity_target_chains = [
        str(chain)
        for chain in (
            candidate.get("biological_target_chains")
            or candidate.get("target_entity_chains")
            or raw.get("biological_target_chains")
            or raw.get("target_entity_chains")
            or []
        )
        if str(chain)
    ]
    if entity_target_chains and not target_chains and not bool(raw.get("target_refolding")):
        target_chains = entity_target_chains
    if target_chains:
        target_set = set(target_chains)
        binder_chains = [chain for chain in binder_chains if chain not in target_set]
    return binder_chains, target_chains


def _collapse_target_chains_for_ipsae(
    structure: Path,
    *,
    binder_chain: str,
    target_chains: list[str],
) -> Path:
    """Create a PDB scoring copy where target fragments form one chain group."""
    target_group = [
        chain for chain in dict.fromkeys(target_chains) if chain and chain != binder_chain
    ]
    if not structure.exists() or structure.suffix.lower() != ".pdb" or len(target_group) <= 1:
        return structure

    merged_chain = target_group[0]
    merged = structure.with_name(
        f"{structure.stem}_targetgroup_{merged_chain}{structure.suffix}"
    )
    lines: list[str] = []
    residue_numbers: dict[tuple[str, str, str, str], int] = {}
    next_residue = 1
    changed = False
    for line in structure.read_text(errors="ignore").splitlines(keepends=True):
        record = line[:6].strip()
        source_chain = line[21].strip() if len(line) > 21 else ""
        if record in {"ATOM", "HETATM"} and source_chain in target_group:
            residue_key = (
                source_chain,
                line[22:26].strip(),
                line[26:27],
                line[17:20],
            )
            if residue_key not in residue_numbers:
                residue_numbers[residue_key] = next_residue
                next_residue += 1
            line = (
                f"{line[:21]}{merged_chain[:1]}"
                f"{residue_numbers[residue_key]:4d}{line[26:]}"
            )
            changed = True
        elif record in {"TER", "ANISOU"} and source_chain in target_group:
            line = f"{line[:21]}{merged_chain[:1]}{line[22:]}"
            changed = True
        lines.append(line)
    if changed:
        merged.write_text("".join(lines), encoding="utf-8")
        return merged
    return structure


def _distance(left: tuple[float, float, float], right: tuple[float, float, float]) -> float:
    return math.sqrt(sum((left[index] - right[index]) ** 2 for index in range(3)))


def _interface_pae_scores(
    structure_path: Path | None,
    pae_path: Path | None,
    binder_chain: str = "A",
    target_chain: str = "B",
    target_chains: list[str] | None = None,
    dist_cutoff: float = 10.0,
) -> dict[str, Any]:
    if not structure_path or not structure_path.exists() or not pae_path or not pae_path.exists() or pae_path.suffix.lower() != ".json":
        return {}
    matrix = _pae_matrix_from_json(pae_path)
    if matrix is None:
        return {}
    chains = _chain_ca_coordinates(structure_path)
    binder = chains.get(binder_chain) or []
    target_group = [
        chain
        for chain in (target_chains or [target_chain])
        if chain != binder_chain and chains.get(chain)
    ]
    if not binder or not target_group:
        return {}
    if len(matrix) < sum(len(coords) for coords in chains.values()):
        return {}

    offsets: dict[str, int] = {}
    offset = 0
    for chain, coords in chains.items():
        offsets[chain] = offset
        offset += len(coords)
    binder_offset = offsets.get(binder_chain, 0)
    pairs: list[tuple[int, int]] = []
    for current_target_chain in target_group:
        target = chains[current_target_chain]
        target_offset = offsets[current_target_chain]
        for binder_index, binder_coord in enumerate(binder):
            for target_index, target_coord in enumerate(target):
                if _distance(binder_coord, target_coord) <= dist_cutoff:
                    pairs.append((binder_offset + binder_index, target_offset + target_index))
    if not pairs:
        return {}

    binder_to_target = [matrix[binder_index][target_index] for binder_index, target_index in pairs]
    target_to_binder = [matrix[target_index][binder_index] for binder_index, target_index in pairs]
    combined = binder_to_target + target_to_binder
    return {
        "ipae": sum(combined) / len(combined),
        "ipae_target_chain_count": len(target_group),
        "ipae_interface_contact_count": len(pairs),
        "ipae_binder_to_target": sum(binder_to_target) / len(binder_to_target),
        "ipae_target_to_binder": sum(target_to_binder) / len(target_to_binder),
        "ipae_contact_pairs": len(pairs),
        "ipae_dist_cutoff": dist_cutoff,
    }


def _run_ipsae(
    source_run_dir: Path,
    run_dir: Path,
    candidates: list[dict[str, Any]],
    pae_cutoff: float = 10.0,
    dist_cutoff: float = 10.0,
    workers: int | None = None,
    ipdae_contact_cutoff: float = 8.0,
    structure_contact_cutoff: float = 4.0,
    hotspot_contact_cutoff: float = 5.0,
) -> list[dict[str, Any]]:
    ipsae_dir = run_dir / "artifacts" / "analysis" / "ipsae"
    input_dir = ipsae_dir / "inputs"
    output_dir = ipsae_dir / "rows"
    logs_dir = ipsae_dir / "logs"
    input_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)
    worker_count = _allocated_cpu_cores(workers)
    steps: list[dict[str, Any]] = []
    enhanced: list[dict[str, Any]] = []
    used_ids: set[str] = set()
    staged_structures: dict[str, Path] = {}
    staged_measurement_structures: dict[str, Path] = {}

    for index, candidate in enumerate(candidates, start=1):
        safe_id = _safe_id(candidate)
        if safe_id in used_ids:
            safe_id = f"{safe_id}_{index}"
        used_ids.add(safe_id)
        metrics = dict(candidate.get("metrics") or {})
        metrics.pop("ipsae_error", None)
        binder_chains, target_chains = _chain_role_groups(candidate)
        if not binder_chains or not target_chains:
            metrics["ipsae_error"] = "binder or target chain roles are missing"
            enhanced.append({**candidate, "metrics": metrics, "_ipsae_batch_id": safe_id})
            continue
        binder_chain, target_chain = binder_chains[0], target_chains[0]
        structure_path = _resolve_path(source_run_dir, candidate.get("complex_pdb"))
        pae_path = _find_pae_path(source_run_dir, candidate, structure_path)
        staged_dir = input_dir / safe_id
        staged_dir.mkdir(parents=True, exist_ok=True)
        row_path = output_dir / f"{safe_id}_ipsae_row.csv"

        if not structure_path or not structure_path.exists():
            metrics["ipsae_error"] = "complex structure missing"
            enhanced.append({**candidate, "metrics": metrics, "_ipsae_batch_id": safe_id})
            continue
        staged_measurement_structure = staged_dir / structure_path.name
        shutil.copy2(structure_path, staged_measurement_structure)
        staged_measurement_structures[safe_id] = staged_measurement_structure
        staged_pae = None
        if pae_path and pae_path.exists():
            staged_pae = staged_dir / pae_path.name
            shutil.copy2(pae_path, staged_pae)
        else:
            metrics["ipsae_error"] = "PAE file missing; structure metrics were still calculated"
        staged_structure = _collapse_target_chains_for_ipsae(
            staged_measurement_structure,
            binder_chain=binder_chain,
            target_chains=target_chains,
        )
        staged_structures[safe_id] = staged_structure
        design_structure, design_binder_chains = _candidate_design_reference(source_run_dir, candidate)
        design_complex, design_complex_binder_chains, design_complex_target_chains = _candidate_design_complex(source_run_dir, candidate)
        staged_design_structure = None
        if design_structure:
            staged_design_structure = staged_dir / f"design_reference{design_structure.suffix.lower()}"
            shutil.copy2(design_structure, staged_design_structure)
        staged_design_complex = None
        if design_complex:
            staged_design_complex = staged_dir / f"design_complex{design_complex.suffix.lower()}"
            shutil.copy2(design_complex, staged_design_complex)
        hotspots = _candidate_hotspots(candidate)
        engine, pae_format = _ipsae_engine(candidate, staged_pae)
        command = [
            "python3",
            str(IPSAE_RUNNER),
            "--id",
            safe_id,
            "--structure",
            f"/work/artifacts/analysis/ipsae/inputs/{safe_id}/{staged_structure.name}",
            "--measurement-structure",
            f"/work/artifacts/analysis/ipsae/inputs/{safe_id}/{staged_measurement_structure.name}",
            "--pae-format",
            pae_format,
            "--engine",
            engine,
            "--binder-chain",
            binder_chain,
            "--target-chain",
            target_chain,
            "--measurement-target-chains",
            ",".join(target_chains),
            "--pae-cutoff",
            str(pae_cutoff),
            "--dist-cutoff",
            str(dist_cutoff),
            "--ipdae-contact-cutoff",
            str(ipdae_contact_cutoff),
            "--structure-contact-cutoff",
            str(structure_contact_cutoff),
            "--hotspot-contact-cutoff",
            str(hotspot_contact_cutoff),
            "--hotspots",
            hotspots,
            "--output",
            f"/work/artifacts/analysis/ipsae/rows/{safe_id}_ipsae_row.csv",
        ]
        if staged_pae:
            command.extend(["--pae", f"/work/artifacts/analysis/ipsae/inputs/{safe_id}/{staged_pae.name}"])
        if staged_design_structure:
            command.extend([
                "--design-structure",
                f"/work/artifacts/analysis/ipsae/inputs/{safe_id}/{staged_design_structure.name}",
                "--design-reference-binder-chains",
                ",".join(design_binder_chains),
            ])
        if staged_design_complex:
            command.extend([
                "--design-complex",
                f"/work/artifacts/analysis/ipsae/inputs/{safe_id}/{staged_design_complex.name}",
                "--design-complex-binder-chains",
                ",".join(design_complex_binder_chains or design_binder_chains),
                "--design-target-chains",
                ",".join(design_complex_target_chains),
            ])
        if staged_pae and staged_pae.suffix.lower() == ".json" and staged_structure.suffix.lower() == ".pdb":
            # Preserve the existing richer AF2/AF3 scorer; the container adapter
            # runs it unchanged, then emits the common PAE/i_pDAE measurements.
            command.extend(["--ipsae-script", str(IPSAE_FULL_RUNNER)])
        steps.append(
            {
                "id": safe_id,
                "name": f"ipsae-{safe_id}",
                "command": command,
                "row_path": str(row_path),
                "log_path": f"/work/artifacts/analysis/ipsae/logs/{safe_id}.log",
            }
        )
        enhanced.append({**candidate, "metrics": metrics, "_ipsae_batch_id": safe_id})

    docker_return_code: int | None = None
    if steps:
        manifest_path = ipsae_dir / "batch_manifest.json"
        status_path = ipsae_dir / "batch_status.json"
        manifest = {
            "work_dir": "/work",
            "status_path": "/work/artifacts/analysis/ipsae/batch_status.json",
            "jobs": [
                {
                    "id": step["id"],
                    "command": step["command"],
                    "log_path": step["log_path"],
                }
                for step in steps
            ],
        }
        write_json(manifest_path, manifest)
        command_payload = {
            "mode": "internal+docker-batch",
            "candidate_count": len(steps),
            "cpu_workers": worker_count,
            "steps": steps,
            "note": "One CPU-limited ipSAE container processes each candidate independently in parallel. Structure contacts, chain-scoped Cα B-factor means, hotspot contacts, and optional design-to-cofold binder RMSD and sequence-mapped designed-contact preservation are calculated from saved structures; ipSAE and PAE-derived measurements require that candidate's saved PAE matrix.",
        }
        docker_command = [
            "docker",
            "run",
            "--rm",
            "-v",
            f"{run_dir}:/work",
            "-v",
            f"{IPSAE_BATCH_RUNNER}:/opt/mds-ipsae-batch-runner.py:ro",
            "-w",
            "/work",
            IPSAE_IMAGE,
            "python3",
            "/opt/mds-ipsae-batch-runner.py",
            "--manifest",
            "/work/artifacts/analysis/ipsae/batch_manifest.json",
            "--workers",
            str(min(worker_count, len(steps))),
        ]
        docker_step = apply_docker_cpu_limits_to_steps(run_dir, [{"name": "ipsae-batch", "command": docker_command}])[0]
        docker_command = docker_step["command"]
        command_payload["docker_command"] = docker_command
        write_json(run_dir / "command.json", command_payload)
        with (run_dir / "stdout.log").open("a") as stdout, (run_dir / "stderr.log").open("a") as stderr:
            stdout.write(f"$ {' '.join(docker_command)}\n")
            stdout.flush()
            try:
                docker_return_code = subprocess.run(docker_command, stdout=stdout, stderr=stderr, check=False).returncode
            except OSError as exc:
                stderr.write(f"Could not launch ipSAE batch container: {exc}\n")
                docker_return_code = 127

    batch_status = read_json(status_path) if steps and status_path.is_file() else {}
    status_by_id = {
        str(row.get("id")): row
        for row in batch_status.get("jobs", [])
        if isinstance(row, dict)
    } if isinstance(batch_status, dict) else {}

    for candidate in enhanced:
        safe_id = str(candidate.pop("_ipsae_batch_id", _safe_id(candidate)))
        row_path = output_dir / f"{safe_id}_ipsae_row.csv"
        binder_chains, target_chains = _chain_role_groups(candidate)
        if not binder_chains or not target_chains:
            continue
        structure_path = _resolve_path(source_run_dir, candidate.get("complex_pdb"))
        pae_path = _find_pae_path(source_run_dir, candidate, structure_path)
        staged_structure = staged_structures.get(safe_id)
        row_metrics = _parse_ipsae_row(row_path) if row_path.exists() else {}
        summary_metrics: dict[str, Any] = {}
        if staged_structure and staged_structure.exists():
            summary_path = _find_ipsae_summary(staged_structure, pae_cutoff, dist_cutoff)
            if summary_path:
                binder_chain, target_chain = _chain_roles(candidate)
                summary_metrics = _parse_ipsae_summary_scores(summary_path, focus_chain=binder_chain)
        if row_metrics or summary_metrics:
            binder_chains, target_chains = _chain_role_groups(candidate)
            pae_metrics = _interface_pae_scores(
                structure_path,
                pae_path,
                binder_chain=binder_chains[0],
                target_chain=target_chains[0],
                target_chains=target_chains,
                dist_cutoff=dist_cutoff,
            ) if pae_path and pae_path.exists() else {}
            merged_metrics = dict(candidate.get("metrics") or {})
            for metric_block in [pae_metrics, row_metrics, summary_metrics]:
                merged_metrics.update(
                    {
                        key: value
                        for key, value in metric_block.items()
                        if value is not None
                    }
                )
            candidate["metrics"] = merged_metrics
        else:
            metrics = dict(candidate.get("metrics") or {})
            status = status_by_id.get(safe_id, {})
            return_code = status.get("return_code")
            if return_code not in (None, 0):
                reason = f"IPSAE command exited with code {return_code}; see {status.get('log_path') or 'candidate log'}."
            elif docker_return_code not in (None, 0):
                reason = f"IPSAE batch container exited with code {docker_return_code}."
            else:
                reason = "IPSAE scorer produced no parsable result"
            metrics.setdefault("ipsae_error", reason)
            candidate["metrics"] = metrics
    return enhanced


def run_ipsae_analysis_contract(
    source_run_dir: Path,
    candidates_jsonl: Path,
    pae_cutoff: float = 10.0,
    dist_cutoff: float = 10.0,
    candidate_ids: list[str] | None = None,
    ipdae_contact_cutoff: float = 8.0,
    structure_contact_cutoff: float = 4.0,
    hotspot_contact_cutoff: float = 5.0,
    selected_structure_path: Path | None = None,
    selected_model_id: str | None = None,
    binder_chains: list[str] | None = None,
    target_chains: list[str] | None = None,
) -> Path:
    """Analyze saved complex PAE and structures without launching a predictor."""
    source_run_dir = Path(source_run_dir)
    candidates_jsonl = Path(candidates_jsonl)
    candidates = _analysis_source_candidates(source_run_dir, candidates_jsonl)
    candidates = _select_complex_analysis_candidates(
        source_run_dir, candidates, candidate_ids=candidate_ids,
        selected_structure_path=selected_structure_path, selected_model_id=selected_model_id,
        binder_chains=binder_chains, target_chains=target_chains,
    )
    if selected_structure_path is not None and (
        len(candidates[0].get("binder_chains") or []) != 1
        or len(candidates[0].get("target_chains") or []) != 1
    ):
        raise ValueError("The current ipSAE adapter requires one selected protein chain in each interface group.")
    selected_ids = {str(candidate.get("candidate_id") or "") for candidate in candidates}
    job = _analysis_job(
        source_run_dir,
        candidates_jsonl,
        tool="ipsae",
        params={"pae_cutoff": pae_cutoff, "dist_cutoff": dist_cutoff,
                "ipdae_contact_cutoff": ipdae_contact_cutoff,
                "structure_contact_cutoff": structure_contact_cutoff,
                "hotspot_contact_cutoff": hotspot_contact_cutoff,
                "candidate_ids": sorted(selected_ids) if candidate_ids is not None else None,
                "selected_structure_path": str(selected_structure_path or ""),
                "selected_model_id": str(selected_model_id or ""),
                "binder_chains": binder_chains, "target_chains": target_chains},
        note="Calculate ipSAE and PAE-derived interface scores when a saved PAE matrix exists; also measure structure contacts, chain-scoped Cα B-factors, configured hotspot contacts, binder-only fitted RMSD, target-aligned binder-pose RMSD, and sequence-mapped designed-interface contact preservation when a separate design complex is available. Native RFD4 refold RMSD remains separately attributed. No predictor is launched.",
        mode="internal+docker",
    )
    if selected_structure_path is not None and len(candidates) == 1:
        _write_complex_analysis_selection(job, source_run_dir, candidates[0], workflow="analysis.ipsae", parameters={
            "pae_cutoff": float(pae_cutoff), "dist_cutoff": float(dist_cutoff),
            "ipdae_contact_cutoff": float(ipdae_contact_cutoff),
            "structure_contact_cutoff": float(structure_contact_cutoff),
            "hotspot_contact_cutoff": float(hotspot_contact_cutoff),
        })
    inputs = [{**copy.deepcopy(candidate), "metrics": {}, "_measurement_source": candidate} for candidate in candidates]
    scored = _run_ipsae(source_run_dir, job.run_dir, inputs, pae_cutoff=float(pae_cutoff),
                        dist_cutoff=float(dist_cutoff), ipdae_contact_cutoff=float(ipdae_contact_cutoff),
                        structure_contact_cutoff=float(structure_contact_cutoff),
                        hotspot_contact_cutoff=float(hotspot_contact_cutoff))
    for candidate in scored:
        raw = candidate.setdefault("raw_metadata", {})
        design_reference, _reference_binder_chains = _candidate_design_reference(source_run_dir, candidate)
        design_complex, _complex_binder_chains, _complex_target_chains = _candidate_design_complex(source_run_dir, candidate)
        raw["analysis_parameters"] = {"pae_cutoff": float(pae_cutoff), "dist_cutoff": float(dist_cutoff),
            "ipdae_contact_cutoff": float(ipdae_contact_cutoff),
            "structure_contact_cutoff": float(structure_contact_cutoff),
            "hotspot_contact_cutoff": float(hotspot_contact_cutoff),
            "image": IPSAE_IMAGE,
            "pae_file": str(_find_pae_path(source_run_dir, candidate, _resolve_path(source_run_dir, candidate.get("complex_pdb"))) or "")}
        raw["measurement_provenance"] = {
            **dict(raw.get("measurement_provenance") or {}),
            "bindcraft2_ipdae_recomputed": {
                "definition": "BindCraft2 anchored_interface_tm_scores; maximum binder- and target-anchored contact scores",
                "contact_cutoff_angstrom": float(ipdae_contact_cutoff),
                "source_pae_matrix": raw["analysis_parameters"]["pae_file"],
                "calculation_location": "mn-protein-scoring-ipsae container adapter",
                "native_bindcraft2_score": False,
            },
            "ipsae_structure_analysis": {
                "calculation_location": "mn-protein-scoring-ipsae container adapter",
                "structure_contact_cutoff_angstrom": float(structure_contact_cutoff),
                "hotspot_contact_cutoff_angstrom": float(hotspot_contact_cutoff),
                "structure_confidence_source": "saved complex C-alpha B-factors; kept distinct from predictor-native confidence summaries and not assumed to be pLDDT",
                "binder_rmsd": "sequence-aligned binder C-alpha RMSD after optimal binder-only superposition; available only when a separate design structure is found",
                "target_aligned_binder_pose_rmsd": "sequence-matched design and cofold complexes; align target C-alpha atoms, then measure binder C-alpha RMSD in that target frame without binder refitting",
                "design_reference_found": bool(design_reference),
                "design_complex_found": bool(design_complex),
                "design_reference_structure": str(design_reference or ""),
                "design_complex_structure": str(design_complex or ""),
                "interaction_preservation": "sequence-mapped binder–target residue contact-pair overlap using the same configured heavy-atom distance cutoff; no claim about specific chemical interactions",
                "hotspots": _candidate_hotspots(candidate),
                "native_prediction_performed": False,
            },
        }

    return _finish_metric_analysis(
        job,
        scored,
        tool="ipsae",
        table_name="ipsae_metrics.csv",
        metric_names=["ipsae", "ipsae_min", "ipsae_max", "ipsae_avg", "lis", "ipae",
                      "ipae_binder_to_target", "ipae_target_to_binder", "ipae_bidirectional",
                      "ipae_binder_to_target_norm31", "ipae_target_to_binder_norm31", "ipae_bidirectional_norm31",
                      "bindcraft2_ipdae_recomputed", "bindcraft2_ipdae_binder_anchor_max",
                      "bindcraft2_ipdae_target_anchor_max", "structure_interface_atom_contact_pair_count",
                      "structure_interface_residue_contact_pair_count", "structure_binder_interface_residue_count",
                      "structure_target_interface_residue_count", "structure_ca_clash_pair_count_2_5a",
                      "structure_binder_ca_bfactor_mean", "structure_target_ca_bfactor_mean",
                      "structure_complex_ca_bfactor_mean", "structure_hotspot_contact_fraction_found",
                      "structure_hotspot_contact_fraction_requested", "structure_hotspot_contacted_count",
                      "structure_hotspot_min_distance", "structure_binder_design_to_cofold_fitted_ca_rmsd",
                      "structure_target_aligned_binder_pose_ca_rmsd", "structure_target_alignment_ca_rmsd",
                      "structure_target_alignment_ca_count", "structure_target_aligned_binder_ca_count",
                      "structure_target_aligned_binder_sequence_coverage",
                      "structure_design_interface_residue_contact_pairs", "structure_cofold_interface_residue_contact_pairs",
                      "structure_interface_contact_pairs_preserved", "structure_design_interaction_preservation_fraction",
                      "structure_interface_contact_jaccard", "structure_binder_interface_residue_preservation_fraction",
                      "structure_target_interface_residue_preservation_fraction", "structure_interaction_sequence_mapping_coverage"],
    )


def run_rosetta_metrics_contract(
    source_run_dir: Path,
    candidates_jsonl: Path,
    nprocs: int | None = None,
    candidate_ids: list[str] | None = None,
    selected_structure_path: Path | None = None,
    selected_model_id: str | None = None,
    binder_chains: list[str] | None = None,
    target_chains: list[str] | None = None,
) -> Path:
    """Score one-binder/one-target complex structures using the installed PyRosetta image."""
    source_run_dir = Path(source_run_dir)
    candidates_jsonl = Path(candidates_jsonl)
    nprocs = _allocated_cpu_cores(nprocs)
    candidates = _analysis_source_candidates(source_run_dir, candidates_jsonl)
    candidates = _select_complex_analysis_candidates(
        source_run_dir, candidates, candidate_ids=candidate_ids,
        selected_structure_path=selected_structure_path, selected_model_id=selected_model_id,
        binder_chains=binder_chains, target_chains=target_chains,
    )
    if selected_structure_path is not None and (
        len(candidates[0].get("binder_chains") or []) != 1
        or len(candidates[0].get("target_chains") or []) != 1
    ):
        raise ValueError("The current Rosetta interface scorer requires one selected protein chain in each interface group.")
    candidates = [{**copy.deepcopy(candidate), "metrics": {}, "_measurement_source": candidate} for candidate in candidates]
    job = _analysis_job(
        source_run_dir,
        candidates_jsonl,
        tool="rosetta_metrics",
        params={"nprocs": nprocs, "candidate_ids": candidate_ids,
                "selected_structure_path": str(selected_structure_path or ""),
                "selected_model_id": str(selected_model_id or ""),
                "binder_chains": binder_chains, "target_chains": target_chains},
        note="Calculate Rosetta interface, packing, SASA, hydrogen-bond, and binder surface metrics for saved complex structures.",
        mode="docker",
    )
    if selected_structure_path is not None and len(candidates) == 1:
        _write_complex_analysis_selection(job, source_run_dir, candidates[0], workflow="analysis.rosetta", parameters={})
    scoring_dir = job.run_dir / "artifacts" / "analysis" / "rosetta" / "inputs"
    source_pdb_dir = job.run_dir / "artifacts" / "analysis" / "rosetta" / "source_pdb"
    relaxed_dir = job.run_dir / "artifacts" / "analysis" / "rosetta" / "relaxed"
    output_csv = job.run_dir / "artifacts" / "analysis" / "rosetta_metrics.csv"
    scoring_dir.mkdir(parents=True, exist_ok=True)
    source_pdb_dir.mkdir(parents=True, exist_ok=True)
    relaxed_dir.mkdir(parents=True, exist_ok=True)

    from mn_molecular_design_studio.engine.workflow_engine.workflows import refolding as refolding_workflow
    from mn_molecular_design_studio.engine.workflow_engine.workflows.sequence_design import (
        _stage_ligandmpnn_role_normalized_input,
    )

    by_safe_id: dict[str, dict[str, Any]] = {}
    for index, candidate in enumerate(candidates, start=1):
        metrics = dict(candidate.get("metrics") or {})
        metrics.pop("rosetta_error", None)
        structure = _resolve_path(source_run_dir, str(candidate.get("complex_pdb") or ""))
        if not structure or not structure.is_file():
            metrics["rosetta_error"] = "complex structure missing"
            candidate["metrics"] = metrics
            continue
        binders, targets = _chain_role_groups(candidate)
        if len(binders) != 1 or len(targets) != 1:
            metrics["rosetta_error"] = "current Rosetta interface scorer needs one binder chain and one target chain; explicit roles required"
            candidate["metrics"] = metrics
            continue
        safe_id = re.sub(r"[^A-Za-z0-9_.-]+", "_", _safe_id(candidate)) or f"candidate_{index}"
        if safe_id in by_safe_id:
            safe_id = f"{safe_id}_{index}"
        source_pdb = structure
        if structure.name.lower().endswith((".cif", ".mmcif", ".cif.gz", ".mmcif.gz")):
            # The PyRosetta batch script scores every PDB in ``scoring_dir``.
            # Keep the converted source there only as an input to chain
            # normalization; otherwise it is also scored with the original
            # chain IDs and fails the normalized binder-chain lookup.
            source_pdb = source_pdb_dir / f"{safe_id}_source.pdb"
            try:
                refolding_workflow._cif_to_pdb(structure, source_pdb)
            except (OSError, ValueError, KeyError) as exc:
                metrics["rosetta_error"] = f"could not convert mmCIF to PDB: {exc}"
                candidate["metrics"] = metrics
                continue
        staged_pdb = scoring_dir / f"{safe_id}.pdb"
        try:
            role_metadata = _stage_ligandmpnn_role_normalized_input(
                source_pdb=source_pdb,
                staged_pdb=staged_pdb,
                candidate=candidate,
                requested_design_chains="",
            )
        except (OSError, ValueError, IndexError, KeyError) as exc:
            metrics["rosetta_error"] = f"could not prepare chain-normalized structure: {exc}"
            candidate["metrics"] = metrics
            continue
        binder_chains = role_metadata.get("binder_chains") or []
        target_chains = role_metadata.get("target_chains") or []
        if len(binder_chains) != 1 or len(target_chains) != 1:
            metrics["rosetta_error"] = (
                "current Rosetta interface scorer needs one binder chain and one target chain; "
                f"found {len(binder_chains)} binder chain(s) and {len(target_chains)} target chain(s)"
            )
            candidate["metrics"] = metrics
            staged_pdb.unlink(missing_ok=True)
            continue
        metadata = dict(candidate.get("raw_metadata") or {})
        metadata["rosetta_chain_role_schema"] = role_metadata.get("chain_role_schema")
        metadata["rosetta_chain_roles"] = role_metadata.get("chain_roles")
        candidate["raw_metadata"] = metadata
        by_safe_id[safe_id] = candidate

    command: list[str] = []
    return_code = None
    if by_safe_id:
        command = [
            "docker",
            "run",
            "--rm",
            "-v",
            f"{job.run_dir}:{job.run_dir}",
            "-w",
            str(PYROSETTA_WORKDIR),
            PYROSETTA_IMAGE,
            "python",
            str(PYROSETTA_SCRIPT),
            "--folder",
            f"rosetta:{scoring_dir}",
            "--out-csv",
            str(output_csv),
            "--relaxed-dir",
            str(relaxed_dir),
            "--binder-chain",
            "Z",
            "--target-chain",
            "A",
            "--nprocs",
            str(nprocs),
            "--dalphaball-path",
            str(PYROSETTA_WORKDIR / "functions" / "DAlphaBall.gcc"),
        ]
        command = apply_docker_cpu_limits_to_steps(job.run_dir, [{"name": "rosetta-batch", "command": command}])[0]["command"]
        write_json(
            job.run_dir / "command.json",
            {
                "mode": "docker",
                "command": command,
                "cpu_workers": nprocs,
                "candidate_count": len(by_safe_id),
                "note": "Standalone PyRosetta analysis. Input structures are copied and chain-normalized; the selected source result is not modified.",
            },
        )
        with (job.run_dir / "stdout.log").open("a") as stdout, (job.run_dir / "stderr.log").open("a") as stderr:
            stdout.write(f"$ {' '.join(command)}\n")
            stdout.flush()
            try:
                return_code = subprocess.run(command, stdout=stdout, stderr=stderr, check=False).returncode
            except OSError as exc:
                stderr.write(f"Could not launch Rosetta scoring container: {exc}\n")
                return_code = 127

    result_rows: dict[str, dict[str, Any]] = {}
    if output_csv.is_file():
        with output_csv.open(newline="", encoding="utf-8", errors="replace") as handle:
            for row in csv.DictReader(handle):
                result_rows[str(row.get("binder_id") or "")] = {
                    key.removeprefix("rosetta_"): _safe_float(value) if _safe_float(value) is not None else value
                    for key, value in row.items()
                    if key != "binder_id" and value not in {None, ""}
                }

    for safe_id, candidate in by_safe_id.items():
        metrics = dict(candidate.get("metrics") or {})
        rosetta_metrics = result_rows.get(safe_id, {})
        if rosetta_metrics:
            metrics.update(rosetta_metrics)
            metrics.pop("rosetta_error", None)
            metrics["rosetta_scoring_status"] = "completed"
        else:
            reason = f"Rosetta scorer produced no row (exit code {return_code})." if return_code is not None else "Rosetta scorer did not run."
            metrics["rosetta_error"] = reason
        candidate["metrics"] = metrics

    for candidate in candidates:
        candidate.setdefault("raw_metadata", {})["analysis_parameters"] = {"image": PYROSETTA_IMAGE, "nprocs": nprocs}
    return _finish_metric_analysis(
        job,
        candidates,
        tool="rosetta_metrics",
        table_name="rosetta_metrics.csv",
        metric_names=["rosetta_binder_score", "rosetta_interface_dG", "rosetta_interface_sc", "rosetta_interface_dSASA"],
    )


def _metric_value(flat_metrics: dict[str, Any], canonical: str) -> float | None:
    aliases = METRIC_ALIASES.get(canonical, [canonical])
    lowered = {key.lower(): value for key, value in flat_metrics.items()}
    for alias in aliases:
        if alias.lower() in lowered:
            value = _safe_float(lowered[alias.lower()])
            if value is not None:
                low, high = RANGES.get(canonical, (None, None))
                if (low is None or value >= low) and (high is None or value <= high):
                    return value
    return None


def _add_hotspot_metrics(source_run_dir: Path, candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    enhanced: list[dict[str, Any]] = []
    for candidate in candidates:
        metrics = dict(candidate.get("metrics") or {})
        structure_path = _resolve_path(source_run_dir, candidate.get("complex_pdb"))
        binder_chains, target_chains = _chain_role_groups(candidate)
        if structure_path and structure_path.exists() and structure_path.suffix.lower() == ".pdb" and binder_chains and target_chains:
            binder_chain, target_chain = _chain_roles(candidate)
            hotspot_metrics = calculate_hotspot_metrics(
                structure_path,
                candidate.get("hotspots") or [],
                binder_chain=str(metrics.get("binder_chain") or binder_chain),
                target_chain=str(metrics.get("target_chain") or target_chain),
                contact_cutoff=8.0,
                contig=str(candidate.get("contig") or ""),
            )
            metrics.update(hotspot_metrics)
        enhanced.append({**candidate, "metrics": metrics})
    return enhanced


def _passes_thresholds(row: dict[str, Any], thresholds: dict[str, float]) -> tuple[bool, list[str]]:
    failures: list[str] = [str(row[key]) for key in ("ipsae_error", "rosetta_error", "validation_error") if row.get(key)]
    checks = [
        ("binder_plddt", ">=", thresholds.get("min_binder_plddt")),
        ("confidence", ">=", thresholds.get("min_confidence")),
        ("iptm", ">=", thresholds.get("min_iptm")),
        ("ipsae", ">=", thresholds.get("min_ipsae")),
        ("ipae", "<=", thresholds.get("max_ipae")),
        ("ipde", "<=", thresholds.get("max_ipde")),
        ("binder_rmsd", "<=", thresholds.get("max_binder_rmsd")),
        ("monomer_rmsd", "<=", thresholds.get("max_binder_rmsd")),
        ("hotspot_contact_fraction", ">=", thresholds.get("min_hotspot_contact_fraction")),
        ("min_binder_to_hotspot_distance", "<=", thresholds.get("max_hotspot_distance")),
    ]
    for metric, op, threshold in checks:
        if threshold is None:
            continue
        value = row.get(metric)
        if _safe_float(value) is None:
            failures.append(f"{metric} missing or non-finite")
            continue
        if op == ">=" and float(value) < float(threshold):
            failures.append(f"{metric} < {threshold:g}")
        if op == "<=" and float(value) > float(threshold):
            failures.append(f"{metric} > {threshold:g}")
    return not failures, failures


def _analysis_attempt_key(candidate: dict[str, Any]) -> str:
    raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
    source_candidate = raw.get("source_candidate") if isinstance(raw.get("source_candidate"), dict) else None
    if source_candidate:
        return _analysis_attempt_key(source_candidate)
    for key in ("source_candidate_id", "sequence_refinement_parent_id"):
        value = raw.get(key)
        if value:
            return str(value)
    parents = [str(parent) for parent in candidate.get("parents") or [] if str(parent)]
    if parents:
        return parents[0]
    return str(candidate.get("candidate_id") or "")


def _ranking_sort_columns(metric: str) -> tuple[list[str], list[bool]]:
    metric = str(metric or "analysis_score")
    if metric in {
        "binder_rmsd",
        "monomer_rmsd",
        "ipae",
        "ipde",
    }:
        return ["passes_filters", metric, "analysis_score"], [False, True, False]
    if metric == "analysis_score":
        return ["passes_filters", "analysis_score"], [False, False]
    return ["passes_filters", metric, "analysis_score"], [False, False, False]


def build_analysis_table(
    candidates: list[dict[str, Any]],
    thresholds: dict[str, float] | None = None,
    *,
    ranking_metric: str = "analysis_score",
) -> pd.DataFrame:
    thresholds = {**DEFAULT_THRESHOLDS, **(thresholds or {})}
    rows: list[dict[str, Any]] = []
    cohort_methods: dict[str, set[str]] = {}
    for source_index, candidate in enumerate(candidates):
        from mn_molecular_design_studio.output_vocabulary import candidate_output
        output = candidate_output(candidate)
        common = output["metrics"]
        flat_metrics = _flatten_metrics(candidate.get("metrics") or {})
        row = {
            "source_index": source_index,
            "candidate_id": candidate.get("candidate_id"),
            "source_tool": candidate.get("source_tool"),
            "complex_pdb": candidate.get("complex_pdb"),
            "binder_pdb": candidate.get("binder_pdb"),
            "binder_sequence": candidate.get("binder_sequence"),
            "parent_count": len(candidate.get("parents") or []),
            "metric_count": len(flat_metrics),
            "attempt_key": _analysis_attempt_key(candidate),
        }
        for metric in METRIC_ALIASES:
            if metric == "ipde":
                backend = str(flat_metrics.get("complex_refolding_backend") or "").lower()
                row[metric] = _metric_value(flat_metrics, metric) if "boltz" in backend else None
            else:
                row[metric] = _metric_value(flat_metrics, metric)
        for metric, field in {"binder_plddt": "confidence.validation.bound_binder.plddt", "iptm": "confidence.validation.complex.iptm",
                              "ipae": "confidence.validation.interface.mean_pae", "ipsae": "interface.ipsae.minimum"}.items():
            if field in common:
                row[metric] = common[field]
            if row.get(metric) is not None:
                cohort_methods.setdefault(metric, set()).add(output.get("provenance", {}).get(field, {}).get("method") or str(candidate.get("source_tool")))
        for metric in BINDCRAFT_METRICS:
            if metric in flat_metrics:
                row[f"bindcraft_{metric}"] = flat_metrics[metric]
        for metric, value in flat_metrics.items():
            if metric in row:
                metric = f"source_{metric}"
            row[metric] = value
        native_pass_filters = flat_metrics.get("pass_filters")
        row["native_pass_filters"] = native_pass_filters if isinstance(native_pass_filters, bool) else None
        row["ipsae_error"] = str(flat_metrics.get("ipsae_error") or "")
        row["analysis_score"] = None
        passed, failures = _passes_thresholds(row, thresholds)
        row["passes_filters"] = passed if thresholds or failures else None
        row["filter_failures"] = "; ".join(failures)
        rows.append(row)
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    score_fields = [field for field in ("binder_plddt", "iptm", "ipsae", "ipae", "ipde", "binder_rmsd", "monomer_rmsd", "hotspot_contact_fraction")
                    if field in df and df[field].notna().all() and len(cohort_methods.get(field, set())) <= 1]
    df["analysis_score"] = None  # legacy column retained; no universal composite score
    df["comparison_metrics"] = ", ".join(score_fields)
    df["metric_coverage"] = df[[field for field in METRIC_ALIASES if field in df]].notna().sum(axis=1)
    df["ranking_status"] = "Choose an explicit ranking measurement" if ranking_metric == "analysis_score" else "Within recorded protocol"
    from mn_molecular_design_studio.candidate_selection import _digest, _protocol
    df["comparison_protocol"] = [_digest(_protocol(candidates[int(index)])) for index in df["source_index"]]
    sort_cols, ascending = _ranking_sort_columns(ranking_metric)
    selected_sort = [(column, direction) for column, direction in zip(sort_cols, ascending) if column in df.columns]
    selected_columns = ["comparison_protocol", *[column for column, _direction in selected_sort], "candidate_id"]
    selected_directions = [True, *[direction for _column, direction in selected_sort], True]
    return df.sort_values(
        selected_columns or ["passes_filters", "analysis_score"],
        ascending=selected_directions or [False, False],
    ).reset_index(drop=True)


def _summary_payload(df: pd.DataFrame, thresholds: dict[str, float]) -> dict[str, Any]:
    summary: dict[str, Any] = {
        "thresholds": thresholds,
        "candidate_count": int(len(df)),
        "passing_count": int(df["passes_filters"].sum()) if "passes_filters" in df else 0,
        "ipsae_error_count": int(df["ipsae_error"].fillna("").astype(bool).sum()) if "ipsae_error" in df else 0,
    }
    numeric_cols = [col for col in METRIC_ALIASES if col in df.columns]
    for col in ["analysis_score", *numeric_cols]:
        if col not in df:
            continue
        series = pd.to_numeric(df[col], errors="coerce").dropna()
        if series.empty:
            continue
        summary[col] = {
            "min": float(series.min()),
            "median": float(series.median()),
            "max": float(series.max()),
        }
    return summary


def run_analysis_contract(
    source_run_dir: Path,
    candidates_jsonl: Path,
    tool: str = "ranking",
    keep_top_n: int = 100,
    keep_per_attempt: int | None = None,
    ranking_metric: str = "analysis_score",
    thresholds: dict[str, float] | None = None,
) -> Path:
    source_run_dir = Path(source_run_dir)
    candidates_jsonl = Path(candidates_jsonl)
    source_metadata = read_json(source_run_dir / "metadata.json")
    source_candidates = _analysis_source_candidates(source_run_dir, candidates_jsonl)

    job = create_job(
        ANALYSIS_GROUP,
        job_type="analysis",
        tool=tool,
        inputs={"source_run_dir": str(source_run_dir), "candidates_jsonl": str(candidates_jsonl)},
        params={
            "keep_top_n": keep_top_n,
            "keep_per_attempt": keep_per_attempt,
            "ranking_metric": ranking_metric,
            "thresholds": {**DEFAULT_THRESHOLDS, **(thresholds or {})},
            "backend": "internal",
        },
    )
    inherited_metadata = {
        "campaign_name": source_metadata.get("campaign_name") or "",
        "campaign_id": source_metadata.get("campaign_id") or source_metadata.get("run_id") or "",
        "upstream_task_group": source_metadata.get("task_group") or source_run_dir.parent.name,
        "upstream_run_id": source_metadata.get("run_id") or source_run_dir.name,
        "upstream_job_code": source_metadata.get("job_code") or "",
    }
    write_json(job.run_dir / "metadata.json", {**read_json(job.run_dir / "metadata.json"), **inherited_metadata})
    write_json(
        job.run_dir / "command.json",
        {
            "mode": "internal",
            "command": [],
            "note": "Internal ranking/filter analysis over saved candidate metrics. The upstream result and its native ranking are preserved.",
        },
    )
    update_status(job.run_dir, "running")
    thresholds = {**DEFAULT_THRESHOLDS, **(thresholds or {})}
    analysis_dir = job.run_dir / "artifacts" / "analysis"
    analysis_dir.mkdir(parents=True, exist_ok=True)
    source_candidates = _add_hotspot_metrics(source_run_dir, source_candidates)
    df = build_analysis_table(source_candidates, thresholds, ranking_metric=ranking_metric)
    df.insert(0, "analysis_rank", range(1, len(df) + 1))
    ranked_csv = analysis_dir / "ranked_candidates.csv"
    pass_fail_csv = analysis_dir / "pass_fail_counts.csv"
    summary_json = analysis_dir / "metrics_summary.json"
    df.to_csv(ranked_csv, index=False)
    if "passes_filters" in df:
        df["passes_filters"].value_counts(dropna=False).rename_axis("passes_filters").reset_index(name="count").to_csv(pass_fail_csv, index=False)
    summary = _summary_payload(df, thresholds)
    summary.update(
        {
            "ranking_metric": ranking_metric,
            "source_run_id": source_metadata.get("run_id") or source_run_dir.name,
            "source_ranking_preserved": True,
            "ranking_scope": "derived_analysis_job",
        }
    )
    write_json(summary_json, summary)

    from mn_molecular_design_studio.candidate_selection import evaluate_selection, write_selection_artifacts
    threshold_fields = {"min_binder_plddt": "binder_plddt", "min_confidence": "confidence", "min_iptm": "iptm",
        "min_ipsae": "ipsae", "max_ipae": "ipae", "max_ipde": "ipde", "max_binder_rmsd": "binder_rmsd",
        "min_hotspot_contact_fraction": "hotspot_contact_fraction", "max_hotspot_distance": "min_binder_to_hotspot_distance"}
    rules = [{"field": threshold_fields[key], "operator": ">=" if key.startswith("min_") else "<=", "value": value}
             for key, value in thresholds.items() if key in threshold_fields and value is not None]
    projections = [{**source_candidates[int(row["source_index"])], "metrics": {**(source_candidates[int(row["source_index"])].get("metrics") or {}),
                    **{key: row[key] for key in METRIC_ALIASES}}} for row in df.to_dict("records")]
    outcomes = evaluate_selection(projections, {"rules": rules, "limit": max(0, int(keep_top_n)),
        "ranking_metric": "" if ranking_metric == "analysis_score" else ranking_metric,
        "direction": "lower" if ranking_metric in {"ipae", "ipde", "binder_rmsd", "monomer_rmsd", "rosetta_interface_dG"} else "higher"})
    if keep_per_attempt is not None:
        attempt_counts = {}
        for candidate in sorted(outcomes, key=lambda item: (item["raw_metadata"]["selection"]["rank"] or 0, item["candidate_id"])):
            decision = candidate["raw_metadata"]["selection"]
            group = (decision["cohort"], _analysis_attempt_key(candidate))
            if decision["selected"]:
                if attempt_counts.get(group, 0) >= max(1, int(keep_per_attempt)):
                    decision["selected"] = False
                    decision["reasons"].append("Per-attempt quota reached")
                else:
                    attempt_counts[group] = attempt_counts.get(group, 0) + 1
    decisions = {row["candidate_id"]: row["raw_metadata"]["selection"] for row in outcomes}
    df["selection_status"] = [decisions[identifier]["status"] for identifier in df["candidate_id"]]
    df["selected"] = [decisions[identifier]["selected"] for identifier in df["candidate_id"]]
    df.to_csv(ranked_csv, index=False)
    ranked = [source_candidates[int(index)] for index in df["source_index"]]
    ranked_rows = df.to_dict("records")
    candidates: list[dict[str, Any]] = []
    for rank, source in enumerate(ranked, start=1):
        row = ranked_rows[rank - 1] if rank - 1 < len(ranked_rows) else {}
        metrics = {}
        metrics.update(
            {
                "analysis_rank": decisions[source["candidate_id"]]["rank"],
                "analysis_score": row.get("analysis_score"),
                "analysis_ranking_metric": ranking_metric,
                "analysis_attempt_key": row.get("attempt_key"),
                "analysis_keep_per_attempt": keep_per_attempt,
                "passes_filters": row.get("passes_filters"),
                "filter_failures": row.get("filter_failures"),
                "analysis_backend": "internal",
                "result_kind": "app_reanalysis",
            }
        )
        candidates.append(
            {
                **source,
                "candidate_id": f"{source.get('candidate_id')}_{tool}_{rank:03d}",
                "stage": STAGE_ANALYSIS,
                "source_tool": tool,
                "tool": tool,
                "metrics": metrics,
                "parents": [str(source.get("candidate_id") or "")],
                "raw_metadata": {
                    **dict(source.get("raw_metadata") or {}),
                    "result_kind": "app_reanalysis",
                    "source_candidate": source,
                    "backend_status": "internal",
                    "selection": decisions[source["candidate_id"]],
                    "measurement_status": {"status": "not_measured", "error": "Ranking reads source-stage measurements; it does not measure new properties."},
                },
            }
        )

    normalized = write_candidates(job.run_dir, tool, candidates)
    write_selection_artifacts(analysis_dir / "selection", candidates)
    if rules or ranking_metric != "analysis_score":
        scored_count = sum(decision["status"] in {"pass", "reject"} for decision in decisions.values())
    else:
        scored_count = sum(int(row["metric_coverage"]) > 0 and not any(pd.notna(row.get(key)) and bool(row.get(key)) for key in
            ("ipsae_error", "rosetta_error", "qc_error", "validation_error")) for row in ranked_rows)
    scientific_status = "failed" if not scored_count else "valid" if scored_count == len(candidates) else "partial"
    summary.update(scored_count=scored_count, scientific_status=scientific_status,
                   selected_count=sum(row["raw_metadata"]["selection"]["selected"] for row in candidates))
    write_json(summary_json, summary)
    artifacts = _candidate_artifacts(job.run_dir)
    finish_job(
        job.run_dir,
        scored_count > 0,
        {
            "outputs": {"artifacts": artifacts, "candidates": normalized},
            "metrics": {
                "candidate_count": len(normalized),
                "passing_count": summary.get("passing_count", 0),
                "backend_status": "internal",
            },
            "downstream_artifacts": {
                "candidates_jsonl": "artifacts/normalized_candidates/candidates.jsonl",
                "campaign_result": "artifacts/normalized_candidates/campaign_result.json",
            },
            "summary": summary,
            "execution_status": "completed",
            "scientific_status": summary["scientific_status"],
        },
    )
    return job.run_dir
