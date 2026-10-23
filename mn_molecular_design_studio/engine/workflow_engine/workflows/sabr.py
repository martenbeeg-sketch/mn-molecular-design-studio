"""Standalone SAbR antibody numbering and CDR annotation workflow."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from mn_molecular_design_studio.engine.workflow_engine.core.artifacts import Artifact, artifact_path
from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import JobPaths, create_job, finish_job, write_json
from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path
from mn_molecular_design_studio.engine.workflow_engine.workflows.mpnn import STAGES
from mn_molecular_design_studio import sabr_runtime


def validate_sabr_request(parameters: dict[str, Any], *, require_files: bool = True) -> None:
    structure_path = parameters.get("structure_path")
    source_run_dir = parameters.get("source_run_dir")
    candidates_jsonl = parameters.get("candidates_jsonl")
    if bool(structure_path) == bool(source_run_dir and candidates_jsonl):
        raise ValueError("Choose either one structure file or an existing candidate collection.")
    if (source_run_dir is None) != (candidates_jsonl is None):
        raise ValueError("Existing results need both source_run_dir and candidates_jsonl.")
    if parameters.get("mode", "sabr") not in {"sabr", "softalign"}:
        raise ValueError("SAbR mode must be 'sabr' or 'softalign'.")
    if parameters.get("scheme", "imgt") not in {"imgt", "chothia", "kabat", "martin", "aho", "wolfguy"}:
        raise ValueError("Choose a supported SAbR numbering scheme.")
    noise = float(parameters.get("noise_level", 0.0))
    if noise not in {0.0, 0.2, 0.5, 1.0, 2.0}:
        raise ValueError("SAbR noise level must be one of 0, 0.2, 0.5, 1, or 2 Å.")
    if require_files:
        if structure_path and not Path(structure_path).expanduser().is_file():
            raise ValueError("The selected antibody structure does not exist.")
        if source_run_dir and not Path(source_run_dir).expanduser().is_dir():
            raise ValueError("The selected source job does not exist.")
        if candidates_jsonl and not Path(candidates_jsonl).expanduser().is_file():
            raise ValueError("The selected candidate collection does not exist.")


def run_sabr_annotation(
    structure_path: Path | None = None,
    source_run_dir: Path | None = None,
    candidates_jsonl: Path | None = None,
    selected_candidate_ids: list[str] | None = None,
    chain_ids: list[str] | None = None,
    chain_type: str = "auto",
    scheme: str = "imgt",
    mode: str = "sabr",
    noise_level: float = 0.0,
    write_numbered_structures: bool = True,
    job_name: str = "sabr-antibody-numbering",
    existing_job: JobPaths | None = None,
) -> Path:
    """Annotate antibody CDRs with SAbR and optionally write renumbered copies."""
    validate_sabr_request({
        "structure_path": structure_path, "source_run_dir": source_run_dir,
        "candidates_jsonl": candidates_jsonl, "mode": mode, "scheme": scheme,
        "noise_level": noise_level,
    })
    if mode not in {"sabr", "softalign"}:
        raise ValueError("SAbR mode must be 'sabr' or 'softalign'.")
    if scheme not in {"imgt", "chothia", "kabat", "martin", "aho", "wolfguy"}:
        raise ValueError("Choose a supported SAbR numbering scheme.")
    if not 0 <= float(noise_level) <= 2 or float(noise_level) not in {0.0, 0.2, 0.5, 1.0, 2.0}:
        raise ValueError("SAbR noise level must be one of 0, 0.2, 0.5, 1, or 2 Å.")
    if bool(structure_path) == bool(source_run_dir and candidates_jsonl):
        raise ValueError("Choose either one structure file or an existing candidate collection.")
    if (source_run_dir is None) != (candidates_jsonl is None):
        raise ValueError("Existing results need both source_run_dir and candidates_jsonl.")

    inputs = {"structure_path": str(structure_path)} if structure_path else {
        "source_run_dir": str(source_run_dir), "candidates_jsonl": str(candidates_jsonl)
    }
    params = {
        "chain_ids": chain_ids or [], "chain_type": chain_type, "scheme": scheme,
        "mode": mode, "noise_level": float(noise_level),
        "write_numbered_structures": bool(write_numbered_structures),
        "selected_candidate_ids": selected_candidate_ids or [],
    }
    job = existing_job or create_job(
        "analysis", job_type="antibody_numbering", tool="sabr",
        inputs=inputs, params=params,
    )
    items: list[dict[str, Any]] = []
    try:
        if structure_path:
            source_rows = [{"candidate_id": "structure", "complex_pdb": str(Path(structure_path).resolve()),
                           "binder_chains": chain_ids or []}]
        else:
            source_root = Path(source_run_dir)
            source_rows = [row for row in read_candidates(Path(candidates_jsonl))
                           if row.get("stage") in STAGES
                           and (not selected_candidate_ids or row.get("candidate_id") in selected_candidate_ids)
                           and (row.get("complex_pdb") or row.get("binder_pdb"))]
            if not source_rows:
                raise ValueError("No selected candidates have a designable structure.")

        annotation_bundle: dict[str, Any] = {
            "tool": "SAbR", "scheme": scheme, "mode": mode,
            "noise_level": float(noise_level), "items": items,
        }
        for candidate in source_rows:
            candidate_id = str(candidate.get("candidate_id") or "structure")
            raw_path = candidate.get("complex_pdb") or candidate.get("binder_pdb")
            path = Path(raw_path).resolve() if structure_path else resolve_stored_path(
                raw_path, run_dir=Path(source_run_dir), must_exist=True
            )
            row: dict[str, Any] = {"candidate_id": candidate_id, "structure": str(path), "status": "failed"}
            try:
                if path is None or not path.is_file():
                    raise ValueError("The candidate structure is missing.")
                chains = list(chain_ids or candidate.get("binder_chains") or [])
                if not chains:
                    raise ValueError("Choose an antibody chain or select a candidate with binder-chain metadata.")
                annotation = sabr_runtime.annotate_file(
                    path, chains, chain_type=chain_type, noise_level=float(noise_level),
                    mode=mode, scheme=scheme,
                )
                row.update(status="completed", annotation=annotation)
                if write_numbered_structures:
                    numbered_dir = job.run_dir / "artifacts" / "sabr" / "numbered"
                    numbered_dir.mkdir(parents=True, exist_ok=True)
                    suffix = path.suffix if path.suffix.lower() in {".pdb", ".ent", ".cif", ".mmcif"} else ".pdb"
                    numbered_path = numbered_dir / f"{candidate_id}_numbered{suffix}"
                    sabr_runtime.write_renumbered(
                        path, numbered_path, chains, scheme=scheme, chain_type=chain_type,
                        noise_level=float(noise_level), mode=mode,
                    )
                    row["numbered_structure"] = str(numbered_path)
            except Exception as exc:
                row["error"] = str(exc)
            items.append(row)

        bundle_path = artifact_path(job.run_dir, "sabr", "annotations.json")
        write_json(bundle_path, annotation_bundle)
        successes = sum(item["status"] == "completed" for item in items)
        artifacts = [Artifact("SAbR annotations", bundle_path, "antibody_annotation").to_json(job.run_dir)]
        for item in items:
            numbered = item.get("numbered_structure")
            if numbered:
                artifacts.append(Artifact(Path(numbered).stem, Path(numbered), "structure").to_json(job.run_dir))
        success = successes > 0
        finish_job(job.run_dir, success, {
            "outputs": {
                "tool_status": {"success": success, "completed": successes, "failed": len(items) - successes},
                "items": items, "artifacts": artifacts,
            },
            "metrics": {"candidate_count": len(items), "annotated_count": successes},
            "downstream_artifacts": {"annotations": "artifacts/sabr/annotations.json"},
            "warning": "Some selected structures could not be numbered." if successes < len(items) else "",
        })
    except Exception as exc:
        finish_job(job.run_dir, False, {"outputs": {"items": items}, "error": str(exc)})
        raise
    return job.run_dir
