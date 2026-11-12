"""Tool-independent structure and provenance workspace for saved complexes."""
from __future__ import annotations

import gzip
import hashlib
import html
import io
import json
import os
import re
import sqlite3
import tempfile
from pathlib import Path, PurePosixPath
import streamlit as st


STRUCTURE_SUFFIXES = {".pdb", ".ent", ".cif", ".mmcif", ".pdbqt", ".sdf", ".mol"}
VIEWER_MAX_BYTES = 50 * 1024 * 1024


def initialize_workspace_view(session_state) -> None:
    """Set the Database-only analysis workspace without replacing Inspector links."""
    if session_state.get("complex-analysis-open") or session_state.get("complex-analysis-active"):
        session_state["complex-analysis-workspace-view"] = "inspector"
        return
    # Old sessions may retain the former run-first dashboard or Results board.
    # They are presentation state only and must not make the Results catalog a
    # dependency of ordinary Complex Analysis page loads.
    session_state["complex-analysis-workspace-view"] = "database"
    session_state.pop("complex-analysis-source-kind", None)
    session_state["complex-analysis-dashboard-version"] = "database-only-v1"


def _structure_suffix(path: Path) -> str:
    name = path.name.casefold()
    if name.endswith(".pdb.gz"):
        return ".pdb"
    if name.endswith((".cif.gz", ".mmcif.gz")):
        return ".cif" if name.endswith(".cif.gz") else ".mmcif"
    return path.suffix.casefold()


def open_database_complex(molecule_id: str, *, evaluation_id: str = "", structure_id: str = "") -> None:
    """Navigate by stable Database, evaluation, and structure identifiers."""
    st.session_state["complex-analysis-molecule-id"] = str(molecule_id)
    st.session_state[f"complex-analysis-evaluation-{molecule_id}"] = str(evaluation_id or "")
    st.session_state[f"complex-analysis-structure-id-{molecule_id}"] = str(structure_id or "")
    st.session_state["complex-analysis-open"] = {
        "kind": "database", "molecule_id": str(molecule_id),
        "evaluation_id": str(evaluation_id or ""), "structure_id": str(structure_id or ""),
    }
    st.session_state["complex-analysis-source-kind"] = "Database"
    st.session_state["complex-analysis-workspace-view"] = "inspector"
    st.session_state["page"] = "Complex Analysis"
    st.rerun()


def _registered_cross_run_structure(user: dict, structure_path: Path, *, source_run_dir: Path,
                                    runs_root: Path) -> tuple[str, str, str]:
    """Authorize and identify a structure artifact in a linked Results job."""
    from mn_molecular_design_studio.artifact_integrity import authorize_result_artifact
    return authorize_result_artifact(user, structure_path, source_run_dir=source_run_dir, runs_root=runs_root)


def linked_result_structure(run_dir: Path | str, candidate: dict, *, user: dict) -> dict | None:
    """Resolve one candidate's cross-run structure through its authorized Results provenance.

    This intentionally reuses the Results candidate resolver and the same source-job
    authorization/manifest checks as ``open_result_complex``. It only describes
    child-run structures; artifacts physically stored in the parent are listed by
    the ordinary Results artifact table.
    """
    from mn_molecular_design_studio import cofolding_results
    from mn_molecular_design_studio.runtime import runs_root as configured_runs_root

    parent = Path(run_dir).expanduser().resolve(strict=True)
    root = Path(configured_runs_root()).expanduser().resolve(strict=True)
    structure = cofolding_results.structure_for(parent, candidate)
    if structure is None:
        return None
    try:
        structure = structure.expanduser().resolve(strict=True)
        structure.relative_to(root)
    except (OSError, ValueError) as exc:
        raise ValueError("The candidate structure is outside Studio's managed Results storage.") from exc
    try:
        structure.relative_to(parent)
    except ValueError:
        pass
    else:
        return None

    task_group, run_id, artifact_reference = _registered_cross_run_structure(
        user, structure, source_run_dir=parent, runs_root=root,
    )
    raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
    source_candidate = raw.get("source_candidate") if isinstance(raw.get("source_candidate"), dict) else {}
    source_raw = source_candidate.get("raw_metadata") if isinstance(source_candidate.get("raw_metadata"), dict) else {}
    component_rows = raw.get("components") if isinstance(raw.get("components"), list) else []
    source_components = source_raw.get("components") if isinstance(source_raw.get("components"), list) else []
    receptor = next((item for item in component_rows if isinstance(item, dict) and item.get("role") == "receptor"), {})
    ligand = next((item for item in component_rows if isinstance(item, dict) and item.get("role") == "ligand"), {})
    if not receptor:
        receptor = next((item for item in source_components if isinstance(item, dict) and item.get("role") == "receptor"), {})
    if not ligand:
        ligand = next((item for item in source_components if isinstance(item, dict) and item.get("role") == "ligand"), {})
    receptor_identity = raw.get("receptor_identity") if isinstance(raw.get("receptor_identity"), dict) else {}
    ligand_identity = raw.get("ligand_identity") if isinstance(raw.get("ligand_identity"), dict) else {}
    hashes = raw.get("source_artifact_sha256") if isinstance(raw.get("source_artifact_sha256"), dict) else {}
    source_hashes = source_raw.get("source_artifact_sha256") if isinstance(source_raw.get("source_artifact_sha256"), dict) else {}
    return {
        "path": structure,
        "candidate": candidate,
        "candidate_id": str(candidate.get("candidate_id") or ""),
        "artifact_reference": artifact_reference,
        "producing_job": f"{task_group}/{run_id}",
        "artifact_task_group": task_group,
        "artifact_run_id": run_id,
        "structure_sha256": hashlib.sha256(structure.read_bytes()).hexdigest(),
        "source_pose_id": str(raw.get("pose_id") or source_candidate.get("candidate_id") or ""),
        "source_job_reference": str(raw.get("source_job_reference") or source_raw.get("source_job_reference") or ""),
        "source_pose_reference": str(
            (source_raw.get("source_artifact_references") or {}).get("pose")
            or source_candidate.get("complex_pdb") or ""
        ),
        "receptor_reference": str(
            (source_raw.get("source_artifact_references") or {}).get("receptor")
            or source_candidate.get("target_pdb") or ""
        ),
        "receptor_sha256": str(receptor_identity.get("structure_sha256") or receptor.get("structure_sha256")
                                or hashes.get("receptor") or source_hashes.get("receptor") or ""),
        "pose_sha256": str(hashes.get("pose") or ligand.get("pose_sha256") or source_hashes.get("pose") or ""),
        "ligand_name": str(ligand_identity.get("name") or ligand_identity.get("compound_id")
                           or ligand.get("compound_id") or raw.get("compound_name") or ""),
        "ligand_id": str(ligand_identity.get("compound_id") or ligand.get("compound_id")
                          or raw.get("compound_id") or ""),
        "ligand_smiles": str(ligand_identity.get("smiles") or ligand.get("smiles") or ""),
    }


def open_result_complex(run_dir: Path | str, structure_path: Path | str, candidate: dict | None = None,
                        *, user: dict) -> None:
    """Navigate to an authorized Results structure, including a linked child-run artifact."""
    from mn_molecular_design_studio.auth import can_access_job
    from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path
    from mn_molecular_design_studio.job_catalog import get_job_access_record
    from mn_molecular_design_studio.runtime import runs_root

    root = Path(runs_root()).expanduser().resolve()
    run_dir = Path(run_dir).expanduser().resolve(strict=True)
    try:
        parent_relative = run_dir.relative_to(root)
    except ValueError as exc:
        raise ValueError("The source Results job is outside Studio's managed Results storage.") from exc
    if len(parent_relative.parts) < 2:
        raise ValueError("The source Results job reference is invalid.")
    source_task_group, source_run_id = parent_relative.parts[:2]
    source_access = get_job_access_record(source_task_group, source_run_id)
    if not source_access or not can_access_job(user, source_access):
        raise PermissionError("The source Results job is unavailable or you do not have access to it.")
    resolved_structure = resolve_stored_path(str(structure_path), run_dir=run_dir, must_exist=True)
    if resolved_structure is None:
        raise ValueError("The selected structure artifact is unavailable.")
    try:
        resolved_structure = resolved_structure.expanduser().resolve(strict=True)
        resolved_structure.relative_to(root)
    except (OSError, ValueError) as exc:
        raise ValueError("The selected structure is outside Studio's managed Results storage.") from exc
    task_group, run_id, relative_path = _registered_cross_run_structure(
        user, resolved_structure, source_run_dir=run_dir, runs_root=root)
    from mn_molecular_design_studio.artifact_integrity import retained_artifact_file
    verified_structure = retained_artifact_file(root / task_group / run_id, resolved_structure)
    candidate = candidate if isinstance(candidate, dict) else {}
    raw_metadata = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
    if candidate:
        from mn_molecular_design_studio import cofolding_results
        candidate_structure = cofolding_results.structure_for(run_dir, candidate)
        if candidate_structure is None or candidate_structure.resolve(strict=True) != resolved_structure:
            raise ValueError("The selected structure does not match this Results candidate.")
    model_id = candidate.get("model_id")
    if model_id in (None, ""):
        model_id = candidate.get("model_index")
    if model_id in (None, ""):
        model_id = raw_metadata.get("model_id")
    if model_id in (None, ""):
        model_id = raw_metadata.get("model_index")
    st.session_state["complex-analysis-results-job"] = f"{source_task_group}/{source_run_id}"
    st.session_state["complex-analysis-results-structure"] = relative_path
    st.session_state["complex-analysis-open"] = {
        "kind": "results",
        "task_group": source_task_group,
        "run_id": source_run_id,
        "artifact_task_group": task_group,
        "artifact_run_id": run_id,
        "candidate_id": str(candidate.get("candidate_id") or ""),
        "candidate": candidate,
        "model_index": str(model_id if model_id is not None else ""),
        "relative_path": relative_path,
        "structure_sha256": hashlib.sha256(verified_structure.read_bytes()).hexdigest(),
    }
    st.session_state["complex-analysis-source-kind"] = "Results"
    st.session_state["complex-analysis-workspace-view"] = "inspector"
    st.session_state["page"] = "Complex Analysis"
    # This opener can be called from Results' periodic fragment. Force a full
    # app rerun so the top-level page router switches away from Job details.
    st.rerun(scope="app")


def _remote_cache_root(user: dict) -> Path:
    from mn_molecular_design_studio.cluster.config import home as cluster_home

    cache_root = Path(cluster_home()).expanduser() / "complex-analysis-artifacts"
    owner = hashlib.sha256(str(user.get("id") or "").encode()).hexdigest()[:24]
    cache_root = (cache_root / owner).resolve()
    cache_root.mkdir(parents=True, exist_ok=True, mode=0o700)
    try:
        cache_root.chmod(0o700)
    except OSError:
        pass
    return cache_root


def _stage_remote_result_file(content: bytes, suffix: str, *, user: dict) -> Path:
    """Cache one authorized peer artifact under an opaque, content-addressed name."""
    cache_root = _remote_cache_root(user)
    suffix = suffix if suffix in STRUCTURE_SUFFIXES else ".pdb"
    digest = hashlib.sha256(content).hexdigest()
    path = cache_root / f"{digest}{suffix}"
    path.resolve().relative_to(cache_root)
    if path.is_file() and hashlib.sha256(path.read_bytes()).hexdigest() == digest:
        return path
    with tempfile.NamedTemporaryFile(prefix=".stage-", dir=cache_root, delete=False) as stream:
        temporary = Path(stream.name)
        stream.write(content)
        stream.flush()
        os.fsync(stream.fileno())
    try:
        temporary.chmod(0o600)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)
    if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
        path.unlink(missing_ok=True)
        raise ValueError("The retrieved remote structure failed its content-integrity check.")
    return path


def _remote_candidate_manifest(candidate: dict, structure: Path, *, user: dict) -> Path:
    """Create the minimal local picker manifest needed by existing Inspector controls."""
    row = dict(candidate)
    row["complex_pdb"] = str(structure)
    row.setdefault("raw_metadata", {})
    encoded = (json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n").encode("utf-8")
    digest = hashlib.sha256(encoded).hexdigest()
    root = (_remote_cache_root(user) / "sources" / digest).resolve()
    root.mkdir(parents=True, exist_ok=True, mode=0o700)
    (root / "artifacts" / "normalized_candidates").mkdir(parents=True, exist_ok=True, mode=0o700)
    path = root / "artifacts" / "normalized_candidates" / "candidates.jsonl"
    if not path.is_file() or path.read_bytes() != encoded:
        path.write_bytes(encoded)
        path.chmod(0o600)
    return path


def open_remote_result_complex(job: dict, candidate: dict, artifact_path: str, *, user: dict,
                               structure_content: bytes, source_run_dir: str,
                               structure_field: str, structure_sha256: str = "") -> None:
    """Open one authorized peer-owned result artifact in the existing Inspector."""
    from mn_molecular_design_studio.cluster import ui as cluster_ui

    job = job if isinstance(job, dict) else {}
    candidate = candidate if isinstance(candidate, dict) else {}
    job_id = str(job.get("cluster_job_id") or "")
    node_id = str(job.get("node_id") or "")
    relative = cluster_ui._safe_run_artifact_path(artifact_path)
    if not job_id or not node_id or not relative:
        raise ValueError("The remote Results job or selected artifact reference is invalid.")
    if structure_field not in {"complex_pdb", "binder_pdb"}:
        raise ValueError("The selected candidate does not identify a supported structure field.")
    candidate_structure = candidate.get(structure_field)
    if not isinstance(candidate_structure, str) or not candidate_structure.strip():
        raise ValueError("The selected candidate does not identify this structure artifact.")
    from pathlib import PurePosixPath
    candidate_path = PurePosixPath(candidate_structure.strip())
    run_path = PurePosixPath(str(source_run_dir or ""))
    if candidate_path.is_absolute():
        try:
            candidate_path = candidate_path.relative_to(run_path)
        except ValueError as exc:
            raise ValueError("The candidate structure is outside its remote source job.") from exc
    if (candidate_path.is_absolute() or ".." in candidate_path.parts
            or cluster_ui._safe_run_artifact_path(candidate_path.as_posix()) != relative):
        raise ValueError("The selected structure does not match this candidate's registered artifact.")
    token = str(st.session_state.get("studio-session-token") or "development-mode")
    authorized_job = cluster_ui._request("/jobs/" + job_id, token)
    if (str(authorized_job.get("node_id") or "") != node_id
            or str(authorized_job.get("local_reference") or "") != str(job.get("local_reference") or "")):
        raise PermissionError("The remote Results job changed or is no longer accessible.")
    local_reference = str(authorized_job.get("local_reference") or "")
    reference = PurePosixPath(local_reference)
    run_parts = run_path.parts
    if (reference.is_absolute() or len(reference.parts) != 2 or ".." in reference.parts
            or tuple(run_parts[-2:]) != reference.parts):
        raise ValueError("The selected candidate run does not match its authorized remote job.")
    # The selected bytes were retrieved through this job's artifact action. The
    # source agent resolves local_reference and verifies the file stays inside
    # that run; coordinator manifests can be incomplete for historical runs.
    if not str(candidate.get("candidate_id") or ""):
        raise ValueError("The remote structure has no candidate identity.")
    if not isinstance(structure_content, bytes) or not structure_content:
        raise ValueError("The selected remote structure has not been retrieved.")
    digest = hashlib.sha256(structure_content).hexdigest()
    if structure_sha256 and digest != str(structure_sha256):
        raise ValueError("The selected remote structure failed its content-integrity check.")
    if relative.casefold().endswith(".gz"):
        raise ValueError("Compressed remote structures are not yet supported in the Scientific Inspector.")
    suffix = _structure_suffix(Path(relative))
    if suffix not in STRUCTURE_SUFFIXES:
        raise ValueError("The selected remote artifact is not a supported structure format.")
    artifact_cache_key = "complex-analysis-remote-artifact-" + hashlib.sha256(
        f"{job_id}:{relative}".encode()
    ).hexdigest()[:20]
    st.session_state[artifact_cache_key] = structure_content
    staged_structure = _stage_remote_result_file(structure_content, suffix, user=user)
    candidate_file = _remote_candidate_manifest(candidate, staged_structure, user=user)
    source_run_dir = candidate_file.parents[2]
    artifact_relative = PurePosixPath(relative)
    if (reference.is_absolute() or len(reference.parts) < 2 or ".." in reference.parts
            or artifact_relative.is_absolute() or ".." in artifact_relative.parts):
        raise ValueError("The remote Results reference is invalid.")
    cluster_uri = f"cluster://{node_id}/" + (reference / artifact_relative).as_posix()
    st.session_state["complex-analysis-open"] = {
        "kind": "remote_results",
        "cluster_job_id": job_id,
        "node_id": node_id,
        "local_reference": str(authorized_job.get("local_reference") or ""),
        "task_group": str(authorized_job.get("task_group") or ""),
        "run_id": str(authorized_job.get("run_id") or ""),
        "project_id": str(authorized_job.get("project_id") or ""),
        "campaign_name": str(authorized_job.get("campaign_name") or authorized_job.get("tool") or "Remote Results"),
        "candidate_id": str(candidate.get("candidate_id") or ""),
        "model_index": str(candidate.get("model_index") or (candidate.get("raw_metadata") or {}).get("model_index") or ""),
        "candidate": candidate,
        "artifact_path": relative,
        "structure_sha256": digest,
        "staged_structure": str(staged_structure),
        "candidate_file": str(candidate_file),
        "staged_run_dir": str(source_run_dir),
        "analysis_structure_uri": cluster_uri,
    }
    st.session_state["complex-analysis-source-kind"] = "Results"
    st.session_state["complex-analysis-workspace-view"] = "inspector"
    st.session_state["page"] = "Complex Analysis"
    # This button is rendered inside the Results fragment; route through the
    # top-level app so the remote candidate opens in the Inspector page.
    st.rerun(scope="app")


def component_inventory(content: bytes, filename: str, stored_components: list[dict] | None = None) -> list[dict]:
    """Read the component inventory through the existing PDB/mmCIF import parser."""
    from mn_molecular_design_studio.molecule_imports import parse_structure

    parsed: list[dict] = []
    parse_error = ""
    try:
        rows = parse_structure(content, filename, include_ligands=True)
        if rows:
            parsed = [dict(item) for item in rows[0].get("components", []) if isinstance(item, dict)]
    except (ImportError, OSError, ValueError, TypeError) as exc:
        parse_error = str(exc)

    stored = [dict(item) for item in (stored_components or []) if isinstance(item, dict)]
    if not parsed:
        parsed = stored
    else:
        # Keep stored annotations such as chain roles while treating the parsed
        # file as authoritative for which components are present in this model.
        by_key = {_component_match_key(item): item for item in stored}
        parsed = [{**by_key.get(_component_match_key(item), {}), **item} for item in parsed]

    result = []
    counts: dict[str, int] = {}
    for item in parsed:
        kind = str(item.get("type") or "other")
        chain = str(item.get("chain") or item.get("auth_asym_id") or "_")
        residue_key = str(item.get("residue_key") or "")
        source_chain = residue_key.split(":", 1)[0] if residue_key else chain
        key = f"{kind}:{residue_key or chain}"
        counts[key] = counts.get(key, 0) + 1
        component_id = key if counts[key] == 1 else f"{key}:{counts[key]}"
        result.append({
            **item,
            "id": component_id,
            "type": kind,
            "chain": chain,
            "source_chain": source_chain,
            "residue_key": residue_key,
            "label": str(item.get("name") or item.get("ccd") or item.get("role") or kind.title()),
            "parse_error": parse_error,
        })
    return result


def _component_match_key(item: dict) -> str:
    kind = str(item.get("type") or "other")
    residue_key = str(item.get("residue_key") or "")
    chain = str(item.get("chain") or item.get("auth_asym_id") or "_")
    return f"{kind}:{residue_key or chain}"


def _viewer_text(path: Path, context: dict | None = None) -> str:
    context = context if isinstance(context, dict) else {}
    suffix = _structure_suffix(path)
    if path.suffix.casefold() == ".gz":
        raw = gzip.decompress(path.read_bytes())
    else:
        raw = path.read_bytes()
    text = raw.decode("utf-8", errors="replace")
    if suffix in {".sdf", ".mol"}:
        from rdkit import Chem
        supplier = Chem.ForwardSDMolSupplier(io.BytesIO(raw), removeHs=False)
        molecule = supplier[int(context.get("molecule_index", 0))]
        if molecule is None:
            raise ValueError("Could not read the selected ligand coordinates.")
        return Chem.MolToPDBBlock(molecule)
    if suffix == ".pdbqt":
        lines = []
        for line in text.splitlines():
            if line.startswith(("ATOM  ", "HETATM")):
                atom_type = line.split()[-1]
                element = {"A": "C", "NA": "N", "OA": "O", "SA": "S", "HD": "H"}.get(atom_type, atom_type.title())
                lines.append(line[:66].ljust(76) + element.rjust(2))
        return "\n".join([*lines, "END", ""])
    return text


def _evaluation_options(record: dict) -> tuple[list[dict], dict[str, tuple[dict, dict]]]:
    from mn_molecular_design_studio.evaluations import evaluation_descriptor

    rows = []
    by_id = {}
    for observation in record.get("observations", []):
        try:
            descriptor = evaluation_descriptor(observation)
        except (KeyError, TypeError, ValueError):
            continue
        by_id[descriptor["id"]] = (observation, descriptor)
        rows.append(descriptor)
    rows.sort(key=lambda row: (not row.get("preferred", False), row.get("created_at") or ""), reverse=True)
    return rows, by_id


def _evaluation_label(item: dict) -> str:
    candidate = item.get("candidate_id") or item.get("id", "")[-8:]
    return " · ".join(value for value in (
        str(item.get("evaluator") or "Unknown evaluator"),
        str(item.get("run_id") or "No run"),
        str(candidate),
    ) if value)


def _has_database_analysis_complex(item: dict) -> bool:
    """Recognize stored complexes and source-linked target-conditioned poses.

    A few design imports register the binder as the molecular identity while
    keeping the target only in the exact saved complex representation and its
    source metadata. The lightweight Database catalog marks that explicit
    provenance without hydrating every record. The selected structure is
    still parsed before the user chooses chains or launches an evaluator.
    """
    if len(item.get("components") or []) >= 2:
        return True
    if any((structure.get("context") or {}).get("database_structure_policy") == "bound_docking_complex"
           for structure in item.get("representations") or item.get("structures") or []):
        return True
    if item.get("source_target_provenance") is True:
        return True
    metadata = item.get("metadata") if isinstance(item.get("metadata"), dict) else {}
    return bool(metadata.get("target_chain_inference") and (metadata.get("target_fasta") or metadata.get("target_pdb")))


def _target_structure_notice(metadata: dict) -> str:
    """Describe sequence-defined target provenance without inferring a fold."""
    validation = metadata.get("validation_acceptance") if isinstance(metadata.get("validation_acceptance"), dict) else {}
    target_type = str(validation.get("target_type") or metadata.get("target_type") or "").strip().casefold()
    if target_type == "fasta_idr":
        return (
            "The source workflow identifies this as a sequence-defined IDR target. Structural interface, confidence, "
            "and energy measurements involving that chain describe this saved predicted conformation; they do not "
            "establish a stable native target fold."
        )
    if target_type in {"sequence", "fasta"}:
        return (
            "The source workflow defines this target by sequence. Structural measurements involving that chain "
            "describe this saved predicted conformation; sequence provenance does not establish a native fold or "
            "a full-length structural reference."
        )
    return ""


def _evaluation_matches_structure(evaluation: dict, structure: dict) -> bool:
    context = structure.get("context") if isinstance(structure.get("context"), dict) else {}
    representation_id = str(evaluation.get("representation_id") or "")
    if representation_id:
        return representation_id in {
            str(structure.get("id") or ""),
            str(context.get("representation_id") or ""),
        }
    representation = str(evaluation.get("representation") or "")
    if representation:
        if representation in {str(structure.get("path") or ""), str(context.get("structure_digest") or "")}:
            return True
        # Older Database evaluations may store a digest/path in the descriptor
        # while their canonical association is the run and candidate context.
    run_id = str(evaluation.get("run_id") or "")
    candidate_id = str(evaluation.get("candidate_id") or "")
    if run_id and str(context.get("run_id") or "") == run_id:
        return not candidate_id or str(context.get("candidate_id") or "") == candidate_id
    return False


def _representation_label(item: dict) -> str:
    context = item.get("context") if isinstance(item.get("context"), dict) else {}
    model = context.get("model_index") or context.get("model") or context.get("sample") or ""
    candidate = context.get("candidate_id") or ""
    prefix = f"Model {model} · " if model else ""
    prefix += f"{candidate} · " if candidate else ""
    return prefix + str(item.get("label") or Path(str(item.get("path") or "structure")).name)


def _accessible_result_jobs(user: dict) -> list[dict]:
    from mn_molecular_design_studio.auth import accessible_project_ids
    from mn_molecular_design_studio.job_catalog import search_jobs

    rows = search_jobs(
        statuses=["completed", "succeeded"],
        access_owner_id=str(user.get("id") or ""),
        accessible_project_ids=accessible_project_ids(user),
        access_admin=user.get("role") == "admin",
        # Older completed Results are still valid analysis sources. Keep the
        # contextual picker aligned with the job catalog's one-page maximum so
        # existing runs remain reachable without making a separate registry.
        limit=1000,
        sort_by="finished_desc",
    )["items"]
    return [row for row in rows if str(row.get("status") or "").casefold() in {"completed", "succeeded"} and row.get("success") is not False]


def _job_structure_options(user: dict, job: dict) -> list[dict]:
    from mn_molecular_design_studio.auth import can_access_job
    from mn_molecular_design_studio.job_catalog import get_job_access_record, job_artifacts
    from mn_molecular_design_studio.runtime import runs_root

    task_group = str(job.get("task_group") or "")
    run_id = str(job.get("run_id") or "")
    if not task_group or not run_id or not can_access_job(user, get_job_access_record(task_group, run_id)):
        return []
    run_dir = Path(str(job.get("run_dir") or "")).expanduser().resolve()
    root = Path(runs_root()).expanduser().resolve()
    candidates: list[dict] = []
    candidate_path = run_dir / "artifacts" / "normalized_candidates" / "candidates.jsonl"
    if candidate_path.is_file():
        try:
            from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
            candidates = [row for row in read_candidates(candidate_path) if isinstance(row, dict)]
        except (OSError, ValueError, TypeError):
            candidates = []
    from mn_molecular_design_studio import cofolding_results
    candidate_by_path = {}
    for candidate in candidates:
        structure = cofolding_results.structure_for(run_dir, candidate)
        if structure is not None:
            try:
                candidate_by_path[str(structure.resolve())] = candidate
            except OSError:
                pass

    rows = []
    for artifact in job_artifacts(task_group, run_id):
        relative = str(artifact.get("relative_path") or "")
        try:
            path = (root / relative).resolve(strict=True)
            path.relative_to(run_dir)
            path.relative_to(root)
        except (KeyError, OSError, ValueError):
            continue
        if not path.is_file() or _structure_suffix(path) not in STRUCTURE_SUFFIXES:
            continue
        try:
            candidate = candidate_by_path.get(str(path.resolve()), {})
            rows.append({
                "path": path,
                "relative_path": path.relative_to(run_dir).as_posix(),
                "artifact_path": relative,
                "artifact_task_group": task_group,
                "artifact_run_id": run_id,
                "candidate": candidate,
                "candidate_id": str(candidate.get("candidate_id") or ""),
                "label": str(artifact.get("label") or path.name),
            })
        except (OSError, ValueError):
            continue
    # A normalized candidate can refer to a saved structure which is absent
    # from the source job's artifact catalog. For connected outputs in another
    # run, require the producing job ACL and its exact registered artifact.
    known = {str(item["path"]) for item in rows}
    for candidate in candidates:
        structure = cofolding_results.structure_for(run_dir, candidate)
        if structure is None:
            continue
        try:
            resolved = structure.resolve(strict=True)
            producer_task_group, producer_run_id, relative_path = _registered_cross_run_structure(
                user, resolved, source_run_dir=run_dir, runs_root=root)
        except (OSError, ValueError, PermissionError):
            continue
        if str(resolved) in known or _structure_suffix(resolved) not in STRUCTURE_SUFFIXES:
            continue
        relative_to_root = resolved.relative_to(root)
        rows.append({
            "path": resolved,
            "relative_path": relative_path,
            "artifact_path": relative_to_root.as_posix(),
            "artifact_task_group": producer_task_group,
            "artifact_run_id": producer_run_id,
            "candidate": candidate,
            "candidate_id": str(candidate.get("candidate_id") or ""),
            "label": resolved.name,
        })
        known.add(str(resolved))
    rows.sort(key=lambda item: (item["relative_path"].casefold(), item["candidate_id"].casefold()))
    return rows


def _results_review_reference(user: dict, job: dict, selected: dict) -> tuple[str, str]:
    """Address a Results review to the authorized job that owns the structure."""
    from mn_molecular_design_studio.runtime import runs_root

    root = Path(runs_root()).expanduser().resolve()
    source_run_dir = Path(str(job.get("run_dir") or "")).expanduser().resolve(strict=True)
    structure_path = Path(str(selected.get("path") or "")).expanduser().resolve(strict=True)
    producer_group, producer_run_id, _ = _registered_cross_run_structure(
        user, structure_path, source_run_dir=source_run_dir, runs_root=root)
    producer_run_dir = (root / producer_group / producer_run_id).resolve(strict=True)
    relative_path = structure_path.relative_to(producer_run_dir).as_posix()
    return f"{producer_group}/{producer_run_id}", relative_path


def _analysis_source_reference(job: dict) -> str:
    """Use the Studio scheduler's authorized peer-job form when applicable."""
    cluster_job_id = str(job.get("cluster_job_id") or "")
    if cluster_job_id and job.get("node_id") and job.get("local_reference"):
        return "cluster:" + cluster_job_id
    return f"{job.get('task_group') or ''}/{job.get('run_id') or ''}"


def _authorized_remote_source_project(job: dict) -> dict | None:
    """Resolve a remote source project's name from the coordinator's authorized organization view."""
    from mn_molecular_design_studio.cluster import ui as cluster_ui

    project_id = str(job.get("project_id") or "")
    job_id = str(job.get("cluster_job_id") or "")
    if not project_id or not job_id:
        return None
    token = str(st.session_state.get("studio-session-token") or "development-mode")
    organization = cluster_ui._request("/organization", token)
    projects = organization.get("projects") if isinstance(organization, dict) else None
    return next((dict(project) for project in projects or []
                 if isinstance(project, dict) and str(project.get("id") or "") == project_id), None)


def _analysis_structure_reference(job: dict, selected: dict, path: Path) -> str:
    if job.get("cluster_job_id") and selected.get("analysis_structure_uri"):
        return str(selected["analysis_structure_uri"])
    return str(path)


def _authorized_remote_structure_bytes(user: dict, opened: dict) -> tuple[dict, bytes] | None:
    """Reauthorize a peer job and fetch only its selected registered artifact."""
    from mn_molecular_design_studio.cluster import ui as cluster_ui

    job_id = str(opened.get("cluster_job_id") or "")
    if not job_id:
        raise ValueError("The remote Results job reference is missing.")
    token = str(st.session_state.get("studio-session-token") or "development-mode")
    remote_job = cluster_ui._request("/jobs/" + job_id, token)
    expected = {
        "node_id": str(opened.get("node_id") or ""),
        "local_reference": str(opened.get("local_reference") or ""),
        "project_id": str(opened.get("project_id") or ""),
    }
    if any(str(remote_job.get(key) or "") != value for key, value in expected.items()):
        raise PermissionError("The remote Results job changed or is no longer accessible.")
    relative = cluster_ui._safe_run_artifact_path(str(opened.get("artifact_path") or ""))
    if not relative:
        raise ValueError("The selected structure is not a valid run-relative artifact of this remote job.")
    reference = PurePosixPath(str(remote_job.get("local_reference") or ""))
    expected_uri = f"cluster://{expected['node_id']}/" + (reference / PurePosixPath(relative)).as_posix()
    if (reference.is_absolute() or len(reference.parts) != 2 or ".." in reference.parts
            or str(opened.get("analysis_structure_uri") or "") != expected_uri):
        raise ValueError("The selected structure reference does not match its authorized source job.")
    key = "complex-analysis-remote-artifact-" + hashlib.sha256(f"{job_id}:{relative}".encode()).hexdigest()[:20]
    content = st.session_state.get(key)
    if not isinstance(content, bytes):
        command_id = st.session_state.get(key + "-command")
        if not command_id:
            operation = cluster_ui._request(
                "/jobs/" + job_id + "/actions", token,
                {"action": "artifact", "parameters": {"path": relative}}, method="POST",
            )
            command_id = str(operation.get("command_id") or "")
            if not command_id:
                raise ValueError("The source workstation did not accept the artifact request.")
            st.session_state[key + "-command"] = command_id
        status = cluster_ui._request("/commands/" + command_id, token)
        if status.get("status") == "failed":
            raise ValueError(str(status.get("error") or "The source workstation could not serve this structure."))
        if status.get("status") != "completed":
            st.info("Waiting for the source workstation to return this structure. Use Refresh to check again.")
            return None
        metadata = status.get("response") if isinstance(status.get("response"), dict) else {}
        size = int(metadata.get("size_bytes") or 0)
        if size > VIEWER_MAX_BYTES:
            raise ValueError("The selected remote structure exceeds the 50 MB Scientific Inspector viewing limit.")
        content = cluster_ui._request("/commands/" + command_id + "/file", token, raw=True, timeout=60)
        if not isinstance(content, bytes):
            raise ValueError("The source workstation returned an invalid structure payload.")
        st.session_state[key] = content
    if len(content) > VIEWER_MAX_BYTES:
        raise ValueError("The selected remote structure exceeds the 50 MB Scientific Inspector viewing limit.")
    expected_digest = str(opened.get("structure_sha256") or "")
    digest = hashlib.sha256(content).hexdigest()
    if expected_digest and digest != expected_digest:
        raise ValueError("The remote structure changed since it was selected. Reopen the current Results artifact.")
    return remote_job, content


def _render_remote_results_source(user: dict, opened: dict, *, workflow_schema=None,
                                  project_selector=None, submit=None) -> None:
    """Render one peer-owned candidate in the existing structure Inspector."""
    fetched = _authorized_remote_structure_bytes(user, opened)
    if fetched is None:
        if st.button("Refresh remote structure", key="complex-analysis-remote-refresh"):
            st.rerun()
        return
    remote_job, content = fetched
    candidate = opened.get("candidate") if isinstance(opened.get("candidate"), dict) else {}
    if str(candidate.get("candidate_id") or "") != str(opened.get("candidate_id") or ""):
        raise ValueError("The selected structure no longer matches its source candidate.")
    relative = str(opened.get("artifact_path") or "")
    digest = hashlib.sha256(content).hexdigest()
    path = Path(str(opened.get("staged_structure") or ""))
    if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != digest:
        path = _stage_remote_result_file(content, _structure_suffix(Path(relative)), user=user)
    candidate_file = Path(str(opened.get("candidate_file") or ""))
    if not candidate_file.is_file():
        candidate_file = _remote_candidate_manifest(candidate, path, user=user)
    staged_run_dir = candidate_file.parents[2]
    local_reference = str(remote_job.get("local_reference") or "")
    from pathlib import PurePosixPath
    uri = "cluster://" + str(remote_job.get("node_id")) + "/" + (
        PurePosixPath(local_reference) / PurePosixPath(relative)
    ).as_posix()
    if uri != str(opened.get("analysis_structure_uri") or ""):
        raise ValueError("The peer artifact reference changed; reopen this candidate from Results.")

    raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
    components = component_inventory(content, path.name, raw.get("components"))
    node_name = str(remote_job.get("node_name") or "source workstation")
    st.caption(f"Source Results: {node_name} · {remote_job.get('task_group')}/{remote_job.get('run_id')} · {relative}")
    st.caption(f"Candidate/model: {candidate.get('candidate_id')} · {opened.get('model_index') or 'model not specified'}")
    st.caption(f"SHA-256: {digest}")
    _show_components(components)
    key = "complex-analysis-remote-viewer-" + hashlib.sha256(f"{opened.get('cluster_job_id')}:{relative}".encode()).hexdigest()[:14]
    _show_structure(path, raw, components, key=key,
                    highlighted_selections=_active_biotite_highlights(key))
    with st.expander("Remote Results provenance", expanded=False):
        st.json({
            "source_workstation": node_name,
            "node_id": remote_job.get("node_id"),
            "cluster_job_id": opened.get("cluster_job_id"),
            "task_group": remote_job.get("task_group"), "run_id": remote_job.get("run_id"),
            "project_id": remote_job.get("project_id"),
            "campaign_name": remote_job.get("campaign_name"),
            "candidate_id": candidate.get("candidate_id"), "model_id": opened.get("model_index"),
            "artifact_uri": uri, "structure_sha256": digest,
            "candidate_source_metadata": candidate.get("raw_metadata") or {},
        })
    st.caption("Read-only Results inspection. To analyze or shortlist this design, explicitly save the selected candidate/model to the Molecular Database from its Results page.")


def _read_json(path: Path) -> dict:
    try:
        value = json.loads(path.read_text())
    except (OSError, ValueError):
        return {}
    return value if isinstance(value, dict) else {}


def _structure_digest(path: Path) -> str:
    try:
        with path.open("rb") as stream:
            return hashlib.file_digest(stream, "sha256").hexdigest()
    except OSError:
        return ""


def _protein_chain_ids(components: list[dict]) -> list[str]:
    return sorted({
        str(item.get("source_chain") or item.get("chain") or "")
        for item in components
        if str(item.get("type") or "").casefold() == "protein"
        and str(item.get("source_chain") or item.get("chain") or "")
    })


def _ligand_selection_for_component(component: dict, ligand_components: list[dict], metadata: dict) -> dict:
    """Build an exact, pose-local selector with an explicit isomeric identity."""
    residue_key = str(component.get("residue_key") or "")
    # Imported component.chain may include an internal residue suffix (for
    # example ``A_201``); residue_key retains the author chain separately.
    chain = str(component.get("source_chain") or residue_key.split(":", 1)[0]
                or component.get("chain") or "")
    ccd = str(component.get("ccd") or "").upper()
    smiles = str(component.get("smiles") or "").strip()
    saved_selection = metadata.get("ligand_selection")
    if isinstance(saved_selection, dict):
        saved_chain = str(saved_selection.get("chain") or "")
        saved_residue = str(saved_selection.get("residue_key") or "")
        if (saved_chain == chain and saved_residue == residue_key
                and str(saved_selection.get("smiles") or "").strip()):
            return {
                "chain": saved_chain, "residue_key": saved_residue,
                "ccd": str(saved_selection.get("ccd") or ccd).upper(),
                "smiles": str(saved_selection["smiles"]).strip(),
            }
    if not smiles:
        stored = next((item for item in (metadata.get("components") or [])
                       if isinstance(item, dict) and str(item.get("residue_key") or "") == residue_key
                       and str(item.get("type") or "") == "small_molecule"), {})
        smiles = str(stored.get("smiles") or "").strip()
        ccd = ccd or str(stored.get("ccd") or "").upper()
    if not smiles:
        values = metadata.get("ligand_smiles")
        if isinstance(values, str) and len(ligand_components) == 1:
            smiles = values.strip()
        elif isinstance(values, list) and len(values) == len(ligand_components):
            value = values[ligand_components.index(component)]
            smiles = str(value.get("smiles") or value.get("value") or "").strip() if isinstance(value, dict) else str(value or "").strip()
    return {"chain": chain, "residue_key": residue_key, "ccd": ccd, "smiles": smiles}


def _biotite_item_matches_source(item: dict, outputs: dict[str, dict], *, candidate_id: str,
                                structure_digest: str, binder_chains: list[str],
                                target_chains: list[str], source_run_id: str) -> tuple[bool, dict]:
    """Match a tool input to native or CandidateSet identity by its exact receipt."""
    if (str(item.get("structure_digest") or "") != structure_digest
            or sorted(map(str, item.get("binder_chains") or [])) != sorted(map(str, binder_chains))
            or sorted(map(str, item.get("target_chains") or [])) != sorted(map(str, target_chains))):
        return False, {}
    if str(item.get("candidate_id") or "") == candidate_id:
        return True, {}
    output = outputs.get(str(item.get("candidate_id") or ""), {})
    raw = output.get("raw_metadata") if isinstance(output.get("raw_metadata"), dict) else {}
    receipt = raw.get("complex_analysis_entry") if isinstance(raw.get("complex_analysis_entry"), dict) else {}
    receipt_candidate = str(receipt.get("native_candidate_id") or receipt.get("candidate_id") or "")
    receipt_run = str(receipt.get("source_run") or Path(str(receipt.get("source_run_dir") or "")).name)
    matches = (receipt_candidate == candidate_id and receipt_run == source_run_id
               and str(receipt.get("structure_sha256") or "") == structure_digest)
    return matches, receipt if matches else {}


def _results_analysis_rows(user, tool, source_task_group, source_run_id):
    from mn_molecular_design_studio.auth import accessible_project_ids
    from mn_molecular_design_studio.job_catalog import search_jobs
    rows, cursor = [], None
    while True:
        page = search_jobs(tools=[tool], task_group="analysis", limit=500, cursor=cursor,
                           upstream_task_group=source_task_group, upstream_run_id=source_run_id,
                           access_owner_id=str(user.get("id") or ""),
                           accessible_project_ids=accessible_project_ids(user),
                           access_admin=user.get("role") == "admin", sort_by="active_desc")
        rows.extend(page["items"])
        next_cursor = page.get("next_cursor")
        if next_cursor is None or next_cursor == cursor:
            return rows
        cursor = next_cursor


def _saved_analysis_protocol(run_dir, settings):
    from mn_molecular_design_studio.scientific_settings import PROTOCOL_FIELDS
    params = _read_json(Path(run_dir) / "input.json").get("params") or {}
    command = _read_json(Path(run_dir) / "command.json")
    return {**{key: value for key, value in params.items() if key in PROTOCOL_FIELDS},
            **(params.get("assessment_protocol") or {}),
            **{key: value for key, value in settings.items() if key in PROTOCOL_FIELDS},
            **{key: value for key, value in command.items() if key in {"image", "image_digest"}}}


def _matching_biotite_jobs(user: dict, *, source_task_group: str, source_run_id: str,
                           candidate_id: str, structure_digest: str,
                           binder_chains: list[str], target_chains: list[str],
                           contact_distance: float, model_id: str = "") -> list[dict]:
    """Find reusable Biotite evaluations for this exact source pose and chain mapping."""
    from mn_molecular_design_studio.auth import accessible_project_ids, can_access_job
    from mn_molecular_design_studio.job_catalog import get_job_access_record, search_jobs

    rows = _results_analysis_rows(user, "biotite_interface", source_task_group, source_run_id)
    matches = []
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
    for row in rows:
        task_group = str(row.get("task_group") or "")
        run_id = str(row.get("run_id") or "")
        if not task_group or not run_id or not can_access_job(user, get_job_access_record(task_group, run_id)):
            continue
        run_dir = Path(str(row.get("run_dir") or "")).expanduser().resolve()
        metadata = _read_json(run_dir / "metadata.json")
        if (str(metadata.get("upstream_task_group") or "") != source_task_group
                or str(metadata.get("upstream_run_id") or "") != source_run_id):
            continue
        manifest = _read_json(run_dir / "artifacts" / "analysis" / "batch_manifest.json")
        if manifest.get("tool") != "biotite_interface":
            continue
        settings = manifest.get("settings") if isinstance(manifest.get("settings"), dict) else {}
        try:
            saved_distance = float(settings.get("contact_distance", 5.0))
        except (TypeError, ValueError):
            continue
        if saved_distance != float(contact_distance):
            continue
        selection = manifest.get("selection") if isinstance(manifest.get("selection"), dict) else {}
        desired_binder, desired_target = sorted(map(str, binder_chains)), sorted(map(str, target_chains))
        output_path = run_dir / "artifacts" / "normalized_candidates" / "candidates.jsonl"
        outputs = {str(candidate.get("candidate_id") or ""): candidate
                   for candidate in read_candidates(output_path)} if output_path.is_file() else {}
        selected_item = None
        matched_receipt = {}
        for item in manifest.get("items", []):
            identity_matches, receipt = _biotite_item_matches_source(
                item, outputs, candidate_id=candidate_id, structure_digest=structure_digest,
                binder_chains=desired_binder, target_chains=desired_target, source_run_id=source_run_id)
            if identity_matches:
                selected_item, matched_receipt = item, receipt
                break
        if selected_item is None:
            continue
        saved_model_id = str(selected_item.get("model_id") or matched_receipt.get("model_id")
                              or selection.get("model_id") or "")
        if saved_model_id != model_id:
            continue
        status = str(row.get("status") or metadata.get("status") or "unknown").casefold()
        from mn_molecular_design_studio.complex_analysis_batch import assessment_protocol, protocol_matches
        if status in {"completed", "succeeded"} and not protocol_matches(
                _saved_analysis_protocol(run_dir, settings), assessment_protocol("analysis.biotite_interface")):
            continue
        result = _read_json(run_dir / "artifacts" / "analysis" / "batch_status.json")
        item_status = next((item for item in result.get("items", []) if item.get("id") == selected_item.get("id")), {})
        error = str(item_status.get("error") or metadata.get("error") or "")
        selected_status = str(item_status.get("status") or "").casefold()
        if status in {"completed", "succeeded"} and selected_status and selected_status not in {"completed", "succeeded"}:
            status = selected_status
        matches.append({**row, "status": status, "run_dir": run_dir,
                        "manifest": manifest, "metadata": metadata,
                        "selected_item": selected_item, "error": error})
    # Prefer a successful exact evaluation. A failed retry does not make an
    # earlier successful matching calculation stale.
    rank = {"completed": 0, "succeeded": 0, "running": 1, "starting": 1, "queued": 1,
            "preparing": 1, "pending": 1, "failed": 2, "cancelled": 3, "canceled": 3}
    return sorted(matches, key=lambda item: (rank.get(item["status"], 4), str(item.get("created_at") or "")), reverse=False)


def _analysis_evaluation_id(run_dir: Path, user: dict, source_candidate_id: str, *,
                            evaluator: str = "biotite_interface") -> str:
    """Resolve the persistent Evaluation descriptor created by normal Results ingestion."""
    try:
        from mn_molecular_design_studio.evaluations import evaluation_descriptor
        from mn_molecular_design_studio.molecular_database_results import list_run_molecules
        records = list_run_molecules(run_dir, user)
    except (ImportError, OSError, PermissionError, ValueError):
        return ""
    for record in records:
        for observation in record.get("observations", []):
            try:
                descriptor = evaluation_descriptor(observation)
            except (KeyError, TypeError, ValueError):
                continue
            candidate_context = (observation.get("context") or {}).get("candidate")
            candidate = candidate_context if isinstance(candidate_context, dict) else {}
            candidate_id = str(descriptor.get("candidate_id") or candidate.get("candidate_id") or "")
            source_id = str(candidate.get("candidate_id") or "")
            candidate_matches = candidate_id.startswith(source_candidate_id + "_") or source_id == source_candidate_id
            if descriptor.get("evaluator") == evaluator and candidate_matches:
                return str(descriptor.get("id") or "")
    return ""


def _matching_contextual_jobs(user: dict, *, workflow_id: str, source_task_group: str,
                              source_run_id: str, candidate_id: str, structure_digest: str,
                              model_id: str, binder_chains: list[str], target_chains: list[str],
                              parameters: dict, ligand_selection: dict | None = None) -> list[dict]:
    """Find an exact reusable evaluation for one selected source pose."""
    from mn_molecular_design_studio.scientific_settings import scientific_parameters
    from mn_molecular_design_studio.auth import accessible_project_ids, can_access_job
    from mn_molecular_design_studio.job_catalog import get_job_access_record, search_jobs

    tool = {"analysis.ipsae": "ipsae", "analysis.rosetta": "rosetta_metrics", "analysis.pep_patch": "pep_patch",
            "analysis.plip": "plip", "analysis.pandamap": "pandamap", "analysis.posebusters": "posebusters"}.get(workflow_id)
    if not tool:
        return []
    rows = _results_analysis_rows(user, tool, source_task_group, source_run_id)
    desired_binder = sorted(map(str, binder_chains))
    desired_target = sorted(map(str, target_chains))
    matches = []
    for row in rows:
        task_group = str(row.get("task_group") or "")
        run_id = str(row.get("run_id") or "")
        if not task_group or not run_id or not can_access_job(user, get_job_access_record(task_group, run_id)):
            continue
        run_dir = Path(str(row.get("run_dir") or "")).expanduser().resolve()
        metadata = _read_json(run_dir / "metadata.json")
        if (str(metadata.get("upstream_task_group") or "") != source_task_group
                or str(metadata.get("upstream_run_id") or "") != source_run_id):
            continue
        if workflow_id in {"analysis.plip", "analysis.pandamap", "analysis.posebusters"}:
            manifest = _read_json(run_dir / "artifacts" / "analysis" / "batch_manifest.json")
            selection = manifest.get("selection") if isinstance(manifest.get("selection"), dict) else {}
            item = next((candidate for candidate in manifest.get("items", [])
                         if str(candidate.get("candidate_id") or "") == candidate_id
                         and str(candidate.get("source_structure_digest") or "") == structure_digest
                         and str(candidate.get("model_id") or (selection.get("model_ids") or {}).get(candidate_id)
                                 or selection.get("model_id") or "") == model_id
                         and candidate.get("ligand_selection") == (ligand_selection or {})), None)
            saved_parameters = {key: value for key, value in (manifest.get("settings") or {}).items()
                                if key not in {"workers", "protocol_revision", "checkpoint_sha256"}}
            saved_model_id = str(selection.get("model_id") or "")
        elif workflow_id == "analysis.pep_patch":
            manifest = _read_json(run_dir / "artifacts" / "analysis" / "batch_manifest.json")
            selection = manifest.get("selection") if isinstance(manifest.get("selection"), dict) else {}
            item = next((candidate for candidate in manifest.get("items", [])
                         if str(candidate.get("candidate_id") or "") == candidate_id
                         and str(candidate.get("structure_digest") or "") == structure_digest
                         and sorted(map(str, candidate.get("binder_chains") or [])) == desired_binder
                         and not desired_target), None)
            saved_parameters = {key: value for key, value in (manifest.get("settings") or {}).items()
                                if key not in {"workers", "protocol_revision", "checkpoint_sha256"}}
            saved_model_id = str(selection.get("model_id") or "")
        else:
            manifest = _read_json(run_dir / "artifacts" / "analysis" / "complex_analysis_selection.json")
            item = manifest if manifest.get("schema") == "studio.complex_analysis_selection.v1" else None
            saved_parameters = manifest.get("parameters") if isinstance(manifest.get("parameters"), dict) else {}
            saved_model_id = str(manifest.get("model_id") or "")
        if not item:
            continue
        saved_model_id = str(item.get("model_id") or (selection.get("model_ids") or {}).get(candidate_id)
                             or saved_model_id or "")
        item_digest = (str(item.get("source_structure_digest") or "")
                       if workflow_id in {"analysis.plip", "analysis.pandamap", "analysis.posebusters"}
                       else str(item.get("structure_digest") or ""))
        if (str(item.get("candidate_id") or "") != candidate_id
                or item_digest != structure_digest
                or saved_model_id != model_id
                or (workflow_id not in {"analysis.plip", "analysis.pandamap", "analysis.posebusters"}
                    and sorted(map(str, item.get("binder_chains") or [])) != desired_binder)
                or (workflow_id not in {"analysis.pep_patch", "analysis.plip", "analysis.pandamap", "analysis.posebusters"}
                    and sorted(map(str, item.get("target_chains") or [])) != desired_target)
                or scientific_parameters(saved_parameters, workflow_id) != scientific_parameters(parameters, workflow_id)
                or (workflow_id in {"analysis.plip", "analysis.pandamap", "analysis.posebusters"}
                    and ((selection.get("ligand_selections") or {}).get(candidate_id)
                         or selection.get("ligand_selection") or item.get("ligand_selection") or {})
                    != (ligand_selection or {}))):
            continue
        status = str(row.get("status") or metadata.get("status") or "unknown").casefold()
        from mn_molecular_design_studio.complex_analysis_batch import assessment_protocol, protocol_matches
        if status in {"completed", "succeeded"} and not protocol_matches(
                _saved_analysis_protocol(run_dir, manifest.get("settings") or manifest.get("parameters") or {}),
                assessment_protocol(workflow_id)):
            continue
        error = str(metadata.get("error") or "")
        if workflow_id in {"analysis.pep_patch", "analysis.plip", "analysis.pandamap", "analysis.posebusters"} and status in {"completed", "succeeded"}:
            batch = _read_json(run_dir / "artifacts" / "analysis" / "batch_status.json")
            batch_item = next((candidate for candidate in batch.get("items", []) if candidate.get("id") == item.get("id")), {})
            item_status = str(batch_item.get("status") or "").casefold()
            error = str(batch_item.get("error") or error)
            if item_status and item_status not in {"completed", "succeeded"}:
                status = item_status
        matches.append({**row, "status": status, "run_dir": run_dir,
                        "manifest": manifest, "metadata": metadata, "selected_item": item, "error": error})
    rank = {"completed": 0, "succeeded": 0, "running": 1, "starting": 1, "queued": 1,
            "preparing": 1, "pending": 1, "failed": 2, "cancelled": 3, "canceled": 3}
    return sorted(matches, key=lambda item: (rank.get(item["status"], 4), str(item.get("created_at") or "")))


def _biotite_contact_residues(run_dir: Path, candidate_id: str) -> tuple[list[str], list[str], list[str]]:
    """Return saved author-chain residue selectors from this exact batch item."""
    manifest = _read_json(run_dir / "artifacts" / "analysis" / "batch_manifest.json")
    item = next((row for row in manifest.get("items", [])
                 if str(row.get("candidate_id") or "") == candidate_id), None)
    if not item:
        return [], [], []
    contacts = _read_json(run_dir / "artifacts" / "analysis" / str(item.get("id") or "") / "contacts.json")
    binder, target = set(), set()
    for pair in contacts.get("residue_pairs", []):
        if not isinstance(pair, dict):
            continue
        for destination, field in ((binder, "binder_residue"), (target, "target_residue")):
            selector = str(pair.get(field) or "").strip()
            # Keep unusual identifiers in the map, but pass only selectors
            # supported by Mol*'s current auth-chain/sequence selection syntax.
            if re.fullmatch(r"[^:\s]+:-?\d+", selector):
                destination.add(selector)
    binder_rows, target_rows = sorted(binder), sorted(target)
    return binder_rows, target_rows, sorted(binder | target)


def _active_biotite_highlights(viewer_key: str) -> list[str]:
    residue_sets = st.session_state.get(f"{viewer_key}-biotite-residue-sets", {})
    mode = st.session_state.get(f"{viewer_key}-biotite-highlight-mode", "No highlight")
    selected = residue_sets.get(mode, []) if isinstance(residue_sets, dict) else []
    selected = list(selected) if isinstance(selected, list) else []
    inspector = st.session_state.get(f"{viewer_key}-inspector-highlights", [])
    if isinstance(inspector, list):
        selected.extend(item for item in inspector if isinstance(item, str))
    return list(dict.fromkeys(selected))


def _select_inspector_evidence(state_key: str, selectors: list[str]) -> None:
    st.session_state[state_key] = list(dict.fromkeys(item for item in selectors if item))


def _authorized_evaluation_run_dir(user: dict, observation: dict) -> Path | None:
    """Resolve a saved Evaluation's artifact directory through the job ACL."""
    from mn_molecular_design_studio.auth import can_access_job
    from mn_molecular_design_studio.complex_analysis_inspector import observation_context
    from mn_molecular_design_studio.job_catalog import get_job_access_record
    from mn_molecular_design_studio.runtime import runs_root

    context = observation_context(observation)
    run_id = str(observation.get("run_id") or context.get("run_id") or "")
    task_group = str(context.get("task_group") or "analysis")
    if not run_id:
        return None
    job_access = get_job_access_record(task_group, run_id)
    if not can_access_job(user, job_access):
        return None
    relative_path = str(job_access.get("run_relative_path") or "") if isinstance(job_access, dict) else ""
    if not relative_path:
        return None
    root = Path(runs_root()).expanduser().resolve()
    relative = Path(relative_path)
    if relative.is_absolute():
        return None
    path = (root / relative).resolve()
    try:
        path.relative_to(root)
    except ValueError:
        return None
    if not can_access_job(user, job_access):
        return None
    return path if path.is_dir() else None


def _csv_rows(path: Path) -> list[dict]:
    import csv

    try:
        with path.open("r", newline="", encoding="utf-8-sig") as handle:
            return [dict(row) for row in csv.DictReader(handle)]
    except (OSError, UnicodeError, csv.Error):
        return []


def _saved_evaluation_artifacts(user: dict, observation: dict, structure: dict,
                                current_digest: str) -> list[dict]:
    """Read only artifacts belonging to an accessible, exact-structure Evaluation."""
    from mn_molecular_design_studio.complex_analysis_inspector import observation_context

    context = observation_context(observation)
    candidate = context.get("candidate") if isinstance(context.get("candidate"), dict) else {}
    raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
    selection = raw.get("complex_analysis_selection") if isinstance(raw.get("complex_analysis_selection"), dict) else {}
    representation_id = str(context.get("representation_id") or "")
    selection_source = selection.get("database_source") if isinstance(selection.get("database_source"), dict) else {}
    source_structure_id = str(context.get("source_structure_id") or selection_source.get("structure_id") or "")
    saved_digest = str(context.get("source_structure_digest") or context.get("structure_digest")
                       or selection.get("structure_digest") or "")
    structure_id = str(structure.get("id") or "")
    if saved_digest and saved_digest != current_digest:
        return []
    if not (representation_id == structure_id or source_structure_id == structure_id
            or saved_digest == current_digest):
        return []
    run_dir = _authorized_evaluation_run_dir(user, observation)
    if not run_dir:
        return []
    analysis_dir = run_dir / "artifacts" / "analysis"
    manifest = _read_json(analysis_dir / "batch_manifest.json")
    candidate_id = str(context.get("candidate_id") or candidate.get("candidate_id") or "")
    item = next((row for row in manifest.get("items", [])
                 if candidate_id and str(row.get("candidate_id") or "") == candidate_id), None)
    if item is None:
        candidate_dirs = sorted(analysis_dir.glob("item_*/"))
        if len(candidate_dirs) == 1:
            item = {"id": candidate_dirs[0].name}
    item_id = str((item or {}).get("id") or "")
    if item_id and Path(item_id).name != item_id:
        return []
    item_dir = analysis_dir / item_id if item_id else analysis_dir
    results = []
    for name in ("tag_fusion_suitability.json", "sequence_liability.json", "binder_interface_chemistry.json", "contacts.json", "interactions.csv",
                 "full_checks.csv", "report.json", "patch_metrics.json"):
        path = item_dir / name
        if not path.is_file():
            continue
        if name.endswith(".json"):
            payload = _read_json(path)
            if name == "contacts.json":
                rows = payload.get("residue_pairs", []) if isinstance(payload, dict) else []
                salt_pairs = payload.get("salt_bridge_pairs", []) if isinstance(payload, dict) else []
                results.append({"name": name, "path": path,
                                "rows": [dict(row) for row in rows if isinstance(row, dict)]})
                if salt_pairs:
                    results.append({"name": "salt_bridges.json", "path": path,
                                    "rows": [{"binder_residue": pair[0], "target_residue": pair[1],
                                              "interaction_type": "Salt bridge"}
                                             for pair in salt_pairs if isinstance(pair, (list, tuple)) and len(pair) == 2]})
                continue
            elif name == "reference_environment.json":
                results.append({"name": name, "path": path,
                                "rows": payload.get("residue_evidence", []) if isinstance(payload, dict) else [],
                                "payload": payload})
                continue
            elif name == "sequence_liability.json":
                results.append({"name": name, "path": path,
                                "rows": payload.get("findings", []) if isinstance(payload, dict) else [],
                                "payload": payload})
                continue
            elif name == "binder_interface_chemistry.json":
                results.append({"name": name, "path": path,
                                "rows": payload.get("residue_annotations", []) if isinstance(payload, dict) else [],
                                "payload": payload})
                continue
            else:
                rows = [payload] if payload else []
        else:
            rows = _csv_rows(path)
        results.append({"name": name, "path": path, "rows": [dict(row) for row in rows if isinstance(row, dict)]})
    environment_path = run_dir / "artifacts" / "analysis" / "reference_environment" / "reference_environment.json"
    if environment_path.is_file():
        payload = _read_json(environment_path)
        if isinstance(payload, dict):
            results.append({"name": "reference_environment.json", "path": environment_path,
                            "rows": payload.get("residue_evidence", []), "payload": payload})
    for path in sorted(item_dir.glob("residue_patches*.csv")):
        rows = _csv_rows(path)
        if rows:
            results.append({"name": path.name, "path": path, "rows": rows})
    liability_path = analysis_dir / "sequence_liability" / "sequence_liability.json"
    if liability_path.is_file():
        payload = _read_json(liability_path)
        if isinstance(payload, dict):
            results.append({"name": "sequence_liability.json", "path": liability_path,
                            "rows": payload.get("findings", []), "payload": payload})
    interface_chemistry_path = analysis_dir / "binder_interface_chemistry.json"
    if interface_chemistry_path.is_file():
        already_loaded = any(
            item.get("name") == "binder_interface_chemistry.json"
            and Path(str(item.get("path") or "")).resolve() == interface_chemistry_path.resolve()
            for item in results
        )
        payload = _read_json(interface_chemistry_path) if not already_loaded else None
        if isinstance(payload, dict):
            results.append({"name": "binder_interface_chemistry.json", "path": interface_chemistry_path,
                            "rows": payload.get("residue_annotations", []), "payload": payload})
    return results


def _database_inspector_data(user: dict, record: dict, structure: dict, path: Path,
                             components: list[dict], associated: list[tuple[dict, dict]]) -> dict:
    from mn_molecular_design_studio.complex_analysis_inspector import (
        author_residue_selector, inspector_evidence_rows, observation_context,
        parse_author_residue_selector, protein_sequence_tracks, resolve_author_residue_selector,
        structure_digest,
    )

    protein_chains = sorted({str(item.get("source_chain") or item.get("chain") or item.get("chain_id") or "")
                             for item in components
                             if item.get("type") == "protein"
                             and (item.get("source_chain") or item.get("chain") or item.get("chain_id"))})
    try:
        tracks = protein_sequence_tracks(path, protein_chains)
    except (ImportError, OSError, ValueError, TypeError):
        tracks = []
    current_digest = structure_digest(path)
    ligand_components = [item for item in components if str(item.get("type") or "").casefold() == "small_molecule"]
    ligand_selection = {}
    for observation, _evaluation in associated:
        context = observation_context(observation)
        candidate = context.get("candidate") if isinstance(context.get("candidate"), dict) else {}
        raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
        selection = raw.get("complex_analysis_selection") if isinstance(raw.get("complex_analysis_selection"), dict) else {}
        candidate_selection = selection.get("ligand_selection") if isinstance(selection.get("ligand_selection"), dict) else {}
        if candidate_selection:
            ligand_selection = candidate_selection
            break
    if not ligand_selection:
        context = structure.get("context") if isinstance(structure.get("context"), dict) else {}
        ligand_selection = context.get("ligand_selection") if isinstance(context.get("ligand_selection"), dict) else {}
    if not ligand_selection and ligand_components:
        ligand = ligand_components[0]
        ligand_selection = {
            "chain": str(ligand.get("source_chain") or ligand.get("chain") or ""),
            "residue_key": str(ligand.get("residue_key") or ""),
            "ccd": str(ligand.get("ccd") or ""),
            "smiles": str(ligand.get("smiles") or ""),
        }
    ligand_selector = ""
    ligand_key = str(ligand_selection.get("residue_key") or "")
    parsed_ligand_key = parse_author_residue_selector(ligand_key)
    if parsed_ligand_key:
        ligand_selector = author_residue_selector(
            str(ligand_selection.get("chain") or parsed_ligand_key["chain_id"] or "_"),
            parsed_ligand_key["number"], parsed_ligand_key["insertion_code"] or "",
        )
    elif ligand_selection.get("chain") and ligand_selection.get("residue_number") is not None:
        ligand_selector = author_residue_selector(
            str(ligand_selection["chain"]), int(ligand_selection["residue_number"]),
            str(ligand_selection.get("insertion_code") or ""),
        )

    evaluations = []
    evidence_rows = []
    for observation, descriptor in associated:
        evaluator = str(observation.get("namespace") or descriptor.get("evaluator") or "unknown").casefold()
        artifacts = _saved_evaluation_artifacts(user, observation, structure, current_digest)
        item = {"observation": observation, "evaluation": descriptor, "evaluator": evaluator,
                "artifacts": artifacts}
        context = observation_context(observation)
        candidate = context.get("candidate") if isinstance(context.get("candidate"), dict) else {}
        raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
        saved_tag_artifact = next((artifact.get("rows", [None])[0]
                                   for artifact in artifacts
                                   if artifact.get("name") == "tag_fusion_suitability.json"
                                   and artifact.get("rows")), None)
        # The access-checked run artifact is the authoritative evaluator output.
        # Candidate normalization may resolve short DSSP assignment codes such
        # as "H" or "-" as paths in older saved Evaluation metadata.
        tag_result = (saved_tag_artifact if isinstance(saved_tag_artifact, dict)
                      else raw.get("tag_fusion_suitability"))
        if evaluator in {"tag_fusion_suitability", "tag_fusion_context"} and isinstance(tag_result, dict):
            tag_source = tag_result.get("source") if isinstance(tag_result.get("source"), dict) else {}
            expected_model = str((structure.get("context") or {}).get("model_id")
                                 or (structure.get("context") or {}).get("model_index")
                                 or structure.get("model_index") or structure.get("id") or "")
            if (str(tag_source.get("molecule_id") or "") == str(record.get("id") or "")
                    and str(tag_source.get("structure_id") or "") == str(structure.get("id") or "")
                    and str(tag_source.get("structure_digest") or "") == current_digest
                    and str(tag_source.get("model_id") or "") == expected_model):
                item["tag_fusion_result"] = tag_result
                verified_terminal_selectors = {}
                for end, terminal in (tag_result.get("terminal_geometry") or {}).items():
                    identity = terminal.get("terminal_residue") if isinstance(terminal, dict) else None
                    if not isinstance(identity, dict):
                        continue
                    selector = resolve_author_residue_selector(
                        tracks, str(identity.get("chain_id") or ""), identity.get("author_number"),
                        str(identity.get("insertion_code") or ""),
                    )
                    verified_terminal_selectors[str(end)] = selector
                    residue_ref = (f"{identity.get('chain_id')}:{identity.get('author_number')}"
                                   f"{identity.get('insertion_code') or ''}")
                    evidence_rows.append({
                        "id": f"{evaluator}:terminus:{end}:{descriptor.get('id')}",
                        "evaluator": evaluator, "artifact": "tag_fusion_suitability.json",
                        "label": f"{end}-terminus geometry · {residue_ref}",
                        "residue_refs": [residue_ref], "selectors": [selector] if selector else [],
                        "raw": terminal,
                    })
                item["tag_fusion_terminal_selectors"] = verified_terminal_selectors
        environment_artifact = next((artifact for artifact in artifacts
                                     if artifact.get("name") == "reference_environment.json"
                                     and isinstance(artifact.get("payload"), dict)), None)
        environment_result = environment_artifact.get("payload") if environment_artifact else None
        if evaluator == "reference_environment" and isinstance(environment_result, dict):
            env_source = environment_result.get("source") if isinstance(environment_result.get("source"), dict) else {}
            env_reference = environment_result.get("reference") if isinstance(environment_result.get("reference"), dict) else {}
            expected_model = str((structure.get("context") or {}).get("model_id")
                                 or (structure.get("context") or {}).get("model_index")
                                 or structure.get("model_index") or structure.get("id") or "")
            reference_is_current = False
            try:
                from mn_molecular_design_studio.molecular_database import MolecularDatabase
                from mn_molecular_design_studio.evaluations import file_digest
                db = MolecularDatabase()
                saved_context = db.get_complex_reference_context(
                    user, source_molecule_id=str(record.get("id") or ""),
                    source_structure_id=str(structure.get("id") or ""),
                    source_chain=str(env_source.get("chain_id") or ""),
                )
                if (saved_context and saved_context.get("status") == "validated"
                        and str(saved_context.get("id") or "") == str(environment_result.get("reference_context_id") or "")
                        and str(saved_context.get("reference_molecule_id") or "") == str(env_reference.get("molecule_id") or "")
                        and str(saved_context.get("reference_structure_id") or "") == str(env_reference.get("structure_id") or "")
                        and str(saved_context.get("reference_structure_digest") or "") == str(env_reference.get("structure_digest") or "")
                        and str(saved_context.get("reference_model_id") or "") == str(env_reference.get("model_id") or "")
                        and str(saved_context.get("reference_component_chain") or "") == str(env_reference.get("chain_id") or "")):
                    current_reference = db.structure_file(user, str(env_reference.get("structure_id") or ""), materialize=True)
                    reference_is_current = file_digest(current_reference) == str(env_reference.get("structure_digest") or "")
            except (OSError, PermissionError, ValueError, TypeError, KeyError):
                reference_is_current = False
            if (reference_is_current
                    and str(env_source.get("molecule_id") or "") == str(record.get("id") or "")
                    and str(env_source.get("structure_id") or "") == str(structure.get("id") or "")
                    and str(env_source.get("structure_digest") or "") == current_digest
                    and str(env_source.get("model_id") or "") == expected_model):
                item["reference_environment_result"] = environment_result
                for row in environment_result.get("residue_evidence", []):
                    if not isinstance(row, dict):
                        continue
                    for role_key, role_label in (
                        ("source_component_geometry", "associated source component"),
                        ("other_source_complex_geometry", "other source-complex component"),
                    ):
                        role_geometry = row.get(role_key) if isinstance(row.get(role_key), dict) else {}
                        if not int(role_geometry.get("atom_pairs_within_proximity_cutoff") or 0):
                            continue
                        pair = role_geometry.get("closest_atom_pair")
                        src = pair.get("source") if isinstance(pair, dict) else None
                        if not isinstance(src, dict):
                            continue
                        selector = resolve_author_residue_selector(
                            tracks, str(src.get("chain_id") or ""), src.get("author_number"),
                            str(src.get("insertion_code") or ""),
                        )
                        if selector:
                            evidence_rows.append({
                                "id": f"reference_environment:{role_key}:{row.get('context_kind')}:{row.get('chain_id')}:{row.get('author_number')}:{row.get('insertion_code', '')}:{descriptor.get('id')}",
                                "evaluator": evaluator, "artifact": "reference_environment.json",
                                "label": f"Reference context {row.get('context_kind')} · {role_label} {src.get('chain_id')}:{src.get('author_residue')} · {role_geometry.get('closest_distance_angstrom')} Å",
                                "residue_refs": [f"{src.get('chain_id')}:{src.get('author_residue')}"],
                                "selectors": [selector], "raw": row,
                            })
        liability_artifact = next((artifact for artifact in artifacts
                                   if artifact.get("name") == "sequence_liability.json"
                                   and isinstance(artifact.get("payload"), dict)), None)
        liability_result = liability_artifact.get("payload") if liability_artifact else raw.get("sequence_liability")
        if evaluator == "sequence_liability" and isinstance(liability_result, dict):
            liability_source = liability_result.get("source") if isinstance(liability_result.get("source"), dict) else {}
            expected_model = str((structure.get("context") or {}).get("model_id")
                                 or (structure.get("context") or {}).get("model_index")
                                 or structure.get("model_index") or structure.get("id") or "")
            chains = [str(value) for value in liability_result.get("selected_component_chains") or []]
            if (str(liability_source.get("molecule_id") or "") == str(record.get("id") or "")
                    and str(liability_source.get("structure_id") or "") == str(structure.get("id") or "")
                    and str(liability_source.get("structure_digest") or "") == current_digest
                    and str(liability_source.get("model_id") or "") == expected_model
                    and chains and set(chains) <= set(protein_chains)):
                item["sequence_liability_result"] = liability_result
                for finding in liability_result.get("findings") or []:
                    if not isinstance(finding, dict):
                        continue
                    chain_id = str(finding.get("chain_id") or "")
                    number = finding.get("author_number")
                    if number in (None, "") or chain_id not in chains:
                        selector = ""
                    else:
                        selector = resolve_author_residue_selector(
                            tracks, chain_id, number, str(finding.get("insertion_code") or ""),
                        )
                    residue_ref = str(finding.get("residue_ref") or f"{chain_id}:sequence#{finding.get('sequence_position') or '?'}")
                    evidence_rows.append({
                        "id": f"sequence_liability:{finding.get('finding_id') or finding.get('kind')}:{descriptor.get('id')}",
                        "evaluator": evaluator, "artifact": "sequence_liability.json",
                        "label": f"{str(finding.get('severity') or 'info').upper()} · {finding.get('kind') or 'liability'} · {residue_ref}",
                        "residue_refs": [residue_ref], "selectors": [selector] if selector else [],
                        "raw": finding,
                    })
        interface_chemistry_artifact = next((artifact for artifact in artifacts
            if artifact.get("name") == "binder_interface_chemistry.json"
            and isinstance(artifact.get("payload"), dict)), None)
        interface_chemistry = (interface_chemistry_artifact.get("payload") if interface_chemistry_artifact
                               else raw.get("binder_interface_chemistry"))
        if evaluator == "binder_interface_chemistry_v1" and isinstance(interface_chemistry, dict):
            source = interface_chemistry.get("source") if isinstance(interface_chemistry.get("source"), dict) else {}
            roles = interface_chemistry.get("selected_chains") if isinstance(interface_chemistry.get("selected_chains"), dict) else {}
            binder = [str(value) for value in roles.get("binder") or []]
            target = [str(value) for value in roles.get("target") or []]
            expected_model = str((structure.get("context") or {}).get("model_id")
                                 or (structure.get("context") or {}).get("model_index")
                                 or structure.get("model_index") or structure.get("id") or "")
            if (str(source.get("molecule_id") or "") == str(record.get("id") or "")
                    and str(source.get("structure_id") or "") == str(structure.get("id") or "")
                    and str(source.get("structure_digest") or "") == current_digest
                    and str(source.get("model_id") or "") == expected_model
                    and binder and target and not (set(binder) & set(target)
                    or (set(binder) | set(target)) - set(protein_chains))):
                item["binder_interface_chemistry_result"] = interface_chemistry
                for residue in interface_chemistry.get("residue_annotations") or []:
                    if not isinstance(residue, dict):
                        continue
                    if not (residue.get("heavy_atom_contact_count") or residue.get("buried_unpaired_polar_atoms")
                            or residue.get("hydrogen_bond_evidence_count")):
                        continue
                    chain_id = str(residue.get("chain_id") or "")
                    selector = resolve_author_residue_selector(
                        tracks, chain_id, residue.get("author_number"), str(residue.get("insertion_code") or ""),
                    )
                    residue_ref = f"{chain_id}:{residue.get('residue_name')}{residue.get('author_residue')}"
                    evidence_rows.append({
                        "id": f"binder_interface_chemistry:{residue_ref}:{descriptor.get('id')}",
                        "evaluator": evaluator, "artifact": "binder_interface_chemistry.json",
                        "label": f"Interface · {residue_ref} · {residue.get('side')}",
                        "residue_refs": [residue_ref], "selectors": [selector] if selector else [],
                        "raw": residue,
                    })
        evaluations.append(item)
        for artifact in artifacts:
            if artifact["name"] in {"tag_fusion_suitability.json", "sequence_liability.json", "binder_interface_chemistry.json", "reference_environment.json", "report.json", "patch_metrics.json"}:
                # These are summary measurements shown in their evaluator panel,
                # not residue-level findings to attach to a Mol* selection.
                continue
            evidence_rows.extend(inspector_evidence_rows(
                evaluator, artifact["name"], artifact["rows"], tracks,
                ligand_selector=ligand_selector,
            ))
    confidence_keys = ("plddt", "pae", "ptm", "iptm", "ipae", "confidence")
    confidence = []
    for item in evaluations:
        observation = item["observation"]
        metrics = observation.get("metrics") if isinstance(observation.get("metrics"), dict) else {}
        values = {key: value for key, value in metrics.items() if any(stem in str(key).casefold() for stem in confidence_keys)}
        context = observation_context(observation)
        candidate = context.get("candidate") if isinstance(context.get("candidate"), dict) else {}
        raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
        outputs = raw.get("confidence_outputs")
        if values or outputs:
            confidence.append({**item, "confidence_metrics": values, "confidence_outputs": outputs or {}})
    return {
        "tracks": tracks,
        "evaluations": evaluations,
        "evidence": evidence_rows,
        "confidence": confidence,
        "ligand_selection": ligand_selection,
        "ligand_components": ligand_components,
        "structure_digest": current_digest,
        "protein_chains": protein_chains,
        "record": record,
        "structure": structure,
        "path": path,
    }


def _render_inspector_measurements(items: list[dict], *, evaluators: set[str] | None = None) -> None:
    import pandas as pd

    selected = [item for item in items if evaluators is None or item["evaluator"] in evaluators]
    if not selected:
        st.info("No saved evaluation measurements are available for this panel.")
        return
    for item in selected:
        observation, evaluation = item["observation"], item["evaluation"]
        st.markdown(f"**{item['evaluator']}** · `{evaluation.get('id')}`")
        st.caption(f"Run {evaluation.get('run_id') or observation.get('run_id') or 'unknown'} · saved measurement values retain the evaluator's native definitions and units.")
        metrics = observation.get("metrics") if isinstance(observation.get("metrics"), dict) else {}
        if metrics:
            st.dataframe(pd.DataFrame([{"Measurement": key, "Stored value": value}
                                       for key, value in metrics.items()]), hide_index=True, use_container_width=True)
        normalized = observation.get("normalized_output") if isinstance(observation.get("normalized_output"), dict) else {}
        normalized_metrics = normalized.get("metrics") if isinstance(normalized.get("metrics"), dict) else {}
        normalized_provenance = normalized.get("provenance") if isinstance(normalized.get("provenance"), dict) else {}
        if normalized_metrics:
            with st.expander("Existing Studio measurement definitions and units"):
                rows = []
                for key, value in normalized_metrics.items():
                    definition = normalized_provenance.get(key) if isinstance(normalized_provenance.get(key), dict) else {}
                    rows.append({
                        "Registered measure": key,
                        "Stored projection": value,
                        "Unit": definition.get("unit") or "not specified",
                        "Scope": definition.get("scope") or "not specified",
                        "Method": definition.get("method") or "not specified",
                        "Source tool": definition.get("source_tool") or item["evaluator"],
                    })
                st.dataframe(pd.DataFrame(rows), hide_index=True, use_container_width=True)
        for artifact in item["artifacts"]:
            rows = artifact["rows"]
            if rows and artifact["name"] not in {"report.json", "patch_metrics.json"}:
                with st.expander(f"Saved {artifact['name']} · {len(rows)} row(s)"):
                    st.dataframe(pd.DataFrame(rows), hide_index=True, use_container_width=True)
            elif rows:
                with st.expander(f"Saved {artifact['name']}"):
                    st.json(rows[0])
        with st.expander("Evaluation provenance", expanded=False):
            context = observation.get("context") if isinstance(observation.get("context"), dict) else {}
            st.json({
                "evaluation_id": evaluation.get("id"), "observation_id": observation.get("id"),
                "evaluator": item["evaluator"], "run_id": evaluation.get("run_id") or observation.get("run_id"),
                "source_structure_id": context.get("source_structure_id") or context.get("representation_id"),
                "source_structure_digest": context.get("source_structure_digest") or context.get("structure_digest"),
                "settings": context.get("settings") or {}, "evaluation_context": context,
            })


def _render_inspector_findings(data: dict, viewer_key: str) -> None:
    evidence = data["evidence"]
    if not evidence:
        st.info("No saved residue-linked artifact evidence is available for this structure. Existing scalar evaluations remain available in their measurement panels.")
        return
    st.caption("Choose a saved evaluator observation and highlight its exact residue references in the structure. Raw evaluator fields, units, and Evaluation provenance are retained below.")
    evidence_ids = [item["id"] for item in evidence]
    labels = {item["id"]: f"{item['evaluator']} · {item['label']}" for item in evidence}
    select_key = f"{viewer_key}-inspector-evidence-select"
    if st.session_state.get(select_key) not in evidence_ids:
        st.session_state[select_key] = evidence_ids[0]
    selected_id = st.selectbox("Saved finding / interaction", evidence_ids,
                               format_func=lambda value: labels[value], key=select_key)
    selected = next(item for item in evidence if item["id"] == selected_id)
    state_key = f"{viewer_key}-inspector-highlights"
    cols = st.columns([1, 1])
    cols[0].button("Highlight evidence in Mol*", key=f"{viewer_key}-inspector-evidence-highlight",
                   disabled=not selected["selectors"],
                   on_click=_select_inspector_evidence, args=(state_key, selected["selectors"]))
    cols[1].button("Clear Inspector highlight", key=f"{viewer_key}-inspector-evidence-clear",
                   on_click=_select_inspector_evidence, args=(state_key, []))
    if selected["selectors"]:
        st.caption("Mol* selection: " + ", ".join(selected["selectors"]))
    else:
        st.info("This output does not contain a uniquely resolvable author residue reference in the selected representation, so no highlight is inferred.")
    st.json({"evaluator": selected["evaluator"], "artifact": selected["artifact"],
             "residue_references": selected["residue_refs"], "stored_evaluator_fields": selected["raw"]})
    descriptor = next((item["evaluation"] for item in data["evaluations"]
                       if item["evaluator"] == selected["evaluator"]), {})
    st.caption(f"Evidence Evaluation ID: {descriptor.get('id') or 'unavailable'}")


def _tag_fusion_terminal_metric_rows(terminal: dict) -> list[dict]:
    """Build Arrow-safe rows for numeric and categorical terminal evidence."""
    definitions = (
        ("Distance to selected interface (Å)", "distance_to_selected_interface_angstrom"),
        ("Terminal exit direction toward partner (degrees)", "aim_toward_partner_degrees"),
        ("Interface-side cosine", "interface_side_cosine"),
        ("Terminal residue SASA (Å²)", "terminal_residue_sasa_angstrom2"),
        ("DSSP relative accessibility", "relative_sasa_sander"),
        ("Secondary structure", "secondary_structure"),
        ("Terminal window mean pLDDT", "terminal_window_mean_plddt"),
        ("Terminal window mean B-factor", "terminal_window_mean_b_factor"),
        ("Terminal contact at cutoff", "contact_at_terminal_residue"),
        ("Loop residues in terminal window", "terminal_window_loop_residue_count"),
        ("Flexibility evidence", "flexibility_evidence"),
    )
    rows = []
    for label, field in definitions:
        value = terminal.get(field)
        if value is None:
            continue
        if field == "secondary_structure":
            value = _display_dssp_assignment(value)
        rows.append({"Evidence": label, "Stored value": str(value)})
    return rows


def _display_dssp_assignment(value: object) -> str:
    """Keep DSSP codes readable when older normalized metadata stored them as paths."""
    assignment = str(value).strip()
    dssp_codes = {"H", "G", "I", "E", "B", "T", "S", "-"}
    if assignment in dssp_codes:
        return assignment
    normalized = assignment.replace("\\", "/")
    if normalized.startswith("/") or normalized.startswith("runs:///"):
        leaf = normalized.rsplit("/", 1)[-1]
        if leaf in dssp_codes:
            return leaf
        return "unavailable (invalid stored assignment)"
    return assignment


def _render_inspector_tag_fusion(items: list[dict], viewer_key: str) -> None:
    import pandas as pd

    results = [item for item in items if isinstance(item.get("tag_fusion_result"), dict)]
    if not results:
        st.info("No matching saved Tag and Fusion Evaluation is available for this exact Database representation.")
        return
    options = [str(item["evaluation"].get("id") or index) for index, item in enumerate(results)]
    labels = {option: f"{item['evaluator']} · {option}" for option, item in zip(options, results)}
    key = f"{viewer_key}-tag-fusion-evaluation"
    if st.session_state.get(key) not in options:
        st.session_state[key] = options[0]
    selected_id = st.selectbox("Tag and Fusion Evaluation", options,
                               format_func=lambda value: labels[value], key=key)
    item = results[options.index(selected_id)]
    result = item["tag_fusion_result"]
    st.caption(f"{result.get('method')} · version {result.get('method_version')} · {result.get('assessment_scope')}")
    st.markdown(f"**Assessment status:** `{result.get('assessment_status', 'unknown')}`")
    if result.get("not_applicable_reason"):
        st.info(result["not_applicable_reason"])
    geometry = result.get("terminal_geometry") if isinstance(result.get("terminal_geometry"), dict) else {}
    terminal_columns = st.columns(2)
    for column, end in zip(terminal_columns, ("N", "C")):
        with column:
            terminal = geometry.get(end) if isinstance(geometry.get(end), dict) else {}
            residue = terminal.get("terminal_residue") if isinstance(terminal.get("terminal_residue"), dict) else {}
            selector = str((item.get("tag_fusion_terminal_selectors") or {}).get(end) or "")
            st.markdown(f"#### {end}-terminus")
            st.caption(f"{residue.get('residue_name') or 'Residue'} {residue.get('author_number', 'unavailable')}{residue.get('insertion_code') or ''} · chain {residue.get('chain_id') or 'unavailable'} · {terminal.get('assessment_status', 'not assessed')}")
            if terminal:
                metrics = _tag_fusion_terminal_metric_rows(terminal)
                if metrics:
                    st.dataframe(pd.DataFrame(metrics), hide_index=True, use_container_width=True)
                if selector:
                    st.button(f"Highlight {end}-terminus in Mol*", key=f"{viewer_key}-tag-fusion-highlight-{end}",
                              on_click=_select_inspector_evidence,
                              args=(f"{viewer_key}-inspector-highlights", [selector]))
                clearances = terminal.get("tag_clearance") or {}
                if clearances:
                    rows = [{"Model": value.get("label"), "Measured geometric status": value.get("status"),
                             "Suitability advisory": (value.get("suitability_advisory") or {}).get("status"),
                             "Advisory basis": (value.get("suitability_advisory") or {}).get("reason"),
                             "Minimum tested linker (residues)": value.get("minimum_linker_residues_tested"),
                             "Path (Å)": value.get("shortest_grid_path_angstrom")}
                            for value in clearances.values()]
                    st.dataframe(pd.DataFrame(rows), hide_index=True, use_container_width=True)
                    st.caption("Measured envelope clearance is geometric evidence. The separate suitability advisory accounts for terminal interface participation and does not predict tag behavior experimentally.")
            else:
                st.caption("This component format has no free, assessable N/C terminus under the selected method.")
    comparison = result.get("tag_comparison") if isinstance(result.get("tag_comparison"), dict) else {}
    if comparison:
        st.markdown("#### Advisory comparison")
        rows = [{"Tag / fusion model": value.get("label", key),
                 "Preferred end advisory": value.get("preferred_end_advisory", value.get("recommendation", "not assessed")),
                 "Basis": value.get("advice_basis", value.get("basis", ""))}
                for key, value in comparison.items() if isinstance(value, dict)]
        if rows:
            st.dataframe(pd.DataFrame(rows), hide_index=True, use_container_width=True)
    st.caption("These are coordinate-based advisories. Flexible polymer conformations, tag–receptor interactions, expression, and experimental SPR/BLI behavior are not simulated.")
    with st.expander("Evaluation provenance, definitions and limitations"):
        st.json({"evaluation_id": item["evaluation"].get("id"), "run_id": item["evaluation"].get("run_id"),
                 "source": result.get("source"), "parameters": result.get("parameters"),
                 "method_parameters": result.get("method_parameters"),
                 "reference_context": result.get("reference_context"),
                 "checks_not_performed": result.get("checks_not_performed"),
                 "interpretation_limitations": result.get("interpretation_limitations"),
                 "adaptation_attribution": result.get("adaptation_attribution")})


def _render_inspector_liabilities(items: list[dict], viewer_key: str, evidence: list[dict]) -> None:
    import pandas as pd

    results = [item for item in items if isinstance(item.get("sequence_liability_result"), dict)]
    if not results:
        st.info("No matching saved Sequence and Developability Liability Evaluation is available for this exact Database representation.")
        return
    options = [str(item["evaluation"].get("id") or index) for index, item in enumerate(results)]
    labels = {option: f"{item['evaluation'].get('run_id') or item['observation'].get('run_id') or 'saved evaluation'} · {option}"
              for option, item in zip(options, results)}
    key = f"{viewer_key}-inspector-liability-evaluation"
    if st.session_state.get(key) not in options:
        st.session_state[key] = options[0]
    selected_id = st.selectbox("Saved liability evaluation", options, format_func=lambda value: labels[value], key=key)
    item = results[options.index(selected_id)]
    result = item["sequence_liability_result"]
    st.caption("Sequence alerts are susceptibility motifs, not observed modifications or experimental stability measurements. Structural findings are tied to the saved model and their own evaluator evidence.")
    st.markdown(f"**Assessment status:** {result.get('assessment_status', 'unknown')} · **Method:** {result.get('method', 'unavailable')} · **Version:** {result.get('method_version', 'unavailable')}")
    properties = []
    for chain_id, record in (result.get("sequence_properties") or {}).items():
        metrics = record.get("metrics") if isinstance(record.get("metrics"), dict) else {}
        properties.append({"Chain": chain_id, "Sequence calculation": record.get("status"),
                           "pI": metrics.get("qc_sequence_pi"), "Net charge at pH 7.4": metrics.get("qc_sequence_charge"),
                           "GRAVY": metrics.get("qc_sequence_gravy"), "Length": metrics.get("qc_sequence_length"),
                           "Method": record.get("method") or record.get("evaluator") or "unavailable"})
    if properties:
        st.markdown("#### Sequence properties")
        st.dataframe(pd.DataFrame(properties), hide_index=True, use_container_width=True)
    coverage = result.get("sequence_coverage") or {}
    if coverage:
        st.markdown("#### Sequence and coordinate coverage")
        rows = []
        for chain_id, value in coverage.items():
            rows.append({"Chain": chain_id, "Sequence source": value.get("sequence_source"),
                         "Coverage": value.get("status"), "Sequence residues": value.get("sequence_length"),
                         "Coordinate residues": value.get("coordinate_residue_count"),
                         "Unresolved/unmapped": value.get("unresolved_or_unmapped_residue_count"),
                         "Termini": (value.get("terminal_assessment") or {}).get("status"),
                         "Explanation": (value.get("terminal_assessment") or {}).get("reason") or value.get("reason")})
        st.dataframe(pd.DataFrame(rows), hide_index=True, use_container_width=True)
    findings = result.get("findings") or []
    if findings:
        st.markdown("#### Residue-linked findings")
        st.dataframe(pd.DataFrame([{
            "Severity": str(row.get("severity") or "info").upper(), "Evidence basis": row.get("basis"),
            "Finding": row.get("kind"), "Residue": row.get("residue_ref"), "Motif": row.get("motif"),
            "Interface status": row.get("interface_status"), "Observation": row.get("message"),
        } for row in findings]), hide_index=True, use_container_width=True)
    else:
        st.info("No supported liability findings were reported for this selected sequence/model. This does not establish absence of experimental liabilities.")
    missing = result.get("missing_data") or []
    if missing:
        st.markdown("#### Unavailable or unsupported inputs")
        st.dataframe(pd.DataFrame([{"Input": row.get("data"), "Status": row.get("status"), "Reason": row.get("reason")}
                                   for row in missing]), hide_index=True, use_container_width=True)
    with st.expander("Methods and provenance", expanded=False):
        st.json({"evaluation_id": item["evaluation"].get("id"), "run_id": item["evaluation"].get("run_id"),
                 "source": result.get("source"), "selected_component_chains": result.get("selected_component_chains"),
                 "selected_partner_chains": result.get("selected_partner_chains"), "parameters": result.get("parameters"),
                 "method_provenance": result.get("method_provenance"),
                 "interpretation_limitations": result.get("interpretation_limitations"),
                 "attribution": result.get("attribution")})
    liability_evidence = [row for row in evidence if row.get("evaluator") == "sequence_liability"
                          and str(item["evaluation"].get("id") or "") in str(row.get("id") or "")]
    if liability_evidence:
        evidence_ids = [row["id"] for row in liability_evidence]
        names = {row["id"]: row.get("label") or row["id"] for row in liability_evidence}
        select_key = f"{viewer_key}-liability-finding-select"
        if st.session_state.get(select_key) not in evidence_ids:
            st.session_state[select_key] = evidence_ids[0]
        selected_evidence_id = st.selectbox("Finding to highlight", evidence_ids,
                                            format_func=lambda value: names[value], key=select_key)
        selected_evidence = next(row for row in liability_evidence if row["id"] == selected_evidence_id)
        st.button("Highlight this liability in Mol*", key=f"{viewer_key}-liability-highlight",
                  disabled=not selected_evidence.get("selectors"),
                  on_click=_select_inspector_evidence,
                  args=(f"{viewer_key}-inspector-highlights", selected_evidence.get("selectors") or []))
        if not selected_evidence.get("selectors"):
            st.info("This finding refers to an unresolved sequence position, so no structural highlight is inferred.")


def _render_inspector_tracks(data: dict, viewer_key: str) -> None:
    tracks = data["tracks"]
    if not tracks:
        return
    with st.expander("Interactive author-numbered sequence tracks", expanded=True):
        st.caption("Each residue is shown in sequence order with the saved author residue number and insertion code. Select a tile to highlight that exact residue in Mol*.")
        chain_ids = [str(item["chain_id"]) for item in tracks]
        chain_key = f"{viewer_key}-inspector-track-chain"
        if st.session_state.get(chain_key) not in chain_ids:
            st.session_state[chain_key] = chain_ids[0]
        chain_id = st.selectbox("Protein component track", chain_ids, key=chain_key)
        track = next(item for item in tracks if item["chain_id"] == chain_id)
        page_size = 24
        page_count = max(1, (len(track["residues"]) + page_size - 1) // page_size)
        page_key = f"{viewer_key}-inspector-track-page-{chain_id}"
        if int(st.session_state.get(page_key, 0) or 0) >= page_count:
            st.session_state[page_key] = 0
        page = st.selectbox("Sequence segment", list(range(page_count)),
                            format_func=lambda index: f"Residues {index * page_size + 1}–{min((index + 1) * page_size, len(track['residues']))} of {len(track['residues'])}",
                            key=page_key)
        start = page * page_size
        residues = track["residues"][start:start + page_size]
        st.caption("Sequence: " + "".join(row["amino_acid"] for row in residues))
        for offset in range(0, len(residues), 12):
            batch = residues[offset:offset + 12]
            columns = st.columns(12)
            for column, residue in zip(columns, batch):
                column.button(
                    f"{residue['amino_acid']}\n{residue['author_residue']}",
                    key=f"{viewer_key}-inspector-track-{chain_id}-{residue['sequence_index']}",
                    help=f"Chain {chain_id} · author residue {residue['author_residue']} · sequence index {residue['sequence_index'] + 1}",
                    on_click=_select_inspector_evidence,
                    args=(f"{viewer_key}-inspector-highlights", [residue["selector"]]),
                    disabled=not residue["selector"],
                )


def _render_inspector_interface_chemistry(items: list[dict], viewer_key: str) -> None:
    """Show versioned interface chemistry values and exact author-numbered residue evidence."""
    import pandas as pd

    results = [(item, item.get("binder_interface_chemistry_result")) for item in items
               if isinstance(item.get("binder_interface_chemistry_result"), dict)]
    if not results:
        st.info("No exact-structure BinderScope-style interface chemistry Evaluation is available.")
        return
    for item, result in results:
        st.markdown(f"**{result.get('method')}** · Evaluation `{item['evaluation'].get('id')}`")
        st.caption(f"Method {result.get('method_version')} · source molecule {result.get('source', {}).get('molecule_id')} · structure {result.get('source', {}).get('structure_id')} · model {result.get('source', {}).get('model_id')} · SHA-256 {str(result.get('source', {}).get('structure_digest') or '')[:16]}")
        metrics = result.get("metrics") if isinstance(result.get("metrics"), dict) else {}
        st.dataframe([{
            "Two-side interface dSASA (Å²)": metrics.get("interface_total_delta_sasa_angstrom2"),
            "Half ΔSASA BSA (Å²)": metrics.get("interface_buried_surface_area_half_delta_sasa_angstrom2"),
            "C/S buried area (Å²)": metrics.get("hydrophobic_cs_buried_area_angstrom2"),
            "C/S share of ΔSASA": metrics.get("hydrophobic_cs_buried_fraction"),
            "Buried unpaired polar atoms": metrics.get("buried_unpaired_polar_atom_count"),
            "Interfacial H-bond geometry": metrics.get("interfacial_hydrogen_bond_count"),
            "Residue-pair contacts": metrics.get("contacting_residue_pair_count"),
        }], hide_index=True, use_container_width=True)
        residue_rows = [row for row in result.get("residue_annotations") or [] if isinstance(row, dict)]
        if residue_rows:
            st.markdown("##### Per-residue burial and contact evidence")
            st.dataframe([{
                "Side": row.get("side"), "Chain": row.get("chain_id"), "Residue": f"{row.get('residue_name')}{row.get('author_residue')}",
                "Sequence AA": row.get("amino_acid"), "ΔSASA (Å²)": row.get("buried_sasa_angstrom2"),
                "Complex SASA (Å²)": row.get("complex_sasa_angstrom2"), "Heavy-atom contacts": row.get("heavy_atom_contact_count"),
                "Partner residues": ", ".join(row.get("partner_residues") or []),
                "Buried unpaired polar atoms": len(row.get("buried_unpaired_polar_atoms") or []),
                "H-bond evidence": row.get("hydrogen_bond_evidence_count"),
            } for row in residue_rows], hide_index=True, use_container_width=True)
            highlight_rows = [row for row in residue_rows if row.get("selector") and
                              (row.get("heavy_atom_contact_count") or row.get("buried_unpaired_polar_atoms")
                               or row.get("hydrogen_bond_evidence_count"))]
            if highlight_rows:
                key = f"{viewer_key}-interface-chemistry-residue"
                labels = {f"{row['chain_id']}:{row['residue_name']}{row['author_residue']} · {row['side']}": row["selector"]
                          for row in highlight_rows}
                selected = st.selectbox("Interface residue to highlight in Mol*", list(labels), key=key)
                st.button("Highlight interface residue", key=f"{viewer_key}-interface-chemistry-highlight",
                          on_click=_select_inspector_evidence,
                          args=(f"{viewer_key}-inspector-highlights", [labels[selected]]))
        hbond_rows = [row for row in result.get("hydrogen_bonds") or [] if isinstance(row, dict)]
        if hbond_rows:
            st.markdown("##### Interfacial hydrogen-bond geometry")
            st.dataframe([{
                "Donor": f"{row.get('donor', {}).get('chain_id')}:{row.get('donor', {}).get('residue_name')}{row.get('donor', {}).get('author_number')} {row.get('donor', {}).get('atom')}",
                "Acceptor": f"{row.get('acceptor', {}).get('chain_id')}:{row.get('acceptor', {}).get('residue_name')}{row.get('acceptor', {}).get('author_number')} {row.get('acceptor', {}).get('atom')}",
                "D–A (Å)": row.get("donor_acceptor_distance_angstrom"),
                "H–A (Å)": row.get("hydrogen_acceptor_distance_angstrom"),
                "D–H–A (°)": row.get("donor_hydrogen_acceptor_angle_degrees"),
                "Evidence basis": row.get("evidence_basis"),
            } for row in hbond_rows], hide_index=True, use_container_width=True)
        with st.expander("Method parameters, evidence basis and limitations"):
            st.json({"parameters": result.get("parameters"), "hydrogen_bond_basis_counts": metrics.get("interfacial_hydrogen_bond_evidence_basis_counts"),
                     "assumptions": result.get("methodological_assumptions"), "attribution": result.get("attribution")})


def _render_database_inspector(user: dict, data: dict, viewer_key: str, *, record: dict,
                               structure: dict, components: list[dict], model_id: str,
                               reference_view_context: dict | None = None) -> None:
    items = data["evaluations"]
    evaluators = {item["evaluator"] for item in items}
    interface_evaluators = {"biotite_interface", "rosetta_metrics", "ipsae", "plip", "pandamap", "binder_interface_chemistry_v1"}
    surface_evaluators = {"pep_patch", "posebusters"}
    tag_fusion_evaluators = {"tag_fusion_suitability", "tag_fusion_context"}
    liability_evaluators = {"sequence_liability"}
    ligand_present = bool(data["ligand_components"])
    panel_specs = []
    if evaluators & interface_evaluators:
        panel_specs.append(("Interface", "interface"))
    if any(isinstance(item.get("binder_interface_chemistry_result"), dict) for item in items):
        panel_specs.append(("Interface chemistry", "interface_chemistry"))
    if data["confidence"]:
        panel_specs.append(("Confidence", "confidence"))
    if evaluators & surface_evaluators or ligand_present:
        panel_specs.append(("Surface / Chemistry", "chemistry"))
    if data["evidence"]:
        panel_specs.append(("Findings", "findings"))
    if data["protein_chains"]:
        panel_specs.append(("Reference Context", "reference"))
    if any(isinstance(item.get("reference_environment_result"), dict) for item in items):
        panel_specs.append(("Spatial Context", "spatial_context"))
    if evaluators & tag_fusion_evaluators:
        panel_specs.append(("Tag & Fusion", "tag_fusion"))
    if evaluators & liability_evaluators:
        panel_specs.append(("Sequence Liabilities", "sequence_liabilities"))
    st.markdown("### Scientific Inspector")
    if not panel_specs:
        st.info("This saved structure has no applicable saved Evaluation data yet. Its structure and existing analysis launch controls remain available.")
        return
    tabs = st.tabs([title for title, _kind in panel_specs])
    for tab, (_title, kind) in zip(tabs, panel_specs):
        with tab:
            if kind == "interface":
                _render_inspector_measurements(items, evaluators=interface_evaluators)
            elif kind == "interface_chemistry":
                _render_inspector_interface_chemistry(items, viewer_key)
            elif kind == "confidence":
                for item in data["confidence"]:
                    st.markdown(f"**{item['evaluator']}** · `{item['evaluation'].get('id')}`")
                    if item["confidence_metrics"]:
                        st.dataframe([{"Measurement": key, "Stored value": value}
                                      for key, value in item["confidence_metrics"].items()],
                                     hide_index=True, use_container_width=True)
                    if item["confidence_outputs"]:
                        st.json(item["confidence_outputs"])
                    normalized = item["observation"].get("normalized_output")
                    provenance = normalized.get("provenance") if isinstance(normalized, dict) else None
                    if isinstance(provenance, dict) and provenance:
                        with st.expander("Confidence measurement definitions and units"):
                            keys = [key for key in (normalized.get("metrics") or {}) if key in provenance]
                            st.dataframe([{
                                "Registered measure": key,
                                "Unit": provenance[key].get("unit") or "not specified",
                                "Scope": provenance[key].get("scope") or "not specified",
                                "Method": provenance[key].get("method") or "not specified",
                                "Source tool": provenance[key].get("source_tool") or item["evaluator"],
                            } for key in keys], hide_index=True, use_container_width=True)
                    context = item["observation"].get("context")
                    if isinstance(context, dict):
                        st.caption(f"Confidence provenance · run {item['evaluation'].get('run_id') or item['observation'].get('run_id')} · values shown as stored by {item['evaluator']}.")
            elif kind == "chemistry":
                ligand = data["ligand_selection"]
                smiles = str(ligand.get("smiles") or "")
                if ligand_present:
                    st.markdown("#### Selected ligand identity")
                    if smiles:
                        st.caption("Saved isomeric SMILES")
                        st.code(smiles, language=None)
                    st.caption(f"Ligand residue {ligand.get('residue_key') or 'unavailable'} · chain {ligand.get('chain') or 'unavailable'} · CCD {ligand.get('ccd') or 'unavailable'}")
                _render_inspector_measurements(items, evaluators=surface_evaluators)
                if "posebusters" in evaluators:
                    st.info("PoseBusters results are structural and chemical pose-validation outputs, separate from ligand–receptor interaction analyses.")
            elif kind == "findings":
                _render_inspector_findings(data, viewer_key)
            elif kind == "reference":
                key = hashlib.sha256(f"database:{record['id']}:{structure['id']}".encode()).hexdigest()[:14]
                _render_optional_reference_context(
                    user, record, structure, components, model_id, key,
                    view_context=reference_view_context,
                )
            elif kind == "tag_fusion":
                _render_inspector_tag_fusion(items, viewer_key)
            elif kind == "sequence_liabilities":
                _render_inspector_liabilities(items, viewer_key, data["evidence"])
            elif kind == "spatial_context":
                for item in items:
                    result = item.get("reference_environment_result")
                    if not isinstance(result, dict):
                        continue
                    st.markdown(f"**Reference environment** · Evaluation `{item['evaluation'].get('id')}`")
                    placement = result.get("placement") or {}
                    st.caption(
                        f"Placement fit: {placement.get('aligned_pair_count')} mapped Cα pairs · "
                        f"source coverage {placement.get('source_coverage')} · "
                        f"reference coverage {placement.get('reference_coverage')} · "
                        f"RMSD {placement.get('rmsd_angstrom')} Å · "
                        f"mapping {placement.get('mapping_method') or 'unavailable'}"
                    )
                    st.caption("Geometry is reported against both the associated source component and the other chains in the selected complex. Same-component contacts can reflect structural continuity; they are not partner-interface measurements.")
                    st.dataframe([result.get("measurements", {})], hide_index=True, use_container_width=True)
                    evidence = result.get("residue_evidence") or []
                    if evidence:
                        st.dataframe([{
                            "Context": row.get("context_kind"), "Residue": f"{row.get('chain_id')}:{row.get('residue_name')}{row.get('author_number')}{row.get('insertion_code') or ''}",
                            "Closest distance (Å)": row.get("closest_distance_angstrom"),
                            "Atom pairs ≤ 8 Å": row.get("atom_pairs_within_proximity_cutoff"),
                            "Potential VDW overlaps": row.get("atom_pairs_with_vdw_overlap_above_epsilon"),
                            "Associated component nearest (Å)": (row.get("source_component_geometry") or {}).get("closest_distance_angstrom"),
                            "Other complex chains nearest (Å)": (row.get("other_source_complex_geometry") or {}).get("closest_distance_angstrom"),
                            "Other chain potential overlaps": (row.get("other_source_complex_geometry") or {}).get("atom_pairs_with_vdw_overlap_above_epsilon"),
                            "Closest source role": ((row.get("closest_atom_pair") or {}).get("source_role")),
                            "Closest source atom": ((row.get("closest_atom_pair") or {}).get("source") or {}).get("atom_name"),
                            "Closest context atom": ((row.get("closest_atom_pair") or {}).get("context") or {}).get("atom_name"),
                        } for row in evidence], hide_index=True, use_container_width=True)
                    with st.expander("Method, exact identities and limitations"):
                        st.json({"method": result.get("method"), "method_version": result.get("method_version"),
                                 "source": result.get("source"), "reference": result.get("reference"),
                                 "placement": result.get("placement"), "parameters": result.get("parameters"),
                                 "reference_coordinate_origin": result.get("reference_coordinate_origin"),
                                 "sequence_only_potential_n_linked_sequons": result.get("sequence_only_potential_n_linked_sequons"),
                                 "limitations": result.get("limitations")})


def _render_biotite_results(user: dict, evaluation: dict, *, source_candidate_id: str, key: str,
                            structure_path: Path | None = None, structure_context: dict | None = None,
                            structure_components: list[dict] | None = None, viewer_key: str = "") -> None:
    run_dir = Path(evaluation["run_dir"])
    from mn_molecular_design_studio.quality_ui import render_contact_maps, render_progress

    status = str(evaluation.get("status") or "unknown").casefold()
    display_status = "Completed" if status in {"completed", "succeeded"} else status.title()
    st.caption(f"Biotite Interface · {display_status} · {evaluation.get('job_code') or evaluation.get('run_id')}")
    metadata_path = run_dir / "metadata.json"

    @st.fragment(run_every=5 if status not in {"completed", "succeeded", "failed", "cancelled", "canceled"} else None)
    def refresh_status():
        latest = str(_read_json(metadata_path).get("status") or status).casefold()
        st.metric("Analysis job status", latest.title())
        if status not in {"completed", "succeeded", "failed", "cancelled", "canceled"} and latest in {"completed", "succeeded", "failed", "cancelled", "canceled"}:
            st.rerun(scope="app")
    refresh_status()
    render_progress(run_dir, status, key + "-progress")
    if status in {"failed", "cancelled", "canceled"}:
        if evaluation.get("error"):
            st.error(evaluation["error"])
        else:
            st.error("This analysis job did not complete. Review its existing Results logs before retrying.")
    metrics_path = run_dir / "artifacts" / "analysis" / "qc_metrics.csv"
    if status in {"completed", "succeeded"} and metrics_path.is_file():
        import pandas as pd
        try:
            frame = pd.read_csv(metrics_path)
        except (OSError, ValueError, pd.errors.ParserError):
            frame = pd.DataFrame()
        if not frame.empty:
            st.markdown("##### Measurements saved by this evaluation")
            st.dataframe(frame, hide_index=True, use_container_width=True)
        render_contact_maps(run_dir, key + "-contact-map")
        binder_residues, target_residues, all_residues = _biotite_contact_residues(run_dir, source_candidate_id)
        if all_residues and structure_path is not None and viewer_key:
            residue_sets = {"No highlight": [], "All interface residues": all_residues,
                            "Binder side": binder_residues, "Target side": target_residues}
            st.session_state[f"{viewer_key}-biotite-residue-sets"] = residue_sets
            mode = st.selectbox("Highlight in structure viewer",
                list(residue_sets), key=f"{viewer_key}-biotite-highlight-mode")
            st.caption(f"{len(residue_sets[mode])} residues selected from the saved contact map; the structure viewer above will highlight them.")
        evaluation_id = _analysis_evaluation_id(run_dir, user, source_candidate_id)
        if evaluation_id:
            st.success(f"Saved Evaluation · {evaluation_id}")
        else:
            st.caption("The workflow Results are complete. The Database evaluation projection will appear when the run is indexed.")
    with st.expander("Analysis provenance", expanded=False):
        manifest = evaluation.get("manifest") or {}
        item = evaluation.get("selected_item") or {}
        metadata = evaluation.get("metadata") or {}
        database_source = metadata.get("database_source") if isinstance(metadata.get("database_source"), dict) else {}
        st.json({
            "evaluation_run": f"{evaluation.get('task_group')}/{evaluation.get('run_id')}",
            "workflow": "analysis.biotite_interface",
            "source_kind": metadata.get("upstream_source_kind") or "results",
            "source_job": (f"{metadata.get('upstream_task_group')}/{metadata.get('upstream_run_id')}"
                           if metadata.get("upstream_source_kind") != "database" else ""),
            "database_source": {key: database_source.get(key) for key in
                                ("molecule_id", "structure_id", "group_id", "evaluation_id", "model_id", "structure_digest")
                                if key in database_source},
            "source_candidate_id": item.get("candidate_id"),
            "selected_structure_sha256": item.get("structure_digest"),
            "binder_chains": item.get("binder_chains") or [],
            "target_chains": item.get("target_chains") or [],
            "method": manifest.get("tool"),
            "parameters": manifest.get("settings") or {},
            "source_candidate_structure": (manifest.get("selection") or {}).get("structure_path") or "source candidate structure",
            "source_job_code": evaluation.get("metadata", {}).get("upstream_job_code"),
        })


def _render_contextual_results(user: dict, evaluation: dict, *, workflow_id: str,
                               source_candidate_id: str, key: str) -> None:
    """Render existing metric tables and saved artifacts for a contextual workflow."""
    run_dir = Path(evaluation["run_dir"])
    status = str(evaluation.get("status") or "unknown").casefold()
    title = {"analysis.ipsae": "ipSAE & Structure Analysis",
             "analysis.rosetta": "Rosetta Interface Analysis",
             "analysis.pep_patch": "PEP-Patch", "analysis.plip": "PLIP interaction analysis",
             "analysis.pandamap": "PandaMap interaction analysis",
             "analysis.posebusters": "PoseBusters pose validation",
             "analysis.reference_placement": "Reference-frame placement",
             "analysis.reference_environment": "Reference-frame steric and environmental context",
             "analysis.tag_fusion_suitability": "Tag and Fusion Suitability",
             "analysis.tag_fusion_context": "Reference-context Tag and Fusion Suitability",
             "analysis.sequence_liability": "Sequence and Developability Liabilities",
             "analysis.binder_interface_chemistry": "BinderScope-style Interface Surface & Chemistry"}.get(workflow_id, workflow_id)
    ligand_workflow = workflow_id in {"analysis.plip", "analysis.pandamap", "analysis.posebusters"}
    if workflow_id in {"analysis.plip", "analysis.pandamap"}:
        st.caption("Interaction analysis · reports contacts and interaction classes for the selected ligand pose and receptor.")
    elif workflow_id == "analysis.posebusters":
        st.caption("Structural/chemical pose validation · reports geometry and chemistry checks; it is not an interaction analysis.")
    st.caption(f"{title} · {status.title()} · {evaluation.get('job_code') or evaluation.get('run_id')}")

    @st.fragment(run_every=5 if status not in {"completed", "succeeded", "failed", "cancelled", "canceled"} else None)
    def refresh_status():
        latest = str(_read_json(run_dir / "metadata.json").get("status") or status).casefold()
        st.metric("Analysis job status", latest.title())
        if status not in {"completed", "succeeded", "failed", "cancelled", "canceled"} and latest in {"completed", "succeeded", "failed", "cancelled", "canceled"}:
            st.rerun(scope="app")

    refresh_status()
    if workflow_id == "analysis.pep_patch":
        from mn_molecular_design_studio.quality_ui import render_progress
        render_progress(run_dir, status, key + "-progress")
    if status in {"failed", "cancelled", "canceled"}:
        st.error(evaluation.get("error") or "The analysis job did not complete. Review the existing Results logs before retrying.")
    else:
        import pandas as pd
        table_name = {"analysis.ipsae": "ipsae_metrics.csv",
                      "analysis.rosetta": "rosetta_metrics.csv",
                      "analysis.pep_patch": "qc_metrics.csv",
                      "analysis.plip": "qc_metrics.csv", "analysis.pandamap": "qc_metrics.csv",
                      "analysis.posebusters": "qc_metrics.csv",
                      "analysis.reference_placement": "reference_placement_metrics.csv",
                      "analysis.reference_environment": "reference_environment_metrics.csv",
                      "analysis.tag_fusion_suitability": "tag_fusion_suitability_metrics.csv",
                      "analysis.tag_fusion_context": "tag_fusion_suitability_metrics.csv",
                      "analysis.sequence_liability": "sequence_liability_metrics.csv",
                      "analysis.binder_interface_chemistry": "binder_interface_chemistry_metrics.csv"}[workflow_id]
        metrics_path = run_dir / "artifacts" / "analysis" / table_name
        if status in {"completed", "succeeded"} and metrics_path.is_file():
            try:
                frame = pd.read_csv(metrics_path)
            except (OSError, ValueError, pd.errors.ParserError):
                frame = pd.DataFrame()
            if not frame.empty:
                source_column = "source_candidate_id" if "source_candidate_id" in frame else "analysis_source_candidate_id"
                if source_column in frame:
                    selected_frame = frame[frame[source_column].astype(str) == source_candidate_id]
                    if not selected_frame.empty:
                        frame = selected_frame
                st.markdown("##### Measurements saved by this evaluation")
                st.dataframe(frame, hide_index=True, use_container_width=True)
        if workflow_id == "analysis.pep_patch":
            for path in sorted((run_dir / "artifacts" / "analysis").glob("item_*/patch_metrics.json")):
                metrics = _read_json(path)
                if metrics:
                    with st.expander("PEP-Patch surface metrics"):
                        st.json(metrics)
            for path in sorted((run_dir / "artifacts" / "analysis").glob("item_*/residue_patches*.csv")):
                try:
                    frame = pd.read_csv(path, index_col=0)
                except (OSError, ValueError, pd.errors.ParserError):
                    continue
                with st.expander(f"Residue patch output · {path.name}"):
                    st.dataframe(frame, hide_index=False, use_container_width=True)
        if ligand_workflow and status in {"completed", "succeeded"}:
            item = evaluation.get("selected_item") or {}
            item_dir = run_dir / "artifacts" / "analysis" / str(item.get("id") or "")
            csv_name = "full_checks.csv" if workflow_id == "analysis.posebusters" else "interactions.csv"
            artifact_path = item_dir / csv_name
            if artifact_path.is_file():
                try:
                    frame = pd.read_csv(artifact_path)
                except (OSError, ValueError, pd.errors.ParserError):
                    frame = pd.DataFrame()
                with st.expander("Pose validation checks" if workflow_id == "analysis.posebusters" else "Ligand–receptor interaction records", expanded=True):
                    if frame.empty:
                        st.caption("The registered evaluator produced no tabular rows for this item.")
                    else:
                        st.dataframe(frame, hide_index=True, use_container_width=True)
        if workflow_id == "analysis.reference_placement" and status in {"completed", "succeeded"}:
            placement = _read_json(run_dir / "artifacts" / "analysis" / "reference_placement" / "reference_placement.json")
            if placement:
                with st.expander("Rigid transformation, exact structures, mapped residues and limitations"):
                    st.json(placement)
            summary = _read_json(run_dir / "artifacts" / "analysis" / "reference_placement_summary.json")
            if summary:
                with st.expander("Registered evaluator summary", expanded=False):
                    st.json(summary)
        if workflow_id == "analysis.reference_environment" and status in {"completed", "succeeded"}:
            result = _read_json(run_dir / "artifacts" / "analysis" / "reference_environment" / "reference_environment.json")
            if result:
                st.dataframe([result.get("measurements", {})], hide_index=True, use_container_width=True)
                with st.expander("Residue/atom evidence, method and limitations", expanded=True):
                    st.json(result)
        if workflow_id in {"analysis.tag_fusion_suitability", "analysis.tag_fusion_context"} and status in {"completed", "succeeded"}:
            result = _read_json(run_dir / "artifacts" / "analysis" / "tag_fusion_suitability.json")
            if result:
                st.markdown(f"**Assessment status:** `{result.get('assessment_status', 'unknown')}` · **Scope:** `{result.get('assessment_scope', 'unknown')}`")
                if result.get("not_applicable_reason"):
                    st.info(result["not_applicable_reason"])
                with st.expander("Terminal evidence, tag/fusion clearance and scientific limitations", expanded=True):
                    st.json(result)
                if result.get("adaptation_attribution"):
                    st.caption(result["adaptation_attribution"])
        if workflow_id == "analysis.sequence_liability" and status in {"completed", "succeeded"}:
            result = _read_json(run_dir / "artifacts" / "analysis" / "sequence_liability" / "sequence_liability.json")
            if result:
                st.markdown(f"**Assessment status:** `{result.get('assessment_status', 'unknown')}` · **Method:** `{result.get('method', 'unavailable')}`")
                with st.expander("Residue-linked findings, measurements and limitations", expanded=True):
                    st.json(result)
        if workflow_id == "analysis.binder_interface_chemistry" and status in {"completed", "succeeded"}:
            result = _read_json(run_dir / "artifacts" / "analysis" / "binder_interface_chemistry.json")
            if result:
                st.dataframe([result.get("metrics", {})], hide_index=True, use_container_width=True)
                with st.expander("Residue annotations, hydrogen-bond geometry, settings and assumptions", expanded=True):
                    st.json(result)
    evaluator = {"analysis.ipsae": "ipsae", "analysis.rosetta": "rosetta_metrics", "analysis.pep_patch": "pep_patch",
                 "analysis.plip": "plip", "analysis.pandamap": "pandamap", "analysis.posebusters": "posebusters",
                 "analysis.reference_environment": "reference_environment",
                 "analysis.tag_fusion_suitability": "tag_fusion_suitability",
                 "analysis.tag_fusion_context": "tag_fusion_context",
                 "analysis.sequence_liability": "sequence_liability",
                 "analysis.binder_interface_chemistry": "binder_interface_chemistry_v1"}.get(workflow_id, workflow_id)
    evaluation_id = _analysis_evaluation_id(run_dir, user, source_candidate_id, evaluator=evaluator)
    if evaluation_id:
        st.success(f"Saved Evaluation · {evaluation_id}")
    with st.expander("Analysis provenance", expanded=False):
        if workflow_id == "analysis.pep_patch":
            provenance = evaluation.get("manifest") or {}
            selection = provenance.get("selection") or {}
            item = evaluation.get("selected_item") or {}
            scientific_parameters = provenance.get("settings") or {}
        elif ligand_workflow:
            provenance = evaluation.get("manifest") or {}
            selection = provenance.get("selection") or {}
            item = evaluation.get("selected_item") or {}
            scientific_parameters = provenance.get("settings") or {}
        else:
            provenance = _read_json(run_dir / "artifacts" / "analysis" / "complex_analysis_selection.json")
            selection = provenance
            item = provenance
            scientific_parameters = provenance.get("parameters") or {}
        metadata = evaluation.get("metadata") or {}
        database_source = metadata.get("database_source") if isinstance(metadata.get("database_source"), dict) else {}
        selection_database_source = selection.get("database_source") if isinstance(selection.get("database_source"), dict) else {}
        st.json({
            "evaluation_run": f"{evaluation.get('task_group')}/{evaluation.get('run_id')}",
            "workflow": workflow_id,
            "source_kind": metadata.get("upstream_source_kind") or "results",
            "source_job": (f"{metadata.get('upstream_task_group')}/{metadata.get('upstream_run_id')}"
                           if metadata.get("upstream_source_kind") != "database" else ""),
            "database_source": {key: (selection_database_source or database_source).get(key) for key in
                                ("molecule_id", "structure_id", "group_id", "evaluation_id", "model_id", "structure_digest")
                                if key in (selection_database_source or database_source)},
            "source_candidate_id": item.get("candidate_id"),
            "selected_structure_sha256": item.get("source_structure_digest") or item.get("structure_digest"),
            "analysis_complex_sha256": item.get("structure_digest") if ligand_workflow else "",
            "model_id": selection.get("model_id") or item.get("model_id"),
            "ligand_selection": selection.get("ligand_selection") or item.get("ligand_selection") or {},
            "ligand_identity": item.get("ligand_identity") or {},
            "ligand_pose_sha256": item.get("ligand_pose_digest") or "",
            "receptor_coordinates_sha256": item.get("receptor_digest") or "",
            "binder_chains": item.get("binder_chains") or [],
            "target_chains": item.get("target_chains") or [],
            "method": provenance.get("tool") or provenance.get("workflow"),
            "parameters": scientific_parameters,
            "source_job_code": metadata.get("upstream_job_code"),
            "evaluation_id": evaluation_id,
        })


def _render_contextual_workflow(user: dict, *, workflow_id: str, schema: dict | None,
                                job: dict, selected: dict, candidate: dict,
                                structure_path: Path, model_id: str, candidate_file: Path,
                                task_group: str, run_id: str, project_id: str,
                                binder_chains: list[str], target_chains: list[str],
                                scientific_parameters: dict, key: str, submit,
                                unavailable: str = "", ligand_selection: dict | None = None) -> None:
    title = {"analysis.ipsae": "ipSAE & Structure Analysis",
             "analysis.rosetta": "Rosetta Interface Analysis",
             "analysis.pep_patch": "PEP-Patch", "analysis.plip": "PLIP interaction analysis",
             "analysis.pandamap": "PandaMap interaction analysis",
             "analysis.posebusters": "PoseBusters pose validation"}.get(workflow_id, workflow_id)
    st.markdown(f"##### {title}")
    with st.container():
        if not schema:
            st.info(unavailable or "The registered workflow is unavailable in this page context.")
            return
        if unavailable:
            st.info(unavailable)
            return
        if not candidate.get("candidate_id") or not candidate_file.is_file():
            st.info("This Results artifact is not linked to a normalized source candidate; the registered analysis cannot preserve its source lineage.")
            return
        if not task_group or not run_id or not str(job.get("project_id") or ""):
            st.info("The source Results job has no complete accessible project identity for a new workflow.")
            return
        if not project_id or project_id != str(job.get("project_id") or ""):
            st.info("Select the source Results job's project so the workflow remains within its access boundary.")
            return

        digest = _structure_digest(structure_path)
        candidate_id = str(candidate.get("candidate_id") or "")
        is_pep_patch = workflow_id == "analysis.pep_patch"
        existing = _matching_contextual_jobs(
            user, workflow_id=workflow_id, source_task_group=task_group, source_run_id=run_id,
            candidate_id=candidate_id, structure_digest=digest, model_id=model_id,
            binder_chains=binder_chains, target_chains=target_chains,
            parameters=scientific_parameters, ligand_selection=ligand_selection,
        ) if digest else []
        current = next((item for item in existing if item["status"] in {"completed", "succeeded"}), None)
        if current is None:
            current = next((item for item in existing if item["status"] not in {"failed", "cancelled", "canceled"}), None)
        if current:
            _render_contextual_results(user, current, workflow_id=workflow_id,
                                       source_candidate_id=candidate_id, key=key)
        elif existing:
            failed = existing[0]
            st.error(f"Previous matching analysis failed: {failed.get('error') or 'No per-item error was recorded.'}")
            _render_contextual_results(user, failed, workflow_id=workflow_id,
                                       source_candidate_id=candidate_id, key=key)
        else:
            st.caption("No matching evaluation is saved for this exact structure, model, chain selection, and method settings.")

        if current and current["status"] in {"completed", "succeeded"}:
            return
        button_label = "Retry analysis" if existing else "Run analysis"
        if st.button(button_label, type="primary", key=f"complex-analysis-run-{workflow_id.replace('.', '-')}-{key}", disabled=submit is None):
            parameters = {
                **scientific_parameters,
                "candidate_ids": [candidate_id],
                "selected_structure_path": _analysis_structure_reference(job, selected, structure_path),
                "selected_model_id": model_id,
            }
            if workflow_id in {"analysis.plip", "analysis.pandamap", "analysis.posebusters"}:
                parameters["ligand_selection"] = dict(ligand_selection or {})
            elif is_pep_patch:
                parameters["binder_chains"] = list(binder_chains)
            else:
                parameters["binder_chains"] = list(binder_chains)
                parameters["target_chains"] = list(target_chains)
            source_job = _analysis_source_reference(job)
            if not job.get("cluster_job_id"):
                st.session_state["complex-analysis-active"] = {
                    "kind": "results", "task_group": task_group, "run_id": run_id,
                    "candidate_id": candidate_id, "model_index": model_id,
                    "relative_path": selected["relative_path"],
                }
            before = st.session_state.get("last_submitted_job")
            try:
                submit(schema, workflow_id, {**parameters, "source_job": source_job},
                       project_id, {}, title)
            except (OSError, RuntimeError, TypeError, ValueError) as exc:
                st.error(f"Could not start {title}: {exc}")
                return
            if st.session_state.get("last_submitted_job") != before:
                st.rerun()


def _matching_database_analysis_jobs(user: dict, *, workflow_id: str, source: dict,
                                    candidate_id: str, model_id: str, binder_chains: list[str],
                                    target_chains: list[str], parameters: dict,
                                    ligand_selection: dict | None = None, indexed_rows=None,
                                    load_evidence=True, expected_protocol=None, allow_reassessment=False) -> list[dict]:
    """Find exact indexed work; full evidence is read only for targeted inspection."""
    from mn_molecular_design_studio.auth import can_access_job
    from mn_molecular_design_studio.job_catalog import database_analysis_jobs, analysis_source_key
    from mn_molecular_design_studio.complex_analysis_batch import scientific_parameters, protocol_matches, assessment_protocol
    require_verified_protocol = expected_protocol is None
    expected_protocol = assessment_protocol(workflow_id) if expected_protocol is None else expected_protocol
    tool = {"analysis.biotite_interface": "biotite_interface", "analysis.ipsae": "ipsae",
            "analysis.rosetta": "rosetta_metrics", "analysis.pep_patch": "pep_patch",
            "analysis.sequence_qc": "sequence_qc", "analysis.proteinsol": "proteinsol",
            "analysis.plip": "plip", "analysis.pandamap": "pandamap", "analysis.posebusters": "posebusters",
            "analysis.reference_placement": "reference_placement",
            "analysis.reference_environment": "reference_environment",
            "analysis.tag_fusion_suitability": "tag_fusion_suitability",
            "analysis.tag_fusion_context": "tag_fusion_context",
            "analysis.sequence_liability": "sequence_liability",
            "analysis.binder_interface_chemistry": "binder_interface_chemistry_v1"}.get(workflow_id)
    if not tool:
        return []
    rows = database_analysis_jobs(user, [source]) if indexed_rows is None else indexed_rows
    key = analysis_source_key(source)
    identity_keys = ("molecule_id", "structure_id", "group_id", "structure_digest", "model_id")
    if workflow_id == "analysis.ipsae" and source.get("pae_digest"):
        identity_keys += ("pae_digest",)
    if workflow_id == "analysis.sequence_liability":
        identity_keys += ("ligand_selection",)
    if workflow_id == "analysis.reference_environment":
        identity_keys += ("reference_context_id", "reference_molecule_id", "reference_structure_id",
                          "reference_structure_digest", "reference_model_id", "reference_chain_id", "placement_run_id")
    matches = []
    for row in rows:
        projection = row.get("database_analysis") or {}
        if row.get("tool") != tool or projection.get("source_key") != key or projection.get("invalid"):
            continue
        if not can_access_job(user, row):
            continue
        saved_source = projection.get("source") or {}
        if any((saved_source.get(field) or {}) != (source.get(field) or {}) if field == "ligand_selection"
               else saved_source.get(field, "") != source.get(field, "") for field in identity_keys):
            continue
        if (projection.get("candidate_id") != candidate_id or projection.get("model_id") != str(model_id)
                or sorted(map(str, projection.get("binder_chains") or [])) != sorted(map(str, binder_chains))
                or sorted(map(str, projection.get("target_chains") or [])) != sorted(map(str, target_chains))):
            continue
        if workflow_id in {"analysis.plip", "analysis.pandamap", "analysis.posebusters"} and projection.get("ligand_selection") != (ligand_selection or {}):
            continue
        saved_parameters = scientific_parameters(projection.get("parameters") or {}, workflow_id)
        requested = scientific_parameters(parameters, workflow_id)
        if workflow_id == "analysis.reference_environment":
            fields = ("source_component_chain", "selected_reference_chains", "assembly_id")
            if any(saved_parameters.get(field) != requested.get(field) for field in fields):
                continue
        elif saved_parameters != requested:
            continue
        status = str(row.get("status") or "unknown").casefold()
        if expected_protocol is not None and status in {"completed", "succeeded", "failed", "cancelled", "canceled", "error"}:
            saved_protocol = projection.get("protocol") or {}
            if not protocol_matches(saved_protocol, expected_protocol):
                known = bool(expected_protocol) and all(v not in (None, "") and saved_protocol.get(k) not in (None, "")
                                                       for k, v in expected_protocol.items())
                if known or allow_reassessment or require_verified_protocol:
                    continue
        if status in {"completed", "succeeded"} and projection.get("item_status"):
            status = str(projection["item_status"]).casefold()
        result = {**row, "status": status, "error": projection.get("error") or ""}
        matches.append(result)
    matches.sort(key=lambda item: (0 if item["status"] in {"completed", "succeeded"} else
                 2 if item["status"] in {"failed", "cancelled", "canceled", "error"} else 1,
                 str(item.get("created_at") or ""), str(item.get("run_id") or "")))
    if load_evidence:
        # These callers have explicitly opened an individual structure/analysis.
        # Coverage planning passes False and never reads any Results files.
        for match in matches:
            run = Path(match["run_dir"])
            batch = _read_json(run / "artifacts/analysis/batch_manifest.json")
            match["manifest"] = batch or _read_json(run / "artifacts/analysis/complex_analysis_selection.json")
            match["metadata"] = _read_json(run / "metadata.json")
            match["selected_item"] = next((item for item in batch.get("items", [])
                if str(item.get("candidate_id") or "") == candidate_id), match["manifest"])
    return matches


def _render_biotite_workspace(user: dict, job: dict, selected: dict, components: list[dict], *,
                              workflow_schema, project_selector, submit, viewer_key: str) -> None:
    st.markdown("#### Available analyses")
    protein_chains = _protein_chain_ids(components)
    protein_ligand = any(item.get("type") == "small_molecule" for item in components) and bool(protein_chains)
    nucleic_acids = any(item.get("type") in {"dna", "rna", "nucleic_acid"} for item in components)
    ppi_applicable = len(protein_chains) >= 2
    st.markdown("##### Protein–protein")
    if ppi_applicable:
        st.write("**Biotite Interface** · existing Studio workflow")
        st.caption("Heavy-atom contacts and residue maps; geometric salt bridges; hydrogen bonds when explicit hydrogens are present. The result is structural evidence, not an experimental binding claim.")
    else:
        st.info("Biotite Interface is not applicable: this selected structure does not contain at least two protein chains.")
    if ppi_applicable:
        st.caption("Rosetta Interface Analysis and ipSAE are available below through their registered Studio workflows. They currently require one selected chain per group; ipSAE-derived scores additionally require a saved pairwise PAE matrix.")
    st.caption("Compare saved complexes and the Collection Map expose compatible saved structures and registered Interface Clustering results. Optional reference mapping is available for Database components; spatial placement and general counter-target specificity remain unavailable.")
    if protein_chains:
        st.caption("PEP-Patch evaluates hydrophobic and pH-specific electrostatic surface patches on a user-selected protein chain group.")
    st.markdown("##### Protein–ligand")
    if protein_ligand:
        st.caption("PLIP and PandaMap measure ligand–receptor interactions. PoseBusters evaluates ligand-pose geometry and chemistry; its checks are not interaction scores.")
    else:
        st.caption("Not applicable to the selected components.")
    st.markdown("##### Protein–DNA/RNA")
    if nucleic_acids:
        st.caption("The selected structure contains nucleic acid. No registered contextual analysis in this Biotite Interface slice evaluates protein–nucleic-acid interactions; capabilities will be listed as unavailable until a compatible existing workflow is connected.")
    else:
        st.caption("Not applicable to the selected components.")
    if not protein_chains:
        return

    candidate = selected.get("candidate") if isinstance(selected.get("candidate"), dict) else {}
    candidate_id = str(candidate.get("candidate_id") or "")
    path = Path(selected["path"]).expanduser().resolve()
    task_group = str(job.get("task_group") or "")
    run_id = str(job.get("run_id") or "")
    source_run_dir = Path(str(job.get("run_dir") or "")).expanduser().resolve()
    candidate_file = source_run_dir / "artifacts" / "normalized_candidates" / "candidates.jsonl"
    source_project = str(job.get("project_id") or "")
    candidate_metadata = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
    entry = candidate_metadata.get("complex_analysis_entry") or {}
    if entry:
        st.caption(f"Analysis entry: {entry.get('entry_mode')} · source candidate {entry.get('native_candidate_id')} · "
                   f"model {entry.get('model_id') or 'unspecified'} · {entry.get('method')}")
        with st.expander("Original candidate and structure provenance", expanded=False):
            st.json({key: entry.get(key) for key in (
                "candidate_id", "native_candidate_id", "variant_id", "parents", "parent_variant_ids", "target",
                "source_run_dir", "source_run", "source_structure", "structure_sha256", "model_id", "prediction_id", "method", "acceptance")})
            origin = Path(str(entry.get("source_run_dir") or ""))
            from mn_molecular_design_studio.auth import can_access_job
            from mn_molecular_design_studio.job_catalog import get_job_access_record
            if origin.name and can_access_job(user, get_job_access_record(origin.parent.name, origin.name)):
                if st.button("Return to original Results", key="complex-analysis-return-original"):
                    st.session_state["selected-job-reference"] = f"{origin.parent.name}/{origin.name}"
                    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
                    original_rows = read_candidates(origin / "artifacts/normalized_candidates/candidates.jsonl")
                    index = next((index for index, item in enumerate(original_rows)
                                  if item.get("candidate_id") == entry.get("native_candidate_id")), None)
                    if index is not None:
                        st.session_state[f"job-candidate-{origin.name}"] = index
                        st.session_state[f"cofolding-prediction-{origin.name}"] = index
                        st.session_state[f"job-design-scope-{origin.name}"] = "All"
                    st.session_state["page"] = "Job details"
                    st.rerun()
    candidate_binder = [str(value) for value in candidate.get("binder_chains") or candidate_metadata.get("binder_chains") or []]
    candidate_target = [str(value) for value in candidate.get("target_chains") or candidate_metadata.get("target_chains") or []]
    candidate_binder = [value for value in candidate_binder if value in protein_chains]
    candidate_target = [value for value in candidate_target if value in protein_chains and value not in candidate_binder]
    key = hashlib.sha256(f"{task_group}/{run_id}:{selected['relative_path']}".encode()).hexdigest()[:14]
    model_id = str(next((value for value in (candidate.get('model_id'), candidate.get('model_index'),
        candidate_metadata.get('model_id'), candidate_metadata.get('model_index'))
        if value is not None and str(value) != ''), ''))

    if ppi_applicable:
        biotite_id = "analysis.biotite_interface"
        biotite_schema = None
        biotite_unavailable = ""
        try:
            if workflow_schema is None:
                raise RuntimeError("The workflow registry is unavailable in this page context.")
            biotite_schema = workflow_schema(biotite_id)
        except Exception as exc:
            biotite_unavailable = f"The registered workflow could not be loaded: {exc}"
        if not candidate_id:
            biotite_unavailable = "This structure artifact has no normalized candidate identity in the selected Results job."
        elif not candidate_file.is_file():
            biotite_unavailable = "The selected Results job has no normalized candidate collection."
        elif _structure_suffix(path) not in {".pdb", ".ent", ".cif", ".mmcif"}:
            biotite_unavailable = "Biotite Interface requires a PDB or mmCIF protein complex artifact."
        elif not task_group or not run_id:
            biotite_unavailable = "The source Results job identity is incomplete."
        elif not source_project:
            biotite_unavailable = "The source Results job has no project assignment for a new workflow."
        st.markdown("##### Interface chain groups")
        st.caption(f"Exact model/artifact: {candidate_id or path.name} · {model_id or 'model not specified'} · {path.name}")
        first, second = st.columns(2)
        a_key, b_key = f"complex-analysis-interface-a-{key}", f"complex-analysis-interface-b-{key}"
        if a_key not in st.session_state:
            st.session_state[a_key] = candidate_binder
        if b_key not in st.session_state:
            st.session_state[b_key] = candidate_target
        binder_chains = first.multiselect("Interface chain group A", protein_chains, key=a_key,
            help="This group name does not assert design or experimental binding.")
        target_chains = second.multiselect("Interface chain group B", protein_chains, key=b_key,
            help="May contain one or more protein chains.")
        schema_parameters = {str(row.get("name")): row for row in (biotite_schema or {}).get("parameters", [])}
        default_distance = float(schema_parameters.get("contact_distance", {}).get("default", 5.0))
        contact_distance = st.number_input("Biotite contact distance (Å)", min_value=0.1, max_value=20.0,
            value=default_distance, step=0.1, key=f"complex-analysis-biotite-cutoff-{key}")
        invalid_groups = (not binder_chains or not target_chains or bool(set(binder_chains) & set(target_chains)))
        if invalid_groups:
            st.info("Choose two nonempty, disjoint protein chain groups to enable the interface analyses.")
        project_prefix = f"complex-analysis-{key}"
        if job.get("cluster_job_id"):
            try:
                source_project_row = _authorized_remote_source_project(job)
            except (ValueError, RuntimeError, OSError) as exc:
                source_project_row = None
                st.info(f"Could not load the authorized source project from the coordinator: {exc}")
            if source_project_row is None:
                analysis_project = None
                st.info("The source project is not present in the coordinator's authorized project list.")
            else:
                project_key = f"complex-analysis-{key}-project-id"
                allowed_project_ids = {"", source_project}
                if st.session_state.get(project_key) not in allowed_project_ids:
                    st.session_state[project_key] = ""
                analysis_project = st.selectbox(
                    "Analysis project *",
                    options=["", source_project],
                    format_func=lambda value: "Select a project" if value == "" else str(source_project_row.get("name") or value),
                    key=project_key,
                    help="Only the authorized project that owns this remote Results source can receive its analysis.",
                ) or None
        elif project_selector:
            project_key = f"{project_prefix}-project-id"
            if project_key not in st.session_state:
                st.session_state[project_key] = source_project
            analysis_project = project_selector(project_prefix, label="Analysis project", allow_create=False)
        else:
            analysis_project = source_project
            st.caption(f"Analysis project: {source_project or 'unavailable'}")
        if analysis_project != source_project:
            st.warning("Select the source Results job's project. Studio keeps analyses in the source project's access boundary.")

        digest = _structure_digest(path)
        contact_distance = round(float(contact_distance), 3)
        remote_measurements = {'entries': [], 'pending': False, 'errors': []}
        if job.get('cluster_job_id') and not invalid_groups and digest:
            from mn_molecular_design_studio.complex_analysis_results_ui import saved_results_measurements
            from mn_molecular_design_studio.cluster.transport import ClusterUnavailable
            try:
                saved = saved_results_measurements(user, job, [{
                    'candidate': {**candidate, 'binder_chains': binder_chains, 'target_chains': target_chains},
                    'digest': digest, 'label': candidate_id}],
                    token=str(st.session_state.get('studio-session-token') or 'development-mode'))
                remote_measurements = {**saved, 'entries': [row for row in saved['entries']
                    if row['method'] == 'biotite_interface'
                    and row['parameters'].get('contact_distance') == contact_distance
                    and row['status'] in {'measured', 'recorded', 'valid', 'completed', 'succeeded'}]}
            except (OSError, ValueError, ClusterUnavailable) as exc:
                remote_measurements['errors'] = [{'error': str(exc)}]
        existing = _matching_biotite_jobs(user, source_task_group=task_group, source_run_id=run_id,
            candidate_id=candidate_id, structure_digest=digest, binder_chains=binder_chains,
            target_chains=target_chains, contact_distance=contact_distance, model_id=model_id) if not job.get('cluster_job_id') and not invalid_groups and digest else []
        current = next((item for item in existing if item["status"] in {"completed", "succeeded"}), None)
        if current is None:
            current = next((item for item in existing if item["status"] not in {"failed", "cancelled", "canceled"}), None)
        if current:
            _render_biotite_results(user, current, source_candidate_id=candidate_id, key=key,
                structure_path=path, structure_context=candidate_metadata,
                structure_components=components,
                viewer_key=viewer_key)
        elif existing:
            failed = existing[0]
            st.error(f"Previous matching analysis failed: {failed.get('error') or 'No per-item error was recorded.'}")
            _render_biotite_results(user, failed, source_candidate_id=candidate_id, key=key,
                structure_path=path, structure_context=candidate_metadata,
                structure_components=components,
                viewer_key=viewer_key)
        elif remote_measurements['entries']:
            st.success('Saved Biotite Interface measurements are available for this exact source structure and chain selection.')
            st.dataframe([{'Evaluation ID': row['evaluation_id'], 'Source workstation': row['source_workstation'],
                           'Contact distance (Å)': row['parameters']['contact_distance'],
                           **row['metrics']} for row in remote_measurements['entries']], hide_index=True, width='stretch')
            st.caption('Open Results Evaluation Board for the saved measurements and full source provenance.')
        elif remote_measurements['pending'] or remote_measurements['errors']:
            st.info('Saved remote measurements are loading or unavailable. Refresh before submitting an additional analysis.')
            for error in remote_measurements['errors']:
                st.warning(error['error'])
            if st.button('Refresh saved measurements', key=f'complex-analysis-remote-measurements-{key}'):
                st.rerun()
        elif not biotite_unavailable:
            st.caption("No matching Biotite Interface evaluation is saved for this exact structure, model, chain mapping, and contact distance.")
        if biotite_unavailable:
            st.warning(f"Biotite Interface unavailable: {biotite_unavailable}")
        elif not (current and current["status"] in {"completed", "succeeded"}) and not remote_measurements['entries']:
            button_label = "Retry analysis" if existing else "Run analysis"
            disabled = bool(invalid_groups or not analysis_project or analysis_project != source_project or submit is None or biotite_schema is None
                            or remote_measurements['pending'] or remote_measurements['errors'])
            if st.button(button_label, type="primary", disabled=disabled, key=f"complex-analysis-run-biotite-{key}"):
                if not job.get("cluster_job_id"):
                    st.session_state["complex-analysis-active"] = {
                        "kind": "results", "task_group": task_group, "run_id": run_id,
                        "candidate_id": candidate_id, "model_index": model_id,
                        "relative_path": selected["relative_path"],
                    }
                workflow_parameters = {
                    "candidate_ids": [candidate_id],
                    "selected_structure_path": _analysis_structure_reference(job, selected, path),
                    "selected_model_id": model_id, "binder_chains": list(binder_chains),
                    "target_chains": list(target_chains), "contact_distance": contact_distance,
                }
                before = st.session_state.get("last_submitted_job")
                try:
                    submit(biotite_schema, biotite_id, {**workflow_parameters, "source_job": _analysis_source_reference(job)},
                           str(analysis_project), {}, "Biotite Interface")
                except (OSError, RuntimeError, TypeError, ValueError) as exc:
                    st.error(f"Could not start Biotite Interface: {exc}")
                    return
                if st.session_state.get("last_submitted_job") != before:
                    st.rerun()
        pep_surface_key = f"complex-analysis-pep-patch-chains-{key}"
        if pep_surface_key not in st.session_state:
            st.session_state[pep_surface_key] = candidate_binder or binder_chains
        pep_surface_chains = st.multiselect("PEP-Patch surface chain group", protein_chains,
            key=pep_surface_key, help="Choose the protein surface whose patches to evaluate.")
    else:
        binder_chains = []
        target_chains = []
        st.markdown("##### Protein surface group")
        surface_key = f"complex-analysis-pep-patch-chains-{key}"
        if surface_key not in st.session_state:
            st.session_state[surface_key] = candidate_binder or protein_chains[:1]
        pep_surface_chains = st.multiselect("Protein chains for PEP-Patch", protein_chains, key=surface_key)
        binder_chains = list(pep_surface_chains)
        analysis_project = source_project
        if project_selector:
            project_prefix = f"complex-analysis-{key}"
            project_key = f"{project_prefix}-project-id"
            if project_key not in st.session_state:
                st.session_state[project_key] = source_project
            analysis_project = project_selector(project_prefix, label="Analysis project", allow_create=False)

    common_unavailable = ""
    if not candidate_id:
        common_unavailable = "This structure has no normalized candidate identity in its Results job."
    elif not candidate_file.is_file():
        common_unavailable = "The source Results job has no normalized candidate collection."
    elif not source_project or not task_group or not run_id:
        common_unavailable = "The source Results job lacks the project or stable job identity required for this workflow."
    elif _structure_suffix(path) not in {".pdb", ".ent", ".cif", ".mmcif"}:
        common_unavailable = "The selected analysis requires a PDB or mmCIF protein structure."

    schemas = {}
    schema_errors = {}
    for workflow_id in ("analysis.rosetta", "analysis.ipsae", "analysis.pep_patch",
                        "analysis.plip", "analysis.pandamap", "analysis.posebusters"):
        try:
            schemas[workflow_id] = workflow_schema(workflow_id) if workflow_schema else None
            if not schemas[workflow_id]:
                schema_errors[workflow_id] = "The workflow registry did not return this analysis schema."
        except Exception as exc:
            schema_errors[workflow_id] = f"The registered workflow is unavailable: {exc}"

    if ppi_applicable:
        role_unavailable = common_unavailable or ("Choose two nonempty, disjoint protein chain groups." if invalid_groups else "")
        single_chain_unavailable = role_unavailable or (
            "The registered workflow supports one selected chain in each group." if len(binder_chains) != 1 or len(target_chains) != 1 else ""
        )
        st.markdown("#### Protein–protein evaluations")
        _render_contextual_workflow(user, workflow_id="analysis.rosetta", schema=schemas.get("analysis.rosetta"),
            job=job, selected=selected, candidate=candidate, structure_path=path, model_id=model_id,
            candidate_file=candidate_file, task_group=task_group, run_id=run_id, project_id=str(analysis_project or ""),
            binder_chains=list(binder_chains), target_chains=list(target_chains), scientific_parameters={},
            key=key, submit=submit, unavailable=schema_errors.get("analysis.rosetta") or single_chain_unavailable)
        ipsae_schema = schemas.get("analysis.ipsae")
        ipsae_parameters = {str(row.get("name")): row for row in (ipsae_schema or {}).get("parameters", [])}
        with st.expander("Configure ipSAE parameters", expanded=False):
            defaults = {name: float(ipsae_parameters.get(name, {}).get("default", value)) for name, value in {
                "pae_cutoff": 10., "dist_cutoff": 10., "ipdae_contact_cutoff": 8.,
                "structure_contact_cutoff": 4., "hotspot_contact_cutoff": 5.,
            }.items()}
            fields = st.columns(3)
            ipsae_values = {}
            for index, (name, default) in enumerate(defaults.items()):
                ipsae_values[name] = round(float(fields[index % len(fields)].number_input(
                    name.replace("_", " ").title() + (" (Å)" if name != "pae_cutoff" else " (Å)"),
                    min_value=0.1, max_value=50.0, value=default, step=0.5,
                    key=f"complex-analysis-ipsae-{name}-{key}")), 3)
        candidate_pae = str(candidate_metadata.get("pae_path") or "")
        pae_path = Path(candidate_pae).expanduser() if candidate_pae else None
        if pae_path is not None and not pae_path.is_absolute():
            pae_path = source_run_dir / pae_path
        if pae_path is None or not pae_path.is_file():
            st.caption("No saved pairwise PAE matrix was found for this model. The existing ipSAE workflow can still report structure metrics, while PAE-derived scores will be marked unavailable.")
        _render_contextual_workflow(user, workflow_id="analysis.ipsae", schema=ipsae_schema,
            job=job, selected=selected, candidate=candidate, structure_path=path, model_id=model_id,
            candidate_file=candidate_file, task_group=task_group, run_id=run_id, project_id=str(analysis_project or ""),
            binder_chains=list(binder_chains), target_chains=list(target_chains), scientific_parameters=ipsae_values,
            key=key, submit=submit, unavailable=schema_errors.get("analysis.ipsae") or single_chain_unavailable)

    if protein_chains:
        pep_schema = schemas.get("analysis.pep_patch")
        pep_parameters = {str(row.get("name")): row for row in (pep_schema or {}).get("parameters", [])}
        default_ph = float(pep_parameters.get("ph", {}).get("default", 7.4))
        ph = st.number_input("PEP-Patch pH", min_value=0.0, max_value=14.0, value=default_ph, step=0.1,
                             key=f"complex-analysis-pep-patch-ph-{key}")
        surface_unavailable = common_unavailable or ("Choose at least one protein chain for the surface calculation." if not pep_surface_chains else "")
        if analysis_project != source_project:
            surface_unavailable = "Select the source Results job's project to keep this workflow in its access boundary."
        _render_contextual_workflow(user, workflow_id="analysis.pep_patch", schema=pep_schema,
            job=job, selected=selected, candidate=candidate, structure_path=path, model_id=model_id,
            candidate_file=candidate_file, task_group=task_group, run_id=run_id, project_id=str(analysis_project or ""),
            binder_chains=list(pep_surface_chains), target_chains=[], scientific_parameters={"ph": round(float(ph), 3)},
            key=key, submit=submit, unavailable=schema_errors.get("analysis.pep_patch") or surface_unavailable)

    if protein_ligand:
        st.markdown("#### Protein–ligand analyses")
        ligand_components = [item for item in components if str(item.get("type") or "").casefold() == "small_molecule"]
        ligand_key = f"complex-analysis-results-ligand-{key}"
        ligand_options = [item.get("id") for item in ligand_components]
        ligand_labels = {item.get("id"): f"{item.get('label') or item.get('ccd') or 'Ligand'} · {item.get('residue_key') or item.get('source_chain')}"
                         for item in ligand_components}
        if st.session_state.get(ligand_key) not in ligand_options:
            st.session_state[ligand_key] = ligand_options[0] if ligand_options else ""
        selected_ligand_id = st.selectbox("Ligand residue", ligand_options,
            format_func=lambda value: ligand_labels[value], key=ligand_key)
        selected_ligand = next((item for item in ligand_components if item.get("id") == selected_ligand_id), {})
        ligand_selection = _ligand_selection_for_component(selected_ligand, ligand_components, candidate_metadata)
        ligand_unavailable = common_unavailable
        if not ligand_selection.get("residue_key") or not ligand_selection.get("chain"):
            ligand_unavailable = ligand_unavailable or "The selected ligand has no stable source chain and residue identity."
        elif not ligand_selection.get("smiles"):
            ligand_unavailable = ligand_unavailable or "The source does not preserve this ligand's isomeric SMILES identity. PLIP, PandaMap, and PoseBusters require a chemically matched ligand template; add a verified identity to the source record first."
        st.caption(f"Selected ligand: {ligand_selection.get('ccd') or 'SMILES-defined'} · {ligand_selection.get('residue_key') or 'residue identity unavailable'} · receptor model {model_id}")
        for workflow_id in ("analysis.plip", "analysis.pandamap"):
            _render_contextual_workflow(user, workflow_id=workflow_id, schema=schemas.get(workflow_id),
                job=job, selected=selected, candidate=candidate, structure_path=path, model_id=model_id,
                candidate_file=candidate_file, task_group=task_group, run_id=run_id, project_id=str(analysis_project or ""),
                binder_chains=[], target_chains=[], scientific_parameters={}, ligand_selection=ligand_selection,
                key=key, submit=submit, unavailable=schema_errors.get(workflow_id) or ligand_unavailable)
        _render_contextual_workflow(user, workflow_id="analysis.posebusters", schema=schemas.get("analysis.posebusters"),
            job=job, selected=selected, candidate=candidate, structure_path=path, model_id=model_id,
            candidate_file=candidate_file, task_group=task_group, run_id=run_id, project_id=str(analysis_project or ""),
            binder_chains=[], target_chains=[], scientific_parameters={}, ligand_selection=ligand_selection,
            key=key, submit=submit, unavailable=schema_errors.get("analysis.posebusters") or ligand_unavailable)


def _show_component_groups(components: list[dict], key: str) -> None:
    if not components:
        return
    st.markdown("#### Component groups")
    st.caption("Create any number of named groups. Group selection is for inspection only; it does not run an analysis.")
    group_ids_key = f"complex-analysis-group-ids-{key}"
    group_ids = list(st.session_state.get(group_ids_key) or [])
    if not group_ids:
        group_ids = ["group-1"]
        st.session_state[group_ids_key] = group_ids
    options = [item["id"] for item in components]
    labels = {item["id"]: f"{item['type']} · chain {item['source_chain']} · {item['label']}" for item in components}
    for group_id in list(group_ids):
        name_col, member_col, remove_col = st.columns([1.1, 3.2, .5], vertical_alignment="bottom")
        name_key = f"complex-analysis-group-name-{key}-{group_id}"
        member_key = f"complex-analysis-group-members-{key}-{group_id}"
        if name_key not in st.session_state:
            st.session_state[name_key] = f"Group {group_ids.index(group_id) + 1}"
        name_col.text_input("Group name", key=name_key, label_visibility="collapsed", placeholder="Group name")
        current = [value for value in st.session_state.get(member_key, []) if value in options]
        st.session_state[member_key] = current
        member_col.multiselect("Components", options, key=member_key, format_func=lambda value: labels[value], placeholder="Choose components")
        if len(group_ids) > 1 and remove_col.button("×", key=f"complex-analysis-group-remove-{key}-{group_id}", help="Remove this group"):
            group_ids.remove(group_id)
            st.session_state[group_ids_key] = group_ids
            st.rerun()
        selected = st.session_state.get(member_key, [])
        if selected:
            selection = [item for item in components if item["id"] in selected]
            st.caption(f"{st.session_state.get(name_key) or 'Group'}: " + ", ".join(
                f"{item['type']} chain {item['source_chain']}" for item in selection
            ))
            with st.expander(f"Inspect {st.session_state.get(name_key) or 'selected group'}", expanded=False):
                _show_components(selection)
    if st.button("Add component group", key=f"complex-analysis-group-add-{key}"):
        next_index = max((int(value.split("-")[-1]) for value in group_ids if value.split("-")[-1].isdigit()), default=0) + 1
        st.session_state[group_ids_key] = [*group_ids, f"group-{next_index}"]
        st.rerun()


def _show_structure(path: Path, context: dict, components: list[dict], *, key: str,
                    highlighted_selections: list[str] | None = None) -> None:
    if path.stat().st_size > VIEWER_MAX_BYTES:
        st.info("This structure is over 50 MB. Download it from its source results or Database record to inspect it locally.")
        return
    try:
        text = _viewer_text(path, context)
        from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer import molstar_custom_component
        from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer.dataclasses import StructureVisualization

        molstar_custom_component(
            [StructureVisualization(pdb=text, color="chain-id", representation_type="cartoon+ball-and-stick",
                                    highlighted_selections=highlighted_selections)],
            key=key, height=650, show_controls=True, show_structure_tools=False,
        )
    except (ImportError, OSError, ValueError, TypeError) as exc:
        st.error(f"Could not display this structure in Mol*: {exc}")


def _show_reference_placement_overlay(reference_path: Path, source_path: Path,
                                      source_context: dict, result: dict, *,
                                      key: str, highlighted_selections: list[str] | None = None,
                                      reference_highlighted_selections: list[str] | None = None) -> None:
    """Show the exact associated reference and transformed source in one Mol* view."""
    if max(reference_path.stat().st_size, source_path.stat().st_size) > VIEWER_MAX_BYTES:
        st.info("One of these structures is over 50 MB; the placement remains saved, but this paired Mol* view is unavailable.")
        return
    try:
        from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer import molstar_custom_component
        from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer.dataclasses import StructureVisualization

        visualizations = [
            StructureVisualization(
                pdb=_viewer_text(reference_path, {}), color="uniform",
                color_params={"value": "0x9AA3AE"},
                representation_type="cartoon+ball-and-stick",
                highlighted_selections=reference_highlighted_selections,
            ),
            StructureVisualization(
                pdb=_viewer_text(source_path, source_context), color="chain-id",
                representation_type="cartoon+ball-and-stick",
                highlighted_selections=highlighted_selections,
                transform_matrix=result["matrix_4x4"],
            ),
        ]
        molstar_custom_component(
            visualizations, key=key, height=650, show_controls=True,
            show_structure_tools=False,
        )
    except (ImportError, OSError, ValueError, TypeError, KeyError) as exc:
        st.error(f"Could not display the validated reference placement in Mol*: {exc}")


def _show_components(components: list[dict]) -> None:
    if not components:
        st.info("No supported components could be identified from this structure or its saved metadata.")
        return
    rows = [{
        "Component": item["label"],
        "Type": item["type"],
        "Structure chain": item["source_chain"],
        "Component ID": item["id"],
        "Role": item.get("role") or "—",
        "Sequence / identifier": item.get("sequence") or item.get("ccd") or item.get("smiles") or "—",
    } for item in components]
    st.dataframe(rows, hide_index=True, use_container_width=True)


def _show_human_review(database, user: dict, *, source_kind: str, subject_id: str,
                       structure_ref: str, evaluation_id: str = "") -> None:
    """Edit one account's durable review note, separate from scientific evaluations."""
    from mn_molecular_design_studio.complex_analysis_review_ui import render_review_editor, mark_review_saved
    identity = {"source_kind": source_kind, "subject_id": subject_id,
                "structure_ref": structure_ref, "evaluation_id": evaluation_id}
    try:
        review = database.get_complex_review(user, **identity)
    except (OSError, PermissionError, ValueError, TypeError) as exc:
        st.warning(f"Cannot load saved human review: {exc}")
        return
    st.markdown("#### Human review")
    st.caption("Your review is saved persistently under your account, separately from calculated evaluations.")
    decision, base_key, dirty = render_review_editor(user, identity, review,
        key_prefix="complex-analysis-review", status_label="Review status", note_label="Review notes",
        placeholder="Record observations or questions for a human reviewer.")
    if st.button("Save review", key=base_key.replace("-base-", "-save-")):
        if decision["status"] in {"Shortlisted for experimental testing", "Not shortlisted for experimental testing"} and not decision["note"].strip():
            st.error("Enter a short rationale for this shortlist decision.")
        elif dirty:
            try:
                saved = database.save_complex_review(user, **decision)
            except (OSError, PermissionError, ValueError, TypeError) as exc:
                st.error(str(exc))
            else:
                mark_review_saved(base_key, saved)
                st.success(f"Review saved · {saved['updated_at']}")
        else:
            st.info("No review changes to save.")


def _explicit_reference_link(structure: dict) -> tuple[str, str, str]:
    """Return only a reference pair explicitly carried by this representation's provenance."""
    context = structure.get("context") if isinstance(structure.get("context"), dict) else {}
    for molecule_key, structure_key in (
        ("design_target_molecule_id", "design_target_structure_id"),
        ("target_molecule_id", "target_structure_id"),
        ("reference_molecule_id", "reference_structure_id"),
    ):
        molecule_id = str(context.get(molecule_key) or "").strip()
        structure_id = str(context.get(structure_key) or "").strip()
        if structure_id:
            return molecule_id, structure_id, f"structure.context.{molecule_key}/{structure_key}"
    return "", "", ""


def _paired_reference_for_viewer(database, user: dict, record: dict, structure: dict,
                                 components: list[dict], key: str) -> dict | None:
    """Resolve the current user's digest-checked optional reference for the viewer.

    This only loads the exact saved Database representation. It does not place,
    align, transform, or score the reference against the selected complex.
    """
    protein_chains = sorted({str(item.get("source_chain") or "") for item in components
                             if item.get("type") == "protein" and item.get("source_chain")})
    if not protein_chains:
        return None
    component_key = f"complex-analysis-reference-component-{key}"
    source_chain = str(st.session_state.get(component_key) or "")
    if source_chain not in protein_chains:
        source_chain = protein_chains[0]
        st.session_state[component_key] = source_chain
    try:
        saved = database.get_complex_reference_context(
            user, source_molecule_id=str(record["id"]),
            source_structure_id=str(structure["id"]), source_chain=source_chain,
        )
        if not saved or not saved.get("configured"):
            return None
        reference_id = str(saved.get("reference_structure_id") or "")
        if not reference_id:
            return None
        reference_path = database.structure_file(user, reference_id)
        reference_components = component_inventory(reference_path.read_bytes(), reference_path.name)
    except (OSError, PermissionError, ValueError, TypeError, KeyError):
        return None
    return {
        "saved": saved,
        "source_chain": source_chain,
        "path": reference_path,
        "components": reference_components,
        "reference_id": reference_id,
    }


def _reference_placement_parameters(saved: dict, model_id: str) -> dict:
    from mn_molecular_design_studio.complex_analysis_placement import PLACEMENT_METHOD, PLACEMENT_VERSION

    return {
        "source_component_chain": str(saved.get("source_component_chain") or ""),
        "selected_model_id": str(model_id),
        "reference_context_id": str(saved.get("id") or ""),
        "reference_molecule_id": str(saved.get("reference_molecule_id") or ""),
        "reference_structure_id": str(saved.get("reference_structure_id") or ""),
        "reference_structure_digest": str(saved.get("reference_structure_digest") or ""),
        "reference_model_id": str(saved.get("reference_model_id") or ""),
        "reference_component_chain": str(saved.get("reference_component_chain") or ""),
        "method": PLACEMENT_METHOD,
        "method_version": PLACEMENT_VERSION,
    }


def _reference_placement_match(user: dict, record: dict, structure: dict, saved: dict,
                              model_id: str) -> dict | None:
    """Load only a completed placement Evaluation matching both exact structures."""
    from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer.dataclasses import validate_rigid_transform_matrix
    from mn_molecular_design_studio.complex_analysis_placement import PLACEMENT_METHOD, PLACEMENT_VERSION

    if saved.get("status") != "validated":
        return None
    source = {
        "molecule_id": str(record.get("id") or ""),
        "structure_id": str(structure.get("id") or ""),
        "group_id": str(structure.get("group_id") or ""),
        "structure_digest": str(saved.get("source_structure_digest") or ""),
        "model_id": str(model_id),
    }
    candidate_id = "database-" + hashlib.sha256(
        (source["molecule_id"] + ":" + source["structure_id"]).encode()
    ).hexdigest()[:24]
    parameters = _reference_placement_parameters(saved, model_id)
    try:
        matches = _matching_database_analysis_jobs(
            user, workflow_id="analysis.reference_placement", source=source,
            candidate_id=candidate_id, model_id=str(model_id), binder_chains=[], target_chains=[],
            parameters=parameters,
        )
    except (OSError, PermissionError, RuntimeError, ValueError, sqlite3.Error) as exc:
        st.info(f"Saved placement job status is unavailable: {exc}")
        return None
    current = next((item for item in matches if item.get("status") in {"completed", "succeeded"}), None)
    if current is None:
        return None
    result_path = Path(current["run_dir"]) / "artifacts" / "analysis" / "reference_placement" / "reference_placement.json"
    result = _read_json(result_path)
    source_result = result.get("source") if isinstance(result.get("source"), dict) else {}
    reference_result = result.get("reference") if isinstance(result.get("reference"), dict) else {}
    if (result.get("schema") != "studio.complex_analysis_reference_placement.v1"
            or result.get("method") != PLACEMENT_METHOD
            or result.get("method_version") != PLACEMENT_VERSION
            or str(result.get("reference_context_id") or "") != str(saved.get("id") or "")
            or str(source_result.get("molecule_id") or "") != source["molecule_id"]
            or str(source_result.get("structure_id") or "") != source["structure_id"]
            or str(source_result.get("group_id") or "") != source["group_id"]
            or str(source_result.get("structure_digest") or "") != source["structure_digest"]
            or str(source_result.get("model_id") or "") != str(model_id)
            or str(source_result.get("chain_id") or "") != str(saved.get("source_component_chain") or "")
            or str(reference_result.get("molecule_id") or "") != str(saved.get("reference_molecule_id") or "")
            or str(reference_result.get("structure_id") or "") != str(saved.get("reference_structure_id") or "")
            or str(reference_result.get("group_id") or "") != str(saved.get("reference_group_id") or "")
            or str(reference_result.get("structure_digest") or "") != str(saved.get("reference_structure_digest") or "")
            or str(reference_result.get("model_id") or "") != str(saved.get("reference_model_id") or "")
            or str(reference_result.get("chain_id") or "") != str(saved.get("reference_component_chain") or "")):
        return None
    try:
        result["matrix_4x4"] = validate_rigid_transform_matrix(result.get("matrix_4x4"))
        rmsd = float(result.get("rmsd_angstrom"))
        pair_count = int(result.get("aligned_pair_count"))
    except (TypeError, ValueError):
        return None
    if not (rmsd >= 0.0 and pair_count >= 3):
        return None
    return {"job": current, "result": result}


def _render_reference_context_track(track: dict, viewer_key: str) -> None:
    residues = list(track.get("residues") or [])
    if not residues:
        st.info(track.get("reason") or "The selected reference does not expose a full polymer sequence track.")
        return
    st.markdown("##### Full reference-chain sequence")
    st.caption(
        "Mapped residues are present in the source construct and mapped to this reference. "
        "Resolved reference-only residues are omitted from the source construct. "
        "Unresolved residues occur in the declared polymer sequence but have no coordinates in this reference representation."
    )
    counts = track.get("counts") or {}
    st.caption(
        f"Chain {track.get('chain_id')} · {len(residues)} declared residues · "
        f"{counts.get('mapped', 0)} mapped/present in source · "
        f"{counts.get('omitted', 0)} reference-resolved but omitted · "
        f"{counts.get('unresolved', 0)} unresolved in reference coordinates"
    )
    page_size = 24
    page_count = max(1, (len(residues) + page_size - 1) // page_size)
    state_key = f"{viewer_key}-reference-highlights"
    page_key = f"{viewer_key}-reference-track-page-{track.get('chain_id')}"
    if int(st.session_state.get(page_key, 0) or 0) >= page_count:
        st.session_state[page_key] = 0
    page = st.selectbox(
        "Reference sequence segment", list(range(page_count)),
        format_func=lambda index: f"Residues {index * page_size + 1}–{min((index + 1) * page_size, len(residues))} of {len(residues)}",
        key=page_key,
    )
    selected = residues[page * page_size:(page + 1) * page_size]
    st.caption("Reference sequence: " + "".join(row["amino_acid"] for row in selected))
    symbols = {"mapped": "●", "omitted": "○", "unresolved": "?"}
    labels = {"mapped": "Mapped and present in source construct",
              "omitted": "Resolved in reference; omitted from source construct",
              "unresolved": "Unresolved in selected reference coordinates"}
    for offset in range(0, len(selected), 12):
        batch = selected[offset:offset + 12]
        columns = st.columns(12)
        for column, residue in zip(columns, batch):
            status = str(residue.get("status") or "unresolved")
            number = str(residue.get("author_residue") or f"seq#{residue['sequence_index'] + 1}")
            selector = str(residue.get("selector") or "")
            column.button(
                f"{symbols.get(status, '?')} {residue['amino_acid']}\n{number}",
                key=f"{viewer_key}-reference-track-{track.get('chain_id')}-{residue['sequence_index']}",
                help=f"Reference chain {track.get('chain_id')} · author residue {number} · {labels.get(status, status)}",
                on_click=_select_inspector_evidence if selector else None,
                args=(state_key, [selector]) if selector else None,
                disabled=not selector,
            )


def _render_optional_reference_context(user: dict, record: dict, structure: dict,
                                       components: list[dict], model_id: str, key: str,
                                       view_context: dict | None = None) -> None:
    """Select and persist an optional Database reference for one protein component."""
    protein_chains = sorted({str(item.get("source_chain") or "") for item in components
                             if item.get("type") == "protein" and item.get("source_chain")})
    if not protein_chains:
        return
    from mn_molecular_design_studio.molecular_database import MolecularDatabase
    from mn_molecular_design_studio.complex_analysis_catalog import authorized_structure_catalog

    database = MolecularDatabase()
    view_context = view_context if isinstance(view_context, dict) else {}
    paired_reference = view_context.get("paired_reference") if isinstance(view_context.get("paired_reference"), dict) else None
    placement = view_context.get("placement") if isinstance(view_context.get("placement"), dict) else None
    with st.expander("Optional component reference context", expanded=bool(placement)):
        st.caption(
            "A component-specific link to another saved Database structure is optional; Studio does not assert that the reference is full-length. Core analyses and compatible comparisons do not require it. "
            "Only exact Database IDs recorded in this structure's provenance are reused automatically; Studio does not infer a reference from sequence similarity or an unrelated cofold."
        )
        component_key = f"complex-analysis-reference-component-{key}"
        if st.session_state.get(component_key) not in protein_chains:
            st.session_state[component_key] = protein_chains[0]
        source_chain = st.selectbox("Source protein component", protein_chains, key=component_key)
        try:
            saved = database.get_complex_reference_context(
                user, source_molecule_id=str(record["id"]),
                source_structure_id=str(structure["id"]), source_chain=source_chain,
            )
        except (PermissionError, ValueError, OSError) as exc:
            saved = None
            st.info(f"Reference context is not configured for this component: {exc}")

        if saved and saved.get("status") in {"validated", "partial", "ambiguous", "unmatched"}:
            mapping = saved.get("mapping") if isinstance(saved.get("mapping"), dict) else {}
            st.markdown(
                f"Saved reference: **{saved.get('reference_record_name') or saved.get('reference_molecule_id')}** · "
                f"{saved.get('reference_structure_label') or saved.get('reference_structure_id')} · chain {saved.get('reference_component_chain')}"
            )
            st.caption(
                f"Source {saved.get('source_molecule_id')} / {saved.get('source_structure_id')} · "
                f"SHA256 {saved.get('source_structure_digest', '')[:16]} · chain {source_chain} · model {saved.get('source_model_id') or model_id}"
            )
            st.caption(
                f"Reference {saved.get('reference_molecule_id')} / {saved.get('reference_structure_id')} · "
                f"SHA256 {saved.get('reference_structure_digest', '')[:16]} · model {saved.get('reference_model_id')} · "
                f"mapping: {saved.get('status')} · {mapping.get('mapping_reason') or 'Residue mapping unavailable.'}"
            )
            st.caption(f"Association source: {mapping.get('association_source') or 'user_selected'}")
            if saved.get("status") == "validated":
                st.info(
                    "Residue correspondence is uniquely mapped and sequence-validated. Optional reference-frame placement is available under Available analyses; "
                    "the reference does not change or gate this complex's existing analyses."
                )
            else:
                st.info(
                    "Reference context is retained for provenance, but reference-dependent evaluations are disabled because residue correspondence is not fully validated. "
                    "Core analyses remain available."
                )
            if (placement and paired_reference
                    and str(saved.get("id") or "") == str(paired_reference.get("saved", {}).get("id") or "")
                    and placement.get("result")):
                from mn_molecular_design_studio.complex_analysis_reference_context import (
                    explicit_reference_components, reference_component_track,
                )
                result = placement["result"]
                st.markdown("#### Reference context in the placed Mol* view")
                st.caption(
                    "The associated reference file is shown in its original coordinates in grey; the selected source complex is transformed into that frame and shown in chain colors. "
                    "The reference file is loaded as saved, including its other resolved chains and components. Its coverage of the biological protein is not inferred."
                )
                st.caption(
                    f"Reference Database record {saved.get('reference_molecule_id')} · structure {saved.get('reference_structure_id')} · "
                    f"model {saved.get('reference_model_id')} · chain {saved.get('reference_component_chain')} · "
                    f"SHA256 {str(saved.get('reference_structure_digest') or '')[:16]}"
                )
                reference_track = reference_component_track(
                    paired_reference["path"], str(saved.get("reference_component_chain") or ""),
                    result.get("mapped_residue_pairs") or [],
                )
                if reference_track.get("status") == "available":
                    _render_reference_context_track(
                        reference_track, f"complex-analysis-db-viewer-{structure['id']}",
                    )
                    st.caption("Residue buttons select the exact author-numbered residue in the reference (first structure) in the paired Mol* viewer. They do not add or infer measurements.")
                else:
                    st.info(reference_track.get("reason") or "Full reference sequence context is unavailable for this file.")
                try:
                    explicit_components = explicit_reference_components(paired_reference["path"])
                except (OSError, ValueError, TypeError):
                    explicit_components = []
                if explicit_components:
                    st.markdown("##### Explicit branched entities declared in the reference file")
                    st.dataframe(explicit_components, hide_index=True, use_container_width=True)
                    st.caption("These are coordinate-file annotations only. A declared glycan entity does not establish occupancy or accessibility; no context measurement is calculated.")
                else:
                    if paired_reference["path"].name.casefold().endswith((".pdb", ".ent")):
                        st.caption("This PDB file has no mmCIF branched-entity annotations; non-polymer coordinates remain visible in Mol*. Studio does not infer missing glycans.")
                    else:
                        st.caption("No explicitly declared branched entities were detected in this reference file; Studio does not infer missing glycans.")
                if paired_reference.get("components"):
                    st.markdown("##### Other parsed reference components")
                    st.dataframe([{
                        "Type": component.get("type"),
                        "Author chain": component.get("source_chain"),
                        "Component": component.get("label"),
                        "Residue / identity": component.get("residue_key") or component.get("smiles") or component.get("sequence"),
                    } for component in paired_reference["components"]], hide_index=True, use_container_width=True)
            residues = mapping.get("residue_mapping") if isinstance(mapping.get("residue_mapping"), list) else []
            if residues:
                rows = []
                for item in residues:
                    source_residue, reference_residue = item.get("source") or {}, item.get("reference") or {}
                    rows.append({
                        "Source chain": source_chain,
                        "Source residue": f"{source_residue.get('number')}{source_residue.get('insertion') or ''}",
                        "Source sequence index": source_residue.get("sequence_index"),
                        "Source amino acid": source_residue.get("amino_acid"),
                        "Reference chain": saved.get("reference_component_chain"),
                        "Reference residue": f"{reference_residue.get('number')}{reference_residue.get('insertion') or ''}",
                        "Reference sequence index": reference_residue.get("sequence_index"),
                        "Reference amino acid": reference_residue.get("amino_acid"),
                        "Exact sequence match": bool(item.get("sequence_match")),
                    })
                st.markdown("##### Residue-numbering provenance")
                st.caption("Author residue numbers and insertion codes come from the exact saved PDB/mmCIF artifacts; sequence indices and amino-acid identities are retained too.")
                st.dataframe(rows, hide_index=True, use_container_width=True)
                st.caption(
                    f"{mapping.get('mapping_method')} · alignments {mapping.get('alignment_count')} · "
                    f"source coverage {float(mapping.get('source_coverage') or 0):.3f} · "
                    f"mapped sequence identity {float(mapping.get('mapped_sequence_identity') or 0):.3f}"
                )
                st.json(mapping.get("mapping_parameters") or {})
            st.caption("The exact associated Database reference is displayed beside the selected complex above in its original coordinate frame.")
        elif saved:
            status_message = {
                "stale": "The saved source or reference file changed; reselect the reference to refresh its digest-bound numbering map.",
                "source_unavailable": "The source representation is no longer accessible.",
                "reference_unavailable": "The reference is no longer accessible.",
            }.get(str(saved.get("status")), "The saved reference context is unavailable.")
            st.warning(status_message + " The complex remains available for its existing analyses.")
        else:
            st.info("Reference-context analyses: not configured. This is not a missing-input condition for the complex.")

        linked_molecule, linked_structure, link_provenance = _explicit_reference_link(structure)
        search = st.text_input("Find a Database reference", key=f"complex-analysis-reference-search-{key}", placeholder="Name or sequence")
        try:
            reference_records = []
            if search.strip():
                catalog_key = f"complex-analysis-reference-catalog-{user.get('id', '')}-{key}-{source_chain}"
                cached_catalog = st.session_state.get(catalog_key)
                same_query = isinstance(cached_catalog, dict) and str(cached_catalog.get("query") or "") == search
                refresh_catalog = st.button("Refresh reference choices", key=f"complex-analysis-reference-refresh-{key}-{source_chain}")
                load_next = st.button("Load next reference page", key=f"complex-analysis-reference-next-{key}-{source_chain}",
                                      disabled=not (same_query and cached_catalog.get("next_offset") is not None))
                catalog = authorized_structure_catalog(
                    database, user, session_state=st.session_state, cache_key=catalog_key,
                    query=search, limit=100, advance=load_next, refresh=refresh_catalog,
                )
                reference_records = catalog["records"]
        except (PermissionError, ValueError, OSError) as exc:
            reference_records = []
            st.warning(f"Could not list reference structures: {exc}")
        choices = {str(row["id"]): row for row in reference_records
                   if row.get("representations") or row.get("structures")}
        saved_molecule = str((saved or {}).get("reference_molecule_id") or "")
        if linked_structure and not linked_molecule:
            try:
                linked_molecule = database.accessible_structure_molecule_id(user, linked_structure)
            except (FileNotFoundError, PermissionError, ValueError, OSError):
                linked_molecule = ""
        for candidate_id in (saved_molecule, linked_molecule):
            if candidate_id and candidate_id not in choices:
                try:
                    extra = database.list_molecules(user, molecule_ids=[candidate_id], collection="all", limit=1)["items"]
                    if extra and (extra[0].get("representations") or extra[0].get("structures")):
                        choices[candidate_id] = extra[0]
                except (PermissionError, ValueError, OSError):
                    pass
        molecule_options = ["", *choices]
        default_molecule = saved_molecule if saved_molecule in choices else linked_molecule if linked_molecule in choices else ""
        molecule_key = f"complex-analysis-reference-molecule-{key}-{source_chain}"
        if st.session_state.get(molecule_key) not in molecule_options:
            st.session_state[molecule_key] = default_molecule
        reference_molecule_id = st.selectbox(
            "Reference Database record", molecule_options,
            format_func=lambda value: "Choose a saved reference…" if not value else f"{choices[value].get('name') or 'Structure'} · {value[-8:]}",
            key=molecule_key,
        )
        if reference_molecule_id:
            reference_record = choices[reference_molecule_id]
            representations = list(reference_record.get("representations") or reference_record.get("structures") or [])
            structure_choices = {str(item["id"]): item for item in representations
                                 if str(item.get("format") or "").casefold() in {"pdb", "ent", "cif", "mmcif"}}
            structure_options = ["", *structure_choices]
            saved_structure = str((saved or {}).get("reference_structure_id") or "")
            default_structure = saved_structure if saved_structure in structure_choices else linked_structure if linked_structure in structure_choices else ""
            structure_key = f"complex-analysis-reference-structure-{key}-{source_chain}"
            if st.session_state.get(structure_key) not in structure_options:
                st.session_state[structure_key] = default_structure
            reference_structure_id = st.selectbox(
                "Reference structure / model", structure_options,
                format_func=lambda value: "Choose a saved structure…" if not value else _representation_label(structure_choices[value]),
                key=structure_key,
            )
            reference_chain = ""
            if reference_structure_id:
                try:
                    reference_path = database.structure_file(user, reference_structure_id)
                    reference_components = component_inventory(reference_path.read_bytes(), reference_path.name, reference_record.get("components"))
                    reference_chains = sorted({str(item.get("source_chain") or "") for item in reference_components
                                               if item.get("type") == "protein" and item.get("source_chain")})
                except (OSError, PermissionError, ValueError) as exc:
                    reference_path, reference_chains = None, []
                    st.warning(f"Could not inspect the selected reference structure: {exc}")
                if reference_path and not reference_chains:
                    st.info("This saved representation contains no parsed protein chain for component mapping.")
                if reference_chains:
                    saved_chain = str((saved or {}).get("reference_component_chain") or "")
                    source_context = structure.get("context") if isinstance(structure.get("context"), dict) else {}
                    linked_chains = [str(value) for value in source_context.get("target_chains") or []]
                    default_chain = saved_chain if saved_chain in reference_chains else (
                        linked_chains[0] if linked_structure == reference_structure_id and len(linked_chains) == 1
                        and linked_chains[0] in reference_chains else ""
                    )
                    chain_options = ["", *reference_chains]
                    chain_key = f"complex-analysis-reference-chain-{key}-{source_chain}"
                    if st.session_state.get(chain_key) not in chain_options:
                        st.session_state[chain_key] = default_chain
                    reference_chain = st.selectbox(
                        "Reference protein chain", chain_options,
                        format_func=lambda value: "Choose a reference component…" if not value else value,
                        key=chain_key,
                    )
                if st.button("Save component reference", key=f"complex-analysis-reference-save-{key}-{source_chain}",
                             type="primary", disabled=not reference_chain):
                    try:
                        database.save_complex_reference_context(
                            user, source_molecule_id=str(record["id"]), source_structure_id=str(structure["id"]),
                            source_chain=source_chain, reference_molecule_id=reference_molecule_id,
                            reference_structure_id=reference_structure_id, reference_chain=reference_chain,
                            association_source=(link_provenance if linked_structure == reference_structure_id
                                                and linked_molecule in {"", reference_molecule_id} else "user_selected"),
                        )
                        st.success("Optional component reference saved with exact structure digests and residue-numbering provenance.")
                        st.rerun()
                    except (OSError, PermissionError, TypeError, ValueError) as exc:
                        st.error(f"Could not save this component reference: {exc}")
        if saved and saved.get("configured"):
            if st.button("Clear component reference", key=f"complex-analysis-reference-clear-{key}-{source_chain}"):
                try:
                    database.clear_complex_reference_context(
                        user, source_molecule_id=str(record["id"]),
                        source_structure_id=str(structure["id"]), source_chain=source_chain,
                    )
                    st.rerun()
                except (OSError, PermissionError, ValueError) as exc:
                    st.error(f"Could not clear this component reference: {exc}")


def _database_chain_defaults(context: dict, protein_chains: list[str], *,
                             stored_components: list[dict] | None = None,
                             parsed_components: list[dict] | None = None) -> tuple[list[str], list[str]]:
    """Use roles recorded for this representation, validated against its chains."""
    binder = context.get("binder_chains")
    target = context.get("target_chains")
    if (isinstance(binder, list) and isinstance(target, list) and binder
            and all(isinstance(chain, str) for chain in binder + target)
            and len(binder) == len(set(binder)) and len(target) == len(set(target))
            and set(binder + target) <= set(protein_chains)
            and not set(binder) & set(target)
            and (target or len(protein_chains) == len(binder))):
        return list(binder), list(target)
    if binder is None and target is None:
        # Older binder identities may predate representation-level roles. Use
        # their explicit role only if chain AND sequence match this actual file.
        parsed = {item.get("source_chain") or item.get("chain"): item.get("sequence")
                  for item in parsed_components or [] if item.get("type") == "protein"}
        matched = [item for item in stored_components or []
                   if item.get("type") == "protein" and item.get("sequence")
                   and parsed.get(item.get("chain")) == item["sequence"]]
        binders = sorted({item["chain"] for item in matched if item.get("role") == "binder"})
        targets = sorted({item["chain"] for item in matched if item.get("role") == "target"})
        if binders:
            return binders, targets or [chain for chain in protein_chains if chain not in binders]
    return protein_chains[:1], protein_chains[1:]


def _render_database_analysis(user: dict, record: dict, structure: dict, path: Path,
                              components: list[dict], evaluation_id: str, *,
                              workflow_schema=None, submit=None) -> None:
    protein_chains = sorted({str(item.get("source_chain") or "") for item in components
                             if item.get("type") == "protein" and item.get("source_chain")})
    kinds = {str(item.get("type") or "other") for item in components}
    group_id = str(structure.get("group_id") or "")
    from mn_molecular_design_studio.molecular_database import MolecularDatabase
    db_group = MolecularDatabase().group(group_id, user)
    project_id = str(db_group.get("project_id") or "")
    context = structure.get("context") if isinstance(structure.get("context"), dict) else {}
    default_binder, default_target = _database_chain_defaults(
        context, protein_chains, stored_components=record.get("components"), parsed_components=components,
    )
    from mn_molecular_design_studio.scientific_settings import model_identity
    model_id = model_identity(context.get("model_id"), context.get("model_index"), context.get("model"), structure.get("model_index"), structure["id"])
    source = {
        "molecule_id": str(record["id"]), "structure_id": str(structure["id"]),
        "group_id": group_id, "evaluation_id": evaluation_id, "model_id": model_id,
        "binder_chains": [], "target_chains": [],
        "structure_context": context, "group_origin": db_group.get("origin") or {},
    }
    digest = _structure_digest(path)
    if not digest:
        st.info("The selected structure cannot be checksummed, so contextual analysis is unavailable.")
        return
    source["structure_digest"] = digest
    from mn_molecular_design_studio.complex_analysis_sources import database_pae_evidence
    pae_evidence = database_pae_evidence(user, source)
    source["pae_digest"] = pae_evidence["sha256"] if pae_evidence else ""
    key = hashlib.sha256(f"database:{record['id']}:{structure['id']}".encode()).hexdigest()[:14]
    st.markdown("### Available analyses")
    st.caption("Choose a registered analysis explicitly. Opening a Database structure does not start calculations.")
    schemas, schema_errors = {}, {}
    for workflow_id in ("analysis.biotite_interface", "analysis.rosetta", "analysis.ipsae", "analysis.pep_patch",
                        "analysis.plip", "analysis.pandamap", "analysis.posebusters",
                        "analysis.reference_environment", "analysis.tag_fusion_suitability", "analysis.tag_fusion_context",
                        "analysis.sequence_liability", "analysis.binder_interface_chemistry"):
        try:
            schemas[workflow_id] = workflow_schema(workflow_id) if workflow_schema else None
            if not schemas[workflow_id]:
                schema_errors[workflow_id] = "The workflow registry did not return this analysis schema."
        except Exception as exc:
            schema_errors[workflow_id] = f"The registered workflow is unavailable: {exc}"

    ppi = kinds == {"protein"} and len(protein_chains) >= 2
    if ppi:
        st.markdown("#### Protein–protein")
        st.caption(f"Selected representation: {structure['id']} · model {model_id} · SHA256 {digest[:16]}")
        left, right = st.columns(2)
        group_a_key, group_b_key = f"complex-analysis-db-a-{key}", f"complex-analysis-db-b-{key}"
        if group_a_key not in st.session_state:
            st.session_state[group_a_key] = default_binder
        if group_b_key not in st.session_state:
            st.session_state[group_b_key] = default_target
        binder_chains = left.multiselect("Interface chain group A", protein_chains, key=group_a_key)
        target_chains = right.multiselect("Interface chain group B", protein_chains, key=group_b_key)
        source["binder_chains"], source["target_chains"] = list(binder_chains), list(target_chains)
        invalid_groups = not binder_chains or not target_chains or bool(set(binder_chains) & set(target_chains))
        if invalid_groups:
            st.info("Choose nonempty, disjoint protein chain groups. Multi-chain groups are supported by Biotite Interface.")
        interface_params = {"contact_distance": float(next((item.get("default", 5.0) for item in (schemas.get("analysis.biotite_interface") or {}).get("parameters", []) if item.get("name") == "contact_distance"), 5.0))}
        interface_params["contact_distance"] = round(float(st.number_input(
            "Biotite contact distance (Å)", min_value=0.1, max_value=20.0,
            value=interface_params["contact_distance"], step=0.1, key=f"complex-analysis-db-biotite-cutoff-{key}")), 3)
        def show_ppi(workflow_id: str, parameters: dict, unavailable: str = "") -> None:
            _render_database_workflow(user, workflow_id, schemas.get(workflow_id), source, model_id,
                                      binder_chains, target_chains, parameters, key, submit,
                                      project_id, unavailable or ("Choose two nonempty, disjoint protein chain groups." if invalid_groups else ""),
                                      structure_path=path, structure_context=context, structure_components=components,
                                      viewer_key=f"complex-analysis-db-viewer-{structure['id']}")
        show_ppi("analysis.biotite_interface", interface_params, schema_errors.get("analysis.biotite_interface", ""))
        st.caption("This fixed-method profile uses 960-point Shrake–Rupley surface sampling and documented geometric polar cutoffs; its method and assumptions are saved with the Evaluation.")
        show_ppi("analysis.binder_interface_chemistry", {}, schema_errors.get("analysis.binder_interface_chemistry", ""))
        show_ppi("analysis.rosetta", {}, schema_errors.get("analysis.rosetta", "") or (
            "The registered Rosetta workflow currently accepts exactly one chain in each group."
            if len(binder_chains) != 1 or len(target_chains) != 1 else ""))
        ipsae_schema = schemas.get("analysis.ipsae")
        ipsae_schema_params = {str(item.get("name")): item for item in (ipsae_schema or {}).get("parameters", [])}
        defaults = {name: float(ipsae_schema_params.get(name, {}).get("default", default)) for name, default in (
            ("pae_cutoff", 10.), ("dist_cutoff", 10.), ("ipdae_contact_cutoff", 8.),
            ("structure_contact_cutoff", 4.), ("hotspot_contact_cutoff", 5.))}
        with st.expander("Configure ipSAE parameters", expanded=False):
            fields = st.columns(3)
            ipsae_params = {name: round(float(fields[index % len(fields)].number_input(
                name.replace("_", " ").title() + " (Å)", min_value=0.1, max_value=50.0,
                value=value, step=0.5, key=f"complex-analysis-db-ipsae-{name}-{key}")), 3)
                for index, (name, value) in enumerate(defaults.items())}
        st.caption("ipSAE can use saved structure geometry here. PAE-derived terms remain unavailable unless the source stores a compatible pairwise PAE matrix.")
        show_ppi("analysis.ipsae", ipsae_params, schema_errors.get("analysis.ipsae", "") or (
            "The registered ipSAE workflow currently accepts exactly one chain in each group."
            if len(binder_chains) != 1 or len(target_chains) != 1 else ""))
    elif {"protein", "small_molecule"} <= kinds:
        st.markdown("#### Protein–ligand")
        st.caption("PLIP and PandaMap report ligand–receptor interactions. PoseBusters checks pose geometry and chemistry; it does not report interaction strength.")
        ligand_components = [item for item in components if str(item.get("type") or "").casefold() == "small_molecule"]
        ligand_key = f"complex-analysis-db-ligand-{key}"
        ligand_options = [item.get("id") for item in ligand_components]
        ligand_labels = {item.get("id"): f"{item.get('label') or item.get('ccd') or 'Ligand'} · {item.get('residue_key') or item.get('source_chain')}"
                         for item in ligand_components}
        if st.session_state.get(ligand_key) not in ligand_options:
            st.session_state[ligand_key] = ligand_options[0] if ligand_options else ""
        selected_ligand_id = st.selectbox("Ligand residue", ligand_options,
            format_func=lambda value: ligand_labels[value], key=ligand_key)
        selected_ligand = next((item for item in ligand_components if item.get("id") == selected_ligand_id), {})
        ligand_selection = _ligand_selection_for_component(selected_ligand, ligand_components, {
            "components": record.get("components") or [],
            "ligand_selection": context.get("ligand_selection"),
        })
        st.caption(f"Selected representation: {structure['id']} · model {model_id} · SHA256 {digest[:16]}")
        unavailable = ""
        if not ligand_selection.get("chain") or not ligand_selection.get("residue_key"):
            unavailable = "The selected ligand has no stable source chain and residue identity."
        elif not ligand_selection.get("smiles"):
            unavailable = "The saved record does not preserve this ligand's isomeric SMILES identity. A verified source identity is required for topology matching."
        for workflow_id in ("analysis.plip", "analysis.pandamap"):
            _render_database_workflow(user, workflow_id, schemas.get(workflow_id), source, model_id, [], [], {},
                key + "-" + workflow_id.rsplit(".", 1)[-1], submit, project_id,
                schema_errors.get(workflow_id, "") or unavailable, structure_path=path,
                structure_context=context, structure_components=components, viewer_key=f"complex-analysis-db-viewer-{structure['id']}",
                ligand_selection=ligand_selection)
        _render_database_workflow(user, "analysis.posebusters", schemas.get("analysis.posebusters"), source, model_id, [], [], {},
            key + "-posebusters", submit, project_id,
            schema_errors.get("analysis.posebusters", "") or unavailable, structure_path=path,
            structure_context=context, structure_components=components, viewer_key=f"complex-analysis-db-viewer-{structure['id']}",
            ligand_selection=ligand_selection)
    elif {"protein", "dna"} <= kinds or {"protein", "rna"} <= kinds:
        st.markdown("#### Protein–DNA/RNA")
        st.info("The selected representation contains protein and nucleic acid. The registered protein–protein interface workflows do not apply; no contextual DNA/RNA interaction evaluator is currently registered here.")
    else:
        st.info("No compatible registered interaction analysis was identified for this component combination.")

    if protein_chains:
        st.markdown("#### Sequence and Developability Liabilities")
        st.caption("Post-hoc sequence and model screen. Exact saved Biotite, PEP-Patch, PLIP, PandaMap and Sequence Composition Evaluations are reused when available; this workflow does not rerun those scientific tools.")
        liability_component_key = f"complex-analysis-liability-components-{key}"
        if st.session_state.get(liability_component_key) is None:
            st.session_state[liability_component_key] = default_binder
        selected_liability_chains = st.multiselect(
            "Protein chain(s) to assess", protein_chains,
            key=liability_component_key,
        )
        liability_partner_options = [chain for chain in protein_chains if chain not in selected_liability_chains]
        liability_partner_key = f"complex-analysis-liability-partners-{key}"
        if liability_partner_key not in st.session_state:
            st.session_state[liability_partner_key] = liability_partner_options
        else:
            st.session_state[liability_partner_key] = [value for value in st.session_state[liability_partner_key]
                                                       if value in liability_partner_options]
        liability_partner_chains = st.multiselect(
            "Optional protein interface partners", liability_partner_options,
            key=liability_partner_key,
            help="Interface-specific findings are included only when matching saved contact/interactions Evaluations exist.",
        )
        liability_ligand_selection = {}
        liability_ligands = [item for item in components if str(item.get("type") or "").casefold() == "small_molecule"]
        if liability_ligands:
            ligand_evidence_key = f"complex-analysis-liability-ligand-evidence-{key}"
            ligand_options = [""] + [str(item.get("id") or "") for item in liability_ligands]
            ligand_labels = {str(item.get("id") or ""): f"{item.get('label') or item.get('ccd') or 'Ligand'} · {item.get('residue_key') or item.get('source_chain')}"
                             for item in liability_ligands}
            chosen_ligand_id = st.selectbox(
                "Optional ligand interaction-evidence selection", ligand_options,
                format_func=lambda value: "No ligand-specific evidence" if not value else ligand_labels.get(value, value),
                key=ligand_evidence_key,
                help="This only links a matching saved PLIP/PandaMap Evaluation to protein interface findings; it does not rerun or score the ligand pose.",
            )
            if chosen_ligand_id:
                chosen_ligand = next(item for item in liability_ligands if str(item.get("id") or "") == chosen_ligand_id)
                liability_ligand_selection = _ligand_selection_for_component(chosen_ligand, liability_ligands, {
                    "components": record.get("components") or [],
                    "ligand_selection": context.get("ligand_selection"),
                })
                if not liability_ligand_selection.get("smiles"):
                    st.warning("The selected ligand has no saved isomeric SMILES, so it cannot identify an exact saved interaction Evaluation.")
                    liability_ligand_selection = {}
                else:
                    st.caption("Only saved PLIP/PandaMap interaction artifacts with this exact ligand identity are eligible for interface-specific findings.")
        format_key = f"complex-analysis-liability-format-{key}"
        liability_format_labels = {
            "linear_protein": "Linear protein chain(s)",
            "cyclic_construct": "Cyclic construct (no free termini)",
            "multichain_construct": "Multichain construct (chain-by-chain screen)",
        }
        if st.session_state.get(format_key) not in liability_format_labels:
            st.session_state[format_key] = "multichain_construct" if len(selected_liability_chains) > 1 else "linear_protein"
        component_format = st.selectbox("Construct format", list(liability_format_labels),
                                        format_func=lambda value: liability_format_labels[value], key=format_key)
        b_factor_declared = st.checkbox(
            "The selected structure stores predictor pLDDT values in Cα B-factors (verified from its source)",
            value=False, key=f"complex-analysis-liability-plddt-{key}",
            help="Otherwise B-factors remain uninterpreted metadata and do not affect a liability severity.",
        )
        liability_source = {**source, "binder_chains": list(selected_liability_chains),
                            "target_chains": list(liability_partner_chains),
                            "ligand_selection": liability_ligand_selection}
        _render_database_workflow(
            user, "analysis.sequence_liability", schemas.get("analysis.sequence_liability"),
            liability_source, model_id, list(selected_liability_chains), list(liability_partner_chains),
            {"component_format": component_format,
             "b_factor_interpretation": "plddt_0_100" if b_factor_declared else "uninterpreted"},
            key + "-sequence-liability", submit, project_id,
            schema_errors.get("analysis.sequence_liability", "") or
            ("Choose one or more protein chains to assess." if not selected_liability_chains else ""),
            structure_path=path, structure_context=context, structure_components=components,
            viewer_key=f"complex-analysis-db-viewer-{structure['id']}",
            ligand_selection=liability_ligand_selection or None,
        )

    if len(protein_chains) >= 2:
        st.markdown("#### Tag and Fusion Suitability")
        st.caption("Post-hoc advisory for one selected protein component against a selected protein interface. Flexible tags use geometric clearance paths; folded fusions use a spherical clearance envelope. Neither model simulates the actual tag or predicts assay performance.")
        component_key = f"complex-analysis-tag-fusion-component-{key}"
        if st.session_state.get(component_key) not in protein_chains:
            st.session_state[component_key] = default_binder[0]
        assessed_chain = st.selectbox("Protein component to assess", protein_chains,
                                      key=component_key)
        partner_options = [chain for chain in protein_chains if chain != assessed_chain]
        partner_key = f"complex-analysis-tag-fusion-partners-{key}-{assessed_chain}"
        if partner_key not in st.session_state:
            st.session_state[partner_key] = partner_options
        partner_chains = st.multiselect("Interface partner chains", partner_options,
                                        key=partner_key)
        format_key = f"complex-analysis-tag-fusion-format-{key}"
        format_labels = {
            "single_chain_linear_protein": "Single-chain linear protein",
            "cyclic_peptide": "Cyclic peptide (not applicable)",
            "multichain_construct": "Multichain construct (not applicable)",
        }
        component_format = st.selectbox("Component format", list(format_labels),
            format_func=lambda value: format_labels[value], key=format_key)
        termini_confirmed = st.checkbox(
            "I confirm the selected representation contains the complete intended construct termini",
            value=False, key=f"complex-analysis-tag-fusion-termini-{key}",
            help="If the coordinate file declares an unresolved terminal sequence position, the evaluator still marks that end incomplete.",
        )
        plddt_declared = st.checkbox(
            "The selected structure stores predictor pLDDT in Cα B-factors (verified from its source)",
            value=False, key=f"complex-analysis-tag-fusion-plddt-{key}",
            help="B-factors are otherwise shown only as uninterpreted file metadata.",
        )
        tag_cutoff = round(float(st.number_input(
            "Interface contact cutoff (Å)", min_value=1.0, max_value=10.0, value=4.5, step=0.25,
            key=f"complex-analysis-tag-fusion-cutoff-{key}")), 3)
        tag_parameters = {
            "interface_contact_cutoff": tag_cutoff,
            "termini_confirmed": bool(termini_confirmed),
            "component_format": component_format,
            "b_factor_interpretation": "plddt_0_100" if plddt_declared else "uninterpreted",
        }
        tag_source = {**source, "binder_chains": [assessed_chain], "target_chains": list(partner_chains)}
        tag_unavailable = "" if partner_chains else "Choose at least one disjoint protein interface partner chain."
        _render_database_workflow(
            user, "analysis.tag_fusion_suitability", schemas.get("analysis.tag_fusion_suitability"),
            tag_source, model_id, [assessed_chain], list(partner_chains), tag_parameters,
            key + "-tag-fusion", submit, project_id,
            schema_errors.get("analysis.tag_fusion_suitability", "") or tag_unavailable,
        )

        st.markdown("##### Optional full-reference surroundings")
        context_unavailable = ""
        contextual_parameters = dict(tag_parameters)
        db = MolecularDatabase()
        reference_associations = {}
        for chain in protein_chains:
            try:
                saved = db.get_complex_reference_context(
                    user, source_molecule_id=str(record["id"]),
                    source_structure_id=str(structure["id"]), source_chain=chain,
                )
            except (OSError, PermissionError, ValueError, TypeError, KeyError):
                continue
            if isinstance(saved, dict) and saved.get("id") and saved.get("status") in {"validated", "stale"}:
                reference_associations[chain] = saved
        if reference_associations:
            reference_key = f"complex-analysis-tag-fusion-reference-chain-{key}-{assessed_chain}"
            reference_chains = sorted(reference_associations)
            if st.session_state.get(reference_key) not in reference_chains:
                st.session_state[reference_key] = reference_chains[0]
            reference_source_chain = st.selectbox(
                "Reference-associated source chain used to place the complex",
                reference_chains, key=reference_key,
                help="The explicit reference association may belong to another component in this same saved complex. Its saved transform places the entire complex.",
            )
            contextual_source = {**tag_source, "source_component_chain": reference_source_chain}
            associated_reference = reference_associations[reference_source_chain]
            if associated_reference.get("status") != "validated":
                context_unavailable = "This reference association is stale or ambiguous. The intrinsic assessment above remains available."
            else:
                matched_placement = _reference_placement_match(
                    user, record, structure, associated_reference, model_id,
                )
                if matched_placement is None:
                    context_unavailable = "Run and save Reference-frame placement for this exact association first. Contextual clearance does not refit or infer a reference."
                else:
                    st.caption(f"Uses the saved placement Evaluation for source chain {reference_source_chain} and includes resolved reference coordinates outside its mapped component.")
                    contextual_parameters["source_component_chain"] = reference_source_chain
        else:
            contextual_source = {**tag_source, "source_component_chain": ""}
            context_unavailable = "No current, explicit reference association is saved for a protein component in this structure. The intrinsic assessment above remains available."
        _render_database_workflow(
            user, "analysis.tag_fusion_context", schemas.get("analysis.tag_fusion_context"),
            contextual_source, model_id, [assessed_chain], list(partner_chains),
            contextual_parameters, key + "-tag-fusion-context", submit, project_id,
            schema_errors.get("analysis.tag_fusion_context", "") or tag_unavailable or context_unavailable,
        )

    if protein_chains:
        placement_workflow = "analysis.reference_placement"
        try:
            schemas[placement_workflow] = workflow_schema(placement_workflow) if workflow_schema else None
            if not schemas[placement_workflow]:
                schema_errors[placement_workflow] = "The workflow registry did not return this analysis schema."
        except Exception as exc:
            schemas[placement_workflow] = None
            schema_errors[placement_workflow] = f"The registered workflow is unavailable: {exc}"
        reference_component_key = f"complex-analysis-reference-component-{key}"
        source_component_chain = str(st.session_state.get(reference_component_key) or "")
        if source_component_chain not in protein_chains:
            source_component_chain = protein_chains[0]
        placement_source = {**source, "binder_chains": [], "target_chains": []}
        placement_parameters: dict = {}
        placement_unavailable = ""
        try:
            saved_reference = MolecularDatabase().get_complex_reference_context(
                user, source_molecule_id=str(record["id"]),
                source_structure_id=str(structure["id"]), source_chain=source_component_chain,
            )
        except (OSError, PermissionError, ValueError) as exc:
            saved_reference = None
            placement_unavailable = f"The saved reference association is unavailable: {exc}"
        if not placement_unavailable:
            if not saved_reference:
                placement_unavailable = "Select and save an explicit component reference in the Inspector's Reference Context panel first."
            elif saved_reference.get("status") != "validated":
                placement_unavailable = (
                    f"Reference placement is unavailable for this association ({saved_reference.get('status')}). "
                    "It requires a current, unique, exact sequence correspondence; core analyses remain available."
                )
            else:
                placement_parameters = _reference_placement_parameters(saved_reference, model_id)
                if not all(placement_parameters.get(name) for name in (
                    "source_component_chain", "reference_context_id", "reference_molecule_id",
                    "reference_structure_id", "reference_structure_digest", "reference_model_id",
                    "reference_component_chain",
                )):
                    placement_unavailable = "The saved reference association is missing required identity provenance."
        st.markdown("#### Optional reference-frame placement")
        st.caption("Fits the selected protein component to its explicitly associated reference chain with mapped Cα atoms. It records a rigid display transform; it does not infer reference identity or calculate biological accessibility.")
        _render_database_workflow(
            user, placement_workflow, schemas.get(placement_workflow), placement_source,
            model_id, [], [], placement_parameters, key + "-reference-placement", submit,
            project_id, schema_errors.get(placement_workflow, "") or placement_unavailable,
            structure_path=path, structure_context=context, structure_components=components,
            viewer_key=f"complex-analysis-db-viewer-{structure['id']}",
        )
        st.markdown("#### Optional full-reference steric and environmental context")
        environment_unavailable = ""
        environment_source = {**source, "binder_chains": [], "target_chains": [],
                              "source_component_chain": source_component_chain}
        environment_parameters = {"source_component_chain": source_component_chain,
                                  "selected_reference_chains": [], "assembly_id": None}
        environment_reference_chains = []
        environment_association = saved_reference if isinstance(saved_reference, dict) else None
        if not environment_association or environment_association.get("status") != "validated":
            environment_unavailable = "A current, explicit, validated component reference association is required. Ordinary complex analyses remain available."
        else:
            placement_match = _reference_placement_match(user, record, structure, environment_association, model_id)
            if placement_match is None:
                environment_unavailable = "Run and save Reference-frame placement for this exact association first."
            else:
                try:
                    database = MolecularDatabase()
                    ref_path = database.structure_file(user, str(environment_association["reference_structure_id"]), materialize=True)
                    from mn_molecular_design_studio.complex_analysis_comparison import _read_biopython_structure
                    ref_models = list(_read_biopython_structure(ref_path).get_models())
                    if len(ref_models) != 1:
                        raise ValueError("The associated reference must have the same single-model representation validated for placement.")
                    environment_reference_chains = sorted(str(chain.id) for chain in ref_models[0]
                                                          if str(chain.id) != str(environment_association.get("reference_component_chain") or ""))
                    environment_source.update({
                        "reference_context_id": str(environment_association.get("id") or ""),
                        "reference_molecule_id": str(environment_association.get("reference_molecule_id") or ""),
                        "reference_structure_id": str(environment_association.get("reference_structure_id") or ""),
                        "reference_structure_digest": str(environment_association.get("reference_structure_digest") or ""),
                        "reference_model_id": str(environment_association.get("reference_model_id") or ""),
                        "reference_chain_id": str(environment_association.get("reference_component_chain") or ""),
                        "placement_run_id": str(placement_match["job"].get("run_id") or ""),
                    })
                except (OSError, PermissionError, ValueError, KeyError, TypeError) as exc:
                    environment_unavailable = f"The selected reference coordinates cannot be used for context analysis: {exc}"
        if not environment_unavailable:
            selected_context_key = f"complex-analysis-reference-environment-chains-{key}-{source_component_chain}"
            if st.session_state.get(selected_context_key) is None:
                st.session_state[selected_context_key] = []
            selected_context_chains = st.multiselect(
                "Additional reference components to compare",
                environment_reference_chains, key=selected_context_key,
                help="The validated mapped reference chain is always excluded from context comparison at mapped residues. Other chains are included only when explicitly selected.",
            )
            environment_parameters["selected_reference_chains"] = list(selected_context_chains)
            assembly_id = st.text_input(
                "Reference assembly ID (optional)", value="",
                key=f"complex-analysis-reference-environment-assembly-{key}-{source_component_chain}",
                help="Select an assembly only when its identity-operator membership is explicit in the associated mmCIF. Non-identity assembly operators are not expanded.",
            ).strip()
            environment_parameters["assembly_id"] = assembly_id or None
            st.caption("Reports resolved omitted-reference residues, closest heavy-atom distances, atom-pair counts and potential van-der-Waals envelope overlaps. Static overlaps are advisory. Explicit modeled sugars are separate from sequence-only motifs; no glycan occupancy or accessibility is inferred.")
        _render_database_workflow(
            user, "analysis.reference_environment", schemas.get("analysis.reference_environment"),
            environment_source, model_id, [], [], environment_parameters,
            key + "-reference-environment", submit, project_id,
            schema_errors.get("analysis.reference_environment", "") or environment_unavailable,
            structure_path=path, structure_context=context, structure_components=components,
            viewer_key=f"complex-analysis-db-viewer-{structure['id']}",
        )

    if protein_chains:
        st.markdown("#### Protein surface analysis")
        pep_key = f"complex-analysis-db-pep-chains-{key}"
        if pep_key not in st.session_state:
            st.session_state[pep_key] = default_binder
        surface_chains = st.multiselect("PEP-Patch surface chain group", protein_chains, key=pep_key)
        st.number_input("PEP-Patch pH", min_value=0.0, max_value=14.0, value=7.4, step=0.1,
                        key=f"complex-analysis-db-pep-ph-{key}")
        surface_source = {**source, "binder_chains": list(surface_chains), "target_chains": []}
        _render_database_workflow(user, "analysis.pep_patch", schemas.get("analysis.pep_patch"),
            surface_source, model_id, list(surface_chains), [],
            {"ph": round(float(st.session_state[f"complex-analysis-db-pep-ph-{key}"]), 3)},
            key + "-pep", submit, project_id,
            schema_errors.get("analysis.pep_patch", "") or ("Choose at least one protein chain for the surface analysis." if not surface_chains else ""))


def _render_database_workflow(user: dict, workflow_id: str, schema: dict | None, source: dict,
                              model_id: str, binder_chains: list[str], target_chains: list[str],
                              parameters: dict, key: str, submit, project_id: str,
                              unavailable: str = "", *, structure_path: Path | None = None,
                              structure_context: dict | None = None,
                              structure_components: list[dict] | None = None,
                              viewer_key: str = "", ligand_selection: dict | None = None) -> None:
    title = {"analysis.biotite_interface": "Biotite Interface", "analysis.rosetta": "Rosetta Interface Analysis",
             "analysis.ipsae": "ipSAE & Structure Analysis", "analysis.pep_patch": "PEP-Patch",
             "analysis.plip": "PLIP interaction analysis", "analysis.pandamap": "PandaMap interaction analysis",
             "analysis.posebusters": "PoseBusters pose validation",
             "analysis.reference_placement": "Reference-frame placement",
             "analysis.reference_environment": "Reference-frame steric and environmental context",
             "analysis.tag_fusion_suitability": "Tag and Fusion Suitability",
             "analysis.tag_fusion_context": "Reference-context Tag and Fusion Suitability",
             "analysis.sequence_liability": "Sequence and Developability Liabilities",
             "analysis.binder_interface_chemistry": "BinderScope-style Interface Surface & Chemistry"}.get(workflow_id, workflow_id)
    st.markdown(f"##### {title}")
    if not schema:
        st.info(unavailable or "The registered workflow schema is unavailable.")
        return
    if unavailable:
        st.info(unavailable)
        return
    candidate_id = "database-" + hashlib.sha256((str(source["molecule_id"]) + ":" + str(source["structure_id"])).encode()).hexdigest()[:24]
    try:
        matches = _matching_database_analysis_jobs(user, workflow_id=workflow_id, source=source,
            candidate_id=candidate_id, model_id=model_id, binder_chains=list(binder_chains),
            target_chains=list(target_chains), parameters=parameters, ligand_selection=ligand_selection,
            expected_protocol=schema.get("assessment_protocol"))
    except (OSError, PermissionError, RuntimeError, ValueError, sqlite3.Error) as exc:
        st.info(f"Cannot verify analysis job status: {exc}")
        if st.button("Refresh analysis job index", key=f"refresh-analysis-index-{workflow_id}-{key}"):
            from mn_molecular_design_studio.job_catalog import reconcile
            reconcile(wait=False)
            st.rerun()
        return
    current = next((item for item in matches if item["status"] in {"completed", "succeeded"}), None)
    if current is None:
        current = next((item for item in matches if item["status"] not in {"failed", "cancelled", "canceled"}), None)
    if current:
        from mn_molecular_design_studio.complex_analysis_batch import protocol_matches
        if (current["status"] in {"completed", "succeeded"}
                and not protocol_matches((current.get("database_analysis") or {}).get("protocol") or {},
                                         schema.get("assessment_protocol") or {})):
            st.info("Historical analysis has unverified protocol compatibility; it is not reused as a current Evaluation.")
            return
        if workflow_id == "analysis.biotite_interface":
            _render_biotite_results(user, current, source_candidate_id=candidate_id, key=key,
                structure_path=structure_path, structure_context=structure_context,
                structure_components=structure_components, viewer_key=viewer_key)
        else:
            _render_contextual_results(user, current, workflow_id=workflow_id,
                                       source_candidate_id=candidate_id, key=key)
    elif matches:
        failed = matches[0]
        st.error(f"Previous matching analysis failed: {failed.get('error') or 'No per-item error was recorded.'}")
        if workflow_id == "analysis.biotite_interface":
            _render_biotite_results(user, failed, source_candidate_id=candidate_id, key=key,
                structure_path=structure_path, structure_context=structure_context,
                structure_components=structure_components, viewer_key=viewer_key)
        else:
            _render_contextual_results(user, failed, workflow_id=workflow_id,
                                       source_candidate_id=candidate_id, key=key)
    else:
        st.caption("No matching evaluation is saved for this exact Database structure, model, chain selection, and method settings.")
    if current and current["status"] in {"completed", "succeeded"}:
        return
    if current and current["status"] not in {"failed", "cancelled", "canceled"}:
        return
    if st.button("Retry analysis" if matches else "Run analysis", type="primary", key=f"complex-analysis-db-run-{workflow_id.replace('.', '-')}-{key}", disabled=submit is None):
        from mn_molecular_design_studio.complex_analysis_batch import analysis_chain_parameters
        st.session_state["complex-analysis-active"] = {
            "kind": "database", "molecule_id": source["molecule_id"],
            "structure_id": source["structure_id"], "evaluation_id": source.get("evaluation_id", ""),
        }
        try:
            submit(schema, workflow_id, {**parameters, "candidate_ids": [candidate_id],
                "selected_structure_path": source.get("selected_structure_path"),
                "selected_model_id": model_id, **analysis_chain_parameters(workflow_id, binder_chains, target_chains),
                "database_source": {
                    "molecule_id": source["molecule_id"], "structure_id": source["structure_id"],
                    "group_id": source["group_id"], "evaluation_id": source.get("evaluation_id", ""),
                    "model_id": model_id, "binder_chains": list(binder_chains),
                    "target_chains": list(target_chains),
                    **({"source_component_chain": str(parameters.get("source_component_chain") or "")}
                       if workflow_id in {"analysis.reference_placement", "analysis.tag_fusion_context"} else {}),
                    **({"ligand_selection": dict(ligand_selection)} if ligand_selection else {}),
                }, **({"ligand_selection": dict(ligand_selection)}
                      if ligand_selection and workflow_id in {"analysis.plip", "analysis.pandamap", "analysis.posebusters"} else {})}, project_id,
                {"_execution": {"placement": "local", "reassessment_of": matches[0]["run_id"]}} if matches else {}, title)
            st.rerun()
        except (OSError, RuntimeError, TypeError, ValueError) as exc:
            st.error(f"Could not start {title}: {exc}")


def _render_database_source(user: dict, opened: dict, *, workflow_schema=None, submit=None) -> None:
    from mn_molecular_design_studio.molecular_database import MolecularDatabase
    from mn_molecular_design_studio.complex_analysis_catalog import authorized_structure_catalog

    database = MolecularDatabase()
    st.caption("Analyses launched from Database use the existing Studio workflow queue on this workstation; the selected structure stays linked to its Database record without creating a prediction Results job.")
    query = st.text_input("Search Database complexes", key="complex-analysis-db-query", placeholder="Name or sequence")
    catalog_key = f"complex-analysis-db-catalog-{user.get('id', '')}"
    cached_catalog = st.session_state.get(catalog_key)
    same_query = isinstance(cached_catalog, dict) and str(cached_catalog.get("query") or "") == query
    refresh_catalog = st.button("Refresh Database structure choices", key="complex-analysis-db-refresh")
    load_next = st.button("Load next Database structure page", key="complex-analysis-db-next",
                          disabled=not (same_query and cached_catalog.get("next_offset") is not None))
    catalog = authorized_structure_catalog(
        database, user, session_state=st.session_state, cache_key=catalog_key,
        query=query, limit=100, advance=load_next, refresh=refresh_catalog,
        eligible=lambda item: bool(item.get("representations") and _has_database_analysis_complex(item)),
    )
    records = catalog["records"]
    st.caption(f"Loaded {len(records):,} authorized Database records with eligible complex structures"
               + ("; more eligible records may be available." if catalog.get("next_offset") is not None else "."))
    by_id = {str(item["id"]): item for item in records}
    requested_id = str(opened.get("molecule_id") or st.session_state.get("complex-analysis-molecule-id") or "")
    if requested_id and requested_id not in by_id:
        try:
            requested = database.list_molecules(user, molecule_ids=[requested_id], collection="all", limit=1)["items"]
        except PermissionError:
            requested = []
        if requested and requested[0].get("representations") and _has_database_analysis_complex(requested[0]):
            records = [requested[0], *records]
            by_id[requested_id] = requested[0]
    if not records:
        st.info("No accessible Database records with a saved multi-component structure or explicit target-conditioned structure provenance are available. Import a PDB/mmCIF structure in Database or open a completed Results structure.")
        return
    ids = list(dict.fromkeys(str(item["id"]) for item in records))
    labels = {
        str(item["id"]): (
            f"{item.get('name') or 'Complex'} · "
            + ("target in saved structure · " if item.get("source_target_provenance") else
               f"{len(item.get('components') or [])} registered components · ")
            + str(item["id"])[-8:]
        )
        for item in records
    }
    if requested_id not in ids:
        requested_id = ids[0]
    molecule_widget_key = "complex-analysis-molecule-id"
    if st.session_state.get(molecule_widget_key) not in ids:
        st.session_state[molecule_widget_key] = requested_id
    molecule_id = st.selectbox("Database record", ids, index=ids.index(requested_id), format_func=lambda value: labels[value], key=molecule_widget_key)
    rows = database.list_molecules(user, molecule_ids=[molecule_id], collection="all", limit=1)["items"]
    if not rows:
        st.warning("This Database record is unavailable or you do not have access to it.")
        return
    all_record = rows[0]
    all_evaluation_rows, all_evaluation_by_id = _evaluation_options(all_record)
    record = all_record
    evaluation_rows, evaluation_by_id = all_evaluation_rows, all_evaluation_by_id
    default_evaluation = str(opened.get("evaluation_id") or (record.get("selected_evaluation") or {}).get("id") or "")
    evaluation_ids = ["", *[row["id"] for row in evaluation_rows]]
    if default_evaluation not in evaluation_ids:
        default_evaluation = ""
    evaluation_widget_key = f"complex-analysis-evaluation-{molecule_id}"
    if st.session_state.get(evaluation_widget_key) not in evaluation_ids:
        st.session_state[evaluation_widget_key] = default_evaluation
    eval_col, structure_col = st.columns([1, 1])
    evaluation_id = eval_col.selectbox(
        "Evaluation / model metadata", evaluation_ids,
        index=evaluation_ids.index(default_evaluation),
        format_func=lambda value: "All saved structures · no evaluation selected" if not value else _evaluation_label(evaluation_by_id[value][1]),
        key=evaluation_widget_key,
    )
    if evaluation_id:
        record_rows = database.list_molecules(user, molecule_ids=[molecule_id], collection="all", limit=1,
                                              evaluation_filters={"evaluation_id": evaluation_id})["items"]
        if not record_rows:
            st.warning("This evaluation is unavailable or you do not have access to its source group.")
            return
        record = record_rows[0]
        evaluation_rows, evaluation_by_id = _evaluation_options(record)
    # list_molecules applies the canonical evaluation-to-representation
    # association and ACL filtering. Use its selected structures when an
    # evaluation is chosen instead of trying to reconstruct that contract here.
    representations = list(record.get("structures") or []) if evaluation_id else list(
        record.get("representations") or record.get("structures") or []
    )
    if not representations:
        if evaluation_id:
            st.info("The selected evaluation has no unambiguous saved structure association. Choose another evaluation or clear the selector to view saved structures separately.")
        else:
            st.info("No structure is attached to this Database record.")
        return
    representation_ids = [str(item["id"]) for item in representations]
    requested_structure = str(opened.get("structure_id") or st.session_state.get(f"complex-analysis-structure-id-{molecule_id}") or "")
    if requested_structure not in representation_ids:
        requested_structure = representation_ids[0]
    structure_widget_key = f"complex-analysis-structure-id-{molecule_id}"
    if st.session_state.get(structure_widget_key) not in representation_ids:
        st.session_state[structure_widget_key] = requested_structure
    structure_id = structure_col.selectbox(
        "Structure / model / pose", representation_ids,
        index=representation_ids.index(requested_structure),
        format_func=lambda value: next(_representation_label(item) for item in representations if str(item["id"]) == value),
        key=structure_widget_key,
    )
    structure = next(item for item in representations if str(item["id"]) == structure_id)
    try:
        path = database.structure_file(user, structure_id)
    except (OSError, ValueError, PermissionError) as exc:
        st.warning(f"The selected structure is unavailable: {exc}")
        return
    st.caption(
        f"Database record: {record.get('name') or record['id']} · {record['id']} · "
        f"Structure {structure_id} · {structure.get('label') or path.name}"
    )
    component_rows = component_inventory(path.read_bytes(), path.name, all_record.get("components"))
    source_metadata = all_record.get("metadata") if isinstance(all_record.get("metadata"), dict) else {}
    target_notice = _target_structure_notice(source_metadata)
    if target_notice:
        st.info(target_notice)
    if source_metadata.get("target_chain_inference") and (source_metadata.get("target_fasta") or source_metadata.get("target_pdb")):
        st.caption(
            "This Database molecular identity stores fewer chains than its target-conditioned source record. "
            "The selected saved structure is parsed below; choose the target chain from its coordinates explicitly. "
            "The source record preserves a target sequence relationship, not a separate linked target molecule."
        )
    _show_components(component_rows)
    viewer_key = f"complex-analysis-db-viewer-{structure_id}"
    structure_context = structure.get("context") if isinstance(structure.get("context"), dict) else {}
    from mn_molecular_design_studio.scientific_settings import model_identity
    model_id = model_identity(structure_context.get("model_id"), structure_context.get("model_index"),
                              structure_context.get("model"), structure.get("model_index"), structure_id)
    associated = []
    for observation in all_record.get("observations", []):
        evaluation = next((item for known_observation, item in all_evaluation_by_id.values()
                           if known_observation.get("id") == observation.get("id")), None)
        if evaluation and _evaluation_matches_structure(evaluation, structure):
            associated.append((observation, evaluation))
    if evaluation_id and not any(evaluation.get("id") == evaluation_id for _, evaluation in associated):
        selected = all_evaluation_by_id.get(evaluation_id)
        if selected:
            associated.append(selected)
    inspector_data = _database_inspector_data(user, all_record, structure, path, component_rows, associated)
    selected_finding_key = f"{viewer_key}-inspector-evidence-select"
    selected_finding_id = st.session_state.get(selected_finding_key)
    selected_finding = next((item for item in inspector_data["evidence"]
                             if item["id"] == selected_finding_id), None)
    if selected_finding:
        st.session_state[f"{viewer_key}-inspector-highlights"] = selected_finding["selectors"]
    reference_key = hashlib.sha256(f"database:{all_record['id']}:{structure['id']}".encode()).hexdigest()[:14]
    paired_reference = _paired_reference_for_viewer(
        database, user, all_record, structure, component_rows, reference_key,
    )
    if paired_reference:
        saved = paired_reference["saved"]
        placement = _reference_placement_match(user, all_record, structure, saved, model_id)
        if placement:
            st.markdown("#### Reference-frame structural placement")
            st.caption(
                f"Cα fit · {placement['result']['aligned_pair_count']} aligned pairs · "
                f"source coverage {placement['result']['source_coverage']:.3f} · "
                f"reference coverage {placement['result']['reference_coverage']:.3f} · "
                f"RMSD {placement['result']['rmsd_angstrom']:.3f} Å"
            )
            st.caption(f"{placement['result']['method']} · {placement['result']['method_version']} · {placement['result']['mapping_method']}")
            with st.expander("Placement limitations"):
                for limitation in placement["result"].get("limitations") or []:
                    st.write(f"• {limitation}")
            _show_reference_placement_overlay(
                paired_reference["path"], path, structure.get("context") or {},
                placement["result"],
                key=f"complex-analysis-reference-placement-viewer-{reference_key}-{paired_reference['source_chain']}",
                highlighted_selections=_active_biotite_highlights(viewer_key),
                reference_highlighted_selections=st.session_state.get(f"{viewer_key}-reference-highlights"),
            )
            st.caption("Mol*: the complete associated reference is shown first in grey; the selected complex is second, transformed into the reference frame and colored by chain.")
            st.caption("The selected complex is transformed for display into the associated reference frame. Its Database coordinates and existing evaluations are unchanged. This fit does not establish experimental binding, biological accessibility, clashes, or glycan context.")
        else:
            source_view, reference_view = st.columns(2)
            with source_view:
                st.markdown("#### Selected complex")
                st.caption(f"Source {all_record.get('name') or all_record['id']} · {structure_id} · component {paired_reference['source_chain']}")
                _show_structure(path, structure.get("context") or {}, component_rows, key=viewer_key,
                                highlighted_selections=_active_biotite_highlights(viewer_key))
            with reference_view:
                st.markdown("#### Associated reference structure")
                st.caption(
                    f"{saved.get('reference_record_name') or saved.get('reference_molecule_id')} · "
                    f"{saved.get('reference_structure_label') or paired_reference['reference_id']} · "
                    f"{paired_reference['reference_id']} · component {saved.get('reference_component_chain')}"
                )
                _show_structure(
                    paired_reference["path"], {}, paired_reference["components"],
                    key=f"complex-analysis-reference-viewer-{reference_key}-{paired_reference['source_chain']}",
                )
            st.caption("No completed placement evaluation matches this exact source/reference association. These structures are displayed in their original, separate coordinate frames.")
        _render_database_inspector(user, inspector_data, viewer_key, record=all_record,
                                   structure=structure, components=component_rows, model_id=model_id,
                                   reference_view_context={"paired_reference": paired_reference, "placement": placement})
    else:
        structure_col, inspector_col = st.columns([1.15, 1])
        with structure_col:
            _show_structure(path, structure.get("context") or {}, component_rows, key=viewer_key,
                            highlighted_selections=_active_biotite_highlights(viewer_key))
        with inspector_col:
            _render_database_inspector(user, inspector_data, viewer_key, record=all_record,
                                       structure=structure, components=component_rows, model_id=model_id,
                                       reference_view_context=None)
    _render_inspector_tracks(inspector_data, viewer_key)
    _show_component_groups(component_rows, structure_id)
    _render_database_analysis(user, record, structure, path, component_rows, evaluation_id,
                              workflow_schema=workflow_schema, submit=submit)
    st.markdown("#### Existing evaluations associated with this structure")
    if associated:
        for observation, evaluation in associated:
            with st.expander(_evaluation_label(evaluation), expanded=evaluation.get("id") == evaluation_id):
                observation_context = observation.get("context") if isinstance(observation.get("context"), dict) else {}
                candidate_context = observation_context.get("candidate") if isinstance(observation_context.get("candidate"), dict) else {}
                raw_candidate_metadata = candidate_context.get("raw_metadata") if isinstance(candidate_context.get("raw_metadata"), dict) else {}
                if not _evaluation_matches_structure(evaluation, structure):
                    st.caption("This evaluation was selected explicitly, but the saved records do not establish an exact structure association.")
                st.json({
                    "evaluation": {key: value for key, value in evaluation.items() if key != "observation_id"},
                    "stored_metrics": observation.get("metrics") or {},
                    "normalized_scores": observation.get("normalized_output", {}).get("metrics", {}),
                    "model_confidence_outputs": raw_candidate_metadata.get("confidence_outputs") or {},
                    "evaluation_context": observation_context,
                })
    else:
        st.info("No saved evaluation is associated with this structure. No confidence or scientific measurements are inferred.")
    with st.expander("Structure provenance", expanded=False):
        st.json({
            "record_id": record["id"],
            "record_origin": record.get("origin") or {},
            "source_groups": record.get("groups") or [],
            "record_sources": record.get("sources") or [],
            "structure_id": structure_id,
            "structure_group_id": structure.get("group_id"),
            "structure_format": structure.get("format"),
            "structure_context": structure.get("context") or {},
        })
    _show_human_review(
        database, user, source_kind="database", subject_id=molecule_id,
        structure_ref=structure_id, evaluation_id=evaluation_id,
    )


def _render_results_source(user: dict, opened: dict, *, workflow_schema=None, project_selector=None, submit=None) -> None:
    """Render one explicitly opened, authorized Results artifact without catalog discovery."""
    from pathlib import PurePosixPath
    from mn_molecular_design_studio.auth import can_access_job
    from mn_molecular_design_studio.job_catalog import get_job_access_record, job_artifacts
    from mn_molecular_design_studio.runtime import runs_root

    source_ref = (str(opened.get("task_group") or ""), str(opened.get("run_id") or ""))
    artifact_ref = (str(opened.get("artifact_task_group") or opened.get("task_group") or ""),
                    str(opened.get("artifact_run_id") or opened.get("run_id") or ""))
    reference = str(opened.get("relative_path") or "")
    relative = PurePosixPath(reference)
    # Cross-run navigation persists a portable runs:/// reference so that the
    # parent campaign and the child structure-producing job remain distinct.
    # Resolve that reference only after both jobs have passed authorization.
    if reference.startswith("runs:///"):
        portable_parts = PurePosixPath(reference.removeprefix("runs:///"))
        if (portable_parts.is_absolute() or len(portable_parts.parts) < 3
                or ".." in portable_parts.parts or tuple(portable_parts.parts[:2]) != artifact_ref):
            raise ValueError("The portable Results Inspector link does not match its producing job.")
        artifact_relative = PurePosixPath(*portable_parts.parts[2:])
    else:
        artifact_relative = relative
    if (not all(source_ref) or not all(artifact_ref) or artifact_relative.is_absolute()
            or not artifact_relative.parts or ".." in artifact_relative.parts):
        raise ValueError("The Results Inspector link is incomplete or invalid.")

    source_access = get_job_access_record(*source_ref)
    artifact_access = get_job_access_record(*artifact_ref)
    if (not source_access or not artifact_access or not can_access_job(user, source_access)
            or not can_access_job(user, artifact_access)):
        raise PermissionError("The source job or structure-producing job is unavailable or not authorized.")

    root = Path(runs_root()).expanduser().resolve()
    source_run_dir = (root / str(source_access.get("run_relative_path") or "")).resolve(strict=True)
    artifact_run_dir = (root / str(artifact_access.get("run_relative_path") or "")).resolve(strict=True)
    source_run_dir.relative_to(root)
    artifact_run_dir.relative_to(root)
    structure_path = (artifact_run_dir / Path(*artifact_relative.parts)).resolve(strict=True)
    structure_path.relative_to(artifact_run_dir)
    structure_path.relative_to(root)

    # Revalidate the artifact against the producing job's registered manifest.
    registered = False
    for artifact in job_artifacts(*artifact_ref):
        registered_relative = PurePosixPath(str(artifact.get("relative_path") or ""))
        if registered_relative.is_absolute() or ".." in registered_relative.parts:
            continue
        try:
            registered_path = (artifact_run_dir / Path(*registered_relative.parts)).resolve(strict=True)
            registered_path.relative_to(artifact_run_dir)
        except (OSError, ValueError):
            continue
        if registered_relative == artifact_relative and registered_path == structure_path:
            registered = True
            break
    if not registered:
        raise PermissionError("The selected structure is not a registered artifact of its authorized Results job.")

    candidate = opened.get("candidate") if isinstance(opened.get("candidate"), dict) else {}
    candidate_id = str(opened.get("candidate_id") or "")
    if not candidate and candidate_id:
        from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
        candidate_path = source_run_dir / "artifacts" / "normalized_candidates" / "candidates.jsonl"
        if candidate_path.is_file():
            candidate = next((row for row in read_candidates(candidate_path)
                              if str(row.get("candidate_id") or "") == candidate_id), {})
    if candidate_id and str(candidate.get("candidate_id") or "") != candidate_id:
        raise ValueError("The Results Inspector candidate no longer matches its saved source reference.")
    if candidate:
        from mn_molecular_design_studio import cofolding_results
        expected = cofolding_results.structure_for(source_run_dir, candidate)
        if expected is None or expected.resolve(strict=True) != structure_path:
            raise ValueError("The selected structure does not match the explicitly selected Results candidate.")

    digest = _structure_digest(structure_path)
    if not digest or (opened.get("structure_sha256") and digest != str(opened["structure_sha256"])):
        raise ValueError("The selected Results structure changed or failed its saved digest check.")
    metadata = _read_json(source_run_dir / "metadata.json")
    candidate_metadata = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
    model_identity = (candidate.get("model_id") or candidate.get("model_index")
                      or candidate_metadata.get("model_id") or candidate_metadata.get("model_index")
                      or opened.get("model_index") or "model not specified")

    st.caption(f"Results: {source_ref[0]}/{source_ref[1]} · {metadata.get('tool') or metadata.get('display_tool') or 'Unknown tool'} · {reference}")
    if artifact_ref != source_ref:
        st.caption(f"Structure produced by linked job {artifact_ref[0]}/{artifact_ref[1]}.")
    st.caption(f"Candidate/model: {candidate_id or 'not recorded'} · {model_identity}")
    st.caption(f"SHA-256: {digest}")
    content = structure_path.read_bytes()
    components = component_inventory(content, structure_path.name, candidate_metadata.get("components"))
    _show_components(components)
    viewer_key = "complex-analysis-targeted-results-" + hashlib.sha256(
        f"{source_ref}:{artifact_ref}:{relative.as_posix()}:{digest}".encode()).hexdigest()[:16]
    _show_structure(structure_path, candidate_metadata, components, key=viewer_key,
                    highlighted_selections=_active_biotite_highlights(viewer_key))
    _show_component_groups(components, viewer_key)
    if candidate:
        st.markdown("#### Saved candidate evidence")
        st.json({
            "candidate_id": candidate_id,
            "stage": candidate.get("stage"),
            "source_tool": candidate.get("source_tool") or metadata.get("tool"),
            "stored_metrics": candidate.get("metrics") or {},
            "confidence_outputs": candidate_metadata.get("confidence_outputs") or {},
        })
    with st.expander("Results provenance", expanded=False):
        st.json({
            "source_job": f"{source_ref[0]}/{source_ref[1]}",
            "source_project_id": source_access.get("project_id"),
            "source_owner_id": source_access.get("owner_id"),
            "structure_job": f"{artifact_ref[0]}/{artifact_ref[1]}",
            "candidate_id": candidate_id,
            "model_id": model_identity,
            "artifact_path": reference,
            "structure_sha256": digest,
            "candidate_source_metadata": candidate_metadata,
        })
    st.caption("Read-only Results inspection. To analyze or shortlist this design, explicitly save the selected candidate/model to the Molecular Database from its Results page.")


def _render_complex_analysis_inspector(user: dict, *, workflow_schema=None, project_selector=None, submit=None) -> None:
    """Render one explicitly selected Database structure or Results artifact."""
    opened = st.session_state.pop("complex-analysis-open", None)
    if opened:
        st.session_state["complex-analysis-active"] = opened
    active = st.session_state.get("complex-analysis-active") or {}
    st.title("Scientific Inspector")
    if st.button("← Scientific Work", key="complex-analysis-return-database"):
        st.session_state["complex-analysis-workspace-view"] = "database"
        st.session_state.pop("complex-analysis-source-kind", None)
        st.session_state.pop("complex-analysis-active", None)
        st.session_state.pop("complex-analysis-detail-visible", None)
        st.rerun()
    kind = str(active.get("kind") or "")
    st.caption("Inspect the exact saved structure and its scientific provenance.")
    try:
        if kind == "database":
            _render_database_source(user, active,
                                    workflow_schema=workflow_schema, submit=submit)
        elif kind == "remote_results":
            _render_remote_results_source(user, active,
                                          workflow_schema=workflow_schema,
                                          project_selector=project_selector, submit=submit)
        elif kind == "results":
            _render_results_source(user, active if active.get("kind") == "results" else {},
                                   workflow_schema=workflow_schema,
                                   project_selector=project_selector, submit=submit)
        else:
            st.info("Open an individual structure from Results, or select a saved Database structure to inspect it here.")
    except PermissionError:
        st.warning("This record or job is unavailable or you do not have access to it.")
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        st.error(f"Could not load this structure source: {html.escape(str(exc))}")


def render_complex_analysis_workspace(user: dict, *, workflow_schema=None, project_selector=None, submit=None) -> None:
    """Render the Database-first Select → Analyze → Compare → Shortlist workspace.

    Results catalog discovery belongs to the Results workspace. Complex Analysis
    only enters a Results path for an explicitly opened artifact deep link.
    """
    opened = bool(st.session_state.get("complex-analysis-open"))
    active = st.session_state.get("complex-analysis-active")
    view = str(st.session_state.get("complex-analysis-workspace-view") or "")
    if opened or (view == "inspector" and isinstance(active, dict)):
        _render_complex_analysis_inspector(user, workflow_schema=workflow_schema,
                                           project_selector=project_selector, submit=submit)
        return
    st.session_state["complex-analysis-workspace-view"] = "database"
    st.title("Complex Analysis")
    st.caption("Open saved work by project and target/system. Review recorded evidence, assess compatible inputs, inspect structures and document selections.")
    from mn_molecular_design_studio.complex_analysis_collection_ui import render_collection_workspace
    render_collection_workspace(user, workflow_schema=workflow_schema, submit=submit)
