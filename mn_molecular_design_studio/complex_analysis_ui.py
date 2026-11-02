"""Tool-independent structure and provenance workspace for saved complexes."""
from __future__ import annotations

import gzip
import hashlib
import html
import io
import json
from pathlib import Path
import streamlit as st


STRUCTURE_SUFFIXES = {".pdb", ".ent", ".cif", ".mmcif", ".pdbqt", ".sdf", ".mol"}
VIEWER_MAX_BYTES = 50 * 1024 * 1024


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
    st.session_state["page"] = "Complex Analysis"
    st.rerun()


def open_result_complex(run_dir: Path | str, structure_path: Path | str, candidate: dict | None = None) -> None:
    """Navigate to a saved Results structure using job and artifact identity."""
    run_dir = Path(run_dir).expanduser().resolve()
    structure_path = Path(structure_path).expanduser().resolve()
    try:
        relative_path = structure_path.relative_to(run_dir).as_posix()
    except ValueError as exc:
        raise ValueError("The selected structure is outside its source job.") from exc
    candidate = candidate if isinstance(candidate, dict) else {}
    raw_metadata = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
    task_group, run_id = run_dir.parent.name, run_dir.name
    st.session_state["complex-analysis-results-job"] = f"{task_group}/{run_id}"
    st.session_state["complex-analysis-results-structure"] = relative_path
    st.session_state["complex-analysis-open"] = {
        "kind": "results",
        "task_group": task_group,
        "run_id": run_id,
        "candidate_id": str(candidate.get("candidate_id") or ""),
        "model_index": str(candidate.get("model_index") or raw_metadata.get("model_index") or ""),
        "relative_path": relative_path,
    }
    st.session_state["complex-analysis-source-kind"] = "Results"
    st.session_state["page"] = "Complex Analysis"
    st.rerun()


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
        limit=300,
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
                "candidate": candidate,
                "candidate_id": str(candidate.get("candidate_id") or ""),
                "label": str(artifact.get("label") or path.name),
            })
        except (OSError, ValueError):
            continue
    # A normalized candidate can refer to a saved structure which has not yet
    # been added to the artifact catalog. Keep it usable if it remains inside
    # the authorized job directory.
    known = {str(item["path"]) for item in rows}
    for candidate in candidates:
        structure = cofolding_results.structure_for(run_dir, candidate)
        if structure is None:
            continue
        try:
            resolved = structure.resolve(strict=True)
            resolved.relative_to(run_dir)
            resolved.relative_to(root)
        except (OSError, ValueError):
            continue
        if str(resolved) in known or _structure_suffix(resolved) not in STRUCTURE_SUFFIXES:
            continue
        rows.append({
            "path": resolved,
            "relative_path": resolved.relative_to(run_dir).as_posix(),
            "artifact_path": resolved.relative_to(root).as_posix(),
            "candidate": candidate,
            "candidate_id": str(candidate.get("candidate_id") or ""),
            "label": resolved.name,
        })
        known.add(str(resolved))
    rows.sort(key=lambda item: (item["relative_path"].casefold(), item["candidate_id"].casefold()))
    return rows


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


def _matching_biotite_jobs(user: dict, *, source_task_group: str, source_run_id: str,
                           candidate_id: str, structure_digest: str,
                           binder_chains: list[str], target_chains: list[str],
                           contact_distance: float, model_id: str = "") -> list[dict]:
    """Find reusable Biotite evaluations for this exact source pose and chain mapping."""
    from mn_molecular_design_studio.auth import accessible_project_ids, can_access_job
    from mn_molecular_design_studio.job_catalog import get_job_access_record, search_jobs

    rows = search_jobs(
        tools=["biotite_interface"], task_group="analysis", limit=1000,
        access_owner_id=str(user.get("id") or ""),
        accessible_project_ids=accessible_project_ids(user),
        access_admin=user.get("role") == "admin", sort_by="active_desc",
    )["items"]
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
        saved_model_id = str(selection.get("model_id") or "")
        if model_id and saved_model_id and model_id != saved_model_id:
            continue
        desired_binder, desired_target = sorted(map(str, binder_chains)), sorted(map(str, target_chains))
        selected_item = next((item for item in manifest.get("items", [])
                              if str(item.get("candidate_id") or "") == candidate_id
                              and str(item.get("structure_digest") or "") == structure_digest
                              and sorted(map(str, item.get("binder_chains") or [])) == desired_binder
                              and sorted(map(str, item.get("target_chains") or [])) == desired_target), None)
        if selected_item is None:
            continue
        status = str(row.get("status") or metadata.get("status") or "unknown").casefold()
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


def _analysis_evaluation_id(run_dir: Path, user: dict, source_candidate_id: str) -> str:
    """Resolve the persistent Evaluation descriptor created by normal Results ingestion."""
    try:
        from mn_molecular_design_studio.evaluations import evaluation_descriptor
        from mn_molecular_design_studio.molecular_database_results import list_run_molecules
        records = list_run_molecules(run_dir, user)
    except (ImportError, OSError, PermissionError, ValueError):
        return ""
    output_candidate_id = source_candidate_id + "_biotite_interface"
    for record in records:
        for observation in record.get("observations", []):
            try:
                descriptor = evaluation_descriptor(observation)
            except (KeyError, TypeError, ValueError):
                continue
            candidate = (observation.get("context") or {}).get("candidate") or {}
            if descriptor.get("evaluator") == "biotite_interface" and (
                descriptor.get("candidate_id") == output_candidate_id
                or candidate.get("candidate_id") == output_candidate_id
            ):
                return str(descriptor.get("id") or "")
    return ""


def _render_biotite_results(user: dict, evaluation: dict, *, source_candidate_id: str, key: str) -> None:
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
        evaluation_id = _analysis_evaluation_id(run_dir, user, source_candidate_id)
        if evaluation_id:
            st.success(f"Saved Evaluation · {evaluation_id}")
        else:
            st.caption("The workflow Results are complete. The Database evaluation projection will appear when the run is indexed.")
    with st.expander("Analysis provenance", expanded=False):
        manifest = evaluation.get("manifest") or {}
        item = evaluation.get("selected_item") or {}
        st.json({
            "evaluation_run": f"{evaluation.get('task_group')}/{evaluation.get('run_id')}",
            "workflow": "analysis.biotite_interface",
            "source_job": f"{evaluation.get('metadata', {}).get('upstream_task_group')}/{evaluation.get('metadata', {}).get('upstream_run_id')}",
            "source_candidate_id": item.get("candidate_id"),
            "selected_structure_sha256": item.get("structure_digest"),
            "binder_chains": item.get("binder_chains") or [],
            "target_chains": item.get("target_chains") or [],
            "method": manifest.get("tool"),
            "parameters": manifest.get("settings") or {},
            "source_candidate_structure": (manifest.get("selection") or {}).get("structure_path") or "source candidate structure",
            "source_job_code": evaluation.get("metadata", {}).get("upstream_job_code"),
        })


def _render_biotite_workspace(user: dict, job: dict, selected: dict, components: list[dict], *,
                              workflow_schema, project_selector, submit) -> None:
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
    st.caption("ipSAE and Interface Clustering are already available in Studio. Contextual selection and launch for those workflows will follow this first Biotite integration; their existing inputs remain available in All Tools.")
    st.markdown("##### Protein–ligand")
    if protein_ligand:
        st.caption("PLIP, PandaMap, and PoseBusters are registered Studio workflows. Contextual launch for the selected pose is a subsequent integration phase.")
    else:
        st.caption("Not applicable to the selected components.")
    st.markdown("##### Protein–DNA/RNA")
    if nucleic_acids:
        st.caption("The selected structure contains nucleic acid. No registered contextual analysis in this Biotite Interface slice evaluates protein–nucleic-acid interactions; capabilities will be listed as unavailable until a compatible existing workflow is connected.")
    else:
        st.caption("Not applicable to the selected components.")
    if not ppi_applicable:
        return

    candidate = selected.get("candidate") if isinstance(selected.get("candidate"), dict) else {}
    candidate_id = str(candidate.get("candidate_id") or "")
    path = Path(selected["path"]).expanduser().resolve()
    task_group = str(job.get("task_group") or "")
    run_id = str(job.get("run_id") or "")
    source_run_dir = Path(str(job.get("run_dir") or "")).expanduser().resolve()
    candidate_file = source_run_dir / "artifacts" / "normalized_candidates" / "candidates.jsonl"
    workflow_id = "analysis.biotite_interface"
    schema = None
    unavailable = ""
    try:
        if workflow_schema is None:
            raise RuntimeError("The workflow registry is unavailable in this page context.")
        schema = workflow_schema(workflow_id)
    except Exception as exc:
        unavailable = f"The existing workflow could not be loaded: {exc}"
    if not candidate_id:
        unavailable = "This structure artifact has no normalized candidate identity in the selected Results job."
    elif not candidate_file.is_file():
        unavailable = "The selected Results job has no normalized candidate collection for the existing workflow."
    elif _structure_suffix(path) not in {".pdb", ".ent", ".cif", ".mmcif"}:
        unavailable = "Biotite Interface requires a PDB or mmCIF protein complex artifact."
    elif not task_group or not run_id:
        unavailable = "The source Results job identity is incomplete."
    elif not str(job.get("project_id") or ""):
        unavailable = "The source Results job has no project assignment; Studio requires an authorized project for a new workflow."

    candidate_metadata = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
    candidate_binder = [str(value) for value in candidate.get("binder_chains") or []]
    candidate_target = [str(value) for value in candidate.get("target_chains") or []]
    if not candidate_binder:
        candidate_binder = [str(value) for value in candidate_metadata.get("binder_chains") or []]
    if not candidate_target:
        candidate_target = [str(value) for value in candidate_metadata.get("target_chains") or []]
    candidate_binder = [value for value in candidate_binder if value in protein_chains]
    candidate_target = [value for value in candidate_target if value in protein_chains and value not in candidate_binder]
    key = hashlib.sha256(f"{task_group}/{run_id}:{selected['relative_path']}".encode()).hexdigest()[:14]
    a_key, b_key = f"complex-analysis-biotite-a-{key}", f"complex-analysis-biotite-b-{key}"
    st.markdown("##### Configure Biotite Interface")
    st.caption(f"Exact model/artifact: {candidate.get('candidate_id') or path.name} · {candidate.get('model_index') or candidate_metadata.get('model_index') or path.name}")
    first, second = st.columns(2)
    if a_key not in st.session_state:
        st.session_state[a_key] = candidate_binder
    if b_key not in st.session_state:
        st.session_state[b_key] = candidate_target
    binder_chains = first.multiselect("Interface chain group A", protein_chains,
        key=a_key, help="The existing workflow calls this binder_chains; this label is only a group name and does not assert design or binding.")
    target_chains = second.multiselect("Interface chain group B", protein_chains,
        key=b_key, help="May contain one or more protein chains.")
    parameters = {str(row.get("name")): row for row in (schema or {}).get("parameters", [])}
    default_distance = float(parameters.get("contact_distance", {}).get("default", 5.0))
    contact_distance = st.number_input("Contact distance (Å)", min_value=0.1, max_value=20.0,
        value=default_distance, step=0.1, key=f"complex-analysis-biotite-cutoff-{key}")

    invalid_selection = (not binder_chains or not target_chains or bool(set(binder_chains) & set(target_chains)))
    if invalid_selection:
        st.info("Choose two nonempty, disjoint protein chain groups to enable the existing interface workflow.")
    if unavailable:
        st.warning(f"Biotite Interface unavailable: {unavailable}")
        return

    project_prefix = f"complex-analysis-biotite-{key}"
    source_project = str(job.get("project_id") or "")
    project_key = f"{project_prefix}-project-id"
    if project_selector:
        if project_key not in st.session_state:
            st.session_state[project_key] = source_project
        analysis_project = project_selector(project_prefix, label="Analysis project", allow_create=False)
    else:
        analysis_project = source_project
        st.caption(f"Analysis project: {source_project}")
    if analysis_project != source_project:
        st.warning("Select the source Results job's project. Studio requires analyses to remain in the source project's access boundary.")

    distance_digest = _structure_digest(path)
    model_id = str(candidate.get("model_index") or candidate_metadata.get("model_index") or path.name)
    existing = _matching_biotite_jobs(user, source_task_group=task_group, source_run_id=run_id,
        candidate_id=candidate_id, structure_digest=distance_digest,
        binder_chains=binder_chains, target_chains=target_chains,
        contact_distance=float(contact_distance), model_id=model_id) if not invalid_selection and distance_digest else []
    current = next((item for item in existing if item["status"] in {"completed", "succeeded"}), None)
    if current is None:
        current = next((item for item in existing if item["status"] not in {"failed", "cancelled", "canceled"}), None)
    if current:
        _render_biotite_results(user, current, source_candidate_id=candidate_id, key=key)
    elif existing:
        failed = existing[0]
        st.error(f"Previous matching analysis failed: {failed.get('error') or 'No per-item error was recorded.'}")
        _render_biotite_results(user, failed, source_candidate_id=candidate_id, key=key)
    else:
        st.caption("No matching Biotite Interface evaluation is saved for this exact structure, model, chain mapping, and contact distance.")

    if current and current["status"] in {"completed", "succeeded"}:
        return
    button_label = "Retry analysis" if existing else "Run analysis"
    disabled = bool(invalid_selection or not analysis_project or analysis_project != source_project or submit is None or schema is None)
    if st.button(button_label, type="primary", disabled=disabled, key=f"complex-analysis-run-biotite-{key}"):
        if not submit:
            st.warning("Workflow submission is unavailable in this page context.")
            return
        st.session_state["complex-analysis-active"] = {
            "kind": "results", "task_group": task_group, "run_id": run_id,
            "candidate_id": candidate_id, "model_index": model_id,
            "relative_path": selected["relative_path"],
        }
        workflow_parameters = {
            "candidate_ids": [candidate_id],
            "selected_structure_path": str(path),
            "selected_model_id": model_id,
            "binder_chains": list(binder_chains),
            "target_chains": list(target_chains),
            "contact_distance": float(contact_distance),
        }
        source_job = f"{task_group}/{run_id}"
        before = st.session_state.get("last_submitted_job")
        try:
            submit(schema, workflow_id, {**workflow_parameters, "source_job": source_job},
                   str(analysis_project), {}, "Biotite Interface")
        except (OSError, RuntimeError, TypeError, ValueError) as exc:
            st.error(f"Could not start Biotite Interface: {exc}")
            return
        if st.session_state.get("last_submitted_job") != before:
            st.rerun()


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


def _show_structure(path: Path, context: dict, components: list[dict], *, key: str) -> None:
    if path.stat().st_size > VIEWER_MAX_BYTES:
        st.info("This structure is over 50 MB. Download it from its source results or Database record to inspect it locally.")
        return
    try:
        text = _viewer_text(path, context)
        from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer import molstar_custom_component
        from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer.dataclasses import StructureVisualization

        molstar_custom_component(
            [StructureVisualization(pdb=text, color="chain-id", representation_type="cartoon+ball-and-stick")],
            key=key, height=650, show_controls=True, show_structure_tools=False,
        )
    except (ImportError, OSError, ValueError, TypeError) as exc:
        st.error(f"Could not display this structure in Mol*: {exc}")


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
    review = database.get_complex_review(
        user, source_kind=source_kind, subject_id=subject_id,
        structure_ref=structure_ref, evaluation_id=evaluation_id,
    )
    token = hashlib.sha256("\0".join((source_kind, subject_id, structure_ref, evaluation_id)).encode()).hexdigest()[:16]
    status_key = f"complex-analysis-review-status-{token}"
    note_key = f"complex-analysis-review-note-{token}"
    if status_key not in st.session_state:
        st.session_state[status_key] = review["status"]
    if note_key not in st.session_state:
        st.session_state[note_key] = review["note"]
    st.markdown("#### Human review")
    st.caption("Your review is saved persistently under your account, separately from calculated evaluations.")
    st.selectbox("Review status", ["Unreviewed", "Needs review", "Reviewed", "Not applicable"], key=status_key)
    st.text_area("Review notes", key=note_key, max_chars=10_000,
                 placeholder="Record observations or questions for a human reviewer.")
    if st.button("Save review", key=f"complex-analysis-review-save-{token}"):
        saved = database.save_complex_review(
            user, source_kind=source_kind, subject_id=subject_id,
            structure_ref=structure_ref, evaluation_id=evaluation_id,
            status=st.session_state[status_key], note=st.session_state[note_key],
        )
        st.success(f"Review saved · {saved['updated_at']}")


def _render_database_source(user: dict, opened: dict) -> None:
    from mn_molecular_design_studio.molecular_database import MolecularDatabase

    database = MolecularDatabase()
    query = st.text_input("Search Database complexes", key="complex-analysis-db-query", placeholder="Name or sequence")
    page = database.list_molecules(user, collection="all", query=query, limit=1000)
    records = [item for item in page["items"] if item.get("representations") and len(item.get("components") or []) >= 2]
    by_id = {str(item["id"]): item for item in records}
    requested_id = str(opened.get("molecule_id") or st.session_state.get("complex-analysis-molecule-id") or "")
    if requested_id and requested_id not in by_id:
        try:
            requested = database.list_molecules(user, molecule_ids=[requested_id], collection="all", limit=1)["items"]
        except PermissionError:
            requested = []
        if requested and len(requested[0].get("components") or []) >= 2 and requested[0].get("representations"):
            records = [requested[0], *records]
            by_id[requested_id] = requested[0]
    if not records:
        st.info("No accessible Database records with a saved multi-component structure are available. Import a PDB/mmCIF structure in Database or open a completed Results structure.")
        return
    ids = list(dict.fromkeys(str(item["id"]) for item in records))
    labels = {str(item["id"]): f"{item.get('name') or 'Complex'} · {len(item.get('components') or [])} components · {str(item['id'])[-8:]}" for item in records}
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
    record = rows[0]
    evaluation_rows, evaluation_by_id = _evaluation_options(record)
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
    component_rows = component_inventory(path.read_bytes(), path.name, record.get("components"))
    _show_components(component_rows)
    _show_structure(path, structure.get("context") or {}, component_rows, key=f"complex-analysis-db-viewer-{structure_id}")
    _show_component_groups(component_rows, structure_id)
    associated = []
    for observation in record.get("observations", []):
        evaluation = next((item for known_observation, item in evaluation_by_id.values()
                           if known_observation.get("id") == observation.get("id")), None)
        if evaluation and _evaluation_matches_structure(evaluation, structure):
            associated.append((observation, evaluation))
    if evaluation_id and not any(evaluation.get("id") == evaluation_id for _, evaluation in associated):
        selected = evaluation_by_id.get(evaluation_id)
        if selected:
            associated.append(selected)
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
    jobs = _accessible_result_jobs(user)
    if not jobs:
        st.info("No accessible completed Results jobs are available.")
        return
    opened_ref = f"{opened.get('task_group') or ''}/{opened.get('run_id') or ''}".strip("/")
    refs = [f"{row.get('task_group')}/{row.get('run_id')}" for row in jobs]
    labels = {
        f"{row.get('task_group')}/{row.get('run_id')}":
            f"{row.get('campaign_name') or row.get('display_tool') or row.get('tool') or 'Result'} · {row.get('job_code') or row.get('run_id')} · {row.get('status')}"
        for row in jobs
    }
    if opened_ref not in refs:
        opened_ref = refs[0]
    job_widget_key = "complex-analysis-results-job"
    if st.session_state.get(job_widget_key) not in refs:
        st.session_state[job_widget_key] = opened_ref
    selected_ref = st.selectbox("Completed Results job", refs, index=refs.index(opened_ref), format_func=lambda value: labels[value], key=job_widget_key)
    job = next(row for row in jobs if f"{row.get('task_group')}/{row.get('run_id')}" == selected_ref)
    options = _job_structure_options(user, job)
    if not options:
        st.info("This completed job has no accessible PDB/mmCIF or other supported structure artifact.")
        return
    relative_paths = [item["relative_path"] for item in options]
    requested_path = str(opened.get("relative_path") or st.session_state.get("complex-analysis-results-structure") or "")
    if requested_path not in relative_paths:
        requested_path = relative_paths[0]
    structure_widget_key = "complex-analysis-results-structure"
    if st.session_state.get(structure_widget_key) not in relative_paths:
        st.session_state[structure_widget_key] = requested_path
    structure_path = st.selectbox(
        "Structure / model / pose", relative_paths,
        index=relative_paths.index(requested_path),
        format_func=lambda value: next(
            (f"{item['candidate_id']} · {item['label']}" if item.get("candidate_id") else item["label"])
            for item in options if item["relative_path"] == value
        ),
        key=structure_widget_key,
    )
    selected = next(item for item in options if item["relative_path"] == structure_path)
    candidate = selected.get("candidate") or {}
    candidate_metadata = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
    st.caption(
        f"Results: {job.get('task_group')}/{job.get('run_id')} · "
        f"{job.get('display_tool') or job.get('tool') or 'Unknown tool'} · {selected['relative_path']}"
    )
    if candidate.get("candidate_id"):
        st.caption(f"Candidate/model: {candidate.get('candidate_id')} · {candidate.get('model_index') or candidate_metadata.get('model_index') or 'model not specified'}")
    try:
        content = selected["path"].read_bytes()
        components = component_inventory(content, selected["path"].name, candidate_metadata.get("components"))
    except (OSError, ValueError) as exc:
        st.error(f"Could not read the selected structure: {exc}")
        return
    _show_components(components)
    _show_structure(selected["path"], candidate_metadata, components,
                    key=f"complex-analysis-results-viewer-{selected_ref}-{structure_path}")
    _show_component_groups(components, f"{selected_ref}-{structure_path}")
    if candidate:
        st.markdown("#### Existing model evaluation metadata")
        st.json({
            "candidate_id": candidate.get("candidate_id"),
            "stage": candidate.get("stage"),
            "source_tool": candidate.get("source_tool") or job.get("tool"),
            "stored_metrics": candidate.get("metrics") or {},
            "model_confidence_outputs": candidate_metadata.get("confidence_outputs") or {},
            "source_metadata": {key: value for key, value in candidate_metadata.items()
                                if key not in {"components", "confidence_outputs"}},
        })
    else:
        st.info("No candidate evaluation metadata is associated with this structure artifact.")
    with st.expander("Results provenance", expanded=False):
        st.json({
            "task_group": job.get("task_group"),
            "run_id": job.get("run_id"),
            "job_code": job.get("job_code"),
            "tool": job.get("display_tool") or job.get("tool"),
            "project_id": job.get("project_id"),
            "completed_at": job.get("completed_at"),
            "candidate_id": candidate.get("candidate_id"),
            "model_index": candidate.get("model_index") or candidate_metadata.get("model_index"),
            "artifact_path": selected.get("artifact_path"),
            "job_relative_path": selected["relative_path"],
        })
    from mn_molecular_design_studio.molecular_database import MolecularDatabase
    _render_biotite_workspace(user, job, selected, components,
                              workflow_schema=workflow_schema,
                              project_selector=project_selector, submit=submit)
    _show_human_review(
        MolecularDatabase(), user, source_kind="results", subject_id=selected_ref,
        structure_ref=selected["relative_path"],
    )


def render_complex_analysis_workspace(user: dict, *, workflow_schema=None, project_selector=None, submit=None) -> None:
    """Render the shared Database/Results structure and provenance workspace."""
    opened = st.session_state.pop("complex-analysis-open", None)
    if opened:
        st.session_state["complex-analysis-active"] = opened
    active = st.session_state.get("complex-analysis-active") or {}
    if active.get("kind") and "complex-analysis-source-kind" not in st.session_state:
        st.session_state["complex-analysis-source-kind"] = "Database" if active["kind"] == "database" else "Results"
    st.title("Complex Analysis")
    st.caption("Inspect saved molecular structures, component identities, model metadata, and source provenance.")
    source = st.radio("Source", ["Database", "Results"], horizontal=True, key="complex-analysis-source-kind")
    try:
        if source == "Database":
            _render_database_source(user, active if active.get("kind") == "database" else {})
        else:
            _render_results_source(user, active if active.get("kind") == "results" else {},
                                   workflow_schema=workflow_schema,
                                   project_selector=project_selector, submit=submit)
    except PermissionError:
        st.warning("This record or job is unavailable or you do not have access to it.")
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        st.error(f"Could not load this structure source: {html.escape(str(exc))}")
