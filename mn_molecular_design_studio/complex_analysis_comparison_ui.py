"""Read-only collection comparison for saved Complex Analysis structures."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import streamlit as st

from mn_molecular_design_studio.complex_analysis_comparison import (
    compatibility_report,
    external_alignment_template,
    normalize_external_alignment_record,
    protein_chain_numbering,
    saved_target_contact_rows,
    target_structure_fingerprint,
    validate_supplied_alignment,
)


def _key(value: str) -> str:
    return hashlib.sha256(value.encode()).hexdigest()[:16]


def _read_json(path: Path) -> dict:
    try:
        value = json.loads(Path(path).read_text())
    except (OSError, ValueError):
        return {}
    return value if isinstance(value, dict) else {}


def _jobs_for_selected_run(jobs: list[dict], run_group: dict) -> list[dict]:
    """Intersect a selected logical run with the independently authorized local jobs."""
    sources = {(str(source.get("task_group") or ""), str(source.get("run_id") or ""))
               for source in run_group.get("run_sources") or [] if source.get("is_local")
               and int((source.get("source_summary") or {}).get("structure_count") or 0) > 0}
    return [job for job in jobs if (str(job.get("task_group") or ""), str(job.get("run_id") or "")) in sources]


def _comparison_slot(user: dict, slot: str, run_group: dict | None = None) -> dict | None:
    """Resolve one independently authorized Results or Database structure."""
    from mn_molecular_design_studio import complex_analysis_ui as workspace
    from mn_molecular_design_studio.complex_analysis_catalog import authorized_structure_catalog
    from mn_molecular_design_studio.molecular_database import MolecularDatabase

    if run_group:
        st.caption("Run: " + str(run_group.get("logical_run_name") or "Selected run"))
        source_kind = "Results"
    else:
        source_kind = st.selectbox("Structure source", ["Results", "Database"], key=f"complex-analysis-compare-{slot}-source")
    source: dict = {"kind": source_kind.casefold()}
    if source_kind == "Results":
        if run_group:
            from mn_molecular_design_studio.auth import can_access_job
            from mn_molecular_design_studio.job_catalog import get_job_access_record

            # A failed/partial generator can retain valid result structures.
            # Keep those poses available while rechecking ordinary job access.
            jobs = [source for source in run_group.get("run_sources") or []
                    if source.get("is_local") and source.get("task_group") and source.get("run_id")
                    and can_access_job(user, get_job_access_record(str(source["task_group"]), str(source["run_id"])))]
            jobs = _jobs_for_selected_run(jobs, run_group)
        else:
            jobs = workspace._accessible_result_jobs(user)
        if not jobs:
            st.info("No accessible Results jobs with local structures are available for this selected run. Remote files can be inspected on their source workstation." if run_group else
                    "No accessible completed Results jobs are available for this comparison slot.")
            return None
        # A BioEmu trajectory topology is not a coordinate-bearing pose, so
        # don't make it the first comparison choice when static predictions
        # are also available. Keep it selectable for workflows that expose a
        # saved static sample alongside the trajectory metadata.
        default_job_index = next(
            (index for index, item in enumerate(jobs)
             if str(item.get("display_tool") or item.get("tool") or "").casefold() != "bioemu2"),
            0,
        )
        refs = [f"{row.get('task_group')}/{row.get('run_id')}" for row in jobs]
        by_ref = {f"{row.get('task_group')}/{row.get('run_id')}": row for row in jobs}
        job_key = f"complex-analysis-compare-{slot}-job"
        if st.session_state.get(job_key) not in refs:
            st.session_state[job_key] = refs[default_job_index]
        selected_ref = st.selectbox(
            "Results source job" if run_group else "Completed Results job", refs,
            format_func=lambda ref: f"{by_ref[ref].get('display_tool') or by_ref[ref].get('tool') or 'Result'} · {by_ref[ref].get('job_code') or ref}",
            index=default_job_index,
            key=job_key,
        )
        job = by_ref[selected_ref]
        if run_group:
            from mn_molecular_design_studio.complex_analysis_dashboard import _discover_structures

            options = _discover_structures(user, job)
        else:
            options = workspace._job_structure_options(user, job)
        # BioEmu stores topology separately from trajectory coordinates. A
        # topology-only artifact is not a selected structural sample/pose.
        if str(job.get("display_tool") or job.get("tool") or "").casefold() == "bioemu2":
            options = [item for item in options if Path(str(item.get("relative_path") or "")).name.casefold() != "topology.cif"]
        if not options:
            st.info("This Results job has no accessible static structure model for comparison. Trajectory topology files without selected coordinates are not treated as poses.")
            return None
        option_keys = [str(item["relative_path"]) for item in options]
        by_path = {str(item["relative_path"]): item for item in options}
        relative = st.selectbox(
            "Structure / model / pose", option_keys,
            format_func=lambda path: f"{by_path[path].get('candidate_id') or 'Structure'} · {by_path[path].get('label') or path}",
            key=f"complex-analysis-compare-{slot}-results-structure",
        )
        selected = by_path[relative]
        candidate = selected.get("candidate") or {}
        metadata = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
        path = Path(selected["path"])
        source.update({
            "identity": selected_ref + ":" + relative,
            "label": f"{job.get('display_tool') or job.get('tool') or 'Results'} · {selected_ref} · {relative}",
            "path": path,
            "context": metadata,
            "record": {},
            "candidate": candidate,
            "job": job,
            "candidate_id": str(candidate.get("candidate_id") or ""),
            "model_id": str(candidate.get("model_index") or metadata.get("model_index") or path.name),
            "source_ref": selected_ref,
            "relative_path": relative,
            "project_id": str(job.get("project_id") or ""),
            "component_metadata": metadata.get("components") or [],
            "evaluation_entries": [],
            "structure_digest": workspace._structure_digest(path),
        })
        metrics = candidate.get("metrics") if isinstance(candidate.get("metrics"), dict) else {}
        if metrics:
            source["evaluation_entries"].append({
                "evaluation_id": str(candidate.get("candidate_id") or ""),
                "evaluator": str(candidate.get("source_tool") or job.get("display_tool") or job.get("tool") or "prediction model"),
                "run_id": str(job.get("run_id") or ""),
                "status": str(job.get("status") or "completed"),
                "metrics": metrics,
                "provenance": {"source": "Results candidate", "task_group": job.get("task_group"),
                               "run_id": job.get("run_id"), "candidate_id": candidate.get("candidate_id"),
                               "model_id": source["model_id"], "artifact": relative},
            })
        raw_target = metadata.get("target_id") or metadata.get("target_identity") or metadata.get("uniprot_id")
        source["target_identity"] = str(raw_target or "")
        source["coordinate_frame_id"] = str(metadata.get("coordinate_frame_id") or "")
    else:
        database = MolecularDatabase()
        query = st.text_input("Search Database", key=f"complex-analysis-compare-{slot}-query")
        catalog_key = f"complex-analysis-compare-catalog-{user.get('id', '')}-{slot}"
        cached_catalog = st.session_state.get(catalog_key)
        same_query = isinstance(cached_catalog, dict) and str(cached_catalog.get("query") or "") == query
        refresh_catalog = st.button("Refresh Database choices", key=f"complex-analysis-compare-{slot}-refresh")
        load_next = st.button("Load next Database page", key=f"complex-analysis-compare-{slot}-next",
                              disabled=not (same_query and cached_catalog.get("next_offset") is not None))
        catalog = authorized_structure_catalog(
            database, user, session_state=st.session_state, cache_key=catalog_key,
            query=query, limit=100, advance=load_next, refresh=refresh_catalog,
        )
        records = [item for item in catalog["records"] if item.get("representations") or item.get("structures")]
        st.caption(f"Showing {len(records):,} authorized choices on this page"
                   + ("; additional pages are available." if catalog.get("next_offset") is not None else "."))
        if not records:
            st.info("No accessible Database structures are available for this comparison slot.")
            return None
        record_ids = [str(item["id"]) for item in records]
        by_id = {str(item["id"]): item for item in records}
        molecule_key = f"complex-analysis-compare-{slot}-molecule"
        if st.session_state.get(molecule_key) not in record_ids:
            st.session_state[molecule_key] = record_ids[0]
        molecule_id = st.selectbox(
            "Database complex", record_ids,
            format_func=lambda value: f"{by_id[value].get('name') or 'Complex'} · {value}",
            key=molecule_key,
        )
        rows = database.list_molecules(user, molecule_ids=[molecule_id], collection="all", limit=1)["items"]
        if not rows:
            st.warning("This Database record is unavailable or you do not have access to it.")
            return None
        record = rows[0]
        structures = list(record.get("representations") or record.get("structures") or [])
        if not structures:
            st.info("This Database record has no accessible structure representations.")
            return None
        structure_ids = [str(item["id"]) for item in structures]
        by_structure = {str(item["id"]): item for item in structures}
        structure_id = st.selectbox(
            "Structure / model / pose", structure_ids,
            format_func=lambda value: workspace._representation_label(by_structure[value]),
            key=f"complex-analysis-compare-{slot}-database-structure",
        )
        structure = by_structure[structure_id]
        try:
            path = database.structure_file(user, structure_id)
        except (OSError, PermissionError, ValueError) as exc:
            st.warning(f"The selected Database structure is unavailable: {exc}")
            return None
        context = structure.get("context") if isinstance(structure.get("context"), dict) else {}
        entries = []
        evaluation_rows, by_evaluation = workspace._evaluation_options(record)
        for observation in record.get("observations", []):
            descriptor = next((desc for known, desc in by_evaluation.values() if known.get("id") == observation.get("id")), None)
            if not descriptor or not workspace._evaluation_matches_structure(descriptor, structure):
                continue
            entries.append({
                "evaluation_id": descriptor.get("id", ""),
                "evaluator": descriptor.get("evaluator", "Unknown evaluator"),
                "run_id": descriptor.get("run_id", ""),
                "status": descriptor.get("status", "completed"),
                "metrics": observation.get("metrics") or {},
                "provenance": {"evaluation": descriptor,
                               "context": observation.get("context") or {},
                               "score_sources": (observation.get("normalized_output") or {}).get("score_sources") or {}},
            })
        source.update({
            "identity": f"{molecule_id}:{structure_id}",
            "label": f"Database · {record.get('name') or molecule_id} · {structure_id}",
            "path": path,
            "context": context,
            "record": record,
            "structure": structure,
            "candidate": {},
            "candidate_id": str(context.get("candidate_id") or ""),
            "model_id": str(context.get("model_id") or context.get("model_index") or context.get("model") or structure.get("model_index") or structure_id),
            "source_ref": str(record.get("id") or molecule_id),
            "component_metadata": record.get("components") or [],
            "evaluation_entries": entries,
            "structure_digest": workspace._structure_digest(path),
            "target_identity": str(context.get("target_id") or context.get("target_identity") or context.get("target_molecule_id") or ""),
            "coordinate_frame_id": str(context.get("coordinate_frame_id") or (
                f"receptor-sha256:{context.get('target_structure_digest')}" if context.get("target_structure_digest") else "")),
        })
        try:
            source["project_id"] = str(database.group(str(structure.get("group_id") or ""), user).get("project_id") or "")
        except (KeyError, PermissionError, ValueError):
            source["project_id"] = ""

    try:
        content = source["path"].read_bytes()
        components = workspace.component_inventory(content, source["path"].name, source.get("component_metadata"))
    except (OSError, ValueError) as exc:
        st.warning(f"Could not read this structure: {exc}")
        return None
    protein_chains = workspace._protein_chain_ids(components)
    if not protein_chains:
        st.info("This structure has no identified protein chains. It remains available in its source page; comparison requires a protein component.")
        return None
    component_types = {str(item.get("type") or "other").casefold() for item in components}
    is_ppi = component_types == {"protein"}
    is_protein_ligand = {"protein", "small_molecule"} <= component_types and not {"dna", "rna"} & component_types
    if not (is_ppi or is_protein_ligand):
        kind = "protein–DNA/RNA" if {"protein", "dna"} <= component_types or {"protein", "rna"} <= component_types else "mixed-component"
        st.info(f"This comparison supports protein–protein and protein–ligand structures. This representation also contains {kind} components; single-structure viewing remains available.")
        return None
    chain_key = _key(source["identity"])
    context = source.get("context") if isinstance(source.get("context"), dict) else {}
    raw = source.get("candidate", {}).get("raw_metadata") if isinstance(source.get("candidate", {}).get("raw_metadata"), dict) else {}
    default_target = [str(value) for value in (source.get("candidate", {}).get("target_chains") or raw.get("target_chains") or context.get("target_chains") or []) if str(value) in protein_chains]
    if is_ppi:
        default_binder = [str(value) for value in (source.get("candidate", {}).get("binder_chains") or raw.get("binder_chains") or context.get("binder_chains") or []) if str(value) in protein_chains]
        if not default_target:
            default_target = protein_chains[:1]
        if not default_binder:
            default_binder = [chain for chain in protein_chains if chain not in default_target]
        a, b = st.columns(2)
        target_chains = a.multiselect("Target chain group", protein_chains, default=default_target,
                                      key=f"complex-analysis-compare-{slot}-target-{chain_key}")
        binder_chains = b.multiselect("Protein partner chain group", protein_chains, default=default_binder,
                                      key=f"complex-analysis-compare-{slot}-partner-{chain_key}")
        if not binder_chains or set(target_chains) & set(binder_chains):
            st.info("Select a nonempty protein-partner group and keep it separate from the target group.")
            return None
    else:
        if not default_target:
            default_target = list(protein_chains)
        target_chains = st.multiselect("Receptor protein chain group", protein_chains, default=default_target,
                                       key=f"complex-analysis-compare-{slot}-target-{chain_key}")
        binder_chains = []
        ligands = [item for item in components if str(item.get("type") or "").casefold() == "small_molecule"]
        if not ligands:
            st.info("No ligand residue was identified in this representation.")
            return None
        ligand_key = f"complex-analysis-compare-{slot}-ligand-{chain_key}"
        ligand_ids = [str(item.get("id") or "") for item in ligands]
        ligand_by_id = {str(item.get("id") or ""): item for item in ligands}
        if st.session_state.get(ligand_key) not in ligand_ids:
            st.session_state[ligand_key] = ligand_ids[0]
        selected_ligand_id = st.selectbox(
            "Ligand residue", ligand_ids,
            format_func=lambda value: f"{ligand_by_id[value].get('label') or ligand_by_id[value].get('ccd') or 'Ligand'} · {ligand_by_id[value].get('residue_key') or ligand_by_id[value].get('source_chain')}",
            key=ligand_key,
        )
        from mn_molecular_design_studio.complex_analysis_ui import _ligand_selection_for_component
        ligand_selection = _ligand_selection_for_component(ligand_by_id[selected_ligand_id], ligands, {
            "components": source.get("component_metadata") or [],
            "ligand_selection": context.get("ligand_selection"),
            "ligand_smiles": context.get("ligand_smiles"),
        })
    try:
        target_numbering = protein_chain_numbering(source["path"], list(target_chains)) if target_chains else []
        binder_numbering = protein_chain_numbering(source["path"], list(binder_chains)) if is_ppi else []
    except ValueError as exc:
        st.info(str(exc))
        return None
    source.update({
        "components": components,
        "target_chains": list(target_chains),
        "binder_chains": list(binder_chains),
        "target_numbering": target_numbering,
        "target_chains_sequence": [row["sequence"] for row in target_numbering],
        "binder_chains_sequence": [row["sequence"] for row in binder_numbering],
        "complex_type": "protein-protein" if is_ppi else "protein-ligand",
        "format": source["path"].suffix.casefold(),
        "viewer_key": f"complex-analysis-compare-viewer-{slot}-{chain_key}",
    })
    if is_protein_ligand:
        source["ligand_selection"] = ligand_selection
    partner = ", ".join(binder_chains) if is_ppi else f"{ligand_selection.get('ccd') or 'Ligand'} {ligand_selection.get('residue_key') or ''}".strip()
    st.caption(f"{source['label']} · model {source['model_id'] or 'not specified'} · target chains {', '.join(target_chains)} · partner {partner}")
    return source


def _render_saved_evaluations(sources: list[dict]) -> None:
    entries = []
    for index, source in enumerate(sources):
        sample_label = f"Structure {chr(ord('A') + index)}"
        for entry in source.get("evaluation_entries", []):
            entries.append((sample_label, entry))
    st.markdown("#### Saved evaluations and measurements")
    if not entries:
        st.info("No saved model measurements or structure-associated evaluations were found for these selected representations.")
        return
    evaluators = sorted({str(entry.get("evaluator") or "Unknown evaluator") for _, entry in entries})
    selected = st.multiselect("Filter by evaluator", evaluators, default=evaluators, key="complex-analysis-compare-evaluator-filter")
    rows = []
    for sample_label, entry in entries:
        if entry.get("evaluator") not in selected:
            continue
        rows.append({
            "Structure": sample_label,
            "Evaluator": entry.get("evaluator"),
            "Status": entry.get("status"),
            "Run": entry.get("run_id"),
            "Evaluation / candidate ID": entry.get("evaluation_id"),
            "Stored measurements": json.dumps(entry.get("metrics") or {}, sort_keys=True, default=str),
            "Measurement provenance": json.dumps(entry.get("provenance") or {}, sort_keys=True, default=str),
        })
    if rows:
        st.dataframe(rows, hide_index=True, use_container_width=True)
    else:
        st.info("No evaluation matches the selected evaluator filter.")
    st.caption("Rows preserve each saved evaluation's own method, context, and values. This view does not normalize measurements into a common score or rank structures.")


def _saved_biotite_contacts(user: dict, source: dict) -> dict | None:
    from mn_molecular_design_studio import complex_analysis_ui as workspace

    # The unified cohort may have already resolved a source-qualified saved contact
    # artifact (including an authorized remote Results artifact). Reuse that exact
    # receipt instead of repeating discovery or fetching it again.
    saved = source.get("saved_biotite_contacts")
    if isinstance(saved, dict) and isinstance(saved.get("contacts"), list):
        return saved

    if source["kind"] == "results":
        job = source.get("job") or {}
        matches = workspace._matching_biotite_jobs(
            user, source_task_group=str(job.get("task_group") or ""), source_run_id=str(job.get("run_id") or ""),
            candidate_id=str(source.get("candidate_id") or ""), structure_digest=str(source.get("structure_digest") or ""),
            binder_chains=list(source["binder_chains"]), target_chains=list(source["target_chains"]),
            contact_distance=5.0, model_id=str(source.get("model_id") or ""),
        )
    else:
        structure = source.get("structure") or {}
        record = source.get("record") or {}
        group_id = str(structure.get("group_id") or "")
        source_selection = {
            "molecule_id": str(record.get("id") or ""), "structure_id": str(structure.get("id") or ""),
            "group_id": group_id, "model_id": str(source.get("model_id") or ""),
            "evaluation_id": "", "structure_digest": source.get("structure_digest"),
        }
        try:
            from mn_molecular_design_studio.molecular_database import MolecularDatabase
            group = MolecularDatabase().group(group_id, user)
            source_selection["project_id"] = str(group.get("project_id") or "")
        except (KeyError, PermissionError, ValueError):
            return None
        candidate_id = "database-" + hashlib.sha256((source_selection["molecule_id"] + ":" + source_selection["structure_id"]).encode()).hexdigest()[:24]
        matches = workspace._matching_database_analysis_jobs(
            user, workflow_id="analysis.biotite_interface", source=source_selection,
            candidate_id=candidate_id, model_id=str(source.get("model_id") or ""),
            binder_chains=list(source["binder_chains"]), target_chains=list(source["target_chains"]),
            parameters={"contact_distance": 5.0},
        )
    completed = next((item for item in matches if str(item.get("status") or "").casefold() in {"completed", "succeeded"}), None)
    if not completed:
        return None
    item = completed.get("selected_item") or {}
    contacts_path = Path(completed["run_dir"]) / "artifacts" / "analysis" / str(item.get("id") or "") / "contacts.json"
    contacts = workspace._read_json(contacts_path)
    if not contacts:
        return None
    return {"evaluation": completed, "contacts": contacts.get("residue_pairs") or [], "selected_item": item}


def _render_supplied_alignment_overlay(sources: list[dict]) -> None:
    """Display two saved structures using a user-supplied transform and map."""
    if len(sources) < 2:
        return
    from mn_molecular_design_studio import complex_analysis_ui as workspace

    st.markdown("#### Overlay with supplied transform")
    st.caption(
        "Paste an externally supplied moving-to-reference 4×4 transform and explicit target-residue map. "
        "Studio validates the matrix form, selected model identities, and residue references; it does not fit, "
        "estimate, or score the alignment. The saved source structures are never modified."
    )
    labels = [f"{chr(ord('A') + index)} · {item['label']} · model {item.get('model_id') or 'unspecified'}" for index, item in enumerate(sources)]
    slots = list(range(len(sources)))
    first, second = st.columns(2)
    reference_index = first.selectbox(
        "Reference structure", slots, format_func=lambda index: labels[index],
        key="complex-analysis-supplied-transform-reference",
    )
    moving_index = second.selectbox(
        "Moving structure", slots, index=1, format_func=lambda index: labels[index],
        key="complex-analysis-supplied-transform-moving",
    )
    if reference_index == moving_index:
        st.info("Choose two different saved structure representations.")
        return
    reference = sources[reference_index]
    moving = sources[moving_index]
    form_key = _key("|".join((
        str(reference.get("identity")), str(reference.get("model_id")), str(reference.get("structure_digest")),
        ",".join(reference.get("target_chains") or []),
        str(moving.get("identity")), str(moving.get("model_id")), str(moving.get("structure_digest")),
        ",".join(moving.get("target_chains") or []),
    )))
    input_mode = st.radio(
        "Transformation input", ["Manual fields", "External JSON record"], horizontal=True,
        key=f"complex-analysis-transform-input-mode-{form_key}",
    )
    matrix_text = mapping_text = transform_source = transform_method = mapping_note = record_json = ""
    if input_mode == "External JSON record":
        template = external_alignment_template(reference, moving)
        st.caption(
            "Import schema: mnds.external-alignment/v1. The record must identify these exact structures, models, "
            "and file digests; residue correspondence and source provenance remain explicit. The downloaded "
            "identity matrix is only a placeholder and must be replaced with the external transform."
        )
        st.download_button(
            "Download bound record template", data=json.dumps(template, indent=2),
            file_name="mnds-external-alignment-v1.json", mime="application/json",
            key=f"complex-analysis-transform-template-{form_key}",
        )
        uploaded = st.file_uploader(
            "Upload external alignment record", type=["json"],
            key=f"complex-analysis-transform-upload-{form_key}",
        )
        pasted_json = st.text_area(
            "Or paste external alignment record (JSON)", placeholder=json.dumps(template, indent=2),
            key=f"complex-analysis-transform-record-{form_key}", height=300,
        )
        try:
            uploaded_json = uploaded.getvalue().decode("utf-8") if uploaded is not None else ""
        except UnicodeDecodeError:
            uploaded_json = "\0invalid-utf8"
        record_json = uploaded_json or pasted_json
        if uploaded_json:
            st.caption(f"Using uploaded record: {uploaded.name}")
    else:
        matrix_text = st.text_area(
            "Moving-to-reference 4×4 transform (JSON rows)",
            placeholder='[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]',
            key=f"complex-analysis-transform-matrix-{form_key}", height=90,
        )
        mapping_text = st.text_area(
            "Explicit target-residue correspondence (JSON array)",
            placeholder='[{"reference":{"chain":"A","number":10,"insertion":"","residue_name":"ALA"},'
                        '"moving":{"chain":"X","number":10,"insertion":"","residue_name":"ALA"}}]',
            key=f"complex-analysis-transform-map-{form_key}", height=130,
        )
        source_col, method_col = st.columns(2)
        transform_source = source_col.text_input(
            "Transformation source", placeholder="Record, file, or external tool and version",
            key=f"complex-analysis-transform-source-{form_key}",
        )
        transform_method = method_col.text_input(
            "Transformation method", placeholder="As reported by its source",
            key=f"complex-analysis-transform-method-{form_key}",
        )
        mapping_note = st.text_input(
            "Mapping uncertainty / provenance note (optional)",
            key=f"complex-analysis-transform-note-{form_key}",
        )
    current_input_digest = _key("\0".join((input_mode, record_json, matrix_text, mapping_text, transform_source, transform_method, mapping_note,
                                           str(reference.get("identity")), str(reference.get("model_id")),
                                           str(reference.get("structure_digest")), ",".join(reference.get("target_chains") or []),
                                           str(moving.get("identity")), str(moving.get("model_id")),
                                           str(moving.get("structure_digest")), ",".join(moving.get("target_chains") or []))))
    state_key = f"complex-analysis-supplied-transform-result-{form_key}"
    if st.button("Validate and show overlay", key=f"complex-analysis-transform-apply-{form_key}"):
        try:
            if input_mode == "External JSON record":
                record = normalize_external_alignment_record(json.loads(record_json), reference, moving)
            else:
                matrix = json.loads(matrix_text)
                correspondence = json.loads(mapping_text)
                record = {
                    "reference_identity": str(reference.get("identity") or ""),
                    "reference_model_id": str(reference.get("model_id") or ""),
                    "reference_digest": str(reference.get("structure_digest") or ""),
                    "moving_identity": str(moving.get("identity") or ""),
                    "moving_model_id": str(moving.get("model_id") or ""),
                    "moving_digest": str(moving.get("structure_digest") or ""),
                    "matrix_4x4": matrix,
                    "residue_correspondence": correspondence,
                    "source": transform_source,
                    "method": transform_method,
                    "mapping_note": mapping_note,
                }
            report = validate_supplied_alignment(reference, moving, record)
            st.session_state[state_key] = {"digest": current_input_digest, "record": record, "report": report}
        except (json.JSONDecodeError, TypeError, ValueError) as exc:
            st.session_state[state_key] = {
                "digest": current_input_digest,
                "record": {},
                "report": {"valid": False, "errors": [f"Could not read the supplied JSON: {exc}"], "warnings": []},
            }

    saved = st.session_state.get(state_key) or {}
    if saved.get("digest") != current_input_digest:
        return
    report = saved.get("report") or {}
    for message in report.get("errors", []):
        st.error(message)
    for message in report.get("warnings", []):
        st.warning(message)
    if not report.get("valid"):
        return

    st.success(
        f"Rigid-transform form and residue references validated · {report['correspondence_count']} explicit residue pairs · "
        f"reference unmapped: {report['reference_unmapped_residues']} · moving unmapped: {report['moving_unmapped_residues']}"
    )
    st.caption(
        f"Source: {report['source']} · method: {report['method']} · "
        f"reference: {reference['identity']} / model {reference['model_id']} / {str(reference.get('structure_digest') or '')[:16]} · "
        f"moving: {moving['identity']} / model {moving['model_id']} / {str(moving.get('structure_digest') or '')[:16]}"
    )
    if report.get("mapping_note"):
        st.caption(f"User-supplied mapping note: {report['mapping_note']}")
    external_provenance = report.get("external_provenance") or {}
    if external_provenance:
        external_source = external_provenance.get("source") or {}
        st.caption(
            f"External record: {external_source.get('tool')} {external_source.get('version')} · "
            f"record ID: {external_source.get('record_id') or 'not supplied'} · "
            "provenance is supplied by the user and is not independently authenticated."
        )
        with st.expander("Full external source metadata"):
            st.json(external_provenance)
    from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer import molstar_custom_component
    from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer.dataclasses import StructureVisualization

    visualizations = [
        StructureVisualization(
            pdb=workspace._viewer_text(reference["path"], reference.get("context") or {}),
            color="chain-id", representation_type="cartoon+ball-and-stick",
        ),
        StructureVisualization(
            pdb=workspace._viewer_text(moving["path"], moving.get("context") or {}),
            color="chain-id", representation_type="cartoon+ball-and-stick",
            transform_matrix=report["matrix_4x4"],
        ),
    ]
    molstar_custom_component(
        visualizations, key=f"complex-analysis-supplied-transform-viewer-{form_key}",
        height=650, show_controls=True, show_structure_tools=False, force_reload=True,
    )


def _saved_structure_selection(sources: list[dict], *, workflow_id: str) -> dict:
    schemas = {
        "analysis.interface_clustering": "studio.interface_clustering_selection.v1",
        "analysis.foldseek": "studio.foldseek_selection.v1",
    }
    try:
        schema = schemas[workflow_id]
    except KeyError as exc:
        raise ValueError(f"Unsupported saved-structure collection workflow: {workflow_id}") from exc
    selected = []
    for sample in sources:
        item = {
            "kind": str(sample.get("kind") or "").casefold(),
            "structure_digest": str(sample.get("structure_digest") or ""),
            "model_id": str(sample.get("model_id") or ""),
            "target_chains": list(sample.get("target_chains") or []),
            "binder_chains": list(sample.get("binder_chains") or []),
            "evaluation_ids": [str(entry.get("evaluation_id") or "")
                               for entry in sample.get("evaluation_entries", [])
                               if str(entry.get("evaluation_id") or "")],
            "label": str(sample.get("label") or sample.get("identity") or "Saved structure"),
        }
        if item["kind"] == "results":
            job = sample.get("job") or {}
            item.update({
                "task_group": str(job.get("task_group") or ""),
                "run_id": str(job.get("run_id") or ""),
                "structure_ref": str(sample.get("relative_path") or ""),
                "candidate_id": str(sample.get("candidate_id") or ""),
            })
        elif item["kind"] == "database":
            record = sample.get("record") or {}
            structure = sample.get("structure") or {}
            item.update({
                "molecule_id": str(record.get("id") or ""),
                "structure_id": str(structure.get("id") or ""),
                "group_id": str(structure.get("group_id") or ""),
            })
        selected.append(item)
    return {"schema": schema, "sources": selected}


def _interface_clustering_selection(sources: list[dict]) -> dict:
    """Keep the original Interface Clustering selector API for callers/tests."""
    return _saved_structure_selection(sources, workflow_id="analysis.interface_clustering")


def _foldseek_selection(sources: list[dict]) -> dict:
    return _saved_structure_selection(sources, workflow_id="analysis.foldseek")


def _matching_interface_clustering_jobs(user: dict, *, selection_key: str,
                                        similarity_threshold: float,
                                        contact_distance: float) -> list[dict]:
    """Find accessible batch results for this exact selected structure cohort."""
    return _matching_collection_jobs(
        user, tool="interface_clustering", selection_key=selection_key,
        settings={"similarity_threshold": similarity_threshold, "contact_distance": contact_distance},
    )


def _matching_foldseek_jobs(user: dict, *, selection_key: str,
                            similarity_threshold: float, coverage: float) -> list[dict]:
    """Find accessible Foldseek results for this exact cohort and settings."""
    return _matching_collection_jobs(
        user, tool="foldseek", selection_key=selection_key,
        settings={"similarity_threshold": similarity_threshold, "coverage": coverage},
    )


def _matching_collection_jobs(user: dict, *, tool: str, selection_key: str,
                              settings: dict[str, float]) -> list[dict]:
    """Find accessible cohort results by the exact saved selection and settings."""
    from mn_molecular_design_studio.auth import accessible_project_ids, can_access_job
    from mn_molecular_design_studio.job_catalog import get_job_access_record, search_jobs

    rows = search_jobs(
        tools=[tool], task_group="analysis", limit=1000,
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
        collection = _read_json(run_dir / "artifacts" / "analysis" / "complex_analysis_collection.json")
        manifest = _read_json(run_dir / "artifacts" / "analysis" / "batch_manifest.json")
        if collection.get("selection_key") != selection_key or manifest.get("tool") != tool:
            continue
        saved_settings = manifest.get("settings") if isinstance(manifest.get("settings"), dict) else {}
        defaults = {"similarity_threshold": 0.5, "contact_distance": 5.0, "coverage": 0.8}
        try:
            if any(float(saved_settings.get(name, defaults.get(name, 0.0))) != float(value)
                   for name, value in settings.items()):
                continue
        except (TypeError, ValueError):
            continue
        metadata = _read_json(run_dir / "metadata.json")
        batch = _read_json(run_dir / "artifacts" / "analysis" / "batch_status.json")
        status = str(row.get("status") or metadata.get("status") or "unknown").casefold()
        matches.append({**row, "status": status, "run_dir": run_dir,
                        "manifest": manifest, "collection": collection,
                        "metadata": metadata, "batch_status": batch})
    rank = {"completed": 0, "succeeded": 0, "running": 1, "starting": 1, "queued": 1,
            "preparing": 1, "pending": 1, "failed": 2, "cancelled": 3, "canceled": 3}
    return sorted(matches, key=lambda item: (rank.get(item["status"], 4), str(item.get("created_at") or "")))


def _render_interface_clustering_results(user: dict, evaluation: dict, *, key: str) -> None:
    run_dir = Path(evaluation["run_dir"])
    status = str(evaluation.get("status") or "unknown").casefold()
    st.caption(f"Interface Clustering · {status.title()} · {evaluation.get('job_code') or evaluation.get('run_id')}")

    @st.fragment(run_every=5 if status not in {"completed", "succeeded", "failed", "cancelled", "canceled"} else None)
    def refresh_status():
        latest = str(_read_json(run_dir / "metadata.json").get("status") or status).casefold()
        st.metric("Clustering job status", latest.title())
        if status not in {"completed", "succeeded", "failed", "cancelled", "canceled"} and latest in {"completed", "succeeded", "failed", "cancelled", "canceled"}:
            st.rerun(scope="app")

    refresh_status()
    from mn_molecular_design_studio.quality_ui import render_progress
    render_progress(run_dir, status, key + "-progress")
    if status in {"failed", "cancelled", "canceled"}:
        st.error(evaluation.get("error") or "Interface Clustering did not complete. Review the existing Results logs before retrying.")
    collection = evaluation.get("collection") or {}
    if status in {"completed", "succeeded"}:
        metrics_path = run_dir / "artifacts" / "analysis" / "qc_metrics.csv"
        if metrics_path.is_file():
            import pandas as pd
            try:
                frame = pd.read_csv(metrics_path)
            except (OSError, ValueError, pd.errors.ParserError):
                frame = pd.DataFrame()
            if not frame.empty:
                items = collection.get("items") if isinstance(collection.get("items"), list) else []
                by_candidate = {str(item.get("candidate_id") or ""): item for item in items}
                source_column = "source_candidate_id" if "source_candidate_id" in frame else "analysis_source_candidate_id"
                if source_column in frame:
                    frame["Source structure"] = frame[source_column].astype(str).map(
                        lambda candidate_id: by_candidate.get(candidate_id, {}).get("source_identity", candidate_id)
                    )
                    from mn_molecular_design_studio import complex_analysis_ui as workspace
                    frame["Saved Evaluation"] = frame[source_column].astype(str).map(
                        lambda candidate_id: workspace._analysis_evaluation_id(
                            run_dir, user, candidate_id, evaluator="interface_clustering",
                        )
                    )
                st.markdown("##### Cluster assignments saved by this evaluation")
                st.dataframe(frame, hide_index=True, use_container_width=True)
    with st.expander("Interface Clustering provenance", expanded=False):
        st.json({
            "workflow": "analysis.interface_clustering",
            "analysis_job": f"{evaluation.get('task_group')}/{evaluation.get('run_id')}",
            "collection_id": collection.get("collection_id"),
            "selection_key": collection.get("selection_key"),
            "settings": (evaluation.get("manifest") or {}).get("settings") or {},
            "sources": collection.get("items") or [],
            "method": "Single-linkage clusters of contacted target-residue sets using Jaccard similarity.",
        })


def _render_interface_clustering(user: dict, sources: list[dict], *, compatible: bool,
                                 workflow_schema=None, submit=None) -> None:
    st.markdown("#### Interface Clustering")
    st.caption("Runs the registered workflow on these saved structures: single-linkage clustering by Jaccard similarity of contacted target residues. The workflow supports up to 2,000 items and requires one identical target structure and author numbering throughout.")
    selection = _interface_clustering_selection(sources)
    fingerprints = []
    fingerprint_errors = []
    for sample in sources:
        try:
            fingerprints.append(target_structure_fingerprint(sample["path"], list(sample["target_chains"])))
        except (KeyError, OSError, ValueError, ImportError) as exc:
            fingerprints.append("")
            fingerprint_errors.append(f"{sample.get('label') or 'Selected structure'}: {exc}")
    same_target = bool(fingerprints) and all(fingerprints) and len(set(fingerprints)) == 1
    if fingerprint_errors:
        st.info("Interface Clustering cannot verify the exact selected target structure. " + "; ".join(fingerprint_errors))
    elif compatible and not same_target:
        st.info("The target sequences and numbering are compatible for side-by-side review, but Interface Clustering also requires identical target coordinates and author chain IDs.")

    project_ids = {str(sample.get("project_id") or "") for sample in sources}
    only_private_database = project_ids == {""} and all(sample.get("kind") == "database" for sample in sources)
    same_project = len(project_ids) == 1 and (bool(next(iter(project_ids))) or only_private_database)
    project_id = next(iter(project_ids)) if same_project else ""
    if not same_project:
        st.info("The selected structures do not share one project access boundary. Choose structures from the same project before launching a shared evaluation.")
    elif project_id:
        st.caption(f"Analysis access: source project {project_id}.")
    else:
        st.caption("Analysis access: private to your account, matching the selected private Database records.")

    schema = None
    schema_error = ""
    try:
        schema = workflow_schema("analysis.interface_clustering") if workflow_schema else None
        if not schema:
            schema_error = "The registered Interface Clustering workflow schema is unavailable."
    except Exception as exc:
        schema_error = f"The registered Interface Clustering workflow could not be loaded: {exc}"
    if schema_error:
        st.info(schema_error)
    parameter_specs = {str(item.get("name")): item for item in (schema or {}).get("parameters", [])}
    first, second = st.columns(2)
    similarity_threshold = float(first.number_input(
        "Jaccard similarity threshold", min_value=0.01, max_value=1.0,
        value=float(parameter_specs.get("similarity_threshold", {}).get("default", 0.5)), step=0.05,
        key="complex-analysis-interface-clustering-threshold",
        help="The registered workflow uses single-linkage clustering of contacted target residues.",
    ))
    contact_distance = float(second.number_input(
        "Contact distance (Å)", min_value=0.1, max_value=20.0,
        value=float(parameter_specs.get("contact_distance", {}).get("default", 5.0)), step=0.1,
        key="complex-analysis-interface-clustering-distance",
    ))
    selection_key = ""
    try:
        from mn_molecular_design_studio.complex_analysis_sources import interface_clustering_selection_key
        selection_key = interface_clustering_selection_key(selection, project_id)
    except (TypeError, ValueError):
        pass
    existing = _matching_interface_clustering_jobs(
        user, selection_key=selection_key, similarity_threshold=similarity_threshold,
        contact_distance=contact_distance,
    ) if selection_key and compatible and same_target and same_project and schema else []
    current = next((item for item in existing if item["status"] in {"completed", "succeeded"}), None)
    if current is None:
        current = next((item for item in existing if item["status"] not in {"failed", "cancelled", "canceled"}), None)
    if current:
        _render_interface_clustering_results(user, current, key=selection_key[:14])
    elif existing:
        failed = existing[0]
        st.error(f"Previous matching Interface Clustering job failed: {failed.get('error') or 'No per-item error was recorded.'}")
        _render_interface_clustering_results(user, failed, key=selection_key[:14])

    launch_ready = (compatible and same_target and same_project and bool(schema)
                    and bool(selection_key) and submit is not None)
    if current and current["status"] in {"completed", "succeeded"}:
        return
    button_label = "Retry Interface Clustering" if existing else "Run Interface Clustering"
    if st.button(button_label, type="primary", disabled=not launch_ready,
                 key=f"complex-analysis-interface-clustering-run-{selection_key[:14]}"):
        parameters = {
            "similarity_threshold": round(similarity_threshold, 4),
            "contact_distance": round(contact_distance, 3),
            "interface_clustering_collection": selection,
        }
        try:
            submit(schema, "analysis.interface_clustering", parameters, project_id, {}, "Interface Clustering")
        except (OSError, RuntimeError, TypeError, ValueError) as exc:
            st.error(f"Could not start Interface Clustering: {exc}")
            return
        st.rerun()


def _render_foldseek_results(user: dict, evaluation: dict, *, key: str) -> None:
    run_dir = Path(evaluation["run_dir"])
    status = str(evaluation.get("status") or "unknown").casefold()
    st.caption(f"Foldseek · {status.title()} · {evaluation.get('job_code') or evaluation.get('run_id')}")

    @st.fragment(run_every=5 if status not in {"completed", "succeeded", "failed", "cancelled", "canceled"} else None)
    def refresh_status():
        latest = str(_read_json(run_dir / "metadata.json").get("status") or status).casefold()
        st.metric("Foldseek job status", latest.title())
        if status not in {"completed", "succeeded", "failed", "cancelled", "canceled"} and latest in {"completed", "succeeded", "failed", "cancelled", "canceled"}:
            st.rerun(scope="app")

    refresh_status()
    from mn_molecular_design_studio.quality_ui import render_progress
    render_progress(run_dir, status, key + "-progress")
    if status in {"failed", "cancelled", "canceled"}:
        st.error(evaluation.get("error") or "Foldseek did not complete. Review the existing Results logs before retrying.")
    collection = evaluation.get("collection") or {}
    if status in {"completed", "succeeded"}:
        metrics_path = run_dir / "artifacts" / "analysis" / "qc_metrics.csv"
        if metrics_path.is_file():
            import pandas as pd
            try:
                frame = pd.read_csv(metrics_path)
            except (OSError, ValueError, pd.errors.ParserError):
                frame = pd.DataFrame()
            if not frame.empty:
                items = collection.get("items") if isinstance(collection.get("items"), list) else []
                by_candidate = {str(item.get("candidate_id") or ""): item for item in items}
                source_column = "source_candidate_id" if "source_candidate_id" in frame else "analysis_source_candidate_id"
                if source_column in frame:
                    frame["Source structure"] = frame[source_column].astype(str).map(
                        lambda candidate_id: by_candidate.get(candidate_id, {}).get("source_identity", candidate_id)
                    )
                    from mn_molecular_design_studio import complex_analysis_ui as workspace
                    frame["Saved Evaluation"] = frame[source_column].astype(str).map(
                        lambda candidate_id: workspace._analysis_evaluation_id(
                            run_dir, user, candidate_id, evaluator="foldseek",
                        )
                    )
                st.markdown("##### Foldseek cluster assignments saved by this evaluation")
                st.dataframe(frame, hide_index=True, use_container_width=True)
    with st.expander("Foldseek provenance", expanded=False):
        st.json({
            "workflow": "analysis.foldseek",
            "analysis_job": f"{evaluation.get('task_group')}/{evaluation.get('run_id')}",
            "collection_id": collection.get("collection_id"),
            "selection_key": collection.get("selection_key"),
            "settings": (evaluation.get("manifest") or {}).get("settings") or {},
            "sources": collection.get("items") or [],
            "method": "Foldseek easy-cluster using TM-score threshold and coverage on each selected protein chain.",
        })


def _render_foldseek(user: dict, sources: list[dict], *, workflow_schema=None, submit=None) -> None:
    st.markdown("#### Foldseek")
    st.caption("Runs the registered Foldseek cohort workflow on the selected protein chains. Each input must resolve to exactly one protein chain; target compatibility, shared target coordinates, and common author numbering are not required.")
    selection = _foldseek_selection(sources)
    valid_chains = all(len(sample.get("binder_chains") or []) == 1 for sample in sources)
    if not valid_chains:
        st.info("Foldseek requires exactly one selected protein chain per structure. Select one partner chain in every slot.")

    project_ids = {str(sample.get("project_id") or "") for sample in sources}
    only_private_database = project_ids == {""} and all(sample.get("kind") == "database" for sample in sources)
    same_project = len(project_ids) == 1 and (bool(next(iter(project_ids))) or only_private_database)
    project_id = next(iter(project_ids)) if same_project else ""
    if not same_project:
        st.info("The selected structures do not share one project access boundary. Choose structures from the same project before launching a shared evaluation.")
    elif project_id:
        st.caption(f"Analysis access: source project {project_id}.")
    else:
        st.caption("Analysis access: private to your account, matching the selected private Database records.")

    schema = None
    schema_error = ""
    try:
        schema = workflow_schema("analysis.foldseek") if workflow_schema else None
        if not schema:
            schema_error = "The registered Foldseek workflow schema is unavailable."
    except Exception as exc:
        schema_error = f"The registered Foldseek workflow could not be loaded: {exc}"
    if schema_error:
        st.info(schema_error)
    parameter_specs = {str(item.get("name")): item for item in (schema or {}).get("parameters", [])}
    first, second = st.columns(2)
    similarity_threshold = float(first.number_input(
        "Foldseek TM-score threshold", min_value=0.01, max_value=1.0,
        value=float(parameter_specs.get("similarity_threshold", {}).get("default", 0.5)), step=0.05,
        key="complex-analysis-foldseek-threshold",
        help="Passed to the registered Foldseek easy-cluster protocol as --tmscore-threshold.",
    ))
    coverage = float(second.number_input(
        "Foldseek coverage", min_value=0.01, max_value=1.0,
        value=float(parameter_specs.get("coverage", {}).get("default", 0.8)), step=0.05,
        key="complex-analysis-foldseek-coverage",
        help="Passed to the registered Foldseek easy-cluster protocol as -c with cov-mode 0.",
    ))
    selection_key = ""
    try:
        from mn_molecular_design_studio.complex_analysis_sources import saved_structure_collection_selection_key
        selection_key = saved_structure_collection_selection_key(selection, project_id)
    except (TypeError, ValueError):
        pass
    existing = _matching_foldseek_jobs(
        user, selection_key=selection_key, similarity_threshold=similarity_threshold, coverage=coverage,
    ) if selection_key and valid_chains and same_project and schema else []
    current = next((item for item in existing if item["status"] in {"completed", "succeeded"}), None)
    if current is None:
        current = next((item for item in existing if item["status"] not in {"failed", "cancelled", "canceled"}), None)
    if current:
        _render_foldseek_results(user, current, key=selection_key[:14])
    elif existing:
        failed = existing[0]
        st.error(f"Previous matching Foldseek job failed: {failed.get('error') or 'No per-item error was recorded.'}")
        _render_foldseek_results(user, failed, key=selection_key[:14])

    launch_ready = (valid_chains and same_project and bool(schema) and bool(selection_key) and submit is not None)
    if current and current["status"] in {"completed", "succeeded"}:
        return
    button_label = "Retry Foldseek" if existing else "Run Foldseek"
    if st.button(button_label, type="primary", disabled=not launch_ready,
                 key=f"complex-analysis-foldseek-run-{selection_key[:14]}"):
        parameters = {
            "similarity_threshold": round(similarity_threshold, 4),
            "coverage": round(coverage, 4),
            "foldseek_collection": selection,
        }
        try:
            submit(schema, "analysis.foldseek", parameters, project_id, {}, "Foldseek")
        except (OSError, RuntimeError, TypeError, ValueError) as exc:
            st.error(f"Could not start Foldseek: {exc}")
            return
        st.rerun()


def render_comparison(user: dict, *, workflow_schema=None, project_selector=None, submit=None,
                      run_groups: list[dict] | None = None, selected_sources: list[dict] | None = None) -> None:
    """Render an additive multi-source comparison panel for existing structures."""
    st.subheader("Compare saved complexes")
    st.caption("Select 2–5 existing representations independently. Registered workflows can analyze a compatible selection while preserving each source structure and evaluation's provenance.")
    n_structures = len(selected_sources) if selected_sources else len(run_groups) if run_groups else int(st.number_input(
        "Complexes to compare", min_value=2, max_value=5, value=2, step=1,
        key="complex-analysis-compare-count",
    ))
    sources: list[dict] = list(selected_sources or [])
    for start in range(0, 0 if selected_sources else n_structures, 2):
        columns = st.columns(min(2, n_structures - start))
        for offset, column in enumerate(columns):
            slot_index = start + offset
            with column:
                if run_groups:
                    source = _comparison_slot(user, chr(ord("a") + slot_index), run_group=run_groups[slot_index])
                else:
                    source = _comparison_slot(user, chr(ord("a") + slot_index))
            if source:
                sources.append(source)
    if len(sources) != n_structures:
        return

    reference = sources[0]
    reports = [
        (chr(ord("B") + index), compatibility_report(reference, source))
        for index, source in enumerate(sources[1:])
    ]
    compatible = all(report["compatible"] for _, report in reports)
    protein_protein = all(sample.get("complex_type") == "protein-protein" for sample in sources)
    st.markdown("#### Target and component compatibility")
    if compatible:
        st.success("Compatible target sequence and author residue numbering; selected components have matching roles and ligand identity.")
    else:
        st.warning("These structures are kept as independent samples and cannot share a target-residue comparison map.")
    for label, report in reports:
        for reason in report["errors"]:
            st.caption(f"Structure A vs {label}: {reason}")
    notes = list(dict.fromkeys(note for _, report in reports for note in report["notes"]))
    for note in notes:
        st.caption(note)
    st.dataframe([
        {"Structure": f"Structure {chr(ord('A') + index)}", "Source": sample["label"], "Model": sample.get("model_id") or "Not recorded",
         "Target chains": ", ".join(sample["target_chains"]),
         "Partner / ligand": (", ".join(sample["binder_chains"]) if sample.get("complex_type") == "protein-protein"
                              else f"{sample.get('ligand_selection', {}).get('ccd') or 'Ligand'} · {sample.get('ligand_selection', {}).get('residue_key') or ''}".strip(" ·")),
         "Target sequence lengths": ", ".join(str(len(row["sequence"])) for row in sample["target_numbering"]),
         "Structure format": sample.get("format") or "Unknown"}
        for index, sample in enumerate(sources)
    ], hide_index=True, use_container_width=True)

    st.markdown("#### Comparative structure inspection")
    from mn_molecular_design_studio import complex_analysis_ui as workspace
    for start in range(0, len(sources), 2):
        columns = st.columns(min(2, len(sources) - start))
        for offset, column in enumerate(columns):
            index = start + offset
            sample = sources[index]
            with column:
                st.caption(f"Structure {chr(ord('A') + index)}")
                workspace._show_structure(sample["path"], sample.get("context") or {}, sample["components"], key=sample["viewer_key"])

    if compatible and protein_protein:
        contact_sets = {
            f"Structure {chr(ord('A') + index)}": _saved_biotite_contacts(user, sample)
            for index, sample in enumerate(sources)
        }
        st.markdown("#### Saved interface maps")
        same_target_chain_labels = all(sample["target_chains"] == reference["target_chains"] for sample in sources[1:])
        all_contact_sets = all(contact_sets.values())
        if all_contact_sets and same_target_chain_labels:
            rows = saved_target_contact_rows({label: saved["contacts"] for label, saved in contact_sets.items()})
            if rows:
                st.dataframe(rows, hide_index=True, use_container_width=True)
            st.caption("Presence is transcribed from each exact saved Biotite contact artifact. No contact distances or between-structure similarity values are recalculated.")
        elif any(contact_sets.values()):
            for start in range(0, len(sources), 2):
                columns = st.columns(min(2, len(sources) - start))
                for offset, column in enumerate(columns):
                    index = start + offset
                    label = f"Structure {chr(ord('A') + index)}"
                    saved = contact_sets[label]
                    with column:
                        st.caption(label)
                        if saved:
                            st.write(sorted({str(pair.get("target_residue")) for pair in saved["contacts"] if isinstance(pair, dict) and pair.get("target_residue")}))
                            st.caption("Saved target residue identifiers retain this structure's author chain numbering.")
                        else:
                            st.info("No matching completed Biotite contact map is saved for this exact structure, model, and chain selection.")
        else:
            st.info("No matching completed Biotite Interface evaluations are saved for these exact structures and selected chain groups. Open a structure's analysis panel to run the registered workflow.")

    elif compatible:
        st.markdown("#### Ligand identity comparison")
        st.dataframe([{
            "Structure": f"Structure {chr(ord('A') + index)}",
            "Ligand residue": sample["ligand_selection"].get("residue_key"),
            "Isomeric SMILES": sample["ligand_selection"].get("smiles"),
        } for index, sample in enumerate(sources)], hide_index=True, use_container_width=True)
        st.caption("This view compares the selected saved structures and their saved evaluations. It does not recalculate ligand interactions or pose-validation metrics.")
    if compatible:
        if reports and all(report["overlay_compatible"] for _, report in reports):
            st.markdown("#### Shared-frame overlay")
            bases = list(dict.fromkeys(report["overlay_basis"] for _, report in reports))
            st.caption("Overlay is enabled because " + "; ".join(bases) + ". Structures from different frames remain in the side-by-side viewers above.")
            if st.checkbox("Show Mol* overlay", key="complex-analysis-compare-overlay"):
                from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer import molstar_custom_component
                from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer.dataclasses import StructureVisualization
                visualizations = [
                    StructureVisualization(pdb=workspace._viewer_text(sample["path"], sample.get("context") or {}),
                                           color="chain-id", representation_type="cartoon+ball-and-stick")
                    for sample in sources
                ]
                molstar_custom_component(visualizations, key="complex-analysis-compatible-overlay", height=650,
                                         show_controls=True, show_structure_tools=False, force_reload=True)
        else:
            st.caption("These structures do not establish a common coordinate frame. Fit identical target residues explicitly below, or supply a reviewed external transform.")
        if protein_protein and len(sources) >= 2:
            st.markdown("#### Target-based multi-design overlay")
            st.caption("Each complex is independently aligned on its validated target C-alpha correspondence. Binder coordinates do not contribute to the fit. Transforms are temporary visualization data and are not saved as Evaluations.")
            if any(len(sample.get("target_chains") or []) != 1 for sample in sources):
                st.info("The selected cohort includes a multi-chain target. Current automatic target fitting requires one explicitly selected target chain; structures remain available in the side-by-side viewers.")
            elif any(len(str(chain)) != 1 for sample in sources
                     for chain in [*(sample.get("target_chains") or []), *(sample.get("binder_chains") or [])]):
                st.info("Per-chain colors in the existing Mol* overlay require single-character author chain IDs. Use the individual viewers or Inspector for this cohort.")
            else:
                labels = [f"{chr(ord('A') + index)} · {sample['label']} · model {sample.get('model_id') or 'unspecified'}"
                          for index, sample in enumerate(sources)]
                reference_index = st.selectbox("Target reference frame", list(range(len(sources))),
                    format_func=lambda index: labels[index], key="complex-analysis-target-fit-reference")
                binder_visible = st.multiselect("Binders visible in overlay", list(range(len(sources))),
                    default=list(range(len(sources))), format_func=lambda index: labels[index],
                    key="complex-analysis-target-fit-visible-binders")
                highlight_contacts = st.checkbox("Highlight saved target-contact residues where available", value=True,
                    key="complex-analysis-target-fit-highlight-contacts",
                    help="Only exact saved Biotite contact maps for these structure/model/chain identities are highlighted. Missing evaluations are not recalculated.")
                if st.checkbox("Align targets and show shared Mol* view", key="complex-analysis-fit-targets-multi"):
                    try:
                        from mn_molecular_design_studio.complex_analysis_comparison import fit_target_cohort
                        fit = fit_target_cohort(sources, reference_index=reference_index)
                        fit_rows = {int(row["index"]): row for row in fit["rows"]}
                        st.dataframe([{
                            "Design / model": labels[index],
                            "Target chain": ", ".join(row["target_chains"]),
                            "Alignment": "Reference frame" if index == reference_index else row["method"],
                            "Mapped target residues": row["mapped_residues"],
                            "Target Cα RMSD (Å)": row["target_ca_rmsd_angstrom"],
                            "Structure SHA-256": row["structure_digest"],
                        } for index, row in fit_rows.items()], hide_index=True, use_container_width=True)
                        st.caption(f"{fit['method']}. {fit['limitations']}")
                        from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer import (
                            ChainVisualization, StructureVisualization, molstar_custom_component,
                        )
                        palette = ["0x377eb8", "0xe41a1c", "0x4daf4a", "0x984ea3", "0xff7f00"]
                        visualizations = []
                        for index, sample in enumerate(sources):
                            chains = [ChainVisualization(chain_id=str(chain), color="uniform",
                                color_params={"value": "0xb8c0cc"}, label=f"Target · {chain}",
                                representation_type="cartoon") for chain in sample["target_chains"]]
                            if index in binder_visible:
                                chains.extend(ChainVisualization(chain_id=str(chain), color="uniform",
                                    color_params={"value": palette[index % len(palette)]},
                                    label=f"Binder {chr(ord('A') + index)} · {sample.get('candidate_id') or sample['label']} · model {sample.get('model_id') or 'unspecified'}",
                                    representation_type="cartoon+ball-and-stick")
                                    for chain in sample["binder_chains"])
                            highlights = []
                            if highlight_contacts:
                                saved = sample.get("saved_biotite_contacts") or {}
                                from mn_molecular_design_studio.complex_analysis_inspector import (
                                    author_residue_selector, parse_author_residue_selector,
                                )
                                for pair in saved.get("contacts") or []:
                                    parsed = parse_author_residue_selector(str(pair.get("target_residue") or "")) if isinstance(pair, dict) else None
                                    if parsed and parsed.get("chain_id") in set(sample["target_chains"]):
                                        selector = author_residue_selector(parsed["chain_id"], parsed["number"], parsed.get("insertion_code") or "")
                                        if selector:
                                            highlights.append(selector)
                            visualizations.append(StructureVisualization(
                                pdb=workspace._viewer_text(sample["path"], sample.get("context") or {}),
                                representation_type=None, chains=chains,
                                highlighted_selections=sorted(set(highlights)),
                                transform_matrix=(None if index == reference_index else fit_rows[index]["matrix_4x4"]),
                            ))
                        molstar_custom_component(visualizations, key="complex-analysis-target-fit-overlay-multi",
                            height=680, show_controls=True, show_structure_tools=False, force_reload=True)
                    except (OSError, ValueError, KeyError, ImportError) as exc:
                        st.info(f"Target-based shared view unavailable: {exc}")
        elif len(sources) == 2 and st.checkbox("Fit identical target C-alpha residues and show overlay", key="complex-analysis-fit-targets"):
            from mn_molecular_design_studio.complex_analysis_comparison import fit_matching_targets
            try:
                fit = fit_matching_targets(sources[0], sources[1])
                st.metric("Target C-alpha RMSD (Å)", f"{fit['target_ca_rmsd_angstrom']:.3f}")
                st.caption(fit["interpretation"])
                from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer import molstar_custom_component
                from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer.dataclasses import StructureVisualization
                molstar_custom_component([
                    StructureVisualization(pdb=workspace._viewer_text(sample["path"], sample.get("context") or {}),
                        color="chain-id", representation_type="cartoon+ball-and-stick",
                        transform_matrix=fit["matrix_4x4"] if index else None)
                    for index, sample in enumerate(sources)], key="complex-analysis-target-fit-overlay", height=650,
                    show_controls=True, show_structure_tools=False, force_reload=True)
                st.download_button("Download target-fit provenance", json.dumps(fit, indent=2),
                                   "target-fit.json", "application/json")
            except (OSError, ValueError, KeyError) as exc:
                st.info(f"Target-fit overlay unavailable: {exc}")
    else:
        st.info("Interface residue maps require compatible target sequence and numbering. The structures remain available for independent visual inspection.")

    if protein_protein:
        _render_interface_clustering(user, sources, compatible=compatible,
                                      workflow_schema=workflow_schema, submit=submit)
        _render_foldseek(user, sources, workflow_schema=workflow_schema, submit=submit)
    _render_saved_evaluations(sources)
    _render_supplied_alignment_overlay(sources)
    st.markdown("#### Human review")
    st.caption("Human review is separate from computational acceptance and scientific Evaluations. Saved annotations do not change the source candidate or Database record.")
    from mn_molecular_design_studio.complex_analysis_shortlist_export import render_shortlist_review_export
    render_shortlist_review_export(user, sources)
