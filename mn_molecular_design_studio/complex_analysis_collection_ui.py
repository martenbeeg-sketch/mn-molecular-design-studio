"""Database-centered Collection Map and saved Evaluation board."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import streamlit as st

from mn_molecular_design_studio.complex_analysis_collection import (
    aggregate_contact_footprints,
    assignment_rows,
    collection_page,
    component_kind,
    evaluation_table_rows,
    filter_evaluation_rows,
    heatmap_selection,
    mapped_contact_positions,
    prepare_inspector_navigation_state,
    sequence_correspondence,
)
from mn_molecular_design_studio.complex_analysis_catalog import authorized_structure_catalog


def _widget_key(prefix: str, value: str) -> str:
    return prefix + hashlib.sha256(str(value).encode()).hexdigest()[:12]


def _structure_label(record: dict, structure: dict) -> str:
    from mn_molecular_design_studio.complex_analysis_ui import _representation_label

    return f"{record.get('name') or 'Complex'} · {_representation_label(structure)} · {str(structure.get('id') or '')[-8:]}"


def _accessible_records(database, user: dict, query: str) -> list[dict]:
    """Use lightweight, ACL-filtered Database discovery for structure choices."""
    page = database.list_structure_records(user, collection="all", query=query, limit=100)
    return list(page.get("items", []))


def _merge_database_record_options(page_records: list[dict], cached_records: dict[str, dict],
                                   selected_ids: list[str]) -> dict[str, dict]:
    """Keep selected lightweight records available while a user pages or searches."""
    selected = {str(value) for value in selected_ids}
    merged = {str(row.get("id") or ""): row for row in page_records if row.get("id")}
    for record_id in selected:
        cached = cached_records.get(record_id)
        if cached is not None:
            merged.setdefault(record_id, cached)
    return merged


def _load_selected_structures(database, user: dict, records: list[dict], selected_ids: list[str]) -> list[dict]:
    """Build selected entries from authorized Database records and structures."""
    from mn_molecular_design_studio import complex_analysis_ui as workspace
    from mn_molecular_design_studio.complex_analysis_inspector import (
        observation_context, protein_sequence_tracks, structure_digest,
    )

    record_by_id = {str(row.get("id") or ""): row for row in records}
    samples = []
    for structure_id in selected_ids:
        owner = next((record for record in record_by_id.values()
                      if any(str(item.get("id") or "") == str(structure_id)
                             for item in record.get("representations") or record.get("structures") or [])), None)
        if owner is None:
            continue
        # Structure discovery intentionally uses a lightweight ACL-filtered
        # catalog. Hydrate only the explicitly selected molecule so its saved
        # Evaluations, source lineage and review annotations are available to
        # the board and Inspector without rebuilding the whole Database index.
        try:
            hydrated = database.list_molecules(
                user, collection="all", molecule_ids=[str(owner.get("id") or "")], limit=1,
            ).get("items", [])
        except (OSError, PermissionError, ValueError, TypeError):
            continue
        if not hydrated:
            continue
        record = hydrated[0]
        structures = list(record.get("representations") or record.get("structures") or [])
        structure = next((row for row in structures if str(row.get("id") or "") == str(structure_id)), None)
        if structure is None:
            continue
        try:
            path = database.structure_file(user, structure_id)
        except (OSError, PermissionError, ValueError):
            continue
        try:
            components = workspace.component_inventory(path.read_bytes(), path.name, record.get("components"))
        except (OSError, ValueError, TypeError):
            components = []

        _evaluation_rows, by_id = workspace._evaluation_options(record)
        descriptor_by_observation = {
            str(observation.get("id") or ""): descriptor
            for observation, descriptor in by_id.values()
        }
        associated = []
        for observation in record.get("observations", []):
            descriptor = descriptor_by_observation.get(str(observation.get("id") or ""))
            if descriptor and workspace._evaluation_matches_structure(descriptor, structure):
                associated.append((observation, descriptor))
        try:
            tracks = protein_sequence_tracks(
                path, sorted({str(item.get("source_chain") or "") for item in components
                              if item.get("type") == "protein" and item.get("source_chain")}),
            )
        except (ImportError, OSError, ValueError, TypeError):
            tracks = []
        try:
            inspector_data = workspace._database_inspector_data(user, record, structure, path, components, associated)
        except (OSError, PermissionError, ValueError, TypeError, KeyError):
            inspector_data = {"evaluations": [], "evidence": [], "tracks": tracks, "structure_digest": ""}

        context = structure.get("context") if isinstance(structure.get("context"), dict) else {}
        model_id = str(context.get("model_id") or context.get("model_index") or context.get("model") or structure_id)
        candidate = _database_candidate_view(record, structure, components, tracks, associated)
        samples.append({
            "identity": str(structure_id),
            "molecule_id": str(record["id"]),
            "record": record,
            "structure": structure,
            "structure_id": str(structure_id),
            "label": _structure_label(record, structure),
            "path": Path(path),
            "components": components,
            "candidate": candidate,
            "kind": component_kind(components),
            "context": context,
            "model_id": model_id,
            "tracks": tracks,
            "associated": associated,
            "inspector_data": inspector_data,
            "digest": structure_digest(path),
            "evaluation_entries": evaluation_table_rows([
                {"molecule_id": record["id"], "structure_id": structure_id,
                 "structure_label": _structure_label(record, structure),
                 "observation": observation, "evaluation": descriptor,
                 "context": observation_context(observation)}
                for observation, descriptor in associated
            ]),
        })
    return samples


def _database_candidate_view(record: dict, structure: dict, components: list[dict], tracks: list[dict],
                            associated: list[tuple[dict, dict]] | None = None) -> dict:
    """Expose saved Database provenance to cohort filters without inventing decisions.

    Database representations store their source acceptance evidence in the
    structure provenance (and the matching molecule source row). This transient
    CandidateSet-shaped view lets the existing acceptance UI read that evidence;
    it is never written back to the Database.
    """
    context = structure.get("context") if isinstance(structure.get("context"), dict) else {}
    candidate_id = str(context.get("candidate_id") or context.get("native_candidate_id") or "")
    source = next((item for item in record.get("sources") or []
                   if candidate_id and str(item.get("candidate_id") or "") == candidate_id), {})
    provenance = source if isinstance(source, dict) else {}
    source_tool = str(context.get("source_tool") or provenance.get("source_tool") or "")
    # When the selected structure was produced by a generator job, its own
    # recorded tool is the generator. A later validation Evaluation may carry
    # the validator's name in its candidate metadata and must not replace it.
    if not source_tool and str(context.get("task_group") or "").casefold() != "analysis":
        source_tool = str(context.get("tool") or "")
    decision = context.get("acceptance_decision") if isinstance(context.get("acceptance_decision"), dict) else {}
    if not decision:
        decision = provenance.get("acceptance_decision") if isinstance(provenance.get("acceptance_decision"), dict) else {}
    evidence = decision.get("evidence") if isinstance(decision.get("evidence"), dict) else {}
    raw = {}
    validation = evidence.get("validation_acceptance")
    if isinstance(validation, dict):
        raw["validation_acceptance"] = validation
        validation_provenance = context.get("validation_acceptance_provenance") or provenance.get("validation_acceptance_provenance")
        if isinstance(validation_provenance, dict):
            raw["validation_acceptance_provenance"] = validation_provenance
    campaign = evidence.get("distributed_campaign")
    if isinstance(campaign, dict):
        raw["distributed_campaign"] = campaign
    parent_candidate = context.get("source_candidate") or provenance.get("source_candidate")
    if isinstance(parent_candidate, dict):
        raw["source_candidate"] = parent_candidate
    metrics = evidence.get("native_filter_metrics")
    metrics = {str(key): value for key, value in metrics.items()
               if key in {"native_pass_filters", "pass_filters", "passes_filters", "passed_filters"}} if isinstance(metrics, dict) else {}
    binder_chains = [str(value) for value in context.get("binder_chains") or provenance.get("binder_chains") or []]
    target_chains = [str(value) for value in context.get("target_chains") or provenance.get("target_chains") or []]
    # Registered validation decisions are sometimes stored as linked Database
    # Evaluations rather than copied into the representation context. Reuse
    # only a decision attached to this exact representation and candidate; a
    # generic accepted metric or an Evaluation for another model is not enough.
    selected_representation = str(structure.get("id") or context.get("representation_id") or "")
    for observation, descriptor in associated or []:
        descriptor = descriptor if isinstance(descriptor, dict) else {}
        observed_representation = str(descriptor.get("representation_id") or "")
        if selected_representation and observed_representation and observed_representation != selected_representation:
            continue
        observation_context = (observation or {}).get("context")
        if isinstance(observation_context, str):
            try:
                observation_context = json.loads(observation_context)
            except (TypeError, ValueError):
                observation_context = {}
        observation_context = observation_context if isinstance(observation_context, dict) else {}
        observed_candidate = observation_context.get("candidate")
        observed_candidate = observed_candidate if isinstance(observed_candidate, dict) else {}
        observed_candidate_id = str(observed_candidate.get("candidate_id") or "")
        if candidate_id and observed_candidate_id and observed_candidate_id != candidate_id:
            continue
        if candidate_id and observed_candidate_id == candidate_id and not source_tool:
            # Database structures can be stored by a downstream validation
            # run. Recover the original generator from the exact linked
            # candidate record rather than treating that evaluator as the
            # generator in cohort filters and presentation.
            source_tool = str(observed_candidate.get("source_tool") or observed_candidate.get("generator") or "")
        if not binder_chains:
            binder_chains = [str(value) for value in observation_context.get("binder_chains")
                             or observed_candidate.get("binder_chains") or []]
        if not target_chains:
            target_chains = [str(value) for value in observation_context.get("target_chains")
                             or observed_candidate.get("target_chains") or []]
        observation_metrics = (observation or {}).get("metrics")
        if isinstance(observation_metrics, str):
            try:
                observation_metrics = json.loads(observation_metrics)
            except (TypeError, ValueError):
                observation_metrics = {}
        if isinstance(observation_metrics, dict):
            for key in ("native_pass_filters", "pass_filters", "passes_filters", "passed_filters",
                        "ppi.pass", "ppi.all_metrics_present"):
                if key in observation_metrics:
                    metrics.setdefault(key, observation_metrics[key])
        if str(descriptor.get("evaluator") or descriptor.get("namespace") or
               (observation or {}).get("namespace") or "").casefold() != "validation_acceptance":
            continue
        observed_raw = observed_candidate.get("raw_metadata")
        observed_raw = observed_raw if isinstance(observed_raw, dict) else {}
        observed_decision = observed_raw.get("validation_acceptance")
        if not isinstance(observed_decision, dict):
            continue
        if not raw.get("validation_acceptance"):
            raw["validation_acceptance"] = observed_decision
        if not raw.get("validation_acceptance_provenance"):
            saved_provenance = observed_raw.get("validation_acceptance_provenance")
            raw["validation_acceptance_provenance"] = saved_provenance if isinstance(saved_provenance, dict) else {
                "evaluation_id": str(descriptor.get("id") or ""),
                "run_id": str((observation or {}).get("run_id") or ""),
                "source": str(descriptor.get("evaluator") or "validation_acceptance"),
            }
    by_chain = {str(row.get("chain_id") or ""): str(row.get("sequence") or "")
                for row in tracks if isinstance(row, dict)}
    binder_sequences = [by_chain.get(chain, "") for chain in binder_chains]
    binder_sequence = ":".join(value for value in binder_sequences if value) if binder_sequences else ""
    # Some older representations store explicit roles in components rather
    # than the top-level provenance. Use only those annotations.
    if not binder_sequence:
        role_components = [item for item in components if item.get("type") == "protein"
                           and str(item.get("role") or "").casefold() == "binder"]
        if role_components:
            binder_sequence = ":".join(str(item.get("sequence") or "") for item in role_components)
    return {
        "candidate_id": candidate_id or str(record.get("id") or structure.get("id") or ""),
        "native_candidate_id": str(context.get("native_candidate_id") or provenance.get("native_candidate_id") or candidate_id),
        "model_id": str(context.get("model_id") or context.get("model_index") or ""),
        "model_index": str(context.get("model_index") or ""),
        "source_tool": source_tool or (
            str(context.get("tool") or "") if str(context.get("task_group") or "").casefold() != "analysis" else ""
        ),
        "binder_sequence": binder_sequence,
        "binder_chains": binder_chains,
        "target_chains": target_chains,
        "metrics": metrics,
        "raw_metadata": raw,
    }


def _contact_snapshot(sample: dict) -> dict:
    saved = sample.get("contact_snapshot")
    if isinstance(saved, dict):
        return saved
    evaluations = [item for item in sample.get("inspector_data", {}).get("evaluations", [])
                   if "biotite_interface" in str(item.get("evaluator") or "").casefold()]
    for item in evaluations:
        artifact = next((row for row in item.get("artifacts", []) if row.get("name") == "contacts.json"), None)
        if artifact is None:
            continue
        context = item.get("observation", {}).get("context")
        if isinstance(context, str):
            try:
                context = json.loads(context)
            except (TypeError, ValueError):
                context = {}
        context = context if isinstance(context, dict) else {}
        settings = context.get("settings") if isinstance(context.get("settings"), dict) else {}
        candidate = context.get("candidate") if isinstance(context.get("candidate"), dict) else {}
        raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
        selection = raw.get("complex_analysis_selection") if isinstance(raw.get("complex_analysis_selection"), dict) else {}
        return {
            "status": "saved",
            "rows": list(artifact.get("rows") or []),
            "evaluation": item.get("evaluation") or {},
            "observation": item.get("observation") or {},
            "settings": settings,
            "selection": selection,
            "contact_distance": settings.get("contact_distance"),
        }
    return {"status": "artifact unavailable" if evaluations else "missing evaluation",
            "rows": None, "evaluation": evaluations[0].get("evaluation", {}) if evaluations else {},
            "observation": evaluations[0].get("observation", {}) if evaluations else {},
            "settings": {}, "selection": {}, "contact_distance": None}


def _contact_profile_label(value: Any) -> str:
    if value is None or value == "":
        return "Parameters not recorded"
    try:
        return f"{float(value):g} Å"
    except (TypeError, ValueError, OverflowError):
        return f"Unparsed contact-distance setting: {value}"


def _style_contact_heatmap(frame, columns: list[str]):
    """Apply status colors across pandas versions before and after Styler.map."""
    def cell_style(value):
        return {"●": "background-color: #7252a3; color: white; text-align: center",
                "·": "background-color: #edf0f4; color: #64748b; text-align: center",
                "?": "background-color: #fff3cd; color: #7c5b00; text-align: center",
                "!": "background-color: #f8d7da; color: #842029; text-align: center",
                "×": "background-color: #f5c2c7; color: #842029; text-align: center",
                "–": "background-color: #e2e3e5; color: #555; text-align: center"}.get(value, "")

    styled = frame.style
    if columns:
        # Pandas 3 removed Styler.applymap; map is its scalar-cell replacement.
        # Keep the older method as a fallback for Studio environments on pandas <2.1.
        map_cells = getattr(styled, "map", None)
        if callable(map_cells):
            styled = map_cells(cell_style, subset=columns)
        else:
            styled = styled.applymap(cell_style, subset=columns)
    return styled


def _ligand_identity(sample: dict) -> dict[str, str]:
    context = sample.get("context") if isinstance(sample.get("context"), dict) else {}
    selection = context.get("ligand_selection") if isinstance(context.get("ligand_selection"), dict) else {}
    ligand = next((item for item in sample.get("components", [])
                   if str(item.get("type") or "").casefold() in {"small_molecule", "ligand"}), {})
    return {
        "isomeric_smiles": str(selection.get("smiles") or ligand.get("smiles") or ""),
        "ccd": str(selection.get("ccd") or ligand.get("ccd") or ""),
        "pose_residue": str(selection.get("residue_key") or ligand.get("residue_key") or ""),
        "pose_chain": str(selection.get("chain") or ligand.get("source_chain") or ligand.get("chain") or ""),
    }


def _ligand_pose_label(ligand: dict[str, str]) -> str:
    """Format the recorded ligand chain and residue selector without duplication."""
    chain = str(ligand.get("pose_chain") or "").strip()
    residue = str(ligand.get("pose_residue") or "").strip()
    if chain and residue and (residue == chain or residue.startswith(f"{chain}:")):
        return residue
    if chain and residue:
        return f"{chain}:{residue}"
    return chain or residue


def _protein_ligand_rows(samples: list[dict]) -> list[dict]:
    groups: dict[tuple, str] = {}
    rows = []
    for sample in samples:
        ligand = _ligand_identity(sample)
        receptor_sequences = tuple(sorted(track.get("sequence", "") for track in sample.get("tracks", [])))
        if ligand["isomeric_smiles"] and receptor_sequences:
            identity = (ligand["isomeric_smiles"], receptor_sequences)
            group = groups.setdefault(identity, f"Protein–ligand group {len(groups) + 1}")
            status = "Compatible ligand identity and exact receptor sequence"
        else:
            group = "Separate · missing ligand isomeric SMILES or receptor sequence"
            status = "Not grouped: chemical identity or receptor sequence is unavailable"
        rows.append({
            "Structure": sample["label"], "Group": group, "Compatibility": status,
            "Ligand isomeric SMILES": ligand["isomeric_smiles"] or "Not recorded",
            "Ligand pose": _ligand_pose_label(ligand),
            "Receptor protein sequence lengths": ", ".join(str(len(value)) for value in receptor_sequences),
            "Saved evaluations": ", ".join(sorted({str(row.get("method")) for row in sample["evaluation_entries"]})) or "None",
            "Structure ID": sample["structure_id"],
            "Structure digest": sample["digest"],
        })
    return rows


def _render_evaluation_board(samples: list[dict]) -> list[dict]:
    entries = []
    seen = set()
    for sample in samples:
        for row in sample.get("evaluation_entries", []):
            identity = row.get('evaluation_id') or (row.get('run_id'), row.get('candidate_id'), row.get('structure_digest'), row.get('method'))
            if identity not in seen:
                entries.append(row)
                seen.add(identity)
    st.markdown("#### Evaluation Board")
    st.caption("Scientific findings and model confidence remain separate from human review. This board filters only saved Evaluations and their recorded provenance; it does not create a combined quality score.")
    if not entries:
        st.info("No saved Evaluation is associated with the selected structure representations. Missing Evaluations remain distinct from zero-valued measurements.")
        return []

    methods = sorted({str(row.get("method") or "unknown") for row in entries})
    engines = sorted({str(row.get("source_engine") or "unknown") for row in entries})
    numeric_metrics = sorted({str(key) for row in entries for key, value in (row.get("metrics") or {}).items()
                              if type(value) in {int, float}})
    first, second, third = st.columns([1, 1, 1.2])
    method = first.selectbox("Method", ["All methods", *methods], key="collection-eval-method")
    engine = second.selectbox("Source engine", ["All engines", *engines], key="collection-eval-engine")
    metric = third.selectbox("Measured value", ["No numeric filter", *numeric_metrics], key="collection-eval-metric")
    minimum = maximum = None
    if metric != "No numeric filter":
        values = [float(row["metrics"][metric]) for row in entries
                  if type((row.get("metrics") or {}).get(metric)) in {int, float}]
        low, high = (min(values), max(values)) if values else (0.0, 0.0)
        range_left, range_right = st.columns(2)
        minimum = float(range_left.number_input("Minimum value", value=low, key=_widget_key("collection-eval-min-", metric)))
        maximum = float(range_right.number_input("Maximum value", value=high, key=_widget_key("collection-eval-max-", metric)))
    filtered = filter_evaluation_rows(
        entries, method="" if method == "All methods" else method,
        source_engine="" if engine == "All engines" else engine,
        metric="" if metric == "No numeric filter" else metric,
        minimum=minimum, maximum=maximum,
    )
    import pandas as pd

    comparison_columns = st.multiselect('Comparison measurement columns', numeric_metrics,
        default=[key for key in ('i_pTM', 'i_pDAE', 'qc_interface_contact_pairs', 'qc_interface_residue_pairs', 'plip_interaction_count') if key in numeric_metrics],
        key='collection-eval-comparison-columns')

    display_rows = [{
        "Structure": row["structure_label"], "Kind": row["category"], "Method": row["method"],
        "Source engine": row["source_engine"], "Status": row["status"], "Created": row["created_at"],
        "Evaluation ID": row["evaluation_id"], "Run ID": row["run_id"] or "Not recorded",
        "Source job": row.get('source_job') or 'Not recorded',
        "Source workstation": row.get('source_workstation') or 'Not recorded',
        "Source candidate": row.get('source_candidate_id') or 'Not recorded',
        "Source model": row.get('source_model_id') or 'Not recorded',
        "Chains": json.dumps({'binder': row.get('binder_chains'), 'target': row.get('target_chains')})
                  if 'binder_chains' in row else 'Not recorded',
        "Candidate / pose ID": row["candidate_id"] or "Not recorded",
        "Structure digest": row["structure_digest"] or "Not recorded",
        **({"Native acceptance": row['native_acceptance'],
            "Exploratory acceptance": row.get('exploratory_acceptance', 'not_evaluated'),
            "Manual selection": row.get('manual_selection', False)} if 'native_acceptance' in row else {}),
        "Parameters": json.dumps(row["parameters"], sort_keys=True, default=str),
        "Recorded measurements": json.dumps(row["metrics"], sort_keys=True, default=str),
        **{key: row['metrics'].get(key) for key in comparison_columns},
    } for row in filtered]
    st.caption(f"{len(filtered)} saved Evaluation row(s) match these filters; {len(entries)} are associated with the selected structures.")
    st.dataframe(pd.DataFrame(display_rows), hide_index=True, use_container_width=True)

    assignments = assignment_rows(entries)
    cluster_evaluations = [row for row in entries
                           if "interface_clustering" in str(row.get("method") or "").casefold()
                           or "foldseek" in str(row.get("method") or "").casefold()]
    if assignments:
        st.markdown("##### Saved collection assignments")
        st.dataframe(pd.DataFrame(assignments), hide_index=True, use_container_width=True)
        sample_by_id = {sample["structure_id"]: sample for sample in samples}
        member_assignments = [row for row in assignments if row.get("Structure ID") in sample_by_id]
        if member_assignments:
            assignment_ids = [str(row.get("Evaluation ID") or "") for row in member_assignments]
            assignment_lookup = {str(row.get("Evaluation ID") or ""): row for row in member_assignments}
            chosen_id = st.selectbox(
                "Saved assignment structure",
                assignment_ids,
                format_func=lambda value: (
                    f"{assignment_lookup[value]['Method']} · {assignment_lookup[value]['Structure']} · "
                    f"cluster {assignment_lookup[value]['Assignment']} · {value[-8:]}"
                ),
                key="collection-assignment-structure",
            )
            chosen = assignment_lookup[chosen_id]
            if st.button("Open saved assignment in Scientific Inspector", key="collection-assignment-open"):
                _open_in_inspector(sample_by_id[chosen["Structure ID"]], "", chosen_id)
        representatives = [row for row in assignments if row.get("Representative structure ID") in sample_by_id]
        if representatives:
            representative_ids = sorted({str(row["Representative structure ID"]) for row in representatives})
            representative_labels = {str(row["Representative structure ID"]): row["Representative structure"]
                                     for row in representatives}
            representative_id = st.selectbox(
                "Saved representative structure", representative_ids,
                format_func=lambda value: representative_labels.get(value, value),
                key="collection-assignment-representative",
            )
            if st.button("Open representative in Scientific Inspector", key="collection-assignment-representative-open"):
                _open_in_inspector(sample_by_id[representative_id])
    elif cluster_evaluations:
        st.markdown("##### Saved collection workflow Evaluations")
        st.dataframe(pd.DataFrame([{
            "Structure": row["structure_label"],
            "Method": ("Interface Clustering · binding-mode similarity" if "interface_clustering" in row["method"].casefold()
                       else "Foldseek · structural similarity"),
            "Recorded outputs": json.dumps(row["metrics"], sort_keys=True, default=str),
            "Evaluation ID": row["evaluation_id"], "Run ID": row["run_id"],
            "Parameters": json.dumps(row["parameters"], sort_keys=True, default=str),
        } for row in cluster_evaluations]), hide_index=True, use_container_width=True)
    st.caption("Interface Clustering assignments describe binding-mode similarity. Foldseek assignments describe structural similarity and are shown separately.")

    evaluation_ids = [str(row["evaluation_id"]) for row in filtered if row.get("evaluation_id")]
    if evaluation_ids:
        labels = {str(row["evaluation_id"]): f"{row['method']} · {row['structure_label']} · {str(row['evaluation_id'])[-8:]}"
                  for row in filtered if row.get("evaluation_id")}
        selected_id = st.selectbox("Evaluation provenance detail", evaluation_ids,
                                   format_func=lambda value: labels.get(value, value),
                                   key="collection-eval-detail")
        selected = next(row for row in filtered if str(row["evaluation_id"]) == selected_id)
        st.json({"evaluation_id": selected["evaluation_id"], "observation_id": selected["observation_id"],
                 "method": selected["method"], "source_engine": selected["source_engine"],
                 "status": selected["status"], "metrics": selected["metrics"],
                 "parameters": selected["parameters"], "source_context": selected["provenance"]})
    return filtered


def _render_batch_assessment(user: dict, samples: list[dict], *, workflow_schema=None, submit=None) -> None:
    """Explicitly queue registered analyses against selected Database structures."""
    database_samples = [sample for sample in samples if sample.get("source_kind") == "Database"]
    if not database_samples:
        return
    from mn_molecular_design_studio.complex_analysis_batch import (
        ASSESSMENT_METHODS, CORE_PROFILE, COMPREHENSIVE_PROFILE,
        assessment_applicability, evaluation_coverage, recorded_chain_roles,
    )
    from mn_molecular_design_studio.complex_analysis_sources import _database_chain_sequence
    from mn_molecular_design_studio.molecular_database import MolecularDatabase

    database = MolecularDatabase()
    st.markdown("#### Analyze selected designs")
    st.caption("This action runs existing registered Studio evaluators only when you submit. Each result is attached to the exact selected Database design and structure through the normal Evaluation ingestion path.")
    profile_name = st.selectbox("Assessment profile", ["Core binder assessment", "Comprehensive binder assessment", "Custom"],
                                key="complex-analysis-batch-profile")
    profile_methods = list(CORE_PROFILE if profile_name == "Core binder assessment" else
                           COMPREHENSIVE_PROFILE if profile_name == "Comprehensive binder assessment" else ())
    default_methods = [method for method in profile_methods if method in ASSESSMENT_METHODS]
    profile_state_key, method_state_key = "complex-analysis-batch-last-profile", "complex-analysis-batch-methods"
    if st.session_state.get(profile_state_key) != profile_name:
        st.session_state[method_state_key] = default_methods
        st.session_state[profile_state_key] = profile_name
    methods = st.multiselect("Registered methods", list(ASSESSMENT_METHODS), default=default_methods,
        format_func=lambda value: ASSESSMENT_METHODS[value]["label"], key=method_state_key,
        help="Interface Clustering and Foldseek remain cohort-level methods in Detailed structural comparison. They do not produce one independent per-design evaluation.")
    st.caption("Biotite and Rosetta interface geometry, ipSAE, PEP-Patch, sequence composition, Protein-Sol, liabilities and tag/fusion suitability keep their own method definitions and units. No universal score is calculated.")
    first, second = st.columns(2)
    contact_distance = float(first.number_input("Biotite contact distance (Å)", min_value=0.1, max_value=20.0,
        value=5.0, step=0.1, key="complex-analysis-batch-biotite-distance"))
    ph = float(second.number_input("Sequence / PEP-Patch pH", min_value=0.0, max_value=14.0,
        value=7.4, step=0.1, key="complex-analysis-batch-ph"))
    liability_format_options = {"": "Choose construct format…", "linear_protein": "Linear protein chain(s)",
                                "cyclic_construct": "Cyclic construct (no free termini)",
                                "multichain_construct": "Multichain construct (chain-by-chain screen)"}
    liability_format = st.selectbox("Sequence-liability construct format", list(liability_format_options),
        format_func=lambda value: liability_format_options[value], key="complex-analysis-batch-liability-format",
        disabled="analysis.sequence_liability" not in methods)
    liability_plddt = st.checkbox("The selected structure stores pLDDT in Cα B-factors (verified from its source)",
        value=False, key="complex-analysis-batch-liability-plddt", disabled="analysis.sequence_liability" not in methods)
    tag_plddt = st.checkbox("The selected structure stores pLDDT in Cα B-factors (verified from its source)",
        value=False, key="complex-analysis-batch-tag-plddt", disabled="analysis.tag_fusion_suitability" not in methods)
    tag_columns = st.columns(2)
    tag_format = tag_columns[0].selectbox("Tag/fusion component format",
        ["", "single_chain_linear_protein", "cyclic_peptide", "multichain_construct"],
        format_func=lambda value: {"": "Choose format…", "single_chain_linear_protein": "Single-chain linear protein",
                                   "cyclic_peptide": "Cyclic peptide", "multichain_construct": "Multichain construct"}[value],
        key="complex-analysis-batch-tag-format", disabled="analysis.tag_fusion_suitability" not in methods)
    tag_termini_confirmed = tag_columns[1].checkbox(
        "Complete intended termini are present (verified)", value=False,
        key="complex-analysis-batch-tag-termini", disabled="analysis.tag_fusion_suitability" not in methods)
    tag_cutoff = float(st.number_input("Tag/fusion interface contact cutoff (Å)", min_value=1.0, max_value=10.0,
        value=4.5, step=0.25, key="complex-analysis-batch-tag-cutoff",
        disabled="analysis.tag_fusion_suitability" not in methods))
    ipsae_defaults = {"pae_cutoff": 10.0, "dist_cutoff": 10.0, "ipdae_contact_cutoff": 8.0,
                     "structure_contact_cutoff": 4.0, "hotspot_contact_cutoff": 5.0}
    ipsae_parameters = {}
    if "analysis.ipsae" in methods:
        with st.expander("ipSAE and structure-analysis cutoffs", expanded=False):
            fields = st.columns(3)
            for index, (name, default) in enumerate(ipsae_defaults.items()):
                ipsae_parameters[name] = round(float(fields[index % len(fields)].number_input(
                    name.replace("_", " ").title() + " (Å)", min_value=0.1, max_value=50.0,
                    value=default, step=0.5, key="complex-analysis-batch-ipsae-" + name)), 3)
    st.caption("Cohort-level Interface Clustering and Foldseek remain available under Detailed structural comparison below. They produce saved per-structure assignments, not an independent scalar per design.")
    include_mismatched = st.checkbox("Allow a new Evaluation when an earlier one used different settings", value=False,
        key="complex-analysis-batch-include-mismatched")
    retry_failed = st.checkbox("Retry previously failed requested analyses", value=False,
        key="complex-analysis-batch-retry-failed")

    method_schemas: dict[str, dict | None] = {}
    schema_errors: dict[str, str] = {}
    for method in methods:
        try:
            method_schemas[method] = workflow_schema(method) if workflow_schema else None
            if not method_schemas[method]:
                schema_errors[method] = "The registered workflow schema is unavailable."
        except Exception as exc:
            method_schemas[method] = None
            schema_errors[method] = f"The registered workflow schema could not be loaded: {exc}"

    role_map: dict[str, tuple[list[str], list[str]]] = {}
    for sample in database_samples:
        component_chains = sorted({str(item.get("source_chain") or item.get("chain") or "")
                                   for item in sample.get("components") or []
                                   if isinstance(item, dict) and str(item.get("type") or "").casefold() == "protein"})
        recorded_binder, recorded_target = recorded_chain_roles(sample)
        token = hashlib.sha256(str(sample.get("identity") or "").encode()).hexdigest()[:12]
        with st.expander(f"Chain assignment · {sample.get('label') or sample.get('structure_id')}",
                         expanded=not bool(recorded_binder)):
            st.caption(f"Model {sample.get('model_id') or 'not recorded'} · structure SHA-256 {sample.get('digest') or 'unavailable'}")
            binder_key, target_key = f"batch-assess-binder-{token}", f"batch-assess-target-{token}"
            if binder_key not in st.session_state:
                st.session_state[binder_key] = recorded_binder
            if target_key not in st.session_state:
                st.session_state[target_key] = recorded_target
            binder = st.multiselect("Protein component(s) to assess", component_chains,
                key=binder_key, help="For sequence methods, choose exactly one chain. For interface methods, these are the designed/binder chains.")
            target_options = [chain for chain in component_chains if chain not in binder]
            st.session_state[target_key] = [chain for chain in st.session_state.get(target_key, []) if chain in target_options]
            target = st.multiselect("Interface partner / target chains", target_options, key=target_key)
            role_map[str(sample.get("identity") or sample.get("structure_id"))] = (list(binder), list(target))
            if not recorded_binder and not binder:
                st.info("No saved binder role was found. Select the component explicitly; chain order is not used to guess it.")

    try:
        from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
    except ImportError:
        read_candidates = None

    rows, submission_plan = [], []
    for sample in database_samples:
        identity = str(sample.get("identity") or sample.get("structure_id"))
        binder, target = role_map.get(identity, ([], []))
        group_id = str((sample.get("structure") or {}).get("group_id") or "")
        molecule_id = str(sample.get("molecule_id") or "")
        structure_id = str(sample.get("structure_id") or "")
        model_id = str(sample.get("model_id") or "")
        digest = str(sample.get("digest") or "")
        try:
            group = database.group(group_id, user)
            project_id = str(group.get("project_id") or "")
            group_origin = group.get("origin") or {}
        except (PermissionError, OSError, ValueError, KeyError) as exc:
            project_id, group_origin = "", {}
            role_error = f"Source group is inaccessible: {exc}"
        else:
            role_error = ""
        source = {"molecule_id": molecule_id, "structure_id": structure_id, "group_id": group_id,
                  "model_id": model_id, "structure_digest": digest, "binder_chains": binder,
                  "target_chains": target, "structure_context": sample.get("context") or {},
                  "group_origin": group_origin}
        for method in methods:
            method_target = [] if method in {"analysis.sequence_qc", "analysis.proteinsol", "analysis.pep_patch"} else target
            parameters: dict[str, Any] = {}
            if method == "analysis.biotite_interface":
                parameters = {"contact_distance": round(contact_distance, 3)}
            elif method == "analysis.pep_patch":
                parameters = {"ph": round(ph, 3)}
            elif method == "analysis.sequence_qc":
                parameters = {"ph": round(ph, 3)}
            elif method == "analysis.sequence_liability":
                parameters = {"component_format": liability_format,
                              "b_factor_interpretation": "plddt_0_100" if liability_plddt else "uninterpreted"}
            elif method == "analysis.tag_fusion_suitability":
                parameters = {"interface_contact_cutoff": tag_cutoff, "termini_confirmed": tag_termini_confirmed,
                              "component_format": tag_format,
                              "b_factor_interpretation": "plddt_0_100" if tag_plddt else "uninterpreted"}
            elif method in {"analysis.rosetta", "analysis.ipsae"}:
                parameters = {str(item.get("name")): item.get("default")
                              for item in (method_schemas.get(method) or {}).get("parameters", [])
                              if isinstance(item, dict) and item.get("name") and item.get("default") is not None}
                if method == "analysis.ipsae":
                    parameters.update(ipsae_parameters)
            elif method == "analysis.proteinsol":
                parameters = {}
            compatible, reason = assessment_applicability(sample, method, binder, target)
            if method == "analysis.sequence_liability" and not liability_format:
                compatible, reason = False, "Choose the construct format before submitting sequence-liability screening."
            if method == "analysis.tag_fusion_suitability" and not tag_format:
                compatible, reason = False, "Choose the component format before submitting tag/fusion assessment."
            if role_error:
                compatible, reason = False, role_error
            if not schema_errors.get(method) and not method_schemas.get(method):
                compatible, reason = False, schema_errors.get(method, "Schema unavailable")
            elif schema_errors.get(method):
                compatible, reason = False, schema_errors[method]
            status, detail = (evaluation_coverage(sample, method, parameters=parameters,
                         binder_chains=binder, target_chains=method_target) if compatible else ("Unsupported input", reason))
            if compatible and status in {"Missing evaluation", "Existing Evaluation · different settings", "Previously failed"}:
                from mn_molecular_design_studio import complex_analysis_ui as workspace
                candidate_id = "database-" + hashlib.sha256((molecule_id + ":" + structure_id).encode()).hexdigest()[:24]
                try:
                    matching_jobs = workspace._matching_database_analysis_jobs(
                        user, workflow_id=method, source=source, candidate_id=candidate_id,
                        model_id=model_id, binder_chains=binder, target_chains=method_target,
                        parameters=parameters,
                    )
                except (OSError, PermissionError, RuntimeError, TypeError, ValueError):
                    matching_jobs = []
                if matching_jobs:
                    existing_job = matching_jobs[0]
                    job_status = str(existing_job.get("status") or "").casefold()
                    if job_status in {"completed", "succeeded"}:
                        status, detail = "Completed job · Evaluation ingestion pending", str(existing_job.get("run_id") or "Completed job")
                    elif job_status in {"queued", "pending", "submitted", "scheduled", "starting", "preparing", "running", "staging", "waiting"}:
                        status, detail = "Queued / running", str(existing_job.get("run_id") or job_status)
                    elif job_status in {"failed", "cancelled", "canceled"}:
                        status, detail = "Previously failed", str(existing_job.get("error") or existing_job.get("run_id") or "Previous job failed")
            rows.append({"Design / structure": sample.get("label") or structure_id,
                         "Method": ASSESSMENT_METHODS[method]["label"], "Status": status,
                         "Reason / Evaluation": detail})
            queue = compatible and status == "Missing evaluation"
            queue = queue or (compatible and status == "Existing Evaluation · different settings" and include_mismatched)
            queue = queue or (compatible and status == "Previously failed" and retry_failed)
            if queue and project_id is not None:
                candidate_id = "database-" + hashlib.sha256((molecule_id + ":" + structure_id).encode()).hexdigest()[:24]
                source_selection = {"molecule_id": molecule_id, "structure_id": structure_id,
                    "group_id": group_id, "model_id": model_id, "binder_chains": binder,
                    "target_chains": method_target}
                parameters.update({"candidate_ids": [candidate_id], "selected_structure_path": str(sample.get("path") or ""),
                                   "selected_model_id": model_id, "binder_chains": binder,
                                   "target_chains": method_target, "database_source": source_selection})
                submission_plan.append({"method": method, "sample": sample, "parameters": parameters,
                                        "project_id": project_id, "source": source})

    if not rows:
        st.info("Select at least one registered assessment method to inspect coverage.")
        return
    import pandas as pd
    st.markdown("##### Analysis coverage")
    st.dataframe(pd.DataFrame(rows), hide_index=True, use_container_width=True)
    if submission_plan:
        st.caption(f"Ready to submit {len(submission_plan)} explicit design × method job(s). Current Database-source adapter submits one exact model at a time; no missing Evaluations are launched automatically.")
    else:
        st.caption("No missing requested Evaluations are ready to submit. Compatible saved work is reused; different settings and previous failures require the explicit controls above.")
    if len(submission_plan) > 50:
        st.warning("This selection would queue more than 50 jobs. Narrow the selected designs or methods before submission.")
    if st.button("Analyze selected designs", type="primary", key="complex-analysis-submit-batch-assessment",
                 disabled=submit is None or not submission_plan or len(submission_plan) > 50):
        submitted, errors = 0, []
        for item in submission_plan:
            method = item["method"]
            try:
                schema = method_schemas[method]
                submit(schema, method, item["parameters"], item["project_id"], {},
                       ASSESSMENT_METHODS[method]["label"])
                submitted += 1
            except (OSError, RuntimeError, TypeError, ValueError, PermissionError) as exc:
                errors.append(f"{ASSESSMENT_METHODS[method]['label']} · {item['sample'].get('label')}: {exc}")
        if submitted:
            st.success(f"Submitted {submitted} analysis job(s). Evaluations will appear here after the normal Results ingestion completes.")
        for error in errors:
            st.error(error)


def _open_in_inspector(sample: dict, selector: str = "", evaluation_id: str = "") -> None:
    from mn_molecular_design_studio import complex_analysis_ui as workspace

    if sample.get("source_kind") == "Results":
        job = sample["job"]
        if sample.get("remote"):
            workspace.open_remote_result_complex(
                job, sample["candidate"], sample["relative_path"], user=sample["user"],
                structure_content=Path(sample["path"]).read_bytes(),
                source_run_dir=str(sample.get("source_run_dir") or job.get("run_dir") or ""),
                structure_field=str(sample.get("structure_field") or "complex_pdb"),
                structure_sha256=str(sample.get("digest") or ""),
            )
        else:
            if selector:
                viewer_key = (f"complex-analysis-results-viewer-{job.get('task_group')}/{job.get('run_id')}"
                              f"-{sample['relative_path']}")
                st.session_state[f"{viewer_key}-inspector-highlights"] = [selector]
            workspace.open_result_complex(
                job["run_dir"], sample["path"], sample["candidate"], user=sample["user"],
            )
        return

    prepare_inspector_navigation_state(
        st.session_state, molecule_id=sample["molecule_id"], structure_id=sample["structure_id"],
        evaluation_id=evaluation_id, selector=selector,
    )
    workspace.open_database_complex(
        sample["molecule_id"], evaluation_id=evaluation_id, structure_id=sample["structure_id"],
    )


def _render_protein_protein_map(samples: list[dict]) -> None:
    from mn_molecular_design_studio.complex_analysis_inspector import protein_sequence_tracks

    st.markdown("#### Protein–protein contact map")
    if len(samples) < 2:
        st.info("Select at least two compatible protein–protein structures from Results or the Database to build a collection map.")
        return
    sample_ids = [str(sample.get("identity") or sample["structure_id"]) for sample in samples]
    sample_labels = {str(sample.get("identity") or sample["structure_id"]): sample["label"] for sample in samples}
    anchor_id = st.selectbox("Residue correspondence anchor", sample_ids,
                             format_func=lambda value: sample_labels.get(value, value),
                             key="collection-map-anchor")
    anchor = next(sample for sample in samples if str(sample.get("identity") or sample["structure_id"]) == anchor_id)
    chain_ids = [str(track.get("chain_id") or "") for track in anchor.get("tracks", [])]
    if not chain_ids:
        st.warning("The selected anchor has no single-model protein sequence track. It remains in the Evaluation Board but cannot define a residue map.")
        return
    anchor_chain = st.selectbox("Component to track", chain_ids,
                                format_func=lambda value: f"Chain {value} · {len(next(track for track in anchor['tracks'] if track['chain_id'] == value)['sequence'])} residues",
                                key=_widget_key("collection-map-anchor-chain-", anchor_id))
    reference_track = next(track for track in anchor["tracks"] if track["chain_id"] == anchor_chain)
    reference_sequence = reference_track["sequence"]
    tracked = []
    compatibility_rows = []
    with st.expander("Component correspondence and group membership", expanded=len(samples) <= 8):
        for sample in samples:
            tracks = sample.get("tracks", [])
            compatible = []
            reports = {}
            for track in tracks:
                report = sequence_correspondence(reference_sequence, track.get("sequence", ""))
                reports[str(track.get("chain_id") or "")] = report
                if report["compatible"]:
                    compatible.append(track)
            sample_identity = str(sample.get("identity") or sample["structure_id"])
            if sample_identity == anchor_id:
                selected_chain = anchor_chain
            elif len(compatible) == 1:
                selected_chain = str(compatible[0]["chain_id"])
            elif len(compatible) > 1:
                options = [str(track["chain_id"]) for track in compatible]
                selected_chain = st.selectbox(
                    f"Map component in {sample['label']}", options,
                    format_func=lambda value: f"Chain {value} · {reports[value]['status']}",
                    key=_widget_key("collection-map-chain-", sample["structure_id"]),
                )
            else:
                selected_chain = ""
            report = reports.get(selected_chain) if selected_chain else None
            track = next((row for row in tracks if str(row.get("chain_id") or "") == selected_chain), None)
            if sample_identity == anchor_id:
                report = {"status": "exact", "compatible": True,
                          "mapping": list(range(len(reference_sequence))), "reason": "Anchor component."}
                track = reference_track
                selected_chain = anchor_chain
            status = report["status"] if report else "unsupported"
            reason = report["reason"] if report else (
                "No protein chain has a unique exact sequence correspondence." if tracks
                else "No single-model protein sequence track is available.")
            compatibility_rows.append({"Structure": sample["label"], "Mapped component": selected_chain or "None",
                                       "Correspondence": status, "Explanation": reason})
            if report and report.get("compatible") and track:
                member = dict(sample)
                member.update({"tracked_chain": selected_chain, "tracked_track": track,
                               "identity": sample_identity,
                               "candidate_to_reference": report["mapping"], "correspondence": report})
                member["contact"] = _contact_snapshot(sample)
                tracked.append(member)
    import pandas as pd
    st.dataframe(pd.DataFrame(compatibility_rows), hide_index=True, use_container_width=True)
    incompatible = [row for row in compatibility_rows if row["Correspondence"] not in {"exact", "unique_exact_subsequence", "reference_is_unique_subsequence"}]
    if incompatible:
        st.info(f"{len(incompatible)} selected protein–protein structure(s) remain outside this residue group because their component correspondence is unsupported or ambiguous. Their saved Evaluations remain in the board.")
    if not tracked:
        return

    profiles = sorted({_contact_profile_label(member["contact"].get("contact_distance"))
                       for member in tracked if member["contact"]["status"] == "saved"})
    mixed_profiles = len(profiles) > 1
    profile_options = (["Choose one contact-distance setting", *profiles] if mixed_profiles
                       else ["All compatible saved maps", *profiles])
    if st.session_state.get("collection-map-contact-profile") not in profile_options:
        st.session_state["collection-map-contact-profile"] = profile_options[0]
    profile = st.selectbox("Saved Biotite map settings", profile_options,
                           key="collection-map-contact-profile")
    map_members = []
    for member in tracked:
        contact = member["contact"]
        distance = contact.get("contact_distance")
        label = _contact_profile_label(distance)
        if profile == "Choose one contact-distance setting" and contact["status"] == "saved":
            status = "different parameters"
            rows = None
        elif profile not in {"All compatible saved maps", "Choose one contact-distance setting"} and contact["status"] == "saved" and label != profile:
            status = "different parameters"
            rows = None
        else:
            status = contact["status"]
            rows = contact.get("rows")
        map_members.append({**member, "contact_status": status, "contact_rows": rows,
                            "contact_profile": label,
                            "contact_positions": (mapped_contact_positions(
                                rows, chain_id=member["tracked_chain"], track=member["tracked_track"],
                                candidate_to_reference=member["candidate_to_reference"],
                            ) if isinstance(rows, list) else None)})
    saved_profiles = {row["contact_profile"] for row in map_members if row["contact_status"] == "saved"}
    if mixed_profiles and profile == "Choose one contact-distance setting":
        st.warning("These saved maps use different contact-distance settings. Select one setting to show a comparable residue-frequency map; the mixed cohort remains visible in the Evaluation Board.")
    footprints = aggregate_contact_footprints(map_members, reference_track)
    evaluable_count = next((row["evaluable_count"] for row in footprints), 0) if footprints else 0
    st.caption(f"Contact frequency is transcribed from saved Biotite Interface artifacts on the selected component. Denominator: {evaluable_count} saved map(s) with these settings. A saved empty map is zero contacts; missing or inaccessible Evaluations are shown separately.")

    st.markdown("##### Contact-frequency sequence track")
    if not footprints:
        st.info("The selected component has no mapped residues.")
    else:
        page_size = 24
        page_count = max(1, (len(footprints) + page_size - 1) // page_size)
        segment = st.selectbox("Sequence segment", list(range(page_count)),
                               format_func=lambda value: f"Residues {value * page_size + 1}–{min((value + 1) * page_size, len(footprints))} of {len(footprints)}",
                               key="collection-map-sequence-segment")
        start = segment * page_size
        visible_track = footprints[start:start + page_size]
        st.caption("Residues are ordered by sequence; labels retain anchor author numbering and insertion codes. Select a tile to carry that component residue into the existing Inspector.")
        for row_start in range(0, len(visible_track), 12):
            columns = st.columns(12)
            for column, residue in zip(columns, visible_track[row_start:row_start + 12]):
                suffix = residue["insertion_code"]
                count = "—" if residue["evaluable_count"] == 0 else f"{residue['contact_count']}/{residue['evaluable_count']}"
                label = f"{residue['amino_acid']}\n{residue['author_number']}{suffix}\n{count}"
                if column.button(label, key=_widget_key("collection-map-residue-", f"{anchor_id}:{anchor_chain}:{residue['reference_index']}"),
                                 help=f"{residue['amino_acid']} · author {residue['author_number']}{suffix} · {count} saved maps"):
                    st.session_state["complex-analysis-collection-residue-index"] = residue["reference_index"]

    st.markdown("##### Residue footprint by complex")
    rows_per_page = int(st.selectbox("Complexes per heatmap page", [10, 25, 50, 100], index=1,
                                     key="collection-map-page-size"))
    max_page = max(0, (len(map_members) - 1) // rows_per_page)
    page = int(st.number_input("Heatmap page (0-based)", min_value=0, max_value=max_page, value=0,
                               key="collection-map-page"))
    visible_members = collection_page(map_members, page=page, page_size=rows_per_page)
    residue_page_size = 100
    residue_page_count = max(1, (len(footprints) + residue_page_size - 1) // residue_page_size)
    residue_page = int(st.selectbox("Residue columns", list(range(residue_page_count)),
                                    format_func=lambda value: f"Residues {value * residue_page_size + 1}–{min((value + 1) * residue_page_size, len(footprints))}",
                                    key="collection-map-heatmap-residue-page"))
    visible_residues = collection_page(footprints, page=residue_page, page_size=residue_page_size)
    residue_columns = {
        f"{residue['amino_acid']}{residue['author_number']}{residue['insertion_code']}": residue["reference_index"]
        for residue in visible_residues
    }
    import pandas as pd
    matrix = []
    for member in visible_members:
        row = {"Complex": member["label"], "Biotite Evaluation": member["contact_status"]}
        positions = member.get("contact_positions")
        for column, reference_index in residue_columns.items():
            if member["contact_status"] == "saved":
                row[column] = "●" if reference_index in (positions or set()) else "·"
            elif member["contact_status"] == "missing evaluation":
                row[column] = "?"
            elif member["contact_status"] == "artifact unavailable":
                row[column] = "!"
            elif member["contact_status"] == "failed evaluation":
                row[column] = "×"
            elif member["contact_status"] == "different parameters":
                row[column] = "–"
            else:
                row[column] = "?"
        matrix.append(row)
    frame = pd.DataFrame(matrix)
    style_columns = list(residue_columns)
    styled = _style_contact_heatmap(frame, style_columns)
    grid = st.dataframe(styled, hide_index=True, use_container_width=True, height=420,
                        on_select="rerun", selection_mode="single-cell",
                        key="complex-analysis-collection-footprint-heatmap")
    st.caption("Legend: ● saved contact · · saved map with no contact · ? no matching saved Evaluation · × Evaluation failed · ! Evaluation artifact unavailable · – saved map uses another selected setting.")
    selection_state = getattr(grid, "selection", {}) or {}
    cells = selection_state.get("cells", []) if hasattr(selection_state, "get") else []
    visible_ids = [str(member.get("identity") or member["structure_id"]) for member in visible_members]
    selected_cell = heatmap_selection(cells, visible_ids, residue_columns)

    st.markdown("##### Open a selected complex in the Scientific Inspector")
    mapped_members = [member for member in tracked if member["correspondence"].get("compatible")]
    mapped_ids = [str(member.get("identity") or member["structure_id"]) for member in mapped_members]
    mapped_by_id = {str(member.get("identity") or member["structure_id"]): member for member in mapped_members}
    default_id = selected_cell["member_id"] if selected_cell else anchor_id
    if default_id not in mapped_ids:
        default_id = mapped_ids[0] if mapped_ids else ""
    if mapped_ids:
        selected_id = st.selectbox("Complex for sequence selection", mapped_ids,
                                   index=mapped_ids.index(default_id),
                                   format_func=lambda value: sample_labels.get(value, value),
                                   key="collection-map-open-member")
        selected_member = mapped_by_id[selected_id]
        residue_index = (selected_cell["reference_index"] if selected_cell else
                         st.session_state.get("complex-analysis-collection-residue-index"))
        if residue_index is not None:
            local_index = next((index for index, value in enumerate(selected_member["candidate_to_reference"])
                                if value == int(residue_index)), None)
            selector = (selected_member["tracked_track"]["residues"][local_index]["selector"]
                        if local_index is not None and local_index < len(selected_member["tracked_track"]["residues"])
                        else "")
            residue = next((item for item in footprints if item["reference_index"] == int(residue_index)), None)
            if residue:
                st.caption(f"Selected anchor residue: {residue['amino_acid']} {residue['author_number']}{residue['insertion_code'] or ''} · mapped Mol* selector {selector or 'unavailable for this component'}")
        else:
            selector = ""
            st.caption("Select a sequence tile or a heatmap cell to carry a residue into Mol*.")
        contact = selected_member["contact"]
        evaluation_id = str(contact.get("evaluation", {}).get("id") or "")
        if st.button("Open selected complex in Scientific Inspector", key="collection-map-open-inspector",
                     disabled=bool(residue_index is not None and not selector)):
            _open_in_inspector(selected_member, selector or "", evaluation_id)

def _database_collection_samples(user: dict) -> list[dict]:
    from mn_molecular_design_studio.molecular_database import MolecularDatabase

    database = MolecularDatabase()
    query = st.text_input("Search Database structures", key="complex-analysis-collection-query",
                          placeholder="Name, sequence, or molecule identifier")
    catalog_key = _widget_key("complex-analysis-collection-catalog-", str(user.get("id") or ""))
    refresh_catalog = st.button("Refresh Database structures", key=_widget_key("complex-analysis-collection-refresh-", str(user.get("id"))))
    catalog = authorized_structure_catalog(database, user, session_state=st.session_state,
        cache_key=catalog_key, query=query, limit=100, refresh=refresh_catalog)
    load_next = st.button("Load next Database structure page",
        key=_widget_key("complex-analysis-collection-next-", f"{user.get('id')}:{query}"),
        disabled=catalog.get("next_offset") is None)
    if load_next:
        catalog = authorized_structure_catalog(database, user, session_state=st.session_state,
            cache_key=catalog_key, query=query, limit=100, advance=True)
    records = list(catalog.get("records") or [])
    records_key = _widget_key("complex-analysis-collection-selected-records-", str(user.get("id") or ""))
    cache_key = _widget_key("complex-analysis-collection-selected-record-cache-", str(user.get("id") or ""))
    structure_key = _widget_key("complex-analysis-collection-selected-structures-", str(user.get("id") or ""))
    cached_records = st.session_state.get(cache_key, {})
    cached_records = dict(cached_records) if isinstance(cached_records, dict) else {}
    selected_before = [str(value) for value in st.session_state.get(records_key, [])]
    # Retain only selected prior-page rows, then refresh the visible page. The
    # later full-record and structure lookups recheck ACLs before any artifact
    # is read, so this cache only keeps picker labels and references available.
    cached_records = {key: cached_records[key] for key in selected_before if key in cached_records}
    for row in records:
        if row.get("id"):
            cached_records[str(row["id"])] = row
    st.session_state[cache_key] = cached_records
    record_by_id = _merge_database_record_options(records, cached_records, selected_before)
    if not records:
        st.info("No accessible Database structures match this search. Previously selected records remain available until removed.")
        if not record_by_id:
            return []
    st.caption(f"Showing {len(records):,} authorized Database structure records"
               + ("; more pages are available." if catalog.get("next_offset") is not None else "."))
    selected_records = st.multiselect("Database records", list(record_by_id),
        format_func=lambda value: f"{record_by_id[value].get('name') or 'Complex'} · {str(value)[-8:]} · {len(record_by_id[value].get('representations') or record_by_id[value].get('structures') or [])} structures",
        key=records_key)
    if not selected_records:
        st.info("Select Database records, then choose the saved structures and models to add to this cohort.")
        return []
    structure_lookup = {}
    for record_id in selected_records:
        record = record_by_id.get(record_id) or {}
        for structure in record.get("representations") or record.get("structures") or []:
            structure_id = str(structure.get("id") or "")
            if structure_id:
                structure_lookup[structure_id] = (record, structure)
    structure_ids = list(structure_lookup)
    if structure_key in st.session_state:
        st.session_state[structure_key] = [value for value in st.session_state[structure_key] if value in structure_lookup]
    selected_structures = st.multiselect("Database structures, models, and poses", structure_ids,
        format_func=lambda value: _structure_label(*structure_lookup[value]),
        key=structure_key)
    if not selected_structures:
        st.info("Choose Database structures to add them to the comparison cohort.")
        return []
    detailed = database.list_molecules(user, molecule_ids=list(dict.fromkeys(selected_records)),
                                       collection="all", limit=1000)["items"]
    samples = _load_selected_structures(database, user, detailed, selected_structures)
    for sample in samples:
        sample.update({"source_kind": "Database", "identity": "database:" + sample["molecule_id"] + ":" + sample["structure_id"],
                       "user": user, "source_reference": "database:" + sample["molecule_id"]})
        try:
            sample["review_annotation"] = database.get_complex_review(
                user, source_kind="database", subject_id=sample["molecule_id"],
                structure_ref=sample["structure_id"], evaluation_id="",
            )
        except (OSError, PermissionError, ValueError, TypeError):
            sample["review_annotation"] = {"status": "Unreviewed", "note": ""}
    return samples


def _results_source_jobs(user: dict) -> list[dict]:
    from mn_molecular_design_studio.complex_analysis_dashboard import _run_catalog
    from mn_molecular_design_studio.run_results import completed_sources

    token = str(st.session_state.get("studio-session-token") or "development-mode")
    return completed_sources(user, token=token, runs=_run_catalog(user))


def _results_structure_choices(user: dict, jobs: list[dict]) -> list[dict]:
    from mn_molecular_design_studio import complex_analysis_ui as workspace, logical_runs, run_results
    from mn_molecular_design_studio.complex_analysis_results_ui import _remote_board_options

    token = str(st.session_state.get("studio-session-token") or "development-mode")
    choices = []
    for job in jobs:
        reference = logical_runs.source_reference(job)
        remote = reference.startswith("cluster:")
        if remote:
            collection = run_results.read_source(job, user, token=token)
            if collection.get("source_errors"):
                st.warning(f"Could not read authorized Results metadata for {job.get('run_id')}: {collection['source_errors'][0]['error']}")
                continue
            options = _remote_board_options(job, collection)
        else:
            options = [item for item in workspace._job_structure_options(user, job) if item.get("candidate_id")]
        for option in options:
            candidate = option.get("candidate") if isinstance(option.get("candidate"), dict) else {}
            candidate_id = str(candidate.get("candidate_id") or option.get("candidate_id") or "")
            if not candidate_id:
                continue
            relative = str(option.get("relative_path") or "")
            raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
            model_value = next((value for value in (candidate.get("model_id"), candidate.get("model_index"),
                raw.get("model_id"), raw.get("model_index")) if value is not None and str(value).strip() != ""), None)
            model = str(model_value) if model_value is not None else "model unspecified"
            identity = "results:" + hashlib.sha256(f"{reference}\0{relative}\0{candidate_id}\0{model}".encode()).hexdigest()[:24]
            source_tool = str(job.get("display_tool") or job.get("tool") or job.get("job_type") or "Results")
            choices.append({"identity": identity, "reference": reference, "job": job, "option": option,
                            "candidate": candidate, "candidate_id": candidate_id, "model_id": model,
                            "relative_path": relative, "remote": remote,
                            "label": f"{source_tool} · {job.get('node_name') or 'This workstation'} · {job.get('run_id')} · {candidate_id} · model {model}"})
    return choices


def _materialize_results_sample(user: dict, choice: dict) -> dict | None:
    from mn_molecular_design_studio import complex_analysis_ui as workspace
    from mn_molecular_design_studio.complex_analysis_results_ui import _materialize_board_structure, board_entries
    from mn_molecular_design_studio.complex_analysis_inspector import protein_sequence_tracks
    from mn_molecular_design_studio.logical_runs import source_reference

    job, option, candidate = choice["job"], choice["option"], choice["candidate"]
    token = str(st.session_state.get("studio-session-token") or "development-mode")
    item = _materialize_board_structure(user, job, option, token=token) if choice["remote"] else option
    if item is None:
        return None
    path = Path(item["path"]).resolve(strict=True)
    content = path.read_bytes()
    digest = hashlib.sha256(content).hexdigest()
    raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
    components = workspace.component_inventory(content, path.name, raw.get("components"))
    chain_ids = sorted({str(row.get("source_chain") or row.get("chain") or "") for row in components
                        if row.get("type") == "protein" and (row.get("source_chain") or row.get("chain"))})
    try:
        tracks = protein_sequence_tracks(path, chain_ids)
    except (ImportError, OSError, TypeError, ValueError):
        tracks = []
    reference = source_reference(job)
    model_id = str(choice["model_id"] if choice["model_id"] != "model unspecified" else "")
    source_tool = str(candidate.get("source_tool") or job.get("display_tool") or job.get("tool") or "Unknown")
    label = choice["label"]
    sample = {
        "source_kind": "Results", "kind": component_kind(components),
        "identity": choice["identity"], "structure_id": choice["identity"],
        "label": label, "path": path, "components": components, "tracks": tracks,
        "context": raw, "candidate": candidate, "candidate_id": choice["candidate_id"],
        "model_id": model_id, "digest": digest, "structure_digest": digest,
        "job": job, "user": user, "source_reference": reference,
        "source_job": f"{job.get('task_group')}/{job.get('run_id')}",
        "source_workstation": str(job.get("node_name") or "This workstation"),
        "remote": bool(choice["remote"]), "relative_path": choice["relative_path"],
        "source_run_dir": str(option.get("source_run_dir") or job.get("run_dir") or ""),
        "structure_field": str(option.get("structure_field") or "complex_pdb"),
        "project_id": str(job.get("project_id") or ""),
        "evaluation_entries": board_entries(candidate, label=label, run_id=str(job.get("run_id") or ""),
                                             digest=digest, job=job),
        "generator_engine": str(raw.get("generator_engine") or raw.get("design_engine")
            or ((raw.get("source_candidate") or {}).get("source_tool") if isinstance(raw.get("source_candidate"), dict) else "")
            or ((raw.get("candidate_set") or {}).get("generator") if isinstance(raw.get("candidate_set"), dict) else "")
            or (job.get("display_tool") or job.get("tool") or "Unknown")),
        "prediction_engine": str(raw.get("prediction_engine") or (source_tool if str(candidate.get("stage") or "").casefold() in {"prediction", "refolding", "structure_prediction"} else "")),
    }
    return sample


def _attach_results_measurements(user: dict, samples: list[dict], jobs: list[dict]) -> None:
    if not samples:
        return
    from mn_molecular_design_studio import logical_runs, run_results
    from mn_molecular_design_studio.complex_analysis_results_ui import _linked_measurements

    token = str(st.session_state.get("studio-session-token") or "development-mode")
    grouped = {}
    for sample in samples:
        grouped.setdefault(sample["source_reference"], []).append(sample)
    by_reference = {logical_runs.source_reference(job): job for job in jobs}
    for reference, selected in grouped.items():
        job = by_reference.get(reference)
        if not job:
            continue
        def load(source):
            return run_results.read_source(source, user, token=token)
        try:
            linked = _linked_measurements(user, job, [
                {"candidate": sample["candidate"], "digest": sample["digest"], "label": sample["label"]}
                for sample in selected
            ], jobs=jobs, load=load)
        except (OSError, PermissionError, ValueError, TypeError):
            linked = []
        for sample in selected:
            candidate_id = sample["candidate_id"]
            model = sample["model_id"]
            sample["evaluation_entries"].extend(row for row in linked
                if str(row.get("source_candidate_id") or "") in {candidate_id, str((sample["candidate"].get("raw_metadata") or {}).get("native_candidate_id") or "")}
                and str(row.get("source_structure_digest") or "") == sample["digest"]
                and (not row.get("source_model_id") or str(row.get("source_model_id")) == model))
            sample["contact_snapshot"] = _results_contact_snapshot(user, sample, jobs)


def _results_contact_snapshot(user: dict, sample: dict, jobs: list[dict]) -> dict:
    """Read a saved Biotite contact artifact for this exact Results model."""
    from mn_molecular_design_studio import complex_analysis_ui as workspace, logical_runs, run_results

    job = sample["job"]
    chains_b = list(sample.get("candidate", {}).get("binder_chains") or [])
    chains_t = list(sample.get("candidate", {}).get("target_chains") or [])
    if not chains_b or not chains_t or set(chains_b) & set(chains_t):
        return {"status": "missing evaluation", "rows": None, "evaluation": {},
                "settings": {}, "selection": {}, "contact_distance": None}
    if not sample.get("remote"):
        try:
            matches = workspace._matching_biotite_jobs(
                user, source_task_group=str(job.get("task_group") or ""),
                source_run_id=str(job.get("run_id") or ""), candidate_id=sample["candidate_id"],
                structure_digest=sample["digest"], binder_chains=chains_b, target_chains=chains_t,
                contact_distance=5.0, model_id=str(sample.get("model_id") or ""))
        except (OSError, PermissionError, ValueError, TypeError):
            matches = []
        for evaluation in matches:
            if str(evaluation.get("status") or "").casefold() not in {"completed", "succeeded"}:
                continue
            item = evaluation.get("selected_item") or {}
            contacts = workspace._read_json(Path(evaluation["run_dir"]) / "artifacts" / "analysis" / str(item.get("id") or "") / "contacts.json")
            if contacts.get("schema") == "studio.contacts.v1" and isinstance(contacts.get("residue_pairs"), list):
                return {"status": "saved", "rows": contacts["residue_pairs"],
                        "evaluation": {"id": workspace._analysis_evaluation_id(Path(evaluation["run_dir"]), user, sample["candidate_id"]),
                                       "run_id": evaluation.get("run_id"), "method": "biotite_interface"},
                        "settings": {"contact_distance": contacts.get("distance_angstrom")},
                        "selection": {"binder_chains": contacts.get("binder_chains"), "target_chains": contacts.get("target_chains")},
                        "contact_distance": contacts.get("distance_angstrom")}
        has_measurement = any("biotite_interface" in str(row.get("method") or row.get("evaluator") or "").casefold()
                              for row in sample.get("evaluation_entries", []))
        failed = any(str(item.get("status") or "").casefold() in {"failed", "error", "cancelled", "canceled"}
                     for item in matches)
        recorded_failure = any("biotite_interface" in str(row.get("method") or row.get("evaluator") or "").casefold()
                              and str(row.get("status") or "").casefold() in {"failed", "error", "cancelled", "canceled"}
                              for row in sample.get("evaluation_entries", []))
        status = ("failed evaluation" if failed or recorded_failure else
                  "artifact unavailable" if has_measurement or matches else "missing evaluation")
        return {"status": status, "rows": None, "evaluation": {}, "settings": {},
                "selection": {}, "contact_distance": None}

    # Remote artifacts use the same Results authorization and per-file artifact endpoint.
    # The small batch manifest identifies the exact generated item; only then is contacts.json requested.
    matching_jobs = []
    source_node = str(job.get("node_id") or "")
    for analysis in jobs:
        if str(analysis.get("task_group") or "") != "analysis":
            continue
        if "biotite" not in str(analysis.get("tool") or analysis.get("display_tool") or "").casefold():
            continue
        if (str(analysis.get("node_id") or "") != source_node
                or str(analysis.get("upstream_task_group") or "") != str(job.get("task_group") or "")
                or str(analysis.get("upstream_run_id") or "") != str(job.get("run_id") or "")):
            continue
        matching_jobs.append(analysis)
    token = str(st.session_state.get("studio-session-token") or "development-mode")
    for analysis in matching_jobs:
        try:
            collection = run_results.read_source(analysis, user, token=token)
            if collection.get("source_errors") or not collection.get("candidates"):
                continue
            outputs = {str(row.get("candidate_id") or ""): row for row in collection["candidates"]}
            from mn_molecular_design_studio.cluster import ui as cluster_ui
            manifest_bytes = cluster_ui._fetch_remote_artifact_bytes(
                str(analysis.get("cluster_job_id") or ""), "artifacts/analysis/batch_manifest.json", token,
                "cohort-biotite-manifest-" + hashlib.sha256(str(analysis.get("cluster_job_id")).encode()).hexdigest()[:16],
                label="saved Biotite batch manifest")
            if not manifest_bytes:
                continue
            manifest = json.loads(manifest_bytes)
            settings = manifest.get("settings") if isinstance(manifest.get("settings"), dict) else {}
            if float(settings.get("contact_distance", 5.0)) != 5.0:
                continue
            selected_item = None
            output_candidate = None
            matched_receipt = {}
            for item in manifest.get("items") or []:
                matched, receipt = workspace._biotite_item_matches_source(
                    item, outputs, candidate_id=sample["candidate_id"], structure_digest=sample["digest"],
                    binder_chains=chains_b, target_chains=chains_t, source_run_id=str(job.get("run_id") or ""))
                if not matched:
                    continue
                output_candidate = outputs.get(str(item.get("candidate_id") or ""), {})
                saved_model = str(item.get("model_id") or receipt.get("model_id")
                                  or (manifest.get("selection") or {}).get("model_id") or "")
                if saved_model != str(sample.get("model_id") or ""):
                    continue
                selected_item, matched_receipt = item, receipt
                break
            if selected_item is None:
                continue
            contact_path = f"artifacts/analysis/{selected_item.get('id')}/contacts.json"
            contact_bytes = cluster_ui._fetch_remote_artifact_bytes(
                str(analysis.get("cluster_job_id") or ""), contact_path, token,
                "cohort-biotite-contacts-" + hashlib.sha256(f"{analysis.get('cluster_job_id')}:{contact_path}".encode()).hexdigest()[:16],
                label="saved Biotite contact artifact")
            if not contact_bytes:
                continue
            contacts = json.loads(contact_bytes)
            if (contacts.get("schema") != "studio.contacts.v1"
                    or sorted(map(str, contacts.get("binder_chains") or [])) != sorted(map(str, chains_b))
                    or sorted(map(str, contacts.get("target_chains") or [])) != sorted(map(str, chains_t))):
                continue
            from mn_molecular_design_studio.evaluations import results_evaluation_descriptor
            descriptor = results_evaluation_descriptor(output_candidate or {}, analysis, structure_digest=sample["digest"])
            return {"status": "saved", "rows": contacts.get("residue_pairs") or [],
                    "evaluation": descriptor, "settings": {"contact_distance": contacts.get("distance_angstrom")},
                    "selection": {"binder_chains": chains_b, "target_chains": chains_t,
                                  "model_id": sample.get("model_id"), "source_receipt": matched_receipt},
                    "contact_distance": contacts.get("distance_angstrom")}
        except (OSError, PermissionError, ValueError, TypeError, KeyError, json.JSONDecodeError):
            continue
    has_measurement = any("biotite_interface" in str(row.get("method") or row.get("evaluator") or "").casefold()
                          for row in sample.get("evaluation_entries", []))
    failed = any(str(item.get("status") or "").casefold() in {"failed", "error", "cancelled", "canceled"}
                 for item in matching_jobs)
    recorded_failure = any("biotite_interface" in str(row.get("method") or row.get("evaluator") or "").casefold()
                          and str(row.get("status") or "").casefold() in {"failed", "error", "cancelled", "canceled"}
                          for row in sample.get("evaluation_entries", []))
    status = ("failed evaluation" if failed or recorded_failure else
              "artifact unavailable" if has_measurement or matching_jobs else "missing evaluation")
    return {"status": status, "rows": None, "evaluation": {}, "settings": {},
            "selection": {}, "contact_distance": None}


def _cohort_inventory_row(sample: dict) -> dict:
    candidate = sample.get("candidate") if isinstance(sample.get("candidate"), dict) else {}
    raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else sample.get("context", {})
    from mn_molecular_design_studio.native_filter_metrics import native_filter_evidence
    from mn_molecular_design_studio.complex_analysis_collection import acceptance_decisions

    native = native_filter_evidence(candidate) if candidate else {"status": "not recorded"}
    studio = raw.get("validation_acceptance") if isinstance(raw.get("validation_acceptance"), dict) else {}
    review = raw.get("human_review") if isinstance(raw.get("human_review"), dict) else {}
    decisions = acceptance_decisions(candidate)
    final = next((row for row in decisions if row.get("stage_id") == "final"), None)
    stage_text = []
    for decision in decisions:
        policy = decision.get("policy") if isinstance(decision.get("policy"), dict) else {}
        policy_name = policy.get("schema") or policy.get("version") or decision.get("source") or "Policy not recorded"
        rules = policy.get("mandatory_filters") or policy.get("rules") or []
        thresholds = "; ".join(
            f"{rule.get('field')} {rule.get('operator')} {rule.get('value')}"
            for rule in rules if isinstance(rule, dict) and rule.get("field")
        )
        scope = f"model {decision['model_id']}" if decision.get("model_id") else f"candidate {decision.get('candidate_id') or 'not recorded'}"
        stage_text.append(f"{decision['stage']}: {decision['status']} · {scope} · {policy_name}"
                          + (f" · {thresholds}" if thresholds else ""))
    experimental = raw.get("experimental_validation") if isinstance(raw.get("experimental_validation"), dict) else {}
    return {
        "Candidate / molecule": str(candidate.get("candidate_id") or sample.get("molecule_id") or sample.get("structure_id") or "Unknown"),
        "Original candidate": str((raw.get("source_candidate") or {}).get("candidate_id") or raw.get("native_candidate_id") or candidate.get("candidate_id") or "Not recorded"),
        "Generator": sample.get("generator_engine") or candidate.get("source_tool") or sample.get("source_kind"),
        "Prediction engine": sample.get("prediction_engine") or "Not recorded",
        "Target": raw.get("target_identity") or raw.get("target_id") or sample.get("target_identity") or _recorded_target_sequence_identity(sample) or "Not recorded",
        "Binder sequence": candidate.get("binder_sequence") or "Not recorded",
        "Model / pose": sample.get("model_id") or "Not recorded",
        "Molecular type": sample.get("kind") or "Not classified",
        "Native acceptance": native.get("status") or "not recorded",
        "Studio filter": studio.get("status") or "not evaluated",
        "Final computational acceptance": final.get("status") if final else "Not recorded",
        "Recorded acceptance evidence": " | ".join(stage_text) or "No explicit decision recorded",
        "Human review": sample.get("review_annotation", {}).get("status")
                        or review.get("status") or ("selected" if review.get("selected") else "Not recorded"),
        "Manual selection": bool(review.get("selected", False)),
        "Experimental validation": experimental.get("status") or "Not recorded",
        "Evaluation coverage": ", ".join(sorted({str(row.get("method") or row.get("evaluator") or "Unknown") for row in sample.get("evaluation_entries", [])})) or "No saved Evaluation",
        "Source job": sample.get("source_job") or sample.get("record", {}).get("source_job") or "Database",
        "Workstation": sample.get("source_workstation") or "Database",
        "Structure SHA-256": sample.get("digest") or "Not recorded",
        "Stored measurements": json.dumps(candidate.get("metrics") or {}, sort_keys=True, default=str),
        "Row identity": sample.get("identity") or sample.get("structure_id"),
    }


def _recorded_target_sequence_identity(sample: dict) -> str:
    """Name an exact saved target construct only when its target-chain roles are recorded."""
    candidate = sample.get("candidate") if isinstance(sample.get("candidate"), dict) else {}
    raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
    context = sample.get("context") if isinstance(sample.get("context"), dict) else {}
    target_chains = candidate.get("target_chains") or raw.get("target_chains") or context.get("target_chains")
    if not isinstance(target_chains, (list, tuple)) or not target_chains:
        return ""
    by_chain = {str(track.get("chain_id") or ""): str(track.get("sequence") or "")
                for track in sample.get("tracks", []) if isinstance(track, dict)}
    sequences = [by_chain.get(str(chain), "") for chain in target_chains]
    if not sequences or any(not sequence for sequence in sequences):
        return ""
    joined = "|".join(sequences)
    digest = hashlib.sha256(joined.encode()).hexdigest()[:12]
    chain_count = len(sequences)
    return f"Target sequence · {sum(map(len, sequences))} aa · {chain_count} chain(s) · {digest}"


def _comparison_sample(user: dict, sample: dict) -> dict:
    from mn_molecular_design_studio.complex_analysis_comparison import protein_chain_numbering, target_structure_fingerprint

    candidate = sample.get("candidate") if isinstance(sample.get("candidate"), dict) else {}
    raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
    context = sample.get("context") if isinstance(sample.get("context"), dict) else {}
    target = [str(value) for value in (candidate.get("target_chains") or raw.get("target_chains")
             or context.get("target_chains") or context.get("receptor_chains") or [])]
    binder = [str(value) for value in (candidate.get("binder_chains") or raw.get("binder_chains")
             or context.get("binder_chains") or [])]
    target_numbering = protein_chain_numbering(sample["path"], target) if target else []
    binder_numbering = protein_chain_numbering(sample["path"], binder) if binder else []
    is_results = sample.get("source_kind") == "Results"
    record, structure = sample.get("record") or {}, sample.get("structure") or {}
    comparison = {
        "kind": "results" if is_results else "database", "identity": sample["identity"],
        "label": sample["label"], "path": sample["path"], "context": context,
        "candidate": candidate, "job": sample.get("job") or {},
        "record": record, "structure": structure,
        "candidate_id": str(sample.get("candidate_id") or ""), "project_id": str(sample.get("project_id") or ""),
        "source_ref": str(sample.get("source_reference") or ""), "relative_path": str(sample.get("relative_path") or ""),
        "components": sample.get("components") or [], "component_metadata": sample.get("components") or [],
        "model_id": str(sample.get("model_id") or ""), "structure_digest": str(sample.get("digest") or ""),
        "target_chains": target, "binder_chains": binder,
        "target_numbering": target_numbering, "binder_numbering": binder_numbering,
        "target_chains_sequence": [row["sequence"] for row in target_numbering],
        "binder_chains_sequence": [row["sequence"] for row in binder_numbering],
        "target_structure_fingerprint": target_structure_fingerprint(sample["path"], target) if target else "",
        "complex_type": "protein-protein" if sample.get("kind") == "protein-protein" else "protein-ligand" if sample.get("kind") == "protein-ligand" else sample.get("kind"),
        "format": sample["path"].suffix.casefold(), "viewer_key": "complex-analysis-cohort-" + hashlib.sha256(sample["identity"].encode()).hexdigest()[:14],
        "evaluation_entries": sample.get("evaluation_entries") or [],
        "saved_biotite_contacts": ({"contacts": sample["contact_snapshot"].get("rows") or [],
            "evaluation": sample["contact_snapshot"].get("evaluation") or {}}
            if sample.get("contact_snapshot", {}).get("status") == "saved" else None),
    }
    return comparison


def render_collection_workspace(user: dict, run_rows: list[dict] | None = None, *,
                                workflow_schema=None, submit=None) -> None:
    """Render the mixed Results/Database design explorer inside Complex Analysis."""
    st.subheader("Unified Design Explorer · Evaluation Board")
    st.caption("Select authorized saved structures from Results and/or the curated Database. By default only explicitly accepted final computational decisions enter the cohort. Native passes, earlier-stage passes, human review and experimental validation remain distinct; no universal score or Database promotion is applied.")
    sources = st.multiselect("Include sources", ["Results", "Database"], default=["Database"],
                             key="complex-analysis-cohort-source-types")
    samples = []
    jobs = []
    if "Database" in sources:
        samples.extend(_database_collection_samples(user))
    if "Results" in sources:
        try:
            jobs = _results_source_jobs(user)
        except (OSError, RuntimeError, ValueError, PermissionError) as exc:
            st.warning(f"Could not load the authorized Results catalog: {exc}")
            jobs = []
        if not jobs:
            st.info("No authorized completed Results collections are available.")
        else:
            from mn_molecular_design_studio.logical_runs import source_reference
            by_ref = {source_reference(job): job for job in jobs}
            references = list(by_ref)
            run_key = "complex-analysis-cohort-result-jobs"
            if run_key in st.session_state:
                st.session_state[run_key] = [value for value in st.session_state[run_key] if value in references]
            selected_refs = st.multiselect("Results jobs (select one or more campaigns, engines or workstations)", references,
                format_func=lambda value: f"{by_ref[value].get('display_tool') or by_ref[value].get('tool') or 'Result'} · {by_ref[value].get('node_name') or 'This workstation'} · {by_ref[value].get('task_group')}/{by_ref[value].get('run_id')}",
                key=run_key)
            selected_jobs = [by_ref[value] for value in selected_refs]
            choices = _results_structure_choices(user, selected_jobs) if selected_jobs else []
            if selected_jobs and not choices:
                st.info("The selected Results jobs have no accessible candidate-associated structure models.")
            elif choices:
                choice_by_id = {row["identity"]: row for row in choices}
                result_key = "complex-analysis-cohort-result-structures"
                if result_key in st.session_state:
                    st.session_state[result_key] = [value for value in st.session_state[result_key] if value in choice_by_id]
                selected_ids = st.multiselect("Candidate structures / prediction models", list(choice_by_id),
                    format_func=lambda value: choice_by_id[value]["label"], key=result_key)
                for identity in selected_ids:
                    try:
                        sample = _materialize_results_sample(user, choice_by_id[identity])
                    except (OSError, PermissionError, ValueError, TypeError) as exc:
                        st.warning(f"Could not open this authorized Results structure: {exc}")
                        continue
                    if sample is not None:
                        samples.append(sample)
                _attach_results_measurements(user, [item for item in samples if item.get("source_kind") == "Results"], jobs)
    if not samples:
        st.info("Choose saved Database structures and/or Results candidate models to populate the Design Explorer.")
        return

    from mn_molecular_design_studio.complex_analysis_collection import (
        acceptance_decisions, acceptance_stage_options, filter_acceptance_cohort,
    )
    st.caption(f"Loaded {len(samples)} selected structure/model row(s). Results structures are retrieved only after selection; original Results remain unmodified.")
    stages = acceptance_stage_options(samples)
    final_rows = [row for sample in samples for row in acceptance_decisions(sample.get("candidate") or {})
                  if row.get("stage_id") == "final"]
    cohort_mode = st.selectbox(
        "Acceptance cohort",
        ["Final accepted designs", "Passed at selected stage", "Rejected at selected stage",
         "Pending / unresolved at selected stage", "No recorded decision", "All designs"],
        key="complex-analysis-cohort-acceptance-mode",
        help="The default includes only an explicit final computational pass attached to this candidate/model. Select a stage to inspect its own decisions; a native generator pass is never treated as final acceptance.",
    )
    stage_id = ""
    if cohort_mode in {"Passed at selected stage", "Rejected at selected stage", "Pending / unresolved at selected stage"}:
        if stages:
            stage_ids = [item[0] for item in stages]
            stage_key = "complex-analysis-cohort-acceptance-stage"
            if st.session_state.get(stage_key) not in stage_ids:
                st.session_state[stage_key] = stage_ids[0]
            stage_id = st.selectbox("Recorded acceptance stage", stage_ids,
                format_func=lambda value: next(label for key, label in stages if key == value), key=stage_key)
            selected_stage = next(label for key, label in stages if key == stage_id)
            st.caption(f"Stage decisions retain their recorded source and policy per candidate. Selected stage: {selected_stage}.")
        else:
            st.info("No recorded acceptance stages are available for these structures. Choose All designs or No recorded decision to inspect them without treating them as passed.")
    cohort_samples = filter_acceptance_cohort(samples, cohort_mode, stage_id)
    if cohort_mode == "Final accepted designs" and not cohort_samples:
        if final_rows:
            st.info("No selected candidate/model has a recorded final computational pass. Rejected, unresolved and other stages remain separate; choose a stage explicitly to inspect another cohort.")
        else:
            st.info("No explicit final computational decisions are recorded for these structures. Choose a recorded stage to view its passed cohort, or deliberately select All designs or No recorded decision.")
    st.caption(f"Acceptance cohort: {len(cohort_samples)} of {len(samples)} selected structure/model row(s). Human review and experimental validation do not change computational acceptance.")
    if not cohort_samples:
        return
    samples = cohort_samples
    import pandas as pd
    inventory = [_cohort_inventory_row(sample) for sample in samples]
    frame = pd.DataFrame(inventory)
    search = st.text_input("Search candidates, targets, sequences, jobs, or workstations", key="complex-analysis-cohort-search").casefold().strip()
    engine_options = sorted({str(row["Generator"]) for row in inventory if row.get("Generator")})
    engine_filter = st.multiselect("Generator filter", engine_options, default=engine_options,
                                   key="complex-analysis-cohort-generator-filter")
    prediction_options = sorted({str(row["Prediction engine"]) for row in inventory if row.get("Prediction engine") and row["Prediction engine"] != "Not recorded"})
    prediction_filter = st.multiselect("Prediction engine filter", prediction_options, default=prediction_options,
                                       key="complex-analysis-cohort-prediction-filter")
    visible = [row for row in inventory if row.get("Generator") in engine_filter
        and (row.get("Prediction engine") == "Not recorded" or row.get("Prediction engine") in prediction_filter)
        and (not search or search in json.dumps(row, sort_keys=True, default=str).casefold())]
    if not visible:
        st.info("No candidate structures match the current filters.")
        return
    visible_ids = [str(row["Row identity"]) for row in visible]
    sample_by_id = {str(sample.get("identity") or sample.get("structure_id")): sample for sample in samples}
    st.dataframe(pd.DataFrame([{key: value for key, value in row.items() if key != "Row identity"} for row in visible]),
                 hide_index=True, use_container_width=True)
    st.caption("Per-engine measurements retain their stored method, model and original units. Missing, not evaluated, and measured-zero values are not combined.")
    select_key = "complex-analysis-cohort-selected"
    if select_key in st.session_state:
        st.session_state[select_key] = [value for value in st.session_state[select_key] if value in visible_ids]
    selected_ids = st.multiselect("Select candidate/model rows for residue map, comparison, or Inspector", visible_ids,
        format_func=lambda value: sample_by_id[value]["label"], key=select_key)
    selected_samples = [sample_by_id[value] for value in selected_ids if value in sample_by_id]
    if not selected_samples:
        return
    _render_batch_assessment(user, selected_samples, workflow_schema=workflow_schema, submit=submit)
    filtered_evaluations = _render_evaluation_board(selected_samples)
    from mn_molecular_design_studio.complex_analysis_shortlist_export import render_shortlist_review_export
    render_shortlist_review_export(user, selected_samples, key_prefix="explorer-shortlist")
    kinds = {sample.get("kind") for sample in selected_samples}
    if "protein-ligand" in kinds or "protein-protein-and-ligand" in kinds:
        ligand_samples = [sample for sample in selected_samples if sample.get("kind") in {"protein-ligand", "protein-protein-and-ligand"}]
        st.markdown("#### Protein–ligand pose group")
        st.caption("Ligand chemistry, receptor and pose remain distinct from protein–protein contacts; no PPI contact frequency is applied.")
        st.dataframe(pd.DataFrame(_protein_ligand_rows(ligand_samples)), hide_index=True, use_container_width=True)
    ppi_samples = [sample for sample in selected_samples if sample.get("kind") == "protein-protein"]
    if ppi_samples:
        _render_protein_protein_map(ppi_samples)
    other = [sample for sample in selected_samples if sample.get("kind") not in {"protein-protein", "protein-ligand"}]
    if other:
        st.markdown("#### Structures kept outside current contact maps")
        st.dataframe([{"Structure": sample["label"], "Molecular type": sample.get("kind"),
                       "Reason": "Current residue contact map applies only to protein–protein complexes; this representation remains visible in the candidate board and Inspector."}
                      for sample in other], hide_index=True, use_container_width=True)

    inspect_id = st.selectbox("Candidate/model to open in Scientific Inspector", selected_ids,
        format_func=lambda value: sample_by_id[value]["label"], key="complex-analysis-cohort-inspector")
    evaluation_id = st.selectbox("Evaluation context (optional)", [""] + [str(row.get("evaluation_id") or "") for row in filtered_evaluations if row.get("evaluation_id")],
        format_func=lambda value: "No specific Evaluation" if not value else next((f"{row.get('method')} · {row.get('structure_label')}" for row in filtered_evaluations if str(row.get("evaluation_id") or "") == value), value),
        key="complex-analysis-cohort-inspector-evaluation")
    if st.button("Open selected model in Scientific Inspector", key="complex-analysis-cohort-open-inspector"):
        sample = sample_by_id[inspect_id]
        if sample.get("source_kind") == "Results":
            job = sample["job"]
            if sample.get("remote"):
                from mn_molecular_design_studio import complex_analysis_ui as workspace
                workspace.open_remote_result_complex(job, sample["candidate"], sample["relative_path"], user=user,
                    structure_content=Path(sample["path"]).read_bytes(), source_run_dir=sample["source_run_dir"],
                    structure_field=sample["structure_field"], structure_sha256=sample["digest"])
            else:
                from mn_molecular_design_studio import complex_analysis_ui as workspace
                workspace.open_result_complex(job["run_dir"], sample["path"], sample["candidate"], user=user)
        else:
            _open_in_inspector(sample, "", evaluation_id)

    compare_ids = st.multiselect("Choose 2–5 structures for detailed comparison", selected_ids,
        default=selected_ids[:2] if len(selected_ids) >= 2 else [], max_selections=5,
        format_func=lambda value: sample_by_id[value]["label"], key="complex-analysis-cohort-compare")
    if len(compare_ids) >= 2:
        with st.expander("Detailed structural comparison", expanded=False):
            compare_sources = [_comparison_sample(user, sample_by_id[value]) for value in compare_ids]
            missing_roles = [sample["label"] for sample in compare_sources
                if (sample.get("complex_type") == "protein-protein"
                    and (not sample.get("target_chains") or not sample.get("binder_chains")
                         or set(sample.get("target_chains") or []) & set(sample.get("binder_chains") or [])))
                or (sample.get("complex_type") == "protein-ligand" and not sample.get("target_chains"))]
            if missing_roles:
                st.info("Comparison needs recorded target/binder or receptor chain assignments. Roles are missing for: "
                        + ", ".join(missing_roles) + ". The structures remain available in the board and Inspector.")
            else:
                from mn_molecular_design_studio.complex_analysis_comparison_ui import render_comparison
                render_comparison(user, selected_sources=compare_sources)
