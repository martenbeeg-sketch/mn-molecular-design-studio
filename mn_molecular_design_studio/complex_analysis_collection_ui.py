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
        # `records` came from the ACL-filtered Database listing for this user.
        # Avoid a second Database-wide index rebuild per selected representation;
        # structure_file below independently rechecks access to the exact file.
        record = owner
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
        samples.append({
            "identity": str(structure_id),
            "molecule_id": str(record["id"]),
            "record": record,
            "structure": structure,
            "structure_id": str(structure_id),
            "label": _structure_label(record, structure),
            "path": Path(path),
            "components": components,
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


def _contact_snapshot(sample: dict) -> dict:
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
    entries = [row for sample in samples for row in sample.get("evaluation_entries", [])]
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

    display_rows = [{
        "Structure": row["structure_label"], "Kind": row["category"], "Method": row["method"],
        "Source engine": row["source_engine"], "Status": row["status"], "Created": row["created_at"],
        "Evaluation ID": row["evaluation_id"], "Run ID": row["run_id"] or "Not recorded",
        "Candidate / pose ID": row["candidate_id"] or "Not recorded",
        "Structure digest": row["structure_digest"] or "Not recorded",
        "Parameters": json.dumps(row["parameters"], sort_keys=True, default=str),
        "Recorded measurements": json.dumps(row["metrics"], sort_keys=True, default=str),
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


def _open_in_inspector(sample: dict, selector: str = "", evaluation_id: str = "") -> None:
    from mn_molecular_design_studio import complex_analysis_ui as workspace

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
        st.info("Select at least two protein–protein Database structures to build a collection map.")
        return
    sample_ids = [sample["structure_id"] for sample in samples]
    sample_labels = {sample["structure_id"]: sample["label"] for sample in samples}
    anchor_id = st.selectbox("Residue correspondence anchor", sample_ids,
                             format_func=lambda value: sample_labels.get(value, value),
                             key="collection-map-anchor")
    anchor = next(sample for sample in samples if sample["structure_id"] == anchor_id)
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
            if sample["structure_id"] == anchor_id:
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
            if sample["structure_id"] == anchor_id:
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
    profile_options = ["All compatible saved maps", *profiles]
    profile = st.selectbox("Saved Biotite map settings", profile_options,
                           key="collection-map-contact-profile")
    map_members = []
    for member in tracked:
        contact = member["contact"]
        distance = contact.get("contact_distance")
        label = _contact_profile_label(distance)
        if profile != "All compatible saved maps" and contact["status"] == "saved" and label != profile:
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
    if profile == "All compatible saved maps" and len(saved_profiles) > 1:
        st.warning("Saved contact maps use more than one recorded contact-distance setting. The table shows each saved map, but the pooled frequency denominator mixes those settings; choose a single setting above for a comparable frequency map.")
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
    st.caption("Legend: ● saved contact · · saved map with no contact · ? no matching saved Evaluation · ! Evaluation artifact unavailable · – saved map uses another selected setting.")
    selection_state = getattr(grid, "selection", {}) or {}
    cells = selection_state.get("cells", []) if hasattr(selection_state, "get") else []
    visible_ids = [member["structure_id"] for member in visible_members]
    selected_cell = heatmap_selection(cells, visible_ids, residue_columns)

    st.markdown("##### Open a selected complex in the Scientific Inspector")
    mapped_members = [member for member in tracked if member["correspondence"].get("compatible")]
    mapped_ids = [member["structure_id"] for member in mapped_members]
    mapped_by_id = {member["structure_id"]: member for member in mapped_members}
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

def render_collection_workspace(user: dict) -> None:
    """Render a scalable, ACL-scoped collection browser over saved Database data."""
    from mn_molecular_design_studio.molecular_database import MolecularDatabase

    st.subheader("Collection Map & Evaluation Board")
    st.caption("Explore saved Database structures as a collection. Component correspondence is explicit and sequence-based; unsupported structures stay separate with a reason.")
    database = MolecularDatabase()
    query = st.text_input("Search Database structures", key="complex-analysis-collection-query",
                          placeholder="Name, sequence, or molecule identifier")
    catalog_key = _widget_key("complex-analysis-collection-catalog-", str(user.get("id") or ""))
    refresh_catalog = st.button("Refresh Database structures", key=_widget_key("complex-analysis-collection-refresh-", str(user.get("id"))))
    catalog = authorized_structure_catalog(
        database, user, session_state=st.session_state, cache_key=catalog_key,
        query=query, limit=100, refresh=refresh_catalog,
    )
    load_next = st.button(
        "Load next Database structure page",
        key=_widget_key("complex-analysis-collection-next-", f"{user.get('id')}:{query}"),
        disabled=catalog.get("next_offset") is None,
    )
    if load_next:
        catalog = authorized_structure_catalog(
            database, user, session_state=st.session_state, cache_key=catalog_key,
            query=query, limit=100, advance=True,
        )
    records = list(catalog.get("records") or [])
    if not records:
        st.info("No accessible Database structures match this search.")
        return
    st.caption(f"Showing {len(records):,} authorized structure-bearing record choices"
               + ("; more pages are available." if catalog.get("next_offset") is not None else "."))
    record_by_id = {str(row["id"]): row for row in records}
    record_ids = list(record_by_id)
    selection_suffix = hashlib.sha256(query.encode()).hexdigest()[:10]
    selection_key = f"complex-analysis-collection-records-{selection_suffix}"
    if selection_key in st.session_state:
        st.session_state[selection_key] = [value for value in st.session_state[selection_key] if value in record_ids]
    selected_records = st.multiselect(
        "Database records", record_ids,
        format_func=lambda value: f"{record_by_id[value].get('name') or 'Complex'} · {str(value)[-8:]} · {len(record_by_id[value].get('representations') or record_by_id[value].get('structures') or [])} structures",
        key=selection_key,
    )
    if not selected_records:
        st.info("Select one or more accessible Database records, then choose their saved structures and models.")
        return
    structure_records = [record_by_id[value] for value in selected_records if value in record_by_id]
    structure_lookup = {}
    for record in structure_records:
        for structure in record.get("representations") or record.get("structures") or []:
            structure_id = str(structure.get("id") or "")
            if structure_id:
                structure_lookup[structure_id] = (record, structure)
    structure_ids = list(structure_lookup)
    selected_structures = st.multiselect(
        "Saved structures, models, and poses", structure_ids,
        format_func=lambda value: _structure_label(*structure_lookup[value]),
        key="complex-analysis-collection-structures",
    )
    if not selected_structures:
        st.info("Choose saved representations to populate the map and Evaluation Board. There is no 2–5 structure limit.")
        return
    if len(selected_structures) > 1000:
        st.warning("The interactive map displays large collections in pages; select at most 1,000 structures per workspace view.")
        selected_structures = selected_structures[:1000]
    with st.spinner("Loading authorized Database representations and their saved Evaluations…"):
        detailed = database.list_molecules(
            user, molecule_ids=list(dict.fromkeys(str(record_id) for record_id in selected_records)),
            collection="all", limit=1000,
        )["items"]
        samples = _load_selected_structures(database, user, detailed, selected_structures)
    if not samples:
        st.warning("The selected structures are unavailable or no longer accessible to this account.")
        return
    st.caption(f"Loaded {len(samples)} authorized structure representation(s) from {len({item['molecule_id'] for item in samples})} Database record(s).")
    filtered = _render_evaluation_board(samples)
    st.markdown("#### Molecular component groups")
    kinds = {sample["kind"] for sample in samples}
    st.dataframe([{"Structure": sample["label"], "Molecular type": sample["kind"],
                   "Protein components": ", ".join(f"{row['chain_id']} ({len(row['sequence'])} aa)" for row in sample["tracks"]) or "None",
                   "Saved Evaluation methods": ", ".join(sorted({str(row.get('method')) for row in sample["evaluation_entries"]})) or "None",
                   "Structure ID": sample["structure_id"], "SHA256": sample["digest"]}
                  for sample in samples], hide_index=True, use_container_width=True)
    st.caption("Database identities, structure IDs, exact structure digests, Evaluation IDs, source run IDs, recorded settings, and raw metric fields stay attached to the rows above.")

    if "protein-ligand" in kinds or "protein-protein-and-ligand" in kinds:
        ligand_samples = [sample for sample in samples if sample["kind"] in {"protein-ligand", "protein-protein-and-ligand"}]
        if ligand_samples:
            st.markdown("#### Protein–ligand pose groups")
            st.caption("Protein–ligand grouping preserves the recorded isomeric ligand identity, receptor sequence, pose residue and structure identity. No protein–protein contact-frequency metric is applied here.")
            import pandas as pd
            st.dataframe(pd.DataFrame(_protein_ligand_rows(ligand_samples)), hide_index=True, use_container_width=True)

    ppi_samples = [sample for sample in samples if sample["kind"] == "protein-protein"]
    if ppi_samples:
        _render_protein_protein_map(ppi_samples)
    other = [sample for sample in samples if sample["kind"] not in {"protein-protein", "protein-ligand"}]
    if other:
        st.markdown("#### Kept in separate groups")
        st.dataframe([{"Structure": sample["label"], "Reason": (
            "Protein–protein with ligand or mixed components; map a saved evaluator to a specific component pair before using PPI-only contact summaries."
            if sample["kind"] == "protein-protein-and-ligand" else
            "This structure does not meet the protein–protein or protein–ligand collection criteria for the current map views."),
            "Molecular type": sample["kind"]} for sample in other], hide_index=True, use_container_width=True)

    if filtered:
        selected_eval_id = str(st.session_state.get("collection-eval-detail") or "")
        chosen = next((row for row in filtered if row["evaluation_id"] == selected_eval_id), None)
        if chosen:
            sample = next((item for item in samples if item["structure_id"] == chosen["structure_id"]), None)
            if sample and st.button("Open Evaluation structure in Scientific Inspector", key="collection-board-open-inspector"):
                _open_in_inspector(sample, "", chosen["evaluation_id"])
