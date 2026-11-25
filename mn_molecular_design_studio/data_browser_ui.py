"""The Database workspace: scientific folders, saved results and linked details."""
from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path

import streamlit as st

from mn_molecular_design_studio.data_browser import (
    ReadOnlyLibrary, admission_measurements, concise_result_label, expand_scope, listed_results, load_saved_record,
    navigation_tree, producing_evaluations, saved_receptor_file, scope_path, source_collection_scope, dataset_results_scope,
)
from mn_molecular_design_studio.database_browser_cache import BrowserCache
from mn_molecular_design_studio.database_browser_labels import calculation_label, run_label, tool_label
from mn_molecular_design_studio.database_design_evaluations import analysis_results, connected_evaluations, design_subject
from mn_molecular_design_studio.folder_tree import render_tree
from mn_molecular_design_studio.scientific_browser_tree import scientific_branches, scientific_scopes, scope_counts, source_jobs, target_constructs
from mn_molecular_design_studio.data_browser_table import (
    COLUMN_PRESETS, column_label, default_columns, evaluation_entries, result_table, table_columns, table_csv, table_profile,
)
from mn_molecular_design_studio.binderscope_table import NUMBER_FORMATS, PREFIX, VIEWER_COLUMN, inspector_choices, inspector_label
from mn_molecular_design_studio.database_actions_ui import read_only_workspace, render_library_actions, render_selection_actions
from mn_molecular_design_studio.navigation import database_destination


def _cache(database):
    cached = st.session_state.get("data-browser-cache")
    if not isinstance(cached, BrowserCache) or cached.path != database.path.resolve():
        if isinstance(cached, BrowserCache):
            cached.close()
        cached = BrowserCache(database.path)
        st.session_state["data-browser-cache"] = cached
    return cached


def _clear_result_search():
    st.session_state["data-browser-result-query"] = ""


def _reset_folder_refinements(scope_id, selected_id):
    _clear_result_search()
    st.session_state["data-browser-type"] = ""
    st.session_state["data-browser-coordinates"] = False
    st.session_state["data-browser-activity-" + scope_id] = ""
    st.session_state["data-browser-construct-" + selected_id] = ""


def _change_workspace():
    layout = st.session_state.get("data-browser-layout") or "Split"
    st.session_state["data-browser-layout"] = layout
    st.session_state["data-browser-pin-folders"] = layout == "Split"


def _workspace():
    controls, pin_control, folder_control = st.columns([3, 1, 1], vertical_alignment="bottom")
    with controls:
        layout = st.segmented_control("Workspace", ["Split", "Table", "Structure"], default="Split",
            key="data-browser-layout", on_change=_change_workspace,
            help="Split keeps the original layout. Table or Structure expands that pane without changing the selected data.") or "Split"
    with pin_control:
        st.session_state.setdefault("data-browser-pin-folders", True)
        pinned = st.toggle("Pin folders", key="data-browser-pin-folders",
            help="Keep the folder tree beside the content, or unpin it to open it as an overlay.")
    weights = {"folders": .27, "table": .43, "details": .30}
    if layout != "Split":
        weights.update(table=.73, details=.73)
    names = (["folders"] if pinned else []) + (["table", "details"] if layout == "Split" else
        ["table"] if layout == "Table" else ["details"])
    containers = st.columns([weights[name] for name in names], gap="medium") if len(names) > 1 else [st.container()]
    panes = dict(zip(names, containers))
    if pinned:
        navigation = panes["folders"]
    else:
        with folder_control:
            navigation = st.popover("Folders", icon=":material/folder_open:", width="stretch",
                help="Browse or refine folders without leaving this focused workspace.")
    return layout, navigation, panes.get("table", panes.get("details")), panes.get("details")


def _page_selection(signature, page, visible):
    context = (signature, page)
    selection = st.session_state.get("data-browser-page-selection")
    if not selection or selection["context"] != context:
        selection = {"context": context, "selected": [], "focused": visible[0]["id"] if visible else ""}
    identifiers = {entry["id"] for entry in visible}
    selection["selected"] = [identifier for identifier in selection["selected"] if identifier in identifiers]
    if selection["focused"] not in identifiers:
        selection["focused"] = visible[0]["id"] if visible else ""
    st.session_state["data-browser-page-selection"] = selection
    return selection


def _listing_window(entries, signature, *, focus_molecule=""):
    options = ["50", "100", "250", "500", *(["All"] if len(entries) <= 1000 else [])]
    key = "data-browser-rows-shown"
    st.session_state.setdefault(key, "All" if len(entries) <= 1000 else "50")
    limited = st.session_state[key] not in options
    if limited:
        st.session_state[key] = "50"
    with st.popover(f"{len(entries):,} results · " + ("all shown" if st.session_state[key] == "All" else "paged"),
            icon=":material/table_rows:", width="stretch"):
        amount = st.selectbox("Rows shown", options, key=key,
            format_func=lambda value: f"All {len(entries):,} · no paging" if value == "All" else value,
            help="Show the whole matching view for up to 1,000 entries. Larger folders retain paging to limit transfer and memory.")
        size = max(1, len(entries)) if amount == "All" else int(amount)
        pages = max(1, math.ceil(len(entries) / size))
        view_signature = hashlib.sha256(json.dumps([signature, amount]).encode()).hexdigest()[:20]
        position = next((position for position, entry in enumerate(entries) if entry["result"]["molecule_id"] == focus_molecule), None) if focus_molecule else None
        if position is not None and pages > 1:
            st.session_state["data-browser-page-" + view_signature] = position // size + 1
        page = int(st.number_input("Page", min_value=1, max_value=pages, value=1, step=1,
            key="data-browser-page-" + view_signature)) if pages > 1 else 1
        visible = entries[(page - 1) * size:page * size]
        st.caption(f"All {len(entries):,} matching entries · no paging" if amount == "All" else
            f"{len(entries):,} matching entries · page {page} of {pages} · {len(visible):,} shown")
        if limited:
            st.caption("This folder has more than 1,000 entries; paging is enabled to keep the table responsive.")
    return visible, page, view_signature


def _table_controls(database, index, entries, visible, snapshot, signature, page, *, cache=None):
    table_key = ("binderscope-summary-v5", snapshot, signature, page)
    cached = st.session_state.get("data-browser-page-table")
    if not cached or cached[0] != table_key:
        cached = (table_key, result_table(database, index, visible, cache=cache))
        st.session_state["data-browser-page-table"] = cached
    rows = cached[1]
    available = table_columns(rows)
    automatic = table_profile(rows)
    st.caption(automatic + " columns · original scores and separately linked post-analysis")
    with st.expander("Table columns & CSV export"):
        preset = st.selectbox("Column preset", COLUMN_PRESETS, key="data-browser-column-preset")
        details = st.toggle("Show optional BinderScope details", key="data-browser-binderscope-details",
            disabled=preset in {"Compact", "All columns"},
            help="Add model, version, confidence, terminal accessibility and target-reference RMSD to activity presets. Compact stays compact; All columns includes everything.") if PREFIX + "Status" in available else False
        columns_key = hashlib.sha256(json.dumps([available, table_profile(rows, preset), details]).encode()).hexdigest()[:16]
        columns = st.multiselect("Visible columns", available, default=default_columns(rows, preset, binderscope_details=details), format_func=column_label,
            key="data-browser-columns-v4-" + columns_key)
        st.caption("Native scores and inputs belong to each row’s displayed evaluation, not a best-score merge. "
            "Analysis numbers require one completed, verified same-structure evaluation (or the exact MD segment). "
            "Other models, derived structures and repeated evaluations are labelled separately. "
            "Missing values stay blank; native scales can differ between tools.")
        if PREFIX + "Status" in available:
            st.caption("BinderScope is post-analysis, not an admission decision. Hotspots show contacted / requested. "
                "Tag accessibility is geometry-based; terminal confidence describes its reliability, not experimental validation.")
        export_mode = st.selectbox("CSV rows", ["Current listing", "All original evaluations"], key="data-browser-csv-mode")
        st.caption(f"Exports all {len(entries):,} matching entries across every page, with all columns, "
            "including recorded input parameters, native results and source identifiers. "
            "Linked analysis has separate columns and complete JSON evidence. The evaluation option expands "
            "original models / replicates without adding analysis-only rows.")
        export_key = ("binderscope-summary-v5", snapshot, signature, export_mode)
        prepared = st.session_state.get("data-browser-prepared-csv")
        if prepared and prepared[0] != export_key:
            st.session_state.pop("data-browser-prepared-csv", None)
            prepared = None
        if st.button("Prepare CSV", key="data-browser-prepare-csv"):
            with st.spinner("Preparing all matching SQL rows…"):
                selected = evaluation_entries(entries) if export_mode == "All original evaluations" else entries
                exported = result_table(database, index, selected, cache=cache)
                prepared = (export_key, table_csv(exported), len(exported))
                st.session_state["data-browser-prepared-csv"] = prepared
        if prepared:
            st.download_button(f"Download CSV ({prepared[2]:,} rows)", prepared[1], file_name="saved-results-" + signature + ".csv",
                mime="text/csv", on_click="ignore", key="data-browser-table-csv")
    return rows, columns


def _render_structure(database, user, record, preferred):
    structures = record["structures"]
    if not structures:
        if (record.get("selected_evaluation") or {}).get("representation_id"):
            st.warning("This evaluation has an unresolved saved-pose reference. Other structures may exist for this molecule, "
                "but a different run, replicate or tool's pose is not substituted for its scores.")
        else:
            st.info("No saved coordinates are attached to this result. This does not mean the molecule has no other structures; "
                "some calculations contain scores only.")
        return
    if record.get("structure_link_status"):
        st.caption("Pose recovered by exact checksum and matching source context; the original legacy reference is retained.")
    lookup = {entry["id"]: entry for entry in structures}
    selected_evaluation = record.get("selected_evaluation") or {}
    model_context = (record["id"], preferred, selected_evaluation.get("id"))
    remembered = st.session_state.get("data-browser-model-choice") or {}
    default = remembered.get("chosen") if remembered.get("context") == model_context else preferred
    chosen = st.selectbox("Saved structure", list(lookup),
        index=list(lookup).index(default) if default in lookup else 0,
        format_func=lambda value: lookup[value]["label"],
        key="data-browser-structure-" + record["id"] + (preferred or "") + str(selected_evaluation.get("id") or ""))
    st.session_state["data-browser-model-choice"] = {"context": model_context, "chosen": chosen}
    viewer_height = 720 if st.session_state.get("data-browser-layout") == "Structure" else 320
    structure = lookup[chosen]
    if selected_evaluation.get("representation_id") and chosen != selected_evaluation["representation_id"]:
        st.caption("Comparison model: the scores below remain attached to the selected evaluation’s own structure, not these comparison coordinates.")
    try:
        path = database.structure_file(user, chosen)
        context = structure.get("context") or {}
        from mn_molecular_design_studio.database_ui import _structure_text
        from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer import molstar_custom_component
        from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer.dataclasses import StructureVisualization

        if path.stat().st_size > 50_000_000:
            st.caption("This structure is too large for the inline viewer. Download it for inspection.")
        elif record["type"] == "small_molecule" and context.get("database_structure_policy") in {"best_docking_pose", "best_prediction"}:
            from mn_molecular_design_studio.docking_results import pose_visualizations
            visualizations = pose_visualizations(path.read_text(errors="replace"), path.name,
                include_pocket_sidechains=False)
            if not getattr(database, "input_only", False) and not any(item.representation_type == "cartoon" for item in visualizations):
                try:
                    receptor = saved_receptor_file(database, user, structure)
                    visualizations = pose_visualizations(path.read_text(errors="replace"), path.name,
                        receptor_text=receptor.read_text(errors="replace"), receptor_name=receptor.name,
                        include_pocket_sidechains=False)
                except (OSError, ValueError, PermissionError) as exc:
                    st.warning("Ligand coordinates only; the receptor frame could not be verified: " + str(exc))
            molstar_custom_component(visualizations, key="data-browser-molstar-native-" + chosen, height=viewer_height,
                show_controls=False, show_structure_tools=False)
        else:
            representation = "ball-and-stick" if record["type"] == "small_molecule" else "cartoon"
            molstar_custom_component([StructureVisualization(pdb=_structure_text(path, context), color="chain-id",
                representation_type=representation)], key="data-browser-molstar-native-" + chosen, height=viewer_height,
                show_controls=False, show_structure_tools=False)
        st.download_button("Download selected structure", path.read_bytes(), file_name=record["name"][:70] + path.suffix,
            key="data-browser-download-" + chosen)
    except (OSError, ValueError, PermissionError) as exc:
        st.warning("The saved coordinates could not be opened: " + str(exc))
        return
    if context.get("rfd4_pae") and st.toggle("Show internal refold PAE", key=f"saved-rfd4-pae-{chosen}"):
        from mn_molecular_design_studio.cofolding_results_ui import render_openfold3_pae

        try:
            model = database.rfd4_pae_file(user, chosen, "model")
            pae = database.rfd4_pae_file(user, chosen, "matrix")
            render_openfold3_pae(model, pae, key=f"saved-pae-{chosen}")
        except (OSError, ValueError, RuntimeError) as exc:
            st.info(f"The saved OpenFold3 confidence could not be opened. {exc}")


def _render_details(database, user, cache, index, entry):
    result = entry["result"]
    record = load_saved_record(database, user, result)
    if record is None:
        st.warning("This result is no longer accessible. Choose another folder or refresh.")
        return
    heading = concise_result_label(entry["outcome"])
    st.markdown("#### " + heading)
    from mn_molecular_design_studio.entity_properties_ui import render_record_properties
    render_record_properties(database, user, record, index=index, outcome=entry["outcome"],
        read_only=read_only_workspace(), key_prefix="data-browser-result")
    st.caption(" · ".join(value for value in [entry["outcome"].get("target_reference_label") or entry["outcome"].get("target_name"),
        ", ".join(tool_label(tool) for tool in entry["outcome"]["origin_tools"]),
        entry["outcome"].get("design_created_at") if entry["outcome"]["outcome_kind"] == "Design" else run_label(result["run_id"])] if value))
    st.caption(entry["outcome"]["activity"] + " · " + entry["outcome"]["eligibility"])
    from mn_molecular_design_studio.rfd4_representations import expects_refolded_model
    if result.get("receipt", {}).get("representation_role") == "refolded_validation":
        st.caption("Refolded validation model · generated coordinates remain attached for comparison. Scores belong to their recorded model.")
    elif expects_refolded_model(record):
        st.warning("Only the generated model is saved for this result. Refolding scores exist, but the refolded coordinates have not been recovered; this model is not a substitute.")
    _render_structure(database, user, record, entry["structure_id"])
    cohort = design_subject(index, result)
    if cohort:
        cohort_results = cache.get_results(database, index, cohort, include_metrics=True)
        connected, unresolved = connected_evaluations(index, result, cohort_results)
    else:
        cohort_results, connected, unresolved = [], [], 0
    original_results = producing_evaluations(index, entry["outcome"], connected)
    originals = {row["id"] for row in original_results}
    analyses = analysis_results(connected)
    validations = [row for row in connected if row not in analyses and (
        row.get("receipt", {}).get("role") == "design_validation" or row.get("receipt", {}).get("design_refold"))]
    other_calculations = [row for row in connected if row["id"] not in originals and row not in analyses and row not in validations]
    scores_tab, parameters_tab, history_tab = st.tabs(["Scores & analyses", "Parameters", "History"])
    from mn_molecular_design_studio.database_ui import _render_evaluation_measurements
    with scores_tab:
        st.markdown("**Original result scores**")
        if not original_results:
            preserved = admission_measurements(record)
            if preserved:
                st.caption("Generator measurements preserved in SQL admission evidence; these are not a newly reconstructed Evaluation.")
                st.dataframe(preserved, hide_index=True, use_container_width=True, column_order=["Value", "Measurement", "Tool"],
                    column_config={"Value": st.column_config.TextColumn(width="small")})
            elif record["observations"]:
                for observation in record["observations"]:
                    st.dataframe([{"Metric": name, "Value": json.dumps(value, ensure_ascii=False)}
                        for name, value in observation["metrics"].items()], hide_index=True, use_container_width=True)
            else:
                st.caption("No original measurements are recorded on this SQL result. Missing values are not zero.")
        for original in original_results:
            label = calculation_label(original)
            with st.expander(label, expanded=len(original_results) == 1 or original["evaluation_id"] == result["evaluation_id"]):
                st.caption(original["structure_relation"])
                _render_evaluation_measurements(original, compact=True)
        st.markdown("**Attached post-analysis**")
        if not analyses:
            st.caption("No post-analysis is recorded for this design / target context.")
        for analysis in analyses:
            with st.expander(tool_label(analysis["tool"]) + " · " + analysis["structure_relation"],
                             expanded=analysis["tool"] == "binderscope" and len(analyses) == 1):
                if analysis["tool"] == "binderscope":
                    from mn_molecular_design_studio.binderscope_ui import render_design_result
                    from mn_molecular_design_studio.binderscope import saved_result
                    assessment = saved_result(database, user, analysis)
                    if assessment:
                        assessment = {**assessment, "_studio_target_label": entry["outcome"].get("target_reference_label") or entry["outcome"].get("target_name") or "Saved target"}
                        render_design_result(assessment, database=database, user=user, key=analysis["observation_id"], compact=True)
                    else:
                        st.caption("This assessment is no longer accessible.")
                else:
                    _render_evaluation_measurements(analysis, compact=True)
        if validations:
            with st.expander(f"Design-internal refolding ({len(validations)})"):
                for validation in validations:
                    st.markdown("**" + tool_label(validation["tool"]) + "**")
                    _render_evaluation_measurements(validation, compact=True)
        if unresolved:
            st.caption(f"{unresolved} calculations have unresolved links and are not presented as attached scores.")
    with parameters_tab:
        validation = record.get("metadata", {}).get("validation_acceptance")
        if validation:
            with st.expander("Saved admission policy and measurement provenance"):
                st.json(validation, expanded=False)
        if record["observations"]:
            for observation in record["observations"]:
                st.caption(observation["namespace"] + " · parameters and context recorded in SQL")
                st.json(observation["context"], expanded=False)
        elif not validation:
            st.caption("No calculation parameters are recorded for this saved molecular identity.")
    with history_tab:
        st.caption("Display labels do not rename scientific identities. Original names and identifiers remain here.")
        st.write({"Original name": entry["outcome"]["name"], "Molecule ID": record["id"],
                  "Admission": entry["outcome"].get("acceptance_level", "Not recorded")})
        naming = record.get("metadata", {}).get("studio_naming") or {}
        if naming.get("native_names"):
            st.write({"Native names": naming["native_names"]})
        creation = record.get("metadata", {}).get("design_creation")
        if creation:
            with st.expander("Design creation time and original output evidence"):
                st.json(creation, expanded=False)
        st.dataframe([{"Tool": tool_label(row["tool"]), "Run": row["run_id"], "Candidate": row.get("candidate", ""),
                       "Model / pose": row["model"], "Collection": row["collection"]} for row in entry["outcome"]["originals"]],
                     hide_index=True, use_container_width=True)
        with st.expander("Source evidence"):
            st.json(record["sources"], expanded=False)
        if other_calculations:
            with st.expander(f"Other calculations for this molecule ({len(other_calculations)})"):
                st.caption("These are not asserted to be post-analysis of the selected result. Open their own saved result to inspect scores.")
                st.dataframe([{"Method": tool_label(row["tool"]), "Run": row["run_id"],
                    "Relationship": row["structure_relation"]} for row in other_calculations],
                    hide_index=True, use_container_width=True)


def render_data_browser(user, projects, *, database=None, workflow_schema=None, submit=None):
    if database is None:
        configured = os.environ.get("MDS_DATA_BROWSER_DATABASE_ROOT", "").strip()
        database = ReadOnlyLibrary(Path(configured) if configured else None)
    st.markdown('<style>[data-testid="stMainBlockContainer"]{padding:1.5rem 1.6rem;max-width:none}.stMarkdown h4{overflow-wrap:anywhere}</style>', unsafe_allow_html=True)
    destination = database_destination(st.session_state)
    dataset_request = st.session_state.pop("dataset-results-request", None)
    if destination and dataset_request and destination.get("group_id") == dataset_request["group_id"]:
        st.session_state["data-browser-dataset-context"] = dataset_request
    st.title("Data browser")
    st.caption("Saved scientific outcomes, structures and attached analyses. Reusable input libraries are managed in Datasets.")
    from mn_molecular_design_studio.entity_properties_ui import render_notice
    render_notice()
    render_library_actions(database, user, projects, browser_only=True)
    if not database.path.is_file():
        st.info("No saved SQL data is available yet. Import reusable inputs in Datasets, or admit passing results from a run.")
        return
    cache = _cache(database)
    with st.spinner("Opening saved data…"):
        index = cache.get_index(database, user, scope="all",
            collection="all" if st.session_state.get("data-browser-indexed-history") else "library",
            hide_smoke_test=False, include_related=True, hide_pipeline_batches=True)
    if st.session_state.get("data-browser-indexed-history"):
        st.info("Historical indexed data is included. Stored does not mean reviewed or passed. Library actions can restore the curated-only view.")
    index["project_names"] = {project["id"]: project["name"] for project in projects}
    if destination is not None:
        st.session_state["data-browser-source-group"] = destination.get("group_id")
    source_group = st.session_state.get("data-browser-source-group")
    dataset_context = st.session_state.get("data-browser-dataset-context")
    scopes_key = ("scientific-tree-v3-cached-navigation", cache.index_key, id(index), source_group,
        tuple(sorted((dataset_context or {}).items())),
        tuple((row["id"], row["name"]) for row in projects))
    if st.session_state.get("data-browser-scopes-key") != scopes_key:
        st.session_state["data-browser-scopes"] = scientific_scopes(index, projects)
        st.session_state["data-browser-scopes-key"] = scopes_key
        st.session_state.pop("data-browser-inventory", None)
        st.session_state.pop("data-browser-inventory-key", None)
    scopes = st.session_state["data-browser-scopes"]
    if source_scope := source_collection_scope(index, source_group):
        scopes[source_scope["id"]] = scopes.get(source_scope["id"], source_scope)
    connected_scope = dataset_results_scope(scopes, index, dataset_context["group_id"],
        dataset_context.get("molecule_id", "")) if dataset_context else None
    if dataset_context and dataset_context.get("molecule_id") and connected_scope is None:
        st.warning("No unique connected-data route was found for this dataset member. Its input data remains available in Datasets; result relationships have not been guessed.")
    if destination is not None:
        desired = "collection:" + str(destination["group_id"]) if destination.get("group_id") else (
            "project:" + (destination.get("project_id") or "personal") if "project_id" in destination else "projects")
        desired = connected_scope or desired
        if desired not in scopes:
            st.warning("The requested collection or project has no accessible saved data in this view. No other collection has been substituted.")
            desired = "collections" if destination.get("group_id") else "projects"
        st.session_state.update({"data-browser-entry-v2": scope_path(scopes, desired)[0],
            "data-browser-scope": desired, "data-browser-folder": desired,
            "data-browser-focus-molecule": "" if connected_scope else destination.get("molecule_id") or ""})
    for ancestor in st.session_state.get("data-browser-scope-path", []):
        if ancestor in scopes:
            expand_scope(scopes, index, ancestor)
    if not index["molecules"] and not index["groups"]:
        st.info("No saved molecules are accessible to this account.")
        return
    layout, left, middle, right = _workspace()
    roots = {"projects": "Projects", "targets": "Targets", "ligands": "Compounds", "collections": "Libraries & datasets"}
    with left:
        st.markdown("### Explore")
        root = st.selectbox("Start at", list(roots), format_func=roots.get, key="data-browser-entry-v2")
    default = next((key for key, scope in scopes.items() if scope["parent"] == root), root)
    scope_id = st.session_state.get("data-browser-scope", default)
    if scope_id not in scopes or scope_path(scopes, scope_id)[0] != root:
        retained = [key for key in st.session_state.get("data-browser-scope-path", []) if key in scopes and scope_path(scopes, key)[0] == root]
        scope_id = retained[-1] if retained else default
        st.session_state["data-browser-scope"] = scope_id
        st.session_state["data-browser-folder"] = scope_id
    st.session_state["data-browser-scope"] = scope_id
    scope = scopes[scope_id]
    expand_scope(scopes, index, scope_id)
    st.session_state["data-browser-scope-path"] = scope_path(scopes, scope_id)
    collection_view = any(scopes[ancestor].get("route") == "collection" for ancestor in scope_path(scopes, scope_id))
    with st.spinner("Reading this folder’s saved results…"):
        inventory = cache.get_saved_results(database, index, scope["subject"], run_view=scope["run_view"],
            include_references=collection_view, fallback_references=(scope.get("entity") or {}).get("kind") == "target" or
                scope.get("route") == "scientific-project")
        inventory = [outcome for outcome in inventory if outcome["outcome_kind"] not in {"Imported molecule", "Reference molecule"}]
    if scope.get("prepared_snapshot") != scopes_key:
        scope_counts(scopes, index, scope, inventory)
        children, selections = scientific_branches(index, scope, inventory)
        scope["result_branches"] = children
        scope["result_selections"] = selections
        scope["prepared_snapshot"] = scopes_key
    children, selections = scope["result_branches"], scope["result_selections"]
    selected_id = st.session_state.get("data-browser-folder", scope_id)
    if selected_id not in selections and selected_id != scope_id:
        selected_id = scope_id
    search_context = (str(database.path.resolve()), str(user["id"]), scope_id, selected_id)
    if st.session_state.get("data-browser-result-query-context") != search_context:
        _clear_result_search()
        st.session_state["data-browser-result-query-context"] = search_context
        st.session_state.pop("binderscope-browser-preflight", None)
    with left:
        query = st.text_input("Find a folder", placeholder="Target, ligand, project, sequence…", key="data-browser-folder-query")
        with st.expander("Refine this level"):
            molecule_type = st.selectbox("Molecule type", ["", "protein", "small_molecule", "dna", "rna"],
                format_func=lambda value: {"": "All molecule types", "small_molecule": "Small molecules", "protein": "Proteins", "dna": "DNA", "rna": "RNA"}[value],
                key="data-browser-type")
            activities = sorted({row["activity"] for row in inventory})
            activity = st.selectbox("Activity", ["", *activities], format_func=lambda value: value or "All activities", key="data-browser-activity-" + scope_id)
            construct_inventory = [row for row in inventory if not activity or row["activity"] == activity]
            constructs = target_constructs(index, [entry["outcome"] for entry in listed_results(
                construct_inventory, selection=selections.get(selected_id, {}))])
            construct = st.selectbox("Target construct", ["", *constructs],
                format_func=lambda value: constructs.get(value, "All constructs"), key="data-browser-construct-" + selected_id) if len(constructs) > 1 else ""
            if len(constructs) > 1:
                st.caption("Constructs refine this activity tree; they are not duplicate folders. "
                    "Results without an exact recorded construct remain visible with All constructs.")
            coordinates_only = st.checkbox("With saved coordinates", key="data-browser-coordinates")
            st.caption("These controls refine the selected folder, not the stored records.")
        limit = st.session_state.get("data-browser-folder-limit", 100)
        nodes, hidden = navigation_tree(scopes, active_scope=scope_id, children=children, query=query, limit=limit, root=root, molecule_type=molecule_type)
        expanded = [*scope_path(scopes, scope_id), *st.session_state.get("data-browser-expanded", [])]
        if selected_id in selections:
            expanded += [node["id"] for node in children if selected_id.startswith(node["id"])]
        event = render_tree(nodes, selected_id=selected_id, expanded_ids=expanded, key="data-browser-folders-v2")
        if event and event.get("token") != st.session_state.get("data-browser-tree-event"):
            st.session_state["data-browser-tree-event"] = event.get("token")
            chosen = event.get("selectedId")
            branch_scope = next((identifier for identifier, entry in scopes.items()
                                 if chosen in entry.get("result_selections", {})), None)
            if chosen in scopes or branch_scope:
                destination = chosen if chosen in scopes else branch_scope
                st.session_state["data-browser-folder"] = chosen
                st.session_state["data-browser-scope"] = destination
                st.session_state["data-browser-scope-path"] = scope_path(scopes, destination)
                st.session_state["data-browser-expanded"] = event.get("expandedIds", [])
                st.rerun()
        if hidden and st.button(f"Show more folders ({hidden:,} remaining)", key="data-browser-more-folders"):
            st.session_state["data-browser-folder-limit"] = limit + 100
            st.rerun()
        if st.button("Refresh SQL view", icon=":material/refresh:", key="data-browser-refresh"):
            cache.close()
            st.session_state.pop("data-browser-cache", None)
            st.rerun()
    md_context = None
    with middle:
        selection = selections.get(selected_id, {})
        breadcrumb = [scopes[key]["label"] for key in scope_path(scopes, scope_id)]
        breadcrumb += selection.get("path", [])
        st.caption(" / ".join(breadcrumb))
        if scope.get("reference_family"):
            st.caption("Includes uniquely matching target crops. Exact construct identities and measurements remain separate for every result.")
        from mn_molecular_design_studio.entity_properties_ui import render_folder_properties
        if selected_id == scope_id:
            render_folder_properties(database, user, index, scope, projects,
                read_only=read_only_workspace(), key_prefix="data-browser-folder")
            if (scope.get("entity") or {}).get("kind") == "dataset":
                if st.button("Open dataset input members", icon=":material/library_books:", key="data-browser-dataset-members"):
                    from mn_molecular_design_studio.navigation import open_datasets
                    open_datasets(st.session_state, group_id=scope["entity"]["group_id"])
                    st.rerun()
        result_controls = st.container() if layout == "Split" else st.popover("Results, filters & analysis", icon=":material/table_view:")
        md_selected = bool(selection.get("md_outcome_id"))
        with result_controls:
            if md_selected:
                listing, query = "Saved results", ""
            else:
                if layout != "Structure":
                    st.markdown("### Saved results")
                listing = st.segmented_control("Listing", ["Saved results", "All saved structures"], default="Saved results", key="data-browser-listing")
                search_column, clear_column = st.columns([5, 1], vertical_alignment="bottom")
                with search_column:
                    query = st.text_input("Search this folder", placeholder="Name, original identifier, tool or campaign…", key="data-browser-result-query")
                clear_column.button("Clear", key="data-browser-clear-result-query", on_click=_clear_result_search,
                    disabled=not bool(query), help="Clear this folder’s search. Choosing another folder also clears it.")
        refined = [row for row in inventory if (not activity or row["activity"] == activity) and (
            not construct or row.get("target_id") == construct) and (
            not coordinates_only or row["structure_count"]) and (not molecule_type or any(
                component["type"] == molecule_type for component in index["molecules"][row["molecule_id"]]["components"]))]
        entries = listed_results(refined, selection=selection, query=query, structures=listing == "All saved structures")
        entries.sort(key=lambda entry: not bool(entry["result"]["structure_ids"]))
        signature = hashlib.sha256(json.dumps([scope_id, selected_id, query, listing, [entry["id"] for entry in entries]]).encode()).hexdigest()[:20]
        if md_selected:
            visible, page, view_signature = entries[:50], 1, signature
            analysis_toolbar = None
        else:
            window_control, analysis_toolbar = st.columns([1, 2], vertical_alignment="bottom")
            with window_control:
                visible, page, view_signature = _listing_window(entries, signature,
                    focus_molecule=st.session_state.get("data-browser-focus-molecule", ""))
        with result_controls:
            if not md_selected:
                if any(entry["outcome"]["outcome_kind"] == "Design" for entry in entries):
                    st.caption("Passed means the recorded native or final acceptance filters passed, not that every later analysis passed. "
                        "Tree counts are whole-folder totals; search and refinements can reduce the list.")
        page_selection = _page_selection(view_signature, page, visible)
        if focus_molecule := st.session_state.pop("data-browser-focus-molecule", ""):
            match = next((entry for entry in visible if entry["result"]["molecule_id"] == focus_molecule), None)
            if match:
                page_selection["focused"] = match["id"]
            else:
                st.info("The requested molecular record has no saved scientific outcome at this level. It may be a reusable input or unreviewed run output; no passing result has been inferred.")
        selected_entries = []
        if visible and md_selected:
            focused = visible[0]
        elif visible:
            with result_controls:
                rows, columns = _table_controls(database, index, entries, visible, scopes_key, view_signature, page, cache=cache)
            if layout == "Structure":
                lookup = {entry["id"]: entry for entry in visible}
                chosen = st.selectbox("Result on this page", list(lookup), index=list(lookup).index(page_selection["focused"]),
                    format_func=lambda value: concise_result_label(lookup[value]["outcome"]),
                    key="data-browser-focused-result-" + view_signature + str(page))
                page_selection["focused"] = chosen
            else:
                from mn_molecular_design_studio.binderscope_ui import request_saved_row

                context = view_signature + ":" + str(page)
                click_key = "data-browser-binderscope-row-" + context
                requests = [{"entry_id": entry["id"], "label": inspector_label(row),
                    "assessments": [(choice["evaluation_id"], choice["observation_id"]) for choice in inspector_choices(row)]
                    if entry["outcome"]["outcome_kind"] == "Design" else []} for entry, row in zip(visible, rows)]
                has_viewers = any(request["assessments"] for request in requests)
                displayed_rows = [{**row, VIEWER_COLUMN: request["label"] if request["assessments"] else None}
                    for row, request in zip(rows, requests)] if has_viewers else rows
                displayed_columns = list(columns)
                if has_viewers:
                    displayed_columns.insert(displayed_columns.index("Name") + 1 if "Name" in displayed_columns else 0, VIEWER_COLUMN)
                displayed_rows = [{field: row.get(field) for field in displayed_columns} for row in displayed_rows]
                retained_rows = [position for position, entry in enumerate(visible) if entry["id"] in page_selection["selected"]]
                event = st.dataframe(displayed_rows, hide_index=True, width="stretch", height=720 if layout == "Table" else 565,
                    column_order=displayed_columns,
                    column_config={**{field: column_label(field) for field in table_columns(rows)},
                        **{field: st.column_config.NumberColumn(column_label(field), format=number_format)
                            for field, number_format in NUMBER_FORMATS.items() if field in table_columns(rows)},
                        "Name": st.column_config.TextColumn(width="small", pinned=has_viewers),
                        **({VIEWER_COLUMN: st.column_config.ButtonColumn("BinderScope", type="tertiary", width=155, pinned=True,
                            help="Open this row’s saved analysis, without selecting it or starting another calculation.",
                            key=click_key, on_click=request_saved_row, args=(click_key, context, requests))} if has_viewers else {}),
                        "Method": st.column_config.TextColumn(width="small"),
                        "Project": st.column_config.TextColumn(width="medium"),
                        "Target": st.column_config.TextColumn(width="medium"),
                        "Run": st.column_config.TextColumn(width="medium"),
                        "Structures": st.column_config.NumberColumn(width="small"),
                        "Admission": st.column_config.TextColumn(width="small")},
                    on_select="rerun", selection_mode="multi-row", selection_default={"selection": {"rows": retained_rows}},
                    key="data-browser-table-" + view_signature + str(page))
                selected_ids = [visible[position]["id"] for position in event.selection.rows if 0 <= position < len(visible)]
                if selected_ids != page_selection["selected"]:
                    page_selection["focused"] = selected_ids[-1] if selected_ids else visible[0]["id"]
                page_selection["selected"] = selected_ids
            selected_entries = [entry for entry in visible if entry["id"] in page_selection["selected"]]
            focused = next(entry for entry in visible if entry["id"] == page_selection["focused"])
            from mn_molecular_design_studio.binderscope_ui import render_saved_browser_action, render_saved_row_request
            if layout == "Structure":
                position = next(position for position, entry in enumerate(visible) if entry["id"] == focused["id"])
                render_saved_browser_action(database, user, index, focused, rows[position], signature=signature)
            else:
                render_saved_row_request(database, user, index, visible, context=context)
        else:
            focused = None
            if inventory:
                st.info(f'No saved results match the search “{query}” in this folder. The search does not remove SQL data.' if query else
                    "No saved results match the current refinements. This does not mean the SQL data was removed.")
                if query or activity or construct or coordinates_only or molecule_type:
                    st.button("Show this folder’s saved data", key="data-browser-reset-result-refinements",
                        on_click=_reset_folder_refinements, args=(scope_id, selected_id),
                        help="Clear the search and result refinements without changing any saved records.")
            else:
                st.info("This folder has no saved scientific outcomes at this level. Reusable inputs remain in Datasets. Designs appear only after a recorded filter pass or explicit later selection; run history is not imported here.")
        from mn_molecular_design_studio.binderscope_ui import render_browser_action
        if analysis_toolbar is not None:
            with analysis_toolbar:
                render_browser_action(database, user, index, entries, selected_entries, signature=signature,
                                      workflow_schema=workflow_schema, submit=submit)
        if not md_selected:
            project_id = (scope.get("entity") or {}).get("project_id", "")
            render_selection_actions(database, user, projects, selected_entries, project_id=project_id)
        if selection.get("experiment"):
            jobs = source_jobs(index, inventory, selection)
            with result_controls.expander(f"Source jobs & provenance ({len(jobs):,})"):
                st.caption("Execution batches are retained here; they are not separate scientific experiments.")
                st.dataframe(jobs, hide_index=True, use_container_width=True)
        if focused and focused["outcome"]["outcome_kind"] == "Dynamics":
            from mn_molecular_design_studio.md_browser_ui import render_md_study
            context_panel = st.popover("MD analysis & snapshot selection", icon=":material/analytics:") if layout == "Structure" else st.container()
            with context_panel:
                md_context = render_md_study(database, user, index, focused, selection)
    if right is not None:
        with right:
            if focused and md_context is not None:
                from mn_molecular_design_studio.md_browser_ui import render_md_structure
                render_md_structure(database, user, focused, md_context)
            elif focused:
                _render_details(database, user, cache, index, focused)
            else:
                st.caption("Select a saved result to open its structure, parameters and attached analyses.")
