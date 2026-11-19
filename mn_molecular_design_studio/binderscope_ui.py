"""Explicit BinderScope launch and inspection surfaces; no calculations on render."""
from __future__ import annotations

import hashlib
import json
import os
from uuid import uuid4
from pathlib import Path

import streamlit as st

from mn_molecular_design_studio.binderscope import (
    METHOD_VERSION, REVISION, TOPOLOGIES, existing_assessments, prepare_browser_selection,
    project_selections, settings_key, validate_assessment_settings,
)


BINDER_FORMAT_LABELS = {"unspecified": "Not recorded", "linear protein": "Protein",
    "linear peptide": "Linear peptide", "cyclic peptide": "Cyclic peptide (head-to-tail)"}


def _read_only_preview():
    return (os.environ.get("MDS_READ_ONLY_PREVIEW") == "1"
            or os.environ.get("MDS_API_BASE_URL", "").rstrip("/") == "http://127.0.0.1:1")


def settings_rows(sources, counts, included=None):
    return [{"Source": str(position), "Include": counts.get((source["molecule_id"], source["structure_id"], source["group_id"]), 0) == 0
            if included is None else str(position) in included, "Design": source["label"],
        "Structure": source.get("structure_label") or source["structure_id"],
        "Topology": source["assessment"]["component_format"],
        "Confidence": "Confirmed pLDDT (0–100)" if source["assessment"]["b_factor_interpretation"] == "plddt_0_100" else "Unknown / uninterpreted",
        "Existing analyses": counts.get((source["molecule_id"], source["structure_id"], source["group_id"]), 0),
        "Hotspots": (str(len(source.get("target_context", {}).get("hotspots", {}).get("requested") or [])) + " requested")
            if source.get("target_context", {}).get("hotspots", {}).get("status") == "recorded" else
            source.get("target_context", {}).get("hotspots", {}).get("status", "not_recorded").replace("_", " "),
        "Target reference": source.get("target_context", {}).get("reference", {}).get("status", "unavailable"),
        "Evidence": source["assessment"].get("topology_basis", "") + " " + source["assessment"].get("confidence_basis", "")}
        for position, source in enumerate(sources)]


def edited_selection(selection, rows, *, included_only=True):
    if len(rows) != len(selection["sources"]) or {row["Source"] for row in rows} != {str(position) for position in range(len(rows))}:
        raise ValueError("The selected structure inventory changed. Reopen BinderScope.")
    rows = {row["Source"]: row for row in rows}
    sources = []
    for position, source in enumerate(selection["sources"]):
        row = rows[str(position)]
        if type(row.get("Include")) is not bool:
            raise ValueError("Choose which designs to include using the checkboxes.")
        if row["Confidence"] not in {"Confirmed pLDDT (0–100)", "Unknown / uninterpreted"}:
            raise ValueError("Choose a supported confidence encoding.")
        assessment = dict(source["assessment"])
        values = {"component_format": row["Topology"], "b_factor_interpretation":
            "plddt_0_100" if row["Confidence"] == "Confirmed pLDDT (0–100)" else "uninterpreted"}
        for field, value in values.items():
            if assessment[field] != value:
                assessment[field] = value
                assessment["topology_basis" if field == "component_format" else "confidence_basis"] = "User declaration for this exact selected structure."
        assessment = validate_assessment_settings(assessment)
        if row["Include"] or not included_only:
            sources.append({**source, "assessment": assessment})
    return {**selection, "sources": sources}


def _set_editor_selection(prepared, selection, included):
    prepared["selection"] = selection
    prepared["included"] = included
    prepared["editor_version"] += 1


def _dismiss_analysis():
    st.session_state.pop("binderscope-browser-preflight", None)


@st.dialog("Analyze designs · BinderScope", width="large", on_dismiss=_dismiss_analysis)
def _browser_dialog(database, user, index, chosen, token, *, workflow_schema, submit):
    if st.session_state.get("studio-session-token"):
        from mn_molecular_design_studio.auth import resolve_session
        current = resolve_session(st.session_state["studio-session-token"])
        if not current or current["id"] != user["id"] or current.get("role") != user.get("role"):
            st.error("Your sign-in or permissions changed. Reopen the browser before submitting analysis.")
            return
    key = "binderscope-browser-preflight"
    prepared = st.session_state.get(key) or {}
    if prepared.get("token") != token:
        with st.spinner("Checking exact saved structures and recorded settings…"):
            selection, excluded = prepare_browser_selection(database, user, chosen)
            counts = existing_assessments(database, user, selection["sources"])
        prepared = {"token": token, "selection": selection, "excluded": excluded, "counts": counts,
                    "editor_version": 0, "submitted": {}}
        st.session_state[key] = prepared
    selection = prepared["selection"]
    sources = selection["sources"]
    st.write(f"**{len(sources):,} eligible designs** · {len(prepared['excluded']):,} excluded")
    st.caption("Only admitted SQL designs are assessed. Original scores and admission decisions stay unchanged. "
               "Each result attaches to its exact design and saved structure; projects keep separate access boundaries.")
    if prepared["excluded"]:
        with st.expander("Excluded designs and missing inputs"):
            st.dataframe(prepared["excluded"], hide_index=True, use_container_width=True)
    if not sources:
        return
    st.caption("Checked = analyze. Designs with a completed assessment on this exact complex are unchecked by default. "
               "Tick them to assess again; use View BinderScope in the main table to inspect existing results instead.")
    rows = settings_rows(sources, prepared["counts"], prepared.get("included"))
    edited = st.data_editor(rows, hide_index=True, use_container_width=True, height=min(350, 80 + 35 * len(rows)),
        column_order=["Include", "Design", "Hotspots", "Target reference", "Topology", "Confidence", "Existing analyses", "Structure", "Evidence"],
        column_config={"Source": None,
            "Include": st.column_config.CheckboxColumn("Analyze", required=True, width="small",
                help="Only checked designs will be queued. You may check an already assessed design to assess it again."),
            "Existing analyses": st.column_config.NumberColumn("Completed assessments", format="%d",
                help="Verified completed reports for this saved complex and chain selection, within the same project. This is not a pass/fail count."),
            "Topology": st.column_config.SelectboxColumn("Binder format", options=list(TOPOLOGIES), required=True,
                format_func=BINDER_FORMAT_LABELS.get,
                help="Protein/peptide describes binder format, not folded versus disordered structure. Linear peptide has unjoined backbone ends. Only head-to-tail cyclization removes both termini; disulfide/side-chain loops do not."),
            "Confidence": st.column_config.SelectboxColumn("Binder confidence", options=["Unknown / uninterpreted", "Confirmed pLDDT (0–100)"], required=True),
            "Design": st.column_config.TextColumn(width="medium"), "Evidence": st.column_config.TextColumn(width="large")},
        disabled=["Source", "Design", "Structure", "Evidence", "Existing analyses", "Hotspots", "Target reference"],
        key="binderscope-settings-" + token + str(prepared["editor_version"]))
    selection = edited_selection(selection, edited)
    all_settings = edited_selection(prepared["selection"], edited, included_only=False)
    included = [row["Source"] for row in edited if row["Include"]]
    select_unassessed, select_all, clear_selection = st.columns(3)
    select_unassessed.button("Unassessed only", key="binderscope-select-unassessed-" + token, on_click=_set_editor_selection,
        args=(prepared, all_settings, [row["Source"] for row in rows if not row["Existing analyses"]]))
    select_all.button("Include all", key="binderscope-select-all-" + token, on_click=_set_editor_selection,
        args=(prepared, all_settings, [row["Source"] for row in edited]))
    clear_selection.button("Clear selection", key="binderscope-clear-selection-" + token, on_click=_set_editor_selection,
        args=(prepared, all_settings, []))
    with st.expander("Adjust unknown settings / execution"):
        topology = st.selectbox("For designs without recorded binder format", ["Keep per-design settings", *TOPOLOGIES[1:]],
            format_func=lambda value: BINDER_FORMAT_LABELS.get(value, value), key="binderscope-unknown-topology-" + token)
        st.caption("Protein includes folded or disordered protein binders. Linear versus cyclic describes peptide backbone connectivity, "
            "not folding. A disulfide loop is not head-to-tail cyclization. Missing format remains Not recorded; length alone is not used to choose it.")
        confidence = st.selectbox("For designs without confirmed confidence encoding", ["Keep per-design settings", "I confirm B-factors are pLDDT (0–100)"],
                                  key="binderscope-unknown-confidence-" + token)
        st.caption("Do not confirm pLDDT from an average score alone. This declaration concerns binder Cα B-factors in these exact structures; target B-factors stay uninterpreted.")
        if st.button("Apply to unknown settings", key="binderscope-apply-" + token):
            for position, source in enumerate(all_settings["sources"]):
                if str(position) not in included:
                    continue
                assessment = source["assessment"]
                if topology != "Keep per-design settings" and assessment["component_format"] == "unspecified":
                    assessment.update(component_format=topology, topology_basis="User declaration for selected structures with unknown topology.")
                if confidence != "Keep per-design settings" and assessment["b_factor_interpretation"] == "uninterpreted":
                    assessment.update(b_factor_interpretation="plddt_0_100", confidence_basis="User declaration for selected structures with unknown encoding.")
            _set_editor_selection(prepared, all_settings, included)
            st.rerun(scope="fragment")
        cpu = st.number_input("CPU cores", min_value=1, max_value=64, value=4, step=1, key="binderscope-browser-cpu-" + token)
        repeat = st.checkbox("Create a new assessment instead of reusing this submission", key="binderscope-browser-repeat-" + token)
    known_confidence = sum(source["assessment"]["b_factor_interpretation"] == "plddt_0_100" for source in selection["sources"])
    known_topology = sum(source["assessment"]["component_format"] != "unspecified" for source in selection["sources"])
    selected_count = len(selection["sources"])
    groups = project_selections(selection, index)
    st.write(f"**{selected_count:,} selected** · {len(groups):,} project job(s)")
    st.caption(f"{known_topology:,}/{selected_count:,} binder formats recorded/declared · {known_confidence:,}/{selected_count:,} confidence encodings confirmed. "
               "Unknown confidence stays uninterpreted: geometry, contacts, patches, liabilities and structural comparisons still run; "
               "confidence-dependent tag/fusion checks are unavailable. Cyclic binders have no free termini.")
    st.caption("Recorded hotspots use each design’s own input target and residue mapping. Missing or ambiguous evidence is not a pass. "
               "For disordered targets, contacts describe this predicted conformation—not preservation of a fixed fold.")
    if _read_only_preview():
        st.info("Read-only development preview: inspect the real selection and settings here; submission is disabled.")
    parameters = {"use_recorded_settings": True, "cpu_cores": cpu}
    groups = project_selections(selection, index)
    identity = settings_key(selection, parameters)
    if st.button(f"Queue BinderScope · {selected_count:,} designs", type="primary", key="binderscope-browser-submit-" + token,
                 disabled=not selected_count or _read_only_preview() or submit is None or workflow_schema is None):
        if repeat:
            prepared.setdefault("repeat_identity", uuid4().hex)
            identity += prepared["repeat_identity"]
        for project_id, project_selection in groups.items():
            request_id = settings_key(project_selection, {**parameters, "submission": identity})
            if request_id in prepared["submitted"]:
                st.success("This project assessment is already queued.")
                continue
            try:
                result = submit(workflow_schema("analysis.binderscope"), "analysis.binderscope",
                    {**parameters, "binderscope_collection": project_selection}, project_id,
                    {"_execution": {"idempotency_key": request_id}, "prioritize": False}, "BinderScope")
                if result and result.get("job"):
                    prepared["submitted"][request_id] = result["job"]
            except (ValueError, PermissionError, OSError, RuntimeError) as exc:
                st.error("Could not queue this project assessment: " + str(exc))
                break
    if prepared["submitted"]:
        st.caption("Queued jobs run in the background. Their source-bound results appear under each design’s attached analyses after completion.")


def request_saved_row(click_key, context, requests):
    click = st.session_state.get(click_key)
    position = click.get("row") if isinstance(click, dict) else None
    if type(position) is int and 0 <= position < len(requests):
        request = requests[position]
        if request["assessments"] and click.get("label") == request["label"]:
            st.session_state["binderscope-row-request"] = {**request, "context": context}


@st.dialog("BinderScope · design inspector", width="large")
def _saved_row_dialog(database, user, index, entry, assessments, key):
    from mn_molecular_design_studio.binderscope import saved_result
    from mn_molecular_design_studio.binderscope_inspector import render_inspector
    from mn_molecular_design_studio.binderscope_table import inspector_choices
    from mn_molecular_design_studio.data_browser_analyses import linked_analysis_columns

    try:
        fresh = linked_analysis_columns(database, index, [entry])[0]
        choices = [choice for choice in inspector_choices(fresh)
            if (choice["evaluation_id"], choice["observation_id"]) in assessments]
        if not choices:
            st.warning("This saved assessment is no longer accessible for the selected model. Refresh the SQL view.")
            return
        chosen = choices[0] if len(assessments) == 1 else None
        if len(assessments) > 1:
            labels = {choice["observation_id"]: choice for choice in choices}
            identifier = st.selectbox("Saved BinderScope assessment", list(labels), index=None,
                placeholder="Choose an assessment; none is selected automatically",
                format_func=lambda value: " · ".join(str(part) for part in (
                    labels[value].get("analysis_version") or "Version not recorded",
                    labels[value].get("run_id") or "Run not recorded", value)), key=key + "-choice")
            chosen = labels.get(identifier)
        if chosen is None:
            return
        result = saved_result(database, user, chosen)
        if not result:
            st.warning("This saved assessment is no longer accessible for the selected model. Refresh the SQL view.")
            return
        result = {**result, "_studio_target_label": entry["outcome"].get("target_reference_label") or entry["outcome"].get("target_name")}
        render_inspector(result, database=database, user=user, key=key)
    except (OSError, ValueError, PermissionError) as exc:
        st.warning("Could not open the saved assessment: " + str(exc))


def render_saved_row_request(database, user, index, entries, *, context):
    request = st.session_state.pop("binderscope-row-request", None)
    if not request or request["context"] != context:
        return
    entry = next((entry for entry in entries if entry["id"] == request["entry_id"]), None)
    if entry is None or entry["outcome"]["outcome_kind"] != "Design":
        return
    token = hashlib.sha256(json.dumps([context, entry["id"], request["assessments"]]).encode()).hexdigest()[:24]
    _saved_row_dialog(database, user, index, entry, request["assessments"], "binderscope-row-" + token)


def render_saved_browser_action(database, user, index, entry, row, *, signature):
    from mn_molecular_design_studio.binderscope import saved_result
    from mn_molecular_design_studio.binderscope_inspector import inspection_dialog
    from mn_molecular_design_studio.binderscope_table import PREFIX, inspector_choices
    from mn_molecular_design_studio.data_browser import concise_result_label
    from mn_molecular_design_studio.data_browser_analyses import linked_analysis_columns

    if entry["outcome"]["outcome_kind"] != "Design":
        return
    choices = inspector_choices(row)
    token = hashlib.sha256(json.dumps([signature, entry["id"], [choice["observation_id"] for choice in choices]]).encode()).hexdigest()[:24]
    st.caption("Inspecting: " + concise_result_label(entry["outcome"]) + " · BinderScope: " + row.get(PREFIX + "Status", "Not assessed"))
    chosen = choices[0] if len(choices) == 1 else None
    if len(choices) > 1:
        labels = {choice["observation_id"]: choice for choice in choices}
        identifier = st.selectbox("Saved BinderScope assessment", list(labels), index=None,
            placeholder="Choose an assessment; none is selected automatically",
            format_func=lambda value: " · ".join(str(part) for part in (
                labels[value].get("analysis_version") or "Version not recorded",
                labels[value].get("run_id") or "Run not recorded", value)), key="binderscope-saved-choice-" + token)
        chosen = labels.get(identifier)
    if not choices and row.get(PREFIX + "Status") != "Not assessed":
        st.caption("No completed, verifiable saved report for this exact model. Other models and incomplete reports are not substituted.")
    if st.button("Open saved BinderScope · 3D & interface", key="binderscope-saved-open-" + token, disabled=chosen is None,
            help="Opens the focused design’s saved assessment. Does not start analysis or change admission."):
        try:
            fresh = linked_analysis_columns(database, index, [entry])[0]
            current = next((choice for choice in inspector_choices(fresh)
                if choice["observation_id"] == chosen["observation_id"] and choice["evaluation_id"] == chosen["evaluation_id"]), None)
            result = saved_result(database, user, current) if current else None
            if not result:
                st.warning("This saved assessment is no longer accessible for the selected model. Refresh the SQL view.")
                return
            result = {**result, "_studio_target_label": entry["outcome"].get("target_reference_label") or entry["outcome"].get("target_name")}
            inspection_dialog(result, database, user, "binderscope-table-" + token)
        except (OSError, ValueError, PermissionError) as exc:
            st.warning("Could not open the saved assessment: " + str(exc))


def render_browser_action(database, user, index, entries, selected, *, signature, workflow_schema=None, submit=None):
    designs = [entry for entry in entries if entry.get("outcome", {}).get("outcome_kind") == "Design"]
    if not designs:
        return
    selected_designs = [entry for entry in selected if entry.get("outcome", {}).get("outcome_kind") == "Design"]
    checked = [entry["id"] for entry in selected_designs]
    checked_key = hashlib.sha256(json.dumps(checked).encode()).hexdigest()[:12]
    scope_key = "binderscope-browser-selection-" + signature + "-" + checked_key
    scope = st.session_state.get(scope_key, "Selected table rows" if checked else "All matching designs")
    scope_control, launch_control = st.columns([1, 1], vertical_alignment="bottom")
    with scope_control.popover(f"{len(selected_designs):,} checked designs" if scope == "Selected table rows" else
            f"Whole view · {len(designs):,} designs", icon=":material/checklist:", width="stretch"):
        scope = st.radio("Which designs?", ["Selected table rows", "All matching designs"], key=scope_key,
            index=0 if selected_designs else 1,
            format_func=lambda value: f"All {len(designs):,} designs · all pages" if value == "All matching designs" else
                f"Only {len(selected_designs):,} checked designs")
        st.caption("The checkboxes at the far left of the results table choose a subset for analysis; they are not pass/fail indicators. "
            "Checking rows switches to that subset. Clear them to return to the whole view.")
        st.caption("The next screen skips already assessed complexes by default and lets you change that. Nothing runs until you confirm Queue.")
    chosen = selected_designs if scope == "Selected table rows" else designs
    if launch_control.button("Analyze · BinderScope", type="primary", disabled=not chosen, width="stretch",
            key="binderscope-browser-open-" + signature,
            help=f"Review {len(chosen):,} designs before queuing. Completed assessments are unchecked by default."):
        token = hashlib.sha256(json.dumps([user["id"], signature, [entry["id"] for entry in chosen]], sort_keys=True).encode()).hexdigest()[:24]
        _browser_dialog(database, user, index, chosen, token, workflow_schema=workflow_schema, submit=submit)


def render_setup(schema, *, user, project_selector, submit):
    st.header("BinderScope")
    st.caption("Integrated analysis source · revision " + REVISION[:7] + " · " + METHOD_VERSION)
    if _read_only_preview():
        st.info("Read-only development preview: configuration is available, but job submission is disabled.")
    st.write("Interface chemistry, liabilities, surface patches, structural diversity, binding modes and tag/fusion assessment. "
             "Reference context, specificity and scaffold evidence are enabled only when their inputs are supplied.")
    st.info("Start in Data browser: select table rows or all matching admitted designs at any tree level, then choose Analyze · BinderScope. "
            "Generator, topology and confidence settings are taken from the selected designs, with per-design overrides when needed.")
    if st.button("Open Data browser", key="binderscope-open-browser"):
        st.session_state["page"] = "Data Browser"
        st.rerun()
    with st.expander("Advanced · original campaign inputs"):
        _render_native_setup(schema, project_selector=project_selector, submit=submit)


def _render_native_setup(schema, *, project_selector, submit):
    parameters = {}
    st.caption("For original campaign funnels, rejected attempts and scaffold/reference evidence. This is not a second set of admitted designs, "
               "and running it does not import outputs into SQL.")
    project_id = project_selector("binderscope")
    directories = st.text_area("Native campaign folders · one absolute path per line", key="binderscope-native-folders")
    parameters["native_run_dirs"] = [value.strip() for value in directories.splitlines() if value.strip()]
    target = st.text_input("Original input target PDB / mmCIF", key="binderscope-native-target")
    if target:
        parameters["target_pdb"] = target
    for name, label in (("reference_cif", "Full-reference structure"), ("reference_chain", "Reference target chain"),
                        ("profile_path", "BinderScope profile JSON"), ("uniprot", "Target UniProt accession"),
                        ("reference_data_dir", "Offline BinderScope reference-data folder")):
        value = st.text_input(label, key="binderscope-native-" + name)
        if value:
            parameters[name] = value
    parameters["partner_chains"] = st.text_input("Reference partner chains (explicit comma-separated IDs, or none)", value="none", key="binderscope-partners")
    st.caption("The runtime is offline. Specificity requires supplied/cached records. ConSurf is not enabled.")
    ready = bool(project_id and parameters["native_run_dirs"] and target)
    parameters["cpu_cores"] = st.number_input("CPU cores", min_value=1, max_value=64, value=4, step=1, key="binderscope-cpu")
    st.warning("These are computational screens, not experimental binding or developability verdicts. "
               "This action does not admit, reject, delete or rename designs.")
    repeat = st.checkbox("Create a new assessment instead of reusing this submission", key="binderscope-repeat",
        help="Use this to retry a failed job or assess changed native campaign/reference files. Each assessment keeps its own provenance.")
    if st.button("Queue campaign assessment", type="primary", disabled=not ready or _read_only_preview(), key="binderscope-submit"):
        selection = parameters.get("binderscope_collection") or parameters.get("native_run_dirs")
        identity = uuid4().hex if repeat else settings_key(selection, parameters)
        submit(schema, "analysis.binderscope", parameters, project_id, {
            "_execution": {"idempotency_key": identity}, "prioritize": False,
        }, "BinderScope")


def render_design_result(result, *, database=None, user=None, key="", compact=False):
    from mn_molecular_design_studio.binderscope_inspector import inspector_key, inspection_dialog, render_inspector, render_summary

    if compact:
        render_summary(result)
        if st.button("Open BinderScope inspector · 3D & interface", key=inspector_key(result, key) + "-open", use_container_width=True):
            inspection_dialog(result, database, user, key)
        return
    render_inspector(result, database=database, user=user, key=key)


def render_job_result(run_dir):
    path = Path(run_dir) / "artifacts/analysis/binderscope/binderscope_report.json"
    if not path.is_file():
        return False
    report = json.loads(path.read_text())
    st.subheader("BinderScope assessment")
    st.caption("Integrated engine " + report["upstream_revision"][:7] + " · " + report["mode"] + " inputs")
    for limitation in report.get("limitations") or []:
        st.caption(limitation)
    results = report.get("results") or []
    if results:
        lookup = {str(index): result for index, result in enumerate(results)}
        chosen = st.selectbox("Assessed design", list(lookup), format_func=lambda key:
            str(lookup[key].get("source", {}).get("label") or lookup[key]["candidate_id"]), key="binderscope-report-design-" + str(run_dir))
        user = st.session_state.get("_molecular_database_user")
        from mn_molecular_design_studio.data_browser import ReadOnlyLibrary
        database = ReadOnlyLibrary() if user else None
        render_design_result(lookup[chosen], database=database, user=user, key=str(run_dir))
    with st.expander("Cohort families, diversity and coverage"):
        st.json(report.get("cohorts") or [], expanded=False)
    archive = Path(run_dir) / "artifacts/analysis/binderscope-results.zip"
    if archive.is_file():
        st.download_button("Download complete analysis archive", archive.read_bytes(),
            file_name="binderscope-results.zip", mime="application/zip", key="binderscope-archive-" + str(run_dir))
    return True
