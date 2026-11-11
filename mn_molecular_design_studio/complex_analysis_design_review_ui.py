"""One decision workbench over saved design-model evidence."""
from __future__ import annotations

from collections import Counter
import hashlib
import json

import streamlit as st

from mn_molecular_design_studio.complex_analysis_design_review import (
    design_review_page, review_measurements, remember_review_selection,
)
from mn_molecular_design_studio.complex_analysis_collection import (
    filter_acceptance_cohort, acceptance_stage_options,
)


def _field_label(field):
    from mn_molecular_design_studio.output_vocabulary import field_label
    known = {"oracle_ppi.overall_iptm": "Overall ipTM", "oracle_ppi.iptm_coi": "Interface ipTM",
             "oracle_ppi.overall_plddt": "Overall pLDDT", "oracle_ppi.plddt_coi": "Interface pLDDT",
             "qc_interface_contact_pairs": "Interface atom contacts", "binder_residue_count": "Binder residues",
             "oracle_ppi.overall_pae": "Overall PAE", "interface_buried_surface_area_half_delta_sasa_angstrom2": "Interface BSA (Å²)",
             "common.confidence.generator.complex.iptm": "Generator ipTM (0–1)",
             "common.interface.native.rfd4_proteina.binder_residue_count_4a": "Binder contacts ≤4 Å"}
    if field in known:
        return known[field]
    label = field_label(field)
    label = label.replace("qc_interface_", "").replace("_", " ")
    return label if len(label) <= 36 else label[:33] + "…"


def _method_label(method):
    return {"rfd4_proteina": "RFD4", "bindcraft2": "BindCraft2", "biotite_interface": "Biotite",
            "validation_acceptance": "Saved validation", "binder_interface_chemistry_v1": "Interface chemistry",
            "sequence_qc": "Sequence Composition", "proteinsol": "Protein-Sol", "rosetta_metrics": "Rosetta",
            "ipsae": "ipSAE", "tag_fusion_suitability": "Tag & fusion suitability",
            "sequence_liability": "Sequence liabilities", "pep_patch": "PEP-Patch"}.get(method, method)


def _metric_priority(item):
    field = item["field"].lower()
    priorities = ("iptm", "i_ptm", "plddt", "contact_pairs", "buried_surface", "unsatisfied", "pae", "hydrophobic")
    return next((index for index, part in enumerate(priorities) if part in field), len(priorities))


def render_design_review(database, user, selection, *, requested_action="", workflow_schema=None, submit=None):
    """Review numeric evidence first; retrieve only explicitly selected coordinates."""
    from mn_molecular_design_studio import complex_analysis_collection_ui as collection
    import pandas as pd
    import plotly.express as px

    token = hashlib.sha256(json.dumps([user["id"], selection["group_ids"], selection.get("system_scope")], sort_keys=True).encode()).hexdigest()[:12]
    assessment_requested = requested_action == "analyze"
    inspection_requested = requested_action == "compare"
    st.markdown("### Shortlist & export" if requested_action == "shortlist" else "### Review & assess")
    st.write(selection["label"])
    st.caption("Review saved results, select designs, then check assessment coverage below or inspect their structures. New calculations require explicit submission; coordinates load only when requested.")
    with st.expander("Campaign totals, input target and source provenance"):
        if st.button("Show campaign context", key="review-context-" + token):
            st.session_state["review-context-open-" + token] = True
        if st.session_state.get("review-context-open-" + token):
            from mn_molecular_design_studio.complex_analysis_evidence_ui import render_work_overview
            render_work_overview(database, user, selection)
    page_number = int(st.number_input("Design page · 250 models per page", min_value=1, value=1,
                                    step=1, key="review-page-" + token))
    page = design_review_page(database, user, selection["group_ids"], offset=(page_number - 1) * 250)
    samples = page["items"]
    if selection.get("system_scope") is not None:
        from mn_molecular_design_studio.complex_analysis_catalog import saved_relationship_page, relationship_structure_ids
        # Resolve only the explicit selected context; no cross-target substitution.
        linked, offset = set(), 0
        while True:
            relationships = saved_relationship_page(database, user, selection["group_ids"], limit=100,
                offset=offset, system_scope=selection["system_scope"])
            for row in relationships["items"]:
                linked.update(relationship_structure_ids(database, user, selection["group_ids"], row["relationship"]))
            if relationships["next_offset"] is None:
                break
            offset = relationships["next_offset"]
        samples = [sample for sample in samples if sample["structure_id"] in linked]
    if not samples:
        st.info("No exact retained models on this page. Choose another page or return to campaign selection; source evidence remains unchanged.")
        return
    if requested_action == "shortlist":
        selected_ids = st.session_state.get("design-review-selected-" + token, [])
        selected = [s for s in samples if s["structure_id"] in selected_ids]
        if selected:
            from mn_molecular_design_studio.complex_analysis_shortlist_export import render_shortlist_review_export
            render_shortlist_review_export(user, selected, key_prefix="review-shortlist-" + token)
            return
        st.info("Select designs in the review table below, then open Shortlist & export to record decisions and rationale.")
    final_count = sum(any(d["stage_id"] == "final" and d["cohort_status"] == "passed" for d in sample["decisions"]) for sample in samples)
    evidence_count = sum(bool(sample["associated"]) for sample in samples)
    shortlisted = sum(sample["review_annotation"]["status"] == "Shortlisted for experimental testing" for sample in samples)
    for col, (label, value) in zip(st.columns(4), [("Saved models", len(samples)), ("Final passed", final_count),
                                                ("With saved evidence", evidence_count), ("Your shortlist", shortlisted)]):
        col.metric(label, value)
    st.caption("Acceptance counts refer to explicitly linked decisions on these exact models. They are not an inferred campaign-wide pass count.")
    controls = st.columns([1.2, 1.4, 1.3])
    with controls[0]:
        mode = collection._stateful_selectbox("Cohort", ["Final accepted designs", "Passed at selected stage", "All designs",
            "Rejected at selected stage", "Pending / unresolved at selected stage", "No recorded decision"],
            key="review-cohort-" + token, state_key="review-cohort-value-" + token, default="Final accepted designs")
    stages = acceptance_stage_options(samples)
    with controls[1]:
        stage = collection._stateful_selectbox("Recorded passing stage", [key for key, _ in stages] or [""],
            format_func=lambda key: dict(stages).get(key, "No decision recorded"),
            disabled="selected stage" not in mode, key="review-stage-" + token,
            state_key="review-stage-value-" + token, default=stages[0][0] if stages else "")
    with controls[2]:
        search = collection._stateful_text_input("Find a design", placeholder="Design, native ID or source run",
            key="review-search-" + token, state_key="review-search-value-" + token)
    if not final_count and mode == "Final accepted designs":
        st.info("No final computational pass is recorded on this page. Choose Passed at selected stage to review an earlier recorded pass, or All designs for diagnostics.")
    filtered = filter_acceptance_cohort(samples, mode, stage)
    filtered = [s for s in filtered if not search or search.casefold() in
                json.dumps([s["label"], s["candidate"]["candidate_id"], s["context"]], default=str).casefold()]
    series, values = review_measurements(filtered)
    with st.expander("Measurement methods and additional fields"):
        include_raw = st.checkbox("Show additional raw numeric fields", value=False, key="review-raw-" + token)
        st.caption("Columns retain their method, version and parameter set. Different settings are separate choices; repeated conflicting measurements stay blank.")
    if not include_raw:
        series = {key: row for key, row in series.items() if any(part in row["field"].lower() for part in (
            "iptm", "i_ptm", "plddt", "ipae", "overall_pae", "contact_pairs", "buried_surface",
            "unsatisfied", "hydrophobic", "hydrogen_bond_count", "solubility", "liability", "binder_residue_count"))}
    ordered = sorted(series, key=lambda key: (-sum(key in value for value in values.values()), _metric_priority(series[key]), series[key]["field"]))
    protocols = sorted({(row["method"], row["version"], row["parameters_digest"]) for row in series.values()})
    protocol_number = {protocol: index + 1 for index, protocol in enumerate(protocols)}
    metric_label = lambda key: (f"{_field_label(series[key]['field'])} · {_method_label(series[key]['method'])}" +
        (f" · settings {protocol_number[(series[key]['method'], series[key]['version'], series[key]['parameters_digest'])]}" if len(protocols) > 1 else ""))
    defaults = []
    for key in sorted(ordered, key=lambda key: (_metric_priority(series[key]), -sum(key in value for value in values.values()))):
        if series[key]["field"] not in {series[other]["field"] for other in defaults}:
            defaults.append(key)
        if len(defaults) == 3:
            break
    metrics = st.multiselect("Measurements to compare", ordered, default=defaults,
                            format_func=metric_label, max_selections=6, key="review-metrics-" + token + hashlib.sha256(str(ordered).encode()).hexdigest()[:6])
    with st.expander("Filter by a saved measurement"):
        metric = st.selectbox("Measurement", ["", *metrics], format_func=lambda key: metric_label(key) if key else "No numeric filter", key="review-filter-" + token)
        if metric:
            recorded = [values[s["structure_id"]][metric] for s in filtered if metric in values.get(s["structure_id"], {})]
            if recorded:
                left, right = st.columns(2)
                low = left.number_input("Minimum", value=float(min(recorded)), key="review-min-" + token + metric)
                high = right.number_input("Maximum", value=float(max(recorded)), key="review-max-" + token + metric)
                keep_missing = st.checkbox("Keep designs without this measurement", value=True, key="review-missing-" + token + metric)
                filtered = [s for s in filtered if (keep_missing if metric not in values.get(s["structure_id"], {}) else
                            low <= values[s["structure_id"]][metric] <= high)]
    allowed = {s["structure_id"] for s in filtered}
    state_key = "design-review-selected-" + token
    selected_ids = remember_review_selection(st.session_state, token, st.session_state.get(state_key, []), allowed)
    st.caption(f"{len(filtered)} models shown · page {page_number} of {(page['total'] + 249) // 250}. Blank measurements mean missing, conflicting or failed evidence, not zero. Parameters and versions remain separate; no combined quality score.")
    if not filtered:
        return
    selected = [s for s in filtered if s["structure_id"] in selected_ids]
    st.markdown(f"#### {len(selected)} selected for review")
    if inspection_requested and not selected:
        st.info("Select designs in the table or plot below, then click Inspect selected in 3D. Use 2–5 compatible models for overlays.")
    if assessment_requested and not selected:
        st.info("Select designs in the table or plot below, then click Calculate missing assessments. The table shows saved evidence; opening it does not calculate new parameters.")
    actions = st.columns(3)
    assess = actions[0].button("Calculate missing assessments", type="primary", disabled=not selected,
                               key="review-assess-" + token)
    inspect = actions[1].button("Inspect selected in 3D", disabled=not selected, key="review-inspect-" + token)
    actions[2].button("Refresh saved results", key="review-refresh-" + token,
                      help="Reload saved Evaluations without resubmitting calculations.")
    if assess or inspect or ((inspection_requested or (assessment_requested and not
            st.session_state.get("complex-analysis-assessment-open"))) and selected):
        if len(selected) > 50:
            st.warning("Select at most 50 models for explicit structure loading; use 2–5 for overlays.")
        else:
            assessment = assess or (assessment_requested and not inspect)
            st.session_state["complex-analysis-assessment-open"] = assessment
            collection._load_collection_structures(database, user, [s["record"] for s in selected], selection["group_ids"],
                selected_structure_ids=selected_ids, passing_cohort=False, next_step="overview" if assessment else "compare",
                display_labels={s["structure_id"]: s["label"] for s in selected})
            # Loading can fail (for example an unavailable exact source file).
            # Keep the review visible without claiming an assessment was opened.
            st.session_state["complex-analysis-assessment-open"] = False
    if st.session_state.get("complex-analysis-assessment-open"):
        st.caption("Assessment settings are open below · [Go to coverage & settings](#assessment-coverage-and-settings). Compatible saved Evaluations are reused.")
    else:
        st.caption("Calculate missing assessments opens coverage and settings on this page for the selected exact models. Compatible saved Evaluations are reused.")
    # Labels only identify rows within this view; persistent IDs remain in samples.
    labels = {key: metric_label(key) for key in metrics}
    columns = {key: f"{_field_label(series[key]['field'])[:28]} [{i+1}]" for i, key in enumerate(metrics)}
    data = []
    for sample in filtered:
        decisions = sample["decisions"]
        decision = next((d for d in decisions if d["stage_id"] == (stage if "selected stage" in mode else "final")), None)
        data.append({"Design": ("✓ " if sample["structure_id"] in selected_ids else "") + sample["label"],
                     **{columns[key]: values.get(sample["structure_id"], {}).get(key) for key in metrics},
                     "Decision": decision["cohort_status"] if decision else "Not recorded",
                     "Engine": _method_label(sample["candidate"].get("source_tool") or "Not recorded"),
                     "Model": sample["model_id"] or "Saved model", "Evidence": len(sample["associated"]),
                     "Review": sample["review_annotation"]["status"].replace(" for experimental testing", "")})
    selection_version = hashlib.sha256(str([s["structure_id"] for s in filtered]).encode()).hexdigest()[:10]
    table_key = "review-table-" + token + selection_version
    def from_table():
        event = st.session_state[table_key]
        remember_review_selection(st.session_state, token,
            [filtered[i]["structure_id"] for i in event["selection"]["rows"] if 0 <= i < len(filtered)], allowed)
    left, right = st.columns([1.6, 1])
    with left:
        st.dataframe(pd.DataFrame(data), hide_index=True, width="stretch", height=380,
                     on_select=from_table, selection_mode="multi-row", key=table_key,
                     selection_default={"selection": {"rows": [i for i, sample in enumerate(filtered)
                         if sample["structure_id"] in selected_ids]}},
                     column_config={"Design": st.column_config.TextColumn(width="medium"),
                                    **{column: st.column_config.NumberColumn(width="small", format="%.3g") for column in columns.values()}},
                     column_order=["Design", *columns.values(), "Decision", "Engine"])
        st.caption("Select rows to review, inspect, assess or shortlist. Long IDs and method details are below, not in the main table.")
    with right:
        if len(metrics) >= 2:
            x = st.selectbox("X axis", metrics, format_func=metric_label, key="review-x-" + token)
            y = st.selectbox("Y axis", metrics, index=1, format_func=metric_label, key="review-y-" + token)
            plotted = [{"x": values[s["structure_id"]][x], "y": values[s["structure_id"]][y],
                        "Design": s["label"], "Engine": _method_label(s["candidate"].get("source_tool") or "Unknown"),
                        "Selected": s["structure_id"] in selected_ids, "Structure": s["structure_id"]}
                       for s in filtered if x in values.get(s["structure_id"], {}) and y in values.get(s["structure_id"], {})]
            if plotted:
                plot_key = "review-plot-" + token + selection_version + x + y
                def from_plot():
                    event = st.session_state[plot_key]
                    remember_review_selection(st.session_state, token,
                        [p["customdata"][0] for p in event["selection"]["points"] if p.get("customdata")], allowed)
                fig = px.scatter(pd.DataFrame(plotted), x="x", y="y", color="Engine", symbol="Selected",
                    hover_name="Design", custom_data=["Structure"], labels={"x": _field_label(series[x]["field"]), "y": _field_label(series[y]["field"])})
                fig.update_layout(height=280, margin=dict(l=5,r=5,t=8,b=5),
                                  showlegend=len({row["Engine"] for row in plotted}) > 1, legend=dict(orientation="h", title_text=""),
                                  clickmode="event+select", dragmode="select")
                st.plotly_chart(fig, on_select=from_plot, selection_mode=("points", "box", "lasso"), key=plot_key, width="stretch")
                st.caption(f"{len(plotted)} / {len(filtered)} models have both measurements. Select points or drag a box; no absent value is plotted as zero.")
            else:
                st.info("No model has both selected measurements with these saved settings. Choose different columns.")
        else:
            counts = Counter(row["Decision"] for row in data)
            st.bar_chart(pd.DataFrame({"Models": counts}), height=280)
            st.caption("Recorded decisions only. Select at least two numeric measurements for an evidence plot.")
    selected_ids = st.session_state.get(state_key, [])
    selected = [s for s in filtered if s["structure_id"] in selected_ids]
    if not selected:
        st.info("Select one or more designs above to inspect their evidence or calculate missing assessments.")
        return
    chosen = st.selectbox("Focused design", selected_ids, format_func=lambda key: next(s["label"] for s in selected if s["structure_id"] == key), key="review-focus-" + token + selection_version)
    focused = next(s for s in selected if s["structure_id"] == chosen)
    if st.button("Open focused Scientific Inspector", key="review-inspector-" + token):
        from mn_molecular_design_studio.complex_analysis_ui import open_database_complex
        open_database_complex(focused["molecule_id"], structure_id=focused["structure_id"])
        st.rerun()
    evidence_tab, shortlist_tab, provenance_tab = st.tabs(["Design evidence", "Shortlist & export", "Exact identities & methods"])
    with evidence_tab:
        st.caption(focused["label"] + " · " + (focused["candidate"].get("source_tool") or "Engine not recorded"))
        st.write("Saved evidence: " + (", ".join(sorted({_method_label(o["namespace"]) for o, _ in focused["associated"]})) or "Unavailable"))
        st.caption("Recorded chains · binder: " + (", ".join(focused["candidate"].get("binder_chains") or []) or "not recorded") +
                   " · target: " + (", ".join(focused["candidate"].get("target_chains") or []) or "not recorded"))
        if focused["candidate"].get("binder_sequence"):
            st.code(focused["candidate"]["binder_sequence"], language=None)
        st.dataframe([{"Measurement": labels[key], "Value": values.get(chosen, {}).get(key)} for key in metrics], hide_index=True, width="stretch")
        st.markdown("#### Saved assessment results")
        saved = sorted(focused["associated"], key=lambda item: (
            item[1].get("stage") == "analysis", str(item[1].get("created_at") or item[0].get("created_at") or "")), reverse=True)
        if saved:
            options = {descriptor["id"]: (observation, descriptor) for observation, descriptor in saved}
            result_id = st.selectbox("Saved result", list(options),
                format_func=lambda value: _method_label(options[value][0]["namespace"]) + " · " +
                    str(options[value][1].get("run_id") or value), key="review-saved-result-" + token + chosen)
            observation, descriptor = options[result_id]
            st.caption("Status: " + str(descriptor.get("status") or "not recorded") + " · Evaluation " + result_id +
                       " · exact structure SHA-256 " + focused["digest"])
            if observation["namespace"] == "tag_fusion_suitability" and not (
                    observation.get("context", {}).get("settings") or {}).get("termini_confirmed"):
                st.info("This saved tag/fusion assessment is provisional: complete intended construct termini were not confirmed. "
                        "Its geometry evidence is available; it is not experimental tag validation.")
            st.dataframe([{"Measurement (saved field)": key, "Saved value": json.dumps(value, ensure_ascii=False)}
                          for key, value in observation["metrics"].items() if isinstance(value, (int, float, str, bool)) or value is None],
                         hide_index=True, width="stretch")
            with st.expander("Selected result settings, units and residue-level evidence"):
                st.json({"evaluation": descriptor, "measurements": observation["metrics"],
                         "source_context": observation["context"]}, expanded=False)
            st.caption("Values retain their saved field names and method definitions. Open Scientific Inspector for residue-level evidence and the exact model.")
        else:
            st.info("No saved assessment is linked to this exact model yet. Use Calculate missing assessments to review coverage and explicitly request an analysis.")
        st.caption("Computational acceptance and experimental shortlisting are separate. Supporting decisions retain their original scope and policy.")
        with st.expander("Recorded acceptance policies and evidence"):
            st.json(focused["decisions"], expanded=False)
        with st.expander("All saved measurements, status and method limitations"):
            for observation, descriptor in focused["associated"]:
                st.caption(_method_label(observation["namespace"]) + " · " + str(descriptor.get("status") or "Not recorded"))
                st.json({"measurements": observation["metrics"], "evaluation": descriptor,
                         "source_context": observation["context"]}, expanded=False)
    with shortlist_tab:
        # Streamlit tabs execute eagerly: hydrate persistent review/export only
        # after an explicit request, avoiding per-design full-record reads.
        if st.button("Review shortlist decisions and exports", key="review-shortlist-" + token):
            st.session_state["review-shortlist-open-" + token] = True
        if st.session_state.get("review-shortlist-open-" + token):
            from mn_molecular_design_studio.complex_analysis_shortlist_export import render_shortlist_review_export
            render_shortlist_review_export(user, selected, key_prefix="review-shortlist-" + token)
    with provenance_tab:
        st.json({"molecule_id": focused["molecule_id"], "structure_id": chosen, "sha256": focused["digest"],
                 "model_id": focused["model_id"], "candidate": focused["candidate"],
                 "source": focused["context"], "evaluations": [descriptor for _, descriptor in focused["associated"]],
                 "displayed_measurements": {columns[key]: series[key] for key in metrics}}, expanded=False)

    if st.session_state.get("complex-analysis-assessment-open"):
        with st.container(border=True):
            st.markdown("#### Assessment coverage & settings")
            if st.button("Close assessment settings", key="review-close-assess-" + token):
                st.session_state["complex-analysis-assessment-open"] = False
                st.rerun()
            # Revalidate access on every rerun; never assess stale coordinates
            # from another campaign, model selection or user account.
            hydrated = collection._load_workspace_samples(user)
            expected = {(sample["molecule_id"], sample["structure_id"], sample["digest"]) for sample in selected}
            actual = {(sample["molecule_id"], sample["structure_id"], sample["digest"]) for sample in hydrated}
            if not hydrated or actual != expected:
                st.session_state["complex-analysis-assessment-open"] = False
                st.warning("Selected models changed or could not be revalidated. Click Calculate missing assessments to reopen their exact coverage.")
            else:
                st.caption(f"{len(hydrated)} selected model(s) · saved results remain above. Choose methods and settings, review coverage, then explicitly submit missing work. Changing the design selection closes these settings.")
                # Coordinates stay cached, but coverage must reflect newly
                # ingested Evaluations read by this page, not an old snapshot.
                evidence = {(sample["molecule_id"], sample["structure_id"]): sample["associated"] for sample in selected}
                assessment_samples = [{**sample, "associated": evidence[(sample["molecule_id"], sample["structure_id"])]}
                                      for sample in hydrated]
                collection._render_selected_assessments(user, assessment_samples, workflow_schema=workflow_schema, submit=submit)
