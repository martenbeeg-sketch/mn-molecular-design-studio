"""Campaign → exact retained models → the existing assessment workspace."""
from __future__ import annotations

import streamlit as st

from mn_molecular_design_studio import scientific_lineage as lineage


def render_campaign_browser(database, user):
    st.markdown("### Scientific campaigns")
    st.caption("Project → target/construct → campaign → designs and models. This view reads indexed metadata; original Results remain on their source workstation.")
    clause, parameters = lineage._access(user)
    with database.connection() as c:
        facets = c.execute("SELECT c.project_id,c.target_key,json_extract(c.summary,'$.targets') targets FROM scientific_campaigns c WHERE " + clause, parameters).fetchall()
    from mn_molecular_design_studio.workspace import load_organization
    project_names = {p["id"]: p.get("name") or p["id"] for p in load_organization().get("projects", [])}
    target_names = {}
    import json
    for row in facets:
        targets = json.loads(row["targets"] or "[]")
        if row["target_key"]:
            label = " · ".join(str(t.get("biological_target") or t.get("label") or t.get("identity") or "Recorded construct") for t in targets)
            target_names[row["target_key"]] = label + " · " + row["target_key"][-10:]
    projects = sorted({row["project_id"] for row in facets if row["project_id"]})
    project = st.selectbox("Project", ["", *projects], format_func=lambda v: project_names.get(v, v) if v else "All accessible projects", key="scientific-project")
    kind = st.selectbox("Scientific activity", ["design", "folding", "docking", "dynamics", "analysis", ""],
                        format_func=lambda v: v.title() if v else "All activities", key="scientific-kind")
    targets = sorted({row["target_key"] for row in facets if row["target_key"] and (not project or row["project_id"] == project)})
    target = st.selectbox("Exact target/construct", ["", *targets], format_func=lambda v: target_names.get(v, v) if v else "All recorded constructs", key="scientific-target")
    query = st.text_input("Find campaign", key="scientific-query")
    filters = (project, kind, target, query)
    if st.session_state.get("scientific-filters") != filters:
        st.session_state["scientific-offset"] = 0
        st.session_state["scientific-filters"] = filters
    offset = int(st.session_state.get("scientific-offset", 0))
    page = lineage.campaigns(database, user, kind=kind, project_id=project, target=target, query=query, offset=offset)
    if not page["items"]:
        st.info("No indexed campaign matches. Historical indexing is an explicit staging/administrator operation; opening this page does not scan Results or contact remote workstations.")
        return []
    st.dataframe([{"Campaign": r["name"], "Type": r["kind"], "Target/construct": target_names.get(r["target_key"], "Not recorded / mixed constructs"), "Jobs": r["summary"]["jobs"],
        "Workstations": ", ".join(r["summary"]["workstations"]), "Methods": ", ".join(r["summary"]["methods"]),
        "Indexed occurrences": r["summary"]["indexed_occurrences"],
        "Authoritative output rows": r["summary"]["authoritative_output_occurrences"],
        "Recorded acceptance (candidate/model rows)": "; ".join(f"{stage}: {count['passed']} pass / {count['rejected']} reject / {count['unresolved']} unresolved / {count['unrecorded']} no decision"
            for stage, count in r["summary"].get("acceptance_counts", {}).items()) or "No decisions recorded",
        "Retained representation references": r["summary"]["retained_representation_refs"],
        "Mapping warnings": len(r["summary"]["mapping_warnings"]), "Updated": r["updated_at"]} for r in page["items"]], hide_index=True)
    left, right = st.columns(2)
    if left.button("Previous campaigns", disabled=not offset, key="scientific-prev"):
        st.session_state["scientific-offset"] = max(0, offset - 25)
        st.rerun()
    if right.button("Next campaigns", disabled=page["next_offset"] is None, key="scientific-next"):
        st.session_state["scientific-offset"] = page["next_offset"]
        st.rerun()
    choices = {row["ref"]: row for row in page["items"]}
    selected = st.selectbox("Open campaign", ["", *choices], format_func=lambda v: choices[v]["name"] if v else "Choose a campaign", key="scientific-campaign")
    if not selected:
        return []
    output_offset = int(st.session_state.get("scientific-output-offset", 0))
    if st.session_state.get("scientific-open-ref") != selected:
        output_offset = 0
        st.session_state["scientific-output-offset"] = 0
        st.session_state["scientific-open-ref"] = selected
    detail = lineage.campaign_detail(database, user, selected, offset=output_offset)
    st.dataframe(detail["jobs"], hide_index=True, column_order=["name", "workstation", "task_group", "run_id", "method", "status"])
    with st.expander("Recorded counts, source links and mapping limitations"):
        st.dataframe(detail["campaign"]["summary"]["counts"], hide_index=True)
        st.json({"dependencies": detail["dependencies"], "targets": detail["campaign"]["summary"]["targets"],
                 "warnings": detail["campaign"]["summary"]["mapping_warnings"]})
    stages = {d.get("stage_id"): d.get("stage") for o in detail["outputs"] for d in o.get("acceptance", []) if d.get("stage_id")}
    for stage in detail["campaign"]["summary"].get("acceptance_counts", {}):
        stages.setdefault(stage, {"native": "Native generator acceptance", "sequence": "Sequence design acceptance", "shared_filter": "Studio filter", "final": "Final computational acceptance"}.get(stage, stage))
    stage = st.selectbox("Recorded acceptance cohort", ["final", *[s for s in stages if s != "final"], "all"],
        format_func=lambda v: "All outputs — diagnostic selection" if v == "all" else stages.get(v) or "Final computational pass", key="scientific-stage-" + selected)
    if stage == "final" and "final" not in stages:
        st.info("No final decision is recorded on this page. Explicitly choose a recorded earlier passing stage, or All outputs for diagnostics. Native passing is not final acceptance.")
    execution_refs = st.checkbox("Include worker/controller representation references (diagnostic)", key="scientific-execution-refs")
    model_scope = (selected, stage, execution_refs)
    if st.session_state.get("scientific-model-scope") != model_scope:
        output_offset = 0
        st.session_state["scientific-output-offset"] = 0
        st.session_state["scientific-model-scope"] = model_scope
    model_page = lineage.retained_models(database, user, selected, stage=stage, offset=output_offset, include_execution_refs=execution_refs)
    st.caption(f"{model_page['output_total']} output occurrence(s) in the selected source/acceptance scope; {len(model_page['items'])} retained references on this page. Owner-local availability is recorded evidence, not a promise of an offline local copy.")
    if st.button("Next output page", disabled=model_page["next_offset"] is None, key="scientific-output-next"):
        st.session_state["scientific-output-offset"] = model_page["next_offset"]
        st.rerun()
    if st.button("Previous output page", disabled=not output_offset, key="scientific-output-prev"):
        st.session_state["scientific-output-offset"] = max(0, output_offset - 100)
        st.rerun()
    with st.expander("Indexed outputs and artifact availability"):
        st.json(detail["outputs"])
    models = {r["structure_id"]: r for r in model_page["items"]}
    if not models:
        st.info("No exact authorized retained representation matches this cohort/page. Source outputs remain inspectable in My results; missing model mappings are not replaced by another structure.")
        return []
    st.dataframe([dict(candidate=r["candidate_id"], model=r["model_id"] or "not recorded", molecule=r["molecule_id"],
                       structure=r["structure_id"], sha256=r["sha256"], saved_evaluations=len(r["evaluations"])) for r in models.values()], hide_index=True)
    ids = st.multiselect("Retained structures/models", list(models), default=list(models),
        format_func=lambda v: f"{models[v]['candidate_id']} · {models[v]['model_id'] or 'model not recorded'} · {v[-10:]}",
        key="scientific-models-" + selected + "-" + str(output_offset) + "-" + stage + str(execution_refs))
    if st.button("Open saved models for analysis", type="primary", disabled=not ids, key="scientific-open-models"):
        # Requery the exact source references and ACL before hydrating coordinates.
        current = {r["structure_id"] for r in lineage.retained_models(database, user, selected, stage=stage, offset=output_offset, include_execution_refs=execution_refs)["items"]}
        if not set(ids) <= current:
            st.error("Selection or permissions changed. Refresh this campaign.")
            return []
        from mn_molecular_design_studio.complex_analysis_collection_ui import _load_selected_structures
        molecules = sorted({models[i]["molecule_id"] for i in ids})
        records = database.list_structure_records(user, collection="all", molecule_ids=molecules, limit=len(molecules))["items"]
        samples = _load_selected_structures(database, user, records, ids)
        if not samples:
            st.error("The selected structures could not be hydrated. Their exact source references remain unchanged.")
            return []
        st.session_state["complex-analysis-workflow-samples"] = samples
        st.session_state["complex-analysis-workflow-selection"] = {"user_id": user["id"], "structure_ids": ids,
            "molecule_ids": molecules, "collection": "all", "project_id": "", "group_id": "", "campaign_ref": selected}
        mode = "All designs" if stage == "all" else "Final accepted designs" if stage == "final" else "Passed at selected stage"
        for suffix in ("widget", "value"):
            st.session_state["complex-analysis-cohort-acceptance-mode-" + suffix] = mode
            st.session_state["complex-analysis-cohort-acceptance-stage-" + suffix] = stage
        st.session_state["complex-analysis-flow-step"] = "analyze"
        st.rerun()
    return []
