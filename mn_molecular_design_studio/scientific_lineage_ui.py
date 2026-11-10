"""Campaign → exact retained models → the existing assessment workspace."""
from __future__ import annotations

import json

import streamlit as st

from mn_molecular_design_studio import scientific_lineage as lineage


def _local_source_status():
    """Read coordinator heartbeat metadata only; never probe owner filesystems."""
    try:
        from mn_molecular_design_studio.cluster.store import Store
        return {row["id"]: {"availability": row.get("availability", "unknown"),
                            "last_seen_at": row.get("last_seen_at")}
                for row in Store().nodes()}
    except Exception:
        return {}


def _campaign_source_status(summary, nodes):
    refs = summary.get("source_job_refs") or []
    source_ids = sorted({ref.split("/", 3)[2] for ref in refs if ref.startswith("cluster://") and len(ref.split("/", 3)) > 2})
    if not source_ids:
        source_ids = list(summary.get("workstations") or [])
    statuses = []
    for node_id in source_ids:
        node = nodes.get(node_id)
        if not node:
            statuses.append(f"{node_id}: availability not known")
        elif node["availability"] == "online":
            statuses.append(f"{node_id}: heartbeat recent; detailed Results remain owner-local")
        else:
            last_seen = f", last seen {node['last_seen_at']}" if node.get("last_seen_at") else ""
            statuses.append(f"{node_id}: source temporarily unavailable ({node['availability']}{last_seen})")
    return "; ".join(statuses) or "Source availability not known"


def render_campaign_browser(database, user):
    st.markdown("### Scientific campaigns")
    st.caption("Project → target/construct → campaign → designs and models. This view reads indexed metadata; original Results remain on their source workstation.")
    clause, parameters = lineage._access(user)
    with database.connection() as c:
        facets = c.execute("SELECT c.project_id,ct.target_key,t.biological_identity_key,t.target_json "
            "FROM scientific_campaigns c LEFT JOIN scientific_campaign_targets ct ON ct.campaign_ref=c.ref "
            "LEFT JOIN scientific_target_constructs t ON t.target_key=ct.target_key WHERE " + clause, parameters).fetchall()
    from mn_molecular_design_studio.workspace import load_organization
    project_names = {p["id"]: p.get("name") or p["id"] for p in load_organization().get("projects", [])}
    target_names = {}
    for row in facets:
        if row["target_key"] and row["target_json"]:
            target_record = json.loads(row["target_json"])
            label = str(target_record.get("biological_target") or target_record.get("label")
                        or target_record.get("uniprot_id") or target_record.get("target_id") or "Recorded construct")
            target_names[row["target_key"]] = label + " · construct " + row["target_key"][-8:]
            if row["biological_identity_key"]:
                target_names[row["biological_identity_key"]] = label + " · all recorded constructs"
    projects = sorted({row["project_id"] for row in facets if row["project_id"]})
    project = st.selectbox("Project", ["", *projects], format_func=lambda v: project_names.get(v, v) if v else "All accessible projects", key="scientific-project")
    activity_labels = {"design": "Protein / molecular design", "folding": "Folding and prediction",
        "ligand_cofolding": "Ligand co-folding", "docking": "Docking and rescoring",
        "dynamics": "Molecular dynamics", "analysis": "Downstream analysis", "unknown": "Unclassified"}
    kind = st.selectbox("Scientific activity", ["", *activity_labels],
                        format_func=lambda v: activity_labels.get(v, "All activities"), key="scientific-kind")
    targets = sorted({key for row in facets if not project or row["project_id"] == project
                      for key in (row["target_key"], row["biological_identity_key"]) if key})
    target = st.selectbox("Biological target or exact construct", ["", *targets],
        format_func=lambda v: target_names.get(v, v) if v else "All recorded targets/constructs", key="scientific-target")
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
    source_nodes = _local_source_status()
    st.caption("Campaign counts and job relationships are synchronized metadata. Detailed Results stay on their owner workstation; an offline owner is shown as temporarily unavailable, and this page never scans or contacts remote workstations.")
    st.dataframe([{"Campaign": r["name"], "Primary activity": activity_labels.get(r["kind"], r["kind"]),
        "Activities": ", ".join(activity_labels.get(v, v) for v in r["summary"].get("activity_kinds", [])) or "Not recorded",
        "Target/construct": "; ".join(
            target_names.get(lineage._target_construct_key(t), str(t.get("biological_target") or t.get("label") or "Recorded construct"))
            for t in r["summary"].get("targets", [])) or "Not recorded / unresolved", "Jobs": r["summary"]["jobs"],
        "Workstations": ", ".join(r["summary"]["workstations"]), "Methods": ", ".join(r["summary"]["methods"]),
        "Source output occurrences": r["summary"].get("source_output_occurrences", 0),
        "Generated designs/backbones (normalized)": r["summary"].get("generated_design_occurrences", 0),
        "Candidate/model/pose occurrences": r["summary"].get("candidate_model_pose_occurrences", 0),
        "Output types": ", ".join(f"{name}: {count}" for name, count in sorted(r["summary"].get("output_kind_counts", {}).items())),
        "Recorded acceptance (candidate/model/pose rows)": "; ".join(f"{stage}: {count['passed']} pass / {count['rejected']} reject / {count['unresolved']} unresolved / {count['unrecorded']} no decision"
            for stage, count in r["summary"].get("acceptance_counts", {}).items()) or "No decisions recorded",
        "Database-retained candidate/models": r["summary"].get("database_retained_candidate_model_occurrences", 0),
        "Database-retained structures": r["summary"].get("database_retained_structure_references", 0),
        "Retained links synchronized": ("Complete" if r["summary"].get("publication", {}).get("complete") else
                             f"Pending · {r['summary'].get('publication', {}).get('indexed_occurrences', 0)}/"
                             f"{r['summary'].get('publication', {}).get('expected_occurrences', 0)} Database links"),
        "Source Results availability": _campaign_source_status(r["summary"], source_nodes),
        "Mapping warning counts": ", ".join(f"{name}: {count}" for name, count in r["summary"].get("mapping_warning_counts", {}).items()) or "None recorded",
        "Updated": r["updated_at"]} for r in page["items"]], hide_index=True)
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
    job_offset = int(st.session_state.get("scientific-job-offset", 0))
    dependency_offset = int(st.session_state.get("scientific-dependency-offset", 0))
    if st.session_state.get("scientific-open-ref") != selected:
        output_offset = 0
        st.session_state["scientific-output-offset"] = 0
        st.session_state["scientific-open-ref"] = selected
    detail = lineage.campaign_detail(database, user, selected, offset=output_offset, job_offset=job_offset,
                                     dependency_offset=dependency_offset)
    campaign_summary = detail["campaign"]["summary"]
    activities = set(campaign_summary.get("activity_kinds") or [detail["campaign"].get("kind", "unknown")])
    if activities & {"ligand_cofolding", "docking"}:
        st.info("Ligand co-folding and docking stay distinct in the source records. Open a retained protein–ligand model below to review PLIP/PandaMap interaction evidence and PoseBusters pose validation. Native docking/co-folding scores remain method-specific in the source Results; they are not combined with these analyses.")
    if detail["incomplete_job_count"]:
        st.warning(f"Retained Database links are still synchronizing for {detail['incomplete_job_count']} job(s). The source-authored campaign counts and acceptance totals remain available; exact retained models become selectable after all link pages arrive.")
    st.caption(f"Showing {min(job_offset + 1, detail['job_total'])}–{min(job_offset + len(detail['jobs']), detail['job_total'])} of {detail['job_total']} distributed execution jobs.")
    job_rows = [dict(row, source_availability=_campaign_source_status(
        {"source_job_refs": [row.get("ref", "")], "workstations": [row.get("workstation", "")]}, source_nodes))
        for row in detail["jobs"]]
    st.dataframe(job_rows, hide_index=True, column_order=["name", "workstation", "source_availability", "task_group", "run_id", "method", "status", "publication_state"])
    source_jobs = {row["ref"]: row for row in detail["jobs"] if row.get("ref") and row.get("task_group") and row.get("run_id")}
    if source_jobs:
        source_ref = st.selectbox("Open detailed source Results", ["", *source_jobs],
            format_func=lambda value: (f"{source_jobs[value].get('method') or source_jobs[value].get('tool') or 'Scientific job'} · "
                f"{source_jobs[value].get('workstation') or source_jobs[value].get('node_id')} · "
                f"{source_jobs[value].get('task_group')}/{source_jobs[value].get('run_id')}" if value else "Choose the exact source job"),
            key="scientific-source-job-" + selected)
        if st.button("Open full Results, native scores and unretained outputs", disabled=not source_ref,
                     key="scientific-open-source-results-" + selected):
            source = source_jobs[source_ref]
            from mn_molecular_design_studio.cluster import config as cluster_config
            local_reference = f"{source['task_group']}/{source['run_id']}"
            if str(source.get("node_id") or "") == cluster_config.identity():
                destination = local_reference
            else:
                cluster_id = str(source.get("cluster_job_id") or cluster_config.job_id(source["node_id"], local_reference))
                destination = "cluster:" + cluster_id
            st.session_state["selected-job-reference"] = destination
            st.session_state["page"] = "Job details"
            st.rerun()
        st.caption("This opens the existing authorized Results view for one selected source job; campaign browsing itself does not retrieve or scan its files.")
    job_left, job_right = st.columns(2)
    if job_left.button("Previous execution jobs", disabled=not job_offset, key="scientific-job-prev"):
        st.session_state["scientific-job-offset"] = max(0, job_offset - 100)
        st.rerun()
    if job_right.button("Next execution jobs", disabled=detail["next_job_offset"] is None, key="scientific-job-next"):
        st.session_state["scientific-job-offset"] = detail["next_job_offset"]
        st.rerun()
    with st.expander("Recorded stage counts, owner-local Results and mapping limitations"):
        st.dataframe(detail["campaign"]["summary"]["counts"], hide_index=True)
        st.json({"source_output_occurrences": detail["campaign"]["summary"].get("source_output_occurrences", 0),
                 "output_kind_counts": detail["campaign"]["summary"].get("output_kind_counts", {}),
                 "stage_occurrence_counts": detail["campaign"]["summary"].get("stage_occurrence_counts", {}),
                 "acceptance_counts": detail["campaign"]["summary"].get("acceptance_counts", {}),
                 "Database_retained_links": detail["campaign"]["summary"].get("database_retained_candidate_model_occurrences", 0),
                 "job_relationships": detail["dependencies"], "targets": detail["campaign"]["summary"]["targets"],
                 "mapping_warning_counts": detail["campaign"]["summary"].get("mapping_warning_counts", {})})
        st.caption(f"Showing up to 100 of {detail['dependency_total']} compact job/controller relationships. Per-candidate files, metrics and trajectory contents remain in source Results.")
        dep_left, dep_right = st.columns(2)
        if dep_left.button("Previous dependencies", disabled=not dependency_offset, key="scientific-dependency-prev"):
            st.session_state["scientific-dependency-offset"] = max(0, dependency_offset - 100)
            st.rerun()
        if dep_right.button("Next dependencies", disabled=detail["next_dependency_offset"] is None, key="scientific-dependency-next"):
            st.session_state["scientific-dependency-offset"] = detail["next_dependency_offset"]
            st.rerun()
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
    st.caption(f"{model_page['output_total']} exact retained Database link(s) match this selected source/acceptance scope; {len(model_page['items'])} authorized structures on this page. Full Results are owner-local and may be temporarily unavailable when that workstation is offline.")
    if st.button("Next output page", disabled=model_page["next_offset"] is None, key="scientific-output-next"):
        st.session_state["scientific-output-offset"] = model_page["next_offset"]
        st.rerun()
    if st.button("Previous output page", disabled=not output_offset, key="scientific-output-prev"):
        st.session_state["scientific-output-offset"] = max(0, output_offset - 100)
        st.rerun()
    with st.expander("Database-retained source links and artifact availability"):
        st.json(detail["outputs"])
    models = {r["structure_id"]: r for r in model_page["items"]}
    if not models:
        st.info("No exact authorized retained representation matches this cohort/page. Source outputs remain inspectable in My results; missing model mappings are not replaced by another structure.")
        return []
    st.dataframe([dict(candidate=r["candidate_id"], model=r["model_id"] or "not recorded", molecule=r["molecule_id"],
                       structure=r["structure_id"], sha256=r["sha256"], source_output=r["output_ref"],
                       source_availability=r["data_location"]["owner_results"],
                       database_representation=r["data_location"]["database_representation"],
                       saved_evaluations=len(r["evaluations"])) for r in models.values()], hide_index=True)
    st.caption("Source availability is the owner’s recorded indexing evidence. Remote reachability and coordinate-byte locality are not probed during campaign browsing; coordinates load only after you open the exact model.")
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
