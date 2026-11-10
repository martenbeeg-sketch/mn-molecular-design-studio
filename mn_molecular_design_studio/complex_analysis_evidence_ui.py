"""Evidence-first presentation of existing Database work; no scientific execution."""
from __future__ import annotations

from collections import Counter
import hashlib
import json

import streamlit as st

from mn_molecular_design_studio.complex_analysis_catalog import (
    saved_evidence_page, saved_evidence_summary, saved_relationship_page,
    saved_execution_groups, relationship_structure_ids,
)


METHOD_LABELS = {
    "md.mmgbsa.gb": "MD · MM/GBSA", "md.mmgbsa.pb": "MD · MM/PBSA",
    "md.analysis": "MD · trajectory analysis", "md.production": "MD · production information",
}


def method_label(namespace: str) -> str:
    return METHOD_LABELS.get(namespace, namespace)


def evidence_row(observation: dict) -> dict:
    """Keep scopes and zero values intact; never collapse method-specific scores."""
    descriptor = observation["evaluation"]
    context = observation.get("context") or {}
    return {
        **{key: value for key, value in observation.get("metrics", {}).items()
           if isinstance(value, (str, int, float, bool)) or value is None},
        "Molecule / design": observation.get("molecule_name") or observation["molecule_id"],
        "Method": method_label(observation["namespace"]),
        "Model / pose": descriptor.get("model_id") or descriptor.get("pose_id") or "Not recorded",
        "Source run": observation["run_id"],
        "Study / replica": " · ".join(str(context[key]) for key in ("study_id", "replica_index")
                                     if context.get(key) is not None),
        "Production segment": context.get("source_production_run_id") or context.get("production_run_id") or "Not recorded",
        "Continuation of": context.get("continuation_of_run_id") or "Not recorded",
        "Evidence history": ("Superseded by " + str(context["endpoint_superseded_by_run_id"]))
                            if context.get("endpoint_superseded_by_run_id") else "No supersession recorded",
        "Target": descriptor.get("target_id") or "Not recorded",
        "Recorded status": descriptor.get("status") or "Not recorded",
        "Saved measurements": "; ".join(f"{key}: {value}" for key, value in observation.get("metrics", {}).items()
                                         if isinstance(value, (int, float, bool)))[:380] or "See recorded evidence",
        "Evaluation": descriptor["id"],
        "Observation": observation["id"],
    }


def open_saved_work(state, *, user_id: str, group_ids: list[str], label: str,
                    project_label: str = "", work_type: str = "", system_label: str = "", system_scope=None) -> None:
    """Navigation state only; source membership/acceptance is never rewritten."""
    state["complex-analysis-evidence-selection"] = {
        "user_id": user_id, "group_ids": list(dict.fromkeys(group_ids)), "label": label,
        "project_label": project_label, "work_type": work_type,
        "system_label": system_label, "system_scope": system_scope,
    }
    state.pop("complex-analysis-workflow-samples", None)
    state.pop("complex-analysis-workflow-selection", None)
    state["complex-analysis-assessment-open"] = False
    for key in list(state):
        if key.startswith("complex-analysis-cohort-"):
            del state[key]
    state["complex-evidence-view"] = "Measurements"
    state["complex-analysis-flow-step"] = "overview"


def render_work_overview(database, user, selection):
    """Explain the saved experiment before loading models or submitting work."""
    from mn_molecular_design_studio import complex_analysis_collection_ui as collection
    from mn_molecular_design_studio.complex_analysis_catalog import collection_source_contexts, collection_system_facets, retained_system_contexts

    groups = [database.group(key, user) for key in selection["group_ids"]]
    contexts = collection_source_contexts(database, user, groups)
    fallback_ids = [group["id"] for group in groups
                    if not collection_system_facets({**group, "scientific_context": contexts.get(group["id"], {})})]
    systems = retained_system_contexts(database, user, fallback_ids) if fallback_ids else {}
    stats = database.cohort_statistics(user, {"selected": selection["group_ids"]})["selected"]
    st.markdown("### Campaign / study overview")
    st.write(selection["label"])
    st.caption(" → ".join(value for value in (selection.get("project_label"), selection.get("system_label"),
                                               selection.get("work_type")) if value))
    for column, (label, value) in zip(st.columns(3), [
        ("Retained molecules / designs", stats["molecule_count"]),
        ("Retained models / poses", stats["structure_count"]),
        ("Saved source memberships", len(groups)),
    ]):
        column.metric(label, f"{value:,}")
    st.caption("These are retained Database inventories. Generated attempts and recorded acceptance stages are shown separately below. Opening this overview reads metadata only.")
    prepared = [{**group, "scientific_context": contexts.get(group["id"], {}),
                 "retained_system_contexts": systems.get(group["id"], []),
                 "summary": collection._design_collection_label(group)} for group in groups]
    cohorts = collection._database_analysis_cohorts(prepared)
    for cohort in cohorts:
        if cohort["activity"] in {"Protein design", "Ligand design"}:
            if len(cohorts) > 1:
                st.markdown("**" + cohort["name"] + "**")
            for column, (label, value) in zip(st.columns(3), collection._cohort_source_counts(cohort).items()):
                with column.container(border=True):
                    st.caption(label)
                    st.write(value)
    st.dataframe([{"Campaign / source": cohort["name"], "Grouping evidence": cohort["cohort_type"],
                   "Source jobs": cohort["source_runs"], "Generator / source engine": cohort["generators"],
                   **collection._cohort_source_counts(cohort)} for cohort in cohorts],
                 hide_index=True, use_container_width=True)
    st.markdown("#### Original target / input system")
    st.caption("The input target or assembly is distinct from generated, cofolded and independently predicted output models. Variants and models retain their own acceptance and Evaluation evidence.")
    displayed = set()
    for group in prepared:
        facets = collection_system_facets(group)
        for facet in facets:
            if facet["key"] in displayed:
                continue
            displayed.add(facet["key"])
            target = facet["target"]
            with st.expander(facet["label"], expanded=len(groups) == 1):
                st.caption("Recorded input: " + facet["input_kind"])
                sequences = target.get("sequences_by_chain") or {}
                for chain, sequence in sorted(sequences.items()):
                    st.caption(f"Recorded target chain {chain} · {len(sequence)} residues")
                    st.code(sequence, language=None)
                st.json(target, expanded=False)
    if not displayed:
        st.info("No exact original input target or assembly is recorded in the source summary. Saved evidence can still expose recorded target/preparation references; names do not supply missing identity.")
    with st.expander("Distributed source jobs, settings and provenance"):
        st.dataframe([{"Source job": (contexts.get(group["id"], {}).get("job_ref") or group["summary"]["run"]),
                       "Database membership": group["id"], "Source label": group["name"]} for group in prepared],
                     hide_index=True, use_container_width=True)
        st.json([{ "group_id": group["id"], "origin": group.get("origin"),
                   "saved_source_context": contexts.get(group["id"], {})} for group in prepared], expanded=False)
    st.caption("Existing saved measurements are available without coordinates. New compatible assessments require explicit model selection and submission; no calculations start here.")
    left, right = st.columns(2)
    left.button("Review evidence & assessments", type="primary", key="overview-open-evidence",
                on_click=lambda: st.session_state.update({"complex-analysis-flow-step": "analyze"}))
    right.button("Choose structures to visualize", key="overview-open-models",
                 on_click=lambda: st.session_state.update({"complex-analysis-flow-step": "compare"}))


def _relationship_label(row: dict) -> str:
    target = " · ".join(str(row[key]) for key in ("target_name", "target_pdb") if row.get(key))
    target = target or row["target_ref"] or "Target not recorded"
    if row["preparation_id"]:
        target += " · prep " + row["preparation_id"][-8:]
    return target


def render_relationship_choices(database, user, selection, retained_stats):
    """Navigate recorded systems/studies; never create a molecular identity."""
    group_ids = selection["group_ids"]
    token = hashlib.sha256(str([group_ids, selection.get("system_scope")]).encode()).hexdigest()[:10]
    st.markdown("#### Molecule–target systems and studies")
    query = st.text_input("Find molecule, target or study", key="relationship-search-" + token)
    page_key = token + hashlib.sha256(query.encode()).hexdigest()[:8]
    number = int(st.number_input("Relationship page", min_value=1, value=1, step=1, key="relationship-page-" + page_key))
    page = saved_relationship_page(database, user, group_ids, query=query, offset=(number - 1) * 25,
                                   system_scope=selection.get("system_scope"))
    totals, rows = page["totals"], page["items"]
    cards = [("Molecular identities", totals["molecules"] if totals["relationships"] or query else retained_stats["molecule_count"]),
             ("Target / preparation contexts", totals["target_contexts"]),
             ("Recorded studies", totals["studies"]) if totals["studies"] else ("Source relationships", totals["relationships"])]
    for column, (label, value) in zip(st.columns(3), cards):
        column.metric(label, f"{value:,}")
    st.caption("One molecule can participate in several target systems and studies. Each row keeps its recorded target/preparation and source membership; matching names do not merge studies. Counts apply to the current search.")
    if not rows:
        selection.pop("relationship", None)
        st.session_state.pop("complex-analysis-workflow-samples", None)
        st.session_state.pop("complex-analysis-workflow-selection", None)
        st.session_state["complex-analysis-assessment-open"] = False
        st.session_state["complex-analysis-evidence-selection"] = dict(selection)
        st.info("No saved evidence relationships on this page. Clear the search, choose another page, or inspect retained molecules and models below.")
        return None, []
    relationship_columns = (["Target / preparation", "Study", "Replicas", "Production segments", "Molecule / design", "Saved evidence"]
                            if any(row["study_id"] for row in rows) else
                            ["Molecule / design", "Target / preparation", "Saved evidence", "Study"])
    st.dataframe([{
        "Molecule / design": row["molecule_name"] or row["molecule_id"],
        "Target / preparation": _relationship_label(row),
        "Study": row["study_id"] or "No study recorded · source membership kept",
        "Replicas": row["replica_count"], "Production segments": row["production_count"],
        "Saved evidence": ", ".join(method_label(method) for method in row["methods"]),
        "Source membership": row["group_id"],
    } for row in rows], hide_index=True, use_container_width=True, column_order=relationship_columns,
        column_config={"Target / preparation": st.column_config.TextColumn(width="medium"),
                       "Study": st.column_config.TextColumn(width="medium"),
                       "Replicas": st.column_config.NumberColumn(width="small"),
                       "Production segments": st.column_config.NumberColumn(width="small")})
    chosen = st.selectbox("Open molecular system / study", list(range(len(rows))),
        format_func=lambda index: f"{rows[index]['molecule_name'] or rows[index]['molecule_id']} · {_relationship_label(rows[index])} · {rows[index]['study_id'] or rows[index]['group_id']}",
        key="relationship-selected-" + page_key + str(number))
    row = rows[chosen]
    relationship = row["relationship"]
    if selection.get("relationship") != relationship:
        st.session_state.pop("complex-analysis-workflow-samples", None)
        st.session_state.pop("complex-analysis-workflow-selection", None)
        st.session_state["complex-analysis-assessment-open"] = False
        for key in list(st.session_state):
            if key.startswith("complex-analysis-cohort-"):
                del st.session_state[key]
        selection["relationship"] = relationship
        st.session_state["complex-analysis-evidence-selection"] = dict(selection)
    with st.expander("Exact molecule, target, preparation and study references"):
        st.json(relationship)
    execution = saved_execution_groups(database, user, group_ids, relationship)
    return row, execution


def render_model_choices(database, user: dict, selection: dict, *, intent: str = "compare") -> None:
    """Read one metadata page. Coordinates are fetched only by the load action."""
    from mn_molecular_design_studio import complex_analysis_collection_ui as collection

    group_ids = selection["group_ids"]
    relationship = selection.get("relationship")
    token = hashlib.sha256(json.dumps([group_ids, relationship], sort_keys=True).encode()).hexdigest()[:12]
    st.markdown("#### Choose exact models to inspect or assess")
    st.caption("Coordinates are loaded only after you click an action. Docking and cofolding models share this view; their source scores remain separate.")
    page_number = int(st.number_input("Model page", min_value=1, value=1, step=1, key="evidence-model-page-" + token))
    page = database.list_structure_records(user, collection="all", group_ids=group_ids,
        molecule_ids=[relationship["molecule_id"]] if relationship else None,
        restrict_structure_groups=True, limit=25, offset=(page_number - 1) * 25)
    records = list(page.get("items") or [])
    options = {str(model["id"]): (record, model) for record in records
               for model in record.get("representations") or record.get("structures") or []}
    if relationship:
        linked = relationship_structure_ids(database, user, group_ids, relationship)
        options = {key: value for key, value in options.items() if key in linked}
    if not options:
        st.info("No exact representation link to this selected system/study is recorded on this page. Saved measurements remain available. This view does not substitute a model from another target or study.")
        return
    chosen = st.multiselect("Models / poses", list(options), default=list(options)[:1], max_selections=50,
        format_func=lambda key: collection._structure_label(*options[key]), key="evidence-models-" + token)
    st.caption(f"{len(records)} molecular records on this page · {len(options)} model / pose references. No coordinates loaded by this list.")
    left, right = st.columns(2)
    visualize = left.button("Inspect & visualize selected models", type="primary", disabled=not chosen,
                            key="evidence-visualize-" + token)
    assess = right.button("Assess selected models…", disabled=not chosen, key="evidence-assess-" + token)
    if visualize or assess:
        groups = [database.group(key, user) for key in group_ids]
        # Only an explicitly saved passing collection carries a collected pass.
        # A source run alone is not evidence that every model passed filtering.
        passed = all((group.get("origin") or {}).get("type") == "results_design_collection"
                     for group in groups)
        st.session_state["complex-analysis-assessment-open"] = bool(assess)
        collection._load_collection_structures(database, user, records, group_ids,
            selected_structure_ids=chosen, passing_cohort=passed,
            next_step="analyze" if assess else intent)


def render_saved_evidence(database, user: dict, selection: dict) -> None:
    group_ids = selection["group_ids"]
    groups = [database.group(key, user) for key in group_ids]
    summaries = saved_evidence_summary(database, user, group_ids)
    stats = database.cohort_statistics(user, {"selected": group_ids})["selected"]
    methods = Counter()
    for summary in summaries.values():
        methods.update(summary["methods"])
    st.markdown("### Saved evidence")
    st.write(selection.get("label") or "Selected scientific work")
    scope = " · ".join(value for value in (selection.get("project_label"), selection.get("work_type")) if value)
    st.caption((scope + " · " if scope else "") + f"Viewing {len(groups)} selected collection(s), not the entire Database.")
    row, execution = render_relationship_choices(database, user, selection, stats)
    relationship = row["relationship"] if row else None
    replica = None
    include_history = True
    if row:
        st.markdown("#### Evidence for the selected system / study")
        if row["historical_count"]:
            include_history = st.checkbox(f"Include explicitly superseded history ({row['historical_count']} saved records)",
                value=False, key="relationship-history-" + hashlib.sha256(str(relationship).encode()).hexdigest()[:10])
            st.caption("Only records explicitly marked as superseded are hidden by default. Other results are not automatically ranked or replaced.")
        replicas = list(dict.fromkeys(item["replica"] for item in execution if item["replica"] is not None))
        if replicas:
            replica = st.selectbox("Recorded replica", [None, *replicas],
                format_func=lambda value: "All recorded replicas" if value is None else str(value),
                key="relationship-replica-" + hashlib.sha256(str(relationship).encode()).hexdigest()[:10])
            with st.expander("Replicas, production segments and continuation references"):
                st.dataframe([{"Replica": item["replica"], "Production segment": item["production_run"] or "Not recorded",
                               "Recorded continuation parent(s)": item["recorded_parents"] or "Not recorded",
                               "Superseded evidence records": item["historical_count"]}
                              for item in execution], hide_index=True, use_container_width=True)
        scoped = saved_evidence_summary(database, user, group_ids, relationship=relationship,
                                        replica=replica, include_history=include_history)
        methods = Counter()
        for summary in scoped.values():
            methods.update(summary["methods"])
    view = st.radio("Evidence view", ["Measurements", "Models & sequences", "Sources & provenance"],
                    horizontal=True, key="complex-evidence-view")
    if view == "Sources & provenance":
        st.caption(f"Selected evidence scope: {sum(methods.values()):,} visible saved records across {len(methods)} method namespaces. Collection inventories below retain their full history.")
        st.json({"selected_relationship": relationship, "replica": replica, "include_superseded_history": include_history}, expanded=False)
        for group in groups:
            with st.expander(group["name"], expanded=len(groups) == 1):
                st.json({"group_id": group["id"], "project_id": group["project_id"],
                         "origin": group["origin"], "methods": summaries[group["id"]]["methods"]})
        return
    if view == "Models & sequences":
        render_model_choices(database, user, selection)
        # This explicit view can inspect sequence/chemical identities even when
        # no coordinate representation exists. Hydration is bounded to one page.
        component_groups = [relationship["group_id"]] if relationship else group_ids
        group_id = st.selectbox("Molecular identities from collection", component_groups,
            format_func=lambda key: next(group["name"] for group in groups if group["id"] == key),
            key="complex-evidence-sequence-group")
        page = int(st.number_input("Molecule page", min_value=1, value=1, step=1, key="complex-evidence-molecule-page"))
        records = database.list_molecules(user, group_id=group_id, collection="all", limit=10, offset=(page - 1) * 10,
                                         molecule_ids=[relationship["molecule_id"]] if relationship else None)
        for record in records.get("items") or []:
            with st.expander(record.get("name") or record["id"]):
                st.json({"molecule_id": record["id"], "components": record.get("components"),
                         "source_identity": record.get("origin")})
        return
    if not methods:
        st.info("No saved measurements in these memberships. Open Models & sequences to inspect retained models and explicitly assess compatible inputs.")
        return
    if not row:
        return
    method = st.selectbox("Recorded method", ["", *sorted(methods)],
        format_func=lambda key: "All recorded methods" if not key else f"{method_label(key)} · {methods[key]} records",
        key="complex-evidence-method-" + hashlib.sha256(str([relationship, replica, include_history]).encode()).hexdigest()[:10])
    st.caption("Scores keep their original method and scope. Missing values are not zero. MM/GBSA and MM/PBSA are separate methods; native confidence and affinity predictions are not interchangeable.")
    token = hashlib.sha256((str([group_ids, relationship, replica, include_history]) + method).encode()).hexdigest()[:10]
    page_number = int(st.number_input("Evidence page", min_value=1, value=1, step=1, key="evidence-page-" + token))
    page = saved_evidence_page(database, user, group_ids, namespace=method, offset=(page_number - 1) * 25,
                               relationship=relationship, replica=replica, include_history=include_history)
    items = page["items"]
    st.caption(f"{page['total']:,} matching saved observations · page {page_number} · up to 25 records per page")
    if not items:
        st.info("No observations on this page.")
        return
    st.dataframe([evidence_row(item) for item in items], hide_index=True, use_container_width=True,
                 column_order=["Method", "Study / replica", "Production segment", "Saved measurements", "Recorded status",
                               "Model / pose", "Source run", "Evidence history"],
                 column_config={"Source run": st.column_config.TextColumn(width="medium")})
    item_by_id = {item["id"]: item for item in items}
    selected = st.selectbox("Inspect saved evidence", list(item_by_id),
        format_func=lambda key: f"{method_label(item_by_id[key]['namespace'])} · {item_by_id[key]['run_id']} · {key[-8:]}",
        key="evidence-observation-" + token + str(page_number))
    item = item_by_id[selected]
    with st.container(border=True):
        st.markdown("#### Exact measurement and provenance")
        st.caption(f"Evaluation {item['evaluation']['id']} · molecule {item['molecule_id']}")
        st.json({"measurements": item["metrics"], "evaluation": item["evaluation"]}, expanded=False)
        with st.expander("Measurement definitions and existing metric mappings"):
            from mn_molecular_design_studio.output_vocabulary import vocabulary
            mapped = item.get("normalized_output") or {}
            fields = set((mapped.get("metrics") or {}).keys())
            definitions = [definition for definition in vocabulary()["measures"] if definition["key"] in fields]
            if definitions:
                st.dataframe([{key: definition[key] for key in ("label", "unit", "scope", "description")}
                              for definition in definitions], hide_index=True, use_container_width=True)
            st.json(mapped, expanded=False)
            st.caption("Unmapped native fields retain their recorded names and source context. No equivalence, unit or affinity calibration is inferred.")
        with st.expander("Recorded parameters, units, artifact references and source context"):
            st.json(item["context"])
    if stats["structure_count"]:
        st.button("Choose models to visualize or assess", key="evidence-go-models",
                  on_click=lambda: st.session_state.update({"complex-evidence-view": "Models & sequences"}))
    else:
        st.info("This saved work has no attached 3D model. Its measurements and provenance are available above; 3D inspection requires a separately retained, explicitly linked representation.")
