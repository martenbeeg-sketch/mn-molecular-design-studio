"""Evidence-first presentation of existing Database work; no scientific execution."""
from __future__ import annotations

from collections import Counter
import hashlib
import json

import streamlit as st

from mn_molecular_design_studio.complex_analysis_catalog import saved_evidence_page, saved_evidence_summary


METHOD_LABELS = {
    "md.mmgbsa.gb": "MD · MM/GBSA", "md.mmgbsa.pb": "MD · MM/PBSA",
    "md.analysis": "MD · trajectory analysis", "md.production": "MD · production information",
}


def method_label(namespace: str) -> str:
    return METHOD_LABELS.get(namespace, namespace)


def evidence_summary_cards(stats: dict, methods: dict, collection_count: int) -> list[tuple[str, int]]:
    """Describe the selected scientific work, rather than a mixed record total."""
    if methods and all(method.startswith("md.") for method in methods):
        return [(method_label(method).removeprefix("MD · ") + " records", count)
                for method, count in sorted(methods.items())]
    return [("Molecules / designs", stats["molecule_count"]),
            ("Saved structures", stats["structure_count"]),
            ("Evidence types", len(methods)), ("Selected collections", collection_count)]


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
        "Target": descriptor.get("target_id") or "Not recorded",
        "Recorded status": descriptor.get("status") or "Not recorded",
        "Saved measurements": "; ".join(f"{key}: {value}" for key, value in observation.get("metrics", {}).items()
                                         if isinstance(value, (int, float, bool)))[:380] or "See recorded evidence",
        "Evaluation": descriptor["id"],
        "Observation": observation["id"],
    }


def open_saved_work(state, *, user_id: str, group_ids: list[str], label: str,
                    project_label: str = "", work_type: str = "") -> None:
    """Navigation state only; source membership/acceptance is never rewritten."""
    state["complex-analysis-evidence-selection"] = {
        "user_id": user_id, "group_ids": list(dict.fromkeys(group_ids)), "label": label,
        "project_label": project_label, "work_type": work_type,
    }
    state.pop("complex-analysis-workflow-samples", None)
    state.pop("complex-analysis-workflow-selection", None)
    state["complex-analysis-assessment-open"] = False
    for key in list(state):
        if key.startswith("complex-analysis-cohort-"):
            del state[key]
    state["complex-evidence-view"] = "Measurements"
    state["complex-analysis-flow-step"] = "analyze"


def render_model_choices(database, user: dict, selection: dict, *, intent: str = "compare") -> None:
    """Read one metadata page. Coordinates are fetched only by the load action."""
    from mn_molecular_design_studio import complex_analysis_collection_ui as collection

    group_ids = selection["group_ids"]
    token = hashlib.sha256(json.dumps(group_ids, sort_keys=True).encode()).hexdigest()[:12]
    st.markdown("#### Choose exact models to inspect or assess")
    st.caption("Coordinates are loaded only after you click an action. Docking and cofolding models share this view; their source scores remain separate.")
    page_number = int(st.number_input("Model page", min_value=1, value=1, step=1, key="evidence-model-page-" + token))
    page = database.list_structure_records(user, collection="all", group_ids=group_ids,
        restrict_structure_groups=True, limit=25, offset=(page_number - 1) * 25)
    records = list(page.get("items") or [])
    options = {str(model["id"]): (record, model) for record in records
               for model in record.get("representations") or record.get("structures") or []}
    if not options:
        st.info("No structural models on this page. Saved measurements remain available without coordinates.")
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
    summaries_to_show = evidence_summary_cards(stats, methods, len(groups))
    for start in range(0, len(summaries_to_show), 4):
        row = summaries_to_show[start:start + 4]
        for card, (title, value) in zip(st.columns(len(row)), row):
            card.metric(title, f"{value:,}")
    st.caption("Counts describe saved evidence, including historical or superseded entries; they are not counts of independent experiments. Opening this view runs no calculations.")
    view = st.radio("Evidence view", ["Measurements", "Models & sequences", "Sources & provenance"],
                    horizontal=True, key="complex-evidence-view")
    if view == "Sources & provenance":
        st.caption(f"Technical inventory: {sum(methods.values()):,} underlying saved records across {len(methods)} method namespaces.")
        for group in groups:
            with st.expander(group["name"], expanded=len(groups) == 1):
                st.json({"group_id": group["id"], "project_id": group["project_id"],
                         "origin": group["origin"], "methods": summaries[group["id"]]["methods"]})
        return
    if view == "Models & sequences":
        render_model_choices(database, user, selection)
        # This explicit view can inspect sequence/chemical identities even when
        # no coordinate representation exists. Hydration is bounded to one page.
        group_id = st.selectbox("Molecular identities from collection", group_ids,
            format_func=lambda key: next(group["name"] for group in groups if group["id"] == key),
            key="complex-evidence-sequence-group")
        page = int(st.number_input("Molecule page", min_value=1, value=1, step=1, key="complex-evidence-molecule-page"))
        records = database.list_molecules(user, group_id=group_id, collection="all", limit=10, offset=(page - 1) * 10)
        for record in records.get("items") or []:
            with st.expander(record.get("name") or record["id"]):
                st.json({"molecule_id": record["id"], "components": record.get("components"),
                         "source_identity": record.get("origin")})
        return
    if not methods:
        st.info("No saved measurements in these memberships. Open Models & sequences to inspect retained models and explicitly assess compatible inputs.")
        return
    method = st.selectbox("Recorded method", ["", *sorted(methods)],
        format_func=lambda key: "All recorded methods" if not key else f"{method_label(key)} · {methods[key]} records",
        key="complex-evidence-method-" + hashlib.sha256(str(group_ids).encode()).hexdigest()[:10])
    st.caption("Scores keep their original method and scope. Missing values are not zero. MM/GBSA and MM/PBSA are separate methods; native confidence and affinity predictions are not interchangeable.")
    token = hashlib.sha256((str(group_ids) + method).encode()).hexdigest()[:10]
    page_number = int(st.number_input("Evidence page", min_value=1, value=1, step=1, key="evidence-page-" + token))
    page = saved_evidence_page(database, user, group_ids, namespace=method, offset=(page_number - 1) * 25)
    items = page["items"]
    st.caption(f"{page['total']:,} matching saved observations · page {page_number} · up to 25 records per page")
    if not items:
        st.info("No observations on this page.")
        return
    st.dataframe([evidence_row(item) for item in items], hide_index=True, use_container_width=True,
                 column_order=["Molecule / design", "Method", "Target", "Saved measurements", "Recorded status", "Model / pose", "Source run"],
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
