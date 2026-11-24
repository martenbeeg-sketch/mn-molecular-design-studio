"""Studio's browsable design capability matrix."""
from __future__ import annotations

import pandas as pd
import streamlit as st

from mn_molecular_design_studio.binder_pipeline import (
    RFD4_MODULAR_ROUTES, RFD4_REFOLD_ENGINES, compatible_routes,
)
from mn_molecular_design_studio.design_capabilities import catalogue


def render_binder_preset(key, routes, *, button_label="Open pipeline draft"):
    route = (st.selectbox("Workflow variant", routes, format_func=lambda item: item.variant, key=key + "-variant")
             if len(routes) > 1 else routes[0])
    st.caption(route.summary)
    options = {}
    if route.id in RFD4_MODULAR_ROUTES:
        options["soluble_mpnn"] = st.checkbox(
            "Add SolubleMPNN sequence redesign", key=key + "-soluble", help="Keep the original designs and protect their interface in redesign variants."
        )
    if route.id == "rfd4-external":
        options["refold_engine"] = st.selectbox(
            "Refolding engine", list(RFD4_REFOLD_ENGINES), format_func=RFD4_REFOLD_ENGINES.get, key=key + "-refold"
        )
    if route.id in {"rfd4-internal", "rfd4-external"}:
        options["interface_analysis"] = st.checkbox("Add ipSAE & Structure Analysis", key=key + "-analysis")
    if st.button(button_label, type="primary", key=key + "-create"):
        return {"route_id": route.id, "target_input": route.target_input, "binder_format": route.format, **options}
    return None


def render(key):
    data = catalogue()
    st.subheader("Design capabilities")
    group = st.segmented_control("Design task group", ["Protein binders", "Other design tasks"],
                                 default="Protein binders", key=key + "-group")
    columns = data["binder_formats"] if group != "Other design tasks" else data["other_tasks"]
    rows = [row for row in data["items"] if row["output"] in columns]
    selected = "all"
    if group != "Other design tasks":
        targets = {"all": "All target inputs", **data["target_inputs"]}
        selected = st.selectbox("Target input", list(targets), index=1, format_func=targets.get, key=key + "-target")
        if selected != "all":
            rows = [row for row in rows if row["target_input"] == selected]
    if not rows:
        st.info("No documented design workflows for this target input.")
        return
    access = data["access"]
    st.markdown(
        "Select a **Build pipeline** cell to open its preset. "
        "**Add canvas step** means the tool is in the canvas library."
    )
    by_engine = {}
    for row in rows:
        cells = by_engine.setdefault(row["engine"], {})
        available = cells.setdefault(row["output"], set())
        available.add(row["access"])
    matrix = pd.DataFrame(
        [{"Engine / checkpoint": engine, **{
            label: next((access[level] for level in access if level in cells.get(code, ())), "—")
            for code, label in columns.items()}}
         for engine, cells in sorted(by_engine.items(), key=lambda entry: entry[0].casefold())]
    )
    cell_styles = {access[level]: "color: #176044; background-color: #e7f5ee" for level in ("pipeline", "canvas")}
    styled = matrix.style.map(lambda value: cell_styles.get(value, ""), subset=list(columns.values()))
    selection = st.dataframe(
        styled, hide_index=True, width="stretch", on_select="rerun", selection_mode="single-cell",
        key=f"{key}-matrix-{group}-{selected}",
        column_config={"Engine / checkpoint": st.column_config.TextColumn(pinned=True)},
    ).selection.cells
    output = None
    if selection:
        row_index, column = selection[0]
        engine = matrix.iloc[row_index]["Engine / checkpoint"]
        output = next((code for code, label in columns.items() if label == column), None)
    else:
        engine = st.selectbox("Engine details", list(matrix["Engine / checkpoint"]),
                              index=None, placeholder="Choose an engine", key=f"{key}-{group}-{selected}-details")
    if engine:
        selected_rows = [row for row in rows if row["engine"] == engine]
        if output is None:
            outputs = list(dict.fromkeys(row["output"] for row in selected_rows))
            output = (st.selectbox("Format or task", outputs, format_func=columns.get, key=f"{key}-output-{engine}")
                      if len(outputs) > 1 else outputs[0])
        selected_rows = [row for row in selected_rows if row["output"] == output]
        if not selected_rows:
            st.caption("Choose a supported cell.")
            return None
        target_inputs = list(dict.fromkeys(row["target_input"] for row in selected_rows if row["access"] == "pipeline"))
        if target_inputs:
            st.markdown(f"**{engine} · {columns[output]}**")
            target_input = (st.selectbox("Target context", target_inputs, format_func=data["target_inputs"].get,
                                         key=f"{key}-{engine}-{output}-context")
                            if len(target_inputs) > 1 else target_inputs[0])
            routes = [route for route in compatible_routes(target_input, output) if route.engine == engine]
            return render_binder_preset(f"{key}-{target_input}-{engine}-{output}", routes)
        for row in selected_rows:
            condition = data["target_inputs"].get(row["target_input"], "")
            st.markdown(" · ".join(filter(None, [f"**{columns[row['output']]}**", condition, access[row["access"]]])))
            if row["note"]:
                st.caption(row["note"] + (f" [Project reference]({row['reference']})" if row["reference"] else ""))
            if row["access"] == "canvas":
                from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG

                st.caption(f"Add **{NODE_CATALOG[row['canvas_node']]['label']}** from the canvas tool library.")
