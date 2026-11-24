"""Studio's browsable design capability matrix."""
from __future__ import annotations

import pandas as pd
import streamlit as st

from mn_molecular_design_studio.design_capabilities import catalogue


def render(key):
    data = catalogue()
    st.subheader("Design capabilities")
    group = st.segmented_control("Design task group", ["Protein binders", "Other design tasks"],
                                 default="Protein binders", key=key + "-group")
    columns = data["binder_formats"] if group != "Other design tasks" else data["other_tasks"]
    rows = [row for row in data["items"] if row["output"] in columns]
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
        "**Build pipeline** creates a connected draft. Configure its steps, then select **Run pipeline**. "
        "**Add canvas step** means the tool is in the canvas library. "
        "**Upstream only** is unavailable in Studio."
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
    cell_styles = {access["pipeline"]: "color: #176044; background-color: #e7f5ee",
                   access["canvas"]: "color: #176044; background-color: #e7f5ee",
                   access["upstream"]: "color: #674d23; background-color: #fdf2df"}
    styled = matrix.style.map(lambda value: cell_styles.get(value, ""), subset=list(columns.values()))
    st.dataframe(styled, hide_index=True, width="stretch",
                 column_config={"Engine / checkpoint": st.column_config.TextColumn(pinned=True)})
    engine = st.selectbox("Engine details", list(matrix["Engine / checkpoint"]),
                          index=None, placeholder="Choose an engine", key=key + "-details")
    if engine:
        selected_rows = [row for row in rows if row["engine"] == engine]
        outputs = list(dict.fromkeys(row["output"] for row in selected_rows))
        if len(outputs) > 1:
            output = st.selectbox("Format or task", outputs, format_func=columns.get, key=f"{key}-output-{engine}")
            selected_rows = [row for row in selected_rows if row["output"] == output]
        for row in selected_rows:
            condition = data["target_inputs"].get(row["target_input"], "")
            st.markdown(" · ".join(filter(None, [f"**{columns[row['output']]}**", condition, access[row["access"]]])))
            if row["note"]:
                st.caption(row["note"] + (f" [Project reference]({row['reference']})" if row["reference"] else ""))
            if row["access"] == "canvas":
                from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG

                st.caption(f"Add **{NODE_CATALOG[row['canvas_node']]['label']}** from the canvas tool library.")
