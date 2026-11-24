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
    st.caption("Pipeline is available in the binder picker. Tool setup uses a Studio workflow. Upstream means the current Studio adapter does not expose that mode. A dash indicates no entry for this combination.")
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
    cell_styles = {"Pipeline": "color: #176044; background-color: #e7f5ee",
                   "Tool setup": "color: #234f81; background-color: #eaf2fc",
                   "Upstream": "color: #674d23; background-color: #fdf2df"}
    styled = matrix.style.map(lambda value: cell_styles.get(value, ""), subset=list(columns.values()))
    st.dataframe(styled, hide_index=True, width="stretch",
                 column_config={"Engine / checkpoint": st.column_config.TextColumn(pinned=True)})
    engine = st.selectbox("Engine details", list(matrix["Engine / checkpoint"]),
                          index=None, placeholder="Choose an engine", key=key + "-details")
    if engine:
        selected_rows = [row for row in rows if row["engine"] == engine]
        for row in selected_rows:
            condition = data["target_inputs"].get(row["target_input"], "")
            st.markdown(" · ".join(filter(None, [f"**{columns[row['output']]}**", condition,
                                                  access[row["access"]], f"`{row['workflow']}`"])))
            if row["note"]:
                st.caption(row["note"] + (f" [Project reference]({row['reference']})" if row["reference"] else ""))
