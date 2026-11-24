"""Browse design tasks and registered canvas operations by capability."""
from __future__ import annotations

from html import escape

import streamlit as st

from mn_molecular_design_studio.binder_pipeline import (
    RFD4_MODULAR_ROUTES, RFD4_REFOLD_ENGINES, compatible_routes,
)
from mn_molecular_design_studio.design_capabilities import catalogue


_OPERATION_CATEGORIES = {
    "Structure prediction & cofolding": {"Structure prediction", "Cofolding input"},
    "Docking & affinity": {"Docking"},
    "Molecular dynamics": {"Molecular dynamics"},
    "Library preparation": {"Ligand library"},
    "Analysis & validation": {"Analysis", "Antibody analysis", "Collection"},
}


_matrix = st.components.v2.component(
    "studio_design_capabilities",
    html='<div class="matrix" role="region" aria-label="Capability table" tabindex="0"></div>',
    css="""
        .matrix {
            overflow: auto;
            border: 1px solid var(--st-border-color, #ddd);
            border-radius: .5rem;
        }
        table { border-spacing: 0; width: 100%; font-size: .875rem; }
        th, td {
            border-bottom: 1px solid var(--st-border-color, #ddd);
            padding: .45rem .65rem;
            text-align: left;
        }
        thead th { min-width: 8rem; }
        th:first-child {
            position: sticky;
            left: 0;
            min-width: 9rem;
            background: var(--st-background-color, white);
        }
        tbody tr:last-child > * { border-bottom: 0; }
        button {
            font: inherit;
            color: inherit;
            text-align: left;
            background: transparent;
            border: 1px solid var(--st-border-color, #999);
            border-radius: .4rem;
            padding: .32rem .5rem;
            cursor: pointer;
        }
        button[data-pipeline] { color: #14543e; background: #e7f5ee; border-color: #8bb8a4; }
        button:hover, button:focus-visible { outline: 2px solid var(--st-primary-color, #267356); }
    """,
    js="""
        export default function({ data, parentElement, setTriggerValue }) {
            const matrix = parentElement.querySelector('.matrix');
            if (matrix.markup !== data) {
                matrix.innerHTML = data;
                matrix.markup = data;
            }
            matrix.onclick = (event) => {
                const button = event.target.closest('button');
                if (button) setTriggerValue('choose', [button.dataset.engine, button.dataset.output]);
            };
        }
    """,
)


def _matrix_markup(by_engine, columns, access, *, action_column=False):
    labels = ["Studio access"] if action_column else columns.values()
    headers = "".join(f'<th scope="col">{escape(label)}</th>' for label in labels)
    table = [f'<table><thead><tr><th scope="col">Engine / checkpoint</th>{headers}</tr></thead><tbody>']
    for engine, cells in sorted(by_engine.items(), key=lambda entry: entry[0].casefold()):
        table.append(f'<tr><th scope="row">{escape(engine)}</th>')
        for code, label in columns.items():
            level = next((level for level in access if level in cells.get(code, ())), None)
            if level is None:
                table.append("<td>—</td>")
                continue
            title = f"{access[level]} · {engine} · {label}"
            style = ' data-pipeline="true"' if level == "pipeline" else ""
            table.append(
                f'<td><button type="button" data-engine="{escape(engine)}" data-output="{escape(code)}"'
                f' aria-label="{escape(title)}"{style}>{escape(access[level])}</button></td>'
            )
        table.append("</tr>")
    return "".join([*table, "</tbody></table>"])


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


def render_design_tasks(key, open_pipeline):
    data = catalogue()
    available_outputs = {row["output"] for row in data["items"]}
    tasks = {"protein_binders": "Protein binders", **{
        code: label for code, label in data["other_tasks"].items() if code in available_outputs
    }}
    task = st.selectbox("Design task", list(tasks), format_func=tasks.get, key=key + "-task")
    other_tasks = task != "protein_binders"
    if other_tasks:
        columns = {task: tasks[task]}
        target = "all"
    else:
        columns = data["binder_formats"]
        targets = {"all": "All target inputs", **data["target_inputs"]}
        target = st.selectbox("Target input", list(targets), index=1, format_func=targets.get, key=key + "-target")
    rows = [row for row in data["items"] if row["output"] in columns
            and (target == "all" or row["target_input"] == target)]
    if not rows:
        st.info("No documented design workflows for this target input.")
        return
    access = data["access"]
    if any(row["access"] == "pipeline" for row in rows):
        st.caption("Choose Build pipeline to configure a workflow.")
    by_engine = {}
    for row in rows:
        cells = by_engine.setdefault(row["engine"], {})
        available = cells.setdefault(row["output"], set())
        available.add(row["access"])
    choice = _matrix(
        data=_matrix_markup(by_engine, columns, access, action_column=other_tasks), key=f"{key}-matrix-{task}-{target}",
        on_choose_change=lambda: None,
    ).get("choose")
    if not choice:
        return
    engine, output = choice
    selected_rows = [row for row in rows if row["engine"] == engine and row["output"] == output]
    if not selected_rows:
        return

    @st.dialog(f"{engine} · {columns[output]}")
    def details():
        target_inputs = list(dict.fromkeys(row["target_input"] for row in selected_rows if row["access"] == "pipeline"))
        if target_inputs:
            target_input = (st.selectbox("Target context", target_inputs, format_func=data["target_inputs"].get,
                                         key=f"{key}-{engine}-{output}-context")
                            if len(target_inputs) > 1 else target_inputs[0])
            routes = [route for route in compatible_routes(target_input, output) if route.engine == engine]
            preset = render_binder_preset(f"{key}-{target_input}-{engine}-{output}", routes)
            if preset:
                try:
                    open_pipeline(preset)
                except ValueError as exc:
                    st.error(str(exc))
            return
        for row in selected_rows:
            condition = data["target_inputs"].get(row["target_input"], "")
            st.markdown(" · ".join(filter(None, [f"**{columns[row['output']]}**", condition, access[row["access"]]])))
            if row["note"]:
                st.caption(row["note"] + (f" [Project reference]({row['reference']})" if row["reference"] else ""))
            if row["access"] == "canvas":
                from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG

                st.caption(f"Add **{NODE_CATALOG[row['canvas_node']]['label']}** from the canvas tool library.")

    details()


def render_capabilities(key, open_binder, open_step):
    from mn_molecular_design_studio.tool_connections import editor_tool_catalog

    st.subheader("Capabilities")
    category = st.selectbox("Capability", ["Design", *_OPERATION_CATEGORIES], key=key + "-capability")
    if category == "Design":
        render_design_tasks(key, open_binder)
        return

    operations = {kind: spec for kind, spec in editor_tool_catalog().items()
                  if spec["category"] in _OPERATION_CATEGORIES[category] and not spec.get("hidden")}
    table = ["<table><thead><tr><th scope='col'>Tool or step</th><th scope='col'>Inputs</th>"
             "<th scope='col'>Outputs</th><th scope='col'>Pipeline</th></tr></thead><tbody>"]
    for kind, spec in sorted(operations.items(), key=lambda item: item[1]["label"].casefold()):
        label = escape(spec["label"])
        inputs = escape(", ".join(port["label"] for port in spec["inputs"]) or "—")
        outputs = escape(", ".join(port["label"] for port in spec["outputs"]) or "—")
        table.append(f'<tr><th scope="row">{label}</th><td>{inputs}</td><td>{outputs}</td>'
                     f'<td><button type="button" data-engine="{escape(kind)}" '
                     f'aria-label="Set up {label}">Set up</button></td></tr>')
    choice = _matrix(data="".join([*table, "</tbody></table>"]), key=f"{key}-operations-{category}",
                     on_choose_change=lambda: None).get("choose")
    if not choice or choice[0] not in operations:
        return

    kind = choice[0]
    spec = operations[kind]

    @st.dialog(spec["label"])
    def details():
        st.write(spec["palette"]["description"])
        if st.button("Open pipeline draft", type="primary", key=f"{key}-open-operation-{kind}"):
            try:
                open_step(kind)
            except ValueError as exc:
                st.error(str(exc))

    details()
