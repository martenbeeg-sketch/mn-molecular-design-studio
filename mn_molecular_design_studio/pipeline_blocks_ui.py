"""A structured editor for the same saved pipelines used by the canvas."""
from pathlib import Path

import streamlit as st
import streamlit.components.v1 as components

from mn_molecular_design_studio.cluster.ui import _job_details_reference, _navigate_to_job_details
from mn_molecular_design_studio.pipeline_blocks import (
    SCREEN, add_block, blank_program, block_action, connect_blocks, next_blocks, program,
)
from mn_molecular_design_studio.pipeline_graph import (
    NODE_CATALOG, PipelineInputError, list_definitions, migrate_template_graph, new_graph, validate_graph,
)
from mn_molecular_design_studio.pipeline_graph_ui import (
    _activate_graph, _open_node_editor, _save_draft, _show_node_editor, _show_pipeline_save_as,
)
from mn_molecular_design_studio.scientific_tasks_ui import render_task_setup


_COMPONENT = components.declare_component(
    "studio_pipeline_program", path=str(Path(__file__).parent / "components/pipeline_blocks"),
)


def _start_or_open(prefix, owner_id, project_id, graph):
    drafts = list_definitions(owner_id, project_id) if project_id else []
    with st.expander("Open or start a pipeline", expanded=not graph["nodes"]):
        recipes, existing = st.tabs(["Starters", "Saved drafts"])
        with recipes:
            if st.button("Large-library screen", key=prefix + "-screen"):
                _activate_graph(prefix, new_graph(SCREEN))
                st.rerun()
            if graph["nodes"] and st.button("Blank program", key=prefix + "-blank"):
                _activate_graph(prefix, blank_program())
                st.rerun()
            created = render_task_setup(prefix + "-starter")
            if created is not None:
                _activate_graph(prefix, created)
                st.rerun()
        with existing:
            if not drafts:
                st.caption("Choose a project with saved drafts." if not project_id else "No saved drafts in this project.")
                return
            selected = st.selectbox("Saved draft", drafts, format_func=lambda draft: draft["graph"]["name"],
                                    key=prefix + "-saved")
            st.button("Open draft", key=prefix + "-open", on_click=_activate_graph,
                      args=(prefix, selected["graph"], selected["revision"]))


def _read_program(prefix, graph, display):
    node_ids = display["blocks"]
    selection_key = prefix + "-block"
    selected = st.session_state.get(selection_key)
    if selected not in node_ids:
        selected = display["screen"]["setup"][0] if "screen" in display else display["stages"][0][0]
    event = _COMPONENT(program=display, selected=selected, key=prefix + "-program-" + graph["id"], default=None)
    event_key = prefix + "-block-event"
    if event and event.get("token") != st.session_state.get(event_key):
        st.session_state[event_key] = event.get("token")
        if event.get("id") in node_ids and event["id"] != selected:
            st.session_state[selection_key] = event["id"]
            st.rerun()
    st.session_state[selection_key] = selected
    return selected


def _add_input(prefix, graph):
    choices = [kind for kind, spec in NODE_CATALOG.items() if spec["outputs"] and spec["category"] != "Active learning"
               and all(port.get("optional") for port in spec["inputs"])]
    choices.sort(key=lambda kind: (kind != "cofold_input", NODE_CATALOG[kind]["label"]))
    kind = st.selectbox("Input or standalone step", choices, format_func=lambda item: NODE_CATALOG[item]["label"],
                        key=prefix + "-input-kind")
    if st.button("Add block", icon=":material/add:", key=prefix + "-add-input"):
        updated, node_id = add_block(graph, kind)
        st.session_state[prefix + "-graph"] = updated
        st.session_state[prefix + "-block"] = node_id
        st.rerun()


def _next_block(prefix, graph, node):
    node_id = node["id"]
    outputs = NODE_CATALOG[node["data"]["kind"]]["outputs"]
    if not outputs:
        st.caption("This block collects the workflow output.")
        return
    names = {port["id"]: port["label"] for port in outputs}
    st.subheader("What comes next")
    output = (st.radio("Output branch", list(names), format_func=names.get, horizontal=True,
                       key=prefix + "-output-" + node_id) if len(outputs) > 1 else outputs[0]["id"])
    options, existing = next_blocks(graph, node_id, output)
    key = f"{prefix}-next-{node_id}-{output}"
    if options:
        actions = sorted({block_action(item[0])[0] for item in options})
        action = (st.selectbox("Next action", actions, key=key + "-action", index=None,
                               placeholder="All compatible actions") if len(actions) > 1 else None)
        options = sorted((item for item in options if action is None or block_action(item[0])[0] == action),
                         key=lambda item: (item[0] != "results", NODE_CATALOG[item[0]]["label"], item[1]))

        def option_label(option):
            kind, port = option
            spec = NODE_CATALOG[kind]
            name = spec["label"]
            if len(spec["inputs"]) > 1:
                name += " · " + next(item["label"] for item in spec["inputs"] if item["id"] == port)
            return name

        kind, target_port = st.selectbox("Compatible new block", options, format_func=option_label, key=key)
        if st.button("Add and connect", icon=":material/add:", key=key + "-add"):
            updated, target_id = add_block(graph, kind, source=node_id, output=output, target_port=target_port)
            st.session_state[prefix + "-graph"] = updated
            st.session_state[prefix + "-block"] = target_id
            st.rerun()
    else:
        st.caption("There are no compatible new blocks for this output.")
    if existing:
        with st.expander("Connect to a block already in this pipeline"):
            lookup = {item["id"]: item["data"] for item in graph["nodes"]}
            def destination_label(item):
                target, port = item
                data = lookup[target]
                label = next(entry["label"] for entry in NODE_CATALOG[data["kind"]]["inputs"] if entry["id"] == port)
                return f"{data['label']} · {label}"

            target_id, port = st.selectbox("Compatible destination", existing, key=key + "-existing",
                                          format_func=destination_label)
            if st.button("Connect blocks", key=key + "-connect"):
                st.session_state[prefix + "-graph"] = connect_blocks(graph, node_id, output, target_id, port)
                st.session_state[prefix + "-block"] = target_id
                st.rerun()


def _inspect_block(prefix, graph, selected):
    nodes = {item["id"]: item for item in graph["nodes"]}
    node = nodes[selected]
    data = node["data"]
    spec = NODE_CATALOG[data["kind"]]
    screen = graph.get("template_id") == SCREEN
    st.subheader(data["label"])
    st.caption("Compare the cumulative ranking with the previous round." if data["kind"] == "active_decision"
               else data.get("description") or spec["description"])
    if spec["inputs"] and not screen:
        st.write("Receives · " + ", ".join(port["label"] for port in spec["inputs"]))
    if spec["outputs"] and not screen:
        st.write("Produces · " + ", ".join(port["label"] for port in spec["outputs"]))
    st.button("Configure block", key=prefix + "-configure", on_click=_open_node_editor, args=(prefix, selected))

    if screen:
        edges = [edge for edge in graph["edges"] if edge["source"] == selected]
        if data["kind"] == "active_decision":
            edges.sort(key=lambda edge: nodes[edge["target"]]["data"]["kind"] != "active_results")
        for edge in edges:
            label = nodes[edge["target"]]["data"]["label"]
            lead = "Repeat at" if edge.get("loop_back") else "Then"
            if data["kind"] == "active_decision":
                lead = "If stopping" if nodes[edge["target"]]["data"]["kind"] == "active_results" else "Otherwise"
            st.write(f"{lead} · {label}")
        return

    _next_block(prefix, graph, node)
    with st.expander("Edit connections or remove block"):
        adjacent = [edge for edge in graph["edges"] if selected in (edge["source"], edge["target"])]
        if adjacent:
            chosen = st.selectbox("Connection", adjacent, key=f"{prefix}-edge-{selected}",
                format_func=lambda edge: f"{nodes[edge['source']]['data']['label']} → {nodes[edge['target']]['data']['label']}")
            if st.button("Disconnect", key=prefix + "-disconnect"):
                st.session_state[prefix + "-graph"] = {
                    **graph, "edges": [edge for edge in graph["edges"] if edge["id"] != chosen["id"]],
                }
                st.rerun()
        if st.button("Remove block", key=prefix + "-remove"):
            st.session_state[prefix + "-graph"] = {
                **graph, "nodes": [item for item in graph["nodes"] if item["id"] != selected],
                "edges": [edge for edge in graph["edges"] if edge not in adjacent],
            }
            st.rerun()


def render_pipeline_blocks(*, owner_id, project_selector, submit_api, search_jobs, read_only_projects=()):
    st.title("Pipeline Builder")
    project_id = project_selector("pipeline-blocks")
    if project_id in read_only_projects:
        st.info("Open archived project drafts in Pipelines.")
        return
    prefix = f"pipeline-blocks-{project_id or 'unassigned'}"
    graph_key = prefix + "-graph"
    if graph_key not in st.session_state:
        _activate_graph(prefix, blank_program())
    graph = migrate_template_graph(st.session_state[graph_key])
    _start_or_open(prefix, owner_id, project_id, graph)
    name, run, save, duplicate = st.columns([3.4, 1.2, 1, 1], vertical_alignment="bottom")
    graph = {**graph, "name": name.text_input("Pipeline name", value=graph["name"], max_chars=160,
        key=prefix + "-name-" + graph["id"]).strip() or "Untitled pipeline"}
    st.session_state[graph_key] = graph
    ready = False
    try:
        validate_graph(graph, require_inputs=True, require_files=False)
        ready = True
    except PipelineInputError as exc:
        notice, action = st.columns([3.4, 1.2], vertical_alignment="center")
        notice.info(str(exc))
        action.button("Configure required input", key=prefix + "-configure-input",
                      on_click=_open_node_editor, args=(prefix, exc.node_id))
    except ValueError as exc:
        if graph["nodes"]:
            st.info(str(exc))

    if graph["nodes"]:
        if run.button("Run pipeline", type="primary", disabled=not (ready and project_id), key=prefix + "-run"):
            try:
                result = submit_api("/api/v1/pipelines/run", {
                    "project_id": project_id, "graph": graph, "reuse_completed": True, "gpu_device": "all",
                })
                st.session_state[prefix + "-last-job"] = result["job"]
                st.session_state[prefix + "-last-job-graph"] = graph["id"]
                st.session_state["last_submitted_job"] = result["job"]
            except (ValueError, OSError, KeyError, TypeError) as exc:
                st.error(str(exc))
        if save.button("Save draft", disabled=not project_id, key=prefix + "-save"):
            try:
                _save_draft(prefix, project_id, graph, submit_api)
                st.success("Draft saved.")
            except (ValueError, OSError, KeyError, TypeError) as exc:
                st.error(str(exc))
        if duplicate.button("Save a copy", disabled=not project_id, key=prefix + "-copy"):
            st.session_state[prefix + "-save-as-open"] = True
    if st.session_state.get(prefix + "-last-job-graph") == graph["id"]:
        job = st.session_state[prefix + "-last-job"]
        reference = _job_details_reference(job)
        if st.button(f"Open run results · {job.get('job_code') or reference}", key=prefix + "-results"):
            _navigate_to_job_details(reference, final_results=graph.get("template_id") == SCREEN)
    if not project_id:
        st.caption("Choose a project to save or run this pipeline.")

    try:
        display = program(graph)
    except ValueError as exc:
        st.error(str(exc))
        labels = {node["id"]: node["data"]["label"] for node in graph["nodes"]}
        if labels:
            selected = st.selectbox("Block to inspect", list(labels), format_func=labels.get, key=prefix + "-repair")
            _inspect_block(prefix, graph, selected)
    else:
        workflow, inspector = st.columns([1.7, 1], gap="large")
        with workflow:
            if graph.get("template_id") != SCREEN:
                with st.expander("Add an input or standalone step", expanded=not graph["nodes"]):
                    try:
                        _add_input(prefix, graph)
                    except ValueError as exc:
                        st.error(str(exc))
            if graph["nodes"]:
                selected = _read_program(prefix, graph, display)
        if graph["nodes"]:
            with inspector:
                try:
                    _inspect_block(prefix, graph, selected)
                except ValueError as exc:
                    st.error(str(exc))
    _show_pipeline_save_as(prefix, project_id, submit_api)
    _show_node_editor(prefix, owner_id=owner_id, project_id=project_id,
                      submit_api=submit_api, search_jobs=search_jobs, stage_edits=True)
