"""Streamlit host for the freely movable pipeline canvas."""
from __future__ import annotations

import copy
import hashlib
import json
import uuid
from pathlib import Path

import streamlit as st
import streamlit.components.v1 as components

from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG, clean_graph, fingerprint, list_definitions, list_versions, migrate_template_graph, new_graph, validate_graph
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json
from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path

_COMPONENTS = {}
_BUILD = Path(__file__).with_name("components") / "pipeline_editor" / "build"


def _component():
    files = sorted(p for p in _BUILD.rglob("*") if p.is_file())
    digest = hashlib.sha256(b"".join(p.read_bytes() for p in files)).hexdigest()[:12]
    if digest not in _COMPONENTS:
        _COMPONENTS[digest] = components.declare_component("pipeline_canvas_" + digest, path=str(_BUILD))
    return _COMPONENTS[digest]


def _matching_execution(graph, snapshot, execution):
    def keys(graph):
        plan = validate_graph(graph)
        nodes = {n["id"]: n for n in plan["graph"]["nodes"]}
        result = {}
        for node_id in plan["order"]:
            data = nodes[node_id]["data"]
            result[node_id] = fingerprint(data["kind"], data["params"], [result[p] for p in plan["incoming"][node_id]])
        return result
    try:
        current, previous = keys(graph), keys(snapshot)
        return {n: state for n, state in execution.get("nodes", {}).items() if current.get(n) == previous.get(n)}
    except (ValueError, KeyError):
        return {}


def _target_ancestor(graph, node_id):
    nodes = {node["id"]: node for node in graph.get("nodes", [])}
    pending = [edge["source"] for edge in graph.get("edges", []) if edge.get("target") == node_id]
    seen = set()
    while pending:
        source_id = pending.pop(0)
        if source_id in seen or source_id not in nodes:
            continue
        seen.add(source_id)
        source = nodes[source_id]
        if source["data"]["kind"] == "target":
            return source
        pending.extend(edge["source"] for edge in graph.get("edges", []) if edge.get("target") == source_id)
    return None


def _render_node_input(selected, owner_id, project_id, search_jobs, prefix, graph=None):
    from mn_molecular_design_studio.target_editor import TargetSource, _hash, _target_widget_keys, render_target_coldspots, render_target_editor, render_target_source
    from mn_molecular_design_studio.design_validation import design_candidates
    from mn_molecular_design_studio.generation_io import generated_candidates
    contexts = {}
    node_id, params = selected["id"], selected["data"]["params"]
    key = f"{prefix}-{node_id}"
    if selected["data"]["kind"] == "target":
        path = resolve_stored_path(params.get("target_pdb"), must_exist=True) if params.get("target_pdb") else None
        source = render_target_source(key_prefix=key, label="Target Structure", example_code="7QF4", allow_cif=True, upload_in_expander=True, owner_id=owner_id)
        if not source.path and path:
            source = TargetSource(path=path, name=path.stem, owner_id=owner_id)
        if source.path:
            chain_key = f"{key}-chains-" + _hash(str(source.path))
            st.session_state.setdefault(chain_key, [c.strip() for c in params.get("target_chains", "A").split(",") if c.strip()])
            st.session_state.setdefault(_target_widget_keys(key, source.path)["hotspots"], params.get("hotspots", ""))
        target = render_target_editor(source, key_prefix=key, label="Target Structure", show_structure_tools=False, rerun_scope="fragment")
        if target.source_path:
            coldspot_key = _target_widget_keys(key, target.source_path)["coldspots"]
            if coldspot_key not in st.session_state:
                st.session_state[coldspot_key] = params.get("coldspots", "")
            with st.expander("Coldspot residues (optional)", expanded=bool(params.get("coldspots"))):
                st.caption("Coldspots are inherited only by generator engines that support them; each generator node shows its support status.")
                render_target_coldspots(target, key_prefix=key, show_structure_tools=False, rerun_scope="fragment")
        if target.design_path and target.ready and target.valid_trim_ranges:
            contexts[node_id] = {"target_pdb": str(target.design_path), "target_chains": ",".join(target.chains), "hotspots": target.hotspots, "coldspots": target.coldspots}
            params.update(contexts[node_id])
    elif selected["data"]["kind"] == "binder_context":
        context_name = selected["data"].get("label", "Binder Design Context")
        context_subject = context_name.removesuffix(" Context")
        target_node = _target_ancestor(graph or {}, node_id)
        target_params = target_node["data"]["params"] if target_node else {}
        path = resolve_stored_path(target_params.get("target_pdb"), must_exist=True) if target_params.get("target_pdb") else None
        if not path:
            st.info("Connect this context to a Target Structure node and load its structure first.")
            return contexts
        source = TargetSource(path=path, name=path.stem, owner_id=owner_id)
        chain_key = f"{key}-chains-" + _hash(str(source.path))
        selected_chains = target_params.get("target_chains", "")
        if isinstance(selected_chains, str):
            selected_chains = [chain.strip() for chain in selected_chains.replace(";", ",").split(",") if chain.strip()]
        st.session_state[chain_key] = list(selected_chains or [])
        widget_keys = _target_widget_keys(key, source.path)
        st.session_state.setdefault(widget_keys["hotspots"], params.get("hotspots", ""))
        st.session_state.setdefault(widget_keys["coldspots"], params.get("coldspots", ""))
        target = render_target_editor(source, key_prefix=key, label=f"{context_subject} target", show_structure_tools=False, allow_trimming=False,
                                      chain_selection_label="Target chains (set in Target Structure)", chain_selection_disabled=True,
                                      residue_selection_title="Hotspot override",
                                      residue_selection_caption="Choose residues the binder should contact. Leave empty to inherit hotspots from Target Structure.",
                                      residue_selection_input_label="Hotspot override", rerun_scope="fragment")
        if target.source_path:
            with st.expander("Coldspot residues (optional)", expanded=bool(params.get("coldspots"))):
                st.caption("Choose residues the binder should avoid. Leave empty to inherit coldspots from Target Structure. Support depends on the connected generator.")
                render_target_coldspots(target, key_prefix=key, show_structure_tools=False, rerun_scope="fragment")
        if target.ready and target.valid_trim_ranges:
            contexts[node_id] = {"hotspots": target.hotspots, "coldspots": target.coldspots}
            params.update(contexts[node_id])
    else:
        project_design=selected['data']['kind']=='project_design_source'
        generated=selected['data']['kind'] in {'generation_source','backbone_source'}
        stage='generation.backbone' if selected['data']['kind']=='backbone_source' else 'generation.backbone_sequence'
        if project_design:
            from mn_molecular_design_studio.engine.workflow_engine.workflows.graph_pipeline import source_project_designs
            candidate_loader=lambda path: source_project_designs(path, project_id)[1]
        else:
            candidate_loader=(lambda path:generated_candidates(path,stage=stage)) if generated else design_candidates
        if not project_id:
            st.info("Select a project to choose its completed design jobs." if project_design else "Select a project to choose its completed generator jobs." if generated else "Select a project to choose its successful design jobs.")
            return contexts
        rows=[];cursor=None
        while True:
            page = search_jobs(project_id=project_id, task_group="design", statuses=["completed", "succeeded"], limit=1000, cursor=cursor)
            rows.extend(page["items"]);cursor=page.get("next_cursor")
            if not cursor:
                break
        eligible = {str(row["run_dir"]): (row, designs) for row in rows if (designs := candidate_loader(Path(row["run_dir"])))}
        if not eligible:
            st.info("This project has no completed jobs with structure-bearing design candidates." if project_design else "This project has no completed generator jobs with the required output type." if generated else "This project has no successful design jobs with filter-passing designs.")
            return contexts
        options = list(eligible)
        default = options.index(params["source_run_dir"]) if params.get("source_run_dir") in options else 0
        source_dir = st.selectbox("Completed design job" if project_design else "Generator job" if generated else "Successful design job", options, index=default, format_func=lambda p: f"{eligible[p][0].get('campaign_name') or eligible[p][0]['tool']} · {eligible[p][0].get('job_code')}", key=key+"-source")
        ids = [c["candidate_id"] for c in eligible[source_dir][1]]
        chosen = st.multiselect("Project design candidates" if project_design else "Generated designs" if generated else "Filter-passing designs", ids, default=[c for c in params.get("selected_candidate_ids", ids) if c in ids], key=key+"-designs-"+Path(source_dir).name)
        if project_design:st.caption("Choose structures from any completed generator or design job. Antibody designs are annotated by SAbR before AbMPNN redesign.")
        elif generated:st.caption("These designs have not been validated. Connect sequence design and structure prediction steps to assess them.")
        contexts[node_id] = {"source_run_dir": source_dir, "selected_candidate_ids": chosen}
        params.update(contexts[node_id])
    return contexts


def _open_node_editor(prefix, node_id):
    graph = st.session_state[prefix+"-graph"]
    node = next((n for n in graph["nodes"] if n["id"] == node_id), None)
    if node:
        st.session_state[prefix+"-editor"] = {"id": uuid.uuid4().hex, "node_id": node_id, "params": copy.deepcopy(node["data"]["params"]), "revision": 0}


def _close_node_editor(prefix):
    st.session_state.pop(prefix+"-editor", None)


def _show_node_editor(prefix, *, owner_id, project_id, submit_api, search_jobs):
    from mn_molecular_design_studio.pipeline_node_ui import render_node_form
    editor = st.session_state.get(prefix+"-editor")
    if not editor:
        return
    graph = st.session_state[prefix+"-graph"]
    selected = next((n for n in graph["nodes"] if n["id"] == editor["node_id"]), None)
    if selected is None:
        _close_node_editor(prefix)
        return

    @st.dialog("Configure pipeline step", width="large", on_dismiss=lambda: _close_node_editor(prefix))
    def dialog():
        working = copy.deepcopy(selected)
        working["data"]["params"] = copy.deepcopy(editor["params"])
        key = f"{prefix}-editor-{editor['id']}-{editor['revision']}"
        def render_input():
            return _render_node_input(working, owner_id, project_id, search_jobs, key, graph)
        values, label, action, replacement_kind = render_node_form(working, graph, key=key, render_input=render_input, owner_id=owner_id)
        if action == "replace":
            try:
                from mn_molecular_design_studio.pipeline_node_ui import replacement_node_data
                updated = copy.deepcopy(graph)
                node = next(n for n in updated["nodes"] if n["id"] == selected["id"])
                replacement_source = copy.deepcopy(selected)
                replacement_source["data"].update(params=values, label=label)
                replacement = replacement_node_data(replacement_source, replacement_kind)
                node["data"].update(replacement)
                updated = clean_graph(updated)
                validate_graph(updated)
                if project_id:
                    submit_api("/api/v1/pipelines/save", {"project_id":project_id, "graph":updated})
                st.session_state[prefix+"-graph"] = updated
                st.session_state[prefix+"-settings-update"] = {
                    "id":uuid.uuid4().hex, "node_id":node["id"], **replacement,
                }
                _close_node_editor(prefix)
                st.rerun(scope="app")
            except (ValueError, OSError, TypeError, KeyError) as exc:
                st.error(str(exc))
        if action == "reset":
            from mn_molecular_design_studio.pipeline_node_ui import node_form_defaults
            editor["params"] = node_form_defaults(working["data"]["kind"])
            editor["revision"] += 1
            st.rerun(scope="fragment")
        if action == "cancel":
            _close_node_editor(prefix)
            st.rerun(scope="app")
        if action == "apply":
            try:
                updated = copy.deepcopy(graph)
                node = next(n for n in updated["nodes"] if n["id"] == selected["id"])
                node["data"].update(params=values, label=label)
                updated = clean_graph(updated)
                if project_id:
                    submit_api("/api/v1/pipelines/save", {"project_id":project_id, "graph":updated})
                st.session_state[prefix+"-graph"] = updated
                st.session_state[prefix+"-settings-update"] = {"id":uuid.uuid4().hex, "node_id":node["id"], "params":values, "label":label}
                _close_node_editor(prefix)
                st.rerun(scope="app")
            except (ValueError, OSError, TypeError, KeyError) as exc:
                st.error(str(exc))
    dialog()


def _activate_graph(prefix, graph):
    """Widget callback: replace the draft before its widgets are instantiated."""
    st.session_state[prefix+"-graph"] = copy.deepcopy(graph)
    st.session_state[prefix+"-reset"] = uuid.uuid4().hex
    st.session_state[prefix+"-name"] = graph["name"]
    _close_node_editor(prefix)


def _activate_selected_template(prefix, templates):
    """Load the selected starter graph before rendering the canvas again."""
    selected = st.session_state.get(prefix + "-template")
    if selected in templates:
        _activate_graph(prefix, new_graph(templates[selected]))


def render_pipeline_builder(*, owner_id, project_selector, submit_api, search_jobs):
    st.markdown("<style>.block-container{max-width:1900px;padding-top:1.35rem} [data-testid=stMainBlockContainer]{padding-left:2rem;padding-right:2rem}</style>", unsafe_allow_html=True)
    project_col, template_col, new_col = st.columns([1.2, 1.1, .55], vertical_alignment="bottom")
    with project_col:
        project_id = project_selector("visual-pipeline")
    prefix = f"visual-pipeline-{project_id or 'unassigned'}"
    state_key = prefix+"-graph"
    definitions = list_definitions(owner_id, project_id) if project_id else []
    if state_key not in st.session_state:
        st.session_state[state_key] = copy.deepcopy(definitions[0]["graph"] if definitions else new_graph())
        st.session_state[prefix+"-reset"] = uuid.uuid4().hex
        st.session_state[prefix+"-name"] = st.session_state[state_key]["name"]
    pending=st.session_state.get('pipeline-pending-draft')
    if pending and pending['project_id']==project_id:
        _activate_graph(prefix,pending['graph'])
        st.session_state.pop('pipeline-pending-draft',None)
    current_graph = st.session_state[state_key]
    migrated_graph = migrate_template_graph(current_graph)
    if migrated_graph is not current_graph:
        _activate_graph(prefix, migrated_graph)
    templates = {"Foundry · RFD3 → MPNN → RF3": "foundry",
                 "Antibody CDRs · Existing designs → SAbR → AbMPNN → sequence analysis": "abmpnn_cdr_redesign",
                 "PXDesign · MPNN → AF2 + Protenix":"pxdesign","RFdiffusion · MPNN → AF2":"rfdiffusion",
                 "Genie3 · MPNN → AF2 template":"genie3","Protpardelle-1C · MPNN → ESMFold":"protpardelle",
                 "ESMFold2 Binder Design → ProteinMPNN → RF3":"esmfold2_binder",
                 "Miniprotein binders · all generators → ProteinMPNN → Boltz-2 + analysis":"all_generators_boltz2",
                 "VHH binders · BindCraft2 + BoltzGen + ESMFold2 → ProteinMPNN → Boltz-2":"vhh_design",
                 "RFD4 Protein Binder → Boltz-2 → interface analysis":"rfd4_protein_binder",
                 "RFD4 Small-Molecule Pocket Design → Boltz-2":"rfd4_small_molecule",
                 "RFD4 Enzyme Motif Scaffolding → Boltz-2":"rfd4_enzyme",
                 "Cofolding · Protein monomer":"cofold_protein",
                 "Cofolding · Protein / protein":"cofold_protein_protein",
                 "Cofolding · Protein / ligand":"cofold_protein_ligand",
                 "Cofolding · Protein / DNA":"cofold_protein_dna",
                 "Cofolding · Protein / RNA":"cofold_protein_rna",
                 "Cofolding · Protein / macrocyclic peptide":"cofold_protein_macrocyclic",
                 "Blank canvas": "blank"}
    template_widget_key = prefix + "-template"
    old_template_label = "All generators → ProteinMPNN → Boltz-2 + analysis"
    if st.session_state.get(template_widget_key) == old_template_label:
        st.session_state[template_widget_key] = "Miniprotein binders · all generators → ProteinMPNN → Boltz-2 + analysis"
    with template_col:
        template = st.selectbox("Start from a template", list(templates), key=prefix+"-template",
                                on_change=_activate_selected_template, args=(prefix, templates))
    with new_col:
        st.button("Reset canvas", key=prefix+"-new", use_container_width=True,
                  on_click=_activate_selected_template, args=(prefix, templates))
    graph = st.session_state[state_key]
    from mn_molecular_design_studio.pipeline_graph import TEMPLATE_REFERENCES
    if graph.get('template_id') in TEMPLATE_REFERENCES:
        url,description=TEMPLATE_REFERENCES[graph['template_id']]
        with st.expander('About this starting template'):
            st.write(description)
            st.markdown('[Author workflow recommendations]('+url+')')
    name_col, reuse_col, gpu_col = st.columns([1.6, 1.25, .75], vertical_alignment="bottom")
    graph["name"] = name_col.text_input("Pipeline name", key=prefix+"-name")
    reuse = reuse_col.checkbox("Reuse unchanged completed steps", value=True, key=prefix+"-reuse", help="New runs can reuse successful steps in this project. Layout changes do not invalidate them.")
    gpu = gpu_col.selectbox("GPU assignment", ["Automatic", "GPU 0", "GPU 1"], key=prefix+"-gpu")
    with st.expander("Saved pipelines & versions", expanded=False):
        if definitions:
            saved = st.selectbox("Saved draft", definitions, format_func=lambda d:f"{d['graph']['name']} · {d['versions']} saved versions", key=prefix+"-saved")
            st.button("Open saved draft", key=prefix+"-open", on_click=_activate_graph, args=(prefix, saved["graph"]))
        versions = list_versions(owner_id, project_id, graph["id"]) if project_id else []
        if versions:
            version = st.selectbox("Version history", versions, format_func=lambda d:f"{d['saved_at'][:19].replace('T',' ')} UTC · {d['revision'][:10]}", key=prefix+"-version-"+graph["id"])
            st.button("Restore version into draft", key=prefix+"-restore", on_click=_activate_graph, args=(prefix, version["graph"]))
        if not definitions and not versions:
            st.caption("Drafts save automatically after you select a project. Save version keeps an immutable copy.")
    contexts = {}
    if not project_id:
        st.info("Choose a project before saving a pipeline version or running it. Each pipeline run is recorded under its selected project.")
    last_job = st.session_state.get(prefix+"-last-job", {})
    execution = {}
    if last_job.get("run_dir"):
        last_path = Path(last_job["run_dir"])
        execution = _matching_execution(graph, read_json(last_path/"graph_snapshot.json"), read_json(last_path/"graph_execution.json"))
    value = _component()(graph=graph, catalog=NODE_CATALOG, reset_token=st.session_state[prefix+"-reset"], input_contexts=contexts, execution=execution, settings_update=st.session_state.get(prefix+"-settings-update"), key=prefix+"-canvas", default=None)
    if value:
        signature = json.dumps(value, sort_keys=True)
        if signature != st.session_state.get(prefix+"-last-value"):
            st.session_state[prefix+"-last-value"] = signature
            try:
                updated = clean_graph(value["graph"])
                # Input paths are authored through the trusted upload/project
                # controls, rather than accepted from component messages.
                known = {n["id"]: n for n in graph["nodes"]}
                for n in updated["nodes"]:
                    if n["data"]["kind"] in {"target", "binder_context", "design_source", "generation_source", "backbone_source", "project_design_source", "complexa_generator", "pxdesign", "rfdiffusion", "genie3", "protpardelle"}:
                        previous = known.get(n["id"])
                        n["data"]["params"] = copy.deepcopy(previous["data"]["params"] if previous and previous["data"]["kind"]==n["data"]["kind"] else {f["name"]:f["default"] for f in NODE_CATALOG[n["data"]["kind"]]["fields"]})
                updated["name"] = graph["name"]
                if project_id:
                    submit_api("/api/v1/pipelines/save", {"project_id": project_id, "graph": updated})
                st.session_state[state_key] = graph = updated
                action=value.get("action") or {}
                if action.get("id") and action["id"] != st.session_state.get(prefix+"-handled-action"):
                    st.session_state[prefix+"-handled-action"] = action["id"]
                    if action["type"] == "edit_node":
                        _open_node_editor(prefix, action.get("node_id"))
                    elif not project_id:
                        st.warning("Select a project before saving a version or running the pipeline.")
                    elif action["type"] == "save_version":
                        result=submit_api("/api/v1/pipelines/save", {"project_id":project_id,"graph":graph,"version":True})
                        st.success(f"Version {result['revision'][:10]} saved. Existing versions and runs are preserved.")
                    elif action["type"] == "run":
                        result=submit_api("/api/v1/pipelines/run", {"project_id":project_id,"graph":graph,"reuse_completed":reuse,"gpu_device":{"Automatic":"all","GPU 0":"0","GPU 1":"1"}[gpu]})
                        st.session_state[prefix+"-last-job"] = result["job"]
                        st.session_state["last_submitted_job"] = result["job"]
                        st.success(f"Pipeline queued · {result['job'].get('job_code') or result['job']['run_id']}")
            except (ValueError, OSError, TypeError, KeyError) as exc:
                st.error(str(exc))
    _show_node_editor(prefix, owner_id=owner_id, project_id=project_id, submit_api=submit_api, search_jobs=search_jobs)
    try:
        plan=validate_graph(graph, require_inputs=True)
        st.caption(f"✓ Connections and inputs ready · {len(plan['order'])} steps. Runs use a separate snapshot of this graph.")
    except ValueError as exc:
        st.caption(f"Before running: {exc}")
    if project_id:
        with st.expander("Pipeline runs", expanded=bool(last_job)):
            if st.button("Refresh run status", key=prefix+"-refresh"):
                st.rerun()
            rows=search_jobs(project_id=project_id, task_group="campaign", tools=["visual_pipeline"], limit=25)["items"]
            if not rows:
                st.caption("No runs yet. Editing and saving the canvas does not start any jobs.")
            for row in rows[:10]:
                run_dir=Path(row["run_dir"])
                snapshot=read_json(run_dir/"graph_snapshot.json")
                if snapshot and snapshot.get("id") != graph["id"]:
                    continue
                columns=st.columns([2,1,1])
                columns[0].write(f"{row.get('job_code')} · {row.get('campaign_name') or graph['name']}")
                columns[1].write(row.get("status", "queued"))
                columns[2].link_button("Open results", "/?job_ref="+row["run_id"])
                execution=read_json(run_dir/"graph_execution.json")
                if execution.get("nodes"):
                    st.dataframe([{"Step": next((n["data"]["label"] for n in snapshot.get("nodes", []) if n["id"]==nid),nid), "Status": state.get("status"), "Job":state.get("job_code", ""), "Candidates":state.get("candidate_count"), "Detail":state.get("error", "")} for nid,state in execution["nodes"].items()], hide_index=True, use_container_width=True)
