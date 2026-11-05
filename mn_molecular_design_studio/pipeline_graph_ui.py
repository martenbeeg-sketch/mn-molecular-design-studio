"""Streamlit host for the freely movable pipeline canvas."""
from __future__ import annotations

import copy
import hashlib
import json
import uuid
from pathlib import Path

import streamlit as st
import streamlit.components.v1 as components

from mn_molecular_design_studio.pipeline_graph import ACTIVE_LEARNING_LAYOUT, NODE_CATALOG, clean_graph, fingerprint, list_definitions, list_versions, migrate_template_graph, new_graph, validate_graph
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


def _restore_component_parameters(component_graph, previous_graph):
    """Keep node forms and trusted paths server-owned across canvas edits.

    The canvas edits graph topology and layout; workflow settings are edited by
    the Studio node form. Preserve them uniformly for every node kind instead
    of maintaining an incomplete per-tool whitelist.
    """
    trusted = clean_graph(previous_graph)
    updated = clean_graph(component_graph)
    known = {node["id"]: node for node in trusted.get("nodes", [])}
    for node in updated.get("nodes", []):
        previous = known.get(node["id"])
        if previous and previous.get("data", {}).get("kind") == node["data"].get("kind"):
            node["data"]["params"] = copy.deepcopy(node["data"].get("params", {})) | copy.deepcopy(
                previous.get("data", {}).get("params", {})
            )
    return updated


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


def _active_library_options():
    from mn_molecular_design_studio.molecular_database import MolecularDatabase
    result = []
    for record in MolecularDatabase().external_library_records():
        preparations = record.get("preparations") or []
        features = [item for item in preparations if item.get("schema_version") == "studio.library-fingerprints.v2"
                    and Path(str(item.get("path") or "")).is_file()]
        centroids = [item for item in preparations if item.get("schema_version") == "studio.library-centroids.v2"
                     and Path(str(item.get("path") or "")).is_file()]
        if features or centroids:
            result.append((record, features, centroids))
    return result


def _external_library_label(record, features=()):
    """Name a source by its registered database file, not by one sample run."""
    sources = record.get("sources") or []
    source_names = list(dict.fromkeys(
        Path(str(source.get("path") or "")).name
        for source in sources if source.get("path")
    ))
    name = ", ".join(source_names) if source_names else str(record.get("name") or record.get("id") or "External library")
    prepared_rows = max((int(item.get("rows") or 0) for item in features), default=0)
    if prepared_rows:
        name += f" · {prepared_rows:,} compounds prepared"
    return name


def _centroid_preparation_label(item):
    """Summarize a centroid output using its source population and method."""
    manifest = {}
    try:
        manifest = json.loads(Path(str(item.get("path") or "")).read_text(encoding="utf-8"))
    except (OSError, ValueError, TypeError):
        pass

    input_rows = int(manifest.get("input_rows") or item.get("input_rows") or 0)
    centroid_rows = int(manifest.get("row_count") or item.get("rows") or 0)
    algorithm = manifest.get("algorithm")
    methods = {
        "exact_greedy_molecular_weight_order": "exact greedy",
        "random_balanced_sampled_best_first": "chunked OpenVS-style, sampled",
        "random_balanced_hierarchical_best_first": "chunked OpenVS-style, hierarchical",
    }
    method = methods.get(algorithm, str(algorithm or "centroid selection").replace("_", " "))
    threshold = manifest.get("similarity_threshold")
    threshold_label = f" · Tanimoto {float(threshold):.2f}" if threshold is not None else ""
    created_at = str(manifest.get("created_at") or item.get("created_at") or "")
    date_label = f" · {created_at[:10]}" if created_at else ""
    return (f"{input_rows:,} input → {centroid_rows:,} centroids · {method}"
            f"{threshold_label}{date_label}")


def _render_active_library(params, key):
    libraries = [row for row in _active_library_options() if row[1]]
    if not libraries:
        st.info("No registered compound library has prepared fingerprints. Prepare a library first in the Library tools.")
        return None
    ids = [record["id"] for record, _features, _centroids in libraries]
    current_id = params.get("library_id")
    selected_id = current_id if current_id in ids else ids[0]
    selected_id = st.selectbox("Source compound database", ids,
        index=ids.index(selected_id),
        format_func=lambda value: next(_external_library_label(record, features)
            for record, features, _centroids in libraries if record["id"] == value),
        key=key + "-library")
    record, features, _centroids = next(row for row in libraries if row[0]["id"] == selected_id)
    feature_paths = [str(item["path"]) for item in features]
    feature_labels = {str(item["path"]): f"{int(item.get('rows') or 0):,} fingerprints · {Path(item['path']).parent.name}" for item in features}
    current_feature = params.get("feature_manifest")
    if current_feature not in feature_paths:
        current_feature = feature_paths[-1]
    feature_path = st.selectbox("Fingerprint dataset", feature_paths, index=feature_paths.index(current_feature),
        format_func=lambda value: feature_labels[value], key=key + "-features")
    st.caption("Prediction and random exploration use this compound library. The centroid input is selected independently in its own node.")
    return {"library_id": selected_id, "feature_manifest": feature_path}


def _render_active_centroids(params, key):
    libraries = [row for row in _active_library_options() if row[2]]
    if not libraries:
        st.info("No registered centroid datasets were found. Create centroids first in the Library tools.")
        return None
    ids = [record["id"] for record, _features, _centroids in libraries]
    selected_id = params.get("library_id") if params.get("library_id") in ids else ids[0]
    selected_id = st.selectbox("Source compound database", ids, index=ids.index(selected_id),
        format_func=lambda value: next(_external_library_label(record, features)
            for record, features, _centroids in libraries if record["id"] == value),
        key=key + "-centroid-library")
    record, _features, centroids = next(row for row in libraries if row[0]["id"] == selected_id)
    centroid_paths = [str(item["path"]) for item in centroids]
    centroid_labels = {str(item["path"]): _centroid_preparation_label(item) for item in centroids}
    current_centroid = params.get("centroid_manifest")
    if current_centroid not in centroid_paths:
        current_centroid = centroid_paths[-1]
    source_names = [Path(str(source.get("path"))).name for source in record.get("sources") or [] if source.get("path")]
    if source_names:
        st.caption(f"Centroid preparations registered for {', '.join(dict.fromkeys(source_names))}. Select the prepared centroid set to use as the random docking seed.")
    else:
        st.caption("Centroid preparations registered for this source database. Select the set to use as the random docking seed.")
    centroid_path = st.selectbox("Centroid preparation", centroid_paths,
        index=centroid_paths.index(current_centroid), format_func=lambda value: centroid_labels[value],
        key=key + "-centroids")
    # The initial sample-size widget is rendered after this input. Clear its
    # cached value when the user switches centroid sets so a value larger than
    # the newly selected set cannot remain in Streamlit's widget state.
    sample_widget_key = key + "-sample_size"
    selected_manifest_key = key + "-selected-centroid-manifest"
    if st.session_state.get(selected_manifest_key) != centroid_path:
        st.session_state.pop(sample_widget_key, None)
        st.session_state[selected_manifest_key] = centroid_path
    st.caption("This centroid source can differ from the compound database used for prediction. Both must use compatible fingerprint settings.")
    return {"library_id": selected_id, "centroid_manifest": centroid_path}


def _active_target_options():
    from mn_molecular_design_studio.workspace import _studio_home
    options = []
    registry = _studio_home() / "workdir" / "targets" / "mn_ligand_targets.json"
    try:
        payload = json.loads(registry.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        payload = {}
    for item in payload.get("targets", []) if isinstance(payload, dict) else []:
        path = Path(str(item.get("path") or "")).expanduser()
        if path.is_file():
            label = f"MN-Ligand · {item.get('pdb_id') or path.stem} · {item.get('name') or item.get('receptor') or 'imported receptor'}"
            options.append((label, str(path), item.get("name") or label))
    campaign = Path("/mnt/data/RESULTS/unidockpro/20251024_THRa_Large_screen_CD/mn_docking_data/THR_medchem_agonist/targets/thra/4lnw_T3_A_501_protein_mut_protein_refined.pdb")
    if campaign.is_file():
        campaign_ligand = campaign.with_name("4lnx_T3_A_502_ligand_refined.sdf")
        box = {
            "center": {"x": 26.616, "y": 38.249, "z": 32.711},
            "size": {"x": 30.0, "y": 20.0, "z": 20.0},
            "units": "angstrom",
        }
        if campaign_ligand.is_file():
            box["align_longest_axis_to_x"] = True
            box["alignment_ligand_sdf"] = str(campaign_ligand)
        options.append(("THRα campaign · 4LNW refined docking receptor", str(campaign), "THRα · 4LNW refined campaign receptor", box))
    user = st.session_state.get("_molecular_database_user")
    if user:
        from mn_molecular_design_studio.molecular_database import MolecularDatabase
        database = MolecularDatabase()
        records = database.list_molecules(
            user, scope="mine", collection="library", molecule_type="protein", limit=100_000
        )["items"]
        for record in records:
            protein_components = [row for row in record.get("components", []) if row.get("type") == "protein"]
            if not protein_components:
                continue
            component = protein_components[0]
            structure = next((row for row in record.get("structures", [])
                              if row.get("format") in {"pdb", "ent", "cif", "mmcif"}
                              and Path(str(row.get("path") or "")).is_file()), None)
            if structure is None:
                continue
            path = str(Path(structure["path"]).expanduser().resolve())
            metadata = record.get("metadata") or {}
            box = metadata.get("docking_box")
            label = f"Database · {record.get('name') or record['id']} · chain {component.get('chain') or '?'}"
            options.append((label, path, record.get("name") or label, box if isinstance(box, dict) else None))
    unique = {}
    for row in options:
        # Preserve docking-box metadata when the same structure appears in the
        # MN-Ligand registry and the reusable molecular database.
        previous = unique.get(row[1])
        if previous is None or (len(previous) < 4 and len(row) >= 4 and row[3]):
            unique[row[1]] = row
    return list(unique.values())


def _render_active_target(params, key):
    options = _active_target_options()
    current = str(params.get("target_pdb") or "")
    previous_target = current
    if current and Path(current).is_file() and current not in {item[1] for item in options}:
        options.insert(0, (params.get("target_name") or Path(current).stem, current, params.get("target_name") or Path(current).stem, None))
    if not options:
        st.info("No saved receptors were found. Import a receptor structure into Studio first.")
        return None
    paths = [item[1] for item in options]
    if current not in paths:
        current = paths[0]
    selected = st.selectbox("Existing receptor", paths, index=paths.index(current),
        format_func=lambda value: next(item[0] for item in options if item[1] == value), key=key + "-target")
    item = next(option for option in options if option[1] == selected)
    st.caption("Using a structure already imported into Studio. No upload or duplicate receptor entry is created.")
    result = {"target_pdb": selected, "target_name": item[2]}
    if selected != previous_target and len(item) >= 4 and isinstance(item[3], dict):
        box = item[3]
        center = box.get("center") or box.get("center_angstrom") or {}
        size = box.get("size") or box.get("size_angstrom") or {}
        if isinstance(center, (list, tuple)) and len(center) == 3:
            center = dict(zip(("x", "y", "z"), center))
        if isinstance(size, (list, tuple)) and len(size) == 3:
            size = dict(zip(("x", "y", "z"), size))
        if isinstance(center, dict) and isinstance(size, dict):
            for axis in "xyz":
                if center.get(axis) is not None:
                    result[f"center_{axis}"] = float(center[axis])
                if size.get(axis) is not None:
                    result[f"size_{axis}"] = float(size[axis])
        if "align_longest_axis_to_x" in box:
            result["align_longest_axis_to_x"] = bool(box["align_longest_axis_to_x"])
        if "alignment_ligand_sdf" in box:
            result["alignment_ligand_sdf"] = str(box["alignment_ligand_sdf"])
            result["reference_ligand_sdf"] = str(box["alignment_ligand_sdf"])
    return result


def _render_node_input(selected, owner_id, project_id, search_jobs, prefix, graph=None):
    if selected["data"]["kind"] == "run_source":
        from mn_molecular_design_studio.run_inputs_ui import select
        choice = select(project_id, prefix + "-run-output", require_complex=False,
                        reference=selected["data"]["params"].get("run_reference"))
        if choice:
            selected["data"]["params"]["run_reference"] = choice["reference"]
        return {}
    from mn_molecular_design_studio.target_editor import TargetSource, _hash, _target_widget_keys, render_target_coldspots, render_target_editor, render_target_source
    from mn_molecular_design_studio.design_validation import design_candidates
    from mn_molecular_design_studio.generation_io import generated_candidates
    contexts = {}
    node_id, params = selected["id"], selected["data"]["params"]
    key = f"{prefix}-{node_id}"
    if selected["data"]["kind"] == "target":
        path = resolve_stored_path(params.get("target_pdb"), must_exist=True) if params.get("target_pdb") else None
        source = render_target_source(
            key_prefix=key, label="Target Structure", example_code="7QF4",
            allow_cif=True, upload_in_expander=True, owner_id=owner_id,
            allow_saved_targets=True, project_id=project_id,
            initial_molecule_id=params.get("target_source_molecule_id", ""),
            initial_structure_id=params.get("target_source_structure_id", ""),
        )
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
            contexts[node_id] = {
                "target_pdb": str(target.design_path), "target_chains": ",".join(target.chains),
                "hotspots": target.hotspots, "coldspots": target.coldspots,
                "target_source_molecule_id": source.molecule_id,
                "target_source_structure_id": source.structure_id,
            }
            params.update(contexts[node_id])
    elif selected["data"]["kind"] == "active_library":
        chosen = _render_active_library(params, key)
        if chosen:
            contexts[node_id] = chosen
            params.update(chosen)
    elif selected["data"]["kind"] == "active_seed":
        chosen = _render_active_centroids(params, key)
        if chosen:
            contexts[node_id] = chosen
            params.update(chosen)
    elif selected["data"]["kind"] == "active_target":
        chosen = _render_active_target(params, key)
        if chosen:
            contexts[node_id] = chosen
            params.update(chosen)
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
            params['source_scope'] = st.selectbox("Source scope", ['designs', 'complex_results'],
                index=1 if params.get('source_scope') == 'complex_results' else 0,
                format_func=lambda value: 'Complex Results · any completed project job' if value == 'complex_results' else 'Design jobs', key=key+'-scope')
            candidate_loader=lambda path: source_project_designs(path, project_id,
                include_complex_results=params['source_scope'] == 'complex_results')[1]
        else:
            candidate_loader=(lambda path:generated_candidates(path,stage=stage)) if generated else design_candidates
        if not project_id:
            st.info("Select a project to choose its completed design jobs." if project_design else "Select a project to choose its completed generator jobs." if generated else "Select a project to choose its successful design jobs.")
            return contexts
        rows=[];cursor=None
        while True:
            page = search_jobs(project_id=project_id, task_group=None if project_design and params.get('source_scope') == 'complex_results' else "design", statuses=["completed", "succeeded"], limit=1000, cursor=cursor)
            rows.extend(page["items"]);cursor=page.get("next_cursor")
            if not cursor:
                break
        eligible = {}
        for row in rows:
            try:
                designs = candidate_loader(Path(row['run_dir']))
            except (ValueError, FileNotFoundError):
                continue
            if designs:
                eligible[str(row['run_dir'])] = (row, designs)
        if not eligible:
            st.info("This project has no completed jobs with structure-bearing design candidates." if project_design else "This project has no completed generator jobs with the required output type." if generated else "This project has no successful design jobs with filter-passing designs.")
            return contexts
        options = list(eligible)
        default = options.index(params["source_run_dir"]) if params.get("source_run_dir") in options else 0
        source_dir = st.selectbox("Completed design job" if project_design else "Generator job" if generated else "Successful design job", options, index=default, format_func=lambda p: f"{eligible[p][0].get('campaign_name') or eligible[p][0]['tool']} · {eligible[p][0].get('job_code')}", key=key+"-source")
        ids = [c["candidate_id"] for c in eligible[source_dir][1]]
        chosen = st.multiselect("Project design candidates" if project_design else "Generated designs" if generated else "Filter-passing designs", ids, default=[c for c in params.get("selected_candidate_ids", ids) if c in ids], key=key+"-designs-"+Path(source_dir).name)
        if project_design:st.caption("Choose saved structures with explicit target chains. Complex Results also includes non-binder predictions. Connect Candidate collection → Complex Analysis with Manual selection to analyze them without an acceptance decision.")
        elif generated:st.caption("These designs have not been validated. Connect sequence design and structure prediction steps to assess them.")
        contexts[node_id] = {"source_run_dir": source_dir, "selected_candidate_ids": chosen}
        if project_design:
            contexts[node_id]['source_scope'] = params['source_scope']
        params.update(contexts[node_id])
    return contexts


def _open_node_editor(prefix, node_id):
    graph = st.session_state[prefix+"-graph"]
    node = next((n for n in graph["nodes"] if n["id"] == node_id), None)
    if node:
        st.session_state[prefix+"-editor"] = {"id": uuid.uuid4().hex, "node_id": node_id, "params": copy.deepcopy(node["data"]["params"]), "revision": 0}


def _close_node_editor(prefix):
    st.session_state.pop(prefix+"-editor", None)


def _filter_observed_candidates(graph, node_id, run_dir):
    """Read the selected run's exact upstream node outputs for Filter choices."""
    try:
        root = resolve_stored_path(run_dir, must_exist=True)
        snapshot = read_json(root / "graph_snapshot.json")
        if snapshot != graph:
            return []
        execution = read_json(root / "graph_execution.json").get("nodes", {})
        from mn_molecular_design_studio import candidate_sets
        rows = []
        for edge in graph.get("edges", []):
            if edge.get("target") != node_id:
                continue
            source_state = execution.get(edge.get("source"), {})
            if source_state.get("status") not in {"completed", "cached", "partial"}:
                continue
            output = source_state.get("output") or {}
            output = (output.get("port_outputs") or {}).get(edge.get("sourceHandle", "out"), output)
            if not output.get("candidates_jsonl"):
                continue
            for candidate in candidate_sets.iter_candidates(output):
                rows.append(candidate)
                if len(rows) >= 5000:
                    return rows
        return rows
    except (OSError, ValueError, TypeError, KeyError):
        return []


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
        last_job = st.session_state.get(prefix + "-last-job") or {}
        observed_candidates = (_filter_observed_candidates(graph, working["id"], last_job.get("run_dir"))
                              if working["data"]["kind"] == "filter" and last_job.get("run_dir") else [])
        def render_input():
            return _render_node_input(working, owner_id, project_id, search_jobs, key, graph)
        values, label, action, replacement_kind, node_updates = render_node_form(
            working, graph, key=key, render_input=render_input, owner_id=owner_id,
            include_node_updates=True, observed_candidates=observed_candidates,
        )
        if action == "replace":
            try:
                from mn_molecular_design_studio.pipeline_node_ui import replacement_node_data
                updated = copy.deepcopy(graph)
                node = next(n for n in updated["nodes"] if n["id"] == selected["id"])
                replacement_source = copy.deepcopy(selected)
                replacement_source["data"].update(params=values, label=label)
                replacement = replacement_node_data(replacement_source, replacement_kind)
                node["data"].update(replacement)
                for name, value in node_updates.items():
                    if name == "execution" and value is None:
                        node["data"].pop(name, None)
                    else:
                        node["data"][name] = copy.deepcopy(value)
                updated = clean_graph(updated)
                validate_graph(updated)
                if project_id:
                    submit_api("/api/v1/pipelines/save", {"project_id":project_id, "graph":updated})
                st.session_state[prefix+"-graph"] = updated
                st.session_state[prefix+"-settings-update"] = {
                    "id":uuid.uuid4().hex, "node_id":node["id"],
                    "data": {**replacement, **node_updates},
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
                for name, value in node_updates.items():
                    if name == "execution" and value is None:
                        node["data"].pop(name, None)
                    else:
                        node["data"][name] = copy.deepcopy(value)
                updated = clean_graph(updated)
                if project_id:
                    submit_api("/api/v1/pipelines/save", {"project_id":project_id, "graph":updated})
                st.session_state[prefix+"-graph"] = updated
                st.session_state[prefix+"-settings-update"] = {
                    "id":uuid.uuid4().hex, "node_id":node["id"],
                    "data": {"params": values, "label": label, **node_updates},
                }
                _close_node_editor(prefix)
                st.rerun(scope="app")
            except (ValueError, OSError, TypeError, KeyError) as exc:
                st.error(str(exc))
    dialog()


def _activate_graph(prefix, graph):
    """Widget callback: replace the draft before its widgets are instantiated."""
    st.session_state[prefix+"-graph"] = copy.deepcopy(graph)
    st.session_state[prefix+"-reset"] = uuid.uuid4().hex
    _close_node_editor(prefix)


def _show_pipeline_save_as(prefix, project_id, submit_api):
    """Create a separately named pipeline draft from the current canvas."""
    if not st.session_state.get(prefix + "-save-as-open"):
        return

    @st.dialog("Save pipeline as…", width="small")
    def dialog():
        state_key = prefix + "-graph"
        source = copy.deepcopy(st.session_state[state_key])
        name = st.text_input("Pipeline name", value=f"{source['name']} copy", key=prefix + "-save-as-name-" + source["id"]).strip()
        save_col, cancel_col = st.columns(2)
        save = save_col.button("Save as…", key=prefix + "-save-as-submit", type="primary", use_container_width=True)
        cancel = cancel_col.button("Cancel", key=prefix + "-save-as-cancel", use_container_width=True)
        if cancel:
            st.session_state[prefix + "-save-as-open"] = False
            st.rerun()
        if save:
            if not project_id:
                st.error("Choose a project before saving a pipeline.")
                return
            if not name:
                st.error("Enter a pipeline name.")
                return
            duplicate = copy.deepcopy(source)
            duplicate["id"] = uuid.uuid4().hex
            duplicate["name"] = name
            try:
                submit_api("/api/v1/pipelines/save", {"project_id": project_id, "graph": duplicate})
                _activate_graph(prefix, duplicate)
                st.session_state[prefix + "-save-as-open"] = False
                st.success(f"Saved as {name}.")
                st.rerun()
            except (ValueError, OSError, TypeError, KeyError) as exc:
                st.error(str(exc))

    dialog()


def _activate_selected_template(prefix, templates):
    """Load the selected starter graph before rendering the canvas again."""
    selected = st.session_state.get(prefix + "-template")
    if selected in templates:
        _activate_graph(prefix, new_graph(templates[selected]))


def _show_pipeline_resources(prefix, owner_id, project_id, submit_api):
    """Canvas toolbar dialog for shared scheduler and workstation defaults."""
    if not st.session_state.get(prefix + "-resources-open"):
        return

    @st.dialog("Pipeline resources", width="large")
    def dialog():
        from mn_molecular_design_studio import candidate_execution
        from mn_molecular_design_studio.cluster import config as cluster_config

        state_key = prefix + "-graph"
        graph = copy.deepcopy(st.session_state[state_key])
        settings = candidate_execution.policy(graph.get("candidate_execution"))
        resource_key = prefix + "-" + graph["id"]
        nodes = []
        try:
            if cluster_config.load()["mode"] in {"peer", "coordinator"}:
                from mn_molecular_design_studio.cluster.store import Store
                nodes = Store().nodes()
        except (OSError, ValueError, KeyError, TypeError):
            nodes = []

        st.caption("These are pipeline-wide scheduler defaults. Open a canvas step's Configure dialog to override resources for that step.")
        sharded = st.checkbox(
            "Distribute compatible candidate batches",
            value=settings["mode"] == "sharded", key=resource_key + "-resources-sharded",
            help="Splits independent candidates into scheduler jobs. Engine-native batch size and design counts remain in each tool's settings.",
        )
        choices = {
            "local": "This workstation",
            "auto": "Automatic · compatible trusted workstations",
            **{item["id"]: cluster_config.display_name(item["name"]) + " · " + item["availability"]
               for item in nodes},
        }
        selected = settings.get("placement", "local")
        settings["placement"] = st.selectbox(
            "Workstation", list(choices), format_func=choices.get,
            index=list(choices).index(selected) if selected in choices else 0,
            key=resource_key + "-resources-placement",
        )
        left, middle, right = st.columns(3)
        settings["gpu_device"] = left.text_input(
            "GPU selector", value=settings.get("gpu_device", "all"), key=resource_key + "-resources-gpu",
            help="Use all, auto, or explicit device indices such as 0 or 0,1.",
        ).strip() or "all"
        settings["cpu_cores"] = int(middle.number_input(
            "CPU cores per job", min_value=1, max_value=256, value=settings["cpu_cores"], key=resource_key + "-resources-cpu",
        ))
        settings["min_vram_gb"] = float(right.number_input(
            "Minimum free VRAM (GiB)", min_value=0., value=float(settings["min_vram_gb"]), key=resource_key + "-resources-vram",
        ))

        resources_ready = True
        if sharded:
            settings["mode"] = "sharded"
            columns = st.columns(4)
            settings["batch_size"] = int(columns[0].number_input(
                "Candidates per scheduler batch", min_value=1, max_value=10000, value=settings["batch_size"], key=resource_key + "-resources-batch",
            ))
            settings["max_in_flight"] = int(columns[1].number_input(
                "Concurrent batches", min_value=1, max_value=32, value=settings["max_in_flight"], key=resource_key + "-resources-concurrency",
            ))
            settings["ram_gb"] = float(columns[2].number_input(
                "RAM per job (GiB)", min_value=0., value=float(settings["ram_gb"]), key=resource_key + "-resources-ram",
            ))
            settings["retries"] = int(columns[3].number_input(
                "Retries per failed batch", min_value=0, max_value=5, value=settings["retries"], key=resource_key + "-resources-retries",
            ))
            settings["incremental_cpu"] = st.checkbox(
                "Process completed batches while prediction continues (CPU)", value=settings["incremental_cpu"], key=resource_key + "-resources-incremental",
            )
            settings["failure_policy"] = "partial" if st.checkbox(
                "Allow explicit partial candidate pools and missing branches", value=settings["failure_policy"] == "partial",
                key=resource_key + "-resources-partial",
            ) else "all"
            from mn_molecular_design_studio.campaign_objectives_ui import render_time_limits
            st.caption("Dispatch time limits apply to each partitioned step. Shared campaign objectives can set limits across selected design branches.")
            limits, ready = render_time_limits(st, settings, key=resource_key + "-resources-time")
            resources_ready = ready
            if ready:
                for name in ("max_runtime_seconds", "deadline_at"):
                    settings.pop(name, None)
                settings.update(limits)
        else:
            settings["mode"] = "native"
            settings.pop("max_runtime_seconds", None)
            settings.pop("deadline_at", None)

        if nodes:
            with st.expander("Connected workstation availability", expanded=False):
                st.dataframe([{
                    "Workstation": cluster_config.display_name(item["name"]),
                    "Status": item["availability"],
                    "Free GPUs": len(item["snapshot"].get("resources", {}).get("free_gpus", [])),
                    "Running allocations": item["snapshot"].get("resources", {}).get("reservations", 0),
                    "Free CPU threads": item["snapshot"].get("resources", {}).get("cpu_free", 0),
                } for item in nodes], hide_index=True, use_container_width=True)

        apply_col, cancel_col = st.columns(2)
        if cancel_col.button("Cancel", key=resource_key + "-resources-cancel", use_container_width=True):
            st.session_state[prefix + "-resources-open"] = False
            st.rerun()
        if apply_col.button("Apply resource defaults", key=resource_key + "-resources-apply", type="primary", disabled=not resources_ready, use_container_width=True):
            try:
                graph["candidate_execution"] = candidate_execution.policy(settings)
                normalized = clean_graph(graph)
                if project_id:
                    submit_api("/api/v1/pipelines/save", {"project_id": project_id, "graph": normalized})
                st.session_state[state_key] = normalized
                st.session_state[prefix + "-resources-open"] = False
                if not project_id:
                    st.info("Select a project to save these resource defaults.")
                st.rerun()
            except (ValueError, OSError, TypeError, KeyError) as exc:
                st.error(str(exc))

    dialog()


def _render_shared_campaign_editor(prefix, graph):
    """Render the optional shared objective editor inside the canvas Campaign dialog."""
    ready = True
    stage_error = None
    time_ready = True
    from mn_molecular_design_studio.workflow_modes import capability as design_capability
    from mn_molecular_design_studio.engine.workflow_engine.workflows.generator_capabilities import generator_capabilities
    eligible_campaign_nodes = []
    passing_campaign_nodes = []
    for candidate_node in graph['nodes']:
        route = design_capability(candidate_node['data']['kind'])
        if not route.get('modular_source') or not route.get('route_registered') or not route.get('design_count_field'):
            continue
        eligible_campaign_nodes.append(candidate_node)
        try:
            campaign = generator_capabilities(route.get('tool'))
        except ValueError:
            campaign = None
        if campaign and campaign.campaign_adapter:
            passing_campaign_nodes.append(candidate_node)
    shared_existing = graph.get('shared_campaign') if isinstance(graph.get('shared_campaign'), dict) else None
    shared_enabled = st.checkbox(
        'Use a shared campaign objective across design branches',
        value=bool(shared_existing), key=prefix + '-shared-campaign-enabled-' + graph['id'],
        help='Count outputs or validated passing designs across compatible native branches. Generation remains in each native workflow; the existing controller, scheduler and ledger coordinate campaign-enabled routes.',
    )
    if shared_enabled:
        from mn_molecular_design_studio.distributed_campaign_config_ui import _default_validation_stages
        existing_shared = shared_existing or {}
        from mn_molecular_design_studio.campaign_objectives import OBJECTIVE_LABELS
        available_objectives = ['generated'] + (['accepted', 'time'] if passing_campaign_nodes else [])
        old_objective = existing_shared.get('objective_type', 'generated')
        default_objective = old_objective if old_objective in available_objectives else 'generated'
        objective_type = st.selectbox(
            'Campaign objective', available_objectives,
            index=available_objectives.index(default_objective),
            format_func=OBJECTIVE_LABELS.get,
            key=prefix + '-shared-objective-type-' + graph['id'],
            help='Generated counts every native output, including rejected outputs. Passing counts only supported campaign adapters that complete the configured validation stages.',
        )
        objective_nodes = passing_campaign_nodes if objective_type in {'accepted', 'time'} else eligible_campaign_nodes
        branch_ids = [n['id'] for n in objective_nodes]
        selected_branches = st.multiselect(
            'Design branches in the shared objective', branch_ids,
            default=[node for node in existing_shared.get('generator_nodes', []) if node in branch_ids],
            format_func=lambda value: next(n['data']['label'] + ' · ' + value for n in objective_nodes if n['id'] == value),
            key=prefix + '-shared-campaign-generators-' + graph['id'],
        )
        planned_by_node = {}
        for branch_id in selected_branches:
            branch = next(n for n in objective_nodes if n['id'] == branch_id)
            count_field = design_capability(branch['data']['kind']).get('design_count_field')
            planned_by_node[branch_id] = max(1, int(branch['data']['params'].get(count_field, 1)))
        default_generated_target = sum(
            int((existing_shared.get('generator_limits') or {}).get(branch, {}).get('target_count') or planned_by_node[branch])
            for branch in selected_branches
        )
        default_target = int(existing_shared.get('accepted_target') or 1) if objective_type == 'accepted' else int(existing_shared.get('generated_target') or default_generated_target or 1)
        controls = st.columns(4)
        objective_target = None if objective_type == 'time' else int(controls[0].number_input(
            'Unique passing designs' if objective_type == 'accepted' else 'Generated-design target', min_value=1, value=default_target,
            step=1, key=prefix + '-shared-target-' + graph['id'],
        ))
        planned_total = sum(planned_by_node.values())
        max_attempts = int(controls[1].number_input(
            'Maximum generation attempts', min_value=1, value=int(existing_shared.get('max_attempts') or max(objective_target or 1000, planned_total, 1)),
            step=1, key=prefix + '-shared-attempts-' + graph['id'],
            help='Upper bound on generation tries across the selected branches. A try may yield no output or a candidate that fails acceptance. Separate from each engine’s GPU inference batch size.',
        ))
        max_jobs = int(controls[2].number_input(
            'Overall concurrent jobs', min_value=1, value=int(existing_shared.get('max_concurrent_jobs') or 2),
            step=1, key=prefix + '-shared-jobs-' + graph['id'],
        ))
        max_generation_jobs = int(controls[3].number_input(
            'Concurrent generation jobs', min_value=1,
            value=int(existing_shared.get('max_concurrent_generation_jobs') or min(max_jobs, 2)),
            step=1, key=prefix + '-shared-generation-jobs-' + graph['id'],
        ))
        from mn_molecular_design_studio.campaign_objectives_ui import render_time_limits
        time_limits, time_ready = render_time_limits(
            st, existing_shared, key=prefix + '-shared-time-' + graph['id'], required=objective_type == 'time',
        )
        max_compute = 0.0
        if objective_type == 'accepted':
            validation_left, validation_right = st.columns(2)
            max_validations = int(validation_left.number_input(
                'Concurrent downstream validations', min_value=1,
                value=int(existing_shared.get('max_concurrent_validations') or min(max_jobs, 2)),
                step=1, key=prefix + '-shared-validation-jobs-' + graph['id'],
            ))
            max_compute = float(validation_right.number_input(
                'Maximum shared GPU seconds (0 = no limit)', min_value=0.0,
                value=float(existing_shared.get('max_compute_seconds') or 0),
                key=prefix + '-shared-compute-' + graph['id'],
            ))
        else:
            max_validations = min(max_jobs, int(existing_shared.get('max_concurrent_validations') or max_jobs))
        previous_placement = existing_shared.get('workstation_selection', 'auto')
        workstation_selection = st.selectbox(
            'Generator workstations', ['auto', 'local'],
            index=['auto', 'local'].index(previous_placement)
            if previous_placement in ('auto', 'local') else 0,
            format_func=lambda value: 'Automatic · compatible workstations' if value == 'auto' else 'This workstation',
            key=prefix + '-shared-workstation-selection-' + graph['id'],
            help='Automatic placement sends attempts to available workstations with compatible GPU runtime, images, and model files.',
        )
        limits = {key: copy.deepcopy(value) for key, value in (existing_shared.get('generator_limits') or {}).items() if key in selected_branches}
        with st.expander('Optional per-generator limits', expanded=False):
            for branch_id in selected_branches:
                old_limit = limits.get(branch_id, {})
                branch_columns = st.columns(3 if objective_type == 'generated' else 2)
                branch_target = 0
                if objective_type == 'generated':
                    branch_target = int(branch_columns[0].number_input(
                        f'{branch_id} generated-design target', min_value=0,
                        value=int(old_limit.get('target_count') or planned_by_node[branch_id]), step=1,
                        key=prefix + '-shared-limit-target-' + graph['id'] + '-' + branch_id,
                    ))
                left, right = branch_columns[-2:]
                attempt_limit = int(left.number_input(
                    f'{branch_id} maximum attempts (0 = shared budget)', min_value=0,
                    value=int(old_limit.get('max_attempts') or 0), step=1,
                    key=prefix + '-shared-limit-attempts-' + graph['id'] + '-' + branch_id,
                ))
                concurrency_limit = int(right.number_input(
                    f'{branch_id} concurrent generation jobs (0 = shared cap)', min_value=0,
                    value=int(old_limit.get('max_concurrent_generation_jobs') or 0), step=1,
                    key=prefix + '-shared-limit-jobs-' + graph['id'] + '-' + branch_id,
                ))
                if attempt_limit or concurrency_limit:
                    limits[branch_id] = {
                        **({'target_count': branch_target} if objective_type == 'generated' and branch_target else {}),
                        **({'max_attempts': attempt_limit} if attempt_limit else {}),
                        **({'max_concurrent_generation_jobs': concurrency_limit} if concurrency_limit else {}),
                    }
                else:
                    if objective_type == 'generated' and branch_target:
                        limits[branch_id] = {'target_count': branch_target}
                    else:
                        limits.pop(branch_id, None)
        stages = existing_shared.get('validation_stages') or _default_validation_stages()
        stage_error = None
        acceptance_node = None
        if objective_type == 'accepted':
            filter_nodes = [node for node in graph['nodes'] if node['data']['kind'] == 'shared_filter']
            choices = [''] + [node['id'] for node in filter_nodes]
            existing_node = existing_shared.get('acceptance_node') or ''
            acceptance_node = st.selectbox(
                'Optional filter node whose passing decisions define the goal', choices,
                index=choices.index(existing_node) if existing_node in choices else 0,
                format_func=lambda value: 'Existing campaign-stage criteria (legacy)' if not value else next(
                    node['data']['label'] + ' · ' + value for node in filter_nodes if node['id'] == value),
                key=prefix + '-acceptance-node-' + graph['id'],
                help='Select the existing format-aware Validation & Acceptance filter. Count/timed campaigns can omit it and apply filtering later.',
            )
            if acceptance_node:
                st.caption('Mandatory filters, binder format, target type, informational scores and ranking come from the selected node. Edit that node to change them.')

            stages_text = st.text_area(
                'Independent prediction settings (JSON)' if acceptance_node else 'Common downstream prediction and acceptance criteria (JSON)',
                value=json.dumps(stages, indent=2), height=240,
                key=prefix + '-shared-validation-stages-' + graph['id'],
            )
            try:
                stages = json.loads(stages_text)
            except json.JSONDecodeError as exc:
                stages = None
                stage_error = f'Common validation stages must be valid JSON (line {exc.lineno}, column {exc.colno}).'
        graph['shared_campaign'] = {
            'enabled': True, 'objective_type': objective_type, 'generator_nodes': selected_branches,
            **({'accepted_target': objective_target} if objective_type == 'accepted' else
               {'generated_target': objective_target} if objective_type == 'generated' else {}),
            'max_attempts': max_attempts,
            'max_concurrent_jobs': max_jobs,
            'max_concurrent_generation_jobs': max_generation_jobs,
            'max_concurrent_validations': max_validations,
            'workstation_selection': workstation_selection,
            'generator_limits': limits, 'validation_stages': stages if objective_type == 'accepted' else [],
            **({'acceptance_node': acceptance_node} if acceptance_node else {}),
            **time_limits,
            **({'max_compute_seconds': max_compute} if objective_type == 'accepted' and max_compute else {}),
        }
        if stage_error:
            st.error(stage_error)
        if objective_type == 'accepted':
            st.caption('Only registered native campaign adapters are eligible. The shared target counts unique sequences that pass the configured predictor and downstream metrics; native rejection, downstream rejection, duplicate and unresolved outcomes remain visible.')
        elif objective_type == 'time':
            st.caption('The shared timer keeps supported branches producing until the duration, deadline, or attempt budget is reached. Acceptance filtering is optional and can run later.')
        else:
            st.caption(f'Planned native design count: {planned_total}. Every emitted candidate is recorded, including native-filter rejections and rows without an available native filter. Shared generated target: {objective_target}.')
        if not selected_branches:
            st.info('Choose at least one compatible design branch to enable the shared campaign objective.')
    else:
        graph.pop('shared_campaign', None)
    return (not shared_enabled) or (bool(selected_branches) and time_ready and stage_error is None)


def _show_pipeline_campaign(prefix, project_id, submit_api):
    if not st.session_state.get(prefix + "-campaign-open"):
        return

    @st.dialog("Campaign objective", width="large")
    def dialog():
        state_key = prefix + "-graph"
        graph = copy.deepcopy(st.session_state[state_key])
        ready = _render_shared_campaign_editor(prefix, graph)
        apply_col, cancel_col = st.columns(2)
        if cancel_col.button("Cancel", key=prefix + "-campaign-cancel", use_container_width=True):
            st.session_state[prefix + "-campaign-open"] = False
            st.rerun()
        if apply_col.button("Apply campaign objective", key=prefix + "-campaign-apply", type="primary", disabled=not ready, use_container_width=True):
            try:
                normalized = clean_graph(graph)
                if project_id:
                    submit_api("/api/v1/pipelines/save", {"project_id": project_id, "graph": normalized})
                st.session_state[state_key] = normalized
                st.session_state[prefix + "-campaign-open"] = False
                st.rerun()
            except (ValueError, OSError, TypeError, KeyError) as exc:
                st.error(str(exc))

    dialog()


def render_pipeline_builder(*, owner_id, project_selector, submit_api, search_jobs):
    st.markdown("<style>.block-container{max-width:1900px;padding-top:1.35rem} [data-testid=stMainBlockContainer]{padding-left:2rem;padding-right:2rem}</style>", unsafe_allow_html=True)
    project_id = project_selector("visual-pipeline")
    prefix = f"visual-pipeline-{project_id or 'unassigned'}"
    state_key = prefix + "-graph"
    definitions = list_definitions(owner_id, project_id) if project_id else []
    if state_key not in st.session_state:
        st.session_state[state_key] = copy.deepcopy(definitions[0]["graph"] if definitions else new_graph("targeted_protein_binder"))
        st.session_state[prefix + "-reset"] = uuid.uuid4().hex
    reuse_key = prefix + "-reuse"
    st.session_state.setdefault(reuse_key, True)
    pending = st.session_state.get("pipeline-pending-draft")
    if pending and pending["project_id"] == project_id:
        _activate_graph(prefix, pending["graph"])
        st.session_state.pop("pipeline-pending-draft", None)
    current_graph = st.session_state[state_key]
    migrated_graph = migrate_template_graph(current_graph)
    if migrated_graph is not current_graph:
        _activate_graph(prefix, migrated_graph)

    graph = st.session_state[state_key]
    with st.expander("Saved pipelines & versions", expanded=True):
        if definitions:
            saved = st.selectbox("Saved draft", definitions, format_func=lambda d:f"{d['graph']['name']} · {d['versions']} saved versions", key=prefix+"-saved")
            st.button("Open saved draft", key=prefix+"-open", on_click=_activate_graph, args=(prefix, saved["graph"]))
        versions = list_versions(owner_id, project_id, graph["id"]) if project_id else []
        if versions:
            version = st.selectbox("Version history", versions, format_func=lambda d:f"{d['saved_at'][:19].replace('T',' ')} UTC · {d['revision'][:10]}", key=prefix+"-version-"+graph["id"])
            st.button("Restore version into draft", key=prefix+"-restore", on_click=_activate_graph, args=(prefix, version["graph"]))
        if not definitions and not versions:
            st.caption("Canvas edits autosave to the current draft. Use Save as… to create a separate pipeline; Save version keeps immutable history.")

    from mn_molecular_design_studio.scientific_tasks_ui import render_legacy_recipe_selector
    created = render_legacy_recipe_selector(prefix)
    if created is not None:
        _activate_graph(prefix, created)
    graph = st.session_state[state_key]
    reuse = bool(st.session_state.get(reuse_key, True))

    from mn_molecular_design_studio import candidate_execution
    from mn_molecular_design_studio.cluster import config as cluster_config
    settings = candidate_execution.policy(graph.get("candidate_execution"))
    graph["candidate_execution"] = settings
    nodes = []
    try:
        if cluster_config.load()["mode"] in {"peer", "coordinator"}:
            from mn_molecular_design_studio.cluster.store import Store
            nodes = Store().nodes()
    except (OSError, ValueError, KeyError, TypeError):
        nodes = []
    execution_controls_ready = True

    contexts = {}
    if not project_id:
        st.info("Choose a project before saving a pipeline version or running it. Each pipeline run is recorded under its selected project.")
    last_job = st.session_state.get(prefix+"-last-job", {})
    execution = {}
    if last_job.get("run_dir"):
        last_path = Path(last_job["run_dir"])
        execution = _matching_execution(graph, read_json(last_path/"graph_snapshot.json"), read_json(last_path/"graph_execution.json"))
    value = _component()(graph=graph, catalog=NODE_CATALOG, active_study_layout=ACTIVE_LEARNING_LAYOUT,
                         active_execution_ready=graph.get("template_id") == "unidock_active_learning",
                         reuse_completed=reuse,
                         reset_token=st.session_state[prefix+"-reset"], input_contexts=contexts,
                         execution=execution, settings_update=st.session_state.get(prefix+"-settings-update"),
                         key=prefix+"-canvas", default=None)
    if value:
        signature = json.dumps(value, sort_keys=True)
        if signature != st.session_state.get(prefix+"-last-value"):
            st.session_state[prefix+"-last-value"] = signature
            try:
                component_reuse = value.get("reuse_completed", reuse)
                if type(component_reuse) is bool:
                    st.session_state[reuse_key] = component_reuse
                reuse = bool(st.session_state.get(reuse_key, True))
                component_graph = copy.deepcopy(value["graph"])
                component_graph["candidate_execution"] = copy.deepcopy(graph["candidate_execution"])
                for component_node in component_graph.get("nodes", []):
                    previous_node = next((n for n in graph["nodes"] if n["id"] == component_node["id"]), {})
                    if previous_node.get("data", {}).get("execution"):
                        component_node["data"]["execution"] = copy.deepcopy(previous_node["data"]["execution"])
                    for branch_option in ('enabled', 'retained_run_dir', 'evidence_namespace'):
                        if branch_option in previous_node.get('data', {}):
                            component_node['data'][branch_option] = copy.deepcopy(previous_node['data'][branch_option])
                for field_name in ("template_id", "workflow_mode", "task_configuration", "shared_campaign"):
                    if field_name in graph:
                        component_graph[field_name] = copy.deepcopy(graph[field_name])
                    else:
                        component_graph.pop(field_name, None)
                updated = _restore_component_parameters(component_graph, graph)
                updated["name"] = component_graph.get("name") or graph["name"]
                if project_id:
                    submit_api("/api/v1/pipelines/save", {"project_id": project_id, "graph": updated})
                st.session_state[state_key] = graph = updated
                action=value.get("action") or {}
                if action.get("id") and action["id"] != st.session_state.get(prefix+"-handled-action"):
                    st.session_state[prefix+"-handled-action"] = action["id"]
                    if action["type"] == "edit_node":
                        _open_node_editor(prefix, action.get("node_id"))
                    elif action["type"] == "campaign":
                        st.session_state[prefix + "-campaign-open"] = True
                    elif action["type"] == "resources":
                        st.session_state[prefix + "-resources-open"] = True
                    elif action["type"] == "save":
                        if not project_id:
                            st.warning("Choose a project before saving a pipeline.")
                        else:
                            submit_api("/api/v1/pipelines/save", {"project_id": project_id, "graph": graph})
                            st.success(f"Saved {graph['name']}.")
                    elif action["type"] == "save_as":
                        if not project_id:
                            st.warning("Choose a project before saving a pipeline.")
                        else:
                            st.session_state[prefix + "-save-as-open"] = True
                    elif not project_id:
                        st.warning("Select a project before saving a version or running the pipeline.")
                    elif action["type"] == "save_version":
                        result=submit_api("/api/v1/pipelines/save", {"project_id":project_id,"graph":graph,"version":True})
                        st.success(f"Version {result['revision'][:10]} saved. Existing versions and runs are preserved.")
                    elif action["type"] == "run":
                        if not execution_controls_ready:
                            raise ValueError('Correct the dispatch time limits before starting this pipeline.')
                        result=submit_api("/api/v1/pipelines/run", {"project_id":project_id,"graph":graph,"reuse_completed":reuse,"gpu_device":"all"})
                        st.session_state[prefix+"-last-job"] = result["job"]
                        st.session_state["last_submitted_job"] = result["job"]
                        st.success(f"Pipeline queued · {result['job'].get('job_code') or result['job']['run_id']}")
            except (ValueError, OSError, TypeError, KeyError) as exc:
                st.error(str(exc))
    _show_pipeline_campaign(prefix, project_id, submit_api)
    _show_pipeline_resources(prefix, owner_id, project_id, submit_api)
    _show_pipeline_save_as(prefix, project_id, submit_api)
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
            from mn_molecular_design_studio.workflow_modes import gpu_monitor_rows
            from mn_molecular_design_studio.cluster.capabilities import resources as scheduler_resources
            online, local_resources = scheduler_resources()
            inventory = list(nodes)
            if not any(n['id'] == cluster_config.identity() for n in inventory):
                inventory.append(dict(id=cluster_config.identity(), name='This workstation', availability='online' if online else 'unknown', snapshot={'resources': local_resources}))
            gpu_rows = gpu_monitor_rows(inventory)
            if gpu_rows:
                st.caption('GPU availability from scheduler/heartbeat snapshots. Model, image, checkpoint and memory eligibility is checked separately for every job. Queued work has no GPU allocation until leased.')
                st.dataframe(gpu_rows, hide_index=True)
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
                if execution.get("comparison_status") == "partial":
                    st.warning("Comparison partially completed: successful branches are available in Results. Failed branches and their errors are listed below.")
                shared_progress = execution.get('shared_campaign')
                if isinstance(shared_progress, dict):
                    objective_type = shared_progress.get('objective_type', 'accepted')
                    target = max(1, int(shared_progress.get('objective_target') or shared_progress.get('accepted_target') or shared_progress.get('generated_target') or 1))
                    count = int(shared_progress.get('generated_designs') or 0) if objective_type == 'generated' else int(shared_progress.get('unique_accepted') or 0)
                    count_label = 'generated designs' if objective_type == 'generated' else 'unique passing designs'
                    if objective_type == 'time':
                        remaining = shared_progress.get('remaining_seconds')
                        st.caption(f"Timed campaign · {shared_progress.get('generated_designs', 0)} generated · "
                                   + (f"{remaining / 60:.1f} minutes remaining" if remaining is not None else str(shared_progress.get('status'))))
                    else:
                        st.progress(min(1.0, count / target), text=f"Shared objective · {count}/{target} {count_label} · {shared_progress.get('status', 'running')}")
                    if objective_type == 'generated':
                        scheduler_text = (f"Scheduler generation jobs {shared_progress.get('scheduler_jobs_active', 0)}/"
                                          f"{shared_progress.get('max_concurrent_generation_jobs', 0)} active · "
                                          f"scientifically passing {shared_progress.get('scientifically_passing_designs', 0)} · "
                                          f"scientific rejects {shared_progress.get('scientifically_rejected_designs', 0)}")
                    else:
                        scheduler_text = (f"Attempts {shared_progress.get('attempts_authorized', 0)}/{shared_progress.get('max_attempts', 0)} · "
                                          f"active {shared_progress.get('attempts_active', 0)} generation and {shared_progress.get('validations_active', 0)} validation jobs · "
                                          f"scientific rejections {shared_progress.get('unique_rejected', 0)} · "
                                          f"duplicate sequences {shared_progress.get('duplicate_candidates', 0)}")
                    compute_seconds = shared_progress.get('compute_seconds')
                    compute_text = f"GPU compute {float(compute_seconds):.1f}s" if compute_seconds is not None else 'GPU compute unavailable for this native route'
                    st.caption(scheduler_text + ' · '
                        f"native rejections {shared_progress.get('native_rejected_unique', 0) + shared_progress.get('generated_native_rejected', 0)} · "
                        f"native pass {shared_progress.get('generated_native_pass', 0)} · "
                        f"native filter unavailable {shared_progress.get('generated_native_filter_unavailable', 0)} · "
                        + compute_text)
                    if shared_progress.get('per_generator'):
                        st.dataframe([
                            {'Generator branch': branch, **counts}
                            for branch, counts in shared_progress['per_generator'].items()
                        ], hide_index=True, use_container_width=True)
                batch_rows = []
                for node_id, step in execution.get('nodes', {}).items():
                    for batch in step.get('batches', []):
                        attempt = (batch.get('attempts') or [{}])[-1]
                        batch_rows.append({'Step': node_id, 'Batch': batch['batch_id'], 'Status': batch['status'],
                            'Worker': attempt.get('worker', batch.get('worker', '')),
                            'GPU': str((attempt.get('allocation') or {}).get('gpu_device', '')),
                            'Completed inputs': batch.get('completed_units', 0), 'Total inputs': batch['unit_count'],
                            'Retries': max(0, len(batch.get('attempts', []))-1), 'Detail': batch.get('error') or batch.get('wait_reason', '')})
                queued_steps = sum(step.get('status') in {'queued', 'pending'} for step in execution.get('nodes', {}).values())
                queued_chunks = sum(batch['Status'] in {'queued', 'pending', 'retrying', 'submitting'} for batch in batch_rows)
                st.caption(f'Queued/pending stages: {queued_steps} · Queued/pending processing jobs: {queued_chunks}')
                if batch_rows:
                    completed = sum(b['Completed inputs'] for b in batch_rows)
                    total = sum(b['Total inputs'] for b in batch_rows)
                    st.progress(completed / total if total else 0, text=f'{completed}/{total} candidate-stage inputs completed')
                    st.dataframe(batch_rows, hide_index=True, use_container_width=True)
                if execution.get("nodes"):
                    st.dataframe([{"Step": next((n["data"]["label"] for n in snapshot.get("nodes", []) if n["id"]==nid),nid), "Status": state.get("status"), "Job":state.get("job_code", ""), "Candidates":state.get("candidate_count"),
                        **{key: state.get('native_execution', {}).get(key) for key in ('collection_chunks', 'native_invocations', 'docker_processes', 'gpu_allocation_events', 'gpu_slots', 'prediction_outputs')},
                        "Detail":state.get("error", "")} for nid,state in execution["nodes"].items()], hide_index=True, use_container_width=True)
                    stage_rows = [{'Step': nid, 'Stage count': name, 'Count': count} for nid, step in execution['nodes'].items()
                                  for name, count in step.get('stage_counts', {}).items()]
                    if stage_rows:
                        st.dataframe(stage_rows, hide_index=True, use_container_width=True)
