"""Streamlit host for the freely movable pipeline canvas."""
from __future__ import annotations

import copy
import hashlib
import json
import uuid
from collections import Counter
from pathlib import Path
from urllib.parse import urlencode

import streamlit as st
import streamlit.components.v1 as components

from mn_molecular_design_studio.pipeline_graph import (
    ACTIVE_LEARNING_LAYOUT, NODE_CATALOG, PipelineInputError, clean_graph, fingerprint, is_screening_graph,
    get_shared_pipeline_template, list_definitions, list_pipeline_templates, list_shared_pipeline_templates,
    list_versions, migrate_template_graph, new_graph, screening_graph_issues,
    shared_template_starting_points, validate_graph,
)
from mn_molecular_design_studio.api_client import StudioAPIConflictError
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


def _render_active_library(params, key, *, library_catalog_available=True):
    if not library_catalog_available:
        st.info("The live library catalog is unavailable in this development preview. Open the pipeline in Studio to choose prepared fingerprints.")
        return None
    libraries = [row for row in _active_library_options() if row[1]]
    if not libraries:
        st.info("No registered fingerprint manifests are accessible on this workstation. Check the library registration and preparation files in Library tools.")
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


def _render_active_centroids(params, key, *, library_catalog_available=True):
    if not library_catalog_available:
        st.info("The live library catalog is unavailable in this development preview. Open the pipeline in Studio to choose an existing centroid set.")
        return None
    libraries = [row for row in _active_library_options() if row[2]]
    if not libraries:
        st.info("No registered centroid manifests are accessible on this workstation. Check the library registration and preparation files in Library tools.")
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


def _render_node_input(selected, owner_id, project_id, search_jobs, prefix, graph=None, *, library_catalog_available=True):
    if selected["data"]["kind"] == "ligand_sdf_input":
        from mn_molecular_design_studio.ligand_tools_ui import _save_reference_ligand_upload

        params = selected["data"]["params"]
        path = resolve_stored_path(params.get("ligand_sdf"), must_exist=True)
        upload = st.file_uploader("Ligand structure", type=["sdf"], key=prefix + "-ligand-sdf")
        if upload is not None:
            try:
                path = _save_reference_ligand_upload(owner_id, upload.name, upload.getvalue())
            except (OSError, ValueError) as exc:
                st.error(str(exc))
                return {}
        if path:
            st.caption(f"Selected {upload.name if upload else path.name}")
            params["ligand_sdf"] = str(path)
            return {selected["id"]: {"ligand_sdf": str(path)}}
        return {}
    if selected["data"]["kind"] == "md_source":
        from mn_molecular_design_studio.md_analysis_ui import _cached_source_options_for_job

        if not project_id:
            st.info("Choose a project before selecting an MD run.")
            return {}
        rows, cursor = [], None
        while True:
            page = search_jobs(project_id=project_id, task_group="ligand-md", statuses=["completed"], limit=1000, cursor=cursor)
            rows.extend(row for row in page["items"] if row.get("run_dir") and Path(row["run_dir"]).is_dir())
            cursor = page.get("next_cursor")
            if not cursor:
                break
        if not rows:
            st.info("No completed MD runs are present on this pipeline workstation. Open the pipeline on the workstation holding the trajectory; it is not copied between computers.")
            return {}
        params = selected["data"]["params"]
        by_path = {str(row["run_dir"]): row for row in rows}
        paths = list(by_path)
        source = st.selectbox("Completed MD job", paths,
            index=paths.index(params["source_run_dir"]) if params.get("source_run_dir") in paths else 0,
            format_func=lambda path: str(by_path[path].get("campaign_name") or by_path[path].get("job_code") or by_path[path].get("run_id")),
            key=prefix + "-md-job")
        options = _cached_source_options_for_job(by_path[source])
        lookup = {row["option_id"]: row for row in options}
        if not lookup:
            st.info("No exact trajectory/topology pair is recorded for this MD job.")
            return {}
        identifiers = list(lookup)
        chosen = st.selectbox("Recorded production segment", identifiers,
            index=identifiers.index(params["option_id"]) if params.get("option_id") in identifiers else 0,
            format_func=lambda value: lookup[value]["label"], key=prefix + "-md-segment-" + Path(source).name)
        topologies = [path.relative_to(Path(source).resolve()).as_posix() for path in lookup[chosen]["topology_paths"]]
        topology = st.selectbox("Matching topology", topologies,
            index=topologies.index(params["topology_relative_path"]) if params.get("topology_relative_path") in topologies else 0,
            key=prefix + "-md-topology-" + hashlib.sha256(chosen.encode()).hexdigest()[:12])
        st.caption("One segment and replica only. Source files stay read-only on this workstation. Configure ligand selection and sampling in MD Analysis.")
        values = {"source_run_dir": source, "option_id": chosen, "topology_relative_path": topology}
        params.update(values)
        return {selected["id"]: values}
    if selected["data"]["kind"] == "run_source":
        from mn_molecular_design_studio.run_inputs_ui import select
        choice = select(project_id, prefix + "-run-output", require_complex=False, allow_docking=True,
                        reference=selected["data"]["params"].get("run_reference"),
                        selected_docking_poses=selected["data"]["params"].get("selected_docking_poses"))
        if choice:
            selected["data"]["params"].update(
                run_reference=choice["reference"],
                selected_docking_poses=choice.get("selected_docking_poses"),
            )
            return {selected["id"]: {
                "run_reference": choice["reference"],
                "selected_docking_poses": choice.get("selected_docking_poses"),
            }}
        return {}
    from mn_molecular_design_studio.target_editor import (
        TargetSource, _hash, _target_widget_keys, render_saved_target_sequence, render_target_coldspots,
        render_target_editor, render_target_source,
    )
    from mn_molecular_design_studio.design_validation import design_candidates
    from mn_molecular_design_studio.generation_io import generated_candidates
    contexts = {}
    node_id, params = selected["id"], selected["data"]["params"]
    key = f"{prefix}-{node_id}"
    if selected["data"]["kind"] == "target" and params.get("target_input") == "idr":
        from mn_molecular_design_studio.binder_pipeline import stage_idr_target

        sequence = st.text_area("Target sequence or FASTA", value=params.get("target_sequence", ""),
                                key=key + "-idr-sequence")
        hotspots = st.text_input("Hotspot residues (optional)", value=params.get("hotspots", ""),
                                key=key + "-idr-hotspots", help="Use 1-based residue positions in the sequence.")
        if sequence.strip():
            try:
                values = stage_idr_target(sequence, owner_id, hotspots=hotspots)
            except (ValueError, OSError) as exc:
                st.error(str(exc))
                values = {"target_pdb": "", "target_sequence": ""}
        else:
            values = {"target_pdb": "", "target_sequence": ""}
        params.update(values)
        return {node_id: values}
    if selected["data"]["kind"] == "target":
        from mn_molecular_design_studio.rfd4_proteina_input import CYCLIC_TASK

        successors = {edge["target"] for edge in (graph or {}).get("edges", []) if edge["source"] == node_id}
        if any(node["id"] in successors and node["data"]["kind"] == "rfd4"
               and node["data"]["params"].get("task") == CYCLIC_TASK for node in (graph or {}).get("nodes", [])):
            from mn_molecular_design_studio import rfd4_proteina_input as inputs

            st.caption("Upload the annotated mmCIF with the receptor, peptide ring and design conditions.")
            upload = st.file_uploader("Annotated cyclic complex", type=["cif", "mmcif"], key=key + "-cyclic-complex")
            try:
                path = resolve_stored_path(params.get("target_pdb"), must_exist=True) if params.get("target_pdb") else None
                if upload is not None:
                    path, = inputs.store_structure_uploads([upload], owner_id)
                if path is None:
                    return {}
                ring, receptors = inputs.cyclic_complex_chains(path)
            except (ValueError, OSError) as exc:
                st.error(str(exc))
                return {}
            st.caption(f"{path.name} · receptor chains {', '.join(receptors)} · cyclic chain {ring}")
            values = {"target_pdb": str(path), "target_chains": ",".join(receptors),
                      "hotspots": "", "coldspots": "", "target_source_molecule_id": "", "target_source_structure_id": ""}
            params.update(values)
            return {node_id: values}
        path = resolve_stored_path(params.get("target_pdb"), must_exist=True) if params.get("target_pdb") else None
        source = render_target_source(
            key_prefix=key, label="Target Structure", example_code="7QF4",
            allow_cif=True, upload_in_expander=True, owner_id=owner_id,
            allow_saved_targets=True, project_id=project_id,
            initial_molecule_id=params.get("target_source_molecule_id", "") if not path else "",
            initial_structure_id=params.get("target_source_structure_id", ""),
        )
        if not source.path and path:
            source = TargetSource(path=path, name=path.stem, owner_id=owner_id,
                                  molecule_id=params.get("target_source_molecule_id", ""),
                                  structure_id=params.get("target_source_structure_id", ""))
        render_saved_target_sequence(owner_id=owner_id, project_id=project_id, key=key)
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
        chosen = _render_active_library(params, key, library_catalog_available=library_catalog_available)
        if chosen:
            contexts[node_id] = chosen
            params.update(chosen)
    elif selected["data"]["kind"] == "active_seed":
        chosen = _render_active_centroids(params, key, library_catalog_available=library_catalog_available)
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


def _show_node_editor(prefix, *, owner_id, project_id, submit_api, search_jobs, stage_edits=False,
                      library_catalog_available=True):
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
            return _render_node_input(
                working, owner_id, project_id, search_jobs, key, graph,
                library_catalog_available=library_catalog_available,
            )
        values, label, action, replacement_kind, node_updates = render_node_form(
            working, graph, key=key, render_input=render_input, owner_id=owner_id,
            project_id=project_id, include_node_updates=True, observed_candidates=observed_candidates,
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
                validate_graph(updated, allow_incomplete=stage_edits)
                if project_id and not stage_edits:
                    _save_draft(prefix, project_id, updated, submit_api)
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
                if stage_edits:
                    validate_graph(updated, allow_incomplete=True)
                if project_id and not stage_edits:
                    _save_draft(prefix, project_id, updated, submit_api)
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


def _draft_revision(prefix, graph):
    identifier, revision = st.session_state.get(prefix + "-draft-base", (None, None))
    return revision if identifier == graph["id"] else None


def _save_draft(prefix, project_id, graph, submit_api, *, version=False, from_canvas=False):
    try:
        saved = submit_api("/api/v1/pipelines/save", {"project_id": project_id, "graph": graph,
            "version": version, "expected_revision": _draft_revision(prefix, graph)})
    except StudioAPIConflictError:
        st.session_state[prefix + "-graph"] = graph
        st.session_state[prefix + "-draft-conflict"] = graph["id"]
        if not from_canvas:
            st.session_state[prefix + "-reset"] = uuid.uuid4().hex
        raise
    if (saved or {}).get("revision"):
        st.session_state[prefix + "-draft-base"] = (graph["id"], saved["revision"])
    st.session_state.pop(prefix + "-draft-conflict", None)
    return saved


def _activate_graph(prefix, graph, revision=None):
    """Widget callback: replace the draft before its widgets are instantiated."""
    st.session_state[prefix+"-graph"] = copy.deepcopy(graph)
    st.session_state[prefix + "-draft-base"] = (graph["id"], revision)
    st.session_state[prefix + "-view"] = "editor"
    st.session_state.pop(prefix + "-restore-viewport", None)
    st.session_state.pop(prefix + "-draft-conflict", None)
    st.session_state.pop(prefix + "-run-check", None)
    st.session_state[prefix+"-reset"] = uuid.uuid4().hex
    st.session_state.pop(prefix + "-guided-open", None)
    st.session_state.pop(prefix + "-binder-open", None)
    _close_node_editor(prefix)


def open_binder_preset(preset, *, prefix=None):
    from mn_molecular_design_studio.binder_pipeline import new_binder_pipeline

    _open_capability_graph(new_binder_pipeline(**preset), prefix)


def _open_capability_graph(graph, prefix):
    if prefix is None:
        st.session_state.pop(_TEMPLATE_VIEW_PREFIX + "-open", None)
        st.session_state["pipeline-pending-draft"] = {"project_id": None, "graph": graph}
        st.session_state.setdefault("project-pending-selection", {})["visual-pipeline-project-id"] = ""
        st.session_state["page"] = "Pipelines"
    else:
        _activate_graph(prefix, graph)
    st.rerun()


def _binder_picker(key, *, connected_target=None):
    from mn_molecular_design_studio.binder_pipeline import (
        FORMATS, MULTIDOMAIN_TARGET, STRUCTURE_TARGET_INPUTS, TARGET_INPUTS, compatible_routes, legacy_sequence_carrier,
    )
    from mn_molecular_design_studio.design_capabilities import ACCESS, picker_tools
    from mn_molecular_design_studio.design_capabilities_ui import render_binder_preset

    if connected_target:
        params = connected_target["data"]["params"]
        if legacy_sequence_carrier(params):
            st.info("This legacy Target carries a sequence. Use New pipeline → Design a binder → Sequence to configure an IDR route.")
            return None
        target_input = params.get("target_input", "structure")
        st.caption(f"Target step · {connected_target['data']['label']}")
    else:
        sources = {TARGET_INPUTS[value]: value for value in STRUCTURE_TARGET_INPUTS}
        source = st.segmented_control("Target input", [*sources, "Sequence"], default="Structure", key=key + "-source")
        if source is None:
            return None
        target_input = sources.get(source)
        if source == "Sequence":
            target_input = st.selectbox("Target context", ["idr", "sequence"], index=None,
                                        format_func=TARGET_INPUTS.get, placeholder="Choose how to interpret the sequence",
                                        key=key + "-context")
            if target_input is None:
                return None
    if target_input == MULTIDOMAIN_TARGET:
        st.caption("These workflows use the supplied arrangement of the target domains. Select the relevant chains in the Target step.")
    binder_format = st.selectbox("Binder format", list(FORMATS), format_func=FORMATS.get, key=key + "-format")
    routes = compatible_routes(target_input, binder_format)
    additional = picker_tools(target_input, binder_format)
    engines = {}
    for route in routes:
        engine_name = route.engine if route.format == binder_format else f"{route.engine} ({FORMATS[route.format]})"
        engines[route.engine] = f"{engine_name} · {ACCESS['pipeline']}"
    for entry in additional:
        engines.setdefault(entry.engine, f"{entry.engine} · {ACCESS[entry.access]}")
    if not engines:
        st.info("No documented binder-design engine for this target input and format.")
        return None
    engine = st.selectbox("Design engine", list(engines), format_func=engines.get,
                          key=f"{key}-{target_input}-{binder_format}-engine")
    variants = [route for route in routes if route.engine == engine]
    if not variants:
        entry = next(entry for entry in additional if entry.engine == engine)
        st.info(entry.note)
        if entry.reference:
            st.markdown(f"[Project reference]({entry.reference})")
        return None
    return render_binder_preset(f"{key}-{target_input}-{binder_format}-{engine}", variants,
                                button_label="Add to canvas" if connected_target else "Create binder pipeline")


def _show_binder_branch(prefix, graph, project_id, submit_api):
    from mn_molecular_design_studio.binder_pipeline import new_binder_pipeline

    open_key = prefix + "-binder-open"
    target = next((node for node in graph["nodes"] if node["id"] == st.session_state.get(open_key)
                   and node["data"]["kind"] == "target"), None)
    if target is None:
        return

    @st.dialog("Design a binder", on_dismiss=lambda: st.session_state.pop(open_key, None))
    def dialog():
        preset = _binder_picker(f"{prefix}-{graph['id']}-{target['id']}-binder", connected_target=target)
        if preset is None:
            return
        try:
            created = new_binder_pipeline(**preset, graph=graph, target_node_id=target["id"])
            revision = _draft_revision(prefix, graph)
            if project_id and revision:
                saved = _save_draft(prefix, project_id, created, submit_api)
                created, revision = saved["graph"], saved["revision"]
        except (ValueError, RuntimeError, OSError, KeyError) as exc:
            st.error(str(exc))
            return
        _activate_graph(prefix, created, revision)
        st.rerun()

    dialog()


def _show_guided_pipeline_setup(prefix, graph, project_id, submit_api):
    open_key = prefix + "-guided-open"
    if st.session_state.get(open_key) != graph["id"] or graph["nodes"]:
        return

    from mn_molecular_design_studio.pipeline_graph import MAIN_STARTERS
    from mn_molecular_design_studio.scientific_tasks_ui import render_task_setup

    @st.dialog("Start a pipeline", width="medium", on_dismiss=lambda: st.session_state.pop(open_key, None))
    def dialog():
        if st.button("Start with an empty canvas", key=open_key + "-blank"):
            st.session_state.pop(open_key, None)
            st.rerun()
        categories = {value: name for name, value in MAIN_STARTERS.items()}
        categories["targeted_protein_binder"] = "Design a binder"
        categories["advanced"] = "Advanced workflows"
        key = prefix + "-guided-" + graph["id"]
        task = st.selectbox("Workflow type", list(categories), format_func=categories.get, index=None, key=key + "-task")
        if task is None:
            return
        if task == "targeted_protein_binder":
            from mn_molecular_design_studio.binder_pipeline import new_binder_pipeline

            preset = _binder_picker(key + "-binder")
            if preset is None:
                return
            try:
                created = new_binder_pipeline(**preset)
            except ValueError as exc:
                st.error(str(exc))
                return
        else:
            created = render_task_setup(key + "-" + task, task=task, compact=True)
        if created is None:
            return
        created["id"] = graph["id"]
        if graph["name"] != "Untitled pipeline":
            created["name"] = graph["name"]
        revision = _draft_revision(prefix, graph)
        if project_id and revision:
            try:
                saved = _save_draft(prefix, project_id, created, submit_api)
                created, revision = saved["graph"], saved["revision"]
            except (ValueError, RuntimeError, OSError, KeyError) as exc:
                st.error(str(exc))
                return
        _activate_graph(prefix, created, revision)
        st.rerun()

    dialog()


def _show_pipeline_save_project(prefix, graph, project_selector, submit_api, read_only_projects):
    open_key = prefix + "-save-project-open"
    if not st.session_state.get(open_key):
        return

    @st.dialog("Save pipeline to a project", on_dismiss=lambda: st.session_state.pop(open_key, None))
    def dialog():
        project_id = project_selector("visual-pipeline-save-destination")
        save, cancel = st.columns(2)
        writable = bool(project_id and project_id not in read_only_projects)
        if save.button("Save pipeline", type="primary", key=open_key + "-confirm", disabled=not writable):
            try:
                saved = submit_api("/api/v1/pipelines/save", {"project_id": project_id, "graph": graph})
                _activate_graph(f"visual-pipeline-{project_id}", saved["graph"], saved["revision"])
                st.session_state.setdefault("project-pending-selection", {})["visual-pipeline-project-id"] = project_id
                st.session_state.pop(prefix + "-graph", None)
                st.session_state.pop(open_key, None)
                st.rerun()
            except (ValueError, OSError, RuntimeError, KeyError, TypeError) as exc:
                st.error(str(exc))
        if cancel.button("Cancel", key=open_key + "-cancel"):
            st.session_state.pop(open_key, None)
            st.rerun()

    dialog()


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
                saved = submit_api("/api/v1/pipelines/save", {"project_id": project_id, "graph": duplicate})
                _activate_graph(prefix, duplicate, (saved or {}).get("revision"))
                st.session_state[prefix + "-save-as-open"] = False
                st.success(f"Saved as {name}.")
                st.rerun()
            except (ValueError, OSError, TypeError, KeyError) as exc:
                st.error(str(exc))

    dialog()


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
                    _save_draft(prefix, project_id, normalized, submit_api)
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
                    _save_draft(prefix, project_id, normalized, submit_api)
                st.session_state[state_key] = normalized
                st.session_state[prefix + "-campaign-open"] = False
                st.rerun()
            except (ValueError, OSError, TypeError, KeyError) as exc:
                st.error(str(exc))

    dialog()


def _change_archive_status(draft, project_id, submit_api, *, archived):
    try:
        submit_api(f"/api/v1/pipelines/{draft['graph']['id']}/archive", {
            "project_id": project_id, "archived": archived, "expected_revision": draft["revision"]})
        return True
    except (ValueError, RuntimeError) as exc:
        st.error(str(exc))
        return False


def _open_saved_pipeline(prefix, draft):
    _activate_graph(prefix, draft["graph"], draft["revision"])


_TEMPLATE_VIEW_PREFIX = "visual-pipeline-shared-template"


def _open_template_view(template_id, return_prefix, project_id, *, working_routes=None):
    working_routes = {} if working_routes is None else working_routes
    template = working_routes.get(template_id) or get_shared_pipeline_template(template_id)
    _activate_graph(_TEMPLATE_VIEW_PREFIX, template["graph"])
    st.session_state[_TEMPLATE_VIEW_PREFIX + "-restore-viewport"] = bool(
        template.get("restore_viewport", template.get("saved_at"))
    )
    st.session_state[_TEMPLATE_VIEW_PREFIX + "-open"] = {
        "id": template_id, "revision": template["revision"], "return_prefix": return_prefix,
        "return_project": project_id, "working_routes": working_routes,
    }


def _show_save_pipeline_template(prefix, submit_api):
    key = prefix + "-template-source"
    source = st.session_state.get(key)
    if not source:
        return

    @st.dialog("Save as shared template", on_dismiss=lambda: st.session_state.pop(key, None))
    def dialog():
        name = st.text_input("Template name", value=source["name"], max_chars=160,
                             key=key + "-name-" + source["id"]).strip()
        st.caption("All signed-in users can use this template. Project inputs are chosen when setting up a pipeline.")
        if st.button("Save template", type="primary", disabled=not name, key=prefix + "-publish-template"):
            try:
                submit_api("/api/v1/pipelines/shared-templates/save", {"graph": source, "name": name})
            except (ValueError, RuntimeError, OSError, LookupError) as exc:
                st.error(str(exc))
                return
            st.session_state.pop(key, None)
            st.session_state[prefix + "-view"] = "library"
            st.session_state[prefix + "-library-next-section"] = "Templates"
            st.session_state[prefix + "-library-notice"] = f"Saved shared template {name}."
            st.rerun()

    dialog()


def _change_pipeline_trash(prefix, pipeline, project_id, submit_api, *, trashed):
    graph = pipeline["graph"]
    try:
        submit_api(f"/api/v1/pipelines/{graph['id']}/trash", {
            "project_id": project_id, "trashed": trashed, "expected_revision": pipeline["revision"],
        })
    except (ValueError, RuntimeError, OSError, LookupError) as exc:
        st.error(str(exc))
        return
    if trashed and (st.session_state.get(prefix + "-graph") or {}).get("id") == graph["id"]:
        for suffix in ("-graph", "-draft-base", "-draft-conflict"):
            st.session_state.pop(prefix + suffix, None)
    message = f"Moved {graph['name']} to Trash." if trashed else f"Restored {graph['name']}."
    st.session_state[prefix + "-library-notice"] = message
    st.session_state[prefix + "-library-next-section"] = "Saved pipelines"
    st.rerun()


def _render_pipeline_library(prefix, owner_id, project_id, definitions, archived, trashed,
                             submit_api, search_jobs, search_cluster_runs, *, can_edit_templates=False, database_user=None):
    from mn_molecular_design_studio.pipeline_blocks import blank_program

    with st.container(key="pipeline-library"):
        section_key = prefix + "-library-section"
        if next_section := st.session_state.pop(prefix + "-library-next-section", None):
            st.session_state[section_key] = next_section
        st.session_state.setdefault(section_key, "Saved pipelines")
        title, create = st.columns([3, 1], vertical_alignment="bottom")
        title.title("Pipeline library")
        if create.button("New pipeline", type="primary", icon=":material/add:", width="stretch", key=prefix + "-new"):
            graph = blank_program()
            _activate_graph(prefix, graph)
            st.session_state[prefix + "-guided-open"] = graph["id"]
            st.rerun()
        section = st.segmented_control(
            "Pipeline library view", ["Saved pipelines", "Templates", "Runs", "Capabilities", "Trash"],
            selection_mode="single", label_visibility="collapsed", key=section_key,
        ) or "Saved pipelines"
        if notice := st.session_state.pop(prefix + "-library-notice", None):
            st.success(notice)
        if section == "Capabilities":
            from mn_molecular_design_studio.design_capabilities_ui import render_capabilities

            render_capabilities(prefix + "-design-capabilities",
                                lambda graph: _open_capability_graph(graph, prefix),
                                database_user=database_user, project_id=project_id)
            return
        if section == "Runs":
            _render_pipeline_runs(prefix, project_id, search_jobs, search_cluster_runs)
            return
        if section == "Trash":
            st.subheader("Trashed pipelines")
            if not trashed:
                st.caption("Trash is empty." if project_id else "Choose a project to open its Trash.")
            for item in trashed:
                graph = item["graph"]
                desc, restore = st.columns([3, 1], vertical_alignment="center")
                desc.write(graph["name"])
                desc.caption(f"Moved to Trash {item['trashed_at'][:10]} · {graph['id'][:8]}")
                if restore.button("Restore pipeline", key=prefix + "-restore-trash-" + graph["id"]):
                    _change_pipeline_trash(prefix, item, project_id, submit_api, trashed=False)
            return
        if section == "Templates":
            shared = list_shared_pipeline_templates()
            cards = st.columns(3)
            for index, item in enumerate(shared):
                template_id = item["id"]
                with cards[index % len(cards)], st.container(border=True):
                    if st.button(item["name"], key=prefix + "-template-" + template_id, width="stretch",
                                 help="Open the template canvas"):
                        try:
                            _open_template_view(template_id, prefix, project_id)
                            st.rerun()
                        except (ValueError, OSError, LookupError) as exc:
                            st.error(str(exc))
                    if item.get("description"):
                        st.caption(item["description"])
                    if item.get("saved_at"):
                        st.caption(f"Updated {item['saved_at'][:10]} UTC")

            saved_templates = list_pipeline_templates(owner_id, project_id) if project_id else []
            if saved_templates:
                st.subheader("Previously saved in this project")
                columns = st.columns(3)
                for index, item in enumerate(saved_templates):
                    with columns[index % len(columns)], st.container(border=True):
                        st.write(item["name"])
                        st.caption(f"{item['steps']} {'step' if item['steps'] == 1 else 'steps'} · Saved {item['created_at'][:10]}")
                        if st.button("Create from template", key=prefix + "-use-saved-template-" + item["id"], width="stretch"):
                            try:
                                result = submit_api(f"/api/v1/pipelines/saved-templates/{item['id']}/use",
                                                    {"project_id": project_id})
                                _activate_graph(prefix, result["graph"])
                                st.session_state[prefix + "-restore-viewport"] = True
                                st.rerun()
                            except (ValueError, RuntimeError, OSError, LookupError, KeyError) as exc:
                                st.error(str(exc))
            return

        if not project_id:
            st.info("Choose a project to access saved pipelines and run history.")
        elif not definitions:
            st.info("No saved pipelines in this project. Start a new canvas or choose a template.")
        current = st.session_state.get(prefix + "-graph")
        if current and not any(item["graph"]["id"] == current["id"] for item in (*definitions, *archived, *trashed)):
            if st.button(f"Continue editing {current['name']}", key=prefix + "-resume"):
                st.session_state[prefix + "-view"] = "editor"
                st.session_state[prefix + "-restore-viewport"] = True
                st.rerun()
        if definitions:
            names = Counter(row["graph"]["name"].casefold() for row in definitions)
            query = st.text_input("Search saved pipelines", placeholder="Search by pipeline name",
                                  key=prefix + "-library-search").strip().casefold()
            visible = [row for row in definitions if query in row["graph"]["name"].casefold()]
            if not visible:
                st.caption("No pipelines match this search.")
            cards = st.columns(3)
            for index, draft in enumerate(visible):
                graph = draft["graph"]
                steps = len(graph["nodes"])
                versions = draft["versions"]
                with cards[index % len(cards)], st.container(border=True):
                    if st.button(graph["name"], key=prefix + "-open-" + graph["id"], width="stretch"):
                        _open_saved_pipeline(prefix, draft)
                        st.rerun()
                    details, actions = st.columns([2, 1], vertical_alignment="center")
                    details.caption(f"{steps} {'step' if steps == 1 else 'steps'} · {versions} saved {'version' if versions == 1 else 'versions'}")
                    updated = draft["saved_at"][:19].replace("T", " ")
                    identifier = f" · {graph['id'][:8]}" if names[graph["name"].casefold()] > 1 else ""
                    details.caption(f"Updated {updated} UTC{identifier}")
                    with actions.popover("Actions", width="stretch"):
                        if can_edit_templates and st.button("Save as shared template", key=prefix + "-save-template-" + graph["id"]):
                            st.session_state[prefix + "-template-source"] = graph
                        if st.button("Move to Trash", key=prefix + "-trash-" + graph["id"]):
                            _change_pipeline_trash(prefix, draft, project_id, submit_api, trashed=True)

        if archived:
            with st.expander("Archived pipelines"):
                selected = st.selectbox("Archived pipeline", archived, format_func=lambda row: row["graph"]["name"],
                                        key=prefix + "-archived-pipeline")
                if st.button("Restore pipeline", key=prefix + "-unarchive"):
                    if _change_archive_status(selected, project_id, submit_api, archived=False):
                        _activate_graph(prefix, selected["graph"], selected["revision"])
                        st.rerun()
        if can_edit_templates:
            _show_save_pipeline_template(prefix, submit_api)


def _show_pipeline_versions(prefix, owner_id, project_id, graph, draft, submit_api):
    open_key = prefix + "-versions-open"
    if not st.session_state.get(open_key):
        return

    @st.dialog("Pipeline versions", width="medium", on_dismiss=lambda: st.session_state.pop(open_key, None))
    def dialog():
        st.write(graph["name"])
        if not project_id:
            st.info("Choose a project before saving a version.")
            return
        if st.button("Save a version", type="primary", key=prefix + "-version-save"):
            try:
                _save_draft(prefix, project_id, graph, submit_api, version=True)
            except StudioAPIConflictError:
                st.session_state.pop(open_key, None)
                st.rerun()
            except (ValueError, OSError, RuntimeError) as exc:
                st.error(str(exc))
            else:
                st.success("Version saved.")
        versions = list_versions(owner_id, project_id, graph["id"])
        if versions:
            version = st.selectbox("Saved version", versions, key=prefix + "-version-" + graph["id"],
                format_func=lambda row: f"{row['saved_at'][:19].replace('T', ' ')} UTC · {row['revision'][:10]}")
            if st.button("Restore to editor", key=prefix + "-restore"):
                try:
                    _save_draft(prefix, project_id, version["graph"], submit_api)
                except StudioAPIConflictError:
                    st.session_state.pop(open_key, None)
                    st.rerun()
                except (ValueError, OSError, RuntimeError) as exc:
                    st.error(str(exc))
                else:
                    _activate_graph(prefix, version["graph"], _draft_revision(prefix, graph))
                    st.session_state.pop(open_key, None)
                    st.rerun()
        else:
            st.caption("No saved versions yet.")
        if draft and st.button("Archive pipeline", key=prefix + "-archive"):
            current = dict(draft, revision=_draft_revision(prefix, graph) or draft["revision"])
            if _change_archive_status(current, project_id, submit_api, archived=True):
                for suffix in ("-graph", "-draft-base", "-draft-conflict", "-versions-open"):
                    st.session_state.pop(prefix + suffix, None)
                st.session_state[prefix + "-view"] = "library"
                st.rerun()

    dialog()


def _screen_run_requirements(graph):
    fields = {
        "active_library": ("Screening library", ("library_id", "feature_manifest"),
                           "Choose a prepared compound library and fingerprint dataset."),
        "active_seed": ("Initial centroids", ("library_id", "centroid_manifest"),
                        "Choose a registered centroid dataset."),
        "active_target": ("Receptor", ("target_pdb",), "Choose a receptor structure and review its docking box."),
        "active_docking": ("Docking workstations", ("worker_nodes",), "Select at least one docking workstation."),
    }
    nodes = {node["data"]["kind"]: node for node in graph["nodes"]}
    for kind, (label, required, missing) in fields.items():
        node = nodes.get(kind)
        if node is None:
            yield label, None, "Restore the " + label.lower() + " step.", True
            continue
        params = node["data"]["params"]
        if not all(params.get(field) for field in required):
            yield label, node["id"], missing, True
            continue
        detail = str(params.get("target_name") or "") if kind == "active_target" else ""
        if kind == "active_docking":
            count = len(params["worker_nodes"])
            detail = f"{params['engine']} · {count} {'workstation' if count == 1 else 'workstations'}"
        if not detail:
            detail = Path(params[required[-1]]).name
        yield label, node["id"], detail, False


def _render_pipeline_run_setup(prefix, project_id, graph, submit_api, execution_available):
    from mn_molecular_design_studio.cluster.ui import _job_details_reference

    view_key = prefix + "-view"
    st.button("Back to canvas", icon=":material/arrow_back:", key=prefix + "-run-back",
              on_click=st.session_state.update, args=({view_key: "editor"},))
    st.title("Run setup")
    st.caption(graph["name"])
    screen = is_screening_graph(graph)
    inputs = list(_screen_run_requirements(graph)) if screen else []
    blocked = any(row[3] for row in inputs)
    graph_problem = ""
    try:
        plan = validate_graph(graph, require_inputs=True, require_files=False)
    except PipelineInputError as exc:
        blocked = True
        if exc.node_id not in {row[1] for row in inputs if row[3]}:
            node = next((node for node in graph["nodes"] if node["id"] == exc.node_id), None)
            inputs.append((node["data"]["label"] if node else "Input", exc.node_id, str(exc), True))
    except ValueError as exc:
        blocked = True
        graph_problem = str(exc)

    input_column, execution_column = st.columns([1.5, 1], gap="large")
    with input_column, st.container(border=True):
        st.subheader("Inputs")
        if graph_problem:
            st.warning(graph_problem)
        if not blocked and not inputs:
            st.write(f"{len(plan['order'])} steps configured")
        for label, node_id, detail, missing in inputs:
            text, action = st.columns([3, 1], vertical_alignment="center")
            text.write(label)
            (text.warning if missing else text.caption)(detail)
            if node_id:
                action.button("Configure", key=prefix + "-run-input-" + node_id,
                              help=f"Configure {label}", on_click=_open_node_editor, args=(prefix, node_id))
    with execution_column, st.container(border=True):
        st.subheader("Execution")
        if screen:
            reuse = False
            st.write("Start a new screening campaign")
            if execution_available:
                from mn_molecular_design_studio.pipeline_node_ui import _available_docking_machines

                docking = next((node["data"]["params"] for node in graph["nodes"]
                                if node["data"]["kind"] == "active_docking"), {})
                gpu_engine = docking.get("engine") in {"UniDock Pro", "GNINA", "AutoDock-GPU"}
                try:
                    machines = _available_docking_machines(gpu_engine)
                except (OSError, ValueError, KeyError, TypeError):
                    st.caption("Could not read the workstation snapshot. Open Docking workstations to check compute.")
                else:
                    if machines:
                        with st.expander("Workstation snapshot"):
                            st.dataframe([{"Workstation": machine["name"], "Free CPUs": machine["cpu_free"],
                                           "Reported GPUs": len(machine["gpus"]),
                                           "Heartbeat": "Stale" if machine["heartbeat_stale"] else "Online"}
                                          for machine in machines], hide_index=True, width="stretch")
                    else:
                        st.caption("The latest snapshot shows no available workstation for this docking engine.")
        else:
            reuse = st.checkbox("Reuse unchanged completed steps", value=st.session_state.get(prefix + "-reuse", True),
                                key=prefix + "-run-reuse-" + graph["id"],
                                help="Clear this option to compute every step fresh.")
            st.session_state[prefix + "-reuse"] = reuse
            st.button("Workstation settings", key=prefix + "-run-resources", on_click=st.session_state.update,
                      args=({prefix + "-resources-open": True},))
        if not project_id:
            st.info("Choose a project to submit this run.")
        if not execution_available:
            st.info("Job submission is unavailable in this development preview.")

    with execution_column:
        request = {"project_id": project_id, "graph": graph, "reuse_completed": reuse, "gpu_device": "all"}
        check_key = prefix + "-run-check"
        signature = hashlib.sha256(json.dumps(request, sort_keys=True).encode()).hexdigest()
        saved_check = st.session_state.get(check_key) or {}
        check = saved_check.get("result", {}) if saved_check.get("signature") == signature else {}
        available = bool(project_id and execution_available and not blocked and not st.session_state.get(prefix + "-draft-conflict"))
        checks, queue = st.columns([1, 1])
        if checks.button("Check readiness", key=prefix + "-preflight", disabled=not available, width="stretch"):
            st.session_state.pop(check_key, None)
            check = {}
            try:
                with st.spinner("Checking saved inputs…"):
                    result = submit_api("/api/v1/pipelines/validate", request)
                if result.get("valid") is not True:
                    raise ValueError("The pipeline did not pass validation.")
                check = result
                st.session_state[check_key] = {"signature": signature, "result": check}
            except (ValueError, OSError, RuntimeError, TypeError) as exc:
                st.error(str(exc))
        if check.get("valid"):
            if check.get("ready_to_run"):
                st.success("Inputs verified. The scheduler will allocate compute when the run is queued.")
            else:
                st.info(check.get("file_validation") or "The assigned workstation will check saved files and compute before starting.")
        if queue.button("Submit run", key=prefix + "-submit-run", type="primary",
                        disabled=not (available and check.get("valid")), width="stretch"):
            try:
                with st.spinner("Submitting pipeline…"):
                    job = submit_api("/api/v1/pipelines/run", request)["job"]
                st.session_state[prefix + "-last-job"] = job
                st.session_state[prefix + "-last-job-graph"] = graph["id"]
                st.session_state["last_submitted_job"] = job
                st.session_state.pop(check_key, None)
                st.session_state[view_key] = "runs"
                st.session_state[prefix + "-run-notice"] = job.get("job_code") or _job_details_reference(job)
                st.rerun()
            except (ValueError, OSError, RuntimeError, TypeError, KeyError) as exc:
                st.error(str(exc))


def render_pipeline_builder(*, owner_id, project_selector, submit_api, search_jobs,
                            read_only_projects=(), execution_available=True, search_cluster_runs=None,
                            library_catalog_available=True, save_project_selector=None, can_edit_templates=False,
                            database_user=None):
    st.markdown("<style>.block-container{max-width:1900px;padding-top:3.75rem} [data-testid=stMainBlockContainer]{padding-left:2rem;padding-right:2rem}</style>", unsafe_allow_html=True)
    template_key = _TEMPLATE_VIEW_PREFIX + "-open"
    opened_template = st.session_state.get(template_key)
    template_read_only = bool(opened_template and not can_edit_templates)
    starting_points = shared_template_starting_points(opened_template["id"]) if opened_template else ()
    project_id = None if opened_template else project_selector("visual-pipeline")
    prefix = _TEMPLATE_VIEW_PREFIX if opened_template else f"visual-pipeline-{project_id or 'unassigned'}"
    state_key = prefix + "-graph"
    saved_drafts = list_definitions(owner_id, project_id, include_archived=True, include_trashed=True) if project_id else []
    definitions = [item for item in saved_drafts if not (item.get("archived_at") or item.get("trashed_at"))]
    archived = [item for item in saved_drafts if item.get("archived_at")]
    trashed = [item for item in saved_drafts if item.get("trashed_at")]
    if project_id in read_only_projects:
        st.title("Pipeline library")
        st.info("This project is read-only. You can inspect and download its saved pipeline graphs.")
        if st.button("Browse shared templates", key=prefix + "-browse-shared-templates"):
            st.session_state.setdefault("project-pending-selection", {})["visual-pipeline-project-id"] = ""
            st.session_state["visual-pipeline-unassigned-view"] = "library"
            st.session_state["visual-pipeline-unassigned-library-next-section"] = "Templates"
            st.rerun()
        if saved_drafts:
            draft = st.selectbox("Saved pipeline", saved_drafts, format_func=lambda item: item["graph"]["name"],
                                 key=prefix + "-archived-draft")
            versions = list_versions(owner_id, project_id, draft["graph"]["id"])
            snapshot = st.selectbox("Revision", [draft, *versions], key=prefix + "-archived-revision",
                format_func=lambda item: ("Saved pipeline" if "versions" in item else "Saved version")
                + " · " + item["saved_at"][:19] + " · " + item["revision"][:10])
            st.download_button("Download pipeline graph", json.dumps(snapshot["graph"], indent=2),
                               file_name=f"pipeline-{snapshot['revision'][:10]}.json", mime="application/json")
            with st.expander("Graph and parameters"):
                st.json(snapshot["graph"], expanded=False)
        else:
            st.caption("No saved pipelines in this project.")
        return
    requested_pipeline = st.query_params.get("pipeline")
    request_key = prefix + "-linked-pipeline"
    if not opened_template and requested_pipeline and st.session_state.get(request_key) != requested_pipeline:
        st.session_state[request_key] = requested_pipeline
        linked = next((draft for draft in definitions if draft["graph"]["id"] == requested_pipeline), None)
        if linked:
            _open_saved_pipeline(prefix, linked)
        else:
            st.info("This saved pipeline is unavailable in the selected project.")
    if state_key not in st.session_state:
        preview_prefix = f"pipeline-blocks-{project_id or 'unassigned'}"
        from_preview = st.session_state.get(preview_prefix + "-graph")
        if from_preview:
            _activate_graph(prefix, from_preview, _draft_revision(preview_prefix, from_preview))
    reuse_key = prefix + "-reuse"
    st.session_state.setdefault(reuse_key, True)
    pending = st.session_state.get("pipeline-pending-draft")
    if not opened_template and pending and pending["project_id"] == project_id:
        _activate_graph(prefix, pending["graph"])
        st.session_state.pop("pipeline-pending-draft", None)
    current_graph = st.session_state.get(state_key)
    view_key = prefix + "-view"
    view = st.session_state.get(view_key, "editor" if current_graph else "library")
    if view == "library" or current_graph is None:
        _render_pipeline_library(prefix, owner_id, project_id, definitions, archived, trashed,
                                 submit_api, search_jobs, search_cluster_runs, can_edit_templates=can_edit_templates, database_user=database_user)
        return
    if view == "runs":
        _render_pipeline_runs(prefix, project_id, search_jobs, search_cluster_runs, graph=current_graph)
        return
    if not opened_template:
        migrated_graph = migrate_template_graph(current_graph)
        if migrated_graph is not current_graph:
            _activate_graph(prefix, migrated_graph, _draft_revision(prefix, current_graph))

    graph = st.session_state[state_key]
    saved_draft = next((item for item in definitions if item["graph"]["id"] == graph["id"]), None)
    in_trash = any(item["graph"]["id"] == graph["id"] for item in trashed)
    if in_trash:
        st.warning("This saved pipeline is in Trash. Restore it to edit or run it.")
        if view == "run-setup":
            if st.button("Back to pipeline library", key=prefix + "-trashed-run-back"):
                st.session_state[view_key] = "library"
                st.session_state[prefix + "-library-next-section"] = "Trash"
                st.rerun()
            return
    if view == "run-setup":
        _render_pipeline_run_setup(prefix, project_id, graph, submit_api, execution_available)
        _show_pipeline_resources(prefix, owner_id, project_id, submit_api)
        _show_node_editor(prefix, owner_id=owner_id, project_id=project_id, submit_api=submit_api,
                          search_jobs=search_jobs, library_catalog_available=library_catalog_available)
        return
    if st.session_state.get(prefix + "-draft-conflict") == graph["id"]:
        st.warning("The saved pipeline has changed. Your edits are still in this tab. Open the latest saved pipeline or save your edits as a copy.")
        if saved_draft:
            st.button("Open latest saved pipeline", key=prefix + "-open-latest", on_click=_activate_graph,
                args=(prefix, saved_draft["graph"], saved_draft["revision"]))
        if st.button("Save edits as a copy", key=prefix + "-copy-edits"):
            st.session_state[prefix + "-save-as-open"] = True

    from mn_molecular_design_studio import candidate_execution
    settings = candidate_execution.policy(graph.get("candidate_execution"))
    graph["candidate_execution"] = settings
    if opened_template:
        if can_edit_templates:
            st.caption("Saving updates the open route for everyone." if starting_points
                       else "Saving updates the shared template for everyone.")
        else:
            st.caption("Choose Use template to create an editable pipeline.")
        if notice := st.session_state.pop(prefix + "-save-notice", None):
            st.success(notice)
    elif not project_id:
        st.info("Save lets you choose a project for this pipeline. Run setup also requires a project.")
    launch_panel = st.container()
    last_job = st.session_state.get(prefix+"-last-job", {})
    if st.session_state.get(prefix + "-last-job-graph") != graph["id"]:
        last_job = {}
    execution = {}
    if last_job.get("run_dir"):
        last_path = Path(last_job["run_dir"])
        execution = _matching_execution(graph, read_json(last_path/"graph_snapshot.json"), read_json(last_path/"graph_execution.json"))
    reset_token = st.session_state[prefix + "-reset"]
    from mn_molecular_design_studio.tool_connections import editor_connections, editor_tool_catalog
    round_outline = None
    screen_issues = screening_graph_issues(graph) if is_screening_graph(graph) else None
    if screen_issues == []:
        from mn_molecular_design_studio.pipeline_blocks import program
        try:
            round_outline = program(graph)
        except (KeyError, ValueError):
            pass  # The runnable-graph check below explains the draft issue.
    elif screen_issues:
        from mn_molecular_design_studio.pipeline_graph import restore_screening_connections
        try:
            restored = restore_screening_connections(graph)
            validate_graph(restored, allow_incomplete=True)
        except ValueError:
            pass  # The checklist below explains the missing stages.
        else:
            if not template_read_only and restored != graph and st.button("Restore screening connections", key=prefix + "-restore-screen"):
                _activate_graph(prefix, restored, _draft_revision(prefix, graph))
                if opened_template:
                    st.session_state[prefix + "-restore-viewport"] = True
                st.rerun()
    value = _component()(graph=graph, catalog=editor_tool_catalog(), screening_positions=ACTIVE_LEARNING_LAYOUT,
                         connections=editor_connections(), round_outline=round_outline,
                         reset_token=reset_token, input_contexts={},
                         draft_saved=bool(_draft_revision(prefix, graph)),
                         shared_template=bool(opened_template), template_read_only=template_read_only,
                         template_id=opened_template["id"] if opened_template else None,
                         starting_points=starting_points,
                         can_edit_templates=can_edit_templates,
                         binder_picker=not opened_template and not in_trash,
                         guided_setup=not opened_template and not graph["nodes"],
                         restore_viewport=bool(_draft_revision(prefix, graph) or st.session_state.get(prefix + "-restore-viewport")),
                         execution=execution, settings_update=st.session_state.get(prefix+"-settings-update"),
                         key=f"{prefix}-canvas-{reset_token}", default=None)
    if value and (value.get("graph") or {}).get("id") == graph["id"]:
        signature = json.dumps(value, sort_keys=True)
        if signature != st.session_state.get(prefix+"-last-value"):
            st.session_state[prefix+"-last-value"] = signature
            try:
                component_graph = copy.deepcopy(graph if template_read_only else value["graph"])
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
                action=value.get("action") or {}
                new_action = action.get("id") and action["id"] != st.session_state.get(prefix+"-handled-action")
                action_type = action.get("type") if new_action else None
                changed = dict(updated, viewport=graph.get("viewport")) != graph
                if project_id and action_type not in {"save_as", "save_template"} and (changed or action_type == "save"):
                    _save_draft(prefix, project_id, updated, submit_api, from_canvas=True)
                st.session_state[state_key] = graph = updated
                if new_action:
                    st.session_state[prefix+"-handled-action"] = action["id"]
                    if action_type == "edit_node" and not template_read_only:
                        _open_node_editor(prefix, action.get("node_id"))
                    elif opened_template:
                        if action_type == "template_route":
                            requested = action.get("template_id")
                            allowed = {route["id"] for point in starting_points for route in point["routes"]}
                            if requested in allowed and requested != opened_template["id"]:
                                working_routes = opened_template["working_routes"]
                                working_routes[opened_template["id"]] = {
                                    "graph": graph, "revision": opened_template["revision"],
                                    "restore_viewport": can_edit_templates or st.session_state.get(prefix + "-restore-viewport"),
                                }
                                _open_template_view(requested, opened_template["return_prefix"],
                                                    opened_template["return_project"], working_routes=working_routes)
                                st.rerun()
                        if action_type == "save" and can_edit_templates:
                            saved = submit_api("/api/v1/pipelines/shared-templates/save", {
                                "graph": graph, "template_id": opened_template["id"],
                                "expected_revision": opened_template["revision"],
                            })
                            opened_template["revision"] = saved["revision"]
                            _activate_graph(prefix, saved["graph"])
                            st.session_state[prefix + "-restore-viewport"] = True
                            st.session_state[prefix + "-save-notice"] = f"Saved shared template {saved['name']}."
                            st.rerun()
                        if action_type in {"library", "use_template"}:
                            destination = opened_template["return_prefix"]
                            st.session_state.pop(template_key, None)
                            st.session_state.setdefault("project-pending-selection", {})["visual-pipeline-project-id"] = (
                                opened_template["return_project"] or ""
                            )
                            if action_type == "use_template":
                                _activate_graph(destination, {**graph, "id": uuid.uuid4().hex})
                                st.session_state[destination + "-restore-viewport"] = (
                                    can_edit_templates or st.session_state.get(prefix + "-restore-viewport", False)
                                )
                            else:
                                st.session_state[destination + "-view"] = "library"
                                st.session_state[destination + "-library-next-section"] = "Templates"
                            st.rerun()
                    elif action_type == "guided_setup" and not graph["nodes"]:
                        st.session_state[prefix + "-guided-open"] = graph["id"]
                    elif action_type == "binder_picker":
                        st.session_state[prefix + "-binder-open"] = action.get("node_id")
                    elif action["type"] == "campaign":
                        st.session_state[prefix + "-campaign-open"] = True
                    elif action["type"] == "resources":
                        st.session_state[prefix + "-resources-open"] = True
                    elif action["type"] == "versions":
                        st.session_state[prefix + "-versions-open"] = True
                    elif action["type"] in {"library", "runs", "run"}:
                        if not st.session_state.get(prefix + "-draft-conflict"):
                            st.session_state[view_key] = "run-setup" if action["type"] == "run" else action["type"]
                            st.rerun()
                    elif action["type"] == "save":
                        if not project_id:
                            st.session_state[prefix + "-save-project-open"] = True
                        else:
                            st.success(f"Saved {graph['name']}.")
                            link = urlencode({"pipeline": graph["id"], "pipeline_project": project_id})
                            st.markdown(f"[Open saved pipeline](?{link})")
                    elif action_type == "save_template" and can_edit_templates:
                        st.session_state[prefix + "-template-source"] = graph
                    elif action["type"] == "save_as":
                        if not project_id:
                            st.warning("Choose a project before saving a pipeline.")
                        else:
                            st.session_state[prefix + "-save-as-open"] = True
                if changed and not action_type and (screen_issues is not None or is_screening_graph(updated)):
                    st.rerun()
            except StudioAPIConflictError:
                st.rerun()
            except (ValueError, RuntimeError, OSError, TypeError, KeyError) as exc:
                launch_panel.error(str(exc))
    if opened_template:
        if can_edit_templates:
            _show_node_editor(prefix, owner_id=owner_id, project_id=None, submit_api=submit_api,
                              search_jobs=search_jobs, stage_edits=True, library_catalog_available=library_catalog_available)
        return
    from mn_molecular_design_studio.cluster.ui import _job_details_reference, _navigate_to_job_details
    with launch_panel:
        if screen_issues is not None:
            if screen_issues:
                with st.expander("Screening graph checklist"):
                    for issue in screen_issues:
                        st.write(issue)
            else:
                st.caption("Screening graph complete · includes next-round feedback")
        try:
            plan = validate_graph(graph, require_inputs=True, require_files=False)
        except PipelineInputError as exc:
            node = next(n for n in graph["nodes"] if n["id"] == exc.node_id)
            st.info(str(exc))
            st.button(f"Configure {node['data']['label']}", key=prefix + "-configure-required",
                      on_click=_open_node_editor, args=(prefix, node["id"]))
        except ValueError as exc:
            st.info(f"Before running · {exc}")
        else:
            if project_id:
                st.caption(f"{len(plan['order'])} steps configured")
        if last_job:
            reference = _job_details_reference(last_job)
            st.caption(f"Last submission · {last_job.get('job_code') or reference}")
            if st.button("Open run results", key=prefix + "-open-latest-result"):
                _navigate_to_job_details(reference, final_results=graph.get("template_id") == "unidock_active_learning")
    _show_guided_pipeline_setup(prefix, graph, project_id, submit_api)
    if not opened_template and not in_trash:
        _show_binder_branch(prefix, graph, project_id, submit_api)
    _show_pipeline_campaign(prefix, project_id, submit_api)
    _show_pipeline_resources(prefix, owner_id, project_id, submit_api)
    _show_pipeline_versions(prefix, owner_id, project_id, graph, saved_draft, submit_api)
    if can_edit_templates:
        _show_save_pipeline_template(prefix, submit_api)
    _show_pipeline_save_as(prefix, project_id, submit_api)
    _show_pipeline_save_project(prefix, graph, save_project_selector or project_selector, submit_api, read_only_projects)
    _show_node_editor(prefix, owner_id=owner_id, project_id=project_id, submit_api=submit_api,
                      search_jobs=search_jobs, library_catalog_available=library_catalog_available)


def _pipeline_run_graph(row):
    if not row.get("run_dir"):
        return {}
    run_dir = Path(row["run_dir"])
    snapshot = read_json(run_dir / "graph_snapshot.json")
    if snapshot:
        return snapshot
    params = read_json(run_dir / "input.json").get("params")
    return (params.get("graph") or {}) if isinstance(params, dict) else {}


def _render_pipeline_runs(prefix, project_id, search_jobs, search_cluster_runs, *, graph=None):
    from mn_molecular_design_studio.cluster import config
    from mn_molecular_design_studio.cluster.ui import _job_details_reference, _navigate_to_job_details

    if graph:
        st.button("Back to canvas", icon=":material/arrow_back:", key=prefix + "-runs-back",
                  on_click=st.session_state.update, args=({prefix + "-view": "editor"},))
        st.title("Pipeline runs")
        st.caption(graph["name"])
    submitted = st.session_state.pop(prefix + "-run-notice", None)
    if submitted:
        st.success(f"Run submitted · {submitted}")
    if not project_id:
        st.info("Choose a project to view pipeline runs.")
        return
    if st.button("Refresh runs", key=prefix + "-refresh-runs", icon=":material/refresh:"):
        st.rerun()

    def pages(search, **filters):
        page = search(project_id=project_id, tools=["visual_pipeline"], limit=100, offset=0, **filters)
        rows = list(page["items"])
        while len(rows) < page.get("total", len(rows)):
            more = search(project_id=project_id, tools=["visual_pipeline"], limit=100, offset=len(rows), **filters)["items"]
            if not more:
                break
            rows.extend(more)
        return rows

    local_rows = pages(search_jobs, task_group="campaign")
    remote_rows = []
    if search_cluster_runs:
        try:
            remote_rows = [{key: value for key, value in row.items() if key != "run_dir"}
                           for row in pages(search_cluster_runs)]
        except (OSError, RuntimeError, ValueError) as exc:
            st.caption(f"Could not load shared run history. {exc}")
    local_id = config.identity()
    local_refs = {_job_details_reference(row) for row in local_rows}
    rows = local_rows + [row for row in remote_rows
                         if row.get("node_id") != local_id or _job_details_reference(row) not in local_refs]
    last_job = st.session_state.get(prefix + "-last-job")
    if last_job and (not graph or st.session_state.get(prefix + "-last-job-graph") == graph["id"]):
        recent = _job_details_reference(last_job)
        if recent not in {_job_details_reference(row) for row in rows}:
            rows.insert(0, {**last_job, "pipeline_graph_id": st.session_state.get(prefix + "-last-job-graph")})

    if graph:
        rows = [row for row in rows
                if (row.get("pipeline_graph_id") or _pipeline_run_graph(row).get("id")) == graph["id"]]
    if not rows:
        st.info("No runs recorded for this pipeline." if graph else "No pipeline runs in this project yet.")
        return
    rows.sort(key=lambda row: str(row.get("created_at") or ""), reverse=True)
    page_size = 15
    page_count = (len(rows) + page_size - 1) // page_size
    page_number = 1
    if page_count > 1:
        page_number = st.number_input("Page", min_value=1, max_value=page_count, value=1,
                                      key=prefix + "-runs-page-" + (graph["id"] if graph else "all"))
    st.caption(f"{len(rows)} {'run' if len(rows) == 1 else 'runs'}")
    rows = rows[(page_number - 1) * page_size:page_number * page_size]
    for row in rows:
        reference = _job_details_reference(row)
        name, status, action = st.columns([3, 1, 1], vertical_alignment="center")
        name.write(config.display_text(row.get("pipeline_name") or row.get("campaign_name") or (graph or {}).get("name") or "Pipeline run"))
        job_code = row.get("job_code") or reference
        name.caption(" · ".join(value for value in (job_code, str(row.get("created_at") or "")[:10], row.get("node_name")) if value))
        status.write(str(row.get("status") or "submitted").replace("_", " ").capitalize())
        if action.button("Open results", key=prefix + "-results-" + reference, width="stretch"):
            snapshot = _pipeline_run_graph(row)
            _navigate_to_job_details(reference, final_results=(snapshot or graph or {}).get("template_id") == "unidock_active_learning")

    inspectable = [row for row in rows if row.get("run_dir")]
    if inspectable:
        selected = st.selectbox("Inspect run progress", inspectable, index=None, placeholder="Choose a run",
                                format_func=lambda row: row.get("job_code") or _job_details_reference(row),
                                key=prefix + "-run-progress")
        if selected:
            _render_pipeline_run_progress(Path(selected["run_dir"]), _pipeline_run_graph(selected))


def _render_pipeline_run_progress(run_dir, snapshot):
    execution = read_json(run_dir / "graph_execution.json")
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
