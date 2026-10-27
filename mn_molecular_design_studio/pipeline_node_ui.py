"""Familiar tool input forms embedded in the visual pipeline workspace."""
from __future__ import annotations

import copy
import json
from functools import lru_cache
from html import escape
from pathlib import Path

import streamlit as st

from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG


@lru_cache(maxsize=32)
def node_form_defaults(kind):
    """Read engine defaults from the same schema used by standalone tool pages."""
    spec = NODE_CATALOG[kind]
    result = {f["name"]: copy.deepcopy(f["default"]) for f in spec["fields"]}
    if kind == "boltz2":
        result["prediction_config_version"] = 2
    if not spec["workflow"]:
        return result
    from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import workflow_schema
    native = {p["name"]: p.get("default") for p in workflow_schema(spec["workflow"])["parameters"]}
    aliases = {"sequences": "num_seq_per_target" if kind in {"foundry_mpnn", "ligandmpnn"} else "num_mpnn_seqs",
               "num_designs": "max_trajectories" if kind == "bindcraft2_generator" else "number_of_final_designs" if kind == "bindcraft2" else "num_designs"}
    for field in spec["fields"]:
        source = aliases.get(field["name"], field["name"])
        if native.get(source) is not None:
            result[field["name"]] = copy.deepcopy(native[source])
    if kind == "rfd3":
        result["binder_length"] = f"{native['length_min']}-{native['length_max']}"
    return result


def compatible_node_replacements(node, graph):
    """List same-category nodes that can retain the current connections."""
    current_kind = node["data"]["kind"]
    current_spec = NODE_CATALOG[current_kind]
    current_group = current_spec.get("replacement_group", current_spec.get("category"))
    nodes = {item["id"]: item for item in graph.get("nodes", [])}
    task_input = next((item for item in nodes.values() if item["data"].get("kind") == "cofold_input"), None)
    task_type = (task_input or {}).get("data", {}).get("params", {}).get("task_type")
    ligand_input_type = (task_input or {}).get("data", {}).get("params", {}).get("ligand_input_type", "smiles")
    if not task_type:
        template_id = graph.get("template_id")
        if template_id in {"all_generators_boltz2", "vhh_design"}:
            task_type = "protein_protein"
    target_structure_available = bool(
        (task_input or {}).get("data", {}).get("params", {}).get("target_pdb")
        or any(item["data"].get("kind") == "target" and item["data"].get("params", {}).get("target_pdb") for item in nodes.values())
    )
    incoming = [nodes[edge["source"]] for edge in graph.get("edges", [])
                if edge.get("target") == node["id"] and edge.get("source") in nodes]
    outgoing = [nodes[edge["target"]] for edge in graph.get("edges", [])
                if edge.get("source") == node["id"] and edge.get("target") in nodes]

    def output_type(item):
        outputs = NODE_CATALOG[item["data"]["kind"]].get("outputs") or []
        return outputs[0].get("type") if outputs else None

    def accepts(output, input_spec):
        accepted = input_spec.get("accepts", [])
        return output == "candidates" if accepted == ["target"] else output in accepted

    choices = []
    for kind, candidate in NODE_CATALOG.items():
        candidate_group = candidate.get("replacement_group", candidate.get("category"))
        if kind == current_kind or candidate.get("hidden") or candidate_group != current_group:
            continue
        if task_type and candidate.get("molecule_tasks") and task_type not in candidate["molecule_tasks"]:
            continue
        if task_type == "protein_ligand" and candidate.get("ligand_input_formats") and ligand_input_type not in candidate["ligand_input_formats"]:
            continue
        if candidate.get("requires_target_structure") and not target_structure_available:
            continue
        inputs = candidate.get("inputs") or []
        outputs = candidate.get("outputs") or []
        if incoming:
            if not inputs or (len(incoming) > 1 and not inputs[0].get("multiple")):
                continue
            if any(not accepts(output_type(parent), inputs[0]) for parent in incoming):
                continue
        if outgoing:
            if not outputs:
                continue
            candidate_output = outputs[0].get("type")
            if any(not (NODE_CATALOG[target["data"]["kind"]].get("inputs")
                        and accepts(candidate_output, NODE_CATALOG[target["data"]["kind"]]["inputs"][0]))
                   for target in outgoing):
                continue
        choices.append(kind)
    return choices


def replacement_node_data(node, new_kind):
    """Build replacement settings, retaining only fields valid for both tools."""
    old_spec = NODE_CATALOG[node["data"]["kind"]]
    new_spec = NODE_CATALOG[new_kind]
    values = node_form_defaults(new_kind)
    old_fields = {item["name"]: item for item in old_spec["fields"]}
    new_fields = {item["name"]: item for item in new_spec["fields"]}
    old_values = node["data"].get("params", {})
    for name, value in old_values.items():
        before, after = old_fields.get(name), new_fields.get(name)
        if not before or not after or before["type"] != after["type"]:
            continue
        field_type = after["type"]
        if field_type == "select" and value not in after.get("options", []):
            continue
        if field_type == "number":
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                continue
            if "min" in after and value < after["min"] or "max" in after and value > after["max"]:
                continue
            if isinstance(after["default"], int) and not float(value).is_integer():
                continue
        elif field_type == "boolean" and type(value) is not bool:
            continue
        elif field_type == "text" and not isinstance(value, str):
            continue
        elif field_type == "json" and not isinstance(value, list):
            continue
        values[name] = copy.deepcopy(value)

    old_label = node["data"].get("label") or old_spec["label"]
    label = new_spec["label"] if old_label == old_spec["label"] else old_label
    return {"kind": new_kind, "label": label, "description": new_spec["description"], "params": values}


def _connected_inputs(graph, node_id):
    nodes = {n["id"]: n for n in graph["nodes"]}
    parents = [nodes[e["source"]] for e in graph["edges"] if e["target"] == node_id and e["source"] in nodes]
    seen, upstream, queue = set(), [], list(parents)
    while queue:
        node = queue.pop(0)
        if node["id"] in seen:
            continue
        seen.add(node["id"])
        if node["data"]["kind"] in {"target", "binder_context", "design_source", "generation_source", "backbone_source", "molecule_source", "project_design_source"}:
            upstream.append(node)
        queue.extend(nodes[e["source"]] for e in graph["edges"] if e["target"] == node["id"] and e["source"] in nodes)
    return parents, upstream


def _render_connection_card(node, graph):
    parents, inputs = _connected_inputs(graph, node["id"])
    with st.container(border=True):
        st.markdown("**Connected inputs**")
        if not parents:
            st.caption("Connect an upstream block to provide this step’s input.")
        else:
            st.write(" → ".join(p["data"]["label"] for p in parents) + " → " + node["data"]["label"])
        for source in inputs:
            params = source["data"]["params"]
            if source["data"]["kind"] == "target":
                path = params.get("target_pdb")
                a, b = st.columns(2)
                a.caption("TARGET STRUCTURE")
                a.write(Path(path).name if path else "Choose a structure in the Target block")
                b.caption("CHAINS · HOTSPOTS · COLDSPOTS")
                b.write(f"{params.get('target_chains') or '—'} · {params.get('hotspots') or 'No hotspots'} · {params.get('coldspots') or 'No coldspots'}")
            elif source["data"]["kind"] == "binder_context":
                empty_state = "No annotation" if params.get("clear_unselected") else "Inherited"
                context_name = source["data"].get("label", "Binder Design Context").upper()
                st.caption(context_name + " · " + f"{('Hotspots: ' + params['hotspots']) if params.get('hotspots') else 'Hotspots ' + empty_state} · {('Coldspots: ' + params['coldspots']) if params.get('coldspots') else 'Coldspots ' + empty_state}")
            else:
                st.caption(source['data']['label']+" · " + str(len(params.get("selected_candidate_ids") or [])) + " selected")
        if inputs:
            st.caption("Structures, sequences, and hotspot selections follow the connections. Click their input block to edit them.")


def _render_field(field, values, key):
    name = field["name"]
    value = values.get(name, field["default"])
    label = {"num_designs":"Number of designs", "timesteps":"Diffusion timesteps", "seed":"Random seed (0 = random)"}.get(name, field["label"])
    kwargs = {"key": key+"-"+name, "help":field.get("help")}
    if name == "scaffold_ids" and field["type"] == "json":
        from mn_molecular_design_studio.peptide_trap_input import SCAFFOLDS
        return st.multiselect(label, list(SCAFFOLDS), default=list(value or []),
                              format_func=lambda item: SCAFFOLDS[item]["label"], **kwargs)
    if field["type"] == "json":
        text = st.text_area(label, value=json.dumps(value, indent=2), **kwargs)
        try:
            parsed = json.loads(text)
            if not isinstance(parsed, list):
                raise ValueError("Expected a JSON list.")
            return parsed
        except ValueError as exc:
            st.error(str(exc))
            return text  # graph validation prevents saving an invalid list
    if field["type"] == "boolean":
        return st.checkbox(label, value=bool(value), **kwargs)
    if field["type"] == "select":
        options = field["options"]
        labels = {"protein_mpnn":"ProteinMPNN", "ligand_mpnn":"LigandMPNN", "soluble_mpnn":"SolubleMPNN",
                  "molscrub_0_1_0":"MolScrub 0.1.0 defaults · fastest benchmark",
                  "molscrub_0_3_0_defaults":"MolScrub 0.3.0 defaults · 3 conformers, 400 FF iterations",
                  "molscrub_0_3_0_matched":"MolScrub 0.3.0 matched · 1 conformer, 200 FF iterations",
                  "classic":"Structure-based (classic)", "ligand_based":"Ligand-based", "hybrid":"Hybrid (receptor + reference ligand)",
                  "fast":"Fast", "balance":"Balanced", "detail":"Detailed",
                  "smiles":"SMILES", "ccd":"CCD component code",
                  "single_sequence":"Single sequence (no MSA)", "reference_then_colabfold":"Reference MSA; ColabFold if missing",
                  "protein_protein":"Protein / protein", "protein_ligand":"Protein / ligand",
                  "protein_dna":"Protein / DNA", "protein_rna":"Protein / RNA",
                  "protein_macrocyclic":"Protein / macrocyclic peptide",
                  "protein":"Protein monomer",
                  "miniprotein":"Miniprotein", "scfv":"Antibody (scFv)", "vhh":"Nanobody (VHH)",
                  "trastuzumab_framework_vhvl":"Trastuzumab framework", "atezolizumab_framework_vhvl":"Atezolizumab framework",
                  "ocankitug_framework_vhvl":"Ocankitug framework", "caplacizumab_vhh":"Caplacizumab VHH",
                  "ozoralizumab_tnf_vhh":"Ozoralizumab TNF VHH", "vobarilizumab_il6r_vhh":"Vobarilizumab IL-6R VHH",
                  "custom":"Custom masked template", "dna_rna":"DNA / RNA", "ion":"Ion", "ligand":"Ligand", "lipid":"Lipid",
                  "protein":"Protein–protein", "interface":"Interface", "epitope":"Epitope", "idp":"Disordered region",
                  "alphafold3":"Official AlphaFold3", "openfold3":"OpenFold3 P2", "openbind0":"OpenBind-0",
                  "boltz2":"Boltz-2 weights", "protenix2":"Protenix V2", "rosettafold3":"RosettaFold3",
                  "intellifold2":"IntelliFold2", "opendde":"OpenDDE"}
        return st.selectbox(label, options, index=options.index(value) if value in options else 0,
                            format_func=lambda v:labels.get(v,v), **kwargs)
    if field["type"] == "number":
        integer = isinstance(field["default"], int)
        convert = int if integer else float
        return st.number_input(label, value=convert(value), step=1 if integer else .1,
                               min_value=convert(field["min"]) if "min" in field else None,
                               max_value=convert(field["max"]) if "max" in field else None,
                               disabled=name=="cfg_scale" and not values.get("classifier_free_guidance"), **kwargs)
    if name in {"target_sequence", "target_protein_sequences", "partner_protein_sequences", "dna_sequences", "rna_sequences", "ligand_smiles", "options_json", "condition_spec"}:
        height = 170 if name == "options_json" else 220 if name == "condition_spec" else 110 if "sequence" in name else 85
        return st.text_area(label, value=str(value), height=height, **kwargs)
    return st.text_input(label, value=str(value), **kwargs)


def _centroid_manifest_row_count(path):
    """Return the available centroid rows for the active-learning seed input."""
    if not path:
        return None
    try:
        manifest = json.loads(Path(str(path)).read_text(encoding="utf-8"))
        count = int(manifest.get("row_count"))
    except (OSError, ValueError, TypeError, AttributeError):
        return None
    return count if count >= 0 else None


def _available_docking_machines(gpu_engine):
    """Show enabled registered machines with a recent usable resource snapshot.

    A stale peer heartbeat is not proof that a remote computer is down. The
    active-learning controller rechecks its scheduler and free GPUs over SSH
    when the node is applied, before it can be run.
    """
    from mn_molecular_design_studio.cluster.store import Store
    from mn_molecular_design_studio.cluster import config as cluster_config
    from mn_molecular_design_studio.cluster.config import display_name

    machines = []
    for row in Store().nodes():
        if row.get("disabled"):
            continue
        snapshot = row.get("snapshot") or {}
        resources = snapshot.get("resources") or {}
        if not snapshot.get("accept_jobs") or not snapshot.get("scheduler_online"):
            continue
        gpu_list = resources.get("free_gpus") if "free_gpus" in resources else resources.get("gpus", [])
        gpu_list = gpu_list or []
        if gpu_engine and not gpu_list:
            continue
        machines.append({"id": row["id"], "name": display_name(row["name"]), "gpus": gpu_list,
                         "cpu_free": int(resources.get("cpu_free") or 0),
                         "local": row["id"] == cluster_config.identity(),
                         "heartbeat_stale": row.get("availability") != "online"})
    return machines


def render_node_form(node, graph, *, key, render_input, owner_id):
    kind = node["data"]["kind"]
    spec = NODE_CATALOG[kind]
    values = {f["name"]:copy.deepcopy(node["data"]["params"].get(f["name"],f["default"])) for f in spec["fields"]}
    values.update({k:copy.deepcopy(v) for k,v in node["data"]["params"].items() if k in {"selected_candidate_ids", "configuration_migration", "prediction_config_version"}})
    if kind == "boltz2":
        values["prediction_config_version"] = 2
    st.html('''<style>
      [role="dialog"]:has(.st-key-pipeline-node-form){width:min(960px,calc(100vw - 48px));max-width:960px}
      .st-key-pipeline-node-form {color:#193e30}
      .st-key-pipeline-node-form .node-form-heading{text-align:center;padding:.1rem 0 1rem}
      .st-key-pipeline-node-form .node-form-heading h2{margin:.3rem 0;color:#173d31;font-size:1.9rem;letter-spacing:-.035em}
      .st-key-pipeline-node-form .node-form-heading p{margin:.5rem auto;color:#718478;font-size:.9rem;max-width:640px}
      .st-key-pipeline-node-form .node-form-eyebrow{color:#688974;text-transform:uppercase;font-size:.68rem;letter-spacing:.12em}
      .st-key-pipeline-node-form [data-testid="stVerticalBlockBorderWrapper"]{border-color:#d7e6dc!important;border-radius:12px!important}
      .st-key-pipeline-node-form [data-baseweb="input"],.st-key-pipeline-node-form [data-baseweb="select"]>div{background:#f3f8f5!important;border-color:#d7e6dc!important}
      .st-key-pipeline-node-form [data-testid="stExpander"]{border-color:#d7e6dc!important;border-radius:8px!important;background:#fcfefc}
      .st-key-pipeline-node-form [data-testid="stButton"] button{border-radius:7px!important;border-color:#bad2c3;color:#27674e}
      .st-key-pipeline-node-form button[kind="primary"]{background:#27674e!important;color:white!important;border-color:#27674e!important}
      .st-key-pipeline-node-form [data-baseweb="tag"]{background:#e0eee5!important;color:#27674e!important}
    </style>''')
    with st.container(key="pipeline-node-form"):
        st.html(f'<div class="node-form-heading"><div class="node-form-eyebrow">{escape(spec["category"])} · Pipeline step</div>'
                f'<h2>{escape(node["data"].get("label") or spec["label"])}</h2><p>{escape(node["data"].get("description") or spec["description"])}</p></div>')
        replacements = compatible_node_replacements(node, graph)
        replacement_kind = None
        if replacements:
            with st.expander("Replace this step", expanded=False):
                st.caption("Choose another tool in this category. The node stays in place, compatible connections are kept, and shared settings carry over when valid.")
                replacement_kind = st.selectbox(
                    "Replacement tool",
                    [""] + replacements,
                    format_func=lambda value: "Choose a tool…" if not value else NODE_CATALOG[value]["label"],
                    key=key + "-replacement-kind",
                ) or None
                if st.button("Replace step", key=key + "-replace", disabled=not replacement_kind,
                             help="Replace this node while retaining compatible connections."):
                    action = "replace"
                else:
                    action = None
        else:
            action = None
        coldspot_support = spec.get("coldspot_support")
        hotspot_support = spec.get("hotspot_support")
        if hotspot_support == "native":
            st.caption("Connected Target or Binder Design Context hotspots are passed to this engine’s hotspot input.")
        elif hotspot_support == "conditional":
            st.caption("Hotspots are used for the protein-binder task; ligand-design tasks use their ligand or motif constraints instead.")
        elif hotspot_support == "unsupported":
            st.caption("This step does not accept binder hotspots. Connected interface annotations are ignored here.")
        if coldspot_support == "native":
            st.caption("Connected Target or Binder Design Context coldspots are passed to this engine’s native coldspot constraint.")
        elif coldspot_support == "mapped":
            st.caption("Connected Target or Binder Design Context coldspots are passed as BoltzGen non-binding residues and mapped to its chain-local numbering.")
        elif coldspot_support == "unsupported":
            st.caption("This engine does not accept coldspots. Connected Target or Binder Design Context coldspots are ignored for this step.")
        label = st.text_input("Step name", value=node["data"]["label"], key=key+"-label").strip() or spec["label"]
        ready = True
        if kind == "binder_context":
            _render_connection_card(node, graph)
            values["clear_unselected"] = st.checkbox(
                "Clear inherited annotations when overrides are empty",
                value=bool(values.get("clear_unselected", False)),
                key=key + "-clear-unselected",
                help="Use this when the branch should receive no hotspots or coldspots, even if the connected Target Structure has them set.",
            )
            with st.container(border=True):
                contexts = render_input()
                ready = node["id"] in contexts
                if ready:
                    values.update(contexts[node["id"]])
        elif kind == "active_library":
            with st.container(border=True):
                contexts = render_input()
                ready = node["id"] in contexts
                if ready:
                    values.update(contexts[node["id"]])
        elif kind == "active_seed":
            with st.container(border=True):
                contexts = render_input()
                ready = node["id"] in contexts
                if ready:
                    values.update(contexts[node["id"]])
            for field_spec in spec["fields"]:
                if field_spec["name"] in {"sample_size", "sampling_seed"}:
                    if field_spec["name"] == "sample_size":
                        centroid_count = _centroid_manifest_row_count(values.get("centroid_manifest"))
                        if centroid_count is not None:
                            if centroid_count < 1:
                                st.error("The selected centroid dataset contains no compounds.")
                                ready = False
                                continue
                            maximum = min(int(field_spec.get("max", centroid_count)), centroid_count)
                            current = int(values.get("sample_size") or field_spec["default"])
                            values["sample_size"] = min(max(current, int(field_spec["min"])), maximum)
                            st.caption(f"Selected centroid set: {centroid_count:,} compounds. The initial docking sample cannot exceed this count.")
                            bounded_field = {**field_spec, "max": maximum}
                            values[field_spec["name"]] = _render_field(bounded_field, values, key)
                            continue
                    values[field_spec["name"]] = _render_field(field_spec, values, key)
        elif kind == "active_target":
            with st.container(border=True):
                contexts = render_input()
                ready = node["id"] in contexts
                if ready:
                    values.update(contexts[node["id"]])
            align_field = next(field for field in spec["fields"] if field["name"] == "align_longest_axis_to_x")
            values["align_longest_axis_to_x"] = _render_field(align_field, values, key)
            alignment_field = next(field for field in spec["fields"] if field["name"] == "alignment_ligand_sdf")
            values["alignment_ligand_sdf"] = _render_field(alignment_field, values, key)
            reference_field = next(field for field in spec["fields"] if field["name"] == "reference_ligand_sdf")
            if not values.get("reference_ligand_sdf"):
                values["reference_ligand_sdf"] = values.get("alignment_ligand_sdf") or ""
            values["reference_ligand_sdf"] = _render_field(reference_field, values, key)
            st.markdown("**Docking box**")
            columns = st.columns(3)
            for index, name in enumerate(("center_x", "center_y", "center_z", "size_x", "size_y", "size_z")):
                field_spec = next(field for field in spec["fields"] if field["name"] == name)
                with columns[index % 3]:
                    values[name] = _render_field(field_spec, values, key)
        elif kind == "active_docking":
            st.info("Intermediate compounds and scores stay in external job files. Database import is disabled for every active-learning batch on every computer. Set the optional final shortlist in Final results.")
            _render_connection_card(node, graph)
            engine_field = next(field for field in spec["fields"] if field["name"] == "engine")
            previous_engine = values.get("engine", "UniDock Pro")
            values["engine"] = _render_field(engine_field, values, key)
            if values["engine"] != previous_engine:
                profile = ("total_score", "REU") if values["engine"] == "RosettaLigand" else ("binding_energy", "kcal/mol") if values["engine"] == "AutoDock-GPU" else ("vina_score", "kcal/mol")
                values.update(score_field=profile[0], score_units=profile[1], score_direction="Lower is better")
            gpu_engine = values["engine"] in {"UniDock Pro", "GNINA", "AutoDock-GPU"}
            st.caption("Each selected GPU processes one batch at a time through its Studio scheduler. The computer's CPU budget is shared between its GPUs. Idle GPUs can take pending batches. Failed compounds are recorded once without rescue attempts; the organizer collects scores and retained poses before round cleanup.")
            machines = []
            try:
                machines = _available_docking_machines(gpu_engine)
            except (OSError, ValueError, KeyError, TypeError):
                machines = []
            if not machines:
                st.warning("No registered Studio computer with compatible resources is available for this docking engine.")
                ready = False
            else:
                labels = {machine["id"]: f"{machine['name']} · {str(len(machine['gpus'])) + ' GPU(s)' if gpu_engine else 'CPU'}"
                          + (" · stale heartbeat; rechecked on Apply" if machine["heartbeat_stale"] else "")
                          for machine in machines}
                available_ids = list(labels)
                current_workers = values.get("worker_nodes") if isinstance(values.get("worker_nodes"), list) else []
                current_ids = [str(item.get("node_id")) for item in current_workers if isinstance(item, dict)]
                default_ids = [item for item in current_ids if item in available_ids] or available_ids
                chosen_ids = st.multiselect("Docking computers", available_ids,
                    default=default_ids, format_func=lambda value: labels[value], key=key + "-docking-machines-" + values["engine"])
                worker_specs = []
                no_cpu_slots = False
                for machine_id in chosen_ids:
                    machine = next(item for item in machines if item["id"] == machine_id)
                    gpu_ids = [str(gpu.get("gpu_id", index)) for index, gpu in enumerate(machine["gpus"])]
                    old = next((item for item in current_workers if isinstance(item, dict) and item.get("node_id") == machine_id), {})
                    cpu_available = int(machine["cpu_free"])
                    if machine["local"]:
                        # The campaign coordinator holds one scheduled CPU slot while it queues the docking shard.
                        cpu_available = max(0, cpu_available - 1)
                    if cpu_available < 1:
                        no_cpu_slots = True
                        continue
                    cpu_default = min(cpu_available, max(1, int(old.get("cpu_cores") or cpu_available)))
                    selected_cpus = st.number_input(
                        f"Total CPU worker slots · {machine['name']}", min_value=1, max_value=cpu_available,
                        value=cpu_default, step=1, key=key + "-docking-cpus-" + machine_id,
                    )
                    selected_gpus = []
                    if gpu_engine:
                        gpu_values = [str(value) for value in (old.get("gpu_devices") or [old.get("gpu_device") or gpu_ids[0]])]
                        selected_gpus = st.multiselect(f"GPUs · {machine['name']}", gpu_ids,
                            default=[value for value in gpu_values if value in gpu_ids],
                            key=key + "-docking-gpus-" + machine_id,
                            format_func=lambda value: f"GPU {value}")
                        if not selected_gpus or len(selected_gpus) > int(selected_cpus):
                            ready = False
                            st.warning("Select at least one GPU and allocate at least one CPU slot per GPU.")
                    batch_share = int(st.number_input(
                        f"Batch share per GPU · {machine['name']}", min_value=1, max_value=100,
                        value=min(100, max(1, int(old.get("batch_share") or 1))), step=1,
                        key=key + "-docking-share-" + machine_id,
                        help="Relative initial share per GPU. Faster GPUs can take batches still pending on other workers, so the final allocation adjusts during the round.",
                    ))
                    worker_specs.append({"node_id": machine_id, "gpu_devices": selected_gpus,
                                         "gpu_device": selected_gpus[0] if selected_gpus else "none",
                                         "cpu_cores": int(selected_cpus), "name": machine["name"],
                                         "local": machine["local"], "batch_share": batch_share})
                values["worker_nodes"] = worker_specs
                fields = [field for field in spec["fields"] if field["name"] not in {"worker_nodes", "engine"}
                          and (field["name"] not in {"search_mode", "docking_mode"} or values["engine"] == "UniDock Pro")]
                with st.container(border=True):
                    st.markdown("**Docking and batch settings**")
                    columns = st.columns(2)
                    for index, field_spec in enumerate(fields):
                        with columns[index % 2]:
                            field_key = key + "-" + values["engine"] if field_spec["name"] in {"score_field", "score_units", "score_direction"} else key
                            values[field_spec["name"]] = _render_field(field_spec, values, field_key)
                if no_cpu_slots:
                    st.warning("A selected computer has no free CPU slots for a docking worker. Free some slots or choose another computer.")
                    ready = False
                if not worker_specs:
                    st.info("Select at least one computer with free CPU slots. The organizer must be included.")
                    ready = False
        elif kind in {"target", "design_source", "generation_source", "backbone_source", "project_design_source"}:
            with st.container(border=True):
                contexts = render_input()
                ready = node["id"] in contexts
                if ready:
                    values.update(contexts[node["id"]])
            if kind in {"design_source","generation_source","backbone_source","project_design_source"}:
                ready = ready and bool(values.get("selected_candidate_ids"))
        elif kind == "cofold_input":
            values["task_type"] = _render_field(next(f for f in spec["fields"] if f["name"] == "task_type"), values, key)
            task_type = values["task_type"]
            st.caption("Provide sequence inputs for this task. A target PDB path is optional for most models and required by AlphaFold2 Initial Guess and BoltzGen Fold.")
            st.caption("Choose compatible prediction nodes for this molecule type and ligand notation. Task setup defaults to one engine; comparisons are explicit.")
            visible = {"target_protein_sequences", "target_pdb", "target_chain_id", "sequence_file"}
            if task_type in {"protein_protein", "protein_macrocyclic"}:
                visible.add("partner_protein_sequences")
            elif task_type == "protein_ligand":
                visible.update({"ligand_input_type", "ligand_ccd_codes", "ligand_smiles", "compound_set_file"})
            elif task_type == "protein_dna":
                visible.add("dna_sequences")
            elif task_type == "protein_rna":
                visible.add("rna_sequences")
            for field in spec["fields"]:
                if field["name"] in visible and field["name"] != "task_type":
                    values[field["name"]] = _render_field(field, values, key)
            target = bool(values.get("target_protein_sequences", "").strip() or values.get("sequence_file"))
            partner = bool(values.get("partner_protein_sequences", "").strip())
            if task_type == "protein":
                ready = target
            elif task_type == "protein_protein":
                ready = target and partner
            elif task_type == "protein_ligand":
                if values.get("ligand_input_type") == "ccd":
                    values["ligand_smiles"] = ""
                    ready = target and bool(values.get("compound_set_file") or values.get("ligand_ccd_codes", "").strip())
                else:
                    values["ligand_ccd_codes"] = ""
                    ready = target and bool(values.get("compound_set_file") or values.get("ligand_smiles", "").strip())
            elif task_type == "protein_dna":
                ready = target and bool(values.get("dna_sequences", "").strip())
            elif task_type == "protein_rna":
                ready = target and bool(values.get("rna_sequences", "").strip())
            elif task_type == "protein_macrocyclic":
                ready = target and partner
                values["cyclic_partner"] = True
            if not ready:
                st.info("Enter the sequence or molecule inputs required for the selected prediction task.")
        elif kind in {"esm_if1", "caliby", "boltzgen_inverse"}:
            from mn_molecular_design_studio.inverse_folding_ui import render_inputs
            _render_connection_card(node, graph)
            with st.container(key="inverse-node"):
                updated, error = render_inputs(kind, owner_id=owner_id, key_prefix=key + "-inverse", initial=values, connected=True)
            values.update({f["name"]: updated[f["name"]] for f in spec["fields"] if f["name"] in updated})
            ready = not error
            if error:
                st.caption(error)
        elif kind in {"ipsae", "rosetta_metrics"}:
            _render_connection_card(node, graph)
            if kind == "ipsae":
                st.caption("Reads the predicted complex and its saved pairwise PAE matrix. Boltz-2 full PAE output is enabled in the all-generators starter template.")
            else:
                st.caption("Adds PyRosetta interface energy, shape complementarity, buried surface area, and related metrics to each connected predicted complex.")
        elif kind == "surface_analysis":
            _render_connection_card(node, graph)
            st.caption("Structure and selected chains only. Surface inspection ignores binder hotspots and coldspots; those annotations remain available to connected design branches.")
            _, upstream = _connected_inputs(graph, node["id"])
            target_node = next((source for source in upstream if source["data"]["kind"] == "target"), None)
            target_chains = target_node["data"]["params"].get("target_chains", "") if target_node else ""
            target_chains = [chain for chain in target_chains.replace(";", ",").replace(" ", ",").split(",") if chain] if isinstance(target_chains, str) else list(target_chains or [])
            if not target_chains:
                st.warning("Connect a Target block and select at least one protein chain.")
                ready = False
            else:
                if values.get("masif_chain_id") not in target_chains:
                    values["masif_chain_id"] = target_chains[0]
                values["masif_chain_id"] = st.selectbox(
                    "MaSIF-seed target chain",
                    target_chains,
                    index=target_chains.index(values["masif_chain_id"]),
                    key=key + "-masif-chain",
                )
            st.caption("PeSTo and ScanNet keep their individual prediction modes. Surf2Spot and MaSIF-seed run on the same selected target; MaSIF-seed uses the chain chosen above.")
            with st.container(border=True):
                st.markdown("**PeSTo and ScanNet modes**")
                columns = st.columns(2)
                for index, field in enumerate(f for f in spec["fields"] if f["name"] != "masif_chain_id"):
                    with columns[index % 2]:
                        values[field["name"]] = _render_field(field, values, key)
        elif kind == "esm_binder":
            _render_connection_card(node, graph)
            st.caption("Uses the selected target chain’s sequence; no MSA is required. The output contains designed sequences and ESMFold2-predicted complexes for downstream ProteinMPNN or structure-prediction steps.")
            _, upstream = _connected_inputs(graph, node["id"])
            target_node = next((source for source in upstream if source["data"]["kind"] == "target"), None)
            if target_node:
                target_chains = target_node["data"]["params"].get("target_chains", "")
                target_chains = [chain for chain in target_chains.replace(";", ",").replace(" ", ",").split(",") if chain] if isinstance(target_chains, str) else list(target_chains or [])
                if len(target_chains) != 1:
                    st.warning("ESMFold2 Binder Design requires exactly one selected target chain. Select one chain in the connected Target block.")
                    ready = False
            from mn_molecular_design_studio.esmfold2_binder_ui import FRAMEWORKS, VHH_FRAMEWORK_OPTIONS, _template_from_text
            modalities = {"Miniprotein": "miniprotein", "Antibody (scFv)": "scfv", "Nanobody (VHH)": "vhh"}
            selected_modality = st.pills("Binder type", list(modalities), selection_mode="single",
                default=next((label for label, value in modalities.items() if value == values.get("binder_type")), "Miniprotein"),
                key=key + "-esm-binder-type") or "Miniprotein"
            values["binder_type"] = modalities[selected_modality]
            values["num_designs"] = st.number_input("Number of designs", min_value=1, max_value=128,
                value=int(values.get("num_designs", 1)), step=1, key=key + "-esm-num-designs")
            if values["binder_type"] == "miniprotein":
                values["binder_length"] = st.number_input("Binder length", min_value=40, max_value=200,
                    value=int(values.get("binder_length", 80)), step=1, key=key + "-esm-binder-length")
            elif values["binder_type"] == "scfv":
                choices = {**FRAMEWORKS}
                current = values.get("binder_framework")
                default_label = next((label for label, framework in choices.items() if framework == current), next(iter(choices)))
                framework_label = st.selectbox("Antibody framework", list(choices), index=list(choices).index(default_label), key=key + "-esm-scfv-framework")
                values["binder_framework"] = choices[framework_label]
            else:
                choices = {**VHH_FRAMEWORK_OPTIONS, "Custom masked template": "custom"}
                current = values.get("binder_framework")
                default_label = next((label for label, framework in choices.items() if framework == current), next(iter(choices)))
                framework_label = st.selectbox("VHH framework", list(choices), index=list(choices).index(default_label), key=key + "-esm-vhh-framework")
                values["binder_framework"] = choices[framework_label]
                if values["binder_framework"] == "custom":
                    values["binder_template"] = st.text_area("VHH framework template", value=str(values.get("binder_template") or ""),
                        height=110, placeholder="Paste a VHH sequence and replace its three CDRs with #…", key=key + "-esm-vhh-template",
                        help="Mark exactly three CDR runs with # characters; keep framework residues as amino-acid letters.")
                    try:
                        values["binder_template"] = _template_from_text(values["binder_template"])
                    except ValueError as exc:
                        st.error(str(exc))
                        ready = False
            with st.expander("Optional settings", expanded=False):
                columns = st.columns(2)
                optional = [f for f in spec["fields"] if f["name"] in {"optimization_steps", "seed", "learning_rate", "compile_model", "checkpoint_lm"}]
                for index, field in enumerate(optional):
                    with columns[index % 2]:
                        values[field["name"]] = _render_field(field, values, key)
        elif kind in {"ligandmpnn", "foundry_mpnn"}:
            from mn_molecular_design_studio.mpnn_ui import render_inputs
            _render_connection_card(node, graph)
            initial = dict(values)
            initial["num_seq_per_target"] = int(initial.pop("sequences")) * int(initial.pop("number_of_batches", 1))
            initial["noise_level"] = float(initial["noise_level"]) if initial.get("noise_level") not in ("", None) else None
            with st.container(key="mpnn-node"):
                updated, error = render_inputs(owner_id=owner_id, key_prefix=key + "-mpnn", initial=initial, connected=True)
            values["sequences"] = updated["num_seq_per_target"]
            if kind == "foundry_mpnn":
                values["number_of_batches"] = 1
            for field in spec["fields"]:
                name = field["name"]
                if name not in updated:
                    continue
                value = updated[name]
                if isinstance(value, dict):
                    value = json.dumps(value) if value else ""
                if name == "noise_level":
                    value = "" if value is None else str(value)
                if name == "seed" and value is None:
                    value = 0
                values[name] = value
            ready = not error
            if error:
                st.caption(error)
        elif kind in {'pxdesign','rfdiffusion','genie3','protpardelle'}:
            from mn_molecular_design_studio.generator_ui import render_inputs
            parents,upstream=_connected_inputs(graph,node['id'])
            connected=next((n['data']['params'] for n in upstream if n['data']['kind']=='target'),None)
            if parents:_render_connection_card(node,graph)
            updated,_,error=render_inputs(kind,owner_id=owner_id,key_prefix='generator-node-'+key,
                                          initial=values,connected_target=connected,include_jobs=False)
            updated['target_chains']=','.join(updated['target_chains'])
            allowed={f['name'] for f in spec['fields']}
            values.update({k:('' if v is None else v) for k,v in updated.items() if k in allowed})
            ready=not error
            if error:st.caption(error)
        elif kind=='complexa_generator':
            from mn_molecular_design_studio.complexa_ui import render_inputs
            parents,upstream=_connected_inputs(graph,node['id'])
            connected=next((n['data']['params'] for n in upstream if n['data']['kind']=='target'),None)
            if parents:_render_connection_card(node,graph)
            # This is the same structure, task and sampling form as the tile.
            updated,_,error=render_inputs(owner_id=owner_id,key_prefix=key+'-complexa',generation_only=True,
                                          initial=values,connected_target=connected,include_jobs=False)
            updated['target_chains']=','.join(updated['target_chains'])
            values.update(updated)
            ready=not error
            if error:st.caption(error)
        else:
            _render_connection_card(node, graph)
            if kind == "rfd3":
                st.caption("Protein binder design · the connected Target block provides structure, chains and hotspots.")
            primary = [f for f in spec["fields"] if f.get("section")=="basic" or spec["category"]=="Analysis"]
            optional = [f for f in spec["fields"] if f not in primary]
            if kind == "bindcraft2_generator":
                primary = sorted(primary, key=lambda item: item["name"] != "modality")
            if primary:
                with st.container(border=True):
                    st.markdown("**Design settings**" if spec["category"].startswith("Design") else "**Prediction settings**" if spec["category"]=="Structure prediction" else "**Step settings**")
                    fields = primary
                    if kind == "rfd3":
                        lengths = str(values["binder_length"]).split("-",1)
                        try:
                            lo, hi = int(lengths[0]), int(lengths[-1])
                        except ValueError:
                            lo,hi=100,150
                        a,b = st.columns(2)
                        minimum = int(a.number_input("Minimum binder length",1,3000,lo,key=key+"-length-min"))
                        maximum = int(b.number_input("Maximum binder length",1,3000,hi,key=key+"-length-max"))
                        values["binder_length"] = f"{minimum}-{maximum}"
                        if maximum < minimum:
                            st.error("Maximum binder length must be at least the minimum.")
                            ready = False
                        fields = [f for f in primary if f["name"]!="binder_length"]
                    columns = st.columns(2)
                    for i,field in enumerate(fields):
                        if kind == "bindcraft2_generator" and field["name"] == "binder_length":
                            if values.get("modality") == "VHH":
                                values["binder_length"] = ""
                                continue
                            if not values.get("binder_length"):
                                values["binder_length"] = "60-100"
                        if kind == "bindcraft2_generator" and field["name"] == "scaffold_edits":
                            if values.get("modality") != "VHH":
                                values["scaffold_edits"] = ""
                                continue
                            if not values.get("scaffold_edits"):
                                from mn_molecular_design_studio.engine.workflow_engine.workflows.bindcraft2 import BC2_SCAFFOLD_EDIT_DEFAULTS
                                values["scaffold_edits"] = BC2_SCAFFOLD_EDIT_DEFAULTS["VHH"]
                        with columns[i%2]:
                            if kind == "bindcraft2_generator" and field["name"] == "scaffold_edits":
                                values[field["name"]] = st.text_area(field["label"], value=str(values[field["name"]]), height=95,
                                    key=key+"-"+field["name"], help=field.get("help"))
                            else:
                                values[field["name"]] = _render_field(field,values,key)
                    if kind == "bindcraft2_generator" and values.get("modality") == "VHH":
                        st.caption("VHH scaffold mode fixes the framework; the edit ranges above define the CDR and loop regions to design.")
            if optional:
                with st.expander("Optional settings", expanded=False):
                    columns = st.columns(2)
                    for i,field in enumerate(optional):
                        with columns[i%2]:
                            values[field["name"]] = _render_field(field,values,key)
            if kind == "rfd3" and values.get("low_memory_mode"):
                st.warning("Low memory mode uses a chunked path with partially missing checkpoint weights and may reduce prediction quality.")
            if kind == "boltz2":
                for notice in values.get("configuration_migration", []):
                    st.warning(notice)
                st.caption("Samples and recycles are the effective Boltz settings. The installed engine always writes full PAE, retained for ipSAE; confidence scores are not experimental binding probabilities.")
            if kind == "af2":
                from mn_molecular_design_studio.prediction_configuration import AF2_MODES
                st.caption(AF2_MODES[values.get("af2_validation_mode", "conditioned_self_consistency")])
            if kind == "results":
                st.info("Collect candidate structures and metrics from the connected branches. Existing results remain available when you edit the pipeline.")
            if kind == "filter":
                st.caption("Choose native fields, common fields, or evidence/<stage>/<tool>/<metric>. Add a candidate ID before the metric to resolve repeated stages. Missing, failed and nonfinite required measurements cannot pass. Ranking and limits apply within recorded comparable protocols. All decisions and selected sequence/structure exports remain in this step's selection files.")
        st.caption("Apply saves this step to the draft. The pipeline runs only when you choose Run pipeline on the canvas.")
        a,b,c = st.columns([1.3,1,1.3])
        if a.button("Reset to tool defaults", key=key+"-defaults", disabled=kind in {"target","binder_context","design_source","generation_source","backbone_source","project_design_source","results"}, use_container_width=True):
            action="reset"
        if b.button("Cancel", key=key+"-cancel", use_container_width=True):
            action="cancel"
        if c.button("Apply to pipeline", type="primary", key=key+"-apply", disabled=not ready, use_container_width=True):
            action="apply"
    return values,label,action,replacement_kind
