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


def _connected_inputs(graph, node_id):
    nodes = {n["id"]: n for n in graph["nodes"]}
    parents = [nodes[e["source"]] for e in graph["edges"] if e["target"] == node_id and e["source"] in nodes]
    seen, upstream, queue = set(), [], list(parents)
    while queue:
        node = queue.pop(0)
        if node["id"] in seen:
            continue
        seen.add(node["id"])
        if node["data"]["kind"] in {"target", "design_source", "generation_source", "backbone_source"}:
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
                b.caption("CHAINS · HOTSPOTS")
                b.write(f"{params.get('target_chains') or '—'} · {params.get('hotspots') or 'No hotspots selected'}")
            else:
                st.caption(source['data']['label']+" · " + str(len(params.get("selected_candidate_ids") or [])) + " selected")
        if inputs:
            st.caption("Structures, sequences, and hotspot selections follow the connections. Click their input block to edit them.")


def _render_field(field, values, key):
    name = field["name"]
    value = values.get(name, field["default"])
    label = {"num_designs":"Number of designs", "timesteps":"Diffusion timesteps", "seed":"Random seed (0 = random)"}.get(name, field["label"])
    kwargs = {"key": key+"-"+name, "help":field.get("help")}
    if field["type"] == "boolean":
        return st.checkbox(label, value=bool(value), **kwargs)
    if field["type"] == "select":
        options = field["options"]
        labels = {"protein_mpnn":"ProteinMPNN", "ligand_mpnn":"LigandMPNN", "soluble_mpnn":"SolubleMPNN",
                  "single_sequence":"Single sequence (no MSA)", "reference_then_colabfold":"Reference MSA; ColabFold if missing",
                  "miniprotein":"Miniprotein", "scfv":"Antibody (scFv)", "vhh":"Nanobody (VHH)",
                  "trastuzumab_framework_vhvl":"Trastuzumab framework", "atezolizumab_framework_vhvl":"Atezolizumab framework",
                  "ocankitug_framework_vhvl":"Ocankitug framework", "caplacizumab_vhh":"Caplacizumab VHH",
                  "ozoralizumab_tnf_vhh":"Ozoralizumab TNF VHH", "vobarilizumab_il6r_vhh":"Vobarilizumab IL-6R VHH",
                  "custom":"Custom masked template", "dna_rna":"DNA / RNA", "ion":"Ion", "ligand":"Ligand", "lipid":"Lipid",
                  "protein":"Protein–protein", "interface":"Interface", "epitope":"Epitope", "idp":"Disordered region"}
        return st.selectbox(label, options, index=options.index(value) if value in options else 0,
                            format_func=lambda v:labels.get(v,v), **kwargs)
    if field["type"] == "number":
        integer = isinstance(field["default"], int)
        convert = int if integer else float
        return st.number_input(label, value=convert(value), step=1 if integer else .1,
                               min_value=convert(field["min"]) if "min" in field else None,
                               max_value=convert(field["max"]) if "max" in field else None,
                               disabled=name=="cfg_scale" and not values.get("classifier_free_guidance"), **kwargs)
    return st.text_input(label, value=str(value), **kwargs)


def render_node_form(node, graph, *, key, render_input, owner_id):
    kind = node["data"]["kind"]
    spec = NODE_CATALOG[kind]
    values = {f["name"]:copy.deepcopy(node["data"]["params"].get(f["name"],f["default"])) for f in spec["fields"]}
    values.update({k:copy.deepcopy(v) for k,v in node["data"]["params"].items() if k=="selected_candidate_ids"})
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
                f'<h2>{escape(spec["label"])}</h2><p>{escape(spec["description"])}</p></div>')
        label = st.text_input("Step name", value=node["data"]["label"], key=key+"-label").strip() or spec["label"]
        ready = True
        if kind in {"target", "design_source", "generation_source", "backbone_source"}:
            with st.container(border=True):
                contexts = render_input()
                ready = node["id"] in contexts
                if ready:
                    values.update(contexts[node["id"]])
            if kind in {"design_source","generation_source","backbone_source"}:
                ready = ready and bool(values.get("selected_candidate_ids"))
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
                        with columns[i%2]:
                            values[field["name"]] = _render_field(field,values,key)
            if optional:
                with st.expander("Optional settings", expanded=False):
                    columns = st.columns(2)
                    for i,field in enumerate(optional):
                        with columns[i%2]:
                            values[field["name"]] = _render_field(field,values,key)
            if kind == "rfd3" and values.get("low_memory_mode"):
                st.warning("Low memory mode uses a chunked path with partially missing checkpoint weights and may reduce prediction quality.")
            if kind == "af2":
                st.caption("Initial guess is enabled for this pipeline step.")
            if kind == "results":
                st.info("Collect candidate structures and metrics from the connected branches. Existing results remain available when you edit the pipeline.")
        st.caption("Apply saves this step to the draft. The pipeline runs only when you choose Run pipeline on the canvas.")
        a,b,c = st.columns([1.3,1,1.3])
        action = None
        if a.button("Reset to tool defaults", key=key+"-defaults", disabled=kind in {"target","design_source","generation_source","backbone_source","results"}, use_container_width=True):
            action="reset"
        if b.button("Cancel", key=key+"-cancel", use_container_width=True):
            action="cancel"
        if c.button("Apply to pipeline", type="primary", key=key+"-apply", disabled=not ready, use_container_width=True):
            action="apply"
    return values,label,action
