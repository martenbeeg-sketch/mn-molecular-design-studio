"""Tamarind-inspired inverse-folding forms, shared with the node editor."""
from __future__ import annotations

import copy
import hashlib
from html import escape
from pathlib import Path

import streamlit as st

from mn_molecular_design_studio import inverse_folding_input as inputs
from mn_molecular_design_studio.auth import user_input_root
from mn_molecular_design_studio.mpnn_ui import STYLE as MPNN_STYLE, _existing_results, _structure_fields

STYLE = MPNN_STYLE.replace("mpnn-surface", "inverse-surface").replace("mpnn-node", "inverse-node") + """<style>
 :is(.st-key-inverse-surface,.st-key-inverse-node) button[data-variant=segmented_control][aria-checked=true]{background:#e4f0e8!important;color:#27674e!important;border-color:#4c926c!important}
 :is(.st-key-inverse-surface,.st-key-inverse-node) button[data-variant=segmented_control][aria-checked=true] p{color:#27674e!important}
 :is(.st-key-inverse-surface,.st-key-inverse-node) [data-baseweb=checkbox] input:checked+div{background-color:#27674e!important;border-color:#27674e!important}
 </style>"""
ASSETS = Path(__file__).with_name("assets") / "inverse_folding"
MODEL_HELP = {
    "caliby": "General-purpose model trained on PDB protein chains. Use ensembles for the full Caliby protocol.",
    "soluble_caliby": "Trained without annotated transmembrane proteins for soluble-expression targets.",
    "soluble_caliby_v1": "SolubleCaliby trained on monomers and interfaces.",
    "caliby_distill": "Distilled ensemble-conditioned model for single-structure screening.",
    "soluble_caliby_distill": "Distilled soluble model for single-structure screening.",
}


def _save_upload(upload, owner, folder):
    contents = upload.getvalue()
    destination = user_input_root(owner) / "inverse-folding" / folder / (hashlib.sha256(contents).hexdigest()[:12] + "-" + Path(upload.name).name)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if not destination.exists():
        destination.write_bytes(contents)
    return destination


def _load_example(engine, prefix, owner):
    filename = "3HTN.pdb" if engine == "esm_if1" else "1CRN.pdb"
    original = ASSETS / filename
    path = user_input_root(owner) / "inverse-folding/examples" / filename
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(original.read_bytes())
    p = inputs.defaults(engine) | dict(structure_path=str(path), design_chains="A")
    if engine == "esm_if1":
        p["designed_residues"] = "A62-70"
    if engine == "caliby":
        p["omit_aas"] = "CG"
    for key in list(st.session_state):
        if key.startswith(prefix + "-field-") or key.startswith(prefix + "-source-"):
            del st.session_state[key]
    st.session_state[prefix + "-preset"] = p


def _provided_structures(owner, prefix, label, folder):
    uploads = st.file_uploader(label, type=["pdb", "ent", "cif", "mmcif", "zip"], accept_multiple_files=True,
        key=prefix + "-" + folder, help="PDB/mmCIF files, or a ZIP containing structure files.")
    paths = []
    for upload in uploads or []:
        path = _save_upload(upload, owner, folder)
        paths += inputs.expand_structure_archive(path, path.with_suffix("")) if path.suffix.lower() == ".zip" else [path]
    return [str(v) for v in paths]


def render_inputs(engine, *, owner_id, key_prefix, initial=None, connected=False, batch=False):
    p = inputs.defaults(engine) | copy.deepcopy(initial or {})
    prefix = key_prefix
    if not connected:
        p.update(st.session_state.get(prefix + "-preset", {}))
    error = None
    st.html(STYLE)
    if engine == "caliby":
        p["task"] = st.pills("Task", ["design", "score"], default=p["task"],
            format_func=lambda v: "Sequence Design" if v == "design" else "Scoring", key=prefix + "-field-task") or "design"
        st.caption("Design sequences for a structure or structural ensemble." if p["task"] == "design" else "Score the global Potts energy of the input structure's sequence.")
    if not connected:
        if batch:
            p.pop("structure_path", None)
            p["structure_paths"] = _provided_structures(owner_id, prefix, "Input structures", "batch")
            if not p["structure_paths"]:
                error = "Upload one or more structures."
            p["design_chains"] = st.text_input("Designed chains", p["design_chains"], key=prefix + "-field-batch-chains", help="Blank designs every protein chain in each structure.")
            p["designed_residues"] = st.text_input("Designed residues", p["designed_residues"], key=prefix + "-field-batch-residues", help="Blank designs all positions in the chosen chains.")
        else:
            st.markdown("**Input structure(s) to design**" if engine == "caliby" else "**PDB File**")
            p, error = _structure_fields(p, prefix, owner_id, single_chain=engine == "esm_if1", allow_non_protein=engine == "boltzgen_inverse")
    else:
        st.caption("Connected structures and binder chain roles are inherited. The target and fixed motifs remain protected.")
        p["design_chains"] = st.text_input("Designed chains (optional)", p["design_chains"], key=prefix + "-field-connected-chains", help="Blank uses the upstream candidate's designed chains.")
        p["designed_residues"] = st.text_input("Designed residues (optional)", p["designed_residues"], key=prefix + "-field-connected-residues", placeholder="e.g. B26-32 B52-56")
    if engine == "caliby":
        p["model_type"] = st.pills("Model", list(inputs.CALIBY_MODELS), default=p["model_type"], format_func=inputs.CALIBY_MODELS.get,
            key=prefix + "-field-model") or "caliby"
        st.caption(MODEL_HELP[p["model_type"]])
        distill = p["model_type"].endswith("_distill")
        options = ["single"] if distill else ["single", "provided"] if p["task"] == "score" else ["single", "generate", "provided"]
        initial_mode = p["ensemble_mode"] if p["ensemble_mode"] in options else "single"
        labels = {"single": "Single structure", "generate": "Generate Ensembles", "provided": "Provide Ensembles"}
        p["ensemble_mode"] = st.pills("Structural ensembles", options, default=initial_mode, format_func=labels.get,
            key=prefix + "-field-ensemble-" + p["task"] + str(distill)) or initial_mode
        if p["ensemble_mode"] == "provided":
            if connected or batch:
                error = "Provided ensembles need one primary structure. Use Single structure or Generate Ensembles for connected batches."
            p["ensemble_paths"] = _provided_structures(owner_id, prefix, "Conformer structures", "conformers")
        if p["ensemble_mode"] == "generate":
            p["num_conformers"] = int(st.number_input("Number of ensembles per PDB", 1, 500, p["num_conformers"], key=prefix + "-field-conformers"))
            st.caption("Protpardelle-1C generates conformers for Caliby's ensemble-conditioned sequence design.")
    if p["task"] == "design":
        p["num_sequences"] = int(st.number_input("Number of sequences per PDB" if engine == "caliby" else "Number of Sequences", 1, 10000, p["num_sequences"], key=prefix + "-field-count"))
    with st.expander("Optional settings", expanded=False):
        p["fixed_residues"] = st.text_input("Residues to fix", p["fixed_residues"], key=prefix + "-field-fixed", placeholder="e.g. A10-20 A25")
        if engine != "esm_if1" and p["task"] == "design":
            p["omit_aas"] = st.text_input("Amino Acids to Omit", p["omit_aas"], key=prefix + "-field-omit", placeholder="e.g. C,G")
        a, b = st.columns(2)
        p["temperature"] = a.number_input("Sampling temperature", .000001, 5., float(p["temperature"]), format="%.4f", key=prefix + "-field-temperature")
        p["seed"] = int(b.number_input("Random seed (0 = random)", 0, 2147483647, p["seed"], key=prefix + "-field-seed"))
        if engine == "esm_if1":
            p["multichain_backbone"] = st.checkbox("Use other chains as structural context", value=p["multichain_backbone"], key=prefix + "-field-multichain")
        if engine == "caliby":
            if p["task"] == "design":
                st.markdown("**Position constraints CSV**")
                st.caption("Use pdb_key plus fixed_pos_seq, fixed_pos_scn, fixed_pos_override_seq, pos_restrict_aatype or symmetry_pos. Positions use the input's chain IDs and author residue numbers.")
                upload = st.file_uploader("Upload position constraints", type=["csv"], key=prefix + "-constraints")
                if upload:
                    p["constraints_csv"] = str(_save_upload(upload, owner_id, "constraints"))
                p["num_sweeps"] = int(st.number_input("Potts sampling sweeps", 1, 100000, p["num_sweeps"], key=prefix + "-field-sweeps"))
            p["batch_size"] = int(st.number_input("Batch size", 1, 64, p["batch_size"], key=prefix + "-field-batch-size"))
            modes = ["off"]
            if p["task"] == "design" and p["ensemble_mode"] == "single":
                modes += ["fast"]
            if p["task"] == "design" and p["ensemble_mode"] == "generate" and p["num_conformers"] == 32:
                modes += ["exact"]
            p["optimization_mode"] = st.selectbox("Acceleration", modes, index=modes.index(p["optimization_mode"]) if p["optimization_mode"] in modes else 0,
                key=prefix + "-field-acceleration-" + p["ensemble_mode"] + p["task"] + str(p["num_conformers"]),
                format_func=lambda v: {"off": "Upstream", "fast": "Fast · small numerical differences", "exact": "Exact · same output as upstream"}[v])
    try:
        inputs.prepare(engine, p, require_input=not connected)
    except (ValueError, TypeError, OSError) as exc:
        error = error or str(exc)
    return p, error


def render_setup(engine, schema, *, owner_id, project_selector, submit, search_jobs):
    meta = inputs.ENGINES[engine]
    prefix = "inverse-" + engine
    st.html(STYLE)
    if st.button("← Back to All tools", key=prefix + "-back"):
        st.session_state["page"] = "All tools"
        st.rerun()
    with st.container(key="inverse-surface"):
        st.html(f'<div class="mpnn-header"><div class="mpnn-eyebrow">Inverse folding</div><h1>{escape(meta["title"])}</h1><p>{escape(meta["description"])}</p></div>')
        a, b, c, d = st.columns(4)
        if a.button("Show Outputs", key=prefix + "-outputs", use_container_width=True):
            st.session_state[prefix + "-show-outputs"] = not st.session_state.get(prefix + "-show-outputs")
        if b.button("Show Settings", key=prefix + "-settings", use_container_width=True):
            st.session_state[prefix + "-show-settings"] = not st.session_state.get(prefix + "-show-settings")
        c.link_button("Source code", meta["source"], use_container_width=True)
        d.link_button("Input reference", meta["reference"], use_container_width=True)
        if st.session_state.get(prefix + "-show-outputs"):
            st.info("FASTA sequences, scores CSV, and matching backbone structures with original chain and residue IDs. Structures retain input geometry; connect a structure-prediction node to evaluate each new sequence.")
        if st.session_state.get(prefix + "-show-settings"):
            st.dataframe([dict(parameter=f["name"], default=str(f.get("default", ""))) for f in schema["parameters"]], hide_index=True)
        with st.expander("About this workflow"):
            st.write(meta["description"] + " Inputs can come from any generator, structure-prediction engine or uploaded structure. Output sequences connect to the same folding and analysis nodes as ProteinMPNN.")
            st.caption(meta["version"])
        a, b = st.columns(2)
        name = a.text_input("Job name", engine.replace("_", "-") + "-design", key=prefix + "-name").strip()
        with b:
            project = project_selector(prefix)
        mode = st.segmented_control("Input source", ["Structure", "Batch structures", "Existing results"], default="Structure", key=prefix + "-input-mode") or "Structure"
        source, source_error = _existing_results(project, search_jobs, prefix) if mode == "Existing results" else ({}, None)
        if mode == "Structure":
            st.button("Load example inputs", key=prefix + "-example", on_click=_load_example, args=(engine, prefix, owner_id))
        p, error = render_inputs(engine, owner_id=owner_id, key_prefix=prefix, connected=mode == "Existing results", batch=mode == "Batch structures")
        p.update(source)
        p["job_name"] = name
        ready = inputs.readiness(engine, p)
        with st.expander("Scheduler & reference files"):
            p["gpu_device"] = st.selectbox("GPU assignment", ["all", "0", "1"], format_func=lambda v: {"all": "Automatic", "0": "GPU 0", "1": "GPU 1"}[v], key=prefix + "-gpu")
            priority = st.checkbox("Prioritize in scheduler queue", key=prefix + "-priority")
            if ready["ready"]:
                st.success("Selected model parameters are available in shared reference files.")
            else:
                st.warning("Missing reference files: " + ", ".join(ready["missing"]))
        error = source_error or error
        if error:
            st.caption(error)
        if not project:
            st.caption("Select a project to submit.")
        with st.columns([1, 1, 1])[1]:
            clicked = st.button("Submit", type="primary", key=prefix + "-submit", use_container_width=True, disabled=not project or bool(error) or not ready["ready"])
        if clicked:
            allowed = {f["name"] for f in schema["parameters"]}
            submit(schema, inputs.WORKFLOWS[engine], {k: v for k, v in p.items() if k in allowed}, project, {"prioritize": priority}, name or meta["title"])
