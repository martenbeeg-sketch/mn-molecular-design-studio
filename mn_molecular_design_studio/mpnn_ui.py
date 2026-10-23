"""ProteinMPNN tool page and the identical MPNN settings on pipeline nodes."""
from __future__ import annotations

import copy
import hashlib
import json
from html import escape
from pathlib import Path

import streamlit as st

from mn_molecular_design_studio import mpnn_input as inputs
from mn_molecular_design_studio.auth import user_input_root
from mn_molecular_design_studio.target_editor import TargetSource, render_target_source
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json
from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer import StructureVisualization, molstar_custom_component

STYLE = """<style>
 .st-key-mpnn-surface{max-width:1000px;margin:auto;color:#233e31;--mpnn-green:#27674e}
 .mpnn-header{text-align:center;margin:1rem auto 1.3rem}.mpnn-header h1{font-size:2rem;letter-spacing:-.035em;color:#183d30;margin-bottom:.3rem}
 .mpnn-header p{font-size:.9rem;color:#819087}.mpnn-eyebrow{font-size:.68rem;letter-spacing:.13em;color:#7d9284;text-transform:uppercase}
 :is(.st-key-mpnn-surface,.st-key-mpnn-node) [data-testid=stButton] button{border:1px solid #b8d0c0;border-radius:6px;color:#29644c;min-height:2.5rem;background:#fff}
 :is(.st-key-mpnn-surface,.st-key-mpnn-node) [data-testid=stButton] button[kind=primary]{background:#27674e;color:white;border-color:#27674e}
 :is(.st-key-mpnn-surface,.st-key-mpnn-node) [data-testid=stButton] button:disabled{background:#a6c1b1;color:white;border-color:#a6c1b1}
 :is(.st-key-mpnn-surface,.st-key-mpnn-node) [data-variant=pills]{border:1px solid #c5d9cc!important;background:white!important;color:#39694f!important;border-radius:5px!important;padding:.35rem .8rem!important}
 :is(.st-key-mpnn-surface,.st-key-mpnn-node) [data-variant=pills]:is([aria-checked=true],[aria-pressed=true]){background:#e4f0e8!important;color:#1f6545!important;border-color:#4c926c!important;box-shadow:inset 0 0 0 1px #4c926c}
 :is(.st-key-mpnn-surface,.st-key-mpnn-node) [data-baseweb=input],:is(.st-key-mpnn-surface,.st-key-mpnn-node) [data-baseweb=textarea],:is(.st-key-mpnn-surface,.st-key-mpnn-node) [data-baseweb=select]>div{background:#f2f6f3!important;border-color:#d8e3dc!important;border-radius:6px!important}
 :is(.st-key-mpnn-surface,.st-key-mpnn-node) :is(input,textarea){background:#f2f6f3!important;color:#243e31!important}
 :is(.st-key-mpnn-surface,.st-key-mpnn-node) [data-baseweb=tag]{background:#dfede4!important;color:#27674e!important}
 :is(.st-key-mpnn-surface,.st-key-mpnn-node) [data-tag]{background:#dfede4!important;color:#27674e!important}
 :is(.st-key-mpnn-surface,.st-key-mpnn-node) [data-tag] :is(span,button){color:#27674e!important}
 :is(.st-key-mpnn-surface,.st-key-mpnn-node) [data-testid=stMultiSelect] .react-aria-ComboBox>[role=group]{background:#f2f6f3;border-color:#d8e3dc}
 :is(.st-key-mpnn-surface,.st-key-mpnn-node) [data-testid=stCheckbox] label>div:first-of-type{border-color:#a4c0af}
 :is(.st-key-mpnn-surface,.st-key-mpnn-node) [data-testid=stCheckbox] label:has(input:checked)>div:first-of-type{background:#27674e!important;border-color:#27674e!important}
 :is(.st-key-mpnn-surface,.st-key-mpnn-node) [data-testid=stExpander]{border-color:#d8e4dc;border-radius:6px;background:#fcfefc}
 :is(.st-key-mpnn-surface,.st-key-mpnn-node) [data-testid=stExpander] summary{background:#f5f9f6!important}
 :is(.st-key-mpnn-surface,.st-key-mpnn-node) [data-testid=stFileUploaderDropzone]{border:1px dashed #b8cfbf;border-radius:7px;background:#f6faf7}
 :is(.st-key-mpnn-surface,.st-key-mpnn-node) [data-testid=stCaptionContainer]{color:#7d8e82}
 .st-key-mpnn-surface [data-testid=stLinkButton] a{border:0;color:#547963;font-size:.77rem;min-height:1.8rem;padding:.2rem}
 .mpnn-divider{border-top:1px solid #dfe8e2;margin:.6rem 0}
 </style>"""

EXAMPLES = {
    "9F6O · Nanobody CDRs": dict(structure_path=str(inputs.ASSETS / "9F6O.pdb"), design_chains="B", designed_residues="B26-32 B52-56 B95-102", model_type="protein_mpnn"),
    "1BC8 · DNA context": dict(structure_path=str(inputs.ASSETS / "1BC8.pdb"), design_chains="C", designed_residues="", model_type="ligand_mpnn"),
}


def _load_example(label, prefix, owner):
    p = inputs.defaults() | EXAMPLES[label]
    source = Path(p["structure_path"])
    destination = user_input_root(owner) / "mpnn-examples" / source.name
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(source.read_bytes())
    p["structure_path"] = str(destination)
    for key in list(st.session_state):
        if key.startswith(prefix + "-field-") or key.startswith(prefix + "-source-"):
            del st.session_state[key]
    st.session_state[prefix + "-preset"] = p


def _structure_fields(p, prefix, owner, *, single_chain=False, allow_non_protein=True, abmpnn_only=False):
    fmt = st.pills("Structure format", ["PDB", "CIF"], default="PDB", selection_mode="single", key=prefix + "-source-format", label_visibility="collapsed") or "PDB"
    st.caption("Protein structure file to be designed")
    source = render_target_source(key_prefix=prefix + "-source", label="Structure", owner_id=owner,
        example_code="9F6O", allow_cif=True, allow_non_protein=allow_non_protein, upload_in_expander=False,
        upload_label="Drag & drop or browse your structure file", upload_file_types=["pdb", "ent"] if fmt == "PDB" else ["cif", "mmcif"],
        upload_help="Keep ligand and nucleic-acid atoms for LigandMPNN. The original structure is saved with the job.")
    if not source.path and p.get("structure_path") and Path(p["structure_path"]).is_file():
        source = TargetSource(path=Path(p["structure_path"]), name=Path(p["structure_path"]).stem, owner_id=owner)
    if not source.path:
        st.multiselect("Designed Chains", [], disabled=True, key=prefix + "-source-unloaded-chains")
        st.caption("Load a structure and choose chains to enable residue selection.")
        if not abmpnn_only:
            p["designed_residues"] = st.text_input("Designed Residues", value=p["designed_residues"], disabled=True, key=prefix + "-field-unloaded-residues")
        return p, "Upload a PDB/mmCIF structure or fetch a PDB accession."
    p["structure_path"] = str(source.path)
    info = inputs.inventory(source.path)
    chains = [c for c, rows in info["chains"].items() if rows["kind"] == "protein"]
    if not chains:
        return p, "The structure needs at least one protein chain."
    digest = hashlib.sha256(str(source.path).encode()).hexdigest()[:10]
    key = prefix + "-field-" + digest
    default = [c for c in inputs.chain_list(p["design_chains"]) if c in chains] or chains[:1]
    if single_chain:
        chosen = [st.selectbox("Chain", chains, index=chains.index(default[0]), key=key + "-chain")]
    else:
        chosen = st.multiselect("Designed Chains", chains, default=default, key=key + "-chains",
                               help="Only these chains are redesigned. All other chain sequences remain fixed.")
    p["design_chains"] = ",".join(chosen)
    if source.path.suffix.lower() in {".pdb", ".ent"} and st.checkbox("Trim PDB", key=key + "-trim", disabled=not chosen):
        from mn_molecular_design_studio.target_editor import (
            TargetDesignInput, _render_target_trimming, _protein_residue_ids_by_chain,
            _trim_target_chain_ranges, _derived_target_storage_path,
        )
        editor = TargetDesignInput(label="designed", key_prefix=key, source_path=source.path, chains=chosen)
        original = source.path.read_text()
        ranges = _render_target_trimming(editor, _protein_residue_ids_by_chain(original))
        if not editor.valid_trim_ranges:
            return p, "Retain at least one residue per designed chain."
        trimmed = _trim_target_chain_ranges(original, ranges)
        path = _derived_target_storage_path(trimmed, source.name + "-trimmed.pdb", owner_id=owner)
        source = TargetSource(path=path, name=source.name + " · trimmed", owner_id=owner)
        p["structure_path"] = str(path)
        info = inputs.inventory(path)
    selection_key = key + "-residues"
    pending = selection_key + "-pending"
    if pending in st.session_state:
        st.session_state[selection_key] = st.session_state.pop(pending)
    if abmpnn_only:
        selected = []
        st.caption(f"{source.name} · {len(chosen)} antibody chain(s) · CDR positions assigned automatically by SAbR")
    else:
        p["designed_residues"] = st.text_input("Designed Residues", value=p["designed_residues"], key=selection_key,
            placeholder="e.g. B26-32 B52-56 B95-102", help="Use author residue numbers, ranges, or individual IDs. Leave blank to design every residue in the selected chains.")
        try:
            selected = inputs.residues(p["designed_residues"], info, chains=chosen)
        except ValueError as exc:
            selected = []
            st.error(str(exc))
        st.caption(f"{source.name} · {len(chosen)} designed chain(s) · " + (f"{len(selected)} selected residues" if selected else "all residues in those chains"))
    with st.expander("Structure & residue selection", expanded=True):
        a, b, _ = st.columns([1, 1, 3])
        if a.button("Select all", key=key + "-all", use_container_width=True):
            st.session_state[pending] = " ".join(r["label"] for c in chosen for r in info["chains"][c]["residues"])
            st.rerun()
        if b.button("Clear", key=key + "-clear", use_container_width=True):
            st.session_state[pending] = ""
            st.rerun()
        value = molstar_custom_component([StructureVisualization(pdb=source.path.read_text(), color="chain-id",
            representation_type="cartoon", highlighted_selections=selected)], key=key + "-viewer",
            height=400, show_controls=True, show_structure_tools=False, selection_mode=True,
            force_reload=int(hashlib.sha256(p["designed_residues"].encode()).hexdigest()[:7], 16))
        signature = json.dumps(value, sort_keys=True, default=str)
        seen = st.session_state.get(key + "-seen")
        st.session_state[key + "-seen"] = signature
        if seen is not None and signature != seen:
            from mn_molecular_design_studio.target_editor import _selected_viewer_residues, _residue_map
            picked = _selected_viewer_residues(value, _residue_map(source.path.read_text(), source.path))
            labels = [f"{c}{n}" for c, n in sorted(picked) if c in chosen and f"{c}{n}" in info["residues"]]
            st.session_state[pending] = " ".join(labels)
            st.rerun()
    return p, None if chosen else "Choose at least one designed chain."


def render_inputs(*, owner_id, key_prefix, initial=None, connected=False, abmpnn_only=False):
    style = STYLE.replace("mpnn-surface", "abmpnn-surface") if abmpnn_only else STYLE
    st.html(style)
    p = inputs.defaults() | copy.deepcopy(initial or {})
    if abmpnn_only:
        p.update(model_type="ab_mpnn", auto_detect_cdrs=True)
    if not connected:
        p.update(st.session_state.get(key_prefix + "-preset", {}))
    error = None
    prefix = key_prefix
    if not connected:
        st.markdown("**PDB File**")
        p, error = _structure_fields(p, prefix, owner_id, abmpnn_only=abmpnn_only)
    else:
        st.caption("Connected structures provide the backbone and chain roles for this step.")
        p["design_chains"] = st.text_input("Designed Chains", value=p["design_chains"], key=prefix + "-field-chains",
            placeholder="Automatic · designed/binder chains from the previous step", help="Leave blank to preserve the upstream target/design chain assignments.")
        if not abmpnn_only:
            p["designed_residues"] = st.text_input("Designed Residues", value=p["designed_residues"], key=prefix + "-field-residues", placeholder="All designed-chain residues, or e.g. B26-32 B52-56")
        else:
            st.caption("SAbR will determine the selected CDR residue positions for each candidate.")
    p["num_seq_per_target"] = int(st.number_input("Number of Sequences", min_value=1, max_value=10000, value=int(p["num_seq_per_target"]), key=prefix + "-field-count",
        help="Number of sequences generated for each input backbone."))
    if abmpnn_only:
        model = "ab_mpnn"
        p["model_type"] = model
        p["cdr_regions"] = st.multiselect(
            "CDR regions to redesign", ["CDR1", "CDR2", "CDR3"],
            default=p.get("cdr_regions") or ["CDR1", "CDR2", "CDR3"],
            key=prefix + "-field-cdr-regions",
            help="SAbR identifies these loops using IMGT numbering. Antibody framework and unselected CDRs remain fixed.",
        )
        p["sabr_chain_type"] = "auto"
        p["sabr_mode"] = "sabr"
        p["sabr_noise_level"] = 0.0
        st.caption("AbMPNN redesigns the selected loops. SAbR identifies CDR positions from the antibody backbone; the target and framework stay fixed.")
        if p.get("structure_path") and Path(p["structure_path"]).is_file() and p.get("design_chains"):
            try:
                from mn_molecular_design_studio.sabr_runtime import cdr_fixed_residues
                _, selected_cdrs = cdr_fixed_residues(
                    Path(p["structure_path"]), inputs.chain_list(p["design_chains"]),
                    p["cdr_regions"], chain_type="auto",
                )
                parts = [f"{name}: {len(selected_cdrs[name])} residues" for name in p["cdr_regions"]]
                st.success("SAbR found " + " · ".join(parts))
            except Exception as exc:
                st.warning(f"SAbR will validate the selected antibody chain when the job starts: {exc}")
    else:
        model = st.pills("Model Type", list(inputs.MODELS), format_func=lambda m: inputs.MODELS[m]["label"],
            default=p["model_type"], selection_mode="single", key=prefix + "-field-model") or p["model_type"]
        p["model_type"] = model
        st.caption(inputs.MODELS[model]["description"])
        if model == "ab_mpnn":
            st.caption("Select antibody chains and the CDR or framework residues to redesign. For paired antibodies choose heavy before light.")
    if not connected:
        p["verify_sequences"] = st.selectbox("Verify Sequences", [False, True], format_func=lambda v: "ESMFold · each designed chain" if v else "None",
            index=int(p.get("verify_sequences", False)), key=prefix + "-field-verify", help="Optionally queue structure verification after sequence design. Other folding engines can be connected on the pipeline canvas.")
    with st.expander("View Optional Settings", expanded=st.session_state.get(prefix + "-show-settings", False)):
        p["sampling_temp"] = float(st.number_input("Temperature", min_value=.000001, max_value=1., value=float(p["sampling_temp"]), step=.05, format="%.4f", key=prefix + "-field-temperature",
            help="Suggested values: 0.1, 0.15, 0.2, 0.25, 0.3. Higher values give more sequence diversity."))
        levels = inputs.MODELS[model]["noise"]
        if levels:
            current = p["noise_level"] if p["noise_level"] in levels else inputs.MODELS[model]["default"]
            p["noise_level"] = st.selectbox("Noise Level", levels, index=levels.index(current), format_func=lambda n: f"{n:g} Å", key=prefix + "-field-noise-" + model,
                help="Select the checkpoint trained with this backbone noise. This does not add noise to your input coordinates.")
        else:
            p["noise_level"] = None
        p["omit_aas"] = st.text_input("Amino Acids to Omit", value=p["omit_aas"], key=prefix + "-field-omit", help="Default C avoids new cysteines. Use 'empty' to omit none.")
        p["seed"] = st.number_input("Seed", min_value=1, max_value=4294967295, value=p["seed"] or None, step=1, key=prefix + "-field-seed", placeholder="Random seed", help="Leave blank to draw a random seed.")
        p["bias_aa"] = st.text_area("Bias AA", value=p["bias_aa"], placeholder="W:3.0,P:3.0,C:3.0,A:-3.0", key=prefix + "-field-bias")
        for field, label, hint in [("bias_aa_per_residue", "Bias AA per Residue", '{"A12": {"G": -0.3, "P": 2.0}}'), ("omit_aa_per_residue", "Omit AA per Residue", '{"A12": "CP", "A13": "CW"}')]:
            raw = p[field]
            p[field] = st.text_area(label, value=json.dumps(raw) if isinstance(raw, dict) and raw else raw or "", placeholder=hint, key=prefix + "-field-" + field)
        p["homo_oligomer"] = int(st.number_input("Homo Oligomer", min_value=0, max_value=62, value=int(p["homo_oligomer"]), key=prefix + "-field-homo",
            help="Number of existing designed-chain copies to tie (2 = dimer, 3 = trimer). Copies must already be in the input and have matching residue counts. 0/1 disables tying."))
        p["fixed_residues"] = st.text_input("Fixed Residues", value=p["fixed_residues"], key=prefix + "-field-fixed", help="These residues and all upstream fixed motifs retain their sequence.")
        p["fix_interface_residues"] = st.checkbox(
            "Fix interface residues from valid generator sequences",
            value=p["fix_interface_residues"],
            key=prefix + "-field-fix-interface",
            help="For upstream generator candidates only, keep binder residues within 4 Å of the target fixed when the candidate has a valid amino-acid sequence matching its structure. Candidates without a valid sequence are redesigned normally.",
        )
        if model == "global_label_membrane_mpnn":
            p["global_transmembrane_label"] = st.checkbox("Transmembrane protein", value=p["global_transmembrane_label"], key=prefix + "-field-global-membrane")
        if model == "per_residue_label_membrane_mpnn":
            p["transmembrane_buried"] = st.text_input("Membrane-buried Residues", value=p["transmembrane_buried"], key=prefix + "-field-buried")
            p["transmembrane_interface"] = st.text_input("Membrane-interface Residues", value=p["transmembrane_interface"], key=prefix + "-field-interface")
        if model == "ligand_mpnn":
            p["use_side_chain_context"] = st.checkbox("Use fixed side chains as atomic context", value=p["use_side_chain_context"], key=prefix + "-field-side-context")
        else:
            p["use_side_chain_context"] = False
        p["pack_side_chains"] = st.checkbox("Pack designed side chains", value=p["pack_side_chains"], key=prefix + "-field-packing")
        if p["pack_side_chains"]:
            p["number_of_packs_per_design"] = int(st.number_input("Packing samples per sequence", 1, 100, max(1, int(p["number_of_packs_per_design"])), key=prefix + "-field-pack-count"))
        p["backend"] = st.selectbox("Execution backend", ["auto", "shared_ligandmpnn", "foundry"], index=["auto", "shared_ligandmpnn", "foundry"].index(p["backend"]),
            format_func=lambda v: {"auto": "Automatic", "shared_ligandmpnn": "Shared MPNN", "foundry": "Foundry MPNN"}[v], key=prefix + "-field-backend",
            help="Automatic uses the native Foundry path for compatible CIF backbones and the shared MPNN image for other model families.")
    try:
        p = inputs.prepare(p)
    except (ValueError, TypeError, OSError) as exc:
        error = str(exc)
    return p, error


def _existing_results(project, search_jobs, prefix, *, antibody_backbones_only=False):
    if not project:
        st.info("Select a project to choose its completed design results.")
        return {}, "Select a project."
    rows = []
    page = {"items": []}
    cursor = None
    while True:
        result = search_jobs(project_id=project, statuses=["completed", "succeeded"], limit=1000, cursor=cursor)
        page["items"].extend(result["items"])
        cursor = result.get("next_cursor")
        if not cursor:
            break
    from mn_molecular_design_studio.engine.workflow_engine.workflows.mpnn import STAGES
    for row in page["items"]:
        run = Path(row["run_dir"])
        if read_json(run / "result.json").get("success") is not True:
            continue
        allowed_stages = {"generation.backbone", "generation.backbone_sequence"} if antibody_backbones_only else set(STAGES)
        candidates = [c for c in read_candidates(run / "artifacts/normalized_candidates/candidates.jsonl") if c.get("stage") in allowed_stages and c.get("complex_pdb")]
        if candidates:
            rows.append((row, candidates))
    if not rows:
        st.info("No completed jobs with designable structures are available in this project.")
        return {}, "Choose existing design results."
    index = st.selectbox("Source job", range(len(rows)), format_func=lambda i: f"{rows[i][0].get('campaign_name') or rows[i][0]['tool']} · {rows[i][0].get('job_code')}", key=prefix + "-source-job-" + project)
    row, candidates = rows[index]
    ids = [c["candidate_id"] for c in candidates]
    selected = st.multiselect("Designs", ids, default=ids, key=prefix + "-source-designs-" + row["run_id"])
    return {"source_run_dir": row["run_dir"], "candidates_jsonl": str(Path(row["run_dir"]) / "artifacts/normalized_candidates/candidates.jsonl"), "selected_candidate_ids": selected}, None if selected else "Select at least one input design."


def render_setup(schema, *, owner_id, project_selector, submit, search_jobs, abmpnn_only=False):
    prefix = "abmpnn" if abmpnn_only else "mpnn"
    style = STYLE.replace("mpnn-surface", "abmpnn-surface") if abmpnn_only else STYLE
    st.html(style)
    if st.button("← Back to All tools", key=prefix + "-back"):
        st.session_state["page"] = "All tools"
        st.rerun()
    with st.container(key=prefix + "-surface"):
        title = "AbMPNN" if abmpnn_only else "ProteinMPNN"
        subtitle = "Antibody structure to CDR sequence redesign · SAbR IMGT loop detection" if abmpnn_only else "Protein structure to sequence design · ProteinMPNN family, LigandMPNN, HyperMPNN & AbMPNN"
        st.html(f'<div class="mpnn-header"><div class="mpnn-eyebrow">Inverse folding</div><h1>{title}</h1><p>{subtitle}</p></div>')
        a, b, c, d = st.columns(4)
        if a.button("Show outputs", key=prefix + "-outputs", use_container_width=True):
            st.session_state[prefix + "-show-outputs"] = not st.session_state.get(prefix + "-show-outputs", False)
        if b.button("Show settings", key=prefix + "-settings", use_container_width=True):
            st.session_state[prefix + "-show-settings"] = not st.session_state.get(prefix + "-show-settings", False)
        c.link_button("Source code", "https://github.com/dauparas/LigandMPNN", use_container_width=True)
        d.link_button("Input reference", "https://app.tamarind.bio/tools/abmpnn" if abmpnn_only else "https://app.tamarind.bio/tools/proteinmpnn", use_container_width=True)
        if st.session_state.get(prefix + "-show-outputs"):
            st.info("Each input backbone produces designed sequences, matching structure files, model scores and parameters. Outputs can be connected to folding and analysis on the pipeline canvas.")
        with st.expander("About this workflow"):
            if abmpnn_only:
                st.write("SAbR assigns antibody numbering and identifies CDR loops from the supplied backbone. AbMPNN redesigns only the selected loops; the target, framework, and unselected loops are kept fixed.")
            else:
                st.write("Design all or selected residues on a fixed protein backbone. Keep the target, other chains and fixed motifs unchanged. Use LigandMPNN to retain atomic context around non-protein molecules.")
        a, b = st.columns(2)
        name = a.text_input("Job name", "abmpnn-cdr-redesign" if abmpnn_only else "proteinmpnn-design", key=prefix + "-name").strip()
        with b:
            project = project_selector(prefix)
        mode = st.segmented_control("Input source", ["Structure", "Existing results"], default="Structure", key=prefix + "-mode") or "Structure"
        source, source_error = ({}, None) if mode == "Structure" else _existing_results(project, search_jobs, prefix, antibody_backbones_only=abmpnn_only)
        if mode == "Structure":
            a, b = st.columns(2)
            examples = {"9F6O · Nanobody CDRs": EXAMPLES["9F6O · Nanobody CDRs"]} if abmpnn_only else EXAMPLES
            for column, label in zip((a, b), examples):
                column.button(label, key=prefix + "-example-" + label, on_click=_load_example, args=(label, prefix, owner_id), use_container_width=True)
        parameters, error = render_inputs(owner_id=owner_id, key_prefix=prefix, connected=mode != "Structure", abmpnn_only=abmpnn_only)
        if mode != "Structure":
            parameters["verify_sequences"] = st.selectbox("Verify Sequences", [False, True],
                format_func=lambda v: "ESMFold · each designed chain" if v else "None", key=prefix + "-field-verify-existing")
        parameters["job_name"] = name or ("abmpnn-cdr-redesign" if abmpnn_only else "proteinmpnn-design")
        parameters.update(source)
        if mode != "Structure":
            parameters.pop("structure_path", None)
        with st.expander("Scheduler & reference files"):
            gpu = st.selectbox("GPU assignment", ["Automatic", "GPU 0", "GPU 1"], key=prefix + "-gpu")
            parameters["gpu_device"] = {"Automatic": "all", "GPU 0": "0", "GPU 1": "1"}[gpu]
            priority = st.checkbox("Prioritize in scheduler queue", key=prefix + "-priority")
            readiness = inputs.readiness(parameters["model_type"], parameters["noise_level"])
            if readiness["ready"]:
                st.success("Selected checkpoint is available in reference files.")
            else:
                st.warning("Selected checkpoint is missing: " + readiness["checkpoint"])
        error = source_error or error
        if error:
            st.caption(error)
        with st.columns([1, 1, 1])[1]:
            clicked = st.button("Submit", type="primary", key=prefix + "-submit", use_container_width=True, disabled=not project or bool(error) or not readiness["ready"])
        if not project:
            st.caption("Select a project to submit.")
        if clicked:
            allowed = {f["name"] for f in schema["parameters"]}
            submit(schema, schema["workflow"], {k: v for k, v in parameters.items() if k in allowed}, project, {"prioritize": priority}, name or ("AbMPNN" if abmpnn_only else "ProteinMPNN"))
