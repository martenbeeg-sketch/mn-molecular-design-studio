"""RFdiffusion3 setup, sharing Studio project, structure and scheduler controls."""
from __future__ import annotations

import hashlib
import io
import json
from pathlib import Path
import shutil
import zipfile

import pandas as pd
import streamlit as st

from mn_molecular_design_studio import rfdiffusion3_input as inputs
from mn_molecular_design_studio.auth import user_input_root
from mn_molecular_design_studio.target_editor import render_target_source, render_target_editor
from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer import StructureVisualization, molstar_custom_component

EXAMPLES = {
    "protein_binder": "5o45_cropped.pdb", "enzyme": "M0255_1mg5.pdb",
    "nucleic_acid_binder": "1bna.pdb", "ligand_binder": "IAI.pdb",
}
ENZYME_ATOMS = {"A108": ["ND2", "CG"], "A139": ["OG", "CB", "CA"],
                "A152": ["OH", "CZ"], "A156": ["NZ", "CE", "CD"]}


def store_json_structures(uploaded, owner_id: str) -> list[Path]:
    """Preserve JSON input basenames; never extract archive paths or symlinks."""
    content = []
    for upload in uploaded:
        if upload.name.lower().endswith(".zip"):
            with zipfile.ZipFile(io.BytesIO(upload.getvalue())) as archive:
                for item in archive.infolist():
                    if item.is_dir():
                        continue
                    if item.file_size > 50_000_000 or item.flag_bits & 1:
                        raise ValueError("ZIP members must be unencrypted and smaller than 50 MB.")
                    if Path(item.filename).suffix.lower() not in {".pdb", ".ent", ".cif", ".mmcif"}:
                        continue
                    content.append((Path(item.filename).name, archive.read(item)))
                    if len(content) > 128 or sum(len(data) for _, data in content) > 100_000_000:
                        raise ValueError("Upload at most 128 structures and 100 MB in total.")
        else:
            content.append((Path(upload.name).name, upload.getvalue()))
    if len({name for name, _ in content}) != len(content):
        raise ValueError("Referenced structures must have unique filenames, including ZIP members.")
    if len(content) > 128 or sum(len(data) for _, data in content) > 100_000_000:
        raise ValueError("Upload at most 128 structures and 100 MB in total.")
    result = []
    for name, data in content:
        if len(data) > 50_000_000:
            raise ValueError("Structure files must be smaller than 50 MB.")
        destination = user_input_root(owner_id) / "rfdiffusion3-json" / hashlib.sha256(data).hexdigest()[:16] / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(data)
        inputs.inventory(destination)
        result.append(destination)
    return result


def render_setup(schema: dict, *, owner_id: str, header, project_selector, submit) -> None:
    st.html("""<style>
    .block-container {max-width:1240px!important;padding-top:1.5rem!important}
    .st-key-rfd3-surface {--rfd-green:#27674e;--rfd-ink:#163d31;--rfd-mint:#edf6f0}
    .st-key-rfd3-surface .rfd3-heading {text-align:center;padding:.5rem 0 1.7rem}
    .st-key-rfd3-surface .rfd3-eyebrow {font-size:.69rem;letter-spacing:.13em;font-weight:650;color:#658677;margin-bottom:.7rem}
    .st-key-rfd3-surface .rfd3-heading h1 {font-size:2.1rem;color:#173d31;margin:0;font-weight:700;letter-spacing:-.04em}
    .st-key-rfd3-surface .rfd3-heading p {font-size:.92rem;color:#6c8276;margin:.65rem auto .6rem;max-width:570px;line-height:1.6}
    .st-key-rfd3-surface .rfd3-meta {font-size:.73rem;color:#789084}
    .st-key-rfd3-surface [data-testid="stVerticalBlockBorderWrapper"] {border-color:#d7e6dc!important;border-radius:12px!important}
    .st-key-rfd3-input-card {padding:1.3rem!important;background:#fff;border:1px solid #d7e6dc;border-radius:12px;margin:.5rem 0 .8rem}
    .st-key-rfd3-sampling-card {padding:1.1rem 1.3rem!important;background:#f3f8f5;border:1px solid #d7e6dc;border-radius:12px;margin:.3rem 0 .6rem}
    .st-key-rfd3-surface .rfd3-section {font-size:.72rem;font-weight:650;letter-spacing:.1em;color:#5c7b68;margin:0 0 .8rem}
    .st-key-rfd3-surface .rfd3-task-description {color:#537260;font-size:.88rem;line-height:1.5;margin:.2rem 0 .4rem}
    .st-key-rfd3-surface [data-variant="pills"] {border-radius:7px!important;border:1px solid #c2d9cc!important;background:#f7faf8;color:#38614e!important;min-height:2.6rem;padding:.45rem .9rem!important;font-size:.84rem!important;box-shadow:none!important}
    .st-key-rfd3-surface [data-variant="pills"][aria-checked="true"] {background:#27674e!important;color:white!important;border-color:#27674e!important}
    .st-key-rfd3-surface [data-testid="stTextInput"] input, .st-key-rfd3-surface [data-testid="stNumberInput"] input, .st-key-rfd3-surface textarea {background:#f5f9f6!important;color:#193e30!important}
    .st-key-rfd3-surface [data-baseweb="input"], .st-key-rfd3-surface [data-baseweb="textarea"] {border-color:#dce7df!important;border-radius:7px!important;background:#f5f9f6!important}
    .st-key-rfd3-surface [data-baseweb="tag"] {background:#e0eee5!important;color:#27674e!important}
    .st-key-rfd3-surface [data-testid="stButton"] button {border-radius:7px!important;color:#27674e;border-color:#bad2c3}
    .st-key-rfd3-surface [data-testid="stButton"] button:hover {background:#edf6f0!important;border-color:#659d7d!important;color:#1d5039!important}
    .st-key-rfd3-surface button[kind="primary"] {background:#27674e!important;border-color:#27674e!important;color:#fff!important;min-height:2.8rem}
    .st-key-rfd3-surface button[kind="primary"]:disabled {background:#99b8a5!important;border-color:#99b8a5!important;color:#fff!important}
    .st-key-rfd3-surface [data-testid="stExpander"] {border-color:#d7e6dc!important;border-radius:8px!important;background:#fcfefc}
    .st-key-rfd3-surface [data-testid="stFileUploaderDropzone"] {background:#f3f8f5!important;border:1px dashed #b9d2c2!important;border-radius:8px}
    .st-key-rfd3-surface [role="radio"][aria-checked="true"] {background:#27674e!important;border-color:#27674e!important;color:#fff!important}
    .st-key-rfd3-surface [data-testid="stCaptionContainer"] {color:#7c8f82}
    </style>""")
    with st.container(key="rfd3-surface"):
        _render_setup_content(schema, owner_id=owner_id, header=header, project_selector=project_selector, submit=submit)


def _render_setup_content(schema: dict, *, owner_id: str, header, project_selector, submit) -> None:
    if st.button("← Back to All tools", key="rfd3-back"):
        st.session_state["page"] = "All tools"
        st.session_state.pop("selected_workflow", None)
        st.rerun()
    st.html('<div class="rfd3-heading"><div class="rfd3-eyebrow">ALL-ATOM PROTEIN DESIGN</div><h1>RFdiffusion3</h1>'
            '<p>Design protein backbones around proteins, catalytic motifs, DNA, RNA and small molecules.</p>'
            f'<span class="rfd3-meta">{inputs.CODE_VERSION}</span></div>')
    with st.expander("About this workflow"):
        st.caption(inputs.CODE_VERSION)
        st.markdown("[RFdiffusion3 documentation](https://github.com/RosettaCommons/foundry/tree/production/models/rfd3) · "
                    "[Tamarind input reference](https://app.tamarind.bio/tools/rfdiffusion3)")
        st.write("This step generates all-atom backbones. Use Studio’s sequence-design tools next, then refold and filter the designed sequences. "
                 "Generated backbones are not classified as filter-passing designs.")
        st.caption("The installed engine uses standard inference. Claude inference-optimization modes shown on Tamarind are not available in this image.")
    a, b = st.columns(2)
    with a:
        job_name = st.text_input("Job name", "rfdiffusion3-design", key="rfd3-job-name").strip()
    with b:
        project = project_selector("rfd3")
    st.html('<div class="rfd3-section">DESIGN TASK</div>')
    task = st.pills("Task", list(inputs.TASKS), format_func=inputs.TASKS.get, selection_mode="single", default="protein_binder", key="rfd3-task-pills", label_visibility="collapsed") or "protein_binder"
    descriptions = {"protein_binder": "Design a new protein that binds a selected surface on your target.",
                    "enzyme": "Build a protein scaffold around catalytic residues and a bound molecule.",
                    "nucleic_acid_binder": "Design a protein binder around a DNA or RNA structure.",
                    "ligand_binder": "Build a protein binding pocket around a small molecule.",
                    "json": "Use the native specification for custom motifs, symmetry and multiple design inputs."}
    st.html(f'<p class="rfd3-task-description">{descriptions[task]}</p>')
    prefix = f"rfd3-{task}"
    p = {"task": task, "job_name": job_name or "rfdiffusion3-design"}
    error = None
    with st.container(key="rfd3-input-card"):
        if task == "json":
            st.caption("Provide named RFdiffusion3 dialect-2 specifications. Input filenames must match the uploaded structures. Each specification receives the number of designs below.")
            config_file = st.file_uploader("JSON configuration", type=["json"], key=prefix+"-json-upload")
            text = st.text_area("Configuration JSON", value='{"unconditional": {"length": "100-150", "is_non_loopy": true}}', height=220, key=prefix+"-json")
            uploaded = st.file_uploader("Referenced structures (optional)", type=["pdb", "ent", "cif", "mmcif", "zip"], accept_multiple_files=True, key=prefix+"-files")
            try:
                p.update(json_config=json.loads(config_file.getvalue() if config_file else text), input_files=[str(f) for f in store_json_structures(uploaded, owner_id)])
            except (ValueError, OSError, zipfile.BadZipFile) as exc:
                error = str(exc)
        else:
            if st.button("Load example inputs", icon=":material/science:", key=prefix+"-example"):
                source = inputs.ASSETS / EXAMPLES[task]
                destination = user_input_root(owner_id) / "structures" / source.name
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source, destination)
                st.session_state[prefix+"-fetched-target"] = str(destination)
                st.session_state[prefix+"-fetched-target-name"] = source.stem
                st.session_state.pop(prefix+"-target-upload", None)
                st.rerun()
            source = render_target_source(key_prefix=prefix, label="Target / motif structure", example_code="1BNA" if task == "nucleic_acid_binder" else "5O45", owner_id=owner_id,
                                          allow_cif=True, allow_non_protein=True, upload_help="PDB or mmCIF. Protein, DNA/RNA and ligand-only structures are supported.", upload_label="Upload structure")
            info = None
            if source.path:
                p["target_structure"] = str(source.path)
                try:
                    info = inputs.inventory(source.path)
                    if task == "protein_binder":
                        target = render_target_editor(source, key_prefix=prefix+"-editor", show_structure_tools=False)
                        p.update(target_structure=str(target.design_path or source.path), target_chains=target.chains, hotspots=target.hotspots)
                    else:
                        st.success(f"{source.name} loaded")
                        if task == "nucleic_acid_binder":
                            chains = [c for c, row in info["chains"].items() if row["kind"] in {"dna", "rna"}]
                            p["target_chains"] = st.multiselect("DNA / RNA target chains", chains, default=chains, key=prefix+"-chains")
                            p["hotspots"] = st.text_input("Hotspot residues (optional)", placeholder="A3,A7", key=prefix+"-hotspots")
                        with st.expander("Structure preview", expanded=False):
                            molstar_custom_component([StructureVisualization(pdb=source.path.read_text(), color="molecule-type", representation_type="cartoon+ball-and-stick")], key=prefix+"-viewer", height=450, show_structure_tools=False)
                    ligand_options = list(info["ligands"])
                    if ligand_options:
                        p["ligands"] = st.multiselect("Ligands to retain", ligand_options, default=ligand_options if task in {"enzyme", "ligand_binder"} else [], key=prefix+"-ligands")
                    elif task == "ligand_binder":
                        st.warning("Load a structure containing a small-molecule ligand.")
                    if task == "enzyme":
                        mode = st.radio("Scaffolding mode", ["atom", "residue"], format_func=lambda x: "Atom-level motif" if x == "atom" else "Residue-level motif", horizontal=True, key=prefix+"-mode")
                        p["scaffolding_mode"] = mode
                        if mode == "atom":
                            options = [label for label, row in info["residues"].items() if row["kind"] == "protein"]
                            selected = st.multiselect("Catalytic motif residues", options, default=[r for r in ENZYME_ATOMS if r in options], key=prefix+"-motifs")
                            fixed = {}
                            for residue in selected:
                                options = info["residues"][residue]["atoms"]
                                atoms = st.multiselect(f"Fixed atoms · {residue} ({info['residues'][residue]['name']})", options,
                                                       default=[a for a in ENZYME_ATOMS.get(residue, options) if a in options], key=prefix+"-atoms-"+residue+"-"+source.path.stem)
                                fixed[residue] = ",".join(atoms)
                            p["fixed_atoms"] = fixed
                            p["unindex"] = st.text_input("Unindexed motif residues", "all", help="all places the selected motif residues without fixed sequence positions. Alternatively enter selected residue IDs, e.g. A108,A139.", key=prefix+"-unindex")
                        else:
                            st.caption("Alternate designed segment lengths with retained motif residues. Example: 20,A108,20,A139,20,A152,20,A156,20. Residue numbers refer to the loaded file.")
                            p["motif_contig"] = st.text_input("Scaffold contig", "20,A108,20,A139,20,A152,20,A156,20", key=prefix+"-contig")
                except (ValueError, OSError, KeyError) as exc:
                    error = str(exc)
            if not (task == "enzyme" and p.get("scaffolding_mode") == "residue"):
                a,b = st.columns(2)
                with a:
                    p["length_min"] = int(st.number_input("Minimum scaffold length" if task == "enzyme" else "Minimum binder length", 1, 3000, 100, key=prefix+"-min"))
                with b:
                    p["length_max"] = int(st.number_input("Maximum scaffold length" if task == "enzyme" else "Maximum binder length", 1, 3000, 150, key=prefix+"-max"))
            if p.get("ligands"):
                st.markdown("**Ligand atom conditions**")
                st.caption("Use comma-separated atom names or ALL. A Fixed condition with empty atoms allows that ligand to diffuse. Other conditions need atoms or ALL.")
                edited = st.data_editor(pd.DataFrame(columns=["condition", "ligand", "atoms"]), num_rows="dynamic", hide_index=True,
                                        column_config={"condition": st.column_config.SelectboxColumn("Condition", options=list(inputs.CONDITIONS), required=True),
                                                       "ligand": st.column_config.SelectboxColumn("Ligand", options=p["ligands"], required=True),
                                                       "atoms": st.column_config.TextColumn("Atom names")}, key=prefix+"-conditions")
                p["atom_conditions"] = edited.fillna("").to_dict("records")
    with st.container(key="rfd3-sampling-card"):
        st.html('<div class="rfd3-section">DESIGN SAMPLING</div>')
        a,b = st.columns(2)
        with a:
            p["num_designs"] = int(st.number_input("Number of designs", 1, 10000, 8, key="rfd3-count"))
        with b:
            p["diffusion_batch_size"] = int(st.number_input("Diffusion batch size", 1, 32, 1, key="rfd3-batch"))
        st.caption("Batch size controls GPU memory use. Each batch shares its sampled length. Extra designs from the final batch remain in the raw files; the results table contains the requested number.")
    with st.expander("Optional settings"):
        a,b = st.columns(2)
        with a:
            p["timesteps"] = int(st.number_input("Diffusion timesteps", 1, 1000, 200, key="rfd3-steps"))
            p["step_scale"] = st.number_input("Step scale", 1., 5., 1.5, .1, key="rfd3-scale")
            p["gamma_0"] = st.number_input("Gamma 0", 0., 1., .6, .05, key="rfd3-gamma")
            if task != "json":
                p["is_non_loopy"] = st.checkbox("Condition on non-loopy structures", True, key="rfd3-nonloopy")
        with b:
            p["seed"] = int(st.number_input("Random seed (0 = random)", 0, 2147483647, 0, key="rfd3-seed"))
            p["cpu_cores"] = int(st.number_input("CPU cores", 1, 256, 4, key="rfd3-cpus"))
            p["low_memory_mode"] = st.checkbox(
                "Low memory mode (experimental)",
                False,
                key="rfd3-lowmem",
                help=(
                    "Uses Foundry's chunked pairwise path to reduce GPU memory. "
                    "The pinned Foundry 0.1.9 checkpoint omits some chunked-only weights, "
                    "so this mode may reduce prediction quality."
                ),
            )
            p["classifier_free_guidance"] = st.checkbox("Classifier-free guidance", False, key="rfd3-cfg")
            p["cfg_scale"] = st.number_input("Guidance scale", 0., 10., 1.5, .1, disabled=not p["classifier_free_guidance"], key="rfd3-cfgscale")
        if task != "json":
            ori = st.text_input("Design center (optional [x,y,z])", value="[24,20,10]" if task == "nucleic_acid_binder" else "", key=prefix+"-ori")
            if ori.strip():
                try:
                    p["ori_token"] = json.loads(ori)
                except ValueError:
                    error = "Design center must be three coordinates, e.g. [24,20,10]."
    prioritize = st.checkbox("Prioritize in scheduler queue", key="rfd3-priority")
    if error:
        st.error(error)
    if st.button("Submit design", type="primary", disabled=not project or bool(error), use_container_width=True, key="rfd3-submit"):
        try:
            inputs.prepare_request(p)
            submit(schema, "design.rfdiffusion3", p, project, {"prioritize": prioritize}, "RFdiffusion3")
        except (ValueError, TypeError, OSError) as exc:
            st.error(str(exc))
    if not project:
        st.caption("Select a project to submit.")


def render_results(run_id: str, run_dir: Path, candidates: list[dict], parameters: dict, result: dict) -> None:
    task = parameters.get("task") or (result.get("metrics") or {}).get("task")
    st.subheader(inputs.TASKS.get(task, "RFdiffusion3 designs"))
    st.caption("Generated backbones · sequence design and refolding are still required. These structures have not passed designability filters.")
    if not candidates:
        st.info("Design structures will appear here when RFdiffusion3 finishes.")
        return
    st.dataframe(pd.DataFrame([{"Design": c.get("candidate_id"), "Input": (c.get("raw_metadata") or {}).get("specification_name"),
                                 "Designed chains": ",".join(c.get("binder_chains") or []),
                                 "Designed residues": (c.get("metrics") or {}).get("generated_residue_count"),
                                 "Stage": "Backbone generation"} for c in candidates]), hide_index=True, height=min(360, 35*(len(candidates)+1)+3))
    selected = st.selectbox("Design", range(len(candidates)), format_func=lambda i: candidates[i]["candidate_id"], key=f"rfd3-result-{run_id}")
    candidate = candidates[selected]
    path = Path(candidate.get("complex_pdb") or "")
    if path.is_file():
        molstar_custom_component([StructureVisualization(pdb=path.read_text(), color="molecule-type", representation_type="cartoon+ball-and-stick")], key=f"rfd3-result-view-{run_id}-{selected}", height=560, show_controls=True, show_structure_tools=False)
        st.download_button("Download generated structure", path.read_bytes(), file_name=path.name, key=f"rfd3-result-download-{run_id}")
    with st.expander("Design settings and native metadata"):
        st.json(candidate.get("raw_metadata") or {})
        st.json(candidate.get("metrics") or {})
    st.caption("Use this job as the source in ProteinMPNN / LigandMPNN to design sequences, then validate them with a cofolding engine.")
