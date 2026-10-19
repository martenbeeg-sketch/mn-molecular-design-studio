"""RFantibody standalone input and result pages for Molecular Design Studio."""
from __future__ import annotations

import hashlib
import shutil
from pathlib import Path

import pandas as pd
import streamlit as st

from mn_molecular_design_studio import rfantibody_input as inputs
from mn_molecular_design_studio.auth import user_input_root
from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer import StructureVisualization, molstar_custom_component
from mn_molecular_design_studio.target_editor import _target_widget_keys, render_target_editor, render_target_source


def _save_custom_framework(upload, owner_id: str) -> Path | None:
    if upload is None:
        return None
    content = upload.getvalue()
    digest = hashlib.sha256(content).hexdigest()[:16]
    destination = user_input_root(owner_id) / "rfantibody" / f"custom-framework-{digest}.pdb"
    destination.parent.mkdir(parents=True, exist_ok=True)
    if not destination.exists():
        destination.write_bytes(content)
    return destination


def _load_example_inputs(owner_id: str) -> None:
    target = inputs.EXAMPLE_INPUTS / "flu_HA.pdb"
    stored_target = user_input_root(owner_id) / "rfantibody" / "flu_HA.pdb"
    stored_target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(target, stored_target)
    st.session_state["rfab-task"] = "nanobody"
    st.session_state["rfab-framework"] = "h-NbBCII10"
    st.session_state["rfab-regions"] = ["H1", "H2", "H3"]
    st.session_state["rfab-example-target-path"] = str(stored_target)
    st.session_state["rfab-example-target-name"] = "flu_HA"
    hotspot_key = _target_widget_keys("rfab-target-editor", stored_target)["hotspots"]
    st.session_state[hotspot_key] = "B146,B170,B177"
    st.session_state.pop("rfab-target-upload", None)
    for suffix in ("H1", "H2", "H3", "L1", "L2", "L3"):
        st.session_state[f"rfab-loop-length-{suffix}"] = "auto"


def render_setup(schema: dict, *, owner_id: str, project_selector, submit) -> None:
    st.html("""<style>
    .block-container {max-width:1240px!important;padding-top:1.5rem!important}
    .st-key-rfantibody-page {--rf-green:#2e6a50;--rf-ink:#173d31;--rf-mint:#eef6f1}
    .st-key-rfantibody-page .rfab-heading {text-align:center;padding:.45rem 0 1.35rem}
    .st-key-rfantibody-page .rfab-eyebrow {font-size:.69rem;letter-spacing:.13em;font-weight:650;color:#688575;margin-bottom:.6rem}
    .st-key-rfantibody-page .rfab-heading h1 {font-size:2.1rem;color:#173d31;margin:0;font-weight:700;letter-spacing:-.04em}
    .st-key-rfantibody-page .rfab-heading p {font-size:.92rem;color:#6c8276;margin:.55rem auto .45rem;max-width:600px;line-height:1.6}
    .st-key-rfantibody-page .rfab-meta {font-size:.73rem;color:#789084}
    .st-key-rfantibody-page [data-testid="stVerticalBlockBorderWrapper"] {border-color:#d7e6dc!important;border-radius:12px!important}
    .st-key-rfab-input-card {padding:1.25rem!important;background:#fff;border:1px solid #d7e6dc;border-radius:12px;margin:.45rem 0 .75rem}
    .st-key-rfab-design-count-card {padding:1.05rem 1.25rem!important;background:#f3f8f5;border:1px solid #d7e6dc;border-radius:12px;margin:.35rem 0 .65rem}
    .st-key-rfantibody-page .rfab-section {font-size:.72rem;font-weight:650;letter-spacing:.1em;color:#5c7b68;margin:0 0 .7rem}
    .st-key-rfantibody-page [data-variant="pills"] {border-radius:7px!important;border:1px solid #c2d9cc!important;background:#f7faf8;color:#38614e!important;min-height:2.55rem;padding:.4rem .85rem!important;font-size:.84rem!important;box-shadow:none!important}
    .st-key-rfantibody-page [data-variant="pills"][aria-checked="true"] {background:#2e6a50!important;color:white!important;border-color:#2e6a50!important}
    .st-key-rfantibody-page [data-testid="stTextInput"] input,.st-key-rfantibody-page [data-testid="stNumberInput"] input,.st-key-rfantibody-page textarea {background:#f5f9f6!important;color:#193e30!important}
    .st-key-rfantibody-page [data-baseweb="input"],.st-key-rfantibody-page [data-baseweb="textarea"] {border-color:#dce7df!important;border-radius:7px!important;background:#f5f9f6!important}
    .st-key-rfantibody-page [data-baseweb="tag"] {background:#e0eee5!important;color:#27674e!important}
    .st-key-rfantibody-page [data-testid="stButton"] button {border-radius:7px!important;color:#27674e;border-color:#bad2c3}
    .st-key-rfantibody-page [data-testid="stButton"] button:hover {background:#edf6f0!important;border-color:#659d7d!important;color:#1d5039!important}
    .st-key-rfantibody-page button[kind="primary"] {background:#2e6a50!important;border-color:#2e6a50!important;color:#fff!important;min-height:2.8rem}
    .st-key-rfantibody-page [data-testid="stExpander"] {border-color:#d7e6dc!important;border-radius:8px!important;background:#fcfefc}
    .st-key-rfantibody-page [data-testid="stFileUploaderDropzone"] {background:#f3f8f5!important;border:1px dashed #b9d2c2!important;border-radius:8px}
    .st-key-rfantibody-page [data-testid="stCaptionContainer"] {color:#7c8f82}
    </style>""")
    with st.container(key="rfantibody-page"):
        _render_setup_content(schema, owner_id=owner_id, project_selector=project_selector, submit=submit)


def _render_setup_content(schema: dict, *, owner_id: str, project_selector, submit) -> None:
    if st.button("← Back to All tools", key="rfab-back"):
        st.session_state["page"] = "All tools"
        st.session_state.pop("selected_workflow", None)
        st.rerun()
    st.html(f'<div class="rfab-heading"><div class="rfab-eyebrow">ANTIBODY DESIGN</div><h1>RFantibody</h1>'
            '<p>Design antibody and nanobody binders with RFdiffusion, ProteinMPNN and RF2.</p>'
            f'<span class="rfab-meta">{inputs.CODE_VERSION} · GPU required</span></div>')
    with st.expander("About this tool", expanded=False):
        st.write("RFantibody runs backbone generation, CDR sequence design and complex prediction as one standalone design job.")
        st.markdown("[RFantibody source](https://github.com/RosettaCommons/RFantibody) · "
                    "[Tamarind input reference](https://app.tamarind.bio/tools/rfantibody) · "
                    "[RFantibody paper](https://www.biorxiv.org/content/10.1101/2024.03.14.585103v1)")
        st.caption("Custom frameworks must use Chothia numbering. The optional AbMPNN checkpoint is read from Studio’s reference weights.")

    left, right = st.columns([1.35, 1], vertical_alignment="bottom")
    with left:
        job_name = st.text_input("Job name", "rfantibody-design", key="rfab-job-name").strip() or "rfantibody-design"
    with right:
        project = project_selector("rfantibody")

    st.html('<div class="rfab-section">DESIGN MODE</div>')
    task = st.pills("Task", list(inputs.TASKS), format_func=inputs.TASKS.get, selection_mode="single", default="antibody", key="rfab-task", label_visibility="collapsed") or "antibody"
    framework_choices = ["hu-4D5-8_Fv", "custom"] if task == "antibody" else ["h-NbBCII10", "custom"]
    framework_default = framework_choices[0]
    selected_framework = st.session_state.get("rfab-framework")
    if selected_framework not in framework_choices:
        st.session_state["rfab-framework"] = framework_default
    framework = st.selectbox("Framework", framework_choices, format_func=inputs.FRAMEWORKS.get, key="rfab-framework",
                             help="Choose an RFantibody framework or upload a custom Chothia-numbered framework.")

    with st.container(key="rfab-input-card"):
        st.html('<div class="rfab-section">TARGET ANTIGEN</div>')
        target_format = st.pills(
            "Target file format", ["PDB", "CIF"], selection_mode="single", default="PDB",
            key="rfab-target-format", label_visibility="collapsed",
        ) or "PDB"
        upload_extensions = ["pdb", "ent"] if target_format == "PDB" else ["cif", "mmcif"]
        previous_upload = st.session_state.get("rfab-target-upload")
        if previous_upload is not None and Path(previous_upload.name).suffix.lower() not in {f".{ext}" for ext in upload_extensions}:
            st.session_state.pop("rfab-target-upload", None)
        source = render_target_source(
            key_prefix="rfab", label="Target", example_code="7JZU", allow_cif=True,
            upload_file_types=upload_extensions,
            upload_expander_label=f"Upload a {target_format} file instead",
            upload_label=f"Drop a {target_format} file here or browse",
            upload_help=f"Protein target structure in {target_format} format.",
            upload_in_expander=False,
            fetched_path_key="rfab-example-target-path", fetched_name_key="rfab-example-target-name", owner_id=owner_id,
        )
        target_chains: list[str] = []
        hotspots = ""
        target_file: Path | None = None
        if source.path:
            target = render_target_editor(
                source, key_prefix="rfab-target-editor", label="Antigen", chain_selection_label="Antigen Chain IDs",
                residue_selection_title="Hotspot residues",
                residue_selection_caption="Please select antigen chains to enable residue selection. Choose hotspots in the MolView sequence panel or structure.",
                residue_selection_input_label="Hotspot residues",
                allow_trimming=False, show_structure_tools=False, viewer_height=480,
            )
            target_chains, hotspots = target.chains, target.hotspots
            target_file = target.design_path or source.path
        else:
            st.markdown("#### Antigen Chain IDs")
            st.caption("Chain selection becomes available after loading a PDB or mmCIF structure.")
            st.markdown("#### Hotspot residues")
            st.caption("Please upload a PDB file and select chains to enable residue selection.")

    custom_framework: Path | None = None
    heavy_chain = "H"
    light_chain = "L"
    if framework == "custom":
        with st.container(key="rfab-framework-card"):
            st.html('<div class="rfab-section">CUSTOM FRAMEWORK</div>')
            st.caption("Upload a Chothia-numbered framework PDB. RFantibody converts selected chains into its H/L/T format.")
            framework_upload = st.file_uploader("Framework PDB", type=["pdb", "ent"], key="rfab-framework-upload")
            if framework_upload is not None:
                try:
                    custom_framework = _save_custom_framework(framework_upload, owner_id)
                except OSError as exc:
                    st.error(f"Could not save custom framework: {exc}")
            heavy_chain = st.text_input("Heavy chain ID", "H", max_chars=1, key="rfab-heavy-chain").strip() or "H"
            if task == "antibody":
                light_chain = st.text_input("Light chain ID", "L", max_chars=1, key="rfab-light-chain").strip() or "L"

    st.html('<div class="rfab-section">REGIONS TO DESIGN</div>')
    available_regions = list(inputs.NANOBODY_REGIONS if task == "nanobody" else inputs.CDR_REGIONS)
    previous_regions = st.session_state.get("rfab-regions")
    if previous_regions is None or not set(previous_regions).issubset(available_regions):
        st.session_state["rfab-regions"] = available_regions
    regions = st.multiselect("CDR regions", available_regions, default=available_regions, key="rfab-regions", label_visibility="collapsed") or []
    st.caption("Choose which antibody loops RFdiffusion should design. Other framework regions stay fixed.")

    if task == "antibody":
        st.caption("Antibody mode designs selected heavy- and/or light-chain CDRs; nanobody mode uses heavy-chain CDRs.")
    error = None
    settings = {
        "task": task, "target_pdb": str(target_file) if target_file else None, "target_chains": target_chains,
        "hotspots": hotspots, "framework": framework, "framework_pdb": str(custom_framework) if custom_framework else None,
        "custom_heavy_chain": heavy_chain, "custom_light_chain": light_chain, "regions": list(regions), "loop_lengths": {},
        "job_name": job_name,
    }
    with st.container(key="rfab-design-count-card"):
        st.html('<div class="rfab-section">NUMBER OF DESIGNS</div>')
        settings["number_of_designs"] = int(st.number_input(
            "Number of designs", min_value=1, max_value=100000, value=1, step=1, key="rfab-count",
            help="Total number of antibody or nanobody designs to generate. Larger runs need more GPU time and storage.",
        ))
        st.caption("Generate from 1 to 100,000 designs. Use the Optional settings section for loop lengths and model parameters.")
    with st.expander("View Optional Settings", expanded=False):
        if regions:
            st.markdown("**CDR loop lengths**")
            st.caption("Use auto for Chothia-defined boundaries, or specify a fixed length or range (1–40 residues).")
            for region in regions:
                chain, number = region[0], region[1]
                settings["loop_lengths"][region] = st.text_input(
                    f"{chain}CDR{number} Length", "auto", key=f"rfab-loop-length-{region}",
                    placeholder="auto or 10-15",
                )
        else:
            st.info("Select at least one CDR region before configuring loop lengths.")

        st.markdown("**Sequence design**")
        settings["number_of_sequences"] = int(st.number_input("Sequences per structure", min_value=1, max_value=100, value=1, step=1, key="rfab-sequences"))
        settings["use_abmpnn"] = st.checkbox("Use AbMPNN weights", False, key="rfab-use-abmpnn", help="Use the installed antibody-specific ProteinMPNN checkpoint from Studio’s reference directory.")
        settings["mpnn_temperature"] = st.number_input("ProteinMPNN temperature", 0.01, 2.0, 0.1, 0.01, key="rfab-temperature")
        settings["num_connections"] = int(st.number_input("ProteinMPNN neighbors", 1, 100, 48, 1, key="rfab-connections"))
        settings["omit_aas"] = st.text_input("Omit amino acids", "CX", key="rfab-omit-aas", help="Use one-letter codes, for example CX.")
        st.caption("The pinned RFantibody image does not yet support global or per-residue amino-acid bias in its ProteinMPNN interface.")
        settings["bias_aas"] = st.text_input("Global amino-acid bias", "", key="rfab-bias-aas", placeholder="W:3.0,P:3.0,C:3.0,A:-3.0", disabled=True)
        settings["bias_by_res_json"] = st.text_area("Per-residue bias JSON", "{}", height=90, key="rfab-bias-by-res", help='Use antibody chain IDs and 1-based sequence positions, e.g. {"H":{"31":{"W":2.0}}}.', disabled=True)

        st.markdown("**Backbone and structure prediction**")
        settings["augment_eps"] = st.number_input("Backbone noise augmentation", 0.0, 1.0, 0.05, 0.01, key="rfab-augment")
        settings["num_recycles"] = int(st.number_input("RF2 recycles", 0, 20, 3, 1, key="rfab-recycles"))
        settings["rf2_repeats"] = int(st.number_input("Number of models (RF2)", 1, 20, 5, 1, key="rfab-repeats", help="Independent RF2 predictions with different random seeds."))
        settings["calculate_epitope_distance"] = st.checkbox("Calculate distance to epitope", True, key="rfab-distance")
        settings["hotspot_show_prop"] = st.number_input("Hotspot visibility", 0.0, 1.0, 0.1, 0.05, key="rfab-hotspot-prop")
        settings["save_trajectories"] = st.checkbox("Save diffusion trajectories", False, key="rfab-trajectories")

        st.markdown("**Compute**")
        settings["cpu_cores"] = int(st.number_input("CPU cores", min_value=1, max_value=256, value=4, step=1, key="rfab-cpus"))
        settings["gpu_device"] = st.selectbox("GPU device", ["0", "1", "all"], index=0, key="rfab-gpu", help="RFantibody requires an NVIDIA GPU.")
        prioritize = st.checkbox("Prioritize in scheduler queue", key="rfab-priority")
    try:
        inputs.validate_request(settings, require_files=bool(target_file))
        if not target_file:
            error = "Load a target structure before submitting."
    except (TypeError, ValueError, OSError) as exc:
        error = str(exc)
    if error:
        st.error(error)

    st.button("Load example inputs", icon=":material/science:", key="rfab-example", on_click=_load_example_inputs, args=(owner_id,))
    can_submit = bool(project and target_file and not error and regions)
    if st.button("Submit design", type="primary", use_container_width=True, disabled=not can_submit, key="rfab-submit"):
        try:
            inputs.validate_request(settings, require_files=True)
            submit(schema, "design.rfantibody", settings, project, {"prioritize": prioritize}, "RFantibody")
        except (TypeError, ValueError, OSError) as exc:
            st.error(str(exc))
    if not project:
        st.caption("Select a project to submit.")


def render_results(run_id: str, run_dir: Path, candidates: list[dict], parameters: dict, result: dict) -> None:
    st.subheader("RFantibody designs")
    metrics = result.get("metrics") if isinstance(result.get("metrics"), dict) else {}
    st.caption(f"{inputs.TASKS.get(parameters.get('task'), 'Antibody')} design · {metrics.get('candidate_count', len(candidates))} predicted structure(s) · RFantibody {inputs.UPSTREAM_COMMIT[:8]}")
    if not candidates:
        st.info("Completed RFantibody structure predictions will appear here.")
        return
    rows = []
    for rank, candidate in enumerate(candidates, 1):
        score = candidate.get("metrics") if isinstance(candidate.get("metrics"), dict) else {}
        rows.append({
            "Rank": rank,
            "Design": candidate.get("candidate_id"),
            "Sequence": candidate.get("binder_sequence"),
            "Chains": ",".join(candidate.get("binder_chains") or []),
            "Length": candidate.get("binder_length"),
            "pLDDT": score.get("pred_lddt"),
            "Interaction PAE": score.get("interaction_pae"),
            "RF2 model": score.get("rf2_repeat"),
        })
    frame = pd.DataFrame(rows)
    st.dataframe(frame, hide_index=True, use_container_width=True, height=min(410, 38 * (len(frame) + 1) + 3), column_config={"Sequence": st.column_config.TextColumn(width="large")})
    selected = st.selectbox("Design to visualize", range(len(candidates)), format_func=lambda index: f"{index + 1}. {candidates[index].get('candidate_id')}", key=f"rfab-result-{run_id}")
    candidate = candidates[selected]
    structure_path = Path(str(candidate.get("complex_pdb") or ""))
    binder_path = Path(str(candidate.get("binder_pdb") or ""))
    if structure_path.is_file():
        molstar_custom_component(
            [StructureVisualization(pdb=structure_path.read_text(errors="ignore"), color="molecule-type", representation_type="cartoon+ball-and-stick")],
            key=f"rfab-result-view-{run_id}-{selected}", height=560, show_controls=True, show_structure_tools=False,
        )
        download_cols = st.columns(2)
        download_cols[0].download_button("Download predicted complex", structure_path.read_bytes(), file_name=structure_path.name, key=f"rfab-complex-download-{run_id}-{selected}")
        if binder_path.is_file():
            download_cols[1].download_button("Download antibody", binder_path.read_bytes(), file_name=binder_path.name, key=f"rfab-binder-download-{run_id}-{selected}")
    st.markdown("**Designed antibody sequence**")
    st.code(candidate.get("binder_sequence") or "", language=None)
    score = candidate.get("metrics") or {}
    shown_metrics = {key.replace("_", " ").title(): value for key, value in score.items() if isinstance(value, (int, float, str, bool)) and key not in {"binder_ca_trace_c3", "binder_ca_trace_segments"}}
    with st.expander("Structure scores and design metadata", expanded=False):
        st.json(shown_metrics)
        st.json(candidate.get("raw_metadata") or {})
    with st.expander("Requested settings", expanded=False):
        st.json({key: value for key, value in parameters.items() if key not in {"target_pdb", "framework_pdb"}})
