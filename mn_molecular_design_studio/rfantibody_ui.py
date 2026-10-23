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
from mn_molecular_design_studio.target_editor import _hash, _target_widget_keys, render_target_editor, render_target_source


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


def _save_mimosa_complex(upload, owner_id: str) -> Path | None:
    if upload is None:
        return None
    content = upload.getvalue()
    digest = hashlib.sha256(content).hexdigest()[:16]
    suffix = Path(upload.name).suffix.lower()
    destination = user_input_root(owner_id) / "rfantibody" / f"mimosa-complex-{digest}{suffix}"
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


def _load_mimosa_paper_example(owner_id: str, example_id: str) -> None:
    from mn_molecular_design_studio.mimosa_ui import (
        contacting_partner_residues,
        fetch_paper_mimosa_example,
    )

    example_id = str(example_id).strip().lower()
    example = fetch_paper_mimosa_example(
        example_id, user_input_root(owner_id) / "rfantibody" / "mimosa" / example_id
    )
    target_path = Path(example["target_path"])
    source_path = Path(example["source_path"])
    target_chain = str(example["target_chain"])
    partner_chain = str(example["source_partner_chain"])
    source_target_chain = str(example["source_target_chain"])
    target_contacts = contacting_partner_residues(
        source_path, source_target_chain, partner_chain, contact_cutoff=5.0
    )

    st.session_state["rfab-example-target-path"] = str(target_path)
    st.session_state["rfab-example-target-name"] = str(example["target_name"])
    st.session_state["rfab-mimosa-job-name"] = f"mimosa-{example_id.replace('_', '-')}-rfantibody"
    st.session_state["rfab-target-format"] = "PDB"
    st.session_state.pop("rfab-target-upload", None)
    st.session_state.pop("rfab-mimosa-complex", None)
    st.session_state[_target_widget_keys("rfab-target-editor", target_path)["viewer_mode"]] = "Whole complex"
    st.session_state[f"rfab-target-editor-chains-{_hash(str(target_path))}"] = [target_chain]
    st.session_state[_target_widget_keys("rfab-target-editor", target_path)["hotspots"]] = ",".join(
        f"{target_chain}{row['position']}" for row in target_contacts
    )
    design_loop = str(example["design_loop"])
    st.session_state["rfab-mimosa-regions"] = list(
        example.get("rfantibody_regions")
        or (["H1", "H2", "H3"] if design_loop.lower() == "all_cdrs"
            else [design_loop.replace("CDR", "H")])
    )
    for key, value in {
        "rfab-mimosa-scale-v2": 0.3,
        "rfab-mimosa-particles-v2": 4,
        "rfab-mimosa-fk-lambda-v2": 4.0,
        "rfab-mimosa-resampling-interval-v2": 3,
    }.items():
        st.session_state[key] = value
    for region, length in (("H1", "7"), ("H2", "6"), ("H3", "5-17")):
        st.session_state[f"rfab-loop-length-{region}"] = length
    st.session_state["rfab-mimosa-loaded-source"] = str(source_path)
    st.session_state["rfab-mimosa-paper-template"] = str(example["source_template"])
    st.session_state["rfab-mimosa-condition-only"] = str(example["condition"])
    st.session_state["rfab-mimosa-binder-chain"] = [
        chain.strip() for chain in partner_chain.split(",") if chain.strip()
    ]
    st.session_state["rfab-mimosa-source-target-chain"] = source_target_chain
    st.session_state["rfab-mimosa-design-target-chain"] = target_chain
    st.session_state["rfab-mimosa-design-loop"] = (
        "all_cdrs" if design_loop.lower() == "all_cdrs"
        else design_loop.replace("CDR", "H")
    )
    st.session_state["rfab-mimosa-contact-cutoff"] = 8.0
    st.session_state["rfab-mimosa-source-positions"] = list(example["motif_positions"])
    st.session_state["rfab-mimosa-source-selection-context"] = (
        str(source_path.resolve()), partner_chain, source_target_chain, 8.0
    )


def _load_mimosa_paper_example_from_button(owner_id: str, example_id: str) -> None:
    try:
        _load_mimosa_paper_example(owner_id, example_id)
        st.session_state.pop("rfab-mimosa-example-load-error", None)
    except (OSError, ValueError) as exc:
        st.session_state["rfab-mimosa-example-load-error"] = str(exc)


def render_setup(schema: dict, *, owner_id: str, project_selector, submit, mimosa_only: bool = False) -> None:
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
        _render_setup_content(schema, owner_id=owner_id, project_selector=project_selector, submit=submit, mimosa_only=mimosa_only)


def _render_setup_content(schema: dict, *, owner_id: str, project_selector, submit, mimosa_only: bool = False) -> None:
    if st.button("← Back to All tools", key="rfab-back"):
        st.session_state["page"] = "All tools"
        st.session_state.pop("selected_workflow", None)
        st.rerun()
    page_title = "RFantibody · MIMOSA" if mimosa_only else "RFantibody"
    st.html(f'<div class="rfab-heading"><div class="rfab-eyebrow">ANTIBODY DESIGN</div><h1>{page_title}</h1>'
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
        default_job_name = "mimosa-keap1-nrf2-rfantibody" if mimosa_only else "rfantibody-design"
        job_name_key = "rfab-mimosa-job-name" if mimosa_only else "rfab-job-name"
        job_name = st.text_input("Job name", default_job_name, key=job_name_key).strip() or default_job_name
    with right:
        project = project_selector("rfantibody")

    st.html('<div class="rfab-section">DESIGN MODE</div>')
    if mimosa_only:
        task = "nanobody"
        st.info("The paper’s MIMOSA validation starts with VHH nanobody design.")
    else:
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
        st.session_state.setdefault("rfab-target-format", "PDB")
        target_format = st.pills(
            "Target file format", ["PDB", "CIF"], selection_mode="single",
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
                residue_selection_caption=(
                    "Optional target-site hotspots that guide RFdiffusion toward the epitope. These are separate from the source-partner motif residues selected below."
                    if mimosa_only
                    else "Please select antigen chains to enable residue selection. Choose hotspots in the MolView sequence panel or structure."
                ),
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
    regions_key = "rfab-mimosa-regions" if mimosa_only else "rfab-regions"
    previous_regions = st.session_state.get(regions_key)
    if previous_regions is None or not set(previous_regions).issubset(available_regions):
        st.session_state[regions_key] = available_regions
    regions = st.multiselect("CDR regions", available_regions, key=regions_key, label_visibility="collapsed") or []
    st.caption("Choose which antibody loops RFdiffusion should design. Other framework regions stay fixed.")

    settings: dict = {}
    with st.expander("MIMOSA: guide a cognate interface motif", expanded=mimosa_only):
        st.caption(
            "Use a source protein–target complex to select the partner amino acids whose interface geometry to mimic. "
            "The partner can be an antibody, receptor, or other protein. MIMOSA dynamically matches those residues "
            "to the generated loop during denoising; only motif residues placed within 2 Å keep their source identity."
        )
        settings["mimosa_condition"] = "unconstrained"
        if mimosa_only:
            st.session_state.setdefault("rfab-mimosa-condition-only", "sequence_fixed")
            condition = st.selectbox(
                "MIMOSA condition",
                ["geometry_only", "sequence_fixed", "unconstrained"],
                format_func={
                    "geometry_only": "Geometry-only · guide motif geometry; redesign all CDR identities",
                    "sequence_fixed": "Sequence-fixed · guide geometry; keep placed motif identities",
                    "unconstrained": "Unconstrained · BoltzGen/RFantibody baseline without motif guidance",
                }.get,
                key="rfab-mimosa-condition-only",
            )
            settings["mimosa_condition"] = condition
            settings["mimosa_enabled"] = condition != "unconstrained"
        else:
            settings["mimosa_enabled"] = st.checkbox(
                "Enable MIMOSA motif guidance", key="rfab-mimosa-enabled"
            )
            if settings["mimosa_enabled"]:
                condition = st.selectbox(
                    "MIMOSA condition",
                    ["geometry_only", "sequence_fixed"],
                    index=1,
                    format_func={
                        "geometry_only": "Geometry-only · redesign all CDR identities",
                        "sequence_fixed": "Sequence-fixed · keep placed motif identities",
                    }.get,
                    key="rfab-mimosa-condition",
                )
                settings["mimosa_condition"] = condition
        settings["mimosa_complex"] = None
        settings["mimosa_binder_chain"] = ""
        settings["mimosa_complex_target_chain"] = ""
        settings["mimosa_target_chain"] = target_chains[0] if target_chains else ""
        settings["mimosa_cdr"] = "CDR3"
        settings["mimosa_design_loop"] = "H3" if task == "nanobody" else "H3"
        settings["mimosa_source_positions"] = ""
        settings["mimosa_contact_cutoff"] = 8.0
        settings["mimosa_guidance_scale"] = 0.3
        settings["mimosa_num_particles"] = 4
        settings["mimosa_fk_lambda"] = 4.0
        settings["mimosa_resampling_interval"] = 3
        if settings["mimosa_enabled"]:
            complex_upload = st.file_uploader(
                "Source interaction complex", type=["pdb", "ent", "cif", "mmcif"],
                key="rfab-mimosa-complex",
                help="A known interaction structure, such as the complex whose motif you want the designed nanobody to mimic.",
            )
            from mn_molecular_design_studio.mimosa_ui import (
                PAPER_INTERACTION_TEMPLATES,
                fetch_paper_interaction_template,
                render_source_motif_viewer,
            )

            st.session_state.setdefault("rfab-mimosa-paper-template", "2FLU")
            paper_template = st.selectbox(
                "Or use a paper interaction template",
                ["custom", *PAPER_INTERACTION_TEMPLATES],
                format_func=lambda value: "Upload my own" if value == "custom" else f"{value} · {PAPER_INTERACTION_TEMPLATES[value]}",
                key="rfab-mimosa-paper-template",
            )
            st.caption("These PDB entries are interaction templates from the paper; load the design target separately above.")
            if paper_template != "custom" and st.button("Load interaction template", key="rfab-mimosa-load-template"):
                template_path = user_input_root(owner_id) / "mimosa" / f"{paper_template}.cif"
                try:
                    loaded_path = fetch_paper_interaction_template(paper_template, template_path)
                    st.session_state["rfab-mimosa-loaded-source"] = str(loaded_path)
                    st.success(f"Loaded {paper_template} interaction complex.")
                except (OSError, ValueError) as exc:
                    st.error(str(exc))
            try:
                uploaded_path = _save_mimosa_complex(complex_upload, owner_id)
                complex_path = uploaded_path or (
                    Path(st.session_state["rfab-mimosa-loaded-source"])
                    if st.session_state.get("rfab-mimosa-loaded-source")
                    else None
                )
                # Workflow payloads cross the Studio JSON API boundary. Keep the
                # selected source file JSON-safe while downstream validation can
                # still normalize it back to Path.
                settings["mimosa_complex"] = str(complex_path) if complex_path else None
            except OSError as exc:
                st.error(f"Could not save the MIMOSA source complex: {exc}")
            settings["mimosa_contact_cutoff"] = st.number_input(
                "Source interface contact cutoff (Å)", min_value=1.0, max_value=20.0, value=8.0, step=0.5,
                key="rfab-mimosa-cutoff",
                help="A source partner amino acid is offered below when any heavy atom is within this distance of the source target.",
            )
            if settings["mimosa_complex"]:
                from mn_molecular_design_studio.mimosa_ui import contacting_partner_residues, source_complex_chains

                try:
                    source_chain_ids = source_complex_chains(settings["mimosa_complex"])
                    if len(source_chain_ids) < 2:
                        raise ValueError("The source complex needs at least two protein chains.")
                    partner_key = "rfab-mimosa-binder-chain"
                    stored_partner_chains = st.session_state.get(partner_key)
                    if isinstance(stored_partner_chains, str):
                        stored_partner_chains = [item.strip() for item in stored_partner_chains.split(",") if item.strip()]
                    stored_partner_chains = [chain for chain in (stored_partner_chains or []) if chain in source_chain_ids]
                    if not stored_partner_chains:
                        stored_partner_chains = ["H" if "H" in source_chain_ids else source_chain_ids[0]]
                    st.session_state[partner_key] = stored_partner_chains
                    selected_partner_chains = st.multiselect(
                        "Source partner chain(s)", source_chain_ids, key=partner_key,
                        help="Select multiple partner chains when the source interaction motif spans them, as in the 6WZM antibody example.",
                    )
                    settings["mimosa_binder_chain"] = ",".join(selected_partner_chains)
                    source_target_options = [chain for chain in source_chain_ids if chain not in selected_partner_chains]
                    if not source_target_options:
                        raise ValueError("Choose at least one source partner chain and leave a different chain as the source target.")
                    source_target_key = "rfab-mimosa-source-target-chain"
                    if st.session_state.get(source_target_key) not in source_target_options:
                        st.session_state[source_target_key] = "A" if "A" in source_target_options else source_target_options[0]
                    settings["mimosa_complex_target_chain"] = st.selectbox(
                        "Source target chain", source_target_options, key=source_target_key
                    )
                    if target_chains:
                        settings["mimosa_target_chain"] = st.selectbox(
                            "Design target chain", target_chains, key="rfab-mimosa-design-target-chain"
                        )
                    residues = contacting_partner_residues(
                        settings["mimosa_complex"], settings["mimosa_binder_chain"],
                        settings["mimosa_complex_target_chain"], settings["mimosa_contact_cutoff"],
                    )
                    position_options = [row["position"] for row in residues]
                    labels = {row["position"]: f"{row['label']} · {row['distance']}" for row in residues}
                    if not position_options:
                        raise ValueError("No contacting amino acids were found on the selected source partner chain(s) at this cutoff.")
                    selection_key = "rfab-mimosa-source-positions"
                    selection_context_key = "rfab-mimosa-source-selection-context"
                    selection_context = (
                        str(Path(settings["mimosa_complex"]).resolve()),
                        settings["mimosa_binder_chain"],
                        settings["mimosa_complex_target_chain"],
                        float(settings["mimosa_contact_cutoff"]),
                    )
                    if st.session_state.get(selection_context_key) != selection_context:
                        st.session_state[selection_key] = position_options
                        st.session_state[selection_context_key] = selection_context
                    else:
                        old_selection = st.session_state.get(selection_key, position_options)
                        st.session_state[selection_key] = [value for value in old_selection if value in position_options]
                    settings["mimosa_source_positions"] = st.multiselect(
                        "Source partner amino acids to keep as the motif",
                        options=position_options,
                        format_func=lambda value: labels.get(value, value),
                        key=selection_key,
                        help="These selected contacting residues supply the motif coordinates and identities. Clear unwanted residues; select at least one.",
                    )
                    st.caption(f"{len(settings['mimosa_source_positions'])} of {len(position_options)} contacting source residues selected.")
                    render_source_motif_viewer(
                        settings["mimosa_complex"],
                        settings["mimosa_binder_chain"],
                        settings["mimosa_source_positions"],
                        key_prefix="rfab-mimosa",
                    )
                except (OSError, ValueError, KeyError) as exc:
                    st.error(f"Could not inspect the source interaction complex: {exc}")
            else:
                st.info("Upload a source interaction complex to select the partner residues that form the motif.")
            loop_col = st.container()
            with loop_col:
                default_loop = "H3" if task == "nanobody" or "H3" in regions else (regions[0] if regions else "H3")
                loop_options = list(regions) or [default_loop]
                if task == "nanobody" and {"H1", "H2", "H3"}.issubset(regions):
                    loop_options.append("all_cdrs")
                loop_key = "rfab-mimosa-design-loop"
                if st.session_state.get(loop_key) not in loop_options:
                    st.session_state[loop_key] = default_loop if default_loop in loop_options else loop_options[0]
                settings["mimosa_design_loop"] = st.selectbox(
                    "Generated CDR to guide", loop_options,
                    format_func=lambda value: "All three CDRs · scattered motif assignment" if value == "all_cdrs" else value,
                    key=loop_key,
                )
            defaults = {
                "rfab-mimosa-scale-v2": 0.3,
                "rfab-mimosa-particles-v2": 4,
                "rfab-mimosa-fk-lambda-v2": 4.0,
                "rfab-mimosa-resampling-interval-v2": 3,
            }
            for key, value in defaults.items():
                st.session_state.setdefault(key, value)
            if st.button("Restore paper defaults", key="rfab-mimosa-restore-defaults"):
                for key, value in defaults.items():
                    st.session_state[key] = value
            c2, c3 = st.columns(2)
            with c2:
                settings["mimosa_guidance_scale"] = st.number_input(
                    "Maximum coordinate guidance", min_value=0.0, max_value=2.0, step=0.05,
                    key="rfab-mimosa-scale-v2",
                )
            with c3:
                settings["mimosa_num_particles"] = int(st.number_input(
                    "Particle trajectories", min_value=1, max_value=64, step=1,
                    key="rfab-mimosa-particles-v2",
                ))
            c4, c5 = st.columns(2)
            with c4:
                settings["mimosa_fk_lambda"] = st.number_input(
                    "Feynman–Kac lambda", min_value=0.0, max_value=20.0, step=0.25,
                    key="rfab-mimosa-fk-lambda-v2",
                )
            with c5:
                settings["mimosa_resampling_interval"] = int(st.number_input(
                    "Resampling interval (steps)", min_value=1, max_value=50, step=1,
                    key="rfab-mimosa-resampling-interval-v2",
                ))
            st.caption("Final sequence identities are fixed only for assigned motif residues within 2 Å of their source Cα positions.")
            st.caption("MIMOSA and SAbR run in Studio's local runtime image.")

    if task == "antibody":
        st.caption("Antibody mode designs selected heavy- and/or light-chain CDRs; nanobody mode uses heavy-chain CDRs.")
    error = None
    settings.update({
        "task": task, "target_pdb": str(target_file) if target_file else None, "target_chains": target_chains,
        "hotspots": hotspots, "framework": framework, "framework_pdb": str(custom_framework) if custom_framework else None,
        "custom_heavy_chain": heavy_chain, "custom_light_chain": light_chain, "regions": list(regions), "loop_lengths": {},
        "job_name": job_name,
    })
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

    if mimosa_only:
        from mn_molecular_design_studio.mimosa_ui import PAPER_MIMOSA_EXAMPLES

        example_ids = list(PAPER_MIMOSA_EXAMPLES)
        example_key = "rfab-mimosa-paper-example"
        st.session_state.setdefault(example_key, "keap1_nrf2")
        selected_example = st.selectbox(
            "Paper interaction example",
            example_ids,
            format_func=lambda value: PAPER_MIMOSA_EXAMPLES[value]["label"],
            key=example_key,
        )
        st.button(
            "Load selected paper example into all fields",
            icon=":material/science:",
            key="rfab-mimosa-load-paper-example",
            on_click=_load_mimosa_paper_example_from_button,
            args=(owner_id, selected_example),
        )
        example_error = st.session_state.pop("rfab-mimosa-example-load-error", None)
        if example_error:
            st.error(f"Could not load the MIMOSA paper example: {example_error}")
        st.caption(
            "Fig. 2 offers contiguous single-loop examples (2FLU, 3PB1); Fig. 3 offers scattered IL-8 motifs from CXCR1 (8IC0) or LY3041658 heavy and light chains (6WZM). "
            "The loader fills the target, source chains, motif residues, target contact hotspots, and MolView structure."
        )
    else:
        st.button("Load example inputs", icon=":material/science:", key="rfab-example", on_click=_load_example_inputs, args=(owner_id,))
    can_submit = bool(project and target_file and not error and regions)
    if st.button("Submit design", type="primary", use_container_width=True, disabled=not can_submit, key="rfab-submit"):
        try:
            inputs.validate_request(settings, require_files=True)
            selected_workflow = "design.rfantibody_mimosa" if mimosa_only else "design.rfantibody"
            submit(schema, selected_workflow, settings, project, {"prioritize": prioritize}, page_title)
        except (TypeError, ValueError, OSError) as exc:
            st.error(str(exc))
    if not project:
        st.caption("Select a project to submit.")


def render_results(run_id: str, run_dir: Path, candidates: list[dict], parameters: dict, result: dict) -> None:
    st.subheader("RFantibody designs")
    metrics = result.get("metrics") if isinstance(result.get("metrics"), dict) else {}
    generated_count = sum(candidate.get("stage") == "generation.backbone" for candidate in candidates)
    rf2_count = sum(
        candidate.get("stage") == "complex_refolding"
        or (candidate.get("metrics") or {}).get("rf2_repeat") is not None
        for candidate in candidates
    )
    stage_counts = []
    if generated_count:
        stage_counts.append(f"{generated_count} generated backbone" + ("s" if generated_count != 1 else ""))
    if rf2_count:
        stage_counts.append(f"{rf2_count} RF2 prediction" + ("s" if rf2_count != 1 else ""))
    summary = " · ".join(stage_counts) if stage_counts else f"{metrics.get('candidate_count', len(candidates))} predicted structure(s)"
    st.caption(f"{inputs.TASKS.get(parameters.get('task'), 'Antibody')} design · {summary} · RFantibody {inputs.UPSTREAM_COMMIT[:8]}")
    if not candidates:
        st.info("Completed RFantibody structure predictions will appear here.")
        return
    rows = []
    rf2_indices = [
        index for index, candidate in enumerate(candidates)
        if candidate.get("stage") == "complex_refolding"
        or (candidate.get("metrics") or {}).get("rf2_repeat") is not None
    ]
    rf2_indices.sort(key=lambda index: (
        -float((candidates[index].get("metrics") or {}).get("pred_lddt") or 0.0),
        float((candidates[index].get("metrics") or {}).get("interaction_pae") or float("inf")),
    ))
    rf2_rank = {candidate_index: rank for rank, candidate_index in enumerate(rf2_indices, 1)}
    for rank, candidate in enumerate(candidates, 1):
        score = candidate.get("metrics") if isinstance(candidate.get("metrics"), dict) else {}
        stage = candidate.get("stage")
        stage_label = {
            "generation.backbone": "Generated backbone",
            "complex_refolding": "RF2 prediction",
        }.get(stage, stage or "Prediction")
        rows.append({
            "Rank": rf2_rank.get(rank - 1, "—"),
            "Design": candidate.get("candidate_id"),
            "Stage": stage_label,
            "Sequence": candidate.get("binder_sequence") or "—",
            "Chains": ",".join(candidate.get("binder_chains") or []),
            "Length": candidate.get("binder_length") or "—",
            "pLDDT": score.get("pred_lddt", "—"),
            "Interaction PAE": score.get("interaction_pae", "—"),
            "RF2 model": score.get("rf2_repeat", "—"),
        })
    frame = pd.DataFrame(rows)
    st.dataframe(frame, hide_index=True, use_container_width=True, height=min(410, 38 * (len(frame) + 1) + 3), column_config={"Sequence": st.column_config.TextColumn(width="large")})
    default_index = rf2_indices[0] if rf2_indices else 0
    selected = st.selectbox(
        "Design to visualize",
        range(len(candidates)),
        index=default_index,
        format_func=lambda index: (
            f"{rf2_rank[index]}. RF2 · {candidates[index].get('candidate_id')}"
            if index in rf2_rank else f"Generated backbone · {candidates[index].get('candidate_id')}"
        ),
        key=f"rfab-result-{run_id}",
    )
    candidate = candidates[selected]
    structure_path = Path(str(candidate.get("complex_pdb") or ""))
    binder_path = Path(str(candidate.get("binder_pdb") or ""))
    if structure_path.is_file():
        molstar_custom_component(
            [StructureVisualization(pdb=structure_path.read_text(errors="ignore"), color="chain-id", representation_type="cartoon")],
            key=f"rfab-result-view-{run_id}-{selected}", height=560, show_controls=True, show_structure_tools=False,
        )
        download_cols = st.columns(2)
        download_label = "Download generated backbone" if candidate.get("stage") == "generation.backbone" else "Download RF2 complex"
        download_cols[0].download_button(download_label, structure_path.read_bytes(), file_name=structure_path.name, key=f"rfab-complex-download-{run_id}-{selected}")
        if binder_path.is_file():
            download_cols[1].download_button("Download antibody", binder_path.read_bytes(), file_name=binder_path.name, key=f"rfab-binder-download-{run_id}-{selected}")
    st.markdown("**Designed antibody sequence**")
    if candidate.get("binder_sequence"):
        st.code(candidate["binder_sequence"], language=None)
    elif candidate.get("stage") == "generation.backbone":
        st.caption("No sequence has been assigned to this generated backbone.")
    else:
        st.code("", language=None)
    score = candidate.get("metrics") or {}
    shown_metrics = {key.replace("_", " ").title(): value for key, value in score.items() if isinstance(value, (int, float, str, bool)) and key not in {"binder_ca_trace_c3", "binder_ca_trace_segments"}}
    with st.expander("Structure scores and design metadata", expanded=False):
        st.json(shown_metrics)
        st.json(candidate.get("raw_metadata") or {})
    with st.expander("Requested settings", expanded=False):
        st.json({key: value for key, value in parameters.items() if key not in {"target_pdb", "framework_pdb"}})
    render_confidence_guide()


def render_confidence_guide() -> None:
    st.markdown("#### How to interpret these scores")
    st.markdown(
        "RFantibody's maintainers suggest **RF2 pAE <10 Å** as a minimal filter and "
        "**design-to-prediction RMSD <2 Å**; **Rosetta ddG <−20** may also help when available. "
        "They caution that current filters are not yet reliable at separating binders from non-binders. "
        "**Interaction PAE/ipAE** (Å) is lower when the predicted interface placement is more confident, "
        "but it is not the same metric as whole-complex pAE and RFantibody's guidance does not establish "
        "a validated ipAE cutoff. **pLDDT** measures local structure confidence, not binding probability; "
        "RFantibody reports it on a 0–1 scale here, so 0.90 corresponds to 90/100. No antibody-specific "
        "pLDDT cutoff is established by that guidance. Treat these as prioritization signals, not proof of binding."
    )
    st.caption(
        "Sources: [RFantibody filtering guidance](https://github.com/RosettaCommons/RFantibody#filtering-strategies) · "
        "[AlphaFold confidence guidance](https://alphafold.ebi.ac.uk/faq)"
    )
