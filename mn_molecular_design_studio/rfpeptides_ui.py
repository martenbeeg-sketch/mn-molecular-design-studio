"""Tamarind-inspired standalone RFpeptides input and result pages."""
from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd
import streamlit as st

from mn_molecular_design_studio import rfpeptides_input as inputs
from mn_molecular_design_studio.auth import user_input_root
from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer import StructureVisualization, molstar_custom_component
from mn_molecular_design_studio.target_editor import _target_widget_keys, render_target_editor, render_target_source

EXAMPLE_TARGET = Path(__file__).resolve().parent / "resources" / "rfpeptides" / "7zkr_GABARAP.pdb"


def _load_example_inputs(owner_id: str) -> None:
    destination = user_input_root(owner_id) / "rfpeptides" / "7zkr_GABARAP.pdb"
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(EXAMPLE_TARGET, destination)
    st.session_state["rfpeptides-target-format"] = "PDB"
    st.session_state["rfpeptides-example-target-path"] = str(destination)
    st.session_state["rfpeptides-example-target-name"] = "7zkr_GABARAP"
    hotspot_key = _target_widget_keys("rfpeptides-target-editor", destination)["hotspots"]
    st.session_state[hotspot_key] = "A51,A52,A50,A48,A62,A65"
    st.session_state.pop("rfpeptides-target-upload", None)


def render_setup(schema: dict, *, owner_id: str, project_selector, submit) -> None:
    st.html("""<style>
    .block-container {max-width:1240px!important;padding-top:1.5rem!important}
    .st-key-rfpeptides-page {--rfp-green:#2e6a50;--rfp-ink:#173d31;--rfp-mint:#eef6f1}
    .st-key-rfpeptides-page .rfp-heading {text-align:center;padding:.45rem 0 1.3rem}
    .st-key-rfpeptides-page .rfp-eyebrow {font-size:.69rem;letter-spacing:.13em;font-weight:650;color:#688575;margin-bottom:.6rem}
    .st-key-rfpeptides-page .rfp-heading h1 {font-size:2.1rem;color:#173d31;margin:0;font-weight:700;letter-spacing:-.04em}
    .st-key-rfpeptides-page .rfp-heading p {font-size:.92rem;color:#6c8276;margin:.55rem auto .45rem;max-width:620px;line-height:1.6}
    .st-key-rfpeptides-page .rfp-meta {font-size:.73rem;color:#789084}
    .st-key-rfp-target-card {padding:1.2rem 1.25rem!important;background:#fff;border:1px solid #d7e6dc;border-radius:12px;margin:.45rem 0 .75rem}
    .st-key-rfp-count-card {padding:1.05rem 1.25rem!important;background:#f3f8f5;border:1px solid #d7e6dc;border-radius:12px;margin:.35rem 0 .65rem}
    .st-key-rfpeptides-page .rfp-section {font-size:.72rem;font-weight:650;letter-spacing:.1em;color:#5c7b68;margin:0 0 .7rem}
    .st-key-rfpeptides-page [data-variant="pills"] {border-radius:7px!important;border:1px solid #c2d9cc!important;background:#f7faf8;color:#38614e!important;min-height:2.55rem;padding:.4rem .85rem!important;font-size:.84rem!important;box-shadow:none!important}
    .st-key-rfpeptides-page [data-variant="pills"][aria-checked="true"] {background:#2e6a50!important;color:white!important;border-color:#2e6a50!important}
    .st-key-rfpeptides-page [data-testid="stTextInput"] input,.st-key-rfpeptides-page [data-testid="stNumberInput"] input,.st-key-rfpeptides-page textarea {background:#f5f9f6!important;color:#193e30!important}
    .st-key-rfpeptides-page [data-baseweb="input"],.st-key-rfpeptides-page [data-baseweb="textarea"] {border-color:#dce7df!important;border-radius:7px!important;background:#f5f9f6!important}
    .st-key-rfpeptides-page [data-baseweb="tag"] {background:#e0eee5!important;color:#27674e!important}
    .st-key-rfpeptides-page [data-testid="stButton"] button {border-radius:7px!important;color:#27674e;border-color:#bad2c3}
    .st-key-rfpeptides-page [data-testid="stButton"] button:hover {background:#edf6f0!important;border-color:#659d7d!important;color:#1d5039!important}
    .st-key-rfpeptides-page button[kind="primary"] {background:#2e6a50!important;border-color:#2e6a50!important;color:#fff!important;min-height:2.8rem}
    .st-key-rfpeptides-page [data-testid="stExpander"] {border-color:#d7e6dc!important;border-radius:8px!important;background:#fcfefc}
    .st-key-rfpeptides-page [data-testid="stFileUploaderDropzone"] {background:#f3f8f5!important;border:1px dashed #b9d2c2!important;border-radius:8px}
    .st-key-rfpeptides-page [data-testid="stCaptionContainer"] {color:#7c8f82}
    </style>""")
    with st.container(key="rfpeptides-page"):
        _render_setup_content(schema, owner_id=owner_id, project_selector=project_selector, submit=submit)


def _render_setup_content(schema: dict, *, owner_id: str, project_selector, submit) -> None:
    if st.button("← Back to All tools", key="rfpeptides-back"):
        st.session_state["page"] = "All tools"
        st.session_state.pop("selected_workflow", None)
        st.rerun()
    st.html(
        '<div class="rfp-heading"><div class="rfp-eyebrow">PROTEIN DESIGN</div><h1>RFpeptides</h1>'
        '<p>Design macrocyclic peptide binder backbones for a selected protein surface.</p>'
        f'<span class="rfp-meta">{inputs.CODE_VERSION} · GPU required</span></div>'
    )
    with st.expander("About RFpeptides", expanded=False):
        st.write("This standalone tool runs the RFpeptides macrocyclic binder protocol from RFdiffusion.")
        st.markdown(
            "[RFdiffusion source](https://github.com/RosettaCommons/RFdiffusion) · "
            "[Tamarind input reference](https://app.tamarind.bio/tools/rfpeptides) · "
            "[RFpeptides paper](https://www.biorxiv.org/content/10.1101/2024.11.18.622547v1)"
        )
        st.caption("This run generates cyclic peptide backbones. It does not design peptide sequences or use Rosetta; use ProteinMPNN or another sequence-design tool next.")

    job_col, project_col = st.columns([1.35, 1], vertical_alignment="bottom")
    with job_col:
        job_name = st.text_input("Job Name", "rfpeptides-design", key="rfpeptides-job-name").strip() or "rfpeptides-design"
    with project_col:
        project = project_selector("rfpeptides")

    with st.container(key="rfp-target-card"):
        st.html('<div class="rfp-section">TARGET</div>')
        target_format = st.pills(
            "PDB File", ["PDB", "CIF"], selection_mode="single", default="PDB",
            key="rfpeptides-target-format", label_visibility="collapsed",
        ) or "PDB"
        upload_extensions = ["pdb", "ent"] if target_format == "PDB" else ["cif", "mmcif"]
        previous_upload = st.session_state.get("rfpeptides-target-upload")
        if previous_upload is not None and Path(previous_upload.name).suffix.lower() not in {f".{ext}" for ext in upload_extensions}:
            st.session_state.pop("rfpeptides-target-upload", None)
        source = render_target_source(
            key_prefix="rfpeptides", label="Target", example_code="7ZKR", allow_cif=True,
            upload_file_types=upload_extensions,
            upload_expander_label=f"Upload a {target_format} file instead",
            upload_label=f"Drop a {target_format} file here or browse",
            upload_help=f"Protein target structure in {target_format} format.",
            upload_in_expander=False,
            fetched_path_key="rfpeptides-example-target-path", fetched_name_key="rfpeptides-example-target-name",
            owner_id=owner_id,
        )
        target_chains: list[str] = []
        hotspots = ""
        target_file: Path | None = None
        binder_length = "12-18"
        if source.path:
            length_value: dict[str, str] = {}

            def render_binder_length(_target) -> None:
                length_value["value"] = st.text_input(
                    "Binder Length", "12-18", key="rfpeptides-binder-length",
                    help="Length of the cyclic peptide region to design. Enter one length or a range sampled uniformly, such as 12-18.",
                )

            target = render_target_editor(
                source, key_prefix="rfpeptides-target-editor", label="Target",
                chain_selection_label="Target Chains", residue_selection_title="Hotspots",
                residue_selection_caption="Select target residues that the cyclic peptide should contact. Choose residues in the MolView structure or enter chain IDs and residue numbers.",
                residue_selection_input_label="Hotspots", allow_trimming=False,
                show_structure_tools=False, viewer_height=480,
                before_residue_selection=render_binder_length,
            )
            target_chains, hotspots = target.chains, target.hotspots
            target_file = target.design_path or source.path
            binder_length = length_value.get("value", "12-18")
        else:
            binder_length = st.text_input(
                "Binder Length", "12-18", key="rfpeptides-binder-length",
                help="Length of the cyclic peptide region to design. Enter one length or a range sampled uniformly, such as 12-18.",
            )
            st.markdown("#### Target Chains")
            st.caption("Load a PDB or CIF structure to select target chains.")
            st.markdown("#### Hotspots")
            st.caption("Please upload a structure and select chains to enable residue selection.")

    error = None
    settings = {
        "target_pdb": str(target_file) if target_file else None,
        "target_chains": target_chains,
        "hotspots": hotspots,
        "binder_length": binder_length,
        "job_name": job_name,
        "cpu_cores": 4,
        "gpu_device": "0",
    }
    with st.container(key="rfp-count-card"):
        st.html('<div class="rfp-section">NUMBER OF DESIGNS</div>')
        settings["number_of_designs"] = int(st.number_input(
            "Number of Designs", min_value=1, max_value=100000, value=1, step=1,
            key="rfpeptides-count", help="Total number of cyclic peptide backbone designs to generate.",
        ))
        st.caption("Range: [1, 100000]")
    with st.expander("View Optional Settings", expanded=False):
        settings["diffusion_steps"] = int(st.number_input(
            "Diffusion Steps", min_value=1, max_value=200, value=50, step=1,
            key="rfpeptides-diffusion-steps", help="Number of RFdiffusion denoising steps (T). The default 50 matches the RFpeptides paper.",
        ))

    st.info("Running with Rosetta disabled. This tool returns cyclic peptide backbones; sequence design is a separate step.")
    try:
        inputs.validate_request(settings, require_files=bool(target_file))
        if not target_file:
            error = "Load a target structure before submitting."
    except (TypeError, ValueError, OSError) as exc:
        error = str(exc)
    if error:
        st.error(error)

    st.button("Load example inputs", icon=":material/science:", key="rfpeptides-example", on_click=_load_example_inputs, args=(owner_id,))
    can_submit = bool(project and target_file and target_chains and not error)
    if st.button("Submit", type="primary", use_container_width=True, disabled=not can_submit, key="rfpeptides-submit"):
        try:
            inputs.validate_request(settings, require_files=True)
            submit(schema, "design.rfpeptides", settings, project, {"prioritize": False}, "RFpeptides")
        except (TypeError, ValueError, OSError) as exc:
            st.error(str(exc))
    if not project:
        st.caption("Select a project to submit.")


def _result_path(value: object, run_dir: Path) -> Path | None:
    if not value:
        return None
    path = Path(str(value))
    if not path.is_absolute():
        path = run_dir / path
    return path if path.is_file() else None


def render_results(run_id: str, run_dir: Path, candidates: list[dict], parameters: dict, result: dict) -> None:
    st.subheader("RFpeptides cyclic peptide backbones")
    metrics = result.get("metrics") if isinstance(result.get("metrics"), dict) else {}
    st.caption(
        f"{metrics.get('candidate_count', len(candidates))} backbone(s) · "
        f"length {parameters.get('binder_length', '—')} · {inputs.CODE_VERSION}"
    )
    st.info("RFpeptides generated cyclic peptide backbone structures only. The displayed structures do not contain designed amino-acid sequences.")
    if not candidates:
        st.warning("No cyclic peptide structures were collected for this run.")
        return
    rows = []
    for rank, candidate in enumerate(candidates, 1):
        rows.append({
            "Rank": rank,
            "Design": candidate.get("candidate_id"),
            "Peptide chain": ",".join(candidate.get("binder_chains") or ["A"]),
            "Length": candidate.get("binder_length"),
            "Topology": "Head-to-tail cyclic",
            "Sequence design": "Not run",
        })
    frame = pd.DataFrame(rows)
    st.dataframe(frame, hide_index=True, use_container_width=True, height=min(410, 38 * (len(frame) + 1) + 3))
    selected = st.selectbox(
        "Design to visualize", range(len(candidates)),
        format_func=lambda index: f"{index + 1}. {candidates[index].get('candidate_id')}",
        key=f"rfpeptides-result-{run_id}",
    )
    candidate = candidates[selected]
    complex_path = _result_path(candidate.get("complex_pdb"), run_dir)
    peptide_path = _result_path(candidate.get("binder_pdb"), run_dir)
    if complex_path:
        molstar_custom_component(
            [StructureVisualization(pdb=complex_path.read_text(errors="ignore"), color="chain-id", representation_type="cartoon")],
            key=f"rfpeptides-result-view-{run_id}-{selected}", height=560,
            show_controls=True, show_structure_tools=False,
        )
        download_cols = st.columns(2)
        download_cols[0].download_button("Download target–peptide complex", complex_path.read_bytes(), file_name=complex_path.name, key=f"rfpeptides-complex-download-{run_id}-{selected}")
        if peptide_path:
            download_cols[1].download_button("Download peptide backbone", peptide_path.read_bytes(), file_name=peptide_path.name, key=f"rfpeptides-peptide-download-{run_id}-{selected}")
    with st.expander("Job parameters and provenance", expanded=False):
        st.json({
            "engine_version": metrics.get("engine_version"),
            "diffusion_steps": metrics.get("diffusion_steps"),
            "target_chain_map": metrics.get("target_chain_map"),
            "selected_hotspots": metrics.get("selected_hotspots"),
            "sequence_design": metrics.get("sequence_design"),
        })
