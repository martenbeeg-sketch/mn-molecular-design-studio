"""Tamarind-style input page for the RFDiffusion4-Proteina generator."""
from __future__ import annotations

import shutil
import uuid
from pathlib import Path

import streamlit as st

from mn_molecular_design_studio import rfd4_proteina_input as inputs
from mn_molecular_design_studio.auth import user_input_root
from mn_molecular_design_studio.target_editor import render_target_editor, render_target_source
from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer import (
    StructureVisualization,
    molstar_custom_component,
)


EXAMPLE_SPEC = """# Protein binder condition spec. Upload target.cif below.
pd_l1_binder:
  input: target.cif
  contig: \"60-120, /0, A17-131\"
  conditions:
    C_CRDGRP:
      1:
        - select: \"sel('A/*/*/[N,CA,C,O]')\"
    C_SEQ:
      true:
        - select: \"sel('A')\"
    C_HOT:
      true:
        - select: \"sel('A/*/[56,115,123]')\"
"""


def _load_example_inputs(example: str, owner_id: str) -> None:
    if example == "pd_l1_binder":
        task, input_mode = "protein_binder", "guided_binder"
        example_target = Path(__file__).resolve().parent / "assets" / "proteina_complexa" / "PD-L1.pdb"
        stored_target = user_input_root(owner_id) / "rfd4-proteina" / "PD-L1.pdb"
        stored_target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(example_target, stored_target)
        st.session_state["rfd4-target-fetched-target"] = str(stored_target)
        st.session_state["rfd4-target-fetched-target-name"] = "PD-L1"
        st.session_state.pop("rfd4-target-target-upload", None)
        st.session_state["rfd4-target-pdb-code"] = "5O45"
        st.session_state["rfd4-job-name"] = f"rfd4-proteina-example-{uuid.uuid4().hex[:5]}"
        st.session_state["rfd4-length-min"] = 60
        st.session_state["rfd4-length-max"] = 120
        notice = "PD-L1 binder example loaded with its target structure and chain A selected. Set a project before submitting."
    else:
        task, input_mode = "unconditional", "unconditional"
        st.session_state["rfd4-job-name"] = "rfdiffusion4-unconditional-example"
        notice = "Unconditional monomer example loaded. Choose a project before submitting."

    st.session_state["rfd4-task-pills"] = task
    st.session_state[f"rfd4-input-mode-pills-{task}"] = input_mode
    st.session_state["rfd4-count"] = 8 if example == "pd_l1_binder" else 1
    st.session_state["rfd4-example-loaded-notice"] = notice


def render_setup(schema: dict, *, owner_id: str, project_selector, submit) -> None:
    st.html("""<style>
    .block-container {max-width:800px!important;padding:1rem 1.15rem 2rem!important;margin:0 auto!important}
    [data-testid="stMainBlockContainer"] {background:#fff;border:1px solid #d9ddda;border-radius:11px;box-shadow:0 2px 8px rgba(16,24,40,.08);padding:1.1rem 1.35rem 1.5rem}
    .st-key-rfd4-surface {--rfd4-green:#28674d;--rfd4-ink:#202426;--rfd4-mint:#edf5f0}
    .st-key-rfd4-surface .rfd4-heading {text-align:center;padding:.1rem 0 .55rem}
    .st-key-rfd4-surface .rfd4-heading h1 {font-size:1.9rem;color:#17191b;margin:0;font-weight:650;letter-spacing:-.025em}
    .st-key-rfd4-surface .rfd4-heading .badge {display:inline-block;vertical-align:middle;margin-left:.35rem;padding:.16rem .48rem;border-radius:1rem;background:#28674d;color:#fff;font-size:.62rem;font-weight:700;letter-spacing:.02em}
    .st-key-rfd4-surface .rfd4-heading .hero-link {display:inline-block;vertical-align:middle;margin-left:.4rem;color:#68756e;text-decoration:none;font-size:.78rem;font-weight:650}
    .st-key-rfd4-surface .rfd4-heading .hero-link:hover {color:#28674d}
    .st-key-rfd4-surface .rfd4-heading p {font-size:.92rem;color:#52616a;margin:.55rem auto .45rem;max-width:700px;line-height:1.5}
    .st-key-rfd4-surface .rfd4-meta {font-size:.8rem;color:#69767f}
    .st-key-rfd4-surface .rfd4-project-reminder {font-size:.78rem;color:#7b858b;margin:-.15rem 0 .65rem}
    .st-key-rfd4-surface [data-testid="stVerticalBlockBorderWrapper"] {border-color:#d9ddda!important;border-radius:9px!important}
    .st-key-rfd4-surface [data-testid="stButton"] button {border-radius:7px!important}
    .st-key-rfd4-surface button[kind="primary"] {background:#28674d!important;border-color:#28674d!important;color:#fff!important;min-height:2.6rem}
    .st-key-rfd4-surface button[kind="primary"]:disabled {background:#91b5a4!important;border-color:#91b5a4!important;color:#fff!important}
    .st-key-rfd4-surface [data-testid="stFileUploaderDropzone"] {background:#f3f8f5!important;border:1px dashed #b9d2c2!important;border-radius:8px}
    .st-key-rfd4-surface [data-testid="stExpander"] {border-color:#d9ddda!important;border-radius:8px!important;background:#fff}
    .st-key-rfd4-surface [data-testid="stTextArea"] textarea {font-family:ui-monospace,SFMono-Regular,Consolas,monospace;font-size:.83rem}
    .st-key-rfd4-surface [data-testid="stTabs"] button[aria-selected="true"] {color:#28674d!important;border-bottom-color:#28674d!important}
    .st-key-rfd4-surface .st-key-rfd4-task-pills button,
    .st-key-rfd4-surface [class*="st-key-rfd4-input-mode-pills"] button {border-radius:7px!important;border:1px solid #d5d9d7!important}
    .st-key-rfd4-surface .st-key-rfd4-task-pills button[aria-pressed="true"],
    .st-key-rfd4-surface [class*="st-key-rfd4-input-mode-pills"] button[aria-pressed="true"] {background:#e8ece9!important;color:#20332a!important;border-color:#cfd8d2!important;box-shadow:none!important}
    </style>""")
    with st.container(key="rfd4-surface"):
        _render_content(schema, owner_id=owner_id, project_selector=project_selector, submit=submit)


def _render_content(schema: dict, *, owner_id: str, project_selector, submit) -> None:
    if st.button("← Back to All tools", key="rfd4-back"):
        st.session_state["page"] = "All tools"
        for key in ("selected_workflow", "selected_workflow_engine", "selected_workflow_group"):
            st.session_state.pop(key, None)
        st.rerun()
    st.html(
        '<div class="rfd4-heading"><h1>RFDiffusion4-Proteina <span class="badge">NEW</span>'
        '<a class="hero-link" href="https://github.com/RosettaCommons/RFD4-Proteina" target="_blank" title="Source code" aria-label="Source code">GH</a>'
        '<a class="hero-link" href="https://rfd4-proteina.org/latest/" target="_blank" title="Paper and documentation" aria-label="Paper and documentation">▤</a></h1>'
        '<p>Generates protein structure and sequence together, from scratch or around a target, motif or ligand.</p>'
        f'<span class="rfd4-meta">Protein Design · ~10 min · {inputs.UPSTREAM_VERSION}</span></div>'
    )
    action_cols = st.columns([1, 1, 1, .8])
    with action_cols[0]:
        with st.popover("▤ Show Outputs"):
            st.write("Generated CIF structures, normalized sequences, and candidate records are saved in Results.")
    with action_cols[1]:
        with st.popover("⚙ Show Settings"):
            st.write("Choose a design task, input style, design count, GPU batch size, and sampling effort. Seed and CPU controls are also available.")
    with action_cols[2]:
        with st.popover("↗ API Docs"):
            st.markdown("[RFDiffusion4-Proteina GitHub](https://github.com/RosettaCommons/RFD4-Proteina) · [Documentation](https://rfd4-proteina.org/latest/) · [Paper](https://rfd4-proteina.org/assets/paper/RFD4-Proteina.pdf)")
    with action_cols[3]:
        with st.popover("⌄ Details"):
            st.write("This tile runs the upstream generator in a CUDA 13 container. Inline refolding and oracle scoring are not enabled; use the downstream prediction and analysis tools for those steps.")
            st.caption("The same image is intended for RTX 4090 and RTX 5090 hosts. Upstream documents a CUDA 13-compatible driver, but does not certify these exact cards in its hardware list.")
    st.divider()

    design_tab, examples_tab = st.tabs(["Design", "Examples"])
    with design_tab:
        if notice := st.session_state.pop("rfd4-example-loaded-notice", None):
            st.success(notice)
        name_col, project_col = st.columns([1.2, 1.2])
        with name_col:
            job_name = st.text_input("Job name", "rfd4-proteina-design", key="rfd4-job-name").strip() or "rfd4-proteina-design"
        with project_col:
            project_id = project_selector("rfd4-proteina")
        st.markdown(
            '<div class="rfd4-project-reminder">Recommended: keep related design, prediction, and analysis jobs in the same project so their inputs and results stay together.</div>',
            unsafe_allow_html=True,
        )

        task = st.pills(
            "Design task",
            list(inputs.TASKS),
            selection_mode="single",
            default="protein_binder",
            format_func=inputs.TASKS.get,
            key="rfd4-task-pills",
            help="Choose the design problem. The available inputs below update to match the selected task.",
            width="stretch",
        ) or "protein_binder"
        st.caption(inputs.TASK_HELP[task])

        input_modes = ["guided_binder", "annotated_cif", "condition_spec"]
        if task == "unconditional":
            input_modes = ["unconditional", "condition_spec", "annotated_cif"]
        selected_default = "guided_binder" if task == "protein_binder" else ("unconditional" if task == "unconditional" else "annotated_cif")
        input_mode = st.pills(
            "Input method",
            input_modes,
            selection_mode="single",
            default=selected_default,
            format_func={
                "guided_binder": "Target structure and hotspot selection",
                "annotated_cif": "RFD4 annotated CIF/mmCIF",
                "condition_spec": "Native YAML/JSON conditions",
                "unconditional": "Bundled monomer template",
            }.get,
            key="rfd4-input-mode-pills-" + task,
            width="stretch",
        ) or selected_default

        parameters = {
            "task": task,
            "input_mode": input_mode,
            "job_name": job_name,
            "input_files": [],
            "target_chains": [],
            "hotspots": "",
            "condition_spec": "",
        }
        ready = True
        if input_mode == "guided_binder":
            source = render_target_source(
                key_prefix="rfd4-target",
                label="Target",
                example_code="5O45",
                owner_id=owner_id,
                allow_cif=True,
                upload_help="Protein target as PDB or mmCIF. The selected chains are converted to an RFD4-compatible mmCIF input.",
            )
            if source.path:
                target = render_target_editor(
                    source,
                    key_prefix="rfd4-target-editor",
                    label="Target",
                    residue_selection_title="Hotspot residues (optional)",
                    show_structure_tools=False,
                    allow_trimming=True,
                    viewer_height=430,
                )
                parameters.update(
                    target_structure=str(target.design_path or source.path),
                    target_chains=target.chains,
                    hotspots=target.hotspots,
                )
                ready = target.ready
                min_col, max_col = st.columns(2)
                with min_col:
                    parameters["length_min"] = int(st.number_input(
                        "Designed chain minimum length", min_value=1, max_value=3000, value=60,
                        key="rfd4-length-min",
                    ))
                with max_col:
                    parameters["length_max"] = int(st.number_input(
                        "Designed chain maximum length", min_value=1, max_value=3000, value=120,
                        key="rfd4-length-max",
                    ))
                if parameters["length_min"] > parameters["length_max"]:
                    st.error("Minimum chain length must not exceed maximum chain length.")
                    ready = False
                with st.expander("How the guided input is translated", expanded=False):
                    st.write("Studio writes an RFD4 condition spec that preserves the selected target chains, conditions their sequence and backbone geometry, applies the optional hotspot residues, and inserts a new binder chain in the selected length range.")
            else:
                ready = False

        elif input_mode == "annotated_cif":
            uploaded = st.file_uploader(
                "Upload annotated CIF/mmCIF structures",
                type=["cif", "mmcif"],
                accept_multiple_files=True,
                key="rfd4-annotated-files-" + task,
                help="Upload RFDiffusion4-Proteina annotated inputs. For multiple structures, each file becomes one named design input.",
            )
            try:
                parameters["input_files"] = [str(path) for path in inputs.store_structure_uploads(uploaded or [], owner_id)]
            except (OSError, ValueError) as exc:
                st.error(str(exc))
                ready = False
            ready = ready and bool(parameters["input_files"])
            if task == "unconditional" and not uploaded:
                ready = True
                parameters["input_mode"] = "unconditional"

        elif input_mode == "condition_spec":
            uploaded = st.file_uploader(
                "Upload condition-spec CIF/mmCIF files",
                type=["cif", "mmcif"],
                accept_multiple_files=True,
                key="rfd4-custom-files-" + task,
                help="Use the uploaded filenames in each design's input field, such as input: target.cif.",
            )
            try:
                parameters["input_files"] = [str(path) for path in inputs.store_structure_uploads(uploaded or [], owner_id)]
            except (OSError, ValueError) as exc:
                st.error(str(exc))
                ready = False
            parameters["condition_spec"] = st.text_area(
                "Condition specification",
                value=st.session_state.get("rfd4-condition-spec", EXAMPLE_SPEC),
                height=260,
                key="rfd4-condition-spec",
                help="Native RFDiffusion4-Proteina condition-spec YAML or JSON. Inputs must refer to CIF files uploaded above.",
            )
            ready = ready and bool(parameters["input_files"]) and bool(parameters["condition_spec"].strip())

        elif input_mode == "unconditional":
            st.info("Uses the bundled RFDiffusion4-Proteina short monomer template. No structure upload is needed.")

        with st.expander("Sampling and compute", expanded=True):
            count_col, batch_col, effort_col = st.columns(3)
            with count_col:
                parameters["num_designs"] = int(st.number_input("Number of designs", 1, 10000, 8, key="rfd4-count"))
            with batch_col:
                parameters["batch_size"] = int(st.number_input("GPU batch size", 1, 16, 1, key="rfd4-batch", help="Start with 1 on 24 GB GPUs; increase only if memory allows."))
            with effort_col:
                effort = st.selectbox("Sampling effort", ["Fast", "Balanced", "Thorough"], index=1, key="rfd4-effort")
            effort_steps = {"Fast": 30, "Balanced": 60, "Thorough": 90}
            advanced = st.checkbox("Set sampling steps manually", key="rfd4-manual-steps")
            parameters["nsteps"] = int(st.number_input("Sampling steps", 1, 200, 60 if advanced else effort_steps[effort], disabled=not advanced, key="rfd4-steps"))
            adv_a, adv_b = st.columns(2)
            with adv_a:
                parameters["seed"] = int(st.number_input("Random seed (0 = choose automatically)", 0, 2147483647, 0, key="rfd4-seed"))
            with adv_b:
                parameters["cpu_cores"] = int(st.number_input("CPU cores", 1, 256, 8, key="rfd4-cpus"))
            parameters["compile_model"] = st.checkbox("Use torch.compile (longer first-run warm-up)", value=False, key="rfd4-compile")
            st.caption("The model weights are downloaded on the first run and cached on the workstation. Batch size controls how many designs are sampled at once; the requested designs are generated in batches on the assigned GPU.")

        prioritize = st.checkbox("Prioritize in scheduler queue", key="rfd4-priority")
        st.button(
            "Load example inputs",
            key="rfd4-load-example-inputs",
            on_click=_load_example_inputs,
            args=("pd_l1_binder", owner_id),
            help="Loads the PD-L1 example structure, selects chain A, and sets eight designs, as on Tamarind.",
        )
        submit_cols = st.columns([1, 1, 1])
        with submit_cols[1]:
            submit_clicked = st.button("Submit", type="primary", disabled=not project_id or not ready, key="rfd4-submit")
        if submit_clicked:
            submit(schema, "design.rfd4_proteina", parameters, project_id, {"prioritize": prioritize}, "RFDiffusion4-Proteina")
        if not project_id:
            st.markdown('<div style="text-align:center;color:#748087;font-size:.82rem;font-style:italic">Select a project to submit.</div>', unsafe_allow_html=True)

    with examples_tab:
        st.subheader("Examples")
        st.write("Choose one of Tamarind's native task patterns, then enter it as an annotated CIF or a condition specification on the Design tab.")
        example_cols = st.columns(2)
        with example_cols[0]:
            st.button(
                "Load PD-L1 binder example",
                key="rfd4-load-pd-l1-example",
                on_click=_load_example_inputs,
                args=("pd_l1_binder", owner_id),
                use_container_width=True,
            )
        with example_cols[1]:
            st.button(
                "Load unconditional monomer example",
                key="rfd4-load-unconditional-example",
                on_click=_load_example_inputs,
                args=("unconditional", owner_id),
                use_container_width=True,
            )
        st.markdown(
            "- **Protein binder:** guided target-chain and hotspot form; the app writes a native RFD4 condition spec.\n"
            "- **Antibody / nanobody:** use the upstream antibody examples and upload the annotated CIF or native specification.\n"
            "- **Enzyme / motif / ligand pocket:** use an annotated CIF when motif and ligand annotations already exist.\n"
            "- **DNA / RNA binder:** supply the nucleic-acid structure and conditions in an annotated CIF or YAML/JSON spec.\n"
            "- **Existing structure:** use an annotated structure as the starting point for partial diffusion.\n"
            "- **Unconditional generation:** use the bundled short monomer template."
        )
        st.code(EXAMPLE_SPEC, language="yaml")
        st.download_button("Download binder condition-spec example", EXAMPLE_SPEC, "rfd4-proteina-binder-example.yaml", "text/yaml", key="rfd4-example-download")
        st.markdown("[Upstream cookbook](https://rfd4-proteina.org/latest/how-to/cookbook.html) · [Antibody and nanobody examples](https://github.com/RosettaCommons/RFD4-Proteina/tree/main/examples/antibody_design) · [Condition-spec reference](https://rfd4-proteina.org/latest/how-to/input_prep/yaml-json.html)")


def render_results(run_id: str, run_dir: Path, candidates: list[dict], parameters: dict, result: dict) -> None:
    st.subheader("RFDiffusion4-Proteina designs")
    st.caption("RFDiffusion4-Proteina co-designs sequence and structure. Inline refolding is not enabled; continue with sequence analysis or a cofolding validation tool.")
    if not candidates:
        st.info("Generated structures will appear here after RFDiffusion4-Proteina finishes.")
        return
    import pandas as pd

    rows = []
    for candidate in candidates:
        rows.append({
            "Design": candidate.get("candidate_id"),
            "Sequence": candidate.get("binder_sequence") or "—",
            "Chains": ", ".join(candidate.get("binder_chains") or []),
            "Residues": (candidate.get("metrics") or {}).get("generated_residue_count"),
        })
    st.dataframe(pd.DataFrame(rows), hide_index=True, use_container_width=True)
    selected = st.selectbox("Design to view", range(len(candidates)), format_func=lambda i: candidates[i].get("candidate_id", f"Design {i + 1}"), key=f"rfd4-result-{run_id}")
    candidate = candidates[selected]
    path = Path(candidate.get("complex_pdb") or "")
    if path.is_file():
        try:
            molstar_custom_component(
                [StructureVisualization(pdb=path.read_text(errors="ignore"), color="molecule-type", representation_type="cartoon")],
                key=f"rfd4-view-{run_id}-{selected}", height=560, show_controls=True, show_structure_tools=False,
            )
            st.download_button("Download generated CIF", path.read_bytes(), file_name=path.name, key=f"rfd4-download-{run_id}-{selected}")
        except (OSError, UnicodeDecodeError):
            st.warning("The structure is available in the run artifacts but could not be displayed in the viewer.")
