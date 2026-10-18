"""Tamarind-inspired input page for Biohub ESMFold2 binder design."""
from __future__ import annotations

import re

import streamlit as st

from mn_molecular_design_studio.esmfold2_binder_frameworks import (
    get_vhh_framework,
    vhh_framework_choices,
)


PD_L1_EXAMPLE = (
    "AFTVTVPKDLYVVEYGSNMTIECKFPVEKQLDLAALIVYWEMEDKNIIQFVHGEEDLKVQHSSYRQRARLLKDQLSLGNAALQITDVKLQDAGVYRCMISYGGADYKRITVKVNA"
)
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
FRAMEWORKS = {
    "Trastuzumab framework": "trastuzumab_framework_vhvl",
    "Atezolizumab framework": "atezolizumab_framework_vhvl",
    "Ocankitug framework": "ocankitug_framework_vhvl",
}
VHH_FRAMEWORKS = vhh_framework_choices()
VHH_FRAMEWORK_OPTIONS = {
    str(record["display_name"]): str(record["id"])
    for record in VHH_FRAMEWORKS
}


def _sequence_from_text(text: str) -> str:
    sequence = "".join(
        line.strip() for line in text.splitlines() if line.strip() and not line.lstrip().startswith(">")
    ).upper()
    sequence = re.sub(r"\s+", "", sequence)
    if not sequence:
        raise ValueError("Enter a target protein sequence.")
    invalid = sorted(set(sequence) - STANDARD_AA)
    if invalid:
        raise ValueError(f"Target sequence contains unsupported characters: {', '.join(invalid)}")
    return sequence


def _template_from_text(text: str) -> str:
    template = re.sub(r"\s+", "", text).upper()
    if not template:
        raise ValueError("Paste a VHH framework template and mark its three CDRs with #.")
    if len(re.findall(r"#+", template)) != 3:
        raise ValueError("Mark exactly three CDR runs with #. Keep framework residues as amino-acid letters.")
    invalid = sorted(set(template) - STANDARD_AA - {"#"})
    if invalid:
        raise ValueError(f"VHH template contains unsupported characters: {', '.join(invalid)}")
    return template


def render_setup(schema: dict, *, project_selector, submit) -> None:
    st.html("""<style>
    .block-container {max-width:1240px!important;padding-top:1.5rem!important}
    .st-key-esm2bd-surface {--esm-green:#27674e;--esm-ink:#173d31;--esm-mint:#edf6f0}
    .st-key-esm2bd-surface .esm2bd-heading {text-align:center;padding:.35rem 0 1.35rem}
    .st-key-esm2bd-surface .esm2bd-eyebrow {font-size:.69rem;letter-spacing:.13em;font-weight:650;color:#658677;margin-bottom:.7rem}
    .st-key-esm2bd-surface .esm2bd-heading h1 {font-size:2rem;color:#173d31;margin:0;font-weight:700;letter-spacing:-.035em}
    .st-key-esm2bd-surface .esm2bd-heading p {font-size:.92rem;color:#6c8276;margin:.65rem auto .4rem;max-width:620px;line-height:1.6}
    .st-key-esm2bd-surface [data-testid="stVerticalBlockBorderWrapper"] {border-color:#d7e6dc!important;border-radius:12px!important}
    .st-key-esm2bd-input-card {padding:1.15rem 1.3rem!important;background:#fff;border:1px solid #d7e6dc;border-radius:12px;margin:.5rem 0 .8rem}
    .st-key-esm2bd-sampling-card {padding:1rem 1.3rem!important;background:#f3f8f5;border:1px solid #d7e6dc;border-radius:12px;margin:.3rem 0 .6rem}
    .st-key-esm2bd-surface .esm2bd-section {font-size:.72rem;font-weight:650;letter-spacing:.1em;color:#5c7b68;margin:0 0 .8rem}
    .st-key-esm2bd-surface [data-variant="pills"] {border-radius:7px!important;border:1px solid #c2d9cc!important;background:#f7faf8;color:#38614e!important;min-height:2.55rem;padding:.42rem .85rem!important;font-size:.84rem!important;box-shadow:none!important}
    .st-key-esm2bd-surface [data-variant="pills"][aria-checked="true"] {background:#27674e!important;color:white!important;border-color:#27674e!important}
    .st-key-esm2bd-surface [data-testid="stTextInput"] input,.st-key-esm2bd-surface [data-testid="stNumberInput"] input,.st-key-esm2bd-surface textarea {background:#f5f9f6!important;color:#193e30!important}
    .st-key-esm2bd-surface [data-baseweb="input"],.st-key-esm2bd-surface [data-baseweb="textarea"] {border-color:#dce7df!important;border-radius:7px!important;background:#f5f9f6!important}
    .st-key-esm2bd-surface [data-testid="stButton"] button {border-radius:7px!important;color:#27674e;border-color:#bad2c3}
    .st-key-esm2bd-surface [data-testid="stButton"] button:hover {background:#edf6f0!important;border-color:#659d7d!important;color:#1d5039!important}
    .st-key-esm2bd-surface button[kind="primary"] {background:#27674e!important;border-color:#27674e!important;color:#fff!important;min-height:2.8rem}
    .st-key-esm2bd-surface button[kind="primary"]:disabled {background:#99b8a5!important;border-color:#99b8a5!important;color:#fff!important}
    .st-key-esm2bd-surface [data-testid="stExpander"] {border-color:#d7e6dc!important;border-radius:8px!important;background:#fcfefc}
    .st-key-esm2bd-surface [data-testid="stCaptionContainer"] {color:#7c8f82}
    </style>""")
    with st.container(key="esm2bd-surface"):
        _render_setup_content(schema, project_selector=project_selector, submit=submit)


def _render_setup_content(schema: dict, *, project_selector, submit) -> None:
    if st.button("← Back to All tools", key="esm2bd-back"):
        st.session_state["page"] = "All tools"
        st.session_state.pop("selected_workflow", None)
        st.rerun()

    st.html(
        '<div class="esm2bd-heading"><div class="esm2bd-eyebrow">DE NOVO BINDER DESIGN</div>'
        '<h1>ESMFold2 Binder Design</h1>'
        '<p>Design protein binders against a target sequence with gradient-guided ESMFold2 optimization.</p></div>'
    )
    with st.expander("About this workflow", expanded=False):
        st.write("The Biohub method optimizes binder sequences using ESMFold2 structure and interface predictions with ESMC sequence regularization. It runs from sequence alone; no target PDB or MSA is required.")
        st.markdown("[Tamarind workflow reference](https://app.tamarind.bio/tools/esmfold2-binder-design) · [Biohub ESM source](https://github.com/Biohub/esm/tree/main/cookbook/tutorials)")
        st.caption("This tile uses the installed Biohub ESMFold2 experimental fast checkpoint. The upstream protocol provides miniprotein and scFv prompts; the VHH presets come from the local ESMFold2 pipeline framework panel and are passed as masked sequence prompts.")

    job_col, project_col = st.columns([1.1, 1.5])
    with job_col:
        job_name = st.text_input("Job name", "esmfold2-binder-design", key="esm2bd-job-name").strip()
    with project_col:
        project_id = project_selector("esmfold2-binder-design")

    st.markdown("**EXAMPLES**")
    if st.button("PD-L1 target sequence", icon=":material/science:", key="esm2bd-example-pdl1"):
        st.session_state["esm2bd-target-sequence"] = PD_L1_EXAMPLE
        st.session_state["esm2bd-hotspots"] = ""
        st.rerun()

    modality_labels = {"Miniprotein": "miniprotein", "Antibody (scFv)": "scfv", "Nanobody (VHH)": "vhh"}
    modality_label = st.pills(
        "Binder type",
        list(modality_labels),
        selection_mode="single",
        default="Miniprotein",
        key="esm2bd-binder-type",
        label_visibility="collapsed",
    ) or "Miniprotein"
    binder_type = modality_labels[modality_label]
    st.caption({
        "miniprotein": "Design a de novo miniprotein binder with a fixed sequence length.",
        "scfv": "Design the six CDR loops on a single-chain VH–VL antibody framework.",
        "vhh": "Design the three CDR loops on a single-domain VHH framework. Framework residues stay fixed.",
    }[binder_type])

    with st.container(key="esm2bd-input-card"):
        st.html('<div class="esm2bd-section">TARGET SEQUENCE</div>')
        st.session_state.setdefault("esm2bd-target-sequence", "")
        target_text = st.text_area(
            "Target protein sequence",
            key="esm2bd-target-sequence",
            height=150,
            placeholder="Paste a protein sequence or a single-record FASTA here…",
            help="Only the 20 standard amino acids are supported. FASTA headers and line breaks are ignored.",
        )
        target_sequence = ""
        input_error = ""
        try:
            target_sequence = _sequence_from_text(target_text)
            st.caption(f"{len(target_sequence)} target residues")
        except ValueError as exc:
            input_error = str(exc)
        hotspots = st.text_input(
            "Target hotspot residues (optional)",
            key="esm2bd-hotspots",
            placeholder="For example, 28, 45",
            help="Enter 1-based positions in the sequence. The residue letter is inferred and checked automatically.",
        )

        binder_length = 80
        binder_framework = ""
        binder_template = ""
        if binder_type == "miniprotein":
            binder_length = int(st.number_input("Binder length", min_value=40, max_value=200, value=80, step=1, key="esm2bd-binder-length"))
        elif binder_type == "scfv":
            framework_label = st.selectbox("Antibody framework", list(FRAMEWORKS), key="esm2bd-scfv-framework")
            binder_framework = FRAMEWORKS[framework_label]
        else:
            framework_options = {**VHH_FRAMEWORK_OPTIONS, "Custom masked template": "custom"}
            framework_label = st.selectbox(
                "VHH framework",
                list(framework_options),
                key="esm2bd-vhh-framework",
            )
            binder_framework = framework_options[framework_label]
            if binder_framework == "custom":
                binder_template = st.text_area(
                    "VHH framework template",
                    key="esm2bd-vhh-template",
                    height=120,
                    placeholder="Paste a single-domain VHH sequence and replace its three CDRs with #…",
                    help="Mark exactly the three CDR runs to redesign with # characters. Keep framework residues as amino-acid letters.",
                )
                try:
                    binder_template = _template_from_text(binder_template)
                except ValueError as exc:
                    input_error = input_error or str(exc)
            else:
                record = get_vhh_framework(binder_framework)
                cdr_ranges = record["cdr_lengths"]
                cdr_summary = ", ".join(
                    f"{name.upper()}: {bounds['min']}–{bounds['max']} aa"
                    for name, bounds in cdr_ranges.items()
                )
                st.caption(
                    f"{record['source_target']} framework · variable CDR lengths: {cdr_summary}. "
                    "The CDRs are redesigned; this is not transfer of the source antibody's binding sequence."
                )

    with st.container(key="esm2bd-sampling-card"):
        st.html('<div class="esm2bd-section">DESIGN SAMPLING</div>')
        num_designs = int(st.number_input("Number of designs", min_value=1, max_value=128, value=1, step=1, key="esm2bd-num-designs"))

    with st.expander("Optional settings"):
        a, b = st.columns(2)
        with a:
            optimization_steps = int(st.number_input("Optimization steps", min_value=10, max_value=300, value=150, step=10, key="esm2bd-optimization-steps"))
        with b:
            seed = int(st.number_input("Random seed", min_value=0, max_value=2_147_483_647, value=0, step=1, key="esm2bd-seed"))
    prioritize = st.checkbox("Prioritize in scheduler queue", key="esm2bd-prioritize")

    if input_error:
        st.error(input_error)
    can_submit = bool(project_id and target_sequence and not input_error)
    if st.button("Submit", type="primary", use_container_width=True, disabled=not can_submit, key="esm2bd-submit"):
        parameters = {
            "target_sequence": target_sequence,
            "binder_type": binder_type,
            "binder_framework": binder_framework,
            "binder_template": binder_template,
            "hotspots": hotspots,
            "binder_length": binder_length,
            "num_designs": num_designs,
            "campaign_name": job_name or "esmfold2-binder-design",
            "optimization_steps": optimization_steps,
            "seed": seed,
        }
        submit(
            schema,
            "design.esmfold2_binder",
            parameters,
            project_id,
            {"prioritize": prioritize},
            "ESMFold2 Binder Design",
        )
    elif not project_id:
        st.caption("Select a project to submit.")
