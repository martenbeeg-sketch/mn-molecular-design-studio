"""Review a saved PDB and submit a focused MODELLER preparation job."""
from __future__ import annotations

import json
import secrets
from pathlib import Path

import streamlit as st

from mn_molecular_design_studio.database_actions_ui import ERRORS
from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import workflow_schema
from mn_molecular_design_studio.engine.workflow_engine.workflows.ligand_tools import inspect_complex_gaps
from mn_molecular_design_studio.protein_modelling import (
    MAX_LOOP_MODELS, MAX_MODELED_GAP, loop_residues, observed_chains, parse_options,
)


@st.cache_data(max_entries=16, show_spinner=False)
def _gaps(source: str):
    return inspect_complex_gaps(Path(source))


def _gap_options(path, key):
    gaps = _gaps(str(path))
    by_gap = {f"Chain {gap['chain']} · {gap['left_code']}{gap['left_resseq']}–{gap['right_code']}{gap['right_resseq']}": gap
              for gap in gaps if gap["numbering_gap_length"] <= MAX_MODELED_GAP}
    if not by_gap:
        st.info(f"No internal numbering gaps of up to {MAX_MODELED_GAP} residues with observed flanks were found.")
        return {}, False
    selected = st.multiselect("Missing segments to rebuild", list(by_gap), key=key + "-gaps",
                              help="Compare each numbering gap with the expected sequence before confirming.")
    definitions = []
    valid = bool(selected)
    for label in selected:
        gap = by_gap[label]
        gap_key = f"{key}-gap-{gap['chain']}-{gap['left_resseq']}-{gap['right_resseq']}"
        st.markdown(f"**{label}**")
        count = gap["numbering_gap_length"]
        sequence = st.text_input("Missing amino-acid sequence", value=gap["suggestion"], max_chars=MAX_MODELED_GAP,
                                 key=gap_key + "-sequence", help=f"Enter {count} one-letter codes in sequence order.")
        confirmed = st.checkbox("I confirm the missing sequence and both flanking residues", key=gap_key + "-confirmed")
        valid = valid and confirmed
        definitions.append({field: gap[field] for field in ("chain", "left_resseq", "right_resseq")})
        definitions[-1].update(sequence=sequence, confirmed=confirmed)
    return {"gaps": definitions}, valid


def _loop_options(source, key):
    chains = {chain: residues for chain, residues in observed_chains(source).items() if len(residues) >= 3}
    if not chains:
        st.info("An existing loop needs a protein chain with observed residues on both sides.")
        return {}, False
    chain = st.selectbox("Chain", list(chains), key=key + "-loop-chain")
    labels = list(chains[chain])[1:-1]
    chain_key = f"{key}-loop-{chain}"
    start = st.selectbox("First residue", labels, index=None, key=chain_key + "-start")
    end = st.selectbox("Last residue", labels, index=None, key=chain_key + "-end")
    models = st.number_input("MODELLER conformations", 1, MAX_LOOP_MODELS, 5, key=chain_key + "-models")
    selection = {"chain": chain, "start": start, "end": end, "models": models}
    if start is None or end is None:
        return {}, False
    chosen = loop_residues(chains, selection)
    st.caption("Selected sequence · " + " ".join(chains[chain][label].resname for label in chosen))
    return {"loop": selection}, True


def render_modelling(path, source, key, projects, preferred_projects, submit):
    panel = st.expander("Model missing residues or an existing loop", key=key + "-modelling", on_change="rerun")
    if not panel.open:
        return
    with panel:
        operation = st.radio("MODELLER operation", ["Rebuild missing residues", "Remodel an existing loop"],
                             key=key + "-operation")
        options, ready = {}, False
        try:
            if operation == "Rebuild missing residues":
                options, ready = _gap_options(path, key)
            else:
                options, ready = _loop_options(source, key)
            if ready:
                parse_options(json.dumps(options))
        except (*ERRORS, RuntimeError) as exc:
            st.error(str(exc))
            ready = False
        active = {project["id"]: project["name"] for project in projects
                  if not project.get("archived_at") and not project.get("deleted_at")}
        ids = list(active)
        preferred = next((project for project in preferred_projects if project in active), None)
        project = st.selectbox("Project for the modelling job", ids, index=ids.index(preferred) if preferred else None,
                               format_func=active.get, key=key + "-modelling-project")
        if not ids:
            st.caption("Create or select an active project to run MODELLER.")
        if st.button("Run MODELLER", help="Uses the MODELLER installation on this workstation.", disabled=not ready or not project, type="primary", key=key + "-run-modeller"):
            try:
                workflow = "ligand.complex.prepare"
                execution = {"placement": "local", "idempotency_key": secrets.token_urlsafe(24)}
                result = submit(workflow_schema(workflow), workflow,
                                {"source_pdb": str(path), "repair": False, "modelling_json": json.dumps(options)},
                                project, {"_execution": execution}, operation)
                if result and (run_id := result["job"].get("run_id")):
                    task_group = result["job"].get("task_group") or "ligand-preparation"
                    st.session_state[key + "-last-modelling-job"] = f"{task_group}/{run_id}"
            except (*ERRORS, RuntimeError) as exc:
                st.error(f"Could not submit the MODELLER job: {exc}")
        if reference := st.session_state.get(key + "-last-modelling-job"):
            if st.button("View modelling job", key=key + "-view-modeller"):
                st.session_state["selected-job-reference"] = reference
                st.session_state["page"] = "Job details"
                st.rerun()
