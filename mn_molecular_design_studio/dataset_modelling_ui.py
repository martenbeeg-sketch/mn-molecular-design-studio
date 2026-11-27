"""Run repair and sequence editing on the coordinates shown in Datasets."""
from __future__ import annotations

import json
import secrets
from pathlib import Path

import streamlit as st

from mn_molecular_design_studio.database_actions_ui import ERRORS
from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import workflow_schema
from mn_molecular_design_studio.engine.workflow_engine.workflows.ligand_tools import inspect_complex_gaps
from mn_molecular_design_studio.molecular_database import MolecularDatabase
from mn_molecular_design_studio.protein_modelling import (
    AMINO_ACIDS, MAX_LOOP_MODELS, MAX_MODELED_GAP, loop_residues, observed_chains, parse_options,
)


@st.cache_data(max_entries=16, show_spinner=False)
def _gaps(source: str):
    return inspect_complex_gaps(Path(source))


def _gap_options(path, chains, key):
    gaps = _gaps(str(path))
    by_gap = {f"Chain {gap['chain']} · {gap['left_code']}{gap['left_resseq']}–{gap['right_code']}{gap['right_resseq']}": gap
              for gap in gaps if gap["numbering_gap_length"] <= MAX_MODELED_GAP
              and all(str(gap[field]) in chains.get(gap["chain"], {}) for field in ("left_resseq", "right_resseq"))}
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


def _loop_options(chains, key):
    chains = {chain: residues for chain, residues in chains.items() if len(residues) >= 3}
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


def _extension_options(chains, key):
    chain = st.selectbox("Chain to extend", list(chains), key=key + "-extension-chain")
    st.caption(f"Append after residue {next(reversed(chains[chain]))} on the selected trim.")
    sequence = "".join(st.text_input("Amino acids to append", placeholder="e.g. GSG", max_chars=30,
                                     key=key + "-extension-sequence").upper().split())
    models = st.number_input("MODELLER conformations", min_value=1, max_value=50, value=10,
                             key=key + "-extension-models")
    if set(sequence) - set(AMINO_ACIDS):
        raise ValueError("Enter canonical one-letter amino-acid codes.")
    return {"extension_chain": chain, "extension_sequence": sequence, "extension_models": models}, bool(sequence)


def _substitution_options(chains, key):
    from Bio.SeqUtils import seq1, seq3

    chain = st.selectbox("Chain", list(chains), key=key + "-substitution-chain")
    residues = {label: residue.resname for label, residue in chains[chain].items()
                if all(atom in residue for atom in ("N", "CA", "C"))}
    if not residues:
        st.info("Substitution needs an observed amino acid with complete backbone coordinates.")
        return {}, False
    selected = st.multiselect("Residues to substitute", list(residues), max_selections=50,
                              format_func=lambda label: f"{residues[label]} {label}", key=f"{key}-{chain}-substitutions")
    substitutions = []
    for label in selected:
        target = st.selectbox(f"New amino acid for {residues[label]} {label}",
            [code for code in AMINO_ACIDS if code != seq1(residues[label])], index=None,
            format_func=lambda code: f"{seq3(code).upper()} ({code})", key=f"{key}-{chain}-{label}-substitution")
        if target:
            substitutions.append({"chain": chain, "residue": label, "amino_acid": target})
    return {"substitutions": substitutions}, bool(selected) and len(substitutions) == len(selected)


def _operation_parameters(path, selected_pdb, operation, key):
    if operation == "PDBFixer repair":
        refine = st.checkbox("Refine rebuilt atoms with OpenMM", value=True, key=key + "-openmm",
                             help="Relaxes rebuilt atoms while holding the observed protein backbone fixed.")
        map_modified = st.checkbox("Map supported modified residues", value=True, key=key + "-map-modified")
        return {"repair": True, "refine_rebuilt_positions": refine, "preserve_original_backbone": True,
                "map_modified_residues": map_modified, "max_internal_gap": 0}, True
    chains = {chain: residues for chain, residues in observed_chains(selected_pdb).items() if residues}
    if not chains:
        st.info("Keep a protein chain in the preview to use MODELLER.")
        return {}, False
    if operation == "Extend the C-terminus":
        parameters, ready = _extension_options(chains, key)
        return {"repair": False, **parameters}, ready
    if operation == "Rebuild missing residues or a loop":
        options, ready = _gap_options(path, chains, key)
    elif operation == "Remodel an existing loop":
        options, ready = _loop_options(chains, key)
    else:
        options, ready = _substitution_options(chains, key)
    return {"repair": False, "modelling_options": parse_options(options) if ready else None}, ready


def render_protein_tools(user, source_path, selected_pdb, input_history, key, projects, preferred_projects, submit):
    st.markdown("**Repair and model the selected trim**")
    operation = st.selectbox("Preparation operation", ["PDBFixer repair", "Extend the C-terminus",
        "Substitute amino acids", "Rebuild missing residues or a loop", "Remodel an existing loop"], key=key + "-operation")
    parameters, ready = {}, False
    try:
        parameters, ready = _operation_parameters(source_path, selected_pdb.decode("utf-8"), operation, key)
    except (*ERRORS, RuntimeError) as exc:
        st.error(str(exc))
    active = {project["id"]: project["name"] for project in projects
              if not project.get("archived_at") and not project.get("deleted_at")}
    ids = list(active)
    preferred = next((project for project in preferred_projects if project in active), None)
    project = st.selectbox("Project for the preparation job", ids, index=ids.index(preferred) if preferred else None,
                           format_func=active.get, key=key + "-modelling-project")
    if not ids:
        st.caption("Create or select an active project to run preparation.")
    if st.button("Run preparation", disabled=not ready or not project, type="primary", key=key + "-run-preparation"):
        try:
            source = MolecularDatabase.materialize_structure_input(user, selected_pdb, ".pdb")
            workflow = "ligand.complex.prepare"
            execution = {"placement": "local", "idempotency_key": secrets.token_urlsafe(24)}
            result = submit(workflow_schema(workflow), workflow,
                            {"source_pdb": str(source), "input_preparation_history_text": json.dumps(input_history), **parameters},
                            project, {"_execution": execution}, operation)
            if result and (run_id := (result.get("job") or {}).get("run_id")):
                task_group = result["job"].get("task_group") or "ligand-preparation"
                st.session_state[key + "-last-preparation-job"] = f"{task_group}/{run_id}"
        except (*ERRORS, RuntimeError) as exc:
            st.error(f"Could not submit preparation: {exc}")
    if reference := st.session_state.get(key + "-last-preparation-job"):
        if st.button("View preparation job", key=key + "-view-preparation"):
            st.session_state.pop("datasets-prepare-request", None)
            st.session_state["selected-job-reference"] = reference
            st.session_state["page"] = "Job details"
            st.rerun()
