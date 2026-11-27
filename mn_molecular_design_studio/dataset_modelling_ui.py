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

REPAIR = "PDBFixer repair"
EXTEND = "Extend the C-terminus"
SUBSTITUTE = "Substitute amino acids"
REBUILD = "Rebuild missing residues or a loop"
REMODEL = "Remodel an existing loop"

PICK_HINTS = {
    EXTEND: "Click a residue to choose the attachment point. An internal residue requires a trim.",
    SUBSTITUTE: "Click residues in the same chain to select or deselect them for substitution.",
    REBUILD: "Click a residue flanking the missing segment to select the gap.",
    REMODEL: "Click the first and last residues of the loop.",
}


@st.cache_data(max_entries=16, show_spinner=False)
def _gaps(source: str):
    return inspect_complex_gaps(Path(source))


def _gap_flanks(gap):
    return tuple((gap["chain"], str(gap[field])) for field in ("left_resseq", "right_resseq"))


def _available_gaps(path, chains):
    return {f"Chain {gap['chain']} · {gap['left_code']}{gap['left_resseq']}–{gap['right_code']}{gap['right_resseq']}": gap
            for gap in _gaps(str(path)) if gap["numbering_gap_length"] <= MAX_MODELED_GAP
            and all(label in chains.get(chain, {}) for chain, label in _gap_flanks(gap))}


def _gap_options(by_gap, key):
    if not by_gap:
        st.info(f"No internal numbering gaps of up to {MAX_MODELED_GAP} residues with observed flanks were found.")
        return {}, False, []
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
    return {"gaps": definitions}, valid, list(dict.fromkeys(site for label in selected for site in _gap_flanks(by_gap[label])))


def _substitutable_residues(residues):
    return {label: residue.resname for label, residue in residues.items()
            if all(atom in residue for atom in ("N", "CA", "C"))}


def _read_pick(viewer_key, state_key, chains):
    value = st.session_state.get(viewer_key) if viewer_key else None
    signature = (viewer_key, value)
    if not value or st.session_state.get(state_key) == signature:
        return None
    st.session_state[state_key] = signature
    try:
        pick = json.loads(value) if isinstance(value, str) else value
        pick = pick.get("residuePick") or {}
        chain = (pick["chainId"] or "_").strip() or "_"
        label = f"{int(pick['residueNumber'])}{pick['insertionCode'].strip()}"
        residue = chains.get(chain, {}).get(label)
        if residue is not None and residue.resname == pick["residueName"]:
            return chain, label
    except (AttributeError, KeyError, TypeError, ValueError):
        return None
    return None


def _apply_pick(site, operation, chains, gaps, key, viewer_key):
    chain, label = site
    residues = chains[chain]
    if operation == EXTEND:
        st.session_state[key + "-extension-chain"] = chain
        st.session_state[key + "-extension-site"] = (viewer_key, chain, label)
    elif operation == SUBSTITUTE:
        if label not in _substitutable_residues(residues):
            return "Choose a residue with observed N, CA and C atoms."
        st.session_state[key + "-substitution-chain"] = chain
        widget = f"{key}-{chain}-substitutions"
        selected = [item for item in st.session_state.get(widget, []) if item in residues]
        if label in selected:
            selected.remove(label)
        elif len(selected) < 50:
            selected.append(label)
        st.session_state[widget] = selected
    elif operation == REBUILD:
        matches = [name for name, gap in gaps.items() if site in _gap_flanks(gap)]
        if len(matches) != 1:
            return "Click the other flank of the intended gap, or select the gap in the list."
        selected = [name for name in st.session_state.get(key + "-gaps", []) if name in gaps]
        st.session_state[key + "-gaps"] = list(dict.fromkeys([*selected, *matches]))
    elif operation == REMODEL:
        labels = list(residues)[1:-1]
        if label not in labels:
            return "Choose a loop endpoint with an observed residue beyond it."
        start_key, end_key = f"{key}-loop-{chain}-start", f"{key}-loop-{chain}-end"
        start = st.session_state.get(start_key)
        if (st.session_state.get(key + "-loop-chain") == chain and start in labels
                and st.session_state.get(end_key) is None):
            start, end = sorted((start, label), key=labels.index)
        else:
            start, end = label, None
        st.session_state[key + "-loop-chain"] = chain
        st.session_state[start_key], st.session_state[end_key] = start, end
    return None


def _loop_options(chains, key):
    chains = {chain: residues for chain, residues in chains.items() if len(residues) >= 3}
    if not chains:
        st.info("An existing loop needs a protein chain with observed residues on both sides.")
        return {}, False, []
    chain = st.selectbox("Chain", list(chains), key=key + "-loop-chain")
    labels = list(chains[chain])[1:-1]
    chain_key = f"{key}-loop-{chain}"
    start = st.selectbox("First residue", labels, index=None, key=chain_key + "-start")
    end = st.selectbox("Last residue", labels, index=None, key=chain_key + "-end")
    models = st.number_input("MODELLER conformations", 1, MAX_LOOP_MODELS, 5, key=chain_key + "-models")
    selection = {"chain": chain, "start": start, "end": end, "models": models}
    sites = [(chain, label) for label in dict.fromkeys((start, end)) if label is not None]
    if start is None or end is None:
        return {}, False, sites
    chosen = loop_residues(chains, selection)
    st.caption("Selected sequence · " + " ".join(chains[chain][label].resname for label in chosen))
    return {"loop": selection}, True, sites


def _trim_to_attachment(key, chain, number):
    widget = f"{key}-{chain}"
    st.session_state[widget] = (st.session_state[widget][0], number)
    st.session_state.pop(key + "-extension-site", None)


def _extension_options(chains, key, viewer_key):
    chain = st.selectbox("Chain to extend", list(chains), key=key + "-extension-chain")
    residues = chains[chain]
    terminus = next(reversed(residues))
    candidate = st.session_state.get(key + "-extension-site")
    attachment = candidate[2] if candidate and candidate[:2] == (viewer_key, chain) and candidate[2] in residues else terminus
    st.caption(f"Append after residue {terminus} on the selected trim.")
    if attachment != terminus:
        number = residues[attachment].id[1]
        other_insertions = any(residue.id[1] == number and label != attachment for label, residue in residues.items())
        if other_insertions or residues[attachment].id[2].strip():
            st.info("Choose an attachment point with an unambiguous residue number.")
        elif f"{key}-{chain}" in st.session_state:
            st.button(f"Trim chain {chain} through residue {attachment}", key=key + "-apply-attachment",
                      help="Makes this residue the end of the chain in the preview.", on_click=_trim_to_attachment,
                      args=(key, chain, number))
            st.caption("Confirm the trim to attach here, or click the current terminus.")
    sequence = "".join(st.text_input("Amino acids to append", placeholder="e.g. GSG", max_chars=30,
                                     key=key + "-extension-sequence").upper().split())
    models = st.number_input("MODELLER conformations", min_value=1, max_value=50, value=10,
                             key=key + "-extension-models")
    if set(sequence) - set(AMINO_ACIDS):
        raise ValueError("Enter canonical one-letter amino-acid codes.")
    return {"extension_chain": chain, "extension_sequence": sequence, "extension_models": models}, (
        bool(sequence) and attachment == terminus), [(chain, attachment)]


def _substitution_options(chains, key):
    from Bio.SeqUtils import seq1, seq3

    chain = st.selectbox("Chain", list(chains), key=key + "-substitution-chain")
    residues = _substitutable_residues(chains[chain])
    if not residues:
        st.info("Substitution needs an observed amino acid with complete backbone coordinates.")
        return {}, False, []
    selected = st.multiselect("Residues to substitute", list(residues), max_selections=50,
                              format_func=lambda label: f"{residues[label]} {label}", key=f"{key}-{chain}-substitutions")
    substitutions = []
    for label in selected:
        target = st.selectbox(f"New amino acid for {residues[label]} {label}",
            [code for code in AMINO_ACIDS if code != seq1(residues[label])], index=None,
            format_func=lambda code: f"{seq3(code).upper()} ({code})", key=f"{key}-{chain}-{label}-substitution")
        if target:
            substitutions.append({"chain": chain, "residue": label, "amino_acid": target})
    return {"substitutions": substitutions}, bool(selected) and len(substitutions) == len(selected), [
        (chain, label) for label in selected]


def _operation_parameters(operation, chains, gaps, key, viewer_key):
    if operation == REPAIR:
        refine = st.checkbox("Refine rebuilt atoms with OpenMM", value=True, key=key + "-openmm",
                             help="Relaxes rebuilt atoms while holding the observed protein backbone fixed.")
        map_modified = st.checkbox("Map supported modified residues", value=True, key=key + "-map-modified")
        return {"repair": True, "refine_rebuilt_positions": refine, "preserve_original_backbone": True,
                "map_modified_residues": map_modified, "max_internal_gap": 0}, True, []
    if not chains:
        st.info("Keep a protein chain in the preview to use MODELLER.")
        return {}, False, []
    if operation == EXTEND:
        parameters, ready, sites = _extension_options(chains, key, viewer_key)
        return {"repair": False, **parameters}, ready, sites
    if operation == REBUILD:
        options, ready, sites = _gap_options(gaps, key)
    elif operation == REMODEL:
        options, ready, sites = _loop_options(chains, key)
    else:
        options, ready, sites = _substitution_options(chains, key)
    return {"repair": False, "modelling_options": parse_options(options) if ready else None}, ready, sites


def render_protein_tools(user, source_path, selected_pdb, input_history, key, projects, preferred_projects, submit,
                         *, viewer_key=""):
    st.markdown("**Repair and model the selected trim**")
    operation = st.selectbox("Preparation operation", [REPAIR, *PICK_HINTS], key=key + "-operation")
    parameters, ready, highlights = {}, False, []
    try:
        chains = ({chain: residues for chain, residues in observed_chains(selected_pdb.decode("utf-8")).items()
                   if residues} if operation in PICK_HINTS else {})
        gaps = _available_gaps(source_path, chains) if operation == REBUILD else {}
        site = _read_pick(viewer_key, key + "-last-pick", chains)
        if site and (message := _apply_pick(site, operation, chains, gaps, key, viewer_key)):
            st.info(message)
        parameters, ready, selected = _operation_parameters(operation, chains, gaps, key, viewer_key)
        for chain, label in selected:
            number, insertion = chains[chain][label].id[1:]
            highlights.append(f"{chain}:{number}^{insertion.strip() or '_'}")
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
    return PICK_HINTS.get(operation), highlights
