"""Validate user-selected regions for MODELLER structure preparation."""
from __future__ import annotations

from copy import deepcopy
import re

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
MAX_MODELED_GAP = 15
MAX_MODELED_LOOP = 30
MAX_LOOP_MODELS = 20


def parse_options(value: dict | None) -> dict:
    if value is None:
        return {}
    if not isinstance(value, dict) or set(value) - {"gaps", "loop", "substitutions"}:
        raise ValueError("MODELLER settings must specify missing residues, a loop, or amino-acid substitutions.")
    options = deepcopy(value)
    gaps = options.get("gaps", [])
    if not isinstance(gaps, list) or len(gaps) > 50:
        raise ValueError("Choose at most 50 missing segments.")
    for gap in gaps:
        if not isinstance(gap, dict) or gap.get("confirmed") is not True:
            raise ValueError("Confirm the missing sequence and both flanking residues.")
        chain, left, right = (gap.get(key) for key in ("chain", "left_resseq", "right_resseq"))
        if not _valid_chain(chain) or type(left) is not int or type(right) is not int:
            raise ValueError("A missing segment requires a chain and two observed residue numbers.")
        sequence = "".join(str(gap.get("sequence") or "").upper().split())
        if not 1 <= len(sequence) <= MAX_MODELED_GAP or set(sequence) - set(AMINO_ACIDS):
            raise ValueError(f"Enter up to {MAX_MODELED_GAP} canonical amino-acid letters for a missing segment.")
        if len(sequence) != right - left - 1:
            raise ValueError("The missing sequence must match the number of positions between the selected flanks.")
        gap["sequence"] = sequence
    if "loop" in options:
        loop = options["loop"]
        if not isinstance(loop, dict) or not _valid_chain(loop.get("chain")):
            raise ValueError("Choose the chain of the loop to remodel.")
        if any(not isinstance(loop.get(key), str) or not re.fullmatch(r"-?\d{1,4}[A-Za-z]?", loop[key])
               for key in ("start", "end")):
            raise ValueError("Choose observed start and end residues for the loop.")
        models = loop.get("models", 5)
        if type(models) is not int or not 1 <= models <= MAX_LOOP_MODELS:
            raise ValueError(f"Choose between 1 and {MAX_LOOP_MODELS} loop conformations.")
        loop["models"] = models
    substitutions = options.get("substitutions", [])
    if not isinstance(substitutions, list) or len(substitutions) > 50:
        raise ValueError("Select at most 50 amino-acid substitutions.")
    selected = set()
    for item in substitutions:
        if not isinstance(item, dict) or not _valid_chain(item.get("chain")):
            raise ValueError("Select a chain and an observed residue for each substitution.")
        residue = item.get("residue")
        if not isinstance(residue, str) or not re.fullmatch(r"-?\d{1,4}[A-Za-z]?", residue):
            raise ValueError("Choose an observed residue number for each substitution.")
        amino_acid = str(item.get("amino_acid") or "").strip().upper()
        if len(amino_acid) != 1 or amino_acid not in AMINO_ACIDS:
            raise ValueError("Choose a canonical amino acid for each substitution.")
        site = (item["chain"], residue)
        if site in selected:
            raise ValueError("Choose a single substitution for each residue.")
        selected.add(site)
        item["amino_acid"] = amino_acid
    if not gaps and "loop" not in options and not substitutions:
        raise ValueError("Select missing residues, an existing loop, or amino-acid substitutions.")
    return options


def _valid_chain(chain):
    return isinstance(chain, str) and re.fullmatch(r"[A-Za-z0-9_]", chain) is not None


def observed_chains(source: str):
    from mn_molecular_design_studio.molecule_imports import read_coordinate_structure
    from Bio.PDB.Polypeptide import is_aa

    models = list(read_coordinate_structure(source.encode("utf-8"), "source.pdb").get_models())
    if len(models) != 1:
        raise ValueError("Select a PDB with one coordinate model for MODELLER.")
    return {chain.id.strip() or "_": {f"{residue.id[1]}{residue.id[2].strip()}": residue
                                     for residue in chain if is_aa(residue, standard=True)}
            for chain in models[0] if _valid_chain(chain.id.strip() or "_")}


def substitution_sites(chains, substitutions):
    from Bio.SeqUtils import seq3

    sites = []
    for substitution in substitutions:
        chain, label = substitution["chain"], substitution["residue"]
        residue = chains.get(chain, {}).get(label)
        if residue is None or any(atom not in residue for atom in ("N", "CA", "C")):
            raise ValueError(f"Chain {chain} residue {label} needs observed backbone coordinates for substitution.")
        resname = residue.resname
        target = seq3(substitution["amino_acid"]).upper()
        if resname == target:
            raise ValueError(f"Chain {chain} residue {label} is already {target}.")
        number, insertion = str(residue.id[1]), residue.id[2].strip()
        sites.append({"key": "|".join((resname, chain, number, insertion or "_")), "chain": chain,
                      "resseq": number, "icode": insertion, "resname": resname, "target": target})
    return sites


def loop_residues(chains, loop):
    residues = chains.get(loop["chain"], {})
    labels = list(residues)
    try:
        start, end = labels.index(loop["start"]), labels.index(loop["end"])
    except ValueError as exc:
        raise ValueError("The loop endpoints must exist in the selected chain.") from exc
    if not 0 < start <= end < len(labels) - 1 or end - start + 1 > MAX_MODELED_LOOP:
        raise ValueError(f"Select an internal loop of up to {MAX_MODELED_LOOP} residues with an observed flank at each end.")
    window = [residues[label] for label in labels[start - 1:end + 2]]
    if any(any(atom not in residue for atom in ("N", "CA", "C")) for residue in window):
        raise ValueError("The loop and its flanks need complete backbone coordinates.")
    if any(left["C"] - right["N"] > 2.0 for left, right in zip(window, window[1:])):
        raise ValueError("The selected region crosses a backbone break. Rebuild the missing residues first.")
    return labels[start:end + 1]
