"""Shared model, checkpoint and input contract for the MPNN tile and canvas."""
from __future__ import annotations

import copy
import json
import math
import re
from pathlib import Path

from mn_molecular_design_studio.runtime import reference_root
from mn_molecular_design_studio.generator_input import inventory

WORKFLOW = "sequence.mpnn"
IMAGE = "mn-ligandmpnn:latest"
ASSETS = Path(__file__).with_name("assets") / "mpnn"
AA = "ACDEFGHIKLMNPQRSTVWYX"
MODELS = {
    "protein_mpnn": {"label": "ProteinMPNN", "description": "Standard model for general protein sequence design.", "native": "protein_mpnn", "noise": [.02, .1, .2, .3], "default": .2},
    "ligand_mpnn": {"label": "LigandMPNN", "description": "Design sequences using the atomic context of ligands, nucleic acids and metals.", "native": "ligand_mpnn", "noise": [.05, .1, .2, .3], "default": .1},
    "soluble_mpnn": {"label": "SolubleMPNN", "description": "ProteinMPNN trained on soluble proteins.", "native": "soluble_mpnn", "noise": [.02, .1, .2, .3], "default": .2},
    "hyper_mpnn": {"label": "HyperMPNN", "description": "ProteinMPNN weights trained on hyperthermophiles for thermostable protein design.", "native": "protein_mpnn", "noise": [.02, .1, .2, .3], "default": .2},
    "ab_mpnn": {"label": "AbMPNN", "description": "ProteinMPNN weights specialized for antibody sequence design.", "native": "protein_mpnn", "noise": [], "default": None},
    "global_label_membrane_mpnn": {"label": "MembraneMPNN · global", "description": "Condition sequence design on whether the whole protein is transmembrane.", "native": "global_label_membrane_mpnn", "noise": [.2], "default": .2},
    "per_residue_label_membrane_mpnn": {"label": "MembraneMPNN · residues", "description": "Label membrane-buried and membrane-interface residues individually.", "native": "per_residue_label_membrane_mpnn", "noise": [.2], "default": .2},
}


def defaults():
    return dict(model_type="protein_mpnn", noise_level=None, num_seq_per_target=2, sampling_temp=.1,
                design_chains="", designed_residues="", fixed_residues="", omit_aas="C", seed=None,
                fix_interface_residues=False,
                bias_aa="", bias_aa_per_residue={}, omit_aa_per_residue={}, homo_oligomer=0,
                global_transmembrane_label=True, transmembrane_buried="",
                transmembrane_interface="", use_side_chain_context=False, pack_side_chains=False,
                number_of_packs_per_design=1, backend="auto")


def checkpoint(model_type="protein_mpnn", noise_level=None, *, root=None):
    if model_type not in MODELS:
        raise ValueError("Choose a supported MPNN model.")
    spec = MODELS[model_type]
    if not spec["noise"] and noise_level is not None:
        raise ValueError(f"{spec['label']} provides one checkpoint and has no noise-level selector.")
    noise = spec["default"] if noise_level is None else float(noise_level)
    if spec["noise"] and noise not in spec["noise"]:
        raise ValueError(f"{spec['label']} has no checkpoint trained at {noise:g} Å noise.")
    root = Path(root) if root is not None else reference_root()
    if model_type == "ab_mpnn":
        return root / "mpnn" / "abmpnn.pt"
    level = round(noise * 100)
    if model_type == "hyper_mpnn":
        epoch = 240 if level == 2 else 300
        return root / "mpnn" / f"v48_{level:03d}_epoch{epoch}_hyper.pt"
    base = root / "protpardelle-1c" / "model_params" / "LigandMPNN"
    prefix = {"protein_mpnn": "proteinmpnn", "ligand_mpnn": "ligandmpnn", "soluble_mpnn": "solublempnn"}.get(model_type)
    if prefix:
        return base / (f"{prefix}_v_{32 if model_type == 'ligand_mpnn' else 48}_{level:03d}" + ("_25" if model_type == "ligand_mpnn" else "") + ".pt")
    return base / f"{model_type}_v_48_020.pt"


def chain_list(value):
    return list(dict.fromkeys(str(v) for v in value)) if isinstance(value, (list, tuple)) else list(dict.fromkeys(re.split(r"[,;\s]+", str(value or "").strip()))) if str(value or "").strip() else []


def residues(value, info=None, *, chains=None):
    """Accept explicit author IDs, insertion codes, and inclusive numeric ranges."""
    result = []
    for token in re.split(r"[,;\s]+", str(value or "").strip()):
        if not token:
            continue
        match = re.fullmatch(r"([A-Za-z0-9]+?)(-?\d+)([A-Za-z]?)(?:-(?:([A-Za-z]+))?(-?\d+))?", token)
        if not match:
            raise ValueError(f"Invalid residue '{token}'; use A12, A12A or A12-20.")
        chain, number, insertion, end_chain, end = match.groups()
        if chains is not None and chain not in chains:
            raise ValueError(f"{token} is outside the designed chains.")
        start = int(number)
        if end is not None:
            if insertion or end_chain and end_chain != chain or not start <= int(end) <= start + 10000:
                raise ValueError(f"Invalid residue range: {token}.")
            labels = [f"{chain}{n}" for n in range(start, int(end) + 1)]
        else:
            labels = [f"{chain}{start}{insertion}"]
        if info and any(label not in info["residues"] or info["residues"][label]["kind"] != "protein" for label in labels):
            raise ValueError(f"{token} includes residues absent from the input protein.")
        result.extend(labels)
    return list(dict.fromkeys(result))


def _json_map(value, label):
    if value in (None, ""):
        return {}
    try:
        result = json.loads(value) if isinstance(value, str) else copy.deepcopy(value)
    except ValueError as exc:
        raise ValueError(f"{label} must be a JSON object.") from exc
    if not isinstance(result, dict):
        raise ValueError(f"{label} must be a JSON object.")
    return result


def prepare(parameters, *, require_input=False):
    p = defaults() | copy.deepcopy(parameters)
    for key in ("structure_path", "source_run_dir", "candidates_jsonl"):
        if isinstance(p.get(key), Path):
            p[key] = str(p[key])
    model = p["model_type"]
    if not isinstance(p["fix_interface_residues"], bool):
        raise ValueError("Fix interface residues must be enabled or disabled.")
    checkpoint(model, p["noise_level"])
    p["noise_level"] = MODELS[model]["default"] if p["noise_level"] is None else p["noise_level"]
    for key, lo, hi in [("num_seq_per_target", 1, 10000), ("homo_oligomer", 0, 62), ("number_of_packs_per_design", 0, 100)]:
        if isinstance(p[key], bool) or not isinstance(p[key], int) or not lo <= p[key] <= hi:
            raise ValueError(f"{key.replace('_', ' ').capitalize()} must be an integer between {lo} and {hi}.")
    if p["seed"] == 0:
        p["seed"] = None
    if p["seed"] is not None and (isinstance(p["seed"], bool) or not isinstance(p["seed"], int) or not 1 <= p["seed"] <= 4294967295):
        raise ValueError("Seed must be between 1 and 4294967295, or blank for a random seed.")
    if not isinstance(p["sampling_temp"], (float, int)) or isinstance(p["sampling_temp"], bool) or not math.isfinite(p["sampling_temp"]) or not 0 < p["sampling_temp"] <= 1:
        raise ValueError("Temperature must be greater than 0 and at most 1.")
    p["design_chains"] = ",".join(chain_list(p["design_chains"]))
    omit = str(p["omit_aas"] or "").strip().upper()
    p["omit_aas"] = "" if omit in {"", "EMPTY"} else omit
    if any(a not in AA for a in p["omit_aas"]) or set(p["omit_aas"]) >= set(AA.replace("X", "")):
        raise ValueError("Amino acids to omit must leave at least one standard amino acid available.")
    bias = {}
    raw_bias = p["bias_aa"]
    if raw_bias:
        if isinstance(raw_bias, str):
            try:
                raw_bias = dict(token.strip().split(":") for token in raw_bias.split(",") if token.strip())
            except ValueError as exc:
                raise ValueError("Bias AA must use entries such as W:3.0,P:-1.0.") from exc
        if not isinstance(raw_bias, dict):
            raise ValueError("Bias AA must be an amino-acid to number map.")
        for aa, score in raw_bias.items():
            aa = aa.strip().upper()
            if len(aa) != 1 or aa not in AA or isinstance(score, bool) or not math.isfinite(float(score)):
                raise ValueError("Bias AA must contain amino-acid letters and finite numeric values.")
            bias[aa] = float(score)
    p["bias_aa"] = ",".join(f"{aa}:{score:g}" for aa, score in bias.items())
    info = inventory(Path(p["structure_path"])) if p.get("structure_path") else None
    chains = chain_list(p["design_chains"])
    if info:
        available = {c: v for c, v in info["chains"].items() if v["kind"] == "protein"}
        if not chains:
            raise ValueError("Select at least one designed chain.")
        if not set(chains) <= set(available):
            raise ValueError("Every designed chain must be a protein chain in the input structure.")
        if p["homo_oligomer"] > 1 and (len(chains) != p["homo_oligomer"] or len({len(available[c]["residues"]) for c in chains}) != 1):
            raise ValueError("Homo oligomer copies must match equally sized designed chains already present in the structure.")
    for key in ("designed_residues", "fixed_residues", "transmembrane_buried", "transmembrane_interface"):
        p[key] = " ".join(residues(p[key], info, chains=chains or None))
    if set(p["designed_residues"].split()) & set(p["fixed_residues"].split()):
        raise ValueError("A residue cannot be both designed and fixed.")
    if set(p["transmembrane_buried"].split()) & set(p["transmembrane_interface"].split()):
        raise ValueError("A residue cannot be both membrane-buried and at the membrane interface.")
    for key in ("bias_aa_per_residue", "omit_aa_per_residue"):
        p[key] = _json_map(p[key], key.replace("_", " "))
        for position, value in p[key].items():
            if len(residues(position, info, chains=chains or None)) != 1:
                raise ValueError("Per-residue settings require individual residue IDs.")
            if key == "bias_aa_per_residue":
                if not isinstance(value, dict) or any(len(aa) != 1 or aa not in AA or isinstance(score, bool) or not math.isfinite(float(score)) for aa, score in value.items()):
                    raise ValueError("Per-residue biases require amino-acid to finite number maps.")
            elif not isinstance(value, str) or any(a not in AA for a in value) or set(value) >= set(AA.replace("X", "")):
                raise ValueError("Per-residue omissions must leave at least one standard amino acid available.")
    if p["backend"] not in {"auto", "shared_ligandmpnn", "foundry"}:
        raise ValueError("Choose Automatic, Shared MPNN or Foundry as the backend.")
    if p["backend"] == "foundry" and (model not in {"protein_mpnn", "ligand_mpnn"} or p["pack_side_chains"]):
        raise ValueError("This model or side-chain packing option requires the Shared MPNN backend; choose Automatic.")
    if p["use_side_chain_context"] and model != "ligand_mpnn":
        raise ValueError("Atomic side-chain context is available with LigandMPNN.")
    if p["pack_side_chains"] and p["number_of_packs_per_design"] < 1:
        raise ValueError("Choose at least one packing sample when side-chain packing is enabled.")
    if require_input and not (p.get("structure_path") or p.get("source_run_dir") and p.get("candidates_jsonl")):
        raise ValueError("Upload a structure, fetch a PDB, or select existing results.")
    if p.get("structure_path") and (p.get("source_run_dir") or p.get("candidates_jsonl")):
        raise ValueError("Choose either a structure input or existing results.")
    if p.get("source_run_dir") and p.get("candidates_jsonl"):
        if not Path(p["candidates_jsonl"]).resolve().is_relative_to(Path(p["source_run_dir"]).resolve()):
            raise ValueError("The candidate collection must belong to the selected source job.")
    return p


def readiness(model_type="protein_mpnn", noise_level=None):
    path = checkpoint(model_type, noise_level)
    return {"ready": path.is_file() and path.stat().st_size > 1000000, "checkpoint": str(path)}
