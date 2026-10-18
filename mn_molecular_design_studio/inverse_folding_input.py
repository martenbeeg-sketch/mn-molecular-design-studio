"""Portable structure-to-sequence contracts shared by tiles and pipeline nodes."""
from __future__ import annotations

import copy
import math
import re
import zipfile
from pathlib import Path

from mn_molecular_design_studio.mpnn_input import chain_list, residues
from mn_molecular_design_studio.generator_input import inventory
from mn_molecular_design_studio.runtime import reference_root

CALIBY_MODELS = {
    "caliby": "Caliby", "soluble_caliby": "SolubleCaliby", "soluble_caliby_v1": "SolubleCaliby v1",
    "caliby_distill": "Caliby-distill", "soluble_caliby_distill": "SolubleCaliby-distill",
}
ENGINES = {
    "esm_if1": dict(title="ESM-IF1", image="mn-esm-if1-opt:2.0.1-2b369911", version="ESM-IF1 · fair-esm 2.0.1 · 2b369911",
        description="Design protein sequences from a backbone with ESM-IF1.", source="https://github.com/facebookresearch/esm", reference="https://app.tamarind.bio/tools/esm-if1"),
    "caliby": dict(title="Caliby", image="mn-caliby-opt:0.1-41d31560", version="Caliby 0.1 · 41d31560",
        description="Potts model sequence design and scoring, optionally conditioned on structural ensembles.", source="https://github.com/ProteinDesignLab/caliby", reference="https://app.tamarind.bio/tools/caliby"),
    "boltzgen_inverse": dict(title="BoltzGen Inverse Folding", image="mn-boltzgen-opt:a07fea1-cu130-batch3", version="BoltzGen 0.3.2 · a07fea1",
        description="Redesign sequences on an existing protein structure from any generator or prediction engine.", source="https://github.com/HannesStark/boltzgen", reference="https://github.com/HannesStark/boltzgen#additional-configuration-options"),
}
WORKFLOWS = {engine: "sequence." + engine for engine in ENGINES}
AA = "ACDEFGHIKLMNPQRSTVWY"


def defaults(engine):
    return dict(design_chains="", designed_residues="", fixed_residues="", num_sequences=4 if engine == "caliby" else 10,
        temperature=.01 if engine == "caliby" else .1 if engine == "boltzgen_inverse" else 1., seed=0, optimization_mode="off", omit_aas="",
        task="design", model_type="caliby", ensemble_mode="generate" if engine == "caliby" else "single",
        num_conformers=32, ensemble_paths=None, constraints_csv=None, batch_size=1, num_sweeps=500,
        multichain_backbone=True)


def weights(engine, p=None, root=None):
    p = defaults(engine) | (p or {})
    root = Path(root) if root else reference_root()
    if engine == "esm_if1":
        return [root / "esm_models/esm_if1_gvp4_t16_142M_UR50.pt"]
    if engine == "caliby":
        paths = [root / "caliby/caliby" / (p["model_type"] + ".ckpt")]
        if p["ensemble_mode"] == "generate":
            paths += [root / "caliby/protpardelle-1c/weights/cc95_epoch3490.pth",
                      root / "caliby/protpardelle-1c/configs/cc95.yaml"]
            paths += list((root / "caliby/proteinmpnn").glob("*.pt")) or [root / "caliby/proteinmpnn"]
        return paths
    snapshots = root / "boltzgen-cache/models--boltzgen--boltzgen-1/snapshots"
    paths = sorted(snapshots.glob("*/boltzgen1_ifold.ckpt"))
    return paths[-1:] or [snapshots / "missing/boltzgen1_ifold.ckpt"]


def readiness(engine, p=None):
    missing = [str(v) for v in weights(engine, p) if not v.is_file()]
    return dict(ready=not missing, missing=missing)


def prepare(engine, parameters, *, require_input=False):
    if engine not in ENGINES:
        raise ValueError("Choose a supported inverse-folding engine.")
    p = defaults(engine) | copy.deepcopy(parameters)
    for name, low, high in [("num_sequences", 1, 10000), ("num_conformers", 1, 500),
                             ("batch_size", 1, 64), ("num_sweeps", 1, 100000), ("seed", 0, 2147483647)]:
        if type(p[name]) is not int or not low <= p[name] <= high:
            raise ValueError(f"{name.replace('_', ' ').capitalize()} must be an integer between {low} and {high}.")
    if isinstance(p["temperature"], bool) or not isinstance(p["temperature"], (float, int)) or not math.isfinite(p["temperature"]) or not 0 < p["temperature"] <= 5:
        raise ValueError("Temperature must be greater than zero and at most 5.")
    if p["task"] not in {"design", "score"} or p["task"] == "score" and engine != "caliby":
        raise ValueError("Scoring is available in Caliby.")
    if p["optimization_mode"] not in {"off", "fast", "exact"}:
        raise ValueError("Choose a supported acceleration mode.")
    if engine != "caliby" and p["optimization_mode"] != "off":
        raise ValueError("This structure-selection adapter currently uses the upstream inference mode.")
    if p["ensemble_mode"] not in {"single", "generate", "provided"}:
        raise ValueError("Choose single structure, generated ensembles, or provided ensembles.")
    if engine == "caliby":
        if p["model_type"] not in CALIBY_MODELS:
            raise ValueError("Choose a supported Caliby checkpoint.")
        if p["model_type"].endswith("_distill") and p["ensemble_mode"] != "single":
            raise ValueError("Distilled Caliby models use a single input structure.")
        if p["task"] == "score" and p["ensemble_mode"] == "generate":
            raise ValueError("For scoring, use the input structure or provide an ensemble.")
        if p["task"] == "score" and p["optimization_mode"] != "off":
            raise ValueError("Caliby scoring uses upstream inference mode.")
        if p["optimization_mode"] == "exact" and (p["ensemble_mode"] != "generate" or p["num_conformers"] != 32):
            raise ValueError("Caliby exact acceleration supports the generated 32-conformer protocol.")
        if p["optimization_mode"] == "fast" and p["ensemble_mode"] != "single":
            raise ValueError("Caliby fast acceleration supports single-structure design.")
    omit = re.sub(r"[,;\s]", "", str(p["omit_aas"] or "")).upper()
    if any(a not in AA for a in omit) or set(omit) == set(AA):
        raise ValueError("Omitted amino acids must be standard one-letter codes and leave at least one choice.")
    if engine == "esm_if1" and omit:
        raise ValueError("ESM-IF1 supports fixed positions; global amino-acid omissions are not exposed by this adapter.")
    p["omit_aas"] = "".join(dict.fromkeys(omit))
    p["design_chains"] = ",".join(chain_list(p["design_chains"]))
    direct = [Path(v) for v in (p.get("structure_paths") or [])]
    if p.get("structure_path"):
        direct.insert(0, Path(p["structure_path"]))
    direct = list(dict.fromkeys(direct))
    source = bool(p.get("source_run_dir") and p.get("candidates_jsonl"))
    if bool(p.get("source_run_dir")) != bool(p.get("candidates_jsonl")):
        raise ValueError("Select a source run and its candidate collection together.")
    if source:
        root = Path(p["source_run_dir"]).resolve()
        candidate_path = Path(p["candidates_jsonl"]).resolve()
        if not candidate_path.is_relative_to(root):
            raise ValueError("The candidate collection must belong to the selected source run.")
        if not candidate_path.is_file():
            raise ValueError("The source candidate collection is missing.")
    if direct and source:
        raise ValueError("Choose structure inputs or an existing candidate collection.")
    if require_input and not direct and not source:
        raise ValueError("Upload a structure or select existing results.")
    if source and p.get("selected_candidate_ids") == []:
        raise ValueError("Select at least one input candidate.")
    if engine == "caliby" and p["ensemble_mode"] == "provided":
        if not p.get("ensemble_paths"):
            raise ValueError("Upload conformer structures for the provided ensemble.")
        if len(direct) != 1 or source:
            raise ValueError("A provided ensemble must be associated with one primary structure.")
    for path in direct + [Path(v) for v in (p.get("ensemble_paths") or [])]:
        if not path.is_file() or path.suffix.lower() not in {".pdb", ".ent", ".cif", ".mmcif"}:
            raise ValueError(f"Structure input must be an existing PDB or mmCIF file: {path}")
    for path in direct:
        info = inventory(path)
        chains = chain_list(p["design_chains"]) or [c for c, v in info["chains"].items() if v["kind"] == "protein"]
        if not chains or any(c not in info["chains"] or info["chains"][c]["kind"] != "protein" for c in chains):
            raise ValueError("Choose protein chains present in every input structure.")
        for key in ("designed_residues", "fixed_residues"):
            residues(p[key], info, chains=chains)
    if p.get("constraints_csv") and (engine != "caliby" or not Path(p["constraints_csv"]).is_file()):
        raise ValueError("Upload a Caliby position-constraints CSV file.")
    p["structure_paths"] = [str(v) for v in direct]
    return p


def expand_structure_archive(archive, destination):
    """Extract only bounded structure entries, preserving no archive paths."""
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive) as handle:
        entries = [v for v in handle.infolist() if not v.is_dir() and Path(v.filename).suffix.lower() in {".pdb", ".ent", ".cif", ".mmcif"}]
        if not entries or len(entries) > 1000 or sum(v.file_size for v in entries) > 200 * 1024 * 1024:
            raise ValueError("The ZIP must contain 1–1000 structure files totaling at most 200 MB.")
        paths = []
        for index, entry in enumerate(entries):
            path = destination / f"{index:04d}-{Path(entry.filename).name}"
            path.write_bytes(handle.read(entry))
            paths.append(path)
    return paths


def canvas_fields(engine):
    from mn_molecular_design_studio.pipeline_graph import field
    p = defaults(engine)
    names = ["num_sequences", "design_chains", "designed_residues", "fixed_residues", "temperature", "seed"]
    if engine == "esm_if1":
        names += ["multichain_backbone"]
    else:
        names += ["omit_aas"]
    if engine == "caliby":
        names += ["task", "model_type", "ensemble_mode", "num_conformers", "num_sweeps", "batch_size", "optimization_mode"]
    options = dict(task=["design", "score"], model_type=list(CALIBY_MODELS), ensemble_mode=["single", "generate"], optimization_mode=["off", "fast", "exact"])
    result = []
    for name in names:
        value = p[name]
        kind = "select" if name in options else "boolean" if isinstance(value, bool) else "text" if isinstance(value, str) else "number"
        result.append(field(name, name.replace("_", " ").capitalize(), value, kind=kind,
            options=options.get(name), section="basic" if name == "num_sequences" else "optional"))
    return result
