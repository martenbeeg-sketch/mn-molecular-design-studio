"""Standalone small-molecule tools exposed as Studio workflow tiles."""

from __future__ import annotations

import csv
import ast
import json
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
import numpy as np

from mn_molecular_design_studio.engine.workflow_engine.core.docker_runner import (
    run_docker_job,
)
from mn_molecular_design_studio.engine.workflow_engine.core.gpu import docker_gpu_args
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import (
    create_job,
    finish_job,
    read_json,
    update_status,
    write_json,
)
from mn_molecular_design_studio.engine.workflow_engine.core.manifests import load_manifest
from mn_molecular_design_studio.engine.workflow_engine.runtime import reference_root
from mn_molecular_design_studio import ligand_compounds


MAX_COMPOUNDS_PER_LAUNCH = 150_000
_MAX_COMPOUNDS_PER_LAUNCH = MAX_COMPOUNDS_PER_LAUNCH
_LDDM_CHECKPOINT = "generation/lddm/lddm_CDBB.ckpt"
_POCKETXMOL_CHECKPOINT = "generation/pocketxmol/data/trained_models/pxm_use/checkpoints/pocketxmol.ckpt"
_POCKETXMOL_CONFIG = "generation/pocketxmol/data/trained_models/pxm_use/train_config/train.yml"
_GENERATOR_CONFIGS: dict[str, dict[str, Any]] = {
    "omtra": {
        "name": "OMTRA", "tool": "omtra_generation", "image": "mn-omtra:cu128-745127d",
        "references": ("generation/omtra/checkpoint.ckpt",),
        "conditions": ("target_pdb",),
    },
    "flowr_root": {
        "name": "FLOWR.root", "tool": "flowr_root_generation", "image": "mn-flowr-root:cu128-b2263e2",
        "references": ("generation/flowr_root/flowr_root_v2.2.ckpt",),
        "conditions": ("target_pdb", "pocket_pdb", "reference_ligand_sdf"),
    },
    "conditar": {
        "name": "conDitar", "tool": "conditar_generation", "image": "mn-conditar:cu128-4294d286",
        "references": ("generation/conditar/Diff.pt", "generation/conditar/PocketAE.pt"),
        "conditions": ("target_pdb", "pocket_pdb"),
    },
    "paopt": {
        "name": "conDitar + paOPT", "tool": "paopt_generation", "image": "mn-paopt:cu128-4294d286",
        "references": ("generation/conditar/Diff.pt", "generation/conditar/PocketAE.pt"),
        "conditions": ("target_pdb", "pocket_pdb", "reference_ligand_sdf"),
    },
    "drugrpg": {
        "name": "DrugRPG", "tool": "drugrpg_generation", "image": "mn-drugrpg:cu128-6fa0e41",
        "references": ("generation/drugrpg/training.pt",),
        "conditions": ("target_pdb", "pocket_pdb"),
    },
    "pfm": {
        "name": "PFM", "tool": "pfm_generation", "image": "mn-pfm:cu128-33be6c1",
        "references": (
            "generation/pfm/pfm.pth", "generation/pfm/x_predictor.pth",
            "generation/pfm/h_predictor.pth", "generation/pfm/training.yml",
        ),
        "conditions": ("target_pdb", "pocket_pdb"),
    },
    "pocketflow": {
        "name": "PocketFlow", "tool": "pocketflow_generation", "image": "mn-pocketflow:cu128-a31a5a0",
        "references": ("generation/pocketflow/ZINC-pretrained-255000.pt",),
        "conditions": ("target_pdb", "pocket_pdb"),
    },
    "pgmg": {
        "name": "PGMG", "tool": "pgmg_generation", "image": "mn-pgmg:cu128-85fb712",
        "references": ("generation/pgmg/model.pth", "generation/pgmg/tokenizer.pkl"),
        "conditions": ("target_pdb", "pharmacophore_posp"),
    },
}
_ROSETTALIGAND_IMAGE = "mn-openvs:local"
_GROMACS_MD_IMAGE = "mn-gromacs:2026.3-cu128"
_OPENMM_MD_IMAGE = "mn-md:cu128"
_GROMACS_TILE_RUNNER = r'''from __future__ import annotations
import json, sys, traceback
from pathlib import Path
from mn_ligand.workflows.gromacs_md import prepare_gromacs_system, run_gromacs_production

def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))

def write_json(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

root = Path("/work/artifacts/md")
config_path = Path("/work/artifacts/input/gromacs_config.json")
result_path = root / "result.json"
config = read_json(config_path)
try:
    prep_path = root / "preparation" / "result.json"
    prep_result = prepare_gromacs_system(config, prep_path)
    prep_files = ((prep_result.get("md_result") or {}).get("output_files") or {})
    if not prep_result.get("success") or not prep_files:
        raise RuntimeError("GROMACS preparation did not return system files")
    production_config = {
        **config,
        "source_topology_path": prep_files["gromacs_topology"],
        "source_coordinates_path": prep_files["gromacs_coordinates"],
        "source_checkpoint_path": prep_files["gromacs_checkpoint"],
        "source_index_path": prep_files["gromacs_index"],
    }
    production_path = root / "production" / "result.json"
    production_result = run_gromacs_production(production_config, production_path)
    if not production_result.get("success"):
        raise RuntimeError("GROMACS production did not complete successfully")
    write_json(result_path, {
        "success": True,
        "engine": "GROMACS",
        "selected_ligand": config.get("ligand_key"),
        "production_ns": config.get("production_length_ns"),
        "preparation": prep_result.get("md_result"),
        "production": production_result.get("md_result"),
        "output_files": (production_result.get("md_result") or {}).get("output_files") or {},
    })
except Exception as exc:
    write_json(result_path, {"success": False, "error": str(exc), "traceback": traceback.format_exc()})
    raise
'''
_OPENMM_TILE_RUNNER = r'''from __future__ import annotations
import json, traceback
from pathlib import Path
from mn_ligand.workflows.bound_ligand_md import run_ligandx_md

root = Path("/work/artifacts/md")
config = json.loads(Path("/work/artifacts/input/openmm_config.json").read_text(encoding="utf-8"))
result_path = root / "result.json"
try:
    result = run_ligandx_md(config, result_path)
    if not result.get("success"):
        raise RuntimeError(str(result.get("error") or (result.get("md_result") or {}).get("error") or "OpenMM MD did not complete successfully"))
except Exception as exc:
    result_path.write_text(json.dumps({"success": False, "error": str(exc), "traceback": traceback.format_exc()}, indent=2) + "\\n", encoding="utf-8")
    raise
'''
_RDKIT_3D_SCRIPT = r'''from __future__ import annotations
import argparse, math
from rdkit import Chem
from rdkit.Chem import AllChem
p = argparse.ArgumentParser()
p.add_argument("--smiles", required=True)
p.add_argument("--name", required=True)
p.add_argument("--output", required=True)
a = p.parse_args()
m = Chem.MolFromSmiles(a.smiles, sanitize=True)
if m is None:
    raise SystemExit("RDKit could not sanitize the input SMILES")
m = Chem.AddHs(m)
params = AllChem.ETKDGv3()
params.randomSeed = 0x4D4E
status = AllChem.EmbedMolecule(m, params)
if status != 0:
    params.useRandomCoords = True
    status = AllChem.EmbedMolecule(m, params)
if status != 0 or m.GetNumConformers() != 1:
    raise SystemExit("RDKit could not generate a 3D conformer")
try:
    if AllChem.MMFFHasAllMoleculeParams(m):
        AllChem.MMFFOptimizeMolecule(m, maxIters=300)
    else:
        AllChem.UFFOptimizeMolecule(m, maxIters=300)
except (RuntimeError, ValueError):
    pass
c = m.GetConformer()
if any(not math.isfinite(v) for i in range(m.GetNumAtoms()) for v in tuple(c.GetAtomPosition(i))):
    raise SystemExit("RDKit generated non-finite coordinates")
m.SetProp("_Name", a.name)
w = Chem.SDWriter(a.output)
w.write(m)
w.close()
'''


def _safe_id(value: object, index: int) -> str:
    text = re.sub(r"[^A-Za-z0-9._-]+", "-", str(value or "")).strip("-._")
    return (text or f"compound-{index:06d}")[:100]


def _validated_compound_set(path: Path, compound_limit: int) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    if type(compound_limit) is not int or not 0 <= compound_limit <= _MAX_COMPOUNDS_PER_LAUNCH:
        raise ValueError(f"compound_limit must be between 0 (all selected) and {_MAX_COMPOUNDS_PER_LAUNCH}.")
    payload = ligand_compounds.load_compound_set(Path(path))
    compounds = payload["compounds"][:compound_limit] if compound_limit else payload["compounds"]
    if len(compounds) > _MAX_COMPOUNDS_PER_LAUNCH:
        raise ValueError(f"Select at most {_MAX_COMPOUNDS_PER_LAUNCH} compounds for this run; no compounds were silently omitted.")
    if not compounds:
        raise ValueError("The selected compound set is empty.")
    return payload, compounds


def validate_target_pdb(path: Path) -> None:
    target = Path(path)
    if target.suffix.lower() not in {".pdb", ".ent"}:
        raise ValueError("Ligand docking and Nesso currently require a PDB or ENT target structure.")
    try:
        text = target.read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        raise ValueError(f"Could not read target structure: {exc}") from exc
    if not any(line.startswith(("ATOM  ", "HETATM")) for line in text.splitlines()):
        raise ValueError("The target PDB contains no atom records.")


def validate_ligand_sdf(path: Path) -> None:
    ligand = Path(path)
    if ligand.suffix.lower() not in {".sdf", ".mol"}:
        raise ValueError("The reference ligand must be an SDF or MOL file.")
    try:
        text = ligand.read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        raise ValueError(f"Could not read reference ligand: {exc}") from exc
    records = [record for record in text.split("$$$$") if record.strip()]
    if not records or not any("M  END" in record for record in records):
        raise ValueError("The reference ligand contains no readable SDF/MOL record.")


def _sdf_heavy_atom_coordinates(path: Path) -> np.ndarray:
    """Read heavy-atom coordinates from a V2000/V3000 SDF without RDKit."""
    validate_ligand_sdf(path)
    lines = Path(path).read_text(encoding="utf-8", errors="replace").split("$$$$", 1)[0].splitlines()
    coordinates: list[tuple[float, float, float]] = []
    if any("V3000" in line for line in lines[:5]):
        in_atoms = False
        for line in lines:
            if "M  V30 BEGIN ATOM" in line:
                in_atoms = True
            elif "M  V30 END ATOM" in line:
                break
            elif in_atoms and "M  V30 " in line:
                fields = line.split()
                if len(fields) < 7 or fields[3].upper() in {"H", "D"}:
                    continue
                try:
                    coordinates.append(tuple(float(value) for value in fields[4:7]))
                except ValueError as exc:
                    raise ValueError("The reference SDF has invalid atom coordinates.") from exc
    else:
        if len(lines) < 4:
            raise ValueError("The reference SDF is missing its atom-count header.")
        try:
            atom_count = int(lines[3][:3])
        except ValueError as exc:
            raise ValueError("The reference SDF has an invalid atom-count header.") from exc
        if atom_count < 1 or len(lines) < 4 + atom_count:
            raise ValueError("The reference SDF atom records are incomplete.")
        for line in lines[4 : 4 + atom_count]:
            element = line[31:34].strip().upper()
            if element in {"H", "D"}:
                continue
            try:
                coordinates.append(
                    (float(line[0:10]), float(line[10:20]), float(line[20:30]))
                )
            except ValueError as exc:
                raise ValueError("The reference SDF has invalid atom coordinates.") from exc
    if len(coordinates) < 2:
        raise ValueError("Axis alignment requires at least two ligand heavy atoms.")
    result = np.asarray(coordinates, dtype=float)
    if not np.isfinite(result).all():
        raise ValueError("The reference ligand contains non-finite coordinates.")
    return result


def _rotation_from_vectors(source: np.ndarray, destination: np.ndarray) -> np.ndarray:
    source = source / np.linalg.norm(source)
    destination = destination / np.linalg.norm(destination)
    cross = np.cross(source, destination)
    sine = float(np.linalg.norm(cross))
    cosine = float(np.dot(source, destination))
    if sine < 1e-12:
        if cosine > 0:
            return np.eye(3)
        perpendicular = np.array([0.0, 1.0, 0.0])
        if abs(float(np.dot(source, perpendicular))) > 0.9:
            perpendicular = np.array([0.0, 0.0, 1.0])
        axis = np.cross(source, perpendicular)
        axis /= np.linalg.norm(axis)
        return 2.0 * np.outer(axis, axis) - np.eye(3)
    skew = np.array(
        [[0.0, -cross[2], cross[1]], [cross[2], 0.0, -cross[0]], [-cross[1], cross[0], 0.0]]
    )
    return np.eye(3) + skew + (skew @ skew) * ((1.0 - cosine) / (sine * sine))


def ligand_longest_axis_transform(path: Path) -> dict[str, Any]:
    """Rotate a coordinate ligand's longest axis onto global X."""
    coordinates = _sdf_heavy_atom_coordinates(path)
    pivot = coordinates.mean(axis=0)
    centered = coordinates - pivot
    eigenvalues, eigenvectors = np.linalg.eigh(centered.T @ centered)
    axis = eigenvectors[:, int(np.argmax(eigenvalues))]
    dominant = int(np.argmax(np.abs(axis)))
    if axis[dominant] < 0:
        axis = -axis
    rotation = _rotation_from_vectors(axis, np.array([1.0, 0.0, 0.0]))
    aligned = centered @ rotation.T + pivot
    minimum, maximum = aligned.min(axis=0), aligned.max(axis=0)
    return {
        "pivot": [float(value) for value in pivot],
        "rotation": [[float(value) for value in row] for row in rotation],
        "source_longest_axis": [float(value) for value in axis],
        "aligned_longest_axis": [1.0, 0.0, 0.0],
        "ligand_box": {
            "center": [float(value) for value in (minimum + maximum) / 2.0],
            "size": [float(value) for value in np.maximum(maximum - minimum, 1.0)],
        },
    }


def _apply_coordinate_transform(coordinates: np.ndarray, transform: dict[str, Any]) -> np.ndarray:
    pivot = np.asarray(transform["pivot"], dtype=float)
    rotation = np.asarray(transform["rotation"], dtype=float)
    return (np.asarray(coordinates, dtype=float) - pivot) @ rotation.T + pivot


def transform_pdb_data(pdb_data: str, transform: dict[str, Any]) -> str:
    output: list[str] = []
    transformed_count = 0
    for line in pdb_data.splitlines():
        if line.startswith(("ATOM  ", "HETATM")) and len(line) >= 54:
            try:
                coordinate = np.asarray([[float(line[30:38]), float(line[38:46]), float(line[46:54])]])
            except ValueError:
                output.append(line)
                continue
            x, y, z = _apply_coordinate_transform(coordinate, transform)[0]
            line = f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}"
            transformed_count += 1
        output.append(line)
    if transformed_count == 0:
        raise ValueError("The target contains no PDB coordinate records to align.")
    return "\n".join(output) + "\n"


def transform_axis_aligned_box(
    center: tuple[float, float, float],
    size: tuple[float, float, float],
    transform: dict[str, Any],
) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    box_center = np.asarray(center, dtype=float)
    half_size = np.asarray(size, dtype=float) / 2.0
    corners = np.asarray(
        [
            box_center + np.asarray((x, y, z), dtype=float) * half_size
            for x in (-1.0, 1.0)
            for y in (-1.0, 1.0)
            for z in (-1.0, 1.0)
        ]
    )
    transformed = _apply_coordinate_transform(corners, transform)
    minimum, maximum = transformed.min(axis=0), transformed.max(axis=0)
    return (
        tuple(float(value) for value in (minimum + maximum) / 2.0),
        tuple(float(value) for value in maximum - minimum),
    )


def _docking_axis_alignment(
    target_pdb: Path,
    center: tuple[float, float, float],
    size: tuple[float, float, float],
    *,
    align_longest_axis_to_x: bool,
    alignment_ligand_sdf: Path | None,
) -> tuple[str, tuple[float, float, float], tuple[float, float, float], dict[str, Any] | None]:
    target_data = Path(target_pdb).read_text(encoding="utf-8", errors="replace")
    if not align_longest_axis_to_x:
        return target_data, center, size, None
    if alignment_ligand_sdf is None:
        raise ValueError("Axis alignment requires a coordinate-bearing reference ligand SDF.")
    transform = ligand_longest_axis_transform(Path(alignment_ligand_sdf))
    return (
        transform_pdb_data(target_data, transform),
        *transform_axis_aligned_box(center, size, transform),
        transform,
    )


def bound_ligand_choices(path: Path) -> list[dict[str, Any]]:
    """List non-polymer HETATM residues that can be selected for ligand MD."""
    validate_target_pdb(path)
    water = {"HOH", "WAT", "DOD", "H2O"}
    ions = {"NA", "CL", "K", "CA", "MG", "ZN", "MN", "FE", "CU", "CO", "NI", "CD", "HG"}
    polymer = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
        "ASX", "GLX", "MSE", "SEC", "PYL", "HID", "HIE", "HIP", "CYX", "CYM",
    }
    groups: dict[str, dict[str, Any]] = {}
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("HETATM") or len(line) < 54:
            continue
        resname = line[17:20].strip().upper()
        if not resname or resname in water | ions | polymer:
            continue
        chain = line[21].strip() or "_"
        resseq = line[22:26].strip()
        icode = line[26].strip() or "_"
        key = "|".join((resname, chain, resseq, icode))
        atom_name = line[12:16].strip()
        element = (line[76:78].strip() or atom_name[:1]).upper()
        item = groups.setdefault(key, {
            "key": key, "resname": resname, "chain": chain, "resseq": resseq,
            "icode": icode, "coordinates": {}, "atom_names": set(),
            "heavy_atom_names": set(), "heavy_atom_elements": {},
        })
        if atom_name:
            item["atom_names"].add(atom_name)
            if element != "H":
                item["heavy_atom_names"].add(atom_name)
                item["heavy_atom_elements"].setdefault(atom_name, element)
        try:
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue
        if atom_name:
            item["coordinates"].setdefault(atom_name, xyz)
    result = []
    for item in groups.values():
        atom_count = len(item["atom_names"])
        heavy_count = len(item["heavy_atom_names"])
        if not heavy_count:
            continue
        coordinates = list(item["coordinates"].values())
        center = [
            round(sum(float(xyz[index]) for xyz in coordinates) / max(1, len(coordinates)), 3)
            for index in range(3)
        ]
        result.append({
            key: value for key, value in item.items()
            if key not in {"atom_names", "heavy_atom_names", "heavy_atom_elements", "coordinates"}
        } | {
            "atom_count": atom_count,
            "heavy_atom_count": heavy_count,
            "heavy_atom_elements": list(item["heavy_atom_elements"].values()),
            "heavy_atom_coordinates": {
                atom_name: xyz
                for atom_name, xyz in item["coordinates"].items()
                if atom_name in item["heavy_atom_names"]
            },
            "center": center,
        })
    return sorted(result, key=lambda item: (-item["heavy_atom_count"], item["resname"], item["chain"], item["resseq"]))


def _sdf_heavy_atom_elements(path: Path) -> list[str]:
    """Read the first V2000 or V3000 SDF record without requiring RDKit in Studio."""
    validate_ligand_sdf(path)
    lines = Path(path).read_text(encoding="utf-8", errors="replace").split("$$$$", 1)[0].splitlines()
    if any("V3000" in line for line in lines[:5]):
        in_atoms = False
        elements = []
        for line in lines:
            if "M  V30 BEGIN ATOM" in line:
                in_atoms = True
            elif "M  V30 END ATOM" in line:
                break
            elif in_atoms and "M  V30 " in line:
                fields = line.split()
                if len(fields) >= 4 and fields[3].upper() not in {"H", "D"}:
                    elements.append(fields[3].upper())
        if elements:
            return elements
        raise ValueError("Could not read atom records from the V3000 ligand SDF.")
    if len(lines) < 4:
        raise ValueError("The ligand SDF is missing its atom-count header.")
    try:
        atom_count = int(lines[3][:3])
    except ValueError as exc:
        raise ValueError("The ligand SDF has an invalid V2000 atom-count header.") from exc
    if atom_count < 1 or len(lines) < 4 + atom_count:
        raise ValueError("The ligand SDF atom records are incomplete.")
    elements = []
    for line in lines[4:4 + atom_count]:
        element = line[31:34].strip().upper() if len(line) >= 34 else ""
        if not element:
            fields = line.split()
            element = fields[3].upper() if len(fields) >= 4 else ""
        if element not in {"H", "D"}:
            elements.append(element)
    if not elements:
        raise ValueError("The ligand SDF contains no heavy atoms.")
    return elements


def _sdf_heavy_atom_count(path: Path) -> int:
    return len(_sdf_heavy_atom_elements(path))


def _validate_generation_request(workflow_id: str, parameters: dict[str, Any], *, require_files: bool) -> None:
    count = parameters.get("requested_count", 64)
    batch_size = parameters.get("batch_size", 8)
    seed = parameters.get("seed", 42)
    if type(count) is not int or not 1 <= count <= 10_000:
        raise ValueError("Generated molecule count must be between 1 and 10,000.")
    if type(batch_size) is not int or not 1 <= batch_size <= 512:
        raise ValueError("Generation batch size must be between 1 and 512.")
    if type(seed) is not int or not 1 <= seed <= 2_147_483_646:
        raise ValueError("Generation seed must be a positive 32-bit integer.")
    gpu_device = str(parameters.get("gpu_device", "0")).strip()
    if not gpu_device or gpu_device.lower() in {"none", "cpu", "off"}:
        raise ValueError("These molecule generators require a CUDA GPU device.")
    cpu_cores = parameters.get("cpu_cores", 8)
    if type(cpu_cores) is not int or not 1 <= cpu_cores <= 256:
        raise ValueError("Generation CPU allocation must be between 1 and 256 cores.")
    if workflow_id == "ligand.generation.lddm":
        steps = parameters.get("n_steps", 100)
        noise = parameters.get("sampling_noise", 5.0)
        sampler = parameters.get("sampler", "ForwardEuler")
        if type(steps) is not int or not 1 <= steps <= 2_000:
            raise ValueError("LDDM integration steps must be between 1 and 2,000.")
        if not math.isfinite(float(noise)) or not 0 <= float(noise) <= 20:
            raise ValueError("LDDM sampling noise must be between 0 and 20.")
        if sampler not in {"ForwardEuler", "HeunSampler"}:
            raise ValueError("LDDM sampler must be ForwardEuler or HeunSampler.")
        if require_files:
            validate_target_pdb(Path(parameters["target_pdb"]))
            validate_ligand_sdf(Path(parameters["reference_ligand_sdf"]))
            if not (reference_root() / _LDDM_CHECKPOINT).is_file():
                raise ValueError("The LDDM CD+BB checkpoint is missing from the configured reference directory.")
    elif workflow_id == "ligand.generation.pocketxmol":
        steps = parameters.get("diffusion_steps", 100)
        mean_atoms = parameters.get("ligand_atoms_mean", 28.0)
        std_atoms = parameters.get("ligand_atoms_std", 2.0)
        if type(steps) is not int or not 1 <= steps <= 2_000:
            raise ValueError("PocketXMol diffusion steps must be between 1 and 2,000.")
        if not math.isfinite(float(mean_atoms)) or not 5 <= float(mean_atoms) <= 100:
            raise ValueError("PocketXMol mean atom count must be between 5 and 100.")
        if not math.isfinite(float(std_atoms)) or not 0 <= float(std_atoms) <= 20:
            raise ValueError("PocketXMol atom-count deviation must be between 0 and 20.")
        if require_files:
            validate_target_pdb(Path(parameters["target_pdb"]))
            validate_target_pdb(Path(parameters["pocket_pdb"]))
            for reference in (_POCKETXMOL_CHECKPOINT, _POCKETXMOL_CONFIG):
                if not (reference_root() / reference).is_file():
                    raise ValueError(f"PocketXMol reference file is missing: {reference}")
    else:
        raise ValueError(f"Unknown ligand generation workflow: {workflow_id}.")


def _validate_generator_tile_request(workflow_id: str, parameters: dict[str, Any], *, require_files: bool) -> None:
    engine = workflow_id.removeprefix("ligand.generation.")
    config = _GENERATOR_CONFIGS.get(engine)
    if config is None:
        raise ValueError(f"Unknown ligand generation workflow: {workflow_id}.")
    count = parameters.get("requested_count", 64)
    batch_size = parameters.get("batch_size", 8)
    seed = parameters.get("seed", 42)
    cores = parameters.get("cpu_cores", 8)
    if type(count) is not int or not 1 <= count <= 10_000:
        raise ValueError("Generated molecule count must be between 1 and 10,000.")
    if type(batch_size) is not int or not 1 <= batch_size <= 512:
        raise ValueError("Generation batch size must be between 1 and 512.")
    if type(seed) is not int or not 1 <= seed <= 2_147_483_646:
        raise ValueError("Generation seed must be a positive 32-bit integer.")
    if type(cores) is not int or not 1 <= cores <= 256:
        raise ValueError("Generation CPU allocation must be between 1 and 256 cores.")
    if str(parameters.get("gpu_device", "0")).strip().lower() in {"", "none", "cpu", "off"}:
        raise ValueError(f"{config['name']} requires a CUDA GPU device.")
    for field in config["conditions"]:
        if not parameters.get(field):
            raise ValueError(f"{config['name']} requires {field.replace('_', ' ')}.")
    settings = parameters.get("engine_settings") or {}
    if not isinstance(settings, dict):
        raise ValueError("engine_settings must be an object.")
    if engine == "paopt":
        objectives = settings.get("optimize_properties") or []
        if not isinstance(objectives, list) or not objectives or any(not str(value).strip() for value in objectives):
            raise ValueError("paOPT requires at least one ADMET endpoint to optimize.")
    if engine == "pgmg":
        posp = Path(parameters["pharmacophore_posp"])
        if require_files:
            _validate_pgmg_posp(posp)
    if engine == "omtra" and parameters.get("pharmacophore_json") and require_files:
        _validate_omtra_pharmacophore_json(Path(parameters["pharmacophore_json"]))
    if require_files:
        validate_target_pdb(Path(parameters["target_pdb"]))
        if parameters.get("pocket_pdb"):
            validate_target_pdb(Path(parameters["pocket_pdb"]))
        if parameters.get("reference_ligand_sdf"):
            validate_ligand_sdf(Path(parameters["reference_ligand_sdf"]))
        missing = [relative for relative in config["references"] if not (reference_root() / relative).is_file()]
        if missing:
            raise ValueError(f"{config['name']} reference files are missing: " + ", ".join(missing))


def _validate_pgmg_posp(path: Path) -> None:
    source = Path(path)
    if source.suffix.lower() != ".posp":
        raise ValueError("PGMG pharmacophore input must be a .posp file.")
    try:
        rows = [line.split() for line in source.read_text(encoding="utf-8", errors="replace").splitlines() if line.strip()]
    except OSError as exc:
        raise ValueError(f"Could not read PGMG pharmacophore input: {exc}") from exc
    supported = {"AROM", "HDON", "HACC", "POSC", "HYBL,LHYBL", "HYBL", "LHYBL"}
    if not 1 <= len(rows) <= 8:
        raise ValueError("PGMG supports between one and eight pharmacophore points.")
    for row in rows:
        if len(row) != 4 or row[0] not in supported:
            raise ValueError("Each PGMG point must contain a supported feature type and three coordinates.")
        try:
            coordinates = [float(value) for value in row[1:]]
        except ValueError as exc:
            raise ValueError("PGMG pharmacophore coordinates must be finite numbers.") from exc
        if any(not math.isfinite(value) for value in coordinates):
            raise ValueError("PGMG pharmacophore coordinates must be finite numbers.")


def _validate_omtra_pharmacophore_json(path: Path) -> None:
    source = Path(path)
    if source.suffix.lower() != ".json":
        raise ValueError("OMTRA pharmacophore input must be a Pharmit JSON file.")
    try:
        payload = json.loads(source.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError("Could not read the OMTRA Pharmit pharmacophore JSON.") from exc
    points = payload.get("points") if isinstance(payload, dict) else None
    if not isinstance(points, list) or not points:
        raise ValueError("The OMTRA pharmacophore JSON must contain a non-empty points list.")
    enabled = 0
    for point in points:
        if not isinstance(point, dict):
            raise ValueError("Each OMTRA pharmacophore point must be an object.")
        if point.get("enabled", True) is False:
            continue
        try:
            coordinates = [float(point[axis]) for axis in ("x", "y", "z")]
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError("OMTRA pharmacophore points require numeric x, y, and z coordinates.") from exc
        if any(not math.isfinite(value) for value in coordinates):
            raise ValueError("OMTRA pharmacophore coordinates must be finite numbers.")
        enabled += 1
    if not enabled:
        raise ValueError("The OMTRA pharmacophore JSON has no enabled points.")


def _validate_gromacs_md_request(parameters: dict[str, Any], *, require_files: bool, engine: str = "gromacs") -> None:
    engine_label = "OpenMM" if engine == "openmm" else "GROMACS"
    duration = float(parameters.get("production_ns", 1.0))
    timestep = float(parameters.get("timestep_fs", 2.0))
    temperature = float(parameters.get("temperature", 300.0))
    pressure = float(parameters.get("pressure", 1.0))
    padding = float(parameters.get("padding_nm", 1.0))
    ionic_strength = float(parameters.get("ionic_strength", 0.15))
    if not math.isfinite(duration) or not 0.001 <= duration <= 1000:
        raise ValueError(f"{engine_label} production duration must be between 0.001 and 1,000 ns.")
    if timestep not in {2.0, 4.0}:
        raise ValueError(f"{engine_label} production timestep must be 2 or 4 fs.")
    if timestep == 4.0 and not bool(parameters.get("hydrogen_mass_repartitioning", True)):
        raise ValueError(f"A 4 fs {engine_label} timestep requires hydrogen-mass repartitioning.")
    if not math.isfinite(temperature) or not 200 <= temperature <= 400:
        raise ValueError(f"{engine_label} temperature must be between 200 and 400 K.")
    if not math.isfinite(pressure) or not 0.1 <= pressure <= 20:
        raise ValueError(f"{engine_label} pressure must be between 0.1 and 20 bar.")
    if not math.isfinite(padding) or not 0.6 <= padding <= 3.0:
        raise ValueError(f"{engine_label} solvent padding must be between 0.6 and 3 nm.")
    if not math.isfinite(ionic_strength) or not 0 <= ionic_strength <= 2:
        raise ValueError(f"{engine_label} ionic strength must be between 0 and 2 M.")
    seed = parameters.get("seed", 42)
    cores = parameters.get("cpu_cores", 8)
    if type(seed) is not int or not 1 <= seed <= 2_147_483_646:
        raise ValueError(f"{engine_label} seed must be a positive 32-bit integer.")
    if type(cores) is not int or not 1 <= cores <= 256:
        raise ValueError(f"{engine_label} CPU allocation must be between 1 and 256 cores.")
    allowed_ligand_ff = {"openff-2.2.0", "gaff2"} if engine == "openmm" else {"gaff", "gaff2"}
    allowed_protein_ff = {"amber14-all", "amber19-all"} if engine == "openmm" else {"ff14SB", "ff19SB", "ff15ipq", "ff03.r1"}
    if parameters.get("forcefield_method", "gaff2" if engine == "gromacs" else "openff-2.2.0") not in allowed_ligand_ff:
        raise ValueError(f"The ligand force field is unsupported by the {engine_label} adapter.")
    if parameters.get("protein_forcefield_method", "ff14SB" if engine == "gromacs" else "amber14-all") not in allowed_protein_ff:
        raise ValueError(f"The protein force field is unsupported by the {engine_label} adapter.")
    if parameters.get("water_model", "tip3p") not in {"tip3p", "opc", "spce", "tip4pew"}:
        raise ValueError(f"The selected water model is not supported by the {engine_label} preparation adapter.")
    if bool(parameters.get("use_gpu", True)) and str(parameters.get("gpu_device", "0")).strip().lower() in {"", "none", "cpu", "off"}:
        raise ValueError("Select a CUDA GPU device or turn off GPU acceleration.")
    if require_files:
        complex_path = Path(parameters["complex_pdb"])
        sdf_path = Path(parameters["ligand_sdf"])
        validate_target_pdb(complex_path)
        selected = next((item for item in bound_ligand_choices(complex_path) if item["key"] == str(parameters.get("ligand_key") or "")), None)
        if selected is None:
            raise ValueError("The selected bound ligand residue is not present in the uploaded complex PDB.")
        sdf_heavy_atoms = _sdf_heavy_atom_count(sdf_path)
        if sdf_heavy_atoms != int(selected["heavy_atom_count"]):
            raise ValueError(
                "The ligand SDF does not match the selected PDB residue: "
                f"{sdf_heavy_atoms} SDF heavy atoms versus {selected['heavy_atom_count']} PDB heavy atoms."
            )
        sdf_elements = _sdf_heavy_atom_elements(sdf_path)
        if sdf_elements != list(selected["heavy_atom_elements"]):
            raise ValueError(
                "The bound-ligand SDF must preserve the selected PDB residue's heavy-atom element order. "
                f"PDB: {selected['heavy_atom_elements']}; SDF: {sdf_elements}."
            )
        readiness = openmm_md_runner_readiness() if engine == "openmm" else gromacs_md_runner_readiness()
        if not readiness["ready"]:
            raise ValueError(f"Studio's bundled {engine_label} runtime is unavailable at {readiness['source']}.")


def validate_ligand_tool_request(workflow_id: str, parameters: dict[str, Any], *, require_files: bool = True) -> None:
    """Validate user-visible tile requests before a worker is queued."""
    if workflow_id == "ligand.complex.prepare":
        _validate_complex_preparation_request(parameters, require_files=require_files)
        return
    if workflow_id == "ligand.md.gromacs":
        _validate_gromacs_md_request(parameters, require_files=require_files)
        return
    if workflow_id == "ligand.md.openmm":
        _validate_gromacs_md_request(parameters, require_files=require_files, engine="openmm")
        return
    if workflow_id in {"ligand.generation.lddm", "ligand.generation.pocketxmol"}:
        _validate_generation_request(workflow_id, parameters, require_files=require_files)
        return
    if workflow_id.startswith("ligand.generation."):
        _validate_generator_tile_request(workflow_id, parameters, require_files=require_files)
        return
    compound_limit = parameters.get("compound_limit", 0)
    if type(compound_limit) is not int or not 0 <= compound_limit <= _MAX_COMPOUNDS_PER_LAUNCH:
        raise ValueError(f"compound_limit must be between 0 (all selected) and {_MAX_COMPOUNDS_PER_LAUNCH}.")
    if require_files:
        validate_target_pdb(Path(parameters["target_pdb"]))
        _validated_compound_set(Path(parameters["compound_set_file"]), compound_limit)
    if workflow_id == "ligand.docking.rosetta_ligand":
        protocol = str(parameters.get("protocol", "vsh")).strip().lower()
        reference_mode = str(parameters.get("reference_mode", "pocket_center")).strip().lower()
        if protocol not in {"vsh", "vsx", "convergence"}:
            raise ValueError("RosettaLigand protocol must be VSH, VSX, or convergence.")
        if reference_mode not in {"reference_guided", "pocket_center"}:
            raise ValueError("RosettaLigand placement must be reference-guided or pocket-centered.")
        if reference_mode == "reference_guided" and not parameters.get("reference_ligand_sdf"):
            raise ValueError("Reference-guided RosettaLigand requires reference ligand coordinates in SDF format.")
        if type(parameters.get("cpu_cores", 8)) is not int or not 1 <= parameters.get("cpu_cores", 8) <= 128:
            raise ValueError("RosettaLigand CPU workers must be between 1 and 128.")
        if type(parameters.get("replicates", 1)) is not int or not 1 <= parameters.get("replicates", 1) <= 100:
            raise ValueError("RosettaLigand replicates must be between 1 and 100.")
        ph = float(parameters.get("ph", 7.4))
        padding = float(parameters.get("padding", 4.0))
        if not math.isfinite(ph) or not 0 <= ph <= 14:
            raise ValueError("RosettaLigand preparation pH must be between 0 and 14.")
        if not math.isfinite(padding) or not 1 <= padding <= 30:
            raise ValueError("RosettaLigand search padding must be between 1 and 30 Å.")
        if require_files:
            if not rosettaligand_runner_readiness()["ready"]:
                raise ValueError(f"Studio's bundled RosettaLigand runner is unavailable at {rosettaligand_runner_readiness()['source']}.")
            validate_target_pdb(Path(parameters["target_pdb"]))
            if reference_mode == "reference_guided":
                validate_ligand_sdf(Path(parameters["reference_ligand_sdf"]))
            _validated_compound_set(Path(parameters["compound_set_file"]), compound_limit)
        return
    if workflow_id in {"ligand.docking.vina", "ligand.docking.gnina", "ligand.docking.unidock_pro", "ligand.docking.autodock_gpu"}:
        center = [parameters.get(f"center_{axis}") for axis in "xyz"]
        if any(value is None or not math.isfinite(float(value)) for value in center):
            raise ValueError("Docking box center coordinates must be finite numbers.")
        sizes = [parameters.get(f"size_{axis}", 20.0) for axis in "xyz"]
        if any(not math.isfinite(float(value)) or not 1 <= float(value) <= 60 for value in sizes):
            raise ValueError("Docking box dimensions must be between 1 and 60 Å.")
        if workflow_id == "ligand.docking.unidock_pro":
            if str(parameters.get("search_mode", "balance")).strip().lower() not in {"fast", "balance", "detail"}:
                raise ValueError("Uni-Dock Pro search mode must be fast, balance, or detail.")
            if str(parameters.get("gpu_device", "0")).strip().lower() in {"", "none", "cpu", "off"}:
                raise ValueError("Uni-Dock Pro requires a CUDA GPU device.")
        elif workflow_id == "ligand.docking.autodock_gpu":
            if str(parameters.get("gpu_device", "0")).strip().lower() in {"", "none", "cpu", "off"}:
                raise ValueError("AutoDock-GPU requires a CUDA GPU device.")
            if not 1 <= int(parameters.get("nrun", 20)) <= 500:
                raise ValueError("AutoDock-GPU LGA runs must be between 1 and 500.")
        elif not 1 <= int(parameters.get("exhaustiveness", 8)) <= 1024:
            raise ValueError("Vina exhaustiveness must be between 1 and 1,024.")
        if workflow_id == "ligand.docking.gnina" and str(parameters.get("cnn_scoring", "rescore")).strip().lower() not in {"none", "rescore", "refinement", "all"}:
            raise ValueError("GNINA CNN scoring must be none, rescore, refinement, or all.")
        use_scrub = parameters.get("use_scrub", True)
        skip_tautomer = parameters.get("scrub_skip_tautomer", True)
        scrub_ph = float(parameters.get("scrub_ph", 7.4))
        if type(use_scrub) is not bool or type(skip_tautomer) is not bool:
            raise ValueError("Scrub and tautomer-enumeration settings must be booleans.")
        if not math.isfinite(scrub_ph) or not 0 <= scrub_ph <= 14:
            raise ValueError("Scrub pH must be between 0 and 14.")
        if not 1 <= int(parameters.get("num_modes", 9)) <= 100:
            raise ValueError("Vina pose count must be between 1 and 100.")
        if not 1 <= int(parameters.get("seed", 42)) <= 2_147_483_646:
            raise ValueError("Vina seed must be a positive 32-bit integer.")
        align = parameters.get("align_longest_axis_to_x", False)
        if type(align) is not bool:
            raise ValueError("align_longest_axis_to_x must be a boolean.")
        if align:
            ligand_path = parameters.get("alignment_ligand_sdf")
            if not ligand_path:
                raise ValueError("Axis alignment requires a coordinate-bearing reference ligand SDF.")
            if require_files:
                validate_ligand_sdf(Path(ligand_path))
                _sdf_heavy_atom_coordinates(Path(ligand_path))
    elif workflow_id == "ligand.affinity.nesso":
        if not 1 <= int(parameters.get("recycling_steps", 5)) <= 64:
            raise ValueError("Nesso recycling steps must be between 1 and 64.")
        if not 1 <= int(parameters.get("seed", 42)) <= 2_147_483_646:
            raise ValueError("Nesso seed must be a positive 32-bit integer.")
        if require_files and not nesso_readiness()["ready"]:
            raise ValueError("Nesso references are incomplete. Install its checkpoint, CCD dictionary, and ESM2 cache under the configured reference directory.")
    else:
        raise ValueError(f"Unknown ligand tool workflow: {workflow_id}.")


def _pdb_polymer_chain_ranges(pdb_data: str) -> dict[str, dict[str, Any]]:
    residues: dict[str, set[int]] = {}
    amino_acids = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
        "ASX", "GLX", "MSE", "SEC", "PYL", "HID", "HIE", "HIP", "CYX", "CYM",
    }
    for line in pdb_data.splitlines():
        if not line.startswith("ATOM  ") or len(line) < 27 or line[17:20].strip().upper() not in amino_acids:
            continue
        try:
            resseq = int(line[22:26])
        except ValueError:
            continue
        residues.setdefault(line[21].strip() or "_", set()).add(resseq)
    return {
        chain: {"start": min(values), "end": max(values), "residues": sorted(values)}
        for chain, values in sorted(residues.items())
        if values
    }


def _validate_complex_preparation_request(parameters: dict[str, Any], *, require_files: bool) -> None:
    max_gap = parameters.get("max_internal_gap", 15)
    if type(max_gap) is not int or not 0 <= max_gap <= 100:
        raise ValueError("Maximum internal gap must be an integer from 0 to 100 residues.")
    if type(parameters.get("repair", True)) is not bool:
        raise ValueError("repair must be a boolean.")
    extension = "".join(str(parameters.get("extension_sequence") or "").split()).upper()
    if len(extension) > 30 or any(residue not in "ACDEFGHIKLMNPQRSTVWY" for residue in extension):
        raise ValueError("C-terminal extension must use up to 30 canonical amino-acid letters.")
    model_count = parameters.get("extension_models", 10)
    if type(model_count) is not int or not 1 <= model_count <= 50:
        raise ValueError("Terminal-extension models must be between 1 and 50.")
    if require_files:
        source = Path(parameters["source_pdb"])
        validate_target_pdb(source)
        if not ligand_complex_preparation_readiness()["ready"]:
            raise ValueError(
                f"Studio's bundled complex-preparation runtime is unavailable at {ligand_complex_preparation_readiness()['source']}."
            )
        pdb_data = source.read_text(encoding="utf-8", errors="replace")
        chains = _pdb_polymer_chain_ranges(pdb_data)
        if not chains:
            raise ValueError("The uploaded complex must contain at least one protein chain.")
        if extension and str(parameters.get("extension_chain") or "") not in chains:
            raise ValueError("Select a protein chain present in the uploaded complex for C-terminal extension.")
        raw_ranges = parameters.get("trim_ranges_json") or ""
        if raw_ranges:
            try:
                trim_ranges = json.loads(str(raw_ranges))
            except json.JSONDecodeError as exc:
                raise ValueError("Trim ranges must be valid JSON.") from exc
            if not isinstance(trim_ranges, dict) or not trim_ranges:
                raise ValueError("Select at least one protein chain range to trim.")
            for chain, bounds in trim_ranges.items():
                if chain not in chains or not isinstance(bounds, (list, tuple)) or len(bounds) != 2:
                    raise ValueError(f"Invalid trim selection for chain {chain}.")
                start, end = int(bounds[0]), int(bounds[1])
                chain_data = chains[chain]
                allowed_end = chain_data["end"] + (len(extension) if chain == str(parameters.get("extension_chain") or "") else 0)
                if start > end or start < chain_data["start"] or end > allowed_end:
                    raise ValueError(f"Trim range for chain {chain} is outside its observed residue range.")


def ligand_complex_preparation_readiness() -> dict[str, Any]:
    root = _studio_ligand_runtime_root()
    required = (
        root / "mn_ligand" / "workflows" / "protein_preparation.py",
        root / "mn_ligand" / "workflows" / "terminal_repair.py",
        root / "mn_ligand" / "workflows" / "target_trimming.py",
        root / "mn_ligand" / "workflows" / "modeller_terminal_extension_runner.py",
    )
    return {"ready": all(path.is_file() for path in required), "source": str(root)}


def run_ligand_complex_prepare(
    source_pdb: Path,
    repair: bool = True,
    map_modified_residues: bool = True,
    skip_terminal_missing_residues: bool = True,
    max_internal_gap: int = 15,
    extension_chain: str = "",
    extension_sequence: str = "",
    extension_models: int = 10,
    trim_ranges_json: str = "",
) -> Path:
    """Import and repair a ligand complex using Studio's bundled ligand runtime."""
    source_pdb = Path(source_pdb).expanduser().resolve()
    parameters = {
        "source_pdb": source_pdb,
        "repair": repair,
        "map_modified_residues": map_modified_residues,
        "skip_terminal_missing_residues": skip_terminal_missing_residues,
        "max_internal_gap": max_internal_gap,
        "extension_chain": extension_chain,
        "extension_sequence": extension_sequence,
        "extension_models": extension_models,
        "trim_ranges_json": trim_ranges_json,
    }
    validate_ligand_tool_request("ligand.complex.prepare", parameters)
    job = create_job(
        "ligand-preparation",
        "ligand_complex_preparation",
        "studio_ligand_complex_preparation",
        {"source_pdb": str(source_pdb)},
        {key: value for key, value in parameters.items() if key != "source_pdb"},
    )
    input_dir = job.run_dir / "artifacts" / "input"
    prepared_dir = job.run_dir / "artifacts" / "prepared"
    input_dir.mkdir(parents=True, exist_ok=True)
    prepared_dir.mkdir(parents=True, exist_ok=True)
    staged_source = input_dir / "source_complex.pdb"
    shutil.copy2(source_pdb, staged_source)
    params_path = input_dir / "preparation_parameters.json"
    write_json(params_path, {key: value for key, value in parameters.items() if key != "source_pdb"})
    helper = Path(__file__).with_name("ligand_complex_prepare_runner.py")
    source_root = _studio_ligand_runtime_root()
    command = [sys.executable, str(helper), str(staged_source), str(prepared_dir), str(params_path)]
    write_json(
        job.run_dir / "command.json",
        {"mode": "host", "command": command, "studio_ligand_runtime": str(source_root)},
    )
    update_status(job.run_dir, "running", current_phase="Importing and preparing protein–ligand complex")
    environment = os.environ.copy()
    existing_pythonpath = environment.get("PYTHONPATH", "")
    environment["PYTHONPATH"] = str(source_root) + (os.pathsep + existing_pythonpath if existing_pythonpath else "")
    environment["MN_LIGAND_RUN_DIR"] = str(job.run_dir / "artifacts" / "mn_ligand_runs")
    try:
        process = subprocess.run(command, capture_output=True, text=True, check=False, env=environment)
        (job.run_dir / "stdout.log").write_text(process.stdout or "", encoding="utf-8")
        (job.run_dir / "stderr.log").write_text(process.stderr or "", encoding="utf-8")
        complex_path = prepared_dir / "prepared_complex.pdb"
        target_path = prepared_dir / "prepared_target.pdb"
        report_path = prepared_dir / "repair_report.json"
        success = process.returncode == 0 and complex_path.is_file() and target_path.is_file() and report_path.is_file()
        report = read_json(report_path) if report_path.is_file() else {}
        artifact_outputs = [
            {"name": "prepared_complex", "type": "protein_ligand_complex/pdb/v1", "path": "artifacts/prepared/prepared_complex.pdb"},
            {"name": "prepared_target", "type": "protein_structure/pdb/v1", "path": "artifacts/prepared/prepared_target.pdb"},
            {"name": "repair_report", "type": "complex_preparation_report/v1", "path": "artifacts/prepared/repair_report.json"},
        ] if success else []
        finish_job(job.run_dir, success, {
            "outputs": {"artifacts": artifact_outputs},
            "metrics": {
                "repair_enabled": bool(repair),
                "terminal_extension_enabled": bool(str(extension_sequence).strip()),
                "trimming_enabled": bool(str(trim_ranges_json).strip()),
                "source": "Studio-bundled complex import, protein repair, optional MODELLER C-terminal extension, and chain-aware trimming",
            },
            "downstream_artifacts": {
                "prepared_complex": "artifacts/prepared/prepared_complex.pdb",
                "prepared_target": "artifacts/prepared/prepared_target.pdb",
                "repair_report": "artifacts/prepared/repair_report.json",
            } if success else {},
            "error": "\n".join((process.stderr or process.stdout or "Complex preparation failed").splitlines()[-40:]) if not success else "",
            "source_import": report.get("source_import", {}),
        })
    except Exception as exc:
        (job.run_dir / "stderr.log").write_text(f"{type(exc).__name__}: {exc}\n", encoding="utf-8")
        finish_job(job.run_dir, False, {"outputs": {}, "metrics": {}, "error": f"{type(exc).__name__}: {exc}"})
        raise
    return job.run_dir


def _studio_ligand_runtime_root() -> Path:
    """Return the ligand runtime vendored inside the Studio source/package."""
    return (Path(__file__).resolve().parents[1] / "vendor" / "ligand_runtime").resolve()


def gromacs_md_runner_readiness() -> dict[str, Any]:
    """Check Studio's bundled source files required by the GROMACS container."""
    root = _studio_ligand_runtime_root()
    required = (
        root / "mn_ligand" / "workflows" / "gromacs_md.py",
        root / "mn_ligand" / "workflows" / "bound_ligand_md.py",
    )
    return {
        "ready": all(path.is_file() for path in required),
        "source": str(root),
        "image": _GROMACS_MD_IMAGE,
        "missing": [str(path) for path in required if not path.is_file()],
    }


def openmm_md_runner_readiness() -> dict[str, Any]:
    """Check Studio's bundled source files mounted into the OpenMM image."""
    root = _studio_ligand_runtime_root()
    required = (
        root / "mn_ligand" / "workflows" / "bound_ligand_md.py",
        root / "mn_ligand" / "ligandx" / "services" / "md" / "service.py",
    )
    return {
        "ready": all(path.is_file() for path in required),
        "source": str(root),
        "image": _OPENMM_MD_IMAGE,
        "missing": [str(path) for path in required if not path.is_file()],
    }


def rosettaligand_runner_readiness() -> dict[str, Any]:
    """Check Studio's bundled RosettaLigand runner source."""
    source = _studio_ligand_runtime_root() / "mn_ligand" / "workflows" / "openvs.py"
    return {"ready": source.is_file(), "source": str(source), "image": _ROSETTALIGAND_IMAGE}


def _openvs_runner_constants() -> dict[str, str]:
    source = Path(rosettaligand_runner_readiness()["source"])
    if not source.is_file():
        raise FileNotFoundError(f"Studio's bundled RosettaLigand runner source was not found: {source}")
    tree = ast.parse(source.read_text(encoding="utf-8"), filename=str(source))
    wanted = {"_RUNNER_SCRIPT", "_CENTER_PDB_SCRIPT", "_ROSETTA_PARAMS_VALIDATOR_SCRIPT"}
    constants: dict[str, str] = {}
    for node in tree.body:
        if isinstance(node, ast.Assign):
            names = {target.id for target in node.targets if isinstance(target, ast.Name)}
            match = names.intersection(wanted)
            if match:
                value = ast.literal_eval(node.value)
                for name in match:
                    constants[name] = str(value)
    missing = wanted - set(constants)
    if missing:
        raise ValueError(f"Studio's bundled RosettaLigand runner is missing code: {', '.join(sorted(missing))}")
    return constants


def nesso_reference_paths() -> tuple[Path, Path, Path]:
    root = reference_root().expanduser().resolve()
    checkpoint = Path(os.getenv("MN_NESSO_CHECKPOINT_DIR", str(root / "nesso" / "v1.0.0"))).expanduser()
    ccd = Path(os.getenv("MN_NESSO_CCD_PATH", str(root / "nesso" / "ccd.pkl"))).expanduser()
    configured_cache = os.getenv("MN_NESSO_ESM_CACHE_DIR", "").strip()
    if configured_cache:
        esm_cache = Path(configured_cache).expanduser()
    else:
        portable_cache = root / "nesso" / "huggingface"
        legacy_cache = root / "nesso" / ".cache" / "huggingface"
        shared_cache = root / "proteina-complexa" / "hf-cache"
        candidates = (portable_cache, legacy_cache, shared_cache)
        esm_cache = next((candidate for candidate in candidates if _nesso_esm_cache_ready(candidate)), None)
        if esm_cache is None:
            esm_cache = next((candidate for candidate in candidates if candidate.exists()), portable_cache)
    return checkpoint.resolve(), ccd.resolve(), esm_cache.resolve()


def _nesso_esm_cache_ready(cache: Path) -> bool:
    model_cache = cache / "models--facebook--esm2_t33_650M_UR50D" / "snapshots"
    return any(
        (snapshot / "model.safetensors").is_file()
        and (snapshot / "config.json").is_file()
        and (snapshot / "vocab.txt").is_file()
        for snapshot in model_cache.glob("*")
    )


def nesso_readiness() -> dict[str, Any]:
    checkpoint, ccd, cache = nesso_reference_paths()
    esm_ready = _nesso_esm_cache_ready(cache)
    checkpoint_ready = (checkpoint / "model.safetensors").is_file() and (checkpoint / "hparams.json").is_file()
    return {
        "checkpoint_dir": str(checkpoint),
        "ccd_path": str(ccd),
        "esm_cache_dir": str(cache),
        "checkpoint_ready": checkpoint_ready,
        "ccd_ready": ccd.is_file(),
        "esm_ready": esm_ready,
        "ready": checkpoint_ready and ccd.is_file() and esm_ready,
    }


def _pdb_protein_sequences(path: Path) -> list[tuple[str, str]]:
    structure = PDBParser(QUIET=True).get_structure("ligand-target", str(path))
    model = next(structure.get_models(), None)
    if model is None:
        raise ValueError("The target structure contains no coordinate model.")
    chains: list[tuple[str, str]] = []
    for chain in model:
        residues = []
        for residue in chain:
            if residue.id[0] != " ":
                continue
            letter = seq1(residue.resname, custom_map={"MSE": "M"}, undef_code="X")
            if letter != "X":
                residues.append(letter)
        sequence = "".join(residues)
        if sequence:
            chain_id = str(chain.id or "_")
            chains.append((chain_id, sequence))
    if not chains:
        raise ValueError("The target PDB contains no standard protein residues for Nesso.")
    return chains


def _nesso_yaml(protein_chains: list[tuple[str, str]], smiles: str) -> str:
    if not smiles or "\n" in smiles or "\r" in smiles:
        raise ValueError("Compound SMILES must be one non-empty line.")
    lines = ["sequences:"]
    for index, (_source_chain, sequence) in enumerate(protein_chains):
        chain_id = chr(ord("A") + index)
        lines.extend([
            "  - protein:",
            f"      id: {chain_id}",
            f"      sequence: {json.dumps(sequence)}",
            "      msa: empty",
        ])
    lines.extend([
        "  - ligand:",
        "      id: L",
        f"      smiles: {json.dumps(smiles)}",
        "properties:",
        "  - affinity:",
        "      binder: L",
    ])
    return "\n".join(lines) + "\n"


def _write_compound_inputs(run_dir: Path, payload: dict[str, Any], compounds: list[dict[str, Any]]) -> list[dict[str, Any]]:
    input_root = run_dir / "artifacts" / "input"
    input_root.mkdir(parents=True, exist_ok=True)
    staged_path = input_root / "compound_set.json"
    staged_payload = ligand_compounds.write_compound_set(
        staged_path,
        name=str(payload.get("name") or "Compound set"),
        compounds=compounds,
        owner_id=str(payload.get("owner_id") or ""),
        project_id=str(payload.get("project_id") or ""),
        created_at=str(payload.get("created_at") or ""),
    )
    staged_payload["source_compound_set_id"] = str(payload.get("compound_set_id") or "")
    staged_payload["source_content_sha256"] = str(payload.get("content_sha256") or "")
    staged_path.write_text(json.dumps(staged_payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return compounds


def _finite_float(payload: dict[str, Any], key: str) -> float | None:
    try:
        value = float(payload[key])
    except (KeyError, TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def _write_csv(path: Path, rows: list[dict[str, Any]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _write_job_contract(
    job_dir: Path,
    workflow: str,
    input_type: str,
    output_type: str,
    compound_set: dict[str, Any],
    compounds: list[dict[str, Any]],
    *,
    target_artifact: str,
    outputs: list[dict[str, Any]],
    extra_inputs: list[dict[str, Any]] | None = None,
) -> None:
    source_id = str(compound_set.get("compound_set_id") or "")
    source_digest = str(compound_set.get("content_sha256") or "")
    metadata_path = Path(job_dir) / "metadata.json"
    metadata = read_json(metadata_path)
    metadata.update({
        "workflow": workflow,
        "engine": {
            "ligand.docking.vina": "AutoDock Vina",
            "ligand.docking.gnina": "GNINA",
            "ligand.docking.unidock_pro": "Uni-Dock Pro",
            "ligand.docking.rosetta_ligand": "RosettaLigand",
            "ligand.docking.autodock_gpu": "AutoDock-GPU",
            "ligand.affinity.nesso": "Nesso-1",
        }.get(workflow, workflow),
        "result_type": output_type,
        "artifact_contract_version": 1,
        "compound_set_name": str(compound_set.get("name") or "Compound set"),
        "compound_set_id": source_id,
        "compound_set_sha256": source_digest,
        "compound_count": len(compounds),
    })
    write_json(metadata_path, metadata)
    write_json(Path(job_dir) / "artifact_contract.json", {
        "schema_version": 1,
        "workflow": workflow,
        "input_type": input_type,
        "output_type": output_type,
        "inputs": [
            {"name": "target_structure", "type": "protein_structure/pdb/v1", "path": target_artifact},
            {
                "name": "compound_set", "type": "compound_set/v1",
                "path": "artifacts/input/compound_set.json",
                "source_record_id": source_id,
                "source_content_sha256": source_digest,
                "consumed_count": len(compounds),
                "consumed_content_sha256": ligand_compounds.compound_content_sha256(compounds),
            },
            *(extra_inputs or []),
        ],
        "outputs": outputs,
        "structure_output": workflow in {"ligand.docking.vina", "ligand.docking.gnina", "ligand.docking.unidock_pro", "ligand.docking.autodock_gpu"},
        "compound_set_name": str(compound_set.get("name") or "Compound set"),
        "compound_set_id": source_id,
        "compound_set_sha256": source_digest,
        "compound_count": len(compounds),
    })


def _run_docking(
    *,
    engine: str,
    target_pdb: Path,
    compound_set_file: Path,
    center_x: float,
    center_y: float,
    center_z: float,
    size_x: float,
    size_y: float,
    size_z: float,
    exhaustiveness: int,
    num_modes: int,
    seed: int,
    compound_limit: int,
    cpu_cores: int,
    gpu_device: str,
    compound_workers: int | None = None,
    align_longest_axis_to_x: bool = False,
    alignment_ligand_sdf: Path | None = None,
    job_name: str | None = None,
    cnn_scoring: str = "rescore",
    use_scrub: bool = True,
    scrub_ph: float = 7.4,
    scrub_skip_tautomer: bool = True,
) -> Path:
    if engine not in {"vina", "gnina"}:
        raise ValueError(f"Unsupported small-molecule docking engine: {engine}.")
    validate_target_pdb(target_pdb)
    if any(not math.isfinite(value) for value in (center_x, center_y, center_z)):
        raise ValueError("Docking box center coordinates must be finite numbers.")
    if any(not math.isfinite(value) or not 1 <= value <= 60 for value in (size_x, size_y, size_z)):
        raise ValueError("Docking box dimensions must be between 1 and 60 Å.")
    if not 1 <= exhaustiveness <= 1024 or not 1 <= num_modes <= 100 or not 1 <= seed <= 2_147_483_646:
        raise ValueError("Exhaustiveness, pose count, or seed is outside the supported range.")
    cnn_scoring = str(cnn_scoring).strip().lower()
    if engine == "gnina" and cnn_scoring not in {"none", "rescore", "refinement", "all"}:
        raise ValueError("GNINA CNN scoring must be none, rescore, refinement, or all.")
    if type(cpu_cores) is not int or not 1 <= cpu_cores <= 256:
        raise ValueError("cpu_cores must be between 1 and 256.")
    if type(use_scrub) is not bool or type(scrub_skip_tautomer) is not bool:
        raise ValueError("Scrub and tautomer-enumeration settings must be booleans.")
    scrub_ph = float(scrub_ph)
    if not math.isfinite(scrub_ph) or not 0 <= scrub_ph <= 14:
        raise ValueError("Scrub pH must be between 0 and 14.")
    original_center = [float(center_x), float(center_y), float(center_z)]
    original_box_size = [float(size_x), float(size_y), float(size_z)]
    alignment_ligand_sdf = Path(alignment_ligand_sdf) if alignment_ligand_sdf else None
    target_data, aligned_center, aligned_size, coordinate_transform = _docking_axis_alignment(
        target_pdb,
        (float(center_x), float(center_y), float(center_z)),
        (float(size_x), float(size_y), float(size_z)),
        align_longest_axis_to_x=align_longest_axis_to_x,
        alignment_ligand_sdf=alignment_ligand_sdf,
    )
    center_x, center_y, center_z = aligned_center
    size_x, size_y, size_z = aligned_size
    payload, compounds = _validated_compound_set(compound_set_file, compound_limit)
    if engine == "gnina":
        if compound_workers is None:
            compound_workers = 2
        if type(compound_workers) is not int or not 1 <= compound_workers <= 4:
            raise ValueError("GNINA compound_workers must be an integer between 1 and 4.")
        compound_workers = min(compound_workers, cpu_cores, len(compounds))
        cpu_threads_per_compound = max(1, cpu_cores // compound_workers)
    else:
        compound_workers = 1
        cpu_threads_per_compound = cpu_cores
    params = {
        "engine": engine,
        "center": [float(center_x), float(center_y), float(center_z)],
        "box_size": [float(size_x), float(size_y), float(size_z)],
        "original_center": original_center,
        "original_box_size": original_box_size,
        "exhaustiveness": int(exhaustiveness), "num_modes": int(num_modes),
        "seed": int(seed), "compound_limit": int(compound_limit), "cpu_cores": cpu_cores,
        **({
            "compound_workers": min(cpu_cores, len(compounds)),
            "cpu_threads_per_compound": 1,
        } if engine == "vina" else {}),
        **({
            "compound_workers": compound_workers,
            "cpu_threads_per_compound": cpu_threads_per_compound,
        } if engine == "gnina" else {}),
        "compound_set_name": payload.get("name", "Compound set"),
        "align_longest_axis_to_x": bool(coordinate_transform),
        "coordinate_transform": coordinate_transform,
        "use_scrub": use_scrub,
        "scrub_ph": scrub_ph,
        "scrub_skip_tautomer": scrub_skip_tautomer,
    }
    if job_name is not None:
        normalized_job_name = str(job_name).strip()[:120]
        if normalized_job_name:
            params["job_name"] = normalized_job_name
    if engine == "gnina":
        params["cnn_scoring"] = cnn_scoring
    workflow = f"ligand.docking.{engine}"
    display_engine = "GNINA" if engine == "gnina" else "AutoDock Vina"
    job = create_job(
        "ligand-docking", "ligand_docking", engine,
        {
            "target_pdb": str(target_pdb), "compound_set_file": str(compound_set_file),
            **({"alignment_ligand_sdf": str(alignment_ligand_sdf)} if coordinate_transform and alignment_ligand_sdf else {}),
        }, params,
    )
    contract_outputs = [
        {"name": "docking_scores", "type": "ligand_docking_scores/v1", "path": "artifacts/docking_scores.csv"},
        {"name": "poses", "type": "ligand_pose_collection/v1", "path": "artifacts/poses"},
    ]
    if coordinate_transform:
        contract_outputs.append({"name": "coordinate_transform", "type": "coordinate_transform/v1", "path": "artifacts/input/alignment_transform.json"})
    _write_job_contract(
        job.run_dir, workflow, "compound_set/v1", "ligand_docking_result/v1",
        payload, compounds, target_artifact="artifacts/input/receptor.pdb",
        outputs=contract_outputs,
        extra_inputs=(
            [{"name": "axis_alignment_ligand", "type": "ligand_coordinates/sdf/v1", "path": "artifacts/input/alignment_ligand.sdf"}]
            if coordinate_transform
            else []
        ),
    )
    artifacts = job.run_dir / "artifacts"
    input_dir = artifacts / "input"
    input_dir.mkdir(parents=True, exist_ok=True)
    receptor = input_dir / "receptor.pdb"
    receptor.write_text(target_data, encoding="utf-8")
    if coordinate_transform and alignment_ligand_sdf is not None:
        shutil.copy2(alignment_ligand_sdf, input_dir / "alignment_ligand.sdf")
        write_json(input_dir / "alignment_transform.json", coordinate_transform)
    staged_compounds = _write_compound_inputs(job.run_dir, payload, compounds)
    (input_dir / "smiles_to_3d.py").write_text(_RDKIT_3D_SCRIPT, encoding="utf-8")
    with (input_dir / "compounds.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["compound_id", "smiles"])
        writer.writerows((item["compound_id"], item["smiles"]) for item in staged_compounds)
    manifest = load_manifest(engine)
    if engine == "vina":
        dock_command = (
            f"vina --receptor artifacts/prepared/receptor.pdbqt --ligand \"artifacts/ligands/${{compound_id}}.pdbqt\" "
            f"--center_x {center_x:.3f} --center_y {center_y:.3f} --center_z {center_z:.3f} "
            f"--size_x {size_x:.3f} --size_y {size_y:.3f} --size_z {size_z:.3f} "
            f"--exhaustiveness {int(exhaustiveness)} --num_modes {int(num_modes)} --seed {int(seed)} --cpu 1 "
            "--out \"artifacts/poses/${compound_id}.pdbqt\" "
        )
        dock_command += "> \"artifacts/logs/${compound_id}.vina.log\" 2>&1"
        score_parse = (
            "score=$(awk '$1 ~ /^[0-9]+$/ && $2 ~ /^-?[0-9]+([.][0-9]+)?$/ {print $2; exit}' "
            "\"artifacts/logs/${compound_id}.vina.log\"); cnn_score=''; cnn_affinity=''"
        )
        failure_message = "Vina returned an error"
        selected_gpu = "none"
    else:
        dock_command = (
            f"gnina --receptor artifacts/prepared/receptor.pdbqt --ligand \"artifacts/ligands/${{compound_id}}.pdbqt\" "
            f"--center_x {center_x:.3f} --center_y {center_y:.3f} --center_z {center_z:.3f} "
            f"--size_x {size_x:.3f} --size_y {size_y:.3f} --size_z {size_z:.3f} "
            f"--exhaustiveness {int(exhaustiveness)} --num_modes {int(num_modes)} --seed {int(seed)} --cpu {cpu_threads_per_compound} "
            f"--cnn_scoring {cnn_scoring} --out \"artifacts/poses/${{compound_id}}.pdbqt\" "
            "--log \"artifacts/logs/${compound_id}.gnina.log\""
        )
        score_parse = (
            "score=$(awk '/^REMARK minimizedAffinity / {print $3; exit}' \"artifacts/poses/${compound_id}.pdbqt\"); "
            "if [[ -z \"$score\" ]]; then score=$(awk '/^REMARK VINA RESULT:/ {print $4; exit}' \"artifacts/poses/${compound_id}.pdbqt\"); fi; "
            "cnn_score=$(awk '/^REMARK CNNscore / {print $3; exit}' \"artifacts/poses/${compound_id}.pdbqt\"); "
            "cnn_affinity=$(awk '/^REMARK CNNaffinity / {print $3; exit}' \"artifacts/poses/${compound_id}.pdbqt\")"
        )
        failure_message = "GNINA returned an error"
        selected_gpu = gpu_device
    scrub_step = (
        "if scrub.py \"artifacts/ligands/${compound_id}_source.sdf\" "
        "-o \"artifacts/ligands/${compound_id}.sdf\" "
        f"--ph {scrub_ph:.2f}"
        + (" --skip_tautomer" if scrub_skip_tautomer else "")
        + " >>\"artifacts/logs/${compound_id}.prepare.log\" 2>&1 "
        + "&& [[ -s \"artifacts/ligands/${compound_id}.sdf\" ]]; then :; else "
        + "cp \"artifacts/ligands/${compound_id}_source.sdf\" \"artifacts/ligands/${compound_id}.sdf\"; "
        + "echo 'Scrub failed; using the RDKit SDF for Meeko preparation.' "
        + ">> \"artifacts/logs/${compound_id}.prepare.log\"; fi; "
        if use_scrub else
        "cp \"artifacts/ligands/${compound_id}_source.sdf\" \"artifacts/ligands/${compound_id}.sdf\"; "
    )
    common_setup = (
        "set -euo pipefail; cd /work; "
        "mkdir -p artifacts/prepared artifacts/ligands artifacts/poses artifacts/logs; "
        "if ! mk_prepare_receptor.py --read_pdb artifacts/input/receptor.pdb -o artifacts/prepared/receptor -p; then "
        "awk 'substr($0,77,2) !~ /^[[:space:]]*H[[:space:]]*$/' artifacts/input/receptor.pdb > artifacts/input/receptor-noh.pdb; "
        "mk_prepare_receptor.py -i artifacts/input/receptor-noh.pdb -o artifacts/prepared/receptor -p; fi; "
        "printf 'compound_id\\tstatus\\tbest_score_kcal_mol\\tcnn_score\\tcnn_affinity\\tpose_path\\terror\\n' > artifacts/docking_scores.csv; "
    )
    preparation_loop = (
        "tail -n +2 artifacts/input/compounds.tsv | while IFS=$'\\t' read -r compound_id smiles; do "
        "[[ -n \"$compound_id\" && -n \"$smiles\" ]] || continue; "
        "if ! python artifacts/input/smiles_to_3d.py --smiles \"$smiles\" --name \"$compound_id\" "
        "--output \"artifacts/ligands/${compound_id}_source.sdf\" >\"artifacts/logs/${compound_id}.prepare.log\" 2>&1; then "
        "printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \"$compound_id\" failed '' '' '' '' 'RDKit 3D preparation failed' > \"artifacts/logs/${compound_id}.result.tsv\"; continue; fi; "
        + scrub_step
        + "if ! mk_prepare_ligand.py -i \"artifacts/ligands/${compound_id}.sdf\" "
        "-o \"artifacts/ligands/${compound_id}.pdbqt\" >>\"artifacts/logs/${compound_id}.prepare.log\" 2>&1; then "
        "printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \"$compound_id\" failed '' '' '' '' 'Meeko ligand preparation failed' > \"artifacts/logs/${compound_id}.result.tsv\"; continue; fi; "
    )
    if engine == "vina":
        vina_workers = min(cpu_cores, len(compounds))
        script = (
            common_setup
            + " : > artifacts/input/ready_compounds.txt; "
            + preparation_loop
            + 'printf "%s\\n" "$compound_id" >> artifacts/input/ready_compounds.txt; done; '
            + f"DOCKING_COMPOUND_WORKERS={vina_workers}; "
            + "dock_vina_compound() { compound_id=\"$1\"; if "
            + dock_command
            + "; then " + score_parse + "; "
            + "printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \"$compound_id\" completed \"$score\" \"$cnn_score\" \"$cnn_affinity\" \"artifacts/poses/${compound_id}.pdbqt\" '' > \"artifacts/logs/${compound_id}.result.tsv.tmp\"; "
            + "else printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \"$compound_id\" failed '' '' '' '' 'Vina returned an error' > \"artifacts/logs/${compound_id}.result.tsv.tmp\"; fi; "
            + "mv \"artifacts/logs/${compound_id}.result.tsv.tmp\" \"artifacts/logs/${compound_id}.result.tsv\"; "
            + "}; export -f dock_vina_compound; "
            + "xargs -r -P \"$DOCKING_COMPOUND_WORKERS\" -n 1 bash -c 'dock_vina_compound \"$1\"' _ "
            + "< artifacts/input/ready_compounds.txt; "
            + "while IFS=$'\\t' read -r compound_id smiles; do [[ -n \"$compound_id\" ]] || continue; "
            + "if [[ -s \"artifacts/logs/${compound_id}.result.tsv\" ]]; then cat \"artifacts/logs/${compound_id}.result.tsv\" >> artifacts/docking_scores.csv; "
            + "else printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \"$compound_id\" failed '' '' '' '' 'Compound task produced no result' >> artifacts/docking_scores.csv; fi; "
            + "done < <(tail -n +2 artifacts/input/compounds.tsv)"
        )
    elif engine == "gnina":
        script = (
            common_setup
            + " : > artifacts/input/ready_compounds.txt; "
            + preparation_loop
            + 'printf "%s\\n" "$compound_id" >> artifacts/input/ready_compounds.txt; done; '
            + f"DOCKING_COMPOUND_WORKERS={compound_workers}; "
            + "dock_gnina_compound() { compound_id=\"$1\"; if "
            + dock_command
            + "; then " + score_parse + "; "
            + "printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \"$compound_id\" completed \"$score\" \"$cnn_score\" \"$cnn_affinity\" \"artifacts/poses/${compound_id}.pdbqt\" '' > \"artifacts/logs/${compound_id}.result.tsv.tmp\"; "
            + "else printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \"$compound_id\" failed '' '' '' '' 'GNINA returned an error' > \"artifacts/logs/${compound_id}.result.tsv.tmp\"; fi; "
            + "mv \"artifacts/logs/${compound_id}.result.tsv.tmp\" \"artifacts/logs/${compound_id}.result.tsv\"; "
            + "}; export -f dock_gnina_compound; "
            + "xargs -r -P \"$DOCKING_COMPOUND_WORKERS\" -n 1 bash -c 'dock_gnina_compound \"$1\"' _ "
            + "< artifacts/input/ready_compounds.txt; "
            + "while IFS=$'\\t' read -r compound_id smiles; do [[ -n \"$compound_id\" ]] || continue; "
            + "if [[ -s \"artifacts/logs/${compound_id}.result.tsv\" ]]; then cat \"artifacts/logs/${compound_id}.result.tsv\" >> artifacts/docking_scores.csv; "
            + "else printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \"$compound_id\" failed '' '' '' '' 'Compound task produced no result' >> artifacts/docking_scores.csv; fi; "
            + "done < <(tail -n +2 artifacts/input/compounds.tsv)"
        )
    else:
        script = (
            common_setup
            + preparation_loop
            + f"if {dock_command}; then {score_parse}; "
            + "printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \"$compound_id\" completed \"$score\" \"$cnn_score\" \"$cnn_affinity\" \"artifacts/poses/${compound_id}.pdbqt\" '' >> artifacts/docking_scores.csv; "
            + f"else printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \"$compound_id\" failed '' '' '' '' '{failure_message}' >> artifacts/docking_scores.csv; fi; done; "
            + "while IFS=$'\\t' read -r compound_id smiles; do [[ -n \"$compound_id\" ]] || continue; "
            + "if [[ -s \"artifacts/logs/${compound_id}.result.tsv\" ]]; then cat \"artifacts/logs/${compound_id}.result.tsv\" >> artifacts/docking_scores.csv; fi; "
            + "done < <(tail -n +2 artifacts/input/compounds.tsv)"
        )
    return_code = run_docker_job(job.run_dir, manifest, ["bash", "-lc", script], gpu_device=selected_gpu)
    score_path = artifacts / "docking_scores.csv"
    rows: list[dict[str, Any]] = []
    if score_path.is_file():
        with score_path.open(newline="", encoding="utf-8", errors="replace") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                score = _finite_float(row, "best_score_kcal_mol") if row.get("best_score_kcal_mol") else None
                cnn_score = _finite_float(row, "cnn_score") if row.get("cnn_score") else None
                cnn_affinity = _finite_float(row, "cnn_affinity") if row.get("cnn_affinity") else None
                rows.append({
                    "compound_id": row.get("compound_id", ""),
                    "engine": display_engine,
                    "status": row.get("status", "failed"),
                    "best_score_kcal_mol": score,
                    "cnn_score": cnn_score,
                    "cnn_affinity": cnn_affinity,
                    "pose_path": row.get("pose_path", ""),
                    "error": row.get("error", ""),
                })
        _write_csv(score_path, rows, ["compound_id", "engine", "status", "best_score_kcal_mol", "cnn_score", "cnn_affinity", "pose_path", "error"])
    successful = [row for row in rows if row.get("status") == "completed" and row.get("best_score_kcal_mol") is not None]
    success = return_code == 0 and bool(successful)
    best = (
        max(successful, key=lambda item: float(item["cnn_score"]) if item.get("cnn_score") is not None else -math.inf)
        if engine == "gnina" and cnn_scoring != "none" and successful
        else min(successful, key=lambda item: float(item["best_score_kcal_mol"])) if successful
        else None
    )
    finish_job(job.run_dir, success, {
        "outputs": {"scores": "artifacts/docking_scores.csv", "pose_directory": "artifacts/poses"},
        "metrics": {
            "engine": display_engine, "compound_count": len(compounds),
            "scored_count": len(successful), "failed_count": len(rows) - len(successful),
            "primary_ranking_metric": "cnn_score" if engine == "gnina" and cnn_scoring != "none" else "best_score_kcal_mol",
            "primary_ranking_direction": "higher" if engine == "gnina" and cnn_scoring != "none" else "lower",
            "best_compound_id": best["compound_id"] if best else "",
            "best_score_kcal_mol": best.get("best_score_kcal_mol") if best else None,
            "best_cnn_score": best.get("cnn_score") if best else None,
            "return_code": return_code,
        },
        "downstream_artifacts": {"docking_scores": "artifacts/docking_scores.csv", "docking_poses": "artifacts/poses"},
    })
    return job.run_dir


def run_ligand_autodock_gpu(
    target_pdb: Path,
    compound_set_file: Path,
    center_x: float,
    center_y: float,
    center_z: float,
    size_x: float = 20.0,
    size_y: float = 20.0,
    size_z: float = 20.0,
    nrun: int = 20,
    seed: int = 42,
    compound_limit: int = 0,
    cpu_cores: int = 4,
    gpu_device: str = "0",
    align_longest_axis_to_x: bool = False,
    alignment_ligand_sdf: Path | None = None,
    use_scrub: bool = True,
    scrub_ph: float = 7.4,
    scrub_skip_tautomer: bool = True,
) -> Path:
    """Dock a Database compound set with AutoDock-GPU and AutoGrid maps."""
    validate_target_pdb(target_pdb)
    center = [float(center_x), float(center_y), float(center_z)]
    box_size = [float(size_x), float(size_y), float(size_z)]
    if any(not math.isfinite(value) for value in center):
        raise ValueError("Docking box center coordinates must be finite numbers.")
    if any(not math.isfinite(value) or not 1 <= value <= 60 for value in box_size):
        raise ValueError("Docking box dimensions must be between 1 and 60 Å.")
    if type(nrun) is not int or not 1 <= nrun <= 500:
        raise ValueError("AutoDock-GPU LGA runs must be between 1 and 500.")
    if type(seed) is not int or not 1 <= seed <= 2_147_483_646:
        raise ValueError("AutoDock-GPU seed must be a positive 32-bit integer.")
    if type(cpu_cores) is not int or not 1 <= cpu_cores <= 256:
        raise ValueError("AutoDock-GPU CPU allocation must be between 1 and 256 cores.")
    if type(use_scrub) is not bool or type(scrub_skip_tautomer) is not bool:
        raise ValueError("Scrub and tautomer-enumeration settings must be booleans.")
    scrub_ph = float(scrub_ph)
    if not math.isfinite(scrub_ph) or not 0 <= scrub_ph <= 14:
        raise ValueError("Scrub pH must be between 0 and 14.")
    if not str(gpu_device).strip() or str(gpu_device).strip().lower() in {"none", "cpu", "off"}:
        raise ValueError("AutoDock-GPU requires a CUDA GPU device.")
    original_center, original_box_size = list(center), list(box_size)
    alignment_ligand_sdf = Path(alignment_ligand_sdf) if alignment_ligand_sdf else None
    target_data, aligned_center, aligned_size, coordinate_transform = _docking_axis_alignment(
        target_pdb,
        tuple(center),
        tuple(box_size),
        align_longest_axis_to_x=align_longest_axis_to_x,
        alignment_ligand_sdf=alignment_ligand_sdf,
    )
    center, box_size = list(aligned_center), list(aligned_size)
    payload, compounds = _validated_compound_set(compound_set_file, compound_limit)
    params = {
        "engine": "AutoDock-GPU",
        "center": center,
        "box_size": box_size,
        "original_center": original_center,
        "original_box_size": original_box_size,
        "align_longest_axis_to_x": bool(coordinate_transform),
        "coordinate_transform": coordinate_transform,
        "nrun": int(nrun),
        "seed": int(seed),
        "compound_limit": int(compound_limit),
        "cpu_cores": int(cpu_cores),
        "gpu_device": str(gpu_device).strip(),
        "compound_set_name": payload.get("name", "Compound set"),
        "compound_count": len(compounds),
        "primary_ranking_metric": "best_score_kcal_mol",
        "primary_ranking_direction": "lower",
        "score_units": "kcal/mol (AutoDock4 scoring function)",
        "use_scrub": use_scrub,
        "scrub_ph": scrub_ph,
        "scrub_skip_tautomer": scrub_skip_tautomer,
    }
    job = create_job(
        "ligand-docking", "ligand_docking", "autodock_gpu",
        {
            "target_pdb": str(target_pdb), "compound_set_file": str(compound_set_file),
            **({"alignment_ligand_sdf": str(alignment_ligand_sdf)} if coordinate_transform and alignment_ligand_sdf else {}),
        },
        params,
    )
    contract_outputs = [
        {"name": "docking_scores", "type": "ligand_docking_scores/v1", "path": "artifacts/docking_scores.csv", "score_units": "kcal/mol (AutoDock4 scoring function)", "ranking_direction": "lower"},
        {"name": "poses", "type": "ligand_pose_collection/v1", "path": "artifacts/poses"},
    ]
    if coordinate_transform:
        contract_outputs.append({"name": "coordinate_transform", "type": "coordinate_transform/v1", "path": "artifacts/input/alignment_transform.json"})
    _write_job_contract(
        job.run_dir, "ligand.docking.autodock_gpu", "compound_set/v1", "ligand_docking_result/v1",
        payload, compounds, target_artifact="artifacts/input/receptor.pdb",
        outputs=contract_outputs,
        extra_inputs=(
            [{"name": "axis_alignment_ligand", "type": "ligand_coordinates/sdf/v1", "path": "artifacts/input/alignment_ligand.sdf"}]
            if coordinate_transform
            else []
        ),
    )
    input_root = job.run_dir / "artifacts" / "input"
    input_root.mkdir(parents=True, exist_ok=True)
    (input_root / "receptor.pdb").write_text(target_data, encoding="utf-8")
    if coordinate_transform and alignment_ligand_sdf is not None:
        shutil.copy2(alignment_ligand_sdf, input_root / "alignment_ligand.sdf")
        write_json(input_root / "alignment_transform.json", coordinate_transform)
    _write_compound_inputs(job.run_dir, payload, compounds)
    with (input_root / "compounds.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["compound_id", "smiles"])
        writer.writerows((item["compound_id"], item["smiles"]) for item in compounds)
    write_json(input_root / "autodock_gpu_config.json", {
        "center": center,
        "box_size": box_size,
        "nrun": int(nrun),
        "seed": int(seed),
        "use_scrub": use_scrub,
        "scrub_ph": scrub_ph,
        "scrub_skip_tautomer": scrub_skip_tautomer,
    })
    (input_root / "smiles_to_3d.py").write_text(_RDKIT_3D_SCRIPT, encoding="utf-8")
    shutil.copy2(Path(__file__).with_name("autodock_gpu_runner.py"), input_root / "autodock_gpu_runner.py")
    return_code = run_docker_job(
        job.run_dir,
        load_manifest("autodock_gpu"),
        ["python", "/work/artifacts/input/autodock_gpu_runner.py"],
        gpu_device=str(gpu_device).strip(),
    )
    score_path = job.run_dir / "artifacts" / "docking_scores.csv"
    rows: list[dict[str, Any]] = []
    if score_path.is_file():
        with score_path.open(newline="", encoding="utf-8", errors="replace") as handle:
            for row in csv.DictReader(handle):
                rows.append({
                    "compound_id": row.get("compound_id", ""),
                    "engine": "AutoDock-GPU",
                    "status": row.get("status", "failed"),
                    "best_score_kcal_mol": _finite_float(row, "best_score_kcal_mol"),
                    "pose_path": row.get("pose_path", ""),
                    "error": row.get("error", ""),
                })
        _write_csv(score_path, rows, ["compound_id", "engine", "status", "best_score_kcal_mol", "pose_path", "error"])
    successful = [row for row in rows if row["status"] == "completed" and row["best_score_kcal_mol"] is not None]
    best = min(successful, key=lambda row: float(row["best_score_kcal_mol"])) if successful else None
    summary_path = job.run_dir / "artifacts" / "autodock_gpu_result.json"
    summary = read_json(summary_path) if summary_path.is_file() else {}
    success = return_code == 0 and summary.get("success") is True and bool(successful)
    finish_job(job.run_dir, success, {
        "outputs": {"scores": "artifacts/docking_scores.csv", "pose_directory": "artifacts/poses", "native_run_summary": "artifacts/autodock_gpu_result.json"},
        "metrics": {
            "engine": "AutoDock-GPU",
            "score_units": "kcal/mol (AutoDock4 scoring function)",
            "compound_count": len(compounds),
            "scored_count": len(successful),
            "failed_count": len(rows) - len(successful),
            "primary_ranking_metric": "best_score_kcal_mol",
            "primary_ranking_direction": "lower",
            "best_compound_id": best["compound_id"] if best else "",
            "best_score_kcal_mol": best["best_score_kcal_mol"] if best else None,
            "return_code": return_code,
        },
        "downstream_artifacts": {"docking_scores": "artifacts/docking_scores.csv", "docking_poses": "artifacts/poses"},
    })
    return job.run_dir


def run_ligand_vina(
    target_pdb: Path,
    compound_set_file: Path,
    center_x: float,
    center_y: float,
    center_z: float,
    size_x: float = 20.0,
    size_y: float = 20.0,
    size_z: float = 20.0,
    exhaustiveness: int = 32,
    num_modes: int = 10,
    seed: int = 42,
    compound_limit: int = 0,
    cpu_cores: int = 4,
    gpu_device: str = "none",
    align_longest_axis_to_x: bool = False,
    alignment_ligand_sdf: Path | None = None,
    job_name: str = "autodock-vina",
    use_scrub: bool = True,
    scrub_ph: float = 7.4,
    scrub_skip_tautomer: bool = True,
) -> Path:
    """Dock a Studio Database compound set with AutoDock Vina."""
    return _run_docking(
        engine="vina",
        target_pdb=target_pdb, compound_set_file=compound_set_file,
        center_x=center_x, center_y=center_y, center_z=center_z,
        size_x=size_x, size_y=size_y, size_z=size_z,
        exhaustiveness=exhaustiveness, num_modes=num_modes, seed=seed,
        compound_limit=compound_limit, cpu_cores=cpu_cores, gpu_device="none",
        align_longest_axis_to_x=align_longest_axis_to_x,
        alignment_ligand_sdf=alignment_ligand_sdf,
        job_name=job_name,
        use_scrub=use_scrub,
        scrub_ph=scrub_ph,
        scrub_skip_tautomer=scrub_skip_tautomer,
    )


def run_ligand_gnina(
    target_pdb: Path,
    compound_set_file: Path,
    center_x: float,
    center_y: float,
    center_z: float,
    size_x: float = 20.0,
    size_y: float = 20.0,
    size_z: float = 20.0,
    exhaustiveness: int = 8,
    num_modes: int = 9,
    seed: int = 42,
    compound_limit: int = 0,
    cpu_cores: int = 4,
    compound_workers: int = 2,
    gpu_device: str = "0",
    align_longest_axis_to_x: bool = False,
    alignment_ligand_sdf: Path | None = None,
    job_name: str | None = None,
    cnn_scoring: str = "rescore",
    use_scrub: bool = True,
    scrub_ph: float = 7.4,
    scrub_skip_tautomer: bool = True,
) -> Path:
    """Run GNINA docking and CNN rescoring on a Studio Database compound set."""
    return _run_docking(
        engine="gnina",
        target_pdb=target_pdb, compound_set_file=compound_set_file,
        center_x=center_x, center_y=center_y, center_z=center_z,
        size_x=size_x, size_y=size_y, size_z=size_z,
        exhaustiveness=exhaustiveness, num_modes=num_modes, seed=seed,
        compound_limit=compound_limit, cpu_cores=cpu_cores, gpu_device=gpu_device,
        compound_workers=compound_workers,
        align_longest_axis_to_x=align_longest_axis_to_x,
        alignment_ligand_sdf=alignment_ligand_sdf,
        job_name=job_name,
        cnn_scoring=cnn_scoring,
        use_scrub=use_scrub,
        scrub_ph=scrub_ph,
        scrub_skip_tautomer=scrub_skip_tautomer,
    )


def run_ligand_unidock_pro(
    target_pdb: Path,
    compound_set_file: Path,
    center_x: float,
    center_y: float,
    center_z: float,
    size_x: float = 20.0,
    size_y: float = 20.0,
    size_z: float = 20.0,
    search_mode: str = "balance",
    num_modes: int = 9,
    seed: int = 42,
    compound_limit: int = 0,
    cpu_cores: int = 4,
    gpu_device: str = "0",
    align_longest_axis_to_x: bool = False,
    alignment_ligand_sdf: Path | None = None,
    job_name: str = "unidock-pro",
    use_scrub: bool = True,
    scrub_ph: float = 7.4,
    scrub_skip_tautomer: bool = True,
) -> Path:
    """Run the classic Uni-Dock Pro receptor-docking mode on a compound set."""
    search_mode = str(search_mode).strip().lower()
    if search_mode not in {"fast", "balance", "detail"}:
        raise ValueError("Uni-Dock Pro search mode must be fast, balance, or detail.")
    if type(num_modes) is not int or not 1 <= num_modes <= 100:
        raise ValueError("Uni-Dock Pro pose count must be between 1 and 100.")
    if type(seed) is not int or not 1 <= seed <= 2_147_483_646:
        raise ValueError("Uni-Dock Pro seed must be a positive 32-bit integer.")
    if not str(gpu_device).strip() or str(gpu_device).strip().lower() in {"none", "cpu", "off"}:
        raise ValueError("Uni-Dock Pro requires a CUDA GPU device.")
    validate_target_pdb(target_pdb)
    if any(not math.isfinite(value) for value in (center_x, center_y, center_z)):
        raise ValueError("Docking box center coordinates must be finite numbers.")
    if any(not math.isfinite(value) or not 1 <= value <= 60 for value in (size_x, size_y, size_z)):
        raise ValueError("Docking box dimensions must be between 1 and 60 Å.")
    if type(cpu_cores) is not int or not 1 <= cpu_cores <= 256:
        raise ValueError("cpu_cores must be between 1 and 256.")
    if type(use_scrub) is not bool or type(scrub_skip_tautomer) is not bool:
        raise ValueError("Scrub and tautomer-enumeration settings must be booleans.")
    scrub_ph = float(scrub_ph)
    if not math.isfinite(scrub_ph) or not 0 <= scrub_ph <= 14:
        raise ValueError("Scrub pH must be between 0 and 14.")

    original_center = [float(center_x), float(center_y), float(center_z)]
    original_box_size = [float(size_x), float(size_y), float(size_z)]
    alignment_ligand_sdf = Path(alignment_ligand_sdf) if alignment_ligand_sdf else None
    target_data, aligned_center, aligned_size, coordinate_transform = _docking_axis_alignment(
        target_pdb,
        tuple(original_center),
        tuple(original_box_size),
        align_longest_axis_to_x=align_longest_axis_to_x,
        alignment_ligand_sdf=alignment_ligand_sdf,
    )
    center_x, center_y, center_z = aligned_center
    size_x, size_y, size_z = aligned_size
    payload, compounds = _validated_compound_set(compound_set_file, compound_limit)
    params = {
        "engine": "unidock_pro", "search_mode": search_mode,
        "center": [float(center_x), float(center_y), float(center_z)],
        "box_size": [float(size_x), float(size_y), float(size_z)],
        "original_center": original_center,
        "original_box_size": original_box_size,
        "align_longest_axis_to_x": bool(coordinate_transform),
        "coordinate_transform": coordinate_transform,
        "num_modes": num_modes, "seed": seed, "compound_limit": compound_limit,
        "cpu_cores": cpu_cores, "gpu_device": str(gpu_device).strip(),
        "compound_set_name": payload.get("name", "Compound set"),
        "job_name": str(job_name).strip()[:120] or "unidock-pro",
        "use_scrub": use_scrub,
        "scrub_ph": scrub_ph,
        "scrub_skip_tautomer": scrub_skip_tautomer,
    }
    job = create_job(
        "ligand-docking", "ligand_docking", "unidock_pro",
        {
            "target_pdb": str(target_pdb), "compound_set_file": str(compound_set_file),
            **({"alignment_ligand_sdf": str(alignment_ligand_sdf)} if coordinate_transform and alignment_ligand_sdf else {}),
        }, params,
    )
    contract_outputs = [
        {"name": "docking_scores", "type": "ligand_docking_scores/v1", "path": "artifacts/docking_scores.csv"},
        {"name": "poses", "type": "ligand_pose_collection/v1", "path": "artifacts/poses"},
    ]
    if coordinate_transform:
        contract_outputs.append({"name": "coordinate_transform", "type": "coordinate_transform/v1", "path": "artifacts/input/alignment_transform.json"})
    _write_job_contract(
        job.run_dir, "ligand.docking.unidock_pro", "compound_set/v1", "ligand_docking_result/v1",
        payload, compounds, target_artifact="artifacts/input/receptor.pdb",
        outputs=contract_outputs,
        extra_inputs=(
            [{"name": "axis_alignment_ligand", "type": "ligand_coordinates/sdf/v1", "path": "artifacts/input/alignment_ligand.sdf"}]
            if coordinate_transform
            else []
        ),
    )
    artifacts = job.run_dir / "artifacts"
    input_dir = artifacts / "input"
    input_dir.mkdir(parents=True, exist_ok=True)
    (input_dir / "receptor.pdb").write_text(target_data, encoding="utf-8")
    if coordinate_transform and alignment_ligand_sdf is not None:
        shutil.copy2(alignment_ligand_sdf, input_dir / "alignment_ligand.sdf")
        write_json(input_dir / "alignment_transform.json", coordinate_transform)
    staged_compounds = _write_compound_inputs(job.run_dir, payload, compounds)
    (input_dir / "smiles_to_3d.py").write_text(_RDKIT_3D_SCRIPT, encoding="utf-8")
    with (input_dir / "compounds.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["compound_id", "smiles"])
        writer.writerows((item["compound_id"], item["smiles"]) for item in staged_compounds)
    config_path = input_dir / "docking.conf"
    config_path.write_text("".join(
        f"{key} = {value:.3f}\n"
        for key, value in zip(
            ("center_x", "center_y", "center_z", "size_x", "size_y", "size_z"),
            (center_x, center_y, center_z, size_x, size_y, size_z),
        )
    ), encoding="utf-8")
    scrub_step = (
        "if scrub.py \"artifacts/ligands/${compound_id}_source.sdf\" "
        "-o \"artifacts/ligands/${compound_id}.sdf\" "
        f"--ph {scrub_ph:.2f}"
        + (" --skip_tautomer" if scrub_skip_tautomer else "")
        + " >> artifacts/logs/preparation.log 2>&1 "
        + "&& [[ -s \"artifacts/ligands/${compound_id}.sdf\" ]]; then :; else "
        + "cp \"artifacts/ligands/${compound_id}_source.sdf\" \"artifacts/ligands/${compound_id}.sdf\"; "
        + "echo 'Scrub failed; using the RDKit SDF for Meeko preparation.' >> artifacts/logs/preparation.log; fi; "
        if use_scrub else
        "cp \"artifacts/ligands/${compound_id}_source.sdf\" \"artifacts/ligands/${compound_id}.sdf\"; "
    )
    script = (
        "set -euo pipefail; cd /work; "
        "mkdir -p artifacts/prepared artifacts/ligands artifacts/poses artifacts/logs; "
        "if ! mk_prepare_receptor.py --read_pdb artifacts/input/receptor.pdb -o artifacts/prepared/receptor -p; then "
        "awk 'substr($0,77,2) !~ /^[[:space:]]*H[[:space:]]*$/' artifacts/input/receptor.pdb > artifacts/input/receptor-noh.pdb; "
        "mk_prepare_receptor.py -i artifacts/input/receptor-noh.pdb -o artifacts/prepared/receptor -p; fi; "
        ": > artifacts/input/ligand_index.txt; "
        "tail -n +2 artifacts/input/compounds.tsv | while IFS=$'\\t' read -r compound_id smiles; do "
        "[[ -n \"$compound_id\" && -n \"$smiles\" ]] || continue; "
        "if ! python artifacts/input/smiles_to_3d.py --smiles \"$smiles\" --name \"$compound_id\" "
        "--output \"artifacts/ligands/${compound_id}_source.sdf\" >> artifacts/logs/preparation.log 2>&1; then continue; fi; "
        + scrub_step
        + "if ! mk_prepare_ligand.py -i \"artifacts/ligands/${compound_id}.sdf\" "
        "-o \"artifacts/ligands/${compound_id}.pdbqt\" >> artifacts/logs/preparation.log 2>&1; then continue; fi; "
        "printf '%s\\n' \"artifacts/ligands/${compound_id}.pdbqt\" >> artifacts/input/ligand_index.txt; done; "
        "dock_status=0; if [[ -s artifacts/input/ligand_index.txt ]]; then "
        f"udp --receptor artifacts/prepared/receptor.pdbqt --ligand_index artifacts/input/ligand_index.txt "
        f"--config artifacts/input/docking.conf --dir artifacts/poses --search_mode {search_mode} "
        f"--num_modes {num_modes} --seed {seed} || dock_status=$?; fi; "
        "printf 'compound_id\\tstatus\\tbest_score_kcal_mol\\tcnn_score\\tcnn_affinity\\tpose_path\\terror\\n' > artifacts/docking_scores.csv; "
        "tail -n +2 artifacts/input/compounds.tsv | while IFS=$'\\t' read -r compound_id smiles; do "
        "[[ -n \"$compound_id\" ]] || continue; pose=\"artifacts/poses/${compound_id}_out.pdbqt\"; "
        "if [[ -s \"$pose\" ]]; then "
        "score=$(awk '/^REMARK VINA RESULT:/ {print $4; exit} /^REMARK minimizedAffinity / {print $3; exit}' \"$pose\"); "
        "if [[ -n \"$score\" ]]; then printf '%s\\tcompleted\\t%s\\t\\t\\t%s\\t\\n' \"$compound_id\" \"$score\" \"$pose\" >> artifacts/docking_scores.csv; "
        "else printf '%s\\tfailed\\t\\t\\t\\t\\tPose file has no readable docking score\\n' \"$compound_id\" >> artifacts/docking_scores.csv; fi; "
        "else printf '%s\\tfailed\\t\\t\\t\\t\\tNo pose returned (dock exit %s)\\n' \"$compound_id\" \"$dock_status\" >> artifacts/docking_scores.csv; fi; done; "
        "exit 0"
    )
    return_code = run_docker_job(job.run_dir, load_manifest("unidock_pro"), ["bash", "-lc", script], gpu_device=str(gpu_device).strip())
    score_path = artifacts / "docking_scores.csv"
    rows: list[dict[str, Any]] = []
    if score_path.is_file():
        with score_path.open(newline="", encoding="utf-8", errors="replace") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                rows.append({
                    "compound_id": row.get("compound_id", ""), "engine": "Uni-Dock Pro",
                    "status": row.get("status", "failed"),
                    "best_score_kcal_mol": _finite_float(row, "best_score_kcal_mol") if row.get("best_score_kcal_mol") else None,
                    "cnn_score": None, "cnn_affinity": None,
                    "pose_path": row.get("pose_path", ""), "error": row.get("error", ""),
                })
        _write_csv(score_path, rows, ["compound_id", "engine", "status", "best_score_kcal_mol", "cnn_score", "cnn_affinity", "pose_path", "error"])
    successful = [row for row in rows if row.get("status") == "completed" and row.get("best_score_kcal_mol") is not None]
    best = min(successful, key=lambda item: float(item["best_score_kcal_mol"])) if successful else None
    finish_job(job.run_dir, return_code == 0 and bool(successful), {
        "outputs": {"scores": "artifacts/docking_scores.csv", "pose_directory": "artifacts/poses"},
        "metrics": {
            "engine": "Uni-Dock Pro", "search_mode": search_mode,
            "compound_count": len(compounds), "scored_count": len(successful),
            "failed_count": len(rows) - len(successful), "primary_ranking_metric": "best_score_kcal_mol",
            "primary_ranking_direction": "lower", "best_compound_id": best["compound_id"] if best else "",
            "best_score_kcal_mol": best.get("best_score_kcal_mol") if best else None,
            "return_code": return_code,
        },
        "downstream_artifacts": {"docking_scores": "artifacts/docking_scores.csv", "docking_poses": "artifacts/poses"},
    })
    return job.run_dir


def _rosettaligand_protocol_xml(protocol: str, padding: float, reference_guided: bool) -> str:
    runmode = "VSH" if protocol == "convergence" else protocol.upper()
    reference_options = "" if reference_guided else ' reference_pool="none" use_pharmacophore="0"'
    if protocol == "convergence":
        protocol_options = ' final_optH_mode="1" pre_optH_relax="1" entropy_method="Simple" nrelax="30" nreport="1" final_exact_minimize="sc"'
        stages = '''
      <Stage repeats="8" npool="150" rmsdthreshold="2.0" pmut="0.2"
        maxiter="40" pack_cycles="75" smoothing="0.75"
        ramp_schedule="0.1,1.0"/>
      <Stage repeats="8" npool="150" rmsdthreshold="2.0" pmut="0.2"
        maxiter="40" pack_cycles="75" smoothing="0.375"
        ramp_schedule="0.1,1.0"/>'''
    elif protocol == "vsh":
        protocol_options = ' final_optH_mode="1" pre_optH_relax="1" entropy_method="Simple"'
        stages = ""
    else:
        protocol_options = ' aligner_fastmode="1" final_exact_minimize="ligandonly"'
        stages = ""
    return f'''<ROSETTASCRIPTS>
  <SCOREFXNS>
    <ScoreFunction name="genpot_soft" weights="beta_cart">
      <Reweight scoretype="fa_rep" weight="0.2"/>
    </ScoreFunction>
    <ScoreFunction name="genpot" weights="beta_cart">
      <Reweight scoretype="coordinate_constraint" weight="0.1"/>
    </ScoreFunction>
  </SCOREFXNS>
  <MOVERS>
    <GALigandDock name="dock" scorefxn="genpot_soft" scorefxn_relax="genpot"
      runmode="{runmode}" premin_ligand="1"
      multiple_ligands_file="%%liglist%%" estimate_dG="1"
      padding="{padding:.3f}"{reference_options}{protocol_options}>{stages}
    </GALigandDock>
  </MOVERS>
  <PROTOCOLS>
    <Add mover="dock"/>
  </PROTOCOLS>
  <OUTPUT scorefxn="genpot"/>
</ROSETTASCRIPTS>
'''


def _parse_rosettaligand_scores(run_dir: Path, protocol: str, compounds: list[dict[str, Any]]) -> list[dict[str, Any]]:
    compound_ids = {str(item["compound_id"]) for item in compounds}
    raw_rows: list[dict[str, Any]] = []
    for path in sorted((run_dir / "native").rglob("run.score.sc")):
        header: list[str] = []
        for line in path.read_text(errors="replace").splitlines():
            if not line.startswith("SCORE:"):
                continue
            values = line.split()[1:]
            if values and values[0] == "total_score":
                header = values
                continue
            if not header or len(values) != len(header):
                continue
            source = dict(zip(header, values))
            compound_id = str(source.get("ligandname") or "")
            if compound_id not in compound_ids:
                continue
            replicate_match = re.search(r"replicate_(\d+)", str(path.relative_to(run_dir)))
            raw_rows.append({
                "compound_id": compound_id,
                "protocol": protocol.upper(),
                "replicate": int(replicate_match.group(1)) if replicate_match else 1,
                "estimated_dg_reu": _finite_float(source, "dG"),
                "total_score_reu": _finite_float(source, "total_score"),
                "ligand_rmsd_angstrom": _finite_float(source, "lig_rms"),
                "description": str(source.get("description") or ""),
                "native_score_file": path.relative_to(run_dir).as_posix(),
            })
    pose_files = sorted((run_dir / "poses").glob("*.pdb")) if (run_dir / "poses").is_dir() else []
    used_poses: set[Path] = set()
    pose_root = run_dir / "artifacts" / "poses"
    pose_root.mkdir(parents=True, exist_ok=True)
    ranks: dict[str, int] = {}
    rows: list[dict[str, Any]] = []
    for index, row in enumerate(raw_rows):
        compound_id = row["compound_id"]
        ranks[compound_id] = ranks.get(compound_id, 0) + 1
        row["pose_rank"] = ranks[compound_id]
        description = row["description"]
        matching = [path for path in pose_files if description and description in path.name]
        source_pose = next((path for path in matching if path not in used_poses), None)
        if source_pose is None and index < len(pose_files) and pose_files[index] not in used_poses:
            source_pose = pose_files[index]
        if source_pose is not None:
            used_poses.add(source_pose)
            target_pose = pose_root / f"{compound_id}.replicate_{row['replicate']:03d}.pose_{row['pose_rank']:02d}.pdb"
            shutil.copy2(source_pose, target_pose)
            row["pose_path"] = target_pose.relative_to(run_dir).as_posix()
        else:
            row["pose_path"] = ""
        rows.append(row)
    return rows


def run_ligand_rosetta_ligand(
    target_pdb: Path,
    compound_set_file: Path,
    center_x: float,
    center_y: float,
    center_z: float,
    protocol: str = "vsh",
    reference_mode: str = "pocket_center",
    reference_ligand_sdf: Path | None = None,
    compound_limit: int = 0,
    cpu_cores: int = 8,
    ph: float = 7.4,
    conformers: int = 20,
    minimization_steps: int = 2000,
    padding: float = 4.0,
    replicates: int = 1,
    seed_start: int = 1001,
    job_name: str = "rosettaligand",
) -> Path:
    """Dock a Database compound set with Studio's bundled RosettaLigand runner."""
    protocol = str(protocol).strip().lower()
    reference_mode = str(reference_mode).strip().lower()
    parameters = {
        "target_pdb": target_pdb, "compound_set_file": compound_set_file,
        "center_x": center_x, "center_y": center_y, "center_z": center_z,
        "protocol": protocol, "reference_mode": reference_mode,
        "reference_ligand_sdf": reference_ligand_sdf,
        "compound_limit": compound_limit, "cpu_cores": cpu_cores,
        "ph": ph, "conformers": conformers,
        "minimization_steps": minimization_steps, "padding": padding,
        "replicates": replicates, "seed_start": seed_start,
        "job_name": str(job_name).strip()[:120] or "rosettaligand",
    }
    validate_ligand_tool_request("ligand.docking.rosetta_ligand", parameters)
    if any(not math.isfinite(float(value)) for value in (center_x, center_y, center_z)):
        raise ValueError("RosettaLigand pocket-center coordinates must be finite numbers.")
    payload, compounds = _validated_compound_set(compound_set_file, compound_limit)
    if type(conformers) is not int or not 1 <= conformers <= 1000:
        raise ValueError("RosettaLigand conformers must be between 1 and 1,000.")
    if type(minimization_steps) is not int or not 1 <= minimization_steps <= 100_000:
        raise ValueError("RosettaLigand minimization steps must be between 1 and 100,000.")
    if type(seed_start) is not int or seed_start < 1 or seed_start + replicates >= 2_147_483_647:
        raise ValueError("RosettaLigand seed range must contain positive 32-bit integers.")
    source = _openvs_runner_constants()
    image = _ROSETTALIGAND_IMAGE
    params = {
        "engine": "rosetta_ligand", "protocol": protocol, "reference_mode": reference_mode,
        "compound_limit": compound_limit, "compound_count": len(compounds),
        "center": [float(center_x), float(center_y), float(center_z)],
        "cpu_cores": cpu_cores, "ph": float(ph), "conformers": conformers,
        "minimization_steps": minimization_steps, "padding_angstrom": float(padding),
        "replicates": replicates, "seed_start": seed_start, "docker_image": image,
    }
    input_payload = {"target_pdb": str(target_pdb), "compound_set_file": str(compound_set_file)}
    if reference_ligand_sdf is not None:
        input_payload["reference_ligand_sdf"] = str(reference_ligand_sdf)
    job = create_job("ligand-docking", "ligand_docking", "rosetta_ligand", input_payload, params)
    _write_job_contract(
        job.run_dir, "ligand.docking.rosetta_ligand", "compound_set/v1", "ligand_docking_result/v1",
        payload, compounds, target_artifact="artifacts/input/receptor.pdb",
        outputs=[
            {"name": "docking_scores", "type": "ligand_docking_scores/v1", "path": "artifacts/docking_scores.csv", "score_units": "REU"},
            {"name": "poses", "type": "ligand_pose_collection/v1", "path": "artifacts/poses"},
        ],
    )
    run_root = job.run_dir
    input_dir = run_root / "input"
    input_dir.mkdir(parents=True, exist_ok=True)
    (input_dir / "receptor.pdb").write_bytes(Path(target_pdb).read_bytes())
    with (input_dir / "compounds.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["compound_id", "smiles"])
        writer.writerows((item["compound_id"], item["smiles"]) for item in compounds)
    (input_dir / "dock.xml").write_text(_rosettaligand_protocol_xml(protocol, float(padding), reference_mode == "reference_guided"), encoding="utf-8")
    (input_dir / "center_pdb.py").write_text(source["_CENTER_PDB_SCRIPT"], encoding="utf-8")
    (input_dir / "validate_rosetta_params.py").write_text(source["_ROSETTA_PARAMS_VALIDATOR_SCRIPT"], encoding="utf-8")
    (input_dir / "smiles_to_3d.py").write_text(_RDKIT_3D_SCRIPT, encoding="utf-8")
    runner = source["_RUNNER_SCRIPT"].replace("cd /workspace", "cd /work", 1)
    (input_dir / "run_openvs.sh").write_text(runner, encoding="utf-8")
    if reference_mode == "reference_guided":
        (input_dir / "reference_ligand.sdf").write_bytes(Path(reference_ligand_sdf).read_bytes())
    environment = {
        "OPENVS_PROTOCOL": protocol.upper(), "OPENVS_REFERENCE_MODE": reference_mode,
        "OPENVS_REFERENCE_SUFFIX": ".sdf" if reference_ligand_sdf else "",
        "OPENVS_CPU_WORKERS": str(cpu_cores), "OPENVS_REPLICATES": str(replicates),
        "OPENVS_SEED_START": str(seed_start), "OPENVS_PH": f"{float(ph):.2f}",
        "OPENVS_PRESERVE_INPUT_PROTONATION": "0", "OPENVS_CONFORMERS": str(conformers),
        "OPENVS_MINIMIZATION_STEPS": str(minimization_steps),
        "OPENVS_CENTER_X": f"{float(center_x):.4f}", "OPENVS_CENTER_Y": f"{float(center_y):.4f}",
        "OPENVS_CENTER_Z": f"{float(center_z):.4f}",
    }
    shell = "set -euo pipefail; cd /work; " + "; ".join(
        f"export {key}={shlex.quote(value)}" for key, value in environment.items()
    ) + "; bash input/run_openvs.sh"
    return_code = run_docker_job(run_root, load_manifest("rosetta_ligand"), ["bash", "-lc", shell], gpu_device="none")
    rows = _parse_rosettaligand_scores(run_root, protocol, compounds)
    _write_csv(
        run_root / "artifacts" / "docking_scores.csv", rows,
        ["compound_id", "protocol", "replicate", "pose_rank", "estimated_dg_reu", "total_score_reu", "ligand_rmsd_angstrom", "pose_path", "native_score_file"],
    )
    scored = [row for row in rows if row.get("estimated_dg_reu") is not None]
    best = min(scored, key=lambda row: float(row["estimated_dg_reu"])) if scored else None
    finish_job(run_root, return_code == 0 and bool(scored), {
        "outputs": {"scores": "artifacts/docking_scores.csv", "pose_directory": "artifacts/poses"},
        "metrics": {
            "engine": "RosettaLigand", "protocol": protocol.upper(),
            "score_units": "Rosetta energy units (REU; relative ranking only)",
            "compound_count": len(compounds), "scored_pose_count": len(scored),
            "primary_ranking_metric": "estimated_dg_reu", "primary_ranking_direction": "lower",
            "best_compound_id": best.get("compound_id") if best else "",
            "best_estimated_dg_reu": best.get("estimated_dg_reu") if best else None,
            "return_code": return_code,
        },
        "downstream_artifacts": {"docking_scores": "artifacts/docking_scores.csv", "docking_poses": "artifacts/poses"},
    })
    return run_root


def _write_generation_contract(
    job_dir: Path,
    *,
    workflow: str,
    engine: str,
    inputs: list[dict[str, Any]],
    checkpoint: str,
) -> None:
    metadata_path = Path(job_dir) / "metadata.json"
    metadata = read_json(metadata_path)
    input_payload = read_json(Path(job_dir) / "input.json")
    params = input_payload.get("params") if isinstance(input_payload.get("params"), dict) else {}
    metadata.update({
        "workflow": workflow,
        "engine": engine,
        "result_type": "generated_compound_set/v1",
        "artifact_contract_version": 1,
        "requested_compound_count": params.get("requested_compound_count"),
        "checkpoint": checkpoint,
    })
    write_json(metadata_path, metadata)
    write_json(Path(job_dir) / "artifact_contract.json", {
        "schema_version": 1,
        "workflow": workflow,
        "input_type": "protein_pocket_and_generation_conditions/v1",
        "output_type": "generated_compound_set/v1",
        "inputs": inputs,
        "outputs": [
            {"name": "compound_table", "type": "compound_table/v1", "path": "artifacts/generated_compounds.csv"},
            {"name": "compound_set", "type": "compound_set/v1", "path": "artifacts/generated_compound_set.json"},
            {"name": "coordinate_set", "type": "ligand_coordinate_set/sdf/v1", "path": "artifacts/generated_compounds.sdf"},
            {"name": "generation_report", "type": "generation_report/v1", "path": "normalized/generation_report.json"},
        ],
        "compound_count": 0,
        "checkpoint": checkpoint,
    })


def _run_ligand_generation(
    *,
    engine: str,
    workflow: str,
    target_pdb: Path,
    pocket_pdb: Path | None,
    reference_ligand_sdf: Path | None,
    requested_count: int,
    batch_size: int,
    seed: int,
    gpu_device: str,
    cpu_cores: int,
    n_steps: int = 100,
    sampler: str = "ForwardEuler",
    sampling_noise: float = 5.0,
    diffusion_steps: int = 100,
    ligand_atoms_mean: float = 28.0,
    ligand_atoms_std: float = 2.0,
) -> Path:
    parameters = {
        "target_pdb": target_pdb,
        "requested_count": requested_count,
        "batch_size": batch_size,
        "seed": seed,
        "gpu_device": gpu_device,
    }
    if engine == "lddm":
        parameters.update({"reference_ligand_sdf": reference_ligand_sdf, "n_steps": n_steps, "sampler": sampler, "sampling_noise": sampling_noise})
        checkpoint_relative = _LDDM_CHECKPOINT
        input_payload = {"target_pdb": str(target_pdb), "reference_ligand_sdf": str(reference_ligand_sdf)}
        tool = "lddm_generation"
        task_group = "ligand-generation"
    elif engine == "pocketxmol":
        parameters.update({"pocket_pdb": pocket_pdb, "diffusion_steps": diffusion_steps, "ligand_atoms_mean": ligand_atoms_mean, "ligand_atoms_std": ligand_atoms_std})
        checkpoint_relative = _POCKETXMOL_CHECKPOINT
        input_payload = {"target_pdb": str(target_pdb), "pocket_pdb": str(pocket_pdb)}
        tool = "pocketxmol_generation"
        task_group = "ligand-generation"
    else:
        raise ValueError(f"Unsupported ligand generation engine: {engine}.")
    if type(cpu_cores) is not int or not 1 <= cpu_cores <= 256:
        raise ValueError("cpu_cores must be between 1 and 256.")
    _validate_generation_request(workflow, parameters, require_files=True)
    reference = reference_root().expanduser().resolve()
    checkpoint = reference / checkpoint_relative
    manifest = load_manifest(engine)
    params = {
        "engine": engine,
        "requested_compound_count": int(requested_count),
        "batch_size": int(batch_size),
        "seed": int(seed),
        "gpu_device": str(gpu_device),
        "cpu_cores": int(cpu_cores),
        "checkpoint": checkpoint_relative,
    }
    if engine == "lddm":
        params.update({"n_steps": int(n_steps), "sampler": sampler, "sampling_noise": float(sampling_noise)})
    else:
        params.update({"diffusion_steps": int(diffusion_steps), "ligand_atoms_mean": float(ligand_atoms_mean), "ligand_atoms_std": float(ligand_atoms_std)})
    job = create_job(task_group, "ligand_generation", tool, input_payload, params)
    input_root = job.run_dir / "artifacts" / "input"
    input_root.mkdir(parents=True, exist_ok=True)
    target_artifact = input_root / "target.pdb"
    target_artifact.write_bytes(Path(target_pdb).read_bytes())
    contract_inputs = [{"name": "target_structure", "type": "protein_structure/pdb/v1", "path": "artifacts/input/target.pdb"}]
    if engine == "lddm":
        reference_artifact = input_root / "reference_ligand.sdf"
        reference_artifact.write_bytes(Path(reference_ligand_sdf).read_bytes())
        contract_inputs.append({"name": "reference_ligand", "type": "ligand_coordinates/sdf/v1", "path": "artifacts/input/reference_ligand.sdf"})
    else:
        pocket_artifact = input_root / "pocket.pdb"
        pocket_artifact.write_bytes(Path(pocket_pdb).read_bytes())
        contract_inputs.append({"name": "pocket", "type": "protein_pocket/pdb/v1", "path": "artifacts/input/pocket.pdb"})
    _write_generation_contract(job.run_dir, workflow=workflow, engine=engine, inputs=contract_inputs, checkpoint=checkpoint_relative)

    output_root = job.run_dir / "normalized"
    native_root = job.run_dir / "native"
    output_root.mkdir(parents=True, exist_ok=True)
    native_root.mkdir(parents=True, exist_ok=True)
    image = str(manifest["image"])
    run_name = f"{engine}_{job.run_id[-8:]}"
    arguments = [
        "--output", "/work/native", "--normalized-output", "/work/normalized",
        "--count", str(int(requested_count)), "--batch-size", str(int(batch_size)),
        "--seed", str(int(seed)), "--run-name", run_name,
        "--checkpoint", f"/references/{checkpoint_relative}",
        "--native-mode", "fixed_protein_ligand_denovo_condensed",
    ]
    if engine == "lddm":
        arguments.extend([
            "--target", "/work/artifacts/input/target.pdb",
            "--reference-ligand", "/work/artifacts/input/reference_ligand.sdf",
            "--n-steps", str(int(n_steps)), "--sampler", sampler,
            "--sampling-noise", str(float(sampling_noise)),
        ])
    else:
        arguments.extend([
            "--target", "/work/artifacts/input/target.pdb",
            "--pocket-structure", "/work/artifacts/input/pocket.pdb",
            "--diffusion-steps", str(int(diffusion_steps)),
            "--ligand-atoms-mean", str(float(ligand_atoms_mean)),
            "--ligand-atoms-std", str(float(ligand_atoms_std)),
        ])
    command = [
        "docker", "run", "--rm", *docker_gpu_args(gpu_device),
        "--cpus", str(int(cpu_cores)), "--shm-size=8g",
        "-v", f"{job.run_dir}:/work", "-v", f"{reference}:/references:ro",
        image, *arguments,
    ]
    write_json(job.run_dir / "command.json", {"mode": "docker", "image": image, "command": command})
    update_status(job.run_dir, "running")
    try:
        completed = subprocess.run(command, text=True, capture_output=True, check=False)
        return_code, stdout, stderr = int(completed.returncode), completed.stdout or "", completed.stderr or ""
    except OSError as exc:
        return_code, stdout, stderr = -1, "", str(exc)
    (job.run_dir / "stdout.log").write_text(stdout, encoding="utf-8")
    (job.run_dir / "stderr.log").write_text(stderr, encoding="utf-8")
    generated_rows: list[dict[str, Any]] = []
    generated_table = output_root / "generated_compounds.csv"
    if generated_table.is_file():
        with generated_table.open(newline="", encoding="utf-8", errors="replace") as handle:
            for row in csv.DictReader(handle):
                smiles = str(row.get("canonical_isomeric_smiles") or row.get("smiles") or "").strip()
                if not smiles:
                    continue
                generated_rows.append({
                    "compound_id": row.get("compound_id") or "",
                    "smiles": smiles,
                    "properties": {
                        key: value for key, value in row.items()
                        if key not in {"compound_id", "canonical_isomeric_smiles", "smiles"} and value not in (None, "")
                    },
                })
    output_set = None
    if generated_rows:
        output_set = ligand_compounds.write_compound_set(
            job.run_dir / "artifacts" / "generated_compound_set.json",
            name=f"{engine.upper()} generated compounds",
            compounds=generated_rows,
            owner_id="",
            project_id="",
            created_at=str(read_json(job.run_dir / "metadata.json").get("created_at") or ""),
            compound_set_id=f"generated-{job.run_id}",
        )
        artifacts_root = job.run_dir / "artifacts"
        shutil.copy2(generated_table, artifacts_root / "generated_compounds.csv")
        generated_sdf = output_root / "generated_compounds.sdf"
        if generated_sdf.is_file():
            shutil.copy2(generated_sdf, artifacts_root / "generated_compounds.sdf")
    contract_path = job.run_dir / "artifact_contract.json"
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    contract["compound_count"] = len(output_set["compounds"]) if output_set else 0
    contract_path.write_text(json.dumps(contract, indent=2) + "\n", encoding="utf-8")
    success = return_code == 0 and bool(output_set)
    finish_job(job.run_dir, success, {
        "outputs": {
            "generated_compounds": "artifacts/generated_compounds.csv",
            "generated_compound_set": "artifacts/generated_compound_set.json",
            "generated_coordinates": "artifacts/generated_compounds.sdf",
            "generation_report": "normalized/generation_report.json",
        },
        "metrics": {
            "engine": "LDDM" if engine == "lddm" else "PocketXMol",
            "requested_compound_count": int(requested_count),
            "generated_compound_count": len(output_set["compounds"]) if output_set else 0,
            "return_code": return_code,
        },
        "downstream_artifacts": {
            "generated_compound_set": "artifacts/generated_compound_set.json",
            "generated_compounds": "artifacts/generated_compounds.csv",
            "generated_coordinates": "artifacts/generated_compounds.sdf",
        },
    })
    return job.run_dir


def run_ligand_lddm(
    target_pdb: Path,
    reference_ligand_sdf: Path,
    requested_count: int = 64,
    batch_size: int = 8,
    seed: int = 42,
    n_steps: int = 100,
    sampler: str = "ForwardEuler",
    sampling_noise: float = 5.0,
    gpu_device: str = "0",
    cpu_cores: int = 8,
) -> Path:
    """Generate pocket-conditioned ligands with the LDDM CD+BB checkpoint."""
    return _run_ligand_generation(
        engine="lddm", workflow="ligand.generation.lddm",
        target_pdb=target_pdb, pocket_pdb=None, reference_ligand_sdf=reference_ligand_sdf,
        requested_count=requested_count, batch_size=batch_size, seed=seed,
        gpu_device=gpu_device, cpu_cores=cpu_cores, n_steps=n_steps,
        sampler=sampler, sampling_noise=sampling_noise,
    )


def run_ligand_pocketxmol(
    target_pdb: Path,
    pocket_pdb: Path,
    requested_count: int = 64,
    batch_size: int = 8,
    seed: int = 42,
    diffusion_steps: int = 100,
    ligand_atoms_mean: float = 28.0,
    ligand_atoms_std: float = 2.0,
    gpu_device: str = "0",
    cpu_cores: int = 8,
) -> Path:
    """Generate pocket-conditioned ligands with the PocketXMol SBDD model."""
    return _run_ligand_generation(
        engine="pocketxmol", workflow="ligand.generation.pocketxmol",
        target_pdb=target_pdb, pocket_pdb=pocket_pdb, reference_ligand_sdf=None,
        requested_count=requested_count, batch_size=batch_size, seed=seed,
        gpu_device=gpu_device, cpu_cores=cpu_cores, diffusion_steps=diffusion_steps,
        ligand_atoms_mean=ligand_atoms_mean, ligand_atoms_std=ligand_atoms_std,
    )


def _run_ligand_generator(
    *, engine: str, target_pdb: Path, pocket_pdb: Path | None,
    reference_ligand_sdf: Path | None, pharmacophore_posp: Path | None,
    requested_count: int, batch_size: int, seed: int, gpu_device: str,
    cpu_cores: int, engine_settings: dict[str, Any] | None = None,
    pharmacophore_json: Path | None = None,
) -> Path:
    """Run a Studio-managed ligand-generation image as a standalone tile."""
    workflow = f"ligand.generation.{engine}"
    parameters: dict[str, Any] = {
        "target_pdb": target_pdb, "requested_count": requested_count,
        "batch_size": batch_size, "seed": seed, "gpu_device": gpu_device,
        "cpu_cores": cpu_cores, "engine_settings": engine_settings or {},
    }
    if pocket_pdb is not None:
        parameters["pocket_pdb"] = pocket_pdb
    if reference_ligand_sdf is not None:
        parameters["reference_ligand_sdf"] = reference_ligand_sdf
    if pharmacophore_posp is not None:
        parameters["pharmacophore_posp"] = pharmacophore_posp
    if pharmacophore_json is not None:
        parameters["pharmacophore_json"] = pharmacophore_json
    _validate_generator_tile_request(workflow, parameters, require_files=True)
    config = _GENERATOR_CONFIGS[engine]
    checkpoint = config["references"][0]
    reference = reference_root().expanduser().resolve()
    manifest = load_manifest(engine)
    params = {
        "engine": engine, "requested_compound_count": int(requested_count),
        "batch_size": int(batch_size), "seed": int(seed),
        "gpu_device": str(gpu_device), "cpu_cores": int(cpu_cores),
        "checkpoint": checkpoint, "engine_settings": engine_settings or {},
    }
    input_payload = {key: str(value) for key, value in parameters.items() if key in {"target_pdb", "pocket_pdb", "reference_ligand_sdf", "pharmacophore_posp", "pharmacophore_json"} and value is not None}
    job = create_job("ligand-generation", "ligand_generation", config["tool"], input_payload, params)
    input_root = job.run_dir / "artifacts" / "input"
    input_root.mkdir(parents=True, exist_ok=True)
    staged: dict[str, Path] = {}
    staged_names = {
        "target_pdb": (target_pdb, "target.pdb"),
        "pocket_pdb": (pocket_pdb, "pocket.pdb"),
        "reference_ligand_sdf": (reference_ligand_sdf, "reference_ligand.sdf"),
        "pharmacophore_posp": (pharmacophore_posp, "pharmacophore.posp"),
        "pharmacophore_json": (pharmacophore_json, "pharmit_pharmacophore.json"),
    }
    for key, (source, name) in staged_names.items():
        if source is None:
            continue
        destination = input_root / name
        shutil.copy2(Path(source), destination)
        staged[key] = destination
    input_contract_types = {
        "target_pdb": ("target_structure", "protein_structure/pdb/v1"),
        "pocket_pdb": ("pocket", "protein_pocket/pdb/v1"),
        "reference_ligand_sdf": ("reference_ligand", "ligand_coordinates/sdf/v1"),
        "pharmacophore_posp": ("pharmacophore", "pharmacophore/pgmg_posp/v1"),
        "pharmacophore_json": ("pharmacophore", "pharmacophore/pharmit_json/v1"),
    }
    contract_inputs = [
        {"name": input_contract_types[key][0], "type": input_contract_types[key][1], "path": path.relative_to(job.run_dir).as_posix()}
        for key, path in staged.items()
    ]
    _write_generation_contract(job.run_dir, workflow=workflow, engine=engine, inputs=contract_inputs, checkpoint=checkpoint)
    native_root = job.run_dir / "native"
    output_root = job.run_dir / "normalized"
    native_root.mkdir(parents=True, exist_ok=True)
    output_root.mkdir(parents=True, exist_ok=True)
    work_path = lambda key: f"/work/{staged[key].relative_to(job.run_dir).as_posix()}"
    settings = dict(engine_settings or {})
    native_mode = "fixed_protein_ligand_denovo_condensed"
    if engine == "omtra" and "pharmacophore_json" in staged:
        native_mode = "fixed_protein_pharmacophore_ligand_denovo_condensed"
    args = [
        "--output", "/work/native", "--normalized-output", "/work/normalized",
        "--count", str(int(requested_count)), "--batch-size", str(int(batch_size)),
        "--seed", str(int(seed)), "--run-name", f"{engine}_{job.run_id[-8:]}",
        "--checkpoint", f"/references/{checkpoint}", "--native-mode", native_mode,
    ]
    if engine == "omtra":
        args.extend([
            "--integration-steps", str(max(25, int(settings.get("integration_steps", 250)))),
            "--noise-scale", str(float(settings.get("noise_scale", 1.0))),
            "--epsilon", str(float(settings.get("epsilon", 0.01))),
            "--ligand-atoms-mean", str(float(settings.get("ligand_atoms_mean", 0.0))),
            "--ligand-atoms-std", str(float(settings.get("ligand_atoms_std", 2.0))),
        ])
        if settings.get("stochastic_sampling"):
            args.append("--stochastic-sampling")
        args.extend(["--target", work_path("target_pdb")])
        for key, flag in (("pocket_pdb", "--pocket-structure"), ("reference_ligand_sdf", "--reference-ligand")):
            if key in staged:
                args.extend([flag, work_path(key)])
        if "pharmacophore_json" in staged:
            args.extend(["--pharmacophore", work_path("pharmacophore_json")])
    elif engine in {"pocketflow", "pfm", "drugrpg"}:
        if engine == "pocketflow":
            args.extend([
                "--atom-temperature", str(float(settings.get("atom_temperature", 1.0))),
                "--bond-temperature", str(float(settings.get("bond_temperature", 1.0))),
                "--max-atoms", str(max(5, int(settings.get("max_atoms", 40)))),
                "--focus-strategy", str(settings.get("focus_strategy", "maximum")),
                "--focus-threshold", str(float(settings.get("focus_threshold", 0.5))),
                "--min-protein-distance", str(float(settings.get("min_protein_distance", 3.0))),
            ])
        elif engine == "drugrpg":
            args.extend(["--max-atoms", str(max(5, int(settings.get("max_atoms", 30))))])
        args.extend(["--pocket-structure", work_path("pocket_pdb")])
    elif engine == "lddm":
        args.extend(["--target", work_path("target_pdb"), "--reference-ligand", work_path("reference_ligand_sdf")])
    elif engine in {"conditar", "paopt"}:
        args.extend(["--pocket-structure", work_path("pocket_pdb"), "--secondary-reference", f"/references/{config['references'][1]}"])
        args.extend(["--diffusion-steps", str(max(100, int(settings.get("diffusion_steps", 1000)))), "--pocket-radius", str(float(settings.get("pocket_radius", 10.0)))])
        if engine == "paopt":
            args.extend(["--reference-ligand", work_path("reference_ligand_sdf")])
            args.extend(["--optimization-steps", str(max(1, int(settings.get("optimization_steps", 1))))])
            args.extend(["--gradient-estimate-pairs", str(max(1, int(settings.get("gradient_estimate_pairs", 4))))])
            args.extend(["--perturbation-size", str(float(settings.get("perturbation_size", 0.03)))])
            for value in settings.get("optimize_properties") or []:
                args.extend(["--optimize-properties", str(value)])
            for value in settings.get("minimize_properties") or []:
                args.extend(["--minimize-properties", str(value)])
    elif engine == "flowr_root":
        args.extend([
            "--integration-steps", str(max(10, int(settings.get("integration_steps", 100)))),
            "--corrector-steps", str(max(0, int(settings.get("corrector_steps", 0)))),
            "--solver", str(settings.get("solver", "euler")),
            "--diversity-threshold", str(float(settings.get("diversity_threshold", 0.9))),
            "--pocket-structure", work_path("pocket_pdb"),
            "--reference-ligand", work_path("reference_ligand_sdf"),
        ])
        for setting, flag in (("use_sde_simulation", "--use-sde-simulation"), ("sample_molecule_sizes", "--sample-molecule-sizes"), ("filter_diversity", "--filter-diversity")):
            if settings.get(setting):
                args.append(flag)
    elif engine == "pgmg":
        args.extend(["--pharmacophore", work_path("pharmacophore_posp"), "--secondary-reference", f"/references/{config['references'][1]}"])
    command = [
        "docker", "run", "--rm", *docker_gpu_args(gpu_device), "--cpus", str(int(cpu_cores)),
        "--shm-size=8g", "-v", f"{job.run_dir}:/work", "-v", f"{reference}:/references:ro",
        str(manifest["image"]), *args,
    ]
    write_json(job.run_dir / "command.json", {"mode": "docker", "image": manifest["image"], "command": command})
    update_status(job.run_dir, "running")
    try:
        completed = subprocess.run(command, text=True, capture_output=True, check=False)
        return_code, stdout, stderr = int(completed.returncode), completed.stdout or "", completed.stderr or ""
    except OSError as exc:
        return_code, stdout, stderr = -1, "", str(exc)
    (job.run_dir / "stdout.log").write_text(stdout, encoding="utf-8")
    (job.run_dir / "stderr.log").write_text(stderr, encoding="utf-8")
    generated_rows: list[dict[str, Any]] = []
    generated_table = output_root / "generated_compounds.csv"
    if generated_table.is_file():
        with generated_table.open(newline="", encoding="utf-8", errors="replace") as handle:
            for row in csv.DictReader(handle):
                smiles = str(row.get("canonical_isomeric_smiles") or row.get("smiles") or "").strip()
                if smiles:
                    generated_rows.append({
                        "compound_id": row.get("compound_id") or "", "smiles": smiles,
                        "properties": {key: value for key, value in row.items() if key not in {"compound_id", "canonical_isomeric_smiles", "smiles"} and value not in (None, "")},
                    })
    output_set = None
    if generated_rows:
        output_set = ligand_compounds.write_compound_set(
            job.run_dir / "artifacts" / "generated_compound_set.json",
            name=f"{config['name']} generated compounds", compounds=generated_rows,
            owner_id="", project_id="", created_at=str(read_json(job.run_dir / "metadata.json").get("created_at") or ""),
            compound_set_id=f"generated-{job.run_id}",
        )
        artifacts_root = job.run_dir / "artifacts"
        shutil.copy2(generated_table, artifacts_root / "generated_compounds.csv")
        generated_sdf = output_root / "generated_compounds.sdf"
        if generated_sdf.is_file():
            shutil.copy2(generated_sdf, artifacts_root / "generated_compounds.sdf")
    contract_path = job.run_dir / "artifact_contract.json"
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    contract["compound_count"] = len(output_set["compounds"]) if output_set else 0
    contract_path.write_text(json.dumps(contract, indent=2) + "\n", encoding="utf-8")
    finish_job(job.run_dir, return_code == 0 and bool(output_set), {
        "outputs": {
            "generated_compounds": "artifacts/generated_compounds.csv",
            "generated_compound_set": "artifacts/generated_compound_set.json",
            "generated_coordinates": "artifacts/generated_compounds.sdf",
            "generation_report": "normalized/generation_report.json",
        },
        "metrics": {"engine": config["name"], "requested_compound_count": int(requested_count), "generated_compound_count": len(output_set["compounds"]) if output_set else 0, "return_code": return_code},
        "downstream_artifacts": {"generated_compound_set": "artifacts/generated_compound_set.json", "generated_compounds": "artifacts/generated_compounds.csv", "generated_coordinates": "artifacts/generated_compounds.sdf"},
    })
    return job.run_dir


def _generator_tile_wrapper(engine: str, target_pdb: Path, pocket_pdb: Path | None = None, reference_ligand_sdf: Path | None = None, pharmacophore_posp: Path | None = None, requested_count: int = 64, batch_size: int = 8, seed: int = 42, gpu_device: str = "0", cpu_cores: int = 8, engine_settings: dict[str, Any] | None = None, pharmacophore_json: Path | None = None) -> Path:
    return _run_ligand_generator(
        engine=engine, target_pdb=target_pdb, pocket_pdb=pocket_pdb,
        reference_ligand_sdf=reference_ligand_sdf, pharmacophore_posp=pharmacophore_posp,
        requested_count=requested_count, batch_size=batch_size, seed=seed,
        gpu_device=gpu_device, cpu_cores=cpu_cores, engine_settings=engine_settings,
        pharmacophore_json=pharmacophore_json,
    )


def run_ligand_omtra(target_pdb: Path, pocket_pdb: Path | None = None, reference_ligand_sdf: Path | None = None, pharmacophore_posp: Path | None = None, requested_count: int = 64, batch_size: int = 8, seed: int = 42, gpu_device: str = "0", cpu_cores: int = 8, engine_settings: dict[str, Any] | None = None, pharmacophore_json: Path | None = None) -> Path:
    return _generator_tile_wrapper("omtra", target_pdb, pocket_pdb, reference_ligand_sdf, pharmacophore_posp, requested_count, batch_size, seed, gpu_device, cpu_cores, engine_settings, pharmacophore_json)


def run_ligand_flowr_root(target_pdb: Path, pocket_pdb: Path | None = None, reference_ligand_sdf: Path | None = None, pharmacophore_posp: Path | None = None, requested_count: int = 64, batch_size: int = 8, seed: int = 42, gpu_device: str = "0", cpu_cores: int = 8, engine_settings: dict[str, Any] | None = None) -> Path:
    return _generator_tile_wrapper("flowr_root", target_pdb, pocket_pdb, reference_ligand_sdf, pharmacophore_posp, requested_count, batch_size, seed, gpu_device, cpu_cores, engine_settings)


def run_ligand_conditar(target_pdb: Path, pocket_pdb: Path | None = None, reference_ligand_sdf: Path | None = None, pharmacophore_posp: Path | None = None, requested_count: int = 64, batch_size: int = 8, seed: int = 42, gpu_device: str = "0", cpu_cores: int = 8, engine_settings: dict[str, Any] | None = None) -> Path:
    return _generator_tile_wrapper("conditar", target_pdb, pocket_pdb, reference_ligand_sdf, pharmacophore_posp, requested_count, batch_size, seed, gpu_device, cpu_cores, engine_settings)


def run_ligand_paopt(target_pdb: Path, pocket_pdb: Path | None = None, reference_ligand_sdf: Path | None = None, pharmacophore_posp: Path | None = None, requested_count: int = 64, batch_size: int = 8, seed: int = 42, gpu_device: str = "0", cpu_cores: int = 8, engine_settings: dict[str, Any] | None = None) -> Path:
    return _generator_tile_wrapper("paopt", target_pdb, pocket_pdb, reference_ligand_sdf, pharmacophore_posp, requested_count, batch_size, seed, gpu_device, cpu_cores, engine_settings)


def run_ligand_drugrpg(target_pdb: Path, pocket_pdb: Path | None = None, reference_ligand_sdf: Path | None = None, pharmacophore_posp: Path | None = None, requested_count: int = 64, batch_size: int = 8, seed: int = 42, gpu_device: str = "0", cpu_cores: int = 8, engine_settings: dict[str, Any] | None = None) -> Path:
    return _generator_tile_wrapper("drugrpg", target_pdb, pocket_pdb, reference_ligand_sdf, pharmacophore_posp, requested_count, batch_size, seed, gpu_device, cpu_cores, engine_settings)


def run_ligand_pfm(target_pdb: Path, pocket_pdb: Path | None = None, reference_ligand_sdf: Path | None = None, pharmacophore_posp: Path | None = None, requested_count: int = 64, batch_size: int = 8, seed: int = 42, gpu_device: str = "0", cpu_cores: int = 8, engine_settings: dict[str, Any] | None = None) -> Path:
    return _generator_tile_wrapper("pfm", target_pdb, pocket_pdb, reference_ligand_sdf, pharmacophore_posp, requested_count, batch_size, seed, gpu_device, cpu_cores, engine_settings)


def run_ligand_pocketflow(target_pdb: Path, pocket_pdb: Path | None = None, reference_ligand_sdf: Path | None = None, pharmacophore_posp: Path | None = None, requested_count: int = 64, batch_size: int = 8, seed: int = 42, gpu_device: str = "0", cpu_cores: int = 8, engine_settings: dict[str, Any] | None = None) -> Path:
    return _generator_tile_wrapper("pocketflow", target_pdb, pocket_pdb, reference_ligand_sdf, pharmacophore_posp, requested_count, batch_size, seed, gpu_device, cpu_cores, engine_settings)


def run_ligand_pgmg(target_pdb: Path, pocket_pdb: Path | None = None, reference_ligand_sdf: Path | None = None, pharmacophore_posp: Path | None = None, requested_count: int = 64, batch_size: int = 8, seed: int = 42, gpu_device: str = "0", cpu_cores: int = 8, engine_settings: dict[str, Any] | None = None) -> Path:
    return _generator_tile_wrapper("pgmg", target_pdb, pocket_pdb, reference_ligand_sdf, pharmacophore_posp, requested_count, batch_size, seed, gpu_device, cpu_cores, engine_settings)


def run_ligand_nesso(
    target_pdb: Path,
    compound_set_file: Path,
    recycling_steps: int = 5,
    seed: int = 42,
    compound_limit: int = 0,
    gpu_device: str = "0",
) -> Path:
    """Predict affinity for a Studio Database compound set with Nesso-1."""
    validate_target_pdb(target_pdb)
    if not 1 <= recycling_steps <= 64 or not 1 <= seed <= 2_147_483_646:
        raise ValueError("Nesso recycling steps or seed is outside the supported range.")
    payload, compounds = _validated_compound_set(compound_set_file, compound_limit)
    readiness = nesso_readiness()
    if not readiness["ready"]:
        missing = [key for key in ("checkpoint_ready", "ccd_ready", "esm_ready") if not readiness[key]]
        raise ValueError("Nesso references are incomplete (" + ", ".join(missing) + ").")
    proteins = _pdb_protein_sequences(target_pdb)
    params = {
        "recycling_steps": recycling_steps, "seed": seed, "compound_limit": compound_limit,
        "gpu_device": str(gpu_device), "compound_set_name": payload.get("name", "Compound set"),
        "compound_count": len(compounds), "structure_output": False,
    }
    job = create_job(
        "ligand-affinity", "ligand_affinity", "nesso_affinity",
        {"target_pdb": str(target_pdb), "compound_set_file": str(compound_set_file)}, params,
    )
    _write_job_contract(
        job.run_dir, "ligand.affinity.nesso", "compound_set/v1", "ligand_affinity_result/v1",
        payload, compounds, target_artifact="artifacts/input/target.pdb",
        outputs=[
            {"name": "affinity_metrics", "type": "ligand_affinity_scores/v1", "path": "artifacts/nesso_affinity.csv"},
            {"name": "native_predictions", "type": "nesso_affinity_predictions/v1", "path": "output"},
        ],
    )
    input_root = job.run_dir / "artifacts" / "input"
    input_root.mkdir(parents=True, exist_ok=True)
    (input_root / "target.pdb").write_bytes(Path(target_pdb).read_bytes())
    _write_compound_inputs(job.run_dir, payload, compounds)
    nesso_inputs = job.run_dir / "inputs"
    nesso_inputs.mkdir(parents=True, exist_ok=True)
    for index, compound in enumerate(compounds, start=1):
        compound_id = _safe_id(compound.get("compound_id"), index)
        (nesso_inputs / f"{compound_id}.yaml").write_text(
            _nesso_yaml(proteins, str(compound["smiles"])), encoding="utf-8"
        )
    image = str(load_manifest("nesso")["image"])
    reference = reference_root().expanduser().resolve()
    checkpoint = Path(readiness["checkpoint_dir"]).resolve()
    ccd = Path(readiness["ccd_path"]).resolve()
    cache = Path(readiness["esm_cache_dir"]).resolve()
    try:
        relative_checkpoint = checkpoint.relative_to(reference)
        relative_ccd = ccd.relative_to(reference)
        relative_cache = cache.relative_to(reference)
    except ValueError as exc:
        raise ValueError("Nesso checkpoint, CCD, and ESM cache must be inside the configured Studio reference directory.") from exc
    native_args = [
        "predict", "/work/inputs", "--out_dir", "/work/output", "--cache", "/cache",
        "--checkpoint", f"/ref/{relative_checkpoint.as_posix()}", "--ccd", f"/ref/{relative_ccd.as_posix()}",
        "--accelerator", "gpu", "--devices", "1", "--precision", "bf16-mixed",
        "--recycling_steps", str(recycling_steps), "--num_workers", "2", "--require_affinity", "--override",
        "--seed", str(seed), "--no_kernels",
    ]
    command = [
        "docker", "run", "--rm", *docker_gpu_args(gpu_device),
        "-v", f"{job.run_dir}:/work", "-v", f"{reference}:/ref:ro",
        "-v", f"{cache}:/cache/huggingface:ro",
        "--env", "HF_HUB_OFFLINE=1", "--env", "TRANSFORMERS_OFFLINE=1", "--shm-size", "8g",
        image, *native_args,
    ]
    write_json(job.run_dir / "command.json", {"mode": "docker", "command": command, "image": image})
    update_status(job.run_dir, "running")
    try:
        completed = subprocess.run(command, text=True, capture_output=True, check=False)
        return_code = int(completed.returncode)
        stdout, stderr = completed.stdout or "", completed.stderr or ""
    except OSError as exc:
        return_code, stdout, stderr = -1, "", str(exc)
    (job.run_dir / "stdout.log").write_text(stdout, encoding="utf-8")
    (job.run_dir / "stderr.log").write_text(stderr, encoding="utf-8")
    rows: list[dict[str, Any]] = []
    for affinity_path in sorted((job.run_dir / "output").glob("**/predictions/*/affinity.json")):
        try:
            native = json.loads(affinity_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        if not isinstance(native, dict):
            continue
        affinity = _finite_float(native, "affinity_pred_value")
        if affinity is None:
            continue
        rows.append({
            "compound_id": affinity_path.parent.name,
            "affinity_log10_ic50_uM": affinity,
            "pIC50": 6.0 - affinity,
            "ic50_uM": 10.0 ** affinity,
            "binder_probability": _finite_float(native, "affinity_probability_binary"),
            "native_file": str(affinity_path.relative_to(job.run_dir)),
        })
    metrics_path = job.run_dir / "artifacts" / "nesso_affinity.csv"
    _write_csv(metrics_path, rows, ["compound_id", "affinity_log10_ic50_uM", "pIC50", "ic50_uM", "binder_probability", "native_file"])
    success = return_code == 0 and len({row["compound_id"] for row in rows}) == len(compounds)
    finish_job(job.run_dir, success, {
        "outputs": {"affinity_metrics": "artifacts/nesso_affinity.csv", "structure_output": False},
        "metrics": {"engine": "Nesso-1", "compound_count": len(compounds), "scored_count": len(rows), "return_code": return_code},
        "downstream_artifacts": {"affinity_metrics": "artifacts/nesso_affinity.csv", "affinity_results": [row["native_file"] for row in rows]},
    })
    return job.run_dir


def run_ligand_gromacs_md(
    complex_pdb: Path,
    ligand_sdf: Path,
    ligand_key: str,
    production_ns: float = 1.0,
    timestep_fs: float = 2.0,
    temperature: float = 300.0,
    pressure: float = 1.0,
    padding_nm: float = 1.0,
    ionic_strength: float = 0.15,
    forcefield_method: str = "gaff2",
    protein_forcefield_method: str = "ff14SB",
    water_model: str = "tip3p",
    hydrogen_mass_repartitioning: bool = True,
    use_gpu: bool = True,
    gpu_device: str = "0",
    cpu_cores: int = 8,
    seed: int = 42,
) -> Path:
    """Run bound protein–ligand preparation and production MD with GROMACS."""
    params = {
        "ligand_key": str(ligand_key),
        "production_ns": float(production_ns),
        "timestep_fs": float(timestep_fs),
        "temperature": float(temperature),
        "pressure": float(pressure),
        "padding_nm": float(padding_nm),
        "ionic_strength": float(ionic_strength),
        "forcefield_method": str(forcefield_method),
        "protein_forcefield_method": str(protein_forcefield_method),
        "water_model": str(water_model),
        "hydrogen_mass_repartitioning": bool(hydrogen_mass_repartitioning),
        "use_gpu": bool(use_gpu),
        "gpu_device": str(gpu_device if use_gpu else "none"),
        "cpu_cores": int(cpu_cores),
        "seed": int(seed),
    }
    _validate_gromacs_md_request(
        {**params, "complex_pdb": str(complex_pdb), "ligand_sdf": str(ligand_sdf)},
        require_files=True,
    )
    readiness = gromacs_md_runner_readiness()
    runtime_root = Path(readiness["source"]).resolve()
    job = create_job(
        "ligand-md", "ligand_md", "gromacs_md",
        {"complex_pdb": str(complex_pdb), "ligand_sdf": str(ligand_sdf)},
        params,
    )
    artifacts = job.run_dir / "artifacts"
    input_root = artifacts / "input"
    input_root.mkdir(parents=True, exist_ok=True)
    staged_complex = input_root / "complex.pdb"
    staged_ligand = input_root / "ligand.sdf"
    shutil.copy2(complex_pdb, staged_complex)
    shutil.copy2(ligand_sdf, staged_ligand)
    md_root = artifacts / "md"
    md_root.mkdir(parents=True, exist_ok=True)
    production_steps = max(1, int(round(float(production_ns) * 1_000_000.0 / float(timestep_fs))))
    config = {
        "pdb_id": "STUDIO",
        "system_id": f"studio_{job.run_id[-8:]}",
        "job_id": job.run_id,
        "ligand_key": str(ligand_key),
        "input_complex_pdb_path": "/work/artifacts/input/complex.pdb",
        "prepared_complex_path": "/work/artifacts/input/complex.pdb",
        "ligand_refined_sdf_path": "/work/artifacts/input/ligand.sdf",
        "ligand_refined_sdf_data": staged_ligand.read_text(encoding="utf-8", errors="replace"),
        "strict_refined_ligand": True,
        "forcefield_method": str(forcefield_method),
        "protein_forcefield_method": str(protein_forcefield_method),
        "water_model": str(water_model),
        "box_shape": "dodecahedron",
        "padding_nm": float(padding_nm),
        "ionic_strength": float(ionic_strength),
        "mmgbsa_backend": "ambertools_mmpbsa",
        "temperature": float(temperature),
        "pressure": float(pressure),
        "production_timestep_fs": float(timestep_fs),
        "hydrogen_mass_amu": 4.0 if timestep_fs == 4.0 and hydrogen_mass_repartitioning else None,
        "mass_repartition_factor": 3.0 if timestep_fs == 4.0 and hydrogen_mass_repartitioning else None,
        "production_steps": production_steps,
        "production_length_ns": float(production_ns),
        "production_report_interval": max(1, int(round(5000 * float(timestep_fs) / 2.0))),
        "replica_seed": int(seed),
        "use_gpu": bool(use_gpu),
        "density_stabilization_min_ns": 0.25,
        "density_stabilization_max_ns": 1.0,
        "density_stabilization_increment_ns": 0.25,
        "density_plateau_required": False,
        "density_sample_interval_ps": 4.0,
        **params,
    }
    # The downstream runner uses production_timestep_fs rather than the UI key.
    config["production_timestep_fs"] = float(timestep_fs)
    config_path = input_root / "gromacs_config.json"
    write_json(config_path, config)
    (input_root / "gromacs_tile_runner.py").write_text(_GROMACS_TILE_RUNNER, encoding="utf-8")
    metadata_path = job.run_dir / "metadata.json"
    metadata = read_json(metadata_path)
    metadata.update({
        "workflow": "ligand.md.gromacs",
        "engine": "GROMACS 2026.3",
        "result_type": "ligand_md_result/v1",
        "artifact_contract_version": 1,
        "selected_ligand_key": str(ligand_key),
        "production_ns": float(production_ns),
    })
    write_json(metadata_path, metadata)
    write_json(job.run_dir / "artifact_contract.json", {
        "schema_version": 1,
        "workflow": "ligand.md.gromacs",
        "input_type": "protein_ligand_complex/v1",
        "output_type": "ligand_md_result/v1",
        "inputs": [
            {"name": "protein_ligand_complex", "type": "protein_ligand_complex/pdb/v1", "path": "artifacts/input/complex.pdb"},
            {"name": "ligand_coordinates", "type": "ligand_coordinates/sdf/v1", "path": "artifacts/input/ligand.sdf", "selected_residue_key": str(ligand_key)},
        ],
        "outputs": [
            {"name": "trajectory_set", "type": "md_trajectory_set/v1", "path": "artifacts/md/production"},
            {"name": "trajectory_analysis", "type": "md_analysis/v1", "path": "artifacts/md/production/trajectory_analysis.json"},
            {"name": "prepared_system", "type": "md_prepared_system/v1", "path": "artifacts/md/preparation"},
        ],
        "structure_output": True,
        "compound_set_id": "",
        "compound_count": 1,
    })
    manifest = load_manifest("gromacs_md")
    selected_gpu = str(gpu_device if use_gpu else "none")
    return_code = run_docker_job(
        job.run_dir,
        manifest,
        ["python", "/work/artifacts/input/gromacs_tile_runner.py"],
        gpu_device=selected_gpu,
        extra_mounts=[(runtime_root, "/studio-ligand-runtime", True)],
        environment={"PYTHONPATH": "/studio-ligand-runtime"},
    )
    result_path = md_root / "result.json"
    result = read_json(result_path) if result_path.is_file() else {}
    output_files = result.get("output_files") or {}
    def _resolve_container_output(value: object) -> Path:
        path = Path(str(value or ""))
        if path.is_absolute() and path.is_relative_to(Path("/work")):
            return job.run_dir / path.relative_to("/work")
        return path

    analysis_path = _resolve_container_output(output_files.get("analysis_report"))
    trajectory_path = _resolve_container_output(output_files.get("production_trajectory"))
    success = bool(
        return_code == 0
        and result.get("success") is True
        and analysis_path.is_file()
        and trajectory_path.is_file()
    )
    if not success and not result.get("error"):
        stderr_path = job.run_dir / "stderr.log"
        result["error"] = stderr_path.read_text(encoding="utf-8", errors="replace")[-4000:] if stderr_path.is_file() else "GROMACS container did not produce a valid trajectory and analysis report."
    finish_job(job.run_dir, success, {
        "outputs": {
            "trajectory_set": "artifacts/md/production",
            "trajectory_analysis": "artifacts/md/production/trajectory_analysis.json",
            "prepared_system": "artifacts/md/preparation",
            "run_summary": "artifacts/md/result.json",
        },
        "metrics": {
            "engine": "GROMACS 2026.3",
            "ligand_key": str(ligand_key),
            "production_ns": float(production_ns),
            "timestep_fs": float(timestep_fs),
            "return_code": return_code,
            "error": result.get("error", ""),
        },
        "downstream_artifacts": {
            "md_trajectory_set": "artifacts/md/production",
            "md_analysis": "artifacts/md/production/trajectory_analysis.json",
            "prepared_md_system": "artifacts/md/preparation",
        },
    })
    return job.run_dir


def run_ligand_openmm_md(
    complex_pdb: Path,
    ligand_sdf: Path,
    ligand_key: str,
    production_ns: float = 1.0,
    timestep_fs: float = 2.0,
    temperature: float = 300.0,
    pressure: float = 1.0,
    padding_nm: float = 1.0,
    ionic_strength: float = 0.15,
    forcefield_method: str = "openff-2.2.0",
    protein_forcefield_method: str = "amber14-all",
    water_model: str = "tip3p",
    hydrogen_mass_repartitioning: bool = True,
    use_gpu: bool = True,
    gpu_device: str = "0",
    cpu_cores: int = 8,
    seed: int = 42,
) -> Path:
    """Run bound protein–ligand system preparation and MD with OpenMM."""
    params = {
        "ligand_key": str(ligand_key), "production_ns": float(production_ns),
        "timestep_fs": float(timestep_fs), "temperature": float(temperature),
        "pressure": float(pressure), "padding_nm": float(padding_nm),
        "ionic_strength": float(ionic_strength), "forcefield_method": str(forcefield_method),
        "protein_forcefield_method": str(protein_forcefield_method), "water_model": str(water_model),
        "hydrogen_mass_repartitioning": bool(hydrogen_mass_repartitioning),
        "use_gpu": bool(use_gpu), "gpu_device": str(gpu_device if use_gpu else "none"),
        "cpu_cores": int(cpu_cores), "seed": int(seed),
    }
    _validate_gromacs_md_request(
        {**params, "complex_pdb": str(complex_pdb), "ligand_sdf": str(ligand_sdf)},
        require_files=True, engine="openmm",
    )
    readiness = openmm_md_runner_readiness()
    runtime_root = Path(readiness["source"]).resolve()
    job = create_job(
        "ligand-md", "ligand_md", "openmm_md",
        {"complex_pdb": str(complex_pdb), "ligand_sdf": str(ligand_sdf)}, params,
    )
    artifacts = job.run_dir / "artifacts"
    input_root = artifacts / "input"
    input_root.mkdir(parents=True, exist_ok=True)
    staged_complex = input_root / "complex.pdb"
    staged_ligand = input_root / "ligand.sdf"
    shutil.copy2(complex_pdb, staged_complex)
    shutil.copy2(ligand_sdf, staged_ligand)
    md_root = artifacts / "md"
    md_root.mkdir(parents=True, exist_ok=True)
    production_steps = max(1, int(round(float(production_ns) * 1_000_000.0 / float(timestep_fs))))
    config = {
        "pdb_id": "STUDIO", "system_id": f"studio_{job.run_id[-8:]}", "job_id": job.run_id,
        "ligand_key": str(ligand_key), "input_complex_pdb_path": "/work/artifacts/input/complex.pdb",
        "prepared_complex_path": "/work/artifacts/input/complex.pdb",
        "ligand_refined_sdf_path": "/work/artifacts/input/ligand.sdf",
        "ligand_refined_sdf_data": staged_ligand.read_text(encoding="utf-8", errors="replace"),
        "strict_refined_ligand": True, "forcefield_method": str(forcefield_method),
        "protein_forcefield_method": str(protein_forcefield_method), "water_model": str(water_model),
        "box_shape": "dodecahedron", "padding_nm": float(padding_nm),
        "ionic_strength": float(ionic_strength), "temperature": float(temperature),
        "pressure": float(pressure), "production_timestep_fs": float(timestep_fs),
        "hydrogen_mass_amu": 4.0 if timestep_fs == 4.0 and hydrogen_mass_repartitioning else None,
        "production_steps": production_steps, "production_length_ns": float(production_ns),
        "production_report_interval": max(1, int(round(5000 * float(timestep_fs) / 2.0))),
        "replica_seed": int(seed), "charge_method": "gasteiger", "generate_conformer": False,
        "preparation_protocol": "current_staged",
        "mmgbsa_enabled": False, "density_stabilization_min_ns": 0.25,
        "density_stabilization_max_ns": 1.0, "density_stabilization_increment_ns": 0.25,
        "density_plateau_required": False, "density_sample_interval_ps": 4.0,
        **params,
    }
    config["production_timestep_fs"] = float(timestep_fs)
    config["ligand_refined_sdf_data"] = staged_ligand.read_text(encoding="utf-8", errors="replace")
    config_path = input_root / "openmm_config.json"
    write_json(config_path, config)
    (input_root / "openmm_tile_runner.py").write_text(_OPENMM_TILE_RUNNER, encoding="utf-8")
    metadata_path = job.run_dir / "metadata.json"
    metadata = read_json(metadata_path)
    metadata.update({
        "workflow": "ligand.md.openmm", "engine": "OpenMM",
        "result_type": "ligand_md_result/v1", "artifact_contract_version": 1,
        "selected_ligand_key": str(ligand_key), "production_ns": float(production_ns),
    })
    write_json(metadata_path, metadata)
    write_json(job.run_dir / "artifact_contract.json", {
        "schema_version": 1, "workflow": "ligand.md.openmm",
        "input_type": "protein_ligand_complex/v1", "output_type": "ligand_md_result/v1",
        "inputs": [
            {"name": "protein_ligand_complex", "type": "protein_ligand_complex/pdb/v1", "path": "artifacts/input/complex.pdb"},
            {"name": "ligand_coordinates", "type": "ligand_coordinates/sdf/v1", "path": "artifacts/input/ligand.sdf", "selected_residue_key": str(ligand_key)},
        ],
        "outputs": [
            {"name": "trajectory_set", "type": "md_trajectory_set/v1", "path": "artifacts/md"},
            {"name": "run_summary", "type": "md_run_summary/v1", "path": "artifacts/md/result.json"},
        ],
        "structure_output": True, "compound_set_id": "", "compound_count": 1,
    })
    manifest = load_manifest("openmm_md")
    return_code = run_docker_job(
        job.run_dir, manifest,
        ["python", "/work/artifacts/input/openmm_tile_runner.py"],
        gpu_device=str(gpu_device if use_gpu else "none"),
        extra_mounts=[(runtime_root, "/studio-ligand-runtime", True)],
        environment={"PYTHONPATH": "/studio-ligand-runtime"},
    )
    result_path = md_root / "result.json"
    result = read_json(result_path) if result_path.is_file() else {}
    nested = result.get("md_result") or {}
    output_files = nested.get("output_files") or {}
    def _resolve_openmm_output(value: object) -> Path:
        path = Path(str(value or ""))
        if path.is_absolute() and path.is_relative_to(Path("/work")):
            return job.run_dir / path.relative_to("/work")
        return path
    trajectory = _resolve_openmm_output(output_files.get("production_trajectory"))
    success = bool(return_code == 0 and result.get("success") is True and trajectory.is_file())
    if not success and not result.get("error"):
        stderr_path = job.run_dir / "stderr.log"
        result["error"] = stderr_path.read_text(encoding="utf-8", errors="replace")[-4000:] if stderr_path.is_file() else "OpenMM did not produce a valid production trajectory."
    finish_job(job.run_dir, success, {
        "outputs": {"trajectory_set": "artifacts/md", "run_summary": "artifacts/md/result.json"},
        "metrics": {
            "engine": "OpenMM", "ligand_key": str(ligand_key),
            "production_ns": float(production_ns), "timestep_fs": float(timestep_fs),
            "return_code": return_code, "error": result.get("error", ""),
        },
        "downstream_artifacts": {"md_trajectory_set": "artifacts/md", "md_run_summary": "artifacts/md/result.json"},
    })
    return job.run_dir
