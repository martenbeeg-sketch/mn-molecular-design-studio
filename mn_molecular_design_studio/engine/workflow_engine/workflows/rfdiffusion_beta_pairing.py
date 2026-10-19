"""Paper-inspired RFdiffusion generation with explicit SS/adjacency tensors."""
from __future__ import annotations

import re
import secrets
import shutil
import subprocess
from pathlib import Path
from typing import Any

from mn_molecular_design_studio import rfdiffusion_beta_pairing_input as inputs
from mn_molecular_design_studio.engine.workflow_engine.core.candidates import write_candidates
from mn_molecular_design_studio.engine.workflow_engine.core.gpu import docker_gpu_args, normalize_gpu_device
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import create_job, update_status, write_json
from mn_molecular_design_studio.runtime import reference_root
from mn_molecular_design_studio.engine.workflow_engine.workflows import design


RF_IMAGE = "mn-rfdiffusion:latest"
TENSOR_IMAGE = "mn-protein-scoring-pyrosetta:latest"
RF_SOURCE_REVISION = "86507b6"


def _container_path(run_dir: Path, path: Path) -> str:
    return "/work/" + str(path.relative_to(run_dir))


def _normalize_target(source: Path, chain: str, selected_numbers: list[int], destination: Path) -> tuple[list[list[int]], dict[str, int]]:
    """Keep one chain and renumber it contiguously, retaining a reversible ID map."""
    if source.suffix.lower() not in {".pdb", ".ent"}:
        raise ValueError("RFdiffusion β-pairing currently needs a PDB target.")
    destination.parent.mkdir(parents=True, exist_ok=True)
    source_lines = source.read_text(errors="replace").splitlines()
    residues: list[tuple[int, str, str]] = []
    seen: set[tuple[int, str]] = set()
    for line in source_lines:
        if not line.startswith("ATOM  ") or line[21:22] != chain:
            continue
        resname = line[17:20].strip().upper()
        if resname not in inputs.AA3:
            continue
        key = (int(line[22:26]), line[26:27].strip())
        if key not in seen:
            residues.append((key[0], key[1], resname))
            seen.add(key)
    if not residues:
        raise ValueError(f"No standard protein atoms were found in target chain {chain}.")
    if len(residues) > 5000:
        raise ValueError("The selected target chain is longer than 5,000 residues; trim it before generation.")
    # The author script addresses target positions by chain and residue number.
    # Reject ambiguous duplicate sequence numbers and preserve explicit insertion-code residues only
    # after assigning each a unique normalized number.
    lookup: dict[tuple[int, str], int] = {
        (old_number, insertion): new_number
        for new_number, (old_number, insertion, _name) in enumerate(residues, start=1)
    }
    old_to_new: dict[str, int] = {}
    for old_number, insertion, _name in residues:
        label = f"{chain}{old_number}"
        if insertion:
            label += insertion
        old_to_new[label] = lookup[(old_number, insertion)]
    missing = [n for n in selected_numbers if f"{chain}{n}" not in old_to_new]
    if missing:
        raise ValueError("Selected β-strand residues could not be mapped in the PDB: " + ", ".join(f"{chain}{n}" for n in missing[:8]))

    normalized_to_old = {new: f"{chain}{old}{insertion}" for (old, insertion), new in lookup.items()}
    selected_positions = sorted({old_to_new[f"{chain}{n}"] for n in selected_numbers})
    groups: list[list[int]] = []
    for position in selected_positions:
        if not groups or position != groups[-1][-1] + 1:
            groups.append([])
        groups[-1].append(position)
    if any(len(group) < 3 for group in groups):
        raise ValueError("Each selected target β-strand segment must contain at least three consecutive protein residues.")
    if len(groups) > 4:
        raise ValueError("Select no more than four target β-strand segments in one job.")

    output_lines: list[str] = []
    serial = 1
    for line in source_lines:
        if not line.startswith("ATOM  ") or line[21:22] != chain:
            continue
        resname = line[17:20].strip().upper()
        if resname not in inputs.AA3 or line[16:17] not in {" ", "A"}:
            continue
        key = (int(line[22:26]), line[26:27].strip())
        new_number = lookup.get(key)
        if new_number is None:
            continue
        atom_line = line
        if resname == "MSE":
            atom_name = line[12:16].strip()
            if atom_name == "SE":
                atom_line = line[:12] + " SD " + line[16:]
                atom_line = atom_line[:76] + " S " + atom_line[79:]
            atom_line = atom_line[:17] + "MET" + atom_line[20:]
        atom_line = atom_line[:6] + f"{serial:5d}" + atom_line[11:21] + "A" + f"{new_number:4d}" + " " + atom_line[27:]
        # Remove alternate-location labels after keeping blank/A atoms.
        atom_line = atom_line[:16] + " " + atom_line[17:]
        output_lines.append(atom_line)
        serial += 1
    if not output_lines:
        raise ValueError("No protein coordinates remained after target-chain normalization.")
    destination.write_text("\n".join([*output_lines, "TER", "END", ""]))
    return groups, normalized_to_old


def _collapse(group: list[int]) -> str:
    return f"A{group[0]}-{group[-1]}"


def _docker_tensor_command(run_dir: Path, target: Path, tensor_dir: Path, p: dict[str, Any]) -> list[str]:
    command = [
        "docker", "run", "--rm", "--shm-size=2G", "--entrypoint", "/bin/bash",
        "-v", f"{run_dir}:/work", "-v", f"{inputs.TENSOR_SCRIPT.parent}:/scripts:ro", "-w", "/work", TENSOR_IMAGE,
        "-c", 'exec "$@"', "studio-beta-tensors", "python", "/scripts/make_interface_tensor.py",
        "--input_pdb", _container_path(run_dir, target),
        "--out_dir", _container_path(run_dir, tensor_dir),
        "--binderlen", str(p["binder_length"]),
        "--target_adj", p["target_adj"],
        "--binder_ss", ",".join("E" for _ in p["target_strand_groups"]),
        "--binder_ss_len", ",".join(str(p["binder_strand_length"]) for _ in p["target_strand_groups"]),
    ]
    return command


def _docker_rfdiffusion_command(run_dir: Path, target: Path, tensor_dir: Path, p: dict[str, Any]) -> list[str]:
    models = reference_root() / "rfdiffusion_models"
    output = run_dir / "artifacts" / "raw" / inputs.TOOL / "output" / "design"
    binders = p["binder_length"]
    target_length = len(p["target_residue_map"])
    contig = f"{binders}/0 A1-{target_length}"
    hotspots = ",".join(p["target_beta_residue_ids"])
    seed = p["seed"] or secrets.randbelow(2147483646) + 1
    p["resolved_seed"] = seed
    command = [
        "docker", "run", "--rm", *docker_gpu_args(p["gpu_device"]), "--shm-size=8G",
        "--entrypoint", "/bin/bash", "-v", f"{run_dir}:/work",
        "-v", f"{models}:/models:ro", "-w", "/opt/RFdiffusion", RF_IMAGE,
        "-c", 'export PYTHONPATH=/opt/RFdiffusion:/opt/RFdiffusion/env/SE3Transformer; exec "$@"',
        "studio-rfdiffusion-beta", "python3", "/opt/RFdiffusion/scripts/run_inference.py",
        f"inference.output_prefix={_container_path(run_dir, output)}",
        f"inference.input_pdb={_container_path(run_dir, target)}",
        f"inference.num_designs={p['num_designs']}",
        f"inference.model_directory_path=/models",
        f"inference.schedule_directory_path=/opt/RFdiffusion/schedules",
        f"inference.model_runner=ScaffoldedSampler",
        f"inference.design_startnum={seed}",
        f"diffuser.T={p['timesteps']}",
        f"contigmap.contigs=[{contig}]",
        # The installed mn-rfdiffusion image predates the upstream config-key
        # rename: it selects ScaffoldedSampler with scaffoldguided_enable.
        "scaffoldguided.scaffoldguided_enable=True",
        f"scaffoldguided.scaffold_dir={_container_path(run_dir, tensor_dir)}",
        f"scaffoldguided.scaffold_list={_container_path(run_dir, tensor_dir.parent / (tensor_dir.name + '.txt'))}",
        f"scaffoldguided.systematic={str(p['systematic_sampling']).lower()}",
        "scaffoldguided.sampled_insertion=0",
        "scaffoldguided.sampled_N=0",
        "scaffoldguided.sampled_C=0",
        f"ppi.hotspot_res=[{hotspots}]",
        "inference.deterministic=True",
    ]
    return command


def run_rfdiffusion_beta_pairing(
    target_pdb: Path,
    target_chain: str,
    target_beta_strands: str,
    binder_length: int = 65,
    binder_strand_length: int = 6,
    num_designs: int = 10,
    timesteps: int = 50,
    systematic_sampling: bool = False,
    seed: int = 0,
    gpu_device: object = "0",
    cpu_cores: int = 4,
    job_name: str = "",
) -> Path:
    p = inputs.prepare_request(locals())
    model_weights = p["model_weights"]
    p["gpu_device"] = normalize_gpu_device(gpu_device)
    selected_numbers = inputs.parse_residue_selection(target_beta_strands, target_chain)
    job = create_job(
        "design", "design_generation", inputs.TOOL,
        {"target_pdb": str(target_pdb), "target_chain": target_chain, "target_beta_strands": target_beta_strands},
        p,
    )
    raw = job.run_dir / "artifacts" / "raw" / inputs.TOOL
    input_dir = raw / "input"
    tensor_dir = raw / "interface_tensors"
    output_dir = raw / "output"
    input_dir.mkdir(parents=True, exist_ok=True)
    tensor_dir.mkdir(parents=True, exist_ok=True)
    target = input_dir / "target.pdb"
    groups, target_residue_map = _normalize_target(Path(target_pdb), target_chain, selected_numbers, target)
    p["target_residue_map"] = target_residue_map
    p["normalized_target_chain"] = "A"
    p["target_strand_groups"] = [[f"A{n}" for n in group] for group in groups]
    p["target_adj"] = ":".join(_collapse(group) for group in groups)
    p["target_beta_strands_normalized"] = p["target_adj"]
    p["target_beta_residue_ids"] = [f"A{n}" for group in groups for n in group]
    p["register_count"] = inputs.estimate_register_count(binder_length, binder_strand_length, len(groups))
    p["tensor_generator"] = "published interface_tensors/make_interface_tensor.py (Studio fixes: inclusive final register and 1-based PDB offset)"
    p["rfdiffusion_source_revision"] = RF_SOURCE_REVISION
    p["engine_version"] = f"{RF_IMAGE} · β-pairing interface tensors"
    p["rfdiffusion_runtime_config"] = "scaffoldguided.scaffoldguided_enable (installed image config)"
    p["campaign_name"] = job_name.strip() or "RFdiffusion β-pairing"
    output_dir.mkdir(parents=True, exist_ok=True)
    write_json(job.run_dir / "input.json", {
        "inputs": {"target_pdb": str(target_pdb), "target_chain": target_chain, "target_beta_strands": target_beta_strands},
        "params": p,
    })
    update_status(job.run_dir, "preparing", campaign_name=p["campaign_name"], current_phase="Preparing β-pairing interface tensors")
    tensor_command = _docker_tensor_command(job.run_dir, target, tensor_dir, p)
    steps = [
        {"name": "RFdiffusion β-pairing · build SS/adjacency tensors", "command": tensor_command},
        {"name": "RFdiffusion β-pairing · generate backbones", "command": _docker_rfdiffusion_command(job.run_dir, target, tensor_dir, p)},
    ]
    write_json(raw / "run_config.json", {"contig": f"{binder_length}/0 A1-{len(target_residue_map)}", "parameters": p})
    rc = design._run_shell_steps(job.run_dir, steps)
    output_paths = sorted(output_dir.glob("design_*.pdb")) if rc == 0 else []
    candidates = []
    for index, structure in enumerate(output_paths, start=1):
        binders, targets, inference = design._infer_generated_chain_roles(target, structure, ["A"])
        candidates.append({
            "candidate_id": f"rfdiffusion_beta_pairing_{index:05d}",
            "source_tool": inputs.TOOL,
            "stage": "generation.backbone",
            "target_pdb": str(target.relative_to(job.run_dir)),
            "complex_pdb": str(structure.relative_to(job.run_dir)),
            "binder_chains": binders,
            "target_chains": targets,
            "binder_sequence": None,
            "binder_length": binder_length,
            "contig": f"{binder_length}/0 A1-{len(target_residue_map)}",
            "hotspots": p["target_beta_residue_ids"],
            "metrics": {"result_kind": "generated", "pass_filters": None},
            "raw_metadata": {
                "generator_only": True,
                "generator_engine": inputs.TOOL,
                "result_kind": "generated",
                "conditioning_method": "RFdiffusion SS/adjacency interface tensor conditioning",
                "conditioned_binder_secondary_structure": "E",
                "target_beta_strands": p["target_beta_strands_normalized"],
                "target_chain_inference": inference,
                "register_count": p["register_count"],
                "systematic_register_sampling": systematic_sampling,
                "sequence_design_required": True,
                "model_weights": model_weights,
                "seed": p["resolved_seed"],
            },
        })
    normalized = write_candidates(job.run_dir, inputs.TOOL, candidates) if candidates else []
    design._finish_design_job(
        job.run_dir, rc == 0 and bool(normalized), rc, normalized,
        [(f"raw/{inputs.TOOL}/**/*.{extension}", extension) for extension in ("pdb", "pt", "txt", "json", "trb", "log")],
    )
    return job.run_dir
