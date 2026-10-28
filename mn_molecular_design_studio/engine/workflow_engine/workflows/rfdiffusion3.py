"""RFdiffusion3 tasks using the existing, pinned Foundry container."""
from __future__ import annotations

import csv
import math
import secrets
import shutil
import subprocess
import traceback
from pathlib import Path
from typing import Literal

from mn_molecular_design_studio.engine.workflow_engine.core.execution_timing import execution_phase
from mn_molecular_design_studio import rfdiffusion3_input as inputs
from mn_molecular_design_studio.design_validation import protein_chains, read_structure, sequence
from mn_molecular_design_studio.engine.workflow_engine.core.candidates import STAGE_GENERATION_BACKBONE, write_candidates
from mn_molecular_design_studio.engine.workflow_engine.core.gpu import docker_gpu_args, normalize_gpu_device
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import create_job, finish_job, read_json, update_status, write_json
from mn_molecular_design_studio.engine.workflow_engine.runtime import reference_root


def _normalize(run: Path, prepared: dict, requested: int) -> list[dict]:
    root = run / "artifacts/raw/rfdiffusion3"
    rows = []
    for identifier, context in prepared["contexts"].items():
        paths = sorted((root / "output").glob(f"studio_{identifier}_*.cif*"))
        for index, path in enumerate(paths[:requested], 1):
            # Keep a standard, uncompressed CIF alongside the native output.
            # Existing Studio viewers and downstream MPNN routes read this path.
            from Bio.PDB import MMCIFIO
            structure = read_structure(path)
            destination = root / "structures" / f"{identifier}_{index:05d}.cif"
            destination.parent.mkdir(parents=True, exist_ok=True)
            writer = MMCIFIO(); writer.set_structure(structure); writer.save(str(destination))
            chains = protein_chains(structure)
            target_ids = []
            if context["source_path"] and context["task"] == "protein_binder":
                original = protein_chains(read_structure(Path(context["source_path"])))
                for source_id in context["target_chains"]:
                    expected = sequence(original[source_id])
                    matches = [c for c, residues in chains.items() if c not in target_ids and sequence(residues) == expected]
                    if source_id in matches:
                        target_ids.append(source_id)
                    elif len(matches) == 1:
                        target_ids.append(matches[0])
                    else:
                        raise ValueError("Could not map the retained protein target into the RFdiffusion3 output.")
            binder_ids = [c for c in chains if c not in target_ids]
            if not binder_ids:
                raise ValueError("RFdiffusion3 returned no generated protein chain.")
            native_json = path.with_name(path.name.removesuffix(".gz").removesuffix(".cif") + ".json")
            metadata = read_json(native_json)
            metrics = {key: value for key, value in (metadata.get("metrics") or {}).items() if isinstance(value, (int, float, bool))}
            metrics.update({key: value for key, value in metadata.items() if isinstance(value, (int, float, bool))})
            fixed_residues = []
            if context["source_path"]:
                source_info = inputs.inventory(Path(context["source_path"]))
                unfixed = prepared["specifications"][identifier].get("select_unfixed_sequence")
                unfixed_labels = set()
                if isinstance(unfixed, str):
                    unfixed_labels = {r["label"] for r in inputs._selected_rows(unfixed, source_info)}
                elif isinstance(unfixed, dict):
                    for selector in unfixed:
                        unfixed_labels.update(r["label"] for r in inputs._selected_rows(selector, source_info))
                for source_label, output_label in (metadata.get("diffused_index_map") or {}).items():
                    if unfixed is True or source_label in unfixed_labels:
                        continue
                    source_row = source_info["residues"].get(source_label, {})
                    if source_row.get("kind") == "protein" and any(output_label == f"{c}{r.id[1]}" for c in binder_ids for r in chains[c]):
                        fixed_residues.append(output_label)
            metrics.update(generated_residue_count=sum(len(chains[c]) for c in binder_ids),
                           sequence_design_required=True, result_kind="generation_only")
            rows.append({
                "candidate_id": f"rfdiffusion3_{identifier}_{index:05d}",
                "source_tool": "rfdiffusion3_foundry", "stage": STAGE_GENERATION_BACKBONE,
                "target_pdb": context["source_path"], "complex_pdb": str(destination), "binder_sequence": None,
                "binder_chains": binder_ids, "target_chains": target_ids, "metrics": metrics,
                "raw_metadata": {"task": context["task"], "specification_name": context["name"],
                                 "result_kind": "generation_only", "native_structure": str(path),
                                 "metadata_path": str(native_json), "engine_version": inputs.CODE_VERSION,
                                 "input_specification": prepared["specifications"][identifier],
                                 "target_entity_chains": context["target_chains"],
                                 "input_target_chains": context["target_chains"],
                                 "fixed_residues": sorted(set(fixed_residues)),
                                 "diffused_index_map": metadata.get("diffused_index_map") or {},
                                 "sequence_design_required": True},
            })
    return rows


def run_rfdiffusion3_design(
    task: Literal["protein_binder", "enzyme", "nucleic_acid_binder", "ligand_binder", "json"] = "protein_binder",
    target_structure: Path | None = None,
    target_chains: list[str] | None = None,
    length_min: int = 100, length_max: int = 150,
    hotspots: str = "", ligands: list[str] | None = None,
    scaffolding_mode: Literal["atom", "residue"] = "atom",
    fixed_atoms: dict | None = None, unindex: str = "all", motif_contig: str = "",
    atom_conditions: list[dict] | None = None, ori_token: list[float] | None = None,
    json_config: dict | None = None, input_files: list[Path] | None = None,
    num_designs: int = 8, diffusion_batch_size: int = 1,
    timesteps: int = 200, step_scale: float = 1.5, gamma_0: float = .6,
    is_non_loopy: bool = True, classifier_free_guidance: bool = False, cfg_scale: float = 1.5,
    low_memory_mode: bool = False, seed: int = 0,
    job_name: str = "rfdiffusion3-design", cpu_cores: int = 4, gpu_device: str = "0",
) -> Path:
    parameters = dict(locals())
    parameters["target_structure"] = str(target_structure) if target_structure else None
    parameters["input_files"] = [str(path) for path in input_files or []]
    prepared = inputs.prepare_request(parameters)
    seed = seed or secrets.randbelow(2147483646) + 1
    parameters.update(seed=seed, gpu_device=normalize_gpu_device(gpu_device), engine_version=inputs.CODE_VERSION, image=inputs.IMAGE)
    job = create_job("design", "rfdiffusion3_design", "rfdiffusion3", {"target_structure": str(target_structure) if target_structure else None}, parameters)
    update_status(job.run_dir, "preparing", campaign_name=job_name, current_phase="Preparing RFdiffusion3 inputs")
    root = job.run_dir / "artifacts/raw/rfdiffusion3"
    output = root / "output"; output.mkdir(parents=True, exist_ok=True)
    try:
        paths = {}
        for name, source in prepared["files"].items():
            destination = root / "inputs" / name
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)
            paths[name] = "/work/" + destination.relative_to(job.run_dir).as_posix()
        for identifier, spec in prepared["specifications"].items():
            if spec.get("input"):
                spec["input"] = paths[spec["input"]]
            context = prepared["contexts"][identifier]
            if context["source_path"]:
                context["source_path"] = str(root / "inputs" / Path(context["source_path"]).name)
        configuration = root / "inputs.json"; write_json(configuration, prepared["specifications"])
        write_json(root / "input_mapping.json", prepared["contexts"])
        weights = Path(reference_root()) / "foundry"
        if not (weights / "rfd3_latest.ckpt").is_file():
            raise FileNotFoundError(f"RFdiffusion3 checkpoint is missing from {weights}.")
        batch_size = min(diffusion_batch_size, num_designs)
        batches = math.ceil(num_designs / batch_size)
        command = [
            "docker", "run", "--rm", "--name", f"mn-studio-rfd3-{job.run_id[-8:]}",
            "--cpus", str(cpu_cores), "--shm-size=16G", *docker_gpu_args(gpu_device),
            "--mount", f"type=bind,src={job.run_dir.resolve()},dst=/work",
            "--mount", f"type=bind,src={weights.resolve()},dst=/weights,readonly",
            "-e", f"OMP_NUM_THREADS={cpu_cores}", "-e", "PYTHONUNBUFFERED=1", "-w", "/work", "--entrypoint", "rfd3", inputs.IMAGE,
            "design", f"inputs=/work/{configuration.relative_to(job.run_dir).as_posix()}",
            f"out_dir=/work/{output.relative_to(job.run_dir).as_posix()}", "ckpt_path=/weights/rfd3_latest.ckpt",
            f"diffusion_batch_size={batch_size}", f"n_batches={batches}",
            f"inference_sampler.num_timesteps={timesteps}", f"inference_sampler.step_scale={step_scale}",
            f"inference_sampler.gamma_0={gamma_0}", f"inference_sampler.use_classifier_free_guidance={classifier_free_guidance}",
            f"inference_sampler.cfg_scale={cfg_scale}", f"low_memory_mode={low_memory_mode}", f"seed={seed}",
            "skip_existing=False", "prevalidate_inputs=True", "global_prefix=studio_",
        ]
        write_json(job.run_dir / "command.json", {"mode": "docker", "image": inputs.IMAGE, "command": command})
        update_status(job.run_dir, "running", current_phase=f"Generating {num_designs} design(s) per task input")
        with (job.run_dir / "stdout.log").open("a") as stdout, (job.run_dir / "stderr.log").open("a") as stderr:
            with execution_phase(job.run_dir, "docker_command") as timing:
                process = subprocess.run(command, stdout=stdout, stderr=stderr, cwd=job.run_dir, check=False)
                timing["return_code"] = int(process.returncode)
        update_status(job.run_dir, "running", current_phase="Collecting RFdiffusion3 designs")
        with execution_phase(job.run_dir, "candidate_normalization"):
            candidates = _normalize(job.run_dir, prepared, num_designs)
            normalized = write_candidates(job.run_dir, "rfdiffusion3_foundry", candidates) if candidates else []
        table = root / "designs.csv"
        with table.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["candidate_id", "task", "input", "designed_chains", "generated_residue_count", "structure"])
            writer.writeheader()
            writer.writerows({"candidate_id": c["candidate_id"], "task": task,
                             "input": c["raw_metadata"]["specification_name"], "designed_chains": ",".join(c["binder_chains"]),
                             "generated_residue_count": c["metrics"]["generated_residue_count"], "structure": c["complex_pdb"]} for c in candidates)
        artifacts = [{"path": str(path), "name": path.name, "type": "rfdiffusion3_output"}
                     for path in root.rglob("*") if path.is_file()]
        success = process.returncode == 0 and len(normalized) == num_designs * len(prepared["specifications"])
        finish_job(job.run_dir, success, {
            "outputs": {"artifacts": artifacts, "candidates": normalized, "design_table": str(table)},
            "metrics": {"return_code": process.returncode, "task": task, "candidate_count": len(normalized),
                        "requested_designs_per_input": num_designs, "input_count": len(prepared["specifications"]),
                        "native_generated_count": len(list(output.glob("*.cif*"))), "engine_version": inputs.CODE_VERSION},
            "downstream_artifacts": {"candidates_jsonl": "artifacts/normalized_candidates/candidates.jsonl",
                                     "campaign_result": "artifacts/normalized_candidates/campaign_result.json"},
        })
    except Exception as exc:
        with (job.run_dir / "stderr.log").open("a") as handle:
            handle.write(traceback.format_exc())
        finish_job(job.run_dir, False, {"error": str(exc), "metrics": {"task": task, "error": str(exc)}})
    return job.run_dir
