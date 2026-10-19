"""Standalone RFpeptides macrocyclic peptide backbone design."""
from __future__ import annotations

import csv
import shutil
import subprocess
import traceback
from pathlib import Path

from Bio.PDB import PDBIO, Select

from mn_molecular_design_studio import rfpeptides_input as inputs
from mn_molecular_design_studio.design_validation import protein_chains, read_structure
from mn_molecular_design_studio.engine.workflow_engine.core.candidates import STAGE_GENERATION_BACKBONE, write_candidates
from mn_molecular_design_studio.engine.workflow_engine.core.docker_runner import build_docker_command
from mn_molecular_design_studio.engine.workflow_engine.core.gpu import normalize_gpu_device
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import create_job, finish_job, update_status, write_json
from mn_molecular_design_studio.engine.workflow_engine.core.manifests import load_manifest
from mn_molecular_design_studio.engine.workflow_engine.runtime import reference_root


def _write_binder_structure(source: Path, destination: Path, binder_chain: str = "A") -> None:
    structure = read_structure(source)
    model = next(structure.get_models(), None)
    if model is None or binder_chain not in model:
        return

    class _BinderOnly(Select):
        def accept_chain(self, chain):
            return chain.id == binder_chain

        def accept_residue(self, residue):
            return residue.id[0] == " " and "CA" in residue

        def accept_atom(self, atom):
            return atom.element != "H" and (not atom.is_disordered() or atom.get_altloc() in {"A", " "})

    destination.parent.mkdir(parents=True, exist_ok=True)
    writer = PDBIO()
    writer.set_structure(structure)
    writer.save(str(destination), _BinderOnly())


def _collect_candidates(run_dir: Path, target_path: Path, mapping: dict, hotspots: list[str]) -> list[dict]:
    output_dir = run_dir / "artifacts" / "raw" / "rfpeptides" / "output"
    binder_dir = run_dir / "artifacts" / "raw" / "rfpeptides" / "peptides"
    candidates = []
    for index, structure_path in enumerate(sorted(output_dir.glob("rfpeptides_*.pdb")), 1):
        try:
            structure = read_structure(structure_path)
            chains = protein_chains(structure)
            if "A" not in chains:
                raise ValueError("RFdiffusion output does not contain the expected cyclic peptide chain A.")
            target_chains = [chain for chain in mapping["target_chains"] if chain in chains]
            if not target_chains:
                raise ValueError("RFdiffusion output does not contain the selected target chains.")
            binder_path = binder_dir / f"{structure_path.stem}_peptide.pdb"
            _write_binder_structure(structure_path, binder_path)
            binder_length = len(chains["A"])
            candidates.append({
                "candidate_id": structure_path.stem,
                "source_tool": "rfpeptides",
                "stage": STAGE_GENERATION_BACKBONE,
                "target_pdb": str(target_path.relative_to(run_dir)),
                "complex_pdb": str(structure_path.relative_to(run_dir)),
                "binder_pdb": str(binder_path.relative_to(run_dir)),
                "binder_chains": ["A"],
                "target_chains": target_chains,
                "hotspots": hotspots,
                "binder_length": binder_length,
                "metrics": {"design_index": index, "binder_length": binder_length, "cyclic": True},
                "raw_metadata": {
                    "upstream_commit": inputs.UPSTREAM_COMMIT,
                    "result_kind": "cyclic_peptide_backbone",
                    "sequence_design_status": "not_run",
                    "target_chain_map": mapping["chain_map"],
                },
            })
        except (OSError, ValueError, KeyError) as exc:
            with (output_dir / "collection_warnings.txt").open("a") as handle:
                handle.write(f"{structure_path.name}: {type(exc).__name__}: {exc}\n")
    return candidates


def run_rfpeptides_design(
    target_pdb: Path | None = None,
    target_chains: list[str] | None = None,
    binder_length: str = "12-18",
    hotspots: str = "",
    number_of_designs: int = 1,
    diffusion_steps: int = 50,
    cpu_cores: int = 4,
    gpu_device: str = "0",
    job_name: str = "rfpeptides-design",
) -> Path:
    parameters = dict(locals())
    normalized = inputs.validate_request(parameters, require_files=True)
    models_dir = reference_root() / "rfdiffusion_models"
    checkpoint = reference_root() / inputs.WEIGHTS_RELATIVE_PATH
    if not checkpoint.is_file():
        raise FileNotFoundError(
            f"RFpeptides requires {checkpoint}. Download Complex_base_ckpt.pt with the official "
            "RFdiffusion model download script into the Studio reference directory."
        )
    gpu_device = normalize_gpu_device(gpu_device)
    parameters.update({
        "target_pdb": str(normalized["target_pdb"]),
        "target_chains": normalized["target_chains"],
        "hotspots": normalized["hotspots"],
        "binder_length": normalized["binder_length"],
        "number_of_designs": normalized["number_of_designs"],
        "diffusion_steps": normalized["diffusion_steps"],
        "cpu_cores": normalized["cpu_cores"],
        "gpu_device": gpu_device,
        "engine_version": inputs.CODE_VERSION,
        "image": inputs.IMAGE,
    })
    job = create_job(
        "design", "rfpeptides_design", "rfpeptides",
        {"target_pdb": str(normalized["target_pdb"]), "target_chains": normalized["target_chains"]},
        parameters,
    )
    update_status(job.run_dir, "preparing", campaign_name=job_name, current_phase="Preparing cyclic peptide target")
    root = job.run_dir / "artifacts" / "raw" / "rfpeptides"
    inputs_dir = root / "inputs"
    output_dir = root / "output"
    inputs_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    target_path = inputs_dir / "target.pdb"
    try:
        mapping = inputs.prepare_target(normalized["target_pdb"], normalized["target_chains"], normalized["hotspots"], target_path)
        shutil.copy2(normalized["target_pdb"], inputs_dir / f"source_target{normalized['target_pdb'].suffix.lower()}")
        contig = inputs.build_contig(normalized["binder_length"], mapping["target_chains"], mapping["chain_lengths"])
        args = [
            "--config-name", "base",
            f"inference.output_prefix=/work/artifacts/raw/rfpeptides/output/rfpeptides",
            f"inference.input_pdb=/work/artifacts/raw/rfpeptides/inputs/target.pdb",
            f"inference.model_directory_path=/ref/rfdiffusion_models",
            f"inference.ckpt_override_path=/ref/rfdiffusion_models/Complex_base_ckpt.pt",
            f"inference.num_designs={normalized['number_of_designs']}",
            f"inference.cyclic=True",
            "inference.cyc_chains=a",
            "inference.write_trajectory=False",
            f"diffuser.T={normalized['diffusion_steps']}",
            f"contigmap.contigs=[{contig}]",
        ]
        if mapping["hotspots"]:
            hotspot_list = ",".join(mapping["hotspots"])
            args.append(f"ppi.hotspot_res=[{hotspot_list}]")
        manifest = load_manifest("rfpeptides")
        command = build_docker_command(manifest, job.run_dir, args, gpu_device=gpu_device)
        image_position = command.index(manifest["image"])
        command[image_position:image_position] = [
            "--name", f"mn-rfpeptides-{job.run_id[-8:]}",
            "--shm-size=16G",
            "-e", f"OMP_NUM_THREADS={normalized['cpu_cores']}",
            "-e", "PYTHONUNBUFFERED=1",
            "-e", "TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD=1",
        ]
        write_json(job.run_dir / "command.json", {"mode": "docker", "command": command, "image": inputs.IMAGE})
        update_status(job.run_dir, "running", current_phase="RFpeptides cyclic backbone generation")
        with (job.run_dir / "stdout.log").open("w") as stdout, (job.run_dir / "stderr.log").open("w") as stderr:
            process = subprocess.run(command, stdout=stdout, stderr=stderr, cwd=job.run_dir, check=False)
        update_status(job.run_dir, "running", current_phase="Collecting cyclic peptide backbones")
        candidates = _collect_candidates(job.run_dir, target_path, mapping, mapping["hotspots"])
        normalized_candidates = write_candidates(job.run_dir, "rfpeptides", candidates) if candidates else []
        table_path = root / "designs.csv"
        with table_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["rank", "candidate_id", "peptide_chain", "length", "cyclic", "structure"])
            writer.writeheader()
            for rank, candidate in enumerate(normalized_candidates, 1):
                writer.writerow({
                    "rank": rank,
                    "candidate_id": candidate["candidate_id"],
                    "peptide_chain": "A",
                    "length": candidate.get("binder_length"),
                    "cyclic": True,
                    "structure": candidate.get("complex_pdb"),
                })
        artifacts = [{"path": str(path), "name": path.name, "type": "rfpeptides_output"} for path in root.rglob("*") if path.is_file()]
        success = process.returncode == 0 and bool(normalized_candidates)
        finish_job(job.run_dir, success, {
            "outputs": {"artifacts": artifacts, "candidates": normalized_candidates, "designs_table": str(table_path)},
            "metrics": {
                "return_code": process.returncode,
                "candidate_count": len(normalized_candidates),
                "requested_designs": normalized["number_of_designs"],
                "binder_length": normalized["binder_length"],
                "diffusion_steps": normalized["diffusion_steps"],
                "cyclic": True,
                "sequence_design": "not included; use ProteinMPNN or a sequence-design workflow",
                "engine_version": inputs.CODE_VERSION,
                "target_chain_map": mapping["chain_map"],
                "selected_hotspots": mapping["hotspots"],
            },
            "downstream_artifacts": {
                "candidates_jsonl": "artifacts/normalized_candidates/candidates.jsonl",
                "campaign_result": "artifacts/normalized_candidates/campaign_result.json",
                "designs_table": str(table_path.relative_to(job.run_dir)),
            },
        })
    except Exception as exc:
        with (job.run_dir / "stderr.log").open("a") as handle:
            handle.write(traceback.format_exc())
        finish_job(job.run_dir, False, {"error": str(exc), "metrics": {"error": str(exc), "engine_version": inputs.CODE_VERSION}})
    return job.run_dir
