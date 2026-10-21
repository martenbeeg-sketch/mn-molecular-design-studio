"""Queued RFDiffusion4-Proteina co-design generation workflow."""
from __future__ import annotations

import csv
import gzip
import shutil
import secrets
import traceback
from pathlib import Path
from typing import Literal

from Bio.PDB import MMCIFIO, Select

from mn_molecular_design_studio import rfd4_proteina_input as inputs
from mn_molecular_design_studio.design_validation import protein_chains, read_structure, sequence
from mn_molecular_design_studio.engine.workflow_engine.core.candidates import (
    STAGE_GENERATION_BACKBONE_SEQUENCE,
    write_candidates,
)
from mn_molecular_design_studio.engine.workflow_engine.core.docker_runner import run_docker_job
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import (
    create_job,
    finish_job,
    update_status,
    write_json,
)
from mn_molecular_design_studio.engine.workflow_engine.core.manifests import load_manifest
from mn_molecular_design_studio.engine.workflow_engine.runtime import reference_root


Task = Literal[
    "protein_binder", "nanobody", "antibody", "tcr", "motif_scaffolding", "enzyme",
    "small_molecule_binder", "dna_binder", "rna_design", "existing_structure",
    "custom_conditions", "unconditional",
]
InputMode = Literal["guided_binder", "annotated_cif", "condition_spec", "unconditional"]


def _target_sequences(paths: list[Path]) -> dict[str, str]:
    found: dict[str, str] = {}
    for path in paths:
        if path.suffix.lower() not in {".cif", ".mmcif", ".pdb", ".ent"}:
            continue
        try:
            for chain, residues in protein_chains(read_structure(path)).items():
                value = sequence(residues)
                if value:
                    found[f"{path.stem}:{chain}"] = value
        except (OSError, ValueError, KeyError, StopIteration):
            continue
    return found


def _materialize_cif(path: Path, output_root: Path) -> Path:
    name = path.name
    if name.endswith(".cif.zst"):
        import zstandard

        destination = output_root / "normalized_structures" / name.removesuffix(".zst")
        destination.parent.mkdir(parents=True, exist_ok=True)
        with path.open("rb") as source, destination.open("wb") as target:
            zstandard.ZstdDecompressor().copy_stream(source, target)
        return destination
    if name.endswith(".cif.gz"):
        destination = output_root / "normalized_structures" / name.removesuffix(".gz")
        destination.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(path, "rb") as source, destination.open("wb") as target:
            shutil.copyfileobj(source, target)
        return destination
    return path


def _normalize_outputs(run_dir: Path, prepared: dict, task: str, requested: int) -> list[dict]:
    output_root = prepared["root"] / "output"
    generated = sorted(
        path for path in output_root.glob("inference_samples/**/*_generated.cif*")
        if path.name.endswith(("_generated.cif", "_generated.cif.zst", "_generated.cif.gz"))
    )
    reference_sequences = _target_sequences(prepared["staged_files"])
    candidates = []
    for index, source_path in enumerate(generated[:10_000], 1):
        path = _materialize_cif(source_path, output_root)
        try:
            structure = read_structure(path)
            chains = protein_chains(structure)
        except (OSError, ValueError, KeyError, StopIteration):
            continue
        target_chains = []
        if reference_sequences:
            reference_values = set(reference_sequences.values())
            target_chains = [chain for chain, residues in chains.items() if sequence(residues) in reference_values]
        binder_chains = [chain for chain in chains if chain not in target_chains]
        if not binder_chains:
            # For partial diffusion or a complex task with ambiguous chain roles,
            # retain the generated structure and leave role assignment unset.
            binder_chains = []
        binder_sequence = "".join(sequence(chains[chain]) for chain in binder_chains) or None
        total_residues = sum(len(residues) for residues in chains.values())
        candidate_name = path.name.removesuffix(".zst").removesuffix("_generated.cif")
        candidate_id = f"rfd4_proteina_{candidate_name or index}"
        binder_path = path
        if binder_chains:
            class BinderChainSelect(Select):
                def accept_chain(self, chain):
                    return chain.id in binder_chains

            binder_path = output_root / "normalized_structures" / f"{candidate_name or index}_binder.cif"
            binder_path.parent.mkdir(parents=True, exist_ok=True)
            writer = MMCIFIO()
            writer.set_structure(structure)
            writer.save(str(binder_path), select=BinderChainSelect())
        candidates.append({
            "candidate_id": candidate_id,
            "source_tool": "rfd4_proteina",
            "stage": STAGE_GENERATION_BACKBONE_SEQUENCE,
            "target_pdb": prepared.get("target_structure"),
            "complex_pdb": str(path),
            "binder_pdb": str(binder_path),
            "binder_sequence": binder_sequence,
            "target_chains": target_chains,
            "binder_chains": binder_chains,
            "metrics": {"generated_residue_count": total_residues, "result_kind": "sequence_structure_codesign"},
            "raw_metadata": {
                "task": task,
                "result_kind": "sequence_structure_codesign",
                "native_structure": str(source_path),
                "display_structure": str(path),
                "input_condition_spec": str(prepared["spec_path"]),
                "target_chain_sequences": prepared.get("target_sequences") or reference_sequences,
                "engine_version": inputs.UPSTREAM_VERSION,
                "upstream_commit": inputs.UPSTREAM_COMMIT,
                "requested_designs": requested,
            },
        })
    return candidates


def run_rfd4_proteina_design(
    task: Task = "protein_binder",
    input_mode: InputMode = "guided_binder",
    target_structure: Path | None = None,
    target_chains: list[str] | None = None,
    hotspots: str = "",
    length_min: int = 60,
    length_max: int = 120,
    condition_spec: str = "",
    input_files: list[Path] | None = None,
    num_designs: int = 8,
    batch_size: int = 1,
    nsteps: int = 60,
    seed: int = 0,
    compile_model: bool = False,
    cpu_cores: int = 8,
    gpu_device: str = "all",
    job_name: str = "rfd4-proteina-design",
) -> Path:
    """Generate protein sequence–structure designs with the pinned RFD4 image."""
    parameters = dict(locals())
    prepared = None
    job = None
    try:
        # Create the job first so all user inputs and native outputs are retained
        # under its result directory and remain available from Results.
        source_inputs = {
            "target_structure": str(target_structure) if target_structure else None,
            "input_files": [str(path) for path in input_files or []],
        }
        job = create_job("design", "rfd4_proteina_design", "rfd4_proteina", source_inputs, parameters)
        update_status(job.run_dir, "preparing", campaign_name=job_name, current_phase="Preparing RFDiffusion4-Proteina inputs")
        prepared = inputs.prepare_request(parameters, run_dir=job.run_dir)
        if not seed:
            seed = secrets.randbelow(2_147_483_646) + 1
        parameters.update(
            seed=seed,
            gpu_device=gpu_device,
            engine_version=inputs.UPSTREAM_VERSION,
            upstream_commit=inputs.UPSTREAM_COMMIT,
            image=inputs.IMAGE,
        )
        write_json(job.run_dir / "parameters.json", parameters)
        output_root = prepared["root"] / "output"
        output_root.mkdir(parents=True, exist_ok=True)
        manifest = load_manifest("rfd4_proteina")
        arguments = [
            "inputs=/work/" + prepared["spec_path"].relative_to(job.run_dir).as_posix(),
            f"num_replicates={num_designs}",
            f"batch_size={batch_size}",
            f"generation.args.nsteps={nsteps}",
            f"seed={seed}",
            "metrics_tags=[]",
            f"out_dir=/work/{output_root.relative_to(job.run_dir).as_posix()}",
        ]
        if not compile_model:
            arguments.append("compile_lengths=null")
        cache_dir = Path(reference_root()) / "rfd4-proteina" / "huggingface"
        cache_dir.mkdir(parents=True, exist_ok=True)
        update_status(job.run_dir, "preparing", current_phase=f"Staged {len(prepared['staged_files'])} structure input(s)")
        return_code = run_docker_job(
            job.run_dir,
            manifest,
            arguments,
            gpu_device=gpu_device,
            extra_mounts=[(cache_dir, "/cache/huggingface", False)],
            environment={
                "HF_HOME": "/cache/huggingface",
                "HF_HUB_CACHE": "/cache/huggingface/hub",
                "OMP_NUM_THREADS": str(cpu_cores),
                "PYTHONUNBUFFERED": "1",
            },
        )
        update_status(job.run_dir, "running", current_phase="Collecting generated RFD4 structures")
        candidates = _normalize_outputs(job.run_dir, prepared, task, num_designs)
        normalized = write_candidates(job.run_dir, "rfd4_proteina", candidates) if candidates else []
        table = prepared["root"] / "designs.csv"
        with table.open("w", newline="") as handle:
            columns = ["candidate_id", "task", "binder_sequence", "binder_chains", "target_chains", "generated_residue_count", "structure"]
            writer = csv.DictWriter(handle, fieldnames=columns)
            writer.writeheader()
            writer.writerows({
                "candidate_id": item["candidate_id"], "task": task,
                "binder_sequence": item.get("binder_sequence") or "",
                "binder_chains": ",".join(item.get("binder_chains") or []),
                "target_chains": ",".join(item.get("target_chains") or []),
                "generated_residue_count": item["metrics"].get("generated_residue_count"),
                "structure": item.get("complex_pdb"),
            } for item in candidates)
        artifacts = [
            {"path": str(path), "name": path.name, "type": "rfd4_proteina_output"}
            for path in prepared["root"].rglob("*") if path.is_file()
        ]
        success = return_code == 0 and bool(normalized)
        finish_job(job.run_dir, success, {
            "outputs": {"artifacts": artifacts, "candidates": normalized, "design_table": str(table)},
            "metrics": {
                "return_code": return_code, "task": task, "input_mode": input_mode,
                "candidate_count": len(normalized), "requested_designs": num_designs,
                "input_design_count": prepared["design_count"],
                "engine_version": inputs.UPSTREAM_VERSION, "upstream_commit": inputs.UPSTREAM_COMMIT,
                "gpu_device": gpu_device,
            },
            "error": None if success else ("RFDiffusion4-Proteina returned no generated CIF structures." if return_code == 0 else f"RFDiffusion4-Proteina exited with status {return_code}."),
            "downstream_artifacts": {"candidates_jsonl": "artifacts/normalized_candidates/candidates.jsonl"},
        })
    except Exception as exc:
        if job is not None:
            with (job.run_dir / "stderr.log").open("a") as handle:
                handle.write(traceback.format_exc())
            finish_job(job.run_dir, False, {"error": str(exc), "metrics": {"task": task, "error": str(exc)}})
    return job.run_dir if job is not None else Path()
