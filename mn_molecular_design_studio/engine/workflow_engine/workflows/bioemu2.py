"""BioEmu-2 structure sampling jobs for protein targets and designs."""
from __future__ import annotations

import hashlib
import json
import secrets
import shutil
from pathlib import Path
from typing import Literal

from mn_molecular_design_studio.engine.workflow_engine.core.jobs import (
    create_job,
    finish_job,
    update_status,
    write_json,
)
from mn_molecular_design_studio.engine.workflow_engine.core.docker_runner import run_docker_job
from mn_molecular_design_studio.engine.workflow_engine.core.manifests import load_manifest
from mn_molecular_design_studio.engine.workflow_engine.runtime import reference_root


Mode = Literal["propagate", "sample-equilibrium", "predict-structure"]
Sampler = Literal["edm", "edm_gd", "edm_fk", "edm_fkc"]
PhysicalityFilter = Literal["high_quality", "high_quality_per_residue", "no_violations", "no_filter"]
Precision = Literal["bf16-mixed", "32-true"]
MatmulPrecision = Literal["highest", "high", "medium"]


def run_bioemu2(
    mode: Mode = "propagate",
    structure_paths: list[Path] | None = None,
    protein_sequences: list[str] | None = None,
    num_samples: int = 8,
    lag_time_ns: float = 1.0,
    sampler: Sampler = "edm",
    physicality_filter: PhysicalityFilter = "high_quality",
    precision: Precision = "bf16-mixed",
    matmul_precision: MatmulPrecision = "high",
    batch_size: int | None = None,
    num_steps: int | None = None,
    seed: int = 0,
    compile_denoiser: bool = False,
    job_name: str = "bioemu2-sampling",
    cpu_cores: int = 8,
    gpu_device: str = "all",
) -> Path:
    """Sample one or more starting complexes, or a sequence-defined system."""
    structures = [Path(path).expanduser().resolve() for path in structure_paths or []]
    sequences = [str(sequence).strip() for sequence in protein_sequences or [] if str(sequence).strip()]
    if num_samples < 1:
        raise ValueError("Number of samples must be at least one.")
    if cpu_cores < 1:
        raise ValueError("CPU cores must be at least one.")
    if mode == "propagate":
        if not structures:
            raise ValueError("Propagate mode requires at least one uploaded starting structure.")
        if sequences:
            raise ValueError("Provide starting structures or sequences, not both.")
        if not 0 < lag_time_ns <= 1000:
            raise ValueError("BioEmu-2 recommends a propagation lag from above 0 to 1000 ns.")
        for path in structures:
            if not path.is_file() or path.suffix.lower() not in {".pdb", ".ent", ".cif", ".mmcif"}:
                raise ValueError(f"Starting structure must be a PDB or mmCIF file: {path}")
    else:
        if structures:
            raise ValueError("Sequence sampling mode accepts sequences, not starting structures.")
        if not sequences:
            raise ValueError("Provide at least one protein sequence for this mode.")
        if len(sequences) > 26:
            raise ValueError("A sequence-defined system can contain at most 26 protein chains.")
    if batch_size is not None and batch_size < 1:
        raise ValueError("Batch size must be at least one.")
    if num_steps is not None and num_steps < 1:
        raise ValueError("Denoising steps must be at least one.")
    if seed == 0:
        seed = secrets.randbelow(2_147_483_646) + 1

    parameters = {
        "mode": mode,
        "structure_paths": [str(path) for path in structures],
        "protein_sequences": sequences,
        "num_samples": num_samples,
        "lag_time_ns": lag_time_ns,
        "sampler": sampler,
        "physicality_filter": physicality_filter,
        "precision": precision,
        "matmul_precision": matmul_precision,
        "batch_size": batch_size,
        "num_steps": num_steps,
        "seed": seed,
        "compile_denoiser": compile_denoiser,
        "job_name": job_name,
        "cpu_cores": cpu_cores,
        "gpu_device": gpu_device,
    }
    inputs = {
        "structure_paths": [str(path) for path in structures],
        "protein_chain_count": len(sequences),
    }
    job = create_job("simulation", "bioemu2_sampling", "bioemu2", inputs, parameters)
    input_root = job.run_dir / "artifacts" / "input"
    output_root = job.run_dir / "artifacts" / "bioemu2"
    input_root.mkdir(parents=True, exist_ok=True)
    output_root.mkdir(parents=True, exist_ok=True)
    update_status(job.run_dir, "preparing", campaign_name=job_name, current_phase="Staging BioEmu-2 inputs")

    staged_structures = []
    try:
        for index, source in enumerate(structures, start=1):
            digest = hashlib.sha256(source.read_bytes()).hexdigest()[:10]
            safe_name = "".join(char if char.isalnum() or char in "._-" else "_" for char in source.name)
            destination = input_root / f"complex_{index:03d}_{digest}_{safe_name}"
            shutil.copy2(source, destination)
            staged_structures.append(destination.relative_to(job.run_dir).as_posix())

        request = {
            "mode": mode,
            "structure_paths": staged_structures,
            "protein_sequences": sequences,
            "num_samples": num_samples,
            "lag_time_ns": lag_time_ns,
            "sampler": sampler,
            "physicality_filter": physicality_filter,
            "precision": precision,
            "matmul_precision": matmul_precision,
            "batch_size": batch_size,
            "num_steps": num_steps,
            "seed": seed,
            "compile_denoiser": compile_denoiser,
            "output_root": output_root.relative_to(job.run_dir).as_posix(),
        }
        request_path = input_root / "bioemu2_request.json"
        write_json(request_path, request)

        cache_dir = reference_root() / "bioemu2" / "huggingface"
        cache_dir.mkdir(parents=True, exist_ok=True)
        manifest = load_manifest("bioemu2")
        update_status(
            job.run_dir,
            "preparing",
            current_phase=(f"Sampling {len(structures)} starting complex(es)" if structures else "Preparing sequence-defined system"),
        )
        return_code = run_docker_job(
            job.run_dir,
            manifest,
            ["--studio-request", f"/work/{request_path.relative_to(job.run_dir).as_posix()}"],
            gpu_device=gpu_device,
            extra_mounts=[(cache_dir, "/cache/huggingface", False)],
            environment={
                "HF_HOME": "/cache/huggingface",
                "HF_HUB_CACHE": "/cache/huggingface/hub",
                "OMP_NUM_THREADS": str(cpu_cores),
            },
        )
        summary_path = output_root / "studio_summary.json"
        success = return_code == 0 and summary_path.is_file()
        if success:
            summary = json.loads(summary_path.read_text())
            total_systems = len(summary.get("records", []))
            finish_job(
                job.run_dir,
                True,
                {
                    "outputs": {
                        "ensemble_root": str(output_root),
                        "summary": str(summary_path),
                    },
                    "metrics": {
                        "mode": mode,
                        "systems": total_systems,
                        "samples_per_system": num_samples,
                        "total_sampled_structures": total_systems * num_samples,
                        "lag_time_ns": lag_time_ns if mode == "propagate" else None,
                    },
                },
            )
        else:
            finish_job(
                job.run_dir,
                False,
                {
                    "error": f"BioEmu-2 exited with code {return_code}.",
                    "outputs": {"ensemble_root": str(output_root)},
                },
            )
    except Exception as exc:
        with (job.run_dir / "stderr.log").open("a") as stderr:
            stderr.write(f"BioEmu-2 workflow error: {type(exc).__name__}: {exc}\n")
        finish_job(job.run_dir, False, {"error": f"{type(exc).__name__}: {exc}"})
    return job.run_dir
