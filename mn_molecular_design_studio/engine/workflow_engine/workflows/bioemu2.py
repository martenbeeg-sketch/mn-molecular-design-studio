"""BioEmu-2 structure sampling jobs for protein targets and designs."""
from __future__ import annotations

import hashlib
import json
import os
import secrets
import shutil
import tempfile
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


def _seed_bioemu2_unpaired_msa_cache(
    sequences: list[str],
    *,
    cache_dir: Path,
    reference_dir: Path,
) -> dict:
    """Seed BioEmu-2's raw ColabFold cache from Studio's exact-query A3Ms.

    BioEmu-2 keys the unpaired search by the complete ordered unique sequence
    list, while Studio stores each exact sequence alignment separately. Build
    the upstream cache record only when every sequence has a validated shared
    A3M. BioEmu-2 can then reuse it directly and still fetch a paired alignment
    when its paired-search cache has no matching entry.
    """
    from mn_molecular_design_studio.cofolding_msa import validated_text

    unique_sequences = list(dict.fromkeys("".join(str(sequence).split()).upper() for sequence in sequences))
    unique_sequences = [sequence for sequence in unique_sequences if sequence]
    if not unique_sequences:
        return {"seeded": False, "alignment_count": 0, "reason": "no protein sequences", "alignments": []}

    request = {
        "sequences": unique_sequences,
        "use_pairing": False,
        "pairing_strategy": "greedy",
        "use_env": True,
    }
    key = hashlib.sha256(json.dumps(request, sort_keys=True).encode("utf-8")).hexdigest()
    destination = cache_dir / f"{key}.json"
    if destination.is_file():
        return {
            "seeded": False,
            "cache_hit": True,
            "reason": "BioEmu-2 cache entry already exists",
            "bioemu_cache_path": str(destination),
            "alignment_count": len(unique_sequences),
            "alignments": [],
        }

    a3ms: list[str] = []
    alignment_rows: list[dict[str, str | int]] = []
    missing: list[str] = []
    for sequence in unique_sequences:
        sequence_digest = hashlib.sha256(sequence.encode("utf-8")).hexdigest()
        source = reference_dir / "msa_cache" / f"{sequence_digest}.a3m"
        if not source.is_file():
            missing.append(sequence_digest)
            continue
        try:
            a3m = validated_text(source.read_text(encoding="utf-8"), sequence)
        except (OSError, UnicodeError, ValueError):
            missing.append(sequence_digest)
            continue
        a3ms.append(a3m)
        alignment_rows.append(
            {
                "sequence_sha256": sequence_digest,
                "source_path": str(source),
                "bytes": source.stat().st_size,
            }
        )

    if missing:
        return {
            "seeded": False,
            "cache_hit": False,
            "reason": "one or more exact-query shared A3Ms are missing or invalid",
            "missing_sequence_sha256": missing,
            "alignment_count": len(alignment_rows),
            "alignments": alignment_rows,
        }

    result_digest = hashlib.sha256(json.dumps(a3ms).encode("utf-8")).hexdigest()
    payload = {"request": request, "a3ms": a3ms, "sha256": result_digest}
    cache_dir.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{key}.", suffix=".tmp", dir=cache_dir)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle)
        os.replace(temporary, destination)
    finally:
        Path(temporary).unlink(missing_ok=True)
    return {
        "seeded": True,
        "cache_hit": False,
        "reason": "seeded BioEmu-2 unpaired ColabFold cache from Studio exact-query A3Ms",
        "bioemu_cache_key": key,
        "bioemu_cache_path": str(destination),
        "alignment_count": len(alignment_rows),
        "alignments": alignment_rows,
        "paired_alignment": "separate BioEmu-2 paired search/cache; not synthesized from unpaired A3Ms",
    }


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
    source_run_dir: Path | None = None,
    candidates_jsonl: Path | None = None,
) -> Path:
    """Sample one or more starting complexes, or a sequence-defined system."""
    if source_run_dir is not None or candidates_jsonl is not None:
        if source_run_dir is None or candidates_jsonl is None or structure_paths:
            raise ValueError("Select a Run output or explicit structures, with both source paths supplied together.")
        from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
        from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path
        source_candidates = read_candidates(Path(candidates_jsonl))
        structure_paths = [resolve_stored_path(c["complex_pdb"], run_dir=Path(source_run_dir), must_exist=True)
                           for c in source_candidates if c.get("complex_pdb")]
        if len(structure_paths) != len(source_candidates):
            raise ValueError("Every candidate in the selected Run output must have a complex structure for BioEmu-2 propagation.")
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
        **({"source_run_dir": str(source_run_dir), "candidates_jsonl": str(candidates_jsonl)} if source_run_dir else {}),
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
            "output_dir": output_root.relative_to(job.run_dir).as_posix(),
        }
        request_path = input_root / "bioemu2_request.json"
        write_json(request_path, request)

        cache_dir = reference_root() / "bioemu2" / "huggingface"
        cache_dir.mkdir(parents=True, exist_ok=True)
        bioemu_msa_cache_dir = reference_root() / "bioemu2" / "colabfold_msa"
        bioemu_msa_cache_dir.mkdir(parents=True, exist_ok=True)
        msa_cache_summary = _seed_bioemu2_unpaired_msa_cache(
            sequences if mode in {"sample-equilibrium", "predict-structure"} else [],
            cache_dir=bioemu_msa_cache_dir,
            reference_dir=reference_root(),
        )
        msa_cache_manifest = input_root / "bioemu2_msa_cache_manifest.json"
        write_json(msa_cache_manifest, msa_cache_summary)
        manifest = load_manifest("bioemu2")
        update_status(
            job.run_dir,
            "preparing",
            current_phase=(
                f"Sampling {len(structures)} starting complex(es)"
                if structures
                else f"Preparing sequence-defined system; {msa_cache_summary.get('alignment_count', 0)} exact shared A3M(s) available to BioEmu-2"
            ),
        )
        return_code = run_docker_job(
            job.run_dir,
            manifest,
            ["--studio-request", f"/work/{request_path.relative_to(job.run_dir).as_posix()}"],
            gpu_device=gpu_device,
            extra_mounts=[
                (cache_dir, "/cache/huggingface", False),
                (bioemu_msa_cache_dir, "/cache/bioemu2_msa", False),
            ],
            environment={
                "HF_HOME": "/cache/huggingface",
                "HF_HUB_CACHE": "/cache/huggingface/hub",
                "BIOEMU2_MSA_CACHE": "/cache/bioemu2_msa",
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
                        "msa_cache_manifest": str(msa_cache_manifest),
                    },
                    "metrics": {
                        "mode": mode,
                        "systems": total_systems,
                        "samples_per_system": num_samples,
                        "total_sampled_structures": total_systems * num_samples,
                        "lag_time_ns": lag_time_ns if mode == "propagate" else None,
                        "shared_msa_cache_seeded": bool(msa_cache_summary.get("seeded")),
                        "shared_msa_alignment_count": msa_cache_summary.get("alignment_count", 0),
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
