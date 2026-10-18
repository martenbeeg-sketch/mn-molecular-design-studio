from __future__ import annotations

import json

from mn_molecular_design_studio import cofolding_msa
import os
import re
import secrets
import subprocess
import traceback
from pathlib import Path
from typing import Literal

from mn_molecular_design_studio.engine.workflow_engine.core.candidates import (
    STAGE_COMPLEX_REFOLDING,
    write_candidates,
)
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import (
    create_job,
    finish_job,
    update_status,
    write_json,
)
from mn_molecular_design_studio.engine.workflow_engine.runtime import reference_root
from mn_molecular_design_studio.engine.workflow_engine.workflows import target_msa as target_msa_workflow
from mn_molecular_design_studio.engine.workflow_engine.workflows.esmfold2_runtime import stage_esmfold2_msa


TOOL_NAME = "esmfold2_cofolding"
ESMFOLD2_IMAGE = os.getenv(
    "MN_MOLECULAR_DESIGN_STUDIO_ESMFOLD2_IMAGE",
    os.getenv("MN_COFOLDING_ESMFOLD2_IMAGE", "mn-esmfold2-opt:0.1.0-cu130"),
)
ESMFOLD2_RUNNER = "/usr/local/bin/mn-esmfold2-runner"
_STANDARD_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
_CHAIN_IDS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"


def parse_protein_sequences(text: str) -> list[dict[str, str]]:
    """Parse Cofolding-style colon-separated protein chains."""
    chains: list[dict[str, str]] = []
    for raw in str(text or "").split(":"):
        sequence = "".join(raw.split()).upper()
        if not sequence:
            continue
        invalid = sorted(set(sequence) - _STANDARD_AMINO_ACIDS)
        if invalid:
            raise ValueError(
                f"Chain {len(chains) + 1} contains unsupported amino-acid letters: {', '.join(invalid)}. "
                "Use the 20 standard amino acids."
            )
        chains.append({"id": _CHAIN_IDS[len(chains)] if len(chains) < len(_CHAIN_IDS) else "", "sequence": sequence})
    if not chains:
        raise ValueError("Enter at least one protein sequence. Separate chains with colons.")
    if len(chains) > len(_CHAIN_IDS):
        raise ValueError(f"ESMFold2 supports at most {len(_CHAIN_IDS)} input chains in this tool.")
    return chains


def _gpu_ids(gpu_device: object) -> list[str]:
    value = str(gpu_device or "0").strip().lower().removeprefix("device=")
    if value in {"all", "auto", "any"}:
        result = subprocess.run(
            ["nvidia-smi", "--query-gpu=index", "--format=csv,noheader"],
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
        ids = [line.strip() for line in result.stdout.splitlines() if line.strip()]
        if result.returncode or not ids:
            raise RuntimeError("Could not determine the allocated GPU IDs for ESMFold2.")
        return ids
    ids = [part.strip() for part in value.split(",") if part.strip()]
    if not ids or any(not item.isdigit() for item in ids):
        raise ValueError("gpu_device must be a GPU index or comma-separated indices such as '0' or '0,1'.")
    return list(dict.fromkeys(ids))


def _kit_config(gpu_id: str) -> str:
    result = subprocess.run(
        ["nvidia-smi", "-i", gpu_id, "--query-gpu=compute_cap", "--format=csv,noheader"],
        capture_output=True,
        text=True,
        timeout=10,
        check=False,
    )
    try:
        capability = float(result.stdout.splitlines()[0].strip()) if result.returncode == 0 else 0.0
    except (ValueError, IndexError):
        capability = 0.0
    if capability >= 12.0:
        return "rtx5090"
    if capability >= 8.9:
        return "rtx4090"
    raise RuntimeError(
        f"The optimized ESMFold2 image supports RTX 4090 (compute capability 8.9) and RTX 5090 (12.0); "
        f"GPU {gpu_id} reports {capability}."
    )


def _mean_plddt_from_cif(path: Path) -> tuple[float | None, dict[str, float]]:
    try:
        from Bio.PDB import MMCIFParser

        structure = MMCIFParser(QUIET=True).get_structure("esmfold2", str(path))
        by_chain: dict[str, list[float]] = {}
        for atom in structure.get_atoms():
            if atom.get_name().strip() != "CA":
                continue
            value = float(atom.get_bfactor())
            if value >= 0:
                chain_id = str(atom.get_parent().get_parent().id)
                by_chain.setdefault(chain_id, []).append(value)
        means = {
            chain_id: round(sum(values) / len(values), 3)
            for chain_id, values in by_chain.items()
            if values
        }
        all_values = [value for values in by_chain.values() for value in values]
        mean = round(sum(all_values) / len(all_values), 3) if all_values else None
        return mean, means
    except Exception:
        return None, {}


def _seed_from_name(path: Path) -> int | None:
    match = re.search(r"(?:seed[_-]?|__s)(\d+)", path.stem, flags=re.IGNORECASE)
    return int(match.group(1)) if match else None


def _read_prediction_rows(path: Path) -> dict[str, dict]:
    """Index the optimized kit's per-prediction confidence metrics by CIF filename."""
    rows: dict[str, dict] = {}
    if not path.is_file():
        return rows
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        try:
            row = json.loads(line)
        except json.JSONDecodeError:
            continue
        if isinstance(row, dict) and row.get("pred_file"):
            rows[Path(str(row["pred_file"])).name] = row
    return rows


def _confidence_summary_from_npz(path: Path, chain_ids: list[str]) -> tuple[float | None, dict[str, float]]:
    """Read ESMFold2's token pLDDT output and summarize it in the UI's 0–100 scale."""
    try:
        import numpy as np

        with np.load(path, allow_pickle=False) as confidence:
            values = np.asarray(confidence["plddt"], dtype=float).reshape(-1)
            asym_ids = np.asarray(confidence["asym_id"]).reshape(-1)
        if values.size and float(np.nanmax(values)) <= 1.5:
            values = values * 100.0
        finite = np.isfinite(values)
        mean = round(float(np.nanmean(values[finite])), 3) if finite.any() else None
        by_chain: dict[str, float] = {}
        if asym_ids.size == values.size:
            ordered_asym_ids = list(dict.fromkeys(asym_ids.tolist()))
            for index, asym_id in enumerate(ordered_asym_ids):
                chain_id = chain_ids[index] if index < len(chain_ids) else str(asym_id)
                mask = (asym_ids == asym_id) & finite
                if mask.any():
                    by_chain[chain_id] = round(float(np.mean(values[mask])), 3)
        return mean, by_chain
    except (ImportError, OSError, ValueError, KeyError, TypeError):
        return None, {}


def _prediction_command(
    *,
    run_dir: Path,
    input_file: Path,
    output_dir: Path,
    model_variant: str,
    msa_source: str,
    optimization_mode: str,
    seeds: list[int],
    num_loops: int,
    num_diffusion_samples: int,
    msa_max_depth: int,
    cpu_threads: int,
    gpu_ids: list[str],
) -> list[str]:
    if optimization_mode == "big" and len(gpu_ids) < 2:
        raise ValueError("ESMFold2 big mode requires GPUs 0 and 1.")
    if optimization_mode != "big" and not gpu_ids:
        raise ValueError("ESMFold2 requires one reserved GPU.")
    used_gpu_ids = gpu_ids[:2] if optimization_mode == "big" else gpu_ids[:1]
    model_variant_arg = "fast" if model_variant == "fast" else "full_msa" if msa_source != "single_sequence" else "full_nomsa"
    config = _kit_config(used_gpu_ids[0])
    cache_dir = Path(
        os.getenv(
            "MN_MOLECULAR_DESIGN_STUDIO_ESMFOLD2_CACHE_DIR",
            os.getenv("MN_COFOLDING_CACHE_DIR", str(reference_root() / "mn-cofolding" / "cache")),
        )
    ).expanduser().resolve()
    cache_dir.mkdir(parents=True, exist_ok=True)
    selected = ",".join(used_gpu_ids)
    docker = os.getenv("MN_MOLECULAR_DESIGN_STUDIO_DOCKER", os.getenv("MN_COFOLDING_DOCKER", "docker"))
    command = [
        docker,
        "run",
        "--rm",
        "--name",
        f"mds-{run_dir.name}-esmfold2",
        "--shm-size=16g",
        "--mount",
        f"type=bind,src={run_dir.resolve()},dst=/work",
        "--mount",
        f"type=bind,src={reference_root().resolve()},dst=/reference,readonly",
        "--mount",
        f"type=bind,src={cache_dir},dst=/cache",
        "-w",
        "/work",
        "--gpus",
        f'"device={selected}"',
        "--env",
        "HOME=/cache/home",
        "--env",
        "XDG_CACHE_HOME=/cache/xdg",
        "--env",
        "HF_HOME=/cache/huggingface/esmfold2",
        "--env",
        "MODEL_OPT_JIT_ROOT=/cache/esmfold2/jit",
        "--env",
        "HF_HUB_ENABLE_HF_TRANSFER=0",
    ]
    for variable in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
        "TF_NUM_INTRAOP_THREADS",
        "TF_NUM_INTEROP_THREADS",
    ):
        command.extend(["--env", f"{variable}={cpu_threads}"])
    command.extend(
        [
            ESMFOLD2_IMAGE,
            ESMFOLD2_RUNNER,
            "pred",
            "--config",
            config,
            "--variant",
            model_variant_arg,
            "--mode",
            optimization_mode,
            "--n_gpu",
            "2" if optimization_mode == "big" else "1",
            "--input",
            f"/work/{input_file.relative_to(run_dir).as_posix()}",
            "--out_dir",
            f"/work/{output_dir.relative_to(run_dir).as_posix()}",
            "--seeds",
            ",".join(str(seed) for seed in seeds),
            "--num_loops",
            str(num_loops),
            "--num_sampling_steps",
            "200",
            "--num_diffusion_samples",
            str(num_diffusion_samples),
            "--msa_max_depth",
            str(msa_max_depth),
            "--det",
            "1",
        ]
    )
    return command


def run_esmfold2_cofolding(
    protein_sequences: str,
    job_name: str = "esmfold2-cofold",
    model_variant: Literal["full", "fast"] = "full",
    msa_source: Literal["auto", "reference_then_colabfold", "single_sequence", "local_mmseqs_gpu", "boltz_msa_server"] = "auto",
    optimization_mode: Literal["off", "exact", "fast", "big"] = "off",
    seeds: list[int] | None = None,
    num_diffusion_samples: int = 5,
    num_loops: int = 10,
    msa_max_depth: int = 1024,
    cpu_threads: int = 4,
    gpu_device: str = "0",
) -> Path:
    """Run the optimized Biohub ESMFold2 image on one multichain sequence input."""
    if msa_source == "auto":
        msa_source = "single_sequence" if model_variant == "fast" else "reference_then_colabfold"
    requested_seeds = [int(seed) for seed in (seeds or [0])]
    used_seeds = {seed for seed in requested_seeds if seed > 0}
    seeds = []
    for seed in requested_seeds:
        if seed == 0:
            seed = secrets.randbelow(2**31 - 1) + 1
            while seed in used_seeds:
                seed = secrets.randbelow(2**31 - 1) + 1
            used_seeds.add(seed)
        seeds.append(seed)
    job_name = str(job_name or "esmfold2-cofold").strip() or "esmfold2-cofold"
    params = {
        "job_name": job_name,
        "campaign_name": job_name,
        "model_variant": model_variant,
        "msa_source": msa_source,
        "optimization_mode": optimization_mode,
        "seeds": seeds,
        "num_diffusion_samples": int(num_diffusion_samples),
        "num_loops": int(num_loops),
        "msa_max_depth": int(msa_max_depth),
        "cpu_threads": int(cpu_threads),
        "gpu_device": str(gpu_device),
        "container_image": ESMFOLD2_IMAGE,
    }
    job = create_job(
        "refolding-validation",
        job_type="complex_refolding",
        tool=TOOL_NAME,
        inputs={"protein_sequences": protein_sequences},
        params=params,
    )
    update_status(job.run_dir, "running", campaign_name=job_name)
    raw_dir = job.run_dir / "artifacts" / "raw" / TOOL_NAME
    input_dir = job.run_dir / "input"
    msa_dir = input_dir / "esmfold2_msas"
    output_dir = raw_dir / "output"
    raw_dir.mkdir(parents=True, exist_ok=True)
    input_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        if model_variant not in {"full", "fast"}:
            raise ValueError("Choose the ESMFold2 Full or ESMFold2 Fast model.")
        if msa_source not in {"reference_then_colabfold", "single_sequence", "local_mmseqs_gpu", "boltz_msa_server"}:
            raise ValueError("Choose single-sequence input, local MMseqs GPU, or the Boltz MSA server.")
        if model_variant == "fast" and msa_source != "single_sequence":
            raise ValueError("ESMFold2 Fast is sequence-only and cannot use an MSA source.")
        if optimization_mode not in {"off", "exact", "fast", "big"}:
            raise ValueError("Optimization mode must be off, exact, fast, or big.")
        if not 1 <= int(num_diffusion_samples) <= 5:
            raise ValueError("Diffusion samples per seed must be between 1 and 5.")
        if not 1 <= int(num_loops) <= 100:
            raise ValueError("Loop count must be between 1 and 100.")
        if not 1 <= int(msa_max_depth) <= 8192:
            raise ValueError("Maximum MSA rows per chain must be between 1 and 8,192.")
        if not 1 <= int(cpu_threads) <= 64:
            raise ValueError("CPU threads must be between 1 and 64.")
        if not seeds or len(set(seeds)) != len(seeds) or any(seed < 0 for seed in seeds):
            raise ValueError("Enter one or more distinct non-negative random seeds.")

        chains = parse_protein_sequences(protein_sequences)
        gpu_ids = _gpu_ids(gpu_device)
        if optimization_mode == "big" and len(gpu_ids) < 2:
            raise ValueError("ESMFold2 big mode requires GPUs 0 and 1; select gpu_device='0,1'.")
        if optimization_mode == "big" and gpu_ids[:2] != ["0", "1"]:
            raise ValueError("ESMFold2 big mode is configured to row-shard across GPUs 0 and 1 in that order.")

        with (job.run_dir / "stdout.log").open("a", encoding="utf-8") as stdout:
            stdout.write(
                f"Preparing {len(chains)} ESMFold2 protein chain(s); model={model_variant}; "
                f"MSA={msa_source}; optimization={optimization_mode}.\n"
            )
            stdout.flush()

        if msa_source != "single_sequence":
            msa_dir.mkdir(parents=True, exist_ok=True)
            for chain in chains:
                supplied = ''
                if msa_source == "local_mmseqs_gpu" and cofolding_msa.find_reference(chain["sequence"]) is None:
                    target_msa_workflow.ensure_boltz_msa_for_sequence(
                        job.run_dir, chain["sequence"], f"esmfold2_chain_{chain['id']}",
                        raw_subdir=f"{TOOL_NAME}/msa", gpu_device=gpu_ids[0], msa_source="alphafast_mmseqs_gpu")
                    source_msa, _ = target_msa_workflow.boltz_msa_paths(chain["sequence"])
                    supplied = source_msa.read_text()
                staged = cofolding_msa.prepare_chain(job.run_dir, chain["sequence"], chain["id"],
                                                     depth=int(msa_max_depth), supplied=supplied, supplied_source="local_mmseqs_gpu")
                destination = msa_dir / f"chain_{chain['id']}.a3m"
                destination.write_text(staged.read_text())
                chain["msa"] = destination.relative_to(input_dir).as_posix()

        input_payload = {
            "id": re.sub(r"[^A-Za-z0-9_.-]+", "_", params["job_name"]).strip("._")[:120] or "esmfold2_cofold",
            "sequences": [{"type": "protein", **chain} for chain in chains],
            "seeds": seeds,
            "num_diffusion_samples": int(num_diffusion_samples),
        }
        input_file = input_dir / "esmfold2-cofold-input.json"
        input_file.write_text(json.dumps(input_payload, indent=2) + "\n", encoding="utf-8")
        (raw_dir / "input_sequences.fasta").write_text(
            "".join(f">chain_{chain['id']}\n{chain['sequence']}\n" for chain in chains),
            encoding="utf-8",
        )
        write_json(
            raw_dir / "input_mapping.json",
            {
                "chains": [{"id": chain["id"], "length": len(chain["sequence"])} for chain in chains],
                "model_variant": model_variant,
                "msa_source": msa_source,
                "msa_manifest": "artifacts/raw/cofolding_msa/manifest.json",
                "msa_files_by_chain": {chain["id"]: str((input_dir / chain["msa"]).relative_to(job.run_dir)) for chain in chains if chain.get("msa")},
                "optimization_mode": optimization_mode,
                "seeds": seeds,
            },
        )

        command = _prediction_command(
            run_dir=job.run_dir,
            input_file=input_file,
            output_dir=output_dir,
            model_variant=model_variant,
            msa_source=msa_source,
            optimization_mode=optimization_mode,
            seeds=seeds,
            num_loops=int(num_loops),
            num_diffusion_samples=int(num_diffusion_samples),
            msa_max_depth=int(msa_max_depth),
            cpu_threads=int(cpu_threads),
            gpu_ids=gpu_ids,
        )
        with (job.run_dir / "stdout.log").open("a", encoding="utf-8") as stdout, (job.run_dir / "stderr.log").open(
            "a", encoding="utf-8"
        ) as stderr:
            stdout.write("$ " + " ".join(command) + "\n")
            stdout.flush()
            completed = subprocess.run(command, stdout=stdout, stderr=stderr, check=False)

        prediction_log = output_dir / "esmfold2.pred.log"
        structures_dir = output_dir / "cif_all"
        structure_paths = sorted(structures_dir.glob("*.cif")) if structures_dir.is_dir() else []
        if not structure_paths:
            structure_paths = sorted(output_dir.rglob("*.cif"))
        if completed.returncode != 0:
            raise RuntimeError(
                f"Optimized ESMFold2 exited with code {completed.returncode}; produced {len(structure_paths)} structure(s)."
            )
        if not structure_paths:
            raise RuntimeError("Optimized ESMFold2 completed without producing a CIF structure.")

        rows_path = output_dir / "pred_rows.jsonl"
        prediction_rows = _read_prediction_rows(rows_path)
        candidates: list[dict] = []
        for index, structure_path in enumerate(structure_paths, start=1):
            plddt_mean, plddt_by_chain = _mean_plddt_from_cif(structure_path)
            confidence_path = structure_path.with_name(f"{structure_path.stem}_pae.npz")
            prediction = prediction_rows.get(structure_path.name, {})
            if confidence_path.is_file():
                confidence_mean, confidence_by_chain = _confidence_summary_from_npz(
                    confidence_path,
                    [chain["id"] for chain in chains],
                )
                plddt_mean = confidence_mean if confidence_mean is not None else plddt_mean
                plddt_by_chain = confidence_by_chain or plddt_by_chain
            seed = prediction.get("seed")
            if seed is None:
                seed = _seed_from_name(structure_path)
            sample_index = prediction.get("sample_index")
            candidate_id = f"esmfold2_cofold_{index:04d}"
            candidates.append(
                {
                    "candidate_id": candidate_id,
                    "stage": STAGE_COMPLEX_REFOLDING,
                    "source_tool": TOOL_NAME,
                    "complex_pdb": str(structure_path.relative_to(job.run_dir)),
                    "metrics": {
                        "plddt_mean": plddt_mean,
                        "plddt_by_chain": plddt_by_chain,
                        "ptm": prediction.get("ptm"),
                        "iptm": prediction.get("iptm"),
                        "prediction_index": index,
                        "sample_index": sample_index,
                        "seed": seed,
                        "model_variant": model_variant,
                        "msa_source": msa_source,
                        "optimization_mode": optimization_mode,
                    },
                    "raw_metadata": {
                        "structure_file": str(structure_path.relative_to(job.run_dir)),
                        "confidence_file": (
                            str(confidence_path.relative_to(job.run_dir)) if confidence_path.is_file() else None
                        ),
                        "prediction_rows_file": (
                            str(rows_path.relative_to(job.run_dir)) if rows_path.is_file() else None
                        ),
                        "chain_sequences": {chain["id"]: chain["sequence"] for chain in chains},
                        "model_variant": model_variant,
                        "optimization_mode": optimization_mode,
                        "container_image": ESMFOLD2_IMAGE,
                    },
                }
            )
        normalized = write_candidates(job.run_dir, TOOL_NAME, candidates)
        artifact_paths = [
            str(input_file.relative_to(job.run_dir)),
            str((raw_dir / "input_mapping.json").relative_to(job.run_dir)),
            str((raw_dir / "input_sequences.fasta").relative_to(job.run_dir)),
            *[str(path.relative_to(job.run_dir)) for path in structure_paths],
            *[
                str(path.relative_to(job.run_dir))
                for path in (item.with_name(f"{item.stem}_pae.npz") for item in structure_paths)
                if path.is_file()
            ],
        ]
        if rows_path.is_file():
            artifact_paths.append(str(rows_path.relative_to(job.run_dir)))
        if prediction_log.is_file():
            artifact_paths.append(str(prediction_log.relative_to(job.run_dir)))
        summary_path = raw_dir / "prediction_summary.json"
        artifact_paths.append(str(summary_path.relative_to(job.run_dir)))
        write_json(summary_path, {"candidate_count": len(normalized), "artifacts": artifact_paths})
        finish_job(
            job.run_dir,
            True,
            {
                "outputs": {"artifacts": artifact_paths, "candidates": normalized},
                "metrics": {
                    "candidate_count": len(normalized),
                    "chain_count": len(chains),
                    "model_variant": model_variant,
                    "msa_source": msa_source,
                    "optimization_mode": optimization_mode,
                    "container_image": ESMFOLD2_IMAGE,
                    "seeds": seeds,
                    "num_diffusion_samples": int(num_diffusion_samples),
                    "num_loops": int(num_loops),
                    "msa_max_depth": int(msa_max_depth),
                },
                "downstream_artifacts": {
                    "candidates_jsonl": "artifacts/normalized_candidates/candidates.jsonl",
                    "campaign_result": "artifacts/normalized_candidates/campaign_result.json",
                },
            },
        )
    except Exception as exc:
        with (job.run_dir / "stderr.log").open("a", encoding="utf-8") as stderr:
            stderr.write("\n" + traceback.format_exc() + "\n")
        finish_job(
            job.run_dir,
            False,
            {"outputs": {}, "metrics": {"error": str(exc), "error_type": type(exc).__name__}},
        )
    return job.run_dir
