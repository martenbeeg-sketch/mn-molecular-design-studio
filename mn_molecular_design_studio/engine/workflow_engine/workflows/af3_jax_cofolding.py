from __future__ import annotations

import json

from mn_molecular_design_studio import cofolding_msa
import hashlib
import math
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


TOOL_NAME = "af3_jax_cofolding"
AF3_IMAGE = os.getenv(
    "MN_MOLECULAR_DESIGN_STUDIO_AF3_IMAGE",
    os.getenv("MN_COFOLDING_AF3_IMAGE", "mn-cofolding-af3-memory:3.1.14-cu13"),
)
AF3_RUNNER = "/usr/local/bin/cofolding-af3-runner.py"
_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWYBXZJUO")
_DNA = frozenset("ACGTNRYKMSWBDHV")
_RNA = frozenset("ACGUNRYKMSWBDHV")
_CHAIN_IDS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"

# Models supported by MN Cofolding's shared AF3-family JAX runner. Native
# model containers remain available through their optimized-prediction tiles.
MODEL_PROFILES: dict[str, dict] = {
    "openbind0": {"label": "AF3-family runner · OpenBind-0 weights", "modes": ("off", "fast", "big"), "two_gpu_big": False, "terms": False},
    "openfold3": {"label": "AlphaFold3 JAX · OpenFold3 P2 weights", "modes": ("off", "exact", "fast", "big"), "two_gpu_big": False, "terms": False},
    "boltz2": {"label": "AF3-family runner · Boltz-2 weights", "modes": ("off", "fast", "big"), "two_gpu_big": False, "terms": False, "supports_affinity": True},
    "protenix2": {"label": "AF3-family runner · Protenix V2 weights", "modes": ("off", "fast", "big"), "two_gpu_big": False, "terms": False},
    "rosettafold3": {"label": "AF3-family runner · RosettaFold3 weights", "modes": ("off", "fast", "big"), "two_gpu_big": False, "terms": False},
    "intellifold2": {"label": "AF3-family runner · IntelliFold 2 weights", "modes": ("off", "fast", "big"), "two_gpu_big": False, "terms": False},
    "opendde": {"label": "AF3-JAX · OpenDDE weights", "modes": ("off", "fast", "big"), "two_gpu_big": False, "terms": False},
    "alphafold3": {"label": "AlphaFold3 JAX · official weights", "modes": ("off", "exact", "fast", "big"), "two_gpu_big": True, "terms": True},
}


def stage_batch_model_patch(run_dir: Path) -> Path:
    """Expose the pinned preview runner's model families to its batch adapter."""
    path = run_dir / "artifacts/raw/optimized_prediction/af3_family_batch_patch.py"
    path.parent.mkdir(parents=True, exist_ok=True)
    allowed = "{" + ", ".join(repr(name) for name in sorted(MODEL_PROFILES)) + "}"
    path.write_text(
        "from pathlib import Path\n"
        "path = Path('/opt/mn_batch_api/adapters.py')\n"
        "source = path.read_text()\n"
        "anchor = 'if model not in {\"alphafold3\", \"openfold3\"}:'\n"
        "if source.count(anchor) != 1:\n"
        "    raise RuntimeError('Pinned AF3 batch adapter model guard changed')\n"
        f"source = source.replace(anchor, {('if model not in ' + allowed + ':')!r})\n"
        "path.write_text(source)\n"
    )
    return path


def _reference_dir() -> Path:
    for variable in ("MN_MOLECULAR_DESIGN_STUDIO_REFERENCE_DIR", "MN_COFOLDING_REFERENCE_DIR"):
        configured = os.getenv(variable, "").strip()
        if configured:
            return Path(configured).expanduser().resolve()
    return reference_root()


def _split_colon_separated(value: str) -> list[str]:
    return [part.strip() for part in re.split(r":+", str(value or "").replace("\r", "").replace("\n", "")) if part.strip()]


def _clean_sequence(value: str, alphabet: frozenset[str], molecule: str) -> str:
    sequence = "".join(str(value or "").split()).upper()
    if not sequence:
        raise ValueError(f"{molecule} sequences cannot be empty.")
    invalid = sorted(set(sequence) - alphabet)
    if invalid:
        raise ValueError(f"Invalid character(s) in {molecule} sequence: {''.join(invalid)}")
    return sequence


def parse_molecule_fields(
    protein_sequences: str,
    dna_sequences: str,
    rna_sequences: str,
    ligand_ccd_codes: str,
    ligand_smiles: str,
) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]], list[str], list[str]]:
    proteins = [
        {"id": "", "sequence": _clean_sequence(value, _AMINO_ACIDS, "Protein")}
        for value in _split_colon_separated(protein_sequences)
    ]
    dna = [
        {"id": "", "sequence": _clean_sequence(value, _DNA, "DNA")}
        for value in _split_colon_separated(dna_sequences)
    ]
    rna = [
        {"id": "", "sequence": _clean_sequence(value, _RNA, "RNA")}
        for value in _split_colon_separated(rna_sequences)
    ]
    ccd_codes: list[str] = []
    for code in _split_colon_separated(ligand_ccd_codes):
        if not re.fullmatch(r"[A-Za-z0-9]{1,8}", code):
            raise ValueError(f"Invalid CCD code {code!r}; use one or more 1–8 character component IDs.")
        ccd_codes.append(code.upper())
    smiles = [line.strip() for line in str(ligand_smiles or "").splitlines() if line.strip()]
    if any(len(value) > 50_000 for value in smiles):
        raise ValueError("A ligand SMILES string exceeds 50,000 characters.")

    total_entities = len(proteins) + len(dna) + len(rna) + len(ccd_codes) + len(smiles)
    if not total_entities:
        raise ValueError("Add at least one protein, DNA, RNA, CCD ligand, or SMILES ligand.")
    if total_entities > len(_CHAIN_IDS):
        raise ValueError(f"AlphaFold3 supports at most {len(_CHAIN_IDS)} input chains or ligands in this tool.")
    for index, entity in enumerate([*proteins, *dna, *rna]):
        entity["id"] = _CHAIN_IDS[index]
    return proteins, dna, rna, ccd_codes, smiles


def _parse_seeds(seeds: list[int] | None) -> list[int]:
    requested = [int(seed) for seed in (seeds or [0])]
    values = list(dict.fromkeys(seed for seed in requested if seed != 0))
    if any(seed == 0 for seed in requested):
        used = set(values)
        for _ in (seed for seed in requested if seed == 0):
            random_seed = secrets.randbelow(2**31 - 1) + 1
            while random_seed in used:
                random_seed = secrets.randbelow(2**31 - 1) + 1
            values.append(random_seed)
            used.add(random_seed)
    if not values or any(seed < 0 or seed > 2**31 - 1 for seed in values):
        raise ValueError("Provide one or more seeds from 0 to 2,147,483,647; 0 means random.")
    if len(values) > 32:
        raise ValueError("AlphaFold3 accepts at most 32 random seeds per job.")
    return values


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
            raise RuntimeError("Could not determine the allocated GPU IDs for MN Cofolding.")
        return ids
    ids = [part.strip() for part in value.split(",") if part.strip()]
    if not ids or any(not item.isdigit() for item in ids):
        raise ValueError("gpu_device must be a GPU index or comma-separated indices such as '0' or '0,1'.")
    return list(dict.fromkeys(ids))


def _gpu_tuning(gpu_id: str) -> tuple[str, str]:
    try:
        result = subprocess.run(
            ["nvidia-smi", "-i", gpu_id, "--query-gpu=compute_cap", "--format=csv,noheader"],
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
        capability = float(result.stdout.splitlines()[0].strip()) if result.returncode == 0 else 0.0
    except (OSError, ValueError, IndexError, subprocess.TimeoutExpired):
        capability = 0.0
    if capability < 8.0:
        return "xla", "--xla_disable_hlo_passes=custom-kernel-fusion-rewriter"
    if capability < 9.0:
        return "xla", "--xla_gpu_enable_triton_gemm=false"
    return "triton", "--xla_gpu_enable_triton_gemm=false"


def _cap_a3m(text: str, maximum_rows: int) -> str:
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    for raw_line in text.replace("\x00", "").replace("\r\n", "\n").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks)))
            header, chunks = line, []
        elif header is not None:
            chunks.append(line)
        else:
            raise ValueError("A3M content appeared before its first FASTA header.")
    if header is not None:
        records.append((header, "".join(chunks)))
    if not records:
        raise ValueError("The generated protein MSA contains no FASTA records.")
    return "\n".join(part for row in records[:maximum_rows] for part in row) + "\n"


def _find_reference_msa(sequence: str) -> Path | None:
    """Reuse exact-sequence A3Ms from the configured reference_files trees."""
    reference_dir = _reference_dir()
    for repository in (
        reference_dir / "msa_cache",
        reference_dir / "boltz_models" / "msa_repository",
    ):
        cached, _source = target_msa_workflow.find_cached_msa_for_sequence(
            sequence, msa_repository_dir=repository,
        )
        if cached is not None and cached.is_file():
            return cached
    return None


def _run_colabfold_server_msa(run_dir: Path, sequence: str, chain_id: str) -> Path:
    """Acquire and stage through the common exact-query reference cache."""
    return cofolding_msa.prepare_chain(run_dir, sequence, chain_id, depth=1_000_000)


def _mean_plddt_from_cif(path: Path) -> float | None:
    try:
        from Bio.PDB import MMCIFParser

        structure = MMCIFParser(QUIET=True).get_structure("af3_jax", str(path))
        values = [
            float(atom.get_bfactor())
            for atom in structure.get_atoms()
            if atom.get_name().strip() == "CA" and float(atom.get_bfactor()) >= 0
        ]
        return round(sum(values) / len(values), 3) if values else None
    except Exception:
        return None


def _summary_confidence(path: Path) -> dict:
    stem = path.stem[:-6] if path.stem.endswith("_model") else path.stem
    summary_path = path.with_name(f"{stem}_summary_confidences.json")
    if not summary_path.is_file():
        return {}
    try:
        value = json.loads(summary_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    return value if isinstance(value, dict) else {}


def _confidence_data(path: Path) -> tuple[dict, Path | None]:
    stem = path.stem[:-6] if path.stem.endswith("_model") else path.stem
    confidence_path = path.with_name(f"{stem}_confidences.json")
    if not confidence_path.is_file():
        confidence: dict = {}
    else:
        try:
            value = json.loads(confidence_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            value = {}
        confidence = value if isinstance(value, dict) else {}

    # Boltz may write its affinity head separately from the structure confidence
    # file. Merge that optional sidecar so the result metrics can expose it too.
    affinity_stems = list(dict.fromkeys((stem, re.sub(r"_model_\d+$", "", stem))))
    for affinity_stem in affinity_stems:
        affinity_path = path.with_name(f"affinity_{affinity_stem}.json")
        if not affinity_path.is_file():
            continue
        try:
            affinity_value = json.loads(affinity_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        if isinstance(affinity_value, dict):
            confidence.update(affinity_value)
    return confidence, confidence_path if confidence_path.is_file() else None


def _numeric_metric(value: object) -> float | None:
    if isinstance(value, (list, tuple)) and len(value) == 1:
        value = value[0]
    if isinstance(value, bool):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return round(number, 5) if math.isfinite(number) else None


def _prediction_command(
    *,
    run_dir: Path,
    input_file: Path,
    output_dir: Path,
    model_name: str,
    optimization_mode: str,
    num_diffusion_samples: int,
    num_recycles: int,
    cpu_cores: int,
    gpu_ids: list[str],
) -> list[str]:
    model_profile = MODEL_PROFILES[model_name]
    uses_two_gpu_big = optimization_mode == "big" and bool(model_profile["two_gpu_big"])
    if uses_two_gpu_big:
        if gpu_ids[:2] != ["0", "1"]:
            raise ValueError(f"{model_profile['label']} big mode needs GPUs 0 and 1 in that order.")
        used_gpu_ids = gpu_ids[:2]
    else:
        if not gpu_ids:
            raise ValueError(f"{model_profile['label']} requires one reserved GPU.")
        used_gpu_ids = gpu_ids[:1]

    cache_dir = Path(
        os.getenv(
            "MN_MOLECULAR_DESIGN_STUDIO_AF3_CACHE_DIR",
            os.getenv("MN_COFOLDING_CACHE_DIR", str(_reference_dir() / "mn-cofolding" / "cache")),
        )
    ).expanduser().resolve()
    cache_dir.mkdir(parents=True, exist_ok=True)
    msa_cache_dir = Path(
        os.getenv("MN_SHARED_MSA_CACHE_DIR", str(_reference_dir() / "msa_cache"))
    ).expanduser().resolve()
    msa_cache_dir.mkdir(parents=True, exist_ok=True)
    docker = os.getenv("MN_MOLECULAR_DESIGN_STUDIO_DOCKER", os.getenv("MN_COFOLDING_DOCKER", "docker"))
    attention, xla_flags = _gpu_tuning(used_gpu_ids[0])
    selected = ",".join(used_gpu_ids)
    command = [
        docker,
        "run",
        "--rm",
        "--name",
        f"mds-{run_dir.name}-af3-jax",
        "--shm-size=16g",
        "--mount",
        f"type=bind,src={run_dir.resolve()},dst=/work",
        "--mount",
        f"type=bind,src={_reference_dir()},dst=/reference,readonly",
        "--mount",
        f"type=bind,src={cache_dir},dst=/cache",
        "--mount",
        f"type=bind,src={msa_cache_dir},dst=/msa_cache",
        "-w",
        "/work",
        "--gpus",
        f'"device={selected}"',
        "--env",
        "HOME=/cache/home",
        "--env",
        "XDG_CACHE_HOME=/cache/xdg",
        "--env",
        "HF_HOME=/cache/huggingface",
        "--env",
        "AF3_WEIGHTS_DIR=/cache/alphafold3/weights",
        "--env",
        "AF3_CCD_CACHE_DIR=/cache/alphafold3/ccd",
        "--env",
        "AF3_MSA_CACHE_DIR=/msa_cache/mn-cofolding/alphafold3",
        "--env",
        "AF3_SHARED_MSA_CACHE_DIR=/msa_cache",
        "--env",
        f"MN_COFOLDING_AF3_OPT_MODE={optimization_mode}",
        "--env",
        f"MN_COFOLDING_AF3_N_GPU={2 if uses_two_gpu_big else 1}",
    ]
    for variable in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
        "TF_NUM_INTRAOP_THREADS",
        "TF_NUM_INTEROP_THREADS",
    ):
        command.extend(["--env", f"{variable}={cpu_cores}"])
    if xla_flags:
        command.extend(["--env", f"XLA_FLAGS={xla_flags}"])
    command.extend(
        [
            AF3_IMAGE,
            "python",
            AF3_RUNNER,
            f"--json_path=/work/{input_file.relative_to(run_dir).as_posix()}",
            f"--model={model_name}",
            "--norun_data_pipeline",
            f"--output_dir=/work/{output_dir.relative_to(run_dir).as_posix()}",
            "--force_output_dir",
            "--cache_dir=/cache/alphafold3",
            f"--num_recycles={num_recycles}",
            f"--num_diffusion_samples={num_diffusion_samples}",
            f"--flash_attention_implementation={attention}",
        ]
    )
    if model_name == "alphafold3":
        command.append("--model_dir=/reference/alphafold3")
    return command


def run_af3_jax_cofolding(
    protein_sequences: str = "",
    dna_sequences: str = "",
    rna_sequences: str = "",
    ligand_ccd_codes: str = "",
    ligand_smiles: str = "",
    model_name: Literal[
        "openbind0", "openfold3", "boltz2", "protenix2", "rosettafold3",
        "intellifold2", "opendde", "alphafold3",
    ] = "alphafold3",
    job_name: str = "af3-jax-cofold",
    msa_source: Literal["single_sequence", "local_mmseqs_gpu", "reference_then_colabfold"] = "reference_then_colabfold",
    optimization_mode: Literal["off", "exact", "fast", "big"] = "off",
    seeds: list[int] | None = None,
    num_diffusion_samples: int = 5,
    num_recycles: int = 10,
    msa_max_depth: int = 1024,
    cpu_cores: int = 4,
    gpu_device: str = "0",
    accept_alphafold3_terms: bool = False,
    structure_templates: list[dict] | None = None,
) -> Path:
    """Run an AF3-family model through MN Cofolding's shared JAX runner."""
    if model_name not in MODEL_PROFILES:
        raise ValueError(f"Unsupported MN Cofolding model: {model_name}.")
    model_profile = MODEL_PROFILES[model_name]
    from mn_molecular_design_studio.structure_templates import validate_uploads, stage_templates
    structure_templates = validate_uploads(structure_templates, supported=True)
    seeds = _parse_seeds(seeds)
    job_name = str(job_name or "af3-jax-cofold").strip() or "af3-jax-cofold"
    params = {
        "model_name": model_name,
        "model_label": model_profile["label"],
        "job_name": job_name,
        "msa_source": msa_source,
        "optimization_mode": optimization_mode,
        "seeds": seeds,
        "num_diffusion_samples": int(num_diffusion_samples),
        "num_recycles": int(num_recycles),
        "msa_max_depth": int(msa_max_depth),
        "cpu_cores": int(cpu_cores),
        "gpu_device": str(gpu_device),
        "accept_alphafold3_terms": bool(accept_alphafold3_terms),
        "container_image": AF3_IMAGE,
    }
    molecule_inputs = {
        "protein_sequences": protein_sequences,
        "dna_sequences": dna_sequences,
        "rna_sequences": rna_sequences,
        "ligand_ccd_codes": ligand_ccd_codes,
        "ligand_smiles": ligand_smiles,
    }
    job = create_job(
        "refolding-validation",
        job_type="complex_refolding",
        tool=TOOL_NAME,
        inputs=molecule_inputs,
        params=params,
    )
    update_status(job.run_dir, "running", campaign_name=job_name)
    raw_dir = job.run_dir / "artifacts" / "raw" / TOOL_NAME
    input_dir = job.run_dir / "input"
    msa_dir = input_dir / "af3_msas"
    output_dir = raw_dir / "output"
    raw_dir.mkdir(parents=True, exist_ok=True)
    input_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        if model_profile["terms"] and not accept_alphafold3_terms:
            raise ValueError("Review and accept AlphaFold3's published weight terms before using the official model.")
        if msa_source not in {"single_sequence", "local_mmseqs_gpu", "reference_then_colabfold"}:
            raise ValueError("Choose single-sequence input, local MMseqs GPU, or reference MSAs with ColabFold server fallback.")
        if optimization_mode not in model_profile["modes"]:
            raise ValueError(
                f"{model_profile['label']} supports these inference modes: {', '.join(model_profile['modes'])}."
            )
        if optimization_mode != "off" and AF3_IMAGE == "mn-cofolding-af3:3.1.14-cu13":
            raise ValueError("Inference optimization requires the optimized mn-cofolding-af3-memory image.")
        if optimization_mode in {"exact", "fast"} and AF3_IMAGE != "mn-cofolding-af3-memory:3.1.14-cu13":
            raise ValueError(
                f"{model_profile['label']} {optimization_mode} mode requires mn-cofolding-af3-memory:3.1.14-cu13, "
                "the source-pinned optimized image."
            )
        if not 1 <= int(num_diffusion_samples) <= 5:
            raise ValueError("Diffusion samples per seed must be between 1 and 5.")
        if not 1 <= int(num_recycles) <= 100:
            raise ValueError("Recycle count must be between 1 and 100.")
        if not 1 <= int(msa_max_depth) <= 65536:
            raise ValueError("Maximum MSA rows per chain must be between 1 and 65,536.")
        if not 1 <= int(cpu_cores) <= 256:
            raise ValueError("CPU cores to reserve must be between 1 and 256.")

        proteins, dna, rna, ccd_codes, smiles = parse_molecule_fields(
            protein_sequences,
            dna_sequences,
            rna_sequences,
            ligand_ccd_codes,
            ligand_smiles,
        )
        gpu_ids = _gpu_ids(gpu_device)
        uses_two_gpu_big = optimization_mode == "big" and bool(model_profile["two_gpu_big"])
        if uses_two_gpu_big and gpu_ids[:2] != ["0", "1"]:
            raise ValueError(f"{model_profile['label']} big mode needs GPUs 0 and 1, in that order; reserve both in the scheduler.")
        if not gpu_ids:
            raise ValueError(f"{model_profile['label']} requires one reserved GPU.")

        with (job.run_dir / "stdout.log").open("a", encoding="utf-8") as stdout:
            stdout.write(
                f"Preparing {model_profile['label']}; proteins={len(proteins)}; DNA={len(dna)}; "
                f"RNA={len(rna)}; ligands={len(ccd_codes) + len(smiles)}; MSA={msa_source}; "
                f"optimization={optimization_mode}.\n"
            )
            stdout.flush()

        msa_gpu_device = gpu_ids[0]
        msa_source_by_chain: dict[str, str] = {}
        if msa_source != "single_sequence":
            msa_dir.mkdir(parents=True, exist_ok=True)
            for chain in proteins:
                supplied = ''
                if msa_source == "local_mmseqs_gpu" and cofolding_msa.find_reference(chain["sequence"]) is None:
                    target_msa_workflow.ensure_boltz_msa_for_sequence(
                        job.run_dir, chain["sequence"], f"af3_jax_chain_{chain['id']}",
                        raw_subdir=f"{TOOL_NAME}/msa", gpu_device=msa_gpu_device,
                        msa_source="alphafast_mmseqs_gpu")
                    source_msa, _ = target_msa_workflow.boltz_msa_paths(chain["sequence"])
                    supplied = source_msa.read_text()
                staged_msa = cofolding_msa.prepare_chain(job.run_dir, chain["sequence"], chain["id"],
                                                        depth=int(msa_max_depth), supplied=supplied, supplied_source="local_mmseqs_gpu")
                # Keep the established AF3 input mapping path, with the exact capped input.
                destination = msa_dir / f"chain_{chain['id']}.a3m"
                destination.write_text(staged_msa.read_text())
                chain["msa"] = destination.read_text()
                msa_source_by_chain[chain["id"]] = next(
                    row["source"] for row in json.loads((job.run_dir / "artifacts/raw/cofolding_msa/manifest.json").read_text())["alignments"]
                    if row["chain_id"] == chain["id"])

        entities: list[dict] = []
        for chain in proteins:
            protein = {
                "id": [chain["id"]],
                "sequence": chain["sequence"],
                "modifications": [],
                "templates": [],
            }
            if msa_source == "single_sequence":
                protein["unpairedMsa"] = f">query\n{chain['sequence']}\n"
                protein["pairedMsa"] = ""
            else:
                protein["unpairedMsa"] = chain["msa"]
                protein["pairedMsa"] = ""
            entities.append({"protein": protein})
        for kind, chains in (("dna", dna), ("rna", rna)):
            entities.extend(
                {kind: {"id": [chain["id"]], "sequence": chain["sequence"]}}
                for chain in chains
            )
        next_id = len(proteins) + len(dna) + len(rna)
        for code in ccd_codes:
            entities.append({"ligand": {"id": [_CHAIN_IDS[next_id]], "ccdCodes": [code]}})
            next_id += 1
        for value in smiles:
            entities.append({"ligand": {"id": [_CHAIN_IDS[next_id]], "smiles": value}})
            next_id += 1

        safe_name = re.sub(r"[^A-Za-z0-9_-]+", "_", job_name).strip("_")[:80] or "af3_jax_cofold"
        input_payload = {
            "name": safe_name,
            "sequences": entities,
            "modelSeeds": seeds,
            "dialect": "alphafold3",
            "version": 1,
        }
        input_file = input_dir / "af3-jax-cofold-input.json"
        input_file.write_text(json.dumps(input_payload, indent=2) + "\n", encoding="utf-8")
        if structure_templates:
            stage_templates(job.run_dir, input_file, item_id="preview", profile={"input_kind": "af3_json"},
                            proteins=proteins, uploads=structure_templates)
        (raw_dir / "input_sequences.fasta").write_text(
            "".join(f">chain_{chain['id']}\n{chain['sequence']}\n" for chain in proteins),
            encoding="utf-8",
        )
        write_json(
            raw_dir / "input_mapping.json",
            {
                "protein_chains": [{"id": chain["id"], "length": len(chain["sequence"])} for chain in proteins],
                "dna_chains": [{"id": chain["id"], "length": len(chain["sequence"])} for chain in dna],
                "rna_chains": [{"id": chain["id"], "length": len(chain["sequence"])} for chain in rna],
                "ligand_ccd_codes": ccd_codes,
                "ligand_smiles_count": len(smiles),
                "msa_source": msa_source,
                "msa_source_by_chain": msa_source_by_chain,
                "msa_files_by_chain": {
                    chain["id"]: str((msa_dir / f"chain_{chain['id']}.a3m").relative_to(job.run_dir))
                    for chain in proteins
                    if (msa_dir / f"chain_{chain['id']}.a3m").is_file()
                },
                "msa_max_depth": int(msa_max_depth),
                "model_name": model_name,
                "model_label": model_profile["label"],
                "optimization_mode": optimization_mode,
                "seeds": seeds,
            },
        )

        command = _prediction_command(
            run_dir=job.run_dir,
            input_file=input_file,
            output_dir=output_dir,
            model_name=model_name,
            optimization_mode=optimization_mode,
            num_diffusion_samples=int(num_diffusion_samples),
            num_recycles=int(num_recycles),
            cpu_cores=int(cpu_cores),
            gpu_ids=gpu_ids,
        )
        with (job.run_dir / "stdout.log").open("a", encoding="utf-8") as stdout, (job.run_dir / "stderr.log").open(
            "a", encoding="utf-8"
        ) as stderr:
            stdout.write("$ " + " ".join(command) + "\n")
            stdout.flush()
            completed = subprocess.run(command, stdout=stdout, stderr=stderr, check=False)

        structure_paths = sorted(
            path
            for path in output_dir.rglob("*.cif")
            if path.is_file() and ".rowpair" not in path.parts
        )
        if completed.returncode != 0:
            raise RuntimeError(
                f"{model_profile['label']} exited with code {completed.returncode}; produced {len(structure_paths)} structure(s)."
            )
        if not structure_paths:
            raise RuntimeError(f"{model_profile['label']} completed without producing a CIF structure.")

        sample_structure_paths = [
            path for path in structure_paths if _is_af3_sample_structure(path)
        ]
        # AF3 also writes a top-level copy of its selected model alongside the
        # per-seed sample files. Use the actual sample files for prediction rows
        # when present so that copied aggregate output is not shown as a sixth
        # prediction. Keep every file in artifact_paths below for download.
        prediction_paths = sample_structure_paths or structure_paths
        candidates: list[dict] = []
        for index, structure_path in enumerate(prediction_paths, start=1):
            summary = _summary_confidence(structure_path)
            confidence, confidence_path = _confidence_data(structure_path)
            confidence_plddt = confidence.get("atom_plddts") or confidence.get("plddt")
            if isinstance(confidence_plddt, list):
                try:
                    confidence_mean = sum(float(value) for value in confidence_plddt) / len(confidence_plddt)
                except (TypeError, ValueError, ZeroDivisionError):
                    confidence_mean = None
            else:
                confidence_mean = None
            plddt_mean = _mean_plddt_from_cif(structure_path)
            metrics = {
                "plddt_mean": plddt_mean if plddt_mean is not None else confidence_mean,
                "ptm": summary.get("ptm", confidence.get("ptm")),
                "iptm": summary.get("iptm", confidence.get("iptm")),
                "ranking_score": summary.get(
                    "ranking_score",
                    summary.get("confidence_score", confidence.get("ranking_score", confidence.get("confidence_score"))),
                ),
                "prediction_index": index,
                "seed": _seed_from_name(structure_path),
                "optimization_mode": optimization_mode,
                "msa_source": msa_source,
                "model_name": model_name,
            }
            for confidence_source in (summary, confidence):
                for metric_name, metric_value in confidence_source.items():
                    normalized_name = str(metric_name).casefold().replace("-", "_")
                    if not normalized_name.startswith(("affinity_", "binding_probability", "complex_", "ligand_", "protein_")):
                        continue
                    numeric_value = _numeric_metric(metric_value)
                    if numeric_value is not None:
                        metrics[normalized_name] = numeric_value
            candidates.append(
                {
                    "candidate_id": f"af3_jax_cofold_{index:04d}",
                    "stage": STAGE_COMPLEX_REFOLDING,
                    "source_tool": TOOL_NAME,
                    "complex_pdb": str(structure_path.relative_to(job.run_dir)),
                    "metrics": metrics,
                    "raw_metadata": {
                        "structure_file": str(structure_path.relative_to(job.run_dir)),
                        "protein_chains": {chain["id"]: chain["sequence"] for chain in proteins},
                        "dna_chains": {chain["id"]: chain["sequence"] for chain in dna},
                        "rna_chains": {chain["id"]: chain["sequence"] for chain in rna},
                        "ligand_ccd_codes": ccd_codes,
                        "ligand_smiles": smiles,
                        "model_name": model_name,
                        "model_label": model_profile["label"],
                        "optimization_mode": optimization_mode,
                        "container_image": AF3_IMAGE,
                        **(
                            {"pae_path": str(confidence_path.relative_to(job.run_dir))}
                            if confidence_path
                            and (confidence.get("pae") or confidence.get("predicted_aligned_error"))
                            else {}
                        ),
                    },
                }
            )
        normalized = write_candidates(job.run_dir, TOOL_NAME, candidates)
        artifact_paths = [
            str(input_file.relative_to(job.run_dir)),
            str((raw_dir / "input_mapping.json").relative_to(job.run_dir)),
            str((raw_dir / "input_sequences.fasta").relative_to(job.run_dir)),
            *[
                str(path.relative_to(job.run_dir))
                for path in sorted(msa_dir.glob("chain_*.a3m"))
                if path.is_file()
            ],
            *[str(path.relative_to(job.run_dir)) for path in structure_paths],
        ]
        for structure_path in structure_paths:
            stem = structure_path.stem[:-6] if structure_path.stem.endswith("_model") else structure_path.stem
            for suffix in ("_summary_confidences.json", "_confidences.json"):
                confidence_path = structure_path.with_name(f"{stem}{suffix}")
                if confidence_path.is_file():
                    artifact_paths.append(str(confidence_path.relative_to(job.run_dir)))
            affinity_stems = list(dict.fromkeys((stem, re.sub(r"_model_\d+$", "", stem))))
            for affinity_stem in affinity_stems:
                affinity_path = structure_path.with_name(f"affinity_{affinity_stem}.json")
                if affinity_path.is_file():
                    artifact_paths.append(str(affinity_path.relative_to(job.run_dir)))
        artifact_paths = list(dict.fromkeys(artifact_paths))
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
                    "protein_chain_count": len(proteins),
                    "dna_chain_count": len(dna),
                    "rna_chain_count": len(rna),
                    "ligand_count": len(ccd_codes) + len(smiles),
                    "model_name": model_name,
                    "model_label": model_profile["label"],
                    "msa_source": msa_source,
                    "msa_source_by_chain": msa_source_by_chain,
                    "optimization_mode": optimization_mode,
                    "container_image": AF3_IMAGE,
                    "seeds": seeds,
                    "num_diffusion_samples": int(num_diffusion_samples),
                    "num_recycles": int(num_recycles),
                    "accept_alphafold3_terms": bool(model_profile["terms"] and accept_alphafold3_terms),
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


def _seed_from_name(path: Path) -> int | None:
    match = re.search(r"seed[_-]?(\d+)", path.stem, flags=re.IGNORECASE)
    return int(match.group(1)) if match else None


def _is_af3_sample_structure(path: Path) -> bool:
    return bool(re.search(r"seed[-_]?\d+_sample[-_]?\d+", path.as_posix(), flags=re.IGNORECASE))
