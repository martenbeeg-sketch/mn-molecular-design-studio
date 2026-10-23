from __future__ import annotations

import glob
import csv
import gzip
import hashlib
import json
import random
import re
import shlex
import shutil
import subprocess
from pathlib import Path

from mn_molecular_design_studio.engine.workflow_engine.core.artifacts import Artifact, artifact_path
from mn_molecular_design_studio.engine.workflow_engine.core.candidates import (
    STAGE_COMPLEX_REFOLDING,
    STAGE_GENERATION_BACKBONE,
    STAGE_GENERATION_BACKBONE_SEQUENCE,
    candidate_stage_counts,
    read_candidates,
    write_candidates,
)
from mn_molecular_design_studio.engine.workflow_engine.core.gpu import docker_gpu_args, normalize_gpu_device
from mn_molecular_design_studio.engine.workflow_engine.core.hotspot_metrics import (
    calculate_hotspot_metrics,
    passes_hotspot_prefilter,
)
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import collect_jobs, create_job, finish_job, read_json, update_status, write_json
from mn_molecular_design_studio.engine.workflow_engine.core.manifests import load_manifest
from mn_molecular_design_studio.engine.workflow_engine.core.scheduler import apply_docker_cpu_limit, apply_docker_cpu_limits_to_steps
from mn_molecular_design_studio.engine.workflow_engine.core.structures import (
    filter_pdb_text,
    mmcif_protein_residues_by_chain,
    pdb_summary,
)
from mn_molecular_design_studio.engine.workflow_engine.runtime import reference_root
from mn_molecular_design_studio.engine.workflow_engine.workflows.target_msa import (
    ensure_boltz_msas_for_target as _shared_ensure_boltz_msas_for_target,
    ensure_pxdesign_msa_dirs_for_target,
)
from mn_molecular_design_studio.mimosa_input import normalize_mimosa_source_positions


DESIGN_GROUP = "design"
MIMOSA_BOLTZGEN_IMAGE = "mn-boltzgen-opt:mimosa-local"
RFDIFFUSION_SCAFFOLD_LIBRARY_CONTAINER_DIR = "/models/ppi_scaffolds"
BINDCRAFT_RESOURCE_DIR = Path(__file__).resolve().parents[1] / "data" / "bindcraft"
RFDIFFUSION_BUNDLED_SCAFFOLD_TAR = (
    Path(__file__).resolve().parents[1] / "data" / "reference" / "ppi_scaffolds_subset.tar.gz"
)
AA3_TO_1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "MSE": "M",
    "SEC": "U",
    "PYL": "O",
    "ASX": "B",
    "GLX": "Z",
}


def rfdiffusion_scaffold_library_status(scaffold_dir: Path | None = None) -> dict:
    scaffold_dir = scaffold_dir or reference_root() / "rfdiffusion_models" / "ppi_scaffolds"
    ss_files = sorted(scaffold_dir.glob("*_ss.pt")) if scaffold_dir.exists() else []
    adj_files = sorted(scaffold_dir.glob("*_adj.pt")) if scaffold_dir.exists() else []
    return {
        "path": str(scaffold_dir),
        "exists": scaffold_dir.exists(),
        "ss_count": len(ss_files),
        "adj_count": len(adj_files),
        "ready": bool(ss_files and adj_files),
        "bundled_tar": str(RFDIFFUSION_BUNDLED_SCAFFOLD_TAR),
        "bundled_tar_exists": RFDIFFUSION_BUNDLED_SCAFFOLD_TAR.exists(),
    }


def prepare_rfdiffusion_scaffold_library(scaffold_dir: Path | None = None) -> dict:
    import tarfile

    scaffold_dir = scaffold_dir or reference_root() / "rfdiffusion_models" / "ppi_scaffolds"
    status = rfdiffusion_scaffold_library_status(scaffold_dir)
    if status["ready"]:
        return status
    if not RFDIFFUSION_BUNDLED_SCAFFOLD_TAR.exists():
        raise FileNotFoundError(
            f"Bundled RFdiffusion scaffold archive not found: {RFDIFFUSION_BUNDLED_SCAFFOLD_TAR}"
        )
    scaffold_dir.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(RFDIFFUSION_BUNDLED_SCAFFOLD_TAR, "r:gz") as archive:
        parent = scaffold_dir.parent.resolve()
        for member in archive.getmembers():
            target = (scaffold_dir.parent / member.name).resolve()
            if parent != target and parent not in target.parents:
                raise ValueError(f"Unsafe path in RFdiffusion scaffold archive: {member.name}")
        archive.extractall(scaffold_dir.parent)
    return rfdiffusion_scaffold_library_status(scaffold_dir)


def rfdiffusion_scaffold_ids(scaffold_dir: Path | None = None) -> list[str]:
    scaffold_dir = scaffold_dir or reference_root() / "rfdiffusion_models" / "ppi_scaffolds"
    if not scaffold_dir.exists():
        return []
    return sorted(path.name[: -len("_ss.pt")] for path in scaffold_dir.glob("*_ss.pt"))


def random_rfdiffusion_scaffold_id(scaffold_dir: Path | None = None) -> str:
    scaffold_ids = rfdiffusion_scaffold_ids(scaffold_dir)
    if not scaffold_ids:
        raise RuntimeError("RFdiffusion scaffold-guided mode has no scaffold *_ss.pt files to choose from.")
    return random.choice(scaffold_ids)


def _has_rfdiffusion_scaffold_list_override(run_parameters: str) -> bool:
    return bool(re.search(r"(?:^|\s)\+{0,2}scaffoldguided\.scaffold_list=", run_parameters))


def prepared_design_targets() -> list[dict]:
    """Return every target structure that can be staged for binder design."""
    from mn_molecular_design_studio.engine.workflow_engine.workflows.detection import ppi_target_jobs

    return ppi_target_jobs()


def target_label(row: dict) -> str:
    chains = ",".join(row.get("chains") or [])
    source = str(row.get("source_category") or row.get("prepared_kind") or "target")
    provenance = str(row.get("source_label") or row.get("tool") or "").strip()
    provenance_text = f", source {provenance}" if provenance else ""
    return f"{row['target_name']} [{source}] ({row['job_code']}, chains {chains or 'unknown'}{provenance_text})"


def design_jobs_for_target(target_pdb: Path) -> list[dict]:
    target_path = str(target_pdb)
    rows: list[dict] = []
    for row in collect_jobs(DESIGN_GROUP):
        run_dir = Path(row["run_dir"])
        payload = read_json(run_dir / "input.json")
        inputs = payload.get("inputs") or {}
        if inputs.get("target_pdb") != target_path:
            continue
        result = read_json(run_dir / "result.json")
        rows.append(
            {
                **row,
                "success": result.get("success"),
                "candidate_count": (result.get("metrics") or {}).get("candidate_count", ""),
            }
        )
    return rows


def default_target_contig(target_pdb: Path, target_chains: list[str], binder_length: str) -> str:
    summary = pdb_summary(target_pdb.read_text(errors="ignore"))
    chain_rows = {row["chain_id"]: row for row in summary.get("chains", [])}
    segments: list[str] = []
    for chain_id in target_chains:
        row = chain_rows.get(chain_id)
        residues = sorted(row.get("residues") or []) if row else []
        if not residues:
            continue
        start = end = residues[0]
        for residue in residues[1:]:
            if residue == end + 1:
                end = residue
                continue
            segments.append(f"{chain_id}{start}-{end}")
            start = end = residue
        segments.append(f"{chain_id}{start}-{end}")
    return "/".join(segments) + f"/0 {binder_length}" if segments else ""


def _target_input_spec(target_pdb: Path, target_chains: list[str]) -> str:
    summary = pdb_summary(target_pdb.read_text(errors="ignore"))
    chain_rows = {row["chain_id"]: row for row in summary.get("chains", [])}
    segments: list[str] = []
    for chain_id in target_chains:
        row = chain_rows.get(chain_id)
        residues = sorted(row.get("residues") or []) if row else []
        if not residues:
            continue
        start = end = residues[0]
        for residue in residues[1:]:
            if residue == end + 1:
                end = residue
                continue
            segments.append(f"{chain_id}{start}-{end}")
            start = end = residue
        segments.append(f"{chain_id}{start}-{end}")
    return ",".join(segments)


def normalize_hotspots(text: str) -> str:
    tokens = []
    for raw in re.split(r"[,;\s]+", text.strip()):
        token = raw.strip()
        if not token:
            continue
        token = token.replace(":", "")
        if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_.-]*[0-9]+", token):
            raise ValueError(f"Hotspot '{raw}' must contain a chain ID followed by a residue number, such as A123.")
        tokens.append(token)
    return ",".join(tokens)


def validate_rfdiffusion_inputs(contig: str, hotspots: str, binder_length: str, num_designs: int, timesteps: int) -> None:
    if not contig.strip():
        raise ValueError("Contig is required.")
    if "/0 " not in contig:
        raise ValueError("Binder contig should look like A10-150/0 55-55.")
    if not re.fullmatch(r"\d+(-\d+)?", binder_length.strip()):
        raise ValueError("Binder length should look like 55 or 55-80.")
    if num_designs < 1:
        raise ValueError("Number of designs must be at least 1.")
    if timesteps < 1:
        raise ValueError("Timesteps must be at least 1.")
    if hotspots:
        normalize_hotspots(hotspots)


def validate_rfdiffusion_optional_params(contigmap_length: str = "", inpaint_seq: str = "") -> None:
    if contigmap_length and not re.fullmatch(r"\d+(-\d+)?", contigmap_length.strip()):
        raise ValueError("contigmap.length should look like 123 or 123-456.")
    if inpaint_seq:
        for segment in inpaint_seq.split("/"):
            if not re.fullmatch(r"[A-Za-z][0-9]+(-[0-9]+)?", segment.strip()):
                raise ValueError("contigmap.inpaint_seq should look like A10-20/A22/B30-40.")


def build_rfdiffusion_scaffoldguided_parameters(
    enabled: bool,
    scaffold_dir: str = RFDIFFUSION_SCAFFOLD_LIBRARY_CONTAINER_DIR,
    target_pdb: bool = True,
    target_path: str = "",
    target_ss: str = "",
    target_adj: str = "",
) -> str:
    if not enabled:
        return ""
    scaffold_dir = str(scaffold_dir or RFDIFFUSION_SCAFFOLD_LIBRARY_CONTAINER_DIR).strip().rstrip("/")
    if not scaffold_dir:
        raise ValueError("RFdiffusion scaffold-guided mode requires a scaffold directory.")
    parts = [
        "++scaffoldguided.scaffoldguided=True",
        f"++scaffoldguided.scaffold_dir={scaffold_dir}/",
    ]
    if target_pdb:
        parts.append("++scaffoldguided.target_pdb=True")
    if target_path.strip():
        parts.append(f"++scaffoldguided.target_path={target_path.strip()}")
    if target_ss.strip():
        parts.append(f"++scaffoldguided.target_ss={target_ss.strip()}")
    if target_adj.strip():
        parts.append(f"++scaffoldguided.target_adj={target_adj.strip()}")
    return " ".join(parts)


def build_rfdiffusion_run_parameters(
    timesteps: int,
    partial_diffusion: bool = False,
    contigmap_length: str = "",
    inpaint_seq: str = "",
    model_weights: str = "Complex_base",
    deterministic: bool = False,
    noise_scale_ca: str = "",
    noise_scale_frame: str = "",
    scaffoldguided: bool = False,
    scaffold_dir: str = RFDIFFUSION_SCAFFOLD_LIBRARY_CONTAINER_DIR,
    scaffold_target_pdb: bool = True,
    scaffold_target_path: str = "",
    scaffold_target_ss: str = "",
    scaffold_target_adj: str = "",
    extra_run_parameters: str = "",
) -> str:
    validate_rfdiffusion_optional_params(contigmap_length, inpaint_seq)
    parts = [f"diffuser.partial_T={timesteps}" if partial_diffusion else f"diffuser.T={timesteps}"]
    if contigmap_length.strip():
        length = contigmap_length.strip()
        if "-" not in length:
            length = f"{length}-{length}"
        parts.append(f"contigmap.length={length}")
    if inpaint_seq.strip():
        parts.append(f"contigmap.inpaint_seq=[{inpaint_seq.strip()}]")
    if model_weights not in {"", "Base", "Complex_base"}:
        parts.append(f"inference.ckpt_override_path=/models/{model_weights}_ckpt.pt")
    if deterministic:
        parts.append("inference.deterministic=True")
    if noise_scale_ca.strip():
        parts.append(f"denoiser.noise_scale_ca={noise_scale_ca.strip()}")
    if noise_scale_frame.strip():
        parts.append(f"denoiser.noise_scale_frame={noise_scale_frame.strip()}")
    scaffoldguided_parameters = build_rfdiffusion_scaffoldguided_parameters(
        scaffoldguided,
        scaffold_dir=scaffold_dir,
        target_pdb=scaffold_target_pdb,
        target_path=scaffold_target_path,
        target_ss=scaffold_target_ss,
        target_adj=scaffold_target_adj,
    )
    if scaffoldguided_parameters:
        parts.append(scaffoldguided_parameters)
    if extra_run_parameters.strip():
        parts.append(extra_run_parameters.strip())
    return " ".join(parts)


def parse_binder_lengths(text: str) -> list[int]:
    value = str(text or "").strip()
    if not value:
        raise ValueError("Binder length is required.")
    match = re.fullmatch(r"(\d+)(?:\s*(?:-|,|\s)\s*(\d+))?", value)
    if not match:
        raise ValueError("Binder length should be a single value or a two-value range.")
    lengths = [int(value) for value in match.groups() if value is not None]
    if any(length < 1 for length in lengths):
        raise ValueError("Binder lengths must be positive.")
    if len(lengths) > 2:
        raise ValueError("Binder length should be a single value or a min-max range.")
    if len(lengths) == 2 and lengths[1] < lengths[0]:
        raise ValueError("Binder length max must be greater than or equal to min.")
    return [lengths[0]] if len(lengths) == 2 and lengths[0] == lengths[1] else lengths


def _bindcraft_resource(*parts: str) -> Path:
    return BINDCRAFT_RESOURCE_DIR.joinpath(*parts)


def _read_json_file(path: Path) -> dict:
    if not path.exists():
        raise FileNotFoundError(f"BindCraft settings file not found: {path}")
    return json.loads(path.read_text())


def _bindcraft_advanced_settings(
    settings_file: str,
    max_trajectories: int,
    num_seqs_override: int | None = None,
    max_mpnn_sequences_override: int | None = None,
    enable_mpnn: bool = True,
) -> dict:
    settings = _read_json_file(_bindcraft_resource("settings_advanced", settings_file))
    settings.update(
        {
            "max_trajectories": max_trajectories,
            "save_design_animations": False,
            "save_design_trajectory_plots": False,
            "save_trajectory_pickle": False,
            "remove_unrelaxed_trajectory": False,
            "remove_unrelaxed_complex": False,
            "remove_binder_monomer": False,
            "enable_mpnn": bool(enable_mpnn),
            "af_params_dir": "alphafold_models_path",
        }
    )
    if num_seqs_override is not None:
        settings["num_seqs"] = num_seqs_override
    if max_mpnn_sequences_override is not None:
        settings["max_mpnn_sequences"] = max_mpnn_sequences_override
    return settings


def _copy_bindcraft_settings(run_dir: Path, params: dict) -> None:
    bindcraft_dir = run_dir / "artifacts" / "raw" / "bindcraft"
    bindcraft_dir.mkdir(parents=True, exist_ok=True)
    binder_name = params.get("campaign_name") or "design"
    binder_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", binder_name.strip()).strip("_") or "design"
    input_dict = {
        "design_path": "output",
        "starting_pdb": "target.pdb",
        "binder_name": binder_name,
        "chains": ",".join(params["target_chains"]),
        "target_hotspot_residues": params.get("hotspots") or "",
        "lengths": parse_binder_lengths(params["binder_length"]),
        "number_of_final_designs": params["number_of_final_designs"],
    }
    write_json(bindcraft_dir / "input.json", input_dict)
    write_json(
        bindcraft_dir / "settings_advanced.json",
        _bindcraft_advanced_settings(
            params["advanced_settings_file"],
            params["max_trajectories"],
            params.get("num_seqs_override"),
            params.get("max_mpnn_sequences_override"),
            bool(params.get("enable_mpnn", True)),
        ),
    )
    filters_path = _bindcraft_resource("settings_filters", params["filter_settings"])
    write_json(bindcraft_dir / "settings_filters.json", _read_json_file(filters_path))


def _read_bindcraft_stats(output_dir: Path) -> dict[str, dict]:
    stats: dict[str, dict] = {}
    for csv_name in ["final_design_stats.csv", "trajectory_stats.csv", "mpnn_design_stats.csv"]:
        path = output_dir / csv_name
        if not path.exists():
            continue
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle):
                design_id = row.get("Design") or row.get("Name")
                if not design_id:
                    continue
                cleaned = {k: v for k, v in row.items() if v not in {"", None}}
                for key, value in list(cleaned.items()):
                    if key.startswith("Average_"):
                        cleaned.setdefault(key.removeprefix("Average_"), value)
                stats.setdefault(design_id, {}).update(cleaned)
    return stats


def _bindcraft_stats_for_design(stats: dict[str, dict], design_name: str) -> dict:
    if design_name in stats:
        return stats[design_name]
    model_stripped = re.sub(r"_model\d+$", "", design_name)
    if model_stripped in stats:
        return stats[model_stripped]
    for stats_name, row in stats.items():
        if design_name.startswith(stats_name) or stats_name.startswith(model_stripped):
            return row
    return {}


def _normalize_bindcraft_candidates(run_dir: Path, params: dict, target_artifact: Path) -> list[dict]:
    output_dir = run_dir / "artifacts" / "raw" / "bindcraft" / "output"
    stats = _read_bindcraft_stats(output_dir)
    pdb_by_design: dict[str, Path] = {}
    patterns = ["Accepted/*.pdb"]
    if not bool(params.get("enable_mpnn", True)):
        patterns.extend(["Trajectory/Relaxed/*.pdb", "Trajectory/*.pdb"])
    for pattern in patterns:
        for pdb_path in sorted(output_dir.glob(pattern)):
            pdb_by_design.setdefault(pdb_path.stem, pdb_path)
    candidates: list[dict] = []
    for index, pdb_path in enumerate(pdb_by_design.values(), start=1):
        design_name = pdb_path.stem
        metrics = _bindcraft_stats_for_design(stats, design_name)
        trajectory_name = re.sub(r"_mpnn\d+(?:_model\d+)?$", "", design_name)
        trajectory_path = output_dir / "Trajectory" / f"{trajectory_name}.pdb"
        output_kind = pdb_path.parent.name
        is_generator_trajectory = output_kind in {"Trajectory", "Relaxed"}
        result_kind = "generator_trajectory" if is_generator_trajectory else "native_pipeline"
        binder_chains, target_chains, chain_inference = _infer_generated_chain_roles(
            target_artifact,
            pdb_path,
            params.get("target_chains", []),
        )
        binder_sequence = metrics.get("Sequence")
        if not binder_sequence and binder_chains:
            sequences = _pdb_sequences_by_chain(pdb_path)
            binder_sequence = "".join(sequences.get(chain, "") for chain in binder_chains) or None
        try:
            complex_path = str(pdb_path.relative_to(run_dir))
            target_path = str(target_artifact.relative_to(run_dir))
        except ValueError:
            complex_path = str(pdb_path)
            target_path = str(target_artifact)
        candidates.append(
            {
                "candidate_id": f"bindcraft_{index:05d}",
                "source_tool": "bindcraft",
                "stage": STAGE_GENERATION_BACKBONE_SEQUENCE
                if is_generator_trajectory
                else STAGE_COMPLEX_REFOLDING,
                "target_pdb": target_path,
                "complex_pdb": complex_path,
                "binder_pdb": None,
                "binder_sequence": binder_sequence,
                "target_chains": target_chains,
                "binder_chains": binder_chains,
                "hotspots": [token for token in params.get("hotspots", "").split(",") if token],
                "binder_length": params.get("binder_length"),
                "contig": None,
                "metrics": {
                    **metrics,
                    "complex_refolding_backend": "bindcraft" if not is_generator_trajectory else "",
                    "result_kind": result_kind,
                    "native_pass_filters": output_kind == "Accepted",
                    "target_chain_inference": chain_inference,
                },
                "raw_metadata": {
                    "result_kind": result_kind,
                    "bindcraft_design": design_name,
                    "input_target_chains": params.get("target_chains", []),
                    "target_chain_inference": chain_inference,
                    "bindcraft_output_kind": output_kind,
                    "design_reference_pdb": _rel_path(run_dir, trajectory_path)
                    if trajectory_path.exists()
                    else None,
                    "design_reference_kind": "bindcraft_trajectory"
                    if trajectory_path.exists()
                    else "",
                },
            }
        )
    return write_candidates(run_dir, "bindcraft", candidates)


def _rel_path(run_dir: Path, path: Path | None) -> str | None:
    if path is None:
        return None
    try:
        return str(path.relative_to(run_dir))
    except ValueError:
        return str(path)


def _write_candidates(run_dir: Path, tool: str, candidates: list[dict]) -> None:
    write_candidates(run_dir, tool, candidates)


def _candidate_structure_key(candidate: dict) -> str:
    return str(candidate.get("complex_pdb") or candidate.get("binder_pdb") or candidate.get("candidate_id") or "")


def _with_pool_metadata(candidate: dict, pool_level: str, coverage: str) -> dict:
    raw_metadata = dict(candidate.get("raw_metadata") or {})
    raw_metadata.setdefault("candidate_pool_level", pool_level)
    raw_metadata.setdefault("candidate_pool_coverage", coverage)
    metrics = dict(candidate.get("metrics") or {})
    metrics.setdefault("candidate_pool_level", pool_level)
    candidate = dict(candidate)
    candidate["raw_metadata"] = raw_metadata
    candidate["metrics"] = metrics
    return candidate


def _merge_candidate_pools(
    *,
    tool: str,
    native_candidates: list[dict],
    raw_candidates: list[dict],
    native_coverage: str,
    raw_coverage: str,
) -> list[dict]:
    merged: list[dict] = []
    seen: set[str] = set()
    for candidate in native_candidates:
        enriched = _with_pool_metadata(candidate, "native_pipeline", native_coverage)
        key = _candidate_structure_key(enriched)
        seen.add(key)
        merged.append(enriched)
    for candidate in raw_candidates:
        key = _candidate_structure_key(candidate)
        if key in seen:
            continue
        seen.add(key)
        index = len(merged) + 1
        enriched = _with_pool_metadata(candidate, "prefilter_generated", raw_coverage)
        enriched["candidate_id"] = f"{tool}_prefilter_{index:05d}"
        merged.append(enriched)
    return merged


def _generic_candidate_records(
    run_dir: Path,
    tool: str,
    stage: str,
    params: dict,
    target_artifact: Path | None,
    structure_patterns: list[str],
    metadata_patterns: list[str] | None = None,
) -> list[dict]:
    metadata_patterns = metadata_patterns or []
    structure_paths: list[Path] = []
    for pattern in structure_patterns:
        structure_paths.extend(sorted((run_dir / "artifacts").glob(pattern)))
    metadata_paths: list[Path] = []
    for pattern in metadata_patterns:
        metadata_paths.extend(sorted((run_dir / "artifacts").glob(pattern)))
    candidates: list[dict] = []
    for index, structure_path in enumerate(structure_paths, start=1):
        metadata_path = metadata_paths[index - 1] if index - 1 < len(metadata_paths) else None
        candidates.append(
            {
                "candidate_id": f"{tool}_{index:05d}",
                "source_tool": tool,
                "stage": stage,
                "target_pdb": _rel_path(run_dir, target_artifact),
                "complex_pdb": _rel_path(run_dir, structure_path),
                "binder_pdb": None,
                "binder_sequence": None,
                "target_chains": params.get("target_chains", []),
                "hotspots": [token for token in params.get("hotspots", "").split(",") if token],
                "binder_length": params.get("binder_length"),
                "contig": params.get("contig"),
                "metrics": {},
                "raw_metadata": {"metadata_path": _rel_path(run_dir, metadata_path)},
            }
        )
    _write_candidates(run_dir, tool, candidates)
    return candidates


def _read_csv_rows(path: Path | None) -> list[dict[str, str]]:
    if path is None or not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _normalize_boltzgen_generation_candidates(
    run_dir: Path,
    params: dict,
    target_artifact: Path | None,
    *,
    write_output: bool = True,
) -> list[dict]:
    source_tool = "boltzgen_generator" if bool(params.get("generator_only_tile")) else "boltzgen"
    structure_paths = sorted(
        [
            *(run_dir / "artifacts").glob("raw/boltzgen/run-generation-only/intermediate_designs/*.cif"),
            *(run_dir / "artifacts").glob("raw/boltzgen/run-vanilla/intermediate_designs/*.cif"),
            *(run_dir / "artifacts").glob("raw/boltzgen/run-vanilla/rejected_backbones/*.cif"),
        ]
    )
    metadata_files = [
        *(run_dir / "artifacts").glob("raw/boltzgen/**/intermediate_designs/*.npz"),
        *(run_dir / "artifacts").glob("raw/boltzgen/**/rejected_backbones/*.npz"),
    ]
    metadata_paths = {path.stem: path for path in sorted(metadata_files)}
    integrity_report_path = run_dir / "artifacts/raw/boltzgen/run-vanilla/mimosa_integrity_results.json"
    integrity_by_file = {}
    if integrity_report_path.is_file():
        try:
            integrity_by_file = {
                row.get("file"): row for row in json.loads(integrity_report_path.read_text()).get("candidates", [])
            }
        except (OSError, ValueError, TypeError):
            integrity_by_file = {}
    candidates: list[dict] = []
    for index, structure_path in enumerate(structure_paths, start=1):
        binder_chains, target_chains, chain_inference = _infer_generated_chain_roles(
            target_artifact,
            structure_path,
            params.get("target_chains", []),
        )
        sequences = _sequences_by_chain(structure_path)
        binder_sequence = "".join(sequences.get(chain, "") for chain in binder_chains) or None
        metadata_path = metadata_paths.get(structure_path.stem)
        mimosa_metrics, mimosa_diagnostics = _mimosa_generation_metadata(metadata_path)
        backbone_integrity = _mimosa_backbone_integrity(structure_path, binder_chains) if params.get("mimosa_enabled") else {}
        integrity_record = integrity_by_file.get(structure_path.name, {})
        if integrity_record:
            backbone_integrity = {
                "backbone_integrity_status": integrity_record.get("status", "not_checked"),
                "backbone_bad_bond_count": int(integrity_record.get("bad_bond_count", 0)),
                "backbone_checked_bond_count": int(integrity_record.get("checked_bond_count", 0)),
            }
        candidate = {
            "candidate_id": f"boltzgen_{index:05d}",
            "source_tool": source_tool,
            "stage": STAGE_GENERATION_BACKBONE_SEQUENCE,
            "target_pdb": _rel_path(run_dir, target_artifact) if target_artifact else None,
            "complex_pdb": _rel_path(run_dir, structure_path),
            "binder_pdb": None,
            "binder_sequence": binder_sequence,
            "target_chains": target_chains,
            "binder_chains": binder_chains,
            "hotspots": [token for token in params.get("hotspots", "").split(",") if token],
            "binder_length": params.get("binder_length"),
            "contig": None,
            "metrics": {
                "target_chain_inference": chain_inference,
                **mimosa_metrics,
                **{key: value for key, value in backbone_integrity.items() if key != "backbone_bad_bonds"},
                **({"result_kind": "structural_rejection"} if backbone_integrity.get("backbone_integrity_status") == "rejected" else {}),
                **({"result_kind": "generated"} if params.get("generator_only_tile") and backbone_integrity.get("backbone_integrity_status") != "rejected" else {}),
            },
            "raw_metadata": {
                "result_kind": "generation_only",
                "metadata_path": _rel_path(run_dir, metadata_path),
                "input_target_chains": params.get("target_chains", []),
                "target_chain_inference": chain_inference,
                **({"mimosa_diagnostics": mimosa_diagnostics} if mimosa_diagnostics else {}),
                **({"backbone_bad_bonds": backbone_integrity.get("backbone_bad_bonds", [])} if backbone_integrity else {}),
                **({"backbone_bad_bonds": integrity_record.get("bad_bonds", [])} if integrity_record else {}),
            },
        }
        candidates.append(candidate)
    return write_candidates(run_dir, source_tool, candidates) if write_output else candidates


def _boltzgen_metric_rows(run_dir: Path, budget: int) -> list[dict[str, str]]:
    ranked_dir = run_dir / "artifacts" / "raw" / "boltzgen" / "run-vanilla" / "final_ranked_designs"
    all_metrics = ranked_dir / "all_designs_metrics.csv"
    if all_metrics.exists():
        return _read_csv_rows(all_metrics)
    preferred = ranked_dir / f"final_designs_metrics_{budget}.csv"
    if preferred.exists():
        return _read_csv_rows(preferred)
    candidates = sorted(ranked_dir.glob("final_designs_metrics_*.csv"))
    return _read_csv_rows(candidates[0] if candidates else None)


def _boltzgen_final_design_path(run_dir: Path, budget: int, row: dict[str, str]) -> Path | None:
    ranked_dir = run_dir / "artifacts" / "raw" / "boltzgen" / "run-vanilla" / "final_ranked_designs"
    final_dirs = [ranked_dir / f"final_{budget}_designs", *sorted(ranked_dir.glob("final_*_designs"))]
    design_id = str(row.get("id") or "").strip()
    file_name = str(row.get("file_name") or "").strip()
    stems = [Path(file_name).stem] if file_name else []
    if design_id:
        stems.append(design_id)
    for final_dir in final_dirs:
        if not final_dir.exists():
            continue
        for stem in stems:
            for path in sorted(final_dir.glob(f"*{stem}*.cif")):
                if path.is_file():
                    return path
        if file_name:
            direct = final_dir / file_name
            if direct.exists():
                return direct
    return None


def _boltzgen_design_pae_index(run_dir: Path) -> dict[str, list[Path]]:
    """Index full BoltzGen PAE archives from the native refolding pipeline."""
    output_dir = run_dir / "artifacts" / "raw" / "boltzgen" / "run-vanilla"
    pae_files: set[Path] = set()
    if not output_dir.is_dir():
        return {}
    for path in output_dir.rglob("*.npz"):
        parts = {part.lower() for part in path.parts}
        if "fold_out_npz" in parts or "pae" in parts or "pae" in path.stem.lower():
            pae_files.add(path)

    index: dict[str, list[Path]] = {}
    for path in sorted(pae_files):
        stem = path.stem
        if stem.lower().startswith("pae_"):
            stem = stem[4:]
        elif stem.lower().startswith("pae-"):
            stem = stem[4:]
        if stem.lower().endswith("_full_pae"):
            stem = stem[:-9]
        elif stem.lower().endswith("_pae"):
            stem = stem[:-4]
        key = stem.casefold()
        if key:
            index.setdefault(key, []).append(path)
    return index


def _boltzgen_design_pae_path(
    index: dict[str, list[Path]], row: dict[str, str], structure_path: Path
) -> Path | None:
    aliases = [
        Path(str(row.get("file_name") or "")).stem,
        Path(str(row.get("id") or "")).stem,
        structure_path.stem,
        re.sub(r"^rank\d+_", "", structure_path.stem, flags=re.IGNORECASE),
    ]
    for alias in aliases:
        key = alias.casefold()
        matches = sorted(
            index.get(key, []),
            key=lambda path: (
                0 if structure_path.parent in path.parents else 1,
                0 if "final_ranked_designs" in path.parts else 1,
                len(path.parts),
            ),
        )
        for candidate in matches:
            try:
                import numpy as np

                with np.load(candidate, allow_pickle=False) as archive:
                    for array_key in archive.files:
                        if "pae" not in array_key.casefold():
                            continue
                        matrix = np.asarray(archive[array_key])
                        if matrix.ndim >= 2 and matrix.shape[-1] == matrix.shape[-2]:
                            return candidate
            except (OSError, ValueError, KeyError):
                continue
    return None


def _numeric_or_text(value: str) -> object:
    if value is None:
        return ""
    text = str(value).strip()
    if text == "":
        return ""
    lowered = text.lower()
    if lowered in {"true", "false"}:
        return lowered == "true"
    try:
        number = float(text)
    except ValueError:
        return text
    if number.is_integer():
        return int(number)
    return number


def _mimosa_generation_metadata(metadata_path: Path | None) -> tuple[dict, dict]:
    """Read the MIMOSA writer's scalar JSON fields without loading object arrays."""
    if metadata_path is None or not metadata_path.is_file():
        return {}, {}
    try:
        import numpy as np

        with np.load(metadata_path, allow_pickle=False) as archive:
            def scalar_json(key):
                if key not in archive.files:
                    return {}
                value = archive[key]
                try:
                    value = value.item()
                except ValueError:
                    pass
                return json.loads(str(value)) if value is not None else {}

            generated = scalar_json("mimosa_generated_metrics")
            diagnostics = scalar_json("mimosa_diagnostics")
            fixed = archive["mimosa_fixed_aa"] if "mimosa_fixed_aa" in archive.files else None
            fixed_count = int((fixed >= 0).sum()) if fixed is not None else 0
        metrics = {
            "generated_motif_rmsd_angstrom": generated.get("generated_motif_rmsd_angstrom"),
            "generated_motif_potential_half_sum_angstrom2": generated.get("potential_half_sum_angstrom_squared"),
            "fixed_motif_residue_count": fixed_count,
        }
        return {key: value for key, value in metrics.items() if value is not None}, diagnostics
    except (OSError, ValueError, KeyError, json.JSONDecodeError):
        return {}, {}


def _mimosa_backbone_integrity(structure_path: Path, binder_chains: list[str]) -> dict:
    """Check C–N bonds between consecutive modeled polymer residues in binder chains."""
    if not binder_chains:
        return {}
    try:
        model = _first_model(_parse_structure(structure_path))
    except Exception:
        return {}
    if model is None:
        return {}
    bad = []
    checked = 0
    source_gaps = 0
    for chain_id in binder_chains:
        if chain_id not in model:
            continue
        residues = [r for r in model[chain_id] if getattr(r, "id", ("",))[0] == " "]
        for left, right in zip(residues, residues[1:]):
            lid, rid = left.id, right.id
            # Chothia/IMGT residue numbering is intentionally discontinuous.
            # Check the actual ordered polymer neighbors, independent of IDs.
            if int(rid[1]) > int(lid[1]) + 1:
                source_gaps += int(rid[1]) - int(lid[1]) - 1
            if "C" not in left or "N" not in right:
                continue
            checked += 1
            distance = float(left["C"] - right["N"])
            if abs(distance - 1.329) > 0.5:
                bad.append({"chain": str(chain_id), "left": f"{lid[1]}{str(lid[2]).strip()}",
                            "right": f"{rid[1]}{str(rid[2]).strip()}", "distance_angstrom": distance})
    return {"backbone_integrity_status": "rejected" if bad else "passed" if checked else "not_checked",
            "backbone_bad_bond_count": len(bad), "backbone_checked_bond_count": checked,
            "backbone_unresolved_source_gap_count": source_gaps,
            "backbone_bad_bonds": bad}


def _float_or_none(value: object) -> float | None:
    if value in {None, ""}:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _normalize_boltzgen_vanilla_candidates(run_dir: Path, params: dict, target_artifact: Path | None) -> list[dict]:
    budget = int(params.get("budget") or 1)
    rows = _boltzgen_metric_rows(run_dir, budget)
    pae_index = _boltzgen_design_pae_index(run_dir)
    candidates: list[dict] = []
    for index, row in enumerate(rows, start=1):
        structure_path = _boltzgen_final_design_path(run_dir, budget, row)
        if structure_path is None:
            continue
        pae_path = _boltzgen_design_pae_path(pae_index, row, structure_path)
        binder_chains, target_chains, chain_inference = _infer_generated_chain_roles(
            target_artifact,
            structure_path,
            params.get("target_chains", []),
        )
        sequences = _sequences_by_chain(structure_path)
        # BoltzGen's ``designed_sequence`` column can contain only the
        # redesigned residues (for example, a CDR), while
        # ``designed_chain_sequence`` and the ranked CIF contain the complete
        # binder. Prefer the sequence parsed from the actual output structure
        # so the normalized candidate and database identity describe the
        # displayed complex.
        binder_sequence = (
            "".join(sequences.get(chain, "") for chain in binder_chains)
            or row.get("designed_chain_sequence")
            or row.get("sequence")
            or row.get("designed_sequence")
            or None
        )
        metrics = {key: _numeric_or_text(value) for key, value in row.items() if value not in {"", None}}
        metrics.update(
            {
                "complex_refolding_backend": "boltzgen_vanilla",
                "result_kind": "native_pipeline",
                "target_chain_inference": chain_inference,
                "boltzgen_has_full_pae": pae_path is not None,
            }
        )
        candidate = {
            "candidate_id": f"boltzgen_{index:05d}",
            "source_tool": "boltzgen",
            "stage": STAGE_COMPLEX_REFOLDING,
            "target_pdb": _rel_path(run_dir, target_artifact) if target_artifact else None,
            "complex_pdb": _rel_path(run_dir, structure_path),
            "binder_pdb": None,
            "binder_sequence": binder_sequence,
            "target_chains": target_chains,
            "binder_chains": binder_chains,
            "hotspots": [token for token in params.get("hotspots", "").split(",") if token],
            "binder_length": params.get("binder_length"),
            "contig": None,
            "metrics": metrics,
            "raw_metadata": {
                "result_kind": "native_pipeline",
                "boltzgen_design_id": row.get("id"),
                "boltzgen_file_name": row.get("file_name"),
                "input_target_chains": params.get("target_chains", []),
                "target_chain_inference": chain_inference,
                "output_kind": "final_ranked_design",
                "pae_source_path": _rel_path(run_dir, pae_path) if pae_path else None,
                "pae_path": _rel_path(run_dir, pae_path) if pae_path else None,
            },
        }
        candidates.append(candidate)
    raw_candidates = _normalize_boltzgen_generation_candidates(
        run_dir,
        params,
        target_artifact,
        write_output=False,
    )
    merged = _merge_candidate_pools(
        tool="boltzgen",
        native_candidates=candidates,
        raw_candidates=raw_candidates,
        native_coverage="final ranked native pipeline designs",
        raw_coverage="intermediate generated BoltzGen structures before final ranking",
    )
    return write_candidates(run_dir, "boltzgen", merged)


def _clean_key(text: str, fallback: str = "mn_app_target") -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(text or "").strip()).strip("_")
    return cleaned or fallback


def _pdb_residue_records_by_chain(path: Path) -> dict[str, list[tuple[int, str, str]]]:
    records: dict[str, list[tuple[int, str, str]]] = {}
    seen: set[tuple[str, str, str]] = set()
    for line in path.read_text(errors="ignore").splitlines():
        if not line.startswith("ATOM  "):
            continue
        chain = line[21].strip() or "_"
        residue_key = (chain, line[22:26].strip(), line[26].strip())
        if residue_key in seen:
            continue
        seen.add(residue_key)
        try:
            residue_number = int(line[22:26])
        except ValueError:
            continue
        records.setdefault(chain, []).append(
            (residue_number, line[26].strip(), AA3_TO_1.get(line[17:20].strip().upper(), "X"))
        )
    return {
        chain: sorted(chain_records, key=lambda item: (item[0], item[1]))
        for chain, chain_records in records.items()
    }


_A3M_SEQUENCE_ALLOWED = set("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz-")
_ALPHAFOLD_MODELS_DIR = Path("/mnt/db/reference_files/alphafold_models")
_BOLTZ_CACHE_DIR = Path("/mnt/db/reference_files/boltz_models")
_BOLTZ_MSA_REPOSITORY_DIR = _BOLTZ_CACHE_DIR / "msa_repository"


def _target_chain_sequences(path: Path, target_chains: list[str]) -> dict[str, str]:
    residue_records = _pdb_residue_records_by_chain(path)
    chains = [chain for chain in target_chains if chain in residue_records] or list(residue_records)
    return {chain: "".join(record[2] for record in residue_records.get(chain, [])) for chain in chains}


def _boltz_msa_paths(sequence: str, msa_repository_dir: Path = _BOLTZ_MSA_REPOSITORY_DIR) -> tuple[Path, str]:
    digest = hashlib.sha256(sequence.encode("utf-8")).hexdigest()
    filename = f"{digest}.a3m"
    return msa_repository_dir / filename, f"/msa_repository/{filename}"


def _validate_a3m_file(path: Path) -> tuple[bool, str]:
    if not path.exists():
        return False, "file does not exist"
    try:
        raw = path.read_bytes()
    except Exception as exc:
        return False, f"read failed: {exc}"
    if not raw:
        return False, "file is empty"
    if b"\x00" in raw:
        return False, "contains NUL byte(s)"
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as exc:
        return False, f"invalid UTF-8: {exc}"
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines:
        return False, "contains no non-empty lines"
    if not lines[0].startswith(">"):
        return False, "first non-empty line is not a FASTA header"
    saw_sequence = False
    for line in lines:
        if line.startswith(">"):
            continue
        saw_sequence = True
        if any(char not in _A3M_SEQUENCE_ALLOWED for char in line):
            return False, "contains invalid sequence characters"
    if not saw_sequence:
        return False, "contains headers only and no sequence"
    return True, ""


def _repair_a3m_file(path: Path) -> bool:
    if not path.exists():
        return False
    cleaned = path.read_bytes().rstrip(b"\x00").replace(b"\r\n", b"\n")
    valid, _reason = _validate_a3m_payload(cleaned)
    if not valid:
        return False
    if cleaned != path.read_bytes():
        path.write_bytes(cleaned)
    return True


def _validate_a3m_payload(raw: bytes) -> tuple[bool, str]:
    temp = Path("/tmp") / f"mn_molecular_design_studio.engine.workflow_engine_a3m_check_{hashlib.sha256(raw).hexdigest()}.a3m"
    try:
        temp.write_bytes(raw)
        return _validate_a3m_file(temp)
    finally:
        temp.unlink(missing_ok=True)


def _first_valid_a3m(root: Path) -> Path | None:
    for candidate in sorted(root.glob("**/*.a3m"), key=lambda path: ("processed" not in str(path), len(str(path)))):
        valid, _reason = _validate_a3m_file(candidate)
        if valid:
            return candidate
        if _repair_a3m_file(candidate):
            return candidate
    return None


def _write_boltz_msa_probe_yaml(path: Path, sequence: str) -> None:
    path.write_text(
        "\n".join(
            [
                "version: 1",
                "sequences:",
                "  - protein:",
                "      id: A",
                f"      sequence: {sequence}",
                "",
            ]
        )
    )


def _ensure_boltz_msa_for_sequence(run_dir: Path, sequence: str, label: str, gpu_device: object = "0") -> str:
    host_path, container_path = _boltz_msa_paths(sequence)
    host_path.parent.mkdir(parents=True, exist_ok=True)
    valid, reason = _validate_a3m_file(host_path)
    if not valid and host_path.exists():
        if _repair_a3m_file(host_path):
            valid, reason = _validate_a3m_file(host_path)
    if valid:
        return container_path

    probe_dir = run_dir / "artifacts" / "raw" / "genie3" / "boltz_msa_cache" / _clean_key(label)
    probe_dir.mkdir(parents=True, exist_ok=True)
    yaml_path = probe_dir / "input.yaml"
    _write_boltz_msa_probe_yaml(yaml_path, sequence)
    command = [
        "docker",
        "run",
        "--rm",
        *docker_gpu_args(gpu_device),
        "--shm-size=64G",
        "-v",
        f"{probe_dir}:/work",
        "-v",
        f"{_BOLTZ_CACHE_DIR}:/cache",
        "-v",
        f"{_BOLTZ_MSA_REPOSITORY_DIR}:/msa_repository",
        "-e",
        "BOLTZ_CACHE=/cache",
        "--ipc=host",
        "--shm-size=48G",
        "mn-boltz2:cu128",
        "predict",
        "/work/input.yaml",
        "--out_dir",
        "/work",
        "--sampling_steps",
        "200",
        "--recycling_steps",
        "3",
        "--diffusion_samples",
        "1",
        "--accelerator",
        "gpu",
        "--override",
        "--use_msa_server",
    ]
    with (run_dir / "stdout.log").open("a") as stdout, (run_dir / "stderr.log").open("a") as stderr:
        stdout.write(f"$ {' '.join(command)}\n")
        stdout.write(f"MSA cache miss for {label}: {host_path} ({reason})\n")
        stdout.flush()
        command = apply_docker_cpu_limit(command, run_dir)
        completed = subprocess.run(command, stdout=stdout, stderr=stderr, check=False)
    if completed.returncode != 0:
        raise RuntimeError(f"Boltz2 MSA preflight failed for {label} with return code {completed.returncode}.")
    generated = _first_valid_a3m(probe_dir)
    if generated is None:
        raise RuntimeError(f"Boltz2 MSA preflight did not produce a valid A3M for {label}.")
    payload = generated.read_bytes().rstrip(b"\x00").replace(b"\r\n", b"\n")
    valid_payload, payload_reason = _validate_a3m_payload(payload)
    if not valid_payload:
        raise RuntimeError(f"Boltz2 MSA preflight produced invalid A3M for {label}: {payload_reason}.")
    host_path.write_bytes(payload)
    valid, reason = _validate_a3m_file(host_path)
    if not valid:
        raise RuntimeError(f"Cached Boltz2 MSA is invalid for {label}: {reason}.")
    with (run_dir / "stdout.log").open("a") as stdout:
        stdout.write(f"MSA cached for {label}: {host_path}\n")
    return container_path


def _ensure_boltz_msas_for_target(run_dir: Path, target_artifact: Path, target_chains: list[str], gpu_device: object = "0") -> dict[str, str]:
    return _shared_ensure_boltz_msas_for_target(
        run_dir,
        target_artifact,
        target_chains,
        raw_subdir="genie3/boltz_msa_cache",
        gpu_device=gpu_device,
    )


def _write_genie3_target_files(
    *,
    run_dir: Path,
    target_artifact: Path,
    target_chains: list[str],
    binder_length: str,
    hotspots: str,
    campaign_name: str,
    boltz_msa_by_source_chain: dict[str, str] | None = None,
) -> tuple[Path, str, dict]:
    lengths = parse_binder_lengths(binder_length)
    min_length = lengths[0]
    max_length = lengths[-1]
    problem_key = _clean_key(campaign_name or target_artifact.stem)
    dataset_dir = run_dir / "artifacts" / "raw" / "genie3" / "dataset" / "mn_app"
    for subdir in ["problems", "targets/pdb", "targets/fasta", "targets/msa"]:
        (dataset_dir / subdir).mkdir(parents=True, exist_ok=True)

    residue_records = _pdb_residue_records_by_chain(target_artifact)
    source_chains = [chain for chain in target_chains if chain in residue_records] or list(residue_records)
    chain_map = {source_chain: chr(ord("B") + index) for index, source_chain in enumerate(source_chains)}
    source_by_new_chain = {new_chain: source_chain for source_chain, new_chain in chain_map.items()}
    residue_map: dict[str, str] = {}
    chain_tags: list[str] = []
    chain_sequences: dict[str, str] = {}
    for source_chain in source_chains:
        new_chain = chain_map[source_chain]
        records = residue_records.get(source_chain, [])
        sequence = "".join(record[2] for record in records)
        chain_sequences[new_chain] = sequence
        if records:
            chain_tags.append(f"{new_chain}1-{len(records)}")
        for new_index, (old_residue, _insertion, _aa) in enumerate(records, start=1):
            residue_map[f"{source_chain}{old_residue}"] = f"{new_chain}{new_index}"

    target_pdb = dataset_dir / "targets" / "pdb" / f"{problem_key}.pdb"
    target_fasta = dataset_dir / "targets" / "fasta" / f"{problem_key}.fasta"
    target_msa = dataset_dir / "targets" / "msa" / f"{problem_key}.a3m"
    target_pdb_by_chain: list[str] = []
    target_fasta_by_chain: list[str] = []
    target_msa_by_chain: list[str] = []

    rewritten_lines = [
        f"REMARK 999 KEY    {problem_key}\n",
        f"REMARK 999 NAME   {campaign_name or problem_key}\n",
    ]
    for tag in chain_tags:
        rewritten_lines.append(f"REMARK 999 TARGET {tag[0]} {tag[1:].split('-')[0].rjust(4)} {tag[1:].split('-')[1].rjust(4)}\n")
    source_text = target_artifact.read_text(errors="ignore")
    for line in source_text.splitlines():
        if not line.startswith("ATOM  "):
            continue
        source_chain = line[21].strip() or "_"
        if source_chain not in chain_map:
            continue
        try:
            old_residue = int(line[22:26])
        except ValueError:
            continue
        mapped = residue_map.get(f"{source_chain}{old_residue}")
        if not mapped:
            continue
        new_chain = mapped[0]
        new_residue = int(mapped[1:])
        # Genie3's problem preparation rejects altlocs/insertion codes, so the
        # app writes a normalized target PDB for the generated problem set.
        rewritten_lines.append(line[:16] + " " + line[17:21] + new_chain + str(new_residue).rjust(4) + " " + line[27:] + "\n")
    target_pdb.write_text("".join(rewritten_lines))

    merged_sequence = ":".join(chain_sequences[chain] for chain in sorted(chain_sequences))
    target_fasta.write_text(f">{problem_key}\n{merged_sequence}\n")
    target_msa.write_text(f">{problem_key}\n{merged_sequence}\n")
    for new_chain, sequence in sorted(chain_sequences.items()):
        chain_pdb = dataset_dir / "targets" / "pdb" / f"{problem_key}-chain_{new_chain}.pdb"
        chain_fasta = dataset_dir / "targets" / "fasta" / f"{problem_key}-chain_{new_chain}.fasta"
        chain_msa = dataset_dir / "targets" / "msa" / f"{problem_key}-chain_{new_chain}.a3m"
        chain_pdb.write_text("".join(line for line in rewritten_lines if line.startswith("ATOM  ") and line[21] == new_chain))
        chain_fasta.write_text(f">{problem_key}-chain_{new_chain}\n{sequence}\n")
        chain_msa.write_text(f">{problem_key}-chain_{new_chain}\n{sequence}\n")
        target_pdb_by_chain.append(_rel_path(run_dir, chain_pdb) or str(chain_pdb))
        target_fasta_by_chain.append(_rel_path(run_dir, chain_fasta) or str(chain_fasta))
        source_chain = source_by_new_chain.get(new_chain, "")
        repository_msa = (boltz_msa_by_source_chain or {}).get(source_chain)
        target_msa_by_chain.append(repository_msa or f"/work/{_rel_path(run_dir, chain_msa)}")

    mapped_hotspots = [residue_map[token] for token in hotspots.split(",") if token and token in residue_map]
    extended = sorted(
        {
            f"{token[0]}{neighbor}"
            for token in mapped_hotspots
            for neighbor in range(max(1, int(token[1:]) - 2), int(token[1:]) + 3)
            if any(tag.startswith(token[0]) and neighbor <= int(tag.split("-")[-1]) for tag in chain_tags)
        },
        key=lambda value: (value[0], int(value[1:])),
    )
    interface_residues = {
        "hotspot": mapped_hotspots,
        "extended": extended or mapped_hotspots,
        "common": mapped_hotspots,
    }
    problem = {
        "key": problem_key,
        "name": campaign_name or problem_key,
        "target_pdb_filepath": f"/work/{_rel_path(run_dir, target_pdb)}",
        "target_fasta_filepath": f"/work/{_rel_path(run_dir, target_fasta)}",
        "target_msa_filepath": f"/work/{_rel_path(run_dir, target_msa)}",
        "target_pdb_filepath_by_chain": [f"/work/{path}" for path in target_pdb_by_chain],
        "target_fasta_filepath_by_chain": [f"/work/{path}" for path in target_fasta_by_chain],
        "target_msa_filepath_by_chain": target_msa_by_chain,
        "target_chain_and_residues": chain_tags,
        "target_interface_residues": interface_residues,
        "binder_min_length": min_length,
        "binder_max_length": max_length,
        "other": {
            "source_target_chains": source_chains,
            "genie3_chain_map": chain_map,
            "genie3_residue_map": residue_map,
        },
    }
    write_json(dataset_dir / "problems" / f"{problem_key}.json", problem)
    write_json(
        run_dir / "artifacts" / "raw" / "genie3" / "input_mapping.json",
        {
            "problem_key": problem_key,
            "source_target_chains": source_chains,
            "genie3_chain_map": chain_map,
            "genie3_residue_map": residue_map,
            "requested_hotspots": [token for token in hotspots.split(",") if token],
            "mapped_hotspots": mapped_hotspots,
        },
    )
    return dataset_dir, problem_key, problem


def _yaml_scalar(value: object) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)):
        return str(value)
    text = str(value)
    if re.fullmatch(r"[A-Za-z0-9_./:-]+", text):
        return text
    return json.dumps(text)


def _write_genie3_experiment_yaml(
    path: Path,
    *,
    experiment_name: str,
    rootdir: str,
    dataset: str,
    problem_key: str,
    num_designs: int,
    seed: int,
    num_devices: int,
    direction_scale: float,
    inverse_folding_num_seq: int,
    folding_mode: str,
    folding_model_name: str,
    folding_num_models: int,
    folding_num_recycles: int,
    compile_generation: bool,
    run_mode: str,
    cond_strategy: str,
    enable_beam_search: bool,
    beam_width: int,
) -> None:
    lines = [
        "experiment:",
        f"  name: {_yaml_scalar(experiment_name)}",
        f"  seed: {seed}",
        "",
        "paths:",
        f"  rootdir: {_yaml_scalar(rootdir)}",
        f"  dataset: {_yaml_scalar(dataset)}",
        "",
        "runtime:",
        f"  num_devices: {num_devices}",
        "",
        "generation:",
        f"  compile: {_yaml_scalar(compile_generation)}",
        "  dataset:",
        "    source: target",
        f"    selections: {_yaml_scalar(problem_key)}",
        f"    n_sample: {num_designs}",
        f"    cond_strategy: {_yaml_scalar(cond_strategy)}",
    ]
    if enable_beam_search:
        lines.extend(
            [
                "  inference:",
                "    sampler:",
                "      sampler:",
                f"        direction_scale: {direction_scale}",
                "    search:",
                "      name: beam",
                "      search:",
                f"        beam_width: {beam_width}",
                "    reward:",
                "      name: colabfold",
            ]
        )
    else:
        lines.extend(
            [
                "  sampler:",
                "    sampler:",
                f"      direction_scale: {direction_scale}",
            ]
        )
    if run_mode == "full_vanilla_pipeline":
        lines.extend(
            [
                "",
                "evaluation:",
                "  version: binder",
                "  inverse_folding:",
                "    model_name: proteinmpnn",
                f"    num_seq: {inverse_folding_num_seq}",
                "  folding:",
                f"    model_name: {_yaml_scalar(folding_model_name)}",
                f"    mode: {_yaml_scalar(folding_mode)}",
                f"    num_models: {folding_num_models}",
                f"    num_recycles: {folding_num_recycles}",
            ]
        )
    if cond_strategy in {"iter_common", "iter_common_prob"}:
        lines.extend(["", "rounds:", "  - id: round_0", f"    cond_strategy: {cond_strategy}"])
    path.write_text("\n".join(lines) + "\n")


def _genie3_results_dir(run_dir: Path, problem_key: str) -> Path:
    return run_dir / "artifacts" / "raw" / "genie3" / "output" / problem_key / "results"


def _genie3_results_dirs(run_dir: Path, problem_key: str) -> list[Path]:
    output_root = run_dir / "artifacts" / "raw" / "genie3" / "output"
    direct = output_root / problem_key / "results"
    dirs = [direct]
    dirs.extend(sorted(output_root.glob(f"round_*/{problem_key}/results")))
    return dirs


def _genie3_row_key(row: dict[str, str]) -> str:
    return str(row.get("name") or row.get("domain") or row.get("sample") or row.get("id") or "").strip()


def _genie3_structure_for_row(results_dir: Path, row: dict[str, str], successful_only: bool) -> Path | None:
    row_key = _genie3_row_key(row)
    candidates: list[Path] = []
    if successful_only:
        success_dir = results_dir / "v0_success" / "successful_complexes"
        candidates.extend([success_dir / f"{row_key}.pdb", success_dir / f"{row_key}.cif"])
    design_path = str(row.get("design_filepath") or "").strip()
    if design_path:
        path = Path(design_path)
        candidates.append(path)
        if path.is_absolute() and str(path).startswith("/work/"):
            candidates.append(results_dir.parents[5] / path.relative_to("/work"))
    for candidate in candidates:
        if candidate.exists():
            return candidate
    search_roots = [
        results_dir / "v0_success" / "successful_complexes",
        results_dir.parent / "structures" / row_key if row_key else results_dir.parent / "structures",
        results_dir.parent / "structures",
        results_dir.parent / "pdbs",
        results_dir.parent,
    ]
    for root in search_roots:
        if not root.exists():
            continue
        if row_key:
            matches = sorted([*root.glob(f"*{row_key}*.pdb"), *root.glob(f"*{row_key}*.cif")])
            if matches:
                return matches[0]
        matches = sorted([*root.glob("*.pdb"), *root.glob("*.cif")])
        if matches:
            return matches[0]
    return None


def _genie3_generated_structure_paths(output_root: Path, problem_key: str) -> list[Path]:
    evaluated = [
        *output_root.glob(f"{problem_key}/eval_shards/*/structures/*/unrelaxed_rank_*.pdb"),
        *output_root.glob(f"round_*/{problem_key}/eval_shards/*/structures/*/unrelaxed_rank_*.pdb"),
    ]
    paths = [
        *output_root.glob(f"{problem_key}/pdbs/*.pdb"),
        *output_root.glob("pdbs/*.pdb"),
        *output_root.glob("**/pdbs/*.pdb"),
    ]
    preferred = evaluated if any(path.is_file() for path in evaluated) else paths
    readable: dict[Path, Path] = {}
    for path in preferred:
        if not path.is_file():
            continue
        resolved = path.resolve()
        readable.setdefault(resolved, path)
    return sorted(readable.values())


def _normalize_genie3_candidates(run_dir: Path, params: dict, target_artifact: Path) -> list[dict]:
    problem_key = str(params.get("problem_key") or "mn_app_target")
    run_mode = str(params.get("run_mode") or "generation_only")
    output_root = run_dir / "artifacts" / "raw" / "genie3" / "output"
    results_dirs = _genie3_results_dirs(run_dir, problem_key)
    results_dir = results_dirs[0]
    rows: list[dict[str, str]] = []
    successful_only = False
    if run_mode == "full_vanilla_pipeline":
        for candidate_results_dir in results_dirs:
            success_csv = candidate_results_dir / "v0_success" / "success_info.csv"
            info_csv = candidate_results_dir / "info.csv"
            if success_csv.exists():
                success_rows = _read_csv_rows(success_csv)
                if success_rows:
                    rows = success_rows
                    results_dir = candidate_results_dir
                    successful_only = True
                    break
            if info_csv.exists():
                rows = _read_csv_rows(info_csv)
                results_dir = candidate_results_dir
                break
    def generated_genie3_candidates() -> list[dict]:
        generated_paths = _genie3_generated_structure_paths(output_root, problem_key)
        generated: list[dict] = []
        for index, structure_path in enumerate(generated_paths, start=1):
            evaluated = "eval_shards" in structure_path.parts
            binder_chains, target_chains, chain_inference = _infer_generated_chain_roles(
                target_artifact,
                structure_path,
                params.get("target_chains", []),
            )
            metrics: dict[str, Any] = {
                "result_kind": "native_evaluation_fallback" if evaluated else "generation_only",
                "target_chain_inference": chain_inference,
            }
            if evaluated:
                score_paths = sorted(structure_path.parent.glob("scores_rank_*.json"))
                if score_paths:
                    scores = read_json(score_paths[0])
                    metrics.update(
                        {
                            "complex_refolding_backend": "genie3_colabfold",
                            "ptm": scores.get("ptm"),
                            "iptm": scores.get("iptm"),
                            "max_pae": scores.get("max_pae"),
                        }
                    )
            generated.append(
                {
                    "candidate_id": f"genie3_{index:05d}",
                    "source_tool": "genie3",
                    "stage": STAGE_COMPLEX_REFOLDING if evaluated else STAGE_GENERATION_BACKBONE_SEQUENCE,
                    "target_pdb": _rel_path(run_dir, target_artifact),
                    "complex_pdb": _rel_path(run_dir, structure_path),
                    "binder_pdb": None,
                    "binder_sequence": "".join(_pdb_sequences_by_chain(structure_path).get(chain, "") for chain in binder_chains) or None,
                    "target_chains": target_chains,
                    "binder_chains": binder_chains,
                    "hotspots": [token for token in params.get("hotspots", "").split(",") if token],
                    "binder_length": params.get("binder_length"),
                    "contig": None,
                    "metrics": metrics,
                    "raw_metadata": {
                        "result_kind": metrics["result_kind"],
                        "problem_key": problem_key,
                        "output_kind": "genie3_evaluated_structure" if evaluated else "genie3_generated_backbone",
                    },
                }
            )
        return generated

    if not rows:
        return write_candidates(run_dir, "genie3", generated_genie3_candidates())

    candidates: list[dict] = []
    for index, row in enumerate(rows, start=1):
        structure_path = _genie3_structure_for_row(results_dir, row, successful_only)
        if structure_path is None:
            continue
        binder_chains, target_chains, chain_inference = _infer_generated_chain_roles(
            target_artifact,
            structure_path,
            params.get("target_chains", []),
        )
        sequences = _sequences_by_chain(structure_path)
        metrics = {key: _numeric_or_text(value) for key, value in row.items() if value not in {"", None}}
        pass_filters = True if successful_only else None
        metrics.update(
            {
                "complex_refolding_backend": "genie3_vanilla",
                "result_kind": "native_pipeline",
                "target_chain_inference": chain_inference,
                "native_final_rank": index,
                "native_pass_filters": pass_filters,
                "pass_filters": pass_filters,
                "binder_plddt": metrics.get("binder_plddt") or metrics.get("avg_binder_plddt") or metrics.get("avg_plddt"),
                "iptm": metrics.get("iptm") or metrics.get("binder_ptm"),
                "ipae": metrics.get("min_interaction_pae") or metrics.get("min_interface_pae"),
                "binder_rmsd": metrics.get("complex_scrmsd") or metrics.get("binder_scrmsd"),
            }
        )
        candidates.append(
            {
                "candidate_id": f"genie3_{index:05d}",
                "source_tool": "genie3",
                "stage": STAGE_COMPLEX_REFOLDING,
                "target_pdb": _rel_path(run_dir, target_artifact),
                "complex_pdb": _rel_path(run_dir, structure_path),
                "binder_pdb": None,
                "binder_sequence": row.get("binder_seq") or "".join(sequences.get(chain, "") for chain in binder_chains) or None,
                "target_chains": target_chains,
                "binder_chains": binder_chains,
                "hotspots": [token for token in params.get("hotspots", "").split(",") if token],
                "binder_length": params.get("binder_length"),
                "contig": None,
                "metrics": metrics,
                "raw_metadata": {
                    "result_kind": "native_pipeline",
                    "problem_key": problem_key,
                    "result_csv": _rel_path(run_dir, results_dir / ("v0_success/success_info.csv" if successful_only else "info.csv")),
                    "input_target_chains": params.get("target_chains", []),
                    "target_chain_inference": chain_inference,
                    "output_kind": "genie3_v0_success" if successful_only else "genie3_info",
                },
            }
        )
    merged = _merge_candidate_pools(
        tool="genie3",
        native_candidates=candidates,
        raw_candidates=generated_genie3_candidates(),
        native_coverage="Genie3 native evaluation rows and success table",
        raw_coverage="generated Genie3 PDBs before native success filtering",
    )
    return write_candidates(run_dir, "genie3", merged)


def _pxdesign_summary_paths(run_dir: Path) -> list[Path]:
    raw_dir = run_dir / "artifacts" / "raw" / "pxdesign"
    paths = [
        *raw_dir.glob("output/design_outputs/*/summary.csv"),
        *raw_dir.glob("design_outputs/*/summary.csv"),
    ]
    return sorted(path for path in paths if path.is_file())


def _pxdesign_structure_path(summary_path: Path, row: dict[str, str]) -> Path | None:
    chosen = str(row.get("chosen_struct_path") or "").strip()
    if chosen:
        candidate = summary_path.parent / chosen
        if candidate.exists():
            return candidate
    rank = str(row.get("rank") or "").strip()
    if rank:
        for subdir in ["passing-Protenix-basic", "passing-AF2-IG-easy", "orig_designed"]:
            candidate = summary_path.parent / subdir / f"rank_{rank}.cif"
            if candidate.exists():
                return candidate
    for path in sorted(summary_path.parent.glob("**/*.cif")):
        if path.is_file():
            return path
    return None


def _pxdesign_bool(value: object) -> bool | None:
    text = str(value or "").strip().lower()
    if text in {"true", "1", "yes"}:
        return True
    if text in {"false", "0", "no"}:
        return False
    return None


def _normalize_pxdesign_candidates(run_dir: Path, params: dict, target_artifact: Path) -> list[dict]:
    candidates: list[dict] = []
    for summary_path in _pxdesign_summary_paths(run_dir):
        rows = _read_csv_rows(summary_path)
        for row in rows:
            structure_path = _pxdesign_structure_path(summary_path, row)
            if structure_path is None:
                continue
            index = len(candidates) + 1
            binder_chains, target_chains, chain_inference = _infer_generated_chain_roles(
                target_artifact,
                structure_path,
                params.get("target_chains", []),
            )
            sequences = _sequences_by_chain(structure_path)
            binder_sequence = row.get("sequence") or "".join(sequences.get(chain, "") for chain in binder_chains) or None
            metrics = {key: _numeric_or_text(value) for key, value in row.items() if value not in {"", None}}
            success_values = [
                _pxdesign_bool(row.get("Protenix-success")),
                _pxdesign_bool(row.get("Protenix-basic-success")),
                _pxdesign_bool(row.get("AF2-IG-success")),
                _pxdesign_bool(row.get("AF2-IG-easy-success")),
            ]
            native_success = any(value is True for value in success_values) if any(value is not None for value in success_values) else None
            metrics.update(
                {
                    "complex_refolding_backend": f"pxdesign_{params.get('preset', 'preview')}",
                    "target_chain_inference": chain_inference,
                    "pass_filters": native_success,
                    "result_kind": "native_pipeline",
                    "binder_plddt": metrics.get("af2_plddt"),
                    "iptm": metrics.get("af2_iptm"),
                    "ipae": metrics.get("af2_ipAE"),
                    "binder_rmsd": metrics.get("af2_bound_unbound_RMSD"),
                }
            )
            candidates.append(
                {
                    "candidate_id": f"pxdesign_{index:05d}",
                    "source_tool": "pxdesign",
                    "stage": STAGE_COMPLEX_REFOLDING,
                    "target_pdb": _rel_path(run_dir, target_artifact),
                    "complex_pdb": _rel_path(run_dir, structure_path),
                    "binder_pdb": None,
                    "binder_sequence": binder_sequence,
                    "target_chains": target_chains,
                    "binder_chains": binder_chains,
                    "hotspots": [token for token in params.get("hotspots", "").split(",") if token],
                    "binder_length": params.get("binder_length"),
                    "contig": None,
                    "metrics": metrics,
                    "raw_metadata": {
                        "result_kind": "native_pipeline",
                        "summary_csv": _rel_path(run_dir, summary_path),
                        "chosen_struct_path": row.get("chosen_struct_path"),
                        "chosen_struct_type": row.get("chosen_struct_type"),
                        "input_target_chains": params.get("target_chains", []),
                        "target_chain_inference": chain_inference,
                    },
                }
            )
    raw_candidates = _normalize_pxdesign_generation_candidates(
        run_dir,
        params,
        target_artifact,
        write_output=False,
    )
    merged = _merge_candidate_pools(
        tool="pxdesign",
        native_candidates=candidates,
        raw_candidates=raw_candidates,
        native_coverage="PXDesign summary/native retained structures",
        raw_coverage="raw PXDesign inference prediction CIFs before native filtering",
    )
    return write_candidates(run_dir, "pxdesign", merged)


def _normalize_pxdesign_generation_candidates(
    run_dir: Path,
    params: dict,
    target_artifact: Path,
    *,
    write_output: bool = True,
) -> list[dict]:
    raw_dir = run_dir / "artifacts" / "raw" / "pxdesign" / "output"
    structure_paths = sorted(
        path
        for path in raw_dir.glob("**/*.cif")
        if path.is_file() and not (path.parent == raw_dir and path.name.lower() in {"target.cif", "target_binder.cif"})
    )
    candidates: list[dict] = []
    for index, structure_path in enumerate(structure_paths, start=1):
        binder_chains, target_chains, chain_inference = _infer_generated_chain_roles(
            target_artifact,
            structure_path,
            params.get("target_chains", []),
        )
        sequences = _sequences_by_chain(structure_path)
        binder_sequence = "".join(sequences.get(chain, "") for chain in binder_chains) or None
        candidates.append(
            {
                "candidate_id": f"pxdesign_{index:05d}",
                "source_tool": "pxdesign",
                "stage": STAGE_GENERATION_BACKBONE,
                "target_pdb": _rel_path(run_dir, target_artifact),
                "complex_pdb": _rel_path(run_dir, structure_path),
                "binder_pdb": None,
                "binder_sequence": binder_sequence,
                "target_chains": target_chains,
                "binder_chains": binder_chains,
                "hotspots": [token for token in params.get("hotspots", "").split(",") if token],
                "binder_length": params.get("binder_length"),
                "contig": None,
                "metrics": {
                    "target_chain_inference": chain_inference,
                    "pxdesign_mode": "generation_only",
                    "result_kind": "generation_only",
                },
                "raw_metadata": {
                    "result_kind": "generation_only",
                    "input_target_chains": params.get("target_chains", []),
                    "target_chain_inference": chain_inference,
                    "output_kind": "infer_prediction",
                },
            }
        )
    return write_candidates(run_dir, "pxdesign", candidates) if write_output else candidates


def _proteina_complexa_result_rows(run_dir: Path) -> tuple[Path | None, list[dict[str, str]]]:
    eval_root = run_dir / "artifacts/raw/proteina_complexa/evaluation_results"
    preferred = sorted(p for p in eval_root.glob("**/RAW_*_results_*_combined.csv") if "transposed" not in p.name)
    if preferred:
        return preferred[0], _read_csv_rows(preferred[0])
    fallback = sorted(p for pattern in ("**/motif_binder_results_*.csv","**/binder_results_*.csv") for p in eval_root.glob(pattern) if "transposed" not in p.name)
    if fallback:
        return fallback[0], _read_csv_rows(fallback[0])
    return None, []


def _proteina_complexa_success_sets(run_dir):
    """Native success tables are authoritative, including joint AME criteria."""
    tables = sorted((run_dir/"artifacts/raw/proteina_complexa/evaluation_results").glob("**/all_successes_*.csv"))
    result={}
    for table in tables:
        sequence=next((seq for seq in ("mpnn_fixed","mpnn","self") if table.stem.endswith("_"+seq)),None)
        if sequence:
            for row in _read_csv_rows(table):
                if row.get("pdb_path"):result.setdefault(row["pdb_path"],set()).add(sequence)
    return bool(tables),result


def _proteina_complexa_selected_row(csv_path, row, sequence, task):
    """Select a passing redesign using the native exported criteria and list cells."""
    import ast
    import math
    import operator

    class NativeConstants(ast.NodeTransformer):
        def visit_Name(self,node):
            constants={"nan":float("nan"),"inf":float("inf")}
            return ast.copy_location(ast.Constant(constants[node.id]),node) if node.id in constants else node

    def values(value):
        if isinstance(value,str) and value.strip().startswith("["):
            try:value=ast.literal_eval(NativeConstants().visit(ast.parse(value,mode="eval")))
            except (ValueError,SyntaxError):return []
        return value if isinstance(value,list) else [value]

    paths=values(row.get(sequence+"_complex_pdb_path_all",row.get(sequence+"_complex_pdb_path","")))
    if len(paths)<=1:
        return row,0,True
    criteria_name={"protein_binder":"protein_binder","ligand_binder":"ligand_binder","motif_ligand":"motif_binder_motif_ligand_binder"}[task]
    criteria_path=csv_path.parent/("success_criteria_"+criteria_name+".json")
    try:criteria=json.loads(criteria_path.read_text())
    except (OSError,ValueError):return row,None,False
    rules=[]
    for metric,spec in criteria.get("thresholds",criteria.get("binder_thresholds",{})).items():
        rules.append((f"{sequence}_{spec.get('column_prefix','complex')}_{metric}_all",spec))
    rules.extend((spec["column"].format(seq_type=sequence),spec) for spec in criteria.get("motif_criteria",[]))
    if not rules:return row,None,False
    operations={"<":operator.lt,"<=":operator.le,">":operator.gt,">=":operator.ge,"==":operator.eq,"!=":operator.ne}
    def passes(index):
        for column,spec in rules:
            cell=values(row.get(column,row.get(column.removesuffix("_all"))))
            # Native scalar motif criteria are broadcast across redesigns.
            value=cell[index] if index<len(cell) else cell[0] if len(cell)==1 else None
            try:
                number=float(value)*float(spec.get("scale",1))
                if not math.isfinite(number) or not operations[spec.get("op","<=")](number,float(spec["threshold"])):return False
            except (TypeError,ValueError,KeyError):return False
        return True
    passing=[index for index in range(len(paths)) if passes(index)]
    if not passing:return row,None,False
    best_path=row.get(sequence+"_complex_pdb_path")
    best_index=paths.index(best_path) if best_path in paths else None
    index=best_index if best_index in passing else passing[0]
    selected=dict(row)
    for key,value in row.items():
        if key.startswith(sequence+"_") and key.endswith("_all"):
            items=values(value)
            if index<len(items):selected[key[:-4]]=json.dumps(items[index]) if isinstance(items[index],(list,dict)) else items[index]
    return selected,index,True


def _proteina_complexa_resolve_path(run_dir: Path, csv_path: Path | None, path_text: str) -> Path | None:
    text = str(path_text or "").strip().strip("'\"")
    if not text:
        return None
    path = Path(text)
    if text.startswith("/work/"):
        path = run_dir / text[len("/work/"):]
        return path if path.is_file() else None
    if path.is_absolute():
        return path if path.is_file() else None
    raw_root = run_dir / "artifacts" / "raw" / "proteina_complexa"
    candidates = [raw_root / path]
    if text.startswith("./"):
        candidates.append(raw_root / text[2:])
    if csv_path is not None:
        candidates.extend([csv_path.parent / path, csv_path.parent / text.lstrip("./")])
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def _proteina_complexa_native_pass(row: dict[str, str]) -> bool | None:
    ipae = _float_or_none(row.get("self_complex_i_pAE"))
    plddt = _float_or_none(row.get("self_complex_pLDDT"))
    rmsd = _float_or_none(row.get("self_binder_scRMSD_ca") or row.get("self_binder_scRMSD"))
    if ipae is None or plddt is None or rmsd is None:
        return None
    return ipae * 31.0 <= 7.0 and plddt >= 0.9 and rmsd < 1.5


def _proteina_complexa_ranked_rows(rows: list[dict[str, str]]) -> list[tuple[int, dict[str, str], str | None, float | None]]:
    if any(_float_or_none(row.get("self_complex_i_pAE")) is not None for row in rows):
        ranked = sorted(
            enumerate(rows),
            key=lambda item: (
                _float_or_none(item[1].get("self_complex_i_pAE")) is None,
                _float_or_none(item[1].get("self_complex_i_pAE")) or 0.0,
                item[0],
            ),
        )
        return [
            (rank, row, "self_complex_i_pAE", _float_or_none(row.get("self_complex_i_pAE")))
            for rank, (_original_index, row) in enumerate(ranked, start=1)
        ]
    if any(_float_or_none(row.get("self_complex_i_pTM")) is not None for row in rows):
        ranked = sorted(
            enumerate(rows),
            key=lambda item: (
                _float_or_none(item[1].get("self_complex_i_pTM")) is None,
                -(_float_or_none(item[1].get("self_complex_i_pTM")) or 0.0),
                item[0],
            ),
        )
        return [
            (rank, row, "self_complex_i_pTM", _float_or_none(row.get("self_complex_i_pTM")))
            for rank, (_original_index, row) in enumerate(ranked, start=1)
        ]
    return [(rank, row, None, None) for rank, row in enumerate(rows, start=1)]


def _normalize_proteina_complexa_candidates(run_dir: Path, params: dict, target_artifact: Path) -> list[dict]:
    csv_path, rows = _proteina_complexa_result_rows(run_dir)
    has_success_tables,success_sets=_proteina_complexa_success_sets(run_dir)
    task=params.get("pipeline_task","protein_binder")
    candidates: list[dict] = []
    def roles(path):
        if task=="protein_binder":return _infer_generated_chain_roles(target_artifact,path,params.get("target_chains",[]))
        from mn_molecular_design_studio.design_validation import read_structure
        protein=list(_sequences_by_chain(path))
        all_chains=[c.id for c in next(read_structure(path).get_models())]
        return protein,[c for c in all_chains if c not in protein],"ligand_target"
    for index, row, rank_metric, rank_value in _proteina_complexa_ranked_rows(rows):
        successful=success_sets.get(row.get("pdb_path",""),set())
        sequence_type=next((seq for seq in ("self","mpnn","mpnn_fixed") if seq in successful),"self")
        selected_index=None;selection_verified=True
        if successful:
            selection_verified=False
            for sequence in (seq for seq in ("self","mpnn","mpnn_fixed") if seq in successful):
                selected,redesign_index,verified=_proteina_complexa_selected_row(csv_path,row,sequence,task)
                if verified:
                    row,sequence_type,selected_index,selection_verified=selected,sequence,redesign_index,True
                    break
        evaluated_path=_proteina_complexa_resolve_path(run_dir,csv_path,row.get(sequence_type+"_complex_pdb_path",""))
        structure_path = (
            evaluated_path
            or _proteina_complexa_resolve_path(run_dir, csv_path, row.get("complex_pdb_path", ""))
            or _proteina_complexa_resolve_path(run_dir, csv_path, row.get("pdb_path", ""))
        )
        if structure_path is None:
            continue
        binder_chains,target_chains,chain_inference=roles(structure_path)
        sequences = _sequences_by_chain(structure_path)
        binder_sequence = row.get(sequence_type+"_binder_sequence") or row.get(sequence_type+"_sequence") or row.get("binder_sequence") or row.get("_res_mpnn_best_sequence") or "".join(
            sequences.get(chain, "") for chain in binder_chains
        ) or None
        metrics = {key: _numeric_or_text(value) for key, value in row.items() if value not in {"", None}}
        native_pass = bool(successful) if has_success_tables else _proteina_complexa_native_pass(row) if task=="protein_binder" else None
        selected_pass=native_pass if not native_pass or selection_verified and evaluated_path else None
        native_ipae = _float_or_none(row.get(sequence_type+"_complex_i_pAE"))
        metrics.update(
            {
                "complex_refolding_backend": "proteina_complexa_native",
                "result_kind": "native_pipeline",
                "target_chain_inference": chain_inference,
                "binder_plddt": metrics.get(sequence_type+"_complex_pLDDT"),
                "iptm": metrics.get(sequence_type+"_complex_i_pTM"),
                "ipae": metrics.get(sequence_type+"_complex_i_pAE"),
                "ipsae": metrics.get(sequence_type+"_complex_min_ipSAE"),
                "ipsae_min": metrics.get(sequence_type+"_complex_min_ipSAE"),
                "ipsae_max": metrics.get(sequence_type+"_complex_max_ipSAE"),
                "ipsae_avg": metrics.get(sequence_type+"_complex_avg_ipSAE"),
                "binder_rmsd": next((metrics[key] for suffix in ("_binder_scRMSD_ca","_binder_scRMSD","_binder_scRMSD_bb3") if (key:=sequence_type+suffix) in metrics),None),
                "native_final_rank": index,
                "native_rank_metric": rank_metric,
                "native_rank_value": rank_value,
                "native_ipae_scaled": native_ipae * 31.0 if native_ipae is not None else None,
                "native_design_pass_filters": native_pass,
                "native_pass_filters": selected_pass,
                "pass_filters": selected_pass,
                "proteina_complexa_success_criteria": "Native task-specific success tables" if has_success_tables else "Protein AlphaProteo criteria" if task=="protein_binder" else "Native success tables unavailable",
                "pipeline_task":task,"evaluated_sequence_type":sequence_type,
                "evaluated_redesign_index":selected_index,
                "passing_redesign_verified":selection_verified if successful else None,
            }
        )
        candidates.append(
            {
                "candidate_id": f"proteina_complexa_{index:05d}",
                "source_tool": "proteina_complexa",
                "stage": STAGE_COMPLEX_REFOLDING,
                "target_pdb": _rel_path(run_dir, target_artifact),
                "complex_pdb": _rel_path(run_dir, structure_path),
                "binder_pdb": None,
                "binder_sequence": binder_sequence,
                "target_chains": target_chains,
                "binder_chains": binder_chains,
                "hotspots": [token for token in params.get("hotspots", "").split(",") if token],
                "binder_length": params.get("binder_length"),
                "contig": None,
                "metrics": metrics,
                "raw_metadata": {
                    "result_kind": "native_pipeline",
                    "result_csv": _rel_path(run_dir, csv_path),
                    "complex_pdb_path": row.get("complex_pdb_path"),
                    "self_complex_pdb_path": row.get("self_complex_pdb_path"),
                    "evaluated_sequence_type": sequence_type,
                    "evaluated_complex_pdb_path": row.get(sequence_type+"_complex_pdb_path"),
                    "input_target_chains": params.get("target_chains", []),
                    "target_chain_inference": chain_inference,
                    "output_kind": "proteina_complexa_native",
                },
            }
        )
    raw_candidates: list[dict] = []
    seen_sequences: set[str] = set()
    inference_root = run_dir / "artifacts" / "raw" / "proteina_complexa" / "inference"
    raw_limit = int(params.get("num_designs") or 0) if bool(params.get("generator_only", False)) else 0
    for structure_path in sorted(inference_root.glob("**/job_*/*.pdb")):
        if raw_limit and len(raw_candidates) >= raw_limit:
            break
        if structure_path.name.endswith("_binder.pdb") or structure_path.name.endswith("_updated.pdb"):
            continue
        index = len(raw_candidates) + 1
        binder_chains,target_chains,chain_inference=roles(structure_path)
        sequences = _sequences_by_chain(structure_path)
        binder_sequence = "".join(sequences.get(chain, "") for chain in binder_chains) or None
        if params.get("generator_only") and params.get("deduplicate_sequences", True) and binder_sequence:
            if binder_sequence in seen_sequences:
                continue
            seen_sequences.add(binder_sequence)
        generation_metadata = read_json(structure_path.with_suffix(".studio.json"))
        raw_candidates.append(
            {
                "candidate_id": f"proteina_complexa_generated_{index:05d}",
                "source_tool": "proteina_complexa",
                "stage": STAGE_GENERATION_BACKBONE_SEQUENCE,
                "target_pdb": _rel_path(run_dir, target_artifact),
                "complex_pdb": _rel_path(run_dir, structure_path),
                "binder_pdb": None,
                "binder_sequence": binder_sequence,
                "target_chains": target_chains,
                "binder_chains": binder_chains,
                "hotspots": [token for token in params.get("hotspots", "").split(",") if token],
                "binder_length": params.get("binder_length"),
                "contig": None,
                "metrics": {
                    "result_kind": "generated" if params.get("generator_only") else "generated_prefilter",
                    "pipeline_task": task,
                    "pass_filters": None,
                    "target_chain_inference": chain_inference,
                },
                "raw_metadata": {
                    **generation_metadata,
                    "result_kind": "generated" if params.get("generator_only") else "generated_prefilter",
                    "generator_only": bool(params.get("generator_only")),
                    "pipeline_task": task,
                    "target_entity_kind": "protein" if task == "protein_binder" else "ligand",
                    "ligand_original_code": (params.get("ligand_residue") or {}).get("name"),
                    "input_target_chains": params.get("target_chains", []),
                    "target_chain_inference": chain_inference,
                    "output_kind": "proteina_complexa_generation" if params.get("generator_only") else "proteina_complexa_inference_prefilter",
                },
            }
        )
    if bool(params.get("generator_only", False)):
        return write_candidates(run_dir, "proteina_complexa", raw_candidates)
    merged = _merge_candidate_pools(
        tool="proteina_complexa",
        native_candidates=candidates,
        raw_candidates=raw_candidates,
        native_coverage="ranked Proteina-Complexa evaluation/filter result rows",
        raw_coverage="generated Proteina-Complexa inference complexes before evaluation filtering",
    )
    return write_candidates(run_dir, "proteina_complexa", merged)


def _normalize_rfdiffusion3_contig(contig: str) -> str:
    normalized = contig.strip()
    if not normalized:
        return normalized
    if ",/0," in normalized:
        return normalized
    return re.sub(r"/0\s+", ",/0,", normalized)


def _cif_chain_lengths(path: Path) -> dict[str, int]:
    return {chain: len(residues) for chain, residues in mmcif_protein_residues_by_chain(path).items()}


def _pdb_sequences_by_chain(path: Path) -> dict[str, str]:
    residues: dict[str, list[tuple[int, str, str]]] = {}
    seen: set[tuple[str, str, str]] = set()
    for line in path.read_text(errors="ignore").splitlines():
        if not line.startswith("ATOM  "):
            continue
        chain = line[21].strip() or "_"
        key = (chain, line[22:26].strip(), line[26].strip())
        if key in seen:
            continue
        seen.add(key)
        try:
            residue_number = int(line[22:26])
        except ValueError:
            residue_number = 0
        residues.setdefault(chain, []).append(
            (residue_number, line[26].strip(), AA3_TO_1.get(line[17:20].strip().upper(), "X"))
        )
    return {
        chain: "".join(item[2] for item in sorted(chain_residues, key=lambda item: (item[0], item[1])))
        for chain, chain_residues in residues.items()
    }


def _cif_sequences_by_chain(path: Path) -> dict[str, str]:
    return {
        chain: "".join(AA3_TO_1.get(residue_name, "X") for _number, _insertion, residue_name in residues)
        for chain, residues in mmcif_protein_residues_by_chain(path).items()
    }


def _sequences_by_chain(path: Path) -> dict[str, str]:
    if path.suffix.lower() in {".cif", ".mmcif"} or path.name.endswith((".cif.gz", ".mmcif.gz")):
        return _cif_sequences_by_chain(path)
    return _pdb_sequences_by_chain(path)


def _target_chains_in_generated_structure(
    target_artifact: Path,
    generated_structure: Path,
    requested_target_chains: list[str],
    output_chains: list[str],
) -> list[str]:
    requested = [chain for chain in requested_target_chains if chain in output_chains]
    if not target_artifact.exists():
        return requested
    target_sequences = _sequences_by_chain(target_artifact)
    if requested_target_chains:
        target_sequences={chain:sequence for chain,sequence in target_sequences.items() if chain in requested_target_chains}
    generated_sequences = _sequences_by_chain(generated_structure)
    target_chains: list[str] = []
    used_target_sequences: set[str] = set()
    for output_chain, output_sequence in generated_sequences.items():
        if not output_sequence:
            continue
        for target_chain, target_sequence in target_sequences.items():
            if target_chain in used_target_sequences or not target_sequence:
                continue
            if output_sequence == target_sequence:
                target_chains.append(output_chain)
                used_target_sequences.add(target_chain)
                break
    return target_chains or requested


def _chain_ids_from_structure(path: Path) -> list[str]:
    if path.suffix.lower() in {".cif", ".mmcif"} or path.name.endswith((".cif.gz", ".mmcif.gz")):
        return list(_cif_chain_lengths(path))
    return list(_pdb_sequences_by_chain(path))


def _infer_generated_chain_roles(target_artifact: Path | None, generated_structure: Path, requested_target_chains: list[str]) -> tuple[list[str], list[str], str]:
    output_chains = _chain_ids_from_structure(generated_structure)
    requested = [chain for chain in requested_target_chains if chain in output_chains]
    target_chains = (
        _target_chains_in_generated_structure(
            target_artifact,
            generated_structure,
            requested_target_chains,
            output_chains,
        )
        if target_artifact is not None
        else []
    )
    inference = "sequence_match" if target_chains and target_artifact is not None and target_artifact.exists() else "requested_chain_ids"
    if not target_chains:
        if "B" in output_chains:
            target_chains = ["B"]
            inference = "fallback_B"
        elif len(output_chains) > 1:
            target_chains = output_chains[1:]
            inference = "fallback_nonfirst"
    target_set = set(target_chains)
    binder_chains = [chain for chain in output_chains if chain not in target_set]
    if not binder_chains and output_chains:
        binder_chains = [output_chains[0]]
        target_chains = [chain for chain in output_chains if chain not in set(binder_chains)]
        inference = "fallback_first_binder"
    return binder_chains, target_chains, inference


def _normalize_rfdiffusion3_candidates(
    run_dir: Path,
    params: dict,
    target_artifact: Path,
    *,
    write_output: bool = True,
) -> list[dict]:
    raw_dir = run_dir / "artifacts" / "raw" / "rfdiffusion3_foundry" / "rfd3"
    structure_paths = sorted([*raw_dir.glob("*.cif"), *raw_dir.glob("*.cif.gz")])
    metadata_paths = {path.stem.replace(".cif", ""): path for path in raw_dir.glob("*.json")}
    target_chains = params.get("target_chains", [])
    candidates: list[dict] = []
    for index, structure_path in enumerate(structure_paths, start=1):
        chain_lengths = _cif_chain_lengths(structure_path)
        output_chains = list(chain_lengths)
        binder_chains, inferred_target_chains, inference = _infer_generated_chain_roles(target_artifact, structure_path, target_chains)
        generated_binder_present = bool(binder_chains)
        metadata_path = metadata_paths.get(structure_path.stem.replace(".cif", ""))
        candidates.append(
            {
                "candidate_id": f"rfdiffusion3_foundry_{index:05d}",
                "source_tool": "rfdiffusion3_foundry",
                "stage": STAGE_GENERATION_BACKBONE,
                "target_pdb": _rel_path(run_dir, target_artifact),
                "complex_pdb": _rel_path(run_dir, structure_path),
                "binder_pdb": None,
                "binder_sequence": None,
                "target_chains": inferred_target_chains,
                "binder_chains": binder_chains,
                "hotspots": [token for token in params.get("hotspots", "").split(",") if token],
                "binder_length": params.get("binder_length"),
                "contig": params.get("contig"),
                "metrics": {
                    "chain_lengths": chain_lengths,
                    "generated_binder_present": generated_binder_present,
                },
                "raw_metadata": {
                    "metadata_path": _rel_path(run_dir, metadata_path),
                    "output_chains": output_chains,
                    "input_target_chains": target_chains,
                    "target_chain_inference": inference,
                    "normalization_warning": "" if generated_binder_present else "RFdiffusion3 output has no non-target binder chain.",
                },
            }
        )
    if write_output:
        _write_candidates(run_dir, "rfdiffusion3_foundry", candidates)
    return candidates


def _rf3_summary_metrics(summary_path: Path | None) -> dict:
    if summary_path is None or not summary_path.exists():
        return {}
    payload = read_json(summary_path)
    metrics: dict = {}
    for key in [
        "ranking_score",
        "confidence",
        "ptm",
        "iptm",
        "has_clash",
        "fraction_disordered",
    ]:
        if key in payload:
            metrics[key] = payload[key]
    pae = payload.get("pae") or payload.get("predicted_aligned_error")
    if isinstance(pae, list):
        metrics["rf3_pae_available"] = True
    return metrics


def _find_rf3_summary(rf3_out_dir: Path) -> Path | None:
    summary_paths = sorted(rf3_out_dir.glob("**/*_summary_confidences.json"))
    return summary_paths[0] if summary_paths else None


def _write_foundry_cif_with_target_msas(
    source_cif: Path,
    staged_cif: Path,
    target_artifact: Path,
    output_target_chains: list[str],
    source_target_chains: list[str],
    msa_by_source_chain: dict[str, str],
) -> dict[str, str]:
    text = gzip.open(source_cif, "rt", errors="ignore").read() if source_cif.name.endswith(".gz") else source_cif.read_text(errors="ignore")
    generated_sequences = _sequences_by_chain(source_cif)
    target_sequences = _pdb_sequences_by_chain(target_artifact)
    msa_by_output_chain: dict[str, str] = {}
    for output_chain in output_target_chains:
        output_sequence = generated_sequences.get(output_chain, "")
        matched_source_chain = ""
        for source_chain in source_target_chains:
            if output_sequence and output_sequence == target_sequences.get(source_chain, ""):
                matched_source_chain = source_chain
                break
        if not matched_source_chain and output_chain in msa_by_source_chain:
            matched_source_chain = output_chain
        msa_path = msa_by_source_chain.get(matched_source_chain)
        if msa_path:
            msa_by_output_chain[output_chain] = msa_path

    header_lines = ["#"]
    for chain_id, msa_path in sorted(msa_by_output_chain.items()):
        header_lines.append(f"_msa_paths_by_chain_id.{chain_id}   {msa_path}")
    header_lines.append("#")
    staged_cif.parent.mkdir(parents=True, exist_ok=True)
    staged_cif.write_text("\n".join(header_lines) + "\n" + text)
    return msa_by_output_chain


def _normalize_rfdiffusion3_foundry_native_candidates(run_dir: Path, params: dict, target_artifact: Path) -> list[dict]:
    raw_root = run_dir / "artifacts" / "raw" / "rfdiffusion3_foundry"
    rf3_root = raw_root / "rf3"
    mapping_rows = _read_csv_rows(raw_root / "foundry_native_mapping.tsv")
    candidates: list[dict] = []
    for index, row in enumerate(mapping_rows, start=1):
        source_id = row.get("rfd3_id") or f"rfd3_{index:05d}"
        mpnn_id = row.get("mpnn_id") or f"mpnn_{index:05d}"
        rf3_out_dir = Path(row.get("rf3_out_dir") or "")
        if not rf3_out_dir.is_absolute():
            rf3_out_dir = run_dir / rf3_out_dir
        model_paths = sorted(rf3_out_dir.glob("*_model.cif"))
        if not model_paths:
            model_paths = sorted(rf3_out_dir.glob("**/*_model.cif"))
        if not model_paths:
            continue
        structure_path = model_paths[0]
        summary_path = _find_rf3_summary(rf3_out_dir)
        rfd3_path = run_dir / row.get("rfd3_cif", "")
        mpnn_path = run_dir / row.get("mpnn_cif", "")
        binder_chains, target_chains, chain_inference = _infer_generated_chain_roles(
            target_artifact,
            structure_path,
            params.get("target_chains", []),
        )
        sequences = _sequences_by_chain(structure_path)
        binder_sequence = "".join(sequences.get(chain, "") for chain in binder_chains) or None
        metrics = {
            **_rf3_summary_metrics(summary_path),
            "complex_refolding_backend": "rf3",
            "sequence_design_backend": "foundry_mpnn",
            "generation_backend": "rfdiffusion3",
            "result_kind": "native_pipeline",
            "native_final_rank": index,
        }
        candidates.append(
            {
                "candidate_id": f"rfdiffusion3_foundry_{index:05d}_rf3",
                "source_tool": "rfdiffusion3_foundry",
                "stage": STAGE_COMPLEX_REFOLDING,
                "target_pdb": _rel_path(run_dir, target_artifact),
                "complex_pdb": _rel_path(run_dir, structure_path),
                "binder_pdb": None,
                "binder_sequence": binder_sequence,
                "target_chains": target_chains,
                "binder_chains": binder_chains,
                "hotspots": [token for token in params.get("hotspots", "").split(",") if token],
                "binder_length": params.get("binder_length"),
                "contig": params.get("contig"),
                "metrics": metrics,
                "raw_metadata": {
                    "result_kind": "native_pipeline",
                    "rfd3_cif": _rel_path(run_dir, rfd3_path) if rfd3_path.exists() else row.get("rfd3_cif"),
                    "mpnn_cif": _rel_path(run_dir, mpnn_path) if mpnn_path.exists() else row.get("mpnn_cif"),
                    "rf3_summary": _rel_path(run_dir, summary_path),
                    "input_target_chains": params.get("target_chains", []),
                    "target_chain_inference": chain_inference,
                    "target_msa_by_rf3_chain": row.get("target_msa_by_rf3_chain", ""),
                },
            }
        )
    raw_candidates = _normalize_rfdiffusion3_candidates(
        run_dir,
        params,
        target_artifact,
        write_output=False,
    )
    mpnn_candidates: list[dict] = []
    for mpnn_path in sorted(raw_root.glob("mpnn/**/*.cif")):
        index = len(mpnn_candidates) + 1
        binder_chains, target_chains, chain_inference = _infer_generated_chain_roles(
            target_artifact,
            mpnn_path,
            params.get("target_chains", []),
        )
        sequences = _sequences_by_chain(mpnn_path)
        mpnn_candidates.append(
            {
                "candidate_id": f"rfdiffusion3_foundry_mpnn_{index:05d}",
                "source_tool": "rfdiffusion3_foundry",
                "stage": STAGE_GENERATION_BACKBONE_SEQUENCE,
                "target_pdb": _rel_path(run_dir, target_artifact),
                "complex_pdb": _rel_path(run_dir, mpnn_path),
                "binder_pdb": None,
                "binder_sequence": "".join(sequences.get(chain, "") for chain in binder_chains) or None,
                "target_chains": target_chains,
                "binder_chains": binder_chains,
                "hotspots": [token for token in params.get("hotspots", "").split(",") if token],
                "binder_length": params.get("binder_length"),
                "contig": params.get("contig"),
                "metrics": {
                    "result_kind": "foundry_mpnn_prefilter",
                    "target_chain_inference": chain_inference,
                },
                "raw_metadata": {
                    "result_kind": "foundry_mpnn_prefilter",
                    "input_target_chains": params.get("target_chains", []),
                    "target_chain_inference": chain_inference,
                    "output_kind": "foundry_mpnn_prefilter",
                },
            }
        )
    merged = _merge_candidate_pools(
        tool="rfdiffusion3_foundry",
        native_candidates=candidates,
        raw_candidates=[*mpnn_candidates, *raw_candidates],
        native_coverage="RF3-folded native Foundry survivors",
        raw_coverage="RFdiffusion3 and Foundry MPNN structures before RF3 survivor filtering",
    )
    _write_candidates(run_dir, "rfdiffusion3_foundry", merged)
    return merged


def _finish_design_job(run_dir: Path, success: bool, rc: int, candidates: list[dict], artifact_patterns: list[tuple[str, str]]) -> None:
    artifacts = _collect_artifacts(
        run_dir,
        run_dir / "artifacts",
        artifact_patterns
        + [
            ("normalized_candidates/*.jsonl", "normalized_candidates"),
            ("normalized_candidates/*.json", "campaign_result"),
        ],
    )
    finish_job(
        run_dir,
        success and bool(candidates),
        {
            "outputs": {"artifacts": artifacts, "candidates": candidates},
            "metrics": {"return_code": rc, "artifact_count": len(artifacts), "candidate_count": len(candidates)},
            "downstream_artifacts": {
                "candidates_jsonl": "artifacts/normalized_candidates/candidates.jsonl",
                "campaign_result": "artifacts/normalized_candidates/campaign_result.json",
            },
        },
    )


def _copy_target_for_design(job_run_dir: Path, target_pdb: Path, target_chains: list[str], tool_dir: str = "input") -> Path:
    target_artifact = artifact_path(job_run_dir, tool_dir, "target.pdb")
    target_artifact.write_text(filter_pdb_text(target_pdb.read_text(errors="ignore"), keep_chains=set(target_chains) or None))
    return target_artifact


def _run_shell_steps(run_dir: Path, steps: list[dict]) -> int:
    steps = apply_docker_cpu_limits_to_steps(run_dir, steps)
    write_json(run_dir / "command.json", {"mode": "docker", "steps": steps})
    update_status(run_dir, "running")
    with (run_dir / "stdout.log").open("a") as stdout, (run_dir / "stderr.log").open("a") as stderr:
        for step in steps:
            stdout.write(f"$ {' '.join(step['command'])}\n")
            stdout.flush()
            proc = subprocess.run(step["command"], stdout=stdout, stderr=stderr, check=False)
            if proc.returncode != 0:
                return int(proc.returncode)
    return 0


def _collect_artifacts(run_dir: Path, root: Path, patterns: list[tuple[str, str]]) -> list[dict]:
    artifacts: list[dict] = []
    for pattern, artifact_type in patterns:
        for path_text in glob.glob(str(root / pattern), recursive=True):
            path = Path(path_text)
            if path.is_file():
                artifacts.append(Artifact(path.stem, path, artifact_type).to_json(run_dir))
    return artifacts


def _candidate_record(
    run_dir: Path,
    candidate_id: str,
    complex_pdb: Path,
    trb_path: Path | None,
    params: dict,
    target_artifact: Path,
) -> dict:
    def rel(path: Path | None) -> str | None:
        if path is None:
            return None
        try:
            return str(path.relative_to(run_dir))
        except ValueError:
            return str(path)

    binder_chains, target_chains, chain_inference = _infer_generated_chain_roles(
        target_artifact,
        complex_pdb,
        params.get("target_chains", []),
    )
    binder_chain_for_metrics = binder_chains[0] if binder_chains else "A"
    target_chain_for_metrics = target_chains[0] if target_chains else "B"
    metrics = calculate_hotspot_metrics(
        complex_pdb,
        params.get("hotspots", ""),
        binder_chain=binder_chain_for_metrics,
        target_chain=target_chain_for_metrics,
        contact_cutoff=float(params.get("hotspot_contact_cutoff") or 8.0),
        contig=str(params.get("contig") or ""),
    )
    metrics["binder_chains"] = ",".join(binder_chains)
    metrics["target_chains"] = ",".join(target_chains)
    metrics["target_chain_inference"] = chain_inference
    if params.get("apply_backbone_hotspot_prefilter"):
        passed, failures = passes_hotspot_prefilter(
            metrics,
            min_contact_fraction=float(params.get("backbone_min_hotspot_contact_fraction") or 0.25),
            max_min_distance=float(params.get("backbone_max_hotspot_distance") or 10.0),
        )
        metrics["passes_backbone_hotspot_filter"] = passed
        metrics["backbone_hotspot_filter_failures"] = "; ".join(failures)
    else:
        metrics["passes_backbone_hotspot_filter"] = None
        metrics["backbone_hotspot_filter_failures"] = ""

    return {
        "candidate_id": candidate_id,
        "source_tool": "rfdiffusion_classic",
        "stage": STAGE_GENERATION_BACKBONE,
        "target_pdb": rel(target_artifact),
        "complex_pdb": rel(complex_pdb),
        "binder_pdb": None,
        "binder_sequence": None,
        "target_chains": target_chains,
        "binder_chains": binder_chains,
        "hotspots": [token for token in params.get("hotspots", "").split(",") if token],
        "binder_length": params.get("binder_length"),
        "contig": params.get("contig"),
        "metrics": metrics,
        "raw_metadata": {
            "trb_path": rel(trb_path),
            "input_target_chains": params.get("target_chains", []),
            "output_chains": _chain_ids_from_structure(complex_pdb),
            "target_chain_inference": chain_inference,
            "chain_convention": "dynamic_target_sequence_match_with_classic_fallback",
        },
    }


def _normalize_rfdiffusion_candidates(run_dir: Path, params: dict, target_artifact: Path) -> list[dict]:
    raw_dir = run_dir / "artifacts" / "raw" / "rfdiffusion" / "output"
    candidates: list[dict] = []
    for index, pdb_path in enumerate(sorted(raw_dir.glob("design_*.pdb")), start=1):
        trb_path = pdb_path.with_suffix(".trb")
        candidates.append(
            _candidate_record(
                run_dir,
                candidate_id=f"rfdiffusion_classic_{index:05d}",
                complex_pdb=pdb_path,
                trb_path=trb_path if trb_path.exists() else None,
                params=params,
                target_artifact=target_artifact,
            )
        )
    return write_candidates(run_dir, "rfdiffusion_classic", candidates)


def _job_source_for_run(run_dir: Path) -> dict:
    metadata = read_json(run_dir / "metadata.json")
    candidates = read_candidates(run_dir)
    return {
        "task_group": metadata.get("task_group") or run_dir.parent.name,
        "run_id": run_dir.name,
        "run_dir": str(run_dir),
        "job_code": metadata.get("job_code") or "",
        "tool": metadata.get("tool") or "",
        "candidates_jsonl": str(run_dir / "artifacts" / "normalized_candidates" / "candidates.jsonl"),
        "candidate_count": len(candidates),
        "stage_counts": candidate_stage_counts(candidates),
    }


def _mark_internal_child_run(child_run_dir: Path, parent_run_dir: Path, step_index: int) -> None:
    metadata = read_json(child_run_dir / "metadata.json")
    metadata.update(
        {
            "hidden": True,
            "internal_parent_run_id": parent_run_dir.name,
            "internal_parent_run_dir": str(parent_run_dir),
            "internal_step_index": step_index,
        }
    )
    write_json(child_run_dir / "metadata.json", metadata)


def _append_internal_child_logs(parent_run_dir: Path, child_run_dir: Path, label: str) -> None:
    for log_name in ("stdout.log", "stderr.log"):
        source_log = child_run_dir / log_name
        if not source_log.exists():
            continue
        with (parent_run_dir / log_name).open("a") as target:
            target.write(f"\n\n===== internal step: {label} ({child_run_dir.name}) =====\n")
            target.write(source_log.read_text(errors="ignore"))


def _copy_pipeline_file(parent_run_dir: Path, source_run_dir: Path, path_text: str | None, folder: Path) -> str | None:
    if not path_text:
        return None
    path = Path(path_text)
    source_path = path if path.is_absolute() else source_run_dir / path
    if not source_path.exists() or not source_path.is_file():
        return path_text
    folder.mkdir(parents=True, exist_ok=True)
    target_path = folder / source_path.name
    if source_path.resolve() != target_path.resolve():
        shutil.copy2(source_path, target_path)
    return _rel_path(parent_run_dir, target_path)


def _adopt_rfdiffusion_native_candidates(
    parent_run_dir: Path,
    complex_run_dir: Path,
    params: dict,
    target_artifact: Path,
) -> list[dict]:
    adopted: list[dict] = []
    source_candidates = [
        candidate
        for candidate in read_candidates(complex_run_dir)
        if candidate.get("stage") == STAGE_COMPLEX_REFOLDING
    ]
    native_dir = parent_run_dir / "artifacts" / "raw" / "rfdiffusion" / "native_pipeline"
    target_rel = _rel_path(parent_run_dir, target_artifact)
    for index, source in enumerate(source_candidates, start=1):
        candidate_id = f"rfdiffusion_classic_{index:05d}_native"
        candidate_dir = native_dir / candidate_id
        raw_metadata = dict(source.get("raw_metadata") or {})
        complex_rel = _copy_pipeline_file(parent_run_dir, complex_run_dir, source.get("complex_pdb"), candidate_dir)
        binder_rel = _copy_pipeline_file(parent_run_dir, complex_run_dir, source.get("binder_pdb"), candidate_dir)
        pae_rel = _copy_pipeline_file(parent_run_dir, complex_run_dir, raw_metadata.get("pae_path"), candidate_dir)
        if pae_rel:
            raw_metadata["pae_path"] = pae_rel
        raw_metadata.update(
            {
                "source_candidate": source,
                "source_complex_run_dir": str(complex_run_dir),
                "result_kind": "native_pipeline",
                "native_pipeline": "rfdiffusion_classic_ligandmpnn_monomer_complex_refolding",
            }
        )
        metrics = dict(source.get("metrics") or {})
        metrics.update(
            {
                "result_kind": "native_pipeline",
                "generation_backend": "rfdiffusion_classic",
                "sequence_design_backend": params.get("sequence_design_method") or "ligandmpnn",
            }
        )
        adopted.append(
            {
                **source,
                "candidate_id": candidate_id,
                "stage": STAGE_COMPLEX_REFOLDING,
                "source_tool": "rfdiffusion_classic",
                "tool": "rfdiffusion_classic",
                "target_pdb": target_rel,
                "complex_pdb": complex_rel or source.get("complex_pdb"),
                "binder_pdb": binder_rel,
                "target_chains": source.get("target_chains") or params.get("target_chains", []),
                "hotspots": [token for token in str(params.get("hotspots") or "").split(",") if token],
                "binder_length": params.get("binder_length"),
                "contig": params.get("contig"),
                "metrics": metrics,
                "parents": [str(source.get("candidate_id") or "")],
                "raw_metadata": raw_metadata,
            }
        )
    return write_candidates(parent_run_dir, "rfdiffusion_classic", adopted) if adopted else []


def _rfdiffusion_pipeline_payload(
    *,
    target_pdb_path: str,
    contig: str,
    hotspots: str,
    num_designs: int,
    final_run_parameters: str,
    rfdiffusion_guidance_preset: str,
    scaffoldguided: bool,
    scaffold_dir: str,
    scaffold_target_pdb: bool,
    scaffold_target_ss: str,
    scaffold_target_adj: str,
    backbone_filters: str,
    apply_backbone_hotspot_prefilter: bool,
    backbone_min_hotspot_contact_fraction: float,
    backbone_max_hotspot_distance: float,
    mpnn_num_sequences: int,
    mpnn_sampling_temp: float,
    mpnn_omit_aa: str,
    mpnn_bias_aa: str,
    mpnn_run_parameters: str,
    sequence_design_method: str,
    mpnn_fastrelax_cycles: int,
    refolding_test: str,
    disable_pyrosetta_scoring: bool = True,
) -> dict:
    if sequence_design_method == "fastrelax":
        mpnn_args = f'-omit_AAs "{mpnn_omit_aa.strip()}" -temperature {mpnn_sampling_temp}'
        if mpnn_bias_aa.strip():
            mpnn_args += f' -bias_AA "{mpnn_bias_aa.strip()}"'
    else:
        mpnn_args = f'--omit_AA "{mpnn_omit_aa.strip()}" --temperature {mpnn_sampling_temp}'
        if mpnn_bias_aa.strip():
            mpnn_args += f' --bias_AA "{mpnn_bias_aa.strip()}"'
    if mpnn_run_parameters.strip():
        mpnn_args += f" {mpnn_run_parameters.strip()}"
    return {
        "pipeline": "mn-molecular-design-studio-rfdiffusion",
        "design_type": "binder",
        "rfdiffusion_input_pdb": target_pdb_path,
        "rfdiffusion_num_designs": num_designs,
        "rfdiffusion_contig": contig,
        "rfdiffusion_guidance_preset": rfdiffusion_guidance_preset,
        "rfdiffusion_run_parameters": final_run_parameters,
        "scaffoldguided": scaffoldguided,
        "scaffold_dir": scaffold_dir,
        "scaffold_target_pdb": scaffold_target_pdb,
        "scaffold_target_ss": scaffold_target_ss.strip(),
        "scaffold_target_adj": scaffold_target_adj.strip(),
        "hotspot": hotspots,
        "backbone_filters": backbone_filters.strip() or "none",
        "backbone_hotspot_prefilter": {
            "enabled": apply_backbone_hotspot_prefilter,
            "min_hotspot_contact_fraction": backbone_min_hotspot_contact_fraction,
            "max_min_binder_to_hotspot_distance": backbone_max_hotspot_distance,
            "contact_cutoff": 8.0,
        },
        "mpnn_num_sequences": 1 if sequence_design_method == "fastrelax" else mpnn_num_sequences,
        "mpnn_fastrelax_cycles": mpnn_fastrelax_cycles if sequence_design_method == "fastrelax" else 0,
        "mpnn_run_parameters": mpnn_args.strip(),
        "disable_pyrosetta_scoring": disable_pyrosetta_scoring,
        "refolding_tests": refolding_test or None,
        "planned_steps": [
            "rfdiffusion_backbone",
            "backbone_metrics",
            "sequence_design",
            "refolding_validation" if refolding_test else "skip_refolding_validation",
        ],
    }


def _hydra_overrides_for_bash(overrides: str) -> str:
    return overrides.replace('"', '\\"')


def run_rfdiffusion_classic(
    target_pdb: Path,
    target_chains: list[str],
    contig: str,
    binder_length: str,
    hotspots: str = "",
    campaign_name: str = "",
    num_designs: int = 1,
    timesteps: int = 50,
    model_weights: str = "Complex_base",
    rfdiffusion_guidance_preset: str = "none",
    extra_run_parameters: str = "",
    partial_diffusion: bool = False,
    contigmap_length: str = "",
    inpaint_seq: str = "",
    backbone_filters: str = "",
    apply_backbone_hotspot_prefilter: bool = False,
    backbone_min_hotspot_contact_fraction: float = 0.25,
    backbone_max_hotspot_distance: float = 10.0,
    deterministic: bool = False,
    noise_scale_ca: str = "",
    noise_scale_frame: str = "",
    scaffoldguided: bool = False,
    scaffold_dir: str = RFDIFFUSION_SCAFFOLD_LIBRARY_CONTAINER_DIR,
    scaffold_target_pdb: bool = True,
    scaffold_target_path: str = "",
    scaffold_target_ss: str = "",
    scaffold_target_adj: str = "",
    save_trajectory: bool = False,
    editable_run_parameters: str = "",
    edited_target_pdb_text: str | None = None,
    mpnn_num_sequences: int = 1,
    mpnn_sampling_temp: float = 0.0001,
    mpnn_omit_aa: str = "CX",
    mpnn_bias_aa: str = "",
    mpnn_run_parameters: str = "",
    sequence_design_method: str = "protein_mpnn",
    mpnn_fastrelax_cycles: int = 0,
    refolding_test: str = "af2_model_1_multimer_tt_3rec",
    execution_backend: str = "docker",
    run_vanilla_pipeline: bool = True,
    monomer_refolding_tool: str = "boltz2_monomer",
    complex_refolding_tool: str = "af2_initial_guess",
    complex_template_mode: str = "target_template",
    complex_multimer: bool = True,
    complex_num_recycles: int = 3,
    analysis_keep_top_n: int = 100,
    analysis_keep_per_attempt: int | None = None,
    analysis_rank_metric: str = "analysis_score",
    analysis_thresholds: dict[str, float] | None = None,
    gpu_device: object = "0",
) -> Path:
    sequence_design_method = {
        "ligandmpnn": "protein_mpnn",
        "proteinmpnn": "protein_mpnn",
        "protein_mpnn": "protein_mpnn",
        "ligand_mpnn": "ligand_mpnn",
        "soluble_mpnn": "soluble_mpnn",
        "solublempnn": "soluble_mpnn",
    }.get(str(sequence_design_method), "protein_mpnn")
    hotspots = normalize_hotspots(hotspots) if hotspots else ""
    validate_rfdiffusion_inputs(contig, hotspots, binder_length, num_designs, timesteps)
    validate_rfdiffusion_optional_params(contigmap_length, inpaint_seq)
    manifest = load_manifest("rfdiffusion_classic")
    final_run_parameters = editable_run_parameters.strip() or build_rfdiffusion_run_parameters(
        timesteps=timesteps,
        partial_diffusion=partial_diffusion,
        contigmap_length=contigmap_length,
        inpaint_seq=inpaint_seq,
        model_weights=model_weights,
        deterministic=deterministic,
        noise_scale_ca=noise_scale_ca,
        noise_scale_frame=noise_scale_frame,
        scaffoldguided=scaffoldguided,
        scaffold_dir=scaffold_dir,
        scaffold_target_pdb=scaffold_target_pdb,
        scaffold_target_path=scaffold_target_path,
        scaffold_target_ss=scaffold_target_ss,
        scaffold_target_adj=scaffold_target_adj,
        extra_run_parameters=extra_run_parameters,
    )
    scaffold_selection: dict[str, object] = {}
    if scaffoldguided:
        scaffold_status = prepare_rfdiffusion_scaffold_library()
        if not scaffold_status.get("ready"):
            raise RuntimeError(
                "RFdiffusion scaffold-guided mode was requested, but the scaffold library is not ready "
                f"at {scaffold_status.get('path')}."
            )
        if "scaffoldguided.scaffoldguided" not in final_run_parameters:
            final_run_parameters = " ".join(
                [
                    final_run_parameters,
                    build_rfdiffusion_scaffoldguided_parameters(
                        True,
                        scaffold_dir=scaffold_dir,
                        target_pdb=scaffold_target_pdb,
                        target_path=scaffold_target_path,
                        target_ss=scaffold_target_ss,
                        target_adj=scaffold_target_adj,
                    ),
                ]
            ).strip()
        if scaffold_target_pdb and not re.search(r"(?:^|\s)\+{0,2}scaffoldguided\.target_path=", final_run_parameters):
            final_run_parameters = " ".join(
                [final_run_parameters, "++scaffoldguided.target_path=/work/artifacts/input/target.pdb"]
            ).strip()
        if _has_rfdiffusion_scaffold_list_override(final_run_parameters):
            scaffold_selection = {
                "mode": "manual_or_editable_scaffold_list",
                "scaffold_dir": scaffold_dir.strip() or RFDIFFUSION_SCAFFOLD_LIBRARY_CONTAINER_DIR,
            }
        else:
            selected_scaffold_id = random_rfdiffusion_scaffold_id(Path(scaffold_status["path"]))
            final_run_parameters = " ".join(
                [final_run_parameters, f"++scaffoldguided.scaffold_list=[{selected_scaffold_id}]"]
            ).strip()
            scaffold_selection = {
                "mode": "random_single_scaffold",
                "scaffold_id": selected_scaffold_id,
                "scaffold_dir": scaffold_dir.strip() or RFDIFFUSION_SCAFFOLD_LIBRARY_CONTAINER_DIR,
                "host_scaffold_dir": scaffold_status.get("path"),
                "ss_file": f"{selected_scaffold_id}_ss.pt",
                "adj_file": f"{selected_scaffold_id}_adj.pt",
            }
    params = {
        "target_chains": target_chains,
        "campaign_name": campaign_name.strip(),
        "contig": contig,
        "binder_length": binder_length,
        "hotspots": hotspots,
        "num_designs": num_designs,
        "timesteps": timesteps,
        "model_weights": model_weights,
        "rfdiffusion_guidance_preset": rfdiffusion_guidance_preset,
        "extra_run_parameters": extra_run_parameters.strip(),
        "partial_diffusion": partial_diffusion,
        "contigmap_length": contigmap_length.strip(),
        "inpaint_seq": inpaint_seq.strip(),
        "backbone_filters": backbone_filters.strip(),
        "apply_backbone_hotspot_prefilter": apply_backbone_hotspot_prefilter,
        "backbone_min_hotspot_contact_fraction": backbone_min_hotspot_contact_fraction,
        "backbone_max_hotspot_distance": backbone_max_hotspot_distance,
        "hotspot_contact_cutoff": 8.0,
        "deterministic": deterministic,
        "noise_scale_ca": noise_scale_ca.strip(),
        "noise_scale_frame": noise_scale_frame.strip(),
        "scaffoldguided": scaffoldguided,
        "scaffold_dir": scaffold_dir.strip(),
        "scaffold_target_pdb": scaffold_target_pdb,
        "scaffold_target_path": scaffold_target_path.strip(),
        "scaffold_target_ss": scaffold_target_ss.strip(),
        "scaffold_target_adj": scaffold_target_adj.strip(),
        "scaffold_selection": scaffold_selection,
        "scaffold_selection_mode": scaffold_selection.get("mode", ""),
        "scaffold_id": scaffold_selection.get("scaffold_id", ""),
        "save_trajectory": save_trajectory,
        "rfdiffusion_run_parameters": final_run_parameters,
        "mpnn_num_sequences": mpnn_num_sequences,
        "mpnn_sampling_temp": mpnn_sampling_temp,
        "mpnn_omit_aa": mpnn_omit_aa.strip(),
        "mpnn_bias_aa": mpnn_bias_aa.strip(),
        "mpnn_run_parameters": mpnn_run_parameters.strip(),
        "sequence_design_method": sequence_design_method,
        "mpnn_fastrelax_cycles": mpnn_fastrelax_cycles,
        "refolding_test": refolding_test,
        "execution_backend": execution_backend,
        "run_vanilla_pipeline": run_vanilla_pipeline,
        "monomer_refolding_tool": monomer_refolding_tool,
        "complex_refolding_tool": complex_refolding_tool,
        "complex_template_mode": complex_template_mode,
        "complex_multimer": complex_multimer,
        "complex_num_recycles": complex_num_recycles,
        "analysis_keep_top_n": analysis_keep_top_n,
        "analysis_keep_per_attempt": analysis_keep_per_attempt,
        "analysis_rank_metric": analysis_rank_metric,
        "analysis_thresholds": analysis_thresholds or {},
        "gpu_device": normalize_gpu_device(gpu_device),
        "pipeline_mode": "mn_molecular_design_studio.engine.workflow_engine_rfdiffusion",
    }
    job = create_job(
        DESIGN_GROUP,
        job_type="design_campaign",
        tool="rfdiffusion_classic",
        inputs={"target_pdb": str(target_pdb), "target_chains": target_chains},
        params=params,
    )
    if campaign_name.strip():
        update_status(job.run_dir, "queued", campaign_name=campaign_name.strip())

    target_artifact = artifact_path(job.run_dir, "input", "target.pdb")
    target_text = edited_target_pdb_text if edited_target_pdb_text is not None else target_pdb.read_text(errors="ignore")
    target_artifact.write_text(filter_pdb_text(target_text, keep_chains=set(target_chains) or None))

    raw_dir = job.run_dir / "artifacts" / "raw" / "rfdiffusion"
    raw_dir.mkdir(parents=True, exist_ok=True)
    pipeline_dir = artifact_path(job.run_dir, "raw", "rfdiffusion", "pipeline_inputs")
    pipeline_dir.mkdir(parents=True, exist_ok=True)
    pipeline_payload = _rfdiffusion_pipeline_payload(
        target_pdb_path="/work/artifacts/input/target.pdb",
        contig=contig,
        hotspots=hotspots,
        num_designs=num_designs,
        final_run_parameters=final_run_parameters,
        rfdiffusion_guidance_preset=rfdiffusion_guidance_preset,
        scaffoldguided=scaffoldguided,
        scaffold_dir=scaffold_dir,
        scaffold_target_pdb=scaffold_target_pdb,
        scaffold_target_ss=scaffold_target_ss,
        scaffold_target_adj=scaffold_target_adj,
        backbone_filters=backbone_filters,
        apply_backbone_hotspot_prefilter=apply_backbone_hotspot_prefilter,
        backbone_min_hotspot_contact_fraction=backbone_min_hotspot_contact_fraction,
        backbone_max_hotspot_distance=backbone_max_hotspot_distance,
        mpnn_num_sequences=mpnn_num_sequences,
        mpnn_sampling_temp=mpnn_sampling_temp,
        mpnn_omit_aa=mpnn_omit_aa,
        mpnn_bias_aa=mpnn_bias_aa,
        mpnn_run_parameters=mpnn_run_parameters,
        sequence_design_method=sequence_design_method,
        mpnn_fastrelax_cycles=mpnn_fastrelax_cycles,
        refolding_test=refolding_test,
    )
    write_json(pipeline_dir / "mn_rfdiffusion_pipeline_params.json", pipeline_payload)
    if scaffold_selection:
        pipeline_payload["scaffold_selection"] = scaffold_selection
        write_json(pipeline_dir / "mn_rfdiffusion_pipeline_params.json", pipeline_payload)
        write_json(pipeline_dir / "rfdiffusion_scaffold_selection.json", scaffold_selection)
    (pipeline_dir / "rfdiffusion_run_parameters.txt").write_text(final_run_parameters + "\n")
    (pipeline_dir / "rfdiffusion_contig.txt").write_text(contig + "\n")
    (pipeline_dir / "hotspots.txt").write_text(hotspots + "\n")

    hotspot_arg = f' "ppi.hotspot_res=[{hotspots}]"' if hotspots else ""
    traj_arg = "true" if save_trajectory else "false"
    shell_run_parameters = _hydra_overrides_for_bash(final_run_parameters)
    rfdiffusion_models_dir = reference_root() / "rfdiffusion_models"
    script = (
        "set -euxo pipefail; "
        "rm -rf /work/artifacts/raw/rfdiffusion/output; "
        "export PYTHONPATH=/opt/RFdiffusion:/opt/RFdiffusion/env/SE3Transformer; "
        "mkdir -p /work/artifacts/raw/rfdiffusion/output; "
        "cd /work/artifacts/raw/rfdiffusion; "
        "python3 /opt/RFdiffusion/scripts/run_inference.py "
        "inference.output_prefix=output/design "
        "inference.model_directory_path=/models "
        "inference.schedule_directory_path=/opt/RFdiffusion/schedules "
        "inference.input_pdb=/work/artifacts/input/target.pdb "
        f"inference.num_designs={num_designs} "
        f'"contigmap.contigs=[{contig}]"'
        f"{hotspot_arg} "
        f"inference.write_trajectory={traj_arg} "
        f"{shell_run_parameters}; "
        'test -n "$(find /work/artifacts/raw/rfdiffusion/output -maxdepth 1 -type f -name \"design_*.pdb\" -print -quit)"'
    )
    if execution_backend == "nextflow":
        if shutil.which("nextflow") is None:
            raise RuntimeError("Nextflow is not on PATH. Update the app env with `conda env update -f environment.yml --prune`.")
        nextflow_dir = raw_dir / "nextflow"
        nextflow_dir.mkdir(parents=True, exist_ok=True)
        pipeline_path = Path(__file__).resolve().parents[1] / "pipelines" / "rfdiffusion_backbone"
        command = [
            "bash",
            "-lc",
            " ".join(
                [
                    f"cd {shlex.quote(str(nextflow_dir))}",
                    "&&",
                    "nextflow",
                    "run",
                    shlex.quote(str(pipeline_path)),
                    "-config",
                    shlex.quote(str(pipeline_path / "nextflow.config")),
                    "-work-dir",
                    shlex.quote(str(nextflow_dir / "work")),
                    "-with-trace",
                    "trace.txt",
                    "-with-report",
                    "report.html",
                    "-with-timeline",
                    "timeline.html",
                    "--input_pdb",
                    shlex.quote(str(target_artifact.resolve())),
                    "--publish_dir",
                    shlex.quote(str(raw_dir / "output")),
                    "--rfdiffusion_image",
                    shlex.quote(str(manifest["image"])),
                    "--rfdiffusion_num_designs",
                    str(num_designs),
                    "--rfdiffusion_contig",
                    shlex.quote(contig),
                    "--hotspot",
                    shlex.quote(hotspots),
                    "--write_trajectory",
                    str(save_trajectory).lower(),
                    "--rfdiffusion_run_parameters",
                    shlex.quote(final_run_parameters),
                    "--rfdiffusion_models_dir",
                    shlex.quote(str(rfdiffusion_models_dir)),
                ]
            ),
        ]
        steps = [{"name": "rfdiffusion-nextflow-backbone", "command": command}]
    else:
        command = [
            "docker",
            "run",
            "--rm",
            *docker_gpu_args(gpu_device),
            "-v",
            f"{job.run_dir}:/work",
            "-v",
            f"{rfdiffusion_models_dir}:/models:ro",
            "-w",
            "/work",
            manifest["image"],
            "bash",
            "-lc",
            script,
        ]
        steps = [{"name": "rfdiffusion-docker-backbone", "command": command}]
    rc = _run_shell_steps(job.run_dir, steps)

    candidates = _normalize_rfdiffusion_candidates(job.run_dir, params, target_artifact) if rc == 0 else []
    child_runs: list[Path] = []
    if rc == 0 and run_vanilla_pipeline and candidates:
        from mn_molecular_design_studio.engine.workflow_engine.workflows.campaigns import run_lineage_steps

        thresholds = analysis_thresholds or {}
        steps = [
            {
                "module": "sequence_design",
                "tool": "ligandmpnn",
                "params": {
                    "model_type": sequence_design_method if sequence_design_method in {"protein_mpnn", "ligand_mpnn", "soluble_mpnn"} else "protein_mpnn",
                    "design_chains": "",
                    "num_seq_per_target": int(mpnn_num_sequences),
                    "sampling_temp": float(mpnn_sampling_temp),
                    "omit_aas": str(mpnn_omit_aa or "CX"),
                    "seed": None,
                    "require_backbone_hotspot_filter_pass": bool(apply_backbone_hotspot_prefilter),
                    "gpu_device": normalize_gpu_device(gpu_device),
                },
            },
            {
                "module": "monomer_refolding",
                "tool": monomer_refolding_tool,
                "params": {
                    "min_plddt": float(thresholds.get("min_binder_plddt", 70.0)),
                    "gpu_device": normalize_gpu_device(gpu_device),
                },
            },
            {
                "module": "complex_refolding",
                "tool": complex_refolding_tool,
                "params": {
                    "require_monomer_success": True,
                    "template_mode": complex_template_mode,
                    "multimer": complex_multimer,
                    "num_recycles": complex_num_recycles,
                    "gpu_device": normalize_gpu_device(gpu_device),
                },
            },
            {
                "module": "analysis",
                "tool": "ranking",
                "params": {
                    "keep_top_n": int(analysis_keep_top_n),
                    "keep_per_attempt": analysis_keep_per_attempt,
                    "ranking_metric": analysis_rank_metric,
                    "thresholds": thresholds,
                },
            },
        ]
        child_runs = run_lineage_steps(
            campaign_name.strip() or f"RFdiffusion classic vanilla from {job.run_dir.name}",
            _job_source_for_run(job.run_dir),
            steps,
        )
        for step_index, child_run in enumerate(child_runs, start=1):
            _mark_internal_child_run(child_run, job.run_dir, step_index)
            _append_internal_child_logs(job.run_dir, child_run, f"vanilla-{step_index}")
        if len(child_runs) >= 3 and read_json(child_runs[2] / "result.json").get("success") is True:
            candidates = _adopt_rfdiffusion_native_candidates(job.run_dir, child_runs[2], params, target_artifact)
        if not child_runs or read_json(child_runs[-1] / "result.json").get("success") is not True:
            rc = 1
    artifacts = _collect_artifacts(
        job.run_dir,
        job.run_dir / "artifacts",
        [
            ("raw/rfdiffusion/output/*.pdb", "designed_complex_pdb"),
            ("raw/rfdiffusion/output/*.trb", "rfdiffusion_trb"),
            ("raw/rfdiffusion/native_pipeline/**/*", "rfdiffusion_native_pipeline"),
            ("raw/rfdiffusion/pipeline_inputs/*", "rfdiffusion_pipeline_input"),
            ("raw/rfdiffusion/nextflow/*.txt", "nextflow_trace"),
            ("raw/rfdiffusion/nextflow/*.html", "nextflow_report"),
            ("normalized_candidates/*.jsonl", "normalized_candidates"),
            ("normalized_candidates/*.json", "campaign_result"),
        ],
    )
    finish_job(
        job.run_dir,
        rc == 0 and bool(candidates),
        {
            "outputs": {"artifacts": artifacts, "candidates": candidates},
            "metrics": {
                "return_code": rc,
                "artifact_count": len(artifacts),
                "candidate_count": len(candidates),
                "internal_pipeline_steps": len(child_runs),
            },
            "downstream_artifacts": {
                "target_pdb": str(target_artifact.relative_to(job.run_dir)),
                "pipeline_inputs": "artifacts/raw/rfdiffusion/pipeline_inputs",
                "candidates_jsonl": "artifacts/normalized_candidates/candidates.jsonl",
                "campaign_result": "artifacts/normalized_candidates/campaign_result.json",
            },
        },
    )
    return job.run_dir


def run_bindcraft(
    target_pdb: Path,
    target_chains: list[str],
    binder_length: str,
    hotspots: str = "",
    campaign_name: str = "",
    number_of_final_designs: int = 1,
    num_seqs_override: int | None = None,
    max_mpnn_sequences_override: int | None = None,
    time_limit_seconds: int | None = None,
    enable_mpnn: bool = True,
    max_trajectories: int = 1,
    filter_settings: str = "default_filters.json",
    advanced_settings_file: str = "default_4stage_multimer_mpnn.json",
    gpu_device: object = "0",
) -> Path:
    if not target_chains:
        raise ValueError("At least one target chain is required.")
    hotspots = normalize_hotspots(hotspots) if hotspots else ""
    parse_binder_lengths(binder_length)
    if number_of_final_designs < 1:
        raise ValueError("Desired final designs must be at least 1.")
    if num_seqs_override is not None and num_seqs_override < 1:
        raise ValueError("num_seqs must be at least 1.")
    if max_mpnn_sequences_override is not None and max_mpnn_sequences_override < 1:
        raise ValueError("max_mpnn_sequences must be at least 1.")
    if time_limit_seconds is not None and time_limit_seconds < 60:
        raise ValueError("Time limit should be at least 60 seconds.")
    if max_trajectories < 1:
        raise ValueError("Max trajectories must be at least 1.")
    advanced_settings_path = _bindcraft_resource("settings_advanced", advanced_settings_file)
    if not advanced_settings_path.exists():
        raise ValueError(f"Unknown BindCraft advanced settings file: {advanced_settings_file}")
    manifest = load_manifest("bindcraft")
    params = {
        "target_chains": target_chains,
        "binder_length": binder_length,
        "hotspots": hotspots,
        "campaign_name": campaign_name.strip(),
        "number_of_final_designs": number_of_final_designs,
        "num_seqs_override": num_seqs_override,
        "max_mpnn_sequences_override": max_mpnn_sequences_override,
        "time_limit_seconds": time_limit_seconds,
        "enable_mpnn": enable_mpnn,
        "max_trajectories": max_trajectories,
        "filter_settings": filter_settings,
        "advanced_settings_file": advanced_settings_file,
        "gpu_device": normalize_gpu_device(gpu_device),
    }
    job = create_job(
        DESIGN_GROUP,
        job_type="design_campaign",
        tool="bindcraft",
        inputs={"target_pdb": str(target_pdb), "target_chains": target_chains},
        params=params,
    )
    target_artifact = artifact_path(job.run_dir, "raw", "bindcraft", "input_pdb.pdb")
    target_text = target_pdb.read_text(errors="ignore")
    target_artifact.write_text(filter_pdb_text(target_text, keep_chains=set(target_chains) or None))
    _copy_bindcraft_settings(job.run_dir, params)

    bindcraft_command = "python /content/bindcraft/bindcraft.py"
    if time_limit_seconds is not None:
        bindcraft_command = f"timeout {int(time_limit_seconds)} {bindcraft_command}"
    script = (
        "set -euxo pipefail; "
        "cd /work/artifacts/raw/bindcraft; "
        "ln -sf input_pdb.pdb target.pdb; "
        "ln -sfn /af_models alphafold_models_path; "
        "cp settings_advanced.json advanced.json; "
        f"{bindcraft_command} "
        "--settings input.json --advanced advanced.json --filters settings_filters.json; "
        'test -n "$(find output -type f -name \"*.pdb\" -print -quit)"'
    )
    command = [
        "docker",
        "run",
        "--rm",
        *docker_gpu_args(gpu_device),
        "--shm-size=64G",
        "-v",
        f"{job.run_dir}:/work",
        "-v",
        "/mnt/db/reference_files/alphafold_models:/af_models:ro",
        "-w",
        "/work",
        manifest["image"],
        "bash",
        "-lc",
        script,
    ]
    rc = _run_shell_steps(job.run_dir, [{"name": "bindcraft", "command": command}])
    candidates = _normalize_bindcraft_candidates(job.run_dir, params, target_artifact)
    artifacts = _collect_artifacts(
        job.run_dir,
        job.run_dir / "artifacts",
        [
            ("raw/bindcraft/output/**/*.pdb", "designed_complex_pdb"),
            ("raw/bindcraft/output/**/*.csv", "bindcraft_metrics_table"),
            ("normalized_candidates/*.jsonl", "normalized_candidates"),
            ("normalized_candidates/*.json", "campaign_result"),
            ("raw/bindcraft/*.json", "bindcraft_input_settings"),
        ],
    )
    finish_job(
        job.run_dir,
        bool(candidates),
        {
            "outputs": {"artifacts": artifacts, "candidates": candidates},
            "metrics": {"return_code": rc, "artifact_count": len(artifacts), "candidate_count": len(candidates)},
            "downstream_artifacts": {
                "target_pdb": str(target_artifact.relative_to(job.run_dir)),
                "candidates_jsonl": "artifacts/normalized_candidates/candidates.jsonl",
                "campaign_result": "artifacts/normalized_candidates/campaign_result.json",
            },
        },
    )
    return job.run_dir


def run_rfdiffusion3_foundry(
    target_pdb: Path,
    target_chains: list[str],
    binder_length: str,
    hotspots: str = "",
    campaign_name: str = "",
    num_designs: int = 1,
    timesteps: int = 50,
    contig: str | None = None,
    is_non_loopy: bool = True,
    infer_ori_strategy: str = "",
    run_vanilla_pipeline: bool = True,
    mpnn_sequences_per_backbone: int = 1,
    mpnn_model_type: str = "protein_mpnn",
    mpnn_checkpoint_path: str = "/weights/proteinmpnn_v_48_020.pt",
    rf3_checkpoint_path: str = "/weights/rf3_foundry_01_24_latest_remapped.ckpt",
    prepare_target_msa: bool = True,
    gpu_device: object = "0",
) -> Path:
    hotspots = normalize_hotspots(hotspots) if hotspots else ""
    manifest = load_manifest("rfdiffusion3_foundry")
    final_contig = _normalize_rfdiffusion3_contig(contig) if contig and contig.strip() else ""
    if mpnn_sequences_per_backbone < 1:
        raise ValueError("MPNN sequences per backbone must be at least 1.")
    params = {
        "target_chains": target_chains,
        "binder_length": binder_length,
        "hotspots": hotspots,
        "campaign_name": campaign_name.strip(),
        "num_designs": num_designs,
        "timesteps": timesteps,
        "contig": final_contig,
        "is_non_loopy": is_non_loopy,
        "infer_ori_strategy": infer_ori_strategy.strip() or ("hotspots" if hotspots else "default"),
        "run_vanilla_pipeline": run_vanilla_pipeline,
        "mpnn_sequences_per_backbone": mpnn_sequences_per_backbone,
        "mpnn_model_type": mpnn_model_type,
        "mpnn_checkpoint_path": mpnn_checkpoint_path,
        "rf3_checkpoint_path": rf3_checkpoint_path,
        "prepare_target_msa": prepare_target_msa,
        "gpu_device": normalize_gpu_device(gpu_device),
    }
    job = create_job(DESIGN_GROUP, "design_campaign", "rfdiffusion3_foundry", {"target_pdb": str(target_pdb), "target_chains": target_chains}, params)
    if campaign_name.strip():
        update_status(job.run_dir, "queued", campaign_name=campaign_name.strip())
    target_artifact = _copy_target_for_design(job.run_dir, target_pdb, target_chains, "raw/rfdiffusion3_foundry")
    if not final_contig:
        final_contig = _normalize_rfdiffusion3_contig(default_target_contig(target_artifact, target_chains, binder_length))
        params["contig"] = final_contig
        write_json(
            job.run_dir / "input.json",
            {
                "job_type": "design_campaign",
                "tool": "rfdiffusion3_foundry",
                "inputs": {"target_pdb": str(target_pdb), "target_chains": target_chains},
                "params": params,
            },
        )
    msa_by_source_chain: dict[str, str] = {}
    if run_vanilla_pipeline and prepare_target_msa:
        try:
            msa_by_source_chain = _shared_ensure_boltz_msas_for_target(
                job.run_dir,
                target_artifact,
                target_chains,
                raw_subdir="rfdiffusion3_foundry/rf3_msa_cache",
                gpu_device=gpu_device,
            )
            params["target_msa_by_chain"] = msa_by_source_chain
            write_json(
                job.run_dir / "input.json",
                {
                    "job_type": "design_campaign",
                    "tool": "rfdiffusion3_foundry",
                    "inputs": {"target_pdb": str(target_pdb), "target_chains": target_chains},
                    "params": params,
                },
            )
        except Exception as exc:
            params["target_msa_error"] = str(exc)
            write_json(
                job.run_dir / "input.json",
                {
                    "job_type": "design_campaign",
                    "tool": "rfdiffusion3_foundry",
                    "inputs": {"target_pdb": str(target_pdb), "target_chains": target_chains},
                    "params": params,
                },
            )
    design_input = {
        "design_1": {
            "dialect": 2,
            "input": "target.pdb",
            "contig": final_contig,
            "is_non_loopy": is_non_loopy,
            "infer_ori_strategy": params["infer_ori_strategy"],
        }
    }
    if hotspots:
        design_input["design_1"]["select_hotspots"] = {token: "CA" for token in hotspots.split(",") if token}
    write_json(job.run_dir / "artifacts" / "raw" / "rfdiffusion3_foundry" / "rfd3_inputs_staged.json", design_input)
    if run_vanilla_pipeline:
        script = (
            "set -euxo pipefail; cd /work/artifacts/raw/rfdiffusion3_foundry; "
            "rm -rf rfd3 mpnn rf3 tmp foundry_native_mapping.tsv; mkdir -p rfd3 mpnn rf3 tmp; "
            f"rfd3 design out_dir=rfd3 inputs=rfd3_inputs_staged.json ckpt_path=/weights/rfd3_latest.ckpt "
            f"diffusion_batch_size=1 n_batches={num_designs} inference_sampler.num_timesteps={timesteps} "
            "skip_existing=False prevalidate_inputs=True; "
            'test -n "$(find rfd3 -maxdepth 1 -type f \\( -name \"*.cif\" -o -name \"*.cif.gz\" \\) -print -quit)"'
        )
    else:
        script = (
            "set -euxo pipefail; cd /work/artifacts/raw/rfdiffusion3_foundry; rm -rf rfd3; "
            f"rfd3 design out_dir=rfd3 inputs=rfd3_inputs_staged.json ckpt_path=/weights/rfd3_latest.ckpt "
            f"diffusion_batch_size=1 n_batches={num_designs} inference_sampler.num_timesteps={timesteps} "
            "skip_existing=False prevalidate_inputs=True; "
            'test -n "$(find rfd3 -type f \\( -name \"*.cif\" -o -name \"*.cif.gz\" \\) -print -quit)"'
        )
    command = [
        "docker",
        "run",
        "--rm",
        *docker_gpu_args(gpu_device),
        "-v",
        f"{job.run_dir}:/work",
        "-v",
        "/mnt/db/reference_files/foundry:/weights:ro",
        "-v",
        "/mnt/db/reference_files/boltz_models/msa_repository:/msa_repository:ro",
        "-w",
        "/work",
        manifest["image"],
        "bash",
        "-lc",
        script,
    ]
    rc = _run_shell_steps(job.run_dir, [{"name": "rfdiffusion3-foundry", "command": command}])
    if rc == 0 and run_vanilla_pipeline:
        raw_root = job.run_dir / "artifacts" / "raw" / "rfdiffusion3_foundry"
        mapping_path = raw_root / "foundry_native_mapping.tsv"
        with mapping_path.open("w", newline="") as handle:
            fieldnames = ["rfd3_cif", "rfd3_id", "mpnn_cif", "mpnn_id", "rf3_out_dir", "target_msa_by_rf3_chain"]
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            for rfd3_index, rfd3_path in enumerate(sorted([*raw_root.glob("rfd3/*.cif"), *raw_root.glob("rfd3/*.cif.gz")]), start=1):
                rfd3_id = rfd3_path.name.removesuffix(".cif.gz").removesuffix(".cif")
                staged_cif = raw_root / "tmp" / f"{rfd3_id}.cif"
                if rfd3_path.name.endswith(".gz"):
                    with gzip.open(rfd3_path, "rt", errors="ignore") as src:
                        staged_cif.write_text(src.read())
                else:
                    shutil.copy2(rfd3_path, staged_cif)
                mpnn_binder_chains, _mpnn_target_chains, _mpnn_role_inference = _infer_generated_chain_roles(
                    target_artifact,
                    staged_cif,
                    target_chains,
                )
                if not mpnn_binder_chains:
                    with (job.run_dir / "stderr.log").open("a") as stderr:
                        stderr.write(f"Could not identify binder chain(s) for Foundry MPNN input {rfd3_id}; target chain will not be redesigned.\n")
                    rc = 1
                    break
                mpnn_out_dir = raw_root / "mpnn" / rfd3_id
                mpnn_out_dir.mkdir(parents=True, exist_ok=True)
                mpnn_command = [
                    "docker",
                    "run",
                    "--rm",
                    *docker_gpu_args(gpu_device),
                    "-v",
                    f"{job.run_dir}:/work",
                    "-v",
                    "/mnt/db/reference_files/foundry:/weights:ro",
                    "-w",
                    "/work",
                    manifest["image"],
                    "mpnn",
                    "--structure_path",
                    f"/work/artifacts/raw/rfdiffusion3_foundry/tmp/{rfd3_id}.cif",
                    "--checkpoint_path",
                    mpnn_checkpoint_path,
                    "--is_legacy_weights",
                    "True",
                    "--model_type",
                    mpnn_model_type,
                    "--designed_chains",
                    ",".join(mpnn_binder_chains),
                    "--batch_size",
                    str(mpnn_sequences_per_backbone),
                    "--number_of_batches",
                    "1",
                    "--remove_waters",
                    "True",
                    "--out_directory",
                    f"/work/artifacts/raw/rfdiffusion3_foundry/mpnn/{rfd3_id}",
                ]
                rc = _run_shell_steps(job.run_dir, [{"name": f"foundry-mpnn-{rfd3_id}", "command": mpnn_command}])
                if rc != 0:
                    break
                for mpnn_path in sorted(mpnn_out_dir.glob("*.cif")):
                    mpnn_id = mpnn_path.stem
                    binder_chains, output_target_chains, _inference = _infer_generated_chain_roles(
                        target_artifact,
                        mpnn_path,
                        target_chains,
                    )
                    rf3_input = raw_root / "rf3_inputs" / rfd3_id / f"{mpnn_id}.cif"
                    target_msa_by_rf3_chain = _write_foundry_cif_with_target_msas(
                        mpnn_path,
                        rf3_input,
                        target_artifact,
                        output_target_chains,
                        target_chains,
                        msa_by_source_chain,
                    )
                    rf3_out_dir = raw_root / "rf3" / rfd3_id / mpnn_id
                    rf3_command = [
                        "docker",
                        "run",
                        "--rm",
                        *docker_gpu_args(gpu_device),
                        "-v",
                        f"{job.run_dir}:/work",
                        "-v",
                        "/mnt/db/reference_files/foundry:/weights:ro",
                        "-v",
                        "/mnt/db/reference_files/boltz_models/msa_repository:/msa_repository:ro",
                        "-w",
                        "/work",
                        manifest["image"],
                        "rf3",
                        "fold",
                        f"inputs=/work/artifacts/raw/rfdiffusion3_foundry/rf3_inputs/{rfd3_id}/{mpnn_id}.cif",
                        f"ckpt_path={rf3_checkpoint_path}",
                        f"out_dir=/work/artifacts/raw/rfdiffusion3_foundry/rf3/{rfd3_id}/{mpnn_id}",
                        "raise_if_missing_msa_for_protein_of_length_n=10000",
                    ]
                    rc = _run_shell_steps(job.run_dir, [{"name": f"rf3-fold-{rfd3_id}-{mpnn_id}", "command": rf3_command}])
                    if rc != 0:
                        break
                    writer.writerow(
                        {
                            "rfd3_cif": _rel_path(job.run_dir, rfd3_path),
                            "rfd3_id": rfd3_id,
                            "mpnn_cif": _rel_path(job.run_dir, mpnn_path),
                            "mpnn_id": mpnn_id,
                            "rf3_out_dir": _rel_path(job.run_dir, rf3_out_dir),
                            "target_msa_by_rf3_chain": json.dumps(target_msa_by_rf3_chain, sort_keys=True),
                        }
                    )
                if rc != 0:
                    break
        candidates = _normalize_rfdiffusion3_foundry_native_candidates(job.run_dir, params, target_artifact) if rc == 0 else []
    else:
        candidates = _normalize_rfdiffusion3_candidates(job.run_dir, params, target_artifact) if rc == 0 else []
    _finish_design_job(
        job.run_dir,
        rc == 0,
        rc,
        candidates,
        [
            ("raw/rfdiffusion3_foundry/rfd3/**/*", "rfdiffusion3_output"),
            ("raw/rfdiffusion3_foundry/mpnn/**/*", "foundry_mpnn_output"),
            ("raw/rfdiffusion3_foundry/rf3/**/*", "rf3_output"),
            ("raw/rfdiffusion3_foundry/rf3_inputs/**/*", "rf3_input"),
            ("raw/rfdiffusion3_foundry/*.json", "rfdiffusion3_input"),
            ("raw/rfdiffusion3_foundry/*.tsv", "foundry_native_mapping"),
        ],
    )
    return job.run_dir


def run_boltzgen(
    target_pdb: Path | None = None,
    target_chains: list[str] | None = None,
    binder_length: str = "60-180",
    hotspots: str = "",
    non_binding_residues: str = "",
    campaign_name: str = "",
    num_designs: int = 1,
    budget: int = 1,
    sampling_steps: int = 20,
    diffusion_batch_size: int | None = None,
    omit_amino_acids: str = "",
    additional_filters: list[str] | None = None,
    filter_binding_site_adherence: bool = False,
    run_vanilla_pipeline: bool = True,
    generator_only_tile: bool = False,
    gpu_device: object = "0",
    task: str = "protein",
    scaffold_name: str = "7eow",
    ligand_type: str = "ccd",
    ligand_value: str = "TSA",
    peptide_sequence: str = "12..20",
    cyclotide_sequence: str = "3C8C6C5C3C1C2",
    redesign_residues: str = "",
    design_yaml: str = "",
    yaml_protocol: str = "protein-anything",
    mimosa_condition: str = "unconstrained",
    mimosa_complex: Path | None = None,
    mimosa_binder_chain: str = "",
    mimosa_complex_target_chain: str = "",
    mimosa_target_chain: str = "",
    mimosa_cdr: str = "CDR3",
    mimosa_design_loop: str = "CDR3",
    mimosa_source_positions: str | list[str] = "",
    mimosa_contact_cutoff: float = 8.0,
    mimosa_guidance_scale: float = 0.3,
    mimosa_num_particles: int = 4,
    mimosa_fk_lambda: float = 4.0,
    mimosa_resampling_interval: int = 3,
) -> Path:
    supported_tasks = {"protein", "nanobody", "small_molecule", "peptide", "cyclotide", "redesign", "yaml"}
    if task not in supported_tasks:
        raise ValueError(f"Unsupported BoltzGen task: {task}")
    if generator_only_tile:
        run_vanilla_pipeline = False
        # This setting is consumed only by BoltzGen's inverse-folding stage.
        omit_amino_acids = ""
    target_chains = list(target_chains or [])
    hotspots = normalize_hotspots(hotspots) if hotspots else ""
    non_binding_residues = normalize_hotspots(non_binding_residues) if non_binding_residues else ""
    redesign_residues = normalize_hotspots(redesign_residues) if redesign_residues else ""
    omit_amino_acids = "".join(dict.fromkeys(str(omit_amino_acids or "").upper().strip()))
    additional_filters = [str(value).strip() for value in (additional_filters or []) if str(value).strip()]
    if task != "small_molecule" and (target_pdb is None or not target_chains):
        raise ValueError("BoltzGen requires a target structure and at least one target chain for this task.")
    if task == "small_molecule" and ligand_type not in {"ccd", "smiles"}:
        raise ValueError("BoltzGen ligand type must be CCD or SMILES.")
    if task == "small_molecule" and not ligand_value.strip():
        raise ValueError("Enter a CCD code or SMILES string for the BoltzGen ligand.")
    if task == "nanobody" and scaffold_name not in {"7eow", "7xl0", "8coh", "8z8v"}:
        raise ValueError("Choose one of BoltzGen's bundled nanobody scaffolds.")
    if task == "redesign" and not redesign_residues:
        raise ValueError("Select at least one residue to redesign.")
    if task == "yaml" and not design_yaml.strip():
        raise ValueError("Enter a BoltzGen design YAML specification.")
    native_protocols = {"protein-anything", "peptide-anything", "protein-small_molecule", "nanobody-anything"}
    if task == "yaml" and yaml_protocol not in native_protocols:
        raise ValueError("Choose a native BoltzGen protocol for the YAML configuration.")
    mimosa_condition = str(mimosa_condition).strip().lower()
    mimosa_source_positions = normalize_mimosa_source_positions(mimosa_source_positions)
    if mimosa_condition not in {"geometry_only", "sequence_fixed", "unconstrained"}:
        raise ValueError("MIMOSA condition must be geometry_only, sequence_fixed, or unconstrained.")
    mimosa_enabled = mimosa_condition != "unconstrained"
    if mimosa_enabled:
        if task != "nanobody":
            raise ValueError("BoltzGen MIMOSA currently supports its nanobody VHH scaffold workflow.")
        if generator_only_tile:
            raise ValueError("MIMOSA needs BoltzGen's inverse-folding stage; use the full pipeline.")
        if mimosa_complex is None or not Path(mimosa_complex).is_file():
            raise ValueError("MIMOSA requires a cognate binder–target complex PDB or mmCIF.")
        if Path(mimosa_complex).suffix.lower() not in {".pdb", ".ent", ".cif", ".mmcif"}:
            raise ValueError("MIMOSA source complexes must be PDB or mmCIF files.")
        if not all((mimosa_binder_chain, mimosa_complex_target_chain, mimosa_target_chain)):
            raise ValueError("Select the source partner, source target, and design-target chains for MIMOSA.")
        if mimosa_target_chain not in target_chains:
            raise ValueError("The MIMOSA design-target chain must be among the selected target chains.")
        if mimosa_design_loop not in {"CDR1", "CDR2", "CDR3", "all_cdrs"}:
            raise ValueError("Choose one nanobody CDR or all CDRs for MIMOSA guidance.")
        if not mimosa_source_positions:
            raise ValueError("Select at least one contacting amino acid from the source partner chain.")
        if not 0 < float(mimosa_contact_cutoff) <= 20:
            raise ValueError("MIMOSA contact cutoff must be greater than 0 and at most 20 Å.")
        if not 0 <= float(mimosa_guidance_scale) <= 2:
            raise ValueError("MIMOSA guidance scale must be between 0 and 2.")
        if not 1 <= int(mimosa_num_particles) <= 64:
            raise ValueError("MIMOSA particle count must be between 1 and 64.")
        if float(mimosa_fk_lambda) < 0:
            raise ValueError("MIMOSA Feynman-Kac lambda must be non-negative.")
        if int(mimosa_resampling_interval) < 1:
            raise ValueError("MIMOSA resampling interval must be at least one step.")
    lengths = parse_binder_lengths(binder_length)
    binder_sequence_spec = str(lengths[0]) if len(lengths) == 1 else f"{lengths[0]}..{lengths[1]}"
    manifest = load_manifest("boltzgen")
    if mimosa_enabled:
        manifest["image"] = MIMOSA_BOLTZGEN_IMAGE
        image_check = subprocess.run(["docker", "image", "inspect", MIMOSA_BOLTZGEN_IMAGE], capture_output=True, text=True, check=False)
        if image_check.returncode != 0:
            raise RuntimeError(
                f"MIMOSA runtime image {MIMOSA_BOLTZGEN_IMAGE!r} is unavailable. "
                "Restore or rebuild the Studio runtime image before submitting this job."
            )
    if diffusion_batch_size is not None and diffusion_batch_size < 1:
        raise ValueError("BoltzGen diffusion batch size must be at least 1.")
    invalid_amino_acids = sorted(set(omit_amino_acids) - set("ACDEFGHIKLMNPQRSTVWY"))
    if invalid_amino_acids:
        raise ValueError("BoltzGen omitted amino acids must use standard one-letter amino-acid codes.")
    if len(omit_amino_acids) == 20:
        raise ValueError("BoltzGen must have at least one amino acid available for inverse folding.")
    if additional_filters and not run_vanilla_pipeline:
        raise ValueError("BoltzGen additional filters require the full refinement and ranking pipeline.")
    if filter_binding_site_adherence and not run_vanilla_pipeline:
        raise ValueError("BoltzGen binding-site adherence filtering requires the full pipeline.")
    if filter_binding_site_adherence and not hotspots:
        raise ValueError("Select at least one binding-site residue for BoltzGen's native adherence filter.")
    if budget < 1:
        raise ValueError("BoltzGen final design budget must be at least 1.")
    if budget > num_designs:
        budget = num_designs
    protocol = {
        "protein": "protein-anything",
        "nanobody": "nanobody-anything",
        "small_molecule": "protein-small_molecule",
        "peptide": "peptide-anything",
        "cyclotide": "peptide-anything",
        "redesign": "protein-anything",
        "yaml": yaml_protocol,
    }[task]
    params = {
        "task": task,
        "protocol": protocol,
        "target_chains": target_chains,
        "binder_length": binder_length,
        "hotspots": hotspots,
        "non_binding_residues": non_binding_residues,
        "campaign_name": campaign_name.strip(),
        "num_designs": num_designs,
        "budget": budget,
        "sampling_steps": sampling_steps,
        "diffusion_batch_size": diffusion_batch_size,
        "omit_amino_acids": omit_amino_acids,
        "additional_filters": additional_filters,
        "filter_binding_site_adherence": filter_binding_site_adherence,
        "run_vanilla_pipeline": run_vanilla_pipeline,
        "generator_only_tile": bool(generator_only_tile),
        "gpu_device": normalize_gpu_device(gpu_device),
        "scaffold_name": scaffold_name,
        "ligand_type": ligand_type,
        "ligand_value": ligand_value,
        "peptide_sequence": peptide_sequence,
        "cyclotide_sequence": cyclotide_sequence,
        "redesign_residues": redesign_residues,
        "design_yaml": design_yaml,
        "yaml_protocol": yaml_protocol,
        "mimosa_enabled": mimosa_enabled,
        "mimosa_condition": mimosa_condition,
        "mimosa_complex": str(mimosa_complex) if mimosa_complex else None,
        "mimosa_binder_chain": mimosa_binder_chain,
        "mimosa_complex_target_chain": mimosa_complex_target_chain,
        "mimosa_target_chain": mimosa_target_chain,
        "mimosa_cdr": mimosa_cdr,
        "mimosa_design_loop": mimosa_design_loop,
        "mimosa_source_positions": mimosa_source_positions,
        "mimosa_contact_cutoff": float(mimosa_contact_cutoff),
        "mimosa_guidance_scale": float(mimosa_guidance_scale),
        "mimosa_num_particles": int(mimosa_num_particles),
        "mimosa_fk_lambda": float(mimosa_fk_lambda),
        "mimosa_resampling_interval": int(mimosa_resampling_interval),
    }
    job = create_job(
        DESIGN_GROUP,
        "design_generation" if generator_only_tile else "design_campaign",
        "boltzgen_generator" if generator_only_tile else "boltzgen",
        {"target_pdb": str(target_pdb) if target_pdb else "", "target_chains": target_chains},
        params,
    )
    raw_dir = job.run_dir / "artifacts" / "raw" / "boltzgen"
    raw_dir.mkdir(parents=True, exist_ok=True)
    target_artifact: Path | None = None
    residue_index_by_chain: dict[str, dict[int, int]] = {}
    if target_pdb is not None:
        is_mmcif = target_pdb.suffix.casefold() in {".cif", ".mmcif"}
        if is_mmcif:
            target_artifact = artifact_path(job.run_dir, "raw/boltzgen/input", "target.cif")
            shutil.copy2(target_pdb, target_artifact)
            mmcif_residues = mmcif_protein_residues_by_chain(target_pdb)
            for chain_id, entries in mmcif_residues.items():
                chain_indices: dict[int, int] = {}
                for index, (residue_number, _insertion_code, _residue_name) in enumerate(entries, start=1):
                    chain_indices.setdefault(residue_number, index)
                residue_index_by_chain[chain_id] = chain_indices
        else:
            target_artifact = _copy_target_for_design(job.run_dir, target_pdb, target_chains, "raw/boltzgen/input")
            summary = pdb_summary(target_artifact.read_text(errors="ignore"))
            residues_by_chain = {row["chain_id"]: sorted(row.get("residues") or []) for row in summary.get("chains", [])}
            residue_index_by_chain = {
                chain_id: {residue: index for index, residue in enumerate(residues, start=1)}
                for chain_id, residues in residues_by_chain.items()
            }
    mimosa_source_artifact: Path | None = None
    if mimosa_enabled:
        source_path = Path(mimosa_complex).expanduser()
        source_suffix = source_path.suffix.lower()
        mimosa_source_artifact = artifact_path(
            job.run_dir, "raw/boltzgen/input", f"mimosa_source_complex{source_suffix}"
        )
        shutil.copy2(source_path, mimosa_source_artifact)
        if run_vanilla_pipeline:
            integrity_script = Path(__file__).resolve().parents[4] / "scripts" / "check_mimosa_backbone.py"
            if not integrity_script.is_file():
                raise ValueError("The MIMOSA peptide-bond integrity checker is missing from the Studio source tree.")
            shutil.copy2(integrity_script, raw_dir / "input" / "check_mimosa_backbone.py")
    binding_by_chain: dict[str, list[str]] = {}
    boltzgen_hotspot_map: dict[str, str] = {}
    unmapped_hotspots: list[str] = []

    def map_residues(text: str) -> tuple[dict[str, list[str]], dict[str, str], list[str]]:
        mapped_by_chain: dict[str, list[str]] = {}
        residue_map: dict[str, str] = {}
        unmapped: list[str] = []
        for token in text.split(",") if text else []:
            if not token:
                continue
            chain_id = next(
                (
                    chain
                    for chain in sorted(target_chains, key=len, reverse=True)
                    if token.startswith(chain) and token[len(chain) :].isdigit()
                ),
                "",
            )
            if not chain_id:
                unmapped.append(token)
                continue
            residue_number = int(token[len(chain_id) :])
            chain_local_index = residue_index_by_chain.get(chain_id, {}).get(residue_number)
            if chain_local_index is None:
                unmapped.append(token)
                continue
            mapped = str(chain_local_index)
            mapped_by_chain.setdefault(chain_id, []).append(mapped)
            residue_map[token] = f"{chain_id}{mapped}"
        return mapped_by_chain, residue_map, unmapped

    if task == "redesign" and not hotspots:
        hotspots = redesign_residues
        params["hotspots"] = hotspots
    binding_by_chain, boltzgen_hotspot_map, unmapped_hotspots = map_residues(hotspots)
    non_binding_by_chain, boltzgen_non_binding_map, unmapped_non_binding = map_residues(
        non_binding_residues,
    )
    for chain_id in set(binding_by_chain) & set(non_binding_by_chain):
        overlap = set(binding_by_chain[chain_id]) & set(non_binding_by_chain[chain_id])
        if overlap:
            raise ValueError(
                f"BoltzGen residues on chain {chain_id} cannot be both binding and non-binding: "
                f"{', '.join(sorted(overlap, key=int))}."
            )
    params["boltzgen_hotspot_map"] = boltzgen_hotspot_map
    params["boltzgen_binding_by_chain"] = binding_by_chain
    params["unmapped_hotspots"] = unmapped_hotspots
    params["boltzgen_non_binding_map"] = boltzgen_non_binding_map
    params["boltzgen_non_binding_by_chain"] = non_binding_by_chain
    params["unmapped_non_binding_residues"] = unmapped_non_binding
    write_json(job.run_dir / "input.json", {"inputs": {"target_pdb": str(target_pdb) if target_pdb else "", "target_chains": target_chains}, "params": params})
    write_json(
        raw_dir / "input" / "hotspot_mapping.json",
        {
            "requested_hotspots": [token for token in hotspots.split(",") if token],
            "mapped_hotspots": boltzgen_hotspot_map,
            "binding_by_chain": binding_by_chain,
            "unmapped_hotspots": unmapped_hotspots,
            "requested_non_binding_residues": [token for token in non_binding_residues.split(",") if token],
            "mapped_non_binding_residues": boltzgen_non_binding_map,
            "non_binding_by_chain": non_binding_by_chain,
            "unmapped_non_binding_residues": unmapped_non_binding,
            "mapping_contract": "PDB/mmCIF author residue number -> BoltzGen chain-local residue index",
        },
    )
    input_dir = raw_dir / "input"
    input_dir.mkdir(parents=True, exist_ok=True)
    suggested_design_chain_id = next(
        (chain_id for chain_id in "ABCDEFGHIJKLMNOPQRSTUVWXYZ" if chain_id not in set(target_chains)),
        "0",
    )
    if task == "yaml":
        binding_lines: list[str] = ["      binding_types:"] if binding_by_chain or non_binding_by_chain else []
        for chain_id in sorted(set(binding_by_chain) | set(non_binding_by_chain)):
            binding_lines.extend(["        - chain:", f"            id: {json.dumps(chain_id)}"])
            if binding_by_chain.get(chain_id):
                binding_lines.append(f"            binding: {','.join(binding_by_chain[chain_id])}")
            if non_binding_by_chain.get(chain_id):
                binding_lines.append(f"            not_binding: {','.join(non_binding_by_chain[chain_id])}")
        binding_spec = "\n".join(binding_lines) + "\n" if binding_lines else ""
        yaml_text = design_yaml.replace("{{TARGET_FILE}}", target_artifact.name if target_artifact else "")
        yaml_text = yaml_text.replace("{{TARGET_CHAIN}}", target_chains[0] if target_chains else "")
        yaml_text = yaml_text.replace("{{DESIGN_CHAIN}}", suggested_design_chain_id)
        yaml_text = yaml_text.replace("{{BINDING_TYPES}}", binding_spec)
    elif task == "small_molecule":
        ligand_key = "ccd" if ligand_type == "ccd" else "smiles"
        ligand_scalar = str(ligand_value).strip().upper() if ligand_type == "ccd" else json.dumps(str(ligand_value).strip())
        yaml_text = (
            "entities:\n"
            "  - protein:\n"
            "      id: A\n"
            f"      sequence: {binder_sequence_spec}\n"
            "  - ligand:\n"
            "      id: B\n"
            f"      {ligand_key}: {ligand_scalar}\n"
        )
    elif task == "nanobody":
        for suffix in (".yaml", ".cif"):
            source = Path(__file__).parent / "../assets" / "boltzgen_nanobody_scaffolds" / f"{scaffold_name}{suffix}"
            if not source.exists():
                raise ValueError(f"Bundled BoltzGen scaffold file is missing: {source.name}")
            shutil.copy2(source, input_dir / source.name)
        yaml_text = "entities:\n  - file:\n      path: " + target_artifact.name + "\n      include:\n"
        for chain_id in target_chains:
            yaml_text += f"        - chain:\n            id: {json.dumps(chain_id)}\n"
        if binding_by_chain or non_binding_by_chain:
            yaml_text += "      binding_types:\n"
            for chain_id in sorted(set(binding_by_chain) | set(non_binding_by_chain)):
                yaml_text += f"        - chain:\n            id: {json.dumps(chain_id)}\n"
                if binding_by_chain.get(chain_id):
                    yaml_text += f"            binding: {','.join(binding_by_chain[chain_id])}\n"
                if non_binding_by_chain.get(chain_id):
                    yaml_text += f"            not_binding: {','.join(non_binding_by_chain[chain_id])}\n"
        yaml_text += f"  - file:\n      path: {scaffold_name}.yaml\n"
    elif task == "redesign":
        yaml_text = "entities:\n  - file:\n      path: " + target_artifact.name + "\n      include:\n"
        for chain_id in target_chains:
            yaml_text += f"        - chain:\n            id: {json.dumps(chain_id)}\n"
        yaml_text += "      design:\n"
        for chain_id, indices in sorted(binding_by_chain.items()):
            if indices:
                yaml_text += (
                    "        - chain:\n"
                    f"            id: {json.dumps(chain_id)}\n"
                    f"            res_index: {','.join(indices)}\n"
                )
        yaml_text += '      structure_groups: "all"\n'
    else:
        designed_chain_id = suggested_design_chain_id
        sequence_spec = (
            cyclotide_sequence.strip()
            if task == "cyclotide"
            else peptide_sequence.strip()
            if task == "peptide"
            else binder_sequence_spec
        )
        cyclic_line = "      cyclic: true\n" if task == "cyclotide" else ""
        yaml_text = (
            "entities:\n"
            "  - protein:\n"
            f"      id: {designed_chain_id}\n"
            f"      sequence: {sequence_spec}\n"
            f"{cyclic_line}"
            "  - file:\n"
            f"      path: {target_artifact.name}\n"
            "      include:\n"
        )
        for chain_id in target_chains:
            yaml_text += f"        - chain:\n            id: {json.dumps(chain_id)}\n"
        if binding_by_chain or non_binding_by_chain:
            yaml_text += "      binding_types:\n"
            for chain_id in sorted(set(binding_by_chain) | set(non_binding_by_chain)):
                yaml_text += f"        - chain:\n            id: {json.dumps(chain_id)}\n"
                if binding_by_chain.get(chain_id):
                    yaml_text += f"            binding: {','.join(binding_by_chain[chain_id])}\n"
                if non_binding_by_chain.get(chain_id):
                    yaml_text += f"            not_binding: {','.join(non_binding_by_chain[chain_id])}\n"
        yaml_text += '      structure_groups: "all"\n'
        if task == "cyclotide":
            yaml_text += (
                "constraints:\n"
                f"  - bond:\n      atom1: [{designed_chain_id}, 4, SG]\n      atom2: [{designed_chain_id}, 26, SG]\n"
                f"  - bond:\n      atom1: [{designed_chain_id}, 13, SG]\n      atom2: [{designed_chain_id}, 30, SG]\n"
                f"  - bond:\n      atom1: [{designed_chain_id}, 20, SG]\n      atom2: [{designed_chain_id}, 32, SG]\n"
            )
    yaml_path = input_dir / "target_binder.yaml"
    yaml_path.write_text(yaml_text.rstrip() + "\n")
    if task == "yaml":
        yaml_normalizer_path = input_dir / "normalize_boltzgen_yaml.py"
        yaml_normalizer_path.write_text(
            "from pathlib import Path\n"
            "import yaml\n"
            "path = Path('/work/artifacts/raw/boltzgen/input/target_binder.yaml')\n"
            "data = yaml.safe_load(path.read_text())\n"
            "if not isinstance(data, dict):\n"
            "    raise ValueError('BoltzGen YAML must contain a top-level mapping')\n"
            "if 'entities' not in data:\n"
            "    aliases = [key for key in data if str(key).rstrip('/').endswith('/entities')]\n"
            "    if len(aliases) != 1:\n"
            "        raise ValueError('BoltzGen YAML must have a top-level entities list')\n"
            "    data['entities'] = data.pop(aliases[0])\n"
            "    path.write_text(yaml.safe_dump(data, sort_keys=False))\n"
            "    print('Normalized BoltzGen entities key from the mounted YAML path')\n"
        )
    if run_vanilla_pipeline:
        filter_patch_path = input_dir / "patch_optional_esmfold2_filter.py"
        filter_patch_path.write_text(
            "from importlib.util import find_spec\n"
            "from pathlib import Path\n"
            "module = find_spec('boltzgen.task.filter.filter')\n"
            "if module is None or module.origin is None:\n"
            "    raise RuntimeError('BoltzGen Filter source is unavailable')\n"
            "path = Path(module.origin)\n"
            "source = path.read_text()\n"
            "replacements = [\n"
            "    ('        use_affinity: bool = False,  # This changes the filtering metrics to metrics more amenable to small molecule binder design\\n',\n"
            "     '        use_affinity: bool = False,  # This changes the filtering metrics to metrics more amenable to small molecule binder design\\n        use_esmfold2: bool | None = None,\\n'),\n"
            "    ('        self.use_esmfold2 = not use_affinity\\n',\n"
            "     '        self.use_esmfold2 = (not use_affinity) if use_esmfold2 is None else bool(use_esmfold2)\\n'),\n"
            "    ('        else:\\n            del self.metrics[\"design_to_target_iptm\"]\\n            del self.metrics[\"neg_min_design_to_target_pae\"]\\n            self.metrics[self.esmfold2_score_key] = 1\\n',\n"
            "     '        elif self.use_esmfold2:\\n            del self.metrics[\"design_to_target_iptm\"]\\n            del self.metrics[\"neg_min_design_to_target_pae\"]\\n            self.metrics[self.esmfold2_score_key] = 1\\n'),\n"
            "]\n"
            "for old, new in replacements:\n"
            "    if source.count(old) != 1:\n"
            "        raise RuntimeError('Pinned BoltzGen Filter patch anchor not found exactly once')\n"
            "    source = source.replace(old, new, 1)\n"
            "path.write_text(source)\n"
            "print('Applied native-metric fallback for missing optional ESMFold2 runtime')\n"
        )
    output_dir = "run-vanilla" if run_vanilla_pipeline else "run-generation-only"
    run_args = (
        f"boltzgen run input/target_binder.yaml --output {output_dir} --protocol {protocol} "
        f"--num_designs {num_designs} --devices 1 --num_workers 1 --cache /cache "
    )
    if diffusion_batch_size is not None:
        run_args += f"--diffusion_batch_size {diffusion_batch_size} "
    if omit_amino_acids:
        run_args += f"--inverse_fold_avoid {shlex.quote(omit_amino_acids)} "
    if run_vanilla_pipeline:
        run_args += f"--budget {budget} "
        # The optimized image does not include BoltzGen's external uv-managed
        # ESMFold2 runtime. Keep the native Boltz2 refolding, metric, and ranking
        # stages, which provide the complex PAE output used by Studio.
        steps = ["design", "inverse_folding", "folding", "design_folding"]
        if protocol == "protein-small_molecule":
            steps.append("affinity")
        steps.extend(["analysis", "filtering"])
        run_args += "--steps " + " ".join(steps) + " "
        # The default analysis config expects ESMFold2 score files even when
        # that optional runtime is absent. Disable only those scores so native
        # BoltzGen analysis and ranking can complete.
        run_args += "--config analysis esmfold2_metrics=false "
        run_args += "--config filtering use_esmfold2=false "
    else:
        run_args += "--steps design "
    if additional_filters:
        run_args += "--additional_filters " + " ".join(
            shlex.quote(value) for value in additional_filters
        ) + " "
    if filter_binding_site_adherence:
        run_args += "--config filtering filter_bindingsite=true "
    run_args += (
        f"--config design sampling_steps={sampling_steps} "
        "compile_pairformer=false compile_structure=false data.num_workers=1"
    )
    if run_vanilla_pipeline:
        # BoltzGen computes the native ipSAE metrics during refolding. This
        # optimized image's writer extension persists the full PAE matrix.
        pae_metric_keys = (
            "pae,design_ipsae_min,design_to_target_ipsae,target_to_design_ipsae,"
            "min_interaction_pae,min_design_to_target_pae,interaction_pae,"
            "ligand_iptm,protein_iptm,iptm,design_iptm,design_iiptm,"
            "design_to_target_iptm,design_residue_iptm,design_ptm,target_ptm,ptm"
        )
        run_args += (
            " --config folding "
            + shlex.quote(f"keys_dict_out=[{pae_metric_keys}]")
            + " "
            + shlex.quote("writer.extra_eval_keys=[pae]")
        )
    if run_vanilla_pipeline:
        pae_archive_check = shlex.quote(
            "import numpy as np, pathlib; "
            "archives = pathlib.Path('run-vanilla').rglob('fold_out_npz/*.npz'); "
            "assert any('pae' in np.load(path).files for path in archives), "
            "'BoltzGen did not write a full PAE matrix archive'"
        )
        output_check = (
            f"test -n \"$(find {output_dir}/final_ranked_designs -type f -path '*/final_*_designs/*.cif' -print -quit)\" "
            f"&& python -c {pae_archive_check}"
        )
    else:
        output_check = f"test -n \"$(find {output_dir}/intermediate_designs -maxdepth 1 -type f -name '*.cif' -print -quit)\""
    mimosa_prepare = ""
    if mimosa_enabled:
        if target_artifact is None or mimosa_source_artifact is None:
            raise ValueError("MIMOSA needs both the target structure and cognate complex artifacts.")
        motif_command = [
            "python", "/workspace/scripts/prepare_mimosa_motif.py",
            "--complex", f"/work/artifacts/raw/boltzgen/{mimosa_source_artifact.relative_to(raw_dir)}",
            "--target", f"/work/artifacts/raw/boltzgen/{target_artifact.relative_to(raw_dir)}",
            "--binder-chain", mimosa_binder_chain,
            "--complex-target-chain", mimosa_complex_target_chain,
            "--target-chain", mimosa_target_chain,
            "--cdr", mimosa_cdr,
            "--contact-cutoff", str(mimosa_contact_cutoff),
            "--guidance-scale", str(mimosa_guidance_scale),
            "--fk-lambda", str(mimosa_fk_lambda),
            "--resampling-interval", str(mimosa_resampling_interval),
            "--num-particles", str(mimosa_num_particles),
            "--condition", mimosa_condition,
            "--design-loop", mimosa_design_loop,
            "--output", "/work/artifacts/raw/boltzgen/input/target_binder.mimosa.json",
        ]
        if mimosa_source_positions:
            motif_command.extend(["--source-residues", *mimosa_source_positions])
        mimosa_prepare = shlex.join(motif_command) + "; "
    mimosa_split_pipeline = ""
    if mimosa_enabled and run_vanilla_pipeline:
        all_steps = "--steps " + " ".join(steps)
        generation_args = run_args.replace(all_steps, "--steps design")
        post_generation_steps = [step for step in steps if step != "design"]
        if not post_generation_steps:
            raise ValueError("BoltzGen MIMOSA needs inverse-folding/refolding stages after the integrity check.")
        downstream_args = run_args.replace(all_steps, "--steps " + " ".join(post_generation_steps))
        integrity_report = "run-vanilla/mimosa_integrity_results.json"
        integrity_check = shlex.quote(
            "import json,sys; d=json.load(open(" + repr(integrity_report) + ")); "
            "sys.exit(0 if d.get('candidate_count',0)>0 and d.get('candidate_count')==d.get('rejected_count') else 1)"
        )
        target_inside = f"/work/artifacts/raw/boltzgen/input/{target_artifact.name}"
        checker_inside = "/work/artifacts/raw/boltzgen/input/check_mimosa_backbone.py"
        mimosa_split_pipeline = (
            f"{generation_args}; "
            f"python {shlex.quote(checker_inside)} run-vanilla/intermediate_designs {shlex.quote(target_inside)} {shlex.quote(integrity_report)}; "
            f"if python -c {integrity_check}; then echo '[MIMOSA] All generated backbones were rejected by the peptide-bond filter'; "
            f"else {downstream_args}; fi; "
        )
        output_check = (
            f"if python -c {integrity_check}; then echo '[MIMOSA] Structural rejection report is the completed result'; "
            f"else test -n \"$(find run-vanilla/final_ranked_designs -type f -path '*/final_*_designs/*.cif' -print -quit)\" "
            f"&& python -c {pae_archive_check}; fi"
        )
    script = (
        "set -euxo pipefail; cd /work/artifacts/raw/boltzgen; "
        + mimosa_prepare
        + (
            "python /work/artifacts/raw/boltzgen/input/normalize_boltzgen_yaml.py; "
            if task == "yaml"
            else ""
        )
        + (
            "python /work/artifacts/raw/boltzgen/input/patch_optional_esmfold2_filter.py; "
            if run_vanilla_pipeline
            else ""
        )
        + f"boltzgen check input/target_binder.yaml --output checked --cache /cache; rm -rf {output_dir}; "
        + (mimosa_split_pipeline if mimosa_split_pipeline else f"{run_args}; ")
        + f"{output_check}"
    )
    command = [
        "docker",
        "run",
        "--rm",
        *docker_gpu_args(gpu_device),
        "--shm-size=64G",
        "--entrypoint",
        "/bin/bash",
        "-v",
        "/mnt/db/reference_files/boltzgen-cache:/cache",
        "-v",
        f"{job.run_dir}:/work",
        "-w",
        "/work",
        manifest["image"],
        "-lc",
        script,
    ]
    rc = _run_shell_steps(job.run_dir, [{"name": "boltzgen", "command": command}])
    if rc == 0 and run_vanilla_pipeline:
        candidates = _normalize_boltzgen_vanilla_candidates(job.run_dir, params, target_artifact)
    elif rc == 0:
        candidates = _normalize_boltzgen_generation_candidates(job.run_dir, params, target_artifact)
    else:
        candidates = []
    _finish_design_job(
        job.run_dir,
        rc == 0 and bool(candidates),
        rc,
        candidates,
        [
            ("raw/boltzgen/**/*.cif", "cif"),
            ("raw/boltzgen/**/*.npz", "npz"),
            ("raw/boltzgen/**/*.json", "json"),
            ("raw/boltzgen/**/*.csv", "csv"),
            ("raw/boltzgen/**/*.pdf", "pdf"),
            ("raw/boltzgen/**/*.yaml", "yaml"),
        ],
    )
    return job.run_dir


def run_genie3(
    target_pdb: Path,
    target_chains: list[str],
    binder_length: str,
    hotspots: str = "",
    campaign_name: str = "",
    num_designs: int = 1,
    seed: int = 7,
    num_devices: int = 1,
    direction_scale: float = 0.0,
    cond_strategy: str = "extended",
    inverse_folding_num_seq: int = 1,
    folding_model_name: str = "colabfold",
    folding_mode: str = "template",
    folding_num_models: int = 5,
    folding_num_recycles: int = 20,
    compile_generation: bool = False,
    run_mode: str = "full_vanilla_pipeline",
    enable_beam_search: bool = False,
    beam_width: int = 4,
    gpu_device: object = "0",
) -> Path:
    if not target_chains:
        raise ValueError("At least one target chain is required.")
    hotspots = normalize_hotspots(hotspots) if hotspots else ""
    parse_binder_lengths(binder_length)
    if num_designs < 1:
        raise ValueError("Design attempts must be at least 1.")
    if num_devices < 1:
        raise ValueError("Number of devices must be at least 1.")
    if inverse_folding_num_seq < 1:
        raise ValueError("inverse_folding.num_seq must be at least 1.")
    if cond_strategy not in {"hotspot", "extended", "common", "iter_common", "iter_common_prob"}:
        raise ValueError("Unknown Genie3 conditioning strategy.")
    if folding_model_name not in {"colabfold", "boltz2"}:
        raise ValueError("Genie3 folding model must be colabfold or boltz2.")
    if folding_mode not in {"template", "msa"}:
        raise ValueError("Genie3 folding mode must be template or msa.")
    if folding_model_name == "boltz2" and folding_mode != "msa":
        raise ValueError("Genie3 Boltz2 binder evaluation requires folding.mode=msa.")
    if run_mode not in {"generation_only", "full_vanilla_pipeline"}:
        raise ValueError("Unknown Genie3 run mode.")
    manifest = load_manifest("genie3")
    params = {
        "target_chains": target_chains,
        "binder_length": binder_length,
        "hotspots": hotspots,
        "campaign_name": campaign_name.strip(),
        "num_designs": num_designs,
        "seed": seed,
        "num_devices": num_devices,
        "direction_scale": direction_scale,
        "cond_strategy": cond_strategy,
        "inverse_folding_num_seq": inverse_folding_num_seq,
        "folding_model_name": folding_model_name,
        "folding_mode": folding_mode,
        "folding_num_models": folding_num_models,
        "folding_num_recycles": folding_num_recycles,
        "compile_generation": compile_generation,
        "run_mode": run_mode,
        "enable_beam_search": enable_beam_search,
        "beam_width": beam_width,
        "gpu_device": normalize_gpu_device(gpu_device),
        "pipeline_mode": "genie3_vanilla_binder",
    }
    job = create_job(
        DESIGN_GROUP,
        "design_campaign",
        "genie3",
        {"target_pdb": str(target_pdb), "target_chains": target_chains},
        params,
    )
    if campaign_name.strip():
        update_status(job.run_dir, "queued", campaign_name=campaign_name.strip())
    target_artifact = _copy_target_for_design(job.run_dir, target_pdb, target_chains, "raw/genie3/input")
    raw_dir = job.run_dir / "artifacts" / "raw" / "genie3"
    raw_dir.mkdir(parents=True, exist_ok=True)
    boltz_msa_by_source_chain: dict[str, str] = {}
    if folding_model_name == "boltz2":
        try:
            update_status(job.run_dir, "running")
            boltz_msa_by_source_chain = _ensure_boltz_msas_for_target(
                job.run_dir,
                target_artifact,
                target_chains,
                gpu_device=gpu_device,
            )
        except Exception as exc:
            with (job.run_dir / "stderr.log").open("a") as stderr:
                stderr.write(f"Boltz2 MSA cache preparation failed: {exc}\n")
            _finish_design_job(
                job.run_dir,
                False,
                1,
                [],
                [
                    ("raw/genie3/**/*.pdb", "pdb"),
                    ("raw/genie3/**/*.json", "json"),
                    ("raw/genie3/**/*.yaml", "yaml"),
                    ("raw/genie3/**/*.a3m", "a3m"),
                    ("raw/genie3/**/*.fasta", "fasta"),
                ],
            )
            return job.run_dir
    dataset_dir, problem_key, problem = _write_genie3_target_files(
        run_dir=job.run_dir,
        target_artifact=target_artifact,
        target_chains=target_chains,
        binder_length=binder_length,
        hotspots=hotspots,
        campaign_name=campaign_name,
        boltz_msa_by_source_chain=boltz_msa_by_source_chain,
    )
    params["problem_key"] = problem_key
    write_json(job.run_dir / "input.json", {"inputs": {"target_pdb": str(target_pdb), "target_chains": target_chains}, "params": params})
    config_path = raw_dir / "experiment.yaml"
    experiment_name = _clean_key(campaign_name or f"genie3_{job.run_dir.name}")
    _write_genie3_experiment_yaml(
        config_path,
        experiment_name=experiment_name,
        rootdir="/work/artifacts/raw/genie3/output",
        dataset=f"/work/{_rel_path(job.run_dir, dataset_dir)}",
        problem_key=problem_key,
        num_designs=num_designs,
        seed=seed,
        num_devices=num_devices,
        direction_scale=direction_scale,
        inverse_folding_num_seq=inverse_folding_num_seq,
        folding_model_name=folding_model_name,
        folding_mode=folding_mode,
        folding_num_models=folding_num_models,
        folding_num_recycles=folding_num_recycles,
        compile_generation=compile_generation,
        run_mode=run_mode,
        cond_strategy=cond_strategy,
        enable_beam_search=enable_beam_search,
        beam_width=beam_width,
    )
    write_json(raw_dir / "problem.json", problem)
    if run_mode == "full_vanilla_pipeline":
        genie3_args = (
            "genie3 generate -c /work/artifacts/raw/genie3/experiment.yaml "
            f"--num-devices {num_devices} --verbose && "
            "genie3 evaluate -c /work/artifacts/raw/genie3/experiment.yaml "
            f"--num-devices {num_devices} --verbose"
        )
        entrypoint_args = ["bash", "-lc", genie3_args]
    else:
        entrypoint_args = [
            "genie3",
            "generate",
            "-c",
            "/work/artifacts/raw/genie3/experiment.yaml",
            "--num-devices",
            str(num_devices),
            "--verbose",
        ]
    command = [
        "docker",
        "run",
        "--rm",
        *docker_gpu_args(gpu_device),
        "--shm-size=64G",
        "-v",
        "/mnt/db/reference_files/genie3/pretrained:/opt/genie3/pretrained:ro",
        "-v",
        f"{job.run_dir}:/work",
        "-v",
        f"{_ALPHAFOLD_MODELS_DIR}:/alphafold_models:ro",
        "-v",
        f"{_BOLTZ_CACHE_DIR}:/cache",
        "-v",
        f"{_BOLTZ_MSA_REPOSITORY_DIR}:/msa_repository:ro",
        "-e",
        "GENIE3_COLABFOLD_DATA_DIR=/alphafold_models",
        "-e",
        "BOLTZ_CACHE=/cache",
        "-w",
        "/opt/genie3",
        manifest["image"],
        *entrypoint_args,
    ]
    rc = _run_shell_steps(job.run_dir, [{"name": "genie3", "command": command}])
    candidates = _normalize_genie3_candidates(job.run_dir, params, target_artifact) if rc == 0 else []
    _finish_design_job(
        job.run_dir,
        rc == 0 and bool(candidates),
        rc,
        candidates,
        [
            ("raw/genie3/**/*.pdb", "pdb"),
            ("raw/genie3/**/*.cif", "cif"),
            ("raw/genie3/**/*.csv", "csv"),
            ("raw/genie3/**/*.json", "json"),
            ("raw/genie3/**/*.yaml", "yaml"),
            ("raw/genie3/**/*.a3m", "a3m"),
            ("raw/genie3/**/*.fasta", "fasta"),
        ],
    )
    return job.run_dir


def run_pxdesign(
    target_pdb: Path,
    target_chains: list[str],
    binder_length: str,
    hotspots: str = "",
    campaign_name: str = "",
    num_designs: int = 1,
    n_steps: int = 400,
    dtype: str = "bf16",
    preset: str = "preview",
    run_mode: str = "pipeline",
    n_max_runs: int = 1,
    use_fast_ln: bool = True,
    use_deepspeed_evo_attention: bool = False,
    prepare_target_msa: bool = True,
    gpu_device: object = "0",
) -> Path:
    hotspots = normalize_hotspots(hotspots) if hotspots else ""
    lengths = parse_binder_lengths(binder_length)
    px_binder_length = int(lengths[0])
    manifest = load_manifest("pxdesign")
    params = {
        "target_chains": target_chains,
        "binder_length": binder_length,
        "pxdesign_binder_length": px_binder_length,
        "hotspots": hotspots,
        "campaign_name": campaign_name.strip(),
        "num_designs": num_designs,
        "n_steps": n_steps,
        "dtype": dtype,
        "preset": preset,
        "run_mode": run_mode,
        "n_max_runs": n_max_runs,
        "use_fast_ln": use_fast_ln,
        "use_deepspeed_evo_attention": use_deepspeed_evo_attention,
        "prepare_target_msa": prepare_target_msa,
        "gpu_device": normalize_gpu_device(gpu_device),
    }
    job = create_job(DESIGN_GROUP, "design_campaign", "pxdesign", {"target_pdb": str(target_pdb), "target_chains": target_chains}, params)
    if campaign_name.strip():
        update_status(job.run_dir, "queued", campaign_name=campaign_name.strip())
    raw_dir = job.run_dir / "artifacts" / "raw" / "pxdesign"
    input_dir = raw_dir / "input"
    input_dir.mkdir(parents=True, exist_ok=True)
    target_artifact = _copy_target_for_design(job.run_dir, target_pdb, target_chains, "raw/pxdesign/input")
    pxdesign_msa_dirs: dict[str, str] = {}
    target_msa_error = ""
    if prepare_target_msa and run_mode != "generation_only":
        try:
            pxdesign_msa_dirs = ensure_pxdesign_msa_dirs_for_target(
                job.run_dir,
                target_artifact,
                target_chains,
                raw_subdir="pxdesign/input/msa",
                yaml_base_dir=raw_dir,
                gpu_device=gpu_device,
            )
            params["pxdesign_msa_dirs"] = pxdesign_msa_dirs
        except Exception as exc:
            target_msa_error = str(exc)
            params["target_msa_error"] = target_msa_error
            with (job.run_dir / "stderr.log").open("a") as stderr:
                stderr.write(f"PXDesign target MSA preparation failed: {type(exc).__name__}: {exc}\n")
    hotspots_by_chain: dict[str, list[int]] = {}
    for token in hotspots.split(",") if hotspots else []:
        if not token:
            continue
        try:
            hotspots_by_chain.setdefault(token[0], []).append(int(token[1:]))
        except ValueError:
            continue
    yaml_lines = [
        "target:",
        "  file: input/target.pdb",
        "  chains:",
    ]
    for chain_id in target_chains or ["A"]:
        chain_hotspots = hotspots_by_chain.get(chain_id, [])
        chain_msa = pxdesign_msa_dirs.get(chain_id)
        if chain_hotspots or chain_msa:
            yaml_lines.extend(
                [
                    f"    {chain_id}:",
                ]
            )
            if chain_hotspots:
                yaml_lines.append(f"      hotspots: [{', '.join(str(item) for item in chain_hotspots)}]")
            if chain_msa:
                yaml_lines.append(f"      msa: {chain_msa}")
        else:
            yaml_lines.append(f"    {chain_id}: all")
    yaml_lines.extend(["", f"binder_length: {px_binder_length}"])
    (input_dir / "target_binder.yaml").write_text("\n".join(yaml_lines) + "\n")
    output_dir = "output"
    common_prefix = (
        "set -euxo pipefail; cd /work/artifacts/raw/pxdesign; "
        "export PROTENIX_DATA_ROOT_DIR=/ref/pxdesign/release_data/ccd_cache; "
        "export TOOL_WEIGHTS_ROOT=/ref/pxdesign/tool_weights; "
        "pxdesign check-input --yaml input/target_binder.yaml; "
        "rm -rf parsed_target output; "
        "pxdesign parse-target --yaml input/target_binder.yaml -o parsed_target; "
    )
    if run_mode == "generation_only":
        script = (
            common_prefix +
            f"pxdesign infer -i input/target_binder.yaml -o {output_dir} "
            f"--N_sample {int(num_designs)} "
            f"--N_step {int(n_steps)} "
            f"--dtype {shlex.quote(dtype)} "
            f"--use_fast_ln {str(bool(use_fast_ln))} "
            f"--use_deepspeed_evo_attention {str(bool(use_deepspeed_evo_attention))} "
            "--load_checkpoint_dir /ref/pxdesign/release_data/checkpoint; "
            'test -n "$(find output -type f -name \"*.cif\" -print -quit)"'
        )
    else:
        script = (
            common_prefix +
            f"pxdesign pipeline --preset {shlex.quote(preset)} "
            "-i input/target_binder.yaml "
            f"-o {output_dir} "
            f"--N_sample {int(num_designs)} "
            f"--N_step {int(n_steps)} "
            f"--N_max_runs {int(n_max_runs)} "
            f"--dtype {shlex.quote(dtype)} "
            f"--use_fast_ln {str(bool(use_fast_ln))} "
            f"--use_deepspeed_evo_attention {str(bool(use_deepspeed_evo_attention))} "
            "--load_checkpoint_dir /ref/pxdesign/release_data/checkpoint; "
            "test -s output/design_outputs/target_binder/summary.csv"
        )
    command = [
        "docker",
        "run",
        "--rm",
        *docker_gpu_args(gpu_device),
        "-v",
        f"{job.run_dir}:/work",
        "-v",
        "/mnt/db/reference_files/pxdesign:/ref/pxdesign",
        "-w",
        "/work",
        manifest["image"],
        "bash",
        "-lc",
        script,
    ]
    rc = _run_shell_steps(job.run_dir, [{"name": "pxdesign", "command": command}])
    if rc == 0 and run_mode == "generation_only":
        candidates = _normalize_pxdesign_generation_candidates(job.run_dir, params, target_artifact)
    elif rc == 0:
        candidates = _normalize_pxdesign_candidates(job.run_dir, params, target_artifact)
    else:
        candidates = []
    _finish_design_job(
        job.run_dir,
        rc == 0 and bool(candidates),
        rc,
        candidates,
        [
            ("raw/pxdesign/**/*.cif", "cif"),
            ("raw/pxdesign/**/*.pdb", "pdb"),
            ("raw/pxdesign/**/*.csv", "csv"),
            ("raw/pxdesign/**/*.json", "json"),
            ("raw/pxdesign/**/*.yaml", "yaml"),
            ("raw/pxdesign/**/*.a3m", "a3m"),
            ("raw/pxdesign/**/*.png", "plot"),
        ],
    )
    return job.run_dir


def _write_protpardelle_target_motif(
    source_pdb: Path,
    target_chains: list[str],
    motif_path: Path,
) -> tuple[dict[str, int], dict[str, str], dict[tuple[str, int], tuple[str, int]]]:
    """Write a Protpardelle motif PDB with compact chain IDs and 1-based residues."""
    chain_ids = [chain for chain in target_chains if chain]
    chain_map = {old_chain: chr(ord("A") + index) for index, old_chain in enumerate(chain_ids)}
    residue_map: dict[tuple[str, int], tuple[str, int]] = {}
    next_residue: dict[str, int] = {old_chain: 0 for old_chain in chain_ids}
    seen_residues: set[tuple[str, str, str]] = set()
    lengths: dict[str, int] = {}
    output_lines: list[str] = []

    for line in source_pdb.read_text(errors="ignore").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        old_chain = line[21].strip() or "_"
        if old_chain not in chain_map:
            continue
        try:
            old_residue = int(line[22:26])
        except ValueError:
            continue
        residue_key = (old_chain, line[22:26], line[26])
        if residue_key not in seen_residues:
            seen_residues.add(residue_key)
            next_residue[old_chain] += 1
            residue_map[(old_chain, old_residue)] = (chain_map[old_chain], next_residue[old_chain])
            lengths[chain_map[old_chain]] = next_residue[old_chain]
        new_chain, new_residue = residue_map[(old_chain, old_residue)]
        output_lines.append(f"{line[:21]}{new_chain}{new_residue:4d} {line[27:]}")

    if not output_lines:
        raise ValueError("Could not write a Protpardelle motif PDB from the selected target chains.")
    motif_path.parent.mkdir(parents=True, exist_ok=True)
    motif_path.write_text("\n".join(output_lines) + "\nEND\n")
    return lengths, chain_map, residue_map


def _map_protpardelle_hotspots(hotspots: str, residue_map: dict[tuple[str, int], tuple[str, int]]) -> str:
    mapped: list[str] = []
    for token in hotspots.split(",") if hotspots else []:
        if not token:
            continue
        try:
            chain = token[0]
            residue = int(token[1:])
        except ValueError:
            continue
        mapped_residue = residue_map.get((chain, residue))
        if mapped_residue:
            mapped.append(f"{mapped_residue[0]}{mapped_residue[1]}")
    return ",".join(mapped)


def _write_protpardelle_sampling_yaml(
    path: Path,
    motif_name: str,
    chain_lengths: dict[str, int],
    binder_length: str,
    hotspots: str,
    model_name: str,
    model_epoch: str,
    sampling_config: str,
    step_scale: float,
    schurn: int,
    crop_cond_start: float,
    translation: tuple[float, float, float],
) -> dict:
    lengths = parse_binder_lengths(binder_length)
    binder_range = [lengths[0], lengths[0]] if len(lengths) == 1 else [lengths[0], lengths[1]]
    target_contig = ";/;".join(f"{chain}1-{length}" for chain, length in chain_lengths.items())
    target_lengths = [[length, length] for length in chain_lengths.values()]
    motif_contig = f"{target_contig};/;{binder_range[0]}-{binder_range[1]}"
    total_lengths = [*target_lengths, binder_range]
    payload = {
        "search_space": {
            "models": [[model_name, model_epoch, sampling_config]],
            "step_scales": [float(step_scale)],
            "schurns": [int(schurn)],
            "crop_cond_starts": [float(crop_cond_start)],
            "translations": [[float(translation[0]), float(translation[1]), float(translation[2])]],
        },
        "motifs": [motif_name],
        "motif_contigs": [motif_contig],
        "total_lengths": [total_lengths],
        "hotspots": [hotspots],
        "ssadj": [None],
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            [
                "search_space:",
                "  models:",
                f"    - {json.dumps([model_name, model_epoch, sampling_config])}",
                f"  step_scales: {json.dumps([float(step_scale)])}",
                f"  schurns: {json.dumps([int(schurn)])}",
                f"  crop_cond_starts: {json.dumps([float(crop_cond_start)])}",
                f"  translations: {json.dumps([[float(translation[0]), float(translation[1]), float(translation[2])]])}",
                "",
                "motifs:",
                f"  - {motif_name}",
                "",
                "motif_contigs:",
                f"  - {motif_contig}",
                "",
                "total_lengths:",
                f"  - {json.dumps(total_lengths)}",
                "",
                "hotspots:",
                f'  - "{hotspots}"',
                "",
                "ssadj:",
                "  - null",
                "",
            ]
        )
    )
    return payload


def _protpardelle_native_pass(metrics: dict) -> bool | None:
    sample_rmsd = _float_or_none(metrics.get("ca_motif_sample_rmsd"))
    pred_rmsd = _float_or_none(metrics.get("allatom_motif_pred_rmsd"))
    scaffold_rmsd = _float_or_none(metrics.get("ca_scaffold_scrmsd"))
    if sample_rmsd is None or pred_rmsd is None or scaffold_rmsd is None:
        return None
    return sample_rmsd < 1.0 and pred_rmsd < 1.0 and scaffold_rmsd < 2.0


def _resolve_protpardelle_structure(run_dir: Path, row: dict[str, str]) -> Path | None:
    save_name = str(row.get("save_name") or "").strip()
    candidates: list[Path] = []
    if save_name:
        save_path = Path(save_name)
        candidates.append(save_path)
        if save_path.is_absolute() and str(save_path).startswith("/work/"):
            candidates.append(run_dir / save_path.relative_to("/work"))
        if save_path.suffix != ".pdb":
            candidates.append(save_path.with_suffix(".pdb"))
    for candidate in candidates:
        if candidate.exists():
            return candidate
    output_dir = run_dir / "artifacts" / "raw" / "protpardelle_1c" / "output"
    for path in sorted(output_dir.glob("**/esmfold/*.pdb")):
        if save_name and Path(save_name).stem not in path.stem:
            continue
        return path
    for path in sorted(output_dir.glob("**/*.pdb")):
        if "sample_" in path.name:
            continue
        return path
    return None


def _normalize_protpardelle_1c_candidates(run_dir: Path, params: dict, target_artifact: Path) -> list[dict]:
    output_dir = run_dir / "artifacts" / "raw" / "protpardelle_1c" / "output"
    metric_paths = sorted(output_dir.glob("**/esm_metrics.csv"))
    candidates: list[dict] = []
    seq_counts_by_structure: dict[str, int] = {}
    for metric_path in metric_paths:
        for row in _read_csv_rows(metric_path):
            structure_path = _resolve_protpardelle_structure(run_dir, row)
            if structure_path is None:
                continue
            structure_index = str(row.get("structure_index") or len(candidates) + 1)
            seq_counts_by_structure[structure_index] = seq_counts_by_structure.get(structure_index, 0) + 1
            design_number = int(float(structure_index)) + 1 if re.fullmatch(r"\d+(\.0)?", structure_index) else len(seq_counts_by_structure)
            seq_number = seq_counts_by_structure[structure_index]
            binder_chains, target_chains, chain_inference = _infer_generated_chain_roles(
                target_artifact,
                structure_path,
                params.get("target_chains", []),
            )
            sequences = _sequences_by_chain(structure_path)
            metrics = {key: _numeric_or_text(value) for key, value in row.items() if value not in {"", None}}
            native_pass = _protpardelle_native_pass(metrics)
            plddt = metrics.get("plddt")
            if isinstance(plddt, (int, float)) and plddt <= 1.0:
                plddt = plddt * 100.0
            metrics.update(
                {
                    "complex_refolding_backend": "protpardelle_1c_esmfold",
                    "result_kind": "native_pipeline",
                    "target_chain_inference": chain_inference,
                    "native_final_rank": len(candidates) + 1,
                    "native_design_index": design_number,
                    "native_sequence_index": seq_number,
                    "native_pass_filters": native_pass,
                    "pass_filters": native_pass,
                    "binder_plddt": plddt,
                    "ipae": metrics.get("pae"),
                    "binder_rmsd": metrics.get("ca_scaffold_scrmsd"),
                    "monomer_rmsd": metrics.get("ca_scaffold_scrmsd"),
                    "protpardelle_success_criteria": "ca_motif_sample_rmsd<1.0; allatom_motif_pred_rmsd<1.0; ca_scaffold_scrmsd<2.0",
                }
            )
            candidates.append(
                {
                    "candidate_id": f"protpardelle_1c_{design_number:05d}_mpnn_{seq_number:03d}",
                    "source_tool": "protpardelle_1c",
                    "stage": STAGE_COMPLEX_REFOLDING,
                    "target_pdb": _rel_path(run_dir, target_artifact),
                    "complex_pdb": _rel_path(run_dir, structure_path),
                    "binder_pdb": None,
                    "binder_sequence": "".join(sequences.get(chain, "") for chain in binder_chains) or None,
                    "target_chains": target_chains,
                    "binder_chains": binder_chains,
                    "hotspots": [token for token in params.get("hotspots", "").split(",") if token],
                    "binder_length": params.get("binder_length"),
                    "contig": params.get("motif_contig"),
                    "metrics": metrics,
                    "raw_metadata": {
                        "result_kind": "native_pipeline",
                        "metric_csv": _rel_path(run_dir, metric_path),
                        "input_target_chains": params.get("target_chains", []),
                        "target_chain_inference": chain_inference,
                        "output_kind": "protpardelle_1c_native",
                    },
                }
            )
    if candidates:
        def _rank_value(candidate: dict, key: str, default: float, invert: bool = False) -> float:
            value = (candidate.get("metrics") or {}).get(key)
            try:
                number = float(value)
            except (TypeError, ValueError):
                number = default
            return -number if invert else number

        candidates.sort(
            key=lambda candidate: (
                0 if (candidate.get("metrics") or {}).get("native_pass_filters") is True else 1,
                _rank_value(candidate, "ca_scaffold_scrmsd", float("inf")),
                _rank_value(candidate, "pae", float("inf")),
                _rank_value(candidate, "binder_plddt", float("-inf"), invert=True),
                _rank_value(candidate, "allatom_motif_pred_rmsd", float("inf")),
            )
        )
        for rank, candidate in enumerate(candidates, start=1):
            metrics = candidate.setdefault("metrics", {})
            metrics["native_final_rank"] = rank
    sample_candidates: list[dict] = []
    sample_paths = sorted(output_dir.glob("**/sample_*.pdb"))
    for index, structure_path in enumerate(sample_paths, start=1):
        binder_chains, target_chains, chain_inference = _infer_generated_chain_roles(
            target_artifact,
            structure_path,
            params.get("target_chains", []),
        )
        sample_candidates.append(
            {
                "candidate_id": f"protpardelle_1c_sample_{index:05d}",
                "source_tool": "protpardelle_1c",
                "stage": STAGE_GENERATION_BACKBONE,
                "target_pdb": _rel_path(run_dir, target_artifact),
                "complex_pdb": _rel_path(run_dir, structure_path),
                "binder_pdb": None,
                "binder_sequence": None,
                "target_chains": target_chains,
                "binder_chains": binder_chains,
                "hotspots": [token for token in params.get("hotspots", "").split(",") if token],
                "binder_length": params.get("binder_length"),
                "contig": params.get("motif_contig"),
                "metrics": {"result_kind": "generation_only", "target_chain_inference": chain_inference},
                "raw_metadata": {"result_kind": "generation_only", "target_chain_inference": chain_inference},
            }
        )
    candidates = _merge_candidate_pools(
        tool="protpardelle_1c",
        native_candidates=candidates,
        raw_candidates=sample_candidates,
        native_coverage="native Protpardelle ESMFold/MPNN metric rows",
        raw_coverage="raw Protpardelle scaffold samples before ESMFold/MPNN filtering",
    )
    return write_candidates(run_dir, "protpardelle_1c", candidates)


def run_protpardelle_1c(
    target_pdb: Path,
    target_chains: list[str],
    binder_length: str,
    hotspots: str = "",
    campaign_name: str = "",
    num_designs: int = 10,
    num_mpnn_seqs: int = 10,
    model_name: str = "cc83",
    model_epoch: str = "2616",
    sampling_config: str = "sampling_sidechain_conditional",
    step_scale: float = 1.2,
    schurn: int = 200,
    crop_cond_start: float = 0.0,
    batch_size: int = 1,
    seed: int = 7,
    gpu_device: object = "0",
) -> Path:
    if not target_chains:
        raise ValueError("At least one target chain is required.")
    if num_designs < 1:
        raise ValueError("Design attempts must be at least 1.")
    if num_mpnn_seqs < 0:
        raise ValueError("MPNN sequences per design must not be negative.")
    if batch_size < 1:
        raise ValueError("Batch size must be at least 1.")
    hotspots = normalize_hotspots(hotspots) if hotspots else ""
    parse_binder_lengths(binder_length)
    manifest = load_manifest("protpardelle_1c")
    params = {
        "target_chains": target_chains,
        "binder_length": binder_length,
        "hotspots": hotspots,
        "campaign_name": campaign_name.strip(),
        "num_designs": num_designs,
        "num_mpnn_seqs": num_mpnn_seqs,
        "model_name": model_name,
        "model_epoch": model_epoch,
        "sampling_config": sampling_config,
        "step_scale": step_scale,
        "schurn": schurn,
        "crop_cond_start": crop_cond_start,
        "batch_size": batch_size,
        "seed": seed,
        "gpu_device": normalize_gpu_device(gpu_device),
        "contract": "protpardelle_1c_native_pipeline",
    }
    job = create_job(DESIGN_GROUP, "design_campaign", "protpardelle_1c", {"target_pdb": str(target_pdb), "target_chains": target_chains}, params)
    if campaign_name.strip():
        update_status(job.run_dir, "queued", campaign_name=campaign_name.strip())
    raw_dir = job.run_dir / "artifacts" / "raw" / "protpardelle_1c"
    input_dir = raw_dir / "input"
    motif_dir = input_dir / "motifs"
    output_dir = raw_dir / "output"
    input_dir.mkdir(parents=True, exist_ok=True)
    target_artifact = _copy_target_for_design(job.run_dir, target_pdb, target_chains, "raw/protpardelle_1c/input")
    motif_name = "mn_app_target"
    chain_lengths, chain_map, residue_map = _write_protpardelle_target_motif(target_artifact, target_chains, motif_dir / f"{motif_name}.pdb")
    mapped_hotspots = _map_protpardelle_hotspots(hotspots, residue_map)
    sampling_payload = _write_protpardelle_sampling_yaml(
        input_dir / "protpardelle_sampling.yaml",
        motif_name,
        chain_lengths,
        binder_length,
        mapped_hotspots,
        model_name,
        model_epoch,
        sampling_config,
        step_scale,
        schurn,
        crop_cond_start,
        (0.0, 0.0, 0.0),
    )
    params.update(
        {
            "protpardelle_chain_map": chain_map,
            "protpardelle_hotspots": mapped_hotspots,
            "motif_contig": sampling_payload["motif_contigs"][0],
            "total_lengths": sampling_payload["total_lengths"][0],
        }
    )
    write_json(input_dir / "protpardelle_sampling_payload.json", sampling_payload)
    write_json(job.run_dir / "input.json", {"inputs": {"target_pdb": str(target_pdb), "target_chains": target_chains}, "params": params})
    command = [
        "docker",
        "run",
        "--rm",
        *docker_gpu_args(gpu_device),
        "-v",
        f"{job.run_dir}:/work",
        "-v",
        "/mnt/db/reference_files/protpardelle-1c:/ref/protpardelle-1c:ro",
        "-e",
        "PROTPARDELLE_OUTPUT_DIR=/work/artifacts/raw/protpardelle_1c/output",
        "-e",
        "PROTPARDELLE_MODEL_PARAMS=/ref/protpardelle-1c/model_params",
        "-e",
        "ESMFOLD_PATH=/ref/protpardelle-1c/model_params/ESMFold",
        "-e",
        "PROTEINMPNN_WEIGHTS=/ref/protpardelle-1c/model_params/ProteinMPNN/vanilla_model_weights",
        "-e",
        "LIGANDMPNN_WEIGHTS=/ref/protpardelle-1c/model_params/LigandMPNN",
        "-e",
        "FOLDSEEK_BIN=/opt/foldseek/bin/foldseek",
        "-w",
        "/work",
        manifest["image"],
        "python",
        "-m",
        "protpardelle.sample",
        "/work/artifacts/raw/protpardelle_1c/input/protpardelle_sampling.yaml",
        "--motif-dir",
        "/work/artifacts/raw/protpardelle_1c/input/motifs",
        "--num-samples",
        str(num_designs),
        "--num-mpnn-seqs",
        str(num_mpnn_seqs),
        "--batch-size",
        str(batch_size),
        "--seed",
        str(seed),
    ]
    rc = _run_shell_steps(job.run_dir, [{"name": "protpardelle-1c", "command": command}])
    candidates = _normalize_protpardelle_1c_candidates(job.run_dir, params, target_artifact) if rc == 0 else []
    _finish_design_job(
        job.run_dir,
        rc == 0 and bool(candidates),
        rc,
        candidates,
        [
            ("raw/protpardelle_1c/**/*.pdb", "pdb"),
            ("raw/protpardelle_1c/**/*.csv", "csv"),
            ("raw/protpardelle_1c/**/*.yaml", "yaml"),
            ("raw/protpardelle_1c/**/*.json", "json"),
            ("raw/protpardelle_1c/**/*.html", "report"),
        ],
    )
    return job.run_dir


def run_proteina_complexa(
    target_pdb: Path,
    target_chains: list[str],
    binder_length: str,
    hotspots: str = "",
    campaign_name: str = "",
    num_designs: int = 1,
    n_steps: int = 400,
    replicas: int = 2,
    seed: int = 0,
    batch_size: int = 1,
    gpu_device: object = "0",
    generator_only: bool = False,
    pipeline_task: str = "protein_binder",
    ligand: str = "",
    motif_atom_spec: str = "",
    search_strategy: str = "auto",
    guidance_weight: float = 1.0,
    sampling_noise: float = 0.1,
    deduplicate_sequences: bool = True,
    cpu_cores: int = 4,
) -> Path:
    parameters=dict(locals())
    from mn_molecular_design_studio import complexa_input as inputs
    import secrets
    p=inputs.prepare_request(parameters)
    resources=inputs.resource_paths(reference_root())
    ready=inputs.readiness(pipeline_task,reference_root(),generator_only=generator_only)
    if not ready["ready"]:
        raise ValueError("Proteina-Complexa reference files are missing: "+", ".join(ready["missing"]))
    p["seed"]=p["seed"] or secrets.randbelow(2147483646)+1
    p["run_name"]=re.sub(r"[^A-Za-z0-9_.-]+","_",campaign_name.strip()).strip("_") or "studio_complexa"
    p["target_pdb"]=str(target_pdb)
    p["gpu_device"]=normalize_gpu_device(gpu_device)
    p["engine_version"]=inputs.CODE_VERSION
    p["contract"]="complexa_native_generation" if generator_only else "complexa_native_pipeline"
    job=create_job(DESIGN_GROUP,"design_generation" if generator_only else "design_campaign",
                   "proteina_complexa_generator" if generator_only else "proteina_complexa",
                   {"target_pdb":str(target_pdb),"target_chains":p["target_chains"]},p)
    update_status(job.run_dir,"preparing",campaign_name=campaign_name or "Proteina-Complexa",current_phase="Preparing "+inputs.TASKS[pipeline_task])
    raw=job.run_dir/"artifacts/raw/proteina_complexa";raw.mkdir(parents=True,exist_ok=True)
    target,source=inputs.stage_target(p,raw/"input")
    target_input=_target_input_spec(target,p["target_chains"]) if pipeline_task=="protein_binder" else ""
    if pipeline_task=="protein_binder" and not target_input:
        raise ValueError("The selected protein target contains no readable residues.")
    p["target_input"]=target_input
    cache=resources["hf_cache"];cache.mkdir(parents=True,exist_ok=True)
    environment={"COMPLEXA_INIT":"1","DATA_PATH":"/workspace/protein-foundation-models/assets","CKPT_PATH":"/weights",
        "AF2_DIR":"/af2","RF3_CKPT_PATH":"/rf3.ckpt","RF3_EXEC_PATH":"/workspace/.venv/bin/rf3",
        "FOLDSEEK_EXEC":"/studio-tools/foldseek","MMSEQS_EXEC":"/studio-tools/mmseqs",
        "SC_EXEC":"/workspace/protein-foundation-models/src/proteinfoundation/result_analysis/sc",
        "LOCAL_CODE_PATH":"/workspace/protein-foundation-models","COMMUNITY_MODELS_PATH":"/workspace/protein-foundation-models/community_models",
        "LOCAL_CACHE_DIR":"/cache","CACHE_DIR":"/cache","HF_HOME":"/cache","HF_HUB_CACHE":"/cache",
        "JAX_COMPILATION_CACHE_DIR":"/cache/jax","JAX_ENABLE_COMPILATION_CACHE":"true","CUDA_CACHE_PATH":"/cache/cuda",
        "PYTORCH_CUDA_ALLOC_CONF":"expandable_segments:True"}
    spec={**p,"source_file":source.name,"base_config":inputs.TASK_CONFIGS[pipeline_task],"overrides":inputs.native_overrides(p,target_input),"environment":environment}
    write_json(raw/"input/target_config.json",spec)
    helper=Path(inputs.__file__).with_name("complexa_container.py")
    shutil.copy2(helper,raw/"run_complexa.py")
    mounts=["-v",f"{job.run_dir}:/work","-v",f"{resources['checkpoints']}:/weights:ro","-v",f"{cache}:/cache"]
    if not generator_only:
        mounts += ["-v",f"{resources['rf3'].resolve()}:/rf3.ckpt:ro","-v",f"{resources['af2']}:/af2:ro",
                   "-v",f"{resources['mpnn']}:/workspace/protein-foundation-models/community_models/LigandMPNN/model_params:ro",
                   "-v",f"{resources['tools']}:/studio-tools:ro"]
    command=["docker","run","--rm",*docker_gpu_args(gpu_device),"--entrypoint","/workspace/.venv/bin/python",
        *mounts,inputs.IMAGE,"/work/artifacts/raw/proteina_complexa/run_complexa.py","/work/artifacts/raw/proteina_complexa/input/target_config.json"]
    rc=_run_shell_steps(job.run_dir,[{"name":"proteina-complexa-"+pipeline_task,"command":command}])
    candidates=_normalize_proteina_complexa_candidates(job.run_dir,p,target) if rc==0 else []
    has_evaluation=bool(_proteina_complexa_result_rows(job.run_dir)[1])
    _finish_design_job(job.run_dir,rc==0 and bool(candidates) and (generator_only or has_evaluation),rc,candidates,
        [("raw/proteina_complexa/**/*."+extension,extension) for extension in ("pdb","cif","csv","json","yaml","log")])
    return job.run_dir


def run_proteina_complexa_generation(
    target_pdb: Path,
    target_chains: list[str],
    binder_length: str,
    hotspots: str = "",
    campaign_name: str = "",
    num_designs: int = 1,
    n_steps: int = 400,
    seed: int = 0,
    batch_size: int = 1,
    gpu_device: object = "0",
    pipeline_task: str = "protein_binder",
    ligand: str = "",
    motif_atom_spec: str = "",
    guidance_weight: float = 1.0,
    sampling_noise: float = 0.1,
    deduplicate_sequences: bool = True,
    cpu_cores: int = 4,
) -> Path:
    """Generate backbone and sequence, ready for independently chosen tools."""
    return run_proteina_complexa(**locals(),generator_only=True,search_strategy="single-pass",replicas=1)
