from __future__ import annotations

import csv
import math
import re
import shlex
import subprocess
from pathlib import Path
from typing import Any

from mn_molecular_design_studio.engine.workflow_engine.core.artifacts import Artifact, artifact_path
from mn_molecular_design_studio.engine.workflow_engine.core.candidates import (
    STAGE_COMPLEX_REFOLDING,
    STAGE_GENERATION_BACKBONE_SEQUENCE,
    write_candidates,
)
from mn_molecular_design_studio.engine.workflow_engine.core.gpu import docker_gpu_args, normalize_gpu_device
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import (
    create_job,
    finish_job,
    read_json,
    update_status,
    write_json,
)
from mn_molecular_design_studio.engine.workflow_engine.core.structures import filter_pdb_text
from mn_molecular_design_studio.engine.workflow_engine.core.local_worker import spawn_worker_for_run
from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import portable_path
from mn_molecular_design_studio.engine.workflow_engine.core.scheduler import apply_docker_cpu_limit, configured_cpu_slots
from mn_molecular_design_studio.engine.workflow_engine.runtime import reference_root
from mn_molecular_design_studio.engine.workflow_engine.workflows.design import (
    _infer_generated_chain_roles,
    _pdb_sequences_by_chain,
    _sequences_by_chain,
    normalize_hotspots,
    parse_binder_lengths,
)


DESIGN_GROUP = "design"
BINDCRAFT2_TOOL = "bindcraft2"
BINDCRAFT2_GENERATOR_TOOL = "bindcraft2_generator"
BINDCRAFT2_JOB_KIND = "bindcraft2_design"
BINDCRAFT2_IMAGE = "mn-bindcraft2:v1.0.4-cu13-pae-v1"
BINDCRAFT2_SOURCE_COMMIT = "d502f19a187a1bc8c05e833ca9cbab7bc1934fe7"
BINDCRAFT2_SOURCE_URL = "https://github.com/PacesaLab/BindCraft2.git"
BINDCRAFT2_REFERENCE_SUBDIR = "alphafold_models"
BINDCRAFT2_MODELS = tuple(
    [f"model_{index}_multimer_v3" for index in range(1, 6)] + ["model_1_ptm", "model_2_ptm"]
)
ALPHAFOLD_MINIMUM_BYTES = 100 * 1024 * 1024
BC2_AMINO_ACIDS = tuple("ACDEFGHIKLMNPQRSTVWY")
BC2_BINDER_FORMATS = (
    "binder",
    "large_binder",
    "peptide",
    "cyclic_peptide",
    "homo_oligomer",
    "multidomain",
    "VHH",
    "ARP",
    "scFv",
    "Fab",
)
BC2_CONFORMATIONAL_OBJECTIVES = ("induced_fit", "fold_switch")
SUPPORTED_MODALITIES = BC2_BINDER_FORMATS + BC2_CONFORMATIONAL_OBJECTIVES
BC2_PROPERTY_PRESETS = {
    "bigbang": "Start gradient design from the coordinates on hand.",
    "disulfide_staple": "Allow cysteine and require a geometrically detected disulfide.",
    "forced_targeting": "Concentrate the paratope on declared hotspots (requires hotspots).",
    "humanize": "Add human sequence preferences and an MHC-II anchor-score objective.",
    "initial_guess": "Use the trajectory pose as the starting guess for validation.",
    "mixed_topology": "Encourage mixed helix and beta-sheet topology.",
    "protease_stable": "Penalize protease motifs, exposed loops, and exposed termini.",
    "termini_accessible": "Orient both chain ends away from the target.",
    "termini_together": "Bring the N and C termini close together.",
}
BC2_SCAFFOLD_FORMATS = frozenset({"VHH", "scFv", "Fab", "ARP"})
BC2_SCAFFOLD_EDIT_DEFAULTS = {
    "VHH": "A27-33,A35*,A37*,A45*,A47*,A50*,A52-57(5-7),A59,A97-99*,A100-114(3-25),A115,A116",
    "scFv": "L30-32,L91-94,L96,H31-33,H50,H52-55(4-6),H57,H59,H99-106(5-18),H72*",
    "Fab": "L30-32,L91-94,L96,H31-33,H50,H52-55(4-6),H57,H59,H99-106(5-18),H72*",
    "ARP": "A16,A19,A20,A23,A44-48(5-10),A53,A56,A77-81(5-10),A86,A89,A110-114(5-10),A119,A122,A143-147(5-10),A152,A156",
}
BC2_PRESET_LENGTHS = {
    "binder": (60, 180),
    "large_binder": (250, 600),
    "peptide": (12, 25),
    "cyclic_peptide": (7, 16),
    "homo_oligomer": (40, 120),
    "multidomain": (120, 300),
}
WORKERS_PER_GPU_OPTIONS: tuple[str | int, ...] = ("auto", 1, 2, 3, 4)


def bindcraft2_cpu_options() -> tuple[list[int], int]:
    capacity = configured_cpu_slots()
    options = sorted({1, min(2, capacity), min(4, capacity), min(8, capacity)})
    return options, min(4, capacity)


def bindcraft2_reference_directory(root: Path | None = None) -> Path:
    return (root or reference_root()) / BINDCRAFT2_REFERENCE_SUBDIR


def missing_bindcraft2_parameters(
    root: Path | None = None,
    *,
    minimum_bytes: int = ALPHAFOLD_MINIMUM_BYTES,
) -> list[str]:
    """Return the required BC2 AF2 checkpoint files that are absent or incomplete."""
    directory = bindcraft2_reference_directory(root)
    missing: list[str] = []
    for model in BINDCRAFT2_MODELS:
        possible = (directory / f"params_{model}.npz", directory / "params" / f"params_{model}.npz")
        if not any(path.is_file() and path.stat().st_size >= minimum_bytes for path in possible):
            missing.append(f"params_{model}.npz")
    return missing


def bindcraft2_image_available(image: str = BINDCRAFT2_IMAGE) -> bool:
    try:
        result = subprocess.run(
            ["docker", "image", "inspect", image],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
            timeout=5,
        )
    except (OSError, subprocess.TimeoutExpired):
        return False
    return result.returncode == 0


def _validate_hotspots(hotspots: str, target_chains: list[str]) -> str:
    normalized = normalize_hotspots(hotspots) if hotspots else ""
    allowed = {str(chain).strip().upper() for chain in target_chains}
    invalid = []
    for token in normalized.split(","):
        match = re.fullmatch(r"(.+?)(\d+)", token)
        if token and (not match or match.group(1).upper() not in allowed):
            invalid.append(token)
    if invalid:
        raise ValueError(f"Hotspots must belong to selected target chains: {', '.join(invalid)}")
    return normalized


def normalize_disordered_target_sequence(text: str) -> str:
    """Accept one raw amino-acid sequence or one FASTA record for a BC2 IDR."""
    lines = [line.strip() for line in str(text or "").lstrip("\ufeff").splitlines() if line.strip()]
    headers = [line for line in lines if line.startswith(">")]
    if len(headers) > 1:
        raise ValueError("Upload or paste one FASTA target sequence at a time.")
    sequence = "".join(line for line in lines if not line.startswith(">"))
    sequence = re.sub(r"\s+", "", sequence).upper()
    if not sequence:
        raise ValueError("Enter a target amino-acid sequence or upload a FASTA file.")
    invalid = sorted(set(sequence) - set(BC2_AMINO_ACIDS))
    if invalid:
        raise ValueError("Use the 20 standard amino-acid letters (A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y).")
    if len(sequence) > 10000:
        raise ValueError("The disordered target sequence must be 10,000 residues or fewer.")
    return sequence


def parse_disordered_target_fasta(text: str) -> tuple[str, str]:
    """Parse a single named FASTA record for an intrinsically disordered target."""
    lines = [line.strip() for line in str(text or "").lstrip("\ufeff").splitlines() if line.strip()]
    headers = [line[1:].strip() for line in lines if line.startswith(">")]
    if len(headers) != 1:
        raise ValueError("Upload one FASTA record with a single >header.")
    target_name = headers[0]
    if not target_name:
        raise ValueError("The FASTA header must include a target name.")
    sequence = normalize_disordered_target_sequence(text)
    return target_name, sequence


def normalize_bindcraft2_modalities(modality: str | list[str] | tuple[str, ...]) -> list[str]:
    """Normalize a BC2 binder format and its optional conformational objective."""
    names = [str(value).strip() for value in ([modality] if isinstance(modality, str) else modality) if str(value).strip()]
    if not names:
        names = ["binder"]
    if len(set(names)) != len(names):
        raise ValueError("BindCraft 2 modalities cannot be repeated.")
    unknown = [name for name in names if name not in SUPPORTED_MODALITIES]
    if unknown:
        raise ValueError(f"Unknown BindCraft 2 modality: {', '.join(unknown)}")
    formats = [name for name in names if name in BC2_BINDER_FORMATS]
    objectives = [name for name in names if name in BC2_CONFORMATIONAL_OBJECTIVES]
    if len(formats) > 1:
        raise ValueError("Choose one BindCraft 2 binder format.")
    if len(objectives) > 1:
        raise ValueError("Choose at most one conformational objective.")
    if objectives and not formats:
        formats = ["binder"]
    return formats + objectives


def build_bindcraft2_settings(
    *,
    target_chains: list[str],
    binder_length: str | None,
    target_type: str = "structure",
    target_sequence: str = "",
    target_name: str = "app_target",
    target_filename: str | None = None,
    scaffold_edits: str | None = None,
    target_crop_lengths: tuple[int, int] | list[int] | None = None,
    target_crop_count: int | None = None,
    hotspots: str = "",
    coldspots: str = "",
    campaign_name: str = "",
    modality: str | list[str] | tuple[str, ...] = "binder",
    design_properties: list[str] | tuple[str, ...] = (),
    number_of_final_designs: int = 1,
    max_trajectories: int = 1,
    campaign_seed: int = 0,
    workers_per_gpu: str | int = "auto",
    redesign_interface: bool = False,
    binder_copies: int | None = None,
    number_of_domains: int | None = None,
    second_target_filename: str | None = None,
    second_target_chains: list[str] | None = None,
    second_target_hotspots: str = "",
    second_target_coldspots: str = "",
    second_target_goal: str = "bind",
    second_target_strength: float = 1.0,
    excluded_amino_acids: list[str] | tuple[str, ...] = (),
    amino_acid_bias: dict[str, float] | None = None,
    min_interface_iptm: float | None = None,
    max_interface_ipae: float | None = None,
    min_complex_ptm: float | None = None,
    min_unbound_binder_plddt: float | None = None,
    min_interface_residues: int | None = None,
    max_binder_rmsd: float | None = None,
    min_plddt_refine: float | None = None,
    min_plddt_mutate: float | None = None,
    quality_filters: dict[str, Any] | None = None,
    engine_overrides: dict[str, Any] | None = None,
    trajectory_only: bool = False,
) -> dict[str, Any]:
    target_type = str(target_type).strip().lower()
    if target_type not in {"structure", "sequence"}:
        raise ValueError("Target type must be structure or sequence.")
    chains = [str(chain).strip() for chain in target_chains if str(chain).strip()]
    normalized_sequence = normalize_disordered_target_sequence(target_sequence) if target_type == "sequence" else ""
    if target_type == "structure" and not chains:
        raise ValueError("At least one target chain is required.")
    if target_type == "sequence" and (chains or hotspots or coldspots):
        raise ValueError("Sequence targets do not use chains, hotspots, or coldspots.")
    lengths = parse_binder_lengths(binder_length) if str(binder_length or "").strip() else None
    modalities = normalize_bindcraft2_modalities(modality)
    objectives = [name for name in modalities if name in BC2_CONFORMATIONAL_OBJECTIVES]
    binder_format = next(name for name in modalities if name in BC2_BINDER_FORMATS)
    if objectives and second_target_filename:
        raise ValueError("BindCraft 2 conformational objectives require a single target.")
    if objectives and target_type == "sequence" and target_crop_count is not None and int(target_crop_count) != 1:
        raise ValueError("BindCraft 2 conformational objectives require exactly one target sequence crop.")
    if "fold_switch" in objectives and binder_format in BC2_SCAFFOLD_FORMATS:
        raise ValueError("The fold_switch objective cannot be combined with a fixed scaffold format.")
    properties = list(dict.fromkeys(str(value).strip() for value in design_properties if str(value).strip()))
    unknown_properties = [name for name in properties if name not in BC2_PROPERTY_PRESETS]
    if unknown_properties:
        raise ValueError(f"Unknown BindCraft 2 design property: {', '.join(unknown_properties)}")
    if "forced_targeting" in properties and not hotspots:
        raise ValueError("BindCraft 2 forced_targeting requires at least one target hotspot.")
    if binder_format in BC2_SCAFFOLD_FORMATS and "mixed_topology" in properties:
        raise ValueError("BindCraft 2 mixed_topology cannot be combined with a fixed scaffold modality.")
    if int(number_of_final_designs) < 1:
        raise ValueError("Desired final designs must be at least 1.")
    if int(max_trajectories) < 1:
        raise ValueError("Design attempts must be at least 1.")
    if int(campaign_seed) < 0:
        raise ValueError("Campaign seed must be zero or greater.")
    if workers_per_gpu != "auto":
        try:
            worker_count = int(workers_per_gpu)
        except (TypeError, ValueError) as exc:
            raise ValueError("Workers per GPU must be auto or a positive integer.") from exc
        if worker_count < 1:
            raise ValueError("Workers per GPU must be auto or a positive integer.")
        workers_per_gpu = worker_count

    hotspot_text = _validate_hotspots(hotspots, chains) if target_type == "structure" else ""
    coldspot_text = _validate_hotspots(coldspots, chains) if target_type == "structure" and coldspots else ""
    target_file = target_filename or ("target.fasta" if target_type == "sequence" else "target.pdb")
    target = {"name": str(target_name or "app_target").strip(), "target_path": target_file}
    if target_type == "structure":
        target["chains"] = ",".join(chains)
    if hotspot_text:
        target["hotspots"] = hotspot_text
    if coldspot_text:
        target["coldspots"] = coldspot_text
    targets = [target]
    second_chains = [str(chain).strip() for chain in (second_target_chains or []) if str(chain).strip()]
    if second_target_filename:
        if target_type != "structure":
            raise ValueError("A second structure target can only be added to a structure-target design.")
        if not second_chains:
            raise ValueError("Select at least one chain in the second target.")
        if second_target_goal not in {"bind", "avoid"}:
            raise ValueError("The second target goal must be bind or avoid.")
        second = {
            "name": "app_second_target",
            "target_path": second_target_filename,
            "chains": ",".join(second_chains),
            "weight": abs(float(second_target_strength)) * (1 if second_target_goal == "bind" else -1),
        }
        if second_target_hotspots:
            second["hotspots"] = _validate_hotspots(second_target_hotspots, second_chains)
        if second_target_coldspots:
            second["coldspots"] = _validate_hotspots(second_target_coldspots, second_chains)
        targets.append(second)
    name = str(campaign_name or "bindcraft2").strip()
    settings: dict[str, Any] = {
        "targets": targets,
        "modality": modalities[0] if len(modalities) == 1 else modalities,
        "campaign_name": name,
        "number_of_final_designs": int(number_of_final_designs),
        "max_trajectories": int(max_trajectories),
        "campaign_seed": int(campaign_seed),
        "workers_per_gpu": workers_per_gpu,
        "project_folder": "output",
        "resume": True,
        "trajectory_only": bool(trajectory_only),
        "save_design_trajectory": True,
        "save_design_animations": False,
        "save_design_frames": False,
    }
    if lengths is not None:
        settings["binder_lengths"] = lengths
    if scaffold_edits is not None:
        if binder_format not in BC2_SCAFFOLD_FORMATS:
            raise ValueError("Scaffold edit regions are only available for VHH, scFv, Fab, and ARP.")
        if not str(scaffold_edits).strip():
            raise ValueError("Scaffold edit regions cannot be empty.")
        settings["mutate_positions"] = str(scaffold_edits).strip()
    if target_type == "sequence":
        if target_crop_lengths is not None:
            crop_min, crop_max = (int(value) for value in target_crop_lengths)
            if crop_min < 1 or crop_max < crop_min:
                raise ValueError("Disordered-target crop lengths must be positive and ordered.")
            settings["crop_fasta_sequence"] = [crop_min, crop_max]
        if target_crop_count is not None:
            if int(target_crop_count) < 1:
                raise ValueError("The number of disordered-target crops must be at least 1.")
            settings["idr_crop_count"] = int(target_crop_count)
    settings.update({name: True for name in properties})
    if redesign_interface:
        settings["redesign_interface"] = True
    if binder_copies is not None:
        settings["copies"] = int(binder_copies)
    if number_of_domains is not None:
        settings["n_domains"] = int(number_of_domains)
    aa_bias = {str(name).upper(): 0.0 for name in excluded_amino_acids}
    aa_bias.update({str(name).upper(): float(value) for name, value in (amino_acid_bias or {}).items()})
    if aa_bias:
        settings["aa_bias"] = aa_bias
    for name, value in (
        ("min_iptm_final", min_interface_iptm),
        ("max_ipae_final", max_interface_ipae),
        ("min_ptm_final", min_complex_ptm),
        ("min_monomer_plddt_final", min_unbound_binder_plddt),
    ):
        if value is not None:
            settings[name] = float(value)
    if min_interface_residues is not None:
        settings.setdefault("filters", {})["Interface_Residues"] = {
            "threshold": int(min_interface_residues),
            "higher": True,
        }
    if max_binder_rmsd is not None:
        settings.setdefault("filters", {})["Binder_RMSD"] = {"threshold": float(max_binder_rmsd), "higher": False}
    if min_plddt_refine is not None:
        settings["min_plddt_refine"] = float(min_plddt_refine)
    if min_plddt_mutate is not None:
        settings["min_plddt_mutate"] = float(min_plddt_mutate)
    if quality_filters:
        filter_overrides = quality_filters.get("filters")
        if isinstance(filter_overrides, dict):
            settings["filters"] = {**settings.get("filters", {}), **filter_overrides}
        for setting_name, value in quality_filters.items():
            if setting_name != "filters":
                settings[setting_name] = value
    if engine_overrides:
        if not isinstance(engine_overrides, dict):
            raise ValueError("Advanced BindCraft 2 settings must be a JSON object.")
        settings.update(engine_overrides)
    # The selected tile owns whether this is a generator or a full campaign;
    # advanced engine overrides must not silently change the workflow type.
    settings["trajectory_only"] = bool(trajectory_only)
    if trajectory_only:
        # The standalone generator contract must not be changed into a full
        # ProteinMPNN campaign by an advanced override.
        settings.update(
            trajectory_only=True,
            number_of_final_designs=1,
            max_trajectories=int(max_trajectories),
            save_design_trajectory=True,
            save_failed_trajectories=True,
        )
    return settings


def bounded_bindcraft2_settings(settings: dict[str, Any], campaign_seed: int) -> dict[str, Any]:
    """Make one isolated, reproducible native trajectory attempt."""
    if type(campaign_seed) is not int or campaign_seed < 0:
        raise ValueError("A bounded BindCraft2 attempt needs a nonnegative campaign seed.")
    bounded = dict(settings)
    bounded.update(
        number_of_final_designs=1,
        max_trajectories=1,
        campaign_seed=campaign_seed,
        resume=True,
        save_design_trajectory=True,
        save_failed_trajectories=True,
    )
    return bounded


def _portable_command(command: list[str], reference_dir: Path) -> list[str]:
    source = str(reference_dir.expanduser().resolve())
    replacement = portable_path(reference_dir)
    return [token.replace(source, replacement) for token in command]


def build_bindcraft2_docker_command(
    *,
    run_dir: Path,
    reference_dir: Path,
    gpu_device: object = "0",
    image: str = BINDCRAFT2_IMAGE,
) -> list[str]:
    gpu = normalize_gpu_device(gpu_device)
    reference = reference_dir.expanduser().resolve()
    workdir = "/work/artifacts/raw/bindcraft2"
    command = [
        "docker",
        "run",
        "--rm",
        *docker_gpu_args(gpu),
        "--shm-size=32G",
        "-e",
        f"BINDCRAFT_AF2_PARAMS=/ref/{BINDCRAFT2_REFERENCE_SUBDIR}",
        "-v",
        f"{run_dir.expanduser().resolve()}:/work:rw",
        "-v",
        f"{reference}:/ref:ro",
        "-w",
        workdir,
        image,
        "design",
        "settings.json",
    ]
    return apply_docker_cpu_limit(command, run_dir)


def _prepare_bindcraft2_run(
    *,
    target_pdb: Path | None = None,
    target_chains: list[str] | None = None,
    binder_length: str | None,
    target_type: str = "structure",
    target_fasta: Path | None = None,
    scaffold_edits: str | None = None,
    target_crop_lengths: tuple[int, int] | list[int] | None = None,
    target_crop_count: int | None = None,
    hotspots: str = "",
    coldspots: str = "",
    campaign_name: str = "",
    modality: str | list[str] | tuple[str, ...] = "binder",
    design_properties: list[str] | tuple[str, ...] = (),
    number_of_final_designs: int = 1,
    max_trajectories: int = 1,
    campaign_seed: int = 0,
    workers_per_gpu: str | int = "auto",
    redesign_interface: bool = False,
    binder_copies: int | None = None,
    number_of_domains: int | None = None,
    second_target_pdb: Path | None = None,
    second_target_chains: list[str] | None = None,
    second_target_hotspots: str = "",
    second_target_coldspots: str = "",
    second_target_goal: str = "bind",
    second_target_strength: float = 1.0,
    excluded_amino_acids: list[str] | tuple[str, ...] = (),
    amino_acid_bias: dict[str, float] | None = None,
    min_interface_iptm: float | None = None,
    max_interface_ipae: float | None = None,
    min_complex_ptm: float | None = None,
    min_unbound_binder_plddt: float | None = None,
    min_interface_residues: int | None = None,
    max_binder_rmsd: float | None = None,
    min_plddt_refine: float | None = None,
    min_plddt_mutate: float | None = None,
    quality_filters: dict[str, Any] | None = None,
    engine_overrides: dict[str, Any] | None = None,
    trajectory_only: bool = False,
    distributed_campaign: dict[str, Any] | None = None,
    cpu_cores: int = 4,
    gpu_device: object = "0",
    dispatch_worker: bool,
    parent_run_dir: Path | None = None,
) -> Path:
    distributed_config = dict(distributed_campaign or {})
    distributed_enabled = bool(distributed_config.get("enabled"))
    if distributed_enabled:
        if trajectory_only:
            raise ValueError("Distributed acceptance campaigns require the full BindCraft2 design tile.")
        from mn_molecular_design_studio.engine.workflow_engine.workflows.distributed_campaign import validate_campaign_config

        distributed_config = validate_campaign_config(distributed_config)
        distributed_config.setdefault("gpu_selection", normalize_gpu_device(gpu_device))
    target_type = str(target_type).strip().lower()
    target_chains = [str(chain).strip() for chain in (target_chains or []) if str(chain).strip()]
    if target_type not in {"structure", "sequence"}:
        raise ValueError("Target type must be structure or sequence.")
    target_name = "app_target"
    normalized_sequence = ""
    if target_type == "structure":
        if target_pdb is None:
            raise ValueError("Choose a PDB or mmCIF structure target.")
        target_pdb = Path(target_pdb).expanduser().resolve()
        if not target_pdb.is_file():
            raise FileNotFoundError(f"Target structure not found: {target_pdb}")
        if target_pdb.suffix.lower() not in {".pdb", ".ent", ".cif", ".mmcif"}:
            raise ValueError("Structure targets must be PDB or mmCIF files.")
        if not target_chains:
            raise ValueError("Select at least one target chain.")
        available_target_chains = _sequences_by_chain(target_pdb)
        missing_chains = [chain for chain in target_chains if chain not in available_target_chains]
        if missing_chains:
            raise ValueError(f"Target structure does not contain selected chain(s): {', '.join(missing_chains)}")
        target_filename = "target.cif" if target_pdb.suffix.lower() in {".cif", ".mmcif"} else "target.pdb"
    else:
        if target_pdb is not None or target_chains or hotspots or coldspots:
            raise ValueError("Sequence targets do not use structure files, chains, hotspots, or coldspots.")
        if target_fasta is None:
            raise ValueError("Upload a named FASTA file for a disordered sequence target.")
        target_fasta = Path(target_fasta).expanduser().resolve()
        if not target_fasta.is_file():
            raise FileNotFoundError(f"Target FASTA not found: {target_fasta}")
        try:
            target_name, normalized_sequence = parse_disordered_target_fasta(
                target_fasta.read_text(encoding="utf-8-sig")
            )
        except UnicodeDecodeError as exc:
            raise ValueError("The target FASTA file must be UTF-8 text.") from exc
        target_filename = "target.fasta"
    second_target_pdb = Path(second_target_pdb).expanduser().resolve() if second_target_pdb else None
    second_target_chains = [str(chain).strip() for chain in (second_target_chains or []) if str(chain).strip()]
    if second_target_pdb is not None:
        if target_type != "structure":
            raise ValueError("A second structure target can only be added to a structure-target design.")
        if not second_target_pdb.is_file():
            raise FileNotFoundError(f"Second target structure not found: {second_target_pdb}")
        if not second_target_chains:
            raise ValueError("Select at least one chain in the second target.")
        missing_second_chains = [chain for chain in second_target_chains if chain not in _sequences_by_chain(second_target_pdb)]
        if missing_second_chains:
            raise ValueError(f"Second target does not contain selected chain(s): {', '.join(missing_second_chains)}")
    if int(cpu_cores) < 1:
        raise ValueError("CPU cores must be at least 1.")
    settings = build_bindcraft2_settings(
        target_chains=target_chains,
        binder_length=binder_length,
        target_type=target_type,
        target_sequence=normalized_sequence,
        target_name=target_name,
        target_filename=target_filename,
        scaffold_edits=scaffold_edits,
        target_crop_lengths=target_crop_lengths,
        target_crop_count=target_crop_count,
        hotspots=hotspots,
        coldspots=coldspots,
        campaign_name=campaign_name,
        modality=modality,
        design_properties=design_properties,
        number_of_final_designs=number_of_final_designs,
        max_trajectories=max_trajectories,
        campaign_seed=campaign_seed,
        workers_per_gpu=workers_per_gpu,
        redesign_interface=redesign_interface,
        binder_copies=binder_copies,
        number_of_domains=number_of_domains,
        second_target_filename="second_target.pdb" if second_target_pdb else None,
        second_target_chains=second_target_chains,
        second_target_hotspots=second_target_hotspots,
        second_target_coldspots=second_target_coldspots,
        second_target_goal=second_target_goal,
        second_target_strength=second_target_strength,
        excluded_amino_acids=excluded_amino_acids,
        amino_acid_bias=amino_acid_bias,
        min_interface_iptm=min_interface_iptm,
        max_interface_ipae=max_interface_ipae,
        min_complex_ptm=min_complex_ptm,
        min_unbound_binder_plddt=min_unbound_binder_plddt,
        min_interface_residues=min_interface_residues,
        max_binder_rmsd=max_binder_rmsd,
        min_plddt_refine=min_plddt_refine,
        min_plddt_mutate=min_plddt_mutate,
        quality_filters=quality_filters,
        engine_overrides=engine_overrides,
        trajectory_only=trajectory_only,
    )
    parameters = missing_bindcraft2_parameters()
    if parameters:
        raise FileNotFoundError(
            "BindCraft 2 AlphaFold parameters are missing or incomplete under "
            f"{bindcraft2_reference_directory()}: {', '.join(parameters)}"
        )
    if not bindcraft2_image_available():
        raise RuntimeError(
            f"Docker image {BINDCRAFT2_IMAGE} is missing. Build it with "
            "bash scripts/build_bindcraft2_pae_image.sh from the Design Studio repository root."
        )

    selected_gpu = "none" if distributed_enabled else normalize_gpu_device(gpu_device)
    params = {
        "target_type": target_type,
        "target_file": target_filename,
        "target_name": target_name if target_type == "sequence" else "",
        "target_chains": target_chains,
        "binder_length": binder_length,
        "target_crop_lengths": list(target_crop_lengths) if target_crop_lengths is not None else None,
        "target_crop_count": target_crop_count,
        "scaffold_edits": scaffold_edits,
        "engine_overrides": engine_overrides or {},
        "hotspots": settings["targets"][0].get("hotspots", ""),
        "coldspots": settings["targets"][0].get("coldspots", ""),
        "campaign_name": settings["campaign_name"],
        "modality": settings["modality"],
        "design_properties": [name for name in BC2_PROPERTY_PRESETS if settings.get(name)],
        "redesign_interface": bool(redesign_interface),
        "binder_copies": binder_copies,
        "number_of_domains": number_of_domains,
        "second_target_chains": second_target_chains,
        "second_target_hotspots": second_target_hotspots,
        "second_target_coldspots": second_target_coldspots,
        "second_target_goal": second_target_goal,
        "second_target_strength": float(second_target_strength),
        "excluded_amino_acids": list(excluded_amino_acids),
        "amino_acid_bias": amino_acid_bias or {},
        "number_of_final_designs": int(number_of_final_designs),
        "max_trajectories": int(max_trajectories),
        "campaign_seed": int(campaign_seed),
        "workers_per_gpu": settings["workers_per_gpu"],
        "cpu_cores": int(cpu_cores),
        "gpu_device": selected_gpu,
        "docker_image": BINDCRAFT2_IMAGE,
        "bindcraft2_source_commit": BINDCRAFT2_SOURCE_COMMIT,
        "bindcraft2_source_url": BINDCRAFT2_SOURCE_URL,
        "reference_key": f"reference:///{BINDCRAFT2_REFERENCE_SUBDIR}",
        "quality_filters": quality_filters or {},
        "trajectory_only": bool(trajectory_only),
        "distributed_campaign": distributed_config if distributed_enabled else {},
    }
    job = create_job(
        DESIGN_GROUP,
        "bindcraft2_campaign" if distributed_enabled else "bindcraft2_generator" if trajectory_only else "bindcraft2_design",
        "bindcraft2_campaign" if distributed_enabled else BINDCRAFT2_GENERATOR_TOOL if trajectory_only else BINDCRAFT2_TOOL,
        inputs={
            "target_type": target_type,
            "target_structure": f"artifacts/raw/bindcraft2/{target_filename}" if target_type == "structure" else None,
            "target_pdb": f"artifacts/raw/bindcraft2/{target_filename}" if target_filename == "target.pdb" else None,
            "target_fasta": "artifacts/raw/bindcraft2/target.fasta" if target_type == "sequence" else None,
            "target_chains": params["target_chains"],
            "second_target_pdb": "artifacts/raw/bindcraft2/second_target.pdb" if second_target_pdb else None,
        },
        params=params,
    )
    raw_dir = artifact_path(job.run_dir, "raw", "bindcraft2")
    raw_dir.mkdir(parents=True, exist_ok=True)
    staged_target = raw_dir / target_filename
    if target_type == "sequence":
        staged_target.write_text(f">{target_name}\n{normalized_sequence}\n", encoding="utf-8")
    elif target_pdb.suffix.lower() in {".cif", ".mmcif"}:
        staged_target.write_bytes(target_pdb.read_bytes())
    else:
        target_text = target_pdb.read_text(errors="ignore")
        staged_target.write_text(filter_pdb_text(target_text, keep_chains=set(params["target_chains"])))
    if second_target_pdb is not None:
        second_target_text = second_target_pdb.read_text(errors="ignore")
        (raw_dir / "second_target.pdb").write_text(
            filter_pdb_text(second_target_text, keep_chains=set(second_target_chains))
        )
    settings_path = raw_dir / "settings.json"
    write_json(settings_path, settings)
    if dispatch_worker:
        write_json(
            job.run_dir / "worker_request.json",
            {"kind": "bindcraft2_campaign_controller" if distributed_enabled else BINDCRAFT2_JOB_KIND, "kwargs": {}, "path_kwargs": []},
        )
        write_json(
            job.run_dir / "command.json",
            {
                "mode": "local_worker",
                "command": ["python", "-m", "mn_molecular_design_studio.engine.workflow_engine.core.local_worker", "--run-dir", "runs:///design/" + job.run_dir.name],
            },
        )
    else:
        parent_metadata = read_json(Path(parent_run_dir) / "metadata.json") if parent_run_dir else {}
        parent_allocation = parent_metadata.get("resource_allocation") if isinstance(parent_metadata.get("resource_allocation"), dict) else {}
        child_metadata = read_json(job.run_dir / "metadata.json")
        child_metadata["resource_allocation"] = {
            **parent_allocation,
            "cpu_cores": int(parent_allocation.get("cpu_cores") or cpu_cores),
            "gpu_device": selected_gpu,
            "scheduler": str(parent_allocation.get("scheduler") or "design-campaign-parent"),
        }
        write_json(job.run_dir / "metadata.json", child_metadata)
        write_json(
            job.run_dir / "command.json",
            {
                "mode": "embedded_campaign",
                "parent_run": portable_path(Path(parent_run_dir)) if parent_run_dir else "",
            },
        )
    if campaign_name.strip():
        update_status(job.run_dir, "queued", campaign_name=campaign_name.strip())
    if dispatch_worker:
        spawn_worker_for_run(job.run_dir)
    return job.run_dir


def enqueue_bindcraft2_design(
    *,
    target_pdb: Path | None = None,
    target_chains: list[str] | None = None,
    binder_length: str | None = None,
    target_type: str = "structure",
    target_fasta: Path | None = None,
    hotspots: str = "",
    coldspots: str = "",
    campaign_name: str = "",
    modality: str | list[str] | tuple[str, ...] = "binder",
    design_properties: list[str] | tuple[str, ...] = (),
    number_of_final_designs: int = 1,
    max_trajectories: int = 1,
    campaign_seed: int = 0,
    workers_per_gpu: str | int = "auto",
    redesign_interface: bool = False,
    binder_copies: int | None = None,
    number_of_domains: int | None = None,
    second_target_pdb: Path | None = None,
    second_target_chains: list[str] | None = None,
    second_target_hotspots: str = "",
    second_target_coldspots: str = "",
    second_target_goal: str = "bind",
    second_target_strength: float = 1.0,
    excluded_amino_acids: list[str] | tuple[str, ...] = (),
    amino_acid_bias: dict[str, float] | None = None,
    min_interface_iptm: float | None = None,
    max_interface_ipae: float | None = None,
    min_complex_ptm: float | None = None,
    min_unbound_binder_plddt: float | None = None,
    min_interface_residues: int | None = None,
    max_binder_rmsd: float | None = None,
    min_plddt_refine: float | None = None,
    min_plddt_mutate: float | None = None,
    quality_filters: dict[str, Any] | None = None,
    scaffold_edits: str | None = None,
    target_crop_lengths: tuple[int, int] | list[int] | None = None,
    target_crop_count: int | None = None,
    engine_overrides: dict[str, Any] | None = None,
    trajectory_only: bool = False,
    distributed_campaign: dict[str, Any] | None = None,
    cpu_cores: int = 4,
    gpu_device: object = "0",
) -> Path:
    campaign_config = dict(distributed_campaign or {})
    selected_gpu = normalize_gpu_device(gpu_device)
    if campaign_config.get("enabled"):
        campaign_config.setdefault("gpu_selection", selected_gpu)
        gpu_device = "none"
    return _prepare_bindcraft2_run(
        target_pdb=target_pdb,
        target_chains=target_chains,
        binder_length=binder_length,
        target_type=target_type,
        target_fasta=target_fasta,
        hotspots=hotspots,
        coldspots=coldspots,
        campaign_name=campaign_name,
        modality=modality,
        design_properties=design_properties,
        number_of_final_designs=number_of_final_designs,
        max_trajectories=max_trajectories,
        campaign_seed=campaign_seed,
        workers_per_gpu=workers_per_gpu,
        redesign_interface=redesign_interface,
        binder_copies=binder_copies,
        number_of_domains=number_of_domains,
        second_target_pdb=second_target_pdb,
        second_target_chains=second_target_chains,
        second_target_hotspots=second_target_hotspots,
        second_target_coldspots=second_target_coldspots,
        second_target_goal=second_target_goal,
        second_target_strength=second_target_strength,
        excluded_amino_acids=excluded_amino_acids,
        amino_acid_bias=amino_acid_bias,
        min_interface_iptm=min_interface_iptm,
        max_interface_ipae=max_interface_ipae,
        min_complex_ptm=min_complex_ptm,
        min_unbound_binder_plddt=min_unbound_binder_plddt,
        min_interface_residues=min_interface_residues,
        max_binder_rmsd=max_binder_rmsd,
        min_plddt_refine=min_plddt_refine,
        min_plddt_mutate=min_plddt_mutate,
        quality_filters=quality_filters,
        scaffold_edits=scaffold_edits,
        target_crop_lengths=target_crop_lengths,
        target_crop_count=target_crop_count,
        engine_overrides=engine_overrides,
        trajectory_only=trajectory_only,
        distributed_campaign=campaign_config,
        cpu_cores=cpu_cores,
        gpu_device=gpu_device,
        dispatch_worker=True,
    )


def run_bindcraft2_campaign(
    *,
    parent_run_dir: Path,
    target_pdb: Path,
    target_chains: list[str],
    binder_length: str | None,
    hotspots: str = "",
    campaign_name: str = "",
    modality: str | list[str] | tuple[str, ...] = "binder",
    design_properties: list[str] | tuple[str, ...] = (),
    number_of_final_designs: int = 1,
    max_trajectories: int = 1,
    campaign_seed: int = 0,
    workers_per_gpu: str | int = "auto",
    cpu_cores: int = 1,
    gpu_device: object = "0",
    scaffold_edits: str | None = None,
    engine_overrides: dict[str, Any] | None = None,
) -> Path:
    """Run BC2 as an inline child of an already scheduled vanilla campaign."""
    child_run = _prepare_bindcraft2_run(
        target_pdb=target_pdb,
        target_chains=target_chains,
        binder_length=binder_length,
        hotspots=hotspots,
        campaign_name=campaign_name,
        modality=modality,
        design_properties=design_properties,
        scaffold_edits=scaffold_edits,
        engine_overrides=engine_overrides,
        number_of_final_designs=number_of_final_designs,
        max_trajectories=max_trajectories,
        campaign_seed=campaign_seed,
        workers_per_gpu=workers_per_gpu,
        cpu_cores=cpu_cores,
        gpu_device=gpu_device,
        dispatch_worker=False,
        parent_run_dir=parent_run_dir,
    )
    run_bindcraft2_job(child_run)
    return child_run


def _coerce_metric(value: str) -> Any:
    value = value.strip()
    if not value:
        return None
    try:
        number = float(value)
    except ValueError:
        return value
    if not math.isfinite(number):
        return value
    return int(number) if number.is_integer() else number


def _normalize_bindcraft2_trajectory_candidates(
    run_dir: Path,
    settings: dict[str, Any],
    raw_dir: Path,
    target_artifact: Path,
) -> list[dict[str, Any]]:
    output_dir = raw_dir / "output"
    trajectory_root = output_dir / "1_Trajectories"
    trajectory_csv = trajectory_root / "!_Trajectories.csv"
    rows: dict[str, dict[str, str]] = {}
    if trajectory_csv.is_file():
        with trajectory_csv.open(newline="", encoding="utf-8-sig") as handle:
            rows = {
                str(row.get("design") or "").strip(): row
                for row in csv.DictReader(handle)
                if str(row.get("design") or "").strip()
            }

    target_type = str(settings.get("target_type") or "structure")
    target_sequence = str(settings.get("target_sequence") or "").strip().upper()
    if target_type == "sequence" and not target_sequence and target_artifact.is_file():
        try:
            _, target_sequence = parse_disordered_target_fasta(target_artifact.read_text(encoding="utf-8-sig"))
        except (OSError, UnicodeDecodeError, ValueError):
            target_sequence = ""
    input_target_chains = [str(chain) for chain in settings.get("target_chains") or []]
    requested_hotspots = [token for token in str(settings.get("hotspots") or "").split(",") if token]
    # BindCraft2 prefixes each structure with the trajectory name, e.g.
    # ``<design>_trajectory.cif``. Rows with a non-empty ``terminated`` stage
    # failed before the normal campaign's ProteinMPNN handoff.
    trajectory_files = sorted(
        path
        for path in trajectory_root.glob("**/*_trajectory*.cif")
        if not path.stem.endswith("_monomer")
    )
    trajectory_groups: dict[str, list[Path]] = {}
    for path in trajectory_files:
        trajectory_groups.setdefault(path.parent.name, []).append(path)
    target_reference_sequences = (
        {
            chain: sequence
            for chain, sequence in _sequences_by_chain(target_artifact).items()
            if not input_target_chains or chain in input_target_chains
        }
        if target_type == "structure" and target_artifact.is_file()
        else {}
    )
    candidates: list[dict[str, Any]] = []
    for index, (trajectory_name, structures) in enumerate(sorted(trajectory_groups.items()), start=1):
        structure_path = next(
            (path for path in structures if path.name == f"{trajectory_name}_trajectory.cif"),
            None,
        )
        if structure_path is None:
            # Multi-target and cropped campaigns add a state suffix. Use the
            # state matching the input target as the downstream MPNN template.
            for possible_path in structures:
                output_sequences = _sequences_by_chain(possible_path)
                if target_type == "sequence" and target_sequence:
                    matches_target = any(sequence and sequence in target_sequence for sequence in output_sequences.values())
                elif target_reference_sequences:
                    target_sequences = set(target_reference_sequences.values())
                    matches_target = any(sequence in target_sequences for sequence in output_sequences.values())
                else:
                    matches_target = False
                if matches_target:
                    structure_path = possible_path
                    break
        if structure_path is None:
            structure_path = structures[0]
        sequences = _sequences_by_chain(structure_path)
        row = rows.get(trajectory_name)
        if row is None or "terminated" not in row or str(row.get("terminated") or "").strip():
            continue
        binder_chains: list[str]
        target_chains: list[str]
        chain_inference: str
        if target_type == "sequence" and target_sequence:
            target_chains = [
                chain for chain, sequence in sequences.items()
                if sequence and sequence in target_sequence
            ]
            binder_chains = [chain for chain in sequences if chain not in target_chains]
            chain_inference = "fasta_sequence_match" if target_chains else "trajectory_roles"
        else:
            binder_chains, target_chains, chain_inference = _infer_generated_chain_roles(
                target_artifact,
                structure_path,
                input_target_chains,
            )
        binder_sequence = str(row.get("Binder_Sequence") or "").strip().upper() or "".join(
            sequences.get(chain, "") for chain in binder_chains
        )
        metrics = {
            key: value
            for key, raw_value in row.items()
            if key not in {"design", "Binder_Sequence"}
            and (value := _coerce_metric(str(raw_value or ""))) is not None
        }
        metrics.update(
            {
                "result_kind": "generated",
                "pass_filters": None,
                "pre_mpnn_filter_pass": True,
                "reached_mpnn_stage": True,
                "target_chain_inference": chain_inference,
            }
        )
        structure_ref = portable_path(structure_path, run_dir=run_dir)
        candidates.append(
            {
                "candidate_id": f"bindcraft2_trajectory_{index:05d}",
                "stage": STAGE_GENERATION_BACKBONE_SEQUENCE,
                "source_tool": BINDCRAFT2_GENERATOR_TOOL,
                "target_pdb": portable_path(target_artifact, run_dir=run_dir)
                if target_type == "structure" and target_artifact.is_file() else None,
                "complex_pdb": structure_ref,
                "binder_pdb": None,
                "binder_sequence": binder_sequence or None,
                "target_chains": target_chains or input_target_chains,
                "binder_chains": binder_chains,
                "hotspots": requested_hotspots,
                "binder_length": str(row.get("length") or settings.get("binder_length") or ""),
                "metrics": metrics,
                "raw_metadata": {
                    "result_kind": "trajectory_only",
                    "generator_only": True,
                    "pre_mpnn_filter_pass": True,
                    "reached_mpnn_stage": True,
                    "mpnn_eligibility_gate": "passed_bindcraft2_gradient_design_stages",
                    "trajectory_name": trajectory_name,
                    "trajectory_csv": trajectory_csv.relative_to(run_dir).as_posix() if trajectory_csv.is_file() else None,
                    "trajectory_structure": structure_ref,
                    "input_target_chains": input_target_chains,
                    "target_type": target_type,
                    "target_fasta": portable_path(target_artifact, run_dir=run_dir)
                    if target_type == "sequence" and target_artifact.is_file() else None,
                    "target_chain_inference": chain_inference,
                    "sequence_source": "gradient-optimized trajectory before ProteinMPNN sequence redesign",
                },
            }
        )
    return write_candidates(run_dir, BINDCRAFT2_GENERATOR_TOOL, candidates)


def normalize_bindcraft2_candidates(run_dir: Path, settings: dict[str, Any]) -> list[dict[str, Any]]:
    raw_dir = run_dir / "artifacts" / "raw" / "bindcraft2"
    target_type = str(settings.get("target_type") or "structure")
    target_filename = Path(
        str(settings.get("target_file") or ("target.fasta" if target_type == "sequence" else "target.pdb"))
    ).name
    target_artifact = raw_dir / target_filename
    output_dir = raw_dir / "output"
    if bool(settings.get("trajectory_only")):
        return _normalize_bindcraft2_trajectory_candidates(run_dir, settings, raw_dir, target_artifact)
    ranked_csv = output_dir / "3_Ranked" / "!_Ranked.csv"
    if not ranked_csv.is_file():
        return write_candidates(run_dir, BINDCRAFT2_TOOL, [])

    # BindCraft 2 can include designs in its ranked table that failed the
    # native refolding filters. Keep the ranked order, but carry the actual
    # pass/fail result from the refolded table instead of assuming every row
    # in 3_Ranked passed.
    refolded_csv = output_dir / "2_Refolded" / "!_Refolded.csv"
    refolded_by_sequence: dict[str, list[dict[str, str]]] = {}
    if refolded_csv.is_file():
        with refolded_csv.open(newline="", encoding="utf-8-sig") as handle:
            for row in csv.DictReader(handle):
                sequence = "".join(str(row.get("Binder_Sequence") or "").split()).upper()
                if sequence:
                    refolded_by_sequence.setdefault(sequence, []).append(row)

    target_sequence = str(settings.get("target_sequence") or "").strip().upper()
    if target_type == "sequence" and not target_sequence and target_artifact.is_file():
        _, target_sequence = parse_disordered_target_fasta(target_artifact.read_text(encoding="utf-8"))
    input_target_chains = [str(chain) for chain in settings.get("target_chains") or []]
    requested_hotspots = [token for token in str(settings.get("hotspots") or "").split(",") if token]
    requested_length = str(settings.get("binder_length") or "")
    candidates: list[dict[str, Any]] = []
    with ranked_csv.open(newline="", encoding="utf-8-sig") as handle:
        for index, row in enumerate(csv.DictReader(handle), start=1):
            design_name = Path(str(row.get("design") or "")).name
            if not design_name:
                continue
            ranked_sequence = "".join(str(row.get("Binder_Sequence") or "").split()).upper()
            refolded_matches = refolded_by_sequence.get(ranked_sequence, [])
            passing_matches = [
                match for match in refolded_matches
                if str(match.get("outcome") or "").strip().casefold() in {"passed", "pass"}
                and not str(match.get("failed_filters") or "").strip()
            ]
            passed_filters = bool(passing_matches)
            refolded_outcome = "passed" if passed_filters else (
                str(refolded_matches[0].get("outcome") or "").strip().casefold() if refolded_matches else ""
            )
            failed_filters = ""
            if not passed_filters and refolded_matches:
                failed_filters = str(refolded_matches[0].get("failed_filters") or "").strip()
            structures = sorted(
                path
                for path in (output_dir / "3_Ranked").glob(f"{design_name}*.cif")
                if not path.stem.endswith("_monomer")
            )
            if not structures:
                continue
            complex_path = next((path for path in structures if path.stem == design_name), structures[0])
            sequences = _sequences_by_chain(complex_path)
            if target_type == "sequence":
                designed_sequence = str(row.get("Binder_Sequence") or "").strip().upper()
                inferred_target_chains = [
                    chain for chain, chain_sequence in sequences.items()
                    if chain_sequence and chain_sequence in target_sequence and chain_sequence != designed_sequence
                ]
                # IDR campaigns can export only sampled target windows. Several
                # windows may be joined, so the native binder sequence is a
                # second authoritative way to identify the remaining target.
                if not inferred_target_chains and designed_sequence:
                    designed_chains = [chain for chain, sequence in sequences.items() if sequence == designed_sequence]
                    if designed_chains:
                        inferred_target_chains = [chain for chain in sequences if chain not in designed_chains]
                if inferred_target_chains:
                    binder_chains = [chain for chain in sequences if chain not in inferred_target_chains]
                    chain_inference = (
                        "fasta_sequence_match"
                        if any(sequence == target_sequence for sequence in sequences.values())
                        else "fasta_crop_or_binder_sequence_match"
                    )
                else:
                    binder_chains, inferred_target_chains, chain_inference = _infer_generated_chain_roles(
                        target_artifact,
                        complex_path,
                        input_target_chains,
                    )
            else:
                binder_chains, inferred_target_chains, chain_inference = _infer_generated_chain_roles(
                    target_artifact,
                    complex_path,
                    input_target_chains,
                )
            sequence = str(row.get("Binder_Sequence") or "").strip() or "".join(
                sequences.get(chain, "") for chain in binder_chains
            )
            metrics = {
                key: value
                for key, raw_value in row.items()
                if key not in {"design", "Binder_Sequence"} and (value := _coerce_metric(str(raw_value or ""))) is not None
            }
            rank = _coerce_metric(str(row.get("rank") or index))
            metrics.update(
                {
                    "bindcraft2_rank": rank,
                    "native_pass_filters": passed_filters,
                    "result_kind": "native_pipeline",
                    "target_chain_inference": chain_inference,
                }
            )
            # Store stable app-managed URIs here; the candidate writer and result.json both retain
            # them without recording this machine's active run-root path.
            complex_ref = portable_path(complex_path, run_dir=run_dir)
            pae_sidecar = complex_path.with_name(f"{complex_path.stem}.pae.json")
            pae_ref = portable_path(pae_sidecar, run_dir=run_dir) if pae_sidecar.is_file() else None
            target_ref = (
                portable_path(target_artifact, run_dir=run_dir)
                if target_type == "structure" and target_artifact.is_file()
                else None
            )
            candidates.append(
                {
                    "candidate_id": f"bindcraft2_{index:05d}",
                    "stage": STAGE_COMPLEX_REFOLDING,
                    "source_tool": BINDCRAFT2_TOOL,
                    "target_pdb": target_ref,
                    "complex_pdb": complex_ref,
                    "binder_pdb": None,
                    "binder_sequence": sequence or None,
                    "target_chains": inferred_target_chains or input_target_chains,
                    "binder_chains": binder_chains,
                    "hotspots": requested_hotspots,
                    "binder_length": str(row.get("length") or requested_length),
                    "metrics": metrics,
                    "raw_metadata": {
                        "bindcraft2_design": design_name,
                        "bindcraft2_ranked_csv": ranked_csv.relative_to(run_dir).as_posix(),
                        "bindcraft2_output_structure": complex_ref,
                        "pae_path": pae_ref,
                        "bindcraft2_pae_matrix": pae_ref,
                        "bindcraft2_row": row,
                        "bindcraft2_refolded_outcome": refolded_outcome or None,
                        "bindcraft2_failed_filters": failed_filters or None,
                        "bindcraft2_refolded_sequence_matches": len(refolded_matches),
                        "input_target_chains": input_target_chains,
                        "target_type": target_type,
                        "target_fasta": portable_path(target_artifact) if target_type == "sequence" and target_artifact.is_file() else None,
                        "target_chain_inference": chain_inference,
                    },
                }
            )
    return write_candidates(run_dir, BINDCRAFT2_TOOL, candidates)


def _bindcraft2_artifacts(run_dir: Path, *, trajectory_only: bool = False) -> list[dict[str, Any]]:
    raw = run_dir / "artifacts" / "raw" / "bindcraft2"
    candidates = run_dir / "artifacts" / "normalized_candidates"
    selected: list[tuple[Path, str]] = []
    patterns = (
        (raw / "output" / "3_Ranked" / "!_Ranked.csv", "bindcraft2_ranked_table"),
        (raw / "output" / "2_Refolded" / "!_Refolded.csv", "bindcraft2_refolded_table"),
        (raw / "output" / "1_Trajectories" / "!_Trajectories.csv", "bindcraft2_trajectory_table"),
        (raw / "output" / "summary.csv", "bindcraft2_summary"),
        (raw / "output" / "campaign_metadata.json", "bindcraft2_campaign_metadata"),
        (raw / "output" / ".campaign_state.json", "bindcraft2_native_campaign_state"),
        (raw / "settings.json", "bindcraft2_input_settings"),
        (candidates / "candidates.jsonl", "normalized_candidates"),
        (candidates / "campaign_result.json", "campaign_result"),
    )
    for path, artifact_type in patterns:
        if path.is_file():
            selected.append((path, artifact_type))
    for target_filename in ("target.pdb", "target.cif", "target.fasta"):
        target_path = raw / target_filename
        if target_path.is_file():
            selected.append((target_path, "bindcraft2_input_target" if target_filename != "target.fasta" else "bindcraft2_input_sequence"))
    metadata_files = sorted((raw / "output").glob("campaign_metadata*.json"))
    selected.extend((path, "bindcraft2_campaign_metadata") for path in metadata_files if path.is_file())
    ranked_dir = raw / "output" / "3_Ranked"
    if ranked_dir.is_dir():
        selected.extend((path, "designed_complex_mmcif") for path in sorted(ranked_dir.glob("*.cif")) if not path.stem.endswith("_monomer"))
    pae_directories = (ranked_dir, raw / "output" / "2_Refolded" / "Complexes")
    for pae_directory in pae_directories:
        if pae_directory.is_dir():
            selected.extend((path, "pairwise_pae_matrix") for path in sorted(pae_directory.glob("*.pae.json")))
    if trajectory_only:
        trajectory_root = raw / "output" / "1_Trajectories"
        if trajectory_root.is_dir():
            selected.extend(
                (path, "bindcraft2_generated_complex")
                for path in sorted(trajectory_root.glob("**/*_trajectory*.cif"))
                if not path.stem.endswith("_monomer")
            )
    return [Artifact(path.stem, path, artifact_type).to_json(run_dir) for path, artifact_type in selected]


def read_bindcraft2_campaign_state(run_dir: Path) -> dict[str, Any] | None:
    """Read and validate native progress; the normalized full-campaign result omits this count."""
    path = Path(run_dir) / "artifacts" / "raw" / "bindcraft2" / "output" / ".campaign_state.json"
    if not path.is_file():
        return None
    import json

    try:
        state = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise RuntimeError(f"BindCraft2 campaign-state file is unreadable: {path}") from exc
    if not isinstance(state, dict):
        raise RuntimeError("BindCraft2 campaign state must be a JSON object.")
    for key in ("trajectories", "accepted"):
        value = state.get(key)
        if type(value) is not int or value < 0:
            raise RuntimeError(f"BindCraft2 campaign state has an invalid {key} count.")
    rejections = state.get("rejections", {})
    if not isinstance(rejections, dict):
        raise RuntimeError("BindCraft2 campaign-state rejections must be an object.")
    return {"trajectories": state["trajectories"], "accepted": state["accepted"], "rejections": rejections}


def run_bindcraft2_job(run_dir: Path) -> Path:
    run_dir = Path(run_dir).expanduser().resolve()
    payload = read_json(run_dir / "input.json")
    params = payload.get("params") if isinstance(payload.get("params"), dict) else {}
    inputs = payload.get("inputs") if isinstance(payload.get("inputs"), dict) else {}
    reference_dir = reference_root()
    missing = missing_bindcraft2_parameters(reference_dir)
    if missing:
        raise FileNotFoundError(
            f"BindCraft 2 AlphaFold parameters are missing under {bindcraft2_reference_directory(reference_dir)}: "
            + ", ".join(missing)
        )
    image = str(params.get("docker_image") or BINDCRAFT2_IMAGE)
    command = build_bindcraft2_docker_command(
        run_dir=run_dir,
        reference_dir=reference_dir,
        gpu_device=params.get("gpu_device", "0"),
        image=image,
    )
    command_for_record = _portable_command(command, reference_dir)
    write_json(
        run_dir / "command.json",
        {"mode": "docker", "steps": [{"name": BINDCRAFT2_TOOL, "command": command_for_record}]},
    )
    update_status(run_dir, "running")
    stdout_path = run_dir / "stdout.log"
    stderr_path = run_dir / "stderr.log"
    with stdout_path.open("a", encoding="utf-8") as stdout, stderr_path.open("a", encoding="utf-8") as stderr:
        stdout.write(f"$ {shlex.join(command_for_record)}\n")
        stdout.flush()
        process = subprocess.run(command, cwd=run_dir, stdout=stdout, stderr=stderr, check=False)

    candidates = normalize_bindcraft2_candidates(run_dir, params)
    trajectory_only = bool(params.get("trajectory_only"))
    artifacts = _bindcraft2_artifacts(run_dir, trajectory_only=trajectory_only)
    trajectory_csv = run_dir / "artifacts" / "raw" / "bindcraft2" / "output" / "1_Trajectories" / "!_Trajectories.csv"
    trajectory_rows: list[dict[str, str]] = []
    if trajectory_csv.is_file():
        with trajectory_csv.open(newline="", encoding="utf-8-sig") as handle:
            trajectory_rows = list(csv.DictReader(handle))
    native_state = read_bindcraft2_campaign_state(run_dir)
    finish_job(
        run_dir,
        process.returncode == 0,
        {
            "outputs": {"artifacts": artifacts, "candidates": candidates},
            "metrics": {
                "return_code": int(process.returncode),
                "candidate_count": len(candidates),
                "artifact_count": len(artifacts),
                "accepted_design_count": 0 if trajectory_only else sum(
                    candidate.get("metrics", {}).get("native_pass_filters") is True for candidate in candidates
                ),
                "trajectory_only": trajectory_only,
                "generated_trajectory_count": len(candidates) if trajectory_only else 0,
                "mpnn_eligible_trajectory_count": len(candidates) if trajectory_only else 0,
                "trajectory_attempt_count": (
                    native_state["trajectories"] if native_state is not None
                    else len(trajectory_rows) if trajectory_only else 0
                ),
                "native_accepted_count": native_state["accepted"] if native_state is not None else None,
                "native_rejections": native_state["rejections"] if native_state is not None else {},
                "native_campaign_state_verified": native_state is not None,
                "bindcraft2_image": image,
                "bindcraft2_source_commit": str(params.get("bindcraft2_source_commit") or BINDCRAFT2_SOURCE_COMMIT),
                "bindcraft2_source_url": str(params.get("bindcraft2_source_url") or BINDCRAFT2_SOURCE_URL),
                "bindcraft2_workers_per_gpu": params.get("workers_per_gpu", "auto"),
                "bindcraft2_cpu_cores": params.get("cpu_cores", 4),
            },
            "downstream_artifacts": {
                "target_structure": inputs.get("target_structure"),
                "target_pdb": inputs.get("target_pdb"),
                "target_fasta": inputs.get("target_fasta"),
                "settings_json": "artifacts/raw/bindcraft2/settings.json",
                "campaign_output": "artifacts/raw/bindcraft2/output",
                "candidates_jsonl": "artifacts/normalized_candidates/candidates.jsonl",
                "campaign_result": "artifacts/normalized_candidates/campaign_result.json",
            },
        },
    )
    return run_dir
