"""Input validation and file preparation rules for the pinned RFantibody engine."""
from __future__ import annotations

import re
import json
import math
from pathlib import Path

from Bio.SeqUtils import seq1
from mn_molecular_design_studio.mimosa_input import normalize_mimosa_source_positions

EXAMPLE_INPUTS = Path(__file__).resolve().parent / "resources" / "rfantibody"
UPSTREAM_COMMIT = "3e5a5bb1cb98befde9fc049ae3ff58cabb491f67"
CODE_VERSION = f"RFantibody · {UPSTREAM_COMMIT[:8]}"
IMAGE = "mn-rfantibody:3e5a5bb-cu128"
MIMOSA_IMAGE = "mn-rfantibody:mimosa-local"
MIMOSA_ENGINE_ROOT = Path(__file__).resolve().parents[2] / "tools-to-implement" / "RFantibody_mimosa"
TASKS = {"antibody": "Antibody", "nanobody": "Nanobody"}
FRAMEWORKS = {
    "hu-4D5-8_Fv": "Humanized antibody · hu-4D5-8 Fv",
    "h-NbBCII10": "Nanobody · h-NbBCII10",
    "custom": "Upload a custom Chothia-numbered framework",
}
CDR_REGIONS = ("H1", "H2", "H3", "L1", "L2", "L3")
NANOBODY_REGIONS = CDR_REGIONS[:3]
AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
PDB_CHAIN_IDS = "ABCDEFGIJKLMNOPQRSTUVWXYZabcdefghijkmnopqrstuvwxyz0123456789"


def normalize_regions(task: str, regions: object) -> list[str]:
    allowed = NANOBODY_REGIONS if task == "nanobody" else CDR_REGIONS
    if task not in TASKS:
        raise ValueError("Choose Antibody or Nanobody.")
    if isinstance(regions, str):
        regions = [item.strip() for item in regions.split(",") if item.strip()]
    if not isinstance(regions, (list, tuple, set)):
        raise ValueError("Choose at least one CDR region to design.")
    chosen = list(dict.fromkeys(str(item).strip().upper() for item in regions))
    invalid = [item for item in chosen if item not in allowed]
    if invalid:
        raise ValueError(f"These CDR regions are not available for {TASKS[task].lower()} design: {', '.join(invalid)}.")
    if not chosen:
        raise ValueError("Choose at least one CDR region to design.")
    return [item for item in allowed if item in chosen]


def normalize_loop_lengths(task: str, regions: list[str], values: dict | None = None) -> dict[str, str]:
    values = values if isinstance(values, dict) else {}
    normalized: dict[str, str] = {}
    for region in regions:
        value = str(values.get(region, "auto") or "auto").strip().lower()
        if value == "auto":
            normalized[region] = ""
            continue
        match = re.fullmatch(r"(\d+)(?:\s*[-–]\s*(\d+))?", value)
        if not match:
            raise ValueError(f"{region} length must be auto, a residue count, or a range such as 10-15.")
        start, end = int(match.group(1)), int(match.group(2) or match.group(1))
        if not 1 <= start <= end <= 40:
            raise ValueError(f"{region} length must be between 1 and 40 residues, with the minimum no greater than the maximum.")
        normalized[region] = str(start) if start == end else f"{start}-{end}"
    return normalized


def build_design_loops(task: str, regions: object, loop_lengths: dict | None = None) -> str:
    selected = normalize_regions(task, regions)
    lengths = normalize_loop_lengths(task, selected, loop_lengths)
    return ",".join(f"{region}:{lengths[region]}" for region in selected)


def normalize_global_bias(value: object) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    parsed: dict[str, float] = {}
    for token in text.split(","):
        match = re.fullmatch(r"\s*([A-Za-z])\s*:\s*(-?(?:\d+(?:\.\d*)?|\.\d+))\s*", token)
        if not match:
            raise ValueError("Amino-acid bias must use comma-separated entries such as W:3.0,P:-1.0.")
        amino_acid, value_text = match.groups()
        amino_acid = amino_acid.upper()
        score = float(value_text)
        if amino_acid not in AMINO_ACIDS or not math.isfinite(score) or abs(score) > 10:
            raise ValueError("Global amino-acid bias uses standard one-letter codes and scores from -10 to 10.")
        if amino_acid in parsed:
            raise ValueError(f"Amino-acid bias repeats {amino_acid}.")
        parsed[amino_acid] = score
    return ",".join(f"{amino_acid}:{score:g}" for amino_acid, score in parsed.items())


def normalize_per_residue_bias(value: object) -> str:
    text = str(value or "{}").strip() or "{}"
    try:
        parsed = json.loads(text)
    except json.JSONDecodeError as exc:
        raise ValueError("Per-residue bias must be valid JSON.") from exc
    if not isinstance(parsed, dict):
        raise ValueError('Per-residue bias JSON must be an object such as {"H":{"31":{"W":2.0}}}.')
    normalized: dict[str, dict[str, dict[str, float]]] = {}
    for chain, positions in parsed.items():
        if chain not in {"H", "L"} or not isinstance(positions, dict):
            raise ValueError("Per-residue bias keys must be H or L antibody chain IDs with residue objects.")
        normalized[chain] = {}
        for residue, biases in positions.items():
            if not str(residue).isdigit() or int(residue) < 1 or not isinstance(biases, dict) or not biases:
                raise ValueError("Per-residue bias positions must be positive residue numbers with amino-acid score objects.")
            normalized[chain][str(int(residue))] = {}
            for amino_acid, score in biases.items():
                amino_acid = str(amino_acid).upper()
                if amino_acid not in AMINO_ACIDS or isinstance(score, bool) or not isinstance(score, (int, float)) or not math.isfinite(float(score)) or abs(float(score)) > 10:
                    raise ValueError("Per-residue bias uses standard one-letter amino-acid codes and scores from -10 to 10.")
                normalized[chain][str(int(residue))][amino_acid] = float(score)
    return json.dumps(normalized, separators=(",", ":"), sort_keys=True)


def normalize_hotspots(hotspots: object, target_chains: list[str]) -> list[str]:
    if isinstance(hotspots, (list, tuple, set)):
        text = ",".join(str(item) for item in hotspots)
    else:
        text = str(hotspots or "")
    allowed = {str(chain).strip() for chain in target_chains}
    result = []
    for token in re.split(r"[,;\s]+", text.strip()):
        if not token:
            continue
        match = re.fullmatch(r"([A-Za-z0-9]+?)(-?\d+)", token)
        if not match:
            raise ValueError(f"Hotspot {token!r} must include a chain ID and residue number, for example B146.")
        chain, number = match.groups()
        if chain not in allowed:
            raise ValueError(f"Hotspot {token!r} is not in the selected antigen chains.")
        value = f"{chain}{int(number)}"
        if value not in result:
            result.append(value)
    return result


def validate_request(parameters: dict, *, require_files: bool = True) -> dict:
    """Validate user/API inputs without creating jobs or running the model."""
    task = str(parameters.get("task", "antibody"))
    if task not in TASKS:
        raise ValueError("Choose Antibody or Nanobody.")
    framework = str(parameters.get("framework", "h-NbBCII10" if task == "nanobody" else "hu-4D5-8_Fv"))
    if framework not in FRAMEWORKS:
        raise ValueError("Choose a supported framework or upload a custom framework.")
    if task == "nanobody" and framework == "hu-4D5-8_Fv":
        raise ValueError("Choose a nanobody framework for Nanobody mode.")
    if task == "antibody" and framework == "h-NbBCII10":
        raise ValueError("Choose an antibody framework for Antibody mode.")
    target_value = parameters.get("target_pdb")
    target = Path(target_value).expanduser() if target_value else None
    if require_files and (target is None or not target.is_file()):
        raise ValueError("Load a target PDB or mmCIF structure first.")
    target_chains = [str(item).strip() for item in parameters.get("target_chains") or [] if str(item).strip()]
    if len(set(target_chains)) != len(target_chains):
        raise ValueError("Selected antigen chains must be unique.")
    if require_files and not target_chains:
        raise ValueError("Select at least one antigen chain.")
    hotspots = normalize_hotspots(parameters.get("hotspots", ""), target_chains)
    if require_files and target is not None:
        from mn_molecular_design_studio.design_validation import protein_chains, read_structure

        if target.suffix.lower() not in {".pdb", ".ent", ".cif", ".mmcif"}:
            raise ValueError("Target structures must be PDB or mmCIF files.")
        try:
            source_chains = protein_chains(read_structure(target))
        except Exception as exc:
            raise ValueError(f"Could not read target structure {target.name}: {exc}") from exc
        missing = [chain for chain in target_chains if chain not in source_chains]
        if missing:
            raise ValueError(f"Selected antigen chain(s) are missing or contain no standard protein residues: {', '.join(missing)}.")
        for hotspot in hotspots:
            match = re.fullmatch(r"([A-Za-z0-9]+?)(-?\d+)", hotspot)
            chain, number = match.groups()
            if not any(int(residue.id[1]) == int(number) and not str(residue.id[2]).strip() for residue in source_chains[chain]):
                raise ValueError(f"Hotspot {hotspot} does not match a residue without an insertion code in the target structure.")
    regions = normalize_regions(task, parameters.get("regions", NANOBODY_REGIONS if task == "nanobody" else CDR_REGIONS))
    loops = normalize_loop_lengths(task, regions, parameters.get("loop_lengths"))
    enabled_value = parameters.get("mimosa_enabled")
    mimosa_enabled = bool(enabled_value) if enabled_value is not None else False
    mimosa_condition = str(
        parameters.get("mimosa_condition", "sequence_fixed" if mimosa_enabled else "unconstrained")
    ).strip().lower()
    if mimosa_condition not in {"geometry_only", "sequence_fixed", "unconstrained"}:
        raise ValueError("MIMOSA condition must be geometry_only, sequence_fixed, or unconstrained.")
    if enabled_value is None:
        mimosa_enabled = mimosa_condition != "unconstrained"
    elif mimosa_enabled and mimosa_condition == "unconstrained":
        raise ValueError("Choose a MIMOSA condition when motif guidance is enabled.")
    elif not mimosa_enabled:
        # The explicit tile toggle is authoritative. Standard RFantibody calls
        # still carry the MIMOSA condition default, but must remain unconstrained.
        mimosa_condition = "unconstrained"
    mimosa_complex = Path(parameters["mimosa_complex"]).expanduser() if parameters.get("mimosa_complex") else None
    mimosa_binder_chain = str(parameters.get("mimosa_binder_chain", "") or "").strip()
    mimosa_complex_target_chain = str(parameters.get("mimosa_complex_target_chain", "") or "").strip()
    mimosa_target_chain = str(parameters.get("mimosa_target_chain", "") or "").strip()
    mimosa_cdr = str(parameters.get("mimosa_cdr", "") or "").strip().upper()
    mimosa_design_loop = str(parameters.get("mimosa_design_loop", "") or "").strip().upper()
    mimosa_source_positions = normalize_mimosa_source_positions(
        parameters.get("mimosa_source_positions", "")
    )
    if mimosa_enabled:
        if not mimosa_source_positions:
            raise ValueError("Select at least one contacting amino acid from the source partner chain to form the MIMOSA motif.")
        if require_files and (mimosa_complex is None or not mimosa_complex.is_file()):
            raise ValueError("MIMOSA requires a cognate binder-target complex PDB or mmCIF.")
        if mimosa_complex is not None and mimosa_complex.suffix.lower() not in {".pdb", ".ent", ".cif", ".mmcif"}:
            raise ValueError("MIMOSA source complexes must be PDB or mmCIF files.")
        if not all((mimosa_binder_chain, mimosa_complex_target_chain, mimosa_target_chain)):
            raise ValueError("Select the source binder, source target, and design-target chains for MIMOSA.")
        if mimosa_target_chain not in target_chains:
            raise ValueError("The MIMOSA design-target chain must be selected as an antigen chain.")
        if mimosa_cdr not in {"CDR1", "CDR2", "CDR3"}:
            raise ValueError("Choose source antibody CDR1, CDR2, or CDR3 for MIMOSA.")
        if mimosa_design_loop == "ALL_CDRS":
            if not {"H1", "H2", "H3"}.issubset(regions):
                raise ValueError("Scattered MIMOSA assignment requires H1, H2, and H3 in Regions to design.")
        elif mimosa_design_loop not in regions:
            raise ValueError("The MIMOSA design loop must also be selected in Regions to design.")
    try:
        mimosa_contact_cutoff = float(parameters.get("mimosa_contact_cutoff", 8.0))
        mimosa_guidance_scale = float(parameters.get("mimosa_guidance_scale", 0.3))
        mimosa_num_particles = int(parameters.get("mimosa_num_particles", 4))
        mimosa_fk_lambda = float(parameters.get("mimosa_fk_lambda", 4.0))
        mimosa_resampling_interval = int(parameters.get("mimosa_resampling_interval", 3))
    except (TypeError, ValueError) as exc:
        raise ValueError("MIMOSA numeric settings must contain valid numbers.") from exc
    if mimosa_contact_cutoff <= 0 or not math.isfinite(mimosa_contact_cutoff):
        raise ValueError("MIMOSA contact cutoff must be a positive finite value.")
    if mimosa_guidance_scale < 0 or not math.isfinite(mimosa_guidance_scale):
        raise ValueError("MIMOSA guidance scale must be a non-negative finite value.")
    if not 1 <= mimosa_num_particles <= 64:
        raise ValueError("MIMOSA particle count must be between 1 and 64.")
    if mimosa_fk_lambda < 0 or not math.isfinite(mimosa_fk_lambda):
        raise ValueError("MIMOSA Feynman–Kac lambda must be a non-negative finite value.")
    if mimosa_resampling_interval < 1:
        raise ValueError("MIMOSA resampling interval must be at least one diffusion step.")
    try:
        number_of_designs = int(parameters.get("number_of_designs", 1))
        number_of_sequences = int(parameters.get("number_of_sequences", 1))
        rf2_repeats = int(parameters.get("rf2_repeats", 5))
        num_recycles = int(parameters.get("num_recycles", 3))
        num_connections = int(parameters.get("num_connections", 48))
        cpu_cores = int(parameters.get("cpu_cores", 4))
        mpnn_temperature = float(parameters.get("mpnn_temperature", 0.1))
        augment_eps = float(parameters.get("augment_eps", 0.05))
        hotspot_show_prop = float(parameters.get("hotspot_show_prop", 0.1))
        rf2_seed = int(parameters["rf2_seed"]) if parameters.get("rf2_seed") is not None else None
        diffusion_deterministic = bool(parameters.get("diffusion_deterministic", False))
    except (TypeError, ValueError) as exc:
        raise ValueError("Numeric RFantibody settings must contain valid numbers.") from exc
    for value, label, low, high in (
        (number_of_designs, "Number of designs", 1, 100000),
        (number_of_sequences, "Sequences per structure", 1, 100),
        (rf2_repeats, "Number of RF2 models", 1, 20),
        (num_recycles, "RF2 recycles", 0, 20),
        (num_connections, "ProteinMPNN neighbors", 1, 100),
        (cpu_cores, "CPU cores", 1, 256),
    ):
        if not low <= value <= high:
            raise ValueError(f"{label} must be between {low} and {high}.")
    if rf2_seed is not None and not 0 <= rf2_seed <= 2_147_483_647:
        raise ValueError("RF2 seed must be an integer between 0 and 2147483647.")
    if not 0.01 <= mpnn_temperature <= 2.0:
        raise ValueError("ProteinMPNN temperature must be between 0.01 and 2.0.")
    if not 0.0 <= augment_eps <= 1.0:
        raise ValueError("ProteinMPNN augmentation epsilon must be between 0 and 1.")
    if not 0.0 <= hotspot_show_prop <= 1.0:
        raise ValueError("Hotspot visibility must be between 0 and 1.")
    omit_aas = str(parameters.get("omit_aas", "CX")).strip().upper()
    if set(omit_aas) - (AMINO_ACIDS | {"X"}):
        raise ValueError("Omit amino acids using standard one-letter codes, optionally including X.")
    bias_aas = normalize_global_bias(parameters.get("bias_aas", ""))
    bias_by_res_json = normalize_per_residue_bias(parameters.get("bias_by_res_json", "{}"))
    if bias_aas or json.loads(bias_by_res_json):
        raise ValueError(
            "The installed RFantibody ProteinMPNN interface does not support global or per-residue amino-acid bias yet."
        )
    use_abmpnn = bool(parameters.get("use_abmpnn", False))
    abmpnn_checkpoint = None
    if use_abmpnn:
        from mn_molecular_design_studio.mpnn_input import checkpoint

        abmpnn_checkpoint = checkpoint("ab_mpnn")
        if require_files and not abmpnn_checkpoint.is_file():
            raise ValueError(f"AbMPNN weights are not installed at {abmpnn_checkpoint}.")
    framework_path = Path(parameters["framework_pdb"]).expanduser() if parameters.get("framework_pdb") else None
    if framework == "custom" and require_files and (framework_path is None or not framework_path.is_file()):
        raise ValueError("Upload a custom Chothia-numbered framework PDB.")
    heavy_chain = str(parameters.get("custom_heavy_chain", "H") or "H").strip()
    light_chain = str(parameters.get("custom_light_chain", "L") or "L").strip()
    if framework == "custom":
        if not re.fullmatch(r"[A-Za-z0-9]", heavy_chain):
            raise ValueError("The custom framework heavy-chain ID must be one character.")
        if task == "antibody" and (not re.fullmatch(r"[A-Za-z0-9]", light_chain) or light_chain == heavy_chain):
            raise ValueError("The custom framework light-chain ID must be one character different from the heavy-chain ID.")
        if require_files and framework_path is not None:
            if framework_path.suffix.lower() not in {".pdb", ".ent"}:
                raise ValueError("Custom antibody frameworks must be Chothia-numbered PDB files.")
            from mn_molecular_design_studio.design_validation import protein_chains, read_structure

            try:
                framework_chains = protein_chains(read_structure(framework_path))
            except Exception as exc:
                raise ValueError(f"Could not read custom framework {framework_path.name}: {exc}") from exc
            required_framework_chains = [heavy_chain] + ([light_chain] if task == "antibody" else [])
            missing = [chain for chain in required_framework_chains if chain not in framework_chains]
            if missing:
                raise ValueError(f"Custom framework is missing selected chain(s): {', '.join(missing)}.")
    return {
        "task": task,
        "framework": framework,
        "target_pdb": target,
        "target_chains": target_chains,
        "hotspots": hotspots,
        "regions": regions,
        "loop_lengths": loops,
        "design_loops": ",".join(f"{region}:{loops[region]}" for region in regions),
        "framework_pdb": framework_path,
        "custom_heavy_chain": heavy_chain,
        "custom_light_chain": light_chain,
        "number_of_designs": number_of_designs,
        "number_of_sequences": number_of_sequences,
        "rf2_repeats": rf2_repeats,
        "rf2_seed": rf2_seed,
        "diffusion_deterministic": diffusion_deterministic,
        "num_recycles": num_recycles,
        "num_connections": num_connections,
        "cpu_cores": cpu_cores,
        "mpnn_temperature": mpnn_temperature,
        "augment_eps": augment_eps,
        "hotspot_show_prop": hotspot_show_prop,
        "omit_aas": omit_aas,
        "bias_aas": bias_aas,
        "bias_by_res_json": bias_by_res_json,
        "use_abmpnn": use_abmpnn,
        "abmpnn_checkpoint": abmpnn_checkpoint,
        "mimosa_enabled": mimosa_enabled,
        "mimosa_condition": mimosa_condition,
        "mimosa_complex": mimosa_complex,
        "mimosa_binder_chain": mimosa_binder_chain,
        "mimosa_complex_target_chain": mimosa_complex_target_chain,
        "mimosa_target_chain": mimosa_target_chain,
        "mimosa_cdr": mimosa_cdr,
        "mimosa_design_loop": mimosa_design_loop,
        "mimosa_source_positions": mimosa_source_positions,
        "mimosa_contact_cutoff": mimosa_contact_cutoff,
        "mimosa_guidance_scale": mimosa_guidance_scale,
        "mimosa_num_particles": mimosa_num_particles,
        "mimosa_fk_lambda": mimosa_fk_lambda,
        "mimosa_resampling_interval": mimosa_resampling_interval,
    }
