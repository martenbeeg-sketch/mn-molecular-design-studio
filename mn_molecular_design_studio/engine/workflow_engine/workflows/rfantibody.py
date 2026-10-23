"""Queued, standalone RFantibody design using the pinned upstream CLI."""
from __future__ import annotations

import csv
import json
import re
import secrets
import shlex
import shutil
import subprocess
import traceback
from pathlib import Path
from typing import Literal

import numpy as np
from Bio.PDB import PDBIO, Select

from mn_molecular_design_studio import rfantibody_input as inputs
from mn_molecular_design_studio.design_validation import protein_chains, read_structure, sequence
from mn_molecular_design_studio.engine.workflow_engine.core.candidates import (
    STAGE_COMPLEX_REFOLDING,
    STAGE_GENERATION_BACKBONE,
    write_candidates,
)
from mn_molecular_design_studio.engine.workflow_engine.core.gpu import docker_gpu_args, normalize_gpu_device
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import create_job, finish_job, update_status, write_json
from mn_molecular_design_studio.runtime import reference_root

STANDARD_RESIDUES = frozenset({"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "MSE"})
PDB_CHAIN_IDS = "".join(char for char in "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789" if char not in "HLT")


class _SelectedProtein(Select):
    def __init__(self, chain_ids: set[str]):
        self.chain_ids = chain_ids

    def accept_chain(self, chain):
        return chain.id in self.chain_ids

    def accept_residue(self, residue):
        return residue.id[0] == " " and residue.resname in STANDARD_RESIDUES and "CA" in residue

    def accept_atom(self, atom):
        return atom.element != "H" and (not atom.is_disordered() or atom.get_altloc() in {"A", " "})


def _write_selected_structure(structure, destination: Path, chain_ids: set[str]) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    writer = PDBIO()
    writer.set_structure(structure)
    writer.save(str(destination), _SelectedProtein(chain_ids))
    if not destination.is_file() or "ATOM" not in destination.read_text(errors="ignore"):
        raise ValueError(f"No standard protein atoms could be written to {destination.name}.")


def _prepare_target(source: Path, selected_chains: list[str], hotspots: list[str], destination: Path) -> dict:
    structure = read_structure(source)
    model = next(structure.get_models(), None)
    if model is None:
        raise ValueError("The selected target structure has no coordinate model.")
    chains = protein_chains(structure)
    missing = [chain for chain in selected_chains if chain not in chains]
    if missing:
        raise ValueError(f"Selected antigen chain(s) are missing from the target: {', '.join(missing)}.")

    used: set[str] = set()
    chain_map: dict[str, str] = {}
    for source_chain in selected_chains:
        if len(source_chain) == 1 and source_chain in PDB_CHAIN_IDS and source_chain not in used:
            output_chain = source_chain
        else:
            output_chain = next((char for char in PDB_CHAIN_IDS if char not in used), "")
        if not output_chain:
            raise ValueError("The selected target uses more chains than the PDB format can represent.")
        used.add(output_chain)
        chain_map[source_chain] = output_chain

    mapped_hotspots = []
    for hotspot in hotspots:
        match = re.fullmatch(r"([A-Za-z0-9]+?)(-?\d+)", hotspot)
        source_chain, residue_number = match.groups()
        if source_chain not in chain_map:
            raise ValueError(f"Hotspot {hotspot} is not part of the selected antigen chains.")
        residues = chains[source_chain]
        if not any(int(residue.id[1]) == int(residue_number) and not str(residue.id[2]).strip() for residue in residues):
            raise ValueError(f"Hotspot {hotspot} is missing or has an insertion code, which RFantibody cannot address.")
        mapped_hotspots.append(f"{chain_map[source_chain]}{int(residue_number)}")

    for chain in list(model):
        if chain.id not in chain_map:
            model.detach_child(chain.id)
            continue
        original_id = chain.id
        chain.id = chain_map[original_id]
        for residue in list(chain):
            if residue.id[0] != " " or residue.resname not in STANDARD_RESIDUES or "CA" not in residue:
                chain.detach_child(residue.id)
    _write_selected_structure(structure, destination, set(chain_map.values()))
    return {"chain_map": chain_map, "target_chains": list(chain_map.values()), "hotspots": mapped_hotspots}


def _score_from_pdb(path: Path) -> dict[str, float]:
    scores: dict[str, float] = {}
    for line in path.read_text(errors="ignore").splitlines():
        match = re.match(r"SCORE\s+([A-Za-z0-9_]+):\s*(-?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)", line)
        if match:
            scores[match.group(1)] = float(match.group(2))
    return scores


def _write_binder_structure(source: Path, destination: Path, chains: list[str]) -> None:
    structure = read_structure(source)
    model = next(structure.get_models(), None)
    if model is None:
        return
    allowed = set(chains)
    for chain in list(model):
        if chain.id not in allowed:
            model.detach_child(chain.id)
    _write_selected_structure(structure, destination, allowed)


def _read_quiver(path: Path) -> list[tuple[str, list[str]]]:
    """Read the pinned RFantibody text Quiver format without importing its CLI package."""
    records: list[tuple[str, list[str]]] = []
    tag: str | None = None
    pdb_lines: list[str] = []
    if not path.is_file():
        return records
    for line in path.read_text(errors="replace").splitlines(keepends=True):
        if line.startswith("QV_TAG "):
            if tag is not None:
                records.append((tag, pdb_lines))
            tag = line.split(maxsplit=1)[1].strip()
            pdb_lines = []
        elif tag is not None and not line.startswith("QV_SCORE "):
            pdb_lines.append(line)
    if tag is not None:
        records.append((tag, pdb_lines))
    return records


def _mimosa_generation_metrics(metadata: dict) -> dict:
    fixed = metadata.get("fixed_residues") if isinstance(metadata.get("fixed_residues"), list) else []
    if metadata.get("condition") == "geometry_only":
        fixed = []
    metrics = {
        "generated_motif_rmsd_angstrom": metadata.get("generated_motif_rmsd_angstrom"),
        "generated_motif_potential_half_sum_angstrom2": metadata.get("generated_motif_potential_half_sum_angstrom_squared"),
        "fixed_motif_residue_count": len(fixed),
        "backbone_integrity_status": metadata.get("backbone_integrity_status"),
        "backbone_bad_bond_count": metadata.get("backbone_bad_bond_count"),
        "backbone_checked_bond_count": metadata.get("backbone_checked_bond_count"),
    }
    return {key: value for key, value in metrics.items() if value is not None}


def _refolded_motif_rmsd(run_dir: Path, generation_tag: str, structure_path: Path) -> float | None:
    """Measure the selected source motif in the target-aligned RF2 structure."""
    stages_dir = run_dir / "artifacts/raw/rfantibody/stages"
    # The sampler writes the tag returned by the QV reader under its stage
    # prefix (for example, ``samples_design_0`` becomes
    # ``1_rfdiffusion_samples_design_0.mimosa.json``). Keep the old unprefixed
    # spelling as a read-only fallback for earlier artifacts.
    safe_tag = Path(str(generation_tag)).name
    if safe_tag != str(generation_tag) or not safe_tag:
        return None
    metadata_path = next(
        (path for path in (
            stages_dir / f"1_rfdiffusion_{safe_tag}.mimosa.json",
            stages_dir / f"{safe_tag}.mimosa.json",
        ) if path.is_file()),
        None,
    )
    if metadata_path is None:
        return None
    try:
        metadata = json.loads(metadata_path.read_text())
        alignment = metadata.get("target_alignment") or []
        assignments = metadata.get("motif_assignments") or []
        guided_positions = metadata.get("guided_loop_positions") or []
        source = metadata.get("source") or {}
        target_sequence = str(source.get("target_sequence") or "")
        if (len(alignment) < 3 or not assignments or not guided_positions or not target_sequence
                or len(assignments) != int(metadata.get("guided_motif_residues", len(assignments)))):
            return None
        structure = read_structure(structure_path)
        chains = protein_chains(structure)
        binder_chains = {str((row.get("generated_position") or {}).get("chain_id", "")) for row in assignments}
        target_candidates = [
            (chain_id, chain) for chain_id, chain in chains.items()
            if chain_id not in binder_chains and len(chain) == len(target_sequence)
        ]
        target_chain_id, target_residues = next(
            ((chain_id, chain) for chain_id, chain in target_candidates
             if sequence(chain) == target_sequence), (None, None)
        )
        if target_chain_id is None:
            return None
        ref_points, current_points = [], []
        target_residues = list(target_residues)
        for row in alignment:
            offset = int(row["target_offset"])
            if offset < 0 or offset >= len(target_residues):
                continue
            residue = target_residues[offset]
            if "CA" not in residue:
                continue
            ref_points.append(row["ca"])
            current_points.append(residue["CA"].coord)
        if len(ref_points) < 3:
            return None
        ref = np.asarray(ref_points, dtype=float)
        current = np.asarray(current_points, dtype=float)
        if len(ref_points) != len(alignment) or len(current_points) != len(alignment):
            return None
        ref_center, current_center = ref.mean(axis=0), current.mean(axis=0)
        u, _, vh = np.linalg.svd((ref - ref_center).T @ (current - current_center))
        correction = np.eye(3)
        correction[-1, -1] = np.sign(np.linalg.det(u @ vh))
        rotation = u @ correction @ vh
        translation = current_center - ref_center @ rotation
        reference_points = np.asarray([row["source_ca"] for row in assignments], dtype=float) @ rotation + translation
        generated_points = []
        # The RFdiffusion sidecar's author residue IDs use the numbering
        # carried by the sampler (often Chothia/IMGT), while the QV writer
        # renumbers output PDB residues consecutively. Prefer explicit pose
        # indices, and otherwise derive the per-chain offset from the matched
        # motif assignments. Only use author numbering as a legacy fallback.
        pose_offsets: dict[str, int] = {}
        offset_values: dict[str, set[int]] = {}
        for row in assignments:
            position = row.get("generated_position") or {}
            chain_id = str(position.get("chain_id", ""))
            try:
                pose_index = int(row.get("generated_pose_index", row.get("zero_based_pose_index")))
                author_index = int(position.get("author_residue_id"))
            except (TypeError, ValueError):
                continue
            if chain_id:
                offset_values.setdefault(chain_id, set()).add(pose_index - author_index)
        pose_offsets = {chain_id: next(iter(offsets)) for chain_id, offsets in offset_values.items()
                        if len(offsets) == 1}
        # Pose indices address the sampler's concatenated C-alpha array. The
        # QV artifact is the authoritative bridge from that order to the
        # writer-renumbered chain-local residue IDs used by RF2 output.
        qv_records = _read_quiver(stages_dir / "1_rfdiffusion.qv")
        qv_ca: list[tuple[str, int]] = []
        qv_chain_counts: dict[str, int] = {}
        for line in next((lines for tag, lines in qv_records if tag == safe_tag), []):
            if not line.startswith("ATOM") or line[12:16].strip() != "CA":
                continue
            chain_id = line[21].strip()
            qv_ca.append((chain_id, qv_chain_counts.get(chain_id, 0)))
            qv_chain_counts[chain_id] = qv_chain_counts.get(chain_id, 0) + 1
        if not qv_ca:
            return None
        if any(qv_chain_counts.get(chain_id) != len(chains.get(chain_id, []))
               for chain_id in qv_chain_counts):
            return None
        for position in guided_positions:
            chain_id = str(position.get("chain_id", ""))
            if chain_id not in chains:
                return None
            residue = None
            explicit_pose_index = position.get("zero_based_pose_index", position.get("generated_pose_index"))
            if explicit_pose_index is not None:
                try:
                    pose_index = int(explicit_pose_index)
                except (TypeError, ValueError):
                    return None
                if not 0 <= pose_index < len(qv_ca) or qv_ca[pose_index][0] != chain_id:
                    return None
                residue = chains[chain_id][qv_ca[pose_index][1]]
            elif chain_id in pose_offsets:
                # Legacy sidecars did not store the index for every loop
                # position. Their author-number fallback is only unambiguous
                # for a single-chain, insertion-free loop with a consistent
                # assignment offset.
                if ({str(item.get("chain_id", "")) for item in guided_positions} != {chain_id}
                        or any(str(item.get("insertion_code", "") or "") for item in guided_positions)):
                    return None
                try:
                    pose_index = int(position.get("author_residue_id")) + pose_offsets[chain_id]
                except (TypeError, ValueError):
                    return None
                if not 0 <= pose_index < len(qv_ca) or qv_ca[pose_index][0] != chain_id:
                    return None
                residue = chains[chain_id][qv_ca[pose_index][1]]
            else:
                try:
                    residue_number = int(position.get("author_residue_id"))
                except (ValueError, TypeError):
                    return None
                insertion_code = str(position.get("insertion_code", " ") or " ")[:1]
                residue = next((candidate for candidate in chains[chain_id]
                                if int(candidate.id[1]) == residue_number and str(candidate.id[2] or " ") == insertion_code), None)
            if residue is None or "CA" not in residue:
                return None
            generated_points.append(np.asarray(residue["CA"].coord, dtype=float))
        generated = np.asarray(generated_points, dtype=float)
        if (not np.all(np.isfinite(reference_points)) or not np.all(np.isfinite(generated))
                or len(generated) < len(reference_points)):
            return None
        cost = np.linalg.norm(reference_points[:, None, :] - generated[None, :, :], axis=-1)
        if metadata.get("assignment_method") == "sliding_window":
            motif_size = len(reference_points)
            window_cost = np.asarray([
                np.diag(cost[:, start:start + motif_size]).sum()
                for start in range(generated.shape[0] - motif_size + 1)
            ])
            selected = np.arange(int(np.argmin(window_cost)), int(np.argmin(window_cost)) + motif_size)
            distances = cost[np.arange(motif_size), selected]
        else:
            from scipy.optimize import linear_sum_assignment
            rows, columns = linear_sum_assignment(cost)
            if len(rows) != len(reference_points):
                return None
            distances = cost[rows, columns]
        return float(np.sqrt(np.mean(np.square(distances))))
    except (OSError, ValueError, KeyError, TypeError, StopIteration):
        return None


def _collect_generation_candidates(
    run_dir: Path,
    target_pdb: Path,
    target_chain_ids: list[str],
    hotspots: list[str],
) -> list[dict]:
    """Preserve accepted and structurally rejected RFdiffusion outputs as generation records."""
    root = run_dir / "artifacts/raw/rfantibody"
    stages = root / "stages"
    structures = root / "structures/generation"
    report: dict[str, dict] = {}
    report_path = stages / "mimosa_integrity_results.json"
    if report_path.is_file():
        try:
            rows = json.loads(report_path.read_text()).get("candidates", [])
            report = {str(row.get("tag")): row for row in rows if isinstance(row, dict) and row.get("tag")}
        except (OSError, ValueError, TypeError):
            report = {}

    generated_files: dict[str, Path] = {}
    qv_path = stages / "1_rfdiffusion.qv"
    for tag, lines in _read_quiver(qv_path):
        if not any(line.startswith(("ATOM  ", "HETATM")) for line in lines):
            continue
        safe_tag = re.sub(r"[^A-Za-z0-9_.-]+", "_", tag)
        destination = structures / f"{safe_tag}.pdb"
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text("".join(lines))
        generated_files[tag] = destination
    rejected_dir = stages / "rejected_backbones"
    if rejected_dir.is_dir():
        for path in sorted(rejected_dir.glob("*.pdb")):
            generated_files.setdefault(path.stem, path)

    candidates: list[dict] = []
    for index, (tag, structure_path) in enumerate(sorted(generated_files.items()), 1):
        try:
            structure = read_structure(structure_path)
            chains = protein_chains(structure)
            binder_chains = [chain for chain in ("H", "L") if chain in chains]
            target_chains = [chain for chain in target_chain_ids if chain in chains]
            if not binder_chains:
                continue
            metadata = report.get(tag, {})
            metrics = _mimosa_generation_metrics(metadata)
            if metadata and metadata.get("backbone_integrity_status") == "rejected":
                metrics["result_kind"] = "structural_rejection"
            else:
                metrics["result_kind"] = "generated"
            diagnostics = metadata.get("sampler_diagnostics") if isinstance(metadata.get("sampler_diagnostics"), dict) else {}
            candidate_id = f"rfdiffusion_{re.sub(r'[^A-Za-z0-9_.-]+', '_', tag)}"
            candidates.append({
                "candidate_id": candidate_id,
                "source_tool": "rfantibody",
                "stage": STAGE_GENERATION_BACKBONE,
                "target_pdb": str(target_pdb),
                "complex_pdb": str(structure_path),
                "binder_pdb": None,
                "binder_sequence": None,
                "target_chains": target_chains or target_chain_ids,
                "binder_chains": binder_chains,
                "hotspots": hotspots,
                "binder_length": sum(len(chains[chain]) for chain in binder_chains),
                "metrics": metrics,
                "raw_metadata": {
                    "result_kind": "generation_backbone",
                    "generation_tag": tag,
                    "mimosa_diagnostics": diagnostics,
                    "backbone_bad_bonds": metadata.get("backbone_bad_bonds", []),
                },
            })
        except (OSError, ValueError, KeyError) as exc:
            (root / "collection_warnings.txt").open("a").write(
                f"generation {tag}: {type(exc).__name__}: {exc}\n"
            )
    return candidates


def _proteinmpnn_command(
    *,
    input_quiver: str,
    output_quiver: str,
    loops: str,
    sequences_per_structure: int,
    temperature: float,
    augment_eps: float,
    num_connections: int,
    omit_aas: str,
    checkpoint: str,
) -> list[str]:
    """Use the pinned RFantibody script's actual ProteinMPNN argument names."""
    return [
        "python", "scripts/proteinmpnn_interface_design.py",
        "-quiver", input_quiver,
        "-outquiver", output_quiver,
        "-loop_string", loops,
        "-seqs_per_struct", str(sequences_per_structure),
        "-temperature", str(temperature),
        "-augment_eps", str(augment_eps),
        "-num_connections", str(num_connections),
        "-omit_AAs", omit_aas,
        "-checkpoint_path", checkpoint,
    ]


def _collect_candidates(
    run_dir: Path,
    target_pdb: Path,
    target_chain_ids: list[str],
    hotspots: list[str],
    *,
    calculate_epitope_distance: bool,
    generation_candidates: list[dict] | None = None,
) -> list[dict]:
    output_root = run_dir / "artifacts" / "raw" / "rfantibody"
    structures_dir = output_root / "structures"
    candidates = []
    generated_by_tag = {
        str(candidate.get("raw_metadata", {}).get("generation_tag")): candidate
        for candidate in (generation_candidates or [])
        if isinstance(candidate.get("raw_metadata"), dict)
    }
    pdb_paths = sorted(output_root.glob("rf2_models/model_*/*_best.pdb"))
    for index, path in enumerate(pdb_paths, 1):
        try:
            structure = read_structure(path)
            chains = protein_chains(structure)
            binder_chains = [chain for chain in ("H", "L") if chain in chains]
            target_chains = [chain for chain in chains if chain not in binder_chains]
            if not binder_chains or not target_chains:
                continue
            model_match = re.search(r"model_(\d+)", path.parent.name)
            rf2_repeat = int(model_match.group(1)) if model_match else index
            tag_match = re.match(r"(.+?)_dldesign_", path.stem)
            generation_tag = tag_match.group(1) if tag_match else ""
            generated_parent = generated_by_tag.get(generation_tag)
            candidate_id = f"rf2_{rf2_repeat:02d}_{path.stem}"
            binder_sequence = "".join(sequence(chains[chain]) for chain in binder_chains)
            binder_path = structures_dir / f"{candidate_id}_binder.pdb"
            _write_binder_structure(path, binder_path, binder_chains)
            metrics = _score_from_pdb(path)
            if generated_parent:
                refolded_motif_rmsd = _refolded_motif_rmsd(run_dir, generation_tag, path)
                if refolded_motif_rmsd is not None:
                    metrics["refolded_motif_rmsd_angstrom"] = refolded_motif_rmsd
            raw_hotspots = hotspots if calculate_epitope_distance else []
            candidates.append({
                "candidate_id": candidate_id,
                "source_tool": "rfantibody",
                "stage": STAGE_COMPLEX_REFOLDING,
                "target_pdb": str(target_pdb),
                "complex_pdb": str(path),
                "binder_pdb": str(binder_path),
                "binder_sequence": binder_sequence,
                "target_chains": target_chains,
                "binder_chains": binder_chains,
                "hotspots": raw_hotspots,
                "binder_length": len(binder_sequence),
                "metrics": {
                    **metrics,
                    "rf2_repeat": rf2_repeat,
                    "binder_length": len(binder_sequence),
                    "hotspot_distance_analysis": bool(calculate_epitope_distance),
                },
                "raw_metadata": {
                    "upstream_commit": inputs.UPSTREAM_COMMIT,
                    "input_target_chains": target_chain_ids,
                    "input_hotspots": hotspots,
                    "rf2_repeat": rf2_repeat,
                    "result_kind": "rf2_complex_prediction",
                },
                "parents": [generated_parent["candidate_id"]] if generated_parent else [],
            })
        except (OSError, ValueError, KeyError) as exc:
            (output_root / "collection_warnings.txt").open("a").write(f"{path.name}: {type(exc).__name__}: {exc}\n")
    return sorted(candidates, key=lambda candidate: (-float(candidate["metrics"].get("pred_lddt", -1)), float(candidate["metrics"].get("interaction_pae", 1e9))))


def run_rfantibody_design(
    task: Literal["antibody", "nanobody"] = "antibody",
    target_pdb: Path | None = None,
    target_chains: list[str] | None = None,
    hotspots: str = "",
    framework: str = "hu-4D5-8_Fv",
    framework_pdb: Path | None = None,
    custom_heavy_chain: str = "H",
    custom_light_chain: str = "L",
    regions: list[str] | None = None,
    loop_lengths: dict[str, str] | None = None,
    number_of_designs: int = 1,
    number_of_sequences: int = 1,
    use_abmpnn: bool = False,
    generation_only: bool = False,
    calculate_epitope_distance: bool = True,
    mpnn_temperature: float = 0.1,
    bias_aas: str = "",
    bias_by_res_json: str = "{}",
    omit_aas: str = "CX",
    augment_eps: float = 0.05,
    num_connections: int = 48,
    num_recycles: int = 3,
    rf2_repeats: int = 5,
    hotspot_show_prop: float = 0.1,
    save_trajectories: bool = False,
    cpu_cores: int = 4,
    gpu_device: str = "0",
    job_name: str = "rfantibody-design",
    mimosa_enabled: bool = False,
    mimosa_condition: str = "sequence_fixed",
    mimosa_complex: Path | None = None,
    mimosa_binder_chain: str = "",
    mimosa_complex_target_chain: str = "",
    mimosa_target_chain: str = "",
    mimosa_cdr: str = "CDR3",
    mimosa_design_loop: str = "H3",
    mimosa_source_positions: str | list[str] = "",
    mimosa_contact_cutoff: float = 8.0,
    mimosa_guidance_scale: float = 0.3,
    mimosa_num_particles: int = 4,
    mimosa_fk_lambda: float = 4.0,
    mimosa_resampling_interval: int = 3,
    rf2_seed: int | None = None,
    diffusion_deterministic: bool = False,
) -> Path:
    parameters = dict(locals())
    normalized = inputs.validate_request(parameters, require_files=True)
    if normalized["target_pdb"] is None:
        raise ValueError("Load a target PDB or mmCIF structure first.")
    weights_dir = reference_root() / "rfantibody"
    required_weights = ["RFdiffusion_Ab.pt"]
    if not normalized["generation_only"]:
        required_weights.append("RF2_ab.pt")
        if not normalized["use_abmpnn"]:
            required_weights.append("ProteinMPNN_v48_noise_0.2.pt")
    missing_weights = [name for name in required_weights if not (weights_dir / name).is_file()]
    if missing_weights:
        raise ValueError(
            f"RFantibody checkpoints are missing from {weights_dir}: "
            + ", ".join(missing_weights)
        )
    image = inputs.MIMOSA_IMAGE if normalized["mimosa_enabled"] else inputs.IMAGE
    if normalized["mimosa_enabled"]:
        image_check = subprocess.run(
            ["docker", "image", "inspect", image],
            capture_output=True,
            text=True,
            check=False,
        )
        if image_check.returncode != 0:
            raise RuntimeError(
                f"MIMOSA runtime image {image!r} is unavailable. "
                "Restore or rebuild the Studio runtime image before submitting this job."
            )
    seed = normalized["rf2_seed"] if normalized["rf2_seed"] is not None else secrets.randbelow(2_000_000_000) + 1
    gpu_device = normalize_gpu_device(gpu_device)
    parameters.update({
        "target_pdb": str(normalized["target_pdb"]),
        "framework_pdb": str(normalized["framework_pdb"]) if normalized["framework_pdb"] else None,
        "mimosa_complex": str(normalized["mimosa_complex"]) if normalized["mimosa_complex"] else None,
        "gpu_device": gpu_device,
        "cpu_cores": normalized["cpu_cores"],
        "rf2_seed": seed,
        "diffusion_deterministic": bool(normalized["diffusion_deterministic"]),
        "use_abmpnn": bool(normalized["use_abmpnn"]),
        "diffusion_seed": 0 if normalized["diffusion_deterministic"] else None,
        "engine_version": inputs.CODE_VERSION,
        "image": image,
        "generation_only": bool(normalized["generation_only"]),
    })
    inputs_payload = {
        "target_pdb": str(normalized["target_pdb"]),
        "framework_pdb": str(normalized["framework_pdb"]) if normalized["framework_pdb"] else None,
        "mimosa_complex": str(normalized["mimosa_complex"]) if normalized["mimosa_complex"] else None,
    }
    job = create_job(
        "design",
        "design_generation" if normalized["generation_only"] else "rfantibody_design",
        "rfantibody",
        inputs_payload,
        parameters,
    )
    update_status(job.run_dir, "preparing", campaign_name=job_name, current_phase="Preparing RFantibody target and framework")
    root = job.run_dir / "artifacts" / "raw" / "rfantibody"
    inputs_dir = root / "inputs"
    inputs_dir.mkdir(parents=True, exist_ok=True)
    target_path = inputs_dir / "target.pdb"
    try:
        mapping = _prepare_target(normalized["target_pdb"], normalized["target_chains"], normalized["hotspots"], target_path)
        shutil.copy2(normalized["target_pdb"], inputs_dir / f"source_target{normalized['target_pdb'].suffix.lower()}")
        mimosa_packet_path = None
        if normalized["mimosa_enabled"]:
            mimosa_packet_path = inputs_dir / "mimosa_motif.json"
            parameters["mimosa_packet"] = str(mimosa_packet_path)
            source_complex_path = inputs_dir / f"mimosa_source_complex{normalized['mimosa_complex'].suffix.lower()}"
            shutil.copy2(normalized["mimosa_complex"], source_complex_path)
            parameters["mimosa_source_complex"] = str(source_complex_path)
            write_json(
                job.run_dir / "input.json",
                {
                    "job_type": "rfantibody_design",
                    "tool": "rfantibody",
                    "inputs": inputs_payload,
                    "params": parameters,
                },
            )
        if normalized["framework"] == "custom":
            shutil.copy2(normalized["framework_pdb"], inputs_dir / "custom_framework.pdb")

        design_loops = normalized["design_loops"]
        hotspot_text = ",".join(mapping["hotspots"])
        diff_qv = root / "stages" / "1_rfdiffusion.qv"
        mpnn_qv = root / "stages" / "2_proteinmpnn.qv"
        predictions_root = root / "rf2_models"
        output_root = root / "structures"
        (diff_qv.parent).mkdir(parents=True, exist_ok=True)
        predictions_root.mkdir(parents=True, exist_ok=True)
        output_root.mkdir(parents=True, exist_ok=True)
        target_in_container = "/work/artifacts/raw/rfantibody/inputs/target.pdb"
        diff_qv_container = "/work/artifacts/raw/rfantibody/stages/1_rfdiffusion.qv"
        mpnn_qv_container = "/work/artifacts/raw/rfantibody/stages/2_proteinmpnn.qv"
        if normalized["framework"] == "custom":
            framework_container = "/work/artifacts/raw/rfantibody/inputs/framework_hlt.pdb"
            framework_command = [
                "python", "scripts/util/chothia2HLT.py", "/work/artifacts/raw/rfantibody/inputs/custom_framework.pdb",
                "--heavy", normalized["custom_heavy_chain"], "--output", framework_container,
            ]
            if normalized["task"] == "antibody":
                framework_command.extend(["--light", normalized["custom_light_chain"]])
        else:
            framework_name = "h-NbBCII10" if normalized["framework"] == "h-NbBCII10" else "hu-4D5-8_Fv"
            framework_container = f"scripts/examples/example_inputs/{framework_name}.pdb"
            framework_command = []

        image_mounts = [
            "--mount", f"type=bind,src={job.run_dir.resolve()},dst=/work",
            "--mount", f"type=bind,src={weights_dir.resolve()},dst=/ref/rfantibody,readonly",
        ]
        mpnn_checkpoint_in_container = "/ref/rfantibody/ProteinMPNN_v48_noise_0.2.pt"
        if normalized["use_abmpnn"]:
            checkpoint = normalized["abmpnn_checkpoint"]
            image_mounts.extend(["--mount", f"type=bind,src={checkpoint.resolve()},dst=/mnt/studio-weights/abmpnn.pt,readonly"])
            mpnn_checkpoint_in_container = "/mnt/studio-weights/abmpnn.pt"

        command = [
            "docker", "run", "--rm", "--name", f"mn-rfantibody-{job.run_id[-8:]}",
            "--cpus", str(normalized["cpu_cores"]), "--shm-size=16G", *docker_gpu_args(gpu_device),
            *image_mounts, "-w", "/opt/rfantibody",
            "-e", f"OMP_NUM_THREADS={normalized['cpu_cores']}", "-e", "PYTHONUNBUFFERED=1",
            "-e", "RFANTIBODY_WEIGHTS=/ref/rfantibody",
            image, "bash", "-euo", "pipefail", "-c",
        ]
        steps: list[str] = []
        if normalized["mimosa_enabled"]:
            mapped_target_chain = mapping["chain_map"][normalized["mimosa_target_chain"]]
            prepare_motif = [
                "python", "scripts/prepare_mimosa_motif.py",
                "--complex", f"/work/artifacts/raw/rfantibody/inputs/mimosa_source_complex{normalized['mimosa_complex'].suffix.lower()}",
                "--target", target_in_container,
                "--binder-chain", normalized["mimosa_binder_chain"],
                "--complex-target-chain", normalized["mimosa_complex_target_chain"],
                "--target-chain", mapped_target_chain,
                "--cdr", normalized["mimosa_cdr"],
                "--contact-cutoff", str(normalized["mimosa_contact_cutoff"]),
                "--guidance-scale", str(normalized["mimosa_guidance_scale"]),
                "--fk-lambda", str(normalized["mimosa_fk_lambda"]),
                "--resampling-interval", str(normalized["mimosa_resampling_interval"]),
                "--condition", normalized["mimosa_condition"],
                "--design-loop", normalized["mimosa_design_loop"],
                "--output", "/work/artifacts/raw/rfantibody/inputs/mimosa_motif.json",
            ]
            if normalized["mimosa_source_positions"]:
                prepare_motif.extend(
                    ["--source-residues", *normalized["mimosa_source_positions"]]
                )
            steps.append("echo '[MIMOSA] Preparing the cognate CDR motif with containerized SAbR'; " + shlex.join(prepare_motif))
        if framework_command:
            steps.append("echo '[RFantibody] Converting the custom antibody framework'; " + shlex.join(framework_command))
        diffusion = [
            "rfdiffusion", "--target", target_in_container, "--framework", framework_container,
            "--output-quiver", diff_qv_container, "--num-designs", str(normalized["number_of_designs"]),
            "--design-loops", design_loops, "--diffuser-t", "50",
        ]
        if normalized["diffusion_deterministic"]:
            diffusion.append("--deterministic")
        if normalized["mimosa_enabled"]:
            mimosa_packet_container = "/work/artifacts/raw/rfantibody/inputs/mimosa_motif.json"
            diffusion.extend(
                [
                    "--mimosa-packet", mimosa_packet_container,
                    "--mimosa-loop", normalized["mimosa_design_loop"],
                    "--mimosa-guidance-scale", str(normalized["mimosa_guidance_scale"]),
                    "--mimosa-particles", str(normalized["mimosa_num_particles"]),
                    "--mimosa-fk-lambda", str(normalized["mimosa_fk_lambda"]),
                    "--mimosa-resampling-interval", str(normalized["mimosa_resampling_interval"]),
                    "--mimosa-condition", normalized["mimosa_condition"],
                ]
            )
        if hotspot_text:
            diffusion.extend(["--hotspots", hotspot_text])
        if not save_trajectories:
            diffusion.append("--no-trajectory")
        steps.extend(["echo '[RFantibody] Stage 1/3 · RFdiffusion'; " + shlex.join(diffusion)])
        integrity_report_container = "/work/artifacts/raw/rfantibody/stages/mimosa_integrity_results.json"
        if normalized["generation_only"]:
            steps.append("echo '[RFantibody] Generation-only run; skipping internal sequence design and RF2 predictions'")
        else:
            mpnn = _proteinmpnn_command(
                input_quiver=diff_qv_container,
                output_quiver=mpnn_qv_container,
                loops=",".join(normalized["regions"]),
                sequences_per_structure=normalized["number_of_sequences"],
                temperature=normalized["mpnn_temperature"],
                augment_eps=normalized["augment_eps"],
                num_connections=normalized["num_connections"],
                omit_aas=normalized["omit_aas"],
                checkpoint=mpnn_checkpoint_in_container,
            )
        if not normalized["generation_only"] and normalized["mimosa_enabled"]:
            all_rejected = (
                "python -c " + shlex.quote(
                    "import json,sys; d=json.load(open(" + repr(integrity_report_container) + ")); "
                    "sys.exit(0 if d.get('candidate_count',0)>0 and d.get('candidate_count')==d.get('rejected_count') else 1)"
                )
            )
            steps.append(
                f"if {all_rejected}; then echo '[RFantibody] All generated backbones failed the MIMOSA peptide-bond integrity check; skipping inverse folding and refolding'; "
                f"else echo '[RFantibody] Stage 2/3 · ProteinMPNN'; {shlex.join(mpnn)}; fi"
            )
        elif not normalized["generation_only"]:
            steps.append("echo '[RFantibody] Stage 2/3 · ProteinMPNN'; " + shlex.join(mpnn))
        if not normalized["generation_only"]:
            for repeat in range(1, normalized["rf2_repeats"] + 1):
                model_output = f"/work/artifacts/raw/rfantibody/rf2_models/model_{repeat:02d}"
                rf2 = [
                    "rf2", "--input-quiver", mpnn_qv_container, "--output-dir", model_output,
                    "--num-recycles", str(normalized["num_recycles"]), "--hotspot-show-prop", str(normalized["hotspot_show_prop"]),
                    "--seed", str(seed + repeat - 1), "--no-cautious",
                ]
                rf2_command = shlex.join(rf2)
                if normalized["mimosa_enabled"]:
                    steps.append(
                        f"if {all_rejected}; then echo '[RFantibody] Skipping RF2 for structurally rejected designs'; "
                        f"else echo '[RFantibody] Stage 3/3 · RF2 prediction {repeat}/{normalized['rf2_repeats']}'; {rf2_command}; fi"
                    )
                else:
                    steps.append(f"echo '[RFantibody] Stage 3/3 · RF2 prediction {repeat}/{normalized['rf2_repeats']}'; " + rf2_command)
        script = "\n".join(steps)
        write_json(job.run_dir / "command.json", {"mode": "docker", "image": image, "command": command + [script]})
        update_status(
            job.run_dir,
            "running",
            current_phase="RFdiffusion backbone generation only"
            if normalized["generation_only"]
            else "RFdiffusion → ProteinMPNN → RF2",
        )
        with (job.run_dir / "stdout.log").open("a") as stdout, (job.run_dir / "stderr.log").open("a") as stderr:
            process = subprocess.run(command + [script], stdout=stdout, stderr=stderr, cwd=job.run_dir, check=False)
        update_status(
            job.run_dir,
            "running",
            current_phase="Collecting generated RFantibody backbones"
            if normalized["generation_only"]
            else "Collecting RFantibody designs and RF2 scores",
        )
        generation_candidates = _collect_generation_candidates(
            job.run_dir, target_path, mapping["target_chains"], mapping["hotspots"]
        )
        candidates = [] if normalized["generation_only"] else _collect_candidates(
            job.run_dir, target_path, mapping["target_chains"], mapping["hotspots"],
            calculate_epitope_distance=calculate_epitope_distance,
            generation_candidates=generation_candidates,
        )
        candidates = [*generation_candidates, *candidates]
        normalized_candidates = write_candidates(job.run_dir, "rfantibody", candidates) if candidates else []
        table_path = root / "designs.csv"
        with table_path.open("w", newline="") as handle:
            fields = ["rank", "candidate_id", "binder_sequence", "binder_chains", "pred_lddt", "interaction_pae", "rf2_repeat", "structure"]
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            for rank, candidate in enumerate(
                [row for row in normalized_candidates if row.get("stage") == STAGE_COMPLEX_REFOLDING], 1
            ):
                metrics = candidate.get("metrics") or {}
                writer.writerow({"rank": rank, "candidate_id": candidate["candidate_id"], "binder_sequence": candidate.get("binder_sequence"),
                                 "binder_chains": ",".join(candidate.get("binder_chains") or []), "pred_lddt": metrics.get("pred_lddt"),
                                 "interaction_pae": metrics.get("interaction_pae"), "rf2_repeat": metrics.get("rf2_repeat"),
                                 "structure": candidate.get("complex_pdb")})
        artifacts = [{"path": str(path), "name": path.name, "type": "rfantibody_output"}
                     for path in root.rglob("*") if path.is_file()]
        integrity_report_path = root / "stages" / "mimosa_integrity_results.json"
        scientific_rejection = False
        if normalized["mimosa_enabled"] and integrity_report_path.is_file():
            try:
                report = json.loads(integrity_report_path.read_text())
                scientific_rejection = bool(report.get("candidate_count")) and report.get("candidate_count") == report.get("rejected_count")
            except (OSError, ValueError, TypeError):
                scientific_rejection = False
        success = process.returncode == 0 and (bool(normalized_candidates) or scientific_rejection)
        finish_job(job.run_dir, success, {
            "outputs": {"artifacts": artifacts, "candidates": normalized_candidates, "design_table": str(table_path)},
            "metrics": {"return_code": process.returncode, "candidate_count": len(normalized_candidates),
                        "result_status": "completed_with_structural_rejection" if scientific_rejection else "completed",
                        "requested_designs": normalized["number_of_designs"], "sequences_per_structure": normalized["number_of_sequences"],
                        "rf2_repeats": normalized["rf2_repeats"], "engine_version": inputs.CODE_VERSION,
                        "generation_only": bool(normalized["generation_only"]),
                        "rf2_predictions_completed": 0 if normalized["generation_only"] else normalized["rf2_repeats"],
                        "selected_target_chain_map": mapping["chain_map"], "selected_hotspots": mapping["hotspots"]},
            "downstream_artifacts": {"candidates_jsonl": "artifacts/normalized_candidates/candidates.jsonl",
                                     "campaign_result": "artifacts/normalized_candidates/campaign_result.json",
                                     "design_table": str(table_path.relative_to(job.run_dir))},
        })
    except Exception as exc:
        with (job.run_dir / "stderr.log").open("a") as handle:
            handle.write(traceback.format_exc())
        finish_job(job.run_dir, False, {"error": str(exc), "metrics": {"error": str(exc), "engine_version": inputs.CODE_VERSION}})
    return job.run_dir
