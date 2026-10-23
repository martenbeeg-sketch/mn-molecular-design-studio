"""Streaming MDTraj analysis for MD trajectories selected in Molecular Design Studio.

This file is copied into the CPU MD container by the workflow wrapper.  Keep it
self-contained apart from the dependencies already bundled in ``mn-md``.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import struct
import traceback
from pathlib import Path
from typing import Any, Iterable

import numpy as np


_AKMA_PS = 0.04888821


def _dcd_timing(path: Path, *, frame_interval_ps: object = None, first_time_ps: object = None) -> dict[str, Any]:
    """Read CHARMM CORD timing written by OpenMM, or validate an explicit override.

    MDTraj exposes frame-number-like times for some DCDs, so those values are
    never used for DCD analysis. The DCD DELTA field is in AKMA time units.
    """
    if frame_interval_ps is not None or first_time_ps is not None:
        if frame_interval_ps is None or first_time_ps is None:
            raise ValueError("DCD timing override requires both frame_interval_ps and first_time_ps.")
        interval = _finite_positive(frame_interval_ps, "DCD frame interval")
        first = _finite_positive(first_time_ps, "DCD first frame time", allow_zero=True)
        return {"source": "explicit_override", "frame_interval_ps": interval, "first_time_ps": first}

    with path.open("rb") as handle:
        header = handle.read(92)
    if len(header) < 92:
        raise ValueError("DCD header is truncated; provide an explicit validated timing override.")
    endian = "<" if struct.unpack("<i", header[:4])[0] == 84 else ">" if struct.unpack(">i", header[:4])[0] == 84 else None
    if endian is None or header[4:8] != b"CORD":
        raise ValueError("Unsupported DCD header; expected CHARMM CORD. Provide explicit validated timing overrides.")
    nset, first_step, interval, _last_step = struct.unpack(endian + "4i", header[8:24])
    delta_akma = struct.unpack(endian + "f", header[44:48])[0]
    charmm_version = struct.unpack(endian + "i", header[84:88])[0]
    if charmm_version <= 0:
        raise ValueError("Unsupported XPLOR-style DCD timing; provide explicit validated timing overrides.")
    if struct.unpack(endian + "i", header[88:92])[0] != 84:
        raise ValueError("DCD CHARMM timing record marker is invalid; timing is ambiguous.")
    if nset <= 0 or interval <= 0 or first_step < 0 or not math.isfinite(delta_akma) or delta_akma <= 0:
        raise ValueError("DCD header contains invalid frame count or timing; timing is ambiguous.")
    return {
        "source": "charmm_cord_header",
        "nset": nset,
        "first_step": first_step,
        "step_interval": interval,
        "delta_akma": float(delta_akma),
        "frame_interval_ps": float(delta_akma * _AKMA_PS * interval),
        "first_time_ps": float(first_step * delta_akma * _AKMA_PS),
        "charmm_version": charmm_version,
        "byte_order": "little" if endian == "<" else "big",
    }


def _finite_positive(value: object, label: str, *, allow_zero: bool = False) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must be a finite number.") from exc
    if not math.isfinite(result) or (result < 0 if allow_zero else result <= 0):
        qualifier = "non-negative" if allow_zero else "greater than zero"
        raise ValueError(f"{label} must be finite and {qualifier}.")
    return result


def kabsch_align(
    coordinates_angstrom: np.ndarray,
    reference_angstrom: np.ndarray,
    align_indices: np.ndarray,
) -> tuple[np.ndarray, float]:
    """Fit ``coordinates`` onto ``reference`` using a least-squares rigid fit.

    Coordinates are row vectors.  The returned RMSD is measured over the
    alignment selection, in Angstroms, after fitting.
    """
    xyz = np.asarray(coordinates_angstrom, dtype=np.float64)
    reference = np.asarray(reference_angstrom, dtype=np.float64)
    indices = np.asarray(align_indices, dtype=np.int64)
    if xyz.ndim != 2 or xyz.shape[1] != 3 or reference.shape != xyz.shape:
        raise ValueError("Frame coordinates and reference must have matching (atoms, 3) shapes.")
    if indices.ndim != 1 or len(indices) < 3 or np.any(indices < 0) or np.any(indices >= len(xyz)):
        raise ValueError("Protein alignment selection must contain at least three valid atoms.")
    if not np.isfinite(xyz).all() or not np.isfinite(reference).all():
        raise ValueError("Trajectory contains non-finite coordinates.")
    mobile = xyz[indices]
    target = reference[indices]
    mobile_center = mobile.mean(axis=0)
    target_center = target.mean(axis=0)
    covariance = (mobile - mobile_center).T @ (target - target_center)
    left, _singular, right_t = np.linalg.svd(covariance)
    correction = np.eye(3)
    correction[-1, -1] = np.linalg.det(left @ right_t)
    rotation = left @ correction @ right_t
    aligned = (xyz - mobile_center) @ rotation + target_center
    residual = aligned[indices] - target
    rmsd = float(np.sqrt(np.mean(np.sum(residual * residual, axis=1))))
    return aligned, rmsd


def residue_contact_metrics(
    protein_coordinates_angstrom: np.ndarray,
    ligand_coordinates_angstrom: np.ndarray,
    atom_residue_indices: np.ndarray,
    cutoff_angstrom: float,
) -> dict[int, float]:
    """Return minimum heavy-atom distance per protein residue for one frame."""
    protein = np.asarray(protein_coordinates_angstrom, dtype=np.float64)
    ligand = np.asarray(ligand_coordinates_angstrom, dtype=np.float64)
    residue_ids = np.asarray(atom_residue_indices, dtype=np.int64)
    cutoff = _finite_positive(cutoff_angstrom, "Contact cutoff")
    if protein.ndim != 2 or protein.shape[1] != 3 or not len(protein):
        raise ValueError("Protein contact selection contains no atoms.")
    if ligand.ndim != 2 or ligand.shape[1] != 3 or not len(ligand):
        raise ValueError("Ligand contact selection contains no atoms.")
    if len(residue_ids) != len(protein):
        raise ValueError("Protein residue mapping does not match selected protein atoms.")
    if not np.isfinite(protein).all() or not np.isfinite(ligand).all():
        raise ValueError("Trajectory contains non-finite coordinates.")
    try:
        from scipy.spatial import cKDTree
        distances, _indices = cKDTree(ligand).query(protein, k=1)
    except ImportError:
        # Simple fallback for minimal test/runtime environments; the bundled MD
        # image includes SciPy, so this is not used on production trajectories.
        distances = np.sqrt(((protein[:, None, :] - ligand[None, :, :]) ** 2).sum(axis=2)).min(axis=1)
    per_residue = {
        int(residue): float(distances[residue_ids == residue].min())
        for residue in np.unique(residue_ids)
    }
    return per_residue


def _sha256(path: Path) -> tuple[str, int]:
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as handle:
        while block := handle.read(4 * 1024 * 1024):
            digest.update(block)
            size += len(block)
    return digest.hexdigest(), size


def _residue_label(residue: Any) -> str:
    chain_id = str(getattr(residue.chain, "chain_id", "") or "")
    insertion = str(getattr(residue, "insertion_code", "") or "")
    return f"{chain_id}:{residue.resSeq}{insertion}:{residue.name}"


def _write_csv(path: Path, headers: list[str], rows: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def analyze(config: dict[str, Any]) -> dict[str, Any]:
    """Analyze a source trajectory in streaming chunks and write tabular outputs."""
    try:
        import mdtraj as md
    except ImportError as exc:
        raise RuntimeError("MDTraj is not available in the bundled MD runtime image.") from exc

    topology_path = Path(str(config["topology_path"]))
    trajectory_path = Path(str(config["trajectory_path"]))
    output_dir = Path(str(config["output_dir"]))
    if not topology_path.is_file() or not trajectory_path.is_file():
        raise FileNotFoundError("The mounted topology or trajectory is missing.")
    ligand_selection = str(config.get("ligand_selection") or "").strip()
    protein_selection = str(config.get("protein_selection") or "protein").strip()
    if not ligand_selection:
        raise ValueError("Enter an MDTraj atom selection for the ligand.")
    stride = config.get("stride", 1)
    if type(stride) is not int or stride < 1:
        raise ValueError("Stride must be a positive integer.")
    start_ns = _finite_positive(config.get("start_ns", 0.0), "Start time", allow_zero=True)
    end_ns = config.get("end_ns")
    end_ns = None if end_ns in (None, "") else _finite_positive(end_ns, "End time", allow_zero=True)
    if end_ns is not None and end_ns <= start_ns:
        raise ValueError("End time must be greater than start time.")
    cutoff = _finite_positive(config.get("contact_cutoff_angstrom", 4.5), "Contact cutoff")
    chunk_size = config.get("chunk_size", 100)
    if type(chunk_size) is not int or chunk_size < 1:
        raise ValueError("Chunk size must be a positive integer.")
    pbc_mode = str(config.get("pbc_mode") or "auto").lower()
    if pbc_mode not in {"auto", "none", "require"}:
        raise ValueError("PBC mode must be 'auto', 'none', or 'require'.")
    is_dcd = trajectory_path.suffix.casefold() == ".dcd"
    dcd_timing = _dcd_timing(
        trajectory_path,
        frame_interval_ps=config.get("dcd_frame_interval_ps"),
        first_time_ps=config.get("dcd_first_time_ps"),
    ) if is_dcd else None

    # MDTraj topology parsing and selection validation happen on a single frame;
    # this keeps memory bounded independently of trajectory length.
    first = md.load_frame(str(trajectory_path), 0, top=str(topology_path))
    atom_count = int(first.n_atoms)
    ligand_indices = np.asarray(first.topology.select(ligand_selection), dtype=np.int64)
    protein_indices = np.asarray(first.topology.select(protein_selection), dtype=np.int64)
    backbone_indices = np.asarray(first.topology.select(f"({protein_selection}) and backbone"), dtype=np.int64)
    if len(ligand_indices) == 0:
        raise ValueError(f"Ligand selection matched no atoms: {ligand_selection!r}.")
    if len(protein_indices) == 0:
        raise ValueError(f"Protein selection matched no atoms: {protein_selection!r}.")
    if np.intersect1d(ligand_indices, protein_indices).size:
        raise ValueError("Ligand and protein selections overlap; choose disjoint atom selections.")
    if len(backbone_indices) < 3:
        raise ValueError("Protein selection must contain at least three backbone atoms for alignment.")
    protein_heavy = np.asarray(
        [index for index in protein_indices if str(first.topology.atom(int(index)).element).lower() not in {"h", "hydrogen"}],
        dtype=np.int64,
    )
    ligand_heavy = np.asarray(
        [index for index in ligand_indices if str(first.topology.atom(int(index)).element).lower() not in {"h", "hydrogen"}],
        dtype=np.int64,
    )
    if len(protein_heavy) == 0 or len(ligand_heavy) == 0:
        raise ValueError("Protein and ligand contact selections must each contain heavy atoms.")
    residue_indices = np.asarray(
        [int(first.topology.atom(int(index)).residue.index) for index in protein_heavy],
        dtype=np.int64,
    )
    residue_labels = {
        int(residue.index): _residue_label(residue)
        for residue in first.topology.residues
        if int(residue.index) in set(residue_indices.tolist())
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    before_hashes = {}
    for name, path in (("topology", topology_path), ("trajectory", trajectory_path)):
        digest, size = _sha256(path)
        before_hashes[name] = {"path": str(config.get(f"{name}_source_path") or path), "sha256": digest, "size_bytes": size}

    reference = None
    first_selected = None
    contacts: dict[int, dict[str, float | int]] = {
        resid: {"minimum_distance_angstrom": math.inf, "distance_sum_angstrom": 0.0, "contact_frames": 0}
        for resid in residue_labels
    }
    backbone_rmsd_sum = ligand_rmsd_sum = contact_count_sum = 0.0
    backbone_rmsd_max = ligand_rmsd_max = 0.0
    total_frame_count = 0
    selected_frame_count = 0
    selected_start_ns = selected_end_ns = None
    pbc_observed = False
    pbc_applied = False
    pbc_notes: list[str] = []
    previous_time_ns: float | None = None
    raw_frame_index = 0
    frame_path = output_dir / "timeseries.csv"
    frame_headers = [
        "frame_index_after_stride", "time_ns", "protein_aligned_backbone_rmsd_angstrom",
        "ligand_rmsd_angstrom", "protein_ligand_minimum_distance_angstrom",
        "geometric_residue_contacts",
    ]
    frame_handle = frame_path.open("w", newline="", encoding="utf-8")
    frame_writer = csv.DictWriter(frame_handle, fieldnames=frame_headers)
    frame_writer.writeheader()

    # Read every frame in bounded chunks and apply stride ourselves. This keeps
    # source frame indices exact across chunk boundaries for DCD header timing.
    for chunk in md.iterload(str(trajectory_path), top=str(topology_path), chunk=chunk_size, stride=1):
        if int(chunk.n_atoms) != atom_count:
            raise ValueError("Trajectory atom count changed between chunks.")
        total_frame_count += int(chunk.n_frames)
        has_cell = chunk.unitcell_vectors is not None
        pbc_observed = pbc_observed or has_cell
        if pbc_mode == "require" and not has_cell:
            raise ValueError("PBC mode is 'require' but the trajectory has no periodic box vectors.")
        if has_cell and pbc_mode != "none":
            try:
                protein_set = set(int(index) for index in protein_indices)
                molecules = chunk.topology.find_molecules()
                protein_molecules = [
                    molecule
                    for molecule in molecules
                    if len(protein_set.intersection(int(atom.index) for atom in molecule)) > 1
                ]
                connected_protein_atoms = set().union(*(
                    protein_set.intersection(int(atom.index) for atom in molecule)
                    for molecule in protein_molecules
                )) if protein_molecules else set()
                if connected_protein_atoms != protein_set:
                    raise ValueError(
                        "Topology lacks bond connectivity for every selected protein atom; "
                        "cannot safely unwrap periodic coordinates."
                    )
                ligand_set = set(int(index) for index in ligand_indices)
                ligand_molecules = [
                    molecule for molecule in molecules
                    if len(ligand_set.intersection(int(atom.index) for atom in molecule)) > 1
                ]
                connected_ligand_atoms = set().union(*(
                    ligand_set.intersection(int(atom.index) for atom in molecule)
                    for molecule in ligand_molecules
                )) if ligand_molecules else set()
                if len(ligand_indices) > 1 and connected_ligand_atoms != ligand_set:
                    raise ValueError(
                        "Topology lacks bond connectivity for every selected ligand atom; "
                        "cannot safely unwrap periodic coordinates."
                    )
                # MDTraj expects each molecule as a set of TopologyAtom objects,
                # not integer indices. Its implementation dereferences atom.index.
                anchors = protein_molecules
                if not anchors:
                    raise ValueError("Topology has no bonded protein molecule for PBC imaging.")
                other_molecules = [molecule for molecule in molecules if molecule not in anchors]
                chunk.image_molecules(
                    inplace=True,
                    anchor_molecules=anchors,
                    other_molecules=other_molecules,
                )
                pbc_applied = True
            except Exception as exc:
                raise ValueError(f"Periodic imaging failed; coordinates were not analyzed: {exc}") from exc
        elif not has_cell and pbc_mode == "auto" and not pbc_notes:
            pbc_notes.append("Trajectory has no unit-cell vectors; used stored coordinates without PBC correction.")

        for frame_index in range(int(chunk.n_frames)):
            if is_dcd:
                assert dcd_timing is not None
                time_ps = dcd_timing["first_time_ps"] + raw_frame_index * dcd_timing["frame_interval_ps"]
            else:
                # XTC/TRR and other supported formats retain native timestamps
                # in MDTraj, in picoseconds.
                time_ps = float(chunk.time[frame_index])
            time_ns = time_ps / 1000.0
            if not math.isfinite(time_ns):
                raise ValueError("Trajectory contains a non-finite frame time.")
            if previous_time_ns is not None and time_ns <= previous_time_ns:
                raise ValueError("Trajectory frame times must be strictly increasing.")
            previous_time_ns = time_ns
            select_by_stride = raw_frame_index % stride == 0
            raw_frame_index += 1
            # DCD DELTA and XTC timestamps use single precision. Include a
            # nominal boundary frame despite its encoding roundoff.
            time_tolerance_ns = max(1e-9, 1e-7 * max(abs(time_ns), abs(start_ns), abs(end_ns or 0.0)))
            if not select_by_stride or time_ns + time_tolerance_ns < start_ns or (end_ns is not None and time_ns > end_ns + time_tolerance_ns):
                continue
            coordinates_angstrom = np.asarray(chunk.xyz[frame_index], dtype=np.float64) * 10.0
            if reference is None:
                reference = coordinates_angstrom.copy()
                selected_start_ns = time_ns
                first_selected = chunk[frame_index]
            aligned, protein_rmsd = kabsch_align(coordinates_angstrom, reference, backbone_indices)
            if selected_start_ns is None:
                selected_start_ns = time_ns
            ligand_delta = aligned[ligand_indices] - reference[ligand_indices]
            ligand_rmsd = float(np.sqrt(np.mean(np.sum(ligand_delta * ligand_delta, axis=1))))
            residue_distances = residue_contact_metrics(
                aligned[protein_heavy], aligned[ligand_heavy], residue_indices, cutoff
            )
            for residue_index, distance in residue_distances.items():
                aggregate = contacts[residue_index]
                aggregate["minimum_distance_angstrom"] = min(float(aggregate["minimum_distance_angstrom"]), distance)
                aggregate["distance_sum_angstrom"] = float(aggregate["distance_sum_angstrom"]) + distance
                if distance <= cutoff:
                    aggregate["contact_frames"] = int(aggregate["contact_frames"]) + 1
            row = {
                "frame_index_after_stride": selected_frame_count,
                "time_ns": time_ns,
                "protein_aligned_backbone_rmsd_angstrom": protein_rmsd,
                "ligand_rmsd_angstrom": ligand_rmsd,
                "protein_ligand_minimum_distance_angstrom": min(residue_distances.values()),
                "geometric_residue_contacts": sum(distance <= cutoff for distance in residue_distances.values()),
            }
            frame_writer.writerow(row)
            backbone_rmsd_sum += protein_rmsd
            ligand_rmsd_sum += ligand_rmsd
            backbone_rmsd_max = max(backbone_rmsd_max, protein_rmsd)
            ligand_rmsd_max = max(ligand_rmsd_max, ligand_rmsd)
            contact_count_sum += row["geometric_residue_contacts"]
            selected_frame_count += 1
            selected_end_ns = time_ns

    frame_handle.flush()
    frame_handle.close()

    if dcd_timing is not None and dcd_timing.get("nset") is not None and raw_frame_index != dcd_timing["nset"]:
        raise ValueError(
            f"DCD header frame count ({dcd_timing['nset']}) does not match decoded trajectory frame count ({raw_frame_index})."
        )

    if selected_frame_count == 0:
        raise ValueError("No trajectory frames fall inside the selected time range.")
    if not all(math.isfinite(value) for value in (
        backbone_rmsd_sum, ligand_rmsd_sum, backbone_rmsd_max, ligand_rmsd_max, contact_count_sum,
    )):
        raise ValueError("Analysis produced a non-finite scalar summary metric.")
    scalar_metrics = {
        "backbone_rmsd_mean_angstrom": backbone_rmsd_sum / selected_frame_count,
        "backbone_rmsd_max_angstrom": backbone_rmsd_max,
        "ligand_rmsd_mean_angstrom": ligand_rmsd_sum / selected_frame_count,
        "ligand_rmsd_max_angstrom": ligand_rmsd_max,
        "selected_frame_count": int(selected_frame_count),
        "contact_count_mean": contact_count_sum / selected_frame_count,
    }
    # Check inputs again after analysis so a source replacement/modification is
    # reported instead of yielding results tied to ambiguous bytes.
    after_hashes = {}
    for name, path in (("topology", topology_path), ("trajectory", trajectory_path)):
        digest, size = _sha256(path)
        after_hashes[name] = {"path": str(config.get(f"{name}_source_path") or path), "sha256": digest, "size_bytes": size}
    if before_hashes != after_hashes:
        raise RuntimeError("A source topology or trajectory changed while MD analysis was running.")

    contact_rows = []
    for residue_index, aggregate in contacts.items():
        contacting_frames = int(aggregate["contact_frames"])
        contact_rows.append({
            "residue": residue_labels[residue_index],
            "residue_index": residue_index,
            "minimum_distance_angstrom": aggregate["minimum_distance_angstrom"],
            "mean_minimum_distance_angstrom": float(aggregate["distance_sum_angstrom"]) / selected_frame_count,
            "contact_frames": contacting_frames,
            "selected_frames": selected_frame_count,
            "contact_occupancy_fraction": contacting_frames / selected_frame_count,
            "contact_occupancy_percent": 100.0 * contacting_frames / selected_frame_count,
            "contact_cutoff_angstrom": cutoff,
            "contact_definition": "minimum heavy-atom distance <= cutoff; geometric contact, not a hydrogen bond",
        })
    contact_path = output_dir / "residue_contacts.csv"
    _write_csv(contact_path, list(contact_rows[0].keys()) if contact_rows else [], contact_rows)

    representative_path = output_dir / "representative_frame.pdb"
    representative_saved = False
    if first_selected is not None and atom_count <= int(config.get("representative_atom_limit", 250_000)):
        try:
            representative = first_selected
            representative.save_pdb(str(representative_path))
            representative_saved = representative_path.is_file()
        except Exception as exc:
            pbc_notes.append(f"Representative PDB could not be saved: {exc}")
    summary = {
        "schema_version": 1,
        "success": True,
        "analysis": "protein-aligned RMSD and geometric ligand-contact occupancy",
        "units": {
            "time": "ns", "distance": "angstrom", "rmsd": "angstrom",
            "selected_frame_count": "frames",
            "contact_count_mean": "geometric heavy-atom residue contacts per analyzed frame",
        },
        "metrics": scalar_metrics,
        "source": {"topology": after_hashes["topology"], "trajectory": after_hashes["trajectory"]},
        "source_segments": config.get("source_segments") or [],
        "source_time": ({
            **dcd_timing,
            "time_formula": "(first_step + frame_index * step_interval) * delta_akma * 0.04888821 ps",
            "time_unit": "ps",
        } if dcd_timing is not None else {
            "source": "trajectory_native_timestamps",
            "time_unit": "ps",
        }),
        "trajectory_atom_count": atom_count,
        "trajectory_frames_seen": total_frame_count,
        "trajectory_frames_seen_after_stride": (total_frame_count + stride - 1) // stride,
        "selected_frame_count": selected_frame_count,
        "selected_time_range_ns": {"start": selected_start_ns, "end": selected_end_ns},
        "selection": {"ligand": ligand_selection, "protein": protein_selection, "alignment": f"({protein_selection}) and backbone"},
        "parameters": {
            "stride": stride, "chunk_size": chunk_size, "start_ns": start_ns,
            "end_ns": end_ns, "contact_cutoff_angstrom": cutoff,
            "pbc_mode_requested": pbc_mode, "periodic_unit_cell_present": pbc_observed,
            "pbc_imaging_applied": pbc_applied,
        },
        "contact_map_semantics": "Per-residue minimum distance among selected heavy atoms; occupancy is the fraction of selected frames with distance <= cutoff. This is a geometric contact measure, not a hydrogen-bond analysis.",
        "artifacts": {
            "timeseries_csv": "artifacts/md_analysis/timeseries.csv",
            "residue_contacts_csv": "artifacts/md_analysis/residue_contacts.csv",
            "representative_frame_pdb": "artifacts/md_analysis/representative_frame.pdb" if representative_saved else None,
        },
        "warnings": pbc_notes,
    }
    summary_path = output_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    args = parser.parse_args(argv)
    config_path = Path(args.config)
    output_path = Path(str(json.loads(config_path.read_text(encoding="utf-8")).get("output_dir") or ".")) / "summary.json"
    try:
        payload = analyze(json.loads(config_path.read_text(encoding="utf-8")))
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return 0
    except Exception as exc:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(json.dumps({"success": False, "error": str(exc), "traceback": traceback.format_exc()}, indent=2) + "\n", encoding="utf-8")
        raise


if __name__ == "__main__":
    raise SystemExit(main())
