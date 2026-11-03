"""Post-hoc, residue-linked tag and fusion suitability evidence.

This is a geometry-based advisory for one selected protein chain and a
user-selected interface partner group. It does not simulate a tag, linker,
expression, folding, or an experimental assay. The terminal exit direction,
clearance-cone, flexible-chain reach, tag-size separation, and evidence-first
presentation are adapted from BinderScope's ``binder_review.termini`` module
and ``docs/termini.md`` (MIT; Copyright (c) 2026 The University of Iowa).
Studio's implementation is independent of BindCraft run/design assumptions.
"""
from __future__ import annotations

import copy
import math
import shutil
import subprocess
import tempfile
import threading
from pathlib import Path
from typing import Any

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra
from scipy.spatial import cKDTree

METHOD = "Studio terminal geometry and clearance assessment (BinderScope termini-inspired)"
METHOD_VERSION = "studio.tag-fusion-suitability.v4"
CONTEXT_METHOD = "Studio tag/fusion clearance with explicitly placed saved-reference surroundings"
CONTEXT_VERSION = "studio.tag-fusion-context.v4"
CONTACT_CUTOFF_ANGSTROM = 4.5
EXIT_WINDOW = 6
TERMINAL_WINDOW = 5
SPHERE_CONE_DEGREES = 40.0
SPHERE_CONE_RINGS = 2
SPHERE_CONE_AZIMUTHS = 12
SPHERE_STANDOFF_ANGSTROM = 3.0
LINKER_REACH_ANGSTROM_PER_RESIDUE = 3.3
LINKER_STEPS = (0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 15)
HEAVY_ATOM_RADIUS_ANGSTROM = 1.7
FUSION_ENVELOPE_SLACK_ANGSTROM = 2.0
PATH_CLEARANCE_ANGSTROM = 3.0
PATH_SAMPLE_STEP_ANGSTROM = 1.0
FLEX_GRID_SPACING_ANGSTROM = 3.0
FLEX_GRID_HALF_WIDTH_ANGSTROM = 60.0
FLEX_START_RADIUS_ANGSTROM = 4.0
FLEX_START_CLEARANCE_ANGSTROM = 1.5
SASA_SPHERE_POINTS = 2400
_SASA_COMPAT_LOCK = threading.RLock()
TAGS = (
    {"key": "his6", "label": "His6", "kind": "flexible_affinity_tag", "residues": 6,
     "motif_residues": 3, "clearance_angstrom": 5.0},
    {"key": "twin_strep", "label": "Twin-Strep", "kind": "flexible_affinity_tag", "residues": 28,
     "motif_residues": 8, "clearance_angstrom": 10.0},
    {"key": "gfp_snap", "label": "GFP / SNAP-sized folded fusion", "kind": "globular_fusion",
     "radius_angstrom": 16.0},
    {"key": "fc_mbp", "label": "Fc / MBP-sized folded fusion", "kind": "globular_fusion",
     "radius_angstrom": 24.0},
)
_AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U", "PYL": "O",
}
_LOOP_SS = {"-", "T", "S", " "}


def _heavy_atoms(residue) -> list[Any]:
    atoms = []
    for atom in residue.get_atoms():
        element = str(getattr(atom, "element", "") or "").strip().upper()
        if element in {"H", "D"}:
            continue
        if not element:
            name = str(atom.get_name() or "").strip().upper()
            if name.startswith(("H", "D")):
                continue
        xyz = np.asarray(atom.get_coord(), dtype=float)
        if xyz.shape == (3,) and np.isfinite(xyz).all():
            atoms.append(atom)
    return atoms


def _xyz(atoms: list[Any]) -> np.ndarray:
    return np.asarray([atom.get_coord() for atom in atoms], dtype=float).reshape((-1, 3))


def _unit(vector: np.ndarray) -> np.ndarray | None:
    norm = float(np.linalg.norm(vector))
    if not math.isfinite(norm) or norm < 1e-8:
        return None
    return vector / norm


def _protein_residues(chain) -> list[Any]:
    from Bio.PDB.Polypeptide import is_aa
    return [residue for residue in chain.get_residues() if is_aa(residue, standard=False)]


def _ca_atom(residue):
    if "CA" not in residue:
        return None
    atom = residue["CA"]
    if getattr(atom, "is_disordered", lambda: 0)() == 2:
        alternatives = list(atom.child_dict.values())
        occupancies = [float(child.get_occupancy() or 0.0) for child in alternatives]
        if not occupancies or not all(math.isfinite(value) and value >= 0 for value in occupancies):
            return None
        maximum = max(occupancies)
        if occupancies.count(maximum) != 1:
            return None
        atom = alternatives[occupancies.index(maximum)]
    xyz = np.asarray(atom.get_coord(), dtype=float)
    return atom if xyz.shape == (3,) and np.isfinite(xyz).all() else None


def _residue_identity(chain_id: str, residue) -> dict[str, Any]:
    number = int(residue.id[1])
    insertion = str(residue.id[2] or "").strip()
    from mn_molecular_design_studio.complex_analysis_inspector import author_residue_selector
    return {
        "chain_id": str(chain_id), "author_number": number, "insertion_code": insertion,
        "residue_name": str(residue.resname).strip().upper(),
        "selector": author_residue_selector(str(chain_id), number, insertion),
    }


def _read_structure(path: Path):
    from mn_molecular_design_studio.complex_analysis_comparison import _read_biopython_structure
    structure = _read_biopython_structure(path)
    models = list(structure.get_models())
    if len(models) != 1:
        raise ValueError("Tag and Fusion assessment requires one selected single-model Database representation.")
    return structure, models[0]


def _dssp_data(path: Path, model) -> tuple[dict[tuple[str, Any], dict[str, Any]], dict[str, Any]]:
    try:
        from Bio.PDB import DSSP
        executable = shutil.which("mkdssp") or shutil.which("dssp")
        if not executable:
            return {}, {"status": "unavailable", "reason": "DSSP executable is not installed."}
        with tempfile.TemporaryDirectory(prefix="studio-dssp-") as temporary_directory:
            dssp_input = path
            suffix = path.name.casefold()
            is_mmcif = suffix.endswith((".cif", ".mmcif", ".cif.gz", ".mmcif.gz"))
            if not is_mmcif:
                try:
                    prefix = path.read_text(encoding="utf-8", errors="replace")[:256]
                except OSError:
                    prefix = ""
                if not prefix.startswith("HEADER"):
                    # DSSP 4 requires a PDB HEADER to distinguish PDB from mmCIF.
                    dssp_input = Path(temporary_directory) / "dssp_input.pdb"
                    dssp_input.write_text("HEADER    STUDIO TEMPORARY DSSP INPUT\n" + path.read_text(encoding="utf-8", errors="replace"),
                                          encoding="utf-8")
            if Path(executable).name.casefold() == "mkdssp":
                # mkdssp 4 defaults to mmCIF output; Bio.PDB 1.x parses DSSP
                # text, so request it explicitly rather than reporting an
                # empty map as a successful measurement.
                completed = subprocess.run(
                    [executable, "--output-format=dssp", str(dssp_input)],
                    capture_output=True, text=True, timeout=120, check=False,
                )
                if completed.returncode != 0 or not completed.stdout.strip():
                    reason = (completed.stderr or completed.stdout or "mkdssp produced no DSSP output.").strip()
                    return {}, {"status": "unavailable", "reason": reason[:500]}
                dssp_file = Path(temporary_directory) / "dssp_output.dssp"
                dssp_file.write_text(completed.stdout, encoding="utf-8")
                dssp = DSSP(model, str(dssp_file), acc_array="Sander", file_type="DSSP")
            else:
                dssp = DSSP(model, str(dssp_input), dssp=executable, acc_array="Sander")
        result = {}
        for key in dssp.keys():
            row = dssp[key]
            # Bio.PDB DSSP tuple: index, aa, secondary structure, relative ASA, ...
            result[(str(key[0]), key[1])] = {
                "secondary_structure": str(row[2]),
                "relative_accessibility": float(row[3]) if row[3] is not None else None,
            }
        version = "unknown"
        try:
            probe = subprocess.run([executable, "--version"], capture_output=True, text=True, timeout=5, check=False)
            version = (probe.stdout or probe.stderr or "unknown").strip().splitlines()[0][:120]
        except (OSError, subprocess.SubprocessError, IndexError):
            pass
        if not result:
            return {}, {"status": "unavailable", "reason": "DSSP completed but returned no residue assignments for the selected structure."}
        return result, {"status": "measured", "method": "DSSP relative accessibility (Sander scale) and secondary structure",
                        "executable": Path(executable).name, "version": version}
    except Exception as exc:
        return {}, {"status": "unavailable", "reason": f"DSSP could not evaluate this coordinate file: {type(exc).__name__}: {exc}"}


def _sasa_data(model) -> tuple[dict[tuple[str, Any], float], dict[str, Any]]:
    try:
        import numpy as np
        from Bio.PDB.SASA import ShrakeRupley
        # Bio.PDB 1.79 still refers to the removed NumPy ``np.int`` alias.
        # Supply that alias only for this locked call, then restore the module.
        # Newer Biopython versions do not enter this compatibility branch.
        with _SASA_COMPAT_LOCK:
            added_numpy_compat_alias = not hasattr(np, "int")
            if added_numpy_compat_alias:
                setattr(np, "int", int)
            try:
                ShrakeRupley(n_points=SASA_SPHERE_POINTS).compute(model, level="R")
            finally:
                if added_numpy_compat_alias:
                    delattr(np, "int")
        out = {}
        for chain in model:
            for residue in chain:
                value = getattr(residue, "sasa", None)
                if value is not None and math.isfinite(float(value)):
                    out[(str(chain.id), residue.id)] = float(value)
        return out, {"status": "measured", "method": "Bio.PDB ShrakeRupley", "sphere_points": SASA_SPHERE_POINTS,
                     "unit": "Å²", "scope": "selected complex coordinates",
                     "sampling_limit": "Finite sphere sampling is an approximation; residual orientation sensitivity is possible."}
    except (ImportError, OSError, ValueError, TypeError, RuntimeError) as exc:
        return {}, {"status": "unavailable", "reason": f"ShrakeRupley SASA could not be calculated: {type(exc).__name__}: {exc}"}


def _declared_terminal_coverage(path: Path, chain_id: str) -> dict[str, Any]:
    try:
        from mn_molecular_design_studio.complex_analysis_reference_context import reference_component_track
        track = reference_component_track(path, chain_id, [])
    except (ImportError, OSError, ValueError, TypeError, KeyError) as exc:
        return {"status": "unavailable", "reason": str(exc)}
    residues = list(track.get("residues") or [])
    if not residues:
        return {"status": "not_declared", "reason": track.get("reason") or "No declared polymer sequence was available."}
    first, last = residues[0], residues[-1]
    return {
        "status": "declared",
        "sequence_length": len(residues),
        "first_position_resolved": bool(first.get("author_number") is not None and first.get("reference_coordinates_present")),
        "last_position_resolved": bool(last.get("author_number") is not None and last.get("reference_coordinates_present")),
        "first_sequence_residue": first.get("amino_acid"),
        "last_sequence_residue": last.get("amino_acid"),
        "observed_sequence_first": residues[0].get("author_residue"),
        "first_status": first.get("status"), "last_status": last.get("status"),
        "numbering_source": track.get("numbering_source"),
    }


def _transform_xyz(coordinates: np.ndarray, matrix: Any) -> np.ndarray:
    transform = np.asarray(matrix, dtype=float)
    if transform.shape != (4, 4) or not np.isfinite(transform).all():
        raise ValueError("The saved reference-placement transform is not a finite 4×4 matrix.")
    if not np.allclose(transform[3], [0, 0, 0, 1], atol=1e-7):
        raise ValueError("The saved reference-placement transform has an invalid homogeneous row.")
    return (transform[:3, :3] @ coordinates.T).T + transform[:3, 3]


def _reference_obstacles(reference_path: Path, reference_chain: str,
                         mapped_pairs: list[dict[str, Any]]) -> tuple[np.ndarray, dict[str, Any]]:
    structure, model = _read_structure(reference_path)
    if reference_chain not in model:
        raise ValueError("The selected reference association chain is absent from its saved representation.")
    excluded = set()
    for pair in mapped_pairs:
        row = pair.get("reference") if isinstance(pair.get("reference"), dict) else {}
        try:
            excluded.add((int(row["number"]), str(row.get("insertion") or "").strip()))
        except (KeyError, TypeError, ValueError):
            raise ValueError("The saved placement result has incomplete reference residue mapping provenance.")
    coords = []
    reference_chain_residue_count = 0
    other_residue_count = 0
    for chain in model:
        for residue in chain:
            if str(chain.id) == str(reference_chain) and (int(residue.id[1]), str(residue.id[2] or "").strip()) in excluded:
                continue
            atoms = _heavy_atoms(residue)
            if not atoms:
                continue
            coords.extend(atom.get_coord() for atom in atoms)
            if str(chain.id) == str(reference_chain):
                reference_chain_residue_count += 1
            else:
                other_residue_count += 1
    array = np.asarray(coords, dtype=float).reshape((-1, 3))
    return array, {
        "reference_chain": str(reference_chain), "mapped_reference_residue_count_excluded": len(excluded),
        "unmapped_reference_chain_residue_count_included": reference_chain_residue_count,
        "other_reference_residue_count_included": other_residue_count,
        "reference_heavy_atom_count_included": len(array),
        "includes_only_resolved_coordinates": True,
    }


def _terminal_local_frame(residue, axis: np.ndarray) -> np.ndarray | None:
    """Build a proper molecular frame that rotates with the terminal residue."""
    u = _unit(axis)
    if u is None:
        return None
    ca = _ca_atom(residue)
    if ca is None:
        return None
    atom_vectors = []
    # Prefer the side-chain direction, then backbone geometry. These vectors
    # define only the azimuthal sampling frame; they do not change the exit axis.
    for name in ("CB", "N", "C", "O"):
        if name in residue:
            atom = residue[name]
            if getattr(atom, "is_disordered", lambda: 0)() == 2:
                alternatives = list(atom.child_dict.values())
                occupancies = [float(child.get_occupancy() or 0.0) for child in alternatives]
                if not occupancies or not all(math.isfinite(value) and value >= 0 for value in occupancies):
                    continue
                maximum = max(occupancies)
                if occupancies.count(maximum) != 1:
                    continue
                atom = alternatives[occupancies.index(maximum)]
            xyz = np.asarray(atom.get_coord(), dtype=float)
            if xyz.shape == (3,) and np.isfinite(xyz).all():
                atom_vectors.append(xyz - np.asarray(ca.get_coord(), dtype=float))
    for vector in atom_vectors:
        transverse = vector - float(np.dot(vector, u)) * u
        e1 = _unit(transverse)
        if e1 is None:
            continue
        e2 = _unit(np.cross(u, e1))
        if e2 is None:
            continue
        e1 = _unit(np.cross(e2, u))
        if e1 is None:
            continue
        frame = np.column_stack((e1, e2, u))
        if np.isfinite(frame).all() and np.linalg.det(frame) > 0.999999:
            return frame
    return None


def _cone_directions(axis: np.ndarray, frame: np.ndarray | None = None) -> np.ndarray:
    u = _unit(axis)
    if u is None:
        return np.empty((0, 3), dtype=float)
    if frame is None:
        # This fallback is retained for callers without residue geometry.
        helper = np.array([1.0, 0.0, 0.0]) if abs(float(u[0])) < 0.9 else np.array([0.0, 1.0, 0.0])
        e1 = _unit(np.cross(u, helper))
        e2 = np.cross(u, e1)
    else:
        frame = np.asarray(frame, dtype=float)
        if frame.shape != (3, 3) or not np.isfinite(frame).all() or not np.allclose(frame[:, 2], u, atol=1e-6):
            return np.empty((0, 3), dtype=float)
        e1, e2 = frame[:, 0], frame[:, 1]
    result = [u]
    for ring in range(1, SPHERE_CONE_RINGS + 1):
        theta = math.radians(SPHERE_CONE_DEGREES * ring / SPHERE_CONE_RINGS)
        for index in range(SPHERE_CONE_AZIMUTHS):
            phi = 2.0 * math.pi * index / SPHERE_CONE_AZIMUTHS
            result.append(math.cos(theta) * u + math.sin(theta) * (math.cos(phi) * e1 + math.sin(phi) * e2))
    return np.asarray(result, dtype=float)


def _flexible_escape(origin: np.ndarray, obstacles: np.ndarray, clearance: float,
                     frame: np.ndarray | None = None) -> float | None:
    """Shortest grid path to open solvent; not a molecular/tag conformation."""
    if not len(obstacles):
        return 0.0
    half = int(round(FLEX_GRID_HALF_WIDTH_ANGSTROM / FLEX_GRID_SPACING_ANGSTROM))
    side = 2 * half + 1
    axis = np.arange(-half, half + 1, dtype=float) * FLEX_GRID_SPACING_ANGSTROM
    relative = np.stack(np.meshgrid(axis, axis, axis, indexing="ij"), axis=-1).reshape(-1, 3)
    if frame is None:
        # Backward-compatible deterministic frame for direct utility callers.
        u = np.array([0.0, 0.0, 1.0])
        frame = np.eye(3)
    else:
        frame = np.asarray(frame, dtype=float)
        if frame.shape != (3, 3) or not np.isfinite(frame).all() or not np.allclose(frame.T @ frame, np.eye(3), atol=1e-6) or np.linalg.det(frame) < 0.999999:
            raise ValueError("The terminal-local clearance frame must be a proper orthonormal rotation matrix.")
    obstacle_local = (np.asarray(obstacles, dtype=float) - np.asarray(origin, dtype=float)) @ frame
    obstacle_tree = cKDTree(obstacle_local)
    nearest = obstacle_tree.query(relative, k=1)[0]
    near_start = np.linalg.norm(relative, axis=1) <= FLEX_START_RADIUS_ANGSTROM
    free = (nearest >= PATH_CLEARANCE_ANGSTROM) | (near_start & (nearest >= FLEX_START_CLEARANCE_ANGSTROM))
    center = half * side * side + half * side + half
    free[center] = True
    indices = np.arange(side ** 3, dtype=np.int64).reshape(side, side, side)
    edges_a, edges_b, weights = [], [], []
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                if (dx, dy, dz) <= (0, 0, 0):
                    continue
                source_slices = tuple(slice(max(0, -d), side - max(0, d)) for d in (dx, dy, dz))
                target_slices = tuple(slice(max(0, d), side - max(0, -d)) for d in (dx, dy, dz))
                valid = free.reshape(side, side, side)[source_slices] & free.reshape(side, side, side)[target_slices]
                source_ids = indices[source_slices][valid]
                target_ids = indices[target_slices][valid]
                edges_a.append(source_ids)
                edges_b.append(target_ids)
                weights.append(np.full(len(source_ids), FLEX_GRID_SPACING_ANGSTROM * math.sqrt(dx * dx + dy * dy + dz * dz)))
    if not edges_a or not any(len(row) for row in edges_a):
        return None
    graph = coo_matrix((np.concatenate(weights), (np.concatenate(edges_a), np.concatenate(edges_b))),
                       shape=(side ** 3, side ** 3)).tocsr()
    distance = dijkstra(graph, directed=False, indices=center, return_predecessors=False,
                        limit=max(tag["residues"] for tag in TAGS if tag["kind"] == "flexible_affinity_tag")
                        * LINKER_REACH_ANGSTROM_PER_RESIDUE + max(LINKER_STEPS) * LINKER_REACH_ANGSTROM_PER_RESIDUE)
    open_nodes = (nearest >= float(clearance)) & np.isfinite(distance)
    return float(np.min(distance[open_nodes])) if open_nodes.any() else None


def _sphere_fits(origin: np.ndarray, directions: np.ndarray, obstacles: np.ndarray,
                 radius: float, linker_residues: int) -> bool:
    if not len(directions) or not len(obstacles):
        return bool(len(directions))
    tree = cKDTree(obstacles)
    reach = radius + SPHERE_STANDOFF_ANGSTROM + LINKER_REACH_ANGSTROM_PER_RESIDUE * linker_residues
    surface_clearance = radius + HEAVY_ATOM_RADIUS_ANGSTROM - FUSION_ENVELOPE_SLACK_ANGSTROM
    step_count = max(1, int(math.ceil(max(0.0, reach - radius) / PATH_SAMPLE_STEP_ANGSTROM)))
    for direction in directions:
        center = origin + direction * reach
        if float(tree.query(center)[0]) < surface_clearance:
            continue
        path_end = max(0.0, reach - radius)
        samples = origin[None, :] + direction[None, :] * np.linspace(2.0, path_end, max(2, step_count))[:, None]
        if np.all(tree.query(samples)[0] >= PATH_CLEARANCE_ANGSTROM):
            return True
    return False


def _tag_clearances(origin: np.ndarray, direction: np.ndarray, obstacles: np.ndarray,
                    frame: np.ndarray | None = None) -> dict[str, Any]:
    directions = _cone_directions(direction, frame)
    result = {}
    for tag in TAGS:
        if tag["kind"] == "flexible_affinity_tag":
            base_reach = (int(tag["residues"]) - int(tag["motif_residues"])) * LINKER_REACH_ANGSTROM_PER_RESIDUE
            escape = _flexible_escape(origin, obstacles, float(tag["clearance_angstrom"]), frame)
            linker = next((step for step in LINKER_STEPS
                           if escape is not None and escape <= base_reach + step * LINKER_REACH_ANGSTROM_PER_RESIDUE), None)
            result[tag["key"]] = {
                "label": tag["label"], "kind": tag["kind"],
                "status": "no_clear_path_within_linker_search" if linker is None else (
                    "clearance_path_without_added_linker" if linker == 0 else "clearance_path_with_linker"),
                "minimum_linker_residues_tested": linker,
                "shortest_grid_path_angstrom": None if escape is None else round(escape, 2),
                "clearance_radius_angstrom": tag["clearance_angstrom"],
                "tag_length_residues": tag["residues"],
                "assumed_motif_residues_requiring_open_solvent": tag["motif_residues"],
                "reach_model": "(tag residues − exposed motif residues + linker residues) × 3.3 Å/residue",
            }
        else:
            linker = next((step for step in LINKER_STEPS
                           if _sphere_fits(origin, directions, obstacles, float(tag["radius_angstrom"]), step)), None)
            result[tag["key"]] = {
                "label": tag["label"], "kind": tag["kind"],
                "status": "clearance_cone_without_added_linker" if linker == 0 else (
                    "clearance_cone_with_linker" if linker is not None else "no_clearance_cone_within_linker_search"),
                "minimum_linker_residues_tested": linker,
                "effective_spherical_radius_angstrom": tag["radius_angstrom"],
                "cone_half_angle_degrees": SPHERE_CONE_DEGREES,
            }
    return result


def _endpoint_assessment(*, end: str, chain_id: str, residues: list[Any], ca_atoms: list[Any],
                         target_tree: cKDTree, target_xyz: np.ndarray, obstacle_xyz: np.ndarray,
                         dssp: dict, sasa: dict, plddt_enabled: bool,
                         boundary_resolved: bool, termini_confirmed: bool,
                         interface_residue_keys: set[Any], binder_centroid: np.ndarray,
                         interface_center: np.ndarray | None) -> dict[str, Any]:
    index = 0 if end == "N" else len(residues) - 1
    residue = residues[index]
    identity = _residue_identity(chain_id, residue)
    row: dict[str, Any] = {
        "terminus": end, "terminal_residue": identity,
        "assessment_status": "measured", "distance_to_selected_interface_angstrom": None,
        "distance_definition": "terminal Cα to nearest heavy atom in the selected interface-partner chain group",
        "direction": None, "direction_definition": "terminal Cα minus mean of up to six inward Cα positions",
        "aim_toward_partner_degrees": None, "interface_side_cosine": None,
        "contact_at_terminal_residue": residue.id in interface_residue_keys,
        "terminal_residue_sasa_angstrom2": sasa.get((str(chain_id), residue.id)),
        "relative_sasa_sander": (dssp.get((str(chain_id), residue.id)) or {}).get("relative_accessibility"),
        "secondary_structure": (dssp.get((str(chain_id), residue.id)) or {}).get("secondary_structure"),
        "terminal_window_secondary_structure": [], "terminal_window_mean_plddt": None,
        "terminal_window_mean_b_factor": None, "flexibility_evidence": "not_assessed",
        "tag_clearance": {},
    }
    if not boundary_resolved:
        row["assessment_status"] = "not_assessed_incomplete_terminal_structure"
        return row
    if not termini_confirmed:
        row["assessment_status"] = "provisional_observed_chain_endpoint"
    if ca_atoms[index] is None:
        row["assessment_status"] = "not_assessed_terminal_ca_unresolved"
        return row
    ca = np.asarray(ca_atoms[index].get_coord(), dtype=float)
    inner = ca_atoms[index + 1:index + 1 + EXIT_WINDOW] if end == "N" else ca_atoms[max(0, index - EXIT_WINDOW):index]
    if len(inner) < 3 or any(atom is None for atom in inner):
        row["assessment_status"] = "not_assessed_insufficient_direction_coordinates"
        return row
    direction = _unit(ca - np.mean([atom.get_coord() for atom in inner], axis=0))
    if direction is None:
        row["assessment_status"] = "not_assessed_degenerate_terminal_direction"
        return row
    local_frame = _terminal_local_frame(residue, direction)
    nearest_distance, nearest_index = target_tree.query(ca)
    toward = _unit(target_xyz[int(nearest_index)] - ca)
    row["distance_to_selected_interface_angstrom"] = round(float(nearest_distance), 3)
    row["direction"] = [round(float(value), 4) for value in direction]
    if toward is not None:
        row["aim_toward_partner_degrees"] = round(math.degrees(math.acos(float(np.clip(np.dot(direction, toward), -1, 1)))), 2)
    if interface_center is not None:
        away_a = _unit(ca - binder_centroid)
        away_b = _unit(binder_centroid - interface_center)
        if away_a is not None and away_b is not None:
            row["interface_side_cosine"] = round(float(np.dot(away_a, away_b)), 4)
    terminal_window = residues[:TERMINAL_WINDOW] if end == "N" else residues[-TERMINAL_WINDOW:][::-1]
    ss_rows = [dssp.get((str(chain_id), item.id)) for item in terminal_window]
    row["terminal_window_secondary_structure"] = [None if not value else value.get("secondary_structure") for value in ss_rows]
    if plddt_enabled:
        b_factors = []
        for terminal_residue in terminal_window:
            atom = _ca_atom(terminal_residue)
            if atom is not None:
                value = float(atom.get_bfactor())
                if not math.isfinite(value) or not 0 <= value <= 100:
                    raise ValueError("B-factor values declared as pLDDT must be finite and within 0–100.")
                b_factors.append(value)
        if b_factors:
            row["terminal_window_mean_plddt"] = round(float(np.mean(b_factors)), 2)
            row["terminal_window_mean_b_factor"] = round(float(np.mean(b_factors)), 2)
    else:
        values = [float(_ca_atom(item).get_bfactor()) for item in terminal_window if _ca_atom(item) is not None]
        if values and all(math.isfinite(value) for value in values):
            row["terminal_window_mean_b_factor"] = round(float(np.mean(values)), 2)
    loop_count = sum(bool(value and value.get("secondary_structure") in _LOOP_SS) for value in ss_rows)
    if all(value is None for value in ss_rows):
        row["flexibility_evidence"] = "not_assessed_secondary_structure_unavailable"
    elif row["terminal_window_mean_plddt"] is not None and row["terminal_window_mean_plddt"] < 70 and loop_count >= 3:
        row["flexibility_evidence"] = "loop_rich_low_declared_plddt_possible_flexibility"
    elif loop_count >= 3:
        row["flexibility_evidence"] = "loop_rich_static_structure; dynamics not measured"
    else:
        row["flexibility_evidence"] = "no_loop_rich_low_confidence_signal_in_available_static_evidence"
    row["terminal_window_loop_residue_count"] = loop_count
    if not termini_confirmed:
        row["clearance_note"] = "Clearance is reported for the observed file-chain endpoint only; tag suitability is not concluded until construct-end completeness is confirmed."
        return row
    if local_frame is None:
        row["clearance_note"] = "Terminal-local clearance sampling frame could not be defined from resolved terminal atoms; tag/fusion clearance was not assessed."
        return row
    # Exclude only the terminal residue from the component's own-fold obstacle
    # set. All other chains and components in the exact saved complex remain.
    row["tag_clearance"] = _tag_clearances(ca, direction, obstacle_xyz, local_frame)
    # Keep sampled geometric clearance distinct from suitability. BinderScope's
    # termini method treats a terminal interface contact as disqualifying for a
    # folded fusion even where a coarse sphere sample finds an open exit cone.
    # Preserve that geometric result, but carry the attachment-site advisory
    # separately so the Inspector cannot present it as an unqualified pass.
    for tag in TAGS:
        clearance = row["tag_clearance"].get(tag["key"])
        if not isinstance(clearance, dict):
            continue
        status = str(clearance.get("status") or "")
        geometric_pass = status in {
            "clearance_path_without_added_linker", "clearance_path_with_linker",
            "clearance_cone_without_added_linker", "clearance_cone_with_linker",
        }
        if tag["kind"] == "globular_fusion" and row["contact_at_terminal_residue"]:
            clearance["suitability_advisory"] = {
                "status": "not_recommended_terminal_interface_contact",
                "reason": "The terminal residue participates in the selected interface. The sampled geometric envelope may clear, but attachment at this residue can perturb the interface.",
            }
        elif geometric_pass:
            clearance["suitability_advisory"] = {
                "status": "clearance_only_no_terminal_interface_contact",
                "reason": "The sampled geometric clearance passed and the terminal residue is not in the selected interface; this remains a static advisory.",
            }
        else:
            clearance["suitability_advisory"] = {
                "status": "not_recommended_clearance_not_demonstrated",
                "reason": "The sampled geometric clearance did not pass within the tested linker search.",
            }
    return row


def calculate_tag_fusion_suitability(
    structure_path: Path | str,
    *,
    component_chain: str,
    interface_chains: list[str],
    model_id: str,
    source_identity: dict[str, Any],
    interface_contact_cutoff: float = CONTACT_CUTOFF_ANGSTROM,
    termini_confirmed: bool = False,
    component_format: str = "single_chain_linear_protein",
    b_factor_interpretation: str = "uninterpreted",
    saved_confidence_context: dict[str, Any] | None = None,
    reference_context: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Evaluate one chain's observed N/C termini against selected interface chains.

    ``reference_context`` is populated only by the authorized Database adapter
    for the separately registered context workflow. It contains a digest-bound
    explicit association and an already-saved Reference Placement Evaluation;
    this function never fits or infers a reference.
    """
    from mn_molecular_design_studio.evaluations import file_digest
    from mn_molecular_design_studio.complex_analysis_placement import PLACEMENT_METHOD, PLACEMENT_VERSION

    path = Path(structure_path).expanduser().resolve(strict=True)
    source_digest = file_digest(path)
    expected_digest = str(source_identity.get("structure_digest") or "")
    if not expected_digest or source_digest != expected_digest:
        raise ValueError("The selected source structure digest does not match the authorized Database representation.")
    if str(model_id or "") != str(source_identity.get("model_id") or ""):
        raise ValueError("The selected model ID does not match the authorized Database representation.")
    component_chain = str(component_chain or "").strip()
    interface_chains = list(dict.fromkeys(str(chain).strip() for chain in interface_chains if str(chain).strip()))
    if not component_chain or not interface_chains or component_chain in interface_chains:
        raise ValueError("Choose one protein component and one or more disjoint interface-partner chains.")
    if not math.isfinite(float(interface_contact_cutoff)) or not 1.0 <= float(interface_contact_cutoff) <= 10.0:
        raise ValueError("The interface contact cutoff must be between 1 and 10 Å.")
    if component_format not in {"single_chain_linear_protein", "cyclic_peptide", "multichain_construct"}:
        raise ValueError("Choose a supported component-format declaration.")
    if b_factor_interpretation not in {"uninterpreted", "plddt_0_100"}:
        raise ValueError("Choose a supported B-factor interpretation.")
    _structure, model = _read_structure(path)
    if component_chain not in model or any(chain not in model for chain in interface_chains):
        raise ValueError("A selected protein chain is absent from the chosen model.")
    component_residues = _protein_residues(model[component_chain])
    partner_residues = [residue for chain in interface_chains for residue in _protein_residues(model[chain])]
    if not component_residues:
        raise ValueError("The selected component chain contains no supported amino-acid residues.")
    if not partner_residues:
        raise ValueError("The selected interface partner group contains no supported amino-acid residues.")

    identity = {key: str(source_identity.get(key) or "") for key in (
        "molecule_id", "structure_id", "group_id", "structure_digest", "model_id", "evaluation_id",
    )}
    params = {
        "component_chain": component_chain, "interface_chains": interface_chains,
        "interface_contact_cutoff_angstrom": float(interface_contact_cutoff),
        "termini_confirmed": bool(termini_confirmed), "component_format": component_format,
        "b_factor_interpretation": b_factor_interpretation,
    }
    scope = "intrinsic_complex" if reference_context is None else "validated_reference_context"
    method, version = (METHOD, METHOD_VERSION) if reference_context is None else (CONTEXT_METHOD, CONTEXT_VERSION)
    base = {
        "schema": "studio.tag_fusion_suitability.v1", "method": method, "method_version": version,
        "adaptation_attribution": "BinderScope binder_review.termini (MIT; The University of Iowa, 2026); see docs/third-party-notices.md.",
        "assessment_scope": scope, "source": identity, "parameters": params,
        "tag_fusion_models": [dict(item) for item in TAGS],
        "method_parameters": {
            "contact_definition": "inter-chain heavy-atom distance cutoff",
            "exit_direction_window_ca_residues": EXIT_WINDOW,
            "terminal_evidence_window_residues": TERMINAL_WINDOW,
            "sphere_cone_half_angle_degrees": SPHERE_CONE_DEGREES,
            "sphere_cone_ring_count": SPHERE_CONE_RINGS,
            "sphere_cone_azimuth_samples_per_ring": SPHERE_CONE_AZIMUTHS,
            "sphere_standoff_angstrom": SPHERE_STANDOFF_ANGSTROM,
            "linker_reach_angstrom_per_residue": LINKER_REACH_ANGSTROM_PER_RESIDUE,
            "linker_lengths_tested_residues": list(LINKER_STEPS),
            "heavy_atom_radius_angstrom": HEAVY_ATOM_RADIUS_ANGSTROM,
            "fusion_envelope_slack_angstrom": FUSION_ENVELOPE_SLACK_ANGSTROM,
            "path_clearance_angstrom": PATH_CLEARANCE_ANGSTROM,
            "path_line_sample_step_angstrom": PATH_SAMPLE_STEP_ANGSTROM,
            "flexible_path_grid_spacing_angstrom": FLEX_GRID_SPACING_ANGSTROM,
            "flexible_path_grid_half_width_angstrom": FLEX_GRID_HALF_WIDTH_ANGSTROM,
            "clearance_sampling_frame": "terminal-residue local molecular frame, covariant to rigid-body rotation",
            "flexible_tag_reach_assumption": "(tag length − exposed motif length + linker length) × 3.3 Å/residue",
            "surface_accessibility": f"Bio.PDB ShrakeRupley, {SASA_SPHERE_POINTS} sphere points, absolute residue SASA in Å²; DSSP Sander relative accessibility when available",
            "dynamics_simulation": False,
        },
        "source_structure_unchanged": True,
        "saved_confidence_context": {},
        "interpretation_limitations": [
            "This is an advisory analysis of one saved coordinate model, not experimental proof of tag accessibility, expression, folding, binding, or assay performance.",
            "Flexible tags are evaluated by a shortest-clearance-path reach estimate; polymer conformations, self-avoidance, entropy, chemical attachment geometry, and tag-receptor binding are not simulated.",
            "GFP/SNAP and Fc/MBP are approximated as spheres in a fixed exit cone; their actual shapes, orientations, linkers, and domain-specific interactions are not modeled.",
            "Static secondary structure and confidence are evidence about the supplied model, not a measurement of molecular dynamics or flexibility.",
            "A reference-context run checks only resolved coordinates in the explicitly associated saved reference; it does not infer missing atoms, occupancy, biological assembly, glycans, or accessibility.",
        ],
    }
    if saved_confidence_context:
        if (str(saved_confidence_context.get("evaluation_id") or "") != identity["evaluation_id"]
                or str(saved_confidence_context.get("structure_digest") or "") != source_digest
                or str(saved_confidence_context.get("model_id") or "") != identity["model_id"]):
            raise ValueError("The saved confidence Evaluation does not match the exact selected Database representation and model.")
        base["saved_confidence_context"] = {
            **saved_confidence_context,
            "used_for_terminal_grading": False,
            "scope": "saved Evaluation-wide confidence outputs; not localized to terminal residues",
        }
    reference_obstacle_xyz = np.empty((0, 3), dtype=float)
    if reference_context is not None:
        base["reference_context"] = {
            key: value for key, value in reference_context.items()
            if key not in {"structure_path", "placement_result"}
        }
        placement = reference_context.get("placement_result") if isinstance(reference_context.get("placement_result"), dict) else {}
        if (placement.get("schema") != "studio.complex_analysis_reference_placement.v1"
                or placement.get("method") != PLACEMENT_METHOD
                or placement.get("method_version") != PLACEMENT_VERSION):
            raise ValueError("The reference context requires a matching saved Reference Placement Evaluation.")
        ref_source = placement.get("source") if isinstance(placement.get("source"), dict) else {}
        ref_reference = placement.get("reference") if isinstance(placement.get("reference"), dict) else {}
        expected_context_id = str(reference_context.get("reference_context_id") or reference_context.get("id") or "")
        if (str(ref_source.get("molecule_id") or "") != identity["molecule_id"]
                or str(ref_source.get("structure_id") or "") != identity["structure_id"]
                or str(ref_source.get("group_id") or "") != identity["group_id"]
                or str(ref_source.get("structure_digest") or "") != source_digest
                or str(ref_source.get("model_id") or "") != identity["model_id"]
                or str(ref_reference.get("molecule_id") or "") != str(reference_context.get("reference_molecule_id") or "")
                or str(ref_reference.get("structure_id") or "") != str(reference_context.get("reference_structure_id") or "")
                or str(ref_reference.get("group_id") or "") != str(reference_context.get("reference_group_id") or "")
                or str(ref_reference.get("structure_digest") or "") != str(reference_context.get("reference_structure_digest") or "")
                or str(ref_reference.get("model_id") or "") != str(reference_context.get("reference_model_id") or "")
                or str(ref_reference.get("chain_id") or "") != str(reference_context.get("reference_component_chain") or "")
                or str(ref_source.get("chain_id") or "") != str(reference_context.get("source_component_chain") or "")
                or placement.get("transformation_direction") != "selected complex/source coordinates -> associated reference coordinates"
                or str(placement.get("reference_context_id") or "") != expected_context_id):
            raise ValueError("The saved placement Evaluation does not match the exact source/reference association and digests.")
        matrix = placement.get("matrix_4x4")
        # Placement stores a validated source -> reference transform. Apply it
        # transiently to the complete complex for the context calculation.
        from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer.dataclasses import validate_rigid_transform_matrix
        matrix = validate_rigid_transform_matrix(matrix)
        for chain in model:
            for atom in chain.get_atoms():
                atom.set_coord(_transform_xyz(np.asarray(atom.get_coord(), dtype=float).reshape((1, 3)), matrix)[0])
        reference_path = Path(str(reference_context.get("structure_path") or "")).expanduser().resolve(strict=True)
        if file_digest(reference_path) != str(reference_context.get("reference_structure_digest") or ""):
            raise ValueError("The associated reference structure digest changed after Database authorization.")
        reference_chain = str(reference_context.get("reference_component_chain") or "")
        mapped_pairs = list(placement.get("mapped_residue_pairs") or [])
        reference_obstacle_xyz, detail = _reference_obstacles(reference_path, reference_chain, mapped_pairs)
        base["reference_context"].update(detail)
        base["reference_context"].update({
            "reference_structure_id": str(reference_context.get("reference_structure_id") or ""),
            "reference_structure_digest": str(reference_context.get("reference_structure_digest") or ""),
            "reference_model_id": str(reference_context.get("reference_model_id") or ""),
            "reference_component_chain": reference_chain,
            "placement_evaluation_id": str(reference_context.get("placement_evaluation_id") or ""),
            "placement_run_id": str(reference_context.get("placement_run_id") or ""),
            "placement_method": PLACEMENT_METHOD,
            "placement_method_version": PLACEMENT_VERSION,
            "mapping_method": placement.get("mapping_method"),
            "mapping_residue_pair_count": placement.get("aligned_pair_count"),
            "placement_rmsd_angstrom": placement.get("rmsd_angstrom"),
            "transform_direction": placement.get("transformation_direction"),
            "explicit_branch_entity_annotations": [],
        })
        try:
            from mn_molecular_design_studio.complex_analysis_reference_context import explicit_reference_components
            base["reference_context"]["explicit_branch_entity_annotations"] = explicit_reference_components(reference_path)
        except (ImportError, OSError, ValueError, TypeError):
            pass

    component_residues = _protein_residues(model[component_chain])
    partner_residues = [residue for chain in interface_chains for residue in _protein_residues(model[chain])]
    target_atoms = [atom for residue in partner_residues for atom in _heavy_atoms(residue)]
    if not target_atoms:
        raise ValueError("The selected interface partner chains contain no usable heavy atoms.")
    target_xyz = _xyz(target_atoms)
    target_tree = cKDTree(target_xyz)
    contacts_by_component_residue: dict[Any, float] = {}
    partner_interface = set()
    contact_cutoff = float(interface_contact_cutoff)
    partner_chain_residues = {chain: _protein_residues(model[chain]) for chain in interface_chains}
    for residue in component_residues:
        atoms = _heavy_atoms(residue)
        if not atoms:
            continue
        distances = target_tree.query(_xyz(atoms))[0]
        minimum = float(np.min(distances)) if len(distances) else math.inf
        if minimum <= contact_cutoff:
            contacts_by_component_residue[residue.id] = minimum
    component_heavy_xyz = _xyz([atom for residue in component_residues for atom in _heavy_atoms(residue)])
    component_heavy_tree = cKDTree(component_heavy_xyz) if len(component_heavy_xyz) else None
    for chain, residues in partner_chain_residues.items():
        if component_heavy_tree is None:
            continue
        for residue in residues:
            atoms = _heavy_atoms(residue)
            if atoms and float(np.min(component_heavy_tree.query(_xyz(atoms))[0])) <= contact_cutoff:
                partner_interface.add((chain, residue.id))
    component_interface_xyz = []
    partner_interface_xyz = []
    for residue in component_residues:
        if residue.id in contacts_by_component_residue:
            component_interface_xyz.extend(atom.get_coord() for atom in _heavy_atoms(residue))
    for chain, residue_id in partner_interface:
        partner_interface_xyz.extend(atom.get_coord() for atom in _heavy_atoms(model[chain][residue_id]))
    interface_center = None
    if component_interface_xyz and partner_interface_xyz:
        interface_center = (np.mean(component_interface_xyz, axis=0) + np.mean(partner_interface_xyz, axis=0)) / 2.0
    ca_atoms = [_ca_atom(residue) for residue in component_residues]
    component_xyz = component_heavy_xyz
    other_xyz = []
    for chain in model:
        for residue in chain:
            is_component_residue = str(chain.id) == component_chain and residue in component_residues
            is_interface_residue = (str(chain.id), residue.id) in partner_interface
            if is_component_residue or (str(chain.id) in interface_chains and residue in partner_residues):
                continue
            other_xyz.extend(atom.get_coord() for atom in _heavy_atoms(residue))
    sasa, sasa_provenance = _sasa_data(model)
    dssp, dssp_provenance = _dssp_data(path, model)
    declared = _declared_terminal_coverage(path, component_chain)
    terminal_flags = {
        "N": declared.get("first_position_resolved") if declared.get("status") == "declared" else None,
        "C": declared.get("last_position_resolved") if declared.get("status") == "declared" else None,
    }
    if component_format in {"cyclic_peptide", "multichain_construct"}:
        reason = ("Cyclic peptides do not have free N/C termini for this analysis." if component_format == "cyclic_peptide"
                  else "A multichain construct has more than two molecular ends; select and analyze a justified single-chain component separately.")
        return {**base, "assessment_status": "not_applicable", "not_applicable_reason": reason,
                "terminal_geometry": {}, "interface_contact_residue_count": len(contacts_by_component_residue),
                "interface_contact_cutoff_angstrom": contact_cutoff,
                "dssp": dssp_provenance, "sasa": sasa_provenance,
                "terminal_completeness": declared,
                "checks_not_performed": ["tag clearance", "terminal comparison", "SPR/BLI orientation"]}
    if reference_context is None:
        base["checks_not_performed"] = [
            "Structural surroundings outside the selected saved complex and its selected components.",
            "Full-length reference context, additional reference subunits, and reference-only resolved regions.",
            "Glycan-specific accessibility; file annotations are not occupancy or accessibility measurements.",
            "Tag/receptor binding, polymer conformations, dynamics, expression, and experimental SPR/BLI behavior.",
        ]
    else:
        base["checks_not_performed"] = [
            "Unresolved or absent reference atoms, unmodeled glycans, biological assembly corrections, and occupancy validation.",
            "Tag/receptor binding, polymer conformations, dynamics, expression, and experimental SPR/BLI behavior.",
        ]
    base["interface"] = {
        "component_chain": component_chain, "partner_chains": interface_chains,
        "contact_cutoff_angstrom": contact_cutoff,
        "component_residue_count_contacting_partner": len(contacts_by_component_residue),
        "partner_residue_count_contacting_component": len(partner_interface),
        "contact_definition": "at least one inter-chain heavy-atom pair within the saved cutoff",
        "observed_interface": bool(contacts_by_component_residue and partner_interface),
    }
    base["terminal_completeness"] = {
        "declared_sequence": declared,
        "user_confirmed_complete_construct_termini": bool(termini_confirmed),
        "confirmation_rule": "A declared polymer sequence that shows unresolved first/last residues overrides user confirmation.",
    }
    base["dssp"] = dssp_provenance
    base["sasa"] = sasa_provenance
    all_source_other = np.asarray(other_xyz, dtype=float).reshape((-1, 3))
    component_full = np.asarray(component_xyz, dtype=float).reshape((-1, 3))
    interface_full = np.asarray(partner_interface_xyz, dtype=float).reshape((-1, 3))
    interface_distance_xyz = interface_full if len(interface_full) else target_xyz
    obstacles_by_end = {}
    for end in ("N", "C"):
        terminal_residue = component_residues[0 if end == "N" else -1]
        terminal_atoms = {id(atom) for atom in _heavy_atoms(terminal_residue)}
        own_body = [atom.get_coord() for residue in component_residues for atom in _heavy_atoms(residue)
                    if id(atom) not in terminal_atoms]
        obstacles_by_end[end] = np.concatenate([
            np.asarray(own_body, dtype=float).reshape((-1, 3)), interface_full, all_source_other,
            *([reference_obstacle_xyz] if len(reference_obstacle_xyz) else []),
        ], axis=0)
    target_tree = cKDTree(interface_distance_xyz)
    resolved_ca = [atom.get_coord() for atom in ca_atoms if atom is not None]
    binder_centroid = np.mean(resolved_ca, axis=0) if resolved_ca else np.full(3, np.nan)
    n_boundary_ok = terminal_flags["N"] is not False
    c_boundary_ok = terminal_flags["C"] is not False
    terminal_geometry = {}
    for end, valid in (("N", n_boundary_ok), ("C", c_boundary_ok)):
        terminal_geometry[end] = _endpoint_assessment(
            end=end, chain_id=component_chain, residues=component_residues, ca_atoms=ca_atoms,
            target_tree=target_tree, target_xyz=interface_full, obstacle_xyz=obstacles_by_end[end],
            dssp=dssp, sasa=sasa, plddt_enabled=b_factor_interpretation == "plddt_0_100",
            boundary_resolved=bool(valid), termini_confirmed=bool(termini_confirmed),
            interface_residue_keys={residue.id for residue in component_residues if residue.id in contacts_by_component_residue},
            binder_centroid=binder_centroid, interface_center=interface_center,
        )
        if not base["interface"]["observed_interface"]:
            terminal_geometry[end]["assessment_status"] = "not_assessed_no_interface_contacts"
            terminal_geometry[end]["tag_clearance"] = {}
    base["terminal_geometry"] = terminal_geometry
    base["interface_contact_cutoff_angstrom"] = contact_cutoff
    base["confidence_interpretation"] = (
        "Cα B-factors are treated as pLDDT only because the caller explicitly declared that this file stores predictor confidence in B-factor fields."
        if b_factor_interpretation == "plddt_0_100" else
        "Cα B-factors are reported only as uninterpreted coordinate-file metadata; pLDDT is unavailable unless separately supplied with verified provenance."
    )
    eligible = [row for row in terminal_geometry.values() if row.get("tag_clearance")]
    if not base["interface"]["observed_interface"]:
        base["assessment_status"] = "no_interface_contacts_at_selected_cutoff"
        base["tag_comparison"] = {}
    elif not termini_confirmed:
        base["assessment_status"] = "provisional_construct_termini_unconfirmed"
        base["tag_comparison"] = {
            tag["key"]: {"label": tag["label"], "preferred_end_advisory": "not_assessed",
                         "advice_basis": "Confirm the intended construct termini before comparing tag/fusion placement."}
            for tag in TAGS
        }
    elif len(eligible) == 2:
        base["assessment_status"] = "measured"
        comparisons = {}
        for tag in TAGS:
            key = tag["key"]
            n_row, c_row = terminal_geometry["N"], terminal_geometry["C"]
            n_status = n_row.get("tag_clearance", {}).get(key, {}).get("status")
            c_status = c_row.get("tag_clearance", {}).get(key, {}).get("status")
            n_ok = n_status in {"clearance_path_without_added_linker", "clearance_path_with_linker",
                                "clearance_cone_without_added_linker", "clearance_cone_with_linker"}
            c_ok = c_status in {"clearance_path_without_added_linker", "clearance_path_with_linker",
                                "clearance_cone_without_added_linker", "clearance_cone_with_linker"}
            # A measured envelope fit remains in ``tag_clearance``. For folded
            # fusions, however, attaching at a residue that itself participates
            # in the selected interface is not a supported suitability choice.
            if tag["kind"] == "globular_fusion":
                n_ok = n_ok and not bool(n_row.get("contact_at_terminal_residue"))
                c_ok = c_ok and not bool(c_row.get("contact_at_terminal_residue"))
            if not termini_confirmed:
                preferred = "not_assessed"
            elif not n_ok and not c_ok:
                preferred = "neither_end_suitable_under_method"
            elif n_ok and not c_ok:
                preferred = "N"
            elif c_ok and not n_ok:
                preferred = "C"
            else:
                n_side, c_side = n_row.get("interface_side_cosine"), c_row.get("interface_side_cosine")
                if n_side is not None and c_side is not None and ((n_side > 0.25) != (c_side > 0.25)):
                    preferred = "N" if n_side > c_side else "C"
                else:
                    nl = n_row.get("tag_clearance", {}).get(key, {}).get("minimum_linker_residues_tested")
                    cl = c_row.get("tag_clearance", {}).get(key, {}).get("minimum_linker_residues_tested")
                    preferred = "either" if nl == cl else "N" if nl is not None and (cl is None or nl < cl) else "C"
            comparisons[key] = {
                "label": tag["label"], "preferred_end_advisory": preferred,
                "advice_basis": "clearance status, tested linker length, and geometric side relative to the selected interface; no universal quality score",
            }
        his = comparisons.get("his6", {}).get("preferred_end_advisory")
        if his in {"N", "C", "either"}:
            comparisons["spr_bli_immobilization_orientation_advisory"] = {
                "recommendation": his,
                "basis": "His6-compatible end(s) with modeled clearance and geometric orientation away from the selected interface, where resolved.",
                "limitation": "This is an orientation heuristic; it does not model surface chemistry, tag capture, mass transport, or experimental SPR/BLI behavior.",
            }
        else:
            comparisons["spr_bli_immobilization_orientation_advisory"] = {
                "recommendation": "no_supported_end_under_method" if his == "neither_end_suitable_under_method" else "not_assessed",
                "basis": "No experimental assay outcome is predicted.",
            }
        base["tag_comparison"] = comparisons
    else:
        base["assessment_status"] = "partial_terminal_evidence"
        base["tag_comparison"] = {tag["key"]: {"label": tag["label"], "preferred_end_advisory": "not_assessed_incomplete_terminal_structure"} for tag in TAGS}
    return base
