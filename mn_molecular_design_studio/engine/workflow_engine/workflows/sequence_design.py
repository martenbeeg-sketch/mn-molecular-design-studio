from __future__ import annotations

import gzip
import json
import math
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any

from mn_molecular_design_studio.engine.workflow_engine.core.artifacts import Artifact, artifact_path
from mn_molecular_design_studio.engine.workflow_engine.core.candidates import STAGE_GENERATION_BACKBONE, STAGE_SEQUENCE_DESIGN, read_candidates, write_candidates
from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path
from mn_molecular_design_studio.engine.workflow_engine.core.gpu import docker_gpu_args, normalize_gpu_device
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import JobPaths, create_job, finish_job, update_status, write_json
from mn_molecular_design_studio.engine.workflow_engine.core.manifests import load_manifest
from mn_molecular_design_studio.engine.workflow_engine.core.scheduler import apply_docker_cpu_limits_to_steps
from mn_molecular_design_studio.engine.workflow_engine.workflows import chain_roles
from mn_molecular_design_studio.engine.workflow_engine.workflows import refolding as refolding_workflow


DESIGN_GROUP = "design"

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
}


def _rel_path(run_dir: Path, path: Path | None) -> str | None:
    if path is None:
        return None
    try:
        return str(path.relative_to(run_dir))
    except ValueError:
        return str(path)


def _normalize_omit_aas(value: object, default: str = "CX") -> str:
    text = str(default if value is None else value).upper()
    if text.strip() in {"", "EMPTY"}:
        return ""
    letters = "".join(char for char in text if "A" <= char <= "Z")
    return letters or default


def _resolve_candidate_path(source_run_dir: Path, path_text: str | None) -> Path | None:
    if not path_text:
        return None
    return resolve_stored_path(path_text, run_dir=source_run_dir)


def _is_cif_path(path: Path) -> bool:
    name = path.name.lower()
    return path.suffix.lower() in {".cif", ".mmcif"} or name.endswith((".cif.gz", ".mmcif.gz"))


def _chain_ids_from_pdb(path: Path) -> list[str]:
    chains: list[str] = []
    seen: set[str] = set()
    for line in path.read_text(errors="ignore").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        chain = line[21].strip() or "_"
        if chain not in seen:
            chains.append(chain)
            seen.add(chain)
    return chains


def _infer_design_chains(candidate: dict[str, Any], pdb_path: Path, requested: str = "") -> list[str]:
    if requested.strip():
        return [token.strip() for token in re.split(r"[,;\s]+", requested) if token.strip()]
    target_chains = set(candidate.get("target_chains") or [])
    chains = _chain_ids_from_pdb(pdb_path)
    binder_chains = [chain for chain in chains if chain not in target_chains]
    return binder_chains or chains[-1:]


def _strict_target_engine_chains(source_chains: list[str]) -> list[str]:
    return list(chain_roles.TARGET_CHAIN_ORDER[: len(source_chains)])


def _strict_binder_engine_chains(source_chains: list[str], *, reserved: list[str]) -> list[str]:
    reserved_set = set(reserved)
    return [chain for chain in chain_roles.BINDER_CHAIN_ORDER if chain not in reserved_set][: len(source_chains)]


def _source_role_chains(candidate: dict[str, Any], pdb_path: Path, requested_design_chains: str = "") -> tuple[list[str], list[str]]:
    chains = _chain_ids_from_pdb(pdb_path)
    chain_set = set(chains)
    requested_binders = [
        token.strip()[:1]
        for token in re.split(r"[,;\s]+", requested_design_chains.strip())
        if token.strip()
    ]
    binder_chains = [chain for chain in requested_binders if chain in chain_set]
    if not binder_chains:
        binder_chains = [str(chain)[:1] for chain in candidate.get("binder_chains") or [] if str(chain)[:1] in chain_set]
    target_chains = [str(chain)[:1] for chain in candidate.get("target_chains") or [] if str(chain)[:1] in chain_set]
    if binder_chains:
        target_chains = [chain for chain in target_chains if chain not in set(binder_chains)]
    if not binder_chains and target_chains:
        binder_chains = [chain for chain in chains if chain not in set(target_chains)]
    if not target_chains and binder_chains:
        target_chains = [chain for chain in chains if chain not in set(binder_chains)]
    if not binder_chains and chains:
        binder_chains = chains[-1:]
        target_chains = chains[:-1]
    return binder_chains, target_chains


def _ligandmpnn_safe_pdb_lines(lines: list[str]) -> list[str]:
    safe_lines: list[str] = []
    known_residues = set(AA3_TO_1)
    for line in lines:
        if not line.startswith(("ATOM  ", "HETATM")):
            safe_lines.append(line)
            continue
        atom_name = line[12:16].strip() or "CA"
        resname = line[17:20].strip().upper() if len(line) >= 20 else "GLY"
        if line.startswith("ATOM  ") and resname not in known_residues and resname not in {"DA", "DC", "DG", "DT", "DI", "DU", "A", "C", "G", "U", "I"}:
            resname = "GLY"
        chain = (line[21].strip() if len(line) > 21 else "") or "_"
        try:
            resseq = int(line[22:26])
        except ValueError:
            resseq = len(safe_lines) + 1
        coord_values = re.findall(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)", line[27:])
        if len(coord_values) < 3:
            safe_lines.append(line)
            continue
        try:
            x, y, z = (float(coord_values[index]) for index in range(3))
            occupancy = float(coord_values[3]) if len(coord_values) > 3 else 1.0
            bfactor = float(coord_values[4]) if len(coord_values) > 4 else 0.0
        except ValueError:
            safe_lines.append(line)
            continue
        element = (line[76:78].strip() if len(line) >= 78 else "") or re.sub(r"[^A-Za-z]", "", atom_name)[:1] or "C"
        record = "HETATM" if line.startswith("HETATM") else "ATOM  "
        safe_lines.append(
            f"{record}{len(safe_lines) + 1:5d} {atom_name[:4]:>4s} {resname:>3s} {chain[:1]}{resseq:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{bfactor:6.2f}          {element[:2]:>2s}"
        )
    return safe_lines


def _unit_vector(vector: tuple[float, float, float]) -> tuple[float, float, float]:
    length = math.sqrt(sum(value * value for value in vector))
    if length < 1e-6:
        return (1.0, 0.0, 0.0)
    return tuple(value / length for value in vector)  # type: ignore[return-value]


def _cross_vector(
    left: tuple[float, float, float],
    right: tuple[float, float, float],
) -> tuple[float, float, float]:
    return (
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    )


def _add_vector(
    origin: tuple[float, float, float],
    direction: tuple[float, float, float],
    scale: float,
) -> tuple[float, float, float]:
    return tuple(origin[index] + direction[index] * scale for index in range(3))  # type: ignore[return-value]


def _ligandmpnn_complete_backbone_lines(lines: list[str]) -> list[str]:
    residues: dict[tuple[str, int, str], dict[str, Any]] = {}
    order: list[tuple[str, int, str]] = []
    for line in lines:
        if not line.startswith("ATOM  ") or len(line) < 54:
            continue
        atom = line[12:16].strip().upper()
        chain = line[21].strip() or "_"
        try:
            resseq = int(line[22:26])
            coords = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue
        key = (chain, resseq, line[26].strip())
        if key not in residues:
            residues[key] = {"resname": line[17:20].strip().upper() or "GLY", "atoms": {}}
            order.append(key)
        residues[key]["atoms"][atom] = coords
    if not residues:
        return lines
    ca_keys = [key for key in order if "CA" in residues[key]["atoms"]]
    if not ca_keys:
        return lines
    if all({"N", "CA", "C", "O"}.issubset(set(residues[key]["atoms"])) for key in ca_keys):
        return lines

    by_chain: dict[str, list[tuple[str, int, str]]] = {}
    for key in ca_keys:
        by_chain.setdefault(key[0], []).append(key)

    completed: list[str] = []
    serial = 1
    previous_chain = ""
    for key in ca_keys:
        chain, resseq, _icode = key
        if previous_chain and chain != previous_chain:
            completed.append("TER")
        previous_chain = chain
        chain_keys = by_chain[chain]
        position = chain_keys.index(key)
        ca = residues[key]["atoms"]["CA"]
        if position + 1 < len(chain_keys):
            next_ca = residues[chain_keys[position + 1]]["atoms"]["CA"]
            direction = _unit_vector(tuple(next_ca[index] - ca[index] for index in range(3)))  # type: ignore[arg-type]
        elif position > 0:
            prev_ca = residues[chain_keys[position - 1]]["atoms"]["CA"]
            direction = _unit_vector(tuple(ca[index] - prev_ca[index] for index in range(3)))  # type: ignore[arg-type]
        else:
            direction = (1.0, 0.0, 0.0)
        normal = _unit_vector(_cross_vector(direction, (0.0, 0.0, 1.0)))
        if normal == (1.0, 0.0, 0.0):
            normal = _unit_vector(_cross_vector(direction, (0.0, 1.0, 0.0)))
        atom_coords = {
            "N": _add_vector(ca, direction, -1.45),
            "CA": ca,
            "C": _add_vector(ca, direction, 1.52),
            "O": _add_vector(_add_vector(ca, direction, 1.52), normal, 1.23),
        }
        resname = residues[key]["resname"] if residues[key]["resname"] in AA3_TO_1 else "GLY"
        for atom_name, coords in atom_coords.items():
            element = atom_name[0]
            completed.append(
                f"ATOM  {serial:5d} {atom_name:>4s} {resname:>3s} {chain[:1]}{resseq:4d}    "
                f"{coords[0]:8.3f}{coords[1]:8.3f}{coords[2]:8.3f}{1.0:6.2f}{0.0:6.2f}          {element:>2s}"
            )
            serial += 1
    # Non-protein atoms are the conditioning context for LigandMPNN. Keep them
    # when a backbone-only generator required completion of protein atoms.
    protein_keys = set(ca_keys)
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")):
            try:
                key = (line[21].strip() or "_", int(line[22:26]), line[26].strip())
            except ValueError:
                continue
            if key not in protein_keys:
                completed.append(line[:6] + f"{serial:5d}" + line[11:])
                serial += 1
    return completed


def _stage_ligandmpnn_role_normalized_input(
    *,
    source_pdb: Path,
    staged_pdb: Path,
    candidate: dict[str, Any],
    requested_design_chains: str,
) -> dict[str, Any]:
    binder_source_chains, target_source_chains = _source_role_chains(candidate, source_pdb, requested_design_chains)
    target_engine_chains = _strict_target_engine_chains(target_source_chains)
    binder_engine_chains = _strict_binder_engine_chains(binder_source_chains, reserved=target_engine_chains)
    source_to_engine = {
        **dict(zip(target_source_chains, target_engine_chains)),
        **dict(zip(binder_source_chains, binder_engine_chains)),
    }
    lines: list[str] = []
    next_atom = 1
    for source_chain in _chain_ids_from_pdb(source_pdb):
        engine_chain = source_to_engine.get(source_chain)
        if not engine_chain:
            continue
        chain_lines, next_atom = refolding_workflow._renumber_structure_chain(
            source_pdb,
            engine_chain,
            next_atom,
            {source_chain},
        )
        if chain_lines:
            lines.extend(chain_lines + ["TER"])
    if not lines:
        shutil.copy2(source_pdb, staged_pdb)
        staged_lines = staged_pdb.read_text(errors="ignore").splitlines()
    else:
        staged_lines = lines + ["END"]
    staged_lines = _ligandmpnn_safe_pdb_lines(staged_lines)
    staged_lines = _ligandmpnn_complete_backbone_lines(staged_lines)
    staged_lines = [line for line in staged_lines if line.strip() != "END"]
    staged_pdb.write_text("\n".join(staged_lines + ["END", ""]))
    from mn_molecular_design_studio.generation_io import has_ligand_target, restore_ligand_context
    if has_ligand_target(candidate):
        restore_ligand_context(source_pdb,staged_pdb,dict(zip(target_source_chains,target_engine_chains)))
    role_map = chain_roles.build_explicit_role_map(
        binder_source_chains=binder_source_chains,
        binder_engine_chains=binder_engine_chains,
        target_source_chains=target_source_chains,
        target_engine_chains=target_engine_chains,
    )
    refolding_workflow._write_engine_chain_map(
        staged_pdb.parent,
        staged_pdb.stem,
        binder_source_chains=binder_source_chains,
        binder_engine_chains=binder_engine_chains,
        target_source_chains=target_source_chains,
        target_engine_chains=target_engine_chains,
    )
    write_json(
        staged_pdb.with_suffix(".role_normalization.json"),
        {
            "source_structure": str(source_pdb),
            "normalized_structure": str(staged_pdb),
            "chain_role_schema": chain_roles.CHAIN_ROLE_SCHEMA_EXPLICIT_V2,
            "chain_roles": role_map.to_dict(),
            "binder_source_chains": binder_source_chains,
            "binder_chains": binder_engine_chains,
            "target_source_chains": target_source_chains,
            "target_chains": target_engine_chains,
        },
    )
    return {
        "chain_role_schema": chain_roles.CHAIN_ROLE_SCHEMA_EXPLICIT_V2,
        "chain_roles": role_map.to_dict(),
        "binder_source_chains": binder_source_chains,
        "binder_chains": binder_engine_chains,
        "target_source_chains": target_source_chains,
        "target_chains": target_engine_chains,
    }


def _fixed_residues_for_candidate(
    candidate: dict[str, Any], requested: list[str], role_metadata: dict[str, Any]
) -> list[str]:
    raw=candidate.get('raw_metadata') or {}
    if raw.get('fixed_motif_required') and not (raw.get('motif_mapping_verified') and raw.get('fixed_residues')):
        raise ValueError('Sequence redesign requires a verified mapping for every fixed motif residue.')
    # Native generation adapters export fixed output residues. Chain-role
    # normalization must remap them before any sequence model runs.
    authoritative = candidate.get('source_tool') == 'rfdiffusion3_foundry' or raw.get('fixed_motif_required')
    automatic = (raw.get("fixed_residues") or []) if authoritative else []
    combined = list(dict.fromkeys([*automatic, *requested]))
    return _remap_fixed_residues_to_engine_chains(combined, role_metadata)


def _remap_fixed_residues_to_engine_chains(fixed_residues: list[str], role_metadata: dict[str, Any]) -> list[str]:
    source_to_engine: dict[str, str] = {}
    chain_map = (role_metadata.get("chain_roles") or {}).get("chain_map") or {}
    for section in ("binder", "targets"):
        for row in chain_map.get(section) or []:
            if not isinstance(row, dict):
                continue
            original = str(row.get("original_chain") or "").strip()[:1]
            engine = str(row.get("engine_chain") or "").strip()[:1]
            if original and engine:
                source_to_engine[original] = engine
    remapped: list[str] = []
    for residue in fixed_residues:
        text = str(residue or "").strip()
        exact = role_metadata.get("source_to_engine_residues", {}).get(text)
        remapped.append(exact or (f"{source_to_engine[text[:1]]}{text[1:]}" if text[:1] in source_to_engine else text))
    return remapped


def interface_residue_contacts(
    pdb_path: Path,
    binder_chains: list[str],
    target_chains: list[str],
    *,
    cutoff: float = 4.0,
) -> list[str]:
    """Return binder residues with a heavy-atom contact to the target."""
    binder_set = set(binder_chains)
    target_set = set(target_chains)
    binder_atoms: list[tuple[str, tuple[float, float, float]]] = []
    target_atoms: list[tuple[float, float, float]] = []
    for line in pdb_path.read_text(errors="ignore").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")) or len(line) < 54:
            continue
        altloc = line[16].strip()
        if altloc not in {"", "A"}:
            continue
        element = line[76:78].strip().upper() if len(line) >= 78 else ""
        atom_name = line[12:16].strip().upper()
        if element == "H" or (not element and atom_name.startswith("H")):
            continue
        chain = line[21].strip() or "_"
        try:
            coords = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue
        if chain in binder_set:
            residue_number = line[22:26].strip()
            if residue_number:
                binder_atoms.append((f"{chain}{residue_number}{line[26].strip()}", coords))
        elif chain in target_set:
            target_atoms.append(coords)
    if not binder_atoms or not target_atoms:
        return []

    cell_size = max(0.1, float(cutoff))
    target_grid: dict[tuple[int, int, int], list[tuple[float, float, float]]] = {}
    for coords in target_atoms:
        cell = tuple(math.floor(value / cell_size) for value in coords)
        target_grid.setdefault(cell, []).append(coords)
    contacts: set[str] = set()
    cutoff_sq = float(cutoff) ** 2
    for residue, coords in binder_atoms:
        cell = tuple(math.floor(value / cell_size) for value in coords)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for target in target_grid.get((cell[0] + dx, cell[1] + dy, cell[2] + dz), []):
                        if sum((coords[index] - target[index]) ** 2 for index in range(3)) <= cutoff_sq:
                            contacts.add(residue)
                            break
                    if residue in contacts:
                        break
                if residue in contacts:
                    break
            if residue in contacts:
                break
    return sorted(contacts, key=lambda tag: (tag[0], int("".join(ch for ch in tag[1:] if ch.isdigit()) or 0), tag))


def _candidate_sequence_matches_structure(
    candidate: dict[str, Any], structure: Path, binder_chains: list[str]
) -> bool:
    """Only lock an interface when the generator declared a canonical matching sequence."""
    declared = re.sub(r"\s+", "", str(candidate.get("binder_sequence") or "")).upper()
    if not declared or not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", declared):
        return False
    try:
        from mn_molecular_design_studio.engine.workflow_engine.workflows import refolding as refolding_workflow
        by_chain = refolding_workflow._sequences_by_chain(Path(structure))
    except (OSError, ValueError, KeyError, IndexError):
        return False
    observed = "".join(by_chain.get(chain, "") for chain in binder_chains).upper()
    return bool(observed) and re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", observed) is not None and observed == declared


def _interface_contacts_for_structure(
    structure: Path, binder_chains: list[str], target_chains: list[str], *, cutoff: float = 4.0
) -> list[str]:
    """Return contact IDs in the source structure's own chain/residue namespace."""
    structure = Path(structure)
    if not _is_cif_path(structure):
        return interface_residue_contacts(structure, binder_chains, target_chains, cutoff=cutoff)
    from mn_molecular_design_studio.engine.workflow_engine.workflows import mpnn as mpnn_workflow
    converted = structure.with_name(structure.name.removesuffix(".gz").removesuffix(".cif").removesuffix(".mmcif") + ".interface_contacts.pdb")
    chain_map = mpnn_workflow._pdb_copy(structure, converted)
    reverse = {pdb_chain: source_chain for source_chain, pdb_chain in chain_map.items()}
    try:
        contacts = interface_residue_contacts(
            converted,
            [chain_map[chain] for chain in binder_chains if chain in chain_map],
            [chain_map[chain] for chain in target_chains if chain in chain_map],
            cutoff=cutoff,
        )
    finally:
        converted.unlink(missing_ok=True)
    return [reverse.get(label[:1], label[:1]) + label[1:] for label in contacts]


def derive_interface_masks(
    source_run_dir: Path,
    candidates: list[dict[str, Any]],
    *,
    design_chains: str = "",
    cutoff: float = 4.0,
) -> tuple[dict[str, list[str]], list[dict[str, Any]]]:
    """Calculate fixed contact masks and preview rows for PDB complex candidates."""
    masks: dict[str, list[str]] = {}
    preview: list[dict[str, Any]] = []
    for candidate in candidates:
        candidate_id = str(candidate.get("candidate_id") or "")
        pdb_path = _resolve_candidate_path(source_run_dir, candidate.get("complex_pdb"))
        if not candidate_id or pdb_path is None or not pdb_path.exists() or pdb_path.suffix.lower() != ".pdb":
            preview.append({"candidate_id": candidate_id, "status": "Skipped: missing PDB complex"})
            continue
        binder_chains = _infer_design_chains(candidate, pdb_path, design_chains)
        target_chains = [str(chain) for chain in candidate.get("target_chains") or []]
        if not target_chains:
            target_chains = [chain for chain in _chain_ids_from_pdb(pdb_path) if chain not in set(binder_chains)]
        residues = interface_residue_contacts(pdb_path, binder_chains, target_chains, cutoff=cutoff)
        if not residues:
            preview.append({"candidate_id": candidate_id, "status": "Skipped: no contacts found"})
            continue
        masks[candidate_id] = residues
        preview.append({
            "candidate_id": candidate_id,
            "binder_chains": ",".join(binder_chains),
            "target_chains": ",".join(target_chains),
            "fixed_residue_count": len(residues),
            "fixed_residues": ", ".join(residues),
            "status": "Ready",
        })
    return masks, preview


def _read_fasta_sequences(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    if not path.exists():
        return rows
    header = ""
    seq: list[str] = []
    for line in path.read_text(errors="ignore").splitlines():
        if line.startswith(">"):
            if header or seq:
                rows.append({"header": header, "sequence": "".join(seq)})
            header = line[1:].strip()
            seq = []
        else:
            seq.append(line.strip())
    if header or seq:
        rows.append({"header": header, "sequence": "".join(seq)})
    parsed: list[dict[str, Any]] = []
    for row in rows:
        from mn_molecular_design_studio.output_vocabulary import fasta_metrics
        metrics = fasta_metrics(row["header"])
        parsed.append({**row, "metrics": metrics})
    return parsed


def _designed_sequence_from_multichain_fasta(sequence: str, input_pdb: Path | None, design_chains: list[str]) -> str:
    sequence = sequence.strip().replace(" ", "")
    if ":" not in sequence:
        return sequence
    parts = sequence.split(":")
    if not input_pdb or not input_pdb.exists() or not design_chains:
        return parts[0]
    chain_order = _chain_ids_from_pdb(input_pdb)
    by_chain = {chain: parts[index] for index, chain in enumerate(chain_order) if index < len(parts)}
    selected = [by_chain[chain] for chain in design_chains if chain in by_chain]
    return "".join(selected) if selected else parts[0]


def _ligandmpnn_design_index(header: str) -> int | None:
    match = re.search(r"(?:^|,\s*)id=(\d+)", header)
    return int(match.group(1)) if match else None


def _run_shell_steps(run_dir: Path, steps: list[dict]) -> int:
    steps = apply_docker_cpu_limits_to_steps(run_dir, steps)
    write_json(run_dir / "command.json", {"mode": "docker", "steps": steps})
    update_status(run_dir, "running")
    with (run_dir / "stdout.log").open("w") as stdout, (run_dir / "stderr.log").open("w") as stderr:
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
        for path in sorted(root.glob(pattern)):
            if path.is_file():
                artifacts.append(Artifact(path.stem, path, artifact_type).to_json(run_dir))
    return artifacts


def _normalize_ligandmpnn_candidates(
    run_dir: Path,
    source_candidates: list[dict[str, Any]],
    candidate_input_paths: dict[str, Path],
    params: dict[str, Any],
) -> list[dict[str, Any]]:
    normalized: list[dict[str, Any]] = []
    raw_root = run_dir / "artifacts" / "raw" / "ligandmpnn"
    for source in source_candidates:
        source_id = source.get("candidate_id")
        if not source_id:
            continue
        seq_dir = raw_root / "output" / source_id / "seqs"
        fasta_files = sorted(seq_dir.glob("*.fa"))
        input_pdb = candidate_input_paths.get(source_id)
        design_chains = params.get("design_chains_by_candidate", {}).get(source_id, [])
        role_metadata = params.get("role_metadata_by_candidate", {}).get(source_id, {})
        target_chains = list(role_metadata.get("target_chains") or source.get("target_chains", []))
        for fasta_path in fasta_files:
            designed_count = 0
            for seq_row in _read_fasta_sequences(fasta_path):
                header = str(seq_row.get("header") or "")
                if " id=" not in f" {header}":
                    continue
                designed_count += 1
                design_index = _ligandmpnn_design_index(header) or designed_count
                sequence_index = designed_count
                backbone_path = raw_root / "output" / source_id / "backbones" / f"{source_id}_{design_index}.pdb"
                if params.get("pack_side_chains"):
                    packed = sorted((raw_root / "output" / source_id / "packed").glob(f"{source_id}_packed_{design_index}_*.pdb"))
                    if packed:
                        backbone_path = packed[0]
                from mn_molecular_design_studio.generation_io import has_ligand_target, restore_ligand_context
                if backbone_path.exists() and input_pdb and has_ligand_target(source):
                    restore_ligand_context(input_pdb,backbone_path,{chain:chain for chain in target_chains})
                binder_sequence = ""
                if backbone_path.exists() and design_chains:
                    pdb_sequences = refolding_workflow._sequences_by_chain(backbone_path)
                    binder_sequence = "".join(pdb_sequences.get(chain, "") for chain in design_chains)
                if not binder_sequence:
                    binder_sequence = _designed_sequence_from_multichain_fasta(
                        str(seq_row.get("sequence") or ""),
                        input_pdb,
                        design_chains,
                    )
                normalized.append(
                    {
                        "candidate_id": f"{source_id}_{params['model_type']}_{sequence_index:03d}",
                        "stage": STAGE_SEQUENCE_DESIGN,
                        "source_tool": params["model_type"],
                        "target_pdb": _rel_path(run_dir, input_pdb) if params.get("unified_mpnn") else source.get("target_pdb"),
                        "complex_pdb": _rel_path(run_dir, backbone_path if backbone_path.exists() else input_pdb),
                        "binder_pdb": None,
                        "binder_sequence": binder_sequence,
                        "target_chains": target_chains,
                        "binder_chains": design_chains,
                        "chain_role_schema": role_metadata.get("chain_role_schema"),
                        "chain_roles": role_metadata.get("chain_roles"),
                        "hotspots": source.get("hotspots", []),
                        "binder_length": str(len(binder_sequence)) if binder_sequence else source.get("binder_length"),
                        "contig": source.get("contig"),
                        "metrics": seq_row.get("metrics") or {},
                        "parents": [source_id],
                        "raw_metadata": {
                            **dict(source.get('raw_metadata') or {}),
                            "source_candidate": source,
                            "fasta": _rel_path(run_dir, fasta_path),
                            "fasta_header": seq_row.get("header"),
                            "sequence_index": sequence_index,
                            "ligandmpnn_design_index": design_index,
                            "fixed_residues": params.get("fixed_residues_by_candidate", {}).get(
                                source_id, []
                            ),
                            "chain_role_schema": role_metadata.get("chain_role_schema"),
                            "chain_roles": role_metadata.get("chain_roles"),
                            "binder_source_chains": role_metadata.get("binder_source_chains", []),
                            "target_source_chains": role_metadata.get("target_source_chains", []),
                        },
                    }
                )
    return write_candidates(run_dir, params["model_type"], normalized)


def run_ligandmpnn_sequence_design(
    source_run_dir: Path,
    candidates_jsonl: Path,
    model_type: str = "protein_mpnn",
    design_chains: str = "",
    num_seq_per_target: int = 1,
    sampling_temp: float = 0.0001,
    omit_aas: str = "CX",
    seed: int | None = None,
    require_backbone_hotspot_filter_pass: bool = False,
    accepted_stages: list[str] | None = None,
    fixed_residues_by_candidate: dict[str, list[str]] | None = None,
    selected_candidate_ids: list[str] | None = None,
    gpu_device: object = "0",
    existing_job: JobPaths | None = None,
    direct_input: dict[str, Any] | None = None,
) -> Path:
    from mn_molecular_design_studio import mpnn_input
    from mn_molecular_design_studio.engine.workflow_engine.workflows import mpnn
    source_run_dir = Path(source_run_dir)
    candidates_jsonl = Path(candidates_jsonl)
    if model_type not in mpnn_input.MODELS:
        raise ValueError("Choose a supported MPNN model.")
    native_model = mpnn_input.MODELS[model_type]["native"]
    if num_seq_per_target < 1:
        raise ValueError("Sequences per backbone must be at least 1.")
    omit_aas = _normalize_omit_aas(omit_aas)
    accepted_stage_set = set(accepted_stages or [STAGE_GENERATION_BACKBONE])
    source_candidates = [
        candidate for candidate in read_candidates(candidates_jsonl) if candidate.get("stage") in accepted_stage_set
    ]
    if selected_candidate_ids is not None:
        selected_ids = {str(candidate_id) for candidate_id in selected_candidate_ids}
        source_candidates = [candidate for candidate in source_candidates if str(candidate.get("candidate_id") or "") in selected_ids]
    if require_backbone_hotspot_filter_pass:
        source_candidates = [
            candidate
            for candidate in source_candidates
            if (candidate.get("metrics") or {}).get("passes_backbone_hotspot_filter") is not False
        ]
    if not source_candidates:
        raise ValueError(
            "No candidates with an accepted stage were found in the selected candidate set after filters."
        )

    manifest = load_manifest("ligandmpnn")
    fixed_residues_by_candidate = {
        str(candidate_id): [str(residue) for residue in residues if str(residue).strip()]
        for candidate_id, residues in (fixed_residues_by_candidate or {}).items()
    }
    params = {
        "source_run_dir": str(source_run_dir),
        "candidates_jsonl": str(candidates_jsonl),
        "model_type": model_type,
        "design_chains": design_chains,
        "num_seq_per_target": num_seq_per_target,
        "sampling_temp": sampling_temp,
        "omit_AAs": omit_aas,
        "seed": seed,
        "require_backbone_hotspot_filter_pass": require_backbone_hotspot_filter_pass,
        "accepted_stages": sorted(accepted_stage_set),
        "fixed_residues_by_candidate": fixed_residues_by_candidate,
        "selected_candidate_ids": sorted(selected_candidate_ids or []),
        "gpu_device": normalize_gpu_device(gpu_device),
    }
    if direct_input:
        params.update(direct_input)
        params["unified_mpnn"] = True
        params["checkpoint_path"] = str(mpnn_input.checkpoint(model_type, direct_input.get("noise_level")))
    job = existing_job or create_job(
        DESIGN_GROUP,
        job_type="sequence_design",
        tool="ligandmpnn",
        inputs={"source_run_dir": str(source_run_dir), "candidates_jsonl": str(candidates_jsonl)},
        params=params,
    )

    raw_root = job.run_dir / "artifacts" / "raw" / "ligandmpnn"
    input_dir = raw_root / "inputs"
    output_dir = raw_root / "output"
    input_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    candidate_input_paths: dict[str, Path] = {}
    design_chains_by_candidate: dict[str, list[str]] = {}
    role_metadata_by_candidate: dict[str, dict[str, Any]] = {}
    steps: list[dict] = []

    for candidate in source_candidates:
        source_id = str(candidate["candidate_id"])
        candidate_input = dict(direct_input or {})
        source_pdb = _resolve_candidate_path(source_run_dir, candidate.get("complex_pdb") or candidate.get("binder_pdb"))
        if source_pdb is None or not source_pdb.exists():
            continue
        if _is_cif_path(source_pdb):
            converted_pdb = artifact_path(job.run_dir, "raw", "ligandmpnn", "inputs", f"{source_id}.source.pdb")
            if direct_input:
                original_info = mpnn_input.inventory(source_pdb)
                mapping = mpnn._pdb_copy(source_pdb, converted_pdb)
                candidate = dict(candidate)
                for role in ("binder_chains", "target_chains"):
                    candidate[role] = [mapping.get(c, c) for c in candidate.get(role, [])]
                labels = {r["label"]: f"{mapping[c]}{r['number']}{r['insertion_code']}" for c, rows in original_info["chains"].items() for r in rows["residues"]}
                candidate_input["design_chains"] = ",".join(mapping.get(c, c) for c in mpnn_input.chain_list(candidate_input.get("design_chains")))
                for field in ("designed_residues", "fixed_residues", "transmembrane_buried", "transmembrane_interface"):
                    candidate_input[field] = " ".join(labels[r] for r in mpnn_input.residues(candidate_input.get(field), original_info))
                for field in ("bias_aa_per_residue", "omit_aa_per_residue"):
                    candidate_input[field] = {labels[r]: value for r, value in candidate_input.get(field, {}).items()}
                raw = dict(candidate.get("raw_metadata") or {})
                if raw.get("fixed_residues"):
                    raw["fixed_residues"] = [labels.get(r, r) for r in raw["fixed_residues"]]
                candidate["raw_metadata"] = raw
            else:
                refolding_workflow._cif_to_pdb(source_pdb, converted_pdb)
            source_pdb = converted_pdb
        if source_pdb.suffix.lower() != ".pdb":
            continue
        if direct_input and candidate_input.get("design_chains"):
            selected = mpnn_input.chain_list(candidate_input["design_chains"])
            available = {c for c, rows in mpnn_input.inventory(source_pdb)["chains"].items() if rows["kind"] == "protein"}
            if not set(selected) <= available:
                raise ValueError("Every requested designed chain must be a protein chain in this input design.")
        staged_pdb = artifact_path(job.run_dir, "raw", "ligandmpnn", "inputs", f"{source_id}.pdb")
        role_metadata = _stage_ligandmpnn_role_normalized_input(
            source_pdb=source_pdb,
            staged_pdb=staged_pdb,
            candidate=candidate,
            requested_design_chains=candidate_input.get("design_chains", design_chains),
        )
        extended = None
        if direct_input:
            role_metadata["source_to_engine_residues"] = mpnn.residue_map(source_pdb, role_metadata)
            extended = mpnn.candidate_options(
                candidate_input,
                source_pdb,
                role_metadata["binder_source_chains"],
                role_metadata,
                normalized_source=staged_pdb,
            )
        auto_interface_fixed: list[str] = []
        if direct_input and direct_input.get("fix_interface_residues") and _candidate_sequence_matches_structure(
            candidate, staged_pdb, role_metadata.get("binder_chains") or []
        ):
            auto_interface_fixed = _interface_contacts_for_structure(
                staged_pdb,
                role_metadata.get("binder_chains") or [],
                role_metadata.get("target_chains") or [],
                cutoff=4.0,
            )
        candidate_input_paths[source_id] = staged_pdb
        chains = list(role_metadata.get("binder_chains") or _infer_design_chains(candidate, staged_pdb, design_chains))
        design_chains_by_candidate[source_id] = chains
        role_metadata_by_candidate[source_id] = role_metadata
        (output_dir / source_id).mkdir(parents=True, exist_ok=True)
        args = [
            "python",
            "/opt/LigandMPNN/run.py",
            "--pdb_path",
            f"/work/artifacts/raw/ligandmpnn/inputs/{source_id}.pdb",
            "--out_folder",
            f"/work/artifacts/raw/ligandmpnn/output/{source_id}",
            "--model_type",
            native_model,
            "--chains_to_design",
            ",".join(chains),
            "--number_of_batches",
            str(num_seq_per_target),
            "--batch_size",
            "1",
            "--temperature",
            str(sampling_temp),
            "--omit_AA",
            omit_aas,
            "--verbose",
            "1",
        ]
        if seed is not None:
            args.extend(["--seed", str(seed)])
        fixed_residues = _fixed_residues_for_candidate(
            candidate, fixed_residues_by_candidate.get(source_id, []), role_metadata
        )
        fixed_residues = list(dict.fromkeys([*fixed_residues, *auto_interface_fixed]))
        if auto_interface_fixed:
            fixed_residues_by_candidate[source_id] = list(dict.fromkeys([
                *fixed_residues_by_candidate.get(source_id, []), *auto_interface_fixed
            ]))
        if extended:
            fixed_residues = list(dict.fromkeys([*fixed_residues, *extended["fixed_residues"]]))
            args += mpnn.shared_arguments(extended, staged_pdb.parent / source_id, "/work/artifacts/raw/ligandmpnn/inputs/" + source_id)
            args += ["--checkpoint_" + native_model, "/mpnn-checkpoint.pt"]
        if fixed_residues:
            fixed_residues_by_candidate[source_id] = fixed_residues
        if fixed_residues:
            args.extend(["--fixed_residues", " ".join(fixed_residues)])
        command = [
            "docker",
            "run",
            "--rm",
            *docker_gpu_args(gpu_device),
            "-v",
            f"{job.run_dir}:/work",
            *(["-v", f"{mpnn_input.checkpoint(model_type, direct_input.get('noise_level')).resolve()}:/mpnn-checkpoint.pt:ro"] if direct_input else []),
            "-w",
            "/opt/LigandMPNN",
            manifest["image"],
            *args,
        ]
        steps.append({"name": f"{model_type}-{source_id}", "command": command})

    if not steps:
        raise ValueError("No selected candidates had an existing PDB structure file. Foundry CIF outputs should use the Foundry-native MPNN path or be converted first.")

    params["design_chains_by_candidate"] = design_chains_by_candidate
    params["role_metadata_by_candidate"] = role_metadata_by_candidate
    write_json(raw_root / "ligandmpnn_params.json", params)
    rc = _run_shell_steps(job.run_dir, steps)
    candidates = _normalize_ligandmpnn_candidates(job.run_dir, source_candidates, candidate_input_paths, params) if rc == 0 else []
    artifacts = _collect_artifacts(
        job.run_dir,
        job.run_dir / "artifacts",
        [
            ("raw/ligandmpnn/**/*.fa", "fasta"),
            ("raw/ligandmpnn/**/*.pdb", "pdb"),
            ("raw/ligandmpnn/**/*.json", "sequence_design_input"),
            ("normalized_candidates/*.jsonl", "normalized_candidates"),
            ("normalized_candidates/*.json", "campaign_result"),
        ],
    )
    finish_job(
        job.run_dir,
        rc == 0 and bool(candidates),
        {
            "outputs": {"artifacts": artifacts, "candidates": candidates},
            "metrics": {"return_code": rc, "artifact_count": len(artifacts), "candidate_count": len(candidates)},
            "downstream_artifacts": {
                "candidates_jsonl": "artifacts/normalized_candidates/candidates.jsonl",
                "campaign_result": "artifacts/normalized_candidates/campaign_result.json",
            },
        },
    )
    return job.run_dir


def _stage_foundry_structure(source_path: Path, staged_path: Path) -> None:
    staged_path.parent.mkdir(parents=True, exist_ok=True)
    if source_path.name.endswith(".gz"):
        with gzip.open(source_path, "rb") as src, staged_path.open("wb") as dst:
            shutil.copyfileobj(src, dst)
    else:
        shutil.copy2(source_path, staged_path)


def _cif_atom_rows(path: Path) -> list[dict[str, str]]:
    text = gzip.open(path, "rt", errors="ignore").read() if path.name.endswith(".gz") else path.read_text(errors="ignore")
    rows: list[dict[str, str]] = []
    atom_headers: list[str] = []
    in_atom_loop = False
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if line == "loop_":
            atom_headers = []
            in_atom_loop = False
            continue
        if line.startswith("_atom_site."):
            atom_headers.append(line.split(".", 1)[1])
            in_atom_loop = True
            continue
        if not in_atom_loop or not line.startswith(("ATOM ", "HETATM ")):
            continue
        parts = line.split()
        if len(parts) >= len(atom_headers):
            rows.append(dict(zip(atom_headers, parts)))
    return rows


def _chains_from_cif(path: Path) -> list[str]:
    chains: list[str] = []
    seen: set[str] = set()
    for row in _cif_atom_rows(path):
        chain = row.get("auth_asym_id") or row.get("label_asym_id") or "_"
        if chain not in seen:
            seen.add(chain)
            chains.append(chain)
    return chains


def _cif_auth_to_label_chain_map(path: Path) -> dict[str, str]:
    """Map input author chain IDs to the label IDs Foundry writes to outputs."""
    mapping: dict[str, str] = {}
    for row in _cif_atom_rows(path):
        author = row.get("auth_asym_id") or ""
        label = row.get("label_asym_id") or ""
        if author not in {"", ".", "?"} and label not in {"", ".", "?"}:
            mapping.setdefault(author, label)
    return mapping


def _sequence_from_cif(path: Path, chains: list[str] | None = None) -> str:
    keep_chains = set(chains or [])
    residues: list[tuple[str, int, str, str]] = []
    seen: set[tuple[str, str, str]] = set()
    for row in _cif_atom_rows(path):
        if (row.get("label_atom_id") or row.get("auth_atom_id")) != "CA":
            continue
        chain = row.get("auth_asym_id") or row.get("label_asym_id") or "_"
        if keep_chains and chain not in keep_chains:
            continue
        residue = row.get("auth_seq_id") or row.get("label_seq_id") or "0"
        insertion = row.get("pdbx_PDB_ins_code") or ""
        key = (chain, residue, insertion)
        if key in seen:
            continue
        seen.add(key)
        try:
            residue_number = int(residue)
        except ValueError:
            residue_number = 0
        residues.append((chain, residue_number, insertion, AA3_TO_1.get((row.get("label_comp_id") or row.get("auth_comp_id") or "").upper(), "X")))
    residues.sort(key=lambda item: (item[0], item[1], item[2]))
    return "".join(row[3] for row in residues)


def _normalize_foundry_mpnn_candidates(
    run_dir: Path,
    source_candidates: list[dict[str, Any]],
    candidate_input_paths: dict[str, Path],
    params: dict[str, Any],
) -> list[dict[str, Any]]:
    normalized: list[dict[str, Any]] = []
    raw_root = run_dir / "artifacts" / "raw" / "foundry_mpnn"
    for source in source_candidates:
        source_id = source.get("candidate_id")
        if not source_id:
            continue
        output_dir = raw_root / "output" / source_id
        cif_files = sorted([*output_dir.glob("*.cif"), *output_dir.glob("*.cif.gz")])
        input_path = candidate_input_paths.get(source_id)
        input_chain_map = _cif_auth_to_label_chain_map(input_path) if input_path and input_path.is_file() else {}

        def output_chain(chain: str) -> str:
            return input_chain_map.get(str(chain), str(chain))

        for index, cif_path in enumerate(cif_files, start=1):
            output_chains = _chains_from_cif(cif_path)
            target_chains = [
                mapped for chain in source.get("target_chains", [])
                if (mapped := output_chain(str(chain))) in output_chains
            ]
            binder_chains = [
                mapped for chain in source.get("binder_chains", [])
                if (mapped := output_chain(str(chain))) in output_chains
            ]
            if not binder_chains:
                # Older results may not include explicit binder roles. Infer
                # only protein chains here; non-protein context such as DNA
                # must never be marked as a designed binder.
                binder_chains = [
                    chain for chain in output_chains
                    if chain not in set(target_chains) and _sequence_from_cif(cif_path, [chain])
                ]
            if not target_chains and len(output_chains) == 1 and source.get("source_tool") != "rfdiffusion3_foundry":
                target_chains = []
                binder_chains = binder_chains or output_chains
            if not binder_chains:
                continue
            role_map = chain_roles.build_explicit_role_map(
                binder_source_chains=binder_chains,
                binder_engine_chains=binder_chains,
                target_source_chains=target_chains or source.get("target_chains", []),
                target_engine_chains=target_chains or source.get("target_chains", []),
            )
            binder_sequence = _sequence_from_cif(cif_path, binder_chains or None)
            normalized.append(
                {
                    "candidate_id": f"{source_id}_foundry_mpnn_{index:03d}",
                    "stage": STAGE_SEQUENCE_DESIGN,
                    "source_tool": "foundry_mpnn",
                    "target_pdb": _rel_path(run_dir, input_path) if params.get("unified_mpnn") else source.get("target_pdb"),
                    "complex_pdb": _rel_path(run_dir, cif_path),
                    "binder_pdb": None,
                    "binder_sequence": binder_sequence,
                    "target_chains": target_chains or source.get("target_chains", []),
                    "binder_chains": binder_chains,
                    "chain_role_schema": chain_roles.CHAIN_ROLE_SCHEMA_EXPLICIT_V2,
                    "chain_roles": role_map.to_dict(),
                    "hotspots": source.get("hotspots", []),
                    "binder_length": source.get("binder_length"),
                    "contig": source.get("contig"),
                    "metrics": {},
                    "parents": [source_id],
                    "raw_metadata": {
                        "source_candidate": source,
                        "input_cif": _rel_path(run_dir, input_path),
                        "mpnn_cif": _rel_path(run_dir, cif_path),
                        "output_chains": output_chains,
                        "chain_role_schema": chain_roles.CHAIN_ROLE_SCHEMA_EXPLICIT_V2,
                        "chain_roles": role_map.to_dict(),
                        "binder_source_chains": binder_chains,
                        "target_source_chains": target_chains or source.get("target_chains", []),
                    },
                }
            )
    return write_candidates(run_dir, "foundry_mpnn", normalized)


def run_foundry_mpnn_sequence_design(
    source_run_dir: Path,
    candidates_jsonl: Path,
    number_of_batches: int = 1,
    batch_size: int = 10,
    model_type: str = "ligand_mpnn",
    checkpoint_path: str = "/weights/ligandmpnn_v_32_010_25.pt",
    selected_candidate_ids: list[str] | None = None,
    gpu_device: object = "0",
    existing_job: JobPaths | None = None,
    direct_input: dict[str, Any] | None = None,
) -> Path:
    from mn_molecular_design_studio import mpnn_input
    from mn_molecular_design_studio.engine.workflow_engine.workflows import mpnn
    source_run_dir = Path(source_run_dir)
    candidates_jsonl = Path(candidates_jsonl)
    if model_type not in {"protein_mpnn", "ligand_mpnn"}:
        raise ValueError("Foundry MPNN supports protein_mpnn or ligand_mpnn.")
    if number_of_batches < 1:
        raise ValueError("Number of batches must be at least 1.")
    if batch_size < 1:
        raise ValueError("Batch size must be at least 1.")
    source_candidates = [
        candidate for candidate in read_candidates(candidates_jsonl) if candidate.get("stage") in (mpnn.STAGES if direct_input else [STAGE_GENERATION_BACKBONE])
    ]
    if selected_candidate_ids is not None:
        selected_ids = {str(candidate_id) for candidate_id in selected_candidate_ids}
        source_candidates = [candidate for candidate in source_candidates if str(candidate.get("candidate_id") or "") in selected_ids]
    if not source_candidates:
        raise ValueError("No generation.backbone candidates were found in the selected candidate set.")

    manifest = load_manifest("foundry_mpnn")
    params = {
        "source_run_dir": str(source_run_dir),
        "candidates_jsonl": str(candidates_jsonl),
        "number_of_batches": number_of_batches,
        "batch_size": batch_size,
        "model_type": model_type,
        "checkpoint_path": checkpoint_path,
        "selected_candidate_ids": sorted(selected_candidate_ids or []),
        "gpu_device": normalize_gpu_device(gpu_device),
    }
    if direct_input:
        params.update(direct_input)
        params["unified_mpnn"] = True
        params["checkpoint_path"] = str(mpnn_input.checkpoint(model_type, direct_input.get("noise_level")))
    job = existing_job or create_job(
        DESIGN_GROUP,
        job_type="sequence_design",
        tool="foundry_mpnn",
        inputs={"source_run_dir": str(source_run_dir), "candidates_jsonl": str(candidates_jsonl)},
        params=params,
    )

    raw_root = job.run_dir / "artifacts" / "raw" / "foundry_mpnn"
    input_dir = raw_root / "inputs"
    output_dir = raw_root / "output"
    input_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    candidate_input_paths: dict[str, Path] = {}
    steps: list[dict] = []

    for candidate in source_candidates:
        if candidate.get("source_tool") == "rfdiffusion3_foundry" and not candidate.get("binder_chains"):
            continue
        source_id = str(candidate["candidate_id"])
        source_structure = _resolve_candidate_path(source_run_dir, candidate.get("complex_pdb") or candidate.get("binder_pdb"))
        if source_structure is None or not source_structure.exists():
            continue
        if not (source_structure.name.endswith(".cif") or source_structure.name.endswith(".cif.gz")):
            continue
        staged_cif = artifact_path(job.run_dir, "raw", "foundry_mpnn", "inputs", f"{source_id}.cif")
        _stage_foundry_structure(source_structure, staged_cif)
        all_chains = _chains_from_cif(staged_cif)
        target_chains = set(candidate.get("target_chains") or [])
        binder_chains = [chain for chain in (candidate.get("binder_chains") or all_chains) if chain in all_chains and chain not in target_chains]
        if not binder_chains:
            raise ValueError(f"No binder chains were identified for {source_id}; refusing to redesign the target.")
        candidate_input_paths[source_id] = staged_cif
        (output_dir / source_id).mkdir(parents=True, exist_ok=True)
        extra_args = []
        extended = None
        if direct_input:
            requested = mpnn_input.chain_list(direct_input.get("design_chains"))
            binder_chains = requested or binder_chains
            extended = mpnn.candidate_options(direct_input, staged_cif, binder_chains)
            candidate["binder_chains"] = binder_chains
            candidate["target_chains"] = [chain for chain in all_chains if chain not in binder_chains]
            fixed = set(_fixed_residues_for_candidate(candidate, [], {})) | set(extended["fixed_residues"])
            auto_interface_fixed: list[str] = []
            if direct_input.get("fix_interface_residues") and _candidate_sequence_matches_structure(
                candidate, staged_cif, binder_chains
            ):
                auto_interface_fixed = _interface_contacts_for_structure(
                    staged_cif, binder_chains, candidate["target_chains"], cutoff=4.0
                )
                fixed.update(auto_interface_fixed)
                if auto_interface_fixed:
                    raw = dict(candidate.get("raw_metadata") or {})
                    raw["auto_interface_fixed_residues"] = auto_interface_fixed
                    candidate["raw_metadata"] = raw
            if fixed:
                info = mpnn_input.inventory(staged_cif)
                designed = [r["label"] for c in binder_chains for r in info["chains"][c]["residues"] if r["label"] not in fixed]
                if not designed:
                    raise ValueError("No residues remain to design after applying fixed residues.")
                selection_args = ["--designed_residues", json.dumps(designed)]
            else:
                selection_args = ["--designed_chains", ",".join(binder_chains)]
            aa3 = {v: k for k, v in AA3_TO_1.items()} | {"X": "UNK"}
            extra_args += ["--temperature", str(extended["sampling_temp"]), "--omit", json.dumps([aa3[a] for a in extended["omit_aas"]])]
            if extended["seed"] is not None:
                extra_args += ["--seed", str(extended["seed"])]
            if extended["bias_aa"]:
                extra_args += ["--bias", json.dumps({aa3[a]: float(score) for a, score in (v.split(":") for v in extended["bias_aa"].split(","))})]
            if extended["bias_aa_per_residue"]:
                extra_args += ["--bias_per_residue", json.dumps({r: {aa3[a]: score for a, score in value.items()} for r, value in extended["bias_aa_per_residue"].items()})]
            if extended["omit_aa_per_residue"]:
                extra_args += ["--omit_per_residue", json.dumps({r: [aa3[a] for a in value] for r, value in extended["omit_aa_per_residue"].items()})]
            if extended["symmetry_groups"]:
                extra_args += ["--symmetry_residues", json.dumps(extended["symmetry_groups"])]
            if extended["use_side_chain_context"]:
                extra_args += ["--atomize_side_chains", "True"]
        else:
            selection_args = ["--designed_chains", ",".join(binder_chains)]
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
            "/mnt/db/reference_files/protpardelle-1c/model_params/LigandMPNN:/ligandmpnn_weights:ro",
            *(["-v", f"{mpnn_input.checkpoint(model_type, direct_input.get('noise_level')).resolve()}:/mpnn-checkpoint.pt:ro"] if direct_input else []),
            "-w",
            "/work",
            manifest["image"],
            "mpnn",
            "--structure_path",
            f"/work/artifacts/raw/foundry_mpnn/inputs/{source_id}.cif",
            "--checkpoint_path",
            "/mpnn-checkpoint.pt" if direct_input else checkpoint_path,
            "--is_legacy_weights",
            "True",
            "--model_type",
            model_type,
            *selection_args,
            "--batch_size",
            str(batch_size),
            "--number_of_batches",
            str(number_of_batches),
            "--remove_waters",
            "True",
            "--out_directory",
            f"/work/artifacts/raw/foundry_mpnn/output/{source_id}",
            *extra_args,
        ]
        steps.append({"name": f"foundry-mpnn-{source_id}", "command": command})

    if not steps:
        raise ValueError("No selected candidates had an existing CIF structure file. RFdiffusion PDB outputs should use the shared LigandMPNN backend.")

    write_json(raw_root / "foundry_mpnn_params.json", params)
    rc = _run_shell_steps(job.run_dir, steps)
    candidates = _normalize_foundry_mpnn_candidates(job.run_dir, source_candidates, candidate_input_paths, params) if rc == 0 else []
    artifacts = _collect_artifacts(
        job.run_dir,
        job.run_dir / "artifacts",
        [
            ("raw/foundry_mpnn/**/*.cif", "cif"),
            ("raw/foundry_mpnn/**/*.json", "sequence_design_input"),
            ("normalized_candidates/*.jsonl", "normalized_candidates"),
            ("normalized_candidates/*.json", "campaign_result"),
        ],
    )
    finish_job(
        job.run_dir,
        rc == 0 and bool(candidates),
        {
            "outputs": {"artifacts": artifacts, "candidates": candidates},
            "metrics": {"return_code": rc, "artifact_count": len(artifacts), "candidate_count": len(candidates)},
            "downstream_artifacts": {
                "candidates_jsonl": "artifacts/normalized_candidates/candidates.jsonl",
                "campaign_result": "artifacts/normalized_candidates/campaign_result.json",
            },
        },
    )
    return job.run_dir


def enqueue_sequence_design_pipeline(
    *,
    backend: str,
    source_run_dir: Path,
    candidates_jsonl: Path,
    design_kwargs: dict[str, Any],
    validation_kwargs: dict[str, Any] | None = None,
) -> Path:
    """Create a durable sequence-design job with an optional validation continuation.

    The local worker owns both phases. Validation is enqueued only after the
    normalized sequence candidates exist, so browser navigation cannot interrupt
    the handoff.
    """
    tool = "ligandmpnn" if backend == "shared_ligandmpnn" else "foundry_mpnn"
    job = create_job(
        DESIGN_GROUP,
        job_type="sequence_design",
        tool=tool,
        inputs={"source_run_dir": str(source_run_dir), "candidates_jsonl": str(candidates_jsonl)},
        params={
            **design_kwargs,
            "backend": backend,
            "validation_requested": validation_kwargs is not None,
            "gpu_device": normalize_gpu_device(design_kwargs.get("gpu_device", "0")),
        },
    )
    write_json(
        job.run_dir / "worker_request.json",
        {
            "kind": "sequence_design_pipeline",
            "kwargs": {
                "backend": backend,
                "source_run_dir": str(source_run_dir),
                "candidates_jsonl": str(candidates_jsonl),
                "design_kwargs": design_kwargs,
                "validation_kwargs": validation_kwargs or {},
            },
        },
    )
    write_json(
        job.run_dir / "command.json",
        {
            "mode": "local_worker",
            "command": ["python", "-m", "mn_molecular_design_studio.engine.workflow_engine.core.local_worker", "--run-dir", str(job.run_dir)],
        },
    )
    return job.run_dir
