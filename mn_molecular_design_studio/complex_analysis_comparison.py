"""Compatibility helpers for browsing saved complex structures side by side.

Validates component selections and exposes saved evidence. Explicit target
superposition reuses Studio's existing Kabsch fitter and records its inputs.
"""
from __future__ import annotations

import gzip
import hashlib
import io
import json
from pathlib import Path
from typing import Any


EXTERNAL_ALIGNMENT_SCHEMA = "mnds.external-alignment/v1"


def external_alignment_template(reference: dict, moving: dict) -> dict[str, Any]:
    """Create a versioned import template bound to the current selections.

    The identity matrix is a placeholder only. The user must replace it with
    an externally produced moving-to-reference transform and provide explicit
    residue correspondence before Studio will display the overlay.
    """
    def selection(item: dict) -> dict[str, str]:
        return {
            "identity": str(item.get("identity") or ""),
            "model_id": str(item.get("model_id") or ""),
            "structure_digest": str(item.get("structure_digest") or ""),
        }

    return {
        "schema": EXTERNAL_ALIGNMENT_SCHEMA,
        "reference": selection(reference),
        "moving": selection(moving),
        "transform": {
            "direction": "moving_to_reference",
            "matrix_4x4": [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]],
            "method": "",
            "source": {"tool": "", "version": "", "record_id": ""},
        },
        "residue_correspondence": [],
        "mapping_note": "",
    }


def normalize_external_alignment_record(payload: Any, reference: dict, moving: dict) -> dict[str, Any]:
    """Validate the versioned exchange envelope and adapt it for validation.

    This checks record structure and exact selected identities only. Detailed
    matrix and residue existence validation remains in
    :func:`validate_supplied_alignment`.
    """
    if not isinstance(payload, dict):
        raise ValueError("The external alignment record must be a JSON object.")
    if payload.get("schema") != EXTERNAL_ALIGNMENT_SCHEMA:
        raise ValueError(f"Unsupported or missing schema; expected {EXTERNAL_ALIGNMENT_SCHEMA}.")
    transform = payload.get("transform")
    if not isinstance(transform, dict):
        raise ValueError("The external alignment record needs a transform object.")
    if transform.get("direction") != "moving_to_reference":
        raise ValueError("Transform direction must be explicitly moving_to_reference.")
    source = transform.get("source")
    if not isinstance(source, dict):
        raise ValueError("Transform source must include tool, version, and optional record_id.")
    tool = str(source.get("tool") or "").strip()
    version = str(source.get("version") or "").strip()
    method = str(transform.get("method") or "").strip()
    if not tool or not version:
        raise ValueError("External source tool and version are required.")
    if not method:
        raise ValueError("External transform method is required.")

    def check_selection(label: str, supplied: Any, selected: dict) -> None:
        if not isinstance(supplied, dict):
            raise ValueError(f"The {label} selection identity is missing.")
        for field, selection_field in (("identity", "identity"), ("model_id", "model_id"),
                                       ("structure_digest", "structure_digest")):
            if str(supplied.get(field) or "") != str(selected.get(selection_field) or ""):
                raise ValueError(f"The record's {label} {field} does not match the current selection.")

    check_selection("reference", payload.get("reference"), reference)
    check_selection("moving", payload.get("moving"), moving)
    correspondence = payload.get("residue_correspondence")
    if not isinstance(correspondence, list):
        raise ValueError("residue_correspondence must be a JSON array.")
    return {
        "reference_identity": str(reference.get("identity") or ""),
        "reference_model_id": str(reference.get("model_id") or ""),
        "reference_digest": str(reference.get("structure_digest") or ""),
        "moving_identity": str(moving.get("identity") or ""),
        "moving_model_id": str(moving.get("model_id") or ""),
        "moving_digest": str(moving.get("structure_digest") or ""),
        "matrix_4x4": transform.get("matrix_4x4"),
        "residue_correspondence": correspondence,
        "source": f"{tool} {version}",
        "method": method,
        "mapping_note": str(payload.get("mapping_note") or ""),
        "external_provenance": {
            "schema": EXTERNAL_ALIGNMENT_SCHEMA,
            # Retain the source object verbatim so optional external fields
            # such as accession, URI, and producer parameters are not lost.
            "source": dict(source),
        },
    }


def _read_biopython_structure(path: Path | str):
    from Bio.PDB import MMCIFParser, PDBParser

    source = Path(path)
    suffix = source.name.casefold()
    is_cif = suffix.endswith((".cif", ".mmcif", ".cif.gz", ".mmcif.gz"))
    if is_cif:
        class _OptionalAtomSiteFieldsMMCIFParser(MMCIFParser):
            def _build_structure(self, structure_id):
                atom_count = len(self._mmcif_dict.get("_atom_site.id", []))
                self._mmcif_dict.setdefault("_atom_site.occupancy", ["1.0"] * atom_count)
                self._mmcif_dict.setdefault("_atom_site.B_iso_or_equiv", ["nan"] * atom_count)
                return super()._build_structure(structure_id)

        parser = _OptionalAtomSiteFieldsMMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    if suffix.endswith(".gz"):
        return parser.get_structure(
            "comparison",
            io.StringIO(gzip.decompress(source.read_bytes()).decode("utf-8", errors="replace")),
        )
    return parser.get_structure("comparison", str(source))


def protein_chain_numbering(path: Path | str, chains: list[str]) -> list[dict[str, Any]]:
    """Return selected protein chain sequences and author residue numbering.

    Chain IDs are omitted from the returned rows so equivalent chains can be
    mapped across structures with different chain labels. Residue IDs retain
    insertion codes, and rows are sorted by sequence and numbering to make a
    multichain group deterministic.
    """
    if not chains:
        return []
    try:
        try:
            from Bio.Data.PDBData import protein_letters_3to1_extended
        except ImportError:
            from Bio.Data.IUPACData import protein_letters_3to1_extended as _letters
            protein_letters_3to1_extended = {key.upper(): value for key, value in _letters.items()}
            protein_letters_3to1_extended.setdefault("MSE", "M")
    except ImportError as exc:
        raise ValueError("Biopython is required to verify residue numbering.") from exc

    try:
        # MMCIFParser defaults to author chain and residue IDs, which are the
        # identifiers exposed in Studio's saved contact maps.
        structure = _read_biopython_structure(path)
        model = next(structure.get_models())
    except (OSError, StopIteration, ValueError) as exc:
        raise ValueError(f"Could not inspect residue numbering: {exc}") from exc

    result = []
    for chain_id in chains:
        chain = model[chain_id] if chain_id in model else None
        if chain is None:
            continue
        sequence = []
        residue_ids = []
        for residue in chain:
            code = str(residue.resname).strip().upper()
            letter = protein_letters_3to1_extended.get(code)
            if not letter:
                continue
            sequence.append(letter)
            residue_ids.append({"number": int(residue.id[1]), "insertion": str(residue.id[2]).strip()})
        if sequence:
            result.append({"sequence": "".join(sequence), "numbering": residue_ids})
    return sorted(result, key=lambda row: (row["sequence"], tuple((r["number"], r["insertion"]) for r in row["numbering"])))


def target_structure_fingerprint(path: Path | str, chains: list[str]) -> str:
    """Fingerprint the exact selected target atoms and author numbering.

    Interface Clustering compares contacted target residues directly, so a
    matching sequence and numbering alone are not enough: every item must use
    the same target coordinates and author chain IDs. This fingerprint is an
    input-compatibility check; it does not calculate contacts or a score.
    Parsed coordinate values are retained at their stored precision so a
    target with even a small coordinate change does not silently pass.
    """
    from Bio.PDB.Polypeptide import is_aa

    selected_chains = [str(chain).strip() for chain in chains if str(chain).strip()]
    if not selected_chains or len(selected_chains) != len(set(selected_chains)):
        raise ValueError("Choose one or more unique target protein chains.")
    try:
        structure = _read_biopython_structure(path)
        models = list(structure.get_models())
    except (OSError, ValueError, ImportError) as exc:
        raise ValueError(f"Could not inspect the selected target structure: {exc}") from exc
    if len(models) != 1:
        raise ValueError(f"Interface Clustering requires one structural model; found {len(models)}.")
    model = models[0]
    rows = []
    for chain_id in sorted(selected_chains):
        if chain_id not in model:
            raise ValueError(f"Selected target chain {chain_id!r} is missing from the chosen model.")
        for residue in model[chain_id]:
            if not is_aa(residue, standard=False):
                continue
            residue_name = str(residue.resname).strip().upper()
            for atom in residue.get_unpacked_list():
                element = str(atom.element or "").strip().upper()
                if element in {"H", "D"}:
                    continue
                xyz = tuple(float(value) for value in atom.coord)
                rows.append((
                    chain_id, int(residue.id[1]), str(residue.id[2]).strip(), residue_name,
                    str(atom.name).strip(), element, str(atom.get_altloc() or "").strip(), xyz,
                ))
    if not rows:
        raise ValueError("No protein heavy atoms were found in the selected target chains.")
    payload = json.dumps(sorted(rows), separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _protein_residue_inventory(path: Path | str, chains: list[str]) -> tuple[set[tuple[str, int, str, str]], int]:
    """Return author-numbered CA-bearing protein residues and model count."""
    from Bio.PDB.Polypeptide import is_aa

    structure = _read_biopython_structure(path)
    models = list(structure.get_models())
    if len(models) != 1:
        raise ValueError(f"The selected artifact contains {len(models)} models; choose a single-model representation.")
    model = models[0]
    residues: set[tuple[str, int, str, str]] = set()
    for chain_id in chains:
        if chain_id not in model:
            raise ValueError(f"Selected target chain {chain_id!r} is missing from the chosen model.")
        for residue in model[chain_id]:
            if not is_aa(residue, standard=False) or "CA" not in residue:
                continue
            residues.add((str(chain_id), int(residue.id[1]), str(residue.id[2]).strip(), str(residue.resname).strip().upper()))
    if not residues:
        raise ValueError("No CA-bearing protein residues were found in the selected target chains.")
    return residues, len(models)


def validate_supplied_alignment(reference: dict, moving: dict, record: dict) -> dict[str, Any]:
    """Validate an external transform and explicit target-residue map.

    This function does not compute correspondence, fit a transform, or assess
    the scientific fit quality. It binds a user-supplied transform to the
    selected structure/model identities and checks that each listed residue
    exists and the mapping is one-to-one.
    """
    from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer.dataclasses import (
        validate_rigid_transform_matrix,
    )

    errors: list[str] = []
    warnings: list[str] = []
    matrix = None
    try:
        matrix = validate_rigid_transform_matrix(record.get("matrix_4x4"))
    except (TypeError, ValueError) as exc:
        errors.append(f"Transformation: {exc}")

    if str(record.get("reference_identity") or "") != str(reference.get("identity") or ""):
        errors.append("The transformation reference identity does not match the selected structure.")
    if str(record.get("reference_model_id") or "") != str(reference.get("model_id") or ""):
        errors.append("The transformation reference model does not match the selected model.")
    if str(record.get("moving_identity") or "") != str(moving.get("identity") or ""):
        errors.append("The transformation moving identity does not match the selected structure.")
    if str(record.get("moving_model_id") or "") != str(moving.get("model_id") or ""):
        errors.append("The transformation moving model does not match the selected model.")
    if str(record.get("reference_digest") or "") != str(reference.get("structure_digest") or ""):
        errors.append("The transformation reference structure digest does not match the selected artifact.")
    if str(record.get("moving_digest") or "") != str(moving.get("structure_digest") or ""):
        errors.append("The transformation moving structure digest does not match the selected artifact.")
    if not str(reference.get("identity") or "").strip() or not str(moving.get("identity") or "").strip():
        errors.append("Both selected structures need stable identities before an overlay can be shown.")
    if not str(reference.get("model_id") or "").strip() or not str(moving.get("model_id") or "").strip():
        errors.append("Both selected structures need stable model identifiers before an overlay can be shown.")
    if not str(record.get("source") or "").strip():
        errors.append("Transformation source is required.")
    if not str(record.get("method") or "").strip():
        errors.append("Transformation method is required.")

    try:
        reference_residues, _ = _protein_residue_inventory(reference["path"], list(reference.get("target_chains") or []))
        moving_residues, _ = _protein_residue_inventory(moving["path"], list(moving.get("target_chains") or []))
    except (KeyError, OSError, ValueError, ImportError) as exc:
        reference_residues, moving_residues = set(), set()
        errors.append(f"Residue mapping could not be verified: {exc}")

    rows = record.get("residue_correspondence")
    if not isinstance(rows, list) or not rows:
        errors.append("At least one explicit target-residue correspondence is required.")
        rows = []

    seen_reference: set[tuple[str, int, str, str]] = set()
    seen_moving: set[tuple[str, int, str, str]] = set()
    valid_pairs = 0

    def residue_key(value: Any, label: str) -> tuple[str, int, str, str] | None:
        if not isinstance(value, dict):
            errors.append(f"{label} residue must be an object with chain, number, insertion, and residue_name fields.")
            return None
        chain = str(value.get("chain") or "").strip()
        residue_name = str(value.get("residue_name") or "").strip().upper()
        try:
            number = int(value.get("number"))
        except (TypeError, ValueError, OverflowError):
            errors.append(f"{label} residue number must be an integer.")
            return None
        insertion = str(value.get("insertion") or "").strip()
        if not chain or not residue_name:
            errors.append(f"{label} residue needs an explicit chain and residue name.")
            return None
        return chain, number, insertion, residue_name

    for index, row in enumerate(rows, start=1):
        if not isinstance(row, dict):
            errors.append(f"Correspondence row {index} is not an object.")
            continue
        ref_key = residue_key(row.get("reference"), f"Reference row {index}")
        mov_key = residue_key(row.get("moving"), f"Moving row {index}")
        if ref_key is None or mov_key is None:
            continue
        if ref_key[0] not in set(reference.get("target_chains") or []):
            errors.append(f"Reference row {index} is outside the selected target chains.")
        if mov_key[0] not in set(moving.get("target_chains") or []):
            errors.append(f"Moving row {index} is outside the selected target chains.")
        if ref_key in seen_reference:
            errors.append(f"Reference row {index} duplicates a mapped residue; correspondence is ambiguous.")
        if mov_key in seen_moving:
            errors.append(f"Moving row {index} duplicates a mapped residue; correspondence is ambiguous.")
        seen_reference.add(ref_key)
        seen_moving.add(mov_key)
        if ref_key not in reference_residues:
            errors.append(f"Reference row {index} does not identify a residue in the selected structure/model.")
            continue
        if mov_key not in moving_residues:
            errors.append(f"Moving row {index} does not identify a residue in the selected structure/model.")
            continue
        valid_pairs += 1

    if reference_residues:
        missing_reference = len(reference_residues - seen_reference)
    else:
        missing_reference = None
    if moving_residues:
        missing_moving = len(moving_residues - seen_moving)
    else:
        missing_moving = None
    if missing_reference or missing_moving:
        warnings.append("The supplied correspondence is partial; unmapped selected target residues remain unpaired.")
    if not str(record.get("mapping_note") or "").strip():
        warnings.append("No mapping uncertainty/provenance note was supplied; Studio does not estimate mapping uncertainty.")
    warnings.append("Transform form was validated; Studio does not verify alignment fit accuracy or calculate RMSD.")
    warnings.append("Transformation source and mapping notes are user-supplied and are not independently authenticated.")

    return {
        "valid": not errors and valid_pairs > 0,
        "errors": list(dict.fromkeys(errors)),
        "warnings": list(dict.fromkeys(warnings)),
        "matrix_4x4": matrix,
        "correspondence_count": valid_pairs,
        "reference_target_residues": len(reference_residues) if reference_residues else None,
        "moving_target_residues": len(moving_residues) if moving_residues else None,
        "reference_unmapped_residues": missing_reference,
        "moving_unmapped_residues": missing_moving,
        "source": str(record.get("source") or "").strip(),
        "method": str(record.get("method") or "").strip(),
        "mapping_note": str(record.get("mapping_note") or "").strip(),
        "external_provenance": record.get("external_provenance") if isinstance(record.get("external_provenance"), dict) else None,
    }


def _ligand_identity(selection: dict) -> str:
    from rdkit import Chem

    smiles = str(selection.get("smiles") or "").strip() if isinstance(selection, dict) else ""
    molecule = Chem.MolFromSmiles(smiles) if smiles else None
    if molecule is None:
        return ""
    molecule = Chem.RemoveHs(molecule)
    for atom in molecule.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(molecule, isomericSmiles=True, canonical=True)


def compatibility_report(left: dict, right: dict) -> dict[str, Any]:
    """Check protein–protein or protein–ligand comparison compatibility."""
    errors: list[str] = []
    notes: list[str] = []
    if left.get("identity") and str(left.get("identity")) == str(right.get("identity")):
        errors.append("Select two distinct saved structure representations to compare.")
    left_type = str(left.get("complex_type") or "")
    right_type = str(right.get("complex_type") or "")
    if left_type not in {"protein-protein", "protein-ligand"}:
        errors.append("The first structure is not a supported protein complex.")
    if right_type not in {"protein-protein", "protein-ligand"}:
        errors.append("The second structure is not a supported protein complex.")
    if left_type and right_type and left_type != right_type:
        errors.append("The selected structures have different complex types.")
    comparison_type = left_type if left_type == right_type else ""
    for label, item in (("First structure", left), ("Second structure", right)):
        if not item.get("target_chains") or not item.get("target_chains_sequence"):
            errors.append(f"{label} has no selected protein target group with a readable sequence.")
        if comparison_type == "protein-protein":
            if not item.get("binder_chains") or not item.get("binder_chains_sequence"):
                errors.append(f"{label} has no selected protein partner group with a readable sequence.")
        elif comparison_type == "protein-ligand":
            selection = item.get("ligand_selection") if isinstance(item.get("ligand_selection"), dict) else {}
            if not selection.get("chain") or not selection.get("residue_key") or not _ligand_identity(selection):
                errors.append(f"{label} has no uniquely selected ligand with a valid isomeric identity.")
        if not item.get("model_id"):
            notes.append(f"{label} has no saved model identifier; its selected artifact is still kept as an independent sample.")

    if left.get("target_chains_sequence") and right.get("target_chains_sequence"):
        if left["target_chains_sequence"] != right["target_chains_sequence"]:
            errors.append("Selected target chain groups have different sequences; target residues cannot be compared directly.")
    if left.get("target_numbering") and right.get("target_numbering"):
        if left["target_numbering"] != right["target_numbering"]:
            errors.append("Selected target chain groups use different author residue numbering or insertion codes.")
    else:
        errors.append("Target residue numbering could not be verified for both structures.")

    if comparison_type == "protein-ligand":
        first_identity = _ligand_identity(left.get("ligand_selection") or {})
        second_identity = _ligand_identity(right.get("ligand_selection") or {})
        if first_identity and second_identity and first_identity != second_identity:
            errors.append("Selected ligands have different isomeric chemical identities.")
        elif not first_identity or not second_identity:
            errors.append("The selected ligand identity could not be verified for both structures.")

    left_target_id = str(left.get("target_identity") or "").strip()
    right_target_id = str(right.get("target_identity") or "").strip()
    if left_target_id and right_target_id and left_target_id != right_target_id:
        if left.get("target_chains_sequence") == right.get("target_chains_sequence"):
            notes.append("Target labels differ, but the selected target sequences and residue numbering match.")
        else:
            errors.append("Saved target identities differ.")
    elif not left_target_id or not right_target_id:
        notes.append("One or both structures lack a stable target identifier; sequence and numbering were used for compatibility.")

    # A Mol* overlay is meaningful only when the source explicitly records a
    # shared coordinate frame or the underlying coordinate files are identical.
    left_frame = str(left.get("coordinate_frame_id") or "")
    right_frame = str(right.get("coordinate_frame_id") or "")
    coordinate_overlay = bool(left.get("structure_digest") and left.get("structure_digest") == right.get("structure_digest"))
    shared_frame = bool(left_frame and left_frame == right_frame)
    return {
        "compatible": not errors,
        "errors": list(dict.fromkeys(errors)),
        "notes": list(dict.fromkeys(notes)),
        "overlay_compatible": not errors and (coordinate_overlay or shared_frame),
        "overlay_basis": "identical structure coordinates" if coordinate_overlay else "matching saved coordinate-frame identifier" if shared_frame else "",
    }


def fit_matching_targets(reference: dict, moving: dict) -> dict[str, Any]:
    """Explicit single-chain fit; flexible-target RMSD is a comparison only."""
    from mn_molecular_design_studio.complex_analysis_placement import fit_ca_atom_superimposer, PLACEMENT_METHOD

    report = compatibility_report(reference, moving)
    if not report["compatible"]:
        raise ValueError(" ".join(report["errors"]))
    selected_atoms = []
    correspondence = []
    for item in (reference, moving):
        chains = item.get("target_chains") or []
        if len(chains) != 1:
            raise ValueError("Target fitting currently requires one explicitly selected target chain per structure.")
        path = Path(item["path"])
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != item.get("structure_digest"):
            raise ValueError("The selected structure bytes changed; refresh the comparison.")
        models = list(_read_biopython_structure(path).get_models())
        if len(models) != 1:
            raise ValueError("Select a single saved model before target fitting.")
        chain = models[0][chains[0]]
        residues = [residue for residue in chain if residue.id[0] == " "]
        if not residues or any("CA" not in residue for residue in residues):
            raise ValueError("Every selected target residue must have a saved C-alpha atom.")
        correspondence.append([(residue.resname, residue.id[1], residue.id[2]) for residue in residues])
        selected_atoms.append([residue["CA"] for residue in residues])
    if correspondence[0] != correspondence[1]:
        raise ValueError("Target fitting requires identical residue identities and author numbering.")
    fitter, matrix = fit_ca_atom_superimposer(*selected_atoms, require_non_collinear=True)
    return {"schema": "studio.target-fit-comparison.v1", "method": PLACEMENT_METHOD,
            "target_ca_rmsd_angstrom": float(fitter.rms), "mapped_residues": len(selected_atoms[0]),
            "matrix_4x4": matrix,
            "reference": {key: reference.get(key) for key in ("identity", "model_id", "structure_digest", "target_chains")},
            "moving": {key: moving.get(key) for key in ("identity", "model_id", "structure_digest", "target_chains")},
            "interpretation": "Rigid target C-alpha superposition of saved models; not native-state accuracy or experimental binding evidence."}


def saved_target_contact_rows(contact_sets: dict[str, list[dict]]) -> list[dict[str, Any]]:
    """Make a residue-by-structure presence table from saved contact artifacts.

    This only transcribes existing contact records; it does not calculate new
    contacts, distances, similarities, or rankings.
    """
    residues: set[str] = set()
    normalized: dict[str, set[str]] = {}
    for label, pairs in contact_sets.items():
        selected = set()
        for pair in pairs:
            if not isinstance(pair, dict):
                continue
            residue = str(pair.get("target_residue") or "").strip()
            if residue:
                selected.add(residue)
                residues.add(residue)
        normalized[str(label)] = selected
    return [{"Target residue": residue, **{label: residue in values for label, values in normalized.items()}}
            for residue in sorted(residues)]
