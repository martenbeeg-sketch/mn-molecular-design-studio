"""Preparation evidence attached to one exact coordinate representation."""
from __future__ import annotations

import hashlib
import io
import json
from pathlib import Path
import re
from zipfile import ZIP_DEFLATED, ZipFile

MANIFEST = "structure-preparation.json"
MAX_MANIFEST_BYTES = 256_000
COORDINATE_SUFFIXES = {".pdb", ".ent", ".cif", ".mmcif"}


def sha256(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def validate_history(history, expected_sha256):
    if not isinstance(history, dict) or history.get("version") != 1:
        raise ValueError("Unsupported structure preparation history.")
    if not re.fullmatch(r"[0-9a-f]{64}", str(expected_sha256)) or history.get("structure_sha256") != expected_sha256:
        raise ValueError("The preparation history checksum does not match this structure.")
    steps = history.get("steps")
    if not isinstance(steps, list) or not steps or any(
            not isinstance(step, dict) or not isinstance(step.get("label"), str) for step in steps):
        raise ValueError("Preparation history must describe the recorded steps.")
    artifacts = history.get("artifacts", [])
    if not isinstance(artifacts, list) or any(not isinstance(item, dict) for item in artifacts):
        raise ValueError("Preparation artifact references must be a list of records.")
    if len(json.dumps(history, allow_nan=False).encode()) > MAX_MANIFEST_BYTES:
        raise ValueError("The preparation history is larger than 256 KB.")
    return history


def parse_history_text(content, expected_sha256=None):
    if not isinstance(content, str) or len(content.encode()) > MAX_MANIFEST_BYTES:
        raise ValueError("The preparation history must be a JSON document of at most 256 KB.")
    history = json.loads(content)
    digest = history.get("structure_sha256") if isinstance(history, dict) else None
    return validate_history(history, expected_sha256 if expected_sha256 is not None else digest)


def extend_history(parent, source_sha256, output_sha256, steps, *, application="Molecular Design Studio", artifact=None):
    previous = validate_history(parent, source_sha256) if parent else {}
    result = {
        "version": 1, "structure_sha256": output_sha256, "source_sha256": source_sha256,
        "verification": "studio_preparation", "application": application,
        "steps": [*previous.get("steps", []), *steps],
        "artifacts": [*previous.get("artifacts", []), *([artifact] if artifact else [])],
    }
    return validate_history(result, output_sha256)


def portable_bundle(filename, content, history):
    filename = Path(filename).name
    if Path(filename).suffix.lower() not in COORDINATE_SUFFIXES:
        raise ValueError("A preparation bundle requires PDB or mmCIF coordinates.")
    validate_history(history, sha256(content))
    stream = io.BytesIO()
    with ZipFile(stream, "w", ZIP_DEFLATED) as archive:
        archive.writestr(filename, content)
        archive.writestr(MANIFEST, json.dumps(history, indent=2, allow_nan=False))
    return stream.getvalue()


def report_steps(report):
    """Summarize executed Studio preparation steps using the native report."""
    result = []

    def add(label, details=None):
        result.append({"label": label, **({"details": details} if details else {})})

    repair = report.get("protein_repair") or {}
    parameters = report.get("parameters") or {}
    if parameters.get("repair") is True and repair.get("enabled") is not False:
        if replacements := (repair.get("modeller_noncanonical_repair") or {}).get("replacements"):
            add("MODELLER residue repair", {"replacements": [
                {key: row.get(key) for key in ("key", "target")} for row in replacements]})
        if gaps := (repair.get("modeller_internal_gap_repair") or {}).get("modeled_gaps"):
            add("MODELLER gap reconstruction", _gap_details(gaps))
        fixer = repair.get("pdbfixer") or {}
        add("PDBFixer cleaning", {f"{key}_count": len(fixer[key]) for key in (
            "missing_residue_segments_added", "missing_atoms") if isinstance(fixer.get(key), list)})
        refinement = fixer.get("refinement") or {}
        if refinement.get("engine") or refinement.get("potential_energy_after_kj_mol") is not None:
            add("OpenMM minimization", {key: refinement[key] for key in (
                "engine", "forcefield", "ligand_context", "restraint_strategy",
                "potential_energy_before_kj_mol", "potential_energy_after_kj_mol")
                                       if key in refinement})
    if extension := report.get("terminal_extension"):
        if extension.get("enabled") is not False and parameters.get("extension_sequence"):
            add("MODELLER C-terminal extension", {"chain": parameters.get("extension_chain"),
                                                   "sequence": parameters["extension_sequence"]})
    if gaps := (report.get("confirmed_gap_repair") or {}).get("modeled_gaps"):
        add("MODELLER confirmed gap reconstruction", _gap_details(gaps))
    loop = report.get("loop_remodelling") or {}
    if loop.get("enabled"):
        add("MODELLER existing-loop remodelling",
            {key: loop[key] for key in ("chain", "start", "end", "models", "selected_model") if key in loop})
    substitutions = report.get("residue_substitutions") or {}
    if substitutions.get("enabled"):
        add("MODELLER amino-acid substitutions", {"replacements": [
            {key: row.get(key) for key in ("chain", "resseq", "icode", "resname", "target")}
            for row in substitutions.get("replacements", [])]})
    trim = report.get("trimming") or {}
    if trim.get("enabled"):
        add("Chain trimming", {"ranges": json.loads(parameters.get("trim_ranges_json") or "{}")})
    if not result:
        add("Complex preparation")
    return result


def _gap_details(gaps):
    return {"gaps": [{key: gap.get(key) for key in ("chain", "author_start", "author_end", "sequence")}
                     for gap in gaps]}
