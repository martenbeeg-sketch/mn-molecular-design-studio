"""Recover documented decisions for the exact source spreadsheet and row."""
from __future__ import annotations

import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

from mn_molecular_design_studio import molecule_imports

_CANDIDATE_FIELDS = (
    "cid", "smiles", "molecular_formula", "inchi_key", "vendor_formula", "formula_comparison",
    "pubchem_record_smiles", "pubchem_record_formula", "pubchem_record_fragment_count",
    "selected_component_smiles", "selected_component_formula", "selected_component_occurrences",
    "component_selection", "validation_engine", "rdkit_version",
)


def pubchem_review(candidate, **evidence):
    return {"method": "PubChem", "candidate": {key: candidate[key] for key in _CANDIDATE_FIELDS if key in candidate},
            **evidence}


def _object(path, root):
    try:
        if path.resolve().is_relative_to(root) and path.is_file():
            data = json.loads(path.read_text(encoding="utf-8"))
            if isinstance(data, dict):
                return data
    except (OSError, ValueError):
        pass
    return {}


def archived_reviews(source_bytes, audit, roots):
    rows = defaultdict(list)
    for row in audit["rows"]:
        if row["status"] in {"invalid", "empty"} and row["source_id"]:
            rows[(row["source_id"], row["raw_smiles"])].append(row["row_number"])
    if not rows:
        return {}
    mapping = audit["mapping"]
    by_number = defaultdict(list)
    digest = hashlib.sha256(source_bytes).hexdigest()
    for root in map(Path.resolve, roots):
        for run in sorted((root / "compound-import").glob("*")):
            reports = run / "artifacts/reports"
            rejected = reports / "rejected_compounds.csv"
            if not rejected.is_file() or not rejected.resolve().is_relative_to(root):
                continue
            validation = _object(reports / "validation_report.json", root)
            if any(validation.get(key) != mapping.get(key) for key in ("sheet_name", "smiles_column", "id_column")):
                continue
            for source in (run / "artifacts/compounds").glob("*"):
                if source.resolve().is_relative_to(root) and source.is_file() and source.stat().st_size == len(source_bytes):
                    with source.open("rb") as handle:
                        if hashlib.file_digest(handle, "sha256").hexdigest() == digest:
                            break
            else:
                continue
            with rejected.open(newline="", encoding="utf-8-sig") as handle:
                rejected_rows = {(row["compound_id"], int(row["source_row"])): row["smiles"].strip()
                                 for row in csv.DictReader(handle)}
            latest = {}
            for review_dir in sorted((root / "compound-review").glob("*")):
                meta = _object(review_dir / "metadata.json", root)
                if meta.get("source_compound_run_id") != run.name or meta.get("status") != "completed":
                    continue
                report = _object(review_dir / "artifacts/reports/pubchem_review.json", root)
                key = (report.get("compound_id"), report.get("source_row"))
                if report.get("source_run_id") != run.name or key not in rejected_rows:
                    continue
                timestamp = str(report.get("reviewed_at") or "")
                if key not in latest or timestamp > latest[key][0]:
                    latest[key] = (timestamp, report, review_dir.name)
            for (vendor_id, legacy_row), (timestamp, report, review_id) in latest.items():
                numbers = rows.get((vendor_id, rejected_rows[(vendor_id, legacy_row)]), [])
                if len(numbers) != 1:
                    continue
                candidate = report.get("candidate") or {}
                if report.get("decision") == "accepted" and candidate.get("smiles") and candidate.get("cid"):
                    by_number[numbers[0]].append(pubchem_review(candidate, source="MN-Ligand",
                        source_run_id=run.name, source_row=legacy_row, review_run_id=review_id, reviewed_at=timestamp,
                        query=report.get("query"), query_type=report.get("query_type")))
                else:
                    by_number[numbers[0]].append(None)
    return {number: evidence[0] for number, evidence in by_number.items()
            if evidence[0] and all(item and item["candidate"] == evidence[0]["candidate"] for item in evidence)}


def apply_reviews(source, audit, reviews, *, manual=None, pasted=False):
    """Revalidate accepted replacements with the same analyzer as the source file."""
    original = {row["row_number"]: row for row in audit["rows"] if row["status"] in {"invalid", "empty"}}
    replacements = {number: evidence["candidate"]["smiles"] for number, evidence in reviews.items() if number in original}
    replacements.update({number: smiles for number, smiles in (manual or {}).items() if number in original})
    if not replacements:
        return None
    if pasted:
        revised = molecule_imports.analyze_smiles_text(source.decode("utf-8-sig"), corrections=replacements)
        revised.update(source_name=audit["source_name"], source_sha256=audit["source_sha256"])
    else:
        mapping = dict(audit["mapping"])
        previous = {int(number): smiles for number, smiles in mapping.pop("corrections", {}).items()}
        revised = molecule_imports.analyze_small_molecule_table(source, audit["source_name"], **mapping,
                                                              corrections={**previous, **replacements})
    for row in revised["rows"]:
        number = row["row_number"]
        if number not in replacements or row["status"] not in {"valid", "duplicate"}:
            continue
        evidence = {"method": "Manual replacement"} if number in (manual or {}) else reviews[number]
        row["review"] = {**evidence, "original_status": original[number]["status"],
                         "original_errors": original[number]["errors"], "replacement_smiles": replacements[number]}
        if evidence.get("method") == "PubChem":
            row["warnings"] = [message for message in row["warnings"] if "user-supplied replacement" not in message]
            row["warnings"].append("Revalidated the approved PubChem structure for this row.")
    revised["summary"]["reviewed_rows"] = sum(bool(row.get("review")) for row in revised["rows"])
    return revised
