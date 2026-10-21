"""Typed, portable compound-set records shared by ligand tiles and pipelines."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import re
from pathlib import Path
from typing import Any, Iterable


COMPOUND_SET_SCHEMA_VERSION = 2
_SMILES_FIELD_NAMES = (
    "smiles", "canonical_smiles", "canonical_isomeric_smiles",
    "isomeric_smiles", "canonical smiles",
)
_ID_FIELD_NAMES = ("compound_id", "id", "pubchem_cid", "zinc_id", "catalog_id", "compound", "name")
_NAME_FIELD_NAMES = ("compound_name", "name", "compound", "preferred_name", "common_name", "title")


def _clean(value: object) -> str:
    return str(value or "").strip()


def _safe_id(value: object, index: int, smiles: str) -> str:
    text = _clean(value)
    if not text:
        text = f"cmp-{hashlib.sha256(smiles.encode('utf-8')).hexdigest()[:12]}-{index:06d}"
    text = re.sub(r"[^A-Za-z0-9._-]+", "-", text).strip("-._")
    return (text or f"compound-{index:06d}")[:100]


def normalize_compounds(rows: Iterable[dict[str, Any]]) -> list[dict[str, Any]]:
    """Validate IDs and SMILES while preserving all input properties."""
    normalized: list[dict[str, Any]] = []
    seen_ids: set[str] = set()
    for index, raw in enumerate(rows, start=1):
        if not isinstance(raw, dict):
            raise ValueError(f"Compound row {index} must be a table-like object.")
        smiles = _clean(
            raw.get("smiles") or raw.get("canonical_smiles")
            or raw.get("canonical_isomeric_smiles") or raw.get("isomeric_smiles")
        )
        if not smiles:
            raise ValueError(f"Compound row {index} has no SMILES value.")
        properties = raw.get("properties")
        if not isinstance(properties, dict):
            properties = {}
        compound_id = _safe_id(raw.get("compound_id") or raw.get("id") or raw.get("name"), index, smiles)
        compound_name = _clean(raw.get("compound_name") or raw.get("name") or properties.get("name")) or compound_id
        folded_id = compound_id.casefold()
        if folded_id in seen_ids:
            raise ValueError(f"Compound ID {compound_id!r} occurs more than once.")
        seen_ids.add(folded_id)
        properties = raw.get("properties")
        if not isinstance(properties, dict):
            properties = {
                str(key): value
                for key, value in raw.items()
                if key not in {"compound_id", "id", "name", "compound_name", "smiles", "canonical_smiles", "canonical_isomeric_smiles", "isomeric_smiles"}
                and value not in (None, "")
            }
        normalized.append({
            "compound_id": compound_id,
            "compound_name": compound_name,
            "smiles": smiles,
            "properties": properties,
        })
    if not normalized:
        raise ValueError("The compound set contains no rows with SMILES.")
    if len(normalized) > 100_000:
        raise ValueError("A compound set may contain at most 100,000 molecules.")
    return normalized


def compound_content_sha256(compounds: list[dict[str, Any]]) -> str:
    """Return a stable digest for a normalized molecule collection."""
    normalized = normalize_compounds(compounds)
    canonical = json.dumps(normalized, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def parse_smi_text(text: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for index, raw_line in enumerate(str(text or "").splitlines(), start=1):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(maxsplit=1)
        smiles = parts[0].strip()
        compound_id = parts[1].strip() if len(parts) > 1 else f"compound-{index:06d}"
        rows.append({"compound_id": compound_id, "compound_name": compound_id, "smiles": smiles})
    return normalize_compounds(rows)


def parse_compound_table(
    data: bytes,
    *,
    filename: str,
    id_column: str = "",
    name_column: str = "",
    smiles_column: str = "",
) -> list[dict[str, Any]]:
    """Read CSV, TSV, SMI, or SDF records into the shared compound contract."""
    name = Path(filename).name
    suffix = Path(name).suffix.lower()
    text = data.decode("utf-8-sig", errors="replace")
    if suffix in {".smi", ".smiles"}:
        return parse_smi_text(text)
    if suffix in {".sdf", ".mol"}:
        return _parse_sdf(text)
    if suffix not in {".csv", ".tsv", ".txt"}:
        raise ValueError("Compound input must be CSV, TSV, SMI, or SDF.")
    delimiter = "\t" if suffix in {".tsv", ".txt"} else ","
    reader = csv.DictReader(io.StringIO(text), delimiter=delimiter)
    fields = [str(value) for value in (reader.fieldnames or [])]
    if not fields:
        raise ValueError("The compound table has no header row.")
    by_normalized = {re.sub(r"[\s-]+", "_", key.strip().lower()): key for key in fields}
    selected_smiles = smiles_column or next(
        (by_normalized[key] for key in _SMILES_FIELD_NAMES if key in by_normalized), ""
    )
    selected_id = id_column or next(
        (by_normalized[key] for key in _ID_FIELD_NAMES if key in by_normalized), ""
    )
    selected_name = name_column or next(
        (by_normalized[key] for key in _NAME_FIELD_NAMES if key in by_normalized), ""
    )
    if selected_smiles not in fields:
        raise ValueError("Select a SMILES column to import this compound table.")
    rows = []
    for raw in reader:
        smiles = _clean(raw.get(selected_smiles))
        if not smiles:
            continue
        row: dict[str, Any] = {"smiles": smiles}
        if selected_id and selected_id in fields:
            row["compound_id"] = _clean(raw.get(selected_id))
        if selected_name and selected_name in fields:
            row["compound_name"] = _clean(raw.get(selected_name))
        row["properties"] = {
            str(key): value
            for key, value in raw.items()
            if key not in {selected_smiles, selected_id, selected_name} and value not in (None, "")
        }
        rows.append(row)
    return normalize_compounds(rows)


def _parse_sdf(text: str) -> list[dict[str, Any]]:
    """Read SDF properties; when available use RDKit for molecules without SMILES tags."""
    records = [record for record in re.split(r"\$\$\$\$\s*(?:\r?\n|$)", text) if record.strip()]
    rows: list[dict[str, Any]] = []
    rdkit_chem = None
    for index, record in enumerate(records, start=1):
        lines = record.splitlines()
        compound_id = _clean(lines[0] if lines else "") or f"compound-{index:06d}"
        properties: dict[str, str] = {}
        for position, line in enumerate(lines):
            match = re.match(r">\s*<([^>]+)>\s*$", line.strip())
            if match and position + 1 < len(lines):
                properties[match.group(1).strip()] = lines[position + 1].strip()
        smiles = next(
            (value for key, value in properties.items() if key.strip().lower() in _SMILES_FIELD_NAMES),
            "",
        )
        if not smiles:
            if rdkit_chem is None:
                try:
                    from rdkit import Chem as rdkit_chem  # type: ignore[assignment]
                except ImportError:
                    rdkit_chem = False
            if rdkit_chem:
                molecule = rdkit_chem.MolFromMolBlock(record, sanitize=True, removeHs=False)
                if molecule is not None:
                    smiles = rdkit_chem.MolToSmiles(molecule, canonical=True)
        if not smiles:
            raise ValueError(
                f"SDF record {index} has no SMILES property and could not be converted. "
                "Add a SMILES field to the SDF or import a CSV/SMI file."
            )
        rows.append({"compound_id": compound_id, "compound_name": compound_id, "smiles": smiles, "properties": properties})
    return normalize_compounds(rows)


def write_compound_set(
    path: Path,
    *,
    name: str,
    compounds: list[dict[str, Any]],
    owner_id: str,
    project_id: str,
    created_at: str,
    compound_set_id: str = "",
) -> dict[str, Any]:
    normalized = normalize_compounds(compounds)
    content_sha256 = compound_content_sha256(normalized)
    set_id = _clean(compound_set_id) or f"cs-{content_sha256[:24]}"
    if not re.fullmatch(r"[A-Za-z0-9._-]{1,100}", set_id):
        raise ValueError("Compound-set ID contains unsupported characters.")
    payload = {
        "schema_version": COMPOUND_SET_SCHEMA_VERSION,
        "record_type": "compound_set",
        "compound_set_id": set_id,
        "content_sha256": content_sha256,
        "name": str(name).strip() or "Compound set",
        "owner_id": str(owner_id),
        "project_id": str(project_id),
        "created_at": str(created_at),
        "compound_count": len(normalized),
        "compounds": normalized,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)
    return payload


def load_compound_set(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"Could not read compound set {path}: {exc}") from exc
    if not isinstance(payload, dict) or payload.get("record_type") != "compound_set":
        raise ValueError("The selected database artifact is not a compound set.")
    schema_version = payload.get("schema_version")
    if schema_version not in (1, COMPOUND_SET_SCHEMA_VERSION):
        raise ValueError(f"Unsupported compound-set schema version: {payload.get('schema_version')!r}.")
    original_compounds = payload.get("compounds") or []
    payload["compounds"] = normalize_compounds(original_compounds)
    content_sha256 = compound_content_sha256(payload["compounds"])
    recorded_digest = str(payload.get("content_sha256") or "")
    if recorded_digest and recorded_digest != content_sha256 and schema_version == 1:
        legacy = [
            {key: row[key] for key in ("compound_id", "smiles", "properties") if key in row}
            for row in payload["compounds"]
        ]
        legacy_canonical = json.dumps(legacy, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
        legacy_digest = hashlib.sha256(legacy_canonical.encode("utf-8")).hexdigest()
        if recorded_digest == legacy_digest:
            recorded_digest = content_sha256
    if recorded_digest and recorded_digest != content_sha256:
        raise ValueError("The compound set failed its content-integrity check.")
    payload["schema_version"] = COMPOUND_SET_SCHEMA_VERSION
    payload["content_sha256"] = content_sha256
    payload.setdefault("compound_set_id", f"cs-{content_sha256[:24]}")
    payload["compound_count"] = len(payload["compounds"])
    return payload
