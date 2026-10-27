"""Validated molecular imports shared by the Database UI and API."""
from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import math
import re
import zipfile
from pathlib import Path
from urllib.request import Request, urlopen

MAX_ROWS = 100_000
MAX_BYTES = 100 * 1024 * 1024
KINDS = ("protein", "small_molecule", "dna", "rna")
ALPHABETS = {
    "protein": set("ACDEFGHIKLMNPQRSTVWYXBZJUO"),
    "dna": set("ACGTRYSWKMBDHVN"),
    "rna": set("ACGURYSWKMBDHVN"),
}


def canonical_component(component: dict) -> dict:
    kind = str(component.get("type") or "protein").lower()
    if kind not in KINDS:
        raise ValueError(f"Unknown chain type: {kind}.")
    result = {"type": kind, "chain": str(component.get("chain") or "A"),
              "role": str(component.get("role") or "")}
    if kind == "small_molecule":
        smiles = str(component.get("smiles") or "").strip()
        if smiles:
            from rdkit import Chem
            molecule = Chem.MolFromSmiles(smiles)
            if molecule is None:
                raise ValueError(f"Invalid SMILES for chain {result['chain']}.")
            result["smiles"] = Chem.MolToSmiles(molecule, isomericSmiles=True)
        elif component.get("ccd") and re.fullmatch(r"[A-Za-z0-9_-]{1,20}", str(component["ccd"])):
            result["ccd"] = str(component["ccd"]).upper()
        else:
            raise ValueError("A small molecule needs valid SMILES or a CCD component code.")
    else:
        sequence = re.sub(r"\s+", "", str(component.get("sequence") or "")).upper().rstrip("*")
        invalid = sorted(set(sequence) - ALPHABETS[kind])
        if not sequence or invalid:
            raise ValueError(f"Invalid {kind} sequence for chain {result['chain']}" +
                             (f": {', '.join(invalid)}." if invalid else "."))
        result["sequence"] = sequence
    if component.get("residue_key"):
        result["residue_key"] = str(component["residue_key"])
    return result


def canonical_components(components: list[dict]) -> tuple[list[dict], str, str]:
    if not components:
        raise ValueError("A molecule must have at least one chain.")
    rows = [canonical_component(row) for row in components]
    if len({row["chain"] for row in rows}) != len(rows):
        raise ValueError("Chain IDs must be unique within a molecule.")
    identities = sorted((row["type"], row.get("sequence") or row.get("smiles") or "ccd:" + row["ccd"])
                        for row in rows)
    fingerprint = hashlib.sha256(json.dumps(identities, separators=(",", ":")).encode()).hexdigest()
    types = {row["type"] for row in rows}
    kind = next(iter(types)) if len(types) == 1 else "complex"
    return rows, fingerprint, kind


def parse_sequences(text: str, *, kind: str = "protein") -> list[dict]:
    records, current, parts = [], None, []
    for raw in str(text).splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current is not None:
                records.append((current, "".join(parts)))
            current, parts = line[1:].strip() or f"Molecule {len(records) + 1}", []
        elif current is None:
            records.append((f"Molecule {len(records) + 1}", line))
        else:
            parts.append(line)
    if current is not None:
        records.append((current, "".join(parts)))
    if not records:
        raise ValueError("Enter at least one sequence.")
    if len(records) > MAX_ROWS:
        raise ValueError("An import can contain at most 100,000 molecules.")
    result = []
    for name, sequence in records:
        components = [canonical_component({"type": kind, "chain": chr(65 + index) if index < 26 else f"C{index+1}",
                                           "sequence": value}) for index, value in enumerate(sequence.split(":"))]
        result.append({"name": name, "components": components, "metadata": {}})
    return result


def parse_smiles(text: str) -> list[dict]:
    rows = []
    for line in str(text).splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        values = line.strip().split(maxsplit=1)
        rows.append({"name": values[1] if len(values) > 1 else f"Compound {len(rows)+1}",
                     "components": [canonical_component({"type": "small_molecule", "chain": "L", "smiles": values[0]})],
                     "metadata": {}})
    if not rows:
        raise ValueError("Enter at least one SMILES string.")
    if len(rows) > MAX_ROWS:
        raise ValueError("An import can contain at most 100,000 molecules.")
    return rows


def read_table(data: bytes, filename: str, *, sheet_name: str = "", string_values=False):
    import pandas as pd
    if len(data) > MAX_BYTES:
        raise ValueError("The table exceeds 100 MB.")
    suffix = Path(filename).suffix.lower()
    if suffix in {".xlsx", ".xls", ".xlsm"}:
        table = pd.read_excel(io.BytesIO(data), sheet_name=sheet_name or 0, keep_default_na=False,
                              dtype=str if string_values else None)
    else:
        table = pd.read_csv(io.BytesIO(data), sep="\t" if suffix in {".tsv", ".txt"} else ",",
                            keep_default_na=False, dtype=str if string_values else None)
    table.columns = [str(value) for value in table.columns]
    if len(set(table.columns)) != len(table.columns):
        raise ValueError("Table columns must have unique names.")
    if len(table) > MAX_ROWS:
        raise ValueError("An import can contain at most 100,000 molecules.")
    return table


def table_sheets(data: bytes, filename: str) -> tuple[str, ...]:
    """List worksheet names for supported Excel workbooks."""
    import pandas as pd
    if len(data) > MAX_BYTES:
        raise ValueError("The table exceeds 100 MB.")
    if Path(filename).suffix.lower() not in {".xlsx", ".xls", ".xlsm"}:
        return ()
    try:
        return tuple(str(name) for name in pd.ExcelFile(io.BytesIO(data)).sheet_names)
    except Exception as exc:
        raise ValueError(f"Excel workbook could not be read: {exc}") from exc


def parse_table(
    data: bytes,
    filename: str,
    *,
    chain_columns: dict[str, dict],
    name_column: str = "",
    id_column: str = "",
    sheet_name: str = "",
) -> list[dict]:
    table = read_table(data, filename, sheet_name=sheet_name)
    if not chain_columns:
        raise ValueError("Map at least one chain column.")
    if (
        any(column not in table.columns for column in chain_columns)
        or (name_column and name_column not in table.columns)
        or (id_column and id_column not in table.columns)
    ):
        raise ValueError("A mapped column is missing from the table.")
    result = []
    for index, source in enumerate(table.to_dict("records"), 1):
        components = []
        for column, mapping in chain_columns.items():
            value = str(source[column]).strip()
            if not value:
                continue
            field = "smiles" if mapping["type"] == "small_molecule" else "sequence"
            try:
                components.append(canonical_component({**mapping, field: value}))
            except ValueError as exc:
                raise ValueError(f"Row {index}, column {column}: {exc}") from exc
        if not components:
            raise ValueError(f"Row {index} has no molecular chains.")
        canonical_components(components)
        molecule_id = str(source.get(id_column) or "").strip() if id_column else ""
        molecule_name = str(source.get(name_column) or "").strip() if name_column else ""
        molecule_name = molecule_name or molecule_id or f"Molecule {index}"
        metadata = {
            key: value for key, value in source.items()
            if key not in chain_columns and key not in {name_column, id_column}
        }
        if molecule_id:
            if any(component["type"] == "small_molecule" for component in components):
                metadata["compound_id"] = molecule_id
            else:
                metadata["external_id"] = molecule_id
        result.append({"name": molecule_name, "components": components, "metadata": metadata})
    if not result:
        raise ValueError("The table contains no molecules.")
    return result


def _small_molecule_analysis_row(molecule):
    """Describe a sanitized RDKit molecule without changing its chemistry."""
    from rdkit import Chem
    from rdkit.Chem import Crippen, Descriptors, Lipinski, QED, rdMolDescriptors

    canonical = Chem.MolToSmiles(molecule, isomericSmiles=True)
    connectivity = Chem.Mol(molecule)
    Chem.RemoveStereochemistry(connectivity)
    connectivity_smiles = Chem.MolToSmiles(connectivity, isomericSmiles=False)
    fragments = len(Chem.GetMolFrags(molecule))
    charge = int(Chem.GetFormalCharge(molecule))
    component_molecules = Chem.GetMolFrags(molecule, asMols=True)
    centers = Chem.FindMolChiralCenters(molecule, includeUnassigned=True)
    warnings = []
    if fragments > 1:
        warnings.append(f"Contains {fragments} disconnected fragments; salts and mixtures are retained.")
    if charge:
        warnings.append(f"Net formal charge is {charge:+d}; no neutralization was applied.")
    unspecified = sum(label == "?" for _, label in centers)
    if unspecified:
        warnings.append(f"{unspecified} stereocenter(s) have unspecified stereochemistry.")
    if any(Chem.GetFormalCharge(component) for component in component_molecules) and not charge:
        warnings.append("One or more disconnected components are charged, although the net molecular charge is zero.")
    descriptors = {
        "formula": rdMolDescriptors.CalcMolFormula(molecule),
        "molecular_weight": round(float(Descriptors.MolWt(molecule)), 4),
        "exact_mass": round(float(rdMolDescriptors.CalcExactMolWt(molecule)), 4),
        "heavy_atoms": int(molecule.GetNumHeavyAtoms()),
        "clogp": round(float(Crippen.MolLogP(molecule)), 4),
        "tpsa": round(float(rdMolDescriptors.CalcTPSA(molecule)), 4),
        "hbd": int(Lipinski.NumHDonors(molecule)),
        "hba": int(Lipinski.NumHAcceptors(molecule)),
        "rotatable_bonds": int(Lipinski.NumRotatableBonds(molecule)),
        "ring_count": int(Lipinski.RingCount(molecule)),
        "formal_charge": charge,
        "fraction_csp3": round(float(rdMolDescriptors.CalcFractionCSP3(molecule)), 4),
        "qed": round(float(QED.qed(molecule)), 4),
        "fragment_count": fragments,
        "components": [{
            "component": index,
            "smiles": Chem.MolToSmiles(component, canonical=True, isomericSmiles=True),
            "formula": rdMolDescriptors.CalcMolFormula(component),
            "molecular_weight": round(float(Descriptors.MolWt(component)), 4),
            "heavy_atoms": int(component.GetNumHeavyAtoms()),
            "formal_charge": int(Chem.GetFormalCharge(component)),
            "contains_carbon": any(atom.GetAtomicNum() == 6 for atom in component.GetAtoms()),
        } for index, component in enumerate(component_molecules, 1)],
    }
    if any(isinstance(value, float) and not math.isfinite(value) for value in descriptors.values()):
        raise ValueError("RDKit produced a non-finite molecular descriptor")
    return canonical, connectivity_smiles, fragments, charge, warnings, descriptors


def _annotation_stripped_candidate(source_smiles: str) -> str:
    """Suggest removal of known vendor dot-annotations without changing input."""
    parts = []
    current = []
    bracket_depth = 0
    for character in source_smiles:
        if character == "[":
            bracket_depth += 1
        elif character == "]" and bracket_depth:
            bracket_depth -= 1
        if character == "." and bracket_depth == 0:
            parts.append("".join(current))
            current = []
        else:
            current.append(character)
    parts.append("".join(current))
    retained = []
    for part in parts:
        token = part.strip()
        if re.fullmatch(
            r"\[(?:Z|E|Rotation|\(\+\)|\(-\)|\+|-|"
            r"(?:\d+(?:\.\d+)?|\d+/\d+)"
            r"(?:\s+[A-Za-z0-9+\-]+)?)\]",
            token,
            flags=re.IGNORECASE,
        ):
            continue
        retained.append(token)
    return ".".join(retained)


def analyze_small_molecule_table(
    data: bytes,
    filename: str,
    *,
    smiles_column: str,
    name_column: str = "",
    id_column: str = "",
    sheet_name: str = "",
) -> dict:
    """Validate every row of a one-SMILES-per-row table and report exact duplicates.

    SMILES are parsed and sanitized with RDKit, then written back in canonical
    isomeric form. The source structure is never desalted, neutralized, or
    stripped of stereochemistry. Exact structure duplicates are represented by
    their first valid source row in ``records``; every source row remains in the
    accompanying audit report.
    """
    from rdkit import Chem, rdBase

    table = read_table(data, filename, sheet_name=sheet_name, string_values=True)
    columns = set(table.columns)
    if smiles_column not in columns or (name_column and name_column not in columns) or (id_column and id_column not in columns):
        raise ValueError("A mapped SMILES, name, or ID column is missing from the table.")
    if not table.shape[0]:
        raise ValueError("The table contains no molecule rows.")

    batch_id = hashlib.sha256(
        Path(filename).name.encode("utf-8") + b"\0" + data + b"\0" +
        json.dumps({"smiles_column": smiles_column, "name_column": name_column,
                    "id_column": id_column, "sheet_name": sheet_name}, sort_keys=True).encode("utf-8")
    ).hexdigest()[:32]
    result_records = []
    report_rows = []
    canonical_first = {}
    connectivity_first = {}
    used_ids = {}
    counts = {"rows": len(table), "valid_rows": 0, "unique_structures": 0, "invalid_rows": 0,
              "empty_rows": 0, "duplicate_rows": 0, "stereo_variants": 0,
              "multi_fragment_rows": 0, "charged_rows": 0, "charged_component_rows": 0,
              "duplicate_id_rows": 0}

    excluded = {smiles_column}
    excluded.update(value for value in (name_column, id_column) if value)
    for row_number, source in enumerate(table.to_dict("records"), 1):
        raw = str(source.get(smiles_column, "") or "").strip()
        source_id = str(source.get(id_column, "") or "").strip() if id_column else ""
        base_id = source_id or f"compound_{row_number:07d}"
        used_ids[base_id] = used_ids.get(base_id, 0) + 1
        compound_id = base_id if used_ids[base_id] == 1 else f"{base_id}__{used_ids[base_id]}"
        if used_ids[base_id] > 1:
            counts["duplicate_id_rows"] += 1
        name = str(source.get(name_column, "") or "").strip() if name_column else ""
        name = name or compound_id
        row = {
            "row_number": row_number,
            "source_id": source_id,
            "compound_id": compound_id,
            "name": name,
            "raw_smiles": raw,
            "canonical_smiles": "",
            "connectivity_smiles": "",
            "status": "invalid",
            "duplicate_of_row": None,
            "stereo_variant_of_row": None,
            "warnings": [],
            "errors": [],
            "descriptors": {},
            "annotation_stripped_candidate": "",
            "candidate_status": "",
            "source_data": {},
            "molecule_id": None,
            "matched_existing": False,
        }
        if not raw:
            row["status"] = "empty"
            row["errors"].append("SMILES is blank.")
            row["source_data"] = source
            counts["empty_rows"] += 1
            report_rows.append(row)
            continue

        try:
            with rdBase.BlockLogs():
                molecule = Chem.MolFromSmiles(raw, sanitize=True)
            if molecule is None:
                raise ValueError("RDKit could not parse and sanitize the SMILES")
            if molecule.GetNumHeavyAtoms() < 1:
                raise ValueError("SMILES contains no heavy atoms")
            if any(atom.GetAtomicNum() == 0 for atom in molecule.GetAtoms()):
                raise ValueError("SMILES contains dummy or query atoms")
            canonical, connectivity, fragment_count, charge, warnings, descriptors = _small_molecule_analysis_row(molecule)
            if not canonical:
                raise ValueError("RDKit could not generate canonical SMILES")
            row.update(canonical_smiles=canonical, connectivity_smiles=connectivity,
                       warnings=warnings, descriptors=descriptors, status="valid")
            counts["valid_rows"] += 1
            if fragment_count > 1:
                counts["multi_fragment_rows"] += 1
            if charge:
                counts["charged_rows"] += 1
            fragment_molecules = Chem.GetMolFrags(molecule, asMols=True)
            if any(Chem.GetFormalCharge(fragment) for fragment in fragment_molecules):
                counts["charged_component_rows"] += 1
            if used_ids[base_id] > 1:
                row["warnings"].append("Duplicate compound ID renamed for a unique source identity.")
            if canonical in canonical_first:
                row["status"] = "duplicate"
                row["duplicate_of_row"] = canonical_first[canonical]
                row["source_data"] = source
                counts["duplicate_rows"] += 1
            else:
                canonical_first[canonical] = row_number
                counts["unique_structures"] += 1
                if connectivity in connectivity_first:
                    row["stereo_variant_of_row"] = connectivity_first[connectivity]
                    row["warnings"].append(
                        f"Same connectivity as row {connectivity_first[connectivity]}, but a different isomeric SMILES; retained as a separate structure."
                    )
                    counts["stereo_variants"] += 1
                else:
                    connectivity_first[connectivity] = row_number
                metadata = {key: value for key, value in source.items() if key not in excluded}
                metadata["compound_id"] = compound_id
                result_records.append({
                    "name": name,
                    "components": [{"type": "small_molecule", "chain": "L", "smiles": canonical}],
                    "metadata": metadata,
                    "source_key": f"import-row:{batch_id}:{row_number}",
                    "source_row_number": row_number,
                })
        except Exception as exc:
            row["errors"].append(str(exc))
            candidate = _annotation_stripped_candidate(raw)
            if candidate and candidate != raw:
                try:
                    with rdBase.BlockLogs():
                        candidate_molecule = Chem.MolFromSmiles(candidate, sanitize=True)
                    candidate_valid = (
                        candidate_molecule is not None
                        and candidate_molecule.GetNumHeavyAtoms() > 0
                        and not any(atom.GetAtomicNum() == 0 for atom in candidate_molecule.GetAtoms())
                    )
                except Exception:
                    candidate_valid = False
                if candidate_valid:
                    row["annotation_stripped_candidate"] = candidate
                    row["candidate_status"] = "Review required; vendor annotation removed"
                    row["warnings"].append("A possible annotation-stripped candidate is available; it was not imported or substituted automatically.")
        if row["errors"]:
            row["status"] = "invalid"
            row["source_data"] = source
            counts["invalid_rows"] += 1
        report_rows.append(row)

    import pandas as pd
    weights = [row["descriptors"]["molecular_weight"] for row in report_rows if row["status"] in {"valid", "duplicate"}]
    extra_fields = [column for column in table.columns if column not in excluded]
    numeric_columns = []
    for field in extra_fields:
        values = table[field].astype(str).str.strip()
        nonempty_count = int(values.ne("").sum())
        numeric_values = pd.to_numeric(values.where(values.ne("")), errors="coerce").dropna()
        numeric_values = numeric_values[numeric_values.map(math.isfinite)]
        if numeric_values.empty:
            continue
        numeric_columns.append({
            "column": field,
            "numeric_count": int(len(numeric_values)),
            "nonempty_count": nonempty_count,
            "numeric_fraction": round(len(numeric_values) / max(nonempty_count, 1), 4),
            "missing_count": int(len(table) - nonempty_count),
            "minimum": float(numeric_values.min()),
            "median": float(numeric_values.median()),
            "mean": float(numeric_values.mean()),
            "maximum": float(numeric_values.max()),
        })
    counts.update(
        source_row_count=counts["rows"],
        valid_count=counts["valid_rows"],
        invalid_count=counts["invalid_rows"] + counts["empty_rows"],
        duplicate_structure_count=counts["duplicate_rows"],
        molecular_weight_min=min(weights) if weights else None,
        molecular_weight_mean=sum(weights) / len(weights) if weights else None,
        molecular_weight_max=max(weights) if weights else None,
        numeric_source_column_count=len(numeric_columns),
        numeric_columns=numeric_columns,
    )
    return {
        "batch_id": batch_id,
        "source_name": Path(filename).name,
        "source_sha256": hashlib.sha256(data).hexdigest(),
        "mapping": {"smiles_column": smiles_column, "name_column": name_column,
                    "id_column": id_column, "sheet_name": sheet_name},
        "summary": counts,
        "rows": report_rows,
        "records": result_records,
    }


def analyze_smiles_text(text: str) -> dict:
    """Row-wise SMILES validation for pasted text or SMI files."""
    import pandas as pd
    if len(str(text).encode("utf-8")) > MAX_BYTES:
        raise ValueError("The SMILES input exceeds 100 MB.")
    lines = []
    for raw_line in str(text).splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        values = line.split(maxsplit=1)
        lines.append({"SMILES": values[0], "Name": values[1] if len(values) > 1 else ""})
        if len(lines) > MAX_ROWS:
            raise ValueError("An import can contain at most 100,000 molecules.")
    if not lines:
        raise ValueError("Enter at least one SMILES string.")
    table_bytes = io.BytesIO()
    pd.DataFrame(lines).to_csv(table_bytes, index=False)
    return analyze_small_molecule_table(
        table_bytes.getvalue(), "pasted-smiles.csv", smiles_column="SMILES", name_column="Name"
    )


def parse_sdf(data: bytes, filename: str = "molecules.sdf") -> list[dict]:
    from rdkit import Chem
    if len(data) > MAX_BYTES:
        raise ValueError("The structure file exceeds 100 MB.")
    supplier = Chem.ForwardSDMolSupplier(io.BytesIO(data), removeHs=False)
    result = []
    for index, molecule in enumerate(supplier, 1):
        if molecule is None:
            raise ValueError(f"Could not read molecule {index} in {filename}.")
        name = molecule.GetProp("_Name").strip() if molecule.HasProp("_Name") else ""
        smiles = Chem.MolToSmiles(Chem.RemoveHs(molecule), isomericSmiles=True)
        block = Chem.MolToMolBlock(molecule) + "\n$$$$\n"
        result.append({"name": name or f"{Path(filename).stem} {index}",
                       "components": [canonical_component({"type": "small_molecule", "chain": "L", "smiles": smiles})],
                       "metadata": molecule.GetPropsAsDict(),
                       "files": [{"name": f"molecule-{index}.sdf", "content": block.encode(), "label": "Imported structure",
                                  "dimensions": 3 if molecule.GetNumConformers() and molecule.GetConformer().Is3D() else 2}]})
        if len(result) > MAX_ROWS:
            raise ValueError("An import can contain at most 100,000 molecules.")
    if not result:
        raise ValueError("No molecules were found in this file.")
    return result


def parse_structure(data: bytes, filename: str, *, include_ligands: bool = True, chains: list[str] | None = None) -> list[dict]:
    from Bio.PDB import PDBParser, MMCIFParser
    from Bio.Data.PDBData import protein_letters_3to1_extended
    if filename.lower().endswith(".gz"):
        with gzip.GzipFile(fileobj=io.BytesIO(data)) as handle:
            data = handle.read(MAX_BYTES + 1)
        filename = filename[:-3]
    if len(data) > MAX_BYTES:
        raise ValueError("The structure exceeds 100 MB.")
    suffix = Path(filename).suffix.lower()
    if suffix in {".sdf", ".mol"}:
        return parse_sdf(data, filename)
    if suffix in {".cif", ".mmcif"}:
        class _OptionalAtomSiteFieldsMMCIFParser(MMCIFParser):
            def _build_structure(self, structure_id):
                # Generated mmCIFs may omit these optional atom-site fields.
                # Occupancy's mmCIF default is 1.0; absent B factors are
                # unknown, so retain that as NaN rather than inventing a value.
                atom_count = len(self._mmcif_dict.get("_atom_site.id", []))
                self._mmcif_dict.setdefault("_atom_site.occupancy", ["1.0"] * atom_count)
                self._mmcif_dict.setdefault("_atom_site.B_iso_or_equiv", ["nan"] * atom_count)
                return super()._build_structure(structure_id)

        parser = _OptionalAtomSiteFieldsMMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    try:
        model = next(parser.get_structure("molecule", io.StringIO(data.decode("utf-8"))).get_models())
    except Exception as exc:
        raise ValueError(f"Could not read {filename}: {exc}") from exc
    components = []
    for chain in model:
        chain_id = chain.id.strip() or "_"
        if chains is not None and chain_id not in chains:
            continue
        protein, dna, rna = [], [], []
        for residue in chain:
            code = residue.resname.strip().upper()
            if code in {"HOH", "WAT", "DOD"}:
                continue
            if code in protein_letters_3to1_extended:
                protein.append(protein_letters_3to1_extended[code])
            elif code in {"DA", "DC", "DG", "DT", "DI", "DN"}:
                dna.append(code[-1].replace("I", "N"))
            elif code in {"A", "C", "G", "U", "I", "N"}:
                rna.append(code.replace("I", "N"))
            elif include_ligands and len(list(residue.get_atoms())) >= 1:
                locator = f"{chain_id}:{residue.id[1]}{residue.id[2].strip()}"
                components.append({"type": "small_molecule", "chain": f"{chain_id}_{residue.id[1]}{residue.id[2].strip()}",
                                   "ccd": code, "residue_key": locator})
        for kind, letters in (("protein", protein), ("dna", dna), ("rna", rna)):
            if letters:
                identifier = chain_id if not any(row["chain"] == chain_id for row in components) else f"{chain_id}_{kind}"
                components.append({"type": kind, "chain": identifier, "sequence": "".join(letters)})
    if not components:
        raise ValueError("No supported molecular chains were found in the structure.")
    components, _, _ = canonical_components(components)
    return [{"name": Path(filename).stem, "components": components, "metadata": {},
             "files": [{"name": Path(filename).name, "content": data, "label": "Imported structure", "dimensions": 3}]}]


def parse_archive(data: bytes) -> list[dict]:
    result, total = [], 0
    try:
        archive_file = zipfile.ZipFile(io.BytesIO(data))
    except zipfile.BadZipFile as exc:
        raise ValueError("This is not a readable ZIP archive.") from exc
    with archive_file as archive:
        entries = [entry for entry in archive.infolist() if not entry.is_dir() and
                   entry.filename.lower().endswith((".pdb", ".ent", ".cif", ".mmcif", ".sdf", ".mol"))]
        if not entries or len(entries) > 500:
            raise ValueError("A ZIP must contain between 1 and 500 structure files.")
        for entry in entries:
            total += entry.file_size
            if total > MAX_BYTES:
                raise ValueError("The uncompressed archive exceeds 100 MB.")
            try:
                result.extend(parse_structure(archive.read(entry), Path(entry.filename).name))
            except (zipfile.BadZipFile, RuntimeError) as exc:
                raise ValueError("A structure in the archive could not be read.") from exc
    return result


def fetch_identifier(identifier: str, source: str) -> tuple[str, bytes]:
    identifier = identifier.strip().upper()
    if source == "PDB" and re.fullmatch(r"[0-9][A-Z0-9]{3}", identifier):
        filename, url = f"{identifier.lower()}.cif", f"https://files.rcsb.org/download/{identifier}.cif"
    elif source == "UniProt" and re.fullmatch(r"[A-Z0-9]{6,10}(?:-[0-9]+)?", identifier):
        filename, url = f"{identifier}.fasta", f"https://rest.uniprot.org/uniprotkb/{identifier}.fasta"
    else:
        raise ValueError("Enter a valid PDB accession or UniProt accession for the selected source.")
    with urlopen(Request(url, headers={"User-Agent": "MolecularDesignStudio/1.0"}), timeout=30) as response:
        data = response.read(MAX_BYTES + 1)
    if len(data) > MAX_BYTES:
        raise ValueError("The fetched file exceeds 100 MB.")
    return filename, data


def csv_bytes(rows: list[dict]) -> bytes:
    columns = list(dict.fromkeys(key for row in rows for key in row))
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=columns)
    writer.writeheader()
    for row in rows:
        writer.writerow({key: json.dumps(value, ensure_ascii=False) if isinstance(value, (list, dict)) else value for key, value in row.items()})
    return buffer.getvalue().encode("utf-8")
