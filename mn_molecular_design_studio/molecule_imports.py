"""Validated molecular imports shared by the Database UI and API."""
from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
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


def read_table(data: bytes, filename: str):
    import pandas as pd
    if len(data) > MAX_BYTES:
        raise ValueError("The table exceeds 100 MB.")
    suffix = Path(filename).suffix.lower()
    if suffix in {".xlsx", ".xls"}:
        table = pd.read_excel(io.BytesIO(data), keep_default_na=False)
    else:
        table = pd.read_csv(io.BytesIO(data), sep="\t" if suffix in {".tsv", ".txt"} else ",", keep_default_na=False)
    table.columns = [str(value) for value in table.columns]
    if len(set(table.columns)) != len(table.columns):
        raise ValueError("Table columns must have unique names.")
    if len(table) > MAX_ROWS:
        raise ValueError("An import can contain at most 100,000 molecules.")
    return table


def parse_table(data: bytes, filename: str, *, chain_columns: dict[str, dict], name_column: str = "") -> list[dict]:
    table = read_table(data, filename)
    if not chain_columns:
        raise ValueError("Map at least one chain column.")
    if any(column not in table.columns for column in chain_columns) or (name_column and name_column not in table.columns):
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
        result.append({"name": str(source.get(name_column) or f"Molecule {index}"), "components": components,
                       "metadata": {key: value for key, value in source.items() if key not in chain_columns and key != name_column}})
    if not result:
        raise ValueError("The table contains no molecules.")
    return result


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
    parser = MMCIFParser(QUIET=True) if suffix in {".cif", ".mmcif"} else PDBParser(QUIET=True)
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
            elif include_ligands and len(list(residue.get_atoms())) > 1:
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
