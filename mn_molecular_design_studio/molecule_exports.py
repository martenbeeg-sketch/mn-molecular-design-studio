"""Exports preserve exact selected membership and available molecular data."""
from __future__ import annotations
import io
import json
import zipfile

from mn_molecular_design_studio.molecule_imports import csv_bytes


def export_molecules(database, user, ids, format):
    database.molecules_by_ids(user, ids)
    items = database.list_molecules(user, molecule_ids=ids, limit=100_000)["items"]
    by_id = {row["id"]: row for row in items}
    records = [by_id[identifier] for identifier in ids]
    if format == "CSV":
        rows = []
        for record in records:
            row = {"name": record["name"], "molecule_id": record["id"], "type": record["type"], **record["metadata"], **record["scores"]}
            for component in record["components"]:
                row[component["chain"]] = component.get("sequence") or component.get("smiles") or "CCD:" + component["ccd"]
            rows.append(row)
        return csv_bytes(rows), "csv", "text/csv"
    if format == "FASTA":
        rows = []
        for record in records:
            if any(not row.get("sequence") for row in record["components"]):
                raise ValueError("FASTA export requires molecules containing only sequences. Use CSV for mixed complexes.")
            rows.extend([">" + record["name"].replace("\n", " "), ":".join(row["sequence"] for row in record["components"])])
        return ("\n".join(rows) + "\n").encode(), "fasta", "text/plain"
    if format in {"SMILES", "SDF"}:
        if any(len(record["components"]) != 1 or not record["components"][0].get("smiles") for record in records):
            raise ValueError("This export requires single compounds with known SMILES.")
        if format == "SMILES":
            return ("\n".join(record["components"][0]["smiles"] + " " + record["name"].replace("\n", " ") for record in records) + "\n").encode(), "smi", "chemical/x-daylight-smiles"
        from rdkit import Chem
        from rdkit.Chem import AllChem
        buffer = io.StringIO()
        writer = Chem.SDWriter(buffer)
        for record in records:
            structure = next((row for row in record["structures"] if row["format"] in {"sdf", "mol"}), None)
            molecule = None
            if structure:
                path = database.structure_file(user, structure["id"])
                molecule = Chem.SDMolSupplier(str(path), removeHs=False)[int(structure["context"].get("molecule_index", 0))]
            if molecule is None:
                molecule = Chem.MolFromSmiles(record["components"][0]["smiles"])
                AllChem.Compute2DCoords(molecule)
            molecule.SetProp("_Name", record["name"])
            molecule.SetProp("studio_molecule_id", record["id"])
            for key, value in {**record["metadata"], **record["scores"]}.items():
                molecule.SetProp(str(key), json.dumps(value) if isinstance(value, (dict, list)) else str(value))
            writer.write(molecule)
        writer.flush()
        value = buffer.getvalue().encode()
        writer.close()
        return value, "sdf", "chemical/x-mdl-sdfile"
    if format == "Structures ZIP":
        buffer = io.BytesIO()
        with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            for record in records:
                if not record["structures"]:
                    raise ValueError(f"{record['name']} has no structure to export.")
                representation = record["structures"][0]
                path = database.structure_file(user, representation["id"])
                archive.writestr(record["id"] + path.suffix, path.read_bytes())
        return buffer.getvalue(), "zip", "application/zip"
    raise ValueError("Choose a supported export format.")
