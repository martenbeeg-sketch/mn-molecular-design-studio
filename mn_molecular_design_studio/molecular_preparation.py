"""Molecular input preparation shared by HTTP and Streamlit clients."""
from __future__ import annotations

import io
from pathlib import Path

from Bio.PDB.Polypeptide import is_aa, is_nucleic


def _polymer_residue(residue):
    return is_aa(residue) or is_nucleic(residue)


def coordinate_chain_ranges(data, filename):
    """Return author residue bounds for polymer chains, plus non-polymer-only chains."""
    from mn_molecular_design_studio.molecule_imports import read_coordinate_structure

    model = next(read_coordinate_structure(data, filename).get_models())
    ranges = {}
    for chain in model:
        if len(chain):
            positions = [residue.id[1] for residue in chain if _polymer_residue(residue)]
            ranges[chain.id.strip() or "_"] = (min(positions), max(positions)) if positions else None
    return ranges


def cropped_record(record, chains):
    from Bio.PDB import PDBIO, PDBParser, MMCIFParser, Select
    components = [row for row in record["components"] if row["chain"] in chains]
    if not components:
        raise ValueError("Select at least one chain.")
    if set(chains) - {row["chain"] for row in record["components"]}:
        raise ValueError("A selected chain is not part of this molecule.")
    files = []
    if record.get("files"):
        source = record["files"][0]
        data = source.get("content") or Path(source["path"]).read_bytes()
        filename = source.get("name") or source["path"]
        parser = MMCIFParser(QUIET=True) if str(filename).endswith((".cif", ".mmcif")) else PDBParser(QUIET=True)
        structure = parser.get_structure("crop", io.StringIO(data.decode()))
        protein_chains = {row["chain"] for row in components if row["type"] != "small_molecule"}
        ligands = {row.get("residue_key") for row in components if row["type"] == "small_molecule"}

        class Selection(Select):
            def accept_model(self, model):
                return model.id == 0

            def accept_residue(self, residue):
                chain = residue.parent.id.strip() or "_"
                locator = f"{chain}:{residue.id[1]}{residue.id[2].strip()}"
                return int((chain in protein_chains and _polymer_residue(residue)) or locator in ligands)

        writer, buffer = PDBIO(), io.StringIO()
        writer.set_structure(structure)
        try:
            writer.save(buffer, Selection())
        except ValueError:
            from Bio.PDB import MMCIFIO
            writer, buffer = MMCIFIO(), io.StringIO()
            writer.set_structure(structure)
            writer.save(buffer, Selection())
            files.append({"name": "cropped.cif", "content": buffer.getvalue().encode(), "label": "Cropped structure", "dimensions": 3})
        else:
            files.append({"name": "cropped.pdb", "content": buffer.getvalue().encode(), "label": "Cropped structure", "dimensions": 3})
        from mn_molecular_design_studio.preparation_history import extend_history, sha256

        files[0]["preparation_history"] = extend_history(
            source.get("preparation_history"), sha256(data), sha256(files[0]["content"]),
            [{"label": "Select chains", "details": {"chains": list(chains)}}],
        )
    return {"name": record["name"] + " · " + "+".join(chains), "components": components,
            "metadata": record.get("metadata") or {}, "files": files}


def prepare_coordinate_subset(data, filename, ranges, *, remove_waters=True, remove_nonpolymer=False):
    """Trim polymer ranges and apply water and ligand choices independently."""
    from Bio.PDB import MMCIFIO, PDBIO, Select

    from mn_molecular_design_studio.molecule_imports import read_coordinate_structure

    cif = Path(filename).suffix.lower() in {".cif", ".mmcif"}
    structure = read_coordinate_structure(data, filename)
    if not ranges or any(bounds is not None and bounds[0] > bounds[1] for bounds in ranges.values()):
        raise ValueError("Select at least one chain with a valid residue range.")

    class Subset(Select):
        def accept_model(self, model):
            return model.id == 0

        def accept_residue(self, residue):
            chain = residue.parent.id.strip() or "_"
            if chain not in ranges:
                return False
            if residue.resname.strip().upper() in {"HOH", "WAT", "DOD"}:
                return not remove_waters
            if _polymer_residue(residue):
                bounds = ranges[chain]
                return bounds is not None and bounds[0] <= residue.id[1] <= bounds[1]
            return not remove_nonpolymer

    writer = MMCIFIO() if cif else PDBIO()
    writer.set_structure(structure)
    buffer = io.StringIO()
    writer.save(buffer, Subset())
    return buffer.getvalue().encode("utf-8")
