"""Molecular input preparation shared by HTTP and Streamlit clients."""
from __future__ import annotations

import io
from pathlib import Path


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
                return int((chain in protein_chains and not residue.id[0].strip()) or locator in ligands)

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
    return {"name": record["name"] + " · " + "+".join(chains), "components": components,
            "metadata": record.get("metadata") or {}, "files": files}
