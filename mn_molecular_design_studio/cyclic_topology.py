"""Read recorded peptide closure bonds from coordinate structures."""
from pathlib import Path


def head_to_tail_chains(path: Path) -> set[str]:
    """Return protein chains with an explicit covalent first-N to last-C bond."""
    import gemmi

    try:
        block = gemmi.cif.read(str(path)).sole_block()
        atoms = block.find_loop("_atom_site.id")
        if atoms:
            # AtomWorks annotated CIFs may omit these two authoring columns.
            for column, default in (("_atom_site.occupancy", "1.0"), ("_atom_site.B_iso_or_equiv", "0.0")):
                if not block.find_values(column):
                    atoms.get_loop().add_columns([column], default)
        structure = gemmi.make_structure_from_block(block)
    except (RuntimeError, ValueError) as exc:
        raise ValueError(f"Could not read the structure {path.name}.") from exc
    if not len(structure):
        raise ValueError(f"{path.name} has no coordinate model.")

    terminals = {}
    for chain in structure[0]:
        residues = [res for res in chain if gemmi.find_tabulated_residue(res.name).is_amino_acid()]
        if len(residues) >= 3:
            terminals[chain.name] = {(str(residues[0].seqid), "N"), (str(residues[-1].seqid), "C")}

    closed = set()
    for connection in structure.connections:
        if connection.type != gemmi.ConnectionType.Covale:
            continue
        left, right = connection.partner1, connection.partner2
        if left.chain_name == right.chain_name and terminals.get(left.chain_name) == {
            (str(left.res_id.seqid), left.atom_name), (str(right.res_id.seqid), right.atom_name)
        }:
            closed.add(left.chain_name)
    return closed
