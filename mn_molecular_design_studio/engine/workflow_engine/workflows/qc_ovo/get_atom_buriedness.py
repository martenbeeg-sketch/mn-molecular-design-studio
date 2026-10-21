import pandas as pd
from Bio.Data.PDBData import protein_letters_3to1
from Bio.PDB import Structure
from scipy.spatial import cKDTree


def get_sequence(structure: Structure) -> str:
    """Get sequence of heavy (and light chain) in one string.

    Assumes only aa residues are in the structure and it contains coords for all the residues,
    e.g. predicted structures.

    First chain is asserted to be the heavy."""
    # FIXME
    chains = list(structure.get_chains())
    # assert 1 <= len(chains) <= 2
    # assert chains[0].id == 'H'  # first chain is Heavy

    # X for ambiguous/unknown (as in pep-patch)
    return "".join([protein_letters_3to1.get(r.resname, "X") for r in structure.get_residues()])


def get_atom_buriedness(structure, radii: list[float]):
    """Calculates the number of atoms within a specified radius for each atom in a PDB file.

    Parameters:
    structure (Structure): Biopython structure.
    radii (List[float]): Radii to calculate buriedness.

    Returns DataFrame with atom number, residue_id (str, possibly with insertion code),
        chain_id, and buriedness (atom count in query atom's vicinity) for each radius.
    """
    # Extract atom coordinates and info
    atoms = [
        (atom, atom.get_coord())
        for chain in structure.get_chains()
        for residue in chain.get_residues()
        for atom in residue.get_atoms()
    ]

    # Create a KDTree
    kdtree = cKDTree([coord for _, coord in atoms])

    # Prepare a dictionary to hold the results
    results = {}

    # Calculate buriedness for each atom
    for atom, coord in atoms:
        _, res_number, insertion_code = atom.get_parent().get_id()
        res_id = str(res_number) if insertion_code in (" ", "") else str(res_number) + insertion_code  # '80' or '80A'

        atom_key = (atom.get_serial_number(), res_id, atom.get_parent().get_parent().get_id())
        results[atom_key] = {}
        for radius in radii:
            # Count atoms within the radius
            count = kdtree.query_ball_point(coord, radius, return_length=True) - 1  # exclude the atom itself
            results[atom_key][f"buriedness_{radius}"] = count

    # Create a DataFrame from the results
    df = pd.DataFrame.from_dict(results, orient="index").reset_index()
    df.columns = ["atom_number", "residue_id", "chain_id"] + [f"buriedness_{r}" for r in radii]

    # Create a residue_serial column that increments for each unique residue_id
    df = df.sort_values(
        "atom_number"
    ).reset_index()  # first sort to get the same order as in PDB file (if biopython atom order is not guaranteed)
    df["residue_serial"] = df.groupby(
        ["residue_id", "chain_id"], sort=False
    ).ngroup()  # group by also chain (sort=False, so should be safe without chain gb even)

    return df
