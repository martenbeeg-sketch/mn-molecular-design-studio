import numpy as np
import pandas as pd
from surface_analyses.patches import assign_patches


# def get_cdrs(traj: Trajectory, cdr_def_scheme: str = 'chothia'):  # remove dependency on traj from pep-patch
def get_cdrs(sequence: str, cdr_def_scheme: str = "chothia"):
    """Adapted from Pep patch.

    Added to the hydrophobic calculation. Fixed the atom's residues to include chain id to not be ambiguous across chains.
    Returns also the number of CDR loop, not just binary CDR flag.
    Todo maybe create public PR for surface_analyses for this and multiple cutoffs?"""
    try:
        from cdr_annotation import CustomAnnotation

        return {
            # str(traj.top.residue(i)): loop_no
            res_serial: loop_no
            # for i, loop_no in CustomAnnotation.from_traj(traj, scheme=cdr_def_scheme).cdr_indices()
            for res_serial, loop_no in CustomAnnotation(sequence, scheme=cdr_def_scheme).cdr_indices()
        }
    except ImportError as e:
        print(
            f"CDR annotation failed with the following error:\n{e}\n"
            "If the error pertains to the annotation tool ANARCI or ANARCI is missing, "
            "a fresh installation of ANARCI ( https://github.com/oxpig/ANARCI ) or its dependencies might help.\n\n"
            "To use pep_patch_hydrophobic without CDR annotation, rerun the script without the '--check_cdrs' flag."
        )
        raise RuntimeError("CDR Annotation failed")


# def get_vertex_df(surf, pdb_traj: Trajectory, annotate_cdrs=True):  # remove dependency on traj from pep-patch
def get_vertex_df(surf, atom_to_residue: pd.DataFrame, sequence: str, annotate_cdrs=True):
    """Annotate vertices with CDR loop number, compute vertex area."""
    # cdrs = get_cdrs(pdb_traj) if annotate_cdrs else {}  # remove dependency on pdb_traj from pep-patch
    cdrs = get_cdrs(sequence) if annotate_cdrs else {}

    # residues = np.array([str(a.residue) for a in pdb_traj.top.atoms])
    # this is the residue serial (same as in pep-patch), first are H residues and then L residue ids start where the H left off
    residues = atom_to_residue["residue_serial"].to_numpy()

    # Calculate the area of each triangle, and split evenly among the vertices
    tri_areas = surf.areas()
    vert_areas = np.zeros(surf.vertices.shape[0])
    # Distribute the triangle areas evenly among their vertices
    np.add.at(
        vert_areas, surf.faces, tri_areas[:, np.newaxis] / 3.0
    )  # add newaxis, so tri_areas gets broadcasted to surf.faces (n,3)
    return pd.DataFrame(
        {
            "area": vert_areas,
            "atom": surf["atom"],
            "cdr": np.vectorize(cdrs.get)(residues, -1)[surf["atom"]],  # assign CDR loop no
            "value": surf["values"],
        }
    )


def assign_patches_electrostatic(vertex_potential: pd.Series, faces: np.array, patch_cutoff: tuple[float]):
    """Aggregate surface in patches using the cutoff.

    Returns a pd.DataFrame of vertices like in PEP Patch.
    Additionally, patch_id column is added, which has unique identified for each patch (unique over the union of positive and negative).
    - Vertices with *positive* charge patches have positive patch_id (1..inf)
    - Vertices with *negative* charge patches have negative patch_id (-1..-inf)
    - Vertices with no patches have patch_id set to 0

    Adapted from pep patch so we can run this with multiple patch_cutoffs while pep patch (APBS) can run just once.

    [done] And it would be great to remove the dependency with pdb_traj.., we would just need the sequence of residues for each atom and total sequence H+L, which both we have in atom_to_residue
    """
    # code from PEP-Patch
    positive_patches = assign_patches(faces, vertex_potential > patch_cutoff[0])
    negative_patches = assign_patches(faces, vertex_potential < patch_cutoff[1])

    # create column patch_id so that patch_id will be unique within positive _and_ negative patches
    patch_id = positive_patches + 1
    # renumber negative patches, so we can discern between positive and negative
    patch_id[negative_patches != -1] = -1 - negative_patches[negative_patches != -1]

    return patch_id
