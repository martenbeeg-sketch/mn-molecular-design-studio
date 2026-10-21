import json
from argparse import Namespace
from functools import partial
from pathlib import Path
from typing import NotRequired, Optional, TypedDict, Union

import pandas as pd
import surface_analyses.commandline_electrostatic
import surface_analyses.commandline_hydrophobic
from Bio.PDB import PDBParser, Structure
from feature_generation import (
    compute_surface_integrals,
    get_feat_name,
    get_global_electrostatic_descriptors_for_patches,
    get_global_hydrophobic_descriptors,
    get_global_hydrophobic_descriptors_for_patches,
    get_vertices_sub_dfs,
)
from get_atom_buriedness import get_atom_buriedness
from peppatch_postprocess import assign_patches_electrostatic, get_vertex_df
from surface_analyses.patches import assign_patches
from surface_analyses.structure import load_trajectory_using_commandline_args
from surface_analyses.surface import Surface


class ElectrostaticPatchConfig(TypedDict):
    patch_cutoff: tuple[float]  # (2., -2.)  default setting in pep patch
    name: str  # identifier/description
    aggregate_per_residue: NotRequired[bool]


class HydrophobicPatchConfig(TypedDict):
    patch_min: float  # threshold value (potential) for inclusion of a vertex in a patch, default of pep patch: 0.12 (works well with WiWh)    name: str
    name: str
    aggregate_per_residue: NotRequired[bool]


class ElectrostaticConfig(TypedDict):
    pH: float  # e.g. 7.4 of blood
    volume_integral_cutoff: NotRequired[tuple[float]]  # (0.3, -0.3)  default setting in pep patch
    patch_configs: list[ElectrostaticPatchConfig]

    # just for computing ML features:
    surface_integral_high_cutoffs: list[float]
    surface_integral_low_cutoffs: list[float]


class HydrophobicConfig(TypedDict):
    hydrophobicity_scale: str  # path to scale, ex: data/hydrophobicity_scales/scale_WiWh.csv
    patch_configs: list[HydrophobicPatchConfig]

    # just for computing ML features:
    surface_integral_cutoffs: list[float]


class PepPatchMethod(TypedDict):
    name: str  # identifier/description e.g. electrostatic_ph7.4
    method: str  # electrostatic/hydrophobic
    config: Union[ElectrostaticConfig, HydrophobicConfig]  # pH, hydrophobicity scale


class PepPatchSet(TypedDict):
    buriedness_radii: list[float]
    methods: list[PepPatchMethod]


class EventType(TypedDict):
    sample_id: str
    structure_prediction_method: str
    structure_prediction_scheme: str

    set_name: str
    set: NotRequired[
        PepPatchSet
    ]  # if not specified, hardcoded settings from codebase will be selected based on 'set_name'


def get_atom_to_residue_mapping(structure: Structure, atom_buriedness: pd.DataFrame = None):
    """Return DataFrame with atom index (0-based contiguous) -> chain_id, residue_id, residue_serial, atom_number (in pdb records)

    residue_id is '80' or '80A', residue_serial is zero-based and starts with H and then L chain residues (if not VHH).

    Residue serial is helpful when computing the CDRS using extended pep-patch code.
    """
    if atom_buriedness is None:
        atom_buriedness = get_atom_buriedness(structure, [0.0])  # 0.0 just as a placeholder

    return atom_buriedness[["residue_id", "chain_id", "residue_serial", "atom_number"]]


def patch_membership_agg(sub_df, patch_col):
    """Get the area of each patch on a residue (or a groupby). Performs an aggregation within a residue group (meant
      to by `apply`ed to a DataFrameGroupBy) to sum vertex 'area' for each value of patch_id (patch_col).
    Returns a dictionary with patch_id: respective area.  ex: {0: 0.7975, 3: 0.7426}
    """
    return sub_df.groupby(patch_col)["area"].sum().to_dict()


def agg_patches_per_residue(patches_vertices):
    """Convert patch vertices into residue-level patches.

    Returns a DataFrame, ex:
                              area      value                                       patch_labels  patch_label
    chain_id residue_id
    H        1           23.310820   2.124569                            {1: 23.310819893187727}            1
             13           9.196304   2.235697                             {19: 9.19630379374818}           19

    Returns a DataFrame with residue SASA, mean vertex value, and patches per residue and the patch with max area on
    the residue.
    """

    patches_gb_residue = patches_vertices.groupby(["chain_id", "residue_id"])

    residue_patches = patches_gb_residue.agg(
        {
            "area": "sum",
            "value": "mean",
            # todo atom buriedness?? of the patch? for position-ml feats
        }
    )

    # all patches the residue belongs in + respective area on the residue (dict path_id: area)
    residue_patches["patch_labels"] = patches_gb_residue.apply(
        partial(patch_membership_agg, patch_col="patch_id"), include_groups=False
    )
    # patch_id of the largest patch
    residue_patches["patch_label"] = residue_patches["patch_labels"].apply(lambda x: max(x, key=x.get))
    return residue_patches


def run_electrostatic(pdb_file: str, config: ElectrostaticConfig, atom_to_residue: pd.DataFrame, sequence: str):
    """Run PEP-Patch electrostatic and return surface vertices annotated with patch labels for each patch cutoff.

    Returns PEP-Patch result an an iterator of tuples patch_config, and DataFrame patch_vertices.

              area  atom   cdr     value  patch_id  residue_id    chain_id
    0     0.144793   410  True -0.226922         0          29           H
    1     0.084391   410  True -0.224813         0          29           H
    2     0.100023   410  True -0.211886         0          29           H

    PEP Patch electrostatic runtime on mac: 10.5s, then 0.15s/(one patch cutoff) takes the agg in patches and per-residue

    For patch_vertices see ` surf_to_patch_vertices_electrostatic`. `patch_id` is 0 if a patch is not assigned,
        negative integer for negative patch and positive int for positive patch
    """
    # pH (not patch_cutoff) determines the output surface potential of PEP-Patch
    mock_trajectory_args = Namespace(parm=pdb_file, trajs=[pdb_file], stride=1, ref=None, protein_ref=None)
    traj = load_trajectory_using_commandline_args(mock_trajectory_args)

    pep_result = surface_analyses.commandline_electrostatic.run_electrostatics(
        traj,
        apbs_dir="apbs",
        pH=config["pH"],
        out="/dev/null",  # we are using the return value in `result`
    )

    # annotated vertices (cdrs, potential value, closest atom number, ...)
    vertex_data = get_vertex_df(pep_result["surface"], atom_to_residue, sequence, annotate_cdrs=False)
    # add residue id to each vertex (based on atom number)
    vertex_data = vertex_data.merge(atom_to_residue, left_on="atom", right_index=True)
    faces = pep_result["surface"].faces
    # todo could save (vertex_data, faces); just from this, patches and features are computed, maybe also `sequence` (or seqH + seqL)

    # return patches for every patch_cutoff
    def patch_vertices_iterator():
        for patch_config in config["patch_configs"]:
            patch_vertices = vertex_data.copy()
            patch_vertices["patch_id"] = assign_patches_electrostatic(
                vertex_data["value"], faces, patch_config["patch_cutoff"]
            )
            yield patch_config, patch_vertices

    return pep_result, patch_vertices_iterator()


def run_hydrophobic(pdb_file: str, config: HydrophobicConfig, atom_to_residue: pd.DataFrame, sequence: str):
    """Run PEP-Patch hydrophobic and return surface vertices annotated with patch labels for each patch cutoff.

    Returns an iterator of tuples patch_cutoff, and DataFrame patch_vertices.

              area  atom   cdr     value  patch_id  residue_id    chain_id
    0     0.144793   410  True -0.226922         0          29           H
    1     0.084391   410  True -0.224813         0          29           H
    2     0.100023   410  True -0.211886         0          29           H

    PEP hydrophobic runtime on mac: 3s, then 0.25s/(one patch cutoff) takes the agg in patches and per-residue

    Residues with hydrophobic patches have positive patch label (1..inf)
    Residues not belonging to a patch have a label 0 (if the residue vertices are mostly unassigned to a patch).
    """
    # hydrophobicity_scale determines the output surface potential of PEP-Patch (not patch_min's)
    mock_trajectory_args = Namespace(parm=pdb_file, trajs=[pdb_file], stride=1, ref=None, protein_ref=None)
    traj = load_trajectory_using_commandline_args(mock_trajectory_args)
    pep_result = surface_analyses.commandline_hydrophobic.run_hydrophobic(
        pdb_file,
        traj,
        scale=str(config["hydrophobicity_scale"]),
        potential=True,
        sap=True,
        surfscore=True,
        sh=True,  # additional per-atom features computed
    )

    surf = Surface.from_dict(pep_result, "hydrophobic_potential:0")
    # we could save surf

    # annotated vertices (cdrs, potential value, closest atom number, ...)
    vertex_data = get_vertex_df(surf, atom_to_residue, sequence, annotate_cdrs=False)
    vertex_data["area"] *= (
        100  # for some reason, pep_patch_hydrophobic computes with nanometers, but pep_patch_electrostatic with angstroms, convert nanometers^2 to angstroms^2
    )
    # add residue id to each vertex (based on atom number)
    vertex_data = vertex_data.merge(atom_to_residue, left_on="atom", right_index=True)
    faces = surf.faces
    # todo could save (vertex_data, faces); just from this, patches and features are computed, maybe also `sequence` (or seqH + seqL)

    # return patches for every patch_min
    def patch_vertices_iterator():
        for patch_config in config["patch_configs"]:
            patch_vertices = vertex_data.copy()
            patch_vertices["patch_id"] = assign_patches(faces, vertex_data["value"] > patch_config["patch_min"])
            patch_vertices["patch_id"] += (
                1  # 0 is not assigned, 1..inf are patches (for consistency with electrostatic)
            )
            yield patch_config, patch_vertices

    return pep_result, patch_vertices_iterator()


def process_pdb(input_pdb: str, peppatch_set: PepPatchSet):
    pdb_structure = PDBParser().get_structure("structure", input_pdb)

    sequence: str = ""  # Not used (logic for antibody CDR annotation)
    atom_buriedness = (
        get_atom_buriedness(pdb_structure, peppatch_set["buriedness_radii"])
        if peppatch_set["buriedness_radii"]
        else None
    )

    # get mapping atom_index -> chain_id, residue_id, so we can map atom numbers
    #   returned from PEP-Patch to residues on a chain
    atom_to_residue = get_atom_to_residue_mapping(pdb_structure, atom_buriedness)

    mab_type = ""  # Not used (logic for antibody CDR annotation)

    # function to aggregate patches per residue for visualization
    per_residue_patches = []

    def maybe_agg_per_residue(method, patch_config, patch_vertices):
        if "aggregate_per_residue" not in patch_config or not patch_config["aggregate_per_residue"]:
            return
        per_residue_patches.append(
            dict(
                **method,  # to identify the parameters the result below was computed with
                patch_config=patch_config,
                per_residue_patches=agg_patches_per_residue(patch_vertices).reset_index().to_dict(orient="records"),
            )
        )

    global_descriptors = {}
    for method in peppatch_set["methods"]:
        method_name = method["name"]

        # assert there are non-empty patch configs (patch_min)
        assert method["config"]["patch_configs"]

        def add_descriptors(
            descriptors: dict[str, Union[float, int]],
            patch_cutoff: Optional[Union[float, tuple[float]]] = None,
            sub_df_prefix: str = "",
        ):
            global_descriptors.update(
                {
                    get_feat_name(method_name, k, patch_cutoff, sub_df_prefix, mab_type): v
                    for k, v in descriptors.items()
                }
            )

        if method["method"] == "electrostatic":
            pep_result, patch_vertices_it = run_electrostatic(input_pdb, method["config"], atom_to_residue, sequence)

            # volume potential integrals computed by pep_patch_electrostatic
            add_descriptors(pep_result["integrals"])

            for patch_config, patch_vertices in patch_vertices_it:
                # optional per-residue aggregation for visualization in biophi
                maybe_agg_per_residue(method, patch_config, patch_vertices)
                # patch-based global descriptors
                add_descriptors(
                    get_global_electrostatic_descriptors_for_patches(patch_vertices, atom_buriedness),
                    patch_config["patch_cutoff"],
                )

            # surface integrals
            # (use last patch_vertices, as only patch definitions differ in the loop above, but a bit hacky)
            for sub_df_name, vertex_df in get_vertices_sub_dfs(patch_vertices):
                # this descriptor is in fact total patch area by the region (H/L cdr1 etc..), weighted by the intensity
                add_descriptors(
                    compute_surface_integrals(
                        vertex_df,
                        high_cutoffs=method["config"]["surface_integral_high_cutoffs"],
                        low_cutoffs=method["config"]["surface_integral_low_cutoffs"],
                    ),
                    sub_df_prefix=sub_df_name,
                )
        elif method["method"] == "hydrophobic":
            # elif method['method'] == 'hydrophobic':
            pep_result, patch_vertices_it = run_hydrophobic(input_pdb, method["config"], atom_to_residue, sequence)

            # additional descriptors computed by pep_patch_hydrophobic
            # e.g. surfscore is defined as saa * propensity of each atom, summed for the structure
            add_descriptors(get_global_hydrophobic_descriptors(pep_result))

            for patch_config, patch_vertices in patch_vertices_it:
                # optional per-residue aggregation for visualization in biophi
                maybe_agg_per_residue(method, patch_config, patch_vertices)

                # patch-based global descriptors
                add_descriptors(
                    get_global_hydrophobic_descriptors_for_patches(patch_vertices, atom_buriedness),
                    patch_config["patch_min"],
                )

            # surface integrals
            # (use last patch_vertices, as only patch definitions differ in the loop above, but a bit hacky)
            for sub_df_name, vertex_df in get_vertices_sub_dfs(patch_vertices):
                add_descriptors(
                    compute_surface_integrals(vertex_df, high_cutoffs=method["config"]["surface_integral_cutoffs"]),
                    sub_df_prefix=sub_df_name,
                )
        else:
            raise ValueError(f"Method {method['method']} not supported")

    if peppatch_set["methods"]:
        # reuse last variable - patch_vertices - for global definitions of cdr123 and chain total surface
        #   areas for the structure (only depends on the predicted structure and no other params)
        for sub_df_name, vertex_df in get_vertices_sub_dfs(patch_vertices):
            global_descriptors[f"{sub_df_name}area"] = vertex_df.area.sum()

    return per_residue_patches, global_descriptors


def get_peppatch_set_definitions(sets_path: str):
    """Our PEP-Patch settings for computing features/residue annotations.

    Can be stored here in this repo for versioning/consistency.
    """
    with Path(sets_path).open() as f:
        return json.load(f)


# CONFIG_STRUCTURE_BUCKET = os.getenv('CONFIG_STRUCTURE_BUCKET')
# CONFIG_DESCRIPTOR_BUCKET = os.getenv('CONFIG_DESCRIPTOR_BUCKET')


# def lambda_handler(event: EventType, context):
#     print('Received event:', event)
#     sample_id = event['sample_id']
#     set_name = event['set_name']
#     input_pdb_path_s3 = f'{sample_id}/{event["structure_prediction_method"]}_{event["structure_prediction_scheme"]}.pdb'
#     output_path_s3_base = f'{sample_id}/{event["structure_prediction_method"].lower()}_{set_name}'

#     try:
#         # get the `methods` from set definition, if there is one for the set_name, otherwise use event['methods']
#         set_definitions = get_peppatch_set_definitions()
#         if set_name in set_definitions:
#             assert 'methods' not in event, f'Set {set_name} methods are already defined in `data/feature_sets.json` but received another `methods` definitions in the Event.'
#             peppatch_set = set_definitions[set_name]
#         else:
#             peppatch_set = event['set']


#         # Use temp dir as working directory,
#         with tempfile.TemporaryDirectory() as tmpdirname, contextlib.chdir(tmpdirname):
#             input_pdb_filename = 'input.pdb'
#             s3_client.download_file(CONFIG_STRUCTURE_BUCKET, input_pdb_path_s3, input_pdb_filename)

#             per_residue_patches, global_descriptors = process_pdb(input_pdb_filename, peppatch_set)

#         s3_client.put_object(
#             Bucket=CONFIG_DESCRIPTOR_BUCKET,
#             Key=f'{output_path_s3_base}.csv',
#             Body=pd.DataFrame([global_descriptors]).to_csv(index=False).encode('utf-8'),
#         )

#         s3_client.put_object(
#             Bucket=CONFIG_DESCRIPTOR_BUCKET,
#             Key=f'{output_path_s3_base}_per_residue_patches.json',  # todo user cannot change numbering scheme after upload, no? If yes, include it here in the name
#             Body=json.dumps(per_residue_patches).encode('utf-8'),
#         )

#     except Exception as e:
#         traceback.print_exc()

#         s3_client.put_object(
#             Bucket=CONFIG_DESCRIPTOR_BUCKET,
#             Key=f'{output_path_s3_base}.err.json',
#             Body=json.dumps({'message': str(e), 'traceback': ''.join(traceback.format_exception(e))}).encode('utf-8'),
#         )
