import fnmatch
import warnings
from typing import Optional, Union

import numpy as np
import pandas as pd


def get_feat_name(
    method_name,
    feat_name,
    cutoff: Union[float, tuple[float]] = None,
    sub_df_prefix: str = "",
    mab_type: str = "",
):
    """Return feature name using the template in this function.

    example
    H_CDR3_surface_integral_high_0.4_fv|electrostatic_7.4    # with sub_df_prefix `H_CDR3_`
    top5_patches_mean_area_-_fv|electrostatic_7.4_c[3, -3]   # with patch cutoff_str `c[3, -3]`
    """
    cutoff_str = ""
    # if cutoff:
    #     cutoff_str = f"_c{cutoff}"

    if mab_type:
        mab_type = f"_{mab_type}"

    return f"{sub_df_prefix}{feat_name}{mab_type},{method_name}{cutoff_str}"


def get_vertices_sub_dfs(vertices: pd.DataFrame):
    """Produces vertex dfs with names:

    'Total_' (all),
    'Chain_{chain_id}_' (chain vertices) for each chain_id present in the structure.
    """
    vs = vertices
    unique_chains = vertices.chain_id.unique()

    dfs = {f"chain_{chain_id}_": vs[vs.chain_id == chain_id] for chain_id in unique_chains}
    dfs["total_"] = vs
    return dfs.items()


def expand_aggrs(filtered_df, col_pattern_aggrs: dict[str, list]):
    """Expanding the unix-style wildcards (using fnmatch) to actual column names."""

    col_aggrs = {}
    for col_pattern, aggrs in col_pattern_aggrs.items():
        # Use fnmatch to find matching columns
        matching_cols = [col for col in filtered_df.columns if fnmatch.fnmatch(col, col_pattern)]
        if not matching_cols:
            warnings.warn(
                f"No columns matched the colname pattern `{col_pattern}` - did you specify the pattern correctly?"
            )
        for col in matching_cols:
            col_aggrs[col] = aggrs
    return col_aggrs


def aggregate_patch_data(patches_df, top_ns=(1, 2, 3, 5)):
    """Create structure-level features (fixed length) from the dataframe of patches (variable length).

    Three types of features are returned:
    - over all patches (total path area, patch intensity etc)
    - over top n (1,2,...) patches
        - the values for n-largest patch (area, intensity,...)
        - values if 1..n-largest patch aggregated (equal weight for each patch disregarding its size)
            - e.g. for top3, the values for the largest three patches are meaned/summed, for five, the largest 5
                patches' values are meaned

    """

    feats = {
        "num_patches": len(patches_df),
    }

    col_aggrs = {
        "area": ["mean", "sum"],
        "intensity_sum": ["mean", "sum"],
        "intensity_mean": ["mean"],
        "mean_buriedness_*": ["mean"],
        "min_buriedness_*": [],
        "max_buriedness_*": [],
    }

    col_aggrs = expand_aggrs(patches_df, col_aggrs)

    fields = list(col_aggrs.keys())

    # struct-level aggregation, all patches
    for field in fields:
        for aggr in col_aggrs[field]:
            feats[f"{aggr}_patch_{field}"] = patches_df[field].aggregate(aggr)

    # 1) Values for top-n largest patch
    # 2) Aggregated values for all patches 1..top-n
    for n in top_ns:
        # if there are less than n patches - nans
        nlargest_patches = (
            patches_df.nlargest(n, "area")
            if len(patches_df) >= n
            else pd.DataFrame({c: [np.nan] * n for c in patches_df.columns})
        )

        top_n_patch = nlargest_patches.iloc[-1]

        for field in fields:
            feats[f"top{n}_patch_{field}"] = top_n_patch[field]

            if n == 1:  # no aggr. needed
                continue

            for aggr in col_aggrs[field]:
                top_n_feat_aggr = nlargest_patches[field].aggregate(aggr)
                feats[f"top{n}_patches_{aggr}_{field}"] = top_n_feat_aggr

    return feats


def compute_surface_integrals(
    vertices: pd.DataFrame, high_cutoffs: list = (), low_cutoffs: list = ()
) -> dict[str, float]:
    """Computes surface integrals of the potential over all vertices. sum_vertices(potential value * vertex area)

    Additionally, high and low integrals -> the summed quantity, value, is only taken as what is above/below a cutoff threshold.
    """
    features = {}

    features["surface_integral"] = np.sum(vertices["value"] * vertices["area"])

    for cutoff in high_cutoffs:
        features[f"surface_integral_high_{cutoff}"] = np.sum(
            np.maximum(vertices["value"] - cutoff, 0) * vertices["area"]
        )

    for cutoff in low_cutoffs:
        features[f"surface_integral_low_{cutoff}"] = np.sum(
            np.minimum(vertices["value"] - cutoff, 0) * vertices["area"]
        )

    return features


def get_buriedness_for_patches(vertices: pd.DataFrame, atom_buriedness: pd.DataFrame) -> pd.DataFrame:
    """Compute buriedness scores of the patches.

    Buriedness := #atoms in a set radius from query atom.
    Mean -> mean the buriedness of the patch's atoms
    Max, min -> take max/min of the patch's atoms, of their buriedness

    Returns a DataFrame with patch_id and buriedness columns.
    """
    # Select only vertices in patches
    vertices = vertices[vertices.patch_id != 0]
    # Merge with atom_buriedness on 'atom'
    vertices = vertices.drop(columns="chain_id").merge(
        atom_buriedness, left_on="atom", right_index=True
    )  # drop the col because atom_buriedness already has it

    by_patch = vertices.groupby("patch_id")

    # Calculate the total area of each patch
    total_area = by_patch["area"].sum()

    # Calculate the weighted mean, max, and min buriedness for each patch (for each buriedness radius)
    patch_buriedness_df = pd.DataFrame(index=by_patch.size().index)
    for radius_col in atom_buriedness.columns:
        if not radius_col.startswith("buriedness_"):
            continue

        weighted_sum_of_buriedness = (vertices[radius_col] * vertices["area"]).groupby(vertices["patch_id"]).sum()
        patch_buriedness_df[f"mean_{radius_col}"] = weighted_sum_of_buriedness / total_area
        patch_buriedness_df[f"max_{radius_col}"] = by_patch[radius_col].max()
        patch_buriedness_df[f"min_{radius_col}"] = by_patch[radius_col].min()

    return patch_buriedness_df


def get_patch_df(patch_vertices, atom_buriedness: Optional[pd.DataFrame]):
    """Aggregate values for patches' vertices into patch-level values.

    Add buriedness scores to each patch, if atom_buriedness is not None.
    Returns a DataFrame with patch_id and the aggregated values/buriedness score.
    """
    patch_df = patch_vertices.groupby("patch_id").aggregate(
        intensity_mean=pd.NamedAgg("value", "mean"),
        intensity_sum=pd.NamedAgg("value", "sum"),
        area=pd.NamedAgg("area", "sum"),
    )

    if atom_buriedness is not None:
        buriedness_df = get_buriedness_for_patches(patch_vertices, atom_buriedness)
        patch_df = patch_df.merge(buriedness_df, on="patch_id")

    return patch_df


def get_patch_area_features(patch_vertices):
    """Return dict of patch area broken down per all/selected chain/CDR1/2/3 vertices."""
    feats = {}
    # include total patch area per all/selected chain/CDR1/2/3
    for sub_df_name, sub_df in get_vertices_sub_dfs(patch_vertices):
        feats[f"{sub_df_name}patch_area"] = sub_df[sub_df["patch_id"] != 0]["area"].sum()
    return feats


def get_global_electrostatic_descriptors_for_patches(
    patch_vertices: pd.DataFrame, atom_buriedness: Optional[pd.DataFrame]
) -> dict[str, float]:
    """Compute scalar structure-level features from patch vertices.

    Returns a dict of scalar features {name: value}.
    - features with _+ in their name for positive patches, _- for negative

    Example feature names (when additionally passed through get_feat_name):
    - sum_patch_area_+|electrostatic_7.4_c[2, -2]
    - top5_patch_area_-|electrostatic_7.4_c[2, -2]       # this is value for the 5th largest patch
    - top5_patches_sum_area_-|electrostatic_7.4_c[2, -2] # this is value summed over _all_ 5 largest patches

    (Features can be nan – e.g. there are less then 5 patches but we're computing a value for top-5 patch.)
    """

    out = {}
    positive_df = patch_vertices[patch_vertices.patch_id > 0]
    negative_df = patch_vertices[patch_vertices.patch_id < 0]

    positive_patch_feats = aggregate_patch_data(get_patch_df(positive_df, atom_buriedness))
    negative_patch_feats = aggregate_patch_data(get_patch_df(negative_df, atom_buriedness))

    def with_suffix(dictionary, key_suffix):
        return {f"{k}{key_suffix}": v for k, v in dictionary.items()}

    out.update(with_suffix(positive_patch_feats, "_+"))
    out.update(with_suffix(negative_patch_feats, "_-"))
    out.update(with_suffix(get_patch_area_features(positive_df), "_+"))
    out.update(with_suffix(get_patch_area_features(negative_df), "_-"))
    return out


def get_global_hydrophobic_descriptors_for_patches(
    patch_vertices: pd.DataFrame, atom_buriedness: pd.DataFrame
) -> dict[str, float]:
    """Compute scalar structure-level features from patch vertices.

    Returns a dict of scalar features {name: value} like:
    - mean_patch_area|hydrophobic_Ja_c0.8
    - top2_patch_intensity_sum|hydrophobic_Ja_c0.8          # this is value for the 2nd largest patch
    - top2_patches_mean_intensity_sum|hydrophobic_Ja_c0.8   # this is value meaned over _both_ 2 largest patches

    (Names when additionally passed through get_feat_name)
    (Features can be nan – e.g. there are less then 5 patches but we're computing a value for top-5 patch.)
    """
    feats = aggregate_patch_data(get_patch_df(patch_vertices[patch_vertices.patch_id > 0], atom_buriedness))
    # include total patch area per all/selected chain/CDR1/2/3
    feats.update(get_patch_area_features(patch_vertices))
    return feats


def get_global_hydrophobic_descriptors(peppatch_hydrophobic_output: dict[str, np.ndarray]):
    """Aggregate per-atom descriptors computed by PEP-Patch hydrophobic; over the whole structure.

    (these outputs are enabled with extra flags to PEP-Patch)"""
    output = peppatch_hydrophobic_output

    # choose trajectory frame 0 because there may be multiple structures (PEP-Patch accepts mol. dynamics trajectory),
    # we have only a single structure
    FRAME = 0
    return dict(
        # sasa_score is similar to the potential integral, however this uses the raw hydrophobicity scale values
        # (but using which the potential is calculated :)
        sasa_score=np.sum(output["surfscore"][FRAME]),
        positive_sasa_score=np.sum(np.maximum(output["surfscore"][FRAME], 0)),
        sap_score=np.sum(output["sap"][FRAME]),
        surrounding_hydrophobicity=np.sum(output["surrounding_hydrophobicity"][FRAME]),
        sasa=np.sum(output["saa"][FRAME]),
    )
