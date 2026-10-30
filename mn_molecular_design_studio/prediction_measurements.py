"""Versioned measurement contract for installed co-folding predictors.

This describes outputs the local adapters can read. It does not promise that a
particular run emitted a score, and it does not set quality thresholds.
"""
from __future__ import annotations

VERSION = "studio.prediction-measurements.v1"

# Matrices are read from the prediction-specific sidecars used by the local
# result adapters. The shared interface calculation only runs when the matrix
# axes match the saved protein C-alpha inventory and binder/target chain roles.
PAE_MATRIX_ENGINES = {
    "alphafold3", "af3_torch", "openfold3_jax", "boltz2", "protenix_v1",
    "protenix_v2", "rf3_foundry", "esmfold2", "esmfold2_fast", "chai1",
    "atlasfold", "openfold3_ob0", "opendde_native", "colabfold",
    "boltzgen_fold", "mol1", "af3_jax_cofolding", "af2_initial_guess",
    "boltz2_initial_guess", "esmfold2_cofolding",
    "esmfold2_complex_validation", "esmfold2_initial_guess_validation",
    "openfold3", "rf3", "protenix",
}

# Exact output names accepted by the installed scalar-summary adapters. Values
# are conditional: a missing key is never replaced by a different confidence
# score. BoltzGen Fold's fields use its design-to-target nomenclature.
SCALAR_SCORE_ALIASES = {
    "alphafold3": {"ptm": ("ptm",), "iptm": ("iptm",)},
    "af3_torch": {"ptm": ("ptm",), "iptm": ("iptm",)},
    "openfold3_jax": {"ptm": ("ptm",), "iptm": ("iptm",)},
    "boltz2": {"ptm": ("ptm",), "iptm": ("iptm",)},
    "protenix_v1": {"ptm": ("ptm", "protenix_v1_ptm"), "iptm": ("iptm", "protenix_v1_iptm")},
    "protenix_v2": {"ptm": ("ptm", "protenix_v2_ptm"), "iptm": ("iptm", "protenix_v2_iptm")},
    "rf3_foundry": {"ptm": ("ptm",), "iptm": ("iptm",)},
    "esmfold2": {"ptm": ("ptm",), "iptm": ("iptm",)},
    "esmfold2_fast": {"ptm": ("ptm",), "iptm": ("iptm",)},
    "chai1": {"ptm": ("ptm",), "iptm": ("iptm",)},
    "atlasfold": {"ptm": ("ptm",), "iptm": ("iptm",)},
    "openfold3_ob0": {"ptm": ("ptm",), "iptm": ("iptm",)},
    "opendde_native": {"ptm": ("ptm",), "iptm": ("iptm",)},
    "colabfold": {"ptm": ("ptm", "colabfold_ptm"), "iptm": ("iptm", "colabfold_iptm")},
    "boltzgen_fold": {
        "ptm": ("design_ptm", "boltzgen_fold_design_ptm"),
        "iptm": ("design_to_target_iptm", "boltzgen_fold_design_to_target_iptm"),
    },
    "mol1": {"ptm": ("ptm",), "iptm": ("iptm",)},
    "af3_jax_cofolding": {"ptm": ("ptm",), "iptm": ("iptm",)},
    "af2_initial_guess": {"ptm": ("ptm",), "iptm": ("iptm",)},
    "boltz2_initial_guess": {
        "ptm": ("ptm", "boltz2_ptm"), "iptm": ("iptm", "boltz2_iptm"),
    },
    "esmfold2_cofolding": {"ptm": ("ptm",), "iptm": ("iptm",)},
    "esmfold2_complex_validation": {"ptm": ("ptm",), "iptm": ("iptm",)},
    "esmfold2_initial_guess_validation": {"ptm": ("ptm",), "iptm": ("iptm",)},
    "openfold3": {"ptm": ("openfold3_ptm",), "iptm": ("openfold3_iptm",)},
    "rf3": {"ptm": ("rf3_ptm",), "iptm": ("rf3_iptm",)},
    "protenix": {"ptm": ("protenix_ptm",), "iptm": ("protenix_iptm",)},
}

_PAE_SIDECAR_KEYS = (
    "predicted_aligned_error", "pae", "token_pair_pae", "predicted_aligned_errors",
)
_NATIVE_LABELS = {
    "af3_jax_cofolding": "AF3-family JAX co-folding route; model_name is retained per candidate",
    "af2_initial_guess": "AlphaFold2 initial-guess evaluator; native iPAE is emitted in Å",
    "boltz2_initial_guess": "Boltz-2 native complex-refolding route",
    "esmfold2_cofolding": "ESMFold2 co-folding route",
    "esmfold2_complex_validation": "ESMFold2 independent complex validation",
    "esmfold2_initial_guess_validation": "ESMFold2 initial-guess validation",
    "openfold3": "OpenFold3 native refolding",
    "rf3": "RosettaFold3 native refolding",
    "protenix": "Protenix native refolding",
    "bindcraft2": "BindCraft2 native AF2 design evaluation",
    "rfd4_proteina": "RFD4-Proteina native AF3 refolding oracle",
}


def _engine_labels() -> dict[str, dict]:
    """Read labels and versions from the actual installed Studio registry."""
    from mn_molecular_design_studio.engine.workflow_engine.workflows.optimized_prediction import ENGINE_PROFILES

    rows = {engine: {"label": profile["label"], "version": profile.get("code_version"),
                     "api_engine": profile.get("api_engine")}
            for engine, profile in ENGINE_PROFILES.items()}
    rows.update({key: {"label": label, "version": None, "api_engine": key}
                 for key, label in _NATIVE_LABELS.items()})
    return rows


def supported_score_engines(metric: str) -> tuple[str, ...]:
    return tuple(engine for engine, metrics in SCALAR_SCORE_ALIASES.items() if metric in metrics)


def interface_pae_from_matrix(matrix, binder_indices, target_indices) -> dict[str, float] | None:
    """Summarize saved PAE blocks in Å and on the conventional /31 scale.

    ``binder_to_target`` is PAE[binder, target], matching the row/column
    convention used by the installed ESMFold2 analysis adapter. The
    bidirectional value is the mean of the two directional means. These are
    arithmetic normalizations, not cross-model calibration.
    """
    try:
        import numpy as np

        values = np.asarray(matrix, dtype=float)
        binder = np.asarray(list(binder_indices), dtype=int)
        target = np.asarray(list(target_indices), dtype=int)
        if (values.ndim != 2 or values.shape[0] != values.shape[1]
                or not binder.size or not target.size
                or np.any(binder < 0) or np.any(target < 0)
                or np.any(binder >= values.shape[0]) or np.any(target >= values.shape[0])
                or not np.isfinite(values).all()):
            return None
        bt = values[np.ix_(binder, target)]
        tb = values[np.ix_(target, binder)]
        bt_mean, tb_mean = float(bt.mean()), float(tb.mean())
        bidirectional = (bt_mean + tb_mean) / 2.0
        result = {
            "prediction_ipae_binder_to_target": bt_mean,
            "prediction_ipae_target_to_binder": tb_mean,
            "prediction_ipae_bidirectional": bidirectional,
            "prediction_ipae_binder_to_target_norm31": bt_mean / 31.0,
            "prediction_ipae_target_to_binder_norm31": tb_mean / 31.0,
            "prediction_ipae_bidirectional_norm31": bidirectional / 31.0,
        }
        return {key: round(value, 6) for key, value in result.items()}
    except (ImportError, TypeError, ValueError, IndexError, OverflowError):
        return None


def prediction_engine_catalog() -> dict:
    """Machine-readable engine/output catalog included in the output vocabulary API."""
    labels = _engine_labels()
    engines = []
    for engine, identity in labels.items():
        aliases = SCALAR_SCORE_ALIASES.get(engine, {})
        matrix_available = engine in PAE_MATRIX_ENGINES
        plddt_note = (
            "Mean Cα B-factor read from the predicted structure; model-native confidence summaries, when present, are separate."
        )
        if engine == "mol1":
            plddt_note += " Mol-1 summary_confidence_sample plddt is retained separately from structure Cα B-factors."
        measurements = [
            {"metric": "pTM", "available": bool(aliases.get("ptm")), "availability": "conditional scalar output",
             "unit": "0–1", "source_fields": list(aliases.get("ptm", ())),
             "note": "Missing native scalar output remains unavailable."},
            {"metric": "iPTM", "available": bool(aliases.get("iptm")), "availability": "conditional scalar output",
             "unit": "0–1", "source_fields": list(aliases.get("iptm", ())),
             "note": "Missing native scalar output remains unavailable; monomer sentinel values are not interface measurements."},
            {"metric": "pLDDT", "available": True,
             "availability": "predicted-structure Cα B-factor; native scalar summary may differ",
             "unit": "0–100", "source_fields": ["structure Cα B-factor", "plddt_mean when emitted"],
             "note": plddt_note},
            {"metric": "PAE matrix", "available": matrix_available,
             "availability": "conditional saved per-model sidecar" if matrix_available else "not exposed by the installed adapter",
             "unit": "Å", "source_fields": list(_PAE_SIDECAR_KEYS) if matrix_available else [],
             "note": ("Studio derives binder→target, target→binder, and bidirectional inter-chain means only when matrix axes match the saved protein Cα inventory. Mol-1 requires atom-confidence output."
                      if matrix_available and engine == "mol1" else
                      "Studio derives binder→target, target→binder, and bidirectional inter-chain means only when matrix axes match the saved protein Cα inventory."
                      if matrix_available else "No PAE value is inferred from pTM, iPTM, or structure B-factors.")},
            {"metric": "ipSAE", "available": matrix_available,
             "availability": "separate Analysis operation from saved PAE" if matrix_available else "requires a compatible PAE matrix",
             "unit": "0–1", "source_fields": ["analysis.ipsae"],
             "note": "ipSAE is not a PAE conversion and is not treated as equivalent to iPTM."},
            {"metric": "interface contacts", "available": True,
             "availability": "Studio structure analysis when binder/target chain roles and coordinates are available",
             "unit": "residues / pairs", "source_fields": ["Studio heavy-atom contacts", "Studio interface residue count"],
             "note": "Geometry-derived; not a model confidence score."},
            {"metric": "RMSD", "available": True,
             "availability": "Studio comparison when a compatible reference structure and residue mapping are available",
             "unit": "Å", "source_fields": ["target_aligned_binder_rmsd", "binder_self_consistency_rmsd"],
             "note": "Target-aligned pose RMSD is unavailable for a sequence-only/IDR target without a fixed structural reference."},
        ]
        if engine == "bindcraft2":
            measurements = [
                {"metric": "pTM", "available": True, "availability": "conditional native AF2 score",
                 "unit": "0–1", "source_fields": ["pTM"], "note": "BindCraft2 design-stage confidence; not an independent refold."},
                {"metric": "iPTM", "available": True, "availability": "conditional native AF2 score",
                 "unit": "0–1", "source_fields": ["i_pTM"], "note": "BindCraft2 design-stage interface score; missing or monomer values are not independent measurements."},
                {"metric": "iPAE", "available": True, "availability": "conditional native scalar score",
                 "unit": "fraction of 31 Å PAE scale", "source_fields": ["i_pAE"],
                 "note": "Binder→target directional mean PAE already divided by 31 Å. Studio can multiply by 31 to show Å; this is a scale conversion, not predictor calibration."},
                {"metric": "i_pDAE", "available": True, "availability": "conditional native ranking score",
                 "unit": "0–1", "source_fields": ["i_pDAE"],
                 "note": "PAE-derived contact score; distinct from iPAE and retained as a separate metric."},
                {"metric": "pLDDT", "available": True, "availability": "conditional native AF2 score",
                 "unit": "0–100", "source_fields": ["pLDDT", "Target_pLDDT", "Unbound_Binder_pLDDT"],
                 "note": "Bound binder, target, and unbound binder confidence are separate scopes."},
                {"metric": "PAE matrix", "available": False, "availability": "not exposed as a common matrix by the installed adapter",
                 "unit": "Å", "source_fields": [], "note": "The native scalar i_pAE remains available with its documented direction and /31 scale."},
                {"metric": "ipSAE", "available": False, "availability": "separate analysis not emitted by native BindCraft2 filter",
                 "unit": "0–1", "source_fields": [], "note": "Not a conversion of iPAE."},
                {"metric": "interface contacts", "available": True,
                 "availability": "Studio structure analysis when binder/target roles and coordinates are available",
                 "unit": "residues / pairs", "source_fields": ["Studio heavy-atom contacts", "Studio interface residue count"],
                 "note": "Geometry-derived; not a model confidence score."},
                {"metric": "RMSD", "available": True,
                 "availability": "Studio comparison when a compatible reference structure and residue mapping are available",
                 "unit": "Å", "source_fields": ["Binder_RMSD", "target_aligned_binder_rmsd"],
                 "note": "BindCraft2 Binder_RMSD and Studio target-aligned pose RMSD use different references and remain distinct."},
            ]
        elif engine == "rfd4_proteina":
            measurements = [
                {"metric": "pTM", "available": True, "availability": "conditional RFD4 native refolding-oracle output",
                 "unit": "0–1", "source_fields": ["oracle_ppi.overall_ptm", "oracle_ppi.ptm_coi"],
                 "note": "Values in the native flat output are percent-scaled and normalized by Studio to 0–1. The selected folding_engine is retained per candidate."},
                {"metric": "iPTM", "available": True, "availability": "conditional RFD4 native refolding-oracle output",
                 "unit": "0–1", "source_fields": ["oracle_ppi.overall_iptm", "oracle_ppi.iptm_coi"],
                 "note": "AF3 oracle output; overall and chain-of-interest scopes remain separate."},
                {"metric": "iPAE", "available": True, "availability": "conditional native scalar score",
                 "unit": "Å", "source_fields": ["oracle_ppi.ipae_coi"],
                 "note": "RFD4 emits this through its native oracle metric contract: average forward/reverse chain-pair PAE for each partner, then average partner chains for the chain of interest. The candidate stores folding_engine; Studio's /31 field is scale conversion, not calibration."},
                {"metric": "pLDDT", "available": True, "availability": "conditional RFD4 native refolding-oracle output",
                 "unit": "0–100", "source_fields": ["oracle_ppi.overall_plddt", "oracle_ppi.plddt_coi"],
                 "note": "Whole-complex and chain-of-interest values are separate."},
                {"metric": "PAE matrix", "available": False, "availability": "not exposed by the current RFD4 candidate adapter",
                 "unit": "Å", "source_fields": [], "note": "The native oracle's scalar ipae_coi remains available."},
                {"metric": "ipSAE", "available": False, "availability": "not emitted by the current RFD4 oracle adapter",
                 "unit": "0–1", "source_fields": [], "note": "Not a conversion of iPAE."},
                {"metric": "interface contacts", "available": True,
                 "availability": "Studio structure analysis when binder/target roles and coordinates are available",
                 "unit": "residues / pairs", "source_fields": ["Studio heavy-atom contacts", "Studio interface residue count"],
                 "note": "Geometry-derived; not an AF3 confidence score."},
                {"metric": "RMSD", "available": True,
                 "availability": "Studio comparison when a compatible reference structure and residue mapping are available",
                 "unit": "Å", "source_fields": ["refolded_target_aligned_binder_rmsd", "target_aligned_binder_rmsd"],
                 "note": "Target-aligned pose RMSD requires a fixed structural reference; sequence-only/IDR targets do not supply one."},
            ]
        engines.append({"engine": engine, **identity, "measurements": measurements})
    return {"version": VERSION, "normalization": {
        "raw_pae": "Å, preserved with producer and aggregation provenance",
        "normalized_pae": "raw PAE / 31 Å; a shared numeric scale, not a model calibration",
        "directional_fields": ["prediction_ipae_binder_to_target", "prediction_ipae_target_to_binder"],
        "bidirectional_field": "prediction_ipae_bidirectional",
        "missing_values": "remain missing and cannot pass a configured mandatory filter",
    }, "engines": engines}
