"""Versioned scientific output vocabulary; native results remain authoritative.

Aliases are explicit and scope aware. No fuzzy suffix matching, magnitude-based
unit guessing, or substitution of one model's score for another model's score.
This module has no workflow/Database dependencies so both can use the contract.
"""
from __future__ import annotations

import json
import math
import re
from dataclasses import asdict, dataclass
from functools import lru_cache

VERSION = "studio.output.v1"
PREFIX = "common."
ENGINE_ALIASES = {
    "rfdiffusion_generator": "rfdiffusion", "rfdiffusion3_foundry": "rfdiffusion3",
    "rfdiffusion_classic": "rfdiffusion", "bindcraft2_generator": "bindcraft2",
    "rfdiffusion_beta_pairing": "rfdiffusion", "protpardelle_1c": "protpardelle",
    "esmfold2_binder_design": "esmfold2_binder",
    "pxdesign_generator": "pxdesign", "genie3_generator": "genie3",
    "protpardelle_generator": "protpardelle", "boltzgen_generator": "boltzgen",
    "protein_mpnn": "proteinmpnn", "ligand_mpnn": "ligandmpnn",
    "soluble_mpnn": "solublempnn", "ab_mpnn": "abmpnn", "hyper_mpnn": "hypermpnn",
    "esmfold2_cofolding": "esmfold2", "esmfold2_complex_validation": "esmfold2",
    "esmfold2_initial_guess_validation": "esmfold2",
}
MPNN_ENGINES = {"proteinmpnn", "ligandmpnn", "solublempnn", "abmpnn", "hypermpnn"}
GENERATOR_WORKFLOWS = {
    "design.rfdiffusion_generate": "rfdiffusion", "design.rfdiffusion_beta_pairing": "rfdiffusion",
    "design.rfdiffusion3": "rfdiffusion3", "design.bindcraft2": "bindcraft2",
    "design.bindcraft2_generator": "bindcraft2", "design.boltzgen": "boltzgen",
    "design.boltzgen_generator": "boltzgen", "design.pxdesign_generate": "pxdesign",
    "design.rfd4_proteina": "rfd4_proteina",
    "design.genie3_generate": "genie3", "design.protpardelle_generate": "protpardelle",
    "design.proteina_complexa": "proteina_complexa", "design.proteina_complexa_generate": "proteina_complexa",
    "design.rfantibody": "rfantibody", "design.rfantibody_mimosa": "rfantibody",
    "design.boltzgen_mimosa": "boltzgen", "design.rfpeptides": "rfpeptides",
    "design.esmfold2_binder": "esmfold2_binder",
    "design.kejia_idr": "kejia_idr", "design.peptide_trap": "peptide_trap",
}
FASTA_METRICS = {"score", "global_score", "seq_recovery", "seq_rec", "overall_confidence", "ligand_confidence"}
_HEADER_NUMBER = re.compile(r"(?:^|[,\s])([A-Za-z_]+)\s*=\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)(?=,|\s|$)")


def engine_name(tool):
    name = str(tool or "").strip().lower().replace("-", "_")
    if name.startswith("optimized_prediction_"):
        name = name.removeprefix("optimized_prediction_")
    return ENGINE_ALIASES.get(name, name)


def fasta_metrics(header):
    """Read exact FASTA header keys, including the current LigandMPNN format."""
    return {key: float(value) for key, value in _HEADER_NUMBER.findall(str(header or ""))
            if key in FASTA_METRICS and math.isfinite(float(value))}


def recovered_metrics(metrics, raw, tool):
    result = dict(metrics or {})
    if engine_name(tool) == "mol1":
        confidence = (raw or {}).get("confidence_outputs")
        if isinstance(confidence, dict):
            plddt = confidence.get("plddt_mean", confidence.get("plddt"))
            if plddt is not None:
                result.setdefault("mol1_model_plddt_mean", plddt)
    if engine_name(tool) in MPNN_ENGINES:
        for key, value in fasta_metrics((raw or {}).get("fasta_header")).items():
            result.setdefault(key, value)
    return result


def _candidate_stage(candidate):
    """Use molecular input composition to correct legacy monomer result labels."""
    stage = candidate.get("stage")
    if stage != "complex_refolding":
        return stage
    raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
    chains = raw.get("protein_chains") or raw.get("chain_sequences") or []
    count = len(chains) if isinstance(chains, (list, dict)) else 0
    for key in ("dna_chains", "rna_chains", "ligand_ccd_codes", "ligand_smiles"):
        values = raw.get(key) or []
        count += len(values) if isinstance(values, (list, dict)) else 0
    if raw.get("has_target_structure"):
        count += 1
    return "monomer_refolding" if count == 1 else stage


@dataclass(frozen=True)
class Measure:
    key: str
    label: str
    unit: str
    scope: str
    method: str
    aliases: tuple[str, ...]
    direction: str = ""
    minimum: float | None = 0
    maximum: float | None = None
    engines: tuple[str, ...] = ()
    transform: str = "identity"
    stages: tuple[str, ...] = ()
    requires: tuple[tuple[str, object], ...] = ()
    description: str = ""
    value_type: str = "number"
    prerequisites: tuple[str, ...] = ()
    aggregation: str = "scalar"


# Shared measurements computed by Studio for every generator. Native output
# names such as RF3 radius_of_gyration / helix_fraction have different scopes
# and are deliberately absent here.
MEASURES = [
    Measure("design.sequence_length", "Designed sequence length", "residues", "designed chains", "sequence", ("$sequence_length",)),
    Measure("geometry.binder.residue_count", "Binder coordinate residues", "residues", "binder", "coordinate inventory", ("binder_residue_count",)),
    Measure("geometry.binder.ca_count", "Binder Cα atoms", "atoms", "binder", "coordinate inventory", ("binder_ca_count",)),
    Measure("geometry.binder.heavy_atom_count", "Binder heavy atoms", "atoms", "binder", "coordinate inventory", ("binder_heavy_atom_count",)),
    Measure("geometry.binder.heavy_atom_rg", "Binder heavy-atom radius of gyration", "Å", "binder", "Studio unweighted heavy atoms", ("binder_radius_of_gyration",), requires=(("binder_geometry_status", "ok"),)),
    Measure("geometry.binder.ca_rg", "Binder Cα radius of gyration", "Å", "binder", "Studio Cα atoms", ("binder_ca_radius_of_gyration",)),
    Measure("geometry.binder.ca_end_to_end", "Binder Cα end-to-end distance", "Å", "binder", "Studio first/last Cα", ("binder_ca_end_to_end_distance", "binder_end_to_end_distance")),
    Measure("geometry.target.alignment_rmsd", "Retained target alignment RMSD", "Å", "target", "Studio matched Cα superposition", ("generation_target_alignment_rmsd",), "lower"),
    Measure("mimosa.generated_motif_rmsd", "Generated motif RMSD", "Å", "guided generated motif", "MIMOSA target-aligned Cα fit", ("generated_motif_rmsd_angstrom",), "lower", engines=("boltzgen", "rfantibody"), stages=("generation.backbone_sequence", "generation.backbone")),
    Measure("mimosa.refolded_motif_rmsd", "Refolded motif RMSD", "Å", "refolded motif", "MIMOSA target-aligned Cα fit", ("refolded_motif_rmsd_angstrom",), "lower", engines=("boltzgen", "rfantibody"), stages=("complex_refolding",)),
    Measure("mimosa.generated_potential", "Generated motif potential", "Å²", "guided generated motif", "MIMOSA half-sum squared Cα displacement", ("generated_motif_potential_half_sum_angstrom2",), "lower", minimum=0, engines=("boltzgen", "rfantibody"), stages=("generation.backbone_sequence", "generation.backbone")),
    Measure("mimosa.fixed_motif_residues", "Fixed motif residue count", "residues", "generated motif", "MIMOSA strict <2 Å identity fixing", ("fixed_motif_residue_count",), minimum=0, engines=("boltzgen", "rfantibody"), stages=("generation.backbone_sequence", "generation.backbone")),
    Measure("mimosa.backbone_integrity", "Backbone integrity status", "status", "generated binder backbone", "Consecutive polymer C–N distance, target 1.329±0.5 Å", ("backbone_integrity_status",), value_type="text", engines=("boltzgen", "rfantibody"), stages=("generation.backbone_sequence", "generation.backbone")),
    Measure("mimosa.backbone_bad_bond_count", "Backbone bad C–N bond count", "bonds", "generated binder backbone", "Consecutive polymer C–N distance, target 1.329±0.5 Å", ("backbone_bad_bond_count",), minimum=0, engines=("boltzgen", "rfantibody"), stages=("generation.backbone_sequence", "generation.backbone")),
]
for component in ("helix", "sheet", "coil", "structured"):
    MEASURES.append(Measure(f"secondary_structure.binder.dssp.{component}_fraction", f"Binder {component} fraction (DSSP)", "fraction", "binder", "Studio pydssp C3", (f"binder_{component}_fraction",), maximum=1))
for component in ("helix_like", "extended", "other"):
    MEASURES.append(Measure(f"secondary_structure.binder.ca_trace.{component}_fraction", f"Binder {component.replace('_', ' ')} fraction (Cα trace)", "fraction", "binder", "Studio Cα trace heuristic", (f"binder_ca_trace_{component}_fraction",), maximum=1))

# A trajectory's prediction is not an independent validation of its design.
# Only audited, explicitly scoped native scores enter these confidence families.
for role, stages in (("generator", ("generation.backbone_sequence",)),
                     ("validation", ("complex_refolding",))):
    for engine, aliases, method in (
        ("bindcraft2", (("ptm", "pTM"), ("iptm", "i_pTM")), "BindCraft2 AF2"),
        ("boltzgen", (("ptm", "ptm"), ("iptm", "iptm")), "BoltzGen native predictor"),
        ("rfdiffusion3", (("ptm", "ptm"), ("iptm", "iptm")), "RFdiffusion3 native refolding"),
        ("genie3", (("ptm", "ptm"), ("iptm", "iptm")), "Genie3 native refolding"),
        ("pxdesign", (("ptm", "af2_ptm"), ("iptm", "iptm")), "PXDesign native AF2-IG"),
        ("proteina_complexa", (("ptm", "ptm"), ("iptm", "iptm")), "Proteina-Complexa native evaluator"),
    ):
        for canonical, alias in aliases:
            MEASURES.append(Measure(f"confidence.{role}.complex.{canonical}", f"{role.title()} complex {canonical}", "0–1", "complex" if canonical == "ptm" else "complex interfaces", method, (alias,), "higher", maximum=1, engines=(engine,), stages=stages))
    for key, label, scope, alias in (
        ("bound_binder.plddt", "bound binder pLDDT", "binder in complex", "pLDDT"),
        ("target.plddt", "target pLDDT", "target in complex", "Target_pLDDT"),
        ("unbound_binder.plddt", "unbound binder pLDDT", "binder alone", "Unbound_Binder_pLDDT"),
    ):
        MEASURES.append(Measure(f"confidence.{role}.{key}", f"{role.title()} {label}", "0–100", scope, "BindCraft2 AF2", (alias,), "higher", maximum=100, engines=("bindcraft2",), transform="fraction_to_percent", stages=stages))
    MEASURES.append(Measure(f"confidence.{role}.interface.mean_pae", f"{role.title()} binder–target mean PAE", "Å", "binder-to-target directional mean", "BindCraft2 AF2", ("i_pAE",), "lower", engines=("bindcraft2",), transform="times_31", stages=stages))
# BindCraft2 stores its binder-to-target mean PAE after division by 31 Å.
# RFD4's native oracle emits an Å-valued COI/partner mean that averages both
# directions. Keep raw producer metrics and expose the numeric /31 scale;
# provenance retains aggregation and the selected oracle engine.
MEASURES.extend([
    Measure("confidence.generator.interface.pae_norm31", "Generator interface PAE normalized by 31 Å", "PAE / 31 Å", "producer-defined interchain PAE aggregation",
        "BindCraft2 AF2 i_pAE: binder-to-target directional mean, already normalized by 31 Å", ("i_pAE",), "lower", minimum=0, maximum=1, engines=("bindcraft2",)),
    Measure("confidence.generator.interface.pae_norm31", "Generator interface PAE normalized by 31 Å", "PAE / 31 Å", "producer-defined interchain PAE aggregation",
        "RFD4 native PPI oracle ipae_coi: COI-to-partner mean averaged across both directions, divided by 31 Å; folding_engine is retained per run", ("oracle_ppi.ipae_coi",), "lower", minimum=0,
        maximum=1, engines=("rfd4_proteina",), transform="divide_31", stages=("generation.backbone_sequence",)),
])
# RFD4's oracle flattens AF3-shaped confidence outputs after scaling pLDDT/pTM by
# 100; PAE remains in Å. These aliases are scoped to RFD4's generator stage,
# since trajectory confidence is not an independent validation run.
for scope, suffix, aliases in (
    ("complex", "", (("oracle.overall_plddt", "oracle_ppi.overall_plddt"),
                       ("oracle.overall_pae", "oracle_ppi.overall_pae"),
                       ("oracle.overall_ptm", "oracle_ppi.overall_ptm"))),
    ("chain of interest", ".chain_of_interest", (("oracle.plddt_coi", "oracle_ppi.plddt_coi"),
                                                    ("oracle.pae_coi", "oracle_ppi.pae_coi"),
                                                    ("oracle.ptm_coi", "oracle_ppi.ptm_coi"))),
):
    for quantity, unit, direction, transform, maximum, alias in zip(
        ("plddt", "pae", "ptm"), ("0–100", "Å", "0–1"),
        ("higher", "lower", "higher"), ("identity", "identity", "percent_to_fraction"),
        (100, None, 1), aliases,
    ):
        for source_alias, method in reversed(tuple(zip(
            alias,
            ("RFD4 native oracle (backend in candidate metadata)", "RFD4 native PPI oracle (backend in candidate metadata)"),
        ))):
            MEASURES.append(Measure(
                f"confidence.generator.{scope.replace(' ', '_')}.{quantity}",
                f"Generator {scope} {quantity}", unit, scope,
                method, (source_alias,), direction, maximum=maximum,
                engines=("rfd4_proteina",), transform=transform,
                stages=("generation.backbone_sequence",),
            ))
# RFD4-PPI emits confidence fields under ``oracle_ppi``. Keep these
# measurements explicitly native and namespaced;
# they are not independent re-prediction results. Confidence values from this
# provider are percentages and are converted only for pTM/iPTM fields.
for quantity, label, unit, direction, maximum, alias, transform in (
    ("iptm", "Generator complex iPTM", "0–1", "higher", 1,
     "oracle_ppi.overall_iptm", "percent_to_fraction"),
    ("ipae_coi", "Generator chain-of-interest interface PAE", "Å", "lower", None,
     "oracle_ppi.ipae_coi", "identity"),
):
    MEASURES.append(Measure(
        f"confidence.generator.complex.{quantity}" if quantity == "iptm" else f"confidence.generator.interface.{quantity}",
        label, unit, "complex interfaces" if quantity == "iptm" else "RFD4 chain of interest",
        "RFD4 native PPI oracle (backend in candidate metadata)", (alias,), direction,
        maximum=maximum, engines=("rfd4_proteina",), transform=transform,
        stages=("generation.backbone_sequence",),
    ))
MEASURES.extend([
    Measure("interface.native.rfd4_proteina.binder_residue_count_4a",
        "RFD4 saved complex · binder residues contacting target within 4 Å", "residues",
        "saved RFD4 generated complex", "Studio heavy-atom distance; each binder residue counted once",
        ("rfd4_native_interface_residue_count_4a",), minimum=0,
        engines=("rfd4_proteina",), stages=("generation.backbone_sequence",)),
    Measure("geometry.native.rfd4_proteina.backbone_ca_clashes_2_5a",
        "RFD4 saved complex · interchain Cα clashes within 2.5 Å", "atom pairs",
        "saved RFD4 generated complex", "Studio interchain Cα pair count at ≤2.5 Å",
        ("rfd4_native_backbone_ca_clashes_2_5a",), minimum=0,
        engines=("rfd4_proteina",), stages=("generation.backbone_sequence",)),
])
for model in sorted(MPNN_ENGINES):
    for key, label, aliases, unit, direction, transform, maximum in (
        ("redesigned_nll", "redesigned-residue NLL", ("score",), "nats/residue", "lower", "identity", None),
        ("redesigned_nll", "redesigned-residue NLL", ("overall_confidence",), "nats/residue", "lower", "negative_log", None),
        ("global_nll", "global NLL", ("global_score",), "nats/residue", "lower", "identity", None),
        ("redesigned_probability", "redesigned-residue probability", ("overall_confidence",), "fraction", "higher", "identity", 1),
    ):
        MEASURES.append(Measure(f"sequence.{model}.{key}", f"{model} {label}", unit, "global" if key == "global_nll" else "redesigned residues near ligand" if key == "ligand_probability" else "redesigned residues", model, aliases, direction, maximum=maximum, engines=(model,), transform=transform))
MEASURES.extend([
    Measure("diversity.representative", "Selected cluster representative", "0/1", "comparison cohort", "declared clustering method", ("diversity_representative",), maximum=1, transform="bool_to_indicator",
            description="Source-order cluster representative, for explicit selection; this is not a quality score."),
    Measure("sequence.recovery_designed", "Sequence recovery over redesigned residues", "fraction", "redesigned residues", "identity to input sequence", ("seq_recovery", "seq_rec"), maximum=1, engines=tuple(sorted(MPNN_ENGINES))),
    Measure("sequence.esm_if1.log_likelihood", "ESM-IF1 log likelihood", "native", "designed sequence", "ESM-IF1", ("esm_if1_log_likelihood",), "higher", minimum=None, engines=("esm_if1",)),
    Measure("sequence.caliby.energy", "Caliby energy", "native", "designed sequence", "Caliby", ("caliby_energy",), "lower", minimum=None, engines=("caliby",)),
])
# Independent predictors retain the same canonical pTM/iPTM families, but each
# source keeps its exact result key and method. See prediction_measurements for
# the locally installed engines and their full matrix/output capabilities.
from mn_molecular_design_studio.prediction_measurements import (
    PAE_MATRIX_ENGINES, SCALAR_SCORE_ALIASES, prediction_engine_catalog,
)

for engine, scores in SCALAR_SCORE_ALIASES.items():
    for quantity in ("ptm", "iptm"):
        aliases = scores.get(quantity, ())
        if aliases:
            MEASURES.append(Measure(
                f"confidence.validation.complex.{quantity}", f"Validation complex {quantity}", "0–1",
                "complex interfaces" if quantity == "iptm" else "complex",
                f"{engine} native prediction confidence", aliases, "higher", maximum=1,
                engines=(engine,), stages=("complex_refolding",),
            ))
    aliases = scores.get("ptm", ())
    if aliases:
        MEASURES.append(Measure(
            "confidence.validation.monomer.ptm", "Validation monomer pTM", "0–1",
            "single modeled component", f"{engine} native prediction confidence", aliases,
            "higher", maximum=1, engines=(engine,), stages=("monomer_refolding",),
        ))
# The shared prediction adapters emit scalar pTM/iPTM for these routes. Their
# result-level key can be generic or prefixed by the native route name.
for engine, aliases in {
    "optimized_prediction": {"ptm": ("ptm",), "iptm": ("iptm",)},
}.items():
    for quantity, keys in aliases.items():
        MEASURES.append(Measure(
            f"confidence.validation.complex.{quantity}", f"Validation complex {quantity}", "0–1",
            "complex interfaces" if quantity == "iptm" else "complex",
            f"{engine} native prediction confidence", keys, "higher", maximum=1,
            engines=(engine,), stages=("complex_refolding",),
        ))
MEASURES.append(Measure("confidence.validation.structure.mean_ca_plddt", "Mean structure Cα pLDDT", "0–100", "all protein Cα atoms", "Studio mean Cα B-factor read from output PDB/mmCIF",
    ("plddt_mean", "complex_plddt"), "higher", maximum=100, stages=("complex_refolding", "monomer_refolding")))
MEASURES.extend([
    Measure("confidence.validation.structure.binder_ca_plddt", "Bound binder mean Cα B-factor pLDDT", "0–100", "binder Cα atoms in predicted complex", "Studio mean Cα B-factor read from output PDB/mmCIF and binder chain roles",
        ("binder_structure_ca_plddt",), "higher", maximum=100, stages=("complex_refolding",)),
    Measure("confidence.validation.structure.target_ca_plddt", "Target mean Cα B-factor pLDDT", "0–100", "target Cα atoms in predicted complex", "Studio mean Cα B-factor read from output PDB/mmCIF and target chain roles",
        ("target_structure_ca_plddt",), "higher", maximum=100, stages=("complex_refolding",)),
])
MEASURES.append(Measure("confidence.validation.model.mol1_plddt_mean", "Mol-1 model summary pLDDT", "0–100", "all modeled protein chains", "Mol-1 native summary_confidence_sample plddt",
    ("mol1_model_plddt_mean",), "higher", maximum=100, engines=("mol1",), stages=("complex_refolding", "monomer_refolding")))
MEASURES.append(Measure("confidence.validation.bound_binder.plddt", "Validation bound binder pLDDT", "0–100", "binder in complex", "Producer model bound-binder confidence; source tool and protocol retained",
    ("binder_plddt", "af2_binder_plddt", "af2_average_binder_plddt"), "higher", maximum=100,
    stages=("complex_refolding",)))
MEASURES.extend([
    Measure("confidence.validation.interface.mean_pae", "Validation interface PAE", "Å", "producer-defined binder–target interface", "Producer-reported interface PAE; source tool and aggregation retained",
        ("ipae", "i_pae", "binder_pae", "boltz2_ipae"), "lower", minimum=0, stages=("complex_refolding",)),
    Measure("geometry.validation.binder.target_aligned_rmsd", "Target-aligned binder RMSD", "Å", "binder backbone in target frame", "Studio matched target Cα alignment and binder Cα RMSD",
        ("target_aligned_binder_rmsd", "binder_rmsd_target_aligned"), "lower", minimum=0, stages=("complex_refolding",)),
    Measure("geometry.validation.binder.self_consistency_rmsd", "Binder self-consistency RMSD", "Å", "binder backbone after binder-only alignment", "Studio matched binder Cα alignment; target pose is not assessed",
        ("binder_self_consistency_rmsd",), "lower", minimum=0, stages=("complex_refolding", "monomer_refolding")),
    Measure("geometry.validation.cyclic_peptide.closure_distance", "Cyclic peptide head-to-tail closure distance", "Å", "terminal backbone C–N atoms", "Studio output structure distance from the first binder N to the final binder C atom",
        ("cyclic_closure_distance",), "lower", minimum=0, engines=("boltz2",), stages=("complex_refolding", "monomer_refolding")),
    Measure("interface.validation.binder_residue_count_4a", "Binder residues contacting target within 4 Å", "residues", "binder–target interface", "Studio heavy-atom distances; count each binder residue once if any heavy atom is within 4 Å of any target heavy atom",
        ("validation_interface_residue_count_4a",), minimum=0, stages=("complex_refolding",)),
    Measure("geometry.validation.interface.backbone_ca_clashes_2_5a", "Interchain Cα clashes within 2.5 Å", "atom pairs", "predicted multichain structure", "Studio interchain Cα pair count at ≤2.5 Å; same distance definition as BindCraft2 Backbone_Clashes",
        ("validation_backbone_ca_clashes_2_5a",), minimum=0, stages=("complex_refolding",)),
            Measure("confidence.native_design.esmfold2.critic.iptm", "ESMFold2 native design critic iPTM", "0–1", "complex interfaces", "ESMFold2 gradient-design critic; native objective, not independent validation",
        ("iptm",), "higher", maximum=1, engines=("esmfold2_binder",), stages=("complex_refolding",)),
])
# Values calculated from a saved, prediction-specific PAE matrix. Keep raw Å
# results, each direction, and the shared /31 numeric scale as separate fields.
for engine in sorted(PAE_MATRIX_ENGINES):
    method = (f"{engine} saved PAE matrix; arithmetic mean of PAE[binder,target] over all inter-chain residues/tokens"
              if engine != "boltzgen_fold" else
              "BoltzGen Fold saved PAE matrix; arithmetic mean of PAE[binder,target] over all inter-chain residues/tokens")
    for key, label, unit, alias in (
        ("confidence.validation.interface.pae_binder_to_target", "Validation binder→target mean PAE", "Å", "prediction_ipae_binder_to_target"),
        ("confidence.validation.interface.pae_target_to_binder", "Validation target→binder mean PAE", "Å", "prediction_ipae_target_to_binder"),
        ("confidence.validation.interface.pae_bidirectional", "Validation bidirectional mean PAE", "Å", "prediction_ipae_bidirectional"),
        ("confidence.validation.interface.pae_binder_to_target_norm31", "Validation binder→target mean PAE / 31 Å", "PAE / 31 Å", "prediction_ipae_binder_to_target_norm31"),
        ("confidence.validation.interface.pae_target_to_binder_norm31", "Validation target→binder mean PAE / 31 Å", "PAE / 31 Å", "prediction_ipae_target_to_binder_norm31"),
        ("confidence.validation.interface.pae_bidirectional_norm31", "Validation bidirectional mean PAE / 31 Å", "PAE / 31 Å", "prediction_ipae_bidirectional_norm31"),
    ):
        measure_method = method + ("; divide by 31 Å" if unit == "PAE / 31 Å" else "; raw PAE in Å")
        MEASURES.append(Measure(key, label, unit, "all binder×target PAE entries" if "binder_to_target" in key or "target_to_binder" in key else "mean of both inter-chain directions",
            measure_method, (alias,), "lower", minimum=0, engines=(engine,), stages=("complex_refolding",)))
MEASURES.extend([
    Measure("confidence.analysis.interface.pae_binder_to_target", "ipSAE analysis contact PAE · binder→target", "Å", "Cα-contacting binder×target residue pairs", "Saved PAE[binder,target] over Cα contacts within the configured ipSAE distance cutoff",
        ("ipae_binder_to_target",), "lower", minimum=0, engines=("ipsae",), stages=("analysis",)),
    Measure("confidence.analysis.interface.pae_target_to_binder", "ipSAE analysis contact PAE · target→binder", "Å", "Cα-contacting target×binder residue pairs", "Saved PAE[target,binder] over Cα contacts within the configured ipSAE distance cutoff",
        ("ipae_target_to_binder",), "lower", minimum=0, engines=("ipsae",), stages=("analysis",)),
    Measure("confidence.analysis.interface.pae_bidirectional", "ipSAE analysis bidirectional contact PAE", "Å", "mean of directional Cα-contact PAE means", "Arithmetic mean of the two contact-scoped directional PAE means",
        ("ipae_bidirectional",), "lower", minimum=0, engines=("ipsae",), stages=("analysis",)),
    Measure("confidence.analysis.interface.pae_binder_to_target_norm31", "ipSAE analysis binder→target contact PAE / 31 Å", "PAE / 31 Å", "Cα-contacting binder×target residue pairs", "Contact PAE[binder,target] mean divided by 31 Å; scale conversion, not model calibration",
        ("ipae_binder_to_target_norm31",), "lower", minimum=0, engines=("ipsae",), stages=("analysis",)),
    Measure("confidence.analysis.interface.pae_target_to_binder_norm31", "ipSAE analysis target→binder contact PAE / 31 Å", "PAE / 31 Å", "Cα-contacting target×binder residue pairs", "Contact PAE[target,binder] mean divided by 31 Å; scale conversion, not model calibration",
        ("ipae_target_to_binder_norm31",), "lower", minimum=0, engines=("ipsae",), stages=("analysis",)),
    Measure("confidence.analysis.interface.pae_bidirectional_norm31", "ipSAE analysis bidirectional contact PAE / 31 Å", "PAE / 31 Å", "mean of directional Cα-contact PAE means", "Bidirectional contact PAE divided by 31 Å; scale conversion, not model calibration",
        ("ipae_bidirectional_norm31",), "lower", minimum=0, engines=("ipsae",), stages=("analysis",)),
    Measure("interface.analysis.bindcraft2_definition.ipdae", "BC2-definition i_pDAE · Studio recomputed", "0–1", "Cα-contacting binder–target interface", "BindCraft2 anchored interface TM formula recomputed from saved PAE and Cα contacts; max over binder- and target-anchored residue scores",
        ("bindcraft2_ipdae_recomputed",), "higher", minimum=0, maximum=1, engines=("ipsae",), stages=("analysis",)),
    Measure("interface.analysis.bindcraft2_definition.ipdae_binder_anchor", "BC2-definition i_pDAE · binder anchored", "0–1", "binder interface residues", "Maximum BindCraft2 anchored interface TM score among contacting binder residues",
        ("bindcraft2_ipdae_binder_anchor_max",), "higher", minimum=0, maximum=1, engines=("ipsae",), stages=("analysis",)),
    Measure("interface.analysis.bindcraft2_definition.ipdae_target_anchor", "BC2-definition i_pDAE · target anchored", "0–1", "target interface residues", "Maximum BindCraft2 anchored interface TM score among contacting target residues",
        ("bindcraft2_ipdae_target_anchor_max",), "higher", minimum=0, maximum=1, engines=("ipsae",), stages=("analysis",)),
    Measure("interface.analysis.ipsae_structure.atom_contact_pairs", "ipSAE analysis · inter-chain heavy-atom contact pairs", "pairs", "selected binder and target chains", "Biopython NeighborSearch pair count under configured heavy-atom distance cutoff",
        ("structure_interface_atom_contact_pair_count",), "higher", minimum=0, engines=("ipsae",), stages=("analysis",)),
    Measure("interface.analysis.ipsae_structure.residue_contact_pairs", "ipSAE analysis · inter-chain residue contact pairs", "pairs", "selected binder and target chains", "Unique residue pairs with any heavy-atom contact under configured distance cutoff",
        ("structure_interface_residue_contact_pair_count",), "higher", minimum=0, engines=("ipsae",), stages=("analysis",)),
    Measure("interface.analysis.ipsae_structure.binder_interface_residues", "ipSAE analysis · binder interface residues", "residues", "selected binder chains", "Unique binder residues with any heavy-atom contact under configured distance cutoff",
        ("structure_binder_interface_residue_count",), "higher", minimum=0, engines=("ipsae",), stages=("analysis",)),
    Measure("interface.analysis.ipsae_structure.target_interface_residues", "ipSAE analysis · target interface residues", "residues", "selected target chains", "Unique target residues with any heavy-atom contact under configured distance cutoff",
        ("structure_target_interface_residue_count",), "higher", minimum=0, engines=("ipsae",), stages=("analysis",)),
    Measure("interface.analysis.ipsae_structure.design_contact_pairs", "Designed complex · residue contact pairs", "pairs", "separate source design complex", "Unique binder–target residue pairs with heavy-atom contacts under the shared configured distance cutoff",
        ("structure_design_interface_residue_contact_pairs",), "higher", minimum=0, engines=("ipsae",), stages=("analysis",)),
    Measure("interface.analysis.ipsae_structure.cofold_contact_pairs", "Cofold complex · residue contact pairs", "pairs", "saved predicted complex", "Unique binder–target residue pairs with heavy-atom contacts under the shared configured distance cutoff",
        ("structure_cofold_interface_residue_contact_pairs",), "higher", minimum=0, engines=("ipsae",), stages=("analysis",)),
    Measure("interface.analysis.ipsae_structure.preserved_contact_pairs", "Designed contact pairs retained in cofold", "pairs", "sequence-mapped binder–target interface", "Designed residue-contact pairs also found after sequence mapping in the cofolded complex",
        ("structure_interface_contact_pairs_preserved",), "higher", minimum=0, engines=("ipsae",), stages=("analysis",)),
    Measure("interface.analysis.ipsae_structure.design_contact_preservation_fraction", "Designed interface contact preservation fraction", "fraction", "sequence-mapped designed residue-contact pairs", "Fraction of designed binder–target residue-contact pairs also contacting in the cofold",
        ("structure_design_interaction_preservation_fraction",), "higher", minimum=0, maximum=1, engines=("ipsae",), stages=("analysis",)),
    Measure("interface.analysis.ipsae_structure.contact_jaccard", "Designed/cofold interface contact Jaccard", "fraction", "sequence-mapped residue-contact pair sets", "Intersection-over-union of designed and cofold residue-contact pairs",
        ("structure_interface_contact_jaccard",), "higher", minimum=0, maximum=1, engines=("ipsae",), stages=("analysis",)),
    Measure("interface.analysis.ipsae_structure.binder_residue_preservation_fraction", "Binder interface residue preservation fraction", "fraction", "designed binder interface residues", "Fraction of designed binder interface residues that retain any target contact in cofold",
        ("structure_binder_interface_residue_preservation_fraction",), "higher", minimum=0, maximum=1, engines=("ipsae",), stages=("analysis",)),
    Measure("interface.analysis.ipsae_structure.target_residue_preservation_fraction", "Target interface residue preservation fraction", "fraction", "designed target interface residues", "Fraction of designed target interface residues that retain any binder contact in cofold",
        ("structure_target_interface_residue_preservation_fraction",), "higher", minimum=0, maximum=1, engines=("ipsae",), stages=("analysis",)),
    Measure("interface.analysis.ipsae_structure.sequence_mapping_coverage", "Designed/cofold sequence mapping coverage", "fraction", "designed binder and target residues", "Minimum aligned-residue coverage across binder and target sequence alignments used for contact preservation",
        ("structure_interaction_sequence_mapping_coverage",), "higher", minimum=0, maximum=1, engines=("ipsae",), stages=("analysis",)),
    Measure("geometry.analysis.ipsae_structure.ca_clashes_2_5a", "ipSAE analysis · inter-chain Cα pairs ≤2.5 Å", "pairs", "selected binder and target chains", "Biopython NeighborSearch inter-chain C-alpha pairs within 2.5 Å",
        ("structure_ca_clash_pair_count_2_5a",), "lower", minimum=0, engines=("ipsae",), stages=("analysis",)),
    Measure("confidence.analysis.structure.binder_ca_bfactor_mean", "ipSAE analysis · binder mean Cα B-factor", "B-factor units", "selected binder chains", "Arithmetic mean of saved structure C-alpha B-factor fields; not assumed to be pLDDT",
        ("structure_binder_ca_bfactor_mean",), engines=("ipsae",), stages=("analysis",)),
    Measure("confidence.analysis.structure.target_ca_bfactor_mean", "ipSAE analysis · target mean Cα B-factor", "B-factor units", "selected target chains", "Arithmetic mean of saved structure C-alpha B-factor fields; not assumed to be pLDDT",
        ("structure_target_ca_bfactor_mean",), engines=("ipsae",), stages=("analysis",)),
    Measure("confidence.analysis.structure.complex_ca_bfactor_mean", "ipSAE analysis · complex mean Cα B-factor", "B-factor units", "selected binder and target chains", "Arithmetic mean of binder and target chain C-alpha B-factor means; not assumed to be pLDDT",
        ("structure_complex_ca_bfactor_mean",), engines=("ipsae",), stages=("analysis",)),
    Measure("geometry.analysis.binder.design_to_cofold_fitted_ca_rmsd", "Binder design-to-cofold fitted Cα RMSD", "Å", "sequence-aligned binder C-alpha atoms", "Design and predicted binder matched by global sequence alignment; optimal binder-only superposition. Does not measure target-relative pose RMSD",
        ("structure_binder_design_to_cofold_fitted_ca_rmsd",), "lower", minimum=0, engines=("ipsae",), stages=("analysis",)),
    Measure("interface.analysis.hotspots.contact_fraction_found", "ipSAE analysis · hotspot contact fraction (found residues)", "fraction", "configured target hotspots found in predicted structure", "Fraction of located configured hotspots within the configured heavy-atom distance cutoff",
        ("structure_hotspot_contact_fraction_found",), "higher", minimum=0, maximum=1, engines=("ipsae",), stages=("analysis",)),
    Measure("interface.analysis.hotspots.contact_fraction_requested", "ipSAE analysis · hotspot contact fraction (requested residues)", "fraction", "all configured target hotspots", "Contacted configured hotspots divided by requested hotspots; missing target residues do not count as contacts",
        ("structure_hotspot_contact_fraction_requested",), "higher", minimum=0, maximum=1, engines=("ipsae",), stages=("analysis",)),
    Measure("interface.analysis.hotspots.contacted_count", "ipSAE analysis · contacted hotspot count", "residues", "configured target hotspots", "Configured target hotspots with a heavy-atom contact under the configured cutoff",
        ("structure_hotspot_contacted_count",), "higher", minimum=0, engines=("ipsae",), stages=("analysis",)),
    Measure("interface.analysis.hotspots.minimum_distance", "ipSAE analysis · nearest hotspot distance", "Å", "configured target hotspots and binder heavy atoms", "Minimum heavy-atom distance from any located hotspot residue to selected binder chains",
        ("structure_hotspot_min_distance",), "lower", minimum=0, engines=("ipsae",), stages=("analysis",)),
    Measure("confidence.native.bindcraft2.interface.ipdae", "BindCraft2 i_pDAE", "0–1", "contact-focused binder–target interface", "BindCraft2 PAE-derived score over Cα contacts within 8 Å; maximum over contacting residues and directions",
        ("i_pDAE", "ipdae"), "higher", minimum=0, maximum=1, engines=("bindcraft2",)),
    Measure("confidence.native.bindcraft2.interface.ipae_normalized", "BindCraft2 normalized i_pAE", "0–1", "binder-to-target directional mean", "BindCraft2 mean interface PAE normalized by 31 Å; distinct from raw PAE in Å",
        ("i_pAE",), "lower", minimum=0, maximum=1, engines=("bindcraft2",)),
    Measure("geometry.native.bindcraft2.bound_unbound_rmsd", "BindCraft2 Binder_RMSD", "Å", "binder backbone between bound and unbound predictions", "BindCraft2 matched binder Cα RMSD; distinct from target-aligned pose RMSD",
        ("Binder_RMSD",), "lower", minimum=0, engines=("bindcraft2",)),
    Measure("sequence.native.bindcraft2.scaffold_sequence_retained_fraction", "BindCraft2 scaffold sequence retained fraction", "fraction", "held scaffold framework", "BindCraft2 exact sequence identity at non-edited scaffold positions",
        ("Scaffold_Sequence_Retained_Fraction",), "higher", minimum=0, maximum=1, engines=("bindcraft2",)),
    Measure("geometry.native.bindcraft2.scaffold_framework_rmsd", "BindCraft2 scaffold framework RMSD", "Å", "held scaffold framework", "BindCraft2 framework Cα RMSD after scaffold alignment",
        ("Scaffold_Framework_RMSD",), "lower", minimum=0, engines=("bindcraft2",)),
    Measure("interface.native.bindcraft2.off_paratope_contact_fraction", "BindCraft2 off-paratope contact fraction", "fraction", "binder residues contacting outside the paratope", "BindCraft2 flagged-paratope contact residue fraction",
        ("Off_Paratope_Contact_Fraction",), "lower", minimum=0, maximum=1, engines=("bindcraft2",)),
    Measure("structure.native.bindcraft2.framework_packing_fraction", "BindCraft2 framework packing fraction", "fraction", "framework residues packed against paratope", "BindCraft2 framework-to-paratope heavy-atom packing fraction",
        ("Framework_Packing_Fraction",), "higher", minimum=0, maximum=1, engines=("bindcraft2",)),
    Measure("geometry.native.bindcraft2.induced_fit_interface_rmsd", "BindCraft2 induced-fit interface RMSD", "Å", "interface residues between free and bound states", "BindCraft2 interface Cα displacement after core alignment",
        ("Induced_Fit_Interface_RMSD",), "higher", minimum=0, engines=("bindcraft2",)),
    Measure("confidence.native.bindcraft2.induced_fit_tm", "BindCraft2 induced-fit TM score", "0–1", "whole binder between free and bound states", "BindCraft2 structural TM score",
        ("Induced_Fit_TM",), "lower", minimum=0, maximum=1, engines=("bindcraft2",)),
    Measure("geometry.native.bindcraft2.interdomain_contact_fraction", "BindCraft2 interdomain contact fraction", "fraction", "binder domain interfaces", "BindCraft2 interdomain contacts divided by total contacts",
        ("Interdomain_Contact_Fraction",), "lower", minimum=0, maximum=1, engines=("bindcraft2",)),
    Measure("geometry.native.bindcraft2.domain_separation_ratio", "BindCraft2 domain separation ratio", "ratio", "binder domains", "BindCraft2 domain separation relative to domain size",
        ("Domain_Separation_Ratio",), "higher", minimum=0, engines=("bindcraft2",)),
    Measure("interface.native.bindcraft2.binder_residue_count_4a", "BindCraft2 binder interface residue count", "residues", "binder–target interface", "BindCraft2 binder residues with any atom within 4 Å of the target",
        ("Interface_Residues",), "higher", minimum=0, engines=("bindcraft2",)),
    Measure("geometry.native.bindcraft2.backbone_ca_clashes_2_5a", "BindCraft2 backbone Cα clash count", "atom pairs", "multichain structure", "BindCraft2 interchain Cα pairs within 2.5 Å",
        ("Backbone_Clashes",), "lower", minimum=0, engines=("bindcraft2",)),
])
for engine, method, rmsd_key, rmsd_label, rmsd_alias in (
    ("pxdesign", "PXDesign native AF2-IG", "geometry.native.pxdesign.bound_unbound_rmsd", "PXDesign bound–unbound binder RMSD", "binder_rmsd"),
    ("genie3", "Genie3 native evaluation", "geometry.native.genie3.binder_self_consistency_rmsd", "Genie3 binder self-consistency RMSD", "binder_rmsd"),
    ("protpardelle", "Protpardelle native ESMFold evaluation", "geometry.native.protpardelle.binder_scaffold_rmsd", "Protpardelle binder scaffold RMSD", "binder_rmsd"),
    ("proteina_complexa", "Proteina-Complexa native evaluation", "geometry.native.proteina_complexa.binder_scaffold_rmsd", "Proteina-Complexa binder scaffold RMSD", "binder_rmsd"),
):
    for key, label, aliases, unit, direction, maximum in (
        ("confidence.validation.bound_binder.plddt", "Validation bound binder pLDDT", ("binder_plddt",), "0–100", "higher", 100),
        ("confidence.validation.interface.mean_pae", "Validation interface PAE", ("ipae",), "Å", "lower", None),
    ):
        MEASURES.append(Measure(key, label, unit, "binder in complex" if key.endswith("plddt") else "binder–target interface",
            method, aliases, direction, minimum=0, maximum=maximum, engines=(engine,), stages=("complex_refolding",)))
    MEASURES.append(Measure(rmsd_key, rmsd_label, "Å", "native evaluator binder RMSD", method,
        (rmsd_alias,), "lower", minimum=0, engines=(engine,), stages=("complex_refolding",)))
for key, label, aliases in (
    ("interface.ipsae.minimum", "Minimum directional ipSAE", ("ipsae_min",)),
    ("interface.ipsae.maximum", "Maximum directional ipSAE", ("ipsae_max",)),
    ("interface.ipsae.mean", "Mean directional ipSAE", ("ipsae_avg",)),
):
    MEASURES.append(Measure(key, label, "0–1", "native selected chain pair", "Proteina-Complexa native evaluation",
        aliases, "higher", minimum=0, maximum=1, engines=("proteina_complexa",), stages=("complex_refolding",)))
for key, label, unit, scope, method, aliases, direction, minimum, maximum in (
    ("interface.ipsae.minimum", "Minimum directional ipSAE", "0–1", "selected chain pair", "ipSAE", ("ipsae_min", "ipSAE_min"), "higher", 0, 1),
    ("interface.ipsae.maximum", "Maximum directional ipSAE", "0–1", "selected chain pair", "ipSAE", ("ipsae_max", "ipSAE_max"), "higher", 0, 1),
    ("interface.ipsae.mean", "Mean directional ipSAE", "0–1", "selected chain pair", "ipSAE", ("ipsae_avg", "ipSAE_avg"), "higher", 0, 1),
    ("energy.rosetta.interface_dg", "Rosetta interface energy", "REU", "binder–target interface", "PyRosetta InterfaceAnalyzer", ("rosetta_interface_dG",), "lower", None, None),
    ("interface.rosetta.buried_area", "Rosetta buried interface area", "Å²", "binder–target interface", "PyRosetta InterfaceAnalyzer", ("rosetta_interface_dSASA",), "", 0, None),
    ("developability.rosetta.sap", "Spatial aggregation propensity", "Rosetta SAP", "binder", "PyRosetta SapScoreMetric", ("rosetta_binder_sap",), "lower", None, None),
    ("sequence.charge", "Sequence charge", "e", "selected protein chains", "Biopython ProteinAnalysis", ("qc_sequence_charge",), "", None, None),
    ("sequence.isoelectric_point", "Isoelectric point", "pH", "selected protein chains", "Biopython ProteinAnalysis", ("qc_sequence_pi",), "", 0, 14),
    ("sequence.entropy", "Composition entropy", "bits", "selected protein chains", "Shannon amino-acid composition", ("qc_sequence_entropy",), "", 0, None),
    ("sequence.gravy", "Mean hydropathy", "Kyte–Doolittle", "selected protein chains", "Biopython ProteinAnalysis", ("qc_sequence_gravy",), "", None, None),
    ("sequence.proteinsol", "Protein-Sol scaled solubility", "Protein-Sol", "selected protein chains", "Protein-Sol", ("proteinsol_scaled_solubility",), "higher", None, None),
    ("interface.contacts", "Interface heavy-atom contact pairs", "pairs", "selected chain groups", "Biotite distance contacts", ("qc_interface_contact_pairs",), "", 0, None),
    ("interface.residue_contacts", "Interface residue contact pairs", "pairs", "selected chain groups", "Biotite distance contacts", ("qc_interface_residue_pairs",), "", 0, None),
    ("interface.hydrogen_bonds", "Interface hydrogen bonds", "bonds", "selected chain groups", "Biotite geometric hydrogen bonds", ("qc_interface_hydrogen_bonds",), "", 0, None),
    ("diversity.cluster", "Diversity cluster", "cluster ID", "comparison cohort", "declared clustering method", ("diversity_cluster",), "", 0, None),
):
    MEASURES.append(Measure(key, label, unit, scope, method, aliases, direction, minimum, maximum,
        description=f"{label}. Scope: {scope}. Method and settings remain part of the producing evaluation."))
MEASURES.append(Measure("surface.peppatch.area", "PEP-Patch surface area", "Å²", "selected protein chains", "PEP-Patch 65da35b", ("peppatch_total_area",),
    prerequisites=("protein structure", "complete heavy atoms")))
MEASURES.extend([
    Measure("sequence.esm1v.native_probability", "ESM-1v native amino-acid probability", "fraction", "selected sequence", "ESM-1v five-model unmasked wild-type marginals", ("esm1v_native_probability",), "higher", maximum=1),
    Measure("sequence.esmif.native_log_likelihood", "ESM-IF1 native conditional log likelihood", "nats/residue", "selected chains in complex backbone", "ESM-IF1", ("esmif_native_log_likelihood",), "higher", minimum=None, maximum=0),
    Measure("ligand.posebusters.passed", "All PoseBusters checks passed", "0/1", "selected ligand pose and receptor", "PoseBusters 1a5f26a", ("posebusters_passed_all",), "higher", maximum=1, transform="bool_to_indicator"),
    Measure("ligand.plip.interactions", "PLIP interaction count", "interactions", "selected ligand pose and receptor", "PLIP", ("plip_interaction_count",)),
    Measure("ligand.pandamap.interactions", "PandaMap interaction count", "interactions", "selected ligand pose and receptor", "PandaMap 0007347", ("pandamap_interaction_count",)),
])
DEFINITIONS = {measure.key: measure for measure in reversed(MEASURES)}


def field_label(field):
    measure = DEFINITIONS.get(str(field).removeprefix(PREFIX)) if str(field).startswith(PREFIX) else None
    return f"{measure.label} ({measure.unit}) · Common" if measure else field


def vocabulary():
    """Machine-readable dictionary of IDs, definitions and conversion rules."""
    return {"version": VERSION, "field_prefix": PREFIX, "engine_aliases": dict(ENGINE_ALIASES),
            "generator_workflows": dict(GENERATOR_WORKFLOWS),
            "prediction_measurements": prediction_engine_catalog(),
            "measures": [{"field": PREFIX + key, "key": key,
                          **{name: getattr(definition, name) for name in ("label", "unit", "scope", "direction", "minimum", "maximum", "description", "value_type", "prerequisites", "aggregation")},
                          "range_semantics": "Physical/value-domain validation only; these limits are not quality thresholds.",
                          "group": key.split(".", 1)[0],
                          "rules": [asdict(rule) for rule in MEASURES if rule.key == key]}
                         for key, definition in DEFINITIONS.items()]}


def measure_capabilities(tool, *, stages=None):
    """Return only engine-scoped measurements whose exact source aliases match.

    This is a description of output normalization, not a claim that a given
    run or candidate produced a value. Generic analysis measurements are
    advertised by the analysis route that computes them.
    """
    engine = engine_name(tool)
    stage_set = set(stages or ())
    rows = []
    for measure in MEASURES:
        if not measure.engines or engine not in measure.engines:
            continue
        if stage_set and measure.stages and not stage_set.intersection(measure.stages):
            continue
        rows.append({
            "field": PREFIX + measure.key,
            "key": measure.key,
            "label": measure.label,
            "unit": measure.unit,
            "scope": measure.scope,
            "method": measure.method,
            "direction": measure.direction,
            "recognized_aliases": list(measure.aliases),
            "stages": list(measure.stages),
            "prerequisites": list(measure.prerequisites),
            "value_domain": {"minimum": measure.minimum, "maximum": measure.maximum},
            "quality_threshold": "not set; configure explicitly in native filters or Shared scientific filters",
            "availability": "recognized when the source field is emitted; missing values remain missing",
        })
    return rows


def normalize_metrics(metrics, tool="", context=None):
    """Return common values, alias provenance and conflicts, without editing raw metrics."""
    context = context if isinstance(context, dict) else {}
    candidate = context.get("candidate") if isinstance(context.get("candidate"), dict) else context
    raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
    tool = candidate.get("source_tool") or tool
    engine = engine_name(tool)
    stage = _candidate_stage(candidate)
    original = metrics or {}
    values = recovered_metrics(original, raw, tool)
    sequence = candidate.get("binder_sequence")
    if sequence and candidate.get("stage") != "generation.backbone":
        chains = re.split(r"[:/]", str(sequence).replace(" ", "").replace("\n", ""))
        if chains and all(re.fullmatch(r"[A-Za-z]+", chain) for chain in chains):
            values["$sequence_length"] = sum(map(len, chains))
    matches = {}
    issues = []
    for measure in MEASURES:
        if measure.engines and engine not in measure.engines:
            continue
        if measure.stages and stage not in measure.stages:
            continue
        if any(values.get(key, candidate.get(key)) != expected for key, expected in measure.requires):
            continue
        for alias in measure.aliases:
            value = values.get(alias)
            if value is None:
                continue
            if measure.value_type == "text":
                if not isinstance(value, str) or not value:
                    issues.append({"field": alias, "reason": "not a non-empty string"})
                    continue
                matches.setdefault(measure.key, []).append({"source_field": alias, "raw_value": value, "value": value, "conversion": "identity", "method": measure.method, "_engine_scoped": bool(measure.engines)})
                continue
            if measure.transform == "bool_to_indicator" and isinstance(value, bool):
                value = float(value)
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                issues.append({"field": alias, "reason": "not a finite number"})
                continue
            try:
                number = float(value)
            except OverflowError:
                number = math.inf
            if not math.isfinite(number):
                issues.append({"field": alias, "reason": "not a finite number"})
                continue
            if measure.transform == "negative_log":
                if not 0 < number <= 1:
                    issues.append({"field": alias, "reason": "probability outside (0, 1]"})
                    continue
                number = -math.log(number)
            elif measure.transform == "fraction_to_percent":
                number *= 100
            elif measure.transform == "percent_to_fraction":
                number /= 100
            elif measure.transform == "times_31":
                number *= 31
            elif measure.transform == "divide_31":
                number /= 31
            if not math.isfinite(number):
                issues.append({"field": alias, "reason": "conversion is not finite"})
                continue
            if (measure.minimum is not None and number < measure.minimum) or (measure.maximum is not None and number > measure.maximum):
                issues.append({"field": alias, "reason": "outside declared range"})
                continue
            if alias == "$sequence_length":
                source = "binder_sequence"
            elif alias in original:
                source = alias
            elif engine == "mol1" and alias == "mol1_model_plddt_mean":
                confidence = raw.get("confidence_outputs") if isinstance(raw.get("confidence_outputs"), dict) else {}
                confidence_field = "plddt_mean" if confidence.get("plddt_mean") is not None else "plddt"
                source = f"raw_metadata.confidence_outputs.{confidence_field}"
            else:
                source = "fasta_header." + alias
            matches.setdefault(measure.key, []).append({"source_field": source, "raw_value": value, "value": number, "conversion": measure.transform, "method": measure.method, "_engine_scoped": bool(measure.engines)})
    normalized, provenance, conflicts = {}, {}, {}
    for key, entries in matches.items():
        # Native FASTA values are often rounded to four decimal places.
        tolerance = 1e-3 if key.endswith("redesigned_nll") else 1e-8
        if any(not math.isclose(entries[0]["value"], entry["value"], rel_tol=tolerance, abs_tol=tolerance) for entry in entries[1:]):
            conflicts[key] = entries
            continue
        definition = DEFINITIONS[key]
        preferred = next((entry for entry in entries if entry["_engine_scoped"]), entries[0])
        public_entries = [{name: value for name, value in entry.items() if name != "_engine_scoped"} for entry in entries]
        normalized[key] = preferred["value"]
        provenance[key] = {"unit": definition.unit, "scope": definition.scope, "method": preferred["method"],
                           "direction": definition.direction, "source_tool": tool, "aliases": public_entries,
                           "stage": stage, "engine_version": raw.get("engine_version"),
                           "conditioning_mode": raw.get("conditioning_mode"),
                           "settings": context.get("settings") or raw.get("prediction_protocol") or {}}
    return {"version": VERSION, "metrics": normalized, "provenance": provenance, "conflicts": conflicts, "issues": issues}


def candidate_output(candidate):
    result = normalize_metrics(candidate.get("metrics") or {}, candidate.get("source_tool"), candidate)
    stage = _candidate_stage(candidate)
    result.update(engine=engine_name(candidate.get("source_tool")), stage=stage,
                  kind={"generation.backbone": "backbone_proposal", "generation.backbone_sequence": "designed_structure", "sequence_design": "sequence_design", "complex_refolding": "refolded_complex", "monomer_refolding": "refolded_monomer"}.get(stage, stage),
                  designed_sequence=candidate.get("binder_sequence") if stage != "generation.backbone" else None,
                  designed_chains=list(candidate.get("binder_chains") or []), target_chains=list(candidate.get("target_chains") or []))
    result["capabilities"] = {"backbone": bool((candidate.get("complex_pdb") or candidate.get("binder_pdb")) and candidate.get("binder_chains")),
                              "sequence": bool(result["designed_sequence"]), "complex_context": bool(candidate.get("target_pdb") or candidate.get("target_chains"))}
    result["parents"] = list(candidate.get("parents") or [])
    return result


def metric_value(candidate, field):
    """Shared selector for workflow filters; native metric keys still work."""
    if str(field).startswith("evidence/"):
        parts = str(field).split("/")
        if len(parts) not in {4, 5}:
            raise ValueError("Use evidence/<stage>/<tool>/<metric>, optionally with a candidate ID before the metric.")
        _, stage, tool, *rest = parts
        matches = [record for record in candidate_evidence(candidate)
                   if record["stage"] == stage and record["source_tool"] == tool
                   and record["same_sequence"] and (len(rest) == 1 or record["candidate_id"] == rest[0])]
        if len(matches) != 1:  # ambiguity is missing evidence, never an implicit latest/best score
            return None
        record = matches[0]
        if record["status"] in {"error", "failed", "missing", "invalid"}:
            return None
        return metric_value(record["candidate"], rest[-1])
    if str(field).startswith(PREFIX):
        return candidate_output(candidate)["metrics"].get(str(field)[len(PREFIX):])
    return (candidate.get("metrics") or {}).get(field)


def candidate_evidence(candidate):
    """Read lineage without relabeling ancestral measurements as current ones.

    Different sequences remain visible as comparators, but cannot supply a
    measurement for the current sequence. No files are loaded by this view.
    """
    records, seen = [], set()
    pending = [(candidate, 0, True)]
    while pending:
        item, depth, compatible = pending.pop()
        if id(item) in seen:
            continue
        if depth > 64:
            raise ValueError("Candidate evidence lineage exceeds 64 stages.")
        seen.add(id(item))
        raw, metrics = item.get("raw_metadata") or {}, item.get("metrics") or {}
        outcome = raw.get("measurement_status") or {}
        errors = [str(metrics[key]) for key in ("ipsae_error", "rosetta_error", "qc_error", "validation_error") if metrics.get(key)]
        sequence = item.get("binder_sequence")
        compatible = compatible and (depth == 0 or bool(sequence and sequence == candidate.get("binder_sequence")))
        records.append({"candidate_id": item.get("candidate_id"), "stage": item.get("stage"),
                        "source_tool": item.get("source_tool") or item.get("tool"),
                        "parents": list(item.get("parents") or []), "depth": depth,
                        "same_sequence": compatible, "role": "candidate" if depth == 0 else "evidence" if compatible else "comparator",
                        "status": "error" if errors else outcome.get("status") or metrics.get("analysis_status") or "unverified",
                        "error": "; ".join(errors) or outcome.get("error", ""), "candidate": item})
        parent = raw.get("source_candidate")
        if isinstance(parent, dict):
            pending.append((parent, depth + 1, compatible))
        additional = raw.get("evidence_candidates")
        if isinstance(additional, list):
            pending.extend((item, depth + 1, compatible) for item in reversed(additional) if isinstance(item, dict))
    return records


def evidence_rows(candidates):
    """Long-form Results projection; sources, units and protocols stay explicit."""
    rows = []
    for candidate in candidates:
        for record in candidate_evidence(candidate):
            item = record["candidate"]
            common = candidate_output(item)
            fields = {**(item.get("metrics") or {}), **{PREFIX + key: value for key, value in common["metrics"].items()}}
            for field, value in fields.items():
                provenance = common.get("provenance", {}).get(field.removeprefix(PREFIX), {}) if field.startswith(PREFIX) else {}
                rows.append({"candidate_id": candidate.get("candidate_id"), "measurement_candidate_id": record["candidate_id"],
                             "stage": record["stage"], "source_tool": record["source_tool"], "role": record["role"],
                             "conditioning_mode": (item.get("raw_metadata") or {}).get("conditioning_mode", ""),
                             "effective_settings": (item.get("raw_metadata") or {}).get("effective_settings", {}),
                             "same_sequence": record["same_sequence"], "measurement_status": record["status"], "error": record["error"],
                             "field": field, "value": value, "unit": provenance.get("unit", "native"),
                             "direction": provenance.get("direction", ""), "parents": ", ".join(record["parents"]),
                             "selector": f"evidence/{record['stage']}/{record['source_tool']}/{record['candidate_id']}/{field}"})
    return rows


@lru_cache(maxsize=2048)
def sql_common_metrics(metrics_json, namespace, context_json):
    """Bounded projection cache; reads old and synchronized observations equally."""
    return json.dumps(normalize_metrics(json.loads(metrics_json), namespace, json.loads(context_json))["metrics"], allow_nan=False)
