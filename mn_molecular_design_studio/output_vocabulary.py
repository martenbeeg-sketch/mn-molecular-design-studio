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
    return ENGINE_ALIASES.get(name, name)


def fasta_metrics(header):
    """Read exact FASTA header keys, including the current LigandMPNN format."""
    return {key: float(value) for key, value in _HEADER_NUMBER.findall(str(header or ""))
            if key in FASTA_METRICS and math.isfinite(float(value))}


def recovered_metrics(metrics, raw, tool):
    result = dict(metrics or {})
    if engine_name(tool) in MPNN_ENGINES:
        for key, value in fasta_metrics((raw or {}).get("fasta_header")).items():
            result.setdefault(key, value)
    return result


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
    ):
        for canonical, alias in aliases:
            MEASURES.append(Measure(f"confidence.{role}.complex.{canonical}", f"{role.title()} complex {canonical}", "0–1", "complex" if canonical == "ptm" else "complex interfaces", method, (alias,), "higher", maximum=1, engines=(engine,), stages=stages))
    for key, label, scope, alias in (
        ("bound_binder.plddt", "bound binder pLDDT", "binder in complex", "pLDDT"),
        ("target.plddt", "target pLDDT", "target in complex", "Target_pLDDT"),
        ("unbound_binder.plddt", "unbound binder pLDDT", "binder alone", "Unbound_Binder_pLDDT"),
    ):
        MEASURES.append(Measure(f"confidence.{role}.{key}", f"{role.title()} {label}", "0–100", scope, "BindCraft2 AF2", (alias,), "higher", maximum=100, engines=("bindcraft2",), transform="fraction_to_percent", stages=stages))
    MEASURES.append(Measure(f"confidence.{role}.interface.mean_pae", f"{role.title()} binder–target mean PAE", "Å", "binder–target symmetric mean", "BindCraft2 AF2", ("i_pAE",), "lower", engines=("bindcraft2",), transform="times_31", stages=stages))
# RFD4's oracle flattens AF3 confidence outputs after scaling pLDDT/pTM by
# 100; PAE remains in Å. These aliases are scoped to RFD4's generator stage,
# since trajectory confidence is not an independent validation run.
for scope, suffix, aliases in (
    ("complex", "", ("oracle.overall_plddt", "oracle.overall_pae", "oracle.overall_ptm")),
    ("chain of interest", ".chain_of_interest", ("oracle.plddt_coi", "oracle.pae_coi", "oracle.ptm_coi")),
):
    for quantity, unit, direction, transform, maximum, alias in zip(
        ("plddt", "pae", "ptm"), ("0–100", "Å", "0–1"),
        ("higher", "lower", "higher"), ("identity", "identity", "percent_to_fraction"),
        (100, None, 1), aliases,
    ):
        MEASURES.append(Measure(
            f"confidence.generator.{scope.replace(' ', '_')}.{quantity}",
            f"Generator {scope} {quantity}", unit, scope,
            "RFD4 AF3 oracle", (alias,), direction, maximum=maximum,
            engines=("rfd4_proteina",), transform=transform,
            stages=("generation.backbone_sequence",),
        ))
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
# Evaluation descriptors are defined here as well, so ranking and UI use the
# same conversions and cannot substitute confidence for an unrelated measure.
for engine in ("af2_initial_guess", "boltz2_initial_guess", "boltz2", "af3_jax_cofolding", "optimized_prediction", "openfold3", "protenix", "protenix_v1"):
    for quantity in ("ptm", "iptm"):
        MEASURES.append(Measure(f"confidence.validation.complex.{quantity}", f"Validation complex {quantity}", "0–1", "complex interfaces" if quantity == "iptm" else "complex", engine,
            (quantity,), "higher", maximum=1, engines=(engine,), stages=("complex_refolding",)))
MEASURES.append(Measure("confidence.validation.bound_binder.plddt", "Validation bound binder pLDDT", "0–100", "binder in complex", "AF2 binder confidence",
    ("af2_binder_plddt", "af2_average_binder_plddt"), "higher", maximum=100))
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
            "measures": [{"field": PREFIX + key, "key": key,
                          **{name: getattr(definition, name) for name in ("label", "unit", "scope", "direction", "minimum", "maximum", "description", "value_type", "prerequisites", "aggregation")},
                          "group": key.split(".", 1)[0],
                          "rules": [asdict(rule) for rule in MEASURES if rule.key == key]}
                         for key, definition in DEFINITIONS.items()]}


def normalize_metrics(metrics, tool="", context=None):
    """Return common values, alias provenance and conflicts, without editing raw metrics."""
    context = context if isinstance(context, dict) else {}
    candidate = context.get("candidate") if isinstance(context.get("candidate"), dict) else context
    raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
    tool = candidate.get("source_tool") or tool
    engine = engine_name(tool)
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
        if measure.stages and candidate.get("stage") not in measure.stages:
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
                matches.setdefault(measure.key, []).append({"source_field": alias, "raw_value": value, "value": value, "conversion": "identity", "method": measure.method})
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
            if not math.isfinite(number):
                issues.append({"field": alias, "reason": "conversion is not finite"})
                continue
            if (measure.minimum is not None and number < measure.minimum) or (measure.maximum is not None and number > measure.maximum):
                issues.append({"field": alias, "reason": "outside declared range"})
                continue
            source = "binder_sequence" if alias == "$sequence_length" else alias if alias in original else "fasta_header." + alias
            matches.setdefault(measure.key, []).append({"source_field": source, "raw_value": value, "value": number, "conversion": measure.transform, "method": measure.method})
    normalized, provenance, conflicts = {}, {}, {}
    for key, entries in matches.items():
        # Native FASTA values are often rounded to four decimal places.
        tolerance = 1e-3 if key.endswith("redesigned_nll") else 1e-8
        if any(not math.isclose(entries[0]["value"], entry["value"], rel_tol=tolerance, abs_tol=tolerance) for entry in entries[1:]):
            conflicts[key] = entries
            continue
        definition = DEFINITIONS[key]
        normalized[key] = entries[0]["value"]
        provenance[key] = {"unit": definition.unit, "scope": definition.scope, "method": entries[0]["method"],
                           "direction": definition.direction, "source_tool": tool, "aliases": entries,
                           "stage": candidate.get("stage"), "engine_version": raw.get("engine_version"),
                           "settings": context.get("settings") or {}}
    return {"version": VERSION, "metrics": normalized, "provenance": provenance, "conflicts": conflicts, "issues": issues}


def candidate_output(candidate):
    result = normalize_metrics(candidate.get("metrics") or {}, candidate.get("source_tool"), candidate)
    stage = candidate.get("stage")
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
    if str(field).startswith(PREFIX):
        return candidate_output(candidate)["metrics"].get(str(field)[len(PREFIX):])
    return (candidate.get("metrics") or {}).get(field)


@lru_cache(maxsize=2048)
def sql_common_metrics(metrics_json, namespace, context_json):
    """Bounded projection cache; reads old and synchronized observations equally."""
    return json.dumps(normalize_metrics(json.loads(metrics_json), namespace, json.loads(context_json))["metrics"], allow_nan=False)
