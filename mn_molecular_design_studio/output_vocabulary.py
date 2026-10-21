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
    "design.genie3_generate": "genie3", "design.protpardelle_generate": "protpardelle",
    "design.proteina_complexa": "proteina_complexa", "design.proteina_complexa_generate": "proteina_complexa",
    "design.rfantibody": "rfantibody", "design.rfpeptides": "rfpeptides",
    "design.esmfold2_binder": "esmfold2_binder",
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
for model in sorted(MPNN_ENGINES):
    for key, label, aliases, unit, direction, transform, maximum in (
        ("redesigned_nll", "redesigned-residue NLL", ("score",), "nats/residue", "lower", "identity", None),
        ("redesigned_nll", "redesigned-residue NLL", ("overall_confidence",), "nats/residue", "lower", "negative_log", None),
        ("global_nll", "global NLL", ("global_score",), "nats/residue", "lower", "identity", None),
        ("redesigned_probability", "redesigned-residue probability", ("overall_confidence",), "fraction", "higher", "identity", 1),
    ):
        MEASURES.append(Measure(f"sequence.{model}.{key}", f"{model} {label}", unit, "global" if key == "global_nll" else "redesigned residues near ligand" if key == "ligand_probability" else "redesigned residues", model, aliases, direction, maximum=maximum, engines=(model,), transform=transform))
MEASURES.extend([
    Measure("sequence.recovery_designed", "Sequence recovery over redesigned residues", "fraction", "redesigned residues", "identity to input sequence", ("seq_recovery", "seq_rec"), maximum=1, engines=tuple(sorted(MPNN_ENGINES))),
    Measure("sequence.esm_if1.log_likelihood", "ESM-IF1 log likelihood", "native", "designed sequence", "ESM-IF1", ("esm_if1_log_likelihood",), "higher", minimum=None, engines=("esm_if1",)),
    Measure("sequence.caliby.energy", "Caliby energy", "native", "designed sequence", "Caliby", ("caliby_energy",), "lower", minimum=None, engines=("caliby",)),
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
                          **{name: getattr(definition, name) for name in ("label", "unit", "scope", "direction", "minimum", "maximum")},
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
        if any(values.get(key) != expected for key, expected in measure.requires):
            continue
        for alias in measure.aliases:
            value = values.get(alias)
            if value is None:
                continue
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
