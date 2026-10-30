"""Shared validation and acceptance policy for native and composed workflows.

Native tools keep their own decisions. This module evaluates only explicitly
configured criteria against supplied measurements; it does not infer missing
scores, convert one metric into another, or run a prediction engine.
"""
from __future__ import annotations

import copy
import csv
import json
import math
import operator
from pathlib import Path

from mn_molecular_design_studio.measurement_status import RANGES, finite_number
from mn_molecular_design_studio.output_vocabulary import DEFINITIONS, PREFIX, candidate_output, metric_value
from mn_molecular_design_studio.prediction_measurements import PAE_MATRIX_ENGINES, SCALAR_SCORE_ALIASES

COMPARISONS = {">=": operator.ge, "<=": operator.le, ">": operator.gt, "<": operator.lt}
BINDER_FORMATS = {
    "miniprotein": "Miniprotein",
    "large_protein": "Large protein",
    "peptide": "Linear peptide",
    "cyclic_peptide": "Cyclic peptide",
    "vhh": "VHH",
    "scfv": "scFv",
    "fab": "Fab",
    "arp": "ARP",
    "homo_oligomer": "Homo-oligomer",
    "multidomain": "Multidomain binder",
    "induced_fit": "Induced-fit binder",
    "fold_switch": "Fold-switch binder",
    "other": "Other protein binder",
}
TARGET_TYPES = {"structured": "Structured target", "fasta_idr": "FASTA / IDR target"}
FORMAT_ALIASES = {
    "miniprotein": "miniprotein", "binder": "miniprotein", "protein_binder": "miniprotein",
    "de novo binder": "miniprotein", "large_binder": "large_protein", "large protein": "large_protein",
    "peptide": "peptide", "cyclic_peptide": "cyclic_peptide", "cyclic peptide": "cyclic_peptide", "cyclotide": "cyclic_peptide",
    "vhh": "vhh", "nanobody": "vhh", "scfv": "scfv", "fab": "fab", "arp": "arp",
    "homo_oligomer": "homo_oligomer", "multidomain": "multidomain",
    "induced_fit": "induced_fit", "fold_switch": "fold_switch", "fold_switching": "fold_switch",
}
TARGET_ALIASES = {
    "sequence": "fasta_idr", "fasta": "fasta_idr", "fasta_idr": "fasta_idr", "idr": "fasta_idr",
    "disordered": "fasta_idr", "unstructured": "fasta_idr",
    "structure": "structured", "structured": "structured", "pdb": "structured", "cif": "structured",
}


def infer_binder_format(candidate: dict, workflow_parameters: dict | None = None) -> str | None:
    """Infer a format only from explicit design provenance, never sequence length."""
    raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
    metrics = candidate.get("metrics") if isinstance(candidate.get("metrics"), dict) else {}

    def values_for(sources, keys):
        values = []
        for source in sources:
            for key in keys:
                value = source.get(key) if isinstance(source, dict) else None
                values.extend(value if isinstance(value, (list, tuple)) else [value] if value is not None else [])
        return {FORMAT_ALIASES[str(value).strip().lower()] for value in values
                if str(value).strip().lower() in FORMAT_ALIASES}

    def nested_values(source, keys):
        values = set()
        pending = [source] if isinstance(source, dict) else []
        while pending:
            item = pending.pop()
            if isinstance(item, (list, tuple)):
                pending.extend(value for value in item if isinstance(value, (dict, list, tuple)))
                continue
            for key in keys:
                value = item.get(key)
                rows = value if isinstance(value, (list, tuple)) else [value] if value is not None else []
                values.update(FORMAT_ALIASES[str(row).strip().lower()] for row in rows
                              if str(row).strip().lower() in FORMAT_ALIASES)
            pending.extend(value for value in item.values() if isinstance(value, (dict, list, tuple)))
        return values

    candidate_sources = [candidate, raw, metrics]
    explicit = values_for(candidate_sources, ("binder_format", "binder_type"))
    if len(explicit) == 1:
        return next(iter(explicit))
    if explicit:
        return None
    # Per-candidate provenance is more specific than the campaign settings,
    # which may intentionally list several modalities for one campaign.
    configured = values_for(candidate_sources, ("modality", "modalities", "task", "design_task", "task_type"))
    if len(configured) == 1:
        return next(iter(configured))
    if configured:
        return None
    parameters = workflow_parameters or {}
    explicit = nested_values(parameters, ("binder_format", "binder_type"))
    if len(explicit) == 1:
        return next(iter(explicit))
    if explicit:
        return None
    configured = nested_values(parameters, ("modality", "modalities", "task", "design_task", "task_type"))
    if len(configured) == 1:
        return next(iter(configured))
    if configured:
        return None
    # Native tool identities supply a format only where that workflow's
    # documented scope is unambiguous.
    tool = str(candidate.get("source_tool") or "").lower()
    if tool in {"rfpeptides", "cyclic_peptide"}:
        return "cyclic_peptide"
    if tool in {"peptide_trap", "kejia_idr"}:
        return "peptide"
    return None


def infer_target_type(candidate: dict, workflow_parameters: dict | None = None) -> str | None:
    """Read the source campaign's target representation when recorded."""
    raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
    values = [candidate.get("target_type"), raw.get("target_type")]
    params = workflow_parameters or {}
    pending = [params]
    fasta_path_found = False
    while pending:
        item = pending.pop()
        if isinstance(item, dict):
            value = item.get("target_type")
            values.extend(value if isinstance(value, list) else [value] if value is not None else [])
            for key in ("target_path", "target_fasta", "target_fasta_path", "target_sequence_file"):
                path = item.get(key)
                if isinstance(path, str) and Path(path).suffix.lower() in {".fa", ".fasta", ".faa"}:
                    fasta_path_found = True
            pending.extend(value for value in item.values() if isinstance(value, (dict, list)))
        elif isinstance(item, list):
            pending.extend(item)
    normalized = {TARGET_ALIASES[str(value).strip().lower()] for value in values
                  if str(value).strip().lower() in TARGET_ALIASES}
    if len(normalized) > 1:
        return None
    if normalized:
        return next(iter(normalized))
    return "fasta_idr" if fasta_path_found else None

# These are selectable definitions, not default quality thresholds. BindCraft2's
# numeric defaults are checkpoint- and modality-specific and remain native
# there; users configure shared thresholds for the exact independent engine.
METRICS = {
    "ptm": {"label": "pTM · predicted complex confidence", "field": PREFIX + "confidence.validation.complex.ptm", "unit": "0–1", "direction": "higher", "modes": ("complex",),
            "engines": tuple(engine for engine, values in SCALAR_SCORE_ALIASES.items() if values.get("ptm"))},
    "iptm": {"label": "iPTM · interface confidence", "field": PREFIX + "confidence.validation.complex.iptm", "unit": "0–1", "direction": "higher", "modes": ("complex",),
              "engines": tuple(engine for engine, values in SCALAR_SCORE_ALIASES.items() if values.get("iptm"))},
    "predicted_pae_binder_to_target": {"label": "Predicted binder→target mean PAE", "field": PREFIX + "confidence.validation.interface.pae_binder_to_target", "unit": "Å", "direction": "lower", "modes": ("complex",),
              "engines": tuple(sorted(PAE_MATRIX_ENGINES))},
    "predicted_pae_bidirectional": {"label": "Predicted bidirectional mean PAE", "field": PREFIX + "confidence.validation.interface.pae_bidirectional", "unit": "Å", "direction": "lower", "modes": ("complex",),
              "engines": tuple(sorted(PAE_MATRIX_ENGINES))},
    "predicted_pae_binder_to_target_norm31": {"label": "Predicted binder→target mean PAE / 31 Å", "field": PREFIX + "confidence.validation.interface.pae_binder_to_target_norm31", "unit": "PAE / 31 Å", "direction": "lower", "modes": ("complex",),
              "engines": tuple(sorted(PAE_MATRIX_ENGINES))},
    "predicted_pae_bidirectional_norm31": {"label": "Predicted bidirectional mean PAE / 31 Å", "field": PREFIX + "confidence.validation.interface.pae_bidirectional_norm31", "unit": "PAE / 31 Å", "direction": "lower", "modes": ("complex",),
              "engines": tuple(sorted(PAE_MATRIX_ENGINES))},
    "rfd4_native_ptm": {"label": "RFD4 native oracle · complex pTM", "field": PREFIX + "confidence.generator.complex.ptm", "unit": "0–1", "direction": "higher", "modes": ("complex",), "source_tools": ("rfd4_proteina",)},
    "rfd4_native_iptm": {"label": "RFD4 native oracle · complex iPTM", "field": PREFIX + "confidence.generator.complex.iptm", "unit": "0–1", "direction": "higher", "modes": ("complex",), "source_tools": ("rfd4_proteina",)},
    "rfd4_native_plddt": {"label": "RFD4 native oracle · complex pLDDT", "field": PREFIX + "confidence.generator.complex.plddt", "unit": "0–100", "direction": "higher", "modes": ("complex",), "source_tools": ("rfd4_proteina",), "informational_only": True},
    "rfd4_native_interface_pae": {"label": "RFD4 native oracle · chain-of-interest iPAE", "field": PREFIX + "confidence.generator.interface.ipae_coi", "unit": "Å", "direction": "lower", "modes": ("complex",), "source_tools": ("rfd4_proteina",), "informational_only": True},
    "interface_pae_norm31": {"label": "Native interface PAE normalized by 31 Å", "field": PREFIX + "confidence.generator.interface.pae_norm31", "unit": "PAE / 31 Å", "direction": "lower", "modes": ("complex",), "source_tools": ("rfd4_proteina", "bindcraft2")},
    "rfd4_interface_residues_4a": {"label": "Saved RFD4 complex · binder residues within 4 Å", "field": PREFIX + "interface.native.rfd4_proteina.binder_residue_count_4a", "unit": "residues", "direction": "higher", "modes": ("complex",), "source_tools": ("rfd4_proteina",)},
    "rfd4_backbone_clashes_2_5a": {"label": "Saved RFD4 complex · interchain Cα clashes ≤2.5 Å", "field": PREFIX + "geometry.native.rfd4_proteina.backbone_ca_clashes_2_5a", "unit": "atom pairs", "direction": "lower", "modes": ("complex",), "source_tools": ("rfd4_proteina",)},
    "ipsae_min": {"label": "ipSAE · minimum direction", "field": PREFIX + "interface.ipsae.minimum", "unit": "0–1", "direction": "higher", "modes": ("complex",), "engines": ("ipsae", "boltz2", "proteina_complexa")},
    "ipsae_contact_pae_binder_to_target": {"label": "ipSAE analysis · contact PAE binder→target", "field": PREFIX + "confidence.analysis.interface.pae_binder_to_target", "unit": "Å", "direction": "lower", "modes": ("complex",), "engines": ("ipsae",)},
    "ipsae_contact_pae_target_to_binder": {"label": "ipSAE analysis · contact PAE target→binder", "field": PREFIX + "confidence.analysis.interface.pae_target_to_binder", "unit": "Å", "direction": "lower", "modes": ("complex",), "engines": ("ipsae",)},
    "ipsae_contact_pae_bidirectional": {"label": "ipSAE analysis · bidirectional contact PAE", "field": PREFIX + "confidence.analysis.interface.pae_bidirectional", "unit": "Å", "direction": "lower", "modes": ("complex",), "engines": ("ipsae",)},
    "ipsae_contact_pae_binder_to_target_norm31": {"label": "ipSAE analysis · contact PAE binder→target /31", "field": PREFIX + "confidence.analysis.interface.pae_binder_to_target_norm31", "unit": "PAE / 31 Å", "direction": "lower", "modes": ("complex",), "engines": ("ipsae",)},
    "ipsae_contact_pae_target_to_binder_norm31": {"label": "ipSAE analysis · contact PAE target→binder /31", "field": PREFIX + "confidence.analysis.interface.pae_target_to_binder_norm31", "unit": "PAE / 31 Å", "direction": "lower", "modes": ("complex",), "engines": ("ipsae",)},
    "ipsae_contact_pae_bidirectional_norm31": {"label": "ipSAE analysis · bidirectional contact PAE / 31 Å", "field": PREFIX + "confidence.analysis.interface.pae_bidirectional_norm31", "unit": "PAE / 31 Å", "direction": "lower", "modes": ("complex",), "engines": ("ipsae",)},
    "ipsae_bindcraft2_ipdae": {"label": "ipSAE analysis · BC2-definition i_pDAE (Studio recomputed)", "field": PREFIX + "interface.analysis.bindcraft2_definition.ipdae", "unit": "0–1", "direction": "higher", "modes": ("complex",), "engines": ("ipsae",)},
    "ipsae_structure_binder_interface_residues": {"label": "ipSAE & Structure Analysis · binder interface residues", "field": PREFIX + "interface.analysis.ipsae_structure.binder_interface_residues", "unit": "residues", "direction": "higher", "modes": ("complex",), "engines": ("ipsae",)},
    "ipsae_structure_residue_contact_pairs": {"label": "ipSAE & Structure Analysis · residue contact pairs", "field": PREFIX + "interface.analysis.ipsae_structure.residue_contact_pairs", "unit": "pairs", "direction": "higher", "modes": ("complex",), "engines": ("ipsae",)},
    "ipsae_design_contact_preservation_fraction": {"label": "Designed interface contact preservation fraction", "field": PREFIX + "interface.analysis.ipsae_structure.design_contact_preservation_fraction", "unit": "fraction", "direction": "higher", "modes": ("complex",), "engines": ("ipsae",)},
    "ipsae_design_cofold_contact_jaccard": {"label": "Designed/cofold interface contact Jaccard", "field": PREFIX + "interface.analysis.ipsae_structure.contact_jaccard", "unit": "fraction", "direction": "higher", "modes": ("complex",), "engines": ("ipsae",)},
    "ipsae_binder_interface_residue_preservation": {"label": "Binder interface residue preservation fraction", "field": PREFIX + "interface.analysis.ipsae_structure.binder_residue_preservation_fraction", "unit": "fraction", "direction": "higher", "modes": ("complex",), "engines": ("ipsae",)},
    "ipsae_target_interface_residue_preservation": {"label": "Target interface residue preservation fraction", "field": PREFIX + "interface.analysis.ipsae_structure.target_residue_preservation_fraction", "unit": "fraction", "direction": "higher", "modes": ("complex",), "engines": ("ipsae",)},
    "ipsae_structure_atom_contact_pairs": {"label": "ipSAE & Structure Analysis · heavy-atom contact pairs", "field": PREFIX + "interface.analysis.ipsae_structure.atom_contact_pairs", "unit": "pairs", "direction": "higher", "modes": ("complex",), "engines": ("ipsae",)},
    "ipsae_structure_ca_clashes_2_5a": {"label": "ipSAE & Structure Analysis · inter-chain Cα clashes", "field": PREFIX + "geometry.analysis.ipsae_structure.ca_clashes_2_5a", "unit": "pairs", "direction": "lower", "modes": ("complex",), "engines": ("ipsae",)},
    "ipsae_structure_binder_ca_bfactor_mean": {"label": "ipSAE & Structure Analysis · binder mean Cα B-factor (informational)", "field": PREFIX + "confidence.analysis.structure.binder_ca_bfactor_mean", "unit": "B-factor units", "direction": "", "modes": ("complex",), "engines": ("ipsae",), "informational_only": True},
    "ipsae_structure_target_ca_bfactor_mean": {"label": "ipSAE & Structure Analysis · target mean Cα B-factor (informational)", "field": PREFIX + "confidence.analysis.structure.target_ca_bfactor_mean", "unit": "B-factor units", "direction": "", "modes": ("complex",), "engines": ("ipsae",), "informational_only": True},
    "ipsae_binder_design_to_cofold_fitted_ca_rmsd": {"label": "Binder design-to-cofold fitted Cα RMSD", "field": PREFIX + "geometry.analysis.binder.design_to_cofold_fitted_ca_rmsd", "unit": "Å", "direction": "lower", "modes": ("complex",), "engines": ("ipsae",)},
    "ipsae_hotspot_contact_fraction_found": {"label": "Configured hotspot contact fraction (found residues)", "field": PREFIX + "interface.analysis.hotspots.contact_fraction_found", "unit": "fraction", "direction": "higher", "modes": ("complex",), "engines": ("ipsae",)},
    "ipsae_hotspot_contact_fraction_requested": {"label": "Configured hotspot contact fraction (requested residues)", "field": PREFIX + "interface.analysis.hotspots.contact_fraction_requested", "unit": "fraction", "direction": "higher", "modes": ("complex",), "engines": ("ipsae",)},
    "ipsae_hotspot_contacted_count": {"label": "Configured hotspot contacted count", "field": PREFIX + "interface.analysis.hotspots.contacted_count", "unit": "residues", "direction": "higher", "modes": ("complex",), "engines": ("ipsae",)},
    "ipsae_hotspot_minimum_distance": {"label": "Nearest configured hotspot distance", "field": PREFIX + "interface.analysis.hotspots.minimum_distance", "unit": "Å", "direction": "lower", "modes": ("complex",), "engines": ("ipsae",)},
    "bound_binder_plddt": {"label": "Bound binder structure mean Cα pLDDT", "field": PREFIX + "confidence.validation.structure.binder_ca_plddt", "unit": "0–100", "direction": "higher", "modes": ("complex",)},
    "target_plddt": {"label": "Target structure mean Cα pLDDT", "field": PREFIX + "confidence.validation.structure.target_ca_plddt", "unit": "0–100", "direction": "higher", "modes": ("complex",), "targets": ("fasta_idr",), "informational_only": True},
    "complex_structure_ca_plddt": {"label": "Structure mean Cα pLDDT", "field": PREFIX + "confidence.validation.structure.mean_ca_plddt", "unit": "0–100", "direction": "higher", "modes": ("complex", "binder")},
    "target_aligned_binder_rmsd": {"label": "Target-aligned binder pose RMSD", "field": PREFIX + "geometry.validation.binder.target_aligned_rmsd", "unit": "Å", "direction": "lower", "modes": ("complex",), "targets": ("structured",)},
    "binder_self_consistency_rmsd": {"label": "Binder self-consistency RMSD", "field": PREFIX + "geometry.validation.binder.self_consistency_rmsd", "unit": "Å", "direction": "lower", "modes": ("complex", "binder")},
    "interface_residues_4a": {"label": "Binder interface residues · any heavy atom ≤4 Å", "field": PREFIX + "interface.validation.binder_residue_count_4a", "unit": "residues", "direction": "higher", "modes": ("complex",)},
    "backbone_clashes_2_5a": {"label": "Interchain Cα clashes ≤2.5 Å", "field": PREFIX + "geometry.validation.interface.backbone_ca_clashes_2_5a", "unit": "atom pairs", "direction": "lower", "modes": ("complex",)},
    "bindcraft2_ipdae": {"label": "Source design · BindCraft2 i_pDAE", "field": PREFIX + "confidence.native.bindcraft2.interface.ipdae", "unit": "0–1", "direction": "higher", "modes": ("complex",), "source_tools": ("bindcraft2",), "informational_only": True},
    "bindcraft2_ipae": {"label": "Source design · BindCraft2 normalized i_pAE", "field": PREFIX + "confidence.native.bindcraft2.interface.ipae_normalized", "unit": "0–1", "direction": "lower", "modes": ("complex",), "source_tools": ("bindcraft2",), "informational_only": True},
    "bindcraft2_bound_unbound_rmsd": {"label": "Source design · BindCraft2 Binder_RMSD", "field": PREFIX + "geometry.native.bindcraft2.bound_unbound_rmsd", "unit": "Å", "direction": "lower", "modes": ("complex",), "source_tools": ("bindcraft2",), "informational_only": True},
    "bindcraft2_unbound_plddt": {"label": "Source design · BindCraft2 unbound binder pLDDT", "field": PREFIX + "confidence.validation.unbound_binder.plddt", "unit": "0–100", "direction": "higher", "modes": ("complex",), "source_tools": ("bindcraft2",), "informational_only": True},
    "bindcraft2_scaffold_sequence_retained": {"label": "Source design · scaffold sequence retained", "field": PREFIX + "sequence.native.bindcraft2.scaffold_sequence_retained_fraction", "unit": "fraction", "direction": "higher", "modes": ("complex",), "source_tools": ("bindcraft2",), "informational_only": True},
    "bindcraft2_scaffold_framework_rmsd": {"label": "Source design · scaffold framework RMSD", "field": PREFIX + "geometry.native.bindcraft2.scaffold_framework_rmsd", "unit": "Å", "direction": "lower", "modes": ("complex",), "source_tools": ("bindcraft2",), "informational_only": True},
    "bindcraft2_off_paratope_contacts": {"label": "Source design · off-paratope contact fraction", "field": PREFIX + "interface.native.bindcraft2.off_paratope_contact_fraction", "unit": "fraction", "direction": "lower", "modes": ("complex",), "source_tools": ("bindcraft2",), "informational_only": True},
    "bindcraft2_framework_packing": {"label": "Source design · framework packing fraction", "field": PREFIX + "structure.native.bindcraft2.framework_packing_fraction", "unit": "fraction", "direction": "higher", "modes": ("complex",), "source_tools": ("bindcraft2",), "informational_only": True},
    "bindcraft2_induced_fit_rmsd": {"label": "Source design · induced-fit interface RMSD", "field": PREFIX + "geometry.native.bindcraft2.induced_fit_interface_rmsd", "unit": "Å", "direction": "higher", "modes": ("complex",), "source_tools": ("bindcraft2",), "informational_only": True},
    "bindcraft2_induced_fit_tm": {"label": "Source design · induced-fit TM", "field": PREFIX + "confidence.native.bindcraft2.induced_fit_tm", "unit": "0–1", "direction": "lower", "modes": ("complex",), "source_tools": ("bindcraft2",), "informational_only": True},
    "bindcraft2_interdomain_contacts": {"label": "Source design · interdomain contact fraction", "field": PREFIX + "geometry.native.bindcraft2.interdomain_contact_fraction", "unit": "fraction", "direction": "lower", "modes": ("complex",), "source_tools": ("bindcraft2",), "informational_only": True},
    "bindcraft2_domain_separation": {"label": "Source design · domain separation ratio", "field": PREFIX + "geometry.native.bindcraft2.domain_separation_ratio", "unit": "ratio", "direction": "higher", "modes": ("complex",), "source_tools": ("bindcraft2",), "informational_only": True},
    "bindcraft2_interface_residues": {"label": "Source design · interface residues", "field": PREFIX + "interface.native.bindcraft2.binder_residue_count_4a", "unit": "residues", "direction": "higher", "modes": ("complex",), "source_tools": ("bindcraft2",), "informational_only": True},
    "bindcraft2_backbone_clashes": {"label": "Source design · Cα backbone clashes", "field": PREFIX + "geometry.native.bindcraft2.backbone_ca_clashes_2_5a", "unit": "atom pairs", "direction": "lower", "modes": ("complex",), "source_tools": ("bindcraft2",), "informational_only": True},
    "cyclic_closure_distance": {"label": "Cyclic peptide head-to-tail closure distance", "field": PREFIX + "geometry.validation.cyclic_peptide.closure_distance", "unit": "Å", "direction": "lower", "modes": ("complex", "binder"), "formats": ("cyclic_peptide",), "requires_cyclic_support": True},
}

# Defaults are copied from the installed BindCraft2 source settings/core/default.json,
# settings/modality/*.json, and bindcraft/settings.py. Only criteria with a direct
# measurement in the selected independent prediction are activated below. For
# example, a complex-only prediction cannot report Unbound_Binder_pLDDT, and the
# target-aligned RMSD here is not BindCraft2's bound-versus-unbound Binder_RMSD.
BINDCRAFT2_REFERENCE = {
    "version": "BindCraft2 source 4a56313",
    "source": "BindCraft2 settings/core/default.json + settings/modality/*.json + bindcraft/settings.py",
    "shared": {
        "ptm": (">=", 0.55), "iptm": (">=", 0.70),
        "interface_residues_4a": (">=", 7), "backbone_clashes_2_5a": ("<=", 0),
    },
    "formats": {
        "miniprotein": {"unbound_plddt": 0.70, "ptm": 0.55, "ipae": 0.35, "binder_rmsd": 3.5},
        "large_protein": {"unbound_plddt": 0.60, "ptm": 0.55, "ipae": 0.35, "binder_rmsd": 3.5},
        "peptide": {"unbound_plddt": None, "ptm": None, "ipae": 0.33},
        "cyclic_peptide": {"unbound_plddt": None, "ptm": None, "ipae": 0.33, "cyclic_closure": None},
        "vhh": {"unbound_plddt": 0.70, "ptm": 0.55, "ipae": 0.35,
                "scaffold_sequence_retained": 1.0, "off_paratope_contact_fraction": 0.20},
        "scfv": {"unbound_plddt": 0.70, "ptm": 0.55, "ipae": 0.35,
                 "scaffold_sequence_retained": 1.0, "scaffold_framework_rmsd": 3.0,
                 "off_paratope_contact_fraction": 0.20},
        "fab": {"unbound_plddt": 0.70, "ptm": 0.55, "ipae": 0.35,
                "scaffold_sequence_retained": 1.0, "scaffold_framework_rmsd": 3.0,
                "off_paratope_contact_fraction": 0.20},
        "arp": {"unbound_plddt": 0.70, "ptm": 0.55, "ipae": 0.35,
                "scaffold_sequence_retained": 1.0, "scaffold_framework_rmsd": 3.0,
                "off_paratope_contact_fraction": 0.20},
        "homo_oligomer": {"unbound_plddt": 0.70, "ptm": 0.55, "ipae": 0.35, "binder_rmsd": 3.5},
        "multidomain": {"unbound_plddt": 0.70, "ptm": 0.40, "ipae": 0.45,
                        "trajectory_plddt": 0.75, "interdomain_contact_fraction": 0.05,
                        "domain_separation_ratio": 1.0},
        "induced_fit": {"unbound_plddt": 0.70, "ptm": 0.55, "ipae": 0.35,
                        "induced_fit_interface_rmsd": 5.0},
        "fold_switch": {"unbound_plddt": 0.70, "ptm": 0.55, "ipae": 0.35,
                        "induced_fit_tm": 0.60},
    },
}


def bindcraft2_reference_policy(binder_format: str, target_type: str, mode: str,
                                available: dict[str, dict]) -> tuple[list[dict], list[str]]:
    """Map BindCraft2's format-aware defaults only to compatible measurements."""
    format_spec = BINDCRAFT2_REFERENCE["formats"].get(binder_format, {})
    shared = dict(BINDCRAFT2_REFERENCE["shared"])
    if "ptm" in format_spec:
        if format_spec["ptm"] is None:
            shared.pop("ptm", None)
        else:
            shared["ptm"] = (">=", format_spec["ptm"])
    criteria = dict(shared)
    if format_spec.get("ipae") is not None:
        criteria["ipae"] = ("<=", format_spec["ipae"])
    inactive = []
    if target_type == "fasta_idr":
        # The target is supplied as sequence because it is intrinsically
        # disordered. Its local structure confidence is not a binder acceptance
        # criterion, and target-pose RMSD has no fixed reference structure.
        inactive.append("Target pLDDT and target-aligned binder RMSD: not used for sequence-defined IDR targets.")
    metric_aliases = {
        "ptm": "rfd4_native_ptm",
        "iptm": "rfd4_native_iptm",
        "ipae": "interface_pae_norm31",
        "interface_residues_4a": "rfd4_interface_residues_4a",
        "backbone_clashes_2_5a": "rfd4_backbone_clashes_2_5a",
    }
    compatible = []
    for metric_id, (comparison, threshold) in criteria.items():
        selected_metric_id = metric_aliases.get(metric_id) if metric_aliases.get(metric_id) in available else metric_id
        spec = available.get(selected_metric_id)
        if spec is None or mode not in spec.get("modes", ()):
            inactive.append(f"{metric_id}: no compatible measurement is available for this source and mode")
            continue
        compatible.append({"field": spec["field"], "operator": comparison,
                           "value": threshold, "metric_id": selected_metric_id})
    if mode == "complex":
        unsupported = [f"Unbound_Binder_pLDDT ≥{format_spec.get('unbound_plddt', 0.8)}"]
        if format_spec.get("ipae") is not None and not any(row.get("metric_id") == "interface_pae_norm31" for row in compatible):
            unsupported.append(f"Interface PAE normalized by 31 Å ≤{format_spec['ipae']}")
        if "rfd4_native_interface_pae" in available and "interface_pae_norm31" not in available:
            inactive.append("RFD4 raw chain-of-interest iPAE is recorded as information only because no verified 31 Å normalization is available.")
        if format_spec.get("binder_rmsd") is not None:
            unsupported.append(f"Binder_RMSD ≤{format_spec['binder_rmsd']} Å")
        format_filter_names = {
            "scaffold_sequence_retained": "Scaffold_Sequence_Retained_Fraction",
            "off_paratope_contact_fraction": "Off_Paratope_Contact_Fraction",
            "scaffold_framework_rmsd": "Scaffold_Framework_RMSD",
            "interdomain_contact_fraction": "Interdomain_Contact_Fraction",
            "domain_separation_ratio": "Domain_Separation_Ratio",
            "induced_fit_interface_rmsd": "Induced_Fit_Interface_RMSD",
            "induced_fit_tm": "Induced_Fit_TM",
        }
        for key, filter_name in format_filter_names.items():
            if format_spec.get(key) is not None:
                value = format_spec[key]
                direction = "≥" if key in {"scaffold_sequence_retained", "domain_separation_ratio", "induced_fit_interface_rmsd"} else "≤"
                unsupported.append(f"{filter_name} {direction}{value}")
        inactive.extend(f"{name}: no compatible measurement in this run" for name in unsupported)
    return compatible, inactive


def normalize_policy(policy: dict | None) -> dict:
    policy = copy.deepcopy(policy or {})
    binder_format = str(policy.get("binder_format") or "miniprotein")
    target_type = str(policy.get("target_type") or "structured")
    if binder_format not in BINDER_FORMATS:
        raise ValueError(f"Unknown binder format {binder_format!r}.")
    if target_type not in TARGET_TYPES:
        raise ValueError(f"Unknown target type {target_type!r}.")
    filters = policy.get("mandatory_filters", policy.get("rules", []))
    scores = policy.get("informational_scores", [])
    if not isinstance(filters, list) or not isinstance(scores, list):
        raise ValueError("Mandatory filters and informational scores must be lists.")
    normalized = []
    seen = set()
    for rule in filters:
        if not isinstance(rule, dict) or not isinstance(rule.get("field"), str) or not rule["field"]:
            raise ValueError("Each mandatory filter needs an exact measurement field.")
        field = rule["field"].strip()
        op = rule.get("operator")
        value = finite_number(rule.get("value"))
        if op not in COMPARISONS or value is None:
            raise ValueError("Each mandatory filter needs a comparison operator and finite threshold.")
        if field in seen:
            raise ValueError(f"Duplicate mandatory filter for {field}.")
        seen.add(field)
        normalized.append({"field": field, "operator": op, "value": value,
                           "metric_id": str(rule.get("metric_id") or "")})
    for field in scores:
        if not isinstance(field, str) or not field.strip():
            raise ValueError("Informational score fields must be nonempty strings.")
    policy.update(schema="studio.validation-acceptance.v1", binder_format=binder_format,
                  target_type=target_type, mandatory_filters=normalized,
                  informational_scores=list(dict.fromkeys(field.strip() for field in scores)),
                  ranking_metric=str(policy.get("ranking_metric") or ""))
    if policy["ranking_metric"] and policy["ranking_metric"] not in policy["informational_scores"]:
        policy["informational_scores"].append(policy["ranking_metric"])
    return policy


def available_metrics(binder_format: str, target_type: str, mode: str, *, engine: str | None = None,
                      source_tools=(), informational: bool = False) -> list[tuple[str, dict]]:
    """Return metrics that make sense to offer for the selected context."""
    if binder_format not in BINDER_FORMATS or target_type not in TARGET_TYPES or mode not in {"complex", "binder"}:
        return []
    rows = []
    for metric_id, spec in METRICS.items():
        if mode not in spec["modes"]:
            continue
        if spec.get("targets") and target_type not in spec["targets"]:
            continue
        if spec.get("formats") and binder_format not in spec["formats"]:
            continue
        if spec.get("informational_only") and not informational:
            continue
        if spec.get("engines") and engine and engine not in spec["engines"]:
            continue
        if spec.get("requires_cyclic_support") and engine:
            from mn_molecular_design_studio.engine.workflow_engine.workflows.optimized_prediction import ENGINE_PROFILES
            if not ENGINE_PROFILES.get(engine, {}).get("supports_cyclic"):
                continue
        if spec.get("source_tools"):
            selected_tools = set(source_tools)
            if not selected_tools or not selected_tools.issubset(set(spec["source_tools"])):
                continue
        # A short free peptide need not fold alone. Keep those measures out of
        # mandatory-filter choices while retaining them for inspection.
        if not informational and binder_format in {"peptide", "cyclic_peptide"} and metric_id in {"binder_self_consistency_rmsd", "bindcraft2_bound_unbound_rmsd"}:
            continue
        rows.append((metric_id, spec))
    return rows


def _observe(candidate: dict, field: str) -> dict:
    value, measurement_candidate = _selected_measurement(candidate, field)
    if value is None and not str(field).startswith(("prediction/", "analysis/")):
        # Native source-design scores remain available for display/ranking, but
        # are explicitly attributed to the source candidate rather than being
        # presented as measurements from the independent prediction engine.
        source_candidate = (candidate.get("raw_metadata") or {}).get("source_candidate")
        if isinstance(source_candidate, dict):
            source_value = metric_value(source_candidate, field)
            if source_value is not None:
                value, measurement_candidate = source_value, source_candidate
    measurement_candidate = measurement_candidate or candidate
    if value is None:
        return {"field": field, "value": None, "status": "missing_evidence", "error": "Measurement is unavailable."}
    normalized = candidate_output(measurement_candidate)
    if field.startswith(PREFIX):
        key = field[len(PREFIX):]
        definition = DEFINITIONS.get(key)
        low, high = (definition.minimum, definition.maximum) if definition else (None, None)
        if any(issue.get("field") == key for issue in normalized.get("issues", [])):
            return {"field": field, "value": None, "status": "invalid_metadata", "error": "Measurement violates its output contract."}
        if any(any(alias.get("source_field") == key for alias in aliases) for aliases in normalized.get("conflicts", {}).values()):
            return {"field": field, "value": None, "status": "invalid_metadata", "error": "Measurement has conflicting source aliases."}
    else:
        low, high = RANGES.get(field.split("/")[-1], (None, None))
    number = finite_number(value)
    if number is None or low is not None and number < low or high is not None and number > high:
        return {"field": field, "value": None, "status": "invalid_metadata", "raw_value": value,
                "error": "Measurement is nonfinite or outside its declared value range."}
    provenance = normalized.get("provenance") or {}
    metric_provenance = provenance.get(field[len(PREFIX):], {}) if field.startswith(PREFIX) else {}
    return {"field": field, "value": number, "status": "measured", "unit": (DEFINITIONS.get(field[len(PREFIX):]).unit if field.startswith(PREFIX) and DEFINITIONS.get(field[len(PREFIX):]) else None),
            "source_tool": measurement_candidate.get("source_tool"), "stage": measurement_candidate.get("stage"),
            "measurement_candidate_id": measurement_candidate.get("candidate_id"),
            "method": metric_provenance.get("method"), "output_provenance": provenance}


def _selected_measurement(candidate: dict, field: str):
    """Read a candidate metric or one exact completed graph measurement.

    Graph measurements keep their branch and prediction identity. An ambiguous
    branch-only selector is missing evidence, never an implicit best score.
    """
    if not str(field).startswith(("prediction/", "analysis/")):
        return metric_value(candidate, field), candidate
    parts = str(field).split("/")
    if len(parts) not in {3, 4} or not parts[-1]:
        raise ValueError("Use prediction/<branch>/<metric> or analysis/<branch>/<metric>, optionally naming a result ID.")
    family, branch, *rest = parts
    info = (candidate.get("raw_metadata") or {}).get("candidate_set") or {}
    key, id_key = (("predictions", "prediction_id") if family == "prediction"
                   else ("analysis_results", "analysis_result_id"))
    records = [row for row in info.get(key, []) if row.get("branch") == branch
               and row.get("status") == "completed"
               and (len(rest) == 1 or row.get(id_key) == rest[0])]
    if len(records) != 1:
        return None, None
    record = records[0]
    metric = rest[-1]
    evidence = record.get("candidate") or {}
    value = metric_value(evidence, metric)
    if value is not None:
        return value, evidence
    observation = (record.get("measurements") or {}).get(metric)
    return (observation.get("value"), evidence) if isinstance(observation, dict) else (None, evidence)


def evaluate_candidate(candidate: dict, policy: dict) -> dict:
    policy = normalize_policy(policy)
    measurements = {}
    reasons = []
    rejection_reasons = []
    unresolved_reasons = []
    for rule in policy["mandatory_filters"]:
        observation = _observe(candidate, rule["field"])
        passed = None
        if observation["status"] == "measured":
            passed = bool(COMPARISONS[rule["operator"]](observation["value"], rule["value"]))
            if not passed:
                rejection_reasons.append(f"{rule['field']} = {observation['value']:g} fails {rule['operator']} {rule['value']:g}.")
        else:
            unresolved_reasons.append(f"{rule['field']}: {observation['error']}")
        measurements[rule["field"]] = {**observation, "mandatory": True, "threshold": rule["value"], "operator": rule["operator"], "passed": passed}
    for field in policy["informational_scores"]:
        measurements.setdefault(field, {**_observe(candidate, field), "mandatory": False, "passed": None})
    if rejection_reasons:
        status = "rejected"
        reasons.extend(rejection_reasons)
        reasons.extend(unresolved_reasons)
    elif unresolved_reasons:
        status = "unresolved"
        reasons.extend(unresolved_reasons)
    elif policy["mandatory_filters"]:
        status = "accepted"
    else:
        status = "not_evaluated"
    rank = measurements.get(policy["ranking_metric"]) if policy["ranking_metric"] else None
    native_metrics = candidate.get("metrics") if isinstance(candidate.get("metrics"), dict) else {}
    native_status = native_metrics.get("native_pass_filters", native_metrics.get("pass_filters"))
    source_candidate = (candidate.get("raw_metadata") or {}).get("source_candidate")
    source_metrics = native_metrics
    if isinstance(source_candidate, dict):
        source_metrics = source_candidate.get("metrics") if isinstance(source_candidate.get("metrics"), dict) else source_metrics
    if native_status is None:
        native_status = source_metrics.get("native_pass_filters", source_metrics.get("pass_filters"))
        if native_status is None and isinstance(source_candidate, dict) and (source_candidate.get("raw_metadata") or {}).get("bindcraft_output_kind"):
            native_status = (source_candidate["raw_metadata"]["bindcraft_output_kind"] == "Accepted")
        if native_status is None and "ppi.pass" in source_metrics:
            all_metrics_present = source_metrics.get("ppi.all_metrics_present")
            if all_metrics_present not in (1, 1.0, True):
                native_status = "unavailable · native PPI metrics incomplete"
            else:
                native_status = "passed" if source_metrics["ppi.pass"] in (1, 1.0, True) else "failed"
    return {"status": status, "reasons": reasons, "measurements": measurements,
            "binder_format": policy["binder_format"], "target_type": policy["target_type"],
            "ranking_metric": policy["ranking_metric"] or None,
            "ranking_value": rank.get("value") if rank and rank.get("status") == "measured" else None,
            "policy": policy,
            "native_filter_status": native_status}


def campaign_policy_from_node(params: dict) -> dict:
    """Use exactly the same policy as the graph Validation & Acceptance node."""
    policy = normalize_policy({
        "binder_format": params.get("binder_format", "miniprotein"),
        "target_type": params.get("target_type", "structured"),
        "mandatory_filters": params.get("rules", []),
        "informational_scores": params.get("informational_scores", []),
        "ranking_metric": params.get("ranking_metric", ""),
    })
    if not policy["mandatory_filters"]:
        raise ValueError("A passing-design campaign requires mandatory filters on its Validation & Acceptance node.")
    return policy


def campaign_stage(policy: dict, *, engine: str = "boltz2", parameters: dict | None = None) -> dict:
    """Bind the independent prediction workflow to the shared acceptance backend."""
    policy = normalize_policy(policy)
    if not policy["mandatory_filters"]:
        raise ValueError("A passing-design objective requires at least one mandatory acceptance filter.")
    return {
        "workflow": "refolding.optimized_prediction",
        "parameters": {"engine": engine, "samples": 1, "optimization_mode": "off",
                       "use_target_template": False, **(parameters or {})},
        "acceptance_policy": policy,
        "required_metrics": [], "rules": [],
    }


def apply_acceptance(candidates: list[dict], policy: dict) -> tuple[list[dict], dict]:
    policy = normalize_policy(policy)
    rows = []
    counts = {key: 0 for key in ("accepted", "rejected", "unresolved", "not_evaluated")}
    for candidate in candidates:
        row = copy.deepcopy(candidate)
        decision = evaluate_candidate(row, policy)
        row.setdefault("raw_metadata", {})["validation_acceptance"] = decision
        counts[decision["status"]] += 1
        rows.append(row)
    if policy["ranking_metric"]:
        rows.sort(key=lambda row: (
            row["raw_metadata"]["validation_acceptance"]["status"] != "accepted",
            row["raw_metadata"]["validation_acceptance"]["ranking_value"] is None,
            -row["raw_metadata"]["validation_acceptance"]["ranking_value"]
            if (row["raw_metadata"]["validation_acceptance"]["ranking_value"] is not None
                and metric_direction(policy["ranking_metric"]) == "higher") else
            row["raw_metadata"]["validation_acceptance"]["ranking_value"] or 0,
            str(row.get("candidate_id") or ""),
        ))
    summary = {"schema": policy["schema"], "candidate_count": len(rows), **counts,
               "mandatory_filter_count": len(policy["mandatory_filters"]),
               "informational_score_fields": policy["informational_scores"],
               "ranking_metric": policy["ranking_metric"] or None,
               "policy": policy}
    return rows, summary


def metric_direction(field: str) -> str:
    for spec in METRICS.values():
        if spec["field"] == field:
            return spec["direction"]
    definition = DEFINITIONS.get(field.removeprefix(PREFIX)) if field.startswith(PREFIX) else None
    return definition.direction if definition and definition.direction in {"higher", "lower"} else "higher"


def write_report(run_dir: Path, candidates: list[dict], summary: dict) -> Path:
    path = Path(run_dir) / "artifacts/validation/acceptance.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["candidate_id", "source_candidate_ids", "native_filter_status", "status", "binder_format", "target_type", "ranking_metric", "ranking_value", "measurements_json", "reasons"])
        writer.writeheader()
        for candidate in candidates:
            decision = (candidate.get("raw_metadata") or {}).get("validation_acceptance") or {}
            provenance = (candidate.get("raw_metadata") or {}).get("validation_acceptance_provenance") or {}
            source_ids = candidate.get("parents") or ([provenance["source_candidate_id"]]
                                                        if provenance.get("source_candidate_id") else [])
            writer.writerow({"candidate_id": candidate.get("candidate_id"),
                             "source_candidate_ids": ";".join(map(str, source_ids)),
                             "native_filter_status": decision.get("native_filter_status"),
                             **{key: decision.get(key) for key in ("status", "binder_format", "target_type", "ranking_metric", "ranking_value", "reasons")},
                             "reasons": "; ".join(map(str, decision.get("reasons") or [])),
                             "measurements_json": json.dumps(decision.get("measurements", {}), sort_keys=True)})
    (path.parent / "acceptance_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True))
    return path
