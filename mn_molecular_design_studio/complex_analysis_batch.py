"""Small helpers for explicit post-design assessments of saved Database models."""
from __future__ import annotations

from typing import Any


ASSESSMENT_METHODS = {
    "analysis.sequence_qc": {
        "label": "Sequence Composition", "evaluators": {"sequence_qc"},
        "scope": "one explicitly selected protein chain", "kind": "sequence",
    },
    "analysis.proteinsol": {
        "label": "Protein-Sol", "evaluators": {"proteinsol"},
        "scope": "one protein sequence of at least 21 residues", "kind": "sequence",
    },
    "analysis.biotite_interface": {
        "label": "Biotite Interface", "evaluators": {"biotite_interface"},
        "scope": "protein–protein complex with explicit binder and target chains", "kind": "ppi",
    },
    "analysis.binder_interface_chemistry": {
        "label": "BinderScope-style Interface Surface & Chemistry",
        "evaluators": {"binder_interface_chemistry_v1"},
        "scope": "single-model protein–protein complex with explicit binder and target chains",
        "kind": "ppi",
    },
    "analysis.rosetta": {
        "label": "Rosetta Interface Analysis", "evaluators": {"rosetta_metrics", "rosetta"},
        "scope": "one explicit binder chain and one target chain", "kind": "single_ppi",
    },
    "analysis.ipsae": {
        "label": "ipSAE & Structure Analysis", "evaluators": {"ipsae"},
        "scope": "one explicit binder chain and one target chain", "kind": "single_ppi",
    },
    "analysis.pep_patch": {
        "label": "PEP-Patch", "evaluators": {"pep_patch"},
        "scope": "one or more explicitly selected protein surface chains", "kind": "protein_surface",
    },
    "analysis.sequence_liability": {
        "label": "Sequence and Developability Liabilities", "evaluators": {"sequence_liability"},
        "scope": "explicit protein chain(s), sequence and selected model", "kind": "protein_structure",
    },
    "analysis.tag_fusion_suitability": {
        "label": "Tag and Fusion Suitability", "evaluators": {"tag_fusion_suitability"},
        "scope": "one component, partner chain(s), and declared construct termini", "kind": "single_ppi",
    },
}

CORE_PROFILE = (
    "analysis.sequence_qc", "analysis.proteinsol", "analysis.biotite_interface",
    "analysis.pep_patch", "analysis.sequence_liability",
)
COMPREHENSIVE_PROFILE = CORE_PROFILE + (
    "analysis.binder_interface_chemistry", "analysis.rosetta", "analysis.ipsae", "analysis.tag_fusion_suitability",
)


def recorded_component_format(sample: dict[str, Any], binder_chains: list[str]) -> tuple[str, str]:
    """Suggest a component format from saved evidence, never from chain order.

    A target chain does not make the assessed binder a multichain construct.
    The design task establishes topology, not complete termini or B-factor
    semantics; those remain separately declared evaluator parameters.
    """
    if not binder_chains:
        return "", "Select the component to assess."
    formats = {"linear_protein": "single_chain_linear_protein",
               "single_chain_linear_protein": "single_chain_linear_protein",
               "cyclic_construct": "cyclic_peptide", "cyclic_peptide": "cyclic_peptide",
               "multichain_construct": "multichain_construct"}
    context = sample.get("context") or {}
    metadata = (sample.get("record") or {}).get("metadata") or {}
    recorded = {formats[value] for item in (context, metadata)
                for value in [item.get("component_format") or item.get("construct_format")]
                if isinstance(value, str) and value in formats}
    if not recorded:
        for observation, descriptor in sample.get("associated") or []:
            if descriptor.get("status", "completed") not in {"completed", "succeeded"}:
                continue
            if descriptor.get("binder_chains") != binder_chains:
                continue
            if descriptor.get("evaluator") not in {"tag_fusion_suitability", "sequence_liability"}:
                continue
            settings = (observation.get("context") or {}).get("settings") or {}
            value = settings.get("component_format")
            if isinstance(value, str) and value in formats:
                recorded.add(formats[value])
    if len(recorded) > 1:
        return "", "Saved component formats conflict; select the format explicitly."
    if recorded:
        value = next(iter(recorded))
        if value == "single_chain_linear_protein" and len(binder_chains) != 1:
            return "", "The saved single-chain format does not apply to this chain selection."
        return value, "Recorded component format / saved assessment settings."
    if (len(binder_chains) == 1 and metadata.get("task") == "protein_binder"
            and not metadata.get("unsupported_chemistry")):
        return "single_chain_linear_protein", "Recorded protein-binder design task; one assessed binder chain."
    return "", "No applicable component topology is recorded; select the format explicitly."


def analysis_chain_parameters(workflow_id: str, binder_chains: list[str], target_chains: list[str]) -> dict:
    """Keep surface-only selection separate from interface workflow arguments.

    PEP-Patch accepts the selected surface chains as binder_chains; it has no
    target_chains argument. Database source provenance retains both roles.
    """
    parameters = {"binder_chains": list(binder_chains)}
    if workflow_id != "analysis.pep_patch":
        parameters["target_chains"] = list(target_chains)
    return parameters


def offered_assessment_methods(samples: list[dict]) -> list[str]:
    """Offer existing protein evaluators by molecular input, independent of engine.

    Detailed chain, sequence and settings checks still happen at submission.
    Protein–ligand methods use the existing Database ligand-analysis adapter.
    """
    kinds = {sample.get("kind") for sample in samples}
    protein = bool(kinds & {"protein", "protein-protein"})
    return [method for method, profile in ASSESSMENT_METHODS.items()
            if protein and (profile["kind"] not in {"ppi", "single_ppi"} or "protein-protein" in kinds)]


def recorded_chain_roles(sample: dict[str, Any]) -> tuple[list[str], list[str]]:
    """Return only roles explicitly stored with this representation or its components."""
    candidate = sample.get("candidate") if isinstance(sample.get("candidate"), dict) else {}
    context = sample.get("context") if isinstance(sample.get("context"), dict) else {}
    components = sample.get("components") if isinstance(sample.get("components"), list) else []
    protein_chains = {str(item.get("source_chain") or item.get("chain") or "") for item in components
                      if isinstance(item, dict) and str(item.get("type") or "").casefold() == "protein"}
    binder = candidate.get("binder_chains") or context.get("binder_chains") or []
    target = candidate.get("target_chains") or context.get("target_chains") or []
    if not binder:
        binder = [str(item.get("source_chain") or item.get("chain") or "") for item in components
                  if isinstance(item, dict) and str(item.get("type") or "").casefold() == "protein"
                  and str(item.get("role") or "").casefold() == "binder"]
    if not target:
        target = [str(item.get("source_chain") or item.get("chain") or "") for item in components
                  if isinstance(item, dict) and str(item.get("type") or "").casefold() == "protein"
                  and str(item.get("role") or "").casefold() == "target"]
    binder = sorted({str(value) for value in binder if str(value) in protein_chains})
    target = sorted({str(value) for value in target if str(value) in protein_chains and str(value) not in binder})
    return binder, target


def assessment_applicability(sample: dict[str, Any], method: str,
                             binder_chains: list[str], target_chains: list[str]) -> tuple[bool, str]:
    """Check data prerequisites without guessing chain roles or missing sequence."""
    profile = ASSESSMENT_METHODS[method]
    candidate = sample.get("candidate") if isinstance(sample.get("candidate"), dict) else {}
    components = sample.get("components") if isinstance(sample.get("components"), list) else []
    chain_sequences = {str(item.get("source_chain") or item.get("chain") or ""): str(item.get("sequence") or "")
                       for item in components if isinstance(item, dict)
                       and str(item.get("type") or "").casefold() == "protein"}
    is_ppi = sample.get("kind") == "protein-protein"
    if profile["kind"] == "sequence":
        if len(binder_chains) != 1 or not chain_sequences.get(binder_chains[0]):
            return False, "Select one protein chain with a parsed sequence."
        sequence = chain_sequences[binder_chains[0]]
        if method == "analysis.proteinsol" and len(sequence) < 21:
            return False, "Protein-Sol requires a sequence of at least 21 residues."
        return True, ""
    if profile["kind"] in {"ppi", "single_ppi"} and not is_ppi:
        return False, "This analysis applies to protein–protein structures."
    if profile["kind"] == "protein_structure" and sample.get("kind") not in {"protein", "protein-protein"}:
        return False, "This analysis requires a protein structure."
    if profile["kind"] == "protein_surface" and not binder_chains:
        return False, "Select at least one protein chain for surface analysis."
    if profile["kind"] in {"ppi", "single_ppi"} and (
        not binder_chains or not target_chains or set(binder_chains) & set(target_chains)
    ):
        return False, "Explicit, non-overlapping binder and target chain groups are required."
    if profile["kind"] == "single_ppi" and (len(binder_chains) != 1 or len(target_chains) != 1):
        return False, "The registered evaluator accepts exactly one chain in each group."
    if profile["kind"] == "protein_structure" and not binder_chains:
        return False, "Select at least one protein chain to assess."
    if profile["kind"] == "single_ppi" and method == "analysis.tag_fusion_suitability" and not target_chains:
        return False, "Choose at least one interface partner chain."
    return True, ""


def evaluation_coverage(sample: dict[str, Any], method: str, *, parameters: dict[str, Any],
                        binder_chains: list[str], target_chains: list[str]) -> tuple[str, str]:
    """Classify saved Evaluation evidence for exact source model and requested settings."""
    profile = ASSESSMENT_METHODS[method]
    matching_method = []
    exact = []
    for observation, descriptor in sample.get("associated") or []:
        evaluator = str((descriptor or {}).get("evaluator") or (observation or {}).get("namespace") or "").casefold()
        if evaluator not in profile["evaluators"]:
            continue
        context = (observation or {}).get("context") if isinstance((observation or {}).get("context"), dict) else {}
        candidate = context.get("candidate") if isinstance(context.get("candidate"), dict) else {}
        raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
        selection = raw.get("complex_analysis_selection") if isinstance(raw.get("complex_analysis_selection"), dict) else {}
        same_structure = (
            str(context.get("source_structure_id") or context.get("representation_id") or descriptor.get("representation_id") or "")
            == str(sample.get("structure_id") or "")
            and str(context.get("source_structure_digest") or context.get("structure_digest") or descriptor.get("source_structure_digest") or "")
            == str(sample.get("digest") or "")
        )
        if not same_structure:
            continue
        matching_method.append((observation, descriptor, context, candidate, selection))
        saved_binder = candidate.get("binder_chains") or selection.get("binder_chains") or context.get("binder_chains") or []
        saved_target = candidate.get("target_chains") or selection.get("target_chains") or context.get("target_chains") or []
        if sorted(map(str, saved_binder)) != sorted(map(str, binder_chains)):
            continue
        if sorted(map(str, saved_target)) != sorted(map(str, target_chains)):
            continue
        if str(selection.get("model_id") or candidate.get("model_id") or candidate.get("model_index") or context.get("model_id") or "") != str(sample.get("model_id") or ""):
            continue
        settings = context.get("settings") if isinstance(context.get("settings"), dict) else {}
        settings_match = all(settings.get(key) == value for key, value in parameters.items())
        if settings_match:
            exact.append((observation, descriptor, context))
    for _observation, descriptor, context in exact:
        status = str(context.get("status") or descriptor.get("status") or "completed").casefold()
        if status in {"completed", "succeeded"} and not context.get("error"):
            return "Already evaluated · compatible settings", str(descriptor.get("id") or "Compatible Evaluation is saved.")
        if status in {"failed", "error", "cancelled", "canceled"}:
            return "Previously failed", str(context.get("error") or "A matching Evaluation is recorded as failed.")
        return "In progress", str(context.get("error") or status.title())
    if matching_method:
        statuses = {str(context.get("status") or descriptor.get("status") or "completed").casefold()
                    for _obs, descriptor, context, _candidate, _selection in matching_method}
        if statuses & {"failed", "error", "cancelled", "canceled"}:
            return "Previously failed", "A prior Evaluation exists for this exact structure with other settings or chain selection."
        return "Existing Evaluation · different settings", "Same evaluator/source model exists, but requested parameters or chain selection differ."
    return "Missing evaluation", "No compatible Evaluation is linked to this exact Database representation."
