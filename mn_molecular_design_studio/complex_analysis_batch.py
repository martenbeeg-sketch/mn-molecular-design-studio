"""Small helpers for explicit post-design assessments of saved Database models."""
from __future__ import annotations

from typing import Any


_PROTOCOL_FIELDS = {"method_version", "engine_version", "protocol_revision", "image", "image_digest", "checkpoint_sha256"}
_PLACEMENT_FIELDS = {"source_run_dir", "candidates_jsonl", "candidate_ids", "selected_structure_path", "selected_model_id",
                     "binder_chains", "target_chains", "source_job", "job_name", "campaign_name", "backend",
                     "cpu_cores", "cpu_threads", "nprocs", "workers", "gpu_device", "output_dir", "method", "assessment_protocol"}


def scientific_parameters(parameters):
    """Keep scientific settings separate from protocol and source placement."""
    return {key: value for key, value in parameters.items()
            if key not in _PLACEMENT_FIELDS | _PROTOCOL_FIELDS and value is not None}


def assessment_protocol(method):
    """Expected existing method contract; no container/network probing on render.

    Mutable :latest engines cannot establish runtime compatibility from a tag.
    Callers with a verified runtime digest may supply it explicitly; otherwise
    reuse remains visibly unverified for those engines.
    """
    from mn_molecular_design_studio.quality_tools import WORKFLOWS
    if method in WORKFLOWS:
        from mn_molecular_design_studio.engine.workflow_engine.workflows.quality import protocol_revision
        return {"protocol_revision": protocol_revision(), "image": WORKFLOWS[method]["image"]}
    if method == "analysis.sequence_liability":
        from mn_molecular_design_studio.complex_analysis_liabilities import METHOD_VERSION
        return {"method_version": METHOD_VERSION}
    if method == "analysis.binder_interface_chemistry":
        from mn_molecular_design_studio.complex_analysis_interface_chemistry import METHOD_VERSION
        return {"method_version": METHOD_VERSION}
    if method == "analysis.tag_fusion_suitability":
        from mn_molecular_design_studio.complex_analysis_tag_fusion import METHOD_VERSION
        return {"method_version": METHOD_VERSION}
    return {"image_digest": ""}  # Unknown is not an implicit version match.


def protocol_matches(saved, expected):
    return bool(expected) and all(value not in (None, "") and saved.get(key) == value for key, value in expected.items())


def analysis_job_coverage(jobs):
    """Treat all nonterminal scheduler states conservatively, including unknown."""
    from mn_molecular_design_studio.cluster.store import TERMINAL
    successes = sorted((j for j in jobs if str(j.get("status") or "").casefold() in {"completed", "succeeded"}),
                       key=lambda j: str(j.get("run_id") or ""))
    active = sorted((j for j in jobs if str(j.get("status") or "unknown").casefold() not in TERMINAL | {"canceled", "error"}),
                    key=lambda j: str(j.get("run_id") or ""))
    if successes:
        return "Completed job · Evaluation ingestion pending", str(successes[0].get("run_id") or "Completed job")
    if active:
        job = active[0]
        return "Queued / running", str(job.get("run_id") or "Matching job") + " · " + str(job.get("status") or "unknown")
    if jobs:
        return "Previously failed", str(jobs[0].get("error") or jobs[0].get("run_id") or "Previous job failed")
    return None


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
        "scope": "one component, partner chain(s), and declared construct termini", "kind": "component_ppi",
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
            if protein and (profile["kind"] not in {"ppi", "single_ppi", "component_ppi"} or "protein-protein" in kinds)]


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
    if profile["kind"] in {"ppi", "single_ppi", "component_ppi"} and not is_ppi:
        return False, "This analysis applies to protein–protein structures."
    if profile["kind"] == "protein_structure" and sample.get("kind") not in {"protein", "protein-protein"}:
        return False, "This analysis requires a protein structure."
    if profile["kind"] == "protein_surface" and not binder_chains:
        return False, "Select at least one protein chain for surface analysis."
    if profile["kind"] in {"ppi", "single_ppi", "component_ppi"} and (
        not binder_chains or not target_chains or set(binder_chains) & set(target_chains)
    ):
        return False, "Explicit, non-overlapping binder and target chain groups are required."
    if profile["kind"] == "single_ppi" and (len(binder_chains) != 1 or len(target_chains) != 1):
        return False, "The registered evaluator accepts exactly one chain in each group."
    if profile["kind"] == "protein_structure" and not binder_chains:
        return False, "Select at least one protein chain to assess."
    if profile["kind"] == "component_ppi" and len(binder_chains) != 1:
        return False, "Select exactly one assessed component and one or more interface partner chains."
    return True, ""


def evaluation_coverage(sample: dict[str, Any], method: str, *, parameters: dict[str, Any],
                        binder_chains: list[str], target_chains: list[str], expected_protocol=None) -> tuple[str, str]:
    """Classify saved Evaluation evidence for exact source model and requested settings."""
    profile = ASSESSMENT_METHODS[method]
    matching_method = []
    exact = []
    expected_protocol = assessment_protocol(method) if expected_protocol is None else expected_protocol
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
        saved_model = next((v for v in (selection.get("model_id"), candidate.get("model_id"), candidate.get("model_index"), context.get("model_id"))
                            if v is not None and str(v) != ""), "")
        if str(saved_model) != str(sample.get("model_id") if sample.get("model_id") is not None else ""):
            continue
        settings = context.get("settings") if isinstance(context.get("settings"), dict) else {}
        saved_protocol = {**{key: descriptor[key] for key in _PROTOCOL_FIELDS if descriptor.get(key)},
                          **{key: raw[key] for key in _PROTOCOL_FIELDS if raw.get(key)},
                          **{key: settings[key] for key in _PROTOCOL_FIELDS if settings.get(key)},
                          **(settings.get("assessment_protocol") or {})}
        settings_match = (scientific_parameters(settings) == scientific_parameters(parameters)
                          and protocol_matches(saved_protocol, expected_protocol))
        if settings_match:
            exact.append((observation, descriptor, context))
    for _observation, descriptor, context in sorted(exact, key=lambda item: (
        0 if str(item[2].get("status") or item[1].get("status") or "completed").casefold() in {"completed", "succeeded"} and not item[2].get("error") else
        2 if str(item[2].get("status") or item[1].get("status") or "completed").casefold() in {"failed", "error", "cancelled", "canceled"} or item[2].get("error") else 1,
        str(item[1].get("id") or ""))):
        status = str(context.get("status") or descriptor.get("status") or "completed").casefold()
        if status in {"completed", "succeeded"} and not context.get("error"):
            return "Already evaluated · compatible settings", str(descriptor.get("id") or "Compatible Evaluation is saved.")
        if status in {"failed", "error", "cancelled", "canceled"} or context.get("error"):
            return "Previously failed", str(context.get("error") or "A matching Evaluation is recorded as failed.")
        return "In progress", str(context.get("error") or status.title())
    if matching_method:
        return "Existing Evaluation · different settings", "Saved evidence differs in settings, chain/model selection or method version; unknown protocol compatibility is not assumed."
    return "Missing evaluation", "No compatible Evaluation is linked to this exact Database representation."
