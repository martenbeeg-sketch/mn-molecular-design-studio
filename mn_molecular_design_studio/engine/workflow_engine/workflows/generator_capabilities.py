"""Explicit attempt, output, and acceptance contracts for design generators.

Campaign adapters are opt-in. An engine is not eligible merely because it can
write a structure: its attempt unit, seed behavior, normalized output, and
acceptance evidence must all be known.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Literal


@dataclass(frozen=True)
class GeneratorCapabilities:
    engine_id: str
    label: str
    campaign_adapter: bool
    independent_attempts: bool
    attempt_unit: str
    attempt_budget_semantics: str
    native_acceptance_available: bool | None
    native_acceptance_semantics: str
    candidate_output: str
    output_normalization: str
    seed_assignment: str
    progress_reporting: str
    cancellation: str
    sequence_structure_mode: Literal["co_design", "backbone_only", "sequence_only", "mixed", "unverified"]
    downstream_validation: tuple[str, ...]
    gpu_required: bool
    minimum_vram_gb: float | None
    supported_stopping_policies: tuple[str, ...]
    deferred_reason: str = ""

    def as_dict(self) -> dict[str, object]:
        value = asdict(self)
        value["downstream_validation"] = list(self.downstream_validation)
        value["supported_stopping_policies"] = list(self.supported_stopping_policies)
        return value


_CAPABILITIES = {
    "bindcraft2": GeneratorCapabilities(
        engine_id="bindcraft2", label="BindCraft2", campaign_adapter=True,
        independent_attempts=True, attempt_unit="one native trajectory",
        attempt_budget_semantics="one child requests one trajectory and one final design",
        native_acceptance_available=True,
        native_acceptance_semantics="authoritative native campaign-state and normalized-row reconciliation",
        candidate_output="normalized candidate rows with native filter decisions",
        output_normalization="BindCraft2 native trajectory and candidate normalizer",
        seed_assignment="one persisted campaign seed per attempt",
        progress_reporting="authoritative .campaign_state.json plus scheduler job status",
        cancellation="stop the scheduler-owned attempt child",
        sequence_structure_mode="mixed",
        downstream_validation=("refolding.optimized_prediction", "analysis.ipsae"),
        gpu_required=True, minimum_vram_gb=None,
        supported_stopping_policies=("native", "validated", "diversity_qualified"),
    ),
    "rfd4_proteina": GeneratorCapabilities(
        engine_id="rfd4_proteina", label="RFD4-Proteina", campaign_adapter=True,
        independent_attempts=True, attempt_unit="one native design request",
        attempt_budget_semantics=(
            "one child requests one design; the child may return zero or one normalized candidate"
        ),
        native_acceptance_available=False,
        native_acceptance_semantics=(
            "RFD4's task-specific condition/evaluation outputs are not an authoritative campaign acceptance flag"
        ),
        candidate_output="normalized co-designed sequence/structure candidate rows",
        output_normalization="RFD4-Proteina result-index normalizer with chain-role and sequence verification",
        seed_assignment="one persisted campaign seed mapped to the native RFD4 request",
        progress_reporting="scheduler job status and RFD4 workflow phase",
        cancellation="stop the scheduler-owned attempt child",
        sequence_structure_mode="co_design",
        downstream_validation=("refolding.optimized_prediction", "analysis.ipsae"),
        gpu_required=True, minimum_vram_gb=16.0,
        supported_stopping_policies=("validated", "diversity_qualified"),
    ),
    "rfdiffusion": GeneratorCapabilities(
        engine_id="rfdiffusion", label="RFdiffusion", campaign_adapter=False,
        independent_attempts=True, attempt_unit="one backbone trajectory",
        attempt_budget_semantics="trajectory count and generated sequence count are distinct",
        native_acceptance_available=False,
        native_acceptance_semantics="backbone generation has no sequence-level binder acceptance",
        candidate_output="backbone-only structures",
        output_normalization="existing generator normalizer; campaign handoff not yet verified",
        seed_assignment="per-trajectory behavior requires adapter verification",
        progress_reporting="native process state; durable attempt progress not yet verified",
        cancellation="scheduler cancellation exists; candidate handoff not yet verified",
        sequence_structure_mode="backbone_only",
        downstream_validation=("protein_mpnn", "refolding.optimized_prediction", "analysis.ipsae"),
        gpu_required=True, minimum_vram_gb=None,
        supported_stopping_policies=(),
        deferred_reason="The campaign adapter must include real ProteinMPNN sequence assignment before validation.",
    ),
    "rfdiffusion3": GeneratorCapabilities(
        engine_id="rfdiffusion3", label="RFdiffusion3 / Foundry", campaign_adapter=False,
        independent_attempts=True, attempt_unit="one backbone trajectory",
        attempt_budget_semantics="native structure attempts and Foundry-MPNN sequence outputs differ",
        native_acceptance_available=False,
        native_acceptance_semantics="generation or MPNN prefilter outputs are not campaign acceptance",
        candidate_output="backbone structures and optional Foundry-MPNN sequences",
        output_normalization="existing generator normalizer; campaign attempt mapping not yet verified",
        seed_assignment="per-trajectory behavior requires adapter verification",
        progress_reporting="native process state; durable attempt progress not yet verified",
        cancellation="scheduler cancellation exists; candidate handoff not yet verified",
        sequence_structure_mode="backbone_only",
        downstream_validation=("protein_mpnn", "refolding.optimized_prediction", "analysis.ipsae"),
        gpu_required=True, minimum_vram_gb=None,
        supported_stopping_policies=(),
        deferred_reason="Attempt-to-sequence accounting and explicit ProteinMPNN handoff need verification.",
    ),
    "boltzgen": GeneratorCapabilities(
        engine_id="boltzgen", label="BoltzGen", campaign_adapter=False,
        independent_attempts=False, attempt_unit="unverified",
        attempt_budget_semantics="native budget and final/intermediate output accounting are not verified",
        native_acceptance_available=None, native_acceptance_semantics="unverified",
        candidate_output="ranked and intermediate structures",
        output_normalization="existing generator normalizer; campaign semantics unverified",
        seed_assignment="unverified", progress_reporting="unverified", cancellation="scheduler-level only",
        sequence_structure_mode="unverified", downstream_validation=(), gpu_required=True,
        minimum_vram_gb=None, supported_stopping_policies=(),
        deferred_reason="Native attempt and candidate semantics remain unverified.",
    ),
    "pxdesign": GeneratorCapabilities(
        engine_id="pxdesign", label="PXDesign", campaign_adapter=False,
        independent_attempts=False, attempt_unit="unverified",
        attempt_budget_semantics="native run and output-candidate accounting are not verified",
        native_acceptance_available=None, native_acceptance_semantics="unverified",
        candidate_output="PXDesign prediction structures",
        output_normalization="existing generator normalizer; campaign semantics unverified",
        seed_assignment="unverified", progress_reporting="unverified", cancellation="scheduler-level only",
        sequence_structure_mode="unverified", downstream_validation=(), gpu_required=True,
        minimum_vram_gb=None, supported_stopping_policies=(),
        deferred_reason="Native attempt and candidate semantics remain unverified.",
    ),
    "genie3": GeneratorCapabilities(
        engine_id="genie3", label="Genie3", campaign_adapter=False,
        independent_attempts=False, attempt_unit="unverified",
        attempt_budget_semantics="native sampling and output-candidate accounting are not verified",
        native_acceptance_available=None, native_acceptance_semantics="unverified",
        candidate_output="generated structures and native evaluation rows",
        output_normalization="existing generator normalizer; campaign semantics unverified",
        seed_assignment="unverified", progress_reporting="unverified", cancellation="scheduler-level only",
        sequence_structure_mode="unverified", downstream_validation=(), gpu_required=True,
        minimum_vram_gb=None, supported_stopping_policies=(),
        deferred_reason="Native attempt, seed, and candidate semantics remain unverified.",
    ),
}


def generator_capabilities(engine_id: str) -> GeneratorCapabilities:
    key = str(engine_id or "").strip().lower()
    try:
        return _CAPABILITIES[key]
    except KeyError as exc:
        raise ValueError(f"Unknown campaign generator: {engine_id!r}.") from exc


def list_generator_capabilities() -> list[dict[str, object]]:
    return [_CAPABILITIES[key].as_dict() for key in sorted(_CAPABILITIES)]


def validate_stopping_policy(
    engine_id: str,
    acceptance_level: str,
    validation_stages: list[dict[str, object]],
    *,
    has_validation_criteria: bool,
) -> None:
    capability = generator_capabilities(engine_id)
    if not capability.campaign_adapter:
        reason = capability.deferred_reason or "No verified distributed-campaign adapter is available."
        raise ValueError(f"{capability.label} cannot run an acceptance-driven campaign: {reason}")
    if acceptance_level not in capability.supported_stopping_policies:
        supported = ", ".join(capability.supported_stopping_policies) or "none"
        raise ValueError(
            f"{capability.label} does not support {acceptance_level!r} stopping; supported policies: {supported}."
        )
    if not capability.native_acceptance_available and (
        acceptance_level == "native" or not validation_stages or not has_validation_criteria
    ):
        raise ValueError(
            f"{capability.label} has no authoritative native acceptance filter; configure downstream "
            "validation stages with explicit measurements or scientific rules."
        )
