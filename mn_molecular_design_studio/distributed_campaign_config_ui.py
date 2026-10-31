"""Pipeline tile controls for the existing distributed campaign adapter."""

from __future__ import annotations

import json
from typing import Any


_GENERATORS = {
    "bindcraft2": "bindcraft2",
    "rfd4": "rfd4_proteina",
}


def _default_validation_stages() -> list[dict[str, Any]]:
    """Use the retained PD-L1 Boltz-2 measurement policy as an editable example."""
    return [{
        "workflow": "refolding.optimized_prediction",
        "parameters": {
            "engine": "boltz2",
            "engine_options": {
                "sampling_steps": 20,
                "use_msa": False,
                "use_msa_server": False,
                "use_potentials": False,
            },
            "num_recycles": 1,
            "optimization_mode": "off",
            "samples": 1,
            "target_use_msa": False,
            "use_target_template": False,
        },
        "required_metrics": ["iptm", "ptm", "plddt_mean"],
        "resources": {"cpu_cores": 8, "gpu_device": "auto", "min_vram_gb": 16.0},
    }]


def _json_editor(st: Any, label: str, value: object, *, key: str, height: int) -> tuple[object, str | None]:
    text = st.text_area(label, value=json.dumps(value, indent=2), key=key, height=height)
    try:
        return json.loads(text), None
    except json.JSONDecodeError as exc:
        return None, f"{label} must be valid JSON (line {exc.lineno}, column {exc.colno})."


def render_campaign_config_editor(
    st: Any,
    node_kind: str,
    current: dict[str, Any] | None,
    *,
    key: str,
    seed: int = 0,
    gpu_inventory: list[dict[str, Any]] | None = None,
    generation_cpu_cores: int | None = None,
) -> tuple[dict[str, Any] | None, bool]:
    """Render supported campaign settings and return a validated API config."""
    generator_id = _GENERATORS.get(node_kind)
    if generator_id is None:
        return current, True

    existing = current if isinstance(current, dict) else {}
    enabled = st.checkbox(
        "Use distributed campaign",
        value=bool(existing.get("enabled", current is not None)),
        key=key + "-enabled",
        help="Use the existing controller to run to a generated count, passing count, or time limit across eligible workers.",
    )
    if not enabled:
        return None, True

    from mn_molecular_design_studio.engine.workflow_engine.workflows.distributed_campaign import validate_campaign_config
    from mn_molecular_design_studio.engine.workflow_engine.workflows.generator_capabilities import generator_capabilities

    capability = generator_capabilities(generator_id)
    default_level = "validated" if not capability.native_acceptance_available else "native"
    from mn_molecular_design_studio.campaign_objectives import OBJECTIVE_LABELS
    objective_options = list(OBJECTIVE_LABELS)
    current_objective = str(existing.get("objective_type") or "generated")
    if current_objective not in objective_options:
        current_objective = "generated"
    objective_type = st.selectbox(
        "Campaign objective", objective_options, index=objective_options.index(current_objective),
        format_func=OBJECTIVE_LABELS.get,
        key=key + "-objective-type",
    )
    accepted_target = None if objective_type == "time" else int(st.number_input(
        "Target designs" if objective_type == "generated" else "Target unique passing designs",
        min_value=1, value=int(existing.get("accepted_target") or 1),
        step=1, key=key + "-accepted-target",
        help=(
            "Studio-level campaign objective: total generated designs, or unique candidates that pass the configured criteria. "
            "It is independent of the generator's internal GPU/dataloader batch size."
        ),
    ))
    max_attempts = int(st.number_input(
        "Maximum total designs", min_value=1, value=int(existing.get("max_attempts") or max(2, accepted_target or 1000)),
        step=1, key=key + "-max-attempts",
        help="Hard cap on requested design outputs over the whole campaign, including rejected designs. This is independent of native jobs and GPU microbatches.",
    ))
    native_designs_per_job = 1
    if generator_id == "rfd4_proteina":
        try:
            requested_default = int(existing.get("native_designs_per_job") or 50)
        except (TypeError, ValueError):
            requested_default = 50
        native_designs_per_job = int(st.number_input(
            "RFD4 designs per distributed job", min_value=1,
            max_value=max(1, max_attempts), value=min(max_attempts, max(1, requested_default)),
            step=1, key=key + "-native-designs-per-job",
            help="Sets RFD4 num_designs for each scheduler job. RFD4's own GPU batch_size stays in the RFD4 settings and controls its internal microbatching.",
        ))
    from mn_molecular_design_studio.campaign_objectives_ui import render_time_limits
    time_limits, time_ready = render_time_limits(st, existing, key=key, required=objective_type == "time")
    compute_seconds = float(st.number_input(
        "Maximum GPU seconds (0 = no limit)", min_value=0.0,
        value=float(existing.get("max_compute_seconds") or 0), key=key + "-compute-budget",
    ))
    policy_options = list(capability.supported_stopping_policies)
    current_level = str(existing.get("acceptance_level") or default_level)
    if current_level not in policy_options:
        current_level = default_level
    if objective_type == "accepted":
        acceptance_level = st.selectbox(
            "Scientific acceptance policy", policy_options,
            index=policy_options.index(current_level), key=key + "-acceptance-level",
            help="Validated acceptance requires explicit valid measurements or scientific rules. RFD4 has no authoritative native acceptance flag.",
        )
    else:
        acceptance_level = current_level
        st.caption("Generated designs are counted from recorded candidate outputs, independent of native or comparative acceptance.")

    try:
        legacy_batch_size = max(1, int(existing.get("batch_size") or 1))
    except (TypeError, ValueError):
        legacy_batch_size = 1
    try:
        legacy_concurrent_batches = max(1, int(existing.get("max_concurrent_batches") or 1))
    except (TypeError, ValueError):
        legacy_concurrent_batches = 1
    if generator_id == "rfd4_proteina":
        # RFD4's real per-job output count is shown explicitly below; retain
        # the legacy Studio grouping field at one to avoid a second "batch size".
        legacy_batch_size = 1
        legacy_concurrent_batches = 1
    legacy_concurrency = legacy_batch_size * legacy_concurrent_batches
    try:
        generation_job_default = max(1, int(
            existing.get("max_concurrent_generation_jobs")
            or (4 if generator_id == "rfd4_proteina" else legacy_concurrency)
        ))
    except (TypeError, ValueError):
        generation_job_default = legacy_concurrency
    max_concurrent_generation_jobs = int(st.number_input(
        "Maximum concurrent generation jobs", min_value=1,
        value=generation_job_default, step=1,
        key=key + "-generation-job-concurrency",
        help="Caps independently scheduled native jobs across compatible workstations. One RFD4 job can contain multiple designs; each job keeps its own seed, receipt and recovery record.",
    ))
    max_concurrent_validations = int(st.number_input(
        "Concurrent downstream validations", min_value=1,
        value=int(existing.get("max_concurrent_validations") or 2), step=1,
        key=key + "-validation-concurrency",
    ))

    placement_options = ["local", "auto", "eligible"]
    placement = str(existing.get("workstation_selection") or "auto")
    if placement not in placement_options:
        placement = "auto"
    workstation_selection = st.selectbox(
        "Worker workstation placement", placement_options,
        index=placement_options.index(placement), key=key + "-workstation-placement",
        format_func={"local": "This workstation", "auto": "Any available compatible workstation", "eligible": "All eligible workstations"}.get,
    )
    gpu_selection = st.text_input(
        "GPU selection per workstation", value=str(existing.get("gpu_selection") or "auto"),
        key=key + "-gpu-selection",
        help="Use auto, eligible, all, or explicit device IDs such as 0 or 0,1.",
    ).strip() or "auto"

    stages = existing.get("validation_stages")
    stages = stages if isinstance(stages, list) else []
    if objective_type in {"generated", "time"}:
        stages = []
    if objective_type == "accepted" and (acceptance_level != "native" or not capability.native_acceptance_available):
        if not stages:
            stages = _default_validation_stages()
        stages, stage_error = _json_editor(
            st, "Validation stages (JSON)", stages, key=key + "-validation-stages", height=260,
        )
        st.caption("The default retained-evidence example accepts finite, valid iPTM, pTM and mean pLDDT measurements; it sets no score cutoffs. Add explicit rules here when the campaign needs score thresholds.")
        st.caption("Each stage may set match to all or any for its threshold rules, and missing_data to reject, unresolved or ignore. Required metrics remain prerequisites. Different stages must all pass.")
    else:
        stage_error = None

    diversity = existing.get("diversity") if isinstance(existing.get("diversity"), dict) and objective_type == "accepted" else {}
    if objective_type == "accepted" and acceptance_level == "diversity_qualified":
        diversity, diversity_error = _json_editor(
            st, "Diversity rules (JSON)", diversity, key=key + "-diversity", height=130,
        )
        st.caption("An empty diversity object uses the controller defaults: exact-sequence deduplication and one candidate per diversity cluster.")
    else:
        diversity_error = None
    if objective_type == "accepted" and acceptance_level == "native":
        st.caption("Global stopping counts the generator's authoritative native acceptance decisions. Prediction measurements remain separate evidence.")

    config = {
        "enabled": True,
        "generator_id": generator_id,
        "objective_type": objective_type,
        "unique_accepted": objective_type == "accepted",
        "accepted_target": accepted_target,
        "max_attempts": max_attempts,
        "native_designs_per_job": native_designs_per_job,
        "acceptance_level": acceptance_level,
        # Keep legacy dispatch grouping values intact for persisted configs;
        # the explicit job cap now controls actual campaign parallelism.
        "batch_size": legacy_batch_size,
        "max_concurrent_batches": legacy_concurrent_batches,
        "max_concurrent_generation_jobs": max_concurrent_generation_jobs,
        "max_concurrent_validations": max_concurrent_validations,
        "workstation_selection": workstation_selection,
        "gpu_selection": gpu_selection,
        "validation_stages": stages,
        "diversity": diversity,
        "seed": int(existing.get("seed", seed) or 0),
    }
    config.update(time_limits)
    if compute_seconds:
        config["max_compute_seconds"] = compute_seconds
    ready = time_ready
    for error in (stage_error, diversity_error):
        if error:
            st.error(error)
            ready = False
    if ready:
        try:
            config = validate_campaign_config(config)
        except (TypeError, ValueError) as exc:
            st.error(f"Campaign settings are invalid: {exc}")
            ready = False

    st.caption("Available GPUs are candidates for scheduling; allocated GPUs are selected only when each attempt or validation starts. Validation concurrency is independent of generation concurrency.")
    if generator_id == "rfd4_proteina":
        st.caption(
            f"Up to {native_designs_per_job * max_concurrent_generation_jobs:,} RFD4 designs can be in flight "
            f"({max_concurrent_generation_jobs} jobs × {native_designs_per_job} designs); "
            "RFD4 GPU batch_size controls the internal microbatch within each job."
        )
    st.caption(
        ("Studio timed campaign; " if objective_type == "time" else
         f"Studio campaign target: {accepted_target} {'generated designs' if objective_type == 'generated' else 'unique passing designs'}; ")
        +
        f"maximum requested designs: {max_attempts}. "
        + (f"Each RFD4 execution job requests {native_designs_per_job} native designs (num_designs); the native GPU batch_size remains independently configured. "
           if generator_id == "rfd4_proteina" else
           "Each execution job requests one native design. ")
        +
        f"Maximum concurrent generation jobs: {max_concurrent_generation_jobs}; "
        f"maximum concurrent validation jobs: {max_concurrent_validations}. "
        "The controller dispatches available jobs toward the Studio objective and hard design cap. "
        "When the passing target is reached, queued jobs stop and already-running jobs finish."
    )
    st.markdown("**Effective campaign settings before submission**")
    st.code(json.dumps(config, indent=2, sort_keys=True), language="json")
    if node_kind == "rfd4":
        st.caption("Finished RFD4 native batches are reconciled candidate-by-candidate and sent through the configured acceptance stages while other batches continue. Reaching the passing target stops new dispatch; already-running native batches finish.")
    with st.expander("Native execution plan"):
        from mn_molecular_design_studio.engine.workflow_engine.workflows.execution_planning import (
            PROPOSED_NATIVE_CHOICES, plan_campaign_execution,
        )
        st.caption("Preview only: these controls are not saved to the tile or sent for execution. Unsupported choices remain disabled for campaign submission.")
        preview_count = st.selectbox(
            "Hypothetical native attempts per job", list(PROPOSED_NATIVE_CHOICES),
            key=key + "-plan-native-count",
        )
        preview_inputs = int(st.number_input(
            "Hypothetical input count", min_value=1, value=1, step=1, key=key + "-plan-input-count",
        ))
        try:
            plan = plan_campaign_execution(
                config, native_attempts_per_job=preview_count, input_count=preview_inputs,
                eligible_gpus=gpu_inventory, generation_cpu_cores=generation_cpu_cores,
            )
            for issue in plan["errors"]:
                st.warning(issue)
            if gpu_inventory is None:
                st.caption("GPU inventory was not supplied to this preview. Workstation/GPU selectors and resource requirements are shown; exact eligibility remains unresolved until authenticated placement checks.")
            st.json(plan)
        except (TypeError, ValueError) as exc:
            st.warning(f"Execution plan unavailable: {exc}")
    return config, ready
