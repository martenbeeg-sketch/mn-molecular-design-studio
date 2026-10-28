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
) -> tuple[dict[str, Any] | None, bool]:
    """Render supported campaign settings and return a validated API config."""
    generator_id = _GENERATORS.get(node_kind)
    if generator_id is None:
        return current, True

    existing = current if isinstance(current, dict) else {}
    enabled = st.checkbox(
        "Use distributed acceptance campaign",
        value=bool(existing.get("enabled", current is not None)),
        key=key + "-enabled",
        help="One controller tracks the global accepted target and attempt budget across eligible workers.",
    )
    if not enabled:
        return None, True

    from mn_molecular_design_studio.engine.workflow_engine.workflows.distributed_campaign import validate_campaign_config
    from mn_molecular_design_studio.engine.workflow_engine.workflows.generator_capabilities import generator_capabilities

    capability = generator_capabilities(generator_id)
    default_level = "validated" if not capability.native_acceptance_available else "native"
    accepted_target = int(st.number_input(
        "Global accepted target", min_value=1, value=int(existing.get("accepted_target") or 1),
        step=1, key=key + "-accepted-target",
    ))
    max_attempts = int(st.number_input(
        "Global maximum native attempts", min_value=1, value=int(existing.get("max_attempts") or max(2, accepted_target)),
        step=1, key=key + "-max-attempts",
    ))
    policy_options = list(capability.supported_stopping_policies)
    current_level = str(existing.get("acceptance_level") or default_level)
    if current_level not in policy_options:
        current_level = default_level
    acceptance_level = st.selectbox(
        "Scientific acceptance policy", policy_options,
        index=policy_options.index(current_level), key=key + "-acceptance-level",
        help="Validated acceptance requires explicit valid measurements or scientific rules. RFD4 has no authoritative native acceptance flag.",
    )

    left, right = st.columns(2)
    controller_batch_size = int(left.number_input(
        "Controller dispatch batch size", min_value=1,
        value=int(existing.get("batch_size") or 1), step=1,
        key=key + "-controller-batch-size",
        help="Legacy dispatch grouping. Each attempt currently launches a separate execution job and container.",
    ))
    max_concurrent_batches = int(right.number_input(
        "Concurrent generation batches", min_value=1,
        value=int(existing.get("max_concurrent_batches") or 1), step=1,
        key=key + "-generation-concurrency",
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
    if acceptance_level != "native" or not capability.native_acceptance_available:
        if not stages:
            stages = _default_validation_stages()
        stages, stage_error = _json_editor(
            st, "Validation stages (JSON)", stages, key=key + "-validation-stages", height=260,
        )
        st.caption("The default retained-evidence example accepts finite, valid iPTM, pTM and mean pLDDT measurements; it sets no score cutoffs. Add explicit rules here when the campaign needs score thresholds.")
    else:
        stage_error = None

    diversity = existing.get("diversity") if isinstance(existing.get("diversity"), dict) else {}
    if acceptance_level == "diversity_qualified":
        diversity, diversity_error = _json_editor(
            st, "Diversity rules (JSON)", diversity, key=key + "-diversity", height=130,
        )
        st.caption("An empty diversity object uses the controller defaults: exact-sequence deduplication and one candidate per diversity cluster.")
    else:
        diversity_error = None
    if acceptance_level == "native":
        st.caption("Global stopping counts the generator's authoritative native acceptance decisions. Prediction measurements remain separate evidence.")

    config = {
        "enabled": True,
        "generator_id": generator_id,
        "accepted_target": accepted_target,
        "max_attempts": max_attempts,
        "acceptance_level": acceptance_level,
        "batch_size": controller_batch_size,
        "max_concurrent_batches": max_concurrent_batches,
        "max_concurrent_validations": max_concurrent_validations,
        "workstation_selection": workstation_selection,
        "gpu_selection": gpu_selection,
        "validation_stages": stages,
        "diversity": diversity,
        "seed": int(existing.get("seed", seed) or 0),
    }
    ready = True
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
    st.caption(
        f"Native attempts per execution job: 1 (the only verified campaign choice). "
        f"Maximum concurrent generation jobs: {controller_batch_size * max_concurrent_batches}; "
        f"maximum concurrent validation jobs: {max_concurrent_validations}. "
        "Generation capacity is the product of the two legacy dispatch controls, subject to the global target, attempt budget and available resources."
    )
    st.markdown("**Effective campaign settings before submission**")
    st.code(json.dumps(config, indent=2, sort_keys=True), language="json")
    if node_kind == "rfd4":
        st.caption("Distributed RFD4 children request one native design with dataloader batch size 1. Controller dispatch size changes only how the campaign groups separate execution jobs.")
    return config, ready
