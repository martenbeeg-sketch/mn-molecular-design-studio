"""Deterministic read-only campaign arithmetic. No scheduler, jobs or GPU probing.

A plan is descriptive data, never a submission request or authorization manifest.
"""
from __future__ import annotations

from copy import deepcopy
import math
from .generator_capabilities import generator_capabilities

PROPOSED_NATIVE_CHOICES = (1, 5, 10, 25, 50, 100)


def _positive(value, name):
    if type(value) is not int or value < 1:
        raise ValueError(f"{name} must be a positive integer.")
    return value


def _ceil(count, size):
    return (count + size - 1) // size


def describe_native_invocation(engine_id, *, attempts_per_input=1, input_count=1, native_batch_size=1):
    """Describe current ordinary native counts, independently of campaign support.

    Requested outputs are not guaranteed usable candidates or acceptances. The
    exact RFD3 remainder alternative is proposed arithmetic, not implemented code.
    """
    capability = generator_capabilities(engine_id)
    count = _positive(attempts_per_input, "attempts_per_input")
    inputs = _positive(input_count, "input_count")
    batch = _positive(native_batch_size, "native_batch_size")
    key = capability.engine_id
    requested = count * inputs
    result = {"requested_native_attempts": requested, "input_count": inputs,
              "attempts_per_input": count, "native_batch_size": batch,
              "native_attempts": None, "nominal_native_outputs": None,
              "native_python_invocations": None, "rounding_extra_attempts": None,
              "exact_remainder_alternative": None, "notes": []}
    if key == "rfd4_proteina":
        chunks, remainder = divmod(count, 100)
        result.update(native_attempts=requested, nominal_native_outputs=requested,
                      native_python_invocations=chunks + bool(remainder), rounding_extra_attempts=0,
                      wrapper_chunks={"full_chunks": chunks, "full_chunk_replicates": 100,
                                      "remainder_replicates": remainder},
                      loader_batches=chunks * _ceil(100 * inputs, batch) + (_ceil(remainder * inputs, batch) if remainder else 0))
        result["notes"].append("Replicates multiply each input; dataloader final batches do not authorize extra designs. Search/oracle multiplicity is outside this simple no-search count.")
    elif key == "rfdiffusion3":
        effective_batch = min(batch, count)
        full, remainder = divmod(count, effective_batch)
        actual = inputs * _ceil(count, effective_batch) * effective_batch
        result.update(native_attempts=actual, nominal_native_outputs=actual,
                      native_python_invocations=1, rounding_extra_attempts=actual-requested,
                      effective_diffusion_batch_size=effective_batch, n_batches=_ceil(count, effective_batch),
                      exact_remainder_alternative={
                          "implemented": False, "full_batches_per_input": full,
                          "full_batch_size": effective_batch, "remainder_per_input": remainder,
                          "native_python_invocations": 1 + bool(remainder), "native_attempts": requested})
        result["notes"].append("Current tile rounds native outputs up, then selects requested outputs. Campaign execution is disabled; no output may be omitted from future attempt accounting.")
    elif key in {"rfdiffusion", "bindcraft2"}:
        result.update(native_attempts=requested, native_python_invocations=inputs, rounding_extra_attempts=0)
        if key == "rfdiffusion":
            result["nominal_native_outputs"] = requested
        else:
            result["notes"].append("Trajectory budget is an upper bound; native accepted-target stopping can end early. Sequence-level candidate multiplicity is unknown.")
        result["notes"].append("One input per invocation; several inputs require independent invocations. Native GPU batch size is not this engine's trajectory count.")
    else:
        result["notes"].append("Native count semantics have not been verified for this engine.")
    return result


def _partition_budget(budget, inputs, size):
    """One input per hypothetical job; allocate a global budget evenly by input.

    Compact histograms avoid building one record per attempt for large campaigns.
    No input is multiplied by the global budget. Extra attempts go to the earliest
    input ordinals. This is a planning convention, not a new campaign adapter.
    """
    base, extras = divmod(budget, inputs)
    groups = {}
    for per_input, number in ((base, inputs-extras), (base+1, extras)):
        full, remainder = divmod(per_input, size)
        if full and number:
            groups[size] = groups.get(size, 0) + full * number
        if remainder and number:
            groups[remainder] = groups.get(remainder, 0) + number
    return {"base_attempts_per_input": base, "first_inputs_with_one_extra": extras,
            "job_groups": [{"native_attempts_requested_per_job": k, "execution_jobs": groups[k]} for k in sorted(groups)],
            "execution_jobs_if_budget_exhausted": sum(groups.values()),
            "total_requested_attempts": sum(k*v for k,v in groups.items())}


def _gpu_snapshot(inventory, config):
    selector = str(config.get("gpu_selection", "auto")).removeprefix("device=")
    device_ids = None if selector in {"all", "auto", "eligible"} else set(selector.replace(" ", "").split(","))
    resources = config.get("resources") or {}
    minimum = resources.get("min_vram_gb", 0)
    if inventory is not None and (not isinstance(inventory, list) or any(not isinstance(x, dict) for x in inventory)):
        raise ValueError("eligible_gpus must be a list of inventory objects.")
    rows = []
    for item in inventory or []:
        row = deepcopy(item)
        excluded, unresolved = [], []
        if row.get("snapshot_current") is not True or not row.get("observed_at"):
            unresolved.append("Snapshot freshness/timestamp is unverified.")
        if config.get("workstation_selection", "local") == "local":
            if row.get("is_local") is False:
                excluded.append("Campaign placement is local.")
            elif row.get("is_local") is not True:
                unresolved.append("Local workstation identity is unknown.")
        if device_ids is not None and str(row.get("gpu_id")) not in device_ids:
            excluded.append("GPU is outside the selected device IDs.")
        for key, label in (("online", "Workstation is online"), ("compatible", "Exact runtime/image/checkpoint compatibility is verified")):
            if row.get(key) is False:
                excluded.append(label + ": no.")
            elif row.get(key) is not True:
                unresolved.append(label + ": unknown.")
        for field, required, label in [("total_vram_gb", minimum, "VRAM"), *[
            (field, resources.get(setting, 0), setting) for field,setting in
            [("cpu_cores_available", "cpu_cores"), ("ram_gb_available", "ram_gb"), ("scratch_gb_available", "scratch_gb")]
        ]]:
            if required:
                reported = row.get(field)
                if type(reported) not in {int, float} or not math.isfinite(reported):
                    unresolved.append(label + " capacity is unknown.")
                elif reported < required:
                    excluded.append(label + " capacity is below the request in this snapshot.")
        if not row.get("node_id") or row.get("gpu_id") is None:
            unresolved.append("Physical workstation/GPU identity is incomplete.")
        free = row.get("free_vram_gb")
        free_known = type(free) in {int, float} and math.isfinite(free)
        immediate = (not excluded and not unresolved and row.get("available") is True
                     and free_known and free >= minimum)
        row.update(eligibility="excluded" if excluded else "unresolved" if unresolved else "eligible_in_supplied_snapshot",
                   available_for_immediate_start_in_snapshot=bool(immediate),
                   availability_explanation="Explicit free-device and VRAM evidence is required for immediate-start availability; eligibility alone can still require queuing.",
                   reasons=excluded+unresolved)
        rows.append(row)
    return rows


def _manifest_preview(config, budget, inputs, size, resources):
    """Bounded deterministic membership examples, with no durable authorization."""
    from .execution_job_ledger import digest, manifest_entries
    seed = config.get("seed", 0)
    if type(seed) is not int or seed < 0:
        return {"available": False, "reason": "Invalid campaign seed.", "authorization": False}
    namespace = "preview-" + digest({"config": config, "inputs": inputs, "size": size})[:16]
    base, extras = divmod(budget, inputs)
    total_jobs = _partition_budget(budget, inputs, size)["execution_jobs_if_budget_exhausted"]
    jobs, ordinal, shown = [], 1, 0
    for input_index in range(min(inputs, budget)):
        remaining = base + int(input_index < extras)
        while remaining:
            if len(jobs) >= 16 or shown >= 100:
                return {"available": True, "authorization": False, "namespace": namespace, "jobs": jobs,
                        "total_jobs": total_jobs, "total_attempts": budget, "shown_attempts": shown,
                        "truncated": True, "omitted_jobs": total_jobs - len(jobs), "omitted_attempts": budget - shown}
            count = min(size, remaining)
            shown_count = min(count, 100 - shown)
            entries = manifest_entries(namespace, seed, ordinal,
                                       [{"input_index": input_index, "identity_verified": False}] * shown_count)
            jobs.append({"execution_job_id": f"{namespace}:j{ordinal:08d}", "authorization": False,
                         "proposed_attempt_count": count, "manifest": entries,
                         "manifest_digest": digest(entries) if shown_count == count else None,
                         "manifest_complete": shown_count == count, "omitted_members": count - shown_count,
                         "configuration_digest": digest({"campaign": config, "native": {}, "resource_requirements": resources}),
                         "native_configuration_resolved": False,
                         "resource_requirements": deepcopy(resources)})
            remaining -= count
            ordinal += count
            shown += shown_count
    return {"available": True, "authorization": False, "namespace": namespace, "jobs": jobs,
            "total_jobs": total_jobs, "total_attempts": budget, "shown_attempts": shown,
            "truncated": shown < budget, "omitted_jobs": 0, "omitted_attempts": budget - shown}


def plan_campaign_execution(campaign, *, engine_id=None, native_attempts_per_job=None,
                            input_count=1, native_batch_size=1, designs_per_input=1, eligible_gpus=None,
                            generation_cpu_cores=None):
    """Plan fresh-campaign budget exhaustion; never predicts acceptance or dispatches.

    Inventory is caller-supplied, timestamped evidence, not live scheduler authority.
    Unsupported settings retain their requested arithmetic and have no supported
    execution plan. Legacy serialized campaign config is never mutated.
    """
    if not isinstance(campaign, dict):
        raise ValueError("campaign must be an object.")
    raw = deepcopy(campaign)
    engine = engine_id or raw.get("generator_id", "bindcraft2")
    capability = generator_capabilities(engine)
    if engine_id is not None and generator_capabilities(raw.get("generator_id", engine_id)).engine_id != capability.engine_id:
        raise ValueError("engine_id conflicts with campaign.generator_id.")
    engine = capability.engine_id
    inputs = _positive(input_count, "input_count")
    batch = _positive(native_batch_size, "native_batch_size")
    requested = raw.get("native_attempts_per_job", 1) if native_attempts_per_job is None else native_attempts_per_job
    if native_attempts_per_job is not None and "native_attempts_per_job" in raw and (type(raw["native_attempts_per_job"]) is not type(requested) or raw["native_attempts_per_job"] != requested):
        raise ValueError("Conflicting requested native attempts per job.")
    requested = _positive(requested, "native_attempts_per_job")
    budget = _positive(raw.get("max_attempts", 1), "max_attempts")
    target = _positive(raw.get("accepted_target", 1), "accepted_target")
    if not isinstance(raw.get("resources", {}), dict):
        raise ValueError("resources must be an object.")
    for key in ("cpu_cores", "min_vram_gb", "ram_gb", "scratch_gb"):
        value = raw.get("resources", {}).get(key)
        if key in raw.get("resources", {}):
            if key == "cpu_cores":
                _positive(value, "resources.cpu_cores")
            elif type(value) not in {int, float} or not math.isfinite(value) or value < 0:
                raise ValueError(f"resources.{key} must be finite and nonnegative.")
    stages = raw.get("validation_stages", [])
    if not isinstance(stages, list) or any(not isinstance(x, dict) for x in stages):
        raise ValueError("validation_stages must be a list of objects.")
    if any(not isinstance(stage.get("resources", {}), dict) for stage in stages):
        raise ValueError("validation stage resources must be objects.")
    errors = []
    # Validate all existing acceptance/resource/retry fields through the unchanged
    # controller validator. Remove only the separately reported hypothetical count.
    from .distributed_campaign import validate_campaign_config
    validation_input = {**raw, "generator_id": engine}
    validation_input.pop("native_attempts_per_job", None)
    try:
        normalized = validate_campaign_config(validation_input)
    except ValueError as exc:
        errors.append(str(exc))
        normalized = validation_input
    choices = [1] if capability.campaign_adapter else []
    if requested not in choices:
        errors.append(f"{capability.label} verified campaign choices are {choices}; {requested} remains descriptive only.")
    if inputs != 1:
        errors.append("Current campaign adapters require one input; multiple-input partitioning is hypothetical only.")
    if batch != 1:
        errors.append("Current campaign jobs do not support changing the native GPU microbatch from 1.")
    if target > budget and not any("accepted_target cannot exceed" in x for x in errors):
        errors.append("accepted_target cannot exceed max_attempts under current campaign semantics.")
    partition = _partition_budget(budget, inputs, requested)
    physical, nominal = 0, 0
    for group in partition["job_groups"]:
        counts = describe_native_invocation(engine, attempts_per_input=group["native_attempts_requested_per_job"], native_batch_size=batch)
        group["native_count_semantics"] = counts
        physical = None if physical is None or counts["native_attempts"] is None else physical + group["execution_jobs"] * counts["native_attempts"]
        nominal = None if nominal is None or counts["nominal_native_outputs"] is None else nominal + group["execution_jobs"] * counts["nominal_native_outputs"]
    if physical is not None and physical > budget:
        errors.append(f"Native rounding would request {physical} physical attempts against a global budget of {budget}.")
    legacy_batch = _positive(normalized.get("batch_size", normalized.get("dispatch_batch_size", 1)), "dispatch_batch_size")
    legacy_concurrency = _positive(normalized.get("max_concurrent_batches", 1), "max_concurrent_batches")
    generation = legacy_batch * legacy_concurrency
    if "max_concurrent_generation_jobs" in raw:
        explicit = _positive(raw["max_concurrent_generation_jobs"], "max_concurrent_generation_jobs")
        errors.append("max_concurrent_generation_jobs is a proposed API field; execution still uses the legacy dispatch product.")
    else:
        explicit = None
    validations = _positive(normalized.get("max_concurrent_validations", 2), "max_concurrent_validations")
    generation_resources = deepcopy(normalized.get("resources", {}))
    if capability.campaign_adapter:
        parent_cpu = _positive(generation_cpu_cores, "generation_cpu_cores") if generation_cpu_cores is not None else (8 if engine == "rfd4_proteina" else 4)
        generation_resources = {"cpu_cores": parent_cpu, "min_vram_gb": capability.minimum_vram_gb or 0.0,
                                "ram_gb": 0.0, "scratch_gb": 0.0, **generation_resources}
    validation_resources = []
    for stage in normalized.get("validation_stages", []):
        gpu = stage.get("workflow") == "refolding.optimized_prediction"
        validation_resources.append({"workflow": stage.get("workflow"), "cpu_cores": 4,
                                     "gpu_required": gpu, "gpu_device": "auto" if gpu else "cpu",
                                     **deepcopy(stage.get("resources", {}))})
    gpus = _gpu_snapshot(eligible_gpus, {**normalized, "resources": generation_resources})
    supported = not errors
    return {
        "schema_version": 1, "read_only": True, "submission_authorized": False,
        "campaign_configuration_supported": supported, "errors": errors,
        "engine": capability.as_dict(),
        "campaign": {"accepted_target": target, "max_native_attempts": budget,
                     "acceptance_level": normalized.get("acceptance_level", "native")},
        "native_attempts_per_job": {"requested": requested, "supported_campaign_choices": choices,
                                    "proposed_descriptive_choices": list(PROPOSED_NATIVE_CHOICES),
                                    "effective": requested if supported else None},
        "execution_job_manifest_preview": _manifest_preview(normalized, budget, inputs, requested, generation_resources),
        "hypothetical_requested_plan": {**partition, "input_count": inputs,
                                        "total_native_attempts_including_rounding": physical,
                                        "nominal_native_output_upper_bound": nominal},
        "supported_execution_plan": deepcopy(partition) if supported else None,
        "ordinary_native_count_example": describe_native_invocation(engine, attempts_per_input=designs_per_input,
                                                                    input_count=inputs, native_batch_size=batch),
        "execution": {"max_concurrent_generation_jobs": generation,
                      "requested_proposed_generation_jobs": explicit,
                      "generation_concurrency_source": "legacy batch_size * max_concurrent_batches",
                      "max_concurrent_validations": validations,
                      "gpus_per_generation_job": 1 if capability.campaign_adapter else None,
                      "generation_resources": generation_resources,
                      "validation_resources": validation_resources,
                      "workstation_selection": normalized.get("workstation_selection", "local"),
                      "gpu_selection": normalized.get("gpu_selection", "auto")},
        "gpu_inventory": gpus,
        "inventory_status": "supplied_snapshot_only" if eligible_gpus is not None else "unknown_not_probed",
        "notes": ["Job counts are upper bounds if the full attempt budget is consumed; acceptance stopping can reduce actual work.",
                  "Manifest previews use synthetic identities and unverified RNG/input mappings; they authorize no work. Multi-member jobs remain accounting-only.",
                  "Hypothetical groups use one input per job and partition the global attempt budget across inputs; this is not an execution adapter.",
                  "Nominal native outputs are not guaranteed normalized, accepted, validated or selected candidates.",
                  "GPU availability and all exact peer compatibility gates must be checked again by the existing authenticated scheduler at submission."]}


if __name__ == "__main__":
    import argparse
    import json
    from pathlib import Path
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("request", help="JSON object with campaign and optional planner keyword arguments")
    args = parser.parse_args()
    print(json.dumps(plan_campaign_execution(**json.loads(Path(args.request).read_text())), indent=2, sort_keys=True))
