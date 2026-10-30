"""Shared workstation controls; remote files are opened only on request."""
from __future__ import annotations

import hashlib
import json
import sqlite3
import sys
import uuid
import csv
import io
import time
import zipfile
from pathlib import Path, PurePosixPath

import pandas as pd
import streamlit as st

from . import config, progress as job_progress, transport


_TOOL_LABELS = {
    "visual_pipeline": "Visual pipeline",
    "optimized_prediction": "Structure Prediction",
    "bindcraft": "BindCraft",
    "bindcraft2": "BindCraft2",
    "bindcraft2_generator": "BindCraft2 Generator",
    "boltzgen": "BoltzGen",
    "boltzgen_generator": "BoltzGen Generator",
    "boltzgen_fold": "BoltzGen Fold",
    "proteina_complexa": "Proteina-Complexa",
    "proteina_complexa_generator": "Proteina-Complexa Generator",
    "pxdesign": "PXDesign",
    "pxdesign_generator": "PXDesign",
    "rfdiffusion": "RFdiffusion",
    "rfdiffusion_classic": "RFdiffusion",
    "rfdiffusion_generator": "RFdiffusion Generator",
    "rfdiffusion_beta_pairing": "RFdiffusion β-pairing",
    "kejia_idr": "Kejia IDR Binder",
    "peptide_trap": "Peptide Trap",
    "rfdiffusion3": "RFdiffusion3",
    "rfdiffusion3_foundry": "RFdiffusion3",
    "genie3": "Genie3",
    "genie3_generator": "Genie3 Generator",
    "protpardelle_1c": "Protpardelle-1C",
    "protpardelle_generator": "Protpardelle-1C Generator",
    "rfantibody": "RFantibody",
    "esmfold2_binder_design": "ESMFold2 Binder Design",
    "esmfold2_binder": "ESMFold2 Binder Design",
    "esmfold2_cofolding": "ESMFold2",
    "esmfold2_monomer": "ESMFold2",
    "esmfold2_complex_validation": "ESMFold2",
    "esmfold": "ESMFold",
    "af2_initial_guess": "AlphaFold2 Initial Guess",
    "alphafold2": "AlphaFold2",
    "af3_jax_cofolding": "ColabFold 2.0 preview",
    "alphafold3": "AlphaFold3 JAX",
    "af3_torch": "AlphaFold3 PyTorch",
    "openfold3": "OpenFold3",
    "rf3": "RosettaFold3",
    "protenix": "Protenix",
    "protenix_v1": "Protenix CLI",
    "boltz2": "Boltz-2",
    "boltz2_complex": "Boltz-2",
    "boltz2_monomer": "Boltz-2",
    "boltz2_initial_guess": "Boltz-2",
    "ligandmpnn": "ProteinMPNN",
    "abmpnn": "AbMPNN",
    "foundry_mpnn": "ProteinMPNN (Foundry)",
    "esm_if1": "ESM-IF1",
    "caliby": "Caliby",
    "boltzgen_inverse": "BoltzGen Inverse Folding",
    "scannet": "ScanNet",
    "pesto": "PeSTo",
    "surf2spot": "Surf2Spot",
    "masif_seed": "MaSIF-seed",
    "ipsae": "ipSAE Scoring",
    "rosetta_metrics": "Rosetta Interface Metrics",
}


def display_tool_label(row):
    tool = str(row.get("display_tool") or row.get("tool") or row.get("job_type") or "").strip()
    optimized_prefix = "optimized_prediction_"
    if tool.startswith(optimized_prefix):
        engine = tool[len(optimized_prefix):]
        try:
            from mn_molecular_design_studio.engine.workflow_engine.workflows.optimized_prediction import ENGINE_PROFILES
            profile = ENGINE_PROFILES.get(engine)
            if profile:
                return str(profile.get("label") or engine)
        except (ImportError, AttributeError):
            pass
    if "." in tool:
        try:
            from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import WORKFLOWS
            spec = WORKFLOWS.get(tool)
            if spec:
                tool = spec.tool
        except (ImportError, AttributeError):
            pass
    if tool in _TOOL_LABELS:
        return _TOOL_LABELS[tool]
    workflow_aliases = {
        "design.bindcraft2": "BindCraft2",
        "design.bindcraft2_generator": "BindCraft2 Generator",
        "design.boltzgen": "BoltzGen",
        "design.boltzgen_generator": "BoltzGen Generator",
        "design.proteina_complexa": "Proteina-Complexa",
        "design.proteina_complexa_generate": "Proteina-Complexa Generator",
        "design.pxdesign_generate": "PXDesign",
        "design.rfdiffusion_generate": "RFdiffusion Generator",
        "design.rfdiffusion_beta_pairing": "RFdiffusion β-pairing",
        "design.kejia_idr": "Kejia IDR Binder",
        "design.peptide_trap": "Peptide Trap",
        "design.rfdiffusion3": "RFdiffusion3",
        "design.genie3_generate": "Genie3 Generator",
        "design.protpardelle_generate": "Protpardelle-1C Generator",
        "sequence.mpnn": "ProteinMPNN",
        "sequence.ligandmpnn": "ProteinMPNN",
        "sequence.abmpnn": "AbMPNN",
        "sequence.foundry_mpnn": "ProteinMPNN (Foundry)",
        "sequence.esm_if1": "ESM-IF1",
        "sequence.caliby": "Caliby",
        "sequence.boltzgen_inverse": "BoltzGen Inverse Folding",
        "refolding.af3_jax_cofolding": "ColabFold 2.0 preview",
        "refolding.af2_initial_guess": "AlphaFold2 Initial Guess",
        "refolding.openfold3": "OpenFold3",
        "refolding.rf3": "RosettaFold3",
        "refolding.protenix": "Protenix",
        "refolding.boltzgen_fold": "BoltzGen Fold",
    }
    if tool in workflow_aliases:
        return workflow_aliases[tool]
    return tool.replace("_", " ").strip().title() or "Unknown tool"


def _request(path, token, payload=None, method="GET", **kwargs):
    return transport.request("/api/v1/cluster" + path, payload, token=token, method=method, **kwargs)


def _is_bindcraft2_job(row):
    identity = " ".join(str(row.get(key) or "") for key in ("tool", "job_type", "workflow")).casefold()
    return "bindcraft2" in identity.replace(" ", "") or "bindcraft 2" in identity


def _row_progress(row):
    progress = row.get("progress")
    if isinstance(progress, dict) and progress.get("engine") == "bindcraft2":
        return progress
    if row.get("node_id") != config.identity() or not _is_bindcraft2_job(row):
        return None
    run_dir = row.get("run_dir")
    if not isinstance(run_dir, str) or not run_dir:
        return None
    try:
        return job_progress.bindcraft2_progress(Path(run_dir))
    except OSError:
        return None


def _progress_summary(progress):
    if not isinstance(progress, dict):
        return ""
    completed = int(progress.get("trajectories_completed") or 0)
    active = int(progress.get("trajectories_in_progress") or 0)
    interrupted = int(progress.get("trajectories_interrupted") or 0)
    accepted = int(progress.get("designs_accepted") or 0)
    trajectory_limit = progress.get("trajectory_limit")
    design_target = progress.get("design_target")
    trajectory_label = f"{completed:,} / {int(trajectory_limit):,}" if trajectory_limit else f"{completed:,} completed"
    design_label = f"{accepted:,} / {int(design_target):,}" if design_target else f"{accepted:,} accepted"
    activity_label = f"{interrupted:,} interrupted" if interrupted else f"{active:,} running"
    return f"Progress: {trajectory_label} trajectories · {activity_label} · {design_label} accepted"


_RESULTS_TERMINAL_STATUSES = {"completed", "succeeded", "failed", "cancelled", "canceled"}


def _results_job_reference(row):
    return str(row.get("cluster_job_id") or (str(row.get("task_group") or "") + "/" + str(row.get("run_id") or "")))


def _is_visual_pipeline_row(row):
    return str(row.get("tool") or "").casefold() == "visual_pipeline" or str(row.get("job_type") or "").casefold() == "visual_pipeline"


def _visual_pipeline_groups(rows):
    """Group visual-pipeline children with their parent run, keyed by its stable run reference."""
    roots = {}
    for row in rows:
        if _is_visual_pipeline_row(row) and row.get("run_id"):
            roots[(str(row.get("node_id") or ""), str(row["run_id"]))] = row

    groups = {_results_job_reference(root): [root] for root in roots.values()}
    row_to_group = {reference: reference for reference in groups}
    for row in rows:
        if _is_visual_pipeline_row(row):
            continue
        parent_id = str(row.get("campaign_id") or "").strip()
        if not parent_id:
            continue
        parent_node = str(row.get("pipeline_parent_node_id") or row.get("node_id") or "")
        root = roots.get((parent_node, parent_id))
        if not root and not row.get("pipeline_parent_node_id") and row.get("pipeline_iteration"):
            # Legacy active-learning children used a parent path without its
            # workstation. Accept only a unique parent in the same access scope.
            candidates = [candidate for (_, run_id), candidate in roots.items() if run_id == parent_id
                          and str(candidate.get("project_id") or "") == str(row.get("project_id") or "")
                          and str(candidate.get("owner_id") or "") == str(row.get("owner_id") or "")]
            root = candidates[0] if len(candidates) == 1 else None
        if not root:
            continue
        row_project = str(row.get("project_id") or "")
        root_project = str(root.get("project_id") or "")
        if row_project and root_project and row_project != root_project:
            continue
        if row.get("owner_id") and root.get("owner_id") and row["owner_id"] != root["owner_id"]:
            continue
        root_ref = _results_job_reference(root)
        row_ref = _results_job_reference(row)
        if row_ref == root_ref:
            continue
        groups[root_ref].append(row)
        row_to_group[row_ref] = root_ref

    for members in groups.values():
        members[1:] = sorted(
            members[1:],
            key=lambda row: (str(row.get("created_at") or ""), str(row.get("job_code") or "")),
        )
    return groups, row_to_group


def _pipeline_iterations(children):
    grouped = {}
    for child in children:
        try:
            iteration = int(child.get("pipeline_iteration") or 0)
        except (TypeError, ValueError):
            iteration = 0
        grouped.setdefault(iteration, []).append(child)
    return dict(sorted(grouped.items())) if any(number > 0 for number in grouped) else {}


def _pipeline_iteration_status(root, number, children, highest):
    statuses = {str(child.get("status") or "unknown").casefold() for child in children}
    if statuses & {"failed"}:
        return "Failed"
    if statuses & {"running", "starting", "preparing"}:
        return "Running"
    if not statuses.issubset(_RESULTS_TERMINAL_STATUSES):
        return "Queued"
    if str(root.get("status") or "").casefold() not in _RESULTS_TERMINAL_STATUSES and number == highest:
        return "In progress"
    if statuses & {"cancelled", "canceled"}:
        return "Cancelled"
    return "Completed"


def _pipeline_step_count_label(root, children):
    iterations = _pipeline_iterations(children)
    if iterations:
        if str(root.get("status") or "").casefold() not in _RESULTS_TERMINAL_STATUSES:
            return f"Iteration {max(iterations)}"
        return f"{len([number for number in iterations if number > 0])} iterations"
    if not children:
        return "No jobs created yet" if str(root.get("status") or "").casefold() not in _RESULTS_TERMINAL_STATUSES else "No linked jobs"
    finished = sum(str(child.get("status") or "unknown").casefold() in _RESULTS_TERMINAL_STATUSES for child in children)
    qualifier = " finished" if str(root.get("status") or "").casefold() in _RESULTS_TERMINAL_STATUSES else " created jobs finished"
    return f"{finished}/{len(children)}{qualifier}"


def _pipeline_step_label(row, pipeline_name):
    name = str(row.get("campaign_step") or row.get("campaign_name") or "").strip()
    prefix = str(pipeline_name or "").strip() + " · "
    if prefix != " · " and name.startswith(prefix):
        name = name[len(prefix):].strip()
    return name or display_tool_label(row)


def _update_results_table_selection(table_key, entry_members):
    selected_jobs = st.session_state.setdefault("cluster-results-selected", {})
    edits = st.session_state.get(table_key, {}).get("edited_rows", {})
    for index, changes in edits.items():
        if "Select" not in changes or not 0 <= int(index) < len(entry_members):
            continue
        for member in entry_members[int(index)]:
            reference = _results_job_reference(member)
            if changes["Select"]:
                selected_jobs[reference] = member
            else:
                selected_jobs.pop(reference, None)
    st.session_state["cluster-results-selected"] = selected_jobs
    st.session_state.pop("cluster-results-pending-remove-and-purge", None)


def _request_results_pipeline_steps(click_key, entry_refs):
    click = st.session_state.get(click_key)
    if not click:
        return
    index = int(click["row"])
    if 0 <= index < len(entry_refs) and entry_refs[index]:
        st.session_state["cluster-results-pipeline-details-request"] = entry_refs[index]
        st.session_state.pop("cluster-results-pipeline-iteration-request", None)


def _request_pipeline_iteration_view(click_key, reference, numbers):
    click = st.session_state.get(click_key)
    if click and 0 <= int(click["row"]) < len(numbers):
        st.session_state["cluster-results-pipeline-details-request"] = reference
        st.session_state["cluster-results-pipeline-iteration-request"] = numbers[int(click["row"])]


def _request_results_pipeline_results(click_key, entry_refs):
    click = st.session_state.get(click_key)
    if click and 0 <= int(click["row"]) < len(entry_refs) and entry_refs[int(click["row"])]:
        _navigate_to_job_details(entry_refs[int(click["row"])], final_results=True)


def _request_results_job_view(click_key, entry_refs):
    click = st.session_state.get(click_key)
    if not click:
        return
    index = int(click["row"])
    if 0 <= index < len(entry_refs) and entry_refs[index]:
        _navigate_to_job_details(entry_refs[index])


def _request_pipeline_dialog_job_view(click_key, entry_refs):
    click = st.session_state.get(click_key)
    if not click:
        return
    index = int(click["row"])
    if 0 <= index < len(entry_refs) and entry_refs[index]:
        _navigate_to_job_details(entry_refs[index])


def _request_results_worker_details(click_key, entry_workers):
    click = st.session_state.get(click_key)
    if not click:
        return
    index = int(click["row"])
    if 0 <= index < len(entry_workers) and entry_workers[index]:
        st.session_state["cluster-results-library-workers-request"] = entry_workers[index]


def _library_worker_details(row, members):
    workers = {}
    for member in members:
        metrics = member.get("metrics") if isinstance(member.get("metrics"), dict) else {}
        plans = member.get("worker_plan") if isinstance(member.get("worker_plan"), list) else []
        actual = metrics.get("workers") if isinstance(metrics.get("workers"), list) else []
        for item in plans + actual:
            if not isinstance(item, dict):
                continue
            node_id = str(item.get("node_id") or item.get("name") or "")
            if not node_id:
                continue
            current = workers.setdefault(node_id, {"node_id": node_id})
            current.update(item)
            if item in actual:
                current["completed"] = str(item.get("status") or "").casefold() == "completed"
    if not workers:
        return None
    return {
        "job_name": config.display_text(row.get("campaign_name") or row.get("tool") or "Library preparation"),
        "job_status": str(row.get("status") or "unknown").replace("_", " ").title(),
        "current_phase": str(row.get("current_phase") or row.get("progress") or ""),
        "task_kind": "centroids" if any(w.get("task_kind") == "centroids" for w in workers.values()) else "fingerprints",
        "workers": sorted(workers.values(), key=lambda item: str(item.get("name") or "").casefold()),
    }


def _navigate_to_job_details(reference, *, final_results=False):
    st.session_state["pipeline-final-results-reference"] = reference if final_results else None
    st.session_state["selected-job-reference"] = reference
    st.session_state["page"] = "Job details"
    st.rerun(scope="app")


def _job_details_reference(row):
    if row.get("node_id") == config.identity() and row.get("local_reference"):
        return str(row["local_reference"])
    if row.get("cluster_job_id"):
        return "cluster:" + str(row["cluster_job_id"])
    return _results_job_reference(row)


@st.dialog("Pipeline steps", width="large")
def _render_results_pipeline_steps(members, elapsed_label=None, user=None, token=None, submit=None, iteration=None):
    root, *children = members
    pipeline_name = config.display_text(root.get("pipeline_name") or root.get("campaign_name") or "Pipeline run")
    st.markdown(f"**{pipeline_name}**")
    st.caption(
        f"{str(root.get('status') or 'unknown').replace('_', ' ').title()} · "
        f"{root.get('job_code') or root.get('run_id') or ''} · "
        f"{config.display_name(root.get('node_name') or 'Not assigned')}"
    )
    st.button(
        "View pipeline run",
        key="cluster-results-pipeline-root-view",
        type="tertiary",
        on_click=_navigate_to_job_details,
        args=(_job_details_reference(root),),
    )
    if st.button("Final results" if str(root.get("status") or "").casefold() in {"completed", "succeeded"} else "Results so far",
                 key="cluster-pipeline-final-results", type="primary"):
        _navigate_to_job_details(_job_details_reference(root), final_results=True)
    if not children:
        st.info("No linked tool-call jobs are available for this pipeline.")
        return
    metrics = root.get("metrics") if isinstance(root.get("metrics"), dict) else {}
    if metrics.get("round_count"):
        docking = sum(str(child.get("tool") or "") == "unidock_pro" for child in children)
        models = sum(str(child.get("tool") or "") == "active_learning_model" for child in children)
        st.caption(f"{int(metrics['round_count'])} rounds · {docking} docking batches · {models} neural-network jobs · "
                   f"{len({child.get('node_id') for child in children if child.get('node_id')})} computers")
    if str(root.get("status") or "").casefold() not in _RESULTS_TERMINAL_STATUSES:
        st.caption("Counts include jobs created so far. A looping pipeline adds jobs until its stopping criteria or iteration limit are reached.")
    iterations = _pipeline_iterations(children)
    if iterations:
        numbers = list(iterations)
        iteration_rows = []
        for number, jobs in iterations.items():
            iteration_rows.append({"Iteration": f"Iteration {number}" if number else "Other jobs",
                "Status": _pipeline_iteration_status(root, number, jobs, max(numbers)),
                "Batches created": sum(str(child.get("tool")) == "unidock_pro" for child in jobs),
                "Neural-network jobs": sum(str(child.get("tool")) == "active_learning_model" for child in jobs),
                "Computers": len({child.get("node_id") for child in jobs if child.get("node_id")}),
                "Inspect": "View jobs"})
        key = "pipeline-iteration-inspect"
        st.dataframe(pd.DataFrame(iteration_rows).style.map(lambda _: "color: #2f6fba; font-weight: 600;", subset=["Inspect"]),
            hide_index=True, width="stretch", column_config={"Inspect": st.column_config.ButtonColumn(
                "Inspect", type="tertiary", key=key, on_click=_request_pipeline_iteration_view,
                args=(key, _results_job_reference(root), numbers))})
        st.caption("Each row represents one active-learning iteration. Open an iteration to inspect its individual batch and neural-network jobs.")
        if iteration not in iterations:
            return
        st.subheader(f"Iteration {iteration} · job details" if iteration else "Other job details")
        children = iterations[iteration]
    child_refs = []
    step_rows = []
    for child in children:
        child_ref = _results_job_reference(child)
        child_refs.append(child_ref)
        step_rows.append({
        "Step": _pipeline_step_label(child, pipeline_name),
        "View": "View",
        "Tool": display_tool_label(child),
        "#": str(child.get("job_code") or ""),
        "Status": str(child.get("status") or "unknown").replace("_", " ").title(),
        "Elapsed": elapsed_label(child) if elapsed_label else "—",
        "Workstation": config.display_name(child.get("node_name") or "Not assigned"),
        "Progress": _progress_summary(child.get("progress")),
        })
    view_key = "cluster-results-pipeline-child-view"
    styled_steps = pd.DataFrame(step_rows).style.map(
        lambda _: "color: #2f6fba; font-weight: 600;",
        subset=["View"],
    )
    st.dataframe(
        styled_steps,
        hide_index=True,
        width="stretch",
        column_config={"View": st.column_config.ButtonColumn(
            "View", width="small", type="tertiary", key=view_key,
            on_click=_request_pipeline_dialog_job_view,
            args=(view_key, [_job_details_reference(child) for child in children]),
        )},
        height=min(420, 36 * len(step_rows) + 42),
    )


@st.dialog("Distributed library workers", width="large")
def _render_library_worker_details(details):
    workers = details.get("workers") or []
    centroids = details.get("task_kind") == "centroids"
    st.markdown(f"**{details.get('job_name') or 'Library preparation'}**")
    caption = f"{details.get('job_status') or 'Unknown'} · One Studio job · {len(workers)} computers"
    if details.get("current_phase"):
        caption += f" · {details['current_phase']}"
    st.caption(caption)
    st.caption("Random fingerprint chunks are sent from the source workstation to CPU workers. Representatives return after all workers in each round finish; the original library stays in place."
               if centroids else "The original library file stays on the organizer. Record batches stream to the selected computers; generated outputs return after all workers finish.")

    selected_key = "cluster-results-library-worker-selected"
    selected_node = st.session_state.get(selected_key)
    by_node = {str(worker.get("node_id") or worker.get("name")): worker for worker in workers}
    if selected_node not in by_node:
        selected_node = None
        st.session_state.pop(selected_key, None)

    if selected_node:
        worker = by_node[selected_node]
        back_col, _ = st.columns([0.2, 0.8])
        if back_col.button("All computers", key="cluster-library-worker-back"):
            selected_node = None
            st.session_state.pop(selected_key, None)
        if selected_node:
            st.subheader(config.display_name(worker.get("name") or selected_node))
            running = details.get("job_status", "").casefold() in {"running", "starting"}
            status = "completed" if worker.get("completed") else "in progress" if running else "planned"
            job_completed = details.get("job_status", "").casefold() in {"completed", "succeeded"}
            summary = [
                {"Step": "Receive source records", "Status": "Complete" if worker.get("completed") else status.title(),
                 "Details": (f"{int(worker.get('input_rows') or 0):,} rows streamed from the organizer"
                             if worker.get("input_rows") is not None else "Source rows are assigned across the selected computers.")},
                {"Step": "Calculate Morgan fingerprints", "Status": "Complete" if worker.get("completed") else status.title(),
                 "Details": (f"Feature Morgan · radius 2 · 1024 bits · {int(worker.get('cpu_worker_count') or 0)} CPU processes · "
                             f"{int(worker.get('row_count') or 0):,} valid compounds"
                             if worker.get("completed") else "CPU processes calculate fingerprints in bounded batches.")},
                {"Step": "Write compact Parquet output", "Status": "Complete" if worker.get("completed") else status.title(),
                 "Details": (f"{int(worker.get('shard_count') or 0):,} fingerprint shards · "
                             f"{int(worker.get('rejected_rows') or 0):,} rejected rows · "
                             f"{int(worker.get('output_size_bytes') or 0):,} bytes"
                             if worker.get("completed") else "Each worker writes compressed output and a rejected-row table.")},
                {"Step": "Return output to organizer", "Status": "Complete" if worker.get("completed") else "After all workers finish",
                 "Details": "Only generated fingerprint outputs are copied back; the original source file stays on the organizer."},
                {"Step": "Organizer merge and validation", "Status": "Complete" if job_completed else "After worker transfer",
                 "Details": "The organizer verifies shard checksums and unique IDs, then writes the final manifest."},
            ]
            if centroids:
                phase_status = "Complete" if worker.get("completed") else status.title()
                summary = [
                    {"Step": "Receive random balanced chunks", "Status": phase_status,
                     "Details": f"{int(worker.get('chunk_count') or 0):,} chunks across all rounds · {int(worker.get('input_rows') or 0):,} input rows"},
                    {"Step": "Cluster on CPUs", "Status": phase_status,
                     "Details": f"{int(worker.get('cpu_worker_count') or 0)} scheduled CPU slots · {int(worker.get('parallel_chunks') or 0)} concurrent chunks · {int(worker.get('cores_per_chunk') or 0)} threads per chunk"},
                    {"Step": "Write representatives", "Status": phase_status,
                     "Details": f"{int(worker.get('row_count') or 0):,} representatives produced across all rounds · compressed Parquet"},
                    {"Step": "Return outputs to coordinator", "Status": phase_status if worker.get("completed") else "After all workers finish the round",
                     "Details": "Generated representatives are collected and their checksums and row counts verified."},
                    {"Step": "Final comparison and retention", "Status": "Complete" if job_completed else "After chunk collection",
                     "Details": "Run the sampled final comparison or further merge rounds. Keep the final centroid table and manifest; remove temporary chunks."},
                ]
            st.dataframe(pd.DataFrame(summary), hide_index=True, use_container_width=True)
    else:
        for index, worker in enumerate(workers):
            node_id = str(worker.get("node_id") or worker.get("name") or index)
            name = config.display_name(worker.get("name") or node_id)
            cpu_count = int(worker.get("cpu_worker_count") or 0)
            completed = bool(worker.get("completed"))
            status = "Completed" if completed else "Assigned" if details.get("job_status", "").casefold() in {"running", "starting"} else "Planned"
            left, middle, right = st.columns([0.48, 0.25, 0.27], vertical_alignment="center")
            left.markdown(f"**{name}**")
            middle.caption(f"{cpu_count} CPU slots")
            if completed:
                right.caption(f"{int(worker.get('row_count') or 0):,} representatives · {int(worker.get('chunk_count') or 0):,} chunks"
                              if centroids else f"{int(worker.get('row_count') or 0):,} fingerprints · {int(worker.get('rejected_rows') or 0):,} rejected")
            else:
                right.caption(status)
            if st.button(f"View {name} steps", key=f"cluster-library-worker-view-{node_id}", type="tertiary"):
                selected_node = node_id
                st.session_state[selected_key] = node_id


def _render_bindcraft2_progress(row, job_id, token):
    terminal = {"completed", "succeeded", "failed", "cancelled"}
    run_every = 15 if row.get("status") not in terminal else None

    @st.fragment(run_every=run_every)
    def live_progress():
        current = row
        try:
            current = _request("/jobs/" + job_id, token)
        except (transport.ClusterUnavailable, ValueError):
            pass
        progress = _row_progress(current)
        st.subheader("BindCraft progress")
        if not progress:
            st.caption("Live campaign counters will appear when BindCraft starts writing progress.")
            return

        completed = int(progress.get("trajectories_completed") or 0)
        active = int(progress.get("trajectories_in_progress") or 0)
        accepted = int(progress.get("designs_accepted") or 0)
        trajectory_limit = progress.get("trajectory_limit")
        design_target = progress.get("design_target")
        attempts_col, accepted_col = st.columns(2)
        attempts_value = f"{completed:,}" + (f" / {int(trajectory_limit):,}" if trajectory_limit else "")
        accepted_value = f"{accepted:,}" + (f" / {int(design_target):,}" if design_target else "")
        attempts_col.metric("Completed trajectories", attempts_value, delta=f"{active:,} running" if active else None)
        accepted_col.metric("Accepted designs", accepted_value)
        if trajectory_limit:
            attempts_col.progress(min(1.0, completed / max(1, int(trajectory_limit))), text="Attempt cap")
        if design_target:
            accepted_col.progress(min(1.0, accepted / max(1, int(design_target))), text="Accepted-design goal")
        if progress.get("updated_at"):
            st.caption("Counters update while the campaign runs; the trajectory count is completed attempts, with active workers shown separately.")

    live_progress()


def _candidate_structure_path(candidate, data):
    """Return a run-relative candidate structure path safe for the artifact API."""
    if not isinstance(candidate, dict) or not isinstance(data, dict):
        return None
    run_dir = (data.get("job") or {}).get("run_dir")
    if not isinstance(run_dir, str) or not run_dir:
        return None
    root = PurePosixPath(run_dir)
    for key in ("complex_pdb", "binder_pdb"):
        value = candidate.get(key)
        if not isinstance(value, str) or not value.strip():
            continue
        path = PurePosixPath(value.strip())
        if path.is_absolute():
            try:
                path = path.relative_to(root)
            except ValueError:
                continue
        if not path.parts or any(part in {".", ".."} for part in path.parts):
            continue
        return str(path), key
    return None


_BINDCRAFT2_RESULT_COLUMNS = (
    "i_pDAE",
    "i_pTM",
    "pLDDT",
    "pTM",
    "i_pAE",
    "Unbound_Binder_pLDDT",
    "Target_pLDDT",
    "Interface_Residues",
    "Binder_RMSD",
    "Target_Crop_Length",
    "Backbone_Clashes",
    "All_Atom_Clashes",
    "Interface_BuriedArea",
    "Interface_BuriedArea_Fraction",
    "Surface_Hydrophobicity",
    "Interface_Hydrophobicity",
    "Binder_Helix_Fraction",
    "Binder_BetaSheet_Fraction",
    "Binder_Loop_Fraction",
    "SS_pLDDT",
    "Interface_Helix_Fraction",
    "Interface_BetaSheet_Fraction",
    "Interface_Loop_Fraction",
    "Binder_Disulfides",
    "Binder_Mass_kDa",
    "Binder_pI",
    "Binder_Net_Charge",
    "Binder_Cysteines",
    "Binder_Extinction",
    "Binder_Free_Cysteines",
    "Interface_Binder_Residues",
    "Interface_Target_Residues",
)


def _candidate_result_rows(candidates, tool=""):
    """Build a cross-tool result table, retaining BindCraft2's native metrics."""
    is_bindcraft2 = "bindcraft2" in str(tool).casefold() or any(
        isinstance(candidate, dict)
        and "bindcraft2" in str(candidate.get("source_tool") or candidate.get("tool") or "").casefold()
        for candidate in candidates
    )
    metric_columns = (
        ("pred_lddt", "pLDDT"),
        ("interaction_pae", "Interaction PAE"),
        ("pae", "PAE"),
        ("hotspot_atom_contact_fraction", "Hotspot atom contacts"),
        ("hotspot_ca_contact_fraction", "Hotspot CA contacts"),
        ("rf2_repeat", "RF2 repeat"),
        ("cyclic", "Cyclic"),
    )
    rows = []
    for rank, candidate in enumerate(candidates, 1):
        if not isinstance(candidate, dict):
            continue
        metrics = candidate.get("metrics") if isinstance(candidate.get("metrics"), dict) else {}
        raw_metadata = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
        native_row = raw_metadata.get("bindcraft2_row") if isinstance(raw_metadata.get("bindcraft2_row"), dict) else {}
        row = {
            "Rank": metrics.get("bindcraft2_rank", native_row.get("rank", rank)) if is_bindcraft2 else rank,
            "Design": (
                native_row.get("design")
                or raw_metadata.get("bindcraft2_design")
                or candidate.get("candidate_id")
                or f"Design {rank}"
            ),
            "Sequence": candidate.get("binder_sequence") or candidate.get("sequence") or "—",
            "Length": candidate.get("binder_length") or metrics.get("binder_length"),
        }
        if is_bindcraft2:
            for key in _BINDCRAFT2_RESULT_COLUMNS:
                value = metrics.get(key, native_row.get(key))
                if value not in (None, ""):
                    row[key] = value
            if "native_pass_filters" in metrics:
                row["Passed filters"] = metrics["native_pass_filters"]
        else:
            for key, label in metric_columns:
                if key in metrics:
                    row[label] = metrics[key]
        rows.append(row)
    return rows


def _bindcraft2_run_parameter_rows(data):
    """Select human-readable BindCraft2 inputs and filter thresholds from inspect data."""
    input_payload = data.get("input") if isinstance(data, dict) else None
    if not isinstance(input_payload, dict):
        return [], []
    params = input_payload.get("params") if isinstance(input_payload.get("params"), dict) else {}
    inputs = input_payload.get("inputs") if isinstance(input_payload.get("inputs"), dict) else {}

    def present(value):
        return value not in (None, "", [], {})

    rows = []
    values = (
        ("Campaign", params.get("campaign_name") or input_payload.get("campaign_name")),
        ("Target", params.get("target_name") or inputs.get("target_name")),
        ("Target type", params.get("target_type") or inputs.get("target_type")),
        ("Target chains", params.get("target_chains") or inputs.get("target_chains")),
        ("Design format", params.get("modality")),
        ("Binder length range", params.get("binder_length")),
        ("Designs requested", params.get("number_of_final_designs")),
        ("Maximum trajectories", params.get("max_trajectories")),
        ("Target crop count", params.get("target_crop_count")),
        ("Target crop lengths", params.get("target_crop_lengths")),
        ("Hotspots", params.get("hotspots")),
        ("Coldspots", params.get("coldspots")),
        ("Random seed", params.get("campaign_seed")),
        ("GPU device", params.get("gpu_device")),
    )
    for label, value in values:
        if present(value):
            if isinstance(value, (list, tuple)):
                value = ", ".join(map(str, value))
            rows.append({"Parameter": label, "Value": value})

    quality = params.get("quality_filters") if isinstance(params.get("quality_filters"), dict) else {}
    filters = quality.get("filters") if isinstance(quality.get("filters"), dict) else {}
    filter_rows = []
    for name, value in quality.items():
        if name == "filters" or value is None:
            continue
        try:
            value = f"{float(value):g}"
        except (TypeError, ValueError):
            value = str(value)
        key = str(name)
        direction = "≤" if key.startswith("max_") else "≥" if key.startswith("min_") else ""
        detail = key.removeprefix("min_").removeprefix("max_")
        parts = detail.split("_")
        if parts[0] in {"iptm", "plddt"} and len(parts) > 1:
            metric = "i_pTM" if parts[0] == "iptm" else "pLDDT"
            label = f"{parts[-1].title()} {metric} stage minimum"
        else:
            label = detail.replace("_", " ").title()
        filter_rows.append({"Filter": label, "Requirement": f"{direction} {value}".strip()})
    for name, config in filters.items():
        if not isinstance(config, dict) or config.get("threshold") is None:
            continue
        threshold = config["threshold"]
        try:
            threshold = float(threshold)
            threshold = f"{threshold:g}"
        except (TypeError, ValueError):
            threshold = str(threshold)
        direction = "≥" if config.get("higher") else "≤"
        filter_rows.append({"Filter": str(name).replace("_", " "), "Requirement": f"{direction} {threshold}"})
    filter_rows.sort(key=lambda item: item["Filter"].casefold())
    return rows, filter_rows


def _safe_run_artifact_path(value):
    if not isinstance(value, str) or not value.strip():
        return None
    path = PurePosixPath(value.strip())
    if path.is_absolute() or not path.parts or any(part in {".", ".."} for part in path.parts):
        return None
    return str(path)


def _bindcraft2_refolded_table_path(data):
    result = data.get("result") if isinstance(data, dict) else None
    outputs = result.get("outputs") if isinstance(result, dict) else None
    artifacts = outputs.get("artifacts") if isinstance(outputs, dict) else None
    for artifact in artifacts if isinstance(artifacts, list) else []:
        if isinstance(artifact, dict) and artifact.get("type") == "bindcraft2_refolded_table":
            path = _safe_run_artifact_path(artifact.get("path"))
            if path:
                return path
    return "artifacts/raw/bindcraft2/output/2_Refolded/!_Refolded.csv"


def _bindcraft2_monomer_path(refolded_csv, candidate):
    raw = candidate.get("raw_metadata") if isinstance(candidate, dict) else None
    raw = raw if isinstance(raw, dict) else {}
    direct_path = _safe_run_artifact_path(raw.get("bindcraft2_binder_monomer"))
    if direct_path:
        return direct_path
    if not refolded_csv:
        return None
    native = raw.get("bindcraft2_row") if isinstance(raw.get("bindcraft2_row"), dict) else {}
    wanted_sequence = "".join(str(candidate.get("binder_sequence") or native.get("Binder_Sequence") or "").split()).upper()
    matches = []
    if wanted_sequence:
        csv_text = str(refolded_csv or "").lstrip("\ufeff")
        for row in csv.DictReader(io.StringIO(csv_text)):
            sequence = "".join(str(row.get("Binder_Sequence") or "").split()).upper()
            if sequence != wanted_sequence:
                continue
            failed = str(row.get("failed_filters") or "").strip()
            outcome = str(row.get("outcome") or "").strip().casefold()
            matches.append((not failed and outcome in {"", "passed"}, row))
    if not matches:
        design = str(raw.get("bindcraft2_design") or native.get("design") or "").strip()
    else:
        _, matched = next((item for item in matches if item[0]), matches[0])
        design = str(matched.get("design") or "").strip()
    design = PurePosixPath(design).name
    if not design or design in {".", ".."}:
        return None
    if design.endswith("_monomer.cif"):
        return _safe_run_artifact_path(
            f"artifacts/raw/bindcraft2/output/2_Refolded/BinderMonomer/{design}"
        )
    if design.endswith((".cif", ".mmcif", ".pdb")):
        design = PurePosixPath(design).stem
    return _safe_run_artifact_path(
        f"artifacts/raw/bindcraft2/output/2_Refolded/BinderMonomer/{design}_monomer.cif"
    )


def _fetch_remote_artifact_text(jid, path, token, state_key, label="saved structure"):
    relative = _safe_run_artifact_path(path)
    if not relative:
        st.warning("The saved structure path is invalid.")
        return None
    command_id = st.session_state.get(state_key)
    if not command_id:
        try:
            operation = _request(
                "/jobs/" + jid + "/actions",
                token,
                {"action": "artifact", "parameters": {"path": relative}},
                method="POST",
            )
            command_id = operation["command_id"]
            st.session_state[state_key] = command_id
        except (ValueError, transport.ClusterUnavailable) as exc:
            st.warning(f"Could not request the saved structure: {exc}")
            return None
    try:
        status = _request("/commands/" + command_id, token)
        if status["status"] == "failed":
            error = status.get("error") or f"The workstation could not serve this {label}."
            st.warning(f"{error} Requested path: {relative}")
            return None
        if status["status"] != "completed":
            st.info(f"Loading the {label} from the owning workstation…")
            return None
        content_key = "cluster-artifact-bytes-" + command_id
        if content_key not in st.session_state:
            st.session_state[content_key] = _request(
                "/commands/" + command_id + "/file", token, raw=True, timeout=60
            )
        return st.session_state[content_key].decode("utf-8", errors="replace")
    except (ValueError, transport.ClusterUnavailable) as exc:
        st.info(str(exc))
        return None


def _render_remote_docking_preview(jid, row, data, token):
    tool = " ".join(
        str(value or "")
        for value in (row.get("tool"), row.get("workflow"), row.get("job_type"), display_tool_label(row))
    ).casefold().replace("-", "_")
    if not any(name in tool for name in ("vina", "autodock_gpu", "gnina", "unidock_pro", "rosetta_ligand", "rosettaligand")):
        return

    st.subheader("Docked protein–ligand complex")
    score_text = _fetch_remote_artifact_text(
        jid,
        "artifacts/docking_scores.csv",
        token,
        f"cluster-docking-scores-{jid}",
        label="docking score table",
    )
    pose_rows = []
    if score_text:
        for item in csv.DictReader(io.StringIO(score_text)):
            pose_path = _safe_run_artifact_path(str(item.get("pose_path") or ""))
            if not pose_path:
                continue
            item["pose_path"] = pose_path
            compound = str(item.get("compound_id") or Path(pose_path).stem)
            details = []
            if item.get("replicate"):
                details.append(f"replicate {item['replicate']}")
            if item.get("pose_rank"):
                details.append(f"pose {item['pose_rank']}")
            for key, unit in (("estimated_dg_reu", "REU"), ("best_score_kcal_mol", "kcal/mol")):
                try:
                    details.append(f"{float(item.get(key)):g} {unit}")
                    break
                except (TypeError, ValueError):
                    continue
            item["_label"] = compound + (" · " + " · ".join(details) if details else "")
            pose_rows.append(item)
    if not pose_rows:
        st.info("No docked pose is recorded in the saved score table.")
        return

    result = data.get("result") if isinstance(data, dict) else {}
    metrics = result.get("metrics") if isinstance(result, dict) and isinstance(result.get("metrics"), dict) else {}
    best_compound = str(metrics.get("best_compound_id") or "").casefold()
    default_index = next(
        (index for index, item in enumerate(pose_rows)
         if str(item.get("compound_id") or "").casefold() == best_compound),
        0,
    )
    selected = st.selectbox(
        "Docked pose",
        range(len(pose_rows)),
        index=default_index,
        format_func=lambda index: pose_rows[index]["_label"],
        key=f"cluster-docking-pose-{jid}",
    )
    pose_path = pose_rows[selected]["pose_path"]
    pose_text = _fetch_remote_artifact_text(
        jid, pose_path, token, f"cluster-docking-pose-text-{jid}-{selected}", label="docked pose"
    )
    if pose_text is None:
        return

    artifacts = [
        _safe_run_artifact_path(str(item.get("relative_path") or ""))
        for item in row.get("artifacts", [])
        if isinstance(item, dict)
    ]
    artifacts = [path for path in artifacts if path]
    receptor_relative = None
    for suffix in (
        "artifacts/prepared/receptor.pdbqt",
        "prepared/receptor.pdb",
        "artifacts/input/receptor.pdb",
        "input/receptor.pdb",
    ):
        receptor_relative = next(
            (suffix for artifact in artifacts if artifact.casefold().endswith("/" + suffix.casefold()) or artifact.casefold() == suffix.casefold()),
            None,
        )
        if receptor_relative:
            break
    receptor_text = None
    if receptor_relative:
        receptor_text = _fetch_remote_artifact_text(
            jid,
            receptor_relative,
            token,
            f"cluster-docking-receptor-text-{jid}",
            label="prepared receptor",
        )
        if receptor_text is None:
            return

    from mn_molecular_design_studio import docking_results
    pose_name = Path(pose_path).name
    structures = docking_results.pose_visualizations(
        pose_text,
        pose_name,
        receptor_text,
        Path(receptor_relative).name if receptor_relative else "",
    )
    if pose_name.casefold().endswith(".pdbqt") and receptor_text is None:
        st.warning("The prepared receptor is unavailable, so this preview can show only the ligand pose.")
    from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer import molstar_custom_component
    molstar_custom_component(
        structures,
        key=f"cluster-docking-complex-{jid}-{selected}",
        height=520,
        show_controls=False,
        show_structure_tools=False,
    )
    st.caption(
        f"Pose: `{pose_name}` · "
        + (f"receptor: `{Path(receptor_relative).name}`" if receptor_relative else "protein and ligand from the docked pose file")
    )


def _align_bindcraft2_monomer_to_complex(complex_text, monomer_text, candidate):
    """Fit the binder-alone CIF onto its binder chain in the predicted complex."""
    from Bio.PDB import MMCIFIO, MMCIFParser, Superimposer
    from mn_molecular_design_studio.design_validation import protein_chains, residue_pairs, sequence
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import _mmcif_text_with_default_occupancy

    def parse(text, name):
        prepared = _mmcif_text_with_default_occupancy(str(text))
        return MMCIFParser(QUIET=True).get_structure(name, io.StringIO(prepared))

    complex_structure = parse(complex_text, "bindcraft2_complex")
    monomer_structure = parse(monomer_text, "bindcraft2_binder_monomer")
    complex_chains = protein_chains(complex_structure)
    monomer_chains = protein_chains(monomer_structure)
    expected = str(candidate.get("binder_sequence") or "").replace(" ", "").upper()
    if not expected:
        raise ValueError("The candidate sequence is missing, so its binder chain cannot be matched.")

    def best_chain(chains, preferred=()):
        ordered = [chain for chain in preferred if chain in chains]
        ordered.extend(chain for chain in chains if chain not in ordered)
        scored = []
        for chain in ordered:
            observed = sequence(chains[chain])
            pairs = residue_pairs(expected, observed)
            identity = sum(expected[i] == observed[j] for i, j in pairs) / max(1, len(pairs))
            scored.append((identity, chain, pairs))
        return max(scored, default=(0.0, None, []))

    binder_chains = candidate.get("binder_chains") or []
    complex_identity, complex_chain, _ = best_chain(complex_chains, binder_chains)
    monomer_identity, monomer_chain, _ = best_chain(monomer_chains)
    if complex_chain is None or monomer_chain is None or min(complex_identity, monomer_identity) < 0.5:
        raise ValueError("Could not identify the binder chain in both structures.")
    complex_residues, monomer_residues = complex_chains[complex_chain], monomer_chains[monomer_chain]
    complex_sequence, monomer_sequence = sequence(complex_residues), sequence(monomer_residues)
    pairs = residue_pairs(complex_sequence, monomer_sequence)
    pairs = [(i, j) for i, j in pairs if "CA" in complex_residues[i] and "CA" in monomer_residues[j]]
    if len(pairs) < 3:
        raise ValueError("Fewer than three matching binder C-alpha atoms are available for alignment.")
    fixed = [complex_residues[i]["CA"] for i, _ in pairs]
    moving = [monomer_residues[j]["CA"] for _, j in pairs]
    superimposer = Superimposer()
    superimposer.set_atoms(fixed, moving)
    for atom in monomer_structure.get_atoms():
        atom.transform(*superimposer.rotran)
    output = io.StringIO()
    writer = MMCIFIO()
    writer.set_structure(monomer_structure)
    writer.save(output)
    return output.getvalue(), float(superimposer.rms), len(pairs)


def _candidate_structure_color(candidate, structure_text=""):
    """Use pLDDT colors only when per-atom confidence values are in the PDB."""
    metrics = candidate.get("metrics") if isinstance(candidate, dict) else None
    metrics = metrics if isinstance(metrics, dict) else {}
    has_confidence_metric = any(
        isinstance(metrics.get(key), (int, float))
        for key in ("pred_lddt", "plddt", "mean_plddt")
    )
    if has_confidence_metric:
        factors = []
        for line in str(structure_text).splitlines():
            if line.startswith(("ATOM  ", "HETATM")):
                try:
                    factors.append(float(line[60:66]))
                except ValueError:
                    continue
        if len(factors) > 1 and max(factors) - min(factors) > 0.01:
            return "plddt"
    return "chain-id"


def _cache(user_id, name):
    safe = hashlib.sha256(str(user_id).encode()).hexdigest()[:24]
    return config.home() / (name + "-" + safe + ".json")


def catalog(user, token, *, include_hidden=False):
    path = _cache(user["id"], "catalog-cache-hidden" if include_hidden else "catalog-cache")
    try:
        query = "?limit=1000" + ("&include_hidden=true" if include_hidden else "")
        result = _request("/jobs" + query, token)
        total = result["total"]
        while len(result["items"]) < total:
            page = _request(f"/jobs?limit=1000&offset={len(result['items'])}" + ("&include_hidden=true" if include_hidden else ""), token)
            if not page["items"]:
                break
            result["items"].extend(page["items"])
        config.atomic_json(path, result)
        return result, False
    except transport.ClusterUnavailable:
        try:
            return json.loads(path.read_text()), True
        except (FileNotFoundError, json.JSONDecodeError):
            return {"items": [], "total": 0}, True


def _refresh_results_snapshot(user, token, local_search):
    """Load remote and local rows once for Results filtering and display."""
    page, stale = catalog(user, token, include_hidden=True)
    rows = list(page["items"])
    local_page = local_search(limit=1000)
    local = list(local_page["items"])
    while len(local) < local_page["total"]:
        more = local_search(limit=1000, offset=len(local))["items"]
        if not more:
            break
        local.extend(more)
    local_name = config.load()["name"]
    known = {(r.get("node_id"), r.get("task_group"), r.get("run_id")): i for i, r in enumerate(rows)}
    for row in local:
        key = (config.identity(), row["task_group"], row["run_id"])
        item = dict(row, node_id=config.identity(), node_name=local_name, availability="online")
        if _is_bindcraft2_job(item):
            item["progress"] = _row_progress(item)
        if key in known:
            rows[known[key]] = rows[known[key]] | item
        else:
            rows.append(item)
    for row in rows:
        if _is_bindcraft2_job(row) and not row.get("progress"):
            row["progress"] = _row_progress(row)
    return {
        "items": rows,
        "stale": stale,
        "user_id": str(user.get("id") or ""),
        "refreshed_at": time.monotonic(),
    }


def _gpu_counts(node):
    snapshot = node.get("snapshot") if isinstance(node.get("snapshot"), dict) else {}
    resources = snapshot.get("resources") if isinstance(snapshot.get("resources"), dict) else {}
    detected = resources.get("gpus")
    free = resources.get("free_gpus")
    if isinstance(detected, list):
        if not isinstance(free, list):
            return len(detected), None, None
        return len(detected), len(free), max(0, len(detected) - len(free))
    if isinstance(free, list):
        return None, len(free), None
    return None, None, None


def _gpu_inventory_rows(nodes, *, snapshot_current=True):
    rows = []
    for node in nodes:
        snapshot = node.get("snapshot") if isinstance(node.get("snapshot"), dict) else {}
        resources = snapshot.get("resources") if isinstance(snapshot.get("resources"), dict) else {}
        detected = resources.get("gpus")
        free = resources.get("free_gpus")
        free_known = isinstance(free, list)
        free = free if free_known else []
        devices = detected if isinstance(detected, list) else free
        free_ids = {
            str(gpu.get("gpu_id", gpu.get("index")))
            for gpu in free if isinstance(gpu, dict) and gpu.get("gpu_id", gpu.get("index")) is not None
        }
        is_online = snapshot_current and str(node.get("availability") or "").casefold() == "online"
        for gpu in devices:
            if not isinstance(gpu, dict):
                continue
            device_id = str(gpu.get("gpu_id") if gpu.get("gpu_id") is not None else gpu.get("index", "—"))
            available = device_id in free_ids
            if not free_known:
                state = "Allocation status not reported" if is_online else "Last reported allocation status unavailable"
            elif is_online:
                state = "Available now" if available else "Allocated now"
            else:
                state = "Last reported available" if available else "Last reported allocated"
            rows.append({
                "Workstation": config.display_name(str(node.get("name") or node.get("id") or "Unknown")),
                "Device": device_id,
                "Model": str(gpu.get("name") or gpu.get("model") or "GPU model not reported"),
                "Scheduler state": state,
                "Connection": str(node.get("availability") or "unknown").title(),
                "Resource snapshot": str(node.get("last_seen_at") or "Not timestamped"),
            })
    return rows


def execution_selector(user, token):
    if not config.enabled():
        return
    settings = config.load()
    snapshot_current = True
    try:
        nodes = _request("/nodes", token)["nodes"]
        config.atomic_json(_cache(user["id"], "nodes-cache"), {"nodes": nodes})
    except (transport.ClusterUnavailable, ValueError):
        snapshot_current = False
        try:
            nodes = json.loads(_cache(user["id"], "nodes-cache").read_text())["nodes"]
        except (FileNotFoundError, json.JSONDecodeError):
            nodes = []
    names = {n["id"]: f"{config.display_name(n['name'])} · {n['availability']}" for n in nodes if n["id"] != config.identity()}
    options = ["auto", "local", *names]
    if st.session_state.get("cluster-placement") not in options:
        st.session_state["cluster-placement"] = settings["default_placement"]
    st.sidebar.selectbox("Run new jobs on", options, format_func=lambda value: {"auto": "Automatic · available workstation", "local": "This workstation"}.get(value, names.get(value, value)), key="cluster-placement")
    st.sidebar.caption(f"This workstation: {config.display_name(settings['name'])}")
    if nodes:
        with st.sidebar.expander("Workstation resources"):
            pool_counts = [_gpu_counts(node) for node in nodes]
            pool_detected = sum(item[0] or 0 for item in pool_counts) if all(item[0] is not None for item in pool_counts) else None
            pool_available = sum(item[1] for item in pool_counts) if all(item[1] is not None for item in pool_counts) else None
            pool_allocated = sum(item[2] or 0 for item in pool_counts) if all(item[2] is not None for item in pool_counts) else None
            pool_age = "current coordinator response" if snapshot_current else "cached workstation snapshots"
            detected_text = str(pool_detected) if pool_detected is not None else "total not reported"
            available_text = str(pool_available) if pool_available is not None else "available count not reported"
            allocated_text = str(pool_allocated) if pool_allocated is not None else "allocated count not reported"
            st.caption(f"Pool snapshot ({pool_age}): {detected_text} detected · {available_text} · {allocated_text} GPUs.")
            for node in nodes:
                snapshot = node.get("snapshot") if isinstance(node.get("snapshot"), dict) else {}
                resource = snapshot.get("resources") if isinstance(snapshot.get("resources"), dict) else {}
                st.write(f"**{config.display_name(node.get('name') or node.get('id') or 'Unknown')}** · {node.get('availability', 'unknown')}")
                detected = resource.get("gpus")
                free = resource.get("free_gpus")
                if isinstance(free, list):
                    gpu_summary = f"{len(free)} free GPUs available"
                else:
                    free = []
                    gpu_summary = "GPU availability/allocations not reported"
                if isinstance(detected, list) and isinstance(resource.get("free_gpus"), list):
                    gpu_summary += f" · {len(detected)} detected · {max(0, len(detected) - len(free))} allocated"
                elif isinstance(detected, list):
                    gpu_summary += f" · {len(detected)} detected"
                else:
                    gpu_summary += " · total inventory not reported"
                if not snapshot_current:
                    gpu_summary += " (cached report)"
                st.caption(f"{gpu_summary} · {resource.get('cpu_free', 0)} free CPU threads · {resource.get('ram_available_gb', 0):.0f} GB RAM available")
                for gpu_row in _gpu_inventory_rows([node], snapshot_current=snapshot_current):
                    st.caption(f"GPU {gpu_row['Device']} · {gpu_row['Model']} · {gpu_row['Scheduler state']}")


def settings(user, token):
    st.subheader("Workstations")
    current = config.load()
    st.caption("Share the job catalog and send new work to available workstations. Full results stay on their execution machine.")
    if user.get("role") != "admin":
        st.write(f"{current['name']} · {current['mode']}")
        return
    modes = ["disabled", "peer", "coordinator", "worker"]
    mode = st.selectbox("Workstation connection", modes, index=modes.index(current["mode"]),
                        format_func=lambda v: {"disabled": "Standalone", "peer": "Independent workstations", "coordinator": "Coordinator (legacy)", "worker": "Connected workstation (legacy)"}[v], key="cluster-settings-mode")
    name_key = "cluster-settings-name"
    if st.session_state.get(name_key) == current["name"]:
        st.session_state[name_key] = config.display_name(current["name"])
    name = st.text_input("Workstation name", value=config.display_name(current["name"]), key=name_key)
    address = st.text_input("Local API address" if mode == "peer" else "Coordinator API address", value=current["coordinator_url"] or "http://127.0.0.1:8766", key="cluster-settings-url")
    accept = st.checkbox("Accept new jobs on this workstation", value=current["accept_jobs"], key="cluster-settings-accept")
    st.caption(f"Permanent workstation ID: {config.identity()}")
    if mode == "peer" and current["mode"] != "peer":
        st.info("Connect trusted peers with the cluster enable-peers command. Each installation keeps its accounts, scheduler and results available independently.")
    elif mode == "worker" and current["mode"] != "worker":
        st.info("Enroll this workstation on the coordinator, then import its private enrollment file with the cluster join command.")
        st.code("mn-molecular-design-studio cluster join workstation-enrollment.json --coordinator-url " + address, language="bash")
    elif st.button("Save workstation settings", key="cluster-settings-save"):
        try:
            changes = {"mode": mode, "name": name, "coordinator_url": address, "accept_jobs": accept}
            if mode == "coordinator" and not current["node_token"]:
                from .store import Store
                enrollment = Store().enroll(name, config.identity())
                changes.update(node_token=enrollment["node_token"], default_placement="auto")
            config.save(changes)
            st.success("Workstation settings saved. Start or restart the API and agent services to apply connection changes.")
        except (ValueError, OSError) as exc:
            st.error(str(exc))
    if mode != "disabled":
        st.code("mn-molecular-design-studio cluster install-services --enable", language="bash")
        st.caption("The background agent publishes resources and accepts work while the browser interface is closed.")
    if current["mode"] == "disabled":
        return
    try:
        nodes = _request("/nodes", token)["nodes"]
        st.dataframe(pd.DataFrame([{"Workstation": config.display_name(n["name"]), "Connection": n["availability"], "Last seen": n["last_seen_at"],
                                   "Detected GPUs": _gpu_counts(n)[0],
                                   "Available GPUs": _gpu_counts(n)[1],
                                   "Allocated GPUs": _gpu_counts(n)[2],
                                   "Free CPU threads": n["snapshot"].get("resources", {}).get("cpu_free"),
                                   "Accepts jobs": n["snapshot"].get("accept_jobs", False)} for n in nodes]), hide_index=True, use_container_width=True)
        with st.expander("GPU inventory", expanded=True):
            counts = [_gpu_counts(node) for node in nodes]
            total_detected = sum(item[0] or 0 for item in counts) if all(item[0] is not None for item in counts) else None
            total_available = sum(item[1] for item in counts) if all(item[1] is not None for item in counts) else None
            total_allocated = sum(item[2] or 0 for item in counts) if all(item[2] is not None for item in counts) else None
            detected_label = str(total_detected) if total_detected is not None else "Not reported"
            available_label = str(total_available) if total_available is not None else "Not reported"
            allocated_label = str(total_allocated) if total_allocated is not None else "Not reported"
            st.caption(f"Latest workstation snapshots: {detected_label} detected · {available_label} available · {allocated_label} allocated GPUs.")
            gpu_rows = _gpu_inventory_rows(nodes)
            if gpu_rows:
                st.caption("Available means unreserved in the latest scheduler snapshot. Offline workstations show their last reported state.")
                st.dataframe(gpu_rows, hide_index=True, use_container_width=True)
            else:
                st.info("GPU inventory has not been reported by the connected workstations.")
    except (transport.ClusterUnavailable, ValueError) as exc:
        st.info(str(exc))
    if current["mode"] == "peer":
        st.caption("Accounts and projects are stored locally and synchronized with trusted workstations. Each machine remains usable when the others are offline.")
        from .replication import Replica
        with Replica().db() as db:
            conflicts = db.execute("SELECT key,error FROM conflicts ORDER BY key").fetchall()
        if conflicts:
            st.warning(f"{len(conflicts)} metadata changes need administrator attention. Existing account and job identities have been preserved.")
            st.dataframe([dict(r) for r in conflicts], hide_index=True)
    if current["mode"] == "coordinator":
        with st.expander("Enroll another workstation"):
            new_name = st.text_input("New workstation name", placeholder="pc004032", key="cluster-enroll-name")
            if st.button("Create enrollment file", key="cluster-enroll-create"):
                try:
                    enrollment = _request("/nodes", token, {"name": new_name}, method="POST")
                    st.session_state["cluster-enrollment"] = enrollment
                except (ValueError, transport.ClusterUnavailable) as exc:
                    st.error(str(exc))
            enrollment = st.session_state.get("cluster-enrollment")
            if enrollment:
                st.download_button("Download private enrollment file", json.dumps(enrollment, indent=2),
                                   file_name=f"studio-{enrollment['node_id']}-enrollment.json", mime="application/json")
                st.caption("This file grants access to one workstation agent. Copy it privately to that workstation.")


def render_runs(user, token, local_search, *, elapsed_label=None, elapsed_since=None, submit=None):
    live_updates = bool(st.session_state.get("cluster-results-live-updates", False))
    has_active_jobs = bool(st.session_state.get("cluster-results-has-active-jobs", False))
    run_every = 2 if st.session_state.get("cluster-results-purge-operations") else (
        15 if live_updates and has_active_jobs else None
    )
    @st.fragment(run_every=run_every)
    def live_results():
        _render_runs_content(
            user, token, local_search, elapsed_label=elapsed_label,
            elapsed_since=elapsed_since, submit=submit,
        )

    live_results()


def _rerun_results_app():
    st.rerun(scope="app")


def _render_runs_content(user, token, local_search, *, elapsed_label=None, elapsed_since=None, submit=None):
    st.markdown(
        """<style>
        button[data-testid="stBaseButton-tertiary"],
        button[data-testid="stBaseButton-tertiary"] *,
        [data-testid="stDataEditor"] button,
        [data-testid="stDataEditor"] button * {
            color: #2f6fba !important;
        }
        </style>""",
        unsafe_allow_html=True,
    )
    if st.session_state.pop("cluster-results-reset-selection", False):
        for key in list(st.session_state):
            if str(key).startswith(("cluster-results-select-", "cluster-results-table-", "cluster-results-pipeline-select-")):
                st.session_state.pop(key, None)
        st.session_state["cluster-results-selected"] = {}
    title_col, refresh_col, live_col = st.columns([0.66, 0.17, 0.17], vertical_alignment="center")
    title_col.title("My results")
    live_updates = live_col.toggle(
        "Live updates",
        key="cluster-results-live-updates",
        help="Refresh the job list every 15 seconds while jobs are active.",
        on_change=_rerun_results_app,
    )
    refresh_requested = refresh_col.button("Refresh", icon=":material/refresh:", key="cluster-results-refresh")
    if message := st.session_state.pop("cluster-results-flash", None):
        st.success(message)
    if failures := st.session_state.pop("cluster-results-purge-errors", None):
        for failure in failures:
            st.error(failure)
    purge_operations = list(st.session_state.get("cluster-results-purge-operations") or [])
    purge_poll_finished = False
    if purge_operations:
        remaining_operations = []
        completed_operations = 0
        for operation in purge_operations:
            try:
                result = _request("/commands/" + operation["command_id"], token)
            except (ValueError, transport.ClusterUnavailable) as exc:
                remaining_operations.append(operation)
                st.info(f"Waiting for purge on {operation.get('workstation', 'the owning workstation')}: {exc}")
                continue
            if result.get("status") == "completed":
                completed_operations += 1
            elif result.get("status") == "failed":
                st.error(
                    f"Could not purge {operation.get('name') or 'a selected run'} on "
                    f"{operation.get('workstation', 'its workstation')}: {result.get('error') or 'operation failed'}"
                )
            else:
                remaining_operations.append(operation)
        if remaining_operations:
            st.session_state["cluster-results-purge-operations"] = remaining_operations
            st.info(f"Physically removing {len(remaining_operations)} run folder(s) on the owning workstation(s)…")
        else:
            st.session_state.pop("cluster-results-purge-operations", None)
            purge_poll_finished = True
            st.session_state["cluster-results-force-refresh"] = True
            if completed_operations:
                st.success(f"Physically removed {completed_operations} remote run folder(s) and their automatic Database records.")
    snapshot = st.session_state.get("cluster-results-catalog-snapshot")
    same_user = isinstance(snapshot, dict) and snapshot.get("user_id") == str(user.get("id") or "")
    if not same_user:
        snapshot = None
    refresh_due = bool(
        live_updates
        and snapshot
        and time.monotonic() - float(snapshot.get("refreshed_at") or 0) >= 14.5
    )
    force_refresh = st.session_state.pop("cluster-results-force-refresh", False)
    if refresh_requested or force_refresh or snapshot is None or refresh_due:
        try:
            snapshot = _refresh_results_snapshot(user, token, local_search)
            st.session_state["cluster-results-catalog-snapshot"] = snapshot
        except (OSError, ValueError, transport.ClusterUnavailable, sqlite3.Error) as exc:
            if snapshot is None:
                st.error(f"Could not load jobs: {exc}")
                return
            st.error(f"Refresh failed; showing the last loaded jobs. {exc}")
    rows = list(snapshot.get("items") or [])
    stale = bool(snapshot.get("stale"))
    if stale:
        st.info("Coordinator unavailable. Showing the last synchronized catalog and current local jobs. Remote status is last known.")
    show_removed = st.toggle("Show removed jobs", key="cluster-results-show-removed")
    old_has_active = bool(st.session_state.get("cluster-results-has-active-jobs", False))
    active_statuses = {
        "queued", "pending", "starting", "preparing", "running", "cancelling",
        "cancel_requested", "submitted", "scheduled", "staging", "waiting", "waiting_for_workstation",
        "paused", "holding",
    }
    has_active_jobs = any(str(row.get("status") or "").casefold() in active_statuses for row in rows)
    st.session_state["cluster-results-has-active-jobs"] = has_active_jobs
    if purge_poll_finished or (live_updates and old_has_active != has_active_jobs):
        st.rerun(scope="app")
    if live_updates and not has_active_jobs:
        st.caption("No active jobs to track. Refresh manually to check for new jobs.")
    if snapshot:
        age_seconds = max(0, int(time.monotonic() - float(snapshot.get("refreshed_at") or time.monotonic())))
        age_label = "just now" if age_seconds < 5 else f"{age_seconds}s ago" if age_seconds < 60 else f"{age_seconds // 60}m ago"
        st.caption(f"Job list refreshed {age_label}.")
    status_options = sorted({str(r.get("status") or "unknown").replace("_", " ").title() for r in rows}, key=str.casefold)
    status_key = "cluster-results-status-filter"
    if not set(st.session_state.get(status_key, [])).issubset(status_options):
        st.session_state[status_key] = []
    group_pipeline_runs = st.toggle(
        "Group pipeline steps",
        value=True,
        key="cluster-results-group-pipelines",
        help="Show each pipeline once, including jobs on other computers. Counts include jobs created so far; looping workflows add jobs until they stop. The Steps button opens their tool-call jobs.",
    )
    search_col, machine_col, status_col, sort_col = st.columns([0.38, 0.21, 0.21, 0.20])
    query = search_col.text_input("Search jobs", placeholder="Name, engine, or run ID", key="cluster-results-query").casefold()
    names = {r.get("node_id"): config.display_name(r.get("node_name") or "Not assigned") for r in rows if r.get("node_id")}
    machine = machine_col.selectbox("Workstation", ["", *names], format_func=lambda v: names.get(v, "All workstations"), key="cluster-results-machine")
    statuses = status_col.multiselect("Status", status_options, key=status_key, placeholder="All statuses")
    sort_by = sort_col.selectbox(
        "Sort jobs",
        ["active_desc", "finished_desc", "name_asc", "status_asc"],
        format_func=lambda value: {
            "active_desc": "Newest",
            "finished_desc": "Finished (newest first)",
            "name_asc": "Name (A–Z)",
            "status_asc": "Status (A–Z)",
        }[value],
        key="cluster-results-sort-by-v2",
    )
    from mn_molecular_design_studio.workspace import load_organization
    project_records = [p for p in load_organization().get("projects", []) if isinstance(p, dict) and p.get("id")]
    known_project_ids = {str(p["id"]) for p in project_records}
    for row in rows:
        project_id = str(row.get("project_id") or "")
        if project_id and project_id not in known_project_ids:
            project_records.append({"id": project_id, "name": f"Project {project_id[:8]}"})
            known_project_ids.add(project_id)
    project_names = {str(p["id"]): str(p.get("name") or "Unnamed project") for p in project_records}
    project_activity = {}
    for row in rows:
        project_id = str(row.get("project_id") or "")
        updated_at = str(row.get("updated_at") or row.get("completed_at") or row.get("created_at") or "")
        if project_id and updated_at > project_activity.get(project_id, ""):
            project_activity[project_id] = updated_at
    project_col, project_order_col = st.columns([0.72, 0.28])
    project_order = project_order_col.selectbox(
        "Project order",
        ["modified_desc", "created_desc", "name_asc"],
        format_func=lambda value: {
            "modified_desc": "Last modified (newest first)",
            "created_desc": "Created (newest first)",
            "name_asc": "Name (A–Z)",
        }[value],
        key="cluster-results-project-order",
        help="Last modified is based on the newest updated job in each project; projects without jobs fall back to their creation date.",
    )
    if project_order == "name_asc":
        project_records.sort(key=lambda p: (str(p.get("name") or "").casefold(), str(p.get("id") or "")))
    elif project_order == "created_desc":
        project_records.sort(key=lambda p: (str(p.get("name") or "").casefold(), str(p.get("id") or "")))
        project_records.sort(key=lambda p: str(p.get("created_at") or ""), reverse=True)
    else:
        project_records.sort(key=lambda p: (str(p.get("name") or "").casefold(), str(p.get("id") or "")))
        project_records.sort(
            key=lambda p: max(
                str(p.get("updated_at") or ""),
                project_activity.get(str(p.get("id") or ""), ""),
                str(p.get("created_at") or ""),
            ),
            reverse=True,
        )
    project_options = ["", *[str(p["id"]) for p in project_records]]
    project = project_col.selectbox(
        "Project",
        project_options,
        format_func=lambda value: project_names.get(value, "All projects"),
        key="cluster-results-project",
    )
    filters_key = (query, machine, project, tuple(sorted(statuses)), sort_by, show_removed, group_pipeline_runs)
    if st.session_state.get("cluster-results-filter-state") != filters_key:
        st.session_state["cluster-results-filter-state"] = filters_key
        st.session_state["cluster-results-page"] = 1
    scope_rows = [r for r in rows if (show_removed or not r.get("hidden"))
                  and (not machine or r.get("node_id") == machine)
                  and (not project or r.get("project_id") == project)]
    chosen = [r for r in scope_rows
              if (not statuses or str(r.get("status") or "unknown").replace("_", " ").title() in statuses)
              and (not query or query in " ".join(f"{r.get(k, '')} {config.display_text(r.get(k, ''))}" for k in ("campaign_name", "pipeline_name", "run_id", "tool", "job_code")).casefold())]
    page_sort_key = "cluster-results-last-sort"
    if st.session_state.get(page_sort_key) != sort_by:
        st.session_state[page_sort_key] = sort_by
        st.session_state["cluster-results-page"] = 1
    if sort_by == "name_asc":
        chosen.sort(
            key=lambda r: (
                str(r.get("campaign_name") or r.get("tool") or r.get("job_type") or "Job").casefold(),
                str(r.get("created_at") or ""),
            )
        )
    elif sort_by == "status_asc":
        chosen.sort(
            key=lambda r: (
                str(r.get("status") or "unknown").replace("_", " ").casefold(),
                str(r.get("created_at") or ""),
            )
        )
    elif sort_by == "active_desc":
        finished_statuses = {"completed", "succeeded", "failed", "cancelled", "canceled"}
        active_rows = [r for r in chosen if str(r.get("status") or "unknown").casefold() not in finished_statuses]
        finished_rows = [r for r in chosen if str(r.get("status") or "unknown").casefold() in finished_statuses]
        def newest_key(row):
            return str(row.get("updated_at") or row.get("completed_at") or row.get("created_at") or "")

        active_rows.sort(key=newest_key, reverse=True)
        finished_rows.sort(key=newest_key, reverse=True)
        chosen = active_rows + finished_rows
    else:
        finished_statuses = {"completed", "succeeded", "failed", "cancelled", "canceled"}
        chosen.sort(
            key=lambda r: (
                str(r.get("status") or "unknown").casefold() in finished_statuses,
                str(r.get("completed_at") or r.get("updated_at") or r.get("created_at") or ""),
                str(r.get("created_at") or ""),
            ),
            reverse=True,
        )
    selected_jobs = st.session_state.setdefault("cluster-results-selected", {})
    row_by_ref = {_results_job_reference(row): row for row in rows}
    for reference in list(selected_jobs):
        if reference in row_by_ref:
            selected_jobs[reference] = row_by_ref[reference]
        else:
            selected_jobs.pop(reference, None)

    pipeline_groups, row_to_group = _visual_pipeline_groups(rows)
    if group_pipeline_runs:
        for root_ref, members in pipeline_groups.items():
            if root_ref in selected_jobs:
                selected_jobs.update((_results_job_reference(member), member) for member in members)
    display_entries = []
    seen_groups = set()
    for row in chosen:
        reference = _results_job_reference(row)
        group_ref = row_to_group.get(reference) if group_pipeline_runs else None
        if group_ref:
            if group_ref in seen_groups:
                continue
            seen_groups.add(group_ref)
            members = pipeline_groups[group_ref]
            display_entries.append({
                "kind": "pipeline",
                "root": members[0],
                "members": members,
            })
        else:
            display_entries.append({"kind": "job", "row": row})

    grouped_count = sum(entry["kind"] == "pipeline" for entry in display_entries)
    caption = f"{len(chosen)} jobs · Job status and file availability are tracked separately."
    if grouped_count:
        caption += f" {grouped_count} pipeline run{'s' if grouped_count != 1 else ''} grouped."
    st.caption(caption)
    page_size_col, page_col, range_col = st.columns([0.24, 0.24, 0.52], vertical_alignment="bottom")
    page_size = page_size_col.selectbox("Entries per page", [25, 50, 100], index=1, key="cluster-results-page-size")
    page_size_key = "cluster-results-page-size-last"
    if st.session_state.get(page_size_key) != page_size:
        st.session_state[page_size_key] = page_size
        st.session_state["cluster-results-page"] = 1
    pages = max(1, (len(display_entries) + page_size - 1) // page_size)
    page_key = "cluster-results-page"
    if st.session_state.get(page_key, 1) > pages:
        st.session_state[page_key] = pages
    selected_page = int(page_col.number_input("Page", min_value=1, max_value=pages, value=1, step=1, key=page_key))
    first_index = (selected_page - 1) * page_size
    last_index = min(first_index + page_size, len(display_entries))
    range_col.caption(f"{first_index + 1 if display_entries else 0}–{last_index} of {len(display_entries)}")
    page_entries = display_entries[first_index:last_index]

    if not page_entries:
        st.info("No jobs match your search and filters.")
    else:
        entry_worker_details = []
        for entry in page_entries:
            row = entry["root"] if entry["kind"] == "pipeline" else entry["row"]
            members = entry["members"] if entry["kind"] == "pipeline" else [row]
            entry_worker_details.append(_library_worker_details(row, members))
        show_computers = any(entry_worker_details)
        table_rows = []
        entry_members = []
        entry_pipeline_refs = []
        entry_view_refs = []
        entry_result_refs = []
        for entry_index, entry in enumerate(page_entries):
            is_pipeline = entry["kind"] == "pipeline"
            row = entry["root"] if is_pipeline else entry["row"]
            ref = _results_job_reference(row)
            members = entry["members"] if is_pipeline else [row]
            entry_members.append(members)
            entry_pipeline_refs.append(ref if is_pipeline else None)
            entry_view_refs.append(_job_details_reference(row))
            entry_result_refs.append(_job_details_reference(row) if is_pipeline else None)
            worker_details = entry_worker_details[entry_index]
            availability = "Last known" if stale and row.get("node_id") != config.identity() else str(row.get("availability", "pending")).title()
            if row.get("last_seen_at"):
                availability += " · " + str(row["last_seen_at"])
            table_row = {
                "_reference": ref,
                "Select": all(_results_job_reference(member) in selected_jobs for member in members),
                "Job": config.display_text(row.get("campaign_name") or row.get("tool") or "Job"),
                "View": "View",
            }
            if group_pipeline_runs:
                table_row["Results"] = ("Final results" if str(row.get("status") or "").casefold() in {"completed", "succeeded"}
                                        else "Results so far") if is_pipeline else None
            if show_computers:
                table_row["Computers"] = f"{len(worker_details['workers'])} computers" if worker_details else None
            if group_pipeline_runs:
                children = members[1:] if is_pipeline else []
                table_row["Steps"] = _pipeline_step_count_label(row, children) if is_pipeline else None
            table_row.update({
                "Project": project_names.get(row.get("project_id"), "Project"),
                "Tool": display_tool_label(row),
            })
            if not group_pipeline_runs:
                table_row["Pipeline"] = config.display_text(str(row.get("pipeline_name") or "").strip() or "—")
            table_row.update({
                "Code": str(row.get("job_code") or ""),
                "Status": str(row.get("status") or "unknown").replace("_", " ").title(),
                "Elapsed": elapsed_label(row) if elapsed_label else "—",
                "Updated": elapsed_since(row) if elapsed_since else "",
                "Workstation": config.display_name(row.get("node_name") or "Not assigned"),
                "Files": availability,
                "Progress": _progress_summary(row.get("progress")),
            })
            table_rows.append(table_row)
        table = pd.DataFrame(table_rows).set_index("_reference")
        table_key = "cluster-results-table-" + hashlib.sha1(
            repr((query, machine, project, tuple(sorted(statuses)), sort_by, show_removed,
                  group_pipeline_runs, selected_page, page_size, tuple(table.index), tuple(table.columns),
                  tuple(sorted(selected_jobs)))).encode("utf-8")
        ).hexdigest()[:12]
        column_config = {
            "Select": st.column_config.CheckboxColumn(
                "Select", width="small",
                help="Selecting a pipeline includes its parent run and every tool-call job. Finish or cancel jobs before purging them.",
            ),
            "Job": st.column_config.TextColumn("Job", width="large"),
            "View": st.column_config.ButtonColumn(
                "View", width="small", type="tertiary", key=table_key + "-view",
                on_click=_request_results_job_view,
                args=(table_key + "-view", entry_view_refs),
            ),
            "Project": st.column_config.TextColumn("Project", width="medium"),
            "Tool": st.column_config.TextColumn("Tool", width="medium"),
            "Pipeline": st.column_config.TextColumn("Pipeline", width="medium"),
            "Code": st.column_config.TextColumn("#", width="small"),
            "Status": st.column_config.TextColumn("Status", width="small"),
            "Elapsed": st.column_config.TextColumn("Elapsed", width="small"),
            "Updated": st.column_config.TextColumn("Updated", width="medium"),
            "Workstation": st.column_config.TextColumn("Workstation", width="medium"),
            "Files": st.column_config.TextColumn("Files", width="medium"),
            "Progress": st.column_config.TextColumn("Progress", width="large"),
        }
        if show_computers:
            worker_click_key = table_key + "-computers"
            column_config["Computers"] = st.column_config.ButtonColumn(
                "Computers", width="small", type="tertiary",
                help="Open the distributed worker computers and inspect each computer's processing steps.",
                key=worker_click_key,
                on_click=_request_results_worker_details,
                args=(worker_click_key, entry_worker_details),
            )
        if group_pipeline_runs:
            results_key = table_key + "-results"
            column_config["Results"] = st.column_config.ButtonColumn("Results", width=120, type="tertiary",
                help="Open the combined campaign results across all iterations and batches.", key=results_key,
                on_click=_request_results_pipeline_results, args=(results_key, entry_result_refs))
            click_key = table_key + "-steps"
            column_config["Steps"] = st.column_config.ButtonColumn(
                "Steps", width=120, type="tertiary",
                help="Open pipeline iterations, then drill into their batch and neural-network jobs. Includes jobs on every computer.",
                key=click_key,
                on_click=_request_results_pipeline_steps,
                args=(click_key, entry_pipeline_refs),
            )
        styled_table = table.style.map(
            lambda _: "color: #2f6fba; font-weight: 600;",
            subset=[column for column in ("View", "Results", "Computers", "Steps") if column in table.columns],
        )
        st.data_editor(
            styled_table,
            key=table_key,
            hide_index=True,
            num_rows="fixed",
            disabled=[column for column in table.columns if column not in {"Select", "Steps", "View", "Results", "Computers"}],
            column_config=column_config,
            on_change=_update_results_table_selection,
            args=(table_key, entry_members),
            width="stretch",
            row_height=36,
            height=min(680, len(table_rows) * 36 + 42),
        )
    details_reference = st.session_state.pop("cluster-results-pipeline-details-request", None)
    if details_reference in pipeline_groups:
        _render_results_pipeline_steps(
            pipeline_groups[details_reference], elapsed_label, user, token, submit,
            iteration=st.session_state.pop("cluster-results-pipeline-iteration-request", None),
        )
    worker_details_request = st.session_state.pop("cluster-results-library-workers-request", None)
    if worker_details_request:
        _render_library_worker_details(worker_details_request)
    st.session_state["cluster-results-selected"] = selected_jobs
    selected_rows = list(selected_jobs.values())
    if selected_rows:
        st.caption(f"{len(selected_rows)} job{'s' if len(selected_rows) != 1 else ''} selected")
        if show_removed:
            available_to_restore = [row for row in selected_rows if row.get("hidden")]
            if available_to_restore and st.button("Restore selected jobs", key="cluster-results-restore-selected"):
                failures = []
                restored = 0
                for row in available_to_restore:
                    try:
                        jid = str(row.get("cluster_job_id") or "")
                        if jid:
                            _request("/jobs/" + jid + "/actions", token, {"action": "restore"}, method="POST")
                        elif row.get("node_id") == config.identity() and row.get("task_group") and row.get("run_id"):
                            from mn_molecular_design_studio.auth import can_access_job
                            from mn_molecular_design_studio.job_catalog import set_job_hidden, get_job_access_record
                            access_record = get_job_access_record(str(row["task_group"]), str(row["run_id"]))
                            if not can_access_job(user, access_record):
                                raise PermissionError("Job not found.")
                            set_job_hidden(str(row["task_group"]), str(row["run_id"]), False)
                        else:
                            raise ValueError("This job cannot be updated while its workstation is unavailable.")
                        restored += 1
                    except (ValueError, PermissionError, transport.ClusterUnavailable) as exc:
                        failures.append(str(exc))
                st.session_state["cluster-results-selected"] = {}
                st.session_state["cluster-results-reset-selection"] = True
                st.session_state["cluster-results-force-refresh"] = True
                if restored:
                    st.session_state["cluster-results-flash"] = f"Restored {restored} job{'s' if restored != 1 else ''}."
                for failure in failures:
                    st.error(failure)
                st.rerun()
        else:
            pending_key = "cluster-results-pending-remove-and-purge"
            active_selected = any(
                str(row.get("status") or "unknown").casefold() not in _RESULTS_TERMINAL_STATUSES
                for row in selected_rows
            )
            pending = st.session_state.get(pending_key)
            if pending and (active_selected or {_results_job_reference(row) for row in pending} != set(selected_jobs)):
                st.session_state.pop(pending_key, None)
            if active_selected:
                st.info("Finish or cancel the selected jobs before removing them. A pipeline selection includes all its steps.")
            if st.button("Remove and purge selected", key="cluster-results-remove-and-purge-request", disabled=active_selected):
                pending_rows = [dict(row) for row in selected_rows]
                pending_rows.sort(key=_is_visual_pipeline_row)
                st.session_state[pending_key] = pending_rows
            pending = st.session_state.get(pending_key)
            if pending:
                st.warning("Confirming permanently deletes each selected finished run folder from its workstation's Results directory and removes its automatic Database records. Queued or running jobs cannot be purged. Molecules explicitly saved to user collections are kept. Remote runs are deleted when their workstation receives the request.")
                confirm, cancel = st.columns(2)
                if confirm.button("Confirm removal and purge", type="primary", key="cluster-results-remove-and-purge-confirm"):
                    removed_sources = removed_members = succeeded = 0
                    queued_operations = []
                    failures = []
                    for row in pending:
                        try:
                            jid = str(row.get("cluster_job_id") or "")
                            if jid:
                                result = _request("/jobs/" + jid + "/purge", token, method="DELETE")
                                if result.get("status") == "purge_queued":
                                    queued_operations.append({
                                        "command_id": result["command_id"],
                                        "workstation": config.display_name(row.get("node_name") or result.get("node_id")),
                                        "name": row.get("campaign_name") or row.get("run_id") or jid,
                                    })
                                else:
                                    database = result.get("database") or {}
                                    removed_sources += int(database.get("removed_sources", 0))
                                    removed_members += int(database.get("removed_members", 0))
                                    succeeded += 1
                            elif row.get("node_id") == config.identity() and row.get("task_group") and row.get("run_id"):
                                from mn_molecular_design_studio.auth import can_access_job
                                from mn_molecular_design_studio.job_catalog import get_job_access_record, reconcile
                                from mn_molecular_design_studio.molecular_database_results import remove_run_records
                                from mn_molecular_design_studio.engine.workflow_engine.core.jobs import delete_job_run, runs_root
                                access_record = get_job_access_record(str(row["task_group"]), str(row["run_id"]))
                                if not can_access_job(user, access_record):
                                    raise PermissionError("Job not found.")
                                if str(row.get("status") or "").casefold() not in {"completed", "succeeded", "cancelled", "failed"}:
                                    raise ValueError("Only finished jobs can be permanently removed.")
                                result = remove_run_records(
                                    row["task_group"], row["run_id"], user,
                                    owner_id=row.get("owner_id"), project_id=row.get("project_id"),
                                )
                                removed_sources += result["removed_sources"]
                                removed_members += result["removed_members"]
                                try:
                                    delete_job_run(row["task_group"], row["run_id"])
                                except FileNotFoundError:
                                    pass
                                reconcile(runs_root())
                                succeeded += 1
                            else:
                                raise ValueError("This job cannot be updated while its workstation is unavailable.")
                        except (ValueError, PermissionError, OSError, sqlite3.Error, transport.ClusterUnavailable) as exc:
                            failures.append(f"{row.get('campaign_name') or row.get('run_id') or 'Selected job'}: {exc}")
                    st.session_state.pop(pending_key, None)
                    st.session_state["cluster-results-selected"] = {}
                    st.session_state["cluster-results-reset-selection"] = True
                    st.session_state["cluster-results-force-refresh"] = True
                    if queued_operations:
                        st.session_state["cluster-results-purge-operations"] = queued_operations
                    messages = []
                    if succeeded:
                        messages.append(
                            f"Physically removed {succeeded} run folder{'s' if succeeded != 1 else ''} "
                            f"and their automatic Database records "
                            f"({removed_sources} source links, {removed_members} run records)."
                        )
                    if queued_operations:
                        messages.append(
                            f"Queued physical deletion of {len(queued_operations)} run folder"
                            f"{'s' if len(queued_operations) != 1 else ''} on the owning workstation(s)."
                        )
                    if messages:
                        st.session_state["cluster-results-flash"] = " ".join(messages)
                    if failures:
                        st.session_state["cluster-results-purge-errors"] = failures
                    st.rerun()
                if cancel.button("Cancel", key="cluster-results-remove-and-purge-cancel"):
                    st.session_state.pop(pending_key, None)
                    st.rerun()
        clear_col, _ = st.columns([1, 4])
        if clear_col.button("Clear selection", key="cluster-results-clear-selection"):
            st.session_state["cluster-results-selected"] = {}
            st.session_state["cluster-results-reset-selection"] = True
            st.rerun()


def render_job(reference, user, token, submit):
    jid = reference.removeprefix("cluster:")
    try:
        row = _request("/jobs/" + jid, token)
    except transport.ClusterUnavailable:
        cached, _ = catalog(user, token)
        row = next((j for j in cached["items"] if j.get("cluster_job_id") == jid), None)
        if not row:
            st.info("The coordinator and this job's cached record are unavailable.")
            return
        row = dict(row, availability="last known")
    if (st.session_state.get("pipeline-final-results-reference") == reference and _is_visual_pipeline_row(row)
            and "unique_docked_compounds" in (row.get("metrics") or {})):
        _render_remote_pipeline_results(jid, row, token)
        return
    st.title(config.display_text(row.get("campaign_name") or row.get("tool") or "Job"))
    st.caption(f"Tool: {display_tool_label(row)} · Pipeline: {config.display_text(str(row.get('pipeline_name') or '').strip() or '—')}")
    st.caption(f"Files stored on {config.display_name(row.get('node_name') or 'Not assigned')} · {row.get('location') or 'Waiting for placement'}")
    a, b, c = st.columns(3)
    a.metric("Job status", row.get("status", "unknown").replace("_", " ").title())
    b.metric("Workstation", config.display_name(row.get("node_name") or "Not assigned"))
    c.metric("Data availability", row.get("availability", "unknown").title())
    if _is_bindcraft2_job(row):
        _render_bindcraft2_progress(row, str(row.get("cluster_job_id") or jid), token)
    if row.get("availability") != "online":
        st.info("The saved job record remains available. Full files and live analysis will be available when the owning workstation reconnects.")
    if row.get("wait_reason"):
        st.info(row["wait_reason"])
    if row.get("error"):
        st.error(row["error"])
    transfer = row.get("library_transfer")
    if transfer:
        transferred = int(transfer.get("transferred_bytes") or 0)
        total = int(transfer.get("total_bytes") or 0)
        state = str(transfer.get("status") or "waiting").replace("_", " ").title()
        st.caption(f"Fingerprint outputs to organizer: {state} · {transferred:,} / {total:,} bytes")
        if total:
            st.progress(min(1.0, transferred / total))
        if transfer.get("error"):
            st.warning(transfer["error"])
    if row.get("status") not in {"completed", "succeeded", "failed", "cancelled"}:
        if st.button("Request cancellation", key="cluster-cancel-" + jid):
            response = _request("/jobs/" + jid + "/actions", token, {"action": "cancel"}, method="POST")
            st.success("Cancellation requested. The workstation will confirm when execution has stopped.")
    if row.get("metrics"):
        st.json(row.get("metrics", {}), expanded=False)
    acceptance_counts = (row.get("candidate_summary") or {}).get("acceptance_counts") or {}
    accepted_count = int(acceptance_counts.get("accepted") or 0)
    if accepted_count and row.get("local_reference") and row.get("availability") == "online":
        st.info(
            f"Validation explicitly accepted {accepted_count:,} design(s). Their structures remain on "
            f"{config.display_name(row.get('node_name') or 'the source workstation')} until you collect them."
        )
        if st.button("Collect accepted designs into shared Database", type="primary", key="cluster-collect-accepted-" + jid):
            try:
                operation = _request(
                    "/jobs/" + jid + "/actions", token,
                    {"action": "collect_accepted_designs"}, method="POST",
                )
                st.session_state["cluster-accepted-collect-op-" + jid] = operation["command_id"]
            except (ValueError, transport.ClusterUnavailable) as exc:
                st.error(f"Could not request accepted-design collection: {exc}")

        @st.fragment(run_every=2)
        def accepted_collection_operation():
            command_id = st.session_state.get("cluster-accepted-collect-op-" + jid)
            if not command_id:
                return
            try:
                operation = _request("/commands/" + command_id, token)
                if operation["status"] == "failed":
                    st.error(operation.get("error") or "Accepted-design collection failed.")
                    return
                if operation["status"] != "completed":
                    st.info("The source workstation is packaging accepted designs and their structures.")
                    return
                response = operation.get("response") or {}
                bytes_key = "cluster-accepted-bundle-" + command_id
                if bytes_key not in st.session_state:
                    st.session_state[bytes_key] = _request(
                        "/commands/" + command_id + "/file", token, raw=True, timeout=180,
                    )
                if "cluster-accepted-promotion-" + command_id not in st.session_state:
                    from mn_molecular_design_studio.accepted_design_import import promote_accepted_design_bundle
                    promoted = promote_accepted_design_bundle(
                        st.session_state[bytes_key], user=user, project_id=str(row.get("project_id") or ""),
                    )
                    st.session_state["cluster-accepted-promotion-" + command_id] = promoted
                promoted = st.session_state["cluster-accepted-promotion-" + command_id]
                duplicate = " (already collected)" if promoted.get("already_imported") else ""
                st.success(
                    f"Collected {promoted['candidate_count']:,} accepted design(s) into the shared Database{duplicate}."
                )
            except (ValueError, OSError, zipfile.BadZipFile, transport.ClusterUnavailable) as exc:
                st.error(f"Could not finish accepted-design collection: {exc}")

        accepted_collection_operation()
    if str(row.get("task_group") or "") in {"ligand-md", "ligand-mmgbsa", "md-analysis"} and row.get("availability") == "online":
        from mn_molecular_design_studio.md_remote_ui import render_remote_md_results
        render_remote_md_results(st, jid, row, token, _fetch_remote_artifact_text)
    if row.get("local_reference"):
        if st.button("Load full results from workstation", key="cluster-inspect-" + jid):
            operation = _request("/jobs/" + jid + "/actions", token, {"action": "inspect"}, method="POST")
            st.session_state["cluster-details-op-" + jid] = operation["command_id"]
        @st.fragment(run_every=2)
        def inspect_operation():
            command_id = st.session_state.get("cluster-details-op-" + jid)
            if not command_id:
                return
            try:
                op = _request("/commands/" + command_id, token)
            except (ValueError, transport.ClusterUnavailable) as exc:
                st.info(str(exc)); return
            if op["status"] == "failed":
                st.error(op["error"]); return
            if op["status"] != "completed":
                st.info("Waiting for the owning workstation to return the results."); return
            data = op["response"]
            from mn_molecular_design_studio.distributed_campaign_ui import render_campaign_context
            render_campaign_context(st, data.get("distributed_campaign_context"))
            candidates = data.get("results", {}).get("candidates", [])
            if candidates:
                tool = " ".join(
                    str(row.get(key) or "")
                    for key in ("tool", "job_type", "workflow", "pipeline_name", "campaign_name")
                ).casefold() + " " + display_tool_label(row).casefold()
                if "rfpeptides" in tool:
                    st.info("RFpeptides generated cyclic peptide backbones. It does not design amino-acid sequences; run ProteinMPNN or another sequence-design tool as a follow-up.")
                rows = _candidate_result_rows(candidates, tool=tool)
                if "bindcraft2" in tool:
                    st.caption("The table includes BindCraft2's ranked confidence, interface, clash, and developability metrics. Scroll horizontally to see all columns.")
                st.dataframe(
                    pd.DataFrame(rows),
                    hide_index=True,
                    use_container_width=True,
                    height=min(410, 38 * (len(rows) + 1) + 3),
                    column_config={"Sequence": st.column_config.TextColumn(width="large")},
                )
                parameter_rows, filter_rows = ([], [])
                if "bindcraft2" in tool:
                    parameter_rows, filter_rows = _bindcraft2_run_parameter_rows(data)
                selected = st.selectbox(
                    "Design to visualize",
                    range(len(candidates)),
                    format_func=lambda index: f"{index + 1}. {candidates[index].get('candidate_id') or 'Design'}",
                    key="cluster-candidate-" + jid,
                )
                candidate = candidates[selected]
                sequence = candidate.get("binder_sequence") or candidate.get("sequence")
                if sequence:
                    st.markdown("**Designed molecule sequence**")
                    st.code(sequence, language=None)
                elif "rfpeptides" in tool:
                    st.caption("Sequence: not applicable yet; this run contains a backbone structure only.")

                structure = _candidate_structure_path(candidate, data)
                if not structure:
                    st.info("No candidate structure file was recorded for this design.")
                elif row.get("availability") != "online":
                    st.info("The structure is saved on the owning workstation and will load when that workstation reconnects.")
                else:
                    relative_path, structure_field = structure
                    candidate_key = str(candidate.get("candidate_id") or selected)
                    file_key = f"cluster-candidate-file-{jid}-{candidate_key}"
                    command_id = st.session_state.get(file_key)
                    if not command_id:
                        try:
                            file_operation = _request(
                                "/jobs/" + jid + "/actions",
                                token,
                                {"action": "artifact", "parameters": {"path": relative_path}},
                                method="POST",
                            )
                            command_id = file_operation["command_id"]
                            st.session_state[file_key] = command_id
                        except (ValueError, transport.ClusterUnavailable) as exc:
                            st.warning(f"Could not request the representative structure: {exc}")
                    if command_id:
                        try:
                            file_status = _request("/commands/" + command_id, token)
                            if file_status["status"] == "failed":
                                st.warning(file_status.get("error") or "The workstation could not serve this structure.")
                            elif file_status["status"] != "completed":
                                st.info("Loading the representative structure from the owning workstation…")
                            else:
                                metadata = file_status.get("response") or {}
                                name = metadata.get("name") or Path(relative_path).name
                                if metadata.get("size_bytes", 0) > 50 * 1024 * 1024:
                                    st.info(f"{name} is too large for the inline structure viewer. Use the saved-file controls below to download it.")
                                else:
                                    content_key = "cluster-candidate-bytes-" + command_id
                                    if content_key not in st.session_state:
                                        st.session_state[content_key] = _request(
                                            "/commands/" + command_id + "/file",
                                            token,
                                            raw=True,
                                            timeout=60,
                                        )
                                    content = st.session_state[content_key]
                                    structure_text = content.decode("utf-8", errors="replace")
                                    st.subheader("Representative structure")
                                    from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer import molstar_custom_component
                                    from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer.dataclasses import StructureVisualization
                                    visualizations = [StructureVisualization(
                                        pdb=structure_text,
                                        color=_candidate_structure_color(candidate, structure_text),
                                        representation_type="cartoon",
                                    )]
                                    overlay_key = None
                                    if "bindcraft2" in tool:
                                        candidate_token = hashlib.sha256(candidate_key.encode()).hexdigest()[:12]
                                        overlay_key = f"cluster-bc2-overlay-{jid}-{candidate_token}"
                                        if st.checkbox(
                                            "Overlay binder-alone monomer",
                                            key=overlay_key + "-enabled",
                                            help="Align the saved binder-only prediction onto the binder chain in this complex.",
                                        ):
                                            monomer_path = _bindcraft2_monomer_path("", candidate)
                                            if not monomer_path:
                                                refolded_csv = _fetch_remote_artifact_text(
                                                    jid,
                                                    _bindcraft2_refolded_table_path(data),
                                                    token,
                                                    overlay_key + "-refolded-csv",
                                                    label="BindCraft2 refolded results table",
                                                )
                                                if refolded_csv:
                                                    monomer_path = _bindcraft2_monomer_path(refolded_csv, candidate)
                                            if monomer_path:
                                                monomer_token = hashlib.sha256(monomer_path.encode()).hexdigest()[:12]
                                                monomer_text = _fetch_remote_artifact_text(
                                                    jid,
                                                    monomer_path,
                                                    token,
                                                    overlay_key + "-monomer-" + monomer_token,
                                                    label="binder-alone monomer",
                                                )
                                                if monomer_text:
                                                    aligned_key = overlay_key + "-aligned-" + monomer_token
                                                    if aligned_key not in st.session_state:
                                                        try:
                                                            aligned_text, overlay_rmsd, matched_residues = _align_bindcraft2_monomer_to_complex(
                                                                structure_text, monomer_text, candidate
                                                            )
                                                            st.session_state[aligned_key] = {
                                                                "text": aligned_text,
                                                                "rmsd": overlay_rmsd,
                                                                "matched_residues": matched_residues,
                                                            }
                                                        except (ImportError, KeyError, ValueError, OSError) as exc:
                                                            st.warning(f"Could not align the binder-alone monomer: {exc}")
                                                    aligned = st.session_state.get(aligned_key)
                                                    if aligned:
                                                        visualizations.append(StructureVisualization(
                                                            pdb=aligned["text"],
                                                            color="uniform",
                                                            color_params={"value": "0xff7d45"},
                                                            representation_type="cartoon",
                                                        ))
                                                        st.caption(
                                                            "Orange: binder-alone monomer, fitted onto the complex binder. "
                                                            f"Cα RMSD: {aligned['rmsd']:.2f} Å over {aligned['matched_residues']} matched residues."
                                                        )
                                            else:
                                                st.warning("No matching binder-alone monomer file was found for this accepted sequence.")
                                    molstar_custom_component(
                                        visualizations,
                                        key=f"cluster-candidate-view-{jid}-{selected}-{structure_field}-overlay-{bool(overlay_key and overlay_key + '-enabled' in st.session_state and st.session_state[overlay_key + '-enabled'])}",
                                        height=520,
                                        show_controls=True,
                                        show_structure_tools=False,
                                    )
                                    st.download_button(
                                        "Download candidate structure",
                                        content,
                                        file_name=name,
                                        key=f"cluster-candidate-download-{command_id}",
                                    )
                        except (ValueError, transport.ClusterUnavailable) as exc:
                            st.info(str(exc))
                if parameter_rows:
                    st.markdown("**BindCraft2 run parameters**")
                    st.dataframe(
                        pd.DataFrame(parameter_rows),
                        hide_index=True,
                        use_container_width=True,
                        height=min(360, 35 * (len(parameter_rows) + 1) + 3),
                    )
                if filter_rows:
                    st.markdown("**Quality filters used**")
                    st.dataframe(
                        pd.DataFrame(filter_rows),
                        hide_index=True,
                        use_container_width=True,
                        height=min(260, 35 * (len(filter_rows) + 1) + 3),
                    )
            _render_remote_docking_preview(jid, row, data, token)
            with st.expander("Run settings and inputs"):
                st.json(data.get("input", {}))
            with st.expander("Run summary"):
                st.json(data.get("result", {}))
            if candidates and "rfantibody" in tool:
                from mn_molecular_design_studio.rfantibody_ui import render_confidence_guide
                render_confidence_guide()
        inspect_operation()
        artifacts = row.get("artifacts", [])
        if artifacts:
            selected = st.selectbox("Saved file", range(len(artifacts)), format_func=lambda i: artifacts[i].get("relative_path", artifacts[i].get("name", "File")), key="cluster-file-" + jid)
            if st.button("Prepare file for viewing or download", key="cluster-file-load-" + jid):
                path = artifacts[selected].get("relative_path")
                operation = _request("/jobs/" + jid + "/actions", token, {"action": "artifact", "parameters": {"path": path}}, method="POST")
                st.session_state["cluster-file-op-" + jid] = operation["command_id"]
            @st.fragment(run_every=2)
            def file_operation():
                command_id = st.session_state.get("cluster-file-op-" + jid)
                if not command_id:
                    return
                try:
                    operation = _request("/commands/" + command_id, token)
                    if operation["status"] == "failed":
                        st.error(operation["error"]); return
                    if operation["status"] != "completed":
                        st.info("Waiting for the owning workstation to serve this file."); return
                    key = "cluster-file-bytes-" + command_id
                    if key not in st.session_state:
                        st.session_state[key] = _request("/commands/" + command_id + "/file", token, raw=True, timeout=60)
                    content = st.session_state[key]
                    name = operation["response"].get("name", "result")
                    st.download_button("Download " + name, content, file_name=name, key="cluster-download-" + command_id)
                    if Path(name).suffix.lower() in {".pdb", ".cif", ".mmcif"} and len(content) < 50 * 1024 * 1024:
                        from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer import molstar_custom_component
                        from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer.dataclasses import StructureVisualization
                        molstar_custom_component([StructureVisualization(pdb=content.decode("utf-8", errors="replace"), color="plddt")], key="cluster-viewer-" + command_id, height=520, show_controls=False, show_structure_tools=False)
                except (ValueError, transport.ClusterUnavailable) as exc:
                    st.info(str(exc))
            file_operation()
        with st.expander("Continue on the workstation that owns these results"):
            workflow = st.selectbox("Next step", ["analysis.ipsae", "analysis.rosetta", "sequence.mpnn", "refolding.boltz2_complex"],
                                    format_func=lambda v: {"analysis.ipsae": "ipSAE interface confidence", "analysis.rosetta": "Rosetta interface analysis", "sequence.mpnn": "ProteinMPNN sequence redesign", "refolding.boltz2_complex": "Boltz-2 complex prediction"}[v], key="cluster-analysis-" + jid)
            if workflow == "analysis.ipsae":
                parameters = {"pae_cutoff": st.number_input("PAE cutoff", 1.0, 40.0, 10.0, key="cluster-pae-" + jid),
                              "dist_cutoff": st.number_input("Distance cutoff", 1.0, 40.0, 10.0, key="cluster-dist-" + jid)}
            elif workflow == "analysis.rosetta":
                parameters = {"nprocs": st.number_input("CPU processes", 1, 64, 4, key="cluster-nprocs-" + jid)}
            elif workflow == "sequence.mpnn":
                parameters = {"model_type": st.selectbox("MPNN model", ["protein_mpnn", "ligand_mpnn", "soluble_mpnn", "hyper_mpnn", "ab_mpnn"], key="cluster-mpnn-" + jid),
                              "num_seq_per_target": st.number_input("Sequences per design", 1, 1000, 2, key="cluster-sequences-" + jid),
                              "fix_interface_residues": st.checkbox("Keep interface residues when a valid source sequence is available", key="cluster-interface-" + jid)}
            else:
                parameters = {"require_monomer_success": False, "recycling_steps": st.number_input("Recycles", 1, 20, 3, key="cluster-recycles-" + jid),
                              "sampling_steps": st.number_input("Sampling steps", 1, 1000, 200, key="cluster-steps-" + jid),
                              "diffusion_samples": st.number_input("Samples", 1, 20, 1, key="cluster-samples-" + jid)}
            if st.button("Queue on " + config.display_name(row.get("node_name") or "owning workstation"), type="primary", key="cluster-next-" + jid):
                result = submit("/api/v1/workflows", {"schema_version": 1, "workflow": workflow, "parameters": parameters,
                                                     "project_id": row["project_id"], "source_job": "cluster:" + jid,
                                                     "execution": {"placement": "auto", "idempotency_key": str(uuid.uuid4())}})
                st.success("Queued with data locality preserved · " + result["job"].get("job_code", "submitted"))


def _render_remote_pipeline_results(jid, row, token):
    from mn_molecular_design_studio.screening_results_ui import render_summary
    st.caption("Results are stored on " + config.display_name(row.get("node_name") or "the organizer"))
    key = "cluster-screening-inspect-" + jid
    if st.button("Refresh combined results", key=key + "-refresh"):
        st.session_state.pop(key, None)
        st.session_state.pop(key + "-response", None)
    cached = st.session_state.get(key + "-response")
    if cached is not None:
        render_summary(cached.get("screening_results"), (cached.get("result") or {}).get("metrics") or row.get("metrics") or {},
                       name=row.get("campaign_name") or "Pipeline", status=row.get("status") or "unknown")
        return
    if not st.session_state.get(key):
        try:
            response = _request("/jobs/" + jid + "/actions", token, {"action": "inspect"}, method="POST")
            st.session_state[key] = response["command_id"]
        except (ValueError, transport.ClusterUnavailable) as exc:
            st.warning(str(exc))
            return
    @st.fragment(run_every=2)
    def completed_inspection():
        try:
            operation = _request("/commands/" + st.session_state[key], token)
        except (ValueError, transport.ClusterUnavailable) as exc:
            st.info(str(exc))
            return
        if operation.get("status") == "failed":
            st.error(operation.get("error") or "Could not retrieve the campaign results.")
            return
        if operation.get("status") != "completed":
            st.info("Loading the combined campaign results from its organizer…")
            return
        results = (operation.get("response") or {}).get("results") or {}
        st.session_state[key + "-response"] = results
        st.rerun(scope="app")
    completed_inspection()
