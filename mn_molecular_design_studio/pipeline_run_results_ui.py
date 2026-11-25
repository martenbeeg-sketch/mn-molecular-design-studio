"""Step navigation and collected outputs for an authorized pipeline run."""

from graphlib import CycleError, TopologicalSorter
import json
from pathlib import Path

import pandas as pd
import streamlit as st

from mn_molecular_design_studio import logical_runs
from mn_molecular_design_studio.cluster import config
from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path
from mn_molecular_design_studio.result_details_ui import render_result_details, result_field_label


def _step_order(graph, states):
    labels = {node["id"]: node["data"].get("label") or node["id"] for node in graph.get("nodes", [])}
    order = dict.fromkeys([*labels, *states])
    sorter = TopologicalSorter({node: [] for node in order})
    for edge in graph.get("edges", []):
        source, target = edge.get("source"), edge.get("target")
        if source in order and target in order:
            sorter.add(target, source)
    try:
        return list(sorter.static_order()), labels
    except CycleError:
        return list(order), labels


def _status(value):
    return str(value or "unknown").replace("_", " ").capitalize()


def _step_result_reference(state, parent_ref, run_dir, remote_node):
    output = state.get("output") or {}
    return _job_reference(
        state.get("job_run_dir") or (output.get("tile_run") or {}).get("run_dir")
        or (output.get("analysis_run") or {}).get("run_dir") or output.get("run_dir"),
        parent_ref, run_dir, remote_node,
    )


def _job_reference(value, parent_ref, run_dir, remote_node):
    if not value or run_dir is None:
        return None
    # Production receipts can contain a portable URI or a legacy path. For
    # remote controllers, use their recorded root without resolving it here.
    root = run_dir.parent.parent
    try:
        if parent_ref.startswith("cluster:"):
            if not remote_node:
                return None
            relative = Path(value).relative_to(root)
            if len(relative.parts) != 2 or ".." in relative.parts:
                return None
            reference = "cluster:" + config.job_id(remote_node, relative.as_posix())
            return reference if reference != parent_ref else None
        path = resolve_stored_path(value, run_dir=run_dir)
        if path is None:
            return None
        relative = path.resolve().relative_to(root.resolve())
        if len(relative.parts) == 2 and (path / "metadata.json").is_file():
            reference = relative.as_posix()
            return reference if reference != parent_ref else None
    except (OSError, ValueError):
        pass
    return None


def _select_step(click_key, state_key, order):
    click = st.session_state.get(click_key) or {}
    index = click.get("row")
    if isinstance(index, int) and 0 <= index < len(order):
        st.session_state[state_key] = order[index]
        st.session_state[state_key + "-open"] = order[index]


def open_step_job(reference, parent_ref, step_label):
    st.session_state["pipeline-step-return"] = {"parent": parent_ref, "child": reference, "step": step_label}
    st.session_state["pipeline-final-results-reference"] = None
    st.session_state["selected-job-reference"] = reference
    st.session_state["page"] = "Job details"


def _batch_rows(batches, parent_ref, run_dir, remote_node, user_id):
    workstations = {node["id"]: node.get("name") or node["id"]
                    for node in st.session_state.get(f"cluster-nodes-snapshot:{user_id}", [])
                    if isinstance(node, dict) and node.get("id")}
    rows, references = [], []
    for index, entry in enumerate(batches, 1):
        batch, worker = entry["receipt"], entry["worker"]
        attempt = (batch.get("attempts") or [{}])[-1]
        cluster_id = entry["cluster_job_id"]
        reference = ("cluster:" + str(cluster_id)) if cluster_id else _job_reference(
            batch.get("job_run_dir") or attempt.get("run_dir"), parent_ref, run_dir, remote_node,
        )
        reference = reference if reference != parent_ref else None
        references.append(reference)
        node = worker.get("node_id") or batch.get("worker") or attempt.get("worker")
        name = worker.get("node_name") or workstations.get(node) or node
        rows.append({
            "Batch": index, "Results": "Open" if reference else "", "Status": _status(entry["status"]),
            "Outputs": batch.get("candidate_count"), "Inputs": batch.get("unit_count"),
            "Workstation": config.display_name(name) if name else "—",
            "Attempts": len(batch.get("attempts") or []),
            "Detail": batch.get("error") or attempt.get("error") or batch.get("wait_reason") or "",
        })
    return rows, references


def _open_batch(click_key, references, parent_ref, step_label):
    click = st.session_state.get(click_key) or {}
    index = click.get("row")
    if isinstance(index, int) and 0 <= index < len(references) and references[index]:
        open_step_job(references[index], parent_ref, step_label)


def _collected_preview(value, run_dir, key):
    if not value or run_dir is None:
        return False
    try:
        path = resolve_stored_path(value, run_dir=run_dir)
        if path is None or not path.resolve().is_relative_to(run_dir.resolve()) or not path.is_file():
            return False
    except (OSError, ValueError):
        return False
    panel = st.expander("Collected candidate preview", key=f"{key}-collected", on_change="rerun")
    if not panel.open:
        return True
    with panel:
        candidates = []
        try:
            with path.open(encoding="utf-8") as stream:
                for line in stream:
                    if not line.strip():
                        continue
                    row = json.loads(line)
                    if isinstance(row, dict):
                        candidates.append(row)
                    if len(candidates) > 100:
                        break
        except (OSError, ValueError):
            st.caption("The collected candidate file could not be read.")
            return True
        if not candidates:
            st.caption("The collected file contains no candidate rows.")
            return True
        rows = []
        for candidate in candidates[:100]:
            metrics = candidate.get("metrics") or {}
            row = {"Candidate": candidate.get("candidate_id"), "Method": candidate.get("source_tool")}
            measurements = [(name, value) for name, value in metrics.items()
                            if isinstance(value, (float, int, bool))]
            row.update({result_field_label(name): value for name, value in measurements[:12]})
            rows.append(row)
        st.dataframe(rows, hide_index=True, width="stretch")
        if len(candidates) > 100:
            st.caption("Showing the first 100 collected candidates.")
    return True


def render_pipeline_steps(parent_ref, parent, graph, execution, linked, *, run_dir=None, user_id=None):
    states = execution.get("nodes") or {}
    order, labels = _step_order(graph, states)
    if not order:
        st.caption("Step execution has not been recorded for this run.")
        return
    summaries = {step["id"]: step for step in logical_runs.pipeline_step_summary(execution, parent, linked)}
    step_rows = []
    for node in order:
        summary = summaries.get(node) or {}
        batches = summary.get("worker_batch_statuses") or {}
        step_rows.append({
            "Step": labels.get(node, summary.get("label") or node),
            "Status": _status(summary.get("status") or "pending"),
            "Outputs": summary.get("collected_candidate_count"),
            "Batches": " · ".join(f"{count} {status.replace('_', ' ')}" for status, count in batches.items()),
        })
    errors = [row["Step"] for row in step_rows if row["Status"] == "Failed"]
    if errors:
        st.warning(f"{len(errors)} {'step failed' if len(errors) == 1 else 'steps failed'}. Select a failed step for its recorded error.")
    key = f"pipeline-run-step-{parent_ref}"
    click_key = key + "-click"
    st.dataframe(pd.DataFrame(step_rows), hide_index=True, width="stretch",
                 column_config={"Step": st.column_config.ButtonColumn(
                     "Step", type="tertiary", alignment="left", key=click_key,
                     on_click=_select_step, args=(click_key, key, order),
                 )})
    st.caption("Select a step to open its tool results or view its batches and execution details.")

    selected = st.session_state.get(key)
    request = st.session_state.pop(key + "-open", None)
    if selected not in order:
        return
    step = step_rows[order.index(selected)]
    state = states.get(selected) or {}
    summary = summaries.get(selected) or {}
    batches = summary.get("batches") or []
    remote_node = parent.get("node_id") if parent_ref.startswith("cluster:") else None
    direct = None if batches else _step_result_reference(state, parent_ref, run_dir, remote_node)
    if not direct and not batches:
        recorded = [job for job in linked if job.get("graph_node_id") == selected
                    and logical_runs.reference(job) != logical_runs.reference(parent)]
        if len(recorded) == 1:
            direct = logical_runs.source_reference(recorded[0])
    if direct and request == selected:
        open_step_job(direct, parent_ref, step["Step"])
        st.rerun(scope="app")

    st.subheader(step["Step"])
    caption = step["Status"]
    if step["Outputs"] is not None:
        caption += f" · {step['Outputs']:,} collected candidate outputs"
    st.caption(caption)
    if state.get("error"):
        st.error(str(state["error"]))
    if direct:
        st.button("Open tool results", key=key + "-tool-" + selected, on_click=open_step_job,
                  args=(direct, parent_ref, step["Step"]))
    if batches:
        rows, references = _batch_rows(batches, parent_ref, run_dir, remote_node, user_id)
        batch_click = key + "-batch-" + selected
        st.dataframe(rows, hide_index=True, width="stretch",
                     column_config={"Results": st.column_config.ButtonColumn(
                         "Results", type="tertiary", key=batch_click, pinned=True,
                         on_click=_open_batch, args=(batch_click, references, parent_ref, step["Step"]),
                     )})
        if not any(references):
            st.caption("No worker job is linked to the recorded batches.")
    preview = not parent_ref.startswith("cluster:") and _collected_preview(
        (state.get("output") or {}).get("candidates_jsonl"), run_dir, key + "-" + selected,
    )
    if not (preview or batches or direct or state.get("error")):
        st.caption("This step has no separate tool result.")
        graph_node = next((node for node in graph.get("nodes", []) if node.get("id") == selected), {})
        params = graph_node.get("data", {}).get("params")
        if params:
            with st.expander("Recorded step settings"):
                render_result_details(st, params)
