"""Display structural clusters and navigate to the exact saved complexes."""
from __future__ import annotations

import csv
import io
import json
from pathlib import Path

import streamlit as st

from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json_raw
from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path


def _collection_file(run_dir, node):
    reference = (node.get("output") or {}).get("candidates_jsonl")
    if not isinstance(reference, str) or reference.startswith("cluster:"):
        return None
    path = resolve_stored_path(reference, run_dir=run_dir, must_exist=True)
    if path is not None and path.is_file() and path.resolve().is_relative_to(run_dir.resolve()):
        return path


@st.cache_data(max_entries=16, show_spinner=False)
def _cluster_rows(assignments, selected, revision):
    """Keep only table values and exact Database references from the retained receipts."""
    representatives = set()
    if selected is not None:
        with selected.open(encoding="utf-8") as stream:
            for line in stream:
                if line.strip():
                    representatives.add(json.loads(line)["candidate_id"])
    rows = []
    with assignments.open(encoding="utf-8") as stream:
        for line in stream:
            if not line.strip():
                continue
            candidate = json.loads(line)
            metrics = candidate.get("metrics") or {}
            identifier = candidate["candidate_id"]
            raw = candidate.get("raw_metadata") or {}
            provenance = raw.get("metadata") or {}
            rows.append({
                "candidate_id": identifier,
                "name": raw.get("name") or identifier,
                "cluster": metrics.get("diversity_cluster"),
                "size": metrics.get("diversity_cluster_size"),
                "selected": identifier in representatives if selected is not None else metrics.get("diversity_representative") is True,
                "status": metrics.get("qc_status") or metrics.get("analysis_status") or "",
                "molecule_id": raw.get("studio_molecule_id") or "",
                "group_id": provenance.get("source_group_id") or "",
                "structure_id": provenance.get("structure_id") or "",
            })
    rows.sort(key=lambda row: (row["cluster"] is None, row["cluster"] or 0, row["name"], row["candidate_id"]))
    return rows


def inspect_clustering_pipeline(job):
    """Use collected analysis artifacts for a pipeline dedicated to Foldseek."""
    if job.get("tool") != "visual_pipeline":
        return None
    run_dir = Path(job["run_dir"])
    graph = read_json_raw(run_dir / "graph_snapshot.json")
    nodes = {row["id"]: row["data"] for row in graph.get("nodes", [])}
    if set(row["kind"] for row in nodes.values()) - {"foldseek", "molecule_source", "filter", "results"}:
        return None
    foldseek = [node_id for node_id, row in nodes.items() if row["kind"] == "foldseek"]
    if len(foldseek) != 1:
        return None
    node_id = foldseek[0]
    execution = read_json_raw(run_dir / "graph_execution.json")
    states = execution.get("nodes") or {}
    assignments = _collection_file(run_dir, states.get(node_id) or {})
    if assignments is None:
        return None
    quotas = [edge["target"] for edge in graph.get("edges", [])
              if edge["source"] == node_id and nodes.get(edge["target"], {}).get("kind") == "filter"
              and nodes[edge["target"]].get("params", {}).get("diversity_field", "").split("/")[-1] == "diversity_cluster"]
    if len(quotas) > 1 or any(row["kind"] == "filter" and key not in quotas for key, row in nodes.items()):
        return None
    shortlist = _collection_file(run_dir, states.get(quotas[0]) or {}) if quotas else None
    if quotas and shortlist is None:
        return None
    try:
        files = [assignments, *([shortlist] if shortlist is not None else [])]
        rows = _cluster_rows(assignments, shortlist, tuple((p.stat().st_mtime_ns, p.stat().st_size) for p in files))
    except (OSError, ValueError, KeyError, TypeError):
        return None
    return {"rows": rows, "graph": graph, "execution": execution,
            "settings": nodes[node_id].get("params") or {}}


def saved_members(rows):
    return [
        {field: row[field] for field in ("molecule_id", "group_id", "structure_id")}
        for row in rows if all(row[field] for field in ("molecule_id", "group_id", "structure_id"))
    ]


def _open_members(project_id, members, *, focus=False):
    destination = {"project_id": project_id, "foldseek_members": members}
    if focus:
        destination["molecule_id"] = members[0]["molecule_id"]
    st.session_state.update({"page": "Data Browser", "data-browser-open": destination})


def _open_row(key, project_id, rows):
    index = (st.session_state.get(key) or {}).get("row")
    if isinstance(index, int) and 0 <= index < len(rows):
        members = saved_members([rows[index]])
        if members:
            _open_members(project_id, members, focus=True)


def render_clustering_results(summary, project_id, run_id):
    rows = summary["rows"]
    representatives = [row for row in rows if row["selected"]]
    cluster_count = len({row["cluster"] for row in rows if row["cluster"] is not None})
    st.subheader("Foldseek structural clusters")
    st.caption(f"{len(rows):,} designs · {cluster_count:,} clusters · {len(representatives):,} selected")
    settings = summary["settings"]
    st.caption(f"TM-score threshold {settings.get('similarity_threshold', 0.5):g} · Alignment coverage {settings.get('coverage', 0.8):.0%}")
    scope = st.segmented_control("Show", ["Representatives", "All designs"], key=f"foldseek-scope-{run_id}",
                                 default="Representatives" if representatives else "All designs")
    visible = rows if scope == "All designs" else representatives
    members = saved_members(representatives)
    controls = st.columns(2)
    if members and project_id:
        controls[0].button(f"Open {len(members):,} representatives in Data browser", type="primary",
                           key=f"foldseek-cohort-{run_id}", on_click=_open_members, args=(project_id, members))
    table = [{"Design": row["name"], "Cluster": row["cluster"], "Cluster size": row["size"],
              "Selected": row["selected"], "Structure": "Open" if project_id and saved_members([row]) else None}
             for row in visible]
    if not table:
        st.info("No representatives were selected in this run.")
        return
    key = f"foldseek-open-row-{run_id}-{scope}"
    st.dataframe(
        table, hide_index=True, width="stretch",
        column_config={"Structure": st.column_config.ButtonColumn(
            "Saved structure", type="tertiary", key=key, on_click=_open_row, args=(key, project_id, visible))},
    )
    buffer = io.StringIO()
    export_fields = ["Design", "Cluster", "Cluster size", "Selected", "Candidate ID"]
    writer = csv.DictWriter(buffer, fieldnames=export_fields)
    writer.writeheader()
    for item, row in zip(table, visible):
        writer.writerow({field: item[field] for field in export_fields[:-1]} | {"Candidate ID": row["candidate_id"]})
    controls[1].download_button("Download table CSV", buffer.getvalue(), file_name="foldseek-clusters.csv",
                                mime="text/csv", key=f"foldseek-download-{run_id}", on_click="ignore")
