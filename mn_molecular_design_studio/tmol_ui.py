"""Shared TMol settings, standalone upload, and result presentation."""
from __future__ import annotations

import hashlib
from pathlib import Path

import streamlit as st

from mn_molecular_design_studio.auth import user_input_root
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json
from mn_molecular_design_studio.engine.workflow_engine.workflows.tmol import FORMATS, MODE_LABELS, WORKFLOW, validate_inputs


def _store_uploads(uploads, owner_id, formats=FORMATS):
    paths = []
    for upload in uploads:
        name = Path(upload.name).name
        content = upload.getvalue()
        if Path(name).suffix.lower() not in formats or not content:
            raise ValueError(f"Choose a supported coordinate file for {name}.")
        digest = hashlib.sha256(content).hexdigest()
        path = user_input_root(owner_id) / "tmol" / digest / name
        path.parent.mkdir(parents=True, exist_ok=True)
        if not path.is_file():
            path.write_bytes(content)
        paths.append(str(path))
    return paths


def render_settings(settings, *, owner_id, key, connected=False, upload=True):
    values = dict(settings)
    modes = list(MODE_LABELS)
    selected = values.get("mode", "fast_relax")
    values["mode"] = st.selectbox("Operation", modes, index=modes.index(selected) if selected in modes else 2,
                                 format_func=MODE_LABELS.get, key=key + "-mode")
    if values["mode"] == "fast_relax":
        values["repeats"] = int(st.number_input("Relaxation repeats", min_value=1,
            value=int(values.get("repeats", 1)), key=key + "-repeats"))
    values["prepare_ligands"] = st.checkbox(
        "Prepare ligands and nonstandard residues in complexes", value=bool(values.get("prepare_ligands")),
        help="Use the supplied structural chemistry. SDF and MOL2 inputs are parameterized automatically. Unsupported components are reported per structure.",
        key=key + "-ligand-preparation",
    )
    if connected:
        for field in ("structure_paths", "source_run_dir", "candidates_jsonl"):
            values.pop(field, None)
        st.caption("Use the connected structures.")
    elif upload:
        uploads = st.file_uploader("Molecular structures", type=[suffix[1:] for suffix in FORMATS],
                                   accept_multiple_files=True, key=key + "-structures")
        if uploads:
            values["structure_paths"] = _store_uploads(uploads, owner_id)
            values.pop("source_run_dir", None)
            values.pop("candidates_jsonl", None)
        elif paths := values.get("structure_paths"):
            st.caption(", ".join(Path(path).name for path in paths))
    with st.expander("Ligand parameter files"):
        parameter_files = st.file_uploader("Prepared ligand parameters", type=["tmol"], accept_multiple_files=True,
            help="Use TMol parameter files that match the ligand residue and atom names in the input structures.",
            key=key + "-ligand-parameters")
        if parameter_files:
            values["ligand_params_files"] = _store_uploads(parameter_files, owner_id, {".tmol"})
        elif paths := values.get("ligand_params_files"):
            st.caption(", ".join(Path(path).name for path in paths))
            if st.button("Clear ligand parameter files", key=key + "-clear-params"):
                values.pop("ligand_params_files")
    return values


def render_pipeline_inputs(settings, *, connected, owner_id, key):
    values = dict(settings)
    try:
        values = render_settings(settings, connected=connected, owner_id=owner_id, key=key)
        validate_inputs(values, connected=connected)
    except (OSError, ValueError) as exc:
        st.info(str(exc))
        return values, False
    return values, True


def render_setup(schema, *, owner_id, project_selector, submit, source_job_options):
    if st.button("← Back to All tools", key="back-to-tmol-tools"):
        st.session_state["page"] = "All tools"
        st.session_state.pop("selected_workflow", None)
        st.rerun()
    st.title("TMol")
    st.write("Score a saved geometry, minimize coordinates, or run FastRelax on existing molecular structures.")
    project_id = project_selector("tmol")
    source_mode = st.radio("Input source", ["Upload structures", "Saved results"], horizontal=True)
    try:
        parameters = render_settings({}, owner_id=owner_id, key="tmol-setup", upload=source_mode == "Upload structures")
        if source_mode == "Saved results" and project_id:
            sources = source_job_options(project_id)
            labels = {item["reference"]: item["label"] for item in sources}
            if labels:
                parameters["source_job"] = st.selectbox(
                    "Structure collection", list(labels), format_func=labels.get, key="tmol-source-job",
                )
            else:
                st.info("This project has no saved structure collections.")
        ready = bool(project_id and (parameters.get("structure_paths") or parameters.get("source_job")))
    except (OSError, RuntimeError, ValueError) as exc:
        st.error(str(exc))
        return
    if st.button("Run TMol", type="primary", disabled=not ready, key="tmol-submit"):
        submit(schema, WORKFLOW, parameters, project_id, {}, "TMol")


def render_results(run_dir, key):
    output = Path(run_dir) / "artifacts/tmol"
    summary = read_json(output / "summary.json")
    if not summary:
        return
    st.subheader("TMol structure results")
    st.caption(f'{MODE_LABELS.get(summary.get("mode"), "TMol")} · TMol {summary.get("tmol_version") or "unavailable"}'
               f' · {summary.get("gpu") or "GPU unrecorded"}')
    st.dataframe([{
        "Structure": row["name"], "Status": row["status"], "Input energy": row.get("tmol_energy_before"),
        "Output energy": row.get("tmol_energy_after"), "Energy change": row.get("tmol_energy_change"),
        "Input waters": row.get("input_waters"), "Output waters": row.get("output_waters"),
        "Error": row.get("error", ""),
    } for row in summary.get("records", [])], hide_index=True, width="stretch")
    st.caption("Energies use the TMol beta2016 score function and are reported in TMol score units.")
    table = output / "measurements.csv"
    if table.is_file():
        st.download_button("Download measurements CSV", table.read_bytes(), file_name="tmol_measurements.csv",
                           mime="text/csv", key=key + "-tmol-csv")
    completed = [row for row in summary.get("records", []) if row.get("status") == "completed"]
    if completed:
        choices = {row["id"]: row for row in completed}
        selected = st.selectbox("Coordinate output", list(choices), format_func=lambda value: choices[value]["name"],
                                key=key + "-tmol-structure")
        path = output / (selected + ".cif")
        if path.is_file():
            st.download_button("Download selected mmCIF", path.read_bytes(), file_name=path.name,
                               mime="chemical/x-mmcif", key=key + "-tmol-download")
