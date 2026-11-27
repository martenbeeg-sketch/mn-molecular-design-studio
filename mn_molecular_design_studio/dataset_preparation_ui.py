"""Coordinate preparation from the exact structure selected in Datasets."""
from __future__ import annotations

import json
from pathlib import Path

import streamlit as st

from mn_molecular_design_studio import datasets, molecule_imports
from mn_molecular_design_studio.database_actions_ui import ERRORS, read_only_workspace
from mn_molecular_design_studio.molecular_database import MolecularDatabase
from mn_molecular_design_studio.molecular_preparation import coordinate_chain_ranges, prepare_coordinate_subset
from mn_molecular_design_studio.navigation import open_datasets, request_dataset_import
from mn_molecular_design_studio.preparation_history import extend_history, report_steps, sha256


COORDINATE_FORMATS = {"pdb", "ent", "cif", "mmcif"}


def render_coordinate_preview(data, *, key, height=350, show_ligands=False):
    if len(data) > 50_000_000:
        st.caption("These coordinates are too large for the inline viewer.")
        return
    from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer import molstar_custom_component
    from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer.dataclasses import StructureVisualization

    molstar_custom_component([StructureVisualization(pdb=data.decode("utf-8", errors="replace"), color="chain-id",
        representation_type="cartoon+ball-and-stick" if show_ligands else "cartoon")],
        key=key, height=height, show_structure_tools=False)


def _dismiss_preparation():
    st.session_state.pop("datasets-prepare-request", None)


@st.dialog("Prepare structure", width="large", on_dismiss=_dismiss_preparation)
def render_preparation_dialog(view, user, record, *, preferred_structure="", projects=(), preferred_projects=(), submit=None):
    st.caption(record["name"])
    render_preparation(view, user, record, projects=projects, preferred_projects=preferred_projects,
                       submit=submit, preferred_structure=preferred_structure)


def render_preparation(view, user, record, *, projects=(), preferred_projects=(), submit=None,
                       preferred_structure=""):
    structures = {item["id"]: item for item in record["structures"] if item["format"] in COORDINATE_FORMATS}
    if not structures or read_only_workspace():
        return
    if not any(c["type"] in {"protein", "dna", "rna"} for c in record["components"]):
        return
    preview, options = st.columns([5, 4])
    with options:
        _render_preparation_options(view, user, record, structures, projects, preferred_projects, submit,
                                    preferred_structure=preferred_structure, preview=preview)


def _render_preparation_options(view, user, record, structures, projects, preferred_projects, submit, *,
                                preferred_structure, preview):
    structure_ids = list(structures)
    selected = st.selectbox("Source structure", structure_ids, index=structure_ids.index(preferred_structure)
                            if preferred_structure in structures else 0,
                            format_func=lambda identifier: structures[identifier]["label"],
                            key=f"dataset-prepare-source-{view.group_id}-{record['id']}-{preferred_structure}")
    key = "dataset-prepare-" + selected
    try:
        path = view.structure_file(user, selected)
        pdb = path.suffix.lower() in {".pdb", ".ent"}
        data = path.read_bytes()
        ranges = coordinate_chain_ranges(data, path.name)
        keep = st.multiselect("Keep chains", list(ranges), default=list(ranges), key=key + "-chains",
                              format_func=lambda chain: chain if ranges[chain] is not None else f"{chain} · non-polymer only")
        chosen = {}
        for chain in keep:
            bounds = ranges[chain]
            if bounds and bounds[0] < bounds[1]:
                bounds = tuple(st.slider(f"Chain {chain} residue range", min_value=bounds[0], max_value=bounds[1],
                                         value=bounds, key=key + "-" + chain))
            chosen[chain] = bounds
        remove_waters = st.checkbox("Remove waters", value=True, key=key + "-waters")
        remove_nonpolymer = st.checkbox("Remove ligands and other non-polymer residues", key=key + "-ligands")
        has_polymer = any(bounds is not None for bounds in chosen.values())
        save = st.button("Save prepared variant", disabled=not has_polymer, key=key + "-save")
        changed = prepare_coordinate_subset(data, path.name, chosen, remove_waters=remove_waters,
                                            remove_nonpolymer=remove_nonpolymer) if has_polymer else b""
        run_tools = submit and pdb and has_polymer and any(
            component["type"] == "protein" and component["chain"] in chosen for component in record["components"])
        if save or run_tools:
            settings = {"ranges": chosen, "remove_waters": remove_waters, "remove_nonpolymer": remove_nonpolymer}
            history = extend_history((structures[selected].get("context") or {}).get("preparation_history"),
                sha256(data), sha256(changed), [{"label": "Trim structure", "details": settings}])
        with preview:
            if has_polymer:
                st.caption("Selected chains and residue ranges")
                render_coordinate_preview(changed, key=key + "-viewer", height=410,
                    show_ligands=not remove_nonpolymer and any(component["type"] == "small_molecule" for component in record["components"]))
            else:
                st.info("Keep at least one polymer chain to preview and save a prepared variant.")
        if save:
            filename = Path(record["name"]).stem + "-prepared" + path.suffix
            prepared, = molecule_imports.parse_structure(changed, filename)
            prepared["files"][0].update(label="Prepared structure", derived_from_structure_id=selected,
                                        preparation=settings, preparation_history=history)
            writable = MolecularDatabase(view.root)
            writable.group(view.group_id, user, write=True)
            if molecule_imports.canonical_components(record["components"])[1] == molecule_imports.canonical_components(prepared["components"])[1]:
                writable.attach_structure(view.group_id, record["id"], prepared["files"][0])
                molecule_id = record["id"]
            else:
                prepared["metadata"]["derived_from_structure_id"] = selected
                result = writable.import_records(user, datasets.input_records([prepared]), group_id=view.group_id)
                molecule_id = result["molecule_ids"][0]
            _dismiss_preparation()
            st.session_state["molecular-database-flash"] = "Prepared structure saved."
            open_datasets(st.session_state, group_id=view.group_id, molecule_id=molecule_id)
            st.rerun()
        if run_tools:
            from mn_molecular_design_studio.dataset_modelling_ui import render_protein_tools

            render_protein_tools(user, path, changed, history, key, projects, preferred_projects, submit)
    except (*ERRORS, UnicodeError, StopIteration) as exc:
        st.error(f"Could not prepare the structure: {exc}")


def render_prepared_outputs(artifact_paths, run_id, project_id, *, user=None):
    names = {"prepared_target.pdb": "Prepared protein", "prepared_complex.pdb": "Prepared complex"}
    prepared = {path.name: path for _, path in artifact_paths if path.name in names}
    if not prepared or read_only_workspace():
        return
    selected = st.selectbox("Prepared output", list(prepared), format_func=names.get, key=f"prepared-output-{run_id}")
    if st.button("Save output in Datasets", key=f"prepared-output-save-{run_id}"):
        path = prepared[selected]
        try:
            with path.open("rb") as handle:
                content = handle.read(molecule_imports.MAX_BYTES + 1)
            if len(content) > molecule_imports.MAX_BYTES:
                raise ValueError("The prepared structure exceeds the 100 MB import limit.")
            report = next((item for _, item in artifact_paths
                           if item.name == "repair_report.json" and item.parent == path.parent), None)
            history = None
            if report is not None:
                raw = report.read_bytes()
                payload = json.loads(raw)
                output_sha = (payload.get("output_sha256") or {}).get(path.name)
                if output_sha and output_sha != sha256(content):
                    raise ValueError("The prepared coordinates no longer match the run's repair report.")
                input_sha = payload.get("source_sha256")
                parent = payload.get("input_preparation_history")
                if parent is None and user and input_sha:
                    parent = _source_history(user, input_sha)
                history = extend_history(parent, input_sha, sha256(content), report_steps(payload), artifact={
                    "run_id": run_id, "path": path.name, "sha256": sha256(content), "report_sha256": sha256(raw),
                })
            request_dataset_import(st.session_state, "prepared", project_id=project_id,
                                   return_to=st.session_state.get("page", ""))
            st.session_state["molecular-import-prepared"] = {
                "content": content, "name": f"{path.stem}-{run_id}.pdb", "run_id": run_id,
                "preparation_history": history,
            }
            st.rerun()
        except (OSError, ValueError) as exc:
            st.error(f"Could not open the prepared output: {exc}")

def _source_history(user, source_sha256):
    from mn_molecular_design_studio.data_browser import ReadOnlyLibrary

    database = ReadOnlyLibrary()
    allowed = [group["id"] for group in database.groups(user, collection="library", scope="all")]
    with database.connection() as connection:
        results = connection.execute("""
            SELECT DISTINCT json_extract(context, '$.preparation_history') AS history
            FROM molecule_structures
            WHERE group_id IN (SELECT value FROM json_each(?))
              AND coalesce(json_extract(context, '$.portable_sha256'), json_extract(context, '$.content_sha256'))=?
              AND json_type(context, '$.preparation_history')='object'
        """, (json.dumps(allowed), source_sha256)).fetchall()
    return json.loads(results[0]["history"]) if len(results) == 1 else None
