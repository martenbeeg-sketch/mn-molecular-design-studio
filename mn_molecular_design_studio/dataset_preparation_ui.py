"""Coordinate preparation from the exact structure selected in Datasets."""
from __future__ import annotations

import json
from pathlib import Path

import streamlit as st

from mn_molecular_design_studio import datasets, molecule_imports
from mn_molecular_design_studio.database_actions_ui import ERRORS, read_only_workspace
from mn_molecular_design_studio.molecular_database import MolecularDatabase
from mn_molecular_design_studio.molecular_preparation import prepare_coordinate_subset
from mn_molecular_design_studio.navigation import open_datasets, open_workflow_tool, request_dataset_import
from mn_molecular_design_studio.preparation_history import extend_history, report_steps, sha256


def render_preparation(view, user, record, *, projects=(), preferred_projects=(), submit=None):
    formats = {"pdb", "ent", "cif", "mmcif"}
    structures = {item["id"]: item for item in record["structures"] if item["format"] in formats}
    if not structures or read_only_workspace():
        return
    if not any(c["type"] in {"protein", "dna", "rna"} for c in record["components"]):
        return
    panel = st.expander("Prepare structure", key=f"dataset-prepare-{record['id']}", on_change="rerun")
    if not panel.open:
        return
    with panel:
        selected = st.selectbox("Source structure", list(structures), format_func=lambda key: structures[key]["label"],
                                key=f"dataset-prepare-source-{view.group_id}-{record['id']}")
        key = "dataset-prepare-" + selected
        try:
            path = view.structure_file(user, selected)
            if path.suffix.lower() in {".pdb", ".ent"} and st.button("Open repair and extension tools", key=key + "-repair"):
                st.session_state["ligand-complex-prepare-pending-source"] = str(path)
                open_workflow_tool(st.session_state, "ligand.complex.prepare")
                st.rerun()
            data = path.read_bytes()
            if submit and path.suffix.lower() in {".pdb", ".ent"} and any(
                    component["type"] == "protein" for component in record["components"]):
                from mn_molecular_design_studio.dataset_modelling_ui import render_modelling

                render_modelling(path, data.decode("utf-8"), key, projects, preferred_projects, submit)
            model = next(molecule_imports.read_coordinate_structure(data, path.name).get_models())
            ranges = {chain.id.strip() or "_": (min(res.id[1] for res in chain), max(res.id[1] for res in chain))
                      for chain in model if len(chain)}
            st.caption("Trim the first coordinate model and save the selection as a prepared variant.")
            keep = st.multiselect("Keep chains", list(ranges), default=list(ranges), key=key + "-chains")
            chosen = {}
            for chain in keep:
                low, high = ranges[chain]
                chosen[chain] = tuple(st.slider(f"Chain {chain} residue range", min_value=low, max_value=high,
                    value=(low, high), key=key + "-" + chain)) if low < high else (low, high)
            remove_waters = st.checkbox("Remove waters", value=True, key=key + "-waters")
            remove_nonpolymer = st.checkbox("Remove ligands and other non-polymer residues", key=key + "-ligands")
            if st.button("Save prepared variant", disabled=not chosen, key=key + "-save"):
                changed = prepare_coordinate_subset(data, path.name, chosen, remove_waters=remove_waters,
                                                    remove_nonpolymer=remove_nonpolymer)
                filename = Path(record["name"]).stem + "-prepared" + path.suffix
                prepared, = molecule_imports.parse_structure(changed, filename)
                settings = {"ranges": chosen, "remove_waters": remove_waters, "remove_nonpolymer": remove_nonpolymer}
                history = extend_history((structures[selected].get("context") or {}).get("preparation_history"),
                                         sha256(data), sha256(changed), [{"label": "Trim structure", "details": settings}])
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
                st.session_state["molecular-database-flash"] = "Prepared structure saved."
                open_datasets(st.session_state, group_id=view.group_id, molecule_id=molecule_id)
                st.rerun()
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
                parent = _source_history(user, input_sha) if user and input_sha else None
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
