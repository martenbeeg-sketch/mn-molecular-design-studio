"""Frozen Database selections passed to the existing tool input pages."""
from __future__ import annotations
import hashlib
import uuid


def freeze_protein_selection(database, user, ids, project_id, *, name="Database selection"):
    from pathlib import Path
    import shutil
    from mn_molecular_design_studio.auth import require_project_access
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import create_job, finish_job, write_json
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import write_candidates
    from mn_molecular_design_studio.workspace import mark_studio_job
    require_project_access(user, project_id)
    records = database.list_molecules(user, molecule_ids=ids, limit=100_000)["items"]
    by_id = {row["id"]: row for row in records}
    if set(ids) != set(by_id):
        raise PermissionError("Some selected molecules are inaccessible.")
    if any(not any(component["type"] == "protein" for component in row["components"]) for row in records):
        raise ValueError("A protein pipeline selection must contain protein molecules or complexes.")
    job = create_job("molecule-selections", "molecule_selection", "database_selection", {},
                     {"job_name": name, "molecule_ids": ids, "selection_type": "protein_candidates/v1"})
    candidates = []
    for identifier in ids:
        record = by_id[identifier]
        proteins = [row for row in record["components"] if row["type"] == "protein"]
        lead = next((row for row in proteins if row.get("role") == "lead"), proteins[0])
        structures = [row for row in record["structures"] if row["format"] in {"pdb", "ent", "cif", "mmcif"}]
        path = None
        if structures:
            source = database.structure_file(user, structures[0]["id"])
            path = job.run_dir / "artifacts/structures" / (identifier + source.suffix)
            path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, path)
        candidates.append({"candidate_id": identifier, "source_tool": "database_selection", "stage": "generation.backbone_sequence",
                           "binder_sequence": lead["sequence"], "binder_chains": [lead["chain"]],
                           "target_chains": [row["chain"] for row in record["components"] if row["chain"] != lead["chain"]],
                           "complex_pdb": str(path) if path else None, "binder_pdb": str(path) if path and len(record["components"]) == 1 else None,
                           "metrics": {"library_selection": True}, "raw_metadata": {
                               "studio_molecule_id": identifier, "name": record["name"],
                               "metadata": record.get("metadata") or {}, "components": record["components"],
                           },
                           "parents": [{"molecule_id": identifier}]})
    write_json(job.run_dir / "selection_snapshot.json", {"schema_version": 1, "owner_id": user["id"], "project_id": project_id,
               "molecule_ids": ids, "molecules": [{"id": row["id"], "name": row["name"],
                                                     "metadata": row.get("metadata") or {},
                                                     "components": row["components"]} for row in records]})
    write_candidates(job.run_dir, "database_selection", candidates)
    finish_job(job.run_dir, True, {"outputs": {"selection": "selection_snapshot.json"}, "metrics": {"molecule_count": len(ids)}})
    mark_studio_job({"run_dir": str(job.run_dir)}, project_id=project_id, owner_id=user["id"])
    return job.run_dir


def open_in_pipeline(database, user, ids, project_id):
    import streamlit as st
    from mn_molecular_design_studio import pipeline_graph
    source = freeze_protein_selection(database, user, ids, project_id)
    graph = {"schema_version": pipeline_graph.SCHEMA, "id": uuid.uuid4().hex, "name": "Database selection", "nodes": [
        {"id": "database-input", "type": "studio", "position": {"x": 100, "y": 100}, "data": {"kind": "molecule_source", "label": "Database selection", "params": {"source_run_dir": str(source)}}},
        {"id": "database-results", "type": "studio", "position": {"x": 600, "y": 100}, "data": {"kind": "results", "label": "Results", "params": {}}}],
        "edges": [{"id": "database-edge", "source": "database-input", "sourceHandle": "out", "target": "database-results", "targetHandle": "in"}],
        "viewport": {"x": 0, "y": 0, "zoom": 1}}
    prefix = f"visual-pipeline-{project_id}"
    st.session_state[prefix + "-graph"] = pipeline_graph.clean_graph(graph)
    st.session_state[prefix + "-reset"] = uuid.uuid4().hex
    st.session_state[prefix + "-name"] = "Database selection"
    st.session_state.setdefault("project-pending-selection", {})["visual-pipeline-project-id"] = project_id
    st.session_state["pipeline-workspace-view"] = "Node canvas"
    st.session_state["page"] = "Pipelines"


def editor_model(record):
    from mn_molecular_design_studio.cofolding_input import default_model
    model = default_model()
    model["compound_name"] = str(record.get("name") or "")
    model["compound_id"] = str((record.get("metadata") or {}).get("compound_id") or record.get("name") or "")
    model["cards"] = []
    for index, component in enumerate(record["components"], 1):
        ligand = component["type"] == "small_molecule"
        model["cards"].append({"id": index, "kind": "ligand" if ligand else component["type"],
                               "mode": "smiles" if component.get("smiles") else "ccd",
                               "value": component.get("sequence") or component.get("smiles") or component.get("ccd"),
                               "copies": 1, "modifications": [], "editing": False})
    return model


def set_editor_selection(prefix, record):
    import streamlit as st
    for key in list(st.session_state):
        if key == prefix + "-editor" or key.startswith(prefix + "-editor-library-"):
            st.session_state.pop(key, None)
    st.session_state[prefix + "-editor-model"] = editor_model(record)
    st.session_state[prefix + "-library-epoch"] = st.session_state.get(prefix + "-library-epoch", 0) + 1


def open_in_tool(database, user, ids, workflow, project_id):
    import streamlit as st
    from mn_molecular_design_studio.auth import require_project_access
    require_project_access(user, project_id)
    rows = database.molecules_by_ids(user, ids)
    if workflow == "pipeline":
        open_in_pipeline(database, user, ids, project_id)
        return
    if workflow == "boltz2":
        if all(len(row["components"]) == 1 and row["components"][0].get("smiles") for row in rows):
            key = f"design-validation-boltz2-{project_id}"
            st.session_state[key + "-kind"] = "Compound collection"
            st.session_state[key + "-collection-selected"] = database.compound_snapshot(user, ids, name="Database selection", project_id=project_id)
        elif len(rows) == 1:
            set_editor_selection("optimized-boltz2-1", rows[0])
        else:
            key = f"design-validation-boltz2-{project_id}"
            st.session_state[key + "-kind"] = "Database molecules"
            st.session_state[key + "-collection-molecules"] = rows
        st.session_state["selected_workflow"] = "refolding.optimized_prediction"
        st.session_state["selected_workflow_engine"] = "boltz2"
        st.session_state.setdefault("project-pending-selection", {})["optimized-prediction-boltz2-project-id"] = project_id
    else:
        receptor = workflow.endswith(":receptor")
        workflow = workflow.removesuffix(":receptor")
        state = {}
        if receptor:
            record = database.list_molecules(user, molecule_ids=ids)["items"][0]
            structure = next((row for row in record["structures"] if row["format"] in {"pdb", "ent", "cif", "mmcif"}), None)
            if structure is None:
                raise ValueError("This molecule has no receptor structure usable by the selected tool.")
            state["vina-target-path"] = str(database.structure_file(user, structure["id"], materialize=True))
        else:
            selection = database.compound_snapshot(user, ids, name="Database selection", project_id=project_id)
            selection["group_id"] = "selection:" + hashlib.sha256(selection["path"].encode()).hexdigest()[:24]
            selection["source_path"] = selection["path"]
            st.session_state["database-pending-compound-selection"] = selection
            state.update({"vina-compound-input": selection, "vina-ligand-format": "Database compound set",
                          "ligand-tool-compound-set-ligand.docking.vina": selection["group_id"]})
            st.session_state["ligand-tool-compound-set-" + workflow] = selection["group_id"]
        if workflow.startswith("ligand.docking."):
            from mn_molecular_design_studio import docking_preparations
            record = docking_preparations.save(user["id"], project_id, workflow, str(uuid.uuid4()), state)
            picker = "docking-preparation-choice-" + hashlib.sha256(f"{user['id']}:{project_id}:{workflow}".encode()).hexdigest()[:16]
            st.session_state[picker] = record["id"]
        st.session_state["selected_workflow"] = workflow
        st.session_state.pop("selected_workflow_engine", None)
        st.session_state.setdefault("project-pending-selection", {})[f"ligand-tool-{workflow}-project-id"] = project_id
    st.session_state.pop("selected_workflow_group", None)
    st.session_state["page"] = "Workflow setup"
