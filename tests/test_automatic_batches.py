"""Selections stay one scheduled job and use native container batches when possible."""
import inspect
import json
from pathlib import Path

import pytest

from mn_molecular_design_studio import batching, ligand_compounds
from mn_molecular_design_studio.engine.workflow_engine.core import jobs, workflow_queue
from mn_molecular_design_studio.engine.workflow_engine.services import workflow_automation as workflows
from mn_molecular_design_studio.engine.workflow_engine.workflows import ligand_tools, optimized_prediction


@pytest.fixture
def compound_file(tmp_path):
    def make(name="compounds", count=3):
        path = tmp_path / (name + ".json")
        ligand_compounds.write_compound_set(path, name=name, compounds=[
            {"compound_id": f"compound_{i}", "smiles": "CCO"} for i in range(count)
        ], owner_id="owner", project_id="project", created_at="2026-10-20T00:00:00Z")
        return path
    return make


def test_sequence_records_batch_independently_but_keep_complex_chains():
    records = batching.sequence_tasks(">two-chain complex\nACD:\nEFG\n>second\nACDEFG\n")
    assert [row["protein_sequences"] for row in records] == ["ACD:EFG", "ACDEFG"]
    assert batching.sequence_tasks("ACD:EFG\nACDEFG")[0]["protein_sequences"] == "ACD:EFG"


@pytest.mark.parametrize("engine", ["boltz2", "alphafold3", "af3_torch", "protenix_v2"])
def test_compound_collection_is_native_prediction_batch(compound_file, engine):
    prepared = optimized_prediction._prepare_prediction(
        engine=engine, protein_sequences="ACDEFGHIK", compound_set_file=compound_file(count=40),
        samples=1, accept_alphafold3_terms=True, resolve_random_seed=False,
    )
    assert len(prepared.normalized_tasks) == 40
    assert {task["name"] for task in prepared.normalized_tasks} == {f"compound_{i}" for i in range(40)}
    assert {task["compound_name"] for task in prepared.normalized_tasks} == {f"compound_{i}" for i in range(40)}
    assert {task["compound_id"] for task in prepared.normalized_tasks} == {f"compound_{i}" for i in range(40)}
    assert all(len(task["proteins"]) == 1 and len(task["ligand_smiles"]) == 1 for task in prepared.normalized_tasks)
    if engine == "boltz2":
        assert prepared.selected_options["use_potentials"] is True
        assert all(task["affinity_binder"] == "B" for task in prepared.normalized_tasks)
        assert [task["source_compound_id"] for task in prepared.normalized_tasks] == [f"compound_{i}" for i in range(40)]


def test_multiple_ligand_chains_in_one_complex_are_not_split():
    prepared = optimized_prediction._prepare_prediction(
        engine="boltz2", protein_sequences="ACDEFGHIK:KLMNPQ", ligand_smiles="CCO\nCCN",
        compound_name="Probe complex",
        samples=1, resolve_random_seed=False,
    )
    assert len(prepared.normalized_tasks) == 1
    assert len(prepared.normalized_tasks[0]["proteins"]) == 2
    assert len(prepared.normalized_tasks[0]["ligand_smiles"]) == 2
    assert prepared.normalized_tasks[0]["compound_id"] == "Probe complex"


def test_manual_ligand_prediction_requires_explicit_compound_name():
    with pytest.raises(ValueError, match="compound name"):
        optimized_prediction._prepare_prediction(
            engine="boltz2", protein_sequences="ACDEFGHIK", ligand_smiles="CCO",
            samples=1, resolve_random_seed=False,
        )


def test_default_compound_limit_keeps_whole_selection(compound_file):
    path = compound_file(count=137)
    for name in ("run_ligand_vina", "run_ligand_gnina", "run_ligand_autodock_gpu", "run_ligand_unidock_pro", "run_ligand_rosetta_ligand", "run_ligand_nesso"):
        default = inspect.signature(getattr(ligand_tools, name)).parameters["compound_limit"].default
        assert len(ligand_tools._validated_compound_set(path, default)[1]) == 137
    assert len(ligand_tools._validated_compound_set(path, 7)[1]) == 7


@pytest.fixture
def queue(monkeypatch, tmp_path):
    root = tmp_path / "runs"
    monkeypatch.setattr(jobs, "runs_root", lambda: root)
    monkeypatch.setattr("mn_molecular_design_studio.molecular_database_ingest.record_completed_run", lambda path: None)
    dispatched = []
    monkeypatch.setattr(workflow_queue, "spawn_worker_for_run", lambda path: dispatched.append(path))
    return root, dispatched


def test_scalar_adapter_is_one_job_continues_after_item_failure_and_keeps_resources(queue, monkeypatch, tmp_path):
    root, dispatched = queue
    calls = []
    def run_fake(value: int, cpu_cores: int = 2, gpu_device: str = "none", job_name: str = "test") -> Path:
        job = jobs.create_job("analysis", "analysis", "ipsae", {}, {"cpu_cores": cpu_cores, "gpu_device": gpu_device})
        calls.append((value, cpu_cores, gpu_device, jobs.read_json(job.run_dir / "metadata.json")))
        if value == 2:
            raise ValueError("bad selected input")
        jobs.finish_job(job.run_dir, True, {"metrics": {"answer": value * 2}})
        return job.run_dir
    monkeypatch.setattr(workflows, "_get_function", lambda spec: run_fake)
    monkeypatch.setattr(workflows, "_validate_workflow_inputs", lambda *args: None)
    request = tmp_path / "request.json"
    request.write_text(json.dumps({"schema_version": 1, "workflow": "analysis.ipsae", "parameters": [
        {"value": i, "cpu_cores": 2} for i in (1, 2, 3)
    ]}))
    response = workflows.submit_workflow_request(request)
    run = Path(response["job"]["run_dir"])
    assert dispatched == [run]
    metadata = jobs.read_json(run / "metadata.json")
    metadata.update(scheduler_backend="mn-compute-scheduler", resource_allocation={"cpu_cores": 2})
    jobs.write_json(run / "metadata.json", metadata)
    workflow_queue.run_workflow_call(run, jobs.read_json(run / "worker_request.json"))
    assert len(list(root.glob("*/*/metadata.json"))) == 1
    assert [call[:3] for call in calls] == [(i, 2, "none") for i in (1, 2, 3)]
    assert all(call[3]["resource_allocation"] == {"cpu_cores": 2} for call in calls)
    result = jobs.read_json(run / "result.json")
    assert [row["status"] for row in result["outputs"]["batch_items"]] == ["completed", "failed", "completed"]
    assert result["metrics"]["completed_items"] == 2
    assert result["success"] is True
    assert "1 of 3" in result["warning"]
    assert jobs.read_json(run / "metadata.json")["batch"]["completed"] == 2


def test_compatible_compound_selections_use_one_native_invocation(queue, compound_file, monkeypatch, tmp_path):
    root, dispatched = queue
    calls = []
    def run_fake(compound_set_file: Path, compound_limit: int = 0, cpu_cores: int = 2, gpu_device: str = "none") -> Path:
        rows = ligand_tools._validated_compound_set(compound_set_file, compound_limit)[1]
        calls.append(rows)
        job = jobs.create_job("ligand-docking", "ligand_docking", "vina", {"compound_set_file": str(compound_set_file)}, {})
        jobs.finish_job(job.run_dir, True, {"metrics": {"compound_count": len(rows)}})
        return job.run_dir
    monkeypatch.setattr(workflows, "_get_function", lambda spec: run_fake)
    monkeypatch.setattr(workflows, "_validate_workflow_inputs", lambda *args: None)
    request = tmp_path / "request.json"
    request.write_text(json.dumps({"schema_version": 1, "workflow": "ligand.docking.vina", "parameters": [
        {"compound_set_file": str(compound_file("first"))}, {"compound_set_file": str(compound_file("second"))}
    ]}))
    result = workflows.submit_workflow_request(request)
    run = Path(result["job"]["run_dir"])
    workflow_queue.run_workflow_call(run, jobs.read_json(run / "worker_request.json"))
    assert dispatched == [run] and len(calls) == 1 and len(calls[0]) == 6
    assert len({row["compound_id"] for row in calls[0]}) == 6
    assert len(list(root.glob("*/*/metadata.json"))) == 1


def test_pipeline_expands_compound_selection_in_one_native_request(compound_file):
    from mn_molecular_design_studio.engine.workflow_engine.workflows.graph_pipeline import _optimized_prediction_tasks
    source = {"molecule_task": {"task_type": "protein_ligand", "target_protein_sequences": "ACDEFGHIK",
                               "compound_set_file": str(compound_file())}}
    tasks = _optimized_prediction_tasks("boltz2", source, {"use_target_msa": True})
    assert len(tasks) == 3
    assert all(task["molecule_chain_options"] == [{"id": "A", "use_msa": True}] for task in tasks)


def test_scheduler_requirements_for_selection_are_maximum_not_sum(monkeypatch):
    from mn_molecular_design_studio.cluster.capabilities import requirements
    monkeypatch.setattr("mn_molecular_design_studio.cluster.capabilities.code_revision", lambda: "code")
    needs = requirements({"workflow": "refolding.optimized_prediction", "parameters": [
        {"engine": "boltz2", "cpu_cores": 4, "gpu_device": "0"},
        {"engine": "boltz2", "cpu_cores": 8, "gpu_device": "0"},
    ]}, installed_images={})
    assert needs["cpu_threads"] == 8
    assert needs["gpu_count"] == 1


def test_cancelled_batch_does_not_start_next_item_or_overwrite_cancellation(queue, monkeypatch, tmp_path):
    root, dispatched = queue
    calls = []
    def run_fake(value: int) -> Path:
        child = jobs.create_job("analysis", "analysis", "ipsae", {}, {})
        calls.append(value)
        jobs.update_status(dispatched[0], "cancelled")
        jobs.finish_job(child.run_dir, True, {})
        return child.run_dir
    monkeypatch.setattr(workflows, "_get_function", lambda spec: run_fake)
    monkeypatch.setattr(workflows, "_validate_workflow_inputs", lambda *args: None)
    request = tmp_path / "request.json"
    request.write_text(json.dumps({"schema_version": 1, "workflow": "analysis.ipsae",
                                   "parameters": [{"value": 1}, {"value": 2}]}))
    response = workflows.submit_workflow_request(request)
    run = Path(response["job"]["run_dir"])
    workflow_queue.run_workflow_call(run, jobs.read_json(run / "worker_request.json"))
    assert calls == [1]
    assert jobs.read_json(run / "metadata.json")["status"] == "cancelled"
    assert {row["status"] for row in jobs.read_json(run / "batch_manifest.json")["items"]} == {"cancelled"}


def test_sequence_file_is_one_native_prediction_collection(tmp_path):
    path = tmp_path / "sequences.fasta"
    path.write_text("".join(f">protein_{i}\nACDEFGHIK\n" for i in range(40)))
    prepared = optimized_prediction._prepare_prediction(engine="esmfold", sequence_file=path,
                                                         samples=1, resolve_random_seed=False)
    assert len(prepared.normalized_tasks) == 40
    assert prepared.normalized_tasks[-1]["name"] == "protein_39"


def test_pipeline_database_protein_selection_needs_no_structure(tmp_path):
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import write_candidates
    from mn_molecular_design_studio.engine.workflow_engine.workflows.graph_pipeline import _optimized_prediction_tasks
    run = tmp_path / "selection"
    run.mkdir()
    write_candidates(run, "database_selection", [
        {"candidate_id": f"protein_{i}", "binder_sequence": "ACDEFGHIK", "source_tool": "database_selection",
         "raw_metadata": {"studio_molecule_id": str(i), "components": [{"type": "protein", "chain": "A", "sequence": "ACDEFGHIK"}]}}
        for i in range(3)
    ])
    tasks = _optimized_prediction_tasks("esmfold", {"run_dir": str(run), "candidates_jsonl": str(run / "artifacts/normalized_candidates/candidates.jsonl")}, {})
    assert len(tasks) == 3
    assert {task["name"] for task in tasks} == {f"protein_{i}" for i in range(3)}


def test_prediction_tile_submits_full_compound_selection_without_batch_switch(monkeypatch, compound_file):
    from types import SimpleNamespace
    from streamlit.testing.v1 import AppTest
    from mn_molecular_design_studio import design_validation_ui, molecular_database
    path = compound_file()
    monkeypatch.setattr(molecular_database, "MolecularDatabase", lambda: SimpleNamespace(compound_options=lambda user: [
        {"path": str(path), "name": "Selected compounds", "compound_count": 3}]))
    monkeypatch.setattr(design_validation_ui, "render_molecule_editor", lambda *args, **kwargs: (
        {"protein_sequences": "ACDEFGHIK", "molecule_chain_options": [{"id": "A", "use_msa": True}]},
        {"optimization_mode": "off", "samples": 1, "cpu_cores": 4, "use_potentials": True}, True))
    app = AppTest.from_string('''
import streamlit as st
from mn_molecular_design_studio.design_validation_ui import render_batch_controls
from mn_molecular_design_studio.engine.workflow_engine.workflows.optimized_prediction import ENGINE_PROFILES
st.session_state['_molecular_database_user'] = {'id': 'owner'}
st.session_state['design-validation-boltz2-project-kind'] = 'Compound collection'
def submit(schema, workflow, parameters, project, settings, label):
    st.session_state['submitted'] = parameters
render_batch_controls('boltz2', ENGINE_PROFILES['boltz2'], 'project', 'prediction', {}, search_jobs=lambda **kwargs: {}, submit=submit)
''').run()
    assert not app.exception
    parameters = app.session_state['submitted']
    assert len(parameters['tasks']) == 3
    assert all(task['molecule_chain_options'] == [{"id": "A", "use_msa": True}] for task in parameters['tasks'])
    assert not any(element.label == "Batch mode" for element in app.expander)


def test_acceptance_tile_lists_native_rejected_designs_by_their_recorded_format(monkeypatch, design_source_factory):
    import json
    from streamlit.testing.v1 import AppTest
    from mn_molecular_design_studio import design_validation_ui
    from mn_molecular_design_studio.engine.workflow_engine.workflows.optimized_prediction import ENGINE_PROFILES

    run, _ = design_source_factory(passed=False)
    input_path = run / "input.json"
    input_data = json.loads(input_path.read_text())
    input_data["params"] = {"modality": "VHH", "target_type": "structured"}
    input_path.write_text(json.dumps(input_data))
    monkeypatch.setattr(design_validation_ui, "render_molecule_editor", lambda *args, **kwargs: ({}, {
        "optimization_mode": "off", "samples": 1, "cpu_cores": 4, "seed": 1,
    }, False))
    app = AppTest.from_string(f'''
import streamlit as st
from mn_molecular_design_studio.design_validation_ui import render_batch_controls
from mn_molecular_design_studio.engine.workflow_engine.workflows.optimized_prediction import ENGINE_PROFILES
def search_jobs(**kwargs):
    return {{"items": [{{"run_id": "source-run", "run_dir": {str(run)!r}, "tool": "bindcraft2", "campaign_name": "VHH campaign", "job_code": "VHH-1"}}], "next_cursor": None}}
render_batch_controls("boltz2", {{**ENGINE_PROFILES["boltz2"], "form_key": "acceptance-boltz2"}},
    "project", "validation", {{}}, search_jobs=search_jobs, submit=lambda *args: None,
    include_acceptance=True, existing_designs_only=True)
''').run()
    assert not app.exception
    campaigns = next(widget for widget in app.multiselect if widget.label == "Design campaigns")
    campaigns.set_value(["source-run"]).run()
    candidates = next(widget for widget in app.multiselect if widget.label == "Designs · VHH-1")
    assert candidates.value == ["passing-design"]
    assert any("native no pass recorded" in str(option) for option in candidates.options)
    binder_format = next(widget for widget in app.selectbox if widget.label == "Binder format")
    assert binder_format.value == "vhh"
