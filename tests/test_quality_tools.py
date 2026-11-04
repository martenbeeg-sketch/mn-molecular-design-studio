import inspect
import json
import os
from pathlib import Path
import shutil
import subprocess
from types import SimpleNamespace

import pytest

from mn_molecular_design_studio.quality_tools import TOOLS
from mn_molecular_design_studio.engine.workflow_engine.workflows import quality, quality_runner
from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import WORKFLOWS, _get_function, validate_workflow_request


def test_sequence_properties_record_ph_and_do_not_discard_unknown_residues():
    metrics = quality_runner.sequence_metrics("ACDEFGHIKLMNPQRSTVWY", 5.5)
    assert metrics["qc_sequence_length"] == 20
    assert metrics["qc_sequence_ph"] == 5.5
    assert sum(value for key, value in metrics.items() if key.startswith("qc_sequence_fraction_")) == pytest.approx(1.)
    with pytest.raises(ValueError, match="unknown residues"):
        quality_runner.sequence_metrics("ACDX")


def test_workers_do_not_exceed_scheduler_grant(monkeypatch):
    from mn_molecular_design_studio.engine.workflow_engine.workflows.analysis import _allocated_cpu_cores
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_JOB_CPU_CORES", "2")
    assert _allocated_cpu_cores(16) == 2


def test_pipeline_filter_can_select_cluster_representatives_without_treating_cluster_ids_as_quality():
    from mn_molecular_design_studio.engine.workflow_engine.workflows.graph_pipeline import filter_candidates
    candidates = [{"candidate_id": "one", "metrics": {"diversity_representative": True, "diversity_cluster": 1}},
                  {"candidate_id": "two", "metrics": {"diversity_representative": False, "diversity_cluster": 1}}]
    selected = filter_candidates(candidates, {"metric": "common.diversity.representative", "operator": ">=", "threshold": 1,
                                               "passing_only": False, "higher_is_better": True, "limit": 0})
    assert [row["candidate_id"] for row in selected] == ["one"]


@pytest.mark.parametrize("tool", TOOLS)
def test_analysis_tools_have_batch_api_schemas_and_pipeline_adapters(tool, tmp_path, design_source_factory):
    from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import workflow_catalog
    assert "analysis." + tool in {row["workflow"] for row in workflow_catalog()}
    from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG
    from mn_molecular_design_studio.engine.workflow_engine.workflows.graph_pipeline import node_request
    source, _ = design_source_factory()
    spec = NODE_CATALOG[tool]
    request = node_request(tool, {f["name"]: f["default"] for f in spec["fields"]},
        {"run_dir": str(source), "candidates_jsonl": str(source / "artifacts/normalized_candidates/candidates.jsonl")}, "QC")
    function = _get_function(WORKFLOWS[request["workflow"]])
    inspect.signature(function).bind(**request["parameters"])
    path = tmp_path / "request.json"
    path.write_text(json.dumps(request))
    validated = validate_workflow_request(path)
    assert validated["valid"] is True


def test_collection_is_one_container_and_keeps_failed_rows(tmp_path, monkeypatch, design_source_factory):
    source, first = design_source_factory()
    candidate_path = source / "artifacts/normalized_candidates/candidates.jsonl"
    candidate_path.write_text(json.dumps(first) + "\n" + json.dumps({**first, "candidate_id": "second"}) + "\n")
    run = tmp_path / "analysis/job"
    run.mkdir(parents=True)
    monkeypatch.setattr(quality, "_analysis_job", lambda *a, **k: SimpleNamespace(run_dir=run))
    monkeypatch.setattr(quality, "_allocated_cpu_cores", lambda *a: 2)
    monkeypatch.setattr(quality, "apply_docker_cpu_limits_to_steps", lambda run, steps: steps)
    monkeypatch.setattr(quality, "docker_user_args", lambda: [])
    finishes, commands = [], []
    monkeypatch.setattr(quality, "finish_job", lambda *a: finishes.append(a))
    def fake(command, **kwargs):
        if command[:3] == ["docker", "image", "inspect"]:
            return SimpleNamespace(returncode=0, stdout="sha256:test\n")
        commands.append(command)
        (run / "artifacts/analysis/batch_status.json").write_text(json.dumps({"items": [
            {"id": "item_000000", "status": "completed", "metrics": {"qc_sequence_length": 5, "qc_sequence_pi": 7., "qc_sequence_gravy": .1}},
            {"id": "item_000001", "status": "failed", "metrics": {}, "error": "unknown residues"}]}))
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(quality.subprocess, "run", fake)
    assert quality.run_sequence_qc(source, candidate_path) == run
    assert len(commands) == 1
    assert "--network" in commands[0]
    assert "sha256:test" in commands[0]  # Run the inspected immutable image, even if its tag moves.
    manifest = json.loads((run / "artifacts/analysis/batch_manifest.json").read_text())
    assert len(manifest["items"]) == 2
    assert manifest["settings"]["workers"] == 2
    summary = json.loads((run / "artifacts/analysis/qc_summary.json").read_text())
    assert summary["completed_count"] == summary["failed_count"] == 1
    assert finishes[0][1] is True
    rows = [json.loads(line) for line in (run / "artifacts/normalized_candidates/candidates.jsonl").read_text().splitlines()]
    assert rows[1]["metrics"]["qc_status"] == "failed"
    assert rows[0]["raw_metadata"]["source_candidate"]["metrics"] == first["metrics"]
    assert "native_pass_filters" not in rows[0]["metrics"]  # inherited metrics are lineage, not new measurements


def test_biotite_interface_runs_exact_candidate_model_and_chain_groups(tmp_path, monkeypatch, design_source_factory):
    source, candidate = design_source_factory()
    candidate_path = source / "artifacts/normalized_candidates/candidates.jsonl"
    second = {**candidate, "candidate_id": "other-model", "complex_pdb": "designed.cif"}
    candidate_path.write_text(json.dumps(candidate) + "\n" + json.dumps(second) + "\n")
    selected_model = source / "selected-model-2.cif"
    selected_model.write_bytes((source / "designed.cif").read_bytes() + b"# exact selected model\n")
    run = tmp_path / "analysis/job"
    run.mkdir(parents=True)
    monkeypatch.setattr(quality, "_analysis_job", lambda *a, **k: SimpleNamespace(run_dir=run))
    monkeypatch.setattr(quality, "_allocated_cpu_cores", lambda *a: 2)
    monkeypatch.setattr(quality, "apply_docker_cpu_limits_to_steps", lambda run, steps: steps)
    monkeypatch.setattr(quality, "docker_user_args", lambda: [])
    monkeypatch.setattr(quality, "finish_job", lambda *a: None)

    def fake(command, **kwargs):
        if command[:3] == ["docker", "image", "inspect"]:
            return SimpleNamespace(returncode=0, stdout="")
        (run / "artifacts/analysis/batch_status.json").write_text(json.dumps({"items": [
            {"id": "item_000000", "status": "completed", "metrics": {"qc_interface_contact_pairs": 12}},
        ]}))
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(quality.subprocess, "run", fake)
    result = quality.run_biotite_interface(
        source, candidate_path, candidate_ids=["passing-design"],
        selected_structure_path=selected_model, selected_model_id="model-2",
        binder_chains=["B"], target_chains=["A"], contact_distance=4.5,
    )

    manifest = json.loads((run / "artifacts/analysis/batch_manifest.json").read_text())
    assert result == run
    assert manifest["selection"] == {
        "candidate_ids": ["passing-design"], "structure_path": str(selected_model.resolve()),
        "model_id": "model-2", "binder_chains": ["B"], "target_chains": ["A"],
    }
    assert len(manifest["items"]) == 1
    assert manifest["items"][0]["candidate_id"] == "passing-design"
    assert manifest["items"][0]["binder_chains"] == ["B"]
    assert manifest["items"][0]["target_chains"] == ["A"]
    assert manifest["items"][0]["structure_digest"] == __import__("hashlib").sha256(selected_model.read_bytes()).hexdigest()
    outputs = [json.loads(line) for line in (run / "artifacts/normalized_candidates/candidates.jsonl").read_text().splitlines()]
    assert len(outputs) == 1
    assert outputs[0]["raw_metadata"]["complex_analysis_selection"]["model_id"] == "model-2"
    assert outputs[0]["raw_metadata"]["analysis_parameters"]["contact_distance"] == 4.5
    assert outputs[0]["binder_chains"] == ["B"] and outputs[0]["target_chains"] == ["A"]


def test_biotite_batch_preserves_per_candidate_model_identity(tmp_path, monkeypatch, design_source_factory):
    source, candidate = design_source_factory()
    candidate_path = source / "artifacts/normalized_candidates/candidates.jsonl"
    second = {**candidate, "candidate_id": "passing-design-2", "model_index": 2}
    candidate_path.write_text(json.dumps(candidate) + "\n" + json.dumps(second) + "\n")
    run = tmp_path / "analysis" / "batch-models"
    run.mkdir(parents=True)
    monkeypatch.setattr(quality, "_analysis_job", lambda *a, **k: SimpleNamespace(run_dir=run))
    monkeypatch.setattr(quality, "_allocated_cpu_cores", lambda *a: 2)
    monkeypatch.setattr(quality, "apply_docker_cpu_limits_to_steps", lambda run, steps: steps)
    monkeypatch.setattr(quality, "docker_user_args", lambda: [])
    monkeypatch.setattr(quality, "finish_job", lambda *a: None)

    def fake(command, **kwargs):
        if command[:3] == ["docker", "image", "inspect"]:
            return SimpleNamespace(returncode=0, stdout="")
        (run / "artifacts/analysis/batch_status.json").write_text(json.dumps({"items": [
            {"id": "item_000000", "status": "completed", "metrics": {}},
            {"id": "item_000001", "status": "completed", "metrics": {}},
        ]}))
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(quality.subprocess, "run", fake)
    quality.run_biotite_interface(
        source, candidate_path, candidate_ids=["passing-design", "passing-design-2"],
        selected_model_ids={"passing-design": "pose-A", "passing-design-2": "pose-B"},
        binder_chains=["B"], target_chains=["A"],
    )
    manifest = json.loads((run / "artifacts/analysis/batch_manifest.json").read_text())
    assert manifest["selection"]["model_ids"] == {"passing-design": "pose-A", "passing-design-2": "pose-B"}
    assert [item["model_id"] for item in manifest["items"]] == ["pose-A", "pose-B"]
    results = [json.loads(line) for line in (run / "artifacts/normalized_candidates/candidates.jsonl").read_text().splitlines()]
    assert [row["raw_metadata"]["complex_analysis_selection"]["model_id"] for row in results] == ["pose-A", "pose-B"]


def test_ligand_batch_requires_an_exact_selection_for_every_candidate(design_source_factory):
    source, candidate = design_source_factory()
    candidate_path = source / "artifacts/normalized_candidates/candidates.jsonl"
    second = {**candidate, "candidate_id": "candidate-two"}
    candidate_path.write_text(json.dumps(candidate) + "\n" + json.dumps(second) + "\n")
    with pytest.raises(ValueError, match="identify every selected candidate exactly once"):
        quality.run_plip(
            source, candidate_path, candidate_ids=[candidate["candidate_id"], "candidate-two"],
            ligand_selections={candidate["candidate_id"]: {
                "chain": "L", "residue_key": "L:1", "smiles": "CC",
            }},
        )


def test_biotite_interface_selection_rejects_unknown_candidates_and_overlapping_groups(tmp_path, design_source_factory):
    source, candidate = design_source_factory()
    candidate_path = source / "artifacts/normalized_candidates/candidates.jsonl"
    with pytest.raises(ValueError, match="unavailable in the source collection"):
        quality.run_biotite_interface(source, candidate_path, candidate_ids=["missing"],
                                      binder_chains=["B"], target_chains=["A"])
    with pytest.raises(ValueError, match="must not overlap"):
        quality.run_biotite_interface(source, candidate_path, candidate_ids=[candidate["candidate_id"]],
                                      binder_chains=["A"], target_chains=["A"])


def test_biotite_schema_exposes_contextual_selection_through_the_existing_workflow():
    from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import workflow_schema

    schema = workflow_schema("analysis.biotite_interface")
    parameters = {item["name"]: item for item in schema["parameters"]}
    assert "source_job" in schema["request_fields"]
    assert parameters["candidate_ids"]["type"] == "list[str] | None"
    assert parameters["selected_structure_path"]["type"] == "pathlib.Path | None"
    assert parameters["selected_model_id"]["type"] == "str | None"
    assert parameters["binder_chains"]["type"] == "list[str] | None"
    assert parameters["target_chains"]["type"] == "list[str] | None"


def test_biotite_workflow_request_validates_selected_model_and_chain_groups(tmp_path, monkeypatch, design_source_factory):
    from mn_molecular_design_studio.engine.workflow_engine.services import workflow_automation

    source, candidate = design_source_factory()
    selected = source / "designed.cif"
    monkeypatch.setattr(workflow_automation, "resolve_job", lambda reference, **kwargs: {
        "run_dir": str(source), "run_id": source.name, "task_group": "design", "project_id": "test-project",
    })
    request = {
        "schema_version": 1, "workflow": "analysis.biotite_interface",
        "source_job": f"design/{source.name}",
        "parameters": {
            "candidate_ids": [candidate["candidate_id"]],
            "selected_structure_path": str(selected), "selected_model_id": "model-1",
            "binder_chains": ["B"], "target_chains": ["A"], "contact_distance": 4.5,
        },
    }
    request_path = tmp_path / "workflow-request.json"
    request_path.write_text(json.dumps(request))

    result = validate_workflow_request(request_path)

    assert result["valid"] is True
    assert result["workflow"] == "analysis.biotite_interface"
    assert result["parameters"]["candidate_ids"] == [candidate["candidate_id"]]
    assert result["parameters"]["selected_structure_path"] == str(selected.resolve())
    assert result["parameters"]["selected_model_id"] == "model-1"
    assert result["parameters"]["binder_chains"] == ["B"]
    assert result["parameters"]["target_chains"] == ["A"]


def test_biotite_source_job_submission_keeps_existing_job_acl_and_project_boundary(tmp_path, monkeypatch):
    from fastapi import HTTPException
    from mn_molecular_design_studio import api, job_catalog

    source = tmp_path / "prediction" / "source"
    source.mkdir(parents=True)
    (source / "metadata.json").write_text(json.dumps({"status": "completed", "tool": "prediction"}))
    resolved = {"run_dir": str(source), "run_id": "source", "task_group": "prediction", "project_id": "source-project"}
    access = {"owner_id": "source-owner", "project_id": "source-project"}
    monkeypatch.setattr(api, "resolve_job", lambda reference, **kwargs: dict(resolved))
    monkeypatch.setattr(job_catalog, "get_job_access_record", lambda group, run: access)

    with pytest.raises(HTTPException) as denied:
        api._resolve_api_job("prediction/source", {"id": "unrelated-user", "role": "user"})
    assert denied.value.status_code == 404

    access["owner_id"] = USER_ID = "analysis-owner"
    user = {"id": USER_ID, "role": "user"}
    row = api._resolve_api_job("prediction/source", user)
    assert row["run_id"] == "source"
    with pytest.raises(HTTPException) as wrong_project:
        api._guard_source_job({"source_job": "prediction/source", "workflow": "analysis.biotite_interface", "parameters": {}},
                              user, project_id="different-project")
    assert wrong_project.value.status_code == 422
    api._guard_source_job({"source_job": "prediction/source", "workflow": "analysis.biotite_interface", "parameters": {}},
                          user, project_id="source-project")


def test_legacy_docking_source_selects_only_best_poses_and_api_allows_missing_candidate_contract(tmp_path, monkeypatch):
    from mn_molecular_design_studio import ligand_compounds
    from mn_molecular_design_studio.engine.workflow_engine.workflows.quality_ligands import docking_candidates
    from mn_molecular_design_studio.engine.workflow_engine.services import workflow_automation
    source = tmp_path / "vina/legacy"
    (source / "artifacts/input").mkdir(parents=True)
    (source / "artifacts/input/receptor.pdb").write_text("HEADER receptor\nEND\n")
    for index in (1, 2):
        (source / f"pose-{index}.pdbqt").write_text("REMARK pose\n")
    ligand_compounds.write_compound_set(source / "artifacts/input/compound_set.json", name="Probes",
        compounds=[{"compound_id": "CMP1", "compound_name": "Ethanol", "smiles": "CCO"}], owner_id="owner", project_id="p", created_at="now")
    (source / "artifacts/docking_scores.csv").write_text("compound_id,status,pose_rank,best_score_kcal_mol,pose_path\nCMP1,completed,1,-7,pose-1.pdbqt\nCMP1,completed,2,-8,pose-2.pdbqt\n")
    selected = docking_candidates(source)
    assert len(selected) == 1
    assert selected[0]["complex_pdb"].endswith("pose-1.pdbqt")
    assert selected[0]["raw_metadata"]["name"] == "Ethanol"
    monkeypatch.setattr(workflow_automation, "resolve_job", lambda ref, **kwargs: {"run_dir": str(source), "run_id": source.name})
    parameters = workflow_automation._resolve_source_job({}, "legacy", quality.run_plip, WORKFLOWS["analysis.plip"])
    assert parameters == {"source_run_dir": source.resolve()}
    with pytest.raises(workflow_automation.AutomationError, match="no normalized candidates"):
        workflow_automation._resolve_source_job({}, "legacy", quality.run_sequence_qc, WORKFLOWS["analysis.sequence_qc"])


def test_ligand_preparation_rejects_ambiguous_copies_and_uses_a_selected_chain(tmp_path):
    from mn_molecular_design_studio.engine.workflow_engine.workflows.quality_ligands import prepare
    from mn_molecular_design_studio.molecular_database_ingest import _record
    from mn_molecular_design_studio.molecular_database import MolecularDatabase
    from mn_molecular_design_studio import molecule_imports
    assets = Path(__file__).resolve().parents[1] / "mn_molecular_design_studio/assets/vina"
    compound = json.loads((assets / "imatinib.json").read_text())[0]
    candidate = {"complex_pdb": str(assets / "1iep.pdb"), "raw_metadata": {"ligand_smiles": compound["smiles"], "source_compound_id": "CMP1", "compound_name": "Imatinib"}}
    item = {"id": "item_000000", "candidate_id": "CMP1"}
    with pytest.raises(ValueError, match="multiple or no unambiguous"):
        prepare(item, candidate, tmp_path)
    candidate["raw_metadata"]["ligand_chain"] = "B"
    ligand, receptor = prepare(item, candidate, tmp_path)
    assert ligand.is_file() and receptor.is_file()
    from rdkit import Chem
    pose = next(iter(Chem.SDMolSupplier(str(ligand))))
    assert pose.GetNumAtoms() == Chem.MolFromSmiles(compound["smiles"]).GetNumAtoms()
    record = _record(candidate, tmp_path)
    ligand_components = [component for component in record["components"] if component["type"] == "small_molecule"]
    assert len(ligand_components) == 1
    # The topology-bearing analyzed complex reuses the imported compound identity.
    database = MolecularDatabase(tmp_path / "db")
    user = {"id": "owner", "role": "user"}
    original = database.import_records(user, molecule_imports.parse_smiles(compound["smiles"] + " CMP1"))
    analyzed = database.import_records(user, [{**record, "components": ligand_components, "files": []}])
    assert analyzed["molecule_ids"] == original["molecule_ids"]


def test_ligand_preparation_matches_long_mmcif_component_after_pdb_export(tmp_path):
    import gemmi
    from mn_molecular_design_studio.engine.workflow_engine.workflows.quality_ligands import prepare
    from mn_molecular_design_studio.molecule_imports import parse_structure
    from mn_molecular_design_studio.evaluations import file_digest

    assets = Path(__file__).resolve().parents[1] / "mn_molecular_design_studio/assets/vina"
    smiles = json.loads((assets / "imatinib.json").read_text())[0]["smiles"]
    structure = gemmi.read_structure(str(assets / "1iep.pdb"))
    for model in structure:
        for chain in model:
            for residue in chain:
                if chain.name == "B" and residue.name == "STI":
                    residue.name = "LIG_LONG"
    source = tmp_path / "long-component-id.cif"
    structure.make_mmcif_document().write_file(str(source))
    components = parse_structure(source.read_bytes(), source.name, include_ligands=True)[0]["components"]
    candidate = {"complex_pdb": str(source), "raw_metadata": {"components": components}}
    item = {"id": "item_000000", "candidate_id": "long-ccd"}

    prepare(item, candidate, tmp_path / "prepared", ligand_selection={
        "chain": "B", "residue_key": "B:202", "ccd": "LIG_LONG", "smiles": smiles,
    })

    assert item["ligand_identity"]["ccd"] == "LIG_LONG"
    assert item["ligand_identity"]["chain"] == "B"
    assert item["ligand_identity"]["residue_key"] == "B:202"
    assert item["source_structure_digest"] == file_digest(source)
    assert item["ligand_pose_digest"] != item["receptor_digest"]


def test_ligand_workflow_runs_exact_pose_through_existing_batch_and_records_provenance(tmp_path, monkeypatch):
    from mn_molecular_design_studio.molecule_imports import parse_structure
    from mn_molecular_design_studio.evaluations import file_digest

    source = tmp_path / "prediction" / "protein-ligand"
    structure_dir = source / "artifacts" / "output"
    structure_dir.mkdir(parents=True)
    assets = Path(__file__).resolve().parents[1] / "mn_molecular_design_studio/assets/vina"
    structure = structure_dir / "1iep.pdb"
    shutil.copyfile(assets / "1iep.pdb", structure)
    ligand_smiles = json.loads((assets / "imatinib.json").read_text())[0]["smiles"]
    components = parse_structure(structure.read_bytes(), structure.name, include_ligands=True)[0]["components"]
    for component in components:
        if component.get("type") == "small_molecule" and component.get("ccd") == "STI":
            component["smiles"] = ligand_smiles
    candidate = {
        "candidate_id": "saved-1iep", "complex_pdb": str(structure),
        "raw_metadata": {"components": components, "ligand_smiles": ligand_smiles,
                         "compound_id": "STI-571", "compound_name": "Imatinib"},
    }
    candidate_path = source / "artifacts" / "normalized_candidates" / "candidates.jsonl"
    candidate_path.parent.mkdir(parents=True)
    candidate_path.write_text(json.dumps(candidate) + "\n")
    run = tmp_path / "analysis" / "plip-run"
    run.mkdir(parents=True)
    monkeypatch.setattr(quality, "_analysis_job", lambda *a, **k: SimpleNamespace(run_dir=run))
    monkeypatch.setattr(quality, "_allocated_cpu_cores", lambda *a: 1)
    monkeypatch.setattr(quality, "apply_docker_cpu_limits_to_steps", lambda run, steps: steps)
    monkeypatch.setattr(quality, "docker_user_args", lambda: [])
    monkeypatch.setattr(quality, "finish_job", lambda *a: None)

    def fake(command, **kwargs):
        if command[:3] == ["docker", "image", "inspect"]:
            return SimpleNamespace(returncode=0, stdout="sha256:test-image\n")
        (run / "artifacts/analysis/batch_status.json").write_text(json.dumps({"items": [
            {"id": "item_000000", "status": "completed", "metrics": {"plip_hbond_count": 2}},
        ]}))
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(quality.subprocess, "run", fake)
    selection = {"chain": "A", "residue_key": "A:201", "ccd": "STI", "smiles": ligand_smiles}
    result = quality.run_plip(
        source, candidate_path, candidate_ids=["saved-1iep"], selected_structure_path=structure,
        selected_model_id="1", ligand_selection=selection,
    )

    manifest = json.loads((run / "artifacts/analysis/batch_manifest.json").read_text())
    item = manifest["items"][0]
    assert result == run
    assert manifest["selection"]["ligand_selection"] == selection
    assert item["source_structure_digest"] == file_digest(structure)
    assert item["model_id"] == "1"
    assert item["ligand_selection"] == selection
    assert item["ligand_identity"] == {"chain": "A", "residue_key": "A:201", "ccd": "STI",
                                        "smiles": item["ligand_identity"]["smiles"]}
    assert item["ligand_provenance"]["source_structure_digest"] == file_digest(structure)
    assert item["ligand_provenance"]["pose_coordinates_sha256"] == item["ligand_pose_digest"]
    assert item["ligand_provenance"]["receptor_coordinates_sha256"] == item["receptor_digest"]
    assert item["receptor_digest"] != item["ligand_pose_digest"]
    analyzed = json.loads((run / "artifacts/normalized_candidates/candidates.jsonl").read_text().splitlines()[0])
    assert analyzed["raw_metadata"]["ligand_provenance"] == item["ligand_provenance"]
    assert analyzed["metrics"]["qc_status"] == "completed"


@pytest.mark.parametrize("runner,tool", [
    (quality.run_plip, "plip"), (quality.run_pandamap, "pandamap"),
    (quality.run_posebusters, "posebusters"),
])
def test_ligand_workflow_adapters_accept_context_selection_and_reject_interface_groups(monkeypatch, tmp_path, runner, tool):
    captured = {}
    monkeypatch.setattr(quality, "run_quality", lambda *args, **kwargs: captured.update(tool=kwargs.pop("tool"), kwargs=kwargs) or tmp_path)
    selection = {"chain": "A", "residue_key": "A:201", "ccd": "STI", "smiles": "CC"}
    path = tmp_path / "structure.pdb"
    result = runner(tmp_path, tmp_path / "candidates.jsonl", candidate_ids=["candidate-1"],
                    selected_structure_path=path, selected_model_id="model-2", ligand_selection=selection)
    assert result == tmp_path
    assert captured == {"tool": tool, "kwargs": {
        "settings": {}, "cpu_cores": 4, "candidate_ids": ["candidate-1"],
        "selected_structure_path": path, "selected_model_id": "model-2", "ligand_selection": selection,
    }}
    with pytest.raises(ValueError, match="does not accept protein interface chain groups"):
        runner(tmp_path, tmp_path / "candidates.jsonl", binder_chains=["A"])


@pytest.mark.skipif(os.environ.get("MDS_QC_DOCKER_TESTS") != "1", reason="Set MDS_QC_DOCKER_TESTS=1 to exercise installed QC containers.")
@pytest.mark.parametrize("tool", [key for key, profile in TOOLS.items() if not profile.get("gpu")])
def test_real_container_batch(tool, tmp_path, design_source_factory):
    source, candidate = design_source_factory()
    work = tmp_path / "work"
    work.mkdir()
    # A complete local example provides a realistic surface/structure input.
    from mn_molecular_design_studio.rfdiffusion3_input import ASSETS
    protein = ASSETS / "5o45_cropped.pdb"
    shutil.copyfile(protein, work / "protein.pdb")
    shutil.copyfile(source / "designed.cif", work / "complex.cif")
    structure = "/work/complex.cif" if tool in {"biotite_interface", "interface_clustering"} else "/work/protein.pdb"
    items = [{"id": f"item_{i:06d}", "candidate_id": f"candidate-{i}", "sequence": "ACDEFGHIKLMNPQRSTVWY" * 3,
              "structure": structure, "binder_chains": ["B"] if structure.endswith(".cif") else ["A"],
              "target_chains": ["A"] if structure.endswith(".cif") else [], "target_key": "same-target"} for i in range(2)]
    if TOOLS[tool]["kind"] == "ligand":
        from mn_molecular_design_studio.engine.workflow_engine.workflows.quality_ligands import prepare
        assets = Path(__file__).resolve().parents[1] / "mn_molecular_design_studio/assets/vina"
        compound = json.loads((assets / "imatinib.json").read_text())[0]
        for item in items:
            prepare(item, {"complex_pdb": str(assets / "1iep.pdb"), "raw_metadata": {"ligand_smiles": compound["smiles"], "ligand_chain": "A"}}, work)
    manifest = {"schema": "studio.qc.v1", "tool": tool, "items": items, "settings": {"workers": 2, "ph": 7.4, "similarity_threshold": .5}, "output": "/work/output"}
    (work / "manifest.json").write_text(json.dumps(manifest))
    module = Path(quality_runner.__file__).parent
    result = subprocess.run(["docker", "run", "--rm", "--network", "none", "--entrypoint", "python3", "--cpus", "2", *quality.docker_user_args(),
        "-e", "OPENBLAS_NUM_THREADS=1", "-e", "OMP_NUM_THREADS=1", "-v", f"{work}:/work", "-v", f"{module}:/opt/studio-qc:ro", "-w", "/work",
        TOOLS[tool]["image"], "/opt/studio-qc/quality_runner.py", "--manifest", "/work/manifest.json"], capture_output=True, text=True, timeout=180)
    assert (work / "output/batch_status.json").is_file(), result.stdout + result.stderr
    status = json.loads((work / "output/batch_status.json").read_text())
    assert result.returncode == 0, result.stdout + result.stderr + json.dumps(status)
    assert status["finished"] == 2
    assert all(item["status"] == "completed" for item in status["items"]), status
