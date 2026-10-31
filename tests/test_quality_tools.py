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
