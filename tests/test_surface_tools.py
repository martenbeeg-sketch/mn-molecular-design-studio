from __future__ import annotations

import pytest
import pandas as pd

from mn_molecular_design_studio import surface_ui
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import JobPaths
from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import workflow_schema
from mn_molecular_design_studio.engine.workflow_engine.workflows import detection
from mn_molecular_design_studio import pipeline_graph
from mn_molecular_design_studio.engine.workflow_engine.workflows.graph_pipeline import node_request


PESTO_INTERFACES = ("protein", "dna_rna", "ion", "ligand", "lipid")
PDB = "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n"
PDB_COMPLEX = (
    "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n"
    "ATOM      2  CA  GLY B   1       1.000   0.000   0.000  1.00  0.00           C\n"
    "HETATM    3  C1  LIG C   1       2.000   0.000   0.000  1.00  0.00           C\nEND\n"
)


def test_pesto_workflow_schema_exposes_supported_interface_classes():
    parameters = {item["name"]: item for item in workflow_schema("detection.pesto")["parameters"]}

    assert parameters["interface"]["default"] == "protein"
    assert all(interface in parameters["interface"]["type"] for interface in PESTO_INTERFACES)


def test_pesto_residue_table_preserves_insertion_codes():
    table = surface_ui._normalized_score_table(
        pd.DataFrame({"chain": ["A"], "residue": ["42B"], "amino_acid": ["G"], "score": [0.9]}),
        "pesto",
        None,
    )

    assert table is not None
    assert table.iloc[0]["Residue"] == "42B"
    assert table.iloc[0]["Residue ID"] == "A42B"


def test_surf2spot_residue_numbers_are_not_assigned_to_an_arbitrary_chain(tmp_path):
    target = tmp_path / "two-chains.pdb"
    target.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      2  CA  ALA B   1       1.000   0.000   0.000  1.00  0.00           C\nEND\n",
        encoding="utf-8",
    )

    table = surface_ui._normalized_score_table(pd.DataFrame({"aa_id": [1], "score": [0.9]}), "surf2spot", target)

    assert table is None


def test_combined_table_counts_duplicate_method_rows_once():
    pesto = pd.DataFrame({
        "Method": ["PeSTo", "PeSTo", "PeSTo"],
        "Chain": ["A", "A", "A"], "Residue": ["1", "1", "2"], "Residue ID": ["A1", "A1", "A2"],
        "Score": [0.7, 0.9, 0.2],
    })
    scannet = pd.DataFrame({
        "Method": ["ScanNet", "ScanNet"],
        "Chain": ["A", "A"], "Residue": ["1", "2"], "Residue ID": ["A1", "A2"],
        "Score": [4, 8],
    })

    consensus = surface_ui._combined_surface_consensus([pesto, scannet])

    assert len(consensus) == 2
    assert set(consensus["Methods"]) == {2}
    assert set(consensus["Scoring methods"]) == {"PeSTo, ScanNet"}


def test_masif_integer_residue_ids_join_other_tool_scores(monkeypatch, tmp_path):
    from mn_molecular_design_studio.engine.workflow_engine.core import detection_annotations

    run_dir = tmp_path / "surface-run"
    surface = run_dir / "artifacts/masif_seed/pred_surfaces/target.ply"
    prediction = run_dir / "artifacts/masif_seed/pred_data/pred_target.npy"
    surface.parent.mkdir(parents=True)
    prediction.parent.mkdir(parents=True)
    surface.write_text("placeholder", encoding="utf-8")
    prediction.write_bytes(b"placeholder")
    target = tmp_path / "target.pdb"
    target.write_text(PDB, encoding="utf-8")
    monkeypatch.setattr(
        detection_annotations,
        "_masif_residue_scores",
        lambda *_args: pd.DataFrame({"chain": ["A"], "residue": [1], "amino_acid": ["A"], "score": [0.8]}),
    )

    masif = surface_ui._masif_score_table(run_dir, target)
    pesto = pd.DataFrame({
        "Method": ["PeSTo"], "Chain": ["A"], "Residue": ["1"], "Residue ID": ["A1"], "Score": [0.9],
    })
    masif["Method"] = "MaSIF-seed"
    consensus = surface_ui._combined_surface_consensus([pesto, masif])

    assert masif.iloc[0]["Residue"] == "1"
    assert len(consensus) == 1
    assert consensus.iloc[0]["Methods"] == 2
    assert consensus.iloc[0]["Scoring methods"] == "MaSIF-seed, PeSTo"


@pytest.mark.parametrize("interface", PESTO_INTERFACES)
def test_pesto_runner_passes_selected_interface_to_container(monkeypatch, tmp_path, interface):
    target = tmp_path / "target.pdb"
    target.write_text(PDB, encoding="utf-8")
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    job = JobPaths(task_group="detection", run_id="surface-test", run_dir=run_dir)
    steps: list[dict] = []
    monkeypatch.setattr(detection, "_run_shell_steps", lambda _run_dir, value, **_kwargs: steps.extend(value) or 0)
    monkeypatch.setattr(detection, "_collect_artifacts", lambda *_args, **_kwargs: [])
    monkeypatch.setattr(detection, "finish_job", lambda *_args, **_kwargs: None)

    detection.run_pesto(target, ["A"], existing_job=job, interface=interface)

    command = steps[0]["command"]
    interface_index = command.index("--interface")
    assert command[interface_index + 1] == interface


def test_pesto_runner_rejects_unknown_interface(tmp_path):
    with pytest.raises(ValueError, match="PeSTo interface"):
        detection.run_pesto(tmp_path / "target.pdb", ["A"], interface="carbohydrate")  # type: ignore[arg-type]


def test_surface_analysis_workflow_schema_exposes_all_method_modes():
    parameters = {item["name"]: item for item in workflow_schema("detection.surface_analysis")["parameters"]}

    assert parameters["target_pdb"]["required"]
    assert parameters["chain_ids"]["required"]
    assert all(mode in parameters["pesto_interface"]["type"] for mode in PESTO_INTERFACES)
    assert all(mode in parameters["scannet_mode"]["type"] for mode in ("interface", "epitope", "idp"))
    assert parameters["masif_chain_id"]["required"] is False


def test_surface_analysis_pipeline_node_connects_target_to_results():
    target, analysis, results = (pipeline_graph.NODE_CATALOG[kind] for kind in ("target", "surface_analysis", "results"))
    assert results["inputs"][0]["label"] == "Results"
    graph = pipeline_graph.new_graph("blank")
    graph["nodes"] = [
        {"id": "target", "type": "studio", "position": {"x": 0, "y": 0}, "data": {"kind": "target", "label": target["label"], "params": {"target_pdb": "/tmp/target.pdb", "target_chains": "A,B", "hotspots": ""}}},
        {"id": "surface", "type": "studio", "position": {"x": 300, "y": 0}, "data": {"kind": "surface_analysis", "label": analysis["label"], "params": {field["name"]: field["default"] for field in analysis["fields"]}}},
        {"id": "results", "type": "studio", "position": {"x": 600, "y": 0}, "data": {"kind": "results", "label": results["label"], "params": {}}},
    ]
    graph["edges"] = [
        {"id": "target-surface", "source": "target", "target": "surface", "sourceHandle": "out", "targetHandle": "in"},
        {"id": "surface-results", "source": "surface", "target": "results", "sourceHandle": "out", "targetHandle": "in"},
    ]
    plan = pipeline_graph.validate_graph(graph, require_inputs=True)
    assert plan["order"] == ["target", "surface", "results"]

    params = {field["name"]: field["default"] for field in analysis["fields"]}
    params.update(pesto_interface="ligand", scannet_mode="epitope", masif_chain_id="B")
    request = node_request("surface_analysis", params, {"target": {"target_pdb": "/tmp/target.pdb", "target_chains": ["A", "B"]}}, "PPI surface pass")
    assert request["workflow"] == "detection.surface_analysis"
    assert request["parameters"]["pesto_interface"] == "ligand"
    assert request["parameters"]["scannet_mode"] == "epitope"
    assert request["parameters"]["masif_chain_id"] == "B"


def test_surface_analysis_input_and_empty_results_views_render(tmp_path):
    from streamlit.testing.v1 import AppTest

    run_dir = tmp_path / "surface-run"
    run_dir.mkdir()
    (run_dir / "result.json").write_text(
        '{"success": false, "outputs": {"tool_status": {"pesto": {"success": false}, '
        '"scannet": {"success": false}, "surf2spot": {"success": false}, '
        '"masif_seed": {"success": false}}}}',
        encoding="utf-8",
    )
    app = AppTest.from_string(
        """
from pathlib import Path
from mn_molecular_design_studio.surface_ui import render_setup, render_surface_analysis_results
from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import workflow_schema
render_setup(workflow_schema('detection.surface_analysis'), owner_id='0123456789abcdef0123456789abcdef', project_selector=lambda _key: None, submit=lambda *args: None)
render_surface_analysis_results('empty-results', Path(%r), [], {'params': {}})
""" % str(run_dir)
    ).run()

    assert not app.exception
    assert {button.label for button in app.button} >= {
        "Protein–protein", "DNA / RNA", "Ion", "Ligand", "Lipid", "Interface", "Epitope", "Disordered region"
    }


def test_surface_analysis_runner_uses_one_job_and_reports_each_method(monkeypatch, tmp_path):
    target = tmp_path / "target.pdb"
    target.write_text(PDB_COMPLEX, encoding="utf-8")
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    (run_dir / "artifacts").mkdir()
    job = JobPaths(task_group="detection", run_id="bundle-test", run_dir=run_dir)
    created = []
    finished = []
    calls = []

    monkeypatch.setattr(detection, "normalize_gpu_device", lambda _device: "0")
    monkeypatch.setattr(detection, "gpu_queue_resource", lambda _device: "gpu:0")
    monkeypatch.setattr(detection, "create_job", lambda *args, **kwargs: created.append((args, kwargs)) or job)

    def tool_stub(tool):
        def invoke(*args, **kwargs):
            calls.append((tool, args, kwargs))
            status_path = run_dir / "artifacts/surface_analysis/tool_status.json"
            statuses = detection.read_json(status_path) if status_path.is_file() else {}
            statuses[tool] = {"success": True, "metrics": {"return_code": 0}}
            detection.write_json(status_path, statuses)
            return run_dir
        return invoke

    for name, tool in (("run_pesto", "pesto"), ("run_scannet", "scannet"), ("run_surf2spot", "surf2spot"), ("run_masif_seed", "masif_seed")):
        monkeypatch.setattr(detection, name, tool_stub(tool))
    monkeypatch.setattr(detection, "finish_job", lambda *args, **kwargs: finished.append((args, kwargs)))

    result = detection.run_surface_analysis(
        target,
        ["A"],
        pesto_interface="ligand",
        scannet_mode="epitope",
        masif_chain_id="A",
    )

    assert result == run_dir
    assert len(created) == len(finished) == 1
    assert [call[0] for call in calls] == ["pesto", "scannet", "surf2spot", "masif_seed"]
    assert all(call[2]["existing_job"] is job and call[2]["_defer_finish"] is True for call in calls)
    assert calls[0][1][0] == calls[1][1][0]
    selected_target = calls[1][1][0]
    selected_text = selected_target.read_text(encoding="utf-8")
    assert " A   1" in selected_text and " B   1" not in selected_text and "HETATM" not in selected_text
    assert finished[0][0][1] is True
    assert set(finished[0][0][2]["outputs"]["tool_status"]) == {"pesto", "scannet", "surf2spot", "masif_seed"}


def _atom(serial, chain, number, insertion="", x=0):
    return f"ATOM  {serial:5d}  CA  ALA {chain}{number:4d}{insertion or ' '}   {x:8.3f}{0:8.3f}{0:8.3f}  1.00  0.00           C\n"


def test_surf2spot_restores_author_ids_for_each_chain(monkeypatch, tmp_path):
    target = tmp_path / "target.pdb"
    target.write_text(_atom(1, "A", 10) + _atom(2, "A", 20) + _atom(3, "B", 42) + _atom(4, "B", 42, "B") + "END\n")
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    job = JobPaths(task_group="detection", run_id="surf", run_dir=run_dir)
    monkeypatch.setattr(detection, "_run_shell_steps", lambda *_args, **_kwargs: 1)
    monkeypatch.setattr(detection, "finish_job", lambda *_args, **_kwargs: None)

    detection.run_surf2spot(target, existing_job=job)

    numbered = (run_dir / "artifacts/surf2spot/input/target.pdb").read_text().splitlines()
    assert [(line[21], line[22:27].strip()) for line in numbered if line.startswith("ATOM")] == [("A", "1"), ("A", "2"), ("B", "1"), ("B", "2")]
    original = surface_ui._target_structure("surf2spot", run_dir, {})
    assert original.read_text() == target.read_text()
    table = surface_ui._normalized_score_table(
        pd.DataFrame({"aa_id": [1, 2], "score": [0.9, 0.3]}), "surf2spot", original,
        source_path=run_dir / "artifacts/surf2spot/predict/target_B.csv",
    )
    assert list(table["Residue ID"]) == ["B42", "B42B"]


@pytest.mark.parametrize("tool", ["pesto", "scannet", "surf2spot", "masif_seed"])
def test_zero_exit_without_predictions_is_a_failed_method(monkeypatch, tmp_path, tool):
    target = tmp_path / "target.pdb"
    target.write_text(PDB)
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    job = JobPaths(task_group="detection", run_id="no-output", run_dir=run_dir)
    finished = []
    monkeypatch.setattr(detection, "_run_shell_steps", lambda *_args, **_kwargs: 0)
    monkeypatch.setattr(detection, "finish_job", lambda *_args: finished.append(_args))
    kwargs = {"existing_job": job}
    if tool in {"pesto", "scannet"}:
        kwargs["chain_ids"] = ["A"]
    elif tool == "masif_seed":
        kwargs["chain_id"] = "A"

    getattr(detection, f"run_{tool}")(target, **kwargs)

    assert finished[0][1] is False
    assert "prediction files" in finished[0][2]["error"]


def test_bundle_keeps_every_command_and_both_output_streams(monkeypatch, tmp_path):
    import sys

    monkeypatch.setattr(detection, "apply_docker_cpu_limits_to_steps", lambda _run, steps: steps)
    monkeypatch.setattr(detection, "update_status", lambda *_args: None)
    for name in ("pesto", "scannet", "surf2spot", "masif_seed"):
        command = [sys.executable, "-c", f"import sys; print('{name} output'); print('{name} error', file=sys.stderr)"]
        assert detection._run_shell_steps(tmp_path, [{"name": name, "command": command}], append=True) == 0

    commands = detection.read_json(tmp_path / "command.json")["steps"]
    assert [step["name"] for step in commands] == ["pesto", "scannet", "surf2spot", "masif_seed"]
    for step in commands:
        assert f"{step['name']} output" in (tmp_path / "stdout.log").read_text()
        assert f"{step['name']} error" in (tmp_path / "stderr.log").read_text()


@pytest.mark.parametrize("interface", PESTO_INTERFACES)
def test_pesto_modes_score_selected_protein_without_partners(monkeypatch, tmp_path, interface):
    target = tmp_path / "complex.pdb"
    target.write_text(PDB_COMPLEX)
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    job = JobPaths(task_group="detection", run_id="pesto-selection", run_dir=run_dir)
    monkeypatch.setattr(detection, "_run_shell_steps", lambda *_args, **_kwargs: 1)
    monkeypatch.setattr(detection, "finish_job", lambda *_args: None)

    detection.run_pesto(target, ["A"], interface=interface, existing_job=job)

    staged = (run_dir / "artifacts/pesto/input/target.pdb").read_text()
    assert " A   1" in staged
    assert " B   1" not in staged and "HETATM" not in staged


@pytest.mark.parametrize("workflow,selection", [
    ("detection.pesto", {"chain_ids": ["B"]}),
    ("detection.scannet", {"chain_ids": ["B"]}),
    ("detection.masif_seed", {"chain_id": "A;touch /tmp/bad"}),
    ("detection.surface_analysis", {"chain_ids": ["A"], "masif_chain_id": "B"}),
])
def test_surface_preflight_rejects_invalid_chains(tmp_path, workflow, selection):
    from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import WORKFLOWS, _validate_workflow_inputs
    from mn_molecular_design_studio.engine.workflow_engine.services.local_automation import AutomationError

    target = tmp_path / "target.pdb"
    target.write_text(PDB)
    with pytest.raises(AutomationError, match="chain"):
        _validate_workflow_inputs(WORKFLOWS[workflow], {"target_pdb": target, **selection})


def test_masif_projection_preserves_insertion_codes(tmp_path):
    import numpy as np
    from mn_molecular_design_studio.engine.workflow_engine.core.detection_annotations import _masif_residue_scores

    target = tmp_path / "target.pdb"
    target.write_text(_atom(1, "A", 42, x=0) + _atom(2, "A", 42, "B", x=10))
    surface = tmp_path / "surface.ply"
    surface.write_text("ply\nformat ascii 1.0\nelement vertex 2\nproperty float x\nproperty float y\nproperty float z\nend_header\n0 0 0\n10 0 0\n")
    prediction = tmp_path / "pred.npy"
    np.save(prediction, np.array([0.2, 0.9]))

    table = _masif_residue_scores(surface, prediction, target)

    assert dict(zip(table["residue"], table["score"])) == {"42": 0.2, "42B": 0.9}


def test_score_profiles_follow_numeric_residue_order_and_drop_infinity():
    table = surface_ui._normalized_score_table(
        pd.DataFrame({"chain": [" A "] * 5, "residue": [1, 2, 10, "2B", 30], "score": [0.2, 0.4, 0.1, 0.7, float("inf")]}),
        "pesto", None,
    )
    assert list(surface_ui._ordered_score_table(table)["Residue ID"]) == ["A1", "A2", "A2B", "A10"]
