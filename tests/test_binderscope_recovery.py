import copy
from pathlib import Path

import pytest

from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json, write_json
from mn_molecular_design_studio.engine.workflow_engine.workflows import binderscope as workflow
from mn_molecular_design_studio.evaluations import file_digest


@pytest.fixture
def completed_report(tmp_path):
    run = tmp_path / "analysis" / "saved-run"
    inputs = run / "input"
    inputs.mkdir(parents=True)
    source = tmp_path / "authorized"
    source.mkdir()
    structure = source / "design.pdb"
    structure.write_text("saved exact coordinates")
    (source / "candidates.jsonl").write_text("")
    settings = {"component_format": "linear protein", "b_factor_interpretation": "plddt_0_100"}
    selector = {"molecule_id": "design", "structure_id": "complex", "group_id": "group",
                "structure_digest": file_digest(structure), "binder_chains": ["B"], "target_chains": ["A"],
                "assessment": settings}
    selected = {"candidate_id": "candidate", "source": {key: selector[key] for key in ("molecule_id", "structure_id", "group_id")},
                "source_structure_digest": file_digest(structure), "input_structure_digest": file_digest(structure),
                "assessment": settings, "binder_chains": ["B"], "target_chains": ["A"]}
    candidate = {"candidate_id": "candidate", "complex_pdb": str(structure), "binder_chains": ["B"], "target_chains": ["A"],
                 "raw_metadata": {"complex_analysis_selection": selected}}
    (inputs / "candidate.pdb").write_bytes(structure.read_bytes())
    parameters = {"cpu_cores": 32, "use_recorded_settings": True, "component_format": "unspecified", "bfactors_are_plddt": False}
    request = {**parameters, "mode": "saved", "adapter_sha256": file_digest(workflow.RUNTIME),
               "context_sha256": file_digest(workflow.RUNTIME.with_name("binderscope_context.py")),
               "items": [workflow._saved_request_item(candidate, inputs, use_recorded_settings=True,
                                                      component_format="unspecified", bfactors_are_plddt=False)]}
    vendor = Path(workflow.__file__).resolve().parents[3] / "vendor/binderscope/SOURCE_MANIFEST.json"
    report = {"schema": "studio.binderscope_report.v1", "mode": "saved", "upstream_revision": workflow.REVISION,
              "source_manifest_sha256": file_digest(vendor), "adapter_sha256": request["adapter_sha256"],
              "target_context_sha256": request["context_sha256"], "results": [{"candidate_id": "candidate",
              "source": {**request["items"][0], "residue_map": {}}, "analysis": {"dsasa": 42}}]}
    write_json(run / "input.json", {"tool": "binderscope", "inputs": {"source_run_dir": str(source), "candidates_jsonl": str(source / "candidates.jsonl")},
        "params": {**parameters, "method_version": workflow.METHOD_VERSION, "upstream_revision": workflow.REVISION,
                   "protocol_revision": "original-protocol-not-new-workflow-hash"}})
    write_json(run / "input/request.json", request)
    write_json(run / "command.json", {"image_digest": "original-container-image"})
    write_json(run / "metadata.json", {"status": "running", "scheduler_status": "failed"})
    write_json(run / "artifacts/analysis/binderscope/binderscope_report.json", report)
    return run, candidate, request, report, {"sources": [selector]}, parameters


def test_completed_report_retains_original_scientific_protocol_and_never_recalculates(completed_report, monkeypatch):
    run, candidate, request, report, selection, parameters = completed_report
    monkeypatch.setattr(workflow, "_saved_inputs", lambda *args, **kwargs: [copy.deepcopy(candidate)])
    monkeypatch.setattr(workflow.subprocess, "run", lambda *args, **kwargs: pytest.fail("Recovery must not execute Docker"))
    def finish(job, staged, recovered, **settings):
        assert recovered == report
        assert settings["recorded_protocol"] == "original-protocol-not-new-workflow-hash"
        assert settings["image_digest"] == "original-container-image"
        assert staged[0]["complex_pdb"] == str(run / "input/candidate.pdb")
        return job.run_dir
    monkeypatch.setattr(workflow, "_finish_report", finish)
    assert workflow._recover_saved_report(run, selection, parameters) == run
    receipt = read_json(run / "binderscope_recovery.json")
    assert receipt["designs"] == 1 and receipt["calculation_repeated"] is False
    assert receipt["original_container_exit_observed"] is False


@pytest.mark.parametrize("damage", ["duplicate", "missing", "extra", "source", "settings", "runtime", "file", "request"])
def test_report_recovery_rejects_incomplete_or_mismatched_science(completed_report, damage):
    run, candidate, request, report, selection, parameters = completed_report
    if damage == "duplicate":
        report["results"].append(copy.deepcopy(report["results"][0]))
    elif damage == "missing":
        report["results"] = []
    elif damage == "extra":
        report["results"].append({**report["results"][0], "candidate_id": "unselected"})
    elif damage == "source":
        report["results"][0]["source"]["source"] = {"molecule_id": "other"}
    elif damage == "settings":
        report["results"][0]["source"]["assessment"] = {"component_format": "cyclic peptide"}
    elif damage == "runtime":
        report["target_context_sha256"] = "wrong-version"
    elif damage == "file":
        (run / "input/candidate.pdb").write_text("changed coordinates")
    elif damage == "request":
        request["items"][0]["binder_chains"] = ["Z"]
    with pytest.raises(ValueError):
        workflow._validate_saved_report(run, report, request, [candidate])


def test_recovery_rejects_selection_replacement_without_writing_results(completed_report, monkeypatch):
    run, candidate, request, report, selection, parameters = completed_report
    monkeypatch.setattr(workflow, "_saved_inputs", lambda *args, **kwargs: [candidate])
    selection["sources"][0]["structure_id"] = "different-complex"
    with pytest.raises(ValueError, match="selection"):
        workflow._recover_saved_report(run, selection, parameters)
    assert not (run / "binderscope_recovery.json").exists()
