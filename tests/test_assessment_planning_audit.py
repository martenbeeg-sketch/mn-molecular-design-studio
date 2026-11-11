"""Regression contracts for audit findings AUD-004/006/011/012."""
import copy
import itertools

import pytest

from mn_molecular_design_studio.complex_analysis_batch import evaluation_coverage


def evidence(status="completed", *, identifier="saved", **settings):
    context = {"source_structure_id": "s", "source_structure_digest": "a" * 64,
               "model_id": "model-1", "binder_chains": ["B"], "target_chains": ["A"],
               "settings": {"contact_distance": 5.0, "method_version": "v2", **settings}, "status": status}
    return {"context": context}, {"id": identifier, "evaluator": "biotite_interface", "status": status}


def sample_with(*records):
    return {"structure_id": "s", "digest": "a" * 64, "model_id": "model-1", "associated": list(records)}


def coverage(sample, **kwargs):
    return evaluation_coverage(sample, "analysis.biotite_interface", parameters={"contact_distance": 5.0},
                               binder_chains=["B"], target_chains=["A"], **kwargs)


def test_failed_attempt_cannot_obscure_success_by_row_order():
    values = [evidence("failed"), evidence("completed", identifier="success")]
    states = [coverage(sample_with(*order))[0] for order in itertools.permutations(values)]
    assert len(set(states)) == 1


def test_unversioned_or_obsolete_evaluation_is_not_silently_reusable():
    observation, descriptor = evidence()
    observation["context"]["settings"] = {"contact_distance": 5.0}
    descriptor["engine_version"] = "obsolete-engine-version"
    assert coverage(sample_with((observation, descriptor)))[0] != "Already evaluated · compatible settings"


@pytest.mark.parametrize("changed", [
    {"method_version": "v1"}, {"protocol_revision": "old"}, {"image_digest": "sha256:old"},
])
def test_requested_protocol_versions_must_match(changed):
    protocol = {"method_version": "v2", "protocol_revision": "current", "image_digest": "sha256:current"}
    record = evidence(**{**protocol, **changed})
    assert coverage(sample_with(record), expected_protocol=protocol)[0] == "Existing Evaluation · different settings"


def test_success_precedes_active_and_failed_with_deterministic_selection():
    protocol = {"method_version": "v2"}
    values = [evidence("failed"), evidence("collecting"), evidence(identifier="b"), evidence(identifier="a")]
    for order in itertools.permutations(values):
        assert coverage(sample_with(*order), expected_protocol=protocol) == ("Already evaluated · compatible settings", "a")


def test_runtime_settings_are_ignored_but_unrequested_scientific_settings_are_not():
    protocol = {"method_version": "v2"}
    runtime = evidence(workers=8, cpu_cores=16, source_run_dir="/elsewhere")
    assert coverage(sample_with(runtime), expected_protocol=protocol)[0] == "Already evaluated · compatible settings"
    unexpected = evidence(contact_distance=5.0, new_scientific_cutoff=10.0)
    assert coverage(sample_with(unexpected), expected_protocol=protocol)[0] == "Existing Evaluation · different settings"


@pytest.mark.parametrize("field,value", [("model_id", "another"), ("source_structure_digest", "b" * 64),
                                        ("binder_chains", ["A"]), ("target_chains", ["B"])])
def test_reuse_retains_exact_source_and_chain_scope(field, value):
    record = evidence()
    record[0]["context"][field] = value
    assert coverage(sample_with(record), expected_protocol={"method_version": "v2"})[0] != "Already evaluated · compatible settings"


def test_protocol_is_published_and_persisted_through_registered_workflow_contract(tmp_path, monkeypatch):
    from types import SimpleNamespace
    from mn_molecular_design_studio.complex_analysis_batch import assessment_protocol
    from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import workflow_schema
    from mn_molecular_design_studio.engine.workflow_engine.workflows import analysis

    expected = assessment_protocol("analysis.biotite_interface")
    assert workflow_schema("analysis.biotite_interface")["assessment_protocol"] == expected
    job = SimpleNamespace(run_dir=tmp_path / "job")
    captured = {}
    def create_job(*args, **kwargs):
        captured.update(kwargs)
        return job
    monkeypatch.setattr(analysis, "create_job", create_job)
    monkeypatch.setattr(analysis, "update_status", lambda *a, **kw: None)
    analysis._analysis_job(tmp_path / "source", tmp_path / "candidates.jsonl", tool="biotite_interface",
                           params={"contact_distance": 5.0}, note="isolated fixture", mode="docker-batch")
    assert captured["params"]["assessment_protocol"] == expected


def test_unknown_protocol_is_visible_and_cannot_hide_failure_behind_other_settings():
    observation, descriptor = evidence("failed", contact_distance=4.0)
    assert coverage(sample_with((observation, descriptor)))[0] == "Existing Evaluation · different settings"
    assert "unknown" in coverage(sample_with(evidence()), expected_protocol={"image_digest": ""})[1]
