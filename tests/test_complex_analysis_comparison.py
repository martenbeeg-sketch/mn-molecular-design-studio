from __future__ import annotations

from pathlib import Path

from streamlit.testing.v1 import AppTest

from mn_molecular_design_studio.complex_analysis_comparison import (
    compatibility_report,
    protein_chain_numbering,
    saved_target_contact_rows,
    target_structure_fingerprint,
)


def _atom(serial: int, residue: str, chain: str, number: int, x: float, insertion: str = " ") -> str:
    return (
        f"ATOM  {serial:5d}  CA  {residue:>3s} {chain}{number:4d}{insertion:1s}   "
        f"{x:8.3f}{1.0:8.3f}{1.0:8.3f}{1.0:6.2f}{20.0:6.2f}           C"
    )


def _sample(*, engine="boltz2", target_numbering=None, target_identity="target-1", digest="same-coordinate"):
    return {
        "source_tool": engine,
        "complex_type": "protein-protein",
        "target_identity": target_identity,
        "target_chains": ["A"],
        "binder_chains": ["B"],
        "target_chains_sequence": ["AG"],
        "binder_chains_sequence": ["SS"],
        "target_numbering": target_numbering if target_numbering is not None else [
            {"sequence": "AG", "numbering": [{"number": 10, "insertion": ""}, {"number": 11, "insertion": ""}]}
        ],
        "structure_digest": digest,
        "model_id": "model-1",
    }


def test_chain_numbering_retains_auth_numbers_and_insertion_codes(tmp_path):
    path = tmp_path / "complex.pdb"
    path.write_text("\n".join([
        _atom(1, "ALA", "A", 10, 1.0),
        _atom(2, "GLY", "A", 10, 2.0, "A"),
        _atom(3, "SER", "B", 1, 8.0),
        "END",
        "",
    ]))

    assert protein_chain_numbering(path, ["A"]) == [{
        "sequence": "AG",
        "numbering": [{"number": 10, "insertion": ""}, {"number": 10, "insertion": "A"}],
    }]


def test_compatibility_allows_different_prediction_engines_with_same_target_frame():
    left = _sample(engine="boltz2")
    right = _sample(engine="chai1", digest="different-coordinates")

    report = compatibility_report(left, right)

    assert report["compatible"] is True
    assert report["overlay_compatible"] is False
    assert not report["errors"]


def test_compatibility_blocks_mismatched_target_sequence_and_residue_numbering():
    left = _sample()
    right = _sample(engine="imported")
    right["target_chains_sequence"] = ["AT"]
    right["target_numbering"] = [{
        "sequence": "AT",
        "numbering": [{"number": 10, "insertion": ""}, {"number": 12, "insertion": ""}],
    }]

    report = compatibility_report(left, right)

    assert report["compatible"] is False
    assert any("different sequences" in reason for reason in report["errors"])
    assert any("numbering" in reason for reason in report["errors"])


def test_compatibility_requires_numbering_and_selected_component_groups():
    left = _sample()
    right = _sample()
    right["target_numbering"] = []
    right["binder_chains"] = []

    report = compatibility_report(left, right)

    assert report["compatible"] is False
    assert any("numbering could not be verified" in reason for reason in report["errors"])
    assert any("no selected protein partner" in reason for reason in report["errors"])


def test_compatibility_rejects_comparing_a_structure_representation_to_itself():
    left = {**_sample(), "identity": "db-molecule:structure-1"}
    right = {**_sample(), "identity": "db-molecule:structure-1"}

    report = compatibility_report(left, right)

    assert report["compatible"] is False
    assert any("two distinct saved structure representations" in reason for reason in report["errors"])


def test_contact_presence_table_only_transcribes_saved_target_residues():
    rows = saved_target_contact_rows({
        "Structure A": [{"target_residue": "A:10"}, {"target_residue": "A:11"}],
        "Structure B": [{"target_residue": "A:11"}],
    })

    assert rows == [
        {"Target residue": "A:10", "Structure A": True, "Structure B": False},
        {"Target residue": "A:11", "Structure A": True, "Structure B": True},
    ]


def test_comparison_view_supports_multiple_molstar_viewers_and_saved_evaluation_rows(tmp_path, monkeypatch):
    from mn_molecular_design_studio import complex_analysis_ui
    from mn_molecular_design_studio import complex_analysis_comparison_ui as view

    path_a = tmp_path / "a.pdb"
    path_b = tmp_path / "b.pdb"
    path_c = tmp_path / "c.pdb"
    for path in (path_a, path_b, path_c):
        path.write_text("END\n")
    target_numbering = [{"sequence": "AG", "numbering": [{"number": 10, "insertion": ""}, {"number": 11, "insertion": ""}]}]
    common = {
        "kind": "results", "label": "Saved structure", "path": path_a, "context": {},
        "components": [], "model_id": "model-1", "target_chains": ["A"], "binder_chains": ["B"],
        "target_numbering": target_numbering, "target_chains_sequence": ["AG"], "binder_chains_sequence": ["SS"],
        "complex_type": "protein-protein", "format": ".pdb", "viewer_key": "compare-a",
        "target_identity": "target-1", "coordinate_frame_id": "", "structure_digest": "a",
        "evaluation_entries": [{"evaluation_id": "ev-a", "evaluator": "boltz2", "run_id": "run-a",
                                "status": "completed", "metrics": {"iptm": 0.8}, "provenance": {"source": "Results"}}],
    }
    first = {**common, "identity": "a", "path": path_a}
    second = {**common, "identity": "b", "label": "Imported structure", "path": path_b,
              "viewer_key": "compare-b", "structure_digest": "b",
              "evaluation_entries": [{"evaluation_id": "ev-b", "evaluator": "chai1", "run_id": "run-b",
                                      "status": "completed", "metrics": {"iptm": 0.7}, "provenance": {"source": "Database"}}]}
    third = {**common, "identity": "c", "label": "Structure from another engine", "path": path_c,
             "viewer_key": "compare-c", "structure_digest": "c",
             "evaluation_entries": [{"evaluation_id": "ev-c", "evaluator": "protenix", "run_id": "run-c",
                                     "status": "completed", "metrics": {"iptm": 0.6}, "provenance": {"source": "Results"}}]}
    slots = {"a": first, "b": second, "c": third}
    shown = []
    tables = []
    monkeypatch.setattr(view, "_comparison_slot", lambda user, slot: slots[slot])
    monkeypatch.setattr(view, "_saved_biotite_contacts", lambda user, sample: None)
    monkeypatch.setattr(complex_analysis_ui, "_show_structure", lambda path, context, components, **kwargs: shown.append(path))
    monkeypatch.setattr(view.st, "dataframe", lambda rows, **kwargs: tables.append(rows))

    app = AppTest.from_string(
        "import streamlit as st\n"
        "st.session_state['complex-analysis-compare-count'] = 3\n"
        "from mn_molecular_design_studio.complex_analysis_comparison_ui import render_comparison\n"
        "render_comparison({'id': 'user-1'})",
        default_timeout=10,
    ).run()

    assert not app.exception
    assert shown == [path_a, path_b, path_c]
    assert any("Compatible target sequence" in item.value for item in app.success)
    evaluations = next(rows for rows in tables if rows and "Stored measurements" in rows[0])
    assert {row["Evaluator"] for row in evaluations} == {"boltz2", "chai1", "protenix"}
    assert {row["Evaluation / candidate ID"] for row in evaluations} == {"ev-a", "ev-b", "ev-c"}


def test_results_comparison_picker_uses_accessible_jobs_and_exact_saved_model(tmp_path, monkeypatch):
    from mn_molecular_design_studio import complex_analysis_ui
    from mn_molecular_design_studio import complex_analysis_comparison_ui as view

    path = tmp_path / "selected-model.pdb"
    path.write_text("\n".join([
        _atom(1, "ALA", "A", 10, 1.0),
        _atom(2, "SER", "B", 1, 4.0),
        "END",
        "",
    ]))
    trajectory_job = {"task_group": "simulation", "run_id": "trajectory-run", "tool": "bioemu2", "status": "completed"}
    accessible_job = {"task_group": "refolding-validation", "run_id": "accessible-run", "tool": "engine-a", "status": "completed"}
    candidate = {"candidate_id": "candidate-model-3", "model_index": 3,
                 "target_chains": ["A"], "binder_chains": ["B"],
                 "metrics": {"stored_metric": 0.6}, "raw_metadata": {}}
    selected = {"path": path, "relative_path": "models/selected-model.pdb", "candidate": candidate,
                "candidate_id": candidate["candidate_id"], "label": path.name}
    calls = []
    monkeypatch.setattr(complex_analysis_ui, "_accessible_result_jobs", lambda user: [trajectory_job, accessible_job])
    monkeypatch.setattr(complex_analysis_ui, "_job_structure_options", lambda user, job: calls.append(job) or [selected])

    script = '''\
import streamlit as st
from mn_molecular_design_studio.complex_analysis_comparison_ui import _comparison_slot
result = _comparison_slot({"id": "user-1"}, "a")
st.write(result["model_id"])
'''
    app = AppTest.from_string(script, default_timeout=10).run()

    assert not app.exception
    job_picker = next(widget for widget in app.selectbox if widget.label == "Completed Results job")
    assert job_picker.options == ["bioemu2 · simulation/trajectory-run", "engine-a · refolding-validation/accessible-run"]
    assert job_picker.value == "refolding-validation/accessible-run"
    assert calls == [accessible_job]
    assert any(item.value == "3" for item in app.markdown)


def _interface_complex(path: Path, *, target_x=1.0, partner_x=8.0) -> Path:
    path.write_text("\n".join([
        _atom(1, "ALA", "A", 10, target_x),
        _atom(2, "SER", "B", 1, partner_x),
        "END", "",
    ]))
    return path


def _interface_sample(path: Path, *, identity: str, run_id: str, candidate_id: str) -> dict:
    import hashlib

    return {
        "kind": "results", "identity": identity, "label": identity, "path": path,
        "context": {}, "components": [], "model_id": "model-1", "project_id": "project-1",
        "target_chains": ["A"], "binder_chains": ["B"],
        "target_numbering": [{"sequence": "A", "numbering": [{"number": 10, "insertion": ""}]}],
        "target_chains_sequence": ["A"], "binder_chains_sequence": ["S"],
        "complex_type": "protein-protein", "format": ".pdb", "viewer_key": identity,
        "structure_digest": hashlib.sha256(path.read_bytes()).hexdigest(),
        "candidate_id": candidate_id, "relative_path": f"models/{identity}.pdb",
        "job": {"task_group": "prediction", "run_id": run_id, "project_id": "project-1"},
        "evaluation_entries": [{"evaluation_id": candidate_id, "evaluator": "predictor",
                                "metrics": {"iptm": 0.7}}],
    }


def test_target_structure_fingerprint_ignores_partner_coordinates_but_tracks_target_coordinates(tmp_path):
    first = _interface_complex(tmp_path / "first.pdb", target_x=1.0, partner_x=8.0)
    second = _interface_complex(tmp_path / "second.pdb", target_x=1.0, partner_x=19.0)
    changed_target = _interface_complex(tmp_path / "changed-target.pdb", target_x=1.25, partner_x=8.0)

    first_digest = target_structure_fingerprint(first, ["A"])

    assert first_digest == target_structure_fingerprint(second, ["A"])
    assert first_digest != target_structure_fingerprint(changed_target, ["A"])


def test_compare_workspace_launches_registered_interface_clustering_for_exact_saved_sources(tmp_path, monkeypatch):
    from mn_molecular_design_studio import complex_analysis_ui
    from mn_molecular_design_studio import complex_analysis_comparison_ui as view

    first = _interface_sample(_interface_complex(tmp_path / "first.pdb"), identity="first", run_id="run-a", candidate_id="candidate-a")
    second = _interface_sample(_interface_complex(tmp_path / "second.pdb", partner_x=16.0), identity="second", run_id="run-b", candidate_id="candidate-b")
    slots = {"a": first, "b": second}
    monkeypatch.setattr(view, "_comparison_slot", lambda user, slot: slots[slot])
    monkeypatch.setattr(view, "_saved_biotite_contacts", lambda user, sample: None)
    monkeypatch.setattr(view, "_matching_interface_clustering_jobs", lambda *args, **kwargs: [])
    monkeypatch.setattr(complex_analysis_ui, "_show_structure", lambda *args, **kwargs: None)

    script = '''\
import streamlit as st
from mn_molecular_design_studio.complex_analysis_comparison_ui import render_comparison
def schema(workflow):
    assert workflow == "analysis.interface_clustering"
    return {"schema_version":1,"parameters":[
        {"name":"similarity_threshold","default":0.5},
        {"name":"contact_distance","default":5.0}]}
def submit(schema, workflow, parameters, project, settings, label):
    st.session_state["submitted"] = {"workflow":workflow,"parameters":parameters,"project":project,"label":label}
st.session_state["complex-analysis-compare-count"] = 2
render_comparison({"id":"user-1"}, workflow_schema=schema, submit=submit)
'''
    app = AppTest.from_string(script, default_timeout=15).run()
    assert not app.exception
    next(item for item in app.button if item.label == "Run Interface Clustering").click().run()

    assert not app.exception
    submitted = app.session_state["submitted"]
    assert submitted["workflow"] == "analysis.interface_clustering"
    assert submitted["project"] == "project-1"
    assert submitted["parameters"]["similarity_threshold"] == 0.5
    assert submitted["parameters"]["contact_distance"] == 5.0
    selection = submitted["parameters"]["interface_clustering_collection"]
    assert selection["schema"] == "studio.interface_clustering_selection.v1"
    assert [item["run_id"] for item in selection["sources"]] == ["run-a", "run-b"]
    assert [item["evaluation_ids"] for item in selection["sources"]] == [["candidate-a"], ["candidate-b"]]


def test_interface_clustering_results_show_saved_assignments_and_evaluation_ids(tmp_path, monkeypatch):
    import json
    from mn_molecular_design_studio import complex_analysis_ui
    from mn_molecular_design_studio.complex_analysis_comparison_ui import _render_interface_clustering_results

    run_dir = tmp_path / "analysis" / "cluster-run"
    artifacts = run_dir / "artifacts" / "analysis"
    artifacts.mkdir(parents=True)
    (run_dir / "metadata.json").write_text(json.dumps({"status": "completed"}))
    (artifacts / "batch_status.json").write_text(json.dumps({"total": 1, "finished": 1, "items": []}))
    (artifacts / "complex_analysis_collection.json").write_text(json.dumps({
        "collection_id": "collection-1", "selection_key": "selection-1",
        "items": [{"candidate_id": "comparison-item", "source_identity": "prediction/run-a:model.pdb"}],
    }))
    (artifacts / "qc_metrics.csv").write_text(
        "candidate_id,source_candidate_id,qc_status,diversity_cluster,diversity_cluster_size,diversity_representative\n"
        "comparison-item_interface_clustering,comparison-item,completed,1,1,True\n"
    )
    monkeypatch.setattr(complex_analysis_ui, "_analysis_evaluation_id", lambda *args, **kwargs: "evaluation-1")
    script = f'''\
from pathlib import Path
from mn_molecular_design_studio.complex_analysis_comparison_ui import _render_interface_clustering_results
_render_interface_clustering_results({{"id":"user-1"}}, {{"status":"completed","run_id":"cluster-run",
    "task_group":"analysis","run_dir":Path({str(run_dir)!r}),"manifest":{{"settings":{{"similarity_threshold":0.5}}}},
    "collection":{{"collection_id":"collection-1","selection_key":"selection-1","items":[
        {{"candidate_id":"comparison-item","source_identity":"prediction/run-a:model.pdb"}}]}}}}, key="test")
'''
    app = AppTest.from_string(script, default_timeout=10).run()

    assert not app.exception
    assert any("Cluster assignments saved" in item.value for item in app.markdown)
    rows = next(frame.value for frame in app.dataframe if frame.value is not None)
    assert rows.iloc[0]["Source structure"] == "prediction/run-a:model.pdb"
    assert rows.iloc[0]["Saved Evaluation"] == "evaluation-1"
    assert rows.iloc[0]["diversity_cluster"] == 1


def test_compare_workspace_launches_foldseek_for_nonmatching_targets(tmp_path, monkeypatch):
    from mn_molecular_design_studio import complex_analysis_ui
    from mn_molecular_design_studio import complex_analysis_comparison_ui as view

    first = _interface_sample(_interface_complex(tmp_path / "first.pdb"), identity="first", run_id="run-a", candidate_id="candidate-a")
    second = _interface_sample(_interface_complex(tmp_path / "second.pdb", partner_x=16.0), identity="second", run_id="run-b", candidate_id="candidate-b")
    # This selection cannot use target residue correspondence, but each
    # declared partner chain remains a valid independent Foldseek input.
    second["target_chains_sequence"] = ["DIFFERENT"]
    slots = {"a": first, "b": second}
    monkeypatch.setattr(view, "_comparison_slot", lambda user, slot: slots[slot])
    monkeypatch.setattr(view, "_saved_biotite_contacts", lambda user, sample: None)
    monkeypatch.setattr(view, "_matching_interface_clustering_jobs", lambda *args, **kwargs: [])
    monkeypatch.setattr(view, "_matching_foldseek_jobs", lambda *args, **kwargs: [])
    monkeypatch.setattr(complex_analysis_ui, "_show_structure", lambda *args, **kwargs: None)

    script = '''\
import streamlit as st
from mn_molecular_design_studio.complex_analysis_comparison_ui import render_comparison
def schema(workflow):
    parameters = {
        "analysis.interface_clustering": [
            {"name":"similarity_threshold","default":0.5},
            {"name":"contact_distance","default":5.0}],
        "analysis.foldseek": [
            {"name":"similarity_threshold","default":0.5},
            {"name":"coverage","default":0.8}],
    }
    return {"schema_version":1,"parameters":parameters[workflow]}
def submit(schema, workflow, parameters, project, settings, label):
    st.session_state["submitted"] = {"workflow":workflow,"parameters":parameters,"project":project,"label":label}
st.session_state["complex-analysis-compare-count"] = 2
render_comparison({"id":"user-1"}, workflow_schema=schema, submit=submit)
'''
    app = AppTest.from_string(script, default_timeout=15).run()
    assert not app.exception
    assert any("cannot share a target-residue comparison map" in item.value for item in app.warning)
    foldseek_button = next(item for item in app.button if item.label == "Run Foldseek")
    assert not foldseek_button.disabled
    foldseek_button.click().run()

    assert not app.exception
    submitted = app.session_state["submitted"]
    assert submitted["workflow"] == "analysis.foldseek"
    assert submitted["project"] == "project-1"
    assert submitted["parameters"]["similarity_threshold"] == 0.5
    assert submitted["parameters"]["coverage"] == 0.8
    selection = submitted["parameters"]["foldseek_collection"]
    assert selection["schema"] == "studio.foldseek_selection.v1"
    assert [item["run_id"] for item in selection["sources"]] == ["run-a", "run-b"]
    assert [item["binder_chains"] for item in selection["sources"]] == [["B"], ["B"]]


def test_foldseek_results_show_cluster_assignments_and_reused_evaluations(tmp_path, monkeypatch):
    import json
    from mn_molecular_design_studio import complex_analysis_ui
    from mn_molecular_design_studio.complex_analysis_comparison_ui import _render_foldseek_results

    run_dir = tmp_path / "analysis" / "foldseek-run"
    artifacts = run_dir / "artifacts" / "analysis"
    artifacts.mkdir(parents=True)
    (run_dir / "metadata.json").write_text(json.dumps({"status": "completed"}))
    (artifacts / "batch_status.json").write_text(json.dumps({"total": 1, "finished": 1, "items": []}))
    (artifacts / "qc_metrics.csv").write_text(
        "candidate_id,source_candidate_id,qc_status,diversity_cluster,diversity_cluster_size,diversity_representative\n"
        "comparison-item_foldseek,comparison-item,completed,2,3,True\n"
    )
    monkeypatch.setattr(complex_analysis_ui, "_analysis_evaluation_id", lambda *args, **kwargs: "foldseek-evaluation-1")
    script = f'''\
from pathlib import Path
from mn_molecular_design_studio.complex_analysis_comparison_ui import _render_foldseek_results
_render_foldseek_results({{"id":"user-1"}}, {{"status":"completed","run_id":"foldseek-run",
    "task_group":"analysis","run_dir":Path({str(run_dir)!r}),"manifest":{{"settings":{{"similarity_threshold":0.55,"coverage":0.75}}}},
    "collection":{{"collection_id":"collection-foldseek","selection_key":"selection-foldseek","items":[
        {{"candidate_id":"comparison-item","source_identity":"prediction/run-a:model.pdb"}}]}}}}, key="test")
'''
    app = AppTest.from_string(script, default_timeout=10).run()

    assert not app.exception
    assert any("Foldseek cluster assignments saved" in item.value for item in app.markdown)
    rows = next(frame.value for frame in app.dataframe if frame.value is not None)
    assert rows.iloc[0]["Source structure"] == "prediction/run-a:model.pdb"
    assert rows.iloc[0]["Saved Evaluation"] == "foldseek-evaluation-1"
    assert rows.iloc[0]["diversity_cluster"] == 2
