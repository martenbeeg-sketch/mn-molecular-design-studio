from __future__ import annotations

import json
import hashlib
from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest

from mn_molecular_design_studio import molecule_imports
from mn_molecular_design_studio.complex_analysis_ui import component_inventory
from mn_molecular_design_studio.molecular_database import MolecularDatabase


USER = {"id": "a" * 32, "role": "user", "_cluster_project_ids": []}
OTHER_USER = {"id": "b" * 32, "role": "user", "_cluster_project_ids": []}


def _atom(serial: int, record: str, atom: str, residue: str, chain: str, number: int, x: float) -> str:
    element = atom[0]
    return (
        f"{record:<6}{serial:5d} {atom:^4s} {residue:>3s} {chain}{number:4d}    "
        f"{x:8.3f}{1.0:8.3f}{1.0:8.3f}{1.0:6.2f}{20.0:6.2f}          {element:>2s}"
    )


def _complex_pdb(second_type: str, *, x: float = 2.0) -> bytes:
    lines = [_atom(1, "ATOM", "CA", "ALA", "A", 1, 1.0)]
    if second_type == "protein":
        lines.append(_atom(2, "ATOM", "CA", "GLY", "B", 1, x))
    elif second_type == "small_molecule":
        lines.append(_atom(2, "HETATM", "C1", "LIG", "L", 1, x))
    elif second_type == "dna":
        lines.append(_atom(2, "ATOM", "P", "DA", "D", 1, x))
    elif second_type == "rna":
        lines.append(_atom(2, "ATOM", "P", "A", "R", 1, x))
    else:
        raise ValueError(second_type)
    return ("\n".join([*lines, "END", ""])).encode()


def _two_protein_chain_pdb(first_sequence: str, second_sequence: str) -> str:
    residue_names = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE",
                     "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU",
                     "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG",
                     "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR"}
    lines, serial = [], 1
    for chain, sequence, offset in (("A", first_sequence, 0.0), ("B", second_sequence, 10.0)):
        for number, residue in enumerate(sequence, start=1):
            lines.append(
                f"ATOM  {serial:5d}  CA  {residue_names[residue]} {chain}{number:4d}    "
                f"{offset + number:8.3f}{1.0:8.3f}{1.0:8.3f}{1.00:6.2f}{20.00:6.2f}           C"
            )
            serial += 1
        lines.append("TER")
    return "\n".join([*lines, "END", ""])


def _isolated_database(tmp_path, monkeypatch) -> MolecularDatabase:
    workspace = tmp_path / "workspace"
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_HOME", str(workspace))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(workspace))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR", str(workspace / "runs"))
    return MolecularDatabase(workspace / "database")


@pytest.mark.parametrize(
    ("second_type", "expected"),
    [
        ("protein", "protein"),
        ("small_molecule", "small_molecule"),
        ("dna", "dna"),
        ("rna", "rna"),
    ],
)
def test_component_inventory_identifies_supported_complex_components(second_type, expected):
    components = component_inventory(_complex_pdb(second_type), "complex.pdb")

    assert {item["type"] for item in components} == {"protein", expected}
    assert {item["source_chain"] for item in components} == {
        "A", {"protein": "B", "small_molecule": "L", "dna": "D", "rna": "R"}[second_type]
    }


def test_structure_display_reuses_molstar_with_all_chains(tmp_path, monkeypatch):
    path = tmp_path / "complex.pdb"
    path.write_bytes(_complex_pdb("rna"))
    captured = []
    from mn_molecular_design_studio.engine.workflow_engine.app.components import molstar_viewer
    from mn_molecular_design_studio.complex_analysis_ui import _show_structure

    monkeypatch.setattr(molstar_viewer, "molstar_custom_component", lambda structures, **kwargs: captured.append((structures, kwargs)))
    script = f'''\
from pathlib import Path
from mn_molecular_design_studio.complex_analysis_ui import _show_structure
_show_structure(Path({str(path)!r}), {{}}, [], key="complex-test", highlighted_selections=["A:1", "R:1"])
'''
    app = AppTest.from_string(script, default_timeout=10).run()

    assert not app.exception
    assert len(captured) == 1
    visualization = captured[0][0][0]
    assert " ALA A" in visualization.pdb and " A R" in visualization.pdb
    assert visualization.color == "chain-id"
    assert visualization.representation_type == "cartoon+ball-and-stick"
    assert visualization.highlighted_selections == ["A:1", "R:1"]


def test_inspector_finding_selection_reaches_molstar_with_insertion_code(tmp_path, monkeypatch):
    path = tmp_path / "complex.pdb"
    path.write_text("\n".join([
        _atom(1, "ATOM", "CA", "GLY", "A", 10, 1.0),
        _atom(2, "ATOM", "CA", "GLY", "A", 10, 2.0),
        "END", "",
    ]))
    # Encode an author insertion code in the saved fixture for the exact Mol* selector.
    lines = path.read_text().splitlines()
    lines[1] = lines[1][:26] + "A" + lines[1][27:]
    path.write_text("\n".join([*lines, ""]))
    captured = []
    from mn_molecular_design_studio.engine.workflow_engine.app.components import molstar_viewer
    from mn_molecular_design_studio.complex_analysis_ui import (
        _active_biotite_highlights, _render_inspector_findings, _show_structure,
    )

    monkeypatch.setattr(molstar_viewer, "molstar_custom_component",
                        lambda structures, **kwargs: captured.append(structures[0].highlighted_selections))
    data = {
        "evidence": [{"id": "finding-1", "evaluator": "plip", "label": "hydrogen bond · A:10A",
                      "artifact": "interactions.csv", "residue_refs": ["A:10^A"],
                      "selectors": ["A:10^A", "L:1^_"], "raw": {"distance_angstrom": "2.7"}}],
        "evaluations": [{"evaluator": "plip", "evaluation": {"id": "evaluation-1"}}],
    }
    script = f'''\
from pathlib import Path
from mn_molecular_design_studio.complex_analysis_ui import _active_biotite_highlights, _render_inspector_findings, _show_structure
key = "inspector-test"
_show_structure(Path({str(path)!r}), {{}}, [], key=key, highlighted_selections=_active_biotite_highlights(key))
_render_inspector_findings({data!r}, key)
'''
    app = AppTest.from_string(script, default_timeout=10).run()
    assert not app.exception
    highlight = next(button for button in app.button if button.label == "Highlight evidence in Mol*")

    app = highlight.click().run()

    assert not app.exception
    assert captured[-1] == ["A:10^A", "L:1^_"]


def test_inspector_sequence_track_selection_reaches_molstar_with_insertion_code(tmp_path, monkeypatch):
    path = tmp_path / "complex.pdb"
    path.write_text("\n".join([
        _atom(1, "ATOM", "CA", "GLY", "A", 10, 1.0),
        _atom(2, "ATOM", "CA", "ALA", "A", 10, 2.0),
        "END", "",
    ]))
    lines = path.read_text().splitlines()
    lines[1] = lines[1][:26] + "A" + lines[1][27:]
    path.write_text("\n".join([*lines, ""]))
    captured = []
    from mn_molecular_design_studio.engine.workflow_engine.app.components import molstar_viewer
    from mn_molecular_design_studio.complex_analysis_ui import _active_biotite_highlights, _render_inspector_tracks, _show_structure

    monkeypatch.setattr(molstar_viewer, "molstar_custom_component",
                        lambda structures, **kwargs: captured.append(structures[0].highlighted_selections))
    tracks = [{
        "chain_id": "A", "sequence": "GA",
        "residues": [
            {"sequence_index": 0, "amino_acid": "G", "author_residue": "10", "selector": "A:10^_"},
            {"sequence_index": 1, "amino_acid": "A", "author_residue": "10A", "selector": "A:10^A"},
        ],
    }]
    script = f'''\
from pathlib import Path
from mn_molecular_design_studio.complex_analysis_ui import _active_biotite_highlights, _render_inspector_tracks, _show_structure
key = "inspector-track-test"
data = {{"tracks": {tracks!r}}}
_show_structure(Path({str(path)!r}), {{}}, [], key=key, highlighted_selections=_active_biotite_highlights(key))
_render_inspector_tracks(data, key)
'''
    app = AppTest.from_string(script, default_timeout=10).run()

    assert not app.exception
    app = next(button for button in app.button if button.label == "A\n10A").click().run()

    assert not app.exception
    assert captured[-1] == ["A:10^A"]


def test_biotite_contact_residues_keep_saved_auth_chain_numbering(tmp_path):
    run_dir = tmp_path / "analysis-run"
    item_dir = run_dir / "artifacts" / "analysis" / "item-1"
    item_dir.mkdir(parents=True)
    (run_dir / "artifacts" / "analysis" / "batch_manifest.json").write_text(json.dumps({
        "items": [
            {"id": "item-1", "candidate_id": "candidate-model-2"},
            {"id": "item-other", "candidate_id": "candidate-model-1"},
        ],
    }))
    (item_dir / "contacts.json").write_text(json.dumps({"residue_pairs": [
        {"binder_residue": "B:42", "target_residue": "A:8"},
        {"binder_residue": "B:42", "target_residue": "A:9"},
        {"binder_residue": "B:42A", "target_residue": "A:10"},
    ]}))

    from mn_molecular_design_studio.complex_analysis_ui import _biotite_contact_residues
    binder, target, combined = _biotite_contact_residues(run_dir, "candidate-model-2")

    assert binder == ["B:42"]
    assert target == ["A:10", "A:8", "A:9"]
    assert combined == ["A:10", "A:8", "A:9", "B:42"]


def test_completed_results_evaluation_reopens_with_same_structure_viewer_context(tmp_path, monkeypatch):
    from mn_molecular_design_studio import complex_analysis_ui as workspace

    run_dir = tmp_path / "source-run"
    (run_dir / "artifacts" / "normalized_candidates").mkdir(parents=True)
    (run_dir / "artifacts" / "normalized_candidates" / "candidates.jsonl").write_text("{}\n")
    path = tmp_path / "complex.pdb"
    path.write_bytes(_complex_pdb("protein"))
    job = {"task_group": "refolding-validation", "run_id": "source-1", "run_dir": str(run_dir),
           "project_id": "project-1", "status": "completed", "tool": "boltz2"}
    candidate = {"candidate_id": "candidate-1", "model_index": "model-3",
                 "binder_chains": ["B"], "target_chains": ["A"], "raw_metadata": {}}
    selected = {"path": path, "relative_path": "models/complex.pdb", "artifact_path": str(path),
                "candidate_id": "candidate-1", "label": "complex model", "candidate": candidate}
    current = {"status": "completed", "run_dir": run_dir, "manifest": {}, "metadata": {},
               "selected_item": {"candidate_id": "candidate-1"}, "run_id": "analysis-1"}
    captured = {}

    monkeypatch.setattr(workspace, "_accessible_result_jobs", lambda user: [job])
    monkeypatch.setattr(workspace, "_job_structure_options", lambda user, row: [selected])
    monkeypatch.setattr(workspace, "_show_structure", lambda *args, **kwargs: None)
    monkeypatch.setattr(workspace, "_show_component_groups", lambda *args, **kwargs: None)
    monkeypatch.setattr(workspace, "_show_human_review", lambda *args, **kwargs: None)
    monkeypatch.setattr(workspace, "_matching_biotite_jobs", lambda *args, **kwargs: [current])
    monkeypatch.setattr(workspace, "_analysis_evaluation_id", lambda *args, **kwargs: "evaluation-saved")
    monkeypatch.setattr(workspace, "_render_biotite_results", lambda *args, **kwargs: captured.update(kwargs))

    script = '''\
from mn_molecular_design_studio.complex_analysis_ui import _render_results_source
_render_results_source({"id": "user-1"}, {},
    workflow_schema=lambda workflow_id: {"parameters": [{"name": "contact_distance", "default": 5.0}]},
    project_selector=lambda *args, **kwargs: "project-1", submit=None)
'''
    app = AppTest.from_string(script, default_timeout=10).run()

    assert not app.exception
    assert captured["viewer_key"] == "complex-analysis-results-viewer-refolding-validation/source-1-models/complex.pdb"
    assert captured["structure_path"] == path


def test_database_workspace_selects_exact_model_evaluation_and_shows_imported_structure(tmp_path, monkeypatch):
    database = _isolated_database(tmp_path, monkeypatch)
    model_a = _complex_pdb("dna", x=2.0)
    model_b = _complex_pdb("dna", x=9.0)
    imported_record = molecule_imports.parse_structure(model_a, "model-a.pdb")[0]
    imported_record["name"] = "Imported protein–DNA complex"
    imported_record["files"][0].update({"model_index": 1, "run_id": "prediction-1", "candidate_id": "complex-1"})
    imported = database.import_records(USER, [imported_record], name="Imported structures")
    molecule_id = imported["molecule_ids"][0]
    initial = database.list_molecules(USER, molecule_ids=[molecule_id])["items"][0]["representations"][0]
    model_b_id = database.attach_structure(
        imported["group_id"], molecule_id,
        {"name": "model-b.pdb", "content": model_b, "label": "Model 2", "model_index": 2,
         "run_id": "prediction-1", "candidate_id": "complex-1"},
    )
    observation_id = database.add_observation(
        imported["group_id"], molecule_id, "boltz2", "prediction-1", {"iptm": 0.87, "pae": 3.2},
        context={"representation_id": model_b_id, "target_id": "dna-target", "candidate_id": "complex-1",
                 "candidate": {"raw_metadata": {"confidence_outputs": {"ptm": 0.93}}}},
    )
    from mn_molecular_design_studio.evaluations import evaluation_descriptor

    record = database.list_molecules(USER, molecule_ids=[molecule_id])["items"][0]
    observation = next(row for row in record["observations"] if row["id"] == observation_id)
    evaluation_id = evaluation_descriptor(observation)["id"]
    shown = []

    from mn_molecular_design_studio import complex_analysis_ui, molecular_database

    monkeypatch.setattr(molecular_database, "MolecularDatabase", lambda: database)
    monkeypatch.setattr(complex_analysis_ui, "_show_structure", lambda path, context, components, **kwargs: shown.append(path.read_bytes()))
    script = f'''\
import streamlit as st
from mn_molecular_design_studio.complex_analysis_ui import render_complex_analysis_workspace
st.session_state["complex-analysis-source-kind"] = "Database"
st.session_state["complex-analysis-open"] = {{"kind": "database", "molecule_id": {molecule_id!r}, "evaluation_id": {evaluation_id!r}, "structure_id": {model_b_id!r}}}
render_complex_analysis_workspace({USER!r})
'''
    app = AppTest.from_string(script, default_timeout=10).run()

    assert not app.exception
    assert len(shown) == 1
    assert shown[0] == model_b
    assert any("No saved evaluation is associated" in item.value for item in app.info) is False
    saved_json = [item.value for item in app.json]
    assert any("0.87" in json.dumps(value) for value in saved_json)
    assert any("0.93" in json.dumps(value) for value in saved_json)
    assert initial["id"] != model_b_id

    app.multiselect(key=f"complex-analysis-group-members-{model_b_id}-group-1").set_value(["dna:D"]).run()
    assert not app.exception
    assert any("Group 1: dna chain D" in item.value for item in app.caption)


def test_results_workspace_displays_generic_predicted_complex_and_saved_confidence(tmp_path, monkeypatch):
    database = _isolated_database(tmp_path, monkeypatch)
    structure_path = tmp_path / "predicted-complex.pdb"
    structure_path.write_bytes(_complex_pdb("small_molecule"))
    candidate = {
        "candidate_id": "prediction-1",
        "model_index": 1,
        "stage": "complex_refolding",
        "source_tool": "boltz2",
        "metrics": {"iptm": 0.76},
        "raw_metadata": {"confidence_outputs": {"ptm": 0.81, "plddt": 88.4}},
    }
    job = {"task_group": "boltz2", "run_id": "predicted-1", "status": "completed", "tool": "boltz2"}
    option = {
        "path": structure_path,
        "relative_path": "predicted-complex.pdb",
        "candidate": candidate,
        "candidate_id": candidate["candidate_id"],
        "label": structure_path.name,
    }
    shown = []

    from mn_molecular_design_studio import complex_analysis_ui, molecular_database

    monkeypatch.setattr(molecular_database, "MolecularDatabase", lambda: database)
    monkeypatch.setattr(complex_analysis_ui, "_accessible_result_jobs", lambda user: [job])
    monkeypatch.setattr(complex_analysis_ui, "_job_structure_options", lambda user, selected_job: [option])
    monkeypatch.setattr(complex_analysis_ui, "_show_structure", lambda path, context, components, **kwargs: shown.append(path))
    script = f'''\
import streamlit as st
from mn_molecular_design_studio.complex_analysis_ui import render_complex_analysis_workspace
st.session_state["complex-analysis-source-kind"] = "Results"
render_complex_analysis_workspace({USER!r})
'''
    app = AppTest.from_string(script, default_timeout=10).run()

    assert not app.exception
    assert shown == [structure_path]
    assert any("0.76" in json.dumps(item.value) and "0.81" in json.dumps(item.value) for item in app.json)


def test_results_workspace_submits_exact_biotite_model_and_component_groups(tmp_path, monkeypatch):
    database = _isolated_database(tmp_path, monkeypatch)
    from mn_molecular_design_studio import complex_analysis_ui, molecular_database

    source_run = tmp_path / "workspace" / "runs" / "prediction" / "source-job"
    source_run.mkdir(parents=True)
    structure_path = source_run / "model-2.pdb"
    structure_path.write_bytes(_complex_pdb("protein"))
    candidate = {
        "candidate_id": "generic-complex-42", "model_index": 2,
        "stage": "complex_refolding", "source_tool": "generic-predictor",
        "complex_pdb": "model-2.pdb", "binder_chains": ["B"], "target_chains": ["A"],
        "metrics": {"confidence": 0.71}, "raw_metadata": {"model_index": 2},
    }
    candidate_file = source_run / "artifacts" / "normalized_candidates" / "candidates.jsonl"
    candidate_file.parent.mkdir(parents=True)
    candidate_file.write_text(json.dumps(candidate) + "\n")
    (source_run / "metadata.json").write_text(json.dumps({"tool": "generic-predictor", "status": "completed"}))
    (source_run / "input.json").write_text(json.dumps({"tool": "generic-predictor", "params": {}}))
    (source_run / "result.json").write_text(json.dumps({"success": True, "metrics": {}}))
    (source_run / "molecular_design_studio.json").write_text(json.dumps({"studio_owner_id": USER["id"], "project_id": "test-project"}))
    job = {"task_group": "prediction", "run_id": "source-job", "status": "completed",
           "tool": "generic-predictor", "run_dir": str(source_run), "project_id": "test-project"}
    option = {"path": structure_path, "relative_path": "model-2.pdb", "candidate": candidate,
              "candidate_id": candidate["candidate_id"], "label": "Model 2"}

    monkeypatch.setattr(molecular_database, "MolecularDatabase", lambda: database)
    monkeypatch.setattr(complex_analysis_ui, "_accessible_result_jobs", lambda user: [job])
    monkeypatch.setattr(complex_analysis_ui, "_job_structure_options", lambda user, selected_job: [option])
    monkeypatch.setattr(complex_analysis_ui, "_show_structure", lambda *args, **kwargs: None)
    monkeypatch.setattr(complex_analysis_ui, "_show_human_review", lambda *args, **kwargs: None)
    monkeypatch.setattr(complex_analysis_ui, "_matching_biotite_jobs", lambda *args, **kwargs: [])
    submitted = []

    def submit(schema, workflow, parameters, project, settings, title):
        submitted.append((schema, workflow, dict(parameters), project, settings, title))
        import streamlit as st
        st.session_state["last_submitted_job"] = {"run_id": "analysis-1"}

    key = hashlib.sha256(b"prediction/source-job:model-2.pdb").hexdigest()[:14]
    script = f'''\
import streamlit as st
from mn_molecular_design_studio.complex_analysis_ui import render_complex_analysis_workspace
st.session_state["complex-analysis-source-kind"] = "Results"
render_complex_analysis_workspace(
    {USER!r},
    workflow_schema=lambda workflow_id: {{"workflow": workflow_id, "schema_version": 1, "parameters": [{{"name": "contact_distance", "default": 5.0}}]}},
    project_selector=lambda prefix, **kwargs: "test-project",
    submit={submit!r},
)
'''
    # AppTest executes the submitted callback in its own script namespace.
    script = script.replace(f"submit={submit!r}", "submit=submit")
    script = script.replace("render_complex_analysis_workspace(\n", "def submit(schema, workflow, parameters, project, settings, title):\n    st.session_state['submitted_parameters'] = parameters\n    st.session_state['submitted_workflow'] = workflow\n    st.session_state['last_submitted_job'] = {'run_id': 'analysis-1'}\nrender_complex_analysis_workspace(\n")
    app = AppTest.from_string(script, default_timeout=10).run()

    assert not app.exception
    run_button = app.button(key=f"complex-analysis-run-biotite-{key}")
    assert run_button is not None and not run_button.disabled
    run_button.click().run()
    assert not app.exception
    assert app.session_state["submitted_workflow"] == "analysis.biotite_interface"
    parameters = app.session_state["submitted_parameters"]
    assert parameters == {
        "candidate_ids": ["generic-complex-42"],
        "selected_structure_path": str(structure_path.resolve()),
        "selected_model_id": "2",
        "binder_chains": ["B"], "target_chains": ["A"], "contact_distance": 5.0,
        "source_job": "prediction/source-job",
    }


def test_biotite_matching_reuses_only_accessible_exact_analysis_jobs(monkeypatch, tmp_path):
    from mn_molecular_design_studio import complex_analysis_ui
    from mn_molecular_design_studio import auth, job_catalog

    source = {"task_group": "prediction", "run_id": "source", "candidate_id": "candidate"}
    accessible_run = tmp_path / "analysis" / "accessible"
    private_run = tmp_path / "analysis" / "private"
    for run_dir, owner in ((accessible_run, USER["id"]), (private_run, OTHER_USER["id"])):
        (run_dir / "artifacts" / "analysis").mkdir(parents=True)
        (run_dir / "metadata.json").write_text(json.dumps({
            "upstream_task_group": "prediction", "upstream_run_id": "source", "status": "completed",
        }))
        manifest = {
            "tool": "biotite_interface", "settings": {"contact_distance": 5.0},
            "items": [{"id": "item_000000", "candidate_id": "candidate", "structure_digest": "sha256-structure",
                       "binder_chains": ["B"], "target_chains": ["A"]}],
        }
        (run_dir / "artifacts" / "analysis" / "batch_manifest.json").write_text(json.dumps(manifest))
        (run_dir / "artifacts" / "analysis" / "batch_status.json").write_text(json.dumps({
            "items": [{"id": "item_000000", "status": "failed", "error": "Selected chains were absent."}],
        }))
    rows = [
        {"task_group": "analysis", "run_id": "accessible", "run_dir": str(accessible_run), "status": "completed", "owner_id": USER["id"]},
        {"task_group": "analysis", "run_id": "private", "run_dir": str(private_run), "status": "completed", "owner_id": OTHER_USER["id"]},
    ]
    monkeypatch.setattr(job_catalog, "search_jobs", lambda **kwargs: {"items": rows})
    monkeypatch.setattr(job_catalog, "get_job_access_record", lambda group, run: {"owner_id": USER["id"] if run == "accessible" else OTHER_USER["id"]})
    monkeypatch.setattr(auth, "can_access_job", lambda user, record: record.get("owner_id") == user["id"])
    monkeypatch.setattr(auth, "accessible_project_ids", lambda user: [])

    matches = complex_analysis_ui._matching_biotite_jobs(
        USER, source_task_group=source["task_group"], source_run_id=source["run_id"],
        candidate_id=source["candidate_id"], structure_digest="sha256-structure",
        binder_chains=["B"], target_chains=["A"], contact_distance=5.0,
    )

    assert [item["run_id"] for item in matches] == ["accessible"]
    assert matches[0]["manifest"]["tool"] == "biotite_interface"
    assert matches[0]["status"] == "failed"
    assert matches[0]["error"] == "Selected chains were absent."


def test_database_workspace_handles_missing_metadata_and_denies_private_structure(tmp_path, monkeypatch):
    database = _isolated_database(tmp_path, monkeypatch)
    imported_record = molecule_imports.parse_structure(_complex_pdb("protein"), "private-complex.pdb")[0]
    imported = database.import_records(OTHER_USER, [imported_record], name="Private PPI")
    row = database.list_molecules(OTHER_USER, molecule_ids=imported["molecule_ids"])["items"][0]
    structure_id = row["representations"][0]["id"]
    from mn_molecular_design_studio import complex_analysis_ui, molecular_database

    monkeypatch.setattr(molecular_database, "MolecularDatabase", lambda: database)
    shown = []
    monkeypatch.setattr(complex_analysis_ui, "_show_structure", lambda *args, **kwargs: shown.append(args[0]))
    script = f'''\
import streamlit as st
from mn_molecular_design_studio.complex_analysis_ui import render_complex_analysis_workspace
st.session_state["complex-analysis-source-kind"] = "Database"
st.session_state["complex-analysis-open"] = {{"kind": "database", "molecule_id": {row['id']!r}}}
render_complex_analysis_workspace({USER!r})
'''
    app = AppTest.from_string(script, default_timeout=10).run()

    assert not app.exception
    assert not shown
    assert any("No accessible Database records" in item.value for item in app.info)
    owner_script = f'''\
import streamlit as st
from mn_molecular_design_studio.complex_analysis_ui import render_complex_analysis_workspace
st.session_state["complex-analysis-source-kind"] = "Database"
st.session_state["complex-analysis-open"] = {{"kind": "database", "molecule_id": {row['id']!r}}}
render_complex_analysis_workspace({OTHER_USER!r})
'''
    owner_app = AppTest.from_string(owner_script, default_timeout=10).run()
    assert not owner_app.exception
    assert any("No saved evaluation is associated" in item.value for item in owner_app.info)
    with pytest.raises(PermissionError):
        database.structure_file(USER, structure_id)


def test_completed_results_structures_follow_existing_job_permissions(tmp_path, monkeypatch):
    database = _isolated_database(tmp_path, monkeypatch)
    run_dir = tmp_path / "workspace" / "runs" / "boltz2" / "private-result"
    run_dir.mkdir(parents=True)
    (run_dir / "molecular_design_studio.json").write_text(json.dumps({"studio_owner_id": OTHER_USER["id"], "project_id": ""}))
    (run_dir / "metadata.json").write_text(json.dumps({"tool": "boltz2", "status": "completed", "campaign_name": "Private prediction"}))
    (run_dir / "input.json").write_text(json.dumps({"tool": "boltz2", "params": {}}))
    (run_dir / "result.json").write_text(json.dumps({"success": True, "metrics": {}}))
    structure_path = run_dir / "complex.pdb"
    structure_path.write_bytes(_complex_pdb("protein"))
    from mn_molecular_design_studio import complex_analysis_ui
    from mn_molecular_design_studio.job_catalog import reconcile

    reconcile()
    assert complex_analysis_ui._accessible_result_jobs(USER) == []
    visible = complex_analysis_ui._accessible_result_jobs(OTHER_USER)
    assert len(visible) == 1
    assert complex_analysis_ui._job_structure_options(USER, visible[0]) == []
    options = complex_analysis_ui._job_structure_options(OTHER_USER, visible[0])
    assert [item["path"] for item in options] == [structure_path]
    saved = database.save_complex_review(
        OTHER_USER, source_kind="results", subject_id="boltz2/private-result",
        structure_ref="complex.pdb", status="Reviewed", note="Manual review note.",
    )
    assert saved["note"] == "Manual review note."
    with pytest.raises(PermissionError):
        database.get_complex_review(
            USER, source_kind="results", subject_id="boltz2/private-result", structure_ref="complex.pdb",
        )


def test_database_review_annotations_persist_separately_and_follow_record_access(tmp_path, monkeypatch):
    database = _isolated_database(tmp_path, monkeypatch)
    imported_record = molecule_imports.parse_structure(_complex_pdb("dna"), "reviewed-complex.pdb")[0]
    imported = database.import_records(USER, [imported_record], name="Reviewed complex")
    record = database.list_molecules(USER, molecule_ids=imported["molecule_ids"])["items"][0]
    structure_id = record["representations"][0]["id"]
    before_observations = list(record["observations"])

    saved = database.save_complex_review(
        USER, source_kind="database", subject_id=record["id"], structure_ref=structure_id,
        status="Needs review", note="Check source numbering before interpreting this model.",
    )
    reloaded = MolecularDatabase(database.root).get_complex_review(
        USER, source_kind="database", subject_id=record["id"], structure_ref=structure_id,
    )

    assert saved["status"] == reloaded["status"] == "Needs review"
    assert saved["note"] == reloaded["note"]
    assert reloaded["updated_at"]
    assert database.list_molecules(USER, molecule_ids=[record["id"]])["items"][0]["observations"] == before_observations
    with pytest.raises(PermissionError):
        database.get_complex_review(
            OTHER_USER, source_kind="database", subject_id=record["id"], structure_ref=structure_id,
        )
    with pytest.raises(ValueError):
        database.save_complex_review(
            USER, source_kind="database", subject_id=record["id"], structure_ref=structure_id,
            status="Predicted binder", note="Not a supported review status.",
        )


@pytest.mark.parametrize(
    ("workflow_id", "button_suffix"),
    [
        ("analysis.rosetta", "analysis-rosetta"),
        ("analysis.ipsae", "analysis-ipsae"),
        ("analysis.pep_patch", "analysis-pep_patch"),
    ],
)
def test_contextual_workflows_submit_selected_pose_through_existing_launcher(
    tmp_path, monkeypatch, workflow_id, button_suffix,
):
    database = _isolated_database(tmp_path, monkeypatch)
    from mn_molecular_design_studio import complex_analysis_ui, molecular_database

    source_run = tmp_path / "workspace" / "runs" / "prediction" / "source-job"
    source_run.mkdir(parents=True)
    structure_path = source_run / "model-2.pdb"
    structure_path.write_bytes(_complex_pdb("protein"))
    candidate = {
        "candidate_id": "generic-complex-42", "model_index": 2,
        "stage": "complex_refolding", "source_tool": "generic-predictor",
        "complex_pdb": "model-2.pdb", "binder_chains": ["B"], "target_chains": ["A"],
        "raw_metadata": {"model_index": 2},
    }
    candidate_file = source_run / "artifacts" / "normalized_candidates" / "candidates.jsonl"
    candidate_file.parent.mkdir(parents=True)
    candidate_file.write_text(json.dumps(candidate) + "\n")
    (source_run / "metadata.json").write_text(json.dumps({"tool": "generic-predictor", "status": "completed"}))
    (source_run / "input.json").write_text(json.dumps({"tool": "generic-predictor", "params": {}}))
    (source_run / "result.json").write_text(json.dumps({"success": True, "metrics": {}}))
    (source_run / "molecular_design_studio.json").write_text(json.dumps({"studio_owner_id": USER["id"], "project_id": "test-project"}))
    job = {"task_group": "prediction", "run_id": "source-job", "status": "completed",
           "tool": "generic-predictor", "run_dir": str(source_run), "project_id": "test-project"}
    option = {"path": structure_path, "relative_path": "model-2.pdb", "candidate": candidate,
              "candidate_id": candidate["candidate_id"], "label": "Model 2"}

    monkeypatch.setattr(molecular_database, "MolecularDatabase", lambda: database)
    monkeypatch.setattr(complex_analysis_ui, "_accessible_result_jobs", lambda user: [job])
    monkeypatch.setattr(complex_analysis_ui, "_job_structure_options", lambda user, selected_job: [option])
    monkeypatch.setattr(complex_analysis_ui, "_show_structure", lambda *args, **kwargs: None)
    monkeypatch.setattr(complex_analysis_ui, "_show_human_review", lambda *args, **kwargs: None)
    monkeypatch.setattr(complex_analysis_ui, "_matching_biotite_jobs", lambda *args, **kwargs: [])
    monkeypatch.setattr(complex_analysis_ui, "_matching_contextual_jobs", lambda *args, **kwargs: [])

    key = hashlib.sha256(b"prediction/source-job:model-2.pdb").hexdigest()[:14]
    script = f'''\
import streamlit as st
from mn_molecular_design_studio.complex_analysis_ui import render_complex_analysis_workspace
def schema_provider(workflow_id):
    names = {{
        "analysis.biotite_interface": [{{"name":"contact_distance", "default":5.0}}],
        "analysis.rosetta": [{{"name":"nprocs", "default":4}}],
        "analysis.ipsae": [{{"name":"pae_cutoff", "default":10.0}}, {{"name":"dist_cutoff", "default":10.0}},
            {{"name":"ipdae_contact_cutoff", "default":8.0}}, {{"name":"structure_contact_cutoff", "default":4.0}},
            {{"name":"hotspot_contact_cutoff", "default":5.0}}],
        "analysis.pep_patch": [{{"name":"ph", "default":7.4}}],
    }}
    return {{"workflow":workflow_id,"schema_version":1,"parameters":names[workflow_id]}}
def submit(schema, workflow, parameters, project, settings, title):
    st.session_state["submitted_workflow"] = workflow
    st.session_state["submitted_parameters"] = parameters
    st.session_state["submitted_project"] = project
    st.session_state["last_submitted_job"] = {{"run_id":"analysis-1"}}
st.session_state["complex-analysis-source-kind"] = "Results"
render_complex_analysis_workspace(
    {USER!r}, workflow_schema=schema_provider,
    project_selector=lambda prefix, **kwargs: "test-project", submit=submit,
)
'''
    app = AppTest.from_string(script, default_timeout=15).run()
    assert not app.exception
    button = app.button(key=f"complex-analysis-run-{button_suffix}-{key}")
    assert button is not None and not button.disabled
    button.click().run()

    assert not app.exception
    assert app.session_state["submitted_workflow"] == workflow_id
    assert app.session_state["submitted_project"] == "test-project"
    parameters = app.session_state["submitted_parameters"]
    assert parameters["candidate_ids"] == [candidate["candidate_id"]]
    assert parameters["selected_structure_path"] == str(structure_path.resolve())
    assert parameters["selected_model_id"] == "2"
    assert parameters["source_job"] == "prediction/source-job"
    if workflow_id == "analysis.pep_patch":
        assert parameters["binder_chains"] == ["B"]
        assert parameters["ph"] == 7.4
        assert "target_chains" not in parameters
    else:
        assert parameters["binder_chains"] == ["B"]
        assert parameters["target_chains"] == ["A"]


@pytest.mark.parametrize("workflow_id", ["analysis.plip", "analysis.pandamap", "analysis.posebusters"])
def test_contextual_ligand_workflows_submit_exact_residue_without_interface_groups(
    tmp_path, monkeypatch, workflow_id,
):
    database = _isolated_database(tmp_path, monkeypatch)
    from mn_molecular_design_studio import complex_analysis_ui, molecular_database
    source_run = tmp_path / "workspace" / "runs" / "prediction" / "ligand-source"
    source_run.mkdir(parents=True)
    structure_path = source_run / "model-2.pdb"
    structure_path.write_bytes(_complex_pdb("small_molecule"))
    smiles = "CCO"
    ligand_component = {"id": "ligand-residue", "type": "small_molecule", "chain": "L_1",
                        "source_chain": "L", "residue_key": "L:1", "ccd": "LIG", "label": "Ethanol"}
    candidate = {
        "candidate_id": "saved-ligand-42", "model_index": 2, "complex_pdb": "model-2.pdb",
        "raw_metadata": {"components": [{**ligand_component, "smiles": smiles}],
                         "ligand_smiles": smiles},
    }
    candidate_file = source_run / "artifacts" / "normalized_candidates" / "candidates.jsonl"
    candidate_file.parent.mkdir(parents=True)
    candidate_file.write_text(json.dumps(candidate) + "\n")
    (source_run / "metadata.json").write_text(json.dumps({"tool": "ligand-predictor", "status": "completed"}))
    (source_run / "input.json").write_text(json.dumps({"tool": "ligand-predictor", "params": {}}))
    (source_run / "result.json").write_text(json.dumps({"success": True, "metrics": {}}))
    (source_run / "molecular_design_studio.json").write_text(json.dumps({
        "studio_owner_id": USER["id"], "project_id": "test-project",
    }))
    job = {"task_group": "prediction", "run_id": "ligand-source", "status": "completed",
           "tool": "ligand-predictor", "run_dir": str(source_run), "project_id": "test-project"}
    option = {"path": structure_path, "relative_path": "model-2.pdb", "candidate": candidate,
              "candidate_id": candidate["candidate_id"], "label": "Model 2"}
    components = [
        {"id": "protein-A", "type": "protein", "source_chain": "A", "chain": "A", "label": "Protein A", "residue_key": "A:1"},
        ligand_component,
    ]
    monkeypatch.setattr(molecular_database, "MolecularDatabase", lambda: database)
    monkeypatch.setattr(complex_analysis_ui, "_accessible_result_jobs", lambda user: [job])
    monkeypatch.setattr(complex_analysis_ui, "_job_structure_options", lambda user, selected_job: [option])
    monkeypatch.setattr(complex_analysis_ui, "component_inventory", lambda *args, **kwargs: components)
    monkeypatch.setattr(complex_analysis_ui, "_show_structure", lambda *args, **kwargs: None)
    monkeypatch.setattr(complex_analysis_ui, "_show_human_review", lambda *args, **kwargs: None)
    monkeypatch.setattr(complex_analysis_ui, "_matching_biotite_jobs", lambda *args, **kwargs: [])
    monkeypatch.setattr(complex_analysis_ui, "_matching_contextual_jobs", lambda *args, **kwargs: [])
    key = hashlib.sha256(b"prediction/ligand-source:model-2.pdb").hexdigest()[:14]
    script = f'''\
import streamlit as st
from mn_molecular_design_studio.complex_analysis_ui import render_complex_analysis_workspace
def schema_provider(workflow):
    return {{"workflow": workflow, "schema_version": 1, "parameters": []}}
def submit(schema, workflow, parameters, project, settings, title):
    st.session_state["submitted_workflow"] = workflow
    st.session_state["submitted_parameters"] = parameters
    st.session_state["submitted_project"] = project
    st.session_state["last_submitted_job"] = {{"run_id":"analysis-ligand-1"}}
st.session_state["complex-analysis-source-kind"] = "Results"
render_complex_analysis_workspace(
    {USER!r}, workflow_schema=schema_provider,
    project_selector=lambda prefix, **kwargs: "test-project", submit=submit,
)
'''
    app = AppTest.from_string(script, default_timeout=15).run()
    assert not app.exception
    button_key = f"complex-analysis-run-{workflow_id.replace('.', '-')}-{key}"
    button = app.button(key=button_key)
    assert button is not None and not button.disabled
    button.click().run()

    assert not app.exception
    parameters = app.session_state["submitted_parameters"]
    assert app.session_state["submitted_workflow"] == workflow_id
    assert app.session_state["submitted_project"] == "test-project"
    assert parameters["candidate_ids"] == [candidate["candidate_id"]]
    assert parameters["selected_structure_path"] == str(structure_path.resolve())
    assert parameters["selected_model_id"] == "2"
    assert parameters["source_job"] == "prediction/ligand-source"
    assert parameters["ligand_selection"] == {
        "chain": "L", "residue_key": "L:1", "ccd": "LIG", "smiles": smiles,
    }
    assert "binder_chains" not in parameters and "target_chains" not in parameters


def test_contextual_matching_reuses_only_exact_accessible_evaluations(tmp_path, monkeypatch):
    from mn_molecular_design_studio import auth, complex_analysis_ui, job_catalog

    rows = []
    owners = {}
    for workflow_id, tool in (("analysis.rosetta", "rosetta_metrics"),
                              ("analysis.ipsae", "ipsae"), ("analysis.pep_patch", "pep_patch")):
        run_dir = tmp_path / tool
        (run_dir / "artifacts" / "analysis").mkdir(parents=True)
        (run_dir / "metadata.json").write_text(json.dumps({
            "upstream_task_group": "prediction", "upstream_run_id": "source", "status": "completed",
        }))
        selection = {
            "schema": "studio.complex_analysis_selection.v1", "candidate_id": "candidate-1",
            "structure_digest": "digest-1", "model_id": "model-2",
            "binder_chains": ["B"], "target_chains": ["A"],
            "parameters": {"pae_cutoff": 10.0} if workflow_id == "analysis.ipsae" else {},
        }
        if workflow_id == "analysis.pep_patch":
            manifest = {"selection": {"model_id": "model-2"}, "settings": {"ph": 7.4, "workers": 2, "protocol_revision": "abc"},
                        "items": [{"id": "item_000000", "candidate_id": "candidate-1", "structure_digest": "digest-1",
                                   "binder_chains": ["B"], "target_chains": ["A"]}]}
            (run_dir / "artifacts" / "analysis" / "batch_manifest.json").write_text(json.dumps(manifest))
            (run_dir / "artifacts" / "analysis" / "batch_status.json").write_text(json.dumps({
                "items": [{"id": "item_000000", "status": "completed"}],
            }))
        else:
            (run_dir / "artifacts" / "analysis" / "complex_analysis_selection.json").write_text(json.dumps(selection))
        owner = USER["id"] if tool != "ipsae" else OTHER_USER["id"]
        owners[tool] = owner
        rows.append({"task_group": "analysis", "run_id": tool, "run_dir": str(run_dir), "status": "completed"})

    monkeypatch.setattr(job_catalog, "search_jobs", lambda **kwargs: {"items": rows})
    monkeypatch.setattr(job_catalog, "get_job_access_record", lambda group, run: {"owner_id": owners[run]})
    monkeypatch.setattr(auth, "can_access_job", lambda user, record: record.get("owner_id") == user["id"])
    monkeypatch.setattr(auth, "accessible_project_ids", lambda user: [])

    for workflow_id, parameters in (("analysis.rosetta", {}), ("analysis.ipsae", {"pae_cutoff": 10.0}),
                                    ("analysis.pep_patch", {"ph": 7.4})):
        matches = complex_analysis_ui._matching_contextual_jobs(
            USER, workflow_id=workflow_id, source_task_group="prediction", source_run_id="source",
            candidate_id="candidate-1", structure_digest="digest-1", model_id="model-2",
            binder_chains=["B"], target_chains=[] if workflow_id == "analysis.pep_patch" else ["A"],
            parameters=parameters,
        )
        expected = {"analysis.rosetta": "rosetta_metrics", "analysis.ipsae": "ipsae", "analysis.pep_patch": "pep_patch"}[workflow_id]
        assert [item["run_id"] for item in matches] == ([] if expected == "ipsae" else [expected])

    assert complex_analysis_ui._matching_contextual_jobs(
        USER, workflow_id="analysis.pep_patch", source_task_group="prediction", source_run_id="source",
        candidate_id="candidate-1", structure_digest="digest-1", model_id="wrong-model",
        binder_chains=["B"], target_chains=[], parameters={"ph": 7.4},
    ) == []


@pytest.mark.parametrize(
    ("workflow_id", "analysis_parameters"),
    [
        ("analysis.rosetta", {}),
        ("analysis.ipsae", {"pae_cutoff": 10.0, "dist_cutoff": 10.0, "ipdae_contact_cutoff": 8.0,
                             "structure_contact_cutoff": 4.0, "hotspot_contact_cutoff": 5.0}),
        ("analysis.pep_patch", {"ph": 7.4}),
    ],
)
def test_contextual_workflows_validate_through_registered_request_schema(
    tmp_path, monkeypatch, workflow_id, analysis_parameters,
):
    from mn_molecular_design_studio.engine.workflow_engine.services import workflow_automation

    source_run = tmp_path / "prediction" / "source"
    source_run.mkdir(parents=True)
    selected_structure = source_run / "model-2.pdb"
    selected_structure.write_bytes(_complex_pdb("protein"))
    candidates = source_run / "artifacts" / "normalized_candidates" / "candidates.jsonl"
    candidates.parent.mkdir(parents=True)
    candidates.write_text(json.dumps({"candidate_id": "candidate-1", "complex_pdb": str(selected_structure)}) + "\n")
    (source_run / "metadata.json").write_text(json.dumps({"tool": "generic-predictor", "status": "completed"}))
    monkeypatch.setattr(workflow_automation, "resolve_job", lambda reference, **kwargs: {
        "run_dir": str(source_run), "run_id": "source", "task_group": "prediction", "project_id": "test-project",
    })
    parameters = {
        **analysis_parameters,
        "candidate_ids": ["candidate-1"], "selected_structure_path": str(selected_structure),
        "selected_model_id": "2", "binder_chains": ["B"],
    }
    if workflow_id != "analysis.pep_patch":
        parameters["target_chains"] = ["A"]
    request_file = tmp_path / f"{workflow_id.replace('.', '-')}-request.json"
    request_file.write_text(json.dumps({
        "schema_version": 1, "workflow": workflow_id, "source_job": "prediction/source",
        "parameters": parameters,
    }))

    result = workflow_automation.validate_workflow_request(request_file)

    assert result["valid"] is True
    assert result["workflow"] == workflow_id
    assert result["parameters"]["source_run_dir"] == str(source_run.resolve())
    assert result["parameters"]["selected_structure_path"] == str(selected_structure.resolve())
    assert result["parameters"]["candidate_ids"] == ["candidate-1"]


def test_contextual_evaluation_projection_uses_source_candidate_identity(tmp_path, monkeypatch):
    from mn_molecular_design_studio import evaluations, molecular_database_results
    from mn_molecular_design_studio.complex_analysis_ui import _analysis_evaluation_id

    monkeypatch.setattr(molecular_database_results, "list_run_molecules", lambda run_dir, user: [{
        "observations": [{"context": {"candidate": {"candidate_id": "source-candidate"}}}],
    }])
    monkeypatch.setattr(evaluations, "evaluation_descriptor", lambda observation: {
        "id": "evaluation-stable-id", "candidate_id": "source-candidate_biotite_interface",
        "evaluator": "biotite_interface",
    })

    assert _analysis_evaluation_id(tmp_path, USER, "source-candidate") == "evaluation-stable-id"


def test_database_open_preserves_stable_evaluation_and_structure_ids(monkeypatch):
    script = '''\
import streamlit as st
from mn_molecular_design_studio.complex_analysis_ui import open_database_complex
if st.button("Open selected Database structure"):
    open_database_complex("record-1", evaluation_id="evaluation-2", structure_id="model-3")
'''
    app = AppTest.from_string(script, default_timeout=10).run()
    app.button[0].click().run()

    assert not app.exception
    assert app.session_state["complex-analysis-open"] == {
        "kind": "database", "molecule_id": "record-1", "evaluation_id": "evaluation-2", "structure_id": "model-3",
    }
    assert app.session_state["page"] == "Complex Analysis"


def test_database_picker_accepts_explicit_target_provenance_for_single_component_identity():
    from mn_molecular_design_studio.complex_analysis_ui import _has_database_analysis_complex, _target_structure_notice

    binder_identity = {
        "components": [{"type": "protein", "chain": "B", "role": "binder", "sequence": "ACDEFG"}],
        "source_target_provenance": True,
        "representations": [{"id": "structure-1", "context": {"candidate_id": "candidate-1"}}],
    }
    assert _has_database_analysis_complex(binder_identity)
    assert not _has_database_analysis_complex({
        "components": binder_identity["components"],
        "representations": binder_identity["representations"],
    })
    assert "saved predicted conformation" in _target_structure_notice({"target_type": "sequence"})
    assert "stable native target fold" in _target_structure_notice({
        "validation_acceptance": {"target_type": "fasta_idr"},
    })
    assert not _target_structure_notice({"target_type": "structured"})


def test_workspace_review_form_saves_persistent_annotation(tmp_path, monkeypatch):
    database = _isolated_database(tmp_path, monkeypatch)
    imported_record = molecule_imports.parse_structure(_complex_pdb("protein"), "review-form.pdb")[0]
    imported = database.import_records(USER, [imported_record], name="Review form complex")
    record = database.list_molecules(USER, molecule_ids=imported["molecule_ids"])["items"][0]
    structure_id = record["representations"][0]["id"]
    token = hashlib.sha256("\0".join(("database", record["id"], structure_id, "")).encode()).hexdigest()[:16]
    from mn_molecular_design_studio import complex_analysis_ui, molecular_database

    monkeypatch.setattr(molecular_database, "MolecularDatabase", lambda: database)
    monkeypatch.setattr(complex_analysis_ui, "_show_structure", lambda *args, **kwargs: None)
    script = f'''\
import streamlit as st
from mn_molecular_design_studio.complex_analysis_ui import render_complex_analysis_workspace
st.session_state["complex-analysis-source-kind"] = "Database"
st.session_state["complex-analysis-open"] = {{"kind": "database", "molecule_id": {record['id']!r}, "structure_id": {structure_id!r}}}
render_complex_analysis_workspace({USER!r})
'''
    app = AppTest.from_string(script, default_timeout=10).run()
    app.selectbox(key=f"complex-analysis-review-status-{token}").set_value("Reviewed").run()
    app.text_area(key=f"complex-analysis-review-note-{token}").set_value("Human review completed.").run()
    app.button(key=f"complex-analysis-review-save-{token}").click().run()

    assert not app.exception
    saved = database.get_complex_review(
        USER, source_kind="database", subject_id=record["id"], structure_ref=structure_id,
    )
    assert saved["status"] == "Reviewed"
    assert saved["note"] == "Human review completed."


def test_database_analysis_source_is_exact_and_does_not_create_a_prediction_job(tmp_path, monkeypatch):
    database = _isolated_database(tmp_path, monkeypatch)
    imported_record = molecule_imports.parse_structure(_complex_pdb("protein"), "imported-ppi.pdb")[0]
    imported = database.import_records(USER, [imported_record], name="Database PPI")
    record = database.list_molecules(USER, molecule_ids=imported["molecule_ids"])["items"][0]
    structure = record["representations"][0]
    selection = {"molecule_id": record["id"], "structure_id": structure["id"],
                 "group_id": imported["group_id"], "model_id": structure["id"],
                 "binder_chains": ["B"], "target_chains": ["A"]}
    from mn_molecular_design_studio.complex_analysis_sources import build_database_analysis_source

    bundle = build_database_analysis_source(USER, selection, workflow_id="analysis.biotite_interface",
                                            project_id=None, database=database)
    source_dir = Path(bundle["source_run_dir"])
    candidate = json.loads(Path(bundle["candidates_jsonl"]).read_text().splitlines()[0])

    assert candidate["source_molecule_id"] == record["id"]
    assert Path(candidate["complex_pdb"]).resolve().parent == source_dir.resolve()
    assert hashlib.sha256(Path(candidate["complex_pdb"]).read_bytes()).hexdigest() == bundle["source"]["structure_digest"]
    assert not (source_dir / "result.json").exists()
    assert not (source_dir / "artifacts/normalized_candidates/candidates.jsonl").exists()

    from mn_molecular_design_studio.complex_analysis_sources import authorize_database_structure
    with pytest.raises(PermissionError):
        authorize_database_structure(OTHER_USER, selection, project_id=None, database=database,
                                     workflow_id="analysis.biotite_interface")


def test_database_ligand_adapter_keeps_exact_ligand_identity_and_drops_interface_groups(tmp_path, monkeypatch):
    database = _isolated_database(tmp_path, monkeypatch)
    from mn_molecular_design_studio.complex_analysis_sources import (
        apply_database_analysis_source, authorize_database_structure, build_database_analysis_source,
    )
    from mn_molecular_design_studio.evaluations import file_digest
    assets = Path(__file__).resolve().parents[1] / "mn_molecular_design_studio/assets/vina"
    smiles = json.loads((assets / "imatinib.json").read_text())[0]["smiles"]
    parsed = molecule_imports.parse_structure((assets / "1iep.pdb").read_bytes(), "1iep.pdb", include_ligands=True)[0]
    for component in parsed["components"]:
        if component.get("type") == "small_molecule" and component.get("ccd") == "STI":
            component["smiles"] = smiles
    imported = database.import_records(USER, [parsed], name="Saved protein-ligand reference")
    record = database.list_molecules(USER, molecule_ids=imported["molecule_ids"])["items"][0]
    structure = record["representations"][0]
    ligand_selection = {"chain": "A", "residue_key": "A:201", "ccd": "STI", "smiles": smiles}
    selector = {"molecule_id": record["id"], "structure_id": structure["id"],
                "group_id": imported["group_id"], "model_id": structure["id"],
                "binder_chains": [], "target_chains": [], "ligand_selection": ligand_selection}

    authorized = authorize_database_structure(USER, selector, project_id=None, database=database,
                                              workflow_id="analysis.plip")
    bundle = build_database_analysis_source(USER, selector, workflow_id="analysis.plip",
                                            project_id=None, database=database)
    adapted = apply_database_analysis_source({
        "schema_version": 1, "workflow": "analysis.plip", "database_source": selector,
        "parameters": {"binder_chains": [], "target_chains": [], "cpu_cores": 2},
    }, USER, None)
    candidate = json.loads(Path(bundle["candidates_jsonl"]).read_text().splitlines()[0])

    assert authorized["ligand_selection"]["chain"] == "A"
    assert authorized["ligand_selection"]["residue_key"] == "A:201"
    assert authorized["ligand_selection"]["ccd"] == "STI"
    assert authorized["ligand_selection"]["smiles"] == smiles
    assert authorized["structure_digest"] == file_digest(authorized["structure_path"])
    assert bundle["source_manifest"]["ligand_selection"] == authorized["ligand_selection"]
    assert candidate["raw_metadata"]["complex_analysis_selection"]["ligand_selection"] == authorized["ligand_selection"]
    assert adapted["workflow"] == "analysis.plip"
    assert adapted["parameters"]["ligand_selection"] == authorized["ligand_selection"]
    assert "binder_chains" not in adapted["parameters"] and "target_chains" not in adapted["parameters"]

    wrong_identity = {**selector, "ligand_selection": {**ligand_selection, "smiles": "CC"}}
    with pytest.raises(ValueError, match="does not match the saved Database component"):
        authorize_database_structure(USER, wrong_identity, project_id=None, database=database,
                                     workflow_id="analysis.plip")


def test_ligand_workflow_reuses_evaluation_with_prepared_digest_and_exact_source_digest(tmp_path, monkeypatch):
    database = _isolated_database(tmp_path, monkeypatch)
    assets = Path(__file__).resolve().parents[1] / "mn_molecular_design_studio/assets/vina"
    smiles = json.loads((assets / "imatinib.json").read_text())[0]["smiles"]
    imported_record = molecule_imports.parse_structure(
        (assets / "1iep.pdb").read_bytes(), "1iep.pdb", include_ligands=True,
    )[0]
    for component in imported_record["components"]:
        if component.get("type") == "small_molecule" and component.get("ccd") == "STI":
            component["smiles"] = smiles
    imported = database.import_records(USER, [imported_record], name="Saved protein-ligand reference")
    record = database.list_molecules(USER, molecule_ids=imported["molecule_ids"])["items"][0]
    structure = record["representations"][0]
    source_digest = hashlib.sha256(Path(structure["path"]).read_bytes()).hexdigest()
    ligand_selection = {"chain": "A", "residue_key": "A:201", "ccd": "STI", "smiles": smiles}
    observation_id = database.add_observation(
        imported["group_id"], record["id"], "plip", "plip-run",
        {"plip_success": True, "plip_interaction_count": 2},
        context={
            "representation_id": structure["id"],
            "structure_digest": "prepared-evaluator-complex-digest",
            "source_structure_digest": source_digest,
            "candidate": {"candidate_id": "selected-ligand-pose", "source_tool": "plip", "stage": "analysis",
                          "raw_metadata": {"complex_analysis_selection": {"ligand_selection": ligand_selection}}},
        },
    )
    from mn_molecular_design_studio.evaluations import evaluation_descriptor
    saved_observation = next(item for item in database.list_molecules(
        USER, molecule_ids=[record["id"]], collection="all", limit=1,
    )["items"][0]["observations"] if item["id"] == observation_id)
    evaluation_id = evaluation_descriptor(saved_observation)["id"]
    selector = {
        "molecule_id": record["id"], "structure_id": structure["id"],
        "group_id": imported["group_id"], "model_id": structure["id"],
        "evaluation_id": evaluation_id, "binder_chains": [], "target_chains": [],
        "ligand_selection": ligand_selection,
    }

    from mn_molecular_design_studio.complex_analysis_sources import authorize_database_structure
    authorized = authorize_database_structure(USER, selector, project_id=None, database=database,
                                              workflow_id="analysis.pandamap")

    assert authorized["structure_id"] == structure["id"]
    assert authorized["structure_digest"] == source_digest
    assert authorized["evaluation_id"] == evaluation_id
    assert authorized["ligand_selection"] == ligand_selection

    bad_context = {**saved_observation["context"], "source_structure_digest": "0" * 64}
    database.add_observation(
        imported["group_id"], record["id"], "plip", "plip-run-wrong-source",
        {"plip_success": True}, context={**bad_context, "candidate": saved_observation["context"]["candidate"]},
    )
    latest = database.list_molecules(USER, molecule_ids=[record["id"]], collection="all", limit=1)["items"][0]
    bad_observation = next(item for item in latest["observations"] if item["run_id"] == "plip-run-wrong-source")
    bad_selector = {**selector, "evaluation_id": evaluation_descriptor(bad_observation)["id"]}
    with pytest.raises(ValueError, match="different structure representation"):
        authorize_database_structure(USER, bad_selector, project_id=None, database=database,
                                     workflow_id="analysis.pandamap")


def test_database_workspace_submits_direct_analysis_with_stable_source_ids(tmp_path, monkeypatch):
    database = _isolated_database(tmp_path, monkeypatch)
    imported_record = molecule_imports.parse_structure(_complex_pdb("protein"), "imported-ppi.pdb")[0]
    imported = database.import_records(USER, [imported_record], name="Database PPI")
    record = database.list_molecules(USER, molecule_ids=imported["molecule_ids"])["items"][0]
    structure = record["representations"][0]
    from mn_molecular_design_studio import complex_analysis_ui, molecular_database
    monkeypatch.setattr(molecular_database, "MolecularDatabase", lambda: database)
    monkeypatch.setattr(complex_analysis_ui, "_show_structure", lambda *args, **kwargs: None)
    key = hashlib.sha256(f"database:{record['id']}:{structure['id']}".encode()).hexdigest()[:14]
    script = f'''\
import streamlit as st
from mn_molecular_design_studio.complex_analysis_ui import render_complex_analysis_workspace
def schema_provider(workflow_id):
    names = {{"analysis.biotite_interface": [{{"name":"contact_distance", "default":5.0}}],
             "analysis.rosetta": [{{"name":"nprocs", "default":4}}],
             "analysis.ipsae": [{{"name":"pae_cutoff", "default":10.0}}, {{"name":"dist_cutoff", "default":10.0}},
                 {{"name":"ipdae_contact_cutoff", "default":8.0}}, {{"name":"structure_contact_cutoff", "default":4.0}},
                 {{"name":"hotspot_contact_cutoff", "default":5.0}}],
             "analysis.pep_patch": [{{"name":"ph", "default":7.4}}]}}
    return {{"workflow":workflow_id,"schema_version":1,"parameters":names[workflow_id]}}
def submit(schema, workflow, parameters, project, settings, title):
    st.session_state["submitted_workflow"] = workflow
    st.session_state["submitted_parameters"] = parameters
    st.session_state["submitted_project"] = project
    st.session_state["last_submitted_job"] = {{"run_id":"analysis-db-1"}}
st.session_state["complex-analysis-source-kind"] = "Database"
st.session_state["complex-analysis-open"] = {{"kind":"database", "molecule_id":{record['id']!r},
    "structure_id":{structure['id']!r}}}
render_complex_analysis_workspace({USER!r}, workflow_schema=schema_provider, submit=submit)
'''
    app = AppTest.from_string(script, default_timeout=20).run()
    assert not app.exception
    app.button(key=f"complex-analysis-db-run-analysis-biotite_interface-{key}").click().run()

    assert not app.exception
    assert app.session_state["submitted_workflow"] == "analysis.biotite_interface"
    assert app.session_state["submitted_project"] == ""
    parameters = app.session_state["submitted_parameters"]
    source = parameters["database_source"]
    assert source["molecule_id"] == record["id"]
    assert source["structure_id"] == structure["id"]
    assert source["group_id"] == imported["group_id"]
    assert source["binder_chains"] == ["A"]
    assert source["target_chains"] == ["B"]


@pytest.mark.parametrize("reverse_roles", [False, True])
def test_database_analysis_ingestion_attaches_evaluation_to_exact_existing_representation(tmp_path, monkeypatch, reverse_roles):
    database = _isolated_database(tmp_path, monkeypatch)
    imported_record = molecule_imports.parse_structure(_complex_pdb("protein"), "imported-ppi.pdb")[0]
    target = database.import_records(USER, [{"name":"Target", "components":[
        component for component in imported_record["components"] if component.get("chain") == "B"
    ]}])
    target_id = target["molecule_ids"][0]
    imported_record["files"][0].update(design_target_molecule_id=target_id,
                                       binder_chains=["A"], target_chains=["B"])
    imported = database.import_records(USER, [imported_record], name="Database PPI")
    record = database.list_molecules(USER, molecule_ids=imported["molecule_ids"])["items"][0]
    structure = record["representations"][0]
    from mn_molecular_design_studio.complex_analysis_sources import build_database_analysis_source
    binders, targets = (["B"], ["A"]) if reverse_roles else (["A"], ["B"])
    adapter = build_database_analysis_source(USER, {
        "molecule_id": record["id"], "structure_id": structure["id"], "group_id": imported["group_id"],
        "model_id": structure["id"], "binder_chains": binders, "target_chains": targets,
    }, workflow_id="analysis.biotite_interface", project_id=None, database=database)
    run = tmp_path / "workspace" / "runs" / "analysis" / "database-interface-evaluation"
    candidate = {
        "candidate_id": adapter["candidate_id"] + "_biotite_interface_001",
        "source_molecule_id": record["id"], "stage": "analysis",
        "source_tool": "biotite_interface", "tool": "biotite_interface",
        "complex_pdb": adapter["selected_structure_path"], "binder_chains": binders, "target_chains": targets,
        "metrics": {"biotite_interface_residue_pair_count": 1},
        "raw_metadata": {"database_source": adapter["source_manifest"], "complex_analysis_selection": {
            "structure_digest": adapter["source"]["structure_digest"], "model_id": structure["id"],
            "binder_chains": binders, "target_chains": targets,
        }},
    }
    (run / "artifacts/normalized_candidates").mkdir(parents=True)
    (run / "artifacts/normalized_candidates/candidates.jsonl").write_text(json.dumps(candidate) + "\n")
    (run / "molecular_design_studio.json").write_text(json.dumps({
        "application": "mn-molecular-design-studio", "studio_owner_id": USER["id"], "project_id": None,
    }))
    (run / "metadata.json").write_text(json.dumps({"status": "completed", "tool": "biotite_interface",
        "database_collection_policy": "automatic"}))
    (run / "input.json").write_text(json.dumps({"params": {"contact_distance": 5.0}}))
    (run / "result.json").write_text(json.dumps({"success": True}))

    from mn_molecular_design_studio.molecular_database_ingest import ingest_run
    assert ingest_run(run, database=database, force=True) == 1

    after = database.list_molecules(USER, molecule_ids=[record["id"]])["items"][0]
    evaluations = [item for item in after["observations"] if item.get("namespace") == "biotite_interface"]
    assert len(evaluations) == 1
    assert evaluations[0]["context"]["representation_id"] == structure["id"]
    assert evaluations[0]["context"]["structure_digest"] == adapter["source"]["structure_digest"]
    assert evaluations[0]["context"]["source_structure_digest"] == adapter["source"]["structure_digest"]
    assert (evaluations[0]["context"]["target_id"] == target_id) is not reverse_roles
    assert len(after["representations"]) == 1


def test_database_ligand_analysis_ingestion_keeps_original_digest_after_evaluator_preparation(tmp_path, monkeypatch):
    database = _isolated_database(tmp_path, monkeypatch)
    from mn_molecular_design_studio.complex_analysis_sources import build_database_analysis_source
    assets = Path(__file__).resolve().parents[1] / "mn_molecular_design_studio/assets/vina"
    smiles = json.loads((assets / "imatinib.json").read_text())[0]["smiles"]
    imported_record = molecule_imports.parse_structure(
        (assets / "1iep.pdb").read_bytes(), "1iep.pdb", include_ligands=True,
    )[0]
    for component in imported_record["components"]:
        if component.get("type") == "small_molecule" and component.get("ccd") == "STI":
            component["smiles"] = smiles
    imported = database.import_records(USER, [imported_record], name="Saved protein-ligand reference")
    record = database.list_molecules(USER, molecule_ids=imported["molecule_ids"])["items"][0]
    structure = record["representations"][0]
    ligand_selection = {"chain": "A", "residue_key": "A:201", "ccd": "STI", "smiles": smiles}
    adapter = build_database_analysis_source(USER, {
        "molecule_id": record["id"], "structure_id": structure["id"], "group_id": imported["group_id"],
        "model_id": structure["id"], "binder_chains": [], "target_chains": [],
        "ligand_selection": ligand_selection,
    }, workflow_id="analysis.plip", project_id=None, database=database)
    source_digest = adapter["source"]["structure_digest"]
    prepared_output = tmp_path / "plip-prepared-complex.pdb"
    selected_structure = Path(adapter["selected_structure_path"])
    prepared_output.write_bytes(b"REMARK PLIP evaluator-prepared complex\n" + selected_structure.read_bytes())
    candidate = json.loads(Path(adapter["candidates_jsonl"]).read_text().splitlines()[0])
    candidate.update({
        "candidate_id": candidate["candidate_id"] + "_plip",
        "source_tool": "plip", "stage": "analysis", "complex_pdb": str(prepared_output),
        "metrics": {"plip_success": True, "plip_interaction_count": 8.0, "qc_status": "completed"},
        "raw_metadata": {
            **candidate["raw_metadata"],
            "analysis_subject": "ligand",
            "complex_analysis_selection": {
                **candidate["raw_metadata"]["complex_analysis_selection"],
                "structure_digest": source_digest,
                "ligand_selection": ligand_selection,
            },
            "structure_digest": hashlib.sha256(prepared_output.read_bytes()).hexdigest(),
        },
    })
    run = tmp_path / "workspace" / "runs" / "analysis" / "database-plip-evaluation"
    candidate_path = run / "artifacts/normalized_candidates/candidates.jsonl"
    candidate_path.parent.mkdir(parents=True)
    candidate_path.write_text(json.dumps(candidate) + "\n")
    (run / "molecular_design_studio.json").write_text(json.dumps({
        "application": "mn-molecular-design-studio", "studio_owner_id": USER["id"], "project_id": None,
    }))
    (run / "metadata.json").write_text(json.dumps({
        "status": "completed", "tool": "plip", "upstream_source_kind": "database",
        "database_collection_policy": "automatic",
    }))
    (run / "input.json").write_text(json.dumps({"tool": "plip", "params": {"ligand_selection": ligand_selection}}))
    (run / "result.json").write_text(json.dumps({"success": True, "execution_status": "completed"}))

    from mn_molecular_design_studio.molecular_database_ingest import ingest_run
    assert ingest_run(run, database=database, force=True) == 1

    after = database.list_molecules(USER, molecule_ids=[record["id"]])["items"][0]
    evaluations = [item for item in after["observations"] if item.get("namespace") == "plip"]
    assert len(evaluations) == 1
    assert evaluations[0]["context"]["representation_id"] == structure["id"]
    assert evaluations[0]["context"]["source_structure_digest"] == source_digest
    assert evaluations[0]["context"]["structure_digest"] != source_digest
    assert evaluations[0]["context"]["candidate"]["raw_metadata"]["database_source"]["structure_id"] == structure["id"]
    assert len(after["representations"]) == 1

    # The selected ligand still matches, but replacing every protein residue
    # with alanine makes the prepared analysis coordinates a different source
    # molecule. A valid ligand identity must not mask that receptor mismatch.
    mismatched_protein = tmp_path / "plip-mismatched-protein-complex.pdb"
    mismatched_protein.write_text("\n".join(
        (line[:17] + "ALA" + line[20:]) if line.startswith("ATOM  ") else line
        for line in prepared_output.read_text().splitlines()
    ) + "\n")
    mismatch_run = tmp_path / "workspace" / "runs" / "analysis" / "database-plip-mismatched-protein"
    mismatch_candidate = json.loads(json.dumps(candidate))
    mismatch_candidate.update(candidate_id=candidate["candidate_id"] + "_wrong-protein",
                              complex_pdb=str(mismatched_protein))
    mismatch_candidate["raw_metadata"]["structure_digest"] = hashlib.sha256(
        mismatched_protein.read_bytes()).hexdigest()
    mismatch_path = mismatch_run / "artifacts/normalized_candidates/candidates.jsonl"
    mismatch_path.parent.mkdir(parents=True)
    mismatch_path.write_text(json.dumps(mismatch_candidate) + "\n")
    for name in ("molecular_design_studio.json", "metadata.json", "input.json", "result.json"):
        (mismatch_run / name).write_text((run / name).read_text())
    assert ingest_run(mismatch_run, database=database, force=True) == 0
    after_mismatch = database.list_molecules(USER, molecule_ids=[record["id"]])["items"][0]
    assert len([item for item in after_mismatch["observations"] if item.get("namespace") == "plip"]) == 1

    bad_run = tmp_path / "workspace" / "runs" / "analysis" / "database-plip-bad-provenance"
    bad_candidate = json.loads(json.dumps(candidate))
    bad_candidate["raw_metadata"]["complex_analysis_selection"]["structure_digest"] = "0" * 64
    bad_candidate_path = bad_run / "artifacts/normalized_candidates/candidates.jsonl"
    bad_candidate_path.parent.mkdir(parents=True)
    bad_candidate_path.write_text(json.dumps(bad_candidate) + "\n")
    for name in ("molecular_design_studio.json", "metadata.json", "input.json", "result.json"):
        (bad_run / name).write_text((run / name).read_text())

    assert ingest_run(bad_run, database=database, force=True) == 0
    after_bad = database.list_molecules(USER, molecule_ids=[record["id"]])["items"][0]
    assert len([item for item in after_bad["observations"] if item.get("namespace") == "plip"]) == 1
    assert len(after_bad["representations"]) == 1


@pytest.mark.parametrize(("calculated_target", "expected_ingestions"), [("EFG", 1), ("KLM", 0)])
def test_binder_owned_complex_analysis_matches_selected_representation_composition(
        tmp_path, monkeypatch, calculated_target, expected_ingestions):
    """A binder identity may own a complex representation, whose target must also match."""
    database = _isolated_database(tmp_path, monkeypatch)
    source = tmp_path / "workspace" / "runs" / "design" / "bindcraft2-binder-complex-source"
    source_candidate = {
        "candidate_id": "bindcraft2_00001", "stage": "complex_refolding",
        "binder_sequence": "ACD", "binder_chains": ["B"], "target_chains": ["A"],
        "complex_pdb": "source-complex.pdb",
        "raw_metadata": {"validation_acceptance": {"status": "accepted"}},
        "metrics": {"native_pass_filters": True},
    }
    source_candidate_path = source / "artifacts/normalized_candidates/candidates.jsonl"
    source_candidate_path.parent.mkdir(parents=True)
    source_candidate_path.write_text(json.dumps(source_candidate) + "\n")
    (source / "source-complex.pdb").write_text(_two_protein_chain_pdb("EFG", "ACD"))
    (source / "molecular_design_studio.json").write_text(json.dumps({
        "application": "mn-molecular-design-studio", "studio_owner_id": USER["id"], "project_id": None,
    }))
    (source / "metadata.json").write_text(json.dumps({
        "status": "completed", "tool": "bindcraft2", "workflow": "design_generation",
        "database_collection_policy": "automatic",
    }))
    (source / "input.json").write_text(json.dumps({"tool": "bindcraft2", "workflow": "design_generation"}))
    (source / "result.json").write_text(json.dumps({"success": True, "execution_status": "completed"}))

    from mn_molecular_design_studio.molecular_database_ingest import ingest_run
    from mn_molecular_design_studio.complex_analysis_sources import build_database_analysis_source
    assert ingest_run(source, database=database, force=True) == 1
    molecule = next(item for item in database.list_molecules(USER, collection="all")["items"]
                    if any(component.get("type") == "protein" and component.get("sequence") == "ACD"
                           for component in item["components"]))
    assert [component.get("sequence") for component in molecule["components"]
            if component.get("type") == "protein"] == ["ACD"]
    molecules_before_analysis = database.list_molecules(USER, collection="all")["total"]
    representation = molecule["representations"][0]
    adapter = build_database_analysis_source(USER, {
        "molecule_id": molecule["id"], "structure_id": representation["id"],
        "group_id": representation["group_id"],
        "model_id": representation.get("context", {}).get("model_id") or representation["id"],
        "binder_chains": ["B"], "target_chains": ["A"],
    }, workflow_id="analysis.biotite_interface", project_id=None, database=database)

    calculated = tmp_path / "calculation-complex.pdb"
    calculated.write_text(_two_protein_chain_pdb(calculated_target, "ACD"))
    run = tmp_path / "workspace" / "runs" / "analysis" / f"target-{calculated_target.lower()}-analysis"
    candidate = {
        "candidate_id": adapter["candidate_id"] + "_biotite_interface",
        "source_molecule_id": molecule["id"], "stage": "analysis", "source_tool": "biotite_interface",
        "tool": "biotite_interface", "complex_pdb": str(calculated),
        "binder_chains": ["B"], "target_chains": ["A"],
        "metrics": {"biotite_interface_residue_pair_count": 1},
        "raw_metadata": {"database_source": adapter["source_manifest"],
                         "complex_analysis_selection": {
                             "source_structure_digest": adapter["source"]["structure_digest"],
                             "model_id": adapter["selected_model_id"],
                             "binder_chains": ["B"], "target_chains": ["A"],
                         }},
    }
    candidate_path = run / "artifacts/normalized_candidates/candidates.jsonl"
    candidate_path.parent.mkdir(parents=True)
    candidate_path.write_text(json.dumps(candidate) + "\n")
    (run / "molecular_design_studio.json").write_text(json.dumps({
        "application": "mn-molecular-design-studio", "studio_owner_id": USER["id"], "project_id": None,
    }))
    (run / "metadata.json").write_text(json.dumps({
        "status": "completed", "tool": "biotite_interface", "workflow": "analysis.biotite_interface",
        "database_collection_policy": "automatic",
    }))
    (run / "input.json").write_text(json.dumps({"tool": "biotite_interface", "workflow": "analysis.biotite_interface"}))
    (run / "result.json").write_text(json.dumps({"success": True, "execution_status": "completed"}))

    ingested = ingest_run(run, database=database, force=True)
    assert ingested == expected_ingestions
    after = database.list_molecules(USER, molecule_ids=[molecule["id"]], collection="all")["items"][0]
    evaluations = [item for item in after["observations"] if item.get("namespace") == "biotite_interface"]
    assert len(evaluations) == expected_ingestions
    assert database.list_molecules(USER, collection="all")["total"] == molecules_before_analysis
    if expected_ingestions:
        assert evaluations[0]["context"]["representation_id"] == representation["id"]


def test_interface_clustering_collection_stages_authorized_results_and_database_sources(tmp_path, monkeypatch):
    from mn_molecular_design_studio import complex_analysis_ui
    from mn_molecular_design_studio.complex_analysis_sources import build_interface_clustering_collection

    user = {**USER, "_cluster_project_ids": ["project-1"]}
    database = _isolated_database(tmp_path, monkeypatch)
    record_data = molecule_imports.parse_structure(_complex_pdb("protein"), "database-source.pdb")[0]
    imported = database.import_records(user, [record_data], name="Project PPI", project_id="project-1")
    record = database.list_molecules(user, molecule_ids=imported["molecule_ids"])['items'][0]
    representation = record["representations"][0]
    database_path = database.structure_file(user, representation["id"])
    import gzip
    results_path = tmp_path / "results-source.pdb.gz"
    results_bytes = gzip.compress(_complex_pdb("protein"))
    results_path.write_bytes(results_bytes)
    candidate = {
        "candidate_id": "result-candidate", "model_index": "model-2",
        "source_tool": "prediction-engine", "binder_chains": ["B"], "target_chains": ["A"],
        "metrics": {"iptm": 0.7}, "raw_metadata": {},
    }
    job = {"task_group": "prediction", "run_id": "result-run", "project_id": "project-1",
           "run_dir": str(tmp_path / "prediction" / "result-run"), "tool": "prediction-engine",
           "status": "completed"}
    selected = {"path": results_path, "relative_path": "models/model-2.pdb", "candidate": candidate,
                "candidate_id": candidate["candidate_id"], "label": "model-2.pdb"}
    monkeypatch.setattr(complex_analysis_ui, "_accessible_result_jobs", lambda current_user: [job])
    monkeypatch.setattr(complex_analysis_ui, "_job_structure_options", lambda current_user, current_job: [selected])
    from mn_molecular_design_studio.evaluations import file_digest
    selection = {
        "schema": "studio.interface_clustering_selection.v1",
        "sources": [
            {"kind": "results", "task_group": "prediction", "run_id": "result-run",
             "structure_ref": "models/model-2.pdb", "candidate_id": "result-candidate",
             "model_id": "model-2", "structure_digest": file_digest(results_path),
             "target_chains": ["A"], "binder_chains": ["B"], "evaluation_ids": ["result-candidate"]},
            {"kind": "database", "molecule_id": record["id"], "structure_id": representation["id"],
             "group_id": imported["group_id"], "model_id": representation["id"],
             "structure_digest": file_digest(database_path), "target_chains": ["A"],
             "binder_chains": ["B"], "evaluation_ids": []},
        ],
    }

    bundle = build_interface_clustering_collection(
        user, selection, project_id="project-1", database=database,
    )
    candidates = [json.loads(line) for line in Path(bundle["candidates_jsonl"]).read_text().splitlines()]
    provenance = bundle["collection_manifest"]["items"]

    assert len(candidates) == 2
    assert len({item["candidate_id"] for item in candidates}) == 2
    assert len({item["raw_metadata"]["target_structure_digest"] for item in candidates}) == 1
    assert {item["source_kind"] for item in provenance} == {"results", "database"}
    result_item = next(item for item in provenance if item["source_kind"] == "results")
    assert result_item["source"]["candidate_id"] == "result-candidate"
    assert result_item["evaluation_ids"] == ["result-candidate"]
    db_candidate = next(item for item in candidates if "database_source" in item["raw_metadata"])
    assert db_candidate["raw_metadata"]["database_source"]["structure_id"] == representation["id"]
    assert all(Path(item["complex_pdb"]).is_relative_to(Path(bundle["source_run_dir"])) for item in candidates)
    results_candidate = next(item for item in candidates if "database_source" not in item["raw_metadata"])
    assert results_candidate["complex_pdb"].endswith(".pdb.gz")
    assert Path(results_candidate["complex_pdb"]).read_bytes() == results_bytes

    changed_target = tmp_path / "changed-target.pdb"
    changed_lines = _complex_pdb("protein").decode().splitlines()
    changed_lines[0] = changed_lines[0][:30] + f"{1.5:8.3f}" + changed_lines[0][38:]
    changed_target.write_text("\n".join(changed_lines) + "\n")
    selected["path"] = changed_target
    selection["sources"][0]["structure_digest"] = file_digest(changed_target)
    with pytest.raises(ValueError, match="identical target structure"):
        build_interface_clustering_collection(user, selection, project_id="project-1", database=database)


def test_interface_clustering_api_request_adapter_uses_registered_input_paths(monkeypatch):
    from mn_molecular_design_studio import api, complex_analysis_sources

    monkeypatch.setattr(complex_analysis_sources, "build_interface_clustering_collection", lambda user, selection, **kwargs: {
        "collection_id": "collection-1", "source_run_dir": "/owned/inputs/collection-1",
        "candidates_jsonl": "/owned/inputs/collection-1/analysis_subjects.jsonl",
    })
    payload = {
        "schema_version": 1, "workflow": "analysis.interface_clustering",
        "interface_clustering_collection": {"schema": "studio.interface_clustering_selection.v1", "sources": [{"kind": "database"}]},
        "parameters": {"similarity_threshold": 0.6, "contact_distance": 4.5},
    }

    adapted = api._prepare_interface_clustering_collection_request(payload, USER, "project-1")

    assert "interface_clustering_collection" not in adapted
    assert adapted["workflow"] == "analysis.interface_clustering"
    assert adapted["parameters"] == {
        "similarity_threshold": 0.6, "contact_distance": 4.5,
        "source_run_dir": "/owned/inputs/collection-1",
        "candidates_jsonl": "/owned/inputs/collection-1/analysis_subjects.jsonl",
    }


def test_foldseek_collection_preserves_sources_without_interface_target_constraints(tmp_path, monkeypatch):
    from mn_molecular_design_studio.evaluations import file_digest
    from mn_molecular_design_studio.complex_analysis_sources import build_foldseek_collection

    database = _isolated_database(tmp_path, monkeypatch)
    user = {**USER, "_cluster_project_ids": []}
    imported = []
    for index in range(2):
        structure_bytes = _complex_pdb("protein", x=8.0 + index)
        if index == 1:
            lines = structure_bytes.decode().splitlines()
            # Change the selected target's sequence and author residue number.
            lines = [
                line[:17] + "SER" + line[20:22] + f"{9:4d}" + line[26:]
                if line.startswith("ATOM  ") and line[21:22] == "A" else line
                for line in lines
            ]
            structure_bytes = ("\n".join(lines) + "\n").encode()
        record_data = molecule_imports.parse_structure(structure_bytes, f"foldseek-{index}.pdb")[0]
        result = database.import_records(user, [record_data], name=f"Foldseek source {index + 1}")
        record = database.list_molecules(user, molecule_ids=result["molecule_ids"])['items'][0]
        representation = record["representations"][0]
        path = database.structure_file(user, representation["id"])
        imported.append({
            "kind": "database", "molecule_id": record["id"], "structure_id": representation["id"],
            "group_id": result["group_id"], "model_id": representation["id"],
            "structure_digest": file_digest(path), "target_chains": ["A"],
            "binder_chains": ["B"], "evaluation_ids": [],
        })

    selection = {"schema": "studio.foldseek_selection.v1", "sources": imported}
    bundle = build_foldseek_collection(user, selection, project_id=None, database=database)
    candidates = [json.loads(line) for line in Path(bundle["candidates_jsonl"]).read_text().splitlines()]
    provenance = bundle["collection_manifest"]

    assert len(candidates) == 2
    assert provenance["schema"] == "studio.foldseek_collection.v1"
    assert all(candidate["binder_chains"] == ["B"] for candidate in candidates)
    assert candidates[0]["raw_metadata"]["complex_analysis_selection"]["target_chains"] == ["A"]
    assert candidates[1]["raw_metadata"]["complex_analysis_selection"]["target_chains"] == ["A"]
    assert len({item["source_structure_digest"] for item in provenance["items"]}) == 2
    assert {item["source"]["structure_id"] for item in provenance["items"]} == {
        item["structure_id"] for item in imported
    }
    assert len({item["target_structure_digest"] for item in provenance["items"]}) == 2

    # Foldseek also accepts a selection that omits target groups entirely.
    without_targets = json.loads(json.dumps(selection))
    without_targets["sources"][1]["target_chains"] = []
    unbound_bundle = build_foldseek_collection(user, without_targets, project_id=None, database=database)
    unbound_candidates = [json.loads(line) for line in Path(unbound_bundle["candidates_jsonl"]).read_text().splitlines()]
    assert unbound_candidates[1]["target_chains"] == []

    invalid = json.loads(json.dumps(selection))
    invalid["sources"][0]["binder_chains"] = ["A", "B"]
    invalid["sources"][0]["target_chains"] = []
    with pytest.raises(ValueError, match="exactly one selected protein chain"):
        build_foldseek_collection(user, invalid, project_id=None, database=database)
