from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import subprocess

import pytest
from test_database_browser import USER, OTHER, database, _import

from mn_molecular_design_studio.binderscope import REVISION, recorded_chains, source_selector
from mn_molecular_design_studio.engine.workflow_engine.tools.binderscope import normalize_complex
from mn_molecular_design_studio.engine.workflow_engine.workflows.binderscope import validate_request


def write_complex(path, *, target_chain="X", binder_chain="Y", target_offset=10, insertion=False):
    lines, serial = [], 1
    for chain, offset, shift, residue in ((target_chain, target_offset, 0, "ALA"), (binder_chain, 101, 6, "SER")):
        for index in range(16):
            angle = index * 1.1
            center = (math.cos(angle) * 2 + shift, math.sin(angle) * 2, index * 1.5)
            for atom, dx, dy, dz, element in (("N", -.7, 0, -.5, "N"), ("CA", 0, 0, 0, "C"),
                ("C", .8, .2, .5, "C"), ("O", .9, .9, .8, "O"), ("CB", 0, 1.2, -.2, "C")):
                code = "A" if insertion and index == 5 else " "
                lines.append(f"ATOM  {serial:5d} {atom:>4s} {residue} {chain}{offset + index:4d}{code}   "
                             f"{center[0] + dx:8.3f}{center[1] + dy:8.3f}{center[2] + dz:8.3f}  1.00 80.00          {element:>2s}\n")
                serial += 1
        lines.append("TER\n")
    path.write_text("".join(lines) + "END\n")


def test_recorded_roles_are_mapped_by_unique_sequence_not_chain_order():
    stored = [{"type": "protein", "role": "target", "sequence": "AAA"},
              {"type": "protein", "role": "binder", "sequence": "SSS"}]
    parsed = [{"type": "protein", "chain": "Y", "sequence": "SSS"},
              {"type": "protein", "chain": "X", "sequence": "AAA"}]
    assert recorded_chains({}, stored, parsed) == (["Y"], ["X"])
    with pytest.raises(ValueError, match="No unambiguous"):
        recorded_chains({}, [], parsed)
    with pytest.raises(ValueError, match="uniquely"):
        recorded_chains({}, stored, parsed + [{"type": "protein", "chain": "Z", "sequence": "SSS"}])
    with pytest.raises(ValueError, match="incomplete"):
        recorded_chains({"binder_chains": ["Z"], "target_chains": ["X"]}, stored, parsed)


def test_source_selector_does_not_accept_non_design_or_missing_structure():
    assert source_selector({"outcome": {"outcome_kind": "Analysis"}})[0] is None
    assert source_selector({"outcome": {"outcome_kind": "Design"}})[0] is None


def test_saved_analysis_requires_confidence_and_topology_not_silent_assumptions(tmp_path):
    parameters = {"source_run_dir": tmp_path, "candidates_jsonl": tmp_path / "subjects.jsonl"}
    with pytest.raises(ValueError, match="pLDDT"):
        validate_request(parameters, require_files=False)
    parameters["bfactors_are_plddt"] = True
    with pytest.raises(ValueError, match="topology"):
        validate_request(parameters, require_files=False)
    parameters["component_format"] = "linear protein"
    validate_request(parameters, require_files=False)
    with pytest.raises(ValueError, match="not both"):
        validate_request({**parameters, "native_run_dirs": [tmp_path]}, require_files=False)
    with pytest.raises(ValueError, match="native campaign"):
        validate_request({**parameters, "reference_cif": tmp_path}, require_files=False)


def test_normalization_preserves_coordinates_author_numbers_and_source(tmp_path):
    import gemmi

    source, destination = tmp_path / "original.pdb", tmp_path / "canonical.pdb"
    write_complex(source, insertion=True)
    original = source.read_bytes()
    mapping, sequences = normalize_complex(source, destination, "Y", "X")
    assert sequences == {"A": "A" * 16, "B": "S" * 16}
    assert mapping["A"][5] == {"canonical_chain": "A", "canonical_number": 6, "chain": "X", "number": 15, "insertion_code": "A"}
    assert source.read_bytes() == original
    before, after = gemmi.read_structure(str(source)), gemmi.read_structure(str(destination))
    assert [(atom.pos.x, atom.pos.y, atom.pos.z) for chain in before[0] for residue in chain for atom in residue] == [
        (atom.pos.x, atom.pos.y, atom.pos.z) for chain in after[0] for residue in chain for atom in residue]
    with pytest.raises(ValueError, match="unique"):
        normalize_complex(source, tmp_path / "invalid.pdb", "A", "B")


@pytest.mark.parametrize("compressed", [False, True])
@pytest.mark.parametrize("missing_confidence", [False, True, "unknown_values"])
def test_native_cif_without_optional_columns_preserves_identity(tmp_path, compressed, missing_confidence):
    import gzip
    import gemmi

    original = tmp_path / "original.pdb"
    write_complex(original, insertion=True)
    document = gemmi.read_structure(str(original)).make_mmcif_document()
    block = document.sole_block()
    atoms = block.get_mmcif_category("_atom_site.")
    del atoms["occupancy"]
    if missing_confidence == "unknown_values":
        atoms["B_iso_or_equiv"] = [False, None] * (len(atoms["id"]) // 2)
    elif missing_confidence:
        del atoms["B_iso_or_equiv"]
    block.set_mmcif_category("_atom_site.", atoms)
    content = document.as_string().encode()
    source = tmp_path / ("native.cif.gz" if compressed else "native.cif")
    source.write_bytes(gzip.compress(content) if compressed else content)
    before = source.read_bytes()
    destination = tmp_path / "canonical.pdb"
    mapping, sequences = normalize_complex(source, destination, "Y", "X")
    assert source.read_bytes() == before
    assert sequences == {"A": "A" * 16, "B": "S" * 16}
    assert mapping["A"][5]["number"] == 15 and mapping["A"][5]["insertion_code"] == "A"
    result = gemmi.read_structure(str(destination))
    expected = gemmi.read_structure(str(original))
    assert [(atom.pos.x, atom.pos.y, atom.pos.z) for chain in result[0] for residue in chain for atom in residue] == [
        (atom.pos.x, atom.pos.y, atom.pos.z) for chain in expected[0] for residue in chain for atom in residue]
    assert all(math.isnan(atom.b_iso) if missing_confidence else atom.b_iso == 80
               for chain in result[0] for residue in chain for atom in residue)


def test_workflow_schema_registers_integrated_tile():
    from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import workflow_schema

    schema = workflow_schema("analysis.binderscope")
    names = {item["name"] for item in schema["parameters"]}
    assert {"native_run_dirs", "candidates_jsonl", "bfactors_are_plddt", "component_format"} <= names


def saved_subject(database, tmp_path):
    source = tmp_path / "design.pdb"
    write_complex(source)
    record = {"name": "Saved design", "components": [
        {"type": "protein", "chain": "X", "role": "target", "sequence": "A" * 16},
        {"type": "protein", "chain": "Y", "role": "binder", "sequence": "S" * 16}],
        "files": [{"name": "design.pdb", "content": source.read_bytes(), "model_id": "0",
                   "binder_chains": ["Y"], "target_chains": ["X"]}]}
    imported = _import(database, records=[record])
    saved = database.list_molecules(USER, molecule_ids=imported["molecule_ids"], collection="all")["items"][0]
    representation = saved["representations"][0]
    selector = {"kind": "database", "molecule_id": saved["id"], "structure_id": representation["id"],
                "group_id": imported["group_id"], "binder_chains": ["Y"], "target_chains": ["X"]}
    return imported, saved, selector


def test_collection_adapter_rechecks_access_and_checksums(database, tmp_path):
    from mn_molecular_design_studio.complex_analysis_sources import build_saved_structure_collection
    imported, saved, selector = saved_subject(database, tmp_path)
    selection = {"schema": "studio.binderscope_selection.v1", "sources": [selector]}
    adapter = build_saved_structure_collection(USER, selection, workflow_id="analysis.binderscope", project_id=None, database=database)
    manifest = json.loads((Path(adapter["source_run_dir"]) / "complex_analysis_collection.json").read_text())
    assert manifest["schema"] == "studio.binderscope_collection.v1"
    assert manifest["items"][0]["source"]["structure_id"] == selector["structure_id"]
    with pytest.raises(PermissionError):
        build_saved_structure_collection(OTHER, selection, workflow_id="analysis.binderscope", project_id=None, database=database)
    with pytest.raises(ValueError, match="changed"):
        build_saved_structure_collection(USER, {**selection, "sources": [{**selector, "structure_digest": "changed"}]},
            workflow_id="analysis.binderscope", project_id=None, database=database)


def test_browser_preparation_uses_saved_roles_and_does_not_admit_more_data(database, tmp_path):
    from mn_molecular_design_studio.binderscope import prepare_browser_selection
    imported, saved, selector = saved_subject(database, tmp_path)
    entry = {"id": "selected", "structure_id": selector["structure_id"],
             "outcome": {"name": "Saved design", "outcome_kind": "Design"},
             "result": {"molecule_id": saved["id"], "group_id": imported["group_id"],
                        "structure_ids": [selector["structure_id"]], "evaluation_id": "", "model": ""}}
    before = database.list_molecules(USER, collection="all")["total"]
    selection, excluded = prepare_browser_selection(database, USER, [entry, entry])
    assert not excluded and len(selection["sources"]) == 1
    assert selection["sources"][0]["binder_chains"] == ["Y"]
    assert selection["sources"][0]["target_chains"] == ["X"]
    assert selection["sources"][0]["structure_digest"]
    assert database.list_molecules(USER, collection="all")["total"] == before
    selection, excluded = prepare_browser_selection(database, USER, [{**entry, "structure_id": ""}])
    assert not excluded and selection["sources"][0]["structure_id"] == selector["structure_id"]
    selection, excluded = prepare_browser_selection(database, USER, [{**entry, "structure_id": "",
        "result": {**entry["result"], "structure_ids": [selector["structure_id"], "another-model"]}}])
    assert not selection["sources"] and "Multiple saved models" in excluded[0]["Reason"]


def test_vendored_source_inventory_is_self_contained_and_unchanged():
    root = Path(__file__).resolve().parents[1] / "mn_molecular_design_studio/vendor/binderscope"
    manifest = json.loads((root / "SOURCE_MANIFEST.json").read_text())
    assert manifest["upstream_revision"] == REVISION
    for relative, expected in manifest["files"].items():
        assert hashlib.sha256((root / relative).read_bytes()).hexdigest() == expected, relative
    dockerfile = (Path(__file__).resolve().parents[1] / "containers/binderscope/Dockerfile").read_text()
    assert "git clone" not in dockerfile and "COPY mn_molecular_design_studio/vendor/binderscope" in dockerfile


def test_tile_configuration_queues_only_on_explicit_action(monkeypatch):
    from streamlit.testing.v1 import AppTest
    monkeypatch.delenv("MDS_READ_ONLY_PREVIEW", raising=False)
    monkeypatch.delenv("MDS_API_BASE_URL", raising=False)
    app = AppTest.from_string('''
import streamlit as st
from mn_molecular_design_studio.binderscope_ui import render_setup
from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import workflow_schema
user = {"id": "a" * 32}
def submit(*args):
    st.session_state["submitted"] = args
render_setup(workflow_schema("analysis.binderscope"), user=user, project_selector=lambda key: "project-1", submit=submit)
''').run()
    assert not app.exception
    assert "submitted" not in app.session_state
    assert not app.radio
    assert "Advanced · original campaign inputs" in [item.label for item in app.expander]
    app.text_area(key="binderscope-native-folders").set_value("/original/campaign")
    app.text_input(key="binderscope-native-target").set_value("/original/target.pdb").run()
    app.button(key="binderscope-submit").click().run()
    assert not app.exception
    assert app.session_state["submitted"][1] == "analysis.binderscope"
    assert app.session_state["submitted"][2]["native_run_dirs"] == ["/original/campaign"]


def test_metadata_suggestions_are_per_design_and_do_not_invent_confidence():
    from mn_molecular_design_studio.binderscope import infer_assessment_settings

    structure = {"id": "structure-1", "context": {}}
    record = {"metadata": {"task": "protein_binder"}, "observations": [
        {"context": {"representation_id": "another-model", "bfactors_are_plddt": True}}]}
    settings = infer_assessment_settings(record, structure, ["B"])
    assert settings["component_format"] == "linear protein"
    assert settings["b_factor_interpretation"] == "uninterpreted"
    record["observations"].append({"context": {"representation_id": "structure-1",
        "binder_chains": ["B"],
        "settings": {"b_factor_interpretation": "plddt_0_100"}}})
    assert infer_assessment_settings(record, structure, ["B"])["b_factor_interpretation"] == "plddt_0_100"
    structure["context"] = {"component_format": "cyclic_peptide"}
    record["metadata"] = {"validation_acceptance": {"binder_format": "cyclic_peptide"}}
    assert infer_assessment_settings(record, structure, ["B"])["component_format"] == "cyclic peptide"
    structure["context"]["component_format"] = "linear_peptide"
    assert infer_assessment_settings(record, structure, ["B"])["component_format"] == "unspecified"
    structure["context"]["bfactors_are_plddt"] = False
    assert infer_assessment_settings(record, structure, ["B"])["b_factor_interpretation"] == "uninterpreted"
    assert infer_assessment_settings({"metadata": {}, "components": [{"sequence": "AAAA"}]},
        {"id": "short", "context": {}}, ["B"])["component_format"] == "unspecified"


def browser_entry(imported, saved, selector):
    return {"id": "selected", "structure_id": selector["structure_id"],
            "outcome": {"name": "Saved design", "outcome_kind": "Design"},
            "result": {"molecule_id": saved["id"], "group_id": imported["group_id"],
                       "structure_ids": [selector["structure_id"]], "evaluation_id": "", "model": ""}}


def test_settings_editor_preserves_mixed_settings_and_exact_source_identity():
    from mn_molecular_design_studio.binderscope_ui import settings_rows, edited_selection

    sources = [{"label": "First design", "molecule_id": "mol-1", "structure_id": "structure-1", "group_id": "group-1",
        "assessment": {"component_format": "linear protein", "b_factor_interpretation": "plddt_0_100"}},
        {"label": "Second design", "molecule_id": "mol-2", "structure_id": "structure-2", "group_id": "group-2",
         "assessment": {"component_format": "cyclic peptide", "b_factor_interpretation": "uninterpreted"}}]
    selection = {"schema": "studio.binderscope_selection.v1", "sources": sources}
    rows = settings_rows(sources, {})
    assert edited_selection(selection, rows) == selection
    rows[1]["Topology"] = "linear peptide"
    edited = edited_selection(selection, rows)
    assert edited["sources"][0] == sources[0]
    assert edited["sources"][1]["assessment"]["component_format"] == "linear peptide"
    assert edited["sources"][1]["structure_id"] == "structure-2"
    assert sources[1]["assessment"]["component_format"] == "cyclic peptide"
    rows[0]["Include"] = False
    selected = edited_selection(selection, rows)
    assert [source["structure_id"] for source in selected["sources"]] == ["structure-2"]
    assert selected["sources"][0]["assessment"]["component_format"] == "linear peptide"
    assert len(edited_selection(selection, rows, included_only=False)["sources"]) == 2
    rows[1]["Include"] = False
    assert not edited_selection(selection, rows)["sources"]
    assert not any(row["Include"] for row in settings_rows(sources, {}, []))
    assert [row["Include"] for row in settings_rows(sources, {}, ["1"])] == [False, True]
    completed = {("mol-1", "structure-1", "group-1"): 1}
    defaults = settings_rows(sources, completed)
    assert [row["Include"] for row in defaults] == [False, True]
    assert [row["Existing analyses"] for row in defaults] == [1, 0]
    assert [source["molecule_id"] for source in edited_selection(selection, defaults)["sources"]] == ["mol-2"]
    assert all(row["Include"] for row in settings_rows(sources, completed, ["0", "1"]))
    rows[0]["Source"] = "forged"
    with pytest.raises(ValueError, match="inventory"):
        edited_selection(selection, rows)


def test_browser_scopes_partition_jobs_without_combining_project_access():
    from mn_molecular_design_studio.binderscope import project_selections

    sources = [{"group_id": "one"}, {"group_id": "two"}, {"group_id": "personal"}]
    selection = {"schema": "studio.binderscope_selection.v1", "sources": sources}
    index = {"groups": {"one": {"project_id": "project-1"}, "two": {"project_id": "project-2"}, "personal": {}}}
    grouped = project_selections(selection, index)
    assert set(grouped) == {"project-1", "project-2", ""}
    assert grouped["project-1"]["sources"] == [sources[0]]
    with pytest.raises(PermissionError):
        project_selections(selection, {"groups": {}})


@pytest.mark.parametrize("read_only,assessed", [(False, False), (True, False), (False, True)])
def test_browser_dialog_prepares_and_queues_without_navigation(database, tmp_path, monkeypatch, read_only, assessed):
    from streamlit.testing.v1 import AppTest

    monkeypatch.delenv("MDS_READ_ONLY_PREVIEW", raising=False)
    monkeypatch.delenv("MDS_API_BASE_URL", raising=False)
    if read_only:
        monkeypatch.setenv("MDS_READ_ONLY_PREVIEW", "1")
    imported, saved, selector = saved_subject(database, tmp_path)
    entry = browser_entry(imported, saved, selector)
    if assessed:
        path = database.structure_file(USER, selector["structure_id"])
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        database.add_observation(imported["group_id"], saved["id"], "binderscope", "completed-review", {}, context={
            "representation_id": selector["structure_id"], "source_structure_digest": digest,
            "binder_chains": ["Y"], "target_chains": ["X"], "status": "completed",
            "candidate": {"raw_metadata": {"binderscope": {"source": {"input_structure_digest": digest}}}}})
    app = AppTest.from_string(f'''
import streamlit as st
from pathlib import Path
from mn_molecular_design_studio.data_browser import ReadOnlyLibrary
from mn_molecular_design_studio.binderscope_ui import _browser_dialog
from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import workflow_schema
database = ReadOnlyLibrary(Path({str(database.root)!r}))
entry = {entry!r}
index = {{"groups": {{{imported["group_id"]!r}: {{"project_id": ""}}}}}}
def submit(*args):
    st.session_state.setdefault("submitted", []).append(args)
    return {{"job": {{"run_id": "test-queued-job"}}}}
_browser_dialog(database, {USER!r}, index, [entry], "test",
                workflow_schema=workflow_schema, submit=submit)
''').run()
    assert not app.exception
    assert "submitted" not in app.session_state
    assert len(app.get("dialog")) == 1
    assert "submitted" not in app.session_state
    settings = next(frame for frame in app.dataframe if "Include" in frame.value.columns)
    column_settings = json.loads(settings.proto.columns)
    assert column_settings["Topology"]["label"] == "Binder format"
    assert column_settings["Topology"]["type_config"]["options"] == [
        {"value": "unspecified", "label": "Not recorded"}, {"value": "linear protein", "label": "Protein"},
        {"value": "linear peptide", "label": "Linear peptide"}, {"value": "cyclic peptide", "label": "Cyclic peptide (head-to-tail)"}]
    assert settings.value["Topology"].tolist() == ["unspecified"]
    submit_button = next(button for button in app.button if button.label.startswith("Queue BinderScope"))
    if assessed:
        assert submit_button.disabled and "0 designs" in submit_button.label
        settings = next(frame for frame in app.dataframe if "Include" in frame.value.columns)
        assert settings.value["Include"].tolist() == [False]
        app.button(key="binderscope-select-all-test").click().run()
        assert not app.exception
        assert not next(button for button in app.button if button.label.startswith("Queue BinderScope")).disabled
        app.button(key="binderscope-select-unassessed-test").click().run()
        assert next(button for button in app.button if button.label.startswith("Queue BinderScope")).disabled
    if read_only:
        assert submit_button.disabled
        return
    app.button(key="binderscope-clear-selection-test").click().run()
    assert not app.exception
    submit_button = next(button for button in app.button if button.label.startswith("Queue BinderScope"))
    assert submit_button.disabled and "0 designs" in submit_button.label
    assert "submitted" not in app.session_state
    app.button(key="binderscope-select-all-test").click().run()
    assert not app.exception
    submit_button = next(button for button in app.button if button.label.startswith("Queue BinderScope"))
    assert not submit_button.disabled and "1 designs" in submit_button.label
    submit_button.click().run()
    assert not app.exception
    submitted = app.session_state["submitted"]
    assert len(submitted) == 1
    assert submitted[0][1] == "analysis.binderscope"
    assert submitted[0][2]["use_recorded_settings"] is True
    source = submitted[0][2]["binderscope_collection"]["sources"][0]
    assert source["structure_id"] == selector["structure_id"]
    assert source["assessment"]["b_factor_interpretation"] == "uninterpreted"
    assert "page" not in app.session_state
    next(button for button in app.button if button.label.startswith("Queue BinderScope")).click().run()
    assert len(app.session_state["submitted"]) == 1


@pytest.mark.parametrize("count", [75, 506])
def test_browser_launcher_uses_table_rows_or_all_matching_pages(monkeypatch, count):
    from streamlit.testing.v1 import AppTest
    from mn_molecular_design_studio import binderscope_ui

    monkeypatch.setattr(binderscope_ui, "_browser_dialog", binderscope_ui._browser_dialog)

    app = AppTest.from_string('''
import streamlit as st
from mn_molecular_design_studio import binderscope_ui as ui
entries = [{"id": str(number), "outcome": {"outcome_kind": "Design"}} for number in range(COUNT)]
def dialog(database, user, index, chosen, token, **kwargs):
    st.session_state["chosen"] = [entry["id"] for entry in chosen]
ui._browser_dialog = dialog
ui.render_browser_action(None, {"id": "test-user"}, {}, entries, entries[:2], signature="scope")
'''.replace("COUNT", str(count))).run()
    assert not app.exception
    assert any(radio.key.startswith("binderscope-browser-selection-scope-")
        for popover in app.get("popover") for radio in popover.get("radio"))
    app.button(key="binderscope-browser-open-scope").click().run()
    assert app.session_state["chosen"] == ["0", "1"]
    app.radio[0].set_value("All matching designs").run()
    app.button(key="binderscope-browser-open-scope").click().run()
    assert app.session_state["chosen"] == [str(number) for number in range(count)]


def test_browser_checked_subset_changes_scope_without_marking_designs_as_assessed(monkeypatch):
    from streamlit.testing.v1 import AppTest
    from mn_molecular_design_studio import binderscope_ui

    monkeypatch.setattr(binderscope_ui, "_browser_dialog", binderscope_ui._browser_dialog)
    app = AppTest.from_string('''
import streamlit as st
from mn_molecular_design_studio import binderscope_ui as ui
entries = [{"id": str(number), "outcome": {"outcome_kind": "Design"}} for number in range(75)]
checked = st.checkbox("Test a checked row")
def dialog(database, user, index, chosen, token, **kwargs):
    st.session_state["chosen"] = [entry["id"] for entry in chosen]
ui._browser_dialog = dialog
ui.render_browser_action(None, {"id": "test-user"}, {}, entries, entries[:1] if checked else [], signature="scope")
''').run()
    app.button(key="binderscope-browser-open-scope").click().run()
    assert len(app.session_state["chosen"]) == 75
    app.checkbox[0].check().run()
    app.button(key="binderscope-browser-open-scope").click().run()
    assert app.session_state["chosen"] == ["0"]
    app.checkbox[0].uncheck().run()
    app.button(key="binderscope-browser-open-scope").click().run()
    assert len(app.session_state["chosen"]) == 75


def test_per_design_settings_survive_authorization_and_change_collection_identity(database, tmp_path):
    from mn_molecular_design_studio.complex_analysis_sources import build_saved_structure_collection
    from mn_molecular_design_studio.engine.workflow_engine.workflows.binderscope import _saved_inputs

    imported, saved, selector = saved_subject(database, tmp_path)
    settings = {"component_format": "cyclic peptide", "b_factor_interpretation": "uninterpreted"}
    selection = {"schema": "studio.binderscope_selection.v1", "sources": [{**selector, "assessment": settings}]}
    adapter = build_saved_structure_collection(USER, selection, workflow_id="analysis.binderscope", project_id=None, database=database)
    candidates = _saved_inputs(Path(adapter["source_run_dir"]), Path(adapter["candidates_jsonl"]), use_recorded_settings=True)
    assert candidates[0]["raw_metadata"]["complex_analysis_selection"]["assessment"] == settings
    changed = {**selection, "sources": [{**selector, "assessment": {**settings, "component_format": "linear peptide"}}]}
    second = build_saved_structure_collection(USER, changed, workflow_id="analysis.binderscope", project_id=None, database=database)
    assert second["source_run_dir"] != adapter["source_run_dir"]
    validate_request({**adapter, "use_recorded_settings": True}, require_files=False)
    with pytest.raises(ValueError, match="topology"):
        build_saved_structure_collection(USER, {**selection, "sources": [{**selector, "assessment": {**settings, "component_format": "guessed"}}]},
                                         workflow_id="analysis.binderscope", project_id=None, database=database)


def test_api_validates_saved_collection_without_running_it(database, tmp_path, monkeypatch):
    from mn_molecular_design_studio import api, binderscope_submission
    imported, saved, selector = saved_subject(database, tmp_path)
    monkeypatch.setattr(binderscope_submission, "ReadOnlyLibrary", lambda: database)
    payload = {"schema_version": 1, "workflow": "analysis.binderscope", "project_id": None,
               "binderscope_collection": {"schema": "studio.binderscope_selection.v1", "sources": [selector]},
               "parameters": {"bfactors_are_plddt": True, "component_format": "linear protein", "cpu_cores": 1}}
    result = api.validate_workflow(payload, base_dir=None, user=USER, credentials=None)
    assert result["valid"]
    payload["parameters"] = {"use_recorded_settings": True, "cpu_cores": 1}
    selector["assessment"] = {"component_format": "linear protein", "b_factor_interpretation": "uninterpreted"}
    assert api.validate_workflow(payload, base_dir=None, user=USER, credentials=None)["valid"]


def test_replica_preflight_accepts_exact_saved_coordinates_without_readmission(database, tmp_path):
    from mn_molecular_design_studio.binderscope import prepare_browser_selection
    imported, saved, selector = saved_subject(database, tmp_path)
    with database.connection() as connection:
        row = connection.execute("SELECT context FROM molecule_structures WHERE id=?", (selector["structure_id"],)).fetchone()
        context = {**json.loads(row["context"]), "_replica_node": "source-workstation"}
        connection.execute("UPDATE molecule_structures SET context=? WHERE id=?", (json.dumps(context), selector["structure_id"]))
    entry = {"structure_id": selector["structure_id"], "outcome": {"name": "Design", "outcome_kind": "Design"},
             "result": {"molecule_id": saved["id"], "group_id": imported["group_id"],
                        "structure_ids": [selector["structure_id"]], "evaluation_id": "", "model": ""}}
    with database.connection() as connection:
        before = [tuple(row) for row in connection.execute("SELECT * FROM molecule_members")]
    selection, excluded = prepare_browser_selection(database, USER, [entry])
    assert not excluded and selection["sources"][0]["structure_id"] == selector["structure_id"]
    with database.connection() as connection:
        assert before == [tuple(row) for row in connection.execute("SELECT * FROM molecule_members")]


def legacy_subject(database, tmp_path):
    import gemmi

    target = _import(database, name="Recorded target", run_id="target", records=[{
        "name": "Recorded target", "components": [{"type": "protein", "sequence": "A" * 16}]}])
    source = tmp_path / "native.pdb"
    write_complex(source)
    document = gemmi.read_structure(str(source)).make_mmcif_document()
    block = document.sole_block()
    atoms = block.get_mmcif_category("_atom_site.")
    del atoms["occupancy"]
    block.set_mmcif_category("_atom_site.", atoms)
    imported = _import(database, records=[{"name": "Historical design", "components": [
        {"type": "protein", "chain": "Y", "role": "binder", "sequence": "S" * 16}],
        "metadata": {"design_target_molecule_id": target["molecule_ids"][0]},
        "files": [{"name": "native.cif", "content": document.as_string().encode()}]}])
    saved = database.list_molecules(USER, molecule_ids=imported["molecule_ids"], collection="all")["items"][0]
    selector = {"kind": "database", "molecule_id": saved["id"], "group_id": imported["group_id"],
                "structure_id": saved["representations"][0]["id"]}
    return imported, saved, selector, target


def test_legacy_roles_use_accessible_sql_target_and_are_rechecked_at_submission(database, tmp_path):
    from mn_molecular_design_studio.binderscope import prepare_browser_selection
    from mn_molecular_design_studio.complex_analysis_sources import build_saved_structure_collection

    imported, saved, selector, target = legacy_subject(database, tmp_path)
    entry = browser_entry(imported, saved, selector)
    selection, excluded = prepare_browser_selection(database, USER, [entry])
    assert not excluded
    source = selection["sources"][0]
    assert (source["binder_chains"], source["target_chains"]) == (["Y"], ["X"])
    assert source["chain_role_evidence"]["target_molecule_id"] == target["molecule_ids"][0]
    assert source["assessment"]["b_factor_interpretation"] == "uninterpreted"
    adapter = build_saved_structure_collection(USER, selection, workflow_id="analysis.binderscope", project_id=None, database=database)
    manifest = json.loads((Path(adapter["source_run_dir"]) / "complex_analysis_collection.json").read_text())
    assert manifest["items"][0]["source"]["chain_role_evidence"] == source["chain_role_evidence"]
    with pytest.raises(ValueError, match="evidence changed"):
        build_saved_structure_collection(USER, {**selection, "sources": [{**source, "chain_role_evidence": {
            **source["chain_role_evidence"], "target_sequence_sha256": "forged"}}]},
            workflow_id="analysis.binderscope", project_id=None, database=database)
    with database.connection() as connection:
        connection.execute("UPDATE molecule_members SET metadata='{}' WHERE group_id=? AND molecule_id=?",
                           (selector["group_id"], selector["molecule_id"]))
    with pytest.raises(ValueError, match="No unambiguous"):
        build_saved_structure_collection(USER, selection, workflow_id="analysis.binderscope", project_id=None, database=database)


def test_legacy_preflight_is_read_only_without_hydrating_evaluation_projections(database, tmp_path, monkeypatch):
    from mn_molecular_design_studio.binderscope import prepare_browser_selection
    from mn_molecular_design_studio.data_browser import ReadOnlyLibrary

    imported, saved, selector, target = legacy_subject(database, tmp_path)
    reader = ReadOnlyLibrary(database.root)
    monkeypatch.setattr(reader, "list_molecules", lambda *args, **kwargs: pytest.fail("Preflight must not rebuild SQL projections"))
    selection, excluded = prepare_browser_selection(reader, USER, [browser_entry(imported, saved, selector)])
    assert not excluded and len(selection["sources"]) == 1


@pytest.mark.parametrize("problem", ["inaccessible", "conflicting", "mismatch", "ambiguous"])
def test_legacy_target_mapping_never_guesses_or_bypasses_access(database, tmp_path, problem):
    from mn_molecular_design_studio.binderscope import prepare_browser_selection

    imported, saved, selector, target = legacy_subject(database, tmp_path)
    with database.connection() as connection:
        if problem == "inaccessible":
            connection.execute("UPDATE molecule_groups SET owner_id=? WHERE id=?", (OTHER["id"], target["group_id"]))
        elif problem == "conflicting":
            metadata = {"design_target_molecule_id": target["molecule_ids"][0], "target_molecule_id": "other"}
            connection.execute("UPDATE molecule_members SET metadata=? WHERE group_id=?", (json.dumps(metadata), imported["group_id"]))
        else:
            components = [{"type": "protein", "sequence": "S" * 16 if problem == "ambiguous" else "G" * 16}]
            connection.execute("UPDATE molecule_members SET components=? WHERE group_id=?", (json.dumps(components), target["group_id"]))
            connection.execute("UPDATE molecules SET components=? WHERE id=?", (json.dumps(components), target["molecule_ids"][0]))
    selection, excluded = prepare_browser_selection(database, USER, [browser_entry(imported, saved, selector)])
    assert not selection["sources"] and excluded


@pytest.mark.skipif(os.environ.get("MDS_TEST_BINDERSCOPE_IMAGE") != "1", reason="Explicit local BinderScope SQL integration test")
@pytest.mark.parametrize("per_design_settings", [False, True])
def test_binderscope_result_attaches_to_original_sql_design(database, tmp_path, per_design_settings):
    from mn_molecular_design_studio.complex_analysis_sources import build_saved_structure_collection
    from mn_molecular_design_studio.engine.workflow_engine.workflows.binderscope import run_binderscope
    from mn_molecular_design_studio.molecular_database_ingest import ingest_run

    imported, saved, selector = saved_subject(database, tmp_path)
    if per_design_settings:
        selector["assessment"] = {"component_format": "linear peptide", "b_factor_interpretation": "uninterpreted"}
    original_molecules = database.list_molecules(USER, collection="all")["total"]
    adapter = build_saved_structure_collection(USER, {"schema": "studio.binderscope_selection.v1", "sources": [selector]},
        workflow_id="analysis.binderscope", project_id=None, database=database)
    settings = {"use_recorded_settings": True} if per_design_settings else {"bfactors_are_plddt": True, "component_format": "linear protein"}
    run = run_binderscope(Path(adapter["source_run_dir"]), Path(adapter["candidates_jsonl"]), cpu_cores=1, **settings)
    result = json.loads((run / "result.json").read_text())
    assert result["success"], (run / "stderr.log").read_text()
    (run / "molecular_design_studio.json").write_text(json.dumps({"application": "mn-molecular-design-studio", "studio_owner_id": USER["id"]}))
    assert ingest_run(run, database=database, force=True) == 1
    assert database.list_molecules(USER, collection="all")["total"] == original_molecules
    with database.connection() as connection:
        rows = connection.execute("SELECT id,molecule_id,group_id,context FROM molecule_observations WHERE namespace='binderscope'").fetchall()
    assert len(rows) == 1 and rows[0]["molecule_id"] == saved["id"]
    context = json.loads(rows[0]["context"])
    assert context["representation_id"] == selector["structure_id"]
    assert context["candidate"]["raw_metadata"]["binderscope"]["analysis"]["dsasa"] > 0
    if per_design_settings:
        attached = context["candidate"]["raw_metadata"]["binderscope"]
        assert attached["source"]["assessment"] == selector["assessment"]
        assert "plddt_res" not in attached["analysis"]
        assert attached["termini"]["status"] == "unavailable"
    from mn_molecular_design_studio.binderscope import saved_result
    selector = {"group_id": rows[0]["group_id"], "molecule_id": saved["id"], "observation_id": rows[0]["id"]}
    assert saved_result(database, USER, selector)["analysis"]["dsasa"] > 0
    assert saved_result(database, OTHER, selector) is None


@pytest.mark.skipif(os.environ.get("MDS_TEST_BINDERSCOPE_IMAGE") != "1", reason="Explicit local BinderScope image smoke test")
@pytest.mark.parametrize("mixed_settings", [False, True])
def test_integrated_engine_runs_offline_without_external_checkout(tmp_path, mixed_settings):
    root = tmp_path / "work"
    inputs = root / "input"
    inputs.mkdir(parents=True)
    source = inputs / "design.pdb"
    write_complex(source)
    request = {"mode": "saved", "input_root": "/work/input", "output_root": "/work/output", "cpu_cores": 1,
        "bfactors_are_plddt": True, "component_format": "linear protein", "items": [{
            "candidate_id": "saved-design-1", "label": "Saved design", "path": "/work/input/design.pdb",
            "input_structure_digest": hashlib.sha256(source.read_bytes()).hexdigest(),
            "binder_chains": ["Y"], "target_chains": ["X"], "source": {"molecule_id": "mol-original", "structure_id": "structure-original"}}]}
    if mixed_settings:
        request["use_recorded_settings"] = True
        request["items"][0]["assessment"] = {"component_format": "cyclic peptide", "b_factor_interpretation": "plddt_0_100"}
        request["items"].append({**request["items"][0], "candidate_id": "unknown-confidence-design",
            "assessment": {"component_format": "linear protein", "b_factor_interpretation": "uninterpreted"}})
    (inputs / "request.json").write_text(json.dumps(request))
    runtime = Path(__file__).resolve().parents[1] / "mn_molecular_design_studio/engine/workflow_engine/tools/binderscope.py"
    request["adapter_sha256"] = hashlib.sha256(runtime.read_bytes()).hexdigest()
    (inputs / "request.json").write_text(json.dumps(request))
    from mn_molecular_design_studio.engine.workflow_engine.workflows.quality import docker_user_args
    process = subprocess.run(["docker", "run", "--rm", "--network", "none", *docker_user_args(),
        "-e", "HOME=/tmp", "-v", f"{root}:/work", "mn-binderscope:c61ef9c",
        "python", "/opt/studio/binderscope_runtime.py", "/work/input/request.json"], capture_output=True, text=True, timeout=180)
    assert process.returncode == 0, process.stdout + process.stderr
    report = json.loads((root / "output/binderscope_report.json").read_text())
    assert report["upstream_revision"] == REVISION
    assert report["results"][0]["candidate_id"] == "saved-design-1"
    assert report["results"][0]["analysis"]["dsasa"] > 0
    assert report["results"][0]["source"]["source"]["structure_id"] == "structure-original"
    if mixed_settings:
        assert len(report["results"]) == 2
        known = next(item for item in report["results"] if item["candidate_id"] == "saved-design-1")
        unknown = next(item for item in report["results"] if item["candidate_id"] == "unknown-confidence-design")
        assert known["analysis"]["plddt_res"] == [80.] * 16
        assert "plddt_t" not in known["analysis"]
        assert known["termini"]["format"] == "cyclic peptide"
        assert known["termini"]["applies"] is False
        assert "plddt_res" not in unknown["analysis"] and "plddt" not in unknown["analysis"]["tracks"]
        assert unknown["analysis"]["tracks"]["rsa"]
        assert unknown["termini"]["status"] == "unavailable"
        assert unknown["analysis"]["dsasa"] > 0


@pytest.mark.skipif(os.environ.get("MDS_TEST_BINDERSCOPE_IMAGE") != "1", reason="Explicit offline BinderScope peer integration test")
def test_replicated_native_cif_analysis_attaches_and_replays_without_modifying_source(database, tmp_path, monkeypatch):
    from mn_molecular_design_studio.binderscope import prepare_browser_selection, saved_result
    from mn_molecular_design_studio.cluster import config, transport
    from mn_molecular_design_studio.cluster.replication import Replica
    from mn_molecular_design_studio.complex_analysis_sources import build_saved_structure_collection
    from mn_molecular_design_studio.engine.workflow_engine.workflows.binderscope import run_binderscope
    from mn_molecular_design_studio.molecular_database import MolecularDatabase, _id
    from mn_molecular_design_studio.molecular_database_ingest import ingest_run
    from mn_molecular_design_studio.molecular_database_peer import capture
    from test_molecular_database_peer import SOURCE, _sync_all, _sync_author_records

    monkeypatch.setenv("MDS_CLUSTER_HOME", str(tmp_path / "cluster"))
    imported, saved, selector, target = legacy_subject(database, tmp_path)
    original_bytes = database.structure_file(USER, selector["structure_id"]).read_bytes()
    original_digest = hashlib.sha256(original_bytes).hexdigest()
    database.make_structure_portable(USER, selector["structure_id"])
    replica = Replica()
    receiver = MolecularDatabase(tmp_path / "receiver")
    third = MolecularDatabase(tmp_path / "third")
    _sync_all(replica, database, receiver)
    _sync_author_records(replica, SOURCE, [third])
    remote_group = _id("peer-", SOURCE + ":" + imported["group_id"])
    remote_structure = _id("peer-structure-", SOURCE + ":" + selector["structure_id"])
    remote_selector = {**selector, "group_id": remote_group, "structure_id": remote_structure}
    monkeypatch.setattr(config, "load", lambda: {"mode": "disabled", "name": "evaluation-node",
                                              "peers": [{"id": SOURCE, "url": "https://owner.invalid"}]})
    monkeypatch.setattr(transport, "node_request", lambda *args, **kwargs: original_bytes)
    assert receiver.structure_file(USER, remote_structure, materialize=True).read_bytes() == original_bytes
    with receiver.connection() as connection:
        before = {table: [tuple(row) for row in connection.execute(f"SELECT * FROM {table} WHERE group_id=? ORDER BY id", (remote_group,))]
                  for table in ("molecule_structures", "molecule_observations")}
    selection, excluded = prepare_browser_selection(receiver, USER, [browser_entry(
        {**imported, "group_id": remote_group}, saved, remote_selector)])
    assert not excluded and len(selection["sources"]) == 1
    adapter = build_saved_structure_collection(USER, selection, workflow_id="analysis.binderscope", project_id=None, database=receiver)
    run = run_binderscope(Path(adapter["source_run_dir"]), Path(adapter["candidates_jsonl"]), cpu_cores=1, use_recorded_settings=True)
    result = json.loads((run / "result.json").read_text())
    assert result["success"], (run / "stderr.log").read_text()
    (run / "molecular_design_studio.json").write_text(json.dumps({"application": "mn-molecular-design-studio", "studio_owner_id": USER["id"]}))
    assert ingest_run(run, database=receiver, force=True) == 1
    assert receiver.list_molecules(USER, collection="all")["total"] == 2
    with receiver.connection() as connection:
        for table, rows in before.items():
            assert rows == [tuple(row) for row in connection.execute(f"SELECT * FROM {table} WHERE group_id=? ORDER BY id", (remote_group,))]
        observation = dict(connection.execute("SELECT * FROM molecule_observations WHERE namespace='binderscope'").fetchone())
    assert observation["molecule_id"] == saved["id"] and observation["group_id"] != remote_group
    context = json.loads(observation["context"])
    assert context["representation_id"] == remote_structure
    assert context["source_structure_digest"] == original_digest
    assert context["source_structure_owner_node"] == SOURCE
    assert context["source_structure_owner_group_id"] == imported["group_id"]
    assert context["source_structure_owner_structure_id"] == selector["structure_id"]
    assert context["target_molecule_id"] == target["molecule_ids"][0]
    assessment = saved_result(receiver, USER, {"group_id": observation["group_id"], "molecule_id": saved["id"], "observation_id": observation["id"]})
    assert assessment["analysis"]["dsasa"] > 0
    assert assessment["source"]["source"]["chain_role_evidence"] == selection["sources"][0]["chain_role_evidence"]
    assert database.structure_file(USER, selector["structure_id"]).read_bytes() == original_bytes
    capture(replica, database=receiver, node_id="evaluation-node")
    for receiver_db in (database, third):
        _sync_author_records(replica, "evaluation-node", [receiver_db])
        _sync_author_records(replica, "evaluation-node", [receiver_db])
        with receiver_db.connection() as connection:
            rows = connection.execute("SELECT molecule_id,context FROM molecule_observations WHERE namespace='binderscope'").fetchall()
        assert len(rows) == 1 and rows[0]["molecule_id"] == saved["id"]
        replayed = json.loads(rows[0]["context"])
        assert replayed["source_structure_reference_status"] == "resolved"
        assert replayed["representation_id"] == (selector["structure_id"] if receiver_db is database else remote_structure)
        assert replayed["source_structure_digest"] == original_digest
        from mn_molecular_design_studio.binderscope_presentation import result_coordinates
        with receiver_db.connection() as connection:
            attached = dict(connection.execute("SELECT id,molecule_id,group_id FROM molecule_observations WHERE namespace='binderscope'").fetchone())
        view = saved_result(receiver_db, USER, {**attached, "observation_id": attached["id"]})
        assert view["_studio_sql_binding"]["structure_id"] == replayed["representation_id"]
        assert result_coordinates(receiver_db, USER, view).encode() == original_bytes
    from mn_molecular_design_studio.molecular_database_ingest import _database_analysis_subject

    source = json.loads((Path(adapter["candidates_jsonl"])).read_text().splitlines()[0])["raw_metadata"]["database_source"]
    assert _database_analysis_subject(receiver, USER, USER["id"], "", source)
    with receiver.connection() as connection:
        connection.execute("UPDATE molecule_members SET metadata='{}' WHERE group_id=? AND molecule_id=?", (remote_group, saved["id"]))
    assert _database_analysis_subject(receiver, USER, USER["id"], "", source) is None
