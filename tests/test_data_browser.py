from __future__ import annotations

from copy import deepcopy
import sqlite3
import json

import pytest
from streamlit.testing.v1 import AppTest

from mn_molecular_design_studio.data_browser import (
    ReadOnlyLibrary, admission_measurements, branch_nodes, concise_result_label, expand_scope, folder_scopes, friendly, listed_results, load_saved_record,
    navigation_tree, producing_evaluations, result_label, saved_receptor_file, saved_results, scope_path,
)
from mn_molecular_design_studio.database_browser import build_browser_index, subject_results
from test_database_browser import OTHER, USER, _import, _record, database
from test_database_design_evaluations import _analysis, design
from test_database_outcomes import _pass_design
from test_bound_docking import bound


def inventory(database):
    index = build_browser_index(database, USER, collection="library", include_related=True)
    subject = {"id": "all", "molecule_ids": set(index["molecules"]), "group_ids": set(index["groups"])}
    raw = subject_results(database, index, subject)
    return index, raw, saved_results(index, raw)


def test_only_admitted_designs_become_design_entries_and_analyses_stay_attached(database, design):
    _analysis(database, design)
    index, raw, before = inventory(database)
    assert not [row for row in before if row["outcome_kind"] == "Design"]
    _pass_design(database, design)
    index, raw, after = inventory(database)
    designs = [row for row in after if row["outcome_kind"] == "Design"]
    assert len(designs) == 1
    assert designs[0]["activity"] == "Protein design"
    assert all(row["activity"] != "Analysis" for row in after)
    listed = listed_results(designs, structures=True)
    assert {row["structure_id"] for row in listed} == {design["structure"]}


def test_folder_levels_share_exact_result_ids_and_preserve_scope(database, design):
    _pass_design(database, design)
    index, raw, outcomes = inventory(database)
    scopes = folder_scopes(index, [])
    target = next(scope for scope in scopes.values() if scope["parent"] == "targets")
    children, selections = branch_nodes(target, outcomes,
        campaigns=[scope for scope in scopes.values() if scope["parent"] == "campaigns"])
    assert children[0]["label"] == "Protein design"
    assert children[0]["count"] == 1
    assert children[0]["children"][0]["count"] == 1
    assert children[0]["children"][0]["children"]
    for selection in selections.values():
        assert {row["id"] for row in listed_results(outcomes, selection=selection)} == {row["id"] for row in listed_results(outcomes)}
    tree, hidden = navigation_tree(scopes, active_scope=target["id"], children=children)
    assert [node["label"] for node in tree] == ["Targets", "Ligands", "Projects", "Campaigns & runs", "Collections"]
    target_node = next(node for node in tree[0]["children"] if node["id"] == target["id"])
    assert target_node["children"] == children
    assert hidden == 0


def test_private_data_never_enters_folders_or_listing(database, design):
    _pass_design(database, design)
    _import(database, user=OTHER, name="Private secret", records=[_record("KKKK", name="Secret")])
    index, raw, outcomes = inventory(database)
    assert not any("Private" in scope["label"] for scope in folder_scopes(index, []).values())
    assert not any(row["name"] == "Secret" for row in outcomes)


def test_library_reader_does_not_migrate_or_write_scientific_data(database, design):
    _pass_design(database, design)
    readonly = ReadOnlyLibrary(database.root)
    with readonly.connection() as connection:
        assert connection.execute("PRAGMA query_only").fetchone()[0] == 1
        before = connection.execute("SELECT count(*) FROM molecule_observations").fetchone()[0]
        with pytest.raises(sqlite3.OperationalError, match="readonly"):
            connection.execute("DELETE FROM molecule_observations")
    index = build_browser_index(readonly, USER, collection="library", include_related=True)
    scope = next(scope for scope in folder_scopes(index, []).values() if scope["parent"] == "projects")
    assert saved_results(index, subject_results(readonly, index, scope["subject"]))
    with database.connection() as connection:
        assert connection.execute("SELECT count(*) FROM molecule_observations").fetchone()[0] == before


def test_technical_names_are_not_presented_as_biological_names(database, design):
    _pass_design(database, design)
    index, raw, outcomes = inventory(database)
    result = next(row for row in outcomes if row["outcome_kind"] == "Design")
    for name in ["cs_" + "a" * 32, "mol-" + "a" * 32]:
        assert not friendly(name)
        assert name not in result_label({**result, "name": name})
    assert friendly("Tau R3 · 31 aa")


def test_structure_listing_never_includes_post_analysis_only_coordinates(database, design):
    _pass_design(database, design)
    index, raw, outcomes = inventory(database)
    original = next(row for row in outcomes if row["outcome_kind"] == "Design")
    second = {**deepcopy(original["primary"]), "id": "second-model", "structure_ids": ["model-2"], "model": "2"}
    original["originals"].append(second)
    original["originals"].append(deepcopy(second))
    rows = listed_results([original], structures=True)
    assert {row["structure_id"] for row in rows} == {design["structure"], "model-2"}
    assert len(rows) == 2


def test_md_results_are_not_mistaken_for_unadmitted_designs(database, design):
    index, raw, outcomes = inventory(database)
    source = next(row for row in raw if row["observation_id"] == design["observation"])
    md = {**deepcopy(source), "tool": "md.production", "kind_label": "md.production",
          "evaluation": {**source["evaluation"], "stage": "analysis", "evaluator": "md.production"}}
    imported = saved_results(index, [md])
    assert len(imported) == 1
    assert imported[0]["activity"] == "Molecular dynamics"
    assert imported[0]["eligibility"] == "Recorded"


def test_folder_search_preserves_active_path_and_limits_large_roots(database, design):
    index, raw, outcomes = inventory(database)
    scopes = folder_scopes(index, [])
    active = next(scope for scope in scopes.values() if scope["parent"] == "targets")
    for position in range(120):
        key = "target-extra-" + str(position)
        scopes[key] = {**active, "id": key, "label": f"A target {position:03}"}
    tree, hidden = navigation_tree(scopes, active_scope=active["id"], children=[], limit=20)
    assert len(tree[0]["children"]) == 20
    assert active["id"] in {row["id"] for row in tree[0]["children"]}
    assert hidden == 101
    tree, hidden = navigation_tree(scopes, active_scope=active["id"], children=[], query="A target 117")
    assert {row["label"] for row in tree[0]["children"]} == {active["label"], "A target 117"}


def test_result_search_finds_native_aliases_without_adding_analysis_rows(database, design):
    _pass_design(database, design)
    index, raw, outcomes = inventory(database)
    assert listed_results(outcomes, query="run-1")
    assert not listed_results(outcomes, query="not recorded anywhere")


def test_new_browser_has_tree_and_details_without_changing_original_browser_state(database, design, monkeypatch):
    from mn_molecular_design_studio import data_browser_ui
    _pass_design(database, design)
    _analysis(database, design, "ipsae", metrics={"ipsae_example": 0.8})
    monkeypatch.setenv("MDS_DATA_BROWSER_DATABASE_ROOT", str(database.root))
    tree_nodes = []
    monkeypatch.setattr(data_browser_ui, "render_tree", lambda nodes, **kwargs: tree_nodes.extend(nodes))
    monkeypatch.setattr(data_browser_ui, "_render_structure", lambda *args: None)
    app = AppTest.from_string(
        "import streamlit as st\nst.session_state['database-browser-view']='Collections'\n"
        "from mn_molecular_design_studio.data_browser_ui import render_data_browser\n"
        f"render_data_browser({USER!r}, [])", default_timeout=30).run()
    assert not app.exception
    assert tree_nodes
    assert app.session_state["database-browser-view"] == "Collections"
    assert not any(element.label in {"Database entries", "Originating tool"} for element in app.selectbox)
    assert any(element.label == "Start at" for element in app.selectbox)
    assert any(element.label == "Refine this level" for element in app.expander)
    assert [tab.label for tab in app.tabs] == ["Scores & analyses", "Parameters", "History"]
    assert any("ipsae_example" in str(frame.value) for frame in app.dataframe)
    assert any("confidence" in str(frame.value) for frame in app.dataframe)


def test_readonly_detail_survives_missing_projection_and_deleted_observation(database, design):
    _pass_design(database, design)
    index, raw, outcomes = inventory(database)
    selected = next(row for row in raw if row["observation_id"] == design["observation"])
    with database.connection() as connection:
        connection.execute("DELETE FROM molecule_metric_values WHERE observation_id=?", (design["observation"],))
        connection.execute("DELETE FROM molecule_evaluations WHERE observation_id=?", (design["observation"],))
    readonly = ReadOnlyLibrary(database.root)
    record = load_saved_record(readonly, USER, selected)
    assert record["observations"][0]["metrics"]["confidence"] == 0.9
    assert load_saved_record(readonly, OTHER, selected) is None
    with database.connection() as connection:
        connection.execute("DELETE FROM molecule_observations WHERE id=?", (design["observation"],))
    assert load_saved_record(readonly, USER, selected) is None


def test_project_route_preserves_project_context(database, design):
    _pass_design(database, design)
    index, raw, outcomes = inventory(database)
    scopes = folder_scopes(index, [])
    project = next(scope for scope in scopes.values() if scope.get("route") == "project")
    expand_scope(scopes, index, project["id"])
    categories = [row for row in scopes.values() if row["parent"] == project["id"]]
    assert {row["label"] for row in categories} == {"Targets & ligands", "Campaigns & runs", "Collections"}
    campaigns = next(row for row in categories if row["route"] == "project-campaigns")
    expand_scope(scopes, index, campaigns["id"])
    nested = [row for row in scopes.values() if row["parent"] == campaigns["id"]]
    assert nested
    assert all(row["subject"]["group_ids"] <= project["context_groups"] for row in nested)
    assert scope_path(scopes, nested[0]["id"])[0] == "projects"


def test_collection_route_opens_molecule_across_accessible_contexts(database):
    ligand = {"type": "small_molecule", "chain": "L", "smiles": "CCO"}
    imported = _import(database, name="Compound library", run_id="", method="user_uploads", kind="curated",
        records=[_record(name="Ethanol", components=[ligand])])
    _import(database, name="Complex", records=[_record(name="Bound ethanol", components=[ligand,
        {"type": "protein", "chain": "A", "sequence": "ACDE", "role": "target"}])])
    index, raw, outcomes = inventory(database)
    scopes = folder_scopes(index, [])
    collection = scopes["collection:" + imported["group_id"]]
    expand_scope(scopes, index, collection["id"])
    molecules = [row for row in scopes.values() if row["parent"] == collection["id"]]
    assert len(molecules) == 1 and molecules[0]["label"] == "Ethanol"
    assert len(molecules[0]["subject"]["molecule_ids"]) == 2
    tree, hidden = navigation_tree(scopes, active_scope=molecules[0]["id"], children=[], root="collections")
    assert tree[0]["label"] == "Collections"
    assert hidden == 0


@pytest.mark.parametrize(("name", "alias"), [
    ("Kaempferol", "HY-14590"),
    ("Design_pc004033_20261022T060227.043741Z", "bindcraft2_00003"),
])
def test_collection_folders_use_current_member_names_and_keep_native_aliases_searchable(database, name, alias):
    imported = _import(database, name="Curated collection", records=[_record(name=name,
        metadata={"studio_naming": {"display_name": name, "native_names": [alias]}})])
    _import(database, name="Other collection", records=[_record(name="Earlier alias")])
    index, raw, outcomes = inventory(database)
    scopes = folder_scopes(index, [])
    collection_id = "collection:" + imported["group_id"]
    expand_scope(scopes, index, collection_id)
    child = next(scope for scope in scopes.values() if scope["parent"] == collection_id)
    assert child["label"] == name
    tree, hidden = navigation_tree(scopes, active_scope=collection_id, children=[], root="collections", query=alias)
    assert tree[0]["children"][0]["children"][0]["label"] == name
    assert hidden == 0


def test_collection_folders_resolve_scoped_target_renames_instead_of_old_aliases(database, monkeypatch):
    from mn_molecular_design_studio import entity_names

    monkeypatch.setattr("mn_molecular_design_studio.auth.accessible_project_ids", lambda user: ["first", "second"])
    for project, name in [("first", "Tau R3"), ("second", "Tau repeat 3")]:
        imported = _import(database, project_id=project, method="user_uploads", kind="reference",
            records=[_record(name="Old receptor", metadata={"target_identity_type": "sequence"})])
        entity_names.rename_target(database, USER, entity_names.component_id({"type": "protein", "sequence": "ACDE"}),
            "project:" + project, name, expected_revision="")
    index, raw, outcomes = inventory(database)
    scopes = folder_scopes(index, [])
    for collection in [scope for scope in scopes.values() if scope["parent"] == "collections"]:
        expand_scope(scopes, index, collection["id"])
        child = next(scope for scope in scopes.values() if scope["parent"] == collection["id"])
        project = index["groups"][collection["entity"]["group_id"]]["project_id"]
        assert child["label"] == {"first": "Tau R3", "second": "Tau repeat 3"}[project]


def test_collection_homonyms_show_project_context_without_renaming_sql(database, monkeypatch):
    monkeypatch.setattr("mn_molecular_design_studio.auth.accessible_project_ids", lambda user: ["first", "second"])
    for project in ["first", "second"]:
        _import(database, name="Protein Library", project_id=project)
    index, raw, outcomes = inventory(database)
    projects = [{"id": "first", "name": "Tau project"}, {"id": "second", "name": "PD-L1 project"}]
    scopes = folder_scopes(index, projects)
    assert {scope["label"] for scope in scopes.values() if scope["parent"] == "collections"} == {
        "Protein Library · Tau project", "Protein Library · PD-L1 project"}
    assert {group["name"] for group in database.groups(USER, scope="all")} == {"Protein Library"}


def test_collection_search_finds_members_before_expansion_without_exposing_private_members(database):
    imported = _import(database, name="Useful library", run_id="", method="user_uploads", kind="curated",
        records=[_record("WQPR", name="Kaempferol", metadata={"studio_naming": {"native_names": ["HY-14590"]}})])
    _import(database, name="Different library", records=[_record("ACDE", name="Unrelated molecule")])
    _import(database, name="Useful library", user=OTHER, records=[_record("KLMN", name="Secret compound")])
    index, raw, outcomes = inventory(database)
    scopes = folder_scopes(index, [])
    for query in ["Kaempferol", "HY-14590", "WQPR"]:
        tree, hidden = navigation_tree(scopes, active_scope="collections", children=[], root="collections", query=query)
        assert [node["id"] for node in tree[0]["children"]] == ["collection:" + imported["group_id"]]
        assert tree[0]["children"][0]["lazy"]
        assert hidden == 0
    tree, hidden = navigation_tree(scopes, active_scope="collections", children=[], root="collections", query="Secret compound")
    assert not tree[0]["children"]


def test_reference_collections_show_saved_identities_without_asserting_a_design_pass(database):
    imported = _import(database, name="Protein Library", method="tool_inputs", kind="reference", run_id="",
        records=[_record("EFG", name="Tau R3", components=[{"type": "protein", "chain": "A", "sequence": "EFG", "role": "target"}])])
    _import(database, name="Unfiltered job inputs", method="tool_inputs", kind="inputs", records=[_record("KLMN")])
    _import(database, name="Unfiltered outputs", records=[_record("WQPR")])
    index = build_browser_index(database, USER, collection="library", include_related=True)
    scope = folder_scopes(index, [])["collection:" + imported["group_id"]]
    raw = subject_results(database, index, scope["subject"], include_inputs=True)
    assert not saved_results(index, raw)
    references = saved_results(index, raw, include_references=True)
    assert len(references) == 1
    assert references[0]["activity"] == "Reference molecules"
    assert references[0]["eligibility"] == "Reference"
    assert result_label(references[0]) == "Tau R3"
    all_scope = {"molecule_ids": set(index["molecules"]), "group_ids": set(index["groups"])}
    assert len(saved_results(index, subject_results(database, index, all_scope, include_inputs=True), include_references=True)) == 1


@pytest.mark.parametrize(("method", "kind", "include_references"), [
    ("user_uploads", "curated", False), ("tool_inputs", "reference", True),
])
def test_collection_keeps_every_saved_structure_of_each_member(database, method, kind, include_references):
    imported = _import(database, name="Compound library", method=method, kind=kind, run_id="",
        records=[_record(name="Ethanol", components=[{"type": "small_molecule", "chain": "L", "smiles": "CCO"}])])
    index, raw, outcomes = inventory(database)
    scope = folder_scopes(index, [])["collection:" + imported["group_id"]]
    raw = subject_results(database, index, scope["subject"], include_inputs=True)
    base = next(row for row in raw if row["group_id"] == imported["group_id"])
    structures = [{**base, "kind": "Structure", "id": "first", "structure_ids": ["pose-1"]},
                  {**base, "kind": "Structure", "id": "second", "structure_ids": ["pose-2"]}]
    saved = saved_results(index, structures, include_references=include_references)
    assert len(saved) == 1 and saved[0]["structure_count"] == 2
    assert {row["structure_id"] for row in listed_results(saved, structures=True)} == {"pose-1", "pose-2"}


def test_reference_only_inputs_are_not_presented_as_scientific_outcomes(database, monkeypatch):
    from mn_molecular_design_studio import data_browser_ui
    from mn_molecular_design_studio.scientific_browser_tree import scientific_scopes

    imported = _import(database, name="Protein Library", method="tool_inputs", kind="reference", run_id="",
        records=[_record("EFG", name="Tau R3", components=[{"type": "protein", "chain": "A", "sequence": "EFG", "role": "target"}])])
    monkeypatch.setenv("MDS_DATA_BROWSER_DATABASE_ROOT", str(database.root))
    monkeypatch.setattr(data_browser_ui, "_render_structure", lambda *args: None)
    monkeypatch.setattr(data_browser_ui, "render_tree", lambda *args, **kwargs: None)
    scopes = scientific_scopes(build_browser_index(database, USER, collection="library", include_related=True), [])
    index = build_browser_index(database, USER, collection="library", include_related=True)
    expand_scope(scopes, index, "project:personal")
    chosen = next(scope["id"] for scope in scopes.values() if scope["parent"] == "project:personal")
    app = AppTest.from_string("import streamlit as st\n"
        f"st.session_state['data-browser-scope'] = {chosen!r}\n"
        f"st.session_state['data-browser-scope-path'] = {['projects', 'project:personal', chosen]!r}\n"
        "from mn_molecular_design_studio.data_browser_ui import render_data_browser\n"
        f"render_data_browser({USER!r}, [])", default_timeout=30).run()
    assert not app.exception
    assert any(element.value == "All 0 matching entries · no paging" for element in app.caption)
    assert not any(element.value == "Reference molecules · Reference" for element in app.caption)
    assert app.button(key="data-browser-open-datasets")
    from mn_molecular_design_studio.datasets import catalog
    assert catalog(database, USER)[0]["id"] == imported["group_id"]


def test_collection_reference_listing_does_not_promote_legacy_docking_inputs_to_results(database):
    ligand = _record(name="Ethanol", components=[{"type": "small_molecule", "chain": "L", "smiles": "CCO"}])
    _import(database, name="Library", run_id="", method="user_uploads", kind="curated", records=[ligand])
    observations = {}
    for role in ["input", "output"]:
        run = "mnligand-" + role
        imported = _import(database, task="docking", run_id=run, kind="evaluation", records=[ligand])
        observations[role] = database.add_observation(imported["group_id"], imported["molecule_ids"][0], "vina", run, {"score": -7},
            context={"task_group": "docking", "workflow_role": "docking",
                     **({"input_role": "ligand"} if role == "input" else {})})
    index, raw, outcomes = inventory(database)
    scope = {"molecule_ids": set(index["molecules"]), "group_ids": set(index["groups"])}
    raw = subject_results(database, index, scope, include_inputs=True)
    saved = saved_results(index, raw, include_references=True)
    docking = [row for row in saved if row["activity"] == "Docking"]
    assert len(docking) == 1
    assert {original["observation_id"] for row in docking for original in row["originals"]} == {observations["output"]}


def test_nested_folder_selection_survives_sql_commits(database, design, monkeypatch):
    from mn_molecular_design_studio import data_browser_ui
    from mn_molecular_design_studio.scientific_browser_tree import scientific_scopes
    _pass_design(database, design)
    monkeypatch.setenv("MDS_DATA_BROWSER_DATABASE_ROOT", str(database.root))
    monkeypatch.setattr(data_browser_ui, "_render_structure", lambda *args: None)
    index, raw, outcomes = inventory(database)
    scopes = scientific_scopes(index, [])
    expand_scope(scopes, index, "project:personal")
    chosen = next(scope["id"] for scope in scopes.values() if scope["parent"] == "project:personal" and scope["parent"] != "collections")
    monkeypatch.setattr(data_browser_ui, "render_tree", lambda *args, **kwargs: {"selectedId": chosen, "expandedIds": [], "token": "choose-targets"})
    app = AppTest.from_string("from mn_molecular_design_studio.data_browser_ui import render_data_browser\n"
        f"render_data_browser({USER!r}, [])", default_timeout=30).run()
    assert not app.exception
    assert app.session_state["data-browser-scope"] == chosen
    database.update_metadata(USER, design["source_group_id"], design["molecule_id"], {"note": "Peer-style metadata refresh"})
    app.run()
    assert not app.exception
    assert app.session_state["data-browser-scope"] == chosen
    assert app.session_state["data-browser-scope-path"] == ["projects", "project:personal", chosen]


def test_folder_search_can_be_cleared_and_does_not_carry_to_another_folder(database, design, monkeypatch):
    from mn_molecular_design_studio import data_browser_ui

    _pass_design(database, design)
    monkeypatch.setenv("MDS_DATA_BROWSER_DATABASE_ROOT", str(database.root))
    monkeypatch.setattr(data_browser_ui, "_render_structure", lambda *args: None)
    monkeypatch.setattr(data_browser_ui, "render_tree", lambda *args, **kwargs: None)
    app = AppTest.from_string("import streamlit as st\n"
        "st.session_state.setdefault('data-browser-result-query', 'stale search before upgrade')\n"
        "from mn_molecular_design_studio.data_browser_ui import render_data_browser\n"
        f"render_data_browser({USER!r}, [])", default_timeout=30).run()
    assert not app.exception
    assert app.text_input(key="data-browser-result-query").value == ""
    assert len(app.dataframe[0].value) == 1
    app.text_input(key="data-browser-result-query").set_value("not this design").run()
    assert not app.exception
    assert any("not this design" in element.value for element in app.info)
    assert app.button(key="data-browser-reset-result-refinements")
    app.button(key="data-browser-clear-result-query").click().run()
    assert not app.exception and len(app.dataframe[0].value) == 1
    app.text_input(key="data-browser-result-query").set_value("not this design").run()
    app.selectbox(key="data-browser-entry-v2").set_value("targets").run()
    assert not app.exception
    assert app.text_input(key="data-browser-result-query").value == ""
    assert len(app.dataframe[0].value) == 1


def test_search_survives_same_folder_sql_refresh_and_empty_state_resets_refinements(database, design, monkeypatch):
    from mn_molecular_design_studio import data_browser_ui

    _pass_design(database, design)
    monkeypatch.setenv("MDS_DATA_BROWSER_DATABASE_ROOT", str(database.root))
    monkeypatch.setattr(data_browser_ui, "_render_structure", lambda *args: None)
    monkeypatch.setattr(data_browser_ui, "render_tree", lambda *args, **kwargs: None)
    app = AppTest.from_string("from mn_molecular_design_studio.data_browser_ui import render_data_browser\n"
        f"render_data_browser({USER!r}, [])", default_timeout=30).run()
    app.text_input(key="data-browser-result-query").set_value("unmatched design").run()
    database.update_metadata(USER, design["source_group_id"], design["molecule_id"], {"note": "Updated SQL metadata"})
    app.run()
    assert not app.exception
    assert app.text_input(key="data-browser-result-query").value == "unmatched design"
    app.button(key="data-browser-refresh").click().run()
    assert not app.exception
    assert app.text_input(key="data-browser-result-query").value == "unmatched design"
    app.checkbox(key="data-browser-coordinates").check().run()
    app.button(key="data-browser-reset-result-refinements").click().run()
    assert not app.exception
    assert app.text_input(key="data-browser-result-query").value == ""
    assert not app.checkbox(key="data-browser-coordinates").value
    assert len(app.dataframe[0].value) == 1


def test_folder_component_sends_the_streamlit_json_value_protocol():
    from pathlib import Path
    source = (Path(__file__).parents[1] / "mn_molecular_design_studio/components/folder_tree/index.html").read_text()
    assert "dataType:'json'" in source


def test_imported_cofolding_without_pass_flag_remains_accessible_but_is_not_passed(database):
    ligand = {"type": "small_molecule", "chain": "L", "smiles": "CCO"}
    _import(database, run_id="", kind="curated", method="user_uploads", records=[_record(name="Ethanol", components=[ligand])])
    imported = _import(database, run_id="mnligand-imported-study", task="refolding", kind="evaluation",
        records=[_record(name="HY-12345", components=[ligand])])
    database.add_observation(imported["group_id"], imported["molecule_ids"][0], "Boltz-2", "mnligand-imported-study", {"confidence": 0.8},
        context={"task_group": "refolding", "candidate": {"stage": "complex_refolding", "source_tool": "Boltz-2"}})
    index, raw, outcomes = inventory(database)
    cofolds = [row for row in outcomes if row["activity"] == "Cofolding"]
    assert len(cofolds) == 1
    assert cofolds[0]["eligibility"] == "Imported"
    assert cofolds[0]["structure_count"] == 0
    assert result_label(cofolds[0]) == "Ethanol"


def test_collection_name_search_preserves_its_molecule_children(database):
    imported = _import(database, name="Useful compound library", run_id="", method="user_uploads", kind="curated",
        records=[_record(name="Ethanol", components=[{"type": "small_molecule", "chain": "L", "smiles": "CCO"}])])
    index, raw, outcomes = inventory(database)
    scopes = folder_scopes(index, [])
    identifier = "collection:" + imported["group_id"]
    expand_scope(scopes, index, identifier)
    tree, hidden = navigation_tree(scopes, active_scope=identifier, children=[], root="collections", query="Useful compound")
    collection = next(row for row in tree[0]["children"] if row["id"] == identifier)
    assert [row["label"] for row in collection["children"]] == ["Ethanol"]


def test_concise_design_label_never_replaces_creation_time_with_run_start():
    label = concise_result_label({"display_name": "Design_pc004033_20261030T125959Z",
        "outcome_kind": "Design", "display_run_id": "20261030-123456-original", "run_id": "20261101-010203-validation"})
    assert label == "Design_pc004033_20261030T125959Z"


def test_saved_receptor_read_does_not_rebuild_evaluation_projections(bound, monkeypatch):
    database, molecule, structure_id, linked, receptor, pose = bound
    readonly = ReadOnlyLibrary(database.root)
    with readonly.connection() as connection:
        row = connection.execute("SELECT * FROM molecule_structures WHERE id=?", (structure_id,)).fetchone()
        structure = {**dict(row), "context": json.loads(row["context"])}
    monkeypatch.setattr(readonly, "list_molecules", lambda *args, **kwargs: pytest.fail("Writing reader was called"))
    admin = {**USER, "role": "admin"}
    assert saved_receptor_file(readonly, admin, structure).read_bytes() == receptor.read_bytes()
    with pytest.raises(ValueError, match="different coordinate-frame digests"):
        saved_receptor_file(readonly, admin, {**structure, "context": {**structure["context"], "target_structure_digest": "0" * 64}})
    with pytest.raises(ValueError, match="not attached to the saved target"):
        saved_receptor_file(readonly, OTHER, structure)


def test_legacy_ligand_pose_remains_visible_and_downloadable_without_guessing_receptor(bound, monkeypatch):
    from mn_molecular_design_studio.engine.workflow_engine.app.components import molstar_viewer
    database, molecule, structure_id, linked, receptor, pose = bound
    admin = {**USER, "role": "admin"}
    record = database.list_molecules(admin, molecule_ids=[molecule], collection="all", limit=1)["items"][0]
    record["structures"] = [{**record["structures"][0], "context": {"database_structure_policy": "best_docking_pose"}}]
    visualizations = []
    monkeypatch.setattr(molstar_viewer, "molstar_custom_component", lambda structures, **kwargs: visualizations.extend(structures))
    app = AppTest.from_string("from pathlib import Path\n"
        "from mn_molecular_design_studio.data_browser import ReadOnlyLibrary\n"
        "from mn_molecular_design_studio.data_browser_ui import _render_structure\n"
        f"_render_structure(ReadOnlyLibrary(Path({str(database.root)!r})), {admin!r}, {record!r}, '')", default_timeout=30).run()
    assert not app.exception
    assert len(visualizations) == 1
    assert any("Ligand coordinates only" in warning.value for warning in app.warning)
    assert app.get("download_button")


@pytest.mark.parametrize("record_type", ["small_molecule", "complex"])
def test_browser_keeps_native_cofolded_complex_without_protein_sticks_or_an_extra_receptor(database, monkeypatch, record_type):
    from test_docking_results import _cif_complex
    from mn_molecular_design_studio import data_browser_ui
    from mn_molecular_design_studio.engine.workflow_engine.app.components import molstar_viewer

    path = database.root / "complex.cif"
    text = _cif_complex()
    path.write_text(text)
    monkeypatch.setattr(ReadOnlyLibrary, "structure_file", lambda *args: path)
    monkeypatch.setattr(data_browser_ui, "saved_receptor_file", lambda *args: pytest.fail("A full cofolded complex does not need another receptor"))
    visualizations = []
    monkeypatch.setattr(molstar_viewer, "molstar_custom_component", lambda structures, **kwargs: visualizations.extend(structures))
    record = {"id": "compound", "name": "Compound", "type": record_type,
        "components": [{"type": "small_molecule"}, {"type": "protein"}],
        "structures": [{"id": "complex", "label": "Cofolded complex", "context": {"database_structure_policy": "best_prediction"}}]}
    app = AppTest.from_string("from pathlib import Path\n"
        "from mn_molecular_design_studio.data_browser import ReadOnlyLibrary\n"
        "from mn_molecular_design_studio.data_browser_ui import _render_structure\n"
        f"_render_structure(ReadOnlyLibrary(Path({str(database.root)!r})), {USER!r}, {record!r}, '')", default_timeout=30).run()
    assert not app.exception and not app.warning
    assert len(visualizations) == 1 and visualizations[0].pdb == text
    assert visualizations[0].representation_type == "cartoon" and visualizations[0].chains is None
    assert path.read_text() == text and app.get("download_button")


def test_missing_evaluation_pose_is_reported_as_a_broken_link_not_an_absent_molecule_structure():
    record = {"structures": [], "selected_evaluation": {"representation_id": "obsolete-pose"}}
    app = AppTest.from_string("from mn_molecular_design_studio.data_browser_ui import _render_structure\n"
        f"_render_structure(None, {USER!r}, {record!r}, '')").run()
    assert not app.exception
    assert any("unresolved saved-pose reference" in warning.value for warning in app.warning)
    assert not app.info


def test_folder_initially_previews_available_coordinates_without_hiding_score_only_results(database, design, monkeypatch):
    from mn_molecular_design_studio import data_browser_ui

    _pass_design(database, design)
    index, raw, outcomes = inventory(database)
    posed = next(row for row in outcomes if row["outcome_kind"] == "Design")
    assert posed["primary"]["structure_ids"]
    score_only = {**deepcopy(posed), "id": "score-only", "structure_count": 0,
                  "primary": {**deepcopy(posed["primary"]), "structure_ids": []}}
    selected = []
    monkeypatch.setenv("MDS_DATA_BROWSER_DATABASE_ROOT", str(database.root))
    monkeypatch.setattr("mn_molecular_design_studio.data_browser.saved_results", lambda *args, **kwargs: [score_only, posed])
    monkeypatch.setattr(data_browser_ui, "_render_details", lambda database, user, cache, index, entry: selected.append(entry))
    monkeypatch.setattr(data_browser_ui, "render_tree", lambda *args, **kwargs: None)
    app = AppTest.from_string("from mn_molecular_design_studio.data_browser_ui import render_data_browser\n"
        f"render_data_browser({USER!r}, [])", default_timeout=30).run()
    assert not app.exception
    assert selected[-1]["result"]["structure_ids"]
    assert any("All 2 matching entries" in caption.value for caption in app.caption)


def test_original_generation_scores_do_not_need_a_second_admission_receipt(database, design):
    _pass_design(database, design)
    index, raw, outcomes = inventory(database)
    outcome = next(row for row in outcomes if row["outcome_kind"] == "Design")
    original = outcome["originals"][0]
    generation = {**deepcopy(original), "id": "original-generation", "evaluation_id": "generation-evaluation",
        "receipt": {"role": "design_generation", "state": "Needs review"},
        "evaluation": {**original["evaluation"], "stage": "generation.backbone_sequence"}}
    analysis = {**deepcopy(generation), "id": "analysis", "receipt": {"role": "post_analysis"}}
    refold = {**deepcopy(generation), "id": "refold", "receipt": {"role": "design_validation"}}
    selected = producing_evaluations(index, outcome, [original, generation, analysis, refold])
    assert [row["id"] for row in selected] == [original["id"], "original-generation"]


def test_preserved_filter_measurements_are_candidate_checked_and_not_fabricated_evaluations():
    measurement = {"value": 0.72, "unit": "0–1", "status": "measured", "source_tool": "rfd4_proteina",
        "measurement_candidate_id": "rfd4_000001", "method": "native oracle"}
    record = {"metadata": {"source_candidate_id": "rfd4_000001", "validation_acceptance": {"measurements": {
        "common.confidence.generator.complex.iptm": measurement,
        "other_candidate": {**measurement, "measurement_candidate_id": "rfd4_000002"},
        "missing_value": {**measurement, "value": None}, "unavailable": {**measurement, "status": "unavailable"}}}}}
    before = deepcopy(record)
    rows = admission_measurements(record)
    assert len(rows) == 1 and rows[0]["Value"] == "0.72 0–1"
    assert rows[0]["Tool"] == "RFD4-Proteina"
    assert record == before


def test_peer_replayed_creation_name_collisions_are_disambiguated_without_merging_designs(database):
    from test_design_creation import evidence
    base = "Design_pc004033_20261101T123456Z"
    for number, sequence in enumerate(("ACDE", "AAAA")):
        result = _import(database, run_id="run-" + str(number), records=[_record(sequence, name="rfd4_000001",
            metadata={"design_creation": evidence(sequence=sequence, source_sha256=str(number) * 64)})])
        database.add_observation(result["group_id"], result["molecule_ids"][0], "generator", "run-" + str(number), {"score": number},
            context={"workflow_role": "design_generation", "validation_acceptance": {"status": "passed"}})
        with database.connection() as connection:
            connection.execute("UPDATE molecule_members SET name=?,metadata=json_set(metadata,'$.studio_naming.display_name',?) WHERE group_id=?",
                (base, base, result["group_id"]))
    index, raw, outcomes = inventory(database)
    designs = [row for row in outcomes if row["outcome_kind"] == "Design"]
    assert len(designs) == 2
    assert len({result_label(row) for row in designs}) == 2
    assert all(result_label(row).startswith(base + "_") for row in designs)


def test_method_folder_selects_its_own_scores_and_structures(database, design):
    _pass_design(database, design)
    index, raw, outcomes = inventory(database)
    outcome = next(row for row in outcomes if row["outcome_kind"] == "Design")
    original = outcome["primary"]
    vina = {**deepcopy(original), "id": "vina-result", "tool": "vina", "structure_ids": ["vina-pose"]}
    gnina = {**deepcopy(original), "id": "gnina-result", "tool": "gnina", "structure_ids": ["gnina-pose"]}
    mixed = {**outcome, "primary": vina, "originals": [vina, gnina], "origin_tools": ["vina", "gnina"]}
    selected = listed_results([mixed], selection={"tool": "gnina"})
    assert len(selected) == 1
    assert selected[0]["result"]["tool"] == "gnina"
    assert selected[0]["outcome"]["origin_tools"] == ["gnina"]
    structures = listed_results([mixed], selection={"tool": "gnina"}, structures=True)
    assert [row["structure_id"] for row in structures] == ["gnina-pose"]
    assert mixed["primary"] is vina and len(mixed["originals"]) == 2
    admitted = {**outcome, "primary": {**original, "tool": "validation_acceptance"},
        "originals": [{**original, "tool": "validation_acceptance"}], "origin_tools": ["rfd4_proteina"]}
    assert listed_results([admitted], selection={"tool": "rfd4_proteina"})[0]["outcome"]["origin_tools"] == ["rfd4_proteina"]
