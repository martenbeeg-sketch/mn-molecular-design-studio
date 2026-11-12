from __future__ import annotations

import pytest
from streamlit.testing.v1 import AppTest

from mn_molecular_design_studio.database_ui import (
    COLLECTION_PAGE_SIZE,
    _browse_collections,
    _collection_section,
)
from mn_molecular_design_studio.molecular_database import MolecularDatabase


USER = {"id": "a" * 32, "role": "user", "_cluster_project_ids": []}


def _group(group_id, name, *, method="tool_output", kind="generated", count=1, **values):
    return {"id": group_id, "name": name, "import_method": method, "collection_kind": kind,
            "molecule_count": count, "origin": {}, "project_id": "", "updated_at": "", **values}


@pytest.mark.parametrize(("method", "kind", "expected"), [
    ("user_uploads", "curated", "Imported libraries"),
    ("user_uploads", "reference", "Imported libraries"),
    ("tool_output", "generated", "Generated results"),
    ("agent_inputs", "curated", "Saved groups"),
    ("tool_output", "reference", "Saved groups"),
    ("tool_inputs", "inputs", "Job history"),
    ("tool_output", "evaluation", "Job history"),
])
def test_collection_categories_use_provenance_not_names_or_size(method, kind, expected):
    assert _collection_section(_group("group", "Important smoke testing", method=method, kind=kind)) == expected


def test_libraries_first_is_deterministic_and_does_not_hide_singletons():
    groups = [_group(str(index), f"Accepted result {index}") for index in range(1000)]
    groups += [_group("library", "Nuclear receptor library", method="user_uploads", kind="curated", count=100),
               _group("small", "Single imported molecule", method="user_uploads", kind="curated"),
               _group("saved", "Campaign shortlist", kind="curated", count=35)]
    ordered = _browse_collections(groups)
    assert [row["id"] for row in ordered[:3]] == ["library", "small", "saved"]
    assert len(ordered) == len(groups)
    assert ordered == _browse_collections(list(reversed(groups)))


def test_search_matches_all_terms_across_project_and_source_without_merging_groups():
    groups = [_group("local", "Tau binders", project_id="project"),
              _group("peer", "Tau binders", project_id="project", origin={"_replica_node": "node", "workstation": "PC004032"})]
    projects = [{"id": "project", "name": "Student research"}]
    assert [row["id"] for row in _browse_collections(groups, query="TAU student 4032", projects=projects)] == ["peer"]
    assert len(_browse_collections(groups, query="tau", projects=projects)) == 2
    assert not _browse_collections(groups, query="tau unknown", projects=projects)


def test_category_and_sort_only_reorder_or_filter_the_supplied_authorized_metadata():
    groups = [_group("old", "Zebra", updated_at="2026-01-01", count=50),
              _group("new", "Alpha", updated_at="2026-11-12", count=2),
              _group("import", "Beta", method="user_uploads", kind="curated")]
    assert [row["id"] for row in _browse_collections(groups, sort="Recently updated")] == ["new", "old", "import"]
    assert [row["id"] for row in _browse_collections(groups, sort="Name A–Z")] == ["new", "import", "old"]
    assert [row["id"] for row in _browse_collections(groups, sort="Largest first")] == ["old", "new", "import"]
    assert [row["id"] for row in _browse_collections(groups, section="Imported libraries")] == ["import"]
    assert groups[0]["id"] == "old"


@pytest.fixture
def browser(tmp_path, monkeypatch):
    from mn_molecular_design_studio import database_ui

    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_HOME", str(tmp_path / "workspace"))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(tmp_path / "workspace"))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR", str(tmp_path / "runs"))
    database = MolecularDatabase(tmp_path / "database")
    library = database.import_records(USER, [{"name": "Ethanol", "components": [{"chain": "L", "type": "small_molecule", "smiles": "CCO"}]}], name="HY-L126 library")
    for index in range(COLLECTION_PAGE_SIZE + 3):
        group_id = database.create_group(USER, f"Tau result {index:02}", import_method="tool_output",
                                         origin={"collection_kind": "generated"})
        database.add_to_group(USER, group_id, library["molecule_ids"])
    database.import_records({**USER, "id": "b" * 32}, [{"name": "Private", "components": [{"chain": "L", "type": "small_molecule", "smiles": "CCN"}]}], name="Private foreign library")
    monkeypatch.setattr(database_ui, "MolecularDatabase", lambda: database)
    monkeypatch.setattr(database_ui, "workspace_sync_running", lambda *args, **kwargs: False)
    monkeypatch.setattr(database_ui, "sync_workspace", lambda *args, **kwargs: pytest.fail("Browsing must not start indexing"))
    monkeypatch.setattr(database_ui, "render_details", lambda *args, **kwargs: None)
    app = AppTest.from_string(f"import streamlit as st\nst.session_state.setdefault('database-browser-view', 'Collections')\nfrom mn_molecular_design_studio.database_ui import render_database\nrender_database({USER!r}, [])", default_timeout=15).run()
    assert not app.exception
    return app, database, library


def test_browser_paginates_and_keeps_selection_while_browsing(browser):
    app, database, library = browser
    assert len(app.radio(key="database-selected-group").options) == COLLECTION_PAGE_SIZE + 1
    assert "HY-L126 library" in app.radio(key="database-selected-group").options[1]
    assert not any("Private foreign" in label for label in app.radio(key="database-selected-group").options)
    assert app.selectbox(key="molecular-project-filter").options == ["All projects"]
    app.radio(key="database-selected-group").set_value(library["group_id"]).run()
    app.number_input(key="molecular-collection-page").set_value(2).run()
    assert not app.exception
    assert app.session_state["database-selected-group"] == library["group_id"]
    assert app.session_state["molecular-collection-page"] == 2
    assert database.list_molecules(USER)["total"] == 1
    assert len(database.groups(USER)) == COLLECTION_PAGE_SIZE + 4


def test_browser_search_category_and_preselection_reset_pagination(browser):
    app, database, library = browser
    result = next(row for row in database.groups(USER) if row["name"] == "Tau result 17")
    app.number_input(key="molecular-collection-page").set_value(2).run()
    app.text_input(key="molecular-group-search").set_value("HY-L126").run()
    assert not app.exception
    assert app.session_state["molecular-collection-page"] == 1
    assert len(app.radio(key="database-selected-group").options) == 2
    app.radio(key="database-selected-group").set_value(library["group_id"]).run()
    app.selectbox(key="molecular-collection-section").set_value("Generated results").run()
    assert not app.exception
    assert app.session_state["database-selected-group"] == ""
    assert any("No matching collections" in item.value for item in app.info)
    app.session_state["molecular-preselect-group"] = result["id"]
    app.session_state["molecular-preselect-collection"] = "library"
    app.run()
    assert not app.exception
    assert app.session_state["molecular-collection-section"] == "All collections"
    assert app.session_state["database-selected-group"] == result["id"]
    assert app.session_state["molecular-collection-page"] == 2
    assert any(result["name"] in item.value for item in app.markdown)


def test_browser_search_never_filters_all_molecules_silently(browser):
    app, database, library = browser
    app.text_input(key="molecular-group-search").set_value("no such collection").run()
    assert not app.exception
    assert app.radio(key="database-selected-group").options == ["All molecules"]
    assert app.session_state["database-selected-group"] == ""
    assert any("1 molecules · page 1" in item.value for item in app.caption)
