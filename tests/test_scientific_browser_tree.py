from copy import deepcopy

import pytest

from mn_molecular_design_studio.data_browser import expand_scope, listed_results, navigation_tree, saved_results
from mn_molecular_design_studio.database_browser import build_browser_index, subject_results
from mn_molecular_design_studio.scientific_browser_tree import (
    scientific_branches, scientific_scopes, scope_counts, source_jobs, target_constructs,
)
from test_database_browser import OTHER, USER, _import, _record, database


def target(database, sequence="EFG", name="Tau R3", project=""):
    return _import(database, name="Protein Library", run_id="", method="tool_inputs", kind="reference", project_id=project,
        records=[_record(sequence, name=name, metadata={"target_identity_type": "sequence"},
                         components=[{"chain": "A", "type": "protein", "sequence": sequence, "role": "target"}])])


def batch(database, sequence, *, target_id="", run="20261101-010101-first", logical="", project="", node="local", user=USER):
    imported = _import(database, name="R3 accepted · " + run, run_id=run, project_id=project, user=user,
        method="user_uploads", kind="curated", origin={"registration_batch": "historical-accepted", "saved_selection": True,
            "studio_run_id": logical, "tool": "rfd4_proteina"},
        records=[_record(sequence, name="Design_" + sequence, components=[{"chain": "B", "type": "protein", "sequence": sequence, "role": "binder"}],
            metadata={"design_target_molecule_id": target_id, "generator_engine": "rfd4_proteina",
                      "workflow_role": "design_generation", "acceptance_status": "passed"})])
    with database.connection() as connection:
        connection.execute("UPDATE molecule_groups SET origin=json_set(origin,'$._replica_node',?) WHERE id=?", (node, imported["group_id"]))
    return imported


def inventory(database):
    index = build_browser_index(database, USER, collection="library", include_related=True)
    raw = subject_results(database, index, {"molecule_ids": set(index["molecules"]), "group_ids": set(index["groups"])})
    return index, saved_results(index, raw)


def pdb_constructs(database):
    identifiers = []
    for sequence, binder, run in [("EFGH", "ACD", "first"), ("EFGHK", "KLM", "second")]:
        reference = target(database, sequence=sequence, name="TRβ")
        target_id = reference["molecule_ids"][0]
        imported = batch(database, binder, target_id=target_id, run=run, logical=run)
        database.add_observation(imported["group_id"], imported["molecule_ids"][0], "rfd4_proteina", run,
            {"confidence": 0.9}, context={"target_molecule_id": target_id, "target_pdb_id": "3GWS",
                "target_receptor": "TRβ", "workflow_role": "design_generation", "generator_engine": "rfd4_proteina"})
        identifiers.append(target_id)
    return identifiers


def test_pdb_construct_routes_are_removed_only_from_scientific_navigation(database):
    from mn_molecular_design_studio.data_browser import folder_scopes

    identifiers = pdb_constructs(database)
    index, outcomes = inventory(database)
    original = folder_scopes(index, [])
    parent = original["targets:pdb:3GWS"]
    previous_children = [scope for scope in original.values() if scope["parent"] == parent["id"]]
    assert len(previous_children) == 2
    before = saved_results(index, subject_results(database, index, parent["subject"]))
    scopes = scientific_scopes(index, [])
    assert scopes[parent["id"]]["subject"] == parent["subject"]
    assert not any(scope["parent"] == parent["id"] for scope in scopes.values())
    after = saved_results(index, subject_results(database, index, scopes[parent["id"]]["subject"]))
    assert {outcome["id"] for outcome in before} == {outcome["id"] for outcome in after}
    assert {outcome["target_id"] for outcome in after} == set(identifiers)
    assert len([scope for scope in folder_scopes(index, []).values() if scope["parent"] == parent["id"]]) == 2
    expand_scope(scopes, index, "project:personal")
    project_target = next(scope for scope in scopes.values() if scope["parent"] == "project:personal")
    assert project_target["subject"]["target_pdb_id"] == "3GWS"
    assert not any(scope["parent"] == project_target["id"] for scope in scopes.values())


def test_construct_choices_keep_equal_length_variants_distinct_without_exposing_unlinked_targets(database):
    identifiers = ["exact-first", "exact-second"]
    index = {"molecules": {identifier: {"name": "TRβ", "target_names": set(), "import_names": set(),
        "components": [{"type": "protein", "sequence": sequence}]} for identifier, sequence in zip(identifiers, ["EFGH", "EFGK"])}}
    outcomes = [{"target_id": identifier} for identifier in [*identifiers, "unknown", ""]]
    choices = target_constructs(index, outcomes)
    assert set(choices) == set(identifiers)
    assert set(choices.values()) == {"TRβ · 4 aa · variant 1", "TRβ · 4 aa · variant 2"}
    assert choices == target_constructs(index, list(reversed(outcomes)))


def test_construct_choices_use_the_selected_receptor_context_not_a_shared_sequence_alias():
    index = {"molecules": {"target": {"name": "TRβ · 6KKB", "target_names": set(), "import_names": set(),
        "components": [{"type": "protein", "sequence": "EFGH"}]}}}
    choices = target_constructs(index, [{"target_id": "target", "target_reference_label": "TRβ · 3GWS · 4 aa"}])
    assert choices == {"target": "TRβ · 3GWS · 4 aa"}


def test_pdb_construct_refinement_filters_exact_targets_without_recreating_construct_folders(database, monkeypatch):
    from streamlit.testing.v1 import AppTest
    from mn_molecular_design_studio import data_browser_ui

    identifiers = pdb_constructs(database)
    monkeypatch.setenv("MDS_DATA_BROWSER_DATABASE_ROOT", str(database.root))
    monkeypatch.setattr(data_browser_ui, "render_tree", lambda *args, **kwargs: None)
    monkeypatch.setattr(data_browser_ui, "_render_structure", lambda *args: None)
    app = AppTest.from_string("import streamlit as st\n"
        "st.session_state.setdefault('data-browser-entry-v2', 'targets')\n"
        "st.session_state.setdefault('data-browser-scope', 'targets:pdb:3GWS')\n"
        "from mn_molecular_design_studio.data_browser_ui import render_data_browser\n"
        f"render_data_browser({USER!r}, [])", default_timeout=30).run()
    assert not app.exception
    assert len(app.dataframe[0].value) == 2
    selector = next(element for element in app.selectbox if element.label == "Target construct")
    assert selector.options == ["All constructs", "TRβ · 3GWS · 4 aa", "TRβ · 3GWS · 5 aa"]
    selector.set_value(identifiers[0]).run()
    assert not app.exception
    app.selectbox(key="data-browser-column-preset").set_value("All columns").run()
    assert not app.exception
    displayed = app.dataframe[0].value
    assert len(displayed) == 1
    assert {row["Target ID"] for row in app.session_state["data-browser-page-table"][1]} == {identifiers[0]}
    assert displayed.iloc[0]["Target construct"] == "TRβ · 3GWS · 4 aa"
    assert displayed.iloc[0]["Target sequence"] == "A: EFGH"
    next(element for element in app.selectbox if element.label == "Target construct").set_value("").run()
    assert not app.exception and len(app.dataframe[0].value) == 2


def test_cached_construct_selection_returns_to_its_pdb_target_without_losing_results(database, monkeypatch):
    from streamlit.testing.v1 import AppTest
    from mn_molecular_design_studio import data_browser_ui
    from mn_molecular_design_studio.data_browser import folder_scopes

    pdb_constructs(database)
    index, outcomes = inventory(database)
    previous = folder_scopes(index, [])
    stale = next(scope["id"] for scope in previous.values() if scope["parent"] == "targets:pdb:3GWS")
    monkeypatch.setenv("MDS_DATA_BROWSER_DATABASE_ROOT", str(database.root))
    monkeypatch.setattr(data_browser_ui, "render_tree", lambda *args, **kwargs: None)
    monkeypatch.setattr(data_browser_ui, "_render_structure", lambda *args: None)
    app = AppTest.from_string("import streamlit as st\n"
        "st.session_state.setdefault('data-browser-entry-v2', 'targets')\n"
        f"st.session_state.setdefault('data-browser-scope', {stale!r})\n"
        f"st.session_state.setdefault('data-browser-folder', {stale!r})\n"
        f"st.session_state.setdefault('data-browser-scope-path', ['targets', 'targets:pdb:3GWS', {stale!r}])\n"
        "from mn_molecular_design_studio.data_browser_ui import render_data_browser\n"
        f"render_data_browser({USER!r}, [])", default_timeout=30).run()
    assert not app.exception
    assert app.session_state["data-browser-scope"] == "targets:pdb:3GWS"
    assert len(app.dataframe[0].value) == 2


def test_construct_refinement_does_not_hide_md_when_switching_activity(database, monkeypatch):
    from streamlit.testing.v1 import AppTest
    from mn_molecular_design_studio import data_browser_ui, md_browser_ui
    from test_md_browser import _add_md

    identifiers = pdb_constructs(database)
    _add_md(database, pdb="3GWS")
    monkeypatch.setenv("MDS_DATA_BROWSER_DATABASE_ROOT", str(database.root))
    event = {}
    monkeypatch.setattr(data_browser_ui, "render_tree", lambda *args, **kwargs: event.copy())
    monkeypatch.setattr(data_browser_ui, "_render_structure", lambda *args: None)
    rendered = []
    monkeypatch.setattr(md_browser_ui, "render_md_study", lambda database, user, index, entry, selection: rendered.append(entry["outcome"]["id"]) or {})
    monkeypatch.setattr(md_browser_ui, "render_md_structure", lambda *args: None)
    app = AppTest.from_string("import streamlit as st\n"
        "st.session_state.setdefault('data-browser-entry-v2', 'targets')\n"
        "st.session_state.setdefault('data-browser-scope', 'targets:pdb:3GWS')\n"
        "from mn_molecular_design_studio.data_browser_ui import render_data_browser\n"
        f"render_data_browser({USER!r}, [])", default_timeout=30).run()
    assert not app.exception
    next(element for element in app.selectbox if element.label == "Target construct").set_value(identifiers[0]).run()
    assert not app.exception and len(app.dataframe[0].value) == 1
    scope = app.session_state["data-browser-scopes"]["targets:pdb:3GWS"]
    dynamics = next(node for node in scope["result_branches"] if node["label"] == "Molecular dynamics")
    event.update(selectedId=dynamics["id"], token="md", expandedIds=[])
    app.run()
    assert not app.exception and rendered
    assert not any(element.label == "Target construct" for element in app.selectbox)
    assert len(app.dataframe[0].value) == 1
    assert {row["Activity"] for row in app.session_state["data-browser-page-table"][1]} == {"Molecular dynamics"}


def test_project_exposes_targets_and_real_datasets_not_automatic_batches(database):
    reference = target(database)
    batch(database, "ACD", target_id=reference["molecule_ids"][0], logical="one-experiment")
    imported = _import(database, name="Compound library", method="user_uploads", kind="curated", task="import", run_id="",
        records=[_record(name="Ethanol", components=[{"chain": "L", "type": "small_molecule", "smiles": "CCO"}])])
    index, outcomes = inventory(database)
    scopes = scientific_scopes(index, [])
    assert [scope["entity"]["group_id"] for scope in scopes.values() if scope["parent"] == "collections"] == [imported["group_id"]]
    expand_scope(scopes, index, "project:personal")
    children = [scope for scope in scopes.values() if scope["parent"] == "project:personal"]
    assert {child["label"] for child in children} == {"Tau R3", "Compound library"}
    assert not scientific_branches(index, scopes["project:personal"], outcomes)[0]


def test_project_tree_children_are_precreated_without_loading_result_inventories(database, monkeypatch):
    from mn_molecular_design_studio import database_browser

    target(database)
    index, outcomes = inventory(database)
    monkeypatch.setattr(database_browser, "subject_results", lambda *args, **kwargs: pytest.fail("Tree outline must use existing metadata only"))
    scopes = scientific_scopes(index, [])
    project = scopes["project:personal"]
    assert project["loaded"] and scopes["projects"]["loaded"]
    children = [scope for scope in scopes.values() if scope["parent"] == project["id"]]
    assert children and all(not child.get("loaded") for child in children)
    tree, hidden = navigation_tree(scopes, active_scope="projects", children=[], root="projects")
    assert not tree[0]["lazy"]
    assert not tree[0]["children"][0]["lazy"] and tree[0]["children"][0]["children"]


def test_distributed_batches_share_one_recorded_experiment_and_retain_provenance(database):
    reference = target(database)
    for sequence, run, node in [("ACD", "20261101-010101-first", "worker-a"), ("KLM", "20261101-010202-second", "worker-b")]:
        batch(database, sequence, target_id=reference["molecule_ids"][0], run=run, logical="persisted-logical-run", node=node)
    index, outcomes = inventory(database)
    scope = {"id": "subject"}
    nodes, selections = scientific_branches(index, scope, outcomes)
    assert len(nodes) == 1
    assert nodes[0]["label"] == "RFD4-Proteina · 1 Nov 2026"
    assert nodes[0]["count"] == 2
    selected = selections[nodes[0]["id"]]
    assert len(listed_results(outcomes, selection=selected)) == 2
    assert len(source_jobs(index, outcomes, selected)) == 2


@pytest.mark.parametrize("shared_campaign", [True, False])
def test_legacy_pipeline_and_generator_copies_use_exact_saved_campaign_not_matching_names(database, shared_campaign):
    reference = target(database)
    for number, run in enumerate(["20261101-010101-parent", "20261101-010103-generator"]):
        imported = batch(database, "ACD", target_id=reference["molecule_ids"][0], run=run)
        with database.connection() as connection:
            connection.execute("UPDATE molecule_members SET metadata=json_set(metadata,'$.distributed_campaign.campaign_id',?) "
                "WHERE group_id=?", ("campaign-one" if shared_campaign or number == 0 else "campaign-two", imported["group_id"]))
    index, outcomes = inventory(database)
    assert len(outcomes) == 1
    nodes, selections = scientific_branches(index, {"id": "subject"}, outcomes)
    assert len(nodes) == (1 if shared_campaign else 2)
    assert all(node["count"] == len(listed_results(outcomes, selection=selections[node["id"]])) == 1 for node in nodes)
    if shared_campaign:
        selection = selections[nodes[0]["id"]]
        assert len(selection["original_ids"]) == 2
        assert len(source_jobs(index, outcomes, selection)) == 2


def test_legacy_campaign_identifiers_do_not_cross_project_boundaries(database, monkeypatch):
    monkeypatch.setattr("mn_molecular_design_studio.auth.accessible_project_ids", lambda user: ["first", "second"])
    for project, sequence in [("first", "ACD"), ("second", "KLM")]:
        imported = batch(database, sequence, project=project)
        with database.connection() as connection:
            connection.execute("UPDATE molecule_members SET metadata=json_set(metadata,'$.distributed_campaign.campaign_id',?) "
                "WHERE group_id=?", ("same-legacy-identifier", imported["group_id"]))
    index, outcomes = inventory(database)
    nodes, selections = scientific_branches(index, {"id": "subject"}, outcomes)
    assert len(nodes) == 2
    assert all(len(listed_results(outcomes, selection=selection)) == 1 for selection in selections.values())


@pytest.mark.parametrize("boundary", ["logical", "node", "project"])
def test_matching_names_and_dates_do_not_merge_unrelated_experiments(database, monkeypatch, boundary):
    monkeypatch.setattr("mn_molecular_design_studio.auth.accessible_project_ids", lambda user: ["first", "second"])
    for number, sequence in enumerate(["ACD", "KLM"]):
        batch(database, sequence, logical=("experiment-" + str(number) if boundary == "logical" else ""),
            project=("first" if not number else "second") if boundary == "project" else "",
            node="worker-" + str(number) if boundary == "node" else "local")
    index, outcomes = inventory(database)
    nodes, selections = scientific_branches(index, {"id": "subject"}, outcomes)
    assert len(nodes) == 2
    assert len({node["label"] for node in nodes}) == 2
    assert {len(listed_results(outcomes, selection=selection)) for selection in selections.values()} == {1}


def test_private_names_do_not_enter_scientific_tree_or_search(database):
    target(database)
    batch(database, "ACD")
    batch(database, "KLM", logical="secret-experiment", user=OTHER)
    index, outcomes = inventory(database)
    scopes = scientific_scopes(index, [])
    expand_scope(scopes, index, "project:personal")
    assert all("KLM" not in scope.get("search", "") for scope in scopes.values())
    assert len(outcomes) == 1


def test_project_counts_match_exact_target_queries_and_unassigned_results_are_accessible(database):
    first = target(database)
    second = target(database, sequence="WQPR", name="Other construct")
    batch(database, "ACD", target_id=first["molecule_ids"][0], logical="first")
    batch(database, "KLM", target_id=second["molecule_ids"][0], logical="second")
    batch(database, "STV", logical="unknown-subject")
    index, outcomes = inventory(database)
    scopes = scientific_scopes(index, [])
    project = scopes["project:personal"]
    expand_scope(scopes, index, project["id"])
    scope_counts(scopes, index, project, outcomes)
    children = [scope for scope in scopes.values() if scope["parent"] == project["id"]]
    assert len(children) == 2
    assert all(child["result_count"] == len(saved_results(index, subject_results(database, index, child["subject"]))) == 1 for child in children)
    nodes, selections = scientific_branches(index, project, outcomes)
    assert nodes[0]["label"] == "Unassigned results" and nodes[0]["count"] == 1
    assert len(listed_results(outcomes, selection=selections[nodes[0]["id"]])) == 1
    tree, hidden = navigation_tree(scopes, active_scope=project["id"], children=nodes, root="projects")
    assert tree[0]["children"][0]["count"] == 3


def test_nested_library_compound_keeps_project_context(database, monkeypatch):
    monkeypatch.setattr("mn_molecular_design_studio.auth.accessible_project_ids", lambda user: ["first", "second"])
    ligand = _record(name="Ethanol", components=[{"chain": "L", "type": "small_molecule", "smiles": "CCO"}])
    imported = _import(database, project_id="first", name="Library", method="user_uploads", kind="curated", task="import", run_id="", records=[ligand])
    _import(database, project_id="second", name="Other library", method="user_uploads", kind="curated", task="import", run_id="", records=[ligand])
    index, outcomes = inventory(database)
    scopes = scientific_scopes(index, [{"id": "first", "name": "First"}, {"id": "second", "name": "Second"}])
    expand_scope(scopes, index, "project:first")
    library = next(scope for scope in scopes.values() if scope["parent"] == "project:first")
    expand_scope(scopes, index, library["id"])
    compound = next(scope for scope in scopes.values() if scope["parent"] == library["id"])
    assert compound["subject"]["group_ids"] == {imported["group_id"]}
    assert compound["label"] == "Ethanol"


def test_experiment_selection_does_not_leak_other_evaluations_or_structures(database):
    batch(database, "ACD")
    index, outcomes = inventory(database)
    original = outcomes[0]
    first = {**original["primary"], "id": "first-evaluation", "structure_ids": ["first-pose"]}
    second = {**deepcopy(first), "id": "second-evaluation", "structure_ids": ["second-pose"]}
    combined = {**original, "primary": second, "originals": [first, second]}
    selected = listed_results([combined], selection={"original_ids": {first["id"]}}, structures=True)
    assert len(selected) == 1 and selected[0]["structure_id"] == "first-pose"
    assert selected[0]["outcome"]["primary"]["id"] == first["id"]
    assert selected[0]["outcome"]["structure_count"] == 1


def docking_outcome(database, tools):
    batch(database, "ACD", logical="shared-pipeline")
    index, outcomes = inventory(database)
    outcome = outcomes[0]
    originals = [{**deepcopy(outcome["primary"]), "id": "result-" + str(number), "tool": tool,
                  "structure_ids": ["pose-" + str(number)], "metrics": {"affinity": -number},
                  "evaluation_id": "evaluation-" + str(number)} for number, tool in enumerate(tools)]
    return index, {**outcome, "outcome_kind": "Binding/docking", "activity": "Docking",
                   "primary": originals[-1], "originals": originals, "origin_tools": sorted(set(tools))}


def test_docking_tree_exposes_engines_and_selects_only_their_evaluations_and_poses(database):
    from mn_molecular_design_studio.data_browser import producing_evaluations

    index, outcome = docking_outcome(database, ["vina", "gnina", "udp"])
    nodes, selections = scientific_branches(index, {"id": "compound"}, [outcome])
    assert {node["label"] for node in nodes} == {"Vina", "GNINA", "Uni-Dock Pro"}
    originals = deepcopy(outcome["originals"])
    for node in nodes:
        assert node["count"] == 1 and not node["children"]
        selection = selections[node["id"]]
        assert selection["experiment"] and len(selection["original_ids"]) == 1
        listing = listed_results([outcome], selection=selection)
        assert len(listing) == 1
        selected = listing[0]["outcome"]
        original = next(row for row in outcome["originals"] if row["id"] in selection["original_ids"])
        assert selected["primary"] == original and selected["origin_tools"] == [original["tool"]]
        assert selected["primary"]["metrics"] == original["metrics"]
        assert selected["structure_count"] == 1
        assert [entry["structure_id"] for entry in listed_results([outcome], selection=selection, structures=True)] == original["structure_ids"]
        assert producing_evaluations(index, selected, outcome["originals"]) == [original]
    assert outcome["originals"] == originals and outcome["origin_tools"] == ["gnina", "udp", "vina"]


def test_docking_tool_aliases_share_navigation_without_losing_their_distinct_evidence(database):
    index, outcome = docking_outcome(database, ["udp", "unidock_pro"])
    nodes, selections = scientific_branches(index, {"id": "compound"}, [outcome])
    assert len(nodes) == 1 and nodes[0]["label"] == "Uni-Dock Pro"
    selection = selections[nodes[0]["id"]]
    assert selection["original_ids"] == {row["id"] for row in outcome["originals"]}
    listing = listed_results([outcome], selection=selection)
    assert listing[0]["outcome"]["origin_tools"] == ["udp", "unidock_pro"]
    assert listing[0]["outcome"]["structure_count"] == 2


def test_docking_keeps_separate_experiments_under_each_engine(database):
    index, outcome = docking_outcome(database, ["vina", "vina", "gnina"])
    outcome["originals"][1].update(run_key="another-recorded-run", run_id="mnligand-another-study")
    cofold = {**deepcopy(outcome), "id": "cofold", "activity": "Cofolding", "outcome_kind": "Cofolding",
              "originals": [{**deepcopy(outcome["primary"]), "tool": "boltz2"}], "origin_tools": ["boltz2"]}
    nodes, selections = scientific_branches(index, {"id": "compound"}, [outcome, cofold])
    docking = next(node for node in nodes if node["label"] == "Docking")
    assert docking["count"] == 1
    vina = next(node for node in docking["children"] if node["label"] == "Vina")
    assert vina["count"] == 1 and len(vina["children"]) == 2
    assert len(selections[vina["id"]]["original_ids"]) == 2
    assert listed_results([outcome], selection=selections[vina["id"]])[0]["outcome"]["origin_tools"] == ["vina"]
    assert len({node["label"] for node in vina["children"]}) == 2
    for experiment in vina["children"]:
        assert "import" not in experiment["label"].lower()
        selection = selections[experiment["id"]]
        assert selection["path"][:2] == ["Docking", "Vina"]
        assert selection["experiment"] and len(selection["original_ids"]) == 1
        assert len(listed_results([outcome], selection=selection, structures=True)) == 1
    again = {**outcome, "originals": list(reversed(outcome["originals"]))}
    assert scientific_branches(index, {"id": "compound"}, [again, cofold]) == (nodes, selections)


def test_undated_docking_experiments_have_short_distinct_labels_without_import_wording(database):
    index, outcome = docking_outcome(database, ["vina", "vina"])
    for number, original in enumerate(outcome["originals"]):
        original.update(run_key="study-" + str(number), run_id="mnligand-study-" + str(number))
    nodes, selections = scientific_branches(index, {"id": "compound"}, [outcome])
    assert nodes[0]["label"] == "Vina"
    assert {node["label"] for node in nodes[0]["children"]} == {"Experiment 1", "Experiment 2"}
    assert len(selections) == 3


def test_unknown_docking_engine_is_not_inferred_from_the_import_or_design_name(database):
    index, outcome = docking_outcome(database, [""])
    nodes, selections = scientific_branches(index, {"id": "compound"}, [outcome])
    assert nodes[0]["label"] == "Unknown tool"
    assert len(listed_results([outcome], selection=selections[nodes[0]["id"]])) == 1


def test_experiment_names_use_only_the_evidence_in_that_experiment(database):
    index, outcome = docking_outcome(database, ["vina", "gnina"])
    outcome.update(activity="Cofolding", outcome_kind="Cofolding")
    outcome["originals"][0].update(tool="boltz2")
    outcome["originals"][1].update(tool="optimized_prediction_alphafold3", run_key="different-study", run_id="mnligand-alphafold")
    outcome["origin_tools"] = ["boltz2", "optimized_prediction_alphafold3"]
    nodes, selections = scientific_branches(index, {"id": "compound"}, [outcome])
    assert {node["label"].split(" · ")[0] for node in nodes} == {"Boltz-2", "AlphaFold 3"}
    assert all("import" not in node["label"].lower() for node in nodes)


def test_md_study_replica_navigation_survives_collapsing_single_activity(database):
    from test_md_browser import _add_md

    _add_md(database, replica=1)
    _add_md(database, replica=2)
    index, outcomes = inventory(database)
    dynamics = [outcome for outcome in outcomes if outcome["activity"] == "Molecular dynamics"]
    nodes, selections = scientific_branches(index, {"id": "subject"}, dynamics)
    assert len(nodes) == 1 and nodes[0]["label"] != "Molecular dynamics"
    replicas = [selection for selection in selections.values() if selection.get("md_replica")]
    assert {selection["md_replica"] for selection in replicas} == {"1", "2"}
    assert all(len(listed_results(dynamics, selection=selection)) == 1 for selection in replicas)


def crop(database, sequence, *, project=""):
    return _import(database, name="Tool inputs", method="tool_inputs", kind="inputs", run_id="", project_id=project,
        records=[_record(sequence, name="Exact crop", metadata={"target_identity_type": "sequence"},
                         components=[{"chain": "A", "type": "protein", "sequence": sequence, "role": "target"}])])


def test_reference_family_navigation_preserves_exact_targets_and_evaluation_contexts(database):
    reference = target(database, sequence="AEFGHIK")
    cropped = crop(database, "EFGH")
    batch(database, "ACD", target_id=reference["molecule_ids"][0], logical="full")
    batch(database, "KLM", target_id=cropped["molecule_ids"][0], logical="crop")
    index, outcomes = inventory(database)
    scopes = scientific_scopes(index, [])
    expand_scope(scopes, index, "project:personal")
    children = [scope for scope in scopes.values() if scope["parent"] == "project:personal"]
    assert len(children) == 1 and children[0]["reference_family"]
    family = children[0]
    assert "identity" not in family["subject"] and "sequence" not in family["subject"]
    retained = saved_results(index, subject_results(database, index, family["subject"]))
    assert len(retained) == 2
    assert {outcome["target_id"] for outcome in retained} == {reference["molecule_ids"][0], cropped["molecule_ids"][0]}
    assert family["entity"]["subject_id"] != family["subject"]["id"]


@pytest.mark.parametrize("ambiguous", ["references", "repeated", "overlapping", "pdb", "project"])
def test_reference_families_do_not_fold_ambiguous_or_incompatible_constructs(database, monkeypatch, ambiguous):
    monkeypatch.setattr("mn_molecular_design_studio.auth.accessible_project_ids", lambda user: ["first", "second"])
    reference = target(database, sequence={"repeated": "AEFGEFG", "overlapping": "AEFEFEG"}.get(ambiguous, "AEFGHIK"), project="first")
    cropped = crop(database, "EFE" if ambiguous == "overlapping" else "EFG", project="second" if ambiguous == "project" else "first")
    if ambiguous == "references":
        target(database, sequence="KEFGNQR", name="Another reference", project="first")
    batch(database, "ACD", target_id=reference["molecule_ids"][0], project="first")
    batch(database, "KLM", target_id=cropped["molecule_ids"][0], project="second" if ambiguous == "project" else "first")
    index, outcomes = inventory(database)
    if ambiguous == "pdb":
        index["target_references"] = {reference["molecule_ids"][0]: {"4LNW"}, cropped["molecule_ids"][0]: {"3GWS"}}
    scopes = scientific_scopes(index, [])
    assert not any(scope.get("reference_family") for scope in scopes.values())


def test_logical_design_membership_is_not_split_by_legacy_evidence_run_aliases(database):
    imported = batch(database, "ACD", logical="persisted-logical-run")
    index, outcomes = inventory(database)
    original = outcomes[0]
    extra = {**deepcopy(original["primary"]), "id": "legacy-evidence", "run_id": "older-observation-reference", "run_key": "legacy-run"}
    original["originals"].append(extra)
    nodes, selections = scientific_branches(index, {"id": "subject"}, outcomes)
    assert len(nodes) == 1 and nodes[0]["count"] == 1
    selected = selections[nodes[0]["id"]]
    jobs = source_jobs(index, outcomes, selected)
    assert len(jobs) == 1 and jobs[0]["Evidence runs"] == "older-observation-reference"


def test_ingestion_preserves_explicit_lineage_without_accepting_malformed_ids():
    from mn_molecular_design_studio.molecular_database_ingest import _run_lineage

    metadata = {"studio_run_id": "logical-run", "campaign_id": "campaign", "organizer_node_id": "coordinator",
        "parent_campaign_run_dir": "runs:///design/parent", "pipeline_parent_node_id": {"invalid": True}, "score": 0.9}
    assert _run_lineage(metadata) == {key: metadata[key] for key in (
        "studio_run_id", "campaign_id", "organizer_node_id", "parent_campaign_run_dir")}


def test_portable_campaign_parent_identity_uses_its_recorded_organizer():
    from mn_molecular_design_studio.database_browser import _campaign_reference

    source = {"parent_campaign_run_dir": "runs:///design/parent", "organizer_node_id": "coordinator"}
    assert _campaign_reference(source, "worker-one", "owner", "project") == (
        "coordinator", "design", "parent", "owner", "project")
    assert _campaign_reference(source, "worker-two", "owner", "project") == (
        "coordinator", "design", "parent", "owner", "project")


def test_loaded_experiment_branches_remain_visible_when_selecting_another_target(database):
    reference = target(database)
    batch(database, "ACD", target_id=reference["molecule_ids"][0], logical="experiment")
    index, outcomes = inventory(database)
    scopes = scientific_scopes(index, [])
    expand_scope(scopes, index, "project:personal")
    subject = next(scope for scope in scopes.values() if scope["parent"] == "project:personal")
    subject["loaded"] = True
    branches, selections = scientific_branches(index, subject, outcomes)
    subject["result_branches"] = branches
    nodes, hidden = navigation_tree(scopes, active_scope="project:personal", children=[], root="projects")
    project = next(node for node in nodes[0]["children"] if node["id"] == "project:personal")
    retained = next(node for node in project["children"] if node["id"] == subject["id"])
    assert retained["children"] == branches and not retained["lazy"]


def test_selecting_a_cached_experiment_restores_its_owning_target(database, monkeypatch):
    from streamlit.testing.v1 import AppTest
    from mn_molecular_design_studio import data_browser_ui

    reference = target(database)
    batch(database, "ACD", target_id=reference["molecule_ids"][0], logical="experiment")
    monkeypatch.setenv("MDS_DATA_BROWSER_DATABASE_ROOT", str(database.root))
    monkeypatch.setattr(data_browser_ui, "_render_structure", lambda *args: None)
    event = {}
    monkeypatch.setattr(data_browser_ui, "render_tree", lambda *args, **kwargs: event.copy())
    app = AppTest.from_string("from mn_molecular_design_studio.data_browser_ui import render_data_browser\n"
        f"render_data_browser({USER!r}, [])", default_timeout=30).run()
    assert not app.exception
    scopes = app.session_state["data-browser-scopes"]
    subject = next(scope for scope in scopes.values() if scope["parent"] == "project:personal")
    event.update(selectedId=subject["id"], token="target", expandedIds=[])
    app.run()
    assert not app.exception
    experiment = next(iter(subject["result_selections"]))
    event.update(selectedId="project:personal", token="project")
    app.run()
    assert not app.exception
    event.update(selectedId=experiment, token="experiment")
    app.run()
    assert not app.exception
    assert app.session_state["data-browser-scope"] == subject["id"]
    assert app.session_state["data-browser-folder"] == experiment
