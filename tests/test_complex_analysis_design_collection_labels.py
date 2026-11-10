from mn_molecular_design_studio.complex_analysis_collection_ui import (
    _database_analysis_cohorts,
    _design_collection_label,
    _is_design_collection,
    _resolve_collection_target_labels,
)


def test_design_collection_label_uses_campaign_target_and_source_host():
    label = _design_collection_label({
        "id": "collection-1",
        "name": "R3 (31 aa) · accepted RFD4",
        "origin": {
            "collection_kind": "curated",
            "campaign_name": "Tau R3 native campaign",
            "run_id": "20261030-070659-2044e3d6",
            "task_group": "design",
            "workstation": "PC004033",
            "tool": "rfd4_proteina",
            "design_target_name": "R3 (31 aa)",
        },
    })

    assert label["target"] == "R3 (31 aa)"
    assert label["run"] == (
        "Tau R3 native campaign · design/20261030-070659-2044e3d6 · PC004033"
    )
    assert label["generator"] == "rfd4_proteina"


def test_design_collection_label_does_not_call_analysis_tool_a_generator():
    label = _design_collection_label({
        "name": "Acceptance evidence",
        "origin": {
            "collection_kind": "generated",
            "tool": "validation_acceptance",
            "run_id": "20261031-085246-eb53803c",
            "task_group": "design",
            "design_target_molecule_id": "mol-target-123",
        },
    })

    assert label["generator"] == "Not recorded"
    assert label["target"] == "Target record mol-target-123"
    assert label["run"].startswith("Campaign not recorded · design/")


def test_job_collection_name_is_used_when_campaign_name_was_not_normalized():
    label = _design_collection_label({
        "name": "PD-L1 distributed complete BindCraft2 acceptance",
        "origin": {
            "type": "job",
            "collection_kind": "generated",
            "tool": "bindcraft2_campaign",
            "run_id": "20261029-094457-fb2fd80a",
            "task_group": "design",
            "workstation": "pc002975",
        },
    })

    assert label["run"] == (
        "PD-L1 distributed complete BindCraft2 acceptance · "
        "design/20261029-094457-fb2fd80a · pc002975"
    )
    assert label["generator"] == "bindcraft2_campaign"


def test_design_collection_picker_keeps_explicit_design_collections():
    assert not _is_design_collection({"collection_kind": "reference", "origin": {}})
    assert _is_design_collection({"collection_kind": "generated", "origin": {"run_id": "design-1"}})
    assert _is_design_collection({"origin": {"collection_kind": "curated", "saved_selection": True}})


def test_design_collection_picker_hides_unlinked_utility_groups_but_keeps_run_cohorts():
    assert not _is_design_collection({
        "name": "Docking targets",
        "origin": {"type": "user", "collection_kind": "curated", "purpose": "Reusable docking receptors"},
    })
    assert _is_design_collection({
        "name": "Distributed BindCraft2",
        "origin": {"type": "job", "collection_kind": "generated", "run_id": "20261029-094457-fb2fd80a"},
    })


def test_analysis_cohort_groups_cross_engine_source_runs_by_exact_target_not_names():
    same_target = "a" * 64
    groups = [
        {"id": "rfd4-run", "project_id": "project-1", "owner_id": "owner",
         "origin": {"type": "job", "collection_kind": "generated", "tool": "rfd4_proteina",
                     "run_id": "run-rfd4", "task_group": "design", "workstation": "PC004033",
                     "design_target_molecule_id": "target-1", "design_target_sequence_sha256": same_target},
         "summary": {"target": "PD-L1", "generator": "RFD4"}},
        {"id": "bindcraft-run", "project_id": "project-1", "owner_id": "owner",
         "origin": {"type": "job", "collection_kind": "generated", "tool": "bindcraft2",
                     "run_id": "run-bindcraft", "task_group": "design", "workstation": "pc002975",
                     "design_target_molecule_id": "target-1", "design_target_sequence_sha256": same_target},
         "summary": {"target": "PD-L1 construct A", "generator": "BindCraft2"}},
        {"id": "different-construct", "project_id": "project-1", "owner_id": "owner",
         "origin": {"type": "job", "collection_kind": "generated", "tool": "bindcraft2",
                     "run_id": "run-other", "task_group": "design", "workstation": "pc002975",
                     "design_target_molecule_id": "target-1", "design_target_sequence_sha256": "b" * 64},
         "summary": {"target": "PD-L1", "generator": "BindCraft2"}},
        # A matching display name alone does not make an unknown target compatible.
        {"id": "unknown-target", "project_id": "project-1", "owner_id": "owner",
         "origin": {"type": "job", "collection_kind": "generated", "tool": "rfd4_proteina",
                     "run_id": "run-unknown", "task_group": "design", "workstation": "PC004032"},
         "summary": {"target": "PD-L1", "generator": "RFD4"}},
    ]

    cohorts = _database_analysis_cohorts(groups)

    assert len(cohorts) == 3
    compatible = next(row for row in cohorts if set(row["group_ids"]) == {"rfd4-run", "bindcraft-run"})
    assert compatible["cohort_type"] == "Cross-run comparison cohort"
    assert compatible["source_runs"] == 2
    assert compatible["workstations"] == "PC004033, pc002975"
    assert compatible["generators"] == "BindCraft2, RFD4"
    assert compatible["campaign_ids"] == []


def test_analysis_cohort_keeps_unlinked_unknown_targets_distinct_and_names_them():
    cohorts = _database_analysis_cohorts([
        {"id": "unknown-a", "name": "Design A", "project_id": "p",
         "origin": {"type": "job", "tool": "rfd4_proteina", "run_id": "run-a"},
         "summary": {"target": "Not recorded", "collection": "Design A", "generator": "RFD4"}},
        {"id": "unknown-b", "name": "Design B", "project_id": "p",
         "origin": {"type": "job", "tool": "rfd4_proteina", "run_id": "run-b"},
         "summary": {"target": "Not recorded", "collection": "Design B", "generator": "RFD4"}},
    ])

    assert len(cohorts) == 2
    assert {row["cohort_type"] for row in cohorts} == {"Source group · target identity missing"}
    assert {row["name"] for row in cohorts} == {"Design A", "Design B"}


def test_target_wide_cohort_is_not_named_as_one_partially_linked_campaign():
    groups = [{"id": key, "project_id": "p", "origin": {
        "tool": "rfd4_proteina", "target_molecule_id": "target", "run_id": key,
        **({"campaign_id": "campaign-a", "campaign_name": "200 passing designs"} if key == "a" else {}),
    }, "summary": {"target": "Tau R3", "generator": "RFD4"}} for key in ("a", "b")]
    cohort, = _database_analysis_cohorts(groups)
    assert cohort["cohort_type"] == "Cross-run comparison cohort"
    assert cohort["name"] != "200 passing designs"
    assert cohort["campaign_ids"] == ["campaign-a"]


def test_target_display_names_use_authorized_metadata_without_changing_identity():
    groups = [{"origin": {"target_molecule_id": "target"},
               "summary": {"target": "Target record target"}},
              {"origin": {"target_molecule_id": "inaccessible"},
               "summary": {"target": "Target record inaccessible"}}]
    user = {"id": "student"}

    class Database:
        def list_structure_records(self, selected_user, **kwargs):
            assert selected_user is user
            assert kwargs["molecule_ids"] == ["target", "inaccessible"]
            return {"items": [{"id": "target", "name": "Tau R3 (31 aa)"}]}

    _resolve_collection_target_labels(Database(), user, groups)
    assert groups[0]["summary"]["target"] == "Tau R3 (31 aa) · target"
    assert groups[0]["origin"]["target_molecule_id"] == "target"
    assert groups[1]["summary"]["target"] == "Target record inaccessible"


def test_scientific_result_task_groups_show_the_stage_of_the_source_work():
    cohorts = _database_analysis_cohorts([
        {"id": "fold", "project_id": "p", "origin": {
            "type": "scientific_results", "task_group": "refolding-validation", "run_id": "fold-run"},
            "summary": {"target": "Not recorded", "collection": "Fold output"}},
        {"id": "dock", "project_id": "p", "origin": {
            "type": "scientific_results", "task_group": "docking", "run_id": "dock-run"},
            "summary": {"target": "Not recorded", "collection": "Dock output"}},
    ])

    assert {tuple(row["group_ids"]): row["activity"] for row in cohorts} == {
        ("fold",): "Folding and prediction",
        ("dock",): "Ligand docking",
    }


def test_cohorts_with_missing_target_names_show_their_exact_target_identity():
    cohorts = _database_analysis_cohorts([
        {"id": "target-a", "name": "RFD4", "project_id": "p",
         "origin": {"type": "job", "tool": "rfd4_proteina", "target_identity": "Not recorded",
                    "target_molecule_id": "mol-a", "target_sequence_sha256": "a" * 64},
         "summary": {"target": "Not recorded", "generator": "RFD4"}},
        {"id": "target-b", "name": "RFD4", "project_id": "p",
         "origin": {"type": "job", "tool": "rfd4_proteina", "target_identity": "Not recorded",
                    "target_molecule_id": "mol-b", "target_sequence_sha256": "b" * 64},
         "summary": {"target": "Not recorded", "generator": "RFD4"}},
    ])

    assert len(cohorts) == 2
    assert {row["target"] for row in cohorts} == {"Target record mol-a", "Target record mol-b"}
    assert len({row["name"] for row in cohorts}) == 2


def test_analysis_cohort_table_has_one_row_for_a_multi_run_target_cohort(monkeypatch):
    from streamlit.testing.v1 import AppTest

    class Database:
        def groups(self, user, *, collection):
            return [
                {"id": "rfd4", "name": "RFD4 candidate", "project_id": "project-1", "owner_id": user["id"],
                 "collection_kind": "generated", "updated_at": "2026-11-10", "origin": {
                     "type": "job", "collection_kind": "generated", "tool": "rfd4_proteina",
                     "run_id": "run-rfd4", "task_group": "design", "workstation": "PC004033",
                     "design_target_molecule_id": "target-1", "design_target_sequence_sha256": "a" * 64}},
                {"id": "bindcraft", "name": "BindCraft candidate", "project_id": "project-1", "owner_id": user["id"],
                 "collection_kind": "generated", "updated_at": "2026-11-09", "origin": {
                     "type": "job", "collection_kind": "generated", "tool": "bindcraft2",
                     "run_id": "run-bindcraft", "task_group": "design", "workstation": "pc002975",
                     "design_target_molecule_id": "target-1", "design_target_sequence_sha256": "a" * 64}},
            ]

        def collection_statistics(self, user, group_ids):
            return {group_id: {"molecule_count": 1, "structure_count": 1, "evaluated_design_count": 0}
                    for group_id in group_ids}

        def cohort_statistics(self, user, cohorts):
            return {key: {"molecule_count": 2, "structure_count": 2, "evaluated_design_count": 0}
                    for key in cohorts}

    monkeypatch.setattr("mn_molecular_design_studio.molecular_database.MolecularDatabase", Database)
    user = {"id": "student", "role": "user"}
    app = AppTest.from_string(f'''
from mn_molecular_design_studio.complex_analysis_collection_ui import _database_collection_samples
_database_collection_samples({user!r})
''', default_timeout=10).run()

    assert not app.exception
    assert len(app.dataframe) == 1
    rows = app.dataframe[0].value
    assert len(rows) == 1
    assert rows.iloc[0]["DB-retained designs"] == 2
    assert rows.iloc[0]["All-run totals"] == "Not part of retained Database cohort"
    assert rows.iloc[0]["Source runs"] == "2 source run(s) · 2 Database group(s)"
    assert not any(widget.label == "Design collection(s)" for widget in app.multiselect)


def test_analysis_cohort_table_renders_every_distinct_work_group(monkeypatch):
    from streamlit.testing.v1 import AppTest

    groups = [
        {"id": "fold", "name": "Fold run", "project_id": "", "owner_id": "student",
         "collection_kind": "generated", "origin": {"type": "scientific_results",
             "task_group": "refolding-validation", "run_id": "fold-run"}},
        {"id": "dock", "name": "Dock run", "project_id": "", "owner_id": "student",
         "collection_kind": "generated", "origin": {"type": "scientific_results",
             "task_group": "docking", "run_id": "dock-run"}},
    ]

    class Database:
        def groups(self, user, *, collection):
            return groups

        def collection_statistics(self, user, group_ids):
            return {group_id: {"molecule_count": 1, "structure_count": 1, "evaluated_design_count": 0}
                    for group_id in group_ids}

        def cohort_statistics(self, user, cohorts):
            return {cohort_id: {"molecule_count": 1, "structure_count": 1,
                                "evaluated_design_count": 0,
                                "validation_acceptance": {"source_memberships": 0}}
                    for cohort_id in cohorts}

    monkeypatch.setattr("mn_molecular_design_studio.molecular_database.MolecularDatabase", Database)
    app = AppTest.from_string('''
from mn_molecular_design_studio.complex_analysis_collection_ui import _database_collection_samples
_database_collection_samples({"id": "student", "role": "user"})
''', default_timeout=10).run()

    assert not app.exception
    rows = app.dataframe[0].value
    assert set(rows["Work group"]) == {"Fold run", "Dock run"}
    assert set(rows["Scientific work"]) == {"Folding and prediction", "Ligand docking"}

    app.session_state["complex-analysis-collections-open"] = True
    app.selectbox(key="complex-analysis-work-type").select("Ligand docking").run()
    assert not app.exception
    assert not app.session_state["complex-analysis-collections-open"]
    assert list(app.dataframe[0].value["Work group"]) == ["Dock run"]
    app.text_input(key="complex-analysis-work-search").input("fold-run").run()
    assert not app.exception
    assert any("No saved work matches" in item.value for item in app.info)


def test_unlinked_analysis_outputs_are_not_listed_as_design_collections():
    assert not _is_design_collection({
        "collection_kind": "generated",
        "origin": {"type": "job", "tool": "biotite_interface", "run_id": "analysis-1"},
    })
    assert _is_design_collection({
        "collection_kind": "curated",
        "origin": {"type": "job", "tool": "plip", "saved_selection": True, "run_id": "pose-cohort"},
    })
    assert not _is_design_collection({
        "collection_kind": "generated",
        "origin": {"type": "scientific_results", "tool": "plip", "run_id": "analysis-2"},
    })
