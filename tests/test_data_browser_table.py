from copy import deepcopy
import csv
import io

from streamlit.testing.v1 import AppTest

from mn_molecular_design_studio.data_browser import ReadOnlyLibrary
from mn_molecular_design_studio.data_browser_table import (
    default_columns, evaluation_entries, result_table, table_csv,
)
from mn_molecular_design_studio.database_browser import build_browser_index, subject_results
from test_database_browser import OTHER, USER, _import, _record, database


def entries(database, *, count=1, user=USER):
    imported = _import(database, name="Ligand study", task="docking", run_id="mnligand-study", user=user,
        records=[_record(name="Compound", components=[{"chain": "L", "type": "small_molecule", "smiles": "CCO"}])])
    for number in range(count):
        database.add_observation(imported["group_id"], imported["molecule_ids"][0], "vina", "mnligand-study",
            {"best_score_kcal_mol": -8.5 - number, "zero": 0, "flag": False, "replicate": number, "seed": 1001 + number},
            context={"candidate_id": "compound-" + str(number), "settings": {"exhaustiveness": 32, "seed": 42,
                "use_potentials": False, "engine_options": {"sampling_steps": 200}},
                "row": {"replicate": number, "attempt_index": 0, "seed": 1001 + number}})
    index = build_browser_index(database, user, collection="library", include_related=True)
    index["project_names"] = {}
    raw = subject_results(database, index, {"molecule_ids": set(index["molecules"]), "group_ids": set(index["groups"])})
    rows = []
    for result in raw:
        if not result["observation_id"]:
            continue
        outcome = {**result, "primary": result, "originals": [result], "origin_tools": [result["tool"]],
            "outcome_kind": "Binding/docking", "activity": "Docking", "eligibility": "Imported", "target_name": "TRα · 4LNW"}
        rows.append({"id": result["id"], "outcome": outcome, "result": result, "structure_id": ""})
    return index, rows


def test_table_and_csv_expose_exact_constructs_without_combining_target_sequences(database):
    index, selected = entries(database, count=2)
    for entry, identifier, sequence in zip(selected, ["target-one", "target-two"], ["EFGH", "EFGHK"]):
        index["molecules"][identifier] = {"name": "TRβ", "target_names": set(), "import_names": set(), "group_ids": set(),
            "components": [{"chain": "A", "type": "protein", "sequence": sequence}]}
        entry["outcome"]["target_id"] = identifier
    rows = result_table(database, index, selected)
    assert [row["Target construct"] for row in rows] == ["TRβ · 4 aa", "TRβ · 5 aa"]
    assert [row["Target sequence"] for row in rows] == ["A: EFGH", "A: EFGHK"]
    assert "Target construct" in default_columns(rows)
    exported = list(csv.DictReader(io.StringIO(table_csv(rows).decode("utf-8-sig"))))
    assert [row["Target ID"] for row in exported] == ["target-one", "target-two"]
    assert [row["Target sequence"] for row in exported] == ["A: EFGH", "A: EFGHK"]


def test_unresolved_target_does_not_borrow_a_construct_from_another_row(database):
    index, selected = entries(database)
    selected[0]["outcome"].update(target_id="not-visible", target_name="TRβ")
    row = result_table(database, index, selected)[0]
    assert row["Target"] == "TRβ"
    assert row["Target construct"] is None and row["Target sequence"] is None


def test_design_table_and_csv_distinguish_native_from_final_acceptance(database):
    index, selected = entries(database, count=2)
    for entry, level in zip(selected, ["Native acceptance", "Final acceptance"]):
        entry["outcome"].update(outcome_kind="Design", activity="Protein design", eligibility="Passed")
        entry["result"]["acceptance_level"] = level
    rows = result_table(database, index, selected)
    assert [row["Admission"] for row in rows] == ["Passed", "Passed"]
    assert [row["Admission basis"] for row in rows] == ["Native acceptance", "Final acceptance"]
    assert "Admission basis" in default_columns(rows)
    exported = list(csv.DictReader(io.StringIO(table_csv(rows).decode("utf-8-sig"))))
    assert [row["Admission basis"] for row in exported] == ["Native acceptance", "Final acceptance"]


def test_table_construct_name_uses_the_evaluated_receptor_not_a_shared_sequence_alias(database):
    index, selected = entries(database)
    index["molecules"]["target"] = {"name": "TRβ · 6KKB", "target_names": set(), "import_names": set(), "group_ids": set(),
        "components": [{"chain": "A", "type": "protein", "sequence": "EFGH"}]}
    selected[0]["outcome"].update(target_id="target", target_reference_label="TRβ · 3GWS · 4 aa")
    row = result_table(database, index, selected)[0]
    assert row["Target construct"] == "TRβ · 3GWS · 4 aa" and row["Target sequence"] == "A: EFGH"


def test_table_hydrates_only_exact_evaluation_preserves_native_values_and_missing(database):
    index, selected = entries(database, count=2)
    assert all(not entry["result"]["metrics"] for entry in selected)
    before = database.path.read_bytes()
    rows = result_table(ReadOnlyLibrary(database.root), index, selected)
    assert len(rows) == 2
    for entry, row in zip(selected, rows):
        replicate = row["Replicate"]
        assert row["Result · best_score_kcal_mol"] == -8.5 - replicate
        assert row["Input · seed"] == 1001 + replicate
        assert row["Input · requested_seed"] == 42
        assert row["Result · zero"] == 0 and row["Result · flag"] is False
        assert row["Input · use_potentials"] is False
        assert row["Input · engine_options.sampling_steps"] == 200
        assert row["Attempt"] == 0 and row["Model / pose"] is None
        assert row["Observation ID"] == entry["result"]["observation_id"]
        assert row["Project"] == "Personal" and row["Method"] == "Vina"
    assert database.path.read_bytes() == before


def test_table_uses_scoped_project_names_and_does_not_read_private_evidence(database):
    index, selected = entries(database)
    _, private = entries(database, user=OTHER)
    group = index["groups"][selected[0]["result"]["group_id"]]
    group["project_id"] = "project"
    index["project_names"]["project"] = "Renamed project"
    rows = result_table(database, index, [*selected, *private])
    assert len(rows) == 1 and rows[0]["Project"] == "Renamed project"


def test_all_evaluations_csv_keeps_replicates_separate_and_never_expands_attached_analysis(database):
    index, selected = entries(database, count=2)
    combined = deepcopy(selected[0])
    combined["outcome"]["originals"] = [entry["result"] for entry in selected]
    combined["outcome"]["attached"] = [{"id": "post-analysis", "tool": "biotite"}]
    expanded = evaluation_entries([combined, combined])
    assert len(expanded) == 2
    rows = result_table(database, index, expanded)
    assert {row["Input · seed"] for row in rows} == {1001, 1002}
    assert {row["Result · best_score_kcal_mol"] for row in rows} == {-8.5, -9.5}
    assert all(row["Evaluations"] == 2 and row["Method"] == "Vina" for row in rows)


def test_csv_covers_every_page_and_contains_hidden_fields_and_safe_strings(database):
    index, selected = entries(database, count=53)
    rows = result_table(database, index, selected)
    rows[0]["Name"] = '=HYPERLINK("bad")'
    rows[0]["Input · unsafe"] = "\t=1+1"
    rows[0]["Result · absent"] = None
    output = table_csv(rows)
    parsed = list(csv.DictReader(io.StringIO(output.decode("utf-8-sig"))))
    assert len(parsed) == 53
    assert len({row["Observation ID"] for row in parsed}) == 53
    assert parsed[0]["Name"].startswith("'=")
    assert parsed[0]["Input · unsafe"].startswith("'\t=")
    assert parsed[0]["Result · absent"] == ""
    assert float(parsed[0]["Result · best_score_kcal_mol"]) < 0
    assert "Observation ID" not in default_columns(rows)
    assert {"Name", "Project", "Method", "Target", "Result · best_score_kcal_mol", "Input · seed"} <= set(default_columns(rows))


def test_missing_observation_does_not_borrow_other_evaluation_metrics(database):
    index, selected = entries(database, count=2)
    selected[0]["result"]["observation_id"] = "missing"
    row = result_table(database, index, selected[:1])[0]
    assert not any(column.startswith("Result · ") for column in row)


def test_admission_measurements_are_not_mislabelled_as_original_evaluation(database):
    metadata = {"source_candidate_id": "candidate", "validation_acceptance": {"measurements": {
        "accepted": {"status": "measured", "value": 0, "source_tool": "generator", "measurement_candidate_id": "candidate"},
        "wrong": {"status": "measured", "value": 99, "source_tool": "generator", "measurement_candidate_id": "other"}}}}
    imported = _import(database, method="user_uploads", kind="curated", records=[_record(metadata=metadata)])
    index = build_browser_index(database, USER, collection="library", include_related=True)
    result = subject_results(database, index, {"molecule_ids": set(index["molecules"]), "group_ids": set(index["groups"])})[0]
    entry = {"outcome": {**result, "id": "saved-design", "originals": [result], "origin_tools": ["generator"],
        "activity": "Protein design", "outcome_kind": "Design", "eligibility": "Passed"}, "result": result, "structure_id": ""}
    row = result_table(database, index, [entry])[0]
    assert row["Admission evidence · generator · accepted"] == "0"
    assert not any("wrong" in column or column.startswith("Result · ") for column in row)


def test_table_controls_export_full_selection_not_visible_page_and_invalidate_stale_csv(database):
    index, selected = entries(database, count=53)
    source = ("from mn_molecular_design_studio.data_browser_ui import _table_controls\n"
        "from mn_molecular_design_studio.data_browser import ReadOnlyLibrary\n"
        f"index={index!r}\nentries={selected!r}\n"
        f"_table_controls(ReadOnlyLibrary(__import__('pathlib').Path({str(database.root)!r})), "
        "index, entries, entries[:50], 'snapshot', 'selection', 1)")
    app = AppTest.from_string(source, default_timeout=30).run()
    assert not app.exception
    app.button(key="data-browser-prepare-csv").click().run()
    assert not app.exception
    prepared = app.session_state["data-browser-prepared-csv"]
    assert prepared[2] == 53
    assert len(list(csv.DictReader(io.StringIO(prepared[1].decode("utf-8-sig"))))) == 53
    app.selectbox(key="data-browser-column-preset").set_value("Compact").run()
    assert not app.exception
    assert app.session_state["data-browser-prepared-csv"][2] == 53
    app.selectbox(key="data-browser-csv-mode").set_value("All original evaluations").run()
    assert not app.exception
    assert "data-browser-prepared-csv" not in app.session_state
