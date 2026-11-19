from copy import deepcopy
import json

import pytest

from mn_molecular_design_studio.data_browser import ReadOnlyLibrary, listed_results, saved_results
from mn_molecular_design_studio.data_browser_table import default_columns, result_table, table_csv, table_profile
from mn_molecular_design_studio.database_browser import build_browser_index, subject_results
from test_database_browser import OTHER, USER, _import, database
from test_database_design_evaluations import _analysis, design
from test_database_outcomes import _pass_design
from test_md_browser import _add_md


def listing(database, design):
    _pass_design(database, design)
    index = build_browser_index(database, USER, collection="library", include_related=True)
    raw = subject_results(database, index, {"molecule_ids": set(index["molecules"]), "group_ids": set(index["groups"])})
    entries = listed_results([outcome for outcome in saved_results(index, raw) if outcome["outcome_kind"] == "Design"])
    assert len(entries) == 1
    return index, entries


def test_original_and_three_same_structure_analyses_have_separate_numeric_columns(database, design):
    identifiers = {}
    for tool, metrics in (("ipsae", {"ipsae_min": 0}), ("biotite_interface", {"interface_contact_pairs": 12}),
                          ("pep_patch", {"peppatch_total_area": 14.2})):
        identifiers[tool] = _analysis(database, design, tool, run=tool, metrics=metrics)
    index, entries = listing(database, design)
    before = database.path.read_bytes()
    row = result_table(ReadOnlyLibrary(database.root), index, entries)[0]
    assert row["Result · confidence"] == 0.9
    assert row["Analysis · ipSAE · Result · ipsae_min"] == 0
    assert row["Analysis · Biotite · Result · interface_contact_pairs"] == 12
    assert row["Analysis · PEP-Patch · Result · peppatch_total_area"] == 14.2
    assert all(row["Analysis · " + tool + " · Status"] == "Recorded" for tool in ("ipSAE", "Biotite", "PEP-Patch"))
    assert row["Analysis · ipSAE · Engine version"] == "v2"
    recorded = json.loads(row["Analysis · ipSAE · Records"])
    assert recorded[0]["observation_id"] == identifiers["ipsae"]
    assert recorded[0]["evaluation_id"] == row["Analysis · ipSAE · Scored evaluation ID"]
    assert database.path.read_bytes() == before
    assert "Analysis · ipSAE · Result · ipsae_min" in default_columns([row])


@pytest.mark.parametrize("context,expected", [
    ({"source_structure_id": "other", "representation_id": "other", "source_structure_digest": "b" * 64,
      "structure_digest": "b" * 64}, "Other model"),
    ({"structure_digest": "b" * 64}, "Derived from selected structure"),
    ({"status": "failed"}, "Failed / unavailable"),
    ({"status": "pending"}, "Status unknown / incomplete"),
    ({"target_molecule_id": ""}, "Target link unverified"),
    ({"source_reference_kind": "sequence", "source_structure_id": "", "representation_id": "",
      "source_structure_digest": "", "structure_digest": ""}, "Sequence only"),
])
def test_unqualified_analysis_is_labelled_without_supplying_numeric_score(database, design, context, expected):
    _analysis(database, design, "ipsae", metrics={"ipsae_min": 0.85}, context=context)
    index, entries = listing(database, design)
    row = result_table(database, index, entries)[0]
    assert row["Analysis · ipSAE · Status"] == expected
    assert "Analysis · ipSAE · Result · ipsae_min" not in row
    assert json.loads(row["Analysis · ipSAE · Records"])[0]["metrics"]["ipsae_min"] == 0.85


def test_conflicting_declared_chain_selection_is_not_scored_as_the_same_selection(database, design):
    _analysis(database, design, "ipsae", context={"binder_chains": ["C"]})
    index, entries = listing(database, design)
    entries[0]["result"]["evaluation"]["binder_chains"] = ["B"]
    row = result_table(database, index, entries)[0]
    assert row["Analysis · ipSAE · Status"] == "Different component selection"
    assert "Analysis · ipSAE · Scored evaluation ID" not in row


def test_multiple_protocols_are_not_silently_ranked_or_replaced_by_the_latest(database, design):
    for number, cutoff in enumerate((8, 12)):
        _analysis(database, design, "ipsae", run="analysis-" + str(number),
            metrics={"ipsae_min": cutoff / 20}, context={"settings": {"cutoff": cutoff}})
    index, entries = listing(database, design)
    row = result_table(database, index, entries)[0]
    assert row["Analysis · ipSAE · Status"] == "Multiple evaluations (2)"
    assert "Analysis · ipSAE · Result · ipsae_min" not in row
    assert len(json.loads(row["Analysis · ipSAE · Records"])) == 2
    assert b"cutoff" in table_csv([row])


def test_one_completed_and_one_failed_analysis_retains_both_with_one_scored_source(database, design):
    _analysis(database, design, "ipsae", run="valid", metrics={"ipsae_min": 0.5})
    _analysis(database, design, "ipsae", run="failed", metrics={"ipsae_min": 0.9}, context={"status": "failed"})
    index, entries = listing(database, design)
    row = result_table(database, index, entries)[0]
    assert row["Analysis · ipSAE · Result · ipsae_min"] == 0.5
    assert row["Analysis · ipSAE · Evaluations"] == 2
    assert row["Analysis · ipSAE · Status"] == "Recorded; Failed / unavailable (1)"


def test_analysis_does_not_leak_between_targets_or_private_collections(database, design):
    _analysis(database, design, "ipsae", context={"target_molecule_id": "other-target"})
    private = _import(database, user=OTHER, task="analysis", kind="evaluation", records=[design["record"]])
    _analysis(database, {**design, "group_id": private["group_id"]}, "plip", metrics={"interaction_count": 99})
    index, entries = listing(database, design)
    row = result_table(database, index, entries)[0]
    assert row["Post-analysis"] == "Not recorded"
    assert {field: value for field, value in row.items() if field.startswith("Analysis · ")} == {"Analysis · BinderScope · Status": "Not assessed"}


def test_validation_refolding_and_filter_decisions_are_not_post_analysis(database, design):
    _analysis(database, design, "ipsae", context={"workflow_role": "design_validation"})
    _analysis(database, design, "validation_acceptance")
    index, entries = listing(database, design)
    assert result_table(database, index, entries)[0]["Post-analysis"] == "Not recorded"


def test_contradictory_source_digest_is_unverified_not_a_scored_attachment(database, design):
    _analysis(database, design, "ipsae", context={"source_structure_digest": "b" * 64, "structure_digest": "b" * 64})
    index, entries = listing(database, design)
    row = result_table(database, index, entries)[0]
    assert row["Post-analysis"] == "Unverified links (1)"
    assert {field: value for field, value in row.items() if field.startswith("Analysis · ")} == {"Analysis · BinderScope · Status": "Not assessed"}


def test_library_input_rows_do_not_acquire_analysis_from_nearby_design_rows(database, design):
    _analysis(database, design, "ipsae")
    index, entries = listing(database, design)
    reference = deepcopy(entries[0])
    reference["outcome"]["outcome_kind"] = "Reference molecule"
    rows = result_table(database, index, [*entries, reference])
    assert rows[0]["Analysis · ipSAE · Status"] == "Recorded"
    assert rows[1]["Post-analysis"] == "Not applicable"
    assert not any(field.startswith("Analysis · ") for field in rows[1])


def test_structure_listing_does_not_attach_primary_pose_analysis_to_an_alternate_pose(database, design):
    _analysis(database, design, "ipsae")
    index, entries = listing(database, design)
    entries[0]["result"]["structure_ids"].append("alternate-pose")
    entries[0]["structure_id"] = "alternate-pose"
    row = result_table(database, index, entries)[0]
    assert row["Analysis · ipSAE · Status"] == "Other model"
    assert "Analysis · ipSAE · Scored evaluation ID" not in row


def test_batch_analysis_loading_uses_one_cohort_query_not_one_query_per_entry(database, design, monkeypatch):
    from mn_molecular_design_studio import data_browser_analyses

    _analysis(database, design, "ipsae")
    index, entries = listing(database, design)
    calls = []
    original = data_browser_analyses.subject_results

    def counted(*args, **kwargs):
        calls.append(True)
        return original(*args, **kwargs)

    monkeypatch.setattr(data_browser_analyses, "subject_results", counted)
    rows = result_table(database, index, entries * 50)
    assert len(rows) == 50 and len(calls) == 1


@pytest.mark.parametrize("activity,preset,metric,setting", [
    ("Protein design", "Design", "oracle_ppi.iptm", "binder_length"),
    ("Docking", "Docking", "best_score_kcal_mol", "exhaustiveness"),
    ("Cofolding", "Structure prediction", "iptm", "num_recycles"),
    ("Structure prediction", "Structure prediction", "ranking_score", "samples"),
    ("Binding affinity", "Binding affinity", "affinity_pred_value", "seed"),
    ("Molecular dynamics", "Molecular dynamics", "production_actual_duration_ns", "production_ns"),
])
def test_activity_profiles_expose_relevant_native_results_and_settings(activity, preset, metric, setting):
    row = {"Name": "Molecule", "Method": "Engine", "Activity": activity, "Project": "Study", "Admission": "Recorded",
           "Result · " + metric: 0, "Input · " + setting: 0, "Result · irrelevant": 9, "Post-analysis": "Not recorded"}
    assert table_profile([row]) == preset
    assert "Result · " + metric in default_columns([row])
    assert "Input · " + setting in default_columns([row])
    assert "Result · irrelevant" not in default_columns([row])
    assert "Admission" in default_columns([row]) if preset == "Design" else True
    assert "Result · " + metric not in default_columns([row], "Compact")
    assert "Result · irrelevant" in default_columns([row], "All columns")


def test_mixed_activity_folders_have_overview_not_a_misleading_single_activity_preset():
    assert table_profile([{"Activity": "Docking"}, {"Activity": "Protein design"}]) == "Overview"


@pytest.mark.parametrize("boundary", [None, "replica", "segment", "study"])
def test_md_analysis_requires_exact_recorded_segment_replica_and_study(database, boundary):
    production = _add_md(database)
    analysis = _add_md(database, namespace="md.analysis", run="analysis", replica=2 if boundary == "replica" else 1,
                       study="other-study" if boundary == "study" else "study-one")
    if boundary == "segment":
        with database.connection() as connection:
            connection.execute("UPDATE molecule_observations SET context=json_set(context,'$.source_production_run_id','other-segment') WHERE id=?", (analysis,))
    index = build_browser_index(database, USER, collection="library", include_related=True)
    raw = subject_results(database, index, {"molecule_ids": set(index["molecules"]), "group_ids": set(index["groups"])})
    original = next(result for result in raw if result["observation_id"] == production)
    outcome = next(outcome for outcome in saved_results(index, raw) if any(row["observation_id"] == production for row in outcome["originals"]))
    entry = {"outcome": outcome, "result": original, "structure_id": ""}
    row = result_table(database, index, [entry])[0]
    if boundary:
        assert row["Post-analysis"] == "Not recorded"
    else:
        assert row["Analysis · MD Analysis · Status"] == "Recorded MD segment"
        assert row["Analysis · MD Analysis · Result · production_actual_duration_ns"] == 99
