import hashlib
import json

import pytest
from streamlit.testing.v1 import AppTest

from mn_molecular_design_studio.data_browser import branch_nodes, folder_scopes, listed_results, saved_results
from mn_molecular_design_studio.database_browser import build_browser_index, subject_results
from mn_molecular_design_studio.md_browser import (
    MAX_ARTIFACT_BYTES, MDArtifactReader, cache_key, ligand_atom_block, matching_study, refresh_requests, report_choices, snapshot_references, study_observations,
)
from mn_molecular_design_studio.md_results_ui import _render_endpoint_details, _render_structural_dynamics
from test_database_browser import OTHER, USER, _import, database
from test_md_results_ui import FakeStreamlit


NODE = "11111111-1111-4111-8111-111111111111"


def test_refresh_polls_pending_transfers_instead_of_starting_duplicate_work():
    state = {"study-pending": "pending", "study-failed": "failed", "study-done": "completed", "another-study": "other"}
    refresh_requests(state, "study-", lambda command: {"status": command})
    assert state == {"study-pending": "pending", "study-done": "completed", "another-study": "other"}


def test_atom_name_source_uses_exact_segment_or_proven_shared_preparation_not_another_replica():
    from test_md_interaction_diagram import _ligand_pdb

    class Reader:
        def __init__(self):
            self.calls = []

        def text(self, path, description, *, expected):
            self.calls.append((path, expected))
            return _ligand_pdb()

    reader = Reader()
    path = "artifacts/legacy/prep-one/prep-one/topology.pdb"
    study = {"production_runs": [{"run_id": "prod-one", "topology_path": path, "topology_sha256": "a" * 64}]}
    report = {"selected_ligand": {"resname": "LIG", "chain": "Z", "atom_count": 0},
        "md_result": {"output_files": {"system_pdb": "/prepared-system/prep-one/topology.pdb"}}}
    block, note = ligand_atom_block(reader, study, report, "continued-prod-two")
    assert "HETATM" in block and "digest-verified" in note
    assert reader.calls == [(path, "a" * 64)]
    report["md_result"]["output_files"]["system_pdb"] = "/prepared-system/other-preparation/topology.pdb"
    assert ligand_atom_block(reader, study, report, "continued-prod-two")[0] == ""
    assert len(reader.calls) == 1
    with pytest.raises(ValueError, match="contradicts"):
        ligand_atom_block(reader, study, report, "prod-one")
    report["md_result"]["output_files"] = {}
    assert ligand_atom_block(reader, study, report, "prod-one")[0]


def test_inline_atom_source_is_not_silently_replaced_and_ambiguous_topology_selection_fails_closed():
    from test_md_interaction_diagram import _ligand_pdb

    class Reader:
        def text(self, *args, **kwargs):
            return _ligand_pdb() + _ligand_pdb().replace("LIG L   1", "LIG B   1")

    report = {"ligand_pdb": "saved atom block"}
    assert ligand_atom_block(Reader(), None, report, "prod")[0] == "saved atom block"
    study = {"production_runs": [{"run_id": "prod", "topology_path": "artifacts/top.pdb", "topology_sha256": "a" * 64}]}
    report = {"selected_ligand": {"resname": "LIG", "atom_count": 0}}
    assert ligand_atom_block(Reader(), study, report, "prod")[0] == ""


def _add_md(database, *, study="study-one", pdb="4LNW", run="import-md", replica=1, namespace="md.production", node=NODE):
    _import(database, name="Library molecule", run_id="library-import", kind="curated")
    group = _import(database, task="ligand-md", run_id=run, kind="intermediate", origin={"source_node": node})
    context = {"run_id": run, "task_group": "ligand-md", "source_node": node,
        "study_id": study, "target_pdb_id": pdb, "target_receptor": "TRα" if pdb == "4LNW" else "TRβ",
        "source_pose_id": study + "-pose", "selected_complex_run_id": study + "-complex", "compound_id": "Ligand A",
        "replica_index": replica, "source_production_run_id": f"{study}-production-{replica}", "md_engine": "OpenMM",
        "source_summary_relative_path": f"artifacts/legacy/{study}-{replica}/result.json"}
    identifier = database.add_observation(group["group_id"], group["molecule_ids"][0], namespace, run,
        {"production_actual_duration_ns": 99}, context=context, source_key=f"{namespace}:{study}:{replica}:{node}")
    return identifier


def _inventory(database, user=USER):
    index = build_browser_index(database, user, collection="library", include_related=True)
    raw = subject_results(database, index, {"id": "all", "molecule_ids": set(index["molecules"]), "group_ids": set(index["groups"])})
    return index, saved_results(index, raw)


def test_md_study_uses_production_not_an_analysis_as_its_representative(database):
    production = _add_md(database)
    endpoint = _add_md(database, namespace="md.mmgbsa.pb", run="endpoint")
    index = build_browser_index(database, USER, collection="library", include_related=True)
    raw = subject_results(database, index, {"molecule_ids": set(index["molecules"]), "group_ids": set(index["groups"])})
    analysis = next(row for row in raw if row["observation_id"] == endpoint)
    analysis.update(id="000000-analysis-first", structure_ids=["endpoint-snapshot"])
    outcome = next(outcome for outcome in saved_results(index, raw) if outcome["outcome_kind"] == "Dynamics")
    assert outcome["primary"]["observation_id"] == production
    assert {row["observation_id"] for row in outcome["originals"]} == {production, endpoint}


def test_md_table_displays_engine_and_production_without_losing_endpoint_provenance(database):
    from mn_molecular_design_studio.data_browser_table import default_columns, evaluation_entries, result_table

    production = _add_md(database)
    endpoint = _add_md(database, namespace="md.mmgbsa.pb", run="endpoint")
    index, outcomes = _inventory(database)
    listing = listed_results([outcome for outcome in outcomes if outcome["outcome_kind"] == "Dynamics"])
    row = result_table(database, index, listing)[0]
    assert row["Method"] == "OpenMM" and row["Production"] == "Recorded"
    assert row["Replicas"] == 1 and row["Production segments"] == 1
    assert row["Saved analyses"] == "MM/PBSA (1)"
    assert {"Production", "Replicas", "Saved analyses"} <= set(default_columns([row]))
    exported = result_table(database, index, evaluation_entries(listing))
    assert {row["Observation ID"] for row in exported} == {production, endpoint}
    assert {row["Method"] for row in exported} == {"MD production", "MM/PBSA"}


def test_endpoint_only_study_does_not_invent_a_production_record(database):
    from mn_molecular_design_studio.data_browser_table import result_table

    _add_md(database, namespace="md.mmgbsa.pb", run="endpoint-only")
    index, outcomes = _inventory(database)
    listing = listed_results([outcome for outcome in outcomes if outcome["outcome_kind"] == "Dynamics"])
    row = result_table(database, index, listing)[0]
    assert row["Method"] == "OpenMM"
    assert row["Production"] == "Not recorded" and row["Production segments"] == 0
    assert row["Production status"] == "Not recorded"


def test_md_table_controls_render_the_study_engine_and_production_summary(database):
    _add_md(database)
    _add_md(database, namespace="md.mmgbsa.pb", run="endpoint")
    app = AppTest.from_string("from pathlib import Path\n"
        "import streamlit as st\n"
        "from mn_molecular_design_studio.data_browser import ReadOnlyLibrary, listed_results, saved_results\n"
        "from mn_molecular_design_studio.database_browser import build_browser_index, subject_results\n"
        "from mn_molecular_design_studio.data_browser_ui import _table_controls\n"
        f"database = ReadOnlyLibrary(Path({str(database.root)!r}))\n"
        f"index = build_browser_index(database, {USER!r}, collection='library', include_related=True)\n"
        "raw = subject_results(database, index, {'molecule_ids': set(index['molecules']), 'group_ids': set(index['groups'])})\n"
        "entries = listed_results([outcome for outcome in saved_results(index, raw) if outcome['outcome_kind'] == 'Dynamics'])\n"
        "rows, columns = _table_controls(database, index, entries, entries, 'test', 'md-study', 1)\n"
        "st.dataframe([{column: row.get(column) for column in columns} for row in rows])\n", default_timeout=30).run()
    assert not app.exception
    displayed = app.dataframe[0].value.iloc[0]
    assert displayed["Method"] == "OpenMM" and displayed["Production"] == "Recorded"
    assert displayed["Replicas"] == 1 and displayed["Production segments"] == 1
    assert displayed["Saved analyses"] == "MM/PBSA (1)"


def test_md_study_counts_replicas_and_segments_without_summing_unverified_continuations(database):
    from mn_molecular_design_studio.md_browser import study_summary

    originals = [{"id": f"{replica}-{segment}", "tool": "md.production", "run_id": f"{replica}-{segment}",
        "md": {"replica_index": replica, "source_production_run_id": f"{replica}-{segment}",
               "production_status": "completed", "md_engine": "OpenMM"},
        "metrics": {"production_actual_duration_ns": duration}}
        for replica in (1, 2, 3) for segment, duration in ((1, 1), (2, 99))]
    summary = study_summary(originals + originals)
    assert summary["Replicas"] == 3 and summary["Production segments"] == 6
    assert summary["Production status"] == "completed"
    assert not any("duration" in field.lower() for field in summary)


def _cached(reader, relative, payload):
    content = payload if isinstance(payload, bytes) else json.dumps(payload).encode()
    key = cache_key(reader.source, relative)
    reader.cache_root.mkdir(parents=True, exist_ok=True)
    (reader.cache_root / (key + ".data")).write_bytes(content)
    (reader.cache_root / (key + ".json")).write_text(json.dumps({"source": list(reader.source),
        "path": relative, "sha256": hashlib.sha256(content).hexdigest()}))
    return content


def _reader(database):
    _add_md(database)
    index, outcomes = _inventory(database)
    outcome = next(row for row in outcomes if row["outcome_kind"] == "Dynamics")
    observation = study_observations(database, USER, outcome)[0]
    run = index["runs"][outcome["primary"]["run_key"]]
    return MDArtifactReader(database, USER, observation, run), index, outcome


def test_md_separates_studies_targets_and_nodes_but_attaches_analysis_across_import_jobs(database):
    _add_md(database, replica=1)
    _add_md(database, replica=2)
    _add_md(database, namespace="md.analysis", run="analysis-import")
    _add_md(database, namespace="md.mmgbsa.gb", run="endpoint-import")
    _add_md(database, study="study-two", pdb="3GWS")
    _add_md(database, node="22222222-2222-4222-8222-222222222222")
    index, outcomes = _inventory(database)
    md = [row for row in outcomes if row["outcome_kind"] == "Dynamics"]
    assert len(md) == 3
    attached = next(row for row in md if len(row["originals"]) == 4)
    assert {row["md"]["study_id"] for row in attached["originals"]} == {"study-one"}
    assert {row["target_pdb_id"] for row in attached["originals"]} == {"4LNW"}
    nodes, selections = branch_nodes({"id": "folder"}, md)
    assert nodes[0]["label"] == "Molecular dynamics"
    assert len(nodes[0]["children"]) == 3
    replica_selections = [value for value in selections.values() if value.get("md_replica") == "2"]
    assert len(replica_selections) == 1
    assert [row["outcome"]["id"] for row in listed_results(md, selection=replica_selections[0])] == [attached["id"]]
    assert "4LNW" in nodes[0]["children"][0]["label"] or "3GWS" in nodes[0]["children"][0]["label"]


def test_md_lightweight_projection_retains_study_identity_without_loading_metrics(database):
    _add_md(database)
    index, outcomes = _inventory(database)
    primary = outcomes[0]["primary"]
    assert primary["metrics"] == {}
    assert primary["md"]["study_id"] == "study-one"
    assert primary["md"]["replica_index"] == 1
    assert primary["md"]["source_summary_relative_path"].endswith("result.json")


def test_legacy_md_with_recorded_pdb_but_no_canonical_target_is_visible_only_in_its_pdb_context(database):
    from test_target_registration import _design
    from mn_molecular_design_studio.molecular_database_ingest import ingest_run

    ingest_run(_design(database), database=database)
    _add_md(database)
    _add_md(database, study="study-two", pdb="3GWS")
    index, outcomes = _inventory(database, {**USER, "_cluster_project_ids": ["project-1"]})
    subject = folder_scopes(index, [])["targets:pdb:4LNW"]["subject"]
    found = saved_results(index, subject_results(database, index, subject))
    md = [row for row in found if row["outcome_kind"] == "Dynamics"]
    assert len(md) == 1 and md[0]["target_pdb_id"] == "4LNW"
    construct = {**subject, "pdb_md_run_keys": set()}
    assert not [row for row in saved_results(index, subject_results(database, index, construct)) if row["outcome_kind"] == "Dynamics"]


def test_report_choices_preserve_replica_and_deduplicate_energy_namespaces(database):
    _add_md(database, replica=1)
    _add_md(database, replica=2)
    _add_md(database, namespace="md.mmgbsa.gb")
    _add_md(database, namespace="md.mmgbsa.pb")
    index, outcomes = _inventory(database)
    observations = study_observations(database, USER, outcomes[0])
    assert len(report_choices(observations, replica="1")) == 1
    assert len(report_choices(observations, replica="2")) == 1
    assert len(report_choices(observations, replica="1", category="endpoint")) == 1
    assert not report_choices(observations, replica="3")
    different = next(row for row in observations if row["namespace"] == "md.mmgbsa.pb")
    different["context"] = {**different["context"], "source_production_run_id": "other-segment", "source_summary_relative_path": "artifacts/other.json"}
    exact = report_choices(observations, replica="1", category="endpoint", production_id="study-one-production-1")
    assert len(exact) == 1 and exact[0]["namespace"] == "md.mmgbsa.gb"


def test_artifact_reader_checks_cache_digest_sql_identity_fresh_acl_and_deletion(database):
    reader, index, outcome = _reader(database)
    relative = reader.observation["context"]["source_summary_relative_path"]
    content = _cached(reader, relative, {"success": True})
    assert reader.report() == {"success": True}
    with pytest.raises(ValueError, match="digest mismatch"):
        reader.text(relative, expected="a" * 64)
    outsider = MDArtifactReader(database, OTHER, reader.observation, reader.run)
    with pytest.raises(PermissionError):
        outsider.report()
    with database.connection() as connection:
        connection.execute("DELETE FROM molecule_observations WHERE id=?", (reader.observation["id"],))
    with pytest.raises(ValueError, match="changed"):
        reader.report()


@pytest.mark.parametrize("relative", ["../../secret.json", "/tmp/secrets.json", "artifacts/../metadata.json", "artifacts/md/trajectory.dcd", "artifacts\\secret.json"])
def test_artifact_reader_rejects_unsafe_paths_and_trajectories(database, relative):
    reader, index, outcome = _reader(database)
    with pytest.raises(ValueError):
        reader.text(relative)


def test_artifact_reader_rejects_cache_tampering_and_symlink_escape(database, tmp_path):
    reader, index, outcome = _reader(database)
    relative = reader.observation["context"]["source_summary_relative_path"]
    _cached(reader, relative, {"success": True})
    path = reader.cache_root / (cache_key(reader.source, relative) + ".data")
    path.write_text("changed")
    with pytest.raises(ValueError, match="digest mismatch"):
        reader.report()
    path.unlink()
    outside = tmp_path / "outside.json"
    outside.write_text("{}")
    path.symlink_to(outside)
    with pytest.raises(ValueError, match="escapes"):
        reader.report()


def test_reader_bounded_download_and_changed_evidence_are_not_presented(database):
    reader, index, outcome = _reader(database)
    reader.fetch_text = lambda *args: "x" * (MAX_ARTIFACT_BYTES + 1)
    with pytest.raises(ValueError, match="size limit"):
        reader.report()
    reader.fetch_text = lambda *args: None
    assert reader.report() is None
    with database.connection() as connection:
        connection.execute("UPDATE molecule_observations SET context=json_set(context,'$.study_id','changed') WHERE id=?", (reader.observation["id"],))
    with pytest.raises(ValueError, match="evidence changed"):
        reader.report()


def test_manifest_requires_exact_pose_and_target_and_snapshot_roles_are_not_guessed(database):
    reader, index, outcome = _reader(database)
    context = reader.observation["context"]
    study = {field: context[field] for field in ("study_id", "compound_id", "source_pose_id", "target_pdb_id", "selected_complex_run_id")}
    study["production_runs"] = [{"run_id": "production", "topology_path": "artifacts/md/topology.pdb", "topology_sha256": "b" * 64}]
    _cached(reader, "artifacts/md/study.json", {"schema_version": 1, "studies": [study]})
    assert matching_study(reader, outcome) == study
    references = snapshot_references(study, "production")
    assert [row["label"] for row in references] == ["Recorded MD topology snapshot"]
    assert "Not asserted" in references[0]["note"]
    assert not snapshot_references(study, "different-production")
    study["source_pose_id"] = "unrelated-pose"
    _cached(reader, "artifacts/md/study.json", {"schema_version": 1, "studies": [study]})
    with pytest.raises(ValueError, match="source_pose_id"):
        matching_study(reader, outcome)


def test_curves_use_their_own_time_axis_and_missing_ligand_series_does_not_hide_protein():
    streamlit = FakeStreamlit()
    _render_structural_dynamics(streamlit, {"time_ps": [5000, 6000],
        "radius_of_gyration": {"protein_angstrom": [10, 11]}},
        {"time_ps": [0, 1000], "backbone_rmsd_angstrom": [0, 1]}, label="one")
    charts = [args[0] for name, args, kwargs in streamlit.calls if name == "line_chart"]
    assert [row["time_ns"] for row in charts[0]] == [0, 1]
    assert [row["time_ns"] for row in charts[1]] == [5, 6]
    assert charts[0][1]["backbone_rmsd_angstrom"] == 1


def test_endpoint_sampling_uses_preserved_sql_parameters_without_inference():
    streamlit = FakeStreamlit()
    _render_endpoint_details(streamlit, {"mmgbsa": {"status": "success"}},
        {"sampling_parameters": {"stride": 10, "start_pct": 20, "end_pct": 100}})
    details = next(args[0] for name, args, kwargs in streamlit.calls if name == "json")
    assert details["Sampling stride"] == 10
    assert details["Sampling window start (%)"] == 20
    assert details["Sampling window end (%)"] == 100
    assert details["Frames analyzed"] == "not recorded"


@pytest.mark.parametrize("malformed", [False, True])
def test_middle_panel_shows_stored_md_graphs_and_summary_even_without_manifest(database, malformed):
    reader, index, outcome = _reader(database)
    _cached(reader, reader.observation["context"]["source_summary_relative_path"], {"success": True, "md_result": {
        "status": "completed", "analytics": {"rmsd": {"time_ps": [0, 1000], "backbone_rmsd_angstrom": [0, 1]},
            "structural_dynamics": {"time_ps": [0, 1000], "radius_of_gyration": {"protein_angstrom": [10, 11]}}}}})
    if malformed:
        _cached(reader, reader.observation["context"]["source_summary_relative_path"], {"success": True, "md_result": ["invalid"]})
    source = f'''from mn_molecular_design_studio.data_browser import ReadOnlyLibrary, saved_results, listed_results
from mn_molecular_design_studio.database_browser import build_browser_index, subject_results
from mn_molecular_design_studio.md_browser_ui import render_md_study
database=ReadOnlyLibrary({str(database.root)!r})
user={USER!r}
index=build_browser_index(database,user,collection="library",include_related=True)
raw=subject_results(database,index,{{"id":"all","molecule_ids":set(index["molecules"]),"group_ids":set(index["groups"])}})
outcome=next(row for row in saved_results(index,raw) if row["outcome_kind"]=="Dynamics")
render_md_study(database,user,index,listed_results([outcome])[0],{{}})
'''
    app = AppTest.from_string(source).run()
    assert not app.exception
    if malformed:
        assert any("must be a JSON object" in row.value for row in app.warning)
    else:
        assert len(app.get("arrow_vega_lite_chart")) + len(app.get("vega_lite_chart")) == 2, str(app)
        assert any("recorded time axis" in row.value for row in app.caption)
    next(element for element in app.selectbox if element.label == "Analysis").set_value("SQL summary").run()
    assert not app.exception
    assert len(app.dataframe) == 1
