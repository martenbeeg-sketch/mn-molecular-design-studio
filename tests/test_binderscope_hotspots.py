from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
from types import SimpleNamespace

import gemmi
import pytest

from test_database_browser import USER, OTHER, database, _import
from test_binderscope import saved_subject, browser_entry, write_complex
from mn_molecular_design_studio.binderscope import prepare_browser_selection
from mn_molecular_design_studio.binderscope_sources import target_context
from mn_molecular_design_studio.complex_analysis_sources import build_saved_structure_collection
from mn_molecular_design_studio.data_browser import load_saved_record
from mn_molecular_design_studio.engine.workflow_engine.tools.binderscope import normalize_complex, read_structure
from mn_molecular_design_studio.engine.workflow_engine.tools.binderscope_context import prepare_context, hotspot_result


def original_target(tmp_path, *, insertion=False):
    path = tmp_path / "original.pdb"
    write_complex(path, insertion=insertion)
    structure = gemmi.read_structure(str(path))
    del structure[0]["Y"]
    structure.write_pdb(str(path))
    return path


def hotspot_subject(database, tmp_path):
    imported, saved, selector = saved_subject(database, tmp_path)
    reference = original_target(tmp_path)
    digest = hashlib.sha256(reference.read_bytes()).hexdigest()
    target = _import(database, run_id="target", records=[{"name": "Original target", "components": [
        {"type": "protein", "chain": "X", "sequence": "A" * 16}], "files": [{
            "name": "target.pdb", "content": reference.read_bytes(), "input_role": "target",
            "structure_role": "source_reference", "source_structure_digest": digest}]}])
    with database.connection() as connection:
        connection.execute("UPDATE molecule_members SET metadata=? WHERE group_id=? AND molecule_id=?", (json.dumps({
            "design_target_molecule_id": target["molecule_ids"][0], "target_type": "structured",
            "candidate_set": {"source_provenance": {"input": {"params": {"hotspots": "X10,X13", "target_chains": ["X"]}}}}}),
            imported["group_id"], saved["id"]))
    entry = browser_entry(imported, saved, selector)
    return imported, saved, selector, entry, target


def test_sql_hotspots_and_original_reference_are_source_scoped_and_revalidated(database, tmp_path):
    imported, saved, selector, entry, target = hotspot_subject(database, tmp_path)
    database.add_observation(imported["group_id"], saved["id"], "unrelated", "other-model", {}, context={
        "representation_id": "not-this-structure", "settings": {"hotspots": "X99"}})
    selection, excluded = prepare_browser_selection(database, USER, [entry])
    assert not excluded
    context = selection["sources"][0]["target_context"]
    assert context["hotspots"]["requested"] == ["X10", "X13"]
    assert context["reference"]["status"] == "available"
    assert context["reference"]["molecule_id"] == target["molecule_ids"][0]
    adapter = build_saved_structure_collection(USER, selection, workflow_id="analysis.binderscope", project_id=None, database=database)
    from mn_molecular_design_studio.engine.workflow_engine.workflows.binderscope import _saved_inputs

    candidate = _saved_inputs(Path(adapter["source_run_dir"]), Path(adapter["candidates_jsonl"]), use_recorded_settings=True)[0]
    reference = candidate["raw_metadata"]["complex_analysis_selection"]["target_reference"]
    assert Path(reference["path"]).is_relative_to(Path(adapter["source_run_dir"]))
    Path(reference["path"]).write_bytes(b"changed")
    with pytest.raises(ValueError, match="original target reference changed"):
        _saved_inputs(Path(adapter["source_run_dir"]), Path(adapter["candidates_jsonl"]), use_recorded_settings=True)
    database.add_observation(imported["group_id"], saved["id"], "generator", "original", {}, context={
        "representation_id": selector["structure_id"], "settings": {"hotspots": "X14"}})
    with pytest.raises(ValueError, match="evidence changed"):
        build_saved_structure_collection(USER, selection, workflow_id="analysis.binderscope", project_id=None, database=database)
    changed, excluded = prepare_browser_selection(database, USER, [entry])
    assert not excluded and changed["sources"][0]["target_context"]["hotspots"]["status"] == "conflicting"


def test_inaccessible_original_reference_does_not_authorize_hotspot_mapping(database, tmp_path):
    imported, saved, selector, entry, target = hotspot_subject(database, tmp_path)
    with database.connection() as connection:
        connection.execute("UPDATE molecule_groups SET owner_id=? WHERE id=?", (OTHER["id"], target["group_id"]))
    selection, excluded = prepare_browser_selection(database, USER, [entry])
    assert not excluded
    assert selection["sources"][0]["target_context"]["reference"]["status"] == "unavailable"


def test_reprediction_renumbering_and_insertion_codes_follow_sequence_identity(tmp_path, monkeypatch):
    reference = original_target(tmp_path, insertion=True)
    source = tmp_path / "repredicted.pdb"
    write_complex(source, target_chain="Q", target_offset=501)
    canonical = tmp_path / "canonical.pdb"
    mapping, sequences = normalize_complex(source, canonical, "Y", "Q")
    context = {"hotspots": {"status": "recorded", "requested": ["X10", "X13", "X15A", "X999", "Y101"],
                           "numbering": "input_target_author"}, "target_kind": "structured", "reference": {"status": "available"}}
    prepared = prepare_context(context, {"path": str(reference), "chain": "X"}, mapping, sequences, canonical, read_structure)
    assert [(row["requested"], row.get("saved_residue")) for row in prepared["rows"][:3]] == [
        ("X10", "Q:501"), ("X13", "Q:504"), ("X15A", "Q:506")]
    assert prepared["target_reference_comparison"]["target_ca_rmsd_angstrom"] < 1e-6
    assert all(row["status"] == "unmapped" for row in prepared["rows"][3:])
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1] / "mn_molecular_design_studio/vendor/binderscope/binderscope/src"))
    result = hotspot_result(prepared, SimpleNamespace(hotspot_dist={1: 3.0, 4: 8.0, 6: 3.5}, hot=[1, 6]))
    assert (result["contacted_count"], result["missed_count"], result["unmapped_count"]) == (2, 1, 2)
    assert result["status"] == "partial" and result["coverage_fraction"] is None


@pytest.mark.parametrize("condition", ["repeated_crop", "sequence_change", "sequence_only", "no_reference"])
def test_uncertain_target_context_is_not_a_pass(tmp_path, condition):
    reference = original_target(tmp_path)
    source = tmp_path / "design.pdb"
    write_complex(source)
    if condition == "repeated_crop":
        structure = gemmi.read_structure(str(source))
        del structure[0]["X"][0]
        structure.write_pdb(str(source))
    elif condition == "sequence_change":
        structure = gemmi.read_structure(str(source))
        structure[0]["X"][0].name = "GLY"
        structure.write_pdb(str(source))
    canonical = tmp_path / "canonical.pdb"
    mapping, sequences = normalize_complex(source, canonical, "Y", "X")
    context = {"hotspots": {"status": "recorded", "requested": ["X10"], "numbering": "input_target_author"},
               "reference": {"status": "available"}, "target_kind": "sequence" if condition == "sequence_only" else "structured"}
    prepared = prepare_context(context, None if condition == "no_reference" else {"path": str(reference), "chain": "X"},
                               mapping, sequences, canonical, read_structure)
    if condition == "sequence_only":
        assert prepared["target_reference_comparison"]["status"] == "not_applicable"
        assert prepared["rows"][0]["status"] == "mapped"
    else:
        assert prepared["rows"][0]["status"] == "unmapped"
        assert prepared["target_reference_comparison"]["status"] == "unavailable"


def test_unique_crop_maps_original_residues_without_assuming_constant_numbers(tmp_path):
    reference = original_target(tmp_path, insertion=True)
    structure = gemmi.read_structure(str(reference))
    structure[0]["X"][5].name = "SER"
    structure.write_pdb(str(reference))
    source = tmp_path / "repredicted.pdb"
    write_complex(source, insertion=True)
    structure = gemmi.read_structure(str(source))
    structure[0]["X"][5].name = "SER"
    for index in range(15, 12, -1):
        del structure[0]["X"][index]
    del structure[0]["X"][0]
    del structure[0]["X"][0]
    for index, residue in enumerate(structure[0]["X"]):
        residue.seqid = gemmi.SeqId(1001 + index, " ")
    structure.write_pdb(str(source))
    canonical = tmp_path / "canonical.pdb"
    mapping, sequences = normalize_complex(source, canonical, "Y", "X")
    context = {"hotspots": {"status": "recorded", "requested": ["X10", "X15A"], "numbering": "input_target_author"},
               "reference": {"status": "available"}, "target_kind": "structured"}
    prepared = prepare_context(context, {"path": str(reference), "chain": "X"}, mapping, sequences, canonical, read_structure)
    assert prepared["rows"][0]["status"] == "unmapped"
    assert prepared["rows"][1]["saved_residue"] == "X:1004"
    assert prepared["target_reference_comparison"]["mapped_residues"] == 11
    assert prepared["target_reference_comparison"]["scope"] == "unique contiguous target segment"


def test_result_viewer_shows_original_and_repredicted_residue_labels():
    from streamlit.testing.v1 import AppTest

    result = {"analysis": {"dsasa": 500.0}, "hotspot_compliance": {"status": "assessed",
        "requested_count": 1, "contacted_count": 1, "missed_count": 0, "unmapped_count": 0, "cutoff_angstrom": 4.5,
        "residues": [{"requested": "X15A", "saved_residue": "Q:506", "status": "contacted", "minimum_distance_angstrom": 3.5}]},
        "target_reference_comparison": {"status": "not_applicable", "reason": "Sequence-only target: no fixed fold assumed."}}
    app = AppTest.from_string(f'''
from mn_molecular_design_studio.binderscope_ui import render_design_result
render_design_result({result!r})
''').run()
    assert not app.exception
    assert app.dataframe[0].value.iloc[0]["Requested residue"] == "X15A"
    assert app.dataframe[0].value.iloc[0]["Saved target residue"] == "Q:506"
    assert any("no fixed fold" in caption.value for caption in app.caption)


@pytest.mark.skipif(os.environ.get("MDS_TEST_BINDERSCOPE_IMAGE") != "1", reason="Explicit offline hotspot SQL integration")
@pytest.mark.parametrize("replicated", [False, True])
def test_offline_engine_assesses_original_hotspots_and_attaches_full_result(database, tmp_path, monkeypatch, replicated):
    from mn_molecular_design_studio.engine.workflow_engine.workflows.binderscope import run_binderscope
    from mn_molecular_design_studio.molecular_database_ingest import ingest_run
    from mn_molecular_design_studio.binderscope import saved_result

    imported, saved, selector, entry, target = hotspot_subject(database, tmp_path)
    if replicated:
        from mn_molecular_design_studio.cluster import config, transport
        from mn_molecular_design_studio.cluster.replication import Replica
        from mn_molecular_design_studio.molecular_database import MolecularDatabase, _id
        from test_molecular_database_peer import SOURCE, _sync_all

        monkeypatch.setenv("MDS_CLUSTER_HOME", str(tmp_path / "cluster"))
        writer = database
        for molecule in writer.list_molecules(USER, collection="all")["items"]:
            for representation in molecule["representations"]:
                writer.make_structure_portable(USER, representation["id"])
        database = MolecularDatabase(tmp_path / "receiver")
        _sync_all(Replica(), writer, database)
        monkeypatch.setattr(config, "load", lambda: {"mode": "disabled", "name": "evaluation-node",
            "peers": [{"id": SOURCE, "url": "https://owner.invalid"}]})
        monkeypatch.setattr(transport, "node_request", lambda path, *args, **kwargs:
            writer.structure_file(USER, path.rsplit("/", 1)[-1]).read_bytes())
        imported = {**imported, "group_id": _id("peer-", SOURCE + ":" + imported["group_id"])}
        selector = {**selector, "group_id": imported["group_id"],
                    "structure_id": _id("peer-structure-", SOURCE + ":" + selector["structure_id"])}
        entry = browser_entry(imported, saved, selector)
    selection, excluded = prepare_browser_selection(database, USER, [entry])
    assert not excluded
    adapter = build_saved_structure_collection(USER, selection, workflow_id="analysis.binderscope", project_id=None, database=database)
    run = run_binderscope(Path(adapter["source_run_dir"]), Path(adapter["candidates_jsonl"]), cpu_cores=1, use_recorded_settings=True)
    assert json.loads((run / "result.json").read_text())["success"], (run / "stderr.log").read_text()
    (run / "molecular_design_studio.json").write_text(json.dumps({"application": "mn-molecular-design-studio", "studio_owner_id": USER["id"]}))
    assert ingest_run(run, database=database, force=True) == 1
    with database.connection() as connection:
        row = dict(connection.execute("SELECT id,group_id FROM molecule_observations WHERE namespace='binderscope'").fetchone())
    result = saved_result(database, USER, {"molecule_id": saved["id"], "group_id": row["group_id"], "observation_id": row["id"]})
    hotspot = result["hotspot_compliance"]
    assert hotspot["status"] == "assessed" and hotspot["coverage_fraction"] == .5
    assert (hotspot["contacted_count"], hotspot["missed_count"], hotspot["unmapped_count"]) == (1, 1, 0)
    assert result["target_reference_comparison"]["target_ca_rmsd_angstrom"] < 1e-6
    assert result["analysis"]["hotspots"] == [1, 4]
    assert result["binding_modes"]["hotspots_touched"]["hot_of"] == 2
    from mn_molecular_design_studio.molecular_database_ingest import _database_analysis_subject

    source = json.loads(Path(adapter["candidates_jsonl"]).read_text().splitlines()[0])["raw_metadata"]["database_source"]
    assert _database_analysis_subject(database, USER, USER["id"], "", source)
    with database.connection() as connection:
        connection.execute("UPDATE molecule_members SET metadata='{}' WHERE group_id=? AND molecule_id=?", (imported["group_id"], saved["id"]))
    assert _database_analysis_subject(database, USER, USER["id"], "", source) is None
