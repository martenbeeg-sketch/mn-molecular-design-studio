from __future__ import annotations

import hashlib
from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest

from mn_molecular_design_studio.complex_analysis_collection import (
    aggregate_contact_footprints,
    assignment_rows,
    collection_page,
    evaluation_table_rows,
    filter_evaluation_rows,
    heatmap_selection,
    prepare_inspector_navigation_state,
    sequence_correspondence,
)
from mn_molecular_design_studio.molecular_database import MolecularDatabase


USER = {"id": "a" * 32, "role": "user", "_cluster_project_ids": []}
OTHER_USER = {"id": "b" * 32, "role": "user", "_cluster_project_ids": []}


def _atom(serial: int, chain: str, number: int, x: float, *, residue: str = "ALA", insertion: str = " ") -> str:
    return (
        f"ATOM  {serial:5d}  CA  {residue:>3s} {chain}{number:4d}{insertion:1s}   "
        f"{x:8.3f}{1.0:8.3f}{1.0:8.3f}{1.0:6.2f}{20.0:6.2f}           C"
    )


def _ppi(sequence: tuple[str, ...] = ("ALA", "GLY")) -> bytes:
    lines = []
    serial = 1
    for chain, residues in zip(("A", "B"), (sequence, ("SER", "THR"))):
        for number, residue in enumerate(residues, start=1):
            lines.append(_atom(serial, chain, number, float(serial), residue=residue))
            serial += 1
    return ("\n".join([*lines, "END", ""])).encode()


def _isolated_database(tmp_path, monkeypatch) -> MolecularDatabase:
    workspace = tmp_path / "workspace"
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_HOME", str(workspace))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(workspace))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR", str(workspace / "runs"))
    return MolecularDatabase(workspace / "database")


def _track(chain: str = "A") -> dict:
    return {
        "chain_id": chain,
        "sequence": "AG",
        "residues": [
            {"sequence_index": 0, "amino_acid": "A", "author_number": 10, "insertion_code": "", "selector": f"{chain}:10^_"},
            {"sequence_index": 1, "amino_acid": "G", "author_number": 10, "insertion_code": "A", "selector": f"{chain}:10^A"},
        ],
    }


def test_sequence_correspondence_accepts_unique_truncated_construct_and_exact_sequence():
    assert sequence_correspondence("MABCDEK", "ABCDE") == {
        "status": "unique_exact_subsequence", "compatible": True,
        "mapping": [1, 2, 3, 4, 5],
        "reason": "Unique exact subsequence at reference index 2; partial coverage is expected.",
    }
    assert sequence_correspondence("AG", "AG")["status"] == "exact"


def test_sequence_correspondence_separates_repeated_and_nonidentical_sequences():
    assert sequence_correspondence("ABCABC", "ABC")["status"] == "ambiguous"
    changed = sequence_correspondence("ABCDE", "ABXDE")
    assert changed["status"] == "incompatible"
    assert not changed["compatible"]


def test_contact_mapping_preserves_author_insertion_codes_and_partial_offset():
    rows = [{"binder_residue": "A:10A", "target_residue": "B:22"}]
    correspondence = sequence_correspondence("MAGE", "AG")
    assert correspondence["mapping"] == [1, 2]
    positions = aggregate_contact_footprints([{
        "identity": "structure-1", "tracked_chain": "A", "tracked_track": _track(),
        "candidate_to_reference": correspondence["mapping"], "contact_rows": rows,
    }], {"residues": [{"sequence_index": i, "amino_acid": "A", "author_number": i + 1,
                        "insertion_code": "", "selector": f"A:{i + 1}^_"} for i in range(4)]})
    assert [row["reference_index"] for row in positions if row["contact_count"]] == [2]


def test_blank_author_number_is_ambiguous_when_inserted_residue_shares_number():
    from mn_molecular_design_studio.complex_analysis_collection import author_residue_index

    track = _track()
    assert author_residue_index(track, "A:10^A", "A") == 1
    assert author_residue_index(track, "A:10", "A") is None
    assert author_residue_index(track, "B:10^A", "A") is None


def test_saved_zero_contact_map_is_distinct_from_missing_or_unavailable_evaluation():
    reference = {"residues": [{"sequence_index": 0, "amino_acid": "A", "author_number": 10,
                               "insertion_code": "", "selector": "A:10^_"}]}
    result = aggregate_contact_footprints([
        {"identity": "zero", "contact_rows": [], "contact_status": "saved"},
        {"identity": "missing", "contact_rows": None, "contact_status": "missing evaluation"},
        {"identity": "unavailable", "contact_rows": None, "contact_status": "artifact unavailable"},
        {"identity": "failed", "contact_rows": None, "contact_status": "failed evaluation"},
    ], reference)
    residue = result[0]
    assert residue["contact_count"] == 0
    assert residue["evaluable_count"] == 1
    assert residue["frequency"] == 0
    assert residue["member_status"] == {
        "zero": "no contact", "missing": "missing evaluation", "unavailable": "artifact unavailable",
        "failed": "failed evaluation",
    }
    assert residue["failed_count"] == 1


def test_contact_heatmap_styling_works_with_current_pandas_styler_api():
    import pandas as pd

    from mn_molecular_design_studio.complex_analysis_collection_ui import _style_contact_heatmap

    rendered = _style_contact_heatmap(pd.DataFrame({"A10": ["●", "?", "·", "×"]}), ["A10"]).to_html()
    assert "#7252a3" in rendered
    assert "#fff3cd" in rendered
    assert "#f5c2c7" in rendered


def test_collection_pagination_scales_beyond_comparison_slots():
    members = [{"identity": f"structure-{index}"} for index in range(137)]
    assert len(collection_page(members, page=0, page_size=25)) == 25
    assert len(collection_page(members, page=5, page_size=25)) == 12
    assert collection_page(members, page=5, page_size=25)[-1]["identity"] == "structure-136"


def test_heatmap_selection_and_inspector_navigation_keep_database_identity_and_residue():
    assert heatmap_selection([(1, "G10A")], ["s0", "s1"], {"G10A": 7}) == {
        "member_id": "s1", "reference_index": 7, "column": "G10A",
    }
    state = {}
    prepare_inspector_navigation_state(
        state, molecule_id="molecule-1", structure_id="s1",
        evaluation_id="evaluation-9", selector="A:10^A",
    )
    assert state["complex-analysis-db-viewer-s1-inspector-highlights"] == ["A:10^A"]
    assert state["complex-analysis-structure-id-molecule-1"] == "s1"
    assert state["complex-analysis-evaluation-molecule-1"] == "evaluation-9"


def test_evaluation_projection_retains_method_metrics_and_source_provenance():
    source_context = {
        "source_structure_digest": "sha256-original",
        "candidate": {"candidate_id": "pose-2", "raw_metadata": {"engine_version": "engine/4"}},
        "settings": {"contact_distance": 5.0, "include_hydrogens": False},
    }
    rows = evaluation_table_rows([{
        "molecule_id": "molecule-1", "structure_id": "structure-1", "structure_label": "Complex A",
        "observation": {"id": "observation-1", "namespace": "biotite_interface", "run_id": "run-1",
                        "metrics": {"contact_count": 12}, "created_at": "2026-10-01"},
        "evaluation": {"id": "evaluation-1", "evaluator": "biotite_interface", "engine_version": "Biotite 1.0",
                       "candidate_id": "candidate-2", "representation_id": "structure-1"},
        "context": source_context,
    }, {
        "molecule_id": "molecule-1", "structure_id": "structure-1", "structure_label": "Complex A",
        "observation": {"id": "observation-2", "namespace": "alphafold", "metrics": {"pLDDT": 88}},
        "evaluation": {"id": "evaluation-2", "evaluator": "alphafold"}, "context": {},
    }])
    assert rows[0]["category"] == "Scientific evaluation"
    assert rows[0]["structure_digest"] == "sha256-original"
    assert rows[0]["parameters"] == source_context["settings"]
    assert rows[0]["candidate_id"] == "candidate-2"
    assert rows[0]["provenance"]["candidate"]["candidate_id"] == "pose-2"
    assert rows[0]["source_engine"] == "Biotite 1.0"
    assert rows[1]["category"] == "Prediction confidence"
    assert filter_evaluation_rows(rows, method="biotite_interface", metric="contact_count", minimum=10) == [rows[0]]


def test_results_cohort_choices_keep_same_native_candidate_separate_by_job_and_model(monkeypatch):
    from mn_molecular_design_studio import complex_analysis_ui, complex_analysis_results_ui, run_results
    from mn_molecular_design_studio.complex_analysis_collection_ui import _results_structure_choices

    local = {"task_group": "design", "run_id": "rfd4-run", "tool": "rfd4", "run_dir": "/runs/design/rfd4-run",
             "project_id": "project-1"}
    remote = {"task_group": "design", "run_id": "bindcraft-run", "tool": "bindcraft2",
              "node_id": "peer-4033", "cluster_job_id": "peerjob-4033", "local_reference": "design/bindcraft-run",
              "project_id": "project-1", "run_dir": ""}
    candidate = {"candidate_id": "design_0", "binder_sequence": "AAAA", "target_chains": ["A"],
                 "binder_chains": ["B"], "raw_metadata": {}}
    monkeypatch.setattr(complex_analysis_ui, "_job_structure_options", lambda _user, _job: [
        {"candidate": {**candidate, "model_index": 0}, "candidate_id": "design_0",
         "relative_path": "models/model0.cif", "path": Path("/tmp/local0.cif"), "label": "model0"},
    ])
    monkeypatch.setattr(run_results, "read_source", lambda _job, _user, **_kwargs: {"candidates": [candidate]})
    monkeypatch.setattr(complex_analysis_results_ui, "_remote_board_options", lambda _job, _collection: [
        {"candidate": {**candidate, "model_index": 0}, "candidate_id": "design_0",
         "relative_path": "artifacts/models/model0.cif", "structure_field": "complex_pdb", "label": "model0"},
        {"candidate": {**candidate, "model_index": 1}, "candidate_id": "design_0",
         "relative_path": "artifacts/models/model1.cif", "structure_field": "complex_pdb", "label": "model1"},
    ])

    rows = _results_structure_choices(USER, [local, remote])

    assert len(rows) == 3
    assert len({row["identity"] for row in rows}) == 3
    assert {row["reference"] for row in rows} == {"design/rfd4-run", "cluster:peerjob-4033"}
    assert {row["model_id"] for row in rows} == {"0", "1"}
    assert all(row["candidate_id"] == "design_0" for row in rows)
    assert sum(row["remote"] for row in rows) == 2


def test_unified_candidate_row_keeps_generator_prediction_and_acceptance_distinct():
    from mn_molecular_design_studio.complex_analysis_collection_ui import _cohort_inventory_row

    sample = {"identity": "results:source-model", "source_kind": "Results", "kind": "protein-protein",
        "candidate_id": "studio-candidate", "model_id": "model-2", "source_job": "design/run-1",
        "source_workstation": "PC004033", "digest": "sha256:pose", "generator_engine": "rfd4",
        "prediction_engine": "boltz-2", "evaluation_entries": [{"method": "biotite_interface"}],
        "tracks": [{"chain_id": "A", "sequence": "AG"}],
        "candidate": {"candidate_id": "studio-candidate", "source_tool": "boltz-2", "stage": "prediction",
            "target_chains": ["A"], "binder_chains": ["B"], "binder_sequence": "AAAA", "metrics": {"plddt": 88.0}, "raw_metadata": {
                "source_candidate": {"candidate_id": "native-design-7", "source_tool": "rfd4",
                                     "metrics": {"native_pass_filters": True}},
                "validation_acceptance": {"status": "unresolved"},
                "human_review": {"selected": True},
            }}}
    row = _cohort_inventory_row(sample)

    assert row["Candidate / molecule"] == "studio-candidate"
    assert row["Original candidate"] == "native-design-7"
    assert row["Generator"] == "rfd4" and row["Prediction engine"] == "boltz-2"
    assert row["Native acceptance"] == "passed"
    assert row["Studio filter"] == "unresolved" and row["Manual selection"] is True
    assert row["Source job"] == "design/run-1" and row["Workstation"] == "PC004033"
    assert row["Structure SHA-256"] == "sha256:pose"
    assert row["Target"] == f"Target sequence · 2 aa · 1 chain(s) · {hashlib.sha256(b'AG').hexdigest()[:12]}"


def test_results_inspector_opener_uses_existing_authorized_local_route_and_residue_highlight(monkeypatch):
    from mn_molecular_design_studio import complex_analysis_ui
    from mn_molecular_design_studio.complex_analysis_collection_ui import _open_in_inspector

    captured = {}
    monkeypatch.setattr(complex_analysis_ui, "open_result_complex",
                        lambda run_dir, path, candidate, *, user: captured.update(
                            run_dir=run_dir, path=path, candidate=candidate, user=user))
    sample = {"source_kind": "Results", "remote": False, "user": USER,
        "job": {"task_group": "design", "run_id": "run-1", "run_dir": "/runs/design/run-1"},
        "path": Path("/runs/design/run-1/model.cif"), "relative_path": "model.cif",
        "candidate": {"candidate_id": "candidate-4"}}

    _open_in_inspector(sample, "A:10^_", "evaluation-2")

    assert captured["candidate"]["candidate_id"] == "candidate-4"
    assert captured["run_dir"] == "/runs/design/run-1"
    st_key = "complex-analysis-results-viewer-design/run-1-model.cif-inspector-highlights"
    assert __import__("streamlit").session_state[st_key] == ["A:10^_"]


def test_cohort_contact_snapshot_keeps_saved_empty_contacts_separate_from_missing():
    from mn_molecular_design_studio.complex_analysis_collection_ui import _contact_snapshot

    assert _contact_snapshot({"contact_snapshot": {"status": "saved", "rows": [], "contact_distance": 5.0}})["rows"] == []
    assert _contact_snapshot({"contact_snapshot": {"status": "missing evaluation", "rows": None}})["rows"] is None


def test_saved_interface_clustering_and_foldseek_remain_distinct():
    rows = assignment_rows([
        {"method": "analysis.interface_clustering", "structure_label": "Complex A", "evaluation_id": "i1",
         "molecule_id": "m1", "structure_id": "s1", "candidate_id": "candidate-1",
         "metrics": {"diversity_cluster": 2, "diversity_representative": False,
                     "diversity_representative_id": "candidate-2"}},
        {"method": "analysis.interface_clustering", "structure_label": "Complex B", "evaluation_id": "i2",
         "molecule_id": "m2", "structure_id": "s2", "candidate_id": "candidate-2",
         "metrics": {"diversity_cluster": 2, "diversity_representative": True,
                     "diversity_representative_id": "candidate-2"}},
        {"method": "analysis.foldseek", "structure_label": "Complex C", "evaluation_id": "f1",
         "metrics": {"cluster_id": "fold-3", "representative": False}},
    ])
    assert rows[0]["Method"] == "Interface Clustering · binding-mode assignment"
    assert rows[2]["Method"] == "Foldseek · structural-similarity assignment"
    assert rows[0]["Assignment"] == 2 and rows[0]["Representative"] is False
    assert rows[0]["Representative structure ID"] == "s2"


def test_cluster_members_resolve_unique_saved_representative_within_same_run():
    rows = assignment_rows([
        {"method": "analysis.interface_clustering", "run_id": "run-a", "structure_label": "A",
         "molecule_id": "m", "structure_id": "s-a", "candidate_id": "c-a",
         "metrics": {"diversity_cluster": 1, "diversity_representative": False}},
        {"method": "analysis.interface_clustering", "run_id": "run-a", "structure_label": "B",
         "molecule_id": "m", "structure_id": "s-b", "candidate_id": "c-b",
         "metrics": {"diversity_cluster": 1, "diversity_representative": True}},
        {"method": "analysis.interface_clustering", "run_id": "run-b", "structure_label": "C",
         "molecule_id": "m", "structure_id": "s-c", "candidate_id": "c-c",
         "metrics": {"diversity_cluster": 1, "diversity_representative": False}},
        {"method": "analysis.interface_clustering", "run_id": "run-b", "structure_label": "D",
         "molecule_id": "m", "structure_id": "s-d", "candidate_id": "c-d",
         "metrics": {"diversity_cluster": 1, "diversity_representative": True}},
    ])
    assert rows[0]["Representative structure ID"] == "s-b"
    assert rows[0]["Representative resolution"] == "unique saved representative in this run and cluster"
    assert rows[2]["Representative structure ID"] == "s-d"
    assert rows[0]["Representative structure ID"] != rows[2]["Representative structure ID"]


def test_database_collection_listing_and_structure_access_respect_group_acl(tmp_path, monkeypatch):
    from mn_molecular_design_studio import molecule_imports
    from mn_molecular_design_studio.complex_analysis_collection_ui import _accessible_records

    database = _isolated_database(tmp_path, monkeypatch)
    imported = molecule_imports.parse_structure(_ppi(), "saved-ppi.pdb")[0]
    result = database.import_records(USER, [imported], name="ACL scoped PPI")
    record = database.list_molecules(USER, molecule_ids=result["molecule_ids"], collection="all")["items"][0]
    structure_id = record["representations"][0]["id"]

    assert [row["id"] for row in _accessible_records(database, USER, "")] == [record["id"]]
    assert _accessible_records(database, OTHER_USER, "") == []
    with pytest.raises(PermissionError):
        database.structure_file(OTHER_USER, structure_id)
    assert hashlib.sha256(database.structure_file(USER, structure_id).read_bytes()).hexdigest()


def test_structure_discovery_is_lightweight_paginated_and_acl_scoped(tmp_path, monkeypatch):
    from mn_molecular_design_studio import evaluations, molecule_imports

    database = _isolated_database(tmp_path, monkeypatch)
    ids = []
    for index, sequence in enumerate((("ALA", "GLY"), ("GLY", "SER"), ("SER", "ALA"))):
        imported = molecule_imports.parse_structure(_ppi(sequence), "saved-ppi.pdb")[0]
        result = database.import_records(USER, [imported], name=f"PPI {index}")
        ids.extend(result["molecule_ids"])
    record = database.list_molecules(USER, molecule_ids=[ids[0]], collection="all")["items"][0]
    database.add_observation(record["groups"][0]["id"], ids[0], "biotite_interface", "saved-run",
                             {"contact_count": 7}, context={"representation_id": record["representations"][0]["id"]})
    monkeypatch.setattr(evaluations, "rebuild_missing", lambda *_args: (_ for _ in ()).throw(
        AssertionError("structure discovery must not rebuild evaluation projections")))

    first = database.list_structure_records(USER, collection="all", limit=2)
    second = database.list_structure_records(USER, collection="all", limit=2,
                                             offset=first["next_offset"])
    denied = database.list_structure_records(OTHER_USER, collection="all", limit=10)

    assert len(first["items"]) == 2
    assert first["next_offset"] == 2
    assert len(second["items"]) == 1
    assert second["next_offset"] is None
    assert {row["id"] for row in first["items"] + second["items"]} == set(ids)
    assert "observations" not in first["items"][0]
    assert "sources" not in first["items"][0]
    assert denied["items"] == []


def test_structure_discovery_respects_search_and_rejects_invalid_pages(tmp_path, monkeypatch):
    from mn_molecular_design_studio import molecule_imports

    database = _isolated_database(tmp_path, monkeypatch)
    imported = molecule_imports.parse_structure(_ppi(), "saved-ppi.pdb")[0]
    imported["name"] = "Unique Discovery Name"
    database.import_records(USER, [imported], name="Unique Discovery Name")

    result = database.list_structure_records(USER, collection="all", query="Unique Discovery", limit=10)

    assert len(result["items"]) == 1
    assert result["items"][0]["name"] == "Unique Discovery Name"
    with pytest.raises(ValueError, match="page size or offset"):
        database.list_structure_records(USER, limit=0)


def test_streamlit_collection_picker_exposes_more_than_five_accessible_records(monkeypatch):
    from mn_molecular_design_studio import molecular_database

    records = [{"id": f"record-{i}", "name": f"Complex {i}",
                "representations": [{"id": f"structure-{i}", "label": "model"}]}
               for i in range(8)]

    class FakeDatabase:
        def list_structure_records(self, user, **_kwargs):
            assert user == USER
            return {"items": records, "next_offset": None}

        def list_molecules(self, user, **_kwargs):
            assert user == USER
            return {"items": records}

    monkeypatch.setattr(molecular_database, "MolecularDatabase", FakeDatabase)
    app = AppTest.from_string(
        "from mn_molecular_design_studio.complex_analysis_collection_ui import render_collection_workspace\n"
        f"render_collection_workspace({USER!r})\n",
        default_timeout=10,
    ).run()
    assert not app.exception
    record_key = "complex-analysis-collection-records-" + hashlib.sha256(b"").hexdigest()[:10]
    assert len(app.multiselect(key=record_key).options) == 8
    app.multiselect(key=record_key).set_value([f"record-{i}" for i in range(8)]).run()
    assert not app.exception
    assert len(app.multiselect(key="complex-analysis-collection-structures").options) == 8


def test_streamlit_collection_picker_loads_next_acl_filtered_page(monkeypatch):
    from mn_molecular_design_studio import molecular_database

    records = [{"id": f"record-{i}", "name": f"Complex {i}",
                "representations": [{"id": f"structure-{i}", "label": "model"}]}
               for i in range(3)]

    class FakeDatabase:
        def list_structure_records(self, user, **kwargs):
            assert user == USER
            if kwargs.get("molecule_ids") is not None:
                selected = set(kwargs["molecule_ids"])
                return {"items": [row for row in records if row["id"] in selected], "next_offset": None}
            if int(kwargs.get("offset", 0)):
                return {"items": records[2:], "next_offset": None}
            return {"items": records[:2], "next_offset": 2}

    monkeypatch.setattr(molecular_database, "MolecularDatabase", FakeDatabase)
    app = AppTest.from_string(
        "from mn_molecular_design_studio.complex_analysis_collection_ui import render_collection_workspace\n"
        f"render_collection_workspace({USER!r})\n",
        default_timeout=10,
    ).run()
    assert not app.exception
    assert len(app.multiselect(key="complex-analysis-collection-records-" + hashlib.sha256(b"").hexdigest()[:10]).options) == 2

    next(button for button in app.button if button.label == "Load next Database structure page").click().run()

    assert not app.exception
    assert len(app.multiselect(key="complex-analysis-collection-records-" + hashlib.sha256(b"").hexdigest()[:10]).options) == 3


def test_ligand_pose_label_preserves_chain_qualified_and_residue_only_selectors():
    from mn_molecular_design_studio.complex_analysis_collection_ui import _ligand_pose_label

    assert _ligand_pose_label({"pose_chain": "L", "pose_residue": "L:1"}) == "L:1"
    assert _ligand_pose_label({"pose_chain": "L", "pose_residue": "1"}) == "L:1"
    assert _ligand_pose_label({"pose_chain": "", "pose_residue": "1"}) == "1"


def test_results_materialization_uses_shared_component_classifier_and_keeps_source_identity(tmp_path, monkeypatch):
    from mn_molecular_design_studio import complex_analysis_results_ui
    from mn_molecular_design_studio.complex_analysis_collection_ui import _materialize_results_sample

    structure = tmp_path / "rfd4-model.pdb"
    structure.write_bytes(_ppi())
    monkeypatch.setattr(complex_analysis_results_ui, "board_entries", lambda *_args, **_kwargs: [])
    monkeypatch.setattr(complex_analysis_results_ui, "_materialize_board_structure",
                        lambda *_args, **_kwargs: {"path": structure})
    candidate = {"candidate_id": "native-rfd4-17", "source_tool": "RFD4", "stage": "design",
                 "binder_chains": ["B"], "target_chains": ["A"],
                 "raw_metadata": {"generator_engine": "RFD4"}}
    job = {"task_group": "design", "run_id": "run-17", "tool": "rfd4", "node_name": "PC004033",
           "project_id": "project-1", "run_dir": str(tmp_path)}
    choice = {"remote": False, "job": job, "option": {"path": structure, "relative_path": "structures/model.pdb",
             "source_run_dir": str(tmp_path), "candidate": candidate}, "candidate": candidate,
             "candidate_id": "native-rfd4-17", "model_id": "0", "identity": "results:row-17",
             "relative_path": "structures/model.pdb", "label": "RFD4 · PC004033 · run-17 · native-rfd4-17"}

    sample = _materialize_results_sample(USER, choice)

    assert sample["kind"] == "protein-protein"
    assert sample["candidate_id"] == "native-rfd4-17"
    assert sample["generator_engine"] == "RFD4"
    assert sample["source_job"] == "design/run-17"
    assert sample["source_workstation"] == "PC004033"
    assert sample["model_id"] == "0"
    assert sample["digest"] == hashlib.sha256(_ppi()).hexdigest()


def test_comparison_projection_does_not_invent_missing_chain_roles(tmp_path):
    from mn_molecular_design_studio.complex_analysis_collection_ui import _comparison_sample

    structure = tmp_path / "unknown-roles.pdb"
    structure.write_bytes(_ppi())
    source = {"kind": "protein-protein", "identity": "results:unknown", "label": "Unknown roles",
              "path": structure, "candidate": {}, "context": {}, "job": {}, "components": [],
              "candidate_id": "candidate-1", "project_id": "project-1", "source_reference": "design/run-1",
              "relative_path": "complex.pdb", "model_id": "model-1", "digest": hashlib.sha256(_ppi()).hexdigest()}

    projected = _comparison_sample(USER, source)

    assert projected["target_chains"] == []
    assert projected["binder_chains"] == []
    assert projected["target_numbering"] == []
    assert projected["binder_numbering"] == []
