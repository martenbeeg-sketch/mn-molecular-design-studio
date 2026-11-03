from __future__ import annotations

import hashlib
import json

from streamlit.testing.v1 import AppTest

from mn_molecular_design_studio.complex_analysis_inspector import (
    author_residue_selector,
    inspector_evidence_rows,
    parse_author_residue_selector,
    protein_sequence_tracks,
    resolve_author_residue_selector,
)


def _atom(serial: int, residue: str, chain: str, number: int, insertion: str = " ") -> str:
    return (
        f"ATOM  {serial:5d}  CA  {residue:>3s} {chain}{number:4d}{insertion:1s}   "
        f"{float(serial):8.3f}{1.0:8.3f}{1.0:8.3f}{1.0:6.2f}{20.0:6.2f}           C"
    )


def test_author_residue_selector_round_trips_blank_and_insertion_codes():
    assert author_residue_selector("A", 10) == "A:10^_"
    assert author_residue_selector("A", 10, "B") == "A:10^B"
    assert parse_author_residue_selector("A:10^B") == {
        "chain_id": "A", "number": 10, "insertion_code": "B", "exact_insertion": True,
    }
    assert parse_author_residue_selector("A:10") == {
        "chain_id": "A", "number": 10, "insertion_code": None, "exact_insertion": False,
    }


def test_sequence_tracks_preserve_author_numbering_and_distinguish_insertions(tmp_path):
    path = tmp_path / "inserted.pdb"
    path.write_text("\n".join([
        _atom(1, "ALA", "A", 10),
        _atom(2, "GLY", "A", 10, "A"),
        _atom(3, "SER", "A", 11),
        "END",
        "",
    ]))

    tracks = protein_sequence_tracks(path, ["A"])

    assert tracks[0]["sequence"] == "AGS"
    assert [row["author_residue"] for row in tracks[0]["residues"]] == ["10", "10A", "11"]
    assert [row["selector"] for row in tracks[0]["residues"]] == ["A:10^_", "A:10^A", "A:11^_"]
    assert resolve_author_residue_selector(tracks, "A", 10) == ""
    assert resolve_author_residue_selector(tracks, "A", 10, "A") == "A:10^A"


def test_saved_interaction_rows_link_exact_residue_and_ligand_evidence(tmp_path):
    path = tmp_path / "complex.pdb"
    path.write_text("\n".join([
        _atom(1, "ALA", "A", 10, "A"),
        _atom(2, "GLY", "A", 11),
        "END",
        "",
    ]))
    tracks = protein_sequence_tracks(path, ["A"])

    rows = inspector_evidence_rows("plip", "interactions.csv", [{
        "interaction_type": "hydrogen bond",
        "protein_chain": "A",
        "protein_residue_number": "10",
        "protein_insertion_code": "A",
        "distance_angstrom": "2.7",
        "angle_degree": "145.0",
        "native_fields_json": '{"distance": 2.7, "angle": 145.0}',
    }], tracks, ligand_selector="L:1^_")

    assert rows[0]["selectors"] == ["A:10^A", "L:1^_"]
    assert rows[0]["label"] == "hydrogen bond · A:10^A"
    assert rows[0]["raw"]["distance_angstrom"] == "2.7"
    assert rows[0]["raw"]["native_fields_json"] == '{"distance": 2.7, "angle": 145.0}'


def test_biotite_contact_and_salt_bridge_artifacts_keep_each_residue_selector(tmp_path):
    path = tmp_path / "complex.pdb"
    path.write_text("\n".join([
        _atom(1, "ALA", "A", 10),
        _atom(2, "GLY", "B", 20, "C"),
        "END",
        "",
    ]))
    tracks = protein_sequence_tracks(path, ["A", "B"])

    contact = inspector_evidence_rows("biotite_interface", "contacts.json", [{
        "binder_residue": "A:10", "target_residue": "B:20C",
    }], tracks)
    bridge = inspector_evidence_rows("biotite_interface", "salt_bridges.json", [{
        "binder_residue": "A:10", "target_residue": "B:20C", "interaction_type": "Salt bridge",
    }], tracks)

    assert contact[0]["selectors"] == ["A:10^_", "B:20^C"]
    assert bridge[0]["selectors"] == ["A:10^_", "B:20^C"]
    assert bridge[0]["label"].startswith("Salt bridge")


def test_posebusters_evidence_highlights_the_saved_ligand_pose_only():
    rows = inspector_evidence_rows("posebusters", "full_checks.csv", [{
        "bond_lengths": True, "internal_steric_clash": True, "minimum_distance_to_protein": 2.1,
    }], [], ligand_selector="L:1^_")
    assert rows[0]["selectors"] == ["L:1^_"]
    assert rows[0]["raw"]["minimum_distance_to_protein"] == 2.1


def test_saved_evaluation_artifacts_are_bound_to_exact_structure_digest(tmp_path, monkeypatch):
    from mn_molecular_design_studio import complex_analysis_ui

    run_dir = tmp_path / "run"
    item_dir = run_dir / "artifacts" / "analysis" / "item_000000"
    item_dir.mkdir(parents=True)
    (run_dir / "artifacts" / "analysis" / "batch_manifest.json").write_text(json.dumps({
        "items": [{"id": "item_000000", "candidate_id": "database-candidate"}],
    }))
    (item_dir / "interactions.csv").write_text(
        "interaction_type,protein_chain,protein_residue_number,distance_angstrom\n"
        "hydrogen bond,A,10,2.7\n"
    )
    source = tmp_path / "structure.cif"
    source.write_text("structure bytes")
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    monkeypatch.setattr(complex_analysis_ui, "_authorized_evaluation_run_dir", lambda *_args: run_dir)
    observation = {
        "run_id": "analysis-run",
        "context": {
            "task_group": "analysis", "candidate_id": "database-candidate",
            "source_structure_id": "structure-id", "source_structure_digest": digest,
            "candidate": {"raw_metadata": {"complex_analysis_selection": {"structure_digest": digest}}},
        },
    }
    structure = {"id": "structure-id"}

    artifacts = complex_analysis_ui._saved_evaluation_artifacts(
        {"id": "user-id"}, observation, structure, digest,
    )
    assert [item["name"] for item in artifacts] == ["interactions.csv"]
    assert artifacts[0]["rows"][0]["protein_residue_number"] == "10"
    assert complex_analysis_ui._saved_evaluation_artifacts(
        {"id": "user-id"}, observation, structure, "changed-digest",
    ) == []


def test_inspector_findings_keep_residue_evidence_and_leave_summary_in_measurements(tmp_path, monkeypatch):
    from mn_molecular_design_studio import complex_analysis_ui

    source = tmp_path / "complex.pdb"
    source.write_text("\n".join([_atom(1, "ALA", "A", 10), "END", ""]))
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    observation = {
        "id": "observation-id", "namespace": "plip", "run_id": "run-id",
        "metrics": {"plip_interaction_count": 1},
        "context": {"representation_id": "structure-id", "structure_digest": digest},
    }
    descriptor = {"id": "evaluation-id", "evaluator": "plip", "representation_id": "structure-id"}
    artifacts = [
        {"name": "report.json", "path": source, "rows": [{"interaction_count": 1}]},
        {"name": "interactions.csv", "path": source, "rows": [{
            "interaction_type": "hydrogen bond", "protein_chain": "A",
            "protein_residue_number": 10, "distance_angstrom": 2.7,
        }]},
    ]
    monkeypatch.setattr(complex_analysis_ui, "_saved_evaluation_artifacts", lambda *_args: artifacts)

    data = complex_analysis_ui._database_inspector_data(
        {"id": "user-id"}, {"id": "molecule-id"},
        {"id": "structure-id", "context": {"structure_digest": digest}},
        source,
        [{"type": "protein", "source_chain": "A"}],
        [(observation, descriptor)],
    )

    assert [row["artifact"] for row in data["evidence"]] == ["interactions.csv"]
    assert data["evidence"][0]["selectors"] == ["A:10^_"]
    assert [artifact["name"] for artifact in data["evaluations"][0]["artifacts"]] == [
        "report.json", "interactions.csv",
    ]


def test_evaluation_artifact_resolution_rechecks_job_access(monkeypatch):
    from mn_molecular_design_studio import auth, job_catalog
    from mn_molecular_design_studio.complex_analysis_ui import _authorized_evaluation_run_dir

    calls = []
    monkeypatch.setattr(job_catalog, "get_job_access_record", lambda *_args: {"owner_id": "another-user"})
    monkeypatch.setattr(auth, "can_access_job", lambda *_args: False)
    monkeypatch.setattr(job_catalog, "search_jobs", lambda **_kwargs: calls.append("searched") or {"items": []})

    assert _authorized_evaluation_run_dir({"id": "user-id"}, {"run_id": "run-id", "context": {"task_group": "analysis"}}) is None
    assert calls == []


def test_measurement_panel_displays_saved_units_methods_and_evaluation_parameters():
    item = {
        "evaluator": "ipsae",
        "evaluation": {"id": "evaluation-1", "run_id": "run-1"},
        "observation": {
            "id": "observation-1",
            "metrics": {"native_score": 0.5},
            "normalized_output": {
                "metrics": {"confidence.interface.ipsae": 0.5},
                "provenance": {"confidence.interface.ipsae": {
                    "unit": "Å", "scope": "directional interface",
                    "method": "saved ipSAE definition", "source_tool": "ipsae",
                }},
            },
            "context": {"settings": {"interface_threshold_angstrom": 12.0}},
        },
        "artifacts": [],
    }
    script = f'''\
from mn_molecular_design_studio.complex_analysis_ui import _render_inspector_measurements
_render_inspector_measurements([{item!r}])
'''

    app = AppTest.from_string(script, default_timeout=10).run()

    assert not app.exception
    projection = app.dataframe[1].value.iloc[0].to_dict()
    assert projection == {
        "Registered measure": "confidence.interface.ipsae",
        "Stored projection": 0.5,
        "Unit": "Å",
        "Scope": "directional interface",
        "Method": "saved ipSAE definition",
        "Source tool": "ipsae",
    }
    provenance = json.loads(app.json[0].value)
    assert provenance["evaluation_context"]["settings"] == {"interface_threshold_angstrom": 12.0}
