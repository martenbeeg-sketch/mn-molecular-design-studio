from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from mn_molecular_design_studio import molecule_imports
from mn_molecular_design_studio.molecular_database import MolecularDatabase


USER = {"id": "a" * 32, "role": "user", "_cluster_project_ids": []}
_FIXTURE = Path(__file__).parents[1] / "mn_molecular_design_studio/assets/generators/pxdesign/5o45.cif"


def _controlled_real_structure_pair(tmp_path: Path, monkeypatch):
    """Use experimental 5O45 coordinates with a documented exact construct.

    Reference is the complete tracked 5O45 coordinate file, which is a partial
    PD-L1 structural construct. Source is a 61-residue
    author-numbered segment of chain A (45–105) plus the original chain-B
    partner, displaced by a known rigid transform. This is a validation fixture
    only; it is not a production Database or biological reference assignment.
    """
    from Bio.PDB import MMCIFParser, PDBIO, Select

    workspace = tmp_path / "isolated-studio-home"
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_HOME", str(workspace))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(workspace))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR", str(workspace / "runs"))

    structure = MMCIFParser(QUIET=True).get_structure("controlled-reference", str(_FIXTURE))
    model = next(structure.get_models())
    rotation = np.asarray([[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    translation = np.asarray([18.0, -7.0, 4.0])
    for atom in model.get_atoms():
        residue = atom.get_parent()
        if residue.get_parent().id in {"A", "B"} and residue.id[0] != "W":
            atom.set_coord(atom.get_coord() @ rotation + translation)

    class SourceConstruct(Select):
        def accept_chain(self, chain):
            return chain.id in {"A", "B"}

        def accept_residue(self, residue):
            if residue.id[0] == "W":
                return False
            if residue.get_parent().id == "A":
                return residue.id[0] == " " and 45 <= int(residue.id[1]) <= 105
            return True

    source_file = tmp_path / "controlled-source.pdb"
    writer = PDBIO()
    writer.set_structure(structure)
    writer.save(str(source_file), SourceConstruct())

    database = MolecularDatabase(workspace / "database")
    reference_record = molecule_imports.parse_structure(_FIXTURE.read_bytes(), _FIXTURE.name)[0]
    source_record = molecule_imports.parse_structure(source_file.read_bytes(), source_file.name)[0]
    reference_import = database.import_records(USER, [reference_record], name="Controlled 5O45 reference")
    source_import = database.import_records(USER, [source_record], name="Controlled 5O45 construct complex")
    reference = database.list_molecules(USER, molecule_ids=reference_import["molecule_ids"], collection="all")["items"][0]
    source = database.list_molecules(USER, molecule_ids=source_import["molecule_ids"], collection="all")["items"][0]
    reference_structure = reference["representations"][0]
    source_structure = source["representations"][0]
    association = database.save_complex_reference_context(
        USER, source_molecule_id=source["id"], source_structure_id=source_structure["id"],
        source_chain="A", reference_molecule_id=reference["id"],
        reference_structure_id=reference_structure["id"], reference_chain="A",
    )
    assert association["status"] == "validated"
    assert len(association["mapping"]["residue_mapping"]) == 61
    selection = {
        "molecule_id": source["id"], "structure_id": source_structure["id"],
        "group_id": source_import["group_id"],
        "model_id": (source_structure.get("context") or {}).get("model_id") or source_structure["id"],
        "binder_chains": [], "target_chains": [], "source_component_chain": "A",
    }
    return database, source, source_structure, reference, reference_structure, association, selection


def test_real_reference_placement_and_reference_sequence_track_reopen(tmp_path, monkeypatch):
    from Bio.PDB import MMCIFParser, PDBParser

    from mn_molecular_design_studio.complex_analysis_sources import build_database_analysis_source
    from mn_molecular_design_studio.complex_analysis_reference_context import reference_component_track
    from mn_molecular_design_studio.engine.workflow_engine.workflows.analysis import run_reference_placement_contract
    from mn_molecular_design_studio.evaluations import evaluation_descriptor, file_digest
    from mn_molecular_design_studio.molecular_database_ingest import ingest_run
    from mn_molecular_design_studio.workspace import mark_studio_job

    database, source, source_representation, reference, reference_representation, association, selection = (
        _controlled_real_structure_pair(tmp_path, monkeypatch)
    )
    adapter = build_database_analysis_source(
        USER, selection, workflow_id="analysis.reference_placement", project_id=None, database=database,
    )
    run_dir = run_reference_placement_contract(
        Path(adapter["source_run_dir"]), Path(adapter["candidates_jsonl"]),
        candidate_ids=[adapter["candidate_id"]],
        selected_structure_path=Path(adapter["selected_structure_path"]),
        selected_model_id=adapter["selected_model_id"], source_component_chain="A",
    )
    result_path = run_dir / "artifacts/analysis/reference_placement/reference_placement.json"
    result = json.loads(result_path.read_text())
    assert result["aligned_pair_count"] == 61
    assert result["source_coverage"] == 1.0
    assert result["reference_coverage"] == 61 / 129
    assert result["transformation_direction"] == "selected complex/source coordinates -> associated reference coordinates"
    assert result["source"]["structure_id"] == source_representation["id"]
    assert result["reference"]["structure_id"] == reference_representation["id"]
    assert result["reference"]["molecule_id"] == reference["id"]
    assert result["reference"]["structure_digest"] == association["reference_structure_digest"]
    assert result["rmsd_angstrom"] < 0.002

    # Independently verify the Mol* column-vector matrix against original
    # target and partner atoms. Database coordinates are not edited.
    source_path = database.structure_file(USER, source_representation["id"])
    reference_path = database.structure_file(USER, reference_representation["id"])
    parser = PDBParser(QUIET=True)
    source_model = next(parser.get_structure("placed-source", str(source_path)).get_models())
    reference_model = next(MMCIFParser(QUIET=True).get_structure("original-reference", str(reference_path)).get_models())
    matrix = np.asarray(result["matrix_4x4"], dtype=float)
    mapped_target_errors = []
    for pair in result["mapped_residue_pairs"]:
        source_info, reference_info = pair["source"], pair["reference"]
        source_id = (" ", int(source_info["number"]), str(source_info.get("insertion") or " "))
        reference_id = (" ", int(reference_info["number"]), str(reference_info.get("insertion") or " "))
        source_ca = np.asarray(source_model["A"][source_id]["CA"].coord, dtype=float)
        reference_ca = np.asarray(reference_model["A"][reference_id]["CA"].coord, dtype=float)
        transformed = matrix[:3, :3] @ source_ca + matrix[:3, 3]
        mapped_target_errors.append(float(np.linalg.norm(transformed - reference_ca)))
    assert np.sqrt(np.mean(np.square(mapped_target_errors))) < 0.002

    def partner_atoms(model):
        return {
            (residue.resname, int(residue.id[1]), str(residue.id[2] or "").strip(), atom.name.strip(), atom.element.strip().upper()): atom.coord
            for residue in model["B"] for atom in residue
            if atom.element.strip().upper() not in {"H", "D"}
        }

    source_partner, reference_partner = partner_atoms(source_model), partner_atoms(reference_model)
    common = sorted(source_partner.keys() & reference_partner.keys())
    transformed_partner_errors = [
        float(np.linalg.norm(matrix[:3, :3] @ source_partner[key] + matrix[:3, 3] - reference_partner[key]))
        for key in common
    ]
    assert len(common) >= 50
    assert np.sqrt(np.mean(np.square(transformed_partner_errors))) < 0.003

    source_digest_before = file_digest(source_path)
    reference_digest_before = file_digest(reference_path)
    mark_studio_job({"run_dir": str(run_dir)}, project_id=None, owner_id=USER["id"], allow_private_source=True)
    assert ingest_run(run_dir, database=database, force=True) == 1
    reopened = database.list_molecules(USER, molecule_ids=[source["id"]], collection="all")["items"][0]
    observations = [item for item in reopened["observations"] if item.get("namespace") == "reference_placement"]
    assert len(observations) == 1
    assert evaluation_descriptor(observations[0])["evaluator"] == "reference_placement"
    assert observations[0]["context"]["representation_id"] == source_representation["id"]
    assert file_digest(source_path) == source_digest_before == result["source"]["structure_digest"]
    assert file_digest(reference_path) == reference_digest_before == result["reference"]["structure_digest"]

    reopened_association = database.get_complex_reference_context(
        USER, source_molecule_id=source["id"], source_structure_id=source_representation["id"], source_chain="A",
    )
    from mn_molecular_design_studio.complex_analysis_ui import _reference_placement_match
    match = _reference_placement_match(USER, source, source_representation, reopened_association, selection["model_id"])
    assert match is not None
    assert match["result"]["matrix_4x4"] == result["matrix_4x4"]

    track = reference_component_track(reference_path, "A", result["mapped_residue_pairs"])
    assert track["status"] == "available"
    assert len(track["sequence"]) == 129
    assert track["counts"] == {"mapped": 61, "present": 129, "omitted": 68, "unresolved": 0}
    assert track["residues"][0]["status"] == "omitted"
    assert track["residues"][44]["status"] == "mapped"
    assert track["residues"][44]["selector"].startswith("A:")
    assert track["residues"][-1]["status"] == "omitted"


def test_full_sequence_track_distinguishes_unresolved_and_insertion_codes(tmp_path):
    from mn_molecular_design_studio.complex_analysis_reference_context import reference_component_track

    path = tmp_path / "partial-reference.pdb"
    path.write_text(
        "SEQRES   1 A    3  ALA GLY SER\n"
        "ATOM      1  CA  ALA A  10       0.000   0.000   0.000  1.00 20.00           C  \n"
        "ATOM      2  CA  SER A  12A      2.000   1.000   0.000  1.00 20.00           C  \n"
        "TER\nEND\n"
    )
    correspondence = [
        {"reference": {"number": 10, "insertion": "", "amino_acid": "A"}},
        {"reference": {"number": 12, "insertion": "A", "amino_acid": "S"}},
    ]
    track = reference_component_track(path, "A", correspondence)

    assert track["status"] == "available"
    assert track["sequence"] == "AGS"
    assert [row["status"] for row in track["residues"]] == ["mapped", "unresolved", "mapped"]
    assert track["residues"][0]["author_residue"] == "10"
    assert track["residues"][1]["author_residue"] == "seq#2"
    assert track["residues"][2]["author_residue"] == "12A"
    assert track["residues"][2]["selector"] == "A:12^A"
    assert track["counts"] == {"mapped": 2, "present": 2, "omitted": 0, "unresolved": 1}


def test_reference_sequence_button_sends_exact_author_selector_to_paired_viewer():
    from streamlit.testing.v1 import AppTest

    track = {
        "chain_id": "R",
        "residues": [{
            "sequence_index": 0, "amino_acid": "G", "author_residue": "10A",
            "status": "omitted", "selector": "R:10^A",
        }],
    }
    script = f'''\
from mn_molecular_design_studio.complex_analysis_ui import _render_reference_context_track
_render_reference_context_track({track!r}, "paired-viewer")
'''
    app = AppTest.from_string(script, default_timeout=20).run()
    assert not app.exception
    app.button[0].click().run()
    assert not app.exception
    assert app.session_state["paired-viewer-reference-highlights"] == ["R:10^A"]
