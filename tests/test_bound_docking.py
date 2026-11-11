import json
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem
from streamlit.testing.v1 import AppTest

from mn_molecular_design_studio import bound_docking
from mn_molecular_design_studio.molecular_database import MolecularDatabase

USER = {"id": "a" * 32, "role": "admin"}


def pdb_atom(serial, name, residue, chain, number, xyz, element="C", record="ATOM"):
    x, y, z = xyz
    return f"{record:6}{serial:5d} {name:4s} {residue:3s} {chain}{number:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2}"


@pytest.fixture
def bound(tmp_path, monkeypatch):
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_HOME", str(tmp_path))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(tmp_path))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR", str(tmp_path / "runs"))
    receptor = tmp_path / "receptor.pdb"
    receptor.write_text("\n".join([pdb_atom(i, " CA ", "ALA", "L", i, (i, 1, 2)) for i in range(1, 6)]) + "\nEND\n")
    pose = tmp_path / "pose.pdbqt"
    atoms = [pdb_atom(i, symbol, "UNL", " ", 1, xyz, record="ATOM")[:66] + "     0.000 " + symbol
             for i, symbol, xyz in [(1, "O", (8, 3, 4)), (2, "C", (7, 3, 4)), (3, "C", (6, 3, 4))]]
    pose.write_text("\n".join(["MODEL 1", "REMARK SMILES CCO", "REMARK SMILES IDX 1 3 2 2 3 1", *atoms,
                                  "ENDMDL", "MODEL 2", *atoms, "ENDMDL"]))
    db = MolecularDatabase()
    linked = bound_docking.register_receptor(db, USER, receptor)
    result = db.import_records(USER, [{"name": "Ethanol", "components": [{"type": "small_molecule", "chain": "L", "smiles": "CCO"}]}])
    structure = db.attach_structure(result["group_id"], result["molecule_ids"][0], {"path": str(pose), "label": "Docked pose"},
        context={**linked, "database_structure_policy": "best_docking_pose", "run_id": "docking-run"})
    return db, result["molecule_ids"][0], structure, linked, receptor, pose


def test_compound_and_target_remain_separate_and_share_portable_receptor(bound):
    db, molecule, _, linked, receptor, _ = bound
    second = bound_docking.register_receptor(db, USER, receptor)
    assert second == linked
    rows = db.list_molecules(USER)["items"]
    assert {r["type"] for r in rows} == {"small_molecule", "protein"}
    assert len(rows) == 2
    assert next(r for r in rows if r["id"] == molecule)["components"] == [
        {"type": "small_molecule", "chain": "L", "role": "", "smiles": "CCO"}]
    with db.connection() as c:
        assert c.execute("SELECT count(*) FROM molecule_structures WHERE molecule_id=?", (linked["target_molecule_id"],)).fetchone()[0] == 1
        context = json.loads(c.execute("SELECT context FROM molecule_structures WHERE id=?", (linked["target_structure_id"],)).fetchone()[0])
        assert context["portable"] is True


def test_md_staging_keeps_mapped_coordinates_atom_order_and_database_identities(bound):
    db, molecule, structure, linked, receptor, _ = bound
    result = bound_docking.prepare_bound_inputs(db, USER, molecule, structure_id=structure)
    sdf = next(Chem.SDMolSupplier(result["ligand_sdf"], removeHs=False))
    assert Chem.MolToSmiles(sdf) == "CCO"
    np.testing.assert_allclose(sdf.GetConformer().GetPositions(), [[6, 3, 4], [7, 3, 4], [8, 3, 4]])
    text = Path(result["complex_pdb"]).read_text()
    atoms = [line for line in text.splitlines() if line.startswith("HETATM")]
    assert [line[76:78].strip() for line in atoms] == ["C", "C", "O"]
    assert all(line[21] == "Z" for line in atoms)  # Receptor already uses L.
    assert result["ligand_key"] == "LIG|Z|1|_"
    assert all(line in text for line in receptor.read_text().splitlines() if line.startswith("ATOM"))
    assert len(db.list_molecules(USER)["items"]) == 2
    assert result == bound_docking.prepare_bound_inputs(db, USER, molecule, structure_id=structure)


def test_md_rejects_ambiguous_mapping_instead_of_inventing_coordinates(bound):
    _, _, _, _, _, pose = bound
    pose.write_text(pose.read_text().replace("1 3 2 2 3 1", "1 3 2 2"))
    with pytest.raises(ValueError, match="atom map"):
        bound_docking.pose_molecule(pose, {}, "CCO")


def test_md_rejects_pose_chemistry_that_does_not_match_saved_ligand_identity(bound):
    _, _, _, _, _, pose = bound
    pose.write_text(pose.read_text().replace("REMARK SMILES CCO", "REMARK SMILES C1CO1"))
    with pytest.raises(ValueError, match="does not match the saved Database ligand identity"):
        bound_docking.pose_molecule(pose, {}, "CCO")


def test_chemically_annotated_pose_records_source_atom_coordinate_mapping(tmp_path):
    molecule = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    conformer = Chem.Conformer(molecule.GetNumAtoms())
    for index in range(molecule.GetNumAtoms()):
        conformer.SetAtomPosition(index, (float(index), 2.0, 3.0))
    molecule.AddConformer(conformer)
    pose = tmp_path / "pose.sdf"
    with Chem.SDWriter(str(pose)) as writer:
        writer.write(molecule)

    reconstructed, provenance = bound_docking.pose_molecule(
        pose, {}, "CCO", return_provenance=True)

    assert reconstructed.GetNumAtoms() == 3
    assert provenance["atom_mapping"] == [
        {"ligand_atom_index": 0, "pose_atom_index": 0, "element": "C"},
        {"ligand_atom_index": 1, "pose_atom_index": 1, "element": "C"},
        {"ligand_atom_index": 2, "pose_atom_index": 2, "element": "O"},
    ]
    np.testing.assert_allclose(reconstructed.GetConformer().GetPositions(), [[0, 2, 3], [1, 2, 3], [2, 2, 3]])


def test_pose_receptor_digest_must_match_before_complex_reconstruction(bound):
    db, molecule, structure_id, _, _, _ = bound
    with db.connection() as connection:
        row = connection.execute("SELECT context FROM molecule_structures WHERE id=?", (structure_id,)).fetchone()
        context = json.loads(row[0])
        context["target_structure_digest"] = "0" * 64
        connection.execute("UPDATE molecule_structures SET context=? WHERE id=?", (json.dumps(context), structure_id))
    with pytest.raises(ValueError, match="different coordinate-frame digests"):
        bound_docking.prepare_bound_inputs(db, USER, molecule, structure_id=structure_id)


def test_pose_without_any_saved_receptor_digest_fails_coordinate_frame_validation(bound):
    db, molecule, structure_id, linked, _, _ = bound
    with db.connection() as connection:
        pose_row = connection.execute("SELECT context FROM molecule_structures WHERE id=?", (structure_id,)).fetchone()
        pose_context = json.loads(pose_row[0])
        pose_context.pop("target_structure_digest", None)
        pose_context.pop("target_variant_sha256", None)
        connection.execute("UPDATE molecule_structures SET context=? WHERE id=?", (json.dumps(pose_context), structure_id))
        target_row = connection.execute("SELECT context FROM molecule_structures WHERE id=?", (linked["target_structure_id"],)).fetchone()
        target_context = json.loads(target_row[0])
        target_context.pop("portable_sha256", None)
        connection.execute("UPDATE molecule_structures SET context=? WHERE id=?",
                           (json.dumps(target_context), linked["target_structure_id"]))

    with pytest.raises(ValueError, match="no digest to verify their coordinate frame"):
        bound_docking.prepare_bound_inputs(db, USER, molecule, structure_id=structure_id)


def test_bound_complex_is_attached_to_same_record_with_chemical_identity_and_provenance(bound):
    db, molecule, source_pose_id, linked, _, pose = bound

    result = bound_docking.register_bound_complex(db, USER, molecule, structure_id=source_pose_id)

    rows = db.list_molecules(USER, molecule_ids=[molecule], collection="all")["items"]
    assert len(rows) == 1
    record = rows[0]
    assert {item["id"] for item in record["representations"]} == {source_pose_id, result["complex_structure_id"]}
    complex_structure = next(item for item in record["representations"] if item["id"] == result["complex_structure_id"])
    context = complex_structure["context"]
    assert context["database_structure_policy"] == "bound_docking_complex"
    assert context["source_pose_structure_id"] == source_pose_id
    assert context["source_pose_digest"] == __import__("hashlib").sha256(pose.read_bytes()).hexdigest()
    assert context["source_receptor_structure_id"] == linked["target_structure_id"]
    assert context["target_structure_digest"] == linked["target_structure_digest"]
    assert context["coordinate_frame_id"] == "receptor-sha256:" + linked["target_structure_digest"]
    assert context["ligand_selection"] == {"chain": "Z", "residue_key": "Z:1", "ccd": "LIG", "smiles": "CCO"}
    assert context["ligand_isomeric_smiles"] == "CCO"
    assert context["pose_atom_mapping"] == [
        {"ligand_atom_index": 0, "pose_atom_serial": 3, "element": "C"},
        {"ligand_atom_index": 1, "pose_atom_serial": 2, "element": "C"},
        {"ligand_atom_index": 2, "pose_atom_serial": 1, "element": "O"},
    ]
    complex_path = db.structure_file(USER, result["complex_structure_id"])
    assert complex_path.read_bytes() == Path(result["complex_pdb"]).read_bytes()
    assert complex_structure["context"]["portable_sha256"] == __import__("hashlib").sha256(complex_path.read_bytes()).hexdigest()

    from mn_molecular_design_studio.complex_analysis_sources import authorize_database_structure
    source = authorize_database_structure(USER, {
        "molecule_id": molecule, "structure_id": result["complex_structure_id"],
        "group_id": result["group_id"], "model_id": result["complex_structure_id"],
        "ligand_selection": context["ligand_selection"],
    }, project_id="", database=db, workflow_id="analysis.plip")
    assert source["structure_digest"] == context["complex_structure_digest"]
    assert source["ligand_selection"] == context["ligand_selection"]


def test_bound_complex_registration_is_idempotent_and_respects_group_write_access(bound):
    db, molecule, source_pose_id, _, _, _ = bound
    first = bound_docking.register_bound_complex(db, USER, molecule, structure_id=source_pose_id)
    second = bound_docking.register_bound_complex(db, USER, molecule, structure_id=source_pose_id)
    assert first["complex_structure_id"] == second["complex_structure_id"]
    assert len(db.list_molecules(USER, molecule_ids=[molecule], collection="all")["items"][0]["representations"]) == 2
    with pytest.raises(PermissionError):
        bound_docking.register_bound_complex(db, {"id": "b" * 32}, molecule, structure_id=source_pose_id)


def test_registered_bound_complex_is_selectable_and_reopens_in_database_analysis(bound, monkeypatch):
    db, molecule, source_pose_id, _, _, _ = bound
    result = bound_docking.register_bound_complex(db, USER, molecule, structure_id=source_pose_id)
    from mn_molecular_design_studio import complex_analysis_ui, molecular_database
    monkeypatch.setattr(molecular_database, "MolecularDatabase", lambda: db)
    monkeypatch.setattr(complex_analysis_ui, "_show_structure", lambda *args, **kwargs: None)

    opened = {"molecule_id": molecule, "structure_id": result["complex_structure_id"]}
    app = AppTest.from_string(
        "from mn_molecular_design_studio.complex_analysis_ui import _render_database_source\n"
        f"_render_database_source({USER!r}, {opened!r})"
    ).run()

    assert not app.exception
    structures = next(widget for widget in app.selectbox if widget.label == "Structure / model / pose")
    assert "Docking pose · receptor–ligand complex" in structures.options
    assert structures.value == result["complex_structure_id"]
    assert any(widget.label == "Ligand residue" for widget in app.selectbox)


@pytest.mark.parametrize("workflow_id", ["analysis.plip", "analysis.pandamap", "analysis.posebusters"])
def test_registered_bound_complex_uses_existing_ligand_workflow_submission_contract(bound, monkeypatch, workflow_id):
    db, molecule, source_pose_id, _, _, _ = bound
    result = bound_docking.register_bound_complex(db, USER, molecule, structure_id=source_pose_id)
    from mn_molecular_design_studio import complex_analysis_ui, molecular_database
    monkeypatch.setattr(molecular_database, "MolecularDatabase", lambda: db)
    monkeypatch.setattr(complex_analysis_ui, "_show_structure", lambda *args, **kwargs: None)
    monkeypatch.setattr(complex_analysis_ui, "_matching_database_analysis_jobs", lambda *args, **kwargs: [])
    opened = {"molecule_id": molecule, "structure_id": result["complex_structure_id"]}
    script = f'''\
import streamlit as st
from mn_molecular_design_studio.complex_analysis_ui import _render_database_source
def submit(schema, workflow, parameters, project_id, settings, title):
    st.session_state["submitted"] = {{"workflow": workflow, "parameters": parameters, "project_id": project_id}}
_render_database_source({USER!r}, {opened!r}, workflow_schema=lambda workflow_id: {{"parameters": []}}, submit=submit)
'''
    app = AppTest.from_string(script).run()
    key = __import__("hashlib").sha256(f"database:{molecule}:{result['complex_structure_id']}".encode()).hexdigest()[:14]
    run_button = next(item for item in app.button
                      if f"run-{workflow_id.replace('.', '-')}" in str(item.key))
    assert key in run_button.key
    run_button.click().run()
    assert not app.exception
    submitted = app.session_state["submitted"]
    assert submitted["workflow"] == workflow_id
    assert submitted["project_id"] == ""
    assert submitted["parameters"]["database_source"]["molecule_id"] == molecule
    assert submitted["parameters"]["database_source"]["structure_id"] == result["complex_structure_id"]
    assert submitted["parameters"]["ligand_selection"] == {
        "chain": "Z", "residue_key": "Z:1", "ccd": "LIG", "smiles": "CCO",
    }


@pytest.mark.parametrize(("case", "expected_ingestions"), [
    ("unannotated", 1),
    ("correct-smiles", 1),
    ("conflicting-smiles", 0),
    ("receptor-substitution", 0),
    ("ligand-element-substitution", 0),
    ("ligand-geometry-substitution", 0),
    ("wrong-ligand-residue", 0),
])
def test_ligand_owned_bound_complex_plip_checks_selected_representation(
        bound, tmp_path, case, expected_ingestions):
    db, molecule, source_pose_id, _, _, _ = bound
    registered = bound_docking.register_bound_complex(db, USER, molecule, structure_id=source_pose_id)
    canonical = db.list_molecules(USER, molecule_ids=[molecule], collection="all")["items"][0]
    assert canonical["components"] == [
        {"type": "small_molecule", "chain": "L", "role": "", "smiles": "CCO"}
    ]
    complex_structure = next(item for item in canonical["representations"]
                             if item["id"] == registered["complex_structure_id"])
    selection = complex_structure["context"]["ligand_selection"]
    assert selection == {"chain": "Z", "residue_key": "Z:1", "ccd": "LIG", "smiles": "CCO"}
    from mn_molecular_design_studio.complex_analysis_sources import build_database_analysis_source
    adapter = build_database_analysis_source(USER, {
        "molecule_id": molecule, "structure_id": complex_structure["id"],
        "group_id": registered["group_id"],
        "model_id": complex_structure.get("context", {}).get("model_id") or complex_structure["id"],
        "ligand_selection": selection,
    }, workflow_id="analysis.plip", project_id="", database=db)

    original = Path(adapter["selected_structure_path"]).read_text()
    receptor_names = {"A": "ALA", "G": "GLY"}
    calculation_lines = ["REMARK PLIP evaluator-prepared complex"]
    geometry_atom_shifted = False
    for line in original.splitlines():
        if case == "receptor-substitution" and line.startswith("ATOM  ") and line[21:22] == "L":
            line = line[:17] + receptor_names["G"] + line[20:]
        elif case == "ligand-element-substitution" and line.startswith("HETATM") and line[21:22] == "Z" \
                and line[76:78].strip() == "O":
            line = line[:76] + " N" + line[78:]
        elif (case == "ligand-geometry-substitution" and not geometry_atom_shifted
              and line.startswith("HETATM") and line[21:22] == "Z"):
            x = float(line[30:38]) + 4.0
            line = line[:30] + f"{x:8.3f}" + line[38:]
            geometry_atom_shifted = True
        calculation_lines.append(line)
    calculation_text = "\n".join(calculation_lines) + "\n"
    if case in {"unannotated", "correct-smiles"}:
        assert calculation_text == (
            "REMARK PLIP evaluator-prepared complex\n"
            + "\n".join(original.splitlines()) + "\n"
        )
    calculation = tmp_path / f"plip-{case}.pdb"
    calculation.write_text(calculation_text)
    candidate = json.loads(Path(adapter["candidates_jsonl"]).read_text().splitlines()[0])
    candidate.update({
        "candidate_id": candidate["candidate_id"] + "_plip_selected_representation",
        "source_tool": "plip", "tool": "plip", "stage": "analysis", "complex_pdb": str(calculation),
        "metrics": {"plip_success": True, "plip_interaction_count": 2},
    })
    raw = candidate["raw_metadata"]
    raw["analysis_subject"] = "ligand"
    raw["complex_analysis_selection"]["source_structure_digest"] = adapter["source"]["structure_digest"]
    raw["complex_analysis_selection"]["structure_digest"] = adapter["source"]["structure_digest"]
    raw["complex_analysis_selection"]["ligand_selection"] = selection
    raw.pop("ligand_smiles", None)
    if case == "correct-smiles":
        raw["ligand_smiles"] = ["CCO"]
    elif case == "conflicting-smiles":
        raw["ligand_smiles"] = ["COC"]
    elif case == "wrong-ligand-residue":
        raw["complex_analysis_selection"]["ligand_selection"] = {
            **selection, "residue_key": "Z:2",
        }
    raw["structure_digest"] = __import__("hashlib").sha256(calculation.read_bytes()).hexdigest()

    run = tmp_path / "runs" / "analysis" / f"plip-ligand-only-{case}"
    candidate_path = run / "artifacts/normalized_candidates/candidates.jsonl"
    candidate_path.parent.mkdir(parents=True)
    candidate_path.write_text(json.dumps(candidate) + "\n")
    (run / "molecular_design_studio.json").write_text(json.dumps({
        "application": "mn-molecular-design-studio", "studio_owner_id": USER["id"], "project_id": "",
    }))
    (run / "metadata.json").write_text(json.dumps({
        "status": "completed", "tool": "plip", "workflow": "analysis.plip",
        "upstream_source_kind": "database", "database_collection_policy": "automatic",
    }))
    (run / "input.json").write_text(json.dumps({"tool": "plip", "workflow": "analysis.plip",
                                                  "params": {"ligand_selection": selection}}))
    (run / "result.json").write_text(json.dumps({"success": True, "execution_status": "completed"}))

    from mn_molecular_design_studio.molecular_database_ingest import ingest_run
    assert ingest_run(run, database=db, force=True) == expected_ingestions
    after = db.list_molecules(USER, molecule_ids=[molecule], collection="all")["items"][0]
    evaluations = [item for item in after["observations"] if item.get("namespace") == "plip"]
    assert len(evaluations) == expected_ingestions
    assert after["components"] == canonical["components"]
    assert after["id"] == molecule
    all_records = db.list_molecules(USER, collection="all")["items"]
    ethanol_records = [item for item in all_records if item["components"] == canonical["components"]]
    assert [item["id"] for item in ethanol_records] == [molecule]
    if expected_ingestions:
        assert evaluations[0]["context"]["representation_id"] == complex_structure["id"]
        assert ingest_run(run, database=db, force=True) == 1
        replayed = db.list_molecules(USER, molecule_ids=[molecule], collection="all")["items"][0]
        assert len([item for item in replayed["observations"] if item.get("namespace") == "plip"]) == 1



def test_compare_saved_complexes_resolves_attached_docking_complex_as_protein_ligand(bound, monkeypatch):
    db, molecule, source_pose_id, _, _, _ = bound
    result = bound_docking.register_bound_complex(db, USER, molecule, structure_id=source_pose_id)
    from mn_molecular_design_studio import molecular_database, complex_analysis_comparison_ui
    monkeypatch.setattr(molecular_database, "MolecularDatabase", lambda: db)
    observed = []
    monkeypatch.setattr(complex_analysis_comparison_ui.st, "write", lambda value, *args, **kwargs: observed.append(value))
    script = f'''\
import streamlit as st
from mn_molecular_design_studio.complex_analysis_comparison_ui import _comparison_slot
st.session_state["complex-analysis-compare-a-source"] = "Database"
st.session_state["complex-analysis-compare-a-molecule"] = {molecule!r}
st.session_state["complex-analysis-compare-a-database-structure"] = {result['complex_structure_id']!r}
selected = _comparison_slot({USER!r}, "a")
st.write(selected["complex_type"])
st.write(selected["ligand_selection"])
st.write(selected["coordinate_frame_id"])
'''
    app = AppTest.from_string(script).run()

    assert not app.exception
    assert "protein-ligand" in observed
    assert {"chain": "Z", "residue_key": "Z:1", "ccd": "LIG", "smiles": "CCO"} in observed
    assert "receptor-sha256:" + result["target_structure_digest"] in observed


def test_receptor_resolves_via_authorized_database_structure_on_peers(bound):
    db, molecule, structure, linked, _, _ = bound
    record = db.list_molecules(USER, molecule_ids=[molecule])["items"][0]
    pose = record["structures"][0]
    pose["context"].pop("target_path")
    resolved = bound_docking.receptor_file(db, USER, pose)
    assert resolved.read_bytes() == Path(linked["target_path"]).read_bytes()
    assert not resolved.stat().st_mode & 0o222
    with pytest.raises(PermissionError):
        bound_docking.prepare_bound_inputs(db, {"id": "b" * 32}, molecule, structure_id=structure)


def test_bulk_shortlist_upgrade_links_separate_target_without_duplicating_poses(bound, tmp_path):
    import pyarrow as pa
    import pyarrow.parquet as pq
    from mn_molecular_design_studio.bulk_screening_results import publish_final_selection
    db, _, _, _, receptor, pose = bound
    run = tmp_path / "runs/campaign/final"
    data = run / "artifacts/active_learning"
    (data / "retained_poses").mkdir(parents=True)
    retained = data / "retained_poses/pose.pdbqt"
    retained.write_bytes(pose.read_bytes())
    (run / "molecular_design_studio.json").write_text(json.dumps({"studio_owner_id": USER["id"]}))
    (run / "graph_snapshot.json").write_text(json.dumps({"nodes": [{"data": {"kind": "active_target",
        "params": {"target_pdb": str(receptor)}}}]}))
    path = data / "docking_scores.parquet"
    pq.write_table(pa.Table.from_pylist([{"compound_id": "X", "compound_key": 1, "library_id": "L", "score": -9,
        "smiles": "CCO", "pose_path": str(retained)}]), path)
    original = path.read_bytes()
    result = publish_final_selection(run, USER, 1, database=db, explicit_selection=True)
    row = db.list_molecules(USER, group_id=result["group_id"])["items"][0]
    assert row["type"] == "small_molecule" and len(row["components"]) == 1
    assert row["structures"][0]["context"]["target_molecule_id"] != row["id"]
    structure_id = row["structures"][0]["id"]
    # Emulate the previous shortlist representation: same portable ligand,
    # without the target relation or the viewer's best-pose policy.
    with db.connection() as c:
        context = dict(row["structures"][0]["context"])
        for key in ("target_structure_id", "target_molecule_id", "target_path", "database_structure_policy"):
            context.pop(key)
        c.execute("UPDATE molecule_structures SET context=? WHERE id=?", (json.dumps(context), structure_id))
    publish_final_selection(run, USER, 1, database=db, explicit_selection=True)
    row = db.list_molecules(USER, group_id=result["group_id"])["items"][0]
    assert len(row["structures"]) == 1 and row["structures"][0]["id"] == structure_id
    assert bound_docking.is_bound_pose(row["structures"][0])
    assert row["observations"][0]["context"]["representation_id"] == structure_id
    assert path.read_bytes() == original


def test_viewer_composes_receptor_and_ligand_without_a_complex_identity(bound, monkeypatch):
    db, molecule, _, _, _, _ = bound
    from mn_molecular_design_studio.engine.workflow_engine.app.components import molstar_viewer
    observed = []
    monkeypatch.setattr(molstar_viewer, "molstar_custom_component", lambda structures, **kw: observed.append(structures))
    script = f'''
from mn_molecular_design_studio.database_ui import render_details
from mn_molecular_design_studio.molecular_database import MolecularDatabase
db = MolecularDatabase()
user = {USER!r}
render_details(db, user, db.list_molecules(user,molecule_ids=[{molecule!r}])['items'][0])
'''
    app = AppTest.from_string(script).run()
    assert not app.exception
    assert len(observed[0]) == 2
    assert observed[0][0].representation_type == "cartoon"
    assert observed[0][1].representation_type == "ball-and-stick"
    assert "REMARK SMILES" not in observed[0][0].pdb
    assert any("separate Database molecules" in text.value for text in app.caption)


@pytest.mark.parametrize("engine", ["openmm", "gromacs"])
def test_md_tiles_use_database_pdb_and_sdf_without_uploads(bound, monkeypatch, engine):
    db, molecule, structure, _, _, _ = bound
    from mn_molecular_design_studio.engine.workflow_engine.workflows import ligand_tools
    monkeypatch.setattr(ligand_tools, engine + "_md_runner_readiness", lambda: {"ready": True, "source": "fixture"})
    prepared = bound_docking.prepare_bound_inputs(db, USER, molecule, structure_id=structure)
    workflow = "ligand.md." + engine
    app = AppTest.from_string(f'''
import streamlit as st
from mn_molecular_design_studio.ligand_tools_ui import render_setup
st.session_state["database-pending-md-selection"] = {dict(prepared, workflow=workflow, project_id="project")!r}
def submit(schema, workflow, parameters, project, settings, title):
    st.session_state['submitted'] = parameters
render_setup({workflow!r}, {{}}, owner_id={USER['id']!r}, compound_sets=[], project_selector=lambda key:'project', submit=submit)
''').run()
    assert not app.exception
    assert app.selectbox(key="ligand-gromacs-bound-residue").value == prepared["ligand_key"]
    assert not app.button(key="ligand-tool-run-" + workflow).disabled
    app.button(key="ligand-tool-run-" + workflow).click().run()
    assert not app.exception
    assert app.session_state["submitted"]["complex_pdb"] == prepared["complex_pdb"]
    assert app.session_state["submitted"]["ligand_sdf"] == prepared["ligand_sdf"]
