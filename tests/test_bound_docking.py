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


def test_receptor_resolves_via_authorized_database_structure_on_peers(bound):
    db, molecule, structure, linked, _, _ = bound
    record = db.list_molecules(USER, molecule_ids=[molecule])["items"][0]
    pose = record["structures"][0]
    pose["context"].pop("target_path")
    assert bound_docking.receptor_file(db, USER, pose) == Path(linked["target_path"])
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
    result = publish_final_selection(run, USER, 1, database=db)
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
    publish_final_selection(run, USER, 1, database=db)
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
