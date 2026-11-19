from __future__ import annotations

import hashlib
import json

import pytest

from mn_molecular_design_studio.data_browser import folder_scopes, saved_results
from mn_molecular_design_studio.database_browser import build_browser_index, subject_results
from mn_molecular_design_studio.database_browser_labels import tool_label
from mn_molecular_design_studio.legacy_target_links import plan_target_links, recover_target, register_target
from mn_molecular_design_studio.library_curation import apply_curation
from test_database_browser import USER, _import, _record, database


ADMIN = {**USER, "role": "admin"}


def _receptor(sequence):
    residues = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE", "G": "GLY", "H": "HIS", "I": "ILE"}
    return ("\n".join(f"ATOM  {number:5d}  CA  {residues[letter]} A{number:4d}    "
        f"{number * 3.8:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00 20.00           C"
        for number, letter in enumerate(sequence, 1)) + "\nEND\n").encode()


def _attempt(root, tool, pdb="4LNW", sequence="ACDE", *, attempt_id="original-attempt", compound="compound-1"):
    task = "docking" if tool in {"vina", "gnina", "udp"} else "refolding"
    run = root / task / ("mnligand-" + tool.replace(" ", "-") + "-" + pdb)
    attempt = run / "artifacts/replicates" / attempt_id
    attempt.mkdir(parents=True)
    (run / "metadata.json").write_text(json.dumps({"pdb_id": pdb, "tool": tool}))
    (attempt / "input.json").write_text(json.dumps({"target": {"kind": "run_artifact", "role": "axis_aligned_receptor",
        "sha256": "f" * 64, "run_id": "prepared-target"}}))
    context = {"task_group": task, "run_id": run.name, "legacy_run_id": attempt_id,
        "candidate_id": attempt_id + "::" + compound + "::legacy-replicate-001"}
    if task == "docking":
        receptor = attempt / "input/receptor.pdb"
        receptor.parent.mkdir()
        receptor.write_bytes(_receptor(sequence))
        context.update(target_id="Legacy receptor label", target_structure_digest=hashlib.sha256(receptor.read_bytes()).hexdigest())
    else:
        inputs = attempt / "inputs"
        inputs.mkdir()
        (inputs / (compound + ".json")).write_text(json.dumps({"sequences": [
            {"protein": {"id": "A", "sequence": sequence}}, {"ligand": {"id": "L", "smiles": "CCO"}}]}))
    return run, context


@pytest.mark.parametrize("tool", ["vina", "gnina", "udp", "Boltz-2", "Nesso-1", "AlphaFold 3"])
def test_exact_inputs_link_all_imported_methods_without_admitting_unreviewed_designs(database, tmp_path, tool):
    root = tmp_path / "imports"
    _import(database, name="Compound library", run_id="", method="user_uploads", kind="curated",
        records=[_record(name="Ethanol", components=[{"type": "small_molecule", "chain": "L", "smiles": "CCO"}])])
    expected = {}
    for pdb, sequence in [("4LNW", "ACDE"), ("3GWS", "FGHI")]:
        run, context = _attempt(root, tool, pdb, sequence)
        imported = _import(database, name="Original evaluation", run_id=run.name, task=run.parent.name, kind="evaluation",
            records=[_record(name="Ethanol", components=[{"type": "small_molecule", "chain": "L", "smiles": "CCO"}])])
        observation = database.add_observation(imported["group_id"], imported["molecule_ids"][0], tool, run.name,
            {"score": -7.1}, context=context, source_key=pdb)
        expected[sequence] = observation
    with database.connection() as connection:
        before = [tuple(row) for row in connection.execute("SELECT id,molecule_id,namespace,run_id,metrics FROM molecule_observations ORDER BY id")]
        members = [tuple(row) for row in connection.execute("SELECT * FROM molecule_members ORDER BY group_id,molecule_id")]
    plan = plan_target_links(database, ADMIN, root)
    assert not plan["unresolved"]
    assert len(plan["observations"]) == 2
    assert len(plan["targets"]) == 2
    apply_curation(database, ADMIN, plan, backup_path=tmp_path / "recovery.sqlite3")
    index = build_browser_index(database, ADMIN, collection="library", include_related=True)
    targets = [row for row in folder_scopes(index, []).values() if row["parent"] == "targets"]
    assert len(targets) == 2
    for target in targets:
        raw = subject_results(database, index, target["subject"])
        assert {row["observation_id"] for row in raw if row["observation_id"]} == {expected[target["subject"]["sequence"]]}
        outcomes = saved_results(index, raw)
        assert len(outcomes) == 1
        assert outcomes[0]["activity"] == ("Docking" if tool in {"vina", "gnina", "udp"} else
            "Binding affinity" if tool == "Nesso-1" else "Cofolding")
        assert outcomes[0]["eligibility"] == "Imported"
        assert outcomes[0]["target_id"]
    with database.connection() as connection:
        assert before == [tuple(row) for row in connection.execute("SELECT id,molecule_id,namespace,run_id,metrics FROM molecule_observations ORDER BY id")]
        assert all(row in [tuple(entry) for entry in connection.execute("SELECT * FROM molecule_members")] for row in members)
    repeated = plan_target_links(database, ADMIN, root)
    assert not repeated["observations"] and not repeated["unresolved"]


def test_missing_input_never_borrows_another_compound_or_uses_a_matching_target_name(tmp_path):
    run, context = _attempt(tmp_path, "Boltz-2")
    wrong = {**context, "candidate_id": "original-attempt::not-this-compound::legacy-replicate-001"}
    with pytest.raises(ValueError, match="missing or ambiguous"):
        recover_target(run, wrong, namespace="Boltz-2")
    with pytest.raises(ValueError, match="conflicting"):
        recover_target(run, {**context, "legacy_run_id": "different-attempt"}, namespace="Boltz-2")


def test_actual_docking_copy_wins_over_reused_artifact_metadata_and_rejects_wrong_digest(tmp_path):
    run, context = _attempt(tmp_path, "vina")
    binding = recover_target(run, context, namespace="vina")
    assert binding["components"][0]["sequence"] == "ACDE"
    assert binding["evidence"]["input_artifact_sha256"] == "f" * 64
    assert binding["evidence"]["snapshots"][1]["sha256"] == context["target_structure_digest"]
    with pytest.raises(ValueError, match="conflicts with"):
        recover_target(run, {**context, "target_structure_digest": "0" * 64}, namespace="vina")


def test_reconciliation_rejects_changed_evidence_and_keeps_original_observation(database, tmp_path):
    run, context = _attempt(tmp_path / "imports", "Nesso-1")
    imported = _import(database, run_id=run.name, task="refolding", kind="evaluation")
    identifier = database.add_observation(imported["group_id"], imported["molecule_ids"][0], "Nesso-1", run.name, {"score": 1}, context=context)
    plan = plan_target_links(database, ADMIN, tmp_path / "imports")
    model = run / "artifacts/replicates/original-attempt/inputs/compound-1.json"
    model.write_text(model.read_text().replace("ACDE", "FGHI"))
    with pytest.raises(ValueError, match="original target input changed"):
        apply_curation(database, ADMIN, plan, backup_path=tmp_path / "backup.sqlite3")
    with database.connection() as connection:
        assert not json.loads(connection.execute("SELECT context FROM molecule_observations WHERE id=?", (identifier,)).fetchone()[0]).get("target_molecule_id")


def test_existing_conflicting_target_is_not_silently_replaced(database, tmp_path):
    run, context = _attempt(tmp_path / "imports", "Boltz-2")
    wrong = _import(database, records=[_record("FGHI")], kind="reference")
    imported = _import(database, run_id=run.name, task="refolding", kind="evaluation")
    database.add_observation(imported["group_id"], imported["molecule_ids"][0], "Boltz-2", run.name, {"score": 1},
        context={**context, "target_molecule_id": wrong["molecule_ids"][0]})
    plan = plan_target_links(database, ADMIN, tmp_path / "imports")
    assert not plan["observations"]
    assert "conflicts" in plan["unresolved"][0]["reason"]


def test_future_registration_is_a_target_input_not_a_passed_design(database, tmp_path):
    run, context = _attempt(tmp_path, "Nesso-1")
    binding = recover_target(run, context, namespace="Nesso-1")
    target = register_target(database, ADMIN, "", binding, source_key="legacy-run:target", origin={"run_id": run.name})
    with database.connection() as connection:
        member = connection.execute("SELECT mm.metadata,g.origin FROM molecule_members mm JOIN molecule_groups g ON g.id=mm.group_id WHERE molecule_id=?", (target,)).fetchone()
        assert json.loads(member["metadata"])["target_identity_type"] == "molecular_identity"
        assert json.loads(member["origin"])["collection_kind"] == "inputs"
        assert connection.execute("SELECT count(*) FROM molecule_observations").fetchone()[0] == 0
    assert tool_label("udp") == tool_label("unidock_pro") == "Uni-Dock Pro"


def test_repaired_target_links_survive_standard_peer_replication(database, tmp_path):
    from mn_molecular_design_studio.cluster.replication import Replica
    from mn_molecular_design_studio.molecular_database import MolecularDatabase
    from test_molecular_database_peer import _sync_all
    root = tmp_path / "imports"
    run, context = _attempt(root, "Nesso-1", pdb="3GWS")
    imported = _import(database, run_id=run.name, task="refolding", kind="curated", method="user_uploads",
        records=[_record(name="Ethanol", components=[{"type": "small_molecule", "chain": "L", "smiles": "CCO"}])])
    database.add_observation(imported["group_id"], imported["molecule_ids"][0], "Nesso-1", run.name,
        {"binder_probability": 0.8}, context=context)
    plan = plan_target_links(database, ADMIN, root)
    apply_curation(database, ADMIN, plan, backup_path=tmp_path / "recovery.sqlite3")
    peer = MolecularDatabase(tmp_path / "peer")
    _sync_all(Replica(), database, peer)
    index = build_browser_index(peer, ADMIN, collection="library", include_related=True)
    target = next(row for row in folder_scopes(index, []).values() if row["parent"] == "targets")
    outcomes = saved_results(index, subject_results(peer, index, target["subject"]))
    assert len(outcomes) == 1
    assert outcomes[0]["activity"] == "Binding affinity"
    assert outcomes[0]["eligibility"] == "Imported"
    from mn_molecular_design_studio.database_design_evaluations import connected_evaluations
    original = outcomes[0]["primary"]
    connected, _ = connected_evaluations(index, original, [original])
    assert connected[0]["structure_relation"] == "Selected calculation (no coordinates)"
    assert "3GWS" in outcomes[0]["target_name"]


def test_new_legacy_affinity_ingestion_records_exact_target_without_filter_pass(database, tmp_path):
    from mn_molecular_design_studio.ligand_compounds import write_compound_set
    from mn_molecular_design_studio.molecular_database_ingest import ingest_run
    from mn_molecular_design_studio.workspace import STUDIO_SOURCE_ID
    run, context = _attempt(tmp_path / "imports", "Nesso-1")
    metadata = json.loads((run / "metadata.json").read_text())
    metadata.update(status="completed", workflow="nesso_affinity", database_collection_policy="automatic")
    (run / "metadata.json").write_text(json.dumps(metadata))
    (run / "input.json").write_text(json.dumps({"tool": "Nesso-1"}))
    (run / "result.json").write_text(json.dumps({"success": True}))
    (run / "mn_ligand_aggregate.json").write_text("{}")
    (run / "molecular_design_studio.json").write_text(json.dumps({"application": STUDIO_SOURCE_ID,
        "studio_owner_id": USER["id"], "project_id": ""}))
    _import(database, name="Library", run_id="", kind="curated", method="user_uploads",
        records=[_record(name="Ethanol", components=[{"type": "small_molecule", "chain": "L", "smiles": "CCO"}])])
    write_compound_set(run / "artifacts/input/compound_set.json", name="Consumed",
        compounds=[{"compound_id": "compound-1", "smiles": "CCO"}], owner_id=USER["id"], project_id="", created_at="now")
    (run / "artifacts/nesso_affinity.csv").write_text(
        "compound_id,legacy_run_id,binder_probability\ncompound-1,original-attempt,0.8\n")
    assert ingest_run(run, database=database)
    assert ingest_run(run, database=database) == 0
    index = build_browser_index(database, ADMIN, collection="library", include_related=True)
    target = next(row for row in folder_scopes(index, []).values() if row["parent"] == "targets")
    outcomes = saved_results(index, subject_results(database, index, target["subject"]))
    assert len(outcomes) == 1
    assert outcomes[0]["activity"] == "Binding affinity"
    assert outcomes[0]["eligibility"] == "Imported"


def test_same_construct_and_compound_do_not_merge_distinct_pdb_studies(database, tmp_path):
    root = tmp_path / "imports"
    imported = _import(database, run_id="", method="user_uploads", kind="curated",
        records=[_record(name="Ethanol", components=[{"type": "small_molecule", "chain": "L", "smiles": "CCO"}])])
    for pdb in ["3GWS", "6KKB"]:
        run, context = _attempt(root, "Nesso-1", pdb=pdb)
        database.add_observation(imported["group_id"], imported["molecule_ids"][0], "Nesso-1", run.name,
            {"binder_probability": 0.8}, context=context)
    plan = plan_target_links(database, ADMIN, root)
    apply_curation(database, ADMIN, plan, backup_path=tmp_path / "recovery.sqlite3")
    index = build_browser_index(database, ADMIN, collection="library", include_related=True)
    roots = [row for row in folder_scopes(index, []).values() if row["parent"] == "targets"]
    assert {row["label"] for row in roots} == {"TRβ · 3GWS", "TRβ · 6KKB"}
    for scope in roots:
        outcomes = saved_results(index, subject_results(database, index, scope["subject"]))
        assert len(outcomes) == 1
        assert outcomes[0]["target_pdb_id"] == scope["subject"]["target_pdb_id"]
        assert scope["subject"]["target_pdb_id"] in outcomes[0]["target_reference_label"]
    all_subject = {"id": "all", "molecule_ids": set(index["molecules"]), "group_ids": set(index["groups"])}
    all_results = saved_results(index, subject_results(database, index, all_subject))
    assert len([row for row in all_results if row["activity"] == "Binding affinity"]) == 2


def test_pdb_family_exposes_distinct_constructs_without_merging_them(database, tmp_path):
    root = tmp_path / "imports"
    imported = _import(database, run_id="", method="user_uploads", kind="curated",
        records=[_record(name="Ethanol", components=[{"type": "small_molecule", "chain": "L", "smiles": "CCO"}])])
    run, context = _attempt(root, "Nesso-1", sequence="ACDE", attempt_id="first")
    _attempt(root, "Nesso-1", sequence="FGHIA", attempt_id="second")
    for attempt in ["first", "second"]:
        calculation = {**context, "legacy_run_id": attempt,
            "candidate_id": attempt + "::compound-1::legacy-replicate-001"}
        database.add_observation(imported["group_id"], imported["molecule_ids"][0], "Nesso-1", run.name,
            {"binder_probability": 0.8}, context=calculation)
    plan = plan_target_links(database, ADMIN, root)
    apply_curation(database, ADMIN, plan, backup_path=tmp_path / "recovery.sqlite3")
    index = build_browser_index(database, ADMIN, collection="library", include_related=True)
    scopes = folder_scopes(index, [])
    family = scopes["targets:pdb:4LNW"]
    constructs = [scope for scope in scopes.values() if scope["parent"] == family["id"]]
    assert {scope["label"] for scope in constructs} == {"Protein construct · 4 aa", "Protein construct · 5 aa"}
    assert family["entity"] is None
    assert len(saved_results(index, subject_results(database, index, family["subject"]))) == 2
    for scope in constructs:
        assert len(saved_results(index, subject_results(database, index, scope["subject"]))) == 1


def test_native_calculation_without_pdb_stays_accessible_without_guessing_its_reference(database, tmp_path):
    root = tmp_path / "imports"
    run, context = _attempt(root, "Nesso-1")
    imported = _import(database, run_id="", method="user_uploads", kind="curated",
        records=[_record(name="Ethanol", components=[{"type": "small_molecule", "chain": "L", "smiles": "CCO"}])])
    database.add_observation(imported["group_id"], imported["molecule_ids"][0], "Nesso-1", run.name, {"score": 1}, context=context)
    plan = plan_target_links(database, ADMIN, root)
    apply_curation(database, ADMIN, plan, backup_path=tmp_path / "recovery.sqlite3")
    target_id = plan["targets"][0]["id"]
    native = database.add_observation(imported["group_id"], imported["molecule_ids"][0], "native-affinity", "native-run", {"score": 2},
        context={"target_molecule_id": target_id, "workflow_role": "binding_affinity", "saved_selection": True})
    index = build_browser_index(database, ADMIN, collection="library", include_related=True)
    scopes = folder_scopes(index, [])
    unspecified = next(row for row in scopes.values() if row["subject"].get("without_pdb_reference"))
    assert "PDB unspecified" in unspecified["label"]
    outcomes = saved_results(index, subject_results(database, index, unspecified["subject"]))
    assert len(outcomes) == 1 and outcomes[0]["observation_id"] == native
    assert not any(row["observation_id"] == native for row in subject_results(database, index, scopes["targets:pdb:4LNW"]["subject"]))


def test_input_only_screening_coordinates_do_not_admit_unreviewed_compounds(database, tmp_path):
    target = _import(database, run_id="", task="", kind="reference", method="user_uploads",
        records=[_record(name="Receptor", components=[{"type": "protein", "chain": "A", "sequence": "ACDE", "role": "target"}],
            metadata={"target_identity_type": "molecular_identity"})])
    inputs = _import(database, run_id="native-screen", task="ligand-docking", kind="inputs", method="tool_inputs",
        records=[_record(name="Unreviewed compound", components=[{"type": "small_molecule", "chain": "L", "smiles": "CCO"}])])
    path = tmp_path / "pose.pdbqt"
    path.write_text("ROOT\nENDROOT\n")
    structure = database.attach_structure(inputs["group_id"], inputs["molecule_ids"][0], {"path": str(path), "label": "Stored screening pose"},
        context={"run_id": "native-screen", "task_group": "ligand-docking"})
    database.add_observation(inputs["group_id"], inputs["molecule_ids"][0], "unidock_pro", "native-screen", {"score": -7},
        context={"target_molecule_id": target["molecule_ids"][0], "representation_id": structure,
            "task_group": "ligand-docking", "workflow_role": "docking"})
    index = build_browser_index(database, ADMIN, collection="library", include_related=True)
    scope = next(row for row in folder_scopes(index, []).values() if row["parent"] == "targets")
    raw = subject_results(database, index, scope["subject"])
    assert any(row["tool"] == "unidock_pro" for row in raw)
    assert not any("unidock_pro" in row["origin_tools"] for row in saved_results(index, raw))
