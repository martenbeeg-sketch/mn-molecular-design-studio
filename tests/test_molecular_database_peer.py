from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from mn_molecular_design_studio.molecular_database import MolecularDatabase, _id, _json
from mn_molecular_design_studio.molecular_database_peer import apply_record, capture, materialize_portable_structures
from mn_molecular_design_studio.cluster.replication import Replica
from mn_molecular_design_studio import molecule_imports


SOURCE = "source-node-a"
OTHER_SOURCE = "source-node-b"
OWNER = {"id": "a" * 32, "role": "user", "_cluster_project_ids": ["project-current", "project-next"]}
OUTSIDER = {"id": "b" * 32, "role": "user", "_cluster_project_ids": []}


@pytest.fixture
def roots(tmp_path, monkeypatch):
    # Both workspaces and the replica metadata store must live under pytest's
    # temporary root; avoid reading or writing a developer's live Studio DB.
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_HOME", str(tmp_path / "studio"))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(tmp_path / "runtime"))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR", str(tmp_path / "runs"))
    monkeypatch.setenv("MDS_CLUSTER_HOME", str(tmp_path / "cluster"))
    return tmp_path


def _capture_rows(replica, kind=None):
    with replica.db() as connection:
        query = "SELECT key,kind,value FROM records"
        rows = connection.execute(query).fetchall()
    return {row["key"]: (row["kind"], json.loads(row["value"])) for row in rows
            if kind is None or row["kind"] == kind}


def _sync_all(replica, source_database, target_database, author=SOURCE):
    capture(replica, database=source_database, node_id=author)
    records = _capture_rows(replica)
    for key, (kind, value) in sorted(records.items(), key=lambda item: (item[1][0] != "molecule_group", item[0])):
        if key.startswith(kind + ":" + author + ":"):
            apply_record(kind, key, value, author, database=target_database)


def test_saved_comparison_keeps_canonical_evaluations_across_sync(roots):
    source, target = MolecularDatabase(roots / "source-evals"), MolecularDatabase(roots / "target-evals")
    imported = source.import_records(OWNER, molecule_imports.parse_sequences("ACDE"))
    molecule, group = imported["molecule_ids"][0], imported["group_id"]
    source.add_observation(group, molecule, "test", "first", {"score": 1}, context={"target_id": "t"}, created_at="2026-01-01")
    source.add_observation(group, molecule, "test", "second", {"score": 2}, context={"target_id": "t"}, created_at="2026-02-01")
    saved = source.save_comparison(OWNER, "Frozen", [molecule], evaluation_filters={"run_id": "first"})
    _sync_all(Replica(), source, target)
    replica_group = _id("peer-", SOURCE + ":" + saved["group_id"])
    restored = target.list_molecules(OWNER, group_id=replica_group)["items"]
    assert len(restored) == 1
    assert restored[0]["scores"]["test.score"] == 1
    assert restored[0]["selected_evaluation"]["id"] == source.list_molecules(OWNER, group_id=saved["group_id"])["items"][0]["selected_evaluation"]["id"]


def test_synchronized_generator_observations_project_same_common_values(roots):
    source = MolecularDatabase(roots / "generator-source")
    target = MolecularDatabase(roots / "generator-target")
    imported = source.import_records(OWNER, molecule_imports.parse_sequences("ACDE"), name="Generated designs")
    source.add_observation(imported["group_id"], imported["molecule_ids"][0], "bindcraft2_generator", "generation-1",
        {"pTM": .8, "Target_pLDDT": .9}, context={"candidate": {"source_tool": "bindcraft2_generator", "stage": "generation.backbone_sequence", "binder_sequence": "ACDE"}})
    replica = Replica()
    _sync_all(replica, source, target)
    row = target.list_molecules(OWNER)["items"][0]
    assert row["scores"]["common.confidence.generator.complex.ptm"] == .8
    assert row["scores"]["common.confidence.generator.target.plddt"] == 90
    assert row["observations"][0]["metrics"]["Target_pLDDT"] == .9
    assert "common.design.sequence_length" in target.metric_fields(OWNER)


def test_composed_comparison_preserves_subject_and_sources_across_peer_sync(roots):
    from mn_molecular_design_studio.comparisons import compare
    source, target = MolecularDatabase(roots / "source-composition"), MolecularDatabase(roots / "target-composition")
    imported = source.import_records(OWNER, molecule_imports.parse_sequences("ACDE"))
    molecule, group = imported["molecule_ids"][0], imported["group_id"]
    context = {"target_id": "t", "structure_digest": "shared-digest", "target_structure_digest": "target-digest"}
    source.add_observation(group, molecule, "prediction", "p", {"score": .9}, context=context)
    source.add_observation(group, molecule, "quality", "q", {"solubility": .5}, context=context)
    settings = {"combine_same_subject": True, "fields": ["prediction.score", "quality.solubility"],
                "directions": {"prediction.score": "higher", "quality.solubility": "higher"}}
    saved = source.save_comparison(OWNER, "Combined", [molecule], settings=settings)
    expected = compare(source.evaluation_rows(OWNER, group_id=saved["group_id"]), **settings)
    _sync_all(Replica(), source, target)
    replica_group = _id("peer-", SOURCE + ":" + saved["group_id"])
    actual = compare(target.evaluation_rows(OWNER, group_id=replica_group), **settings)
    assert len(actual) == 1
    assert actual.iloc[0]["Source evaluations"] == expected.iloc[0]["Source evaluations"]
    assert actual.iloc[0]["Evaluation ID"] == expected.iloc[0]["Evaluation ID"]
    assert actual.iloc[0]["Comparison rank"] == 1


def test_portable_structure_is_copied_to_peer_database_and_survives_source_offline(roots, monkeypatch):
    source = MolecularDatabase(roots / "portable-source")
    target = MolecularDatabase(roots / "portable-target")
    content = b"HEADER    portable protein structure\nEND\n"
    imported = source.import_records(
        OWNER,
        [{"name": "portable", "components": [{"chain": "A", "type": "protein", "sequence": "ACDE"}],
          "files": [{"name": "model.pdb", "content": content}]}],
        name="Portable structures", project_id="project-current",
    )
    molecule_id = imported["molecule_ids"][0]
    with source.connection() as connection:
        source_structure_id = connection.execute(
            "SELECT id FROM molecule_structures WHERE group_id=? AND molecule_id=?",
            (imported["group_id"], molecule_id),
        ).fetchone()[0]
    portable = source.make_structure_portable(OWNER, source_structure_id)
    assert Path(portable["path"]).is_file()
    assert portable["sha256"] == hashlib.sha256(content).hexdigest()

    replica = Replica()
    _sync_all(replica, source, target)
    peer_group_id = _id("peer-", SOURCE + ":" + imported["group_id"])
    with target.connection() as connection:
        target_structure = dict(connection.execute(
            "SELECT * FROM molecule_structures WHERE group_id=? AND molecule_id=?",
            (peer_group_id, molecule_id),
        ).fetchone())
    assert target_structure["path"].startswith("cluster-library://")

    from mn_molecular_design_studio.cluster import transport
    monkeypatch.setattr(transport, "node_request", lambda *args, **kwargs: content)
    assert materialize_portable_structures(database=target, peer={"id": SOURCE, "url": "https://source"}) == 1
    with target.connection() as connection:
        target_structure = dict(connection.execute(
            "SELECT * FROM molecule_structures WHERE id=?", (target_structure["id"],)
        ).fetchone())
    assert not target_structure["path"].startswith("cluster-library://")
    assert target.structure_file(OWNER, target_structure["id"]).read_bytes() == content

    # A later metadata refresh must not replace the local portable copy with a
    # remote URI; that copy remains usable if the source workstation is offline.
    with source.connection() as connection:
        connection.execute("UPDATE molecule_groups SET name=? WHERE id=?", ("Renamed portable structures", imported["group_id"]))
    _sync_all(replica, source, target)
    with target.connection() as connection:
        refreshed = dict(connection.execute(
            "SELECT * FROM molecule_structures WHERE id=?", (target_structure["id"],)
        ).fetchone())
    assert not refreshed["path"].startswith("cluster-library://")
    assert target.structure_file(OWNER, refreshed["id"]).read_bytes() == content


def _make_source_with_entry(root):
    database = MolecularDatabase(root)
    imported = database.import_records(
        OWNER,
        [{"name": "ethanol", "components": molecule_imports.parse_smiles("CCO")[0]["components"],
          "files": [{"name": "model.pdb", "content": b"HEADER    source model\nEND\n"},],
          "metadata": {"tag": "original"}}],
        name="Reusable molecules", project_id="project-current", import_method="tool_output",
        origin={"type": "curated_collection", "curated": True, "collection_kind": "curated"},
    )
    group_id, molecule_id = imported["group_id"], imported["molecule_ids"][0]
    observation_id = database.add_observation(
        group_id, molecule_id, "vina", "run-1", {"score": -7.2},
        context={"settings": {"exhaustiveness": 8}}, source_key="score-source-1",
    )
    return database, group_id, molecule_id, observation_id


def _group_record(database, group_id):
    with database.connection() as connection:
        return dict(connection.execute("SELECT * FROM molecule_groups WHERE id=?", (group_id,)).fetchone())


def _entry_record(database, group_id, molecule_id, *, digest=True):
    with database.connection() as connection:
        molecule = dict(connection.execute("SELECT * FROM molecules WHERE id=?", (molecule_id,)).fetchone())
        member = dict(connection.execute("SELECT * FROM molecule_members WHERE group_id=? AND molecule_id=?", (group_id, molecule_id)).fetchone())
        structures = [dict(row) for row in connection.execute(
            "SELECT * FROM molecule_structures WHERE group_id=? AND molecule_id=? ORDER BY id", (group_id, molecule_id))]
        observations = [dict(row) for row in connection.execute(
            "SELECT * FROM molecule_observations WHERE group_id=? AND molecule_id=? ORDER BY id", (group_id, molecule_id))]
        sources = [dict(row) for row in connection.execute(
            "SELECT * FROM molecule_sources WHERE group_id=? AND molecule_id=? ORDER BY source_key", (group_id, molecule_id))]
    for structure in structures:
        path = structure["path"]
        data = Path(path).read_bytes() if path else b""
        structure["path"] = ""
        structure["sha256"] = hashlib.sha256(data).hexdigest() if digest and data else ""
        context = json.loads(structure["context"])
        target_path = context.pop("target_path", None)
        if target_path:
            with database.connection() as connection:
                target = connection.execute("SELECT id FROM molecule_structures WHERE path=? LIMIT 1", (target_path,)).fetchone()
            if target:
                context["target_structure_id"] = target[0]
        structure["context"] = _json(context)
    return {"node_id": SOURCE, "group_id": group_id, "molecule": molecule, "member": member,
            "structures": structures, "observations": observations, "sources": sources}


def _peer_structure(node):
    residues = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU"}
    lines = [f"HEADER    authored by {node}"]
    for serial, amino_acid in enumerate("ACDE", start=1):
        lines.append(f"ATOM  {serial:5d}  CA  {residues[amino_acid]} A{serial:4d}    "
                     f"{serial:8.3f}{1.0:8.3f}{1.0:8.3f}{1.00:6.2f}{20.00:6.2f}           C")
    return ("\n".join([*lines, "TER", "END", ""])).encode()


def _author_peer_entry(database, node, *, project="project-current"):
    content = _peer_structure(node)
    imported = database.import_records(
        OWNER,
        [{"name": f"ACDE from {node}", "components": molecule_imports.parse_sequences("ACDE")[0]["components"],
          "files": [{"name": f"{node}.pdb", "content": content}],
          "metadata": {"source_workstation": node}}],
        name=f"Local library {node}", project_id=project, import_method="user_uploads",
        origin={"type": "curated_collection", "collection_kind": "curated", "owner_id": OWNER["id"]},
    )
    molecule_id, group_id = imported["molecule_ids"][0], imported["group_id"]
    with database.connection() as connection:
        structure_id = connection.execute(
            "SELECT id FROM molecule_structures WHERE group_id=? AND molecule_id=?", (group_id, molecule_id)
        ).fetchone()[0]
    database.make_structure_portable(OWNER, structure_id)
    database.add_observation(group_id, molecule_id, f"evaluation.{node}", f"run-{node}",
                             {"quality": 0.7 + len(node) / 100},
                             context={"source_workstation": node}, source_key=f"{node}-evaluation")
    return group_id, molecule_id, structure_id, content


def _sync_author_records(replica, author, targets):
    rows = _capture_rows(replica)
    records = [(key, kind, value) for key, (kind, value) in rows.items()
               if key.startswith(kind + ":" + author + ":")]
    records.sort(key=lambda item: (item[1] != "molecule_group", item[0]))
    for key, kind, value in records:
        for target in targets:
            apply_record(kind, key, value, author, database=target)


def _evaluation_run(root, node, source, digest):
    from mn_molecular_design_studio.workspace import STUDIO_SOURCE_ID

    run_dir = root / "runs" / "analysis" / f"remote-evaluation-{node}"
    candidate_dir = run_dir / "artifacts/normalized_candidates"
    candidate_dir.mkdir(parents=True)
    (run_dir / "molecular_design_studio.json").write_text(json.dumps({
        "application": STUDIO_SOURCE_ID, "studio_owner_id": OWNER["id"], "project_id": "project-current",
    }))
    (run_dir / "metadata.json").write_text(json.dumps({
        "studio_owner_id": OWNER["id"], "project_id": "project-current", "task_group": "analysis",
        "tool": "complex_analysis", "workflow": "analysis.complex_analysis", "status": "completed",
        "database_collection_policy": "automatic",
    }))
    (run_dir / "input.json").write_text(json.dumps({"tool": "complex_analysis", "workflow": "analysis.complex_analysis", "params": {}}))
    (run_dir / "result.json").write_text(json.dumps({"success": True, "metrics": {}}))
    candidate = {
        "candidate_id": f"analysis-{node}", "stage": "analysis", "binder_sequence": "ACDE",
        "complex_pdb": "analysis.pdb", "metrics": {"quality": 0.95},
        "raw_metadata": {"database_source": source,
                         "complex_analysis_selection": {"source_structure_digest": digest}},
    }
    (run_dir / "analysis.pdb").write_bytes(_peer_structure("analysis"))
    (candidate_dir / "candidates.jsonl").write_text(json.dumps(candidate) + "\n")
    return run_dir


def _results_source_directory(root, node, run_id, task_group="design"):
    source = root / "runs" / task_group / run_id
    source.mkdir(parents=True)
    (source / "metadata.json").write_text(json.dumps({"run_id": run_id, "task_group": task_group}))
    (source / "molecular_design_studio.json").write_text(json.dumps({"workstation": node}))
    return source


def _results_analysis_run(root, node, source_dir, digest):
    from mn_molecular_design_studio.workspace import STUDIO_SOURCE_ID

    run = root / "runs" / "analysis" / f"results-analysis-{node}"
    candidates = run / "artifacts/normalized_candidates"
    candidates.mkdir(parents=True)
    (run / "molecular_design_studio.json").write_text(json.dumps({
        "application": STUDIO_SOURCE_ID, "studio_owner_id": OWNER["id"], "project_id": "project-current",
    }))
    (run / "metadata.json").write_text(json.dumps({
        "studio_owner_id": OWNER["id"], "project_id": "project-current", "task_group": "analysis",
        "tool": "plip", "workflow": "analysis.plip", "status": "completed",
        "database_collection_policy": "automatic",
    }))
    (run / "input.json").write_text(json.dumps({"tool": "plip", "workflow": "analysis.plip",
                                                 "source_run_dir": str(source_dir)}))
    (run / "result.json").write_text(json.dumps({"success": True, "metrics": {}}))
    candidate = {
        "candidate_id": f"plip-{node}", "stage": "analysis", "source_tool": "plip",
        "binder_sequence": "ACDE", "binder_chains": ["A"], "complex_pdb": "analysis.pdb",
        "metrics": {"interaction_count": 4},
        "raw_metadata": {"source_candidate": {"candidate_id": "source-candidate", "binder_sequence": "ACDE",
                                                   "binder_chains": ["A"]},
                         "complex_analysis_selection": {"source_structure_digest": digest}},
    }
    (run / "analysis.pdb").write_bytes(_peer_structure("analysis"))
    (candidates / "candidates.jsonl").write_text(json.dumps(candidate) + "\n")
    return run


def _records_for_author(replica, author):
    return [(key, kind, value) for key, (kind, value) in _capture_rows(replica).items()
            if key.startswith(kind + ":" + author + ":")]


def _apply_author_records(replica, author, target):
    for key, kind, value in sorted(_records_for_author(replica, author),
                                   key=lambda item: (item[1] != "molecule_group", item[0])):
        apply_record(kind, key, value, author, database=target)


def test_small_molecule_replication_preserves_source_smiles_across_rdkit_versions(roots, monkeypatch):
    from mn_molecular_design_studio import molecular_database_peer

    source = MolecularDatabase(roots / "source")
    components = molecule_imports.parse_smiles("C[C@H](O)C(=O)N")[0]["components"]
    imported = source.import_records(OWNER, [{"name": "chiral input", "components": components}], name="Inputs")
    group_id, molecule_id = imported["group_id"], imported["molecule_ids"][0]
    group_value = {"node_id": SOURCE, "group": _group_record(source, group_id)}
    entry_value = _entry_record(source, group_id, molecule_id)
    target = MolecularDatabase(roots / "target")
    apply_record("molecule_group", f"molecule_group:{SOURCE}:{group_id}", group_value, SOURCE, database=target)

    # Receiving peers must not use their local RDKit canonicalizer to rewrite a
    # source's stored SMILES. That can change the string and stereochemistry.
    monkeypatch.setattr(
        molecular_database_peer, "_canonical_components",
        lambda *_args, **_kwargs: pytest.fail("peer replication recanonicalized a source molecule"),
    )
    apply_record(
        "molecule_entry", f"molecule_entry:{SOURCE}:{group_id}:{molecule_id}",
        entry_value, SOURCE, database=target,
    )
    with target.connection() as connection:
        stored = connection.execute("SELECT fingerprint,components FROM molecules WHERE id=?", (molecule_id,)).fetchone()
        assert stored["fingerprint"] == entry_value["molecule"]["fingerprint"]
        assert stored["components"] == entry_value["molecule"]["components"]


def test_small_molecule_replication_rejects_fingerprint_mismatch_atomically(roots):
    source = MolecularDatabase(roots / "source")
    components = molecule_imports.parse_smiles("C[C@H](O)C(=O)N")[0]["components"]
    imported = source.import_records(OWNER, [{"name": "chiral input", "components": components}], name="Inputs")
    group_id, molecule_id = imported["group_id"], imported["molecule_ids"][0]
    target = MolecularDatabase(roots / "target")
    apply_record(
        "molecule_group", f"molecule_group:{SOURCE}:{group_id}",
        {"node_id": SOURCE, "group": _group_record(source, group_id)}, SOURCE, database=target,
    )
    malformed = _entry_record(source, group_id, molecule_id)
    changed = json.loads(malformed["molecule"]["components"])
    changed[0]["smiles"] = "C[C@@H](O)C(=O)N"
    malformed["molecule"]["components"] = _json(changed)
    malformed["member"]["components"] = _json(changed)

    with pytest.raises(ValueError, match="fingerprint"):
        apply_record(
            "molecule_entry", f"molecule_entry:{SOURCE}:{group_id}:{molecule_id}",
            malformed, SOURCE, database=target,
        )
    with target.connection() as connection:
        assert connection.execute("SELECT count(*) FROM molecules WHERE id=?", (molecule_id,)).fetchone()[0] == 0
        assert connection.execute("SELECT count(*) FROM molecule_members WHERE molecule_id=?", (molecule_id,)).fetchone()[0] == 0


def test_group_snapshot_round_trips_schema_curated_flags_and_project_acl(roots):
    source = MolecularDatabase(roots / "source")
    replica = MolecularDatabase(roots / "replica")
    group_id = source.create_group(
        OWNER, "Curated collection", project_id="project-current", import_method="tool_output",
        origin={"type": "curated_collection", "curated": True, "collection_kind": "curated", "published": True},
        schema={"fields": [{"name": "confidence", "type": "number", "required": True}]},
    )
    transport = Replica()
    capture(transport, database=source, node_id=SOURCE)
    key = f"molecule_group:{SOURCE}:{group_id}"
    kind, value = _capture_rows(transport)[key]
    apply_record(kind, key, value, SOURCE, database=replica)

    peer_group_id = _id("peer-", SOURCE + ":" + group_id)
    group = replica.group(peer_group_id, OWNER)
    assert group["project_id"] == "project-current"
    assert group["import_method"] == "tool_output"
    assert group["schema"] == {"fields": [{"name": "confidence", "type": "number", "required": True}]}
    assert group["origin"]["curated"] is True
    assert group["origin"]["collection_kind"] == "curated"
    assert replica.groups(OUTSIDER, scope="organization") == []

    # Moving the source collection updates the replica's ACL fields and kind.
    with source.connection() as connection:
        connection.execute(
            "UPDATE molecule_groups SET project_id=?,import_method=?,updated_at=? WHERE id=?",
            ("project-next", "user_uploads", "2026-10-20T12:00:00Z", group_id),
        )
    capture(transport, database=source, node_id=SOURCE)
    kind, value = _capture_rows(transport)[key]
    apply_record(kind, key, value, SOURCE, database=replica)
    updated = replica.group(peer_group_id, OWNER)
    assert updated["project_id"] == "project-next"
    assert updated["import_method"] == "user_uploads"
    assert updated["origin"]["curated"] is True
    assert replica.groups(OUTSIDER, scope="organization") == []

    with pytest.raises(PermissionError):
        apply_record(kind, key, {**value, "node_id": OTHER_SOURCE}, OTHER_SOURCE, database=replica)
    assert replica.group(peer_group_id, OWNER)["project_id"] == "project-next"


def test_entry_update_reconciles_deleted_relations_and_keeps_other_records(roots):
    source, group_id, molecule_id, observation_id = _make_source_with_entry(roots / "source")
    replica = MolecularDatabase(roots / "replica")
    transport = Replica()
    _sync_all(transport, source, replica)
    peer_group_id = _id("peer-", SOURCE + ":" + group_id)

    with replica.connection() as connection:
        # Preserve a local relation in the same pair and a second peer group's
        # membership/score for the globally shared molecule identity.
        connection.execute("INSERT INTO molecule_observations VALUES (?,?,?,?,?,?,?,?)",
                           ("local-score", molecule_id, peer_group_id, "local", "local-run", "{}", "{}", "now"))
        other_group_id = "other-source-group"
        other_peer_group_id = _id("peer-", OTHER_SOURCE + ":" + other_group_id)
        connection.execute("INSERT INTO molecule_groups VALUES (?,?,?,?,?,?,?,?,?)",
                           (other_peer_group_id, "Other source", OUTSIDER["id"], "", "tool_output",
                            _json({"_replica_node": OTHER_SOURCE, "_replica_group": other_group_id}), "{}", "now", "now"))
        connection.execute("INSERT INTO molecule_members VALUES (?,?,?,?,?,?)",
                           (other_peer_group_id, molecule_id, "Other peer's member", "[]", "{}", "now"))
        connection.execute("INSERT INTO molecule_observations VALUES (?,?,?,?,?,?,?,?)",
                           ("other-peer-score", molecule_id, other_peer_group_id, "other", "other-run", "{}",
                            _json({"_replica_node": OTHER_SOURCE, "source_node": OTHER_SOURCE}), "now"))

    # An update to existing relation IDs changes their values without leaving
    # stale duplicates. This also updates the member's source-authoritative name.
    with source.connection() as connection:
        connection.execute("UPDATE molecule_members SET name=?,metadata=? WHERE group_id=? AND molecule_id=?",
                           ("ethanol revised", _json({"tag": "revised"}), group_id, molecule_id))
        connection.execute("UPDATE molecule_structures SET label=? WHERE group_id=? AND molecule_id=?",
                           ("revised model", group_id, molecule_id))
        connection.execute("UPDATE molecule_observations SET metrics=?,context=? WHERE id=?",
                           (_json({"score": -8.1}), _json({"settings": {"exhaustiveness": 16}}), observation_id))
        connection.execute("UPDATE molecule_sources SET provenance=? WHERE group_id=? AND molecule_id=?",
                           (_json({"curated": True}), group_id, molecule_id))
    capture(transport, database=source, node_id=SOURCE)
    entry_key = f"molecule_entry:{SOURCE}:{group_id}:{molecule_id}"
    kind, value = _capture_rows(transport)[entry_key]
    apply_record(kind, entry_key, value, SOURCE, database=replica)
    with replica.connection() as connection:
        member = connection.execute("SELECT name,metadata FROM molecule_members WHERE group_id=? AND molecule_id=?", (peer_group_id, molecule_id)).fetchone()
        structure = connection.execute("SELECT label FROM molecule_structures WHERE group_id=? AND molecule_id=?", (peer_group_id, molecule_id)).fetchone()
        score = connection.execute("SELECT metrics FROM molecule_observations WHERE group_id=? AND molecule_id=? AND namespace='vina'", (peer_group_id, molecule_id)).fetchone()
        provenance = connection.execute("SELECT provenance FROM molecule_sources WHERE group_id=? AND molecule_id=?", (peer_group_id, molecule_id)).fetchone()
    assert member["name"] == "ethanol revised"
    assert json.loads(member["metadata"]) == {"tag": "revised"}
    assert structure["label"] == "revised model"
    assert json.loads(score["metrics"]) == {"score": -8.1}
    assert json.loads(provenance["provenance"])["curated"] is True

    # Source now removes every copied relation. A local relation in the same
    # group/molecule and the other peer group's rows must survive reconciliation.
    with source.connection() as connection:
        for table in ("molecule_sources", "molecule_structures", "molecule_observations"):
            connection.execute(f"DELETE FROM {table} WHERE group_id=? AND molecule_id=?", (group_id, molecule_id))
    capture(transport, database=source, node_id=SOURCE)
    kind, value = _capture_rows(transport)[entry_key]
    apply_record(kind, entry_key, value, SOURCE, database=replica)
    with replica.connection() as connection:
        assert connection.execute("SELECT count(*) FROM molecule_structures WHERE group_id=? AND molecule_id=?", (peer_group_id, molecule_id)).fetchone()[0] == 0
        assert connection.execute("SELECT count(*) FROM molecule_sources WHERE group_id=? AND molecule_id=?", (peer_group_id, molecule_id)).fetchone()[0] == 0
        assert connection.execute("SELECT count(*) FROM molecule_observations WHERE id='local-score'").fetchone()[0] == 1
        assert connection.execute("SELECT count(*) FROM molecule_observations WHERE id='other-peer-score'").fetchone()[0] == 1
        assert connection.execute("SELECT count(*) FROM molecule_members WHERE group_id=? AND molecule_id=?", (other_peer_group_id, molecule_id)).fetchone()[0] == 1
        assert connection.execute("SELECT count(*) FROM molecule_members WHERE group_id=? AND molecule_id=?", (peer_group_id, molecule_id)).fetchone()[0] == 1


def test_entry_tombstone_and_malformed_snapshot_are_atomic_and_author_scoped(roots):
    source, group_id, molecule_id, _ = _make_source_with_entry(roots / "source")
    replica = MolecularDatabase(roots / "replica")
    transport = Replica()
    _sync_all(transport, source, replica)
    peer_group_id = _id("peer-", SOURCE + ":" + group_id)
    entry_key = f"molecule_entry:{SOURCE}:{group_id}:{molecule_id}"
    snapshot = _entry_record(source, group_id, molecule_id)
    with replica.connection() as connection:
        # Local records in a replicated pair survive a source tombstone.
        connection.execute(
            "INSERT INTO molecule_observations VALUES (?,?,?,?,?,?,?,?)",
            ("local-score", molecule_id, peer_group_id, "local", "local-run", "{}", "{}", "now"),
        )

    malformed = json.loads(json.dumps(snapshot))
    malformed["member"]["name"] = "should roll back"
    malformed["structures"][0]["context"] = {"not_finite": float("nan")}
    before = replica.list_molecules(OWNER, group_id=peer_group_id)["items"][0]
    with pytest.raises(ValueError):
        apply_record("molecule_entry", entry_key, malformed, SOURCE, database=replica)
    after = replica.list_molecules(OWNER, group_id=peer_group_id)["items"][0]
    assert after["name"] == before["name"]

    apply_record("molecule_entry", entry_key, None, SOURCE, database=replica)
    assert replica.list_molecules(OWNER, group_id=peer_group_id)["total"] == 0
    with replica.connection() as connection:
        for table in ("molecule_members", "molecule_structures", "molecule_sources"):
            assert connection.execute(f"SELECT count(*) FROM {table} WHERE group_id=? AND molecule_id=?", (peer_group_id, molecule_id)).fetchone()[0] == 0
        assert connection.execute("SELECT count(*) FROM molecule_observations WHERE id='local-score'").fetchone()[0] == 1
        assert connection.execute(
            "SELECT count(*) FROM molecule_observations WHERE group_id=? AND molecule_id=? AND id!='local-score'",
            (peer_group_id, molecule_id),
        ).fetchone()[0] == 0

    # A source's tombstone cannot delete a group's row with a different author.
    other_group = "other-source-group"
    other_peer_group = _id("peer-", OTHER_SOURCE + ":" + other_group)
    group_value = {"node_id": OTHER_SOURCE, "group": {
        "id": other_group, "name": "Other", "owner_id": OUTSIDER["id"], "project_id": "", "import_method": "user_uploads",
        "origin": "{}", "schema_json": "{}", "created_at": "now", "updated_at": "now"}}
    apply_record("molecule_group", f"molecule_group:{OTHER_SOURCE}:{other_group}", group_value, OTHER_SOURCE, database=replica)
    with pytest.raises(PermissionError):
        apply_record("molecule_entry", entry_key, None, OTHER_SOURCE, database=replica)
    assert replica.group(other_peer_group, OUTSIDER)["name"] == "Other"


def test_same_identity_with_group_specific_chains_roles_and_order_replicates(roots):
    source = MolecularDatabase(roots / "source")
    target = MolecularDatabase(roots / "target")
    first = source.import_records(OWNER, [{
        "name": "first representation",
        "components": [
            {"chain": "A", "type": "protein", "sequence": "ACDE", "role": "binder"},
            {"chain": "B", "type": "protein", "sequence": "FGHI", "role": "target"},
        ],
    }], name="Source one")
    second = source.import_records(OWNER, [{
        "name": "second representation",
        "components": [
            {"chain": "H", "type": "protein", "sequence": "FGHI", "role": "target"},
            {"chain": "L", "type": "protein", "sequence": "ACDE", "role": "binder"},
        ],
    }], name="Source two")
    molecule_id = first["molecule_ids"][0]
    assert second["molecule_ids"] == [molecule_id]

    # A destination's existing copy may use its own chain names and order.
    target.import_records(OUTSIDER, [{
        "name": "local representation",
        "components": [
            {"chain": "X", "type": "protein", "sequence": "ACDE", "role": "binder"},
            {"chain": "Y", "type": "protein", "sequence": "FGHI", "role": "target"},
        ],
    }], name="Local copy")
    transport = Replica()
    _sync_all(transport, source, target)

    first_peer_group = _id("peer-", SOURCE + ":" + first["group_id"])
    second_peer_group = _id("peer-", SOURCE + ":" + second["group_id"])
    first_row = target.list_molecules(OWNER, group_id=first_peer_group)["items"][0]
    second_row = target.list_molecules(OWNER, group_id=second_peer_group)["items"][0]
    assert first_row["id"] == second_row["id"] == molecule_id
    assert first_row["components"][0]["chain"] == "A"
    assert second_row["components"][0]["chain"] == "H"
    assert [component["sequence"] for component in second_row["components"]] == ["FGHI", "ACDE"]
    assert target.list_molecules(OUTSIDER, group_id=target.groups(OUTSIDER)[0]["id"])["total"] == 1


def test_group_tombstone_removes_only_that_peer_group(roots):
    source, group_id, molecule_id, _ = _make_source_with_entry(roots / "source")
    replica = MolecularDatabase(roots / "replica")
    transport = Replica()
    _sync_all(transport, source, replica)
    peer_group_id = _id("peer-", SOURCE + ":" + group_id)
    with replica.connection() as connection:
        local_group = "local-group"
        connection.execute("INSERT INTO molecule_groups VALUES (?,?,?,?,?,?,?,?,?)",
                           (local_group, "Local", OWNER["id"], "", "user_uploads", "{}", "{}", "now", "now"))
        connection.execute("INSERT INTO molecule_members VALUES (?,?,?,?,?,?)",
                           (local_group, molecule_id, "Local copy", "[]", "{}", "now"))
    apply_record("molecule_group", f"molecule_group:{SOURCE}:{group_id}", None, SOURCE, database=replica)
    with replica.connection() as connection:
        assert connection.execute("SELECT count(*) FROM molecule_groups WHERE id=?", (peer_group_id,)).fetchone()[0] == 0
        assert connection.execute("SELECT count(*) FROM molecule_members WHERE group_id=? AND molecule_id=?", (local_group, molecule_id)).fetchone()[0] == 1
        assert connection.execute("SELECT count(*) FROM molecules WHERE id=?", (molecule_id,)).fetchone()[0] == 1


def test_three_workstations_ingest_offline_converge_retrieve_evaluate_and_tombstone(roots, monkeypatch):
    from mn_molecular_design_studio.cluster import config, transport
    from mn_molecular_design_studio.molecular_database_ingest import ingest_run

    nodes = ("pc002975", "PC004032", "PC004033")
    databases = {node: MolecularDatabase(roots / f"database-{node}") for node in nodes}
    authored = {node: _author_peer_entry(databases[node], node) for node in nodes}
    replica = Replica()

    # PC004033 has authored its own local record but is offline during the
    # first exchange between pc002975 and PC004032.
    for node in nodes[:2]:
        capture(replica, database=databases[node], node_id=node)
        _sync_author_records(replica, node, [databases[target] for target in nodes[:2] if target != node])
    assert databases[nodes[0]].list_molecules(OWNER, collection="all")["total"] == 1
    with databases[nodes[0]].connection() as connection:
        assert connection.execute("SELECT count(*) FROM molecule_groups").fetchone()[0] == 2

    # The offline writer retains its local record and Evaluation. Reconnection
    # exchanges all three source-author snapshots; a second replay is idempotent.
    assert databases[nodes[2]].list_molecules(OWNER, collection="all")["total"] == 1
    for node in nodes:
        capture(replica, database=databases[node], node_id=node)
    for node in nodes:
        _sync_author_records(replica, node, [databases[target] for target in nodes if target != node])
    for node in nodes:
        _sync_author_records(replica, node, [databases[target] for target in nodes if target != node])

    canonical_ids = {authored[node][1] for node in nodes}
    assert len(canonical_ids) == 1
    molecule_id = next(iter(canonical_ids))
    for node in nodes:
        database = databases[node]
        with database.connection() as connection:
            assert connection.execute("SELECT count(*) FROM molecules WHERE id=?", (molecule_id,)).fetchone()[0] == 1
            assert connection.execute("SELECT count(*) FROM molecule_groups").fetchone()[0] == 3
            assert connection.execute("SELECT count(*) FROM molecule_structures WHERE molecule_id=?", (molecule_id,)).fetchone()[0] == 3
            assert connection.execute("SELECT count(*) FROM molecule_observations WHERE molecule_id=?", (molecule_id,)).fetchone()[0] == 3
            assert connection.execute("SELECT count(*) FROM molecule_sources WHERE molecule_id=?", (molecule_id,)).fetchone()[0] == 3

    # pc002975 authored the representation; PC004032 can retrieve its bytes
    # through owner-authorized content retrieval and attach a new local
    # Evaluation to the same canonical identity and remote structure ID.
    source_node = nodes[0]
    owner_group, _, source_structure, source_content = authored[source_node]
    target_database = databases[nodes[1]]
    remote_group = _id("peer-", source_node + ":" + owner_group)
    with target_database.connection() as connection:
        remote_row = dict(connection.execute(
            "SELECT id,context FROM molecule_structures WHERE group_id=? AND molecule_id=?",
            (remote_group, molecule_id),
        ).fetchone())
    remote_context = json.loads(remote_row["context"])
    remote_structure_id = remote_row["id"]
    remote_digest = remote_context.get("portable_sha256") or remote_context.get("sha256")
    assert remote_context["_replica_node"] == source_node
    assert remote_context["_replica_structure"] == source_structure
    source_spec = {"schema": "studio.complex_analysis_database_source.v1", "molecule_id": molecule_id,
                   "structure_id": remote_structure_id, "group_id": remote_group,
                   "structure_digest": remote_digest}
    peer_urls = {node: f"https://{node}.isolated.test" for node in nodes}
    monkeypatch.setattr(config, "load", lambda: {"mode": "disabled", "peers": [
        {"id": node, "url": peer_urls[node]} for node in nodes]})

    def fetch_from_owner(path, payload, *, method, base_url, raw, timeout):
        node = next(node for node, url in peer_urls.items() if url == base_url)
        requested = path.rsplit("/", 1)[-1]
        source_db = databases[node]
        return source_db.structure_file(OWNER, requested).read_bytes()

    monkeypatch.setattr(transport, "node_request", fetch_from_owner)
    assert target_database.structure_file(OWNER, remote_structure_id).read_bytes() == source_content
    with pytest.raises(PermissionError):
        target_database.structure_file(OUTSIDER, remote_structure_id)

    # The simulated remote evaluation output is written only to PC004032's
    # local Evaluation group. Its source group remains source-authoritative.
    run_dir = _evaluation_run(roots / "analysis-PC004032", nodes[1], source_spec, remote_digest)
    assert ingest_run(run_dir, database=target_database) == 1
    evaluation_groups = [group for group in target_database.groups(OWNER, collection="all")
                         if group["origin"].get("collection_kind") == "evaluation"]
    assert len(evaluation_groups) == 1
    evaluated = target_database.list_molecules(OWNER, group_id=evaluation_groups[0]["id"])["items"]
    assert len(evaluated) == 1 and evaluated[0]["id"] == molecule_id
    evaluation = next(observation for observation in evaluated[0]["observations"]
                      if observation["run_id"] == run_dir.name)
    assert evaluation["context"]["source_group_id"] == remote_group
    assert evaluation["context"]["source_structure_id"] == remote_structure_id
    assert evaluation["context"]["source_structure_digest"] == remote_digest
    assert evaluation["context"]["source_node"] == source_node
    with target_database.connection() as connection:
        assert connection.execute("SELECT count(*) FROM molecule_observations WHERE group_id=? AND molecule_id=?",
                                  (remote_group, molecule_id)).fetchone()[0] == 1
        assert connection.execute("SELECT count(*) FROM molecule_observations WHERE group_id=? AND molecule_id=?",
                                  (evaluation_groups[0]["id"], molecule_id)).fetchone()[0] == 1

    # PC004032's Evaluation is authored in its own collection but names
    # pc002975's structure. A third workstation and the structure owner must
    # resolve the same owner-qualified reference when that Evaluation is
    # replayed onward.
    capture(replica, database=target_database, node_id=nodes[1])
    _sync_author_records(replica, nodes[1], [databases[nodes[0]], databases[nodes[2]]])
    expected_owner_group = owner_group
    expected_owner_structure = source_structure
    for receiver_node in (nodes[0], nodes[2]):
        receiver = databases[receiver_node]
        replicated_evaluation_group = _id("peer-", nodes[1] + ":" + evaluation_groups[0]["id"])
        with receiver.connection() as connection:
            row = connection.execute(
                "SELECT context FROM molecule_observations WHERE group_id=? AND run_id=?",
                (replicated_evaluation_group, run_dir.name),
            ).fetchone()
            assert row is not None
            context = json.loads(row["context"])
            assert context["source_structure_owner_node"] == source_node
            assert context["source_structure_owner_group_id"] == expected_owner_group
            assert context["source_structure_owner_structure_id"] == expected_owner_structure
            assert context["source_structure_reference_status"] == "resolved"
            expected_group = (expected_owner_group if receiver_node == source_node else
                              _id("peer-", source_node + ":" + expected_owner_group))
            expected_structure = (expected_owner_structure if receiver_node == source_node else
                                  _id("peer-structure-", source_node + ":" + expected_owner_structure))
            assert context["source_group_id"] == expected_group
            assert context["source_structure_id"] == expected_structure
            assert context["representation_id"] == expected_structure
            assert connection.execute(
                "SELECT count(*) FROM molecule_structures WHERE id=? AND group_id=? AND molecule_id=?",
                (expected_structure, expected_group, molecule_id),
            ).fetchone()[0] == 1

    # Materialize an authorized portable copy, then prove it remains openable
    # when the source workstation is unavailable.
    assert materialize_portable_structures(database=target_database,
        peer={"id": source_node, "url": peer_urls[source_node]}) >= 1
    monkeypatch.setattr(transport, "node_request",
        lambda *args, **kwargs: (_ for _ in ()).throw(transport.ClusterUnavailable("offline")))
    assert target_database.structure_file(OWNER, remote_structure_id).read_bytes() == source_content

    # A source-authorized deletion publishes a tombstone for only PC004032's
    # own imported member; the peer Evaluation and the other authors survive.
    own_group, own_molecule, _, _ = authored[nodes[1]]
    databases[nodes[1]].remove_from_group(OWNER, own_group, [own_molecule])
    capture(replica, database=databases[nodes[1]], node_id=nodes[1])
    tombstone_key = f"molecule_entry:{nodes[1]}:{own_group}:{own_molecule}"
    captured = _capture_rows(replica)
    assert captured[tombstone_key] == ("molecule_entry", None)
    _sync_author_records(replica, nodes[1], [databases[nodes[0]], databases[nodes[2]]])
    for node in (nodes[0], nodes[2]):
        database = databases[node]
        removed_group = _id("peer-", nodes[1] + ":" + own_group)
        with database.connection() as connection:
            assert connection.execute("SELECT count(*) FROM molecule_members WHERE group_id=? AND molecule_id=?",
                                      (removed_group, molecule_id)).fetchone()[0] == 0
            assert connection.execute("SELECT count(*) FROM molecule_members WHERE molecule_id=?", (molecule_id,)).fetchone()[0] >= 2
            original_group = (authored[nodes[0]][0] if node == nodes[0]
                              else _id("peer-", nodes[0] + ":" + authored[nodes[0]][0]))
            assert connection.execute("SELECT count(*) FROM molecule_groups WHERE id=?",
                                      (original_group,)).fetchone()[0] == 1


def test_results_source_evaluation_preserves_owner_reference_across_three_peers_and_delivery_orders(roots):
    from mn_molecular_design_studio.molecular_database_ingest import ingest_run

    node_a, node_b = "results-owner-a", "results-evaluator-b"
    database_a = MolecularDatabase(roots / "results-source-a")
    database_b = MolecularDatabase(roots / "results-source-b")
    database_c_early = MolecularDatabase(roots / "results-source-c-evaluation-first")
    database_c_ordered = MolecularDatabase(roots / "results-source-c-source-first")
    content = _peer_structure("source")
    imported = database_a.import_records(OWNER, [{
        "name": "ACDE source", "components": molecule_imports.parse_sequences("ACDE")[0]["components"],
        "files": [{"name": "source.pdb", "content": content}],
    }], name="Results-backed source", project_id="project-current")
    molecule_id, group_a = imported["molecule_ids"][0], imported["group_id"]
    with database_a.connection() as connection:
        structure_a = connection.execute(
            "SELECT id FROM molecule_structures WHERE group_id=? AND molecule_id=?", (group_a, molecule_id)
        ).fetchone()[0]
        connection.execute("INSERT INTO molecule_sources VALUES (?,?,?,?)", (
            "results-source-entry", molecule_id, group_a,
            _json({"run_id": "source-results-run", "task_group": "design", "candidate_id": "source-candidate",
                   "workstation": node_a}),
        ))
    database_a.make_structure_portable(OWNER, structure_a)
    digest = __import__("hashlib").sha256(content).hexdigest()
    source_dir = _results_source_directory(roots / "results-run-files", node_a, "source-results-run")
    replica = Replica()
    capture(replica, database=database_a, node_id=node_a)
    _apply_author_records(replica, node_a, database_b)

    remote_group_a = _id("peer-", node_a + ":" + group_a)
    with database_b.connection() as connection:
        structure_b = connection.execute(
            "SELECT id FROM molecule_structures WHERE group_id=? AND molecule_id=?", (remote_group_a, molecule_id)
        ).fetchone()[0]
    analysis = _results_analysis_run(roots / "results-analysis-files", node_b, source_dir, digest)
    assert ingest_run(analysis, database=database_b, force=True) == 1
    evaluation_group_b = next(group["id"] for group in database_b.groups(OWNER, collection="all")
                              if group["origin"].get("collection_kind") == "evaluation")
    evaluation = next(item for item in database_b.list_molecules(OWNER, group_id=evaluation_group_b)["items"][0]["observations"]
                      if item["run_id"] == analysis.name)
    context_b = evaluation["context"]
    assert context_b["source_kind"] == "results"
    assert context_b["source_structure_owner_node"] == node_a
    assert context_b["source_structure_owner_group_id"] == group_a
    assert context_b["source_structure_owner_structure_id"] == structure_a
    assert context_b["source_node"] == node_a
    assert context_b["workstation"] != node_a  # Evaluator and source author remain distinct.

    capture(replica, database=database_b, node_id=node_b)
    # Ordered delivery resolves immediately; delayed delivery must converge
    # after A's source entry arrives without a new evaluation from B.
    _apply_author_records(replica, node_a, database_c_ordered)
    _apply_author_records(replica, node_b, database_c_ordered)
    _apply_author_records(replica, node_b, database_c_early)
    early_group_b = _id("peer-", node_b + ":" + evaluation_group_b)
    with database_c_early.connection() as connection:
        pending = connection.execute(
            "SELECT context FROM molecule_observations WHERE group_id=? AND run_id=?",
            (early_group_b, analysis.name),
        ).fetchone()
        assert json.loads(pending["context"])["source_structure_reference_status"] == "owner_record_unavailable"
    _apply_author_records(replica, node_a, database_c_early)

    expected_group_c = _id("peer-", node_a + ":" + group_a)
    expected_structure_c = _id("peer-structure-", node_a + ":" + structure_a)
    for receiver in (database_c_ordered, database_c_early):
        _apply_author_records(replica, node_a, receiver)  # Duplicate source replay is idempotent.
        _apply_author_records(replica, node_b, receiver)  # Duplicate evaluation replay is idempotent.
        receiver_group_b = _id("peer-", node_b + ":" + evaluation_group_b)
        with receiver.connection() as connection:
            row = connection.execute(
                "SELECT id,context FROM molecule_observations WHERE group_id=? AND run_id=?",
                (receiver_group_b, analysis.name),
            ).fetchone()
            context = json.loads(row["context"])
            assert context["source_structure_reference_status"] == "resolved"
            assert context["source_group_id"] == expected_group_c
            assert context["source_structure_id"] == expected_structure_c
            assert context["representation_id"] == expected_structure_c
            assert connection.execute(
                "SELECT count(*) FROM molecule_structures WHERE id=? AND group_id=? AND molecule_id=?",
                (expected_structure_c, expected_group_c, molecule_id),
            ).fetchone()[0] == 1

    # A -> B -> C -> A: evaluation authorship stays with B while the exact
    # source relationship resolves back to A's original representation.
    _apply_author_records(replica, node_b, database_a)
    with database_a.connection() as connection:
        row = connection.execute(
            "SELECT context FROM molecule_observations WHERE run_id=? AND group_id!=?",
            (analysis.name, group_a),
        ).fetchone()
        returned = json.loads(row["context"])
        assert returned["source_structure_reference_status"] == "resolved"
        assert returned["source_group_id"] == group_a
        assert returned["source_structure_id"] == structure_a
        assert returned["source_structure_owner_node"] == node_a
        assert returned["_replica_node"] == node_b


def test_database_source_evaluation_reconciles_when_owner_records_arrive_late(roots):
    source = MolecularDatabase(roots / "late-source-a")
    evaluator = MolecularDatabase(roots / "late-source-b")
    late = MolecularDatabase(roots / "late-source-c-evaluation-first")
    ordered = MolecularDatabase(roots / "late-source-c-source-first")
    owner_group, molecule_id, owner_structure, source_content = _author_peer_entry(source, SOURCE)
    source_digest = hashlib.sha256(source_content).hexdigest()
    replica = Replica()
    capture(replica, database=source, node_id=SOURCE)
    _apply_author_records(replica, SOURCE, evaluator)
    local_source_group = _id("peer-", SOURCE + ":" + owner_group)
    local_source_structure = _id("peer-structure-", SOURCE + ":" + owner_structure)
    evaluation_group = evaluator.create_group(OWNER, "Local delayed-source Evaluation",
        project_id="project-current", import_method="tool_output",
        origin={"type": "job", "collection_kind": "evaluation", "workstation": "evaluator-node-b"})
    evaluator.add_to_group(OWNER, evaluation_group, [molecule_id])
    evaluator.add_observation(evaluation_group, molecule_id, "plip", "delayed-source-analysis", {"interactions": 4},
        context={"source_molecule_id": molecule_id, "source_group_id": local_source_group,
                 "source_structure_id": local_source_structure, "representation_id": local_source_structure,
                 "source_structure_digest": source_digest, "source_node": SOURCE,
                 "source_structure_owner_node": SOURCE,
                 "source_structure_owner_group_id": owner_group,
                 "source_structure_owner_structure_id": owner_structure},
        source_key="delayed-source-evaluation")
    capture(replica, database=evaluator, node_id=OTHER_SOURCE)

    _apply_author_records(replica, SOURCE, ordered)
    _apply_author_records(replica, OTHER_SOURCE, ordered)
    _apply_author_records(replica, OTHER_SOURCE, late)
    late_group = _id("peer-", OTHER_SOURCE + ":" + evaluation_group)
    with late.connection() as connection:
        pending = connection.execute(
            "SELECT context FROM molecule_observations WHERE group_id=? AND run_id='delayed-source-analysis'",
            (late_group,),
        ).fetchone()
        assert json.loads(pending["context"])["source_structure_reference_status"] == "owner_record_unavailable"
    _apply_author_records(replica, SOURCE, late)
    with late.connection() as connection:
        arrived = json.loads(connection.execute(
            "SELECT context FROM molecule_observations WHERE group_id=? AND run_id='delayed-source-analysis'",
            (late_group,),
        ).fetchone()[0])
    assert arrived["source_structure_reference_status"] == "resolved"
    assert arrived["source_structure_id"] == _id("peer-structure-", SOURCE + ":" + owner_structure)

    for receiver in (ordered, late):
        _apply_author_records(replica, SOURCE, receiver)
        _apply_author_records(replica, OTHER_SOURCE, receiver)
        with receiver.connection() as connection:
            context = json.loads(connection.execute(
                "SELECT context FROM molecule_observations WHERE group_id=? AND run_id='delayed-source-analysis'",
                (_id("peer-", OTHER_SOURCE + ":" + evaluation_group),),
            ).fetchone()[0])
        assert context["source_structure_reference_status"] == "resolved"
        assert context["source_structure_id"] == _id("peer-structure-", SOURCE + ":" + owner_structure)

    _apply_author_records(replica, OTHER_SOURCE, source)
    with source.connection() as connection:
        returned = json.loads(connection.execute(
            "SELECT context FROM molecule_observations WHERE run_id='delayed-source-analysis'"
        ).fetchone()[0])
    assert returned["_replica_node"] == OTHER_SOURCE
    assert returned["source_structure_reference_status"] == "resolved"
    assert returned["source_structure_id"] == owner_structure
