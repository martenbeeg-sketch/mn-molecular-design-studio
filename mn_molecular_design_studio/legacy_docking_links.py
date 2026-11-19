"""Plan and atomically repair docking links from verified original MN-Ligand runs."""
from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import os
from pathlib import Path
import sqlite3

from mn_molecular_design_studio.library_curation import _digest, encoded
from mn_molecular_design_studio.molecular_database import _id
from mn_molecular_design_studio.molecule_imports import parse_structure
from mn_molecular_design_studio.target_registration import contents_identity


POLICY = "studio.legacy-docking-links.v1"


def _engine(value):
    value = str(value or "").casefold().replace("-", "_")
    return "udp" if value in {"udp", "unidock_pro", "uni_dock_pro"} else value


def _inside(root, value):
    path = Path(value).resolve()
    if not path.is_relative_to(root.resolve()):
        raise ValueError("An original artifact escapes its selected source directory.")
    return path


def _snapshot(path):
    content = Path(path).read_bytes()
    return content, {"path": str(Path(path).resolve()), "sha256": hashlib.sha256(content).hexdigest()}


def _basename(value):
    value = str(value or "")
    if not value or Path(value).name != value or value in {".", ".."}:
        raise ValueError("An original run or compound selector is missing or invalid.")
    return value


def _number(value):
    try:
        number = float(value)
        return number if math.isfinite(number) else None
    except (ValueError, TypeError):
        return None


def _source_evidence(root, observation, cache):
    context = json.loads(observation["context"])
    row = context.get("row")
    if not isinstance(row, dict):
        raise ValueError("The original calculation row is unavailable.")
    attempt = _basename(row.get("legacy_run_id"))
    compound = _basename(row.get("compound_id"))
    if context.get("legacy_run_id") != attempt or not row.get("replicate"):
        raise ValueError("The original attempt or replicate selectors are missing or conflicting.")
    engine = _engine(observation["namespace"])
    if engine not in {"vina", "gnina", "udp"} or _engine(row.get("engine")) != engine:
        raise ValueError("The recorded docking engine is missing or conflicting.")
    directory = _inside(root, root / "docking" / attempt)
    if directory not in cache:
        metadata, metadata_snapshot = _snapshot(_inside(root, directory / "metadata.json"))
        inputs, input_snapshot = _snapshot(_inside(root, directory / "input.json"))
        scores, score_snapshot = _snapshot(_inside(root, directory / "scores.csv"))
        metadata = json.loads(metadata)
        if metadata.get("run_id") != attempt or metadata.get("status") != "completed":
            raise ValueError("The original source is not this completed attempt.")
        scores = scores.decode("utf-8-sig")
        reader = csv.DictReader(io.StringIO(scores), delimiter="\t" if "\t" in scores.partition("\n")[0] else ",")
        indexed = {}
        for source_row in reader:
            key = (source_row.get("compound_id"), source_row.get("replicate"))
            indexed.setdefault(key, []).append(source_row)
        cache[directory] = (metadata, json.loads(inputs), indexed, [metadata_snapshot, input_snapshot, score_snapshot])
    metadata, inputs, indexed, snapshots = cache[directory]
    if any(_engine(metadata[key]) != engine for key in ("tool", "engine") if metadata.get(key)):
        raise ValueError("The original attempt's engine conflicts with the SQL calculation.")
    matches = indexed.get((compound, str(row["replicate"])), [])
    if len(matches) != 1:
        raise ValueError("The exact original compound/replicate score row is missing or ambiguous.")
    source_row = matches[0]
    if _engine(source_row.get("engine")) != engine:
        raise ValueError("The original score row belongs to a different engine.")
    metrics = json.loads(observation["metrics"])
    score_fields = {key for key in ("best_score_kcal_mol", "cnn_score", "cnn_affinity", "score")
                    if _number(source_row.get(key)) is not None}
    if not score_fields:
        raise ValueError("The original row contains no verifiable docking score.")
    for key, value in source_row.items():
        numeric = _number(value)
        if numeric is not None:
            saved = _number(metrics.get(key))
            if saved is None or not math.isclose(saved, numeric, rel_tol=0, abs_tol=1e-10):
                raise ValueError("The original numeric score or replicate settings conflict with SQL: " + key)
        if key in row and value and row[key] != value and key not in {"pose_path", "pose_file"}:
            raise ValueError("The original calculation row conflicts with SQL: " + key)
    relative = source_row.get("pose_file") or source_row.get("pose_path")
    if not relative or Path(relative).is_absolute():
        raise ValueError("The original pose is not an explicit source-run-relative path.")
    expected = (Path("artifacts/replicates") / attempt / relative).as_posix()
    if row.get("pose_path") != expected:
        raise ValueError("The SQL pose path does not identify the exact original attempt artifact.")
    pose = _inside(directory, directory / relative)
    receptor = _inside(directory, directory / "input/receptor.pdb")
    pose_content, pose_snapshot = _snapshot(pose)
    receptor_content, receptor_snapshot = _snapshot(receptor)
    if not pose_content or not any(line.startswith((b"ATOM  ", b"HETATM")) for line in pose_content.splitlines()):
        raise ValueError("The original pose contains no atom coordinates.")
    reference = inputs.get("target_artifact") or inputs.get("prepared_target") or {}
    if not isinstance(reference, dict) or reference.get("role") not in {"receptor", "axis_aligned_receptor", "prepared_target"}:
        raise ValueError("The original input does not declare its prepared receptor.")
    for key in ("target_structure_digest", "target_variant_sha256"):
        if context.get(key) != receptor_snapshot["sha256"]:
            raise ValueError("The original receptor frame conflicts with the SQL target digest.")
    return context, source_row, pose_snapshot, receptor_snapshot, receptor_content, [*snapshots, pose_snapshot, receptor_snapshot]


def _structure(database, group, molecule, snapshot, context, label, stamp, identity):
    suffix = Path(snapshot["path"]).suffix.lower()
    if suffix not in {".pdb", ".pdbqt", ".sdf", ".mol", ".cif", ".mmcif"}:
        raise ValueError("The original coordinate format is unsupported.")
    digest = snapshot["sha256"]
    return {"id": _id("structure-", encoded([POLICY, group, molecule, identity, digest])),
        "molecule_id": molecule, "group_id": group,
        "path": str(database.root / "molecular_files" / digest[:2] / (digest + suffix)),
        "label": label, "format": suffix[1:], "created_at": stamp,
        "context": encoded({**context, "portable": True, "portable_sha256": digest, "content_sha256": digest}),
        "source": snapshot}


def plan_docking_links(database, user, source_runs_root, *, run_ids=None):
    """Read-only, fail-closed source reconciliation; never import extra molecules or scores."""
    if user.get("role") != "admin":
        raise PermissionError("Docking reconciliation requires an administrator.")
    root = Path(source_runs_root).resolve()
    with database.connection() as connection:
        connection.execute("BEGIN")
        groups = {row["id"]: dict(row) for row in connection.execute("SELECT * FROM molecule_groups")
                  if database._can_read(row, user) and not json.loads(row["origin"]).get("_replica_node")}
        observations = [dict(row) for row in connection.execute(
            "SELECT * FROM molecule_observations WHERE run_id LIKE 'mnligand-%' "
            "AND json_extract(context,'$.task_group')='docking' ORDER BY id")
            if row["group_id"] in groups and (run_ids is None or row["run_id"] in run_ids)]
        members = {(row["group_id"], row["molecule_id"]): dict(row) for row in connection.execute("SELECT * FROM molecule_members")}
        existing = {row["id"]: dict(row) for row in connection.execute("SELECT * FROM molecule_structures")}
    cache, parsed, structures, changes, unresolved, dependencies = {}, {}, {}, [], [], {}
    for observation in observations:
        try:
            context, row, pose, receptor, receptor_content, snapshots = _source_evidence(root, observation, cache)
            group = groups[observation["group_id"]]
            target = str(context.get("target_molecule_id") or "")
            candidates = [(key, member) for key, member in members.items() if key[1] == target and key[0] in groups
                and groups[key[0]]["owner_id"] == group["owner_id"] and groups[key[0]]["project_id"] == group["project_id"]]
            if not candidates:
                raise ValueError("The exact target has no source-owned membership in this project.")
            target_key, target_member = sorted(candidates)[0]
            digest = receptor["sha256"]
            if digest not in parsed:
                parsed[digest] = contents_identity(parse_structure(receptor_content, "receptor.pdb", include_ligands=False,
                    require_model_selection=True)[0]["components"])
            if parsed[digest] != contents_identity(json.loads(target_member["components"])):
                raise ValueError("The original receptor contents conflict with the registered target identity.")
            receptor_record = _structure(database, target_key[0], target, receptor,
                {"database_structure_policy": "docking_receptor", "dimensions": 3,
                 "target_molecule_id": target, "target_structure_digest": digest},
                "Docking receptor · verified original coordinate frame", observation["created_at"], "receptor")
            previous = {key: context[key] for key in ("representation_id", "structure_digest", "pose_digest",
                "calculation_structure_digest", "source_structure_digest", "target_structure_id", "target_path", "pose_id", "pose_identity") if key in context}
            prior_receipt = context.get("docking_link_reconciliation") or {}
            if prior_receipt.get("policy") == POLICY:
                previous = prior_receipt.get("previous", previous)
            receipt = {"policy": POLICY, "attempt_id": context["legacy_run_id"], "compound_id": row["compound_id"],
                "replicate": row["replicate"], "engine": _engine(observation["namespace"]), "snapshots": snapshots}
            base_context = {key: value for key, value in context.items()
                            if key not in {"representation_id", "docking_link_reconciliation"}}
            pose_context = {**base_context, "structure_digest": pose["sha256"], "pose_digest": pose["sha256"],
                "target_structure_id": receptor_record["id"], "target_path": receptor_record["path"],
                "database_structure_policy": "best_docking_pose", "docking_link_reconciliation": receipt}
            for key in ("source_structure_digest", "calculation_structure_digest"):
                if key in pose_context:
                    raise ValueError("This calculation contains an additional coordinate identity requiring manual review.")
            if context.get("pose_id") or context.get("pose_identity"):
                raise ValueError("This calculation has an explicit pose identity requiring manual review.")
            structure_fields = {"type", "run_id", "task_group", "tool", "workstation", "source_node",
                "legacy_run_id", "attempt_index", "candidate_id", "model_id", "pose_id", "pose_rank",
                "replicate", "compound_id", "structure_digest", "pose_digest", "target_molecule_id",
                "target_structure_id", "target_structure_digest", "target_variant_sha256", "target_path",
                "target_pdb_id", "target_receptor", "target_reference", "database_structure_policy"}
            structure_context = {key: value for key, value in pose_context.items() if key in structure_fields}
            structure_context["row"] = {key: value for key, value in context["row"].items()
                if value not in (None, "") and key in {"compound_id", "legacy_run_id", "replicate", "engine",
                    "pose_path", "pose_file", "pose_rank", "pose_index", "seed", "model_id", "attempt_index"}}
            structure_context["docking_link_reconciliation"] = {
                key: value for key, value in receipt.items() if key != "snapshots"}
            structure_context["docking_link_reconciliation"]["evidence_sha256"] = _digest(receipt)
            pose_record = _structure(database, group["id"], observation["molecule_id"], pose, structure_context,
                f"{observation['namespace']} · attempt {context['legacy_run_id'][:8]} · replicate {row['replicate']}",
                observation["created_at"], [observation["run_id"], context["legacy_run_id"], row["replicate"], observation["namespace"]])
            updated = {**pose_context, "representation_id": pose_record["id"],
                "docking_link_reconciliation": {**receipt, "previous": previous}}
            for record in (pose_record, receptor_record):
                prior = existing.get(record["id"])
                if prior and any(prior[key] != record[key] for key in ("molecule_id", "group_id", "path", "format", "context")):
                    raise ValueError("A verified structure identity already has different saved evidence.")
            if (prior_receipt.get("policy") == POLICY and context.get("representation_id") == pose_record["id"]
                    and context.get("target_structure_id") == receptor_record["id"]
                    and all(record["id"] in existing and _snapshot(record["path"])[1]["sha256"] == record["source"]["sha256"]
                            for record in (pose_record, receptor_record))):
                continue
            for record in (pose_record, receptor_record):
                structures[record["id"]] = record
            changes.append({"before": observation, "after_context": encoded(updated)})
            for key in ((group["id"], observation["molecule_id"]), target_key):
                dependencies[key] = {"group": groups[key[0]], "member": members[key]}
        except (OSError, ValueError, TypeError, KeyError) as exc:
            unresolved.append({"observation_id": observation["id"], "run_id": observation["run_id"], "reason": str(exc)})
    plan = {"schema": POLICY, "database": str(database.path.resolve()), "source_runs_root": str(root),
        "observations": changes, "structures": list(structures.values()), "dependencies": list(dependencies.values()),
        "unresolved": unresolved, "examined": len(observations)}
    plan["digest"] = _digest(plan)
    return plan


def apply_docking_links(database, user, plan, *, backup_path):
    """Apply source-owned links only, with a recovery backup and stale-plan checks."""
    if (user.get("role") != "admin" or plan.get("schema") != POLICY
            or plan.get("digest") != _digest({key: value for key, value in plan.items() if key != "digest"})
            or plan.get("database") != str(database.path.resolve())):
        raise ValueError("An administrator and an unchanged plan for this database are required.")
    backup = Path(backup_path)
    if backup.exists() or backup.resolve() == database.path.resolve():
        raise ValueError("Choose a new, separate recovery backup.")
    backup.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    descriptor = os.open(backup, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    os.close(descriptor)
    with database.connection() as connection, sqlite3.connect(backup) as destination:
        connection.backup(destination)
        if destination.execute("PRAGMA quick_check").fetchone()[0] != "ok":
            raise ValueError("The recovery backup failed its integrity check.")
    backup.chmod(0o600)
    snapshots = {snapshot["path"]: snapshot for change in plan["observations"]
        for snapshot in json.loads(change["after_context"])["docking_link_reconciliation"]["snapshots"]}
    for snapshot in snapshots.values():
        if _snapshot(snapshot["path"])[1] != snapshot:
            raise ValueError("An original artifact changed after planning. Rebuild the plan.")
    for record in plan["structures"]:
        content, snapshot = _snapshot(record["source"]["path"])
        if snapshot != record["source"]:
            raise ValueError("An original coordinate artifact changed after planning.")
        if str(database._store_file(Path(snapshot["path"]).name, content)) != record["path"]:
            raise ValueError("The coordinate store disagrees with the reviewed plan.")
    with database.connection() as connection:
        connection.execute("BEGIN IMMEDIATE")
        for dependency in plan["dependencies"]:
            group, member = dependency["group"], dependency["member"]
            actual = connection.execute("SELECT * FROM molecule_groups WHERE id=?", (group["id"],)).fetchone()
            saved = connection.execute("SELECT * FROM molecule_members WHERE group_id=? AND molecule_id=?",
                (member["group_id"], member["molecule_id"])).fetchone()
            if (not actual or dict(actual) != group or json.loads(actual["origin"]).get("_replica_node")
                    or not database._can_read(actual, user) or not saved or dict(saved) != member):
                raise ValueError("A source-owned group or target membership changed. Rebuild the plan.")
        for change in plan["observations"]:
            actual = connection.execute("SELECT * FROM molecule_observations WHERE id=?", (change["before"]["id"],)).fetchone()
            if not actual or dict(actual) != change["before"]:
                raise ValueError("A reviewed observation changed. Rebuild the plan.")
        fields = ("id", "molecule_id", "group_id", "path", "label", "format", "context", "created_at")
        for record in plan["structures"]:
            actual = connection.execute("SELECT * FROM molecule_structures WHERE id=?", (record["id"],)).fetchone()
            if actual and any(actual[key] != record[key] for key in fields if key not in {"created_at", "label"}):
                raise ValueError("A reviewed coordinate identity changed. Rebuild the plan.")
            connection.execute("INSERT OR IGNORE INTO molecule_structures VALUES (?,?,?,?,?,?,?,?)", tuple(record[key] for key in fields))
        for change in plan["observations"]:
            connection.execute("UPDATE molecule_observations SET context=? WHERE id=?", (change["after_context"], change["before"]["id"]))
    return {"plan": plan["digest"], "backup": str(backup), "observations": len(plan["observations"]),
        "structures": len(plan["structures"]), "unresolved": len(plan["unresolved"])}
