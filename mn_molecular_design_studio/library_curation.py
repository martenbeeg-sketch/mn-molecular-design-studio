"""Reviewable, backed-up naming and exact-target corrections for existing SQL data."""
from __future__ import annotations

from collections import defaultdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sqlite3

from mn_molecular_design_studio.molecule_imports import canonical_components, parse_structure
from mn_molecular_design_studio.scientific_naming import design_name, meaningful_name, named_metadata, technical_design_name
from mn_molecular_design_studio.target_registration import contents_identity, registration_evidence, target_components


def encoded(value):
    return json.dumps(value, sort_keys=True, ensure_ascii=False, separators=(",", ":"))


def _digest(value):
    return hashlib.sha256(encoded(value).encode()).hexdigest()


def _design_group(group):
    return group["collection_kind"] == "generated" or (group["collection_kind"] == "curated"
        and group["origin"].get("task_group") == "design" and bool(group["origin"].get("run_id")))


def plan_curation(database, user, *, reference_names=None, reference_constructs=()):
    if user.get("role") != "admin":
        raise PermissionError("Library maintenance requires an administrator.")
    reference_names = reference_names or {}
    groups = {row["id"]: row for row in database.groups(user, scope="all", collection="all") if not row["read_only"]}
    with database.connection() as connection:
        connection.execute("BEGIN")
        molecules = {row["id"]: json.loads(row["components"]) for row in connection.execute("SELECT id,components FROM molecules")}
        members = [dict(row) for row in connection.execute("SELECT * FROM molecule_members WHERE group_id IN (SELECT value FROM json_each(?))", (encoded(list(groups)),))]
        observations = [dict(row) for row in connection.execute("SELECT * FROM molecule_observations WHERE group_id IN (SELECT value FROM json_each(?))", (encoded(list(groups)),))]
    by_member = {(row["group_id"], row["molecule_id"]): row for row in members}
    by_observation = defaultdict(list)
    for observation in observations:
        by_observation[(observation["group_id"], observation["molecule_id"])].append(observation)
    resolutions, unresolved, targets = {}, [], {}
    parsed = {}
    for observation in observations:
        group = groups[observation["group_id"]]
        if not _design_group(group):
            continue
        member = by_member.get((observation["group_id"], observation["molecule_id"]))
        if not member:
            continue
        context, metadata = json.loads(observation["context"]), json.loads(member["metadata"])
        candidate = context.get("candidate") or {}
        original_target = context.get("target_molecule_id") or context.get("target_id") or metadata.get("design_target_molecule_id") or ""
        components = target_components(json.loads(member["components"]), metadata, candidate)
        if not components and original_target in molecules:
            components = molecules[original_target]
        if not components:
            source = str(candidate.get("complex_pdb") or "")
            model = candidate.get("model_id") or context.get("model_id")
            key = (source, str(model))
            if source and Path(source).is_file():
                if key not in parsed:
                    try:
                        parsed[key] = parse_structure(Path(source).read_bytes(), Path(source).name,
                            model_id=model, require_model_selection=True)[0]["components"]
                    except (OSError, ValueError, KeyError):
                        parsed[key] = []
                components = target_components(parsed[key], metadata, candidate)
        if not components:
            unresolved.append({"observation_id": observation["id"], "target_id": original_target, "reason": "No exact, explicitly selected target contents"})
            continue
        normalized, fingerprint, _ = canonical_components(components)
        target_id = "mol-" + fingerprint[:32]
        if original_target in molecules and contents_identity(molecules[original_target]) != contents_identity(normalized):
            unresolved.append({"observation_id": observation["id"], "target_id": original_target, "reason": "Recorded target conflicts with selected target contents"})
            continue
        names = [reference_names.get(_digest(contents_identity(normalized))),
                 metadata.get("design_target_name"), (context.get("settings") or {}).get("target_name")]
        if len(normalized) == 1 and normalized[0]["type"] == "protein":
            sequence = normalized[0].get("sequence") or ""
            matches = [(reference, reference["sequence"].find(sequence)) for reference in reference_constructs
                       if sequence and reference["sequence"].count(sequence) == 1]
            if len(matches) == 1:
                reference, offset = matches[0]
                if sequence == reference["sequence"]:
                    names.insert(0, reference["name"])
                else:
                    start = int(reference.get("start_residue", 1)) + offset
                    names.insert(0, reference["name"] + f" · crop {start}–{start + len(sequence) - 1} ({len(sequence)} aa)")
        names += [row["name"] for row in members if row["molecule_id"] == target_id and (
            json.loads(row["metadata"]).get("target_identity_type") or groups[row["group_id"]]["collection_kind"] in {"reference", "inputs"})]
        label = next((name for name in names if meaningful_name(name)), "")
        if not label:
            unresolved.append({"observation_id": observation["id"], "target_id": original_target, "reason": "Exact target contents found, but no verified descriptive name"})
            continue
        evidence = registration_evidence(normalized, (context.get("target_registration") or {}).get("original_target_id", original_target))
        resolutions[observation["id"]] = {"target_id": target_id, "target_name": label, "evidence": evidence}
        targets[(group["project_id"], group["owner_id"], target_id)] = {"id": target_id, "name": label,
            "components": normalized, "project_id": group["project_id"], "owner_id": group["owner_id"]}
    member_changes, observation_changes = [], []
    for observation in observations:
        resolution = resolutions.get(observation["id"])
        if not resolution:
            continue
        context = json.loads(observation["context"])
        updated = {**context, "target_molecule_id": resolution["target_id"], "target_registration": resolution["evidence"]}
        if context != updated:
            observation_changes.append({"id": observation["id"], "before": observation["context"], "after": encoded(updated)})
    for member in sorted(members, key=lambda row: (row["added_at"], row["name"], row["molecule_id"])):
        group = groups[member["group_id"]]
        metadata = json.loads(member["metadata"])
        own = by_observation[(member["group_id"], member["molecule_id"])]
        resolved = [resolutions[row["id"]] for row in own if row["id"] in resolutions]
        if not _design_group(group) or not resolved or len({row["target_id"] for row in resolved}) != 1:
            continue
        if any(row["id"] not in resolutions for row in own):
            continue
        resolution = resolved[0]
        updated = {**metadata, "design_target_molecule_id": resolution["target_id"], "design_target_name": resolution["target_name"]}
        name = member["name"]
        naming = metadata.get("studio_naming") or {}
        if (technical_design_name(name) or name == naming.get("display_name")) and naming.get("mode") != "manual":
            origin = group["origin"]
            run_ids = sorted({row["run_id"] for row in own})
            tool = metadata.get("generator_engine") or metadata.get("source_tool") or origin.get("tool") or ""
            run_id = origin.get("run_id") or (run_ids[0] if run_ids else "")
            label = design_name(name, target_name=resolution["target_name"], tool=tool, run_id=run_id, metadata=metadata)
            if label:
                updated = named_metadata(name, label, updated, target_name=resolution["target_name"], target_id=resolution["target_id"], tool=tool, run_id=run_id)
                name = label
        if name != member["name"] or updated != metadata:
            member_changes.append({"group_id": member["group_id"], "molecule_id": member["molecule_id"],
                "before_name": member["name"], "after_name": name, "before_metadata": member["metadata"], "after_metadata": encoded(updated)})
    target_names = {target["id"]: target["name"] for target in targets.values()}
    changed_members = {(change["group_id"], change["molecule_id"]) for change in member_changes}
    for member in members:
        name = target_names.get(member["molecule_id"])
        if not name or meaningful_name(member["name"]) or (member["group_id"], member["molecule_id"]) in changed_members:
            continue
        metadata = json.loads(member["metadata"])
        components = json.loads(member["components"])
        if metadata.get("output_kind") or any(row.get("role") in {"binder", "designed_binder", "design"} for row in components):
            continue
        updated = named_metadata(member["name"], name, metadata, target_name=name, target_id=member["molecule_id"], tool="", run_id="")
        member_changes.append({"group_id": member["group_id"], "molecule_id": member["molecule_id"],
            "before_name": member["name"], "after_name": name, "before_metadata": member["metadata"], "after_metadata": encoded(updated)})
    plan = {"schema": "studio.library-curation.v1", "database": str(database.path.resolve()), "members": member_changes,
        "observations": observation_changes, "targets": list(targets.values()), "unresolved": unresolved}
    plan["digest"] = _digest(plan)
    return plan


def apply_curation(database, user, plan, *, backup_path):
    if user.get("role") != "admin" or plan.get("digest") != _digest({key: value for key, value in plan.items() if key != "digest"}):
        raise ValueError("An administrator and an unchanged reviewed curation plan are required.")
    if str(database.path.resolve()) != plan["database"]:
        raise ValueError("The reviewed plan belongs to a different database.")
    backup_path = Path(backup_path)
    backup_path.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    if backup_path.exists() or backup_path.resolve() == database.path.resolve():
        raise ValueError("Choose a new, separate recovery backup path.")
    with database.connection() as source, sqlite3.connect(backup_path) as destination:
        source.backup(destination)
        if destination.execute("PRAGMA quick_check").fetchone()[0] != "ok":
            raise ValueError("The recovery backup failed its integrity check.")
    backup_path.chmod(0o600)
    if plan.get("schema") == "studio.library-design-creation.v1":
        verified = set()
        for change in plan["members"]:
            evidence = json.loads(change["after_metadata"])["design_creation"]
            source = Path(evidence["source_path"])
            key = (str(source), evidence["source_sha256"], evidence.get("source_mtime_ns"))
            if key in verified:
                continue
            if not source.is_file() or hashlib.sha256(source.read_bytes()).hexdigest() != evidence["source_sha256"] or (
                evidence.get("source_mtime_ns") is not None and source.stat().st_mtime_ns != evidence["source_mtime_ns"]
            ):
                raise ValueError("A reviewed original design artifact changed or became unavailable. Rebuild the plan.")
            verified.add(key)
    if plan.get("schema") == "studio.library-target-links.v1":
        verified = set()
        for change in plan["observations"]:
            evidence = json.loads(change["after"])["target_registration"]
            for snapshot in evidence["snapshots"]:
                source = Path(snapshot["path"])
                key = (str(source), snapshot["sha256"])
                if key in verified:
                    continue
                if not source.is_file() or hashlib.sha256(source.read_bytes()).hexdigest() != snapshot["sha256"]:
                    raise ValueError("An original target input changed or became unavailable. Rebuild the plan.")
                verified.add(key)
    with database.connection() as connection:
        connection.execute("BEGIN IMMEDIATE")
        for change in plan["members"]:
            actual = connection.execute("SELECT mm.name,mm.metadata,g.origin FROM molecule_members mm JOIN molecule_groups g ON g.id=mm.group_id WHERE mm.group_id=? AND mm.molecule_id=?",
                (change["group_id"], change["molecule_id"])).fetchone()
            if actual is None or actual["name"] != change["before_name"] or actual["metadata"] != change["before_metadata"] or json.loads(actual["origin"]).get("_replica_node"):
                raise ValueError("A reviewed molecular member changed or is owned by another workstation. Rebuild the plan.")
        for change in plan["observations"]:
            actual = connection.execute("SELECT o.context,g.origin FROM molecule_observations o JOIN molecule_groups g ON g.id=o.group_id WHERE o.id=?", (change["id"],)).fetchone()
            if actual is None or actual["context"] != change["before"] or json.loads(actual["origin"]).get("_replica_node"):
                raise ValueError("A reviewed observation changed or is owned by another workstation. Rebuild the plan.")
        stamp = datetime.now(timezone.utc).isoformat(timespec="seconds")
        for target in plan["targets"]:
            group_id = "inputs-" + hashlib.sha256(f"{target['owner_id']}:{target['project_id']}".encode()).hexdigest()[:32]
            existing = connection.execute("SELECT owner_id,project_id,origin FROM molecule_groups WHERE id=?", (group_id,)).fetchone()
            if existing and (existing["owner_id"] != target["owner_id"] or existing["project_id"] != target["project_id"] or json.loads(existing["origin"]).get("_replica_node")):
                raise ValueError("The target input collection has conflicting ownership.")
            connection.execute("INSERT OR IGNORE INTO molecule_groups VALUES (?,?,?,?,?,?,?,?,?)", (group_id, "Tool inputs", target["owner_id"], target["project_id"],
                "tool_inputs", encoded({"type": "job", "collection_kind": "inputs"}), "{}", stamp, stamp))
            database._put(connection, group_id, {"name": target["name"], "components": target["components"],
                "metadata": {"target_identity_type": "molecular_identity", "structural_reference": False,
                             "target_registration": {"policy": "studio.exact-target-registration.v1"}}},
                source_key="registered-target:" + target["project_id"] + ":" + target["id"],
                origin={"input_role": "modeled_target", "registration_plan": plan["digest"]}, preserve_existing=True)
        for change in plan["members"]:
            connection.execute("UPDATE molecule_members SET name=?,metadata=? WHERE group_id=? AND molecule_id=?",
                (change["after_name"], change["after_metadata"], change["group_id"], change["molecule_id"]))
        for change in plan["observations"]:
            connection.execute("UPDATE molecule_observations SET context=? WHERE id=?", (change["after"], change["id"]))
    return {"plan": plan["digest"], "backup": str(backup_path), "members": len(plan["members"]),
            "observations": len(plan["observations"]), "targets": len(plan["targets"]), "unresolved": len(plan["unresolved"])}
