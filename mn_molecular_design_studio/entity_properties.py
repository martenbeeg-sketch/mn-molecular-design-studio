"""Typed, revision-checked descriptive metadata; never a scientific record editor."""
from __future__ import annotations

import hashlib
import json

from mn_molecular_design_studio import entity_names as names

FIELDS = {
    "project": {"description": "Description", "tags": "Tags", "notes": "Notes"},
    "dataset": {"description": "Description", "source": "Source / supplier", "citation": "Citation / reference", "tags": "Tags", "notes": "Notes"},
    "target": {"description": "Description / construct annotation", "identifiers": "Reference identifiers (annotations only)", "tags": "Tags", "notes": "Notes"},
    "entry": {"dataset_identifier": "Dataset-specific identifier", "supplier": "Supplier", "batch": "Batch / sample", "tags": "Tags", "notes": "Notes"},
    "result": {"description": "Description", "tags": "Tags", "notes": "Notes"},
}


def validate_properties(kind, properties):
    if kind not in FIELDS or not isinstance(properties, dict) or set(properties) - set(FIELDS[kind]):
        raise ValueError("Only the listed descriptive fields can be edited. Scientific identities, links, scores and admission decisions are protected.")
    for field, value in properties.items():
        if not isinstance(value, str) or len(value) > (12000 if field == "notes" else 4000) or "\x00" in value:
            raise ValueError("Property values must be text within the displayed length limit.")
    return {field: value.strip() for field, value in properties.items()}


def group_subject(group):
    if group.get("properties_subject_id"):
        return group["properties_subject_id"]
    from mn_molecular_design_studio.cluster import config
    origin = group.get("origin") or {}
    return json.dumps([origin.get("_replica_node") or config.identity(), origin.get("_replica_group") or group["id"]], separators=(",", ":"))


def entry_subject(group, molecule_id):
    return hashlib.sha256(json.dumps([group_subject(group), molecule_id], separators=(",", ":")).encode()).hexdigest()


def label_for(labels, kind, subject_id, scope):
    return next((row for row in labels if row["kind"] == kind and row["subject_id"] == subject_id and row["scope"] == scope), None)


def project_snapshot(user, project):
    if not user.get("id") or user.get("disabled"):
        raise PermissionError("An active signed-in account is required.")
    if user.get("role") != "admin" and project.get("created_by") != user["id"]:
        from mn_molecular_design_studio.auth import require_project_access
        require_project_access(user, project["id"])
    writable = bool(user.get("id") and not user.get("disabled") and (
        user.get("role") == "admin" or project.get("created_by") == user["id"]))
    if project.get("deleted_at") or (project.get("lifecycle") or {}).get("status") in {"preparing", "running", "attention"}:
        writable = False
    source_revision = hashlib.sha256(json.dumps([project.get("name"), project.get("description"), project.get("entity_properties", {})], sort_keys=True).encode()).hexdigest()
    return {"kind": "project", "subject_id": project["id"], "scope": "project:" + project["id"], "source_revision": source_revision,
        "name": project["name"], "revision": project.get("name_revision", ""), "aliases": project.get("name_aliases", []),
        "properties": {"description": project.get("description", ""), **project.get("entity_properties", {})},
        "authority_node": project.get("name_authority", ""), "writable": writable,
        "type_label": "Project", "scope_label": "This project, wherever it is referenced.",
        "source": {"Project ID": project["id"], "Created": project.get("created_at", "Not recorded")}}


def _snapshot(database, user, group_id, kind, subject_id, name, defaults, source):
    if not user.get("id") or user.get("disabled"):
        raise PermissionError("An active signed-in account is required.")
    group = database.group(group_id, user)
    with database.connection() as connection:
        label = label_for(names.records(connection), kind, subject_id, names.group_scope(group)) or {}
    writable = bool(user.get("id") and not user.get("disabled"))
    try:
        database.group(group_id, user, write=True)
    except PermissionError:
        writable = False
    return {"kind": kind, "subject_id": subject_id, "scope": names.group_scope(group), "name": label.get("name", name),
        "source_revision": hashlib.sha256(json.dumps([name, defaults, source], sort_keys=True).encode()).hexdigest(),
        "revision": label.get("revision", ""), "aliases": label.get("aliases", []),
        "properties": {**defaults, **label.get("properties", {})}, "authority_node": label.get("authority_node", ""),
        "writable": writable, "source": source, "group_id": group_id}


def dataset_snapshot(database, user, group_id):
    group = database.group(group_id, user)
    snapshot = _snapshot(database, user, group_id, "dataset", group_subject(group), group["name"], {},
        {"Collection ID": (group.get("origin") or {}).get("_replica_group") or group["id"],
         "Import method": group["import_method"], "Created": group["created_at"]})
    snapshot.update(type_label="Dataset / collection", scope_label="This dataset's descriptive properties only. Its members, project assignment and original imports do not change.")
    with database.connection() as connection:
        snapshot["source"]["Members (read-only)"] = connection.execute("SELECT count(*) FROM molecule_members WHERE group_id=?", (group_id,)).fetchone()[0]
        snapshot["source"]["Original import files"] = [dict(row) for row in connection.execute(
            "SELECT source_name,source_sha256,created_at FROM molecule_import_batches WHERE group_id=? ORDER BY created_at DESC LIMIT 10", (group_id,))]
    return snapshot


def entry_snapshot(database, user, group_id, molecule_id):
    group = database.group(group_id, user)
    with database.connection() as connection:
        member = connection.execute("SELECT name,metadata FROM molecule_members WHERE group_id=? AND molecule_id=?", (group_id, molecule_id)).fetchone()
    if member is None:
        raise PermissionError("This molecule is no longer in the selected dataset.")
    metadata = json.loads(member["metadata"])
    snapshot = _snapshot(database, user, group_id, "entry", entry_subject(group, molecule_id), member["name"],
        {field: str(metadata.get(field) or "") for field in ("tags", "notes")},
        {"Molecule ID": molecule_id, "Original entry name": member["name"],
         "Original identifiers": {field: metadata[field] for field in ("compound_id", "pubchem_cid", "source_id") if field in metadata}})
    snapshot.update(type_label="Dataset entry", scope_label="Only this molecule's entry in “" + group["name"] + "”. Other datasets and shared compound/target details are unchanged.", molecule_id=molecule_id)
    return snapshot


def target_snapshot(database, user, subject_id, scope):
    from mn_molecular_design_studio.database_browser import build_browser_index, component_subjects
    if not user.get("id") or user.get("disabled"):
        raise PermissionError("An active signed-in account is required.")
    index = build_browser_index(database, user, collection="all", hide_smoke_test=False, include_related=True)
    subject = next((row for row in component_subjects(index, role="targets") if row["id"] == subject_id), None)
    if subject is None:
        raise PermissionError("This target / compound is no longer accessible.")
    groups = [index["groups"][group_id] for group_id in subject["group_ids"] if names.group_scope(index["groups"][group_id]) == scope]
    if not groups:
        raise PermissionError("This naming scope is no longer accessible.")
    label = label_for(index.get("entity_names", []), "target", subject_id, scope) or {}
    return {"kind": "target", "subject_id": subject_id, "scope": scope, "name": label.get("name", subject["name"]),
        "revision": label.get("revision", ""), "aliases": label.get("aliases", []), "properties": label.get("properties", {}),
        "authority_node": label.get("authority_node", ""), "writable": scope in names.target_options(index, subject, user),
        "type_label": "Compound" if subject["type"] == "small_molecule" else "Target / construct",
        "scope_label": "Shared descriptive metadata for this exact molecular identity within " + (
            "this personal library." if scope.startswith("owner:") else "this project. Other projects are unchanged."),
        "source": {"Molecular type": subject["type"], "Identity type": subject["identity"][1], "Identity": subject["identity"][2]}}


def result_snapshot(database, user, outcome):
    from mn_molecular_design_studio.database_browser import build_browser_index, subject_results
    from mn_molecular_design_studio.data_browser import saved_results
    group = database.group(outcome["group_id"], user)
    index = build_browser_index(database, user, collection="all", project_id=group["project_id"], hide_smoke_test=False, include_related=True)
    subject = {"id": "properties-verification", "group_ids": set(index["groups"]), "molecule_ids": {outcome["molecule_id"]}}
    current = next((row for row in saved_results(index, subject_results(database, index, subject)) if row["id"] == outcome["id"]), None)
    if current is None or not any(row["group_id"] == outcome["group_id"] for row in current["originals"]):
        raise PermissionError("This saved result is no longer accessible in the selected context.")
    from mn_molecular_design_studio.data_browser import concise_result_label
    snapshot = _snapshot(database, user, group["id"], "result", outcome["id"], concise_result_label(current), {},
        {"Result ID": outcome["id"], "Molecule ID": outcome["molecule_id"], "Activity": current["activity"],
         "Admission (read-only)": current["eligibility"], "Methods (read-only)": current["origin_tools"]})
    snapshot.update(type_label="Saved result / design", scope_label="This saved scientific outcome's display label and annotations. Structures, models, calculations and filter decisions stay unchanged.")
    return snapshot


def save_group_properties(database, user, snapshot, name, properties, aliases, *, expected_revision, outcome=None):
    from mn_molecular_design_studio.cluster.replication import locked
    from mn_molecular_design_studio.molecular_database import MolecularDatabase
    if not user.get("id") or user.get("disabled"):
        raise PermissionError("Sign in with an active account before editing metadata.")
    writer = MolecularDatabase(database.root)
    with locked():
        kind = snapshot["kind"]
        if kind == "dataset":
            current = dataset_snapshot(writer, user, snapshot["group_id"])
        elif kind == "entry":
            current = entry_snapshot(writer, user, snapshot["group_id"], snapshot["molecule_id"])
        elif kind == "result" and outcome is not None:
            current = result_snapshot(writer, user, outcome)
        else:
            raise ValueError("Unsupported metadata object.")
        writer.group(current["group_id"], user, write=True)
        if (current["subject_id"], current["scope"]) != (snapshot["subject_id"], snapshot["scope"]):
            raise names.NameConflictError("This object moved or changed. Reopen its properties before editing.")
        if current["source_revision"] != snapshot.get("source_revision"):
            raise names.NameConflictError("Original source annotations changed in another session. Reload current values before saving.")
        with writer.connection() as connection:
            connection.execute("BEGIN IMMEDIATE")
            value = names._write(connection, kind, current["subject_id"], current["scope"], name, user,
                                 expected_revision, [current["name"], *aliases], properties=properties)
        names.publish(database=writer)
        return value


def project_groups(groups, labels):
    lookup = {(row["kind"], row["subject_id"], row["scope"]): row for row in labels}
    for group in groups:
        group["properties_subject_id"] = group_subject(group)
        label = lookup.get(("dataset", group_subject(group), names.group_scope(group)))
        if label:
            group.update(original_name=group["name"], name=label["name"], name_aliases=label["aliases"],
                         entity_properties=label.get("properties", {}))
    return groups


def project_results(index, outcomes):
    labels = index.get("entity_names", [])
    if not labels:
        return outcomes
    lookup = {(row["kind"], row["subject_id"], row["scope"]): row for row in labels}
    for outcome in outcomes:
        group = index["groups"][outcome["group_id"]]
        label = lookup.get(("result", outcome["id"], names.group_scope(group)))
        if label:
            outcome["display_name"] = label["name"]
            outcome.setdefault("aliases", []).extend([outcome["name"], *label["aliases"], *label.get("properties", {}).values()])
        compounds = []
        for component in index["molecules"][outcome["molecule_id"]]["components"]:
            profile = lookup.get(("target", names.component_id(component), names.group_scope(group)))
            if profile:
                outcome.setdefault("aliases", []).extend([profile["name"], *profile["aliases"], *profile.get("properties", {}).values()])
                if component["type"] == "small_molecule":
                    compounds.append(profile)
        if not label and len(compounds) == 1 and outcome["outcome_kind"] != "Design":
            outcome["display_name"] = compounds[0]["name"]
        entry = lookup.get(("entry", entry_subject(group, outcome["molecule_id"]), names.group_scope(group)))
        if entry:
            outcome.setdefault("aliases", []).extend([entry["name"], *entry["aliases"], *entry.get("properties", {}).values()])
        if "aliases" in outcome:
            outcome["aliases"] = list(dict.fromkeys(outcome["aliases"]))
    return outcomes
