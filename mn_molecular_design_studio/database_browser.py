"""Read-only navigation across exact molecular contents and recorded run lineage."""
from __future__ import annotations

import hashlib
import json
import socket
from pathlib import PurePosixPath


def _key(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _run_reference(provenance):
    run_id = str(provenance.get("run_id") or "")
    task_group = str(provenance.get("task_group") or "")
    directory = str(provenance.get("run_dir") or "")
    if directory:
        path = PurePosixPath(directory)
        if not run_id:
            run_id = path.name
        if not task_group and path.name == run_id:
            task_group = path.parent.name
    return task_group, run_id


def _campaign_reference(provenance, node, owner, project):
    for field in ("parent_campaign_run_dir", "pipeline_parent_run_dir"):
        if directory := str(provenance.get(field) or ""):
            path = PurePosixPath(directory)
            if path.is_absolute() and path.name and path.parent.name:
                return (node, path.parent.name, path.name, owner, project)
    if campaign_id := str(provenance.get("campaign_id") or ""):
        return (node, "campaign", campaign_id, owner, project)
    return None


def build_browser_index(database, user, *, project_id="", scope="all", methods=None,
                        collection="library", hide_smoke_test=True):
    visible = database.groups(user, project_id=project_id, scope=scope, methods=methods,
                              collection=collection, hide_smoke_test=hide_smoke_test)
    accessible = database.groups(user, project_id=project_id, scope=scope, collection="all",
                                 hide_smoke_test=hide_smoke_test)
    groups = {row["id"]: row for row in accessible}
    visible_json = json.dumps([row["id"] for row in visible])
    accessible_json = json.dumps(list(groups))
    if not visible:
        return {"molecules": {}, "runs": {}, "groups": groups}
    node_file = database.root.parent / "cluster/node-id"
    local_node = node_file.read_text().strip() if node_file.is_file() else "local"
    molecules = {}
    runs = {}
    with database.connection() as connection:
        connection.execute("BEGIN")
        members = connection.execute(
            "SELECT m.id,m.components,m.created_at,mm.group_id,mm.name FROM molecules m "
            "JOIN molecule_members mm ON mm.molecule_id=m.id "
            "WHERE mm.group_id IN (SELECT value FROM json_each(?)) AND EXISTS "
            "(SELECT 1 FROM molecule_members selected WHERE selected.molecule_id=m.id "
            "AND selected.group_id IN (SELECT value FROM json_each(?))) "
            "ORDER BY m.id,mm.added_at DESC,mm.group_id", (accessible_json, visible_json),
        ).fetchall()
        for row in members:
            molecule = molecules.setdefault(row["id"], {
                "id": row["id"], "name": row["name"], "names": set(),
                "components": json.loads(row["components"]), "created_at": row["created_at"],
                "group_ids": set(), "run_keys": set(),
            })
            molecule["names"].add(row["name"])
            molecule["group_ids"].add(row["group_id"])
        ids_json = json.dumps(list(molecules))
        projections = []
        for table, column, run_column in (("molecule_observations", "context", "run_id"),
                                           ("molecule_sources", "provenance", "NULL")):
            projections.append(connection.execute(
                "SELECT DISTINCT molecule_id,group_id," + run_column + " AS recorded_run_id,"
                "json_extract(" + column + ",'$.run_id') AS run_id,"
                "json_extract(" + column + ",'$.task_group') AS task_group,"
                "json_extract(" + column + ",'$.run_dir') AS run_dir,"
                "json_extract(" + column + ",'$._replica_node') AS _replica_node,"
                "json_extract(" + column + ",'$.source_node') AS source_node,"
                "json_extract(" + column + ",'$.source_node_id') AS source_node_id,"
                "json_extract(" + column + ",'$.workstation') AS workstation,"
                "json_extract(" + column + ",'$.source_workstation') AS source_workstation,"
                "json_extract(" + column + ",'$.tool') AS tool,"
                "json_extract(" + column + ",'$.source_tool') AS source_tool,"
                "json_extract(" + column + ",'$.campaign_id') AS campaign_id,"
                "json_extract(" + column + ",'$.campaign_name') AS campaign_name,"
                "json_extract(" + column + ",'$.parent_campaign_run_dir') AS parent_campaign_run_dir,"
                "json_extract(" + column + ",'$.pipeline_parent_run_dir') AS pipeline_parent_run_dir "
                "FROM " + table + " WHERE molecule_id IN (SELECT value FROM json_each(?)) "
                "AND group_id IN (SELECT value FROM json_each(?))", (ids_json, accessible_json),
            ).fetchall())

    recorded_pairs = set()

    def add_run(molecule_id, group_id, provenance):
        group = groups[group_id]
        origin = group["origin"]
        task_group, run_id = _run_reference(provenance)
        if not run_id:
            return
        origin_reference = _run_reference(origin)
        if not task_group and origin_reference[1] == run_id:
            task_group = origin_reference[0]
        workstation = str(provenance.get("workstation") or provenance.get("source_workstation") or origin.get("workstation") or "")
        fallback_node = ("workstation:" + workstation.casefold() if workstation
                         and workstation.split(".")[0].casefold() != socket.gethostname().split(".")[0].casefold()
                         else local_node)
        node = str(provenance.get("_replica_node") or provenance.get("source_node") or provenance.get("source_node_id")
                   or origin.get("_replica_node") or origin.get("source_node") or fallback_node)
        reference = (node, task_group, run_id, group["owner_id"], group["project_id"])
        run_key = _key(reference if task_group else (*reference, group_id))
        run = runs.setdefault(run_key, {
            "id": run_key, "reference": reference, "name": run_id, "run_id": run_id,
            "task_group": task_group, "node": node, "owner_id": group["owner_id"],
            "project_id": group["project_id"], "molecule_ids": set(), "group_ids": set(),
            "tools": set(), "workstations": set(), "campaign_references": set(),
            "campaign_names": set(), "updated_at": "",
        })
        run["molecule_ids"].add(molecule_id)
        run["group_ids"].add(group_id)
        run["updated_at"] = max(run["updated_at"], group.get("updated_at") or "")
        if tool := str(provenance.get("tool") or provenance.get("source_tool") or ""):
            run["tools"].add(tool)
        run["workstations"].add(workstation or node)
        lineage_sources = [provenance]
        if origin_reference == (task_group, run_id):
            lineage_sources.append(origin)
            if group.get("collection_kind") not in {"evaluation", "inputs"}:
                run["name"] = min(run["name"], group["name"]) if run["name"] != run_id else group["name"]
        for source in lineage_sources:
            if campaign := _campaign_reference(source, node, group["owner_id"], group["project_id"]):
                run["campaign_references"].add(campaign)
                if name := str(source.get("campaign_name") or ""):
                    run["campaign_names"].add(name)
        molecules[molecule_id]["run_keys"].add(run_key)
        recorded_pairs.add((molecule_id, group_id))

    for rows in projections:
        for row in rows:
            provenance = dict(row)
            provenance["run_id"] = provenance.pop("recorded_run_id") or provenance["run_id"]
            add_run(row["molecule_id"], row["group_id"], provenance)
    for molecule_id, molecule in molecules.items():
        for group_id in sorted(molecule["group_ids"]):
            group = groups[group_id]
            if ((molecule_id, group_id) in recorded_pairs or group.get("collection_kind") in {"evaluation", "inputs"}
                    or group.get("import_method") == "tool_inputs"):
                continue
            origin = group["origin"]
            add_run(molecule_id, group_id, origin.get("source_run") or origin)
    return {"molecules": molecules, "runs": runs, "groups": groups}


def molecular_subjects(index, *, individual_sequences=False):
    subjects = {}
    for molecule_id, molecule in index["molecules"].items():
        components = molecule["components"]
        contents = [(row["type"], "sequence" if row.get("sequence") else "smiles" if row.get("smiles") else "ccd",
                     row.get("sequence") or row.get("smiles") or row.get("ccd") or "") for row in components]
        if individual_sequences:
            identities = sorted({(item,) for item in contents if item[1] == "sequence" and item[2]})
        else:
            identities = [tuple(sorted(contents))]
        for identity in identities:
            complete = bool(identity) and all(item[2] for item in identity)
            key = _key(identity if complete else (identity, molecule_id))
            subject = subjects.setdefault(key, {
                "id": key, "name": molecule["name"], "names": set(), "molecule_ids": set(),
                "group_ids": set(), "run_keys": set(), "updated_at": "", "identity": identity,
                "type": " / ".join(sorted({item[0] for item in identity})),
                "sequence": " · ".join(item[2] for item in identity),
            })
            subject["names"].update(molecule["names"])
            subject["molecule_ids"].add(molecule_id)
            subject["group_ids"].update(molecule["group_ids"])
            subject["run_keys"].update(molecule["run_keys"])
            subject["updated_at"] = max(subject["updated_at"], molecule["created_at"])
    for subject in subjects.values():
        subject["name"] = min(subject["names"], key=lambda value: (value.casefold(), value))
    return list(subjects.values())


def run_subjects(index, *, combine_campaigns=True):
    references = {tuple(run["reference"]): key for key, run in index["runs"].items() if run["task_group"]}

    def campaign_root(run_key):
        visited = set()
        current = run_key
        while current not in visited:
            visited.add(current)
            run = index["runs"][current]
            candidates = run["campaign_references"]
            if len(candidates) != 1:
                return run_key
            reference = next(iter(candidates))
            parent = references.get(reference)
            if parent is None or parent == current:
                return _key(reference)
            current = parent
            if not index["runs"][current]["campaign_references"]:
                return current
        return run_key

    subjects = {}
    for run_key, run in index["runs"].items():
        key = campaign_root(run_key) if combine_campaigns else run_key
        subject = subjects.setdefault(key, {
            "id": key, "name": run["name"], "names": set(), "molecule_ids": set(),
            "group_ids": set(), "run_keys": set(), "updated_at": "", "tools": set(),
            "workstations": set(),
        })
        subject["names"].add(run["name"])
        subject["names"].update(run["campaign_names"])
        subject["molecule_ids"].update(run["molecule_ids"])
        subject["group_ids"].update(run["group_ids"])
        subject["run_keys"].add(run_key)
        subject["tools"].update(run["tools"])
        subject["workstations"].update(run["workstations"])
        subject["updated_at"] = max(subject["updated_at"], run["updated_at"])
    for subject in subjects.values():
        parent = index["runs"].get(subject["id"])
        subject["name"] = parent["name"] if parent else min(subject["names"], key=lambda value: (value.casefold(), value))
    return list(subjects.values())


def browse_subjects(subjects, index, *, query="", sort="Most jobs"):
    terms = query.casefold().split()
    result = []
    for subject in subjects:
        run_rows = [index["runs"][key] for key in subject["run_keys"]]
        text = " ".join([*subject["names"], subject.get("sequence", ""),
                         *subject["molecule_ids"], *subject["group_ids"],
                         *[index["groups"][key]["name"] for key in subject["group_ids"]],
                         *[str(value) for row in run_rows for value in
                           (row["run_id"], row["task_group"], *row["tools"], *row["workstations"])]]).casefold()
        if all(term in text for term in terms):
            result.append(subject)
    result.sort(key=lambda row: (row["name"].casefold(), row["id"]))
    if sort == "Most jobs":
        result.sort(key=lambda row: len(row["run_keys"]), reverse=True)
    elif sort == "Most molecules":
        result.sort(key=lambda row: len(row["molecule_ids"]), reverse=True)
    elif sort == "Recently added":
        result.sort(key=lambda row: row["updated_at"], reverse=True)
    return result
