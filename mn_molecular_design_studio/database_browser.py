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
                        collection="library", hide_smoke_test=True, include_related=False):
    visible = database.groups(user, project_id=project_id, scope=scope, methods=methods,
                              collection=collection, hide_smoke_test=hide_smoke_test)
    accessible = database.groups(user, project_id=project_id, scope=scope, collection="all",
                                 hide_smoke_test=hide_smoke_test)
    groups = {row["id"]: row for row in accessible}
    visible_json = json.dumps([row["id"] for row in visible])
    accessible_json = json.dumps(list(groups))
    if not visible:
        return {"molecules": {}, "runs": {}, "groups": groups, "visible_ids": set(), "evidence_runs": {}}
    node_file = database.root.parent / "cluster/node-id"
    local_node = node_file.read_text().strip() if node_file.is_file() else "local"
    molecules = {}
    runs = {}
    visible_ids = set()
    evidence_runs = {}
    with database.connection() as connection:
        connection.execute("BEGIN")
        members = connection.execute(
            "SELECT m.id,m.components,m.created_at,mm.group_id,mm.name,"
            "json_extract(mm.metadata,'$.design_target_molecule_id') AS design_target,"
            "json_extract(mm.metadata,'$.target_molecule_id') AS target,"
            "mm.group_id IN (SELECT value FROM json_each(?)) AS visible FROM molecules m "
            "JOIN molecule_members mm ON mm.molecule_id=m.id "
            "WHERE mm.group_id IN (SELECT value FROM json_each(?)) " + ("" if include_related else "AND EXISTS "
            "(SELECT 1 FROM molecule_members selected WHERE selected.molecule_id=m.id "
            "AND selected.group_id IN (SELECT value FROM json_each(?))) ") +
            "ORDER BY m.id,mm.added_at DESC,mm.group_id",
            [visible_json, accessible_json, *([] if include_related else [visible_json])],
        ).fetchall()
        for row in members:
            molecule = molecules.setdefault(row["id"], {
                "id": row["id"], "name": row["name"], "names": set(),
                "components": json.loads(row["components"]), "created_at": row["created_at"],
                "group_ids": set(), "run_keys": set(), "linked_to": set(), "metadata_links": set(),
            })
            molecule["names"].add(row["name"])
            molecule["group_ids"].add(row["group_id"])
            molecule["metadata_links"].update(value for value in (row["design_target"], row["target"]) if value)
            molecule["linked_to"].update(molecule["metadata_links"])
            if row["visible"]:
                visible_ids.add(row["id"])
        ids_json = json.dumps(list(molecules))
        projections = []
        for table, column, run_column, id_column in (("molecule_observations", "context", "run_id", "id"),
                                                    ("molecule_sources", "provenance", "NULL", "source_key"),
                                                    ("molecule_structures", "context", "NULL", "id")):
            projections.append((table, connection.execute(
                "SELECT " + id_column + " AS evidence_id,molecule_id,group_id," + run_column + " AS recorded_run_id,"
                "json_extract(" + column + ",'$.target_molecule_id') AS target_molecule_id,"
                "json_extract(" + column + ",'$.target_id') AS target_id,"
                "json_extract(" + column + ",'$.design_target_molecule_id') AS design_target_molecule_id,"
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
            ).fetchall()))

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
        return run_key

    for table, rows in projections:
        for row in rows:
            provenance = dict(row)
            provenance["run_id"] = provenance.pop("recorded_run_id") or provenance["run_id"]
            evidence_runs[(table, row["evidence_id"])] = add_run(row["molecule_id"], row["group_id"], provenance)
            molecules[row["molecule_id"]]["linked_to"].update(value for value in
                (row["target_molecule_id"], row["target_id"], row["design_target_molecule_id"]) if value in molecules)
    for molecule_id, molecule in molecules.items():
        for group_id in sorted(molecule["group_ids"]):
            group = groups[group_id]
            if ((molecule_id, group_id) in recorded_pairs or group.get("collection_kind") in {"evaluation", "inputs"}
                    or group.get("import_method") == "tool_inputs"):
                continue
            origin = group["origin"]
            evidence_runs[("membership", molecule_id, group_id)] = add_run(molecule_id, group_id, origin.get("source_run") or origin)
    return {"molecules": molecules, "runs": runs, "groups": groups,
            "visible_ids": visible_ids, "evidence_runs": evidence_runs}


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
    visible = index.get("visible_ids", set(index["molecules"]))
    return [subject for subject in subjects.values() if subject["molecule_ids"] & visible]


def component_subjects(index):
    """Select an exact ligand/polymer across complexes and explicitly linked designs."""
    subjects = {}
    linked_by_target = {}
    for molecule_id, molecule in index["molecules"].items():
        for target_id in molecule.get("linked_to", set()):
            linked_by_target.setdefault(target_id, set()).add(molecule_id)
        for component in molecule["components"]:
            field = next((field for field in ("sequence", "smiles", "ccd") if component.get(field)), "")
            if not field:
                continue
            identity = (component["type"], field, component[field])
            subject = subjects.setdefault(identity, {
                "id": _key(["component", identity]), "name": "", "names": set(), "preferred_names": set(),
                "molecule_ids": set(), "direct_molecule_ids": set(), "group_ids": set(), "run_keys": set(),
                "updated_at": "", "type": identity[0], "sequence": identity[2],
            })
            subject["direct_molecule_ids"].add(molecule_id)
            subject["names"].update(molecule["names"])
            if len(molecule["components"]) == 1:
                subject["preferred_names"].update(molecule["names"])
            elif component.get("name"):
                subject["preferred_names"].add(component["name"])
    visible = index.get("visible_ids", set(index["molecules"]))
    result = []
    for identity, subject in subjects.items():
        direct = subject["direct_molecule_ids"]
        if not direct & visible:
            continue
        subject["target_record_ids"] = set(direct)
        linked = set().union(*(linked_by_target.get(target_id, set()) for target_id in subject["target_record_ids"]))
        subject["molecule_ids"] = direct | linked
        for molecule_id in subject["molecule_ids"]:
            molecule = index["molecules"][molecule_id]
            subject["group_ids"].update(molecule["group_ids"])
            subject["run_keys"].update(molecule["run_keys"])
            subject["updated_at"] = max(subject["updated_at"], molecule["created_at"])
        subject["name"] = (min(subject["preferred_names"], key=lambda value: (value.casefold(), value))
                           if subject["preferred_names"] else identity[0].replace("_", " ").title() + " · " + identity[2][:20])
        subject["names"].add(subject["name"])
        result.append(subject)
    return result


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
    visible = index.get("visible_ids", set(index["molecules"]))
    return [subject for subject in subjects.values() if subject["molecule_ids"] & visible]


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


def subject_results(database, index, subject, *, run_view=False):
    """Read individual evidence rows; never combine scores from different evaluations."""
    from mn_molecular_design_studio.evaluations import evaluation_descriptor

    molecule_ids = json.dumps(sorted(subject["molecule_ids"]))
    group_ids = json.dumps(sorted(subject["group_ids"] & index["groups"].keys()))
    with database.connection() as connection:
        connection.execute("BEGIN")
        observations = connection.execute(
            "SELECT id,molecule_id,group_id,namespace,run_id,context,created_at FROM molecule_observations "
            "WHERE molecule_id IN (SELECT value FROM json_each(?)) "
            "AND group_id IN (SELECT value FROM json_each(?)) ORDER BY created_at DESC,id", (molecule_ids, group_ids),
        ).fetchall()
        structures = connection.execute(
            "SELECT id,molecule_id,group_id,label,context,created_at FROM molecule_structures "
            "WHERE molecule_id IN (SELECT value FROM json_each(?)) AND group_id IN (SELECT value FROM json_each(?))",
            (molecule_ids, group_ids),
        ).fetchall()
        sources = connection.execute(
            "SELECT source_key,molecule_id,group_id,provenance FROM molecule_sources "
            "WHERE molecule_id IN (SELECT value FROM json_each(?)) AND group_id IN (SELECT value FROM json_each(?))",
            (molecule_ids, group_ids),
        ).fetchall()
    evidence_runs = index["evidence_runs"]
    direct = subject.get("direct_molecule_ids", subject["molecule_ids"])
    targets = subject.get("target_record_ids", set())

    def included(molecule_id, group_id, run_key, context):
        if run_view:
            return run_key in subject["run_keys"]
        if molecule_id in direct:
            return True
        target = context.get("target_molecule_id") or context.get("target_id") or context.get("design_target_molecule_id")
        if target:
            return target in targets
        metadata_links = index["molecules"][molecule_id].get("metadata_links", set())
        return bool(metadata_links & targets)

    def make_result(kind, evidence_id, molecule_id, group_id, run_key, context, created_at, **values):
        molecule = index["molecules"][molecule_id]
        run = index["runs"].get(run_key, {})
        return {"id": _key([kind, evidence_id]), "kind": kind, "molecule_id": molecule_id,
                "group_id": group_id, "name": molecule["name"], "run_key": run_key,
                "run_id": run.get("run_id") or context.get("run_id") or "",
                "tool": context.get("tool") or context.get("source_tool") or ", ".join(sorted(run.get("tools", set()))),
                "workstation": ", ".join(sorted(run.get("workstations", set()))),
                "collection": index["groups"][group_id]["name"], "group_ids": {group_id}, "created_at": created_at,
                "evaluation_id": "", "observation_id": "", "structure_ids": [],
                "model": "", "status": "", "relation": "Molecule / complex" if molecule_id in direct else "Recorded target link",
                **values}

    result = []
    covered = set()
    used_structures = set()
    prepared_structures = [{**dict(row), "context": json.loads(row["context"])} for row in structures]
    for row in observations:
        observation = {**dict(row), "context": json.loads(row["context"])}
        run_key = evidence_runs.get(("molecule_observations", row["id"]))
        if not included(row["molecule_id"], row["group_id"], run_key, observation["context"]):
            continue
        evaluation = evaluation_descriptor(observation)
        matches = [structure for structure in prepared_structures if structure["molecule_id"] == row["molecule_id"]
                   and (structure["id"] == evaluation["representation_id"] if evaluation["representation_id"] else
                        structure["group_id"] == row["group_id"] and evaluation["run_id"]
                        and structure["context"].get("run_id") == evaluation["run_id"]
                        and (not evaluation["candidate_id"] or structure["context"].get("candidate_id") == evaluation["candidate_id"])
                        and (not evaluation["model_id"] or str(structure["context"].get("model_id", "")) == evaluation["model_id"]))]
        structure_ids = [structure["id"] for structure in matches]
        used_structures.update(structure_ids)
        covered.add((row["molecule_id"], row["group_id"], run_key))
        result.append(make_result("Evaluation", row["id"], row["molecule_id"], row["group_id"], run_key,
            observation["context"], row["created_at"], evaluation_id=evaluation["id"], observation_id=row["id"],
            structure_ids=structure_ids, tool=evaluation["evaluator"],
            model=evaluation["model_id"] or evaluation["pose_id"], candidate=evaluation["candidate_id"], status=evaluation["status"],
            kind_label=evaluation["stage"] or "Evaluation"))
    for row in prepared_structures:
        run_key = evidence_runs.get(("molecule_structures", row["id"]))
        if row["id"] in used_structures or not included(row["molecule_id"], row["group_id"], run_key, row["context"]):
            continue
        covered.add((row["molecule_id"], row["group_id"], run_key))
        result.append(make_result("Structure", row["id"], row["molecule_id"], row["group_id"], run_key,
                                 row["context"], row["created_at"], structure_ids=[row["id"]], kind_label="Structure",
                                 model=str(row["context"].get("model_id", ""))))
    for row in sources:
        context = json.loads(row["provenance"])
        run_key = evidence_runs.get(("molecule_sources", row["source_key"]))
        pair = (row["molecule_id"], row["group_id"], run_key)
        if pair in covered or not included(row["molecule_id"], row["group_id"], run_key, context):
            continue
        covered.add(pair)
        group = index["groups"][row["group_id"]]
        result.append(make_result("Molecule", row["source_key"], row["molecule_id"], row["group_id"], run_key,
                                 context, group.get("updated_at", ""), kind_label="Design / output" if run_key else "Imported molecule"))
    covered_members = {pair[:2] for pair in covered}
    for molecule_id in sorted(subject["molecule_ids"]):
        molecule = index["molecules"][molecule_id]
        for group_id in sorted(molecule["group_ids"] & subject["group_ids"]):
            if (molecule_id, group_id) in covered_members:
                continue
            run_key = evidence_runs.get(("membership", molecule_id, group_id))
            if included(molecule_id, group_id, run_key, {}):
                result.append(make_result("Molecule", [molecule_id, group_id], molecule_id, group_id, run_key,
                    {}, molecule["created_at"], kind_label="Design / output" if run_key else "Saved molecule"))
    evidence_pairs = {(row["molecule_id"], row["run_key"]) for row in result if row["kind"] != "Molecule" and row["run_key"]}
    combined = {}
    for row in sorted(result, key=lambda row: row["id"]):
        if row["kind"] != "Molecule":
            combined[row["id"]] = row
            continue
        pair = (row["molecule_id"], row["run_key"])
        if pair in evidence_pairs:
            continue
        row["id"] = _key(["molecule-result", *pair])
        existing = combined.setdefault(row["id"], row)
        existing["group_ids"].update(row["group_ids"])
        existing["collection"] = ", ".join(sorted(index["groups"][group_id]["name"] for group_id in existing["group_ids"]))
    result = list(combined.values())
    result.sort(key=lambda row: (row["name"].casefold(), row["id"]))
    result.sort(key=lambda row: row["created_at"], reverse=True)
    return result


def load_result_record(database, user, result, *, project_id="", scope="all", hide_smoke_test=True):
    """Hydrate only the selected original result, with no unrelated latest evaluation."""
    records = database.list_molecules(user, molecule_ids=[result["molecule_id"]], group_id=result["group_id"],
        project_id=project_id, scope=scope, collection="all", hide_smoke_test=hide_smoke_test, limit=1,
        evaluation_filters={"evaluation_id": result["evaluation_id"]} if result["evaluation_id"] else {})["items"]
    if not records:
        return None
    record = records[0]
    if result["evaluation_id"]:
        selected = record.get("selected_evaluation") or {}
        if selected.get("id") != result["evaluation_id"]:
            return None
        record["observations"] = [row for row in record["observations"] if row["id"] == result["observation_id"]]
    else:
        record["selected_evaluation"] = None
        record["observations"] = []
        record["scores"] = {}
        record["score_sources"] = {}
    record["structures"] = [row for row in record["representations"] if row["id"] in result["structure_ids"]]
    record["representations"] = list(record["structures"])
    return record
