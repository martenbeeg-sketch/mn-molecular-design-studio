"""Scientific navigation over saved results, independent of SQL batch folders."""
from collections import defaultdict
from datetime import datetime
import re

from mn_molecular_design_studio.data_browser import (
    _identifier, branch_nodes, folder_scopes, listed_results,
)
from mn_molecular_design_studio.database_browser_labels import target_construct_label, tool_label


def _reference_families(index, scopes):
    targets = [scope for scope in scopes.values() if scope["parent"] == "targets" and scope["subject"].get("type") == "protein"
               and scope["subject"].get("sequence") and not scope["subject"].get("target_pdb_id")
               and not scope["subject"].get("without_pdb_reference")]
    anchors = {}
    for scope in targets:
        references = {group_id for identifier in scope["subject"].get("target_record_ids", set())
                      for group_id in index["molecules"][identifier]["group_ids"]
                      if index["groups"][group_id].get("collection_kind") == "reference"}
        if references:
            anchors[scope["id"]] = {(index["groups"][group]["owner_id"], index["groups"][group]["project_id"]) for group in references}
    for child in targets:
        if child["id"] in anchors:
            continue
        groups = {group_id for identifier in child["subject"].get("target_record_ids", set())
                  for group_id in index["molecules"][identifier]["group_ids"]}
        context = {(index["groups"][group]["owner_id"], index["groups"][group]["project_id"]) for group in groups}
        sequence = child["subject"]["sequence"]
        parents = [parent for parent in targets if parent["id"] in anchors and anchors[parent["id"]] & context
                   and len(sequence) < len(parent["subject"]["sequence"]) and sequence in parent["subject"]["sequence"]
                   and parent["subject"]["sequence"].find(sequence) == parent["subject"]["sequence"].rfind(sequence)]
        if len(parents) != 1:
            continue
        parent = parents[0]
        parent.setdefault("reference_family", []).append(child["label"])
        parent["search"] += " " + child.get("search", "")
        subject = {**parent["subject"]}
        for field in ("group_ids", "molecule_ids", "direct_molecule_ids", "target_record_ids", "target_run_keys", "run_keys"):
            subject[field] = subject.get(field, set()) | child["subject"].get(field, set())
        parent["subject"] = subject
        del scopes[child["id"]]
    families = {}
    for scope in targets:
        if not scope.get("reference_family"):
            continue
        subject = scope["subject"]
        reference_id = subject["id"]
        subject["id"] = _identifier("reference-family", reference_id, sorted(subject["target_record_ids"]))
        subject.pop("sequence", None)
        subject.pop("identity", None)
        families[reference_id] = subject
    index["navigation_families"] = families


def is_dataset(group):
    origin = group.get("origin") or {}
    return group.get("collection_kind") in {"curated", "reference"} and (
        group.get("import_method") == "user_uploads" and not origin.get("registration_batch")
        and not origin.get("saved_selection") and not (origin.get("task_group") == "design" and origin.get("run_id")))


def scientific_scopes(index, projects):
    scopes = folder_scopes(index, projects)
    scopes["collections"]["label"] = "Libraries & datasets"
    scopes["ligands"]["label"] = "Compounds"
    for identifier, scope in list(scopes.items()):
        if scope.get("route") == "project":
            scope["route"] = "scientific-project"
            scope["search"] += " " + " ".join(sorted({name for molecule in index["molecules"].values()
                if molecule["group_ids"] & scope["context_groups"] for name in molecule["names"]}))
        elif scope.get("route") == "collection" and not is_dataset(index["groups"][scope["entity"]["group_id"]]):
            del scopes[identifier]
        elif scope["parent"] in scopes and scopes[scope["parent"]]["parent"] == "targets" and (
            scopes[scope["parent"]]["subject"].get("target_pdb_id")
        ):
            del scopes[identifier]
    _reference_families(index, scopes)
    for scope in list(scopes.values()):
        if scope.get("route") == "scientific-project":
            expand_project(scopes, index, scope)
            scope["loaded"] = True
        elif not scope["parent"]:
            scope["loaded"] = True
    return scopes


def target_constructs(index, outcomes):
    references = {}
    for outcome in outcomes:
        identifier = outcome.get("target_id")
        molecule = index["molecules"].get(identifier, {})
        if not any(component.get("type") == "protein" for component in molecule.get("components", [])):
            continue
        labels = references.setdefault(identifier, set())
        if outcome.get("target_reference_label"):
            labels.add(outcome["target_reference_label"])
    labels = {identifier: target_construct_label(index, identifier, reference_label=" / ".join(sorted(names)))
              for identifier, names in references.items()}
    duplicates = defaultdict(list)
    for identifier, label in labels.items():
        duplicates[label].append(identifier)
    for label, matches in duplicates.items():
        if len(matches) > 1:
            for number, identifier in enumerate(sorted(matches), 1):
                labels[identifier] = label + " · variant " + str(number)
    return dict(sorted(labels.items(), key=lambda entry: (entry[1].casefold(), entry[0])))


def expand_project(scopes, index, scope):
    from mn_molecular_design_studio.entity_names import target_names

    context = scope["context_groups"]
    datasets = [entry for entry in scopes.values() if entry["parent"] == "collections"
                and entry["subject"]["group_ids"] & context]
    covered = set().union(*(entry["subject"]["molecule_ids"] for entry in datasets)) if datasets else set()
    subjects = []
    seen = set()
    for category in ("targets", "ligands"):
        for entry in list(scopes.values()):
            if entry["parent"] != category or entry["subject"]["id"] in seen:
                continue
            if category == "ligands" and entry["subject"]["direct_molecule_ids"] & covered:
                continue
            if not entry["subject"]["group_ids"] & context:
                continue
            seen.add(entry["subject"]["id"])
            subjects.append(entry)
    originals = list(scopes.values())

    def clone(entry, parent):
        subject = entry["subject"]
        groups = subject["group_ids"] & context
        molecules = {identifier for identifier in subject["molecule_ids"] if index["molecules"][identifier]["group_ids"] & groups}
        if not groups or not molecules:
            return
        key = parent + ":" + entry["id"]
        label = entry["label"]
        entity = entry.get("entity") or {}
        if entity.get("kind") == "target":
            names = target_names(index, entity["subject_id"], groups)
            if names:
                label = " / ".join(sorted({name["name"] for name in names}))
                pdb = subject.get("target_pdb_id")
                if pdb and pdb not in label:
                    label += " · " + pdb
        scopes[key] = {**entry, "id": key, "label": label, "parent": parent, "context_groups": context,
            "subject": {**subject, "id": key, "group_ids": groups, "molecule_ids": molecules,
                        "direct_molecule_ids": subject.get("direct_molecule_ids", molecules) & molecules}, "loaded": False}
        for child in originals:
            if child["parent"] == entry["id"]:
                clone(child, key)

    for entry in [*subjects, *datasets]:
        clone(entry, scope["id"])


def _experiment_key(index, row, *, design=False):
    group = index["groups"][row["group_id"]]
    origin = group.get("origin") or {}
    run = index["runs"].get(row["run_key"], {})
    scope = (group.get("owner_id"), group.get("project_id"))
    logical = run.get("logical_run_ids", set())
    if len(logical) == 1:
        return (*scope, "logical", next(iter(logical)))
    if origin.get("run_id") == row.get("run_id") or design:
        if isinstance(origin.get("studio_run_id"), str) and origin["studio_run_id"].strip():
            return (*scope, "logical", origin["studio_run_id"])
        if origin.get("campaign_id") and origin.get("organizer_node_id"):
            return (*scope, "campaign", origin["organizer_node_id"], origin["campaign_id"])
    if design:
        member = index["molecules"].get(row["molecule_id"], {}).get("member_receipts", {}).get(row["group_id"], {})
        campaign_id = (member.get("distributed_campaign") or {}).get("campaign_id")
        if isinstance(campaign_id, str) and campaign_id.strip():
            return (*scope, "distributed", campaign_id)
    references = run.get("campaign_references", set())
    if len(references) == 1:
        return (*scope, "parent", *next(iter(references)))
    return (*scope, "run", row["run_key"] or row["group_id"])


def _experiment_label(outcomes, originals):
    if all(outcome["outcome_kind"] == "Design" for outcome in outcomes):
        tools = sorted({tool_label(tool) for outcome in outcomes for tool in outcome["origin_tools"]})
    else:
        tools = sorted({tool_label(row["tool"]) for row in originals})
    method = tools[0] if len(tools) == 1 else outcomes[0]["activity"] + " study"
    dates = sorted({row["run_id"][:8] for row in originals if re.match(r"^\d{8}[-_]\d{6}", row["run_id"])})
    try:
        parsed = datetime.strptime(dates[0], "%Y%m%d") if dates else None
        date = str(parsed.day) + " " + parsed.strftime("%b %Y") if parsed else ""
    except ValueError:
        date = "Date not recorded"
    return " · ".join(value for value in (method, date) if value)


def _experiments(index, scope, members, activity):
    grouped = {}
    for outcome in members:
        for original in outcome["originals"]:
            identity = _experiment_key(index, original, design=outcome["outcome_kind"] == "Design")
            bucket = grouped.setdefault(identity, {"outcomes": {}, "originals": {}})
            bucket["outcomes"][outcome["id"]] = outcome
            bucket["originals"][original["id"]] = original
    entries = []
    for identity, bucket in grouped.items():
        outcomes = list(bucket["outcomes"].values())
        originals = list(bucket["originals"].values())
        entries.append({"id": scope["id"] + ":activity:" + _identifier(activity) + ":experiment:" + _identifier(identity),
            "label": _experiment_label(outcomes, originals), "count": len(outcomes), "children": [],
            "selection": {"activity": activity, "original_ids": set(bucket["originals"]), "experiment": True}})
    names = defaultdict(list)
    for entry in entries:
        names[entry["label"]].append(entry)
    for peers in names.values():
        if len(peers) > 1:
            for number, entry in enumerate(sorted(peers, key=lambda entry: entry["id"]), 1):
                entry["label"] += " · experiment " + str(number)
    return sorted(entries, key=lambda entry: (entry["label"].casefold(), entry["id"]))


def _docking_methods(index, scope, members, identifier, prefix):
    methods = defaultdict(set)
    for outcome in members:
        for original in outcome["originals"]:
            methods[tool_label(original["tool"])].add(original["id"])
    nodes, selections = [], {}
    for method, original_ids in sorted(methods.items()):
        key = identifier + ":tool:" + _identifier(method)
        selected = [entry["outcome"] for entry in listed_results(members, selection={"original_ids": original_ids})]
        experiments = _experiments(index, {**scope, "id": key}, selected, "Docking")
        if len(experiments) == 1:
            selection = experiments[0].pop("selection")
            children = []
        else:
            selection = {"activity": "Docking", "original_ids": original_ids}
            children = experiments
            for experiment in children:
                experiment["label"] = experiment["label"].removeprefix(method + " · ")
                if experiment["label"].startswith("experiment "):
                    experiment["label"] = experiment["label"].capitalize()
                selections[experiment["id"]] = {**experiment.pop("selection"), "path": [*prefix, method, experiment["label"]]}
        nodes.append({"id": key, "label": method, "count": len(selected), "children": children})
        selections[key] = {**selection, "path": [*prefix, method]}
    return nodes, selections


def scientific_branches(index, scope, results):
    if scope.get("route") == "scientific-project":
        identifiers = scope.get("unassigned_outcomes", set())
        if identifiers:
            key = scope["id"] + ":unassigned"
            return [{"id": key, "label": "Unassigned results", "count": len(identifiers), "children": []}], {
                key: {"outcome_ids": identifiers, "path": ["Unassigned results"]}}
        return [], {}
    if scope["id"] in {"projects", "targets", "ligands", "collections"}:
        return [], {}
    if scope.get("route") == "collection":
        return [], {}
    activities = sorted({row["activity"] for row in results} - {"Imported molecules", "Reference molecules"})
    nodes, selections = [], {}
    for activity in activities:
        members = [row for row in results if row["activity"] == activity]
        identifier = scope["id"] + ":activity:" + _identifier(activity)
        prefix = [activity] if len(activities) > 1 else []
        if activity == "Molecular dynamics":
            old_nodes, old_selections = branch_nodes(scope, members)
            branches = old_nodes[0]["children"]

            def remember(node, path):
                if node["id"] in old_selections:
                    selections[node["id"]] = {**old_selections[node["id"]], "path": [*path, node["label"]]}
                for child in node.get("children", []):
                    remember(child, [*path, node["label"]])

            for branch in branches:
                remember(branch, prefix)
        elif activity == "Docking":
            branches, method_selections = _docking_methods(index, scope, members, identifier, prefix)
            selections.update(method_selections)
        else:
            branches = _experiments(index, scope, members, activity)
            for branch in branches:
                selections[branch["id"]] = {**branch.pop("selection"), "path": [*prefix, branch["label"]]}
        if len(activities) == 1:
            nodes.extend(branches)
        else:
            nodes.append({"id": identifier, "label": activity, "count": len(members), "children": branches})
            selections[identifier] = {"activity": activity, "path": [activity]}
    return nodes, selections


def scope_counts(scopes, index, scope, results):
    from mn_molecular_design_studio.database_browser import _component_identity, includes_record

    scope["result_count"] = len(results)
    covered = set()
    for child in scopes.values():
        if child["parent"] != scope["id"]:
            continue
        if child.get("route") == "collection":
            child["member_count"] = len(child["subject"]["molecule_ids"])
            identities = {_component_identity(component) for identifier in child["subject"]["molecule_ids"]
                          for component in index["molecules"][identifier]["components"]}
            covered.update(outcome["id"] for outcome in results if any(_component_identity(component) in identities
                for original in outcome["originals"] for component in index["molecules"][original["molecule_id"]]["components"]))
            continue
        matches = {outcome["id"] for outcome in results if any(original["group_id"] in child["subject"]["group_ids"] and includes_record(
            index, child["subject"], original["molecule_id"], original["group_id"], original["run_key"],
            {"target_molecule_id": original.get("receipt", {}).get("target_id") or original.get("evaluation", {}).get("target_id"),
             "target_pdb_id": original.get("target_pdb_id"), **original.get("md", {})})
            for original in outcome["originals"])}
        child["result_count"] = len(matches)
        covered.update(matches)
    if scope.get("route") == "scientific-project":
        scope["unassigned_outcomes"] = {outcome["id"] for outcome in results} - covered


def source_jobs(index, outcomes, selection):
    originals = {row["id"]: row for outcome in outcomes for row in outcome["originals"]
                 if row["id"] in selection.get("original_ids", set())}
    jobs = {}
    for row in originals.values():
        group = index["groups"][row["group_id"]]
        origin = group.get("origin") or {}
        run = index["runs"].get(row["run_key"], {})
        run_id = origin.get("run_id") if origin.get("studio_run_id") else row["run_id"]
        key = (run_id, row["group_id"])
        job = jobs.setdefault(key, {"Run": run_id or "No job recorded",
            "Workstation": origin.get("workstation") or row.get("workstation") or ", ".join(sorted(run.get("workstations", set()))),
            "Source collection": group["name"], "Evidence runs": set()})
        if row["run_id"] and row["run_id"] != run_id:
            job["Evidence runs"].add(row["run_id"])
    for job in jobs.values():
        job["Evidence runs"] = ", ".join(sorted(job["Evidence runs"]))
    return sorted(jobs.values(), key=lambda row: (row["Run"], row["Workstation"], row["Source collection"]))
