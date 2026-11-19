"""Read-only folder navigation over existing SQL molecules and scientific results."""
from __future__ import annotations

from collections import defaultdict
from contextlib import contextmanager
import hashlib
import json
import re
import sqlite3

from mn_molecular_design_studio.database_browser import component_subjects, molecular_subjects, run_subjects
from mn_molecular_design_studio.database_browser_labels import run_label, tool_label, target_label
from mn_molecular_design_studio.database_outcomes import filter_outcomes, scientific_outcomes
from mn_molecular_design_studio.molecular_database import MolecularDatabase
from mn_molecular_design_studio.output_vocabulary import sql_common_metrics
from mn_molecular_design_studio.scientific_naming import design_name, meaningful_name


class ReadOnlyLibrary(MolecularDatabase):
    @contextmanager
    def connection(self):
        connection = sqlite3.connect(self.path.resolve().as_uri() + "?mode=ro", uri=True, timeout=30)
        connection.row_factory = sqlite3.Row
        connection.create_function("studio_common_metrics", 3, sql_common_metrics, deterministic=True)
        connection.execute("PRAGMA query_only=ON")
        try:
            yield connection
        finally:
            connection.close()


def saved_receptor_file(database, user, structure):
    from mn_molecular_design_studio.bound_docking import receptor_file

    class ReceptorReader:
        def list_molecules(self, user, *, molecule_ids, collection, limit):
            accessible = database._visible_ids(user, scope="all", collection=collection)
            with database.connection() as connection:
                molecules = connection.execute(
                    "SELECT DISTINCT molecule_id FROM molecule_members WHERE molecule_id IN "
                    "(SELECT value FROM json_each(?)) AND group_id IN (SELECT value FROM json_each(?)) LIMIT ?",
                    (json.dumps(molecule_ids), json.dumps(accessible), limit)).fetchall()
                items = []
                for molecule in molecules:
                    representations = [{**dict(row), "context": json.loads(row["context"])} for row in connection.execute(
                        "SELECT * FROM molecule_structures WHERE molecule_id=? "
                        "AND group_id IN (SELECT value FROM json_each(?))",
                        (molecule["molecule_id"], json.dumps(accessible)))]
                    items.append({"id": molecule["molecule_id"], "representations": representations})
                return {"items": items}

        def structure_file(self, user, structure_id):
            return database.structure_file(user, structure_id)

    return receptor_file(ReceptorReader(), user, structure)


def load_saved_record(database, user, result):
    from mn_molecular_design_studio.evaluations import evaluation_descriptor, project_observation
    from mn_molecular_design_studio.output_vocabulary import normalize_metrics

    accessible = database._visible_ids(user, scope="all", collection="all")
    if result["group_id"] not in accessible:
        return None
    with database.connection() as connection:
        connection.execute("BEGIN")
        member = connection.execute(
            "SELECT m.*,mm.name,mm.metadata,mm.components AS member_components,g.origin,g.schema_json "
            "FROM molecules m JOIN molecule_members mm ON mm.molecule_id=m.id "
            "JOIN molecule_groups g ON g.id=mm.group_id WHERE m.id=? AND mm.group_id=?",
            (result["molecule_id"], result["group_id"])).fetchone()
        if member is None:
            return None
        record = {**dict(member), "components": json.loads(member["member_components"]),
            "metadata": json.loads(member["metadata"]), "origin": json.loads(member["origin"]),
            "schema": json.loads(member["schema_json"]), "scores": {}, "score_sources": {},
            "observations": [], "selected_evaluation": None}
        parameters = (result["molecule_id"], json.dumps(accessible))
        record["groups"] = [dict(row) for row in connection.execute(
            "SELECT g.id,g.name FROM molecule_members mm JOIN molecule_groups g ON g.id=mm.group_id "
            "WHERE mm.molecule_id=? AND mm.group_id IN (SELECT value FROM json_each(?))", parameters)]
        structure_ids = list(result["structure_ids"])
        if result.get("receipt", {}).get("representation_role") == "refolded_validation" and result.get("source_structure_id"):
            structure_ids.append(result["source_structure_id"])
        record["structures"] = [{**dict(row), "context": json.loads(row["context"])} for row in connection.execute(
            "SELECT * FROM molecule_structures WHERE molecule_id=? AND group_id IN (SELECT value FROM json_each(?)) "
            "AND id IN (SELECT value FROM json_each(?)) ORDER BY created_at DESC,id",
            (*parameters, json.dumps(structure_ids)))]
        record["structures"].sort(key=lambda row: row["id"] not in result["structure_ids"])
        if result.get("receipt", {}).get("representation_role") == "refolded_validation":
            for structure in record["structures"]:
                if structure["id"] == result.get("source_structure_id") and structure["id"] not in result["structure_ids"]:
                    structure["label"] = "Generated design complex"
        record["sources"] = [{"source_key": row["source_key"], **json.loads(row["provenance"])} for row in connection.execute(
            "SELECT source_key,provenance FROM molecule_sources WHERE molecule_id=? "
            "AND group_id IN (SELECT value FROM json_each(?))", parameters)]
        if result["evaluation_id"]:
            row = connection.execute("SELECT * FROM molecule_observations WHERE id=? AND molecule_id=? AND group_id=?",
                (result["observation_id"], result["molecule_id"], result["group_id"])).fetchone()
            if row is None:
                return None
            observation = {**dict(row), "metrics": json.loads(row["metrics"]), "context": json.loads(row["context"])}
            descriptor = evaluation_descriptor(observation)
            if descriptor["id"] != result["evaluation_id"]:
                return None
            observation["normalized_output"] = normalize_metrics(observation["metrics"], row["namespace"], observation["context"])
            record["observations"] = [observation]
            record["scores"], record["score_sources"], _ = project_observation(observation)
            record["selected_evaluation"] = descriptor
        record["representations"] = list(record["structures"])
        record["structure_link_status"] = result.get("structure_link_status", "")
        record["evaluation_count"] = len(record["observations"])
        return record


def _identifier(*values):
    return hashlib.sha256(json.dumps(values, sort_keys=True).encode()).hexdigest()[:24]


def friendly(value):
    return meaningful_name(value)


def _subject(index, groups):
    molecules = {identifier for identifier, molecule in index["molecules"].items() if molecule["group_ids"] & groups}
    return {"id": _identifier("folder-scope", sorted(groups)), "group_ids": groups,
            "molecule_ids": molecules, "direct_molecule_ids": molecules}


def source_collection_scope(index, group_id):
    group = index["groups"].get(group_id)
    if group is None:
        return None
    identifier = "collection:" + group_id
    return {"id": identifier, "label": group["name"], "parent": "collections",
        "subject": _subject(index, {group_id}), "run_view": False, "route": "collection",
        "entity": {"kind": "dataset", "group_id": group_id}, "search": group["name"]}


def _scope_label(name, fallback, identifier):
    return str(name) if friendly(name) else fallback


def folder_scopes(index, projects):
    scopes = {}
    references = defaultdict(list)
    groups = set(index["groups"])
    all_subject = _subject(index, groups)
    for category, label in [("targets", "Targets"), ("ligands", "Ligands"), ("projects", "Projects"),
                            ("campaigns", "Campaigns & runs"), ("collections", "Collections")]:
        scopes[category] = {"id": category, "label": label, "parent": "", "subject": all_subject, "run_view": False}
    for entry in component_subjects(index, role="targets"):
        category = "targets" if entry["role"] == "Target" else "ligands"
        if entry["type"] == "small_molecule":
            categories = ["ligands", *(["targets"] if entry["role"] == "Target" else [])]
        else:
            categories = [category] if category == "targets" else []
        for category in categories:
            identifier = category + ":" + entry["id"]
            fallback = "Protein target · " + str(len(entry["sequence"])) + " aa" if entry["type"] == "protein" else "Unnamed ligand"
            scope = {"id": identifier, "label": _scope_label(entry["name"], fallback, entry["id"]),
                "parent": category, "subject": entry, "run_view": False,
                "entity": {"kind": "target", "subject_id": entry["id"]},
                "search": " ".join([*entry["names"], entry.get("sequence", ""),
                    *(value for label in entry.get("canonical_names", []) for value in label.get("properties", {}).values())])}
            pdbs = set().union(*(index.get("target_references", {}).get(target, set()) for target in entry["target_record_ids"]))
            if category == "targets" and pdbs:
                for pdb in sorted(pdbs):
                    references[pdb].append(scope)
                if entry["target_record_ids"] & index.get("unreferenced_targets", set()):
                    scopes[identifier] = {**scope, "label": scope["label"] + " · PDB unspecified",
                        "subject": {**entry, "without_pdb_reference": True}}
            else:
                scopes[identifier] = scope
    for pdb, constructs in sorted(references.items()):
        receptor = " / ".join(sorted(index.get("target_reference_names", {}).get(pdb, set()))) or "Protein target"
        identifier = "targets:pdb:" + pdb
        subject = {**constructs[0]["subject"], "id": identifier, "target_pdb_id": pdb}
        for field in ("molecule_ids", "direct_molecule_ids", "group_ids", "run_keys", "target_record_ids", "target_run_keys"):
            subject[field] = set().union(*(row["subject"].get(field, set()) for row in constructs))
        md_context = index.get("pdb_md_contexts", {}).get(pdb, {})
        for field in ("molecule_ids", "group_ids", "run_keys"):
            subject[field].update(md_context.get(field, set()))
        subject["pdb_md_run_keys"] = md_context.get("run_keys", set())
        if len(constructs) > 1:
            subject.pop("sequence", None)
            subject.pop("identity", None)
        label = receptor + " · " + pdb
        if len(constructs) == 1 and constructs[0]["subject"].get("canonical_names"):
            label = constructs[0]["label"] + " · " + pdb
        scopes[identifier] = {"id": identifier, "label": label, "parent": "targets",
            "subject": subject, "run_view": False, "entity": constructs[0]["entity"] if len(constructs) == 1 else None,
            "search": " ".join(row["search"] for row in constructs)}
        if len(constructs) > 1:
            lengths = defaultdict(int)
            for row in constructs:
                lengths[len(row["subject"].get("sequence", ""))] += 1
            for row in constructs:
                key = identifier + ":" + row["subject"]["id"]
                length = len(row["subject"].get("sequence", ""))
                label = "Protein construct · " + str(length) + " aa"
                if row["subject"].get("canonical_names") or lengths[length] > 1:
                    label = row["label"]
                scopes[key] = {**row, "id": key, "parent": identifier, "label": label,
                    "subject": {**row["subject"], "id": key, "target_pdb_id": pdb}}
    project_names = {project["id"]: project["name"] for project in projects}
    project_aliases = {project["id"]: project.get("name_aliases", []) for project in projects}
    by_project = defaultdict(set)
    for identifier, group in index["groups"].items():
        by_project[group.get("project_id") or ""].add(identifier)
    for project, project_groups in by_project.items():
        identifier = "project:" + (project or "personal")
        scopes[identifier] = {"id": identifier, "label": project_names.get(project, "Personal" if not project else "Unnamed project"),
            "parent": "projects", "subject": _subject(index, project_groups), "run_view": False,
            "route": "project", "context_groups": project_groups,
            "entity": {"kind": "project", "project_id": project} if project else None,
            "search": " ".join(project_aliases.get(project, []))}
    for entry in run_subjects(index, combine_campaigns=True):
        identifier = "campaign:" + entry["id"]
        methods = sorted({tool_label(tool) for key in entry["run_keys"] for tool in index["runs"][key]["tools"]})
        fallback = ", ".join(methods[:2]) or "Imported run"
        name = entry["name"]
        if not friendly(name) or re.match(r"^\d{8}[-_]", name):
            dates = sorted(index["runs"][key]["run_id"] for key in entry["run_keys"] if index["runs"][key]["run_id"])
            name = fallback + " · " + run_label(dates[0] if dates else "")
        scopes[identifier] = {"id": identifier, "label": name, "parent": "campaigns", "subject": entry,
            "run_view": True, "search": " ".join(index["runs"][key]["run_id"] for key in entry["run_keys"])}
    collection_names = defaultdict(int)
    collection_search = defaultdict(set)
    for group in index["groups"].values():
        if group.get("collection_kind") in {"curated", "reference", "generated"}:
            collection_names[group["name"].casefold()] += 1
    for molecule in index["molecules"].values():
        terms = {*molecule["names"], *(component.get("sequence") or component.get("smiles") or component.get("ccd") or ""
                  for component in molecule["components"])}
        for group_id in molecule["group_ids"]:
            collection_search[group_id].update(terms)
    for identifier, group in index["groups"].items():
        if group.get("collection_kind") not in {"curated", "reference", "generated"}:
            continue
        key = "collection:" + identifier
        label = group["name"]
        if collection_names[label.casefold()] > 1:
            project = group.get("project_id") or ""
            label += " · " + project_names.get(project, "Personal" if not project else "Unnamed project")
        scopes[key] = {"id": key, "label": label, "parent": "collections",
            "subject": _subject(index, {identifier}), "run_view": False, "route": "collection",
            "entity": {"kind": "dataset", "group_id": identifier},
            "search": " ".join([*group.get("name_aliases", []), *group.get("entity_properties", {}).values(),
                                *sorted(collection_search[identifier])])}
    return scopes


def scope_path(scopes, identifier):
    path = []
    while identifier in scopes:
        path.append(identifier)
        identifier = scopes[identifier]["parent"]
    return list(reversed(path))


def _collection_molecule_label(index, scope, subject, linked):
    from mn_molecular_design_studio.entity_names import target_names

    canonical = target_names(index, linked["id"], scope["subject"]["group_ids"])
    if canonical:
        return " / ".join(sorted({label["name"] for label in canonical}))
    names = {name for identifier in subject["molecule_ids"] & scope["subject"]["molecule_ids"]
             for group, name in index["molecules"][identifier].get("member_names", {}).items()
             if group in scope["subject"]["group_ids"] and friendly(name)}
    if names:
        return min(names, key=lambda value: (value.casefold(), value))
    return linked["name"] if friendly(linked["name"]) else "Unnamed " + subject["type"].replace("_", " ")


def expand_scope(scopes, index, identifier):
    scope = scopes[identifier]
    if scope.get("loaded"):
        return
    scope["loaded"] = True
    route = scope.get("route")
    if route == "scientific-project":
        from mn_molecular_design_studio.scientific_browser_tree import expand_project
        expand_project(scopes, index, scope)
        return
    if route == "project":
        for kind, label in [("targets", "Targets & ligands"), ("campaigns", "Campaigns & runs"), ("collections", "Collections")]:
            key = identifier + ":" + kind
            scopes[key] = {**scope, "id": key, "parent": identifier, "label": label, "route": "project-" + kind, "loaded": False, "entity": None}
        return
    if route and route.startswith("project-"):
        category = route.removeprefix("project-")
        candidates = [entry for entry in scopes.values() if entry["parent"] == category or (
            category == "targets" and entry["parent"] == "ligands")]
    elif route == "collection":
        if "folder_molecular_subjects" not in index:
            index["folder_molecular_subjects"] = molecular_subjects(index)
            index["folder_component_subjects"] = component_subjects(index)
        subjects = index["folder_molecular_subjects"]
        components_by_identity = {(entry["identity"],): entry for entry in index["folder_component_subjects"]}
        candidates = []
        for subject in subjects:
            if not subject["molecule_ids"] & scope["subject"]["molecule_ids"]:
                continue
            linked = components_by_identity.get(subject["identity"], subject)
            linked = index.get("navigation_families", {}).get(linked["id"], linked)
            candidates.append({"id": "molecule:" + subject["id"], "label": _collection_molecule_label(index, scope, subject, linked),
                "subject": linked, "run_view": False, "search": " ".join([*subject["names"], subject["sequence"]]), "route": "molecule",
                "entity": {"kind": "target", "subject_id": linked["id"]} if linked.get("role") in {"Target", "Ligand"} else None})
    else:
        return
    for entry in candidates:
        subject = entry["subject"]
        groups = subject["group_ids"] & scope.get("context_groups", subject["group_ids"])
        molecules = {key for key in subject["molecule_ids"] if index["molecules"][key]["group_ids"] & groups}
        if not groups or not molecules:
            continue
        key = identifier + ":" + entry["id"]
        scopes[key] = {**entry, "id": key, "parent": identifier, "loaded": False,
            "subject": {**subject, "id": key, "molecule_ids": molecules, "group_ids": groups,
                        "direct_molecule_ids": subject.get("direct_molecule_ids", molecules) & molecules},
            **({"context_groups": scope["context_groups"]} if "context_groups" in scope else {})}


def activity(index, result):
    kind = result.get("outcome_kind") or ""
    if kind == "Design":
        components = index["molecules"][result["molecule_id"]]["components"]
        types = {component["type"] for component in components}
        if types == {"small_molecule"}:
            return "Small-molecule design"
        if types & {"dna", "rna"} and "protein" not in types:
            return "Nucleic-acid design"
        return "Protein design"
    return {"Binding/docking": "Docking", "Cofolding": "Cofolding", "Binding affinity": "Binding affinity", "Structure prediction": "Structure prediction",
            "Dynamics": "Molecular dynamics", "Imported molecule": "Imported molecules",
            "Reference molecule": "Reference molecules"}.get(kind, kind or "Saved structures")


def saved_results(index, results, *, include_references=False):
    prepared = []
    references = []
    for result in results:
        group = index["groups"][result["group_id"]]
        is_input = result.get("kind_label") == "Input / reference" or (
            "molecule_observations", result.get("observation_id")) in index.get("evidence_inputs", set())
        if is_input:
            if include_references and group.get("collection_kind") in {"curated", "reference"}:
                references.append(result)
            continue
        run = index["runs"].get(result["run_key"], {})
        legacy = str(result.get("run_id") or "").startswith("mnligand-")
        role = ("docking" if run.get("task_group") == "docking" else "independent_prediction" if (
            run.get("task_group") == "refolding" and result.get("evaluation", {}).get("stage") == "complex_refolding") else "")
        receipt = result.get("receipt", {})
        if legacy and not receipt.get("design_refold") and (role or receipt.get("role") in {
            "docking", "independent_prediction", "binding_affinity"}):
            result = {**result, "receipt": {**receipt, "role": receipt.get("role") or role}, "legacy_import": True}
        md_evaluation = result["kind"] == "Evaluation" and (
            result["tool"].startswith("md.") or index["runs"].get(result["run_key"], {}).get("task_group") in {"ligand-md", "md-analysis"})
        if group.get("collection_kind") in {"inputs", "intermediate"} and not (
            md_evaluation or result.get("legacy_import") or receipt.get("saved_selection") or receipt.get("role") in {"post_analysis", "design_validation"}
        ):
            continue
        prepared.append(result)
    results = prepared
    projected = scientific_outcomes(index, results)
    for result in projected:
        if result["outcome_kind"] == "Structure prediction" and any(row.get("legacy_import") and (
            row.get("evaluation", {}).get("stage") == "complex_refolding") for row in result["originals"]):
            result.update(outcome_kind="Cofolding", kind_label="Cofolding")
    outcomes = filter_outcomes(projected)
    retained = {row["id"] for row in outcomes}
    for result in projected:
        if result["id"] not in retained and result["outcome_kind"] in {"Binding/docking", "Cofolding", "Structure prediction", "Binding affinity"} and (
            any(row.get("legacy_import") for row in result["originals"]) and result["eligibility"] == "Needs review"):
            outcomes.append({**result, "eligibility": "Imported", "acceptance_level": "Imported calculation; no filter pass asserted"})
    md = defaultdict(list)
    from mn_molecular_design_studio.md_browser import production_result, study_identity, study_label
    for result in results:
        if result["tool"].startswith("md.") or result.get("evaluation", {}).get("evaluator", "").startswith("md.") or (
            index["runs"].get(result["run_key"], {}).get("task_group") in {"ligand-md", "md-analysis"}
        ):
            md[study_identity(index, result)].append(result)
    for key, originals in md.items():
        originals.sort(key=lambda row: (not production_result(row), not bool(row["structure_ids"]), row["id"]))
        primary = originals[0]
        target_id = primary.get("receipt", {}).get("target_id") or primary.get("evaluation", {}).get("target_id") or ""
        outcomes.append({**primary, "id": _identifier("saved-md", key), "primary": primary,
            "originals": originals, "outcome_kind": "Dynamics", "kind_label": "Dynamics", "eligibility": "Recorded",
            "origin_tools": sorted({row["tool"] for row in originals}), "source_runs": len({row["run_key"] for row in originals}), "target_name": target_label(index, target_id),
            "target_id": target_id, "acceptance_level": "Recorded dynamics; no design pass asserted",
            "md_label": study_label(primary),
            "structure_count": len({value for row in originals for value in row["structure_ids"]})})
    represented = {row["molecule_id"] for outcome in outcomes for row in outcome["originals"]}
    members = defaultdict(list)
    for result in [*results, *references]:
        group = index["groups"][result["group_id"]]
        if result["molecule_id"] in represented or group.get("collection_kind") not in {"curated", "reference"}:
            continue
        if result["kind"] not in {"Molecule", "Structure"}:
            continue
        reference = include_references and (group.get("collection_kind") == "reference" or result["kind_label"] == "Input / reference")
        if reference or group["import_method"] == "user_uploads":
            members[result["molecule_id"]].append((result, reference))
    for entries in members.values():
        reference = all(entry[1] for entry in entries)
        originals = [entry[0] for entry in entries]
        originals.sort(key=lambda row: (not bool(row["structure_ids"]), row["id"]))
        primary = originals[0]
        outcomes.append({**primary, "primary": primary, "originals": originals,
            "outcome_kind": "Reference molecule" if reference else "Imported molecule",
            "eligibility": "Reference" if reference else "Imported", "origin_tools": [], "source_runs": 0, "target_name": "", "target_id": "",
            "acceptance_level": "Saved reference; no design pass asserted" if reference else "User import",
            "structure_count": len({value for row in originals for value in row["structure_ids"]})})
    for result in outcomes:
        result["activity"] = activity(index, result)
        pdb = result.get("target_pdb_id")
        if pdb:
            target = index["molecules"].get(result.get("target_id"), {})
            length = sum(len(component.get("sequence", "")) for component in target.get("components", []))
            result["target_reference_label"] = " · ".join(value for value in (
                " / ".join(sorted(index.get("target_reference_names", {}).get(pdb, set()))), pdb,
                str(length) + " aa" if length else "") if value)
        molecule = index["molecules"][result["molecule_id"]]
        naming = molecule.get("member_naming", {}).get(result["group_id"], {})
        result["display_name"] = naming.get("display_name") or ""
        if result["outcome_kind"] == "Reference molecule":
            from mn_molecular_design_studio.entity_names import component_id, target_names
            labels = target_names(index, component_id(molecule["components"][0]), {result["group_id"]}) if len(molecule["components"]) == 1 else []
            result["display_name"] = " / ".join(sorted({label["name"] for label in labels})) if labels else result["name"]
            if not friendly(result["display_name"]):
                result["display_name"] = " / ".join(component["type"].replace("_", " ").title() + (
                    " · " + str(len(component["sequence"])) + (" aa" if component["type"] == "protein" else " nt") if component.get("sequence") else "")
                    for component in molecule["components"])
        imported_names = sorted(name for name in molecule.get("import_names", []) if meaningful_name(name))
        if not result["display_name"] and result["outcome_kind"] != "Design" and (
            len(molecule["components"]) == 1 and molecule["components"][0]["type"] == "small_molecule" and len(imported_names) == 1 and
            re.fullmatch(r"HY-\d+[A-Za-z0-9-]*", str(result["name"]))
        ):
            result["display_name"] = imported_names[0]
        result["display_run_id"] = naming.get("run_id") or result.get("run_id") or ""
        result["design_created_at"] = naming.get("created_at") or ""
        result["generating_workstation"] = naming.get("generating_workstation") or ""
        result["creation_key"] = naming.get("creation_key") or ""
        if result["outcome_kind"] == "Design" and naming.get("tool"):
            result["origin_tools"] = [naming["tool"]]
        target_id = result.get("target_id")
        if target_id and target_id not in index["molecules"] and meaningful_name(target_id):
            result["target_name"] = target_id
        result["aliases"] = sorted(molecule["names"])
        result["native_candidate_id"] = naming.get("native_candidate_id") or result.get("candidate", "")
    names = defaultdict(set)
    for result in outcomes:
        if result["display_name"] and result["creation_key"]:
            names[result["display_name"]].add(result["creation_key"])
    for result in outcomes:
        if len(names.get(result["display_name"], ())) > 1:
            result["display_name"] += "_" + result["creation_key"][:10]
    from mn_molecular_design_studio.entity_properties import project_results
    return sorted(project_results(index, outcomes), key=lambda row: (row["activity"], row["name"].casefold(), row["id"]))


def result_label(result):
    if result.get("outcome_kind") == "Dynamics" and result.get("md_label"):
        return result["md_label"]
    if result.get("display_name"):
        return result["display_name"]
    name = str(result["name"] or "")
    tools = " / ".join(tool_label(value) for value in result.get("origin_tools", []))
    if friendly(name) and not re.match(r"^(?:rfd\d|bindcraft\d*|boltzgen)_\d", name, re.I):
        return name
    proposed = design_name(name, target_name=result.get("target_name"), tool=tools, run_id=result.get("run_id"),
        metadata={"native_candidate_id": result.get("native_candidate_id") or name})
    return proposed or " · ".join(value for value in [result.get("target_name"), tools,
        result.get("outcome_kind", "Result"), "Creation time not recorded"] if value)


def concise_result_label(result):
    return result_label(result)


def producing_evaluations(index, outcome, connected):
    originals = {row["id"] for row in outcome["originals"]}
    producers = {tool_label(tool) for tool in outcome["origin_tools"]}
    selected = []
    for row in connected:
        if not row["evaluation_id"]:
            continue
        receipt = row.get("receipt") or {}
        stage = row.get("evaluation", {}).get("stage") or ""
        generated = outcome["outcome_kind"] == "Design" and tool_label(row["tool"]) in producers and (
            receipt.get("role") in {"design_generation", "sequence_redesign"} or
            stage.startswith(("generation.", "sequence_design")) or
            index["runs"].get(row["run_key"], {}).get("task_group") == "design")
        if row["id"] in originals or generated and receipt.get("role") not in {"post_analysis", "design_validation"} and (
            not receipt.get("design_refold") and stage not in {"analysis", "benchmark"}
        ):
            selected.append(row)
    return selected


def admission_measurements(record):
    from mn_molecular_design_studio.output_vocabulary import field_label

    metadata = record.get("metadata") or {}
    validation = metadata.get("validation_acceptance") or {}
    naming = metadata.get("studio_naming") or {}
    expected = {str(value) for value in [metadata.get("source_candidate_id"), metadata.get("native_candidate_id"),
        metadata.get("candidate_id"), naming.get("native_candidate_id")] if value}
    rows = []
    for field, measurement in (validation.get("measurements") or {}).items():
        if not isinstance(measurement, dict) or measurement.get("status") != "measured" or measurement.get("value") is None:
            continue
        candidate = str(measurement.get("measurement_candidate_id") or "")
        if expected and candidate and candidate not in expected:
            continue
        rows.append({"Measurement": field_label(field), "Value": str(measurement["value"]) + (
            " " + str(measurement["unit"]) if measurement.get("unit") else ""),
            "Tool": tool_label(measurement.get("source_tool")), "Method": measurement.get("method") or "Not recorded"})
    return rows


def branch_nodes(scope, results, *, campaigns=()):
    nodes, selections = [], {}
    activities = sorted({row["activity"] for row in results})
    for name in activities:
        members = [row for row in results if row["activity"] == name]
        identifier = scope["id"] + ":activity:" + _identifier(name)
        node = {"id": identifier, "label": name, "count": len(members), "children": []}
        selections[identifier] = {"activity": name, "tool": ""}
        if name == "Molecular dynamics":
            for study in sorted(members, key=lambda row: (result_label(row).casefold(), row["id"])):
                key = identifier + ":study:" + study["id"]
                replicas = sorted({str(row.get("md", {}).get("replica_index")) for row in study["originals"]
                    if row.get("md", {}).get("replica_index") is not None}, key=lambda value: (len(value), value))
                replica_nodes = []
                selections[key] = {"activity": name, "md_outcome_id": study["id"]}
                for replica in replicas:
                    replica_key = key + ":replica:" + replica
                    replica_nodes.append({"id": replica_key, "label": "Replica " + replica, "children": []})
                    selections[replica_key] = {"activity": name, "md_outcome_id": study["id"], "md_replica": replica}
                node["children"].append({"id": key, "label": result_label(study), "children": replica_nodes})
            nodes.append(node)
            continue
        for tool in sorted({tool for row in members for tool in row["origin_tools"]}):
            key = identifier + ":tool:" + _identifier(tool)
            tool_members = [row for row in members if tool in row["origin_tools"]]
            campaign_nodes = []
            for campaign in campaigns:
                run_keys = campaign["subject"]["run_keys"]
                count = sum(any(original["run_key"] in run_keys for original in row["originals"]) for row in tool_members)
                if count:
                    campaign_key = key + ":" + campaign["id"]
                    campaign_nodes.append({"id": campaign_key, "label": campaign["label"], "count": count, "children": []})
                    selections[campaign_key] = {"activity": name, "tool": tool, "run_keys": run_keys}
            node["children"].append({"id": key, "label": tool_label(tool),
                "count": len(tool_members), "children": campaign_nodes})
            selections[key] = {"activity": name, "tool": tool}
        nodes.append(node)
    return nodes, selections


def navigation_tree(scopes, *, active_scope, children, query="", limit=100, root="", molecule_type=""):
    hidden = 0
    terms = query.casefold().split()
    by_parent = defaultdict(list)
    active_path = set(scope_path(scopes, active_scope))
    matches = set(active_path)
    for scope in scopes.values():
        text = " ".join([*(scopes[key]["label"] for key in scope_path(scopes, scope["id"])), scope.get("search", "")]).casefold()
        if all(term in text for term in terms):
            matches.update(scope_path(scopes, scope["id"]))
        by_parent[scope["parent"]].append(scope)

    def render(identifier):
        nonlocal hidden
        members = [scope for scope in by_parent[identifier] if scope["id"] in matches and (
            not molecule_type or not scope["subject"].get("type") or scope["subject"]["type"] == molecule_type or scope["id"] in active_path)]
        members.sort(key=lambda row: (row["label"].casefold(), row["id"]))
        visible = members[:limit]
        for active in members:
            if active["id"] in active_path and active not in visible:
                visible = [*visible[:-1], active]
        hidden += max(0, len(members) - len(visible))
        nested = [render(row["id"]) for row in visible]
        nested = [*nested, *(children if identifier == active_scope else scopes[identifier].get("result_branches", []))]
        node = {"id": identifier, "label": scopes[identifier]["label"], "children": nested,
                "lazy": identifier != active_scope and not scopes[identifier].get("loaded")}
        if members:
            node.update(count=len(members), countLabel="folders")
        if "result_count" in scopes[identifier] and scopes[identifier]["parent"]:
            node.update(count=scopes[identifier]["result_count"], countLabel="saved results")
        elif "member_count" in scopes[identifier]:
            node.update(count=scopes[identifier]["member_count"], countLabel="molecules")
        return node

    categories = [root] if root else ["targets", "ligands", "projects", "campaigns", "collections"]
    return [render(category) for category in categories], hidden


def listed_results(results, *, selection=None, query="", structures=False):
    selection = selection or {}
    terms = query.casefold().split()
    matched = []
    for result in results:
        if selection.get("outcome_ids") and result["id"] not in selection["outcome_ids"]:
            continue
        if selection.get("md_outcome_id") and result["id"] != selection["md_outcome_id"]:
            continue
        if selection.get("activity") and result["activity"] != selection["activity"]:
            continue
        if selection.get("tool") and selection["tool"] not in result["origin_tools"]:
            continue
        if selection.get("run_keys") and not any(row["run_key"] in selection["run_keys"] for row in result["originals"]):
            continue
        if selection.get("original_ids") and not any(row["id"] in selection["original_ids"] for row in result["originals"]):
            continue
        text = " ".join([result_label(result), result["name"], result.get("target_name", ""), *result.get("aliases", []), *result["origin_tools"],
                         *(str(row.get(key) or "") for row in result["originals"] for key in ("candidate", "run_id", "collection"))]).casefold()
        if not all(term in text for term in terms):
            continue
        originals = [row for row in result["originals"] if (not selection.get("original_ids") or row["id"] in selection["original_ids"]) and (
            not selection.get("tool") or row["tool"] == selection["tool"] or
            row.get("receipt", {}).get("generator_engine") == selection["tool"]) and (
            not selection.get("run_keys") or row["run_key"] in selection["run_keys"])]
        if not originals and result["outcome_kind"] == "Design" and selection.get("tool") in result["origin_tools"]:
            originals = [row for row in result["originals"] if not selection.get("run_keys") or row["run_key"] in selection["run_keys"]]
        if selection.get("tool") or selection.get("run_keys") or selection.get("original_ids"):
            if not originals:
                continue
            primary = result["primary"] if result["primary"] in originals else min(originals,
                key=lambda row: (not bool(row["structure_ids"]), not bool(row["evaluation_id"]), row["id"]))
            result = {**result, "primary": primary, "originals": originals,
                "origin_tools": ([selection["tool"]] if selection.get("tool") else
                    result["origin_tools"] if result["outcome_kind"] == "Design" else sorted({row["tool"] for row in originals})),
                "structure_count": len({value for row in originals for value in row["structure_ids"]})}
        if not structures:
            matched.append({"id": result["id"], "label": result_label(result), "outcome": result, "result": result["primary"], "structure_id": ""})
            continue
        seen = set()
        for original in originals:
            for identifier in [*original["structure_ids"], *([original["source_structure_id"]] if
                    original.get("receipt", {}).get("representation_role") == "refolded_validation" and original.get("source_structure_id") else [])]:
                if identifier in seen:
                    continue
                seen.add(identifier)
                displayed = original
                if identifier not in original["structure_ids"]:
                    displayed = {**original, "structure_ids": [identifier], "evaluation_id": "", "observation_id": "", "evaluation": {},
                        "metrics": {}, "source_structure_id": "", "source_structure_digest": "",
                        "structure_digest": original.get("source_structure_digest") or "",
                        "receipt": {**original.get("receipt", {}), "representation_role": "generated_design"}}
                matched.append({"id": _identifier(result["id"], identifier), "label": result_label(result), "outcome": result,
                    "result": displayed, "structure_id": identifier})
    return matched
