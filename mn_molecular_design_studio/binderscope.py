"""Versioned contracts for the independently packaged BinderScope runtime."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path


REVISION = "c61ef9ca3dfc64dbde018ae7329d7538d8c583f1"
IMAGE = "mn-binderscope:c61ef9c"
METHOD_VERSION = "studio.binderscope.saved-cohort.v4"
SELECTION_SCHEMA = "studio.binderscope_selection.v1"
COLLECTION_SCHEMA = "studio.binderscope_collection.v1"
TOPOLOGIES = ("unspecified", "linear protein", "linear peptide", "cyclic peptide")
CONFIDENCE_ENCODINGS = ("uninterpreted", "plddt_0_100")


def validate_assessment_settings(settings):
    if not isinstance(settings, dict) or set(settings) - {
        "component_format", "b_factor_interpretation", "topology_basis", "confidence_basis",
    }:
        raise ValueError("Invalid per-design BinderScope settings.")
    if settings.get("component_format") not in TOPOLOGIES:
        raise ValueError("Invalid saved binder topology.")
    if settings.get("b_factor_interpretation") not in CONFIDENCE_ENCODINGS:
        raise ValueError("Invalid saved confidence encoding.")
    for field in ("topology_basis", "confidence_basis"):
        if field in settings and (not isinstance(settings[field], str) or len(settings[field]) > 1000):
            raise ValueError("Invalid BinderScope settings evidence.")
    return dict(settings)


def infer_assessment_settings(record, structure, binder_chains):
    from mn_molecular_design_studio.complex_analysis_batch import recorded_component_format

    context = structure.get("context") or {}
    metadata = record.get("metadata") or {}
    evidence = [context]
    for observation in record.get("observations") or []:
        saved = observation.get("context") or {}
        candidate = saved.get("candidate") or {}
        declared_chains = saved.get("binder_chains") or candidate.get("binder_chains") or []
        if saved.get("representation_id") == structure["id"] and declared_chains == binder_chains:
            evidence.append(saved)
    formats, confidence = set(), set()
    aliases = {"single_chain_linear_protein": "linear protein", "linear_protein": "linear protein",
               "linear protein": "linear protein", "miniprotein": "linear protein",
               "linear_peptide": "linear peptide", "linear peptide": "linear peptide",
               "cyclic_construct": "cyclic peptide", "cyclic_peptide": "cyclic peptide", "cyclic peptide": "cyclic peptide"}
    for item in [metadata, *evidence]:
        scopes = [item, item.get("settings") or {}, (item.get("validation_acceptance") or {})]
        candidate = item.get("candidate") or {}
        scopes.append(candidate.get("raw_metadata") or {})
        for scope in scopes:
            for field in ("component_format", "construct_format", "binder_format", "modality"):
                value = scope.get(field)
                if isinstance(value, str) and value in aliases:
                    formats.add(aliases[value])
            if scope.get("cyclize_peptide") is True:
                formats.add("cyclic peptide")
    if not formats:
        recorded, basis = recorded_component_format({"context": context, "record": record}, binder_chains)
        if recorded in aliases:
            formats.add(aliases[recorded])
    else:
        basis = "Recorded binder format / generator settings."
    for item in evidence:
        for scope in (item, item.get("settings") or {}):
            encoding = scope.get("b_factor_interpretation") or scope.get("confidence_encoding")
            if encoding in {"plddt_0_100", "uninterpreted", "experimental_b_factor"}:
                confidence.add("plddt_0_100" if encoding == "plddt_0_100" else "uninterpreted")
            if type(scope.get("bfactors_are_plddt")) is bool:
                confidence.add("plddt_0_100" if scope["bfactors_are_plddt"] else "uninterpreted")
    return {
        "component_format": next(iter(formats)) if len(formats) == 1 else "unspecified",
        "b_factor_interpretation": "plddt_0_100" if confidence == {"plddt_0_100"} else "uninterpreted",
        "topology_basis": basis if len(formats) == 1 else "Conflicting saved formats." if formats else "No recorded topology; sequence length is not used.",
        "confidence_basis": "Explicit encoding recorded for this exact representation." if confidence == {"plddt_0_100"} else
            "Conflicting encoding declarations; confidence checks disabled." if len(confidence) > 1 else
            "No confirmed per-residue pLDDT encoding; native scalar scores are not substituted.",
    }


def project_selections(selection, index):
    grouped = {}
    for source in selection["sources"]:
        group = index["groups"].get(source["group_id"])
        if group is None:
            raise PermissionError("A selected design's project is no longer accessible. Reopen the selection.")
        project_id = str(group.get("project_id") or "")
        grouped.setdefault(project_id, {"schema": SELECTION_SCHEMA, "sources": []})["sources"].append(source)
    return grouped


def existing_assessments(database, user, sources):
    from mn_molecular_design_studio.evaluations import evaluation_descriptor

    visible = database._visible_ids(user, scope="all", collection="all")
    wanted = {(source["molecule_id"], source["structure_id"], source["group_id"]): source for source in sources
        if source["group_id"] in visible}
    counts = dict.fromkeys(wanted, 0)
    if not wanted:
        return counts
    with database.connection() as connection:
        projects = {row["id"]: row["project_id"] for row in connection.execute(
            "SELECT id,project_id FROM molecule_groups WHERE id IN (SELECT value FROM json_each(?))", (json.dumps(visible),))}
        for row in connection.execute(
            "SELECT id,molecule_id,group_id,namespace,run_id,"
            "json_remove(context,'$.candidate.raw_metadata.binderscope') AS context,"
            "json_extract(context,'$.candidate.raw_metadata.binderscope.source.input_structure_digest') AS report_digest "
            "FROM molecule_observations WHERE namespace='binderscope' "
            "AND group_id IN (SELECT value FROM json_each(?)) "
            "AND molecule_id IN (SELECT value FROM json_each(?))",
            (json.dumps(visible), json.dumps(sorted({key[0] for key in wanted}))),
        ):
            observation = {**dict(row), "context": json.loads(row["context"])}
            evaluation = evaluation_descriptor(observation)
            if evaluation["status"].casefold() not in {"completed", "complete", "success", "succeeded", "ok"}:
                continue
            for identity, source in wanted.items():
                if (row["molecule_id"], evaluation["representation_id"]) != identity[:2]:
                    continue
                digest = source.get("structure_digest")
                if not digest or row["report_digest"] != digest or evaluation["source_structure_digest"] != digest:
                    continue
                if evaluation["calculation_structure_digest"] and evaluation["calculation_structure_digest"] != digest:
                    continue
                if source["group_id"] not in projects or row["group_id"] not in projects or projects[row["group_id"]] != projects[source["group_id"]]:
                    continue
                if any(evaluation[field] != sorted(source.get(field) or []) for field in ("binder_chains", "target_chains")):
                    continue
                target = (source.get("chain_role_evidence") or {}).get("target_molecule_id") or (
                    source.get("target_context", {}).get("reference") or {}).get("molecule_id")
                if target and evaluation["target_id"] != target:
                    continue
                counts[identity] += 1
    return counts


def protocol_revision():
    digest = hashlib.sha256(REVISION.encode())
    for relative in ("binderscope.py", "engine/workflow_engine/tools/binderscope.py",
                     "binderscope_sources.py", "engine/workflow_engine/tools/binderscope_context.py",
                     "engine/workflow_engine/workflows/binderscope.py", "vendor/binderscope/SOURCE_MANIFEST.json"):
        path = Path(__file__).parent / relative
        digest.update(relative.encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()


def source_selector(entry):
    outcome = entry.get("outcome") or {}
    result = entry.get("result") or {}
    if outcome.get("outcome_kind") != "Design":
        return None, "Not a saved protein design"
    structure_id = str(entry.get("structure_id") or "")
    if not structure_id:
        return None, "No exact saved structure"
    evaluation = result.get("evaluation") or {}
    binder = list(result.get("binder_chains") or evaluation.get("binder_chains") or [])
    target = list(result.get("target_chains") or evaluation.get("target_chains") or [])
    context = result.get("context") or {}
    binder = binder or list(context.get("binder_chains") or [])
    target = target or list(context.get("target_chains") or [])
    if len(binder) != 1 or len(target) != 1 or set(binder) & set(target):
        return None, "Requires one explicitly recorded binder chain and one target chain"
    return {
        "kind": "database", "molecule_id": result["molecule_id"],
        "structure_id": structure_id, "group_id": result["group_id"],
        "model_id": str(result.get("model") or ""),
        "binder_chains": binder, "target_chains": target,
        "label": str(outcome.get("name") or result.get("name") or structure_id),
    }, ""


def selection_from_entries(entries):
    sources, excluded, seen = [], [], set()
    for entry in entries:
        source, reason = source_selector(entry)
        if source is None:
            excluded.append({"Name": str((entry.get("outcome") or {}).get("name") or "Result"), "Reason": reason})
            continue
        identity = (source["molecule_id"], source["structure_id"], source["group_id"])
        if identity not in seen:
            sources.append(source)
            seen.add(identity)
    return {"schema": SELECTION_SCHEMA, "sources": sources}, excluded


def prepare_browser_selection(database, user, entries):
    from mn_molecular_design_studio.data_browser import load_saved_record
    from mn_molecular_design_studio.molecule_imports import parse_structure

    prepared, excluded = [], []
    for entry in entries:
        automatic_model = not entry.get("structure_id")
        if not entry.get("structure_id"):
            structures = list(dict.fromkeys(entry.get("result", {}).get("structure_ids") or []))
            if len(structures) == 1:
                entry = {**entry, "structure_id": structures[0]}
            elif len(structures) > 1:
                excluded.append({"Name": str(entry.get("outcome", {}).get("name") or "Design"),
                    "Reason": "Multiple saved models: choose All saved structures and select exact models"})
                continue
        if (entry.get("outcome") or {}).get("outcome_kind") != "Design" or not entry.get("structure_id"):
            excluded.append({"Name": str(entry.get("outcome", {}).get("name") or "Result"), "Reason": "Not a saved design with exact coordinates"})
            continue
        try:
            record = load_saved_record(database, user, entry["result"])
            if record is None:
                raise ValueError("Saved design is no longer accessible")
            from mn_molecular_design_studio.rfd4_representations import expects_refolded_model
            if automatic_model and expects_refolded_model(record) and entry["result"].get("receipt", {}).get("representation_role") != "refolded_validation":
                raise ValueError("Only the generated RFD4 model is saved here; its refolded validation complex must be recovered before default analysis. Select a generated model explicitly only for design-stage geometry.")
            structure = next(item for item in record["structures"] if item["id"] == entry["structure_id"])
            path = database.structure_file(user, structure["id"])
            parsed = parse_structure(path.read_bytes(), path.name, include_ligands=True)
            components = parsed[0]["components"] if parsed else []
            binder, target, role_evidence = saved_chain_roles(database, user, record, structure, components)
            selected = {**entry, "result": {**entry["result"], "binder_chains": binder, "target_chains": target}}
            source, reason = source_selector(selected)
            if source is None:
                raise ValueError(reason)
            source.pop("model_id")
            from mn_molecular_design_studio.evaluations import file_digest
            source["structure_digest"] = file_digest(path)
            source["structure_label"] = str(structure.get("label") or path.name)
            source["assessment"] = infer_assessment_settings(record, structure, binder)
            source["chain_role_evidence"] = role_evidence
            from mn_molecular_design_studio.binderscope_sources import target_context
            source["target_context"], _ = target_context(database, user, record, structure)
            prepared.append(source)
        except (ValueError, OSError, PermissionError, StopIteration) as exc:
            excluded.append({"Name": str(entry.get("outcome", {}).get("name") or "Result"), "Reason": str(exc)})
    unique = {(source["molecule_id"], source["structure_id"], source["group_id"]): source for source in prepared}
    return {"schema": SELECTION_SCHEMA, "sources": list(unique.values())}, excluded


def saved_chain_roles(database, user, record, structure, parsed_components):
    context = structure.get("context") or {}
    stored = list(record.get("components") or [])
    evidence = {"method": "recorded_chain_roles"}
    if ("binder_chains" not in context and "target_chains" not in context
            and not any(item.get("type") == "protein" and item.get("role") == "target" for item in stored)):
        identifiers = {str(item[field]) for item in (context, record.get("metadata") or {}, record.get("origin") or {})
                       for field in ("design_target_molecule_id", "target_molecule_id") if item.get(field)}
        if len(identifiers) > 1:
            raise ValueError("Conflicting saved target associations; chain roles cannot be reconciled automatically.")
        if identifiers:
            identifier = next(iter(identifiers))
            accessible = database._visible_ids(user, scope="all", collection="all")
            with database.connection() as connection:
                target_record = connection.execute(
                    "SELECT components FROM molecules m WHERE id=? AND EXISTS "
                    "(SELECT 1 FROM molecule_members mm WHERE mm.molecule_id=m.id "
                    "AND mm.group_id IN (SELECT value FROM json_each(?)))",
                    (identifier, json.dumps(accessible))).fetchone()
            if target_record is None:
                raise PermissionError("The recorded target is unavailable or not accessible; its sequence cannot authorize chain matching.")
            target_components = [item for item in json.loads(target_record["components"])
                                 if item.get("type") == "protein" and item.get("sequence")]
            if len(target_components) != 1:
                raise ValueError("The recorded target must identify one protein sequence for automatic BinderScope chain matching.")
            stored.append({"type": "protein", "role": "target", "sequence": target_components[0]["sequence"]})
            evidence = {"method": "saved_binder_and_sql_target_sequence_match", "target_molecule_id": identifier}
    binder, target = recorded_chains(context, stored, parsed_components)
    sequences = {str(item.get("source_chain") or item.get("auth_asym_id") or item.get("chain") or ""):
                 str(item.get("sequence") or "") for item in parsed_components if item.get("type") == "protein"}
    evidence.update(binder_sequence_sha256=hashlib.sha256(sequences[binder[0]].encode()).hexdigest(),
                    target_sequence_sha256=hashlib.sha256(sequences[target[0]].encode()).hexdigest())
    return binder, target, evidence


def recorded_chains(context, stored_components, parsed_components):
    proteins = [(str(item.get("source_chain") or item.get("auth_asym_id") or item.get("chain") or ""),
                 str(item.get("sequence") or "")) for item in parsed_components if item.get("type") == "protein"]
    available = {chain for chain, sequence in proteins if chain and sequence}
    if "binder_chains" in context or "target_chains" in context:
        binder, target = context.get("binder_chains"), context.get("target_chains")
        if (isinstance(binder, list) and isinstance(target, list) and len(binder) == len(target) == 1
                and all(isinstance(chain, str) for chain in binder + target)
                and set(binder + target) <= available and not set(binder) & set(target)):
            return binder, target
        raise ValueError("Recorded binder/target chain selections are incomplete or incompatible")
    roles = {"binder": [], "target": []}
    for item in stored_components:
        if item.get("type") != "protein" or item.get("role") not in roles or not item.get("sequence"):
            continue
        matches = [chain for chain, sequence in proteins if sequence == item["sequence"]
                   and (not item.get("chain") or chain == item["chain"])]
        if len(matches) != 1:
            raise ValueError("Recorded sequence roles cannot be uniquely mapped onto this structure")
        roles[item["role"]].append(matches[0])
    if len(roles["binder"]) != 1 or len(roles["target"]) != 1 or set(roles["binder"]) & set(roles["target"]):
        raise ValueError("No unambiguous recorded binder and target roles; chain order is not inferred")
    return roles["binder"], roles["target"]


def settings_key(selection, settings):
    return hashlib.sha256(json.dumps([selection, settings, REVISION, METHOD_VERSION, protocol_revision()], sort_keys=True).encode()).hexdigest()


def saved_result(database, user, analysis):
    allowed = database._visible_ids(user, scope="all", collection="all")
    if analysis["group_id"] not in allowed:
        return None
    with database.connection() as connection:
        row = connection.execute("SELECT context FROM molecule_observations WHERE id=? AND molecule_id=? AND group_id=?",
            (analysis["observation_id"], analysis["molecule_id"], analysis["group_id"])).fetchone()
    if row is None:
        return None
    context = json.loads(row["context"])
    result = context.get("candidate", {}).get("raw_metadata", {}).get("binderscope")
    if not isinstance(result, dict):
        return None
    return {**result, "_studio_sql_binding": {
        "structure_id": context.get("representation_id"), "structure_digest": context.get("source_structure_digest"),
        "molecule_id": analysis["molecule_id"], "observation_id": analysis["observation_id"], "group_id": analysis["group_id"],
    }}
