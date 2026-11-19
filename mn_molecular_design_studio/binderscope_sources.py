"""Read-only, source-specific target evidence for saved BinderScope designs."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re


SCHEMA = "studio.binderscope_target_context.v1"


def target_context(database, user, record, structure):
    metadata = record.get("metadata") or {}
    context = structure.get("context") or {}
    candidate_set = metadata.get("candidate_set") or {}
    native_input = (candidate_set.get("source_provenance") or {}).get("input") or {}
    scopes = [("member", metadata), ("structure", context),
              ("native_input.params", native_input.get("params") or {})]
    with database.connection() as connection:
        observations = connection.execute(
            "SELECT namespace,context FROM molecule_observations WHERE group_id=? AND molecule_id=? ORDER BY id",
            (structure["group_id"], record["id"])).fetchall()
    for observation in observations:
        saved = json.loads(observation["context"])
        if saved.get("representation_id") != structure["id"] or saved.get("source_kind") in {"database", "results"}:
            continue
        scopes.extend([(f"evaluation.{observation['namespace']}", saved),
                       (f"evaluation.{observation['namespace']}.settings", saved.get("settings") or {}),
                       (f"evaluation.{observation['namespace']}.candidate", (saved.get("candidate") or {}).get("raw_metadata") or {})])
    declarations, evidence = [], []
    for label, scope in scopes:
        for field in ("hotspots", "target_hotspot_residues", "hotspot_residues"):
            if field not in scope or scope[field] is None:
                continue
            value = scope[field]
            if isinstance(value, str):
                tokens = [token for token in re.split(r"[,;\s]+", value.strip()) if token]
            elif isinstance(value, list) and all(type(item) in {str, int} for item in value):
                tokens = [str(item).strip() for item in value if str(item).strip()]
            else:
                tokens = ["Unsupported hotspot declaration"]
            tokens = sorted(set(tokens))
            declarations.append(tokens)
            evidence.append({"source": label, "field": field, "value": tokens})
    unique = {tuple(tokens) for tokens in declarations}
    hotspots = {"status": "conflicting" if len(unique) > 1 else "recorded" if unique and next(iter(unique)) else
                "not_specified" if declarations else "not_recorded",
                "requested": sorted(set(token for declaration in declarations for token in declaration)), "evidence": evidence,
                "numbering": "saved_target_author" if context.get("hotspot_numbering") == "saved_target_author" else "input_target_author"}
    kinds = set()
    for label, scope in scopes:
        for settings in (scope, scope.get("validation_acceptance") or {}):
            kind = settings.get("target_type")
            if kind in {"sequence", "idr", "disordered", "unstructured"}:
                kinds.add("sequence")
            elif kind in {"structured", "structure"}:
                kinds.add("structured")
    kind = next(iter(kinds)) if len(kinds) == 1 else "unknown"
    identifiers = {str(scope[field]) for label, scope in scopes
                   for field in ("design_target_molecule_id", "target_molecule_id") if scope.get(field)}
    result = {"schema": SCHEMA, "hotspots": hotspots, "target_kind": kind,
              "reference": {"status": "unavailable", "reason": "No unique explicit SQL target association."}}
    if len(identifiers) != 1:
        return result, None
    target_id = next(iter(identifiers))
    visible = database._visible_ids(user, scope="all", collection="all")
    group = database.group(structure["group_id"], user)
    with database.connection() as connection:
        references = connection.execute(
            "SELECT s.* FROM molecule_structures s JOIN molecule_groups g ON g.id=s.group_id "
            "WHERE s.molecule_id=? AND s.group_id IN (SELECT value FROM json_each(?)) AND g.project_id=? ORDER BY s.id",
            (target_id, json.dumps(visible), str(group.get("project_id") or ""))).fetchall()
    paths = []
    for scope in [native_input.get("inputs") or {}, *[scope for label, scope in scopes]]:
        paths.extend(str(scope[field]) for field in ("target_structure", "target_pdb", "target_path") if scope.get(field))
    prefixes = {match.group(1) for path in paths
                if (match := re.search(r"(?:target-|^)([0-9a-f]{16,64})\.(?:pdb|cif|mmcif)(?:\.gz)?$", Path(path).name))}
    digests = {str(scope["target_structure_digest"]) for label, scope in scopes if scope.get("target_structure_digest")}
    eligible = []
    for row in references:
        saved = json.loads(row["context"])
        if saved.get("structure_role") != "source_reference" and saved.get("input_role") not in {"target", "receptor"}:
            continue
        digest = saved.get("portable_sha256") or saved.get("sha256") or saved.get("source_structure_digest")
        if not isinstance(digest, str) or not re.fullmatch(r"[0-9a-f]{64}", digest):
            continue
        if prefixes and not any(digest.startswith(prefix) for prefix in prefixes):
            continue
        if not prefixes and digests and digest not in digests:
            continue
        eligible.append((row, digest))
    if len({digest for row, digest in eligible}) != 1:
        result["reference"] = {"status": "unavailable", "reason": "Original input target is missing, inaccessible or has ambiguous coordinate versions."}
        return result, None
    row, digest = eligible[0]
    try:
        path = database.structure_file(user, row["id"])
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError("Original target digest changed.")
        from mn_molecular_design_studio.molecule_imports import parse_structure

        parsed = parse_structure(path.read_bytes(), path.name, include_ligands=True, require_model_selection=True)[0]
        proteins = [component for component in parsed["components"] if component.get("type") == "protein"]
        chains = {str(component.get("source_chain") or component.get("chain")) for component in proteins}
        declared = set()
        for scope in [native_input.get("params") or {}, context]:
            value = scope.get("input_target_chains") or scope.get("target_chains")
            if isinstance(value, str):
                declared.update(token for token in re.split(r"[,\s]+", value) if token)
            elif isinstance(value, list):
                declared.update(str(token) for token in value)
        selected = chains & declared if declared else chains
        if len(selected) != 1:
            raise ValueError("Original target chain cannot be selected unambiguously.")
        result["reference"] = {"status": "available", "molecule_id": target_id, "group_id": row["group_id"],
                               "structure_id": row["id"], "structure_digest": digest, "chain": next(iter(selected)),
                               "basis": "Explicit SQL target association and verified original input coordinates."}
        return result, path
    except (OSError, ValueError, PermissionError) as exc:
        result["reference"] = {"status": "unavailable", "reason": str(exc)}
        return result, None
