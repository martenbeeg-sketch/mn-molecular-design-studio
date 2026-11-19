"""Recover exact targets from preserved MN-Ligand calculation inputs."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re

from mn_molecular_design_studio.molecule_imports import canonical_components, parse_structure
from mn_molecular_design_studio.target_registration import contents_identity


POLICY = "studio.legacy-input-target.v1"


def _read(path, cache):
    path = Path(path).resolve()
    if path not in cache:
        if path.stat().st_size > 10_000_000:
            raise ValueError("The target input exceeds the supported snapshot size.")
        content = path.read_bytes()
        if len(content) > 10_000_000:
            raise ValueError("The target input exceeds the supported snapshot size.")
        cache[path] = (content, hashlib.sha256(content).hexdigest())
    return cache[path]


def _json(path, cache):
    key = (Path(path).resolve(), "json")
    if key not in cache:
        cache[key] = json.loads(_read(path, cache)[0])
    return cache[key]


def _inside(root, path):
    path = Path(path).resolve()
    if not path.is_relative_to(root.resolve()):
        raise ValueError("The target input escapes its preserved run.")
    return path


def _selector(context):
    candidate = context.get("candidate") or {}
    row = context.get("row") or {}
    candidate_id = str(context.get("candidate_id") or candidate.get("candidate_id") or "")
    parts = candidate_id.split("::")
    attempts = {str(value) for value in (context.get("legacy_run_id"), row.get("legacy_run_id"),
        parts[0] if len(parts) > 1 else "") if value}
    compounds = {str(value) for value in (row.get("compound_id"), context.get("compound_id"),
        parts[1] if len(parts) > 1 else "") if value}
    if len(attempts) != 1 or len(compounds) > 1:
        raise ValueError("The calculation has missing or conflicting original attempt selectors.")
    attempt = next(iter(attempts))
    compound = next(iter(compounds), "")
    if any(Path(value).name != value or value in {".", ".."} for value in (attempt, compound) if value):
        raise ValueError("The original attempt or compound selector is not a basename.")
    return attempt, compound


def _model_components(document):
    components = []
    if not isinstance(document, dict):
        raise ValueError("The prediction input is not an object.")
    for entry in document.get("sequences", []):
        protein = (entry.get("protein") or entry.get("proteinChain")) if isinstance(entry, dict) else None
        if not isinstance(protein, dict):
            continue
        chains = protein.get("id") or protein.get("chain") or "A"
        chains = chains if isinstance(chains, list) else [chains]
        if not isinstance(protein.get("sequence"), str):
            raise ValueError("The protein input does not contain a literal sequence.")
        components.extend({"type": "protein", "chain": str(chain), "sequence": protein["sequence"], "role": "target"}
            for chain in chains)
    if not components:
        raise ValueError("No explicit protein target is recorded in the model input.")
    return canonical_components(components)[0]


def recover_target(run_dir, context, *, namespace="", cache=None):
    cache = cache if cache is not None else {}
    run_dir = Path(run_dir).resolve()
    if not run_dir.name.startswith("mnligand-") or run_dir.parent.name not in {"docking", "refolding"}:
        raise ValueError("Only preserved MN-Ligand docking/refolding imports are supported.")
    attempt_id, compound = _selector(context)
    attempt = _inside(run_dir, run_dir / "artifacts/replicates" / attempt_id)
    input_path = _inside(run_dir, attempt / "input.json")
    inputs = _json(input_path, cache)
    reference = inputs.get("target_artifact") or inputs.get("prepared_target") or inputs.get("target")
    if not isinstance(reference, dict) or reference.get("kind") != "run_artifact":
        raise ValueError("The original input has no explicit prepared-target reference.")
    if reference.get("role") not in {"receptor", "axis_aligned_receptor", "prepared_target"}:
        raise ValueError("The original artifact is not declared to be a receptor.")
    snapshots = [{"path": str(input_path), "sha256": _read(input_path, cache)[1]}]
    target_path = None
    if run_dir.parent.name == "docking":
        candidates = [attempt / relative for relative in ("input/receptor.pdb", "artifacts/input/receptor.pdb",
            "input/target.pdb", "artifacts/input/target.pdb", "artifacts/prepared/receptor.pdb")]
        candidates = [_inside(run_dir, path) for path in candidates if path.is_file()]
        if not candidates:
            raise ValueError("The actual copied docking receptor is unavailable.")
        if len({_read(path, cache)[1] for path in candidates}) != 1:
            raise ValueError("The attempt contains conflicting receptor snapshots.")
        target_path = candidates[0]
        content, digest = _read(target_path, cache)
        key = (digest, "receptor-components")
        if key not in cache:
            cache[key] = parse_structure(content, target_path.name, include_ligands=False,
                require_model_selection=True)[0]["components"]
        components = cache[key]
        components = canonical_components([{**component, "role": "target"} for component in components])[0]
        basis = "actual copied docking receptor"
        snapshots.append({"path": str(target_path), "sha256": digest})
        recorded = context.get("target_structure_digest") or context.get("target_variant_sha256")
        if recorded and recorded != digest:
            raise ValueError("The copied receptor conflicts with the calculation's recorded target digest.")
        workflow_role = "docking"
    else:
        if not compound:
            raise ValueError("The prediction has no exact original compound selector.")
        candidates = [attempt / "inputs" / (compound + suffix) for suffix in (".yaml", ".yml", ".json")]
        candidates = [_inside(run_dir, path) for path in candidates if path.is_file()]
        if len(candidates) != 1:
            raise ValueError("The exact compound's prediction input is missing or ambiguous.")
        model_path = candidates[0]
        content, digest = _read(model_path, cache)
        key = (digest, "model-components")
        if key not in cache:
            if model_path.suffix == ".json":
                document = json.loads(content)
            else:
                import yaml
                try:
                    document = yaml.safe_load(content)
                except yaml.YAMLError as exc:
                    raise ValueError("The preserved model input is not valid YAML.") from exc
            cache[key] = _model_components(document)
        components = cache[key]
        basis = "literal target chains in the exact compound's model input"
        snapshots.append({"path": str(model_path), "sha256": digest})
        workflow_role = "binding_affinity" if namespace.casefold().replace("-", "_") in {
            "nesso_1", "nesso_affinity"} else "independent_prediction"
    metadata_path = _inside(run_dir, run_dir / "metadata.json")
    metadata = _json(metadata_path, cache)
    snapshots.append({"path": str(metadata_path), "sha256": _read(metadata_path, cache)[1]})
    migration = metadata.get("legacy_migration") or {}
    pdb = str(metadata.get("pdb_id") or migration.get("pdb_id") or "").upper()
    from mn_molecular_design_studio.legacy_migration import THR_TARGETS
    receptor = str(metadata.get("receptor") or THR_TARGETS.get(pdb, {}).get("receptor") or "")
    if not re.fullmatch(r"[A-Z0-9]{4}", pdb):
        raise ValueError("The imported study has no recorded PDB identity for its target label.")
    length = sum(len(component.get("sequence", "")) for component in components)
    name = " · ".join(value for value in (receptor, pdb, str(length) + " aa") if value)
    evidence = {"policy": POLICY, "basis": basis, "attempt_id": attempt_id, "compound_id": compound,
        "input_artifact_sha256": str(reference.get("sha256") or ""),
        "contents_sha256": hashlib.sha256(json.dumps(contents_identity(components)).encode()).hexdigest(),
        "snapshots": snapshots}
    return {"name": name, "components": components, "evidence": evidence, "workflow_role": workflow_role,
        "target_path": str(target_path) if target_path else "", "pdb_id": pdb, "receptor": receptor}


def register_target(database, user, project_id, binding, *, source_key, origin):
    from mn_molecular_design_studio.molecular_database_ingest import _find_or_create_input
    record = {"name": binding["name"], "components": binding["components"], "files": [],
        "metadata": {"target_identity_type": "molecular_identity", "structural_reference": False,
            "pdb_id": binding["pdb_id"], "target_registration": binding["evidence"]}}
    molecule, _, _ = _find_or_create_input(database, user, project_id, binding["name"], source_key, record,
        provenance={**origin, "input_role": "modeled_target", "target_registration": binding["evidence"]},
        source_key_suffix=":original-target")
    return molecule


def plan_target_links(database, user, runs_root, *, run_ids=None):
    if user.get("role") != "admin":
        raise PermissionError("Target reconciliation requires an administrator.")
    from mn_molecular_design_studio.library_curation import _digest, encoded
    with database.connection() as connection:
        groups = {row["id"]: dict(row) for row in connection.execute("SELECT * FROM molecule_groups")
                  if database._can_read(row, user) and not json.loads(row["origin"]).get("_replica_node")}
        observations = [dict(row) for row in connection.execute(
            "SELECT * FROM molecule_observations WHERE group_id IN (SELECT value FROM json_each(?)) "
            "AND run_id LIKE 'mnligand-%'", (encoded(list(groups)),))
            if run_ids is None or row["run_id"] in run_ids]
        molecules = {row["id"]: json.loads(row["components"]) for row in connection.execute("SELECT id,components FROM molecules")}
    cache, targets, changes, unresolved = {}, {}, [], []
    for observation in observations:
        context = json.loads(observation["context"])
        task = context.get("task_group")
        if task not in {"docking", "refolding"}:
            continue
        run_dir = Path(runs_root) / task / observation["run_id"]
        if not run_dir.is_dir():
            unresolved.append({"observation_id": observation["id"], "reason": "The original run is unavailable on this workstation."})
            continue
        try:
            binding = recover_target(run_dir, context, namespace=observation["namespace"], cache=cache)
            components, fingerprint, _ = canonical_components(binding["components"])
            target_id = "mol-" + fingerprint[:32]
            original = context.get("target_molecule_id")
            if original and (original not in molecules or contents_identity(molecules[original]) != contents_identity(components)):
                raise ValueError("The recovered input conflicts with the existing target molecule link.")
            if context.get("target_pdb_id") and context["target_pdb_id"] != binding["pdb_id"]:
                raise ValueError("The imported study's PDB label conflicts with the recorded calculation target.")
            updated = {**context, "target_molecule_id": original or target_id, "target_pdb_id": binding["pdb_id"],
                "target_receptor": context.get("target_receptor") or binding["receptor"],
                "target_registration": binding["evidence"]}
            if not updated.get("workflow_role"):
                updated["workflow_role"] = binding["workflow_role"]
            if updated == context:
                continue
            changes.append({"id": observation["id"], "before": observation["context"], "after": encoded(updated)})
            group = groups[observation["group_id"]]
            targets[(group["owner_id"], group["project_id"], original or target_id)] = {
                "id": original or target_id, "name": binding["name"], "components": components,
                "owner_id": group["owner_id"], "project_id": group["project_id"]}
        except (OSError, ValueError, KeyError, TypeError) as exc:
            unresolved.append({"observation_id": observation["id"], "reason": str(exc)})
    plan = {"schema": "studio.library-target-links.v1", "database": str(database.path.resolve()), "members": [],
        "observations": changes, "targets": list(targets.values()), "unresolved": unresolved}
    plan["digest"] = _digest(plan)
    return plan
