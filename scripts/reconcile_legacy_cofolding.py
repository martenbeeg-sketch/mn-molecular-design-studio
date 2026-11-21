"""Reconcile selected MN-Ligand cofolding runs with their exact native model evidence."""
from __future__ import annotations

import argparse
from collections import Counter
from contextlib import closing
import hashlib
import json
import os
from pathlib import Path
import re
import sqlite3

from mn_molecular_design_studio.cofolding_results import legacy_ranked_prediction, structure_for
from mn_molecular_design_studio.data_browser import ReadOnlyLibrary
from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
from mn_molecular_design_studio.evaluations import file_digest, identifier
from mn_molecular_design_studio.molecular_database import MolecularDatabase, _now


POLICY = "studio.legacy-cofolding-model.v1"


def encoded(value):
    return json.dumps(value, sort_keys=True, allow_nan=False)


def _private_file(path):
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    return os.fdopen(descriptor, "w", encoding="utf-8")


def _source_rows(connection, runs):
    return [dict(row) for row in connection.execute(
        """SELECT o.* FROM molecule_observations o JOIN molecule_groups g ON g.id=o.group_id
           WHERE o.run_id IN (SELECT value FROM json_each(?))
           AND json_extract(o.context,'$.task_group')='refolding'
           AND o.namespace IN ('AlphaFold 3','Boltz-2')
           AND json_extract(g.origin,'$._replica_node') IS NULL
           AND g.import_method IN ('tool_inputs','tool_output') ORDER BY o.run_id,o.id""", (encoded(runs),)
    )]


def _protected_structures(connection):
    keys = {row[0] for row in connection.execute(
        "SELECT source_structure_id FROM complex_reference_contexts "
        "UNION SELECT reference_structure_id FROM complex_reference_contexts "
        "UNION SELECT structure_ref FROM complex_review_annotations "
        "UNION SELECT json_extract(o.context,'$.representation_id') "
        "FROM complex_review_annotations a JOIN molecule_observations o ON o.id=a.evaluation_id"
    )}
    references = Counter()
    for raw, in connection.execute("SELECT context FROM molecule_observations"):
        context = json.loads(raw)
        if representation := context.get("representation_id"):
            references[representation] += 1
        keys.update(value for key, value in context.items()
                    if key in {"source_structure_id", "target_structure_id", "source_structure_owner_structure_id"}
                    and isinstance(value, str))
    return keys, references


def plan_corrections(database, root, run_ids):
    root = root.resolve()
    manifests = {}
    for run_id in run_ids:
        if not re.fullmatch(r"mnligand-[a-f0-9]{16}", run_id):
            raise ValueError("Choose an existing MN-Ligand aggregate run ID.")
        run = (root / "refolding" / run_id).resolve()
        if not run.is_relative_to(root) or not (run / "mn_ligand_aggregate.json").is_file():
            raise ValueError(f"The original aggregate is unavailable: {run_id}")
        manifests[run_id] = (run, {row["candidate_id"]: row for row in read_candidates(run)})
    with database.connection() as connection:
        connection.execute("BEGIN")
        observations = _source_rows(connection, run_ids)
        structures = {row["id"]: dict(row) for row in connection.execute("SELECT * FROM molecule_structures")}
        protected, references = _protected_structures(connection)
    changes, unresolved = [], []
    ranking = Counter()
    for observation in observations:
        run, candidates = manifests[observation["run_id"]]
        context = json.loads(observation["context"])
        candidate = candidates.get(context.get("candidate_id"))
        selected = legacy_ranked_prediction(run, candidate) if candidate else None
        original = structure_for(run, candidate) if candidate else None
        identity = [observation["run_id"], context.get("candidate_id")]
        if selected is None or original is None:
            unresolved.append([*identity, "Selected model evidence unavailable"])
            continue
        source = Path(selected["complex_pdb"])
        old_candidate = context.get("candidate") or {}
        stored_source = Path(old_candidate.get("complex_pdb") or "")
        if not stored_source.is_absolute():
            stored_source = run / stored_source
        if stored_source.resolve() not in {original.resolve(), source.resolve()}:
            unresolved.append([*identity, "Saved candidate differs from the archived selection"])
            continue
        receipt = {"policy": POLICY, **selected["raw_metadata"]["legacy_model_selection"]}
        ranking[receipt["rank_status"]] += 1
        representation_id = context.get("representation_id")
        previous = structures.get(representation_id)
        if representation_id and (previous is None or previous["molecule_id"] != observation["molecule_id"]
                                  or previous["group_id"] != observation["group_id"]):
            unresolved.append([*identity, "Saved structure reference is inconsistent"])
            continue
        replacement = None
        digest = file_digest(source)
        if not digest:
            unresolved.append([*identity, "Selected structure file unavailable"])
            continue
        if previous and stored_source.resolve() != source.resolve():
            if representation_id in protected or references[representation_id] != 1:
                unresolved.append([*identity, "The previous structure has another saved reference"])
                continue
            if context.get("structure_digest") != file_digest(original):
                unresolved.append([*identity, "The previous representation digest differs from the original model"])
                continue
            structure_context = json.loads(previous["context"])
            structure_context.update(content_sha256=digest, portable_sha256=digest, validated_content_sha256=digest,
                                     portable=True, legacy_model_selection=receipt)
            replacement_id = identifier("structure-", [representation_id, POLICY, digest])
            replacement = {**previous, "id": replacement_id, "format": source.suffix.lower().lstrip("."),
                           "label": f"{observation['namespace']} best complex",
                           "context": encoded(structure_context)}
            context["representation_id"] = replacement_id
        stored_metrics = selected["metrics"]
        candidate_context = {**old_candidate, "complex_pdb": str(source), "model_index": receipt["model_index"]}
        candidate_context["raw_metadata"] = {**(old_candidate.get("raw_metadata") or {}), "legacy_model_selection": receipt}
        context.update(candidate=candidate_context, legacy_model_selection=receipt,
                       model_id=str(receipt["model_index"]), structure_digest=digest)
        for field in ("source_structure_digest", "calculation_structure_digest"):
            if field in context:
                context[field] = digest
        if json.loads(observation["metrics"]) == stored_metrics and json.loads(observation["context"]) == context:
            continue
        changes.append({"before": observation, "after_context": encoded(context), "after_metrics": encoded(stored_metrics),
                        "replacement": replacement, "source": str(source), "sha256": digest,
                        "previous_structure": previous if replacement else None})
    return {"policy": POLICY, "database": str(database.path.resolve()), "root": str(root),
            "runs": run_ids, "examined": len(observations), "ranking": ranking,
            "changes": changes, "unresolved": unresolved}


def apply_corrections(database, plan, backup):
    if plan["database"] != str(database.path.resolve()):
        raise ValueError("The plan belongs to a different Database.")
    reviewed = plan_corrections(ReadOnlyLibrary(database.root), Path(plan["root"]), plan["runs"])
    if reviewed != plan:
        raise ValueError("Saved records or ranked model evidence changed after planning.")
    with _private_file(backup):
        pass
    with database.connection() as connection, closing(sqlite3.connect(backup)) as destination:
        connection.backup(destination)
        if destination.execute("PRAGMA quick_check").fetchone()[0] != "ok":
            raise ValueError("The recovery backup failed its integrity check.")
    for change in plan["changes"]:
        if change["replacement"]:
            source = Path(change["source"])
            content = source.read_bytes()
            if hashlib.sha256(content).hexdigest() != change["sha256"]:
                raise ValueError("The selected model changed after review.")
            change["replacement"]["path"] = str(database._store_file(source.name, content))
    with database.connection() as connection:
        connection.execute("BEGIN IMMEDIATE")
        protected, references = _protected_structures(connection)
        for change in plan["changes"]:
            old = change["before"]
            current = connection.execute(
                "SELECT o.* FROM molecule_observations o JOIN molecule_groups g ON g.id=o.group_id "
                "WHERE o.id=? AND json_extract(g.origin,'$._replica_node') IS NULL", (old["id"],)
            ).fetchone()
            if current is None or dict(current) != old:
                raise ValueError("A reviewed observation changed before application.")
            if change["replacement"]:
                previous = change["previous_structure"]
                saved = connection.execute("SELECT * FROM molecule_structures WHERE id=?", (previous["id"],)).fetchone()
                if (saved is None or dict(saved) != previous or previous["id"] in protected
                        or references[previous["id"]] != 1):
                    raise ValueError("A previous structure or its links changed before application.")
                fields = ("id", "molecule_id", "group_id", "path", "label", "format", "context", "created_at")
                connection.execute("INSERT INTO molecule_structures VALUES (?,?,?,?,?,?,?,?)",
                                   tuple(change["replacement"][key] for key in fields))
            connection.execute("UPDATE molecule_observations SET metrics=?,context=? WHERE id=?",
                               (change["after_metrics"], change["after_context"], old["id"]))
            if change["replacement"]:
                connection.execute("DELETE FROM molecule_structures WHERE id=?", (previous["id"],))
        stamp = _now()
        connection.executemany("UPDATE molecule_groups SET updated_at=? WHERE id=?",
                               [(stamp, group) for group in sorted({row["before"]["group_id"] for row in plan["changes"]})])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("plan", "apply"))
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--source-runs-root", type=Path)
    parser.add_argument("--run-id", action="append")
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--backup", type=Path)
    parser.add_argument("--admin-id", required=True)
    auth = parser.parse_args()
    with closing(sqlite3.connect((auth.root / "user_access.sqlite3").resolve().as_uri() + "?mode=ro", uri=True)) as connection:
        actor = connection.execute("SELECT role,disabled FROM users WHERE id=?", (auth.admin_id,)).fetchone()
    if actor != ("admin", 0):
        raise PermissionError("Choose an enabled Studio administrator.")
    if auth.mode == "plan":
        if not auth.source_runs_root or not auth.run_id:
            raise ValueError("Choose the source runs root and explicit run IDs.")
        plan = plan_corrections(ReadOnlyLibrary(auth.root), auth.source_runs_root, sorted(set(auth.run_id)))
        with _private_file(auth.plan) as handle:
            handle.write(encoded(plan))
    else:
        if not auth.backup:
            raise ValueError("Choose a separate recovery backup path.")
        plan = json.loads(auth.plan.read_text(encoding="utf-8"))
        apply_corrections(MolecularDatabase(auth.root), plan, auth.backup)
    print(encoded({"examined": plan["examined"], "observations": len(plan["changes"]),
                   "structures": sum(row["replacement"] is not None for row in plan["changes"]),
                   "ranking": plan["ranking"], "unresolved": plan["unresolved"]}))


if __name__ == "__main__":
    main()
