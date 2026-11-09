"""Run-scoped result access and exact-subset promotion into the molecular library."""
from __future__ import annotations

import json
import hashlib
import socket
from pathlib import Path

from mn_molecular_design_studio.molecular_database import MolecularDatabase, _id, _json, _now


def _authorized_run(run_dir, user):
    from mn_molecular_design_studio.auth import can_access_job
    from mn_molecular_design_studio.job_catalog import get_job_access_record
    from mn_molecular_design_studio.runtime import runs_root
    from mn_molecular_design_studio.workspace import is_studio_job

    root = Path(runs_root()).expanduser().resolve()
    selected = Path(run_dir).expanduser().resolve()
    try:
        relative = selected.relative_to(root)
    except ValueError as exc:
        raise FileNotFoundError("Studio run not found.") from exc
    if len(relative.parts) != 2 or any(part in {"", ".", ".."} for part in relative.parts):
        raise FileNotFoundError("Studio run not found.")
    task_group, run_id = relative.parts
    if not selected.is_dir() or not is_studio_job(selected):
        raise FileNotFoundError("Studio run not found.")
    job = get_job_access_record(task_group, run_id)
    if not can_access_job(user, job):
        raise PermissionError("This Studio run is not accessible.")
    if job:
        catalog_path = Path(job.get("run_relative_path") or "")
        if catalog_path.parts != (task_group, run_id):
            raise FileNotFoundError("Studio run not found.")
    marker_path = selected / "molecular_design_studio.json"
    try:
        marker = json.loads(marker_path.read_text())
    except (OSError, ValueError):
        marker = {}
    owner_id = str(marker.get("studio_owner_id") or "")
    project_id = str(marker.get("project_id") or "")
    if job:
        if owner_id and job.get("owner_id") and owner_id != str(job["owner_id"]):
            raise PermissionError("This Studio run is not accessible.")
        if project_id and job.get("project_id") and project_id != str(job["project_id"]):
            raise PermissionError("This Studio run is not accessible.")
        owner_id = str(job.get("owner_id") or owner_id)
        project_id = str(job.get("project_id") or project_id or "")
    if not owner_id:
        metadata_path = selected / "metadata.json"
        try:
            metadata = json.loads(metadata_path.read_text())
        except (OSError, ValueError):
            metadata = {}
        owner_id = str(metadata.get("owner_id") or (metadata.get("molecular_design_studio") or {}).get("studio_owner_id") or "")
        project_id = project_id or str(metadata.get("project_id") or (metadata.get("molecular_design_studio") or {}).get("project_id") or "")
    if not owner_id:
        raise FileNotFoundError("Studio run not found.")
    return selected, task_group, run_id, owner_id, project_id


def _run_sources(database, task_group, run_id):
    prefix = f"{socket.gethostname()}/{task_group}/{run_id}:"
    with database.connection() as connection:
        rows = connection.execute(
            "SELECT source_key,molecule_id,group_id,provenance FROM molecule_sources "
            "WHERE substr(source_key,1,?)=? AND json_extract(provenance,'$.saved_selection') IS NOT 1 ORDER BY source_key",
            (len(prefix), prefix),
        ).fetchall()
        groups = connection.execute(
            "SELECT id,name,owner_id,project_id,import_method,origin FROM molecule_groups "
            "WHERE id IN (SELECT DISTINCT group_id FROM molecule_sources WHERE substr(source_key,1,?)=? AND json_extract(provenance,'$.saved_selection') IS NOT 1)",
            (len(prefix), prefix),
        ).fetchall()
    return [dict(row) for row in rows], [dict(row) for row in groups]


def _result_candidate_token(task_group, run_id, candidate_id, model_id=""):
    return _id("result-", _json([socket.gethostname(), str(task_group), str(run_id), str(candidate_id), str(model_id or "")]))


def _candidate_for_explicit_selection(candidate, run_dir):
    """Resolve one authorized Results candidate into its selected library identity."""
    from mn_molecular_design_studio.molecular_database_ingest import _candidate_acceptance_state, _record

    record = _record(candidate, run_dir)
    if record is None:
        return None
    identity = dict(record)
    identity["files"] = list(record.get("files") or [])
    # Distinguish the prediction's model identity from the coordinate file's
    # internal MODEL number, which _record stores as ``model_id``.
    candidate_id = str(candidate.get("candidate_id") or "")
    candidate_model_id = _candidate_model_id(candidate)
    raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
    target_identity, target_sequences = _candidate_target_provenance(candidate, record.get("components"))
    generator_engine, prediction_engine = _candidate_engine_provenance(candidate)
    for item in identity["files"]:
        item.update({
            "candidate_id": candidate_id,
            "candidate_model_id": candidate_model_id,
            "source_model_id": candidate_model_id,
            "source_tool": str(candidate.get("source_tool") or raw.get("generator") or ""),
            "generator_engine": generator_engine,
            "prediction_engine": prediction_engine,
            "target_identity": target_identity,
            "target_sequences_by_chain": target_sequences,
            "target_chains": list(candidate.get("target_chains") or []),
            "binder_chains": list(candidate.get("binder_chains") or []),
            "parent_ids": list(candidate.get("parent_ids") or candidate.get("parents") or []),
        })
    sequence = candidate.get("binder_sequence") or candidate.get("sequence")
    binder_chains = {str(value) for value in candidate.get("binder_chains") or []}
    if sequence:
        binders = [component for component in record.get("components") or []
                   if component.get("role") == "binder" or str(component.get("chain") or "") in binder_chains]
        if binders:
            identity["components"] = binders
    state = _candidate_acceptance_state(candidate)
    raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
    identity["metadata"] = {
        **(identity.get("metadata") or {}),
        "result_selection": True,
        "original_acceptance_status": state,
        "acceptance_evidence": raw.get("validation_acceptance") or raw.get("distributed_campaign") or {},
        "target_identity": target_identity or "Not recorded",
        "target_sequences_by_chain": target_sequences,
    }
    return identity


def _candidate_model_id(candidate: dict) -> str:
    raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
    selection = raw.get("complex_analysis_selection") if isinstance(raw.get("complex_analysis_selection"), dict) else {}
    for source in (candidate, raw, selection):
        for key in ("model_id", "model_index"):
            value = source.get(key)
            if value is not None and str(value).strip():
                return str(value).strip()
    return ""


def _candidate_engine_provenance(candidate: dict) -> tuple[str, str]:
    """Return explicit generator and prediction identities without guessing."""
    raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
    source_candidate = raw.get("source_candidate") if isinstance(raw.get("source_candidate"), dict) else {}
    candidate_set = raw.get("candidate_set") if isinstance(raw.get("candidate_set"), dict) else {}
    generator = next((str(value).strip() for value in (
        raw.get("generator_engine"), raw.get("design_engine"), source_candidate.get("generator_engine"),
        source_candidate.get("design_engine"), source_candidate.get("source_tool"),
        candidate_set.get("generator"), candidate_set.get("source_tool"),
    ) if value is not None and str(value).strip()), "")
    if not generator and not source_candidate:
        try:
            from mn_molecular_design_studio.native_filter_metrics import native_filter_evidence
            is_generator_output = str(candidate.get("stage") or "").startswith("generation.")
            has_native_gate = bool(native_filter_evidence(candidate).get("source"))
        except (ImportError, TypeError, ValueError):
            is_generator_output = has_native_gate = False
        if is_generator_output or has_native_gate:
            generator = str(candidate.get("source_tool") or "").strip()
    prediction = next((str(value).strip() for value in (
        raw.get("prediction_engine"), raw.get("folding_engine"),
        source_candidate.get("prediction_engine"), source_candidate.get("folding_engine"),
    ) if value is not None and str(value).strip()), "")
    return generator, prediction


def _candidate_target_provenance(candidate: dict, components=None) -> tuple[str, dict[str, str]]:
    """Preserve the exact recorded target-chain sequences without naming inference."""
    raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
    target_chains = [str(value) for value in candidate.get("target_chains") or [] if str(value)]
    sequences = raw.get("sequences_by_chain") if isinstance(raw.get("sequences_by_chain"), dict) else {}
    if not sequences:
        raw_components = raw.get("components") if isinstance(raw.get("components"), list) else []
        sequences = {str(row.get("chain") or ""): str(row.get("sequence") or "")
                     for row in raw_components if isinstance(row, dict) and row.get("chain") and row.get("sequence")}
    if not any(sequences.get(chain) for chain in target_chains) and isinstance(components, list):
        sequences = {str(row.get("chain") or ""): str(row.get("sequence") or "")
                     for row in components if isinstance(row, dict) and row.get("chain") and row.get("sequence")}
    target_sequences = {chain: str(sequences.get(chain) or "") for chain in target_chains
                        if str(sequences.get(chain) or "").strip()}
    target_identity = next((str(value).strip() for value in (
        raw.get("target_identity"), raw.get("target_id"), candidate.get("target_identity"), candidate.get("target_id")
    ) if value is not None and str(value).strip()), "")
    if not target_identity and target_sequences:
        descriptions = [f"chain {chain} · {len(sequence)} aa" for chain, sequence in sorted(target_sequences.items())]
        target_identity = "Target sequence (" + ", ".join(descriptions) + ")"
    return target_identity, target_sequences


def run_acceptance_inventory(run_dir, user) -> dict:
    """Summarize only explicit per-candidate acceptance decisions for one authorized run."""
    selected, task_group, run_id, owner_id, project_id = _authorized_run(run_dir, user)
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
    from mn_molecular_design_studio.complex_analysis_collection import acceptance_decisions

    try:
        metadata = json.loads((selected / "metadata.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        metadata = {}
    candidates = [row for row in read_candidates(selected) if isinstance(row, dict) and row.get("candidate_id")]
    rows = []
    seen = set()
    for candidate in candidates:
        candidate_id = str(candidate.get("candidate_id") or "")
        model_id = _candidate_model_id(candidate)
        key = (candidate_id, model_id)
        if key in seen:
            continue
        seen.add(key)
        raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
        generator_engine, prediction_engine = _candidate_engine_provenance(candidate)
        target_identity, target_sequences = _candidate_target_provenance(candidate)
        rows.append({
            "candidate_id": candidate_id,
            "model_id": model_id,
            "source_tool": str(candidate.get("source_tool") or raw.get("generator") or metadata.get("tool") or ""),
            "generator_engine": generator_engine,
            "prediction_engine": prediction_engine,
            "target_identity": target_identity,
            "target_sequences_by_chain": target_sequences,
            "decisions": acceptance_decisions(candidate),
        })
    # Summaries are accumulated from the shared explicit-decision parser above;
    # status scores are never promoted into acceptance evidence here.
    stage_labels = {}
    order = ("final", "sequence", "shared_filter", "native")
    for row in rows:
        for decision in row["decisions"]:
            stage_id = str(decision.get("stage_id") or "")
            if stage_id:
                stage_labels.setdefault(stage_id, str(decision.get("stage") or stage_id))
    stage_ids = [stage_id for stage_id in order if stage_id in stage_labels]
    stage_summaries = []
    for stage_id in stage_ids:
        passing, failed, pending, unrecorded = [], [], [], []
        for row in rows:
            decisions = [item for item in row["decisions"] if item.get("stage_id") == stage_id]
            if not decisions:
                unrecorded.append(row)
                continue
            # Multiple contradictory decisions on one candidate/model are not
            # collapsed into a pass; they remain unresolved for collection.
            statuses = {str(item.get("cohort_status") or "pending") for item in decisions}
            if statuses == {"passed"}:
                passing.append({**row, "decision": decisions[0]})
            elif statuses == {"rejected"}:
                failed.extend({**row, "decision": item} for item in decisions)
            else:
                pending.extend({**row, "decision": item} for item in decisions)
        stage_summaries.append({
            "stage_id": stage_id, "label": stage_labels[stage_id],
            "passed": passing, "failed": failed, "pending": pending, "unrecorded": unrecorded,
            "passed_design_count": len({item["candidate_id"] for item in passing}),
            "passed_model_count": len(passing),
            "failed_model_count": len(failed), "pending_model_count": len(pending),
            "unrecorded_model_count": len(unrecorded),
        })
    return {
        "task_group": task_group, "run_id": run_id, "owner_id": owner_id, "project_id": project_id,
        "candidate_count": len({row["candidate_id"] for row in rows}), "model_count": len(rows),
        "stages": stage_summaries, "metadata": metadata,
        "source_tool": str(metadata.get("tool") or metadata.get("display_tool") or ""),
        "candidates": rows,
    }


def _unindexed_result_rows(run_dir, task_group, run_id, owner_id, project_id, indexed_sources, items):
    """Expose valid candidates in Results without promoting or ingesting them."""
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
    from mn_molecular_design_studio.molecular_database_ingest import _candidate_acceptance_state, _record
    from mn_molecular_design_studio.molecule_imports import canonical_components

    indexed_candidates = set()
    for source in indexed_sources:
        try:
            provenance = json.loads(source["provenance"])
        except (TypeError, ValueError):
            continue
        if provenance.get("candidate_id"):
            indexed_candidates.add((str(provenance["candidate_id"]), str(provenance.get("model_id") or "")))
    for item in items:
        for structure in item.get("structures") or []:
            context = structure.get("context") if isinstance(structure.get("context"), dict) else {}
            candidate_id = str(context.get("candidate_id") or context.get("native_candidate_id") or "")
            model_id = str(context.get("model_id") or context.get("model_index") or context.get("validated_model_id") or "")
            if candidate_id:
                indexed_candidates.add((candidate_id, model_id))
    for candidate in read_candidates(run_dir):
        candidate_id = str(candidate.get("candidate_id") or "")
        model_id = _candidate_model_id(candidate)
        if not candidate_id or (candidate_id, model_id) in indexed_candidates:
            continue
        record = _record(candidate, run_dir)
        if record is None:
            continue
        identity = _candidate_for_explicit_selection(candidate, run_dir)
        if identity is None:
            continue
        components, fingerprint, molecule_type = canonical_components(identity["components"])
        if identity.get("identity_context"):
            fingerprint = hashlib.sha256((fingerprint + _json(identity["identity_context"])).encode()).hexdigest()
        molecule_id = "mol-" + fingerprint[:32]
        files = identity.get("files") or []
        representations = []
        for index, item in enumerate(files):
            path = Path(item["path"])
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            representation_id = _id("result-structure-", _json([candidate_id, model_id, digest, index]))
            representations.append({
                "id": representation_id,
                "label": str(item.get("label") or path.name),
                "format": path.suffix.lower().lstrip("."),
                "digest": digest,
                "model_id": model_id,
                "available_in_results": True,
            })
        raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
        items.append({
            "id": _result_candidate_token(task_group, run_id, candidate_id, model_id),
            "molecule_id": molecule_id,
            "name": str(identity.get("name") or candidate_id),
            "type": molecule_type,
            "components": components,
            "metadata": identity.get("metadata") or {},
            "groups": [],
            "representations": representations,
            "structures": representations,
            "observations": [],
            "scores": dict(candidate.get("metrics") or {}),
            "run_group_ids": [],
            "run_groups": [],
            "source_group_id": None,
            "result_only": True,
            "candidate_id": candidate_id,
            "model_id": model_id,
            "acceptance_status": _candidate_acceptance_state(candidate),
            "source_run": {"task_group": task_group, "run_id": run_id,
                           "owner_id": owner_id, "project_id": project_id},
            "has_structure": bool(representations),
            "native_metadata": {key: raw[key] for key in (
                "validation_acceptance", "distributed_campaign", "engine_version", "structure_kind") if key in raw},
        })


def _list_authorized_run(run_dir, user, database):
    selected, task_group, run_id, owner_id, project_id = _authorized_run(run_dir, user)
    # This is a read path. Candidate previews remain in Results until a user
    # explicitly promotes them or a separate completion/index operation runs.
    sources, source_groups = _run_sources(database, task_group, run_id)
    molecule_ids = list(dict.fromkeys(row["molecule_id"] for row in sources))
    result = (database.list_molecules(user, molecule_ids=molecule_ids, collection="all",
                                      limit=min(100_000, max(1, len(molecule_ids))))
              if molecule_ids else {"items": []})
    accessible_groups = {}
    for group in source_groups:
        if database._can_read(group, user):
            try:
                origin = json.loads(group["origin"])
            except (TypeError, ValueError):
                origin = {}
            accessible_groups[group["id"]] = {
                "id": group["id"], "name": group["name"], "project_id": group["project_id"],
                "collection_kind": database._collection_kind({"origin": origin}),
                "origin": origin,
            }
    by_molecule = {}
    for source in sources:
        if source["group_id"] in accessible_groups:
            by_molecule.setdefault(source["molecule_id"], set()).add(source["group_id"])
    items = []
    for record in result["items"]:
        from mn_molecular_design_studio.evaluations import project_observation, evaluation_descriptor
        record["structures"] = [s for s in record.get("representations", [])
                                 if (s.get("context") or {}).get("run_id") == run_id]
        run_observations = [observation for observation in record["observations"] if observation["run_id"] == run_id]
        if run_observations:
            run_observations.sort(key=lambda obs: (not evaluation_descriptor(obs)["preferred"], obs["id"]))
            chosen = run_observations[0]
            record["scores"], record["score_sources"], _ = project_observation(chosen)
            record["selected_evaluation"] = evaluation_descriptor(chosen)
            record["observations"] = run_observations
        run_group_ids = sorted(by_molecule.get(record["id"], set()))
        run_group_ids.sort(key=lambda group_id: (
            {"generated": 0, "curated": 1, "inputs": 2, "evaluation": 3, "intermediate": 4}.get(
                accessible_groups[group_id]["collection_kind"], 5), group_id))
        record["run_group_ids"] = run_group_ids
        record["run_groups"] = [accessible_groups[group_id] for group_id in run_group_ids]
        record["source_group_id"] = run_group_ids[0] if run_group_ids else None
        items.append(record)
    _unindexed_result_rows(selected, task_group, run_id, owner_id, project_id, sources, items)
    return selected, task_group, run_id, owner_id, project_id, items, accessible_groups


def list_run_molecules(run_dir, user, *, database=None) -> list[dict]:
    """Return normal molecular records associated with this accessible Studio run."""
    database = database or MolecularDatabase()
    return _list_authorized_run(run_dir, user, database)[5]


def list_saved_run_collections(run_dir, user, *, database=None) -> list[dict]:
    """Find accessible explicit selections linked to an authorized source run.

    Saved selections remain separate from automatic run outputs: a Results
    cleanup must preserve them, while its navigation may still show the link.
    """
    database = database or MolecularDatabase()
    _, task_group, run_id, _, _ = _authorized_run(run_dir, user)
    from mn_molecular_design_studio.cluster import config

    with database.connection() as connection:
        rows = connection.execute(
            "SELECT g.*,s.molecule_id,s.provenance FROM molecule_groups g "
            "JOIN molecule_sources s ON s.group_id=g.id "
            "WHERE json_extract(s.provenance,'$.task_group')=? "
            "AND json_extract(s.provenance,'$.run_id')=? "
            "AND json_extract(s.provenance,'$.saved_selection')=1 "
            "ORDER BY g.name,g.id,s.molecule_id", (task_group, run_id),
        ).fetchall()
    groups = {}
    for row in rows:
        if not database._can_read(row, user):
            continue
        provenance = json.loads(row["provenance"])
        source_node = provenance.get("source_node")
        if source_node and source_node != config.identity():
            continue
        if not source_node and provenance.get("workstation") not in (None, "", socket.gethostname()):
            continue
        item = groups.setdefault(row["id"], {
            "id": row["id"], "name": row["name"], "project_id": row["project_id"],
            "origin": json.loads(row["origin"]), "molecule_ids": [],
        })
        if row["molecule_id"] not in item["molecule_ids"]:
            item["molecule_ids"].append(row["molecule_id"])
    return list(groups.values())


def remove_run_records(task_group, run_id, user, *, owner_id=None, project_id=None, database=None) -> dict:
    """Remove only automatically indexed Database records originating in a run.

    Run files and explicitly promoted/saved molecule selections are preserved.
    The job identity is authorized before any Database rows are changed.
    """
    from mn_molecular_design_studio.auth import can_access_job

    task_group, run_id = str(task_group or ""), str(run_id or "")
    if not task_group or not run_id or any(value in {".", ".."} or "/" in value or "\\" in value for value in (task_group, run_id)):
        raise ValueError("Invalid job reference.")
    if not can_access_job(user, {"owner_id": owner_id, "project_id": project_id}):
        raise PermissionError("This job is not accessible.")
    database = database or MolecularDatabase()
    with database.connection() as connection:
        source_rows = connection.execute(
            "SELECT source_key,molecule_id,group_id FROM molecule_sources "
            "WHERE json_extract(provenance,'$.task_group')=? AND json_extract(provenance,'$.run_id')=? "
            "AND json_extract(provenance,'$.saved_selection') IS NOT 1",
            (task_group, run_id),
        ).fetchall()
        group_rows = {
            str(row["id"]): dict(row)
            for row in connection.execute(
                "SELECT * FROM molecule_groups WHERE id IN (SELECT DISTINCT group_id FROM molecule_sources "
                "WHERE json_extract(provenance,'$.task_group')=? AND json_extract(provenance,'$.run_id')=? "
                "AND json_extract(provenance,'$.saved_selection') IS NOT 1)",
                (task_group, run_id),
            ).fetchall()
        }
    # Never remove a collaborator's Database rows merely because the job is
    # visible through a shared project.
    for group_id in group_rows:
        database.group(group_id, user, write=True)

    pairs = {(str(row["group_id"]), str(row["molecule_id"])) for row in source_rows}
    removed_members = 0
    with database.connection() as connection:
        for row in source_rows:
            connection.execute("DELETE FROM molecule_sources WHERE source_key=?", (row["source_key"],))
        for group_id, molecule_id in pairs:
            for table in ("molecule_structures", "molecule_observations"):
                connection.execute(
                    f"DELETE FROM {table} WHERE group_id=? AND molecule_id=? "
                    "AND json_extract(context,'$.task_group')=? AND json_extract(context,'$.run_id')=?",
                    (group_id, molecule_id, task_group, run_id),
                )
            origin_row = group_rows.get(group_id)
            if origin_row:
                try:
                    origin = json.loads(origin_row["origin"])
                except (TypeError, ValueError):
                    origin = {}
                if origin.get("type") == "job" and origin.get("task_group") == task_group and origin.get("run_id") == run_id:
                    relations_remain = any(
                        connection.execute(f"SELECT 1 FROM {table} WHERE group_id=? AND molecule_id=? LIMIT 1", (group_id, molecule_id)).fetchone()
                        for table in ("molecule_sources", "molecule_structures", "molecule_observations")
                    )
                    if not relations_remain:
                        cursor = connection.execute("DELETE FROM molecule_members WHERE group_id=? AND molecule_id=?", (group_id, molecule_id))
                        removed_members += max(0, cursor.rowcount)
        removable_groups = []
        for group_id, group in group_rows.items():
            try:
                origin = json.loads(group["origin"])
            except (TypeError, ValueError):
                origin = {}
            if (origin.get("type") == "job" and origin.get("task_group") == task_group
                    and origin.get("run_id") == run_id):
                count = connection.execute("SELECT count(*) FROM molecule_members WHERE group_id=?", (group_id,)).fetchone()[0]
                if not count:
                    connection.execute("DELETE FROM molecule_groups WHERE id=?", (group_id,))
                    removable_groups.append(group_id)
    return {
        "removed_sources": len(source_rows),
        "removed_members": removed_members,
        "removed_groups": len(removable_groups),
        "preserved_saved_selections": True,
    }


def _candidate_model_from_context(context: dict) -> tuple[str, str]:
    nested = context.get("candidate") if isinstance(context.get("candidate"), dict) else {}
    candidate_id = str(context.get("candidate_id") or context.get("native_candidate_id")
                       or nested.get("candidate_id") or "")
    # Prefer the prediction model identity over the coordinate file's internal
    # MODEL selector when both are present.
    model_id = str(context.get("candidate_model_id") or context.get("source_model_id")
                   or context.get("model_id") or context.get("model_index") or context.get("validated_model_id")
                   or nested.get("model_id") or nested.get("model_index") or "")
    return candidate_id, model_id


def _matches_candidate_models(candidate_id: str, model_id: str, candidate_models) -> bool:
    if candidate_models is None:
        return True
    matches = [row for row in candidate_models if str(row.get("candidate_id") or "") == str(candidate_id or "")]
    if not matches:
        return False
    selected_models = {str(row.get("model_id") or "") for row in matches}
    if "" in selected_models:
        return True
    return bool(model_id) and str(model_id) in selected_models


def _models_on_run_record(record: dict) -> set[tuple[str, str]]:
    if record.get("result_only"):
        return {(str(record.get("candidate_id") or ""), str(record.get("model_id") or ""))}
    found = set()
    for structure in record.get("structures") or []:
        context = structure.get("context") if isinstance(structure.get("context"), dict) else {}
        pair = _candidate_model_from_context(context)
        if pair[0]:
            found.add(pair)
    return found


def _copy_run_associations(database, task_group, run_id, source_keys, source_group_ids,
                           destination_group_id, molecule_ids, *, candidate_models=None):
    if not molecule_ids:
        return
    id_json = _json(molecule_ids)
    group_json = _json(source_group_ids)
    run_context = {"run_id": run_id, "task_group": task_group}
    with database.connection() as connection:
        structures = connection.execute(
            "SELECT * FROM molecule_structures WHERE molecule_id IN (SELECT value FROM json_each(?)) "
            "AND group_id IN (SELECT value FROM json_each(?)) ORDER BY created_at,id", (id_json, group_json)
        ).fetchall()
        observations = connection.execute(
            "SELECT * FROM molecule_observations WHERE molecule_id IN (SELECT value FROM json_each(?)) "
            "AND group_id IN (SELECT value FROM json_each(?)) ORDER BY created_at,id", (id_json, group_json)
        ).fetchall()
        sources = connection.execute(
            "SELECT * FROM molecule_sources WHERE source_key IN (SELECT value FROM json_each(?)) "
            "AND molecule_id IN (SELECT value FROM json_each(?))", (_json(source_keys), id_json)
        ).fetchall()
        for row in structures:
            context = json.loads(row["context"])
            if context.get("run_id") != run_id or context.get("task_group") != task_group:
                continue
            candidate_id, model_id = _candidate_model_from_context(context)
            if not _matches_candidate_models(candidate_id, model_id, candidate_models):
                continue
            identifier = _id("structure-", _json([row["id"], destination_group_id, run_id, task_group]))
            connection.execute("INSERT OR IGNORE INTO molecule_structures VALUES (?,?,?,?,?,?,?,?)", (
                identifier, row["molecule_id"], destination_group_id, row["path"], row["label"], row["format"], row["context"], row["created_at"]))
        for row in observations:
            context = json.loads(row["context"])
            if context.get("run_id") != run_id or context.get("task_group") != task_group:
                continue
            candidate_id, model_id = _candidate_model_from_context(context)
            if not _matches_candidate_models(candidate_id, model_id, candidate_models):
                continue
            identifier = _id("score-", f"saved:{destination_group_id}:{row['id']}")
            connection.execute("INSERT OR IGNORE INTO molecule_observations VALUES (?,?,?,?,?,?,?,?)", (
                identifier, row["molecule_id"], destination_group_id, row["namespace"], row["run_id"],
                row["metrics"], row["context"], row["created_at"]))
        for row in sources:
            provenance = json.loads(row["provenance"])
            if candidate_models is not None and not _matches_candidate_models(
                str(provenance.get("candidate_id") or ""), str(provenance.get("model_id") or ""), candidate_models
            ):
                continue
            identifier = row["source_key"] + ":saved:" + destination_group_id
            provenance["saved_selection"] = True
            provenance["library_source_key"] = row["source_key"]
            connection.execute("INSERT OR IGNORE INTO molecule_sources VALUES (?,?,?,?)", (
                identifier, row["molecule_id"], destination_group_id, _json(provenance)))


def save_run_selection(run_dir, user, molecule_ids, *, database=None, group_id=None,
                       name="Saved results", project_id=None, candidate_models=None,
                       collection_metadata=None) -> dict:
    """Promote exactly the selected run records into a curated user collection."""
    database = database or MolecularDatabase()
    selected, task_group, run_id, owner_id, run_project_id, items, source_groups = _list_authorized_run(run_dir, user, database)
    requested = list(dict.fromkeys(str(value) for value in molecule_ids or []))
    if not requested:
        raise ValueError("Select at least one run molecule to save.")
    records_by_id = {row["id"]: row for row in items}
    if set(requested) - set(records_by_id):
        raise PermissionError("Some selected molecules are not associated with this run or are inaccessible.")
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
    candidates_by_key = {}
    candidates_by_id = {}
    for candidate in read_candidates(selected):
        candidate_id = str(candidate.get("candidate_id") or "")
        if not candidate_id:
            continue
        key = (candidate_id, _candidate_model_id(candidate))
        if key in candidates_by_key:
            if candidate_models is not None:
                raise ValueError("The Results run contains duplicate candidate/model identities.")
            continue
        candidates_by_key[key] = candidate
        candidates_by_id.setdefault(candidate_id, []).append(candidate)
    selected_candidate_models = None if candidate_models is None else [
        {"candidate_id": str(row.get("candidate_id") or ""), "model_id": str(row.get("model_id") or ""),
         "decision": dict(row.get("decision") or {})}
        for row in candidate_models if isinstance(row, dict) and row.get("candidate_id")
    ]
    if selected_candidate_models is not None:
        if not selected_candidate_models:
            raise ValueError("No passing candidate/model rows were selected for this collection.")
        if isinstance(collection_metadata, dict):
            required_stage = str(collection_metadata.get("acceptance_stage_id") or "")
            for row in selected_candidate_models:
                decision = row.get("decision") or {}
                if (str(decision.get("cohort_status") or "") != "passed"
                        or (required_stage and str(decision.get("stage_id") or "") != required_stage)):
                    raise ValueError("A collected candidate/model does not have a recorded pass at the selected acceptance stage.")
        available = set(candidates_by_key)
        missing = [(row["candidate_id"], row["model_id"]) for row in selected_candidate_models
                   if (row["candidate_id"], row["model_id"]) not in available]
        if missing:
            raise FileNotFoundError("A passing candidate/model is no longer available in the Results run.")
    result_selections = {}
    for selection_id in requested:
        preview = records_by_id[selection_id]
        if not preview.get("result_only"):
            continue
        candidate_id = str(preview.get("candidate_id") or "")
        model_id = str(preview.get("model_id") or "")
        candidate = candidates_by_key.get((candidate_id, model_id))
        if candidate is None and not model_id:
            candidates = candidates_by_id.get(candidate_id) or []
            candidate = candidates[0] if len(candidates) == 1 else None
        if candidate is None:
            raise FileNotFoundError("The selected Results candidate is no longer available.")
        identity_record = _candidate_for_explicit_selection(candidate, selected)
        if identity_record is None:
            raise ValueError("The selected Results candidate has no valid molecular identity.")
        result_selections[selection_id] = (candidate, identity_record, preview)
    if selected_candidate_models is not None:
        collection_metadata = dict(collection_metadata) if isinstance(collection_metadata, dict) else {}
        result_by_pair = {
            (str(candidate.get("candidate_id") or ""), _candidate_model_id(candidate)): identity
            for candidate, identity, _preview in result_selections.values()
        }
        for row in selected_candidate_models:
            identity = result_by_pair.get((row["candidate_id"], row["model_id"]))
            identity_metadata = (identity or {}).get("metadata") if isinstance(identity, dict) else {}
            identity_metadata = identity_metadata if isinstance(identity_metadata, dict) else {}
            target_sequences = identity_metadata.get("target_sequences_by_chain")
            target_sequences = target_sequences if isinstance(target_sequences, dict) else {}
            target_identity = str(identity_metadata.get("target_identity") or "")
            if target_sequences:
                row["target_sequences_by_chain"] = target_sequences
            if target_identity and target_identity != "Not recorded":
                row["target_identity"] = target_identity
        target_labels = list(dict.fromkeys(str(row.get("target_identity") or "").strip()
                                           for row in selected_candidate_models
                                           if str(row.get("target_identity") or "").strip()
                                           and str(row.get("target_identity")) != "Not recorded"))
        if target_labels and collection_metadata.get("target_identity") in (None, "", "Not recorded"):
            collection_metadata["target_identity"] = target_labels[0] if len(target_labels) == 1 else target_labels
            collection_metadata["target_identities"] = target_labels
        target_maps = {
            _json(row.get("target_sequences_by_chain") or {})
            for row in selected_candidate_models
        }
        if len(target_maps) == 1:
            collection_metadata["target_sequences_by_chain"] = json.loads(next(iter(target_maps)))
        collection_metadata["accepted_candidate_models"] = [
            {**row, "target_identity": row.get("target_identity") or "",
             "target_sequences_by_chain": row.get("target_sequences_by_chain") or {}}
            for row in selected_candidate_models
        ]
    from mn_molecular_design_studio.molecule_schema import validate_record

    requested_project = None if project_id is None else str(project_id or "")
    destination_project = requested_project
    from mn_molecular_design_studio.auth import require_project_access
    if destination_project:
        require_project_access(user, destination_project)
    if group_id:
        group = database.group(group_id, user, write=True)
        if requested_project is not None and str(group["project_id"] or "") != requested_project:
            raise ValueError("The destination group belongs to a different project.")
        require_project_access(user, group["project_id"] or None)
        for identifier in requested:
            record = records_by_id[identifier]
            validate_record({"components": record["components"], "metadata": record.get("metadata") or {}}, group["schema"])
        origin = group["origin"]
        editable_curated = (origin.get("collection_kind") == "curated" or origin.get("type") in {"user", "legacy"}
                            or group["import_method"] in {"user_uploads", "agent_inputs"})
        if not editable_curated:
            raise ValueError("Choose a curated or imported destination group.")
        destination_project = str(group["project_id"] or "")
        if origin.get("collection_kind") != "curated":
            updated_origin = {**origin, "collection_kind": "curated"}
            with database.connection() as connection:
                connection.execute("UPDATE molecule_groups SET origin=?,updated_at=? WHERE id=?",
                                   (_json(updated_origin), _now(), group_id))
    else:
        if destination_project is None:
            destination_project = run_project_id
        require_project_access(user, destination_project or None)
        for identifier in requested:
            record = records_by_id[identifier]
            validate_record({"components": record["components"], "metadata": record.get("metadata") or {}}, {})
        origin = {"type": "user", "user_id": user["id"], "collection_kind": "curated"}
        if isinstance(collection_metadata, dict):
            pairs = sorted({(row["candidate_id"], row["model_id"]) for row in selected_candidate_models or []})
            collection_identity = _json([str(user.get("id") or ""), task_group, run_id,
                                         str(collection_metadata.get("acceptance_stage") or ""), pairs,
                                         collection_metadata.get("acceptance_policies") or []])
            group_id = _id("design-collection-", collection_identity)
            allowed_metadata = {
                key: collection_metadata[key] for key in (
                    "acceptance_stage", "acceptance_stage_id", "acceptance_policies",
                    "accepted_candidate_models", "target_identities", "target_identity", "target_sequences_by_chain",
                    "generators", "generator", "source_tools", "source_workstation", "source_node_id",
                    "generated_candidate_count", "collected_design_count", "collected_model_count",
                ) if key in collection_metadata
            }
            origin.update({
                "type": "results_design_collection",
                "design_collection_version": 1,
                **allowed_metadata,
                "collection_kind": "curated",
                "source_run": {"task_group": task_group, "run_id": run_id,
                               "owner_id": owner_id, "project_id": run_project_id,
                               "source_workstation": str(collection_metadata.get("source_workstation") or socket.gethostname()),
                               "source_node_id": str(collection_metadata.get("source_node_id") or "")},
            })
        group_id = database.create_group(
            user, name, project_id=destination_project, import_method="user_uploads",
            origin=origin, group_id=group_id,
        )

    promoted_ids = []
    result_source_keys = []
    selected_model_by_key = {(row["candidate_id"], row["model_id"]): row
                             for row in selected_candidate_models or []}
    for selection_id, (candidate, identity_record, preview) in result_selections.items():
        candidate_id = str(candidate.get("candidate_id") or "")
        model_id = _candidate_model_id(candidate)
        accepted_row = selected_model_by_key.get((candidate_id, model_id), {})
        generator_engine, prediction_engine = _candidate_engine_provenance(candidate)
        provenance = {
            "type": "explicit_results_selection",
            "saved_selection": True,
            "collection_kind": "curated",
            "owner_id": owner_id,
            "project_id": run_project_id,
            "task_group": task_group,
            "run_id": run_id,
            "candidate_id": candidate_id,
            "model_id": _candidate_model_id(candidate),
            "source_model_id": _candidate_model_id(candidate),
            "source_tool": str(candidate.get("source_tool") or ""),
            "source_stage": str(candidate.get("stage") or ""),
            "generator_engine": generator_engine,
            "prediction_engine": prediction_engine,
            "parent_ids": list(candidate.get("parent_ids") or candidate.get("parents") or []),
            "binder_sequence": str(candidate.get("binder_sequence") or candidate.get("sequence") or ""),
            "binder_chains": list(candidate.get("binder_chains") or []),
            "target_chains": list(candidate.get("target_chains") or []),
            "workstation": str((collection_metadata or {}).get("source_workstation") or socket.gethostname()),
            "source_workstation": str((collection_metadata or {}).get("source_workstation") or socket.gethostname()),
            "source_node": str((collection_metadata or {}).get("source_node_id") or ""),
            "acceptance_status": preview.get("acceptance_status", "unknown"),
            "acceptance_evidence": (candidate.get("raw_metadata") or {}).get("validation_acceptance")
                                   or (candidate.get("raw_metadata") or {}).get("distributed_campaign") or {},
            "native_metrics_preserved": True,
        }
        if accepted_row:
            provenance["acceptance_stage"] = (collection_metadata or {}).get("acceptance_stage")
            provenance["acceptance_decision"] = accepted_row.get("decision") or {}
            provenance["target_identity"] = (accepted_row.get("target_identity")
                or _candidate_target_provenance(candidate)[0])
            provenance["target_sequences_by_chain"] = accepted_row.get("target_sequences_by_chain") or _candidate_target_provenance(candidate)[1]
        source_key = _id("results-selection-", _json([group_id, owner_id, task_group, run_id, candidate_id, model_id]))
        with database.connection() as connection:
            molecule_id = database._put(connection, group_id, identity_record, source_key=source_key,
                                         origin=provenance, stamp=_now(), preserve_existing=True)
        promoted_ids.append(molecule_id)
        result_source_keys.append(source_key)
        with database.connection() as connection:
            structure_ids = [str(row[0]) for row in connection.execute(
                "SELECT id FROM molecule_structures WHERE group_id=? AND molecule_id=? "
                "AND json_extract(context,'$.candidate_id')=? AND COALESCE(json_extract(context,'$.candidate_model_id'),json_extract(context,'$.source_model_id'),json_extract(context,'$.model_id'),json_extract(context,'$.model_index'),json_extract(context,'$.validated_model_id'),'')=?",
                (group_id, molecule_id, candidate_id, model_id),
            ).fetchall()]
        for structure_id in structure_ids:
            database.make_structure_portable(user, structure_id)
        metrics = candidate.get("metrics") if isinstance(candidate.get("metrics"), dict) else {}
        if metrics:
            database.add_observation(
                group_id, molecule_id, str(candidate.get("source_tool") or "results_candidate"), run_id, metrics,
                context={**provenance, "candidate": candidate, "settings": {},
                         "native_metrics": metrics, "is_representative": True},
                source_key=source_key + ":metrics",
            )

    existing_ids = [selection_id for selection_id in requested if not records_by_id[selection_id].get("result_only")]
    if existing_ids:
        database.add_to_group(user, group_id, existing_ids)
    source_keys = []
    source_group_ids = set()
    sources, _ = _run_sources(database, task_group, run_id)
    for source in sources:
        if source["molecule_id"] in existing_ids:
            if selected_candidate_models is not None:
                provenance = json.loads(source["provenance"])
                if not _matches_candidate_models(
                    str(provenance.get("candidate_id") or ""), str(provenance.get("model_id") or ""), selected_candidate_models
                ):
                    continue
            source_keys.append(source["source_key"])
            if source["group_id"] in source_groups:
                source_group_ids.add(source["group_id"])
    # The helper's accessible group map is intentionally limited to this run.
    _copy_run_associations(database, task_group, run_id, source_keys, sorted(source_group_ids), group_id, existing_ids,
                           candidate_models=selected_candidate_models)
    if selected_candidate_models is not None:
        # Attach a stable source/provenance row for every explicitly passing
        # candidate/model, including records that were already indexed before
        # this Results collection was created.
        molecule_by_model = {}
        for selection_id, (candidate, _identity, _preview) in result_selections.items():
            molecule_id = next((value for value in promoted_ids if value == records_by_id[selection_id].get("molecule_id")), "")
            if not molecule_id:
                molecule_id = str(records_by_id[selection_id].get("molecule_id") or "")
            molecule_by_model[(str(candidate.get("candidate_id") or ""), _candidate_model_id(candidate))] = molecule_id
        for selection_id in existing_ids:
            record = records_by_id[selection_id]
            for candidate_id, model_id in _models_on_run_record(record):
                if not _matches_candidate_models(candidate_id, model_id, selected_candidate_models):
                    continue
                molecule_by_model[(candidate_id, model_id)] = str(record.get("molecule_id") or record.get("id") or "")
        with database.connection() as connection:
            for (candidate_id, model_id), accepted_row in selected_model_by_key.items():
                molecule_id = molecule_by_model.get((candidate_id, model_id))
                if not molecule_id:
                    # A candidate may have more than one Results presentation
                    # row but the same canonical identity; resolve by candidate
                    # model within the authorized run records only.
                    molecule_id = next((str(records_by_id[item].get("molecule_id") or records_by_id[item].get("id") or "")
                                        for item in requested if (candidate_id, model_id) in _models_on_run_record(records_by_id[item])), "")
                if not molecule_id:
                    continue
                source_key = _id("results-selection-", _json([group_id, owner_id, task_group, run_id, candidate_id, model_id]))
                candidate = candidates_by_key.get((candidate_id, model_id), {})
                generator_engine, prediction_engine = _candidate_engine_provenance(candidate)
                provenance = {
                    "type": "explicit_results_selection", "saved_selection": True,
                    "collection_kind": "curated", "owner_id": owner_id, "project_id": run_project_id,
                    "task_group": task_group, "run_id": run_id, "candidate_id": candidate_id,
                    "model_id": model_id, "source_model_id": model_id,
                    "source_tool": str(candidate.get("source_tool") or ""),
                    "source_stage": str(candidate.get("stage") or ""),
                    "generator_engine": generator_engine,
                    "prediction_engine": prediction_engine,
                    "parent_ids": list(candidate.get("parent_ids") or candidate.get("parents") or []),
                    "binder_sequence": str(candidate.get("binder_sequence") or candidate.get("sequence") or ""),
                    "binder_chains": list(candidate.get("binder_chains") or []),
                    "target_chains": list(candidate.get("target_chains") or []),
                    "workstation": (collection_metadata or {}).get("source_workstation") or socket.gethostname(),
                    "source_workstation": (collection_metadata or {}).get("source_workstation") or socket.gethostname(),
                    "source_node": str((collection_metadata or {}).get("source_node_id") or ""),
                    "acceptance_stage": (collection_metadata or {}).get("acceptance_stage"),
                    "acceptance_decision": accepted_row.get("decision") or {},
                    "acceptance_evidence": (candidate.get("raw_metadata") or {}).get("validation_acceptance")
                                           or (candidate.get("raw_metadata") or {}).get("distributed_campaign") or {},
                }
                provenance["target_identity"] = (accepted_row.get("target_identity")
                    or _candidate_target_provenance(candidate)[0])
                provenance["target_sequences_by_chain"] = accepted_row.get("target_sequences_by_chain") or _candidate_target_provenance(candidate)[1]
                connection.execute(
                    "INSERT INTO molecule_sources VALUES (?,?,?,?) ON CONFLICT(source_key) DO UPDATE SET "
                    "molecule_id=excluded.molecule_id,group_id=excluded.group_id,provenance=excluded.provenance",
                    (source_key, molecule_id, group_id, _json(provenance)),
                )
    promoted_ids.extend(existing_ids)
    promoted_ids = list(dict.fromkeys(promoted_ids))
    return {"group_id": group_id, "molecule_ids": promoted_ids, "saved_count": len(promoted_ids)}
