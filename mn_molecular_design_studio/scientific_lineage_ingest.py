"""Explicit/completion-time metadata capture. No directory discovery on reads."""
from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
from urllib.parse import urlsplit

from mn_molecular_design_studio import scientific_lineage as lineage


def _read(path):
    if not path.is_file():
        return {}
    value = json.loads(path.read_text())
    if not isinstance(value, dict):
        raise ValueError(f"Expected object: {path.name}")
    return value


def experiment_kind(job):
    text = " ".join(str(job.get(k) or "") for k in ("task_group", "tool", "job_type", "workflow")).lower()
    if any(k in text for k in ("ligand-md", "md-analysis", "mmgbsa", "molecular_dynamics", "simulation")):
        return "dynamics"
    if any(k in text for k in ("docking", "active_learning")):
        return "docking"
    if "analysis" in text:
        return "analysis"
    if any(k in text for k in ("folding", "prediction")):
        return "folding"
    return "design" if "design" in text and "visual_pipeline" not in text else "unknown"


def artifact(value, job, *, run_dir=None, sha256="", verify=False, required=True):
    """Describe a source artifact; never follow a path outside the selected job."""
    original = str(value or "")
    uri = urlsplit(original)
    node = job["node_id"]
    availability, reference, checked = "unresolved", "", False
    path = None
    if uri.scheme in {"runs", "cluster"}:
        parts = uri.path.strip("/").split("/")
        if len(parts) < 3 or any(p in {"", ".", ".."} for p in parts) or uri.query or uri.fragment or (uri.scheme == "runs" and uri.netloc):
            raise ValueError("Invalid source artifact reference.")
        node = uri.netloc if uri.scheme == "cluster" else node
        lineage.job_reference(node, *parts[:2])
        reference = f"cluster://{node}/" + "/".join(parts)
        if run_dir and node == job["node_id"] and parts[:2] == [job["task_group"], job["run_id"]]:
            path = Path(run_dir).joinpath(*parts[2:])
        else:
            availability = "owner_local_unverified"
    elif original and run_dir:
        candidate = Path(original)
        path = candidate if candidate.is_absolute() else Path(run_dir) / candidate
        try:
            relative = path.resolve().relative_to(Path(run_dir).resolve())
            reference = job["ref"] + "/" + relative.as_posix()
        except ValueError:
            path = None  # Preserve provenance; never grant access to an arbitrary file.
    if path is not None:
        try:
            path.resolve().relative_to(Path(run_dir).resolve())
            availability = "owner_present" if path.is_file() else "missing"
            if verify and path.is_file() and path.stat().st_size <= 10 * 1024 * 1024:
                actual = hashlib.sha256(path.read_bytes()).hexdigest()
                if sha256 and actual != sha256:
                    availability = "inaccessible"
                else:
                    sha256, checked = actual, True
        except (OSError, ValueError):
            availability = "inaccessible"
    return {"reference": reference, "original_reference": original, "node_id": node,
            "storage": "results", "availability": availability, "sha256": sha256,
            "digest_verified": checked, "verification_node": node if checked else "", "required": required}


def snapshot_from_candidates(metadata, candidates, *, node_id, run_dir=None, verify=False):
    from mn_molecular_design_studio.complex_analysis_collection import acceptance_decisions
    from mn_molecular_design_studio.output_vocabulary import normalize_metrics
    from mn_molecular_design_studio.molecular_database_results import _candidate_target_provenance
    group = str(metadata.get("task_group") or (Path(run_dir).parent.name if run_dir else ""))
    run = str(metadata.get("run_id") or (Path(run_dir).name if run_dir else ""))
    ref = lineage.job_reference(node_id, group, run)
    job = {k: copy.deepcopy(metadata[k]) for k in (
        "owner_id", "project_id", "cluster_job_id", "created_at", "updated_at", "status",
        "studio_run_id", "studio_run_name", "studio_run_output_id", "studio_run_output_name",
        "parent_role", "pipeline_iteration", "campaign_step", "campaign_step_index", "execution_proxy",
        "distributed_campaign_id", "graph_node_id", "pipeline_parent_node_id",
    ) if k in metadata}
    job.update(ref=ref, node_id=node_id, task_group=group, run_id=run, kind=experiment_kind(metadata),
               owner_id=str(metadata.get("owner_id") or ""), project_id=str(metadata.get("project_id") or ""),
               name=metadata.get("campaign_name") or metadata.get("name") or run,
               workstation=metadata.get("workstation") or node_id, method=str(metadata.get("tool") or ""),
               tool=str(metadata.get("tool") or ""), has_candidate_collection=bool(candidates))
    snapshot = {"schema": lineage.VERSION, "job": job, "outputs": [], "dependencies": [], "counts": [], "warnings": []}
    parent_id = metadata.get("parent_run_id") or metadata.get("campaign_id")
    parent_node = metadata.get("parent_node_id") or metadata.get("pipeline_parent_node_id") or metadata.get("organizer_node_id")
    if parent_id and parent_node:
        # A task group cannot be inferred from a run ID. Resolve via the recorded
        # node/run alias at projection time, or retain an unresolved dependency.
        snapshot["dependencies"].append({"source_ref": ref, "target_ref": f"node-run:{parent_node}:{parent_id}",
            "role": "execution_child_of", "status": "recorded", "evidence": {"source": "job metadata"}})
    elif parent_id:
        snapshot["warnings"].append("Parent run recorded without its source workstation; linkage unresolved.")
    target = metadata.get("target_construct") or {}
    constructs = {lineage.encoded(target["sequences_by_chain"])} if target.get("sequences_by_chain") else set()
    duplicate_keys = {}
    for index, candidate in enumerate(candidates):
        raw = candidate.get("raw_metadata") or {}
        candidate_id = str(candidate.get("candidate_id") or "")
        model = next((str(v) for v in (candidate.get("model_id"), candidate.get("model_index"), raw.get("model_id"), raw.get("model_index")) if v is not None and str(v) != ""), "")
        key = (candidate_id, model)
        duplicate_keys[key] = duplicate_keys.get(key, 0) + 1
        output_ref = ref + "#" + lineage.identifier("output-", candidate_id, model, duplicate_keys[key])
        stage = str(candidate.get("stage") or "")
        kind = ("docking_pose" if job["kind"] == "docking" else "prediction_model" if
                stage in {"complex_refolding", "monomer_refolding"} or job["kind"] == "folding" else
                "sequence_variant" if stage == "sequence_design" else "backbone" if stage == "generation.backbone" else "design")
        info = raw.get("candidate_set") or {}
        metrics = candidate.get("metrics") or {}
        row = {"ref": output_ref, "job_ref": ref, "kind": kind, "candidate_id": candidate_id,
               "native_id": str(raw.get("native_candidate_id") or candidate_id), "model_id": model, "stage": stage,
               "method": str(candidate.get("source_tool") or job["method"]),
               "lineage": {k: info[k] for k in ("target_id", "backbone_id", "variant_id", "prediction_id") if k in info},
               "binder_sequence": candidate.get("binder_sequence"), "sequences": raw.get("chain_sequences") or raw.get("sequences_by_chain") or {},
               "binder_chains": candidate.get("binder_chains") or [], "target_chains": candidate.get("target_chains") or [],
               "metrics": metrics, "normalized_measurements": normalize_metrics(metrics, context=candidate),
               "source_metadata": copy.deepcopy(raw),
               "acceptance": acceptance_decisions(candidate), "artifacts": [], "retained": [],
               "molecular_context": {k: raw[k] for k in ("components", "pose_identity", "ligand_selection", "ligand_smiles",
                    "atom_mapping", "preparation", "receptor_preparation", "compound_id", "target_id", "target_sequences_by_chain") if k in raw}}
        for field in ("complex_pdb", "binder_pdb", "target_pdb"):
            if candidate.get(field):
                row["artifacts"].append(dict(artifact(candidate[field], job, run_dir=run_dir,
                    sha256=str(raw.get("structure_digest") or "") if field == "complex_pdb" else "", verify=verify), role=field))
        for field, role in (("pae_path", "pae_matrix"), ("reference_structure", "reference_structure")):
            if raw.get(field):
                dependency = artifact(raw[field], job, run_dir=run_dir, verify=verify)
                row["artifacts"].append(dict(dependency, role=role))
                snapshot["dependencies"].append({"source_ref": output_ref, "target_ref": dependency["reference"] or dependency["original_reference"],
                    "role": role, "status": "recorded" if dependency["reference"] else "unresolved", "evidence": {"field": field}})
        for parent in candidate.get("parents") or []:
            snapshot["dependencies"].append({"source_ref": output_ref, "target_ref": str(parent), "role": "derived_from",
                "status": "recorded" if str(parent).startswith(("cluster:", "runs:", "database:")) else "unresolved",
                "evidence": {"native_parent_id": parent}})
        snapshot["outputs"].append(row)
        target_identity, sequences = _candidate_target_provenance(candidate)
        sequences = raw.get("target_sequences_by_chain") or sequences
        if sequences:
            constructs.add(lineage.encoded(sequences))
            recorded_target = {"sequences_by_chain": sequences}
            if not target:
                target = recorded_target
            elif target.get("sequences_by_chain") and target["sequences_by_chain"] != sequences:
                snapshot["warnings"].append("Multiple target constructs in this job; no cross-construct metric pooling is permitted.")
            row["target_construct"] = dict(recorded_target, identity=target_identity)
        if duplicate_keys[key] > 1:
            snapshot["warnings"].append(f"Repeated candidate/model key {candidate_id}/{model}; occurrences retained separately.")
    snapshot["target"] = copy.deepcopy(target)
    if len(constructs) > 1:
        snapshot["target"] = {"mapping_status": "multiple_constructs", "recorded_constructs": [json.loads(s) for s in sorted(constructs)]}
    from mn_molecular_design_studio.candidate_sets import json_safe
    return json_safe(snapshot)


def build_snapshot(run_dir, *, node_id, metadata=None, verify=False):
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
    root = Path(run_dir)
    meta = {**_read(root / "metadata.json"), **(metadata or {})}
    marker = _read(root / "molecular_design_studio.json")
    meta.setdefault("owner_id", marker.get("studio_owner_id") or "")
    meta.setdefault("project_id", marker.get("project_id") or "")
    candidates = read_candidates(root)
    snapshot = snapshot_from_candidates(meta, candidates, node_id=node_id, run_dir=root, verify=verify)
    result = _read(root / "result.json")
    # Retain method-native count definitions; never sum native variants and
    # normalized candidate rows or relabel native acceptance as final.
    for key, value in (result.get("metrics") or {}).items():
        if isinstance(value, (int, float)) and not isinstance(value, bool) and any(k in key for k in ("candidates", "count", "attempts", "scored", "replicas")):
            snapshot["counts"].append({"name": key, "value": value, "unit": "native " + key,
                "stage": snapshot["job"]["method"], "source": "result.json/metrics", "scope": "producer_reported"})
    snapshot["counts"].append({"name": "normalized_output_occurrences", "value": len(candidates),
        "unit": "candidate/model rows", "stage": "saved output", "source": "candidates.jsonl", "scope": "this_job"})
    ledger = _read(root / "artifacts/distributed_campaign/ledger.json")
    for attempt in (ledger.get("attempts") or {}).values():
        if attempt.get("cluster_job_id"):
            snapshot["dependencies"].append({"source_ref": "cluster-job:" + attempt["cluster_job_id"],
                "target_ref": snapshot["job"]["ref"], "role": "execution_child_of", "status": "recorded",
                "evidence": {"source": "distributed campaign ledger", "attempt_id": attempt.get("attempt_id")}})
    study = _read(root / "artifacts/md/study.json")
    if study:
        snapshot["job"]["kind"] = "dynamics"
        for item in study.get("studies", [study]):
            add_md_study(snapshot, item)
    return snapshot


def add_md_study(snapshot, study):
    """Index MD manifests, not trajectories; preserve continuation and GB/PB subjects."""
    job = snapshot["job"]
    study_ref = job["ref"] + "#study:" + str(study.get("study_id") or "unresolved")
    snapshot["dependencies"].append({"source_ref": study_ref,
        "target_ref": str(study.get("source_pose_id") or "unresolved-pose"), "role": "selected_pose_input",
        "status": "unresolved", "evidence": {k: study[k] for k in ("molecule_id", "source_pose_id", "target_provenance", "ligand_smiles", "selection_status") if k in study}})
    replicas = copy.deepcopy(study.get("replicas") or [])
    if not replicas:
        by_replica = {}
        for production in study.get("production_runs") or []:
            # repeat_index is an independent replica; continuation is a segment.
            repeat = production.get("repeat_index")
            key = str(repeat) if repeat is not None else "unknown:" + str(production.get("run_id") or "")
            by_replica.setdefault(key, {"replica_id": key, "segments": []})["segments"].append(production)
        replicas = list(by_replica.values())
    for replica in replicas:
        rid = str(replica.get("replica_id") or replica.get("id") or "")
        ref = job["ref"] + "#" + lineage.identifier("replica-", study_ref, rid)
        snapshot["outputs"].append({"ref": ref, "kind": "md_replica", "candidate_id": rid,
            "model_id": "", "stage": "dynamics", "study_ref": study_ref, "replica": copy.deepcopy(replica), "artifacts": [], "retained": []})
        for segment in replica.get("segments") or replica.get("production_runs") or []:
            sid = str(segment.get("segment_id") or segment.get("run_id") or "")
            segment_ref = job["ref"] + "#" + lineage.identifier("segment-", rid, sid)
            snapshot["outputs"].append({"ref": segment_ref, "kind": "md_segment", "candidate_id": sid,
                "model_id": "", "stage": "production", "replica_ref": ref, "segment": copy.deepcopy(segment), "retained": [], "artifacts": []})
            snapshot["dependencies"].append({"source_ref": segment_ref, "target_ref": ref, "role": "segment_of", "status": "recorded"})
            if segment.get("trajectory_path"):
                snapshot["outputs"][-1]["artifacts"] = [dict(artifact(segment["trajectory_path"], job), role="trajectory")]
            if segment.get("continuation_of_run_id"):
                snapshot["dependencies"].append({"source_ref": segment_ref,
                    "target_ref": job["ref"] + "#" + lineage.identifier("segment-", rid, str(segment["continuation_of_run_id"])),
                    "role": "continues", "status": "recorded", "evidence": {"source": "MD study manifest"}})
    for entry in [*(study.get("mmgbsa_runs") or []), *(study.get("analysis_runs") or [])]:
        snapshot["dependencies"].append({"source_ref": str(entry.get("run_id") or ""),
            "target_ref": str(entry.get("source_trajectory_ref") or study_ref), "role": "trajectory_analysis_input",
            "status": "recorded" if entry.get("source_trajectory_ref") else "unresolved", "evidence": copy.deepcopy(entry)})


def link_retained(database, snapshot, *, local_node=None):
    """Match explicit source/candidate/model evidence; never promote or match by hash alone."""
    job = snapshot["job"]
    snapshot["outputs"] = [o for o in snapshot["outputs"] if o.get("kind") not in {"evaluation", "retained_representation"}]
    for output in snapshot["outputs"]:
        output["retained"] = []
    if local_node is None:
        from mn_molecular_design_studio.cluster import config
        local_node = config.identity()
    with database.connection() as c:
        target_evidence = set()
        structures = c.execute("SELECT s.*,g.origin FROM molecule_structures s JOIN molecule_groups g ON g.id=s.group_id "
            "WHERE json_extract(s.context,'$.run_id')=? AND json_extract(s.context,'$.task_group')=? AND g.owner_id=? AND g.project_id=?",
            (job["run_id"], job["task_group"], job["owner_id"], job["project_id"])).fetchall()
        for row in snapshot["outputs"]:
            for s in structures:
                context, origin = json.loads(s["context"]), json.loads(s["origin"])
                source_node = context.get("_replica_node") or origin.get("_replica_node")
                if (source_node or local_node) != job["node_id"]:
                    continue
                scientific_target = {k: context[k] for k in ("target_molecule_id", "design_target_molecule_id", "design_target_sequence_sha256", "target_structure_digest", "target_variant_sha256", "target_receptor", "target_pdb_id", "target_id") if context.get(k)}
                if scientific_target:
                    target_evidence.add(lineage.encoded(scientific_target))
                ids = {row.get("candidate_id"), row.get("native_id")}
                if not context.get("candidate_id") or context["candidate_id"] not in ids:
                    continue
                model = next((str(v) for v in (context.get("model_id"), context.get("model_index")) if v is not None and str(v) != ""), "")
                if model != str(row.get("model_id") or ""):
                    continue
                digest = context.get("portable_sha256") or context.get("content_sha256") or ""
                known = {a["sha256"] for a in row.get("artifacts", []) if a.get("sha256") and a.get("role") != "target_pdb"}
                if known and digest not in known:
                    continue
                row["retained"].append({"molecule_id": s["molecule_id"], "structure_id": s["id"], "group_id": s["group_id"],
                    "sha256": digest, "source_node": job["node_id"], "local_reference": not bool(source_node),
                    "owner_structure_id": context.get("_replica_structure") or s["id"]})
            if not row.get("retained") and row.get("acceptance"):
                snapshot["warnings"].append(f"No exact retained structure mapping for {row.get('candidate_id')}/{row.get('model_id', '')}.")
        mapped = {r["structure_id"] for o in snapshot["outputs"] for r in o.get("retained", [])}
        for s in structures:
            if s["id"] in mapped:
                continue
            context, origin = json.loads(s["context"]), json.loads(s["origin"])
            if (context.get("_replica_node") or origin.get("_replica_node") or local_node) != job["node_id"]:
                continue
            # A retained reference is evidence in its own right. It does not
            # establish which unknown prediction sample generated its bytes.
            source_id = context.get("_replica_structure") or s["id"]
            digest = context.get("portable_sha256") or context.get("content_sha256") or ""
            snapshot["outputs"].append({"ref": job["ref"] + "#representation:" + source_id,
                "kind": "retained_representation", "candidate_id": str(context.get("candidate_id") or ""),
                "model_id": next((str(v) for v in (context.get("model_id"), context.get("model_index")) if v is not None and str(v) != ""), ""), "stage": str(context.get("stage") or ""),
                "source_context": context, "artifacts": [], "acceptance": [],
                "retained": [{"molecule_id": s["molecule_id"], "structure_id": s["id"], "group_id": s["group_id"],
                              "sha256": digest, "source_node": job["node_id"], "owner_structure_id": source_id}]})
            snapshot["warnings"].append(f"Retained representation {source_id} preserved; exact output/model association unresolved.")
        observations = c.execute("SELECT e.descriptor,o.* FROM molecule_observations o LEFT JOIN molecule_evaluations e ON o.id=e.observation_id "
            "JOIN molecule_groups g ON g.id=o.group_id WHERE o.run_id=? AND g.owner_id=? AND g.project_id=?",
            (job["run_id"], job["owner_id"], job["project_id"])).fetchall()
        for observation in observations:
            context = json.loads(observation["context"])
            if (context.get("_replica_node") or local_node) != job["node_id"]:
                continue
            from mn_molecular_design_studio.evaluations import evaluation_descriptor
            descriptor = (json.loads(observation["descriptor"]) if observation["descriptor"] else
                          evaluation_descriptor(dict(observation, context=context)))
            scientific_target = {k: context[k] for k in ("target_molecule_id", "target_structure_digest", "target_variant_sha256", "target_receptor", "target_pdb_id", "target_id") if context.get(k)}
            if scientific_target:
                target_evidence.add(lineage.encoded(scientific_target))
            ref = job["ref"] + "#evaluation:" + descriptor["id"]
            snapshot["outputs"].append({"ref": ref, "kind": "evaluation", "candidate_id": descriptor.get("candidate_id", ""),
                "model_id": descriptor.get("model_id", ""), "stage": descriptor.get("stage", ""),
                "evaluation_id": descriptor["id"], "observation_id": observation["id"],
                "molecule_id": observation["molecule_id"], "descriptor": descriptor,
                "raw_measurements": json.loads(observation["metrics"]), "units": copy.deepcopy(context.get("units")),
                "method_context": {k: copy.deepcopy(context[k]) for k in ("settings", "source_units", "source_summary_relative_path", "source_summary_sha256", "target_provenance", "selection", "method_version") if k in context},
                "source_context": context,
                "artifacts": [], "retained": []})
            inputs = list(context.get("analysis_inputs") or [])
            if context.get("source_trajectory_ref"):
                inputs.append({"reference": context["source_trajectory_ref"], "role": "trajectory_snapshot",
                    "sha256": context.get("recorded_trajectory_sha256") or "",
                    "snapshot_matches": context.get("trajectory_snapshot_matches")})
            representation = context.get("representation_id") or descriptor.get("representation_id")
            if representation:
                inputs.append({"reference": "database:" + observation["molecule_id"] + ":" + str(representation), "role": "analyzed_structure"})
            if not inputs:
                inputs = [{"reference": "", "role": "analyzed_input", "sha256": descriptor.get("structure_digest", "")}]
            for entry in inputs:
                snapshot["dependencies"].append({"source_ref": ref, "target_ref": str(entry.get("reference") or "unresolved-input"),
                    "role": str(entry.get("role") or "scientific_input"), "status": "recorded" if entry.get("reference") else "unresolved",
                    "evidence": copy.deepcopy(entry)})
    if not snapshot.get("target") and len(target_evidence) == 1:
        snapshot["target"] = json.loads(next(iter(target_evidence)))
    elif not snapshot.get("target") and len(target_evidence) > 1:
        snapshot["warnings"].append("Recorded target contexts differ; construct/receptor correspondence unresolved.")
    snapshot["warnings"] = list(dict.fromkeys(snapshot["warnings"]))
    return snapshot


def capture_completed_run(run_dir, *, database=None, node_id=None):
    from mn_molecular_design_studio.molecular_database import MolecularDatabase
    from mn_molecular_design_studio.cluster import config
    database = database or MolecularDatabase()
    snapshot = build_snapshot(run_dir, node_id=node_id or config.identity(), verify=True)
    link_retained(database, snapshot)
    lineage.put_snapshot(database, snapshot, author=snapshot["job"]["node_id"])
    return snapshot


def chaperone_records(metadata, rows, *, node_id):
    """v1 native-schema translation; no Chaperone/Polars import in Python 3.11.

    Ordered structural history is explicit. Its elements are representations,
    not independent predictions unless upstream metadata records that meaning.
    """
    candidates = []
    for row in rows:
        if not row.get("design_id") or row.get("artifact_kind") not in {"sequence", "structure"}:
            raise ValueError("Unsupported Chaperone row; design_id and artifact_kind required.")
        sequences = row.get("sequences") or {}
        metrics = row.get("metrics") or {}
        if isinstance(sequences, str):
            sequences = json.loads(sequences)
        if isinstance(metrics, str):
            metrics = json.loads(metrics)
        candidates.append({"candidate_id": row["design_id"], "stage": "sequence_design" if row["artifact_kind"] == "sequence" else "generation.backbone_sequence",
            "source_tool": metadata.get("tool", "chaperone"), "metrics": metrics,
            "parents": [row["parent_design_id"]] if row.get("parent_design_id") else [],
            "raw_metadata": {"chain_sequences": sequences}})
    result = snapshot_from_candidates(metadata, candidates, node_id=node_id)
    parent_refs = {}
    for output in result["outputs"]:
        parent_refs.setdefault(output["candidate_id"], []).append(output["ref"])
    for edge in result["dependencies"]:
        choices = parent_refs.get(edge.get("evidence", {}).get("native_parent_id"), [])
        if len(choices) == 1:
            edge.update(target_ref=choices[0], status="recorded")
    for output, row in zip(result["outputs"], rows):
        output["chaperone"] = {"adapter": "studio.chaperone-lineage.v1", "artifact_kind": row["artifact_kind"],
                               "structure_history": copy.deepcopy(row.get("structures") or [])}
        for index, structure in enumerate(row.get("structures") or []):
            path = structure.get("store_path") or structure.get("path") or ""
            output["artifacts"].append({"reference": path, "node_id": node_id, "storage": "results",
                "availability": "unresolved", "sha256": "", "history_index": index,
                "step": structure.get("name"), "path_handle": structure.get("path_handle"),
                "storage_format": structure.get("storage_format"), "metadata": copy.deepcopy(structure.get("metadata") or {})})
    return result


def main():
    """Controlled, explicit historical manifest; no recursive Results scan."""
    import argparse
    from mn_molecular_design_studio.molecular_database import MolecularDatabase
    from mn_molecular_design_studio.runtime import app_home
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--staging-database-root", required=True)
    parser.add_argument("--manifest", required=True, help="JSON list of source node, metadata, and optional exact run_dir")
    args = parser.parse_args()
    destination = Path(args.staging_database_root).resolve()
    if destination == (app_home() / "database").resolve() or destination == (Path.home() / ".local/share/mn-molecular-design-studio/database").resolve():
        parser.error("This backfill command is staging-only; production destination refused.")
    marker = destination / "scientific-lineage-staging.json"
    if (destination / "molecular_library.sqlite3").exists() and not marker.is_file():
        parser.error("An existing destination requires an explicit scientific-lineage-staging.json marker.")
    destination.mkdir(parents=True, exist_ok=True)
    marker.write_text(json.dumps({"purpose": "isolated scientific lineage validation", "schema": lineage.VERSION}))
    database = MolecularDatabase(destination)
    entries = json.loads(Path(args.manifest).read_text())
    reports = []
    for entry in entries:
        node = entry["node_id"]
        snapshot = (build_snapshot(entry["run_dir"], node_id=node, metadata=entry.get("metadata"), verify=True)
                    if entry.get("run_dir") else snapshot_from_candidates(entry["metadata"], entry.get("candidates", []), node_id=node))
        link_retained(database, snapshot, local_node=entry.get("local_node_id", node))
        lineage.put_snapshot(database, snapshot, author=node, rebuild=False)
        reports.append({"job": snapshot["job"]["ref"], "outputs": len(snapshot["outputs"]), "warnings": snapshot["warnings"]})
    with database.connection() as c:
        lineage.rebuild_campaigns(c)
    print(json.dumps({"schema": lineage.VERSION, "jobs": reports}, indent=2))


if __name__ == "__main__":
    main()
