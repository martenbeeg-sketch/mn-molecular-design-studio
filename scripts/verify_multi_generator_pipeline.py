"""Replay retained PD-L1 evidence through the actual Studio graph, CPU only.

No scientific rows, sequences or structures are fabricated. Imported source
snapshots preserve absolute artifact references, original run IDs and hashes.
Fixtures and live submission workers are not used in this retained-evidence run.
"""
from __future__ import annotations

import argparse
import json
import os
import hashlib
from pathlib import Path
from zipfile import ZipFile

from mn_molecular_design_studio import candidate_pipeline as cp, candidate_sets as cs
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import create_job, finish_job, read_json, write_json
from mn_molecular_design_studio.engine.workflow_engine.core.artifacts import build_run_zip
from mn_molecular_design_studio.collection_policy import collection_scope


def snapshot(sources, tool, project, owner):
    from mn_molecular_design_studio.workspace import mark_studio_job
    job = create_job("design", "retained_evidence_import", tool, {}, {"sources": [str(p) for p in sources]})
    destination = job.run_dir / "artifacts/normalized_candidates/candidates.jsonl"
    destination.parent.mkdir(parents=True, exist_ok=True)
    origins, counts = [], {}
    saved_root = os.environ.get("MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR")
    try:
        with destination.open("w") as handle:
            for source in sources:
                os.environ["MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR"] = str(source.parent.parent)
                provenance = dict(run_id=source.name, run_dir=str(source), input=read_json(source / "input.json"),
                                  metadata=read_json(source / "metadata.json"), result_metrics=read_json(source / "result.json").get("metrics", {}),
                                  hashes={str(p): cs.file_digest(p) for name in ("input.json", "metadata.json", "result.json", "artifacts/distributed_campaign/ledger.json", "artifacts/normalized_candidates/candidates.jsonl") if (p := source / name).is_file()})
                ledger = read_json(source / "artifacts/distributed_campaign/ledger.json")
                if ledger:
                    provenance['native_attempt_records'] = ledger.get('attempts', {})
                    rejection_count = sum(sum(int(value) for value in attempt.get('native_state', {}).get('rejections', {}).get('terminated', {}).values())
                                          for attempt in ledger.get('attempts', {}).values())
                    if rejection_count:
                        provenance['result_metrics']['native_rejected_attempts'] = rejection_count
                origins.append(provenance)
                for key, value in cs.native_counts(provenance).items():
                    counts[key] = counts.get(key, 0) + value
                for row in cs.iter_candidates(dict(run_dir=str(source), candidates_jsonl=str(source / "artifacts/normalized_candidates/candidates.jsonl"))):
                    row["raw_metadata"]["retained_origin"] = provenance
                    handle.write(json.dumps(row, sort_keys=True) + "\n")
    finally:
        if saved_root is None:
            os.environ.pop("MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR", None)
        else:
            os.environ["MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR"] = saved_root
    write_json(job.run_dir / "artifacts/retained_sources.json", dict(origins=origins, stage_counts=counts))
    finish_job(job.run_dir, True, {"outputs": {"candidates_jsonl": str(destination)}, "metrics": counts})
    mark_studio_job({"job": {"run_dir": str(job.run_dir)}}, project_id=project, owner_id=owner)
    return job.run_dir


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", required=True)
    parser.add_argument("--owner", required=True)
    parser.add_argument("--rfd4-run", type=Path, action="append", required=True)
    parser.add_argument("--bindcraft-run", type=Path, required=True)
    parser.add_argument("--target", type=Path, required=True)
    parser.add_argument("--evidence", type=Path, required=True)
    args = parser.parse_args()
    originals = [*args.rfd4_run, args.bindcraft_run]
    before = {str(p): cs.file_digest(p) for root in originals for p in root.rglob("*") if p.is_file() and p.suffix in {".jsonl", ".pdb", ".cif"} or p.is_file() and p.name == "ledger.json"}
    with collection_scope("manual"):
        rfd4 = snapshot(args.rfd4_run, "rfd4_proteina", args.project, args.owner)
        bindcraft = snapshot([args.bindcraft_run], "bindcraft2", args.project, args.owner)
        graph = cp.new_pipeline()
        for node in graph["nodes"]:
            if node["id"] == "target":
                node["data"]["params"].update(target_pdb=str(args.target), target_chains="A", hotspots="A40,A99,A107")
            if node["id"] == "generator-0":
                node["data"]["retained_run_dir"] = str(rfd4)
            if node["id"] == "generator-1":
                node["data"]["retained_run_dir"] = str(bindcraft)
            if node["id"] == "predictor-0":
                node["data"]["retained_run_dir"] = str(rfd4)
            if node["id"] == "ranking":
                node["data"]["params"]["ranking_rules"] = [{"field": "prediction/predictor-0/iptm", "direction": "higher"}]
        from mn_molecular_design_studio.pipeline_graph import clean_graph, save_definition
        graph = clean_graph(graph)
        save_definition(graph, args.owner, args.project, version=True)
        from mn_molecular_design_studio.engine.workflow_engine.workflows.graph_pipeline import run_graph_pipeline
        run = run_graph_pipeline(graph, "Phase 7A retained PD-L1 graph integration", args.owner, args.project)
    result = read_json(run / "result.json")
    rows = list(cs.iter_candidates(dict(run_dir=str(run), candidates_jsonl=str(run / "artifacts/normalized_candidates/candidates.jsonl"))))
    assert result["success"] is True, result
    assert sorted(len(row["binder_sequence"]) for row in rows) == [75, 94, 116]
    assert all(len(cs.metadata(row)["predictions"]) == 1 for row in rows)
    assert all(cs.metadata(row)["predictor_matrix"]["predictor-1"] == ["disabled"] for row in rows)
    after = {p: cs.file_digest(p) for p in before}
    assert before == after, "Original retained scientific artifacts changed"
    archive = run / "phase7a-evidence.zip"
    archive.write_bytes(build_run_zip(run))
    with ZipFile(archive) as bundle:
        assert "artifacts/candidate_set/predictions.jsonl" in bundle.namelist()
    evidence = dict(schema="studio.phase7a-verification.v1", run_dir=str(run), graph_id=graph["id"], project_id=args.project,
        scientific_execution="retained outputs only; zero new inference and zero scheduler submissions",
        real_predictor="retained Boltz-2", additional_predictor="disabled here; deterministic adapter worker tested separately",
        source_hashes_unchanged=before == after, source_files=before,
        candidates=[dict(candidate_id=row["candidate_id"], generator=cs.metadata(row)["generator"], source_run=cs.metadata(row)["source_run"],
            binder_length=len(row["binder_sequence"]), sequence_sha256=hashlib.sha256(row["binder_sequence"].encode()).hexdigest(),
            variant_id=cs.metadata(row)["variant_id"], predictor_matrix=cs.metadata(row)["predictor_matrix"],
            iptm=next(iter(cs.metadata(row)["predictions"]))["candidate"]["metrics"].get("iptm")) for row in rows],
        stage_counts={key: value.get("stage_counts", {}) for key, value in read_json(run / "graph_execution.json")["nodes"].items()},
        export_archive=str(archive), export_sha256=cs.file_digest(archive))
    args.evidence.parent.mkdir(parents=True, exist_ok=True)
    args.evidence.write_text(json.dumps(evidence, indent=2) + "\n")
    print(json.dumps({"run_dir": str(run), "evidence": str(args.evidence), "candidate_count": len(rows)}, indent=2))


if __name__ == "__main__":
    main()
