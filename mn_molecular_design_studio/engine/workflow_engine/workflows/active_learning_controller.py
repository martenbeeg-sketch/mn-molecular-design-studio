"""Resumable active-learning docking campaign for the visual pipeline canvas."""
from __future__ import annotations

import csv
import heapq
import json
import math
import os
import shlex
import shutil
import signal
import sqlite3
import subprocess
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from mn_molecular_design_studio.active_learning_settings import ActiveLearningSettings
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import (
    create_job, finish_job, read_json, update_status, utc_now, write_json,
)
from mn_molecular_design_studio.pipeline_graph import validate_active_library_sources, validate_graph


_DOCKING_WORKER_SCHEMA = "studio.active-learning-docking-shard.v1"
_PACKED_FINGERPRINT_BYTES = 128


def _node_by_kind(graph: dict) -> dict[str, dict]:
    return {str(node["data"]["kind"]): node for node in graph["nodes"]}


def _next_round_selection_counts(*, initial_seed_count: int, remaining_compounds: int,
                                 max_batch_size: int, seed_fraction_percent: float = 50.0,
                                 predicted_fraction_percent: float = 80.0) -> dict[str, int]:
    """Size a later docking batch relative to the actual initial centroid sample."""
    if initial_seed_count < 1:
        raise ValueError("The initial centroid sample must contain at least one compound.")
    if not 1 <= float(seed_fraction_percent) <= 100:
        raise ValueError("The next-round fraction of the initial centroid sample must be between 1 and 100 percent.")
    if not 0 <= float(predicted_fraction_percent) <= 100:
        raise ValueError("The predicted-compound fraction must be between 0 and 100 percent.")
    if max_batch_size < 0:
        raise ValueError("The maximum docking batch size cannot be negative.")

    target_count = max(1, int(math.floor(initial_seed_count * float(seed_fraction_percent) / 100.0 + 0.5)))
    total = min(target_count, max(0, int(remaining_compounds)))
    if max_batch_size:
        total = min(total, int(max_batch_size))
    predicted = int(math.floor(total * float(predicted_fraction_percent) / 100.0 + 0.5))
    return {"total": total, "predicted": predicted, "random": total - predicted,
            "target_total": target_count}


def _compound_batches(rows: list[dict], batch_size: int) -> list[list[dict]]:
    """Split one computer's assignments into launches accepted by the docking adapter."""
    if batch_size < 1:
        raise ValueError("A docking launch must contain at least one compound.")
    return [rows[start:start + batch_size] for start in range(0, len(rows), batch_size)]


def _write_graph_state(path: Path, state: dict, *, node: str | None = None,
                       status: str | None = None, **values) -> None:
    if node is not None:
        item = state.setdefault("nodes", {}).setdefault(node, {})
        if status is not None:
            item["status"] = status
        item.update(values)
        item["updated_at"] = utc_now()
    state["updated_at"] = utc_now()
    write_json(path, state)


def _library_sources(params: dict, seed_params: dict):
    from mn_molecular_design_studio.external_libraries import LibraryCatalog, read_feature_manifest
    from mn_molecular_design_studio.engine.workflow_engine.workflows.library_preparation import _centroid_seed_inputs

    features, centroids = validate_active_library_sources(params, seed_params)
    feature_path, feature_manifest = read_feature_manifest(params["feature_manifest"])
    seed_record, seed_manifest, seed_sample = _centroid_seed_inputs(
        str(seed_params["library_id"]), Path(seed_params["centroid_manifest"]),
        sample_size=int(seed_params.get("sample_size") or 10_000),
        sampling_seed=int(seed_params.get("sampling_seed", 42)),
    )
    screening_record = LibraryCatalog().get(str(params["library_id"]))
    if feature_manifest.get("library_id") != screening_record.get("id"):
        raise ValueError("The selected screening features no longer match the registered library.")
    return {
        "feature_path": feature_path,
        "feature_manifest_path": Path(params["feature_manifest"]).resolve(strict=True),
        "feature_manifest": feature_manifest,
        "screening_record": screening_record,
        "centroid_manifest": Path(seed_params["centroid_manifest"]).resolve(strict=True),
        "centroids": centroids,
        "seed_record": seed_record,
        "seed_manifest": seed_manifest,
        "seed_sample": seed_sample,
    }


def _selected_worker_specs(node: dict) -> list[dict]:
    from mn_molecular_design_studio.cluster import config as cluster_config
    from mn_molecular_design_studio.cluster import distributed_library
    from mn_molecular_design_studio.cluster.store import Store

    raw_specs = node.get("worker_nodes")
    if not isinstance(raw_specs, list) or not raw_specs:
        raise ValueError("Select at least one online computer in the Docking step.")
    specs = distributed_library._node_specs(raw_specs)
    ssh_user = str(node.get("ssh_user") or "compute")
    if not distributed_library._SSH_USER.fullmatch(ssh_user):
        raise ValueError("The SSH user must be a valid account name such as compute.")
    by_id = {item["id"]: item for item in Store().nodes()}
    for spec, configured in zip(specs, raw_specs):
        machine = by_id[spec["node_id"]]
        resources = (machine.get("snapshot") or {}).get("resources") or {}
        available_cpu = int(resources.get("cpu_free") or 0)
        if spec["local"] and spec["cpu_cores"] > available_cpu:
            raise ValueError(
                f"The organizer has {available_cpu:,} CPU slots free after reserving the campaign coordinator, "
                f"but its docking worker requests {spec['cpu_cores']:,}. Lower the organizer's worker count or wait for other jobs to finish."
            )
        gpu_id = str(configured.get("gpu_device") or "none")
        if gpu_id.lower() in {"none", "cpu", "off"}:
            raise ValueError(f"UniDock Pro requires a selected GPU on {spec['name']}.")
        gpu_records = resources.get("free_gpus") if "free_gpus" in resources else resources.get("gpus", [])
        free = {str(item.get("gpu_id")) for item in (gpu_records or [])}
        if gpu_id not in free:
            raise ValueError(f"GPU {gpu_id} on {spec['name']} is no longer available.")
        spec.update(gpu_device=gpu_id, node_name=spec["name"], local=spec["node_id"] == cluster_config.identity())
        if not spec["local"]:
            from mn_molecular_design_studio.cluster.distributed_library import _remote_environment
            remote_python, _repo, _state, _slots = _remote_environment(spec, str(node.get("ssh_user") or "compute"))
            spec["remote_python"] = remote_python
    if not any(spec["local"] for spec in specs):
        raise ValueError("Include the organizer computer; it owns the screening fingerprint files and campaign results.")
    return specs


def _docking_params(params: dict, target: dict, batch_file: Path, *, round_number: int,
                   worker: dict, limit: int) -> dict:
    if str(params.get("engine") or "UniDock Pro") != "UniDock Pro":
        raise ValueError("This first active-learning controller currently runs UniDock Pro. Other docking adapters will be connected next.")
    if str(params.get("score_field") or "vina_score") not in {"vina_score", "best_score_kcal_mol"}:
        raise ValueError("UniDock Pro exposes the primary score as best_score_kcal_mol (vina_score alias).")
    if str(params.get("score_direction") or "Lower is better") != "Lower is better":
        raise ValueError("UniDock Pro Vina scores rank lower values first.")
    alignment_enabled = bool(target.get("align_longest_axis_to_x", False))
    alignment_ligand = str(target.get("alignment_ligand_sdf") or "")
    if alignment_enabled and (not alignment_ligand or not Path(alignment_ligand).is_file()):
        raise ValueError("The selected receptor requests axis alignment, but its reference ligand SDF is missing.")
    return {
        "target_pdb": str(target["target_pdb"]),
        "compound_set_file": str(batch_file),
        "center_x": float(target["center_x"]), "center_y": float(target["center_y"]), "center_z": float(target["center_z"]),
        "size_x": float(target["size_x"]), "size_y": float(target["size_y"]), "size_z": float(target["size_z"]),
        "search_mode": str(params.get("search_mode") or "balance"),
        "num_modes": int(params.get("num_modes") or 1),
        "seed": (42 + int(round_number)) % 2_147_483_646 or 42,
        "compound_limit": int(limit),
        "cpu_cores": int(worker["cpu_cores"]),
        "gpu_device": str(worker["gpu_device"]),
        "job_name": f"Active learning round {round_number:03d} · {worker['node_name']}",
        "align_longest_axis_to_x": alignment_enabled,
        "alignment_ligand_sdf": alignment_ligand if alignment_enabled else None,
        # The dimensions configured in the receptor node describe the aligned
        # frame (for example the requested 30 × 20 × 20 Å THRα box).
        "preserve_box_size_after_alignment": alignment_enabled,
    }


def _remote_request_rewrite_command(python: str, request_path: str, remote_dir: str,
                                    alignment_path: str = "") -> str:
    """Build the remote request rewrite with the worker's configured Python."""
    script = (
        "import json,sys; p=sys.argv[1]; d=json.load(open(p)); "
        "d['parameters']['target_pdb']=sys.argv[2]; "
        "d['parameters']['compound_set_file']=sys.argv[3]; "
        "(d['parameters'].__setitem__('alignment_ligand_sdf',sys.argv[4]) "
        "if len(sys.argv)>4 and sys.argv[4] else None); json.dump(d,open(p,'w'))"
    )
    return shlex.join([
        str(python), "-c", script, request_path,
        remote_dir + "/receptor.pdb", remote_dir + "/compounds.json", alignment_path,
    ])


def _wait_local_job(run_dir: Path, *, cancel_check) -> dict:
    cancellation_sent = False
    while True:
        if cancel_check() and not cancellation_sent:
            from mn_molecular_design_studio.engine.workflow_engine.core.jobs import stop_job
            stop_job(run_dir, reason="Active-learning campaign was cancelled")
            cancellation_sent = True
        metadata = read_json(run_dir / "metadata.json")
        status = str(metadata.get("status") or "queued")
        if status in {"completed", "failed", "cancelled"}:
            return {"status": status, "run_dir": str(run_dir), "job_code": str(metadata.get("job_code") or ""),
                    "error": str(read_json(run_dir / "result.json").get("error") or metadata.get("worker_error") or "")}
        time.sleep(2)


def _dispatch_docking_shard(*, request: dict, worker: dict, campaign: Path,
                            ssh_user: str, cancel_check) -> dict:
    from mn_molecular_design_studio.engine.workflow_engine.core.workflow_queue import enqueue_workflow_call
    from mn_molecular_design_studio.engine.workflow_engine.workflows.ligand_tools import run_ligand_unidock_pro
    from mn_molecular_design_studio.workspace import mark_studio_job

    if worker["local"]:
        parameters = {key: value for key, value in request["parameters"].items()}
        run_dir = enqueue_workflow_call(
            run_ligand_unidock_pro, task_group="ligand-docking", tool="unidock_pro",
            job_type="ligand_docking", kwargs=parameters,
            cpu_cores=int(parameters["cpu_cores"]), gpu_device=str(parameters["gpu_device"]),
        )
        mark_studio_job({"job": {"run_dir": str(run_dir)}}, owner_id=request["owner_id"], project_id=request["project_id"])
        metadata = read_json(run_dir / "metadata.json")
        metadata.update(active_learning_parent_run=str(campaign), active_learning_round=int(request["round_number"]),
                        active_learning_worker_node=worker["node_id"], active_learning_worker_name=worker["node_name"])
        write_json(run_dir / "metadata.json", metadata)
        result = _wait_local_job(run_dir, cancel_check=cancel_check)
        result.update(node_id=worker["node_id"], node_name=worker["node_name"], host="", worker_run_dir=str(run_dir),
                      cpu_cores=int(worker["cpu_cores"]), gpu_device=str(worker["gpu_device"]))
        return result

    host = str(worker["host"])
    destination = f"{ssh_user}@{host}"
    remote_dir = f"/tmp/mn-active-docking-{os.urandom(12).hex()}"
    remote_request = remote_dir + "/request.json"
    remote_worker = remote_dir + "/worker.py"
    local_dir = Path(request["local_work_dir"])
    local_dir.mkdir(parents=True, exist_ok=True)
    local_request = local_dir / f"{worker['node_id']}.request.json"
    local_request.write_text(json.dumps({
        "schema_version": _DOCKING_WORKER_SCHEMA,
        "parameters": request["parameters"], "round_number": request["round_number"],
        "parent_run_dir": str(campaign), "owner_id": request["owner_id"], "project_id": request["project_id"],
        "node_id": worker["node_id"], "node_name": worker["node_name"],
    }, indent=2) + "\n", encoding="utf-8")
    worker_script = Path(__file__).resolve().parents[3] / "cluster" / "active_learning_worker.py"
    python = str(worker["remote_python"])
    try:
        subprocess.run(["ssh", "-T", "-o", "BatchMode=yes", "-o", "ConnectTimeout=10", destination,
                        "mkdir", "-m", "700", "-p", remote_dir], check=True, capture_output=True, text=True, timeout=30)
        for source, remote_path in ((Path(request["parameters"]["target_pdb"]), remote_dir + "/receptor.pdb"),
                                    (Path(request["parameters"]["compound_set_file"]), remote_dir + "/compounds.json"),
                                    (local_request, remote_request), (worker_script, remote_worker)):
            subprocess.run(["scp", "-q", str(source), f"{destination}:{remote_path}"],
                           check=True, capture_output=True, text=True, timeout=None)
        alignment_path = request["parameters"].get("alignment_ligand_sdf")
        if alignment_path:
            subprocess.run(["scp", "-q", str(Path(alignment_path)), f"{destination}:{remote_dir}/alignment_ligand.sdf"],
                           check=True, capture_output=True, text=True, timeout=None)
        # The shard parameters point to the local staged files. Rewrite those
        # two paths in the remote request after transfer.
        rewrite_request = _remote_request_rewrite_command(
            python, remote_request, remote_dir,
            remote_dir + "/alignment_ligand.sdf" if alignment_path else "",
        )
        subprocess.run(
            ["ssh", "-T", "-o", "BatchMode=yes", destination,
             "bash -lc " + shlex.quote(rewrite_request)],
            check=True, capture_output=True, text=True, timeout=30,
        )
        command = ["ssh", "-T", "-o", "BatchMode=yes", "-o", "ConnectTimeout=10", destination,
                   "bash -lc " + shlex.quote(shlex.join([python, remote_worker, remote_request]))]
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        cancel_sent = False
        while True:
            try:
                stdout, stderr = process.communicate(timeout=2)
                break
            except subprocess.TimeoutExpired:
                if cancel_check() and not cancel_sent:
                    subprocess.run(["ssh", "-T", "-o", "BatchMode=yes", "-o", "ConnectTimeout=5", destination,
                                    "touch", remote_dir + "/cancel"],
                                   check=False, capture_output=True, text=True, timeout=10)
                    cancel_sent = True
        marker = next((line.removeprefix("MDS_ACTIVE_DOCKING_RESULT=")
                       for line in reversed((stdout or "").splitlines())
                       if line.startswith("MDS_ACTIVE_DOCKING_RESULT=")), None)
        if marker is None:
            raise RuntimeError(f"Docking worker on {worker['node_name']} did not return a completion record: {(stderr or stdout)[-3000:]}")
        result = json.loads(marker)
        if process.returncode or result.get("status") != "completed":
            raise RuntimeError(f"Docking worker on {worker['node_name']} failed: {result.get('error') or (stderr or '')[-2000:]}")
        remote_run = str(result["run_dir"])
        local_scores = Path(request["local_work_dir"]) / f"{worker['node_id']}.docking_scores.csv"
        remote_scores = remote_run + "/artifacts/docking_scores.csv"
        subprocess.run(["scp", "-q", f"{destination}:{shlex.quote(remote_scores)}", str(local_scores)],
                       check=True, capture_output=True, text=True, timeout=None)
        result.update(node_id=worker["node_id"], node_name=worker["node_name"], host=host,
                      worker_run_dir=remote_run, score_path=str(local_scores), remote_staging_dir=remote_dir,
                      cpu_cores=int(worker["cpu_cores"]), gpu_device=str(worker["gpu_device"]))
        return result
    finally:
        local_request.unlink(missing_ok=True)
        subprocess.run(["ssh", "-T", "-o", "BatchMode=yes", "-o", "ConnectTimeout=5", destination,
                        "rm", "-rf", "--", remote_dir], check=False, capture_output=True, text=True, timeout=10)


def _read_worker_scores(result: dict, rows_by_engine_id: dict) -> tuple[list[dict], list[dict]]:
    path = Path(result.get("score_path") or (Path(result["run_dir"]) / "artifacts/docking_scores.csv"))
    if not path.is_file():
        raise ValueError(f"Docking child job {result.get('job_code') or path.parent.parent.name} produced no score table.")
    parsed, failures, seen = [], [], set()
    worker_fields = {
        "worker_node_id": result["node_id"], "worker_node_name": result["node_name"],
        "worker_host": result.get("host", ""), "worker_run_dir": result["worker_run_dir"],
        "worker_job_code": result.get("job_code", ""),
        "worker_cpu_cores": result.get("cpu_cores"),
        "worker_gpu_device": result.get("gpu_device", "none"),
    }
    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        delimiter = "\t" if sample.splitlines() and sample.splitlines()[0].count("\t") > sample.splitlines()[0].count(",") else ","
        for raw in csv.DictReader(handle, delimiter=delimiter):
            engine_id = str(raw.get("compound_id") or "")
            source = rows_by_engine_id.get(engine_id)
            if source is None:
                continue
            seen.add(engine_id)
            value = raw.get("best_score_kcal_mol") or raw.get("vina_score") or raw.get("score")
            try:
                score = float(value)
            except (TypeError, ValueError):
                score = math.nan
            if str(raw.get("status") or "").lower() != "completed" or not math.isfinite(score):
                failures.append({**source, **worker_fields,
                                 "error": str(raw.get("error") or "Docking did not return a valid score.")})
                continue
            parsed.append({**source, "score": score,
                           "pose_path": str(raw.get("pose_path") or f"artifacts/poses/{engine_id}_out.pdbqt"),
                           **worker_fields})
    for engine_id, source in rows_by_engine_id.items():
        if engine_id not in seen:
            failures.append({**source, **worker_fields,
                             "error": "Docking worker did not return a score-table row."})
    return parsed, failures


def _score_database(path: Path) -> sqlite3.Connection:
    path.parent.mkdir(parents=True, exist_ok=True)
    db = sqlite3.connect(path, timeout=60)
    db.execute("PRAGMA journal_mode=WAL")
    db.execute("PRAGMA synchronous=FULL")
    db.execute("""CREATE TABLE IF NOT EXISTS compounds (
        library_id TEXT NOT NULL, compound_key INTEGER NOT NULL, compound_id TEXT NOT NULL,
        compound_name TEXT NOT NULL, smiles TEXT NOT NULL, fingerprint BLOB NOT NULL,
        best_score REAL NOT NULL, first_round INTEGER NOT NULL, best_round INTEGER NOT NULL,
        pose_path TEXT, worker_node_id TEXT, worker_node_name TEXT, worker_host TEXT,
        worker_run_dir TEXT, worker_job_code TEXT,
        PRIMARY KEY(library_id, compound_key))""")
    db.execute("""CREATE TABLE IF NOT EXISTS failed_compounds (
        library_id TEXT NOT NULL, compound_key INTEGER NOT NULL, compound_id TEXT NOT NULL,
        compound_name TEXT NOT NULL, smiles TEXT NOT NULL, first_failed_round INTEGER NOT NULL,
        last_failed_round INTEGER NOT NULL, error TEXT NOT NULL, worker_node_id TEXT,
        worker_node_name TEXT, worker_host TEXT, worker_run_dir TEXT, worker_job_code TEXT,
        PRIMARY KEY(library_id, compound_key))""")
    db.execute("CREATE INDEX IF NOT EXISTS compounds_score ON compounds(best_score)")
    db.commit()
    return db


def _persist_round_scores(db: sqlite3.Connection, rows: list[dict], *, round_number: int,
                          score_direction: str) -> None:
    ordering = "lower" if score_direction == "Lower is better" else "higher"
    for row in rows:
        existing = db.execute("SELECT best_score,first_round,best_round FROM compounds WHERE library_id=? AND compound_key=?",
                              (row["library_id"], int(row["compound_key"]))).fetchone()
        if existing is None:
            db.execute("DELETE FROM failed_compounds WHERE library_id=? AND compound_key=?",
                       (row["library_id"], int(row["compound_key"])))
            db.execute("""INSERT INTO compounds VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
                       (row["library_id"], int(row["compound_key"]), str(row["compound_id"]),
                        str(row.get("compound_name") or row["compound_id"]),
                        str(row.get("source_smiles") or row.get("smiles") or ""),
                        sqlite3.Binary(row["fingerprint"]), float(row["score"]), int(round_number), int(round_number),
                        row.get("pose_path"), row.get("worker_node_id"), row.get("worker_node_name"),
                        row.get("worker_host"), row.get("worker_run_dir"), row.get("worker_job_code")))
        else:
            better = float(row["score"]) < existing[0] if ordering == "lower" else float(row["score"]) > existing[0]
            if better:
                db.execute("DELETE FROM failed_compounds WHERE library_id=? AND compound_key=?",
                           (row["library_id"], int(row["compound_key"])))
                db.execute("""UPDATE compounds SET best_score=?,best_round=?,pose_path=?,worker_node_id=?,
                    worker_node_name=?,worker_host=?,worker_run_dir=?,worker_job_code=?
                    WHERE library_id=? AND compound_key=?""",
                    (float(row["score"]), int(round_number), row.get("pose_path"), row.get("worker_node_id"),
                     row.get("worker_node_name"), row.get("worker_host"), row.get("worker_run_dir"), row.get("worker_job_code"),
                     row["library_id"], int(row["compound_key"])))
    db.commit()


def _persist_round_failures(db: sqlite3.Connection, rows: list[dict], *, round_number: int) -> None:
    for row in rows:
        library_id, compound_key = str(row["library_id"]), int(row["compound_key"])
        if db.execute("SELECT 1 FROM compounds WHERE library_id=? AND compound_key=?",
                      (library_id, compound_key)).fetchone():
            continue
        db.execute("""INSERT INTO failed_compounds VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)
            ON CONFLICT(library_id,compound_key) DO UPDATE SET
            last_failed_round=excluded.last_failed_round,error=excluded.error,
            worker_node_id=excluded.worker_node_id,worker_node_name=excluded.worker_node_name,
            worker_host=excluded.worker_host,worker_run_dir=excluded.worker_run_dir,
            worker_job_code=excluded.worker_job_code""",
            (library_id, compound_key, str(row["compound_id"]),
             str(row.get("compound_name") or row["compound_id"]),
             str(row.get("source_smiles") or row.get("smiles") or ""),
             int(round_number), int(round_number), str(row.get("error") or "Docking failed."),
             row.get("worker_node_id"), row.get("worker_node_name"), row.get("worker_host"),
             row.get("worker_run_dir"), row.get("worker_job_code")))
    db.commit()


def _ranked_rows(db: sqlite3.Connection, direction: str) -> list[dict]:
    order = "ASC" if direction == "Lower is better" else "DESC"
    cursor = db.execute(f"SELECT library_id,compound_key,compound_id,compound_name,smiles,best_score,first_round,best_round,pose_path,worker_node_id,worker_node_name,worker_host,worker_run_dir,worker_job_code FROM compounds ORDER BY best_score {order},library_id,compound_key")
    names = [column[0] for column in cursor.description]
    return [dict(zip(names, row)) for row in cursor]


def _top_vector(scores: list[dict], top_n: int, previous: list[dict] | None,
                tolerance: float) -> tuple[list[dict], bool, float | None, bool]:
    current = [{"library_id": row["library_id"], "compound_key": int(row["compound_key"]),
                "compound_id": row["compound_id"], "score": float(row["best_score"])}
               for row in scores[:top_n]]
    if previous is None or len(current) < top_n or len(previous) < top_n or not all(isinstance(item, dict) for item in previous):
        return current, False, None, False
    same_order = all((old.get("library_id"), int(old.get("compound_key", -1))) ==
                     (new["library_id"], new["compound_key"])
                     for old, new in zip(previous, current))
    deltas = [abs(float(old["score"]) - new["score"]) for old, new in zip(previous, current)]
    max_delta = max(deltas, default=0.0)
    return current, same_order and all(value <= tolerance for value in deltas), max_delta, same_order


def _count_library_compounds(db: sqlite3.Connection, library_id: str) -> int:
    return int(db.execute("""SELECT COUNT(*) FROM (
        SELECT compound_key FROM compounds WHERE library_id=?
        UNION SELECT compound_key FROM failed_compounds WHERE library_id=?)""",
        (str(library_id), str(library_id))).fetchone()[0])


def _save_score_table(path: Path, scores: list[dict]) -> None:
    import pyarrow as pa
    import pyarrow.parquet as pq
    schema = pa.schema([
        ("library_id", pa.string()), ("compound_key", pa.uint64()), ("compound_id", pa.string()),
        ("compound_name", pa.string()), ("smiles", pa.string()), ("score", pa.float64()),
        ("score_field", pa.string()), ("score_units", pa.string()), ("rank", pa.uint64()),
        ("first_round", pa.uint32()), ("best_round", pa.uint32()),
        ("pose_path", pa.string()), ("worker_node_id", pa.string()), ("worker_node_name", pa.string()),
        ("worker_host", pa.string()), ("worker_job_code", pa.string()),
    ])
    records = [{**row, "rank": rank} for rank, row in enumerate(scores, start=1)]
    pq.write_table(pa.Table.from_pylist(records, schema=schema), path, compression="zstd")


def _save_failed_compound_table(path: Path, failures: list[dict]) -> None:
    import pyarrow as pa
    import pyarrow.parquet as pq
    schema = pa.schema([
        ("library_id", pa.string()), ("compound_key", pa.uint64()), ("compound_id", pa.string()),
        ("compound_name", pa.string()), ("smiles", pa.string()),
        ("first_failed_round", pa.uint32()), ("last_failed_round", pa.uint32()),
        ("error", pa.string()), ("worker_node_id", pa.string()), ("worker_node_name", pa.string()),
        ("worker_host", pa.string()), ("worker_job_code", pa.string()),
    ])
    pq.write_table(pa.Table.from_pylist(failures, schema=schema), path, compression="zstd")


def _collect_best_poses(campaign: Path, scores: list[dict], count: int, *, ssh_user: str) -> list[str]:
    selected = [row for row in scores[:count] if row.get("pose_path") and row.get("worker_run_dir")]
    destination = campaign / "artifacts" / "active_learning" / "retained_poses"
    destination.mkdir(parents=True, exist_ok=True)
    manifest = []
    for index, row in enumerate(selected, start=1):
        filename = f"{index:05d}_{row['library_id'][-8:]}_{row['compound_key']}.pdbqt"
        local_path = destination / filename
        remote = row.get("worker_host")
        relative = str(row["pose_path"])
        remote_path = str((Path(row["worker_run_dir"]) / relative).resolve())
        if remote:
            from mn_molecular_design_studio.cluster.distributed_library import _HOST
            if not _HOST.fullmatch(str(remote)):
                raise ValueError("A docking worker returned an invalid workstation host name.")
            command = ["scp", "-q", f"{ssh_user}@{remote}:{shlex.quote(remote_path)}", str(local_path)]
            subprocess.run(command, check=True, capture_output=True, text=True, timeout=None)
        else:
            source = Path(row["worker_run_dir"]) / relative
            if not source.is_file():
                continue
            shutil.copy2(source, local_path)
        manifest.append({"path": str(local_path.relative_to(campaign)),
                         "library_id": row["library_id"], "compound_key": row["compound_key"],
                         "compound_id": row["compound_id"], "score": row["best_score"],
                         "worker_node_id": row.get("worker_node_id"), "worker_job_code": row.get("worker_job_code")})
    manifest_path = campaign / "artifacts" / "active_learning" / "retained_poses.json"
    write_json(manifest_path, {"retained_count": len(manifest), "poses": manifest})
    return [str(item["path"]) for item in manifest]


def _cleanup_docking_artifacts(worker_runs: list[dict], *, ssh_user: str) -> None:
    """Keep scores and job records, and prune bulky per-ligand work after pose retention."""
    from mn_molecular_design_studio.cluster.distributed_library import _HOST
    for item in worker_runs:
        host = str(item.get("host") or "")
        run_dirs = item.get("run_dirs") or [item.get("run_dir")]
        for run_dir in run_dirs:
            run_dir = str(run_dir or "")
            if not run_dir:
                continue
            bulky = [f"{run_dir}/artifacts/{name}" for name in ("poses", "ligands", "prepared")]
            if host:
                if not _HOST.fullmatch(host):
                    continue
                cleanup = shlex.join(["rm", "-rf", "--", *bulky])
                subprocess.run(["ssh", "-T", "-o", "BatchMode=yes", "-o", "ConnectTimeout=5",
                                f"{ssh_user}@{host}", "bash -lc " + shlex.quote(cleanup)],
                               check=False, capture_output=True, text=True, timeout=30)
            else:
                for path in bulky:
                    shutil.rmtree(path, ignore_errors=True)


def _queue_model_step(*, campaign: Path, feature_manifest: Path, database: Path,
                      settings: dict, docking_node: dict, round_number: int,
                      binder_percent: float, screening_library_id: str,
                      predicted_count: int, random_count: int, selection_seed: int,
                      progress=None, cancel_check=None) -> tuple[Path, Path]:
    from mn_molecular_design_studio.engine.workflow_engine.core.workflow_queue import enqueue_workflow_call

    workers = docking_node.get("worker_nodes") or []
    local = next((item for item in workers if item.get("local")), None)
    if local and str(local.get("gpu_device") or "none").lower() not in {"none", "cpu", "off"}:
        import torch
        if torch.cuda.is_available():
            gpu_device = str(local["gpu_device"])
            cpu_cores = int(local.get("cpu_cores") or 4)
        else:
            gpu_device = "none"
            try:
                from mn_molecular_design_studio.engine.workflow_engine.core.scheduler import configured_cpu_slots
                cpu_cores = min(16, max(0, configured_cpu_slots() - 1))
            except Exception:
                cpu_cores = 3
    else:
        gpu_device = "none"
        try:
            from mn_molecular_design_studio.engine.workflow_engine.core.scheduler import configured_cpu_slots
            cpu_cores = min(16, max(0, configured_cpu_slots() - 1))
        except Exception:
            cpu_cores = 3
    if cpu_cores < 1:
        raise ValueError("No scheduled CPU slots are available for model training and prediction.")
    kwargs = {
        "campaign_dir": str(campaign), "feature_manifest": str(feature_manifest),
        "scores_database": str(database), "model_path": str(campaign / "artifacts" / "active_learning" / "model.pt"),
        "round_number": int(round_number), "binder_percent": float(binder_percent),
        "hidden_nodes": int(settings["hidden_nodes"]), "dropout": float(settings["dropout"]),
        "epochs": int(settings["epochs"]), "training_batch_size": int(settings["training_batch_size"]),
        "prediction_batch_size": int(settings["prediction_batch_size"]),
        "predicted_count": int(predicted_count), "random_count": int(random_count),
        "screening_library_id": str(screening_library_id), "selection_seed": int(selection_seed),
        "cpu_cores": cpu_cores, "gpu_device": gpu_device,
    }
    run_dir = enqueue_workflow_call(
        run_active_learning_model_step, task_group="campaign", tool="active_learning_model",
        job_type="active_learning_model", kwargs=kwargs, cpu_cores=cpu_cores, gpu_device=gpu_device,
    )
    from mn_molecular_design_studio.workspace import mark_studio_job
    marker = read_json(campaign / "molecular_design_studio.json")
    mark_studio_job({"job": {"run_dir": str(run_dir)}}, owner_id=str(marker.get("studio_owner_id") or ""),
                    project_id=str(marker.get("project_id") or ""))
    if progress:
        progress(run_dir)
    cancellation_sent = False
    while True:
        if cancel_check and cancel_check() and not cancellation_sent:
            from mn_molecular_design_studio.engine.workflow_engine.core.jobs import stop_job
            stop_job(run_dir, reason="Active-learning campaign was cancelled")
            cancellation_sent = True
        metadata = read_json(run_dir / "metadata.json")
        status = str(metadata.get("status") or "queued")
        if status in {"completed", "failed", "cancelled"}:
            if status != "completed" or not read_json(run_dir / "result.json").get("success"):
                raise RuntimeError(f"Active-learning model step {metadata.get('job_code') or run_dir.name} {status}: {metadata.get('worker_error') or read_json(run_dir / 'result.json').get('error') or 'no result'}")
            selected = Path(str(read_json(run_dir / "result.json").get("outputs", {}).get("selected_compounds") or ""))
            if not selected.is_file():
                raise RuntimeError("The model step completed without its next docking batch.")
            return run_dir, selected
        time.sleep(2)


def run_active_learning_pipeline(graph: dict, pipeline_name: str, owner_id: str, project_id: str,
                                 reuse_completed: bool = True, gpu_device: str = "all") -> Path:
    """Run an active-learning graph as a resumable campaign on the organizer."""
    graph = validate_graph(graph, require_inputs=True)["graph"]
    if graph.get("template_id") != "unidock_active_learning":
        raise ValueError("The active-learning controller only accepts the active-learning graph template.")
    from mn_molecular_design_studio.engine.workflow_engine.workflows.ligand_tools import validate_target_pdb
    from mn_molecular_design_studio.workspace import mark_studio_job

    job = create_job("campaign", "visual_pipeline", "visual_pipeline", {},
                     {"pipeline_name": pipeline_name, "graph": graph, "reuse_completed": False,
                      "template_id": graph["template_id"], "active_learning": True})
    campaign = job.run_dir
    mark_studio_job({"job": {"run_dir": str(campaign)}}, owner_id=owner_id, project_id=project_id)
    graph_path = campaign / "graph_snapshot.json"
    if graph_path.exists() and read_json(graph_path) != graph:
        raise ValueError("The submitted active-learning graph snapshot cannot be changed.")
    if not graph_path.exists():
        graph_path.write_text(json.dumps(graph, indent=2) + "\n", encoding="utf-8")
    nodes = _node_by_kind(graph)
    params = {kind: node["data"]["params"] for kind, node in nodes.items()}
    screening_node, seed_node, target_node = (params["active_library"], params["active_seed"], params["active_target"])
    features = _library_sources(screening_node, seed_node)
    docking = params["active_docking"]
    if str(docking.get("engine") or "UniDock Pro") != "UniDock Pro":
        raise ValueError("The active-learning round controller currently uses the selected UniDock Pro engine.")
    target_path = Path(str(target_node.get("target_pdb") or "")).expanduser().resolve(strict=True)
    validate_target_pdb(target_path)
    target = {**target_node, "target_pdb": str(target_path)}
    workers = _selected_worker_specs(docking)
    settings = ActiveLearningSettings(
        binder_percent_start=float(params["active_binder_rule"].get("binder_percent_start", 10.0)),
        binder_percent_target=float(params["active_binder_rule"].get("binder_percent_target", 0.01)),
        tightening_rounds=int(params["active_binder_rule"].get("tightening_rounds", 10)),
        max_iterations=int(params["active_decision"].get("max_iterations", 0)),
    )
    batch_size = int(docking.get("batch_size") or 0)
    top_n = int(params["active_decision"].get("top_n") or 100)
    stable_limit = int(params["active_decision"].get("stable_rounds") or 3)
    tolerance = float(params["active_decision"].get("score_tolerance") or 0.0)
    score_field = str(docking.get("score_field") or "vina_score")
    score_units = str(docking.get("score_units") or "kcal/mol")
    direction = str(docking.get("score_direction") or "Lower is better")
    retain_poses = int(params["active_accumulate"].get("retained_pose_count") or 0)
    selection_node = params["active_select"]
    seed_fraction_percent = float(selection_node.get("next_round_seed_fraction_percent", 50.0))
    predicted_fraction = float(selection_node.get("predicted_fraction_percent", 80.0)) / 100.0
    training_settings = params["active_train"]
    prediction_settings = params["active_predict"]
    rounds_dir = campaign / "artifacts" / "active_learning" / "rounds"
    rounds_dir.mkdir(parents=True, exist_ok=True)
    data_root = campaign / "artifacts" / "active_learning"
    state_path = data_root / "state.json"
    state = read_json(state_path) or {"schema_version": "studio.active-learning-state.v1", "rounds": [],
                                      "docking_runs": [], "stable_rounds": 0, "previous_top_scores": None}
    state.setdefault("docking_runs", [])
    state.setdefault("initial_seed_count", len(features["seed_sample"]))
    initial_seed_count = int(state["initial_seed_count"])
    score_db = _score_database(data_root / "scores.sqlite3")
    graph_state_path = campaign / "graph_execution.json"
    graph_state = read_json(graph_state_path) or {"graph_id": graph["id"], "template_id": graph["template_id"], "nodes": {}}
    cancel = threading.Event()
    old_signal = signal.getsignal(signal.SIGTERM)

    def handle_term(*_):
        cancel.set()

    signal.signal(signal.SIGTERM, handle_term)

    def cancelled():
        return cancel.is_set() or str(read_json(campaign / "metadata.json").get("status") or "") == "cancelled"

    def checkpoint():
        write_json(state_path, state)
        write_json(graph_state_path, graph_state)

    def phase(node_id: str, label: str, status="running", **details):
        _write_graph_state(graph_state_path, graph_state, node=node_id, status=status, **details)
        update_status(campaign, "running", current_phase=label, progress_label=label)

    try:
        update_status(campaign, "running", campaign_name=pipeline_name, current_phase="Validating independent compound sources")
        for kind in ("active_library", "active_seed", "active_target"):
            phase(nodes[kind]["id"], nodes[kind]["data"].get("label") or kind, "completed",
                  source_library_id=(screening_node.get("library_id") if kind == "active_library" else seed_node.get("library_id") if kind == "active_seed" else None))
        if not state["rounds"]:
            sample_size = int(seed_node.get("sample_size") or 10_000)
            sample_seed = int(seed_node.get("sampling_seed") or 0)
            phase(nodes["active_seed"]["id"], "Sampling the centroid source")
            seed_rows = features["seed_sample"]
            state["initial_seed_count"] = initial_seed_count = len(seed_rows)
            selection = rounds_dir / "round-0001-selection.parquet"
            _write_selection_table(selection, seed_rows, str(seed_node["library_id"]), reason="centroid_seed")
            state["pending_round"] = {"number": 1, "selection_path": str(selection), "selection_kind": "centroid_seed"}
            phase(nodes["active_seed"]["id"], f"Selected {len(seed_rows):,} centroid seed compounds",
                  status="completed", selection_count=len(seed_rows))
            checkpoint()

        previous_top = state.get("previous_top_scores")
        while True:
            if cancelled():
                cancel.set()
                update_status(campaign, "cancelled", progress_label="Active-learning campaign cancelled; completed rounds are saved")
                _write_graph_state(graph_state_path, graph_state, node=nodes["active_results"]["id"], status="cancelled")
                checkpoint()
                return campaign
            pending = state.get("pending_round")
            if not pending:
                break
            round_number = int(pending["number"])
            selection_path = Path(pending["selection_path"])
            selection_rows = _read_selection_table(selection_path)
            if not selection_rows:
                state.pop("pending_round", None)
                break
            round_folder = rounds_dir / f"round-{round_number:04d}"
            round_folder.mkdir(parents=True, exist_ok=True)
            worker_states = [{"node_id": w["node_id"], "name": w["node_name"], "cpu_cores": w["cpu_cores"],
                              "gpu_device": w["gpu_device"], "status": "queued"} for w in workers]
            phase(nodes["active_docking"]["id"], f"Docking round {round_number} · {len(selection_rows):,} compounds",
                  round_number=round_number, compound_count=len(selection_rows), worker_nodes=worker_states)

            def worker_progress(updated):
                nonlocal worker_states
                worker_states = updated
                phase(nodes["active_docking"]["id"], f"Docking round {round_number} · {len(selection_rows):,} compounds",
                      round_number=round_number, compound_count=len(selection_rows), worker_nodes=worker_states)

            all_results, failed_results = _dock_round(
                campaign=campaign, round_number=round_number, rows=selection_rows, workers=workers,
                docking=docking, target=target, owner_id=owner_id, project_id=project_id,
                ssh_user=str(docking.get("ssh_user") or "compute"), work_directory=round_folder,
                cancel_check=cancelled, worker_progress=worker_progress,
            )
            accounted_count = len(all_results) + len(failed_results)
            if accounted_count != len(selection_rows):
                raise RuntimeError(f"Round {round_number} accounted for {accounted_count:,} of {len(selection_rows):,} attempted compounds.")
            if not all_results:
                raise RuntimeError(f"Round {round_number} produced no usable docking scores from {len(selection_rows):,} attempts.")
            worker_summary = _summarize_round_workers(all_results, failed_results)
            phase(nodes["active_docking"]["id"], f"Docking round {round_number} complete",
                  status="completed", round_number=round_number, compound_count=len(selection_rows),
                  scored_count=len(all_results), failed_count=len(failed_results),
                  worker_nodes=[{**item, "status": "completed"} for item in worker_summary])
            _persist_round_scores(score_db, all_results, round_number=round_number, score_direction=direction)
            _persist_round_failures(score_db, failed_results, round_number=round_number)
            ranked = _ranked_rows(score_db, direction)
            current_top, stable_now, max_delta, same_order = _top_vector(ranked, top_n, previous_top, tolerance)
            state["docking_runs"].extend(worker_summary)
            if round_number > 1 and stable_now:
                state["stable_rounds"] = int(state.get("stable_rounds") or 0) + 1
            else:
                state["stable_rounds"] = 0
            history_item = {
                "round": round_number, "attempted_count": len(selection_rows),
                "docked_count": len(all_results), "successful_count": len(all_results),
                "failed_count": len(failed_results),
                "unique_tested_count": len(ranked), "best_score": ranked[0]["best_score"] if ranked else None,
                "binder_percent": settings.percentile(round_number), "top_n": top_n,
                "top_n_stable": stable_now, "stable_round_count": state["stable_rounds"],
                "maximum_rank_score_change": max_delta, "same_top_order": same_order,
                "top_scores": current_top, "workers": worker_summary,
                "selection_kind": pending.get("selection_kind"), "completed_at": utc_now(),
            }
            state["rounds"].append(history_item)
            state["previous_top_scores"] = current_top
            state["pending_round"] = None
            previous_top = current_top
            write_json(round_folder / "round.json", history_item)
            checkpoint()
            phase(nodes["active_accumulate"]["id"], f"Accumulated and ranked {len(ranked):,} unique compounds",
                  status="completed", round_number=round_number, unique_compound_count=len(ranked), score_units=score_units)
            phase(nodes["active_decision"]["id"], "Comparing cumulative top-ranked scores",
                  round_number=round_number, top_n=top_n, maximum_score_change=max_delta,
                  same_top_order=same_order, stable_round_count=state["stable_rounds"])

            screening_count = int(features["feature_manifest"].get("row_count") or 0)
            tested_screening_count = _count_library_compounds(score_db, str(screening_node["library_id"]))
            exhausted = tested_screening_count >= screening_count
            max_rounds = settings.iteration_limit_reached(round_number)
            converged = state["stable_rounds"] >= stable_limit
            if converged or exhausted or max_rounds:
                reason = "stable_top_scores" if converged else "screening_library_exhausted" if exhausted else "maximum_iterations"
                state["stop_reason"] = reason
                state["rounds"][-1]["stop_reason"] = reason
                checkpoint()
                phase(nodes["active_decision"]["id"], f"Stopping: {reason.replace('_', ' ')}",
                      status="completed", round_number=round_number, stop_reason=reason,
                      same_top_order=same_order, maximum_score_change=max_delta)
                break

            binder_percent = settings.percentile(round_number)
            phase(nodes["active_binder_rule"]["id"], f"Relabeling accumulated results at {binder_percent:g}%",
                  round_number=round_number, binder_percent=binder_percent)
            phase(nodes["active_train"]["id"], f"Training the OpenVS-style model after round {round_number}",
                  round_number=round_number, status="running")
            next_round = round_number + 1
            selection_counts = _next_round_selection_counts(
                initial_seed_count=initial_seed_count,
                remaining_compounds=screening_count - tested_screening_count,
                max_batch_size=batch_size,
                seed_fraction_percent=seed_fraction_percent,
                predicted_fraction_percent=predicted_fraction * 100.0,
            )
            remaining = selection_counts["total"]
            if remaining <= 0:
                state["stop_reason"] = "screening_library_exhausted"
                checkpoint()
                break
            predicted_count = selection_counts["predicted"]
            random_count = selection_counts["random"]
            model_run, next_selection = _queue_model_step(
                campaign=campaign, feature_manifest=features["feature_manifest_path"],
                database=data_root / "scores.sqlite3", settings={**training_settings, **prediction_settings},
                docking_node=docking, round_number=round_number, binder_percent=binder_percent,
                screening_library_id=str(screening_node["library_id"]), predicted_count=predicted_count,
                random_count=random_count, selection_seed=int(seed_node.get("sampling_seed") or 0) + next_round,
                cancel_check=cancelled,
            )
            model_params = read_json(model_run / "input.json").get("params") or {}
            local_worker = next(worker for worker in workers if worker["local"])
            model_worker = {"node_id": local_worker["node_id"], "name": local_worker["node_name"],
                            "job_code": read_json(model_run / "metadata.json").get("job_code"),
                            "run_dir": str(model_run), "cpu_cores": model_params.get("cpu_cores"),
                            "gpu_device": model_params.get("gpu_device", "none"), "status": "completed"}
            phase(nodes["active_binder_rule"]["id"], f"Binder labels set to {binder_percent:g}%",
                  status="completed", round_number=round_number, binder_percent=binder_percent)
            phase(nodes["active_train"]["id"], f"Trained model after round {round_number}", status="completed",
                  model_job_code=read_json(model_run / "metadata.json").get("job_code"), round_number=round_number,
                  worker_nodes=[model_worker])
            phase(nodes["active_predict"]["id"], "Predicted across the screening compound library", status="completed",
                  round_number=round_number, prediction_count=int(features["feature_manifest"].get("row_count") or 0),
                  worker_nodes=[model_worker])
            selected_rows = _read_selection_table(next_selection)
            if not selected_rows:
                state["stop_reason"] = "screening_library_exhausted"
                checkpoint()
                break
            next_selection_parent = rounds_dir / f"round-{next_round:04d}-selection.parquet"
            shutil.copy2(next_selection, next_selection_parent)
            state["pending_round"] = {"number": next_round, "selection_path": str(next_selection_parent),
                                       "selection_kind": "predicted_plus_random", "model_run_dir": str(model_run)}
            checkpoint()
            phase(nodes["active_select"]["id"], f"Selected {len(selected_rows):,} compounds for round {next_round}",
                  status="completed", round_number=next_round, predicted_count=predicted_count,
                  random_count=random_count, target_batch_count=selection_counts["target_total"],
                  initial_seed_count=initial_seed_count, seed_fraction_percent=seed_fraction_percent,
                  model_job_code=read_json(model_run / "metadata.json").get("job_code"),
                  worker_nodes=[model_worker])

        ranked = _ranked_rows(score_db, direction)
        ssh_user = str(docking.get("ssh_user") or "compute")
        pose_paths = _collect_best_poses(campaign, ranked, retain_poses, ssh_user=ssh_user)
        pose_manifest = read_json(data_root / "retained_poses.json")
        retained = {(str(item["library_id"]), int(item["compound_key"])): str(item["path"])
                    for item in pose_manifest.get("poses", [])}
        for row in ranked:
            row["pose_path"] = retained.get((str(row["library_id"]), int(row["compound_key"])), "")
        table_path = data_root / "docking_scores.parquet"
        _save_score_table(table_path, [{**row, "score": row["best_score"], "score_field": score_field,
                                       "score_units": score_units} for row in ranked])
        failures = [dict(zip([column[1] for column in score_db.execute("PRAGMA table_info(failed_compounds)").fetchall()], row))
                    for row in score_db.execute("SELECT * FROM failed_compounds ORDER BY library_id,compound_key")]
        failed_path = data_root / "failed_compounds.parquet"
        _save_failed_compound_table(failed_path, failures)
        _cleanup_docking_artifacts(state.get("docking_runs") or [], ssh_user=ssh_user)
        for selection in rounds_dir.glob("round-*-selection.parquet"):
            selection.unlink(missing_ok=True)
        for round_folder in rounds_dir.glob("round-[0-9][0-9][0-9][0-9]"):
            for worker_folder in round_folder.glob("worker-*"):
                shutil.rmtree(worker_folder, ignore_errors=True)
        history_path = data_root / "round_history.json"
        write_json(history_path, {"schema_version": "studio.active-learning-round-history.v1",
                                  "stop_reason": state.get("stop_reason", "completed"), "rounds": state["rounds"],
                                  "score_field": score_field, "score_units": score_units, "score_direction": direction,
                                  "screening_library_id": screening_node["library_id"],
                                  "screening_feature_manifest": str(features["feature_manifest_path"]),
                                  "centroid_library_id": seed_node["library_id"],
                                  "centroid_manifest": str(features["centroid_manifest"]),
                                  "retained_poses": pose_paths})
        state["status"] = "completed"
        checkpoint()
        for kind in ("active_binder_rule", "active_train", "active_predict", "active_select"):
            if graph_state.get("nodes", {}).get(nodes[kind]["id"], {}).get("status") != "completed":
                _write_graph_state(graph_state_path, graph_state, node=nodes[kind]["id"], status="skipped", reason=state.get("stop_reason"))
        _write_graph_state(graph_state_path, graph_state, node=nodes["active_results"]["id"], status="completed",
                           unique_compound_count=len(ranked), round_count=len(state["rounds"]), stop_reason=state.get("stop_reason"))
        score_db.close()
        (data_root / "scores.sqlite3-wal").unlink(missing_ok=True)
        (data_root / "scores.sqlite3-shm").unlink(missing_ok=True)
        # The packed fingerprints are needed during model fitting, but the
        # compact final score table and model checkpoint are the durable outputs.
        (data_root / "scores.sqlite3").unlink(missing_ok=True)
        outputs = {"docking_scores": "artifacts/active_learning/docking_scores.parquet",
                   "failed_compounds": "artifacts/active_learning/failed_compounds.parquet",
                   "round_history": "artifacts/active_learning/round_history.json",
                   "retained_poses": "artifacts/active_learning/retained_poses"}
        downstream = {"docking_scores": "artifacts/active_learning/docking_scores.parquet",
                      "failed_compounds": "artifacts/active_learning/failed_compounds.parquet",
                      "round_history": "artifacts/active_learning/round_history.json"}
        model_path = data_root / "model.pt"
        if model_path.is_file():
            outputs["model"] = "artifacts/active_learning/model.pt"
            downstream["model"] = "artifacts/active_learning/model.pt"
        finish_job(campaign, True, {
            "outputs": outputs,
            "metrics": {"round_count": len(state["rounds"]), "unique_docked_compounds": len(ranked),
                        "failed_compounds": len(failures),
                        "best_compound_id": ranked[0]["compound_id"] if ranked else "",
                        "best_score": ranked[0]["best_score"] if ranked else None,
                        "stop_reason": state.get("stop_reason"), "score_field": score_field,
                        "score_units": score_units, "score_direction": direction,
                        "screening_library_id": screening_node["library_id"],
                        "centroid_library_id": seed_node["library_id"]},
            "downstream_artifacts": downstream,
        })
        return campaign
    except BaseException as exc:
        try:
            score_db.close()
        except Exception:
            pass
        if cancelled():
            update_status(campaign, "cancelled", progress_label="Active-learning campaign cancelled; completed rounds are saved")
            raise
        state["status"] = "failed"
        state["error"] = f"{type(exc).__name__}: {exc}"
        checkpoint()
        _write_graph_state(graph_state_path, graph_state, node=nodes["active_results"]["id"], status="failed", error=str(exc))
        raise
    finally:
        signal.signal(signal.SIGTERM, old_signal)


def _dock_round(*, campaign: Path, round_number: int, rows: list[dict], workers: list[dict],
                docking: dict, target: dict, owner_id: str, project_id: str,
                ssh_user: str, work_directory: Path, cancel_check, worker_progress=None) -> tuple[list[dict], list[dict]]:
    import numpy as np
    import pyarrow as pa
    import pyarrow.parquet as pq
    from mn_molecular_design_studio import ligand_compounds
    from mn_molecular_design_studio.engine.workflow_engine.workflows.ligand_tools import MAX_COMPOUNDS_PER_LAUNCH

    if not rows:
        return []
    if not workers:
        raise ValueError("Select at least one computer for docking.")
    launch_batch_size = int(docking.get("worker_batch_size") or 2_000)
    if not 1 <= launch_batch_size <= MAX_COMPOUNDS_PER_LAUNCH:
        raise ValueError(f"Each computer's UniDock Pro batch must be between 1 and {MAX_COMPOUNDS_PER_LAUNCH:,} compounds.")
    random = np.random.default_rng(int(docking.get("sampling_seed") or 42) + round_number)
    assignments = [[] for _ in workers]
    indexes = random.permutation(len(rows))
    for rank, row_index in enumerate(indexes):
        assignments[rank % len(workers)].append(rows[int(row_index)])
    request_info = []
    marker = read_json(campaign / "molecular_design_studio.json")
    schema = pa.schema([("engine_id", pa.string()), ("library_id", pa.string()), ("compound_key", pa.uint64()),
                        ("compound_id", pa.string()), ("compound_name", pa.string()), ("smiles", pa.string()),
                        ("fingerprint", pa.binary(_PACKED_FINGERPRINT_BYTES)), ("selection_reason", pa.string()),
                        ("predicted_probability", pa.float64())])
    for worker_index, (worker, selected) in enumerate(zip(workers, assignments)):
        if not selected:
            continue
        batches = []
        for batch_index, batch_rows in enumerate(_compound_batches(selected, launch_batch_size), start=1):
            node_directory = work_directory / f"worker-{worker_index:02d}" / f"batch-{batch_index:04d}"
            node_directory.mkdir(parents=True, exist_ok=True)
            compounds, mapping = [], {}
            for index, row in enumerate(batch_rows, start=1):
                engine_id = f"L{index:07d}"
                mapping[engine_id] = row
                compounds.append({"compound_id": engine_id, "compound_name": str(row.get("compound_name") or row["compound_id"]),
                                  "smiles": str(row["source_smiles"] if "source_smiles" in row else row["smiles"])})
            compound_file = node_directory / "compounds.json"
            ligand_compounds.write_compound_set(compound_file,
                name=f"Active-learning round {round_number} · {worker['node_name']} · batch {batch_index}",
                compounds=compounds, owner_id=str(marker.get("studio_owner_id") or ""),
                project_id=str(marker.get("project_id") or ""), created_at=utc_now())
            rows_file = node_directory / "selection.parquet"
            table_rows = [{"engine_id": engine_id, "library_id": row["library_id"], "compound_key": int(row["compound_key"]),
                           "compound_id": str(row["compound_id"]), "compound_name": str(row.get("compound_name") or row["compound_id"]),
                           "smiles": str(row.get("source_smiles") or row.get("smiles") or ""),
                           "fingerprint": bytes(row["fingerprint"]), "selection_reason": str(row.get("selection_reason") or "centroid_seed"),
                           "predicted_probability": row.get("predicted_probability")}
                          for engine_id, row in mapping.items()]
            pq.write_table(pa.Table.from_pylist(table_rows, schema=schema), rows_file, compression="zstd")
            parameters = _docking_params(docking, target, compound_file, round_number=round_number,
                                         worker=worker, limit=len(batch_rows))
            request = {"schema_version": _DOCKING_WORKER_SCHEMA, "round_number": round_number,
                       "local_work_dir": str(node_directory), "parameters": parameters,
                       "owner_id": owner_id, "project_id": project_id}
            batches.append((mapping, request))
        request_info.append((worker, batches))

    results, failures = [], []
    errors = []
    status_by_worker = {worker["node_id"]: {"node_id": worker["node_id"], "name": worker["node_name"],
        "cpu_cores": worker["cpu_cores"], "gpu_device": worker["gpu_device"], "status": "queued",
        "batch_count": len(batches), "completed_batches": 0, "scored_count": 0,
        "failed_count": 0}
        for worker, batches in request_info}

    def publish_progress():
        if worker_progress:
            worker_progress([dict(item) for item in status_by_worker.values()])

    def run_worker_batches(worker, batches):
        worker_rows, worker_failures = [], []
        child_runs = []
        for batch_index, (mapping, request) in enumerate(batches, start=1):
            if cancel_check():
                raise RuntimeError("Active-learning campaign was cancelled before this docking batch.")
            status_by_worker[worker["node_id"]].update(status="running", current_batch=batch_index)
            publish_progress()
            job_result = _dispatch_docking_shard(request=request, worker=worker,
                campaign=campaign, ssh_user=ssh_user, cancel_check=cancel_check)
            if job_result.get("status") != "completed":
                raise RuntimeError(f"Docking batch {batch_index}/{len(batches)} finished as {job_result.get('status')}.")
            score_rows, failed_rows = _read_worker_scores(job_result, mapping)
            worker_rows.extend(score_rows)
            worker_failures.extend(failed_rows)
            child_runs.append({"job_code": job_result.get("job_code"),
                               "run_dir": job_result.get("worker_run_dir") or job_result.get("run_dir"),
                               "host": job_result.get("host", "")})
            status_by_worker[worker["node_id"]].update(completed_batches=batch_index,
                scored_count=len(worker_rows), failed_count=len(worker_failures),
                job_code=job_result.get("job_code"))
            publish_progress()
        return worker_rows, worker_failures, child_runs

    with ThreadPoolExecutor(max_workers=max(1, len(request_info))) as pool:
        futures = {pool.submit(run_worker_batches, worker, batches): worker
                   for worker, batches in request_info}
        for future in as_completed(futures):
            worker = futures[future]
            try:
                score_rows, failed_rows, child_runs = future.result()
                status_by_worker[worker["node_id"]].update(
                    status="completed", child_runs=child_runs,
                    scored_count=len(score_rows),
                    failed_count=len(failed_rows),
                )
                results.extend(score_rows)
                failures.extend(failed_rows)
                publish_progress()
            except Exception as exc:
                errors.append(f"{worker['node_name']}: {exc}")
                status_by_worker[worker["node_id"]].update(status="failed", error=str(exc))
                publish_progress()
    if errors:
        raise RuntimeError("The docking round did not finish on every computer: " + "; ".join(errors))
    return results, failures


def _write_selection_table(path: Path, rows: list[dict], library_id: str, *, reason: str) -> None:
    import pyarrow as pa
    import pyarrow.parquet as pq
    schema = pa.schema([("library_id", pa.string()), ("compound_key", pa.uint64()),
                        ("compound_id", pa.string()), ("compound_name", pa.string()),
                        ("source_smiles", pa.string()), ("smiles", pa.string()),
                        ("fingerprint", pa.binary(_PACKED_FINGERPRINT_BYTES)),
                        ("selection_reason", pa.string()), ("predicted_probability", pa.float64())])
    materialized = []
    for row in rows:
        smiles = str(row.get("source_smiles") or row.get("smiles") or "")
        materialized.append({"library_id": library_id, "compound_key": int(row["compound_key"]),
                             "compound_id": str(row["compound_id"]),
                             "compound_name": str(row.get("compound_name") or row["compound_id"]),
                             "source_smiles": smiles, "smiles": smiles,
                             "fingerprint": bytes(row["fingerprint"]),
                             "selection_reason": str(row.get("selection_reason") or reason),
                             "predicted_probability": row.get("predicted_probability")})
    pq.write_table(pa.Table.from_pylist(materialized, schema=schema), path, compression="zstd")


def _read_selection_table(path: Path) -> list[dict]:
    import pyarrow.parquet as pq
    if not path.is_file():
        raise ValueError(f"Active-learning selection table is missing: {path.name}.")
    return pq.read_table(path).to_pylist()


def _summarize_round_workers(rows: list[dict], failures: list[dict] | None = None) -> list[dict]:
    grouped = {}
    for row, failed in [(row, False) for row in rows] + [(row, True) for row in (failures or [])]:
        item = grouped.setdefault(row["worker_node_id"], {"node_id": row["worker_node_id"],
            "name": row["worker_node_name"], "job_code": row.get("worker_job_code"),
            "host": row.get("worker_host"), "run_dir": row.get("worker_run_dir"),
            "run_dirs": [], "job_codes": [], "cpu_cores": row.get("worker_cpu_cores"),
            "gpu_device": row.get("worker_gpu_device"), "scored_count": 0, "failed_count": 0})
        run_dir, job_code = row.get("worker_run_dir"), row.get("worker_job_code")
        if run_dir and run_dir not in item["run_dirs"]:
            item["run_dirs"].append(run_dir)
        if job_code and job_code not in item["job_codes"]:
            item["job_codes"].append(job_code)
        item["failed_count" if failed else "scored_count"] += 1
    return list(grouped.values())


class _OpenVSDenseModel:
    def __new__(cls, *, hidden_nodes: int, dropout: float):
        import torch.nn as nn
        import torch.nn.functional as F

        class Network(nn.Module):
            def __init__(self):
                super().__init__()
                self.fc1 = nn.Linear(1024, hidden_nodes)
                self.fc2 = nn.Linear(hidden_nodes, hidden_nodes)
                self.fc3 = nn.Linear(hidden_nodes, hidden_nodes)
                self.out = nn.Linear(hidden_nodes, 1)
                self.bn1 = nn.BatchNorm1d(hidden_nodes)
                self.bn2 = nn.BatchNorm1d(hidden_nodes)
                self.bn3 = nn.BatchNorm1d(hidden_nodes)
                self.drop = nn.Dropout(float(dropout))

            def forward(self, x):
                x = self.drop(F.relu(self.bn1(self.fc1(x))))
                x = self.drop(F.relu(self.bn2(self.fc2(x))))
                x = F.relu(self.bn3(self.fc3(x)))
                return self.out(x)

        return Network()


def run_active_learning_model_step(campaign_dir: str, feature_manifest: str,
                                   scores_database: str, model_path: str,
                                   round_number: int, binder_percent: float,
                                   hidden_nodes: int = 3000, dropout: float = 0.5,
                                   epochs: int = 20, training_batch_size: int = 5000,
                                   prediction_batch_size: int = 10000,
                                   predicted_count: int = 8000, random_count: int = 2000,
                                   screening_library_id: str = "", selection_seed: int = 42,
                                   cpu_cores: int = 4, gpu_device: str = "none") -> Path:
    """Retrain on cumulative scores and stream a full-library prediction pass."""
    import numpy as np
    import torch
    from mn_molecular_design_studio.external_libraries import iter_feature_batches, read_feature_manifest

    job = create_job("campaign", "active_learning_model", "active_learning_model",
                     {"feature_manifest": feature_manifest, "scores_database": scores_database},
                     {"round_number": round_number, "binder_percent": binder_percent, "hidden_nodes": hidden_nodes,
                      "dropout": dropout, "epochs": epochs, "training_batch_size": training_batch_size,
                      "prediction_batch_size": prediction_batch_size, "cpu_cores": cpu_cores, "gpu_device": gpu_device})
    work_dir = job.run_dir
    try:
        update_status(work_dir, "running", current_phase="Preparing cumulative OpenVS training data")
        score_db = sqlite3.connect(f"file:{Path(scores_database).resolve()}?mode=ro", uri=True)
        cursor = score_db.execute("SELECT fingerprint,best_score FROM compounds ORDER BY library_id,compound_key")
        packed, values = [], []
        for fingerprint, score in cursor:
            packed.append(bytes(fingerprint))
            values.append(float(score))
        score_db.close()
        if len(packed) < 4:
            raise ValueError("OpenVS training needs at least four uniquely docked compounds.")
        score_values = np.asarray(values, dtype=np.float64)
        percentile = float(binder_percent)
        if not math.isfinite(percentile) or not 0 < percentile <= 100:
            raise ValueError("Binder percentage must be between 0 and 100.")
        cutoff = float(np.percentile(score_values, percentile))
        labels = (score_values <= cutoff).astype(np.float32)
        positive, negative = np.flatnonzero(labels == 1), np.flatnonzero(labels == 0)
        if len(positive) < 2 or len(negative) < 2:
            raise ValueError(f"Cumulative binder labels need at least two compounds in each class (positive={len(positive)}, negative={len(negative)}).")
        packed_matrix = np.frombuffer(b"".join(packed), dtype=np.uint8).reshape(len(packed), _PACKED_FINGERPRINT_BYTES)
        features = np.unpackbits(packed_matrix, axis=1, bitorder="little").astype(np.float32, copy=False)
        generator = np.random.default_rng(42 + int(round_number))
        torch.manual_seed(42 + int(round_number))
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(42 + int(round_number))
        positive, negative = generator.permutation(positive), generator.permutation(negative)
        positive_validation = min(len(positive) - 1, max(1, int(round(len(positive) * 0.2))))
        negative_validation = min(len(negative) - 1, max(1, int(round(len(negative) * 0.2))))
        validation_ids = np.concatenate((positive[:positive_validation], negative[:negative_validation]))
        training_ids = np.concatenate((positive[positive_validation:], negative[negative_validation:]))
        generator.shuffle(training_ids)
        device = torch.device("cuda" if str(gpu_device).lower() not in {"none", "cpu", "off"} and torch.cuda.is_available() else "cpu")
        if device.type == "cpu":
            torch.set_num_threads(max(1, int(cpu_cores)))
        model = _OpenVSDenseModel(hidden_nodes=int(hidden_nodes), dropout=float(dropout)).to(device)
        previous_checkpoint = Path(model_path)
        if previous_checkpoint.is_file():
            try:
                saved = torch.load(previous_checkpoint, map_location="cpu", weights_only=False)
                if saved.get("architecture", {}).get("hidden_nodes") == int(hidden_nodes):
                    model.load_state_dict(saved["model_state"])
            except (OSError, RuntimeError, ValueError, KeyError, TypeError):
                pass
        optimizer = torch.optim.Adam(model.parameters(), lr=1e-4)
        pos_weight = torch.tensor([len(negative) / len(positive)], dtype=torch.float32, device=device)
        loss_fn = torch.nn.BCEWithLogitsLoss(pos_weight=pos_weight)
        x_tensor = torch.from_numpy(features)
        y_tensor = torch.from_numpy(labels.reshape(-1, 1))
        train_batch = max(2, min(int(training_batch_size), len(training_ids)))
        best_state = None
        best_loss = math.inf
        patience = 5
        stale = 0
        history = []
        update_status(work_dir, "running", current_phase=f"Training round {round_number} · {len(labels):,} cumulative docked compounds")
        for epoch in range(max(1, int(epochs))):
            model.train()
            generator.shuffle(training_ids)
            total_loss = 0.0
            steps = 0
            for start in range(0, len(training_ids), train_batch):
                indexes = training_ids[start:start + train_batch]
                if len(indexes) < 2:
                    continue
                xb = x_tensor[indexes].to(device)
                yb = y_tensor[indexes].to(device)
                optimizer.zero_grad(set_to_none=True)
                loss = loss_fn(model(xb), yb)
                loss.backward()
                optimizer.step()
                total_loss += float(loss.detach().cpu())
                steps += 1
            model.eval()
            val_loss_sum = 0.0
            val_steps = 0
            with torch.no_grad():
                for start in range(0, len(validation_ids), train_batch):
                    indexes = validation_ids[start:start + train_batch]
                    loss = torch.nn.functional.binary_cross_entropy_with_logits(
                        model(x_tensor[indexes].to(device)), y_tensor[indexes].to(device))
                    val_loss_sum += float(loss.detach().cpu())
                    val_steps += 1
            val_loss = val_loss_sum / max(1, val_steps)
            history.append({"epoch": epoch + 1, "train_loss": total_loss / max(1, steps), "validation_loss": val_loss})
            if val_loss < best_loss - 1e-4:
                best_loss, stale = val_loss, 0
                best_state = {key: value.detach().cpu().clone() for key, value in model.state_dict().items()}
            else:
                stale += 1
                if stale >= patience:
                    break
        if best_state is not None:
            model.load_state_dict(best_state)
        model.eval()
        final_model = Path(model_path)
        final_model.parent.mkdir(parents=True, exist_ok=True)
        temporary_model = final_model.with_suffix(".pt.tmp")
        torch.save({"architecture": {"name": "OpenVS VanillaNet", "input_bits": 1024,
                                      "hidden_nodes": int(hidden_nodes), "hidden_layers": 3,
                                      "dropout": float(dropout)},
                    "round_number": int(round_number), "binder_percent": percentile, "cutoff": cutoff,
                    "score_count": len(labels), "model_state": {key: value.detach().cpu() for key, value in model.state_dict().items()},
                    "training_history": history}, temporary_model)
        temporary_model.replace(final_model)

        feature_path, feature_manifest_payload = read_feature_manifest(feature_manifest)
        if feature_manifest_payload.get("library_id") != screening_library_id:
            raise ValueError("The prediction fingerprint manifest belongs to another screening library.")
        already = sqlite3.connect(f"file:{Path(scores_database).resolve()}?mode=ro", uri=True)
        docked = {int(row[0]) for row in already.execute(
            "SELECT compound_key FROM compounds WHERE library_id=? UNION "
            "SELECT compound_key FROM failed_compounds WHERE library_id=?",
            (screening_library_id, screening_library_id))}
        already.close()
        top_count = max(0, int(predicted_count))
        top_heap: list[tuple[float, int, dict]] = []
        seen = 0
        update_status(work_dir, "running", current_phase=f"Predicting across {int(feature_manifest_payload.get('row_count') or 0):,} screening compounds")
        for batch in iter_feature_batches(feature_path, batch_size=max(1, int(prediction_batch_size))):
            matrix = np.frombuffer(b"".join(bytes(value) for value in batch.column("fingerprint").to_pylist()), dtype=np.uint8)
            if not matrix.size:
                continue
            matrix = matrix.reshape(batch.num_rows, _PACKED_FINGERPRINT_BYTES)
            bit_values = np.unpackbits(matrix, axis=1, bitorder="little").astype(np.float32, copy=False)
            with torch.no_grad():
                probabilities = torch.sigmoid(model(torch.from_numpy(bit_values).to(device))).reshape(-1).cpu().numpy()
            keys = batch.column("compound_key").to_pylist()
            ids = batch.column("compound_id").to_pylist()
            names = batch.column("compound_name").to_pylist()
            smiles = batch.column("source_smiles").to_pylist()
            fps = batch.column("fingerprint").to_pylist()
            for offset, probability in enumerate(probabilities):
                key = int(keys[offset])
                if key in docked or top_count == 0:
                    continue
                row = {"library_id": screening_library_id, "compound_key": key, "compound_id": ids[offset],
                       "compound_name": names[offset] or ids[offset], "source_smiles": smiles[offset],
                       "fingerprint": bytes(fps[offset]), "selection_reason": "predicted_hit",
                       "predicted_probability": float(probability)}
                entry = (float(probability), -seen, row)
                if len(top_heap) < top_count:
                    heapq.heappush(top_heap, entry)
                elif entry[0] > top_heap[0][0]:
                    heapq.heapreplace(top_heap, entry)
                seen += 1
        top_rows = [item[2] for item in sorted(top_heap, key=lambda item: (-item[0], -item[1]))]
        top_ids = {int(row["compound_key"]) for row in top_rows}
        rng = np.random.default_rng(int(selection_seed))
        sample: list[dict] = []
        available = 0
        if random_count > 0:
            update_status(work_dir, "running", current_phase="Sampling random exploration compounds")
            for batch in iter_feature_batches(feature_path, batch_size=max(1, int(prediction_batch_size))):
                keys = batch.column("compound_key").to_pylist()
                ids = batch.column("compound_id").to_pylist()
                names = batch.column("compound_name").to_pylist()
                smiles = batch.column("source_smiles").to_pylist()
                fps = batch.column("fingerprint").to_pylist()
                for offset, raw_key in enumerate(keys):
                    key = int(raw_key)
                    if key in docked or key in top_ids:
                        continue
                    available += 1
                    row = {"library_id": screening_library_id, "compound_key": key, "compound_id": ids[offset],
                           "compound_name": names[offset] or ids[offset], "source_smiles": smiles[offset],
                           "fingerprint": bytes(fps[offset]), "selection_reason": "random_exploration",
                           "predicted_probability": None}
                    if len(sample) < random_count:
                        sample.append(row)
                    else:
                        choice = int(rng.integers(0, available))
                        if choice < random_count:
                            sample[choice] = row
        combined = top_rows + sample
        selection_path = work_dir / "artifacts" / "selected_compounds.parquet"
        selection_path.parent.mkdir(parents=True, exist_ok=True)
        _write_selection_table(selection_path, combined, screening_library_id, reason="model_selected")
        if device.type == "cuda":
            torch.cuda.empty_cache()
        finish_job(work_dir, True, {
            "outputs": {"selected_compounds": str(selection_path), "model": str(final_model)},
            "metrics": {"round_number": int(round_number), "binder_percent": percentile, "binder_score_cutoff": cutoff,
                        "training_compounds": len(labels), "positive_count": int(labels.sum()),
                        "negative_count": int(len(labels) - labels.sum()), "epochs_completed": len(history),
                        "best_validation_loss": best_loss, "training_device": str(device),
                        "predicted_selection_count": len(top_rows),
                        "random_selection_count": len(sample), "undocked_candidates_scanned": seen},
            "downstream_artifacts": {"selected_compounds": str(selection_path), "model": str(final_model)},
        })
        return work_dir
    except Exception as exc:
        finish_job(work_dir, False, {"error": str(exc), "outputs": {}})
        raise
