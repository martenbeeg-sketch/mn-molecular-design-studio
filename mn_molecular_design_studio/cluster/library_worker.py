"""Internal streaming worker used by one Studio library-preparation job."""
from __future__ import annotations

import argparse
import importlib
import importlib.util
import json
import os
import sys
import threading
from collections import deque
from pathlib import Path


def _load_library_module():
    """Load the organizer's exact fingerprint implementation on remote workers."""
    module_name = "mn_molecular_design_studio.external_libraries"
    module_path = Path(os.environ["MDS_EXTERNAL_LIBRARIES_MODULE"]).resolve(strict=True)
    package = importlib.import_module("mn_molecular_design_studio")
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load the Studio fingerprint implementation.")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    setattr(package, "external_libraries", module)
    spec.loader.exec_module(module)
    return module


# ProcessPoolExecutor uses spawn. Loading the same module at import time makes
# its worker processes able to resolve the fingerprint function by module name.
library = _load_library_module()


def run(input_stream, output: Path, workers: int, batch_size: int, shard_rows: int):
    import pyarrow as pa
    import pyarrow.parquet as pq
    from concurrent.futures import ProcessPoolExecutor
    from multiprocessing import get_context

    output.mkdir(parents=True, exist_ok=False)
    spec = library.FingerprintSpec()
    schema = library.feature_arrow_schema(spec)
    rejection_schema = pa.schema([
        ("row_number", pa.uint64()), ("compound_id", pa.string()),
        ("smiles", pa.string()), ("reason", pa.string()),
    ])
    reject_path = output / "rejected.parquet"
    reject_writer = pq.ParquetWriter(reject_path, rejection_schema, compression="zstd")
    reject_rows = []
    result_rows = []
    feature_writer = None
    feature_count = total_features = total_rejected = input_rows = 0
    shard_index = 0
    shards = []
    pending = deque()
    pool = ProcessPoolExecutor(max_workers=workers, mp_context=get_context("spawn")) if workers > 1 else None

    def close_shard():
        nonlocal feature_writer, feature_count, shard_index
        if feature_writer is None:
            return
        feature_writer.close()
        path = output / f"features-{shard_index:06d}.parquet"
        shards.append({"path": path.name, "rows": feature_count,
                       "sha256": library.digest_file(path), "size_bytes": path.stat().st_size})
        shard_index += 1
        feature_count = 0
        feature_writer = None

    def write_results(items):
        nonlocal total_features, total_rejected, feature_count, feature_writer
        for item in items:
            if item.get("reason"):
                total_rejected += 1
                reject_rows.append({"row_number": item["row_number"], "compound_id": item.get("compound_id"),
                                    "smiles": item.get("source_smiles"), "reason": item["reason"]})
                if len(reject_rows) >= batch_size:
                    reject_writer.write_table(pa.Table.from_pylist(reject_rows, schema=rejection_schema))
                    reject_rows.clear()
                continue
            result_rows.append({"compound_key": int(item["row_number"]) - 1,
                                "compound_id": item["compound_id"], "compound_name": item["compound_name"],
                                "source_smiles": item["source_smiles"], "fingerprint": item["fingerprint"]})
            total_features += 1
            if len(result_rows) >= batch_size:
                if feature_writer is None:
                    feature_writer = pq.ParquetWriter(output / f"features-{shard_index:06d}.parquet", schema,
                                                      compression="zstd")
                feature_writer.write_table(pa.Table.from_pylist(result_rows, schema=schema))
                feature_count += len(result_rows)
                result_rows.clear()
                if feature_count >= shard_rows:
                    close_shard()

    def dispatch(rows, *, drain=False):
        if not rows:
            return
        records = tuple((item["row_number"], item.get("compound_id"), item.get("compound_name"),
                         item.get("source_smiles") or "") for item in rows)
        if pool:
            pending.append(pool.submit(library._fingerprint_batch, records, spec))
            if drain or len(pending) >= workers * 2:
                write_results(pending.popleft().result())
        else:
            write_results(library._fingerprint_batch(records, spec))

    try:
        batch = []
        for line in input_stream:
            item = json.loads(line)
            input_rows += 1
            if item.get("source_error"):
                total_rejected += 1
                reject_rows.append({"row_number": item["row_number"], "compound_id": item.get("compound_id"),
                                    "smiles": item.get("source_smiles"), "reason": item["source_error"]})
                if len(reject_rows) >= batch_size:
                    reject_writer.write_table(pa.Table.from_pylist(reject_rows, schema=rejection_schema))
                    reject_rows.clear()
                continue
            batch.append(item)
            if len(batch) >= batch_size:
                dispatch(batch)
                batch = []
        if batch:
            dispatch(batch, drain=True)
        while pending:
            write_results(pending.popleft().result())
        if result_rows:
            if feature_writer is None:
                feature_writer = pq.ParquetWriter(output / f"features-{shard_index:06d}.parquet", schema,
                                                  compression="zstd")
            feature_writer.write_table(pa.Table.from_pylist(result_rows, schema=schema))
            feature_count += len(result_rows)
            result_rows.clear()
        close_shard()
        if reject_rows:
            reject_writer.write_table(pa.Table.from_pylist(reject_rows, schema=rejection_schema))
            reject_rows.clear()
        reject_writer.close()
        manifest = {"status": "completed", "input_rows": input_rows, "row_count": total_features,
                    "rejected_rows": total_rejected, "shards": shards,
                    "rejections": {"path": "rejected.parquet", "sha256": library.digest_file(reject_path)}}
        library.atomic_json(output / "worker.json", manifest)
    except BaseException:
        reject_writer.close()
        if feature_writer is not None:
            feature_writer.close()
        if pool:
            pool.shutdown(wait=True, cancel_futures=True)
        raise
    if pool:
        pool.shutdown(wait=True, cancel_futures=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--workers", type=int, required=True)
    parser.add_argument("--batch-size", type=int, default=2048)
    parser.add_argument("--shard-rows", type=int, default=1_000_000)
    parser.add_argument("--run-id", default="")
    args = parser.parse_args()
    cpu_lease = None
    heartbeat_stop = threading.Event()
    heartbeat_thread = None
    if os.environ.get("MDS_LIBRARY_SCHEDULED_WORKER") == "1":
        from mn_molecular_design_studio.engine.mn_compute_scheduler.resources import (
            acquire_cpu_lease, cpu_pool_capacity, shared_state_dir,
        )
        cpu_lease = acquire_cpu_lease(
            args.workers, run_id=args.run_id or "library-worker", worker_id=f"library-{os.getpid()}",
            workflow="library-fingerprint-worker", app_id="mn-molecular-design-studio",
            state_root=shared_state_dir(), capacity=cpu_pool_capacity(), stale_after_seconds=120,
        )
        if cpu_lease is None:
            raise RuntimeError("The workstation scheduler has no free CPU allocation for this worker group.")

        def heartbeat():
            while not heartbeat_stop.wait(30):
                cpu_lease.heartbeat()

        heartbeat_thread = threading.Thread(target=heartbeat, name="library-cpu-lease-heartbeat", daemon=True)
        heartbeat_thread.start()
    try:
        print("MDS_LIBRARY_READY", flush=True)
        run(sys.stdin, Path(args.output), args.workers, args.batch_size, args.shard_rows)
    finally:
        heartbeat_stop.set()
        if heartbeat_thread:
            heartbeat_thread.join(timeout=2)
        if cpu_lease:
            cpu_lease.release()


if __name__ == "__main__":
    main()
