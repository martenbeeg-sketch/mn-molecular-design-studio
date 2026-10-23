"""Internal scheduled CPU worker for bounded centroid chunks."""
from __future__ import annotations

import argparse
from collections import deque
from concurrent.futures import ProcessPoolExecutor
import importlib.util
import json
from multiprocessing import get_context
import os
from pathlib import Path
import signal
import sys
import threading


def _backend():
    path = Path(os.environ["MDS_EXTERNAL_LIBRARIES_MODULE"]).with_name("centroid_cpu.py")
    spec = importlib.util.spec_from_file_location("studio_centroid_backend", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _task(input_path, output_path, threshold, cores):
    return _backend().cluster_chunk_file(input_path, output_path, threshold, cores)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--workers", required=True, type=int)
    parser.add_argument("--run-id", default="")
    parser.add_argument("--tasks", required=True, type=int)
    parser.add_argument("--threshold", required=True, type=float)
    args = parser.parse_args()
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=False)
    lease = None
    stop = threading.Event()
    heart = None
    pool = None
    pidfile = output.parent / "worker.pid"
    pidfile.write_text(str(os.getpid()))

    def terminate(signum, frame):
        if pool:
            for process in list(pool._processes.values()):
                process.terminate()
        raise SystemExit("Centroid workers cancelled.")

    signal.signal(signal.SIGTERM, terminate)
    try:
        if os.environ.get("MDS_LIBRARY_SCHEDULED_WORKER") == "1":
            from mn_molecular_design_studio.engine.mn_compute_scheduler.resources import (
                acquire_cpu_lease, cpu_pool_capacity, shared_state_dir,
            )
            lease = acquire_cpu_lease(
                args.workers, run_id=args.run_id, worker_id=f"centroid-{os.getpid()}",
                workflow="library-centroid-worker", app_id="mn-molecular-design-studio",
                state_root=shared_state_dir(), capacity=cpu_pool_capacity(), stale_after_seconds=120,
            )
            if lease is None:
                raise RuntimeError("The workstation scheduler has no free CPU slots for centroid workers.")

            def heartbeat():
                while not stop.wait(30):
                    lease.heartbeat()

            heart = threading.Thread(target=heartbeat, daemon=True)
            heart.start()
        _backend()._native_library()  # Build once before child processes start.
        print("MDS_LIBRARY_READY", flush=True)
        parallel = max(1, min(args.workers, args.tasks))
        cores_per_chunk = max(1, args.workers // parallel)
        pool = ProcessPoolExecutor(max_workers=parallel, mp_context=get_context("spawn"))
        pending = deque()
        results = []
        for line in sys.stdin:
            command = json.loads(line)
            index = int(command["chunk"])
            pending.append(pool.submit(_task, command["input"], str(output / f"chunk-{index:08d}.parquet"),
                                       args.threshold, cores_per_chunk))
            if len(pending) >= parallel * 2:
                results.append(pending.popleft().result())
        while pending:
            results.append(pending.popleft().result())
        pool.shutdown(wait=True)
        pool = None
        (output / "worker-manifest.json").write_text(json.dumps({
            "chunks": results, "cpu_worker_count": args.workers,
            "parallel_chunks": parallel, "cores_per_chunk": cores_per_chunk,
        }, indent=2))
    except BaseException:
        if pool and pool._processes:
            for process in list(pool._processes.values()):
                process.terminate()
        raise
    finally:
        if pool:
            pool.shutdown(wait=True, cancel_futures=True)
        stop.set()
        if heart:
            heart.join(timeout=2)
        if lease:
            lease.release()
        pidfile.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
