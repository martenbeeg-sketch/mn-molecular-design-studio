from __future__ import annotations

import argparse
import json
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any


def _run_one(job: dict[str, Any], work_dir: Path) -> dict[str, Any]:
    job_id = str(job.get("id") or "candidate")
    command = job.get("command")
    if not isinstance(command, list) or not command:
        return {"id": job_id, "return_code": 2, "error": "missing command"}

    log_path = Path(str(job.get("log_path") or f"logs/{job_id}.log"))
    if not log_path.is_absolute():
        log_path = work_dir / log_path
    log_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        completed = subprocess.run(
            [str(part) for part in command],
            cwd=work_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
        log_path.write_text(
            f"$ {' '.join(str(part) for part in command)}\n{completed.stdout or ''}",
            encoding="utf-8",
        )
        return {
            "id": job_id,
            "return_code": int(completed.returncode),
            "log_path": str(log_path),
        }
    except OSError as exc:
        log_path.write_text(f"Could not run candidate scorer: {exc}\n", encoding="utf-8")
        return {"id": job_id, "return_code": 127, "log_path": str(log_path), "error": str(exc)}


def main() -> int:
    parser = argparse.ArgumentParser(description="Run a batch of ipSAE candidate scoring commands.")
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()

    manifest_path = Path(args.manifest)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    work_dir = Path(str(manifest.get("work_dir") or "/work"))
    jobs = manifest.get("jobs")
    if not isinstance(jobs, list):
        raise ValueError("Batch manifest must contain a jobs list.")
    workers = max(1, min(int(args.workers), len(jobs))) if jobs else 1
    results: list[dict[str, Any]] = []

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(_run_one, job, work_dir) for job in jobs]
        for future in as_completed(futures):
            results.append(future.result())

    order = {str(job.get("id") or "candidate"): index for index, job in enumerate(jobs)}
    results.sort(key=lambda result: order.get(str(result.get("id")), len(order)))
    status_path = Path(str(manifest.get("status_path") or "batch_status.json"))
    if not status_path.is_absolute():
        status_path = work_dir / status_path
    status_path.parent.mkdir(parents=True, exist_ok=True)
    status_path.write_text(
        json.dumps({"workers": workers, "jobs": results}, indent=2) + "\n",
        encoding="utf-8",
    )
    failures = sum(1 for result in results if int(result.get("return_code") or 0) != 0)
    print(f"ipSAE batch finished: {len(results)} candidate(s), {failures} command failure(s), {workers} worker(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
