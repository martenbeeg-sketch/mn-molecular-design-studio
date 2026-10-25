"""Run one streamed library fingerprint job across selected Studio computers."""
from __future__ import annotations

import json
import os
import queue
import random
import re
import shlex
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import threading
import uuid
from pathlib import Path

from . import config


_HOST = re.compile(r"[A-Za-z0-9_.-]{1,120}\Z")
_SSH_USER = re.compile(r"[A-Za-z_][A-Za-z0-9_-]{0,31}\Z")


def _live_local_resources():
    """Read local scheduler and storage state without relying on peer heartbeats."""
    from .agent import Agent
    from mn_molecular_design_studio.cluster.capabilities import resources

    scheduler_online, snapshot = resources()
    settings = config.load()
    accepting = bool(settings.get("accept_jobs") and Agent().storage_online())
    return bool(scheduler_online), accepting, snapshot


def _live_remote_resources(host, ssh_user):
    """Verify a remote CPU worker through its own Studio and scheduler runtimes."""
    spec = {"host": host, "name": host}
    python, _repository, _scheduler_state, _scheduler_cpu_slots = _remote_environment(spec, ssh_user)
    probe = (
        "import json; from mn_molecular_design_studio.cluster import config; "
        "from mn_molecular_design_studio.cluster.agent import Agent; "
        "from mn_molecular_design_studio.cluster.capabilities import resources; "
        "online, snapshot = resources(); settings = config.load(); "
        "accepting = bool(settings.get('accept_jobs') and Agent().storage_online()); "
        "print('MDS_RESOURCE_JSON=' + json.dumps({'scheduler_online':bool(online),"
        "'accept_jobs':accepting,'resources':snapshot},separators=(',',':')))"
    )
    remote_command = "bash -lc " + shlex.quote(shlex.join([python, "-c", probe]))
    result = subprocess.run(
        ["ssh", "-T", "-o", "BatchMode=yes", "-o", "ConnectTimeout=10",
         f"{ssh_user}@{host}", remote_command],
        capture_output=True, text=True, timeout=20, check=False,
    )
    if result.returncode:
        detail = (result.stderr or result.stdout).strip()
        raise RuntimeError(f"Could not check the live scheduler on {host}: {detail[-800:] or 'SSH probe failed'}")
    marker = next((line[len("MDS_RESOURCE_JSON="):]
                   for line in reversed(result.stdout.splitlines())
                   if line.startswith("MDS_RESOURCE_JSON=")), None)
    if marker is None:
        raise RuntimeError(f"The live scheduler check on {host} returned no resource snapshot.")
    try:
        payload = json.loads(marker)
    except json.JSONDecodeError as exc:
        raise RuntimeError(f"The live scheduler check on {host} returned invalid resource data.") from exc
    return bool(payload.get("scheduler_online")), bool(payload.get("accept_jobs")), payload.get("resources") or {}


def _node_specs(worker_nodes, *, ssh_user="compute"):
    from .store import Store

    nodes = Store().nodes()
    by_id = {row["id"]: row for row in nodes}
    specs = []
    seen = set()
    for raw in worker_nodes:
        node_id = str(raw.get("node_id") or "")
        if node_id in seen:
            raise ValueError("Choose each worker computer only once.")
        seen.add(node_id)
        node = by_id.get(node_id)
        if not node or node.get("disabled"):
            raise ValueError("Every selected worker computer must be registered and enabled.")
        count = raw.get("cpu_cores")
        chunk_weight = raw.get("chunk_weight", count)
        max_parallel_chunks = raw.get("max_parallel_chunks", count)
        local = node_id == config.identity()
        if local:
            scheduler_online, accepting, resources = _live_local_resources()
        else:
            host = config.display_name(node["name"])
            if not _HOST.fullmatch(host):
                raise ValueError(f"Workstation {node['name']} does not have a safe SSH host name.")
            scheduler_online, accepting, resources = _live_remote_resources(host, ssh_user)
        if not scheduler_online or not accepting:
            raise ValueError(f"The live scheduler on {node['name']} is unavailable or not accepting jobs.")
        available = int(resources.get("cpu_free") or 0)
        limit = int(resources.get("cpu_total") or available or 0)
        # The parent Studio job already owns the organizer's CPU allocation
        # while it runs. Remote worker groups obtain their own scheduler lease.
        allowed = min(256, limit, limit if local else available)
        if type(count) is not int or count < 1 or count > allowed:
            raise ValueError(
                f"Worker count for {node['name']} must fit its live CPU allocation "
                f"({allowed} slot{'s' if allowed != 1 else ''} available)."
            )
        if type(chunk_weight) is not int or chunk_weight < 1 or chunk_weight > 1000:
            raise ValueError(f"Chunk weight for {node['name']} must be a whole number from 1 to 1,000.")
        if type(max_parallel_chunks) is not int or max_parallel_chunks < 1 or max_parallel_chunks > count:
            raise ValueError(f"Concurrent chunk limit for {node['name']} must be between 1 and its CPU slot count ({count}).")
        host = config.display_name(node["name"])
        specs.append({"node_id": node_id, "name": config.display_name(node["name"]),
                      "host": host, "cpu_cores": count, "local": local,
                      "available_cpu": available, "chunk_weight": chunk_weight,
                      "max_parallel_chunks": max_parallel_chunks,
                      "free_gpus": resources.get("free_gpus") or resources.get("gpus") or []})
    if not specs:
        raise ValueError("Select at least one worker computer.")
    if not any(row["local"] for row in specs):
        raise ValueError("Include this organizer among the workers because it owns the source library.")
    return specs


def _remote_environment(spec, ssh_user, timeout=20):
    destination = f"{ssh_user}@{spec['host']}"
    probe = r'''for p in "$HOME/programs/miniforge3/envs/mn-molecular-design-studio/bin/python3.11" "$HOME/programs/miniforge3/envs/mn-molecular-design-studio/bin/python"; do
  if [ -x "$p" ]; then
    "$p" -c 'import os,pathlib,sys,rdkit,pyarrow,mn_molecular_design_studio; from mn_molecular_design_studio.engine.mn_compute_scheduler.resources import shared_state_dir; print(sys.executable); print(pathlib.Path(mn_molecular_design_studio.__file__).resolve().parent.parent); print(shared_state_dir()); print(os.getenv("MN_COMPUTE_SCHEDULER_CPU_SLOTS", ""))' 2>/dev/null && exit 0
  fi
done
exit 4'''
    result = subprocess.run(["ssh", "-T", "-o", "BatchMode=yes", "-o", "ConnectTimeout=10",
                             destination, "bash -lc " + shlex.quote(probe)], capture_output=True, text=True,
                            timeout=timeout, check=False)
    if result.returncode:
        raise RuntimeError(f"Could not verify the Studio Conda environment on {spec['name']}: "
                           f"{(result.stderr or result.stdout).strip() or 'RDKit/PyArrow environment not found'}")
    lines = [line.strip() for line in result.stdout.splitlines()]
    lines = [line for line in lines if line]
    if len(lines) < 3 or not lines[0].endswith(("/python", "/python3.11")):
        raise RuntimeError(f"Could not identify the Studio Conda environment on {spec['name']}.")
    cpu_slots = lines[3] if len(lines) > 3 else ""
    return lines[0], lines[1], lines[2], cpu_slots


class _Worker:
    def __init__(self, spec, *, worker_script, library_module, temp_root, ssh_user,
                 worker_args=(), extra_files=()):
        self.spec = spec
        self.queue = queue.Queue(maxsize=2)
        self.error = None
        self.thread = None
        self.stderr = tempfile.TemporaryFile(mode="w+b")
        self.remote_dir = None
        self.remote_python = None
        self.destination = None
        self.output_dir = Path(temp_root) / f"worker-{spec['node_id']}"
        if spec["local"]:
            # The local worker creates its output directory itself. Remote
            # workers write elsewhere, so this is only their local fetch target.
            self.output_dir.parent.mkdir(parents=True, exist_ok=True)
        else:
            self.output_dir.mkdir(parents=True)
        self.process = None

        if spec["local"]:
            repository = str(Path(__file__).resolve().parents[2])
            env = dict(os.environ)
            env["PYTHONPATH"] = repository + os.pathsep + env.get("PYTHONPATH", "")
            env["MDS_EXTERNAL_LIBRARIES_MODULE"] = str(library_module)
            command = [sys.executable, str(worker_script), "--output", str(self.output_dir),
                       "--workers", str(spec["cpu_cores"]), "--run-id", uuid.uuid4().hex, *worker_args]
            self.process = subprocess.Popen(command, cwd=repository, env=env, stdin=subprocess.PIPE,
                                            stdout=subprocess.PIPE, stderr=self.stderr, text=True, bufsize=1)
        else:
            self.destination = f"{ssh_user}@{spec['host']}"
            python, repository, scheduler_state, scheduler_cpu_slots = _remote_environment(spec, ssh_user)
            self.remote_python = python
            self.remote_dir = f"/tmp/mn-mds-library-{uuid.uuid4().hex}"
            subprocess.run(["ssh", "-T", "-o", "BatchMode=yes", self.destination,
                            "mkdir", "-m", "700", "-p", self.remote_dir],
                           capture_output=True, text=True, timeout=20, check=True)
            subprocess.run(["scp", "-q", str(worker_script), f"{self.destination}:{self.remote_dir}/worker.py"],
                           capture_output=True, text=True, timeout=30, check=True)
            subprocess.run(["scp", "-q", str(library_module), f"{self.destination}:{self.remote_dir}/external_libraries.py"],
                           capture_output=True, text=True, timeout=30, check=True)
            for source, relative in extra_files:
                relative = Path(relative)
                if relative.is_absolute() or ".." in relative.parts:
                    raise ValueError("Worker support files must use relative paths.")
                destination = f"{self.remote_dir}/{relative.as_posix()}"
                subprocess.run(["ssh", "-T", "-o", "BatchMode=yes", self.destination,
                                "mkdir", "-p", str(Path(destination).parent)],
                               capture_output=True, text=True, timeout=30, check=True)
                subprocess.run(["scp", "-q", str(source), f"{self.destination}:{destination}"],
                               capture_output=True, text=True, timeout=30, check=True)
            remote_worker = f"{self.remote_dir}/worker.py"
            remote_module = f"{self.remote_dir}/external_libraries.py"
            remote_output = f"{self.remote_dir}/output"
            command = ["env", f"PYTHONPATH={repository}",
                       f"MDS_EXTERNAL_LIBRARIES_MODULE={remote_module}",
                       "MDS_LIBRARY_SCHEDULED_WORKER=1", python,
                       remote_worker, "--output", remote_output,
                       "--workers", str(spec["cpu_cores"]), "--run-id", uuid.uuid4().hex, *worker_args]
            command[1:1] = [f"MN_COMPUTE_SCHEDULER_STATE_DIR={scheduler_state}"]
            if scheduler_cpu_slots:
                command[2:2] = [f"MN_COMPUTE_SCHEDULER_CPU_SLOTS={scheduler_cpu_slots}"]
            self.process = subprocess.Popen(
                ["ssh", "-T", "-o", "BatchMode=yes", "-o", "ConnectTimeout=10",
                 self.destination, "bash -lc " + shlex.quote(shlex.join(command))],
                stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=self.stderr, text=True, bufsize=1,
            )
        ready = self.process.stdout.readline().strip()
        if ready != "MDS_LIBRARY_READY":
            try:
                self.process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                self.process.terminate()
                self.process.wait(timeout=5)
            self.stderr.seek(0)
            stderr = self.stderr.read().decode("utf-8", errors="replace").strip()
            self.stderr.close()
            raise RuntimeError(f"Worker on {spec['name']} did not obtain a scheduler CPU allocation: {stderr[-2000:]}")
        self.thread = threading.Thread(target=self._write_input, name=f"library-stream-{spec['name']}", daemon=True)
        self.thread.start()

    def _write_input(self):
        try:
            while True:
                batch = self.queue.get()
                if batch is None:
                    break
                payload = "".join(json.dumps(item, ensure_ascii=False, separators=(",", ":")) + "\n"
                                  for item in batch)
                self.process.stdin.write(payload)
                self.process.stdin.flush()
        except BaseException as exc:
            self.error = exc
        finally:
            try:
                self.process.stdin.close()
            except (AttributeError, OSError):
                pass

    def send(self, batch):
        while True:
            if self.error:
                raise RuntimeError(f"Could not stream records to {self.spec['name']}: {self.error}") from self.error
            if self.process.poll() is not None:
                raise RuntimeError(f"The library worker on {self.spec['name']} exited before input finished.")
            try:
                self.queue.put(batch, timeout=0.2)
                return
            except queue.Full:
                continue

    def close_input(self):
        while self.thread.is_alive():
            if self.error:
                raise RuntimeError(f"Could not stream records to {self.spec['name']}: {self.error}") from self.error
            try:
                self.queue.put(None, timeout=0.2)
                break
            except queue.Full:
                continue
        self.thread.join()

    def wait(self):
        return_code = self.process.wait()
        self.stderr.seek(0)
        stderr = self.stderr.read().decode("utf-8", errors="replace").strip()
        self.stderr.close()
        if return_code:
            raise RuntimeError(f"Library worker on {self.spec['name']} failed (exit {return_code}): {stderr[-4000:]}")
        if self.error:
            raise RuntimeError(f"Could not stream records to {self.spec['name']}: {self.error}")

    def fetch(self):
        if self.spec["local"]:
            return self.output_dir
        subprocess.run(["scp", "-q", "-r", f"{self.destination}:{self.remote_dir}/output/.", str(self.output_dir)],
                       capture_output=True, text=True, timeout=None, check=True)
        return self.output_dir

    def abort(self):
        if self.remote_dir and self.remote_python and self.process and self.process.poll() is None:
            # Centroid workers handle SIGTERM by terminating their child pool
            # and releasing the remote scheduler lease before SSH is closed.
            script = ("import os,pathlib,signal,sys; root=pathlib.Path(sys.argv[1]); "
                      "f=root/'worker.pid'; pid=int(f.read_text()) if f.exists() else 0; "
                      "cmd=pathlib.Path('/proc')/str(pid)/'cmdline'; "
                      "os.kill(pid,signal.SIGTERM) if pid>1 and cmd.exists() "
                      "and str(root).encode() in cmd.read_bytes() else None")
            try:
                subprocess.run(["ssh", "-T", "-o", "BatchMode=yes", "-o", "ConnectTimeout=5",
                                self.destination, shlex.join([self.remote_python, "-c", script, self.remote_dir])],
                               capture_output=True, text=True, timeout=10, check=False)
                self.process.wait(timeout=5)
            except (OSError, subprocess.SubprocessError):
                pass
        if self.process and self.process.poll() is None:
            self.process.terminate()
            try:
                self.process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                self.process.kill()
                self.process.wait()
        if self.thread and self.thread.is_alive():
            try:
                self.queue.put_nowait(None)
            except queue.Full:
                pass
            self.thread.join(timeout=2)
        if self.stderr and not self.stderr.closed:
            self.stderr.close()


def _merge(worker_dirs, output, *, record, worker_specs, input_rows, selection):
    import pyarrow as pa
    import pyarrow.parquet as pq
    from mn_molecular_design_studio.external_libraries import (
        FEATURE_SCHEMA, FingerprintSpec, atomic_json, check_sources, digest_file, digest_json, _now,
    )

    spec = FingerprintSpec()
    schema = pa.schema([("row_number", pa.uint64()), ("compound_id", pa.string()),
                        ("smiles", pa.string()), ("reason", pa.string())])
    stage = output.with_name(output.name + f".staging-{uuid.uuid4().hex}")
    stage.mkdir(parents=True)
    reject_writer = pq.ParquetWriter(stage / "rejected.parquet", schema, compression="zstd")
    shards = []
    row_count = rejected_rows = 0
    shard_index = 0
    worker_summaries = []
    try:
        for worker_index, worker_dir in enumerate(worker_dirs):
            worker_manifest = json.loads((worker_dir / "worker.json").read_text(encoding="utf-8"))
            if worker_manifest.get("status") != "completed":
                raise RuntimeError(f"Worker output is incomplete: {worker_dir}")
            worker_rows = int(worker_manifest["row_count"])
            worker_rejected = int(worker_manifest["rejected_rows"])
            row_count += worker_rows
            rejected_rows += worker_rejected
            worker_shards = worker_manifest["shards"]
            worker_output_bytes = sum(int(item.get("size_bytes") or 0) for item in worker_shards)
            for item in worker_manifest["shards"]:
                source = (worker_dir / item["path"]).resolve(strict=True)
                if not source.is_relative_to(worker_dir.resolve()) or digest_file(source) != item["sha256"]:
                    raise ValueError("A worker fingerprint shard failed its integrity check.")
                table = pq.read_table(source)
                if table.schema != pa.schema([("compound_key", pa.uint64()), ("compound_id", pa.string()),
                                              ("compound_name", pa.string()), ("source_smiles", pa.string()),
                                              ("fingerprint", pa.binary(spec.bits // 8))]):
                    raise ValueError("A worker returned an unexpected fingerprint table schema.")
                destination = stage / f"features-{shard_index:06d}.parquet"
                pq.write_table(table, destination, compression="zstd")
                shards.append({"path": destination.name, "rows": table.num_rows,
                               "sha256": digest_file(destination), "size_bytes": destination.stat().st_size})
                shard_index += 1
            reject_path = worker_dir / "rejected.parquet"
            if digest_file(reject_path) != worker_manifest["rejections"]["sha256"]:
                raise ValueError("A worker rejected-row table failed its integrity check.")
            worker_output_bytes += reject_path.stat().st_size
            parquet = pq.ParquetFile(reject_path)
            for batch in parquet.iter_batches(batch_size=50_000, use_threads=False):
                reject_writer.write_batch(batch)
            worker = worker_specs[worker_index]
            worker_summaries.append({
                "node_id": worker["node_id"], "name": worker["name"],
                "cpu_worker_count": worker["cpu_cores"], "status": "completed",
                "input_rows": int(worker_manifest["input_rows"]), "row_count": worker_rows,
                "rejected_rows": worker_rejected, "shard_count": len(worker_shards),
                "output_size_bytes": worker_output_bytes,
            })
        reject_writer.close()

        # Validate IDs across every worker without holding the library in RAM.
        id_db = sqlite3.connect(stage / ".id-index.sqlite3")
        id_db.execute("CREATE TABLE ids (id TEXT PRIMARY KEY)")
        try:
            for item in shards:
                parquet = pq.ParquetFile(stage / item["path"])
                for batch in parquet.iter_batches(batch_size=50_000, columns=["compound_id"], use_threads=False):
                    for compound_id in batch.column(0).to_pylist():
                        try:
                            id_db.execute("INSERT INTO ids VALUES (?)", (compound_id,))
                        except sqlite3.IntegrityError as exc:
                            raise ValueError(f"Duplicate compound ID across workers: {compound_id}") from exc
                id_db.commit()
            # Invalid SMILES retain IDs in the rejection table; check those too.
            rejected_table = pq.ParquetFile(stage / "rejected.parquet")
            for batch in rejected_table.iter_batches(batch_size=50_000, columns=["compound_id"], use_threads=False):
                for compound_id in batch.column(0).to_pylist():
                    if not compound_id:
                        continue
                    try:
                        id_db.execute("INSERT INTO ids VALUES (?)", (compound_id,))
                    except sqlite3.IntegrityError as exc:
                        raise ValueError(f"Duplicate compound ID across workers: {compound_id}") from exc
                id_db.commit()
        finally:
            id_db.close()
        (stage / ".id-index.sqlite3").unlink(missing_ok=True)

        if row_count == 0:
            raise ValueError("The selected library rows contain no valid compounds.")
        check_sources(record)
        manifest = {
            "schema_version": FEATURE_SCHEMA, "status": "completed", "library_id": record["id"],
            "created_at": _now(), "source_descriptor_sha256": digest_json({
                "sources": record["sources"], "id_column": record["id_column"],
                "name_column": record.get("name_column", ""), "smiles_column": record.get("smiles_column", "")}),
            "source_verification": "size_and_mtime_only", "fingerprint_spec": spec.metadata(),
            "input_rows": input_rows, "row_count": row_count, "rejected_rows": rejected_rows,
            "cpu_parallelism": "distributed_process_pool",
            "cpu_worker_count": sum(item["cpu_cores"] for item in worker_specs),
            "workers": worker_summaries,
            "selection": {**selection, "rows_processed": input_rows},
            "rejections": {"path": "rejected.parquet", "sha256": digest_file(stage / "rejected.parquet")},
            "compression": "zstd", "table_format": "parquet",
            "row_fields": ["compound_id", "compound_name", "source_smiles", "fingerprint"],
            "shards": shards,
        }
        atomic_json(stage / "manifest.json", manifest)
        if output.exists():
            if not output.is_dir() or any(output.iterdir()):
                raise ValueError("The fingerprint output folder must be new or empty.")
            output.rmdir()
        os.replace(stage, output)
    except BaseException:
        reject_writer.close()
        shutil.rmtree(stage, ignore_errors=True)
        raise
    return output / "manifest.json", row_count, rejected_rows


def prepare_distributed_fingerprints(record, output_directory, worker_nodes, *, max_records=0,
                                     progress_callback=None, sample_parquet=None, ssh_user="compute",
                                     selection_mode="first_rows", sample_seed=42):
    """Stream source rows to one long-lived CPU worker on each selected computer.

    Worker processes remain internal to a single scheduled Studio job. Source
    rows stream to workers while they compute; result files are downloaded only
    after every worker has completed.
    """
    import csv
    import pyarrow as pa
    import pyarrow.parquet as pq
    from mn_molecular_design_studio.external_libraries import check_sources, iter_source_records

    if type(max_records) is not int or max_records < 0:
        raise ValueError("Rows to process must be zero (all rows) or a positive integer.")
    if selection_mode not in {"first_rows", "random_sample"}:
        raise ValueError("Choose first rows or a random sample across all source rows.")
    if type(sample_seed) is not int or not 0 <= sample_seed <= 2**31 - 1:
        raise ValueError("Random sample seed must be a whole number from 0 to 2,147,483,647.")
    if selection_mode == "random_sample" and max_records < 1:
        raise ValueError("A random sample requires a positive sample size.")
    if not _SSH_USER.fullmatch(ssh_user):
        raise ValueError("Use a valid SSH account name for the compute workstations.")
    specs = _node_specs(worker_nodes, ssh_user=ssh_user)
    check_sources(record)
    output = Path(output_directory).expanduser().resolve()
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise ValueError("The fingerprint output folder must be new or empty.")
    if sample_parquet and Path(sample_parquet).exists():
        raise ValueError("The selected-row Parquet sample already exists; choose a new input folder.")
    output.parent.mkdir(parents=True, exist_ok=True)
    worker_script = Path(__file__).with_name("library_worker.py").resolve(strict=True)
    library_module = Path(__file__).parents[1].joinpath("external_libraries.py").resolve(strict=True)
    tmp_root = Path(tempfile.mkdtemp(prefix=".library-workers-", dir=output.parent))
    workers = []
    sample_temp = tmp_root / "first-rows.parquet"
    sample_writer = None
    sample_rows = []
    sample_schema = None
    sample_count = 0
    input_rows = 0
    source_rows_scanned = None
    id_index = sqlite3.connect(tmp_root / "source-id-index.sqlite3")
    id_index.execute("CREATE TABLE ids (id TEXT PRIMARY KEY)")
    try:
        selected_rows = None
        if selection_mode == "random_sample":
            # Reservoir sampling chooses an exact, uniform sample without
            # replacement while keeping only row numbers in memory. A second
            # sequential source pass retrieves the selected records for the
            # distributed fingerprint workers.
            rng = random.Random(sample_seed)
            reservoir = []
            source_rows_scanned = 0
            for source_rows_scanned, _ in enumerate(iter_source_records(record, batch_size=10_000), 1):
                if source_rows_scanned <= max_records:
                    reservoir.append(source_rows_scanned)
                else:
                    replacement = rng.randrange(source_rows_scanned)
                    if replacement < max_records:
                        reservoir[replacement] = source_rows_scanned
                if source_rows_scanned % 100_000 == 0 and progress_callback:
                    progress_callback(source_rows_scanned, {"phase": "random_sampling"})
            selected_rows = sorted(reservoir)
            if not selected_rows:
                raise ValueError("The selected library source contains no records.")
            check_sources(record)

        for spec in specs:
            workers.append(_Worker(spec, worker_script=worker_script, library_module=library_module,
                                   temp_root=tmp_root, ssh_user=ssh_user))
        buffers = [[] for _ in workers]
        batch_limit = 2048
        selected_index = 0
        for row_number, (raw_id, raw_name, raw_smiles, source_error) in enumerate(
                iter_source_records(record, batch_size=10_000), 1):
            if selection_mode == "first_rows" and max_records and row_number > max_records:
                break
            if selected_rows is not None:
                if selected_index >= len(selected_rows):
                    break
                selected_row = selected_rows[selected_index]
                if row_number < selected_row:
                    continue
                if row_number > selected_row:
                    raise RuntimeError("The source changed while its random sample was being read.")
                selected_index += 1
            input_rows += 1
            compound_id = None if raw_id is None or not str(raw_id).strip() else str(raw_id).strip()
            if not source_error:
                if not compound_id:
                    raise ValueError(f"Missing compound ID at source row {row_number}.")
                try:
                    id_index.execute("INSERT INTO ids VALUES (?)", (compound_id,))
                except sqlite3.IntegrityError as exc:
                    raise ValueError(f"Duplicate compound ID at source row {row_number}: {compound_id}") from exc
            smiles = None if raw_smiles is None else str(raw_smiles).strip()
            name = None if raw_name is None else str(raw_name)
            item = {"row_number": row_number, "compound_id": compound_id, "compound_name": name,
                    "source_smiles": smiles, "source_error": source_error}
            worker_index = (input_rows - 1) % len(workers)
            buffers[worker_index].append(item)
            if len(buffers[worker_index]) >= batch_limit:
                workers[worker_index].send(buffers[worker_index])
                buffers[worker_index] = []

            if sample_parquet:
                if sample_schema is None:
                    columns = [record["id_column"]]
                    if record.get("name_column"):
                        columns.append(record["name_column"])
                    columns.append(record["smiles_column"])
                    sample_schema = pa.schema([(column, pa.string()) for column in columns])
                    sample_writer = pq.ParquetWriter(sample_temp, sample_schema, compression="zstd")
                sample = {record["id_column"]: compound_id, record["smiles_column"]: smiles}
                if record.get("name_column"):
                    sample[record["name_column"]] = name
                sample_rows.append(sample)
                sample_count += 1
                if len(sample_rows) >= 10_000:
                    sample_writer.write_table(pa.Table.from_pylist(sample_rows, schema=sample_schema))
                    sample_rows.clear()
            if input_rows % 10_000 == 0:
                id_index.commit()
                if progress_callback:
                    progress_callback(input_rows, {"phase": "fingerprinting"})
        for worker, pending in zip(workers, buffers):
            if pending:
                worker.send(pending)
        id_index.commit()
        if sample_rows:
            sample_writer.write_table(pa.Table.from_pylist(sample_rows, schema=sample_schema))
        if sample_writer:
            sample_writer.close()
        if input_rows == 0:
            raise ValueError("The selected library source contains no records.")
        if selected_rows is not None and selected_index != len(selected_rows):
            raise RuntimeError("The source changed while its random sample was being read.")

        # Close all input streams and wait for every worker before moving any
        # generated result files back to the organizer.
        for worker in workers:
            worker.close_input()
        for worker in workers:
            worker.wait()
        worker_dirs = [worker.fetch() for worker in workers]
        if progress_callback:
            progress_callback(input_rows, 0)
        if selection_mode == "random_sample":
            selection = {"kind": "uniform_random_sample_without_replacement", "requested_size": max_records,
                         "seed": sample_seed, "source_rows_scanned": source_rows_scanned,
                         "rows_selected": input_rows}
        else:
            selection = {"kind": "first_rows" if max_records else "all_rows", "max_records": max_records}
        manifest_path, count, rejects = _merge(worker_dirs, output, record=record, worker_specs=specs,
                                               input_rows=input_rows, selection=selection)
        if sample_parquet:
            destination = Path(sample_parquet).expanduser().resolve()
            destination.parent.mkdir(parents=True, exist_ok=True)
            if destination.exists():
                raise ValueError("The selected-row Parquet sample already exists; refusing to overwrite it.")
            os.replace(sample_temp, destination)
        return manifest_path, count, rejects
    except BaseException:
        for worker in workers:
            worker.abort()
        raise
    finally:
        id_index.close()
        for worker in workers:
            if worker.remote_dir and worker.destination:
                subprocess.run(["ssh", "-T", "-o", "BatchMode=yes", worker.destination,
                                "rm", "-rf", "--", worker.remote_dir],
                               capture_output=True, text=True, timeout=30, check=False)
        shutil.rmtree(tmp_root, ignore_errors=True)
