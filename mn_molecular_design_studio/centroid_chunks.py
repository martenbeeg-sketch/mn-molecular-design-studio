"""Random, balanced CPU chunk clustering with population-preserving merges."""
from __future__ import annotations

from collections import Counter
import json
import math
from pathlib import Path
import shutil
import subprocess
import time

from mn_molecular_design_studio.external_libraries import (
    _now, atomic_json, digest_file, read_feature_manifest, new_output_directory,
)


def balanced_chunk_sizes(rows, target_size=300_000, chunk_count=0):
    if type(rows) is not int or rows < 1:
        raise ValueError("A prepared library must contain at least one compound.")
    if type(target_size) is not int or target_size < 1:
        raise ValueError("Target chunk size must be a positive integer.")
    if type(chunk_count) is not int or chunk_count < 0 or chunk_count > rows:
        raise ValueError("Chunk count must be zero (automatic) or between one and the input row count.")
    count = chunk_count or (rows + target_size - 1) // target_size
    size, extra = divmod(rows, count)
    return [size + int(i < extra) for i in range(count)]


def _worker_cycle(specs):
    """Build a repeatable chunk assignment using measured relative machine weights."""
    weights = [spec.get("chunk_weight", spec["cpu_cores"]) for spec in specs]
    scores = [0] * len(weights)
    result = []
    for _ in range(sum(weights)):
        scores = [value + weight for value, weight in zip(scores, weights)]
        chosen = max(range(len(scores)), key=scores.__getitem__)
        result.append(chosen)
        scores[chosen] -= sum(weights)
    return result


def representative_quotas(counts, budget, seed=42):
    """Allocate a representative cap in proportion to each chunk's cluster count.

    Hamilton's largest-remainder apportionment keeps the sampling fraction
    close to uniform across local clusters. When the cap can cover every
    nonempty chunk, a chunk that would otherwise receive zero gets one
    representative, borrowed from a chunk with more than one.
    """
    import numpy as np

    if type(budget) is not int or budget < 1 or any(type(count) is not int or count < 0 for count in counts):
        raise ValueError("The final comparison budget must be a positive integer.")
    quotas = [0] * len(counts)
    target = min(budget, sum(counts))
    active = [index for index, count in enumerate(counts) if count]
    if target == sum(counts):
        return list(counts)

    total = sum(counts)
    remainders = []
    allocated = 0
    for index, count in enumerate(counts):
        whole, remainder = divmod(target * count, total) if total else (0, 0)
        quotas[index] = whole
        allocated += whole
        remainders.append((remainder, index))
    tie_order = np.random.default_rng(seed).permutation(len(counts))
    tie_rank = {int(index): rank for rank, index in enumerate(tie_order)}
    leftover = target - allocated
    for _, index in sorted(remainders, key=lambda item: (-item[0], tie_rank[item[1]]))[:leftover]:
        quotas[index] += 1

    if target >= len(active):
        missing = [index for index in active if quotas[index] == 0]
        donors = sorted(active, key=lambda index: (-quotas[index], tie_rank[index]))
        for index in missing:
            donor = next(candidate for candidate in donors if quotas[candidate] > 1)
            quotas[donor] -= 1
            quotas[index] = 1
            donors = sorted(active, key=lambda item: (-quotas[item], tie_rank[item]))
    return quotas


def _merge_stop_decision(input_rows, representative_rows, source_rows, *,
                         target_reduction_factor, plateau_reduction_percent,
                         plateau_consecutive_rounds, plateau_streak):
    """Return the next plateau streak and an optional hierarchical stop reason."""
    if representative_rows == input_rows:
        return plateau_streak, "no_reduction_in_round"

    reduction_percent = 100.0 * (input_rows - representative_rows) / input_rows
    if plateau_reduction_percent > 0 and reduction_percent < plateau_reduction_percent:
        plateau_streak += 1
    else:
        plateau_streak = 0

    if target_reduction_factor and representative_rows * target_reduction_factor <= source_rows:
        return plateau_streak, "target_reduction_factor_reached"
    if plateau_reduction_percent > 0 and plateau_streak >= plateau_consecutive_rounds:
        return plateau_streak, "plateau_reduction_below_threshold"
    return plateau_streak, None


def _sample_candidates(files, budget, seed, destination):
    import numpy as np
    import pyarrow as pa
    import pyarrow.parquet as pq

    counts = [pq.ParquetFile(path).metadata.num_rows for path in files]
    quotas = representative_quotas(counts, budget, seed)
    writer = None
    selected_population = 0
    try:
        for index, (path, count, quota) in enumerate(zip(files, counts, quotas)):
            if not quota:
                continue
            rng = np.random.default_rng(np.random.SeedSequence([seed, index, 1]))
            selected = np.sort(rng.choice(count, size=quota, replace=False))
            offset = 0
            for batch in pq.ParquetFile(path).iter_batches(batch_size=8192):
                start, end = np.searchsorted(selected, [offset, offset + batch.num_rows])
                if end > start:
                    table = pa.Table.from_batches([batch]).take(pa.array(selected[start:end] - offset))
                    selected_population += sum(table["cluster_size"].to_pylist())
                    column = table.schema.get_field_index("cluster_size")
                    # Final cluster sizes count sampled representatives. Discarded
                    # local clusters are never reported as covered source data.
                    table = table.set_column(column, "cluster_size", pa.array(np.ones(end - start, dtype=np.uint64)))
                    if writer is None:
                        writer = pq.ParquetWriter(destination, table.schema, compression="zstd")
                    writer.write_table(table)
                offset += batch.num_rows
    finally:
        if writer:
            writer.close()
    return sum(quotas), selected_population, quotas


def _shuffle_resources(requested_threads):
    """Choose a bounded DuckDB memory budget and matching sort parallelism.

    The randomization query is external-sort capable, but too many concurrent
    sort threads can exhaust a small per-job DuckDB memory cap before spilling.
    Use available host/cgroup memory, cap the sort at 64 GiB, and scale its
    thread count to the selected memory budget.
    """
    import os

    available = None
    try:
        with open("/proc/meminfo", encoding="ascii") as handle:
            for line in handle:
                if line.startswith("MemAvailable:"):
                    available = int(line.split()[1]) * 1024
                    break
    except (OSError, ValueError):
        pass
    if available is None:
        try:
            available = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_AVPHYS_PAGES")
        except (AttributeError, OSError, ValueError):
            available = 8 * 1024**3

    # Respect cgroup v2 or v1 limits when Studio runs inside a container.
    for limit_path, current_path in (
        ("/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/memory.current"),
        ("/sys/fs/cgroup/memory/memory.limit_in_bytes",
         "/sys/fs/cgroup/memory/memory.usage_in_bytes"),
    ):
        try:
            raw_limit = Path(limit_path).read_text(encoding="ascii").strip()
            if raw_limit == "max":
                continue
            limit = int(raw_limit)
            if 0 < limit < 1 << 60:
                current = int(Path(current_path).read_text(encoding="ascii").strip())
                available = min(available, max(0, limit - current))
        except (OSError, ValueError):
            continue

    minimum = 512 * 1024**2
    if available < minimum:
        raise MemoryError("At least 512 MiB of available memory is required to randomize this library.")
    memory_bytes = min(64 * 1024**3, int(available * 0.70))
    threads = min(requested_threads, max(1, memory_bytes // (1024**3)))
    return max(1, int(threads)), memory_bytes


def _shuffle_chunks(files, sizes, directory, *, seed, cpu_cores, fingerprint_bytes, progress,
                    sample_rows=None):
    """An external sort makes globally randomized, exactly balanced chunks.

    DuckDB spills the shuffle to this job's temporary disk instead of storing
    a full-library permutation or the library in the Studio Database.
    """
    import duckdb
    import pyarrow as pa
    import pyarrow.parquet as pq

    directory.mkdir(parents=True)
    schema = pa.schema([
        ("compound_key", pa.uint64()), ("compound_id", pa.string()), ("compound_name", pa.string()),
        ("source_smiles", pa.string()), ("fingerprint", pa.binary(fingerprint_bytes)),
        ("molecular_weight", pa.float64()), ("cluster_size", pa.uint64()),
    ])
    columns = pq.ParquetFile(files[0]).schema_arrow.names
    mass = "molecular_weight" if "molecular_weight" in columns else "NULL::DOUBLE AS molecular_weight"
    population = "cluster_size" if "cluster_size" in columns else "1::UBIGINT AS cluster_size"
    connection = duckdb.connect()
    writer = None
    chunks = []
    seen = offset = index = 0
    try:
        sort_threads, memory_bytes = _shuffle_resources(cpu_cores)
        connection.execute("SET threads = ?", [sort_threads])
        connection.execute("SET memory_limit = ?", [f"{memory_bytes / 1024**3:.2f}GB"])
        connection.execute("SET preserve_insertion_order = false")
        connection.execute("SET temp_directory = ?", [str(directory / "sort-spill")])
        progress(0, 0, "Randomizing compounds on temporary disk")
        query = (
            f"SELECT compound_key,compound_id,compound_name,source_smiles,fingerprint,{mass},{population} "
            "FROM read_parquet(?) ORDER BY hash(compound_key, ?), compound_key"
        )
        parameters = [[str(path) for path in files], seed]
        if sample_rows is not None:
            query += " LIMIT ?"
            parameters.append(sample_rows)
        reader = connection.execute(query, parameters).fetch_record_batch(8192)
        for batch in reader:
            table = pa.Table.from_batches([batch]).cast(schema)
            position = 0
            while position < table.num_rows:
                if index >= len(sizes):
                    raise ValueError("Prepared file row count exceeds its manifest.")
                if writer is None:
                    path = directory / f"input-{index:08d}.parquet"
                    chunks.append(path)
                    writer = pq.ParquetWriter(path, schema, compression="zstd")
                take = min(sizes[index] - offset, table.num_rows - position)
                writer.write_table(table.slice(position, take))
                seen += take
                offset += take
                position += take
                if offset == sizes[index]:
                    writer.close()
                    writer = None
                    index += 1
                    offset = 0
            progress(seen, 0, f"Preparing random chunks: {seen:,} compounds")
        if seen != sum(sizes) or writer is not None:
            raise ValueError("Prepared file row count disagrees with its manifest.")
        return chunks
    finally:
        if writer:
            writer.close()
        connection.close()
        shutil.rmtree(directory / "sort-spill", ignore_errors=True)


def _run_round(chunks, specs, directory, threshold, ssh_user, progress):
    import pyarrow.parquet as pq
    from mn_molecular_design_studio.cluster.distributed_library import _Worker

    cycle = _worker_cycle(specs)
    assignment = [cycle[i % len(cycle)] for i in range(len(chunks))]
    tasks = Counter(assignment)
    root = Path(__file__).parent
    workers = {}
    summaries, files, metrics = [], [], []
    try:
        for index, spec in enumerate(specs):
            if not tasks[index]:
                continue
            workers[index] = _Worker(
                spec, worker_script=root / "cluster/library_centroid_worker.py",
                library_module=root / "external_libraries.py", temp_root=directory, ssh_user=ssh_user,
                worker_args=["--tasks", str(tasks[index]), "--threshold", str(threshold),
                             "--max-parallel-chunks", str(spec.get("max_parallel_chunks", spec["cpu_cores"]))],
                extra_files=[(root / "centroid_cpu.py", "centroid_cpu.py"),
                             (root / "native/centroid_index.cpp", "native/centroid_index.cpp")],
            )
        for index, (path, selected) in enumerate(zip(chunks, assignment)):
            worker = workers[selected]
            input_path = str(path)
            if not worker.spec["local"]:
                input_path = f"{worker.remote_dir}/input-{index:08d}.parquet"
                subprocess.run(["scp", "-q", str(path), f"{worker.destination}:{input_path}"],
                               capture_output=True, text=True, check=True)
            worker.send([{"chunk": index, "input": input_path}])
            progress(index + 1, 0, f"Distributing chunk {index + 1}/{len(chunks)}")
        for worker in workers.values():
            worker.close_input()
        while any(worker.process.poll() is None for worker in workers.values()):
            progress(0, 0, f"Clustering {len(chunks)} chunks on {len(workers)} computers")
            time.sleep(0.5)
        for worker in workers.values():
            worker.wait()
        progress(0, 0, "Collecting completed chunk representatives on the organizer")
        for index, worker in workers.items():
            result_dir = worker.fetch()
            manifest = json.loads((result_dir / "worker-manifest.json").read_text())
            expected = {f"chunk-{chunk:08d}.parquet": pq.ParquetFile(chunks[chunk]).metadata.num_rows
                        for chunk, selected in enumerate(assignment) if selected == index}
            if (len(manifest["chunks"]) != len(expected)
                    or {item["path"] for item in manifest["chunks"]} != set(expected)):
                raise ValueError("A worker returned missing, duplicate, or unexpected centroid chunks.")
            row_count = input_rows = population = 0
            for item in manifest["chunks"]:
                path = (result_dir / item["path"]).resolve(strict=True)
                if (not path.is_relative_to(result_dir.resolve()) or path.suffix != ".parquet"
                        or digest_file(path) != item["sha256"]):
                    raise ValueError("A centroid chunk failed its output integrity check.")
                if pq.ParquetFile(path).metadata.num_rows != item["row_count"]:
                    raise ValueError("A centroid chunk has an unexpected row count.")
                if item["input_rows"] != expected[item["path"]]:
                    raise ValueError("A centroid chunk did not process every assigned input row.")
                files.append(path)
                row_count += item["row_count"]
                input_rows += item["input_rows"]
                population += item["population"]
                metrics.append(item)
            summaries.append({"node_id": worker.spec["node_id"], "name": worker.spec["name"],
                              "task_kind": "centroids",
                              "cpu_worker_count": worker.spec["cpu_cores"], "status": "completed",
                              "chunk_weight": worker.spec.get("chunk_weight", worker.spec["cpu_cores"]),
                              "max_parallel_chunks": worker.spec.get("max_parallel_chunks", worker.spec["cpu_cores"]),
                              "input_rows": input_rows, "row_count": row_count, "population": population,
                              "chunk_count": len(manifest["chunks"]),
                              "parallel_chunks": manifest["parallel_chunks"],
                              "cores_per_chunk": manifest["cores_per_chunk"]})
        if len(files) != len(chunks):
            raise ValueError("A worker did not return every assigned chunk.")
        return sorted(files, key=lambda path: path.name), summaries, metrics
    except BaseException:
        for worker in workers.values():
            worker.abort()
        raise
    finally:
        for worker in workers.values():
            if worker.remote_dir and worker.destination:
                subprocess.run(["ssh", "-T", "-o", "BatchMode=yes", worker.destination,
                                "rm", "-rf", "--", worker.remote_dir],
                               capture_output=True, text=True, timeout=30, check=False)


def prepare_chunked_centroids(feature_manifest, output, *, worker_nodes, similarity_threshold=0.6,
                              chunk_size=300_000, chunk_count=0, max_merge_rounds=0,
                              random_seed=42, ssh_user="compute", progress_callback=None,
                              method="sampled", final_candidate_count=300_000,
                              source_sample_size=0, target_reduction_factor=10,
                              plateau_reduction_percent=1.0, plateau_consecutive_rounds=3):
    import duckdb
    import pyarrow.parquet as pq
    from mn_molecular_design_studio.cluster.distributed_library import _node_specs

    started = time.monotonic()
    source_path, source = read_feature_manifest(feature_manifest)
    source_rows = source["row_count"]
    if type(source_sample_size) is not int or source_sample_size < 0:
        raise ValueError("The source sample size must be zero (use the full library) or a positive integer.")
    if source_sample_size > source_rows:
        raise ValueError("The source sample size cannot exceed the prepared library size.")
    n = source_sample_size or source_rows
    balanced_chunk_sizes(n, chunk_size, chunk_count)
    if method not in {"sampled", "hierarchical"}:
        raise ValueError("Choose sampled final clustering or hierarchical representative merging.")
    if type(final_candidate_count) is not int or final_candidate_count < 1:
        raise ValueError("Final comparison size must be a positive integer.")
    if not math.isfinite(similarity_threshold) or not 0 < similarity_threshold <= 1:
        raise ValueError("Tanimoto similarity must be in (0, 1].")
    if type(max_merge_rounds) is not int or max_merge_rounds < 0:
        raise ValueError("Merge rounds must be zero (until stable) or a positive integer.")
    if type(target_reduction_factor) is not int or target_reduction_factor < 0 or target_reduction_factor == 1:
        raise ValueError("Target reduction factor must be zero (disabled) or a whole number of at least 2.")
    if (isinstance(plateau_reduction_percent, bool)
            or not isinstance(plateau_reduction_percent, (int, float))
            or not math.isfinite(plateau_reduction_percent)
            or not 0 <= plateau_reduction_percent <= 100):
        raise ValueError("Plateau reduction threshold must be between 0 and 100 percent.")
    if type(plateau_consecutive_rounds) is not int or plateau_consecutive_rounds < 1:
        raise ValueError("Plateau rounds must be a positive whole number.")
    if type(random_seed) is not int or not 0 <= random_seed < 2**63:
        raise ValueError("Use a nonnegative random seed below 2^63.")
    specs = _node_specs(worker_nodes, ssh_user=ssh_user)
    original_files = []
    for shard in source["shards"]:
        path = (source_path.parent / shard["path"]).resolve(strict=True)
        if not path.is_relative_to(source_path.parent) or digest_file(path) != shard["sha256"]:
            raise ValueError("A prepared fingerprint shard has changed or escapes its manifest directory.")
        original_files.append(path)
    output = new_output_directory(output)
    work = output / ".centroid-work"
    work.mkdir()
    last_progress = 0.0

    def progress(rows, representatives, phase):
        nonlocal last_progress
        now = time.monotonic()
        if progress_callback and now - last_progress >= 0.5:
            progress_callback(rows, representatives, phase)
            last_progress = now

    files = original_files
    rounds = []
    workers = {}
    comparisons = 0
    previous_round = None
    writer = None
    selected_population = n
    candidate_count = None
    quotas = None
    plateau_streak = 0
    try:
        round_index = 0
        while True:
            sizes = balanced_chunk_sizes(n, chunk_size, chunk_count if not round_index else 0)
            round_dir = work / f"round-{round_index:03d}"
            round_dir.mkdir()
            chunks = _shuffle_chunks(files, sizes, round_dir / "inputs", seed=(random_seed + round_index) % 2**63,
                                     cpu_cores=next(s["cpu_cores"] for s in specs if s["local"]),
                                     fingerprint_bytes=source["fingerprint_spec"]["bits"] // 8,
                                     progress=progress,
                                     sample_rows=n if round_index == 0 and n < source_rows else None)
            outputs, summaries, metrics = _run_round(chunks, specs, round_dir / "workers",
                                                     similarity_threshold, ssh_user, progress)
            reduced = sum(item["row_count"] for item in metrics)
            if sum(item["input_rows"] for item in metrics) != n:
                raise ValueError("Chunk workers did not process every input representative.")
            if sum(item["population"] for item in metrics) != (source_sample_size or source_rows):
                raise ValueError("Centroid merging lost original compound populations.")
            round_checks = sum(item["exact_similarity_evaluations"] for item in metrics)
            comparisons += round_checks
            rounds.append({"round": round_index, "input_rows": n, "representatives": reduced,
                           "chunk_count": len(chunks), "chunk_rows_min": min(sizes), "chunk_rows_max": max(sizes),
                           "reduction_percent": 100.0 * (n - reduced) / n,
                           "exact_similarity_evaluations": round_checks, "workers": summaries})
            for item in summaries:
                aggregate = workers.setdefault(item["node_id"], {**item, "input_rows": 0, "row_count": 0, "chunk_count": 0})
                for key in ("input_rows", "row_count", "chunk_count"):
                    aggregate[key] += item[key]
            files = outputs
            shutil.rmtree(round_dir / "inputs")
            if previous_round:
                shutil.rmtree(previous_round)
            previous_round = round_dir
            if method == "sampled":
                progress(0, reduced, "Sampling representatives in proportion to chunk cluster counts")
                candidate_path = work / "final-candidates.parquet"
                candidate_count, selected_population, quotas = _sample_candidates(
                    files, final_candidate_count, random_seed, candidate_path,
                )
                final_dir = work / "final-comparison"
                final_files, final_workers, final_metrics = _run_round(
                    [candidate_path], specs, final_dir, similarity_threshold, ssh_user, progress,
                )
                final_checks = sum(item["exact_similarity_evaluations"] for item in final_metrics)
                comparisons += final_checks
                if sum(item["population"] for item in final_metrics) != candidate_count:
                    raise ValueError("Final clustering lost sampled representatives.")
                rounds.append({"round": 1, "kind": "final_sampled_comparison", "input_rows": candidate_count,
                               "representatives": sum(item["row_count"] for item in final_metrics),
                               "chunk_count": 1, "chunk_rows_min": candidate_count, "chunk_rows_max": candidate_count,
                               "exact_similarity_evaluations": final_checks, "workers": final_workers})
                for item in final_workers:
                    aggregate = workers.setdefault(item["node_id"], {**item, "input_rows": 0, "row_count": 0, "chunk_count": 0})
                    for key in ("input_rows", "row_count", "chunk_count"):
                        aggregate[key] += item[key]
                files = final_files
                stop_reason = "final_sampled_comparison_completed"
                break
            if len(chunks) == 1:
                stop_reason = "one_final_chunk"
                break
            plateau_streak, criterion_stop = _merge_stop_decision(
                n, reduced, source_sample_size or source_rows,
                target_reduction_factor=target_reduction_factor,
                plateau_reduction_percent=plateau_reduction_percent,
                plateau_consecutive_rounds=plateau_consecutive_rounds,
                plateau_streak=plateau_streak,
            )
            if criterion_stop:
                stop_reason = criterion_stop
                break
            if max_merge_rounds and round_index >= max_merge_rounds:
                stop_reason = "merge_round_limit"
                break
            n = reduced
            round_index += 1
        destination = output / "centroids.parquet"
        final_rows = population = 0
        for path in files:
            for batch in pq.ParquetFile(path).iter_batches(batch_size=8192):
                if writer is None:
                    writer = pq.ParquetWriter(destination, batch.schema, compression="zstd")
                writer.write_batch(batch)
                final_rows += batch.num_rows
                population += sum(batch.column(batch.schema.get_field_index("cluster_size")).to_pylist())
        writer.close()
        writer = None
        expected_population = candidate_count if method == "sampled" else (source_sample_size or source_rows)
        if final_rows != rounds[-1]["representatives"] or population != expected_population:
            raise ValueError("The final centroid table failed its population or row-count check.")
        manifest = {
            "schema_version": "studio.library-centroids.v2", "status": "completed", "created_at": _now(),
            "library_id": source["library_id"], "source_feature_manifest": str(source_path),
            "source_feature_manifest_sha256": digest_file(source_path), "fingerprint_spec": source["fingerprint_spec"],
            "algorithm": "random_balanced_sampled_best_first" if method == "sampled" else "random_balanced_hierarchical_best_first",
            "implementation": "packed_prefix_index_openmp_v1",
            "algorithm_provenance": {"rule": "OpenVS best-first sphere exclusion; smallest molecular weight first",
                                     "reference": "https://github.com/gfzhou/OpenVS/blob/3d46061/openvs/utils/cluster.py"},
            "similarity_threshold": similarity_threshold, "similarity_scope": "each_chunk_at_each_merge_round",
            "approximate": True if method == "sampled" else len(rounds) > 1 or rounds[0]["chunk_count"] > 1,
            "representative_order": ["molecular_weight_ascending_within_chunk", "source_compound_key_ascending"],
            "input_rows": source_sample_size or source_rows, "source_library_rows": source_rows,
            "source_sample_size": source_sample_size or None,
            "source_sample_method": "seeded_hash_order_prefix" if source_sample_size else None,
            "row_count": final_rows, "population": population,
            "cluster_size_scope": "sampled_chunk_representatives" if method == "sampled" else "original_source_compounds",
            "final_candidate_budget": final_candidate_count if method == "sampled" else None,
            "candidate_count": candidate_count, "chunk_representative_quotas": quotas,
            "chunk_quota_policy": "proportional_to_local_cluster_count_with_one_per_nonempty_chunk_when_budget_allows" if method == "sampled" else None,
            "selected_source_cluster_population": selected_population,
            "coverage_guarantee": False,
            "exact_similarity_evaluations": comparisons, "cpu_worker_count": sum(s["cpu_cores"] for s in specs),
            "workers": list(workers.values()), "random_seed": random_seed,
            "randomization": "DuckDB hash(compound_key, seed + round) external sort",
            "duckdb_version": duckdb.__version__, "target_chunk_size": chunk_size,
            "requested_chunk_count": chunk_count, "max_merge_rounds": max_merge_rounds,
            "stopping_criteria": ({
                "semantics": "stop_when_any_enabled_criterion_is_met",
                "target_reduction_factor": target_reduction_factor or None,
                "plateau_reduction_percent": plateau_reduction_percent or None,
                "plateau_consecutive_rounds": plateau_consecutive_rounds if plateau_reduction_percent else None,
                "final_reduction_factor": (source_sample_size or source_rows) / final_rows,
            } if method == "hierarchical" else None),
            "rounds": rounds, "stop_reason": stop_reason, "limits": None,
            "elapsed_seconds": time.monotonic() - started,
            "retention": "final_centroids_and_manifest_only; temporary_chunks_removed",
            "files": [{"path": "centroids.parquet", "sha256": digest_file(destination),
                       "size_bytes": destination.stat().st_size}],
        }
        atomic_json(output / "manifest.json", manifest)
        if progress_callback:
            progress_callback(source_sample_size or source_rows, final_rows, "Centroid library completed")
        return output / "manifest.json"
    except BaseException as exc:
        atomic_json(output / "failure.json", {"status": "failed", "error": str(exc), "rounds": rounds})
        raise
    finally:
        if writer:
            writer.close()
        shutil.rmtree(work, ignore_errors=True)
