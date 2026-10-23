"""Packed CPU implementation of ordered, exact best-first Tanimoto clustering."""
from __future__ import annotations

import ctypes
from fractions import Fraction
import hashlib
import os
from pathlib import Path
import platform
import shutil
import subprocess
import tempfile
import time


def _native_library():
    """Build once per source revision on this workstation; no Python extension ABI."""
    from mn_molecular_design_studio.runtime import app_home

    source = Path(__file__).with_name("native") / "centroid_index.cpp"
    revision = hashlib.sha256(source.read_bytes() + platform.machine().encode()).hexdigest()[:20]
    cache = app_home() / "cache" / "centroid-cpu" / revision
    cache.mkdir(parents=True, exist_ok=True)
    target = cache / "centroid_index.so"
    if not target.exists():
        compiler = shutil.which("c++") or shutil.which("g++")
        if not compiler:
            raise RuntimeError("The CPU centroid backend needs a C++ compiler with OpenMP on this workstation.")
        fd, temporary = tempfile.mkstemp(prefix="build-", suffix=".so", dir=cache)
        os.close(fd)
        try:
            result = subprocess.run(
                [compiler, "-std=c++17", "-O3", "-march=native", "-fopenmp", "-shared", "-fPIC",
                 str(source), "-o", temporary], capture_output=True, text=True, timeout=180,
            )
            if result.returncode:
                raise RuntimeError("Could not build the CPU centroid backend: " + result.stderr[-3000:])
            os.replace(temporary, target)
        finally:
            Path(temporary).unlink(missing_ok=True)
    library = ctypes.CDLL(str(target))
    u64p = ctypes.POINTER(ctypes.c_uint64)
    callback_type = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_uint64, ctypes.c_uint64, ctypes.c_uint64)
    library.studio_centroid_cluster.argtypes = [
        u64p, ctypes.c_uint64, ctypes.c_uint32, u64p, u64p, ctypes.c_uint64, ctypes.c_uint64,
        ctypes.c_uint32, ctypes.c_uint32, ctypes.c_uint32, u64p, u64p, u64p, callback_type,
    ]
    library.studio_centroid_cluster.restype = ctypes.c_int
    library.studio_centroid_error.restype = ctypes.c_char_p
    return library, callback_type, revision


def _digest_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def cluster_packed(fingerprints, order, work_directory, *, similarity_threshold, cpu_cores,
                   progress_callback=None, block_size=1024, weights=None):
    """Use a lossless rare-bit prefix index and parallel immutable query batches.

    The final decisions retain input order, including first-match assignment.
    Only representative indices/counts are retained, never a pairwise matrix.
    """
    import numpy as np

    library, callback_type, revision = _native_library()
    from rdkit import DataStructs
    zero = DataStructs.ExplicitBitVect(64)
    zero_matches = int(DataStructs.TanimotoSimilarity(zero, zero) == 1.0)
    if fingerprints.dtype != np.uint64 or not fingerprints.flags.c_contiguous:
        raise ValueError("The CPU centroid backend requires contiguous packed uint64 fingerprints.")
    if order.dtype != np.uint64 or not order.flags.c_contiguous or len(order) != len(fingerprints):
        raise ValueError("The CPU centroid backend requires one uint64 ordering entry per input row.")
    if type(cpu_cores) is not int or cpu_cores < 1 or type(block_size) is not int or block_size < 1:
        raise ValueError("CPU cores and query batch size must be positive integers.")
    n = len(fingerprints)
    if not n:
        raise ValueError("The prepared fingerprint library is empty.")
    if weights is not None and (weights.dtype != np.uint64 or not weights.flags.c_contiguous
                                or len(weights) != n or np.any(weights == 0)):
        raise ValueError("Cluster population weights must be positive uint64 values, one per input row.")
    work = Path(work_directory)
    representatives = np.memmap(work / "representatives.u64", mode="w+", dtype=np.uint64, shape=(n,))
    sizes = np.memmap(work / "cluster-sizes.u64", mode="w+", dtype=np.uint64, shape=(n,))
    summary = np.zeros(4, dtype=np.uint64)
    threshold = Fraction(str(similarity_threshold)).limit_denominator(2**63 - 1)
    if not threshold.numerator:
        threshold = Fraction(1, 2**63 - 1)
    failures = []
    last_progress = 0.0

    @callback_type
    def progress(rows, count, comparisons):
        nonlocal last_progress
        now = time.monotonic()
        if progress_callback and (rows == n or now - last_progress >= 0.5):
            try:
                progress_callback(rows, count)
                last_progress = now
            except BaseException as exc:
                failures.append(exc)
                return 1
        return 0

    pointer = lambda array: array.ctypes.data_as(ctypes.POINTER(ctypes.c_uint64))
    result = library.studio_centroid_cluster(
        pointer(fingerprints), n, fingerprints.shape[1], pointer(order), pointer(weights) if weights is not None else None,
        threshold.numerator, threshold.denominator, cpu_cores, block_size, zero_matches,
        pointer(representatives), pointer(sizes), pointer(summary), progress,
    )
    if failures:
        raise failures[0]
    if result:
        error = library.studio_centroid_error()
        raise RuntimeError(error.decode() if error else "CPU centroid clustering was interrupted.")
    count = int(summary[0])
    return representatives[:count], sizes[:count], {
        "implementation": "packed_prefix_index_openmp_v1", "native_source_revision": revision,
        "cpu_worker_count": cpu_cores, "query_batch_size": block_size,
        "exact_similarity_evaluations": int(summary[1]), "prefix_postings": int(summary[2]),
        "approximate": False, "zero_fingerprint_similarity": zero_matches,
    }


def cluster_chunk_file(input_path, output_path, similarity_threshold, cpu_cores=1):
    """Cluster one bounded chunk; population weights survive every merge round."""
    import numpy as np
    import pyarrow as pa
    import pyarrow.parquet as pq
    from rdkit import Chem
    from rdkit.Chem import Descriptors

    started = time.monotonic()
    table = pq.read_table(input_path)
    n = table.num_rows
    if not n:
        raise ValueError("Cannot cluster an empty chunk.")
    if "molecular_weight" in table.column_names and table["molecular_weight"].null_count == 0:
        molecular_weights = table["molecular_weight"].to_numpy().astype(np.float64)
    else:
        masses = []
        for compound_id, smiles in zip(table["compound_id"].to_pylist(), table["source_smiles"].to_pylist()):
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError(f"Prepared SMILES became invalid for compound {compound_id}.")
            masses.append(Descriptors.MolWt(mol))
        molecular_weights = np.asarray(masses, dtype=np.float64)
    if not np.all(np.isfinite(molecular_weights)):
        raise ValueError("A chunk contains a non-finite molecular weight.")
    weights = (table["cluster_size"].to_numpy().astype(np.uint64)
               if "cluster_size" in table.column_names else np.ones(n, dtype=np.uint64))
    keys = table["compound_key"].to_numpy()
    ordering = np.asarray(np.lexsort((keys, molecular_weights)), dtype=np.uint64)
    fingerprint_column = table["fingerprint"].combine_chunks()
    width = fingerprint_column.type.byte_width
    packed = np.frombuffer(fingerprint_column.buffers()[1], dtype=np.uint8).reshape(n, width)
    words = (width + 7) // 8
    if width % 8:
        fingerprints = np.zeros((n, words), dtype=np.uint64)
        fingerprints.view(np.uint8)[:, :width] = packed
    else:
        fingerprints = np.ascontiguousarray(packed).view(np.uint64)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".chunk-", dir=output_path.parent) as work:
        reps, populations, metrics = cluster_packed(
            fingerprints, ordering, work, similarity_threshold=similarity_threshold,
            cpu_cores=cpu_cores, weights=weights,
        )
        selected = np.array(reps, copy=True)
        counts = np.array(populations, copy=True)
    if int(counts.sum()) != int(weights.sum()):
        raise RuntimeError("Chunk clustering lost population weights.")
    columns = ["compound_key", "compound_id", "compound_name", "source_smiles", "fingerprint"]
    result = table.select(columns).take(pa.array(selected, type=pa.uint64()))
    result = result.append_column("molecular_weight", pa.array(molecular_weights[selected], type=pa.float64()))
    result = result.append_column("cluster_size", pa.array(counts, type=pa.uint64()))
    pq.write_table(result, output_path, compression="zstd", row_group_size=8192)
    import rdkit
    return {**metrics, "sha256": _digest_file(output_path), "rdkit_version": rdkit.__version__,
            "input_rows": n, "row_count": len(selected),
            "population": int(counts.sum()), "elapsed_seconds": time.monotonic() - started,
            "path": output_path.name}
