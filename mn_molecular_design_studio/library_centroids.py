"""Disk-backed, exact greedy Tanimoto representative selection."""
from __future__ import annotations

import math
import sqlite3

from mn_molecular_design_studio.external_libraries import (
    _now, atomic_json, digest_file, iter_feature_batches, new_output_directory,
    read_feature_manifest,
)


def prepare_centroids(feature_manifest, output, *, similarity_threshold=0.6,
                      progress_callback=None, comparison_batch_size=8192):
    """Select smallest-molecular-weight greedy Tanimoto representatives.

    The full library and representatives stay on disk. Similarity
    checks are exact; a fingerprint-popcount bound skips pairs that cannot reach
    the requested threshold. There are no row, representative, or comparison
    budgets. The output contains centroid compounds and cluster sizes, not a
    per-input assignment table. Runtime and disk use depend on input and hardware.
    """
    import pyarrow as pa
    import pyarrow.parquet as pq
    from rdkit import Chem, DataStructs
    from rdkit.Chem import Descriptors

    path, source = read_feature_manifest(feature_manifest)
    if isinstance(similarity_threshold, bool) or not math.isfinite(similarity_threshold) or not 0 < similarity_threshold <= 1:
        raise ValueError("Tanimoto similarity must be in (0, 1].")
    if type(comparison_batch_size) is not int or comparison_batch_size < 1:
        raise ValueError("Comparison batch size must be a positive integer.")
    output = new_output_directory(output)
    index = sqlite3.connect(output / ".centroid-sort.sqlite3")
    index.execute("PRAGMA temp_store=FILE")
    index.execute("PRAGMA cache_size=-65536")
    index.execute("CREATE TABLE compounds (compound_key INTEGER PRIMARY KEY, compound_id TEXT, compound_name TEXT, smiles TEXT, fp BLOB, popcount INTEGER, molecular_weight REAL)")
    index.execute("CREATE TABLE representatives (representative_key INTEGER PRIMARY KEY, compound_key INTEGER, compound_id TEXT, compound_name TEXT, smiles TEXT, fp BLOB, popcount INTEGER, molecular_weight REAL, cluster_size INTEGER)")
    index.execute("CREATE INDEX representatives_popcount ON representatives (popcount, representative_key)")

    centroid_schema = pa.schema([
        ("compound_key", pa.uint64()), ("compound_id", pa.string()),
        ("compound_name", pa.string()), ("source_smiles", pa.string()),
        ("fingerprint", pa.binary(source["fingerprint_spec"]["bits"] // 8)),
        ("molecular_weight", pa.float64()), ("cluster_size", pa.uint64()),
    ])
    centroids_path = output / "centroids.parquet"
    centroid_writer = None
    loaded = clustered = representative_count = comparisons = 0
    try:
        # The external SQLite index lets fingerprint batches and sorting spill
        # to the data-owning disk instead of keeping a whole library in RAM.
        for batch in iter_feature_batches(path):
            compound_rows = []
            for row in batch.to_pylist():
                mol = Chem.MolFromSmiles(row["source_smiles"])
                if mol is None:
                    raise ValueError(f"Prepared SMILES became invalid for compound {row['compound_id']}.")
                fp = DataStructs.CreateFromBinaryText(row["fingerprint"])
                compound_rows.append((
                    int(row["compound_key"]), row["compound_id"], row["compound_name"],
                    row["source_smiles"], row["fingerprint"], fp.GetNumOnBits(),
                    float(Descriptors.MolWt(mol)),
                ))
            index.executemany("INSERT INTO compounds VALUES (?,?,?,?,?,?,?)", compound_rows)
            index.commit()
            loaded += len(compound_rows)
            if progress_callback:
                progress_callback(loaded, 0)
        if loaded != source["row_count"]:
            raise ValueError("Fingerprint manifest row count disagrees with its shards.")
        index.execute("CREATE INDEX compound_order ON compounds (molecular_weight, compound_key)")
        index.commit()

        for key, compound_id, compound_name, smiles, packed, popcount, molecular_weight in index.execute(
            "SELECT compound_key, compound_id, compound_name, smiles, fp, popcount, molecular_weight "
            "FROM compounds ORDER BY molecular_weight, compound_key"
        ):
            clustered += 1
            fp = DataStructs.CreateFromBinaryText(packed)
            if fp.GetNumBits() != source["fingerprint_spec"]["bits"]:
                raise ValueError("Fingerprint length disagrees with the feature manifest.")

            # If one fingerprint's on-bit count is too large or too small, even
            # perfect overlap cannot reach the requested Tanimoto threshold.
            if popcount:
                min_popcount = math.ceil(similarity_threshold * popcount - 1e-12)
                max_popcount = math.floor(popcount / similarity_threshold + 1e-12)
            else:
                min_popcount = max_popcount = 0

            match = None
            candidates = index.execute(
                "SELECT representative_key, fp FROM representatives "
                "WHERE popcount BETWEEN ? AND ? ORDER BY representative_key",
                (min_popcount, max_popcount),
            )
            while match is None:
                candidate_rows = candidates.fetchmany(comparison_batch_size)
                if not candidate_rows:
                    break
                candidate_fps = [DataStructs.CreateFromBinaryText(row[1]) for row in candidate_rows]
                values = DataStructs.BulkTanimotoSimilarity(fp, candidate_fps)
                comparisons += len(candidate_rows)
                found = next((i for i, value in enumerate(values) if value >= similarity_threshold), None)
                if found is not None:
                    match = candidate_rows[found][0]

            if match is None:
                match = representative_count
                representative_count += 1
                index.execute(
                    "INSERT INTO representatives VALUES (?,?,?,?,?,?,?,?,1)",
                    (match, key, compound_id, compound_name, smiles, packed, popcount, molecular_weight),
                )
            else:
                index.execute("UPDATE representatives SET cluster_size=cluster_size+1 WHERE representative_key=?", (match,))
            if clustered % 10_000 == 0:
                index.commit()
                if progress_callback:
                    progress_callback(clustered, representative_count)

        index.commit()
        centroid_writer = pq.ParquetWriter(centroids_path, centroid_schema, compression="zstd")
        cursor = index.execute(
            "SELECT compound_key, compound_id, compound_name, smiles, fp, molecular_weight, cluster_size "
            "FROM representatives ORDER BY representative_key"
        )
        while True:
            rows = cursor.fetchmany(10_000)
            if not rows:
                break
            centroid_writer.write_table(pa.Table.from_pylist([
                {"compound_key": key, "compound_id": compound_id, "compound_name": name,
                 "source_smiles": smiles, "fingerprint": packed, "molecular_weight": weight,
                 "cluster_size": cluster_size}
                for key, compound_id, name, smiles, packed, weight, cluster_size in rows
            ], schema=centroid_schema))
        centroid_writer.close()
        centroid_writer = None
        if clustered != source["row_count"]:
            raise ValueError("Clustered row count disagrees with the fingerprint manifest.")
        if progress_callback:
            progress_callback(clustered, representative_count)
        manifest = {
            "schema_version": "studio.library-centroids.v2", "status": "completed", "created_at": _now(),
            "library_id": source["library_id"], "source_feature_manifest": str(path),
            "source_feature_manifest_sha256": digest_file(path), "fingerprint_spec": source["fingerprint_spec"],
            "algorithm": "exact_greedy_molecular_weight_order",
            "similarity_threshold": similarity_threshold,
            "representative_order": ["molecular_weight_ascending", "source_compound_key_ascending"],
            "input_rows": clustered, "row_count": representative_count,
            "exact_similarity_evaluations": comparisons,
            "limits": None,
            "files": [{"path": "centroids.parquet", "sha256": digest_file(centroids_path),
                       "size_bytes": centroids_path.stat().st_size}],
        }
        atomic_json(output / "manifest.json", manifest)
        return output / "manifest.json"
    except Exception as exc:
        atomic_json(output / "failure.json", {"status": "failed", "error": str(exc),
                                               "input_rows_loaded": loaded,
                                               "rows_clustered": clustered,
                                               "representatives": representative_count})
        raise
    finally:
        if centroid_writer is not None:
            centroid_writer.close()
        index.close()
        (output / ".centroid-sort.sqlite3").unlink(missing_ok=True)
