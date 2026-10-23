"""Metadata-only library catalog and streaming compound fingerprint preparation.

All molecule records and preparation indexes are written to the explicitly
selected external output directory, never to Studio's molecular Database.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import sqlite3
import uuid
from concurrent.futures import ProcessPoolExecutor
from collections import deque
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from multiprocessing import get_context
from pathlib import Path

from mn_molecular_design_studio.runtime import app_home

LIBRARY_SCHEMA = "studio.external-library.v1"
FEATURE_SCHEMA = "studio.library-fingerprints.v2"
CENTROID_SCHEMA = "studio.library-centroids.v2"


def _now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def digest_json(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def digest_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + "." + uuid.uuid4().hex + ".tmp")
    try:
        temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


class LibraryCatalog:
    def __init__(self, root=None):
        from mn_molecular_design_studio.molecular_database import MolecularDatabase
        self.database = MolecularDatabase(Path(root) if root is not None else app_home() / "database")

    def _path(self, library_id):
        if not isinstance(library_id, str) or len(library_id) != 28 or not library_id.startswith("lib-") or any(c not in "0123456789abcdef" for c in library_id[4:]):
            raise ValueError("Invalid external library ID.")
        return library_id

    def get(self, library_id):
        try:
            return self.database.external_library_record(self._path(library_id))
        except KeyError as exc:
            raise ValueError("External library metadata record not found in the Studio Database.") from exc

    def list(self):
        return self.database.external_library_records()

    def register(self, paths, *, name, id_column="id", smiles_column="smiles", name_column=""):
        id_column, smiles_column, name_column = str(id_column or "").strip(), str(smiles_column or "").strip(), str(name_column or "").strip()
        if not str(name).strip() or (smiles_column and id_column == smiles_column):
            raise ValueError("Supply a library name and distinct mapped ID/SMILES fields when using a table.")
        paths = sorted({Path(p).expanduser().resolve(strict=True) for p in paths})
        if not paths:
            raise ValueError("Supply at least one library file.")
        sources = []
        for path in paths:
            if not path.is_file() or path.suffix.lower() not in {".csv", ".tsv", ".parquet", ".feather", ".arrow", ".sdf"}:
                raise ValueError("Supported inputs are headed CSV/TSV, Parquet, Arrow/Feather, and SDF files.")
            columns = source_columns(path)
            is_sdf = path.suffix.lower() == ".sdf"
            if (not is_sdf and (not id_column or id_column not in columns or not smiles_column or smiles_column not in columns)
                    or is_sdf and id_column and id_column not in columns):
                raise ValueError(f"Missing required ID/SMILES fields in {path.name}; available: {columns}")
            if name_column and name_column not in columns:
                raise ValueError(f"Missing optional name field {name_column!r} in {path.name}; available: {columns}")
            stat = path.stat()
            sources.append({"path": str(path), "format": path.suffix.lower().lstrip("."),
                            "size_bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns})
        descriptor = {"schema_version": LIBRARY_SCHEMA, "sources": sources,
                      "id_column": id_column, "name_column": name_column or "", "smiles_column": smiles_column or ""}
        library_id = "lib-" + digest_json(descriptor)[:24]
        try:
            record = self.get(library_id)
            if not record.get("owner_node_id"):
                from mn_molecular_design_studio.cluster.config import identity
                record["owner_node_id"] = identity()
                self.database.save_external_library_record(record)
            return record
        except ValueError:
            from mn_molecular_design_studio.cluster.config import identity
            record = dict(descriptor, id=library_id, name=str(name).strip(), owner_node_id=identity(), created_at=_now(),
                          source_verification="size_and_mtime_only", preparations=[])
            return self.database.save_external_library_record(record)

    def attach_preparation(self, library_id, manifest_path):
        path = Path(manifest_path).resolve(strict=True)
        manifest = json.loads(path.read_text())
        if manifest.get("library_id") != library_id or manifest.get("status") != "completed":
            raise ValueError("Only a completed preparation for this library can be cataloged.")
        reference = {"schema_version": manifest["schema_version"], "path": str(path),
                     "sha256": digest_file(path), "rows": manifest.get("row_count"),
                     "input_rows": manifest.get("input_rows"), "rejected_rows": manifest.get("rejected_rows"),
                     "fingerprint_spec": manifest.get("fingerprint_spec"), "created_at": manifest["created_at"]}
        record = self.get(library_id)
        if not record.get("owner_node_id"):
            from mn_molecular_design_studio.cluster.config import identity
            record["owner_node_id"] = identity()
        reference["owner_node_id"] = record["owner_node_id"]
        if reference not in record["preparations"]:
            record["preparations"].append(reference)
            if manifest.get("schema_version") == FEATURE_SCHEMA:
                record["prepared_compound_count"] = manifest.get("row_count")
                record["rejected_row_count"] = manifest.get("rejected_rows")
                record["fingerprint_spec"] = manifest.get("fingerprint_spec")
            elif manifest.get("schema_version") == CENTROID_SCHEMA:
                record["centroid_count"] = manifest.get("row_count")
            self.database.save_external_library_record(record)
        return reference

    def attach_transferred_preparation(self, library_id, manifest_path, *, source_node_id):
        """Catalog organizer-side copies of completed remote fingerprint outputs."""
        path = Path(manifest_path).resolve(strict=True)
        manifest = json.loads(path.read_text(encoding="utf-8"))
        if (manifest.get("schema_version") not in {FEATURE_SCHEMA, CENTROID_SCHEMA}
                or manifest.get("status") != "completed" or manifest.get("library_id") != library_id):
            raise ValueError("The transferred fingerprint files belong to another library.")
        from mn_molecular_design_studio.cluster.config import identity
        record = self.get(library_id)
        reference = {"schema_version": manifest["schema_version"], "path": str(path),
                     "sha256": digest_file(path), "rows": manifest.get("row_count"),
                     "input_rows": manifest.get("input_rows"), "rejected_rows": manifest.get("rejected_rows"),
                     "fingerprint_spec": manifest.get("fingerprint_spec"), "created_at": manifest["created_at"],
                     "owner_node_id": identity(), "computed_on_node_id": str(source_node_id)}
        if not any(item.get("sha256") == reference["sha256"] and item.get("owner_node_id") == identity()
                   for item in record.setdefault("preparations", [])):
            record["preparations"].append(reference)
        if manifest.get("schema_version") == FEATURE_SCHEMA:
            record["prepared_compound_count"] = manifest.get("row_count")
            record["rejected_row_count"] = manifest.get("rejected_rows")
            record["fingerprint_spec"] = manifest.get("fingerprint_spec")
        else:
            record["centroid_count"] = manifest.get("row_count")
        self.database.save_external_library_record(record)
        return reference

    def import_remote_record(self, record, *, node_id):
        """Save a remote node's metadata and locations without copying its files."""
        if not isinstance(record, dict) or record.get("schema_version") != LIBRARY_SCHEMA:
            raise ValueError("Invalid external library metadata from a compute workstation.")
        library_id = self._path(record.get("id"))
        sources = record.get("sources")
        preparations = record.get("preparations", [])
        if not isinstance(sources, list) or not isinstance(preparations, list):
            raise ValueError("External library locations are incomplete.")
        if (not str(record.get("name") or "").strip()
                or any(not isinstance(item, dict) or not isinstance(item.get("path"), str) for item in sources)
                or any(not isinstance(item, dict) or not isinstance(item.get("path"), str) for item in preparations)):
            raise ValueError("External library metadata contains invalid file locations.")
        if len(json.dumps(record, ensure_ascii=False, allow_nan=False, default=str)) > 1_000_000:
            raise ValueError("External library metadata exceeds the synchronization limit.")
        owner_node_id = str(uuid.UUID(str(node_id)))
        imported = dict(record, owner_node_id=owner_node_id)
        imported["sources"] = [dict(item, owner_node_id=owner_node_id) for item in sources]
        imported["preparations"] = [dict(item, owner_node_id=owner_node_id) for item in preparations]
        try:
            existing = self.get(library_id)
        except ValueError:
            existing = None
        if existing:
            from mn_molecular_design_studio.cluster.config import identity
            local_preparations = [item for item in existing.get("preparations", [])
                                  if item.get("owner_node_id") == identity()]
            known_hashes = {item.get("sha256") for item in imported["preparations"]}
            imported["preparations"].extend(item for item in local_preparations
                                             if item.get("sha256") not in known_hashes)
            if local_preparations:
                latest = local_preparations[-1]
                if latest.get("schema_version") == FEATURE_SCHEMA:
                    imported["prepared_compound_count"] = latest.get("rows")
                    imported["rejected_row_count"] = latest.get("rejected_rows")
                    imported["fingerprint_spec"] = latest.get("fingerprint_spec")
                elif latest.get("schema_version") == CENTROID_SCHEMA:
                    imported["centroid_count"] = latest.get("rows")
        # Coordinator imports are local metadata records; peer imports carry an
        # authority marker so that only the owning workstation republishes them.
        return self.database.save_external_library_record(imported)


def source_columns(path):
    path = Path(path)
    if path.suffix.lower() == ".sdf":
        from rdkit import Chem
        supplier = Chem.SDMolSupplier(str(path), sanitize=False, removeHs=False, strictParsing=False)
        for mol in supplier:
            if mol is not None:
                return sorted({"_Name", *mol.GetPropNames(includePrivate=True, includeComputed=False)})
        return ["_Name"]
    if path.suffix.lower() in {".csv", ".tsv"}:
        with path.open(newline="", encoding="utf-8-sig") as handle:
            return next(csv.reader(handle, delimiter="\t" if path.suffix.lower() == ".tsv" else ","), [])
    import pyarrow as pa
    import pyarrow.ipc as ipc
    import pyarrow.parquet as pq
    if path.suffix.lower() == ".parquet":
        return pq.ParquetFile(path).schema_arrow.names
    with pa.memory_map(str(path), "r") as handle:
        return ipc.open_file(handle).schema.names


def check_sources(record):
    for source in record["sources"]:
        stat = Path(source["path"]).stat()
        if stat.st_size != source["size_bytes"] or stat.st_mtime_ns != source["mtime_ns"]:
            raise ValueError("Library source changed; register a new edition before preparing it.")


def iter_source_records(record, batch_size=10_000):
    """Read mapped ID, optional name, and SMILES/structure; CSV IDs retain zeroes."""
    from rdkit import Chem
    id_col, name_col, smiles_col = record["id_column"], record.get("name_column") or "", record.get("smiles_column") or ""
    for source_index, source in enumerate(record["sources"], 1):
        path = Path(source["path"])
        if path.suffix.lower() == ".sdf":
            supplier = Chem.SDMolSupplier(str(path), sanitize=True, removeHs=False, strictParsing=False)
            for record_index, mol in enumerate(supplier, 1):
                if mol is None:
                    yield None, None, None, "invalid_sdf_record"
                    continue
                raw_id = mol.GetProp(id_col).strip() if id_col and mol.HasProp(id_col) else ""
                if not raw_id:
                    raw_id = f"generated-{source_index:04d}-{record_index:012d}"
                raw_name = mol.GetProp(name_col) if name_col and mol.HasProp(name_col) else None
                if not raw_name and mol.HasProp("_Name"):
                    raw_name = mol.GetProp("_Name")
                smiles = Chem.MolToSmiles(Chem.RemoveHs(mol), canonical=False, isomericSmiles=True)
                yield raw_id, raw_name, smiles, None
            continue
        if path.suffix.lower() in {".csv", ".tsv"}:
            with path.open(newline="", encoding="utf-8-sig") as handle:
                for row in csv.DictReader(handle, delimiter="\t" if path.suffix.lower() == ".tsv" else ","):
                    yield row.get(id_col), row.get(name_col) if name_col else None, row.get(smiles_col), None
        else:
            import pyarrow as pa
            import pyarrow.ipc as ipc
            import pyarrow.parquet as pq
            requested_columns = [id_col, smiles_col] + ([name_col] if name_col else [])
            if path.suffix.lower() == ".parquet":
                for batch in pq.ParquetFile(path).iter_batches(batch_size=batch_size, columns=requested_columns, use_threads=False):
                    columns = batch.to_pydict()
                    for index, (raw_id, raw_smiles) in enumerate(zip(columns[id_col], columns[smiles_col])):
                        yield raw_id, columns[name_col][index] if name_col else None, raw_smiles, None
            else:
                with pa.memory_map(str(path), "r") as handle:
                    reader = ipc.open_file(handle)
                    for index in range(reader.num_record_batches):
                        batch = reader.get_batch(index)
                        for start in range(0, batch.num_rows, batch_size):
                            columns = batch.slice(start, batch_size).select(requested_columns).to_pydict()
                            for row_index, (raw_id, raw_smiles) in enumerate(zip(columns[id_col], columns[smiles_col])):
                                yield raw_id, columns[name_col][row_index] if name_col else None, raw_smiles, None


@dataclass(frozen=True)
class FingerprintSpec:
    bits: int = 1024
    radius: int = 2
    use_chirality: bool = True
    # Match OpenVS's library-preparation scripts: feature atom invariants and
    # chirality are both enabled for 1024-bit, radius-2 Morgan fingerprints.
    use_features: bool = True

    def __post_init__(self):
        if type(self.bits) is not int or self.bits <= 0 or self.bits % 8 or type(self.radius) is not int or self.radius < 0:
            raise ValueError("Fingerprint size must be a positive multiple of eight; radius must be nonnegative.")
        if type(self.use_chirality) is not bool or type(self.use_features) is not bool:
            raise ValueError("Fingerprint flags must be boolean.")

    def metadata(self):
        import rdkit
        return dict(asdict(self), algorithm="morgan", encoding="rdkit_binary_text",
                    hydrogen_policy="sanitize_remove_explicit_hydrogens",
                    standardization="retain_charge_fragments_tautomers_and_stereochemistry",
                    rdkit_version=rdkit.__version__)


def new_output_directory(path):
    path = Path(path).expanduser().resolve()
    if path.is_relative_to(app_home()):
        raise ValueError("Choose an external output directory outside the Studio workspace.")
    path.mkdir(parents=True, exist_ok=False)
    return path


def feature_arrow_schema(spec):
    import pyarrow as pa
    return pa.schema([("compound_key", pa.uint64()), ("compound_id", pa.string()),
                      ("compound_name", pa.string()), ("source_smiles", pa.string()),
                      ("fingerprint", pa.binary(spec.bits // 8))])


def _fingerprint_batch(batch, spec):
    """CPU-bound, OpenVS-compatible RDKit work for one batch; writer stays in parent."""
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem

    output = []
    for row_number, compound_id, raw_name, smiles in batch:
        mol = Chem.MolFromSmiles(smiles) if smiles else None
        if mol is None or mol.GetNumHeavyAtoms() == 0:
            output.append({"row_number": row_number, "compound_id": compound_id,
                           "compound_name": raw_name, "source_smiles": smiles, "reason": "invalid_or_empty_smiles"})
            continue
        mol = Chem.RemoveHs(mol)
        output.append({"row_number": row_number, "compound_id": compound_id,
                       "compound_name": str(raw_name).strip() if raw_name is not None and str(raw_name).strip() else None,
                       "source_smiles": smiles,
                       "fingerprint": DataStructs.BitVectToBinaryText(AllChem.GetMorganFingerprintAsBitVect(
                           mol, radius=spec.radius, nBits=spec.bits, useFeatures=spec.use_features,
                           useChirality=spec.use_chirality))})
    return output


def prepare_fingerprints(record, output, *, spec=None, batch_size=10_000, shard_rows=1_000_000,
                         progress_callback=None, worker_count=1):
    """Stream source records through a bounded CPU process pool into Parquet shards.

    Invalid SMILES are recorded and excluded; missing or duplicate table IDs fail
    the preparation. SDF rows without IDs get stable file/record-based IDs. An
    external disk index enforces ID uniqueness across all shards.
    Partial/failed outputs never have a manifest.json and require a new output folder.
    """
    import pyarrow as pa
    import pyarrow.parquet as pq

    spec = spec or FingerprintSpec()
    if type(batch_size) is not int or type(shard_rows) is not int or batch_size < 1 or shard_rows < batch_size:
        raise ValueError("Require positive batch_size and shard_rows >= batch_size.")
    if type(worker_count) is not int or not 1 <= worker_count <= 256:
        raise ValueError("CPU worker count must be an integer from 1 to 256.")
    from mn_molecular_design_studio.engine.workflow_engine.vendor.ligand_runtime.mn_ligand.runtime import cpu_process_limit
    cpu_limit = cpu_process_limit()
    if worker_count > cpu_limit:
        raise ValueError(f"CPU worker count exceeds this workstation's configured limit of {cpu_limit}.")
    check_sources(record)
    output = new_output_directory(output)
    schema = feature_arrow_schema(spec)
    shards, rows, rejects = [], [], []
    writer = None
    shard_count = count = rejected = read_count = 0
    rejection_schema = pa.schema([("row_number", pa.uint64()), ("compound_id", pa.string()),
                                 ("smiles", pa.string()), ("reason", pa.string())])
    reject_writer = pq.ParquetWriter(output / "rejected.parquet", rejection_schema, compression="zstd")
    index = sqlite3.connect(output / ".id-index.sqlite3")
    index.execute("CREATE TABLE ids (id TEXT PRIMARY KEY)")
    executor = ProcessPoolExecutor(max_workers=worker_count, mp_context=get_context("spawn")) if worker_count > 1 else None
    pending = deque()
    input_batch = []

    def flush():
        nonlocal writer, shard_count
        if not rows:
            return
        if writer is None:
            writer = pq.ParquetWriter(output / f"features-{len(shards):06d}.parquet", schema, compression="zstd")
        writer.write_table(pa.Table.from_pylist(rows, schema=schema))
        shard_count += len(rows)
        rows.clear()
        if shard_count >= shard_rows:
            close_shard()

    def close_shard():
        nonlocal writer, shard_count
        if writer is not None:
            writer.close()
            path = output / f"features-{len(shards):06d}.parquet"
            shards.append({"path": path.name, "rows": shard_count, "sha256": digest_file(path), "size_bytes": path.stat().st_size})
            writer, shard_count = None, 0

    def consume_fingerprint_results(results):
        nonlocal count, rejected
        for item in results:
            if item.get("reason"):
                rejected += 1
                rejects.append({"row_number": item["row_number"], "compound_id": item["compound_id"],
                                "smiles": item["source_smiles"], "reason": item["reason"]})
                if len(rejects) >= batch_size:
                    reject_writer.write_table(pa.Table.from_pylist(rejects, schema=rejection_schema)); rejects.clear()
                continue
            rows.append({"compound_key": count, "compound_id": item["compound_id"],
                         "compound_name": item["compound_name"], "source_smiles": item["source_smiles"],
                         "fingerprint": item["fingerprint"]})
            count += 1
            if len(rows) >= batch_size:
                flush()

    def dispatch_batch(*, drain=False):
        if input_batch:
            batch = tuple(input_batch)
            input_batch.clear()
            if executor:
                pending.append(executor.submit(_fingerprint_batch, batch, spec))
            else:
                consume_fingerprint_results(_fingerprint_batch(batch, spec))
        if pending and (drain or len(pending) >= worker_count * 2):
            consume_fingerprint_results(pending.popleft().result())

    try:
        for read_count, (raw_id, raw_name, raw_smiles, source_error) in enumerate(iter_source_records(record, batch_size), 1):
            if read_count % batch_size == 0:
                index.commit()
                if progress_callback:
                    progress_callback(read_count, rejected)
            if source_error:
                rejected += 1
                rejects.append({"row_number": read_count, "compound_id": None, "smiles": None, "reason": source_error})
                if len(rejects) >= batch_size:
                    reject_writer.write_table(pa.Table.from_pylist(rejects, schema=rejection_schema)); rejects.clear()
                continue
            if raw_id is None or not str(raw_id).strip():
                raise ValueError(f"Missing compound ID at source row {read_count}.")
            compound_id = str(raw_id).strip()
            try:
                index.execute("INSERT INTO ids VALUES (?)", (compound_id,))
            except sqlite3.IntegrityError as exc:
                raise ValueError(f"Duplicate compound ID at source row {read_count}: {compound_id}") from exc
            smiles = "" if raw_smiles is None else str(raw_smiles).strip()
            raw_name = None if raw_name is None else str(raw_name)
            input_batch.append((read_count, compound_id, raw_name, smiles))
            if len(input_batch) >= batch_size:
                dispatch_batch()
        dispatch_batch(drain=True)
        while pending:
            consume_fingerprint_results(pending.popleft().result())
        flush(); close_shard()
        if rejects:
            reject_writer.write_table(pa.Table.from_pylist(rejects, schema=rejection_schema))
        reject_writer.close()
        check_sources(record)
        if not count:
            raise ValueError("The library has no valid compounds.")
        if progress_callback:
            progress_callback(read_count, rejected)
        manifest = {"schema_version": FEATURE_SCHEMA, "status": "completed", "library_id": record["id"],
                    "created_at": _now(), "source_descriptor_sha256": digest_json({
                        "sources": record["sources"], "id_column": record["id_column"],
                        "name_column": record.get("name_column", ""), "smiles_column": record.get("smiles_column", "")}),
                    "source_verification": "size_and_mtime_only", "fingerprint_spec": spec.metadata(),
                    "input_rows": read_count, "row_count": count, "rejected_rows": rejected,
                    "cpu_parallelism": "process_pool" if worker_count > 1 else "serial",
                    "cpu_worker_count": worker_count,
                    "rejections": {"path": "rejected.parquet", "sha256": digest_file(output / "rejected.parquet")},
                    "compression": "zstd", "table_format": "parquet",
                    "row_fields": ["compound_id", "compound_name", "source_smiles", "fingerprint"],
                    "shards": shards}
        atomic_json(output / "manifest.json", manifest)
        return output / "manifest.json"
    except Exception as exc:
        atomic_json(output / "failure.json", {"status": "failed", "error": str(exc), "input_rows_seen": read_count})
        raise
    finally:
        if executor:
            executor.shutdown(wait=True, cancel_futures=True)
        if writer is not None:
            writer.close()
        reject_writer.close()
        index.close()
        (output / ".id-index.sqlite3").unlink(missing_ok=True)


def read_feature_manifest(path):
    path = Path(path).resolve(strict=True)
    manifest = json.loads(path.read_text())
    if manifest.get("schema_version") != FEATURE_SCHEMA or manifest.get("status") != "completed":
        raise ValueError("Use a completed library fingerprint manifest.")
    return path, manifest


def tanimoto_similarity(packed_a: bytes, packed_b: bytes) -> float:
    """Calculate RDKit Tanimoto for fingerprints stored by this module."""
    from rdkit import DataStructs

    packed_a, packed_b = bytes(packed_a), bytes(packed_b)
    if not packed_a or len(packed_a) != len(packed_b):
        raise ValueError("Packed fingerprints must have the same nonzero byte length.")
    return float(DataStructs.TanimotoSimilarity(
        DataStructs.CreateFromBinaryText(packed_a), DataStructs.CreateFromBinaryText(packed_b)))


def iter_feature_batches(path, *, batch_size=10_000):
    import pyarrow.parquet as pq
    path, manifest = read_feature_manifest(path)
    total = 0
    for shard in manifest["shards"]:
        source = (path.parent / shard["path"]).resolve(strict=True)
        if not source.is_relative_to(path.parent) or digest_file(source) != shard["sha256"]:
            raise ValueError("Fingerprint shard escapes its manifest directory or has changed.")
        parquet = pq.ParquetFile(source)
        if parquet.metadata.num_rows != shard["rows"]:
            raise ValueError("Fingerprint shard row count disagrees with its manifest.")
        for batch in parquet.iter_batches(batch_size=batch_size, use_threads=False):
            total += batch.num_rows
            yield batch
    if total != manifest["row_count"]:
        raise ValueError("Fingerprint manifest row count disagrees with its shards.")
