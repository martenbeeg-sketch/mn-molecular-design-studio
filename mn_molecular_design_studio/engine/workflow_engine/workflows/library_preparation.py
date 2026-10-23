"""Studio's scheduled CPU workflow for external library preparation."""
from __future__ import annotations

import json
from pathlib import Path

from mn_molecular_design_studio.collection_policy import collection_scope, set_run_policy
from mn_molecular_design_studio.external_libraries import LibraryCatalog, digest_file
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import create_job, finish_job, read_json, update_status, write_json


def validate_library_request(parameters, *, require_files=True):
    if not str(parameters.get("library_name", "")).strip():
        raise ValueError("Enter a library name.")
    files = parameters.get("source_files") or []
    if not isinstance(files, list) or not files:
        raise ValueError("Choose at least one ID/SMILES table or SDF file.")
    id_column = str(parameters.get("id_column", "id")).strip()
    smiles_column = str(parameters.get("smiles_column", "smiles")).strip()
    name_column = str(parameters.get("name_column", "")).strip()
    all_sdf = all(Path(path).suffix.lower() == ".sdf" for path in files)
    if (not id_column and not all_sdf) or (smiles_column and id_column == smiles_column) or (name_column and name_column in {id_column, smiles_column}):
        raise ValueError("Choose distinct ID, Name, and SMILES fields; SDF-only input can generate missing IDs.")
    if parameters.get("create_centroids", False):
        raise ValueError("Centroid generation is a separate tile. Queue fingerprint preparation only here.")
    # Fingerprinting uses the distributed worker's own process pool. The live
    # scheduler allocation is checked per selected computer below; it is not
    # constrained by the general ligand-runtime process setting.
    cpu_limit = 256
    if type(parameters.get("cpu_cores", 1)) is not int or not 1 <= parameters.get("cpu_cores", 1) <= cpu_limit:
        raise ValueError(f"CPU worker count must be between 1 and {cpu_limit}.")
    max_records = parameters.get("max_records", 0)
    if type(max_records) is not int or max_records < 0:
        raise ValueError("Rows to process must be zero (all rows) or a positive integer.")
    selection_mode = str(parameters.get("selection_mode", "first_rows"))
    if selection_mode not in {"first_rows", "random_sample"}:
        raise ValueError("Choose first rows or a random sample across all source rows.")
    sample_seed = parameters.get("sample_seed", 42)
    if type(sample_seed) is not int or not 0 <= sample_seed <= 2**31 - 1:
        raise ValueError("Random sample seed must be a whole number from 0 to 2,147,483,647.")
    if selection_mode == "random_sample" and max_records < 1:
        raise ValueError("A random sample requires a positive sample size.")
    ssh_user = str(parameters.get("ssh_user", "compute")).strip()
    if not ssh_user or any(not (char.isalnum() or char in "_-") for char in ssh_user) or ssh_user[0].isdigit():
        raise ValueError("Enter a valid SSH account name for compute workstations.")
    worker_nodes = parameters.get("worker_nodes") or []
    if not isinstance(worker_nodes, list):
        raise ValueError("Select one or more worker computers.")
    if not worker_nodes:
        from mn_molecular_design_studio.cluster import config as cluster_config
        worker_nodes = [{"node_id": cluster_config.identity(), "cpu_cores": parameters.get("cpu_cores", 1)}]
        parameters["worker_nodes"] = worker_nodes
    if require_files:
        from mn_molecular_design_studio.cluster.distributed_library import _node_specs
        try:
            selected_nodes = _node_specs(worker_nodes)
        except (KeyError, TypeError, ValueError, OSError) as exc:
            raise ValueError(str(exc)) from exc
        local_workers = next(row["cpu_cores"] for row in selected_nodes if row["local"])
        if local_workers != parameters.get("cpu_cores", 1):
            raise ValueError("CPU workers on the organizer must match the scheduled CPU allocation.")
    from mn_molecular_design_studio.runtime import app_home
    output = Path(parameters.get("output_directory", "")).expanduser()
    if (not str(parameters.get("output_directory", "")).strip() or not output.is_absolute()
            or require_files and output.resolve().is_relative_to(app_home())):
        raise ValueError("Choose an absolute external project directory outside the Studio workspace.")
    if require_files and output.exists():
        fingerprint_dir = output / "fingerprints"
        if not output.is_dir() or fingerprint_dir.exists() and (not fingerprint_dir.is_dir() or any(fingerprint_dir.iterdir())):
            raise ValueError("The project folder may exist, but its fingerprints folder must be new or empty.")
    if require_files and max_records:
        sample_name = (f"random-sample-{max_records}-seed-{sample_seed}.parquet"
                       if selection_mode == "random_sample" else f"first-{max_records}-rows.parquet")
        sample = output / "input" / sample_name
        if sample.exists():
            raise ValueError(f"The sample file already exists: {sample}. Choose a new output folder.")
    if require_files:
        from mn_molecular_design_studio.external_libraries import source_columns
        for path in files:
            path = Path(path)
            if not path.is_file() or path.suffix.lower() not in {".csv", ".tsv", ".parquet", ".feather", ".arrow", ".sdf"}:
                raise ValueError("Choose CSV/TSV, Parquet, Arrow/Feather, or SDF library files.")
            columns = source_columns(path)
            required = {id_column} if id_column else set()
            if path.suffix.lower() != ".sdf":
                if not id_column:
                    raise ValueError(f"Set the ID column for table file {path.name}.")
                if not smiles_column:
                    raise ValueError(f"Set the SMILES column for table file {path.name}.")
                required.add(smiles_column)
            if name_column:
                required.add(name_column)
            if not required.issubset(columns):
                raise ValueError(f"Mapped fields are missing in {path.name}; available: {columns}.")


def run_library_preparation(source_files: list[Path], library_name: str, output_directory: str = "",
                            id_column: str = "id", smiles_column: str = "smiles", name_column: str = "",
                            create_centroids: bool = False, cpu_cores: int = 1,
                            worker_nodes: list[dict] | None = None, max_records: int = 0,
                            selection_mode: str = "first_rows", sample_seed: int = 42,
                            ssh_user: str = "compute") -> Path:
    parameters = dict(source_files=source_files, library_name=library_name, output_directory=output_directory,
                      id_column=id_column, smiles_column=smiles_column, name_column=name_column,
                      create_centroids=create_centroids, cpu_cores=cpu_cores,
                      worker_nodes=worker_nodes or [], max_records=max_records,
                      selection_mode=selection_mode, sample_seed=sample_seed, ssh_user=ssh_user)
    validate_library_request(parameters)
    with collection_scope("manual"):
        job = create_job("library-preparation", "library_preparation", "library_preparation",
                         {"source_files": [str(p) for p in source_files]},
                         {k: [str(p) for p in v] if k == "source_files" else str(v) if isinstance(v, Path) else v for k, v in parameters.items()})
    metadata = read_json(job.run_dir / "metadata.json")
    metadata["campaign_name"] = library_name.strip()
    write_json(job.run_dir / "metadata.json", metadata)
    set_run_policy(job.run_dir, "manual")
    catalog = LibraryCatalog()
    record = None
    try:
        update_status(job.run_dir, "running", current_phase="Registering external library")
        record = catalog.register(source_files, name=library_name, id_column=id_column,
                                  smiles_column=smiles_column, name_column=name_column)
        outputs = {"library_id": record["id"], "source_files": [str(p) for p in source_files], "database_collection_policy": "manual"}
        metrics = {"source_file_count": len(source_files)}
        base = Path(output_directory).expanduser().resolve() if output_directory else None

        def progress(count, extra):
            if read_json(job.run_dir / "metadata.json").get("status") in {"cancelled", "stopped", "paused"}:
                raise RuntimeError("Library preparation was cancelled or paused.")
            phase = extra.get("phase") if isinstance(extra, dict) else None
            if phase == "random_sampling":
                label = f"Selecting random sample: scanned {count:,} source rows"
            elif phase == "fingerprinting":
                label = f"Fingerprinting sample: {count:,} rows processed"
            else:
                label = f"Library preparation: {count:,} records processed"
            update_status(job.run_dir, "running", current_phase=label)

        base.mkdir(parents=True, exist_ok=True)
        fingerprint_dir = base / "fingerprints"
        fingerprint_dir.mkdir(parents=True, exist_ok=True)
        update_status(job.run_dir, "running", current_phase="Preparing packed Morgan fingerprints")
        from mn_molecular_design_studio.cluster import config as cluster_config
        from mn_molecular_design_studio.cluster.distributed_library import prepare_distributed_fingerprints
        selected_workers = worker_nodes or [{"node_id": cluster_config.identity(), "cpu_cores": cpu_cores}]
        from mn_molecular_design_studio.cluster.distributed_library import _node_specs
        worker_plan = [{"node_id": item["node_id"], "name": item["name"],
                        "cpu_worker_count": item["cpu_cores"], "status": "starting"}
                       for item in _node_specs(selected_workers)]
        update_status(job.run_dir, "running", current_phase="Starting distributed fingerprint workers",
                      worker_plan=worker_plan)
        sample_name = (f"random-sample-{max_records}-seed-{sample_seed}.parquet"
                       if selection_mode == "random_sample" else f"first-{max_records}-rows.parquet")
        sample_parquet = base / "input" / sample_name if max_records else None
        if sample_parquet:
            (base / "input").mkdir(parents=True, exist_ok=True)
        feature_manifest, _, _ = prepare_distributed_fingerprints(
            record, fingerprint_dir, selected_workers, max_records=max_records,
            progress_callback=progress, sample_parquet=sample_parquet, ssh_user=ssh_user,
            selection_mode=selection_mode, sample_seed=sample_seed,
        )
        catalog.attach_preparation(record["id"], feature_manifest)
        from mn_molecular_design_studio.external_libraries import read_feature_manifest
        _, features = read_feature_manifest(feature_manifest)
        outputs["feature_manifest"] = str(feature_manifest)
        metrics.update(compound_count=features["row_count"], rejected_count=features["rejected_rows"],
                       input_rows=features["input_rows"], worker_count=features.get("cpu_worker_count"),
                       worker_computers=len(features.get("workers", [])),
                       workers=features.get("workers", []), selection=features.get("selection", {}))
        finish_job(job.run_dir, True, {"outputs": outputs, "metrics": metrics})
    except Exception as exc:
        failure_outputs = {"library_name": library_name}
        if record:
            failure_outputs["library_id"] = record["id"]
        finish_job(job.run_dir, False, {"error": str(exc), "outputs": failure_outputs})
        raise
    return job.run_dir


def validate_library_centroid_request(parameters, *, require_files=True):
    import math

    library_id = str(parameters.get("library_id") or "")
    feature_manifest = Path(str(parameters.get("feature_manifest") or "")).expanduser()
    output_directory = Path(str(parameters.get("output_directory") or "")).expanduser()
    if not library_id:
        raise ValueError("Select a registered external library.")
    if not feature_manifest.is_absolute():
        raise ValueError("Use an absolute prepared fingerprint manifest path on the selected workstation.")
    if not output_directory.is_absolute():
        raise ValueError("Choose a new absolute centroid output directory on the selected workstation.")
    from mn_molecular_design_studio.runtime import app_home
    if require_files and output_directory.resolve().is_relative_to(app_home()):
        raise ValueError("Choose an external centroid output directory outside the Studio workspace.")
    if require_files and output_directory.exists():
        raise ValueError("Choose a new centroid output directory; existing preparations are never overwritten.")
    threshold = parameters.get("similarity_threshold", 0.6)
    if isinstance(threshold, bool) or not isinstance(threshold, (int, float)) or not math.isfinite(threshold) or not 0 < threshold <= 1:
        raise ValueError("Tanimoto similarity threshold must be in (0, 1].")
    method = parameters.get("clustering_method", "sampled")
    if method not in {"sampled", "hierarchical", "exact"}:
        raise ValueError("Choose sampled, hierarchical, or exact centroid clustering.")
    for name, default, minimum in (("cpu_cores", 1, 1), ("chunk_size", 300_000, 1),
                                    ("chunk_count", 0, 0), ("final_candidate_count", 300_000, 1),
                                    ("max_merge_rounds", 3, 0), ("random_seed", 42, 0)):
        value = parameters.get(name, default)
        if type(value) is not int or value < minimum:
            raise ValueError(f"{name.replace('_', ' ')} must be a whole number of {minimum} or greater.")
    if parameters.get("cpu_cores", 1) > 256 or parameters.get("random_seed", 42) >= 2**63:
        raise ValueError("Use at most 256 CPU slots per computer and a random seed below 2^63.")
    ssh_user = str(parameters.get("ssh_user", "compute"))
    from mn_molecular_design_studio.cluster.distributed_library import _SSH_USER
    if not _SSH_USER.fullmatch(ssh_user):
        raise ValueError("Enter a valid SSH account for the selected computers.")
    if require_files:
        from mn_molecular_design_studio.external_libraries import LibraryCatalog, digest_file, read_feature_manifest
        record = LibraryCatalog().get(library_id)
        _, manifest = read_feature_manifest(feature_manifest)
        if manifest.get("library_id") != library_id:
            raise ValueError("The prepared fingerprints belong to another library edition.")
        registered = any(
            item.get("schema_version") == "studio.library-fingerprints.v2"
            and Path(item.get("path", "")).resolve() == feature_manifest.resolve()
            and item.get("sha256") == digest_file(feature_manifest)
            for item in record.get("preparations", [])
        )
        if not registered:
            raise ValueError("Register the completed fingerprint manifest in this library's Database record first.")
        if method != "exact":
            from mn_molecular_design_studio.cluster import config
            from mn_molecular_design_studio.cluster.distributed_library import _node_specs
            from mn_molecular_design_studio.centroid_chunks import balanced_chunk_sizes
            balanced_chunk_sizes(manifest["row_count"], parameters.get("chunk_size", 300_000),
                                 parameters.get("chunk_count", 0))
            nodes = parameters.get("worker_nodes") or [{"node_id": config.identity(),
                                                       "cpu_cores": parameters.get("cpu_cores", 1)}]
            if not isinstance(nodes, list):
                raise ValueError("Choose worker computers for chunk clustering.")
            specs = _node_specs(nodes)
            if next(s["cpu_cores"] for s in specs if s["local"]) != parameters.get("cpu_cores", 1):
                raise ValueError("CPU workers on the source workstation must match its scheduled CPU allocation.")


def run_library_centroid_preparation(library_id: str, feature_manifest: Path, output_directory: str,
                                     similarity_threshold: float = 0.6, clustering_method: str = "sampled",
                                     cpu_cores: int = 1, chunk_size: int = 300_000, chunk_count: int = 0,
                                     final_candidate_count: int = 300_000, random_seed: int = 42,
                                     max_merge_rounds: int = 3, worker_nodes: list[dict] | None = None,
                                     ssh_user: str = "compute") -> Path:
    parameters = dict(library_id=library_id, feature_manifest=str(feature_manifest),
                      output_directory=output_directory, similarity_threshold=similarity_threshold,
                      clustering_method=clustering_method, cpu_cores=cpu_cores, chunk_size=chunk_size,
                      chunk_count=chunk_count, final_candidate_count=final_candidate_count,
                      random_seed=random_seed, max_merge_rounds=max_merge_rounds,
                      worker_nodes=worker_nodes or [], ssh_user=ssh_user)
    validate_library_centroid_request(parameters)
    feature_manifest = Path(feature_manifest).expanduser().resolve(strict=True)
    base = Path(output_directory).expanduser().resolve()
    with collection_scope("manual"):
        job = create_job("library-preparation", "library_centroids", "library_centroids",
                         {"feature_manifest": str(feature_manifest)}, parameters)
    set_run_policy(job.run_dir, "manual")
    try:
        update_status(job.run_dir, "running", current_phase="Creating Tanimoto centroids")

        def progress(rows, representatives, phase=None):
            if read_json(job.run_dir / "metadata.json").get("status") in {"cancelled", "stopped", "paused"}:
                raise RuntimeError("Centroid generation was cancelled or paused.")
            update_status(job.run_dir, "running", current_phase=phase or f"Centroid generation: {rows:,} rows · {representatives:,} representatives")

        base.mkdir(parents=True, exist_ok=False)
        if clustering_method == "exact":
            from mn_molecular_design_studio.library_centroids import prepare_centroids
            centroid_manifest = prepare_centroids(
                feature_manifest, base / "centroids", similarity_threshold=similarity_threshold,
                progress_callback=progress,
            )
        else:
            from mn_molecular_design_studio.cluster import config
            from mn_molecular_design_studio.cluster.distributed_library import _node_specs
            from mn_molecular_design_studio.centroid_chunks import prepare_chunked_centroids
            selected_workers = worker_nodes or [{"node_id": config.identity(), "cpu_cores": cpu_cores}]
            update_status(job.run_dir, "running", worker_plan=[
                {"node_id": item["node_id"], "name": item["name"], "cpu_worker_count": item["cpu_cores"],
                 "status": "starting", "task_kind": "centroids"} for item in _node_specs(selected_workers)])
            centroid_manifest = prepare_chunked_centroids(
                feature_manifest, base / "centroids", similarity_threshold=similarity_threshold,
                worker_nodes=selected_workers, chunk_size=chunk_size, chunk_count=chunk_count,
                max_merge_rounds=max_merge_rounds, random_seed=random_seed, ssh_user=ssh_user,
                method=clustering_method, final_candidate_count=final_candidate_count,
                progress_callback=progress,
            )
        LibraryCatalog().attach_preparation(library_id, centroid_manifest)
        manifest = json.loads(centroid_manifest.read_text(encoding="utf-8"))
        finish_job(job.run_dir, True, {"outputs": {"library_id": library_id,
                                                     "centroid_manifest": str(centroid_manifest),
                                                     "feature_manifest": str(feature_manifest)},
                                       "metrics": {"input_compound_count": manifest["input_rows"],
                                                   "centroid_count": manifest["row_count"],
                                                   "exact_similarity_evaluations": manifest["exact_similarity_evaluations"],
                                                   "algorithm": manifest["algorithm"],
                                                   "candidate_count": manifest.get("candidate_count"),
                                                   "workers": manifest.get("workers", []),
                                                   "rounds": len(manifest.get("rounds", [])),
                                                   "elapsed_seconds": manifest.get("elapsed_seconds")}})
    except Exception as exc:
        finish_job(job.run_dir, False, {"error": str(exc), "outputs": {"library_id": library_id}})
        raise
    return job.run_dir


def _centroid_seed_inputs(library_id: str, centroid_manifest: Path, *, sample_size: int, sampling_seed: int):
    import json
    from mn_molecular_design_studio.external_libraries import LibraryCatalog

    centroid_manifest = Path(centroid_manifest).expanduser().resolve(strict=True)
    manifest = json.loads(centroid_manifest.read_text(encoding="utf-8"))
    if manifest.get("schema_version") != "studio.library-centroids.v2" or manifest.get("status") != "completed":
        raise ValueError("Choose a completed centroid manifest.")
    if manifest.get("library_id") != library_id:
        raise ValueError("The centroid set belongs to another registered library edition.")
    record = LibraryCatalog().get(library_id)
    catalog_reference = next((p for p in record["preparations"]
                              if p.get("schema_version") == "studio.library-centroids.v2"
                              and Path(p["path"]).resolve() == centroid_manifest), None)
    if not catalog_reference or catalog_reference.get("sha256") != digest_file(centroid_manifest):
        raise ValueError("Register this completed centroid manifest in the library catalog before docking it.")
    centroid_file_info = next((item for item in manifest.get("files", []) if item.get("path") == "centroids.parquet"), None)
    if not centroid_file_info:
        raise ValueError("The centroid manifest does not reference its centroid table.")
    centroid_file = (centroid_manifest.parent / centroid_file_info["path"]).resolve(strict=True)
    if not centroid_file.is_relative_to(centroid_manifest.parent) or digest_file(centroid_file) != centroid_file_info.get("sha256"):
        raise ValueError("The centroid table is outside its manifest directory or has changed.")
    if type(sample_size) is not int or sample_size < 1 or sample_size > manifest.get("row_count", 0):
        raise ValueError(f"Choose between 1 and {manifest.get('row_count', 0):,} available centroid compounds.")
    import numpy as np
    import pyarrow.parquet as pq
    parquet = pq.ParquetFile(centroid_file)
    if parquet.metadata.num_rows != manifest.get("row_count"):
        raise ValueError("The centroid table row count does not match its manifest.")
    rng = np.random.default_rng(sampling_seed)
    selected = []
    seen = 0
    for batch in parquet.iter_batches(batch_size=10_000,
                                      columns=["compound_key", "compound_id", "compound_name", "source_smiles", "fingerprint"]):
        for row in batch.to_pylist():
            if seen < sample_size:
                selected.append(row)
            else:
                index = int(rng.integers(0, seen + 1))
                if index < sample_size:
                    selected[index] = row
            seen += 1
    if seen != manifest.get("row_count"):
        raise ValueError("The centroid table row count changed while sampling.")
    return record, manifest, selected


def validate_library_seed_dock_request(parameters, *, require_files=True):
    import math
    from mn_molecular_design_studio.engine.workflow_engine.workflows.ligand_tools import MAX_COMPOUNDS_PER_LAUNCH

    library_id = parameters.get("library_id")
    centroid_manifest = parameters.get("centroid_manifest")
    if not isinstance(library_id, str) or not library_id:
        raise ValueError("Select a registered external library.")
    if not centroid_manifest:
        raise ValueError("Select a completed centroid set for this library.")
    if type(parameters.get("sample_size")) is not int or not 1 <= parameters["sample_size"] <= MAX_COMPOUNDS_PER_LAUNCH:
        raise ValueError(f"Choose between 1 and {MAX_COMPOUNDS_PER_LAUNCH:,} centroid compounds for one Uni-Dock Pro launch.")
    if type(parameters.get("sampling_seed", 0)) is not int or parameters.get("sampling_seed", 0) < 0:
        raise ValueError("The sampling seed must be a nonnegative integer.")
    if type(parameters.get("seed", 42)) is not int or not 1 <= parameters.get("seed", 42) <= 2_147_483_646:
        raise ValueError("The docking seed must be a positive 32-bit integer.")
    if not str(parameters.get("target_pdb", "")).strip():
        raise ValueError("Select a receptor PDB file.")
    center = [parameters.get(f"center_{axis}") for axis in "xyz"]
    sizes = [parameters.get(f"size_{axis}", 20.0) for axis in "xyz"]
    if any(isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) for value in center):
        raise ValueError("Docking box center coordinates must be finite numbers.")
    if any(isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or not 1 <= value <= 60 for value in sizes):
        raise ValueError("Docking box dimensions must be between 1 and 60 Å.")
    if parameters.get("search_mode", "balance") not in {"fast", "balance", "detail"}:
        raise ValueError("Choose a supported Uni-Dock Pro search mode.")
    if type(parameters.get("num_modes", 9)) is not int or not 1 <= parameters.get("num_modes", 9) <= 100:
        raise ValueError("Pose modes per compound must be between 1 and 100.")
    if type(parameters.get("cpu_cores", 4)) is not int or not 1 <= parameters.get("cpu_cores", 4) <= 256:
        raise ValueError("CPU cores must be between 1 and 256.")
    if not str(parameters.get("gpu_device", "0")).strip() or str(parameters.get("gpu_device", "0")).strip().lower() in {"none", "cpu", "off"}:
        raise ValueError("Uni-Dock Pro requires a CUDA GPU device.")
    if require_files:
        from mn_molecular_design_studio.engine.workflow_engine.workflows.ligand_tools import validate_target_pdb
        validate_target_pdb(Path(parameters["target_pdb"]))
        _centroid_seed_inputs(library_id, Path(centroid_manifest),
                              sample_size=parameters["sample_size"],
                              sampling_seed=parameters.get("sampling_seed", 0))


def run_library_seed_docking(library_id: str, centroid_manifest: Path, target_pdb: Path,
                             center_x: float, center_y: float, center_z: float,
                             size_x: float = 20.0, size_y: float = 20.0, size_z: float = 20.0,
                             sample_size: int = 1000, sampling_seed: int = 0,
                             search_mode: str = "balance", num_modes: int = 9,
                             seed: int = 42, cpu_cores: int = 4, gpu_device: str = "0") -> Path:
    """Uniformly sample real centroid compounds and run the existing Uni-Dock Pro adapter."""
    from datetime import datetime, timezone
    import numpy as np
    from mn_molecular_design_studio import ligand_compounds
    from mn_molecular_design_studio.collection_policy import collection_scope, set_run_policy
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import (
        pending_job_creation_dir, read_json, write_json,
    )
    from mn_molecular_design_studio.engine.workflow_engine.workflows.ligand_tools import run_ligand_unidock_pro

    parameters = locals().copy()
    validate_library_seed_dock_request(parameters, require_files=False)
    from mn_molecular_design_studio.engine.workflow_engine.workflows.ligand_tools import validate_target_pdb
    validate_target_pdb(Path(target_pdb))
    queued_run_dir = pending_job_creation_dir()
    if queued_run_dir is None:
        raise RuntimeError("Seed docking must be queued from the Studio Library preparation screen.")
    run_dir = queued_run_dir.resolve()
    set_run_policy(run_dir, "manual")
    record, manifest, selected = _centroid_seed_inputs(
        library_id, centroid_manifest, sample_size=sample_size, sampling_seed=sampling_seed)
    selected_rows = [{
        "compound_id": row["compound_id"],
        "compound_name": row["compound_name"] or row["compound_id"],
        "smiles": row["source_smiles"],
        "properties": {"source_compound_key": int(row["compound_key"]), "library_id": library_id},
    } for row in selected]
    input_dir = run_dir / "artifacts" / "input"
    input_dir.mkdir(parents=True, exist_ok=True)
    compound_set_path = input_dir / "centroid_seed_compounds.json"
    marker = read_json(run_dir / "molecular_design_studio.json")
    created_at = str(read_json(run_dir / "metadata.json").get("created_at") or datetime.now(timezone.utc).isoformat())
    compound_set = ligand_compounds.write_compound_set(
        compound_set_path, name=f"{record['name']} centroid seed ({sample_size})", compounds=selected_rows,
        owner_id=str(marker.get("studio_owner_id") or ""), project_id=str(marker.get("project_id") or ""),
        created_at=created_at,
    )
    selection = {
        "schema_version": "studio.library-centroid-selection.v1", "library_id": library_id,
        "library_name": record["name"], "centroid_manifest": str(Path(centroid_manifest).resolve()),
        "centroid_manifest_sha256": digest_file(centroid_manifest),
        "centroid_table_sha256": next(item["sha256"] for item in manifest["files"] if item["path"] == "centroids.parquet"),
        "sampling_method": "reservoir_without_replacement_from_centroid_representatives",
        "sampling_rng": "NumPy PCG64", "numpy_version": np.__version__,
        "sampling_seed": sampling_seed, "sample_size": sample_size,
        "compound_set_id": compound_set["compound_set_id"], "compound_set_sha256": compound_set["content_sha256"],
        "selected_compound_keys": [int(row["compound_key"]) for row in selected],
        "selected_compound_ids": [row["compound_id"] for row in selected],
    }
    selection_path = input_dir / "centroid_seed_selection.json"
    write_json(selection_path, selection)
    with collection_scope("manual"):
        run_dir = run_ligand_unidock_pro(
            target_pdb=Path(target_pdb), compound_set_file=compound_set_path,
            center_x=center_x, center_y=center_y, center_z=center_z,
            size_x=size_x, size_y=size_y, size_z=size_z, search_mode=search_mode,
            num_modes=num_modes, seed=seed, compound_limit=sample_size,
            cpu_cores=cpu_cores, gpu_device=gpu_device,
            job_name=f"{record['name']} centroid seed docking",
        )
    metadata = read_json(run_dir / "metadata.json")
    metadata["library_seed_docking"] = {key: value for key, value in selection.items() if key not in {"selected_compound_keys", "selected_compound_ids"}}
    write_json(run_dir / "metadata.json", metadata)
    result = read_json(run_dir / "result.json")
    result.setdefault("outputs", {})["centroid_seed_selection"] = "artifacts/input/centroid_seed_selection.json"
    result.setdefault("metrics", {}).update(library_id=library_id, library_name=record["name"], centroid_sample_size=sample_size,
        centroid_sample_seed=sampling_seed, centroid_manifest_sha256=selection["centroid_manifest_sha256"])
    write_json(run_dir / "result.json", result)
    contract_path = run_dir / "artifact_contract.json"
    contract = read_json(contract_path)
    contract.setdefault("inputs", []).append({"name": "centroid_seed_selection", "type": "library_centroid_selection/v1",
        "path": "artifacts/input/centroid_seed_selection.json", "source_manifest_sha256": selection["centroid_manifest_sha256"]})
    write_json(contract_path, contract)
    return run_dir
