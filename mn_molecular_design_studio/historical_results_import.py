"""Read-only preview of native BindCraft campaign results.

This module deliberately has no production import/commit operation. It reads a
campaign folder or a bounded ZIP archive, reconciles native evidence and emits
a deterministic preview for later explicit review.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import shutil
import sqlite3
import stat
import tempfile
import zipfile
from pathlib import Path, PurePosixPath
from typing import Any
from urllib.parse import quote

from mn_molecular_design_studio.molecule_imports import canonical_components
from mn_molecular_design_studio.output_vocabulary import normalize_metrics

MAX_ARCHIVE_BYTES = 512 * 1024 * 1024
MAX_SOURCE_FILES = 20_000
MAX_TABLE_ROWS = 100_000
MAX_TABLE_BYTES = 64 * 1024 * 1024
_BC1_TABLES = ("final_design_stats.csv", "mpnn_design_stats.csv", "trajectory_stats.csv")
_BC2_STAGES = {
    "1_trajectories": ("!_Trajectories.csv", "trajectory"),
    "2_refolded": ("!_Refolded.csv", "refolded"),
    "3_ranked": ("!_Ranked.csv", "ranked"),
}
_HASH_IN_NAME = re.compile(r"(?:^|_)([0-9a-fA-F]{16})(?=_|$)")
_RANK_PREFIX = re.compile(r"^\d+_")


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _content_identity(root: Path) -> tuple[str, list[dict]]:
    files, rows = [], []
    for current, directories, names in os.walk(root, followlinks=False):
        current_path = Path(current)
        directories[:] = sorted(name for name in directories if not (current_path / name).is_symlink())
        for name in sorted(names):
            path = current_path / name
            if path.is_symlink() or not path.is_file():
                continue
            relative = path.relative_to(root).as_posix()
            digest = _sha256_file(path)
            rows.append({"path": relative, "sha256": digest, "size_bytes": path.stat().st_size})
            files.append(relative.encode("utf-8") + b"\0" + digest.encode("ascii") + b"\n")
            if len(rows) > MAX_SOURCE_FILES:
                raise ValueError(f"Campaign contains more than {MAX_SOURCE_FILES} files.")
    return hashlib.sha256(b"".join(files)).hexdigest(), rows


def _safe_extract_zip(archive_path: Path, destination: Path) -> None:
    if archive_path.stat().st_size > MAX_ARCHIVE_BYTES:
        raise ValueError("Campaign archive exceeds the 512 MB limit.")
    with zipfile.ZipFile(archive_path) as archive:
        members = archive.infolist()
        if len(members) > MAX_SOURCE_FILES:
            raise ValueError(f"Campaign archive contains more than {MAX_SOURCE_FILES} entries.")
        total = sum(item.file_size for item in members)
        if total > MAX_ARCHIVE_BYTES:
            raise ValueError("Campaign archive expands beyond the 512 MB limit.")
        seen_paths = set()
        for item in members:
            raw_name = item.filename.replace("\\", "/")
            name = PurePosixPath(raw_name)
            mode = item.external_attr >> 16
            if (not raw_name or name.is_absolute() or ".." in name.parts
                    or re.match(r"^[A-Za-z]:", raw_name)
                    or stat.S_ISLNK(mode)):
                raise ValueError(f"Campaign archive contains an unsafe path: {item.filename!r}.")
            normalized_name = name.as_posix().rstrip("/")
            if normalized_name in seen_paths:
                raise ValueError(f"Campaign archive contains a duplicate path: {item.filename!r}.")
            seen_paths.add(normalized_name)
            if item.is_dir():
                continue
            target = destination.joinpath(*name.parts)
            target.parent.mkdir(parents=True, exist_ok=True)
            written = 0
            with archive.open(item) as source, target.open("wb") as output:
                while block := source.read(1024 * 1024):
                    written += len(block)
                    if written > item.file_size or written > MAX_ARCHIVE_BYTES:
                        raise ValueError("Campaign archive entry exceeded its declared size.")
                    output.write(block)
            if written != item.file_size:
                raise ValueError("Campaign archive entry size did not match its directory record.")


def _source_root(source: str | Path):
    selected = Path(source).expanduser()
    if selected.is_symlink():
        raise ValueError("Campaign source paths may not be symlinks.")
    path = selected.resolve(strict=True)
    if path.is_dir():
        for current, directories, names in os.walk(path, followlinks=False):
            current_path = Path(current)
            symlink_names = [name for name in [*directories, *names] if (current_path / name).is_symlink()]
            if symlink_names:
                relative = (current_path / sorted(symlink_names)[0]).relative_to(path).as_posix()
                raise ValueError(f"Campaign folders may not contain symlinks: {relative}.")
        return path, None
    if path.is_file() and path.suffix.lower() == ".zip":
        temporary = tempfile.TemporaryDirectory(prefix="studio-bindcraft-preview-")
        try:
            _safe_extract_zip(path, Path(temporary.name))
        except Exception:
            temporary.cleanup()
            raise
        return Path(temporary.name), temporary
    raise ValueError("Choose a BindCraft campaign folder or ZIP archive.")


def _safe_walk_dirs(root: Path):
    for current, directories, _ in os.walk(root, followlinks=False):
        path = Path(current)
        directories[:] = sorted(name for name in directories if not (path / name).is_symlink())
        yield path


def _find_campaign_root(root: Path) -> tuple[Path, str, list[str]]:
    issues: list[str] = []
    dirs = list(_safe_walk_dirs(root))
    # BindCraft 1's three native statistics files are its unambiguous signature.
    for directory in dirs:
        names = {path.name.casefold(): path for path in directory.iterdir() if path.is_file()}
        if all(name in names for name in _BC1_TABLES):
            return directory, "bindcraft1", issues
    candidates = []
    for directory in dirs:
        children = {path.name.casefold(): path for path in directory.iterdir() if path.is_dir()}
        stages = [stage for stage in _BC2_STAGES if stage in children]
        metadata = any(path.is_file() and path.name.casefold().startswith("campaign_metadata")
                       and path.suffix.casefold() == ".json" for path in directory.iterdir())
        tables = any((children.get(stage) and any(p.is_file() and p.name.startswith("!_") and p.suffix.casefold() == ".csv"
                                                   for p in children[stage].iterdir())) for stage in _BC2_STAGES)
        if len(stages) >= 2 or (stages and (metadata or tables)):
            candidates.append((directory, len(stages)))
    if not candidates:
        raise ValueError("No supported BindCraft1 statistics set or BindCraft2 stage layout was found.")
    candidates.sort(key=lambda item: (len(item[0].parts), -item[1], item[0].as_posix()))
    if len(candidates) > 1 and candidates[0][1] == candidates[1][1] and len(candidates[0][0].parts) == len(candidates[1][0].parts):
        raise ValueError("The source contains multiple plausible BindCraft2 campaign roots.")
    if len(candidates) > 1:
        issues.append("Multiple BindCraft2-like folders were found; the shallowest campaign root was selected.")
    return candidates[0][0], "bindcraft2", issues


def _read_csv(path: Path | None) -> tuple[list[dict[str, str]], list[str]]:
    if path is None:
        return [], []
    if path.stat().st_size > MAX_TABLE_BYTES:
        raise ValueError(f"Native table exceeds 64 MB: {path.name}.")
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            return [], []
        headers = [str(name or "").strip() for name in reader.fieldnames]
        folded = [name.casefold() for name in headers]
        if len(set(folded)) != len(folded):
            raise ValueError(f"Native table has duplicate case-insensitive column names: {path.name}.")
        rows = []
        for index, row in enumerate(reader, 1):
            if index > MAX_TABLE_ROWS:
                raise ValueError(f"Native table exceeds {MAX_TABLE_ROWS} rows: {path.name}.")
            rows.append({headers[n]: str(row.get(reader.fieldnames[n]) or "") for n in range(len(headers))})
    return rows, headers


def _value(row: dict, *names: str) -> str:
    folded = {str(key).casefold(): value for key, value in row.items()}
    for name in names:
        value = folded.get(name.casefold())
        if value is not None and str(value).strip():
            return str(value).strip()
    return ""


def _json_objects(root: Path) -> list[tuple[Path, dict]]:
    files = sorted((path for path in root.iterdir() if path.is_file() and path.suffix.casefold() == ".json"
                    and (path.name.casefold() in {"campaign.json", "campaign_metadata.json"}
                         or path.name.casefold().startswith("campaign_metadata"))), key=lambda item: item.name.casefold())
    result = []
    for path in files:
        if path.stat().st_size > MAX_TABLE_BYTES:
            continue
        try:
            value = json.loads(path.read_text(encoding="utf-8-sig"))
        except (OSError, UnicodeDecodeError, json.JSONDecodeError):
            continue
        if isinstance(value, dict):
            result.append((path, value))
    return result


def _walk_find(value: Any, names: set[str]):
    if isinstance(value, dict):
        for key, item in value.items():
            if str(key).casefold() in names and item not in (None, "", [], {}):
                yield str(key), item
            yield from _walk_find(item, names)
    elif isinstance(value, list):
        for item in value:
            yield from _walk_find(item, names)


def _first_metadata_value(objects, names: set[str]):
    for _, obj in objects:
        found = next(_walk_find(obj, names), None)
        if found:
            return found[1]
    return None


def _campaign_metadata(root: Path) -> dict:
    objects = _json_objects(root)
    return {
        "files": [{"path": path.name, "sha256": _sha256_file(path)} for path, _ in objects],
        "objects": objects,
    }


def _numeric_metrics(row: dict, excluded: set[str]) -> dict:
    result = {}
    for key, value in row.items():
        if key.casefold() in excluded or not str(value).strip():
            continue
        try:
            number = float(value)
        except (TypeError, ValueError, OverflowError):
            continue
        if number == number and abs(number) != float("inf"):
            result[key] = int(number) if number.is_integer() else number
    return result


def _sequence_identity(sequence: str) -> tuple[list[dict], str, str] | None:
    pieces = [re.sub(r"\s+", "", part).upper().rstrip("*") for part in str(sequence or "").split("/")]
    if not pieces or any(not part for part in pieces):
        return None
    chains = []
    for index, part in enumerate(pieces):
        chain = chr(65 + index) if index < 26 else f"C{index + 1}"
        chains.append({"type": "protein", "chain": chain, "role": "binder", "sequence": part})
    try:
        return canonical_components(chains)
    except ValueError:
        return None


def _native_status(row: dict) -> tuple[str, str, str]:
    outcome = _value(row, "outcome", "status", "filter_status", "result")
    failed = _value(row, "failed_filters", "filters_failed", "filter_failures")
    normalized = outcome.casefold()
    if normalized in {"passed", "pass", "accepted"} and not failed:
        return "passed", outcome, failed
    if normalized in {"failed", "fail", "rejected"} or failed:
        return "rejected", outcome, failed
    return "unknown", outcome, failed


def _exact_table(stage_dir: Path | None, expected: str, warnings: list[str]) -> Path | None:
    if stage_dir is None:
        return None
    files = sorted((path for path in stage_dir.iterdir() if path.is_file() and path.suffix.casefold() == ".csv"),
                   key=lambda item: item.name.casefold())
    exact = [path for path in files if path.name.casefold() == expected.casefold()]
    if len(exact) == 1:
        return exact[0]
    fallback = [path for path in files if path.name.startswith("!_")]
    if len(fallback) == 1:
        warnings.append(f"{expected} was absent; used the only native stage table {fallback[0].name}.")
        return fallback[0]
    if len(fallback) > 1:
        warnings.append(f"{stage_dir.name} has multiple !_*.csv tables; its table was not selected.")
    elif files:
        warnings.append(f"No recognized native stage table was found in {stage_dir.name}.")
    return None


def _stage_dirs(root: Path) -> dict[str, Path | None]:
    children = {path.name.casefold(): path for path in root.iterdir() if path.is_dir()}
    return {stage: children.get(stage) for stage in _BC2_STAGES}


def _native_hash(row: dict) -> str:
    value = _value(row, "hash", "trajectory_hash")
    if value:
        return value.casefold()
    design = _value(row, "design")
    match = _HASH_IN_NAME.search(design)
    return match.group(1).casefold() if match else ""


def _find_exact_structure(root: Path, stage: str, design: str, *, allow_rank_prefix=False) -> tuple[str, str, list[str]]:
    if not design:
        return "", "", []
    suffixes = {".cif", ".mmcif", ".pdb", ".ent"}
    stage_dir = next((path for path in root.iterdir() if path.is_dir() and path.name.casefold() == stage.casefold()), None)
    if not stage_dir.is_dir():
        return "", "", []
    children = {path.name.casefold(): path for path in stage_dir.iterdir() if path.is_dir()}
    bases = [children.get("complexes", stage_dir / "Complexes"), stage_dir]
    candidates = []
    wanted = {design.casefold()}
    if stage == "3_Ranked":
        wanted.add((design + "_seq0").casefold())
    for base in bases:
        if not base.is_dir():
            continue
        for path in base.rglob("*"):
            if not path.is_file() or path.suffix.casefold() not in suffixes:
                continue
            stem = path.stem.casefold()
            if stem in wanted:
                candidates.append(path)
            elif allow_rank_prefix:
                normalized = _RANK_PREFIX.sub("", path.stem, count=1).casefold()
                if normalized in wanted:
                    candidates.append(path)
    candidates = sorted(set(candidates), key=lambda item: item.relative_to(root).as_posix().casefold())
    if len(candidates) > 1:
        return "", "", [f"Ambiguous exact structure files for native design {design}."]
    if not candidates:
        return "", "", []
    path = candidates[0]
    return path.relative_to(root).as_posix(), _sha256_file(path), []


def _inspect_native_structure(root: Path, relative_path: str, binder_sequence: str, target: dict,
                              *, model_id: str = "") -> dict:
    """Validate native coordinates, model count, and exact chain sequences."""
    path = (root / relative_path).resolve(strict=True)
    if not path.is_relative_to(root.resolve()):
        return {"path": relative_path, "sha256": _sha256_file(path), "status": "unsafe_path"}
    digest = _sha256_file(path)
    try:
        from mn_molecular_design_studio.molecule_imports import parse_structure
        parsed = parse_structure(path.read_bytes(), path.name, require_model_selection=True)[0]
    except (OSError, ValueError) as exc:
        return {"path": relative_path, "sha256": digest, "status": "invalid_coordinates",
                "validation_error": str(exc)}
    coordinate_model_id = str(parsed.get("model_id") or "")
    selected_model_id = str(model_id or coordinate_model_id)
    proteins = {str(item.get("chain") or ""): str(item.get("sequence") or "")
                for item in parsed["components"] if item.get("type") == "protein"}
    expected_binders = [re.sub(r"\s+", "", part).upper() for part in re.split(r"[/ :]", binder_sequence) if part]
    if not expected_binders:
        return {"path": relative_path, "sha256": digest, "status": "binder_sequence_unavailable",
                "model_id": selected_model_id, "coordinate_model_id": coordinate_model_id,
                "chain_sequences": proteins}
    observed_sequences = list(proteins.values())
    unmatched = list(observed_sequences)
    binder_chains = []
    for expected in expected_binders:
        matches = [chain for chain, sequence in proteins.items() if sequence == expected and chain not in binder_chains]
        if len(matches) != 1:
            return {"path": relative_path, "sha256": digest,
                    "status": "binder_sequence_mismatch" if not matches else "ambiguous_binder_chain",
                    "model_id": selected_model_id, "coordinate_model_id": coordinate_model_id,
                    "chain_sequences": proteins}
        binder_chains.append(matches[0])
        unmatched.remove(expected)
    target_identity = target.get("identity") if isinstance(target, dict) else None
    target_components = target_identity.get("components") if isinstance(target_identity, dict) else []
    expected_targets = [str(item.get("sequence") or "") for item in target_components
                        if isinstance(item, dict) and item.get("type") == "protein"]
    if unmatched:
        if not expected_targets:
            return {"path": relative_path, "sha256": digest, "status": "target_correspondence_unresolved",
                    "model_id": selected_model_id, "coordinate_model_id": coordinate_model_id,
                    "chain_sequences": proteins,
                    "binder_chains": binder_chains}
        remaining = list(unmatched)
        target_chains = []
        for expected in expected_targets:
            matches = [chain for chain, sequence in proteins.items()
                       if sequence == expected and chain not in binder_chains and chain not in target_chains]
            if len(matches) != 1:
                return {"path": relative_path, "sha256": digest, "status": "target_sequence_mismatch",
                        "model_id": selected_model_id, "coordinate_model_id": coordinate_model_id,
                        "chain_sequences": proteins,
                        "binder_chains": binder_chains}
            target_chains.append(matches[0])
            remaining.remove(expected)
        if remaining:
            return {"path": relative_path, "sha256": digest, "status": "unexpected_protein_chains",
                    "model_id": selected_model_id, "coordinate_model_id": coordinate_model_id,
                    "chain_sequences": proteins,
                    "binder_chains": binder_chains, "target_chains": target_chains}
    else:
        target_chains = []
    return {"path": relative_path, "sha256": digest, "status": "available",
            "model_id": selected_model_id, "coordinate_model_id": coordinate_model_id,
            "chain_sequences": proteins,
            "binder_chains": binder_chains, "target_chains": target_chains,
            "validation": "exact_protein_sequences_and_single_model"}


def _row_candidate(row: dict, *, stage: str, row_index: int, engine: str, root: Path,
                   candidate_key: str, warnings: list[str], target: dict | None = None) -> dict:
    sequence = _value(row, "Binder_Sequence", "binder_sequence", "sequence", "designed_sequence")
    normalized_sequence = "/".join(re.sub(r"\s+", "", part).upper() for part in sequence.split("/"))
    identity = _sequence_identity(normalized_sequence)
    design = _value(row, "design", "binder_name", "name", "candidate_id", "id")
    trajectory_hash = _native_hash(row)
    if not trajectory_hash and engine == "bindcraft2":
        trajectory_hash = _value(row, "trajectory")
    status, outcome, failed = _native_status(row)
    metric_values = _numeric_metrics(row, {"design", "binder_sequence", "sequence", "designed_sequence", "outcome", "status", "failed_filters", "filters_failed", "filter_failures"})
    metric_projection = normalize_metrics(metric_values, tool=engine, context={
        "candidate": {"source_tool": engine, "stage": "complex_refolding", "binder_sequence": normalized_sequence,
                      "raw_metadata": {"engine_version": ""}}
    })
    structure_stage = ("3_Ranked" if stage == "ranked" else "2_Refolded") if engine == "bindcraft2" else stage
    structure_path, structure_digest, structure_warnings = _find_exact_structure(
        root, structure_stage, design, allow_rank_prefix=(stage == "ranked")) if engine == "bindcraft2" else ("", "", [])
    warnings.extend(structure_warnings)
    structure = (_inspect_native_structure(root, structure_path, normalized_sequence, target or {})
                 if structure_path else {"path": None, "sha256": None, "status": "missing_or_unmatched"})
    structure["stage"] = "ranked" if stage == "ranked" else "refolded"
    unsupported_fields = sorted(set(metric_values) - {
                         str(item.get("source_field", "")).split(".")[-1]
                         for entries in metric_projection["provenance"].values()
                         for item in entries.get("aliases", [])
                     })
    compatibility = {"status": "partial" if unsupported_fields or metric_projection["issues"] else "normalized_subset",
                     "normalized_fields": sorted(metric_projection["metrics"]),
                     "unsupported_native_fields": unsupported_fields}
    return {
        "candidate_id": candidate_key,
        "native_design_id": design,
        "trajectory_hash": trajectory_hash,
        "native_stage": stage,
        "sequence": normalized_sequence,
        "identity_id": "mol-" + identity[1][:32] if identity else None,
        "identity_fingerprint": identity[1] if identity else None,
        "native_filter_status": status,
        "native_outcome": outcome or None,
        "native_failed_filters": failed or None,
        "ranked": stage == "ranked",
        "rank": _value(row, "rank", "ranking") or None,
        "native_metrics": dict(row),
        "numeric_metrics_for_normalization": metric_values,
        "normalized_metrics": metric_projection["metrics"],
        "metric_provenance": metric_projection["provenance"],
        "metric_conflicts": metric_projection["conflicts"],
        "metric_issues": metric_projection["issues"],
        "unrecognized_metrics": sorted(set(metric_values) - {
            str(item.get("source_field", "")).split(".")[-1]
            for entries in metric_projection["provenance"].values() for item in entries.get("aliases", [])
        }),
        "structure": structure,
        "structures": [structure] if structure_path else [],
        "native_metric_compatibility": compatibility,
        "retention_selected": False,
    }


def _bindcraft2(root: Path, source_id: str, warnings: list[str], target: dict | None = None) -> dict:
    stage_dirs = _stage_dirs(root)
    tables, headers = {}, {}
    for stage, (filename, _) in _BC2_STAGES.items():
        path = _exact_table(stage_dirs[stage], filename, warnings)
        rows, columns = _read_csv(path)
        tables[stage] = rows
        headers[stage] = columns
        if stage_dirs[stage] is None:
            warnings.append(f"BindCraft2 stage folder {stage} is missing.")
        elif path is None:
            warnings.append(f"BindCraft2 stage table for {stage} is unavailable.")

    trajectory_rows, refolded_rows, ranked_rows = (tables[key] for key in _BC2_STAGES)
    # The scientific crosswalk is the native pair (trajectory hash, binder sequence).
    # Duplicate pairs are deliberately left ambiguous instead of picking a row.
    refolds_by_pair: dict[tuple[str, str], list[dict]] = {}
    candidates = []
    for index, row in enumerate(refolded_rows, 1):
        digest = _native_hash(row)
        if not digest:
            digest = _value(row, "trajectory")
        sequence = _value(row, "Binder_Sequence", "binder_sequence", "sequence")
        normalized = "/".join(re.sub(r"\s+", "", part).upper() for part in sequence.split("/"))
        pair = (digest, normalized)
        candidate = _row_candidate(row, stage="refolded", row_index=index, engine="bindcraft2", root=root,
                                   candidate_key="historical-" + hashlib.sha256(
                                       f"{source_id}\0refolded\0{index}\0{digest}\0{normalized}".encode()).hexdigest()[:24],
                                   warnings=warnings, target=target)
        if digest and normalized:
            refolds_by_pair.setdefault(pair, []).append(candidate)
        else:
            candidate["association_status"] = "missing_hash_or_sequence"
            warnings.append(f"Refolded row {index} has no exact trajectory hash and binder sequence key.")
        candidates.append(candidate)

    ranked_candidates = []
    for index, row in enumerate(ranked_rows, 1):
        digest = _native_hash(row)
        sequence = _value(row, "Binder_Sequence", "binder_sequence", "sequence")
        normalized = "/".join(re.sub(r"\s+", "", part).upper() for part in sequence.split("/"))
        matches = refolds_by_pair.get((digest, normalized), []) if digest and normalized else []
        if len(matches) == 1:
            candidate = matches[0]
            candidate["ranked"] = True
            candidate["rank"] = _value(row, "rank", "ranking") or str(index)
            candidate["ranked_design_id"] = _value(row, "design") or None
            candidate["ranked_native_metrics"] = dict(row)
            ranked_structure, ranked_digest, ranked_warnings = _find_exact_structure(
                root, "3_Ranked", _value(row, "design"), allow_rank_prefix=True)
            warnings.extend(ranked_warnings)
            ranked_structure_record = (_inspect_native_structure(root, ranked_structure, normalized, target or {})
                                       if ranked_structure else {"path": None, "sha256": None,
                                                                  "status": "missing_or_unmatched"})
            ranked_structure_record["stage"] = "ranked"
            ranked_values = _numeric_metrics(row, {"design", "binder_sequence", "sequence", "outcome", "status", "failed_filters"})
            ranked_projection = normalize_metrics(ranked_values, tool="bindcraft2", context={
                "candidate": {"source_tool": "bindcraft2", "stage": "complex_refolding", "binder_sequence": normalized,
                              "raw_metadata": {"engine_version": ""}}
            })
            candidate["ranked_metrics"] = ranked_values
            candidate["ranked_normalized_metrics"] = ranked_projection["metrics"]
            candidate["ranked_metric_provenance"] = ranked_projection["provenance"]
            candidate["ranked_metric_conflicts"] = ranked_projection["conflicts"]
            candidate["ranked_structure"] = ranked_structure_record
            if ranked_structure:
                candidate.setdefault("structures", []).append(ranked_structure_record)
            if ranked_structure_record["status"] != "available":
                warnings.append(f"Ranked structure for {candidate['candidate_id']} is unresolved: {ranked_structure_record['status']}.")
            candidate["association_status"] = "exact_hash_and_sequence"
            continue
        candidate = _row_candidate(row, stage="ranked", row_index=index, engine="bindcraft2", root=root,
                                   candidate_key="historical-" + hashlib.sha256(
                                       f"{source_id}\0ranked\0{index}\0{digest}\0{normalized}".encode()).hexdigest()[:24],
                                   warnings=warnings, target=target)
        candidate["association_status"] = "unmatched_refolded" if not matches else "ambiguous_refolded_matches"
        candidate["native_filter_status"] = "unknown"
        candidate["native_outcome"] = None
        candidate["native_failed_filters"] = None
        if len(matches) > 1:
            warnings.append(f"Ranked row {index} matches multiple refolds by exact hash and sequence; no crosswalk was chosen.")
        else:
            warnings.append(f"Ranked row {index} has no refolded row matching its exact hash and sequence.")
        ranked_candidates.append(candidate)
        candidates.append(candidate)
    for pair, rows in refolds_by_pair.items():
        if len(rows) > 1:
            for row in rows:
                row["association_status"] = "duplicate_refold_key"
            warnings.append(f"Multiple refold rows share exact key (hash={pair[0]}, sequence); ranked joins remain unresolved.")

    trajectory_index = {}
    for index, row in enumerate(trajectory_rows, 1):
        digest = _native_hash(row)
        trajectory_index.setdefault(digest, []).append({"row_index": index, "design": _value(row, "design"), "native": dict(row)})
    for candidate in candidates:
        digest = candidate.get("trajectory_hash") or ""
        matches = trajectory_index.get(digest, []) if digest else []
        candidate["trajectory_association"] = "exact_hash" if len(matches) == 1 else (
            "ambiguous_hash" if len(matches) > 1 else "missing_hash_match")
        if len(matches) > 1:
            warnings.append(f"Candidate {candidate['candidate_id']} matches multiple trajectory rows by hash.")
    counts = {
        "trajectories_started": len(trajectory_rows),
        "refolded_candidates": len(refolded_rows),
        "native_filter_passes": sum(_native_status(row)[0] == "passed" for row in refolded_rows),
        "native_filter_rejected": sum(_native_status(row)[0] == "rejected" for row in refolded_rows),
        "native_filter_unknown": sum(_native_status(row)[0] == "unknown" for row in refolded_rows),
        "ranked_designs": len(ranked_rows),
    }
    if not any("outcome" in name.casefold() for name in headers["2_refolded"]):
        warnings.append("The refolded table has no outcome column; native pass/reject status remains unknown.")
    return {"counts": counts, "candidates": candidates, "trajectory_rows": trajectory_rows,
            "tables": {"trajectories": headers["1_trajectories"], "refolded": headers["2_refolded"], "ranked": headers["3_ranked"]}}


def _bindcraft1_structure_paths(root: Path, design: str) -> list[dict]:
    """Recognize native <Design>_modelN and ranked-prefix model layouts."""
    rows = []
    exact = re.compile(r"^" + re.escape(design) + r"_model_?(\d+)$", re.IGNORECASE)
    for path in root.rglob("*"):
        if not path.is_file() or path.suffix.casefold() not in {".pdb", ".cif", ".mmcif", ".ent"}:
            continue
        stem = path.stem
        rank_prefix = _RANK_PREFIX.match(stem)
        normalized = stem[rank_prefix.end():] if rank_prefix else stem
        model_match = exact.fullmatch(normalized)
        layout = "native_model_suffix" if model_match else "exact_design_stem" if normalized.casefold() == design.casefold() else ""
        if not layout:
            continue
        rows.append({"path": path.relative_to(root).as_posix(),
                     "model_id": model_match.group(1) if model_match else "",
                     "layout": layout, "rank_prefix": rank_prefix.group(0) if rank_prefix else ""})
    return sorted(rows, key=lambda item: (item["model_id"], item["path"].casefold()))


def _bindcraft1(root: Path, source_id: str, warnings: list[str], target: dict | None = None) -> dict:
    table_paths = {path.name.casefold(): path for path in root.iterdir() if path.is_file()}
    parsed, headers = {}, {}
    for filename in _BC1_TABLES:
        rows, columns = _read_csv(table_paths.get(filename))
        parsed[filename] = rows
        headers[filename] = columns
    final_rows = parsed["final_design_stats.csv"]
    candidates = []
    for index, row in enumerate(final_rows, 1):
        design = _value(row, "design", "binder_name", "name", "candidate_id", "id")
        sequence = _value(row, "Binder_Sequence", "binder_sequence", "sequence", "designed_sequence")
        if not sequence:
            warnings.append(f"BindCraft1 final-design row {index} has no recognizable binder sequence.")
        candidate = _row_candidate(row, stage="final_design", row_index=index, engine="bindcraft1", root=root,
                                   candidate_key="historical-" + hashlib.sha256(
                                       f"{source_id}\0final_design\0{index}\0{design}\0{sequence}".encode()).hexdigest()[:24],
                                   warnings=warnings, target=target)
        candidate["native_filter_status"] = _native_status(row)[0]
        candidate["association_status"] = "native_final_design_table"
        matched_paths = _bindcraft1_structure_paths(root, design)
        structures = []
        for item in matched_paths:
            structure = _inspect_native_structure(root, item["path"], sequence or "", target or {},
                                                   model_id=item["model_id"])
            structure.update(model_id=item["model_id"] or structure.get("model_id"),
                             layout=item["layout"], rank_prefix=item["rank_prefix"], stage="final_design")
            structures.append(structure)
        candidate["structures"] = structures
        if len(structures) == 1:
            candidate["structure"] = structures[0]
        elif len(structures) > 1:
            candidate["structure"] = {"path": None, "sha256": None, "status": "multiple_models",
                                       "model_ids": [item.get("model_id") for item in structures]}
        else:
            candidate["structure"] = {"path": None, "sha256": None, "status": "missing_or_unmatched"}
            if sequence:
                warnings.append(f"No supported native BindCraft1 exact-design or _modelN structure matches {design}.")
        metric_names = set(candidate.get("unrecognized_metrics") or [])
        if any("average_iptm" in key.casefold() for key in metric_names | set(row)):
            warnings.append(f"BindCraft1 {design} has averaged iPTM evidence preserved as native data; no verified normalized mapping is available.")
        if structures and all(item.get("status") == "available" for item in structures):
            candidate["structure_layout_status"] = "supported_native_models"
        elif structures:
            candidate["structure_layout_status"] = "recognized_layout_unvalidated_coordinates"
        else:
            candidate["structure_layout_status"] = "unsupported_or_missing_layout"
        candidates.append(candidate)
    if not final_rows:
        warnings.append("BindCraft1 final_design_stats.csv contains no design rows.")
    counts = {"trajectories_started": len(parsed["trajectory_stats.csv"]),
              "refolded_candidates": len(parsed["mpnn_design_stats.csv"]),
              "native_filter_passes": sum(_native_status(row)[0] == "passed" for row in final_rows),
              "native_filter_rejected": sum(_native_status(row)[0] == "rejected" for row in final_rows),
              "native_filter_unknown": sum(_native_status(row)[0] == "unknown" for row in final_rows),
              "ranked_designs": 0,
              "final_designs": len(final_rows)}
    return {"counts": counts, "candidates": candidates, "trajectory_rows": parsed["trajectory_stats.csv"],
            "tables": {"final_design_stats": headers["final_design_stats.csv"],
                       "mpnn_design_stats": headers["mpnn_design_stats.csv"],
                       "trajectory_stats": headers["trajectory_stats.csv"]}}


def _target_summary(root: Path, campaign_metadata: dict, source_files: list[dict], warnings: list[str]) -> dict:
    objects = campaign_metadata["objects"]
    target_records = []
    for _, obj in objects:
        for _, value in _walk_find(obj, {"targets"}):
            if isinstance(value, list):
                target_records.extend(item for item in value if isinstance(item, dict))
            elif isinstance(value, dict):
                target_records.append(value)
    target_record = target_records[0] if target_records else {}
    target_sequence = next((target_record.get(key) for key in ("target_sequence", "sequence", "fasta_sequence")
                            if target_record.get(key)), None)
    if target_sequence is None:
        target_sequence = _first_metadata_value(objects, {"target_sequence", "target_fasta_sequence"})
    if isinstance(target_sequence, str):
        target_sequence = re.sub(r"\s+", "", target_sequence).upper().rstrip("*")
    elif isinstance(target_sequence, dict):
        target_sequence = {str(chain): re.sub(r"\s+", "", str(sequence)).upper().rstrip("*")
                           for chain, sequence in target_sequence.items()}
    else:
        target_sequence = None
    target_value = next((target_record.get(key) for key in ("target_path", "target_file", "target_structure", "construct_path")
                         if target_record.get(key)), None) or _first_metadata_value(objects, {"target_path", "target_file", "target_structure", "construct_path"})
    chain_value = next((target_record.get(key) for key in ("chains", "target_chains", "chain_ids") if target_record.get(key)), None)
    chain_value = chain_value or _first_metadata_value(objects, {"target_chains", "chain_ids"})
    hotspot_value = next((target_record.get(key) for key in ("hotspots", "hotspot_residues") if target_record.get(key)), None)
    hotspot_value = hotspot_value or _first_metadata_value(objects, {"hotspots", "hotspot_residues"})
    target_path = None
    source_status = "not_declared"
    digest = None
    coordinate_validation = "not_applicable"
    sequence_source = "campaign_metadata" if target_sequence else None
    if isinstance(target_value, str) and target_value.strip():
        raw_path = Path(target_value).expanduser()
        possible = raw_path if raw_path.is_absolute() else root / raw_path
        try:
            resolved = possible.resolve(strict=True)
            if resolved.is_relative_to(root.resolve()) and resolved.is_file():
                target_path = resolved
                source_status = "available_inside_campaign"
                digest = _sha256_file(resolved)
                if resolved.suffix.casefold() in {".pdb", ".ent", ".cif", ".mmcif"}:
                    try:
                        from mn_molecular_design_studio.molecule_imports import parse_structure
                        parsed_target = parse_structure(resolved.read_bytes(), resolved.name,
                                                        require_model_selection=True)[0]
                        sequences = {str(component.get("chain") or ""): str(component.get("sequence") or "")
                                     for component in parsed_target.get("components") or []
                                     if component.get("type") == "protein"}
                        selected_chains = ([str(value) for value in chain_value]
                                           if isinstance(chain_value, (list, tuple))
                                           else [value.strip() for value in str(chain_value or "").split(",") if value.strip()])
                        if selected_chains:
                            if set(selected_chains) - set(sequences):
                                raise ValueError("Declared target chains are absent from the source coordinates.")
                            observed = {chain: sequences[chain] for chain in selected_chains}
                        elif isinstance(target_sequence, str) and target_sequence:
                            matches = {chain: sequence for chain, sequence in sequences.items()
                                       if sequence == target_sequence}
                            if len(matches) != 1:
                                raise ValueError("Target sequence does not map to exactly one source coordinate chain.")
                            observed = matches
                        else:
                            observed = sequences
                        expected = ({str(chain): str(sequence).upper() for chain, sequence in target_sequence.items()}
                                    if isinstance(target_sequence, dict) else None)
                        if expected:
                            if any(observed.get(chain) != sequence for chain, sequence in expected.items()):
                                raise ValueError("Target sequence metadata does not match the declared source chains.")
                        elif isinstance(target_sequence, str) and target_sequence:
                            if sorted(observed.values()) != [target_sequence]:
                                raise ValueError("Target sequence metadata does not match source coordinates.")
                        elif observed:
                            target_sequence = observed
                            sequence_source = "source_structure_chain_sequences"
                        if not observed:
                            raise ValueError("No protein chains were found in the target coordinates.")
                        coordinate_validation = "exact_chain_sequence_single_model"
                    except (OSError, ValueError, TypeError) as exc:
                        coordinate_validation = "invalid_or_mismatched"
                        source_status = "invalid_or_mismatched_coordinates"
                        warnings.append(f"Target source coordinates were not validated: {exc}")
                elif not target_sequence and resolved.suffix.casefold() in {".fa", ".fasta", ".faa"}:
                    text = resolved.read_text(encoding="utf-8-sig")
                    target_sequence = "".join(line.strip() for line in text.splitlines() if line.strip() and not line.startswith(">" )).upper()
                    sequence_source = "source_fasta"
            else:
                source_status = "external_or_outside_selected_campaign"
                warnings.append("The declared target source is outside the selected campaign; it was not opened.")
        except (OSError, RuntimeError):
            source_status = "missing_or_inaccessible"
            warnings.append("The declared target source is unavailable; no external path was followed.")
    if isinstance(target_sequence, str) and target_sequence and not re.fullmatch(r"[A-Z]+", target_sequence):
        warnings.append("Target sequence metadata contains non-sequence characters and was retained as unresolved.")
        target_sequence = None
    target_components = []
    if isinstance(target_sequence, str) and target_sequence:
        if isinstance(chain_value, (list, tuple)):
            chain_id = str(chain_value[0]) if chain_value else "A"
        else:
            chain_id = str(chain_value or "A").split(",")[0].strip() or "A"
        target_components = [{"type": "protein", "chain": chain_id,
                              "role": "target", "sequence": target_sequence}]
    elif isinstance(target_sequence, dict):
        target_components = [{"type": "protein", "chain": str(chain), "role": "target", "sequence": sequence}
                             for chain, sequence in target_sequence.items()]
    target_identity = None
    if target_components:
        try:
            components, fingerprint, _ = canonical_components(target_components)
            target_identity = {"identity_id": "mol-" + fingerprint[:32], "fingerprint": fingerprint,
                              "components": components}
        except ValueError:
            warnings.append("Target sequence is present but could not be canonicalized as a protein identity.")
    return {"sequence_or_chain_sequences": target_sequence,
            "identity": target_identity,
            "sequence_source": sequence_source,
            "identity_components": target_components,
            "chain_mapping": chain_value,
            "hotspots": hotspot_value,
            "source_structure": {"path": target_path.relative_to(root).as_posix() if target_path else None,
                                 "sha256": digest,
                                 "status": source_status,
                                 "validation": coordinate_validation,
                                 "coordinate_type": "source_coordinates_unspecified" if target_path and target_path.suffix.casefold() in {".pdb", ".ent", ".cif", ".mmcif"} else None,
                                 "experimental_or_native_claim": False}}


def _existing_identity_matches(database_path: str | Path | None, user: dict | None,
                               identities: list[dict]) -> dict[str, list[dict]]:
    if not database_path or not user:
        return {}
    from mn_molecular_design_studio.molecular_database import MolecularDatabase
    path = Path(database_path).expanduser().resolve()
    if not path.is_file():
        return {}
    uri = "file:" + quote(path.as_posix(), safe="/") + "?mode=ro"
    connection = sqlite3.connect(uri, uri=True)
    connection.row_factory = sqlite3.Row
    try:
        available = {row[0] for row in connection.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        if not {"molecules", "molecule_members", "molecule_groups"} <= available:
            return {}
        result = {}
        for identity in identities:
            matches = []
            rows = connection.execute("""SELECT DISTINCT m.id,m.fingerprint,mm.name,mm.group_id AS member_group_id,g.*
                FROM molecules m JOIN molecule_members mm ON mm.molecule_id=m.id
                JOIN molecule_groups g ON g.id=mm.group_id WHERE m.fingerprint=?""",
                                      (identity["fingerprint"],)).fetchall()
            for row in rows:
                if MolecularDatabase._can_read(row, user):
                    matches.append({"molecule_id": row["id"], "name": row["name"], "group_id": row["member_group_id"]})
            if matches:
                result[identity["identity_id"]] = matches
        return result
    finally:
        connection.close()


def preview_historical_campaign(source: str | Path, *, selected_candidate_ids=None,
                                database_path: str | Path | None = None,
                                user: dict | None = None) -> dict:
    """Return a deterministic, non-mutating BindCraft1/2 database preview.

    ``selected_candidate_ids`` is an explicit retention selection for the
    preview only. This Phase 1 reader has no Database write path.
    """
    root_container, temporary = _source_root(source)
    try:
        campaign_root, engine, warnings = _find_campaign_root(root_container)
        source_id, source_files = _content_identity(campaign_root)
        metadata = _campaign_metadata(campaign_root)
        target = _target_summary(campaign_root, metadata, source_files, warnings)
        result = (_bindcraft2(campaign_root, source_id, warnings, target) if engine == "bindcraft2"
                  else _bindcraft1(campaign_root, source_id, warnings, target))
        candidates = result["candidates"]
        candidate_ids = [row["candidate_id"] for row in candidates]
        if len(candidate_ids) != len(set(candidate_ids)):
            raise ValueError("The campaign produced duplicate preview candidate identifiers.")
        selected = set(str(value) for value in (selected_candidate_ids or []))
        unknown = selected - set(candidate_ids)
        if unknown:
            raise ValueError("Selected candidate IDs are not present in this campaign preview.")
        for row in candidates:
            row["retention_selected"] = row["candidate_id"] in selected
        unique = {}
        for row in candidates:
            identity_id = row.get("identity_id")
            fingerprint = row.get("identity_fingerprint")
            if not identity_id or not fingerprint:
                continue
            identity = unique.setdefault(identity_id, {"identity_id": identity_id, "fingerprint": fingerprint,
                                                       "components": _sequence_identity(row["sequence"])[0],
                                                       "candidate_ids": [], "retention_selected": False})
            identity["candidate_ids"].append(row["candidate_id"])
            identity["retention_selected"] = identity["retention_selected"] or row["retention_selected"]
        identities = sorted(unique.values(), key=lambda item: item["identity_id"])
        matches = _existing_identity_matches(database_path, user, identities)
        for identity in identities:
            identity["existing_matches"] = matches.get(identity["identity_id"], [])
            identity["database_action"] = "reuse_existing_identity" if identity["existing_matches"] else "propose_new_identity"
        for row in candidates:
            row["existing_identity_matches"] = matches.get(row.get("identity_id"), [])
        selected_rows = [row for row in candidates if row["retention_selected"]]
        proposed_eval = []
        proposed_structures = []
        unresolved_structures = []
        for row in selected_rows:
            structures = row.get("structures") or ([row.get("structure")] if row.get("structure") else [])
            for structure in structures:
                if not isinstance(structure, dict) or not structure.get("path"):
                    continue
                if structure.get("status") != "available":
                    unresolved_structures.append({"candidate_id": row["candidate_id"], **structure})
                    continue
                proposed_structures.append({"candidate_id": row["candidate_id"], **structure})
                if row.get("identity_id") and structure.get("sha256"):
                    proposed_eval.append({"candidate_id": row["candidate_id"], "molecule_id": row["identity_id"],
                                          "structure_digest": structure["sha256"], "model_id": structure.get("model_id"),
                                          "stage": structure.get("stage"), "evaluator": engine,
                                          "status": "preview_only_not_attached", "native_metrics_preserved": True})
        return {
            "status": "preview_only",
            "production_import_enabled": False,
            "engine": engine,
            "engine_version": _first_metadata_value(metadata["objects"], {"version", "bindcraft_version"}) or "unknown",
            "source_content_id": source_id,
            "import_key": f"bindcraft-import:{source_id}",
            "source_file_count": len(source_files),
            "source_files": source_files,
            "campaign_root_name": campaign_root.name,
            "campaign_metadata": {"files": metadata["files"],
                                  "resolved_settings": next((obj for path, obj in metadata["objects"] if path.name.casefold().startswith("campaign_metadata")), {}),
                                  "requested_settings": next((obj for path, obj in metadata["objects"] if path.name.casefold() == "campaign.json"), {})},
            "target": target,
            "counts": result["counts"],
            "tables": result["tables"],
            "candidates": candidates,
            "proposed_identities": identities,
            "proposed_structures": proposed_structures,
            "unresolved_structures": unresolved_structures,
            "proposed_evaluations": proposed_eval,
            "evaluations_reused": [],
            "proposed_collection": {"name": f"{engine} · {campaign_root.name}", "selected_candidate_ids": sorted(selected),
                                    "status": "preview_only_not_created"},
            "duplicates": [identity for identity in identities if len(identity["candidate_ids"]) > 1],
            "conflicts": [row for row in candidates if row.get("association_status") in {
                "missing_hash_or_sequence", "unmatched_refolded", "ambiguous_refolded_matches", "duplicate_refold_key"}],
            "warnings": sorted(set(warnings)),
            "source_provenance": {"source_content_id": source_id, "import_key": f"bindcraft-import:{source_id}", "engine": engine,
                                  "selected_candidate_ids": sorted(selected),
                                  "native_campaign_metrics_unchanged": True},
        }
    finally:
        if temporary is not None:
            temporary.cleanup()
