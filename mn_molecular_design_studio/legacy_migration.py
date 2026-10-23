"""Read-only planning and safe staging for MN-Ligand run migration.

MN-Ligand jobs are filesystem-backed and use several generations of metadata.
This adapter follows recorded run references from each compound import, then
selects only thyroid-receptor result jobs and their prerequisite run folders.
The planner is deliberately read-only; copying/registration is a separate,
explicit operation after its manifest has been reviewed.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import sys
import tempfile
from datetime import datetime, timezone
from collections import defaultdict, deque
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Callable, Iterable


DEFAULT_SOURCE_ROOT = Path("/mnt/data/RESULTS/mn-ligand-workdir/workdir/runs")


@dataclass(frozen=True)
class DatasetProject:
    source_name: str
    short_name: str
    project_name: str
    group_names: tuple[str, ...] = ()


DATASET_PROJECTS = (
    DatasetProject(
        "HY-L126 Nuclear Receptor Compound Library",
        "HY-L126",
        "Thyroid Hormon Recptor interaction HY-L126",
        ("HY-L126 Nuclear Receptor Compound Library",),
    ),
    DatasetProject(
        "Medchem THR agonist list",
        "Medchem",
        "Thyroid Hormon Recptor interaction Medchem",
        ("Medchem THR agonist list",),
    ),
    DatasetProject(
        "screening ultra-large library top 10 structure-based",
        "Large library screen",
        "Thyroid Hormon Recptor interaction Large library screen",
        (
            "Screening ultra-large library · top 10 structure-based",
            "screening ultra-large library top 10 structure-based",
        ),
    ),
)


# PDB entries used in the source runs, grouped by receptor subtype. Multiple
# structures are kept distinct because their preparation and ligand pockets
# differ even when they represent the same protein.
THR_TARGETS = {
    "2H79": {"receptor": "TRα", "name": "Thyroid hormone receptor alpha"},
    "4LNW": {"receptor": "TRα", "name": "Thyroid hormone receptor alpha"},
    "3GWS": {"receptor": "TRβ", "name": "Thyroid hormone receptor beta"},
    "6KKB": {"receptor": "TRβ", "name": "Thyroid hormone receptor beta"},
}


PRIMARY_RESULT_TYPES = {
    "docking_campaign",
    "openvs_docking",
    "complex_refolding",
    "learned_binding_affinity",
}
FOLLOWUP_RESULT_TYPES = {
    "interaction_analysis",
    "pose_validation",
    "pose_rescoring",
}
RESULT_TYPES = PRIMARY_RESULT_TYPES | FOLLOWUP_RESULT_TYPES
_SAFE_COMPONENT = re.compile(r"^[A-Za-z0-9_.-]{1,180}$")


@dataclass
class LegacyRun:
    run_id: str
    task_group: str
    path: Path
    metadata: dict[str, Any]
    inputs: dict[str, Any]
    result: dict[str, Any]
    references: set[str] = field(default_factory=set)


@dataclass
class PlannedRun:
    run_id: str
    task_group: str
    source_path: str
    role: str
    project_names: list[str]
    result_type: str = ""
    tool: str = ""
    status: str = "unknown"
    pdb_id: str = ""
    receptor: str = ""
    compound_count: int | None = None
    replicates: int | None = None
    campaign_name: str = ""
    bytes: int = 0


@dataclass
class MigrationPlan:
    source_root: str
    datasets: list[dict[str, str]]
    projects: dict[str, list[dict[str, Any]]]
    runs: list[PlannedRun]
    targets: list[dict[str, Any]]
    warnings: list[str]
    source_file_count: int
    total_bytes: int

    def to_json(self) -> dict[str, Any]:
        return {
            "schema_version": 1,
            "source_root": self.source_root,
            "datasets": self.datasets,
            "projects": self.projects,
            "runs": [asdict(run) for run in self.runs],
            "targets": self.targets,
            "warnings": self.warnings,
            "source_file_count": self.source_file_count,
            "total_bytes": self.total_bytes,
        }


def _read_json(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError):
        return {}
    return payload if isinstance(payload, dict) else {}


def _references(value: Any) -> set[str]:
    """Collect run IDs from old and new nested artifact-reference formats."""
    found: set[str] = set()
    if isinstance(value, dict):
        for key, item in value.items():
            # These are descriptive links in an edit log, not run inputs.
            if key in {"modification_history", "legacy_preparation_evidence"}:
                continue
            if key.endswith("run_id") and isinstance(item, str) and item:
                found.add(item)
            elif key.endswith("run_ids") and isinstance(item, list):
                found.update(entry for entry in item if isinstance(entry, str) and entry)
            found.update(_references(item))
    elif isinstance(value, list):
        for item in value:
            found.update(_references(item))
    return found


def scan_runs(source_root: Path | str) -> tuple[dict[str, LegacyRun], list[str]]:
    root = Path(source_root).expanduser().resolve()
    if not root.is_dir():
        raise FileNotFoundError(f"MN-Ligand run directory not found: {root}")
    runs: dict[str, LegacyRun] = {}
    warnings: list[str] = []
    for metadata_path in sorted(root.glob("*/*/metadata.json")):
        run_dir = metadata_path.parent.resolve()
        try:
            run_dir.relative_to(root)
        except ValueError:
            warnings.append(f"Skipped run outside the selected source root: {metadata_path}")
            continue
        metadata = _read_json(metadata_path)
        if not metadata:
            warnings.append(f"Skipped unreadable metadata: {metadata_path}")
            continue
        run_id = str(metadata.get("run_id") or run_dir.name)
        task_group = str(metadata.get("task_group") or run_dir.parent.name)
        if not _SAFE_COMPONENT.fullmatch(run_id) or not _SAFE_COMPONENT.fullmatch(task_group):
            warnings.append(f"Skipped unsafe run path components: {task_path(run_dir)}")
            continue
        inputs = _read_json(run_dir / "input.json")
        run = LegacyRun(
            run_id=run_id,
            task_group=task_group,
            path=run_dir,
            metadata=metadata,
            inputs=inputs,
            # Large Boltz/AlphaFold result files can contain thousands of
            # embedded structures. The planner never needs to deserialize
            # them; target import reports are read lazily below.
            result={},
        )
        run.references = (_references(metadata) | _references(inputs)) - {run_id}
        if run_id in runs:
            warnings.append(
                f"Duplicate run ID {run_id!r} occurs in {runs[run_id].path} and {run_dir}; references are ambiguous."
            )
            continue
        runs[run_id] = run
    return runs, warnings


def task_path(path: Path) -> str:
    return f"{path.parent.name}/{path.name}"


def _descendants(seed_id: str, children: dict[str, set[str]]) -> set[str]:
    visited = {seed_id}
    todo = deque([seed_id])
    while todo:
        parent_id = todo.popleft()
        for child_id in children.get(parent_id, ()):
            if child_id not in visited:
                visited.add(child_id)
                todo.append(child_id)
    return visited


def _tool(run: LegacyRun) -> str:
    return str(run.metadata.get("tool") or run.metadata.get("engine") or run.metadata.get("workflow") or run.metadata.get("job_type") or "unknown")


def _target_source(run: LegacyRun) -> str:
    result = run.result or _read_json(run.path / "result.json")
    for key in ("imported_target", "prepared_target", "target"):
        relative = result.get(key)
        if isinstance(relative, str) and relative:
            candidate = (run.path / relative).resolve()
            try:
                candidate.relative_to(run.path)
            except ValueError:
                continue
            if candidate.is_file():
                return str(candidate)
    for candidate in (
        run.path / "artifacts" / "imported" / f"{str(run.metadata.get('pdb_id') or '').lower()}.pdb",
        run.path / "artifacts" / "prepared" / "prepared_target.pdb",
    ):
        if candidate.is_file():
            return str(candidate.resolve())
    return ""


def _tree_stats(paths: Iterable[Path]) -> tuple[int, int]:
    count = size = 0
    for root in paths:
        for current, dirs, files in os.walk(root, followlinks=False):
            current_path = Path(current)
            dirs[:] = [name for name in dirs if not (current_path / name).is_symlink()]
            for name in files:
                path = current_path / name
                if path.is_symlink():
                    continue
                try:
                    stat = path.stat()
                except OSError:
                    continue
                count += 1
                size += stat.st_size
    return count, size


def build_migration_plan(
    source_root: Path | str = DEFAULT_SOURCE_ROOT,
    *,
    dataset_projects: Iterable[DatasetProject] = DATASET_PROJECTS,
) -> MigrationPlan:
    """Build a dataset-aware, target-aware copy plan without changing files."""
    root = Path(source_root).expanduser().resolve()
    runs, warnings = scan_runs(root)
    ancestor_cache: dict[str, frozenset[str]] = {}

    def ancestors(run_id: str, visiting: frozenset[str] = frozenset()) -> frozenset[str]:
        cached = ancestor_cache.get(run_id)
        if cached is not None:
            return cached
        if run_id in visiting or run_id not in runs:
            return frozenset()
        next_visiting = visiting | {run_id}
        values: set[str] = set()
        for parent_id in runs[run_id].references:
            if parent_id not in runs:
                continue
            values.add(parent_id)
            values.update(ancestors(parent_id, next_visiting))
        frozen = frozenset(values)
        ancestor_cache[run_id] = frozen
        return frozen

    target_cache: dict[str, frozenset[str]] = {}

    def target_pdb_ids(run: LegacyRun) -> frozenset[str]:
        cached = target_cache.get(run.run_id)
        if cached is not None:
            return cached
        pdb_ids = set()
        own_pdb = str(run.metadata.get("pdb_id") or "").upper()
        if own_pdb in THR_TARGETS:
            pdb_ids.add(own_pdb)
        for ancestor_id in ancestors(run.run_id):
            pdb_id = str(runs[ancestor_id].metadata.get("pdb_id") or "").upper()
            if pdb_id in THR_TARGETS:
                pdb_ids.add(pdb_id)
        frozen = frozenset(pdb_ids)
        target_cache[run.run_id] = frozen
        return frozen

    children: dict[str, set[str]] = defaultdict(set)
    for run in runs.values():
        for parent_id in run.references:
            if parent_id in runs:
                children[parent_id].add(run.run_id)

    plans: list[PlannedRun] = []
    projects: dict[str, list[dict[str, Any]]] = {}
    datasets: list[dict[str, str]] = []
    targets_by_pdb: dict[str, dict[str, Any]] = {}

    for dataset in dataset_projects:
        imports = [
            run for run in runs.values()
            if run.metadata.get("job_type") == "compound_import"
            and str(run.metadata.get("dataset_name") or "").casefold() == dataset.source_name.casefold()
        ]
        if len(imports) != 1:
            warnings.append(
                f"Expected one compound import for {dataset.source_name!r}; found {len(imports)}."
            )
            if not imports:
                continue
        if len(imports) > 1:
            warnings.append(f"Using the most recent matching compound import for {dataset.source_name!r}.")
            imports.sort(key=lambda item: str(item.metadata.get("created_at") or ""))
        source_import = imports[-1]
        descendants = _descendants(source_import.run_id, children)
        candidates = [runs[run_id] for run_id in descendants if run_id in runs and runs[run_id].metadata.get("job_type") in RESULT_TYPES]

        primary_ids = {
            run.run_id for run in candidates
            if run.metadata.get("job_type") in PRIMARY_RESULT_TYPES
            and len(target_pdb_ids(run)) == 1
        }
        selected: dict[str, tuple[LegacyRun, str, set[str]]] = {}
        for run in candidates:
            job_type = str(run.metadata.get("job_type") or "")
            pdb_ids = target_pdb_ids(run)
            if len(pdb_ids) != 1:
                continue
            if job_type in FOLLOWUP_RESULT_TYPES:
                run_ancestors = ancestors(run.run_id)
                if not run_ancestors.intersection(primary_ids):
                    continue
            elif job_type not in PRIMARY_RESULT_TYPES:
                continue
            pdb_id = next(iter(pdb_ids))
            selected[run.run_id] = (run, pdb_id, set(ancestors(run.run_id)))

        # Every selected run keeps its complete recorded prerequisite chain.
        dependency_ids: set[str] = set()
        for run_id, (_, _, ancestor_ids) in selected.items():
            dependency_ids.update(ancestor_ids)
            dependency_ids.add(run_id)
        missing_refs: set[str] = set()
        for run_id in dependency_ids:
            missing_refs.update(runs[run_id].references - runs.keys())
        if missing_refs:
            warnings.append(
                f"{dataset.short_name}: {len(missing_refs)} referenced prerequisite run(s) are absent from the source root."
            )

        project_jobs: list[dict[str, Any]] = []
        for run_id, (run, pdb_id, ancestor_ids) in selected.items():
            target = THR_TARGETS[pdb_id]
            tool = _tool(run)
            compound_count = run.metadata.get("compound_count")
            replicates = run.metadata.get("replicates")
            label = f"{target['receptor']} · {pdb_id} · {dataset.short_name} · {tool} · {run_id[:8]}"
            plans.append(PlannedRun(
                run_id=run_id,
                task_group=run.task_group,
                source_path=str(run.path),
                role="result",
                project_names=[dataset.project_name],
                result_type=str(run.metadata.get("job_type") or ""),
                tool=tool,
                status=str(run.metadata.get("status") or "unknown"),
                pdb_id=pdb_id,
                receptor=target["receptor"],
                compound_count=compound_count if isinstance(compound_count, int) else None,
                replicates=replicates if isinstance(replicates, int) else None,
                campaign_name=label,
            ))
            project_jobs.append({
                "run_id": run_id,
                "task_group": run.task_group,
                "name": label,
                "pdb_id": pdb_id,
                "receptor": target["receptor"],
                "tool": tool,
                "job_type": run.metadata.get("job_type"),
                "status": run.metadata.get("status") or "unknown",
                "compound_count": compound_count,
                "replicates": replicates,
                "created_at": run.metadata.get("created_at"),
                "dependency_run_ids": sorted(ancestor_ids),
            })
            targets_by_pdb.setdefault(pdb_id, {
                "pdb_id": pdb_id,
                "receptor": target["receptor"],
                "name": target["name"],
                "source_run_id": "",
                "source_file": "",
            })
        projects[dataset.project_name] = project_jobs
        aggregate_group_count = len({(str(job.get("pdb_id") or "").upper(), str(job.get("tool") or "")) for job in project_jobs})
        datasets.append({
            "source_name": dataset.source_name,
            "short_name": dataset.short_name,
            "project_name": dataset.project_name,
            "source_run_id": source_import.run_id,
            "source_compound_count": str(source_import.metadata.get("compound_count") or ""),
            "result_job_count": str(aggregate_group_count),
            "source_result_run_count": str(len(project_jobs)),
            "prerequisite_run_count": str(len(dependency_ids)),
        })

    # Shared prerequisites are copied once; result runs are tagged to exactly
    # one dataset project. Unassigned dependencies retain their old IDs.
    result_ids = {plan.run_id for plan in plans}
    dependency_to_projects: dict[str, set[str]] = defaultdict(set)
    # Determine full dependency ownership from the plan's ancestry lists.
    for project_name, jobs in projects.items():
        for job in jobs:
            for dependency_id in job["dependency_run_ids"]:
                dependency_to_projects[dependency_id].add(project_name)

    dependency_ids = set(dependency_to_projects) - result_ids
    for run_id in sorted(dependency_ids):
        run = runs.get(run_id)
        if run is None:
            continue
        projects_for_dependency = sorted(dependency_to_projects[run_id])
        plans.append(PlannedRun(
            run_id=run.run_id,
            task_group=run.task_group,
            source_path=str(run.path),
            role="dependency",
            project_names=projects_for_dependency,
            result_type=str(run.metadata.get("job_type") or ""),
            tool=_tool(run),
            status=str(run.metadata.get("status") or "unknown"),
            pdb_id=str(run.metadata.get("pdb_id") or ""),
            receptor=THR_TARGETS.get(str(run.metadata.get("pdb_id") or "").upper(), {}).get("receptor", ""),
            compound_count=run.metadata.get("compound_count") if isinstance(run.metadata.get("compound_count"), int) else None,
            replicates=run.metadata.get("replicates") if isinstance(run.metadata.get("replicates"), int) else None,
        ))

    # Resolve clean source proteins from their MN-Ligand import runs. Keep a
    # separate protein entry per PDB structure, even for the same subtype.
    for pdb_id, target in targets_by_pdb.items():
        choices = [
            run for run in runs.values()
            if str(run.metadata.get("pdb_id") or "").upper() == pdb_id
            and run.metadata.get("job_type") == "protein_import"
            and str(run.metadata.get("status") or "").casefold() == "completed"
        ]
        choices.sort(key=lambda item: str(item.metadata.get("created_at") or ""))
        for run in choices:
            source_file = _target_source(run)
            if source_file:
                target["source_run_id"] = run.run_id
                target["source_file"] = source_file
                break
        if not target["source_file"]:
            warnings.append(f"No completed protein-import structure found for {pdb_id} ({target['receptor']}).")

    planned_paths = [Path(plan.source_path) for plan in plans]
    file_count, total_bytes = _tree_stats(planned_paths)
    plans.sort(key=lambda item: (item.role != "result", item.project_names, item.task_group, item.run_id))
    return MigrationPlan(
        source_root=str(root),
        datasets=datasets,
        projects=projects,
        runs=plans,
        targets=[targets_by_pdb[key] for key in sorted(targets_by_pdb)],
        warnings=warnings,
        source_file_count=file_count,
        total_bytes=total_bytes,
    )


def format_plan_summary(plan: MigrationPlan) -> str:
    lines = [f"MN-Ligand migration plan from {plan.source_root}"]
    for dataset in plan.datasets:
        jobs = plan.projects.get(dataset["project_name"], [])
        statuses: dict[str, int] = defaultdict(int)
        tools: dict[str, set[tuple[str, str]]] = defaultdict(set)
        for job in jobs:
            statuses[str(job["status"])] += 1
            tools[str(job["tool"])].add((str(job.get("pdb_id") or ""), str(job.get("tool") or "")))
        aggregate_jobs = len({(str(job.get("pdb_id") or "").upper(), str(job.get("tool") or "")) for job in jobs})
        lines.append(
            f"- {dataset['short_name']} → {dataset['project_name']}: {aggregate_jobs} grouped jobs from "
            f"{len(jobs)} source runs, "
            f"{dataset['prerequisite_run_count']} run folders in dependency closure"
        )
        if tools:
            lines.append("  Tools: " + ", ".join(f"{key} ({len(tools[key])})" for key in sorted(tools)))
        if statuses:
            lines.append("  Statuses: " + ", ".join(f"{key} ({statuses[key]})" for key in sorted(statuses)))
    if plan.targets:
        lines.append("Targets: " + ", ".join(f"{item['receptor']} {item['pdb_id']}" for item in plan.targets))
    lines.append(
        f"Copy footprint: {len(plan.runs)} unique run folders, {plan.source_file_count} files, "
        f"{plan.total_bytes / (1024 ** 3):.2f} GiB"
    )
    if plan.warnings:
        lines.append("Warnings:")
        lines.extend(f"- {warning}" for warning in plan.warnings)
    return "\n".join(lines)


def _write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _dataset_compound_lookup(database, groups: dict[str, dict[str, Any]]) -> dict[str, dict[str, dict[str, str]]]:
    """Index imported dataset identities by their source compound IDs."""
    result: dict[str, dict[str, dict[str, str]]] = {}
    with database.connection() as connection:
        for project_name, group in groups.items():
            rows = connection.execute(
                "SELECT name,metadata,components FROM molecule_members WHERE group_id=?",
                (str(group["id"]),),
            ).fetchall()
            compounds: dict[str, dict[str, str]] = {}
            for row in rows:
                try:
                    metadata = json.loads(row["metadata"] or "{}")
                    components = json.loads(row["components"] or "[]")
                except (TypeError, ValueError):
                    continue
                compound_id = str(metadata.get("compound_id") or "").strip()
                smiles = next((str(item.get("smiles") or "").strip() for item in components
                               if isinstance(item, dict) and item.get("type") == "small_molecule"), "")
                if compound_id and smiles:
                    compounds[compound_id] = {
                        "compound_id": compound_id,
                        "compound_name": str(row["name"] or compound_id),
                        "smiles": smiles,
                    }
            result[project_name] = compounds
    return result


def _read_csv_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    try:
        with path.open(newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle)
            return list(reader.fieldnames or []), [dict(row) for row in reader]
    except (OSError, UnicodeError, csv.Error):
        return [], []


def _write_bytes_if_changed(path: Path, content: bytes) -> bool:
    if path.is_file() and path.read_bytes() == content:
        return False
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_bytes(content)
    temporary.replace(path)
    return True


def _write_csv_if_changed(
    path: Path,
    columns: list[str],
    rows: list[dict[str, Any]],
    *,
    delimiter: str = ",",
) -> bool:
    from io import StringIO

    buffer = StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=columns, extrasaction="ignore", delimiter=delimiter, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return _write_bytes_if_changed(path, buffer.getvalue().encode("utf-8"))


def _float_or_none(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if number == number and abs(number) != float("inf") else None


def _metric_means(rows: list[dict[str, str]], excluded: set[str]) -> dict[str, float]:
    values: dict[str, list[float]] = defaultdict(list)
    for row in rows:
        for key, value in row.items():
            if key in excluded:
                continue
            number = _float_or_none(value)
            if number is not None:
                values[key].append(number)
    return {key: sum(numbers) / len(numbers) for key, numbers in values.items() if numbers}


def _safe_run_relative_file(run_dir: Path, value: str) -> str:
    if not value:
        return ""
    path = (run_dir / value).resolve()
    try:
        return path.relative_to(run_dir.resolve()).as_posix() if path.is_file() else ""
    except ValueError:
        return ""


def _legacy_prediction_candidates(
    run_dir: Path,
    *,
    csv_path: Path,
    kind: str,
    compounds: dict[str, dict[str, str]],
) -> tuple[list[dict[str, Any]], int]:
    """Convert MN-Ligand's Boltz-2/AF3 score tables into Studio candidates.

    One record is emitted per compound and replicate. For engines that return
    multiple models per replicate, metrics are averaged and the strongest
    confidence/ranking model supplies the attached complex structure. The
    original detailed score table remains in place beside this normalized
    manifest.
    """
    _columns, rows = _read_csv_rows(csv_path)
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    unresolved = 0
    for row in rows:
        compound_id = str(row.get("compound_id") or row.get("candidate_id") or "").strip()
        if not compound_id:
            continue
        if compound_id not in compounds:
            unresolved += 1
            continue
        replicate = str(row.get("replicate") or row.get("model_seed") or row.get("seed") or "1").strip()
        grouped[(compound_id, replicate)].append(row)

    seed_values = sorted({replicate for _compound_id, replicate in grouped}, key=lambda value: (not value.isdigit(), int(value) if value.isdigit() else value))
    seed_to_replicate = {seed: index + 1 for index, seed in enumerate(seed_values)}
    candidates: list[dict[str, Any]] = []
    for (compound_id, replicate), replicate_rows in sorted(grouped.items()):
        if kind == "boltz2":
            structure_key = "structure_file"
            rank_key = "confidence_score"
            excluded = {"candidate_id", "compound_id", "model_id", "replicate", "seed", "structure_file"}
        else:
            structure_key = ""
            rank_key = "ranking_score"
            excluded = {"candidate_id", "compound_id", "model_seed", "prediction_id", "sample"}

        def rank(row: dict[str, str]) -> float:
            value = _float_or_none(row.get(rank_key))
            return value if value is not None else float("-inf")

        representative = max(replicate_rows, key=rank)
        if kind == "boltz2":
            relative_structure = _safe_run_relative_file(run_dir, representative.get(structure_key) or "")
        else:
            prediction_id = str(representative.get("prediction_id") or "").strip()
            expected = f"output/{compound_id}/{prediction_id}/{compound_id}_{prediction_id}_model.cif" if prediction_id else ""
            relative_structure = _safe_run_relative_file(run_dir, expected)
        compound = compounds[compound_id]
        metrics = _metric_means(replicate_rows, excluded)
        metrics.update({
            "replicate": float(seed_to_replicate[replicate]),
            "seed": _float_or_none(representative.get("seed") or representative.get("model_seed")) or 0.0,
            "model_count": float(len(replicate_rows)),
        })
        raw_metadata = {
            "source_compound_id": compound_id,
            "compound_id": compound_id,
            "compound_name": compound["compound_name"],
            "replicate": replicate,
            "model_count": len(replicate_rows),
            "ligand_smiles": [compound["smiles"]],
            "components": [{"type": "small_molecule", "chain": "L", "smiles": compound["smiles"]}],
        }
        if kind == "boltz2":
            raw_metadata["representative_model_id"] = representative.get("model_id") or ""
        else:
            raw_metadata["prediction_id"] = representative.get("prediction_id") or ""
            raw_metadata["model_seed"] = representative.get("model_seed") or ""
            raw_metadata["sample_count"] = len(replicate_rows)
        candidate: dict[str, Any] = {
            "candidate_id": f"{compound_id}::legacy-replicate-{seed_to_replicate[replicate]:03d}",
            "stage": "complex_refolding",
            "source_tool": "Boltz-2" if kind == "boltz2" else "AlphaFold 3",
            "metrics": metrics,
            "raw_metadata": raw_metadata,
        }
        if relative_structure:
            candidate["complex_pdb"] = relative_structure
        candidates.append(candidate)
    return candidates, unresolved


def normalize_mn_ligand_run(
    run_dir: Path | str,
    *,
    project_name: str,
    compounds: dict[str, dict[str, str]],
) -> dict[str, Any]:
    """Normalize result artifacts from one safely copied MN-Ligand run."""
    run_dir = Path(run_dir).expanduser().resolve()
    migration = _read_json(run_dir / "mn_ligand_migration.json")
    if migration.get("source_application") != "MN-Ligand" or migration.get("role") != "result":
        return {"run_id": run_dir.name, "normalized_files": [], "unresolved_compound_ids": 0, "skipped": True}

    normalized_files: list[str] = []
    unresolved_ids: set[str] = set()
    score_rows: list[dict[str, str]] = []
    for source in (run_dir / "scores.csv", run_dir / "artifacts" / "scores.csv"):
        if not source.is_file():
            continue
        columns, score_rows = _read_csv_rows(source)
        if not columns or not any("compound_id" in row for row in score_rows):
            continue
        output_columns = list(columns)
        if "pose_path" not in output_columns:
            output_columns.append("pose_path")
        normalized_rows = []
        for source_row in score_rows:
            row = dict(source_row)
            if not row.get("pose_path") and row.get("pose_file"):
                row["pose_path"] = row["pose_file"]
            normalized_rows.append(row)
        destination = run_dir / "artifacts" / "docking_scores.csv"
        if _write_csv_if_changed(destination, output_columns, normalized_rows):
            normalized_files.append(str(destination.relative_to(run_dir)))
        break

    affinity_path = run_dir / "artifacts" / "nesso_affinity.csv"
    if affinity_path.is_file():
        columns, rows = _read_csv_rows(affinity_path)
        if columns and "compound_id" not in columns and "candidate_id" in columns:
            backup = run_dir / "artifacts" / "mn_ligand_source" / "nesso_affinity.csv"
            if not backup.exists():
                backup.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(affinity_path, backup)
            output_columns = ["compound_id", *columns]
            normalized_rows = [{"compound_id": row.get("candidate_id") or "", **row} for row in rows]
            if _write_csv_if_changed(affinity_path, output_columns, normalized_rows):
                normalized_files.append(str(affinity_path.relative_to(run_dir)))
        if "compound_id" in columns or "candidate_id" in columns:
            for row in rows:
                compound_id = str(row.get("compound_id") or row.get("candidate_id") or "").strip()
                if compound_id:
                    score_rows.append({"compound_id": compound_id})

    # Make input identities explicit for legacy affinity tables. This allows
    # Studio to attach measurements to the already imported molecule records.
    compound_ids = {
        str(row.get("compound_id") or row.get("candidate_id") or "").strip()
        for row in score_rows if str(row.get("compound_id") or row.get("candidate_id") or "").strip()
    }
    candidates_path = run_dir / "artifacts" / "normalized_candidates" / "candidates.jsonl"
    candidate_count = 0
    candidate_csv = next((path for path, kind in (
        (run_dir / "artifacts" / "boltz2_metrics.csv", "boltz2"),
        (run_dir / "artifacts" / "alphafold3_metrics.csv", "af3"),
    ) if path.is_file()), None)
    if candidate_csv and not candidates_path.exists():
        kind = "boltz2" if candidate_csv.name == "boltz2_metrics.csv" else "af3"
        candidates, unresolved = _legacy_prediction_candidates(
            run_dir, csv_path=candidate_csv, kind=kind, compounds=compounds
        )
        unresolved_ids.update(
            str(row.get("compound_id") or row.get("candidate_id") or "").strip()
            for _columns, rows in [_read_csv_rows(candidate_csv)] for row in rows
            if str(row.get("compound_id") or row.get("candidate_id") or "").strip() not in compounds
        )
        candidate_count = len(candidates)
        manifest = "".join(json.dumps(candidate, sort_keys=True, ensure_ascii=False) + "\n" for candidate in candidates)
        if _write_bytes_if_changed(candidates_path, manifest.encode("utf-8")):
            normalized_files.append(str(candidates_path.relative_to(run_dir)))
        for candidate in candidates:
            source_id = str((candidate.get("raw_metadata") or {}).get("source_compound_id") or "").strip()
            if source_id:
                compound_ids.add(source_id)

    input_compounds_path = run_dir / "artifacts" / "input" / "compounds.tsv"
    adapter_report = _read_json(run_dir / "legacy_adapter_report.json")
    adapter_owned_input = "artifacts/input/compounds.tsv" in set(adapter_report.get("normalized_files") or [])
    if compound_ids and (not input_compounds_path.is_file() or adapter_owned_input) and not (run_dir / "input" / "compounds.tsv").is_file():
        known = [compounds[compound_id] for compound_id in sorted(compound_ids) if compound_id in compounds]
        unresolved_ids.update(compound_id for compound_id in compound_ids if compound_id not in compounds)
        if known:
            from io import StringIO

            buffer = StringIO(newline="")
            writer = csv.DictWriter(buffer, fieldnames=["compound_id", "compound_name", "smiles"], delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(known)
            if _write_bytes_if_changed(input_compounds_path, buffer.getvalue().encode("utf-8")):
                normalized_files.append(str(input_compounds_path.relative_to(run_dir)))

    if normalized_files:
        existing = _read_json(run_dir / "legacy_adapter_report.json")
        report = {
            "schema_version": 1,
            "source_application": "MN-Ligand",
            "project_name": project_name,
            "normalized_files": sorted(set(existing.get("normalized_files", [])) | set(normalized_files)),
            "cofold_candidate_count": candidate_count or int(existing.get("cofold_candidate_count") or 0),
            "unresolved_compound_ids": sorted(unresolved_ids),
        }
        _write_json_atomic(run_dir / "legacy_adapter_report.json", report)
    return {
        "run_id": run_dir.name,
        "normalized_files": sorted(set(normalized_files)),
        "cofold_candidate_count": candidate_count,
        "unresolved_compound_ids": len(unresolved_ids),
        "skipped": False,
    }


def _aggregate_descriptor(run_dir: Path, marker: dict[str, Any], metadata: dict[str, Any]) -> dict[str, str]:
    campaign = str(metadata.get("campaign_name") or "")
    parts = [part.strip() for part in campaign.split(" · ") if part.strip()]
    project_names = marker.get("project_names") or []
    embedded_migration = metadata.get("legacy_migration")
    embedded_migration = embedded_migration if isinstance(embedded_migration, dict) else {}
    project_name = str(project_names[0] if project_names else embedded_migration.get("project_name") or "")
    tool = str(metadata.get("tool") or metadata.get("engine") or "unknown")
    pdb_id = str(metadata.get("pdb_id") or "").upper()
    if len(parts) >= 4:
        pdb_id = pdb_id or (parts[1].upper() if re.fullmatch(r"[0-9][A-Z0-9]{3}", parts[1], re.I) else "")
    if not project_name or not tool or not pdb_id:
        raise ValueError(f"Cannot group migrated run without project, target PDB, and tool metadata: {run_dir}")
    target_label = parts[0] if len(parts) >= 4 else pdb_id
    dataset_label = parts[2] if len(parts) >= 4 else project_name
    label = " · ".join((target_label, pdb_id, dataset_label, tool))
    return {
        "project_name": project_name,
        "pdb_id": pdb_id,
        "tool": tool,
        "campaign_name": label,
        "task_group": run_dir.parent.name,
    }


def _aggregate_csv(
    children: list[dict[str, Any]],
    destination: Path,
    filename: str,
    *,
    pose_paths: bool = False,
) -> bool:
    columns: list[str] = []
    merged: list[dict[str, Any]] = []
    for child in children:
        source = child["run_dir"] / "artifacts" / filename
        if not source.is_file():
            continue
        source_columns, rows = _read_csv_rows(source)
        for column in source_columns:
            if column not in columns:
                columns.append(column)
        for source_row in rows:
            row: dict[str, Any] = dict(source_row)
            if pose_paths:
                raw_path = str(row.get("pose_path") or row.get("pose_file") or "").strip()
                if raw_path:
                    candidate = Path(raw_path)
                    relative = ""
                    if candidate.is_absolute():
                        source_run = Path(str(child["marker"].get("source_path") or ""))
                        try:
                            relative = candidate.relative_to(source_run).as_posix()
                        except ValueError:
                            pass
                    else:
                        relative = raw_path
                    if relative and (child["run_dir"] / relative).is_file():
                        row["pose_path"] = (Path("artifacts") / "replicates" / child["run_id"] / relative).as_posix()
            row["legacy_run_id"] = child["run_id"]
            row["attempt_index"] = child["attempt_index"]
            row["attempt_status"] = child["metadata"].get("status") or "unknown"
            merged.append(row)
    if not merged:
        return False
    for column in ("legacy_run_id", "attempt_index", "attempt_status"):
        if column not in columns:
            columns.append(column)
    if pose_paths and "pose_path" not in columns:
        columns.append("pose_path")
    return _write_csv_if_changed(destination / "artifacts" / filename, columns, merged)


def _collect_aggregate_compounds(replicate_dirs: Iterable[Path]) -> dict[str, dict[str, str]]:
    """Read compound identities from either MN-Ligand input-table location."""
    compound_rows: dict[str, dict[str, str]] = {}
    for run_dir in replicate_dirs:
        for compounds_path in (
            run_dir / "artifacts" / "input" / "compounds.tsv",
            run_dir / "input" / "compounds.tsv",
        ):
            if not compounds_path.is_file():
                continue
            try:
                with compounds_path.open(newline="", encoding="utf-8") as handle:
                    reader = csv.DictReader(handle, delimiter="\t")
                    for row in reader:
                        compound_id = str(row.get("compound_id") or "").strip()
                        if compound_id:
                            compound_rows.setdefault(compound_id, dict(row))
            except (OSError, UnicodeError, csv.Error):
                continue
    return compound_rows


def _write_aggregate_compounds(job_dir: Path, compound_rows: dict[str, dict[str, str]]) -> bool:
    if not compound_rows:
        return False
    return _write_csv_if_changed(
        job_dir / "artifacts" / "input" / "compounds.tsv",
        ["compound_id", "compound_name", "smiles"],
        [compound_rows[key] for key in sorted(compound_rows)],
        delimiter="\t",
    )


def aggregate_migrated_result_jobs(
    results_root: Path | str,
    *,
    target_workstation: str | None = None,
) -> list[Path]:
    """Bundle source attempts into one visible Studio job per project/target/tool.

    Original per-attempt directories remain intact beneath each aggregate
    job's ``artifacts/replicates`` folder. Prerequisite runs are retained under
    an ignored ``_mn_ligand_migration`` directory so they do not appear as
    additional jobs in Studio's Results list.
    """
    root = Path(results_root).expanduser().resolve()
    index_path = root / "_mn_ligand_migration" / "aggregate_index.json"
    index = _read_json(index_path)
    indexed_jobs = index.get("jobs") if isinstance(index.get("jobs"), list) else []
    source_to_job = index.get("source_run_to_job") if isinstance(index.get("source_run_to_job"), dict) else {}

    raw_results: list[Path] = []
    dependencies: list[Path] = []
    for group_dir in root.iterdir():
        if not group_dir.is_dir() or group_dir.is_symlink() or group_dir.name.startswith("_"):
            continue
        for run_dir in group_dir.iterdir():
            if not run_dir.is_dir() or run_dir.is_symlink():
                continue
            marker = _read_json(run_dir / "mn_ligand_migration.json")
            if marker.get("source_application") != "MN-Ligand":
                continue
            if marker.get("role") == "result":
                raw_results.append(run_dir)
            elif marker.get("role") == "dependency":
                dependencies.append(run_dir)

    if index_path.is_file() and raw_results:
        raise ValueError("This Results directory has aggregate jobs plus unbundled MN-Ligand runs; refusing to mix them.")
    if index_path.is_file() and not raw_results:
        jobs = []
        for item in indexed_jobs:
            job_path = (root / str(item.get("path") or "")).resolve()
            if job_path.is_dir():
                replicate_root = job_path / "artifacts" / "replicates"
                replicate_dirs = [path for path in sorted(replicate_root.iterdir()) if path.is_dir()]
                _write_aggregate_compounds(
                    job_path,
                    _collect_aggregate_compounds(replicate_dirs),
                )
                jobs.append(job_path)
        return jobs

    grouped: dict[tuple[str, str, str], list[dict[str, Any]]] = defaultdict(list)
    for run_dir in sorted(raw_results):
        marker = _read_json(run_dir / "mn_ligand_migration.json")
        metadata = _read_json(run_dir / "metadata.json")
        descriptor = _aggregate_descriptor(run_dir, marker, metadata)
        key = (descriptor["project_name"], descriptor["pdb_id"], descriptor["tool"])
        grouped[key].append({"path": run_dir, "marker": marker, "metadata": metadata, "descriptor": descriptor})

    aggregate_jobs: list[dict[str, Any]] = []
    source_run_to_job: dict[str, dict[str, str]] = {}
    for key, entries in sorted(grouped.items()):
        entries.sort(key=lambda item: (str(item["metadata"].get("created_at") or ""), item["path"].name))
        project_name, pdb_id, tool = key
        task_groups = {item["path"].parent.name for item in entries}
        task_group = sorted(task_groups)[0]
        aggregate_digest = hashlib.sha256("\0".join(key).encode()).hexdigest()[:16]
        aggregate_id = f"mnligand-{aggregate_digest}"
        final_dir = root / task_group / aggregate_id
        source_ids = [item["path"].name for item in entries]
        if final_dir.exists():
            existing = _read_json(final_dir / "mn_ligand_aggregate.json")
            if existing.get("source_run_ids") != source_ids:
                raise FileExistsError(f"Aggregate destination already contains different data: {final_dir}")
            aggregate_jobs.append({"job_id": aggregate_id, "task_group": task_group,
                                   "path": final_dir.relative_to(root).as_posix(), "source_run_ids": source_ids,
                                   "project_name": project_name, "pdb_id": pdb_id, "tool": tool,
                                   "campaign_name": entries[0]["descriptor"]["campaign_name"]})
            for source_id in source_ids:
                source_run_to_job[source_id] = {"job_id": aggregate_id, "path": final_dir.relative_to(root).as_posix()}
            continue

        stage_parent = Path(tempfile.mkdtemp(prefix=f".{aggregate_id}.stage-", dir=final_dir.parent))
        stage_dir = stage_parent / aggregate_id
        stage_dir.mkdir()
        replicate_root = stage_dir / "artifacts" / "replicates"
        replicate_root.mkdir(parents=True)
        children: list[dict[str, Any]] = []
        moved: list[tuple[Path, Path]] = []
        try:
            for attempt_index, entry in enumerate(entries, start=1):
                source = entry["path"]
                replicate_dir = replicate_root / source.name
                os.replace(source, replicate_dir)
                moved.append((replicate_dir, source))
                children.append({
                    "run_id": source.name,
                    "attempt_index": attempt_index,
                    "run_dir": replicate_dir,
                    "marker": entry["marker"],
                    "metadata": entry["metadata"],
                    "input": _read_json(replicate_dir / "input.json"),
                    "result": _read_json(replicate_dir / "result.json"),
                })

            status_counts: dict[str, int] = defaultdict(int)
            for child in children:
                status_counts[str(child["metadata"].get("status") or "unknown")] += 1
            first = children[0]
            first_metadata = dict(first["metadata"])
            first_marker = _read_json(stage_dir / "artifacts" / "replicates" / first["run_id"] / "molecular_design_studio.json")
            owner_id = str(first_marker.get("studio_owner_id") or first_metadata.get("studio_owner_id") or "")
            project_id = str(first_marker.get("project_id") or first_metadata.get("project_id") or "")
            campaign_name = entries[0]["descriptor"]["campaign_name"]

            compound_rows: dict[str, dict[str, str]] = {}
            candidates: list[dict[str, Any]] = []
            for child in children:
                candidate_file = child["run_dir"] / "artifacts" / "normalized_candidates" / "candidates.jsonl"
                if candidate_file.is_file():
                    with candidate_file.open(encoding="utf-8") as handle:
                        for line in handle:
                            try:
                                candidate = json.loads(line)
                            except json.JSONDecodeError:
                                continue
                            if not isinstance(candidate, dict):
                                continue
                            candidate["candidate_id"] = f"{child['run_id']}::{candidate.get('candidate_id') or 'candidate'}"
                            candidate["raw_metadata"] = dict(candidate.get("raw_metadata") or {})
                            candidate["raw_metadata"]["legacy_run_id"] = child["run_id"]
                            candidate["metrics"] = dict(candidate.get("metrics") or {})
                            candidate["metrics"]["legacy_attempt_index"] = child["attempt_index"]
                            for field_name in ("complex_pdb", "binder_pdb", "target_pdb"):
                                value = str(candidate.get(field_name) or "")
                                if value and not Path(value).is_absolute():
                                    candidate[field_name] = (Path("artifacts") / "replicates" / child["run_id"] / value).as_posix()
                            candidates.append(candidate)

            compound_rows = _collect_aggregate_compounds(child["run_dir"] for child in children)
            _write_aggregate_compounds(stage_dir, compound_rows)
            if candidates:
                candidates_path = stage_dir / "artifacts" / "normalized_candidates" / "candidates.jsonl"
                _write_bytes_if_changed(
                    candidates_path,
                    "".join(json.dumps(candidate, sort_keys=True, ensure_ascii=False) + "\n" for candidate in candidates).encode("utf-8"),
                )
            _aggregate_csv(children, stage_dir, "docking_scores.csv", pose_paths=True)
            _aggregate_csv(children, stage_dir, "nesso_affinity.csv")

            replicate_manifest = {
                "schema_version": 1,
                "aggregate_run_id": aggregate_id,
                "source_run_count": len(children),
                "status_counts": dict(sorted(status_counts.items())),
                "attempts": [
                    {
                        "run_id": child["run_id"],
                        "attempt_index": child["attempt_index"],
                        "status": child["metadata"].get("status") or "unknown",
                        "created_at": child["metadata"].get("created_at"),
                        "completed_at": child["metadata"].get("completed_at"),
                        "source_path": child["marker"].get("source_path"),
                    }
                    for child in children
                ],
            }
            _write_json_atomic(replicate_root / "manifest.json", replicate_manifest)

            metadata = first_metadata
            metadata.update({
                "run_id": aggregate_id,
                "task_group": task_group,
                "tool": tool,
                "pdb_id": pdb_id,
                "campaign_name": campaign_name,
                "project_id": project_id,
                "studio_owner_id": owner_id,
                "status": "completed",
                "created_at": min((str(child["metadata"].get("created_at") or "") for child in children if child["metadata"].get("created_at")), default=""),
                "completed_at": max((str(child["metadata"].get("completed_at") or "") for child in children if child["metadata"].get("completed_at")), default=""),
                "legacy_migration": {
                    "source_application": "MN-Ligand",
                    "role": "aggregate_job",
                    "project_name": project_name,
                    "pdb_id": pdb_id,
                    "tool": tool,
                    "source_run_ids": source_ids,
                    "source_status_counts": dict(sorted(status_counts.items())),
                    "target_workstation": str(target_workstation or ""),
                },
            })
            inputs = {
                "tool": tool,
                "target_pdb": pdb_id,
                "legacy_aggregate": {
                    "project_name": project_name,
                    "source_run_ids": source_ids,
                    "replicate_manifest": "artifacts/replicates/manifest.json",
                },
            }
            result = {
                "success": True,
                "metrics": {"source_run_count": len(children), "completed_run_count": status_counts.get("completed", 0)},
                "legacy_aggregate": {"status_counts": dict(sorted(status_counts.items()))},
            }
            _write_json_atomic(stage_dir / "metadata.json", metadata)
            _write_json_atomic(stage_dir / "input.json", inputs)
            _write_json_atomic(stage_dir / "result.json", result)
            _write_json_atomic(stage_dir / "mn_ligand_aggregate.json", {
                "schema_version": 1,
                "source_run_ids": source_ids,
                "project_name": project_name,
                "pdb_id": pdb_id,
                "tool": tool,
                "campaign_name": campaign_name,
            })
            _write_json_atomic(stage_dir / "mn_ligand_migration.json", {
                "schema_version": 1,
                "source_application": "MN-Ligand",
                "role": "aggregate",
                "project_names": [project_name],
                "source_run_ids": source_ids,
                "target_workstation": str(target_workstation or ""),
            })
            _write_json_atomic(stage_dir / "molecular_design_studio.json", {
                "application": "mn-molecular-design-studio",
                "workspace": first_marker.get("workspace") or "",
                "marked_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                "project_id": project_id,
                "studio_owner_id": owner_id,
                "migration_source": "MN-Ligand",
                "workstation_name": str(target_workstation or first_marker.get("workstation_name") or ""),
            })
            final_dir.parent.mkdir(parents=True, exist_ok=True)
            os.replace(stage_dir, final_dir)
            aggregate_job = {
                "job_id": aggregate_id,
                "task_group": task_group,
                "path": final_dir.relative_to(root).as_posix(),
                "source_run_ids": source_ids,
                "project_name": project_name,
                "pdb_id": pdb_id,
                "tool": tool,
                "campaign_name": campaign_name,
            }
            aggregate_jobs.append(aggregate_job)
            for source_id in source_ids:
                source_run_to_job[source_id] = {"job_id": aggregate_id, "path": aggregate_job["path"]}
        except Exception:
            for current, original in reversed(moved):
                if current.exists() and not original.exists():
                    original.parent.mkdir(parents=True, exist_ok=True)
                    os.replace(current, original)
            raise
        finally:
            shutil.rmtree(stage_parent, ignore_errors=True)

    dependency_rows = []
    hidden_dependency_root = root / "_mn_ligand_migration" / "dependencies"
    for dependency in dependencies:
        marker = _read_json(dependency / "mn_ligand_migration.json")
        run_id = str(marker.get("source_run_id") or dependency.name)
        destination = hidden_dependency_root / dependency.parent.name / run_id
        destination.parent.mkdir(parents=True, exist_ok=True)
        if dependency.exists():
            if destination.exists():
                raise FileExistsError(f"Migration dependency destination already exists: {destination}")
            os.replace(dependency, destination)
        dependency_rows.append({"run_id": run_id, "task_group": dependency.parent.name,
                                "path": destination.relative_to(root).as_posix()})

    index_payload = {
        "schema_version": 1,
        "source_run_count": len(source_run_to_job),
        "aggregate_job_count": len(aggregate_jobs),
        "source_run_to_job": source_run_to_job,
        "jobs": aggregate_jobs,
        "dependencies": dependency_rows,
    }
    _write_json_atomic(index_path, index_payload)
    return [root / job["path"] for job in aggregate_jobs]


def _name_key(value: str) -> str:
    return "".join(character for character in value.casefold() if character.isalnum())


def apply_migration_plan(
    plan: MigrationPlan,
    *,
    user: dict[str, Any],
    results_root: Path | str,
    workspace_home: Path | str | None = None,
    database=None,
    target_workstation: str | None = None,
    progress: Callable[[dict[str, Any]], None] | None = None,
) -> dict[str, Any]:
    """Copy planned data into Studio and register project jobs/targets.

    The source stays untouched. Run folders are staged beside their final
    destination and atomically renamed after metadata and ownership markers
    have been written. Calling this repeatedly skips folders already marked as
    imported from the same source run.
    """
    from mn_molecular_design_studio.auth import grant_project_access
    from mn_molecular_design_studio.molecular_database import MolecularDatabase
    from mn_molecular_design_studio.workspace import create_project, load_organization

    owner_id = str(user.get("id") or "")
    if not owner_id:
        raise ValueError("An authenticated Studio user is required to apply the migration.")
    destination_root = Path(results_root).expanduser().resolve()
    destination_root.mkdir(parents=True, exist_ok=True)
    source_plan_root = Path(plan.source_root).expanduser().resolve()
    if destination_root == source_plan_root:
        raise ValueError("The destination Results directory cannot be the MN-Ligand source directory.")
    required_bytes = int(plan.total_bytes * 1.05)
    if shutil.disk_usage(destination_root).free < required_bytes:
        raise OSError(
            f"Not enough free space for the planned copy: need about {required_bytes / (1024 ** 3):.1f} GiB."
        )
    home = Path(workspace_home).expanduser().resolve() if workspace_home else None
    database = database or MolecularDatabase()

    # Resolve source groups and destination conflicts before changing projects
    # or group memberships.
    accessible_groups = database.groups(user, scope="all")
    group_index: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for group in accessible_groups:
        group_index[_name_key(str(group.get("name") or ""))].append(group)
    source_groups: dict[str, dict[str, Any]] = {}
    for dataset in DATASET_PROJECTS:
        matches: dict[str, dict[str, Any]] = {}
        for name in dataset.group_names:
            for group in group_index.get(_name_key(name), []):
                matches[str(group["id"])] = group
        if len(matches) != 1:
            raise ValueError(
                f"Expected one existing Studio molecule group for {dataset.short_name}; found {len(matches)}."
            )
        source_groups[dataset.project_name] = next(iter(matches.values()))

    source_root = Path(plan.source_root).resolve()
    migration_index = _read_json(destination_root / "_mn_ligand_migration" / "aggregate_index.json")
    previously_aggregated = migration_index.get("source_run_to_job")
    previously_aggregated = previously_aggregated if isinstance(previously_aggregated, dict) else {}
    for run in plan.runs:
        source = Path(run.source_path).resolve()
        source.relative_to(source_root)
        if not source.is_dir():
            raise FileNotFoundError(f"Planned source run disappeared: {source}")
        if run.role == "result" and run.run_id in previously_aggregated:
            aggregate_path = destination_root / str(previously_aggregated[run.run_id].get("path") or "")
            if not aggregate_path.is_dir():
                raise FileNotFoundError(f"Previously imported aggregate job disappeared: {aggregate_path}")
            continue
        if run.role == "dependency":
            hidden_dependency = destination_root / "_mn_ligand_migration" / "dependencies" / run.task_group / run.run_id
            if hidden_dependency.is_dir():
                continue
        destination = destination_root / run.task_group / run.run_id
        provenance_path = destination / "mn_ligand_migration.json"
        if destination.is_dir() and provenance_path.is_file():
            provenance = _read_json(provenance_path)
            if provenance.get("source_run_id") == run.run_id and provenance.get("source_root") == str(source_root):
                continue
            raise FileExistsError(f"Destination already contains different data: {destination}")
        if destination.exists():
            raise FileExistsError(f"Destination already exists without matching migration provenance: {destination}")

    # Resolve or create dataset projects, then make them visible to the caller.
    projects = load_organization(home).get("projects", [])
    project_ids: dict[str, str] = {}
    project_owner_ids = {
        project_name: str(group.get("owner_id") or owner_id)
        for project_name, group in source_groups.items()
    }
    for dataset in DATASET_PROJECTS:
        existing = next(
            (project for project in projects
             if isinstance(project, dict)
             and str(project.get("name") or "").strip().casefold() == dataset.project_name.casefold()),
            None,
        )
        if existing:
            project = existing
        else:
            project = create_project(
                dataset.project_name,
                description=f"MN-Ligand results migrated for {dataset.source_name}.",
                home=home,
                owner_id=project_owner_ids[dataset.project_name],
            )
            projects.append(project)
        project_id = str(project["id"])
        project_owner_id = project_owner_ids[dataset.project_name]
        if user.get("role") == "admin" and project_owner_id != owner_id:
            grant_project_access(project_id, project_owner_id, actor=user)
        elif user.get("role") != "admin":
            if project_owner_id != owner_id:
                raise PermissionError("Only an administrator can migrate a group owned by another user.")
            grant_project_access(project_id, owner_id, actor=user, allow_self_grant=True)
        project_ids[dataset.project_name] = project_id

    # Move the already imported dataset groups under their matching projects,
    # preserving the existing Studio molecule IDs and group contents.
    group_moves: list[dict[str, str]] = []
    for dataset in DATASET_PROJECTS:
        group = source_groups[dataset.project_name]
        project_id = project_ids[dataset.project_name]
        if str(group.get("project_id") or "") != project_id:
            database.set_group_project(user, group["id"], project_id)
        group_moves.append({"group_id": str(group["id"]), "name": str(group["name"]), "project_id": project_id})

    # Verify each source folder still exists and belongs under the planned root.
    planned = {run.run_id: run for run in plan.runs}
    copied = skipped = 0
    copied_results: list[Path] = []
    for index, run in enumerate(plan.runs, start=1):
        source = Path(run.source_path).resolve()
        try:
            source.relative_to(source_root)
        except ValueError as exc:
            raise ValueError(f"Planned run escapes the source root: {source}") from exc
        if not source.is_dir():
            raise FileNotFoundError(f"Planned source run disappeared: {source}")
        destination_group = destination_root / run.task_group
        destination = destination_group / run.run_id
        if run.role == "result" and run.run_id in previously_aggregated:
            skipped += 1
            if progress:
                progress({"phase": "copy", "completed": index, "total": len(plan.runs), "run_id": run.run_id, "status": "grouped"})
            continue
        if run.role == "dependency":
            hidden_dependency = destination_root / "_mn_ligand_migration" / "dependencies" / run.task_group / run.run_id
            if hidden_dependency.is_dir():
                skipped += 1
                if progress:
                    progress({"phase": "copy", "completed": index, "total": len(plan.runs), "run_id": run.run_id, "status": "present"})
                continue
        destination_group.mkdir(parents=True, exist_ok=True)
        provenance_path = destination / "mn_ligand_migration.json"
        if destination.is_dir() and provenance_path.is_file():
            provenance = _read_json(provenance_path)
            if provenance.get("source_run_id") == run.run_id and provenance.get("source_root") == str(source_root):
                skipped += 1
                if run.role == "result":
                    copied_results.append(destination)
                if progress:
                    progress({"phase": "copy", "completed": index, "total": len(plan.runs), "run_id": run.run_id, "status": "present"})
                continue
            raise FileExistsError(f"Destination already contains different data: {destination}")
        if destination.exists():
            raise FileExistsError(f"Destination already exists without matching migration provenance: {destination}")

        staging_parent = Path(tempfile.mkdtemp(prefix=f".{run.run_id}.import-", dir=destination_group))
        staging = staging_parent / run.run_id
        try:
            # Preserve symlinks as links. Some MN-Ligand containers leave
            # dangling absolute links to Docker-only mount points; following
            # them on the host would escape the source run and cannot recover
            # the container data anyway.
            shutil.copytree(source, staging, symlinks=True)
            copied_metadata = _read_json(staging / "metadata.json")
            original_campaign_name = str(copied_metadata.get("campaign_name") or "")
            if run.role == "result":
                project_name = run.project_names[0]
                project_id = project_ids[project_name]
                project_owner_id = project_owner_ids[project_name]
                copied_metadata["campaign_name"] = run.campaign_name
                copied_metadata["project_id"] = project_id
                copied_metadata["studio_owner_id"] = project_owner_id
            copied_metadata["legacy_migration"] = {
                "source_application": "MN-Ligand",
                "source_run_id": run.run_id,
                "source_task_group": run.task_group,
                "source_path": str(source),
                "role": run.role,
                "project_names": run.project_names,
                "source_campaign_name": original_campaign_name,
                "symlinks_preserved": True,
            }
            if run.role == "result":
                copied_metadata["legacy_migration"].update({
                    "pdb_id": run.pdb_id,
                    "receptor": run.receptor,
                    "tool": run.tool,
                    "source_run_id": run.run_id,
                })
            if target_workstation:
                copied_metadata["legacy_migration"]["target_workstation"] = str(target_workstation)
            _write_json_atomic(staging / "metadata.json", copied_metadata)
            _write_json_atomic(staging / "mn_ligand_migration.json", {
                "source_application": "MN-Ligand",
                "source_run_id": run.run_id,
                "source_root": str(source_root),
                "source_path": str(source),
                "migrated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                "role": run.role,
                "project_names": run.project_names,
                "symlinks_preserved": True,
            })
            if run.role == "result":
                _write_json_atomic(staging / "molecular_design_studio.json", {
                    "application": "mn-molecular-design-studio",
                    "workspace": str(home) if home else "",
                    "marked_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                    "project_id": project_ids[run.project_names[0]],
                    "studio_owner_id": project_owner_ids[run.project_names[0]],
                    "migration_source": "MN-Ligand",
                    "workstation_name": str(target_workstation or ""),
                })
            os.replace(staging, destination)
            copied += 1
            if run.role == "result":
                copied_results.append(destination)
            if progress:
                progress({"phase": "copy", "completed": index, "total": len(plan.runs), "run_id": run.run_id, "status": "copied"})
        finally:
            shutil.rmtree(staging_parent, ignore_errors=True)

    # Normalize the old score and prediction formats in the copied runs before
    # indexing them. The original MN-Ligand tree remains untouched.
    compound_lookup = _dataset_compound_lookup(database, source_groups)
    run_by_id = {run.run_id: run for run in plan.runs if run.role == "result"}
    normalization_rows = []
    for run_dir in copied_results:
        run = run_by_id.get(run_dir.name)
        if run is None or not run.project_names:
            continue
        normalization_rows.append(normalize_mn_ligand_run(
            run_dir,
            project_name=run.project_names[0],
            compounds=compound_lookup.get(run.project_names[0], {}),
        ))
        if target_workstation:
            metadata_path = run_dir / "metadata.json"
            copied_metadata = _read_json(metadata_path)
            copied_metadata.setdefault("legacy_migration", {})["target_workstation"] = str(target_workstation)
            _write_json_atomic(metadata_path, copied_metadata)
            marker_path = run_dir / "molecular_design_studio.json"
            marker = _read_json(marker_path)
            marker["workstation_name"] = str(target_workstation)
            _write_json_atomic(marker_path, marker)

    # Persist source protein coordinates in Studio's targets folder and in the
    # authenticated user's reusable structure storage.
    target_registry: list[dict[str, Any]] = []
    if home is None:
        from mn_molecular_design_studio.workspace import _studio_home

        home = _studio_home()
    from mn_molecular_design_studio.target_editor import _target_storage_path

    targets_directory = home / "workdir" / "targets"
    targets_directory.mkdir(parents=True, exist_ok=True)
    for target in plan.targets:
        source = Path(str(target.get("source_file") or ""))
        if not source.is_file():
            raise FileNotFoundError(f"Target structure is unavailable for {target.get('pdb_id')}: {source}")
        content = source.read_bytes()
        digest = hashlib.sha256(content).hexdigest()
        target_filename = f"{target['receptor'].replace('α', 'alpha').replace('β', 'beta')}-{target['pdb_id']}.pdb"
        stored_path = _target_storage_path(content, target_filename, owner_id=owner_id)
        durable_path = targets_directory / target_filename
        if durable_path.exists() and _sha256(durable_path) != digest:
            raise FileExistsError(f"A different target structure already exists at {durable_path}")
        if not durable_path.exists():
            durable_path.write_bytes(content)
        target_registry.append({
            "pdb_id": target["pdb_id"],
            "receptor": target["receptor"],
            "name": target["name"],
            "path": str(durable_path),
            "user_structure_path": str(stored_path),
            "sha256": digest,
            "source_run_id": target.get("source_run_id"),
        })
    registry_path = targets_directory / "mn_ligand_targets.json"
    existing_registry = _read_json(registry_path)
    registry_rows = existing_registry.get("targets", []) if isinstance(existing_registry.get("targets"), list) else []
    registry_by_pdb = {str(item.get("pdb_id")): item for item in registry_rows if isinstance(item, dict)}
    for target in target_registry:
        registry_by_pdb[target["pdb_id"]] = target
    _write_json_atomic(registry_path, {"schema_version": 1, "targets": list(registry_by_pdb.values())})

    # Group the legacy attempts into one Results entry per dataset, target and
    # tool while retaining every source attempt beneath the aggregate job.
    aggregate_jobs = aggregate_migrated_result_jobs(
        destination_root,
        target_workstation=target_workstation,
    )

    # Index the aggregate jobs, then ingest their merged score and prediction
    # artifacts so each attempt's best structures link to Studio compounds.
    from mn_molecular_design_studio.job_catalog import reconcile
    from mn_molecular_design_studio.molecular_database_ingest import record_completed_run

    reconcile(destination_root)
    current_workstation = ""
    try:
        from mn_molecular_design_studio.cluster import config as cluster_config

        current_workstation = cluster_config.display_name(cluster_config.load().get("name") or "")
        if target_workstation:
            target_workstation = cluster_config.display_name(target_workstation)
    except (ImportError, OSError, ValueError):
        pass
    deferred_ingestion = bool(target_workstation and current_workstation and
                              target_workstation.casefold() != current_workstation.casefold())
    if not deferred_ingestion:
        for index, run_dir in enumerate(aggregate_jobs, start=1):
            record_completed_run(run_dir)
            if progress:
                progress({"phase": "register", "completed": index, "total": len(aggregate_jobs), "run_id": run_dir.name, "status": "registered"})
    return {
        "copied_run_count": copied,
        "already_present_run_count": skipped,
        "source_result_run_count": sum(1 for run in planned.values() if run.role == "result"),
        "result_job_count": len(aggregate_jobs),
        "dataset_group_assignments": group_moves,
        "projects": project_ids,
        "targets": target_registry,
        "normalized_run_count": sum(bool(row.get("normalized_files")) for row in normalization_rows),
        "normalized_artifact_count": sum(len(row.get("normalized_files") or []) for row in normalization_rows),
        "unresolved_compound_id_count": sum(int(row.get("unresolved_compound_ids") or 0) for row in normalization_rows),
        "ingestion_deferred_to_workstation": target_workstation if deferred_ingestion else "",
        "results_root": str(destination_root),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, default=DEFAULT_SOURCE_ROOT)
    parser.add_argument("--manifest", type=Path, help="Write the read-only plan as JSON.")
    parser.add_argument("--smoke", action="store_true", help="Validate required datasets and TRα/TRβ coverage.")
    args = parser.parse_args(argv)
    try:
        plan = build_migration_plan(args.source_root)
    except (FileNotFoundError, OSError, ValueError) as exc:
        print(f"Migration planning failed: {exc}", file=sys.stderr)
        return 2
    print(format_plan_summary(plan))
    if args.manifest:
        args.manifest.parent.mkdir(parents=True, exist_ok=True)
        args.manifest.write_text(json.dumps(plan.to_json(), indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print(f"Manifest: {args.manifest.resolve()}")
    if args.smoke:
        present_datasets = {item["source_name"] for item in plan.datasets}
        expected = {item.source_name for item in DATASET_PROJECTS}
        receptors = {target["receptor"] for target in plan.targets if target.get("source_file")}
        errors = []
        if present_datasets != expected:
            errors.append(f"Expected all 3 source datasets; found {len(present_datasets)}.")
        if not {"TRα", "TRβ"}.issubset(receptors):
            errors.append(f"Expected source structures for TRα and TRβ; found {sorted(receptors)}.")
        if not any(job["tool"] == "vina" for jobs in plan.projects.values() for job in jobs):
            errors.append("No Vina jobs were recognized.")
        if not any(job["tool"] == "Nesso-1" for jobs in plan.projects.values() for job in jobs):
            errors.append("No Nesso-1 affinity jobs were recognized.")
        if not any(job["tool"] in {"Boltz-2", "AlphaFold 3"} for jobs in plan.projects.values() for job in jobs):
            errors.append("No cofolding/refolding jobs were recognized.")
        if errors:
            print("Smoke check FAILED:\n- " + "\n- ".join(errors), file=sys.stderr)
            return 1
        print("Smoke check passed: dataset lineage, project mapping, targets, and representative tools recognized.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
