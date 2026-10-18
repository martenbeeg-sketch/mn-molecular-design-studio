from __future__ import annotations

import ast
import csv
import glob
from pathlib import Path
import re
import shlex
import subprocess
from typing import Literal

from mn_molecular_design_studio.engine.workflow_engine.core.artifacts import Artifact, artifact_path
from mn_molecular_design_studio.engine.workflow_engine.core.gpu import docker_gpu_args, gpu_queue_resource, normalize_gpu_device
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import JobPaths, collect_jobs, create_job, finish_job, read_json, update_status, write_json
from mn_molecular_design_studio.engine.workflow_engine.core.manifests import load_manifest
from mn_molecular_design_studio.engine.workflow_engine.core.scheduler import apply_docker_cpu_limits_to_steps
from mn_molecular_design_studio.engine.workflow_engine.core.structures import STANDARD_PROTEIN_RESIDUES, filter_pdb_text, pdb_summary, residue_inventory, sanitize_pdb_for_surface_tools
from mn_molecular_design_studio.engine.workflow_engine.runtime import PROJECT_DIR, app_home


DETECTION_GROUP = "detection"
HOTSPOT_GROUP = "hotspot-detection"
KNOWN_BENCHMARK_ROOT = Path("/mnt/db/reference_files/de_novo_binder_scoring_overath_2025")
KNOWN_BENCHMARK_CSV = KNOWN_BENCHMARK_ROOT / "final_dataset.csv"
KNOWN_BENCHMARK_PDB_DIR = KNOWN_BENCHMARK_ROOT / "input_pdbs"


def _parse_chain_list(value: object) -> list[str]:
    if isinstance(value, (list, tuple, set)):
        chains: list[str] = []
        for item in value:
            if isinstance(item, dict):
                chain_id = str(item.get("chain_id") or item.get("chain") or "").strip()
            else:
                chain_id = str(item).strip()
            if chain_id:
                chains.append(chain_id)
        return chains
    text = str(value or "").strip()
    if not text:
        return []
    try:
        parsed = ast.literal_eval(text)
    except (SyntaxError, ValueError):
        parsed = None
    if isinstance(parsed, (list, tuple, set)):
        return _parse_chain_list(parsed)
    return [item.strip() for item in text.split(",") if item.strip()]


def _resolve_run_path(run_dir: Path, value: object) -> Path | None:
    text = str(value or "").strip()
    if not text:
        return None
    path = Path(text).expanduser()
    if not path.is_absolute():
        path = run_dir / path
    try:
        return path.resolve()
    except OSError:
        return path


def _target_key(path: Path, chains: list[str]) -> tuple[str, tuple[str, ...]]:
    try:
        path_text = str(path.resolve())
    except OSError:
        path_text = str(path)
    return path_text, tuple(chains)


def _safe_target_token(value: object, fallback: str = "target") -> str:
    token = "".join(ch.lower() if ch.isalnum() else "_" for ch in str(value or fallback))
    token = "_".join(part for part in token.split("_") if part)
    return token or fallback


def _normalize_target_category(value: object) -> str:
    category = str(value or "").strip() or "prepared"
    if category == "masked":
        return "mutated"
    return category


def _add_target_row(rows: list[dict], seen: set[tuple[str, tuple[str, ...]]], row: dict) -> None:
    target_pdb = Path(str(row.get("target_pdb") or "")).expanduser()
    if not target_pdb.exists():
        return
    chains = _parse_chain_list(row.get("chains") or [])
    if not chains:
        try:
            chains = chains_for_pdb(target_pdb)
        except Exception:
            chains = []
    key = _target_key(target_pdb, chains)
    existing = next((item for item in rows if _target_key(Path(str(item.get("target_pdb") or "")), item.get("chains") or []) == key), None)
    if existing is not None:
        existing["records"] = int(existing.get("records") or 1) + int(row.get("records") or 1)
        return
    if key in seen:
        return
    seen.add(key)
    row["target_pdb"] = str(target_pdb)
    row["chains"] = chains
    row.setdefault("records", 1)
    row["prepared_kind"] = _normalize_target_category(row.get("prepared_kind") or row.get("source_category") or "prepared")
    row["source_category"] = _normalize_target_category(row.get("source_category") or row.get("prepared_kind") or "prepared")
    rows.append(row)


def _prepared_target_rows(rows: list[dict], seen: set[tuple[str, tuple[str, ...]]]) -> None:
    for row in collect_jobs("target-prep"):
        if row.get("status") != "completed":
            continue
        run_dir = Path(row["run_dir"])
        result = read_json(run_dir / "result.json")
        target = (result.get("outputs") or {}).get("target") or {}
        downstream = result.get("downstream_artifacts", {})
        target_pdb = run_dir / downstream.get("target_trimmed_pdb", "artifacts/target_trimmed.pdb")
        prepared_kind = _normalize_target_category(target.get("prepared_kind") or "trimmed")
        if not target_pdb.exists():
            target_pdb = run_dir / downstream.get("target_clean_pdb", "artifacts/target_clean.pdb")
            prepared_kind = _normalize_target_category(target.get("prepared_kind") or "cleaned")
        source_label = target.get("source_label")
        if not source_label:
            source_label = (
                "Target mutation"
                if prepared_kind == "mutated"
                else "Target fragment split"
                if prepared_kind == "split_fragments"
                else "Target preparation"
            )
        _add_target_row(
            rows,
            seen,
            {
                **row,
                "target_name": target.get("target_name") or row["job_code"],
                "target_pdb": str(target_pdb),
                "chains": target.get("chains") or [],
                "residue_count": target.get("residue_count"),
                "source_category": prepared_kind,
                "prepared_kind": prepared_kind,
                "source_label": source_label,
                "target_entities": target.get("target_entities") or [],
                "fragment_split": target.get("fragment_split") or {},
            },
        )

    for row in collect_jobs("target-crop"):
        if row.get("status") != "completed":
            continue
        run_dir = Path(row["run_dir"])
        result = read_json(run_dir / "result.json")
        crop = (result.get("outputs") or {}).get("crop") or {}
        target_pdb = run_dir / result.get("downstream_artifacts", {}).get("target_cropped_pdb", "artifacts/target_cropped.pdb")
        if not target_pdb.exists():
            continue
        summary = pdb_summary(target_pdb.read_text(errors="ignore"))
        chains = [chain["chain_id"] for chain in summary["chains"]]
        source_name = crop.get("target_name") or (result.get("inputs") or {}).get("target_name") or row["job_code"]
        _add_target_row(
            rows,
            seen,
            {
                **row,
                "target_name": f"{source_name} crop",
                "target_pdb": str(target_pdb),
                "chains": chains,
                "residue_count": sum(chain["residue_count"] for chain in summary["chains"]),
                "source_category": "cropped",
                "prepared_kind": "cropped",
                "source_label": "Target crop",
            },
        )


def _imported_target_rows(rows: list[dict], seen: set[tuple[str, tuple[str, ...]]]) -> None:
    for row in collect_jobs("target-prep"):
        if row.get("status") != "completed":
            continue
        run_dir = Path(row["run_dir"])
        result = read_json(run_dir / "result.json")
        inputs = result.get("inputs") if isinstance(result.get("inputs"), dict) else {}
        target = (result.get("outputs") or {}).get("target") or {}
        prepared_kind = _normalize_target_category(target.get("prepared_kind") or "")
        if row.get("tool") == "target_fragment_splitter" or prepared_kind == "split_fragments":
            continue
        source_pdb = Path(str(inputs.get("source_pdb") or target.get("source_pdb") or "")).expanduser()
        if not source_pdb.exists():
            target_input = run_dir / "artifacts" / "target_input.pdb"
            source_pdb = target_input if target_input.exists() else source_pdb
        if not source_pdb.exists():
            continue
        _add_target_row(
            rows,
            seen,
            {
                **row,
                "target_name": target.get("target_name") or inputs.get("target_name") or source_pdb.stem,
                "target_pdb": str(source_pdb),
                "chains": [],
                "source_category": "imported",
                "prepared_kind": "imported",
                "source_label": "Imported target protein",
            },
        )


def _benchmark_target_rows(rows: list[dict], seen: set[tuple[str, tuple[str, ...]]]) -> None:
    if not KNOWN_BENCHMARK_CSV.exists() or not KNOWN_BENCHMARK_PDB_DIR.exists():
        return
    grouped: dict[tuple[str, str], dict] = {}
    try:
        with KNOWN_BENCHMARK_CSV.open(newline="") as handle:
            for record in csv.DictReader(handle):
                target_id = str(record.get("target_id") or "").strip()
                source_name = str(record.get("source") or "").strip()
                binder_id = str(record.get("binder_id") or "").strip()
                if not target_id or not binder_id:
                    continue
                key = (target_id, source_name)
                grouped.setdefault(key, {"record": record, "records": 0})
                grouped[key]["records"] += 1
    except Exception:
        return
    for (target_id, source_name), payload in sorted(grouped.items()):
        record = payload["record"]
        binder_id = str(record.get("binder_id") or "").strip()
        source_complex_pdb = KNOWN_BENCHMARK_PDB_DIR / f"{binder_id}.pdb"
        chains = _parse_chain_list(record.get("target_chains"))
        if not source_complex_pdb.exists() or not chains:
            continue
        target_cache_dir = app_home() / "workdir" / "targets" / "benchmark_target_only"
        target_cache_dir.mkdir(parents=True, exist_ok=True)
        target_file = (
            f"{_safe_target_token(target_id)}__"
            f"{_safe_target_token(source_name, 'source')}__"
            f"{'_'.join(_safe_target_token(chain, 'chain') for chain in chains)}.pdb"
        )
        target_pdb = target_cache_dir / target_file
        if not target_pdb.exists():
            target_text = filter_pdb_text(
                source_complex_pdb.read_text(errors="ignore"),
                keep_chains=set(chains),
                remove_waters=True,
                remove_hetero=True,
            )
            if "ATOM" not in target_text:
                continue
            target_pdb.write_text(target_text)
        _add_target_row(
            rows,
            seen,
            {
                "run_id": f"overath-{target_id}-{source_name}",
                "job_code": "BENCH",
                "task_group": "benchmark-reference",
                "status": "available",
                "created_at": "",
                "updated_at": "",
                "tool": "Overath 2025",
                "target_name": target_id,
                "target_pdb": str(target_pdb),
                "chains": chains,
                "source_category": "benchmark",
                "prepared_kind": "benchmark",
                "source_label": source_name or "Overath 2025 benchmark",
                "source_complex_pdb": str(source_complex_pdb),
                "records": int(payload["records"]),
            },
        )


def completed_target_jobs() -> list[dict]:
    rows: list[dict] = []
    seen: set[tuple[str, tuple[str, ...]]] = set()
    _prepared_target_rows(rows, seen)

    rows.sort(key=lambda row: row.get("updated_at") or row.get("created_at") or "", reverse=True)
    return rows


def ppi_target_jobs() -> list[dict]:
    rows: list[dict] = []
    seen: set[tuple[str, tuple[str, ...]]] = set()
    _prepared_target_rows(rows, seen)
    _imported_target_rows(rows, seen)
    _benchmark_target_rows(rows, seen)

    rows.sort(key=lambda row: row.get("updated_at") or row.get("created_at") or "", reverse=True)
    return rows


def detection_jobs_for_target(target_pdb: Path) -> list[dict]:
    target_path = str(target_pdb)
    target_resolved = str(target_pdb.expanduser().resolve()) if target_pdb.exists() else target_path
    rows: list[dict] = []
    for group in [DETECTION_GROUP, HOTSPOT_GROUP]:
        for row in collect_jobs(group):
            run_dir = Path(row["run_dir"])
            payload = read_json(run_dir / "input.json")
            inputs = payload.get("inputs") or {}
            input_target = str(inputs.get("target_pdb") or "")
            input_path = Path(input_target).expanduser()
            input_resolved = str(input_path.resolve()) if input_path.exists() else input_target
            if input_target != target_path and input_resolved != target_resolved:
                continue
            result = read_json(run_dir / "result.json")
            rows.append(
                {
                    **row,
                    "task_group": group,
                    "success": result.get("success"),
                    "artifact_count": (result.get("metrics") or {}).get("artifact_count", ""),
                }
            )
    return rows


def target_label(row: dict) -> str:
    chains = ",".join(row.get("chains") or [])
    category = str(row.get("source_category") or row.get("prepared_kind") or "target")
    job_code = str(row.get("job_code") or "")
    records = int(row.get("records") or 0)
    record_text = f", {records} records" if records > 1 else ""
    return f"{row['target_name']} [{category}] ({job_code}, chains {chains or 'unknown'}{record_text})"


def _run_shell_steps(run_dir: Path, steps: list[dict], *, append: bool = False) -> int:
    steps = apply_docker_cpu_limits_to_steps(run_dir, steps)
    previous_steps = read_json(run_dir / "command.json").get("steps", []) if append else []
    write_json(run_dir / "command.json", {"mode": "docker", "steps": [*previous_steps, *steps]})
    update_status(run_dir, "running")
    log_mode = "a" if append else "w"
    with (run_dir / "stdout.log").open(log_mode) as stdout, (run_dir / "stderr.log").open(log_mode) as stderr:
        for step in steps:
            stdout.write(f"$ {' '.join(step['command'])}\n")
            stdout.flush()
            stderr.write(f"\n--- {step['name']} ---\n")
            stderr.flush()
            proc = subprocess.run(step["command"], stdout=stdout, stderr=stderr, check=False)
            if proc.returncode != 0:
                return int(proc.returncode)
    return 0


def _docker_base(image: str, run_dir: Path, gpu: bool = False, gpu_device: object = "0") -> list[str]:
    command = ["docker", "run", "--rm"]
    if gpu:
        command.extend(docker_gpu_args(gpu_device))
    command.extend(["-v", f"{run_dir}:/work"])
    return command + [image]


def _collect_artifacts(run_dir: Path, root: Path, patterns: list[tuple[str, str]]) -> list[dict]:
    artifacts: list[dict] = []
    for pattern, artifact_type in patterns:
        for path_text in glob.glob(str(root / pattern), recursive=True):
            path = Path(path_text)
            if path.is_file() and path.stat().st_size:
                artifacts.append(Artifact(path.stem, path, artifact_type).to_json(run_dir))
    return artifacts


def _finish_detection_tool(job: JobPaths, tool: str, success: bool, result: dict, *, defer_finish: bool) -> None:
    """Finish an independent detection job or record its result inside a bundle."""
    status = {"success": bool(success), **result}
    if not success and not status.get("error"):
        metrics = result.get("metrics") if isinstance(result.get("metrics"), dict) else {}
        return_code = metrics.get("return_code")
        status["error"] = (
            f"The method exited with return code {return_code}."
            if return_code not in (None, 0)
            else "The method completed without producing its required prediction files."
        )
    if not defer_finish:
        finish_job(job.run_dir, success, status)
        return
    status_path = job.run_dir / "artifacts" / "surface_analysis" / "tool_status.json"
    statuses = read_json(status_path) if status_path.is_file() else {}
    statuses[tool] = status
    write_json(status_path, statuses)


def validate_surface_target(target_pdb: Path, chain_ids: list[str] | None = None) -> list[str]:
    """Reject unusable target selections before a container job is queued."""
    if target_pdb.suffix.casefold() not in {".pdb", ".ent"}:
        raise ValueError("Protein surface tools require a PDB (.pdb/.ent) structure.")
    text = target_pdb.read_text(errors="replace")
    available = {
        row["chain"] for row in residue_inventory(text)
        if row["resname"] in STANDARD_PROTEIN_RESIDUES
    }
    chains = list(dict.fromkeys(str(chain).strip() for chain in chain_ids)) if chain_ids is not None else sorted(available)
    if not chains or not set(chains).issubset(available):
        raise ValueError("Select existing protein chains from the target PDB structure.")
    if any(not re.fullmatch(r"[A-Za-z0-9]", chain) for chain in chains):
        raise ValueError("Surface tools require single letter or digit chain IDs. Remap blank or special chain IDs first.")
    return chains


def _surf2spot_numbered_input(pdb_text: str) -> str:
    """Give Surf2Spot contiguous IDs while keeping the original structure separately."""
    residues = [row for row in residue_inventory(pdb_text) if row["resname"] in STANDARD_PROTEIN_RESIDUES]
    counters: dict[str, int] = {}
    numbering = {}
    for row in residues:
        chain = row["chain"]
        counters[chain] = counters.get(chain, 0) + 1
        if counters[chain] > 9999:
            raise ValueError("Surf2Spot accepts at most 9999 residues per PDB chain.")
        numbering[(chain, row["residue_number"], row["insertion_code"])] = counters[chain]
    atoms: dict[tuple, list[str]] = {}
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM  ") or len(line) < 27 or line[17:20].strip() not in STANDARD_PROTEIN_RESIDUES:
            continue
        key = (line[21].strip(), int(line[22:26]), line[26].strip())
        number = numbering[key]
        atoms.setdefault(key, []).append(f"{line[:22]}{number:4d} {line[27:]}")
    return "\n".join([*(line for key in sorted(atoms) for line in atoms[key]), "END", ""])


def _write_detection_worker_request(job: JobPaths) -> None:
    write_json(
        job.run_dir / "worker_request.json",
        {
            "kind": "detection_tool",
            "path_kwargs": ["target_pdb"],
        },
    )
    write_json(
        job.run_dir / "command.json",
        {
            "mode": "local_worker",
            "command": ["python", "-m", "mn_molecular_design_studio.engine.workflow_engine.core.local_worker", "--run-dir", str(job.run_dir)],
        },
    )


def enqueue_scannet(target_pdb: Path, chain_ids: list[str], mode: Literal["interface", "epitope", "idp"] = "interface", use_msa: bool = False, gpu_device: object = "0") -> Path:
    validate_surface_target(target_pdb, chain_ids or None)
    manifest = load_manifest("scannet")
    normalized_gpu = normalize_gpu_device(gpu_device)
    job = create_job(
        DETECTION_GROUP,
        job_type="ppi_detection",
        tool="scannet",
        inputs={"target_pdb": str(target_pdb), "chains": chain_ids},
        params={"mode": mode, "use_msa": use_msa, "gpu_device": normalized_gpu, "queue_resource": gpu_queue_resource(normalized_gpu)},
    )
    _write_detection_worker_request(job)
    return job.run_dir


def run_scannet(
    target_pdb: Path,
    chain_ids: list[str],
    mode: Literal["interface", "epitope", "idp"] = "interface",
    use_msa: bool = False,
    gpu_device: object = "0",
    existing_job: JobPaths | None = None,
    _defer_finish: bool = False,
) -> Path:
    chain_ids = validate_surface_target(target_pdb, chain_ids or None)
    if mode not in {"interface", "epitope", "idp"}:
        raise ValueError("ScanNet mode must be interface, epitope, or idp.")
    manifest = load_manifest("scannet")
    normalized_gpu = normalize_gpu_device(gpu_device)
    job = existing_job or create_job(
        DETECTION_GROUP,
        job_type="ppi_detection",
        tool="scannet",
        inputs={"target_pdb": str(target_pdb), "chains": chain_ids},
        params={"mode": mode, "use_msa": use_msa, "gpu_device": normalized_gpu, "queue_resource": gpu_queue_resource(normalized_gpu)},
    )
    input_pdb = artifact_path(job.run_dir, "input", "target.pdb")
    input_pdb.write_text(filter_pdb_text(target_pdb.read_text(errors="ignore"), keep_chains=set(chain_ids) or None))
    out_dir = artifact_path(job.run_dir, "scannet")
    out_dir.mkdir(parents=True, exist_ok=True)
    query = "/work/artifacts/input/target.pdb"
    args = ["python", "predict_bindingsites.py", query, "--name", "target", "--predictions_folder", "/work/artifacts/scannet", "--mode", mode]
    if not use_msa:
        args.append("--noMSA")
    steps = [{"name": "scannet", "command": _docker_base(manifest["image"], job.run_dir, gpu=True, gpu_device=normalized_gpu) + args}]
    rc = _run_shell_steps(job.run_dir, steps, append=_defer_finish)
    artifacts = _collect_artifacts(job.run_dir, out_dir, [("**/*.csv", "residue_score_table"), ("**/*.pdb", "pdb"), ("**/*.cxc", "chimera_script"), ("**/*.py", "script")])
    _finish_detection_tool(
        job,
        "scannet",
        rc == 0 and any(a["type"] == "residue_score_table" for a in artifacts),
        {
            "outputs": {"artifacts": artifacts},
            "metrics": {"return_code": rc, "artifact_count": len(artifacts)},
            "downstream_artifacts": {"residue_scores": [a["path"] for a in artifacts if a["type"] == "residue_score_table"]},
        },
        defer_finish=_defer_finish,
    )
    return job.run_dir


def enqueue_pesto(
    target_pdb: Path,
    chain_ids: list[str],
    gpu_device: object = "0",
    interface: Literal["protein", "dna_rna", "ion", "ligand", "lipid"] = "protein",
) -> Path:
    validate_surface_target(target_pdb, chain_ids or None)
    if interface not in {"protein", "dna_rna", "ion", "ligand", "lipid"}:
        raise ValueError("PeSTo interface must be protein, dna_rna, ion, ligand, or lipid.")
    load_manifest("pesto")
    normalized_gpu = normalize_gpu_device(gpu_device)
    job = create_job(
        DETECTION_GROUP,
        job_type="ppi_detection",
        tool="pesto",
        inputs={"target_pdb": str(target_pdb), "chains": chain_ids},
        params={
            "mode": interface,
            "interface": interface,
            "gpu_device": normalized_gpu,
            "queue_resource": gpu_queue_resource(normalized_gpu),
        },
    )
    _write_detection_worker_request(job)
    return job.run_dir


def run_pesto(
    target_pdb: Path,
    chain_ids: list[str],
    gpu_device: object = "0",
    existing_job: JobPaths | None = None,
    interface: Literal["protein", "dna_rna", "ion", "ligand", "lipid"] = "protein",
    _defer_finish: bool = False,
) -> Path:
    if interface not in {"protein", "dna_rna", "ion", "ligand", "lipid"}:
        raise ValueError("PeSTo interface must be protein, dna_rna, ion, ligand, or lipid.")
    chain_ids = validate_surface_target(target_pdb, chain_ids or None)
    manifest = load_manifest("pesto")
    normalized_gpu = normalize_gpu_device(gpu_device)
    job = existing_job or create_job(
        DETECTION_GROUP,
        job_type="ppi_detection",
        tool="pesto",
        inputs={"target_pdb": str(target_pdb), "chains": chain_ids},
        params={
            "mode": interface,
            "interface": interface,
            "gpu_device": normalized_gpu,
            "queue_resource": gpu_queue_resource(normalized_gpu),
        },
    )
    root = job.run_dir / "artifacts" / "pesto"
    input_pdb = artifact_path(job.run_dir, "pesto", "input", "target.pdb")
    source_text = target_pdb.read_text(errors="ignore")
    input_pdb.write_text(
        filter_pdb_text(
            source_text,
            keep_chains=set(chain_ids),
            remove_waters=True,
            remove_hetero=True,
        )
    )
    output_dir = artifact_path(job.run_dir, "pesto", "output")
    output_dir.mkdir(parents=True, exist_ok=True)
    steps = [
        {
            "name": "pesto",
            "command": [
                "docker",
                "run",
                "--rm",
                *docker_gpu_args(normalized_gpu),
                "-v",
                f"{job.run_dir}:/work",
                "-v",
                "/mnt/db/reference_files/pesto:/models:ro",
                manifest["image"],
            ]
            + [
                "--input",
                "/work/artifacts/pesto/input/target.pdb",
                "--output-dir",
                "/work/artifacts/pesto/output",
                "--interface",
                interface,
                "--device",
                "cuda",
                "--checkpoint",
                "/models/i_v4_1/model_ckpt.pt",
            ],
        }
    ]
    rc = _run_shell_steps(job.run_dir, steps, append=_defer_finish)
    artifacts = _collect_artifacts(
        job.run_dir,
        root,
        [("output/*.csv", "residue_score_table"), ("output/*.pdb", "pdb")],
    )
    residue_score_artifacts = [a["path"] for a in artifacts if a["type"] == "residue_score_table"]
    scored_pdb_artifacts = [a["path"] for a in artifacts if a["type"] == "pdb"]
    _finish_detection_tool(
        job,
        "pesto",
        rc == 0 and bool(residue_score_artifacts) and bool(scored_pdb_artifacts),
        {
            "outputs": {"artifacts": artifacts},
            "metrics": {"return_code": rc, "artifact_count": len(artifacts), "interface": interface},
            "downstream_artifacts": {
                "residue_scores": residue_score_artifacts,
                "scored_pdbs": scored_pdb_artifacts,
            },
        },
        defer_finish=_defer_finish,
    )
    return job.run_dir


def enqueue_surf2spot(target_pdb: Path, gpu_device: object = "0") -> Path:
    validate_surface_target(target_pdb)
    manifest = load_manifest("surf2spot")
    normalized_gpu = normalize_gpu_device(gpu_device)
    job = create_job(
        HOTSPOT_GROUP,
        job_type="hotspot_detection",
        tool="surf2spot",
        inputs={"target_pdb": str(target_pdb)},
        params={"mode": "HS", "gpu_device": normalized_gpu, "queue_resource": gpu_queue_resource(normalized_gpu)},
    )
    _write_detection_worker_request(job)
    return job.run_dir


def run_surf2spot(target_pdb: Path, gpu_device: object = "0", existing_job: JobPaths | None = None, _defer_finish: bool = False) -> Path:
    validate_surface_target(target_pdb)
    manifest = load_manifest("surf2spot")
    normalized_gpu = normalize_gpu_device(gpu_device)
    job = existing_job or create_job(
        HOTSPOT_GROUP,
        job_type="hotspot_detection",
        tool="surf2spot",
        inputs={"target_pdb": str(target_pdb)},
        params={"mode": "HS", "gpu_device": normalized_gpu, "queue_resource": gpu_queue_resource(normalized_gpu)},
    )
    input_pdb = artifact_path(job.run_dir, "surf2spot", "input", "target.pdb")
    sanitized_text, sanitize_metrics = sanitize_pdb_for_surface_tools(target_pdb.read_text(errors="ignore"))
    artifact_path(job.run_dir, "surf2spot", "source", "target.pdb").write_text(sanitized_text)
    input_pdb.write_text(_surf2spot_numbered_input(sanitized_text))
    base = ["docker", "run", "--rm", *docker_gpu_args(normalized_gpu), "-v", f"{job.run_dir}/artifacts/surf2spot:/work"]
    chain_ref = "/mnt/db/reference_files/surf2spot/chainsaw/saved_models"
    emb_ref = "/mnt/db/reference_files/surf2spot/model_emb/prot_t5_xl_half_uniref50-enc"
    common_mounts = ["-v", f"{chain_ref}:/opt/Surf2Spot/Surf2Spot/data/chainsaw/saved_models:ro"]
    emb_mount = ["-v", f"{emb_ref}:/opt/Surf2Spot/Surf2Spot/data/model_emb/prot_t5_xl_half_uniref50-enc:ro"]
    image = manifest["image"]
    steps = [
        {"name": "HS-preprocess", "command": base + common_mounts + [image, "Surf2Spot", "HS-preprocess", "-i", "/work/input", "-o", "/work/preprocess"]},
        {"name": "HS-craft", "command": base + common_mounts + emb_mount + [image, "Surf2Spot", "HS-craft", "-i", "/work/preprocess"]},
        {"name": "HS-predict", "command": base + common_mounts + emb_mount + [image, "Surf2Spot", "HS-predict", "-i", "/work/preprocess", "-o", "/work/predict", "--model", "/opt/Surf2Spot/model/HS/model.pt"]},
        {"name": "HS-draw", "command": base + common_mounts + emb_mount + [image, "Surf2Spot", "HS-draw", "-i", "/work/preprocess", "-o", "/work/predict"]},
    ]
    rc = _run_shell_steps(job.run_dir, steps, append=_defer_finish)
    root = job.run_dir / "artifacts" / "surf2spot"
    artifacts = _collect_artifacts(job.run_dir, root, [("predict/*.csv", "hotspot_table"), ("predict/*_pred.ply", "surface_mesh"), ("predict/*.pse", "pymol_session")])
    _finish_detection_tool(
        job,
        "surf2spot",
        rc == 0 and any(a["type"] == "hotspot_table" for a in artifacts) and any(a["type"] == "surface_mesh" for a in artifacts),
        {
            "outputs": {"artifacts": artifacts},
            "metrics": {"return_code": rc, "artifact_count": len(artifacts), **sanitize_metrics},
            "downstream_artifacts": {
                "hotspot_tables": [a["path"] for a in artifacts if a["type"] == "hotspot_table"],
                "hotspot_surfaces": [a["path"] for a in artifacts if a["type"] == "surface_mesh"],
            },
        },
        defer_finish=_defer_finish,
    )
    return job.run_dir


def enqueue_masif_seed(target_pdb: Path, chain_id: str, gpu_device: object = "0") -> Path:
    validate_surface_target(target_pdb, [chain_id])
    manifest = load_manifest("masif_seed")
    normalized_gpu = normalize_gpu_device(gpu_device)
    job = create_job(
        DETECTION_GROUP,
        job_type="surface_ppi_detection",
        tool="masif_seed",
        inputs={"target_pdb": str(target_pdb), "chain_id": chain_id},
        params={"mode": "target_site_surface", "gpu_device": normalized_gpu, "queue_resource": gpu_queue_resource(normalized_gpu)},
    )
    _write_detection_worker_request(job)
    return job.run_dir


def run_masif_seed(target_pdb: Path, chain_id: str, gpu_device: object = "0", existing_job: JobPaths | None = None, _defer_finish: bool = False) -> Path:
    chain_id = validate_surface_target(target_pdb, [chain_id])[0]
    manifest = load_manifest("masif_seed")
    normalized_gpu = normalize_gpu_device(gpu_device)
    job = existing_job or create_job(
        DETECTION_GROUP,
        job_type="surface_ppi_detection",
        tool="masif_seed",
        inputs={"target_pdb": str(target_pdb), "chain_id": chain_id},
        params={"mode": "target_site_surface", "gpu_device": normalized_gpu, "queue_resource": gpu_queue_resource(normalized_gpu)},
    )
    input_pdb = artifact_path(job.run_dir, "masif_seed", "input", "target.pdb")
    input_pdb.write_text(filter_pdb_text(target_pdb.read_text(errors="ignore"), keep_chains={chain_id}, remove_waters=True, remove_hetero=True))
    target_id = f"TGT1_{chain_id}"
    shell_target_id = shlex.quote(target_id)
    script = (
        "set -euxo pipefail; "
        "cd /opt/masif_seed/masif_seed_search/data/masif_targets; "
        f"./data_prepare_one.sh --file /job/artifacts/masif_seed/input/target.pdb {shell_target_id}; "
        f"./predict_site.sh {shell_target_id}; "
        f"./color_site.sh {shell_target_id}; "
        f"./compute_descriptors.sh {shell_target_id}; "
        "mkdir -p /job/artifacts/masif_seed/pred_surfaces /job/artifacts/masif_seed/pred_data /job/artifacts/masif_seed/target_run; "
        f"cp -av /opt/masif_seed/masif_seed_search/data/masif_targets/output/all_feat_3l/pred_surfaces/{target_id}.ply /job/artifacts/masif_seed/pred_surfaces/ 2>/dev/null || true; "
        f"cp -av /opt/masif_seed/masif_seed_search/data/masif_targets/output/all_feat_3l/pred_data/pred_{target_id}.npy /job/artifacts/masif_seed/pred_data/ 2>/dev/null || true; "
        f"cp -a /opt/masif_seed/masif_seed_search/data/masif_targets/targets/{target_id} /job/artifacts/masif_seed/target_run/ 2>/dev/null || true"
    )
    steps = [
        {
            "name": "masif-target-site",
            "command": [
                "docker",
                "run",
                "--rm",
                *docker_gpu_args(normalized_gpu),
                "-v",
                f"{job.run_dir}:/job",
                "-w",
                "/opt/masif_seed",
                manifest["image"],
                "bash",
                "-lc",
                script,
            ],
        }
    ]
    rc = _run_shell_steps(job.run_dir, steps, append=_defer_finish)
    root = job.run_dir / "artifacts" / "masif_seed"
    artifacts = _collect_artifacts(job.run_dir, root, [("pred_surfaces/*.ply", "surface_mesh"), ("pred_data/*.npy", "masif_site_prediction"), ("target_run/**/*", "masif_output")])
    _finish_detection_tool(
        job,
        "masif_seed",
        rc == 0 and any(a["type"] == "surface_mesh" for a in artifacts) and any(a["type"] == "masif_site_prediction" for a in artifacts),
        {
            "outputs": {"artifacts": artifacts, "masif_target_id": target_id},
            "metrics": {"return_code": rc, "artifact_count": len(artifacts)},
            "downstream_artifacts": {
                "predicted_surfaces": [a["path"] for a in artifacts if a["type"] == "surface_mesh"],
                "site_predictions": [a["path"] for a in artifacts if a["type"] == "masif_site_prediction"],
            },
        },
        defer_finish=_defer_finish,
    )
    return job.run_dir


def run_surface_analysis(
    target_pdb: Path,
    chain_ids: list[str],
    pesto_interface: Literal["protein", "dna_rna", "ion", "ligand", "lipid"] = "protein",
    scannet_mode: Literal["interface", "epitope", "idp"] = "interface",
    scannet_use_msa: bool = False,
    masif_chain_id: str | None = None,
    gpu_device: object = "0",
    existing_job: JobPaths | None = None,
) -> Path:
    """Run the four protein-surface predictors as one grouped target analysis."""
    chains = list(dict.fromkeys(str(chain).strip() for chain in chain_ids if str(chain).strip()))
    if not chains:
        raise ValueError("Select at least one protein chain for surface analysis.")
    masif_chain = str(masif_chain_id or chains[0]).strip()
    if masif_chain not in chains:
        raise ValueError("The MaSIF-seed chain must be one of the selected target chains.")
    if pesto_interface not in {"protein", "dna_rna", "ion", "ligand", "lipid"}:
        raise ValueError("PeSTo interface must be protein, dna_rna, ion, ligand, or lipid.")
    if scannet_mode not in {"interface", "epitope", "idp"}:
        raise ValueError("ScanNet mode must be interface, epitope, or idp.")
    validate_surface_target(target_pdb, chains)

    normalized_gpu = normalize_gpu_device(gpu_device)
    job = existing_job or create_job(
        DETECTION_GROUP,
        job_type="surface_analysis",
        tool="surface_analysis",
        inputs={"target_pdb": str(target_pdb), "chains": chains, "masif_chain_id": masif_chain},
        params={
            "pesto_interface": pesto_interface,
            "scannet_mode": scannet_mode,
            "scannet_use_msa": bool(scannet_use_msa),
            "gpu_device": normalized_gpu,
            "queue_resource": gpu_queue_resource(normalized_gpu),
        },
    )
    selected_target_pdb = artifact_path(job.run_dir, "surface_analysis", "input", "target.pdb")
    selected_target_pdb.write_text(
        filter_pdb_text(
            target_pdb.read_text(errors="ignore"),
            keep_chains=set(chains),
            remove_waters=True,
            remove_hetero=True,
        )
    )
    status_path = job.run_dir / "artifacts" / "surface_analysis" / "tool_status.json"
    write_json(status_path, {})
    methods = (
        ("pesto", lambda: run_pesto(selected_target_pdb, chains, normalized_gpu, existing_job=job, interface=pesto_interface, _defer_finish=True)),
        ("scannet", lambda: run_scannet(selected_target_pdb, chains, scannet_mode, scannet_use_msa, normalized_gpu, existing_job=job, _defer_finish=True)),
        ("surf2spot", lambda: run_surf2spot(selected_target_pdb, normalized_gpu, existing_job=job, _defer_finish=True)),
        ("masif_seed", lambda: run_masif_seed(selected_target_pdb, masif_chain, normalized_gpu, existing_job=job, _defer_finish=True)),
    )
    for tool, invoke in methods:
        try:
            invoke()
        except Exception as exc:
            statuses = read_json(status_path) if status_path.is_file() else {}
            statuses[tool] = {"success": False, "error": str(exc)}
            write_json(status_path, statuses)

    statuses = read_json(status_path) if status_path.is_file() else {}
    for tool, _invoke in methods:
        statuses.setdefault(tool, {"success": False, "error": "The tool did not report a result."})
    success = all(bool(statuses[tool].get("success")) for tool, _invoke in methods)

    artifact_specs = {
        "pesto": [("output/*.csv", "residue_score_table"), ("output/*.pdb", "pdb")],
        "scannet": [("**/*.csv", "residue_score_table"), ("**/*.pdb", "pdb"), ("**/*.cxc", "chimera_script"), ("**/*.py", "script")],
        "surf2spot": [("predict/*.csv", "hotspot_table"), ("predict/*_pred.ply", "surface_mesh"), ("predict/*.pse", "pymol_session")],
        "masif_seed": [("pred_surfaces/*.ply", "surface_mesh"), ("pred_data/*.npy", "masif_site_prediction"), ("target_run/**/*", "masif_output")],
    }
    artifacts: list[dict] = []
    for tool, patterns in artifact_specs.items():
        artifacts.extend(_collect_artifacts(job.run_dir, job.run_dir / "artifacts" / tool, patterns))
    summary_path = artifact_path(job.run_dir, "surface_analysis", "summary.json")
    write_json(
        summary_path,
        {
            "target": {"pdb": str(target_pdb), "chains": chains, "masif_chain_id": masif_chain},
            "settings": {"pesto_interface": pesto_interface, "scannet_mode": scannet_mode, "scannet_use_msa": bool(scannet_use_msa)},
            "tools": statuses,
        },
    )
    artifacts.append(Artifact("surface_analysis_summary", summary_path, "analysis_summary").to_json(job.run_dir))
    score_paths = [item["path"] for item in artifacts if item["type"] in {"residue_score_table", "hotspot_table"}]
    finish_job(
        job.run_dir,
        success,
        {
            "outputs": {"artifacts": artifacts, "tool_status": statuses},
            "metrics": {
                "tool_count": len(statuses),
                "successful_tools": sum(bool(statuses[tool].get("success")) for tool, _invoke in methods),
                "artifact_count": len(artifacts),
                "score_tables": len(score_paths),
            },
            "downstream_artifacts": {"surface_score_tables": score_paths, "summary": "artifacts/surface_analysis/summary.json"},
        },
    )
    return job.run_dir


def chains_for_pdb(path: Path) -> list[str]:
    return [chain["chain_id"] for chain in pdb_summary(path.read_text(errors="ignore"))["chains"]]
