"""Host-side adapter for Studio's bundled ligand-complex preparation runtime."""

from __future__ import annotations

import hashlib
import json
import math
import subprocess
import sys
import traceback
from pathlib import Path
from typing import Any


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _merge_repaired_protein_with_source_heterogens(
    repaired_protein: str,
    source_complex: str,
) -> str:
    protein_lines = [line for line in repaired_protein.splitlines() if line.startswith("ATOM  ")]
    if not protein_lines:
        raise ValueError("Complex repair returned no protein ATOM records")
    heterogen_lines = [
        line
        for line in source_complex.splitlines()
        if line.startswith("HETATM")
        and line[17:20].strip().upper() not in {"HOH", "WAT", "DOD", "H2O"}
    ]
    headers = [
        line
        for line in source_complex.splitlines()
        if line.startswith(("HEADER", "TITLE ", "COMPND", "SOURCE", "CRYST1", "REMARK"))
    ]
    output = [*headers, *protein_lines, "TER", *heterogen_lines, "END", ""]
    serial = 0
    for index, line in enumerate(output):
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        serial += 1
        if serial > 99999:
            raise ValueError("Prepared complex exceeds the PDB 99,999-atom serial limit")
        output[index] = f"{line[:6]}{serial:5d}{line[11:]}"
    return "\n".join(output)


def _protein_only(complex_data: str) -> str:
    lines = [line for line in complex_data.splitlines() if line.startswith("ATOM  ")]
    if not lines:
        raise ValueError("Prepared complex contains no protein ATOM records")
    return "\n".join((*lines, "TER", "END", ""))


def _parse_ranges(value: object) -> dict[str, tuple[int, int]]:
    if not value:
        return {}
    payload = json.loads(str(value))
    if not isinstance(payload, dict):
        raise ValueError("Trim ranges must be a JSON object keyed by protein chain")
    result: dict[str, tuple[int, int]] = {}
    for chain, bounds in payload.items():
        if not isinstance(bounds, (list, tuple)) or len(bounds) != 2:
            raise ValueError(f"Trim range for chain {chain} must contain a start and end residue")
        start, end = int(bounds[0]), int(bounds[1])
        if start > end:
            raise ValueError(f"Trim range for chain {chain} ends before it starts")
        result[str(chain)] = (start, end)
    return result


def describe_gaps(source: str) -> list[dict]:
    from mn_ligand.workflows.protein_preparation import coordinate_gap_candidates, detect_modeller_internal_gaps

    candidates = coordinate_gap_candidates(source)
    detected = detect_modeller_internal_gaps(source, max_internal_gap=15)
    for candidate in candidates:
        known = next((
            gap for gap in detected.get(candidate["chain"], {}).get("gaps", [])
            if not gap.get("terminal") and (
                (gap.get("author_start"), gap.get("author_end")) == (candidate["author_start"], candidate["author_end"])
                or (gap.get("author_start") is None and gap.get("insertion_index") == candidate["insertion_index"]
                    and gap.get("length") == candidate["numbering_gap_length"])
            )
        ), {})
        candidate["suggestion"] = known.get("sequence", "")
    return candidates



def rebuild_confirmed_gaps(source: str, gaps: list[dict], work_dir: Path):
    from mn_ligand.workflows.protein_preparation import coordinate_gap_candidates, model_internal_gaps_with_modeller
    from mn_molecular_design_studio.protein_modelling import MAX_MODELED_GAP, observed_chains

    observed_chains(source)
    fields = ("chain", "left_resseq", "right_resseq")
    candidates = {tuple(row[key] for key in fields): row for row in coordinate_gap_candidates(source)}
    definitions = []
    for gap in gaps:
        key = tuple(gap[field] for field in fields)
        candidate = candidates.pop(key, None)
        if candidate is None:
            raise ValueError("A confirmed gap must match distinct observed flanks in the source PDB.")
        definitions.append({**candidate, "sequence": gap["sequence"],
                            "confirmed": True, "evidence": "User-confirmed missing sequence and coordinate flanks"})
    modeled, report = model_internal_gaps_with_modeller(
        source, work_dir=work_dir, max_internal_gap=MAX_MODELED_GAP,
        user_gap_definitions=definitions, confirmed_only=True,
    )
    if len(report["modeled_gaps"]) != len(definitions):
        raise ValueError("MODELLER did not rebuild every confirmed gap.")
    return modeled, report


def remodel_observed_loop(source: str, options: dict, work_dir: Path):
    from mn_ligand.modeller_runtime import modeller_python
    from mn_ligand.ligandx.lib.chemistry.preparation.target_validation import audit_target_geometry
    from mn_molecular_design_studio.protein_modelling import loop_residues, observed_chains

    labels = set(loop_residues(observed_chains(source), options))
    chain = options["chain"]

    def selected(line):
        return (line.startswith("ATOM  ") and len(line) >= 54 and (line[21].strip() or "_") == chain
                and line[22:26].strip() + line[26].strip() in labels)

    lines = source.splitlines()
    if any(line[16].strip() for line in lines if selected(line)):
        raise ValueError("Choose a loop with one atom conformation per residue.")
    work_dir.mkdir(parents=True, exist_ok=True)
    template = work_dir / "source.pdb"
    template.write_text(_protein_only(source), encoding="utf-8")
    native = Path(__file__).resolve().parents[1] / "vendor/ligand_runtime/mn_ligand/workflows/modeller_existing_loop_runner.py"
    process = subprocess.run(
        [str(modeller_python()), str(native), template.name, json.dumps(options)],
        cwd=work_dir, capture_output=True, text=True, check=False,
    )
    if process.returncode:
        raise RuntimeError(f"MODELLER loop refinement failed: {(process.stderr or process.stdout)[-2500:]}")
    results = json.loads((work_dir / "models.json").read_text(encoding="utf-8"))

    def score(result):
        value = result.get("dope")
        return value if isinstance(value, (int, float)) and math.isfinite(value) else float("inf")

    for model in sorted(results, key=score):
        candidate = (work_dir / Path(model["name"]).name).read_text(encoding="utf-8")
        coordinates = {line[12:27]: line[30:54] for line in candidate.splitlines() if selected(line)}
        if any(line[12:27] not in coordinates for line in lines if selected(line)):
            continue
        merged = "\n".join(line[:30] + coordinates[line[12:27]] + line[54:] if selected(line) else line
                           for line in lines) + "\n"
        try:
            loop_residues(observed_chains(merged), options)
        except ValueError:
            continue
        if audit_target_geometry(merged)["valid"]:
            return merged, {"enabled": True, "engine": "MODELLER LoopModel", **options,
                            "selected_model": model, "successful_models": results}
    raise ValueError("MODELLER did not produce a loop with matching atom identifiers and valid protein geometry.")


def run(input_path: Path, output_dir: Path, params_path: Path) -> None:
    from mn_ligand.core.artifacts import ArtifactRef, write_artifact_manifest
    from mn_ligand.core.jobs import JobRecord
    from mn_ligand.ligandx.lib.chemistry.preparation.target_validation import (
        prepare_target_for_publication,
    )
    from mn_ligand.workflows.protein_preparation import (
        create_protein_import_job,
        imported_target,
        prepared_target,
        run_protein_cleaning_job,
    )
    from mn_ligand.workflows.target_trimming import trim_pdb_data
    from mn_ligand.workflows.terminal_repair import create_terminal_repair_job

    from mn_molecular_design_studio.protein_modelling import parse_options

    params = json.loads(params_path.read_text(encoding="utf-8"))
    modelling = parse_options(params.get("modelling_options"))
    source_data = input_path.read_text(encoding="utf-8", errors="replace")
    import_job = create_protein_import_job(
        source_data,
        filename=input_path.name,
        source="studio_complex_upload",
    )
    source_ref = imported_target(import_job)
    repair_job = None
    repair_report: dict[str, Any] = {"enabled": False, "reason": "Repair was disabled"}
    if bool(params.get("repair", True)):
        repair_job, native_payload = run_protein_cleaning_job(
            import_job.run_id,
            image="mn-md:cu128",
            use_gpu=False,
            map_modified_residues=bool(params.get("map_modified_residues", True)),
            skip_terminal_missing_residues=bool(params.get("skip_terminal_missing_residues", True)),
            max_internal_gap=int(params.get("max_internal_gap", 15)),
            refine_rebuilt_positions=True,
            preserve_nonwater_heterogens=False,
        )
        if repair_job.status != "completed" or not native_payload.get("success"):
            raise RuntimeError(
                str(native_payload.get("error") or repair_job.result.get("error") or "Protein repair failed")
            )
        repaired_ref = prepared_target(repair_job)
        repaired_path = repaired_ref.resolve(repair_job.run_dir, must_exist=True)
        if repaired_path is None:
            raise FileNotFoundError(repaired_ref.path)
        complex_path = repair_job.run_dir / "artifacts" / "prepared" / "studio_prepared_complex.pdb"
        complex_path.parent.mkdir(parents=True, exist_ok=True)
        complex_path.write_text(
            _merge_repaired_protein_with_source_heterogens(
                repaired_path.read_text(encoding="utf-8", errors="replace"),
                source_data,
            ),
            encoding="utf-8",
        )
        complex_ref = ArtifactRef.from_path(
            repair_job.run_dir,
            complex_path,
            "prepared_complex",
            role="repaired_complex",
            metadata={"source_run_id": import_job.run_id, "ligand_retained": True},
        )
        write_artifact_manifest(
            repair_job.run_dir,
            [*(repair_job.artifact_manifest.artifacts if repair_job.artifact_manifest else ()), complex_ref],
        )
        repair_job = JobRecord.load(repair_job.run_dir, task_group="protein-cleaning")
        repair_report_path = next(
            (
                artifact.resolve(repair_job.run_dir, must_exist=True)
                for artifact in (repair_job.artifact_manifest.artifacts if repair_job.artifact_manifest else ())
                if artifact.artifact_type == "repair_report"
            ),
            None,
        )
        repair_report = (
            json.loads(repair_report_path.read_text(encoding="utf-8"))
            if repair_report_path is not None
            else {}
        )
        complex_data = complex_path.read_text(encoding="utf-8", errors="replace")
    else:
        complex_data = source_data if modelling else _merge_repaired_protein_with_source_heterogens(source_data, source_data)

    extension_report: dict[str, Any] = {"enabled": False}
    extension_sequence = str(params.get("extension_sequence") or "").strip()
    if extension_sequence:
        from mn_ligand.workflows.terminal_repair import normalize_extension_sequence

        extension_sequence = normalize_extension_sequence(extension_sequence)
        extension_source_job = repair_job or import_job
        extension_source_artifact = (
            next(
                artifact
                for artifact in (repair_job.artifact_manifest.artifacts if repair_job and repair_job.artifact_manifest else ())
                if artifact.artifact_type == "prepared_complex"
            )
            if repair_job
            else source_ref
        )
        extension_source_path = (
            extension_source_artifact.resolve(extension_source_job.run_dir, must_exist=True)
            if repair_job
            else source_ref.resolve(import_job.run_dir, must_exist=True)
        )
        if extension_source_path is None:
            raise FileNotFoundError("The complex selected for terminal extension is unavailable")
        extension_job = create_terminal_repair_job(
            source_job=extension_source_job,
            source_artifact=extension_source_artifact,
            source_path=extension_source_path,
            chain=str(params.get("extension_chain") or ""),
            extension_sequence=extension_sequence,
            model_count=int(params.get("extension_models", 10)),
            sequence_origin="user_entered",
        )
        if extension_job.status != "completed":
            raise RuntimeError(str(extension_job.result.get("error") or "C-terminal MODELLER extension failed"))
        extension_complex = next(
            (
                artifact.resolve(extension_job.run_dir, must_exist=True)
                for artifact in (extension_job.artifact_manifest.artifacts if extension_job.artifact_manifest else ())
                if artifact.artifact_type == "prepared_complex"
            ),
            None,
        )
        if extension_complex is None:
            raise FileNotFoundError("MODELLER extension did not publish a prepared complex")
        complex_data = extension_complex.read_text(encoding="utf-8", errors="replace")
        extension_report = dict(extension_job.result)

    gap_report = {"enabled": False}
    if modelling.get("gaps"):
        complex_data, gap_report = rebuild_confirmed_gaps(
            complex_data, modelling["gaps"], output_dir / "confirmed_gaps",
        )
    loop_report = {"enabled": False}
    if "loop" in modelling:
        complex_data, loop_report = remodel_observed_loop(
            complex_data, modelling["loop"], output_dir / "loop_models",
        )

    trim_ranges = _parse_ranges(params.get("trim_ranges_json"))
    trim_report: dict[str, Any] = {"enabled": False}
    if trim_ranges:
        complex_data, trim_report = trim_pdb_data(complex_data, trim_ranges)
        trim_report["enabled"] = True

    prepared_target_data, target_validation = prepare_target_for_publication(_protein_only(complex_data))
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "prepared_complex.pdb").write_text(complex_data, encoding="utf-8")
    (output_dir / "prepared_target.pdb").write_text(prepared_target_data, encoding="utf-8")
    report = {
        "source_sha256": hashlib.sha256(input_path.read_bytes()).hexdigest(),
        "output_sha256": {name: hashlib.sha256((output_dir / name).read_bytes()).hexdigest()
                          for name in ("prepared_complex.pdb", "prepared_target.pdb")},
        "source_import": {
            "run_id": import_job.run_id,
            "source_artifact": source_ref.to_dict(),
            "ligands": import_job.metadata.get("ligands") or [],
            "receptor": import_job.metadata.get("receptor") or {},
        },
        "protein_repair": repair_report,
        "terminal_extension": extension_report,
        "confirmed_gap_repair": gap_report,
        "loop_remodelling": loop_report,
        "trimming": trim_report,
        "target_validation": target_validation,
        "parameters": params,
    }
    _write_json(output_dir / "repair_report.json", report)


def main(argv: list[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if len(args) == 2 and args[0] == "--inspect":
        print(json.dumps(describe_gaps(Path(args[1]).read_text(encoding="utf-8"))))
        return 0
    if len(args) != 3:
        raise SystemExit("usage: ligand_complex_prepare_runner.py INPUT_PDB OUTPUT_DIR PARAMS_JSON")
    input_path, output_dir, params_path = (Path(item).expanduser().resolve() for item in args)
    try:
        run(input_path, output_dir, params_path)
    except Exception as exc:
        output_dir.mkdir(parents=True, exist_ok=True)
        _write_json(output_dir / "runner_error.json", {"error": str(exc), "traceback": traceback.format_exc()})
        print(traceback.format_exc(), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
