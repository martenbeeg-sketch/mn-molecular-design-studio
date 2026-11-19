"""Scheduled, offline BinderScope analysis of exact saved designs or native campaigns."""
from __future__ import annotations

import copy
import json
import shutil
import subprocess
from pathlib import Path
from typing import Literal

from mn_molecular_design_studio.binderscope import COLLECTION_SCHEMA, IMAGE, METHOD_VERSION, REVISION, protocol_revision, validate_assessment_settings
from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import create_job, finish_job, read_json, update_status, write_json
from mn_molecular_design_studio.engine.workflow_engine.core.docker_runner import build_docker_command
from mn_molecular_design_studio.engine.workflow_engine.workflows.analysis import _analysis_job, _finish_metric_analysis
from mn_molecular_design_studio.evaluations import file_digest


RUNTIME = Path(__file__).resolve().parents[1] / "tools" / "binderscope.py"
METRICS = {
    "dsasa": "buried_area_angstrom2", "n_hbonds": "hydrogen_bond_count",
    "n_salt": "salt_bridge_count", "n_unsat": "buried_unpaired_polar_atom_count",
    "unsat_density": "unpaired_polar_atoms_per_1000_angstrom2",
    "hyd_frac": "hydrophobic_burial_fraction", "hphob_patch": "largest_surface_patch_angstrom2",
    "helix_frac": "binder_helix_fraction", "strand_frac": "binder_strand_fraction",
}


def validate_request(parameters, *, require_files=True):
    native = parameters.get("native_run_dirs") or []
    saved = parameters.get("source_run_dir") or parameters.get("candidates_jsonl")
    if bool(native) == bool(saved):
        raise ValueError("Select either saved SQL designs or native BindCraft campaigns, not both.")
    if not 1 <= int(parameters.get("cpu_cores", 4)) <= 64:
        raise ValueError("BinderScope requires between 1 and 64 CPU cores.")
    if native:
        if not parameters.get("target_pdb"):
            raise ValueError("Native campaigns require their input target structure.")
    else:
        if not parameters.get("source_run_dir") or not parameters.get("candidates_jsonl"):
            raise ValueError("Saved BinderScope analysis requires an authorized structure collection.")
        if not parameters.get("use_recorded_settings") and parameters.get("bfactors_are_plddt") is not True:
            raise ValueError("BinderScope assumes coordinate B-factors encode pLDDT. Confirm this for the selected models; otherwise use the geometry-only analysis tools.")
        if not parameters.get("use_recorded_settings") and parameters.get("component_format") not in {"linear protein", "linear peptide", "cyclic peptide"}:
            raise ValueError("Declare the saved binder topology before requesting BinderScope terminal analysis.")
        if parameters.get("reference_cif") or parameters.get("target_pdb") or parameters.get("profile_path"):
            raise ValueError("Reference/scaffold/specificity profiles require native campaign mode; saved cohorts do not infer those associations.")
    if require_files:
        for path in [*native, *[parameters[key] for key in ("target_pdb", "reference_cif", "reference_data_dir", "profile_path", "source_run_dir", "candidates_jsonl") if parameters.get(key)]]:
            if not Path(path).exists():
                raise ValueError("A selected BinderScope input is unavailable on this workstation.")
        if saved and parameters.get("use_recorded_settings"):
            _saved_inputs(Path(parameters["source_run_dir"]).resolve(), Path(parameters["candidates_jsonl"]).resolve(),
                          use_recorded_settings=True)


def _saved_inputs(source, candidates_path, *, use_recorded_settings=False):
    manifest = read_json(source / "complex_analysis_collection.json")
    if manifest.get("schema") != COLLECTION_SCHEMA:
        raise ValueError("BinderScope saved inputs must come from its authorized collection adapter.")
    candidates = read_candidates(candidates_path)
    items = {item["candidate_id"]: item for item in manifest["items"]}
    if not candidates or len(items) != len(candidates):
        raise ValueError("The saved collection inventory is inconsistent.")
    for candidate in candidates:
        selection = candidate.get("raw_metadata", {}).get("complex_analysis_selection") or {}
        if selection != items.get(candidate["candidate_id"]):
            raise ValueError("The selected design does not match its authorized collection inventory.")
        if use_recorded_settings:
            validate_assessment_settings(selection.get("assessment"))
        path = Path(candidate["complex_pdb"]).resolve(strict=True)
        path.relative_to(source)
        if file_digest(path) != selection["input_structure_digest"]:
            raise ValueError("A selected structure changed after authorization.")
        reference = selection.get("target_reference")
        if reference:
            reference_path = Path(reference["path"]).resolve(strict=True)
            reference_path.relative_to(source)
            if file_digest(reference_path) != reference["input_structure_digest"]:
                raise ValueError("An original target reference changed after authorization.")
        if len(candidate.get("binder_chains") or []) != 1 or len(candidate.get("target_chains") or []) != 1:
            raise ValueError("Saved BinderScope designs require one binder and one target chain.")
    return candidates


def run_binderscope(
    source_run_dir: Path | None = None,
    candidates_jsonl: Path | None = None,
    native_run_dirs: list[Path] | None = None,
    target_pdb: Path | None = None,
    reference_cif: Path | None = None,
    reference_data_dir: Path | None = None,
    profile_path: Path | None = None,
    uniprot: str = "",
    reference_chain: str = "",
    partner_chains: str = "none",
    component_format: Literal["unspecified", "linear protein", "linear peptide", "cyclic peptide"] = "unspecified",
    bfactors_are_plddt: bool = False,
    use_recorded_settings: bool = False,
    cpu_cores: int = 4,
) -> Path:
    parameters = dict(locals())
    validate_request(parameters)
    saved = not native_run_dirs
    candidates = []
    if saved:
        source = Path(source_run_dir).resolve(strict=True)
        candidates_path = Path(candidates_jsonl).resolve(strict=True)
        candidates_path.relative_to(source)
        candidates = _saved_inputs(source, candidates_path, use_recorded_settings=use_recorded_settings)
        job = _analysis_job(source, candidates_path, tool="binderscope", params={
            "bfactors_are_plddt": bfactors_are_plddt, "component_format": component_format,
            "cpu_cores": cpu_cores, "method_version": METHOD_VERSION, "upstream_revision": REVISION,
            "use_recorded_settings": use_recorded_settings,
            "design_settings": {candidate["candidate_id"]: candidate["raw_metadata"]["complex_analysis_selection"].get("assessment")
                                for candidate in candidates} if use_recorded_settings else {},
            "protocol_revision": protocol_revision(),
        }, note="Integrated BinderScope analysis; original design admission is unchanged.", mode="docker")
    else:
        job = create_job("analysis", job_type="analysis", tool="binderscope", inputs={
            "native_run_dirs": [str(path) for path in native_run_dirs], "target_pdb": str(target_pdb),
            "reference_cif": str(reference_cif) if reference_cif else None,
            "profile_path": str(profile_path) if profile_path else None,
            "reference_data_dir": str(reference_data_dir) if reference_data_dir else None,
            "reference_chain": reference_chain, "partner_chains": partner_chains, "uniprot": uniprot,
        }, params={"cpu_cores": cpu_cores, "method_version": METHOD_VERSION, "upstream_revision": REVISION,
                   "protocol_revision": protocol_revision()})
    run_dir = job.run_dir
    inputs = run_dir / "input"
    inputs.mkdir(parents=True, exist_ok=True)
    request = {"mode": "saved" if saved else "native", "input_root": "/work/input",
               "output_root": "/work/artifacts/analysis/binderscope", "cpu_cores": cpu_cores,
               "bfactors_are_plddt": bfactors_are_plddt, "component_format": component_format,
               "use_recorded_settings": use_recorded_settings,
               "adapter_sha256": file_digest(RUNTIME)}
    request["context_sha256"] = file_digest(RUNTIME.with_name("binderscope_context.py"))
    mounts = []
    staged = []
    for candidate in candidates:
        candidate = copy.deepcopy(candidate)
        selection = candidate["raw_metadata"]["complex_analysis_selection"]
        original = Path(candidate["complex_pdb"])
        destination = inputs / (candidate["candidate_id"] + original.suffix)
        shutil.copyfile(original, destination)
        candidate["complex_pdb"] = str(destination)
        staged.append(candidate)
        request_item = {**selection, "path": "/work/input/" + destination.name,
            "assessment": selection["assessment"] if use_recorded_settings else {
                "component_format": component_format,
                "b_factor_interpretation": "plddt_0_100" if bfactors_are_plddt else "uninterpreted",
                "topology_basis": "Explicit legacy job setting.", "confidence_basis": "Explicit legacy job setting.",
            },
            "label": str(selection.get("source", {}).get("record_name") or candidate["candidate_id"])}
        if selection.get("target_reference"):
            reference = selection["target_reference"]
            reference_destination = inputs / (candidate["candidate_id"] + "-target" + Path(reference["path"]).suffix)
            shutil.copyfile(reference["path"], reference_destination)
            request_item["target_reference"] = {**reference, "path": "/work/input/" + reference_destination.name}
        request.setdefault("items", []).append(request_item)
    if not saved:
        request["native_run_dirs"] = []
        for index, path in enumerate(native_run_dirs):
            destination = f"/source/campaign{index + 1}"
            mounts.append((Path(path), destination, True))
            request["native_run_dirs"].append(destination)
        for field, path in (("target_pdb", target_pdb), ("reference_cif", reference_cif), ("profile", profile_path)):
            if path is not None:
                destination = inputs / (field + Path(path).suffix)
                shutil.copyfile(path, destination)
                request[field] = "/work/input/" + destination.name
                request.setdefault("input_sha256", {})[field] = file_digest(destination)
        request.update({"uniprot": uniprot, "reference_chain": reference_chain, "partner_chains": partner_chains})
    write_json(inputs / "request.json", request)
    environment = {"HOME": "/tmp", "OMP_NUM_THREADS": str(cpu_cores), "OPENBLAS_NUM_THREADS": "1"}
    if reference_data_dir is not None:
        mounts.append((reference_data_dir, "/reference-data", True))
        environment["BINDERSCOPE_REFDATA_DIR"] = "/reference-data"
    command = build_docker_command({"image": IMAGE, "gpu": False}, run_dir,
        ["python", "/opt/studio/binderscope_runtime.py", "/work/input/request.json"], extra_mounts=mounts,
        use_host_user=True, environment=environment, working_dir="/work")
    command[2:2] = ["--network", "none"]
    inspection = subprocess.run(["docker", "image", "inspect", IMAGE, "--format", "{{.Id}}"], capture_output=True, text=True, check=False)
    if inspection.returncode:
        finish_job(run_dir, False, {"error": "BinderScope runtime is not built. See containers/binderscope/Dockerfile."})
        return run_dir
    write_json(run_dir / "command.json", {"mode": "docker", "command": command, "image": IMAGE,
                                          "image_digest": inspection.stdout.strip(), "upstream_revision": REVISION})
    update_status(run_dir, "running")
    from mn_molecular_design_studio.native_execution import record_process
    with (run_dir / "stdout.log").open("w") as stdout, (run_dir / "stderr.log").open("w") as stderr:
        process = record_process(run_dir, command, native_invocations=None, runner=subprocess.run,
                                 stdout=stdout, stderr=stderr, check=False)
    report_path = run_dir / "artifacts/analysis/binderscope/binderscope_report.json"
    if process.returncode or not report_path.is_file():
        finish_job(run_dir, False, {"error": "BinderScope analysis failed; inspect the job log. No design admission was changed."})
        return run_dir
    report = json.loads(report_path.read_text())
    vendor_root = Path(__file__).resolve().parents[3] / "vendor/binderscope"
    shutil.copyfile(vendor_root / "LICENSE", report_path.parent / "BINDERSCOPE-LICENSE.txt")
    shutil.make_archive(str(run_dir / "artifacts/analysis/binderscope-results"), "zip", report_path.parent)
    if not saved:
        finish_job(run_dir, True, {"summary": {"designs": len(report["results"]), "upstream_revision": REVISION},
            "outputs": {"artifacts": [{"name": "BinderScope report", "path": "artifacts/analysis/binderscope-results.zip", "type": "analysis_summary"}]}})
        return run_dir
    measured = {item["candidate_id"]: item for item in report["results"]}
    if set(measured) != {candidate["candidate_id"] for candidate in staged}:
        finish_job(run_dir, False, {"error": "BinderScope did not return exactly the authorized structure inventory."})
        return run_dir
    for candidate in staged:
        result = measured[candidate["candidate_id"]]
        candidate["metrics"] = {"binderscope_" + label: result["analysis"].get(name) for name, label in METRICS.items()}
        candidate["raw_metadata"]["binderscope"] = result
        hotspot = result.get("hotspot_compliance") or {}
        for name in ("requested_count", "contacted_count", "missed_count", "unmapped_count", "coverage_fraction"):
            if name in hotspot:
                candidate["metrics"]["binderscope_hotspot_" + name] = hotspot[name]
        candidate["metrics"]["binderscope_hotspot_status"] = hotspot.get("status", "unavailable")
        comparison = result.get("target_reference_comparison") or {}
        candidate["metrics"]["binderscope_target_reference_ca_rmsd_angstrom"] = comparison.get("target_ca_rmsd_angstrom")
        candidate["raw_metadata"]["analysis_parameters"] = {
            **result.get("source", {}).get("assessment", {}),
            "use_recorded_settings": use_recorded_settings,
            "method_version": METHOD_VERSION, "upstream_revision": REVISION, "protocol_revision": protocol_revision(),
            "image_digest": inspection.stdout.strip(), "source_manifest_sha256": report["source_manifest_sha256"],
            "target_context_sha256": report["target_context_sha256"],
        }
    return _finish_metric_analysis(job, staged, tool="binderscope", table_name="binderscope_metrics.csv",
        metric_names=["binderscope_" + name for name in METRICS.values()],
        required_metrics=["binderscope_buried_area_angstrom2"], error_key="binderscope_error")
