"""Run RFD4's native input preparation without submitting an inference job."""
from __future__ import annotations

import json
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any

from mn_molecular_design_studio import rfd4_proteina_input as inputs
from mn_molecular_design_studio.engine.workflow_engine.core.docker_runner import build_docker_command
from mn_molecular_design_studio.engine.workflow_engine.core.manifests import load_manifest


def validate_inputs(parameters: dict[str, Any], *, timeout_seconds: int = 300) -> dict[str, Any]:
    """Validate and prepare the current form inputs inside the pinned image.

    This creates only temporary files. The Docker command uses the runtime's
    prepare-only flag and explicitly disables GPU access, so it cannot launch
    RFD4 inference or scoring and does not create a Studio job.
    """
    if type(timeout_seconds) is not int or not 1 <= timeout_seconds <= 1800:
        raise ValueError("Validation timeout must be between 1 and 1800 seconds.")

    payload = json.loads(json.dumps(parameters, default=str))
    with tempfile.TemporaryDirectory(prefix="rfd4-proteina-prepare-") as temporary:
        run_dir = Path(temporary)
        protein_sequences = str(payload.pop("protein_sequences", "") or "").strip()
        if protein_sequences:
            if payload.get("target_structure"):
                raise ValueError("Provide either protein sequences or a structure, not both.")
            from mn_molecular_design_studio.engine.workflow_engine.workflows.rfd4_proteina import sequence_structure

            payload["target_structure"] = str(sequence_structure(protein_sequences, run_dir / "sequence_input.cif"))

        prepared = inputs.prepare_request(payload, run_dir=run_dir)
        runtime_dir = prepared["root"] / "runtime"
        runtime_dir.mkdir(parents=True, exist_ok=True)
        runtime_source = Path(__file__).with_name("rfd4_runtime.py")
        runtime_path = runtime_dir / runtime_source.name
        shutil.copy2(runtime_source, runtime_path)

        manifest = {**load_manifest("rfd4_proteina"), "entrypoint": "python", "gpu": False}
        args = [
            "/work/" + runtime_path.relative_to(run_dir).as_posix(),
            "/work/" + prepared["request_path"].relative_to(run_dir).as_posix(),
            "--prepare-only",
        ]
        command = build_docker_command(
            manifest,
            run_dir,
            args,
            environment={"PYTHONUNBUFFERED": "1"},
        )
        try:
            result = subprocess.run(
                command,
                capture_output=True,
                text=True,
                check=False,
                timeout=timeout_seconds,
            )
        except subprocess.TimeoutExpired as exc:
            raise RuntimeError(f"RFD4 input preparation exceeded {timeout_seconds} seconds.") from exc

        if result.returncode:
            detail = (result.stderr or result.stdout or "Container exited without diagnostic output.").strip()
            raise RuntimeError("RFD4 input preparation failed:\n" + detail[-4000:])

        inventory_path = prepared["root"] / "prepared_inputs.json"
        if not inventory_path.is_file():
            raise RuntimeError("RFD4 container exited without writing its prepared-input report.")
        inventory = json.loads(inventory_path.read_text())
        if not isinstance(inventory, list) or not inventory:
            raise RuntimeError("RFD4 prepared-input report contains no validated inputs.")
        if any(item.get("validation") != "passed" for item in inventory):
            raise RuntimeError("RFD4 prepared-input report contains an input that did not pass validation.")

        return {
            "count": len(inventory),
            "items": [
                {
                    "name": item.get("name", "input"),
                    "atoms": item.get("atoms"),
                    "designed_chains": item.get("designed_chains", []),
                    "target_chains": item.get("target_chains", []),
                    "warnings": item.get("warnings", []),
                    "conditions": item.get("conditions", ""),
                }
                for item in inventory
            ],
        }
