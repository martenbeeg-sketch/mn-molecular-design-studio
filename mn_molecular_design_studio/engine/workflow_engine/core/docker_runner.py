from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Mapping, Sequence

from mn_molecular_design_studio.engine.workflow_engine.core.gpu import docker_gpu_args
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import update_status, write_json
from mn_molecular_design_studio.engine.workflow_engine.core.scheduler import apply_docker_cpu_limit
from mn_molecular_design_studio.engine.workflow_engine.runtime import reference_root


def build_docker_command(
    manifest: dict,
    run_dir: Path,
    args: list[str],
    reference_root_path: Path | None = None,
    gpu_device: object = "all",
    extra_mounts: Sequence[tuple[Path, str, bool]] | None = None,
    environment: Mapping[str, str] | None = None,
) -> list[str]:
    command = ["docker", "run", "--rm"]
    if manifest.get("gpu"):
        command.extend(docker_gpu_args(gpu_device))
    command.extend(["-v", f"{run_dir}:/work"])
    reference_path = reference_root_path or reference_root()
    if reference_path.exists():
        command.extend(["-v", f"{reference_path}:/ref:ro"])
    for host_path, container_path, read_only in extra_mounts or ():
        source = Path(host_path).expanduser().resolve()
        if not source.exists():
            raise FileNotFoundError(f"Docker mount source does not exist: {source}")
        destination = str(container_path).strip()
        if not destination.startswith("/") or ":" in destination or "\n" in destination:
            raise ValueError(f"Docker mount destination must be an absolute container path: {destination!r}")
        mount = f"{source}:{destination}" + (":ro" if read_only else "")
        command.extend(["-v", mount])
    for key, value in (environment or {}).items():
        name = str(key)
        setting = str(value)
        if not name or not name.replace("_", "").isalnum() or not name[0].isalpha() or "\n" in setting or "\r" in setting:
            raise ValueError(f"Invalid Docker environment entry: {name!r}")
        command.extend(["-e", f"{name}={setting}"])
    command.append(str(manifest["image"]))
    command.extend(args)
    return apply_docker_cpu_limit(command, run_dir)


def run_docker_job(
    run_dir: Path,
    manifest: dict,
    args: list[str],
    gpu_device: object = "all",
    *,
    extra_mounts: Sequence[tuple[Path, str, bool]] | None = None,
    environment: Mapping[str, str] | None = None,
) -> int:
    command = build_docker_command(
        manifest,
        run_dir,
        args,
        gpu_device=gpu_device,
        extra_mounts=extra_mounts,
        environment=environment,
    )
    write_json(run_dir / "command.json", {"mode": "docker", "command": command, "image": manifest.get("image")})
    update_status(run_dir, "running")
    with (run_dir / "stdout.log").open("w") as stdout, (run_dir / "stderr.log").open("w") as stderr:
        proc = subprocess.run(command, stdout=stdout, stderr=stderr, check=False)
    return int(proc.returncode)
