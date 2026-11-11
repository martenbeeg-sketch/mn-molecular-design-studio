"""Run CPU release tests without exposing production storage or the Docker socket."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import subprocess


def build_command(*, source, runtime_root, python, image, output, arguments, user):
    source, runtime_root, python, output = map(Path, (source, runtime_root, python, output))
    if not re.fullmatch(r"sha256:[0-9a-f]{64}", image):
        raise ValueError("Use an inspected immutable local image ID, not a mutable tag.")
    if not python.is_relative_to(runtime_root) or runtime_root == Path("/"):
        raise ValueError("The Python executable must be inside a narrowly mounted runtime root.")
    if output == source or output.is_relative_to(source) or source.is_relative_to(output):
        raise ValueError("Test evidence must be separate from the read-only source.")
    environment = {
        "HOME": "/tmp/teaching-user", "XDG_CACHE_HOME": "/tmp/cache", "TMPDIR": "/tmp",
        "PYTHONDONTWRITEBYTECODE": "1", "PYTHONPATH": f"{source}:{source}/tests",
        "MN_MOLECULAR_DESIGN_STUDIO_HOME": "/tmp/studio-home",
        "MN_MOLECULAR_DESIGN_STUDIO_APP_HOME": "/tmp/studio-home",
        "MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR": "/tmp/studio-runs",
        "MN_MOLECULAR_DESIGN_STUDIO_REFERENCE_DIR": "/tmp/studio-references",
        "MDS_CLUSTER_HOME": "/tmp/studio-cluster",
        "MN_MOLECULAR_DESIGN_STUDIO_AUTH_MODE": "production",
        "MN_MOLECULAR_DESIGN_STUDIO_SCHEDULER_MODE": "local",
        "MDS_GPU_TESTS": "0", "MDS_QC_DOCKER_TESTS": "0", "MDS_RFD4_DOCKER_TESTS": "0",
        "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1",
    }
    command = ["docker", "run", "--rm", "--pull=never", "--runtime=runc", "--network=none",
               "--read-only", "--cap-drop=ALL", "--security-opt=no-new-privileges",
               "--pids-limit=512", "--memory=8g", "--cpus=4", "--user", user,
               "--tmpfs", "/tmp:rw,nosuid,nodev,size=4g", "--workdir", str(source),
               "--mount", f"type=bind,src={source},dst={source},readonly",
               "--mount", f"type=bind,src={runtime_root},dst={runtime_root},readonly",
               "--mount", f"type=bind,src={output},dst=/evidence"]
    for name, value in environment.items():
        command.extend(["--env", f"{name}={value}"])
    command.extend(["--entrypoint", str(python), image, "-m", "pytest",
                    "-o", "cache_dir=/tmp/pytest-cache", "--junitxml=/evidence/junit.xml", *arguments])
    return command


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--runtime-root", type=Path, required=True)
    parser.add_argument("--python", type=Path, required=True)
    parser.add_argument("--image", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("arguments", nargs=argparse.REMAINDER)
    options = parser.parse_args()
    source, runtime, python, output = [path.resolve() for path in
        (options.source, options.runtime_root, options.python, options.output)]
    if not (source / "mn_molecular_design_studio/api.py").is_file() or not python.is_file():
        parser.error("The source or Python runtime is absent.")
    arguments = options.arguments[1:] if options.arguments[:1] == ["--"] else options.arguments
    security = json.loads(subprocess.check_output(["docker", "info", "--format", "{{json .SecurityOptions}}"], text=True))
    rootless = "name=rootless" in security
    container_user = "0:0" if rootless else f"{os.getuid()}:{os.getgid()}"
    command = build_command(source=source, runtime_root=runtime, python=python, image=options.image,
                            output=output, arguments=arguments, user=container_user)
    inspected = subprocess.check_output(["docker", "image", "inspect", "--format", "{{.Id}}", options.image], text=True).strip()
    if inspected != options.image:
        parser.error("Local image identity differs from the requested immutable ID.")
    output.mkdir(parents=True, mode=0o700, exist_ok=False)
    started = datetime.now(timezone.utc).isoformat()
    with (output / "pytest.log").open("w") as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
    receipt = {"started_at": started, "finished_at": datetime.now(timezone.utc).isoformat(),
               "hostname": os.uname().nodename, "source": str(source), "python": str(python),
               "image_id": inspected, "command": command, "exit_code": result.returncode,
               "network": "none", "production_storage_mounted": False, "docker_socket_mounted": False,
               "scientific_inference": False, "rootless_daemon": rootless, "container_user": container_user}
    (output / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"hostname": receipt["hostname"], "output": str(output), "exit_code": result.returncode}))
    raise SystemExit(result.returncode)


if __name__ == "__main__":
    main()
