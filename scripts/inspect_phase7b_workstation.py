"""Read-only Phase 7B inventory. Never emits workstation credentials or environment."""
from __future__ import annotations

import hashlib
import argparse
import importlib.metadata
import json
from pathlib import Path
import platform
import socket
import subprocess
import sys
from datetime import datetime, timezone

from mn_molecular_design_studio.cluster import capabilities, config
from mn_molecular_design_studio.cluster.store import Store
from mn_molecular_design_studio.runtime import reference_root, runs_root


def command(*args):
    p = subprocess.run(args, capture_output=True, text=True, timeout=20)
    return {"returncode": p.returncode, "stdout": p.stdout.strip(), "stderr": p.stderr.strip()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--hash-models', action='store_true')
    args = parser.parse_args()
    from mn_molecular_design_studio.workspace import ensure_studio_workspace
    ensure_studio_workspace()
    import mn_molecular_design_studio as studio
    root = Path(studio.__file__).resolve().parent.parent
    settings = config.load()
    tags = ["mn-bindcraft2:v1.0.4-cu13-pae-v1", "mn-rfd4-proteina:0.0.3-cu130",
            "mn-boltz2-opt:2.2.1-cu130-batch2", "mn-protenix-opt:2.0.0-cu128-batch2",
            "mn-ligandmpnn:latest"]
    images = {}
    for tag in tags:
        p = command("docker", "image", "inspect", tag)
        if p["returncode"]:
            images[tag] = {"available": False, "error": p["stderr"]}
        else:
            row = json.loads(p["stdout"])[0]
            images[tag] = {"available": True, "id": row["Id"],
                           "content_fingerprint": capabilities.image_content_fingerprint(row)}
    references = []
    model_paths = [
        'rfd4-proteina/huggingface/hub/models--nvidia--RFD4-Proteina/snapshots/a06e2509c7033b010fc43614e9fd840303588502/config.json',
        'rfd4-proteina/huggingface/hub/models--nvidia--RFD4-Proteina/snapshots/a06e2509c7033b010fc43614e9fd840303588502/model.safetensors',
        'boltz_models/boltz2_conf.ckpt', 'boltz_models/boltz2_aff.ckpt', 'boltz_models/ccd.pkl',
        'protenix/checkpoint/protenix-v2.pt', 'protenix/common/components.cif.rdkit_mol.pkl',
        'alphafold_models/params_model_1_multimer_v3.npz', 'alphafold_models/params_model_1_ptm.npz',
        'mpnn/v48_020_epoch300_hyper.pt',
    ]
    for directory in ("rfd4-proteina", "boltz_models", "protenix", "alphafold_models", "mpnn"):
        folder = reference_root() / directory
        files = [reference_root() / p for p in model_paths if p.startswith(directory + '/')]
        references.append({"directory": directory, "files": [
            {"path": str(p.relative_to(reference_root())), "size": p.stat().st_size if p.is_file() else None}
            for p in files]})
    from mn_molecular_design_studio.cluster.gpu_runtime import profiles, compatible
    gpu_profiles = profiles(tags[:4]) if args.hash_models else {}
    checkpoint_hashes = capabilities.reference_fingerprints(model_paths) if args.hash_models else {}
    online, resources = capabilities.resources()
    nodes = Store().nodes()
    # Keep the published compatibility evidence; omit unrelated image inventories.
    for node in nodes:
        snap = node.get("snapshot", {})
        snap["images"] = {tag: snap.get("images", {}).get(tag) for tag in tags}
    modified = {}
    for path in (root / "mn_molecular_design_studio/run_app.py", root / "mn_molecular_design_studio/tool_references.py"):
        if path.is_file():
            modified[str(path.relative_to(root))] = hashlib.sha256(path.read_bytes()).hexdigest()
    packages = {}
    for name in ("streamlit", "fastapi", "numpy", "gemmi", "mn-compute-scheduler"):
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = None
    result = dict(evidence_class="Source inspection", observed_at=datetime.now(timezone.utc).isoformat(),
                  hostname=socket.gethostname(), python=sys.executable, python_version=platform.python_version(),
                  package_root=str(root), packages=packages, source_fingerprint=capabilities.code_revision(),
                  git_head=command("git", "-C", str(root), "rev-parse", "HEAD"),
                  git_status=command("git", "-C", str(root), "status", "--short"), preserved_files=modified,
                  docker_version=command("docker", "version", "--format", "{{.Server.Version}}"),
                  containers=command("docker", "ps", "--format", "{{.ID}} {{.Image}} {{.Names}}"),
                  services=command("systemctl", "--user", "show", "studio-api.service", "studio-agent.service",
                                   "-p", "Id", "-p", "ActiveState", "-p", "ExecStart", "-p", "WorkingDirectory"),
                  images=images, reference_root=str(reference_root()), references=references,
                  checkpoint_hashes=checkpoint_hashes, gpu_runtime_profiles=gpu_profiles,
                  gpu_image_compatibility={gpu: {tag: compatible(info, {tag: profile}) for tag, profile in gpu_profiles.items()}
                                           for gpu, info in capabilities.gpu_details().items()},
                  runs_root=str(runs_root()), gpu_details=capabilities.gpu_details(),
                  scheduler_online=online, scheduler_resources=resources,
                  cluster={k: settings[k] for k in ("mode", "name", "accept_jobs", "default_placement")},
                  nodes=nodes)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
