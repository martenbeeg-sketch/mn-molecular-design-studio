from __future__ import annotations

import hashlib
import inspect
import json
import subprocess
from pathlib import Path


def code_revision():
    root = Path(__file__).resolve().parents[1]
    digest = hashlib.sha256()
    for path in sorted(root.rglob("*.py")):
        digest.update(path.relative_to(root).as_posix().encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()


def images():
    try:
        result = subprocess.run(["docker", "image", "ls", "--format", "{{.Repository}}:{{.Tag}}"], capture_output=True, text=True, timeout=10, check=True)
        tags = sorted({line.strip() for line in result.stdout.splitlines() if line.strip() and "<none>" not in line})
        inventory = {}
        # Docker's image ID is a config digest with the classic image store and
        # a manifest digest with the containerd image store. Inspect the image
        # contents instead so equivalent images compare equally across hosts.
        for start in range(0, len(tags), 32):
            batch = tags[start:start + 32]
            inspected = subprocess.run(
                ["docker", "image", "inspect", "--format", "{{json .}}", *batch],
                capture_output=True, text=True, timeout=20, check=True,
            )
            rows = [json.loads(line) for line in inspected.stdout.splitlines() if line.strip()]
            for requested, image in zip(batch, rows):
                fingerprint = image_content_fingerprint(image)
                if not fingerprint:
                    continue
                for tag in image.get("RepoTags") or [requested]:
                    if tag and "<none>" not in tag:
                        inventory[tag] = fingerprint
        return inventory
    except (OSError, subprocess.SubprocessError, ValueError):
        return {}


def image_content_fingerprint(image):
    """Return a portable identity for an image's executable contents."""
    stable = {
        "os": image.get("Os"),
        "architecture": image.get("Architecture"),
        "variant": image.get("Variant"),
        "os_version": image.get("OsVersion"),
        "config": image.get("Config"),
        "rootfs": image.get("RootFS"),
    }
    if not isinstance(stable["config"], dict) or not isinstance(stable["rootfs"], dict):
        return ""
    encoded = json.dumps(stable, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return "sha256:" + hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def reference_inventory(probes=()):
    from mn_molecular_design_studio.runtime import reference_root
    root = reference_root()
    result = {}
    # Probe named resources rather than repeatedly walking large model databases.
    paths = list(root.iterdir()) if root.is_dir() else []
    paths.extend(root / str(p) for p in probes if not Path(str(p)).is_absolute() and ".." not in Path(str(p)).parts)
    for path in paths:
        try:
            if path.is_file() and path.stat().st_size:
                result[path.relative_to(root).as_posix()] = path.stat().st_size
            elif path.is_dir() and next(path.iterdir(), None) is not None:
                result[path.relative_to(root).as_posix()] = True
        except OSError:
            continue
    return result


def resources():
    from mn_molecular_design_studio.engine.workflow_engine.core.central_scheduler import _client
    try:
        raw = _client(timeout=3).resources()
        reservations = raw.get("reservations", [])
        busy = {str(i) for r in reservations for i in r.get("gpu_devices", [])}
        return True, {"cpu_free": max(0, raw.get("cpu_pool_capacity", 0) - sum(int(r.get("cpu_threads", 0)) for r in reservations)),
                      "cpu_total": raw.get("cpu_pool_capacity", 0), "ram_available_gb": max(0, raw.get("ram_available_gb", 0) - sum(float(r.get("ram_gb", 0)) for r in reservations)),
                      "scratch_free_gb": max(0, raw.get("scratch_free_gb", 0) - sum(float(r.get("scratch_gb", 0)) for r in reservations)),
                      "free_gpus": [g for g in raw.get("gpus", []) if str(g["gpu_id"]) not in busy],
                      "gpus": raw.get("gpus", []), "reservations": len(reservations)}
    except Exception:
        return False, {}


def requirements(payload, *, references=(), installed_images=None):
    from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import WORKFLOWS
    from mn_molecular_design_studio.engine.workflow_engine.workflows.optimized_prediction import ENGINE_PROFILES
    from mn_molecular_design_studio.runtime import reference_root
    import yaml
    needed_images, needed_refs = set(), set(references)
    workflows = []
    if payload.get("graph"):
        from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG
        for node in payload["graph"]["nodes"]:
            spec = NODE_CATALOG[node["data"]["kind"]]
            if spec.get("workflow"):
                workflows.append((spec["workflow"], node["data"].get("params", {})))
    else:
        from mn_molecular_design_studio.batching import parameter_items
        workflows.extend((payload["workflow"], item) for item in parameter_items(payload.get("parameters", {})))
    cpu_only = True
    minimum_gpus = 0
    for workflow, p in workflows:
        if workflow not in WORKFLOWS:
            raise ValueError(f"Unknown workflow: {workflow}")
        spec = WORKFLOWS[workflow]
        manifest_path = Path(__file__).resolve().parents[1] / "engine/workflow_engine/manifests" / (spec.tool + ".yaml")
        if manifest_path.exists():
            manifest = yaml.safe_load(manifest_path.read_text())
            if manifest.get("image"):
                needed_images.add(manifest["image"])
            cpu_only = cpu_only and not manifest.get("gpu", False)
        if workflow == "refolding.optimized_prediction":
            profile = ENGINE_PROFILES[p.get("engine", "boltz2")]
            needed_images.add(profile["image"])
            directories = {"boltz2": "boltz_models", "esmfold": "esm_models", "esmfold2": "biohub-esm", "esmfold2_fast": "biohub-esm",
                           "protenix_v2": "protenix", "protenix_v1": "protenix", "alphafold3": "alphafold3", "openfold3_jax": "openfold3",
                           "af3_torch": "mn-cofolding/af3_torch", "rf3_foundry": "foundry", "chai1": "chai1", "openfold3_ob0": "openfold3",
                           "atlasfold": "atlasfold", "opendde_native": "opendde", "colabfold": "alphafold_models", "boltzgen_fold": "boltzgen-cache"}
            if p.get("engine", "boltz2") in directories:
                needed_refs.add(directories[p.get("engine", "boltz2")])
            if p.get("engine") == "boltz2":
                needed_refs.update("boltz_models/" + name for name in ("boltz2_conf.ckpt", "boltz2_aff.ckpt", "ccd.pkl"))
            elif p.get("engine") == "rf3_foundry":
                needed_refs.add("foundry/rf3_foundry_01_24_latest_remapped.ckpt")
            elif p.get("engine") == "openfold3_ob0":
                needed_refs.add("openfold3/of3-ob-2025-06-30-174k.pt")
            cpu_only = False
            minimum_gpus = max(minimum_gpus, int(profile.get("gpu_count", {}).get(p.get("optimization_mode", "off"), 1)))
        elif workflow == "sequence.mpnn":
            from mn_molecular_design_studio.mpnn_input import checkpoint
            needed_images.add("mn-foundry:cu128" if p.get("backend") == "foundry" else "mn-ligandmpnn:latest")
            if p.get("backend", "auto") == "auto":
                needed_images.add("mn-foundry:cu128")
            level = p.get("noise_level")
            path = checkpoint(p.get("model_type", "protein_mpnn"), None if level in (None, "") else float(level))
            needed_refs.add(path.relative_to(reference_root()).as_posix())
            cpu_only = False
        elif workflow.startswith("analysis."):
            from mn_molecular_design_studio.engine.workflow_engine.workflows.analysis import IPSAE_IMAGE, PYROSETTA_IMAGE
            from mn_molecular_design_studio.quality_tools import WORKFLOWS as quality_workflows
            if "ipsae" in workflow:
                needed_images.add(IPSAE_IMAGE)
            elif "rosetta" in workflow:
                needed_images.add(PYROSETTA_IMAGE)
            elif workflow in quality_workflows:
                needed_images.add(quality_workflows[workflow]["image"])
                needed_refs.update(quality_workflows[workflow].get("references", []))
                cpu_only = cpu_only and not quality_workflows[workflow].get("gpu", False)
        elif workflow.startswith("design.bindcraft2"):
            from mn_molecular_design_studio.engine.workflow_engine.workflows.bindcraft2 import BINDCRAFT2_IMAGE
            needed_images.add(p.get("docker_image") or BINDCRAFT2_IMAGE)
            needed_refs.add("alphafold_models")
            cpu_only = False
        elif workflow in {"sequence.esm_if1", "sequence.caliby", "sequence.boltzgen_inverse"}:
            from mn_molecular_design_studio.inverse_folding_input import ENGINES, weights
            engine = workflow.removeprefix("sequence.")
            needed_images.add(ENGINES[engine]["image"])
            needed_refs.update(path.relative_to(reference_root()).as_posix() for path in weights(engine, p))
            cpu_only = False
        elif workflow in {"design.rfantibody", "design.rfantibody_mimosa"}:
            from mn_molecular_design_studio.rfantibody_input import IMAGE, MIMOSA_IMAGE
            needed_images.add(MIMOSA_IMAGE if p.get("mimosa_enabled", False) else IMAGE)
            required_weights = ["RFdiffusion_Ab.pt"]
            if not p.get("generation_only", False):
                required_weights.append("RF2_ab.pt")
                if not p.get("use_abmpnn", False):
                    required_weights.append("ProteinMPNN_v48_noise_0.2.pt")
            needed_refs.update("rfantibody/" + name for name in required_weights)
            cpu_only = False
        elif workflow == "design.rfpeptides":
            from mn_molecular_design_studio.rfpeptides_input import WEIGHTS_RELATIVE_PATH
            needed_refs.add(WEIGHTS_RELATIVE_PATH.as_posix())
            cpu_only = False
        elif workflow == "design.rfdiffusion_beta_pairing":
            from mn_molecular_design_studio.engine.workflow_engine.workflows.rfdiffusion_beta_pairing import RF_IMAGE, TENSOR_IMAGE
            needed_images.update({RF_IMAGE, TENSOR_IMAGE})
            needed_refs.add("rfdiffusion_models/Complex_Fold_base_ckpt.pt")
            cpu_only = False
        else:
            if workflow.startswith(("design.", "refolding.", "sequence.")):
                cpu_only = False
            if "rfdiffusion3" in workflow:
                needed_images.add("mn-foundry:cu128"); needed_refs.add("foundry/rfd3_latest.ckpt")
            if workflow in {"design.rfd4_proteina","refolding.rfd4_proteina"}:
                from mn_molecular_design_studio.rfd4_proteina_input import IMAGE,HF_REVISION
                needed_images.add(IMAGE)
                bundle="rfd4-proteina/huggingface/hub/models--nvidia--RFD4-Proteina/snapshots/"+HF_REVISION
                needed_refs.update(bundle+"/"+name for name in ("config.json","model.safetensors"))
                o=p.get("options") or {}
                needs_oracle=workflow=="design.rfd4_proteina" and (o.get("refold",True) or o.get("apo") or o.get("reward_search")=="Best by refolding" or p.get("avoid_structures") or p.get("also_bind_structures"))
                if needs_oracle:
                    engine=o.get("folding_engine","openfold3_jax")
                    from mn_molecular_design_studio.engine.workflow_engine.workflows.rfd4_oracle import ORACLE_IMAGE
                    needed_images.add(ORACLE_IMAGE if engine=="openfold3_jax" else ENGINE_PROFILES[engine]["image"])
                    if engine=="openfold3_jax":
                        needed_refs.update(bundle+"/openfold3/weights/"+name for name in ("of3_ported_weights.bin.zst","of3_ported_weights.provenance.json"))
                    else:
                        other=requirements({"workflow":"refolding.optimized_prediction","parameters":{"engine":engine}},installed_images=installed_images)
                        needed_refs.update(other["references"])
                if o.get("cluster"): needed_images.add("mn-foldseek:10-941cd33-qc1")
            from mn_molecular_design_studio import generator_input
            if workflow in generator_input.WORKFLOWS.values():
                engine = next(k for k, v in generator_input.WORKFLOWS.items() if v == workflow)
                paths = generator_input.readiness(engine, p)["required"]
                needed_refs.update(Path(path).relative_to(reference_root()).as_posix() for path in paths)
            if "complexa" in workflow:
                from mn_molecular_design_studio.complexa_input import CHECKPOINTS
                needed_refs.update("proteina-complexa/ckpts/" + name for name in CHECKPOINTS[p.get("task", "protein_binder")])
            for name, directory in (("rfdiffusion", "rfdiffusion_models"), ("genie3", "genie3"), ("pxdesign", "pxdesign"), ("protpardelle", "protpardelle-1c"), ("complexa", "proteina-complexa"), ("boltzgen", "boltzgen-cache")):
                if name in workflow:
                    needed_refs.add(directory)
    installed = installed_images if installed_images is not None else images()
    requested = payload.get("resources", {})
    selected = payload.get("parameters", {})
    p = selected[0] if isinstance(selected, list) else selected
    default_device = "auto"
    if payload.get("graph"):
        default_device = "all"
    elif workflows:
        from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import _get_function
        parameter = inspect.signature(_get_function(WORKFLOWS[workflows[0][0]])).parameters.get("gpu_device")
        if parameter is not None and parameter.default is not inspect.Parameter.empty:
            default_device = parameter.default
    devices = {str(item.get("gpu_device", default_device)) for _, item in workflows}
    if len(devices) > 1 and not requested.get("gpu_device") and not payload.get("graph"):
        raise ValueError("Selected inputs must share one GPU allocation; set resources.gpu_device for the batch.")
    device = str(requested.get("gpu_device", payload.get("gpu_device", p.get("gpu_device", default_device)))).removeprefix("device=")
    gpu_ids = [int(v) for v in device.split(",")] if device and all(v.isdigit() for v in device.split(",")) else []
    if gpu_ids and len(gpu_ids) < minimum_gpus:
        raise ValueError(f"This engine mode requires at least {minimum_gpus} GPU devices.")
    gpu_count = len(gpu_ids) or (0 if cpu_only or device in {"cpu", "none"} else max(1, minimum_gpus))
    return {"images": {image: installed.get(image, "") for image in sorted(needed_images)}, "references": sorted(needed_refs),
            "code_revision": code_revision(), "gpu_count": gpu_count, "gpu_ids": gpu_ids,
            "all_gpus": not cpu_only and device == "all",
            "min_vram_gb": 8 if gpu_count else 0, "cpu_threads": requested.get("cpu_cores") or max((item.get("cpu_cores") or item.get("cpu_threads") or 4 for _, item in workflows), default=4),
            "ram_gb": 0, "scratch_gb": 1}
