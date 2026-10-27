from __future__ import annotations

import json

from mn_molecular_design_studio.cluster import capabilities


def _image(*, image_id: str, layer: str = "layer-a", command: list[str] | None = None, size: int = 10):
    return {
        "Id": image_id,
        "Size": size,
        "Os": "linux",
        "Architecture": "amd64",
        "Variant": "",
        "OsVersion": "",
        "Config": {
            "Entrypoint": command or ["/usr/bin/run"],
            "Env": ["PATH=/usr/bin"],
        },
        "RootFS": {"Type": "layers", "Layers": [layer]},
        "RepoTags": ["mn-ligandmpnn:latest"],
    }


def test_image_fingerprint_ignores_store_specific_id_and_size():
    classic = _image(image_id="sha256:classic-manifest", size=100)
    containerd = _image(image_id="sha256:containerd-config", size=200)

    assert capabilities.image_content_fingerprint(classic) == capabilities.image_content_fingerprint(containerd)


def test_image_fingerprint_changes_when_config_or_layer_changes():
    baseline = _image(image_id="sha256:one")

    assert capabilities.image_content_fingerprint(baseline) != capabilities.image_content_fingerprint(
        _image(image_id="sha256:two", layer="layer-b")
    )
    assert capabilities.image_content_fingerprint(baseline) != capabilities.image_content_fingerprint(
        _image(image_id="sha256:three", command=["/usr/bin/other"])
    )


def test_images_returns_portable_fingerprints_for_each_tag(monkeypatch):
    first = _image(image_id="sha256:classic-manifest")
    second = _image(image_id="sha256:containerd-config", layer="layer-b")
    first["RepoTags"] = ["mn-ligandmpnn:latest", "mn-ligandmpnn:backup"]
    second["RepoTags"] = ["mn-foundry:cu128"]

    def fake_run(command, **_kwargs):
        if command[:4] == ["docker", "image", "ls", "--format"]:
            return type("Completed", (), {"stdout": "mn-ligandmpnn:latest\nmn-foundry:cu128\n"})()
        assert command[:4] == ["docker", "image", "inspect", "--format"]
        return type(
            "Completed",
            (),
            {"stdout": "\n".join(json.dumps(item) for item in (first, second))},
        )()

    monkeypatch.setattr(capabilities.subprocess, "run", fake_run)

    inventory = capabilities.images()

    assert inventory == {
        "mn-ligandmpnn:latest": capabilities.image_content_fingerprint(first),
        "mn-ligandmpnn:backup": capabilities.image_content_fingerprint(first),
        "mn-foundry:cu128": capabilities.image_content_fingerprint(second),
    }


def test_rfantibody_generation_only_requires_only_diffusion_weights():
    from mn_molecular_design_studio.rfantibody_input import IMAGE, MIMOSA_IMAGE

    installed = {IMAGE: "native-image", MIMOSA_IMAGE: "mimosa-image"}
    for workflow, mimosa_enabled, image in (
        ("design.rfantibody", False, IMAGE),
        ("design.rfantibody_mimosa", True, MIMOSA_IMAGE),
    ):
        needs = capabilities.requirements(
            {
                "workflow": workflow,
                "parameters": {
                    "generation_only": True,
                    "mimosa_enabled": mimosa_enabled,
                    "gpu_device": "0",
                },
            },
            installed_images=installed,
        )

        assert needs["images"] == {image: installed[image]}
        assert "rfantibody/RFdiffusion_Ab.pt" in needs["references"]
        assert "rfantibody/ProteinMPNN_v48_noise_0.2.pt" not in needs["references"]
        assert "rfantibody/RF2_ab.pt" not in needs["references"]


def test_rfantibody_full_pipeline_retains_required_sequence_and_refold_weights():
    from mn_molecular_design_studio.rfantibody_input import IMAGE

    needs = capabilities.requirements(
        {
            "workflow": "design.rfantibody",
            "parameters": {"generation_only": False, "use_abmpnn": False, "gpu_device": "0"},
        },
        installed_images={IMAGE: "native-image"},
    )

    assert {
        "rfantibody/RFdiffusion_Ab.pt",
        "rfantibody/ProteinMPNN_v48_noise_0.2.pt",
        "rfantibody/RF2_ab.pt",
    }.issubset(needs["references"])


def test_reference_fingerprints_detect_changed_model_contents(tmp_path, monkeypatch):
    monkeypatch.setenv('MN_MOLECULAR_DESIGN_STUDIO_REFERENCE_DIR', str(tmp_path))
    path = tmp_path / 'model.ckpt'; path.write_bytes(b'weights-v1')
    before = capabilities.reference_fingerprints(['model.ckpt'])
    path.write_bytes(b'weights-v2')
    after = capabilities.reference_fingerprints(['model.ckpt'])
    assert before != after
    from mn_molecular_design_studio.cluster.store import Store
    node = {'availability': 'online', 'snapshot': {'accept_jobs': True, 'scheduler_online': True, 'reference_hashes': before}}
    assert Store.compatible(node, {'reference_hashes': before})
    assert not Store.compatible(node, {'reference_hashes': after})
