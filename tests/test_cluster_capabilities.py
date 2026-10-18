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
