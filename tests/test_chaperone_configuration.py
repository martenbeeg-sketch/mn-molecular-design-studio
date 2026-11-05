from __future__ import annotations

import json

from mn_molecular_design_studio import chaperone_adapter


def test_managed_worker_configuration_needs_no_python_path_environment(monkeypatch, tmp_path):
    path = tmp_path / "chaperone-worker.json"
    path.write_text(json.dumps({
        "python": "/opt/studio/chaperone-worker/bin/python",
        "source": "/opt/chaperone/src",
        "blocks": "/opt/chaperone/src/blocks",
        "source_revision": "abc123",
    }), encoding="utf-8")
    monkeypatch.setenv("MDS_CHAPERONE_CONFIG", str(path))
    for name in ("MDS_CHAPERONE_PYTHON", "MDS_CHAPERONE_SRC", "MDS_CHAPERONE_BLOCKS"):
        monkeypatch.delenv(name, raising=False)

    assert chaperone_adapter.worker_configuration() == {
        "python": "/opt/studio/chaperone-worker/bin/python",
        "source": "/opt/chaperone/src",
        "blocks": "/opt/chaperone/src/blocks",
        "source_revision": "abc123",
    }


def test_explicit_worker_environment_overrides_managed_configuration(monkeypatch, tmp_path):
    path = tmp_path / "chaperone-worker.json"
    path.write_text(json.dumps({"python": "/saved/python", "source_revision": "abc123"}), encoding="utf-8")
    monkeypatch.setenv("MDS_CHAPERONE_CONFIG", str(path))
    monkeypatch.setenv("MDS_CHAPERONE_PYTHON", "/override/python")

    assert chaperone_adapter.worker_configuration()["python"] == "/override/python"
