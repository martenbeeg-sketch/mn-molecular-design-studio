import hashlib
import json
from types import SimpleNamespace

import pytest

from mn_molecular_design_studio.engine.workflow_engine.workflows import ligand_preparation_benchmark as prep


def test_preparation_rejects_a_different_sample_before_creating_a_job(tmp_path, monkeypatch):
    source = tmp_path / "input.smi"
    source.write_text("CCO compound1\n")
    monkeypatch.setattr(prep, "create_job", lambda *a, **k: pytest.fail("A mismatched sample must not create a job"))
    with pytest.raises(ValueError, match="checksum"):
        prep.run_ligand_preparation_benchmark(source, expected_sha256="0" * 64, expected_count=1)


def test_preparation_rejects_duplicate_compound_ids(tmp_path, monkeypatch):
    source = tmp_path / "input.smi"
    source.write_text("CCO duplicate\nCC duplicate\n")
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    monkeypatch.setattr(prep, "create_job", lambda *a, **k: pytest.fail("Duplicate IDs must not create a job"))
    with pytest.raises(ValueError, match="unique"):
        prep.run_ligand_preparation_benchmark(source, expected_sha256=digest, expected_count=2)


def test_preparation_keeps_partial_success_without_retrying(tmp_path, monkeypatch):
    source = tmp_path / "input.smi"
    source.write_text("CCO compound1\nCC compound2\n")
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    run_dir = tmp_path / "job"
    (run_dir / "artifacts").mkdir(parents=True)
    monkeypatch.setattr(prep, "create_job", lambda *a, **k: SimpleNamespace(run_dir=run_dir))
    finishes = []
    monkeypatch.setattr(prep, "finish_job", lambda *args: finishes.append(args))

    def fake_docker(*args, **kwargs):
        assert kwargs["gpu_device"] == "none"
        assert set(kwargs["environment"].values()) == {"1"}
        scratch = kwargs["extra_mounts"][0][0]
        config = json.loads((scratch / "config.json").read_text())
        route = config["route"]["arguments"]
        assert route[route.index("--cpu") + 1] == "23"
        (scratch / "summary.json").write_text(json.dumps({
            "input_count": 2, "output_pdbqt_count": 1, "failed_input_count": 1, "retry_count": 0,
        }))
        (scratch / "prepared-pdbqt.tar.gz").write_bytes(b"test archive")
        return 0

    monkeypatch.setattr(prep, "run_docker_job", fake_docker)
    assert prep.run_ligand_preparation_benchmark(source, digest, expected_count=2) == run_dir
    assert finishes[-1][1] is True
    assert finishes[-1][2]["metrics"]["failed_input_count"] == 1
    assert finishes[-1][2]["metrics"]["retry_count"] == 0
