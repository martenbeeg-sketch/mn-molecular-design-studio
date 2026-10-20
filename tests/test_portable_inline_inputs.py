from __future__ import annotations

from pathlib import Path

import pytest
from fastapi import HTTPException

from mn_molecular_design_studio import api
from mn_molecular_design_studio.engine.workflow_engine.core import jobs, portable_paths, workflow_queue


@pytest.mark.parametrize(
    ("key", "content"),
    [
        ("design_yaml", "entities:\n  - protein:\n      id: A\n      sequence: 60..60\n"),
        ("target_msa_a3m", ">query\nACDE\n>hit\nAC-E\n"),
        ("target_pdb_text", "HEADER    target.pdb\nATOM      1\nEND\n"),
        ("msa_a3m_text", ">query\nACDE\n"),
        ("target_sequence", "ACDEFG"),
        ("bias_by_res_json", "{}"),
        ("input_smiles", "C/C=C\\C"),
        # Inline fields must remain literal even when they resemble a managed URI.
        ("design_yaml", "runs:///literal.yaml"),
        ("target_pdb_text", "/reference_files/literal.pdb"),
    ],
)
def test_inline_inputs_survive_job_and_worker_serialization(monkeypatch, tmp_path, key, content):
    runs_root = tmp_path / "runs"
    run_dir = runs_root / "design" / "test-run"
    target = run_dir / "target.pdb"
    target.parent.mkdir(parents=True)
    target.write_text("END\n")
    monkeypatch.setattr(jobs, "runs_root", lambda: runs_root)
    monkeypatch.setattr(portable_paths, "runs_root", lambda: runs_root)

    kwargs = {key: content, "target_pdb": target}
    request_path = run_dir / "worker_request.json"
    jobs.write_json(request_path, {"kwargs": workflow_queue._invocation_value(kwargs)})
    restored = workflow_queue._restore_invocation_value(jobs.read_json(request_path)["kwargs"])

    assert restored == kwargs
    assert not portable_paths.is_operational_path_key(key)
    assert not portable_paths.is_operational_path_key("output_format")
    assert not portable_paths.is_operational_path_key("input_format")
    assert portable_paths.is_operational_path_key("target_pdb")
    assert portable_paths.is_operational_path_key("output_directory")
    assert portable_paths.is_operational_path_key("design_yaml_file")


def test_non_admin_path_guard_distinguishes_inline_yaml_and_real_paths(monkeypatch, tmp_path):
    user = {"id": "test-user", "role": "user"}
    monkeypatch.setattr(api, "user_input_root", lambda _user_id: tmp_path / "uploads")
    monkeypatch.setattr(api, "reference_root", lambda: tmp_path / "references")
    monkeypatch.setattr(api, "accessible_project_ids", lambda _user: set())
    monkeypatch.setattr(api, "path_is_within_accessible_job", lambda *args, **kwargs: False)

    api._validate_user_path_values(
        {
            "design_yaml": "entities:\n  - file:\n      path: {{TARGET_FILE}}\n",
            "target_msa_a3m": ">target/1-4\nACDE\n",
            "target_pdb_text": "HEADER /inline/target.pdb\nEND\n",
            "bias_by_res_json": "{}",
        },
        user,
        tmp_path,
    )
    with pytest.raises(HTTPException) as error:
        api._validate_user_path_values({"target_pdb": "/not-allowed/target.pdb"}, user, tmp_path)
    assert error.value.status_code == 403
