from __future__ import annotations

import json

from mn_molecular_design_studio.job_catalog import set_job_hidden


def test_results_visibility_changes_without_removing_run_files(tmp_path, monkeypatch):
    run_root = tmp_path / "runs"
    run_dir = run_root / "vina" / "run-1"
    run_dir.mkdir(parents=True)
    metadata_path = run_dir / "metadata.json"
    metadata_path.write_text(json.dumps({"tool": "vina", "status": "completed"}))
    output_path = run_dir / "result.json"
    output_path.write_text('{"success": true}')
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR", str(run_root))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_HOME", str(tmp_path / "studio"))

    set_job_hidden("vina", "run-1", True)
    assert json.loads(metadata_path.read_text())["hidden"] is True
    assert output_path.read_text() == '{"success": true}'

    set_job_hidden("vina", "run-1", False)
    assert json.loads(metadata_path.read_text())["hidden"] is False
