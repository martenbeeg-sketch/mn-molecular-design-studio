import json
import subprocess

import pytest

from mn_molecular_design_studio.engine.workflow_engine.core import execution_timing, docker_runner


def test_spans_pair_across_retries_and_record_nonzero_exit(tmp_path, monkeypatch):
    monkeypatch.setattr(docker_runner, "build_docker_command", lambda *a, **k: ["docker", "run"])
    monkeypatch.setattr(docker_runner, "update_status", lambda *a, **k: None)
    monkeypatch.setattr(docker_runner.subprocess, "run", lambda *a, **k: subprocess.CompletedProcess(a, 7))
    for _ in range(2):
        assert docker_runner.run_docker_job(tmp_path, {}, []) == 7
    rows = [json.loads(line) for line in (tmp_path / "execution_lifecycle.jsonl").read_text().splitlines()]
    assert len(rows) == 4
    assert rows[0]["span_id"] == rows[1]["span_id"] != rows[2]["span_id"]
    assert [row["event"] for row in rows] == ["started", "finished"] * 2
    for row in rows[1::2]:
        assert row["return_code"] == 7
        assert row["elapsed_seconds"] >= 0
        assert row["phase"] == "docker_command"


def test_launch_exception_is_preserved_and_has_terminal_timing(tmp_path):
    with pytest.raises(FileNotFoundError):
        with execution_timing.execution_phase(tmp_path, "docker_command"):
            raise FileNotFoundError("private command contents")
    output = (tmp_path / "execution_lifecycle.jsonl").read_text()
    row = json.loads(output.splitlines()[-1])
    assert row["exception_type"] == "FileNotFoundError"
    assert "return_code" not in row
    assert "private command contents" not in output


def test_telemetry_failure_does_not_fail_work_or_mask_exception(tmp_path, monkeypatch):
    def unavailable(*a, **k):
        raise OSError("disk full")
    monkeypatch.setattr(execution_timing.os, "open", unavailable)
    with execution_timing.execution_phase(tmp_path, "candidate_normalization"):
        pass
    with pytest.raises(ValueError, match="science"):
        with execution_timing.execution_phase(tmp_path, "candidate_normalization"):
            raise ValueError("science")


def test_shared_steps_keep_early_failure_and_distinguish_shell_wrappers(tmp_path, monkeypatch):
    from mn_molecular_design_studio.engine.workflow_engine.workflows import design

    monkeypatch.setattr(design, "apply_docker_cpu_limits_to_steps", lambda run, steps: steps)
    monkeypatch.setattr(design, "update_status", lambda *a, **k: None)
    calls = []
    def run(command, **kwargs):
        calls.append(command)
        return subprocess.CompletedProcess(command, 0 if len(calls) == 1 else 3)
    monkeypatch.setattr(design.subprocess, "run", run)
    assert design._run_shell_steps(tmp_path, [
        {"command": ["docker", "run", "image"]},
        {"command": ["bash", "-c", "false"]},
        {"command": ["docker", "run", "never-started"]},
    ]) == 3
    assert len(calls) == 2
    rows = [json.loads(line) for line in (tmp_path / "execution_lifecycle.jsonl").read_text().splitlines()]
    assert [row["phase"] for row in rows[1::2]] == ["docker_command", "execution_command"]
    assert [row["return_code"] for row in rows[1::2]] == [0, 3]
