"""Scorer outputs must stay writable under rootful and rootless Docker."""
from types import SimpleNamespace
import pytest
from mn_molecular_design_studio.engine.workflow_engine.workflows import analysis, quality


@pytest.mark.parametrize('workflow', ['rosetta', 'ipsae'])
def test_saved_structure_scorers_honor_shared_container_user_mapping(tmp_path, monkeypatch, design_source_factory, workflow):
    run, _ = design_source_factory(missing_occupancy=False)
    monkeypatch.setenv('MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR', str(tmp_path / 'runs'))
    monkeypatch.setattr(quality, 'docker_user_args', lambda: ['--user','123:456'])
    commands = []
    def execute(command, **kwargs):
        commands.append(list(map(str, command)))
        return SimpleNamespace(returncode=0, stdout='', stderr='')
    monkeypatch.setattr(analysis.subprocess, 'run', execute)
    function = analysis.run_rosetta_metrics_contract if workflow == 'rosetta' else analysis.run_ipsae_analysis_contract
    function(run, run/'artifacts/normalized_candidates/candidates.jsonl')
    command = next(c for c in commands if c[:2] == ['docker','run'])
    assert command[command.index('--user')+1] == '123:456'
