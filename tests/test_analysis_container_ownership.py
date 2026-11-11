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


@pytest.mark.parametrize('legacy_percent', [70.0, 171.0])
def test_rosetta_excludes_mislabeled_sasa_ratio_even_when_within_percentage_range(tmp_path, monkeypatch, design_source_factory, legacy_percent):
    import csv
    import json
    from pathlib import Path
    run, _ = design_source_factory(missing_occupancy=False)
    monkeypatch.setenv('MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR', str(tmp_path/'runs'))
    monkeypatch.setattr(quality, 'docker_user_args', lambda: [])
    def execute(command, **kwargs):
        folder = Path(command[command.index('--folder')+1].removeprefix('rosetta:'))
        output = Path(command[command.index('--out-csv')+1])
        with output.open('w',newline='') as handle:
            writer = csv.writer(handle)
            writer.writerow(['binder_id','rosetta_rosetta_interface_dG','rosetta_rosetta_interface_percent_of_binder_in_intf'])
            writer.writerow([next(folder.glob('*.pdb')).stem,-12.5,legacy_percent])
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(analysis.subprocess, 'run', execute)
    job = analysis.run_rosetta_metrics_contract(run, run/'artifacts/normalized_candidates/candidates.jsonl')
    row = json.loads((job/'artifacts/normalized_candidates/candidates.jsonl').read_text().splitlines()[0])
    assert row['metrics']['rosetta_interface_dG'] == -12.5
    assert not any('percent_of_binder_in_intf' in key for key in row['metrics'])
    excluded = row['raw_metadata']['measurement_provenance']['excluded_rosetta_interface_fraction']
    assert excluded['raw_value'] == legacy_percent
    assert excluded['status'] == 'excluded'


def test_historical_rosetta_sasa_ratio_is_preserved_raw_but_not_projected_or_filterable(tmp_path, monkeypatch):
    from mn_molecular_design_studio import molecule_imports
    from mn_molecular_design_studio import evaluations
    from mn_molecular_design_studio.molecular_database import MolecularDatabase
    monkeypatch.setenv('MN_MOLECULAR_DESIGN_STUDIO_HOME', str(tmp_path))
    monkeypatch.setenv('MN_MOLECULAR_DESIGN_STUDIO_APP_HOME', str(tmp_path))
    database = MolecularDatabase(tmp_path/'database')
    A = {'id':'a'*32, 'role':'user', '_cluster_project_ids':[]}
    imported = database.import_records(A, molecule_imports.parse_sequences('ACDEFG'))
    molecule_id = imported['molecule_ids'][0]
    bad_key = 'rosetta_interface_percent_of_binder_in_intf'
    observation_id = database.add_observation(imported['group_id'], molecule_id, 'rosetta_metrics', 'old-rosetta',
                                             {bad_key:171.0, 'rosetta_interface_dG':-12.5})
    for _ in range(2):
        row = database.list_molecules(A, group_id=imported['group_id'])['items'][0]
        assert row['scores']['rosetta_metrics.rosetta_interface_dG'] == -12.5
        assert not any(bad_key in field for field in row['scores'])
        assert row['observations'][0]['metrics'][bad_key] == 171.0
        assert not any(bad_key in field for field in database.metric_fields(A, group_id=imported['group_id']))
        # Simulate an index produced by the older projection version.
        with database.connection() as connection:
            connection.execute('UPDATE molecule_evaluations SET projection_version=? WHERE observation_id=?',
                               ('evaluation.v4:'+evaluations.VERSION, observation_id))
            connection.execute('INSERT OR REPLACE INTO molecule_metric_values VALUES (?,?,?,?,?,?)',
                               (observation_id, 'rosetta_metrics.'+bad_key, 'number', 171.0, '171.0', '{}'))
        with monkeypatch.context() as isolated:
            normalize = evaluations.normalize_metrics
            calls = []
            def reproject_metrics(*args):
                calls.append(args)
                return normalize(*args)
            isolated.setattr(evaluations, 'normalize_metrics', reproject_metrics)
            with database.connection() as connection:
                evaluations.rebuild_missing(connection)
                assert connection.execute('SELECT projection_version FROM molecule_evaluations WHERE observation_id=?',
                                          (observation_id,)).fetchone()[0] == evaluations.PROJECTION_VERSION
                assert len(calls) == 1
                assert bad_key not in calls[0][0]
