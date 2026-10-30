"""MD metadata actions preserve the normal source-job permissions and receipts."""
import json

import pytest

from mn_molecular_design_studio.cluster import config
from mn_molecular_design_studio.cluster.agent import Agent
from mn_molecular_design_studio.cluster.store import Store


def test_md_source_action_requires_visibility_and_keeps_requesting_actor(monkeypatch, tmp_path):
    monkeypatch.setenv('MDS_CLUSTER_HOME', str(tmp_path / 'cluster'))
    store = Store()
    node = store.enroll('owner workstation')
    owner = {'id': 'actor', 'role': 'user', 'project_ids': ['project']}
    job = {'task_group': 'ligand-md', 'run_id': 'md-1', 'owner_id': 'actor', 'project_id': 'project',
           'status': 'completed', 'tool': 'openmm_md'}
    store.heartbeat(node['node_id'], {}, [job])
    jid = config.job_id(node['node_id'], 'ligand-md/md-1')
    with pytest.raises(LookupError):
        store.action(jid, 'md_source_options', {}, {'id': 'outsider', 'role': 'user', 'project_ids': []})
    operation = store.action(jid, 'md_source_options', {}, owner)
    with store.connect() as db:
        row = db.execute('SELECT * FROM commands WHERE id=?', (operation['command_id'],)).fetchone()
    assert row['kind'] == 'md_source_options'
    assert row['owner_id'] == owner['id']
    assert json.loads(row['payload'])['local_reference'] == 'ligand-md/md-1'


def test_accepted_design_collection_action_requires_an_explicit_acceptance(monkeypatch, tmp_path):
    monkeypatch.setenv('MDS_CLUSTER_HOME', str(tmp_path / 'cluster'))
    store = Store()
    node = store.enroll('design workstation')
    owner = {'id': 'actor', 'role': 'user', 'project_ids': ['project']}
    job = {
        'task_group': 'refolding-validation', 'run_id': 'accepted-validation',
        'owner_id': 'actor', 'project_id': 'project', 'status': 'completed',
        'tool': 'optimized_prediction_boltz2',
        'candidate_summary': {'acceptance_counts': {'accepted': 1, 'rejected': 2, 'unresolved': 0, 'not_evaluated': 0}},
    }
    store.heartbeat(node['node_id'], {}, [job])
    jid = config.job_id(node['node_id'], 'refolding-validation/accepted-validation')

    operation = store.action(jid, 'collect_accepted_designs', {}, owner)
    with store.connect() as db:
        command = db.execute('SELECT * FROM commands WHERE id=?', (operation['command_id'],)).fetchone()
    assert command['kind'] == 'collect_accepted_designs'
    assert json.loads(command['payload'])['local_reference'] == 'refolding-validation/accepted-validation'

    no_acceptance = {**job, 'run_id': 'no-acceptance',
                     'candidate_summary': {'acceptance_counts': {'accepted': 0}}}
    store.heartbeat(node['node_id'], {}, [no_acceptance])
    no_acceptance_id = config.job_id(node['node_id'], 'refolding-validation/no-acceptance')
    with pytest.raises(ValueError, match='no explicitly accepted'):
        store.action(no_acceptance_id, 'collect_accepted_designs', {}, owner)


def test_agent_packages_only_accepted_results_and_replays_collection_receipt(monkeypatch, tmp_path):
    import hashlib
    from mn_molecular_design_studio import accepted_design_import
    from mn_molecular_design_studio.engine.workflow_engine.services import local_automation

    monkeypatch.setenv('MDS_CLUSTER_HOME', str(tmp_path / 'agent-cluster'))
    agent = object.__new__(Agent)
    agent.node_id = 'worker-node'
    payload = b'accepted-only-archive'
    calls = []
    monkeypatch.setattr(
        local_automation, 'resolve_job',
        lambda reference: calls.append(('resolve', reference)) or {'run_dir': str(tmp_path / 'validation'), 'run_id': 'validation'},
    )
    monkeypatch.setattr(
        accepted_design_import, 'build_accepted_design_bundle',
        lambda run_dir, **kwargs: (payload, {'candidate_count': 1, 'accepted_candidate_ids': ['accepted-1']}),
    )
    monkeypatch.setattr(
        'mn_molecular_design_studio.cluster.transport.node_request',
        lambda path, **kwargs: {'blob_id': hashlib.sha256(payload).hexdigest(), 'size_bytes': len(payload)},
    )
    command = {
        'id': 'accepted-collection-command', 'kind': 'collect_accepted_designs',
        'job_id': 'validation-job', 'payload': {'local_reference': 'refolding-validation/validation'},
    }

    first = agent.execute(command)
    second = agent.execute(command)

    assert first == second and first[1] is None
    assert first[0]['candidate_count'] == 1
    assert calls == [('resolve', 'refolding-validation/validation')]


def test_md_source_metadata_receipt_is_replayed_without_repeating_work(monkeypatch, tmp_path):
    from mn_molecular_design_studio import auth, md_analysis_ui
    from mn_molecular_design_studio.engine.workflow_engine.services import local_automation
    monkeypatch.setenv('MDS_CLUSTER_HOME', str(tmp_path / 'cluster'))
    actor = {'id': 'actor', 'role': 'user', 'disabled': False}
    monkeypatch.setattr(auth, 'list_users', lambda: [actor])
    monkeypatch.setattr(local_automation, 'resolve_job', lambda reference: {'run_id': 'md-1'})
    calls = []
    def options(job, user):
        calls.append((job, user))
        return {'source_options': [{'trajectory_relative_path': 'artifacts/legacy/one/production.dcd'}]}
    monkeypatch.setattr(md_analysis_ui, 'md_source_options_for_job', options)
    agent = object.__new__(Agent)
    command = {'id': 'metadata-operation', 'kind': 'md_source_options', 'job_id': 'owner/job',
               'owner_id': 'actor', 'payload': {'local_reference': 'ligand-md/md-1'}}
    first = agent.execute(command)
    second = agent.execute(command)
    assert first == second
    assert first[1] is None
    assert calls == [({'run_id': 'md-1'}, actor)]


def test_disabled_actor_cannot_resolve_new_source_metadata(monkeypatch, tmp_path):
    from mn_molecular_design_studio import auth, md_analysis_ui
    monkeypatch.setenv('MDS_CLUSTER_HOME', str(tmp_path / 'cluster'))
    monkeypatch.setattr(auth, 'list_users', lambda: [{'id': 'actor', 'disabled': True}])
    monkeypatch.setattr(md_analysis_ui, 'md_source_options_for_job', lambda *args: pytest.fail('Disabled actor reached source loader'))
    response, error = object.__new__(Agent).execute({'id': 'disabled-operation', 'kind': 'md_source_options',
        'job_id': 'owner/job', 'owner_id': 'actor', 'payload': {'local_reference': 'ligand-md/md-1'}})
    assert response == {}
    assert 'account is unavailable' in error
