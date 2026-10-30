from datetime import datetime, timedelta, timezone

import pytest

from mn_molecular_design_studio.campaign_objectives import normalize_objective, objective_progress


def test_timer_has_no_invented_count_target_and_recovery_keeps_original_start():
    config = normalize_objective({'objective_type': 'time', 'max_runtime_seconds': 600})
    created = datetime(2026, 10, 30, 7, tzinfo=timezone.utc)
    first = objective_progress(config, 2000, created.isoformat(), now=created + timedelta(seconds=599))
    assert first['objective_target'] is None
    assert first['objective_count'] == 2000
    assert first['objective_met'] is False
    assert first['remaining_seconds'] == 1
    recovered = objective_progress(config, 2001, created.isoformat(), now=created + timedelta(seconds=600))
    assert recovered['objective_met'] is True
    assert recovered['count_target_met'] is False


def test_duration_and_rome_offset_deadline_use_the_earliest_stop():
    config = normalize_objective({'objective_type': 'time', 'max_runtime_seconds': 3600,
                                  'deadline_at': '2026-10-30T08:30:00+01:00'})
    assert config['deadline_at'] == '2026-10-30T07:30:00+00:00'
    progress = objective_progress(config, 0, '2026-10-30T07:00:00Z',
                                 now=datetime(2026, 10, 30, 7, 30, tzinfo=timezone.utc))
    assert progress['time_expired']
    assert progress['remaining_seconds'] == 0


@pytest.mark.parametrize('bad', [
    {'objective_type': 'time'},
    {'objective_type': 'time', 'max_runtime_seconds': 10, 'accepted_target': 1},
    {'objective_type': 'time', 'max_runtime_seconds': float('nan')},
    {'objective_type': 'generated', 'deadline_at': '2026-10-30 08:00'},
    {'objective_type': 'time', 'max_runtime_seconds': True},
])
def test_invalid_or_ambiguous_objectives_are_rejected(bad):
    with pytest.raises(ValueError):
        normalize_objective(bad)


def test_invalid_deadline_draft_stays_editable_and_cannot_be_submitted():
    from streamlit.testing.v1 import AppTest
    app = AppTest.from_string(
        "import streamlit as st\n"
        "from mn_molecular_design_studio.campaign_objectives_ui import render_time_limits\n"
        "limits, ready = render_time_limits(st, {'deadline_at': 'bad date'}, key='test', required=True)\n"
        "st.session_state['ready'] = ready\n"
        "st.session_state['limits'] = limits\n"
    ).run()
    assert not app.exception
    assert app.session_state['ready'] is False
    assert app.session_state['limits']['deadline_at'] == 'bad date'
    app.text_input(key='test-deadline').set_value('2026-10-31 08:00').run()
    assert not app.exception
    assert app.session_state['ready'] is True
    assert app.session_state['limits']['deadline_at'] == '2026-10-31T07:00:00+00:00'


def test_expired_shared_scope_reconciles_pending_batches_without_hanging_or_dispatch(tmp_path):
    from types import SimpleNamespace
    from mn_molecular_design_studio import candidate_execution
    from mn_molecular_design_studio.engine.workflow_engine.workflows.shared_campaign import SharedCampaignLedger
    shared = SharedCampaignLedger.create(tmp_path / 'shared.json', {
        'objective_type': 'time', 'deadline_at': '2000-01-01T00:00:00Z', 'max_attempts': 1000,
        'generator_nodes': ['generator'], 'max_concurrent_jobs': 4, 'max_concurrent_generation_jobs': 4,
    })
    def local(_parent, _node_id, rows, target=None):
        return {'candidate_count': len(list(rows)), 'target': target or {}}
    graph = SimpleNamespace(_local_source=local, _absolute_candidates=lambda s: [],
                            node_request=lambda k, p, s, n: {'workflow': 'design.rfd4_proteina', 'parameters': p},
                            NODE_CATALOG={'rfd4': {'category': 'Design engine'}})
    node = {'id': 'generator', 'data': {'kind': 'rfd4', 'label': 'RFD4'}}
    state = {}
    output, status = candidate_execution.tick(
        tmp_path, node, state, {}, {'num_designs': 1000}, 'a' * 64,
        candidate_execution.policy({'mode': 'sharded', 'batch_size': 10, 'max_in_flight': 4}),
        graph=graph, owner_id='owner', project_id='project', gpu_device='auto', persist=lambda: None,
        shared_campaign_ledger=shared, generator_node='generator',
    )
    assert len(state['batches']) == 100
    assert all(batch['unit_count'] == 10 for batch in state['batches'])
    assert all(batch['objective_stop'] and not batch['attempts'] for batch in state['batches'])
    assert status == 'completed' and output['candidate_count'] == 0
