from copy import deepcopy
import json

import pytest

from mn_molecular_design_studio.engine.workflow_engine.workflows.execution_planning import (
    PROPOSED_NATIVE_CHOICES, describe_native_invocation, plan_campaign_execution,
)
from mn_molecular_design_studio.engine.workflow_engine.workflows.distributed_campaign import validate_campaign_config


def test_legacy_config_is_read_only_and_retains_dispatch_product(monkeypatch, tmp_path):
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs, workflow_queue
    def forbidden(*a, **k):
        pytest.fail("Planning must not create, reserve or dispatch jobs")
    monkeypatch.setattr(jobs, "create_job", forbidden)
    monkeypatch.setattr(workflow_queue, "enqueue_workflow_call", forbidden)
    monkeypatch.setattr(workflow_queue, "spawn_worker_for_run", forbidden)
    config = {"max_attempts": 7, "accepted_target": 2, "batch_size": 2, "max_concurrent_batches": 3}
    original = deepcopy(config)
    plan = plan_campaign_execution(config)
    assert plan == plan_campaign_execution(config)
    assert config == original
    assert not list(tmp_path.iterdir())
    assert plan['read_only'] is True and plan['submission_authorized'] is False
    assert plan['execution']['max_concurrent_generation_jobs'] == 6
    assert plan['execution']['max_concurrent_validations'] == 2
    assert plan['supported_execution_plan']['execution_jobs_if_budget_exhausted'] == 7
    assert plan['native_attempts_per_job']['effective'] == 1
    assert validate_campaign_config(config) == validate_campaign_config(original)
    alias = {**config, 'dispatch_batch_size': config['batch_size']}
    assert plan_campaign_execution(alias) == plan
    json.dumps(plan, allow_nan=False)


@pytest.mark.parametrize('size', PROPOSED_NATIVE_CHOICES[1:])
def test_proposed_choices_have_no_supported_fallback(size):
    plan = plan_campaign_execution({'max_attempts': 12}, native_attempts_per_job=size)
    assert not plan['campaign_configuration_supported']
    assert plan['supported_execution_plan'] is None
    assert plan['native_attempts_per_job']['requested'] == size
    assert plan['native_attempts_per_job']['effective'] is None
    assert plan['hypothetical_requested_plan']['total_requested_attempts'] == 12
    with pytest.raises(ValueError, match='supports only 1'):
        validate_campaign_config({'native_attempts_per_job': size})


def test_global_budget_is_not_multiplied_by_input_count():
    for budget in range(1, 37):
        for inputs in range(1, 7):
            for size in PROPOSED_NATIVE_CHOICES:
                p = plan_campaign_execution({'max_attempts': budget}, input_count=inputs,
                                            native_attempts_per_job=size)['hypothetical_requested_plan']
                base, extra = divmod(budget, inputs)
                assert p['base_attempts_per_input'] == base
                assert p['first_inputs_with_one_extra'] == extra
                assert sum(g['execution_jobs']*g['native_attempts_requested_per_job'] for g in p['job_groups']) == budget
                assert p['execution_jobs_if_budget_exhausted'] == sum((base+(i<extra)+size-1)//size for i in range(inputs))
    p = plan_campaign_execution({'max_attempts': 12}, native_attempts_per_job=5)['hypothetical_requested_plan']
    assert [(g['native_attempts_requested_per_job'],g['execution_jobs']) for g in p['job_groups']] == [(2,1),(5,2)]


@pytest.mark.parametrize('replicates,inputs,batch,invocations,loader_batches', [(5,3,2,1,8),(100,2,16,1,13),(101,2,16,2,14),(200,1,16,2,14)])
def test_rfd4_replicates_and_wrapper_chunks_are_distinct_from_gpu_batch(replicates,inputs,batch,invocations,loader_batches):
    p = describe_native_invocation('rfd4_proteina', attempts_per_input=replicates, input_count=inputs, native_batch_size=batch)
    assert p['native_attempts'] == replicates*inputs
    assert p['rounding_extra_attempts'] == 0
    assert p['native_python_invocations'] == invocations
    assert p['loader_batches'] == loader_batches


@pytest.mark.parametrize('designs,batch,inputs,actual,remainder', [(5,2,3,18,1),(4,2,2,8,0),(1,16,2,2,0),(7,3,1,9,1)])
def test_rfd3_rounding_is_visible_and_exact_remainder_is_only_proposed(designs,batch,inputs,actual,remainder):
    p = describe_native_invocation('rfdiffusion3', attempts_per_input=designs, input_count=inputs, native_batch_size=batch)
    assert p['native_attempts'] == actual
    assert p['rounding_extra_attempts'] == actual-designs*inputs
    exact = p['exact_remainder_alternative']
    assert exact['native_attempts'] == designs*inputs
    assert exact['remainder_per_input'] == remainder
    assert exact['implemented'] is False
    assert exact['native_python_invocations'] == 1+bool(remainder)


def test_rfd3_budget_overrun_cannot_be_a_supported_plan():
    p = plan_campaign_execution({'max_attempts':7}, engine_id='rfdiffusion3', native_attempts_per_job=5, native_batch_size=2)
    assert p['hypothetical_requested_plan']['total_native_attempts_including_rounding'] == 8
    assert any('8 physical attempts' in issue for issue in p['errors'])
    assert p['supported_execution_plan'] is None


def test_other_engines_never_invent_native_acceptance_or_output_counts():
    for engine in ('rfdiffusion','rfdiffusion3','boltzgen','pxdesign','genie3'):
        p = plan_campaign_execution({}, engine_id=engine)
        assert not p['campaign_configuration_supported']
        assert p['native_attempts_per_job']['supported_campaign_choices'] == []
    assert describe_native_invocation('boltzgen')['native_attempts'] is None
    assert describe_native_invocation('bindcraft2')['nominal_native_outputs'] is None
    assert describe_native_invocation('rfdiffusion', attempts_per_input=5, input_count=2)['native_python_invocations'] == 2


def test_gpu_eligibility_needs_current_physical_and_compatibility_evidence():
    base = {'node_id':'local','gpu_id':'0','gpu_uuid':'physical-0','is_local':True,
            'online':True,'compatible':True,'total_vram_gb':24,'cpu_cores_available':64,'observed_at':'fixture', 'snapshot_current':True}
    inventory = [base, {**base,'gpu_id':'1','online':False}, {**base,'gpu_id':'2','compatible':None},
                 {**base,'gpu_id':'3','snapshot_current':False}, {**base,'gpu_id':'4','total_vram_gb':8}]
    before = deepcopy(inventory)
    p = plan_campaign_execution({'resources':{'min_vram_gb':16}}, eligible_gpus=inventory)
    assert inventory == before
    assert [g['eligibility'] for g in p['gpu_inventory']] == ['eligible_in_supplied_snapshot','excluded','unresolved','unresolved','excluded']
    assert p['inventory_status'] == 'supplied_snapshot_only'
    assert plan_campaign_execution({})['inventory_status'] == 'unknown_not_probed'


@pytest.mark.parametrize('kwargs', [{'input_count':True},{'native_attempts_per_job':'auto'},{'native_batch_size':0},{'native_attempts_per_job':1.0},{'eligible_gpus':{}}])
def test_invalid_plan_shapes_raise_without_dispatch(kwargs):
    with pytest.raises(ValueError):
        plan_campaign_execution({}, **kwargs)


def test_proposed_concurrency_cannot_silently_override_execution():
    p = plan_campaign_execution({'batch_size':2,'max_concurrent_batches':3,'max_concurrent_generation_jobs':4})
    assert p['execution']['max_concurrent_generation_jobs'] == 6
    assert p['execution']['requested_proposed_generation_jobs'] == 4
    assert p['supported_execution_plan'] is None


def test_preview_controls_cannot_change_saved_campaign_or_dispatch():
    from streamlit.testing.v1 import AppTest
    app = AppTest.from_string(
        "import streamlit as st\n"
        "from mn_molecular_design_studio.distributed_campaign_config_ui import render_campaign_config_editor\n"
        "config, ready = render_campaign_config_editor(st,'bindcraft2',{'enabled':True,'max_attempts':7},key='preview')\n"
        "st.session_state['actual_campaign'] = config\n"
        "st.session_state['ready'] = ready\n"
    ).run()
    assert not app.exception
    original = deepcopy(app.session_state['actual_campaign'])
    app.selectbox(key='preview-plan-native-count').set_value(25).run()
    assert not app.exception
    assert app.session_state['ready'] is True
    assert app.session_state['actual_campaign'] == original
    assert 'native_attempts_per_job' not in original
    plan = json.loads(app.json[-1].value)
    assert plan['native_attempts_per_job']['requested'] == 25
    assert plan['supported_execution_plan'] is None


def test_resources_resolve_adapter_parent_and_validation_defaults():
    config = {'generator_id':'rfd4_proteina','acceptance_level':'validated',
              'validation_stages':[{'workflow':'refolding.optimized_prediction','required_metrics':['iptm']}]}
    plan = plan_campaign_execution(config)
    assert plan['campaign_configuration_supported']
    assert plan['execution']['generation_resources']['cpu_cores'] == 8
    assert plan['execution']['generation_resources']['min_vram_gb'] == 16
    assert plan['execution']['validation_resources'][0]['cpu_cores'] == 4
    assert plan_campaign_execution(config,generation_cpu_cores=12)['execution']['generation_resources']['cpu_cores'] == 12
    assert plan_campaign_execution({})['execution']['generation_resources']['cpu_cores'] == 4


@pytest.mark.parametrize('config',[{'resources':{'min_vram_gb':'16'}},{'resources':{'min_vram_gb':None}},{'resources':{'cpu_cores':True}},
                                   {'validation_stages':[{'resources':[]}]}])
def test_malformed_resource_requests_have_explanations(config):
    with pytest.raises(ValueError,match='resources'):
        plan_campaign_execution(config,eligible_gpus=[])
