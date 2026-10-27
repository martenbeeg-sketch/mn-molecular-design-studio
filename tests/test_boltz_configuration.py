"""F4/F6 regressions at graph, request, staging and native invocation boundaries."""
import copy
from pathlib import Path

import pytest

from mn_molecular_design_studio import pipeline_graph as catalog
from mn_molecular_design_studio.engine.workflow_engine.workflows import graph_pipeline as graph
from mn_molecular_design_studio.engine.workflow_engine.workflows import optimized_prediction as prediction


def boltz_request(settings):
    return prediction.boltz2_prediction_parameters([{'protein_sequences': 'ACDEFG'}], settings=settings, job_name='test')


def test_boltz_alias_only_config_changes_effective_request():
    request = boltz_request({'diffusion_samples': 4, 'recycling_steps': 7})
    assert (request['samples'], request['num_recycles']) == (4, 7)


@pytest.mark.parametrize('settings', [
    {'samples': 1, 'diffusion_samples': 4}, {'num_recycles': 3, 'recycling_steps': 7},
])
def test_boltz_conflicts_are_rejected(settings):
    with pytest.raises(ValueError, match='[Cc]onflict'):
        boltz_request(settings)


def test_boltz_active_ui_has_one_control_per_quantity():
    fields = {f['name'] for f in catalog.NODE_CATALOG['boltz2']['fields']}
    assert {'samples', 'num_recycles'} <= fields
    assert not {'diffusion_samples', 'recycling_steps', 'write_full_pae'} & fields


def test_boltz_pae_requirement_and_effective_native_options():
    profile = prediction.ENGINE_PROFILES['boltz2']
    with pytest.raises(ValueError, match='PAE'):
        prediction._prepare_prediction('boltz2', protein_sequences='ACDEFG',
                                       engine_options={'write_full_pae': False, 'generate_ipsae': True})
    with pytest.raises(ValueError, match='PAE'):
        prediction._prepare_prediction('boltz2', protein_sequences='ACDEFG',
                                       engine_options={'write_full_pae': False, 'generate_ipsae': False})
    options = prediction._batch_options(profile, seed=1, num_recycles=7,
                                        engine_options={'write_full_pae': True, 'generate_ipsae': False})
    assert options['write_full_pae'] is True
    assert options['recycling_steps'] == 7


def test_saved_graph_shipped_alias_defaults_keep_historical_execution():
    from mn_molecular_design_studio.prediction_configuration import boltz_settings
    old = {'samples': 1, 'num_recycles': 3, 'diffusion_samples': 3, 'recycling_steps': 10, 'write_full_pae': False}
    migrated = boltz_settings(old, saved_graph=True)
    assert (migrated['samples'], migrated['num_recycles'], migrated['write_full_pae']) == (1, 3, True)
    assert len(migrated['configuration_migration']) == 3
    assert boltz_settings(migrated, saved_graph=True) == migrated
    assert old['write_full_pae'] is False


def test_saved_graph_changed_conflicting_alias_requires_resolution():
    from mn_molecular_design_studio.prediction_configuration import boltz_settings
    with pytest.raises(ValueError, match='Conflicting'):
        boltz_settings({'samples': 1, 'diffusion_samples': 5}, saved_graph=True)


def test_effective_settings_change_manifest_and_cache_fingerprint():
    from mn_molecular_design_studio.pipeline_graph import fingerprint
    a, b = boltz_request({'samples': 1, 'num_recycles': 1}), boltz_request({'samples': 2, 'num_recycles': 2})
    assert fingerprint('boltz2', a, []) != fingerprint('boltz2', b, [])
    options = [prediction._batch_options(prediction.ENGINE_PROFILES['boltz2'], seed=5,
               num_recycles=r['num_recycles'], engine_options=r['engine_options']) for r in (a, b)]
    assert options[0]['recycling_steps'] == 1 and options[1]['recycling_steps'] == 2


def test_graph_pae_requirement_follows_downstream_analysis():
    graph = catalog.new_graph('all_generators_boltz2')
    params = next(n['data']['params'] for n in graph['nodes'] if n['data']['kind'] == 'boltz2')
    params.update(write_full_pae=False, generate_ipsae=False, prediction_config_version=2)
    with pytest.raises(ValueError, match='PAE'):
        catalog.validate_graph(graph)


@pytest.mark.parametrize('settings', [{'samples': 0}, {'num_recycles': 0}])
def test_zero_sampling_or_recycles_cannot_silently_restore_defaults(settings):
    with pytest.raises(ValueError):
        prediction.validate_prediction_request(boltz_request(settings))
