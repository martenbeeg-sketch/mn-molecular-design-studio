import hashlib
import json
import sys
from types import ModuleType, SimpleNamespace

import pytest

from mn_molecular_design_studio.engine.workflow_engine.workflows import rfd4_native_diagnostics as diagnostics
from mn_molecular_design_studio import rfd4_proteina_input


@pytest.fixture
def fake_native(monkeypatch):
    calls = []
    result = object()
    def passthrough(*args, **kwargs):
        calls.append((args, kwargs))
        return result
    class Model:
        predict_step = passthrough
        _save_and_log_csv = passthrough
        reserve_shape_cache = passthrough
    class Trainer:
        def __init__(self, *args, **kwargs):
            calls.append((args, kwargs))
        validate = passthrough
    class MetricsPool:
        drain_generated_structure = passthrough
    metrics_module = ModuleType('rfproteina.metrics.utils.metrics_utils')
    metrics_module.MetricsPool = MetricsPool
    monkeypatch.setitem(sys.modules, metrics_module.__name__, metrics_module)
    native = SimpleNamespace(load_ckpt_n_configure_inference=passthrough, get_dataset_from_cfg=passthrough,
                             setup_run_compile=passthrough, Proteina=Model, L=SimpleNamespace(Trainer=Trainer),
                             torch=SimpleNamespace(cuda=SimpleNamespace(_lazy_init=passthrough, is_initialized=lambda:False)))
    return native, result, calls


def test_hooks_preserve_native_argument_identity_results_and_restore(fake_native, tmp_path):
    native, expected, calls = fake_native
    original = native.Proteina.predict_step
    path = tmp_path/'timings.jsonl'
    payload = {'native':'unchanged'}
    with diagnostics.installed_hooks(native, path):
        model = native.Proteina()
        assert model.predict_step(payload, batch_idx=7) is expected
        assert native.load_ckpt_n_configure_inference(payload) is expected
        assert native.torch.cuda._lazy_init() is expected
    assert native.Proteina.predict_step is original
    assert calls[0] == ((model,payload),{'batch_idx':7})
    rows = [json.loads(x) for x in path.read_text().splitlines()]
    starts = {x['span_id']:x for x in rows if x['event']=='started'}
    ends = [x for x in rows if x['event']=='finished']
    assert len(starts)==len(ends)==3
    assert all(x['span_id'] in starts and x['elapsed_seconds']>=0 for x in ends)
    assert payload == {'native':'unchanged'}


def test_native_exception_is_identical_and_restores_hooks(fake_native, tmp_path):
    native, _, _ = fake_native
    error = ValueError('native scientific failure')
    def fail(*a, **k):
        raise error
    native.get_dataset_from_cfg = fail
    with pytest.raises(ValueError) as caught:
        with diagnostics.installed_hooks(native, tmp_path/'timings.jsonl'):
            native.get_dataset_from_cfg({})
    assert caught.value is error
    assert native.get_dataset_from_cfg is fail
    assert 'native scientific failure' not in (tmp_path/'timings.jsonl').read_text()


def test_source_changes_refuse_unreviewed_diagnostics_before_generation(tmp_path, monkeypatch):
    original=b'known native source'
    monkeypatch.setattr(diagnostics,'SOURCE_HASHES',{'generate.py':hashlib.sha256(original).hexdigest()})
    (tmp_path/'generate.py').write_bytes(original)
    diagnostics.verify_sources(tmp_path)
    (tmp_path/'generate.py').write_bytes(b'new upstream code')
    with pytest.raises(RuntimeError,match='requires|require reviewed'):
        diagnostics.verify_sources(tmp_path)


def test_opt_in_is_explicit_and_diagnostic_write_failure_preserves_work(tmp_path, monkeypatch):
    monkeypatch.delenv('MN_STUDIO_RFD4_NATIVE_TIMING_FILE',raising=False)
    with pytest.raises(RuntimeError,match='explicit timing output'):
        diagnostics.main()
    with diagnostics.phase(tmp_path/'missing'/'timing.jsonl','diagnostic_write_failure'):
        pass
    with pytest.raises(ValueError):
        with diagnostics.phase(tmp_path/'missing'/'timing.jsonl','diagnostic_write_failure'):
            raise ValueError('original')


@pytest.mark.parametrize('value',['yes',1,None])
def test_native_diagnostics_option_rejects_non_boolean(value):
    with pytest.raises(ValueError,match='native_phase_timing'):
        rfd4_proteina_input.validate_options({'options':{'native_phase_timing':value}})


@pytest.mark.parametrize('kwargs',[{'distributed_campaign':{'enabled':True}},{'task':'score_structures'}])
def test_diagnostics_refuse_unsupported_campaign_or_scoring_before_job_creation(monkeypatch,kwargs):
    from mn_molecular_design_studio.engine.workflow_engine.workflows import rfd4_proteina
    monkeypatch.setattr(rfd4_proteina,'create_job',lambda *a,**k:pytest.fail('Unexpected job creation'))
    with pytest.raises(ValueError,match='ordinary local generation only'):
        rfd4_proteina.run_rfd4_proteina_design(options={'native_phase_timing':True},**kwargs)


def test_opt_in_command_preserves_hydra_argv_and_cannot_redirect_scoring(tmp_path):
    command=['rfd4-generate','--config-path','/work/path with spaces','--config-name','effective_config']
    wrapped=diagnostics.instrumented_command(command,tmp_path/'launcher.py')
    assert wrapped[2:]==command[1:]
    assert command[0]=='rfd4-generate'
    with pytest.raises(ValueError,match='only rfd4-generate'):
        diagnostics.instrumented_command(['rfd4-evaluate','inputs=/work/input'],tmp_path/'launcher.py')
