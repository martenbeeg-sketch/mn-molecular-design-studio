"""Generator request validation is CPU-only; never launches Docker or GPU work."""
import pytest

from mn_molecular_design_studio import generator_input as inputs


@pytest.mark.parametrize('engine', list(inputs.ENGINES))
def test_all_bundled_generator_examples_validate_and_preserve_inputs(engine):
    examples=inputs.examples(engine)
    assert examples
    for example in examples:
        p=example['parameters']
        original=open(p['target_pdb'],'rb').read() if p.get('target_pdb') else None
        try:
            validated=inputs.prepare_request(engine,p)
        except Exception as exc:
            pytest.fail(f"{engine} / {example['label']}: {exc}")
        assert validated['generator_only'] is True
        if original is not None:assert open(p['target_pdb'],'rb').read()==original


@pytest.mark.parametrize('engine', list(inputs.ENGINES))
@pytest.mark.parametrize('change',[{'num_designs':0},{'seed':-1},{'num_designs':1.5},{'binder_length':'100-50'},{'cpu_cores':float('nan')}])
def test_invalid_generator_requests_fail_before_job_creation(engine,change):
    example=inputs.examples(engine)[0]['parameters']
    with pytest.raises(ValueError):inputs.prepare_request(engine,example|change)
