"""CPU contract checks: command construction, forms and handoffs; no inference."""
import copy
import inspect
import json
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

from mn_molecular_design_studio import generator_input as inputs, generator_ui as ui, pipeline_graph as graph
from mn_molecular_design_studio.engine.workflow_engine.workflows import generators, graph_pipeline as runner
from mn_molecular_design_studio.engine.workflow_engine.services import workflow_automation as automation


def build_run(tmp_path, engine, example):
    parameters=inputs.prepare_request(engine,example['parameters'])
    parameters['gpu_device']='0'
    root=tmp_path/engine/example['id']
    raw=root/'artifacts/raw'/engine
    (raw/'input').mkdir(parents=True)
    target=generators._stage_structure(parameters['target_pdb'],raw/'input/target.pdb') if parameters.get('target_pdb') else None
    if target:parameters['staged_target']=str(target)
    return {'engine':engine,'params':parameters,'path':root,'raw':raw},target


@pytest.mark.parametrize('engine',list(inputs.ENGINES))
def test_all_bundled_native_commands_are_generation_only(engine,tmp_path,monkeypatch):
    refs=tmp_path/'refs';scaffolds=refs/'rfdiffusion_models/ppi_scaffolds';scaffolds.mkdir(parents=True)
    (scaffolds/'fixture_ss.pt').write_bytes(b'not loaded by command construction')
    monkeypatch.setattr(generators,'reference_root',lambda:refs)
    for example in inputs.examples(engine):
        run,target=build_run(tmp_path,engine,example)
        original=target.read_bytes() if target else None
        steps=getattr(generators,'_'+engine)(run,target)
        assert steps
        for step in steps:
            command=step['command']
            assert isinstance(command,list) and command[:3]==['docker','run','--rm']
            assert inputs.ENGINES[engine]['image'] in command
            assert '--entrypoint' in command
            assert not any('esmfold' in str(v).lower() or 'colabfold' in str(v).lower() for v in command)
        if target:assert target.read_bytes()==original
        if engine=='pxdesign':
            command=steps[-1]['command']
            assert 'infer' in command and 'pipeline' not in command
            assert command[command.index('--N_sample')+1]==str(run['params']['num_designs'])
        elif engine=='genie3':
            command=steps[-1]['command'];config=yaml.safe_load((run['raw']/'experiment.yaml').read_text())
            assert 'generate' in command and 'evaluate' not in command
            assert 'evaluation' not in config and 'inference' not in config['generation']
            assert config['generation']['dataset']['n_sample']==run['params']['num_designs']
        elif engine=='protpardelle':
            command=steps[-1]['command'];config=yaml.safe_load((run['raw']/'sampling.yaml').read_text())
            assert command[command.index('--num-mpnn-seqs')+1]=='0'
            assert 'protpardelle.sample' in command
            if run['params']['task']=='multichain':
                assert len(config['motif_contigs'][0].split(';/;'))==len(config['total_lengths'][0])
        else:
            overrides=json.loads((run['raw']/'native_overrides.json').read_text())
            assert any(v.startswith('inference.ckpt_override_path=/models/') for v in overrides)
            if run['params']['symmetry']=='T' and run['params']['task']=='symmetry':assert 'inference.symmetry=tetrahedral' in overrides


@pytest.mark.parametrize('positions',['79','0-80','5-2','-1','1;2'])
def test_partial_diffusion_rejects_invalid_fixed_sequence_positions(positions):
    example=next(e for e in inputs.examples('rfdiffusion') if e['id']=='partial')
    with pytest.raises(ValueError,match='Fixed sequence positions'):
        inputs.prepare_request('rfdiffusion',example['parameters']|{'fixed_sequence_positions':positions})


def test_rf_seed_and_provide_sequence_are_reproducible(tmp_path):
    example=copy.deepcopy(next(e for e in inputs.examples('rfdiffusion') if e['id']=='partial'))
    example['parameters'].update(seed=42,deterministic=True,fixed_sequence_positions='0-5')
    run,target=build_run(tmp_path,'rfdiffusion',example)
    steps=generators._rfdiffusion(run,target)
    overrides=json.loads((run['raw']/'native_overrides.json').read_text())
    assert 'inference.design_startnum=42' in overrides
    assert 'contigmap.provide_seq=[0-5]' in overrides
    assert run['params']['resolved_model']=='InpaintSeq'
    assert all(v in steps[-1]['command'] for v in overrides)


def test_native_multichain_without_fixed_target_is_supported(tmp_path):
    example={'id':'unconditioned','parameters':inputs.default_parameters('protpardelle','multichain')|{
        'motif_contig':'50;/;60','chain_lengths':'50;60','model_name':'cc95'}}
    run,target=build_run(tmp_path,'protpardelle',example)
    generators._protpardelle(run,target)
    config=yaml.safe_load((run['raw']/'sampling.yaml').read_text())
    assert config['motifs']==[None]
    assert config['motif_contigs']==['50-50;/;60-60']


def test_redesign_exports_native_fixed_identity_mapping(tmp_path):
    example=next(e for e in inputs.examples('rfdiffusion') if e['id']=='redesign')
    run,target=build_run(tmp_path,'rfdiffusion',example)
    selected={v for b in inputs.residue_ranges('A20-40',inputs.inventory(target)) for v in b}
    expected=set(generators._residues(target))-selected
    assert set(run['params']['fixed_identity_labels'])==expected
    output=run['raw']/'output';output.mkdir()
    lines=[]
    for line in target.read_text().splitlines():
        if line.startswith('ATOM'):
            label=line[21:22]+line[22:26].strip()
            line=line[:60]+f'{float(label in expected):6.2f}'+line[66:]
        lines.append(line)
    structure=output/'design_0.pdb';structure.write_text('\n'.join(lines)+'\n')
    metadata=generators._fixed_residues('rfdiffusion',structure,run['params'])
    assert metadata['motif_mapping_verified'] and set(metadata['fixed_residues'])==expected
    structure.write_text(structure.read_text().replace('  1.00','  0.00'))
    assert not generators._fixed_residues('rfdiffusion',structure,run['params'])['motif_mapping_verified']


@pytest.mark.parametrize('engine',list(inputs.ENGINES))
def test_runner_creates_one_correct_job_and_never_executes_downstream_models(engine,monkeypatch,tmp_path):
    example=inputs.examples(engine)[0]['parameters']
    path=tmp_path/engine;path.mkdir()
    jobs=[];steps=[];finished=[]
    monkeypatch.setattr(inputs,'readiness',lambda *a,**k:{'ready':True,'missing':[]})
    monkeypatch.setattr(generators,'create_job',lambda *a:jobs.append(a) or SimpleNamespace(run_dir=path))
    monkeypatch.setattr(generators,'update_status',lambda *a,**k:None)
    monkeypatch.setattr(generators.design,'_run_shell_steps',lambda p,s:steps.extend(s) or 1)
    monkeypatch.setattr(generators.design,'_finish_design_job',lambda *a,**k:finished.append(a))
    function=getattr(generators,'run_'+engine+'_generation')
    public=ui.public_parameters(engine,example)
    inspect.signature(function).bind(**public)
    assert function(**public)==path
    assert len(jobs)==1 and jobs[0][1:3]==('design_generation',inputs.TOOLS[engine])
    assert steps and len(finished)==1 and finished[0][1] is False
    assert json.loads((path/'input.json').read_text())['params']['generator_only'] is True


@pytest.mark.parametrize('engine',list(inputs.ENGINES))
def test_shared_page_renders_every_task_and_loads_example_without_submitting(engine,monkeypatch,tmp_path):
    from streamlit.testing.v1 import AppTest
    from mn_molecular_design_studio import target_editor
    monkeypatch.setattr(ui,'user_input_root',lambda _:tmp_path/'inputs')
    monkeypatch.setattr(target_editor,'user_input_root',lambda _:tmp_path/'inputs')
    app=AppTest.from_string(f'''
from mn_molecular_design_studio.generator_ui import render_setup
from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import workflow_schema
render_setup(workflow_schema({inputs.WORKFLOWS[engine]!r}),owner_id='test',project_selector=lambda _:None,
    submit=lambda *a:(_ for _ in ()).throw(AssertionError('Submission is forbidden in this test')))
''').run(timeout=20)
    assert not app.exception
    prefix='generator-'+engine
    for task in inputs.ENGINES[engine]['tasks']:
        app.button_group(key=prefix+'-task').set_value(task).run(timeout=20)
        assert not app.exception
        choices=inputs.examples(engine,task)
        if choices:
            app.button(key=prefix+'-fields-'+task+'-load-example').click().run(timeout=20)
            assert not app.exception
            assert not any('absent residues' in c.value for c in app.caption)
        assert app.button(key=prefix+'-submit').disabled


@pytest.mark.parametrize('engine',list(inputs.ENGINES))
def test_shared_canvas_node_form_and_template_bind_to_public_generator(engine):
    from streamlit.testing.v1 import AppTest
    app=AppTest.from_string(f'''
from mn_molecular_design_studio.pipeline_graph import new_graph
from mn_molecular_design_studio.pipeline_node_ui import render_node_form
g=new_graph({engine!r})
node=next(n for n in g['nodes'] if n['data']['kind']=={engine!r})
render_node_form(node,g,key='test-generator',owner_id='test',render_input=lambda:{{}})
''').run(timeout=20)
    assert not app.exception
    assert any(b.label=='Design task' for b in app.button_group)
    assert not any(n.label=='Number of independent jobs' for n in app.number_input)
    g=graph.new_graph(engine)
    generator_node=next(n for n in g['nodes'] if n['data']['kind']==engine)
    if engine == 'protpardelle':
        pass  # This starter template is intentionally an unconditional design.
    else:
        example=inputs.examples(engine)[0]['parameters']
        target_node=next(n for n in g['nodes'] if n['data']['kind']=='target')
        target_node['data']['params'].update(target_pdb=example['target_pdb'],
            target_chains=','.join(example.get('target_chains') or ['A']),
            hotspots=example.get('hotspots',''))
        generator_node['data']['params'].update({k:('' if v is None else v) for k,v in example.items()
            if k in generator_node['data']['params'] and k not in {'target_pdb','target_chains','hotspots'}})
    plan=graph.validate_graph(g,require_inputs=True)
    assert plan['graph']['template_id']==engine
    request=runner.node_request(engine,plan['graph']['nodes'][0]['data']['params'],{},'test')
    automation._validate_workflow_inputs(automation.WORKFLOWS[request['workflow']],request['parameters'])
    assert graph.new_graph(engine+'_generated')['nodes'][0]['data']['kind']=='backbone_source'


def test_hidden_legacy_pipelines_remain_available():
    for kind in ('pxdesign_pipeline','genie3_pipeline','protpardelle_pipeline'):
        assert graph.NODE_CATALOG[kind]['hidden']
        assert graph.NODE_CATALOG[kind]['workflow'] in automation.WORKFLOWS


@pytest.mark.parametrize('engine',list(inputs.ENGINES))
def test_loaded_task_examples_remain_valid_after_shared_form_editors(engine,monkeypatch,tmp_path):
    from streamlit.testing.v1 import AppTest
    from mn_molecular_design_studio import target_editor
    monkeypatch.setattr(ui,'user_input_root',lambda _:tmp_path/'inputs')
    monkeypatch.setattr(target_editor,'user_input_root',lambda _:tmp_path/'inputs')
    app=AppTest.from_string(f'''
import streamlit as st
from mn_molecular_design_studio.generator_ui import render_inputs
with st.container(key='generator-surface'):
    params,jobs,error=render_inputs({engine!r},owner_id='test',key_prefix='example-test')
    st.session_state['returned_parameters']=params
    st.session_state['returned_error']=error
''').run(timeout=20)
    for task in inputs.ENGINES[engine]['tasks']:
        app.button_group(key='example-test-task').set_value(task).run(timeout=20)
        if not inputs.examples(engine,task):continue
        app.button(key='example-test-fields-'+task+'-load-example').click().run(timeout=20)
        assert not app.exception
        assert app.session_state['returned_error'] is None,(engine,task,app.session_state['returned_error'])
        parameters=inputs.prepare_request(engine,app.session_state['returned_parameters'])
        example=inputs.prepare_request(engine,inputs.examples(engine,task)[0]['parameters'])
        for field in ('hotspots','motif_labels','contig','binder_length'):
            if field in example:assert parameters[field]==example[field],(engine,task,field)


def test_residue_picker_accepts_bounded_ranges_and_long_chain_ids():
    from mn_molecular_design_studio.target_editor import _hotspot_pairs
    assert _hotspot_pairs('A30-32,A40,AA10-11')=={('A',30),('A',31),('A',32),('A',40),('AA',10),('AA',11)}
    assert not _hotspot_pairs('A40-30,A1-10000000')


def test_checkpoint_names_are_preserved_in_portable_requests(tmp_path):
    from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_managed_paths
    for name in inputs.RF_MODELS:
        assert resolve_managed_paths({'model_weights':name},run_dir=tmp_path)=={'model_weights':name}


@pytest.mark.parametrize('engine',list(inputs.ENGINES))
def test_generator_api_preflight_enforces_input_access_without_enqueue(engine,monkeypatch,tmp_path):
    import shutil
    from fastapi.testclient import TestClient
    from mn_molecular_design_studio import api
    monkeypatch.setenv('MDS_CLUSTER_HOME',str(tmp_path/'cluster'))
    uploaded=tmp_path/'inputs';uploaded.mkdir()
    example=inputs.examples(engine,'binder')[0]['parameters']
    target=uploaded/Path(example['target_pdb']).name;shutil.copy2(example['target_pdb'],target)
    monkeypatch.setattr(api,'user_input_root',lambda _:uploaded)
    monkeypatch.setattr(api,'reference_root',lambda:tmp_path/'refs')
    monkeypatch.setattr(api,'path_is_within_accessible_job',lambda *a,**k:False)
    monkeypatch.setattr(api,'accessible_project_ids',lambda _:[])
    monkeypatch.setattr(api,'_request_payload_and_project',lambda payload,user:({k:v for k,v in payload.items() if k!='project_id'},payload['project_id']))
    monkeypatch.setattr(runner,'enqueue_graph_pipeline',lambda *a,**k:pytest.fail('No job may be queued'))
    g=graph.new_graph(engine)
    g['nodes'][0]['data']['params'].update({k:('' if v is None else v) for k,v in example.items() if k in g['nodes'][0]['data']['params']})
    g['nodes'][0]['data']['params'].update(target_pdb=str(target),target_chains=','.join(example['target_chains']))
    api.app.dependency_overrides[api._current_user]=lambda:{'id':'test-user','role':'user'}
    try:
        client=TestClient(api.app)
        response=client.post('/api/v1/pipelines/validate',json={'project_id':'project','graph':g})
        assert response.status_code==200,response.text
        g['nodes'][0]['data']['params']['target_pdb']='/etc/passwd'
        assert client.post('/api/v1/pipelines/validate',json={'project_id':'project','graph':g}).status_code==403
    finally:api.app.dependency_overrides.clear()


@pytest.mark.parametrize('engine',list(inputs.ENGINES))
def test_generated_source_is_separate_from_filter_passing_designs(engine,design_source_factory):
    source,candidate=design_source_factory(passed=False)
    metadata=json.loads((source/'metadata.json').read_text());metadata['tool']=inputs.TOOLS[engine]
    (source/'metadata.json').write_text(json.dumps(metadata))
    (source/'input.json').write_text(json.dumps({'params':{'generator_only':True}}))
    candidate.update(stage='generation.backbone',source_tool=inputs.TOOLS[engine],metrics={'result_kind':'generated','pass_filters':None})
    (source/'artifacts/normalized_candidates/candidates.jsonl').write_text(json.dumps(candidate)+'\n')
    assert len(runner.source_generations(str(source),'test-project',stage='generation.backbone')[1])==1
    with pytest.raises(ValueError,match='project'):runner.source_generations(str(source),'other-project')
    with pytest.raises(ValueError,match='filter-passing'):runner.source_designs(str(source),'test-project')
    with pytest.raises(ValueError):runner.source_generations(str(source),'test-project',stage='generation.backbone_sequence')


def test_partial_complex_normalization_preserves_chain_roles_and_generated_status(tmp_path,design_source_factory):
    source,_=design_source_factory(missing_occupancy=False)
    example={'id':'complex-partial','parameters':inputs.default_parameters('rfdiffusion','partial')|{
        'target_pdb':str(source/'designed.cif'),'designed_chain':'B','designed_residues':'B1-5'}}
    run,target=build_run(tmp_path,'rfdiffusion',example)
    output=run['raw']/'output';output.mkdir()
    lines=[]
    for line in target.read_text().splitlines():
        if line.startswith('ATOM'):line=line[:60]+f'{float(line[21:22]=="A"):6.2f}'+line[66:]
        lines.append(line)
    (output/'design_0.pdb').write_text('\n'.join(lines)+'\n')
    candidates=generators.normalize_outputs(run,target)
    assert len(candidates)==1
    candidate=candidates[0]
    assert candidate['target_chains']==['A'] and candidate['binder_chains']==['B']
    assert candidate['stage']=='generation.backbone' and candidate['binder_sequence'] is None
    assert candidate['metrics']['pass_filters'] is None
    assert candidate['raw_metadata']['motif_mapping_verified']


def test_pxdesign_normalization_ignores_target_cif_and_keeps_xpb_binder_chain(tmp_path):
    import gemmi

    example=next(row for row in inputs.examples('pxdesign') if row['id']=='tamarind')
    params=inputs.prepare_request('pxdesign',example['parameters'])
    run_dir=tmp_path/'pxdesign-run';raw=run_dir/'artifacts/raw/pxdesign'
    predictions=raw/'output/target/seed_1/predictions'
    (raw/'input').mkdir(parents=True);predictions.mkdir(parents=True)
    target=generators._stage_structure(params['target_pdb'],raw/'input/target.pdb')
    params['staged_target']=str(target)
    # PXDesign's converted input is also named target.cif, but is not a
    # generated model and Gemmi cannot parse it as one.
    (raw/'output/target.cif').write_text('data_target\n#\n')

    structure=gemmi.read_structure(params['target_pdb'])
    structure[0][0].name='A0'
    binder=gemmi.Chain('B0')
    for index in range(1,11):
        residue=gemmi.Residue();residue.name='XPB';residue.seqid=gemmi.SeqId(index,' ')
        for atom_name,element,x in [('N','N',0.),('CA','C',1.),('C','C',2.),('O','O',3.)]:
            atom=gemmi.Atom();atom.name=atom_name;atom.element=gemmi.Element(element)
            atom.pos=gemmi.Position(x,index*4.,0.);residue.add_atom(atom)
        binder.add_residue(residue)
    structure[0].add_chain(binder)
    structure.make_mmcif_document().write_file(str(predictions/'target_sample_0.cif'))

    candidates=generators.normalize_outputs(
        {'engine':'pxdesign','params':params,'path':run_dir,'raw':raw},target
    )
    assert len(candidates)==1
    assert candidates[0]['target_chains']==['A0']
    assert candidates[0]['binder_chains']==['B0']
    assert candidates[0]['raw_metadata']['target_chain_inference']=='sequence_match'
    assert candidates[0]['stage']=='generation.backbone'
