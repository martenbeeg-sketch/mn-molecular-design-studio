"""RFD4 request, tile, and pipeline contracts; model checks are recorded separately."""
import json
from pathlib import Path
import pytest
from mn_molecular_design_studio import rfd4_proteina_input as inputs, pipeline_graph
from mn_molecular_design_studio.engine.workflow_engine.workflows import rfd4_proteina as workflow
from mn_molecular_design_studio.engine.workflow_engine.workflows.graph_pipeline import node_request
from mn_molecular_design_studio.engine.workflow_engine.services import workflow_automation as automation


def test_all_thirteen_tasks_and_native_operations():
    assert len(inputs.TASKS)==13
    assert 'score_structures' in inputs.TASKS
    assert set(inputs.OPERATIONS)=={'diversify','inverse_folding','structure_prediction','sidechain_packing'}


@pytest.mark.parametrize('seed',[None,0,4294967295])
def test_unconditional_range_seed_and_effort(tmp_path,seed):
    p=inputs.prepare_request({'task':'unconditional','length_min':25,'length_max':37,
                             'seed':seed,'nsteps':60,'options':{'sampling_effort':'Fast'}},run_dir=tmp_path)
    request=json.loads(p['request_path'].read_text())
    assert request['parameters']['seed']==seed
    assert request['parameters']['nsteps']==60
    assert request['parameters']['length_max']==37


@pytest.mark.parametrize('changes',[{'nsteps':30},{'length_min':0},{'seed':-1},{'num_designs':0},
                                    {'options':{'operation':'invalid'}},{'options':{'folding_engine':'invented'}}])
def test_invalid_requests_fail_before_running(tmp_path,changes):
    with pytest.raises(ValueError):
        inputs.prepare_request({'task':'unconditional',**changes},run_dir=tmp_path)


def test_sequence_carrier_and_folding_registration(tmp_path):
    from mn_molecular_design_studio.design_validation import read_structure,protein_chains,sequence
    path=workflow.sequence_structure('>a\nACDEFGHIK\n>b\nLMNPQRSTV\n',tmp_path/'sequence.cif')
    assert [sequence(chain) for chain in protein_chains(read_structure(path)).values()]==['ACDEFGHIK','LMNPQRSTV']
    schema=automation.workflow_schema('refolding.rfd4_proteina')
    assert schema['workflow']=='refolding.rfd4_proteina'
    assert 'protein_sequences' in json.dumps(schema)


def defaults(kind):
    return {f['name']:f['default'] for f in pipeline_graph.NODE_CATALOG[kind]['fields']}


def test_pipeline_generation_and_native_folding_requests():
    source={'target':{'target_pdb':'/tmp/target.cif','target_chains':['A'],'hotspots':'A25'}}
    request=node_request('rfd4',defaults('rfd4'),source,'design')
    assert request['workflow']=='design.rfd4_proteina'
    assert request['parameters']['options']['refold'] is False
    assert request['parameters']['target_structure']=='/tmp/target.cif'
    request=node_request('rfd4_fold',defaults('rfd4_fold'),{'molecule_task':{'target_protein_sequences':'ACDE','partner_protein_sequences':'FGHI'}},'fold')
    assert request['parameters']['protein_sequences']=='ACDE:FGHI'
    assert request['workflow']=='refolding.rfd4_proteina'
    assert pipeline_graph.NODE_CATALOG['rfd4_fold']['inputs'][0]['optional']


def test_folding_node_rejects_molecules_it_cannot_fold():
    with pytest.raises(ValueError,match='protein sequences'):
        node_request('rfd4_fold',defaults('rfd4_fold'),{'molecule_task':{'dna_sequences':'ACGT'}},'fold')


def test_task_compatible_engine_choices():
    from mn_molecular_design_studio.engine.workflow_engine.workflows.rfd4_oracle import compatible_engines
    assert 'boltz2' in compatible_engines('protein_binder')
    assert 'colabfold' in compatible_engines('protein_binder')
    assert 'colabfold' not in compatible_engines('dna_binder')


@pytest.mark.parametrize('options',[
    {'cdr_lengths':[{'lower_bound':8,'upper_bound':4}]},
    {'cdr_lengths':[{'lower_bound':1,'upper_bound':41}]},
    {'scaffold_lengths':[{'lower_bound':0,'upper_bound':5}]},
    {'ligand_interactions':[{'interaction':'invented','ligand_atoms':'O1'}]},
    {'backbone_level':1.}, {'sequence_level':0.}, {'redesign_shell':31.},
    {'framework_atoms':'invented'}, {'hold_kept_in_place':'yes'},
])
def test_task_specific_invalid_controls(options):
    with pytest.raises(ValueError):inputs.validate_options({'options':options})


def test_task_specific_controls_and_numeric_chain_mapping():
    options={'cdr_lengths':[{'lower_bound':5,'upper_bound':12}],
             'framework_atoms':'backbone','dock_independently':True,'keep_constant_domains':False,
             'ligand_interactions':[{'interaction':name,'ligand_atoms':'SAM:O2'} for name in inputs.LIGAND_INTERACTIONS],
             'noise_mode':'levels','backbone_level':.8,'sequence_level':.7,
             'refold_backbone_rmsd_limit':2.,'motif_rmsd_limit':1.}
    assert inputs.validate_options({'options':options})==options
    assert inputs._translate_residues('1:102,A2:5-6',{'1':'A','A2':'B'},
                                     {('1',102):2,('A2',5):3,('A2',6):4})=='A:2,B:3,B:4'


def test_native_fold_preserves_parent_chain_roles(tmp_path):
    from mn_molecular_design_studio.design_validation import read_structure,protein_chains
    root=tmp_path/'artifacts/rfd4_proteina';root.mkdir(parents=True)
    structure=workflow.sequence_structure('ACDEFG:LMNPQR',root/'folded.cif')
    (root/'metrics.csv').write_text('example_id,score\nfolded,1\n')
    (root/'result_index.json').write_text(json.dumps({'structures':[{
        'structure':str(structure),'designed_chains':['A','B'],'target_chains':[],
        'sequences':{'A':'ACDEFG','B':'LMNPQR'},'residue_count':12,
        'source_candidate':{'candidate_id':'parent','binder_chains':['B'],'target_chains':['A']}
    }],'metrics_files':[str(root/'metrics.csv')]}))
    row=workflow.normalize_outputs(tmp_path,{'root':root,'parameters':{'options':{'operation':'structure_prediction','refold':False}}},'existing_structure',1)[0]
    assert row['parents']==['parent']
    assert row['binder_chains']==['B'] and row['target_chains']==['A']
    assert row['binder_sequence']=='LMNPQR'
    assert list(protein_chains(read_structure(Path(row['binder_pdb']))))==['B']


def test_clicking_selected_pill_keeps_its_choice(monkeypatch):
    from mn_molecular_design_studio import rfd4_proteina_ui as ui
    state={'task':'nanobody'};monkeypatch.setattr(ui.st,'session_state',state)
    ui._retain_pill_choice('task','protein_binder')
    state['task']=None
    ui._retain_pill_choice('task','protein_binder')
    assert state['task']=='nanobody'


def test_pd_l1_example_replaces_stale_form_state_and_loads_complete_inputs(monkeypatch,tmp_path):
    from mn_molecular_design_studio import rfd4_proteina_ui as ui

    state={
        'rfd4-proteina-project-id':'project-42',
        'rfd4-task-pills':'existing_structure',
        'rfd4-target-chains-protein_binder':[],
        'rfd4-hotspots-protein_binder':'A123',
        'rfd4-target-target-upload':'stale upload',
        'rfd4-refold-protein_binder':False,
        'other-app-selection':'keep me',
    }
    monkeypatch.setattr(ui.st,'session_state',state)
    monkeypatch.setattr(ui,'user_input_root',lambda owner:tmp_path/owner/'inputs')

    ui._load_example_inputs('pd_l1_binder','test-user')

    assert state['rfd4-proteina-project-id']=='project-42'
    assert state['other-app-selection']=='keep me'
    assert state['rfd4-task-pills']=='protein_binder'
    assert state['rfd4-target-chains-protein_binder']==['A']
    assert state['rfd4-hotspots-protein_binder']=='A40,A99,A107'
    assert Path(state['rfd4-target-fetched-target']).is_file()
    assert state['rfd4-target-format']=='PDB'
    assert state['rfd4-min-protein_binder']==60
    assert state['rfd4-max-protein_binder']==120
    assert state['rfd4-count']==8
    assert 'rfd4-target-target-upload' not in state
    assert state['rfd4-prioritize'] is False


def test_pd_l1_example_button_renders_preselected_chain_and_hotspots(tmp_path):
    from streamlit.testing.v1 import AppTest

    app=AppTest.from_string(f'''
from pathlib import Path
from mn_molecular_design_studio import rfd4_proteina_ui as ui
ui.user_input_root=lambda owner:Path({str(tmp_path)!r})/owner/'inputs'
ui.render_setup({{}},owner_id='test-user',project_selector=lambda _:None,submit=lambda *args:None)
''').run()
    assert not app.exception

    app.button(key='rfd4-load-example').click().run()

    assert not app.exception
    assert app.multiselect(key='rfd4-target-chains-protein_binder').value==['A']
    assert app.text_input(key='rfd4-hotspots-protein_binder').value=='A40,A99,A107'


def test_prepare_only_validation_never_requests_gpu_or_generation(monkeypatch):
    import subprocess
    from mn_molecular_design_studio.engine.workflow_engine.workflows import rfd4_proteina_prepare as prepare

    captured={}

    def fake_prepare_request(parameters, *, run_dir):
        root=run_dir/'artifacts/rfd4_proteina'
        root.mkdir(parents=True)
        request_path=root/'request.json'
        request_path.write_text(json.dumps({'root':'/work/artifacts/rfd4_proteina'}))
        return {'root':root,'request_path':request_path}

    def fake_build(manifest, run_dir, args, **kwargs):
        captured.update(manifest=manifest,args=args,kwargs=kwargs,run_dir=run_dir)
        return ['docker','run','--rm','--entrypoint','python',manifest['image'],*args]

    def fake_run(command, **kwargs):
        captured['command']=command
        request=json.loads((captured['run_dir']/'artifacts/rfd4_proteina/request.json').read_text())
        output=captured['run_dir']/request['root'].removeprefix('/work').lstrip('/')/'prepared_inputs.json'
        output.write_text(json.dumps([{'name':'test','atoms':12,'designed_chains':['B'],
                                       'target_chains':['A'],'warnings':[],'conditions':'valid',
                                       'validation':'passed'}]))
        return subprocess.CompletedProcess(command,0,'','')

    monkeypatch.setattr(prepare.inputs,'prepare_request',fake_prepare_request)
    monkeypatch.setattr(prepare,'build_docker_command',fake_build)
    monkeypatch.setattr(prepare.subprocess,'run',fake_run)
    report=prepare.validate_inputs({'task':'protein_binder'})

    assert report['count']==1
    assert captured['manifest']['gpu'] is False
    assert captured['args'][-1]=='--prepare-only'
    assert 'rfd4-generate' not in captured['command']
    assert 'rfd4-evaluate' not in captured['command']


def test_prepare_only_sequence_validation_builds_sequence_carrier(monkeypatch):
    from mn_molecular_design_studio.engine.workflow_engine.workflows import rfd4_proteina_prepare as prepare

    captured={}

    def fake_prepare_request(parameters, *, run_dir):
        captured['parameters']=parameters
        root=run_dir/'artifacts/rfd4_proteina'
        root.mkdir(parents=True)
        request_path=root/'request.json'
        request_path.write_text('{}')
        return {'root':root,'request_path':request_path}

    def fake_build(manifest, run_dir, args, **kwargs):
        captured['manifest']=manifest
        output=run_dir/'artifacts/rfd4_proteina/prepared_inputs.json'
        output.write_text(json.dumps([{'validation':'passed','name':'fold','atoms':8,'warnings':[]}]))
        return ['docker','run',manifest['image'],*args]

    monkeypatch.setattr(prepare.inputs,'prepare_request',fake_prepare_request)
    monkeypatch.setattr(prepare,'build_docker_command',fake_build)
    monkeypatch.setattr(prepare.subprocess,'run',lambda *a,**k:type('Result',(),{'returncode':0,'stderr':'','stdout':''})())
    report=prepare.validate_inputs({'task':'existing_structure','protein_sequences':'ACDE',
                                    'options':{'operation':'structure_prediction'}})

    assert report['count']==1
    assert Path(captured['parameters']['target_structure']).suffix=='.cif'
    assert 'protein_sequences' not in captured['parameters']
    assert captured['manifest']['gpu'] is False


def test_queued_rfd4_folding_adopts_refolding_job_and_completes(monkeypatch, tmp_path):
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json, write_json
    from mn_molecular_design_studio.engine.workflow_engine.core.workflow_queue import run_workflow_call

    run_dir=tmp_path/'refolding'/'queued-fold'
    (run_dir/'artifacts').mkdir(parents=True)
    write_json(run_dir/'metadata.json',{'run_id':'queued-fold','task_group':'refolding',
                                        'tool':'rfd4_proteina','job_type':'rfd4_structure_prediction',
                                        'status':'running'})

    captured={}
    monkeypatch.setattr(workflow,'reference_root',lambda:tmp_path/'references')
    monkeypatch.setattr(workflow,'run_docker_job',lambda *args,**kwargs:captured.update(docker=(args,kwargs)) or 0)
    monkeypatch.setattr(workflow,'normalize_outputs',lambda run,prepared,task,requested:
                        captured.update(normalized=(prepared,task,requested)) or
                        [{'candidate_id':'folded','complex_pdb':str(tmp_path/'folded.cif'),
                          'binder_sequence':'ACDEFGHIK'}])
    monkeypatch.setattr(workflow,'write_candidates',lambda *args,**kwargs:['candidate-row'])
    from mn_molecular_design_studio import molecular_database_ingest
    monkeypatch.setattr(molecular_database_ingest,'record_completed_run',lambda *_args,**_kwargs:None)

    result=run_workflow_call(run_dir,{'module':workflow.__name__,'function':'run_rfd4_structure_prediction',
                                     'kwargs':{'protein_sequences':'ACDEFGHIK','num_designs':1}})

    assert result==run_dir
    assert captured['normalized'][1:] == ('existing_structure',1)
    payload=read_json(run_dir/'input.json')
    assert payload['job_type']=='rfd4_structure_prediction'
    assert payload['tool']=='rfd4_proteina'
    assert payload['params']['task_group']=='refolding'
    metadata=read_json(run_dir/'metadata.json')
    assert metadata['task_group']=='refolding'
    assert metadata['job_type']=='rfd4_structure_prediction'
    assert metadata['status']=='completed'
    assert read_json(run_dir/'result.json')['success'] is True

    from mn_molecular_design_studio.engine.workflow_engine.services import workflow_automation as automation
    schema=automation.workflow_schema('refolding.rfd4_proteina')
    assert 'task_group' not in {field['name'] for field in schema['parameters']}


def test_native_acceptance_campaign_requires_rfd4_refold():
    from mn_molecular_design_studio.validation_acceptance import METRICS

    campaign = {
        "enabled": True, "acceptance_level": "validated", "accepted_target": 1,
        "max_attempts": 1,
        "validation_stages": [{
            "workflow": "analysis.validation_acceptance",
            "acceptance_policy": {
                "binder_format": "miniprotein", "target_type": "fasta_idr",
                "mandatory_filters": [{
                    "field": METRICS["rfd4_native_iptm"]["field"],
                    "operator": ">=", "value": 0.7,
                }],
            },
        }],
    }
    with pytest.raises(ValueError, match="requires RFD4 options.refold=true"):
        inputs.validate_campaign_inputs({
            "task": "protein_binder", "input_mode": "guided_binder",
            "target_structure": "/tmp/target.cif", "target_chains": ["A"],
            "options": {"refold": False}, "distributed_campaign": campaign,
        })
