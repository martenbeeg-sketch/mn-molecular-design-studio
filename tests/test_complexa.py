"""Proteina-Complexa task, command, output and UI regressions; no GPU inference."""
import copy
import csv
import json
import shutil
import subprocess
from pathlib import Path
from types import SimpleNamespace

import pytest

from mn_molecular_design_studio import complexa_input as inputs, complexa_ui as ui
from mn_molecular_design_studio.engine.workflow_engine.workflows import design
from mn_molecular_design_studio.engine.workflow_engine.services import workflow_automation as automation


@pytest.mark.parametrize('task',list(inputs.TASKS))
def test_all_task_examples_compile_and_originals_survive(task,tmp_path):
    parameters=ui.example_parameters(task)
    original=Path(parameters['target_pdb']).read_bytes()
    p=inputs.prepare_request(parameters);p.update(run_name='test',seed=17)
    target,source=inputs.stage_target(p,tmp_path)
    assert target.is_file() and source.read_bytes()==original
    assert Path(parameters['target_pdb']).read_bytes()==original
    cfg=inputs.native_overrides(p,'A1-115')
    assert cfg['generation']['args']['nsteps']==400
    assert cfg['metric']['keep_folding_outputs'] is True
    assert cfg['ckpt_name']==inputs.CHECKPOINTS[task][0]
    group='motif_target_dict_cfg' if task=='motif_ligand' else 'target_dict_cfg'
    custom=cfg['generation'][group]['MN_APP_TARGET']
    if task=='motif_ligand':assert custom['contig_atoms'].startswith('B1: [O, C]')
    if task!='protein_binder':assert custom['ligand']=='L:0'


@pytest.mark.parametrize('changes',[{'num_designs':0},{'n_steps':49},{'replicas':9},{'guidance_weight':1.1},{'sampling_noise':float('nan')},{'seed':-1},{'target_chains':['Z']},{'hotspots':'A9999'},{'binder_length':'150-70'},{'search_strategy':'shell-command'},{'deduplicate_sequences':'false'}])
def test_invalid_options_rejected_before_job(changes):
    with pytest.raises(ValueError):inputs.prepare_request(ui.example_parameters('protein_binder')|changes)


def test_ligand_instance_and_motif_atoms_are_validated():
    with pytest.raises(ValueError,match='ligand instance'):inputs.prepare_request(ui.example_parameters('ligand_binder')|{'ligand':'MISSING'})
    for text in ['B64: [MISSING]','B9999: [CA]','B64: [O]; B64: [C]','B64: []','__import__(x)']:
        with pytest.raises(ValueError):inputs.prepare_request(ui.example_parameters('motif_ligand')|{'motif_atom_spec':text})
    with pytest.raises(ValueError,match='Single Pass'):inputs.prepare_request(ui.example_parameters('motif_ligand')|{'search_strategy':'best-of-n'})


def test_api_preflight_validates_all_native_tasks():
    spec=automation.WORKFLOWS['design.proteina_complexa']
    for task in inputs.TASKS:automation._validate_workflow_inputs(spec,ui.example_parameters(task))
    with pytest.raises(automation.AutomationError):automation._validate_workflow_inputs(spec,ui.example_parameters('ligand_binder')|{'ligand':'MISSING'})
    fields={p['name'] for p in automation.workflow_schema(spec.workflow_id)['parameters']}
    assert {'pipeline_task','ligand','motif_atom_spec','search_strategy','guidance_weight','sampling_noise','deduplicate_sequences'}<=fields


def test_independent_jobs_preserve_counts_names_and_seeds():
    p=ui.example_parameters('protein_binder')|{'campaign_name':'complexa','seed':17}
    jobs=ui.job_parameters(p,3)
    assert [j['seed'] for j in jobs]==[17,18,19]
    assert [j['num_designs'] for j in jobs]==[4,4,4]
    assert len({j['campaign_name'] for j in jobs})==3
    assert p['campaign_name']=='complexa'


def test_missing_weights_reported_before_creating_job(monkeypatch,tmp_path):
    monkeypatch.setattr(design,'reference_root',lambda:tmp_path)
    monkeypatch.setattr(design,'create_job',lambda *a:pytest.fail('No job should be created'))
    with pytest.raises(ValueError,match='reference files are missing'):design.run_proteina_complexa(**ui.example_parameters('protein_binder'))


def _reference_fixture(root):
    paths=inputs.resource_paths(root)
    for names in inputs.CHECKPOINTS.values():
        for name in names:
            p=paths['checkpoints']/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(b'weights')
    paths['rf3'].parent.mkdir(parents=True);paths['rf3'].write_bytes(b'weights')
    paths['af2'].mkdir();paths['mpnn'].mkdir(parents=True);(paths['mpnn']/'ligandmpnn_v_32_010_25.pt').write_bytes(b'weights')
    paths['tools'].mkdir(parents=True)
    for name in ['foldseek','mmseqs']:(paths['tools']/name).write_bytes(b'#!/bin/sh\n')
    for model in ['esmfold_v1','esm2_t33_650M_UR50D']:
        folder=paths['hf_cache']/f'models--facebook--{model}'/'snapshots/test';folder.mkdir(parents=True)
        (folder/'model.safetensors').write_bytes(b'weights')


@pytest.mark.parametrize('task',list(inputs.TASKS))
def test_runner_keeps_scheduler_contract_and_stages_correct_native_config(monkeypatch,tmp_path,task):
    refs=tmp_path/'refs';_reference_fixture(refs)
    run=tmp_path/'run';run.mkdir()
    monkeypatch.setattr(design,'reference_root',lambda:refs)
    monkeypatch.setattr(design,'create_job',lambda *a:SimpleNamespace(run_dir=run))
    monkeypatch.setattr(design,'update_status',lambda *a,**k:None)
    captured=[]
    monkeypatch.setattr(design,'_run_shell_steps',lambda run,steps:captured.extend(steps) or 1)
    monkeypatch.setattr(design,'_finish_design_job',lambda *a:None)
    design.run_proteina_complexa(**ui.example_parameters(task))
    command=captured[0]['command']
    assert command[command.index('--entrypoint')+1]=='/workspace/.venv/bin/python'
    assert inputs.IMAGE in command and any('/rf3.ckpt:ro' in arg for arg in command)
    assert any('/model_params:ro' in arg for arg in command)
    assert any('/studio-tools:ro' in arg for arg in command)
    spec=json.loads((run/'artifacts/raw/proteina_complexa/input/target_config.json').read_text())
    assert spec['base_config']==inputs.TASK_CONFIGS[task]
    assert spec['overrides']['ckpt_name']==inputs.CHECKPOINTS[task][0]
    assert spec['overrides']['gen_njobs']==spec['overrides']['eval_njobs']==1
    assert spec['environment']['HF_HUB_CACHE']=='/cache'
    assert (run/'artifacts/raw/proteina_complexa/run_complexa.py').is_file()


def test_native_environment_file_contains_only_explicit_paths(monkeypatch,tmp_path):
    from mn_molecular_design_studio.complexa_container import configure_environment
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv('STUDIO_TEST_CREDENTIAL','must-not-be-written')
    monkeypatch.setenv('PATH','/usr/bin')
    environment={'DATA_PATH':'/workspace/assets','SC_EXEC':'/workspace/sc'}
    for key in environment:monkeypatch.setenv(key,'')
    configure_environment(environment)
    text=(tmp_path/'.env').read_text()
    assert 'DATA_PATH="/workspace/assets"' in text
    assert 'STUDIO_TEST_CREDENTIAL' not in text and 'must-not-be-written' not in text


def _write_csv(path,rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def test_native_ligand_success_table_uses_matching_redesigned_structure(tmp_path):
    run=tmp_path/'run';raw=run/'artifacts/raw/proteina_complexa';raw.mkdir(parents=True)
    source=inputs.ASSETS/'PD-L1.pdb';shutil.copy(source,raw/'prediction.pdb')
    row={'pdb_path':'/work/artifacts/raw/proteina_complexa/generated.pdb','self_complex_pdb_path':'','mpnn_complex_pdb_path':'/work/artifacts/raw/proteina_complexa/prediction.pdb','mpnn_complex_pLDDT':'.95','mpnn_complex_i_pTM':'.81'}
    table=raw/'evaluation_results/test/RAW_ligand_binder_results_studio_complexa_combined.csv'
    _write_csv(table,[row]);_write_csv(table.parent/'all_successes_ligand_success_mpnn.csv',[row])
    assert design._proteina_complexa_result_rows(run)[0]==table
    candidates=design._normalize_proteina_complexa_candidates(run,{'pipeline_task':'ligand_binder','target_chains':[],'hotspots':'','binder_length':'80-140'},source)
    candidate=next(c for c in candidates if c['metrics'].get('result_kind')=='native_pipeline')
    assert candidate['metrics']['pass_filters'] is True
    assert candidate['metrics']['evaluated_sequence_type']=='mpnn'
    assert candidate['metrics']['binder_plddt']==.95
    assert candidate['binder_chains']==['A'] # The mock prediction contains one protein chain.


def test_raw_structures_never_become_filter_passing(tmp_path):
    run=tmp_path/'run';folder=run/'artifacts/raw/proteina_complexa/inference/test/job_0';folder.mkdir(parents=True)
    shutil.copy(inputs.ASSETS/'PD-L1.pdb',folder/'sample.pdb')
    candidates=design._normalize_proteina_complexa_candidates(run,{'target_chains':['A'],'hotspots':'','binder_length':'70-150'},inputs.ASSETS/'PD-L1.pdb')
    assert candidates and all(c['metrics'].get('pass_filters') is not True for c in candidates)


def test_ligand_chain_is_not_treated_as_a_binder(tmp_path):
    run=tmp_path/'run';raw=run/'artifacts/raw/proteina_complexa';raw.mkdir(parents=True)
    # The bundled AME structure contains protein on B and a ligand on A.
    source=inputs.ASSETS/'M0024_1nzy_v3.pdb';shutil.copy(source,raw/'prediction.pdb')
    row={'pdb_path':'/work/artifacts/raw/proteina_complexa/prediction.pdb'}
    table=raw/'evaluation_results/test/RAW_motif_binder_results_studio_complexa_combined.csv'
    _write_csv(table,[row])
    candidates=design._normalize_proteina_complexa_candidates(run,{'pipeline_task':'motif_ligand','target_chains':[],'hotspots':'','binder_length':'120-220'},source)
    candidate=next(c for c in candidates if c['metrics'].get('result_kind')=='native_pipeline')
    assert candidate['binder_chains']==['B']
    assert candidate['target_chains']==['A']
    assert candidate['metrics']['pass_filters'] is None


def test_native_passing_redesign_replaces_nonpassing_best_structure(tmp_path):
    run=tmp_path/'run';raw=run/'artifacts/raw/proteina_complexa';raw.mkdir(parents=True)
    source=inputs.ASSETS/'PD-L1.pdb'
    for name in ['best.pdb','passing.pdb']:shutil.copy(source,raw/name)
    row={'pdb_path':'/work/artifacts/raw/proteina_complexa/generated.pdb',
         'mpnn_complex_pdb_path':'/work/artifacts/raw/proteina_complexa/best.pdb',
         'mpnn_complex_pdb_path_all':str(['/work/artifacts/raw/proteina_complexa/best.pdb','/work/artifacts/raw/proteina_complexa/passing.pdb']),
         'mpnn_complex_min_ipAE_all':'[0.2, 0.01]','mpnn_binder_scRMSD_ca_all':'[3, 1]',
         'mpnn_ligand_scRMSD_aligned_allatom_all':'[9, 1]','mpnn_complex_pLDDT_all':'[0.7, 0.95]',
         'mpnn_sequence_all':str(['AAA','GGG'])}
    table=raw/'evaluation_results/test/RAW_ligand_binder_results_studio_complexa_combined.csv'
    _write_csv(table,[row]);_write_csv(table.parent/'all_successes_ligand_binder_mpnn.csv',[row])
    criteria={'thresholds':{
        'min_ipAE':{'threshold':2,'op':'<','scale':31,'column_prefix':'complex'},
        'scRMSD_ca':{'threshold':2,'op':'<','scale':1,'column_prefix':'binder'},
        'scRMSD_aligned_allatom':{'threshold':5,'op':'<','scale':1,'column_prefix':'ligand'}}}
    (table.parent/'success_criteria_ligand_binder.json').write_text(json.dumps(criteria))
    candidates=design._normalize_proteina_complexa_candidates(run,{'pipeline_task':'ligand_binder','target_chains':[]},source)
    candidate=next(c for c in candidates if c['metrics'].get('result_kind')=='native_pipeline')
    assert candidate['complex_pdb'].endswith('passing.pdb')
    assert candidate['binder_sequence']=='GGG'
    assert candidate['metrics']['binder_plddt']==.95
    assert candidate['metrics']['evaluated_redesign_index']==1
    assert candidate['metrics']['pass_filters'] is True
    row['mpnn_complex_min_ipAE_all']='[nan, 0.01]'
    selected,index,verified=design._proteina_complexa_selected_row(table,row,'mpnn','ligand_binder')
    assert verified and index==1
    (table.parent/'success_criteria_ligand_binder.json').unlink()
    candidates=design._normalize_proteina_complexa_candidates(run,{'pipeline_task':'ligand_binder','target_chains':[]},source)
    candidate=next(c for c in candidates if c['metrics'].get('result_kind')=='native_pipeline')
    from mn_molecular_design_studio.design_validation import candidate_passes_filters
    assert candidate['metrics']['native_design_pass_filters'] is True
    assert candidate_passes_filters(candidate) is False


def test_motif_success_requires_joint_criteria_on_same_redesign(tmp_path):
    table=tmp_path/'RAW_results.csv'
    row={'self_complex_pdb_path_all':str(['first.pdb','second.pdb']),
         'self_binder_scRMSD_bb3_all':'[1,3]','self_motif_rmsd_pred_all':'[2,1]',
         'self_correct_motif_sequence_all':'[1,1]','self_has_ligand_clashes_all':'[0,0]'}
    criteria={'binder_thresholds':{'scRMSD_bb3':{'threshold':2,'op':'<=','scale':1,'column_prefix':'binder'}},
              'motif_criteria':[{'column':'{seq_type}_motif_rmsd_pred_all','threshold':1.5,'op':'<='},
                                {'column':'{seq_type}_correct_motif_sequence_all','threshold':1,'op':'>='},
                                {'column':'{seq_type}_has_ligand_clashes_all','threshold':.5,'op':'<'}]}
    (tmp_path/'success_criteria_motif_binder_motif_ligand_binder.json').write_text(json.dumps(criteria))
    assert design._proteina_complexa_selected_row(table,row,'self','motif_ligand')[2] is False
    row['self_binder_scRMSD_bb3_all']='[1,1]'
    selected,index,verified=design._proteina_complexa_selected_row(table,row,'self','motif_ligand')
    assert verified and index==1 and selected['self_complex_pdb_path']=='second.pdb'


def test_multiple_redesigns_without_native_criteria_are_not_verified(tmp_path):
    row={'self_complex_pdb_path_all':str(['first.pdb','second.pdb'])}
    assert design._proteina_complexa_selected_row(tmp_path/'RAW.csv',row,'self','protein_binder')[2] is False


def test_task_form_renders_all_three_tasks_without_submission():
    from streamlit.testing.v1 import AppTest
    app=AppTest.from_string('''
from mn_molecular_design_studio.complexa_ui import render_setup
from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import workflow_schema
render_setup(workflow_schema("design.proteina_complexa"),owner_id="test",project_selector=lambda _:None,submit=lambda *a:(_ for _ in ()).throw(AssertionError("No job expected")))
''').run()
    assert not app.exception
    for task in inputs.TASKS:
        app.button_group(key='complexa-task').set_value(task).run()
        assert not app.exception
        assert any(b.label=='Submit 1 job' and b.disabled for b in app.button)
