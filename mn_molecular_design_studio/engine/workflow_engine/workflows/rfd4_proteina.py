"""RFD4 task workflow with native input validation and selectable refolding."""
from __future__ import annotations

from contextlib import nullcontext
import csv
import json
import os
import secrets
import re
import shutil
import traceback
from pathlib import Path
from typing import Literal

from Bio.PDB import MMCIFIO, Select

from mn_molecular_design_studio.engine.workflow_engine.core.execution_timing import execution_phase
from mn_molecular_design_studio import rfd4_proteina_input as inputs
from mn_molecular_design_studio.design_validation import read_structure
from mn_molecular_design_studio.engine.workflow_engine.core.candidates import (
    STAGE_GENERATION_BACKBONE_SEQUENCE, STAGE_COMPLEX_REFOLDING,
    STAGE_MONOMER_REFOLDING, STAGE_ANALYSIS, write_candidates,
)
from mn_molecular_design_studio.engine.workflow_engine.core.docker_runner import run_docker_job
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import create_job, finish_job, update_status, write_json
from mn_molecular_design_studio.engine.workflow_engine.core.manifests import load_manifest
from mn_molecular_design_studio.engine.workflow_engine.runtime import reference_root
from mn_molecular_design_studio.engine.workflow_engine.workflows.rfd4_oracle import OpenFold3Bridge

Task = Literal['protein_binder','nanobody','antibody','tcr','motif_scaffolding','enzyme',
               'small_molecule_binder','dna_binder','rna_design','existing_structure',
               'custom_conditions','score_structures','unconditional']
InputMode = Literal['guided','guided_binder','annotated_cif','condition_spec','unconditional']

def sequence_structure(text: str, destination: Path) -> Path:
    """Create a sequence carrier; none of its placeholder coordinates are conditioned."""
    import numpy as np
    from Bio.PDB import Structure, Model, Chain, Residue, Atom
    from Bio.SeqUtils import seq3
    lines=text.strip().splitlines()
    if any(line.startswith('>') for line in lines):
        from Bio import SeqIO
        from io import StringIO
        sequences=[str(record.seq).upper() for record in SeqIO.parse(StringIO(text),'fasta')]
    else:
        sequences=[re.sub(r'\s+','',part).upper() for part in text.split(':')]
    if not sequences or len(sequences)>26 or any(not re.fullmatch('[ACDEFGHIKLMNPQRSTVWY]{1,1024}',s) for s in sequences):
        raise ValueError('Provide 1–26 protein chains with 1–1024 standard amino acids each, as FASTA or colon-separated sequences.')
    structure=Structure.Structure('sequence'); model=Model.Model(0); structure.add(model)
    serial=1
    for c,sequence in zip('ABCDEFGHIJKLMNOPQRSTUVWXYZ',sequences):
        chain=Chain.Chain(c); model.add(chain)
        for i,letter in enumerate(sequence,1):
            residue=Residue.Residue((' ',i,' '),seq3(letter).upper(),' '); chain.add(residue)
            for atom,offset,element in [('N',0.,'N'),('CA',1.4,'C'),('C',2.8,'C'),('O',3.3,'O')]:
                residue.add(Atom.Atom(atom,np.array([i*3.8+offset,0.,0.]),0.,1.,' ',atom,serial,element=element)); serial+=1
    destination.parent.mkdir(parents=True,exist_ok=True)
    writer=MMCIFIO(); writer.set_structure(structure); writer.save(str(destination))
    return destination


def local_path(run_dir, value):
    path=Path(str(value))
    if path.is_relative_to('/work'): return run_dir/path.relative_to('/work')
    if not path.is_absolute(): return run_dir/path
    path.resolve().relative_to(run_dir.resolve())
    return path


def normalize_outputs(run_dir, prepared, task, requested):
    root=prepared['root']; index_path=root/'result_index.json'
    if not index_path.is_file(): return []
    index=json.loads(index_path.read_text()); records=[]
    metric_rows=[]
    for value in index.get('metrics_files',[]):
        with local_path(run_dir,value).open() as handle: metric_rows.extend(csv.DictReader(handle))
    for i,item in enumerate(index.get('structures',[]),1):
        path=local_path(run_dir,item['structure'])
        if not path.is_file(): raise ValueError('A reported RFD4 output structure is missing.')
        binder=list(item.get('designed_chains') or []); targets=list(item.get('target_chains') or [])
        parent=item.get('source_candidate') or {}
        if parent:
            binder=parent.get('binder_chains') or binder;targets=parent.get('target_chains') or []
            if set(binder+targets)-set(item.get('sequences',{})):raise ValueError('RFD4 folding changed input chain identities; cannot preserve target/design roles.')
        sequence=''.join(item.get('sequences',{}).get(c,'') for c in binder) or None
        name=path.stem.removesuffix('_generated')
        ident=f'rfd4_{i:06d}_{name}'
        binder_path=path
        if binder:
            try:
                class Designed(Select):
                    def accept_chain(self,chain): return chain.id in binder
                binder_path=root/'normalized_structures'/f'{i:06d}_designed.cif'
                binder_path.parent.mkdir(exist_ok=True)
                writer=MMCIFIO(); writer.set_structure(read_structure(path)); writer.save(str(binder_path),Designed())
            except (ValueError,KeyError,StopIteration):
                binder_path=path
        native=[r for r in metric_rows if str(r.get('designed_structure_path',''))==item['structure']]
        if not native: native=[r for r in metric_rows if str(r.get('example_id','')) and str(r['example_id']) in str(path.parent)]
        if len(metric_rows)==1: native=metric_rows
        folding=task=='existing_structure' and prepared['parameters'].get('options',{}).get('operation')=='structure_prediction'
        stage=STAGE_ANALYSIS if task=='score_structures' else (STAGE_COMPLEX_REFOLDING if targets else STAGE_MONOMER_REFOLDING) if folding else STAGE_GENERATION_BACKBONE_SEQUENCE
        metrics={'generated_residue_count':item.get('residue_count'),'result_kind':'structure_prediction' if folding else 'structure_scoring' if task=='score_structures' else 'sequence_structure_codesign'}
        if len(native)==1:
            for key,value in native[0].items():
                if value in ('','nan','NaN'): continue
                try: metrics[key]=float(value)
                except (ValueError,TypeError): metrics[key]=value
        from mn_molecular_design_studio.molecule_imports import parse_structure
        from mn_molecular_design_studio.scientific_contracts import structure_chemistry
        unsupported_chemistry = structure_chemistry(path)
        components=parse_structure(path.read_bytes(),path.name,include_ligands=True)[0]['components']
        if task in {'unconditional','protein_binder','small_molecule_binder','dna_binder'} and sequence:
            if not prepared['parameters'].get('length_min',60)<=len(sequence)<=prepared['parameters'].get('length_max',120):
                raise ValueError('RFD4 returned a designed chain outside the requested length range.')
        target_reference=prepared.get('target_structure')
        target_artifact=str(local_path(run_dir,target_reference)) if target_reference else None
        records.append({'candidate_id':ident,'source_tool':'rfd4_proteina','stage':stage,
                        'parents':[parent['candidate_id']] if parent.get('candidate_id') else [],
                        'complex_pdb':str(path),'binder_pdb':str(binder_path),'binder_sequence':sequence,
                        'target_pdb':target_artifact,'target_chains':targets,'binder_chains':binder,
                        'metrics':metrics,'raw_metadata':{'task':task,'native_structure':str(path),
                        'source_candidate':parent,
                        'sequences_by_chain':item.get('sequences',{}),'engine_version':inputs.UPSTREAM_VERSION,
                        'upstream_commit':inputs.UPSTREAM_COMMIT,'requested_designs':requested,
                        'studio_molecule_id':ident,'components':components,'unsupported_chemistry':unsupported_chemistry,
                        'reference_constraints':{key:prepared['parameters'].get(key) for key in ('task','condition_spec','options')},
                        'folding_engine':prepared['parameters'].get('options',{}).get('folding_engine','openfold3_jax') if prepared['parameters'].get('options',{}).get('refold',True) else None}})
    return records


def run_rfd4_proteina_design(
    task: Task='protein_binder', input_mode: InputMode='guided', target_structure: Path|None=None,
    target_chains: list[str]|None=None, hotspots: str='', length_min: int=60, length_max: int=120,
    condition_spec: str='', input_files: list[Path]|None=None,
    avoid_structures: list[Path]|None=None, also_bind_structures: list[Path]|None=None,
    num_designs: int=8, batch_size: int=1, nsteps: int=60, seed: int|None=None,
    options: dict|None=None, compile_model: bool=False, cpu_cores: int=8,
    gpu_device: str='0', job_name: str='rfd4-proteina-design',
    protein_sequences: str='',
    task_group: Literal['design','refolding']='design',
    distributed_campaign: dict|None=None,
) -> Path:
    if isinstance(options,dict) and options.get('native_phase_timing') is True:
        if (isinstance(distributed_campaign,dict) and distributed_campaign.get('enabled')) or task=='score_structures':
            raise ValueError('Native phase timing supports ordinary local generation only; campaign and scoring execution are unchanged.')
    parameters=json.loads(json.dumps(dict(locals()),default=str)); job=None
    try:
        if seed is None: parameters['seed']=secrets.randbelow(2**32)
        parameters['options']=dict(options or {})
        campaign_config=dict(distributed_campaign or {})
        distributed_enabled=bool(campaign_config.get('enabled'))
        if distributed_enabled:
            if task_group!='design' or task!='protein_binder' or input_mode not in {'guided','guided_binder'}:
                raise ValueError('RFD4 acceptance-driven campaigns currently support guided protein-binder co-design only.')
            if not target_structure or not target_chains or input_files or avoid_structures or also_bind_structures or protein_sequences:
                raise ValueError('RFD4 campaigns require one selected target structure and chain set without secondary inputs.')
            campaign_config['generator_id']='rfd4_proteina'
            campaign_config.setdefault('seed',int(parameters['seed']))
            campaign_resources=dict(campaign_config.get('resources') or {})
            campaign_resources.setdefault('cpu_cores',int(cpu_cores or 8))
            campaign_config['resources']=campaign_resources
            from mn_molecular_design_studio.engine.workflow_engine.workflows.distributed_campaign import validate_campaign_config
            campaign_config=validate_campaign_config(campaign_config)
            parameters['distributed_campaign']=campaign_config
        job_type=('rfd4_proteina_campaign' if distributed_enabled else
                  'rfd4_structure_prediction' if task_group=='refolding' else 'rfd4_proteina_design')
        job=create_job(task_group,job_type,'rfd4_proteina',
                       {'target_structure':str(target_structure) if target_structure else None,
                        'input_files':[str(p) for p in input_files or []]},parameters)
        update_status(job.run_dir,'preparing',campaign_name=job_name,current_phase='Validating RFD4 task and inputs')
        if protein_sequences:
            if target_structure: raise ValueError('Provide either sequence input or a structure.')
            if task!='existing_structure' or parameters['options'].get('operation')!='structure_prediction':
                raise ValueError('Sequence input is available for structure prediction.')
            parameters['target_structure']=str(sequence_structure(protein_sequences,job.run_dir/'artifacts/sequence_input.cif'))
        design_parameters={key:value for key,value in parameters.items() if key!='distributed_campaign'}
        prepared=inputs.prepare_request(design_parameters,run_dir=job.run_dir)
        if distributed_enabled:
            write_json(job.run_dir/'parameters.json',{**parameters,'engine_version':inputs.UPSTREAM_VERSION,
                       'upstream_commit':inputs.UPSTREAM_COMMIT,'checkpoint_revision':inputs.HF_REVISION,'image':inputs.IMAGE})
            from mn_molecular_design_studio.engine.workflow_engine.workflows.distributed_campaign import run_rfd4_proteina_distributed_campaign
            return run_rfd4_proteina_distributed_campaign(job.run_dir)
        # Stage the exact adapter with every run for reproducibility and remote execution.
        runtime=prepared['root']/'runtime'; runtime.mkdir(exist_ok=True)
        runtime_files=['rfd4_runtime.py','rfd4_studio_oracle.py']
        if parameters['options'].get('native_phase_timing'):
            runtime_files.append('rfd4_native_diagnostics.py')
        for name in runtime_files:
            source=Path(__file__).with_name(name)
            if source.exists(): shutil.copy2(source,runtime/name)
        write_json(job.run_dir/'parameters.json',{**parameters,'engine_version':inputs.UPSTREAM_VERSION,
                   'upstream_commit':inputs.UPSTREAM_COMMIT,'checkpoint_revision':inputs.HF_REVISION,'image':inputs.IMAGE})
        cache=Path(reference_root())/'rfd4-proteina/huggingface'; cache.mkdir(parents=True,exist_ok=True)
        temp=prepared['root']/'tmp'; temp.mkdir(exist_ok=True)
        o=parameters['options']
        needs_oracle=o.get('refold',True) or o.get('apo') or o.get('reward_search')=='Best by refolding' or avoid_structures or also_bind_structures
        bridge=OpenFold3Bridge(job.run_dir,gpu_device,cpu_cores,options=o) if needs_oracle else nullcontext()
        with bridge as oracle:
            env={'HF_HOME':'/cache/huggingface','HF_HUB_CACHE':'/cache/huggingface/hub',
                 'OMP_NUM_THREADS':str(cpu_cores),'PYTHONUNBUFFERED':'1',
                 'TMPDIR':'/work/'+temp.relative_to(job.run_dir).as_posix(),
                 'PYTHONPATH':'/work/'+runtime.relative_to(job.run_dir).as_posix()}
            if o.get('native_phase_timing'):
                env['MN_STUDIO_RFD4_NATIVE_TIMING_FILE']='/work/'+(prepared['root']/'native_lifecycle.jsonl').relative_to(job.run_dir).as_posix()
            if os.environ.get('HF_HUB_OFFLINE'):
                env['HF_HUB_OFFLINE'] = os.environ['HF_HUB_OFFLINE']
            mounts=[(cache,'/cache/huggingface',False)]
            if oracle:
                mounts.append(oracle.mount)
                env.update(AF3_LAUNCHER_SOCKET_TEMPLATE='/oracle/oracle.sock',
                           AF3_MODEL_DIR=('/cache/huggingface/hub/models--nvidia--RFD4-Proteina/snapshots/'+inputs.HF_REVISION+'/openfold3/weights') if oracle.engine=='openfold3_jax' else '/studio-prediction/'+oracle.engine)
            manifest={**load_manifest('rfd4_proteina'),'entrypoint':'python'}
            args=['/work/'+(runtime/'rfd4_runtime.py').relative_to(job.run_dir).as_posix(),
                  '/work/'+prepared['request_path'].relative_to(job.run_dir).as_posix()]
            code=run_docker_job(job.run_dir,manifest,args,gpu_device=gpu_device,extra_mounts=mounts,environment=env)
        native_commands=prepared['root']/'native_command.json'
        receipt_path=job.run_dir/'native_execution.json'
        if receipt_path.is_file() and native_commands.is_file():
            receipt=json.loads(receipt_path.read_text())
            commands=json.loads(native_commands.read_text())
            receipt['processes'][-1]['native_invocations']=len(commands) if all(c.get('status')=='returned' for c in commands) else None
            native_counts=[p.get('native_invocations') for p in receipt['processes']]
            receipt.update(native_invocations=sum(native_counts) if all(n is not None for n in native_counts) else None,
                           native_commands=commands, prediction_outputs=0,
                           native_count_basis='RFD4 wrapper process receipts; oracle processes recorded separately')
            write_json(receipt_path,receipt)
        if code:
            detail=''
            oracle_log=prepared['root']/'oracle.log'
            if oracle_log.is_file():
                lines=[line for line in oracle_log.read_text(errors='replace').splitlines() if line.strip()]
                if lines: detail=' '+lines[-1][:500]
            raise RuntimeError(f'RFD4 exited with status {code}.{detail} See the run logs for details.')
        with execution_phase(job.run_dir, "candidate_normalization"):
            candidates=normalize_outputs(job.run_dir,prepared,task,num_designs)
            normalized=write_candidates(job.run_dir,'rfd4_proteina',candidates) if candidates else []
        if receipt_path.is_file():
            receipt=json.loads(receipt_path.read_text())
            receipt['prediction_outputs']=len(normalized) if o.get('operation')=='structure_prediction' else 0
            write_json(receipt_path,receipt)
        table=prepared['root']/'designs.csv'
        with table.open('w',newline='') as handle:
            writer=csv.DictWriter(handle,fieldnames=['candidate_id','task','binder_sequence','structure'])
            writer.writeheader()
            writer.writerows({'candidate_id':x['candidate_id'],'task':task,'binder_sequence':x.get('binder_sequence'),
                              'structure':x['complex_pdb']} for x in candidates)
        if o.get('cluster') and normalized:
            from .quality import run_foldseek
            cluster_run=run_foldseek(job.run_dir,job.run_dir/'artifacts/normalized_candidates/candidates.jsonl',cpu_cores=cpu_cores)
            cluster_result=json.loads((cluster_run/'result.json').read_text())
            if not cluster_result.get('success'): raise RuntimeError('Design clustering failed: '+str(cluster_result.get('error') or cluster_run))
            write_json(prepared['root']/'clustering_run.json',{'run_dir':str(cluster_run)})
        if not normalized and not o.get('filter_conditions',True): raise RuntimeError('RFD4 returned no usable structures.')
        artifacts=[{'path':str(p),'name':p.name,'type':'rfd4_proteina_output'} for p in prepared['root'].rglob('*') if p.is_file() and not p.is_socket()]
        finish_job(job.run_dir,True,{'outputs':{'artifacts':artifacts,'candidates':normalized,'design_table':str(table)},
                   'metrics':{'task':task,'candidate_count':len(normalized),'requested_designs':num_designs,
                              'input_design_count':prepared['design_count'],'filter_conditions':o.get('filter_conditions',True),
                              'checkpoint_revision':inputs.HF_REVISION,'folding_engine':o.get('folding_engine','openfold3_jax') if needs_oracle else 'rfd4_proteina' if o.get('operation')=='structure_prediction' else None},
                   'downstream_artifacts':{'candidates_jsonl':'artifacts/normalized_candidates/candidates.jsonl'}})
    except Exception as exc:
        if job is not None:
            with (job.run_dir/'stderr.log').open('a') as handle: handle.write(traceback.format_exc())
            finish_job(job.run_dir,False,{'error':str(exc),'metrics':{'task':task,'error':str(exc)}})
    return job.run_dir if job is not None else Path()


def run_rfd4_structure_prediction(target_structure: Path|None=None, protein_sequences: str='', input_files: list[Path]|None=None, target_chains: list[str]|None=None,
                                   num_designs: int=1, seed: int|None=None, nsteps: int=60,
                                   options: dict|None=None, cpu_cores: int=8, batch_size: int=1, compile_model: bool=False, gpu_device: str='0',
                                   job_name: str='rfd4-structure-prediction') -> Path:
    """RFD4's own sequence-conditioned folding, without an external predictor."""
    settings={**(options or {}),'operation':'structure_prediction','refold':False,'filter_conditions':False}
    if target_chains: settings['design_chains']=target_chains
    return run_rfd4_proteina_design(task='existing_structure',target_structure=target_structure,
                                   num_designs=num_designs,seed=seed,nsteps=nsteps,options=settings,
                                   cpu_cores=cpu_cores,gpu_device=gpu_device,job_name=job_name,
                                   protein_sequences=protein_sequences,input_files=input_files,batch_size=batch_size,compile_model=compile_model,
                                   task_group='refolding')
