"""Native generator runners. Downstream sequence design and folding are graph steps."""
from __future__ import annotations

import collections
import csv
import json
import re
import secrets
import shutil
from pathlib import Path

import yaml

from mn_molecular_design_studio.engine.workflow_engine.core.execution_timing import execution_phase
from mn_molecular_design_studio import generator_input as inputs
from mn_molecular_design_studio.engine.workflow_engine.core.candidates import write_candidates
from mn_molecular_design_studio.engine.workflow_engine.core.gpu import docker_gpu_args, normalize_gpu_device
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import create_job, update_status, write_json
from mn_molecular_design_studio.runtime import reference_root
from mn_molecular_design_studio.engine.workflow_engine.workflows import design


def _container_path(run, path):
    return '/work/'+str(Path(path).relative_to(run))


def _stage_structure(source, destination):
    """Keep residue identifiers and only protein atoms in a job-owned copy."""
    if Path(source).suffix.lower() in {'.cif','.mmcif'}:
        from mn_molecular_design_studio.engine.workflow_engine.workflows.refolding import _cif_to_pdb
        _cif_to_pdb(Path(source),destination)
    else:shutil.copy2(source,destination)
    # Keep motif REMARK records (Genie3 uses them), but remove other molecules.
    lines=[line for line in destination.read_text().splitlines() if line.startswith(('ATOM  ','TER','REMARK'))]
    if not any(line.startswith('ATOM  ') for line in lines):raise ValueError('The input has no readable protein coordinates.')
    destination.write_text('\n'.join([*lines,'END','']))
    return destination


def _docker(engine, run, args, mounts=(), environment=None, cwd='/work'):
    p=inputs.ENGINES[engine]
    command=['docker','run','--rm',*docker_gpu_args(run['params']['gpu_device']),'--shm-size=8G','--entrypoint','/bin/bash',
             '-v',str(run['path'])+':/work']
    for source,destination in mounts:command+=['-v',f'{source}:{destination}:ro']
    for name,value in (environment or {}).items():command+=['-e',f'{name}={value}']
    # A fixed shell wrapper activates the image environment and execs argv;
    # user paths, selections and native overrides are always separate arguments.
    command+=['-w',cwd,p['image'],'-c','exec "$@"','studio-generator',*map(str,args)]
    return command


def _pxdesign(run, target):
    p=run['params'];root=run['raw'];path=run['path'];ref=reference_root()/'pxdesign'
    target_config={'target':{'file':_container_path(path,target),'chains':{}},'binder_length':p['length_range'][0]}
    for chain in p['target_chains']:
        hotspots=[int(v[1:]) for v in p['hotspots'].split(',') if v.startswith(chain)]
        target_config['target']['chains'][chain]={'hotspots':hotspots} if hotspots else 'all'
    config=root/'input/target.yaml';config.write_text(yaml.safe_dump(target_config,sort_keys=False))
    environment={'PROTENIX_DATA_ROOT_DIR':'/ref/pxdesign/release_data/ccd_cache','OMP_NUM_THREADS':str(p['cpu_cores'])}
    mounts=[(ref,'/ref/pxdesign')]
    common=dict(environment=environment,mounts=mounts,cwd='/opt/PXDesign')
    native=_container_path(path,config)
    return [
        {'name':'PXDesign · check input','command':_docker('pxdesign',run,['pxdesign','check-input','--yaml',native],**common)},
        {'name':'PXDesign · prepare target','command':_docker('pxdesign',run,['pxdesign','parse-target','--yaml',native,'-o',_container_path(path,root/'parsed_target')],**common)},
        {'name':'PXDesign · generate backbones','command':_docker('pxdesign',run,['pxdesign','infer','-i',native,'-o',_container_path(path,root/'output'),
            '--N_sample',p['num_designs'],'--N_step',p['n_steps'],'--dtype',p['dtype'],'--use_fast_ln',str(p['use_fast_ln']),
            '--use_deepspeed_evo_attention','False','--seeds',str(p['seed']),
            '--load_checkpoint_dir','/ref/pxdesign/release_data/checkpoint'],**common)}]


def _genie3(run, target):
    p=run['params'];root=run['raw'];path=run['path'];task=p['task'];dataset=root/'dataset'
    generation={'compile':p['compile_generation'],'dataset':{'n_sample':p['num_designs']},
                'sampler':{'sampler':{'direction_scale':p['direction_scale']}}}
    config={'experiment':{'name':'studio_generator','seed':p['seed']},'paths':{'rootdir':_container_path(path,root/'output')},
            'runtime':{'num_devices':1},'generation':generation}
    if task=='binder':
        dataset,key,problem=design._write_genie3_target_files(run_dir=path,target_artifact=target,target_chains=p['target_chains'],
            binder_length=p['binder_length'],hotspots=p['hotspots'],campaign_name=p['campaign_name'])
        config['paths']['dataset']=_container_path(path,dataset)
        generation['dataset'].update(source='target',selections=key,cond_strategy=p['cond_strategy'])
        p['problem_key']=key
    elif task=='unconditional':
        lo,hi=p['length_range'];generation['dataset'].update(source='unconditional',min_length=lo,max_length=hi,length_step=p['length_step'])
    else:
        (dataset/'motifs').mkdir(parents=True,exist_ok=True);(dataset/'problems').mkdir(exist_ok=True)
        motif=dataset/'motifs/studio_motif.pdb';shutil.copy2(target,motif)
        if p.get('motif_problem'):
            problem=copy_dict(p['motif_problem'])
        else:
            blocks=inputs.residue_ranges(p['motif_residues'],inputs.inventory(target))
            headers=['REMARK 999 NAME   studio_motif','REMARK 999 PDB    STUDIO']
            for block in blocks:headers.append(f"REMARK 999 INPUT  {block[0][0]} {int(block[0][1:]):3d} {int(block[-1][1:]):3d}")
            motif.write_text('\n'.join(headers)+'\n'+''.join(line+'\n' for line in target.read_text().splitlines() if not line.startswith('REMARK 999')))
            segments=[]
            for i,region in enumerate(p['region_lengths']):
                segments.append(f"{region['min']}-{region['max']}")
                if i<len(blocks):segments.append(f'A{i+1}')
            problem={'segment_config_str':','.join(segments),'minimum_total_length':p['scaffold_length_range'][0],
                     'maximum_total_length':p['scaffold_length_range'][1]}
        problem['motif_filepaths']=[_container_path(path,motif)]
        write_json(dataset/'problems/studio_motif.json',problem)
        config['paths']['dataset']=_container_path(path,dataset)
        generation['dataset'].update(source='motif',selections='studio_motif')
    # Never add generation.inference.reward/search or an evaluation section.
    config_path=root/'experiment.yaml';config_path.write_text(yaml.safe_dump(config,sort_keys=False))
    return [{'name':'Genie3 · generate structures','command':_docker('genie3',run,
        ['genie3','generate','-c',_container_path(path,config_path),'--num-devices','1','--verbose'],
        mounts=[(reference_root()/'genie3/pretrained','/opt/genie3/pretrained')],
        environment={'OMP_NUM_THREADS':str(p['cpu_cores'])},cwd='/opt/genie3')}]


def copy_dict(value):
    return json.loads(json.dumps(value))


def _protpardelle(run, target):
    p=run['params'];path=run['path'];root=run['raw'];task=p['task'];motifs=root/'input/motifs';motifs.mkdir(exist_ok=True)
    if task=='binder':
        lengths,chain_map,residue_map=design._write_protpardelle_target_motif(target,p['target_chains'],motifs/'target.pdb')
        hotspots=design._map_protpardelle_hotspots(p['hotspots'],residue_map)
        target_contig=';/;'.join(f'{chain}1-{length}' for chain,length in lengths.items())
        contig=target_contig+f";/;{p['length_range'][0]}-{p['length_range'][1]}"
        total=[[length,length] for length in lengths.values()]+[p['length_range']]
        p.update(protpardelle_chain_map=chain_map,motif_contig=contig)
    elif task=='partial':
        shutil.copy2(target,motifs/'target.pdb');info=inputs.inventory(target)
        contig='partial_diffusion';total=[[len(c['residues']),len(c['residues'])] for c in info['chains'].values() if c['kind']=='protein'];hotspots=''
        p['target_chains']=[]
    elif task=='unconditional':contig=None;total=p['chain_length_ranges'];hotspots=None
    else:
        if target:shutil.copy2(target,motifs/'target.pdb')
        contig=p['motif_contig'];total=p['chain_length_ranges'];hotspots=p['hotspots'] or None
    payload={'search_space':{'models':[[p['model_name'],p['model_epoch'],p['sampling_config']]],
        'step_scales':[p['step_scale']],'schurns':[p['schurn']],'crop_cond_starts':[p['crop_cond_start']],'translations':[p['translation']]},
        'motifs':['target' if target else None],'motif_contigs':[contig],'total_lengths':[total],'hotspots':[hotspots],'ssadj':[None]}
    if task=='partial':payload['partial_diffusion']={'enabled':True,'rewind_steps':[p['partial_steps']]}
    config=root/'sampling.yaml';config.write_text(yaml.safe_dump(payload,sort_keys=False))
    return [{'name':'Protpardelle-1C · generate structures','command':_docker('protpardelle',run,
        ['python','-m','protpardelle.sample',_container_path(path,config),'--motif-dir',_container_path(path,motifs),
         '--num-samples',p['num_designs'],'--num-mpnn-seqs','0','--batch-size',p['batch_size'],'--seed',p['seed']],
        mounts=[(reference_root()/'protpardelle-1c/model_params','/ref/protpardelle-1c/model_params')],
        environment={'PROTPARDELLE_OUTPUT_DIR':_container_path(path,root/'output'),'PROTPARDELLE_MODEL_PARAMS':'/ref/protpardelle-1c/model_params',
                     'OMP_NUM_THREADS':str(p['cpu_cores'])},cwd='/opt/protpardelle-1c')}]


def _rfdiffusion(run, target):
    p=run['params'];root=run['raw'];path=run['path'];task=p['task'];input_dir=root/'input'
    mounts=[(reference_root()/'rfdiffusion_models','/models')]
    env={'PYTHONPATH':'/opt/RFdiffusion:/opt/RFdiffusion/env/SE3Transformer','OMP_NUM_THREADS':str(p['cpu_cores'])}
    overrides=[f'inference.output_prefix={_container_path(path,root/"output/design")}',f'inference.num_designs={p["num_designs"]}',
        'inference.model_directory_path=/models','inference.schedule_directory_path=/opt/RFdiffusion/schedules',
        f'inference.ckpt_override_path=/models/{p["resolved_model"]}_ckpt.pt',f'diffuser.T={p["timesteps"]}',
        f'denoiser.noise_scale_ca={p["noise_scale_ca"]}',f'denoiser.noise_scale_frame={p["noise_scale_frame"]}',
        f'inference.write_trajectory={str(p["save_trajectory"]).lower()}',f'inference.deterministic={str(p["deterministic"]).lower()}']
    if target:overrides.append('inference.input_pdb='+_container_path(path,target))
    if p['hotspots']:overrides.append('ppi.hotspot_res=['+p['hotspots']+']')
    if task!='fold':overrides.append('contigmap.contigs=['+p['contig']+']')
    if task=='partial':
        overrides.append(f'diffuser.partial_T={p["partial_steps"]}')
        if p['fixed_sequence_positions']:overrides.append('contigmap.provide_seq=['+p['fixed_sequence_positions']+']')
    if task in {'symmetry','symmetric_motif'}:
        symmetry='tetrahedral' if p['symmetry']=='T' else p['symmetry'].lower()
        overrides.append('inference.symmetry='+symmetry)
    potentials=[p['potential']] if p['potential'] else []
    if task in {'symmetry','symmetric_motif'} and p['oligomer_contacts']:
        potentials.append('type:olig_contacts,weight_intra:1,weight_inter:0.1')
        overrides+=['potentials.olig_intra_all=true','potentials.olig_inter_all=true']
    if potentials:overrides+=['potentials.guiding_potentials='+json.dumps(potentials,separators=(',',':')),f'potentials.guide_scale={p["guide_scale"]}','potentials.guide_decay='+p['guide_decay']]
    steps=[]
    if task=='fold':
        scaffold_dir='/models/ppi_scaffolds'
        if p['scaffold_preset']=='tim_barrel':
            local=input_dir/'scaffolds';shutil.copytree(inputs.ASSETS/'rfdiffusion/tim_barrel_scaffold',local,dirs_exist_ok=True)
            scaffold_dir=_container_path(path,local)
        elif p['scaffold_preset']=='custom':
            local=input_dir/'scaffolds';local.mkdir(exist_ok=True);scaffold_dir=_container_path(path,local)
            if p.get('scaffold_pdb'):
                source=_stage_structure(p['scaffold_pdb'],input_dir/'scaffold.pdb')
                steps.append({'name':'RFdiffusion · prepare binder fold','command':_docker('rfdiffusion',run,
                    ['python3','/opt/RFdiffusion/helper_scripts/make_secstruc_adj.py','--input_pdb',_container_path(path,source),'--out_dir',scaffold_dir],mounts=mounts,environment=env)})
            else:
                shutil.copy2(p['scaffold_ss_file'],local/'scaffold_ss.pt');shutil.copy2(p['scaffold_adj_file'],local/'scaffold_adj.pt')
        else:
            files=sorted((reference_root()/'rfdiffusion_models/ppi_scaffolds').glob('*_ss.pt'))
            if not files:raise ValueError('The RFdiffusion scaffold library is empty.')
            selected=files[p['seed']%len(files)].name.removesuffix('_ss.pt')
            overrides.append('scaffoldguided.scaffold_list=['+selected+']');p['scaffold_id']=selected
        overrides+=['scaffoldguided.scaffoldguided=true',f'scaffoldguided.target_pdb={str(bool(target)).lower()}',
            'scaffoldguided.scaffold_dir='+scaffold_dir,f'scaffoldguided.sampled_insertion=0-{p["sampled_insertion"]}',
            f'scaffoldguided.sampled_N=0-{p["sampled_n"]}',f'scaffoldguided.sampled_C=0-{p["sampled_c"]}']
        if target:overrides.append('scaffoldguided.target_path='+_container_path(path,target))
        for key,native in [('target_ss_file','target_ss'),('target_adj_file','target_adj')]:
            if p.get(key):
                copied=input_dir/(native+'.pt');shutil.copy2(p[key],copied)
                overrides.append('scaffoldguided.'+native+'='+_container_path(path,copied))
    # RFdiffusion has no general --seed flag. Seed all three native RNGs before
    # running its entrypoint; deterministic mode also uses this design index.
    python_seed="import random, numpy as np, torch, runpy, sys; seed=int(sys.argv.pop(1)); random.seed(seed); np.random.seed(seed); torch.manual_seed(seed); sys.argv[0]='/opt/RFdiffusion/scripts/run_inference.py'; runpy.run_path(sys.argv[0],run_name='__main__')"
    if p['deterministic']:overrides.append(f'inference.design_startnum={p["seed"]}')
    write_json(root/'native_overrides.json',overrides)
    steps.append({'name':'RFdiffusion · generate backbones','command':_docker('rfdiffusion',run,['python3','-c',python_seed,p['seed'],*overrides],mounts=mounts,environment=env,cwd='/opt/RFdiffusion')})
    return steps


def _residues(path):
    result=collections.OrderedDict()
    for line in Path(path).read_text().splitlines():
        if line.startswith('ATOM'):
            label=line[21:22]+line[22:26].strip()
            result.setdefault(label,line[17:20].strip())
    return result


def _fixed_residues(engine, structure, p):
    if not p['fixed_motif_required']:return {}
    fixed=[];source_labels=p.get('fixed_identity_labels') or [];expected=[]
    source=_residues(p['staged_target'])
    if engine=='rfdiffusion':
        fixed=list(dict.fromkeys(line[21:22]+line[22:26].strip() for line in structure.read_text().splitlines()
                               if line.startswith('ATOM') and float(line[60:66] or 0)>.5))
    elif engine=='genie3':
        fixed=list(dict.fromkeys(line[21:22]+line[22:26].strip() for line in structure.read_text().splitlines()
                               if line.startswith('ATOM') and line[72:76].strip()))
    else:
        table=structure.parent/'scaffold_info.csv'
        if table.is_file():
            index=int(structure.stem.rsplit('_',1)[-1])
            row=next((r for r in csv.DictReader(table.open()) if int(r['sample_num'])==index),{})
            placement=row.get('motif_placements','')
            observed=_residues(structure);chains=list(dict.fromkeys(k[0] for k in observed))
            for ci,(original,placed) in enumerate(zip(p['motif_contig'].split(';/;'),placement.split(';/;'))):
                motif_blocks=[part for part in original.split(';') if part and part[0].isalpha()]
                cursor=0;mi=0;chain_residues=[k for k in observed if k[0]==chains[ci]]
                for token in placed.split('/'):
                    if token.isdigit():cursor+=int(token)
                    elif token:
                        block=inputs.residue_ranges(motif_blocks[mi],inputs.inventory(Path(p['staged_target'])))[0];mi+=1
                        fixed+=chain_residues[cursor:cursor+len(block)];expected += [source[label] for label in block];cursor+=len(block)
    observed=_residues(structure)
    fixed=[label for label in fixed if observed.get(label)!='UNK']
    expected=[name for name in expected if name!='UNK']
    if not expected:expected=[source[label] for label in source_labels if label in source]
    valid=bool(fixed) and len(fixed)==len(set(fixed)) and all(label in observed for label in fixed)
    if expected:valid=valid and collections.Counter(expected)==collections.Counter(observed[label] for label in fixed)
    metadata={'fixed_motif_required':True,'fixed_residues':fixed,'motif_mapping_verified':valid}
    if not valid:metadata['motif_mapping_warning']='The native output did not provide a complete motif identity mapping. Sequence redesign is blocked until the motif mapping is available.'
    return metadata


def _pxdesign_chain_roles(target, structure, requested_target_chains):
    """Keep PXDesign's ``xpb`` placeholder binder chain in candidate roles.

    The generic protein CIF reader intentionally filters non-standard residue
    names, which also hides PXDesign's generated ``xpb`` backbone residues.
    Read all model chains here, then identify the target by sequence so the
    placeholder chain remains available for later sequence redesign.
    """
    import gemmi

    parsed = gemmi.read_structure(str(structure))
    if not len(parsed):
        raise ValueError(f'PXDesign output has no structural models: {structure}')
    output_chains = [str(chain.name).strip() or '_' for chain in parsed[0]]
    target_chains = design._target_chains_in_generated_structure(
        target, structure, requested_target_chains, output_chains
    ) if target is not None else []
    inference = 'sequence_match' if target_chains else 'requested_chain_ids'
    if not target_chains:
        target_chains = [chain for chain in requested_target_chains if chain in output_chains]
    if not target_chains and target is not None and len(output_chains) > 1:
        # PXDesign preserves target-first chain ordering even when it rewrites
        # chain identifiers. Keep that ordering as a last-resort fallback.
        target_chains = output_chains[:1]
        inference = 'fallback_first_generated_chain'
    target_set = set(target_chains)
    return [chain for chain in output_chains if chain not in target_set], target_chains, inference


def normalize_outputs(run, target):
    engine=run['engine'];p=run['params'];root=run['raw'];path=run['path']
    # PXDesign writes the converted input as ``output/target.cif``; it has no
    # coordinate model and is not a generated design. Only CIFs in a native
    # predictions directory are model outputs.
    if engine=='pxdesign':paths=sorted((root/'output').glob('**/predictions/*.cif'))
    elif engine=='genie3':paths=sorted((root/'output').glob('**/pdbs/*.pdb'))
    elif engine=='rfdiffusion':paths=sorted((root/'output').glob('design_*.pdb'))
    else:paths=sorted((root/'output').glob('**/sample_*.pdb'))
    rows=[]
    for index,structure in enumerate(paths,1):
        if engine=='pxdesign' and target and p['task']=='binder':
            binders,targets,chain_inference=_pxdesign_chain_roles(target,structure,p['target_chains'])
        elif target and p['task'] in {'binder','fold','redesign','partial','multichain'} and p['target_chains']:
            binders,targets,chain_inference=design._infer_generated_chain_roles(target,structure,p['target_chains'])
        else:
            targets=[];binders=design._chain_ids_from_structure(structure);chain_inference='generated_chains'
        try:motif_metadata=_fixed_residues(engine,structure,p)
        except (ValueError,KeyError,IndexError,OSError) as exc:
            motif_metadata={'fixed_motif_required':True,'fixed_residues':[],'motif_mapping_verified':False,
                'motif_mapping_warning':f'Could not read the native motif mapping: {exc}. Sequence redesign is blocked for this structure.'}
        raw={'generator_only':True,'generator_engine':engine,'generation_task':p['task'],'result_kind':'generated',
             'target_entity_kind':'protein' if targets else None,'fixed_motif_required':False,
             'input_target_chains':p['target_chains'],'target_chain_inference':chain_inference,
             'conditioned_motif_residues':p.get('motif_labels',[]),**motif_metadata}
        rows.append({'candidate_id':engine+f'_generated_{index:05d}','source_tool':inputs.TOOLS[engine],
            'stage':'generation.backbone','target_pdb':str(target.relative_to(path)) if targets and target else None,
            'complex_pdb':str(structure.relative_to(path)),'binder_chains':binders,'target_chains':targets,'binder_sequence':None,
            'binder_length':p['binder_length'],'contig':p.get('contig') or p.get('motif_contig'),
            'hotspots':p['hotspots'].split(',') if p['hotspots'] else [],
            'metrics':{'result_kind':'generated','pass_filters':None},'raw_metadata':raw})
    return write_candidates(path,inputs.TOOLS[engine],rows)


def _run(engine, parameters):
    p=inputs.prepare_request(engine,parameters)
    ready=inputs.readiness(engine,p)
    if not ready['ready']:raise ValueError('Missing generator reference files: '+', '.join(ready['missing']))
    p['seed']=p['seed'] or secrets.randbelow(2147483646)+1
    p['gpu_device']=normalize_gpu_device(p['gpu_device'])
    p['engine_version']=inputs.ENGINES[engine]['version']
    job=create_job('design','design_generation',inputs.TOOLS[engine],{'target_pdb':p.get('target_pdb'),'target_chains':p['target_chains']},p)
    root=job.run_dir/'artifacts/raw'/engine;(root/'input').mkdir(parents=True,exist_ok=True)
    target=_stage_structure(p['target_pdb'],root/'input/target.pdb') if p.get('target_pdb') else None
    if target:p['staged_target']=str(target)
    run={'engine':engine,'params':p,'path':job.run_dir,'raw':root}
    update_status(job.run_dir,'preparing',campaign_name=p.get('campaign_name') or inputs.ENGINES[engine]['title'],current_phase='Preparing generator inputs')
    steps={'pxdesign':_pxdesign,'genie3':_genie3,'protpardelle':_protpardelle,'rfdiffusion':_rfdiffusion}[engine](run,target)
    write_json(job.run_dir/'input.json',{'inputs':{'target_pdb':p.get('target_pdb'),'target_chains':p['target_chains']},'params':p})
    rc=design._run_shell_steps(job.run_dir,steps)
    candidates=[]
    if rc==0:
        with execution_phase(job.run_dir, "candidate_normalization"):
            candidates=normalize_outputs(run,target)
    design._finish_design_job(job.run_dir,rc==0 and bool(candidates),rc,candidates,
        [(f'raw/{engine}/**/*.{ext}',ext) for ext in ('pdb','cif','yaml','json','csv','trb','pt','log')])
    return job.run_dir


def run_pxdesign_generation(target_pdb:Path|None=None,target_chains:list[str]|None=None,binder_length:str='80',hotspots:str='',campaign_name:str='',
    num_designs:int=10,n_steps:int=400,dtype:str='bf16',use_fast_ln:bool=True,seed:int=0,gpu_device:object='all',cpu_cores:int=4,task:str='binder')->Path:
    return _run('pxdesign',dict(locals()))


def run_genie3_generation(target_pdb:Path|None=None,target_chains:list[str]|None=None,binder_length:str='80-120',hotspots:str='',campaign_name:str='',
    num_designs:int=5,task:str='binder',direction_scale:float=0.,cond_strategy:str='extended',compile_generation:bool=False,
    length_step:int=1,motif_residues:str='',region_lengths:list[dict]|None=None,motif_example:str='',seed:int=0,gpu_device:object='all',cpu_cores:int=4)->Path:
    return _run('genie3',dict(locals()))


def run_protpardelle_generation(target_pdb:Path|None=None,target_chains:list[str]|None=None,binder_length:str='70-100',hotspots:str='',campaign_name:str='',
    num_designs:int=10,task:str='binder',model_name:str='cc83',step_scale:float=1.2,schurn:int=200,crop_cond_start:float=0.,batch_size:int=1,
    partial_steps:int=100,motif_contig:str='',chain_lengths:str='128',motif_residues:str='',region_lengths:list[dict]|None=None,
    translation:list[float]|None=None,seed:int=0,gpu_device:object='all',cpu_cores:int=4)->Path:
    p=dict(locals());p['translation']=translation if translation is not None else [0.,0.,0.]
    return _run('protpardelle',p)


def run_rfdiffusion_generation(target_pdb:Path|None=None,target_chains:list[str]|None=None,binder_length:str='70-100',hotspots:str='',campaign_name:str='',
    num_designs:int=10,task:str='binder',contig:str='',motif_residues:str='',designed_residues:str='',designed_chain:str='',region_lengths:list[dict]|None=None,
    timesteps:int=50,partial_steps:int=20,fixed_sequence_positions:str='',model_weights:str='auto',noise_scale_ca:float=1.,noise_scale_frame:float=1.,
    potential:str='',guide_scale:float=2.,guide_decay:str='quadratic',symmetry:str='C4',total_length:int=320,flank_before:int=50,flank_after:int=50,
    oligomer_contacts:bool=True,deterministic:bool=False,save_trajectory:bool=False,scaffold_pdb:Path|None=None,
    scaffold_ss_file:Path|None=None,scaffold_adj_file:Path|None=None,target_ss_file:Path|None=None,target_adj_file:Path|None=None,
    scaffold_preset:str='library',sampled_insertion:int=0,sampled_n:int=0,sampled_c:int=0,seed:int=0,gpu_device:object='all',cpu_cores:int=4)->Path:
    return _run('rfdiffusion',dict(locals()))
