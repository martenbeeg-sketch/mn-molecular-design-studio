"""Public, generation-only task contracts shared by tool pages and canvas nodes."""
from __future__ import annotations

import copy
import json
import math
import re
from pathlib import Path

from mn_molecular_design_studio.rfdiffusion3_input import inventory as _inventory
from mn_molecular_design_studio.runtime import reference_root

ASSETS = Path(__file__).with_name('assets') / 'generators'
ENGINES = {
    'pxdesign': {'title':'PXDesign', 'version':'PXDesign f788441', 'image':'mn-pxdesign:cu128', 'repo':'https://github.com/bytedance/PXDesign',
                 'tasks':{'binder':'Binder Design'}, 'description':'Generate protein binder backbones with diffusion models.'},
    'rfdiffusion': {'title':'RFdiffusion', 'version':'RFdiffusion · version not pinned', 'image':'mn-rfdiffusion:latest', 'repo':'https://github.com/RosettaCommons/RFdiffusion',
                   'tasks':{'custom':'Custom Contigs','binder':'Binder Design','redesign':'Binder Redesign','motif':'Motif Scaffolding','partial':'Partial Diffusion','fold':'Fold Conditioning','symmetry':'Symmetric Oligomer','symmetric_motif':'Symmetric Motif Scaffolding'}, 'description':'Design or diversify proteins with target, motif, fold, or symmetry conditioning.'},
    'genie3': {'title':'Genie3', 'version':'Genie3 d77ae5a', 'image':'mn-genie3:cu128', 'repo':'https://github.com/aqlaboratory/genie3',
               'tasks':{'binder':'Binder Design','unconditional':'Unconditional Generation','motif':'Motif Scaffolding'}, 'description':'Generate protein structures with fast all-atom diffusion.'},
    'protpardelle': {'title':'Protpardelle-1C', 'version':'Protpardelle-1C ee37840', 'image':'mn-protpardelle-1c:cu128', 'repo':'https://github.com/ProteinDesignLab/protpardelle-1c',
                    'tasks':{'binder':'Binder Design','unconditional':'Unconditional Generation','motif':'Motif Scaffolding','partial':'Partial Diffusion','multichain':'Multichain Generation'}, 'description':'Generate protein structures with sidechain, motif, and multichain conditioning.'},
}
WORKFLOWS = {key:'design.'+key+'_generate' for key in ENGINES}
TOOLS = {key:key+'_generator' for key in ENGINES}
RF_MODELS = ['auto','Base','Complex_base','Complex_beta','ActiveSite','InpaintSeq','Complex_Fold_base','Base_epoch8']
PROT_MODELS = {'cc58':'416','cc78':'1431','cc83':'2616','cc89':'415','cc91':'383','cc94':'3100','cc95':'3490'}
SYMMETRIES = {'C2':2,'C3':3,'C4':4,'C5':5,'C6':6,'D2':4,'D3':6,'D4':8,'T':12}


def inventory(path):
    info=_inventory(path)
    # Official MotifBench inputs include backbone-only UNK residues. Their
    # geometry is conditioned, while their amino-acid identity remains free.
    for row in info['residues'].values():
        if row['name']=='UNK' and {'N','CA','C'}<=set(row['atoms']):row['kind']='protein'
    for row in info['residues'].values():
        if row['kind']=='protein':
            chain=row['label'][:-len(str(row['number']))] if not row['insertion_code'] else row['label'][:-(len(str(row['number']))+len(row['insertion_code']))]
            info['chains'].setdefault(chain,{'kind':'protein','residues':[]})
    for chain in list(info['chains']):
        if info['chains'][chain]['kind']=='protein':
            info['chains'][chain]['residues']=[r for r in info['residues'].values() if r['label']==chain+str(r['number'])+r['insertion_code'] and r['kind']=='protein']
    return info


def default_parameters(engine, task='binder'):
    p = {'task':task,'target_pdb':None,'target_chains':[],'hotspots':'','campaign_name':'',
         'binder_length':'80' if engine=='pxdesign' else '80-120' if engine=='genie3' else '70-100',
         'num_designs':5 if engine=='genie3' else 10,'seed':0,'gpu_device':'all','cpu_cores':4}
    if engine=='pxdesign':p.update(n_steps=400,dtype='bf16',use_fast_ln=True)
    elif engine=='genie3':p.update(direction_scale=.8 if task=='unconditional' else .1 if task=='motif' else 0.,compile_generation=False,cond_strategy='extended',length_step=1,motif_residues='',motif_example='',region_lengths=[])
    elif engine=='rfdiffusion':p.update(contig='',motif_residues='',designed_residues='',designed_chain='',region_lengths=[],timesteps=50,partial_steps=20,
        fixed_sequence_positions='',model_weights='auto',noise_scale_ca=1.,noise_scale_frame=1.,potential='',guide_scale=2.,guide_decay='quadratic',
        symmetry='C4',total_length=320,flank_before=50,flank_after=50,oligomer_contacts=True,deterministic=False,save_trajectory=False,
        scaffold_pdb=None,scaffold_ss_file=None,scaffold_adj_file=None,target_ss_file=None,target_adj_file=None,scaffold_preset='library',sampled_insertion=0,sampled_n=0,sampled_c=0)
    else:p.update(model_name='cc58' if task in {'unconditional','motif','partial'} else 'cc83',step_scale=1.0 if task=='partial' else 1.2,
        schurn=0 if task=='partial' else 200,crop_cond_start=0.,batch_size=1,partial_steps=100,motif_contig='',chain_lengths='128',
        motif_residues='',region_lengths=[],translation=[0.,0.,0.])
    if task=='unconditional':p['binder_length']='128'
    return p


def _integer(value, label, lower, upper):
    if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or int(value)!=value or not lower<=value<=upper:
        raise ValueError(f'{label} must be a whole number from {lower} to {upper}.')
    return int(value)


def _number(value, label, lower, upper):
    if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or not lower<=value<=upper:
        raise ValueError(f'{label} must be between {lower} and {upper}.')
    return float(value)


def length_range(value, *, zero=False):
    match=re.fullmatch(r'\s*(\d+)(?:\s*-\s*(\d+))?\s*',str(value))
    if not match or not (0 if zero else 1)<=int(match[1])<=int(match[2] or match[1])<=3000:
        raise ValueError('Use a length or increasing range, such as 80 or 70-100.')
    return [int(match[1]),int(match[2] or match[1])]


def residue_ranges(text, info, *, chains=None):
    """Expand only addressable residues; never silently skip absent selections."""
    result=[]
    for token in re.split(r'[,;\s]+',str(text or '').strip()):
        if not token:continue
        match=re.fullmatch(r'([A-Za-z])(\d+)(?:-(\d+))?',token)
        if not match:raise ValueError(f'Invalid selection {token!r}; use A30-40 or A59,A83,A91.')
        chain,start,end=match[1],int(match[2]),int(match[3] or match[2])
        if end<start or end-start>10000 or chains is not None and chain not in chains:
            raise ValueError('Residue selections must belong to the chosen chains and have increasing bounds.')
        labels=[f'{chain}{i}' for i in range(start,end+1)]
        if any(label not in info['residues'] or info['residues'][label]['kind']!='protein' or info['residues'][label]['insertion_code'] for label in labels):
            raise ValueError(f'{token} includes absent residues or insertion codes. Trim or remap the input first.')
        result.append(labels)
    return result


def _fixed_block(labels):
    return labels[0]+'-'+labels[-1][1:]


def _chain_blocks(info, chain):
    groups=[]
    for row in info['chains'][chain]['residues']:
        if row['insertion_code']:raise ValueError('Remap insertion-coded residues before generation.')
        label=row['label']
        if not groups or int(label[1:])!=int(groups[-1][-1][1:])+1:groups.append([])
        groups[-1].append(label)
    return [_fixed_block(g) for g in groups]


def _regions(value, count, default=(10,40)):
    value=value or [{'min':default[0],'max':default[1]} for _ in range(count)]
    if not isinstance(value,list) or len(value)!=count:raise ValueError(f'Specify lengths for all {count} designed regions.')
    result=[]
    for row in value:
        if not isinstance(row,dict) or set(row)-{'min','max'}:raise ValueError('Each region needs minimum and maximum lengths.')
        lo=_integer(row.get('min'),'Minimum region length',0,3000);hi=_integer(row.get('max'),'Maximum region length',lo,3000)
        result.append({'min':lo,'max':hi})
    return result


def _motif_contig(p, info, separator='/'):
    selections=residue_ranges(p.get('motif_residues'),info)
    if not selections:raise ValueError('Choose motif residues to preserve.')
    flat=[v for block in selections for v in block]
    if len(flat)!=len(set(flat)):raise ValueError('Select each motif residue once.')
    p['motif_labels']=flat
    p['region_lengths']=_regions(p.get('region_lengths'),len(selections)+1)
    pieces=[]
    for index,block in enumerate(selections):
        region=p['region_lengths'][index]
        if region['max']:pieces.append(f"{region['min']}-{region['max']}")
        pieces.append(_fixed_block(block))
    region=p['region_lengths'][-1]
    if region['max']:pieces.append(f"{region['min']}-{region['max']}")
    minimum=len(flat)+sum(r['min'] for r in p['region_lengths']);maximum=len(flat)+sum(r['max'] for r in p['region_lengths'])
    if not 1<=minimum<=maximum<=3000:raise ValueError('The scaffold length must be from 1 to 3000 residues.')
    p['scaffold_length_range']=[minimum,maximum]
    return separator.join(pieces)


def _validate_contig(text, info):
    text=str(text).strip().strip('[]').replace(',/0,','/0 ')
    if not text or not re.fullmatch(r'[A-Za-z0-9/ \-]+',text):raise ValueError('Use native RFdiffusion contigs, e.g. A1-150/0 70-100.')
    minimum=maximum=0;fixed=[]
    for token in re.split(r'[/\s]+',text):
        if not token or token=='0':continue
        if token[0].isalpha():
            blocks=residue_ranges(token,info)
            fixed.extend(blocks[0]);minimum+=len(blocks[0]);maximum+=len(blocks[0])
        else:
            lo,hi=length_range(token,zero=True);minimum+=lo;maximum+=hi
    if not 1<=minimum<=maximum<=10000:raise ValueError('Contig length must be from 1 to 10000 residues.')
    return text,minimum,maximum,fixed


def prepare_request(engine, parameters):
    if engine not in ENGINES:raise ValueError('Unknown generator.')
    task=parameters.get('task','binder')
    if task not in ENGINES[engine]['tasks']:raise ValueError('Choose a supported generation task.')
    p=default_parameters(engine,task);p.update(copy.deepcopy(parameters))
    p.update(generator_only=True,generator_engine=engine,contract='studio_native_generator_v1')
    p['num_designs']=_integer(p['num_designs'],'Number of designs',1,10000)
    p['seed']=_integer(p['seed'],'Random seed',0,2147483647)
    p['cpu_cores']=_integer(p['cpu_cores'],'CPU cores',1,256)
    p['length_range']=length_range(p['binder_length'])
    chains=p.get('target_chains') or []
    if isinstance(chains,str):chains=[c for c in re.split(r'[,;\s]+',chains) if c]
    if not isinstance(chains,list) or len(chains)!=len(set(chains)):raise ValueError('Select each target chain once.')
    p['target_chains']=chains
    info={'chains':{},'residues':{}}
    if p.get('target_pdb'):
        path=Path(p['target_pdb'])
        if not path.is_file():raise ValueError('The selected structure is unavailable.')
        info=inventory(path);p['target_pdb']=str(path)
        if any(c not in info['chains'] or info['chains'][c]['kind']!='protein' or not re.fullmatch(r'[A-Za-z]',c) for c in chains):
            raise ValueError('Choose protein chains with single-letter IDs. Use the Studio structure importer to remap long mmCIF chain IDs.')
    required=task not in {'unconditional','symmetry','custom','fold','multichain'}
    if required and not info['residues']:raise ValueError('Load a protein structure before submitting.')
    if task=='binder' and not chains:raise ValueError('Select target chains before submitting.')
    p['hotspots']=','.join(dict.fromkeys(label for block in residue_ranges(p.get('hotspots'),info,chains=chains or None) for label in block))
    if engine=='pxdesign':
        if p['length_range'][0]!=p['length_range'][1]:raise ValueError('PXDesign takes one binder length per job; choose a single length.')
        p['n_steps']=_integer(p['n_steps'],'Diffusion steps',1,2000)
        if p['dtype'] not in {'bf16','fp32'}:raise ValueError('Choose BF16 or FP32 precision.')
        if type(p['use_fast_ln']) is not bool:raise ValueError('Fast layer normalization must be true or false.')
    elif engine=='genie3':
        p['direction_scale']=_number(p['direction_scale'],'Direction scale',0,2)
        p['length_step']=_integer(p['length_step'],'Length increment',1,3000)
        if type(p['compile_generation']) is not bool:raise ValueError('Compile generation must be true or false.')
        if p['cond_strategy'] not in {'hotspot','extended','common'}:raise ValueError('Choose hotspot, extended, or common conditioning.')
        if task=='motif':
            if p['motif_example']:
                if not re.fullmatch(r'[A-Za-z0-9_-]+',p['motif_example']):raise ValueError('Choose a bundled motif example.')
                config_path=ASSETS/'genie3/motifbench/problems'/(p['motif_example']+'.json')
                if not config_path.is_file():raise ValueError('This motif example is unavailable.')
                p['motif_problem']=json.loads(config_path.read_text())
                p['motif_labels']=[]
                for line in Path(p['target_pdb']).read_text().splitlines():
                    if line.startswith('REMARK 999 INPUT'):
                        fields=line.split();p['motif_labels'].extend(v for b in residue_ranges(f'{fields[3]}{fields[4]}-{fields[5]}',info) for v in b)
                if not p['motif_labels']:raise ValueError('The bundled motif structure lost its segment headers.')
            else:p['motif_contig']=_motif_contig(p,info,',')
        # The target-conditioned native model expects an explicit interface.
        if task=='binder' and not p['hotspots']:raise ValueError('Select hotspot residues to condition binder generation.')
    elif engine=='rfdiffusion':
        for key,lower,upper in [('timesteps',1,200),('partial_steps',1,200),('total_length',1,10000),('flank_before',0,3000),('flank_after',0,3000),('sampled_insertion',0,100),('sampled_n',0,100),('sampled_c',0,100)]:
            p[key]=_integer(p[key],key.replace('_',' ').capitalize(),lower,upper)
        for key in ('noise_scale_ca','noise_scale_frame'):p[key]=_number(p[key],key,0,1)
        p['guide_scale']=_number(p['guide_scale'],'Guide scale',0,100)
        if p['guide_decay'] not in {'constant','linear','quadratic','cubic'}:raise ValueError('Choose a supported potential decay.')
        for key in ('oligomer_contacts','deterministic','save_trajectory'):
            if type(p[key]) is not bool:raise ValueError(f'{key} must be true or false.')
        if p['model_weights'] not in RF_MODELS:raise ValueError('Choose an installed RFdiffusion checkpoint.')
        if p['potential'] and not re.fullmatch(r'[A-Za-z0-9_.,:\-]+',p['potential']):raise ValueError('Enter a native potential definition, not a command or Hydra override.')
        if task=='binder':
            target='/0 '.join('/'.join(_chain_blocks(info,c)) for c in chains)
            p['contig']=target+'/0 '+p['binder_length']
        elif task=='motif':p['contig']=_motif_contig(p,info)
        elif task in {'redesign','partial'}:
            chain=p.get('designed_chain') or next((c for c in info['chains'] if c not in chains),next(iter(info['chains']),''))
            if chain not in info['chains']:raise ValueError('Choose the chain to redesign or diffuse.')
            p['designed_chain']=chain
            selected=residue_ranges(p.get('designed_residues'),info,chains=[chain])
            if not selected:raise ValueError('Select residues to redesign or diffuse.')
            selected_set={v for block in selected for v in block}
            rows=info['chains'][chain]['residues'];pieces=[];pending=[];groups=[]
            for row in rows:
                kind=row['label'] in selected_set
                if not groups or groups[-1][0]!=kind or row['number']!=int(groups[-1][1][-1][1:])+1:groups.append((kind,[]))
                groups[-1][1].append(row['label'])
            lengths=_regions(p.get('region_lengths'),sum(kind for kind,block in groups),default=(10,20)) if task=='redesign' else []
            li=0
            for designed,block in groups:
                if designed:
                    region=lengths[li] if task=='redesign' else {'min':len(block),'max':len(block)}
                    pieces.append(f"{region['min']}-{region['max']}");li+=1
                else:pieces.append(_fixed_block(block))
            p['region_lengths']=lengths
            # Partial diffusion must keep original chain order and total size.
            p['target_chains']=[c for c in info['chains'] if c!=chain and info['chains'][c]['kind']=='protein']
            p['motif_labels']=[r['label'] for c in info['chains'].values() if c['kind']=='protein' for r in c['residues'] if r['label'] not in selected_set]
            p['contig']='/0 '.join('/'.join(pieces) if c==chain else '/'.join(_chain_blocks(info,c)) for c in info['chains'] if info['chains'][c]['kind']=='protein')
        elif task=='symmetry':
            if p['symmetry'] not in SYMMETRIES or p['total_length']%SYMMETRIES[p['symmetry']]:raise ValueError('Total oligomer length must be divisible by the number of symmetric chains.')
            p['contig']=str(p['total_length'])
        elif task=='symmetric_motif':
            if p['symmetry'] not in SYMMETRIES:raise ValueError('Choose a supported symmetry.')
            blocks=residue_ranges(p['motif_residues'],info)
            if len(blocks)!=SYMMETRIES[p['symmetry']] or len({len(b) for b in blocks})!=1:
                raise ValueError('Provide one motif range of the same length for each symmetric copy.')
            p['motif_labels']=[label for block in blocks for label in block]
            p['contig']='/0 '.join(f"{p['flank_before']}/{_fixed_block(b)}/{p['flank_after']}" for b in blocks)
        elif task=='fold':
            if p['scaffold_preset'] not in {'library','tim_barrel','custom'}:raise ValueError('Choose a fold scaffold source.')
            if p['scaffold_preset']=='custom' and not p.get('scaffold_pdb') and not (p.get('scaffold_ss_file') and p.get('scaffold_adj_file')):
                raise ValueError('Load a binder scaffold structure or its secondary structure and adjacency files.')
            if info['residues'] and not chains:raise ValueError('Select target chains, or clear the target for monomer fold conditioning.')
        if task!='fold':
            p['contig'],minimum,maximum,fixed=_validate_contig(p['contig'],info)
            p['scaffold_length_range']=[minimum,maximum]
            if task=='partial':
                count=sum(len(c['residues']) for c in info['chains'].values() if c['kind']=='protein')
                if minimum!=maximum or maximum!=count:raise ValueError('Partial diffusion must preserve the input structure length.')
                if p['partial_steps']>p['timesteps']:raise ValueError('Partial diffusion steps cannot exceed total diffusion steps.')
            if task=='custom':p['motif_labels']=fixed
        fixed_positions=str(p.get('fixed_sequence_positions') or '').replace(' ','')
        if fixed_positions and (task!='partial' or not re.fullmatch(r'\d+(?:-\d+)?(?:,\d+(?:-\d+)?)*',fixed_positions)):
            raise ValueError('Fixed sequence positions use 0-based index ranges and are available for Partial Diffusion.')
        for token in fixed_positions.split(',') if fixed_positions else []:
            start,*end=map(int,token.split('-'))
            if not 0<=start<=(end[0] if end else start)<p['scaffold_length_range'][0]:
                raise ValueError('Fixed sequence positions must be increasing and within the 0-based input length.')
            ordered=[r['label'] for c in info['chains'].values() if c['kind']=='protein' for r in c['residues']]
            p['motif_labels']=list(dict.fromkeys([*p.get('motif_labels',[]),*ordered[start:(end[0] if end else start)+1]]))
        p['fixed_sequence_positions']=fixed_positions
        p['resolved_model']='InpaintSeq' if fixed_positions else 'Complex_Fold_base' if task=='fold' else 'Base_epoch8' if task=='symmetric_motif' else 'Complex_base' if task in {'binder','redesign'} else 'Base'
        if p['model_weights']!='auto':p['resolved_model']=p['model_weights']
    else:
        if p['model_name'] not in PROT_MODELS:raise ValueError('Choose a published Protpardelle checkpoint.')
        p['model_epoch']=PROT_MODELS[p['model_name']]
        p['step_scale']=_number(p['step_scale'],'Step scale',0.01,10)
        p['schurn']=_integer(p['schurn'],'Sampling churn',0,1000)
        p['batch_size']=_integer(p['batch_size'],'GPU batch size',1,128)
        p['crop_cond_start']=_number(p['crop_cond_start'],'Crop conditioning start',0,1)
        p['partial_steps']=_integer(p['partial_steps'],'Partial diffusion rewind steps',1,200)
        if not isinstance(p['translation'],list) or len(p['translation'])!=3:raise ValueError('Translation requires X, Y and Z coordinates.')
        p['translation']=[_number(v,'Translation',-10000,10000) for v in p['translation']]
        if task=='motif':
            p['motif_contig']=p.get('motif_contig') or _motif_contig(p,info,';')
            # Both the simple motif editor and native semicolon contigs are supported.
            labels=[]
            for token in p['motif_contig'].split(';'):
                if token and token[0].isalpha():labels.extend(v for b in residue_ranges(token,info) for v in b)
                elif token and token!='/':length_range(token,zero=True)
            if not labels:raise ValueError('The motif contig must preserve at least one residue.')
            p['motif_labels']=labels
        if task in {'motif','multichain','unconditional'}:
            p['chain_length_ranges']=[length_range(v) for v in str(p['chain_lengths']).split(';')]
            if task=='motif':
                groups=p['motif_contig'].split(';/;')
                if len(groups)!=len(p['chain_length_ranges']):raise ValueError('Set one total length per scaffold chain.')
                for group,total in zip(groups,p['chain_length_ranges']):
                    _,minimum,maximum,_=_validate_contig(group.replace(';','/'),info)
                    if maximum<total[0] or minimum>total[1]:raise ValueError('Scaffold chain length does not fit its motif and designed regions.')
            if task=='unconditional' and len(p['chain_length_ranges'])!=1:raise ValueError('Use Multichain Generation for more than one chain.')
        if task=='multichain':
            if not p['motif_contig']:raise ValueError('Provide the multichain motif contig, with ;/; between chains.')
            if len(p['motif_contig'].split(';/;'))!=len(p['chain_length_ranges']):raise ValueError('Set one total length per generated chain.')
            p['motif_labels']=[];fixed_chains=[]
            for group,total in zip(p['motif_contig'].split(';/;'),p['chain_length_ranges']):
                _,minimum,maximum,labels=_validate_contig(group.replace(';','/'),info)
                if maximum<total[0] or minimum>total[1]:raise ValueError('Generated chain length does not fit its motif and designed regions.')
                p['motif_labels'].extend(labels)
                for chain,row in info['chains'].items():
                    if labels and set(labels)=={r['label'] for r in row['residues']} and minimum==maximum==len(labels):fixed_chains.append(chain)
            p['target_chains']=list(dict.fromkeys([*chains,*fixed_chains]))
        if task in {'motif','multichain'}:
            # The native parser requires explicit start-end ranges, including
            # fixed single lengths and single motif residues.
            groups=[]
            for group in p['motif_contig'].split(';/;'):
                tokens=[]
                for token in group.split(';'):
                    if token[0].isalpha():
                        block=residue_ranges(token,info)[0];tokens.append(_fixed_block(block))
                    else:
                        lo,hi=length_range(token,zero=True);tokens.append(f'{lo}-{hi}')
                groups.append(';'.join(tokens))
            p['motif_contig']=';/;'.join(groups)
        allatom=p['model_name'] in {'cc89','cc91','cc94'}
        if p['model_name']=='cc89' and task!='partial':raise ValueError('The cc89 model requires an input sequence and is offered for partial diffusion.')
        if task=='unconditional':sampling='sampling_unconditional_allatom_s1' if allatom else 'sampling_unconditional'
        elif task=='partial':sampling='sampling_partial_diffusion_allatom' if allatom else 'sampling_partial_diffusion'
        else:sampling='sampling_sidechain_conditional_allatom_s1' if allatom else 'sampling_sidechain_conditional'
        p['sampling_config']=sampling
    for key in ('scaffold_pdb','scaffold_ss_file','scaffold_adj_file','target_ss_file','target_adj_file'):
        if p.get(key):
            path=Path(p[key])
            if not path.is_file():raise ValueError(f'{key.replace("_"," ")} is unavailable.')
            p[key]=str(path)
    p['fixed_identity_labels']=[label for label in p.get('motif_labels',[]) if info['residues'][label]['name']!='UNK']
    p['fixed_motif_required']=bool(p['fixed_identity_labels'])
    return p


def readiness(engine, parameters, root=None):
    root=Path(root or reference_root());task=parameters.get('task','binder')
    if engine=='rfdiffusion':
        model=parameters.get('resolved_model') or ('Complex_Fold_base' if task=='fold' else 'Base_epoch8' if task=='symmetric_motif' else 'Complex_base' if task in {'binder','redesign'} else 'Base')
        if parameters.get('model_weights','auto')!='auto':model=parameters['model_weights']
        paths=[root/'rfdiffusion_models'/f'{model}_ckpt.pt']
        if task=='fold' and parameters.get('scaffold_preset','library')=='library':paths.append(root/'rfdiffusion_models/ppi_scaffolds')
    elif engine=='pxdesign':paths=[root/'pxdesign/release_data/checkpoint',root/'pxdesign/release_data/ccd_cache']
    elif engine=='genie3':paths=[root/'genie3/pretrained/v1/checkpoints/step=600000.ckpt',root/'genie3/pretrained/v1/config.yaml']
    else:
        model=parameters.get('model_name','cc83');paths=[root/'protpardelle-1c/model_params/weights'/f'{model}_epoch{PROT_MODELS.get(model, "")}.pth']
        if model in {'cc89','cc91','cc94'}:paths.append(root/'protpardelle-1c/model_params/weights/cc58_epoch97_minimpnn.pth')
    missing=[str(p) for p in paths if not p.exists() or p.is_file() and not p.stat().st_size or p.is_dir() and not any(p.iterdir())]
    return {'ready':not missing,'missing':missing,'required':[str(p) for p in paths]}


def examples(engine, task=None):
    """Bundled coordinates are offline, immutable copies with source attribution."""
    out=[]
    def add(key,label,which,**parameters):
        p=default_parameters(engine,which);p.update(parameters)
        out.append({'id':key,'label':label,'task':which,'parameters':p})
    if engine=='pxdesign':
        add('tamarind','1CRN · Tamarind example','binder',target_pdb=str(ASSETS/'1CRN.pdb'),target_chains=['A'],hotspots='A12-33',binder_length='10')
        add('pdl1','PD-L1 · GitHub quick start','binder',target_pdb=str(ASSETS/'pxdesign/5o45.cif'),target_chains=['A'],hotspots='A40,A99,A107',binder_length='80')
    elif engine=='genie3':
        add('tamarind','1IEP · Tamarind example','binder',target_pdb=str(ASSETS/'1IEP.pdb'),target_chains=['A'],hotspots='A261,A263,A264')
        for file in sorted((ASSETS/'genie3/binderbench/problems').glob('*.json')):
            data=json.loads(file.read_text());name=file.stem
            target=ASSETS/'genie3/binderbench/targets/pdb'/f'{name}.pdb'
            if target.is_file():add(name,name.replace('_',' ',1)+' · BinderBench','binder',target_pdb=str(target),target_chains=list(dict.fromkeys(s[0] for s in data['target_chain_and_residues'])),hotspots=','.join(data['target_interface_residues']['hotspot']),binder_length=f"{data['binder_min_length']}-{data['binder_max_length']}")
        add('unconditional','128 residues · GitHub example','unconditional',binder_length='128')
        for file in sorted((ASSETS/'genie3/motifbench/problems').glob('*.json')):
            data=json.loads(file.read_text());pdb=ASSETS/'genie3/motifbench/motifs'/f'{file.stem}.pdb'
            # Preserve the official segment configuration, including multi-motif ordering.
            if pdb.is_file():add(file.stem,file.stem+' · MotifBench','motif',target_pdb=str(pdb),motif_example=file.stem,binder_length=str(data.get('minimum_total_length',125)))
    elif engine=='rfdiffusion':
        base=ASSETS/'rfdiffusion/input_pdbs'
        add('tamarind','5TPN · Tamarind motif example','motif',target_pdb=str(base/'5TPN.pdb'),motif_residues='A30-40,A60-70',region_lengths=[{'min':0,'max':0},{'min':19,'max':19},{'min':0,'max':0}])
        add('rsv','RSV-F site 5 · GitHub','motif',target_pdb=str(base/'5TPN.pdb'),motif_residues='A163-181',region_lengths=[{'min':10,'max':40},{'min':10,'max':40}])
        add('insulin','Insulin receptor · GitHub','binder',target_pdb=str(base/'insulin_target.pdb'),target_chains=['A'],hotspots='A59,A83,A91',noise_scale_ca=0.,noise_scale_frame=0.)
        add('partial','2KL8 · GitHub partial diffusion','partial',target_pdb=str(base/'2KL8.pdb'),designed_chain='A',designed_residues='A1-79',partial_steps=10)
        add('redesign','2KL8 · Redesign a region','redesign',target_pdb=str(base/'2KL8.pdb'),designed_chain='A',designed_residues='A20-40',region_lengths=[{'min':21,'max':21}])
        add('custom','100 residues · GitHub unconditional','custom',contig='100',target_pdb=None)
        add('cyclic','C6 · GitHub oligomer','symmetry',symmetry='C6',total_length=480)
        add('dihedral','D2 · GitHub oligomer','symmetry',symmetry='D2',total_length=400)
        add('tetrahedral','Tetrahedral · GitHub oligomer','symmetry',symmetry='T',total_length=1200)
        add('nickel','Nickel site · GitHub symmetric motif','symmetric_motif',target_pdb=str(base/'nickel_symmetric_motif.pdb'),motif_residues='A2-4,A7-9,A12-14,A17-19')
        add('timbarrel','TIM barrel · GitHub fold conditioning','fold',scaffold_preset='tim_barrel',target_pdb=None,noise_scale_ca=.5,noise_scale_frame=.5,sampled_insertion=5,sampled_n=5,sampled_c=5)
        add('insulin_fold','Insulin receptor · GitHub fold conditioning','fold',target_pdb=str(base/'insulin_target.pdb'),target_chains=['A'],hotspots='A59,A83,A91',scaffold_preset='library',target_ss_file=str(ASSETS/'rfdiffusion/target_folds/insulin_target_ss.pt'),target_adj_file=str(ASSETS/'rfdiffusion/target_folds/insulin_target_adj.pt'),noise_scale_ca=0.,noise_scale_frame=0.)
    else:
        base=ASSETS/'protpardelle/motifs/nanobody'
        pdb=base/'7eow_B_atom_rot_128_A.pdb'
        add('binder','Nanobody target · GitHub multichain example','binder',target_pdb=str(pdb),target_chains=['A'],hotspots='A19,A76,A80,A82',binder_length='120')
        add('unconditional','128 residues · GitHub unconditional','unconditional')
        add('partial','7EOW · GitHub partial diffusion','partial',target_pdb=str(base/'7eow_B_atom.pdb'),target_chains=[],binder_length='128')
        add('cdr3','7EOW CDR3 · GitHub motif scaffolding','motif',target_pdb=str(base/'7eow_CDR3_atom_rot_128.pdb'),motif_contig='50-150;A1-17;10-150',chain_lengths='128')
        add('multichain','Target + 120 residues · GitHub','multichain',target_pdb=str(pdb),target_chains=['A'],hotspots='A19,A76,A80,A82',motif_contig='A1-128;/;120-120',chain_lengths='128;120')
        import yaml
        for group,config in [('motif_bench','03_motifbench'),('bindcraft','04_bindcraft')]:
            data=yaml.safe_load((ASSETS/f'protpardelle/sampling/{config}.yaml').read_text())
            for name,contig,lengths,hotspots in zip(data['motifs'],data['motif_contigs'],data['total_lengths'],data['hotspots']):
                file=ASSETS/'protpardelle/motifs'/group/(str(name)+'.pdb')
                if file.is_file():
                    which='motif' if group=='motif_bench' else 'multichain'
                    add(group+'-'+str(name),str(name)+' · '+('MotifBench' if group=='motif_bench' else 'BindCraft targets'),which,target_pdb=str(file),motif_contig=contig,chain_lengths=';'.join(str(a) if a==b else f'{a}-{b}' for a,b in lengths),hotspots=hotspots or '')
    return [e for e in out if task is None or e['task']==task]


def canvas_fields(engine):
    """Persist the same public form values, without per-job scheduler controls."""
    result=[]
    choices={'task':list(ENGINES[engine]['tasks']),'model_weights':RF_MODELS,'model_name':list(PROT_MODELS),
             'dtype':['bf16','fp32'],'cond_strategy':['extended','hotspot','common'],
             'guide_decay':['constant','linear','quadratic','cubic'],'symmetry':list(SYMMETRIES),'scaffold_preset':['library','custom','tim_barrel']}
    for name,value in default_parameters(engine).items():
        if name in {'campaign_name','gpu_device','cpu_cores'}:continue
        if name=='target_chains':value=''
        if value is None:value=''
        kind='boolean' if isinstance(value,bool) else 'number' if isinstance(value,(int,float)) else 'json' if isinstance(value,list) else 'text'
        if name in choices:kind='select'
        field={'name':name,'label':name.replace('_',' ').capitalize(),'default':value,'type':kind,'section':'basic' if name in {'num_designs','binder_length','task'} else 'optional'}
        if name in choices:field['options']=choices[name]
        if kind=='number':
            field['min']=0
            field['max']=2147483647 if name=='seed' else 10000
        result.append(field)
    return result
