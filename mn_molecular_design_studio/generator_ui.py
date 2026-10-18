"""Tamarind-inspired generator pages, shared unchanged with their canvas forms."""
from __future__ import annotations

import copy
import hashlib
import inspect
from html import escape
from pathlib import Path
import re
import shutil

import streamlit as st

from mn_molecular_design_studio import generator_input as inputs
from mn_molecular_design_studio.auth import user_input_root
from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer import StructureVisualization, molstar_custom_component
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json
from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path
from mn_molecular_design_studio.target_editor import TargetSource, _hash, _target_widget_keys, render_target_editor, render_target_source

STYLE='''<style>
 .st-key-generator-surface{--generator-green:#27674e;color:#193e30}
 .st-key-generator-surface .gen-heading{text-align:center;padding:.5rem 0 1.25rem}
 .gen-heading h1{color:#173f31;font-size:2.15rem!important;letter-spacing:-.035em;margin:.25rem 0 .5rem!important}
 .gen-heading p{color:#75877c;font-size:.93rem;margin:.5rem auto;max-width:660px;line-height:1.6}
 .gen-eyebrow,.gen-section{font-size:.7rem;font-weight:600;letter-spacing:.11em;color:#658572;text-transform:uppercase}
 .gen-section{margin:.9rem 0 .7rem}.gen-version{font-size:.73rem;color:#82958a}
 .gen-mode{display:inline-block;border:1px solid #bdd4c6;color:#27674e;background:#edf6f0;border-radius:6px;padding:.35rem .8rem;font-size:.78rem}
 :is(.st-key-generator-surface,.st-key-pipeline-node-form) [data-testid=stButton] button{border-radius:6px!important;border-color:#b7d0c0!important;color:#27674e;box-shadow:none!important}
 :is(.st-key-generator-surface,.st-key-pipeline-node-form) [data-testid=stButton] button:hover{border-color:#69977c!important;background:#f0f7f2!important;color:#174d32!important}
 :is(.st-key-generator-surface,.st-key-pipeline-node-form) button[kind=primary]{background:#27674e!important;border-color:#27674e!important;color:white!important}
 :is(.st-key-generator-surface,.st-key-pipeline-node-form) button[kind=primary]:disabled{background:#a4c0ae!important;color:white!important}
 :is(.st-key-generator-surface,.st-key-pipeline-node-form) [data-variant=pills]{border:1px solid #bfd5c7!important;border-radius:6px!important;background:#fbfdfb;color:#38634c!important;padding:.45rem .85rem!important;min-height:2.5rem!important}
 :is(.st-key-generator-surface,.st-key-pipeline-node-form) [data-variant=pills]:is([aria-checked=true],[aria-pressed=true]){background:#27674e!important;color:#fff!important;border-color:#27674e!important}
 :is(.st-key-generator-surface,.st-key-pipeline-node-form) [data-baseweb=input],:is(.st-key-generator-surface,.st-key-pipeline-node-form) [data-baseweb=textarea],:is(.st-key-generator-surface,.st-key-pipeline-node-form) [data-baseweb=select]>div{border-color:#d5e2d8!important;background:#f3f7f4!important;border-radius:6px!important}
 :is(.st-key-generator-surface,.st-key-pipeline-node-form) :is(input,textarea){background:#f3f7f4!important;color:#243e30!important}
 :is(.st-key-generator-surface,.st-key-pipeline-node-form) [data-baseweb=tag]{background:#e0eee5!important;color:#27674e!important}
 :is(.st-key-generator-surface,.st-key-pipeline-node-form) [data-testid=stExpander]{border-color:#d6e3da!important;border-radius:7px!important;background:#fcfefc}
 :is(.st-key-generator-surface,.st-key-pipeline-node-form) [data-testid=stFileUploaderDropzone]{border:1px dashed #bad1c1!important;background:#f5f9f6!important;border-radius:8px}
 [class*=st-key-generator][class*=-sampling-card]{padding:1rem 1.2rem!important;background:#f3f8f5;border:1px solid #d8e6dc;border-radius:10px;margin:.7rem 0}
 [class*=st-key-generator][class*=-input-card]{padding:1rem 1.2rem!important;border:1px solid #dbe6de;border-radius:10px;background:white;margin:.6rem 0}
 :is(.st-key-generator-surface,.st-key-pipeline-node-form) [data-testid=stCaptionContainer]{color:#809186}
 .st-key-generator-surface [data-testid=stLinkButton] a{border:none;color:#597766;font-size:.76rem;padding:.15rem .3rem;min-height:1.8rem}
 .st-key-generator-surface [data-testid=stButtonGroup]{gap:.45rem}
 </style>'''


def public_parameters(engine,p):
    from mn_molecular_design_studio.engine.workflow_engine.workflows import generators
    function=getattr(generators,'run_'+engine+'_generation')
    return {k:v for k,v in p.items() if k in inspect.signature(function).parameters}


def _private_example(engine, example, owner_id):
    p=copy.deepcopy(example['parameters'])
    for key in ('target_pdb','scaffold_pdb','target_ss_file','target_adj_file','scaffold_ss_file','scaffold_adj_file'):
        if p.get(key):
            source=Path(p[key]);content=source.read_bytes()
            destination=user_input_root(owner_id)/'generator-examples'/engine/hashlib.sha256(content).hexdigest()[:12]/source.name
            destination.parent.mkdir(parents=True,exist_ok=True)
            if not destination.exists():destination.write_bytes(content)
            p[key]=str(destination)
    return p


def _load_example(engine,example,owner_id,prefix):
    p=_private_example(engine,example,owner_id)
    # Widget callbacks run before rerendering, so task changes and resets keep
    # other forms, project choices and saved pipeline versions intact.
    for key in list(st.session_state):
        if key.startswith(prefix+'-fields-'):del st.session_state[key]
    st.session_state[prefix+'-preset']=p
    st.session_state[prefix+'-task']=p['task']


def _compact(text,info):
    groups=[]
    for label in dict.fromkeys(v for block in inputs.residue_ranges(text,info) for v in block):
        if not groups or groups[-1][-1][0]!=label[0] or int(groups[-1][-1][1:])+1!=int(label[1:]):groups.append([])
        groups[-1].append(label)
    return ','.join(block[0] if len(block)==1 else block[0]+'-'+block[-1][1:] for block in groups)


def _picker(source, p, prefix, field='hotspots', title='Hotspots', allow_trimming=True):
    info=inputs.inventory(source.path)
    chosen=p.get('target_chains') or []
    if isinstance(chosen,str):chosen=[v for v in re.split(r'[,;\s]+',chosen) if v]
    if field in {'motif_residues','designed_residues'}:
        try:chosen=list(dict.fromkeys(block[0][0] for block in inputs.residue_ranges(p.get(field,''),info)))
        except ValueError:chosen=[]
        if not chosen and p.get('designed_chain'):chosen=[p['designed_chain']]
    if chosen:st.session_state.setdefault(prefix+'-chains-'+_hash(str(source.path)),chosen)
    st.session_state.setdefault(_target_widget_keys(prefix,source.path)['hotspots'],p.get(field,'') or '')
    selected=render_target_editor(source,key_prefix=prefix,label='Target' if field=='hotspots' else 'Motif' if field=='motif_residues' else 'Binder',
        residue_selection_title=title,residue_selection_input_label='Chosen residues',allow_trimming=allow_trimming,viewer_height=430,show_structure_tools=False)
    if field=='hotspots':
        p.update(target_pdb=str(selected.design_path or source.path),target_chains=selected.chains,hotspots=selected.hotspots)
    else:
        p[field]=_compact(selected.hotspots,info) if selected.hotspots else ''
        p['target_pdb']=str(selected.design_path or source.path) if allow_trimming else str(source.path)
        if field=='designed_residues':
            if len(selected.chains)!=1:raise ValueError('Choose exactly one chain to redesign or diffuse.')
            p['designed_chain']=selected.chains[0]
    return selected


def _regions(p, count, prefix, *, label='Lengths to replace each designed region', default=(10,40)):
    st.markdown('##### '+label)
    st.caption('Choose a fixed length or a range to sample for each designed region.')
    existing=p.get('region_lengths') or []
    values=[]
    for index in range(count):
        value=existing[index] if len(existing)==count else {'min':default[0],'max':default[1]}
        with st.container(border=True):
            a,b,c=st.columns([1.3,1,1],vertical_alignment='bottom')
            a.caption(('Before motif' if index==0 else 'After motif' if index==count-1 else 'Between motifs') if label.startswith('Lengths') else f'Region {index+1}')
            lo=int(b.number_input('Min',0,3000,int(value['min']),key=f'{prefix}-region-{count}-{index}-min'))
            hi=int(c.number_input('Max',0,3000,int(value['max']),key=f'{prefix}-region-{count}-{index}-max'))
            values.append({'min':lo,'max':hi})
    p['region_lengths']=values


def _preview(source,prefix):
    with st.expander('Structure preview',expanded=False):
        molstar_custom_component([StructureVisualization(pdb=source.path.read_text(),color='chain-id',representation_type='cartoon')],
            key=prefix+'-preview-'+_hash(str(source.path)),height=400,show_structure_tools=False)


def _file(label,current,prefix,owner_id):
    uploaded=st.file_uploader(label,type=['pt'],key=prefix)
    if uploaded:
        content=uploaded.getvalue();path=user_input_root(owner_id)/'generator-folds'/hashlib.sha256(content).hexdigest()[:12]/Path(uploaded.name).name
        path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(content)
        return str(path)
    if current:st.caption('Loaded: '+Path(current).name)
    return current


def render_inputs(engine, *,owner_id,key_prefix,initial=None,connected_target=None,include_jobs=True):
    st.html(STYLE)
    spec=inputs.ENGINES[engine];initial=copy.deepcopy(initial or {})
    st.html('<div class="gen-section">CHOOSE TASK</div>')
    task=st.pills('Design task',list(spec['tasks']),format_func=spec['tasks'].get,selection_mode='single',
        default=initial.get('task','binder'),key=key_prefix+'-task',label_visibility='collapsed') or 'binder'
    p=inputs.default_parameters(engine,task)
    if connected_target:p.update({k:v for k,v in connected_target.items() if k in p and v not in ('',None,[])})
    if initial.get('task','binder')==task:p.update({k:v for k,v in initial.items() if v not in ('',None,[])})
    preset=st.session_state.get(key_prefix+'-preset',{})
    if preset.get('task')==task:p.update(preset)
    p['task']=task
    if isinstance(p.get('target_chains'),str):p['target_chains']=[c for c in re.split(r'[,;\s]+',p['target_chains']) if c]
    prefix=key_prefix+'-fields-'+task;error=None
    descriptions={'binder':'Given a target structure, design a binder to bind to your target.',
        'motif':'Scaffold selected functional residues with newly designed regions.',
        'redesign':'Replace selected regions of an existing binder while preserving its remaining structure.',
        'partial':'Diffuse selected regions around an existing fold while preserving its length.',
        'fold':'Condition generation on the secondary structure and adjacency of a binder fold.',
        'symmetry':'Design a symmetric cyclic, dihedral, or tetrahedral oligomer.',
        'symmetric_motif':'Scaffold one spatially arranged motif for every symmetric copy.',
        'custom':'Write contigs directly using native RFdiffusion notation.',
        'unconditional':'Generate proteins without an input structure.',
        'multichain':'Generate multiple chains with optional fixed chains or functional motifs.'}
    st.caption(descriptions[task])
    examples=inputs.examples(engine,task)
    if examples:
        a,b=st.columns([2.2,1],vertical_alignment='bottom')
        example=a.selectbox('Examples',examples,format_func=lambda e:e['label'],key=prefix+'-example-choice')
        b.button('Load example inputs',icon=':material/science:',key=prefix+'-load-example',use_container_width=True,
            on_click=_load_example,args=(engine,example,owner_id,key_prefix))
    needs_structure=task not in {'unconditional','symmetry'}
    source=TargetSource()
    with st.container(key=key_prefix+'-input-card'):
        if needs_structure:
            format_value=st.pills('Structure format',['PDB','CIF'],selection_mode='single',default='CIF' if Path(p.get('target_pdb') or '').suffix.lower() in {'.cif','.mmcif'} else 'PDB',key=prefix+'-format',label_visibility='collapsed') or 'PDB'
            source=render_target_source(key_prefix=prefix,label='Target' if task in {'binder','fold'} else 'Input',example_code='1CRN' if engine=='pxdesign' else '1IEP' if engine=='genie3' else '5TPN',
                owner_id=owner_id,allow_cif=True,allow_non_protein=task=='motif',upload_file_types=['pdb','ent'] if format_value=='PDB' else ['cif','mmcif'],
                upload_help='Upload coordinates for the target or fixed motif. mmCIF is converted to a PDB-compatible copy when required.')
            if not source.path and p.get('target_pdb'):
                resolved=resolve_stored_path(p['target_pdb'],must_exist=True)
                if resolved:source=TargetSource(path=resolved,name=resolved.stem,owner_id=owner_id)
            p['target_pdb']=str(source.path) if source.path else None
        else:p.update(target_pdb=None,target_chains=[],hotspots='')
        if source.path:
            try:
                if task in {'binder','fold'}:
                    _picker(source,p,prefix)
                elif engine=='protpardelle' and task=='partial':
                    st.success(source.name+' loaded')
                    st.caption('Partial diffusion uses the entire input structure and preserves its chain lengths.')
                    _preview(source,prefix)
                elif task in {'motif','redesign','partial'}:
                    official=engine=='genie3' and bool(p.get('motif_example')) and st.checkbox('Use official MotifBench layout',value=True,key=prefix+'-official-layout')
                    if engine=='genie3' and not official:p['motif_example']=''
                    native=engine=='protpardelle' and task=='motif' and st.checkbox('Use native motif contig',value=bool(p.get('motif_contig')),key=prefix+'-native-contig')
                    if official:
                        st.success(source.name+' loaded')
                        st.caption('Using the official MotifBench layout: '+p['motif_example'])
                        _preview(source,prefix)
                    elif native:
                        st.success(source.name+' loaded')
                        p['motif_contig']=st.text_input('Motif contig',str(p.get('motif_contig') or ''),key=prefix+'-motif-contig',help='Use semicolons within a chain and ;/; between chains.')
                        _preview(source,prefix)
                    else:
                        field='motif_residues' if task=='motif' else 'designed_residues'
                        _picker(source,p,prefix,field,'Motif residues' if task=='motif' else 'Designed residues' if task=='redesign' else 'Diffused residues',allow_trimming=task=='motif')
                        if task=='motif':
                            p['motif_contig']=''
                            count=len(inputs.residue_ranges(p[field],inputs.inventory(source.path)))+1
                            _regions(p,count,prefix)
                        elif task=='redesign':
                            count=max(1,len(inputs.residue_ranges(p[field],inputs.inventory(source.path))))
                            _regions(p,count,prefix,label='Replacement lengths for designed regions',default=(10,20))
                else:
                    st.success(source.name+' loaded')
                    info=inputs.inventory(source.path)
                    if task=='multichain':
                        chain_options=[c for c,row in info['chains'].items() if row['kind']=='protein']
                        p['target_chains']=st.multiselect('Fixed target chains (optional)',chain_options,default=[c for c in p['target_chains'] if c in chain_options],key=prefix+'-target-chains')
                    _preview(source,prefix)
            except (ValueError,KeyError,OSError) as exc:error=str(exc)
        if task=='custom':p['contig']=st.text_input('Contigs',str(p.get('contig') or ''),key=prefix+'-contig',placeholder='A1-150/0 70-100',help='Fixed PDB residue ranges and designed lengths, separated by /; use /0 and a space for chain breaks.')
        if task=='symmetric_motif':
            p['motif_residues']=st.text_input('Motif residues',p['motif_residues'],key=prefix+'-motif-residues',placeholder='A2-4,A7-9,A12-14,A17-19')
            a,b=st.columns(2)
            p['flank_before']=int(a.number_input('Scaffold length before motif',0,3000,int(p['flank_before']),key=prefix+'-flank-before'))
            p['flank_after']=int(b.number_input('Scaffold length after motif',0,3000,int(p['flank_after']),key=prefix+'-flank-after'))
        if task in {'symmetry','symmetric_motif'}:
            p['symmetry']=st.selectbox('Symmetry type',list(inputs.SYMMETRIES),index=list(inputs.SYMMETRIES).index(p['symmetry']),format_func=lambda s:'Tetrahedral (12 chains)' if s=='T' else s,key=prefix+'-symmetry')
            if task=='symmetry':p['total_length']=int(st.number_input('Total oligomer length',1,10000,int(p['total_length']),key=prefix+'-oligomer-length',help='The sum of residues across every chain; it must be divisible by the symmetry order.'))
        if task=='multichain':p['motif_contig']=st.text_input('Multichain contig',p.get('motif_contig',''),key=prefix+'-multichain-contig',placeholder='A1-128;/;120-120')
        if task in {'binder','unconditional'}:
            if engine=='pxdesign':p['binder_length']=str(int(st.number_input('Binder length',1,3000,inputs.length_range(p['binder_length'])[0],key=prefix+'-length')))
            elif engine=='protpardelle' and task=='unconditional':p['chain_lengths']=st.text_input('Protein length',p['chain_lengths'],key=prefix+'-chain-lengths',help='One length or a minimum–maximum range.')
            else:
                lo,hi=inputs.length_range(p['binder_length']);a,b=st.columns(2)
                minimum=int(a.number_input('Min binder length' if task=='binder' else 'Minimum protein length',1,3000,lo,key=prefix+'-length-min'))
                maximum=int(b.number_input('Max binder length' if task=='binder' else 'Maximum protein length',1,3000,hi,key=prefix+'-length-max'))
                p['binder_length']=str(minimum) if minimum==maximum else f'{minimum}-{maximum}'
        elif engine=='protpardelle' and task in {'motif','multichain'}:
            p['chain_lengths']=st.text_input('Total lengths for scaffold chains',p['chain_lengths'],key=prefix+'-chain-lengths',help='One length or range per chain, separated by semicolons, e.g. 128;100-120.')
        if engine=='rfdiffusion' and task=='fold':
            choices=['library','custom','tim_barrel']
            p['scaffold_preset']=st.pills('Binder fold',choices,format_func={'library':'Scaffold library','custom':'Upload / fetch binder','tim_barrel':'TIM barrel'}.get,selection_mode='single',
                default=p['scaffold_preset'],key=prefix+'-scaffold-preset') or 'library'
            if p['scaffold_preset']=='custom':
                scaffold=render_target_source(key_prefix=prefix+'-scaffold',label='Binder',owner_id=owner_id,allow_cif=True,upload_help='A protein structure defining the binder fold.')
                p['scaffold_pdb']=str(scaffold.path) if scaffold.path else p.get('scaffold_pdb')
    with st.container(key=key_prefix+'-sampling-card'):
        st.html('<div class="gen-section">DESIGN SAMPLING</div>')
        p['num_designs']=int(st.number_input('Number of designs'+(' per length' if engine=='genie3' and task=='unconditional' else ''),1,10000,int(p['num_designs']),key=prefix+'-designs'))
        st.caption('Generation only · sequence redesign and structure prediction are selected on the pipeline canvas.')
    settings_open=st.session_state.get(key_prefix+'-show-settings',False)
    with st.expander('View optional settings',expanded=settings_open):
        a,b=st.columns(2)
        with a:
            p['seed']=int(st.number_input('Random seed (0 = random)',0,2147483647,int(p['seed']),key=prefix+'-seed'))
            if engine=='pxdesign':
                p['n_steps']=int(st.number_input('Diffusion steps',1,2000,int(p['n_steps']),key=prefix+'-steps'))
                p['dtype']=st.selectbox('Precision',['bf16','fp32'],index=['bf16','fp32'].index(p['dtype']),key=prefix+'-dtype')
            elif engine=='genie3':
                p['direction_scale']=float(st.number_input('Direction scale',0.,2.,float(p['direction_scale']),.1,key=prefix+'-direction',help='Quality–diversity trade-off. Authors suggest 0.0 for binders, 0.1 for motifs, and 0.8 for short unconditional monomers.'))
                if task=='unconditional':p['length_step']=int(st.number_input('Length increment',1,3000,int(p['length_step']),key=prefix+'-length-step'))
            elif engine=='protpardelle':
                labels={'cc58':'cc58 · single-chain backbone / motif','cc78':'cc78 · homodimer','cc83':'cc83 · binder / multichain','cc89':'cc89 · all-atom partial diffusion','cc91':'cc91 · single-chain all-atom','cc94':'cc94 · multichain all-atom','cc95':'cc95 · multichain backbone'}
                p['model_name']=st.selectbox('Model',list(labels),index=list(labels).index(p['model_name']),format_func=labels.get,key=prefix+'-model')
                p['step_scale']=float(st.number_input('Step scale',.01,10.,float(p['step_scale']),.1,key=prefix+'-step-scale'))
                p['schurn']=int(st.number_input('Sampling churn',0,1000,int(p['schurn']),key=prefix+'-churn'))
            else:
                p['timesteps']=int(st.number_input('Diffusion steps',1,200,int(p['timesteps']),key=prefix+'-steps'))
                p['noise_scale_ca']=float(st.number_input('Noise scale (C alpha)',0.,1.,float(p['noise_scale_ca']),.1,key=prefix+'-noise-ca'))
                p['noise_scale_frame']=float(st.number_input('Noise scale (frame)',0.,1.,float(p['noise_scale_frame']),.1,key=prefix+'-noise-frame'))
                p['model_weights']=st.selectbox('Model',inputs.RF_MODELS,index=inputs.RF_MODELS.index(p['model_weights']),format_func=lambda m:'Automatic (task-specific)' if m=='auto' else m,key=prefix+'-model')
        with b:
            if engine=='pxdesign':
                p['use_fast_ln']=st.checkbox('Fast layer normalization',bool(p['use_fast_ln']),key=prefix+'-fast-ln')
                st.caption('MSAs and Protenix evaluation are used by downstream folding steps. The PXDesign generation stage does not need an MSA.')
            elif engine=='genie3':
                p['compile_generation']=st.checkbox('Compile generation model',bool(p['compile_generation']),key=prefix+'-compile')
                if task=='binder':p['cond_strategy']=st.selectbox('Hotspot conditioning',['extended','hotspot','common'],index=['extended','hotspot','common'].index(p['cond_strategy']),key=prefix+'-conditioning')
            elif engine=='protpardelle':
                p['crop_cond_start']=float(st.number_input('Crop conditioning start',0.,1.,float(p['crop_cond_start']),.1,key=prefix+'-crop-start'))
                p['batch_size']=int(st.number_input('GPU batch size',1,128,int(p['batch_size']),key=prefix+'-batch'))
                if task=='partial':p['partial_steps']=int(st.number_input('Rewind diffusion steps',1,200,int(p['partial_steps']),key=prefix+'-rewind'))
                with st.expander('Initial translation'):
                    p['translation']=[float(st.number_input(axis,-10000.,10000.,float(value),key=prefix+'-translation-'+axis)) for axis,value in zip('XYZ',p['translation'])]
            else:
                if task=='partial':
                    p['partial_steps']=int(st.number_input('Partial diffusion temperature',1,200,int(p['partial_steps']),key=prefix+'-partial'))
                    p['fixed_sequence_positions']=st.text_input('Fix sequence positions (0-based)',p['fixed_sequence_positions'],key=prefix+'-fix-sequence',placeholder='172-205')
                if task in {'symmetry','symmetric_motif'}:p['oligomer_contacts']=st.checkbox('Use oligomer contacts potential',bool(p['oligomer_contacts']),key=prefix+'-contacts')
                p['potential']=st.text_input('Potential',p['potential'],key=prefix+'-potential',placeholder='type:binder_ROG,binderlen:100,weight:5,min_dist:15')
                p['guide_scale']=float(st.number_input('Guide scale',0.,100.,float(p['guide_scale']),.5,key=prefix+'-guide-scale'))
                p['guide_decay']=st.selectbox('Guide decay',['constant','linear','quadratic','cubic'],index=['constant','linear','quadratic','cubic'].index(p['guide_decay']),key=prefix+'-decay')
                p['deterministic']=st.checkbox('Deterministic sampling',bool(p['deterministic']),key=prefix+'-deterministic')
                p['save_trajectory']=st.checkbox('Save diffusion trajectory',bool(p['save_trajectory']),key=prefix+'-trajectory')
        if engine=='rfdiffusion' and task=='fold':
            a,b=st.columns(2)
            with a:p['target_ss_file']=_file('Target secondary structure file',p.get('target_ss_file'),prefix+'-target-ss',owner_id)
            with b:p['target_adj_file']=_file('Target adjacency matrix file',p.get('target_adj_file'),prefix+'-target-adj',owner_id)
            if p['scaffold_preset']=='custom':
                a,b=st.columns(2)
                with a:p['scaffold_ss_file']=_file('Binder secondary structure file',p.get('scaffold_ss_file'),prefix+'-binder-ss',owner_id)
                with b:p['scaffold_adj_file']=_file('Binder adjacency matrix file',p.get('scaffold_adj_file'),prefix+'-binder-adj',owner_id)
            a,b,c=st.columns(3)
            for column,key,label in [(a,'sampled_insertion','Extra loop residues'),(b,'sampled_n','Extra N-terminal residues'),(c,'sampled_c','Extra C-terminal residues')]:
                p[key]=int(column.number_input(label,0,100,int(p[key]),key=prefix+'-'+key))
    jobs=1
    if include_jobs:
        with st.expander('Batch mode',expanded=st.session_state.get(key_prefix+'-show-batch',False)):
            jobs=int(st.number_input('Number of independent jobs',1,100,1,key=key_prefix+'-jobs'))
            st.caption(f"{jobs*p['num_designs']:,} requested designs. Each job is scheduled independently; fixed seeds advance between jobs.")
    p=public_parameters(engine,p)
    if not error:
        try:
            prepared=inputs.prepare_request(engine,p)
            if engine=='rfdiffusion' and prepared.get('contig'):st.caption('Contigs preview: ['+prepared['contig']+']')
        except (ValueError,TypeError,OSError) as exc:error=str(exc)
    return p,jobs,error


def render_setup(schema,*,owner_id,project_selector,submit):
    engine=next(k for k,v in inputs.WORKFLOWS.items() if v==schema['workflow']);spec=inputs.ENGINES[engine]
    prefix='generator-'+engine
    st.html(STYLE)
    st.html('<style>.block-container{max-width:1120px!important;padding-top:1.4rem!important}</style>')
    with st.container(key='generator-surface'):
        if st.button('← Back to All tools',key=prefix+'-back'):
            st.session_state['page']='All tools';st.session_state.pop('selected_workflow',None);st.rerun()
        st.html('<div class="gen-heading"><div class="gen-eyebrow">DESIGN ENGINE</div><h1>'+escape(spec['title'])+'</h1><p>'+escape(spec['description'])+'</p><span class="gen-version">'+escape(spec['version'])+'</span></div>')
        a,b,c,d=st.columns(4)
        if a.button('Show settings',icon=':material/tune:',key=prefix+'-settings-link',use_container_width=True):st.session_state[prefix+'-show-settings']=True
        if b.button('Batch mode',icon=':material/stacks:',key=prefix+'-batch-link',use_container_width=True):st.session_state[prefix+'-show-batch']=True
        c.link_button('Source code',spec['repo'],icon=':material/open_in_new:',use_container_width=True)
        if engine!='protpardelle':d.link_button('Input reference','https://app.tamarind.bio/tools/'+engine,icon=':material/open_in_new:',use_container_width=True)
        else:d.link_button('Examples',spec['repo']+'/tree/main/examples/sampling',icon=':material/science:',use_container_width=True)
        with st.expander('About this workflow'):
            st.write('This tile generates structures. Choose sequence redesign, structure prediction and result analysis on the pipeline canvas.')
            st.markdown('[Author documentation]('+spec['repo']+')')
            st.caption('Only inference controls supported by the installed image are shown. The original native pipeline remains available to existing API jobs.')
        a,b=st.columns(2)
        name=a.text_input('Job name',engine+'-generation',key=prefix+'-name').strip()
        with b:project=project_selector(prefix)
        parameters,jobs,error=render_inputs(engine,owner_id=owner_id,key_prefix=prefix)
        parameters['campaign_name']=name or engine+'-generation'
        with st.expander('Scheduler & local readiness'):
            a,b=st.columns(2)
            gpu=a.selectbox('GPU assignment',['Automatic','GPU 0','GPU 1'],key=prefix+'-gpu')
            parameters['gpu_device']={'Automatic':'all','GPU 0':'0','GPU 1':'1'}[gpu]
            parameters['cpu_cores']=int(b.number_input('CPU cores per job',1,256,4,key=prefix+'-cpus'))
            priority=st.checkbox('Prioritize in scheduler queue',key=prefix+'-priority')
            ready=inputs.readiness(engine,parameters)
            if ready['ready']:st.success('Generator reference files are available.')
            else:st.warning('Missing reference files: '+', '.join(ready['missing']))
        if error:st.caption(error)
        if st.button('Submit '+str(jobs)+' job'+('s' if jobs!=1 else ''),key=prefix+'-submit',type='primary',use_container_width=True,
            disabled=not project or bool(error) or not ready['ready']):
            for index in range(jobs):
                p=copy.deepcopy(parameters)
                if jobs>1:p['campaign_name']+=f'-{index+1:03d}'
                if p.get('seed'):p['seed']=(p['seed']+index)%2147483647 or 1
                submit(schema,schema['workflow'],p,project,{'prioritize':priority},spec['title'])
        if not project:st.caption('Select a project to submit.')


def render_results(run_id,run_dir,candidates,parameters,result):
    engine=parameters.get('generator_engine','rfdiffusion');spec=inputs.ENGINES[engine]
    st.subheader(spec['title']+' · generated structures')
    st.metric('Generated designs',len(candidates))
    st.caption('Generation outputs are ready for sequence redesign and evaluation. They have no folding confidence or filter pass yet.')
    marker=read_json(Path(run_dir)/'molecular_design_studio.json')
    if candidates and result.get('success') is True and marker.get('project_id'):
        if st.button('Continue on pipeline canvas',key='generator-continue-'+run_id,type='primary'):
            from mn_molecular_design_studio.pipeline_graph import new_graph
            graph=new_graph(engine+'_generated')
            graph['nodes'][0]['data']['params'].update(source_run_dir=str(run_dir),selected_candidate_ids=[c['candidate_id'] for c in candidates])
            project=marker['project_id'];st.session_state['pipeline-pending-draft']={'project_id':project,'graph':graph}
            st.session_state.setdefault('project-pending-selection',{})['visual-pipeline-project-id']=project
            st.session_state['page']='Pipelines';st.rerun()
    if not candidates:st.info('Structures will appear when generation finishes.');return
    st.dataframe([{'Design':c['candidate_id'],'Designed chains':', '.join(c.get('binder_chains') or []),
                   'Target chains':', '.join(c.get('target_chains') or []),'Stage':'Generated structure'} for c in candidates],hide_index=True,use_container_width=True)
    index=st.selectbox('Design',range(len(candidates)),format_func=lambda i:candidates[i]['candidate_id'],key='generator-result-'+run_id)
    candidate=candidates[index];path=resolve_stored_path(candidate.get('complex_pdb'),run_dir=Path(run_dir),must_exist=True)
    if warning:=candidate.get('raw_metadata',{}).get('motif_mapping_warning'):st.warning(warning)
    if path:
        molstar_custom_component([StructureVisualization(pdb=path.read_text(),color='chain-id',representation_type='cartoon')],
            key='generator-view-'+run_id+'-'+str(index),height=520,show_structure_tools=False)
        st.download_button('Download generated structure',path.read_bytes(),file_name=path.name,key='generator-download-'+run_id)
    with st.expander('Generation details'):st.json(candidate.get('raw_metadata') or {})
