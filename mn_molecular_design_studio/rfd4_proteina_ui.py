"""Tamarind-style input page for the RFD4-Proteina generator."""
from __future__ import annotations

import shutil
import uuid
from pathlib import Path

import streamlit as st

from mn_molecular_design_studio import rfd4_proteina_input as inputs
from mn_molecular_design_studio.auth import user_input_root
from mn_molecular_design_studio.target_editor import render_target_editor, render_target_source
from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer import (
    StructureVisualization,
    molstar_custom_component,
)


EXAMPLE_SPEC = """# Protein binder condition spec. Upload target.cif below.
pd_l1_binder:
  input: target.cif
  contig: \"60-120, /0, A17-131\"
  conditions:
    C_CRDGRP:
      1:
        - select: \"sel('A/*/*/[N,CA,C,O]')\"
    C_SEQ:
      true:
        - select: \"sel('A')\"
    C_HOT:
      true:
        - select: \"sel('A/*/[56,115,123]')\"
"""

PREVIEW_TASKS = {
    "idr_peptide_binder": {
        "label": "IDR Peptide Binder (experimental)",
        "description": "Co-generate a binder and a sequence-fixed, coordinate-unconditioned IDR target.",
    },
    "cyclic_peptide_binder": {
        "label": "Cyclic Peptide Binder (preview)",
        "description": "Plan a head-to-tail cyclic peptide binder against a protein target.",
    },
}


def _load_example_inputs(example: str, owner_id: str) -> None:
    # The callback runs before the next render. Clear stale widget state first
    # so this example replaces prior choices instead of only changing defaults.
    for key in list(st.session_state):
        if key.startswith('rfd4-') and not key.endswith('-project-id'):
            del st.session_state[key]

    if example == "pd_l1_binder":
        task, input_mode = "protein_binder", "guided_binder"
        example_target = Path(__file__).resolve().parent / "assets" / "proteina_complexa" / "PD-L1.pdb"
        stored_target = user_input_root(owner_id) / "rfd4-proteina" / "PD-L1.pdb"
        stored_target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(example_target, stored_target)
        from mn_molecular_design_studio.design_validation import read_structure, protein_chains

        residues_by_chain = protein_chains(read_structure(stored_target))
        chain = "A" if "A" in residues_by_chain else next(iter(residues_by_chain), None)
        if chain is None:
            raise ValueError("The bundled PD-L1 example contains no protein chain.")
        available = {residue.id[1] for residue in residues_by_chain[chain]}
        hotspot_numbers = (40, 99, 107) if chain == "A" else ()
        hotspots = ",".join(f"{chain}{number}" for number in hotspot_numbers if number in available)
        st.session_state.update({
            "rfd4-task-pills": task,
            "rfd4-input-mode-pills-protein_binder": input_mode,
            "rfd4-target-format": "PDB",
            "rfd4-target-pdb-code": "5O45",
            "rfd4-target-fetched-target": str(stored_target),
            "rfd4-target-fetched-target-name": "PD-L1",
            "rfd4-target-trim": False,
            "rfd4-target-view": False,
            "rfd4-target-chains-protein_binder": [chain],
            "rfd4-min-protein_binder": 60,
            "rfd4-max-protein_binder": 120,
            "rfd4-hotspots-protein_binder": hotspots,
            "rfd4-job-name": f"rfd4-proteina-example-{uuid.uuid4().hex[:5]}",
            "rfd4-count": 8,
            "rfd4-batch-size": 1,
            "rfd4-cpu-cores": 8,
            "rfd4-compile-model": False,
            "rfd4-prioritize": False,
        })
        notice = "PD-L1 example loaded: target structure, chain, hotspots, and design settings are filled. Review the inputs before continuing."
    else:
        task, input_mode = "unconditional", "unconditional"
        st.session_state.update({
            "rfd4-task-pills": task,
            "rfd4-input-mode-pills-unconditional": input_mode,
            "rfd4-job-name": "rfdiffusion4-unconditional-example",
            "rfd4-count": 1,
            "rfd4-min-unconditional": 100,
            "rfd4-max-unconditional": 100,
            "rfd4-batch-size": 1,
            "rfd4-cpu-cores": 8,
            "rfd4-compile-model": False,
            "rfd4-prioritize": False,
        })
        notice = "Unconditional monomer example loaded. Choose a project before submitting."

    st.session_state["rfd4-example-loaded-notice"] = notice


def render_setup(schema: dict, *, owner_id: str, project_selector, submit, folding=False) -> None:
    st.html("""<style>
    .block-container {max-width:732px!important;padding:1rem 1.5rem 2rem!important;margin:0 auto!important}
    [data-testid="stAppViewContainer"] {background:#f6f6f6!important} [data-testid="stMainBlockContainer"] {background:transparent;border:0;border-radius:0;box-shadow:none;padding:1.5rem 1.5rem 2rem}
    .st-key-rfd4-surface {--rfd4-green:#28674d;--rfd4-ink:#202426;--rfd4-mint:#edf5f0}
    .st-key-rfd4-surface .rfd4-heading {text-align:center;padding:.5rem 0 .55rem}
    .st-key-rfd4-surface .rfd4-heading h1 {font-size:24px;color:#17191b;margin:0;font-weight:700;letter-spacing:-.025em}
    .st-key-rfd4-surface .rfd4-heading .badge {display:inline-block;vertical-align:middle;margin-left:.35rem;padding:.16rem .48rem;border-radius:1rem;background:#28674d;color:#fff;font-size:.62rem;font-weight:700;letter-spacing:.02em}
    .st-key-rfd4-surface .rfd4-heading .hero-link {display:inline-block;vertical-align:middle;margin-left:.4rem;color:#68756e;text-decoration:none;font-size:.78rem;font-weight:650}
    .st-key-rfd4-surface .rfd4-heading .hero-link:hover {color:#28674d}
    .st-key-rfd4-surface .rfd4-heading p {font-size:14px;color:#4c4c4c;margin:1rem auto .45rem;max-width:700px;line-height:1.5}
    .st-key-rfd4-surface .rfd4-meta {font-size:12px;color:#777;margin-top:1.15rem}
    .st-key-rfd4-surface .rfd4-project-reminder {font-size:.78rem;color:#7b858b;margin:-.15rem 0 .65rem}
    .st-key-rfd4-surface [data-testid="stVerticalBlockBorderWrapper"] {border-color:#d9ddda!important;border-radius:9px!important}
    .st-key-rfd4-surface [data-testid="stButton"] button {border-radius:7px!important}
    .st-key-rfd4-surface button[kind="primary"] {background:#28674d!important;border-color:#28674d!important;color:#fff!important;min-height:2.6rem}
    .st-key-rfd4-surface button[kind="primary"]:disabled {background:#91b5a4!important;border-color:#91b5a4!important;color:#fff!important}
    .st-key-rfd4-surface [data-testid="stFileUploaderDropzone"] {background:#f3f8f5!important;border:1px dashed #b9d2c2!important;border-radius:8px}
    .st-key-rfd4-surface [data-testid="stExpander"] {border-color:#d9ddda!important;border-radius:8px!important;background:#fff}
    .st-key-rfd4-surface [data-testid="stTextArea"] textarea {font-family:ui-monospace,SFMono-Regular,Consolas,monospace;font-size:.83rem}
    .st-key-rfd4-surface [data-testid="stTabs"] button[aria-selected="true"] {color:#28674d!important;border-bottom-color:#28674d!important}
    .st-key-rfd4-surface button[data-variant="pills"] {border-radius:0!important;border:1px solid #d5d5d5!important;background:transparent!important;color:#686868!important;padding:7px 12px!important;min-height:38px!important;font-weight:400!important}
    .st-key-rfd4-surface button[data-variant="pills"][data-selected="true"] {background:#28674d!important;color:white!important;border-color:#28674d!important}
    .st-key-rfd4-surface [data-testid="stWidgetLabel"] p {font-size:16px!important;font-weight:600!important;color:#202426!important}
    .st-key-rfd4-surface [data-testid="stPills"] [role="radiogroup"] {gap:4px!important}
    .st-key-rfd4-surface [data-testid="stExpander"] {background:transparent!important}
    .st-key-rfd4-surface [data-testid="stPopover"] button {background:transparent;border:0;font-size:12px;color:#686868}
    .st-key-rfd4-surface input {font-size:14px!important}
    .st-key-rfd4-surface [data-baseweb="input"], .st-key-rfd4-surface [data-baseweb="textarea"], .st-key-rfd4-surface [data-baseweb="select"] > div {background:transparent!important;border:1px solid #d2d2d2!important;border-radius:5px!important}
    .st-key-rfd4-surface [data-baseweb="base-input"], .st-key-rfd4-surface input {background:transparent!important}
    .st-key-rfd4-surface [data-testid="stNumberInput"] button {display:none}
    .st-key-rfd4-surface [data-testid="stFileUploaderDropzone"] {min-height:100px;justify-content:center;background:transparent!important;border:1px dashed #c8c8c8!important}
    .st-key-rfd4-surface button[data-variant="pills"] p {font-size:16px!important}
    .st-key-rfd4-surface [data-testid="stTextInputRootElement"], .st-key-rfd4-surface [data-testid="stTextAreaRootElement"], .st-key-rfd4-surface [data-testid="stNumberInputContainer"], .st-key-rfd4-surface [data-testid="stMultiSelect"] [role="group"]:has(> [data-testid="stMultiSelectTagsContainer"]), .st-key-rfd4-surface [data-testid="stSelectbox"] [role="combobox"] {background:transparent!important;border:1px solid #d2d2d2!important;border-radius:5px!important}
    .st-key-rfd4-surface [data-tag] {background:#28674d!important}
    </style>""")
    with st.container(key="rfd4-surface"):
        _render_content(schema, owner_id=owner_id, project_selector=project_selector, submit=submit, folding=folding)


def _retain_pill_choice(key,default):
    if st.session_state.get(key) is None:st.session_state[key]=st.session_state.get(key+'-last-choice',default)
    else:st.session_state[key+'-last-choice']=st.session_state[key]


def _pills(label, choices, default, *, key, labels=None, help=None):
    if st.session_state.get(key) in choices:st.session_state[key+'-last-choice']=st.session_state[key]
    return st.pills(label,choices,default=default,key=key,format_func=(labels or {}).get if labels else str,help=help,label_visibility='collapsed' if label=='File format' else 'visible',on_change=_retain_pill_choice,args=(key,default)) or default


def _remove_object(key,index):
    count=st.session_state.get(key+'-rows',0)
    for j in range(index,count-1):
        for suffix in ('low','high','type','atoms','aa','probability'):
            nxt=f'{key}-{j+1}-{suffix}'
            if nxt in st.session_state:st.session_state[f'{key}-{j}-{suffix}']=st.session_state[nxt]
    st.session_state[key+'-rows']=max(0,count-1)


def _objects(label,key,kind='range'):
    st.markdown(f'**{label}**')
    count_key=key+'-rows'; st.session_state.setdefault(count_key,0)
    values=[]
    for i in range(st.session_state[count_key]):
        with st.container(border=True):
            if kind=='range':
                a,b=st.columns(2)
                maximum=40 if 'cdr' in key else 1024
                low=a.number_input('Lower Bound',min_value=1,max_value=maximum,value=None,key=f'{key}-{i}-low')
                high=b.number_input('Upper Bound',min_value=1,max_value=maximum,value=None,key=f'{key}-{i}-high')
                if low is not None or high is not None: values.append({'lower_bound':low if low is not None else high,'upper_bound':high if high is not None else low})
            elif kind=='interaction':
                interaction=st.selectbox('Interaction',list(inputs.LIGAND_INTERACTIONS),key=f'{key}-{i}-type')
                atoms=st.text_input('Ligand Atoms',key=f'{key}-{i}-atoms',placeholder='e.g. SAM:N3, SAM:O4')
                values.append({'interaction':interaction,'ligand_atoms':atoms})
            else:
                aa=st.selectbox('Amino Acid',list('ACDEFGHIKLMNPQRSTVWY'),key=f'{key}-{i}-aa')
                probability=st.number_input('Probability',0.,1.,.1,key=f'{key}-{i}-probability')
                values.append({'amino_acid':aa,'probability':probability})
            st.button('Delete Object',key=f'{key}-{i}-delete',on_click=_remove_object,args=(key,i))
    if st.button('Add New Object',key=key+'-add'):
        st.session_state[count_key]+=1; st.rerun()
    return values


def _source(owner_id,key='rfd4-target',label='Input Structure',multiple=False):
    import re
    from urllib.request import Request,urlopen
    from types import SimpleNamespace
    st.markdown(f'**{label}**')
    if label=='Input Structure': st.caption('The target, motif or backbone to design against or around. PDB or mmCIF (.cif); use mmCIF with ligands or a motif, since a PDB cannot hold bond orders or charges.')
    extension=_pills('File format',['PDB','CIF'],'PDB',key=key+'-format')
    a,b=st.columns([3,1],vertical_alignment='bottom')
    code=a.text_input('Enter PDB Code',key=key+'-pdb-code',placeholder='Enter PDB Code').strip().upper()
    if b.button('Fetch PDB',disabled=not code,key=key+'-fetch',width='stretch'):
        try:
            if not re.fullmatch(r'[A-Z0-9]{4}',code): raise ValueError('Enter a four-character PDB code.')
            suffix='cif' if extension=='CIF' else 'pdb'
            with urlopen(Request(f'https://files.rcsb.org/download/{code}.{suffix}',headers={'User-Agent':'Molecular-Design-Studio'}),timeout=30) as response: data=response.read(50_000_001)
            paths=inputs.store_structure_uploads([SimpleNamespace(name=code+'.'+suffix,getvalue=lambda:data)],owner_id)
            st.session_state[key+'-fetched-target']=str(paths[0])
        except (ValueError,OSError) as exc:st.error(str(exc))
    uploaded=st.file_uploader(f'Drag & Drop or Browse your {extension.lower()} file',type=['pdb','ent'] if extension=='PDB' else ['cif','mmcif'],accept_multiple_files=multiple,key=key+'-target-upload',label_visibility='collapsed',max_upload_size=50)
    files=(uploaded or []) if multiple else ([uploaded] if uploaded else [])
    paths=[str(p) for p in inputs.store_structure_uploads(files,owner_id)]
    if not paths:
        fetched=st.session_state.get(key+'-fetched-target')
        if fetched and Path(fetched).is_file():paths=[str(fetched)]; st.caption(Path(fetched).name)
    if not multiple and paths:
        a,b=st.columns([1,1])
        trim=a.checkbox('Trim '+extension.lower(),key=key+'-trim')
        if b.toggle('View',key=key+'-view'):
            molstar_custom_component([StructureVisualization(pdb=Path(paths[0]).read_text(),color='molecule-type',representation_type='cartoon')],key=key+'-viewer',height=340)
        if trim: paths=[_trim_source(Path(paths[0]),owner_id,key)]
    return paths


def _trim_source(path,owner_id,key):
    import gemmi,hashlib
    from mn_molecular_design_studio.design_validation import read_structure,protein_chains
    from Bio.PDB import MMCIFIO,Select
    chains=protein_chains(read_structure(path)); ranges={}
    for chain,residues in chains.items():
        a,b=st.columns(2)
        low=a.number_input(f'Chain {chain} start',value=int(residues[0].id[1]),key=key+'-trim-'+chain+'-low')
        high=b.number_input(f'Chain {chain} end',value=int(residues[-1].id[1]),key=key+'-trim-'+chain+'-high')
        if low>high: raise ValueError('Trim start must precede trim end.')
        ranges[chain]=(low,high)
    if path.suffix.lower() in {'.pdb','.ent'}:
        class Retained(Select):
            def accept_residue(self,residue):
                bounds=ranges.get(residue.parent.id)
                return residue.id[0]!=' ' or bounds is None or bounds[0]<=residue.id[1]<=bounds[1]
        import io
        output=io.StringIO(); writer=MMCIFIO();writer.set_structure(read_structure(path));writer.save(output,Retained());data=output.getvalue().encode()
    else:
        document=gemmi.cif.read_file(str(path));block=document.sole_block();atoms=block.get_mmcif_category('_atom_site.')
        keep=[]
        for i,chain in enumerate(atoms['auth_asym_id']):
            number=atoms['auth_seq_id'][i]; bounds=ranges.get(chain)
            keep.append(bounds is None or not str(number).lstrip('-').isdigit() or bounds[0]<=int(number)<=bounds[1])
        block.set_mmcif_category('_atom_site.',{name:[v for v,k in zip(values,keep) if k] for name,values in atoms.items()})
        # Full sequence records would make parsers restore trimmed residues as
        # unresolved positions. Keep chemistry tables and polymer types intact.
        for category in ('_entity_poly_seq.','_pdbx_poly_seq_scheme.'):
            block.set_mmcif_category(category,{})
        polymer=block.get_mmcif_category('_entity_poly.')
        for name in ('pdbx_seq_one_letter_code','pdbx_seq_one_letter_code_can'):polymer.pop(name,None)
        if polymer:block.set_mmcif_category('_entity_poly.',polymer)
        data=document.as_string().encode()
    destination=user_input_root(owner_id)/'rfd4-proteina'/hashlib.sha256(data).hexdigest()[:16]/'trimmed.cif'
    destination.parent.mkdir(parents=True,exist_ok=True);destination.write_bytes(data)
    return str(destination)


def _chain_choices(path):
    if not path:return []
    from mn_molecular_design_studio.design_validation import read_structure,protein_chains
    return list(protein_chains(read_structure(Path(path))))


def _residues(title,path,chains,key):
    st.markdown(f'**{title}**')
    if not path or not chains:
        st.caption('Please upload a structure and select chains to enable residue selection')
        return ''
    from mn_molecular_design_studio.target_editor import TargetDesignInput,_render_residue_viewer,_increment_session_counter
    from mn_molecular_design_studio.design_validation import read_structure,protein_chains
    rows=protein_chains(read_structure(Path(path))); valid={(c,r.id[1]) for c in chains for r in rows.get(c,[])}
    st.session_state.setdefault(key,'')
    revision=key+'-revision';st.session_state.setdefault(revision,0)
    pending=key+'-pending-viewer-selection'
    if pending in st.session_state: st.session_state[key]=st.session_state.pop(pending)
    a,b,c=st.columns([5,1,1],vertical_alignment='bottom')
    if b.button('Select all',key=key+'-all'):
        st.session_state[key]=','.join(f'{chain}:{number}' if not chain.isalpha() else f'{chain}{number}' for chain,number in sorted(valid));_increment_session_counter(revision);st.rerun()
    if c.button('Clear',key=key+'-clear',disabled=not st.session_state[key]): st.session_state[key]='';_increment_session_counter(revision);st.rerun()
    value=a.text_input('Chosen Residues',key=key,placeholder='e.g. A25, A30-35',on_change=_increment_session_counter,args=(revision,))
    target=TargetDesignInput(label=title,key_prefix=key,source_path=Path(path),design_path=Path(path),chains=chains,
                            source_text=Path(path).read_text(),design_text=Path(path).read_text(),hotspot_state_key=key,valid_residues=valid)
    _render_residue_viewer(target,kind='hotspots',height=280,show_structure_tools=False)
    return value


def _optional(task,o,p,path,chains,folding,*,condition_spec=False):
    antibody=task in ('nanobody','antibody','tcr')
    design_type=o.get('score_design_type',task)
    if task=='score_structures':
        o['refold']=st.checkbox('Check Designs by Refolding',value=True,key='rfd4-score-refold')
    else:
        if task in ('protein_binder','small_molecule_binder','dna_binder','unconditional'):o['cath']=st.text_input('Fold Class (CATH)',placeholder='e.g. 2.60.40',key='rfd4-cath-'+task)
        if antibody:o['cdr_lengths']=_objects('CDR Lengths','rfd4-cdr-lengths-'+task)
        if task in ('protein_binder','nanobody','antibody','tcr') and not condition_spec:
            o['hotspot_atoms']=_pills('Hotspot Atoms',['whole','sidechain','tip','single'],'whole',key='rfd4-hotspot-atoms-'+task,labels={'whole':'Whole residue','sidechain':'Side chain','tip':'Side-chain tip','single':'Single atom'})
            st.caption({'whole':'Mark all resolved atoms in each chosen residue.','sidechain':'Mark side-chain atoms; use CA for glycine.','tip':'Mark atoms beyond CB; use CB for alanine and CA for glycine.','single':'Mark the resolved side-chain atom farthest from CA.'}[o['hotspot_atoms']])
        if task in ('protein_binder','nanobody','antibody','tcr','small_molecule_binder','dna_binder') and not condition_spec:
            o['repack_target']=st.checkbox('Repack Target Side Chains',value=True,key='rfd4-repack-'+task)
        if antibody:
            o['framework_atoms']=_pills('Framework Atoms To Hold',['backbone','all'],'backbone',key='rfd4-framework-atoms-'+task,labels={'backbone':'Backbone only','all':'All atoms'})
            if task in ('antibody','tcr'):o['dock_independently']=st.checkbox('Dock Chains Independently',key='rfd4-dock-'+task)
            if task!='nanobody' or o.get('nanobody_framework')=='Your own nanobody':o['keep_constant_domains']=st.checkbox('Keep Constant Domains',key='rfd4-constant-'+task)
        if task in ('motif_scaffolding','enzyme'):
            o['scaffold_lengths']=_objects('Scaffold Lengths','rfd4-scaffold-lengths-'+task)
            if task=='enzyme':
                p['length_min']=st.number_input('Minimum Scaffold Length',20,1024,value=None,key='rfd4-scaffold-min')
                p['length_max']=st.number_input('Maximum Scaffold Length',20,1024,value=None,key='rfd4-scaffold-max')
                o['scaffold_total_explicit']=p['length_min'] is not None or p['length_max'] is not None
                if p['length_min'] is None:p['length_min']=p['length_max'] or 180
                if p['length_max'] is None:p['length_max']=p['length_min']
                o['active_site_atoms']=st.text_input('Active Site Atoms By Residue',placeholder='A102:CG,OD2 A195:NE2')
        if task in ('enzyme','small_molecule_binder'):
            o['ligands_to_keep']=st.text_input('Ligands To Keep',placeholder='SAM, ATP, MG',key='rfd4-ligands-'+task)
            o['ligand_atoms_to_hold']=st.text_input('Ligand Atoms To Hold',placeholder='SAM:O2, SAM:C5',key='rfd4-ligandhold-'+task)
            if task=='enzyme':
                o['ligand_hbond_donors']=st.text_input('Ligand H-Bond Donors',placeholder='SAM:N3')
                o['ligand_hbond_acceptors']=st.text_input('Ligand H-Bond Acceptors',placeholder='SAM:O2')
            else:
                o['ligand_burial']=_pills('Ligand Burial',['model','buried','intermediate','exposed'],'model',key='rfd4-burial',labels={'model':'Model decides','buried':'Buried','intermediate':'Partly buried','exposed':'Exposed'})
                o['exposed_ligand_atoms']=st.text_input('Exposed Ligand Atoms',placeholder='SAM:N6')
                o['ligand_interactions']=_objects('Ligand Interactions','rfd4-ligand-interactions','interaction')
        if task=='existing_structure':
            if o.get('operation')=='diversify':
                o['noise_mode']=_pills('Noise Mode',['fraction','levels'],'fraction',key='rfd4-noise-mode',labels={'fraction':'Fraction','levels':'Levels'})
                if o['noise_mode']=='levels':
                    o['backbone_level']=st.number_input('Backbone Level',.01,.99,.8)
                    o['sequence_level']=st.number_input('Sequence Level',.01,.99,.8)
                o['redesign_shell']=st.number_input('Redesign Shell (Angstroms)',1.,30.,value=None)
            if o.get('operation') in ('diversify','inverse_folding'):
                o['hold_kept_in_place']=st.checkbox('Hold Kept Residues In Place',value=True)
            if o.get('operation')=='diversify':
                for key,title in [('helix_residues','Helix Residues'),('strand_residues','Strand Residues'),('loop_residues','Loop Residues')]:o[key]=_residues(title,path,chains,'rfd4-'+key)
        o['sampling_effort']=_pills('Sampling Effort',['Fast','Balanced','Thorough'],'Balanced',key='rfd4-effort-'+task)
        p['nsteps']=st.number_input('Sampling Steps',60,120,60,key='rfd4-steps-'+task)
        if task!='rna_design' and not folding:o['reward_search']=_pills('Reward Search',['Off','Interface contacts','Low clashes','Best by refolding'],'Off',key='rfd4-reward-'+task)
        o['deterministic']=st.checkbox('Deterministic Sampling',key='rfd4-deterministic-'+task)
        o['backbone_noise']=st.number_input('Backbone Noise',.01,1.5,1.,key='rfd4-bbnoise-'+task)
        o['latent_noise']=st.number_input('Latent Noise',.01,5.,1.,key='rfd4-latnoise-'+task)
        p['seed']=st.number_input('Random Seed',min_value=0,max_value=4294967295,value=None,placeholder='Automatic',key='rfd4-seed-'+task)
        for name,default,minimum in [('backbone_tsr_k',16.,1.),('backbone_tsr_sigma',3.,.01),('latent_tsr_k',4.,1.),('latent_tsr_sigma',12.,.01)]:
            label=name.replace('_',' ').title().replace('Tsr','TSR').replace(' K',' k')
            o[name]=st.number_input(label,min_value=minimum,value=default,key='rfd4-'+name+'-'+task)
        if task not in ('rna_design','custom_conditions') and not folding:
            o['exclude_amino_acids']=st.multiselect('Exclude Amino Acids',list('ACDEFGHIKLMNPQRSTVWY'),key='rfd4-exclude-'+task)
        if task in ('protein_binder','small_molecule_binder','dna_binder','unconditional'):
            o['favor_amino_acids']={r['amino_acid']:r['probability'] for r in _objects('Favor Amino Acids','rfd4-favor-'+task,'favor')}
        if not folding:
            o['cluster']=st.checkbox('Cluster Designs',key='rfd4-cluster-'+task)
            o['filter_conditions']=st.checkbox('Drop Designs That Break Conditions',value=True,key='rfd4-filter-'+task)
            if task not in ('antibody','tcr','rna_design'):o['refold']=st.checkbox('Check Designs by Refolding',value=True,key='rfd4-refold-'+task)
            else:o['refold']=False
    if o.get('refold'):
        fields=[]
        if design_type in ('protein_binder','nanobody','antibody','tcr','small_molecule_binder','dna_binder','enzyme'):fields.append(('interface_error_limit','Interface Error Limit'))
        if design_type in ('protein_binder','nanobody','antibody','tcr','dna_binder'):fields.append(('binder_rmsd_limit','Binder RMSD Limit'))
        if design_type in ('unconditional','small_molecule_binder','enzyme','motif_scaffolding'):fields.append(('refold_backbone_rmsd_limit','Refold Backbone RMSD Limit'))
        if design_type in ('enzyme','motif_scaffolding'):fields.append(('motif_rmsd_limit','Motif RMSD Limit'))
        for key,label in fields:o[key]=st.number_input(label,.1,30.,value=None,placeholder="Authors' limit",key='rfd4-'+key+'-'+task)
    if task in ('protein_binder','small_molecule_binder','dna_binder'):
        o['apo']=st.checkbox('Check Binders Fold Alone',key='rfd4-apo-'+task)
    if task in ('protein_binder','small_molecule_binder'):
        o['energy']=st.checkbox('Score Interface Energy',key='rfd4-energy-'+task)
    if task not in ('antibody','tcr','rna_design','score_structures') and not folding:o['developability']=st.checkbox('Score Developability',key='rfd4-developability-'+task)
    if task!='score_structures':o['save_trajectory']=st.checkbox('Save Trajectory',key='rfd4-trajectory-'+task)
    if o.get('refold') or o.get('apo') or o.get('reward_search')=='Best by refolding':
        from mn_molecular_design_studio.engine.workflow_engine.workflows.rfd4_oracle import compatible_engines
        engines=compatible_engines(design_type)
        o['folding_engine']=st.selectbox('Folding Engine',list(engines),format_func=engines.get,key='rfd4-engine-'+task)
        o['folding_recycles']=st.number_input('Folding Recycles',1,100,5,key='rfd4-recycles-'+task)
        if o['folding_engine']=='alphafold3':o['accept_alphafold3_terms']=st.checkbox('Accept official AlphaFold3 weight terms')


def _render_content(schema:dict,*,owner_id,project_selector,submit,folding=False):
    if st.button('← Back to All tools',key='rfd4-back'):
        st.session_state['page']='All tools';st.session_state.pop('selected_workflow',None);st.rerun()
    title='RFD4-Proteina Structure Prediction' if folding else 'RFD4-Proteina'
    st.html(f'<div class="rfd4-heading"><h1>{title} <span class="badge">NEW</span> <a class="hero-link" title="Source code" href="https://github.com/RosettaCommons/RFD4-Proteina">⌘</a> <a class="hero-link" title="Paper" href="https://rfd4-proteina.org/">▤</a></h1>'
            f'<div class="rfd4-meta">{"Structure Prediction" if folding else "Protein Design"}　◷ ~10m　<i>v0.0.3</i></div>'
            '<p>Generates protein structure and sequence together, from scratch or around a target, motif or ligand.</p></div>')
    for col,label,text in zip(st.columns([1,1,1,.8]),['Show Outputs','Show Settings','API Docs','Details'],[
        'CIF structures, sequences and per-design metrics.','The task selects its matching controls.','[RFD4 API and documentation](https://rfd4-proteina.org/latest/)',
        'RFD4 native generation and sequence-conditioned folding. Select a separate folding engine for independent checks.']):
        with col:
            with st.popover(label):st.markdown(text)
    st.divider()
    if notice:=st.session_state.pop('rfd4-example-loaded-notice',None):st.success(notice)
    a,b=st.columns([1.1,1])
    job_name=a.text_input('Job Name',value='rfd4-structure-prediction' if folding else 'rfd4-proteina-design',key='rfd4-job-name')
    with b:project_id=project_selector('rfd4-fold' if folding else 'rfd4-proteina')
    task_labels={**inputs.TASKS,**{key:value['label'] for key,value in PREVIEW_TASKS.items()}}
    task='existing_structure' if folding else _pills('Design Task',list(task_labels),'protein_binder',key='rfd4-task-pills',labels=task_labels)
    preview_task=task in PREVIEW_TASKS
    if not folding:
        st.markdown(PREVIEW_TASKS[task]['description'] if preview_task else inputs.TASK_HELP[task])
    p={'task':'protein_binder' if task=='idr_peptide_binder' else task,'input_mode':'guided','job_name':job_name,'target_chains':[],'input_files':[],'hotspots':''};o={};p['options']=o
    ready=True;path=None;chains=[]
    if task=='existing_structure':
        o['operation']='structure_prediction' if folding else _pills('Operation',list(inputs.OPERATIONS),'diversify',key='rfd4-operation',labels=inputs.OPERATIONS)
    if task=='custom_conditions':
        p['input_mode']=_pills('Conditions Input',['condition_spec','annotated_cif'],'condition_spec',key='rfd4-custom-mode',labels={'condition_spec':'Condition File','annotated_cif':'Annotated CIF'})
    sequence_input=folding and st.radio('Input',['Protein sequence','Structure'],horizontal=True)=='Protein sequence'
    if task=='idr_peptide_binder':
        p['input_mode']='condition_spec'
        fasta_text=st.text_area('Target FASTA (IDR region)',placeholder='>target_idr\nPaste one IDR sequence',key='rfd4-idr-preview-fasta')
        hotspot_text=st.text_input('Optional hotspot residues',placeholder='e.g. A5, A8-12',key='rfd4-idr-preview-hotspots')
        p['length_min']=st.number_input('Minimum Binder Length',20,1024,70,key='rfd4-idr-min-length')
        p['length_max']=st.number_input('Maximum Binder Length',20,1024,110,key='rfd4-idr-max-length')
        p['target_chains']=['A']
        p['hotspots']=hotspot_text
        ready=False
        if fasta_text.strip():
            try:
                target_sequence=inputs.parse_idr_fasta(fasta_text)
                target_path=inputs.store_idr_sequence_carrier(target_sequence,user_input_root(owner_id))
                p['target_structure']=str(target_path)
                p['input_files']=[str(target_path)]
                p['condition_spec']=inputs.idr_binder_condition_spec(
                    target_sequence,target_filename=target_path.name,
                    binder_length_min=p['length_min'],binder_length_max=p['length_max'],hotspots=hotspot_text,
                )
                ready=p['length_min']<=p['length_max']
                with st.expander('RFD4 condition specification'):
                    st.code(p['condition_spec'],language='yaml')
            except ValueError as exc:
                st.error(str(exc))
        st.info('Experimental sequence-only target mode: the CIF contains an ideal extended backbone only to carry residue identity. C_SEQ fixes the target sequence; no C_CRD or C_CRDGRP condition is added, so target coordinates are free. Hotspots add C_HOT contact guidance without fixing coordinates. This does not predict or validate the IDR ensemble.')
        st.caption('Hotspot positions are 1-based positions in the pasted FASTA sequence. The native RFD4 protein-binder workflow and its configured Results are used.')
    elif preview_task:
        st.file_uploader('Target structure (PDB or mmCIF)',type=['pdb','ent','cif','mmcif'],disabled=True,key='rfd4-cyclic-preview-structure')
        st.multiselect('Target Chains',[],disabled=True,key='rfd4-cyclic-preview-chains')
        st.text_input('Optional hotspot residues',placeholder='e.g. A25, A30-35',disabled=True,key='rfd4-cyclic-preview-hotspots')
        st.number_input('Cyclic peptide length',min_value=6,max_value=24,value=12,disabled=True,key='rfd4-cyclic-preview-length')
        st.markdown('RFD4 has a [cyclic peptide binder preset](https://github.com/RosettaCommons/RFD4-Proteina/blob/main/configs/design_tasks/presets/cyclic_peptide_binder.yaml) and [tests for head-to-tail ring topology](https://github.com/RosettaCommons/RFD4-Proteina/blob/main/tests/pipeline/test_canonical_polymer_ordering.py). Studio still needs a ring-aware input and output path.')
        st.info('Preview only: job submission is disabled until Studio can validate these inputs and preserve the intended design topology in its outputs.')
        ready=False
    elif sequence_input:
        p['protein_sequences']=st.text_area('Protein sequence(s)',placeholder='Paste sequence or FASTA; use separate records for chains.')
        ready=bool(p['protein_sequences'].strip())
    elif task=='rna_design':
        o['rna_secondary_structure']=st.text_area('Secondary Structure (Dot-Bracket)',placeholder='(((....)))')
        ready=bool(o['rna_secondary_structure'].strip())
        st.caption('Use paired brackets for stems, . for unpaired positions, and ? for unspecified positions. Lengths above 256 are beyond the reported training range.')
    elif task!='unconditional':
        multiple=task in ('score_structures','custom_conditions')
        files=_source(owner_id,label='Structures to Score' if task=='score_structures' else 'Input Structure',multiple=multiple)
        ready=bool(files);path=files[0] if files else None;chains=_chain_choices(path)
        if multiple:p['input_files']=files
        else:p['target_structure']=path
        if task=='custom_conditions' and p['input_mode']=='condition_spec':
            upload=st.file_uploader('Condition File',type=['yaml','yml','json'],key='rfd4-custom-yaml')
            p['condition_spec']=st.text_area('Condition specification',value=upload.getvalue().decode() if upload else EXAMPLE_SPEC,height=220)
    if task=='nanobody':
        o['nanobody_framework']=_pills('Nanobody Framework',inputs.FRAMEWORKS,'Vobarilizumab',key='rfd4-framework')
        st.caption('A built-in framework. Upload the target alone.' if o['nanobody_framework']!='Your own nanobody' else 'Upload one structure containing your nanobody and its target.')
    if task=='score_structures':
        labels={
            'unconditional':'Single Protein','protein_binder':'Protein Binder','nanobody':'Nanobody',
            'antibody':'Antibody','tcr':'TCR','motif_scaffolding':'Motif Scaffold','enzyme':'Enzyme',
            'small_molecule_binder':'Small Molecule Binder','dna_binder':'DNA Binder',
        }
        o['score_design_type']=_pills('Design Type',list(labels),'protein_binder',key='rfd4-score-type',labels=labels)
        o['designed_chain']=st.selectbox('Designed Chain',chains,index=None,placeholder='Select a chain',key='rfd4-score-chain') or ''
        ready=ready and bool(o['designed_chain'])
    if task in ('protein_binder','nanobody','antibody','tcr','small_molecule_binder','dna_binder'):
        p['target_chains']=st.multiselect('Target Chains',chains,default=chains[:1] if task in ('protein_binder','nanobody','antibody','tcr') else [],key='rfd4-target-chains-'+task)
    if task in ('antibody','tcr') or (task=='nanobody' and o['nanobody_framework']=='Your own nanobody'):
        o['design_chains']=st.multiselect('Design Chains',chains,key='rfd4-design-chains-'+task)
        o['cdr_residues']=_residues('CDR Residues',path,o['design_chains'],'rfd4-cdr-'+task)
        ready=ready and bool(o['design_chains']) and bool(o['cdr_residues'])
        if set(o['design_chains'])&set(p['target_chains']):st.error('Target and design chains must be separate.');ready=False
    if task in ('protein_binder','small_molecule_binder','dna_binder','unconditional'):
        low,high={'small_molecule_binder':(100,155),'dna_binder':(120,150),'unconditional':(100,100)}.get(task,(60,120))
        suffix='' if task=='unconditional' else ' Binder'
        p['length_min']=st.number_input('Minimum'+suffix+' Length',20,1024,low,key='rfd4-min-'+task)
        p['length_max']=st.number_input('Maximum'+suffix+' Length',20,1024,high,key='rfd4-max-'+task)
        ready=ready and p['length_min']<=p['length_max']
    if task in ('protein_binder','nanobody','antibody','tcr'):
        p['hotspots']=_residues('Hotspot Residues',path,p['target_chains'],'rfd4-hotspots-'+task)
    if task in ('motif_scaffolding','enzyme'):
        o['motif_residues']=_residues('Motif Residues',path,chains,'rfd4-motif-'+task)
        labels={'all':'All atoms','backbone':'Backbone only','sidechain_beyond_cb':'Side chain beyond CB','tip':'Functional tips'}
        order=['sidechain_beyond_cb','tip','all','backbone'] if task=='enzyme' else ['all','backbone','sidechain_beyond_cb']
        if task=='enzyme':labels['all']='All resolved atoms'
        o['motif_atoms']=_pills('Active Site Atoms To Hold' if task=='enzyme' else 'Motif Atoms To Hold',order,'tip' if task=='enzyme' else 'all',key='rfd4-motif-atoms-'+task,labels=labels)
        o['unindex_motif']=st.checkbox('Let The Model Place The Motif',value=task=='enzyme',key='rfd4-motif-place-'+task)
        ready=ready and bool(o['motif_residues'])
    if task in ('enzyme','small_molecule_binder'):
        o['ligand_pose']=_pills('Ligand Pose',['generated','fixed'],'generated',key='rfd4-ligand-pose-'+task,labels={'generated':'Generated','fixed':'Held as uploaded'})
    if task=='dna_binder':o['dna_conditioning']=_pills('DNA Conditioning',['sequence','structure'],'sequence',key='rfd4-dna-conditioning',labels={'sequence':'Sequence','structure':'Structure'})
    if task=='existing_structure' and not sequence_input:
        if o['operation']!='structure_prediction':
            o['fixed_chains']=st.multiselect('Chains To Keep Fixed' if o['operation']=='diversify' else 'Target Chains',chains,key='rfd4-fixed-chains')
        if o['operation']=='sidechain_packing':o['ligand_pose']=_pills('Ligand Pose',['generated','fixed'],'generated',key='rfd4-ligand-pose-'+task,labels={'generated':'Generated','fixed':'Held as uploaded'})
        if o['operation']=='diversify' and st.session_state.get('rfd4-noise-mode','fraction')=='fraction':o['noise_fraction']=st.number_input('Noise Fraction',.05,1.,.3)
        if o['operation'] in ('diversify','inverse_folding'):o['keep_residues']=_residues('Residues To Keep',path,chains,'rfd4-keep')
    p['num_designs']=1 if task=='score_structures' else st.number_input(
        'Predictions per input (Studio count)' if folding else 'Number of Designs (Studio count)',
        1, 500000, 1 if folding else 8, key='rfd4-count',
        help=(
            'Studio count of independent designs/predictions per input record. For RFD4 generation this maps to '
            'num_replicates; with multiple input records, each receives this many replicates. It is not RFD4’s '
            'internal GPU dataloader batch size.'
        ),
    )
    if task=='protein_binder':
        p['avoid_structures']=_source(owner_id,key='rfd4-avoid',label='Structures to Avoid',multiple=True)
        p['also_bind_structures']=_source(owner_id,key='rfd4-alsobind',label='Structures to Also Bind',multiple=True)
    if task=='idr_peptide_binder':
        with st.expander('View Optional Settings'):_optional('protein_binder',o,p,path,chains,folding,condition_spec=True)
    elif not preview_task:
        with st.expander('View Optional Settings'):_optional(task,o,p,path,chains,folding)
    if folding:o.update(refold=False,filter_conditions=False)
    with st.popover('Compute settings'):
        p['batch_size']=st.number_input(
            'RFD4 GPU Batch Size', 1, 16, 1, key='rfd4-batch-size',
            help=(
                'RFD4’s internal dataloader/inference batch size: how many designs it processes together on the GPU. '
                'This does not change Studio’s requested total; increasing it may improve GPU utilization but uses more '
                'GPU memory and can cause out-of-memory errors.'
            ),
        )
        p['cpu_cores']=st.number_input('CPU Cores',1,256,8,key='rfd4-cpu-cores')
        p['compile_model']=st.checkbox('Compile Model',key='rfd4-compile-model')
        prioritize=st.checkbox('Prioritize in Scheduler Queue',key='rfd4-prioritize')
    campaign_ready = True
    if not folding and task in {'protein_binder', 'idr_peptide_binder'}:
        with st.expander('Campaign objective (optional)'):
            from mn_molecular_design_studio.distributed_campaign_config_ui import render_campaign_config_editor
            campaign, campaign_ready = render_campaign_config_editor(
                st, 'rfd4', None, key='rfd4-native-campaign', seed=int(p.get('seed') or 0),
                generation_cpu_cores=p.get('cpu_cores'),
            )
            if campaign:
                p['distributed_campaign'] = campaign
    if not folding:
        st.button('Load example inputs',on_click=_load_example_inputs,args=('pd_l1_binder',owner_id),key='rfd4-load-example')
    validate_clicked, submit_column, _ = st.columns([1.25,1,1])
    placeholder_preview=task=='cyclic_peptide_binder'
    validate_clicked=validate_clicked.button('Validate / Prepare',disabled=placeholder_preview or not ready,key='rfd4-validate-prepare')
    clicked=submit_column.button('Submit',type='primary',disabled=placeholder_preview or not ready or not campaign_ready or not project_id,key='rfd4-submit')
    st.caption('Validate / Prepare checks the inputs in the RFD4 container. It does not create a job or run inference.')
    if validate_clicked:
        try:
            with st.spinner('Checking inputs with RFD4…'):
                from mn_molecular_design_studio.engine.workflow_engine.workflows.rfd4_proteina_prepare import validate_inputs
                report=validate_inputs(p)
            st.success(f"RFD4 input preparation passed for {report['count']} input(s). No job was created.")
            for item in report['items']:
                st.write(f"**{item['name']}** · {item['atoms']} atoms · designed chains: {', '.join(item['designed_chains']) or 'none'} · target chains: {', '.join(item['target_chains']) or 'none'}")
                for warning in item['warnings']:
                    st.warning(str(warning))
        except (ValueError,OSError,RuntimeError) as exc:
            st.error(str(exc))
    if clicked:
        try:
            inputs.validate_options(p)
            if folding:
                params={k:v for k,v in p.items() if k in {'target_structure','protein_sequences','num_designs','nsteps','seed','options','cpu_cores','job_name','batch_size','compile_model'}}
                submit(schema,'refolding.rfd4_proteina',params,project_id,{'prioritize':prioritize},title)
            else:submit(schema,'design.rfd4_proteina',p,project_id,{'prioritize':prioritize},title)
        except (ValueError,OSError) as exc:st.error(str(exc))
    if not project_id:st.caption('Select a project to submit.')


def render_folding_setup(schema,**kwargs):render_setup(schema,folding=True,**kwargs)


def render_results(run_id: str, run_dir: Path, candidates: list[dict], parameters: dict, result: dict) -> None:
    folding=(parameters.get('options') or {}).get('operation')=='structure_prediction' or bool(parameters.get('protein_sequences'))
    st.subheader("RFD4-Proteina structure predictions" if folding else "RFD4-Proteina designs")
    st.caption("RFDiffusion4-Proteina co-designs sequence and structure. Selected refolding checks and native metrics are recorded with each run. Continue with other prediction or analysis nodes in a pipeline.")
    if not candidates:
        st.info("No retained structures. Check the native metrics and run logs for filtering decisions or errors.")
        return
    import pandas as pd

    rows = []
    for candidate in candidates:
        metrics = candidate.get("metrics") or {}
        rows.append({
            "Design": candidate.get("candidate_id"),
            "Sequence": candidate.get("binder_sequence") or "—",
            "Designed Chains": ", ".join(candidate.get("binder_chains") or []),
            "Designed Residues": len(candidate.get("binder_sequence") or "") or "—",
            "Target Chains": ", ".join(candidate.get("target_chains") or []) or "—",
            "Total Residues": metrics.get("generated_residue_count"),
        })
    st.dataframe(pd.DataFrame(rows), hide_index=True, use_container_width=True)
    selected = st.selectbox("Design to view", range(len(candidates)), format_func=lambda i: candidates[i].get("candidate_id", f"Design {i + 1}"), key=f"rfd4-result-{run_id}")
    candidate = candidates[selected]
    path = Path(candidate.get("complex_pdb") or "")
    if path.is_file():
        try:
            molstar_custom_component(
                [StructureVisualization(pdb=path.read_text(errors="ignore"), color="molecule-type", representation_type="cartoon")],
                key=f"rfd4-view-{run_id}-{selected}", height=560, show_controls=True, show_structure_tools=False,
            )
            st.download_button("Download generated CIF", path.read_bytes(), file_name=path.name, key=f"rfd4-download-{run_id}-{selected}")
        except (OSError, UnicodeDecodeError):
            st.warning("The structure is available in the run artifacts but could not be displayed in the viewer.")
