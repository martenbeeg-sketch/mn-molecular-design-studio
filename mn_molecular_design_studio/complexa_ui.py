"""Proteina-Complexa task page with Studio projects, targets and native results."""
from __future__ import annotations

import copy
from html import escape
from pathlib import Path
import shutil

import pandas as pd
import streamlit as st

from mn_molecular_design_studio import complexa_input as inputs
from mn_molecular_design_studio.auth import user_input_root
from mn_molecular_design_studio.target_editor import render_target_source, render_target_editor, _target_widget_keys
from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer import StructureVisualization, molstar_custom_component
from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path

STYLE = '''<style>
.block-container{max-width:1240px!important;padding-top:3.25rem!important}
.st-key-complexa-surface{--cx-green:#27674e;--cx-mint:#edf6f0}
.st-key-complexa-surface .cx-heading{text-align:center;padding:.5rem 0 1.3rem}
.st-key-complexa-surface .cx-eyebrow{font-size:.68rem;letter-spacing:.12em;font-weight:650;color:#698575}
.st-key-complexa-surface h1{font-size:2.1rem;color:#173d31;letter-spacing:-.035em;margin:.6rem 0}
.st-key-complexa-surface .cx-heading p{font-size:.92rem;color:#718478;margin:.7rem auto;max-width:580px}
.st-key-complexa-surface .cx-version{color:#7c9284;font-size:.75rem}
.st-key-complexa-surface .cx-flow{display:flex;align-items:center;justify-content:center;gap:.6rem;flex-wrap:wrap;margin:1.1rem 0 1.4rem;font-size:.74rem;color:#67816e}
.st-key-complexa-surface .cx-flow span{background:#edf6f0;border:1px solid #d6e5db;border-radius:6px;padding:.5rem .75rem;color:#38624b}
.st-key-complexa-surface [data-testid="stVerticalBlockBorderWrapper"]{border-color:#d7e6dc!important;border-radius:12px!important}
[class*="st-key-complexa"][class*="-input-card"],[class*="st-key-complexa"][class*="-sampling-card"]{padding:1.2rem!important;background:#fff;border:1px solid #d7e6dc;border-radius:12px;margin:.5rem 0 .7rem}
[class*="st-key-complexa"][class*="-sampling-card"]{background:#f4f9f5}
.st-key-complexa-surface .cx-section{font-size:.72rem;font-weight:650;letter-spacing:.09em;color:#5b7b67;margin:.1rem 0 .8rem}
.st-key-complexa-surface [data-variant="pills"]{border:1px solid #bfd6c7!important;border-radius:7px!important;color:#38614e!important;background:#f6faf7;min-height:2.6rem;padding:.45rem .9rem!important;box-shadow:none!important}
.st-key-complexa-surface [data-variant="pills"][aria-checked="true"]{background:#27674e!important;color:white!important;border-color:#27674e!important}
.st-key-complexa-surface [data-baseweb="input"],.st-key-complexa-surface [data-baseweb="textarea"]{background:#f4f8f5!important;border-color:#d7e6dc!important;border-radius:7px!important}
.st-key-complexa-surface [data-testid="stExpander"]{border-color:#d7e6dc!important;border-radius:8px!important;background:#fcfefc}
.st-key-complexa-surface [data-testid="stFileUploaderDropzone"]{background:#f3f8f5!important;border:1px dashed #b9d2c2!important;border-radius:8px}
.st-key-complexa-surface [data-testid="stButton"] button{border-radius:7px!important;border-color:#bad2c3;color:#27674e}
.st-key-complexa-surface button[kind="primary"]{background:#27674e!important;border-color:#27674e!important;color:white!important;min-height:2.8rem}
.st-key-complexa-surface button[kind="primary"]:disabled{background:#9bb9a7!important;border-color:#9bb9a7!important}
.st-key-complexa-surface [data-baseweb="tag"]{background:#e0eee5!important;color:#27674e!important}
</style>'''


def example_parameters(task):
    lo,hi,count,strategy=inputs.DEFAULTS[task]
    p={"pipeline_task":task,"target_pdb":str(inputs.ASSETS/inputs.EXAMPLES[task]),"target_chains":[],"binder_length":f"{lo}-{hi}","num_designs":count,"n_steps":400,"replicas":2,"search_strategy":strategy}
    if task=="protein_binder":p.update(target_chains=["A"],hotspots="A37,A39,A49,A98")
    elif task=="ligand_binder":p["ligand"]="IAI"
    else:p.update(ligand="L:1",motif_atom_spec=inputs.MOTIF_EXAMPLE)
    return p


def job_parameters(parameters, jobs):
    result=[]
    if not 1<=jobs<=100:raise ValueError("Choose from 1 to 100 jobs.")
    for index in range(jobs):
        p=copy.deepcopy(parameters)
        if jobs>1:p["campaign_name"]=f"{p.get('campaign_name') or 'complexa'}-{index+1:03d}"
        if p.get("seed"):p["seed"]=(int(p["seed"])+index)%2147483647 or 1
        result.append(p)
    return result


def render_inputs(*, owner_id, key_prefix="complexa", generation_only=False, initial=None, connected_target=None, include_jobs=True):
    """The same task/structure/sampling form in a tool page or canvas node."""
    initial = copy.deepcopy(initial or {})
    st.html('<div class="cx-section">DESIGN TASK</div>')
    task = st.pills("Design task" if generation_only else "Pipeline", list(inputs.TASKS),
                    format_func=inputs.TASKS.get, selection_mode="single",
                    default=initial.get("pipeline_task", "protein_binder"),
                    key=key_prefix+"-task", label_visibility="collapsed") or "protein_binder"
    defaults = initial if initial.get("pipeline_task", task) == task else {}
    descriptions = {"protein_binder":"Design a protein binder against selected target chains and optional hotspots.",
                    "ligand_binder":"Design a protein binding pocket around one ligand from the input structure.",
                    "motif_ligand":"Preserve selected motif atoms around a ligand while designing the surrounding protein."}
    st.write(descriptions[task])
    prefix=key_prefix+"-"+task
    lo,hi,default_count,strategy=inputs.DEFAULTS[task]
    p={"pipeline_task":task,"target_chains":[]}
    error=None
    fallback={**(connected_target or {}), **{k:v for k,v in defaults.items() if v not in (None,"",[])}}
    with st.container(key=key_prefix+"-input-card"):
        if st.button("Load example inputs",icon=":material/science:",key=prefix+"-example"):
            example=example_parameters(task)
            source=Path(example["target_pdb"]);destination=user_input_root(owner_id)/"structures"/source.name
            destination.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,destination)
            st.session_state[prefix+"-fetched-target"]=str(destination)
            st.session_state[prefix+"-fetched-target-name"]=source.stem
            st.session_state.pop(prefix+"-target-upload",None)
            st.session_state[prefix+"-example-hotspots"]=example.get("hotspots","")
            st.session_state[_target_widget_keys(prefix,destination)["hotspots"]]=example.get("hotspots","")
            st.session_state[prefix+"-motif-spec"]=example.get("motif_atom_spec","")
            st.rerun()
        source=render_target_source(key_prefix=prefix,label="Target",example_code="5SDV" if task=="ligand_binder" else "1NZY" if task=="motif_ligand" else "7QF4",
                                    allow_cif=True,allow_non_protein=task!="protein_binder",owner_id=owner_id,upload_in_expander=True,
                                    upload_help="PDB or mmCIF. Include ligand coordinates for ligand-aware tasks.")
        if not source.path and fallback.get("target_pdb"):
            from mn_molecular_design_studio.target_editor import TargetSource
            path=resolve_stored_path(fallback["target_pdb"],must_exist=True)
            if path:source=TargetSource(path=path,name=path.stem,owner_id=owner_id)
        if source.path:
            p["target_pdb"]=str(source.path)
            try:
                if task=="protein_binder":
                    from mn_molecular_design_studio.target_editor import _hash
                    chains=fallback.get("target_chains",[])
                    if isinstance(chains,str):chains=[v.strip() for v in chains.split(",") if v.strip()]
                    if chains:st.session_state.setdefault(prefix+"-chains-"+_hash(str(source.path)),chains)
                    st.session_state.setdefault(_target_widget_keys(prefix,source.path)["hotspots"],fallback.get("hotspots",""))
                    target=render_target_editor(source,key_prefix=prefix,legacy_hotspots_key=prefix+"-example-hotspots",show_structure_tools=False)
                    p.update(target_pdb=str(target.design_path or source.path),target_chains=target.chains,hotspots=target.hotspots)
                    if not target.ready:error="Select target chains before submitting."
                else:
                    st.success(f"{source.name or source.path.name} loaded")
                    info=inputs.inventory(source.path);choices=inputs.ligand_choices(info)
                    if choices:
                        selected=defaults.get("ligand","")
                        selected=next((k for k,v in choices.items() if k==selected or v['name']==selected),next(iter(choices)))
                        p["ligand"]=st.selectbox("Ligand",list(choices),index=list(choices).index(selected),
                                                format_func=lambda value:choices[value]["name"]+" · residue "+choices[value]["label"],key=prefix+"-ligand-"+source.path.stem)
                    else:error="This structure contains no selectable ligand. Load an example or choose a structure containing a ligand."
                    if task=="motif_ligand":
                        st.session_state.setdefault(prefix+"-motif-spec",defaults.get("motif_atom_spec",""))
                        p["motif_atom_spec"]=st.text_area("Motif atom specification",placeholder="B64: [O, C]; B86: [CB, CA, N, C]",key=prefix+"-motif-spec",height=110,
                                                         help="Use residue IDs and atom names from your original structure. A prepared copy maps ligand and motif chains for the native engine.")
                        with st.expander("Available motif atoms"):
                            st.dataframe([{"Residue":label,"Name":row["name"],"Atoms":", ".join(row["atoms"])} for label,row in info["residues"].items() if row["kind"]=="protein"],hide_index=True,height=230)
                    with st.expander("Structure preview",expanded=False):
                        molstar_custom_component([StructureVisualization(pdb=source.path.read_text(),color="molecule-type",representation_type="cartoon+ball-and-stick")],key=prefix+"-preview-"+source.path.stem,height=430,show_structure_tools=False)
            except (ValueError,OSError,KeyError) as exc:error=str(exc)
        lengths=str(defaults.get("binder_length",f"{lo}-{hi}")).split("-",1)
        a,b=st.columns(2)
        minimum=int(a.number_input("Binder length (min)",1,3000,int(lengths[0]),key=prefix+"-length-min"))
        maximum=int(b.number_input("Binder length (max)",1,3000,int(lengths[-1]),key=prefix+"-length-max"))
        p["binder_length"]=f"{minimum}-{maximum}"
    with st.container(key=key_prefix+"-sampling-card"):
        st.html('<div class="cx-section">DESIGN SAMPLING</div>')
        a,b=st.columns(2)
        jobs=int(a.number_input("Number of jobs",1,100,1,key=prefix+"-jobs",help="Independent jobs, each scheduled on an available GPU.")) if include_jobs else 1
        p["num_designs"]=int(b.number_input("Designs per job" if include_jobs else "Number of designs",1,16,int(defaults.get("num_designs",default_count)),key=prefix+"-designs"))
        st.caption(f"{jobs*p['num_designs']} requested designs across {jobs} job(s). Removing duplicate sequences can reduce the result count." if generation_only else f"{jobs*p['num_designs']} requested designs across {jobs} job(s). Best of N can sample additional candidates per design.")
    with st.expander("Optional settings"):
        a,b=st.columns(2)
        with a:
            if not generation_only:
                modes=["single-pass"] if task=="motif_ligand" else ["best-of-n","single-pass"]
                p["search_strategy"]=st.selectbox("Search strategy",modes,format_func=lambda s:"Best of N" if s=="best-of-n" else "Single Pass",key=prefix+"-search")
                p["replicas"]=int(st.number_input("Candidates per design",1,8,2,disabled=p["search_strategy"]=="single-pass",key=prefix+"-replicas"))
            p["n_steps"]=int(st.number_input("Sampling steps",50,1000,int(defaults.get("n_steps",400)),key=prefix+"-steps"))
            p["batch_size"]=int(st.number_input("GPU batch size",1,16,int(defaults.get("batch_size",1)),key=prefix+"-batch"))
        with b:
            p["guidance_weight"]=st.number_input("Guidance weight",0.,1.,float(defaults.get("guidance_weight",1.)),.1,key=prefix+"-guidance")
            p["sampling_noise"]=st.number_input("Sampling noise",0.,1.,float(defaults.get("sampling_noise",.1)),.05,key=prefix+"-noise")
            p["deduplicate_sequences"]=st.checkbox("Remove duplicate sequences",bool(defaults.get("deduplicate_sequences",True)),key=prefix+"-dedup")
            p["seed"]=int(st.number_input("Random seed (0 = random)",0,2147483647,int(defaults.get("seed",0)),key=prefix+"-seed"))
        if generation_only:st.caption("Single Pass sampling generates backbone and sequence. Choose sequence redesign and structure prediction steps on the pipeline canvas.")
        elif task=="motif_ligand":st.caption("Motif + Ligand uses Single Pass because its native pipeline has no search reward configured.")
    if not error:
        try:inputs.prepare_request(p | ({"generator_only":True,"search_strategy":"single-pass","replicas":1} if generation_only else {}))
        except (ValueError,TypeError,OSError) as exc:error=str(exc)
    return p,jobs,error


def render_setup(schema, *, owner_id, project_selector, submit):
    generation_only=schema["workflow"]=="design.proteina_complexa_generate"
    prefix="complexa-generator" if generation_only else "complexa"
    title="Proteina-Complexa Generator" if generation_only else "Proteina-Complexa"
    st.html(STYLE)
    with st.container(key="complexa-surface"):
        if st.button("← Back to All tools",key=prefix+"-back"):
            st.session_state["page"]="All tools";st.session_state.pop("selected_workflow",None);st.rerun()
        st.html('<div class="cx-heading"><div class="cx-eyebrow">'+("DESIGN ENGINE" if generation_only else "DESIGN PIPELINE")+' · ALL-ATOM BINDER DESIGN</div><h1>'+title+'</h1>'
                '<p>Design protein binders, ligand-binding proteins, and scaffolds around functional motifs.</p>'
                f'<span class="cx-version">{inputs.CODE_VERSION}</span></div>')
        steps=["Generate backbone + sequence","Choose downstream tools on canvas"] if generation_only else ["Generate","Filter","Sequence design & refold","Analyze"]
        st.html('<div class="cx-flow">'+"→".join('<span>'+s+'</span>' for s in steps)+'</div>')
        with st.expander("About this workflow"):
            st.write("Generate structures and sequences with the native design engine. Generated designs can feed sequence design and structure prediction nodes in Pipelines." if generation_only else "Each job runs the installed engine’s complete native pipeline: generation, filtering, refolding, task-specific success analysis, and sequence redesign where configured.")
            st.markdown("[Source code](https://github.com/NVIDIA-BioNeMo/Proteina-Complexa) · [Tamarind input reference](https://app.tamarind.bio/tools/proteina-complexa)")
            if not generation_only:st.caption("Protein binders use AF2 evaluation; ligand-aware pipelines use RF3.")
        a,b=st.columns(2)
        name=a.text_input("Job name","proteina-complexa-generator" if generation_only else "proteina-complexa",key=prefix+"-job-name").strip()
        with b:project=project_selector(prefix)
        p,jobs,error=render_inputs(owner_id=owner_id,key_prefix=prefix,generation_only=generation_only)
        p["campaign_name"]=name or prefix
        with st.expander("Scheduler & local readiness"):
            a,b=st.columns(2)
            gpu=a.selectbox("GPU assignment",["Automatic","GPU 0","GPU 1"],key=prefix+"-gpu")
            p["gpu_device"]={"Automatic":"all","GPU 0":"0","GPU 1":"1"}[gpu]
            p["cpu_cores"]=int(b.number_input("CPU cores per job",1,256,4,key=prefix+"-cpus"))
            prioritize=st.checkbox("Prioritize in scheduler queue",key=prefix+"-priority")
            readiness=inputs.readiness(p["pipeline_task"],generator_only=generation_only)
            if readiness["ready"]:st.success("Generator and autoencoder checkpoints are available." if generation_only else "Reference checkpoints and native folding/sequence-design weights are available.")
            else:st.warning("Missing reference files: "+", ".join(readiness["missing"]))
        if error:st.caption(error)
        if st.button(f"Submit {jobs} job"+("s" if jobs!=1 else ""),type="primary",disabled=not project or bool(error) or not readiness["ready"],use_container_width=True,key=prefix+"-submit"):
            try:
                for parameters in job_parameters(p,jobs):submit(schema,schema["workflow"],parameters,project,{"prioritize":prioritize},title)
            except (ValueError,TypeError,OSError) as exc:st.error(str(exc))
        if not project:st.caption("Select a project to submit.")


def render_results(run_id,run_dir,candidates,parameters,result):
    generation_only=bool(parameters.get('generator_only'))
    st.subheader(inputs.TASKS.get(parameters.get("pipeline_task"),"Protein Binder")+(" · generated designs" if generation_only else " · pipeline results"))
    native=[c for c in candidates if c.get("metrics",{}).get("candidate_pool_level")=="native_pipeline" or c.get("metrics",{}).get("result_kind")=="native_pipeline"]
    passed=[c for c in native if c.get("metrics",{}).get("pass_filters") is True]
    if generation_only:
        st.metric('Generated designs',len(candidates))
        st.caption('Backbone and sequence are ready for sequence redesign or prediction. Confidence scores and filter passes become available after evaluation.')
        from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json
        marker=read_json(Path(run_dir)/'molecular_design_studio.json')
        if marker.get('project_id') and result.get('success') is True and candidates:
            if st.button('Continue on pipeline canvas',key='complexa-continue-'+run_id,type='primary'):
                from mn_molecular_design_studio.pipeline_graph import new_graph
                graph=new_graph('complexa_generated')
                graph['nodes'][0]['data']['params'].update(source_run_dir=str(run_dir),selected_candidate_ids=[c['candidate_id'] for c in candidates])
                project=marker['project_id']
                st.session_state['pipeline-pending-draft']={'project_id':project,'graph':graph}
                st.session_state.setdefault('project-pending-selection',{})['visual-pipeline-project-id']=project
                st.session_state['page']='Pipelines'
                st.rerun()
    else:
        a,b,c=st.columns(3);a.metric("Evaluated designs",len(native));b.metric("Filter-passing designs",len(passed));c.metric("Generated structures",len(candidates)-len(native))
        st.caption("Native pipeline and generated structures are kept separately. Only explicit filter-passing designs are offered to project design-validation batches.")
    rows=native or candidates
    if not rows:
        st.info("Designs will appear after generation finishes." if generation_only else "Designs and evaluation metrics will appear after the pipeline finishes.");return
    def percentage(value):
        return value*100 if isinstance(value,(int,float)) and 0<=value<=1 else value
    table=[{"Design":c['candidate_id'],'Residues':len(c.get('binder_sequence') or ''),'Sequence':c.get('binder_sequence')} for c in rows] if generation_only else [{"Design":c["candidate_id"],"Filter pass":c.get("metrics",{}).get("pass_filters"),"Complex pLDDT (%)":percentage(c.get("metrics",{}).get("binder_plddt")),"ipTM":c.get("metrics",{}).get("iptm"),"Binder RMSD (Å)":c.get("metrics",{}).get("binder_rmsd"),"Sequence":c.get("binder_sequence")} for c in rows]
    st.dataframe(table,hide_index=True,use_container_width=True)
    chosen=st.selectbox("Design",range(len(rows)),format_func=lambda i:rows[i]["candidate_id"],key="complexa-result-"+run_id)
    candidate=rows[chosen];path=resolve_stored_path(candidate.get("complex_pdb"),run_dir=Path(run_dir),must_exist=True)
    if (warning:=candidate.get('raw_metadata',{}).get('motif_mapping_warning')):st.warning(warning)
    if path:
        molstar_custom_component([StructureVisualization(pdb=path.read_text(),color="molecule-type",representation_type="cartoon+ball-and-stick")],key="complexa-result-view-"+run_id+"-"+str(chosen),height=530,show_structure_tools=False)
        st.download_button("Download selected structure",path.read_bytes(),file_name=path.name,key="complexa-result-download-"+run_id)
    with st.expander("Generation details" if generation_only else "Native metrics and result tables"):
        st.json(candidate.get("metrics",{}))
        if generation_only:st.json(candidate.get('raw_metadata',{}))
        for index,csv in enumerate(sorted((Path(run_dir)/"artifacts/raw/proteina_complexa/evaluation_results").glob("**/RAW_*combined.csv"))):
            st.download_button("Download "+csv.name,csv.read_bytes(),file_name=csv.name,key=f"complexa-csv-{run_id}-{index}")
