"""Adapt Studio prediction results to RFD4's native refolding protocol.

Confidence fields are transferred only when the selected engine produced them.
Missing PAE or incompatible chain/atom inventories fail explicitly.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from Bio.PDB import MMCIFIO
from mn_molecular_design_studio import cofolding_results, molecule_imports
from mn_molecular_design_studio.design_validation import read_structure
from mn_molecular_design_studio.engine.workflow_engine.workflows.graph_pipeline import _optimized_prediction_request
from mn_molecular_design_studio.engine.workflow_engine.workflows.optimized_prediction import run_optimized_prediction


def fold_entry(entry, bridge):
    root=bridge.run_dir/'artifacts/rfd4_proteina/pipeline_oracles'/entry['name']
    root.mkdir(parents=True,exist_ok=True)
    source=root/'input.cif'; source.write_text(entry['cif'])
    components=molecule_imports.parse_structure(source.read_bytes(),source.name,include_ligands=True)[0]['components']
    if not components: raise ValueError('No supported molecules in refolding input.')
    if any(row['type']=='small_molecule' and not row.get('ccd') for row in components):
        raise ValueError('This refolding adapter needs a known CCD ligand. Use OpenFold3 JAX for custom CIF chemistry.')
    molecules={key:':'.join(row['sequence'] for row in components if row['type']==kind) for kind,key in [('protein','target_protein_sequences'),('dna','dna_sequences'),('rna','rna_sequences')]}
    molecules['ligand_ccd_codes']=':'.join(row['ccd'] for row in components if row['type']=='small_molecule')
    task_type='protein_ligand' if molecules['ligand_ccd_codes'] else 'protein_dna' if molecules['dna_sequences'] else 'protein_rna' if molecules['rna_sequences'] else 'protein'
    molecules['task_type']=task_type
    p={'samples':1,'seed':entry.get('seed',0),'num_recycles':bridge.options.get('folding_recycles',5),
       'target_use_msa':False,'use_target_template':False,'accept_alphafold3_terms':bridge.options.get('accept_alphafold3_terms',False)}
    request=_optimized_prediction_request(bridge.engine,p,{'molecule_task':molecules},'rfd4-'+entry['name'])
    request['parameters']['cpu_cores']=bridge.cpu_cores
    request['parameters']['gpu_device']=bridge.gpu_device
    child=run_optimized_prediction(**request['parameters'])
    with (root/'prediction_run.json').open('w') as out: json.dump({'run_dir':str(child),'request':request},out,indent=2)
    result=json.loads((child/'result.json').read_text())
    if not result.get('success',result.get('status')=='completed'):
        # Some workflow results keep success in the status file, so rely on outputs too.
        if not result.get('outputs',{}).get('candidates'): raise RuntimeError(f'Refolding failed; see {child}')
    candidates=result.get('outputs',{}).get('candidates',[])
    if not candidates: raise RuntimeError(f'Selected refolding engine returned no structure; see {child}')
    structure,files=cofolding_results.candidate_files(child,candidates[0])
    if structure is None: raise RuntimeError('The refolding structure is missing.')
    confidence=cofolding_results.read_confidence(structure,tuple(files))
    model=next(read_structure(structure).get_models())
    chains=list(model.get_chains())
    # Predictors receive chains grouped by molecule type through the shared adapter.
    ordered=[row for kind in ('protein','dna','rna','small_molecule') for row in components if row['type']==kind]
    if len(chains)!=len(ordered): raise ValueError('Refolding changed the number of chains; cannot safely map coordinates or confidence.')
    original_ids=[row['chain'] for row in ordered]
    if len(set(original_ids))!=len(original_ids): raise ValueError('The refolding input has ambiguous chain IDs.')
    for chain,identity in zip(chains,original_ids): chain.id='STUDIO_'+identity
    for chain,identity in zip(chains,original_ids): chain.id=identity
    predicted=[atom for atom in model.get_atoms()]
    values=confidence['plddt']
    if not values.size: raise ValueError('The selected engine produced no pLDDT confidence.')
    residues=list(model.get_residues())
    if len(values)==len(predicted):
        for atom,value in zip(predicted,values): atom.bfactor=float(value)
    elif len(values)==len(residues):
        for residue,value in zip(residues,values):
            for atom in residue: atom.bfactor=float(value)
    else: raise ValueError('Confidence array cannot be mapped to refolded atoms/residues.')
    summary={'chain_ids':original_ids,'overall_plddt':float(values.mean()/100),
             'chain_plddt':[float(np.mean([a.bfactor for a in c.get_atoms()])/100) for c in chains]}
    for key in ('ptm','iptm','ranking_score'):
        if key in confidence['metrics']: summary[key]=float(confidence['metrics'][key])
    details=confidence['details']
    for key in ('chain_ptm','chain_pair_iptm','chain_pair_pae_min'):
        if key in details: summary[key]=np.asarray(details[key]).tolist()
    pae=confidence['pae']
    if pae is None: raise ValueError('Selected engine did not produce PAE needed for the RFD4 refolding checks.')
    sizes=[len(list(c.get_residues())) for c in chains]
    if sum(sizes)!=len(pae):
        raise ValueError('PAE tokens cannot be mapped safely to these molecule chains. Use the native OpenFold3 oracle for this input.')
    boundaries=np.cumsum([0,*sizes]); slices=[slice(a,b) for a,b in zip(boundaries[:-1],boundaries[1:])]
    summary.update(overall_pae=float(pae.mean()),chain_pae=[float(pae[s,s].mean()) for s in slices],
                   chain_pair_pae=[[float(pae[a,b].mean()) for b in slices] for a in slices])
    summary.setdefault('chain_pair_pae_min',[[float(pae[a,b].min()) for b in slices] for a in slices])
    # Do not replace a missing per-chain pTM with the whole-complex pTM.
    out=(bridge.run_dir/Path(entry['output_dir']).relative_to('/work')).resolve(); out.relative_to(bridge.run_dir)
    out.mkdir(parents=True,exist_ok=True)
    writer=MMCIFIO(); writer.set_structure(model); writer.save(str(out/(entry['name']+'_model.cif')))
    (out/(entry['name']+'_summary_confidences.json')).write_text(json.dumps(summary,indent=2))
