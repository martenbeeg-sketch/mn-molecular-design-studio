"""CPU preparation and native CLI entrypoint, copied into the existing image per job.

This file is also invokable with --prepare-only to validate all task configs
without loading models or starting GPU inference.
"""
import json
import os
from pathlib import Path
import subprocess
import sys


def configure_environment(environment):
    os.environ.update(environment)
    os.environ['PATH']='/workspace/.venv/bin:/studio-tools:'+os.environ.get('PATH','')
    # The native validator expects a local dotenv file. It contains only the
    # explicit runtime paths, never credentials from the Studio environment.
    Path('.env').write_text(''.join(f'{key}={json.dumps(value)}\n' for key,value in environment.items()))


def ordered_motifs(spec):
    motifs=spec.get('motif_atoms') or {}
    order=spec.get('motif_residue_order',list(motifs))
    if len(order)!=len(set(order)) or set(order)!=set(motifs):
        raise ValueError('The motif residue order does not match its atom selections.')
    return order


def prepare(spec):
    from hydra import compose, initialize_config_dir
    from omegaconf import OmegaConf
    root=Path('/work/artifacts/raw/proteina_complexa')
    if spec['pipeline_task']!='protein_binder':
        import numpy as np
        from atomworks.io.utils.io_utils import load_any
        from biotite.structure.io.pdb import PDBFile
        atoms=load_any(str(root/'input'/spec['source_file']))[0]
        if not hasattr(atoms,'chain_id'):
            raise ValueError('Input structure did not produce an atom array.')
        ligand=spec['ligand_residue']
        ligand_label=ligand['label']
        ligand_mask=(atoms.chain_id==ligand_label[0])&(atoms.res_id==int(ligand_label[1:]))&(atoms.res_name==ligand['name'])
        motifs=spec.get('motif_atoms') or {}
        mask=ligand_mask.copy()
        for label in motifs:
            mask|=(atoms.chain_id==label[0])&(atoms.res_id==int(label[1:]))&~atoms.hetero
        atoms=atoms[mask & ~np.isin(atoms.element,['H','D'])].copy()
        ligand_mask=(atoms.chain_id==ligand_label[0])&(atoms.res_id==int(ligand_label[1:]))&(atoms.res_name==ligand['name'])
        if not ligand_mask.any():raise ValueError('Selected ligand is missing from the prepared structure.')
        # Calculate every mask before changing any chain/residue identifiers.
        order=ordered_motifs(spec)
        motif_masks=[(atoms.chain_id==label[0])&(atoms.res_id==int(label[1:]))&~atoms.hetero for label in order]
        for i,mask in enumerate(motif_masks,1):
            if not mask.any():raise ValueError('Selected motif residue is missing.')
            atoms.chain_id[mask]='B';atoms.res_id[mask]=i
        atoms.chain_id[ligand_mask]='A';atoms.res_id[ligand_mask]=1;atoms.res_name[ligand_mask]='L:0'
        pdb=PDBFile();pdb.set_structure(atoms);pdb.write(str(root/'input/target.pdb'))
    with initialize_config_dir(version_base=None,config_dir='/workspace/protein-foundation-models/configs'):
        cfg=compose(config_name=spec['base_config'])
    OmegaConf.set_struct(cfg,False)
    cfg=OmegaConf.merge(cfg,OmegaConf.create(spec['overrides']))
    # Resolve now so no installed task defaults leak into this run's input.
    OmegaConf.resolve(cfg)
    config=root/'studio_complexa.yaml'
    OmegaConf.save(cfg,config)
    return config


def export_generation_metadata(spec, config):
    """Save native motif matches for later sequence design, without refolding."""
    root=Path('/work/artifacts/raw/proteina_complexa')
    matcher=None
    error=None
    if spec['pipeline_task']=='motif_ligand':
        try:
            from omegaconf import OmegaConf
            from proteinfoundation.evaluation.motif_binder_eval import _load_motif_info, get_motif_fixed_residues
            from proteinfoundation.evaluation.motif_binder_eval_utils import _parse_task_config
            cfg=OmegaConf.load(config)
            task=cfg.generation.motif_target_dict_cfg.MN_APP_TARGET
            info=_load_motif_info(_parse_task_config('MN_APP_TARGET',task),task)
            matcher=lambda path:get_motif_fixed_residues(info,str(path),'B')
        except Exception as exc:
            error=str(exc)
    for path in sorted((root/'inference').glob('**/job_*/*.pdb')):
        if path.name.endswith(('_binder.pdb','_updated.pdb')):
            continue
        metadata={'pipeline_task':spec['pipeline_task'],
                  'target_entity_kind':'protein' if spec['pipeline_task']=='protein_binder' else 'ligand',
                  'ligand_original_code':(spec.get('ligand_residue') or {}).get('name'),
                  'fixed_motif_required':spec['pipeline_task']=='motif_ligand'}
        if metadata['fixed_motif_required']:
            fixed=matcher(path) if matcher else None
            def residue_names(pdb):
                return {line[21:22]+line[22:26].strip():line[17:20].strip()
                        for line in Path(pdb).read_text().splitlines() if line.startswith('ATOM  ')}
            expected=residue_names(root/'input/target.pdb')
            observed=residue_names(path)
            count=len(spec['motif_atoms'])
            verified=(bool(fixed) and len(fixed)==len(set(fixed))==count
                      and all(expected.get(f'B{i}') and expected[f'B{i}']==observed.get(label)
                              for i,label in enumerate(fixed,1)))
            metadata.update(fixed_residues=fixed or [],
                            motif_mapping_verified=verified)
            if not metadata['motif_mapping_verified']:
                metadata['motif_mapping_warning']=error or 'The native matcher could not identify every motif residue. Sequence redesign requires a verified motif map.'
        path.with_suffix('.studio.json').write_text(json.dumps(metadata,indent=2))


def main():
    spec=json.loads(Path(sys.argv[1]).read_text())
    os.chdir('/work/artifacts/raw/proteina_complexa')
    configure_environment(spec['environment'])
    community=Path('community_models')
    if not community.exists():community.symlink_to('/workspace/protein-foundation-models/community_models',target_is_directory=True)
    config=prepare(spec)
    if '--check-features' in sys.argv:
        from hydra.utils import instantiate
        from omegaconf import OmegaConf
        cfg=OmegaConf.load(config)
        for feature in instantiate(cfg.generation.dataloader.dataset.conditional_features):
            feature.setup([spec['length_range'][0]])
        print('Validated conditional input features:',spec['pipeline_task'])
        return 0
    if '--validate-only' in sys.argv:
        mode='generate' if spec.get('generator_only') else 'design'
        return subprocess.run(['/workspace/.venv/bin/complexa','validate',mode,str(config)],check=False).returncode
    if '--prepare-only' in sys.argv:
        print('Prepared configuration:',config)
        return 0
    mode='generate' if spec.get('generator_only') else 'design'
    result=subprocess.run(['/workspace/.venv/bin/complexa',mode,str(config)],check=False).returncode
    if result==0 and spec.get('generator_only'):
        export_generation_metadata(spec,config)
    return result


if __name__=='__main__':
    raise SystemExit(main())
