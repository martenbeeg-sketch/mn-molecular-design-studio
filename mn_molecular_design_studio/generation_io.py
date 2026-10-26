"""Generation outputs passed to independently selected design and folding tools."""
from pathlib import Path

from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json
from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates


def generated_candidates(run_dir, *, stage=None):
    """Unvalidated designs have their own source; successful-design filters stay strict."""
    run_dir=Path(run_dir)
    metadata=read_json(run_dir/'metadata.json')
    params=read_json(run_dir/'input.json').get('params',{})
    from mn_molecular_design_studio.generator_input import TOOLS
    supported={
        'kejia_idr',
        'peptide_trap',
        'proteina_complexa_generator',
        'bindcraft2_generator',
        'boltzgen_generator',
        *TOOLS.values(),
    }
    generator_mode = bool(
        params.get('generator_only')
        or params.get('generator_only_tile')
        or params.get('trajectory_only')
        # Earlier Kejia jobs predate the shared generator flag.
        or metadata.get('tool') == 'kejia_idr'
    )
    if (metadata.get('status') not in {'completed','succeeded'}
            or read_json(run_dir/'result.json').get('success') is not True
            or not generator_mode
            or metadata.get('tool') not in supported):
        return []
    return [c for c in read_candidates(run_dir)
            if c.get('stage') in ({stage} if stage else {'generation.backbone_sequence','generation.backbone'})
            and (c.get('metrics',{}).get('result_kind')=='generated'
                 or metadata.get('tool')=='kejia_idr' and c.get('metrics',{}).get('result_kind')=='generation_only')]


def has_ligand_target(candidate):
    raw=candidate.get('raw_metadata') or {}
    return raw.get('target_entity_kind')=='ligand'


def _connections(line):
    try:
        return [int(line[i:i+5]) for i in range(6,len(line),5) if line[i:i+5].strip()]
    except ValueError:
        return []


def restore_ligand_context(source, destination, chain_map):
    """Keep exact ligand atoms, coordinates and CONECT bonds through MPNN I/O."""
    source_lines=Path(source).read_text().splitlines()
    ligand_lines=[line for line in source_lines if line.startswith('HETATM') and line[21:22] in chain_map]
    if not ligand_lines:
        raise ValueError('The ligand target is missing from the source structure.')
    target_chains=set(chain_map.values())
    lines=[line for line in Path(destination).read_text().splitlines()
           if line[:6].strip() not in {'END','ENDMDL','MODEL','CONECT'}
           and not (line.startswith('HETATM') and line[21:22] in target_chains)]
    serial=max((int(line[6:11]) for line in lines if line.startswith(('ATOM  ','HETATM'))),default=0)+1
    serial_map={}
    for line in ligand_lines:
        serial_map[int(line[6:11])]=serial
        lines.append(line[:6]+f'{serial:5d}'+line[11:21]+chain_map[line[21:22]]+line[22:])
        serial+=1
    for line in source_lines:
        if not line.startswith('CONECT'):continue
        ids=_connections(line)
        if ids and ids[0] in serial_map:
            neighbors=[serial_map[v] for v in ids[1:] if v in serial_map]
            if neighbors:lines.append('CONECT'+''.join(f'{v:5d}' for v in [serial_map[ids[0]],*neighbors]))
    Path(destination).write_text('\n'.join([*lines,'END','']))


def extract_chain_with_bonds(source, destination, chain):
    lines=Path(source).read_text().splitlines()
    atoms=[line for line in lines if line.startswith(('ATOM  ','HETATM')) and line[21:22]==chain]
    if not atoms:return None
    serials={int(line[6:11]) for line in atoms}
    bonds=[]
    for line in lines:
        if not line.startswith('CONECT'):continue
        ids=_connections(line)
        if ids and ids[0] in serials:
            neighbors=[v for v in ids[1:] if v in serials]
            if neighbors:bonds.append('CONECT'+''.join(f'{v:5d}' for v in [ids[0],*neighbors]))
    Path(destination).write_text('\n'.join([*atoms,*bonds,'END','']))
    return Path(destination)
