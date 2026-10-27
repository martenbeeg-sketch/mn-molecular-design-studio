"""Identity checks at scientific workflow boundaries; no confidence thresholds."""
from __future__ import annotations

import copy
from pathlib import Path

from mn_molecular_design_studio import design_validation as structures

CHAIN_IDS = structures.CHAIN_IDS


def structure_chemistry(path):
    """Record unsupported chemistry instead of flattening it into standard residues."""
    import gemmi
    from Bio.Data.PDBData import protein_letters_3to1
    structure = gemmi.read_structure(str(path))
    issues = []
    for chain in structure[0]:
        for residue in chain:
            info = gemmi.find_tabulated_residue(residue.name)
            if info.is_nucleic_acid() and residue.name not in {'DA', 'DC', 'DG', 'DT', 'A', 'C', 'G', 'U'}:
                issues.append(f'modified/unsupported nucleic acid residue {chain.name}:{residue.seqid}:{residue.name}')
            if info.is_amino_acid() and residue.name not in protein_letters_3to1:
                issues.append(f"modified/unsupported protein residue {chain.name}:{residue.seqid}:{residue.name}")
    for connection in structure.connections:
        # Canonical disulfides are reference information, not prediction restraints.
        if connection.type != gemmi.ConnectionType.Disulf:
            issues.append(f"unsupported explicit inter-residue bond {connection.name}: {connection.partner1} -> {connection.partner2}")
    if Path(path).suffix.lower() == '.pdb':
        # PDB CONECT cannot distinguish custom topology from a canonical CCD graph.
        if any(line.startswith('CONECT') for line in Path(path).read_text().splitlines()):
            issues.append('explicit PDB CONECT topology requires an atom/bond-aware handoff')
    return issues


def component_prediction_task(candidate, run_dir, *, target_use_msa=True):
    raw = candidate['raw_metadata']
    components = copy.deepcopy(raw['components'])
    kinds = {'protein', 'dna', 'rna', 'small_molecule'}
    if any(row.get('type') not in kinds for row in components):
        raise ValueError('Unsupported molecular component in prediction handoff.')
    source_ids = [row['chain'] for row in components]
    if len(set(source_ids)) != len(source_ids):
        raise ValueError('Ambiguous duplicate component chain identities.')
    targets, binders = list(candidate.get('target_chains') or []), list(candidate.get('binder_chains') or [])
    proteins = {row['chain']: row for row in components if row['type'] == 'protein'}
    if set(targets) & set(binders) or set(targets + binders) - proteins.keys():
        raise ValueError('Target/binder roles do not match protein components.')
    ordered_proteins = [proteins[key] for key in targets + binders]
    ordered_proteins += [row for key, row in proteins.items() if key not in targets + binders]
    ordered = ordered_proteins + [row for kind in ('dna', 'rna') for row in components if row['type'] == kind]
    ordered += [row for row in components if row['type'] == 'small_molecule' and not row.get('smiles')]
    ordered += [row for row in components if row['type'] == 'small_molecule' and row.get('smiles')]
    if len(ordered) > len(CHAIN_IDS):
        raise ValueError('Too many molecular components for prediction chain mapping.')
    chain_map = {row['chain']: CHAIN_IDS[i] for i, row in enumerate(ordered)}
    path_value = candidate.get('complex_pdb') or candidate.get('binder_pdb')
    path = None
    if candidate.get('source_tool') == 'rfd4_proteina' and not path_value:
        raise ValueError('RFD4 candidate structure is missing; sequence correspondence cannot be verified.')
    if path_value:
        path = Path(path_value)
        if not path.is_absolute():
            path = Path(run_dir) / path
        issues = structure_chemistry(path) + list(raw.get('unsupported_chemistry') or [])
        if raw.get('bonds') or raw.get('covalent_bonds'):
            issues.append('explicit component bonds are unsupported by this handoff')
        if issues:
            raise ValueError('Unsupported chemistry: ' + '; '.join(issues))
        actual = structures.protein_chains(structures.read_structure(path))
        if set(actual) != set(proteins):
            raise ValueError('Protein component mapping does not match structure chains.')
        for key, row in proteins.items():
            if structures.sequence(actual[key]) != row['sequence']:
                raise ValueError(f'Component sequence differs from structure chain {key}.')
            if raw.get('sequences_by_chain', {}).get(key, row['sequence']) != row['sequence']:
                raise ValueError(f'Generator sequence differs from component chain {key}.')
    binder_sequence = ':'.join(proteins[key]['sequence'] for key in binders)
    declared = str(candidate.get('binder_sequence') or '').replace(':', '')
    if binders and declared and declared != binder_sequence.replace(':', ''):
        raise ValueError('Declared binder sequence differs from molecular components.')
    ligands = [row for row in ordered if row['type'] == 'small_molecule']
    for row in ligands:
        if not (row.get('ccd') or row.get('smiles')):
            raise ValueError('Unsupported ligand without a CCD or SMILES identity.')
    contract = {
        'schema_version': 1, 'source_candidate_id': candidate['candidate_id'],
        'source_molecule_id': raw.get('studio_molecule_id'),
        'chain_map': chain_map, 'components': ordered,
        'target_protein_ids': [chain_map[key] for key in targets],
        'binder_protein_ids': [chain_map[key] for key in binders],
        'binder_sequence': binder_sequence,
        'source_structure': str(path) if path else '',
        # Reference constraints are retained, never imposed as predictor restraints.
        'reference_metadata': {key: copy.deepcopy(value) for key, value in raw.items()
                               if any(token in key for token in ('fixed', 'constraint', 'condition', 'motif', 'cdr', 'framework'))},
    }
    policies = raw.get('molecule_chain_options') or {}
    if isinstance(policies, list):
        policies = {item['id']: item for item in policies}
    chain_options = []
    for row in ordered_proteins:
        option = copy.deepcopy(policies.get(row['chain'], {}))
        option.update(id=chain_map[row['chain']])
        option.setdefault('use_msa', bool(target_use_msa and row['chain'] in targets))
        chain_options.append(option)
    task = {
        'name': candidate['candidate_id'], 'source_candidate_id': candidate['candidate_id'],
        'source_run_dir': str(run_dir), 'scientific_contract': contract,
        'protein_sequences': ':'.join(row['sequence'] for row in ordered_proteins),
        'dna_sequences': ':'.join(row['sequence'] for row in ordered if row['type'] == 'dna'),
        'rna_sequences': ':'.join(row['sequence'] for row in ordered if row['type'] == 'rna'),
        'ligand_ccd_codes': ':'.join(row['ccd'] for row in ligands if not row.get('smiles')),
        'ligand_smiles': '\n'.join(row['smiles'] for row in ligands if row.get('smiles')),
        'molecule_chain_options': chain_options,
    }
    if ligands:
        identity = raw.get('metadata') or {}
        task.update(compound_id=str(identity.get('compound_id') or identity.get('pubchem_cid') or candidate['candidate_id']),
                    compound_name=str(identity.get('compound_name') or raw.get('name') or candidate['candidate_id']))
    return task


def verify_prediction_structure(path, task):
    """Reject chain/sequence corruption before publishing prediction candidates."""
    structure = structures.read_structure(Path(path))
    actual = {key: structures.sequence(value) for key, value in structures.protein_chains(structure).items()}
    expected = {row['id']: row['sequence'] for row in task['proteins']}
    if (expected or task.get("scientific_contract")) and actual != expected:
        raise ValueError(f'Prediction sequence/chain mismatch: expected {expected}, observed {actual}.')
    contract = task.get('scientific_contract') or {}
    if contract:
        # Validate all polymer identities and CCD residue identities in the mapped chains.
        import gemmi
        model = gemmi.read_structure(str(path))[0]
        if {chain.name for chain in model} != set(contract['chain_map'].values()):
            raise ValueError('Prediction molecular component chain mapping differs from source.')
        for row in contract['components']:
            chain_id = contract['chain_map'][row['chain']]
            residues = [res for chain in model if chain.name == chain_id for res in chain]
            if not residues:
                raise ValueError(f'Prediction lost molecular component {chain_id}.')
            if row['type'] == 'small_molecule':
                if row.get('heavy_atoms'):
                    atoms = {atom.name: atom.element.name.upper() for res in residues for atom in res if not atom.is_hydrogen()}
                    if atoms != row['heavy_atoms']:
                        raise ValueError(f'Prediction changed ligand atom identities on {chain_id}.')
                if row.get('ccd') and {res.name for res in residues} != {row['ccd']}:
                    raise ValueError(f'Prediction changed CCD ligand identity on {chain_id}.')
            elif row['type'] in {'dna', 'rna'}:
                sequence = ''.join(res.name[-1] for res in residues)
                if sequence != row['sequence']:
                    raise ValueError(f'Prediction changed nucleic acid sequence on {chain_id}.')
    return {'sequence_chain_identity': 'pass' if expected or contract else 'not_checked',
            'source_candidate_id': task.get('source_candidate_id', '')}


def validate_component_chemistry(task, reference_root):
    """CCD conversion is supported only when its source atom identities match the installed CCD."""
    import pickle
    import gemmi
    contract = task.get('scientific_contract') or {}
    if not contract or not contract.get('source_structure'):
        return
    source = gemmi.read_structure(contract['source_structure'])[0]
    for row in contract['components']:
        if row['type'] != 'small_molecule':
            continue
        if row.get('smiles'):
            raise ValueError('Unsupported SMILES structure handoff without a verified source atom/bond mapping.')
        code = str(row['ccd'])
        if not code.isalnum():
            raise ValueError(f'Unsupported CCD identity {code!r}.')
        molecule_path = Path(reference_root) / 'boltz_models/mols' / (code + '.pkl')
        if not molecule_path.is_file():
            raise ValueError(f'Unsupported or unavailable CCD component {code}; use an installed, validated chemical definition.')
        molecule = pickle.loads(molecule_path.read_bytes())
        expected = {atom.GetProp('name'): atom.GetSymbol().upper()
                    for atom in molecule.GetAtoms() if atom.GetAtomicNum() != 1
                    and not (atom.HasProp('leaving_atom') and atom.GetBoolProp('leaving_atom'))}
        locator = row.get('residue_key')
        if not locator:
            raise ValueError(f'Unsupported CCD {code}: source residue/atom mapping is missing.')
        residues = [residue for chain in source for residue in chain
                    if f'{chain.name}:{residue.seqid}' == locator]
        observed = {atom.name: atom.element.name.upper() for residue in residues for atom in residue if not atom.is_hydrogen()}
        if observed != expected:
            raise ValueError(f'Unsupported CCD {code} source atom graph: missing {sorted(expected.keys() - observed.keys())}, extra {sorted(observed.keys() - expected.keys())}.')
        row['heavy_atoms'] = expected


def restore_prediction_ccd_names(path, task):
    """Restore Boltz's generic LIG labels using verified CCD atom identities.

    Native outputs remain untouched. This restores serialization metadata; it
    neither changes coordinates nor claims validation of inferred bond geometry.
    """
    import gemmi
    path = Path(path)
    contract = task.get('scientific_contract') or {}
    components = [row for row in contract.get('components', []) if row.get('ccd') and row.get('heavy_atoms')]
    if not components:
        return path
    structure = gemmi.read_structure(str(path))
    restored = {}
    for row in components:
        chain_id = contract['chain_map'][row['chain']]
        residues = [res for chain in structure[0] if chain.name == chain_id for res in chain]
        atoms = {atom.name: atom.element.name.upper() for res in residues for atom in res if not atom.is_hydrogen()}
        if len(residues) != 1 or atoms != row['heavy_atoms']:
            raise ValueError(f'Prediction changed CCD ligand atom identities on {chain_id}.')
        residue = residues[0]
        if residue.name == row['ccd']:
            continue
        if residue.name not in {'LIG', 'UNL'}:
            raise ValueError(f'Unexpected predicted ligand label {residue.name} on {chain_id}.')
        restored[chain_id] = row['ccd']
        row['native_output_residue_name'] = residue.name
        residue.name = row['ccd']
    if not restored:
        return path
    destination = path.with_name(path.stem + '_identity' + path.suffix)
    if path.suffix.lower() == '.pdb':
        lines = []
        for line in path.read_text().splitlines(keepends=True):
            if line.startswith(('HETATM', 'ATOM  ')) and line[21] in restored:
                line = line[:17] + f'{restored[line[21]]:>3}' + line[20:]
            lines.append(line)
        destination.write_text(''.join(lines))
    else:
        structure.make_mmcif_document().write_file(str(destination))
    return destination
