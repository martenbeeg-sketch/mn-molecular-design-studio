from mn_molecular_design_studio.pipeline_graph import new_graph, validate_graph

def test_vhh_scaffold_graph_preserves_generator_sequences():
    plan = validate_graph(new_graph("vhh_design"))
    kinds = [n["data"]["kind"] for n in plan["graph"]["nodes"]]
    assert "ligandmpnn" not in kinds
    assert {"bindcraft2_generator", "boltzgen_generator", "esm_binder", "boltz2"} <= set(kinds)


import pytest
from mn_molecular_design_studio import mpnn_input, sabr_runtime
from mn_molecular_design_studio.antibody_constraints import verify_masked_sequence, verify_fixed_sequence


@pytest.mark.parametrize('modality', ['VHH', 'scFv', 'Fab', 'ARP'])
def test_only_scaffold_mode_allows_empty_binder_length(modality):
    g = new_graph('vhh_design')
    node = next(n for n in g['nodes'] if n['data']['kind'] == 'bindcraft2_generator')
    node['data']['params'].update(modality=modality, binder_length='')
    validate_graph(g)
    node['data']['params']['modality'] = 'binder'
    with pytest.raises(ValueError, match='Binder length'):
        validate_graph(g)


def test_selected_cdr_cysteine_is_fixed():
    from mn_molecular_design_studio import design_validation as dv
    path = mpnn_input.ASSETS / '9F6O.pdb'
    residues = dv.protein_chains(dv.read_structure(path))['B']
    cysteines = [f'B{r.id[1]}{r.id[2].strip()}' for r in residues if r.resname == 'CYS']
    selected = [f'B{r.id[1]}{r.id[2].strip()}' for r in residues[-20:]] + cysteines
    fixed, cdrs = sabr_runtime.fixed_residues_from_annotation(path, ['B'], ['CDR3'], {'CDR3': selected})
    assert set(cysteines) <= set(fixed)
    assert not set(cysteines) & set(cdrs['CDR3'])
    assert not set(cdrs['CDR3']) & set(fixed)


def test_framework_and_unselected_positions_cannot_change():
    assert verify_masked_sequence('ACWGCT', 'AC##CT')['fixed_position_count'] == 4
    with pytest.raises(ValueError, match='framework'):
        verify_masked_sequence('AAWGCT', 'AC##CT')
    with pytest.raises(ValueError, match='length'):
        verify_masked_sequence('ACWGCTT', 'AC##CT')
    assert verify_fixed_sequence('ACWGCT', ['V1', 'V2', 'V3', 'V4', 'V5', 'V6'], 'ACAACT', ['V1', 'V2', 'V5', 'V6'])['status'] == 'verified'
    with pytest.raises(ValueError, match='fixed residues'):
        verify_fixed_sequence('AAWGCT', ['V1', 'V2', 'V3', 'V4', 'V5', 'V6'], 'ACAACT', ['V1', 'V2', 'V5', 'V6'])


def test_saved_default_vhh_drafts_migrate_without_generic_redesign():
    from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG
    g = new_graph('vhh_design'); g.pop('template_version', None)
    node = next(node for node in g['nodes'] if node['data']['kind'] == 'sabr')
    node['data'] = {'kind': 'ligandmpnn', 'label': 'ProteinMPNN', 'params': {
        item['name']: item['default'] for item in NODE_CATALOG['ligandmpnn']['fields']}}
    node['data']['params']['sequences'] = 1
    plan = validate_graph(g)
    assert plan['graph']['template_version'] == 2
    assert next(item for item in plan['graph']['nodes'] if item['id'] == node['id'])['data']['kind'] == 'sabr'
    node['data']['params']['model_type'] = 'ab_mpnn'
    plan = validate_graph(g)
    assert next(item for item in plan['graph']['nodes'] if item['id'] == node['id'])['data']['kind'] == 'ligandmpnn'


def test_abmpnn_missing_matching_structure_is_not_accepted(design_source_factory, tmp_path):
    from mn_molecular_design_studio.engine.workflow_engine.workflows.sequence_design import _normalize_ligandmpnn_candidates
    run, source = design_source_factory(missing_occupancy=False)
    destination = tmp_path / 'redesign'
    sequences = destination / 'artifacts/raw/ligandmpnn/output/passing-design/seqs'
    sequences.mkdir(parents=True)
    (sequences / 'sample.fa').write_text('>sample id=1\nACDEFG/KLMNP\n')
    with pytest.raises(ValueError, match='matching output structure'):
        _normalize_ligandmpnn_candidates(destination, [source], {'passing-design': run / 'designed.cif'},
            {'model_type': 'ab_mpnn', 'design_chains_by_candidate': {'passing-design': ['B']}})
