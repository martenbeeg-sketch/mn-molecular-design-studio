"""The saved-candidate → acceptance → analysis → Results vertical slice."""
import copy
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pytest
from streamlit.testing.v1 import AppTest

from mn_molecular_design_studio import candidate_sets as sets, pipeline_graph as graph
from mn_molecular_design_studio.complex_analysis_handoff import build_collection, inherit_review
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json, write_json
from mn_molecular_design_studio.engine.workflow_engine.workflows import graph_pipeline as runner


def node(kind):
    spec = graph.NODE_CATALOG[kind]
    return {'id': kind, 'type': 'studio', 'position': {'x': 0, 'y': 0},
            'data': {'kind': kind, 'label': spec['label'],
                     'params': {f['name']: copy.deepcopy(f['default']) for f in spec['fields']}}}


def fixture_source(tmp_path, design_source_factory):
    source, base = design_source_factory(tmp_path / 'runs/design/original')
    rows = []
    for identifier, score in [('accepted', .85), ('rejected', .3), ('unknown', None)]:
        row = copy.deepcopy(base)
        row.update(candidate_id=identifier, complex_pdb=str(source / 'designed.cif'),
                   parents=['native-parent'], metrics={'native_pass_filters': True})
        if score is not None:
            row['metrics']['iptm'] = score
        rows.append(row)
    path = source / 'artifacts/normalized_candidates/candidates.jsonl'
    path.write_text(''.join(json.dumps(row) + '\n' for row in rows))
    return source, rows


def merged_source(tmp_path, source, rows):
    native = runner._local_source(tmp_path, 'input', rows,
        {'target_pdb': str(source / 'designed.cif'), 'target_chains': ['A']})
    native.update(source_run_id=source.name, source_run_dir=str(source))
    return sets.merge(tmp_path, 'merge', [('input', native)], {'deduplication': 'keep_all'})


def filtered(tmp_path, source):
    return sets.shared_filter(tmp_path, 'acceptance', source,
        {'binder_format': 'miniprotein', 'target_type': 'structured', 'policy_profile': 'custom',
         'rules': [{'field': 'iptm', 'operator': '>=', 'value': .7}], 'exclude_native_rejected': True})


def test_accepted_route_keeps_all_results_and_exact_provenance(tmp_path, design_source_factory):
    source, rows = fixture_source(tmp_path, design_source_factory)
    before = hashlib.sha256((source / 'designed.cif').read_bytes()).hexdigest()
    merged = merged_source(tmp_path, source, rows)
    accepted = filtered(tmp_path, merged)
    output = build_collection(tmp_path, 'analysis', accepted, {})
    active = list(sets.iter_candidates(output))
    review = list(sets.iter_candidates(output, review=True))
    assert len(active) == 1 and len(review) == 3
    assert {row['raw_metadata']['validation_acceptance']['status'] for row in review} == {'accepted', 'rejected', 'unresolved'}
    entry = active[0]['raw_metadata']['complex_analysis_entry']
    assert entry['native_candidate_id'] == 'accepted'
    assert entry['parents'] == ['native-parent']
    assert entry['source_run_dir'] == str(source)
    assert entry['structure_sha256'] == before
    assert entry['source_structure'] == str(source / 'designed.cif')
    assert entry['target']['target_chains'] == ['A']
    assert Path(active[0]['complex_pdb']).is_relative_to(tmp_path / 'artifacts')
    assert not active[0]['raw_metadata'].get('selection', {}).get('selected')
    # Stable view identities and byte hashes on reopen, with no Database API.
    again = build_collection(tmp_path, 'analysis', accepted, {})
    assert list(sets.iter_candidates(again)) == active
    assert hashlib.sha256((source / 'designed.cif').read_bytes()).hexdigest() == before
    assert not list(tmp_path.rglob('*.sqlite*'))


def test_no_policy_or_stale_accepted_flag_never_passes(tmp_path, design_source_factory):
    source, rows = fixture_source(tmp_path, design_source_factory)
    merged = merged_source(tmp_path, source, rows)
    raw_rows = list(sets.iter_candidates(merged))
    for row in raw_rows:
        row['raw_metadata']['validation_acceptance'] = {'status': 'accepted'}
    fake = sets.write_collection(tmp_path, 'stale', raw_rows, merged['target'])
    output = build_collection(tmp_path, 'analysis', fake, {})
    assert output['candidate_count'] == 0
    assert len(list(sets.iter_candidates(output, review=True))) == 3


def test_distinct_predictions_keep_method_model_and_candidate_ids(tmp_path, design_source_factory):
    source, rows = fixture_source(tmp_path, design_source_factory)
    merged = merged_source(tmp_path, source, rows[:1])
    accepted = filtered(tmp_path, merged)
    row = next(sets.iter_candidates(accepted))
    info = sets.metadata(row)
    info['predictions'] = []
    for model_id in (0, 1):
        prediction = copy.deepcopy(rows[0])
        prediction['candidate_id'] = f'native-model-{model_id}'
        prediction['source_tool'] = 'boltz2'
        prediction['raw_metadata']['model_index'] = model_id
        info['predictions'].append({'prediction_id': f'prediction-{model_id}', 'status': 'completed',
                                    'branch': 'boltz', 'candidate': prediction})
    enriched = sets.write_collection(tmp_path, 'enriched', [row], merged['target'])
    output = build_collection(tmp_path, 'analysis', enriched, {})
    entries = list(sets.iter_candidates(output))
    assert len(entries) == 3
    receipts = [row['raw_metadata']['complex_analysis_entry'] for row in entries]
    assert {r['prediction_id'] for r in receipts} == {'', 'prediction-0', 'prediction-1'}
    assert {r['model_id'] for r in receipts} == {'', '0', '1'}
    assert len({row['candidate_id'] for row in entries}) == 3
    assert len({row['complex_pdb'] for row in entries}) == 3
    assert len({r['variant_id'] for r in receipts}) == 1


def test_manual_complex_entry_does_not_change_acceptance_or_selection(tmp_path):
    import gemmi
    from mn_molecular_design_studio.complex_analysis_ui import component_inventory
    from mn_molecular_design_studio.engine.workflow_engine.workflows.quality_ligands import prepare
    assets = Path(__file__).resolve().parents[1] / 'mn_molecular_design_studio/assets/vina'
    structure = gemmi.read_structure(str(assets / '1iep.pdb'))
    for model in structure:
        for chain in list(model):
            if chain.name != 'A':
                model.remove_chain(chain.name)
    path = tmp_path / 'protein-imatinib.pdb'
    structure.write_pdb(str(path))
    components = component_inventory(path.read_bytes(), path.name)
    smiles = json.loads((assets / 'imatinib.json').read_text())[0]['smiles']
    row = {'candidate_id': 'original-1iep-A', 'source_tool': 'historical_structure', 'stage': 'analysis',
           'complex_pdb': str(path), 'target_chains': ['A'], 'binder_chains': [],
           'raw_metadata': {'components': components, 'ligand_smiles': smiles}}
    original = tmp_path / 'prediction/historical-ligand'
    normalized = original / 'artifacts/normalized_candidates/candidates.jsonl'
    normalized.parent.mkdir(parents=True)
    normalized.write_text(json.dumps(row) + '\n')
    write_json(original / 'metadata.json', {'task_group': 'prediction', 'status': 'completed'})
    write_json(original / 'result.json', {'success': True})
    write_json(original / 'molecular_design_studio.json', {'project_id': 'project'})
    loaded_source, loaded = runner.source_project_designs(str(original), 'project', ['original-1iep-A'], include_complex_results=True)
    assert loaded_source == original and len(loaded) == 1
    with pytest.raises(ValueError, match='completed design job'):
        runner.source_project_designs(str(original), 'project')
    native = runner._local_source(tmp_path, 'native', loaded, {'target_pdb': str(path), 'target_chains': ['A']})
    merged = sets.merge(tmp_path, 'merge', [('native', native)], {'deduplication': 'keep_all'})
    identifier = next(sets.iter_candidates(merged))['candidate_id']
    output = build_collection(tmp_path, 'analysis', merged,
        {'entry_mode': 'manual_selection', 'selected_candidate_ids': [identifier]})
    item = next(sets.iter_candidates(output))
    assert item['raw_metadata']['ligand_smiles'] == smiles
    assert item['raw_metadata']['complex_analysis_entry']['acceptance'] == {}
    assert not item['raw_metadata'].get('selection', {}).get('selected')
    assert len([c for c in components if c['type'] == 'protein']) == 1
    assert any(c['type'] == 'small_molecule' for c in components)
    prepared = {'id': 'ligand', 'candidate_id': identifier}
    prepare(prepared, item, tmp_path / 'plip', ligand_selection={
        'chain': 'A', 'residue_key': 'A:201', 'ccd': 'STI', 'smiles': smiles})
    assert prepared['ligand_identity']['residue_key'] == 'A:201'
    with pytest.raises(ValueError, match='explicit candidate IDs'):
        build_collection(tmp_path, 'invalid', merged, {'entry_mode': 'manual_selection'})


def test_existing_analysis_tool_preserves_rejected_and_unknown_review(tmp_path, design_source_factory):
    source, rows = fixture_source(tmp_path, design_source_factory)
    output = build_collection(tmp_path, 'analysis', filtered(tmp_path, merged_source(tmp_path, source, rows)), {})
    measured = sets.write_collection(tmp_path, 'biotite', list(sets.iter_candidates(output)), output['target'])
    result = inherit_review(tmp_path, 'biotite', measured, output)
    assert len(list(sets.iter_candidates(result))) == 1
    assert len(list(sets.iter_candidates(result, review=True))) == 3


@pytest.mark.parametrize('include_tools', [False, True])
def test_full_graph_acceptance_analysis_results_and_reopen(tmp_path, monkeypatch, design_source_factory, include_tools):
    from mn_molecular_design_studio import workspace
    source, _ = fixture_source(tmp_path, design_source_factory)
    parent = tmp_path / 'runs/campaign/vertical-slice'
    parent.mkdir(parents=True)
    write_json(parent / 'metadata.json', {'status': 'queued', 'job_type': 'visual_pipeline'})
    write_json(parent / 'input.json', {'params': {}})
    monkeypatch.setattr(graph, 'app_home', lambda: tmp_path / 'home')
    monkeypatch.setattr(workspace, 'mark_studio_job', lambda *a, **k: None)
    monkeypatch.setattr(runner, 'create_job', lambda *a, **k: SimpleNamespace(run_dir=parent))
    g = graph.new_graph('accepted_complex_analysis')
    if not include_tools:
        g['nodes'] = [n for n in g['nodes'] if n['id'] != 'biotite_interface']
        g['edges'] = [e for e in g['edges'] if 'biotite_interface' not in (e['source'], e['target'])]
        g['edges'].append({'id': 'analysis-results', 'source': 'complex_analysis', 'target': 'results'})
    def submit(request, parent, node, **kwargs):
        from mn_molecular_design_studio.engine.workflow_engine.core.candidates import write_candidates
        child = tmp_path / 'runs/analysis/contacts'
        child.mkdir(parents=True, exist_ok=True)
        input_source = request['parameters']
        measured = []
        for original in sets.iter_candidates({'run_dir': input_source['source_run_dir'], 'candidates_jsonl': input_source['candidates_jsonl']}):
            row = copy.deepcopy(original)
            row.update(candidate_id=original['candidate_id'] + '_biotite', parents=[original['candidate_id']],
                       source_tool='biotite_interface', stage='analysis', metrics={'qc_interface_contact_pairs': 12})
            measured.append(row)
        write_candidates(child, 'biotite_interface', measured)
        write_json(child / 'metadata.json', {'status': 'completed', 'tool': 'biotite_interface'})
        write_json(child / 'result.json', {'success': True, 'outputs': {'candidates': measured}})
        return child
    monkeypatch.setattr(runner, '_submit', submit)
    next(n for n in g['nodes'] if n['id'] == 'project_design_source')['data']['params'].update(source_run_dir=str(source), selected_candidate_ids=['accepted', 'rejected', 'unknown'])
    next(n for n in g['nodes'] if n['id'] == 'shared_filter')['data']['params'].update(
        target_type='structured', policy_profile='custom', rules=[{'field': 'iptm', 'operator': '>=', 'value': .7}])
    runner.run_graph_pipeline(g, 'Vertical slice', 'test-user', 'test-project')
    result = read_json(parent / 'result.json')
    assert result['success'] is True
    assert len(result['outputs']['candidates']) == 3
    state = read_json(parent / 'graph_execution.json')
    assert state['nodes']['complex_analysis']['candidate_count'] == 1
    candidates_before = result['outputs']['candidates']
    runner.run_graph_pipeline(g, 'Vertical slice', 'test-user', 'test-project')
    assert read_json(parent / 'result.json')['outputs']['candidates'] == candidates_before
    assert len(list((parent / 'artifacts/pipeline_nodes/complex_analysis/artifacts/structures').iterdir())) == 1


def test_analysis_node_renders_in_canvas():
    app = AppTest.from_string('''
import streamlit as st
from mn_molecular_design_studio.pipeline_graph import new_graph
from mn_molecular_design_studio.pipeline_node_ui import render_node_form
g = new_graph('accepted_complex_analysis')
n = next(n for n in g['nodes'] if n['id'] == 'complex_analysis')
render_node_form(n, g, key='analysis-node', render_input=lambda *a, **k: None, owner_id='test-user')
''').run()
    assert not app.exception
    assert any(item.label == 'Analysis entry' for item in app.selectbox)
    app.selectbox(key='analysis-node-entry_mode').select('manual_selection').run()
    assert not app.exception


def test_target_superposition_reuses_fitter_and_binds_exact_bytes(tmp_path, design_source_factory):
    import gemmi
    from mn_molecular_design_studio.complex_analysis_comparison import fit_matching_targets, protein_chain_numbering
    source, _ = design_source_factory(missing_occupancy=False)
    reference_path = source / 'designed.cif'
    structure = gemmi.read_structure(str(reference_path))
    for atom in structure[0].all():
        atom.atom.pos += gemmi.Position(2, 3, -1)
    moving_path = tmp_path / 'translated.cif'
    structure.make_mmcif_document().write_file(str(moving_path))
    def sample(path, identifier):
        return {'identity': identifier, 'model_id': '1', 'path': path,
                'structure_digest': hashlib.sha256(path.read_bytes()).hexdigest(),
                'complex_type': 'protein-protein', 'target_chains': ['A'], 'binder_chains': ['B'],
                'target_chains_sequence': ['ACDEFG'], 'binder_chains_sequence': ['KLMNP'],
                'target_numbering': protein_chain_numbering(path, ['A'])}
    reference, moving = sample(reference_path, 'reference'), sample(moving_path, 'moving')
    fit = fit_matching_targets(reference, moving)
    assert fit['mapped_residues'] == 6
    assert fit['target_ca_rmsd_angstrom'] < 1e-5
    assert [fit['matrix_4x4'][i][3] for i in range(3)] == pytest.approx([-2, -3, 1])
    moving_path.write_text(moving_path.read_text() + '\n')
    with pytest.raises(ValueError, match='bytes changed'):
        fit_matching_targets(reference, moving)
