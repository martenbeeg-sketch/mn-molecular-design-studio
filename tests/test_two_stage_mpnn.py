"""Two-stage routing and preservation regressions; inference is captured by pytest."""
import json

import pytest

from test_mpnn import capture_runner
from mn_molecular_design_studio import pipeline_graph as graph
from mn_molecular_design_studio.engine.workflow_engine.workflows import mpnn, sequence_design as native


def test_comparison_has_initial_and_optional_soluble_stages():
    g = graph.new_graph('all_generators_boltz2')
    stages = [n for n in g['nodes'] if n['data']['kind'] == 'ligandmpnn']
    assert len(stages) == 2
    assert stages[0]['data']['params']['sequence_stage'] == 'initial'
    assert stages[1]['data']['params']['sequence_stage'] == 'soluble_optimization'
    assert stages[1]['data']['params']['optimization_enabled'] is False


def test_initial_codesign_preserves_original_without_running_mpnn(capture_runner):
    source, candidate, collection, run, commands, _ = capture_runner
    candidate['stage'] = 'generation.backbone_sequence'
    collection.write_text(json.dumps(candidate) + '\n')
    mpnn.run_mpnn_sequence_design(source, collection, sequence_stage='initial')
    assert not commands
    rows = json.loads((run / 'artifacts/normalized_candidates/candidates.jsonl').read_text())
    assert rows['binder_sequence'] == candidate['binder_sequence']
    assert rows['candidate_id'] == candidate['candidate_id']


def test_initial_backbone_runs_proteinmpnn(capture_runner):
    source, _, collection, _, commands, _ = capture_runner
    mpnn.run_mpnn_sequence_design(source, collection, sequence_stage='initial')
    assert commands[0][commands[0].index('--model_type') + 1] == 'protein_mpnn'


def test_interface_fixing_rejects_mismatched_sequence(capture_runner):
    source, candidate, collection, _, commands, _ = capture_runner
    candidate['binder_sequence'] = 'WRONG'
    collection.write_text(json.dumps(candidate) + '\n')
    with pytest.raises(ValueError, match='sequence.*structure'):
        mpnn.run_mpnn_sequence_design(source, collection, model_type='soluble_mpnn', fix_interface_residues=True)
    assert not commands


def test_configurable_contact_cutoff_reaches_existing_implementation(capture_runner, monkeypatch):
    source, _, collection, _, commands, _ = capture_runner
    observed = []
    monkeypatch.setattr(native, '_interface_contacts_for_structure', lambda *a, **k: observed.append(k['cutoff']) or ['Z1'])
    mpnn.run_mpnn_sequence_design(source, collection, model_type='soluble_mpnn', fix_interface_residues=True, interface_cutoff=6.5)
    assert observed == [6.5]
    assert '--fixed_residues' in commands[0]




def _fake_valid_design(capture_runner, monkeypatch, *, damage=None):
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import finish_job
    source, candidate, collection, run, commands, _ = capture_runner
    monkeypatch.setattr(native, 'finish_job', finish_job)
    def execute(directory, steps):
        commands.extend(s['command'] for s in steps)
        staged = run / 'artifacts/raw/ligandmpnn/inputs' / (candidate['candidate_id'] + '.pdb')
        output = run / 'artifacts/raw/ligandmpnn/output' / candidate['candidate_id']
        (output / 'seqs').mkdir(parents=True)
        (output / 'backbones').mkdir()
        lines = []
        for line in staged.read_text().splitlines():
            if line.startswith('ATOM'):
                chain, number = line[21], int(line[22:26])
                if (chain == 'Z' and number == (1 if damage == 'interface' else 2 if damage == 'motif' else 3)) or (damage == 'target' and chain == 'A' and number == 1):
                    line = line[:17] + ('VAL' if damage == 'target' and chain == 'A' else 'ALA') + line[20:]
            lines.append(line)
        backbone = output / 'backbones' / (candidate['candidate_id'] + '_1.pdb')
        backbone.write_text('\n'.join(lines) + '\n')
        seqs = native.refolding_workflow._sequences_by_chain(backbone)
        (output / 'seqs/design.fa').write_text('>design, id=1\n' + ':'.join(seqs[c] for c in sorted(seqs)) + '\n')
        return 0
    monkeypatch.setattr(native, '_run_shell_steps', execute)
    monkeypatch.setattr(native, '_interface_contacts_for_structure', lambda *a, **k: ['B1'])
    candidate['stage'] = 'generation.backbone_sequence'
    candidate['raw_metadata'] = {'fixed_residues': ['B2']}
    collection.write_text(json.dumps(candidate) + '\n')
    return source, collection, run


def test_soluble_variants_coexist_with_parent_in_results(capture_runner, monkeypatch, tmp_path):
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
    from mn_molecular_design_studio.engine.workflow_engine.workflows.graph_pipeline import _local_source, _absolute_candidates
    source, collection, run = _fake_valid_design(capture_runner, monkeypatch)
    mpnn.run_mpnn_sequence_design(source, collection, sequence_stage='soluble_optimization', optimization_enabled=True, seed=1)
    rows = read_candidates(run)
    assert [r['binder_sequence'] for r in rows] == ['KLMNP', 'KLANP']
    assert rows[1]['parents'] == [rows[0]['candidate_id']]
    assert rows[1]['raw_metadata']['fixed_sequence_validation']['status'] == 'verified'
    results = _local_source(tmp_path, 'results', rows, final_evidence=True)
    assert [r['binder_sequence'] for r in _absolute_candidates(results)] == ['KLMNP', 'KLANP']
    from mn_molecular_design_studio.output_vocabulary import candidate_evidence
    evidence = candidate_evidence(rows[1])
    assert evidence[1]['role'] == 'comparator'
    assert evidence[1]['candidate']['binder_sequence'] == 'KLMNP'


@pytest.mark.parametrize('damage', ['interface', 'motif', 'target'])
def test_invalid_redesign_is_rejected_after_inference(capture_runner, monkeypatch, damage):
    source, collection, run = _fake_valid_design(capture_runner, monkeypatch, damage=damage)
    with pytest.raises(ValueError, match='fixed residues|target-chain'):
        mpnn.run_mpnn_sequence_design(source, collection, sequence_stage='soluble_optimization', optimization_enabled=True)


def test_disabled_soluble_stage_does_not_run_inference(capture_runner):
    source, candidate, collection, run, commands, _ = capture_runner
    candidate['stage'] = 'sequence_design'
    collection.write_text(json.dumps(candidate) + '\n')
    mpnn.run_mpnn_sequence_design(source, collection, sequence_stage='soluble_optimization', optimization_enabled=False)
    assert not commands
    assert json.loads((run / 'artifacts/normalized_candidates/candidates.jsonl').read_text())['binder_sequence'] == 'KLMNP'




def test_fasta_chain_order_matches_installed_ligandmpnn(capture_runner):
    source, candidate, collection, run, _, _ = capture_runner
    # Structure order B,A must not be confused with native FASTA order A,B.
    path = source / 'reordered.pdb'
    lines = (source / 'designed.pdb').read_text().splitlines()
    path.write_text('\n'.join([*[line for line in lines if line.startswith('ATOM') and line[21] == 'B'],
        *[line for line in lines if line.startswith('ATOM') and line[21] == 'A'], 'END', '']))
    assert native._designed_sequence_from_multichain_fasta('ACDEFG:KLMNP', path, ['B']) == 'KLMNP'
    with pytest.raises(ValueError, match='chain count'):
        native._designed_sequence_from_multichain_fasta('ACDEFG:KLMNP:WRONG', path, ['B'])


def test_legacy_default_comparison_migrates_without_changing_explicit_redesign():
    import copy
    g = graph.new_graph('all_generators_boltz2')
    optional = next(n for n in g['nodes'] if n['id'] == 'optional-soluble-mpnn')
    initial = next(n for n in g['nodes'] if n['data']['kind'] == 'ligandmpnn')
    g['nodes'].remove(optional)
    g['edges'] = [e for e in g['edges'] if e['target'] != optional['id']]
    for e in g['edges']:
        if e['source'] == optional['id']:
            e['source'] = initial['id']
    for key in ('sequence_stage', 'optimization_enabled', 'interface_cutoff'):
        initial['data']['params'].pop(key)
    g['template_version'] = 3
    explicit = copy.deepcopy(g)
    next(n for n in explicit['nodes'] if n['id'] == initial['id'])['data']['params']['sampling_temp'] = .3
    migrated = graph.migrate_template_graph(g)
    assert migrated['template_version'] == 4
    assert len([n for n in migrated['nodes'] if n['data']['kind'] == 'ligandmpnn']) == 2
    assert graph.migrate_template_graph(explicit) == explicit


def test_inherited_motif_and_framework_survive_chain_renaming(capture_runner):
    import copy
    source, original, _, _, _, _ = capture_runner
    original['raw_metadata'] = {'fixed_residues': ['B2'], 'sampled_binder_prompt': 'K####'}
    child = copy.deepcopy(original)
    child['raw_metadata'] = {'source_candidate': original}
    child['binder_chains'] = ['Z']
    path = source / 'renamed.pdb'
    path.write_text((source / 'designed.pdb').read_text().replace(' B ', ' Z '))
    child['complex_pdb'] = str(path)
    assert mpnn.inherited_fixed_residues(child, path, ['Z']) == ['Z1', 'Z2']


def test_native_sequence_conditions_cannot_be_silently_dropped(capture_runner):
    source, candidate, _, _, _, _ = capture_runner
    candidate['raw_metadata'] = {'reference_constraints': {'condition_spec': "binder:\n  conditions:\n    C_SEQ:\n      true:\n        - select: \"sel('B/2')\"\n"}}
    with pytest.raises(ValueError, match='verified fixed-residue mapping'):
        mpnn.inherited_fixed_residues(candidate, source / 'designed.pdb', ['B'])


def test_downstream_filters_use_variant_evidence_without_borrowing_parent_scores(capture_runner, monkeypatch):
    import copy
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
    from mn_molecular_design_studio.output_vocabulary import metric_value
    from mn_molecular_design_studio.candidate_selection import evaluate_selection, selected_candidates
    source, collection, run = _fake_valid_design(capture_runner, monkeypatch)
    mpnn.run_mpnn_sequence_design(source, collection, sequence_stage='soluble_optimization', optimization_enabled=True)
    parent, variant = read_candidates(run)
    variant['metrics']['overall_confidence'] = .8
    predictions = []
    for row in (parent, variant):
        predictions.append({'candidate_id': 'predicted-' + row['candidate_id'], 'stage': 'complex_refolding',
            'source_tool': 'boltz2', 'binder_sequence': row['binder_sequence'], 'metrics': {'iptm': .7},
            'parents': [row['candidate_id']], 'raw_metadata': {'source_candidate': copy.deepcopy(row)}})
    field = 'evidence/sequence_design/soluble_mpnn/overall_confidence'
    assert metric_value(predictions[0], field) is None
    assert metric_value(predictions[1], field) == .8
    selection = evaluate_selection(predictions, {'rules': [{'field': field, 'operator': '>=', 'value': .5}]})
    assert [c['candidate_id'] for c in selected_candidates(selection)] == [predictions[1]['candidate_id']]
