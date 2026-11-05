#!/usr/bin/env python3
"""Run a bounded saved-campaign slice in a fresh, isolated Studio home.

No generation, GPU inference or source writes. CPU analysis uses the existing
pinned QC image. Dispatch is synchronous for this demonstration; the saved graph
and all native node requests are the same ones used by the canvas controller.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True, help='Fresh directory, never a production Studio home')
    parser.add_argument('--count', type=int, default=8)
    args = parser.parse_args()
    if not 2 <= args.count <= 20:
        parser.error('--count must be between 2 and 20')
    source = args.source.resolve(strict=True)
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    os.environ['MN_MOLECULAR_DESIGN_STUDIO_HOME'] = str(output / 'studio')
    os.environ['MN_MOLECULAR_DESIGN_STUDIO_APP_HOME'] = str(output / 'studio')
    os.environ['MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR'] = str(output / 'studio/workdir/runs')

    from mn_molecular_design_studio import candidate_sets as sets, pipeline_graph
    from mn_molecular_design_studio.complex_analysis_comparison import (
        fit_matching_targets, protein_chain_numbering, compatibility_report)
    from mn_molecular_design_studio.design_validation import read_structure, protein_chains, sequence
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json, write_json
    from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import WORKFLOWS, _get_function
    from mn_molecular_design_studio.engine.workflow_engine.workflows import graph_pipeline as runner
    from mn_molecular_design_studio.workspace import mark_studio_job
    from mn_molecular_design_studio.pipeline_jobs import link_pipeline_child

    verification = json.loads((source / 'source-verification.json').read_text())
    hashes = {item['path']: item['sha256'] for item in verification.get('files', [])}
    for relative, expected in hashes.items():
        assert hashlib.sha256((source / relative).read_bytes()).hexdigest() == expected, relative
    rows = read_candidates(source / 'artifacts/normalized_candidates/candidates.jsonl')[:args.count]
    selected = [row['candidate_id'] for row in rows]
    marker = read_json(source / 'molecular_design_studio.json')
    project = marker['project_id']
    owner = marker['studio_owner_id']
    studio = output / 'studio'
    (studio / 'organization').mkdir(parents=True, exist_ok=True)
    (studio / 'settings.json').write_text(json.dumps({'results_directory': str(studio / 'workdir/runs'), 'reference_directory': str(studio / 'reference_files')}))
    (studio / 'organization/organization.json').write_text(json.dumps({'id': 'isolated-demo', 'projects': [{'id': project, 'name': 'R3 isolated demonstration'}], 'members': []}))
    graph = pipeline_graph.new_graph('accepted_complex_analysis')
    for node in graph['nodes']:
        if node['id'] == 'project_design_source':
            node['data']['params'].update(source_run_dir=str(source), selected_candidate_ids=selected)
        elif node['id'] == 'shared_filter':
            # Explicit stricter exploratory IDR selection, distinct from the
            # original BindCraft2 pass decisions. Native i_pDAE is unitless.
            node['data']['params'].update(policy_profile='custom', rules=[
                {'field': 'common.confidence.validation.complex.iptm', 'operator': '>=', 'value': .75},
                {'field': 'i_pDAE', 'operator': '>=', 'value': .5}])
    write_json(output / 'graph.json', graph)
    calls = []

    def synchronous_dispatch(request, parent, node, **kwargs):
        from mn_molecular_design_studio.runtime import runs_root
        assert runs_root().is_relative_to(output), 'workspace escaped the isolated output'
        assert request['workflow'] == 'analysis.biotite_interface', request['workflow']
        child = _get_function(WORKFLOWS[request['workflow']])(**request['parameters'])
        mark_studio_job({'job': {'run_dir': str(child)}}, owner_id=owner, project_id=project)
        link_pipeline_child(child, parent_run_dir=parent, parent_node_id=node['id'],
                            pipeline_name='R3 accepted-only demonstration', step=node['data']['label'])
        calls.append(str(child))
        return child

    with patch.object(runner, '_submit', synchronous_dispatch):
        run = runner.run_graph_pipeline(graph, 'R3 accepted-only demonstration', owner, project)
    state = read_json(run / 'graph_execution.json')
    result = read_json(run / 'result.json')
    assert result['success'], {key: value.get('error') for key, value in state['nodes'].items() if value.get('error')}
    entries = list(sets.iter_candidates(state['nodes']['complex_analysis']['output']))
    assert 2 <= len(entries) < len(rows), (len(entries), len(rows))
    measurements = read_candidates(Path(calls[0]) / 'artifacts/normalized_candidates/candidates.jsonl')
    assert len(measurements) == len(entries)
    assert all(row['metrics'].get('analysis_status') == 'measured' for row in measurements), [(r['candidate_id'], r['metrics'].get('analysis_status')) for r in measurements]
    assert len(result['outputs']['candidates']) == len(rows)
    before = json.dumps(result['outputs']['candidates'], sort_keys=True)
    # Reopening/recovering the same saved execution must reuse completed nodes.
    with patch.object(runner, 'create_job', return_value=SimpleNamespace(run_dir=run)), patch.object(
            runner, '_submit', side_effect=AssertionError('reopen must not resubmit scientific analyses')):
        runner.run_graph_pipeline(graph, 'R3 accepted-only demonstration', owner, project)
    assert json.dumps(read_json(run / 'result.json')['outputs']['candidates'], sort_keys=True) == before

    samples = []
    for row in entries[:2]:
        path = Path(row['complex_pdb'])
        chains = protein_chains(read_structure(path))
        samples.append({'identity': row['candidate_id'], 'model_id': path.name, 'path': path,
            'structure_digest': hashlib.sha256(path.read_bytes()).hexdigest(), 'complex_type': 'protein-protein',
            'target_chains': row['target_chains'], 'binder_chains': row['binder_chains'],
            'target_chains_sequence': ':'.join(sequence(chains[c]) for c in row['target_chains']),
            'binder_chains_sequence': ':'.join(sequence(chains[c]) for c in row['binder_chains']),
            'target_numbering': protein_chain_numbering(path, row['target_chains'])})
    fit = fit_matching_targets(*samples)
    write_json(output / 'target-fit.json', fit)
    for relative, expected in hashes.items():
        assert hashlib.sha256((source / relative).read_bytes()).hexdigest() == expected
    # The same collection handoff accepts a manually selected, experimentally
    # sourced protein-ligand complex. Keep only the explicit A/STI assembly.
    import gemmi
    import shutil
    from mn_molecular_design_studio.complex_analysis_handoff import build_collection
    from mn_molecular_design_studio.complex_analysis_ui import component_inventory
    from mn_molecular_design_studio.engine.workflow_engine.workflows.quality import run_plip
    assets = Path(__file__).resolve().parents[1] / 'mn_molecular_design_studio/assets/vina'
    ligand_source = output / 'protein-ligand-source'
    ligand_source.mkdir()
    crystal = gemmi.read_structure(str(assets / '1iep.pdb'))
    for model in crystal:
        for chain in list(model):
            if chain.name != 'A':
                model.remove_chain(chain.name)
    ligand_path = ligand_source / '1iep-chain-A.pdb'
    crystal.write_pdb(str(ligand_path))
    components = component_inventory(ligand_path.read_bytes(), ligand_path.name)
    smiles = json.loads((assets / 'imatinib.json').read_text())[0]['smiles']
    ligand = {'candidate_id': '1iep-A-STI-201', 'source_tool': 'historical_structure', 'stage': 'analysis',
        'complex_pdb': str(ligand_path), 'target_chains': ['A'], 'binder_chains': [],
        'raw_metadata': {'components': components, 'ligand_smiles': smiles,
                         'original_structure_sha256': hashlib.sha256((assets / '1iep.pdb').read_bytes()).hexdigest()}}
    native = runner._local_source(output, 'ligand-input', [ligand],
        {'target_pdb': str(ligand_path), 'target_chains': ['A']})
    collection = sets.merge(output, 'ligand-collection', [('ligand-input', native)], {'deduplication': 'keep_all'})
    candidate_id = next(sets.iter_candidates(collection))['candidate_id']
    handoff = build_collection(output, 'ligand-analysis', collection,
        {'entry_mode': 'manual_selection', 'selected_candidate_ids': [candidate_id]})
    ligand_run = run_plip(Path(handoff['run_dir']), Path(handoff['candidates_jsonl']), cpu_cores=1,
        ligand_selection={'chain': 'A', 'residue_key': 'A:201', 'ccd': 'STI', 'smiles': smiles})
    ligand_measured = read_candidates(ligand_run / 'artifacts/normalized_candidates/candidates.jsonl')
    assert len(ligand_measured) == 1 and ligand_measured[0]['metrics']['analysis_status'] == 'measured'
    assert ligand_measured[0]['raw_metadata']['complex_analysis_entry']['acceptance'] == {}
    write_json(output / 'protein-ligand-demonstration.json', {
        'source': 'tracked 1iep.pdb, explicit chain A / STI A:201', 'handoff': handoff,
        'analysis_run': str(ligand_run), 'metrics': ligand_measured[0]['metrics'],
        'ppi_interface_analysis': 'unavailable: only one selected protein chain',
        'entry_mode': 'manual_selection'})
    from mn_molecular_design_studio.engine.workflow_engine.workflows.shared_campaign import SharedCampaignLedger
    stop_config = {'enabled': True, 'objective_type': 'accepted', 'generator_nodes': ['bindcraft', 'rfd4'],
        'target_node': 'target', 'unique_accepted': True, 'accepted_target': 2, 'max_attempts': 8,
        'max_concurrent_jobs': 2, 'max_concurrent_generation_jobs': 2, 'max_concurrent_validations': 1,
        'generator_limits': {}, 'validation_stages': []}
    ledger_path = output / 'accepted-stop-ledger.json'
    ledger = SharedCampaignLedger.create(ledger_path, stop_config)
    for index, row in enumerate(entries):
        branch = 'bindcraft' if index % 2 == 0 else 'rfd4'
        claim = ledger.claim_candidate('bounded-demo', branch, row, native_acceptance='accepted')
        if claim['role'] == 'owner':
            ledger.finish_candidate('bounded-demo', row['candidate_id'], 'accepted')
        if ledger.summary()['unique_accepted'] == 2:
            break
    recovered = SharedCampaignLedger.create(ledger_path, stop_config)
    assert recovered.summary()['unique_accepted'] == 2
    assert not recovered.reserve_attempt('bounded-demo', 'after-target', 'bindcraft', 'bindcraft2')
    import sqlite3
    database = output / 'studio/database/molecular_library.sqlite3'
    with sqlite3.connect(database.as_uri() + '?mode=ro', uri=True) as connection:
        scientific_counts = {table: connection.execute('SELECT count(*) FROM ' + table).fetchone()[0]
                             for table in ('molecules', 'molecule_observations', 'molecule_structures')}
    assert not any(scientific_counts.values()), scientific_counts
    summary = {'source_run_dir': str(source), 'original_source': verification, 'selected_candidate_ids': selected,
        'run_dir': str(run), 'accepted_count': len(entries), 'review_count': len(result['outputs']['candidates']),
        'native_analysis_runs': calls, 'analysis_metrics': [{'candidate_id': row['candidate_id'], 'metrics': row['metrics']} for row in measurements],
        'comparison': compatibility_report(*samples), 'target_fit': fit,
        'reopen_reused_scientific_records': True, 'source_checksums_unchanged': True,
        'curated_promotion_requested': False, 'isolated_database_counts': scientific_counts,
        'stop_at_two_accepted': recovered.summary(), 'protein_ligand_analysis_run': str(ligand_run)}
    write_json(output / 'demonstration.json', summary)
    print(json.dumps({'report': str(output / 'demonstration.json'), 'accepted': len(entries),
                      'retained_in_results': len(result['outputs']['candidates']),
                      'target_ca_rmsd_angstrom': fit['target_ca_rmsd_angstrom']}))


if __name__ == '__main__':
    main()
