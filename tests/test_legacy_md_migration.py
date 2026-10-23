from mn_molecular_design_studio.legacy_md_migration import build_md_migration_plan


def _row(id, group, metadata=None, inputs=None, result=None):
    return {'run_id': id, 'group': group, 'path': '/source/' + group + '/' + id,
            'metadata.json': {'status': 'completed', **(metadata or {})},
            'input.json': inputs or {}, 'result.json': result or {}}


def _inventory():
    rows = [
        _row('target', 'target-trimming', {'pdb_id': '4LNW'}),
        _row('set', 'compound-selection', {'dataset_name': 'Medchem THR agonist list'}),
        _row('predict', 'refolding', {'prepared_target_run_id': 'target'},
             {'compound_sets': [{'run_id': 'set'}]}),
        _row('pose', 'selected-complexes', {'compound_id': 'HY-1', 'ligand_smiles': 'CCO',
            'prepared_target_run_id': 'target', 'parent_run_id': 'predict',
            'source_pose_id': 'HY-1__replicate_002__pose_008', 'selection_status': 'User changed automatic selection'}),
        _row('prep', 'md-system-prep', {'target_run_id': 'pose'}),
        _row('one', 'bound-ligand-md', {'target_run_id': 'pose', 'md_engine': 'openmm'},
             {'source_md_system_prep_run_id': 'prep'}, {'md_result': {'equilibration_stats': {'production': {'duration_ns': 1}}}}),
        _row('ninety-nine', 'bound-ligand-md', {'target_run_id': 'pose', 'md_engine': 'openmm',
             'continuation_of_run_id': 'one', 'target_duration_ns': 100, 'previous_duration_ns': 1},
             {'source_md_system_prep_run_id': 'prep'}, {'md_result': {'equilibration_stats': {'production': {'duration_ns': 99, 'trajectory_appended': False}}}}),
        _row('endpoint', 'md-mmgbsa', inputs={'source_production_run_id': 'ninety-nine'}),
        _row('failed', 'bound-ligand-md', {'status': 'failed', 'target_run_id': 'pose', 'md_engine': 'openmm'}),
    ]
    return {'source_root': '/source', 'runs': rows,
            'md_files': [{'run_id': r['run_id'], 'bytes': 10} for r in rows if r['group'].startswith('md-') or r['group'] == 'bound-ligand-md']}


def test_preserves_nonbest_pose_and_segment_boundary_without_inventing_best_rank():
    p = build_md_migration_plan(_inventory())
    job = next(j for j in p['jobs'] if j['stage'] == 'md')
    study = job['studies'][0]
    assert study['source_pose_rank'] == 8
    assert study['is_best_pose'] is None
    assert study['selection_status'] == 'User changed automatic selection'
    assert {r['segment_duration_ns'] for r in study['production_runs']} == {1, 99}
    extension = next(r for r in study['production_runs'] if r['run_id'] == 'ninety-nine')
    assert extension['continuation_of_run_id'] == 'one'
    assert extension['target_duration_ns'] == 100
    assert 'failed' not in job['source_run_ids']
    assert {'prep', 'pose', 'one', 'ninety-nine'} == set(job['source_run_ids'])
    assert study['dataset_name'] == 'Medchem THR agonist list'


def test_failed_dependency_is_preserved_if_a_completed_extension_requires_it():
    v = _inventory()
    next(r for r in v['runs'] if r['run_id'] == 'one')['metadata.json']['status'] = 'failed'
    p = build_md_migration_plan(v)
    md = next(j for j in p['jobs'] if j['stage'] == 'md')
    assert next(s for s in md['studies'][0]['production_runs'] if s['run_id'] == 'one')['status'] == 'failed'


def test_unidentified_dataset_is_withheld_not_guessed_from_compound_name():
    v = _inventory()
    next(r for r in v['runs'] if r['run_id'] == 'set')['metadata.json']['dataset_name'] = 'unknown'
    assert not build_md_migration_plan(v)['jobs']


def test_hundred_ns_filter_keeps_initial_segment_but_excludes_unextended_short_branch():
    v = _inventory()
    v['runs'] += [
        _row('short-only', 'bound-ligand-md', {'target_run_id': 'pose', 'md_engine': 'openmm'},
             result={'success': True, 'md_result': {'equilibration_stats': {'production': {'duration_ns': 1}}}}),
        _row('short-endpoint', 'md-mmgbsa', inputs={'source_production_run_id': 'short-only'}),
        _row('mixed-analysis', 'md-analysis', inputs={'production_run_ids': ['short-only', 'ninety-nine']}),
    ]
    p = build_md_migration_plan(v, completed_duration_ns=100)
    md = next(j for j in p['jobs'] if j['stage'] == 'md')
    assert set(md['source_run_ids']) == {'one', 'ninety-nine', 'pose', 'prep'}
    assert md['studies'][0]['qualifying_production_run_ids'] == ['ninety-nine']
    assert {r['run_id'] for r in p['excluded']} == {'short-only', 'short-endpoint', 'mixed-analysis'}
    assert not p['withheld']


def test_gromacs_step_duration_and_previous_duration_establish_hundred_ns_completion():
    v = _inventory()
    for row in v['runs']:
        if row['run_id'] in {'one', 'ninety-nine'}:
            row['metadata.json']['md_engine'] = 'gromacs'
            row['input.json'].update({'production_steps': 12500000, 'production_timestep_fs': 4})
            row['result.json'] = {'success': True, 'md_result': {'success': True}}
    next(r for r in v['runs'] if r['run_id'] == 'ninety-nine')['metadata.json']['previous_duration_ns'] = 50
    p = build_md_migration_plan(v, completed_duration_ns=100)
    md = next(j for j in p['jobs'] if j['stage'] == 'md')
    assert md['tool'] == 'gromacs_md'
    assert md['studies'][0]['qualifying_production_run_ids'] == ['ninety-nine']


def test_requested_total_alone_does_not_establish_completed_hundred_ns():
    v = _inventory()
    last = next(r for r in v['runs'] if r['run_id'] == 'ninety-nine')
    last['result.json']['md_result']['equilibration_stats']['production']['duration_ns'] = 10
    p = build_md_migration_plan(v, completed_duration_ns=100)
    assert p['jobs'] == []
