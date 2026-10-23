"""Plan legacy MD imports without replacing the recorded pose or segment lineage.

Input is a read-only inventory; no scientific inference or trajectory stitching is
performed here. Every completed endpoint/analysis retains its exact source runs.
"""
from __future__ import annotations

import copy
import hashlib
import math
from collections import defaultdict
from pathlib import Path

from mn_molecular_design_studio.legacy_migration import DATASET_PROJECTS, THR_TARGETS

MD_GROUPS = frozenset({'bound-ligand-md', 'md-system-prep', 'md-mmgbsa', 'md-analysis'})
REFERENCE_PROJECT = 'Thyroid Hormon Recptor interaction MD references'


def _meta(row):
    return row.get('metadata.json') or {}


def _inputs(row):
    return row.get('input.json') or {}


def _digest(value):
    return hashlib.sha256(value.encode()).hexdigest()[:16]


def build_md_migration_plan(inventory, *, include_references=True, completed_duration_ns=None):
    """Select completed thyroid studies plus all required earlier MD segments.

    Dataset identity is taken from the original compound-selection/import runs,
    never from a compound ID prefix or its membership in another library.
    """
    runs = {r['run_id']: r for r in inventory['runs']}
    file_rows = {r['run_id']: r for r in inventory['md_files']}
    datasets = {d.source_name.casefold(): (d.project_name, d.source_name) for d in DATASET_PROJECTS}
    withheld = []
    excluded = []
    if completed_duration_ns is not None:
        completed_duration_ns = float(completed_duration_ns)
        if not math.isfinite(completed_duration_ns) or completed_duration_ns <= 0:
            raise ValueError('Completed study duration must be finite and positive')

    def completed_at_duration(row):
        """Check recorded segment length and prior duration, never a filename."""
        result = row.get('result.json') or {}
        md = result.get('md_result') or {}
        if result.get('success') is False or md.get('success') is False:
            return False
        production = (md.get('equilibration_stats') or {}).get('production') or {}
        segment = production.get('duration_ns')
        if segment is None:
            inp = _inputs(row)
            try:
                segment = float(inp['production_steps']) * float(inp['production_timestep_fs']) / 1e6
            except (KeyError, ValueError, TypeError):
                return False
        try:
            previous = float(_meta(row).get('previous_duration_ns') or 0.0)
            total = float(segment) + previous
        except (ValueError, TypeError):
            return False
        return math.isfinite(total) and math.isclose(total, completed_duration_ns, abs_tol=1e-3)

    def target_id(row):
        m, inp = _meta(row), _inputs(row)
        return m.get('target_run_id') or m.get('structure_run_id') or target_for_prep(
            m.get('md_system_prep_run_id') or inp.get('source_md_system_prep_run_id'))

    def target_for_prep(run_id):
        r = runs.get(run_id, {})
        m = _meta(r)
        return m.get('target_run_id') or m.get('structure_run_id') or m.get('parent_run_id')

    def pdb_for(run_id, visited=None):
        visited = set(visited or ())
        if not run_id or run_id in visited:
            return ''
        visited.add(run_id)
        m = _meta(runs.get(run_id, {}))
        pdb = str(m.get('pdb_id') or '').upper()
        if pdb in THR_TARGETS:
            return pdb
        found = set()
        for key in ('source_target_run_id', 'prepared_target_run_id', 'parent_run_id', 'source_structure_run_id'):
            value = pdb_for(m.get(key), visited)
            if value:
                found.add(value)
        return next(iter(found)) if len(found) == 1 else ''

    def dataset_for(run_id, visited=None):
        visited = set(visited or ())
        if not run_id or run_id in visited:
            return None
        visited.add(run_id)
        r = runs.get(run_id, {})
        m = _meta(r)
        name = str(m.get('dataset_name') or m.get('name') or '').casefold()
        if name in datasets:
            return datasets[name]
        ids = list(m.get('compound_run_ids') or [])
        for key in ('compound_sets', 'compound_artifacts'):
            for artifact in _inputs(r).get(key) or []:
                if isinstance(artifact, dict) and artifact.get('run_id'):
                    ids.append(artifact['run_id'])
        if r.get('group') in {'compound-selection', 'compound-import'}:
            ids += [m.get('source_compound_run_id'), m.get('parent_run_id')]
        found = {value for value in (dataset_for(x, visited) for x in ids if isinstance(x, str)) if value}
        return next(iter(found)) if len(found) == 1 else None

    study_cache = {}

    def describe(production_id):
        row = runs.get(production_id, {})
        tid = target_id(row)
        m, tm = _meta(row), _meta(runs.get(tid, {}))
        engine = str(m.get('md_engine') or _inputs(row).get('md_engine') or '').lower()
        engine_evidence = 'recorded engine'
        if engine not in {'openmm', 'gromacs'}:
            outputs = ((row.get('result.json') or {}).get('md_result') or {}).get('output_files') or {}
            checkpoint = outputs.get('production_checkpoint') or (_inputs(row).get('restart_contract') or {}).get('npt_checkpoint')
            if (str(outputs.get('production_trajectory') or '').endswith('.dcd')
                    and str(checkpoint or '').endswith('.chk')):
                engine, engine_evidence = 'openmm', 'recorded DCD trajectory and OpenMM CHK checkpoint'
        if engine not in {'openmm', 'gromacs'}:
            return None
        key = (tid, engine)
        if key in study_cache:
            return study_cache[key]
        pdb = pdb_for(tid) or str(m.get('pdb_id') or '').upper()
        if pdb not in THR_TARGETS:
            return None
        selected = runs.get(tid, {}).get('group') == 'selected-complexes'
        dataset = dataset_for(tm.get('source_structure_run_id') or tm.get('parent_run_id')) if selected else None
        if selected and (not dataset or not tm.get('compound_id') or not tm.get('ligand_smiles')):
            return None
        if not selected:
            if not include_references:
                return None
            dataset = (REFERENCE_PROJECT, 'MD reference controls')
        project, source_dataset = dataset
        pose = str(tm.get('source_pose_id') or '')
        study = {
            'study_id': 'study-' + _digest(str(tid) + ':' + engine),
            'project_name': project, 'dataset_name': source_dataset,
            'compound_id': str(tm.get('compound_id') or ''),
            'compound_name': str(tm.get('compound_name') or ''),
            'ligand_smiles': str(tm.get('ligand_smiles') or ''),
            'identity_status': 'source_metadata_only' if selected else 'reference_identity_unverified',
            'target_pdb_id': pdb, 'target_receptor': THR_TARGETS[pdb]['receptor'],
            'md_engine': engine, 'selected_complex_run_id': str(tid),
            'engine_evidence': engine_evidence,
            'prepared_target_run_id': tm.get('prepared_target_run_id') or str(tid),
            'source_prediction_run_id': tm.get('source_structure_run_id') or tm.get('parent_run_id'),
            'source_pose_id': pose, 'source_pose_rank': int(pose.rsplit('__pose_', 1)[1]) if '__pose_' in pose and pose.rsplit('__pose_', 1)[1].isdigit() else None,
            'source_pose_sdf_index': tm.get('source_pose_sdf_index'),
            'source_tool': tm.get('engine') or tm.get('tool') or 'crystal reference',
            'selection_status': tm.get('selection_status') or 'native/reference input',
            'selection_origin': tm.get('selection_origin'),
            'manual_override_reasons': tm.get('manual_override_reasons'),
            'source_complex_dataset_name': tm.get('complex_dataset_name'),
            'is_best_pose': None,  # Selection rank does not establish energy/confidence rank.
            'production_runs': [], 'mmgbsa_runs': [], 'analysis_runs': [], 'preparation_runs': [],
        }
        study_cache[key] = study
        return study

    required = set()

    def add_production(run_id, stack=()):
        if run_id in stack:
            raise ValueError('Cycle in MD continuation lineage: ' + run_id)
        row = runs.get(run_id)
        if not row or row['group'] != 'bound-ligand-md':
            raise ValueError('Missing MD production dependency: ' + str(run_id))
        study = describe(run_id)
        if not study:
            return None
        if run_id in required:
            return study
        m, inp = _meta(row), _inputs(row)
        previous = m.get('continuation_of_run_id') or inp.get('continuation_of_run_id')
        if previous:
            parent = add_production(previous, (*stack, run_id))
            if not parent or parent['study_id'] != study['study_id']:
                raise ValueError('Continuation crosses compound/target/engine identity: ' + run_id)
        prep = m.get('md_system_prep_run_id') or inp.get('source_md_system_prep_run_id')
        if prep:
            if prep not in runs or runs[prep]['group'] != 'md-system-prep':
                raise ValueError('Missing preparation dependency: ' + str(prep))
            if prep not in required:
                required.add(prep)
                study['preparation_runs'].append({'run_id': prep, 'status': _meta(runs[prep]).get('status')})
        md = (row.get('result.json') or {}).get('md_result') or {}
        production = (md.get('equilibration_stats') or {}).get('production') or {}
        study['production_runs'].append({
            'run_id': run_id, 'status': m.get('status'), 'repeat_index': m.get('repeat_index'),
            'repeat_group_id': m.get('repeat_group_id'), 'workflow_id': m.get('workflow_id'),
            'continuation_of_run_id': previous, 'target_duration_ns': m.get('target_duration_ns'),
            'previous_duration_ns': m.get('previous_duration_ns'),
            'segment_duration_ns': production.get('duration_ns'),
            'trajectory_appended': production.get('trajectory_appended'),
            'production_steps': inp.get('production_steps'),
            'production_timestep_fs': inp.get('production_timestep_fs'),
            'source_artifacts': inp.get('source_artifacts'),
            'result_path': f'artifacts/legacy/{run_id}/result.json',
        })
        required.add(run_id)
        required.add(study['selected_complex_run_id'])
        return study

    # Start with completed scientific results, then retain required earlier segments.
    qualifying = set()
    for r in runs.values():
        if r['group'] == 'bound-ligand-md' and _meta(r).get('status') == 'completed':
            if completed_duration_ns is not None and not completed_at_duration(r):
                excluded.append({'run_id': r['run_id'], 'reason': 'Not a completed study at the selected duration; included only if required by a qualifying continuation'})
                continue
            if add_production(r['run_id']) is None:
                withheld.append({'run_id': r['run_id'], 'reason': 'Unresolved thyroid/compound/engine lineage'})
            else:
                qualifying.add(r['run_id'])
    for r in runs.values():
        if _meta(r).get('status') != 'completed' or r['group'] not in {'md-mmgbsa', 'md-analysis'}:
            continue
        m, inp = _meta(r), _inputs(r)
        source_ids = [m.get('source_production_run_id') or inp.get('source_production_run_id')] if r['group'] == 'md-mmgbsa' else inp.get('production_run_ids', [])
        if completed_duration_ns is not None and (not source_ids or any(id not in required for id in source_ids)):
            excluded.append({'run_id': r['run_id'], 'reason': 'Endpoint/analysis references production outside the selected completed study chains'})
            continue
        sources = [add_production(id) for id in source_ids if isinstance(id, str)]
        if not sources or any(s is None for s in sources) or len({s['study_id'] for s in sources}) != 1:
            withheld.append({'run_id': r['run_id'], 'reason': 'Analysis source studies cannot be uniquely identified'})
            continue
        study = sources[0]
        entry = {'run_id': r['run_id'], 'status': m.get('status'),
                 'result_path': f'artifacts/legacy/{r["run_id"]}/result.json',
                 'parameters': inp.get('parameters') or {},
                 'superseded_by_run_id': m.get('superseded_by_run_id')}
        if r['group'] == 'md-mmgbsa':
            entry['source_production_run_id'] = source_ids[0]
            study['mmgbsa_runs'].append(entry)
        else:
            entry['source_production_run_ids'] = source_ids
            study['analysis_runs'].append(entry)
        required.add(r['run_id'])

    jobs = {}
    for study in study_cache.values():
        if not study['production_runs']:
            continue
        study['qualifying_production_run_ids'] = [r['run_id'] for r in study['production_runs'] if r['run_id'] in qualifying]
        study['selected_completed_duration_ns'] = completed_duration_ns
        for stage, field, group, tool in (
            ('md', 'production_runs', 'ligand-md', study['md_engine'] + '_md'),
            ('mmgbsa', 'mmgbsa_runs', 'ligand-mmgbsa', 'legacy_mmgbsa'),
            ('analysis', 'analysis_runs', 'md-analysis', 'legacy_md_analysis'),
        ):
            if not study[field]:
                continue
            # One visible retrospective job per dataset and tool; targets remain per study.
            key = (study['project_name'], tool)
            job = jobs.setdefault(key, {'project_name': key[0], 'dataset_name': study['dataset_name'],
                'tool': tool, 'task_group': group, 'run_id': 'mnligand-md-' + _digest(':'.join(key)),
                'stage': stage, 'studies': [], 'source_run_ids': []})
            item = copy.deepcopy(study)
            # Each stage owns its source bytes. Other stages carry provenance only.
            for other in ('production_runs', 'mmgbsa_runs', 'analysis_runs', 'preparation_runs'):
                if other != field and not (stage == 'md' and other == 'preparation_runs'):
                    for entry in item[other]:
                        entry.pop('result_path', None)
            job['studies'].append(item)
            job['source_run_ids'].extend(x['run_id'] for x in study[field])
            if stage == 'md':
                job['source_run_ids'].extend(x['run_id'] for x in study['preparation_runs'])
                job['source_run_ids'].append(study['selected_complex_run_id'])
    for job in jobs.values():
        job['source_run_ids'] = sorted(set(job['source_run_ids']))
        job['sources'] = [{'run_id': id, 'path': runs[id]['path'], 'group': runs[id]['group'],
            'size_bytes': file_rows.get(id, {}).get('bytes'), 'status': _meta(runs[id]).get('status')}
            for id in job['source_run_ids']]
        job['size_bytes'] = sum(s['size_bytes'] or 0 for s in job['sources'])
    policy = ('completed results plus preparation and every continuation predecessor' if completed_duration_ns is None
              else f'completed {completed_duration_ns:g} ns studies plus preparation and every continuation predecessor; endpoints/analyses limited to those chains')
    return {'schema_version': 1, 'source_root': inventory['source_root'], 'source_workstation': 'pc004033',
            'selection_policy': policy, 'completed_duration_ns': completed_duration_ns,
            'jobs': sorted(jobs.values(), key=lambda j: (j['project_name'], j['stage'], j['tool'])),
            'withheld': withheld, 'excluded': [r for r in excluded if r['run_id'] not in required],
            'total_size_bytes': sum(j['size_bytes'] for j in jobs.values())}
