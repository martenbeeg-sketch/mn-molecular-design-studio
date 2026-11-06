"""Results entry into the existing Evaluation Board, canvas and Inspector.

Exploratory policies and manual selections are session state. Source files and
Database records are never edited by this workspace.
"""
from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path

import streamlit as st


def open_results_board(run_dir, candidate_ids=None):
    run = Path(run_dir)
    st.session_state['results-analysis-source'] = {
        'reference': f'{run.parent.name}/{run.name}', 'candidate_ids': candidate_ids or []}
    st.session_state['complex-analysis-workspace-view'] = 'results_board'
    st.session_state['page'] = 'Complex Analysis'
    st.rerun()


def _first_present(*values):
    for value in values:
        if value is not None and str(value).strip() != "":
            return value
    return ""


def _source_candidate_id(candidate):
    raw = candidate.get('raw_metadata') if isinstance(candidate.get('raw_metadata'), dict) else {}
    route = raw.get('complex_analysis_entry') if isinstance(raw.get('complex_analysis_entry'), dict) else {}
    candidate_set = raw.get('candidate_set') if isinstance(raw.get('candidate_set'), dict) else {}
    source = raw.get('source_candidate') if isinstance(raw.get('source_candidate'), dict) else {}
    return str(_first_present(route.get('native_candidate_id'), candidate_set.get('native_candidate_id'),
                              candidate_set.get('native_backbone_identity'), raw.get('native_candidate_id'),
                              source.get('candidate_id'), candidate.get('candidate_id')))


def _measurement_label(candidate):
    raw = candidate.get('raw_metadata') if isinstance(candidate.get('raw_metadata'), dict) else {}
    measurement = raw.get('measurement_status') if isinstance(raw.get('measurement_status'), dict) else {}
    status = str(measurement.get('status') or 'not_measured')
    labels = {
        'measured': 'Measured by a connected tool',
        'recorded': 'Measurements recorded',
        'not_measured': 'No independent measurement recorded',
        'missing': 'Required measurement missing',
        'invalid': 'Measurement invalid',
        'error': 'Measurement failed',
    }
    return labels.get(status, status.replace('_', ' ').capitalize())


def _candidate_measurement_label(candidate, linked_entries=(), *, digest='', model_id=''):
    """Include an exact linked Results analysis when presenting candidate evidence."""
    candidate_id = str(candidate.get('candidate_id') or '')
    source_id = _source_candidate_id(candidate)
    raw = candidate.get('raw_metadata') if isinstance(candidate.get('raw_metadata'), dict) else {}
    model = str(_first_present(candidate.get('model_id'), candidate.get('model_index'),
                               raw.get('model_id'), raw.get('model_index')))
    for entry in linked_entries:
        if not entry.get('metrics'):
            continue
        if str(entry.get('source_candidate_id') or '') not in {candidate_id, source_id}:
            continue
        if digest and str(entry.get('source_structure_digest') or '') != str(digest):
            continue
        if model_id and str(entry.get('source_model_id') or '') != str(model_id):
            continue
        if model and str(entry.get('source_model_id') or '') != model:
            continue
        return 'Measured by a connected tool'
    return _measurement_label(candidate)


def board_entries(candidate, *, label, run_id, digest, policy=None, manual_selection=None):
    """Read saved Results measurements without inventing Database Evaluation IDs."""
    from mn_molecular_design_studio.native_filter_metrics import native_filter_evidence
    raw = candidate.get('raw_metadata') or {}
    route = raw.get('complex_analysis_entry') or {}
    native = native_filter_evidence(candidate)
    saved = raw.get('validation_acceptance') or route.get('acceptance') or {}
    context = {'candidate_id': candidate['candidate_id'], 'source_structure_digest': digest,
               'candidate': candidate, 'entry': route, 'acceptance': saved}
    source_candidate_id = _source_candidate_id(candidate)
    measurement = raw.get('measurement_status') if isinstance(raw.get('measurement_status'), dict) else {}
    measurement_status = str(measurement.get('status') or 'not_measured')
    source_model_id = _first_present(route.get('model_id'), candidate.get('model_id'), candidate.get('model_index'),
                                     raw.get('model_id'), raw.get('model_index'))
    entry = {
        'molecule_id': '', 'structure_id': '', 'structure_label': label,
        'evaluation_id': raw.get('evaluation_id', ''), 'observation_id': raw.get('observation_id', ''),
        'method': candidate.get('source_tool') or 'Saved Results', 'category': 'Saved Results measurements',
        'source_engine': candidate.get('source_tool') or 'Unknown', 'run_id': run_id,
        'candidate_id': candidate['candidate_id'], 'structure_digest': digest,
        'parameters': raw.get('analysis_parameters') or {}, 'provenance': context,
        'metrics': candidate.get('metrics') or {}, 'created_at': '',
        'status': measurement_status,
        'measurement_label': _measurement_label(candidate),
        'native_acceptance': native.get('status', 'unavailable'),
        'exploratory_acceptance': (policy or saved).get('status', 'not_evaluated'),
        'manual_selection': ((raw.get('human_review') or {}).get('selected', False)
                             if manual_selection is None else bool(manual_selection)),
        'source_candidate_id': str(source_candidate_id),
        'source_model_id': str(source_model_id),
        'source_structure_digest': str(route.get('structure_sha256') or digest),
        'source_run_id': str(route.get('source_run') or run_id),
    }
    return [entry]


def canvas_graph(job, candidate_ids, policy, *, entry_mode='accepted_only', manual_ids=None, analysis='biotite_interface'):
    """Use ordinary registered source, policy, handoff and analysis nodes."""
    from mn_molecular_design_studio.pipeline_graph import new_graph, NODE_CATALOG
    graph = new_graph('accepted_complex_analysis')
    nodes = {node['id']: node for node in graph['nodes']}
    nodes['project_design_source']['data']['params'].update(
        source_run_dir=str(job['run_dir']), selected_candidate_ids=candidate_ids, source_scope='complex_results')
    if entry_mode == 'accepted_only':
        nodes['shared_filter']['data']['params'].update(
            policy_profile='custom', binder_format=policy['binder_format'], target_type=policy['target_type'],
            rules=[{key: rule[key] for key in ('field', 'operator', 'value')} for rule in policy['mandatory_filters']])
    else:
        nodes.pop('shared_filter')
    nodes['complex_analysis']['data']['params'].update(entry_mode=entry_mode, selected_candidate_ids=manual_ids or [])
    old = nodes.pop('biotite_interface')
    if analysis:
        spec = NODE_CATALOG[analysis]
        old['id'] = analysis
        old['data'] = {'kind': analysis, 'label': spec['label'],
                       'params': {item['name']: copy.deepcopy(item['default']) for item in spec['fields']}}
        nodes[analysis] = old
    order = ['project_design_source', 'candidate_merge'] + (['shared_filter'] if entry_mode == 'accepted_only' else []) + ['complex_analysis'] + ([analysis] if analysis else []) + ['results']
    graph['nodes'] = [nodes[identifier] for identifier in order]
    graph['edges'] = [{'id': f'handoff-{index}', 'source': a, 'target': b,
                      'sourceHandle': 'out', 'targetHandle': 'in'} for index, (a, b) in enumerate(zip(order, order[1:]))]
    graph['name'] = f"{job.get('campaign_name') or job['run_id']} · exploratory complex analysis"
    return graph


def _linked_measurements(user, job, selected):
    """Join exact source candidate, model, digest and role-qualified receipts."""
    from mn_molecular_design_studio import complex_analysis_ui as workspace
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json
    aliases: dict[str, list[dict]] = {}
    for item in selected:
        candidate = item['candidate']
        raw = candidate.get('raw_metadata') or {}
        candidate_set = raw.get('candidate_set') or {}
        route = raw.get('complex_analysis_entry') or {}
        for value in (candidate.get('candidate_id'), candidate_set.get('native_candidate_id'),
                      route.get('native_candidate_id'), route.get('candidate_id')):
            if value:
                aliases.setdefault(str(value), []).append(item)
    samples = {identifier: matches[0] for identifier, matches in aliases.items() if len(matches) == 1}
    result = []
    for analysis in workspace._accessible_result_jobs(user):
        if analysis.get('task_group') != 'analysis':
            continue
        run = Path(analysis['run_dir'])
        path = run / 'artifacts/normalized_candidates/candidates.jsonl'
        if not path.is_file():
            continue
        metadata = read_json(run / 'metadata.json')
        # Candidate-aware graph analyses carry an exact source receipt. Direct
        # contextual analyses carry the original candidate in source_candidate.
        for candidate in read_candidates(path):
            raw = candidate.get('raw_metadata') or {}
            entry = raw.get('complex_analysis_entry') or {}
            original = raw.get('source_candidate') or {}
            source_id = entry.get('native_candidate_id') or original.get('candidate_id')
            source = samples.get(str(source_id or ''))
            if source is None:
                for alias in (entry.get('candidate_id'), entry.get('structure_candidate_id'),
                              original.get('candidate_id'), candidate.get('candidate_id')):
                    if alias and str(alias) in samples:
                        source = samples[str(alias)]
                        break
            if not source:
                continue
            saved_digest = entry.get('structure_sha256') or raw.get('structure_digest')
            if saved_digest != source['digest']:
                continue
            if entry:
                source_run = Path(entry.get('source_run_dir') or '').name or str(entry.get('source_run') or '')
                if source_run != job['run_id']:
                    continue
            elif str(original.get('candidate_id') or '') != str(source['candidate'].get('candidate_id') or ''):
                # Some registered analysis tools preserve their exact source
                # candidate but only expose the intermediate handoff as the
                # analysis job's upstream_run_id. The embedded source ID,
                # structure digest, chains and model below identify the exact
                # Results candidate in that case.
                continue
            if list(candidate.get('target_chains') or []) != list(source['candidate'].get('target_chains') or []):
                continue
            if list(candidate.get('binder_chains') or []) != list(source['candidate'].get('binder_chains') or []):
                continue
            source_candidate = source['candidate']
            source_raw = source_candidate.get('raw_metadata') or {}
            analysis_model = _first_present(entry.get('model_id'), candidate.get('model_id'), candidate.get('model_index'),
                                            raw.get('model_id'), raw.get('model_index'))
            source_model = _first_present(source_candidate.get('model_id'), source_candidate.get('model_index'),
                                          source_raw.get('model_id'), source_raw.get('model_index'))
            if str(analysis_model) != str(source_model):
                continue
            result.extend(board_entries(candidate, label=source['label'], run_id=analysis['run_id'], digest=saved_digest))
    return result


def render_results_board(user, *, workflow_schema=None, submit=None):
    from mn_molecular_design_studio import complex_analysis_ui as workspace
    from mn_molecular_design_studio.validation_acceptance import apply_acceptance, infer_target_type
    from mn_molecular_design_studio.native_filter_metrics import native_filter_evidence
    from mn_molecular_design_studio.complex_analysis_collection_ui import _render_evaluation_board
    from mn_molecular_design_studio.output_vocabulary import candidate_output, PREFIX

    st.title('Results → Complex Analysis')
    st.caption('Native acceptance, exploratory policy and manual selection remain separate. This workspace reads Results; it does not promote molecules into Database.')
    context = st.session_state.get('results-analysis-source') or {}
    jobs = workspace._accessible_result_jobs(user)
    by_ref = {f"{job['task_group']}/{job['run_id']}": job for job in jobs}
    if not by_ref:
        st.info('No accessible completed Results are available.')
        return
    reference = st.selectbox('Source Results campaign', list(by_ref),
        index=list(by_ref).index(context['reference']) if context.get('reference') in by_ref else 0,
        format_func=lambda ref: f"{by_ref[ref].get('campaign_name') or by_ref[ref].get('tool')} · {ref}", key='results-board-run')
    job = by_ref[reference]
    board_states = st.session_state.setdefault('results-board-state', {})
    saved_state = board_states.get(reference, {})
    options = [item for item in workspace._job_structure_options(user, job) if item.get('candidate_id')]
    by_id = {item['candidate_id']: item for item in options}
    if not by_id:
        st.info('This Results job has no local saved candidate structures. Open it on its source workstation to inspect remote files.')
        return
    wanted = [identifier for identifier in context.get('candidate_ids', []) if identifier in by_id] if context.get('reference') == reference else []
    subset_key = 'results-board-subset-'+job['run_id']
    if subset_key not in st.session_state:
        st.session_state[subset_key] = [item for item in saved_state.get('subset', wanted or list(by_id)[:8]) if item in by_id]
    candidate_file = Path(job['run_dir']) / 'artifacts/normalized_candidates/candidates.jsonl'
    try:
        from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
        saved_candidate_count = len(read_candidates(candidate_file))
    except (OSError, ValueError, TypeError):
        saved_candidate_count = len(options)
    st.caption(f"{len(options)} candidate structures are locally available for this Results source (of {saved_candidate_count} saved candidates).")
    subset = st.multiselect('Candidates to review', list(by_id), key=subset_key)
    if not subset:
        st.info('Select source candidates to inspect or compare.')
        return
    originals = [by_id[identifier]['candidate'] for identifier in subset]
    recorded_policy = next(((row.get('raw_metadata') or {}).get('validation_acceptance', {}).get('policy') for row in originals
                            if (row.get('raw_metadata') or {}).get('validation_acceptance', {}).get('policy')), None)
    numeric = sorted({key for row in originals for key, value in (row.get('metrics') or {}).items() if type(value) in {int, float}}
                     | {PREFIX+key for row in originals for key, value in candidate_output(row).get('metrics', {}).items() if type(value) in {int, float}})
    with st.expander('Exploratory acceptance criteria', expanded=True):
        st.caption('These criteria are an exploratory view of the saved measurements. The original generator decisions remain unchanged. Missing required measurements are unresolved.')
        count_key = 'results-board-rule-count-'+job['run_id']
        stored_rules = saved_state.get('rules')
        if count_key not in st.session_state:
            st.session_state[count_key] = len(stored_rules) if stored_rules is not None else len(recorded_policy.get('mandatory_filters', [])) if recorded_policy else 0
        count = int(st.number_input('Number of mandatory criteria', min_value=0, max_value=8, key=count_key))
        rules = []
        for index in range(count):
            saved = (recorded_policy.get('mandatory_filters', []) if recorded_policy else [])
            rule = (stored_rules[index] if stored_rules is not None and index < len(stored_rules)
                    else saved[index] if index < len(saved) else {})
            previous_fields = {item['field'] for item in rules}
            fields = [key for key in numeric if key not in previous_fields]
            if not fields:
                st.info('No numeric measurements are recorded for these candidates.')
                return
            first, second, third = st.columns([3, 1, 1])
            field_key = f'results-board-field-{job["run_id"]}-{index}'
            preferred = rule.get('field') if rule.get('field') in fields else fields[0]
            if field_key not in st.session_state or st.session_state[field_key] not in fields:
                st.session_state[field_key] = preferred
            field = first.selectbox(f'Measurement {index+1}', fields, key=field_key)
            operators = ['>=', '<=', '>', '<']
            op = second.selectbox(f'Comparison {index+1}', operators, index=operators.index(rule.get('operator', '>=')), key=f'results-board-op-{job["run_id"]}-{index}')
            value = third.number_input(f'Threshold {index+1}', value=float(rule.get('value', 0.)), key=f'results-board-value-{job["run_id"]}-{index}')
            rules.append({'field': field, 'operator': op, 'value': value})
        policy = {'binder_format': 'miniprotein', 'target_type': infer_target_type(originals[0]) or 'structured', 'mandatory_filters': rules}
    board_states[reference] = {**saved_state, 'subset': list(subset), 'rules': copy.deepcopy(rules)}
    try:
        evaluated, summary = apply_acceptance(originals, policy)
    except ValueError as exc:
        st.warning(f'Cannot apply this exploratory policy: {exc}')
        return
    source_items = []
    for identifier in subset:
        item = by_id[identifier]
        digest = hashlib.sha256(item['path'].read_bytes()).hexdigest()
        candidate = next(row for row in evaluated if row['candidate_id'] == identifier)
        source_items.append({**item, 'candidate': candidate, 'digest': digest,
                             'label': f"{identifier} · {candidate.get('source_tool')}"})
    linked = _linked_measurements(user, job, source_items)
    st.caption(f"Exploratory policy: {summary['accepted']} pass · {summary['rejected']} rejected · {summary['unresolved']} unresolved · {summary['not_evaluated']} not evaluated")
    st.dataframe([{'Studio candidate': row['candidate_id'],
        'Source candidate': _source_candidate_id(row),
        'Independent measurement': _candidate_measurement_label(
            row, linked, digest=next(item['digest'] for item in source_items if item['candidate']['candidate_id'] == row['candidate_id'])),
        'Native acceptance': native_filter_evidence(row)['status'],
        'Exploratory acceptance': row['raw_metadata']['validation_acceptance']['status'],
        'Reasons': '; '.join(row['raw_metadata']['validation_acceptance']['reasons']),
        **{rule['field']: row['raw_metadata']['validation_acceptance']['measurements'][rule['field']].get('value') for rule in rules}}
        for row in evaluated], hide_index=True, use_container_width=True)
    mode_key = 'results-board-entry-mode'
    if mode_key not in st.session_state:
        st.session_state[mode_key] = saved_state.get('mode', 'Accepted only')
    mode_label = st.radio('Analysis route', ['Accepted only', 'Manual selection'], horizontal=True, key=mode_key)
    manual_key = 'results-board-manual-'+job['run_id']
    if manual_key not in st.session_state:
        st.session_state[manual_key] = [item for item in saved_state.get('manual', []) if item in subset]
    manual = st.multiselect('Explicit manual selection', subset, key=manual_key)
    board_states[reference].update(mode=mode_label, manual=list(manual))
    if mode_label == 'Accepted only':
        chosen = [row['candidate_id'] for row in evaluated if row['raw_metadata']['validation_acceptance']['status'] == 'accepted'
                  and native_filter_evidence(row)['status'] not in {'failed', 'rejected'}]
    else:
        chosen = manual
    st.caption(f'{len(chosen)} candidate(s) enter this analysis route. Manual selection does not change acceptance.')
    if not chosen:
        st.info('Configure required measurements to obtain passing candidates, or explicitly select candidates in Manual selection mode.')
        return
    selected = [item for item in source_items if item['candidate']['candidate_id'] in chosen]
    samples = []
    for item in selected:
        identifier = item['candidate']['candidate_id']
        candidate = item['candidate']
        digest = item['digest']
        label = item['label']
        samples.append({'evaluation_entries': board_entries(candidate, label=label, run_id=job['run_id'], digest=digest,
                                                           policy=candidate['raw_metadata']['validation_acceptance'],
                                                           manual_selection=(mode_label == 'Manual selection' and identifier in manual))})
    if linked:
        samples.append({'evaluation_entries': linked})
    _render_evaluation_board(samples)
    comparable = [item for item in selected if item['candidate'].get('binder_chains') and item['candidate'].get('target_chains')]
    if len(comparable) >= 2:
        with st.expander('Compare selected structures and interface maps', expanded=False):
            from mn_molecular_design_studio.complex_analysis_comparison import protein_chain_numbering, target_structure_fingerprint
            from mn_molecular_design_studio.complex_analysis_comparison_ui import render_comparison
            compare_key = 'results-board-compare-'+job['run_id']
            compare_ids = [item['candidate']['candidate_id'] for item in comparable]
            if compare_key not in st.session_state:
                st.session_state[compare_key] = compare_ids[:2]
            selected_compare_ids = st.multiselect('Structures for detailed comparison', compare_ids,
                key=compare_key, max_selections=5,
                help='Choose two to five structures. The Evaluation Board above continues to show every selected candidate.')
            selected_for_comparison = [item for item in comparable if item['candidate']['candidate_id'] in selected_compare_ids]
            compared = []
            for item in selected_for_comparison:
                candidate = item['candidate']
                target = candidate.get('target_chains') or []
                binder = candidate.get('binder_chains') or []
                target_numbering = protein_chain_numbering(item['path'], target)
                binder_numbering = protein_chain_numbering(item['path'], binder)
                compared.append({'kind': 'results', 'identity': reference + ':' + item['relative_path'],
                    'label': item['label'], 'job': job, 'candidate': candidate, 'path': item['path'],
                    'candidate_id': candidate['candidate_id'], 'project_id': str(job.get('project_id') or ''),
                    'components': workspace.component_inventory(item['path'].read_bytes(), item['path'].name),
                    'model_id': str(_first_present(candidate.get('model_id'), candidate.get('model_index'),
                                                   (candidate.get('raw_metadata') or {}).get('model_id'),
                                                   (candidate.get('raw_metadata') or {}).get('model_index'))),
                    'structure_digest': item['digest'], 'target_structure_fingerprint': target_structure_fingerprint(item['path'], target),
                    'target_chains': target, 'binder_chains': binder, 'target_numbering': target_numbering,
                    'target_chains_sequence': [row['sequence'] for row in target_numbering],
                    'binder_chains_sequence': [row['sequence'] for row in binder_numbering], 'complex_type': 'protein-protein',
                    'context': candidate.get('raw_metadata') or {}, 'viewer_key': 'results-board-compare-'+candidate['candidate_id'],
                    'evaluation_entries': [
                        *board_entries(candidate, label=item['label'], run_id=job['run_id'], digest=item['digest']),
                        *[entry for entry in linked
                          if entry.get('source_candidate_id') == candidate['candidate_id']
                          and entry.get('source_structure_digest') == item['digest']
                          and entry.get('source_model_id') == str(_first_present(
                              candidate.get('model_id'), candidate.get('model_index'),
                              (candidate.get('raw_metadata') or {}).get('model_id'),
                              (candidate.get('raw_metadata') or {}).get('model_index')))],
                    ]})
            if len(compared) >= 2:
                render_comparison(user, workflow_schema=workflow_schema, submit=submit, selected_sources=compared)
    first, second = st.columns(2)
    with first:
        inspect_id = st.selectbox('Candidate to inspect', chosen, key='results-board-inspect')
        if st.button('Open in Scientific Inspector', key='results-board-inspect-open'):
            item = by_id[inspect_id]
            workspace.open_result_complex(job['run_dir'], item['path'], item['candidate'])
    with second:
        ppi = all(by_id[identifier]['candidate'].get('binder_chains') and by_id[identifier]['candidate'].get('binder_sequence') for identifier in chosen)
        analysis = st.selectbox('Connected analysis tool', ['biotite_interface', 'interface_clustering', 'ipsae', 'rosetta_metrics'] if ppi else ['plip', 'pandamap', 'posebusters'],
            format_func=lambda key: __import__('mn_molecular_design_studio.pipeline_graph', fromlist=['NODE_CATALOG']).NODE_CATALOG[key]['label'], key='results-board-tool')
        if st.button('Continue on canvas', key='results-board-canvas'):
            graph = canvas_graph(job, subset, policy, entry_mode='accepted_only' if mode_label == 'Accepted only' else 'manual_selection', manual_ids=manual, analysis=analysis)
            project = job.get('project_id')
            prefix = f'visual-pipeline-{project}'
            st.session_state[prefix+'-graph'] = graph
            st.session_state[prefix+'-reset'] = graph['id']
            st.session_state[prefix+'-name'] = graph['name']
            st.session_state.setdefault('project-pending-selection', {})['visual-pipeline-project-id'] = project
            st.session_state['pipeline-workspace-view'] = 'Node canvas'
            st.session_state['page'] = 'Pipelines'
            st.rerun()
    st.markdown('#### Selected structures')
    for identifier in chosen[:5]:
        item = by_id[identifier]
        with st.expander(identifier, expanded=len(chosen) == 1):
            workspace._show_structure(item['path'], item['candidate'].get('raw_metadata') or {},
                workspace.component_inventory(item['path'].read_bytes(), item['path'].name), key='results-board-view-'+identifier)
            st.json({'candidate_id': identifier, 'run_id': job['run_id'], 'structure': item['relative_path'],
                     'structure_sha256': hashlib.sha256(item['path'].read_bytes()).hexdigest(),
                     'source_metadata': item['candidate'].get('raw_metadata') or {}})
    if st.button('Return to exact source Results', key='results-board-return'):
        st.session_state['selected-job-reference'] = reference
        st.session_state['page'] = 'Job details'
        st.rerun()
