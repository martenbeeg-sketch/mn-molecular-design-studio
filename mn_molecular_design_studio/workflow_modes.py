"""Scientific presets for the existing graph; independent of campaign eligibility."""
from __future__ import annotations

MODES = {
    'native_complete': 'Complete Native Workflow',
    'native_comparative': 'Complete Native Workflow + Comparative Evaluation',
    'modular_design': 'Modular Design Pipeline',
}
# Reviewed ordinary adapter contracts, including hidden saved native forms.
COMPLETE = {
    'bindcraft2': 'Native design, MPNN, validation and filters; final rank only when recorded.',
    'boltzgen': 'Native generation, inverse folding, refolding, filtering and selection.',
    'pxdesign_pipeline': 'Native ProteinMPNN and AF2-IG; extended preset includes Protenix.',
    'genie3_pipeline': 'Native inverse folding and structure evaluation.',
    'protpardelle_pipeline': 'Native MPNN and binder ESMFold evaluation; no independent interface confidence implied.',
    'complexa': 'Native generation, filtering and evaluation; chosen-sequence evidence can be unresolved.',
}
GENERATORS = {
    'rfd4': ('co_design', 'Guided co-design; refolding/evaluation depends on task options. No validated final native ranking.'),
    'rfdiffusion': ('backbone_only', 'Backbone generation requires initial sequence design.'),
    'rfd3': ('backbone_only', 'RFdiffusion3 generation requires initial sequence design.'),
    'genie3': ('backbone_only', 'Generation only; native evaluated pipeline is a separate operation.'),
    'pxdesign': ('backbone_only', 'Generation only; native evaluated pipeline is a separate operation.'),
    'protpardelle': ('backbone_only', 'Generation only; native evaluated pipeline is a separate operation.'),
    'bindcraft2_generator': ('sequence_ready', 'Gradient design through the pre-MPNN gate only.'),
    'boltzgen_generator': ('sequence_ready', 'Native design stage only.'),
    'complexa_generator': ('sequence_ready', 'Native generation only, with explicit target/binder roles.'),
    'esm_binder': ('sequence_ready', 'Gradient sequence design and native complex; no complete native evaluation workflow.'),
    'rfdiffusion_beta_pairing': ('backbone_only', 'Requires compatible strand conditioning and initial sequence design.'),
    'kejia_idr': ('backbone_only', 'Requires its existing template preparation and initial sequence design.'),
}
PREDICTORS = ('boltz2', 'prediction_protenix_v2')


def capability(kind):
    complete = kind in COMPLETE
    sequence, note = GENERATORS.get(kind, ('native_evaluated' if complete else 'unclassified', COMPLETE.get(kind, 'No reviewed binder-mode operation.')))
    return dict(complete_native_workflow=complete, generator_only=kind in GENERATORS,
                sequence_output=sequence, note=note,
                native_ranking='Only explicit finite native rank evidence is displayed; unavailable otherwise.')


def methods(mode):
    if mode not in MODES:
        raise ValueError('Unknown scientific workflow mode.')
    from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG
    choices = GENERATORS if mode == 'modular_design' else COMPLETE
    return [kind for kind in choices if kind in NODE_CATALOG
            and any('target' in socket['accepts'] for socket in NODE_CATALOG[kind]['inputs'])]


def describe(graph):
    """Actual enabled nodes after editing, separately from the starting preset."""
    from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG
    operations, stages = [], []
    for node in graph['nodes']:
        data = node['data']
        if data.get('enabled') is False:
            continue
        kind = data['kind']
        item = dict(node=node['id'], kind=kind, label=data['label'], retained=bool(data.get('retained_run_dir')))
        if NODE_CATALOG[kind]['category'] in {'Design engine', 'Design pipeline'}:
            operations.append(item | capability(kind))
        elif NODE_CATALOG[kind]['category'] in {'Sequence design', 'Structure prediction', 'Collection'}:
            stages.append(item)
    complete = bool(operations) and all(o['complete_native_workflow'] for o in operations)
    sequence = any(NODE_CATALOG[s['kind']]['category'] == 'Sequence design' for s in stages)
    prediction = any(NODE_CATALOG[s['kind']]['category'] == 'Structure prediction' for s in stages)
    comparative = any(s['kind'] == 'candidate_analyzer' for s in stages)
    effective = ('native_complete' if complete and not sequence and not prediction and not comparative else
                 'native_comparative' if complete and not sequence and prediction and comparative else
                 'modular_design' if operations and all(o['generator_only'] for o in operations) and comparative else 'custom')
    preset = graph.get('workflow_mode')
    return dict(preset_mode=preset, preset_label=MODES.get(preset), effective_mode=effective,
                effective_label=MODES.get(effective, 'Custom edited graph'), generator_operations=operations,
                downstream_stages=stages,
                ranking='Explicit measurements within recorded comparable protocols' if any(s['kind'] == 'candidate_rank' for s in stages)
                        else 'Native evidence only; no common numerical rank')


def native_ranking(row):
    """Never synthesize a native order from row position or a generic score."""
    from mn_molecular_design_studio.measurement_status import finite_number
    metrics = row.get('metrics', {})
    for name in ('native_final_rank', 'native_rank', 'bindcraft2_rank'):
        value = finite_number(metrics.get(name))
        if value is not None and value >= 1:
            return dict(status='recorded_native_order', field=name, rank=value,
                        metric=metrics.get('native_rank_metric'), value=metrics.get('native_rank_value'),
                        scope='Within original engine and source run only; recorded order is not independent validation')
    return dict(status='unavailable', rank=None, scope='No validated final native ranking recorded')


def observed_execution_plan(graph, node_states):
    from mn_molecular_design_studio.candidate_pipeline import execution_plan
    inputs = {}
    for node in graph['nodes']:
        parents = [edge['source'] for edge in graph['edges'] if edge['target'] == node['id']]
        if len(parents) == 1:
            inputs[node['id']] = node_states.get(parents[0], {}).get('output', {}).get('candidate_count')
    plan = execution_plan(graph, inputs)
    for row in plan:
        state = node_states.get(row['node'], {})
        row.update(status=state.get('status', 'pending'), actual_output_candidates=state.get('candidate_count'),
                   actual_processing_jobs=len(state['batches']) if 'batches' in state else int(bool(state.get('job_run_dir'))))
        row['note'] += ' Forecast uses upstream counts; conditional sequence eligibility can reduce work. Actual jobs are reported separately.'
    return plan


def gpu_monitor_rows(nodes):
    """Heartbeat/scheduler availability, never a claim of per-model eligibility."""
    rows = []
    for node in nodes:
        resources = node.get('snapshot', {}).get('resources', {})
        free = {str(g['gpu_id']) for g in resources.get('free_gpus', [])}
        for gpu in resources.get('gpus', []):
            device = str(gpu['gpu_id'])
            rows.append({'Workstation': node.get('name', node['id']), 'GPU': device, 'Model': gpu.get('name'),
                         'State': 'unknown / workstation ' + node.get('availability', 'unknown') if node.get('availability') != 'online' else
                                  'available' if device in free else 'active / reserved',
                         'Free VRAM GiB': gpu.get('free_vram_gb'), 'Snapshot': node.get('last_seen_at', 'local scheduler')})
    return rows
