"""Scientific presets for the existing graph; independent of campaign eligibility."""
from __future__ import annotations

MODES = {
    'native_complete': 'Complete Native Workflow',
    'native_comparative': 'Complete Native Workflow + Comparative Evaluation',
    'modular_design': 'Modular Design Pipeline',
}
# Compatibility alias for older callers. New UI choices are resolved from the
# registered structure-prediction nodes by ``cofolder_methods`` below.
PREDICTORS = ('boltz2', 'prediction_protenix_v2')


def _design_spec(kind):
    from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG
    from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import WORKFLOWS
    specification = NODE_CATALOG.get(kind, {})
    workflow = specification.get('workflow')
    registered = WORKFLOWS.get(workflow)
    category = specification.get('category')
    inputs = specification.get('inputs') or []
    outputs = specification.get('outputs') or []
    output_kind = outputs[0].get('type') if outputs else None
    target_connected = any('target' in socket.get('accepts', []) for socket in inputs)
    design_route = category in {'Design engine', 'Design pipeline'} and registered is not None
    complete = design_route and category == 'Design pipeline'
    generator = design_route and category == 'Design engine' and target_connected and output_kind in {'backbone', 'sequence'}
    return specification, registered, output_kind, target_connected, complete, generator


def capability(kind):
    specification, registered, output_kind, target_connected, complete, generator = _design_spec(kind)
    from mn_molecular_design_studio.output_vocabulary import measure_capabilities
    workflow = specification.get('workflow')
    tool = registered.tool if registered else None
    graph_category = specification.get('category')
    sequence = ('native_evaluated' if complete else 'backbone_only' if output_kind == 'backbone'
                else 'sequence_ready' if output_kind == 'sequence' else 'unclassified')
    note = (specification.get('description') or 'No registered design operation.')
    output_stages = (('complex_refolding',) if complete or kind == 'esm_binder' else
                     ('generation.backbone',) if sequence == 'backbone_only' else
                     ('generation.backbone_sequence',))
    return dict(complete_native_workflow=complete, generator_only=generator,
                modular_source=complete or generator, output_kind=output_kind,
                sequence_output=sequence, note=note, workflow=workflow, tool=tool,
                route_registered=registered is not None, graph_category=graph_category,
                target_connected=target_connected,
                design_count_field=next((field['name'] for field in specification.get('fields', [])
                                         if field.get('name') == 'num_designs'), None),
                native_measurements=measure_capabilities(tool, stages=output_stages) if tool else [],
                native_output_stages=list(output_stages),
                native_measurement_contract='Recognized source fields only; actual values and missing status are shown per candidate in Results.' if tool else 'No registered workflow route.',
                native_ranking='Only explicit finite native rank evidence is displayed; unavailable otherwise.')


def methods(mode):
    if mode not in MODES:
        raise ValueError('Unknown scientific workflow mode.')
    from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG
    selected = []
    for kind, specification in NODE_CATALOG.items():
        if specification['category'] not in {'Design engine', 'Design pipeline'}:
            continue
        route = capability(kind)
        if not route['target_connected'] or not route['route_registered']:
            continue
        if route['complete_native_workflow'] or mode == 'modular_design' and route['modular_source']:
            selected.append(kind)
    return selected


def cofolder_methods(task='protein_protein'):
    """Registered graph cofolders that accept the requested molecular task."""
    from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG
    from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import WORKFLOWS
    from mn_molecular_design_studio.engine.workflow_engine.workflows import optimized_prediction

    direct_tasks = {
        'refolding.rfd4_proteina': {'protein', 'protein_protein'},
        'refolding.protenix': {'protein', 'protein_protein'},
        'refolding.openfold3': {'protein', 'protein_protein'},
        'refolding.protenix_cli': {'protein', 'protein_protein'},
        'refolding.boltzgen_fold': {'protein', 'protein_protein'},
        'refolding.af2_initial_guess': {'protein_protein'},
        'refolding.af3': {'protein', 'protein_protein', 'protein_ligand', 'protein_dna', 'protein_rna'},
        'refolding.rf3': {'protein', 'protein_protein', 'protein_ligand', 'protein_dna', 'protein_rna'},
        'refolding.boltz2': {'protein', 'protein_protein', 'protein_ligand', 'protein_dna', 'protein_rna'},
    }
    result = []
    for kind, specification in NODE_CATALOG.items():
        if specification.get('category') != 'Structure prediction':
            continue
        if specification.get('cofolder_selector_visible') is False:
            continue
        workflow = specification.get('workflow')
        if workflow not in WORKFLOWS:
            continue
        tasks = set(specification.get('molecule_tasks') or direct_tasks.get(workflow, ()))
        if workflow == 'refolding.optimized_prediction':
            from mn_molecular_design_studio.pipeline_graph import OPTIMIZED_PREDICTION_NODE_ENGINES
            engine = OPTIMIZED_PREDICTION_NODE_ENGINES.get(kind, kind)
            profile = optimized_prediction.ENGINE_PROFILES.get(engine, {})
            tasks = set(specification.get('molecule_tasks') or profile.get('molecule_tasks', ()))
        if task in tasks and any(
                {'sequence', 'prediction'}.intersection(socket.get('accepts', []))
                for socket in specification.get('inputs', [])):
            result.append(kind)
    return result


def analysis_methods():
    """Registered candidate-level Analysis tiles with a compatible collection input."""
    from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG
    from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import WORKFLOWS
    return [kind for kind, specification in NODE_CATALOG.items()
            if specification.get('category') == 'Analysis'
            and specification.get('workflow') in WORKFLOWS
            and any('candidates' in socket.get('accepts', []) for socket in specification.get('inputs', []))]


def sequence_methods():
    """Registered sequence designers that accept candidate collections."""
    from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG
    from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import WORKFLOWS
    candidates = []
    for kind, specification in NODE_CATALOG.items():
        if specification.get('category') != 'Sequence design' or specification.get('hidden'):
            continue
        workflow = specification.get('workflow')
        accepts = set().union(*(set(socket.get('accepts', [])) for socket in specification.get('inputs', [])))
        output = (specification.get('outputs') or [{}])[0].get('type')
        # AbMPNN depends on an explicit SAbR annotation stage. It remains
        # available on the canvas and is not presented as a drop-in initial
        # designer for arbitrary backbones.
        if workflow in WORKFLOWS and output == 'sequence' and 'candidates' in accepts and kind != 'abmpnn':
            candidates.append(kind)
    return candidates


def sequence_count_field(kind):
    """Return the registered per-parent output-count control for a designer."""
    from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG
    specification = NODE_CATALOG.get(kind)
    if not specification or specification.get('category') != 'Sequence design':
        raise ValueError(f"{kind} is not a registered sequence-design operation.")
    names = {field.get('name') for field in specification.get('fields', [])}
    for name in ('num_sequences', 'sequences', 'num_seq_per_target'):
        if name in names:
            return name
    raise ValueError(f"{specification.get('label', kind)} has no registered sequence-count control.")


def _shared_analysis_routes():
    """Measurements produced by explicit downstream cofold/analysis nodes."""
    from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG, OPTIMIZED_PREDICTION_NODE_ENGINES
    from mn_molecular_design_studio.output_vocabulary import DEFINITIONS, MEASURES, PREFIX, measure_capabilities
    from mn_molecular_design_studio.quality_tools import WORKFLOWS as QUALITY_WORKFLOWS

    routes = []
    for kind, specification in NODE_CATALOG.items():
        category, workflow = specification.get('category'), specification.get('workflow')
        if category == 'Structure prediction' and workflow == 'refolding.optimized_prediction':
            engine = OPTIMIZED_PREDICTION_NODE_ENGINES.get(kind, kind)
            measurements = measure_capabilities(engine, stages=('complex_refolding',))
            # These Studio-computed/model-reported fields are stage-specific,
            # not engine-scoped. Preserve their distinct definitions.
            measurements.extend({
                'field': PREFIX + measure.key, 'key': measure.key, 'label': measure.label,
                'unit': measure.unit, 'scope': measure.scope, 'method': measure.method,
                'direction': measure.direction, 'recognized_aliases': list(measure.aliases),
                'stages': list(measure.stages), 'prerequisites': list(measure.prerequisites),
                'value_domain': {'minimum': measure.minimum, 'maximum': measure.maximum},
                'quality_threshold': 'not set; configure explicitly in Shared scientific filters',
                'availability': 'recognized when the source field is emitted; missing values remain missing',
            } for measure in MEASURES if not measure.engines and measure.stages == ('complex_refolding',))
        elif category == 'Analysis' and workflow:
            aliases = set(QUALITY_WORKFLOWS.get(workflow, {}).get('metrics', []))
            if workflow == 'analysis.ipsae':
                keys = [key for key in DEFINITIONS if key.startswith('interface.ipsae.')]
            elif workflow == 'analysis.rosetta':
                keys = [key for key in DEFINITIONS if key.startswith(('energy.rosetta.', 'interface.rosetta.', 'developability.rosetta.'))]
            else:
                keys = [measure.key for measure in MEASURES if aliases.intersection(measure.aliases)]
            measurements = [{
                'field': PREFIX + key, 'key': key, 'label': DEFINITIONS[key].label,
                'unit': DEFINITIONS[key].unit, 'scope': DEFINITIONS[key].scope,
                'method': DEFINITIONS[key].method, 'direction': DEFINITIONS[key].direction,
                'recognized_aliases': list(DEFINITIONS[key].aliases),
                'quality_threshold': 'not set; configure explicitly in Shared scientific filters',
                'availability': 'measured by this route when its prerequisites are satisfied',
            } for key in dict.fromkeys(keys)]
        else:
            continue
        routes.append({'kind': kind, 'workflow': workflow, 'label': specification['label'],
                       'route_type': 'independent cofold' if category == 'Structure prediction' else 'shared analysis',
                       'measurements': measurements})
    return routes


def compatibility_catalog():
    """Capability view of the actual registered design workflows and graph routes."""
    from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG
    from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import WORKFLOWS, workflow_catalog
    from mn_molecular_design_studio.engine.workflow_engine.workflows.generator_capabilities import generator_capabilities
    from mn_molecular_design_studio.output_vocabulary import measure_capabilities

    # workflow_catalog resolves each registered implementation to catch stale
    # module/function routes without submitting a job or invoking a model.
    registered_rows = workflow_catalog(category='design')
    graph_by_workflow = {}
    for kind, specification in NODE_CATALOG.items():
        if specification.get('category') in {'Design engine', 'Design pipeline'} and specification.get('workflow'):
            graph_by_workflow.setdefault(specification['workflow'], []).append(kind)

    workflows = []
    for public_row in registered_rows:
        workflow = public_row['workflow']
        registered = WORKFLOWS[workflow]
        tiles = graph_by_workflow.get(workflow, [])
        role = next((capability(kind) for kind in tiles
                     if capability(kind)['complete_native_workflow'] or capability(kind)['generator_only']), None)
        if role and role['complete_native_workflow']:
            route_status = 'complete_native_graph'
        elif role and role['generator_only']:
            route_status = 'modular_graph_generator'
        elif tiles:
            route_status = 'specialized_graph_route'
        else:
            route_status = 'registered_standalone_route'
        try:
            campaign = generator_capabilities(registered.tool).as_dict()
            campaign_status = 'verified campaign adapter' if campaign['campaign_adapter'] else 'campaign adapter deferred'
            campaign_reason = campaign['deferred_reason']
        except ValueError:
            campaign_status, campaign_reason = 'no campaign capability record', ''
        workflows.append({
            'workflow': workflow, 'tool': registered.tool, 'label': public_row['description'],
            'route_status': route_status, 'graph_tiles': tiles,
            'native_measurements': role['native_measurements'] if role else measure_capabilities(registered.tool),
            'campaign_status': campaign_status, 'campaign_note': campaign_reason,
            'registration': 'workflow implementation resolves from the installed Studio registry; GPU runtime is not invoked by this catalog',
        })
    return {
        'version': 'studio.design-workflow-capabilities.v1',
        'complete_native_workflows': [row for row in workflows if row['route_status'] == 'complete_native_graph'],
        'composable_generators': [row for row in workflows if row['route_status'] == 'modular_graph_generator'],
        'specialized_graph_routes': [row for row in workflows if row['route_status'] == 'specialized_graph_route'],
        'standalone_registered_routes': [row for row in workflows if row['route_status'] == 'registered_standalone_route'],
        'shared_analysis_routes': _shared_analysis_routes(),
        'workflows': workflows,
    }


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
            details = capability(kind)
            item.update({key: details[key] for key in (
                 'complete_native_workflow', 'generator_only', 'modular_source', 'sequence_output', 'output_kind', 'note',
                 'workflow', 'tool', 'route_registered', 'target_connected',
            )})
            item['native_measurement_count'] = len(details['native_measurements'])
            operations.append(item)
        elif NODE_CATALOG[kind]['category'] in {'Sequence design', 'Structure prediction', 'Collection'}:
            stages.append(item)
    complete = bool(operations) and all(o['complete_native_workflow'] for o in operations)
    sequence = any(NODE_CATALOG[s['kind']]['category'] == 'Sequence design' for s in stages)
    prediction = any(NODE_CATALOG[s['kind']]['category'] == 'Structure prediction' for s in stages)
    comparative = any(s['kind'] == 'candidate_analyzer' for s in stages)
    effective = ('native_complete' if complete and not sequence and not prediction and not comparative else
                 'native_comparative' if complete and not sequence and prediction and comparative else
                 'modular_design' if operations and all(o['modular_source'] for o in operations) and comparative else 'custom')
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
