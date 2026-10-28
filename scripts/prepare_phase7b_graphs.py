"""Save bounded PD-L1 acceptance graphs through Studio's normal graph API.

This prepares definitions only. It never submits a job or invokes a native tool.
Live runs must use the regular Studio graph submission and scheduler.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from mn_molecular_design_studio import candidate_pipeline as cp, pipeline_graph as pg
from mn_molecular_design_studio.auth import development_user


def prepare(evidence: Path):
    from mn_molecular_design_studio.workspace import ensure_studio_workspace
    ensure_studio_workspace()
    results = Path('/mnt/data/RESULTS/mn-molecular-design-studio')
    source = results / 'workdir/runs/design/20261028-124448-350f7294'
    native_source = results / 'workdir/runs/design/20261028-085821-36356a14'
    target = source / 'artifacts/rfd4_proteina/input/PD-L1.pdb.cif'
    original = json.loads((native_source / 'input.json').read_text())['params']
    retained = json.loads((source / 'input.json').read_text())['params']
    marker = json.loads((source / 'molecular_design_studio.json').read_text())
    owner, project = development_user()['id'], marker['project_id']
    definitions = []
    for mode in ('native_complete', 'native_comparative', 'modular_design'):
        graph = cp.new_pipeline(['rfd4', 'rfd4'] if mode == 'modular_design' else ['bindcraft2'], mode=mode)
        labels = {'native_complete': 'Mode 1 native', 'native_comparative': 'Mode 2 retained native comparison',
                  'modular_design': 'Mode 3 modular + retained P5 collection'}
        graph['name'] = 'Phase 7B PD-L1 · ' + labels[mode]
        for node in graph['nodes']:
            data, nid = node['data'], node['id']
            if nid == 'target':
                data['params'].update(target_pdb=str(target), target_chains='A',
                                      hotspots='A40,A99,A107' if mode == 'modular_design' else '')
            if nid == 'generator-0' and mode != 'modular_design':
                # Preserve the retained native smoke-test scientific settings;
                # its one-trajectory cap is stricter than Phase 7B's maximum two.
                data['params'].update(num_designs=1, binder_length=original['binder_length'],
                    campaign_seed=original['campaign_seed'], engine_overrides_json=json.dumps(original['engine_overrides']),
                    workers_per_gpu='1', max_trajectories=1)
                data['execution'].update(cpu_cores=4, min_vram_gb=8, retries=0, gpu_device='auto')
                if mode == 'native_comparative':
                    snapshot = json.loads((evidence.parent / 'retained/native-snapshot.json').read_text())
                    data['retained_run_dir'] = snapshot['run_dir']
            elif nid.startswith('generator-') and mode == 'modular_design':
                data['params'].update(input_mode='guided', num_designs=1,
                    length_min=retained['length_min'], length_max=retained['length_max'], nsteps=retained['nsteps'],
                    seed=retained['seed'],
                    options_json=json.dumps(retained['options']))
                data['execution'].update(retries=0, gpu_device='auto')
                if nid.endswith('-g1'):
                    data['retained_run_dir'] = str(source)
                    data['label'] = 'Retained P5 RFD4 co-design invocation · two outputs'
            if nid.startswith('predictor-'):
                data['enabled'] = True
                data['params'].update(samples=1, seed=42, target_use_msa=False, use_target_template=False)
                data['execution'].update(retries=0, cpu_cores=8, min_vram_gb=16, gpu_device='auto')
                if data['kind'] == 'boltz2':
                    data['params'].update(num_recycles=1, sampling_steps=20, seed_count=1, generate_ipsae=False)
                    data['execution']['batch_size'] = 2
                else:
                    # Two retained inputs give two useful independent jobs while
                    # Boltz demonstrates one actual multi-input invocation.
                    data['execution']['batch_size'] = 1
            if nid == 'ranking':
                data['params']['ranking_rules'] = [{'field': 'prediction/predictor-0/iptm', 'direction': 'higher'}]
        graph = pg.validate_graph(graph, require_inputs=True)['graph']
        cp.preflight(graph, pg.validate_graph(graph, require_inputs=True), project)
        saved = pg.save_definition(graph, owner, project, version=True)
        evidence.mkdir(parents=True, exist_ok=True)
        path = evidence / (mode + '.json')
        path.write_text(json.dumps(saved, indent=2))
        definitions.append(dict(mode=mode, graph_id=graph['id'], name=graph['name'], revision=saved['revision'],
                                owner_id=owner, project_id=project, path=str(path)))
    (evidence / 'definitions.json').write_text(json.dumps(definitions, indent=2))
    print(json.dumps(definitions, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence', type=Path, required=True)
    prepare(parser.parse_args().evidence)
