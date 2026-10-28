"""Replay retained PD-L1 outcomes through the normal graph; no fresh inference.

All generators and enabled predictors explicitly use retained_run_dir. Artifacts
are checksum-verified before and after. This is never live scientific acceptance.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
from zipfile import ZipFile

from mn_molecular_design_studio import candidate_pipeline as cp, candidate_sets as cs, pipeline_graph as pg
from mn_molecular_design_studio.engine.workflow_engine.core.artifacts import build_run_zip
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json
from mn_molecular_design_studio.engine.workflow_engine.workflows.graph_pipeline import run_graph_pipeline


def main(evidence):
    from mn_molecular_design_studio.workspace import ensure_studio_workspace
    ensure_studio_workspace()
    root = Path('/mnt/data/RESULTS/mn-molecular-design-studio/workdir/runs/design')
    native = root / '20261028-085821-36356a14'
    from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path
    p4 = resolve_stored_path('app:///workdir/runs/design/20261028-101057-a07ef0a4', must_exist=True)
    assert p4 is not None
    p5 = root / '20261028-124448-350f7294'
    target = p5 / 'artifacts/rfd4_proteina/input/PD-L1.pdb.cif'
    tracked = [p for directory in (p4, p5, native) for p in directory.rglob('*')
               if p.is_file() and (p.suffix in {'.jsonl', '.pdb', '.cif'} or p.name in {'ledger.json', 'input.json', 'result.json'})]
    tracked.append(target)
    target_context = dict(target_pdb=str(target), target_chains='A', hotspots='A40,A99,A107')
    target_context['identity'] = cs.target_identity(target_context)
    lineage = []
    for directory in (p4, p5):
        # Older P4 runs predate the external Results-root migration. Resolve
        # their run:/// links against the recorded source root, as snapshot()
        # does, rather than rewriting any original scientific artifact.
        os.environ['MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR'] = str(directory.parent.parent)
        for row in cs.iter_candidates(dict(run_dir=str(directory), candidates_jsonl=str(directory / 'artifacts/normalized_candidates/candidates.jsonl'))):
            cp.verify_target(row, target_context)
            evidence_rows = [row, *row.get('raw_metadata', {}).get('evidence_candidates', [])]
            for item in evidence_rows:
                for key in ('complex_pdb', 'binder_pdb', 'target_pdb'):
                    if item.get(key):
                        path = resolve_stored_path(item[key], run_dir=directory, must_exist=True)
                        assert path is not None, (directory, row['candidate_id'], key)
                        tracked.append(path)
            lineage.append(dict(source_run=str(directory), candidate_id=row['candidate_id'],
                binder_length=len(row['binder_sequence']), binder_sha256=hashlib.sha256(row['binder_sequence'].encode()).hexdigest(),
                chain_and_sequence_identity='pass'))
    before = {str(p): cs.file_digest(p) for p in tracked}
    assert before[str(p4 / 'artifacts/distributed_campaign/ledger.json')] == '0f575ca3e563181423524207c34ebc5fbcf1c6c9b6ff8c5bf71409ae10bc1460'
    assert before[str(p5 / 'artifacts/distributed_campaign/ledger.json')] == '205efddd4c18cea8f8ca9d9ee22a9b5ad193964f5e059854631bf2465e891b39'
    assert before[str(target)] == '7da1c888011cdbd0937e3ff6b73429aa069e4c496cfa0cb8d8b53cc06ca73b57'
    evidence.mkdir(parents=True, exist_ok=True)
    (evidence / 'source-hashes-before.json').write_text(json.dumps(before, indent=2))
    (evidence / 'lineage.json').write_text(json.dumps(dict(target=target_context['identity'], candidates=lineage), indent=2))
    marker = read_json(p5 / 'molecular_design_studio.json')
    records = []
    os.environ['MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR'] = str(evidence / 'runs')
    from verify_multi_generator_pipeline import snapshot
    native_import = snapshot([native], 'bindcraft2', marker['project_id'], marker['studio_owner_id'])
    (evidence / 'native-snapshot.json').write_text(json.dumps({'run_dir': str(native_import)}, indent=2))
    for mode in ('native_complete', 'native_comparative', 'modular_design'):
        graph = cp.new_pipeline(['rfd4'] if mode == 'modular_design' else ['bindcraft2'], mode=mode)
        graph['name'] = 'Phase 7B retained PD-L1 · ' + mode
        for node in graph['nodes']:
            if node['id'] == 'target':
                node['data']['params'].update(target_pdb=str(target), target_chains='A',
                    hotspots='A40,A99,A107' if mode == 'modular_design' else '')
            if node['id'] == 'generator-0':
                node['data']['retained_run_dir'] = str(p5 if mode == 'modular_design' else native_import)
            if node['id'] == 'predictor-0' and mode == 'modular_design':
                node['data']['retained_run_dir'] = str(p5)
            if node['id'].startswith('predictor-') and mode == 'native_comparative':
                node['data']['enabled'] = True  # Empty native collection must submit zero predictor jobs.
            if node['id'] == 'ranking':
                node['data']['params']['ranking_rules'] = [{'field': 'prediction/predictor-0/iptm', 'direction': 'higher'}]
        run = run_graph_pipeline(graph, graph['name'], marker['studio_owner_id'], marker['project_id'],
                                 reuse_completed=False, gpu_device='none')
        result = read_json(run / 'result.json')
        assert result['success'] is True, result
        assert not list((run / 'candidate_batches').rglob('request.json'))
        rows = list(cs.iter_candidates(dict(run_dir=str(run), candidates_jsonl=str(run / 'artifacts/normalized_candidates/candidates.jsonl'))))
        assert len(rows) == (2 if mode == 'modular_design' else 0)
        if mode == 'modular_design':
            assert sorted(len(row['binder_sequence']) for row in rows) == [94, 116]
            assert all(len(cs.metadata(row)['predictions']) == 1 for row in rows)
        archive = run / 'phase7b-retained-results.zip'
        archive.write_bytes(build_run_zip(run))
        with ZipFile(archive) as z:
            assert 'artifacts/candidate_set/workflow_summary.json' in z.namelist()
            assert 'artifacts/candidate_set/native_sources.json' in z.namelist()
        records.append(dict(evidence_class='Retained scientific artifact', mode=mode, run_dir=str(run),
                            candidate_count=len(rows), workflow_summary=result['outputs']['workflow_summary'],
                            archive=str(archive), archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest()))
    after = {p: cs.file_digest(Path(p)) for p in before}
    assert before == after, 'Original retained artifacts changed'
    (evidence / 'source-hashes-after.json').write_text(json.dumps(after, indent=2))
    (evidence / 'summary.json').write_text(json.dumps(dict(source_hashes_unchanged=True, runs=records,
        new_generator_invocations=0, new_predictor_inputs=0), indent=2))
    print(json.dumps([dict(mode=r['mode'], candidates=r['candidate_count'], run_dir=r['run_dir']) for r in records], indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence', type=Path, required=True)
    main(parser.parse_args().evidence)
