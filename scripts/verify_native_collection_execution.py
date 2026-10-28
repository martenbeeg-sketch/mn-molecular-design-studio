"""Verify native adapter launch grouping without GPU inference or Docker runs.

Pass --image-source with copied /opt/mn_batch_api sources (boltz/ and protenix/
subdirectories). Tests replace only process execution and library model calls.
Retained coordinates are copied to temporary fixtures, never modified.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import uuid

from mn_molecular_design_studio import candidate_pipeline
from mn_molecular_design_studio.native_execution import tool_inventory


def summarize(directory):
    return [dict(case=p.stem, **json.loads(p.read_text())) for p in sorted(Path(directory).glob('*.json'))]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--image-source', type=Path)
    parser.add_argument('--evidence-dir', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--retained', action='store_true')
    parser.add_argument('--retained-evidence-dir', type=Path)
    parser.add_argument('--source-snapshots', type=Path)
    parser.add_argument('--summarize-only', action='store_true')
    args = parser.parse_args()
    observations_dir = args.evidence_dir if args.summarize_only else args.evidence_dir / uuid.uuid4().hex
    root = Path(__file__).resolve().parents[1]
    retained = json.loads((root / 'docs/testing/evidence/phase7a-retained-graph.json').read_text())
    checks = {name: dict(expected=checksum, actual=hashlib.sha256(Path(name).read_bytes()).hexdigest())
              for name, checksum in retained['source_files'].items()}
    assert all(v['expected'] == v['actual'] for v in checks.values()), 'Retained source changed'
    if not args.summarize_only:
        env = os.environ.copy()
        env['MDS_NATIVE_EXECUTION_EVIDENCE_DIR'] = str(observations_dir.resolve())
        if args.image_source: env['MDS_BATCH_API_SOURCE'] = str(args.image_source.resolve())
        if args.retained: env['MDS_RETAINED_PDL1'] = '1'
        else: env.pop('MDS_RETAINED_PDL1', None)
        subprocess.run([sys.executable, '-m', 'pytest', 'tests/test_native_collection_execution.py', '-q'], cwd=root, env=env, check=True)
    assert all(hashlib.sha256(Path(name).read_bytes()).hexdigest() == v['expected'] for name, v in checks.items()), 'Retained source changed during verification'
    graph = candidate_pipeline.new_pipeline()
    forecasts = {str(n): candidate_pipeline.execution_plan(graph, {node['id']: n for node in graph['nodes']
                  if node['id'] in {'sequence-preparation', 'optional-redesign', 'predictor-0', 'predictor-1'}})
                 for n in (1, 5, 10, 25, 50, 100)}
    observations = summarize(observations_dir)
    expected = {f'{engine}-{n}' for engine in ('shared_ligandmpnn', 'foundry') for n in (1,5,10,25,50,100)} | {f'{engine}-seeds{seeds}-{n}' for engine, seeds in (('boltz2',1),('boltz2',3),('protenix_v2',1)) for n in (1,5,10,25,50,100)}
    missing = sorted(expected - {row['case'] for row in observations})
    value = dict(schema='studio.phase7a-native-execution-verification.v1',
        scientific_execution='deterministic process and library fixtures; no new inference',
        observations=observations, native_launch_matrix_complete=not missing, missing_native_launch_fixtures=missing,
        tool_inventory=tool_inventory(), read_only_forecasts=forecasts,
        retained_process_observations=summarize(args.retained_evidence_dir) if args.retained_evidence_dir else [],
        inspected_source_sha256={str(p.relative_to(args.source_snapshots)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in args.source_snapshots.rglob('*') if p.is_file() and (p.suffix == '.py' or p.name == 'mpnn')} if args.source_snapshots else {},
        retained_source_checks=checks,
        limitations=['GPU allocation metadata is simulated in launch fixtures.',
                     'Source-verified model reuse is not measured throughput or numerical equivalence.',
                     'Unverified selectable native lifecycles remain unknown, never counted as one model batch.',
                     'RFD4 inverse-folding exists in ordinary multi-file workflows; dedicated CandidateSet sequence-stage handoff remains unverified.',
                     'Phase 6B campaigns still require one verified native attempt per child.'])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(value, indent=2) + '\n')
    print(json.dumps(dict(observations=len(value['observations']), registered_inventory_entries=len(value['tool_inventory']),
                         retained_sources_unchanged=True, output=str(args.output))))


if __name__ == '__main__':
    main()
