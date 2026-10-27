"""Summarize retained real candidate benchmark evidence without running models."""
from __future__ import annotations
import argparse
import json
from pathlib import Path


def summarize(path):
    evidence = json.loads(Path(path).read_text())
    result = {k: evidence[k] for k in ('kind', 'retained_source', 'input_ids', 'effective_parameters', 'gpu_inventory') if k in evidence}
    result['artifact_file'] = str(Path(path).resolve())
    result['trials'] = []
    for trial in evidence.get('trials', []):
        devices = {}
        for sample in trial['telemetry']:
            for line in sample['gpus'].splitlines():
                index, utilization, used = (value.strip() for value in line.split(','))
                record = devices.setdefault(index, {'utilization_samples': [], 'memory_samples_mib': []})
                record['utilization_samples'].append(float(utilization))
                record['memory_samples_mib'].append(float(used))
        telemetry = {index: {'mean_utilization_percent': sum(v['utilization_samples'])/len(v['utilization_samples']),
                            'peak_utilization_percent': max(v['utilization_samples']),
                            'peak_device_vram_mib': max(v['memory_samples_mib']),
                            'minimum_device_vram_mib': min(v['memory_samples_mib']),
                            'sample_count': len(v['utilization_samples'])} for index, v in devices.items()}
        attempts = [a for b in trial['state']['batches'] for a in b['attempts']]
        result['trials'].append(dict(label=trial['label'], wall_seconds=trial['wall_seconds'],
            candidates_per_hour=trial['candidates_per_hour'], candidates=trial['candidate_count'],
            allocations=trial['allocations'], telemetry=telemetry,
            controller_seconds=trial.get('controller_seconds'), scheduler_seconds=trial.get('scheduler_seconds'),
            transfer_bytes=trial.get('transfer_bytes', 0), transfer_seconds=trial.get('transfer_seconds', 0),
            retries=trial.get('retries', sum(max(0,len(b['attempts'])-1) for b in trial['state']['batches'])),
            failure_retry_cost_seconds=0 if len(attempts)==len(trial['state']['batches']) else None,
            model_initialization_seconds=None, native_cpu_preparation_seconds=None,
            jobs=[a['run_dir'] for a in attempts]))
    if evidence.get('cpu_analysis'):
        result['cpu_analysis'] = dict(wall_seconds=evidence['cpu_analysis']['wall_seconds'],
            allocations=[b['attempts'][-1]['allocation'] for b in evidence['cpu_analysis']['state']['batches']],
            output=evidence['cpu_analysis']['output'])
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('evidence', type=Path, nargs='+')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    content = json.dumps([summarize(p) for p in args.evidence], indent=2)
    if args.output:
        args.output.write_text(content + '\n')
    else:
        print(content)
