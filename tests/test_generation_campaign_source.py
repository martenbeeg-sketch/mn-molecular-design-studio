import json

from mn_molecular_design_studio.generation_io import generated_candidates
from mn_molecular_design_studio.engine.workflow_engine.workflows.graph_pipeline import source_generations


def test_rfd4_campaign_outputs_are_reusable_generation_sources(tmp_path):
    for name, value in {
        'metadata.json': {'status': 'completed', 'tool': 'rfd4_proteina'},
        'result.json': {'success': True},
        'input.json': {'params': {'distributed_campaign': {'enabled': True}}},
        'molecular_design_studio.json': {'project_id': 'project'},
    }.items():
        (tmp_path / name).write_text(json.dumps(value))
    path = tmp_path / 'artifacts/normalized_candidates/candidates.jsonl'
    path.parent.mkdir(parents=True)
    row = {'candidate_id': 'design-1', 'stage': 'generation.backbone_sequence',
           'source_tool': 'rfd4_proteina', 'binder_sequence': 'AAA',
           'metrics': {'result_kind': 'sequence_structure_codesign'}}
    path.write_text(json.dumps(row) + '\n')
    assert [r['candidate_id'] for r in generated_candidates(tmp_path)] == ['design-1']
    assert source_generations(str(tmp_path), 'project')[1][0]['candidate_id'] == 'design-1'
    assert generated_candidates(tmp_path, stage='generation.backbone') == []
    (tmp_path / 'input.json').write_text(json.dumps({'params': {}}))
    assert generated_candidates(tmp_path) == []
    (tmp_path / 'input.json').write_text(json.dumps({'params': {'distributed_campaign': {'enabled': True}}}))
    (tmp_path / 'result.json').write_text(json.dumps({'success': False}))
    assert generated_candidates(tmp_path) == []
