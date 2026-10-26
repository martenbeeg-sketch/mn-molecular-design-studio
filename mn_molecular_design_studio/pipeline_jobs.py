"""Durable parent identity for pipeline jobs scheduled on different computers."""
from pathlib import Path


def link_pipeline_child(run_dir, *, parent_run_dir, parent_node_id, pipeline_name,
                        step, step_index="", round_number=None):
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json, write_json
    parent_id = Path(str(parent_run_dir).rstrip("/")).name
    if not parent_id or not str(parent_node_id).strip():
        raise ValueError("A pipeline child requires its parent run and workstation identity.")
    path = Path(run_dir) / "metadata.json"
    metadata = read_json(path)
    metadata.update(pipeline_parent_run_dir=str(parent_run_dir), pipeline_parent_node_id=str(parent_node_id),
                    campaign_id=parent_id, pipeline_name=str(pipeline_name), campaign_step=str(step),
                    campaign_step_index=str(step_index))
    if round_number is not None:
        metadata.update(active_learning_parent_run=str(parent_run_dir), active_learning_round=int(round_number),
                        pipeline_iteration=int(round_number))
    write_json(path, metadata)
