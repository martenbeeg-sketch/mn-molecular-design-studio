from __future__ import annotations

from contextlib import contextmanager
from contextvars import ContextVar

_submission = ContextVar("studio_cluster_submission", default=None)


@contextmanager
def submission(command):
    token = _submission.set(command)
    try:
        yield
    finally:
        _submission.reset(token)


def tag_created_job(run_dir, metadata):
    """Set ownership before the ordinary local scheduler can launch the job."""
    command = _submission.get()
    if command is None:
        return
    import json
    from . import config
    metadata.update(cluster_submission_id=command["job_id"], cluster_job_id=command["job_id"],
                    execution_node_id=config.identity(),
                    project_id=command["project_id"], studio_owner_id=command["owner_id"])
    marker = {"application": "mn-molecular-design-studio", "project_id": command["project_id"],
              "studio_owner_id": command["owner_id"], "execution_node_id": config.identity(),
              "cluster_submission_id": command["job_id"], "cluster_job_id": command["job_id"]}
    config.atomic_json(run_dir / "molecular_design_studio.json", marker)
