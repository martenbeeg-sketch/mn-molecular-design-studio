"""Shared scientific settings, independent of subject identity and placement."""
from __future__ import annotations


PROTOCOL_FIELDS = {"method_version", "engine_version", "protocol_revision", "image", "image_digest", "checkpoint_sha256"}
PLACEMENT_FIELDS = {"source_run_dir", "candidates_jsonl", "candidate_ids", "selected_structure_path", "selected_model_id",
                    "model_id", "structure_digest", "binder_chains", "target_chains", "source_job", "job_name",
                    "campaign_name", "backend", "cpu_cores", "cpu_threads", "nprocs", "workers", "gpu_device",
                    "output_dir", "method", "assessment_protocol", "database_source"}


def _canonical_value(value):
    if isinstance(value, dict):
        return {key: _canonical_value(child) for key, child in value.items()}
    if isinstance(value, (list, tuple)):
        return [_canonical_value(child) for child in value]
    if type(value) is float and value.is_integer():
        return int(value)
    return value


def scientific_parameters(parameters, method=None):
    """Compare effective scientific settings without weakening exact equality."""
    settings = {key: value for key, value in parameters.items()
                if key not in PLACEMENT_FIELDS | PROTOCOL_FIELDS and value is not None
                and not (key == "interface_ligand_selection" and value == {})}
    if method in {"analysis.binder_interface_chemistry", "binder_interface_chemistry_v1"}:
        from mn_molecular_design_studio.complex_analysis_interface_chemistry import fixed_parameters
        settings = {**fixed_parameters(), **settings}
    return _canonical_value(settings)


def model_identity(*values):
    return str(next((value for value in values if value is not None and str(value) != ""), ""))
