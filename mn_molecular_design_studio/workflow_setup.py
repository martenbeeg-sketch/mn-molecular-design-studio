"""Client-independent preparation and discovery for scientific workflow forms."""
from __future__ import annotations

import copy
from pathlib import Path


def setup_contract(workflow_id):
    from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import workflow_schema
    from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG
    schema = workflow_schema(workflow_id)
    nodes = [copy.deepcopy(node) for node in NODE_CATALOG.values() if node.get("workflow") == workflow_id]
    result = {"workflow": workflow_id, "schema": schema, "pipeline_nodes": nodes,
              "helpers": [], "examples": []}
    if workflow_id == "design.kejia_idr":
        from mn_molecular_design_studio import kejia_idr_input as inputs
        result.update(template_libraries=[{"id": key, "count": count} for key, count in inputs.library_counts().items()],
                      helpers=["preview", "templates", "template-file", "stage-template", "example", "pipeline"],
                      examples=[{"id": "tamarind", "target_sequence": inputs.EXAMPLE_SEQUENCE,
                                 "template_library": "upload", "binder_chain": "A", "peptide_chain": "B",
                                 "template_available": inputs.example_template().is_file()}])
    elif workflow_id == "design.peptide_trap":
        from mn_molecular_design_studio import peptide_trap_input as inputs
        result.update(templates=[{"id": key, **value, "available": inputs.template_path(key).is_file()}
                                 for key, value in inputs.SCAFFOLDS.items()],
                      helpers=["preview", "templates", "template-file", "stage-template", "example", "pipeline"],
                      examples=[{"id": "c37", "target_sequence": inputs.EXAMPLE_SEQUENCE, "scaffold_ids": ["C37"]}])
    if workflow_id in {"design.kejia_idr", "design.peptide_trap"}:
        result["default_resources"] = {"gpu_device": "none", "cpu_cores": 4}
    return result


def templates(workflow_id, library="", *, offset=0, limit=100):
    if workflow_id == "design.kejia_idr":
        from mn_molecular_design_studio import kejia_idr_input as inputs
        if library not in inputs.LIBRARIES:
            raise ValueError("Choose walle, walle2, mini or strand.")
        paths = sorted((inputs.source_root() / "templates" / library).glob("*.pdb"))
        rows = [{"id": p.name, "library": library, "binder_chain": "A", "peptide_chain": "B"} for p in paths]
    elif workflow_id == "design.peptide_trap":
        from mn_molecular_design_studio import peptide_trap_input as inputs
        rows = [{"id": key, **value, "available": inputs.template_path(key).is_file(),
                 "binder_chain": "B", "peptide_chain": "A"} for key, value in inputs.SCAFFOLDS.items()]
    else:
        raise ValueError("This workflow does not provide a template library helper.")
    return {"items": rows[offset:offset + limit], "total": len(rows), "offset": offset, "limit": limit}


def template_file(workflow_id, template_id, library=""):
    if workflow_id == "design.kejia_idr":
        from mn_molecular_design_studio import kejia_idr_input as inputs
        if library not in inputs.LIBRARIES or Path(template_id).name != template_id or "\\" in template_id:
            raise ValueError("Choose an installed library and its template filename.")
        root = inputs.source_root() / "templates" / library
        path = root / template_id
        if path.suffix != ".pdb" or not path.resolve().is_relative_to(root.resolve()):
            raise ValueError("Invalid template filename.")
    elif workflow_id == "design.peptide_trap":
        from mn_molecular_design_studio import peptide_trap_input as inputs
        if template_id not in inputs.SCAFFOLDS:
            raise ValueError("Choose an available peptide-trap scaffold.")
        path = inputs.template_path(template_id)
    else:
        raise ValueError("This workflow does not provide a template library helper.")
    if not path.is_file():
        raise FileNotFoundError("Install the template library on this workstation first.")
    return path


def preview(workflow_id, parameters, *, limit=100, offset=0):
    if workflow_id == "design.kejia_idr":
        from mn_molecular_design_studio import kejia_idr_input as inputs
        prepared = inputs.prepare_request(parameters)
        rows = prepared.pop("templates")
        normalized = {k: v for k, v in prepared.items() if k != "template_structure"}
        return {"parameters": normalized, "template_count": len(rows),
                "templates": [{"id": Path(row["path"]).name, "binder_chain": row["binder_chain"],
                               "peptide_chain": row["peptide_chain"]} for row in rows[offset:offset + limit]],
                "offset": offset, "limit": limit}
    if workflow_id == "design.peptide_trap":
        from mn_molecular_design_studio import peptide_trap_input as inputs
        prepared = inputs.prepare_request(parameters, require_files=False)
        all_matches = prepared.pop("all_matches")
        registers = prepared.pop("registers")
        try:
            inputs.validate_library(prepared["scaffold_ids"])
            installed, error = True, None
        except (ValueError, OSError) as exc:
            installed, error = False, str(exc)
        return {"parameters": prepared, "matches": all_matches[offset:offset + limit],
                "match_count": len(all_matches), "selected_count": len(registers),
                "offset": offset, "limit": limit, "assets_ready": installed, "asset_error": error}
    raise ValueError("This workflow has no preparation preview helper; use workflow validation.")
