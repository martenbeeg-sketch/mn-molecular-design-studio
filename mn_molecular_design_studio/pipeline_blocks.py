"""Readable blocks backed by the existing Studio graph and port contracts."""
from functools import lru_cache
import uuid

from mn_molecular_design_studio.pipeline_graph import (
    NODE_CATALOG, connection_compatibility, new_graph, validate_graph,
)

SCREEN = "unidock_active_learning"
SCREEN_INPUTS = ("active_library", "active_seed", "active_target")
SCREEN_ROUND = ("active_docking", "active_accumulate")
SCREEN_CONTINUE = ("active_binder_rule", "active_train", "active_predict", "active_select")


def blank_program():
    return {**new_graph("blank"), "name": "Untitled pipeline", "nodes": [], "edges": []}


def block_action(kind):
    spec = NODE_CATALOG[kind]
    category = spec["category"]
    if kind == "active_decision" or category == "Campaign control":
        return "Decide", "control"
    if kind in {"filter", "shared_filter", "candidate_rank", "active_binder_rule", "active_select"}:
        return "Select", "control"
    if not spec["outputs"]:
        return "Collect", "result"
    if (not spec["inputs"] and not spec["workflow"]) or kind in SCREEN_INPUTS:
        return "Use", "source"
    if category == "Analysis" or kind == "md_analysis":
        return "Measure", "analysis"
    if category == "Collection" or kind == "active_accumulate":
        return "Combine", "analysis"
    if category.startswith("Design") or category in {"Antibody and peptide design", "Ligand generation"}:
        return "Generate", "tool"
    if category == "Sequence design":
        return "Design", "tool"
    if category == "Structure prediction" or kind == "active_predict":
        return "Predict", "tool"
    if kind.startswith("tile_docking_") or kind == "active_docking":
        return "Dock", "tool"
    if kind == "active_train":
        return "Train", "tool"
    return "Run", "tool"


def program(graph):
    """Describe each block once, with named branches and actual predecessors."""
    plan = validate_graph(graph, allow_incomplete=True)
    nodes = {node["id"]: node for node in plan["graph"]["nodes"]}
    depth = {}
    stages = []
    blocks = {}
    for node_id in plan["order"]:
        data = nodes[node_id]["data"]
        action, role = block_action(data["kind"])
        parents = plan["incoming_edges"][node_id]
        depth[node_id] = max((depth[edge["source"]] + 1 for edge in parents), default=0)
        if depth[node_id] == len(stages):
            stages.append([])
        stages[depth[node_id]].append(node_id)
        sources = []
        for edge in parents:
            previous = nodes[edge["source"]]["data"]
            output = next(item for item in NODE_CATALOG[previous["kind"]]["outputs"]
                          if item["id"] == edge.get("sourceHandle", "out"))
            source = previous["label"]
            if len(NODE_CATALOG[previous["kind"]]["outputs"]) > 1:
                source += " · " + output["label"]
            sources.append(source)
        blocks[node_id] = {
            "label": data["label"], "action": action, "role": role,
            "sources": sources, "enabled": data.get("enabled", True),
        }
    result = {"blocks": blocks, "stages": stages}
    if graph.get("template_id") == SCREEN:
        by_kind = {node["data"]["kind"]: node for node in nodes.values()}
        def ids(kinds):
            return [by_kind[kind]["id"] for kind in kinds]
        docking = by_kind["active_docking"]["data"]["params"]
        decision = by_kind["active_decision"]["data"]["params"]
        accumulation = by_kind["active_accumulate"]["data"]["params"]
        tolerance = f"{decision['score_tolerance']:g} {docking['score_units']}".strip()
        conditions = [
            f"Top {decision['top_n']:,} ranks stay in order and scores change by at most "
            f"{tolerance} for {decision['stable_rounds']:,} comparisons",
            "The screening library is exhausted",
        ]
        if limit := decision["max_iterations"]:
            conditions.append(f"The campaign reaches {limit:,} rounds")
        result["screen"] = dict(
            setup=ids(SCREEN_INPUTS), rounds=ids(SCREEN_ROUND), decision=by_kind["active_decision"]["id"],
            finish=ids(("active_results",)), next_round=ids(SCREEN_CONTINUE),
            stop_conditions=conditions, retention=accumulation["retained_pose_count"],
        )
    return result


def add_block(graph, kind, *, source=None, output="out", target_port="in"):
    if graph.get("template_id") == SCREEN:
        raise ValueError("The screening controller manages these blocks.")
    spec = NODE_CATALOG[kind]
    node_id = uuid.uuid4().hex
    position = {"x": 40, "y": 50 + 150 * len(graph["nodes"])}
    if source:
        parent = next(node for node in graph["nodes"] if node["id"] == source)
        branches = sum(edge["source"] == source for edge in graph["edges"])
        position = {"x": parent["position"]["x"] + 310, "y": parent["position"]["y"] + 155 * branches}
    node = {"id": node_id, "type": "studio", "position": position,
            "data": {"kind": kind, "label": spec["label"], "params": {}}}
    updated = {**graph, "nodes": [*graph["nodes"], node]}
    if source:
        return connect_blocks(updated, source, output, node_id, target_port), node_id
    return validate_graph(updated, allow_incomplete=True)["graph"], node_id


def connect_blocks(graph, source, output, target, target_port):
    if graph.get("template_id") == SCREEN:
        raise ValueError("The screening controller manages these connections.")
    connection = {"source": source, "sourceHandle": output, "target": target, "targetHandle": target_port}
    if any(edge["source"] == source and edge["target"] == target
           and edge.get("sourceHandle", "out") == output and edge.get("targetHandle", "in") == target_port
           for edge in graph["edges"]):
        raise ValueError("These ports are already connected.")
    updated = {**graph, "edges": [*graph["edges"], {"id": uuid.uuid4().hex, **connection}]}
    return validate_graph(updated, allow_incomplete=True)["graph"]


@lru_cache(maxsize=None)
def _compatible_ports(kind, output):
    return tuple((target_kind, port["id"]) for target_kind, spec in NODE_CATALOG.items()
                 if spec["category"] != "Active learning" and target_kind != "active_results"
                 for port in spec["inputs"] if connection_compatibility(
                     kind, target_kind, source_handle=output, target_handle=port["id"],
                 )["compatible"])


def next_blocks(graph, source, output):
    """Offer only blocks whose connection the shared editor validator accepts."""
    kind = next(node["data"]["kind"] for node in graph["nodes"] if node["id"] == source)
    compatible = _compatible_ports(kind, output)
    new_blocks = []
    for target_kind, port in compatible:
        try:
            add_block(graph, target_kind, source=source, output=output, target_port=port)
        except ValueError:
            continue
        new_blocks.append((target_kind, port))
    existing = []
    for node in graph["nodes"]:
        target_id = node["id"]
        if target_id == source:
            continue
        for target_kind, port in compatible:
            if target_kind == node["data"]["kind"]:
                try:
                    connect_blocks(graph, source, output, target_id, port)
                except ValueError:
                    continue
                existing.append((target_id, port))
    return new_blocks, existing
