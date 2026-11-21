"""Readable blocks backed by the existing Studio graph and port contracts."""
from functools import lru_cache
import math
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


def port_label(port):
    label = port["label"]
    if label == "molecule_task":
        return "Prediction task"
    return label.replace("_", " ").capitalize() if "_" in label else label


def canvas_program(graph):
    """Expose editable layout, named sockets and actual wires without tool parameters."""
    nodes = {}
    compatibility = {}
    managed = graph.get("template_id") == SCREEN
    for node in graph["nodes"]:
        kind = node["data"]["kind"]
        spec = NODE_CATALOG[kind]
        nodes[node["id"]] = {
            "kind": kind, "position": node.get("position") or {"x": 0, "y": 0},
            **{direction: [{"id": port["id"], "label": port_label(port), "multiple": bool(port.get("multiple"))}
                           for port in spec[direction]] for direction in ("inputs", "outputs")},
        }
        if not managed and kind not in compatibility:
            compatibility[kind] = {port["id"]: _compatible_ports(kind, port["id"])
                                   for port in spec["outputs"]}
    return {"nodes": nodes, "edges": graph["edges"], "compatible": compatibility, "managed": managed}


def canvas_cards(graph):
    if graph.get("template_id") == SCREEN:
        return {}
    cards = {}
    for kind, spec in NODE_CATALOG.items():
        if spec["category"] == "Active learning" or kind == "active_results":
            continue
        action, role = block_action(kind)
        cards[kind] = {"label": spec["label"], "action": action, "role": role, "targets": []}
    if not graph["nodes"]:
        cards[SCREEN] = {"label": "Large-library screen", "action": "Repeat", "role": "control", "targets": []}
    return cards


def _position(value):
    if not isinstance(value, dict) or any(type(value.get(axis)) not in (int, float)
            or not -1e7 < value[axis] < 1e7 or not math.isfinite(value[axis]) for axis in ("x", "y")):
        raise ValueError("Choose a valid position on the canvas.")
    return {axis: value[axis] for axis in ("x", "y")}


def move_canvas_block(graph, node_id, position):
    position = _position(position)
    if graph.get("template_id") == SCREEN:
        raise ValueError("The screening controller manages this layout.")
    if not any(node["id"] == node_id for node in graph["nodes"]):
        raise ValueError("This block is no longer on the canvas.")
    return {**graph, "nodes": [dict(node, position=position) if node["id"] == node_id else node
                               for node in graph["nodes"]]}


def place_canvas_block(graph, kind, position, *, source=None, output=None):
    if not isinstance(kind, str) or kind not in canvas_cards(graph):
        raise ValueError("This block is no longer available in the palette.")
    if kind == SCREEN:
        screen = new_graph(SCREEN)
        return screen, screen["nodes"][0]["id"]
    position = _position(position)
    target_port = "in"
    if source is not None:
        parent = next((item for item in graph["nodes"] if item["id"] == source), None)
        if not parent or not isinstance(output, str):
            raise ValueError("This output is no longer on the canvas.")
        choices = [port for target, port in _compatible_ports(parent["data"]["kind"], output) if target == kind]
        if len(choices) != 1:
            raise ValueError("Place this block on the canvas and connect to a specific input.")
        target_port = choices[0]
    updated, node_id = add_block(graph, kind, source=source, output=output, target_port=target_port)
    return move_canvas_block(updated, node_id, position), node_id


def _new_block(graph, kind, source):
    if graph.get("template_id") == SCREEN:
        raise ValueError("The screening controller manages these blocks.")
    spec = NODE_CATALOG[kind]
    node_id = uuid.uuid4().hex
    position = {"x": 40, "y": 50 + 150 * len(graph["nodes"])}
    if source:
        parent = next(node for node in graph["nodes"] if node["id"] == source)
        branches = sum(edge["source"] == source for edge in graph["edges"])
        position = {"x": parent["position"]["x"] + 310, "y": parent["position"]["y"] + 155 * branches}
    return {"id": node_id, "type": "studio", "position": position,
            "data": {"kind": kind, "label": spec["label"], "params": {}}}


def _edge(source, output, target, target_port):
    return {"id": uuid.uuid4().hex, "source": source, "sourceHandle": output,
            "target": target, "targetHandle": target_port}


def add_block(graph, kind, *, source=None, output="out", target_port="in"):
    node = _new_block(graph, kind, source)
    node_id = node["id"]
    updated = {**graph, "nodes": [*graph["nodes"], node]}
    if source:
        return connect_blocks(updated, source, output, node_id, target_port), node_id
    return validate_graph(updated, allow_incomplete=True)["graph"], node_id


def connect_blocks(graph, source, output, target, target_port):
    if graph.get("template_id") == SCREEN:
        raise ValueError("The screening controller manages these connections.")
    if any(edge["source"] == source and edge["target"] == target
           and edge.get("sourceHandle", "out") == output and edge.get("targetHandle", "in") == target_port
           for edge in graph["edges"]):
        raise ValueError("These ports are already connected.")
    updated = {**graph, "edges": [*graph["edges"], _edge(source, output, target, target_port)]}
    return validate_graph(updated, allow_incomplete=True)["graph"]


def insert_block(graph, edge_id, kind, target_port):
    """Replace one wire with both new wires as a single validated edit."""
    edge = next(item for item in graph["edges"] if item["id"] == edge_id)
    node = _new_block(graph, kind, edge["source"])
    before = _edge(edge["source"], edge.get("sourceHandle", "out"), node["id"], target_port)
    after = _edge(node["id"], "out", edge["target"], edge.get("targetHandle", "in"))
    updated = {
        **graph, "nodes": [*graph["nodes"], node],
        "edges": [item for item in graph["edges"] if item["id"] != edge_id] + [before, after],
    }
    return validate_graph(updated, allow_incomplete=True)["graph"], node["id"]


@lru_cache(maxsize=None)
def _compatible_ports(kind, output):
    return tuple((target_kind, port["id"]) for target_kind, spec in NODE_CATALOG.items()
                 if spec["category"] != "Active learning" and target_kind != "active_results"
                 for port in spec["inputs"] if connection_compatibility(
                     kind, target_kind, source_handle=output, target_handle=port["id"],
                 )["compatible"])


def next_blocks(graph, source, output, *, include_new=True, include_existing=True):
    """Offer only blocks whose connection the shared editor validator accepts."""
    kind = next(node["data"]["kind"] for node in graph["nodes"] if node["id"] == source)
    compatible = _compatible_ports(kind, output)
    new_blocks = []
    for target_kind, port in compatible if include_new else ():
        try:
            add_block(graph, target_kind, source=source, output=output, target_port=port)
        except ValueError:
            continue
        new_blocks.append((target_kind, port))
    existing = []
    for node in graph["nodes"] if include_existing else ():
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


def placements(graph, selected):
    """Build the palette and exact allowed sockets for the current selection."""
    slots, cards = {}, {}
    if graph.get("template_id") == SCREEN:
        return slots, cards

    def offer(socket, kind, port=""):
        if kind == SCREEN:
            label, action, role = "Large-library screen", "Repeat", "control"
        else:
            spec = NODE_CATALOG[kind]
            label, (action, role) = spec["label"], block_action(kind)
            if port and len(spec["inputs"]) > 1:
                label += " · " + next(item["label"] for item in spec["inputs"] if item["id"] == port)
        card = cards.setdefault(f"{kind}/{port}", {
            "kind": kind, "port": port, "label": label, "action": action, "role": role, "targets": [],
        })
        card["targets"].append(socket)

    slots["start"] = {"label": "Start a stack", "action": "start"}
    for kind, spec in NODE_CATALOG.items():
        if (spec["outputs"] and spec["category"] != "Active learning"
                and all(port.get("optional") for port in spec["inputs"])):
            offer("start", kind)
    if not graph["nodes"]:
        offer("start", SCREEN)

    nodes = {node["id"]: node["data"] for node in graph["nodes"]}
    if selected in nodes:
        spec = NODE_CATALOG[nodes[selected]["kind"]]
        for output in spec["outputs"]:
            handle = output["id"]
            choices, _ = next_blocks(graph, selected, handle, include_existing=False)
            outgoing = [edge for edge in graph["edges"]
                        if edge["source"] == selected and edge.get("sourceHandle", "out") == handle]
            socket = f"after/{handle}"
            label = "Add a branch" if outgoing else "Snap in the next block"
            if len(spec["outputs"]) > 1:
                label = f'{output["label"]} · {label}'
            if choices:
                slots[socket] = {"label": label, "action": "after", "source": selected, "output": handle}
                for kind, port in choices:
                    offer(socket, kind, port)
            for edge in outgoing:
                socket = f"between/{edge['id']}"
                destination = nodes[edge["target"]]
                for kind, port in choices:
                    if (destination["kind"], edge.get("targetHandle", "in")) not in _compatible_ports(kind, "out"):
                        continue
                    try:
                        insert_block(graph, edge["id"], kind, port)
                    except ValueError:
                        continue
                    slots[socket] = {"label": f'Insert before {destination["label"]}',
                                     "action": "insert", "source": selected, "edge": edge["id"]}
                    offer(socket, kind, port)
    return slots, cards


def drop_block(graph, slots, cards, socket, card_id):
    """Apply a placement offered by the editor to this graph revision."""
    card = cards.get(card_id) if isinstance(card_id, str) else None
    if not card or socket not in card["targets"]:
        raise ValueError("That block does not fit this socket.")
    slot = slots[socket]
    kind, port = card["kind"], card["port"]
    if slot["action"] == "start":
        if kind == SCREEN:
            screen = new_graph(SCREEN)
            return screen, screen["nodes"][0]["id"]
        return add_block(graph, kind)
    if slot["action"] == "insert":
        return insert_block(graph, slot["edge"], kind, port)
    return add_block(graph, kind, source=slot["source"], output=slot["output"], target_port=port)
