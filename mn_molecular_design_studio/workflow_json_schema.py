"""JSON Schema projection of the existing public workflow descriptors.

Only syntax trees are inspected: type descriptor text is never evaluated.
Runtime workflow validation remains authoritative.
"""
from __future__ import annotations

import ast


def type_schema(descriptor):
    def convert(node):
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.BitOr):
            choices = []
            for side in (node.left, node.right):
                value = convert(side)
                choices.extend(value.get("anyOf", [value]))
            return {"anyOf": choices}
        if isinstance(node, ast.Constant):
            return {"type": "null"} if node.value is None else {}
        if isinstance(node, (ast.Name, ast.Attribute)):
            name = ast.unparse(node).split(".")[-1]
            return {"str": {"type": "string"}, "Path": {"type": "string"},
                    "int": {"type": "integer"}, "float": {"type": "number"},
                    "bool": {"type": "boolean"}, "NoneType": {"type": "null"},
                    "dict": {"type": "object"}, "list": {"type": "array"},
                    "tuple": {"type": "array"}}.get(name, {})
        if isinstance(node, ast.Subscript):
            name = ast.unparse(node.value).split(".")[-1]
            arguments = list(node.slice.elts) if isinstance(node.slice, ast.Tuple) else [node.slice]
            if name == "Literal":
                return {"enum": [ast.literal_eval(a) for a in arguments]}
            if name in {"Union", "Optional"}:
                values = [convert(a) for a in arguments]
                if name == "Optional":
                    values.append({"type": "null"})
                return {"anyOf": values}
            if name in {"list", "List", "Sequence", "Iterable", "set", "Set"}:
                return {"type": "array", "items": convert(arguments[0])}
            if name in {"dict", "Dict", "Mapping"}:
                return {"type": "object", "additionalProperties": convert(arguments[-1])}
            if name in {"tuple", "Tuple"}:
                if len(arguments) == 2 and isinstance(arguments[-1], ast.Constant) and arguments[-1].value is Ellipsis:
                    return {"type": "array", "items": convert(arguments[0])}
                return {"type": "array", "prefixItems": [convert(a) for a in arguments],
                        "minItems": len(arguments), "maxItems": len(arguments)}
        return {}
    try:
        return convert(ast.parse(str(descriptor), mode="eval").body)
    except (SyntaxError, TypeError, ValueError):
        return {}


def parameters_schema(parameters):
    properties, required = {}, []
    for parameter in parameters:
        spec = type_schema(parameter.get("type", "Any"))
        if "default" in parameter:
            spec["default"] = parameter["default"]
        properties[parameter["name"]] = spec
        if parameter.get("required"):
            required.append(parameter["name"])
    return {"type": "object", "properties": properties, "required": required,
            "additionalProperties": False}
