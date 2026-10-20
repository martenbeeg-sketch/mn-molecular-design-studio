"""Typed molecular properties and a small, non-executable formula language."""
from __future__ import annotations
import ast
import math
import operator
import re

MAX_FORMULA_STRING = 100_000
MAX_FORMULA_NUMBER = 1e308
MAX_FORMULA_BITS = 1024


def _bounded_value(value):
    if value is not None and not isinstance(value, (str, int, float, bool)):
        raise ValueError("Formula values must be scalar properties.")
    if isinstance(value, str) and len(value) > MAX_FORMULA_STRING:
        raise ValueError("Formula text exceeds the result limit.")
    if isinstance(value, int) and not isinstance(value, bool) and (value.bit_length() > MAX_FORMULA_BITS or abs(value) > MAX_FORMULA_NUMBER):
        raise ValueError("Formula number exceeds the result limit.")
    if isinstance(value, float) and (not math.isfinite(value) or abs(value) > MAX_FORMULA_NUMBER):
        raise ValueError("Formula number exceeds the result limit.")
    return value


def _bounded_binary(op, left, right):
    if isinstance(op, ast.Add) and isinstance(left, str) and isinstance(right, str):
        if len(left) + len(right) > MAX_FORMULA_STRING:
            raise ValueError("Formula text exceeds the result limit.")
    elif isinstance(op, ast.Mult) and isinstance(left, str) and isinstance(right, int):
        if len(left) * max(0, right) > MAX_FORMULA_STRING:
            raise ValueError("Formula text exceeds the result limit.")
    elif isinstance(op, ast.Mult) and isinstance(right, str) and isinstance(left, int):
        if len(right) * max(0, left) > MAX_FORMULA_STRING:
            raise ValueError("Formula text exceeds the result limit.")
    elif isinstance(op, ast.Pow):
        if not isinstance(left, (int, float)) or not isinstance(right, (int, float)):
            raise ValueError("Formula powers require numbers.")
        if abs(right) > 100 or abs(left) > 1e100:
            raise ValueError("Exponent exceeds formula limits.")
        if isinstance(left, int) and isinstance(right, int) and right > 0 and abs(left) > 1:
            # Estimate result bits before exponentiation to prevent a nested
            # integer power from allocating an enormous Python integer.
            if left.bit_length() * right > MAX_FORMULA_BITS:
                raise ValueError("Formula number exceeds the result limit.")
    elif isinstance(op, ast.Mult) and all(isinstance(value, int) and not isinstance(value, bool) for value in (left, right)):
        if left and right and left.bit_length() + right.bit_length() > MAX_FORMULA_BITS + 1:
            raise ValueError("Formula number exceeds the result limit.")
    return _bounded_value(({ast.Add: operator.add, ast.Sub: operator.sub, ast.Mult: operator.mul,
                           ast.Div: operator.truediv, ast.Mod: operator.mod, ast.Pow: operator.pow}[type(op)])(left, right))


def formula_tree(formula):
    text = str(formula).strip().removeprefix("=").strip()
    if len(text) > 1000:
        raise ValueError("A formula may contain at most 1,000 characters.")
    text = re.sub(r"\[([^\[\]]+)\]", lambda match: "FIELD(" + repr(match.group(1)) + ")", text)
    text = re.sub(r"(?<![<>=!])=(?!=)", "==", text.replace("<>", "!=").replace("^", "**"))
    try:
        tree = ast.parse(text, mode="eval")
    except SyntaxError as exc:
        raise ValueError("The formula syntax is invalid.") from exc
    nodes = list(ast.walk(tree))
    allowed = (ast.Expression, ast.Constant, ast.Name, ast.Attribute, ast.Call, ast.BinOp, ast.UnaryOp, ast.Compare,
               ast.Add, ast.Sub, ast.Mult, ast.Div, ast.Pow, ast.Mod, ast.USub, ast.UAdd, ast.Eq, ast.NotEq,
               ast.Gt, ast.GtE, ast.Lt, ast.LtE, ast.Load)
    if len(nodes) > 256 or any(not isinstance(node, allowed) for node in nodes):
        raise ValueError("This formula uses unsupported expressions.")
    for node in nodes:
        if isinstance(node, ast.Call) and (not isinstance(node.func, ast.Name) or node.keywords):
            raise ValueError("Call a supported formula function.")
        if isinstance(node, (ast.Name, ast.Attribute)) and any(part.startswith("_") for part in _name(node).split(".")):
            raise ValueError("Private names are not supported in formulas.")
    return tree


def _name(node):
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        return _name(node.value) + "." + node.attr
    raise ValueError("Invalid column reference.")


def evaluate(formula, values):
    lookup = {str(key).casefold(): value for key, value in values.items()}
    comparisons = {ast.Eq: operator.eq, ast.NotEq: operator.ne, ast.Gt: operator.gt, ast.GtE: operator.ge, ast.Lt: operator.lt, ast.LtE: operator.le}
    functions = {"ABS": abs, "ROUND": round, "MIN": lambda *args: min(args), "MAX": lambda *args: max(args),
                 "SUM": lambda *args: sum(args), "LOG10": math.log10, "LOG": math.log, "SQRT": math.sqrt,
                 "AND": lambda *args: all(args), "OR": lambda *args: any(args), "NOT": operator.not_,
                 "LEN": len, "UPPER": lambda value: str(value).upper(), "LOWER": lambda value: str(value).lower(),
                 "ISBLANK": lambda value: value in (None, ""), "CONCAT": lambda *args: "".join(map(str, args))}
    def visit(node):
        if isinstance(node, ast.Expression):
            return visit(node.body)
        if isinstance(node, ast.Constant):
            if isinstance(node.value, (str, int, float, bool)) or node.value is None:
                return _bounded_value(node.value)
        elif isinstance(node, (ast.Name, ast.Attribute)):
            name = _name(node).casefold()
            if name in {"true", "false"}:
                return name == "true"
            return _bounded_value(lookup.get(name))
        elif isinstance(node, ast.UnaryOp):
            return _bounded_value(-visit(node.operand) if isinstance(node.op, ast.USub) else +visit(node.operand))
        elif isinstance(node, ast.BinOp):
            left, right = visit(node.left), visit(node.right)
            return _bounded_binary(node.op, left, right)
        elif isinstance(node, ast.Compare):
            left = visit(node.left)
            for op, right_node in zip(node.ops, node.comparators):
                right = visit(right_node)
                if not comparisons[type(op)](left, right):
                    return False
                left = right
            return True
        elif isinstance(node, ast.Call):
            name = node.func.id.upper()
            if name == "IF" and len(node.args) == 3:
                return visit(node.args[1]) if visit(node.args[0]) else visit(node.args[2])
            if name == "IFERROR" and len(node.args) == 2:
                try:
                    return visit(node.args[0])
                except (TypeError, ValueError, ZeroDivisionError, OverflowError):
                    return visit(node.args[1])
            if name == "FIELD" and len(node.args) == 1:
                return _bounded_value(lookup.get(str(visit(node.args[0])).casefold()))
            if name not in functions or len(node.args) > 64:
                raise ValueError("Unsupported formula function or argument count.")
            arguments = [visit(arg) for arg in node.args]
            if name == "CONCAT":
                lengths = [len(value) if isinstance(value, str) else len(str(value)) for value in arguments]
                if sum(lengths) > MAX_FORMULA_STRING:
                    raise ValueError("Formula text exceeds the result limit.")
            if name == "ROUND" and len(arguments) == 2 and abs(arguments[1]) > 100:
                raise ValueError("ROUND precision exceeds formula limits.")
            return _bounded_value(functions[name](*arguments))
        raise ValueError("Unsupported formula value.")
    try:
        value = visit(formula_tree(formula))
        if isinstance(value, float) and not math.isfinite(value):
            return "#NUM!"
        if isinstance(value, (int, float)) and abs(value) > 1e308:
            return "#NUM!"
        return value
    except ZeroDivisionError:
        return "#DIV/0!"
    except (TypeError, ValueError, OverflowError, RecursionError):
        return "#VALUE!"


def validate_schema(schema):
    fields = schema.get("fields", [])
    if not isinstance(fields, list) or len(fields) > 100:
        raise ValueError("A group supports at most 100 property fields.")
    seen = set()
    for field in fields:
        name = str(field.get("name") or "").strip()
        if not name or name in seen or len(name) > 100:
            raise ValueError("Property names must be non-empty and unique.")
        seen.add(name)
        if field.get("type", "text") not in {"text", "number", "boolean", "category"}:
            raise ValueError("Choose a supported property type.")
        if field.get("formula"):
            formula_tree(field["formula"])
            if field.get("required"):
                raise ValueError("A computed column cannot be required.")
    return schema


def validate_record(record, schema):
    for field in schema.get("fields", []):
        value = (record.get("metadata") or {}).get(field["name"])
        if field.get("formula"):
            if value is not None:
                raise ValueError(f"{field['name']} is computed; omit it from supplied properties.")
            continue
        if value in (None, ""):
            if field.get("required"):
                raise ValueError(f"Required property is missing: {field['name']}.")
            continue
        kind = field.get("type", "text")
        if kind == "number" and (type(value) not in {int, float} or not math.isfinite(value)):
            raise ValueError(f"{field['name']} must be a finite number.")
        if kind == "boolean" and type(value) is not bool:
            raise ValueError(f"{field['name']} must be true or false.")
        if kind == "text" and not isinstance(value, str):
            raise ValueError(f"{field['name']} must be text.")
        if kind == "category" and field.get("choices") and value not in field["choices"]:
            raise ValueError(f"{field['name']} is not one of the allowed categories.")
    for chain in schema.get("chains", []):
        if chain.get("required") and not any(row["chain"] == chain["chain"] and row["type"] == chain["type"] for row in record["components"]):
            raise ValueError(f"Required {chain['type']} chain {chain['chain']} is missing.")
