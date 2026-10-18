"""Sequence-prompt framework records for the Biohub ESMFold2 binder tile."""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Any

import yaml


FRAMEWORK_FILE = Path(__file__).with_name("esmfold2_vhh_frameworks.yaml")


@lru_cache(maxsize=1)
def _framework_records() -> dict[str, dict[str, Any]]:
    payload = yaml.safe_load(FRAMEWORK_FILE.read_text()) or {}
    records = payload.get("frameworks")
    if not isinstance(records, list) or not records:
        raise ValueError(f"No ESMFold2 VHH frameworks found in {FRAMEWORK_FILE}.")
    by_name: dict[str, dict[str, Any]] = {}
    for record in records:
        if not isinstance(record, dict):
            raise ValueError(f"Invalid VHH framework record in {FRAMEWORK_FILE}.")
        names = [record.get("id"), record.get("canonical_name"), *(record.get("aliases") or [])]
        for name in names:
            if not isinstance(name, str) or not name.strip():
                continue
            key = name.strip().lower()
            previous = by_name.get(key)
            if previous is not None and previous is not record:
                raise ValueError(f"Duplicate VHH framework alias: {name}")
            by_name[key] = record
    return by_name


def vhh_framework_choices() -> list[dict[str, Any]]:
    """Return one record per named bundled framework, in panel order."""
    records = _framework_records()
    choices: list[dict[str, Any]] = []
    seen: set[str] = set()
    for record in records.values():
        framework_id = str(record["id"])
        if framework_id in seen:
            continue
        seen.add(framework_id)
        choices.append(record)
    return choices


def get_vhh_framework(name: str) -> dict[str, Any]:
    """Resolve an id, canonical name, or alias to its packaged VHH record."""
    try:
        return dict(_framework_records()[name.strip().lower()])
    except KeyError as exc:
        choices = ", ".join(str(record["id"]) for record in vhh_framework_choices())
        raise ValueError(f"Unknown VHH framework {name!r}; choose one of: {choices}.") from exc
