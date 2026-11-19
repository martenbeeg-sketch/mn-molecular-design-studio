"""Readers turn a run folder of one design tool into a model.RunData.

A reader module provides NAME (str), detect(path) -> bool and
read_run(path, run_id, label=None, settings_path=None, **options) -> RunData.
To add a reader, append (NAME, "binder_review.readers.<module>") to _REGISTRY; modules that fail to
import are skipped, so a broken optional reader never blocks the others.
"""

from __future__ import annotations

import importlib
from pathlib import Path
from types import ModuleType
from typing import Optional

class ReaderError(ValueError):
    """A run folder that cannot be read, said in words for the person who chose it.

    Distinct from a bug: the caller prints it and stops, rather than showing a
    traceback. A ValueError so that code written against the older readers,
    which raised a plain one, still catches it.
    """


_REGISTRY: list[tuple[str, str]] = [
    ("bindcraft1", "binder_review.readers.bindcraft1"),
    ("bindcraft2", "binder_review.readers.bindcraft2"),
]


def register_reader(name: str, module: str) -> None:
    """Add a reader at run time (after the built-in ones). Replaces an entry of the same name."""
    for i, (n, _) in enumerate(_REGISTRY):
        if n == name:
            _REGISTRY[i] = (name, module)
            return
    _REGISTRY.append((name, module))


def available_readers() -> dict[str, ModuleType]:
    """Reader name -> imported module, for every registered reader that imports cleanly."""
    out: dict[str, ModuleType] = {}
    for name, mod in _REGISTRY:
        try:
            out[name] = importlib.import_module(mod)
        except ImportError:
            continue
    return out


def detect_reader(path: str | Path) -> Optional[ModuleType]:
    """The first registered reader whose detect(path) is true, or None."""
    for module in available_readers().values():
        try:
            if module.detect(path):
                return module
        except OSError:
            continue
    return None


def get_reader(name: str) -> ModuleType:
    readers = available_readers()
    if name not in readers:
        raise KeyError(f"unknown reader {name!r}; available: {sorted(readers)}")
    return readers[name]
