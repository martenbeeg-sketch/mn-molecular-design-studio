"""Best-effort host wall timings; these do not measure GPU kernel/model time."""
from __future__ import annotations

from contextlib import contextmanager
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import time
from typing import Iterator
from uuid import uuid4


def _emit(run_dir: Path, record: dict) -> None:
    # A single append per small record avoids rewriting another process's events.
    # Timing diagnostics must never change a scientific job's outcome.
    try:
        payload = (json.dumps({"schema_version": 1, "at": datetime.now(timezone.utc).isoformat(),
                               "pid": os.getpid(), **record}) + "\n").encode()
        descriptor = os.open(Path(run_dir) / "execution_lifecycle.jsonl",
                             os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o600)
        try:
            os.write(descriptor, payload)
        finally:
            os.close(descriptor)
    except OSError:
        pass


@contextmanager
def execution_phase(run_dir: Path, phase: str) -> Iterator[dict]:
    """Pair events by span_id, including exceptions and nonzero process exits.

    An unpaired start means the host was interrupted or telemetry is incomplete;
    it must not be interpreted as a completed interval. No command, seed, input
    sequence, or exception message is copied into telemetry.
    """
    span_id = str(uuid4())
    start = time.perf_counter()
    record = {"span_id": span_id, "phase": phase}
    _emit(run_dir, {**record, "event": "started"})
    outcome = {}
    try:
        yield outcome
    except BaseException as exc:
        outcome = {"exception_type": type(exc).__name__}
        raise
    finally:
        _emit(run_dir, {**record, **outcome, "event": "finished",
                        "elapsed_seconds": max(0.0, time.perf_counter() - start)})
