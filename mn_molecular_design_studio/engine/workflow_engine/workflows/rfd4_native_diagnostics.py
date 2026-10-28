"""Opt-in, mounted RFD4 timing launcher; never edits installed engine files.

All durations are host wall intervals. No CUDA synchronization, RNG operation,
model argument, result, worker count or native sampling setting is changed.
"""
from __future__ import annotations

from contextlib import contextmanager
from datetime import datetime, timezone
from functools import wraps
import hashlib
import json
import os
from pathlib import Path
import sys
import time
from uuid import uuid4
import warnings


# Refuse to pretend hooks remain trustworthy after an unreviewed engine update.
SOURCE_HASHES = {
    "generate.py": "aba4eae926d62ddb22caa0f162224fc21e1e65e18175e5980c9b752f6bdcaa34",
    "proteina.py": "74dae63c857e1e325546161ff1aa66367de191d7b1de9ec8e37312290d7dd716",
}


def emit(path, **record):
    try:
        payload = (json.dumps({"schema_version": 1, "at": datetime.now(timezone.utc).isoformat(),
                               "pid": os.getpid(), **record}) + "\n").encode()
        fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o600)
        try:
            os.write(fd, payload)
        finally:
            os.close(fd)
    except OSError:
        pass  # Diagnostics cannot turn a scientific success into failure.


@contextmanager
def phase(path, name):
    span = str(uuid4())
    start = time.perf_counter()
    emit(path, phase=name, span_id=span, event="started")
    outcome = {}
    try:
        yield
    except BaseException as exc:
        outcome["exception_type"] = type(exc).__name__
        raise
    finally:
        emit(path, phase=name, span_id=span, event="finished",
             elapsed_seconds=time.perf_counter() - start, **outcome)


def instrumented_command(command, launcher):
    """Preserve native Hydra argv; never redirect evaluation/scoring to generation."""
    if not command or Path(command[0]).name != "rfd4-generate":
        raise ValueError("Native timing supports only rfd4-generate.")
    return [sys.executable, str(launcher), *command[1:]]


def verify_sources(directory):
    for name, expected in SOURCE_HASHES.items():
        if hashlib.sha256((Path(directory) / name).read_bytes()).hexdigest() != expected:
            raise RuntimeError(f"Native diagnostics require reviewed RFD4 source: {name} differs.")


@contextmanager
def installed_hooks(native, path):
    """Wrap known calls transparently and restore even if native execution fails."""
    restorations = []
    def hook(owner, name, label, condition=None):
        original = getattr(owner, name)
        @wraps(original)
        def timed(*args, **kwargs):
            if condition is not None and not condition():
                return original(*args, **kwargs)
            with phase(path, label):
                return original(*args, **kwargs)
        setattr(owner, name, timed)
        restorations.append((owner, name, original))
    try:
        hook(native, "load_ckpt_n_configure_inference", "model_checkpoint_and_configuration")
        hook(native, "get_dataset_from_cfg", "dataset_preparation")
        hook(native, "setup_run_compile", "compile_setup")
        hook(native.Proteina, "reserve_shape_cache", "compile_shape_cache_setup")
        hook(native.L.Trainer, "__init__", "trainer_construction")
        hook(native.L.Trainer, "validate", "validation_loop_including_setup_inference_metrics_outputs")
        hook(native.Proteina, "predict_step", "native_predict_step_including_decoding")
        hook(native.Proteina, "_save_and_log_csv", "final_metrics_csv_serialization")
        hook(native.torch.cuda, "_lazy_init", "cuda_lazy_initialization",
             condition=lambda: not native.torch.cuda.is_initialized())
        # Futures may perform metrics and CIF writes in separate processes. Timing
        # drain measures waiting, not total worker serialization/CPU time.
        from rfproteina.metrics.utils.metrics_utils import MetricsPool
        hook(MetricsPool, "drain_generated_structure", "metrics_and_structure_output_drain")
        emit(path, event="hooks_installed", cuda_initialized_after_imports=native.torch.cuda.is_initialized(),
             source_hashes=SOURCE_HASHES, gpu_timing="host_wall_without_added_synchronization")
        yield
    finally:
        for owner, name, original in reversed(restorations):
            setattr(owner, name, original)


def main():
    path = os.environ.get("MN_STUDIO_RFD4_NATIVE_TIMING_FILE")
    if not path:
        raise RuntimeError("Native diagnostics require an explicit timing output path.")
    with phase(path, "native_entrypoint_main"):
        # Match script_utils.entrypoints.generate's warning policy and entrypoint.
        with phase(path, "native_entrypoint_imports"):
            with warnings.catch_warnings(record=True) as caught:
                import rfproteina.generate as native
                from rfproteina.utils.verbosity import queue_startup_warnings
            queue_startup_warnings(caught)
        verify_sources(Path(native.__file__).parent)
        with installed_hooks(native, path):
            native.main_hydra()
    emit(path, event="native_entrypoint_returned")
    # Interpreter startup and atexit/process teardown are outside this scope.


if __name__ == "__main__":
    main()
