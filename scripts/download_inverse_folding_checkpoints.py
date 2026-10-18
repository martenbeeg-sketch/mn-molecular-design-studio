#!/usr/bin/env python3
"""Fetch pinned Caliby assets into shared references; never run model inference."""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import urllib.request
from pathlib import Path

REPO = "ProteinDesignLab/caliby-weights"
REVISION = "a51f011f4ec7ffe2daad3ab9b2bfb67ff628a096"


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def install(root):
    destination = root / "caliby"
    destination.mkdir(parents=True, exist_ok=True)
    url = f"https://huggingface.co/api/models/{REPO}/tree/{REVISION}?recursive=true&expand=false"
    files = json.load(urllib.request.urlopen(url, timeout=60))
    records = []
    for row in files:
        name = row["path"]
        if row["type"] != "file" or not name.startswith(("caliby/", "proteinmpnn/", "protpardelle-1c/")):
            continue
        target = destination / name
        expected = row.get("lfs", {}).get("oid")
        target.parent.mkdir(parents=True, exist_ok=True)
        if not target.exists() and name.startswith("protpardelle-1c/"):
            cached = root / "protpardelle-1c/model_params" / name.removeprefix("protpardelle-1c/")
            if cached.is_file() and (not expected or sha256(cached) == expected):
                shutil.copy2(cached, target)
        if not target.exists() or expected and sha256(target) != expected:
            temporary = target.with_suffix(target.suffix + ".part")
            source = f"https://huggingface.co/{REPO}/resolve/{REVISION}/{name}"
            print(f"Downloading {name}", flush=True)
            with urllib.request.urlopen(source, timeout=120) as response, temporary.open("wb") as output:
                shutil.copyfileobj(response, output, length=1024 * 1024)
            if expected and sha256(temporary) != expected:
                raise ValueError(f"Checksum mismatch: {name}")
            temporary.replace(target)
        digest = sha256(target)
        records.append(dict(path=name, bytes=target.stat().st_size, sha256=digest))
        print(f"Verified {name}", flush=True)
    (destination / "checkpoints.json").write_text(json.dumps(dict(repo=REPO, revision=REVISION, files=records), indent=2) + "\n")
    return records


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-dir", type=Path, default=Path("/mnt/db/reference_files"))
    install(parser.parse_args().reference_dir)
