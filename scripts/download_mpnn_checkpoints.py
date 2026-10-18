#!/usr/bin/env python3
"""Install published HyperMPNN/AbMPNN checkpoints without running inference."""
import argparse
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

HYPER_REV = "1877ea181d23fe84fb33147748e1367e5a02b6a9"
HYPER = {
    "v48_002_epoch240_hyper.pt": "e3daa9c556a61bfe3e79bab62b9341a9d9763364",
    "v48_010_epoch300_hyper.pt": "8036866da3be5dddc596c31b4ae52fc8ec6550dc",
    "v48_020_epoch300_hyper.pt": "d0ac2834151197e4ff5ebfdc66dec1f4668ef07f",
    "v48_030_epoch300_hyper.pt": "61d5ab97ac1ba7dc1f89277e40acb12ba2dab24c",
}


def install(root):
    root = Path(root) / "mpnn"
    root.mkdir(parents=True, exist_ok=True)
    sources = [(name, f"https://raw.githubusercontent.com/meilerlab/HyperMPNN/{HYPER_REV}/retrained_models/{name}", "git", digest) for name, digest in HYPER.items()]
    sources.append(("abmpnn.pt", "https://zenodo.org/api/records/8164693/files/abmpnn.pt/content", "md5", "911df5d74f2955373286c3bb3df69409"))
    manifest = {}
    for name, url, algorithm, expected in sources:
        destination = root / name
        content = destination.read_bytes() if destination.is_file() else urlopen(url, timeout=120).read()
        actual = hashlib.sha1(f"blob {len(content)}\0".encode() + content).hexdigest() if algorithm == "git" else hashlib.md5(content).hexdigest()
        if actual != expected:
            raise ValueError(f"Checksum mismatch for {name}; existing file was not overwritten.")
        if not destination.exists():
            temporary = destination.with_suffix(".download")
            temporary.write_bytes(content)
            temporary.replace(destination)
        manifest[name] = {"url": url, "bytes": len(content), "sha256": hashlib.sha256(content).hexdigest(), "source_checksum": f"{algorithm}:{expected}"}
        print(f"Verified {name} · {len(content):,} bytes", flush=True)
    (root / "checkpoints.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-dir", type=Path, required=True)
    install(parser.parse_args().reference_dir)
