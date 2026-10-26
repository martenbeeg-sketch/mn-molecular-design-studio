"""Install the checksum-pinned authors' deposit without executing its scripts."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import subprocess
import sys
import tarfile
import urllib.request

RECORD = "10391229"
SHA256 = "72e228dda952b82a24c1c9c4403bcc2d4340d6a25f4f06cb148954352d9b67c3"
URL = f"https://zenodo.org/api/records/{RECORD}/files/data.tar.gz/content"
SCRIPTS = {"apol.py", "design_after_seqthread.xml", "design_scaffold_peptide.xml",
           "get_apol_positions.py", "seqthread.py", "instructions.txt"}


def digest(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def install(root: Path, *, build: bool = True) -> None:
    destination = root / "peptide-traps" / f"sahtoe-{RECORD}"
    destination.mkdir(parents=True, exist_ok=True)
    archive = destination / "data.tar.gz"
    if not archive.exists():
        temporary = archive.with_suffix(".download")
        with urllib.request.urlopen(URL) as response, temporary.open("wb") as handle:
            while chunk := response.read(1024 * 1024):
                handle.write(chunk)
        if digest(temporary) != SHA256:
            raise ValueError("The downloaded paper archive does not match the pinned checksum.")
        temporary.replace(archive)
    if digest(archive) != SHA256:
        raise ValueError(f"Preserving unexpected archive at {archive}; its checksum does not match.")
    files = {}
    with tarfile.open(archive, "r:gz") as bundle:
        for member in bundle:
            path = PurePosixPath(member.name)
            if not member.isfile() or len(path.parts) < 2 or path.parts[0] != "supplementary_files" or ".." in path.parts:
                continue
            relative = PurePosixPath(*path.parts[1:])
            selected = len(relative.parts) == 1 and relative.name in SCRIPTS
            selected |= len(relative.parts) == 2 and relative.parts[0] == "pdb_models" and relative.suffix == ".pdb"
            if not selected:
                continue
            if member.size > 10_000_000:
                raise ValueError(f"Unexpectedly large source member: {relative}")
            data = bundle.extractfile(member).read()
            output = destination / "source" / relative
            output.parent.mkdir(parents=True, exist_ok=True)
            if output.exists() and output.read_bytes() != data:
                raise ValueError(f"Preserving locally modified paper source: {output}")
            output.write_bytes(data)
            files[str(relative)] = hashlib.sha256(data).hexdigest()
    manifest = {"doi": f"10.5281/zenodo.{RECORD}", "archive_sha256": SHA256,
                "source_url": URL, "license": "CC-BY-4.0", "files": files}
    (destination / "source_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Installed {len(files)} script/model files in {destination / 'source'}", flush=True)
    if build:
        studio = Path(__file__).resolve().parents[2]
        subprocess.run(["docker", "build", "-f", str(studio / "scripts/peptide-trap/Dockerfile"),
                        "-t", "mn-peptide-trap:sahtoe-10391229", str(studio)], check=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-dir", type=Path)
    parser.add_argument("--skip-build", action="store_true")
    arguments = parser.parse_args()
    if arguments.reference_dir is None:
        sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
        from mn_molecular_design_studio.runtime import reference_root
        arguments.reference_dir = reference_root()
    install(arguments.reference_dir.resolve(), build=not arguments.skip_build)
