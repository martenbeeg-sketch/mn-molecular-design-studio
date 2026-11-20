from __future__ import annotations

import hashlib
import re
import shutil
from pathlib import Path

from mn_molecular_design_studio import cofolding_msa
from mn_molecular_design_studio.engine.workflow_engine.runtime import reference_root


AA3_TO_1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}

_A3M_SEQUENCE_ALLOWED = set("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz-")
BOLTZ_CACHE_DIR = Path("/mnt/db/reference_files/boltz_models")
BOLTZ_MSA_REPOSITORY_DIR = BOLTZ_CACHE_DIR / "msa_repository"


def _clean_key(text: str, fallback: str = "target") -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(text or "").strip()).strip("_")
    return cleaned or fallback


def pdb_residue_records_by_chain(path: Path) -> dict[str, list[tuple[int, str, str]]]:
    records: dict[str, list[tuple[int, str, str]]] = {}
    seen: set[tuple[str, str, str]] = set()
    for line in path.read_text(errors="ignore").splitlines():
        if not line.startswith("ATOM  "):
            continue
        chain = line[21].strip() or "_"
        residue_key = (chain, line[22:26].strip(), line[26].strip())
        if residue_key in seen:
            continue
        seen.add(residue_key)
        try:
            residue_number = int(line[22:26])
        except ValueError:
            continue
        records.setdefault(chain, []).append(
            (residue_number, line[26].strip(), AA3_TO_1.get(line[17:20].strip().upper(), "X"))
        )
    return {
        chain: sorted(chain_records, key=lambda item: (item[0], item[1]))
        for chain, chain_records in records.items()
    }


def target_chain_sequences(path: Path, target_chains: list[str] | None = None) -> dict[str, str]:
    residue_records = pdb_residue_records_by_chain(path)
    requested = [chain for chain in target_chains or [] if chain in residue_records]
    chains = requested or list(residue_records)
    return {chain: "".join(record[2] for record in residue_records.get(chain, [])) for chain in chains}


def boltz_msa_repository() -> Path:
    return reference_root() / "boltz_models" / "msa_repository"


def boltz_msa_paths(sequence: str, msa_repository_dir: Path | None = None) -> tuple[Path, str]:
    msa_repository_dir = msa_repository_dir or boltz_msa_repository()
    digest = hashlib.sha256(sequence.encode("utf-8")).hexdigest()
    filename = f"{digest}.a3m"
    return msa_repository_dir / filename, f"/msa_repository/{filename}"


def _normalize_msa_query_sequence(sequence: str) -> str:
    return re.sub(r"[^A-Za-z]+", "", str(sequence or "")).upper()


def _first_a3m_query_sequence(path: Path) -> str:
    sequence_lines: list[str] = []
    in_first_record = False
    for raw_line in path.read_text(errors="ignore").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if in_first_record:
                break
            in_first_record = True
            continue
        if in_first_record:
            sequence_lines.append(line.replace("-", ""))
    return _normalize_msa_query_sequence("".join(sequence_lines))


def _read_a3m_records(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    sequence_lines: list[str] = []
    for raw_line in path.read_text(errors="ignore").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(sequence_lines)))
            header = line
            sequence_lines = []
        else:
            sequence_lines.append(line)
    if header is not None:
        records.append((header, "".join(sequence_lines)))
    return records


def _a3m_query_match_mask(query_sequence: str, residue_count: int) -> list[bool] | None:
    keep: list[bool] = []
    consumed = 0
    for char in query_sequence:
        is_match_column = char.isalpha() and char.isupper()
        if is_match_column and consumed >= residue_count:
            keep.append(False)
            continue
        keep.append(True)
        if is_match_column:
            consumed += 1
    if consumed < residue_count:
        return None
    return keep


def _write_prefix_cropped_a3m(source_path: Path, destination_path: Path, sequence: str) -> bool:
    records = _read_a3m_records(source_path)
    if not records:
        return False
    query = _normalize_msa_query_sequence(records[0][1].replace("-", ""))
    expected = _normalize_msa_query_sequence(sequence)
    if not query.startswith(expected) or len(query) == len(expected):
        return False
    keep_mask = _a3m_query_match_mask(records[0][1], len(expected))
    if keep_mask is None:
        return False
    lines: list[str] = []
    for header, aligned_sequence in records:
        cropped = "".join(char for index, char in enumerate(aligned_sequence) if index < len(keep_mask) and keep_mask[index])
        if not cropped:
            continue
        lines.append(header)
        lines.append(cropped)
    if len(lines) < 2:
        return False
    destination_path.parent.mkdir(parents=True, exist_ok=True)
    destination_path.write_text("\n".join(lines) + "\n")
    valid, _reason = validate_a3m_file(destination_path)
    if not valid or _first_a3m_query_sequence(destination_path) != expected:
        destination_path.unlink(missing_ok=True)
        return False
    return True


def validate_a3m_file(path: Path) -> tuple[bool, str]:
    if not path.exists():
        return False, "file does not exist"
    try:
        raw = path.read_bytes()
    except Exception as exc:
        return False, f"read failed: {exc}"
    if not raw:
        return False, "file is empty"
    if b"\x00" in raw:
        return False, "contains NUL byte(s)"
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as exc:
        return False, f"invalid UTF-8: {exc}"
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines:
        return False, "contains no non-empty lines"
    if not lines[0].startswith(">"):
        return False, "first non-empty line is not a FASTA header"
    saw_sequence = False
    for line in lines:
        if line.startswith(">"):
            continue
        saw_sequence = True
        if any(char not in _A3M_SEQUENCE_ALLOWED for char in line):
            return False, "contains invalid sequence characters"
    if not saw_sequence:
        return False, "contains headers only and no sequence"
    return True, ""


def find_cached_msa_for_sequence(
    sequence: str,
    msa_repository_dir: Path = BOLTZ_MSA_REPOSITORY_DIR,
) -> tuple[Path | None, str]:
    """Find an exact-sequence A3M, preferring the shared MSA cache.

    The shared ``reference_files/msa_cache`` is checked before the legacy Boltz
    repository. The optional repository is retained as a compatibility fallback
    for older reference files and benchmark-specific repositories.
    """
    cleaned_sequence = "".join(str(sequence or "").split()).upper()
    if not cleaned_sequence:
        return None, "empty sequence"
    shared_msa = cofolding_msa.find_reference(cleaned_sequence, reference=Path(reference_root()))
    if shared_msa is not None and shared_msa.is_file():
        shared_cache = (Path(reference_root()) / "msa_cache").resolve()
        try:
            shared_msa.resolve().relative_to(shared_cache)
            source = "shared_msa_cache"
        except ValueError:
            source = "legacy_reference_repository"
        return shared_msa, source

    host_path, _container_path = boltz_msa_paths(cleaned_sequence, msa_repository_dir=msa_repository_dir)
    valid, reason = validate_a3m_file(host_path)
    if valid:
        return host_path, "sequence_hash"
    if host_path.exists():
        repaired = _repair_a3m_file(host_path)
        if repaired:
            return host_path, "sequence_hash_repaired"
    if not Path(msa_repository_dir).exists():
        return None, f"repository missing; hash path: {reason}"

    expected_query = _normalize_msa_query_sequence(cleaned_sequence)
    prefix_sources: list[Path] = []
    for candidate in sorted(Path(msa_repository_dir).glob("**/*.a3m")):
        if candidate == host_path:
            continue
        valid, _candidate_reason = validate_a3m_file(candidate)
        if not valid and not _repair_a3m_file(candidate):
            continue
        try:
            candidate_query = _first_a3m_query_sequence(candidate)
            if candidate_query == expected_query:
                return candidate, "sequence_scan"
            if candidate_query.startswith(expected_query):
                prefix_sources.append(candidate)
        except Exception:
            continue
    for candidate in sorted(prefix_sources, key=lambda path: (len(_first_a3m_query_sequence(path)), str(path))):
        if _write_prefix_cropped_a3m(candidate, host_path, expected_query):
            return host_path, f"sequence_prefix_crop:{candidate.name}"
    return None, reason or "no sequence-matched A3M in repository"


def _validate_a3m_payload(raw: bytes) -> tuple[bool, str]:
    temp = Path("/tmp") / f"mn_molecular_design_studio.engine.workflow_engine_a3m_check_{hashlib.sha256(raw).hexdigest()}.a3m"
    try:
        temp.write_bytes(raw)
        return validate_a3m_file(temp)
    finally:
        temp.unlink(missing_ok=True)


def _repair_a3m_file(path: Path) -> bool:
    if not path.exists():
        return False
    cleaned = path.read_bytes().rstrip(b"\x00").replace(b"\r\n", b"\n")
    valid, _reason = _validate_a3m_payload(cleaned)
    if not valid:
        return False
    if cleaned != path.read_bytes():
        path.write_bytes(cleaned)
    return True


def ensure_boltz_msa_for_sequence(
    run_dir: Path,
    sequence: str,
    label: str,
    raw_subdir: str = "target_msa",
    gpu_device: object = "0",
) -> str:
    host_path, container_path = boltz_msa_paths(sequence)
    # The client stages the full alignment so a later legacy run can reuse it.
    alignment = cofolding_msa.prepare_chain(
        run_dir, sequence, _clean_key(label), item_id=raw_subdir, depth=None,
    )
    text = alignment.read_text()
    if not host_path.is_file() or host_path.read_text(errors="replace") != text:
        cofolding_msa._atomic_text(host_path, text)
    return container_path


def ensure_boltz_msas_for_target(
    run_dir: Path,
    target_artifact: Path,
    target_chains: list[str] | None = None,
    raw_subdir: str = "target_msa",
    gpu_device: object = "0",
) -> dict[str, str]:
    sequences = target_chain_sequences(target_artifact, target_chains)
    msa_paths: dict[str, str] = {}
    for chain, sequence in sequences.items():
        if not sequence:
            continue
        msa_paths[chain] = ensure_boltz_msa_for_sequence(
            run_dir,
            sequence,
            f"target_chain_{chain}",
            raw_subdir=raw_subdir,
            gpu_device=gpu_device,
        )
    return msa_paths


def ensure_pxdesign_msa_dirs_for_target(
    run_dir: Path,
    target_artifact: Path,
    target_chains: list[str] | None = None,
    raw_subdir: str = "pxdesign/input/msa",
    yaml_base_dir: Path | None = None,
    gpu_device: object = "0",
) -> dict[str, str]:
    """Create PXDesign-compatible MSA directories from the shared A3M cache.

    PXDesign/Protenix expects each chain MSA as a directory containing
    ``pairing.a3m`` and ``non_pairing.a3m``. The app keeps a single
    sequence-hashed A3M archive for reuse across workflows, then materializes
    per-run directories next to the PXDesign input YAML.
    """

    sequences = target_chain_sequences(target_artifact, target_chains)
    msa_dirs: dict[str, str] = {}
    yaml_base_dir = yaml_base_dir or run_dir
    for chain, sequence in sequences.items():
        if not sequence:
            continue
        ensure_boltz_msa_for_sequence(
            run_dir,
            sequence,
            f"target_chain_{chain}",
            raw_subdir=raw_subdir,
            gpu_device=gpu_device,
        )
        host_msa, _container_msa = boltz_msa_paths(sequence)
        valid, reason = validate_a3m_file(host_msa)
        if not valid:
            raise RuntimeError(f"Cached target MSA for chain {chain} is invalid: {reason}.")
        chain_dir = run_dir / "artifacts" / "raw" / raw_subdir / _clean_key(chain, fallback="chain")
        chain_dir.mkdir(parents=True, exist_ok=True)
        for filename in ["pairing.a3m", "non_pairing.a3m"]:
            shutil.copyfile(host_msa, chain_dir / filename)
        try:
            msa_dirs[chain] = str(chain_dir.relative_to(yaml_base_dir))
        except ValueError:
            msa_dirs[chain] = str(chain_dir)
    return msa_dirs
