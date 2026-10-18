"""Shared, exact-query MSA acquisition and job-local provenance for cofolding."""
from __future__ import annotations

import fcntl
import hashlib
import json
import os
import re
import subprocess
import tempfile
from pathlib import Path

from .engine.workflow_engine.runtime import reference_root


def records(text: str) -> list[tuple[str, str]]:
    result = []
    header, sequence = '', []
    for line in text.replace('\x00', '').splitlines():
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        if line.startswith('>'):
            if header:
                result.append((header, ''.join(sequence)))
            header, sequence = line, []
        elif header:
            sequence.append(line)
        else:
            raise ValueError('A3M must begin with a FASTA header.')
    if header:
        result.append((header, ''.join(sequence)))
    return result


def validated_text(text: str, sequence: str, depth: int | None = None) -> str:
    if '\x00' in text:
        raise ValueError('A3M contains NUL bytes.')
    rows = records(text)
    if not rows or re.sub(r'[a-z.]', '', rows[0][1]) != sequence:
        raise ValueError('The A3M query does not match the protein sequence exactly.')
    for _, row in rows:
        aligned = re.sub(r'[a-z.]', '', row)
        if len(aligned) != len(sequence) or re.search(r'[^A-Za-z.\-]', row):
            raise ValueError('A3M rows must have the same aligned length as the query.')
    return ''.join(f'{header}\n{row}\n' for header, row in rows[:depth])


def _atomic_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f'.{path.name}.', dir=path.parent)
    try:
        with os.fdopen(fd, 'w', encoding='utf-8') as handle:
            handle.write(text)
        os.replace(temporary, path)
    finally:
        Path(temporary).unlink(missing_ok=True)


def _query(path: Path) -> str:
    # Read only the first record. Scanning every complete cached alignment made
    # a cache miss unnecessarily expensive on the shared reference database.
    query, started = [], False
    with path.open() as handle:
        for line in handle:
            if line.startswith('>'):
                if started:
                    break
                started = True
            elif started:
                query.append(line.strip())
    return ''.join(query)


def find_reference(sequence: str, reference: Path | None = None) -> Path | None:
    root = reference or Path(reference_root())
    digest = hashlib.sha256(sequence.encode()).hexdigest()
    repositories = (root / 'msa_cache', root / 'boltz_models/msa_repository')
    prefixes = []
    for repository in repositories:
        canonical = repository / f'{digest}.a3m'
        candidates = [canonical, *sorted(repository.rglob('*.a3m'))] if repository.is_dir() else []
        for path in dict.fromkeys(candidates):
            try:
                query = _query(path)
                if query == sequence:
                    validated_text(path.read_text(), sequence)
                    return path
                if query.startswith(sequence):
                    prefixes.append((len(query), path, query))
            except (OSError, UnicodeError, ValueError):
                continue
    # Preserve the existing target-prefix reuse behavior for cropped targets.
    # Crop by aligned columns, retaining lowercase insertions correctly.
    for _, source, query in sorted(prefixes):
        try:
            text = validated_text(source.read_text(), query)
            cropped = []
            for header, row in records(text):
                columns, letters = 0, []
                for char in row:
                    if char.isupper() or char == '-':
                        if columns == len(sequence):
                            break
                        columns += 1
                    letters.append(char)
                cropped.append((header, ''.join(letters)))
            text = validated_text(''.join(f'{header}\n{row}\n' for header, row in cropped), sequence)
            destination = root / 'msa_cache' / f'{digest}.a3m'
            _atomic_text(destination, text)
            return destination
        except (OSError, UnicodeError, ValueError):
            continue
    return None


def fetch_colabfold(run_dir: Path, sequence: str, chain_id: str) -> str:
    """Run only the CPU/server MSA step; no prediction or GPU reservation."""
    digest = hashlib.sha256(sequence.encode()).hexdigest()
    root = run_dir / 'artifacts/raw/cofolding_msa/colabfold_server' / digest
    output = root / 'output'
    output.mkdir(parents=True, exist_ok=True)
    query = root / 'query.fasta'
    query.write_text(f'>chain_{chain_id}\n{sequence}\n')
    docker = os.getenv('MN_MOLECULAR_DESIGN_STUDIO_DOCKER', os.getenv('MN_COFOLDING_DOCKER', 'docker'))
    image = os.getenv('MN_COFOLDING_COLABFOLD_MSA_IMAGE', 'mn-colabfold-opt:1.6.1-jax0.5.3-cu124')
    command = [docker, 'run', '--rm', '--mount', f'type=bind,src={run_dir.resolve()},dst=/work',
               '-w', '/work', image, 'colabfold_batch',
               f'/work/{query.relative_to(run_dir).as_posix()}', f'/work/{output.relative_to(run_dir).as_posix()}',
               '--msa-only', '--af3-json', '--msa-mode', 'mmseqs2_uniref_env',
               '--pair-mode', 'unpaired_paired', '--pair-strategy', 'greedy']
    with (run_dir / 'stdout.log').open('a') as stdout, (run_dir / 'stderr.log').open('a') as stderr:
        stdout.write(f'ColabFold server MSA cache miss for chain {chain_id}.\n$ {" ".join(command)}\n')
        stdout.flush()
        completed = subprocess.run(command, stdout=stdout, stderr=stderr, timeout=3600, check=False)
    if completed.returncode:
        raise RuntimeError(f'ColabFold MSA search failed for chain {chain_id} (exit {completed.returncode}). See job logs.')
    for path in sorted(output.glob('*.json')):
        payload = json.loads(path.read_text())
        for entity in payload.get('sequences', []):
            protein = entity.get('protein', {})
            query_sequence = protein.get('sequence', '')
            if isinstance(query_sequence, list):
                query_sequence = ''.join(query_sequence)
            if query_sequence == sequence and protein.get('unpairedMsa'):
                return validated_text(protein['unpairedMsa'], sequence)
    raise RuntimeError(f'ColabFold returned no matching unpaired alignment for chain {chain_id}.')


def prepare_chain(run_dir: Path, sequence: str, chain_id: str, *, item_id: str = '',
                  depth: int = 1024, supplied: str = '', supplied_source: str = 'custom', allow_server: bool = True,
                  reference: Path | None = None, staged_name: str | None = None) -> Path:
    """Resolve/cache a full alignment, then stage the exact capped input in the job."""
    if depth < 1:
        raise ValueError('MSA depth must be positive.')
    root = reference or Path(reference_root())
    digest = hashlib.sha256(sequence.encode()).hexdigest()
    cache = root / 'msa_cache' / f'{digest}.a3m'
    cache.parent.mkdir(parents=True, exist_ok=True)
    with (cache.parent / f'.{digest}.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        existing = find_reference(sequence, root)
        source = supplied_source if supplied else 'reference_cache' if existing else 'colabfold_server'
        if supplied:
            text = validated_text(supplied, sequence)
        elif existing:
            text = validated_text(existing.read_text(), sequence)
        elif allow_server:
            text = validated_text(fetch_colabfold(run_dir, sequence, chain_id), sequence)
        else:
            raise ValueError(f'No reference MSA for chain {chain_id}; enable MSA search or supply an A3M.')
        # Keep the full reusable alignment; never replace a deeper cache with a capped upload.
        previous = existing.read_text() if existing else ''
        if not existing or len(records(text)) > len(records(previous)):
            _atomic_text(cache, text)
        elif existing != cache:
            _atomic_text(cache, previous)
    safe_item = re.sub(r'[^A-Za-z0-9_-]', '_', item_id or 'prediction')
    safe_chain = re.sub(r'[^A-Za-z0-9_-]', '_', chain_id)
    staged_file = staged_name or f'chain_{safe_chain}.a3m'
    staged = run_dir / 'artifacts/raw/cofolding_msa/inputs' / safe_item / staged_file
    capped = validated_text(text, sequence, depth)
    _atomic_text(staged, capped)
    manifest = run_dir / 'artifacts/raw/cofolding_msa/manifest.json'
    # A workflow prepares its chains sequentially. The reference-cache lock above
    # separately serializes acquisition across scheduler jobs.
    payload = json.loads(manifest.read_text()) if manifest.is_file() else {'version': 1, 'alignments': []}
    rows = payload['alignments']
    rows[:] = [row for row in rows if (row['item_id'], row['chain_id']) != (item_id, chain_id)]
    rows.append({'item_id': item_id, 'chain_id': chain_id, 'sequence_sha256': digest,
                 'job_path': staged.relative_to(run_dir).as_posix(), 'cache_path': str(cache),
                 'source': source, 'source_path': str(existing) if existing and not supplied else '',
                 'rows': len(records(capped)), 'cached_rows': len(records(cache.read_text()))})
    _atomic_text(manifest, json.dumps(payload, indent=2) + '\n')
    return staged


def complex_a3m(chains: list[dict], msas: dict[str, Path]) -> str:
    """ColabFold multimer format with unpaired rows and correct homomer cardinality."""
    unique = list(dict.fromkeys(chain['sequence'] for chain in chains))
    counts = [sum(chain['sequence'] == sequence for chain in chains) for sequence in unique]
    texts = {chain['sequence']: msas[chain['id']].read_text() for chain in chains}
    # The body contains each unique sequence once; cardinalities live only in
    # the header (ColabFold input.msa_to_str uses the same convention).
    expanded = unique
    lines = ['#' + ','.join(str(len(s)) for s in unique) + '\t' + ','.join(map(str, counts)),
             '>query', ''.join(expanded)]
    for index, sequence in enumerate(expanded):
        for header, row in records(texts[sequence]):
            lines.extend([header, ''.join(row if i == index else '-' * len(s) for i, s in enumerate(expanded))])
    return '\n'.join(lines) + '\n'
