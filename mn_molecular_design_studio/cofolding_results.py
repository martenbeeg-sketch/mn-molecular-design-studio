"""Engine output adapters for the shared cofolding results page.

Artifact names follow the version-pinned batch images and MN Cofolding's readers.
Only files associated with the selected structure are read; nearby models' scores
must never be substituted when a confidence file is absent.
"""
from __future__ import annotations

import gzip
import json
import math
import re
from pathlib import Path
from typing import Any

import numpy as np


def local_artifact(run_dir: Path, value: object) -> Path | None:
    if not value:
        return None
    path = Path(str(value))
    if str(path).startswith('/work/'):
        path = Path(str(path)[6:])
    if not path.is_absolute():
        path = run_dir / path
    try:
        path = path.resolve(strict=True)
        path.relative_to(run_dir.resolve())
    except (OSError, ValueError):
        return None
    return path if path.is_file() else None


def structure_for(run_dir: Path, candidate: dict) -> Path | None:
    raw = candidate.get('raw_metadata') or {}
    return local_artifact(run_dir, candidate.get('complex_pdb') or raw.get('structure_file'))


def confidence_files(structure: Path) -> list[Path]:
    """Exact per-model sidecars for AF3, Boltz, EF2, Chai, Atlas, PX/DDE, RF3 and AF2."""
    stem = structure.stem
    base = stem.removesuffix('_model')
    names = [f'{base}_summary_confidences.json', f'{base}_confidences_aggregated.json',
             f'confidence_{stem}.json', f'{base}_confidence.json', f'{base}_scores.json',
             f'{base}_confidences.json', f'{base}_confidence.npz', f'{stem}_pae.npz',
             f'pae_{stem}.npz', f'plddt_{stem}.npz', f'{stem}.npz', f'{stem}.json']
    sample = re.search(r'_sample_(\d+)$', stem)
    if sample:
        prefix = stem[:sample.start()]
        names += [f'{prefix}_summary_confidence_sample_{sample[1]}.json', f'{prefix}_full_data_sample_{sample[1]}.json']
    chai = re.search(r'model_idx_(\d+)$', stem)
    if chai:
        names.append(f'scores.model_idx_{chai[1]}.npz')
    for marker in ('_unrelaxed_rank_', '_relaxed_rank_'):
        if marker in stem:
            names.append(stem.replace(marker, '_scores_rank_', 1) + '.json')
    boltz_base = re.sub(r'_model_\d+$', '', stem)
    names.append(f'affinity_{boltz_base}.json')
    paths = [structure.with_name(name) for name in names]
    boltzgen_sample = re.fullmatch(r'(.+)_sample(\d+)', stem)
    if boltzgen_sample and structure.parent.name == boltzgen_sample[1]:
        paths.append(structure.parent.parent.parent / 'boltzgen_data' / 'fold_out_npz' / f'{boltzgen_sample[1]}.npz')
    # Foundry can compress its confidence JSONs.
    paths += [Path(str(path) + '.gz') for path in paths if path.suffix == '.json']
    return list(dict.fromkeys(path for path in paths if path.is_file()))


def candidate_files(run_dir: Path, candidate: dict) -> tuple[Path | None, list[Path]]:
    structure = structure_for(run_dir, candidate)
    raw = candidate.get('raw_metadata') or {}
    explicit = [local_artifact(run_dir, raw.get(key)) for key in ('pae_path', 'confidence_file', 'confidence_arrays', 'confidence_json')]
    files = list(dict.fromkeys([*(path for path in explicit if path), *(confidence_files(structure) if structure else [])]))
    return structure, [path for path in files if local_artifact(run_dir, path)]


def read_payload(path: Path) -> dict[str, Any]:
    try:
        if path.suffix == '.npz':
            with np.load(path, allow_pickle=False) as archive:
                return {key: np.asarray(archive[key]) for key in archive.files
                        if any(token in key for token in ('plddt', 'pae', 'ptm', 'score', 'affinity', 'complex', 'clash', 'pde', 'disorder', 'interaction', 'ipsae'))
                        if np.issubdtype(archive[key].dtype, np.number) or archive[key].dtype == np.bool_}
        if path.suffix == '.gz':
            with gzip.open(path, 'rt', encoding='utf-8') as handle:
                value = json.load(handle)
        else:
            value = json.loads(path.read_text(encoding='utf-8'))
        if isinstance(value, list) and len(value) == 1:
            value = value[0]
        return value if isinstance(value, dict) else {}
    except (OSError, UnicodeError, ValueError, EOFError):
        return {}


def number(value: object) -> float | None:
    try:
        if isinstance(value, dict):
            value = value.get('value')
        array = np.asarray(value)
        if array.size != 1 or array.dtype == np.bool_:
            return None
        result = float(array.reshape(-1)[0])
        return result if math.isfinite(result) else None
    except (TypeError, ValueError, OverflowError):
        return None


def summary_metrics(payload: dict) -> dict:
    values = {**payload, **(payload.get('complex') if isinstance(payload.get('complex'), dict) else {})}
    metrics = {key: numeric for key, value in values.items() if (numeric := number(value)) is not None}
    for key, aliases in {
        'plddt_mean': ('plddt_mean', 'avg_plddt', 'overall_plddt', 'complex_plddt', 'plddt'),
        'ranking_score': ('ranking_score', 'sample_ranking_score', 'confidence_score', 'aggregate_score'),
    }.items():
        for alias in aliases:
            if alias in metrics:
                metrics[key] = metrics[alias]
                break
    if 'plddt_mean' in metrics and metrics['plddt_mean'] <= 1:
        metrics['plddt_mean'] *= 100
    return metrics


def read_confidence(structure: Path, files: tuple[Path, ...]) -> dict:
    plddt = np.array([], dtype=float)
    pae = None
    plddt_kind = 'Atom'
    summary = prediction_row(structure)
    details = {}
    for path in files:
        payload = sample_payload(structure, path, read_payload(path))
        summary.update(summary_metrics(payload))
        for key in ('chain_ids', 'chain_ptm', 'chains_ptm', 'chain_iptm', 'chain_pair_iptm',
                    'pair_chains_iptm', 'per_chain_ptm', 'per_chain_pair_iptm', 'chain_pair_pae_min', 'chains', 'interfaces'):
            if key in payload:
                value = payload[key]
                details[key] = value.tolist() if isinstance(value, np.ndarray) else value
        for key in ('atom_plddts', 'atom_plddt', 'plddt', 'predicted_lddt'):
            value = payload.get(key)
            if value is None:
                continue
            try:
                array = np.asarray(value, dtype=float).reshape(-1)
            except (TypeError, ValueError):
                continue
            if array.size > 1 and np.isfinite(array).all() and array.size > plddt.size:
                plddt = array
                plddt_kind = 'Atom' if key.startswith('atom') else 'Residue / token'
        for key in ('pae', 'token_pair_pae', 'predicted_aligned_error', 'predicted_aligned_errors'):
            value = payload.get(key)
            if isinstance(value, dict):
                value = value.get('data', value.get('values'))
            if value is None:
                continue
            try:
                array = np.asarray(value, dtype=float)
                # Remove batch dimensions, retaining the square token-pair axes.
                while array.ndim > 2 and array.shape[0] == 1:
                    array = array[0]
                if array.ndim == 2 and array.shape[0] == array.shape[1] and np.isfinite(array).all():
                    pae = array
            except (ValueError, TypeError):
                continue
    # Preserve the established Boltz viewer/plot convention: its structure
    # carries per-atom confidence, while its separate pLDDT NPZ is per token.
    use_structure_plddt = any(path.name == f'confidence_{structure.stem}.json' for path in files)
    if plddt.size == 0 or use_structure_plddt:
        from Bio.PDB import MMCIFParser, PDBParser
        parser = PDBParser(QUIET=True) if structure.suffix.lower() in {'.pdb', '.ent'} else MMCIFParser(QUIET=True)
        model = parser.get_structure('prediction', str(structure))
        plddt = np.array([float(atom.bfactor) for atom in model.get_atoms() if float(atom.bfactor) >= 0])
        plddt_kind = 'Atom'
    if plddt.size and plddt.max() <= 1:
        plddt = plddt * 100
    if plddt.size and not np.any(plddt > 0):
        plddt = np.array([], dtype=float)
    if plddt.size:
        summary.setdefault('plddt_mean', float(plddt.mean()))
    sampled = pae is not None and pae.shape[0] > 512
    if sampled:
        indices = np.linspace(0, pae.shape[0] - 1, 512).round().astype(int)
        preview = pae[np.ix_(indices, indices)]
    else:
        preview = pae
    return {'plddt': plddt, 'plddt_kind': plddt_kind, 'pae': pae, 'pae_preview': preview, 'pae_sampled': sampled, 'metrics': summary, 'details': details}


def sample_payload(structure: Path, path: Path, payload: dict) -> dict:
    """Match BoltzGen's displayed CIF to its sample in the native archive."""
    if 'fold_out_npz' not in path.parts:
        return payload
    sample = re.fullmatch(r'(.+)_sample(\d+)', structure.stem)
    if sample and path.stem == sample[1]:
        # The native writer uses zero-based names; the cofolding adapter's
        # exported sample names are one-based.
        index = int(sample[2]) - (0 if structure.parent.name == 'refold_cif' else 1)
        return {key: value[index] for key, value in payload.items()
                if isinstance(value, np.ndarray) and value.ndim and 0 <= index < value.shape[0]}
    if re.sub(r'^rank\d+_', '', structure.stem) == path.stem:
        iptm = np.asarray(payload.get('design_to_target_iptm', []))
        ptm = np.asarray(payload.get('design_ptm', []))
        if iptm.ndim == 1 and iptm.size and ptm.shape == iptm.shape:
            # Pinned BoltzGen writer.py and analyze_utils.get_best_folding_sample
            # use this exact score to choose the exported refolded complex.
            index = int(np.argmax(0.8 * iptm + 0.2 * ptm))
            return {key: value[index]
                    if isinstance(value, np.ndarray) and value.ndim and value.shape[0] == iptm.size
                    else value for key, value in payload.items()}
    return payload


def prediction_row(structure: Path) -> dict:
    """ESMFold2 stores scalar scores in pred_rows.jsonl beside cif_all/."""
    for parent in list(structure.parents)[:2]:
        path = parent / 'pred_rows.jsonl'
        if not path.is_file():
            continue
        with path.open(errors='replace') as handle:
            for line in handle:
                try:
                    row = json.loads(line)
                except ValueError:
                    continue
                if isinstance(row, dict) and Path(str(row.get('pred_file') or '')).name == structure.name:
                    return summary_metrics(row)
    return {}


def enrich_candidate(run_dir: Path, candidate: dict) -> dict:
    """Read small score summaries for older runs without rewriting their outputs."""
    structure, files = candidate_files(run_dir, candidate)
    raw = dict(candidate.get('raw_metadata') or {})
    metrics = dict(candidate.get('metrics') or {})
    scores = {**(prediction_row(structure) if structure else {}), **dict(raw.get('confidence_outputs') or {})}
    for path in files:
        name = path.name
        # Large full-data and PAE files are loaded lazily only for the figures.
        if any(token in name for token in ('summary_confidence', 'confidences_aggregated', 'confidence_', '_confidence.json', '_scores_', 'scores.model_idx_', 'affinity_')) or 'fold_out_npz' in path.parts:
            payload = read_payload(path)
            scores.update(summary_metrics(sample_payload(structure, path, payload) if structure else payload))
    merged = {**summary_metrics(scores), **{key: value for key, value in metrics.items() if value is not None}}
    # Fill aliases even if the old adapter only saved the engine's native name.
    for key, value in summary_metrics(merged).items():
        merged[key] = value
    native_rank = re.search(r'_rank_(\d+)_', structure.stem) if structure else None
    if native_rank:
        merged['native_rank'] = int(native_rank[1])
    raw['confidence_outputs'] = scores
    return {**candidate, 'metrics': merged, 'raw_metadata': raw}


def ranked_candidates(candidates: list[dict]) -> list[dict]:
    # Sort within each complex; unrelated complexes must not compete for "best".
    groups: dict[str, list[dict]] = {}
    for candidate in candidates:
        groups.setdefault(complex_group_key(candidate), []).append(candidate)
    result = []
    for group in groups.values():
        key = next((name for name in ('native_rank', 'ranking_score', 'confidence_score', 'iptm', 'ptm', 'plddt_mean')
                    if any(number((item.get('metrics') or {}).get(name)) is not None for item in group)), '')
        def rank(candidate: dict) -> tuple:
            value = number((candidate.get('metrics') or {}).get(key))
            return (1, 0) if value is None else (0, value if key == 'native_rank' else -value)
        result.extend(sorted(group, key=rank))
    return result


def prediction_input_labels(run_dir: Path) -> dict[str, str]:
    """Read stable task names from the saved prediction input mapping."""
    for relative in (
        'artifacts/raw/optimized_prediction/input_mapping.json',
        'artifacts/raw/af3_jax_cofolding/input_mapping.json',
    ):
        path = local_artifact(run_dir, Path(run_dir) / relative)
        if path is None:
            continue
        try:
            payload = json.loads(path.read_text(encoding='utf-8'))
        except (OSError, UnicodeError, ValueError):
            continue
        rows = payload.get('items') if isinstance(payload, dict) else None
        labels = {}
        for row in rows if isinstance(rows, list) else []:
            if not isinstance(row, dict) or not row.get('id'):
                continue
            label = row.get('source_input_name') or row.get('name')
            if label:
                labels[str(row['id'])] = str(label)
        if labels:
            return labels
    return {}


def complex_label(
    candidate: dict,
    fallback: str = 'Complex',
    input_labels: dict[str, str] | None = None,
) -> str:
    """Return a user-facing identity label, showing both ID and name when distinct."""
    identity = compound_identity(candidate, input_labels)
    if identity["compound_id"] and identity["compound_name"] and identity["compound_id"] != identity["compound_name"]:
        return f"{identity['compound_id']} · {identity['compound_name']}"
    if identity["compound_name"] or identity["compound_id"]:
        return identity["compound_name"] or identity["compound_id"]
    raw = candidate.get('raw_metadata') if isinstance(candidate.get('raw_metadata'), dict) else {}
    reference = raw.get('validation_reference') if isinstance(raw.get('validation_reference'), dict) else {}
    return str(
        raw.get('input_name') or (input_labels or {}).get(str(raw.get('batch_item_id') or ''))
        or reference.get('candidate_id')
        or fallback
        or candidate.get('candidate_id')
    )


def compound_identity(candidate: dict, input_labels: dict[str, str] | None = None) -> dict[str, str]:
    """Resolve stable compound ID, display name, and SMILES across engine output contracts."""
    raw = candidate.get("raw_metadata") if isinstance(candidate.get("raw_metadata"), dict) else {}
    fallback = str(raw.get("source_input_name") or raw.get("input_name") or
                   (input_labels or {}).get(str(raw.get("batch_item_id") or "")) or "").strip()
    source_smiles = raw.get("compound_smiles") or raw.get("ligand_smiles") or []
    if isinstance(source_smiles, (str, dict)):
        source_smiles = [source_smiles]
    smiles = []
    if isinstance(source_smiles, list):
        for row in source_smiles:
            value = row.get("smiles") or row.get("value") or "" if isinstance(row, dict) else row
            text = str(value or "").strip()
            if text and text not in smiles:
                smiles.append(text)
    has_ligand = bool(smiles or raw.get("ligand_ccd_codes") or raw.get("compound_id") or raw.get("source_compound_id"))
    compound_id = str(raw.get("compound_id") or raw.get("source_compound_id") or (fallback if has_ligand else "")).strip()
    compound_name = str(raw.get("compound_name") or raw.get("source_compound_name") or (fallback if has_ligand else "")).strip()
    return {"compound_id": compound_id, "compound_name": compound_name, "smiles": "; ".join(smiles)}


def complex_group_key(candidate: dict) -> str:
    raw = candidate.get('raw_metadata') if isinstance(candidate.get('raw_metadata'), dict) else {}
    return str(raw.get('batch_item_id') or raw.get('source_compound_id') or raw.get('input_name') or '__single_complex__')


def best_candidates_by_complex(candidates: list[dict]) -> list[dict]:
    """Return one top-ranked structure for each independent input, in input order."""
    ranked = ranked_candidates(candidates)
    best: dict[str, dict] = {}
    for candidate in ranked:
        best.setdefault(complex_group_key(candidate), candidate)
    return list(best.values())


def model_index(candidate: dict, fallback: int) -> int:
    stem = Path(str(candidate.get('complex_pdb') or '')).stem
    boltzgen_sample = re.search(r'_sample(\d+)$', stem)
    if boltzgen_sample:
        return max(0, int(boltzgen_sample[1]) - 1)
    matches = re.findall(r'(?:_model_|sample[-_]|model_idx_|_x)(\d+)', stem)
    if matches:
        return int(matches[-1])
    metrics = candidate.get('metrics') or {}
    return max(0, int(metrics.get('prediction_sample') or metrics.get('prediction_index') or fallback + 1) - 1)


def alignment_files(run_dir: Path, candidate: dict) -> list[tuple[str, Path]]:
    raw = candidate.get('raw_metadata') or {}
    proteins = raw.get('protein_chains') or []
    if isinstance(proteins, dict):
        proteins = [{'id': key, 'sequence': value} for key, value in proteins.items()]
    if isinstance(raw.get('chain_sequences'), dict):
        proteins = [{'id': key, 'sequence': value} for key, value in raw['chain_sequences'].items()]
    queries = {str(row.get('id')): str(row.get('sequence', '')).upper() for row in proteins if isinstance(row, dict)}
    direct = []
    if raw.get('engine') == 'boltzgen_fold':
        mapping_path = run_dir / 'artifacts/raw/optimized_prediction/input_mapping.json'
        for item in read_payload(mapping_path).get('items', []):
            if item.get('id') != raw.get('batch_item_id'):
                continue
            details = item.get('input_details') or {}
            msa_path = local_artifact(run_dir, details.get('target_msa'))
            if msa_path:
                direct.append((f"{details.get('target_chain_id') or 'Target'} · A3M alignment", msa_path))
    if not queries:
        return direct
    found: dict[str, tuple[int, Path]] = {}
    for path in run_dir.glob('artifacts/raw/**/*.a3m'):
        if not local_artifact(run_dir, path):
            continue
        query = ''
        started = False
        with path.open(errors='replace') as handle:
            for line in handle:
                if line.startswith('>'):
                    if started:
                        break
                    started = True
                elif started:
                    query += ''.join(char for char in line.strip() if char.isupper())
        for chain, sequence in queries.items():
            if query != sequence:
                continue
            priority = 0 if 'output' not in path.parts else 1 if path.name == 'uniref.a3m' else 2
            if chain not in found or priority < found[chain][0]:
                found[chain] = (priority, path)
    return direct + [(f'{chain} · A3M alignment', path) for chain, (_, path) in sorted(found.items())
                     if path not in {path for _, path in direct}]


def manifest_alignment_files(run_dir: Path, candidate: dict) -> list[tuple[str, Path]]:
    """Prefer the exact capped inputs from this complex, with acquisition provenance."""
    import hashlib
    from .cofolding_msa import records

    raw = candidate.get('raw_metadata') or {}
    item_id = str(raw.get('batch_item_id') or '')
    manifest = read_payload(run_dir / 'artifacts/raw/cofolding_msa/manifest.json')
    found = []
    seen = set()
    for row in manifest.get('alignments', []):
        if row.get('item_id', '') != item_id:
            continue
        path = local_artifact(run_dir, row.get('job_path'))
        if path is None or row.get('chain_id') in seen:
            continue
        try:
            entries = records(path.read_text())
            if not entries or hashlib.sha256(entries[0][1].encode()).hexdigest() != row.get('sequence_sha256'):
                continue
        except (OSError, UnicodeError, ValueError):
            continue
        source = {'reference_cache': 'reference folder', 'custom': 'custom input',
                  'colabfold_server': 'ColabFold server', 'local_mmseqs_gpu': 'local MMseqs GPU'}.get(row.get('source'), 'job input')
        found.append((f"{row['chain_id']} · {source} · {len(entries):,} rows · saved for reuse", path))
        seen.add(row['chain_id'])
    return found


def alignment_status(run_dir: Path, candidate: dict) -> str:
    """Explain disabled MSAs without implying that an alignment went missing."""
    raw = candidate.get('raw_metadata') or {}
    engine = str(raw.get('engine') or '')
    metrics = candidate.get('metrics') or {}
    if engine in {'atlasfold', 'esmfold2_fast'} or metrics.get('model_variant') == 'fast':
        return 'This engine uses single-sequence input and does not consume an MSA.'
    for workflow in ('af3_jax_cofolding', 'esmfold2_cofolding'):
        mapping = read_payload(run_dir / f'artifacts/raw/{workflow}/input_mapping.json')
        if mapping.get('msa_source') == 'single_sequence':
            return 'This prediction was run in single-sequence mode, without an MSA. Enable multiple sequence alignment when submitting a new prediction to reuse reference alignments or fetch missing ones from ColabFold.'
    mapping = read_payload(run_dir / 'artifacts/raw/optimized_prediction/input_mapping.json')
    for item in mapping.get('items', []):
        if item.get('id') == raw.get('batch_item_id') and (item.get('input_details') or {}).get('msa_source') == 'single_sequence':
            return 'MSA use was disabled for this prediction. Enable multiple sequence alignment when submitting a new prediction.'
    if metrics.get('msa_source') == 'single_sequence':
        return 'This prediction was run in single-sequence mode, without an MSA.'
    return 'No saved A3M alignment was found in this older job. New MSA-enabled jobs save the exact per-chain input with the results and retain a reusable copy in the reference MSA cache.'
