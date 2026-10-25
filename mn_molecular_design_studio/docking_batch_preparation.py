"""Container-side batch preparation used by the UniDock Studio adapter.

The process pool belongs to MolScrub. Meeko subsequently runs one CLI per SDF
chunk, with the same CPU budget. Valid partial output is kept without retries.
"""

PREPARATION_RUNNER = r'''
import concurrent.futures
import json
import math
import subprocess
import time
from pathlib import Path
from rdkit import Chem, rdBase

root = Path('artifacts')
config = json.loads((root / 'input/preparation_config.json').read_text())
route, workers = config['route'], config['workers']
installed = subprocess.check_output([route['python_executable'], '-c',
    'from importlib.metadata import version; print(version("molscrub"))'], text=True).strip()
if installed != route['version']:
    raise RuntimeError(f'Expected MolScrub {route["version"]}, found {installed}')
prepared = root / 'prepared/batch.sdf'
started = time.perf_counter()
print(f'Preparing the complete batch with MolScrub {installed}, {workers} CPU workers; no rescue attempts.', flush=True)
with (root / 'logs/scrub.log').open('w') as log:
    scrub = subprocess.run([route['executable'], str(root / 'input/library.smi'),
        '-o', str(prepared), *route['arguments']], stdout=log, stderr=subprocess.STDOUT, check=False)
scrub_seconds = time.perf_counter() - started
if not prepared.is_file() or not prepared.stat().st_size:
    raise RuntimeError(f'MolScrub produced no usable SDF (exit {scrub.returncode}); inspect scrub.log')
split_started = time.perf_counter()
count = len(Chem.SDMolSupplier(str(prepared), removeHs=False))
chunk_size = max(1, math.ceil(count / workers))
paths, writer, in_chunk, invalid = [], None, 0, 0
for mol in Chem.SDMolSupplier(str(prepared), removeHs=False, sanitize=True, strictParsing=False):
    if mol is None:
        invalid += 1
        continue
    if writer is None or in_chunk >= chunk_size:
        if writer is not None:
            writer.close()
        path = root / 'prepared' / f'meeko-{len(paths):04d}.sdf'
        paths.append(path)
        writer = Chem.SDWriter(str(path))
        in_chunk = 0
    writer.write(mol)
    in_chunk += 1
if writer is not None:
    writer.close()
split_seconds = time.perf_counter() - split_started

def convert(item):
    index, path = item
    output = root / 'prepared' / f'pdbqt-{index:04d}'
    output.mkdir()
    with (root / 'logs' / f'meeko-{index:04d}.log').open('w') as log:
        result = subprocess.run(['mk_prepare_ligand.py', '-i', str(path),
            '--multimol_outdir', str(output)], stdout=log, stderr=subprocess.STDOUT, check=False)
    return result.returncode

meeko_started = time.perf_counter()
with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
    return_codes = list(pool.map(convert, enumerate(paths)))
# Only the first valid output for an ID is used; no state is prepared twice.
for path in sorted((root / 'prepared').glob('pdbqt-*/*.pdbqt')):
    destination = root / 'ligands' / path.name
    if path.stat().st_size and not destination.exists():
        path.rename(destination)
meeko_seconds = time.perf_counter() - meeko_started
ids = [line.split()[1] for line in (root / 'input/library.smi').read_text().splitlines() if line.strip()]
successful = [name for name in ids if (root / 'ligands' / f'{name}.pdbqt').is_file()]
summary = dict(input_count=len(ids), output_pdbqt_count=len(successful),
    failed_input_count=len(ids)-len(successful), cpu_workers=workers,
    scrub_version=installed, rdkit_version=rdBase.rdkitVersion,
    scrub_return_code=scrub.returncode, meeko_return_codes=return_codes,
    scrub_seconds=scrub_seconds, meeko_split_seconds=split_seconds,
    meeko_seconds=meeko_seconds, end_to_end_seconds=time.perf_counter()-started,
    invalid_sdf_records=invalid, retry_count=0, failure_policy='skip',
    end_to_end_pdbqts_per_second=len(successful)/(time.perf_counter()-started))
(root / 'preparation_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
print(json.dumps(summary), flush=True)
if not successful:
    raise RuntimeError('No compounds produced PDBQT output')
'''


SCORE_RUNNER = r'''
import csv
import math
import sys
from pathlib import Path

root = Path('artifacts')
with (root / 'input/compounds.tsv').open() as source, (root / 'docking_scores.csv').open('w') as output:
    writer = csv.DictWriter(output, delimiter='\t', lineterminator='\n', fieldnames=[
        'compound_id', 'status', 'best_score_kcal_mol', 'cnn_score', 'cnn_affinity', 'pose_path', 'error'])
    writer.writeheader()
    for row in csv.DictReader(source, delimiter='\t'):
        identifier = row['compound_id']
        pose = root / 'poses' / f'{identifier}_out.pdbqt'
        score = None
        if pose.is_file():
            with pose.open() as handle:
                for line in handle:
                    if line.startswith('REMARK VINA RESULT:'):
                        score = line.split()[3]
                        break
                    if line.startswith('REMARK minimizedAffinity '):
                        score = line.split()[2]
                        break
        try:
            score = float(score)
            if not math.isfinite(score):
                score = None
        except (TypeError, ValueError):
            score = None
        writer.writerow(dict(compound_id=identifier, status='completed' if score is not None else 'failed',
            best_score_kcal_mol=score if score is not None else '',
            pose_path=str(pose) if score is not None else '',
            error='' if score is not None else f'No scored pose; preparation or docking failed (dock exit {sys.argv[1]})'))
'''
