"""Scheduler-owned SMILES-to-PDBQT batch measurements, without docking."""
from __future__ import annotations

import hashlib
import json
import shutil
import tempfile
from pathlib import Path

from mn_molecular_design_studio.docking_preparation_protocols import scrub_invocation
from mn_molecular_design_studio.engine.workflow_engine.core.docker_runner import run_docker_job
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import create_job, finish_job, read_json
from mn_molecular_design_studio.engine.workflow_engine.core.manifests import load_manifest


_RUNNER = r'''
import concurrent.futures
import hashlib
import json
import math
import subprocess
import tarfile
import time
from importlib.metadata import version
from pathlib import Path
from rdkit import Chem, rdBase

root = Path('/bench')
config = json.loads((root / 'config.json').read_text())
workers = config['cpu_cores']
route = config['route']
version_check = subprocess.check_output([
    route['python_executable'], '-c',
    'from importlib.metadata import version; print(version("molscrub"))',
], text=True).strip()
if version_check != route['version']:
    raise RuntimeError(f'Expected MolScrub {route["version"]}, found {version_check}')

inputs = (root / 'input.smi').read_text().splitlines()
input_ids = [line.split()[1] for line in inputs if line.strip()]
started = time.perf_counter()
command = [route['executable'], str(root / 'input.smi'), '-o', str(root / 'prepared.sdf'),
           *route['arguments']]
with (root / 'scrub.log').open('w') as log:
    scrub = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=False)
scrub_seconds = time.perf_counter() - started
prepared = root / 'prepared.sdf'
if not prepared.is_file() or not prepared.stat().st_size:
    raise RuntimeError(f'MolScrub produced no SDF (exit {scrub.returncode})')

split_start = time.perf_counter()
count = len(Chem.SDMolSupplier(str(prepared), removeHs=False, sanitize=True))
chunk_size = max(1, math.ceil(count / workers))
chunk_root = root / 'chunks'
chunk_root.mkdir()
paths, writer, in_chunk, invalid = [], None, 0, 0
for mol in Chem.SDMolSupplier(str(prepared), removeHs=False, sanitize=True, strictParsing=False):
    if mol is None:
        invalid += 1
        continue
    if writer is None or in_chunk >= chunk_size:
        if writer is not None:
            writer.close()
        path = chunk_root / f'chunk-{len(paths):04d}.sdf'
        paths.append(path)
        writer = Chem.SDWriter(str(path))
        in_chunk = 0
    writer.write(mol)
    in_chunk += 1
if writer is not None:
    writer.close()
split_seconds = time.perf_counter() - split_start

def meeko_chunk(item):
    index, path = item
    output = root / 'pdbqt' / f'chunk-{index:04d}'
    output.mkdir(parents=True)
    with (root / f'meeko-{index:04d}.log').open('w') as log:
        result = subprocess.run(['mk_prepare_ligand.py', '-i', str(path),
                                 '--multimol_outdir', str(output)],
                                stdout=log, stderr=subprocess.STDOUT, check=False)
    return {'chunk': index, 'return_code': result.returncode,
            'pdbqt_count': sum(1 for _ in output.glob('*.pdbqt'))}

meeko_start = time.perf_counter()
with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
    chunks = list(pool.map(meeko_chunk, enumerate(paths)))
meeko_seconds = time.perf_counter() - meeko_start
elapsed = time.perf_counter() - started
outputs = list((root / 'pdbqt').rglob('*.pdbqt'))
output_ids = {path.stem for path in outputs if path.stat().st_size}
failed_ids = [name for name in input_ids if name not in output_ids]
summary = {
    'input_count': len(input_ids), 'input_sha256': hashlib.sha256((root / 'input.smi').read_bytes()).hexdigest(),
    'cpu_workers': workers, 'scrub_version': version_check, 'rdkit_version': rdBase.rdkitVersion,
    'meeko_version': version('meeko'), 'scrub_command': command,
    'scrub_return_code': scrub.returncode, 'scrub_seconds': scrub_seconds,
    'scrub_sdf_records': count, 'invalid_sdf_records': invalid,
    'meeko_split_seconds': split_seconds, 'meeko_seconds': meeko_seconds,
    'meeko_split_and_parallel_seconds': split_seconds + meeko_seconds,
    'end_to_end_seconds': elapsed, 'output_pdbqt_count': len(output_ids),
    'failed_input_count': len(failed_ids), 'failed_ids': failed_ids,
    'scrub_records_per_second': count / scrub_seconds,
    'meeko_pdbqts_per_second': len(output_ids) / (split_seconds + meeko_seconds),
    'end_to_end_pdbqts_per_second': len(output_ids) / elapsed,
    'retry_count': 0, 'chunks': chunks,
    'timing_scope': 'MolScrub + SDF splitting + parallel Meeko; excludes Docker startup, validation and archiving',
    'working_directory': str(root),
}
archive_start = time.perf_counter()
with tarfile.open(root / 'prepared-pdbqt.tar.gz', 'w:gz', compresslevel=1) as archive:
    archive.add(root / 'pdbqt', arcname='pdbqt')
summary['archive_seconds'] = time.perf_counter() - archive_start
(root / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
print(json.dumps({key: value for key, value in summary.items() if key not in ('failed_ids', 'chunks')}), flush=True)
if not output_ids:
    raise RuntimeError('No compounds produced PDBQT output')
'''


def run_ligand_preparation_benchmark(
    smiles_file: Path,
    expected_sha256: str,
    cpu_cores: int = 23,
    expected_count: int = 10_000,
    scrub_protocol: str = "molscrub_0_1_0",
    job_name: str = "10,000-compound preparation benchmark",
    gpu_device: str = "none",
) -> Path:
    """Use one bounded Scrub pool, then parallel Meeko chunks; count failures once."""
    if type(cpu_cores) is not int or cpu_cores < 1:
        raise ValueError("The preparation job needs a positive CPU allocation.")
    source = Path(smiles_file).resolve(strict=True)
    data = source.read_bytes()
    if hashlib.sha256(data).hexdigest() != expected_sha256:
        raise ValueError("The preparation input checksum does not match the persisted benchmark sample.")
    rows = [line.split() for line in data.decode().splitlines() if line.strip()]
    if len(rows) != expected_count or any(len(row) != 2 for row in rows):
        raise ValueError("The benchmark requires the expected number of SMILES/ID rows.")
    if len({row[1] for row in rows}) != len(rows):
        raise ValueError("Benchmark compound IDs must be unique.")
    route = scrub_invocation(scrub_protocol, ph=7.4, skip_tautomers=True, cpu=cpu_cores)
    job = create_job("ligand-preparation", "ligand_preparation", "ligand_preparation_benchmark",
                     {"smiles_file": str(source)},
                     {"cpu_cores": cpu_cores, "gpu_device": "none", "compound_count": len(rows),
                      "scrub_protocol": scrub_protocol, "job_name": job_name, "retry_failed": False})
    artifacts = job.run_dir / "artifacts"
    (artifacts / "input.smi").write_bytes(data)
    try:
        with tempfile.TemporaryDirectory(prefix="mn-pdbqt-benchmark-", dir="/tmp") as folder:
            scratch = Path(folder)
            (scratch / "input.smi").write_bytes(data)
            (scratch / "runner.py").write_text(_RUNNER)
            (scratch / "config.json").write_text(json.dumps({"cpu_cores": cpu_cores, "route": route}))
            return_code = run_docker_job(
                job.run_dir, load_manifest("unidock_pro"), ["python", "/bench/runner.py"],
                gpu_device="none", extra_mounts=[(scratch, "/bench", False)],
                environment={name: "1" for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                                                       "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS")},
                use_host_user=True,
            )
            for file in scratch.iterdir():
                if file.suffix == ".log" or file.name in {"summary.json", "config.json", "runner.py", "prepared-pdbqt.tar.gz"}:
                    shutil.copy2(file, artifacts / file.name)
            summary = read_json(artifacts / "summary.json")
        if return_code or not summary.get("output_pdbqt_count"):
            raise RuntimeError(f"Preparation benchmark did not finish successfully (exit {return_code}); inspect job logs.")
        finish_job(job.run_dir, True, {
            "metrics": {key: value for key, value in summary.items() if key not in ("chunks", "failed_ids", "scrub_command")},
            "outputs": {"summary": "artifacts/summary.json", "prepared_compounds": "artifacts/prepared-pdbqt.tar.gz"},
        })
        return job.run_dir
    except Exception as exc:
        finish_job(job.run_dir, False, {"error": str(exc), "outputs": {}})
        raise
