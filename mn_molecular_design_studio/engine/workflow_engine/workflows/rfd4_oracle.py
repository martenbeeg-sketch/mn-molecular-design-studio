"""Job-scoped OpenFold3 bridge for the upstream RFD4 Unix-socket protocol.

The generator never receives a Docker socket. The host launches only the pinned
oracle command, for manifests inside this job's own mounted work directory.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import socketserver
import subprocess
import tempfile
import threading

from mn_molecular_design_studio.engine.workflow_engine.core.gpu import docker_gpu_args
from mn_molecular_design_studio.engine.workflow_engine.runtime import reference_root

ORACLE_IMAGE = 'mn-rfd4-oracle:598acab-cu13'
HF_REVISION = 'a06e2509c7033b010fc43614e9fd840303588502'


def model_bundle(reference=None):
    return Path(reference or reference_root()) / 'rfd4-proteina/huggingface/hub/models--nvidia--RFD4-Proteina/snapshots' / HF_REVISION


def compatible_engines(task):
    from mn_molecular_design_studio.engine.workflow_engine.workflows.optimized_prediction import ENGINE_PROFILES
    required = {'dna_binder':'dna_sequences','rna_design':'rna_sequences','enzyme':'ligand_ccd_codes','small_molecule_binder':'ligand_ccd_codes'}.get(task)
    # Models with sequence-only output or a mandatory template cannot serve as
    # independent all-atom confidence oracles. They remain available on the canvas.
    eligible={key:profile['label'] for key,profile in ENGINE_PROFILES.items()
              if not profile.get('requires_target_template') and key not in {'esmfold','esmfold2_fast','esmfold2_full','boltzgen_fold'}
              and (not required or required in profile.get('molecule_inputs',()))}
    return {'openfold3_jax':eligible.pop('openfold3_jax'),**eligible}


class OpenFold3Bridge:
    def __init__(self, run_dir, gpu_device='0', cpu_cores=8, options=None):
        self.run_dir = Path(run_dir).resolve()
        self.gpu_device = gpu_device
        self.cpu_cores = cpu_cores
        self.options = dict(options or {})
        self.engine = self.options.get("folding_engine", "openfold3_jax")
        if self.engine not in compatible_engines("existing_structure"):
            raise ValueError("Unsupported folding engine for inline RFD4 checks.")
        self.temporary = None
        self.server = None
        self.thread = None
        self.process = None
        self.log = None

    def __enter__(self):
        weights = model_bundle() / 'openfold3/weights'
        if self.engine=='openfold3_jax' and not (weights/'of3_ported_weights.provenance.json').is_file():
            raise ValueError('RFD4 OpenFold3 oracle weights are missing from the pinned model bundle.')
        from mn_molecular_design_studio.engine.workflow_engine.workflows.optimized_prediction import ENGINE_PROFILES
        image=ORACLE_IMAGE if self.engine=='openfold3_jax' else ENGINE_PROFILES[self.engine]['image']
        found = subprocess.run(['docker','image','inspect',image],capture_output=True,check=False)
        if found.returncode:
            raise ValueError('The RFD4 OpenFold3 oracle image is not installed: '+image)
        self.temporary = tempfile.TemporaryDirectory(prefix='mds-rfd4-oracle-')
        root = Path(self.temporary.name)
        runtime = self.run_dir/'artifacts/rfd4_proteina/oracle_runtime'
        shutil.copytree(Path(__file__).with_name('rfd4_oracle_runtime'),runtime,dirs_exist_ok=True,
                        ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
        bridge = self
        class Handler(socketserver.StreamRequestHandler):
            def handle(self):
                try:
                    data = self.rfile.readline(1024*1024+1)
                    if len(data)>1024*1024: raise ValueError('Oracle request too large.')
                    request = json.loads(data)
                    cmd=request.get('cmd')
                    if cmd in ('ping','prewarm','release'):
                        response={'status':'ok','released':False}
                    elif cmd=='fold':
                        bridge.fold(request.get('manifest'))
                        response={'status':'ok'}
                    else: raise ValueError('Unsupported oracle command.')
                except Exception as exc:
                    response={'status':'error','message':str(exc)}
                self.wfile.write((json.dumps(response)+'\n').encode())
        self.server=socketserver.UnixStreamServer(str(root/'oracle.sock'),Handler)
        os.chmod(root,0o700)
        self.thread=threading.Thread(target=self.server.serve_forever,daemon=True)
        self.thread.start()
        return self

    @property
    def mount(self):
        return (Path(self.temporary.name),'/oracle',False)

    def fold(self, manifest):
        path=Path(str(manifest))
        if not path.is_relative_to('/work'):
            raise ValueError('Oracle manifest must be inside this job.')
        local=(self.run_dir/path.relative_to('/work')).resolve()
        local.relative_to(self.run_dir)
        if not local.is_file(): raise ValueError('Oracle manifest is missing.')
        entries=json.loads(local.read_text())
        if not isinstance(entries,list) or not 1<=len(entries)<=10000:
            raise ValueError('Invalid oracle manifest batch.')
        for entry in entries:
            for field in ('output_dir','out_dir'):
                if entry.get(field) and not Path(entry[field]).is_relative_to('/work'):
                    raise ValueError('Oracle outputs must remain inside this job.')
        if self.engine!='openfold3_jax':
            from .rfd4_pipeline_oracle import fold_entry
            for entry in entries: fold_entry(entry,self)
            return
        cache=Path(reference_root())/'rfd4-proteina/jax-cache'; cache.mkdir(parents=True,exist_ok=True)
        root=self.run_dir/'artifacts/rfd4_proteina'
        weights=model_bundle()/'openfold3/weights'
        command=['docker','run','--rm',*docker_gpu_args(self.gpu_device),'--cpus',str(self.cpu_cores),
                 '--shm-size=8g','-v',f'{self.run_dir}:/work','-v',f'{Path(reference_root())/"rfd4-proteina/huggingface"}:/cache/huggingface:ro',
                 '-v',f'{cache}:/cache','-e','XLA_PYTHON_CLIENT_PREALLOCATE=false','-e','JAX_PLATFORMS=cuda',
                 '-e','XLA_FLAGS=--xla_gpu_enable_triton_gemm=false',
                 '--entrypoint','python',ORACLE_IMAGE,
                 '/work/artifacts/rfd4_proteina/oracle_runtime/run_alphafold_batch.py',
                 f'--manifest={path}',f'--model_dir=/cache/huggingface/hub/models--nvidia--RFD4-Proteina/snapshots/{HF_REVISION}/openfold3/weights',f'--num_recycles={int(self.options.get("folding_recycles",5))}','--num_diffusion_samples=1',
                 '--flash_attention_implementation=xla','--jax_compilation_cache_dir=/cache','--save_full_outputs=false']
        with (root/'oracle_commands.jsonl').open('a') as out:
            out.write(json.dumps(command)+'\n')
        with (root/'oracle.log').open('a') as out:
            self.process=subprocess.Popen(command,stdout=out,stderr=subprocess.STDOUT)
            try:
                code=self.process.wait(timeout=max(3600,len(entries)*300))
                if code: raise RuntimeError(f'OpenFold3 oracle exited with status {code}; see artifacts/rfd4_proteina/oracle.log')
            except BaseException:
                self.process.terminate()
                raise
            finally: self.process=None

    def __exit__(self,*exc):
        if self.process is not None: self.process.terminate()
        if self.server is not None:
            self.server.shutdown(); self.server.server_close()
        if self.thread is not None: self.thread.join(timeout=5)
        if self.temporary is not None: self.temporary.cleanup()
