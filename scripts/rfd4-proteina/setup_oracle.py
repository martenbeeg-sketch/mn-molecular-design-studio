#!/usr/bin/env python3
"""Import NVIDIA's checksum-pinned RFD4 OpenFold3 runtime into local Docker.

Requires the base Studio RFD4 image, Docker, and unsquashfs. No host Python
packages are installed. Run with --reference-dir matching Studio Settings.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

REVISION='598acab93d5d0aeacc4362f6e88dc20e7db85ef6'
SHA256='7bbd35478ed5e95131efa7756d03a2e754f3508e603ad08bb8f8db0bd89484a1'
FILE='openfold3/runtime/of3_af3_refolding_env-cuda13.sif'
IMAGE='mn-rfd4-oracle:598acab-cu13'
ENVIRONMENT={
    'PATH':'/hmmer/bin:/alphafold3_venv/bin:/usr/local/nvidia/bin:/usr/local/cuda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin',
    'LD_LIBRARY_PATH':'/usr/local/nvidia/lib:/usr/local/nvidia/lib64:/usr/local/cuda/lib64',
    'NVIDIA_VISIBLE_DEVICES':'all','NVIDIA_DRIVER_CAPABILITIES':'compute,utility',
    'XLA_FLAGS':'--xla_gpu_enable_triton_gemm=false',
}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference-dir',required=True,type=Path)
    parser.add_argument('--keep-rootfs',action='store_true')
    args=parser.parse_args()
    cache=args.reference_dir.resolve()/'rfd4-proteina/huggingface';cache.mkdir(parents=True,exist_ok=True)
    sif=cache/'hub/models--nvidia--RFD4-Proteina/snapshots'/REVISION/FILE
    if not sif.is_file():
        code=f'from huggingface_hub import hf_hub_download; hf_hub_download("nvidia/RFD4-Proteina",{FILE!r},revision={REVISION!r})'
        subprocess.run(['docker','run','--rm','-v',f'{cache}:/cache/huggingface','--entrypoint','python',
                        'mn-rfd4-proteina:0.0.3-cu130','-c',code],check=True)
    with sif.open('rb') as handle: digest=hashlib.file_digest(handle,'sha256').hexdigest()
    if digest!=SHA256: raise ValueError('Runtime does not match the published SHA256; refusing import.')
    with sif.open('rb') as handle: offset=handle.read(4*1024*1024).find(b'hsqs')
    if offset<0: raise ValueError('SquashFS payload not found.')
    build=args.reference_dir.resolve()/'rfd4-proteina/oracle-runtime-598acab';root=build/'rootfs'
    build.mkdir(exist_ok=True)
    if not root.exists():
        subprocess.run(['unsquashfs','-no-xattrs','-processors','4','-o',str(offset),'-d',str(root),str(sif)],check=True)
    command=['docker','import']
    for key,value in ENVIRONMENT.items():command.extend(['--change',f'ENV {key}={value}'])
    for key,value in {'source':'https://huggingface.co/nvidia/RFD4-Proteina','revision':REVISION}.items():
        command.extend(['--change',f'LABEL org.opencontainers.image.{key}={value}'])
    command.extend(['--message',f'Official NVIDIA RFD4 OF3 runtime; SHA256 {digest}','-',IMAGE])
    with subprocess.Popen(['tar','--numeric-owner','-C',str(root),'-cf','-','.'],stdout=subprocess.PIPE) as tar:
        imported=subprocess.run(command,stdin=tar.stdout,capture_output=True,text=True)
        tar.stdout.close();tar_code=tar.wait()
        if tar_code or imported.returncode: raise RuntimeError(imported.stderr or 'Root filesystem export failed.')
    metadata={'source':'https://huggingface.co/nvidia/RFD4-Proteina','revision':REVISION,'file':FILE,
              'sha256':digest,'image':IMAGE,'image_id':imported.stdout.strip(),'environment':ENVIRONMENT}
    (build/'provenance.json').write_text(json.dumps(metadata,indent=2))
    if not args.keep_rootfs:shutil.rmtree(root)
    print(json.dumps(metadata,indent=2))


if __name__=='__main__':main()
