# RFpeptides

Studio's RFpeptides tile wraps the macrocyclic binder-backbone protocol in
[RosettaCommons/RFdiffusion](https://github.com/RosettaCommons/RFdiffusion),
pinned at `86507b6538f51fce57b5a72477165f03999ed7ae`.

The standalone job accepts a PDB/mmCIF target, selected protein chains, optional
hotspots, a fixed or ranged peptide length, a design count, and RFdiffusion's
denoising step count. It writes cyclic peptide backbone candidates and a Mol*
complex viewer. Sequence design is a separate downstream step.

## Build the CUDA 12.8 image

The RFdiffusion source checkout is kept at
`/home/user/programs/git-projects/tools-to-implement/RFdiffusion`. Build its
Blackwell-ready inference image from the checkout root:

```bash
docker build -f docker/Dockerfile.cuda128 -t mn-rfpeptides:86507b6-cu128 .
```

The image uses PyTorch 2.7.1 with CUDA 12.8 and checks for both `sm_89` (RTX
4090) and `sm_120` (RTX 5090) kernels during the build. It compiles the upstream
DGL 1.1.3 source against CUDA 12.8, with native kernels for both GPU generations.

## Install the model checkpoint

Model weights are stored outside the image in Studio's reference directory. The
standalone tool requires `rfdiffusion_models/Complex_base_ckpt.pt`; download it
from the URL published by the upstream RFdiffusion repository:

```bash
REF_DIR=/path/to/Studio/reference_files/rfdiffusion_models
mkdir -p "$REF_DIR"
wget -P "$REF_DIR" \
  https://files.ipd.uw.edu/pub/RFdiffusion/e29311f6f1bf1af907f9ef9f44b8328b/Complex_base_ckpt.pt
```

Replace `/path/to/Studio/reference_files` with the active `reference_root()`
directory configured for the Studio installation.

The scheduler reports the image and checkpoint as required resources. Each GPU
workstation must have this image tag and the shared reference checkpoint before
it can run RFpeptides jobs.
