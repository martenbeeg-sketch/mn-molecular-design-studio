# RFD4-Proteina 0.0.3 CUDA 13 image.
# One cu130 image is used on RTX 4090 (Ada) and RTX 5090 (Blackwell) hosts;
# the host needs a CUDA-13-compatible NVIDIA driver and NVIDIA Container Toolkit.
# The Docker build downloads the upstream repository at the pinned commit below.
FROM python:3.12.12-slim-bookworm@sha256:593bd06efe90efa80dc4eee3948be7c0fde4134606dd40d8dd8dbcade98e669c

ARG RFD4_COMMIT=9d35ee4f6328fd9b17f5b84756ee5db99d6d89ac
ARG UV_VERSION=0.12.13

ENV DATA_PATH=/opt/rfd4/data \
    HF_HOME=/cache/huggingface \
    HF_HUB_CACHE=/cache/huggingface/hub \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PYTHONNOUSERSITE=1 \
    TMOL_DISABLE_WHEEL_FETCH=1 \
    TMOL_JIT_FALLBACK=0 \
    TMOL_USE_JIT=0 \
    NVIDIA_VISIBLE_DEVICES=all \
    NVIDIA_DRIVER_CAPABILITIES=compute,utility

RUN apt-get update && \
    apt-get install -y --no-install-recommends build-essential ca-certificates git libxrender1 libx11-6 libxext6 && \
    rm -rf /var/lib/apt/lists/* && \
    python -m pip install --no-cache-dir "uv==${UV_VERSION}"

WORKDIR /opt/rfd4
RUN git init . && \
    git remote add origin https://github.com/RosettaCommons/RFD4-Proteina.git && \
    git fetch --depth 1 origin "${RFD4_COMMIT}" && \
    git checkout --detach FETCH_HEAD

# Install the exact upstream CUDA-13 extra, including Torch cu130 and the fused kernels.
RUN uv pip install --system --no-cache --no-sources \
      --only-binary tmol,tmol-ops-cu12,tmol-ops-cu13 \
      --extra all-cuda13 --torch-backend cu130 -r pyproject.toml && \
    uv pip install --system --no-cache --no-sources --no-deps . && \
    python -m pip check

RUN python -c "import torch, atomworks, rfproteina; from openbabel import pybel; from tmol.score.ljlk.potentials._compiled import ljlk_pose_scores; assert pybel.readstring('smi', 'CCO').OBMol.NumAtoms() == 3; assert callable(ljlk_pose_scores)" && \
    rfd4-generate --help >/dev/null && \
    mkdir -p /cache/huggingface /workspace

WORKDIR /workspace
ENTRYPOINT ["rfd4-generate"]
CMD ["--help"]
