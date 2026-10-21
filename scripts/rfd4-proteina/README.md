# RFD4-Proteina image

Build the CUDA 13 image used by the RFD4-Proteina tile:

```bash
bash scripts/rfd4-proteina/build_image.sh
```

This is one image for RTX 4090 and RTX 5090 workstations. The image is pinned to upstream RFD4-Proteina commit `9d35ee4f6328fd9b17f5b84756ee5db99d6d89ac` (package 0.0.3). Upstream documents a CUDA 13-compatible driver requirement, but its verified-hardware table does not list these two cards; the card compatibility is expected from the cu130 runtime, not upstream-certified.

The first run downloads about 1.7 GB of public model weights to the reference cache. Studio mounts that cache at `/cache/huggingface`, so later runs on that workstation reuse it. The cache is local to that workstation and is not copied by the app's structure/database sync.

Smoke-check the image without starting a design:

```bash
docker run --rm --gpus all mn-rfd4-proteina:0.0.3-cu130 --help
```
