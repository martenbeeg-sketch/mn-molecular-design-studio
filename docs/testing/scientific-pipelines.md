# Scientific pipeline regression tests

Use the existing `mn-molecular-design-studio` Conda Python. CPU pytest never launches GPU tests unless `MDS_GPU_TESTS=1` is set. Mark real single-engine inference with `gpu`, and real multi-stage inference with `gpu_e2e`. Neither marker authorizes downloads or large campaigns.

```bash
python -m pytest tests/test_scientific_contracts.py -q
MDS_GPU_TESTS=1 MDS_TEST_GPU=0 python -m pytest tests/test_scientific_gpu.py -v -s
```

`MDS_GPU_ARTIFACTS` selects an artifact directory (default `.tmp/gpu-tests`). Each test retains runtime/device/image metadata and production job logs, effective parameters, normalized candidates and structures. `MDS_TEST_GPU` selects a device from `nvidia-smi`. Tests use installed reference caches and images. A missing prerequisite is SKIP; an inference or scientific-contract error after launch is FAIL. Preparation-only checks are not GPU inference. No comparison requires identical coordinates or numerical scores across devices.

These checks establish computational integration and scientific data handling, not experimental predictive validity.
