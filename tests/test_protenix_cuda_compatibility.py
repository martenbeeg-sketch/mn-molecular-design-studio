"""Regression for stock Protenix admitting Blackwell with an Ada-only extension."""
from mn_molecular_design_studio.cluster import gpu_runtime


def test_torch_blackwell_support_does_not_override_missing_protenix_kernel():
    profile = {'cuda': '12.8', 'targets': ['sm_89', 'sm_120', 'compute_120'],
               'required_cuda_targets': {'protenix_fused_layer_norm': ['sm_89', 'sm_100']}}
    assert gpu_runtime.compatible({'compute_capability': '8.9', 'driver_cuda_api': 13000}, {'image': profile})
    assert not gpu_runtime.compatible({'compute_capability': '12.0', 'driver_cuda_api': 13000}, {'image': profile})
    profile['required_cuda_targets']['protenix_fused_layer_norm'].append('sm_120')
    assert gpu_runtime.compatible({'compute_capability': '12.0', 'driver_cuda_api': 13000}, {'image': profile})


def test_missing_required_extension_inventory_fails_closed():
    profile = {'cuda': '12.8', 'targets': ['sm_120'],
               'required_cuda_targets': {'protenix_fused_layer_norm': []}}
    assert not gpu_runtime.compatible({'compute_capability': '12.0', 'driver_cuda_api': 13000}, {'image': profile})


def test_required_extension_ptx_forward_compatibility_retains_driver_check():
    profile = {'cuda': '12.8', 'targets': ['sm_120'],
               'required_cuda_targets': {'protenix_fused_layer_norm': ['compute_89']}}
    assert gpu_runtime.compatible({'compute_capability': '12.0', 'driver_cuda_api': 13000}, {'image': profile})
    assert not gpu_runtime.compatible({'compute_capability': '12.0', 'driver_cuda_api': 12070}, {'image': profile})


def test_failed_protenix_probe_is_retained_as_incompatible(monkeypatch):
    import subprocess
    def inspect(*args, **kwargs):
        return subprocess.CompletedProcess(args[0], 0, stdout='sha256:installed\n')
    monkeypatch.setattr(gpu_runtime.subprocess, 'run', inspect)
    monkeypatch.setattr(gpu_runtime, '_probe', lambda *args: {'cuda': None, 'targets': []})
    result = gpu_runtime.profiles(['mn-protenix-opt:broken'])
    assert not gpu_runtime.compatible({'compute_capability': '12.0', 'driver_cuda_api': 13000}, result)


def test_architecture_specific_ptx_is_not_treated_as_portable(monkeypatch):
    def dump(command, **kwargs):
        if '--list-elf' in command:
            return 'ELF file 1: kernel.sm_89.cubin\nELF file 2: kernel.sm_100a.cubin'
        return '.target sm_89\n.target sm_100a\n.target sm_120f\n'
    monkeypatch.setattr(gpu_runtime.subprocess, 'check_output', dump)
    scope = {'inventory': {}}
    exec(gpu_runtime._PROTENIX_KERNEL_PROBE, scope)
    assert scope['inventory']['required_cuda_targets']['protenix_fused_layer_norm'] == ['sm_89', 'compute_89']
