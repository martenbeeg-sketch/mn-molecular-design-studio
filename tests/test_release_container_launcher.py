from pathlib import Path

import pytest

from scripts.test_release_container import build_command


def command(**changes):
    options = dict(source=Path("/candidate"), runtime_root=Path("/runtime"), python=Path("/runtime/bin/python"),
                   image="sha256:" + "a" * 64, output=Path("/private/evidence"),
                   arguments=["-q", "tests/test_post_phase2_release.py"], user="1000:1000")
    return build_command(**{**options, **changes})


def test_container_tests_have_no_network_host_data_or_daemon_access():
    arguments = command()
    assert {"--read-only", "--network=none", "--cap-drop=ALL", "--pull=never"} <= set(arguments)
    mounts = [arguments[index + 1] for index, value in enumerate(arguments) if value == "--mount"]
    assert mounts == ["type=bind,src=/candidate,dst=/candidate,readonly",
                      "type=bind,src=/runtime,dst=/runtime,readonly",
                      "type=bind,src=/private/evidence,dst=/evidence"]
    assert "MN_MOLECULAR_DESIGN_STUDIO_AUTH_MODE=production" in arguments
    assert "MDS_GPU_TESTS=0" in arguments
    assert arguments[-3:] == ["--junitxml=/evidence/junit.xml", "-q", "tests/test_post_phase2_release.py"]


@pytest.mark.parametrize("changes", [
    {"image": "mutable:latest"}, {"runtime_root": Path("/")},
    {"python": Path("/outside/python")}, {"output": Path("/candidate/results")},
    {"output": Path("/")},
])
def test_container_test_configuration_fails_closed(changes):
    with pytest.raises(ValueError):
        command(**changes)
