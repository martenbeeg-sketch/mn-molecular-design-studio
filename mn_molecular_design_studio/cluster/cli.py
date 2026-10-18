from __future__ import annotations

import json
import os
import shlex
import subprocess
import sys
from pathlib import Path

import typer

from . import config
from .store import Store

app = typer.Typer(help="Connect Studio workstations while keeping results on their execution machine.", no_args_is_help=True)


@app.command("identity")
def identity():
    typer.echo(config.identity())


@app.command("init")
def initialize(name: str = typer.Option(...), coordinator_url: str = "http://127.0.0.1:8766"):
    """Make this installation the coordinator and enroll its local agent."""
    settings = config.load()
    if settings["mode"] == "worker":
        raise typer.BadParameter("Disable the worker connection before creating a new coordinator.")
    if settings["mode"] == "coordinator" and settings["node_token"]:
        config.save({"name": name, "coordinator_url": coordinator_url})
    else:
        enrollment = Store().enroll(name, config.identity())
        config.save({"mode": "coordinator", "name": name, "coordinator_url": coordinator_url,
                     "node_token": enrollment["node_token"], "default_placement": "auto"})
    typer.echo(f"Coordinator configured: {name} ({config.identity()}). Start the Studio API and cluster agent.")


@app.command("enroll")
def enroll(name: str = typer.Option(...), output: Path = typer.Option(...), coordinator_url: str = ""):
    """Create a private enrollment file for one additional workstation."""
    if config.load()["mode"] != "coordinator":
        raise typer.BadParameter("Run this command on the coordinator.")
    if output.exists():
        raise typer.BadParameter("Choose a new enrollment filename; existing credentials are not overwritten.")
    payload = Store().enroll(name) | {"coordinator_url": coordinator_url or config.load()["coordinator_url"]}
    output.parent.mkdir(parents=True, exist_ok=True)
    config.atomic_json(output, payload)
    typer.echo(f"Enrollment saved to {output.resolve()}. Copy it privately to {name} and run cluster join.")


@app.command("join")
def join(enrollment: Path = typer.Argument(..., exists=True), coordinator_url: str = ""):
    """Enroll this installation with the shared coordinator."""
    import uuid
    payload = json.loads(enrollment.read_text())
    node_id = str(uuid.UUID(payload["node_id"]))
    current = config.load()
    if current["mode"] != "disabled" and config.identity() != node_id:
        raise typer.BadParameter("This installation already belongs to another workstation identity.")
    if len(payload.get("node_token", "")) < 32:
        raise typer.BadParameter("The enrollment file has no valid workstation credential.")
    config.save({"mode": "worker", "name": payload["name"], "coordinator_url": coordinator_url or payload["coordinator_url"],
                 "node_token": payload["node_token"], "default_placement": "auto"})
    path = config.home() / "node-id"
    path.write_text(node_id + "\n")
    path.chmod(0o600)
    typer.echo(f"Connected as {payload['name']}. Start cluster agent, then sign in with your coordinator account.")


@app.command("agent")
def agent(once: bool = False):
    """Publish local resources/results and accept assigned workstation commands."""
    import logging
    from .agent import Agent
    from mn_molecular_design_studio.workspace import ensure_studio_workspace
    ensure_studio_workspace()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    Agent().run(once=once)


@app.command("status")
def status():
    settings = config.load()
    result = {k: settings[k] for k in ("mode", "name", "coordinator_url", "accept_jobs", "default_placement")}
    result["node_id"] = config.identity()
    if settings["mode"] == "coordinator":
        result["nodes"] = Store().nodes()
    typer.echo(json.dumps(result, indent=2))


@app.command("disable")
def disable():
    """Stop cluster participation; already running local jobs continue."""
    config.save({"mode": "disabled"})
    typer.echo("Workstation sharing disabled. Local scientific jobs retain their existing scheduler ownership.")


@app.command("install-services")
def install_services(enable: bool = False, api_port: int = 8766, api_host: str = "127.0.0.1"):
    """Install user-level systemd units for the API and workstation agent."""
    from mn_molecular_design_studio.workspace import _studio_home
    if not config.enabled():
        raise typer.BadParameter("Configure this workstation with cluster init or cluster join first.")
    if not 1 <= api_port <= 65535:
        raise typer.BadParameter("Invalid API port.")
    directory = Path.home() / ".config/systemd/user"
    directory.mkdir(parents=True, exist_ok=True)
    executable = shlex.quote(sys.executable)
    environment = "\n".join("Environment=" + shlex.quote(key + "=" + str(value)) for key, value in {
        "MN_MOLECULAR_DESIGN_STUDIO_HOME": _studio_home(), "MDS_CLUSTER_HOME": config.home(),
    }.items())
    for name, arguments in (("studio-api", f"api --port {api_port} --host {shlex.quote(api_host)}"), ("studio-agent", "cluster agent")):
        unit = ("[Unit]\nDescription=MN Molecular Design Studio " + name + "\nAfter=network-online.target\n\n[Service]\nType=simple\n"
                + environment + "\nExecStart=" + executable + " -m mn_molecular_design_studio.cli " + arguments
                + "\nRestart=on-failure\nRestartSec=10\n\n[Install]\nWantedBy=default.target\n")
        (directory / (name + ".service")).write_text(unit)
    if enable:
        subprocess.run(["systemctl", "--user", "daemon-reload"], check=True)
        subprocess.run(["systemctl", "--user", "enable", "--now", "studio-api.service", "studio-agent.service"], check=True)
    typer.echo("User services installed. Enable with: systemctl --user enable --now studio-api.service studio-agent.service")


@app.command("install-tunnel")
def install_tunnel(ssh_host: str = typer.Option(..., help="Existing SSH destination, such as user@coordinator-host."),
                   local_port: int = 18766, remote_port: int = 8766, enable: bool = False):
    """Keep an SSH connection from a worker to the coordinator's loopback API."""
    if not ssh_host or ssh_host.startswith("-") or any(c.isspace() for c in ssh_host):
        raise typer.BadParameter("Enter one SSH hostname or user@hostname destination.")
    if not 1 <= local_port <= 65535 or not 1 <= remote_port <= 65535:
        raise typer.BadParameter("Invalid tunnel port.")
    directory = Path.home() / ".config/systemd/user"
    directory.mkdir(parents=True, exist_ok=True)
    arguments = ["/usr/bin/ssh", "-N", "-o", "BatchMode=yes", "-o", "ExitOnForwardFailure=yes", "-o", "ServerAliveInterval=30",
                 "-o", "ServerAliveCountMax=3", "-L", f"127.0.0.1:{local_port}:127.0.0.1:{remote_port}", ssh_host]
    unit = ("[Unit]\nDescription=Studio coordinator SSH tunnel\nAfter=network-online.target\n\n[Service]\nType=simple\nExecStart="
            + shlex.join(arguments) + "\nRestart=on-failure\nRestartSec=10\n\n[Install]\nWantedBy=default.target\n")
    (directory / "studio-coordinator-tunnel.service").write_text(unit)
    if enable:
        subprocess.run(["systemctl", "--user", "daemon-reload"], check=True)
        subprocess.run(["systemctl", "--user", "enable", "--now", "studio-coordinator-tunnel.service"], check=True)
    typer.echo(f"Tunnel service installed. Join using coordinator URL http://127.0.0.1:{local_port}")
