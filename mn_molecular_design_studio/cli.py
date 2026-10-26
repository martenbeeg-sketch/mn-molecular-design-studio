from __future__ import annotations

import typer
import uvicorn

from mn_molecular_design_studio.auth import create_initial_admin, reset_admin_password
from mn_molecular_design_studio.workspace import ensure_studio_workspace


app = typer.Typer(
    help="Launch the local Molecular Design Studio API.",
    no_args_is_help=True,
    context_settings={"help_option_names": ["-h", "--help"]},
)
from mn_molecular_design_studio.cluster.cli import app as cluster_commands
app.add_typer(cluster_commands, name="cluster")


@app.callback()
def cli() -> None:
    """Local application commands for Molecular Design Studio."""


@app.command()
def api(port: int = typer.Option(8766, min=1, max=65535, help="Local API port."),
        host: str = typer.Option("127.0.0.1", help="Listen address; use SSH tunnels or authenticated private-network access for workers.")) -> None:
    """Start the Studio API and optional coordinator endpoints."""
    from mn_molecular_design_studio.auth import is_development_mode
    if host not in {"127.0.0.1", "localhost", "::1"} and is_development_mode():
        raise typer.BadParameter("Configure production sign-in before exposing the API beyond this computer. SSH tunnels can use the loopback address.")
    uvicorn.run(
        "mn_molecular_design_studio.api:app",
        host=host,
        port=port,
        reload=False,
        access_log=False,
    )


@app.command("create-admin")
def create_admin(
    username: str = typer.Option(..., prompt=True, help="Admin username or email address."),
    display_name: str = typer.Option("", prompt="Display name", help="Name shown in the Studio."),
) -> None:
    """Create the first Studio admin account on this computer."""
    ensure_studio_workspace()
    password = typer.prompt("Admin password", hide_input=True, confirmation_prompt=True)
    try:
        user = create_initial_admin(username, display_name, password)
    except (OSError, TypeError, ValueError) as exc:
        raise typer.BadParameter(str(exc)) from exc
    typer.echo(f"Admin account created for {user['username']}. You can now sign in to Molecular Design Studio.")


@app.command("mcp")
def mcp_server(
    transport: str = typer.Option("stdio", help="stdio or streamable-http."),
    host: str = typer.Option("127.0.0.1", help="HTTP listen address."),
    port: int = typer.Option(8767, min=1, max=65535, help="HTTP MCP port."),
) -> None:
    """Start the optional MCP adapter, using the existing Studio REST API."""
    try:
        from mn_molecular_design_studio.mcp_server import run
        run(transport, host, port)
    except ImportError as exc:
        raise typer.BadParameter("Install the Studio mcp extra in this environment first.") from exc


@app.command("reset-admin-password")
def reset_admin(
    username: str = typer.Option(..., prompt=True, help="Existing admin username or email address."),
) -> None:
    """Recover the admin password from the Studio host computer."""
    ensure_studio_workspace()
    password = typer.prompt("New admin password", hide_input=True, confirmation_prompt=True)
    try:
        reset_admin_password(username, password)
    except (OSError, TypeError, ValueError) as exc:
        raise typer.BadParameter(str(exc)) from exc
    typer.echo(f"Admin password reset for {username.strip()}. All previous sessions were signed out.")


def main() -> None:
    app()


if __name__ == "__main__":
    main()
