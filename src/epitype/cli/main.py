"""Main CLI entry point using Typer."""
from __future__ import annotations

import typer
from rich.console import Console

from epitype.cli.energy import energy_command
from epitype.cli.hbonds import hbonds_command
from epitype.cli.packstat import packstat_command
from epitype.cli.run import run_command
from epitype.cli.sasa import sasa_command
from epitype.cli.shape import shape_command
from epitype.cli.smoke import smoke_command

app = typer.Typer(
    name="epitype",
    help="Protein-protein interface analysis toolkit.",
    no_args_is_help=True,
)

console = Console()

# Register commands
app.command(name="run")(run_command)
app.command(name="sasa")(sasa_command)
app.command(name="energy")(energy_command)
app.command(name="hbonds")(hbonds_command)
app.command(name="shape")(shape_command)
app.command(name="packstat")(packstat_command)
app.command(name="smoke")(smoke_command)


@app.command()
def version():
    """Show version information."""
    from epitype import __version__

    console.print(f"epitype v{__version__}")


@app.command()
def info():
    """Show system information and external tool status."""
    from rich.table import Table

    from epitype import __version__
    from epitype.bin import get_nanoshaper_info

    console.print(f"[bold]epitype[/bold] v{__version__}")
    console.print()

    table = Table(title="External Tools")
    table.add_column("Tool", style="cyan")
    table.add_column("Status")
    table.add_column("Path")
    table.add_column("Source")

    # NanoShaper info
    ns_info = get_nanoshaper_info()
    ns_status = "[green]available[/green]" if ns_info["available"] else "[red]not found[/red]"
    table.add_row("NanoShaper", ns_status, ns_info["path"], ns_info["source"])

    console.print(table)
    console.print()
    console.print(f"[dim]Platform: {ns_info['platform']}[/dim]")
    console.print("[dim]Set EPITYPE_NANOSHAPER_PATH to use a custom NanoShaper binary[/dim]")


@app.callback()
def main_callback():
    """Interface Analyzer: Compute protein-protein interface metrics."""
    pass


def main():
    """Main entry point for the CLI."""
    app()


if __name__ == "__main__":
    main()
