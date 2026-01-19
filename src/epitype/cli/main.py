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


@app.command()
def version():
    """Show version information."""
    from epitype import __version__

    console.print(f"epitype v{__version__}")


@app.callback()
def main_callback():
    """Interface Analyzer: Compute protein-protein interface metrics."""
    pass


def main():
    """Main entry point for the CLI."""
    app()


if __name__ == "__main__":
    main()
