"""Main analysis run command."""
from __future__ import annotations

from pathlib import Path

import typer
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn

from epitype.analysis.pipeline import AnalysisConfig, analyze_interface
from epitype.io.writers import format_metrics_table, write_csv, write_json

console = Console()


def run_command(
    structure: Path = typer.Argument(
        ...,
        help="Path to PDB or mmCIF file",
        exists=True,
        readable=True,
    ),
    chains: str = typer.Option(
        ...,
        "--chains",
        "-c",
        help="Chain specification (e.g., 'HL_A' for chains H+L vs A)",
    ),
    output: Path | None = typer.Option(
        None,
        "--output",
        "-o",
        help="Output file path (JSON or CSV based on extension)",
    ),
    cutoff: float = typer.Option(
        8.0,
        "--cutoff",
        help="Interface detection cutoff in Angstroms",
    ),
    no_energy: bool = typer.Option(
        False,
        "--no-energy",
        help="Skip energy calculation (faster)",
    ),
    no_shape: bool = typer.Option(
        False,
        "--no-shape",
        help="Skip shape complementarity (faster, no NanoShaper required)",
    ),
    no_packstat: bool = typer.Option(
        False,
        "--no-packstat",
        help="Skip PackStat calculation (faster)",
    ),
    quiet: bool = typer.Option(
        False,
        "--quiet",
        "-q",
        help="Suppress progress output",
    ),
):
    """
    Analyze protein-protein interface.

    Computes binding energy, buried surface area, shape complementarity,
    packing quality, and hydrogen bond metrics.

    Example:
        epitype run structure.pdb --chains HL_A -o results.json
    """
    config = AnalysisConfig(
        chain_spec=chains,
        interface_cutoff=cutoff,
        compute_energy=not no_energy,
        compute_shape=not no_shape,
        compute_packstat=not no_packstat,
    )

    try:
        if not quiet:
            with Progress(
                SpinnerColumn(),
                TextColumn("[progress.description]{task.description}"),
                console=console,
            ) as progress:
                task = progress.add_task("Analyzing...", total=None)

                def progress_callback(step: str, frac: float):
                    progress.update(task, description=f"{step}...")

                metrics = analyze_interface(structure, config, progress_callback)
        else:
            metrics = analyze_interface(structure, config)
    except ValueError as e:
        console.print(f"[red]Error: {e}[/red]")
        raise typer.Exit(1) from e
    except FileNotFoundError as e:
        console.print(f"[red]Error: {e}[/red]")
        raise typer.Exit(1) from e

    # Output results
    if output:
        if output.suffix == ".json":
            write_json(metrics, output, structure_name=structure.stem)
            console.print(f"Results written to {output}")
        elif output.suffix == ".csv":
            write_csv(metrics, output, structure_name=structure.stem)
            console.print(f"Results written to {output}")
        else:
            # Default to JSON
            write_json(metrics, output, structure_name=structure.stem)
            console.print(f"Results written to {output}")
    else:
        # Print to console
        console.print()
        console.print(format_metrics_table(metrics, structure_name=structure.name))
