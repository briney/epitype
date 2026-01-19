"""Shape complementarity subcommand."""
from __future__ import annotations

import json
from pathlib import Path

import typer
from rich.console import Console
from rich.table import Table

from epitype.cli.common import StructureResolutionError, resolve_structure

console = Console()


def shape_command(
    structure: str = typer.Argument(
        ...,
        help="Path to PDB/mmCIF file or PDB ID (e.g., 1yy9)",
    ),
    chains: str = typer.Option(
        ...,
        "--chains",
        "-c",
        help="Chain specification (e.g., 'HL_A')",
    ),
    cutoff: float = typer.Option(
        8.0,
        "--cutoff",
        help="Interface detection cutoff in Angstroms",
    ),
    output: Path | None = typer.Option(
        None,
        "--output",
        "-o",
        help="Output file (JSON)",
    ),
):
    """
    Calculate shape complementarity (Sc) for a protein-protein interface.

    Requires NanoShaper to be installed and available in PATH.

    Examples:
        epitype shape structure.pdb --chains HL_A
        epitype shape 1yy9 --chains CD_A
    """
    from epitype.surface.nanoshaper import check_nanoshaper_available

    if not check_nanoshaper_available():
        console.print("[red]Error: NanoShaper is not available.[/red]")
        console.print("Download from: https://gitlab.iit.it/SDecherchi/nanoshaper")
        raise typer.Exit(1)

    # Resolve structure input (file path or PDB ID)
    try:
        structure_path = resolve_structure(structure)
    except StructureResolutionError as e:
        console.print(f"[red]Error: {e}[/red]")
        raise typer.Exit(1) from e

    from epitype.core.interface import detect_interface, parse_chain_groups
    from epitype.io.parsers import parse_structure
    from epitype.metrics.shape import calculate_shape_complementarity

    # Parse structure
    struct = parse_structure(structure_path)

    # Parse chain specification
    group1, group2 = parse_chain_groups(chains)

    # Detect interface
    interface = detect_interface(struct, group1, group2, cutoff=cutoff)

    if interface.num_interface_residues == 0:
        console.print("[red]Error: No interface detected between specified chains.[/red]")
        raise typer.Exit(1)

    # Calculate shape complementarity
    console.print("Calculating shape complementarity (this may take a moment)...")
    result = calculate_shape_complementarity(struct, interface)

    if output:
        data = {
            "sc_value": result.sc_value,
            "sc_median": result.sc_median,
            "sc_group1": result.sc_group1,
            "sc_group2": result.sc_group2,
            "interface_area": result.interface_area,
        }

        with open(output, "w") as f:
            json.dump(data, f, indent=2)
        console.print(f"Results written to {output}")
    else:
        table = Table(title=f"Shape Complementarity: {structure_path.name}")
        table.add_column("Metric", style="cyan")
        table.add_column("Value", style="green", justify="right")

        # Determine quality
        if result.sc_value > 0.65:
            quality = "[green]good[/green]"
        elif result.sc_value > 0.5:
            quality = "[yellow]moderate[/yellow]"
        else:
            quality = "[red]poor[/red]"

        table.add_row("Sc (overall)", f"{result.sc_value:.3f} ({quality})")
        table.add_row("Sc (median)", f"{result.sc_median:.3f}")
        table.add_row("Sc (group 1)", f"{result.sc_group1:.3f}")
        table.add_row("Sc (group 2)", f"{result.sc_group2:.3f}")
        table.add_row("Interface area", f"{result.interface_area:.1f}")

        console.print(table)
