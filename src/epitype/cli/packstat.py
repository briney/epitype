"""PackStat calculation subcommand."""
from __future__ import annotations

import json
from pathlib import Path

import typer
from rich.console import Console
from rich.table import Table

from epitype.cli.common import StructureResolutionError, resolve_structure
from epitype.core.interface import detect_interface, parse_chain_groups
from epitype.io.parsers import parse_structure
from epitype.metrics.packstat import calculate_packstat

console = Console()


def packstat_command(
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
    Calculate packing quality (PackStat) for a protein-protein interface.

    PackStat uses a multi-probe algorithm and measures how well
    the interface is packed.

    Examples:
        epitype packstat structure.pdb --chains HL_A
        epitype packstat 1yy9 --chains CD_A
    """
    # Resolve structure input (file path or PDB ID)
    try:
        structure_path = resolve_structure(structure)
    except StructureResolutionError as e:
        console.print(f"[red]Error: {e}[/red]")
        raise typer.Exit(1) from e

    # Parse structure
    struct = parse_structure(structure_path)

    # Parse chain specification
    group1, group2 = parse_chain_groups(chains)

    # Detect interface
    interface = detect_interface(struct, group1, group2, cutoff=cutoff)

    if interface.num_interface_residues == 0:
        console.print("[red]Error: No interface detected between specified chains.[/red]")
        raise typer.Exit(1)

    # Calculate packstat
    console.print("Calculating PackStat...")
    result = calculate_packstat(struct, interface)

    if output:
        data = {
            "packstat": result.packstat,
            "packstat_interface": result.packstat_interface,
            "cavity_volume": result.cavity_volume,
        }

        with open(output, "w") as f:
            json.dump(data, f, indent=2)
        console.print(f"Results written to {output}")
    else:
        table = Table(title=f"PackStat: {structure_path.name}")
        table.add_column("Metric", style="cyan")
        table.add_column("Value", style="green", justify="right")

        # Determine quality
        if result.packstat > 0.65:
            quality = "[green]good[/green]"
        elif result.packstat > 0.5:
            quality = "[yellow]moderate[/yellow]"
        else:
            quality = "[red]poor[/red]"

        table.add_row("PackStat (overall)", f"{result.packstat:.3f} ({quality})")
        table.add_row("PackStat (interface)", f"{result.packstat_interface:.3f}")
        table.add_row("Cavity volume (A^3)", f"{result.cavity_volume:.1f}")

        console.print(table)
