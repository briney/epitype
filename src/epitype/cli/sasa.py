"""SASA calculation subcommand."""
from __future__ import annotations

import json
from pathlib import Path

import typer
from rich.console import Console
from rich.table import Table

from epitype.cli.common import StructureResolutionError, resolve_structure
from epitype.io.parsers import parse_structure
from epitype.metrics.sasa import calculate_sasa

console = Console()


def sasa_command(
    structure: str = typer.Argument(
        ...,
        help="Path to PDB/mmCIF file or PDB ID (e.g., 1yy9)",
    ),
    probe_radius: float = typer.Option(
        1.4,
        "--probe-radius",
        "-p",
        help="Probe radius in Angstroms",
    ),
    per_residue: bool = typer.Option(
        False,
        "--per-residue",
        help="Show per-residue SASA breakdown",
    ),
    output: Path | None = typer.Option(
        None,
        "--output",
        "-o",
        help="Output file (JSON)",
    ),
):
    """
    Calculate SASA for a structure.

    Examples:
        epitype sasa structure.pdb --per-residue
        epitype sasa 1yy9 --per-residue
    """
    # Resolve structure input (file path or PDB ID)
    try:
        structure_path = resolve_structure(structure)
    except StructureResolutionError as e:
        console.print(f"[red]Error: {e}[/red]")
        raise typer.Exit(1) from e

    struct = parse_structure(structure_path)
    result = calculate_sasa(struct, probe_radius)

    if output:
        data = {
            "total": result.total,
            "polar": result.polar,
            "apolar": result.apolar,
        }
        if per_residue:
            data["per_residue"] = result.per_residue

        with open(output, "w") as f:
            json.dump(data, f, indent=2)
        console.print(f"Results written to {output}")
    else:
        table = Table(title=f"SASA: {structure_path.name}")
        table.add_column("Metric", style="cyan")
        table.add_column("Value (A^2)", style="green", justify="right")

        table.add_row("Total SASA", f"{result.total:.1f}")
        table.add_row("Polar SASA", f"{result.polar:.1f}")
        table.add_row("Apolar SASA", f"{result.apolar:.1f}")

        console.print(table)

        if per_residue:
            console.print("\nPer-residue SASA:")
            for res_id, sasa_val in sorted(result.per_residue.items())[:20]:
                console.print(f"  {res_id}: {sasa_val:.1f} A^2")
            if len(result.per_residue) > 20:
                console.print(f"  ... and {len(result.per_residue) - 20} more")
