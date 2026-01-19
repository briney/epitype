"""Energy calculation subcommand."""
from __future__ import annotations

import json
from pathlib import Path

import typer
from rich.console import Console
from rich.table import Table

console = Console()


def energy_command(
    structure: Path = typer.Argument(
        ...,
        help="Path to PDB or mmCIF file",
        exists=True,
    ),
    chains: str = typer.Option(
        ...,
        "--chains",
        "-c",
        help="Chain specification (e.g., 'HL_A')",
    ),
    minimize: bool = typer.Option(
        True,
        "--minimize/--no-minimize",
        help="Minimize energy before calculation",
    ),
    output: Path | None = typer.Option(
        None,
        "--output",
        "-o",
        help="Output file (JSON)",
    ),
):
    """
    Calculate binding energy for a protein-protein interface.

    Requires OpenMM and PDBFixer to be installed.

    Example:
        epitype energy structure.pdb --chains HL_A
    """
    from epitype.metrics.energy import (
        calculate_binding_energy,
        check_openmm_available,
    )

    if not check_openmm_available():
        console.print("[red]Error: OpenMM is not available.[/red]")
        console.print("Install with: conda install -c conda-forge openmm pdbfixer")
        raise typer.Exit(1)

    from epitype.core.interface import detect_interface, parse_chain_groups
    from epitype.core.separation import separate_chains
    from epitype.io.parsers import parse_structure

    # Parse structure
    struct = parse_structure(structure)

    # Parse chain specification
    group1, group2 = parse_chain_groups(chains)

    # Detect interface
    interface = detect_interface(struct, group1, group2)

    if interface.num_interface_residues == 0:
        console.print("[red]Error: No interface detected between specified chains.[/red]")
        raise typer.Exit(1)

    # Separate chains
    separated = separate_chains(struct, interface)

    # Calculate binding energy
    console.print("Calculating binding energy (this may take a moment)...")
    result = calculate_binding_energy(struct, separated, minimize=minimize)

    if output:
        data = {
            "dG_separated": result.dG_separated,
            "dG_per_dSASA": result.dG_per_dSASA,
            "energy_complex": result.energy_complex,
            "energy_separated": result.energy_separated,
        }

        with open(output, "w") as f:
            json.dump(data, f, indent=2)
        console.print(f"Results written to {output}")
    else:
        table = Table(title=f"Binding Energy: {structure.name}")
        table.add_column("Metric", style="cyan")
        table.add_column("Value", style="green", justify="right")

        table.add_row("dG_separated (REU)", f"{result.dG_separated:.2f}")
        table.add_row("Complex energy (kJ/mol)", f"{result.energy_complex:.2f}")
        table.add_row("Separated energy (kJ/mol)", f"{result.energy_separated:.2f}")

        console.print(table)
