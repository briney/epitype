"""Hydrogen bond analysis subcommand."""
from __future__ import annotations

import json
from pathlib import Path

import typer
from rich.console import Console
from rich.table import Table

from epitype.core.interface import detect_interface, parse_chain_groups
from epitype.io.parsers import parse_structure
from epitype.metrics.hbonds import analyze_hbonds

console = Console()


def hbonds_command(
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
    cutoff: float = typer.Option(
        8.0,
        "--cutoff",
        help="Interface detection cutoff in Angstroms",
    ),
    show_details: bool = typer.Option(
        False,
        "--details",
        help="Show individual hydrogen bonds",
    ),
    output: Path | None = typer.Option(
        None,
        "--output",
        "-o",
        help="Output file (JSON)",
    ),
):
    """
    Analyze hydrogen bonds at a protein-protein interface.

    Example:
        epitype hbonds structure.pdb --chains HL_A --details
    """
    # Parse structure
    struct = parse_structure(structure)

    # Parse chain specification
    group1, group2 = parse_chain_groups(chains)

    # Detect interface
    interface = detect_interface(struct, group1, group2, cutoff=cutoff)

    if interface.num_interface_residues == 0:
        console.print("[red]Error: No interface detected between specified chains.[/red]")
        raise typer.Exit(1)

    # Analyze hydrogen bonds
    result = analyze_hbonds(struct, interface)

    if output:
        data = {
            "interface_hbonds": result.interface_hbonds,
            "total_hbonds": result.total_hbonds,
            "unsatisfied_donors": result.buried_unsatisfied_donors,
            "unsatisfied_acceptors": result.buried_unsatisfied_acceptors,
        }
        if show_details:
            data["hbond_details"] = [
                {
                    "donor_residue": hb.donor_residue,
                    "donor_atom": hb.donor_atom,
                    "acceptor_residue": hb.acceptor_residue,
                    "acceptor_atom": hb.acceptor_atom,
                    "distance": hb.distance,
                    "angle": hb.angle,
                }
                for hb in result.hbonds
            ]

        with open(output, "w") as f:
            json.dump(data, f, indent=2)
        console.print(f"Results written to {output}")
    else:
        table = Table(title=f"Hydrogen Bonds: {structure.name}")
        table.add_column("Metric", style="cyan")
        table.add_column("Value", style="green", justify="right")

        table.add_row("Interface H-bonds", str(result.interface_hbonds))
        table.add_row("Total H-bonds", str(result.total_hbonds))
        table.add_row("Unsatisfied donors", str(result.buried_unsatisfied_donors))
        table.add_row("Unsatisfied acceptors", str(result.buried_unsatisfied_acceptors))

        console.print(table)

        if show_details and result.hbonds:
            console.print("\nInterface hydrogen bonds:")
            for hb in result.hbonds[:20]:
                console.print(
                    f"  {hb.donor_residue}:{hb.donor_atom} -> "
                    f"{hb.acceptor_residue}:{hb.acceptor_atom} "
                    f"({hb.distance:.2f} A, {hb.angle:.1f} deg)"
                )
            if len(result.hbonds) > 20:
                console.print(f"  ... and {len(result.hbonds) - 20} more")
