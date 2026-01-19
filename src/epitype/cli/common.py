"""Shared CLI utilities and options."""
from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console

from epitype.io.download import PDBDownloadError, download_pdb, is_pdb_id

console = Console()


class StructureResolutionError(Exception):
    """Raised when structure cannot be resolved from input."""

    pass


def resolve_structure(structure_input: str) -> Path:
    """
    Resolve structure input to a valid file path.

    This function handles three cases:
    1. Existing file path -> return as Path
    2. Valid PDB ID -> download from RCSB and return path
    3. Neither -> raise error

    Args:
        structure_input: File path or PDB ID string

    Returns:
        Path to structure file (existing or downloaded)

    Raises:
        StructureResolutionError: If input cannot be resolved to a file
    """
    path = Path(structure_input)

    # Case 1: Existing file
    if path.exists():
        if not path.is_file():
            raise StructureResolutionError(
                f"'{structure_input}' is a directory, not a file"
            )
        return path

    # Case 2: Check if it's a PDB ID
    if is_pdb_id(structure_input):
        console.print(f"[dim]Downloading {structure_input} from RCSB...[/dim]")
        try:
            return download_pdb(structure_input)
        except PDBDownloadError as e:
            raise StructureResolutionError(str(e)) from e

    # Case 3: Neither file nor valid PDB ID
    raise StructureResolutionError(
        f"'{structure_input}' is not an existing file or valid PDB ID.\n"
        f"PDB IDs must be 4 alphanumeric characters (e.g., '1yy9') "
        f"or extended format (e.g., 'pdb_00001yy9')."
    )


# Common option types
StructurePath = Annotated[
    Path,
    typer.Argument(
        help="Path to PDB or mmCIF file",
        exists=True,
        readable=True,
    ),
]

ChainSpec = Annotated[
    str,
    typer.Option(
        "--chains",
        "-c",
        help="Chain specification (e.g., 'HL_A')",
    ),
]

OutputPath = Annotated[
    Path | None,
    typer.Option(
        "--output",
        "-o",
        help="Output file path",
    ),
]

Verbose = Annotated[
    bool,
    typer.Option(
        "--verbose",
        "-v",
        help="Enable verbose output",
    ),
]


def validate_chain_spec(chain_spec: str) -> tuple[list[str], list[str]]:
    """Validate and parse chain specification."""
    from epitype.core.interface import parse_chain_groups

    try:
        return parse_chain_groups(chain_spec)
    except ValueError as e:
        raise typer.BadParameter(str(e)) from e
