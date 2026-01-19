"""Shared CLI utilities and options."""
from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer

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
