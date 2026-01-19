"""Import utilities with helpful error messages for optional dependencies."""
from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pdbfixer import PDBFixer as PDBFixerType


class PDBFixerNotFoundError(ImportError):
    """Raised when pdbfixer is required but not installed."""

    pass


def import_pdbfixer() -> type[PDBFixerType]:
    """
    Import pdbfixer with helpful error if not available.

    pdbfixer is only available via conda/mamba, not PyPI.

    Returns:
        The PDBFixer class

    Raises:
        PDBFixerNotFoundError: If pdbfixer is not installed
    """
    try:
        from pdbfixer import PDBFixer

        return PDBFixer
    except ImportError:
        raise PDBFixerNotFoundError(
            "pdbfixer is required for energy calculations but is not installed.\n"
            "\n"
            "pdbfixer is only available via conda/mamba:\n"
            "\n"
            "    conda install -c conda-forge pdbfixer\n"
            "    # or\n"
            "    mamba install -c conda-forge pdbfixer\n"
            "\n"
            "For more information, see: https://github.com/openmm/pdbfixer"
        ) from None
