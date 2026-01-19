"""Core type definitions and constants."""
from dataclasses import dataclass, field
from enum import Enum

import numpy as np
from numpy.typing import NDArray

# Type aliases
Coords = NDArray[np.float64]  # Shape: (N, 3)
AtomIndex = int
ResidueIndex = int
ChainID = str

# Physical constants
PROBE_RADIUS = 1.4  # Water probe radius in Angstroms
DEFAULT_INTERFACE_CUTOFF = 8.0  # Cb-Cb distance cutoff
SEPARATION_DISTANCE = 500.0  # Distance to separate chains


class AtomType(Enum):
    """Atom classification for SASA calculations."""

    POLAR = "polar"
    NONPOLAR = "nonpolar"
    CHARGED = "charged"


@dataclass(frozen=True)
class VdWRadii:
    """Van der Waals radii for common atoms (Angstroms)."""

    C: float = 1.70
    N: float = 1.55
    O: float = 1.52
    S: float = 1.80
    H: float = 1.20
    P: float = 1.80

    def get(self, element: str, default: float = 1.70) -> float:
        """Get radius for element."""
        return getattr(self, element.upper(), default)


VDW_RADII = VdWRadii()


@dataclass
class InterfaceMetrics:
    """Complete interface analysis results."""

    dG_separated: float  # Binding energy in REU
    dSASA_int: float  # Buried surface area in Angstrom^2
    dG_per_dSASA: float  # Energy density
    sc_value: float  # Shape complementarity 0-1
    packstat: float  # Packing quality 0-1
    delta_unsat_hbonds: int  # Buried unsatisfied H-bonds
    hbonds_int: int  # Cross-interface H-bonds

    # Additional details
    interface_residues_1: list[str] = field(default_factory=list)
    interface_residues_2: list[str] = field(default_factory=list)
    sasa_complex: float = 0.0
    sasa_separated: float = 0.0

    def to_dict(self) -> dict:
        """Convert to dictionary for serialization."""
        return {
            "dG_separated": self.dG_separated,
            "dSASA_int": self.dSASA_int,
            "dG_per_dSASA": self.dG_per_dSASA,
            "sc_value": self.sc_value,
            "packstat": self.packstat,
            "delta_unsat_hbonds": self.delta_unsat_hbonds,
            "hbonds_int": self.hbonds_int,
            "interface_residues_1": self.interface_residues_1,
            "interface_residues_2": self.interface_residues_2,
            "sasa_complex": self.sasa_complex,
            "sasa_separated": self.sasa_separated,
        }
