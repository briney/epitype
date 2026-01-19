"""Core data structures and types."""
from .interface import InterfaceRegion, detect_interface, parse_chain_groups
from .separation import (
    compute_interface_normal,
    compute_separation_distance,
    separate_chains,
)
from .structure import Atom, Chain, Residue, Structure
from .types import (
    DEFAULT_INTERFACE_CUTOFF,
    PROBE_RADIUS,
    SEPARATION_DISTANCE,
    VDW_RADII,
    AtomIndex,
    AtomType,
    ChainID,
    Coords,
    InterfaceMetrics,
    ResidueIndex,
    VdWRadii,
)

__all__ = [
    # Types
    "AtomIndex",
    "AtomType",
    "ChainID",
    "Coords",
    "DEFAULT_INTERFACE_CUTOFF",
    "InterfaceMetrics",
    "PROBE_RADIUS",
    "ResidueIndex",
    "SEPARATION_DISTANCE",
    "VDW_RADII",
    "VdWRadii",
    # Structure
    "Atom",
    "Chain",
    "Residue",
    "Structure",
    # Interface
    "InterfaceRegion",
    "detect_interface",
    "parse_chain_groups",
    # Separation
    "compute_interface_normal",
    "compute_separation_distance",
    "separate_chains",
]
