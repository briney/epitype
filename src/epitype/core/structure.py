"""Molecular structure representation."""
from __future__ import annotations

from collections.abc import Iterator, Sequence
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

from .types import VDW_RADII


@dataclass
class Atom:
    """Single atom representation."""

    index: int
    name: str  # e.g., "CA", "CB", "N"
    element: str  # e.g., "C", "N", "O"
    coords: NDArray[np.float64]  # Shape: (3,)
    residue_index: int
    chain_id: str
    b_factor: float = 0.0
    occupancy: float = 1.0

    @property
    def vdw_radius(self) -> float:
        """Van der Waals radius."""
        return VDW_RADII.get(self.element)

    @property
    def is_backbone(self) -> bool:
        """Check if backbone atom."""
        return self.name in {"N", "CA", "C", "O"}

    @property
    def is_hydrogen(self) -> bool:
        """Check if hydrogen atom."""
        return self.element == "H"

    def distance_to(self, other: Atom) -> float:
        """Euclidean distance to another atom."""
        return float(np.linalg.norm(self.coords - other.coords))


@dataclass
class Residue:
    """Single residue representation."""

    index: int
    name: str  # Three-letter code, e.g., "ALA"
    chain_id: str
    seq_num: int  # PDB residue number
    insertion_code: str = ""
    atoms: list[Atom] = field(default_factory=list)

    @property
    def ca(self) -> Atom | None:
        """Get alpha carbon."""
        for atom in self.atoms:
            if atom.name == "CA":
                return atom
        return None

    @property
    def cb(self) -> Atom | None:
        """Get beta carbon (CA for glycine)."""
        for atom in self.atoms:
            if atom.name == "CB":
                return atom
        # Fall back to CA for glycine
        return self.ca

    @property
    def center_of_mass(self) -> NDArray[np.float64]:
        """Calculate residue center of mass."""
        if not self.atoms:
            return np.zeros(3)
        coords = np.array([a.coords for a in self.atoms])
        return coords.mean(axis=0)

    @property
    def full_id(self) -> str:
        """Full residue identifier: chain:name:seqnum."""
        return f"{self.chain_id}:{self.name}:{self.seq_num}{self.insertion_code}"

    def get_atom(self, name: str) -> Atom | None:
        """Get atom by name."""
        for atom in self.atoms:
            if atom.name == name:
                return atom
        return None


@dataclass
class Chain:
    """Single chain representation."""

    chain_id: str
    residues: list[Residue] = field(default_factory=list)

    @property
    def atoms(self) -> Iterator[Atom]:
        """Iterate over all atoms in chain."""
        for residue in self.residues:
            yield from residue.atoms

    @property
    def coords(self) -> NDArray[np.float64]:
        """All atom coordinates as array."""
        atom_list = list(self.atoms)
        if not atom_list:
            return np.empty((0, 3))
        return np.array([a.coords for a in atom_list])

    def __len__(self) -> int:
        return len(self.residues)


@dataclass
class Structure:
    """Complete molecular structure."""

    name: str
    chains: dict[str, Chain] = field(default_factory=dict)

    @property
    def atoms(self) -> Iterator[Atom]:
        """Iterate over all atoms."""
        for chain in self.chains.values():
            yield from chain.atoms

    @property
    def residues(self) -> Iterator[Residue]:
        """Iterate over all residues."""
        for chain in self.chains.values():
            yield from chain.residues

    @property
    def coords(self) -> NDArray[np.float64]:
        """All atom coordinates as array."""
        atom_list = list(self.atoms)
        if not atom_list:
            return np.empty((0, 3))
        return np.array([a.coords for a in atom_list])

    @property
    def num_atoms(self) -> int:
        """Total number of atoms."""
        return sum(1 for _ in self.atoms)

    @property
    def num_residues(self) -> int:
        """Total number of residues."""
        return sum(1 for _ in self.residues)

    def get_chain(self, chain_id: str) -> Chain | None:
        """Get chain by ID."""
        return self.chains.get(chain_id)

    def get_chains(self, chain_ids: Sequence[str]) -> list[Chain]:
        """Get multiple chains by ID."""
        return [self.chains[cid] for cid in chain_ids if cid in self.chains]

    def subset(self, chain_ids: Sequence[str]) -> Structure:
        """Create new structure with only specified chains."""
        new_structure = Structure(name=f"{self.name}_subset")
        for cid in chain_ids:
            if cid in self.chains:
                new_structure.chains[cid] = self.chains[cid]
        return new_structure

    def copy(self) -> Structure:
        """Deep copy of structure."""
        import copy

        return copy.deepcopy(self)

    def translate(self, vector: NDArray[np.float64]) -> None:
        """Translate all atoms by vector (in-place)."""
        for atom in self.atoms:
            atom.coords = atom.coords + vector
