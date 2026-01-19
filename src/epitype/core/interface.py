"""Interface detection and representation."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.spatial import KDTree

from .structure import Residue, Structure
from .types import DEFAULT_INTERFACE_CUTOFF


@dataclass
class InterfaceRegion:
    """Represents the interface between two chain groups."""

    # Chain groups (e.g., ["H", "L"] vs ["A"])
    group1_chains: list[str]
    group2_chains: list[str]

    # Interface residues from each group
    residues_group1: list[Residue] = field(default_factory=list)
    residues_group2: list[Residue] = field(default_factory=list)

    # Interface atoms
    atoms_group1_indices: list[int] = field(default_factory=list)
    atoms_group2_indices: list[int] = field(default_factory=list)

    # Contact pairs: (residue1_idx, residue2_idx)
    contact_pairs: list[tuple[int, int]] = field(default_factory=list)

    @property
    def num_interface_residues(self) -> int:
        """Total interface residues."""
        return len(self.residues_group1) + len(self.residues_group2)

    @property
    def residue_ids_group1(self) -> list[str]:
        """Full IDs of group1 interface residues."""
        return [r.full_id for r in self.residues_group1]

    @property
    def residue_ids_group2(self) -> list[str]:
        """Full IDs of group2 interface residues."""
        return [r.full_id for r in self.residues_group2]


def detect_interface(
    structure: Structure,
    group1_chains: list[str],
    group2_chains: list[str],
    cutoff: float = DEFAULT_INTERFACE_CUTOFF,
    use_cb: bool = True,
) -> InterfaceRegion:
    """
    Detect interface residues between two chain groups.

    Args:
        structure: Input structure
        group1_chains: Chain IDs for first group (e.g., ["H", "L"])
        group2_chains: Chain IDs for second group (e.g., ["A"])
        cutoff: Distance cutoff in Angstroms (Cb-Cb or CA-CA)
        use_cb: If True, use Cb distances; otherwise CA

    Returns:
        InterfaceRegion with detected interface residues
    """
    interface = InterfaceRegion(
        group1_chains=group1_chains,
        group2_chains=group2_chains,
    )

    # Collect residues and representative coordinates for each group
    residues_g1: list[Residue] = []
    coords_g1: list[np.ndarray] = []

    residues_g2: list[Residue] = []
    coords_g2: list[np.ndarray] = []

    for chain_id in group1_chains:
        chain = structure.get_chain(chain_id)
        if chain is None:
            continue
        for residue in chain.residues:
            rep_atom = residue.cb if use_cb else residue.ca
            if rep_atom is not None:
                residues_g1.append(residue)
                coords_g1.append(rep_atom.coords)

    for chain_id in group2_chains:
        chain = structure.get_chain(chain_id)
        if chain is None:
            continue
        for residue in chain.residues:
            rep_atom = residue.cb if use_cb else residue.ca
            if rep_atom is not None:
                residues_g2.append(residue)
                coords_g2.append(rep_atom.coords)

    if not coords_g1 or not coords_g2:
        return interface

    coords_g1_arr = np.array(coords_g1)
    coords_g2_arr = np.array(coords_g2)

    # Build KD-tree for group2 and query with group1
    tree = KDTree(coords_g2_arr)

    interface_g1_indices: set[int] = set()
    interface_g2_indices: set[int] = set()
    contact_pairs: list[tuple[int, int]] = []

    # Query all group1 points against group2
    for i, coord in enumerate(coords_g1_arr):
        neighbors = tree.query_ball_point(coord, cutoff)
        if neighbors:
            interface_g1_indices.add(i)
            for j in neighbors:
                interface_g2_indices.add(j)
                contact_pairs.append((i, j))

    # Populate interface
    interface.residues_group1 = [residues_g1[i] for i in sorted(interface_g1_indices)]
    interface.residues_group2 = [residues_g2[i] for i in sorted(interface_g2_indices)]
    interface.contact_pairs = contact_pairs

    # Get atom indices for interface residues
    for residue in interface.residues_group1:
        interface.atoms_group1_indices.extend([a.index for a in residue.atoms])
    for residue in interface.residues_group2:
        interface.atoms_group2_indices.extend([a.index for a in residue.atoms])

    return interface


def parse_chain_groups(chain_spec: str) -> tuple[list[str], list[str]]:
    """
    Parse chain group specification.

    Format: "HL_A" means chains H+L vs chain A
    Format: "AB_CD" means chains A+B vs chains C+D

    Args:
        chain_spec: Chain specification string

    Returns:
        Tuple of (group1_chains, group2_chains)

    Raises:
        ValueError: If chain_spec is invalid
    """
    if "_" not in chain_spec:
        raise ValueError(f"Invalid chain spec '{chain_spec}'. Expected format: 'HL_A'")

    parts = chain_spec.split("_")
    if len(parts) != 2:
        raise ValueError(
            f"Invalid chain spec '{chain_spec}'. Expected exactly one underscore."
        )

    group1 = list(parts[0])
    group2 = list(parts[1])

    if not group1 or not group2:
        raise ValueError(
            f"Invalid chain spec '{chain_spec}'. Both groups must have at least one chain."
        )

    return group1, group2
