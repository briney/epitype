"""Rigid body separation for binding energy calculations."""
from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from .interface import InterfaceRegion
from .structure import Structure
from .types import SEPARATION_DISTANCE


def compute_interface_normal(
    structure: Structure,
    interface: InterfaceRegion,
) -> NDArray[np.float64]:
    """
    Compute the interface normal vector.

    The normal points from group1 centroid to group2 centroid.

    Args:
        structure: Input structure
        interface: Detected interface region

    Returns:
        Unit normal vector (3,)
    """
    # Get centroids of interface residues
    coords_g1 = []
    for residue in interface.residues_group1:
        coords_g1.append(residue.center_of_mass)

    coords_g2 = []
    for residue in interface.residues_group2:
        coords_g2.append(residue.center_of_mass)

    if not coords_g1 or not coords_g2:
        # Default to z-axis if no interface
        return np.array([0.0, 0.0, 1.0])

    centroid_g1 = np.mean(coords_g1, axis=0)
    centroid_g2 = np.mean(coords_g2, axis=0)

    normal = centroid_g2 - centroid_g1
    norm = np.linalg.norm(normal)

    if norm < 1e-6:
        return np.array([0.0, 0.0, 1.0])

    return normal / norm


def separate_chains(
    structure: Structure,
    interface: InterfaceRegion,
    distance: float = SEPARATION_DISTANCE,
) -> Structure:
    """
    Create separated structure by translating group2 chains.

    Args:
        structure: Input structure (will be copied)
        interface: Interface region defining chain groups
        distance: Separation distance in Angstroms

    Returns:
        New structure with group2 chains translated
    """
    # Deep copy the structure
    separated = structure.copy()

    # Compute interface normal
    normal = compute_interface_normal(structure, interface)

    # Translation vector
    translation = normal * distance

    # Translate group2 chains
    for chain_id in interface.group2_chains:
        chain = separated.get_chain(chain_id)
        if chain is not None:
            for residue in chain.residues:
                for atom in residue.atoms:
                    atom.coords = atom.coords + translation

    return separated


def compute_separation_distance(
    structure: Structure,
    interface: InterfaceRegion,
    target_min_distance: float = 100.0,
) -> float:
    """
    Compute minimum separation distance to eliminate all contacts.

    Args:
        structure: Input structure
        interface: Interface region
        target_min_distance: Target minimum inter-group distance

    Returns:
        Required separation distance
    """
    from scipy.spatial.distance import cdist

    # Get all atom coordinates for each group
    coords_g1 = []
    coords_g2 = []

    for chain_id in interface.group1_chains:
        chain = structure.get_chain(chain_id)
        if chain:
            for atom in chain.atoms:
                coords_g1.append(atom.coords)

    for chain_id in interface.group2_chains:
        chain = structure.get_chain(chain_id)
        if chain:
            for atom in chain.atoms:
                coords_g2.append(atom.coords)

    if not coords_g1 or not coords_g2:
        return SEPARATION_DISTANCE

    # Compute current minimum distance
    distances = cdist(coords_g1, coords_g2)
    current_min = distances.min()

    # Required additional separation
    return max(target_min_distance - current_min, 0) + SEPARATION_DISTANCE
