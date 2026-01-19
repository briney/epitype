"""Shape complementarity calculation using NanoShaper surfaces."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.spatial import KDTree

from epitype.core.interface import InterfaceRegion
from epitype.core.structure import Structure
from epitype.surface.nanoshaper import generate_surface
from epitype.surface.types import SurfacePoint


@dataclass
class ShapeComplementarityResult:
    """Shape complementarity analysis results."""

    sc_value: float  # Overall Sc (0-1)
    sc_median: float  # Median Sc
    sc_group1: float  # Sc for group1 surface
    sc_group2: float  # Sc for group2 surface
    interface_area: float  # Interface surface area


def calculate_shape_complementarity(
    structure: Structure,
    interface: InterfaceRegion,
    distance_cutoff: float = 1.5,
    trim_distance: float = 10.0,
) -> ShapeComplementarityResult:
    """
    Calculate shape complementarity (Sc) between interface surfaces.

    Algorithm (Lawrence & Colman, 1993):
    1. Generate molecular surfaces for each chain group
    2. For each surface point on one side, find nearest point on other side
    3. Calculate Sc = dot(n1, n2) * exp(-distance^2 / 0.5)
    4. Average over all interface points

    Args:
        structure: Input structure
        interface: Interface region
        distance_cutoff: Max distance between surface points to consider
        trim_distance: Distance from interface to include surface points

    Returns:
        ShapeComplementarityResult
    """
    # Generate surfaces for each group
    group1_structure = structure.subset(interface.group1_chains)
    group2_structure = structure.subset(interface.group2_chains)

    surface1 = generate_surface(group1_structure)
    surface2 = generate_surface(group2_structure)

    if not surface1 or not surface2:
        return ShapeComplementarityResult(
            sc_value=0.0, sc_median=0.0, sc_group1=0.0, sc_group2=0.0, interface_area=0.0
        )

    # Trim to interface region
    interface_surface1 = _trim_surface_to_interface(surface1, surface2, trim_distance)
    interface_surface2 = _trim_surface_to_interface(surface2, surface1, trim_distance)

    if not interface_surface1 or not interface_surface2:
        return ShapeComplementarityResult(
            sc_value=0.0, sc_median=0.0, sc_group1=0.0, sc_group2=0.0, interface_area=0.0
        )

    # Calculate Sc for each direction
    sc_values_1to2 = _calculate_sc_one_direction(
        interface_surface1, interface_surface2, distance_cutoff
    )
    sc_values_2to1 = _calculate_sc_one_direction(
        interface_surface2, interface_surface1, distance_cutoff
    )

    # Combine
    sc_group1 = float(np.mean(sc_values_1to2)) if len(sc_values_1to2) else 0.0
    sc_group2 = float(np.mean(sc_values_2to1)) if len(sc_values_2to1) else 0.0

    if len(sc_values_1to2) and len(sc_values_2to1):
        all_sc = np.concatenate([sc_values_1to2, sc_values_2to1])
    else:
        all_sc = np.array([0.0])

    return ShapeComplementarityResult(
        sc_value=float(np.mean(all_sc)),
        sc_median=float(np.median(all_sc)),
        sc_group1=float(sc_group1),
        sc_group2=float(sc_group2),
        interface_area=float(len(interface_surface1) + len(interface_surface2)),
    )


def _trim_surface_to_interface(
    surface: list[SurfacePoint],
    other_surface: list[SurfacePoint],
    trim_distance: float,
) -> list[SurfacePoint]:
    """Keep only surface points near the other surface."""
    if not other_surface:
        return []

    other_coords = np.array([p.coords for p in other_surface])
    tree = KDTree(other_coords)

    trimmed = []
    for point in surface:
        dist, _ = tree.query(point.coords)
        if dist < trim_distance:
            trimmed.append(point)

    return trimmed


def _calculate_sc_one_direction(
    from_surface: list[SurfacePoint],
    to_surface: list[SurfacePoint],
    distance_cutoff: float,
) -> NDArray[np.float64]:
    """
    Calculate Sc values from one surface to another.

    For each point in from_surface, find nearest in to_surface
    and compute Sc contribution.
    """
    if not from_surface or not to_surface:
        return np.array([])

    to_coords = np.array([p.coords for p in to_surface])
    to_normals = np.array([p.normal for p in to_surface])
    tree = KDTree(to_coords)

    sc_values = []

    for point in from_surface:
        dist, idx = tree.query(point.coords)

        if dist > distance_cutoff:
            continue

        # Get normals
        n1 = point.normal
        n2 = to_normals[idx]

        # Sc = -dot(n1, n2) * exp(-d^2 / sigma^2)
        # Negative because normals should point away from each other
        # sigma^2 = 0.5 is typical
        dot_product = -np.dot(n1, n2)
        distance_weight = np.exp(-dist * dist / 0.5)

        sc = dot_product * distance_weight
        sc_values.append(sc)

    return np.array(sc_values)
