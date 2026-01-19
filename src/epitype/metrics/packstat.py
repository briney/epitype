"""PackStat packing quality metric using multi-probe algorithm.

This implementation uses FreeSASA for fast SASA calculations at multiple
probe radii, providing ~100x speedup over custom Python implementations.
"""
from __future__ import annotations

from dataclasses import dataclass

import freesasa
import numpy as np
from numpy.typing import NDArray

from epitype.core.interface import InterfaceRegion
from epitype.core.structure import Structure

# Probe radii for multi-probe SASA (30 values from 0.1 to 3.0 Å)
PROBE_RADII = np.linspace(0.1, 3.0, 30)


@dataclass
class PackStatResult:
    """PackStat analysis results."""
    packstat: float  # Overall packing score (0-1)
    packstat_interface: float  # Interface-specific packing
    cavity_volume: float  # Estimated cavity volume (ų)
    per_probe_scores: dict[float, float]  # Score at each probe radius


def calculate_packstat(
    structure: Structure,
    interface: InterfaceRegion | None = None,
) -> PackStatResult:
    """
    Calculate PackStat using multi-probe algorithm.

    Uses FreeSASA for fast SASA calculations at multiple probe radii.

    Algorithm:
    1. For each probe radius (0.1 to 3.0 Å, 30 values):
       a. Compute SASA with that probe using FreeSASA
       b. Calculate packing fraction (buried / reference surface)
    2. Weight and combine across probe radii (favoring smaller probes)

    Args:
        structure: Input structure
        interface: Optional interface region (for interface-specific score)

    Returns:
        PackStatResult with packing metrics
    """
    # Get all atom data
    atoms = list(structure.atoms)

    if len(atoms) == 0:
        return PackStatResult(
            packstat=0.0, packstat_interface=0.0, cavity_volume=0.0, per_probe_scores={}
        )

    # Prepare data for FreeSASA
    coords: list[float] = []
    radii: list[float] = []

    for atom in atoms:
        coords.extend(atom.coords.tolist())
        radii.append(atom.vdw_radius)

    radii_arr = np.array(radii)

    # Calculate SASA at each probe radius
    per_probe_scores: dict[float, float] = {}
    contact_areas: list[float] = []

    for probe_r in PROBE_RADII:
        # SASA at this probe radius using FreeSASA
        total_sasa = _freesasa_total(coords, radii, probe_r)

        # Reference SASA (isolated atoms)
        ref_sasa = _reference_sasa(radii_arr, probe_r)

        # Packing fraction: how much surface is buried
        if ref_sasa > 0:
            packing = 1.0 - (total_sasa / ref_sasa)
        else:
            packing = 0.0

        per_probe_scores[float(probe_r)] = max(0.0, min(1.0, packing))
        contact_areas.append(packing)

    # Overall packstat: weighted average favoring smaller probes
    # (smaller probes detect finer packing details)
    weights = np.exp(-PROBE_RADII)  # Exponential decay with probe size
    weights /= weights.sum()

    packstat = float(np.average(contact_areas, weights=weights))

    # Interface-specific packing
    packstat_interface = packstat
    if interface is not None:
        packstat_interface = _interface_packstat(structure, interface)

    # Estimate cavity volume from intermediate probe sizes
    cavity_volume = _estimate_cavity_volume(radii_arr, contact_areas)

    return PackStatResult(
        packstat=packstat,
        packstat_interface=packstat_interface,
        cavity_volume=cavity_volume,
        per_probe_scores=per_probe_scores,
    )


def _freesasa_total(
    coords: list[float],
    radii: list[float],
    probe_radius: float,
) -> float:
    """
    Calculate total SASA using FreeSASA.

    Args:
        coords: Flattened coordinates [x1, y1, z1, x2, y2, z2, ...]
        radii: Van der Waals radii for each atom
        probe_radius: Probe radius in Angstroms

    Returns:
        Total SASA in Ų
    """
    params = freesasa.Parameters()
    params.setProbeRadius(probe_radius)
    params.setAlgorithm(freesasa.ShrakeRupley)
    result = freesasa.calcCoord(coords, radii, params)
    return result.totalArea()


def _reference_sasa(radii: NDArray[np.float64], probe_radius: float) -> float:
    """Calculate reference SASA for isolated atoms."""
    eff_radii = radii + probe_radius
    return float(np.sum(4 * np.pi * eff_radii ** 2))


def _interface_packstat(
    structure: Structure,
    interface: InterfaceRegion,
) -> float:
    """Calculate packstat for interface atoms only."""
    # Get interface atom indices
    interface_indices = set(interface.atoms_group1_indices) | set(
        interface.atoms_group2_indices
    )

    atoms = [a for a in structure.atoms if a.index in interface_indices]
    if not atoms:
        return 0.0

    # Prepare data for FreeSASA
    coords: list[float] = []
    radii: list[float] = []

    for atom in atoms:
        coords.extend(atom.coords.tolist())
        radii.append(atom.vdw_radius)

    radii_arr = np.array(radii)

    scores = []
    weights = []

    for probe_r in PROBE_RADII:
        total_sasa = _freesasa_total(coords, radii, probe_r)
        ref_sasa = _reference_sasa(radii_arr, probe_r)

        if ref_sasa > 0:
            packing = 1.0 - (total_sasa / ref_sasa)
        else:
            packing = 0.0

        scores.append(max(0.0, min(1.0, packing)))
        weights.append(np.exp(-probe_r))

    weights_arr = np.array(weights)
    weights_arr /= weights_arr.sum()

    return float(np.average(scores, weights=weights_arr))


def _estimate_cavity_volume(
    radii: NDArray[np.float64],
    packing_scores: list[float],
) -> float:
    """Estimate cavity volume from packing scores."""
    # Rough estimate: volume deficit based on packing scores
    total_vdw_volume = float(np.sum((4 / 3) * np.pi * radii ** 3))

    # Average packing deficit at intermediate probe sizes
    mid_scores = packing_scores[10:20]  # Probe radii ~1.0-2.0 Å
    if mid_scores:
        avg_deficit = 1.0 - np.mean(mid_scores)
    else:
        avg_deficit = 0.0

    return total_vdw_volume * avg_deficit * 0.1  # Scaling factor
