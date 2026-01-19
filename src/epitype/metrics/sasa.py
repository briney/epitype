"""Solvent Accessible Surface Area calculations using FreeSASA."""
from __future__ import annotations

from dataclasses import dataclass

import freesasa

from epitype.core.interface import InterfaceRegion
from epitype.core.structure import Atom, Structure
from epitype.core.types import PROBE_RADIUS


@dataclass
class SASAResult:
    """SASA calculation results."""
    total: float  # Total SASA in Ų
    polar: float  # Polar atom SASA
    apolar: float  # Apolar atom SASA
    per_atom: dict[int, float]  # SASA per atom index
    per_residue: dict[str, float]  # SASA per residue ID


def calculate_sasa(
    structure: Structure,
    probe_radius: float = PROBE_RADIUS,
    algorithm: str = "ShrakeRupley",
) -> SASAResult:
    """
    Calculate SASA for a structure using FreeSASA.

    Args:
        structure: Input structure
        probe_radius: Probe radius in Angstroms
        algorithm: "ShrakeRupley" or "LeeRichards"

    Returns:
        SASAResult with total and per-atom SASA
    """
    # Prepare data for FreeSASA
    coords: list[float] = []
    radii: list[float] = []
    atom_list: list[Atom] = []

    for atom in structure.atoms:
        coords.extend(atom.coords.tolist())
        radii.append(atom.vdw_radius)
        atom_list.append(atom)

    if not atom_list:
        return SASAResult(
            total=0.0,
            polar=0.0,
            apolar=0.0,
            per_atom={},
            per_residue={},
        )

    # Configure FreeSASA
    params = freesasa.Parameters()
    params.setProbeRadius(probe_radius)

    if algorithm == "LeeRichards":
        params.setAlgorithm(freesasa.LeeRichards)
    else:
        params.setAlgorithm(freesasa.ShrakeRupley)

    # Calculate
    result = freesasa.calcCoord(coords, radii, params)

    # Extract per-atom SASA
    per_atom: dict[int, float] = {}
    per_residue: dict[str, float] = {}
    polar_sasa = 0.0
    apolar_sasa = 0.0

    polar_elements = {"N", "O", "S"}

    for i, atom in enumerate(atom_list):
        sasa = result.atomArea(i)
        per_atom[atom.index] = sasa

        # Accumulate per-residue
        res_id = f"{atom.chain_id}:{atom.residue_index}"
        per_residue[res_id] = per_residue.get(res_id, 0.0) + sasa

        # Polar vs apolar
        if atom.element in polar_elements:
            polar_sasa += sasa
        else:
            apolar_sasa += sasa

    return SASAResult(
        total=result.totalArea(),
        polar=polar_sasa,
        apolar=apolar_sasa,
        per_atom=per_atom,
        per_residue=per_residue,
    )


def calculate_dsasa(
    complex_structure: Structure,
    separated_structure: Structure,
    interface: InterfaceRegion,
    probe_radius: float = PROBE_RADIUS,
) -> float:
    """
    Calculate change in SASA upon complex formation (dSASA).

    dSASA = SASA_separated - SASA_complex

    Args:
        complex_structure: Bound complex
        separated_structure: Separated chains
        interface: Interface region
        probe_radius: Probe radius in Angstroms

    Returns:
        dSASA (buried surface area) in Ų
    """
    sasa_complex = calculate_sasa(complex_structure, probe_radius)
    sasa_separated = calculate_sasa(separated_structure, probe_radius)

    return sasa_separated.total - sasa_complex.total


def calculate_interface_sasa(
    structure: Structure,
    interface: InterfaceRegion,
    probe_radius: float = PROBE_RADIUS,
) -> tuple[float, float]:
    """
    Calculate SASA for interface residues only.

    Args:
        structure: Input structure
        interface: Interface region
        probe_radius: Probe radius

    Returns:
        Tuple of (group1_interface_sasa, group2_interface_sasa)
    """
    sasa_result = calculate_sasa(structure, probe_radius)

    g1_sasa = 0.0
    for residue in interface.residues_group1:
        res_id = f"{residue.chain_id}:{residue.index}"
        g1_sasa += sasa_result.per_residue.get(res_id, 0.0)

    g2_sasa = 0.0
    for residue in interface.residues_group2:
        res_id = f"{residue.chain_id}:{residue.index}"
        g2_sasa += sasa_result.per_residue.get(res_id, 0.0)

    return g1_sasa, g2_sasa
