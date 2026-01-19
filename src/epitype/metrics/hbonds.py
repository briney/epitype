"""Hydrogen bond detection and analysis."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.spatial import KDTree

from epitype.core.interface import InterfaceRegion
from epitype.core.structure import Atom, Structure

# Hydrogen bond geometry criteria
HBOND_DISTANCE_MAX = 3.5  # D-A distance in Angstroms
HBOND_ANGLE_MIN = 120.0  # D-H-A angle minimum in degrees


@dataclass
class HydrogenBond:
    """Single hydrogen bond."""
    donor_atom: Atom
    acceptor_atom: Atom
    hydrogen_atom: Atom | None
    distance: float  # D-A distance
    angle: float | None  # D-H-A angle if H available


@dataclass
class HBondResult:
    """Hydrogen bond analysis results."""
    total_hbonds: int
    interface_hbonds: int  # Cross-interface H-bonds
    hbonds: list[HydrogenBond]
    buried_unsatisfied_donors: int
    buried_unsatisfied_acceptors: int


# Donor atoms: N-H groups
DONOR_ATOMS = {
    # Backbone
    "N",
    # Sidechains
    "NE", "NH1", "NH2",  # Arg
    "ND2",  # Asn
    "NE2",  # Gln
    "NZ",  # Lys
    "NE1",  # Trp
    "ND1",  # His (can be either)
    "OG", "OG1",  # Ser, Thr (hydroxyl)
    "OH",  # Tyr
}

# Acceptor atoms: O, N with lone pairs
ACCEPTOR_ATOMS = {
    # Backbone
    "O",
    # Sidechains
    "OD1", "OD2",  # Asp
    "OE1", "OE2",  # Glu
    "ND1", "NE2",  # His
    "OG", "OG1",  # Ser, Thr
    "OH",  # Tyr
    "SD",  # Met (weak)
}


def detect_hydrogen_bonds(
    structure: Structure,
    distance_cutoff: float = HBOND_DISTANCE_MAX,
    angle_cutoff: float = HBOND_ANGLE_MIN,
) -> list[HydrogenBond]:
    """
    Detect all hydrogen bonds in structure.

    Uses geometric criteria:
    - D-A distance < 3.5 Å
    - D-H-A angle > 120° (if H available)

    Args:
        structure: Input structure
        distance_cutoff: Maximum D-A distance
        angle_cutoff: Minimum D-H-A angle

    Returns:
        List of detected hydrogen bonds
    """
    # Collect donors and acceptors
    donors: list[Atom] = []
    acceptors: list[Atom] = []

    for atom in structure.atoms:
        if atom.name in DONOR_ATOMS:
            donors.append(atom)
        if atom.name in ACCEPTOR_ATOMS:
            acceptors.append(atom)

    if not donors or not acceptors:
        return []

    # Build KD-tree for acceptors
    acceptor_coords = np.array([a.coords for a in acceptors])
    tree = KDTree(acceptor_coords)

    hbonds: list[HydrogenBond] = []

    for donor in donors:
        # Find nearby acceptors
        indices = tree.query_ball_point(donor.coords, distance_cutoff)

        for idx in indices:
            acceptor = acceptors[idx]

            # Skip same residue
            if donor.residue_index == acceptor.residue_index:
                continue

            distance = donor.distance_to(acceptor)

            # Create H-bond (angle check would require finding H atom)
            hbond = HydrogenBond(
                donor_atom=donor,
                acceptor_atom=acceptor,
                hydrogen_atom=None,  # Could search for attached H
                distance=distance,
                angle=None,
            )
            hbonds.append(hbond)

    return hbonds


def count_interface_hbonds(
    hbonds: list[HydrogenBond],
    interface: InterfaceRegion,
) -> int:
    """
    Count hydrogen bonds crossing the interface.

    Args:
        hbonds: All detected H-bonds
        interface: Interface region

    Returns:
        Number of cross-interface H-bonds
    """
    g1_chains = set(interface.group1_chains)
    g2_chains = set(interface.group2_chains)

    count = 0
    for hbond in hbonds:
        donor_in_g1 = hbond.donor_atom.chain_id in g1_chains
        acceptor_in_g1 = hbond.acceptor_atom.chain_id in g1_chains
        donor_in_g2 = hbond.donor_atom.chain_id in g2_chains
        acceptor_in_g2 = hbond.acceptor_atom.chain_id in g2_chains

        # Cross-interface if one atom in g1 and other in g2
        if (donor_in_g1 and acceptor_in_g2) or (donor_in_g2 and acceptor_in_g1):
            count += 1

    return count


def count_unsatisfied_hbonds(
    structure: Structure,
    interface: InterfaceRegion,
    hbonds: list[HydrogenBond],
    burial_threshold: float = 0.1,  # SASA threshold for "buried"
) -> tuple[int, int]:
    """
    Count buried unsatisfied H-bond donors and acceptors.

    An atom is unsatisfied if:
    - It's buried (low SASA)
    - It's a potential donor/acceptor
    - It's not participating in an H-bond

    Args:
        structure: Input structure
        interface: Interface region
        hbonds: Detected H-bonds
        burial_threshold: SASA threshold for burial

    Returns:
        Tuple of (unsatisfied_donors, unsatisfied_acceptors)
    """
    from epitype.metrics.sasa import calculate_sasa

    # Get per-atom SASA
    sasa_result = calculate_sasa(structure)

    # Find atoms in H-bonds
    satisfied_donors = {hb.donor_atom.index for hb in hbonds}
    satisfied_acceptors = {hb.acceptor_atom.index for hb in hbonds}

    # Get interface atom indices
    interface_atoms = set(interface.atoms_group1_indices) | set(interface.atoms_group2_indices)

    unsatisfied_donors = 0
    unsatisfied_acceptors = 0

    for atom in structure.atoms:
        # Only check interface atoms
        if atom.index not in interface_atoms:
            continue

        # Check if buried
        atom_sasa = sasa_result.per_atom.get(atom.index, 0.0)
        if atom_sasa > burial_threshold:
            continue

        # Check donors
        if atom.name in DONOR_ATOMS and atom.index not in satisfied_donors:
            unsatisfied_donors += 1

        # Check acceptors
        if atom.name in ACCEPTOR_ATOMS and atom.index not in satisfied_acceptors:
            unsatisfied_acceptors += 1

    return unsatisfied_donors, unsatisfied_acceptors


def analyze_hbonds(
    structure: Structure,
    interface: InterfaceRegion,
) -> HBondResult:
    """
    Complete hydrogen bond analysis.

    Args:
        structure: Input structure
        interface: Interface region

    Returns:
        HBondResult with all H-bond metrics
    """
    hbonds = detect_hydrogen_bonds(structure)
    interface_count = count_interface_hbonds(hbonds, interface)
    unsat_donors, unsat_acceptors = count_unsatisfied_hbonds(
        structure, interface, hbonds
    )

    return HBondResult(
        total_hbonds=len(hbonds),
        interface_hbonds=interface_count,
        hbonds=hbonds,
        buried_unsatisfied_donors=unsat_donors,
        buried_unsatisfied_acceptors=unsat_acceptors,
    )
