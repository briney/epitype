"""Binding energy calculations using OpenMM with AMBER ff14SB."""
from __future__ import annotations

import tempfile
from pathlib import Path

from epitype.core.interface import InterfaceRegion
from epitype.core.structure import Structure

# Lazy imports for optional dependencies
_openmm_available = None


def _check_openmm_available() -> bool:
    """Check if OpenMM is available."""
    global _openmm_available
    if _openmm_available is None:
        try:
            import openmm
            from openmm import app, unit
            from pdbfixer import PDBFixer
            _openmm_available = True
        except ImportError:
            _openmm_available = False
    return _openmm_available


# Conversion: kJ/mol to kcal/mol
KJ_TO_REU = 0.239  # 1 kJ/mol ~ 0.239 kcal/mol

# Standard amino acids recognized by AMBER force field
STANDARD_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    # Common variants
    "HIE", "HID", "HIP",  # Histidine protonation states
    "CYX",  # Disulfide-bonded cysteine
}


def prepare_structure_for_openmm(structure: Structure) -> tuple:
    """
    Prepare structure for OpenMM energy calculation.

    Includes fixing missing atoms/residues via PDBFixer.

    Args:
        structure: Input structure

    Returns:
        Tuple of (Modeller, ForceField)
    """
    if not _check_openmm_available():
        raise RuntimeError(
            "OpenMM is not available. Please install it with: "
            "pip install openmm pdbfixer"
        )

    from openmm.app import ForceField, Modeller
    from pdbfixer import PDBFixer

    # Write to temporary PDB
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False, mode='w') as f:
        temp_path = Path(f.name)
        _write_temp_pdb(structure, temp_path)

    # Fix structure
    fixer = PDBFixer(filename=str(temp_path))
    fixer.findMissingResidues()
    # Clear missing residues - we only want to fix atoms within existing residues,
    # not add new residues at termini or internal gaps
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.0)

    # Clean up temp file
    temp_path.unlink()

    # Create force field
    forcefield = ForceField("amber14-all.xml", "amber14/tip3p.xml")

    # Create modeller (no solvent for interface calculations)
    modeller = Modeller(fixer.topology, fixer.positions)

    return modeller, forcefield


def calculate_potential_energy(
    modeller,
    forcefield,
    minimize: bool = False,
) -> float:
    """
    Calculate potential energy of a structure.

    Args:
        modeller: OpenMM Modeller with structure
        forcefield: OpenMM ForceField
        minimize: If True, minimize before calculating energy

    Returns:
        Potential energy in kJ/mol
    """
    if not _check_openmm_available():
        raise RuntimeError(
            "OpenMM is not available. Please install it with: "
            "pip install openmm pdbfixer"
        )

    from openmm import LangevinMiddleIntegrator, Platform, app, unit

    # Create system
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.NoCutoff,  # No cutoff for accuracy
        constraints=None,
        rigidWater=False,
    )

    # Create integrator (not used for single-point energy)
    integrator = LangevinMiddleIntegrator(
        300 * unit.kelvin,
        1.0 / unit.picosecond,
        0.002 * unit.picoseconds,
    )

    # Create simulation
    platform = Platform.getPlatformByName("CPU")
    simulation = app.Simulation(modeller.topology, system, integrator, platform)
    simulation.context.setPositions(modeller.positions)

    if minimize:
        simulation.minimizeEnergy(maxIterations=100)

    # Get energy
    state = simulation.context.getState(getEnergy=True)
    energy = state.getPotentialEnergy()

    return energy.value_in_unit(unit.kilojoule_per_mole)


def calculate_binding_energy(
    complex_structure: Structure,
    separated_structure: Structure,
    interface: InterfaceRegion,
    minimize: bool = True,
) -> float:
    """
    Calculate binding energy (dG_separated).

    dG = E_complex - E_separated

    Args:
        complex_structure: Bound complex
        separated_structure: Separated chains
        interface: Interface region (for reference)
        minimize: If True, minimize structures before energy calculation

    Returns:
        Binding energy in REU (negative = favorable)
    """
    # Prepare complex
    complex_modeller, forcefield = prepare_structure_for_openmm(complex_structure)
    energy_complex = calculate_potential_energy(complex_modeller, forcefield, minimize)

    # Prepare separated
    sep_modeller, forcefield = prepare_structure_for_openmm(separated_structure)
    energy_separated = calculate_potential_energy(sep_modeller, forcefield, minimize)

    # Calculate delta (complex - separated)
    dG_kj = energy_complex - energy_separated

    # Convert to REU
    return dG_kj * KJ_TO_REU


def check_openmm_available() -> bool:
    """Check if OpenMM is available for energy calculations."""
    return _check_openmm_available()


def _write_temp_pdb(structure: Structure, path: Path) -> None:
    """Write structure to temporary PDB file.

    Only writes standard amino acid residues to ensure compatibility
    with AMBER force field. Glycans, ligands, and other non-standard
    residues are excluded.
    """
    with open(path, "w") as f:
        atom_num = 1
        for chain in structure.chains.values():
            last_protein_residue = None
            has_protein = False

            for residue in chain.residues:
                # Skip non-standard residues (glycans, ligands, etc.)
                if residue.name not in STANDARD_RESIDUES:
                    continue

                has_protein = True
                for atom in residue.atoms:
                    # Handle atom name formatting (left-justify if 4 chars)
                    atom_name = atom.name
                    if len(atom_name) < 4:
                        atom_name = f" {atom_name:<3s}"
                    else:
                        atom_name = f"{atom_name:<4s}"

                    # PDB format columns:
                    # 1-6: record, 7-11: serial, 12: space, 13-16: atom name,
                    # 17: altLoc, 18-20: resName, 21: space, 22: chain,
                    # 23-26: resSeq, 27: iCode, 28-30: spaces, 31-54: coords,
                    # 55-60: occupancy, 61-66: B-factor, 77-78: element
                    f.write(
                        f"ATOM  {atom_num:5d} {atom_name} "
                        f"{residue.name:3s} {chain.chain_id:1s}"
                        f"{residue.seq_num:4d}{residue.insertion_code:1s}   "
                        f"{atom.coords[0]:8.3f}{atom.coords[1]:8.3f}{atom.coords[2]:8.3f}"
                        f"{atom.occupancy:6.2f}{atom.b_factor:6.2f}          "
                        f"{atom.element:>2s}\n"
                    )
                    atom_num += 1
                last_protein_residue = residue

            # Write TER record after protein portion of each chain
            if has_protein and last_protein_residue is not None:
                f.write(
                    f"TER   {atom_num:5d}      "
                    f"{last_protein_residue.name:3s} {chain.chain_id:1s}"
                    f"{last_protein_residue.seq_num:4d}{last_protein_residue.insertion_code:1s}\n"
                )
                atom_num += 1
        f.write("END\n")
