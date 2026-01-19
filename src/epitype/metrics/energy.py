"""Binding energy calculations using OpenMM with AMBER ff14SB."""
from __future__ import annotations

import tempfile
from dataclasses import dataclass
from enum import Enum
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


def _get_best_platform() -> tuple:
    """Select the best available OpenMM platform with fallback.

    Tries platforms in order of performance: CUDA > OpenCL > CPU > Reference.
    GPU platforms (CUDA/OpenCL) use mixed precision for best speed/accuracy balance.

    Returns:
        Tuple of (Platform, properties_dict)
    """
    from openmm import Platform

    # Platform preference order and their properties
    platforms = [
        ("CUDA", {"Precision": "mixed"}),
        ("OpenCL", {"Precision": "mixed"}),
        ("CPU", {}),
        ("Reference", {}),
    ]

    for name, properties in platforms:
        try:
            platform = Platform.getPlatformByName(name)
            return platform, properties
        except Exception:
            continue

    # Should never reach here (Reference is always available)
    raise RuntimeError("No OpenMM platform available")


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


class SolventModel(str, Enum):
    """Implicit solvent model for energy calculations.

    Models:
        VACUUM: No solvent (original behavior, fastest)
        OBC2: Onufriev-Bashford-Case II - recommended for protein-protein binding
        GBN: GB with Neck corrections - can be more accurate for some systems
        HCT: Hawkins-Cramer-Truhlar - fastest implicit solvent model
    """

    VACUUM = "vacuum"
    OBC2 = "obc2"
    GBN = "gbn"
    HCT = "hct"


@dataclass
class SolventConfig:
    """Configuration for implicit solvent calculations.

    Attributes:
        model: Implicit solvent model to use
        salt_concentration: Ionic strength in Molar (default 0.15 M for PBS-like conditions)
        solvent_dielectric: Dielectric constant of solvent (default 80.0 for water at 25C)
        solute_dielectric: Dielectric constant of solute interior (default 1.0)
    """

    model: SolventModel = SolventModel.OBC2
    salt_concentration: float = 0.15
    solvent_dielectric: float = 80.0
    solute_dielectric: float = 1.0


@dataclass
class BindingEnergyResult:
    """Result from binding energy calculation.

    Attributes:
        dG_separated: Binding energy in REU (negative = favorable)
        energy_complex: Potential energy of bound complex in kJ/mol
        energy_separated: Potential energy of separated chains in kJ/mol
        solvent_model: Name of solvent model used
        salt_concentration: Salt concentration used in Molar
    """

    dG_separated: float
    energy_complex: float
    energy_separated: float
    solvent_model: str
    salt_concentration: float


def prepare_structure_for_openmm(
    structure: Structure,
    solvent_model: SolventModel = SolventModel.VACUUM,
) -> tuple:
    """
    Prepare structure for OpenMM energy calculation.

    Includes fixing missing atoms/residues via PDBFixer.

    Args:
        structure: Input structure
        solvent_model: Implicit solvent model to use (affects force field selection)

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

    # Create force field with appropriate solvent model
    if solvent_model == SolventModel.VACUUM:
        # Vacuum: just protein force field (no explicit water model needed)
        forcefield = ForceField("amber14-all.xml")
    else:
        # Implicit solvent: protein force field + GB model
        solvent_xml_map = {
            SolventModel.OBC2: "implicit/obc2.xml",
            SolventModel.GBN: "implicit/gbn.xml",
            SolventModel.HCT: "implicit/hct.xml",
        }
        forcefield = ForceField("amber14-all.xml", solvent_xml_map[solvent_model])

    # Create modeller (no explicit solvent for interface calculations)
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
        forcefield: OpenMM ForceField (may include implicit solvent parameters)
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

    # Create system (implicit solvent is handled by the force field if loaded)
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.NoCutoff,
        constraints=None,
        rigidWater=False,
    )

    # Create integrator (not used for single-point energy)
    integrator = LangevinMiddleIntegrator(
        300 * unit.kelvin,
        1.0 / unit.picosecond,
        0.002 * unit.picoseconds,
    )

    # Create simulation with best available platform (GPU if available)
    platform, platform_properties = _get_best_platform()
    simulation = app.Simulation(
        modeller.topology, system, integrator, platform, platform_properties
    )
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
    solvent_config: SolventConfig | None = None,
) -> BindingEnergyResult:
    """
    Calculate binding energy (dG_separated).

    dG = E_complex - E_separated

    Args:
        complex_structure: Bound complex
        separated_structure: Separated chains
        interface: Interface region (for reference)
        minimize: If True, minimize structures before energy calculation
        solvent_config: Implicit solvent configuration (None uses vacuum)

    Returns:
        BindingEnergyResult with binding energy and component energies
    """
    # Use vacuum if no solvent config provided
    solvent_config = solvent_config or SolventConfig()

    # Prepare complex (force field includes solvent model if specified)
    complex_modeller, forcefield = prepare_structure_for_openmm(
        complex_structure, solvent_config.model
    )
    energy_complex = calculate_potential_energy(
        complex_modeller, forcefield, minimize
    )

    # Prepare separated (use same solvent model)
    sep_modeller, forcefield = prepare_structure_for_openmm(
        separated_structure, solvent_config.model
    )
    energy_separated = calculate_potential_energy(
        sep_modeller, forcefield, minimize
    )

    # Calculate delta (complex - separated)
    dG_kj = energy_complex - energy_separated

    # Return result with all details
    return BindingEnergyResult(
        dG_separated=dG_kj * KJ_TO_REU,
        energy_complex=energy_complex,
        energy_separated=energy_separated,
        solvent_model=solvent_config.model.value,
        salt_concentration=solvent_config.salt_concentration,
    )


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
