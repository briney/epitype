"""
Individual metric calculations.

Each metric can be computed independently:

    from epitype import load_structure, detect_interface
    from epitype.metrics import sasa, hbonds, shape, packstat, energy

    structure = load_structure("complex.pdb")
    interface = detect_interface(structure, ["H", "L"], ["A"])

    # SASA
    sasa_result = sasa.calculate_sasa(structure)
    dsasa = sasa.calculate_dsasa(complex_struct, separated_struct, interface)

    # Hydrogen bonds
    hbond_result = hbonds.analyze_hbonds(structure, interface)

    # Shape complementarity
    sc_result = shape.calculate_shape_complementarity(structure, interface)

    # PackStat
    packstat_result = packstat.calculate_packstat(structure, interface)
"""

# Import submodules for direct access
from epitype.metrics import energy, hbonds, packstat, sasa, shape
from epitype.metrics.energy import (
    KJ_TO_REU,
    calculate_binding_energy,
    calculate_potential_energy,
    check_openmm_available,
    prepare_structure_for_openmm,
)
from epitype.metrics.hbonds import (
    HBondResult,
    HydrogenBond,
    analyze_hbonds,
    count_interface_hbonds,
    count_unsatisfied_hbonds,
    detect_hydrogen_bonds,
)
from epitype.metrics.packstat import (
    PackStatResult,
    calculate_packstat,
)

# Import commonly used classes and functions
from epitype.metrics.sasa import (
    SASAResult,
    calculate_dsasa,
    calculate_interface_sasa,
    calculate_sasa,
)
from epitype.metrics.shape import (
    ShapeComplementarityResult,
    calculate_shape_complementarity,
)

__all__ = [
    # Submodules
    "sasa",
    "energy",
    "hbonds",
    "shape",
    "packstat",
    # SASA
    "SASAResult",
    "calculate_sasa",
    "calculate_dsasa",
    "calculate_interface_sasa",
    # H-bonds
    "HydrogenBond",
    "HBondResult",
    "detect_hydrogen_bonds",
    "count_interface_hbonds",
    "count_unsatisfied_hbonds",
    "analyze_hbonds",
    # PackStat
    "PackStatResult",
    "calculate_packstat",
    # Shape complementarity
    "ShapeComplementarityResult",
    "calculate_shape_complementarity",
    # Energy
    "KJ_TO_REU",
    "prepare_structure_for_openmm",
    "calculate_potential_energy",
    "calculate_binding_energy",
    "check_openmm_available",
]
