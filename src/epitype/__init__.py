"""
Interface Analyzer: Protein-protein interface analysis toolkit.

Example usage:

    from epitype import analyze, load_structure

    # Full analysis
    metrics = analyze("complex.pdb", chains="HL_A")
    print(f"Binding energy: {metrics.dG_separated:.2f} REU")
    print(f"Buried surface: {metrics.dSASA_int:.1f} Angstrom^2")

    # Individual metrics
    from epitype.metrics import sasa, hbonds

    structure = load_structure("complex.pdb")
    sasa_result = sasa.calculate_sasa(structure)
"""

__version__ = "0.1.0"

# Core data structures
# Metric modules (for individual access)
from epitype import metrics
from epitype.analysis.pipeline import (
    AnalysisConfig,
    analyze,
    analyze_interface,
    analyze_structure,
)
from epitype.core.interface import InterfaceRegion, detect_interface, parse_chain_groups
from epitype.core.structure import Atom, Chain, Residue, Structure
from epitype.core.types import InterfaceMetrics

# High-level API
from epitype.io.parsers import parse_structure as load_structure

__all__ = [
    # Version
    "__version__",
    # Data structures
    "Structure",
    "Chain",
    "Residue",
    "Atom",
    "InterfaceRegion",
    "InterfaceMetrics",
    # Functions
    "load_structure",
    "analyze",
    "analyze_interface",
    "analyze_structure",
    "detect_interface",
    "parse_chain_groups",
    # Config
    "AnalysisConfig",
    # Submodules
    "metrics",
]
