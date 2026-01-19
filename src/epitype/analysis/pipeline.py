"""Full interface analysis pipeline."""
from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path

from epitype.core.interface import detect_interface, parse_chain_groups
from epitype.core.separation import separate_chains
from epitype.core.structure import Structure
from epitype.core.types import InterfaceMetrics
from epitype.io.parsers import parse_structure
from epitype.metrics.energy import SolventModel
from epitype.metrics.hbonds import analyze_hbonds
from epitype.metrics.packstat import calculate_packstat
from epitype.metrics.sasa import calculate_sasa


@dataclass
class AnalysisConfig:
    """Configuration for interface analysis."""

    chain_spec: str  # e.g., "HL_A"
    interface_cutoff: float = 8.0
    probe_radius: float = 1.4
    separation_distance: float = 500.0
    minimize_energy: bool = True
    compute_energy: bool = True
    compute_shape: bool = True
    compute_packstat: bool = True
    # Solvent configuration for energy calculations
    solvent_model: SolventModel = field(default=SolventModel.OBC2)
    salt_concentration: float = 0.15  # Molar (PBS-like conditions)


ProgressCallback = Callable[[str, float], None]


def analyze_interface(
    structure_path: str | Path,
    config: AnalysisConfig,
    progress_callback: ProgressCallback | None = None,
) -> InterfaceMetrics:
    """
    Run complete interface analysis pipeline.

    Args:
        structure_path: Path to PDB/mmCIF file
        config: Analysis configuration
        progress_callback: Optional callback for progress updates (step_name, fraction)

    Returns:
        InterfaceMetrics with all computed values
    """
    def report(step: str, progress: float) -> None:
        if progress_callback:
            progress_callback(step, progress)

    # Step 1: Parse structure
    report("Parsing structure", 0.0)
    structure = parse_structure(structure_path)

    # Step 2: Detect interface
    report("Detecting interface", 0.1)
    group1, group2 = parse_chain_groups(config.chain_spec)
    interface = detect_interface(
        structure,
        group1,
        group2,
        cutoff=config.interface_cutoff,
    )

    if not interface.residues_group1 or not interface.residues_group2:
        raise ValueError(
            f"No interface detected between {group1} and {group2}. "
            "Check chain specification and cutoff distance."
        )

    # Step 3: Create separated structure
    report("Separating chains", 0.2)
    separated = separate_chains(
        structure,
        interface,
        distance=config.separation_distance,
    )

    # Step 4: Calculate SASA
    report("Calculating SASA", 0.3)
    sasa_complex = calculate_sasa(structure, config.probe_radius)
    sasa_separated = calculate_sasa(separated, config.probe_radius)
    dsasa = sasa_separated.total - sasa_complex.total

    # Step 5: Calculate binding energy (optional, requires OpenMM)
    report("Calculating binding energy", 0.4)
    dG = 0.0
    if config.compute_energy:
        try:
            from epitype.metrics.energy import (
                SolventConfig,
                calculate_binding_energy,
                check_openmm_available,
            )
            if check_openmm_available():
                solvent_config = SolventConfig(
                    model=config.solvent_model,
                    salt_concentration=config.salt_concentration,
                )
                energy_result = calculate_binding_energy(
                    structure,
                    separated,
                    interface,
                    minimize=config.minimize_energy,
                    solvent_config=solvent_config,
                )
                dG = energy_result.dG_separated
        except ImportError:
            pass  # OpenMM not available

    # Step 6: Analyze hydrogen bonds
    report("Analyzing hydrogen bonds", 0.6)
    hbond_result = analyze_hbonds(structure, interface)

    # Step 7: Calculate shape complementarity (optional, requires NanoShaper)
    report("Calculating shape complementarity", 0.7)
    sc_value = 0.0
    if config.compute_shape:
        try:
            from epitype.metrics.shape import calculate_shape_complementarity
            from epitype.surface.nanoshaper import check_nanoshaper_available
            if check_nanoshaper_available():
                sc_result = calculate_shape_complementarity(structure, interface)
                sc_value = sc_result.sc_value
        except (ImportError, RuntimeError):
            pass  # NanoShaper not available

    # Step 8: Calculate PackStat
    report("Calculating PackStat", 0.85)
    packstat = 0.0
    if config.compute_packstat:
        packstat_result = calculate_packstat(structure, interface)
        packstat = packstat_result.packstat_interface

    report("Complete", 1.0)

    # Compute derived metrics
    dG_per_dSASA = dG / dsasa if dsasa > 0 else 0.0

    return InterfaceMetrics(
        dG_separated=dG,
        dSASA_int=dsasa,
        dG_per_dSASA=dG_per_dSASA,
        sc_value=sc_value,
        packstat=packstat,
        delta_unsat_hbonds=hbond_result.buried_unsatisfied_donors + hbond_result.buried_unsatisfied_acceptors,
        hbonds_int=hbond_result.interface_hbonds,
        interface_residues_1=interface.residue_ids_group1,
        interface_residues_2=interface.residue_ids_group2,
        sasa_complex=sasa_complex.total,
        sasa_separated=sasa_separated.total,
    )


def analyze_structure(
    structure: Structure,
    chain_spec: str,
    **kwargs,
) -> InterfaceMetrics:
    """
    Analyze a pre-loaded Structure object.

    Convenience function for API use.

    Args:
        structure: Loaded Structure object
        chain_spec: Chain specification (e.g., "HL_A")
        **kwargs: Additional config options

    Returns:
        InterfaceMetrics
    """
    config = AnalysisConfig(chain_spec=chain_spec, **kwargs)

    group1, group2 = parse_chain_groups(config.chain_spec)
    interface = detect_interface(structure, group1, group2, cutoff=config.interface_cutoff)

    if not interface.residues_group1 or not interface.residues_group2:
        raise ValueError(
            f"No interface detected between {group1} and {group2}. "
            "Check chain specification and cutoff distance."
        )

    separated = separate_chains(structure, interface, distance=config.separation_distance)

    sasa_complex = calculate_sasa(structure, config.probe_radius)
    sasa_separated = calculate_sasa(separated, config.probe_radius)
    dsasa = sasa_separated.total - sasa_complex.total

    dG = 0.0
    if config.compute_energy:
        try:
            from epitype.metrics.energy import (
                SolventConfig,
                calculate_binding_energy,
                check_openmm_available,
            )
            if check_openmm_available():
                solvent_config = SolventConfig(
                    model=config.solvent_model,
                    salt_concentration=config.salt_concentration,
                )
                energy_result = calculate_binding_energy(
                    structure, separated, interface,
                    minimize=config.minimize_energy,
                    solvent_config=solvent_config,
                )
                dG = energy_result.dG_separated
        except ImportError:
            pass

    hbond_result = analyze_hbonds(structure, interface)

    sc_value = 0.0
    if config.compute_shape:
        try:
            from epitype.metrics.shape import calculate_shape_complementarity
            from epitype.surface.nanoshaper import check_nanoshaper_available
            if check_nanoshaper_available():
                sc_result = calculate_shape_complementarity(structure, interface)
                sc_value = sc_result.sc_value
        except (ImportError, RuntimeError):
            pass

    packstat = 0.0
    if config.compute_packstat:
        packstat_result = calculate_packstat(structure, interface)
        packstat = packstat_result.packstat_interface

    dG_per_dSASA = dG / dsasa if dsasa > 0 else 0.0

    return InterfaceMetrics(
        dG_separated=dG,
        dSASA_int=dsasa,
        dG_per_dSASA=dG_per_dSASA,
        sc_value=sc_value,
        packstat=packstat,
        delta_unsat_hbonds=hbond_result.buried_unsatisfied_donors + hbond_result.buried_unsatisfied_acceptors,
        hbonds_int=hbond_result.interface_hbonds,
        interface_residues_1=interface.residue_ids_group1,
        interface_residues_2=interface.residue_ids_group2,
        sasa_complex=sasa_complex.total,
        sasa_separated=sasa_separated.total,
    )


def analyze(
    structure_path: str | Path,
    chains: str,
    **kwargs,
) -> InterfaceMetrics:
    """
    Convenience function for simple analysis.

    Args:
        structure_path: Path to PDB/mmCIF file
        chains: Chain specification (e.g., "HL_A")
        **kwargs: Additional config options

    Returns:
        InterfaceMetrics
    """
    config = AnalysisConfig(chain_spec=chains, **kwargs)
    return analyze_interface(structure_path, config)
