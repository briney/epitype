"""Unit tests for SASA calculations."""
from pathlib import Path

from epitype.core.interface import InterfaceRegion, detect_interface
from epitype.core.separation import separate_chains
from epitype.core.structure import Structure
from epitype.io.parsers import parse_structure
from epitype.metrics.sasa import (
    calculate_dsasa,
    calculate_interface_sasa,
    calculate_sasa,
)


class TestCalculateSASA:
    """Tests for calculate_sasa function."""

    def test_sasa_is_positive(self, minimal_structure: Structure):
        """Test that SASA is non-negative."""
        result = calculate_sasa(minimal_structure)

        assert result.total >= 0
        assert result.polar >= 0
        assert result.apolar >= 0

    def test_polar_plus_apolar_equals_total(self, minimal_structure: Structure):
        """Test that polar + apolar approximately equals total."""
        result = calculate_sasa(minimal_structure)

        # Allow small numerical tolerance
        assert abs(result.polar + result.apolar - result.total) < 0.1

    def test_per_atom_sasa_computed(self, minimal_structure: Structure):
        """Test that per-atom SASA values are computed."""
        result = calculate_sasa(minimal_structure)

        assert len(result.per_atom) > 0
        assert all(v >= 0 for v in result.per_atom.values())

    def test_per_residue_sasa_computed(self, minimal_structure: Structure):
        """Test that per-residue SASA values are computed."""
        result = calculate_sasa(minimal_structure)

        assert len(result.per_residue) > 0
        assert all(v >= 0 for v in result.per_residue.values())

    def test_sasa_real_structure(self, pdb_1yy9: Path):
        """Test SASA calculation on real PDB structure."""
        structure = parse_structure(pdb_1yy9)
        result = calculate_sasa(structure)

        # Real protein should have significant SASA
        assert result.total > 1000  # Ų

    def test_empty_structure_returns_zero(self):
        """Test that empty structure returns zero SASA."""
        structure = Structure(name="empty")
        result = calculate_sasa(structure)

        assert result.total == 0.0
        assert result.polar == 0.0
        assert result.apolar == 0.0

    def test_lee_richards_algorithm(self, minimal_structure: Structure):
        """Test using Lee-Richards algorithm."""
        result = calculate_sasa(minimal_structure, algorithm="LeeRichards")

        assert result.total > 0

    def test_custom_probe_radius(self, minimal_structure: Structure):
        """Test using custom probe radius."""
        result_small = calculate_sasa(minimal_structure, probe_radius=1.0)
        result_large = calculate_sasa(minimal_structure, probe_radius=2.0)

        # Larger probe should generally find less accessible surface
        # (can get into fewer crevices)
        assert result_small.total != result_large.total


class TestCalculateDSASA:
    """Tests for calculate_dsasa function."""

    def test_dsasa_is_positive_for_interface(self, pdb_1yy9: Path):
        """Test that dSASA is positive for complexed interface."""
        structure = parse_structure(pdb_1yy9)

        # Use H,L chains as antibody vs A chain as antigen
        interface = detect_interface(structure, ["C", "D"], ["A"], cutoff=8.0)
        separated = separate_chains(structure, interface)

        dsasa = calculate_dsasa(structure, separated, interface)

        # dSASA should be positive (buried surface area)
        assert dsasa > 0

    def test_dsasa_reasonable_range(self, pdb_1yy9: Path):
        """Test that dSASA is in reasonable range for antibody-antigen."""
        structure = parse_structure(pdb_1yy9)
        interface = detect_interface(structure, ["C", "D"], ["A"], cutoff=8.0)
        separated = separate_chains(structure, interface)

        dsasa = calculate_dsasa(structure, separated, interface)

        # Antibody-antigen interfaces typically 1500-2500 Ų
        assert 500 < dsasa < 5000


class TestCalculateInterfaceSASA:
    """Tests for calculate_interface_sasa function."""

    def test_interface_sasa_both_groups(self, pdb_1yy9: Path):
        """Test that both group SASAs are computed."""
        structure = parse_structure(pdb_1yy9)
        interface = detect_interface(structure, ["C", "D"], ["A"], cutoff=8.0)

        g1_sasa, g2_sasa = calculate_interface_sasa(structure, interface)

        # Both groups should have positive SASA
        assert g1_sasa >= 0
        assert g2_sasa >= 0

    def test_interface_sasa_minimal(self, minimal_structure: Structure, interface_ab: InterfaceRegion):
        """Test interface SASA on minimal structure."""
        g1_sasa, g2_sasa = calculate_interface_sasa(minimal_structure, interface_ab)

        assert g1_sasa >= 0
        assert g2_sasa >= 0
