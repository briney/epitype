"""Unit tests for PackStat calculations."""
from pathlib import Path

import numpy as np

from epitype.core.interface import InterfaceRegion, detect_interface
from epitype.core.structure import Structure
from epitype.io.parsers import parse_structure
from epitype.metrics.packstat import (
    _reference_sasa,
    calculate_packstat,
)


class TestCalculatePackstat:
    """Tests for calculate_packstat function."""

    def test_packstat_returns_0_to_1(self, minimal_structure: Structure):
        """Test that packstat value is in [0, 1] range."""
        result = calculate_packstat(minimal_structure)

        assert 0.0 <= result.packstat <= 1.0

    def test_per_probe_scores_computed(self, minimal_structure: Structure):
        """Test that per-probe scores are computed."""
        result = calculate_packstat(minimal_structure)

        # Should have 30 probe scores
        assert len(result.per_probe_scores) == 30

        # All scores in valid range
        for score in result.per_probe_scores.values():
            assert 0.0 <= score <= 1.0

    def test_packstat_interface_computed(self, minimal_structure: Structure, interface_ab: InterfaceRegion):
        """Test that interface-specific packstat is computed."""
        result = calculate_packstat(minimal_structure, interface=interface_ab)

        assert 0.0 <= result.packstat_interface <= 1.0

    def test_cavity_volume_non_negative(self, minimal_structure: Structure):
        """Test that cavity volume is non-negative."""
        result = calculate_packstat(minimal_structure)

        assert result.cavity_volume >= 0.0

    def test_packstat_real_structure(self, pdb_1yy9: Path):
        """Test packstat on real structure."""
        structure = parse_structure(pdb_1yy9)
        result = calculate_packstat(structure)

        # Real protein should have reasonable packing
        assert 0.3 <= result.packstat <= 1.0

    def test_empty_structure_returns_zero(self):
        """Test that empty structure returns zero packstat."""
        structure = Structure(name="empty")
        result = calculate_packstat(structure)

        assert result.packstat == 0.0
        assert result.packstat_interface == 0.0
        assert result.cavity_volume == 0.0
        assert result.per_probe_scores == {}


class TestPackstatHelpers:
    """Tests for helper functions."""

    def test_reference_sasa_positive(self):
        """Test that reference SASA is positive."""
        radii = np.array([1.7, 1.55, 1.52])
        ref_sasa = _reference_sasa(radii, probe_radius=1.4)

        assert ref_sasa > 0

    def test_reference_sasa_increases_with_probe(self):
        """Test that reference SASA increases with probe radius."""
        radii = np.array([1.7, 1.55, 1.52])

        ref_small = _reference_sasa(radii, probe_radius=1.0)
        ref_large = _reference_sasa(radii, probe_radius=2.0)

        assert ref_large > ref_small


class TestPackstatInterface:
    """Tests for interface-specific packstat."""

    def test_interface_packstat_computed(self, pdb_1yy9: Path):
        """Test interface packstat on real structure."""
        structure = parse_structure(pdb_1yy9)
        interface = detect_interface(structure, ["C", "D"], ["A"], cutoff=8.0)

        result = calculate_packstat(structure, interface=interface)

        assert 0.0 <= result.packstat_interface <= 1.0

    def test_interface_packstat_reasonable_range(self, pdb_1yy9: Path):
        """Test that interface packstat is in reasonable range."""
        structure = parse_structure(pdb_1yy9)
        interface = detect_interface(structure, ["C", "D"], ["A"], cutoff=8.0)

        result = calculate_packstat(structure, interface=interface)

        # Interface should have decent packing for antibody-antigen
        assert 0.2 <= result.packstat_interface <= 1.0


class TestPackstatPerformance:
    """Tests to verify performance improvements."""

    def test_packstat_completes_quickly(self, pdb_1yy9: Path):
        """Test that packstat calculation completes in reasonable time."""
        import time

        structure = parse_structure(pdb_1yy9)

        start = time.time()
        result = calculate_packstat(structure)
        elapsed = time.time() - start

        # Should complete in under 30 seconds (was 12+ minutes before)
        assert elapsed < 30.0
        assert result.packstat > 0
