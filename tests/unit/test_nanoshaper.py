"""Unit tests for NanoShaper surface generation."""
from pathlib import Path

import numpy as np
import pytest

from epitype.core.structure import Structure
from epitype.io.parsers import parse_structure
from epitype.surface.nanoshaper import (
    _write_xyzr,
    check_nanoshaper_available,
    generate_surface,
)
from epitype.surface.types import SurfacePoint

# Skip all tests if NanoShaper is not available
pytestmark = pytest.mark.skipif(
    not check_nanoshaper_available(),
    reason="NanoShaper not available"
)


class TestCheckNanoShaperAvailable:
    """Tests for NanoShaper availability check."""

    def test_check_returns_bool(self):
        """Test that check_nanoshaper_available returns a boolean."""
        result = check_nanoshaper_available()
        assert isinstance(result, bool)

    def test_check_with_invalid_path(self):
        """Test check with invalid path returns False."""
        result = check_nanoshaper_available(nanoshaper_path="/nonexistent/NanoShaper")
        assert result is False


class TestWriteXYZR:
    """Tests for XYZR file writing."""

    def test_write_xyzr_creates_file(self, minimal_structure: Structure, temp_dir: Path):
        """Test that XYZR file is created."""
        xyzr_path = temp_dir / "test.xyzr"
        _write_xyzr(minimal_structure, xyzr_path)

        assert xyzr_path.exists()

    def test_write_xyzr_correct_format(self, minimal_structure: Structure, temp_dir: Path):
        """Test that XYZR file has correct format."""
        xyzr_path = temp_dir / "test.xyzr"
        _write_xyzr(minimal_structure, xyzr_path)

        with open(xyzr_path) as f:
            lines = f.readlines()

        # Should have one line per atom
        assert len(lines) == minimal_structure.num_atoms

        # Each line should have 4 values (x, y, z, radius)
        for line in lines:
            parts = line.strip().split()
            assert len(parts) == 4
            # All values should be parseable as floats
            for part in parts:
                float(part)


class TestGenerateSurface:
    """Tests for surface generation."""

    def test_surface_generation_returns_list(self, minimal_structure: Structure):
        """Test that generate_surface returns a list."""
        surface = generate_surface(minimal_structure)
        assert isinstance(surface, list)

    def test_surface_points_have_correct_attributes(self, minimal_structure: Structure):
        """Test that surface points have correct attributes."""
        surface = generate_surface(minimal_structure)

        if surface:  # May be empty if structure is too small
            point = surface[0]
            assert isinstance(point, SurfacePoint)
            assert point.coords.shape == (3,)
            assert point.normal.shape == (3,)
            assert isinstance(point.atom_index, int)
            assert isinstance(point.area, float)

    def test_surface_normals_are_unit_vectors(self, minimal_structure: Structure):
        """Test that surface normals are approximately unit vectors."""
        surface = generate_surface(minimal_structure)

        for point in surface:
            norm = np.linalg.norm(point.normal)
            # NanoShaper normals should be close to unit length
            assert 0.9 < norm < 1.1

    def test_surface_on_real_structure(self, pdb_1yy9: Path):
        """Test surface generation on real PDB structure."""
        structure = parse_structure(pdb_1yy9)
        surface = generate_surface(structure)

        # Real protein should have many surface points
        assert len(surface) > 100

    def test_custom_probe_radius(self, minimal_structure: Structure):
        """Test surface generation with custom probe radius."""
        surface_small = generate_surface(minimal_structure, probe_radius=1.0)
        surface_large = generate_surface(minimal_structure, probe_radius=2.0)

        # Different probe radii should give different results
        # (though both may be empty for minimal structure)
        assert isinstance(surface_small, list)
        assert isinstance(surface_large, list)

    def test_custom_grid_scale(self, pdb_1yy9: Path):
        """Test surface generation with custom grid scale."""
        structure = parse_structure(pdb_1yy9)

        surface_low = generate_surface(structure, grid_scale=1.0)
        surface_high = generate_surface(structure, grid_scale=3.0)

        # Higher grid scale should have more surface points
        assert len(surface_high) > len(surface_low)

    def test_surface_invalid_nanoshaper_path(self, minimal_structure: Structure):
        """Test that invalid NanoShaper path raises RuntimeError."""
        with pytest.raises(RuntimeError, match="NanoShaper executable not found"):
            generate_surface(minimal_structure, nanoshaper_path="/nonexistent/NanoShaper")

    def test_empty_structure_returns_empty_list(self):
        """Test that empty structure returns empty surface."""
        structure = Structure(name="empty")
        surface = generate_surface(structure)

        assert surface == []
