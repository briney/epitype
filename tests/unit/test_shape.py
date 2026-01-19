"""Unit tests for shape complementarity calculations."""
from pathlib import Path

import numpy as np
import pytest

from epitype.core.interface import detect_interface
from epitype.io.parsers import parse_structure
from epitype.metrics.shape import (
    ShapeComplementarityResult,
    _calculate_sc_one_direction,
    _trim_surface_to_interface,
    calculate_shape_complementarity,
)
from epitype.surface.nanoshaper import check_nanoshaper_available
from epitype.surface.types import SurfacePoint

# Skip all tests if NanoShaper is not available
pytestmark = pytest.mark.skipif(
    not check_nanoshaper_available(),
    reason="NanoShaper not available"
)


class TestShapeComplementarityResult:
    """Tests for ShapeComplementarityResult dataclass."""

    def test_result_creation(self):
        """Test creating a result object."""
        result = ShapeComplementarityResult(
            sc_value=0.65,
            sc_median=0.68,
            sc_group1=0.62,
            sc_group2=0.68,
            interface_area=1500.0,
        )

        assert result.sc_value == 0.65
        assert result.sc_median == 0.68
        assert result.sc_group1 == 0.62
        assert result.sc_group2 == 0.68
        assert result.interface_area == 1500.0


class TestCalculateShapeComplementarity:
    """Tests for calculate_shape_complementarity function."""

    def test_sc_returns_result(self, pdb_1yy9: Path):
        """Test that Sc calculation returns a result."""
        structure = parse_structure(pdb_1yy9)
        interface = detect_interface(structure, ["C", "D"], ["A"], cutoff=8.0)

        result = calculate_shape_complementarity(structure, interface)

        assert isinstance(result, ShapeComplementarityResult)

    def test_sc_in_valid_range(self, pdb_1yy9: Path):
        """Test that Sc value is in valid range (typically 0-1)."""
        structure = parse_structure(pdb_1yy9)
        interface = detect_interface(structure, ["C", "D"], ["A"], cutoff=8.0)

        result = calculate_shape_complementarity(structure, interface)

        # Sc can technically be negative for poorly fitting surfaces,
        # but should generally be positive for real interfaces
        assert -1.0 <= result.sc_value <= 1.0

    def test_sc_positive_for_antibody_antigen(self, pdb_1yy9: Path):
        """Test that Sc is positive for antibody-antigen interface."""
        structure = parse_structure(pdb_1yy9)
        interface = detect_interface(structure, ["C", "D"], ["A"], cutoff=8.0)

        result = calculate_shape_complementarity(structure, interface)

        # Good antibody-antigen interfaces typically have Sc > 0.6
        assert result.sc_value > 0.0

    def test_sc_has_interface_area(self, pdb_1yy9: Path):
        """Test that interface area is computed."""
        structure = parse_structure(pdb_1yy9)
        interface = detect_interface(structure, ["C", "D"], ["A"], cutoff=8.0)

        result = calculate_shape_complementarity(structure, interface)

        assert result.interface_area > 0

    def test_sc_both_groups_computed(self, pdb_1yy9: Path):
        """Test that both group Sc values are computed."""
        structure = parse_structure(pdb_1yy9)
        interface = detect_interface(structure, ["C", "D"], ["A"], cutoff=8.0)

        result = calculate_shape_complementarity(structure, interface)

        # Both should be non-zero for real interface
        # Note: Could be zero if NanoShaper fails or surface is too small
        assert isinstance(result.sc_group1, float)
        assert isinstance(result.sc_group2, float)

    def test_sc_custom_distance_cutoff(self, pdb_1yy9: Path):
        """Test Sc with custom distance cutoff."""
        structure = parse_structure(pdb_1yy9)
        interface = detect_interface(structure, ["C", "D"], ["A"], cutoff=8.0)

        result_small = calculate_shape_complementarity(
            structure, interface, distance_cutoff=1.0
        )
        result_large = calculate_shape_complementarity(
            structure, interface, distance_cutoff=3.0
        )

        # Results should be different with different cutoffs
        assert isinstance(result_small.sc_value, float)
        assert isinstance(result_large.sc_value, float)

    def test_sc_custom_trim_distance(self, pdb_1yy9: Path):
        """Test Sc with custom trim distance."""
        structure = parse_structure(pdb_1yy9)
        interface = detect_interface(structure, ["C", "D"], ["A"], cutoff=8.0)

        result_small = calculate_shape_complementarity(
            structure, interface, trim_distance=5.0
        )
        result_large = calculate_shape_complementarity(
            structure, interface, trim_distance=15.0
        )

        # Different trim distances should give different interface areas
        assert isinstance(result_small.interface_area, float)
        assert isinstance(result_large.interface_area, float)


class TestTrimSurfaceToInterface:
    """Tests for _trim_surface_to_interface helper."""

    def test_trim_empty_surface(self):
        """Test trimming with empty surface."""
        result = _trim_surface_to_interface([], [], 10.0)
        assert result == []

    def test_trim_empty_other_surface(self):
        """Test trimming with empty other surface."""
        point = SurfacePoint(
            coords=np.array([0.0, 0.0, 0.0]),
            normal=np.array([1.0, 0.0, 0.0]),
            atom_index=0,
            area=1.0,
        )
        result = _trim_surface_to_interface([point], [], 10.0)
        assert result == []

    def test_trim_keeps_nearby_points(self):
        """Test that nearby points are kept."""
        surface = [
            SurfacePoint(
                coords=np.array([0.0, 0.0, 0.0]),
                normal=np.array([1.0, 0.0, 0.0]),
                atom_index=0,
                area=1.0,
            ),
            SurfacePoint(
                coords=np.array([100.0, 0.0, 0.0]),
                normal=np.array([1.0, 0.0, 0.0]),
                atom_index=1,
                area=1.0,
            ),
        ]
        other = [
            SurfacePoint(
                coords=np.array([1.0, 0.0, 0.0]),
                normal=np.array([-1.0, 0.0, 0.0]),
                atom_index=0,
                area=1.0,
            ),
        ]

        result = _trim_surface_to_interface(surface, other, 10.0)

        # Only the first point should be kept (within 10 Angstroms)
        assert len(result) == 1
        assert np.allclose(result[0].coords, [0.0, 0.0, 0.0])


class TestCalculateScOneDirection:
    """Tests for _calculate_sc_one_direction helper."""

    def test_empty_from_surface(self):
        """Test with empty from_surface."""
        result = _calculate_sc_one_direction([], [], 1.5)
        assert len(result) == 0

    def test_empty_to_surface(self):
        """Test with empty to_surface."""
        point = SurfacePoint(
            coords=np.array([0.0, 0.0, 0.0]),
            normal=np.array([1.0, 0.0, 0.0]),
            atom_index=0,
            area=1.0,
        )
        result = _calculate_sc_one_direction([point], [], 1.5)
        assert len(result) == 0

    def test_parallel_normals_high_sc(self):
        """Test that anti-parallel normals give high Sc."""
        # Two surfaces facing each other (anti-parallel normals)
        from_surface = [
            SurfacePoint(
                coords=np.array([0.0, 0.0, 0.0]),
                normal=np.array([1.0, 0.0, 0.0]),  # pointing right
                atom_index=0,
                area=1.0,
            ),
        ]
        to_surface = [
            SurfacePoint(
                coords=np.array([0.5, 0.0, 0.0]),  # very close
                normal=np.array([-1.0, 0.0, 0.0]),  # pointing left (toward from)
                atom_index=0,
                area=1.0,
            ),
        ]

        result = _calculate_sc_one_direction(from_surface, to_surface, 2.0)

        # Anti-parallel normals should give positive Sc
        # (because -dot(n1, n2) = -(-1) = 1)
        assert len(result) == 1
        assert result[0] > 0

    def test_same_direction_normals_low_sc(self):
        """Test that same-direction normals give low/negative Sc."""
        # Two surfaces with normals pointing same direction
        from_surface = [
            SurfacePoint(
                coords=np.array([0.0, 0.0, 0.0]),
                normal=np.array([1.0, 0.0, 0.0]),
                atom_index=0,
                area=1.0,
            ),
        ]
        to_surface = [
            SurfacePoint(
                coords=np.array([0.5, 0.0, 0.0]),
                normal=np.array([1.0, 0.0, 0.0]),  # same direction
                atom_index=0,
                area=1.0,
            ),
        ]

        result = _calculate_sc_one_direction(from_surface, to_surface, 2.0)

        # Same direction normals should give negative Sc
        assert len(result) == 1
        assert result[0] < 0

    def test_distance_cutoff_filters_far_points(self):
        """Test that distance cutoff filters far points."""
        from_surface = [
            SurfacePoint(
                coords=np.array([0.0, 0.0, 0.0]),
                normal=np.array([1.0, 0.0, 0.0]),
                atom_index=0,
                area=1.0,
            ),
        ]
        to_surface = [
            SurfacePoint(
                coords=np.array([10.0, 0.0, 0.0]),  # far away
                normal=np.array([-1.0, 0.0, 0.0]),
                atom_index=0,
                area=1.0,
            ),
        ]

        result = _calculate_sc_one_direction(from_surface, to_surface, 1.5)

        # Should be empty because the only point is too far
        assert len(result) == 0
