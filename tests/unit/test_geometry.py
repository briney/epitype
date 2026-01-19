"""Tests for geometry utility functions."""
import numpy as np
import pytest

from epitype.utils.geometry import (
    angle_between,
    centroid,
    dihedral_angle,
    normalize,
    rmsd,
    rotation_matrix,
)


class TestNormalize:
    """Tests for normalize function."""

    def test_unit_vector_unchanged(self):
        """Test that unit vectors are unchanged."""
        v = np.array([1.0, 0.0, 0.0])
        result = normalize(v)
        np.testing.assert_array_almost_equal(result, v)

    def test_normalizes_to_unit_length(self):
        """Test that result has unit length."""
        v = np.array([3.0, 4.0, 0.0])
        result = normalize(v)
        assert np.linalg.norm(result) == pytest.approx(1.0)
        np.testing.assert_array_almost_equal(result, [0.6, 0.8, 0.0])

    def test_preserves_direction(self):
        """Test that direction is preserved."""
        v = np.array([2.0, 2.0, 2.0])
        result = normalize(v)
        # All components should be equal for this input
        assert result[0] == pytest.approx(result[1])
        assert result[1] == pytest.approx(result[2])

    def test_zero_vector_returns_zero(self):
        """Test that zero vector returns zero."""
        v = np.array([0.0, 0.0, 0.0])
        result = normalize(v)
        np.testing.assert_array_almost_equal(result, [0.0, 0.0, 0.0])

    def test_near_zero_vector_returns_zero(self):
        """Test that near-zero vector returns zero."""
        v = np.array([1e-12, 1e-12, 1e-12])
        result = normalize(v)
        np.testing.assert_array_almost_equal(result, [0.0, 0.0, 0.0])


class TestAngleBetween:
    """Tests for angle_between function."""

    def test_parallel_vectors(self):
        """Test angle between parallel vectors is 0."""
        v1 = np.array([1.0, 0.0, 0.0])
        v2 = np.array([2.0, 0.0, 0.0])
        assert angle_between(v1, v2) == pytest.approx(0.0)

    def test_antiparallel_vectors(self):
        """Test angle between antiparallel vectors is 180."""
        v1 = np.array([1.0, 0.0, 0.0])
        v2 = np.array([-1.0, 0.0, 0.0])
        assert angle_between(v1, v2) == pytest.approx(180.0)

    def test_perpendicular_vectors(self):
        """Test angle between perpendicular vectors is 90."""
        v1 = np.array([1.0, 0.0, 0.0])
        v2 = np.array([0.0, 1.0, 0.0])
        assert angle_between(v1, v2) == pytest.approx(90.0)

    def test_45_degree_angle(self):
        """Test 45 degree angle."""
        v1 = np.array([1.0, 0.0, 0.0])
        v2 = np.array([1.0, 1.0, 0.0])
        assert angle_between(v1, v2) == pytest.approx(45.0)

    def test_60_degree_angle(self):
        """Test 60 degree angle."""
        v1 = np.array([1.0, 0.0, 0.0])
        v2 = np.array([0.5, np.sqrt(3) / 2, 0.0])
        assert angle_between(v1, v2) == pytest.approx(60.0)

    def test_order_independent(self):
        """Test that angle is same regardless of order."""
        v1 = np.array([1.0, 2.0, 3.0])
        v2 = np.array([4.0, 5.0, 6.0])
        assert angle_between(v1, v2) == pytest.approx(angle_between(v2, v1))


class TestDihedralAngle:
    """Tests for dihedral_angle function."""

    def test_planar_dihedral_zero(self):
        """Test dihedral angle for cis configuration (0 degrees)."""
        # Cis configuration: p1 and p4 on same side of p2-p3 axis
        p1 = np.array([1.0, 0.0, 0.0])
        p2 = np.array([0.0, 0.0, 0.0])
        p3 = np.array([0.0, 1.0, 0.0])
        p4 = np.array([1.0, 1.0, 0.0])
        assert dihedral_angle(p1, p2, p3, p4) == pytest.approx(0.0, abs=1e-5)

    def test_dihedral_180(self):
        """Test dihedral angle of 180 degrees."""
        p1 = np.array([1.0, 0.0, 0.0])
        p2 = np.array([0.0, 0.0, 0.0])
        p3 = np.array([0.0, 1.0, 0.0])
        p4 = np.array([-1.0, 1.0, 0.0])
        result = dihedral_angle(p1, p2, p3, p4)
        # Should be close to 180 or -180
        assert abs(result) == pytest.approx(180.0, abs=1e-5)

    def test_dihedral_90(self):
        """Test dihedral angle of 90 degrees."""
        p1 = np.array([1.0, 0.0, 0.0])
        p2 = np.array([0.0, 0.0, 0.0])
        p3 = np.array([0.0, 1.0, 0.0])
        p4 = np.array([0.0, 1.0, 1.0])
        result = dihedral_angle(p1, p2, p3, p4)
        assert result == pytest.approx(90.0, abs=1e-5)


class TestRotationMatrix:
    """Tests for rotation_matrix function."""

    def test_identity_rotation(self):
        """Test rotation by 0 gives identity."""
        axis = np.array([0.0, 0.0, 1.0])
        R = rotation_matrix(axis, 0.0)
        np.testing.assert_array_almost_equal(R, np.eye(3))

    def test_90_degree_z_rotation(self):
        """Test 90 degree rotation around z-axis."""
        axis = np.array([0.0, 0.0, 1.0])
        R = rotation_matrix(axis, np.pi / 2)

        # Rotate x-axis should give y-axis
        x = np.array([1.0, 0.0, 0.0])
        rotated = R @ x
        np.testing.assert_array_almost_equal(rotated, [0.0, 1.0, 0.0])

    def test_180_degree_rotation(self):
        """Test 180 degree rotation."""
        axis = np.array([0.0, 0.0, 1.0])
        R = rotation_matrix(axis, np.pi)

        # Rotate x-axis should give -x-axis
        x = np.array([1.0, 0.0, 0.0])
        rotated = R @ x
        np.testing.assert_array_almost_equal(rotated, [-1.0, 0.0, 0.0])

    def test_rotation_preserves_length(self):
        """Test that rotation preserves vector length."""
        axis = np.array([1.0, 1.0, 1.0])
        R = rotation_matrix(axis, np.pi / 3)

        v = np.array([3.0, 4.0, 5.0])
        rotated = R @ v

        assert np.linalg.norm(rotated) == pytest.approx(np.linalg.norm(v))

    def test_orthogonal_matrix(self):
        """Test that result is orthogonal (R @ R.T = I)."""
        axis = np.array([1.0, 2.0, 3.0])
        R = rotation_matrix(axis, 1.234)

        np.testing.assert_array_almost_equal(R @ R.T, np.eye(3))


class TestCentroid:
    """Tests for centroid function."""

    def test_single_point(self):
        """Test centroid of single point is that point."""
        coords = np.array([[1.0, 2.0, 3.0]])
        result = centroid(coords)
        np.testing.assert_array_almost_equal(result, [1.0, 2.0, 3.0])

    def test_two_points(self):
        """Test centroid of two points is midpoint."""
        coords = np.array([[0.0, 0.0, 0.0], [2.0, 4.0, 6.0]])
        result = centroid(coords)
        np.testing.assert_array_almost_equal(result, [1.0, 2.0, 3.0])

    def test_symmetric_points(self):
        """Test centroid of symmetric points is origin."""
        coords = np.array(
            [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, -1.0, 0.0]]
        )
        result = centroid(coords)
        np.testing.assert_array_almost_equal(result, [0.0, 0.0, 0.0])


class TestRmsd:
    """Tests for rmsd function."""

    def test_identical_coords(self):
        """Test RMSD of identical coordinates is 0."""
        coords = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        assert rmsd(coords, coords) == pytest.approx(0.0)

    def test_known_rmsd(self):
        """Test RMSD with known value."""
        coords1 = np.array([[0.0, 0.0, 0.0]])
        coords2 = np.array([[1.0, 0.0, 0.0]])
        # Distance is 1.0, RMSD = sqrt(1) = 1.0
        assert rmsd(coords1, coords2) == pytest.approx(1.0)

    def test_multiple_points(self):
        """Test RMSD with multiple points."""
        coords1 = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        coords2 = np.array([[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
        # Each point moved by 1.0, RMSD = sqrt((1 + 1) / 2) = 1.0
        assert rmsd(coords1, coords2) == pytest.approx(1.0)

    def test_order_independent(self):
        """Test that RMSD is symmetric."""
        coords1 = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        coords2 = np.array([[1.5, 2.5, 3.5], [4.5, 5.5, 6.5]])
        assert rmsd(coords1, coords2) == pytest.approx(rmsd(coords2, coords1))
