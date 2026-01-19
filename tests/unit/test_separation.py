"""Tests for chain separation module."""
import numpy as np
import pytest

from epitype.core.interface import InterfaceRegion, detect_interface
from epitype.core.separation import (
    compute_interface_normal,
    compute_separation_distance,
    separate_chains,
)
from epitype.core.structure import Atom, Chain, Residue, Structure


def _make_two_chain_structure() -> Structure:
    """Create a simple two-chain structure for testing."""
    structure = Structure(name="test")

    # Chain A at origin
    chain_a = Chain(chain_id="A")
    for i in range(3):
        residue = Residue(
            index=i,
            name="ALA",
            chain_id="A",
            seq_num=i + 1,
        )
        residue.atoms.append(
            Atom(
                index=i * 2,
                name="CA",
                element="C",
                coords=np.array([0.0, float(i) * 3.0, 0.0]),
                residue_index=i,
                chain_id="A",
            )
        )
        residue.atoms.append(
            Atom(
                index=i * 2 + 1,
                name="CB",
                element="C",
                coords=np.array([1.0, float(i) * 3.0, 0.0]),
                residue_index=i,
                chain_id="A",
            )
        )
        chain_a.residues.append(residue)
    structure.chains["A"] = chain_a

    # Chain B offset in x direction
    chain_b = Chain(chain_id="B")
    for i in range(3):
        residue = Residue(
            index=i + 3,
            name="GLY",
            chain_id="B",
            seq_num=i + 1,
        )
        residue.atoms.append(
            Atom(
                index=(i + 3) * 2,
                name="CA",
                element="C",
                coords=np.array([5.0, float(i) * 3.0, 0.0]),
                residue_index=i + 3,
                chain_id="B",
            )
        )
        residue.atoms.append(
            Atom(
                index=(i + 3) * 2 + 1,
                name="CB",
                element="C",
                coords=np.array([4.0, float(i) * 3.0, 0.0]),
                residue_index=i + 3,
                chain_id="B",
            )
        )
        chain_b.residues.append(residue)
    structure.chains["B"] = chain_b

    return structure


class TestComputeInterfaceNormal:
    """Tests for compute_interface_normal function."""

    def test_normal_points_from_g1_to_g2(self):
        """Test that normal points from group1 to group2."""
        structure = _make_two_chain_structure()
        interface = detect_interface(
            structure,
            group1_chains=["A"],
            group2_chains=["B"],
            cutoff=10.0,
        )

        normal = compute_interface_normal(structure, interface)

        # Chain A is at x=0-1, Chain B is at x=4-5
        # Normal should point in positive x direction
        assert normal[0] > 0
        # Should be unit vector
        assert np.linalg.norm(normal) == pytest.approx(1.0)

    def test_empty_interface_default_normal(self):
        """Test default normal when no interface residues."""
        structure = _make_two_chain_structure()
        interface = InterfaceRegion(
            group1_chains=["A"],
            group2_chains=["B"],
        )
        # Empty interface (no residues)

        normal = compute_interface_normal(structure, interface)

        # Should return default z-axis
        np.testing.assert_array_almost_equal(normal, [0.0, 0.0, 1.0])

    def test_normal_is_normalized(self):
        """Test that returned normal is unit length."""
        structure = _make_two_chain_structure()
        interface = detect_interface(
            structure,
            group1_chains=["A"],
            group2_chains=["B"],
            cutoff=10.0,
        )

        normal = compute_interface_normal(structure, interface)

        assert np.linalg.norm(normal) == pytest.approx(1.0)


class TestSeparateChains:
    """Tests for separate_chains function."""

    def test_separation_increases_distance(self):
        """Test that separation increases inter-chain distance."""
        structure = _make_two_chain_structure()
        interface = detect_interface(
            structure,
            group1_chains=["A"],
            group2_chains=["B"],
            cutoff=10.0,
        )

        # Get original minimum distance
        chain_a_coords = structure.get_chain("A").coords
        chain_b_coords = structure.get_chain("B").coords
        from scipy.spatial.distance import cdist

        original_min_dist = cdist(chain_a_coords, chain_b_coords).min()

        # Separate
        separated = separate_chains(structure, interface, distance=100.0)

        # Get new minimum distance
        sep_chain_a = separated.get_chain("A").coords
        sep_chain_b = separated.get_chain("B").coords
        new_min_dist = cdist(sep_chain_a, sep_chain_b).min()

        assert new_min_dist > original_min_dist
        assert new_min_dist > original_min_dist + 50  # Should be significantly larger

    def test_original_unchanged(self):
        """Test that original structure is not modified."""
        structure = _make_two_chain_structure()
        interface = detect_interface(
            structure,
            group1_chains=["A"],
            group2_chains=["B"],
            cutoff=10.0,
        )

        # Get original position
        original_pos = list(structure.get_chain("B").atoms)[0].coords.copy()

        # Separate
        separate_chains(structure, interface, distance=100.0)

        # Check original is unchanged
        new_pos = list(structure.get_chain("B").atoms)[0].coords
        np.testing.assert_array_equal(original_pos, new_pos)

    def test_group1_unchanged(self):
        """Test that group1 chains are not translated."""
        structure = _make_two_chain_structure()
        interface = detect_interface(
            structure,
            group1_chains=["A"],
            group2_chains=["B"],
            cutoff=10.0,
        )

        original_a_coords = structure.get_chain("A").coords.copy()

        separated = separate_chains(structure, interface, distance=100.0)

        new_a_coords = separated.get_chain("A").coords
        np.testing.assert_array_almost_equal(original_a_coords, new_a_coords)

    def test_group2_translated(self):
        """Test that group2 chains are translated."""
        structure = _make_two_chain_structure()
        interface = detect_interface(
            structure,
            group1_chains=["A"],
            group2_chains=["B"],
            cutoff=10.0,
        )

        original_b_coords = structure.get_chain("B").coords.copy()

        separated = separate_chains(structure, interface, distance=100.0)

        new_b_coords = separated.get_chain("B").coords

        # Should be different
        assert not np.allclose(original_b_coords, new_b_coords)

        # Translation should be roughly 100 Angstroms
        displacement = np.linalg.norm(new_b_coords[0] - original_b_coords[0])
        assert displacement == pytest.approx(100.0, rel=0.1)

    def test_custom_separation_distance(self):
        """Test using custom separation distance."""
        structure = _make_two_chain_structure()
        interface = detect_interface(
            structure,
            group1_chains=["A"],
            group2_chains=["B"],
            cutoff=10.0,
        )

        original_b_pos = list(structure.get_chain("B").atoms)[0].coords.copy()

        # Use small distance
        separated = separate_chains(structure, interface, distance=10.0)

        new_b_pos = list(separated.get_chain("B").atoms)[0].coords
        displacement = np.linalg.norm(new_b_pos - original_b_pos)

        assert displacement == pytest.approx(10.0, rel=0.1)


class TestComputeSeparationDistance:
    """Tests for compute_separation_distance function."""

    def test_returns_appropriate_distance(self):
        """Test that computed distance is reasonable."""
        structure = _make_two_chain_structure()
        interface = detect_interface(
            structure,
            group1_chains=["A"],
            group2_chains=["B"],
            cutoff=10.0,
        )

        dist = compute_separation_distance(structure, interface, target_min_distance=100.0)

        # Should be at least SEPARATION_DISTANCE (500)
        assert dist >= 500.0

    def test_includes_target_min_distance(self):
        """Test that result accounts for target minimum distance."""
        structure = _make_two_chain_structure()
        interface = detect_interface(
            structure,
            group1_chains=["A"],
            group2_chains=["B"],
            cutoff=10.0,
        )

        # Get current minimum distance
        chain_a_coords = structure.get_chain("A").coords
        chain_b_coords = structure.get_chain("B").coords
        from scipy.spatial.distance import cdist

        current_min = cdist(chain_a_coords, chain_b_coords).min()

        dist = compute_separation_distance(structure, interface, target_min_distance=100.0)

        # Distance should ensure we reach target_min_distance
        assert dist >= 100.0 - current_min + 500.0

    def test_empty_groups_returns_default(self):
        """Test that empty groups return default separation distance."""
        structure = _make_two_chain_structure()
        interface = InterfaceRegion(
            group1_chains=["C"],  # Non-existent
            group2_chains=["D"],  # Non-existent
        )

        dist = compute_separation_distance(structure, interface)

        assert dist == 500.0  # SEPARATION_DISTANCE


class TestSeparationIntegration:
    """Integration tests for the separation workflow."""

    def test_full_separation_workflow(self):
        """Test complete workflow: detect interface, separate, verify."""
        structure = _make_two_chain_structure()

        # Detect interface
        interface = detect_interface(
            structure,
            group1_chains=["A"],
            group2_chains=["B"],
            cutoff=10.0,
        )

        # Verify interface detected
        assert interface.num_interface_residues > 0

        # Compute and apply separation
        dist = compute_separation_distance(structure, interface)
        separated = separate_chains(structure, interface, distance=dist)

        # Verify separation
        from scipy.spatial.distance import cdist

        sep_chain_a = separated.get_chain("A").coords
        sep_chain_b = separated.get_chain("B").coords
        min_dist = cdist(sep_chain_a, sep_chain_b).min()

        # Should be very large now
        assert min_dist > 100.0
