"""Tests for interface detection module."""
import numpy as np
import pytest

from epitype.core.interface import (
    InterfaceRegion,
    detect_interface,
    parse_chain_groups,
)
from epitype.core.structure import Atom, Chain, Residue, Structure


class TestParseChainGroups:
    """Tests for parse_chain_groups function."""

    def test_simple_chain_spec(self):
        """Test simple two-chain spec."""
        group1, group2 = parse_chain_groups("A_B")
        assert group1 == ["A"]
        assert group2 == ["B"]

    def test_multi_chain_spec(self):
        """Test multi-chain spec."""
        group1, group2 = parse_chain_groups("HL_A")
        assert group1 == ["H", "L"]
        assert group2 == ["A"]

    def test_both_multi_chain(self):
        """Test both groups with multiple chains."""
        group1, group2 = parse_chain_groups("AB_CD")
        assert group1 == ["A", "B"]
        assert group2 == ["C", "D"]

    def test_missing_underscore(self):
        """Test error on missing underscore."""
        with pytest.raises(ValueError, match="Expected format"):
            parse_chain_groups("AB")

    def test_multiple_underscores(self):
        """Test error on multiple underscores."""
        with pytest.raises(ValueError, match="exactly one underscore"):
            parse_chain_groups("A_B_C")

    def test_empty_group1(self):
        """Test error on empty first group."""
        with pytest.raises(ValueError, match="at least one chain"):
            parse_chain_groups("_A")

    def test_empty_group2(self):
        """Test error on empty second group."""
        with pytest.raises(ValueError, match="at least one chain"):
            parse_chain_groups("A_")


class TestInterfaceRegion:
    """Tests for InterfaceRegion dataclass."""

    def test_num_interface_residues(self):
        """Test interface residue counting."""
        interface = InterfaceRegion(
            group1_chains=["A"],
            group2_chains=["B"],
        )
        # Add mock residues
        interface.residues_group1 = [
            Residue(index=0, name="ALA", chain_id="A", seq_num=1),
            Residue(index=1, name="GLY", chain_id="A", seq_num=2),
        ]
        interface.residues_group2 = [
            Residue(index=2, name="VAL", chain_id="B", seq_num=1),
        ]
        assert interface.num_interface_residues == 3

    def test_residue_ids(self):
        """Test residue ID extraction."""
        interface = InterfaceRegion(
            group1_chains=["A"],
            group2_chains=["B"],
        )
        interface.residues_group1 = [
            Residue(index=0, name="ALA", chain_id="A", seq_num=1),
        ]
        interface.residues_group2 = [
            Residue(index=1, name="GLY", chain_id="B", seq_num=5),
        ]
        assert interface.residue_ids_group1 == ["A:ALA:1"]
        assert interface.residue_ids_group2 == ["B:GLY:5"]


class TestDetectInterface:
    """Tests for detect_interface function."""

    def _make_two_chain_structure(
        self, chain_a_coords: list[tuple[float, float, float]],
        chain_b_coords: list[tuple[float, float, float]],
    ) -> Structure:
        """Create structure with two chains at specified CA positions."""
        structure = Structure(name="test")

        # Chain A
        chain_a = Chain(chain_id="A")
        for i, (x, y, z) in enumerate(chain_a_coords):
            residue = Residue(
                index=i,
                name="ALA",
                chain_id="A",
                seq_num=i + 1,
            )
            # Add CA atom
            residue.atoms.append(
                Atom(
                    index=i * 2,
                    name="CA",
                    element="C",
                    coords=np.array([x, y, z]),
                    residue_index=i,
                    chain_id="A",
                )
            )
            # Add CB atom (offset from CA)
            residue.atoms.append(
                Atom(
                    index=i * 2 + 1,
                    name="CB",
                    element="C",
                    coords=np.array([x + 0.5, y, z]),
                    residue_index=i,
                    chain_id="A",
                )
            )
            chain_a.residues.append(residue)
        structure.chains["A"] = chain_a

        # Chain B
        chain_b = Chain(chain_id="B")
        offset = len(chain_a_coords)
        for i, (x, y, z) in enumerate(chain_b_coords):
            residue = Residue(
                index=i + offset,
                name="GLY",
                chain_id="B",
                seq_num=i + 1,
            )
            residue.atoms.append(
                Atom(
                    index=(i + offset) * 2,
                    name="CA",
                    element="C",
                    coords=np.array([x, y, z]),
                    residue_index=i + offset,
                    chain_id="B",
                )
            )
            residue.atoms.append(
                Atom(
                    index=(i + offset) * 2 + 1,
                    name="CB",
                    element="C",
                    coords=np.array([x + 0.5, y, z]),
                    residue_index=i + offset,
                    chain_id="B",
                )
            )
            chain_b.residues.append(residue)
        structure.chains["B"] = chain_b

        return structure

    def test_interface_detection_close_residues(self):
        """Test detection of close residues as interface."""
        # Chain A at origin, Chain B close by
        structure = self._make_two_chain_structure(
            chain_a_coords=[(0, 0, 0), (0, 3, 0), (0, 6, 0)],
            chain_b_coords=[(5, 0, 0), (5, 3, 0), (100, 0, 0)],  # Last one is far
        )

        interface = detect_interface(
            structure,
            group1_chains=["A"],
            group2_chains=["B"],
            cutoff=8.0,
        )

        # First two residues of each chain should be in interface
        assert len(interface.residues_group1) >= 2
        assert len(interface.residues_group2) >= 2
        # Contact pairs should exist
        assert len(interface.contact_pairs) > 0

    def test_interface_detection_no_interface(self):
        """Test when chains are too far apart."""
        # Chains far apart
        structure = self._make_two_chain_structure(
            chain_a_coords=[(0, 0, 0)],
            chain_b_coords=[(100, 0, 0)],
        )

        interface = detect_interface(
            structure,
            group1_chains=["A"],
            group2_chains=["B"],
            cutoff=8.0,
        )

        assert len(interface.residues_group1) == 0
        assert len(interface.residues_group2) == 0
        assert len(interface.contact_pairs) == 0

    def test_interface_detection_cutoff_sensitivity(self):
        """Test that cutoff affects detection."""
        structure = self._make_two_chain_structure(
            chain_a_coords=[(0, 0, 0)],
            chain_b_coords=[(6, 0, 0)],  # ~5.5 Angstrom CB-CB distance
        )

        # With small cutoff, no interface
        interface_small = detect_interface(
            structure,
            group1_chains=["A"],
            group2_chains=["B"],
            cutoff=5.0,
        )
        assert interface_small.num_interface_residues == 0

        # With larger cutoff, should detect interface
        interface_large = detect_interface(
            structure,
            group1_chains=["A"],
            group2_chains=["B"],
            cutoff=10.0,
        )
        assert interface_large.num_interface_residues == 2

    def test_contact_pairs(self):
        """Test contact pair detection."""
        structure = self._make_two_chain_structure(
            chain_a_coords=[(0, 0, 0), (0, 5, 0)],
            chain_b_coords=[(3, 0, 0), (3, 5, 0)],
        )

        interface = detect_interface(
            structure,
            group1_chains=["A"],
            group2_chains=["B"],
            cutoff=10.0,
        )

        # Should have contact pairs
        assert len(interface.contact_pairs) > 0
        # Each pair should be (group1_idx, group2_idx)
        for g1_idx, g2_idx in interface.contact_pairs:
            assert isinstance(g1_idx, int)
            assert isinstance(g2_idx, int)

    def test_atom_indices_populated(self):
        """Test that atom indices are populated for interface residues."""
        structure = self._make_two_chain_structure(
            chain_a_coords=[(0, 0, 0)],
            chain_b_coords=[(3, 0, 0)],
        )

        interface = detect_interface(
            structure,
            group1_chains=["A"],
            group2_chains=["B"],
            cutoff=10.0,
        )

        # Atom indices should be populated
        assert len(interface.atoms_group1_indices) > 0
        assert len(interface.atoms_group2_indices) > 0

    def test_missing_chain(self):
        """Test handling of missing chain."""
        structure = self._make_two_chain_structure(
            chain_a_coords=[(0, 0, 0)],
            chain_b_coords=[(3, 0, 0)],
        )

        # Query for non-existent chain
        interface = detect_interface(
            structure,
            group1_chains=["A"],
            group2_chains=["C"],  # Doesn't exist
            cutoff=10.0,
        )

        assert interface.num_interface_residues == 0

    def test_use_ca_instead_of_cb(self):
        """Test using CA instead of CB for distance."""
        structure = self._make_two_chain_structure(
            chain_a_coords=[(0, 0, 0)],
            chain_b_coords=[(5, 0, 0)],
        )

        # Using CA (at exact positions)
        interface_ca = detect_interface(
            structure,
            group1_chains=["A"],
            group2_chains=["B"],
            cutoff=6.0,
            use_cb=False,
        )
        assert interface_ca.num_interface_residues == 2

        # Using CB (offset by 0.5)
        interface_cb = detect_interface(
            structure,
            group1_chains=["A"],
            group2_chains=["B"],
            cutoff=6.0,
            use_cb=True,
        )
        # CB-CB distance is 4.5 (5 - 0.5 - 0.5 = 4, but CB at x+0.5, so 5-0.5-0.5=4)
        assert interface_cb.num_interface_residues == 2

    def test_multi_chain_groups(self):
        """Test interface detection with multi-chain groups."""
        structure = Structure(name="test")

        # Create chains H, L (group 1) and A (group 2)
        for chain_id, x_offset in [("H", 0), ("L", 2), ("A", 5)]:
            chain = Chain(chain_id=chain_id)
            residue = Residue(
                index=0,
                name="ALA",
                chain_id=chain_id,
                seq_num=1,
            )
            residue.atoms.append(
                Atom(
                    index=0,
                    name="CA",
                    element="C",
                    coords=np.array([float(x_offset), 0.0, 0.0]),
                    residue_index=0,
                    chain_id=chain_id,
                )
            )
            residue.atoms.append(
                Atom(
                    index=1,
                    name="CB",
                    element="C",
                    coords=np.array([float(x_offset) + 0.5, 0.0, 0.0]),
                    residue_index=0,
                    chain_id=chain_id,
                )
            )
            chain.residues.append(residue)
            structure.chains[chain_id] = chain

        interface = detect_interface(
            structure,
            group1_chains=["H", "L"],
            group2_chains=["A"],
            cutoff=10.0,
        )

        # Both H and L should have interface residues with A
        assert interface.num_interface_residues >= 2
        assert interface.group1_chains == ["H", "L"]
        assert interface.group2_chains == ["A"]
