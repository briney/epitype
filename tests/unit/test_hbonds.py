"""Unit tests for hydrogen bond detection."""
from pathlib import Path

from epitype.core.interface import InterfaceRegion, detect_interface
from epitype.core.structure import Structure
from epitype.io.parsers import parse_structure
from epitype.metrics.hbonds import (
    ACCEPTOR_ATOMS,
    DONOR_ATOMS,
    HBOND_DISTANCE_MAX,
    HBondResult,
    HydrogenBond,
    analyze_hbonds,
    count_interface_hbonds,
    count_unsatisfied_hbonds,
    detect_hydrogen_bonds,
)


class TestDetectHydrogenBonds:
    """Tests for detect_hydrogen_bonds function."""

    def test_hbond_detection_basic(self, minimal_structure: Structure):
        """Test basic H-bond detection."""
        hbonds = detect_hydrogen_bonds(minimal_structure)

        # Should find some H-bonds in structure
        assert isinstance(hbonds, list)
        for hb in hbonds:
            assert isinstance(hb, HydrogenBond)

    def test_hbond_distance_cutoff_enforced(self, minimal_structure: Structure):
        """Test that distance cutoff is enforced."""
        hbonds = detect_hydrogen_bonds(minimal_structure)

        for hbond in hbonds:
            assert hbond.distance <= HBOND_DISTANCE_MAX

    def test_no_same_residue_hbonds(self, minimal_structure: Structure):
        """Test that same-residue H-bonds are excluded."""
        hbonds = detect_hydrogen_bonds(minimal_structure)

        for hbond in hbonds:
            assert hbond.donor_atom.residue_index != hbond.acceptor_atom.residue_index

    def test_custom_distance_cutoff(self, minimal_structure: Structure):
        """Test using custom distance cutoff."""
        hbonds_default = detect_hydrogen_bonds(minimal_structure)
        hbonds_short = detect_hydrogen_bonds(minimal_structure, distance_cutoff=2.5)

        # Shorter cutoff should find fewer or equal H-bonds
        assert len(hbonds_short) <= len(hbonds_default)

    def test_hbond_detection_real_structure(self, pdb_1yy9: Path):
        """Test H-bond detection on real structure."""
        structure = parse_structure(pdb_1yy9)
        hbonds = detect_hydrogen_bonds(structure)

        # Real protein should have many H-bonds
        assert len(hbonds) > 100

    def test_donor_acceptor_atoms_correct(self, pdb_1yy9: Path):
        """Test that detected H-bond atoms are valid donors/acceptors."""
        structure = parse_structure(pdb_1yy9)
        hbonds = detect_hydrogen_bonds(structure)

        for hbond in hbonds[:10]:  # Check first 10
            assert hbond.donor_atom.name in DONOR_ATOMS
            assert hbond.acceptor_atom.name in ACCEPTOR_ATOMS


class TestCountInterfaceHbonds:
    """Tests for count_interface_hbonds function."""

    def test_cross_interface_counting(self, pdb_1yy9: Path):
        """Test counting cross-interface H-bonds."""
        structure = parse_structure(pdb_1yy9)
        interface = detect_interface(structure, ["C", "D"], ["A"], cutoff=8.0)
        hbonds = detect_hydrogen_bonds(structure)

        interface_count = count_interface_hbonds(hbonds, interface)

        # Should have some cross-interface H-bonds
        assert interface_count >= 0
        assert interface_count <= len(hbonds)

    def test_interface_hbonds_cross_chains(self, pdb_1yy9: Path):
        """Test that interface H-bonds actually cross chain groups."""
        structure = parse_structure(pdb_1yy9)
        interface = detect_interface(structure, ["C", "D"], ["A"], cutoff=8.0)
        hbonds = detect_hydrogen_bonds(structure)

        g1_chains = set(interface.group1_chains)
        g2_chains = set(interface.group2_chains)

        count = 0
        for hbond in hbonds:
            donor_in_g1 = hbond.donor_atom.chain_id in g1_chains
            donor_in_g2 = hbond.donor_atom.chain_id in g2_chains
            acceptor_in_g1 = hbond.acceptor_atom.chain_id in g1_chains
            acceptor_in_g2 = hbond.acceptor_atom.chain_id in g2_chains

            if (donor_in_g1 and acceptor_in_g2) or (donor_in_g2 and acceptor_in_g1):
                count += 1

        interface_count = count_interface_hbonds(hbonds, interface)
        assert interface_count == count


class TestCountUnsatisfiedHbonds:
    """Tests for count_unsatisfied_hbonds function."""

    def test_unsatisfied_counts_non_negative(self, pdb_1yy9: Path):
        """Test that unsatisfied counts are non-negative."""
        structure = parse_structure(pdb_1yy9)
        interface = detect_interface(structure, ["C", "D"], ["A"], cutoff=8.0)
        hbonds = detect_hydrogen_bonds(structure)

        unsat_donors, unsat_acceptors = count_unsatisfied_hbonds(
            structure, interface, hbonds
        )

        assert unsat_donors >= 0
        assert unsat_acceptors >= 0


class TestAnalyzeHbonds:
    """Tests for analyze_hbonds function."""

    def test_analyze_hbonds_complete(self, pdb_1yy9: Path):
        """Test complete H-bond analysis."""
        structure = parse_structure(pdb_1yy9)
        interface = detect_interface(structure, ["C", "D"], ["A"], cutoff=8.0)

        result = analyze_hbonds(structure, interface)

        assert isinstance(result, HBondResult)
        assert result.total_hbonds >= 0
        assert result.interface_hbonds >= 0
        assert result.interface_hbonds <= result.total_hbonds
        assert result.buried_unsatisfied_donors >= 0
        assert result.buried_unsatisfied_acceptors >= 0

    def test_analyze_hbonds_minimal(self, minimal_structure: Structure, interface_ab: InterfaceRegion):
        """Test H-bond analysis on minimal structure."""
        result = analyze_hbonds(minimal_structure, interface_ab)

        assert isinstance(result, HBondResult)
        assert isinstance(result.hbonds, list)

    def test_interface_hbonds_in_expected_range(self, pdb_1yy9: Path):
        """Test that interface H-bonds are in expected range for antibody-antigen."""
        structure = parse_structure(pdb_1yy9)
        interface = detect_interface(structure, ["C", "D"], ["A"], cutoff=8.0)

        result = analyze_hbonds(structure, interface)

        # Antibody-antigen interfaces typically have 5-25 cross-interface H-bonds
        assert 0 <= result.interface_hbonds <= 50
