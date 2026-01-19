"""Unit tests for structure parsers."""
from pathlib import Path

import pytest

from epitype.io.parsers import (
    parse_pdb,
    parse_structure,
)


class TestParsePDB:
    """Tests for parse_pdb function."""

    def test_parse_1yy9(self, pdb_1yy9: Path):
        """Test parsing 1YY9 PDB file."""
        structure = parse_pdb(pdb_1yy9)

        assert structure.name == "1YY9"
        assert len(structure.chains) > 0
        assert structure.num_atoms > 0
        assert structure.num_residues > 0

    def test_chain_ids_preserved(self, pdb_1yy9: Path):
        """Test that chain IDs are correctly preserved."""
        structure = parse_pdb(pdb_1yy9)

        # 1YY9 has chains H, L (antibody) and A (antigen)
        chain_ids = set(structure.chains.keys())
        assert "H" in chain_ids or "A" in chain_ids  # At least one expected chain

    def test_nonexistent_file_raises(self):
        """Test that nonexistent file raises error."""
        with pytest.raises(Exception):
            parse_pdb("/nonexistent/file.pdb")

    def test_atoms_have_coordinates(self, pdb_1yy9: Path):
        """Test that parsed atoms have valid coordinates."""
        structure = parse_pdb(pdb_1yy9)

        for atom in structure.atoms:
            assert atom.coords is not None
            assert len(atom.coords) == 3
            # Coordinates should be finite
            import numpy as np
            assert np.all(np.isfinite(atom.coords))
            break  # Just check first atom

    def test_residues_have_atoms(self, pdb_1yy9: Path):
        """Test that residues contain atoms."""
        structure = parse_pdb(pdb_1yy9)

        for residue in structure.residues:
            assert len(residue.atoms) > 0
            break  # Just check first residue


class TestParseStructure:
    """Tests for auto-detect parse_structure function."""

    def test_auto_detect_pdb(self, pdb_1yy9: Path):
        """Test that PDB files are auto-detected."""
        structure = parse_structure(pdb_1yy9)
        assert structure.num_atoms > 0

    def test_parse_4fqi(self, pdb_4fqi: Path):
        """Test parsing 4FQI structure."""
        structure = parse_structure(pdb_4fqi)
        assert structure.name == "4FQI"
        assert structure.num_atoms > 0

    def test_structure_has_chains(self, pdb_1yy9: Path):
        """Test that parsed structure has chains."""
        structure = parse_structure(pdb_1yy9)
        assert len(structure.chains) > 0

    def test_residue_names_valid(self, pdb_1yy9: Path):
        """Test that residue names are valid 3-letter codes."""
        structure = parse_structure(pdb_1yy9)

        for residue in structure.residues:
            assert len(residue.name) == 3
            assert residue.name.isupper() or residue.name[0].isupper()
            break  # Just check first residue


class TestParserEdgeCases:
    """Tests for edge cases in parsing."""

    def test_parse_preserves_b_factors(self, pdb_1yy9: Path):
        """Test that B-factors are preserved."""
        structure = parse_pdb(pdb_1yy9)

        has_nonzero_bfactor = False
        for atom in structure.atoms:
            if atom.b_factor != 0.0:
                has_nonzero_bfactor = True
                break

        assert has_nonzero_bfactor

    def test_parse_preserves_occupancy(self, pdb_1yy9: Path):
        """Test that occupancy values are preserved."""
        structure = parse_pdb(pdb_1yy9)

        for atom in structure.atoms:
            # Most atoms have occupancy 1.0
            assert 0.0 <= atom.occupancy <= 1.0
            break

    def test_parse_preserves_element(self, pdb_1yy9: Path):
        """Test that element types are preserved."""
        structure = parse_pdb(pdb_1yy9)

        elements = set()
        for atom in structure.atoms:
            elements.add(atom.element)

        # Should have at least C, N, O in a protein
        assert "C" in elements
        assert "N" in elements
        assert "O" in elements
