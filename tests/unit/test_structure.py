"""Tests for structure module."""
import numpy as np
import pytest

from epitype.core.structure import Atom, Chain, Residue, Structure


class TestAtom:
    """Tests for Atom dataclass."""

    def test_atom_creation(self):
        """Test basic atom creation."""
        atom = Atom(
            index=0,
            name="CA",
            element="C",
            coords=np.array([1.0, 2.0, 3.0]),
            residue_index=0,
            chain_id="A",
        )
        assert atom.index == 0
        assert atom.name == "CA"
        assert atom.element == "C"
        assert atom.chain_id == "A"
        np.testing.assert_array_equal(atom.coords, [1.0, 2.0, 3.0])

    def test_vdw_radius(self):
        """Test VdW radius lookup."""
        carbon = Atom(
            index=0,
            name="CA",
            element="C",
            coords=np.array([0.0, 0.0, 0.0]),
            residue_index=0,
            chain_id="A",
        )
        nitrogen = Atom(
            index=1,
            name="N",
            element="N",
            coords=np.array([0.0, 0.0, 0.0]),
            residue_index=0,
            chain_id="A",
        )
        assert carbon.vdw_radius == 1.70
        assert nitrogen.vdw_radius == 1.55

    def test_is_backbone(self):
        """Test backbone atom identification."""
        ca = Atom(
            index=0,
            name="CA",
            element="C",
            coords=np.array([0.0, 0.0, 0.0]),
            residue_index=0,
            chain_id="A",
        )
        cb = Atom(
            index=1,
            name="CB",
            element="C",
            coords=np.array([0.0, 0.0, 0.0]),
            residue_index=0,
            chain_id="A",
        )
        n = Atom(
            index=2,
            name="N",
            element="N",
            coords=np.array([0.0, 0.0, 0.0]),
            residue_index=0,
            chain_id="A",
        )
        assert ca.is_backbone is True
        assert cb.is_backbone is False
        assert n.is_backbone is True

    def test_is_hydrogen(self):
        """Test hydrogen identification."""
        h = Atom(
            index=0,
            name="H",
            element="H",
            coords=np.array([0.0, 0.0, 0.0]),
            residue_index=0,
            chain_id="A",
        )
        c = Atom(
            index=1,
            name="CA",
            element="C",
            coords=np.array([0.0, 0.0, 0.0]),
            residue_index=0,
            chain_id="A",
        )
        assert h.is_hydrogen is True
        assert c.is_hydrogen is False

    def test_distance_to(self):
        """Test distance calculation between atoms."""
        atom1 = Atom(
            index=0,
            name="CA",
            element="C",
            coords=np.array([0.0, 0.0, 0.0]),
            residue_index=0,
            chain_id="A",
        )
        atom2 = Atom(
            index=1,
            name="CB",
            element="C",
            coords=np.array([3.0, 4.0, 0.0]),
            residue_index=0,
            chain_id="A",
        )
        assert atom1.distance_to(atom2) == pytest.approx(5.0)


class TestResidue:
    """Tests for Residue dataclass."""

    def _make_residue(self, name: str = "ALA", has_cb: bool = True) -> Residue:
        """Create a test residue with atoms."""
        residue = Residue(
            index=0,
            name=name,
            chain_id="A",
            seq_num=1,
        )
        # Add backbone atoms
        residue.atoms.append(
            Atom(
                index=0,
                name="N",
                element="N",
                coords=np.array([0.0, 0.0, 0.0]),
                residue_index=0,
                chain_id="A",
            )
        )
        residue.atoms.append(
            Atom(
                index=1,
                name="CA",
                element="C",
                coords=np.array([1.0, 0.0, 0.0]),
                residue_index=0,
                chain_id="A",
            )
        )
        residue.atoms.append(
            Atom(
                index=2,
                name="C",
                element="C",
                coords=np.array([2.0, 0.0, 0.0]),
                residue_index=0,
                chain_id="A",
            )
        )
        residue.atoms.append(
            Atom(
                index=3,
                name="O",
                element="O",
                coords=np.array([2.5, 1.0, 0.0]),
                residue_index=0,
                chain_id="A",
            )
        )
        if has_cb:
            residue.atoms.append(
                Atom(
                    index=4,
                    name="CB",
                    element="C",
                    coords=np.array([1.0, 1.0, 1.0]),
                    residue_index=0,
                    chain_id="A",
                )
            )
        return residue

    def test_ca_access(self):
        """Test CA atom access."""
        residue = self._make_residue()
        ca = residue.ca
        assert ca is not None
        assert ca.name == "CA"
        np.testing.assert_array_equal(ca.coords, [1.0, 0.0, 0.0])

    def test_cb_access(self):
        """Test CB atom access."""
        residue = self._make_residue()
        cb = residue.cb
        assert cb is not None
        assert cb.name == "CB"
        np.testing.assert_array_equal(cb.coords, [1.0, 1.0, 1.0])

    def test_cb_fallback_for_glycine(self):
        """Test CB falls back to CA for glycine (no CB)."""
        residue = self._make_residue(name="GLY", has_cb=False)
        cb = residue.cb
        assert cb is not None
        assert cb.name == "CA"  # Should fall back to CA

    def test_center_of_mass(self):
        """Test center of mass calculation."""
        residue = self._make_residue()
        com = residue.center_of_mass
        # Expected: mean of all atom coords
        expected = np.array([0.0, 0.0, 0.0]) + np.array([1.0, 0.0, 0.0]) + np.array(
            [2.0, 0.0, 0.0]
        ) + np.array([2.5, 1.0, 0.0]) + np.array([1.0, 1.0, 1.0])
        expected = expected / 5.0
        np.testing.assert_array_almost_equal(com, expected)

    def test_full_id(self):
        """Test full residue identifier."""
        residue = self._make_residue()
        assert residue.full_id == "A:ALA:1"

    def test_full_id_with_insertion_code(self):
        """Test full residue identifier with insertion code."""
        residue = self._make_residue()
        residue.insertion_code = "A"
        assert residue.full_id == "A:ALA:1A"

    def test_get_atom(self):
        """Test get_atom method."""
        residue = self._make_residue()
        ca = residue.get_atom("CA")
        assert ca is not None
        assert ca.name == "CA"

        missing = residue.get_atom("XYZ")
        assert missing is None


class TestChain:
    """Tests for Chain dataclass."""

    def _make_chain(self) -> Chain:
        """Create a test chain with residues."""
        chain = Chain(chain_id="A")
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
                    coords=np.array([float(i), 0.0, 0.0]),
                    residue_index=i,
                    chain_id="A",
                )
            )
            residue.atoms.append(
                Atom(
                    index=i * 2 + 1,
                    name="CB",
                    element="C",
                    coords=np.array([float(i), 1.0, 0.0]),
                    residue_index=i,
                    chain_id="A",
                )
            )
            chain.residues.append(residue)
        return chain

    def test_atoms_iterator(self):
        """Test iteration over all atoms in chain."""
        chain = self._make_chain()
        atoms = list(chain.atoms)
        assert len(atoms) == 6  # 3 residues * 2 atoms

    def test_coords(self):
        """Test coordinate extraction."""
        chain = self._make_chain()
        coords = chain.coords
        assert coords.shape == (6, 3)
        # First atom should be at (0, 0, 0)
        np.testing.assert_array_equal(coords[0], [0.0, 0.0, 0.0])

    def test_len(self):
        """Test chain length (number of residues)."""
        chain = self._make_chain()
        assert len(chain) == 3

    def test_empty_chain_coords(self):
        """Test coords for empty chain."""
        chain = Chain(chain_id="A")
        coords = chain.coords
        assert coords.shape == (0, 3)


class TestStructure:
    """Tests for Structure dataclass."""

    def _make_structure(self) -> Structure:
        """Create a test structure with two chains."""
        structure = Structure(name="test")

        # Chain A
        chain_a = Chain(chain_id="A")
        for i in range(2):
            residue = Residue(
                index=i,
                name="ALA",
                chain_id="A",
                seq_num=i + 1,
            )
            residue.atoms.append(
                Atom(
                    index=i,
                    name="CA",
                    element="C",
                    coords=np.array([float(i), 0.0, 0.0]),
                    residue_index=i,
                    chain_id="A",
                )
            )
            chain_a.residues.append(residue)
        structure.chains["A"] = chain_a

        # Chain B
        chain_b = Chain(chain_id="B")
        for i in range(2):
            residue = Residue(
                index=i + 2,
                name="GLY",
                chain_id="B",
                seq_num=i + 1,
            )
            residue.atoms.append(
                Atom(
                    index=i + 2,
                    name="CA",
                    element="C",
                    coords=np.array([float(i) + 10.0, 0.0, 0.0]),
                    residue_index=i + 2,
                    chain_id="B",
                )
            )
            chain_b.residues.append(residue)
        structure.chains["B"] = chain_b

        return structure

    def test_atoms_iterator(self):
        """Test iteration over all atoms."""
        structure = self._make_structure()
        atoms = list(structure.atoms)
        assert len(atoms) == 4

    def test_residues_iterator(self):
        """Test iteration over all residues."""
        structure = self._make_structure()
        residues = list(structure.residues)
        assert len(residues) == 4

    def test_coords(self):
        """Test coordinate extraction."""
        structure = self._make_structure()
        coords = structure.coords
        assert coords.shape == (4, 3)

    def test_num_atoms(self):
        """Test atom count."""
        structure = self._make_structure()
        assert structure.num_atoms == 4

    def test_num_residues(self):
        """Test residue count."""
        structure = self._make_structure()
        assert structure.num_residues == 4

    def test_get_chain(self):
        """Test get_chain method."""
        structure = self._make_structure()
        chain_a = structure.get_chain("A")
        assert chain_a is not None
        assert chain_a.chain_id == "A"

        chain_c = structure.get_chain("C")
        assert chain_c is None

    def test_get_chains(self):
        """Test get_chains method."""
        structure = self._make_structure()
        chains = structure.get_chains(["A", "B"])
        assert len(chains) == 2

        chains = structure.get_chains(["A", "C"])
        assert len(chains) == 1  # Only A exists

    def test_subset(self):
        """Test structure subsetting."""
        structure = self._make_structure()
        subset = structure.subset(["A"])
        assert "A" in subset.chains
        assert "B" not in subset.chains
        assert subset.num_atoms == 2
        assert subset.name == "test_subset"

    def test_copy(self):
        """Test deep copy."""
        structure = self._make_structure()
        copied = structure.copy()

        # Modify original
        list(structure.atoms)[0].coords[0] = 999.0

        # Check copy is unchanged
        assert list(copied.atoms)[0].coords[0] != 999.0

    def test_translate(self):
        """Test in-place translation."""
        structure = self._make_structure()
        original_first = list(structure.atoms)[0].coords.copy()

        translation = np.array([10.0, 20.0, 30.0])
        structure.translate(translation)

        new_first = list(structure.atoms)[0].coords
        np.testing.assert_array_almost_equal(
            new_first, original_first + translation
        )

    def test_translate_does_not_affect_copy(self):
        """Test that translation doesn't affect a previous copy."""
        structure = self._make_structure()
        copied = structure.copy()

        translation = np.array([10.0, 20.0, 30.0])
        structure.translate(translation)

        # Copy should be unchanged
        assert list(copied.atoms)[0].coords[0] == pytest.approx(0.0)

    def test_empty_structure(self):
        """Test empty structure properties."""
        structure = Structure(name="empty")
        assert structure.num_atoms == 0
        assert structure.num_residues == 0
        assert structure.coords.shape == (0, 3)
