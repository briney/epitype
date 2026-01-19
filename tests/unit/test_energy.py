"""Unit tests for binding energy calculations."""
from pathlib import Path

import numpy as np
import pytest

from epitype.core.interface import detect_interface
from epitype.core.separation import separate_chains
from epitype.core.structure import Structure
from epitype.io.parsers import parse_structure
from epitype.metrics.energy import (
    KJ_TO_REU,
    BindingEnergyResult,
    SolventConfig,
    SolventModel,
    _write_temp_pdb,
    calculate_binding_energy,
    calculate_potential_energy,
    check_openmm_available,
    prepare_structure_for_openmm,
)

# Skip all tests if OpenMM is not available
pytestmark = pytest.mark.skipif(
    not check_openmm_available(),
    reason="OpenMM not available"
)


class TestCheckOpenMMAvailable:
    """Tests for OpenMM availability check."""

    def test_check_returns_bool(self):
        """Test that check returns a boolean."""
        result = check_openmm_available()
        assert isinstance(result, bool)


class TestWriteTempPDB:
    """Tests for temporary PDB file writing."""

    def test_write_pdb_creates_file(self, minimal_structure: Structure, temp_dir: Path):
        """Test that PDB file is created."""
        pdb_path = temp_dir / "test.pdb"
        _write_temp_pdb(minimal_structure, pdb_path)

        assert pdb_path.exists()

    def test_write_pdb_has_atoms(self, minimal_structure: Structure, temp_dir: Path):
        """Test that PDB file contains ATOM records."""
        pdb_path = temp_dir / "test.pdb"
        _write_temp_pdb(minimal_structure, pdb_path)

        with open(pdb_path) as f:
            content = f.read()

        assert "ATOM" in content
        assert "END" in content

    def test_write_pdb_preserves_chains(self, minimal_structure: Structure, temp_dir: Path):
        """Test that chain IDs are preserved."""
        pdb_path = temp_dir / "test.pdb"
        _write_temp_pdb(minimal_structure, pdb_path)

        with open(pdb_path) as f:
            content = f.read()

        # Should have both chains A and B
        assert " A " in content or " A" in content
        assert " B " in content or " B" in content

    def test_write_pdb_column_alignment(self, minimal_structure: Structure, temp_dir: Path):
        """Test that PDB columns are correctly aligned per PDB specification.

        PDB format specification:
        - Columns 1-6: Record name (ATOM)
        - Columns 7-11: Atom serial number (right-justified)
        - Column 12: Space
        - Columns 13-16: Atom name
        - Column 17: Alternate location indicator
        - Columns 18-20: Residue name (right-justified)
        - Column 21: Space
        - Column 22: Chain ID
        - Columns 23-26: Residue sequence number (right-justified)
        """
        pdb_path = temp_dir / "test.pdb"
        _write_temp_pdb(minimal_structure, pdb_path)

        with open(pdb_path) as f:
            lines = f.readlines()

        for line in lines:
            if not line.startswith("ATOM"):
                continue

            # Record name in columns 1-6
            assert line[0:6] == "ATOM  ", f"Record name wrong: '{line[0:6]}'"

            # Column 12 should be space (between serial and atom name)
            assert line[11] == " ", f"Column 12 not space: '{line[11]}'"

            # Column 17 should be space (alternate location indicator)
            assert line[16] == " ", f"Column 17 (altLoc) not space: '{line[16]}'"

            # Column 21 should be space (between residue name and chain)
            assert line[20] == " ", f"Column 21 not space: '{line[20]}'"

            # Chain ID in column 22
            chain_id = line[21]
            assert chain_id in ("A", "B"), f"Chain ID wrong: '{chain_id}'"

    def test_write_pdb_residue_name_position(self, minimal_structure: Structure, temp_dir: Path):
        """Test that residue name is in correct columns 18-20."""
        pdb_path = temp_dir / "test.pdb"
        _write_temp_pdb(minimal_structure, pdb_path)

        with open(pdb_path) as f:
            lines = f.readlines()

        for line in lines:
            if not line.startswith("ATOM"):
                continue

            # Residue name should be in columns 18-20 (indices 17-20)
            res_name = line[17:20]
            assert res_name == "ALA", f"Residue name in wrong position: '{res_name}'"

    def test_write_pdb_parseable_by_pdbfixer(self, minimal_structure: Structure, temp_dir: Path):
        """Test that written PDB can be parsed by PDBFixer without errors.

        This is a regression test for the 'Misaligned residue name' error
        that occurs when PDB columns are incorrectly formatted.
        """
        pdb_path = temp_dir / "test.pdb"
        _write_temp_pdb(minimal_structure, pdb_path)

        # Import PDBFixer and try to parse
        from pdbfixer import PDBFixer

        # This should not raise "Misaligned residue name" error
        fixer = PDBFixer(filename=str(pdb_path))

        # Verify it parsed correctly - should have atoms
        assert fixer.topology.getNumAtoms() > 0


class TestPrepareStructureForOpenMM:
    """Tests for structure preparation."""

    def test_prepare_returns_modeller_and_forcefield(self, pdb_1yy9: Path):
        """Test that prepare returns modeller and forcefield."""
        structure = parse_structure(pdb_1yy9)
        # Use a subset to make test faster
        structure = structure.subset(["A"])

        modeller, forcefield = prepare_structure_for_openmm(structure)

        # Should have topology and positions
        assert hasattr(modeller, "topology")
        assert hasattr(modeller, "positions")
        assert forcefield is not None


class TestCalculatePotentialEnergy:
    """Tests for potential energy calculation."""

    def test_energy_is_float(self, pdb_1yy9: Path):
        """Test that energy calculation returns a float."""
        structure = parse_structure(pdb_1yy9)
        # Use subset to make test faster
        structure = structure.subset(["A"])

        modeller, forcefield = prepare_structure_for_openmm(structure)
        energy = calculate_potential_energy(modeller, forcefield)

        assert isinstance(energy, float)

    def test_energy_with_minimization(self, pdb_1yy9: Path):
        """Test energy calculation with minimization."""
        structure = parse_structure(pdb_1yy9)
        structure = structure.subset(["A"])

        modeller, forcefield = prepare_structure_for_openmm(structure)
        energy_no_min = calculate_potential_energy(modeller, forcefield, minimize=False)

        # Reset positions (minimization modifies them)
        modeller2, forcefield2 = prepare_structure_for_openmm(structure)
        energy_min = calculate_potential_energy(modeller2, forcefield2, minimize=True)

        # Minimized energy should be lower (or equal)
        assert energy_min <= energy_no_min + 1e-3  # Small tolerance


class TestCalculateBindingEnergy:
    """Tests for binding energy calculation."""

    @pytest.mark.slow
    def test_binding_energy_returns_result(self, pdb_1yy9: Path):
        """Test that binding energy calculation returns BindingEnergyResult."""
        structure = parse_structure(pdb_1yy9)
        interface = detect_interface(structure, ["C", "D"], ["A"], cutoff=8.0)
        separated = separate_chains(structure, interface)

        # This test is slow, so we use a simplified calculation
        result = calculate_binding_energy(
            structure, separated, interface, minimize=False
        )

        assert isinstance(result, BindingEnergyResult)
        assert isinstance(result.dG_separated, float)
        assert isinstance(result.energy_complex, float)
        assert isinstance(result.energy_separated, float)

    @pytest.mark.slow
    def test_binding_energy_with_minimization(self, pdb_1yy9: Path):
        """Test that binding energy calculation with minimization returns valid result.

        Note: In vacuum calculations without solvation, binding energy may be
        positive or negative depending on the system and minimization results.
        The sign is not a reliable indicator of binding favorability without
        proper solvation terms.
        """
        structure = parse_structure(pdb_1yy9)
        interface = detect_interface(structure, ["C", "D"], ["A"], cutoff=8.0)
        separated = separate_chains(structure, interface)

        result = calculate_binding_energy(
            structure, separated, interface, minimize=True
        )

        # Just verify we get a valid result
        assert isinstance(result, BindingEnergyResult)
        assert not np.isnan(result.dG_separated)
        assert not np.isinf(result.dG_separated)

    @pytest.mark.slow
    def test_binding_energy_result_includes_metadata(self, pdb_1yy9: Path):
        """Test that result includes solvent metadata."""
        structure = parse_structure(pdb_1yy9)
        interface = detect_interface(structure, ["C", "D"], ["A"], cutoff=8.0)
        separated = separate_chains(structure, interface)

        result = calculate_binding_energy(
            structure, separated, interface, minimize=False
        )

        # Default should be vacuum with 0.15 M salt
        assert result.solvent_model == "vacuum"
        assert result.salt_concentration == 0.15


class TestKJtoREU:
    """Tests for energy conversion constant."""

    def test_conversion_constant_positive(self):
        """Test that conversion constant is positive."""
        assert KJ_TO_REU > 0

    def test_conversion_constant_reasonable(self):
        """Test that conversion constant is in reasonable range."""
        # 1 kJ/mol ~ 0.239 kcal/mol
        assert 0.2 < KJ_TO_REU < 0.3


class TestSolventModel:
    """Tests for SolventModel enum."""

    def test_all_models_defined(self):
        """Test that all expected solvent models are defined."""
        models = {m.value for m in SolventModel}
        assert models == {"vacuum", "obc2", "gbn", "hct"}

    def test_vacuum_is_string_enum(self):
        """Test that SolventModel values are strings."""
        assert SolventModel.VACUUM.value == "vacuum"
        assert SolventModel.OBC2.value == "obc2"
        assert SolventModel.GBN.value == "gbn"
        assert SolventModel.HCT.value == "hct"


class TestSolventConfig:
    """Tests for SolventConfig dataclass."""

    def test_default_is_vacuum(self):
        """Test that default configuration uses vacuum."""
        config = SolventConfig()
        assert config.model == SolventModel.VACUUM
        assert config.salt_concentration == 0.15
        assert config.solvent_dielectric == 80.0
        assert config.solute_dielectric == 1.0

    def test_custom_config(self):
        """Test creating custom solvent configuration."""
        config = SolventConfig(
            model=SolventModel.OBC2,
            salt_concentration=0.3,
            solvent_dielectric=78.0,
        )
        assert config.model == SolventModel.OBC2
        assert config.salt_concentration == 0.3
        assert config.solvent_dielectric == 78.0


class TestBindingEnergyResultDataclass:
    """Tests for BindingEnergyResult dataclass."""

    def test_result_has_expected_fields(self):
        """Test that BindingEnergyResult has all expected fields."""
        result = BindingEnergyResult(
            dG_separated=-10.5,
            energy_complex=-1000.0,
            energy_separated=-950.0,
            solvent_model="obc2",
            salt_concentration=0.15,
        )
        assert result.dG_separated == -10.5
        assert result.energy_complex == -1000.0
        assert result.energy_separated == -950.0
        assert result.solvent_model == "obc2"
        assert result.salt_concentration == 0.15


class TestImplicitSolventEnergy:
    """Tests for implicit solvent energy calculations."""

    @pytest.mark.slow
    def test_obc2_produces_different_energy_than_vacuum(self, pdb_1yy9: Path):
        """Test that OBC2 implicit solvent produces different energy than vacuum."""
        structure = parse_structure(pdb_1yy9)
        # Use single chain for faster test
        structure = structure.subset(["A"])

        # Vacuum calculation
        modeller1, forcefield1 = prepare_structure_for_openmm(
            structure, SolventModel.VACUUM
        )
        energy_vacuum = calculate_potential_energy(modeller1, forcefield1)

        # OBC2 calculation
        modeller2, forcefield2 = prepare_structure_for_openmm(
            structure, SolventModel.OBC2
        )
        energy_obc2 = calculate_potential_energy(modeller2, forcefield2)

        # Energies should be different (solvation adds significant contribution)
        assert energy_vacuum != energy_obc2

    @pytest.mark.slow
    def test_binding_energy_with_obc2(self, pdb_1yy9: Path):
        """Test binding energy calculation with OBC2 implicit solvent."""
        structure = parse_structure(pdb_1yy9)
        interface = detect_interface(structure, ["C", "D"], ["A"], cutoff=8.0)
        separated = separate_chains(structure, interface)

        solvent_config = SolventConfig(model=SolventModel.OBC2)
        result = calculate_binding_energy(
            structure, separated, interface,
            minimize=False,
            solvent_config=solvent_config,
        )

        assert isinstance(result, BindingEnergyResult)
        assert not np.isnan(result.dG_separated)
        assert not np.isinf(result.dG_separated)
        assert result.solvent_model == "obc2"

    @pytest.mark.slow
    def test_different_gb_models_produce_different_energies(self, pdb_1yy9: Path):
        """Test that different GB models produce different energies."""
        structure = parse_structure(pdb_1yy9)
        structure = structure.subset(["A"])

        # OBC2
        modeller1, forcefield1 = prepare_structure_for_openmm(
            structure, SolventModel.OBC2
        )
        energy_obc2 = calculate_potential_energy(modeller1, forcefield1)

        # GBN
        modeller2, forcefield2 = prepare_structure_for_openmm(
            structure, SolventModel.GBN
        )
        energy_gbn = calculate_potential_energy(modeller2, forcefield2)

        # Different models should produce different energies
        assert energy_obc2 != energy_gbn


class TestBackwardCompatibility:
    """Tests ensuring backward compatibility."""

    def test_default_solvent_config_is_vacuum(self):
        """Test that default solvent config uses vacuum model."""
        config = SolventConfig()
        assert config.model == SolventModel.VACUUM

    @pytest.mark.slow
    def test_binding_energy_default_uses_vacuum_model(self, pdb_1yy9: Path):
        """Test that binding energy without solvent_config reports vacuum model."""
        structure = parse_structure(pdb_1yy9)
        interface = detect_interface(structure, ["C", "D"], ["A"], cutoff=8.0)
        separated = separate_chains(structure, interface)

        # Without solvent_config
        result = calculate_binding_energy(
            structure, separated, interface, minimize=False
        )

        # Should report vacuum model
        assert result.solvent_model == "vacuum"
        assert result.salt_concentration == 0.15  # default

    @pytest.mark.slow
    def test_prepare_structure_default_uses_vacuum(self, pdb_1yy9: Path):
        """Test that prepare_structure_for_openmm without solvent_model uses vacuum."""
        structure = parse_structure(pdb_1yy9)
        structure = structure.subset(["A"])

        # Default call should work (vacuum)
        modeller, forcefield = prepare_structure_for_openmm(structure)

        # Should be able to calculate energy
        energy = calculate_potential_energy(modeller, forcefield)
        assert isinstance(energy, float)
        assert not np.isnan(energy)
