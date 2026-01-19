"""API integration tests."""
from pathlib import Path

import pytest

import epitype
from epitype import (
    AnalysisConfig,
    InterfaceMetrics,
    InterfaceRegion,
    Structure,
    __version__,
    analyze,
    analyze_interface,
    analyze_structure,
    detect_interface,
    load_structure,
    metrics,
    parse_chain_groups,
)


class TestVersionAccessibility:
    """Tests for version accessibility."""

    def test_version_is_string(self):
        """Test that version is a string."""
        assert isinstance(__version__, str)

    def test_version_format(self):
        """Test that version follows semver format."""
        parts = __version__.split(".")
        assert len(parts) >= 2
        # First two parts should be numeric
        assert parts[0].isdigit()
        assert parts[1].isdigit()

    def test_version_accessible_from_module(self):
        """Test that version is accessible from module."""
        assert epitype.__version__ == __version__


class TestLoadStructure:
    """Tests for load_structure function."""

    def test_load_structure_returns_structure(self, pdb_1yy9: Path):
        """Test that load_structure returns a Structure object."""
        structure = load_structure(pdb_1yy9)
        assert isinstance(structure, Structure)

    def test_load_structure_has_chains(self, pdb_1yy9: Path):
        """Test that loaded structure has chains."""
        structure = load_structure(pdb_1yy9)
        assert len(structure.chains) > 0

    def test_load_structure_has_atoms(self, pdb_1yy9: Path):
        """Test that loaded structure has atoms."""
        structure = load_structure(pdb_1yy9)
        assert structure.num_atoms > 0

    def test_load_structure_has_residues(self, pdb_1yy9: Path):
        """Test that loaded structure has residues."""
        structure = load_structure(pdb_1yy9)
        assert structure.num_residues > 0


class TestDetectInterface:
    """Tests for detect_interface function."""

    def test_detect_interface_returns_region(self, pdb_1yy9: Path):
        """Test that detect_interface returns an InterfaceRegion."""
        structure = load_structure(pdb_1yy9)
        interface = detect_interface(structure, ["C", "D"], ["A"])
        assert isinstance(interface, InterfaceRegion)

    def test_detect_interface_has_residues(self, pdb_1yy9: Path):
        """Test that detected interface has residues."""
        structure = load_structure(pdb_1yy9)
        interface = detect_interface(structure, ["C", "D"], ["A"], cutoff=8.0)
        assert len(interface.residues_group1) > 0
        assert len(interface.residues_group2) > 0

    def test_detect_interface_with_custom_cutoff(self, pdb_1yy9: Path):
        """Test detect_interface with custom cutoff."""
        structure = load_structure(pdb_1yy9)
        interface_tight = detect_interface(structure, ["C", "D"], ["A"], cutoff=5.0)
        interface_loose = detect_interface(structure, ["C", "D"], ["A"], cutoff=12.0)

        # Larger cutoff should find more interface residues
        assert interface_loose.num_interface_residues >= interface_tight.num_interface_residues


class TestParseChainGroups:
    """Tests for parse_chain_groups function."""

    def test_parse_simple_spec(self):
        """Test parsing simple chain spec."""
        g1, g2 = parse_chain_groups("A_B")
        assert g1 == ["A"]
        assert g2 == ["B"]

    def test_parse_multi_chain_spec(self):
        """Test parsing multi-chain spec."""
        g1, g2 = parse_chain_groups("HL_A")
        assert g1 == ["H", "L"]
        assert g2 == ["A"]

    def test_parse_invalid_spec_raises(self):
        """Test that invalid spec raises ValueError."""
        with pytest.raises(ValueError):
            parse_chain_groups("invalid")


class TestAnalyze:
    """Tests for analyze function."""

    def test_analyze_returns_metrics(self, pdb_1yy9: Path):
        """Test that analyze returns InterfaceMetrics."""
        metrics_result = analyze(
            pdb_1yy9,
            chains="CD_A",
            compute_energy=False,
            compute_shape=False,
        )
        assert isinstance(metrics_result, InterfaceMetrics)

    def test_analyze_has_dsasa(self, pdb_1yy9: Path):
        """Test that analyze computes dSASA."""
        metrics_result = analyze(
            pdb_1yy9,
            chains="CD_A",
            compute_energy=False,
            compute_shape=False,
        )
        assert metrics_result.dSASA_int > 0

    def test_analyze_has_hbonds(self, pdb_1yy9: Path):
        """Test that analyze computes H-bonds."""
        metrics_result = analyze(
            pdb_1yy9,
            chains="CD_A",
            compute_energy=False,
            compute_shape=False,
        )
        assert metrics_result.hbonds_int >= 0

    def test_analyze_has_interface_residues(self, pdb_1yy9: Path):
        """Test that analyze returns interface residues."""
        metrics_result = analyze(
            pdb_1yy9,
            chains="CD_A",
            compute_energy=False,
            compute_shape=False,
        )
        assert len(metrics_result.interface_residues_1) > 0
        assert len(metrics_result.interface_residues_2) > 0


class TestAnalyzeStructure:
    """Tests for analyze_structure function."""

    def test_analyze_structure_returns_metrics(self, pdb_1yy9: Path):
        """Test that analyze_structure returns InterfaceMetrics."""
        structure = load_structure(pdb_1yy9)
        metrics_result = analyze_structure(
            structure,
            chain_spec="CD_A",
            compute_energy=False,
            compute_shape=False,
        )
        assert isinstance(metrics_result, InterfaceMetrics)

    def test_analyze_structure_same_as_analyze(self, pdb_1yy9: Path):
        """Test that analyze_structure gives same results as analyze."""
        structure = load_structure(pdb_1yy9)

        metrics_file = analyze(
            pdb_1yy9,
            chains="CD_A",
            compute_energy=False,
            compute_shape=False,
            compute_packstat=False,
        )
        metrics_struct = analyze_structure(
            structure,
            chain_spec="CD_A",
            compute_energy=False,
            compute_shape=False,
            compute_packstat=False,
        )

        assert abs(metrics_file.dSASA_int - metrics_struct.dSASA_int) < 0.1
        assert metrics_file.hbonds_int == metrics_struct.hbonds_int


class TestAnalyzeInterface:
    """Tests for analyze_interface function."""

    def test_analyze_interface_returns_metrics(self, pdb_1yy9: Path):
        """Test that analyze_interface returns InterfaceMetrics."""
        config = AnalysisConfig(
            chain_spec="CD_A",
            compute_energy=False,
            compute_shape=False,
        )
        metrics_result = analyze_interface(pdb_1yy9, config)
        assert isinstance(metrics_result, InterfaceMetrics)

    def test_analyze_interface_with_progress_callback(self, pdb_1yy9: Path):
        """Test analyze_interface with progress callback."""
        progress_steps = []

        def callback(step: str, progress: float):
            progress_steps.append((step, progress))

        config = AnalysisConfig(
            chain_spec="CD_A",
            compute_energy=False,
            compute_shape=False,
            compute_packstat=False,
        )
        analyze_interface(pdb_1yy9, config, progress_callback=callback)

        # Should have recorded progress
        assert len(progress_steps) > 0
        # Last step should be complete
        assert progress_steps[-1][0] == "Complete"
        assert progress_steps[-1][1] == 1.0


class TestAnalysisConfig:
    """Tests for AnalysisConfig dataclass."""

    def test_config_defaults(self):
        """Test AnalysisConfig default values."""
        config = AnalysisConfig(chain_spec="A_B")

        assert config.chain_spec == "A_B"
        assert config.interface_cutoff == 8.0
        assert config.probe_radius == 1.4
        assert config.separation_distance == 500.0
        assert config.minimize_energy is True
        assert config.compute_energy is True
        assert config.compute_shape is True
        assert config.compute_packstat is True

    def test_config_custom_values(self):
        """Test AnalysisConfig with custom values."""
        config = AnalysisConfig(
            chain_spec="HL_A",
            interface_cutoff=10.0,
            compute_energy=False,
        )

        assert config.chain_spec == "HL_A"
        assert config.interface_cutoff == 10.0
        assert config.compute_energy is False


class TestMetricsSubmoduleAccess:
    """Tests for metrics submodule access."""

    def test_sasa_module_accessible(self):
        """Test that sasa module is accessible."""
        assert hasattr(metrics, "sasa")
        assert hasattr(metrics.sasa, "calculate_sasa")

    def test_hbonds_module_accessible(self):
        """Test that hbonds module is accessible."""
        assert hasattr(metrics, "hbonds")
        assert hasattr(metrics.hbonds, "analyze_hbonds")

    def test_packstat_module_accessible(self):
        """Test that packstat module is accessible."""
        assert hasattr(metrics, "packstat")
        assert hasattr(metrics.packstat, "calculate_packstat")

    def test_shape_module_accessible(self):
        """Test that shape module is accessible."""
        assert hasattr(metrics, "shape")
        assert hasattr(metrics.shape, "calculate_shape_complementarity")

    def test_energy_module_accessible(self):
        """Test that energy module is accessible."""
        assert hasattr(metrics, "energy")
        assert hasattr(metrics.energy, "check_openmm_available")

    def test_direct_function_import(self):
        """Test that functions can be imported directly."""
        from epitype.metrics import analyze_hbonds, calculate_sasa

        assert callable(calculate_sasa)
        assert callable(analyze_hbonds)


class TestInterfaceMetrics:
    """Tests for InterfaceMetrics dataclass."""

    def test_metrics_to_dict(self, pdb_1yy9: Path):
        """Test that InterfaceMetrics can be converted to dict."""
        metrics_result = analyze(
            pdb_1yy9,
            chains="CD_A",
            compute_energy=False,
            compute_shape=False,
            compute_packstat=False,
        )
        result_dict = metrics_result.to_dict()

        assert isinstance(result_dict, dict)
        assert "dSASA_int" in result_dict
        assert "hbonds_int" in result_dict
        assert "interface_residues_1" in result_dict

    def test_metrics_attributes(self, pdb_1yy9: Path):
        """Test InterfaceMetrics attributes."""
        metrics_result = analyze(
            pdb_1yy9,
            chains="CD_A",
            compute_energy=False,
            compute_shape=False,
            compute_packstat=False,
        )

        # Check all expected attributes exist
        assert hasattr(metrics_result, "dG_separated")
        assert hasattr(metrics_result, "dSASA_int")
        assert hasattr(metrics_result, "dG_per_dSASA")
        assert hasattr(metrics_result, "sc_value")
        assert hasattr(metrics_result, "packstat")
        assert hasattr(metrics_result, "delta_unsat_hbonds")
        assert hasattr(metrics_result, "hbonds_int")
        assert hasattr(metrics_result, "interface_residues_1")
        assert hasattr(metrics_result, "interface_residues_2")
        assert hasattr(metrics_result, "sasa_complex")
        assert hasattr(metrics_result, "sasa_separated")
