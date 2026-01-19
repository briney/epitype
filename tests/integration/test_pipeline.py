"""Integration tests for analysis pipeline."""
from pathlib import Path

import pytest

from epitype import InterfaceMetrics, load_structure
from epitype.analysis.pipeline import (
    AnalysisConfig,
    analyze,
    analyze_interface,
    analyze_structure,
)


class TestAnalyzePipeline:
    """Tests for the analyze() function."""

    def test_analyze_basic(self, pdb_1yy9: Path):
        """Test basic analysis with minimal options."""
        metrics = analyze(
            pdb_1yy9,
            chains="CD_A",
            compute_energy=False,
            compute_shape=False,
            compute_packstat=False,
        )

        assert isinstance(metrics, InterfaceMetrics)
        assert metrics.dSASA_int > 0
        assert metrics.hbonds_int >= 0

    def test_analyze_returns_interface_residues(self, pdb_1yy9: Path):
        """Test that interface residues are returned."""
        metrics = analyze(
            pdb_1yy9,
            chains="CD_A",
            compute_energy=False,
            compute_shape=False,
            compute_packstat=False,
        )

        assert len(metrics.interface_residues_1) > 0
        assert len(metrics.interface_residues_2) > 0

    def test_analyze_with_packstat(self, pdb_1yy9: Path):
        """Test analysis with packstat enabled."""
        metrics = analyze(
            pdb_1yy9,
            chains="CD_A",
            compute_energy=False,
            compute_shape=False,
            compute_packstat=True,
        )

        assert 0.0 <= metrics.packstat <= 1.0

    def test_analyze_different_cutoffs(self, pdb_1yy9: Path):
        """Test that different cutoffs produce different results."""
        metrics_tight = analyze(
            pdb_1yy9,
            chains="CD_A",
            interface_cutoff=5.0,
            compute_energy=False,
            compute_shape=False,
            compute_packstat=False,
        )
        metrics_loose = analyze(
            pdb_1yy9,
            chains="CD_A",
            interface_cutoff=12.0,
            compute_energy=False,
            compute_shape=False,
            compute_packstat=False,
        )

        # Looser cutoff should find more buried surface
        assert metrics_loose.dSASA_int >= metrics_tight.dSASA_int


class TestAnalyzeInterface:
    """Tests for the analyze_interface() function with config."""

    def test_analyze_interface_with_config(self, pdb_1yy9: Path):
        """Test analysis with AnalysisConfig."""
        config = AnalysisConfig(
            chain_spec="CD_A",
            interface_cutoff=8.0,
            compute_energy=False,
            compute_shape=False,
            compute_packstat=False,
        )

        metrics = analyze_interface(pdb_1yy9, config)

        assert isinstance(metrics, InterfaceMetrics)
        assert metrics.dSASA_int > 0

    def test_analyze_interface_progress_callback(self, pdb_1yy9: Path):
        """Test that progress callback is called."""
        progress_calls = []

        def progress_callback(step: str, fraction: float):
            progress_calls.append((step, fraction))

        config = AnalysisConfig(
            chain_spec="CD_A",
            compute_energy=False,
            compute_shape=False,
            compute_packstat=False,
        )

        analyze_interface(pdb_1yy9, config, progress_callback=progress_callback)

        # Should have multiple progress updates
        assert len(progress_calls) > 0

        # Last call should indicate completion
        assert progress_calls[-1][0] == "Complete"
        assert progress_calls[-1][1] == 1.0

    def test_analyze_interface_progress_steps(self, pdb_1yy9: Path):
        """Test that expected progress steps are reported."""
        steps_seen = set()

        def progress_callback(step: str, fraction: float):
            steps_seen.add(step)

        config = AnalysisConfig(
            chain_spec="CD_A",
            compute_energy=False,
            compute_shape=False,
            compute_packstat=False,
        )

        analyze_interface(pdb_1yy9, config, progress_callback=progress_callback)

        # Should see these steps
        assert "Parsing structure" in steps_seen
        assert "Detecting interface" in steps_seen
        assert "Complete" in steps_seen


class TestAnalyzeStructure:
    """Tests for the analyze_structure() function."""

    def test_analyze_structure_basic(self, pdb_1yy9: Path):
        """Test analysis with pre-loaded structure."""
        structure = load_structure(pdb_1yy9)

        metrics = analyze_structure(
            structure,
            chain_spec="CD_A",
            compute_energy=False,
            compute_shape=False,
            compute_packstat=False,
        )

        assert isinstance(metrics, InterfaceMetrics)
        assert metrics.dSASA_int > 0

    def test_analyze_structure_reusable(self, pdb_1yy9: Path):
        """Test that structure can be analyzed multiple times."""
        structure = load_structure(pdb_1yy9)

        # Analyze twice with same structure
        metrics1 = analyze_structure(
            structure,
            chain_spec="CD_A",
            compute_energy=False,
            compute_shape=False,
            compute_packstat=False,
        )
        metrics2 = analyze_structure(
            structure,
            chain_spec="CD_A",
            compute_energy=False,
            compute_shape=False,
            compute_packstat=False,
        )

        # Results should be identical
        assert metrics1.dSASA_int == metrics2.dSASA_int
        assert metrics1.hbonds_int == metrics2.hbonds_int

    def test_analyze_structure_matches_file_analysis(self, pdb_1yy9: Path):
        """Test that structure analysis matches file analysis."""
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

        # Results should match within floating point tolerance
        assert abs(metrics_file.dSASA_int - metrics_struct.dSASA_int) < 0.1
        assert metrics_file.hbonds_int == metrics_struct.hbonds_int


class TestAnalysisConfig:
    """Tests for AnalysisConfig dataclass."""

    def test_config_defaults(self):
        """Test default configuration values."""
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
        """Test custom configuration values."""
        config = AnalysisConfig(
            chain_spec="HL_A",
            interface_cutoff=10.0,
            probe_radius=1.5,
            separation_distance=1000.0,
            minimize_energy=False,
            compute_energy=False,
            compute_shape=False,
            compute_packstat=False,
        )

        assert config.chain_spec == "HL_A"
        assert config.interface_cutoff == 10.0
        assert config.probe_radius == 1.5
        assert config.separation_distance == 1000.0
        assert config.minimize_energy is False
        assert config.compute_energy is False
        assert config.compute_shape is False
        assert config.compute_packstat is False


class TestErrorHandling:
    """Tests for error handling in pipeline."""

    def test_invalid_chain_spec_raises(self, pdb_1yy9: Path):
        """Test that invalid chain spec raises ValueError."""
        with pytest.raises(ValueError):
            analyze(pdb_1yy9, chains="invalid_spec")

    def test_nonexistent_chains_raises(self, pdb_1yy9: Path):
        """Test that non-existent chains raise ValueError."""
        with pytest.raises(ValueError):
            analyze(pdb_1yy9, chains="XY_Z")

    def test_no_interface_raises(self, pdb_1yy9: Path):
        """Test that no interface detection raises ValueError."""
        # Using very small cutoff should find no interface
        with pytest.raises(ValueError, match="No interface detected"):
            analyze(
                pdb_1yy9,
                chains="CD_A",
                interface_cutoff=0.1,  # Very small, no contacts
            )

    def test_nonexistent_file_raises(self, tmp_path: Path):
        """Test that non-existent file raises FileNotFoundError."""
        fake_path = tmp_path / "nonexistent.pdb"

        with pytest.raises(FileNotFoundError):
            analyze(fake_path, chains="A_B")


class TestSelectiveMetrics:
    """Tests for running analysis with selective metrics."""

    def test_disable_all_optional_metrics(self, pdb_1yy9: Path):
        """Test analysis with all optional metrics disabled."""
        metrics = analyze(
            pdb_1yy9,
            chains="CD_A",
            compute_energy=False,
            compute_shape=False,
            compute_packstat=False,
        )

        # Should still have core metrics
        assert metrics.dSASA_int > 0
        assert metrics.hbonds_int >= 0

        # Optional metrics should be default/zero
        assert metrics.dG_separated == 0.0
        assert metrics.sc_value == 0.0
        assert metrics.packstat == 0.0
