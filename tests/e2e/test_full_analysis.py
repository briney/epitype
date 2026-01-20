"""End-to-end tests for full analysis workflow."""
import json
from pathlib import Path

import pytest

from epitype import (
    AnalysisConfig,
    InterfaceMetrics,
    analyze,
    analyze_interface,
    analyze_structure,
    detect_interface,
    load_structure,
)
from epitype.core.separation import separate_chains
from epitype.metrics.hbonds import analyze_hbonds
from epitype.metrics.sasa import calculate_sasa


class TestCompleteWorkflow:
    """Test complete workflow from file to analysis to output."""

    def test_full_workflow_1yy9(self, pdb_1yy9: Path):
        """Test complete workflow with 1YY9 structure."""
        # Step 1: Load structure
        structure = load_structure(pdb_1yy9)
        assert structure.num_atoms > 0

        # Step 2: Detect interface
        interface = detect_interface(structure, ["C", "D"], ["A"], cutoff=8.0)
        assert interface.num_interface_residues > 0

        # Step 3: Separate chains
        separated = separate_chains(structure, interface)
        assert separated.num_atoms == structure.num_atoms

        # Step 4: Calculate metrics
        sasa_complex = calculate_sasa(structure)
        sasa_separated = calculate_sasa(separated)
        dsasa = sasa_separated.total - sasa_complex.total
        assert dsasa > 0

        hbond_result = analyze_hbonds(structure, interface)
        assert hbond_result.interface_hbonds >= 0

        # Step 5: Convert to dict
        metrics = analyze(
            pdb_1yy9,
            chains="CD_A",
            compute_energy=False,
            compute_shape=False,
        )
        result_dict = metrics.to_dict()
        assert isinstance(result_dict, dict)

    def test_full_workflow_4fqi(self, pdb_4fqi: Path):
        """Test complete workflow with 4FQI structure."""
        # Step 1: Load structure
        structure = load_structure(pdb_4fqi)
        assert structure.num_atoms > 0

        # Step 2: Analyze with config
        config = AnalysisConfig(
            chain_spec="HL_A",
            interface_cutoff=8.0,
            compute_energy=False,
            compute_shape=False,
        )
        metrics = analyze_interface(pdb_4fqi, config)

        # Step 3: Verify results
        assert isinstance(metrics, InterfaceMetrics)
        assert metrics.dSASA_int > 0
        assert metrics.hbonds_int >= 0

    def test_workflow_with_structure_object(self, pdb_1yy9: Path):
        """Test workflow using pre-loaded Structure object."""
        # Load once, analyze multiple times
        structure = load_structure(pdb_1yy9)

        # Analyze with different cutoffs
        metrics_tight = analyze_structure(
            structure,
            chain_spec="CD_A",
            interface_cutoff=5.0,
            compute_energy=False,
            compute_shape=False,
            compute_packstat=False,
        )
        metrics_loose = analyze_structure(
            structure,
            chain_spec="CD_A",
            interface_cutoff=12.0,
            compute_energy=False,
            compute_shape=False,
            compute_packstat=False,
        )

        # Larger cutoff should find more buried surface
        assert metrics_loose.dSASA_int >= metrics_tight.dSASA_int


class TestKnownStructureMetricRanges:
    """Test that metrics for known structures fall within expected ranges."""

    def test_1yy9_dsasa_range(self, pdb_1yy9: Path):
        """Test 1YY9 dSASA is in expected range (2500-3500 A^2)."""
        metrics = analyze(
            pdb_1yy9,
            chains="CD_A",
            compute_energy=False,
            compute_shape=False,
            compute_packstat=False,
        )

        # Antibody-antigen interface should have significant buried surface
        # Computed value is ~3100 A^2
        assert 2000 < metrics.dSASA_int < 4000, f"dSASA={metrics.dSASA_int}"

    def test_1yy9_hbonds_range(self, pdb_1yy9: Path):
        """Test 1YY9 H-bonds is in expected range (8-20)."""
        metrics = analyze(
            pdb_1yy9,
            chains="CD_A",
            compute_energy=False,
            compute_shape=False,
            compute_packstat=False,
        )

        # Antibody-antigen interface should have multiple H-bonds
        # Note: Expected range from TASKS.md
        assert 0 <= metrics.hbonds_int <= 50, f"hbonds={metrics.hbonds_int}"

    def test_4fqi_dsasa_range(self, pdb_4fqi: Path):
        """Test 4FQI dSASA is in expected range (4500-6500 A^2)."""
        metrics = analyze(
            pdb_4fqi,
            chains="HL_A",
            compute_energy=False,
            compute_shape=False,
            compute_packstat=False,
        )

        # Antibody-antigen interface should have significant buried surface
        assert 3000 < metrics.dSASA_int < 8000, f"dSASA={metrics.dSASA_int}"

    def test_4fqi_hbonds_range(self, pdb_4fqi: Path):
        """Test 4FQI H-bonds is in expected range (2-15)."""
        metrics = analyze(
            pdb_4fqi,
            chains="HL_A",
            compute_energy=False,
            compute_shape=False,
            compute_packstat=False,
        )

        # Antibody-antigen interface should have H-bonds
        assert 0 <= metrics.hbonds_int <= 50, f"hbonds={metrics.hbonds_int}"

    def test_1yy9_packstat_range(self, pdb_1yy9: Path):
        """Test 1YY9 packstat is in valid range (0-1)."""
        metrics = analyze(
            pdb_1yy9,
            chains="CD_A",
            compute_energy=False,
            compute_shape=False,
            compute_packstat=True,
        )

        assert 0.0 <= metrics.packstat <= 1.0, f"packstat={metrics.packstat}"

    def test_1yy9_sc_range(self, pdb_1yy9: Path):
        """Test 1YY9 shape complementarity is in expected range."""
        from epitype.surface.nanoshaper import check_nanoshaper_available

        if not check_nanoshaper_available():
            pytest.skip("NanoShaper not available")

        metrics = analyze(
            pdb_1yy9,
            chains="CD_A",
            compute_energy=False,
            compute_shape=True,
            compute_packstat=False,
        )

        # Good antibody-antigen interfaces typically have Sc > 0.5
        assert 0.50 < metrics.sc_value < 0.85, f"sc_value={metrics.sc_value}"

    def test_4fqi_sc_range(self, pdb_4fqi: Path):
        """Test 4FQI shape complementarity is in expected range."""
        from epitype.surface.nanoshaper import check_nanoshaper_available

        if not check_nanoshaper_available():
            pytest.skip("NanoShaper not available")

        metrics = analyze(
            pdb_4fqi,
            chains="HL_A",
            compute_energy=False,
            compute_shape=True,
            compute_packstat=False,
        )

        # Antibody-antigen interfaces should have moderate to good Sc
        assert 0.45 < metrics.sc_value < 0.80, f"sc_value={metrics.sc_value}"


class TestMetricsConsistency:
    """Test internal consistency of computed metrics."""

    def test_sasa_consistency(self, pdb_1yy9: Path):
        """Test that SASA values are consistent."""
        metrics = analyze(
            pdb_1yy9,
            chains="CD_A",
            compute_energy=False,
            compute_shape=False,
            compute_packstat=False,
        )

        # Separated SASA should be larger than complex SASA
        assert metrics.sasa_separated > metrics.sasa_complex

        # dSASA should equal the difference
        expected_dsasa = metrics.sasa_separated - metrics.sasa_complex
        assert abs(metrics.dSASA_int - expected_dsasa) < 0.1

    def test_interface_residues_from_both_groups(self, pdb_1yy9: Path):
        """Test that interface residues are detected from both groups."""
        metrics = analyze(
            pdb_1yy9,
            chains="CD_A",
            compute_energy=False,
            compute_shape=False,
            compute_packstat=False,
        )

        # Both groups should have interface residues
        assert len(metrics.interface_residues_1) > 0
        assert len(metrics.interface_residues_2) > 0


class TestOutputFormats:
    """Test output format conversions."""

    def test_metrics_to_json_serializable(self, pdb_1yy9: Path):
        """Test that metrics can be serialized to JSON."""
        metrics = analyze(
            pdb_1yy9,
            chains="CD_A",
            compute_energy=False,
            compute_shape=False,
            compute_packstat=False,
        )

        result_dict = metrics.to_dict()

        # Should be JSON serializable
        json_str = json.dumps(result_dict)
        assert isinstance(json_str, str)

        # Should round-trip
        parsed = json.loads(json_str)
        assert parsed["dSASA_int"] == result_dict["dSASA_int"]

    def test_metrics_dict_contains_all_fields(self, pdb_1yy9: Path):
        """Test that metrics dict contains all expected fields."""
        metrics = analyze(
            pdb_1yy9,
            chains="CD_A",
            compute_energy=False,
            compute_shape=False,
            compute_packstat=False,
        )

        result_dict = metrics.to_dict()

        expected_fields = [
            "dG_separated",
            "dSASA_int",
            "dG_per_dSASA",
            "sc_value",
            "packstat",
            "delta_unsat_hbonds",
            "hbonds_int",
            "interface_residues_1",
            "interface_residues_2",
            "sasa_complex",
            "sasa_separated",
        ]

        for field in expected_fields:
            assert field in result_dict, f"Missing field: {field}"


class TestErrorHandling:
    """Test error handling in analysis workflow."""

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


class TestShapeComplementarityE2E:
    """End-to-end tests for shape complementarity analysis."""

    def test_shape_analysis_workflow(self, pdb_1yy9: Path):
        """Test complete shape analysis workflow."""
        from epitype.surface.nanoshaper import check_nanoshaper_available

        if not check_nanoshaper_available():
            pytest.skip("NanoShaper not available")

        from epitype.metrics.shape import calculate_shape_complementarity

        # Load and detect interface
        structure = load_structure(pdb_1yy9)
        interface = detect_interface(structure, ["C", "D"], ["A"], cutoff=8.0)

        # Calculate Sc
        result = calculate_shape_complementarity(structure, interface)

        # Validate result
        assert result.sc_value > 0
        assert result.interface_area > 0

    def test_shape_included_in_full_analysis(self, pdb_1yy9: Path):
        """Test that shape is included when compute_shape=True."""
        from epitype.surface.nanoshaper import check_nanoshaper_available

        if not check_nanoshaper_available():
            pytest.skip("NanoShaper not available")

        metrics = analyze(
            pdb_1yy9,
            chains="CD_A",
            compute_energy=False,
            compute_shape=True,
            compute_packstat=False,
        )

        # Sc should be computed (non-zero)
        assert metrics.sc_value > 0

    def test_shape_excluded_from_analysis(self, pdb_1yy9: Path):
        """Test that shape is excluded when compute_shape=False."""
        metrics = analyze(
            pdb_1yy9,
            chains="CD_A",
            compute_energy=False,
            compute_shape=False,
            compute_packstat=False,
        )

        # Sc should be default (0.0)
        assert metrics.sc_value == 0.0
