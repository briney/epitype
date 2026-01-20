"""Integration tests for CLI."""
import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

from epitype.cli.main import app

runner = CliRunner()


class TestVersionCommand:
    """Tests for the version command."""

    def test_version_command(self):
        """Test that version command outputs version."""
        result = runner.invoke(app, ["version"])

        assert result.exit_code == 0
        assert "epitype" in result.stdout
        assert "v" in result.stdout


class TestRunCommand:
    """Tests for the run command."""

    def test_run_basic(self, pdb_1yy9: Path):
        """Test basic run command."""
        result = runner.invoke(
            app,
            [
                "run",
                str(pdb_1yy9),
                "--chains",
                "CD_A",
                "--no-energy",
                "--no-shape",
                "--no-packstat",
                "--quiet",
            ],
        )

        assert result.exit_code == 0
        assert "dSASA_int" in result.stdout

    def test_run_json_output(self, pdb_1yy9: Path, tmp_path: Path):
        """Test run command with JSON output."""
        output_path = tmp_path / "output.json"

        result = runner.invoke(
            app,
            [
                "run",
                str(pdb_1yy9),
                "--chains",
                "CD_A",
                "--no-energy",
                "--no-shape",
                "--no-packstat",
                "--quiet",
                "--output",
                str(output_path),
            ],
        )

        assert result.exit_code == 0
        assert output_path.exists()

        with open(output_path) as f:
            data = json.load(f)

        assert "metrics" in data
        assert "dSASA_int" in data["metrics"]

    def test_run_csv_output(self, pdb_1yy9: Path, tmp_path: Path):
        """Test run command with CSV output."""
        output_path = tmp_path / "output.csv"

        result = runner.invoke(
            app,
            [
                "run",
                str(pdb_1yy9),
                "--chains",
                "CD_A",
                "--no-energy",
                "--no-shape",
                "--no-packstat",
                "--quiet",
                "--output",
                str(output_path),
            ],
        )

        assert result.exit_code == 0
        assert output_path.exists()

        content = output_path.read_text()
        assert "dSASA_int" in content

    def test_run_invalid_chain_spec(self, pdb_1yy9: Path):
        """Test run command with invalid chain spec."""
        result = runner.invoke(
            app,
            ["run", str(pdb_1yy9), "--chains", "invalid", "--quiet"],
        )

        assert result.exit_code == 1
        assert "Error" in result.stdout

    def test_run_nonexistent_chains(self, pdb_1yy9: Path):
        """Test run command with non-existent chains."""
        result = runner.invoke(
            app,
            ["run", str(pdb_1yy9), "--chains", "XY_Z", "--quiet"],
        )

        assert result.exit_code == 1


class TestSasaCommand:
    """Tests for the sasa command."""

    def test_sasa_basic(self, pdb_1yy9: Path):
        """Test basic SASA command."""
        result = runner.invoke(app, ["sasa", str(pdb_1yy9)])

        assert result.exit_code == 0
        assert "Total SASA" in result.stdout
        assert "Polar SASA" in result.stdout
        assert "Apolar SASA" in result.stdout

    def test_sasa_json_output(self, pdb_1yy9: Path, tmp_path: Path):
        """Test SASA command with JSON output."""
        output_path = tmp_path / "sasa.json"

        result = runner.invoke(
            app,
            ["sasa", str(pdb_1yy9), "--output", str(output_path)],
        )

        assert result.exit_code == 0
        assert output_path.exists()

        with open(output_path) as f:
            data = json.load(f)

        assert "total" in data
        assert "polar" in data
        assert "apolar" in data

    def test_sasa_per_residue(self, pdb_1yy9: Path):
        """Test SASA command with per-residue breakdown."""
        result = runner.invoke(app, ["sasa", str(pdb_1yy9), "--per-residue"])

        assert result.exit_code == 0
        assert "Per-residue SASA" in result.stdout


class TestHbondsCommand:
    """Tests for the hbonds command."""

    def test_hbonds_basic(self, pdb_1yy9: Path):
        """Test basic H-bonds command."""
        result = runner.invoke(
            app,
            ["hbonds", str(pdb_1yy9), "--chains", "CD_A"],
        )

        assert result.exit_code == 0
        assert "Interface H-bonds" in result.stdout

    def test_hbonds_json_output(self, pdb_1yy9: Path, tmp_path: Path):
        """Test H-bonds command with JSON output."""
        output_path = tmp_path / "hbonds.json"

        result = runner.invoke(
            app,
            ["hbonds", str(pdb_1yy9), "--chains", "CD_A", "--output", str(output_path)],
        )

        assert result.exit_code == 0
        assert output_path.exists()

        with open(output_path) as f:
            data = json.load(f)

        assert "interface_hbonds" in data
        assert "total_hbonds" in data


class TestPackstatCommand:
    """Tests for the packstat command."""

    def test_packstat_basic(self, pdb_1yy9: Path):
        """Test basic PackStat command."""
        result = runner.invoke(
            app,
            ["packstat", str(pdb_1yy9), "--chains", "CD_A"],
        )

        assert result.exit_code == 0
        assert "PackStat" in result.stdout

    def test_packstat_json_output(self, pdb_1yy9: Path, tmp_path: Path):
        """Test PackStat command with JSON output."""
        output_path = tmp_path / "packstat.json"

        result = runner.invoke(
            app,
            ["packstat", str(pdb_1yy9), "--chains", "CD_A", "--output", str(output_path)],
        )

        assert result.exit_code == 0
        assert output_path.exists()

        with open(output_path) as f:
            data = json.load(f)

        assert "packstat" in data
        assert "packstat_interface" in data


class TestHelpText:
    """Tests for help text."""

    def test_main_help(self):
        """Test main help text."""
        result = runner.invoke(app, ["--help"])

        assert result.exit_code == 0
        assert "run" in result.stdout
        assert "sasa" in result.stdout
        assert "hbonds" in result.stdout
        assert "packstat" in result.stdout

    def test_run_help(self):
        """Test run command help."""
        result = runner.invoke(app, ["run", "--help"])

        assert result.exit_code == 0
        assert "--chains" in result.stdout
        assert "--output" in result.stdout
        assert "--no-energy" in result.stdout

    def test_sasa_help(self):
        """Test sasa command help."""
        result = runner.invoke(app, ["sasa", "--help"])

        assert result.exit_code == 0
        assert "--probe-radius" in result.stdout
        assert "--per-residue" in result.stdout


class TestShapeCommand:
    """Tests for the shape command."""

    def test_shape_basic(self, pdb_1yy9: Path):
        """Test basic shape command."""
        from epitype.surface.nanoshaper import check_nanoshaper_available

        if not check_nanoshaper_available():
            pytest.skip("NanoShaper not available")

        result = runner.invoke(
            app,
            ["shape", str(pdb_1yy9), "--chains", "CD_A"],
        )

        assert result.exit_code == 0
        assert "Shape Complementarity" in result.stdout
        assert "Sc (overall)" in result.stdout

    def test_shape_json_output(self, pdb_1yy9: Path, tmp_path: Path):
        """Test shape command with JSON output."""
        from epitype.surface.nanoshaper import check_nanoshaper_available

        if not check_nanoshaper_available():
            pytest.skip("NanoShaper not available")

        output_path = tmp_path / "shape.json"

        result = runner.invoke(
            app,
            ["shape", str(pdb_1yy9), "--chains", "CD_A", "--output", str(output_path)],
        )

        assert result.exit_code == 0
        assert output_path.exists()

        with open(output_path) as f:
            data = json.load(f)

        assert "sc_value" in data
        assert "sc_median" in data
        assert "sc_group1" in data
        assert "sc_group2" in data
        assert "interface_area" in data

    def test_shape_custom_cutoff(self, pdb_1yy9: Path):
        """Test shape command with custom cutoff."""
        from epitype.surface.nanoshaper import check_nanoshaper_available

        if not check_nanoshaper_available():
            pytest.skip("NanoShaper not available")

        result = runner.invoke(
            app,
            ["shape", str(pdb_1yy9), "--chains", "CD_A", "--cutoff", "10.0"],
        )

        assert result.exit_code == 0
        assert "Sc (overall)" in result.stdout

    def test_shape_invalid_chains(self, pdb_1yy9: Path):
        """Test shape command with invalid chains."""
        from epitype.surface.nanoshaper import check_nanoshaper_available

        if not check_nanoshaper_available():
            pytest.skip("NanoShaper not available")

        result = runner.invoke(
            app,
            ["shape", str(pdb_1yy9), "--chains", "invalid"],
        )

        assert result.exit_code == 1

    def test_shape_nonexistent_chains(self, pdb_1yy9: Path):
        """Test shape command with non-existent chains."""
        from epitype.surface.nanoshaper import check_nanoshaper_available

        if not check_nanoshaper_available():
            pytest.skip("NanoShaper not available")

        result = runner.invoke(
            app,
            ["shape", str(pdb_1yy9), "--chains", "XY_Z"],
        )

        assert result.exit_code == 1

    def test_shape_no_interface(self, pdb_1yy9: Path):
        """Test shape command when no interface detected."""
        from epitype.surface.nanoshaper import check_nanoshaper_available

        if not check_nanoshaper_available():
            pytest.skip("NanoShaper not available")

        result = runner.invoke(
            app,
            ["shape", str(pdb_1yy9), "--chains", "CD_A", "--cutoff", "0.1"],
        )

        assert result.exit_code == 1
        assert "No interface detected" in result.stdout

    def test_shape_displays_quality_indicator(self, pdb_1yy9: Path):
        """Test that shape command displays quality indicator."""
        from epitype.surface.nanoshaper import check_nanoshaper_available

        if not check_nanoshaper_available():
            pytest.skip("NanoShaper not available")

        result = runner.invoke(
            app,
            ["shape", str(pdb_1yy9), "--chains", "CD_A"],
        )

        assert result.exit_code == 0
        # Should have one of the quality indicators
        assert any(q in result.stdout for q in ["good", "moderate", "poor"])

    def test_shape_shows_nanoshaper_source(self, pdb_1yy9: Path):
        """Test that shape command shows NanoShaper source."""
        from epitype.surface.nanoshaper import check_nanoshaper_available

        if not check_nanoshaper_available():
            pytest.skip("NanoShaper not available")

        result = runner.invoke(
            app,
            ["shape", str(pdb_1yy9), "--chains", "CD_A"],
        )

        assert result.exit_code == 0
        assert "Using NanoShaper" in result.stdout
        # Should show source (bundled, environment, or system)
        assert any(s in result.stdout for s in ["bundled", "environment", "system"])

    def test_shape_help(self):
        """Test shape command help."""
        result = runner.invoke(app, ["shape", "--help"])

        assert result.exit_code == 0
        assert "--chains" in result.stdout
        assert "--cutoff" in result.stdout
        assert "--output" in result.stdout


class TestShapeCommandWithoutNanoshaper:
    """Tests for shape command when NanoShaper is unavailable."""

    def test_shape_unavailable_message(self, pdb_1yy9: Path):
        """Test appropriate error when NanoShaper unavailable."""
        from unittest.mock import patch

        with patch("epitype.bin.get_nanoshaper_info") as mock_info:
            mock_info.return_value = {
                "path": "NanoShaper",
                "source": "system",
                "available": False,
                "platform": "linux-x86_64",
            }

            result = runner.invoke(
                app,
                ["shape", str(pdb_1yy9), "--chains", "CD_A"],
            )

            assert result.exit_code == 1
            assert "NanoShaper is not available" in result.stdout
