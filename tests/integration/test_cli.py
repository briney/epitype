"""Integration tests for CLI."""
import json
from pathlib import Path

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
