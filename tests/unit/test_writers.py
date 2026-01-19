"""Tests for io/writers.py."""
import csv
import json
from pathlib import Path

import pytest

from epitype.core.types import InterfaceMetrics
from epitype.io.writers import format_metrics_table, write_csv, write_json


@pytest.fixture
def sample_metrics() -> InterfaceMetrics:
    """Create sample metrics for testing."""
    return InterfaceMetrics(
        dG_separated=-10.5,
        dSASA_int=1500.0,
        dG_per_dSASA=-0.007,
        sc_value=0.72,
        packstat=0.68,
        delta_unsat_hbonds=3,
        hbonds_int=12,
        interface_residues_1=["A:TYR:100", "A:SER:101"],
        interface_residues_2=["B:GLU:50", "B:ARG:51"],
        sasa_complex=8500.0,
        sasa_separated=10000.0,
    )


class TestWriteJson:
    """Tests for write_json function."""

    def test_write_json_creates_file(self, sample_metrics: InterfaceMetrics, tmp_path: Path):
        """Test that write_json creates a JSON file."""
        output_path = tmp_path / "output.json"
        write_json(sample_metrics, output_path)

        assert output_path.exists()

    def test_write_json_valid_json(self, sample_metrics: InterfaceMetrics, tmp_path: Path):
        """Test that output is valid JSON."""
        output_path = tmp_path / "output.json"
        write_json(sample_metrics, output_path)

        with open(output_path) as f:
            data = json.load(f)

        assert isinstance(data, dict)

    def test_write_json_contains_metrics(self, sample_metrics: InterfaceMetrics, tmp_path: Path):
        """Test that JSON contains metrics data."""
        output_path = tmp_path / "output.json"
        write_json(sample_metrics, output_path)

        with open(output_path) as f:
            data = json.load(f)

        assert "metrics" in data
        assert data["metrics"]["dSASA_int"] == 1500.0
        assert data["metrics"]["hbonds_int"] == 12

    def test_write_json_with_structure_name(
        self, sample_metrics: InterfaceMetrics, tmp_path: Path
    ):
        """Test that structure name is included."""
        output_path = tmp_path / "output.json"
        write_json(sample_metrics, output_path, structure_name="1YY9")

        with open(output_path) as f:
            data = json.load(f)

        assert data["structure"] == "1YY9"

    def test_write_json_without_structure_name(
        self, sample_metrics: InterfaceMetrics, tmp_path: Path
    ):
        """Test that structure is None when not provided."""
        output_path = tmp_path / "output.json"
        write_json(sample_metrics, output_path)

        with open(output_path) as f:
            data = json.load(f)

        assert data["structure"] is None

    def test_write_json_custom_indent(self, sample_metrics: InterfaceMetrics, tmp_path: Path):
        """Test custom indentation."""
        output_path = tmp_path / "output.json"
        write_json(sample_metrics, output_path, indent=4)

        content = output_path.read_text()
        # Check that indentation is 4 spaces (not 2)
        assert "    " in content


class TestWriteCsv:
    """Tests for write_csv function."""

    def test_write_csv_creates_file(self, sample_metrics: InterfaceMetrics, tmp_path: Path):
        """Test that write_csv creates a CSV file."""
        output_path = tmp_path / "output.csv"
        write_csv(sample_metrics, output_path)

        assert output_path.exists()

    def test_write_csv_has_header(self, sample_metrics: InterfaceMetrics, tmp_path: Path):
        """Test that CSV has header row."""
        output_path = tmp_path / "output.csv"
        write_csv(sample_metrics, output_path)

        with open(output_path) as f:
            reader = csv.reader(f)
            header = next(reader)

        assert "structure" in header
        assert "dSASA_int" in header
        assert "hbonds_int" in header

    def test_write_csv_has_data(self, sample_metrics: InterfaceMetrics, tmp_path: Path):
        """Test that CSV has data row."""
        output_path = tmp_path / "output.csv"
        write_csv(sample_metrics, output_path, structure_name="test")

        with open(output_path) as f:
            reader = csv.DictReader(f)
            row = next(reader)

        assert row["structure"] == "test"
        assert float(row["dSASA_int"]) == 1500.0

    def test_write_csv_append_mode(self, sample_metrics: InterfaceMetrics, tmp_path: Path):
        """Test that append mode works."""
        output_path = tmp_path / "output.csv"

        # Write first row
        write_csv(sample_metrics, output_path, structure_name="struct1")

        # Append second row
        write_csv(sample_metrics, output_path, structure_name="struct2", append=True)

        with open(output_path) as f:
            reader = csv.DictReader(f)
            rows = list(reader)

        assert len(rows) == 2
        assert rows[0]["structure"] == "struct1"
        assert rows[1]["structure"] == "struct2"

    def test_write_csv_append_no_duplicate_header(
        self, sample_metrics: InterfaceMetrics, tmp_path: Path
    ):
        """Test that append mode doesn't add duplicate header."""
        output_path = tmp_path / "output.csv"

        write_csv(sample_metrics, output_path, structure_name="struct1")
        write_csv(sample_metrics, output_path, structure_name="struct2", append=True)

        lines = output_path.read_text().strip().split("\n")
        # Should have header + 2 data rows = 3 lines
        assert len(lines) == 3


class TestFormatMetricsTable:
    """Tests for format_metrics_table function."""

    def test_format_metrics_table_returns_string(self, sample_metrics: InterfaceMetrics):
        """Test that format_metrics_table returns a string."""
        result = format_metrics_table(sample_metrics)
        assert isinstance(result, str)

    def test_format_metrics_table_contains_metrics(self, sample_metrics: InterfaceMetrics):
        """Test that table contains metric values."""
        result = format_metrics_table(sample_metrics)

        assert "dG_separated" in result
        assert "dSASA_int" in result
        assert "Shape complementarity" in result
        assert "Packstat" in result
        assert "H-bonds" in result

    def test_format_metrics_table_with_structure_name(self, sample_metrics: InterfaceMetrics):
        """Test that structure name is displayed."""
        result = format_metrics_table(sample_metrics, structure_name="1YY9")

        assert "Structure: 1YY9" in result
        assert "=" * 50 in result

    def test_format_metrics_table_quality_indicators(self, sample_metrics: InterfaceMetrics):
        """Test that quality indicators are present."""
        result = format_metrics_table(sample_metrics)

        # With dG_separated = -10.5, should be "good"
        # With sc_value = 0.72, should be "good"
        # With packstat = 0.68, should be "good"
        assert "good" in result

    def test_format_metrics_table_poor_quality(self):
        """Test quality indicators for poor metrics."""
        poor_metrics = InterfaceMetrics(
            dG_separated=10.0,  # Positive = poor
            dSASA_int=1500.0,
            dG_per_dSASA=0.5,  # Positive = poor
            sc_value=0.3,  # Low = poor
            packstat=0.3,  # Low = poor
            delta_unsat_hbonds=3,
            hbonds_int=12,
            interface_residues_1=[],
            interface_residues_2=[],
            sasa_complex=8500.0,
            sasa_separated=10000.0,
        )

        result = format_metrics_table(poor_metrics)
        assert "poor" in result

    def test_format_metrics_table_multiline(self, sample_metrics: InterfaceMetrics):
        """Test that output is multiline."""
        result = format_metrics_table(sample_metrics)
        lines = result.split("\n")

        assert len(lines) >= 7  # Header row + separator + at least 5 metric rows
