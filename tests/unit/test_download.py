"""Tests for PDB download functionality."""
from __future__ import annotations

import urllib.error
from pathlib import Path
from unittest.mock import patch

import pytest

from epitype.io.download import (
    PDBDownloadError,
    download_pdb,
    is_pdb_id,
    normalize_pdb_id,
)


class TestIsPdbId:
    """Tests for is_pdb_id function."""

    @pytest.mark.parametrize(
        "pdb_id",
        [
            "1yy9",
            "4FQI",
            "1ABC",
            "9xyz",
            "1a2b",
            "2ABC",
            "3def",
        ],
    )
    def test_valid_classic_ids(self, pdb_id: str):
        """Valid classic 4-character PDB IDs."""
        assert is_pdb_id(pdb_id) is True

    @pytest.mark.parametrize(
        "pdb_id",
        [
            "pdb_00001yy9",
            "pdb_12345678",
            "pdb_abcdefgh",
            "PDB_00001YY9",  # Case insensitive
        ],
    )
    def test_valid_extended_ids(self, pdb_id: str):
        """Valid extended PDB IDs."""
        assert is_pdb_id(pdb_id) is True

    @pytest.mark.parametrize(
        "invalid",
        [
            "abc",  # Too short
            "12345",  # Too long (classic)
            "abcd",  # First char must be digit (classic)
            "pdb_1234567",  # Too short (extended)
            "pdb_123456789",  # Too long (extended)
            "pbd_00001yy9",  # Wrong prefix
            "/path/to/1yy9.pdb",  # File path
            "structure.pdb",  # File name
            "",  # Empty string
            "1yy9.pdb",  # With extension
        ],
    )
    def test_invalid_ids(self, invalid: str):
        """Invalid PDB ID formats."""
        assert is_pdb_id(invalid) is False


class TestNormalizePdbId:
    """Tests for normalize_pdb_id function."""

    def test_classic_lowercase_to_upper(self):
        """Classic IDs should be uppercased."""
        assert normalize_pdb_id("1yy9") == "1YY9"

    def test_classic_mixed_case_to_upper(self):
        """Classic IDs with mixed case should be uppercased."""
        assert normalize_pdb_id("4FqI") == "4FQI"

    def test_classic_uppercase_unchanged(self):
        """Classic IDs already uppercase should stay unchanged."""
        assert normalize_pdb_id("4FQI") == "4FQI"

    def test_extended_to_lowercase(self):
        """Extended IDs should be lowercased."""
        assert normalize_pdb_id("PDB_00001YY9") == "pdb_00001yy9"

    def test_extended_lowercase_unchanged(self):
        """Extended IDs already lowercase should stay unchanged."""
        assert normalize_pdb_id("pdb_00001yy9") == "pdb_00001yy9"


class TestDownloadPdb:
    """Tests for download_pdb function."""

    @patch("epitype.io.download.urllib.request.urlretrieve")
    def test_successful_download(self, mock_urlretrieve, tmp_path: Path):
        """Test successful PDB download."""
        # Setup mock
        mock_urlretrieve.return_value = (str(tmp_path / "1YY9.pdb"), None)

        result = download_pdb("1yy9", tmp_path)

        assert result == tmp_path / "1YY9.pdb"
        mock_urlretrieve.assert_called_once_with(
            "https://files.rcsb.org/download/1YY9.pdb", tmp_path / "1YY9.pdb"
        )

    @patch("epitype.io.download.urllib.request.urlretrieve")
    def test_download_extended_id(self, mock_urlretrieve, tmp_path: Path):
        """Test download with extended PDB ID."""
        mock_urlretrieve.return_value = (str(tmp_path / "pdb_00001yy9.pdb"), None)

        result = download_pdb("pdb_00001yy9", tmp_path)

        assert result == tmp_path / "pdb_00001yy9.pdb"
        mock_urlretrieve.assert_called_once_with(
            "https://files.rcsb.org/download/pdb_00001yy9.pdb",
            tmp_path / "pdb_00001yy9.pdb",
        )

    @patch("epitype.io.download.urllib.request.urlretrieve")
    def test_download_404_error(self, mock_urlretrieve, tmp_path: Path):
        """Test handling of 404 (PDB not found)."""
        mock_urlretrieve.side_effect = urllib.error.HTTPError(
            "url", 404, "Not Found", {}, None
        )

        with pytest.raises(PDBDownloadError, match="not found in RCSB"):
            download_pdb("9zzz", tmp_path)

    @patch("epitype.io.download.urllib.request.urlretrieve")
    def test_download_other_http_error(self, mock_urlretrieve, tmp_path: Path):
        """Test handling of other HTTP errors."""
        mock_urlretrieve.side_effect = urllib.error.HTTPError(
            "url", 500, "Server Error", {}, None
        )

        with pytest.raises(PDBDownloadError, match="HTTP 500"):
            download_pdb("1yy9", tmp_path)

    @patch("epitype.io.download.urllib.request.urlretrieve")
    def test_download_network_error(self, mock_urlretrieve, tmp_path: Path):
        """Test handling of network errors."""
        mock_urlretrieve.side_effect = urllib.error.URLError("Connection refused")

        with pytest.raises(PDBDownloadError, match="Network error"):
            download_pdb("1yy9", tmp_path)

    @patch("epitype.io.download.urllib.request.urlretrieve")
    @patch("epitype.io.download.Path.cwd")
    def test_download_default_output_dir(self, mock_cwd, mock_urlretrieve, tmp_path: Path):
        """Test download with default output directory (current working dir)."""
        mock_cwd.return_value = tmp_path
        mock_urlretrieve.return_value = (str(tmp_path / "1YY9.pdb"), None)

        result = download_pdb("1yy9")

        assert result == tmp_path / "1YY9.pdb"
        mock_urlretrieve.assert_called_once()
