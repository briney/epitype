"""PDB download utilities."""
from __future__ import annotations

import re
import urllib.error
import urllib.request
from pathlib import Path


# PDB ID patterns
CLASSIC_PDB_PATTERN = re.compile(r"^[0-9][a-zA-Z0-9]{3}$")  # e.g., 1yy9, 4FQI
EXTENDED_PDB_PATTERN = re.compile(r"^pdb_[a-zA-Z0-9]{8}$", re.IGNORECASE)  # e.g., pdb_00001yy9

RCSB_DOWNLOAD_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"


class PDBDownloadError(Exception):
    """Raised when PDB download fails."""

    pass


def is_pdb_id(value: str) -> bool:
    """
    Check if a string is a valid PDB ID.

    Args:
        value: String to check

    Returns:
        True if value matches classic (4-char) or extended (pdb_XXXXXXXX) PDB ID format
    """
    return bool(CLASSIC_PDB_PATTERN.match(value) or EXTENDED_PDB_PATTERN.match(value))


def normalize_pdb_id(pdb_id: str) -> str:
    """
    Normalize PDB ID to uppercase (classic) or lowercase (extended).

    Args:
        pdb_id: Raw PDB ID string

    Returns:
        Normalized PDB ID (e.g., "1YY9" or "pdb_00001yy9")
    """
    if EXTENDED_PDB_PATTERN.match(pdb_id):
        # Extended IDs: lowercase
        return pdb_id.lower()
    # Classic IDs: uppercase
    return pdb_id.upper()


def download_pdb(pdb_id: str, output_dir: Path | None = None) -> Path:
    """
    Download PDB file from RCSB.

    Args:
        pdb_id: Valid PDB ID (classic or extended)
        output_dir: Directory to save file (defaults to current directory)

    Returns:
        Path to downloaded PDB file

    Raises:
        PDBDownloadError: If download fails
    """
    if output_dir is None:
        output_dir = Path.cwd()

    normalized_id = normalize_pdb_id(pdb_id)
    url = RCSB_DOWNLOAD_URL.format(pdb_id=normalized_id)
    output_path = output_dir / f"{normalized_id}.pdb"

    try:
        urllib.request.urlretrieve(url, output_path)
    except urllib.error.HTTPError as e:
        if e.code == 404:
            raise PDBDownloadError(
                f"PDB ID '{normalized_id}' not found in RCSB database"
            ) from e
        raise PDBDownloadError(
            f"Failed to download PDB '{normalized_id}': HTTP {e.code}"
        ) from e
    except urllib.error.URLError as e:
        raise PDBDownloadError(
            f"Network error downloading PDB '{normalized_id}': {e.reason}"
        ) from e

    return output_path
