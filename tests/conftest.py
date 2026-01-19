"""Pytest fixtures and test configuration."""
import shutil
import tempfile
from collections.abc import Generator
from pathlib import Path

import numpy as np
import pytest

from epitype.core.interface import InterfaceRegion
from epitype.core.structure import Atom, Chain, Residue, Structure

# Test data directory
TEST_DATA_DIR = Path(__file__).parent / "data"


@pytest.fixture(scope="session")
def test_data_dir() -> Path:
    """Path to test data directory."""
    return TEST_DATA_DIR


@pytest.fixture(scope="session")
def pdb_1yy9(test_data_dir: Path) -> Path:
    """Path to 1YY9 PDB file (antibody-antigen complex)."""
    return test_data_dir / "1YY9.pdb"


@pytest.fixture(scope="session")
def pdb_4fqi(test_data_dir: Path) -> Path:
    """Path to 4FQI PDB file (antibody-antigen complex)."""
    return test_data_dir / "4FQI.pdb"


@pytest.fixture
def temp_dir() -> Generator[Path, None, None]:
    """Temporary directory for test outputs."""
    tmpdir = Path(tempfile.mkdtemp())
    yield tmpdir
    shutil.rmtree(tmpdir, ignore_errors=True)


@pytest.fixture
def minimal_structure() -> Structure:
    """Create minimal structure for unit tests."""
    structure = Structure(name="minimal")

    # Chain A - 3 residues
    chain_a = Chain(chain_id="A")
    for i in range(3):
        residue = Residue(
            index=i,
            name="ALA",
            chain_id="A",
            seq_num=i + 1,
        )
        # Add backbone atoms
        for j, (name, element, offset) in enumerate([
            ("N", "N", np.array([0.0, 0.0, 0.0])),
            ("CA", "C", np.array([1.5, 0.0, 0.0])),
            ("C", "C", np.array([2.5, 1.2, 0.0])),
            ("O", "O", np.array([2.5, 2.4, 0.0])),
            ("CB", "C", np.array([1.5, -1.5, 0.0])),
        ]):
            atom = Atom(
                index=i * 5 + j,
                name=name,
                element=element,
                coords=offset + np.array([i * 3.8, 0.0, 0.0]),
                residue_index=i,
                chain_id="A",
            )
            residue.atoms.append(atom)
        chain_a.residues.append(residue)

    # Chain B - 3 residues (offset in z)
    chain_b = Chain(chain_id="B")
    for i in range(3):
        residue = Residue(
            index=i + 3,
            name="ALA",
            chain_id="B",
            seq_num=i + 1,
        )
        for j, (name, element, offset) in enumerate([
            ("N", "N", np.array([0.0, 0.0, 0.0])),
            ("CA", "C", np.array([1.5, 0.0, 0.0])),
            ("C", "C", np.array([2.5, 1.2, 0.0])),
            ("O", "O", np.array([2.5, 2.4, 0.0])),
            ("CB", "C", np.array([1.5, -1.5, 0.0])),
        ]):
            atom = Atom(
                index=15 + i * 5 + j,
                name=name,
                element=element,
                coords=offset + np.array([i * 3.8, 0.0, 5.0]),  # 5 Å away in z
                residue_index=i + 3,
                chain_id="B",
            )
            residue.atoms.append(atom)
        chain_b.residues.append(residue)

    structure.chains["A"] = chain_a
    structure.chains["B"] = chain_b

    return structure


@pytest.fixture
def interface_ab(minimal_structure: Structure) -> InterfaceRegion:
    """Interface between chains A and B."""
    from epitype.core.interface import detect_interface
    return detect_interface(minimal_structure, ["A"], ["B"], cutoff=10.0)


# Expected values for known structures
EXPECTED_METRICS_1YY9 = {
    "dSASA_int": (1800.0, 2200.0),  # Expected range
    "sc_value": (0.60, 0.75),
    "hbonds_int": (8, 20),
}

EXPECTED_METRICS_4FQI = {
    "dSASA_int": (1500.0, 2000.0),
    "sc_value": (0.55, 0.70),
    "hbonds_int": (5, 15),
}
