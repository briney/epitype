# Epitype: Technical Implementation Plan

## Overview

A Python toolkit for computing protein-protein interface metrics. This package provides both a CLI (`epitype`) and a Python API for analyzing binding interfaces.

### Core Metrics

| Metric | Description | Range/Units |
|--------|-------------|-------------|
| `dG_separated` | Binding energy (complex - separated) | REU (lower = better) |
| `dSASA_int` | Buried surface area at interface | Ų |
| `dG_separated/dSASA` | Binding energy density | REU/Ų (< -1.5 good) |
| `sc_value` | Shape complementarity | 0-1 (> 0.65 good) |
| `packstat` | Packing quality | 0-1 (> 0.65 good) |
| `delta_unsatHbonds` | Buried unsatisfied H-bond donors/acceptors | count (lower = better) |
| `hbonds_int` | Cross-interface hydrogen bonds | count |

### Technical Decisions

- **Python Version**: 3.10+
- **Energy Function**: OpenMM with AMBER ff14SB force field
- **Surface Generation**: MSMS (external tool)
- **Output Formats**: JSON and CSV
- **PackStat**: Multi-probe packing quality algorithm (30 probe radii)
- **Parallelization**: Single-structure focus (no built-in parallel processing)

---

## 1. Directory Structure

```
epitype/
├── pyproject.toml              # Package configuration, dependencies, build system
├── README.md                   # Project documentation
├── LICENSE                     # MIT License
├── CHANGELOG.md                # Version history
├── WORKPLAN.md                 # This file
│
├── src/
│   └── epitype/
│       ├── __init__.py         # Package exports, version
│       ├── py.typed            # PEP 561 marker for type hints
│       │
│       ├── cli/
│       │   ├── __init__.py
│       │   ├── main.py         # Main CLI entry point (typer app)
│       │   ├── run.py          # `epitype run` command
│       │   ├── sasa.py         # `epitype sasa` subcommand
│       │   ├── energy.py       # `epitype energy` subcommand
│       │   ├── hbonds.py       # `epitype hbonds` subcommand
│       │   ├── shape.py        # `epitype shape` subcommand
│       │   ├── packstat.py     # `epitype packstat` subcommand
│       │   └── common.py       # Shared CLI utilities, options
│       │
│       ├── core/
│       │   ├── __init__.py
│       │   ├── structure.py    # Structure, Chain, Residue, Atom classes
│       │   ├── interface.py    # InterfaceRegion, interface detection
│       │   ├── separation.py   # Rigid body separation logic
│       │   └── types.py        # Type definitions, enums, constants
│       │
│       ├── io/
│       │   ├── __init__.py
│       │   ├── parsers.py      # PDB/mmCIF parsing
│       │   ├── writers.py      # Output writers (JSON, CSV)
│       │   └── validation.py   # Input validation
│       │
│       ├── metrics/
│       │   ├── __init__.py
│       │   ├── sasa.py         # SASA calculation (FreeSASA wrapper)
│       │   ├── energy.py       # Binding energy (OpenMM)
│       │   ├── hbonds.py       # Hydrogen bond detection
│       │   ├── shape.py        # Shape complementarity
│       │   └── packstat.py     # PackStat packing quality
│       │
│       ├── surface/
│       │   ├── __init__.py
│       │   ├── msms.py         # MSMS wrapper
│       │   └── types.py        # Surface point representations
│       │
│       ├── analysis/
│       │   ├── __init__.py
│       │   └── pipeline.py     # Full analysis pipeline orchestration
│       │
│       └── utils/
│           ├── __init__.py
│           ├── geometry.py     # Vector math, transformations
│           ├── kdtree.py       # KD-tree utilities
│           └── units.py        # Unit conversions
│
├── tests/
│   ├── __init__.py
│   ├── conftest.py             # Pytest fixtures, test data paths
│   │
│   ├── unit/
│   │   ├── __init__.py
│   │   ├── test_structure.py
│   │   ├── test_interface.py
│   │   ├── test_separation.py
│   │   ├── test_sasa.py
│   │   ├── test_energy.py
│   │   ├── test_hbonds.py
│   │   ├── test_shape.py
│   │   ├── test_packstat.py
│   │   ├── test_geometry.py
│   │   └── test_parsers.py
│   │
│   ├── integration/
│   │   ├── __init__.py
│   │   ├── test_pipeline.py
│   │   ├── test_cli.py
│   │   └── test_api.py
│   │
│   ├── e2e/
│   │   ├── __init__.py
│   │   ├── test_full_analysis.py
│   │   └── test_known_structures.py
│   │
│   └── data/
│       ├── README.md           # Description of test structures
│       ├── 1YY9.pdb            # Antibody-antigen complex
│       ├── 4FQI.pdb            # Antibody-antigen complex
│       ├── minimal_dimer.pdb   # Small synthetic test case
│       └── expected/
│           ├── 1YY9_metrics.json
│           └── 4FQI_metrics.json
│
└── docs/
    ├── index.md
    ├── installation.md
    ├── cli.md
    ├── api.md
    └── metrics.md
```

---

## 2. Dependencies

### pyproject.toml

```toml
[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[project]
name = "interface-analyzer"
version = "0.1.0"
description = "Protein-protein interface analysis toolkit"
readme = "README.md"
license = "MIT"
requires-python = ">=3.10"
authors = [
    { name = "Bryan Briney", email = "bryan@example.com" }
]
classifiers = [
    "Development Status :: 3 - Alpha",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: MIT License",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3.11",
    "Programming Language :: Python :: 3.12",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]
keywords = ["protein", "interface", "binding", "structural biology"]

dependencies = [
    "numpy>=1.24.0",
    "scipy>=1.10.0",
    "biopython>=1.81",
    "freesasa>=2.2.0",
    "openmm>=8.0.0",
    "pdbfixer>=1.9",
    "typer>=0.9.0",
    "rich>=13.0.0",
    "pydantic>=2.0.0",
]

[project.optional-dependencies]
dev = [
    "pytest>=7.4.0",
    "pytest-cov>=4.1.0",
    "pytest-mock>=3.11.0",
    "mypy>=1.5.0",
    "ruff>=0.1.0",
    "pre-commit>=3.4.0",
]
docs = [
    "mkdocs>=1.5.0",
    "mkdocs-material>=9.4.0",
    "mkdocstrings[python]>=0.23.0",
]

[project.scripts]
epitype = "epitype.cli.main:app"

[project.urls]
Homepage = "https://github.com/briney/epitype"
Documentation = "https://interface-analyzer.readthedocs.io"
Repository = "https://github.com/briney/epitype"

[tool.hatch.build.targets.wheel]
packages = ["src/epitype"]

[tool.pytest.ini_options]
testpaths = ["tests"]
python_files = ["test_*.py"]
addopts = "-v --cov=epitype --cov-report=term-missing"

[tool.mypy]
python_version = "3.10"
strict = true
warn_return_any = true
warn_unused_ignores = true

[tool.ruff]
line-length = 100
target-version = "py310"

[tool.ruff.lint]
select = ["E", "F", "I", "N", "W", "UP", "B", "C4", "SIM"]
```

---

## 3. Core Infrastructure

### 3.1 Type Definitions (`src/epitype/core/types.py`)

```python
"""Core type definitions and constants."""
from dataclasses import dataclass, field
from enum import Enum
from typing import Literal

import numpy as np
from numpy.typing import NDArray

# Type aliases
Coords = NDArray[np.float64]  # Shape: (N, 3)
AtomIndex = int
ResidueIndex = int
ChainID = str

# Physical constants
PROBE_RADIUS = 1.4  # Water probe radius in Angstroms
DEFAULT_INTERFACE_CUTOFF = 8.0  # Cβ-Cβ distance cutoff
SEPARATION_DISTANCE = 500.0  # Distance to separate chains

class AtomType(Enum):
    """Atom classification for SASA calculations."""
    POLAR = "polar"
    NONPOLAR = "nonpolar"
    CHARGED = "charged"

@dataclass(frozen=True)
class VdWRadii:
    """Van der Waals radii for common atoms (Angstroms)."""
    C: float = 1.70
    N: float = 1.55
    O: float = 1.52
    S: float = 1.80
    H: float = 1.20
    P: float = 1.80

    def get(self, element: str, default: float = 1.70) -> float:
        """Get radius for element."""
        return getattr(self, element.upper(), default)

VDW_RADII = VdWRadii()

@dataclass
class InterfaceMetrics:
    """Complete interface analysis results."""
    dG_separated: float  # Binding energy in REU
    dSASA_int: float  # Buried surface area in Ų
    dG_per_dSASA: float  # Energy density
    sc_value: float  # Shape complementarity 0-1
    packstat: float  # Packing quality 0-1
    delta_unsat_hbonds: int  # Buried unsatisfied H-bonds
    hbonds_int: int  # Cross-interface H-bonds

    # Additional details
    interface_residues_1: list[str] = field(default_factory=list)
    interface_residues_2: list[str] = field(default_factory=list)
    sasa_complex: float = 0.0
    sasa_separated: float = 0.0

    def to_dict(self) -> dict:
        """Convert to dictionary for serialization."""
        return {
            "dG_separated": self.dG_separated,
            "dSASA_int": self.dSASA_int,
            "dG_per_dSASA": self.dG_per_dSASA,
            "sc_value": self.sc_value,
            "packstat": self.packstat,
            "delta_unsat_hbonds": self.delta_unsat_hbonds,
            "hbonds_int": self.hbonds_int,
            "interface_residues_1": self.interface_residues_1,
            "interface_residues_2": self.interface_residues_2,
            "sasa_complex": self.sasa_complex,
            "sasa_separated": self.sasa_separated,
        }
```

### 3.2 Structure Representation (`src/epitype/core/structure.py`)

```python
"""Molecular structure representation."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterator, Sequence

import numpy as np
from numpy.typing import NDArray

from .types import AtomType, Coords, VDW_RADII


@dataclass
class Atom:
    """Single atom representation."""
    index: int
    name: str  # e.g., "CA", "CB", "N"
    element: str  # e.g., "C", "N", "O"
    coords: NDArray[np.float64]  # Shape: (3,)
    residue_index: int
    chain_id: str
    b_factor: float = 0.0
    occupancy: float = 1.0

    @property
    def vdw_radius(self) -> float:
        """Van der Waals radius."""
        return VDW_RADII.get(self.element)

    @property
    def is_backbone(self) -> bool:
        """Check if backbone atom."""
        return self.name in {"N", "CA", "C", "O"}

    @property
    def is_hydrogen(self) -> bool:
        """Check if hydrogen atom."""
        return self.element == "H"

    def distance_to(self, other: Atom) -> float:
        """Euclidean distance to another atom."""
        return float(np.linalg.norm(self.coords - other.coords))


@dataclass
class Residue:
    """Single residue representation."""
    index: int
    name: str  # Three-letter code, e.g., "ALA"
    chain_id: str
    seq_num: int  # PDB residue number
    insertion_code: str = ""
    atoms: list[Atom] = field(default_factory=list)

    @property
    def ca(self) -> Atom | None:
        """Get alpha carbon."""
        for atom in self.atoms:
            if atom.name == "CA":
                return atom
        return None

    @property
    def cb(self) -> Atom | None:
        """Get beta carbon (CA for glycine)."""
        for atom in self.atoms:
            if atom.name == "CB":
                return atom
        # Fall back to CA for glycine
        return self.ca

    @property
    def center_of_mass(self) -> NDArray[np.float64]:
        """Calculate residue center of mass."""
        coords = np.array([a.coords for a in self.atoms])
        return coords.mean(axis=0)

    @property
    def full_id(self) -> str:
        """Full residue identifier: chain:name:seqnum."""
        return f"{self.chain_id}:{self.name}:{self.seq_num}{self.insertion_code}"

    def get_atom(self, name: str) -> Atom | None:
        """Get atom by name."""
        for atom in self.atoms:
            if atom.name == name:
                return atom
        return None


@dataclass
class Chain:
    """Single chain representation."""
    chain_id: str
    residues: list[Residue] = field(default_factory=list)

    @property
    def atoms(self) -> Iterator[Atom]:
        """Iterate over all atoms in chain."""
        for residue in self.residues:
            yield from residue.atoms

    @property
    def coords(self) -> NDArray[np.float64]:
        """All atom coordinates as array."""
        return np.array([a.coords for a in self.atoms])

    def __len__(self) -> int:
        return len(self.residues)


@dataclass
class Structure:
    """Complete molecular structure."""
    name: str
    chains: dict[str, Chain] = field(default_factory=dict)

    @property
    def atoms(self) -> Iterator[Atom]:
        """Iterate over all atoms."""
        for chain in self.chains.values():
            yield from chain.atoms

    @property
    def residues(self) -> Iterator[Residue]:
        """Iterate over all residues."""
        for chain in self.chains.values():
            yield from chain.residues

    @property
    def coords(self) -> NDArray[np.float64]:
        """All atom coordinates as array."""
        return np.array([a.coords for a in self.atoms])

    @property
    def num_atoms(self) -> int:
        """Total number of atoms."""
        return sum(1 for _ in self.atoms)

    @property
    def num_residues(self) -> int:
        """Total number of residues."""
        return sum(1 for _ in self.residues)

    def get_chain(self, chain_id: str) -> Chain | None:
        """Get chain by ID."""
        return self.chains.get(chain_id)

    def get_chains(self, chain_ids: Sequence[str]) -> list[Chain]:
        """Get multiple chains by ID."""
        return [self.chains[cid] for cid in chain_ids if cid in self.chains]

    def subset(self, chain_ids: Sequence[str]) -> Structure:
        """Create new structure with only specified chains."""
        new_structure = Structure(name=f"{self.name}_subset")
        for cid in chain_ids:
            if cid in self.chains:
                new_structure.chains[cid] = self.chains[cid]
        return new_structure

    def copy(self) -> Structure:
        """Deep copy of structure."""
        import copy
        return copy.deepcopy(self)

    def translate(self, vector: NDArray[np.float64]) -> None:
        """Translate all atoms by vector (in-place)."""
        for atom in self.atoms:
            atom.coords = atom.coords + vector
```

### 3.3 Interface Detection (`src/epitype/core/interface.py`)

```python
"""Interface detection and representation."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.spatial import KDTree

from .structure import Structure, Residue, Chain
from .types import DEFAULT_INTERFACE_CUTOFF


@dataclass
class InterfaceRegion:
    """Represents the interface between two chain groups."""

    # Chain groups (e.g., ["H", "L"] vs ["A"])
    group1_chains: list[str]
    group2_chains: list[str]

    # Interface residues from each group
    residues_group1: list[Residue] = field(default_factory=list)
    residues_group2: list[Residue] = field(default_factory=list)

    # Interface atoms
    atoms_group1_indices: list[int] = field(default_factory=list)
    atoms_group2_indices: list[int] = field(default_factory=list)

    # Contact pairs: (residue1_idx, residue2_idx)
    contact_pairs: list[tuple[int, int]] = field(default_factory=list)

    @property
    def num_interface_residues(self) -> int:
        """Total interface residues."""
        return len(self.residues_group1) + len(self.residues_group2)

    @property
    def residue_ids_group1(self) -> list[str]:
        """Full IDs of group1 interface residues."""
        return [r.full_id for r in self.residues_group1]

    @property
    def residue_ids_group2(self) -> list[str]:
        """Full IDs of group2 interface residues."""
        return [r.full_id for r in self.residues_group2]


def detect_interface(
    structure: Structure,
    group1_chains: list[str],
    group2_chains: list[str],
    cutoff: float = DEFAULT_INTERFACE_CUTOFF,
    use_cb: bool = True,
) -> InterfaceRegion:
    """
    Detect interface residues between two chain groups.

    Args:
        structure: Input structure
        group1_chains: Chain IDs for first group (e.g., ["H", "L"])
        group2_chains: Chain IDs for second group (e.g., ["A"])
        cutoff: Distance cutoff in Angstroms (Cβ-Cβ or CA-CA)
        use_cb: If True, use Cβ distances; otherwise CA

    Returns:
        InterfaceRegion with detected interface residues
    """
    interface = InterfaceRegion(
        group1_chains=group1_chains,
        group2_chains=group2_chains,
    )

    # Collect residues and representative coordinates for each group
    residues_g1: list[Residue] = []
    coords_g1: list[np.ndarray] = []

    residues_g2: list[Residue] = []
    coords_g2: list[np.ndarray] = []

    for chain_id in group1_chains:
        chain = structure.get_chain(chain_id)
        if chain is None:
            continue
        for residue in chain.residues:
            rep_atom = residue.cb if use_cb else residue.ca
            if rep_atom is not None:
                residues_g1.append(residue)
                coords_g1.append(rep_atom.coords)

    for chain_id in group2_chains:
        chain = structure.get_chain(chain_id)
        if chain is None:
            continue
        for residue in chain.residues:
            rep_atom = residue.cb if use_cb else residue.ca
            if rep_atom is not None:
                residues_g2.append(residue)
                coords_g2.append(rep_atom.coords)

    if not coords_g1 or not coords_g2:
        return interface

    coords_g1_arr = np.array(coords_g1)
    coords_g2_arr = np.array(coords_g2)

    # Build KD-tree for group2 and query with group1
    tree = KDTree(coords_g2_arr)

    interface_g1_indices: set[int] = set()
    interface_g2_indices: set[int] = set()
    contact_pairs: list[tuple[int, int]] = []

    # Query all group1 points against group2
    for i, coord in enumerate(coords_g1_arr):
        neighbors = tree.query_ball_point(coord, cutoff)
        if neighbors:
            interface_g1_indices.add(i)
            for j in neighbors:
                interface_g2_indices.add(j)
                contact_pairs.append((i, j))

    # Populate interface
    interface.residues_group1 = [residues_g1[i] for i in sorted(interface_g1_indices)]
    interface.residues_group2 = [residues_g2[i] for i in sorted(interface_g2_indices)]
    interface.contact_pairs = contact_pairs

    # Get atom indices for interface residues
    for residue in interface.residues_group1:
        interface.atoms_group1_indices.extend([a.index for a in residue.atoms])
    for residue in interface.residues_group2:
        interface.atoms_group2_indices.extend([a.index for a in residue.atoms])

    return interface


def parse_chain_groups(chain_spec: str) -> tuple[list[str], list[str]]:
    """
    Parse chain group specification.

    Format: "HL_A" means chains H+L vs chain A
    Format: "AB_CD" means chains A+B vs chains C+D

    Args:
        chain_spec: Chain specification string

    Returns:
        Tuple of (group1_chains, group2_chains)
    """
    if "_" not in chain_spec:
        raise ValueError(f"Invalid chain spec '{chain_spec}'. Expected format: 'HL_A'")

    parts = chain_spec.split("_")
    if len(parts) != 2:
        raise ValueError(f"Invalid chain spec '{chain_spec}'. Expected exactly one underscore.")

    group1 = list(parts[0])
    group2 = list(parts[1])

    return group1, group2
```

### 3.4 Rigid Body Separation (`src/epitype/core/separation.py`)

```python
"""Rigid body separation for binding energy calculations."""
from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from .structure import Structure
from .interface import InterfaceRegion
from .types import SEPARATION_DISTANCE


def compute_interface_normal(
    structure: Structure,
    interface: InterfaceRegion,
) -> NDArray[np.float64]:
    """
    Compute the interface normal vector.

    The normal points from group1 centroid to group2 centroid.

    Args:
        structure: Input structure
        interface: Detected interface region

    Returns:
        Unit normal vector (3,)
    """
    # Get centroids of interface residues
    coords_g1 = []
    for residue in interface.residues_group1:
        coords_g1.append(residue.center_of_mass)

    coords_g2 = []
    for residue in interface.residues_group2:
        coords_g2.append(residue.center_of_mass)

    if not coords_g1 or not coords_g2:
        # Default to z-axis if no interface
        return np.array([0.0, 0.0, 1.0])

    centroid_g1 = np.mean(coords_g1, axis=0)
    centroid_g2 = np.mean(coords_g2, axis=0)

    normal = centroid_g2 - centroid_g1
    norm = np.linalg.norm(normal)

    if norm < 1e-6:
        return np.array([0.0, 0.0, 1.0])

    return normal / norm


def separate_chains(
    structure: Structure,
    interface: InterfaceRegion,
    distance: float = SEPARATION_DISTANCE,
) -> Structure:
    """
    Create separated structure by translating group2 chains.

    Args:
        structure: Input structure (will be copied)
        interface: Interface region defining chain groups
        distance: Separation distance in Angstroms

    Returns:
        New structure with group2 chains translated
    """
    # Deep copy the structure
    separated = structure.copy()

    # Compute interface normal
    normal = compute_interface_normal(structure, interface)

    # Translation vector
    translation = normal * distance

    # Translate group2 chains
    for chain_id in interface.group2_chains:
        chain = separated.get_chain(chain_id)
        if chain is not None:
            for residue in chain.residues:
                for atom in residue.atoms:
                    atom.coords = atom.coords + translation

    return separated


def compute_separation_distance(
    structure: Structure,
    interface: InterfaceRegion,
    target_min_distance: float = 100.0,
) -> float:
    """
    Compute minimum separation distance to eliminate all contacts.

    Args:
        structure: Input structure
        interface: Interface region
        target_min_distance: Target minimum inter-group distance

    Returns:
        Required separation distance
    """
    from scipy.spatial.distance import cdist

    # Get all atom coordinates for each group
    coords_g1 = []
    coords_g2 = []

    for chain_id in interface.group1_chains:
        chain = structure.get_chain(chain_id)
        if chain:
            for atom in chain.atoms:
                coords_g1.append(atom.coords)

    for chain_id in interface.group2_chains:
        chain = structure.get_chain(chain_id)
        if chain:
            for atom in chain.atoms:
                coords_g2.append(atom.coords)

    if not coords_g1 or not coords_g2:
        return SEPARATION_DISTANCE

    # Compute current minimum distance
    distances = cdist(coords_g1, coords_g2)
    current_min = distances.min()

    # Required additional separation
    return max(target_min_distance - current_min, 0) + SEPARATION_DISTANCE
```

---

## 4. I/O Module

### 4.1 Parsers (`src/epitype/io/parsers.py`)

```python
"""Structure file parsers."""
from __future__ import annotations

from pathlib import Path
from typing import TextIO

import numpy as np
from Bio.PDB import PDBParser as BioPDBParser
from Bio.PDB import MMCIFParser
from Bio.PDB.Structure import Structure as BioStructure

from epitype.core.structure import Structure, Chain, Residue, Atom


def parse_pdb(filepath: str | Path) -> Structure:
    """
    Parse PDB file into Structure.

    Args:
        filepath: Path to PDB file

    Returns:
        Parsed Structure object
    """
    filepath = Path(filepath)
    parser = BioPDBParser(QUIET=True)
    bio_structure = parser.get_structure(filepath.stem, str(filepath))

    return _convert_biopython_structure(bio_structure, filepath.stem)


def parse_mmcif(filepath: str | Path) -> Structure:
    """
    Parse mmCIF file into Structure.

    Args:
        filepath: Path to mmCIF file

    Returns:
        Parsed Structure object
    """
    filepath = Path(filepath)
    parser = MMCIFParser(QUIET=True)
    bio_structure = parser.get_structure(filepath.stem, str(filepath))

    return _convert_biopython_structure(bio_structure, filepath.stem)


def parse_structure(filepath: str | Path) -> Structure:
    """
    Parse structure file (auto-detect format).

    Args:
        filepath: Path to PDB or mmCIF file

    Returns:
        Parsed Structure object
    """
    filepath = Path(filepath)
    suffix = filepath.suffix.lower()

    if suffix in {".pdb", ".ent"}:
        return parse_pdb(filepath)
    elif suffix in {".cif", ".mmcif"}:
        return parse_mmcif(filepath)
    else:
        # Try PDB first, then mmCIF
        try:
            return parse_pdb(filepath)
        except Exception:
            return parse_mmcif(filepath)


def _convert_biopython_structure(bio_struct: BioStructure, name: str) -> Structure:
    """Convert BioPython structure to our Structure class."""
    structure = Structure(name=name)

    atom_index = 0
    residue_index = 0

    # Take first model only
    model = bio_struct[0]

    for bio_chain in model:
        chain_id = bio_chain.id
        chain = Chain(chain_id=chain_id)

        for bio_residue in bio_chain:
            # Skip water and heteroatoms (unless modified amino acid)
            if bio_residue.id[0] != " ":
                # Keep modified residues (e.g., MSE)
                if bio_residue.id[0] not in {"W", "H_HOH"}:
                    pass  # Include modified residues
                else:
                    continue

            residue = Residue(
                index=residue_index,
                name=bio_residue.resname,
                chain_id=chain_id,
                seq_num=bio_residue.id[1],
                insertion_code=bio_residue.id[2].strip(),
            )

            for bio_atom in bio_residue:
                element = bio_atom.element.strip() if bio_atom.element else bio_atom.name[0]

                atom = Atom(
                    index=atom_index,
                    name=bio_atom.name,
                    element=element,
                    coords=np.array(bio_atom.coord, dtype=np.float64),
                    residue_index=residue_index,
                    chain_id=chain_id,
                    b_factor=bio_atom.bfactor,
                    occupancy=bio_atom.occupancy,
                )
                residue.atoms.append(atom)
                atom_index += 1

            chain.residues.append(residue)
            residue_index += 1

        if chain.residues:
            structure.chains[chain_id] = chain

    return structure
```

### 4.2 Writers (`src/epitype/io/writers.py`)

```python
"""Output writers for analysis results."""
from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from epitype.core.types import InterfaceMetrics


def write_json(
    metrics: InterfaceMetrics,
    filepath: str | Path,
    structure_name: str | None = None,
    indent: int = 2,
) -> None:
    """
    Write metrics to JSON file.

    Args:
        metrics: Analysis results
        filepath: Output file path
        structure_name: Optional structure identifier
        indent: JSON indentation
    """
    filepath = Path(filepath)

    data: dict[str, Any] = {
        "structure": structure_name,
        "metrics": metrics.to_dict(),
    }

    with open(filepath, "w") as f:
        json.dump(data, f, indent=indent)


def write_csv(
    metrics: InterfaceMetrics,
    filepath: str | Path,
    structure_name: str | None = None,
    append: bool = False,
) -> None:
    """
    Write metrics to CSV file.

    Args:
        metrics: Analysis results
        filepath: Output file path
        structure_name: Optional structure identifier
        append: If True, append to existing file
    """
    filepath = Path(filepath)

    fieldnames = [
        "structure",
        "dG_separated",
        "dSASA_int",
        "dG_per_dSASA",
        "sc_value",
        "packstat",
        "delta_unsat_hbonds",
        "hbonds_int",
        "sasa_complex",
        "sasa_separated",
    ]

    row = {
        "structure": structure_name or "",
        **{k: v for k, v in metrics.to_dict().items() if k in fieldnames},
    }

    mode = "a" if append and filepath.exists() else "w"
    write_header = not (append and filepath.exists())

    with open(filepath, mode, newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if write_header:
            writer.writeheader()
        writer.writerow(row)


def format_metrics_table(
    metrics: InterfaceMetrics,
    structure_name: str | None = None,
) -> str:
    """
    Format metrics as human-readable table.

    Args:
        metrics: Analysis results
        structure_name: Optional structure identifier

    Returns:
        Formatted string
    """
    lines = []

    if structure_name:
        lines.append(f"Structure: {structure_name}")
        lines.append("=" * 50)

    lines.append(f"{'Metric':<25} {'Value':>15} {'Quality':>10}")
    lines.append("-" * 50)

    # dG_separated
    quality = "good" if metrics.dG_separated < 0 else "poor"
    lines.append(f"{'dG_separated (REU)':<25} {metrics.dG_separated:>15.2f} {quality:>10}")

    # dSASA
    lines.append(f"{'dSASA_int (Ų)':<25} {metrics.dSASA_int:>15.1f}")

    # dG/dSASA
    quality = "good" if metrics.dG_per_dSASA < -1.5 else "moderate" if metrics.dG_per_dSASA < 0 else "poor"
    lines.append(f"{'dG/dSASA (REU/Ų)':<25} {metrics.dG_per_dSASA:>15.3f} {quality:>10}")

    # Shape complementarity
    quality = "good" if metrics.sc_value > 0.65 else "moderate" if metrics.sc_value > 0.5 else "poor"
    lines.append(f"{'Shape complementarity':<25} {metrics.sc_value:>15.3f} {quality:>10}")

    # Packstat
    quality = "good" if metrics.packstat > 0.65 else "moderate" if metrics.packstat > 0.5 else "poor"
    lines.append(f"{'Packstat':<25} {metrics.packstat:>15.3f} {quality:>10}")

    # H-bonds
    lines.append(f"{'Unsat. H-bonds':<25} {metrics.delta_unsat_hbonds:>15d}")
    lines.append(f"{'Interface H-bonds':<25} {metrics.hbonds_int:>15d}")

    return "\n".join(lines)
```

---

## 5. Metrics Implementation

### 5.1 SASA Calculation (`src/epitype/metrics/sasa.py`)

```python
"""Solvent Accessible Surface Area calculations using FreeSASA."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import freesasa
import numpy as np

from epitype.core.structure import Structure, Atom
from epitype.core.interface import InterfaceRegion
from epitype.core.types import PROBE_RADIUS


@dataclass
class SASAResult:
    """SASA calculation results."""
    total: float  # Total SASA in Ų
    polar: float  # Polar atom SASA
    apolar: float  # Apolar atom SASA
    per_atom: dict[int, float]  # SASA per atom index
    per_residue: dict[str, float]  # SASA per residue ID


def calculate_sasa(
    structure: Structure,
    probe_radius: float = PROBE_RADIUS,
    algorithm: str = "ShrakeRupley",
) -> SASAResult:
    """
    Calculate SASA for a structure using FreeSASA.

    Args:
        structure: Input structure
        probe_radius: Probe radius in Angstroms
        algorithm: "ShrakeRupley" or "LeeRichards"

    Returns:
        SASAResult with total and per-atom SASA
    """
    # Prepare data for FreeSASA
    coords: list[float] = []
    radii: list[float] = []
    atom_list: list[Atom] = []

    for atom in structure.atoms:
        coords.extend(atom.coords.tolist())
        radii.append(atom.vdw_radius)
        atom_list.append(atom)

    # Configure FreeSASA
    params = freesasa.Parameters()
    params.setProbeRadius(probe_radius)

    if algorithm == "LeeRichards":
        params.setAlgorithm(freesasa.LeeRichards)
    else:
        params.setAlgorithm(freesasa.ShrakeRupley)

    # Calculate
    result = freesasa.calcCoord(coords, radii, params)

    # Extract per-atom SASA
    per_atom: dict[int, float] = {}
    per_residue: dict[str, float] = {}
    polar_sasa = 0.0
    apolar_sasa = 0.0

    polar_elements = {"N", "O", "S"}

    for i, atom in enumerate(atom_list):
        sasa = result.atomArea(i)
        per_atom[atom.index] = sasa

        # Accumulate per-residue
        res_id = f"{atom.chain_id}:{atom.residue_index}"
        per_residue[res_id] = per_residue.get(res_id, 0.0) + sasa

        # Polar vs apolar
        if atom.element in polar_elements:
            polar_sasa += sasa
        else:
            apolar_sasa += sasa

    return SASAResult(
        total=result.totalArea(),
        polar=polar_sasa,
        apolar=apolar_sasa,
        per_atom=per_atom,
        per_residue=per_residue,
    )


def calculate_dsasa(
    complex_structure: Structure,
    separated_structure: Structure,
    interface: InterfaceRegion,
    probe_radius: float = PROBE_RADIUS,
) -> float:
    """
    Calculate change in SASA upon complex formation (dSASA).

    dSASA = SASA_separated - SASA_complex

    Args:
        complex_structure: Bound complex
        separated_structure: Separated chains
        interface: Interface region
        probe_radius: Probe radius in Angstroms

    Returns:
        dSASA (buried surface area) in Ų
    """
    sasa_complex = calculate_sasa(complex_structure, probe_radius)
    sasa_separated = calculate_sasa(separated_structure, probe_radius)

    return sasa_separated.total - sasa_complex.total


def calculate_interface_sasa(
    structure: Structure,
    interface: InterfaceRegion,
    probe_radius: float = PROBE_RADIUS,
) -> tuple[float, float]:
    """
    Calculate SASA for interface residues only.

    Args:
        structure: Input structure
        interface: Interface region
        probe_radius: Probe radius

    Returns:
        Tuple of (group1_interface_sasa, group2_interface_sasa)
    """
    sasa_result = calculate_sasa(structure, probe_radius)

    g1_sasa = 0.0
    for residue in interface.residues_group1:
        res_id = f"{residue.chain_id}:{residue.index}"
        g1_sasa += sasa_result.per_residue.get(res_id, 0.0)

    g2_sasa = 0.0
    for residue in interface.residues_group2:
        res_id = f"{residue.chain_id}:{residue.index}"
        g2_sasa += sasa_result.per_residue.get(res_id, 0.0)

    return g1_sasa, g2_sasa
```

### 5.2 Binding Energy (`src/epitype/metrics/energy.py`)

```python
"""Binding energy calculations using OpenMM with AMBER ff14SB."""
from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
from openmm import app, unit, Platform, LangevinMiddleIntegrator
from openmm.app import PDBFile, ForceField, Modeller
from pdbfixer import PDBFixer

from epitype.core.structure import Structure
from epitype.core.interface import InterfaceRegion


# Conversion: kJ/mol to kcal/mol
KJ_TO_REU = 0.239  # 1 kJ/mol ≈ 0.239 kcal/mol


def prepare_structure_for_openmm(structure: Structure) -> tuple[Modeller, ForceField]:
    """
    Prepare structure for OpenMM energy calculation.

    Includes fixing missing atoms/residues via PDBFixer.

    Args:
        structure: Input structure

    Returns:
        Tuple of (Modeller, ForceField)
    """
    # Write to temporary PDB
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as f:
        temp_path = Path(f.name)
        _write_temp_pdb(structure, temp_path)

    # Fix structure
    fixer = PDBFixer(filename=str(temp_path))
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.0)

    # Clean up temp file
    temp_path.unlink()

    # Create force field
    forcefield = ForceField("amber14-all.xml", "amber14/tip3p.xml")

    # Create modeller (no solvent for interface calculations)
    modeller = Modeller(fixer.topology, fixer.positions)

    return modeller, forcefield


def calculate_potential_energy(
    modeller: Modeller,
    forcefield: ForceField,
    minimize: bool = False,
) -> float:
    """
    Calculate potential energy of a structure.

    Args:
        modeller: OpenMM Modeller with structure
        forcefield: OpenMM ForceField
        minimize: If True, minimize before calculating energy

    Returns:
        Potential energy in kJ/mol
    """
    # Create system
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.NoCutoff,  # No cutoff for accuracy
        constraints=None,
        rigidWater=False,
    )

    # Create integrator (not used for single-point energy)
    integrator = LangevinMiddleIntegrator(
        300 * unit.kelvin,
        1.0 / unit.picosecond,
        0.002 * unit.picoseconds,
    )

    # Create simulation
    platform = Platform.getPlatformByName("CPU")
    simulation = app.Simulation(modeller.topology, system, integrator, platform)
    simulation.context.setPositions(modeller.positions)

    if minimize:
        simulation.minimizeEnergy(maxIterations=100)

    # Get energy
    state = simulation.context.getState(getEnergy=True)
    energy = state.getPotentialEnergy()

    return energy.value_in_unit(unit.kilojoule_per_mole)


def calculate_binding_energy(
    complex_structure: Structure,
    separated_structure: Structure,
    interface: InterfaceRegion,
    minimize: bool = True,
) -> float:
    """
    Calculate binding energy (dG_separated).

    dG = E_complex - E_separated

    Args:
        complex_structure: Bound complex
        separated_structure: Separated chains
        interface: Interface region (for reference)
        minimize: If True, minimize structures before energy calculation

    Returns:
        Binding energy in REU (negative = favorable)
    """
    # Prepare complex
    complex_modeller, forcefield = prepare_structure_for_openmm(complex_structure)
    energy_complex = calculate_potential_energy(complex_modeller, forcefield, minimize)

    # Prepare separated
    sep_modeller, forcefield = prepare_structure_for_openmm(separated_structure)
    energy_separated = calculate_potential_energy(sep_modeller, forcefield, minimize)

    # Calculate delta (complex - separated)
    dG_kj = energy_complex - energy_separated

    # Convert to REU
    return dG_kj * KJ_TO_REU


def _write_temp_pdb(structure: Structure, path: Path) -> None:
    """Write structure to temporary PDB file."""
    with open(path, "w") as f:
        atom_num = 1
        for chain in structure.chains.values():
            for residue in chain.residues:
                for atom in residue.atoms:
                    f.write(
                        f"ATOM  {atom_num:5d} {atom.name:<4s} "
                        f"{residue.name:3s} {chain.chain_id:1s}"
                        f"{residue.seq_num:4d}{residue.insertion_code:1s}   "
                        f"{atom.coords[0]:8.3f}{atom.coords[1]:8.3f}{atom.coords[2]:8.3f}"
                        f"{atom.occupancy:6.2f}{atom.b_factor:6.2f}          "
                        f"{atom.element:>2s}\n"
                    )
                    atom_num += 1
        f.write("END\n")
```

### 5.3 Hydrogen Bond Detection (`src/epitype/metrics/hbonds.py`)

```python
"""Hydrogen bond detection and analysis."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from scipy.spatial import KDTree

from epitype.core.structure import Structure, Atom, Residue
from epitype.core.interface import InterfaceRegion


# Hydrogen bond geometry criteria
HBOND_DISTANCE_MAX = 3.5  # D-A distance in Angstroms
HBOND_ANGLE_MIN = 120.0  # D-H-A angle minimum in degrees


@dataclass
class HydrogenBond:
    """Single hydrogen bond."""
    donor_atom: Atom
    acceptor_atom: Atom
    hydrogen_atom: Atom | None
    distance: float  # D-A distance
    angle: float | None  # D-H-A angle if H available


@dataclass
class HBondResult:
    """Hydrogen bond analysis results."""
    total_hbonds: int
    interface_hbonds: int  # Cross-interface H-bonds
    hbonds: list[HydrogenBond]
    buried_unsatisfied_donors: int
    buried_unsatisfied_acceptors: int


# Donor atoms: N-H groups
DONOR_ATOMS = {
    # Backbone
    "N",
    # Sidechains
    "NE", "NH1", "NH2",  # Arg
    "ND2",  # Asn
    "NE2",  # Gln
    "NZ",  # Lys
    "NE1",  # Trp
    "ND1", "NE2",  # His (can be either)
    "OG", "OG1",  # Ser, Thr (hydroxyl)
    "OH",  # Tyr
}

# Acceptor atoms: O, N with lone pairs
ACCEPTOR_ATOMS = {
    # Backbone
    "O",
    # Sidechains
    "OD1", "OD2",  # Asp
    "OE1", "OE2",  # Glu
    "OD1",  # Asn
    "OE1",  # Gln
    "ND1", "NE2",  # His
    "OG", "OG1",  # Ser, Thr
    "OH",  # Tyr
    "SD",  # Met (weak)
}


def detect_hydrogen_bonds(
    structure: Structure,
    distance_cutoff: float = HBOND_DISTANCE_MAX,
    angle_cutoff: float = HBOND_ANGLE_MIN,
) -> list[HydrogenBond]:
    """
    Detect all hydrogen bonds in structure.

    Uses geometric criteria:
    - D-A distance < 3.5 Å
    - D-H-A angle > 120° (if H available)

    Args:
        structure: Input structure
        distance_cutoff: Maximum D-A distance
        angle_cutoff: Minimum D-H-A angle

    Returns:
        List of detected hydrogen bonds
    """
    # Collect donors and acceptors
    donors: list[Atom] = []
    acceptors: list[Atom] = []

    for atom in structure.atoms:
        if atom.name in DONOR_ATOMS:
            donors.append(atom)
        if atom.name in ACCEPTOR_ATOMS:
            acceptors.append(atom)

    if not donors or not acceptors:
        return []

    # Build KD-tree for acceptors
    acceptor_coords = np.array([a.coords for a in acceptors])
    tree = KDTree(acceptor_coords)

    hbonds: list[HydrogenBond] = []

    for donor in donors:
        # Find nearby acceptors
        indices = tree.query_ball_point(donor.coords, distance_cutoff)

        for idx in indices:
            acceptor = acceptors[idx]

            # Skip same residue
            if donor.residue_index == acceptor.residue_index:
                continue

            distance = donor.distance_to(acceptor)

            # Create H-bond (angle check would require finding H atom)
            hbond = HydrogenBond(
                donor_atom=donor,
                acceptor_atom=acceptor,
                hydrogen_atom=None,  # Could search for attached H
                distance=distance,
                angle=None,
            )
            hbonds.append(hbond)

    return hbonds


def count_interface_hbonds(
    hbonds: list[HydrogenBond],
    interface: InterfaceRegion,
) -> int:
    """
    Count hydrogen bonds crossing the interface.

    Args:
        hbonds: All detected H-bonds
        interface: Interface region

    Returns:
        Number of cross-interface H-bonds
    """
    g1_chains = set(interface.group1_chains)
    g2_chains = set(interface.group2_chains)

    count = 0
    for hbond in hbonds:
        donor_in_g1 = hbond.donor_atom.chain_id in g1_chains
        acceptor_in_g1 = hbond.acceptor_atom.chain_id in g1_chains
        donor_in_g2 = hbond.donor_atom.chain_id in g2_chains
        acceptor_in_g2 = hbond.acceptor_atom.chain_id in g2_chains

        # Cross-interface if one atom in g1 and other in g2
        if (donor_in_g1 and acceptor_in_g2) or (donor_in_g2 and acceptor_in_g1):
            count += 1

    return count


def count_unsatisfied_hbonds(
    structure: Structure,
    interface: InterfaceRegion,
    hbonds: list[HydrogenBond],
    burial_threshold: float = 0.1,  # SASA threshold for "buried"
) -> tuple[int, int]:
    """
    Count buried unsatisfied H-bond donors and acceptors.

    An atom is unsatisfied if:
    - It's buried (low SASA)
    - It's a potential donor/acceptor
    - It's not participating in an H-bond

    Args:
        structure: Input structure
        interface: Interface region
        hbonds: Detected H-bonds
        burial_threshold: SASA threshold for burial

    Returns:
        Tuple of (unsatisfied_donors, unsatisfied_acceptors)
    """
    from epitype.metrics.sasa import calculate_sasa

    # Get per-atom SASA
    sasa_result = calculate_sasa(structure)

    # Find atoms in H-bonds
    satisfied_donors = {hb.donor_atom.index for hb in hbonds}
    satisfied_acceptors = {hb.acceptor_atom.index for hb in hbonds}

    # Get interface atom indices
    interface_atoms = set(interface.atoms_group1_indices) | set(interface.atoms_group2_indices)

    unsatisfied_donors = 0
    unsatisfied_acceptors = 0

    for atom in structure.atoms:
        # Only check interface atoms
        if atom.index not in interface_atoms:
            continue

        # Check if buried
        atom_sasa = sasa_result.per_atom.get(atom.index, 0.0)
        if atom_sasa > burial_threshold:
            continue

        # Check donors
        if atom.name in DONOR_ATOMS and atom.index not in satisfied_donors:
            unsatisfied_donors += 1

        # Check acceptors
        if atom.name in ACCEPTOR_ATOMS and atom.index not in satisfied_acceptors:
            unsatisfied_acceptors += 1

    return unsatisfied_donors, unsatisfied_acceptors


def analyze_hbonds(
    structure: Structure,
    interface: InterfaceRegion,
) -> HBondResult:
    """
    Complete hydrogen bond analysis.

    Args:
        structure: Input structure
        interface: Interface region

    Returns:
        HBondResult with all H-bond metrics
    """
    hbonds = detect_hydrogen_bonds(structure)
    interface_count = count_interface_hbonds(hbonds, interface)
    unsat_donors, unsat_acceptors = count_unsatisfied_hbonds(
        structure, interface, hbonds
    )

    return HBondResult(
        total_hbonds=len(hbonds),
        interface_hbonds=interface_count,
        hbonds=hbonds,
        buried_unsatisfied_donors=unsat_donors,
        buried_unsatisfied_acceptors=unsat_acceptors,
    )
```

### 5.4 Shape Complementarity (`src/epitype/metrics/shape.py`)

```python
"""Shape complementarity calculation using MSMS surfaces."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from scipy.spatial import KDTree
from numpy.typing import NDArray

from epitype.core.structure import Structure
from epitype.core.interface import InterfaceRegion
from epitype.surface.msms import generate_surface, SurfacePoint


@dataclass
class ShapeComplementarityResult:
    """Shape complementarity analysis results."""
    sc_value: float  # Overall Sc (0-1)
    sc_median: float  # Median Sc
    sc_group1: float  # Sc for group1 surface
    sc_group2: float  # Sc for group2 surface
    interface_area: float  # Interface surface area


def calculate_shape_complementarity(
    structure: Structure,
    interface: InterfaceRegion,
    distance_cutoff: float = 1.5,
    trim_distance: float = 10.0,
) -> ShapeComplementarityResult:
    """
    Calculate shape complementarity (Sc) between interface surfaces.

    Algorithm (Lawrence & Colman, 1993):
    1. Generate molecular surfaces for each chain group
    2. For each surface point on one side, find nearest point on other side
    3. Calculate Sc = dot(n1, n2) * exp(-distance² / 0.5)
    4. Average over all interface points

    Args:
        structure: Input structure
        interface: Interface region
        distance_cutoff: Max distance between surface points to consider
        trim_distance: Distance from interface to include surface points

    Returns:
        ShapeComplementarityResult
    """
    # Generate surfaces for each group
    group1_structure = structure.subset(interface.group1_chains)
    group2_structure = structure.subset(interface.group2_chains)

    surface1 = generate_surface(group1_structure)
    surface2 = generate_surface(group2_structure)

    if not surface1 or not surface2:
        return ShapeComplementarityResult(
            sc_value=0.0, sc_median=0.0, sc_group1=0.0, sc_group2=0.0, interface_area=0.0
        )

    # Trim to interface region
    interface_surface1 = _trim_surface_to_interface(surface1, surface2, trim_distance)
    interface_surface2 = _trim_surface_to_interface(surface2, surface1, trim_distance)

    if not interface_surface1 or not interface_surface2:
        return ShapeComplementarityResult(
            sc_value=0.0, sc_median=0.0, sc_group1=0.0, sc_group2=0.0, interface_area=0.0
        )

    # Calculate Sc for each direction
    sc_values_1to2 = _calculate_sc_one_direction(
        interface_surface1, interface_surface2, distance_cutoff
    )
    sc_values_2to1 = _calculate_sc_one_direction(
        interface_surface2, interface_surface1, distance_cutoff
    )

    # Combine
    sc_group1 = np.mean(sc_values_1to2) if sc_values_1to2 else 0.0
    sc_group2 = np.mean(sc_values_2to1) if sc_values_2to1 else 0.0

    all_sc = np.concatenate([sc_values_1to2, sc_values_2to1]) if (len(sc_values_1to2) and len(sc_values_2to1)) else np.array([0.0])

    return ShapeComplementarityResult(
        sc_value=float(np.mean(all_sc)),
        sc_median=float(np.median(all_sc)),
        sc_group1=float(sc_group1),
        sc_group2=float(sc_group2),
        interface_area=len(interface_surface1) + len(interface_surface2),
    )


def _trim_surface_to_interface(
    surface: list[SurfacePoint],
    other_surface: list[SurfacePoint],
    trim_distance: float,
) -> list[SurfacePoint]:
    """Keep only surface points near the other surface."""
    if not other_surface:
        return []

    other_coords = np.array([p.coords for p in other_surface])
    tree = KDTree(other_coords)

    trimmed = []
    for point in surface:
        dist, _ = tree.query(point.coords)
        if dist < trim_distance:
            trimmed.append(point)

    return trimmed


def _calculate_sc_one_direction(
    from_surface: list[SurfacePoint],
    to_surface: list[SurfacePoint],
    distance_cutoff: float,
) -> NDArray[np.float64]:
    """
    Calculate Sc values from one surface to another.

    For each point in from_surface, find nearest in to_surface
    and compute Sc contribution.
    """
    if not from_surface or not to_surface:
        return np.array([])

    to_coords = np.array([p.coords for p in to_surface])
    to_normals = np.array([p.normal for p in to_surface])
    tree = KDTree(to_coords)

    sc_values = []

    for point in from_surface:
        dist, idx = tree.query(point.coords)

        if dist > distance_cutoff:
            continue

        # Get normals
        n1 = point.normal
        n2 = to_normals[idx]

        # Sc = -dot(n1, n2) * exp(-d² / σ²)
        # Negative because normals should point away from each other
        # σ² = 0.5 is typical
        dot_product = -np.dot(n1, n2)
        distance_weight = np.exp(-dist * dist / 0.5)

        sc = dot_product * distance_weight
        sc_values.append(sc)

    return np.array(sc_values)
```

### 5.5 PackStat (`src/epitype/metrics/packstat.py`)

> **Note**: This implementation uses FreeSASA for fast SASA calculations at multiple
> probe radii, providing ~100x speedup over custom Python implementations.

```python
"""PackStat packing quality metric using multi-probe algorithm.

This implementation uses FreeSASA for fast SASA calculations at multiple
probe radii, providing ~100x speedup over custom Python implementations.
"""
from __future__ import annotations

from dataclasses import dataclass

import freesasa
import numpy as np
from numpy.typing import NDArray

from epitype.core.structure import Structure
from epitype.core.interface import InterfaceRegion


# Probe radii for multi-probe SASA (30 values from 0.1 to 3.0 Å)
PROBE_RADII = np.linspace(0.1, 3.0, 30)


@dataclass
class PackStatResult:
    """PackStat analysis results."""
    packstat: float  # Overall packing score (0-1)
    packstat_interface: float  # Interface-specific packing
    cavity_volume: float  # Estimated cavity volume (ų)
    per_probe_scores: dict[float, float]  # Score at each probe radius


def calculate_packstat(
    structure: Structure,
    interface: InterfaceRegion | None = None,
) -> PackStatResult:
    """
    Calculate PackStat using multi-probe algorithm.

    Uses FreeSASA for fast SASA calculations at multiple probe radii.

    Algorithm:
    1. For each probe radius (0.1 to 3.0 Å, 30 values):
       a. Compute SASA with that probe using FreeSASA
       b. Calculate packing fraction (buried / reference surface)
    2. Weight and combine across probe radii (favoring smaller probes)

    Args:
        structure: Input structure
        interface: Optional interface region (for interface-specific score)

    Returns:
        PackStatResult with packing metrics
    """
    # Get all atom data
    atoms = list(structure.atoms)

    if len(atoms) == 0:
        return PackStatResult(
            packstat=0.0, packstat_interface=0.0, cavity_volume=0.0, per_probe_scores={}
        )

    # Prepare data for FreeSASA
    coords: list[float] = []
    radii: list[float] = []

    for atom in atoms:
        coords.extend(atom.coords.tolist())
        radii.append(atom.vdw_radius)

    radii_arr = np.array(radii)

    # Calculate SASA at each probe radius
    per_probe_scores: dict[float, float] = {}
    contact_areas: list[float] = []

    for probe_r in PROBE_RADII:
        # SASA at this probe radius using FreeSASA
        total_sasa = _freesasa_total(coords, radii, probe_r)

        # Reference SASA (isolated atoms)
        ref_sasa = _reference_sasa(radii_arr, probe_r)

        # Packing fraction: how much surface is buried
        if ref_sasa > 0:
            packing = 1.0 - (total_sasa / ref_sasa)
        else:
            packing = 0.0

        per_probe_scores[float(probe_r)] = max(0.0, min(1.0, packing))
        contact_areas.append(packing)

    # Overall packstat: weighted average favoring smaller probes
    # (smaller probes detect finer packing details)
    weights = np.exp(-PROBE_RADII)  # Exponential decay with probe size
    weights /= weights.sum()

    packstat = float(np.average(contact_areas, weights=weights))

    # Interface-specific packing
    packstat_interface = packstat
    if interface is not None:
        packstat_interface = _interface_packstat(structure, interface)

    # Estimate cavity volume from intermediate probe sizes
    cavity_volume = _estimate_cavity_volume(radii_arr, contact_areas)

    return PackStatResult(
        packstat=packstat,
        packstat_interface=packstat_interface,
        cavity_volume=cavity_volume,
        per_probe_scores=per_probe_scores,
    )


def _freesasa_total(
    coords: list[float],
    radii: list[float],
    probe_radius: float,
) -> float:
    """
    Calculate total SASA using FreeSASA.

    Args:
        coords: Flattened coordinates [x1, y1, z1, x2, y2, z2, ...]
        radii: Van der Waals radii for each atom
        probe_radius: Probe radius in Angstroms

    Returns:
        Total SASA in Ų
    """
    params = freesasa.Parameters()
    params.setProbeRadius(probe_radius)
    params.setAlgorithm(freesasa.ShrakeRupley)
    result = freesasa.calcCoord(coords, radii, params)
    return result.totalArea()


def _reference_sasa(radii: NDArray[np.float64], probe_radius: float) -> float:
    """Calculate reference SASA for isolated atoms."""
    eff_radii = radii + probe_radius
    return float(np.sum(4 * np.pi * eff_radii ** 2))


def _interface_packstat(
    structure: Structure,
    interface: InterfaceRegion,
) -> float:
    """Calculate packstat for interface atoms only."""
    # Get interface atom indices
    interface_indices = set(interface.atoms_group1_indices) | set(
        interface.atoms_group2_indices
    )

    atoms = [a for a in structure.atoms if a.index in interface_indices]
    if not atoms:
        return 0.0

    # Prepare data for FreeSASA
    coords: list[float] = []
    radii: list[float] = []

    for atom in atoms:
        coords.extend(atom.coords.tolist())
        radii.append(atom.vdw_radius)

    radii_arr = np.array(radii)

    scores = []
    weights = []

    for probe_r in PROBE_RADII:
        total_sasa = _freesasa_total(coords, radii, probe_r)
        ref_sasa = _reference_sasa(radii_arr, probe_r)

        if ref_sasa > 0:
            packing = 1.0 - (total_sasa / ref_sasa)
        else:
            packing = 0.0

        scores.append(max(0.0, min(1.0, packing)))
        weights.append(np.exp(-probe_r))

    weights_arr = np.array(weights)
    weights_arr /= weights_arr.sum()

    return float(np.average(scores, weights=weights_arr))


def _estimate_cavity_volume(
    radii: NDArray[np.float64],
    packing_scores: list[float],
) -> float:
    """Estimate cavity volume from packing scores."""
    # Rough estimate: volume deficit based on packing scores
    total_vdw_volume = float(np.sum((4 / 3) * np.pi * radii ** 3))

    # Average packing deficit at intermediate probe sizes
    mid_scores = packing_scores[10:20]  # Probe radii ~1.0-2.0 Å
    if mid_scores:
        avg_deficit = 1.0 - np.mean(mid_scores)
    else:
        avg_deficit = 0.0

    return total_vdw_volume * avg_deficit * 0.1  # Scaling factor
```

### 5.6 MSMS Surface Wrapper (`src/epitype/surface/msms.py`)

```python
"""MSMS molecular surface generation wrapper."""
from __future__ import annotations

import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from epitype.core.structure import Structure


@dataclass
class SurfacePoint:
    """Single point on molecular surface."""
    coords: NDArray[np.float64]  # (3,) position
    normal: NDArray[np.float64]  # (3,) outward normal
    atom_index: int  # Closest atom
    area: float  # Associated surface area


def generate_surface(
    structure: Structure,
    probe_radius: float = 1.4,
    density: float = 3.0,
    msms_path: str = "msms",
) -> list[SurfacePoint]:
    """
    Generate molecular surface using MSMS.

    Args:
        structure: Input structure
        probe_radius: Probe radius in Angstroms
        density: Surface point density (points per Ų)
        msms_path: Path to MSMS executable

    Returns:
        List of surface points with coordinates and normals
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # Write XYZR file (x, y, z, radius for each atom)
        xyzr_path = tmpdir / "structure.xyzr"
        _write_xyzr(structure, xyzr_path)

        # Output paths
        output_prefix = tmpdir / "surface"

        # Run MSMS
        cmd = [
            msms_path,
            "-if", str(xyzr_path),
            "-of", str(output_prefix),
            "-probe_radius", str(probe_radius),
            "-density", str(density),
            "-no_header",
        ]

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=60,
                check=True,
            )
        except FileNotFoundError:
            raise RuntimeError(
                f"MSMS executable not found at '{msms_path}'. "
                "Please install MSMS and ensure it's in your PATH."
            )
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"MSMS failed: {e.stderr}")

        # Parse output files
        vert_path = tmpdir / "surface.vert"
        face_path = tmpdir / "surface.face"

        if not vert_path.exists():
            return []

        surface_points = _parse_msms_output(vert_path, face_path)

    return surface_points


def _write_xyzr(structure: Structure, path: Path) -> None:
    """Write XYZR file for MSMS input."""
    with open(path, "w") as f:
        for atom in structure.atoms:
            f.write(
                f"{atom.coords[0]:.3f} {atom.coords[1]:.3f} "
                f"{atom.coords[2]:.3f} {atom.vdw_radius:.2f}\n"
            )


def _parse_msms_output(
    vert_path: Path,
    face_path: Path,
) -> list[SurfacePoint]:
    """Parse MSMS vertex and face files."""
    surface_points: list[SurfacePoint] = []

    # Parse vertex file
    # Format: x y z nx ny nz atom_index ses_area
    with open(vert_path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 8:
                continue

            try:
                x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
                nx, ny, nz = float(parts[3]), float(parts[4]), float(parts[5])
                atom_idx = int(parts[6])
                area = float(parts[7])

                point = SurfacePoint(
                    coords=np.array([x, y, z]),
                    normal=np.array([nx, ny, nz]),
                    atom_index=atom_idx,
                    area=area,
                )
                surface_points.append(point)
            except (ValueError, IndexError):
                continue

    return surface_points


def check_msms_available(msms_path: str = "msms") -> bool:
    """Check if MSMS is available."""
    try:
        result = subprocess.run(
            [msms_path, "-h"],
            capture_output=True,
            timeout=5,
        )
        return True
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False
```

---

## 6. Analysis Pipeline

### 6.1 Pipeline Orchestration (`src/epitype/analysis/pipeline.py`)

```python
"""Full interface analysis pipeline."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable

from epitype.core.structure import Structure
from epitype.core.interface import InterfaceRegion, detect_interface, parse_chain_groups
from epitype.core.separation import separate_chains
from epitype.core.types import InterfaceMetrics
from epitype.io.parsers import parse_structure
from epitype.metrics.sasa import calculate_sasa, calculate_dsasa
from epitype.metrics.energy import calculate_binding_energy
from epitype.metrics.hbonds import analyze_hbonds
from epitype.metrics.shape import calculate_shape_complementarity
from epitype.metrics.packstat import calculate_packstat


@dataclass
class AnalysisConfig:
    """Configuration for interface analysis."""
    chain_spec: str  # e.g., "HL_A"
    interface_cutoff: float = 8.0
    probe_radius: float = 1.4
    separation_distance: float = 500.0
    minimize_energy: bool = True
    compute_energy: bool = True
    compute_shape: bool = True
    compute_packstat: bool = True


ProgressCallback = Callable[[str, float], None]


def analyze_interface(
    structure_path: str | Path,
    config: AnalysisConfig,
    progress_callback: ProgressCallback | None = None,
) -> InterfaceMetrics:
    """
    Run complete interface analysis pipeline.

    Args:
        structure_path: Path to PDB/mmCIF file
        config: Analysis configuration
        progress_callback: Optional callback for progress updates (step_name, fraction)

    Returns:
        InterfaceMetrics with all computed values
    """
    def report(step: str, progress: float) -> None:
        if progress_callback:
            progress_callback(step, progress)

    # Step 1: Parse structure
    report("Parsing structure", 0.0)
    structure = parse_structure(structure_path)

    # Step 2: Detect interface
    report("Detecting interface", 0.1)
    group1, group2 = parse_chain_groups(config.chain_spec)
    interface = detect_interface(
        structure,
        group1,
        group2,
        cutoff=config.interface_cutoff,
    )

    if not interface.residues_group1 or not interface.residues_group2:
        raise ValueError(
            f"No interface detected between {group1} and {group2}. "
            "Check chain specification and cutoff distance."
        )

    # Step 3: Create separated structure
    report("Separating chains", 0.2)
    separated = separate_chains(
        structure,
        interface,
        distance=config.separation_distance,
    )

    # Step 4: Calculate SASA
    report("Calculating SASA", 0.3)
    sasa_complex = calculate_sasa(structure, config.probe_radius)
    sasa_separated = calculate_sasa(separated, config.probe_radius)
    dsasa = sasa_separated.total - sasa_complex.total

    # Step 5: Calculate binding energy
    report("Calculating binding energy", 0.4)
    if config.compute_energy:
        dG = calculate_binding_energy(
            structure,
            separated,
            interface,
            minimize=config.minimize_energy,
        )
    else:
        dG = 0.0

    # Step 6: Analyze hydrogen bonds
    report("Analyzing hydrogen bonds", 0.6)
    hbond_result = analyze_hbonds(structure, interface)

    # Step 7: Calculate shape complementarity
    report("Calculating shape complementarity", 0.7)
    if config.compute_shape:
        sc_result = calculate_shape_complementarity(structure, interface)
        sc_value = sc_result.sc_value
    else:
        sc_value = 0.0

    # Step 8: Calculate PackStat
    report("Calculating PackStat", 0.85)
    if config.compute_packstat:
        packstat_result = calculate_packstat(structure, interface)
        packstat = packstat_result.packstat_interface
    else:
        packstat = 0.0

    report("Complete", 1.0)

    # Compute derived metrics
    dG_per_dSASA = dG / dsasa if dsasa > 0 else 0.0

    return InterfaceMetrics(
        dG_separated=dG,
        dSASA_int=dsasa,
        dG_per_dSASA=dG_per_dSASA,
        sc_value=sc_value,
        packstat=packstat,
        delta_unsat_hbonds=hbond_result.buried_unsatisfied_donors + hbond_result.buried_unsatisfied_acceptors,
        hbonds_int=hbond_result.interface_hbonds,
        interface_residues_1=interface.residue_ids_group1,
        interface_residues_2=interface.residue_ids_group2,
        sasa_complex=sasa_complex.total,
        sasa_separated=sasa_separated.total,
    )


def analyze_structure(
    structure: Structure,
    chain_spec: str,
    **kwargs,
) -> InterfaceMetrics:
    """
    Analyze a pre-loaded Structure object.

    Convenience function for API use.

    Args:
        structure: Loaded Structure object
        chain_spec: Chain specification (e.g., "HL_A")
        **kwargs: Additional config options

    Returns:
        InterfaceMetrics
    """
    config = AnalysisConfig(chain_spec=chain_spec, **kwargs)

    group1, group2 = parse_chain_groups(config.chain_spec)
    interface = detect_interface(structure, group1, group2, cutoff=config.interface_cutoff)

    separated = separate_chains(structure, interface, distance=config.separation_distance)

    sasa_complex = calculate_sasa(structure, config.probe_radius)
    sasa_separated = calculate_sasa(separated, config.probe_radius)
    dsasa = sasa_separated.total - sasa_complex.total

    if config.compute_energy:
        dG = calculate_binding_energy(structure, separated, interface, minimize=config.minimize_energy)
    else:
        dG = 0.0

    hbond_result = analyze_hbonds(structure, interface)

    if config.compute_shape:
        sc_result = calculate_shape_complementarity(structure, interface)
        sc_value = sc_result.sc_value
    else:
        sc_value = 0.0

    if config.compute_packstat:
        packstat_result = calculate_packstat(structure, interface)
        packstat = packstat_result.packstat_interface
    else:
        packstat = 0.0

    dG_per_dSASA = dG / dsasa if dsasa > 0 else 0.0

    return InterfaceMetrics(
        dG_separated=dG,
        dSASA_int=dsasa,
        dG_per_dSASA=dG_per_dSASA,
        sc_value=sc_value,
        packstat=packstat,
        delta_unsat_hbonds=hbond_result.buried_unsatisfied_donors + hbond_result.buried_unsatisfied_acceptors,
        hbonds_int=hbond_result.interface_hbonds,
        interface_residues_1=interface.residue_ids_group1,
        interface_residues_2=interface.residue_ids_group2,
        sasa_complex=sasa_complex.total,
        sasa_separated=sasa_separated.total,
    )
```

---

## 7. CLI Interface

### 7.1 Main CLI Entry Point (`src/epitype/cli/main.py`)

```python
"""Main CLI entry point using Typer."""
import typer
from rich.console import Console

from epitype.cli.run import run_app
from epitype.cli.sasa import sasa_app
from epitype.cli.energy import energy_app
from epitype.cli.hbonds import hbonds_app
from epitype.cli.shape import shape_app
from epitype.cli.packstat import packstat_app

app = typer.Typer(
    name="epitype",
    help="Protein-protein interface analysis toolkit.",
    no_args_is_help=True,
)

console = Console()

# Add subcommands
app.add_typer(run_app, name="run")
app.add_typer(sasa_app, name="sasa")
app.add_typer(energy_app, name="energy")
app.add_typer(hbonds_app, name="hbonds")
app.add_typer(shape_app, name="shape")
app.add_typer(packstat_app, name="packstat")


@app.command()
def version():
    """Show version information."""
    from epitype import __version__
    console.print(f"epitype v{__version__}")


@app.callback()
def main_callback():
    """Interface Analyzer: Compute protein-protein interface metrics."""
    pass


if __name__ == "__main__":
    app()
```

### 7.2 Run Command (`src/epitype/cli/run.py`)

```python
"""Main analysis run command."""
from pathlib import Path
from typing import Optional

import typer
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn

from epitype.analysis.pipeline import analyze_interface, AnalysisConfig
from epitype.io.writers import write_json, write_csv, format_metrics_table

run_app = typer.Typer(help="Run complete interface analysis.")
console = Console()


@run_app.callback(invoke_without_command=True)
def run(
    structure: Path = typer.Argument(
        ...,
        help="Path to PDB or mmCIF file",
        exists=True,
        readable=True,
    ),
    chains: str = typer.Option(
        ...,
        "--chains", "-c",
        help="Chain specification (e.g., 'HL_A' for chains H+L vs A)",
    ),
    output: Optional[Path] = typer.Option(
        None,
        "--output", "-o",
        help="Output file path (JSON or CSV based on extension)",
    ),
    cutoff: float = typer.Option(
        8.0,
        "--cutoff",
        help="Interface detection cutoff in Angstroms",
    ),
    no_energy: bool = typer.Option(
        False,
        "--no-energy",
        help="Skip energy calculation (faster)",
    ),
    no_shape: bool = typer.Option(
        False,
        "--no-shape",
        help="Skip shape complementarity (faster, no MSMS required)",
    ),
    no_packstat: bool = typer.Option(
        False,
        "--no-packstat",
        help="Skip PackStat calculation (faster)",
    ),
    quiet: bool = typer.Option(
        False,
        "--quiet", "-q",
        help="Suppress progress output",
    ),
):
    """
    Analyze protein-protein interface.

    Computes binding energy, buried surface area, shape complementarity,
    packing quality, and hydrogen bond metrics.

    Example:
        epitype run structure.pdb --chains HL_A -o results.json
    """
    config = AnalysisConfig(
        chain_spec=chains,
        interface_cutoff=cutoff,
        compute_energy=not no_energy,
        compute_shape=not no_shape,
        compute_packstat=not no_packstat,
    )

    if not quiet:
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            console=console,
        ) as progress:
            task = progress.add_task("Analyzing...", total=None)

            def progress_callback(step: str, frac: float):
                progress.update(task, description=f"{step}...")

            metrics = analyze_interface(structure, config, progress_callback)
    else:
        metrics = analyze_interface(structure, config)

    # Output results
    if output:
        if output.suffix == ".json":
            write_json(metrics, output, structure_name=structure.stem)
            console.print(f"Results written to {output}")
        elif output.suffix == ".csv":
            write_csv(metrics, output, structure_name=structure.stem)
            console.print(f"Results written to {output}")
        else:
            # Default to JSON
            write_json(metrics, output, structure_name=structure.stem)
            console.print(f"Results written to {output}")
    else:
        # Print to console
        console.print()
        console.print(format_metrics_table(metrics, structure_name=structure.name))
```

### 7.3 Subcommand Examples (`src/epitype/cli/sasa.py`)

```python
"""SASA calculation subcommand."""
from pathlib import Path
from typing import Optional

import typer
from rich.console import Console
from rich.table import Table

from epitype.io.parsers import parse_structure
from epitype.metrics.sasa import calculate_sasa

sasa_app = typer.Typer(help="Calculate solvent accessible surface area.")
console = Console()


@sasa_app.callback(invoke_without_command=True)
def sasa(
    structure: Path = typer.Argument(
        ...,
        help="Path to PDB or mmCIF file",
        exists=True,
    ),
    probe_radius: float = typer.Option(
        1.4,
        "--probe-radius", "-p",
        help="Probe radius in Angstroms",
    ),
    per_residue: bool = typer.Option(
        False,
        "--per-residue",
        help="Show per-residue SASA breakdown",
    ),
    output: Optional[Path] = typer.Option(
        None,
        "--output", "-o",
        help="Output file (JSON)",
    ),
):
    """
    Calculate SASA for a structure.

    Example:
        epitype sasa structure.pdb --per-residue
    """
    struct = parse_structure(structure)
    result = calculate_sasa(struct, probe_radius)

    if output:
        import json
        data = {
            "total": result.total,
            "polar": result.polar,
            "apolar": result.apolar,
        }
        if per_residue:
            data["per_residue"] = result.per_residue

        with open(output, "w") as f:
            json.dump(data, f, indent=2)
        console.print(f"Results written to {output}")
    else:
        table = Table(title=f"SASA: {structure.name}")
        table.add_column("Metric", style="cyan")
        table.add_column("Value (Ų)", style="green", justify="right")

        table.add_row("Total SASA", f"{result.total:.1f}")
        table.add_row("Polar SASA", f"{result.polar:.1f}")
        table.add_row("Apolar SASA", f"{result.apolar:.1f}")

        console.print(table)

        if per_residue:
            console.print("\nPer-residue SASA:")
            for res_id, sasa in sorted(result.per_residue.items())[:20]:
                console.print(f"  {res_id}: {sasa:.1f} Ų")
            if len(result.per_residue) > 20:
                console.print(f"  ... and {len(result.per_residue) - 20} more")
```

### 7.4 Common CLI Utilities (`src/epitype/cli/common.py`)

```python
"""Shared CLI utilities and options."""
from pathlib import Path
from typing import Annotated

import typer


# Common option types
StructurePath = Annotated[
    Path,
    typer.Argument(
        help="Path to PDB or mmCIF file",
        exists=True,
        readable=True,
    ),
]

ChainSpec = Annotated[
    str,
    typer.Option(
        "--chains", "-c",
        help="Chain specification (e.g., 'HL_A')",
    ),
]

OutputPath = Annotated[
    Path | None,
    typer.Option(
        "--output", "-o",
        help="Output file path",
    ),
]

Verbose = Annotated[
    bool,
    typer.Option(
        "--verbose", "-v",
        help="Enable verbose output",
    ),
]


def validate_chain_spec(chain_spec: str) -> tuple[list[str], list[str]]:
    """Validate and parse chain specification."""
    from epitype.core.interface import parse_chain_groups

    try:
        return parse_chain_groups(chain_spec)
    except ValueError as e:
        raise typer.BadParameter(str(e))
```

---

## 8. Python API

### 8.1 Package Exports (`src/epitype/__init__.py`)

```python
"""
Interface Analyzer: Protein-protein interface analysis toolkit.

Example usage:

    from epitype import analyze, load_structure

    # Full analysis
    metrics = analyze("complex.pdb", chains="HL_A")
    print(f"Binding energy: {metrics.dG_separated:.2f} REU")
    print(f"Buried surface: {metrics.dSASA_int:.1f} Ų")

    # Individual metrics
    from epitype import sasa, hbonds, shape

    structure = load_structure("complex.pdb")
    sasa_result = sasa.calculate(structure)
"""

__version__ = "0.1.0"

# Core data structures
from epitype.core.structure import Structure, Chain, Residue, Atom
from epitype.core.interface import InterfaceRegion, detect_interface, parse_chain_groups
from epitype.core.types import InterfaceMetrics

# High-level API
from epitype.io.parsers import parse_structure as load_structure
from epitype.analysis.pipeline import analyze_interface as analyze
from epitype.analysis.pipeline import analyze_structure, AnalysisConfig

# Metric modules (for individual access)
from epitype import metrics

__all__ = [
    # Version
    "__version__",
    # Data structures
    "Structure",
    "Chain",
    "Residue",
    "Atom",
    "InterfaceRegion",
    "InterfaceMetrics",
    # Functions
    "load_structure",
    "analyze",
    "analyze_structure",
    "detect_interface",
    "parse_chain_groups",
    # Config
    "AnalysisConfig",
    # Submodules
    "metrics",
]
```

### 8.2 Metrics Submodule (`src/epitype/metrics/__init__.py`)

```python
"""
Individual metric calculations.

Each metric can be computed independently:

    from epitype import load_structure, detect_interface
    from epitype.metrics import sasa, hbonds, shape, packstat, energy

    structure = load_structure("complex.pdb")
    interface = detect_interface(structure, ["H", "L"], ["A"])

    # SASA
    sasa_result = sasa.calculate_sasa(structure)
    dsasa = sasa.calculate_dsasa(complex_struct, separated_struct, interface)

    # Hydrogen bonds
    hbond_result = hbonds.analyze_hbonds(structure, interface)

    # Shape complementarity
    sc_result = shape.calculate_shape_complementarity(structure, interface)

    # PackStat
    packstat_result = packstat.calculate_packstat(structure, interface)
"""

from epitype.metrics import sasa
from epitype.metrics import energy
from epitype.metrics import hbonds
from epitype.metrics import shape
from epitype.metrics import packstat

__all__ = ["sasa", "energy", "hbonds", "shape", "packstat"]
```

### 8.3 API Usage Examples

```python
# Example 1: Full analysis pipeline
from epitype import analyze

metrics = analyze("1YY9.pdb", chains="HL_A")
print(f"dG_separated: {metrics.dG_separated:.2f} REU")
print(f"dSASA: {metrics.dSASA_int:.1f} Ų")
print(f"Shape complementarity: {metrics.sc_value:.3f}")


# Example 2: Custom configuration
from epitype import analyze, AnalysisConfig

config = AnalysisConfig(
    chain_spec="AB_C",
    interface_cutoff=10.0,
    compute_energy=False,  # Skip slow energy calculation
    compute_packstat=True,
)
metrics = analyze("complex.pdb", config)


# Example 3: Working with Structure objects
from epitype import (
    load_structure,
    detect_interface,
    analyze_structure,
)
from epitype.core.separation import separate_chains

structure = load_structure("complex.pdb")
interface = detect_interface(structure, ["H", "L"], ["A"])

print(f"Interface residues (Ab): {len(interface.residues_group1)}")
print(f"Interface residues (Ag): {len(interface.residues_group2)}")

# Analyze with the loaded structure
metrics = analyze_structure(structure, chain_spec="HL_A")


# Example 4: Individual metrics
from epitype import load_structure, detect_interface
from epitype.core.separation import separate_chains
from epitype.metrics.sasa import calculate_sasa, calculate_dsasa
from epitype.metrics.hbonds import analyze_hbonds

structure = load_structure("complex.pdb")
interface = detect_interface(structure, ["H", "L"], ["A"])
separated = separate_chains(structure, interface)

# Just SASA
sasa = calculate_sasa(structure)
print(f"Total SASA: {sasa.total:.1f} Ų")

dsasa = calculate_dsasa(structure, separated, interface)
print(f"Buried surface: {dsasa:.1f} Ų")

# Just H-bonds
hbonds = analyze_hbonds(structure, interface)
print(f"Interface H-bonds: {hbonds.interface_hbonds}")
```

---

## 9. Test Suite

### 9.1 Test Configuration (`tests/conftest.py`)

```python
"""Pytest fixtures and test configuration."""
from pathlib import Path
from typing import Generator
import tempfile
import shutil

import numpy as np
import pytest

from epitype.core.structure import Structure, Chain, Residue, Atom
from epitype.core.interface import InterfaceRegion


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
```

### 9.2 Unit Tests

#### 9.2.1 Structure Tests (`tests/unit/test_structure.py`)

```python
"""Unit tests for structure representation."""
import numpy as np
import pytest

from epitype.core.structure import Structure, Chain, Residue, Atom


class TestAtom:
    """Tests for Atom class."""

    def test_atom_creation(self):
        """Test basic atom creation."""
        atom = Atom(
            index=0,
            name="CA",
            element="C",
            coords=np.array([1.0, 2.0, 3.0]),
            residue_index=0,
            chain_id="A",
        )
        assert atom.name == "CA"
        assert atom.element == "C"
        np.testing.assert_array_equal(atom.coords, [1.0, 2.0, 3.0])

    def test_vdw_radius(self):
        """Test van der Waals radius lookup."""
        carbon = Atom(0, "CA", "C", np.zeros(3), 0, "A")
        nitrogen = Atom(1, "N", "N", np.zeros(3), 0, "A")
        oxygen = Atom(2, "O", "O", np.zeros(3), 0, "A")

        assert carbon.vdw_radius == 1.70
        assert nitrogen.vdw_radius == 1.55
        assert oxygen.vdw_radius == 1.52

    def test_is_backbone(self):
        """Test backbone atom detection."""
        ca = Atom(0, "CA", "C", np.zeros(3), 0, "A")
        cb = Atom(1, "CB", "C", np.zeros(3), 0, "A")

        assert ca.is_backbone is True
        assert cb.is_backbone is False

    def test_distance_to(self):
        """Test distance calculation."""
        a1 = Atom(0, "CA", "C", np.array([0.0, 0.0, 0.0]), 0, "A")
        a2 = Atom(1, "CA", "C", np.array([3.0, 4.0, 0.0]), 0, "A")

        assert a1.distance_to(a2) == pytest.approx(5.0)


class TestResidue:
    """Tests for Residue class."""

    def test_residue_creation(self):
        """Test basic residue creation."""
        residue = Residue(
            index=0,
            name="ALA",
            chain_id="A",
            seq_num=1,
        )
        assert residue.name == "ALA"
        assert residue.chain_id == "A"
        assert len(residue.atoms) == 0

    def test_ca_cb_access(self):
        """Test CA/CB atom access."""
        residue = Residue(0, "ALA", "A", 1)
        ca = Atom(0, "CA", "C", np.zeros(3), 0, "A")
        cb = Atom(1, "CB", "C", np.ones(3), 0, "A")
        residue.atoms = [ca, cb]

        assert residue.ca is ca
        assert residue.cb is cb

    def test_glycine_cb_fallback(self):
        """Test that CB falls back to CA for glycine."""
        residue = Residue(0, "GLY", "A", 1)
        ca = Atom(0, "CA", "C", np.zeros(3), 0, "A")
        residue.atoms = [ca]

        assert residue.cb is ca  # Should fall back to CA

    def test_full_id(self):
        """Test full residue identifier."""
        residue = Residue(0, "ALA", "A", 42, insertion_code="A")
        assert residue.full_id == "A:ALA:42A"


class TestChain:
    """Tests for Chain class."""

    def test_chain_atoms_iterator(self, minimal_structure: Structure):
        """Test atom iteration over chain."""
        chain = minimal_structure.chains["A"]
        atoms = list(chain.atoms)

        assert len(atoms) == 15  # 3 residues × 5 atoms

    def test_chain_coords(self, minimal_structure: Structure):
        """Test coordinate extraction."""
        chain = minimal_structure.chains["A"]
        coords = chain.coords

        assert coords.shape == (15, 3)


class TestStructure:
    """Tests for Structure class."""

    def test_structure_creation(self, minimal_structure: Structure):
        """Test structure creation."""
        assert minimal_structure.name == "minimal"
        assert len(minimal_structure.chains) == 2

    def test_num_atoms(self, minimal_structure: Structure):
        """Test atom counting."""
        assert minimal_structure.num_atoms == 30

    def test_num_residues(self, minimal_structure: Structure):
        """Test residue counting."""
        assert minimal_structure.num_residues == 6

    def test_subset(self, minimal_structure: Structure):
        """Test chain subsetting."""
        subset = minimal_structure.subset(["A"])

        assert len(subset.chains) == 1
        assert "A" in subset.chains
        assert "B" not in subset.chains

    def test_translate(self, minimal_structure: Structure):
        """Test structure translation."""
        original_coords = minimal_structure.coords.copy()
        translation = np.array([10.0, 20.0, 30.0])

        structure_copy = minimal_structure.copy()
        structure_copy.translate(translation)

        new_coords = structure_copy.coords
        expected = original_coords + translation

        np.testing.assert_array_almost_equal(new_coords, expected)

    def test_copy_is_deep(self, minimal_structure: Structure):
        """Test that copy is deep."""
        copy = minimal_structure.copy()

        # Modify original
        for atom in minimal_structure.atoms:
            atom.coords = np.array([999.0, 999.0, 999.0])
            break

        # Copy should be unaffected
        first_atom = next(iter(copy.atoms))
        assert first_atom.coords[0] != 999.0
```

#### 9.2.2 Interface Tests (`tests/unit/test_interface.py`)

```python
"""Unit tests for interface detection."""
import numpy as np
import pytest

from epitype.core.interface import (
    detect_interface,
    parse_chain_groups,
    InterfaceRegion,
)
from epitype.core.structure import Structure


class TestParseChainGroups:
    """Tests for chain group parsing."""

    def test_simple_spec(self):
        """Test simple chain specification."""
        g1, g2 = parse_chain_groups("A_B")
        assert g1 == ["A"]
        assert g2 == ["B"]

    def test_multi_chain_spec(self):
        """Test multi-chain specification."""
        g1, g2 = parse_chain_groups("HL_A")
        assert g1 == ["H", "L"]
        assert g2 == ["A"]

    def test_both_multi(self):
        """Test both groups multi-chain."""
        g1, g2 = parse_chain_groups("AB_CD")
        assert g1 == ["A", "B"]
        assert g2 == ["C", "D"]

    def test_invalid_no_underscore(self):
        """Test invalid spec without underscore."""
        with pytest.raises(ValueError, match="underscore"):
            parse_chain_groups("ABC")

    def test_invalid_multiple_underscores(self):
        """Test invalid spec with multiple underscores."""
        with pytest.raises(ValueError, match="one underscore"):
            parse_chain_groups("A_B_C")


class TestDetectInterface:
    """Tests for interface detection."""

    def test_detect_basic_interface(self, minimal_structure: Structure):
        """Test basic interface detection."""
        interface = detect_interface(
            minimal_structure,
            ["A"],
            ["B"],
            cutoff=10.0,
        )

        assert len(interface.residues_group1) > 0
        assert len(interface.residues_group2) > 0
        assert interface.group1_chains == ["A"]
        assert interface.group2_chains == ["B"]

    def test_no_interface_large_cutoff(self, minimal_structure: Structure):
        """Test that small cutoff finds no interface."""
        interface = detect_interface(
            minimal_structure,
            ["A"],
            ["B"],
            cutoff=1.0,  # Too small
        )

        # With 5 Å separation and 1 Å cutoff, should find nothing
        assert len(interface.residues_group1) == 0
        assert len(interface.residues_group2) == 0

    def test_contact_pairs(self, minimal_structure: Structure):
        """Test contact pair detection."""
        interface = detect_interface(
            minimal_structure,
            ["A"],
            ["B"],
            cutoff=10.0,
        )

        assert len(interface.contact_pairs) > 0
        for i, j in interface.contact_pairs:
            assert isinstance(i, int)
            assert isinstance(j, int)

    def test_atom_indices_populated(self, minimal_structure: Structure):
        """Test that atom indices are populated."""
        interface = detect_interface(
            minimal_structure,
            ["A"],
            ["B"],
            cutoff=10.0,
        )

        assert len(interface.atoms_group1_indices) > 0
        assert len(interface.atoms_group2_indices) > 0

    def test_use_ca_instead_of_cb(self, minimal_structure: Structure):
        """Test using CA for distance calculation."""
        interface_cb = detect_interface(
            minimal_structure, ["A"], ["B"], cutoff=10.0, use_cb=True
        )
        interface_ca = detect_interface(
            minimal_structure, ["A"], ["B"], cutoff=10.0, use_cb=False
        )

        # Both should work
        assert interface_cb.num_interface_residues >= 0
        assert interface_ca.num_interface_residues >= 0
```

#### 9.2.3 Separation Tests (`tests/unit/test_separation.py`)

```python
"""Unit tests for chain separation."""
import numpy as np
import pytest

from epitype.core.separation import (
    compute_interface_normal,
    separate_chains,
    compute_separation_distance,
)
from epitype.core.structure import Structure
from epitype.core.interface import InterfaceRegion, detect_interface


class TestComputeInterfaceNormal:
    """Tests for interface normal computation."""

    def test_normal_is_unit_vector(
        self, minimal_structure: Structure, interface_ab: InterfaceRegion
    ):
        """Test that normal is a unit vector."""
        normal = compute_interface_normal(minimal_structure, interface_ab)

        assert normal.shape == (3,)
        np.testing.assert_almost_equal(np.linalg.norm(normal), 1.0)

    def test_normal_direction(
        self, minimal_structure: Structure, interface_ab: InterfaceRegion
    ):
        """Test that normal points from group1 to group2."""
        normal = compute_interface_normal(minimal_structure, interface_ab)

        # In our minimal structure, B is offset in z direction
        # So normal should have positive z component
        assert normal[2] > 0


class TestSeparateChains:
    """Tests for chain separation."""

    def test_separation_increases_distance(
        self, minimal_structure: Structure, interface_ab: InterfaceRegion
    ):
        """Test that separation increases distance between groups."""
        from scipy.spatial.distance import cdist

        # Get original distances
        coords_a = np.array([a.coords for a in minimal_structure.chains["A"].atoms])
        coords_b = np.array([a.coords for a in minimal_structure.chains["B"].atoms])
        original_min = cdist(coords_a, coords_b).min()

        # Separate
        separated = separate_chains(minimal_structure, interface_ab, distance=100.0)

        # Get new distances
        coords_a_new = np.array([a.coords for a in separated.chains["A"].atoms])
        coords_b_new = np.array([a.coords for a in separated.chains["B"].atoms])
        new_min = cdist(coords_a_new, coords_b_new).min()

        assert new_min > original_min
        assert new_min > 100.0

    def test_separation_preserves_group1(
        self, minimal_structure: Structure, interface_ab: InterfaceRegion
    ):
        """Test that group1 coordinates are preserved."""
        original_coords = minimal_structure.chains["A"].coords.copy()

        separated = separate_chains(minimal_structure, interface_ab, distance=100.0)
        new_coords = separated.chains["A"].coords

        np.testing.assert_array_almost_equal(original_coords, new_coords)

    def test_separation_translates_group2(
        self, minimal_structure: Structure, interface_ab: InterfaceRegion
    ):
        """Test that group2 is translated."""
        original_coords = minimal_structure.chains["B"].coords.copy()

        separated = separate_chains(minimal_structure, interface_ab, distance=100.0)
        new_coords = separated.chains["B"].coords

        # Should not be equal
        assert not np.allclose(original_coords, new_coords)

    def test_original_unchanged(
        self, minimal_structure: Structure, interface_ab: InterfaceRegion
    ):
        """Test that original structure is unchanged."""
        original_coords = minimal_structure.coords.copy()

        _ = separate_chains(minimal_structure, interface_ab, distance=100.0)

        np.testing.assert_array_almost_equal(
            minimal_structure.coords, original_coords
        )
```

#### 9.2.4 SASA Tests (`tests/unit/test_sasa.py`)

```python
"""Unit tests for SASA calculations."""
import numpy as np
import pytest

from epitype.metrics.sasa import (
    calculate_sasa,
    calculate_dsasa,
    calculate_interface_sasa,
    SASAResult,
)
from epitype.core.structure import Structure
from epitype.core.interface import InterfaceRegion


class TestCalculateSASA:
    """Tests for SASA calculation."""

    def test_sasa_positive(self, minimal_structure: Structure):
        """Test that SASA is positive."""
        result = calculate_sasa(minimal_structure)

        assert result.total > 0
        assert result.polar >= 0
        assert result.apolar >= 0

    def test_sasa_components_sum(self, minimal_structure: Structure):
        """Test that polar + apolar approximately equals total."""
        result = calculate_sasa(minimal_structure)

        assert result.polar + result.apolar == pytest.approx(result.total, rel=0.01)

    def test_per_atom_sasa(self, minimal_structure: Structure):
        """Test per-atom SASA values."""
        result = calculate_sasa(minimal_structure)

        assert len(result.per_atom) == minimal_structure.num_atoms
        for atom_idx, sasa in result.per_atom.items():
            assert sasa >= 0

    def test_per_residue_sasa(self, minimal_structure: Structure):
        """Test per-residue SASA values."""
        result = calculate_sasa(minimal_structure)

        assert len(result.per_residue) > 0
        for res_id, sasa in result.per_residue.items():
            assert sasa >= 0

    def test_different_probe_radius(self, minimal_structure: Structure):
        """Test SASA changes with probe radius."""
        result_small = calculate_sasa(minimal_structure, probe_radius=1.0)
        result_large = calculate_sasa(minimal_structure, probe_radius=2.0)

        # Larger probe should give smaller SASA (can't access tight spaces)
        assert result_large.total < result_small.total


class TestCalculateDSASA:
    """Tests for dSASA (buried surface area) calculation."""

    def test_dsasa_positive(
        self,
        minimal_structure: Structure,
        interface_ab: InterfaceRegion,
    ):
        """Test that dSASA is positive for interacting chains."""
        from epitype.core.separation import separate_chains

        separated = separate_chains(minimal_structure, interface_ab, distance=500.0)
        dsasa = calculate_dsasa(
            minimal_structure, separated, interface_ab
        )

        # dSASA should be positive (surface is buried in complex)
        assert dsasa > 0

    def test_dsasa_zero_for_separated(
        self,
        minimal_structure: Structure,
        interface_ab: InterfaceRegion,
    ):
        """Test that dSASA is near zero for already separated chains."""
        from epitype.core.separation import separate_chains

        # Create very separated structure
        separated = separate_chains(minimal_structure, interface_ab, distance=1000.0)

        # dSASA between two already-separated structures should be ~0
        dsasa = calculate_dsasa(separated, separated, interface_ab)
        assert abs(dsasa) < 1.0  # Allow small numerical error
```

#### 9.2.5 Hydrogen Bond Tests (`tests/unit/test_hbonds.py`)

```python
"""Unit tests for hydrogen bond detection."""
import numpy as np
import pytest

from epitype.metrics.hbonds import (
    detect_hydrogen_bonds,
    count_interface_hbonds,
    analyze_hbonds,
    HydrogenBond,
    DONOR_ATOMS,
    ACCEPTOR_ATOMS,
)
from epitype.core.structure import Structure, Chain, Residue, Atom


@pytest.fixture
def hbond_structure() -> Structure:
    """Create structure with potential H-bonds."""
    structure = Structure(name="hbond_test")

    # Create two residues with N-H...O=C hydrogen bond geometry
    chain = Chain(chain_id="A")

    # Residue 1: donor (N)
    res1 = Residue(0, "ALA", "A", 1)
    res1.atoms = [
        Atom(0, "N", "N", np.array([0.0, 0.0, 0.0]), 0, "A"),
        Atom(1, "CA", "C", np.array([1.5, 0.0, 0.0]), 0, "A"),
        Atom(2, "C", "C", np.array([2.0, 1.5, 0.0]), 0, "A"),
        Atom(3, "O", "O", np.array([1.5, 2.5, 0.0]), 0, "A"),
    ]

    # Residue 2: acceptor (O) - positioned for H-bond
    res2 = Residue(1, "ALA", "A", 2)
    res2.atoms = [
        Atom(4, "N", "N", np.array([3.0, 1.5, 0.0]), 1, "A"),  # ~3Å from res1 O
        Atom(5, "CA", "C", np.array([4.5, 1.5, 0.0]), 1, "A"),
        Atom(6, "C", "C", np.array([5.0, 3.0, 0.0]), 1, "A"),
        Atom(7, "O", "O", np.array([4.5, 4.0, 0.0]), 1, "A"),
    ]

    chain.residues = [res1, res2]
    structure.chains["A"] = chain

    return structure


class TestDetectHydrogenBonds:
    """Tests for H-bond detection."""

    def test_detect_hbonds(self, hbond_structure: Structure):
        """Test basic H-bond detection."""
        hbonds = detect_hydrogen_bonds(hbond_structure)

        # Should find at least one H-bond
        assert len(hbonds) > 0

    def test_hbond_distance_cutoff(self, hbond_structure: Structure):
        """Test that distance cutoff is respected."""
        hbonds_tight = detect_hydrogen_bonds(hbond_structure, distance_cutoff=2.0)
        hbonds_loose = detect_hydrogen_bonds(hbond_structure, distance_cutoff=5.0)

        assert len(hbonds_tight) <= len(hbonds_loose)

    def test_hbond_structure(self, hbond_structure: Structure):
        """Test H-bond data structure."""
        hbonds = detect_hydrogen_bonds(hbond_structure)

        if hbonds:
            hb = hbonds[0]
            assert isinstance(hb, HydrogenBond)
            assert hb.donor_atom is not None
            assert hb.acceptor_atom is not None
            assert hb.distance > 0

    def test_no_same_residue_hbonds(self, hbond_structure: Structure):
        """Test that H-bonds within same residue are excluded."""
        hbonds = detect_hydrogen_bonds(hbond_structure)

        for hb in hbonds:
            assert hb.donor_atom.residue_index != hb.acceptor_atom.residue_index


class TestCountInterfaceHbonds:
    """Tests for interface H-bond counting."""

    def test_cross_interface_counting(
        self, minimal_structure: Structure, interface_ab
    ):
        """Test counting of cross-interface H-bonds."""
        hbonds = detect_hydrogen_bonds(minimal_structure)
        count = count_interface_hbonds(hbonds, interface_ab)

        # Count should be non-negative
        assert count >= 0
        assert count <= len(hbonds)


class TestAnalyzeHbonds:
    """Tests for complete H-bond analysis."""

    def test_analyze_returns_result(
        self, minimal_structure: Structure, interface_ab
    ):
        """Test that analysis returns proper result."""
        result = analyze_hbonds(minimal_structure, interface_ab)

        assert result.total_hbonds >= 0
        assert result.interface_hbonds >= 0
        assert result.buried_unsatisfied_donors >= 0
        assert result.buried_unsatisfied_acceptors >= 0
```

#### 9.2.6 Geometry Tests (`tests/unit/test_geometry.py`)

```python
"""Unit tests for geometry utilities."""
import numpy as np
import pytest

from epitype.utils.geometry import (
    normalize,
    angle_between,
    dihedral_angle,
    rotation_matrix,
)


class TestNormalize:
    """Tests for vector normalization."""

    def test_normalize_unit_vector(self):
        """Test normalizing already unit vector."""
        v = np.array([1.0, 0.0, 0.0])
        result = normalize(v)
        np.testing.assert_array_almost_equal(result, v)

    def test_normalize_arbitrary(self):
        """Test normalizing arbitrary vector."""
        v = np.array([3.0, 4.0, 0.0])
        result = normalize(v)

        assert np.linalg.norm(result) == pytest.approx(1.0)
        np.testing.assert_array_almost_equal(result, [0.6, 0.8, 0.0])

    def test_normalize_zero_vector(self):
        """Test normalizing zero vector."""
        v = np.array([0.0, 0.0, 0.0])
        result = normalize(v)

        np.testing.assert_array_equal(result, [0.0, 0.0, 0.0])


class TestAngleBetween:
    """Tests for angle calculation."""

    def test_parallel_vectors(self):
        """Test angle between parallel vectors."""
        v1 = np.array([1.0, 0.0, 0.0])
        v2 = np.array([2.0, 0.0, 0.0])

        angle = angle_between(v1, v2)
        assert angle == pytest.approx(0.0)

    def test_perpendicular_vectors(self):
        """Test angle between perpendicular vectors."""
        v1 = np.array([1.0, 0.0, 0.0])
        v2 = np.array([0.0, 1.0, 0.0])

        angle = angle_between(v1, v2)
        assert angle == pytest.approx(90.0)

    def test_antiparallel_vectors(self):
        """Test angle between antiparallel vectors."""
        v1 = np.array([1.0, 0.0, 0.0])
        v2 = np.array([-1.0, 0.0, 0.0])

        angle = angle_between(v1, v2)
        assert angle == pytest.approx(180.0)
```

#### 9.2.7 Parser Tests (`tests/unit/test_parsers.py`)

```python
"""Unit tests for file parsers."""
from pathlib import Path

import pytest

from epitype.io.parsers import (
    parse_pdb,
    parse_mmcif,
    parse_structure,
)


class TestParsePDB:
    """Tests for PDB parsing."""

    @pytest.mark.skipif(
        not Path("tests/data/1YY9.pdb").exists(),
        reason="Test PDB file not available"
    )
    def test_parse_1yy9(self, pdb_1yy9: Path):
        """Test parsing 1YY9 structure."""
        structure = parse_pdb(pdb_1yy9)

        assert structure.name == "1YY9"
        assert len(structure.chains) > 0
        assert structure.num_atoms > 0
        assert structure.num_residues > 0

    @pytest.mark.skipif(
        not Path("tests/data/1YY9.pdb").exists(),
        reason="Test PDB file not available"
    )
    def test_chain_ids_preserved(self, pdb_1yy9: Path):
        """Test that chain IDs are correctly parsed."""
        structure = parse_pdb(pdb_1yy9)

        # 1YY9 has chains H, L (antibody) and A (antigen)
        chain_ids = set(structure.chains.keys())
        assert "H" in chain_ids or "A" in chain_ids

    def test_parse_nonexistent_file(self):
        """Test parsing non-existent file raises error."""
        with pytest.raises(Exception):
            parse_pdb("/nonexistent/file.pdb")


class TestParseStructure:
    """Tests for auto-format detection."""

    @pytest.mark.skipif(
        not Path("tests/data/1YY9.pdb").exists(),
        reason="Test PDB file not available"
    )
    def test_auto_detect_pdb(self, pdb_1yy9: Path):
        """Test automatic PDB format detection."""
        structure = parse_structure(pdb_1yy9)
        assert structure.num_atoms > 0
```

### 9.3 Integration Tests

#### 9.3.1 Pipeline Integration (`tests/integration/test_pipeline.py`)

```python
"""Integration tests for analysis pipeline."""
from pathlib import Path

import pytest

from epitype.analysis.pipeline import (
    analyze_interface,
    analyze_structure,
    AnalysisConfig,
)
from epitype.core.types import InterfaceMetrics
from epitype.io.parsers import parse_structure


class TestAnalyzeInterface:
    """Tests for full analysis pipeline."""

    @pytest.mark.skipif(
        not Path("tests/data/1YY9.pdb").exists(),
        reason="Test PDB file not available"
    )
    def test_full_analysis_1yy9(self, pdb_1yy9: Path):
        """Test complete analysis of 1YY9."""
        config = AnalysisConfig(
            chain_spec="HL_A",
            compute_energy=False,  # Skip for speed
            compute_shape=False,
            compute_packstat=False,
        )

        metrics = analyze_interface(pdb_1yy9, config)

        assert isinstance(metrics, InterfaceMetrics)
        assert metrics.dSASA_int > 0
        assert metrics.hbonds_int >= 0

    @pytest.mark.skipif(
        not Path("tests/data/1YY9.pdb").exists(),
        reason="Test PDB file not available"
    )
    def test_analysis_with_all_metrics(self, pdb_1yy9: Path):
        """Test analysis with all metrics enabled."""
        config = AnalysisConfig(
            chain_spec="HL_A",
            compute_energy=False,  # Still skip energy (requires OpenMM setup)
            compute_shape=False,  # Skip (requires MSMS)
            compute_packstat=True,
        )

        metrics = analyze_interface(pdb_1yy9, config)

        assert metrics.packstat >= 0
        assert metrics.packstat <= 1

    def test_invalid_chain_spec(self, pdb_1yy9: Path):
        """Test error handling for invalid chain spec."""
        config = AnalysisConfig(chain_spec="XY_Z")

        # Should raise error (chains don't exist)
        with pytest.raises(ValueError):
            analyze_interface(pdb_1yy9, config)

    def test_progress_callback(self, pdb_1yy9: Path):
        """Test that progress callback is called."""
        config = AnalysisConfig(
            chain_spec="HL_A",
            compute_energy=False,
            compute_shape=False,
            compute_packstat=False,
        )

        steps_received = []

        def callback(step: str, progress: float):
            steps_received.append(step)

        analyze_interface(pdb_1yy9, config, progress_callback=callback)

        assert len(steps_received) > 0
        assert "Complete" in steps_received


class TestAnalyzeStructure:
    """Tests for Structure-based analysis."""

    @pytest.mark.skipif(
        not Path("tests/data/1YY9.pdb").exists(),
        reason="Test PDB file not available"
    )
    def test_analyze_preloaded_structure(self, pdb_1yy9: Path):
        """Test analysis with pre-loaded structure."""
        structure = parse_structure(pdb_1yy9)

        metrics = analyze_structure(
            structure,
            chain_spec="HL_A",
            compute_energy=False,
            compute_shape=False,
            compute_packstat=False,
        )

        assert isinstance(metrics, InterfaceMetrics)
        assert metrics.dSASA_int > 0
```

#### 9.3.2 CLI Integration (`tests/integration/test_cli.py`)

```python
"""Integration tests for CLI."""
from pathlib import Path
import json

import pytest
from typer.testing import CliRunner

from epitype.cli.main import app


runner = CliRunner()


class TestCLIRun:
    """Tests for 'epitype run' command."""

    @pytest.mark.skipif(
        not Path("tests/data/1YY9.pdb").exists(),
        reason="Test PDB file not available"
    )
    def test_run_basic(self, pdb_1yy9: Path):
        """Test basic run command."""
        result = runner.invoke(
            app,
            [
                "run",
                str(pdb_1yy9),
                "--chains", "HL_A",
                "--no-energy",
                "--no-shape",
                "--no-packstat",
            ],
        )

        assert result.exit_code == 0
        assert "dG_separated" in result.stdout or "dSASA" in result.stdout

    @pytest.mark.skipif(
        not Path("tests/data/1YY9.pdb").exists(),
        reason="Test PDB file not available"
    )
    def test_run_json_output(self, pdb_1yy9: Path, temp_dir: Path):
        """Test JSON output."""
        output_file = temp_dir / "results.json"

        result = runner.invoke(
            app,
            [
                "run",
                str(pdb_1yy9),
                "--chains", "HL_A",
                "--output", str(output_file),
                "--no-energy",
                "--no-shape",
                "--no-packstat",
                "--quiet",
            ],
        )

        assert result.exit_code == 0
        assert output_file.exists()

        with open(output_file) as f:
            data = json.load(f)

        assert "metrics" in data
        assert "dSASA_int" in data["metrics"]

    @pytest.mark.skipif(
        not Path("tests/data/1YY9.pdb").exists(),
        reason="Test PDB file not available"
    )
    def test_run_csv_output(self, pdb_1yy9: Path, temp_dir: Path):
        """Test CSV output."""
        output_file = temp_dir / "results.csv"

        result = runner.invoke(
            app,
            [
                "run",
                str(pdb_1yy9),
                "--chains", "HL_A",
                "--output", str(output_file),
                "--no-energy",
                "--no-shape",
                "--no-packstat",
                "--quiet",
            ],
        )

        assert result.exit_code == 0
        assert output_file.exists()

        content = output_file.read_text()
        assert "dSASA_int" in content


class TestCLISASA:
    """Tests for 'epitype sasa' command."""

    @pytest.mark.skipif(
        not Path("tests/data/1YY9.pdb").exists(),
        reason="Test PDB file not available"
    )
    def test_sasa_basic(self, pdb_1yy9: Path):
        """Test basic SASA command."""
        result = runner.invoke(app, ["sasa", str(pdb_1yy9)])

        assert result.exit_code == 0
        assert "Total SASA" in result.stdout


class TestCLIVersion:
    """Tests for version command."""

    def test_version(self):
        """Test version command."""
        result = runner.invoke(app, ["version"])

        assert result.exit_code == 0
        assert "epitype" in result.stdout
```

#### 9.3.3 API Integration (`tests/integration/test_api.py`)

```python
"""Integration tests for Python API."""
from pathlib import Path

import pytest

import epitype
from epitype import (
    analyze,
    load_structure,
    detect_interface,
    AnalysisConfig,
)


class TestHighLevelAPI:
    """Tests for high-level Python API."""

    def test_version_available(self):
        """Test version is accessible."""
        assert hasattr(epitype, "__version__")
        assert isinstance(epitype.__version__, str)

    @pytest.mark.skipif(
        not Path("tests/data/1YY9.pdb").exists(),
        reason="Test PDB file not available"
    )
    def test_analyze_function(self, pdb_1yy9: Path):
        """Test analyze() convenience function."""
        config = AnalysisConfig(
            chain_spec="HL_A",
            compute_energy=False,
            compute_shape=False,
            compute_packstat=False,
        )

        metrics = analyze(pdb_1yy9, config)

        assert metrics.dSASA_int > 0

    @pytest.mark.skipif(
        not Path("tests/data/1YY9.pdb").exists(),
        reason="Test PDB file not available"
    )
    def test_load_structure_function(self, pdb_1yy9: Path):
        """Test load_structure() function."""
        structure = load_structure(pdb_1yy9)

        assert structure.num_atoms > 0
        assert structure.num_residues > 0

    @pytest.mark.skipif(
        not Path("tests/data/1YY9.pdb").exists(),
        reason="Test PDB file not available"
    )
    def test_detect_interface_function(self, pdb_1yy9: Path):
        """Test detect_interface() function."""
        structure = load_structure(pdb_1yy9)
        interface = detect_interface(structure, ["H", "L"], ["A"])

        assert len(interface.residues_group1) > 0
        assert len(interface.residues_group2) > 0


class TestMetricsSubmodule:
    """Tests for metrics submodule access."""

    @pytest.mark.skipif(
        not Path("tests/data/1YY9.pdb").exists(),
        reason="Test PDB file not available"
    )
    def test_sasa_submodule(self, pdb_1yy9: Path):
        """Test accessing SASA through metrics submodule."""
        from epitype.metrics import sasa

        structure = load_structure(pdb_1yy9)
        result = sasa.calculate_sasa(structure)

        assert result.total > 0

    @pytest.mark.skipif(
        not Path("tests/data/1YY9.pdb").exists(),
        reason="Test PDB file not available"
    )
    def test_hbonds_submodule(self, pdb_1yy9: Path):
        """Test accessing H-bonds through metrics submodule."""
        from epitype.metrics import hbonds

        structure = load_structure(pdb_1yy9)
        interface = detect_interface(structure, ["H", "L"], ["A"])

        result = hbonds.analyze_hbonds(structure, interface)

        assert result.total_hbonds >= 0
```

### 9.4 End-to-End Tests

#### 9.4.1 Full Analysis Tests (`tests/e2e/test_full_analysis.py`)

```python
"""End-to-end tests for complete analysis workflow."""
from pathlib import Path
import json

import pytest

from epitype import analyze, AnalysisConfig
from epitype.io.writers import write_json, write_csv


class TestFullAnalysisWorkflow:
    """E2E tests for complete analysis workflow."""

    @pytest.mark.skipif(
        not Path("tests/data/1YY9.pdb").exists(),
        reason="Test PDB file not available"
    )
    def test_complete_workflow(self, pdb_1yy9: Path, temp_dir: Path):
        """Test complete analysis workflow from file to output."""
        # Step 1: Run analysis
        config = AnalysisConfig(
            chain_spec="HL_A",
            interface_cutoff=8.0,
            compute_energy=False,
            compute_shape=False,
            compute_packstat=True,
        )

        metrics = analyze(pdb_1yy9, config)

        # Step 2: Verify metrics are reasonable
        assert metrics.dSASA_int > 500  # Should have significant buried surface
        assert metrics.dSASA_int < 5000  # But not unreasonably large
        assert metrics.hbonds_int >= 0
        assert 0 <= metrics.packstat <= 1

        # Step 3: Write outputs
        json_path = temp_dir / "results.json"
        csv_path = temp_dir / "results.csv"

        write_json(metrics, json_path, structure_name="1YY9")
        write_csv(metrics, csv_path, structure_name="1YY9")

        # Step 4: Verify outputs are readable
        with open(json_path) as f:
            data = json.load(f)

        assert data["structure"] == "1YY9"
        assert data["metrics"]["dSASA_int"] == metrics.dSASA_int

        csv_content = csv_path.read_text()
        assert "1YY9" in csv_content
        assert str(int(metrics.dSASA_int)) in csv_content or f"{metrics.dSASA_int:.1f}" in csv_content


class TestKnownStructureMetrics:
    """E2E tests validating against known expected values."""

    @pytest.mark.skipif(
        not Path("tests/data/1YY9.pdb").exists(),
        reason="Test PDB file not available"
    )
    def test_1yy9_dsasa_range(self, pdb_1yy9: Path):
        """Test that 1YY9 dSASA is in expected range."""
        config = AnalysisConfig(
            chain_spec="HL_A",
            compute_energy=False,
            compute_shape=False,
            compute_packstat=False,
        )

        metrics = analyze(pdb_1yy9, config)

        # Expected range for antibody-antigen interface
        assert 1800 < metrics.dSASA_int < 2200

    @pytest.mark.skipif(
        not Path("tests/data/1YY9.pdb").exists(),
        reason="Test PDB file not available"
    )
    def test_1yy9_hbonds_range(self, pdb_1yy9: Path):
        """Test that 1YY9 interface H-bonds are in expected range."""
        config = AnalysisConfig(
            chain_spec="HL_A",
            compute_energy=False,
            compute_shape=False,
            compute_packstat=False,
        )

        metrics = analyze(pdb_1yy9, config)

        # Typical antibody-antigen has 8-20 interface H-bonds
        assert 5 <= metrics.hbonds_int <= 25

    @pytest.mark.skipif(
        not Path("tests/data/4FQI.pdb").exists(),
        reason="Test PDB file not available"
    )
    def test_4fqi_comparison(self, pdb_4fqi: Path):
        """Test 4FQI metrics."""
        config = AnalysisConfig(
            chain_spec="HL_A",
            compute_energy=False,
            compute_shape=False,
            compute_packstat=False,
        )

        metrics = analyze(pdb_4fqi, config)

        assert metrics.dSASA_int > 0
        assert len(metrics.interface_residues_1) > 0
        assert len(metrics.interface_residues_2) > 0
```

### 9.5 Test Data Setup

#### 9.5.1 Test Data README (`tests/data/README.md`)

```markdown
# Test Data

This directory contains test structures for the epitype test suite.

## Required Files

Download the following PDB files from RCSB:

1. **1YY9.pdb** - Antibody-antigen complex (hen egg lysozyme with HyHEL-63 Fab)
   - Download: https://files.rcsb.org/download/1YY9.pdb
   - Chains: H, L (antibody), A (antigen)
   - Interface: ~2000 Ų buried surface area

2. **4FQI.pdb** - Antibody-antigen complex
   - Download: https://files.rcsb.org/download/4FQI.pdb
   - Alternative test case

## Setup Script

Run the following to download test data:

```bash
cd tests/data
curl -O https://files.rcsb.org/download/1YY9.pdb
curl -O https://files.rcsb.org/download/4FQI.pdb
```

## Expected Values

Expected metric ranges for validation:

| Structure | dSASA (Ų) | Sc | H-bonds |
|-----------|-----------|-----|---------|
| 1YY9 | 1800-2200 | 0.60-0.75 | 8-20 |
| 4FQI | 1500-2000 | 0.55-0.70 | 5-15 |
```

---

## 10. Implementation Order

### Phase 1: Project Setup and Core Infrastructure

1. **Create directory structure and configuration**
   - Create all directories as shown in Section 1
   - Set up `pyproject.toml` with dependencies
   - Create `py.typed` marker file
   - Initialize `__init__.py` files

2. **Implement core types and constants**
   - `src/epitype/core/types.py`

3. **Implement structure representation**
   - `src/epitype/core/structure.py`
   - Write unit tests: `tests/unit/test_structure.py`

4. **Implement I/O parsers**
   - `src/epitype/io/parsers.py`
   - Write unit tests: `tests/unit/test_parsers.py`
   - Download test PDB files

### Phase 2: Interface Detection and Separation

5. **Implement interface detection**
   - `src/epitype/core/interface.py`
   - Write unit tests: `tests/unit/test_interface.py`

6. **Implement geometry utilities**
   - `src/epitype/utils/geometry.py`
   - Write unit tests: `tests/unit/test_geometry.py`

7. **Implement chain separation**
   - `src/epitype/core/separation.py`
   - Write unit tests: `tests/unit/test_separation.py`

### Phase 3: Core Metrics

8. **Implement SASA calculation**
   - `src/epitype/metrics/sasa.py`
   - Write unit tests: `tests/unit/test_sasa.py`

9. **Implement hydrogen bond detection**
   - `src/epitype/metrics/hbonds.py`
   - Write unit tests: `tests/unit/test_hbonds.py`

10. **Implement PackStat**
    - `src/epitype/metrics/packstat.py`
    - Write unit tests: `tests/unit/test_packstat.py`

### Phase 4: Advanced Metrics (External Dependencies)

11. **Implement MSMS wrapper**
    - `src/epitype/surface/msms.py`
    - `src/epitype/surface/types.py`

12. **Implement shape complementarity**
    - `src/epitype/metrics/shape.py`
    - Write unit tests: `tests/unit/test_shape.py`

13. **Implement binding energy (OpenMM)**
    - `src/epitype/metrics/energy.py`
    - Write unit tests: `tests/unit/test_energy.py`

### Phase 5: Pipeline and CLI

14. **Implement analysis pipeline**
    - `src/epitype/analysis/pipeline.py`
    - Write integration tests: `tests/integration/test_pipeline.py`

15. **Implement output writers**
    - `src/epitype/io/writers.py`

16. **Implement CLI**
    - `src/epitype/cli/main.py`
    - `src/epitype/cli/run.py`
    - `src/epitype/cli/sasa.py`
    - `src/epitype/cli/energy.py`
    - `src/epitype/cli/hbonds.py`
    - `src/epitype/cli/shape.py`
    - `src/epitype/cli/packstat.py`
    - `src/epitype/cli/common.py`
    - Write CLI tests: `tests/integration/test_cli.py`

### Phase 6: API and Documentation

17. **Finalize Python API exports**
    - `src/epitype/__init__.py`
    - `src/epitype/metrics/__init__.py`
    - Write API tests: `tests/integration/test_api.py`

18. **End-to-end tests**
    - `tests/e2e/test_full_analysis.py`
    - `tests/e2e/test_known_structures.py`

19. **Documentation**
    - `README.md`
    - `docs/` directory

---

## 11. External Dependencies Setup

### 11.1 MSMS Installation

MSMS (Michel Sanner's Molecular Surface) is required for shape complementarity calculations.

**Linux:**
```bash
# Download from https://ccsb.scripps.edu/msms/downloads/
wget https://ccsb.scripps.edu/msms/download/933/msms_i86_64Linux2_2.6.1.tar.gz
tar -xzf msms_i86_64Linux2_2.6.1.tar.gz
sudo mv msms.x86_64Linux2.2.6.1 /usr/local/bin/msms
sudo chmod +x /usr/local/bin/msms
```

**macOS:**
```bash
# Download from https://ccsb.scripps.edu/msms/downloads/
# Choose appropriate version for your Mac architecture
```

**Verify installation:**
```bash
msms -h
```

### 11.2 OpenMM Installation

OpenMM is required for binding energy calculations.

```bash
# Via conda (recommended)
conda install -c conda-forge openmm

# Or via pip (may require additional setup)
pip install openmm
```

**Verify installation:**
```python
import openmm
print(openmm.__version__)
```

### 11.3 FreeSASA Installation

FreeSASA is required for SASA calculations.

```bash
pip install freesasa
```

---

## 12. Verification and Testing

### 12.1 Running Tests

```bash
# Install dev dependencies
pip install -e ".[dev]"

# Download test data
cd tests/data
curl -O https://files.rcsb.org/download/1YY9.pdb
curl -O https://files.rcsb.org/download/4FQI.pdb
cd ../..

# Run all tests
pytest

# Run with coverage
pytest --cov=epitype --cov-report=html

# Run specific test categories
pytest tests/unit/          # Unit tests only
pytest tests/integration/   # Integration tests only
pytest tests/e2e/           # End-to-end tests only

# Run tests for specific module
pytest tests/unit/test_sasa.py -v
```

### 12.2 Manual Verification

**CLI verification:**
```bash
# Install package
pip install -e .

# Check version
epitype version

# Run full analysis
epitype run tests/data/1YY9.pdb --chains HL_A -o results.json

# Run individual metrics
epitype sasa tests/data/1YY9.pdb
epitype hbonds tests/data/1YY9.pdb --chains HL_A
```

**Python API verification:**
```python
from epitype import analyze, load_structure, AnalysisConfig

# Full analysis
config = AnalysisConfig(
    chain_spec="HL_A",
    compute_energy=False,
    compute_shape=False,
)
metrics = analyze("tests/data/1YY9.pdb", config)
print(f"dSASA: {metrics.dSASA_int:.1f} Ų")
print(f"H-bonds: {metrics.hbonds_int}")

# Expected output for 1YY9:
# dSASA: ~2000 Ų
# H-bonds: ~10-15
```

### 12.3 Validation Against Expected Values

| Test Case | Metric | Expected Range | Validation Command |
|-----------|--------|----------------|-------------------|
| 1YY9 | dSASA | 1800-2200 Ų | `pytest tests/e2e/test_full_analysis.py::TestKnownStructureMetrics::test_1yy9_dsasa_range` |
| 1YY9 | H-bonds | 8-20 | `pytest tests/e2e/test_full_analysis.py::TestKnownStructureMetrics::test_1yy9_hbonds_range` |
| 4FQI | dSASA | 1500-2000 Ų | `pytest tests/e2e/test_full_analysis.py::TestKnownStructureMetrics::test_4fqi_comparison` |

### 12.4 Type Checking

```bash
# Run mypy
mypy src/epitype

# Expected: no errors
```

### 12.5 Linting

```bash
# Run ruff
ruff check src/epitype tests

# Auto-fix issues
ruff check --fix src/epitype tests
```

---

## 13. Additional CLI Subcommands

The following CLI subcommand files need to be created following the pattern shown in `sasa.py`:

### 13.1 Energy Subcommand (`src/epitype/cli/energy.py`)

```python
"""Binding energy calculation subcommand."""
from pathlib import Path
from typing import Optional

import typer
from rich.console import Console

from epitype.io.parsers import parse_structure
from epitype.core.interface import detect_interface, parse_chain_groups
from epitype.core.separation import separate_chains
from epitype.metrics.energy import calculate_binding_energy

energy_app = typer.Typer(help="Calculate binding energy.")
console = Console()


@energy_app.callback(invoke_without_command=True)
def energy(
    structure: Path = typer.Argument(..., help="Path to PDB/mmCIF file", exists=True),
    chains: str = typer.Option(..., "--chains", "-c", help="Chain spec (e.g., HL_A)"),
    minimize: bool = typer.Option(True, "--minimize/--no-minimize", help="Energy minimize"),
    output: Optional[Path] = typer.Option(None, "--output", "-o"),
):
    """Calculate binding energy using OpenMM."""
    struct = parse_structure(structure)
    g1, g2 = parse_chain_groups(chains)
    interface = detect_interface(struct, g1, g2)
    separated = separate_chains(struct, interface)

    dG = calculate_binding_energy(struct, separated, interface, minimize=minimize)

    if output:
        import json
        with open(output, "w") as f:
            json.dump({"dG_separated": dG, "unit": "REU"}, f, indent=2)
    else:
        quality = "favorable" if dG < 0 else "unfavorable"
        console.print(f"Binding energy: {dG:.2f} REU ({quality})")
```

### 13.2 H-bonds Subcommand (`src/epitype/cli/hbonds.py`)

```python
"""Hydrogen bond analysis subcommand."""
from pathlib import Path
from typing import Optional

import typer
from rich.console import Console
from rich.table import Table

from epitype.io.parsers import parse_structure
from epitype.core.interface import detect_interface, parse_chain_groups
from epitype.metrics.hbonds import analyze_hbonds

hbonds_app = typer.Typer(help="Analyze hydrogen bonds.")
console = Console()


@hbonds_app.callback(invoke_without_command=True)
def hbonds(
    structure: Path = typer.Argument(..., help="Path to PDB/mmCIF file", exists=True),
    chains: str = typer.Option(..., "--chains", "-c", help="Chain spec (e.g., HL_A)"),
    output: Optional[Path] = typer.Option(None, "--output", "-o"),
):
    """Analyze hydrogen bonds at the interface."""
    struct = parse_structure(structure)
    g1, g2 = parse_chain_groups(chains)
    interface = detect_interface(struct, g1, g2)

    result = analyze_hbonds(struct, interface)

    if output:
        import json
        with open(output, "w") as f:
            json.dump({
                "total_hbonds": result.total_hbonds,
                "interface_hbonds": result.interface_hbonds,
                "buried_unsatisfied_donors": result.buried_unsatisfied_donors,
                "buried_unsatisfied_acceptors": result.buried_unsatisfied_acceptors,
            }, f, indent=2)
    else:
        table = Table(title="Hydrogen Bond Analysis")
        table.add_column("Metric", style="cyan")
        table.add_column("Value", justify="right")

        table.add_row("Total H-bonds", str(result.total_hbonds))
        table.add_row("Interface H-bonds", str(result.interface_hbonds))
        table.add_row("Unsatisfied donors", str(result.buried_unsatisfied_donors))
        table.add_row("Unsatisfied acceptors", str(result.buried_unsatisfied_acceptors))

        console.print(table)
```

### 13.3 Shape Subcommand (`src/epitype/cli/shape.py`)

```python
"""Shape complementarity subcommand."""
from pathlib import Path
from typing import Optional

import typer
from rich.console import Console

from epitype.io.parsers import parse_structure
from epitype.core.interface import detect_interface, parse_chain_groups
from epitype.metrics.shape import calculate_shape_complementarity

shape_app = typer.Typer(help="Calculate shape complementarity.")
console = Console()


@shape_app.callback(invoke_without_command=True)
def shape(
    structure: Path = typer.Argument(..., help="Path to PDB/mmCIF file", exists=True),
    chains: str = typer.Option(..., "--chains", "-c", help="Chain spec (e.g., HL_A)"),
    output: Optional[Path] = typer.Option(None, "--output", "-o"),
):
    """Calculate shape complementarity (requires MSMS)."""
    struct = parse_structure(structure)
    g1, g2 = parse_chain_groups(chains)
    interface = detect_interface(struct, g1, g2)

    result = calculate_shape_complementarity(struct, interface)

    if output:
        import json
        with open(output, "w") as f:
            json.dump({
                "sc_value": result.sc_value,
                "sc_median": result.sc_median,
                "sc_group1": result.sc_group1,
                "sc_group2": result.sc_group2,
            }, f, indent=2)
    else:
        quality = "good" if result.sc_value > 0.65 else "moderate" if result.sc_value > 0.5 else "poor"
        console.print(f"Shape complementarity: {result.sc_value:.3f} ({quality})")
```

### 13.4 PackStat Subcommand (`src/epitype/cli/packstat.py`)

```python
"""PackStat packing quality subcommand."""
from pathlib import Path
from typing import Optional

import typer
from rich.console import Console

from epitype.io.parsers import parse_structure
from epitype.core.interface import detect_interface, parse_chain_groups
from epitype.metrics.packstat import calculate_packstat

packstat_app = typer.Typer(help="Calculate packing quality (PackStat).")
console = Console()


@packstat_app.callback(invoke_without_command=True)
def packstat(
    structure: Path = typer.Argument(..., help="Path to PDB/mmCIF file", exists=True),
    chains: Optional[str] = typer.Option(None, "--chains", "-c", help="Chain spec for interface"),
    output: Optional[Path] = typer.Option(None, "--output", "-o"),
):
    """Calculate PackStat packing quality score."""
    struct = parse_structure(structure)

    interface = None
    if chains:
        g1, g2 = parse_chain_groups(chains)
        interface = detect_interface(struct, g1, g2)

    result = calculate_packstat(struct, interface)

    if output:
        import json
        with open(output, "w") as f:
            json.dump({
                "packstat": result.packstat,
                "packstat_interface": result.packstat_interface,
                "cavity_volume": result.cavity_volume,
            }, f, indent=2)
    else:
        quality = "good" if result.packstat > 0.65 else "moderate" if result.packstat > 0.5 else "poor"
        console.print(f"PackStat: {result.packstat:.3f} ({quality})")
        if interface:
            console.print(f"Interface PackStat: {result.packstat_interface:.3f}")
```

---

## 14. Geometry Utilities

### 14.1 Geometry Module (`src/epitype/utils/geometry.py`)

```python
"""Geometry utility functions."""
from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def normalize(v: NDArray[np.float64]) -> NDArray[np.float64]:
    """Normalize vector to unit length."""
    norm = np.linalg.norm(v)
    if norm < 1e-10:
        return np.zeros_like(v)
    return v / norm


def angle_between(v1: NDArray[np.float64], v2: NDArray[np.float64]) -> float:
    """Calculate angle between two vectors in degrees."""
    v1_n = normalize(v1)
    v2_n = normalize(v2)

    dot = np.clip(np.dot(v1_n, v2_n), -1.0, 1.0)
    return float(np.degrees(np.arccos(dot)))


def dihedral_angle(
    p1: NDArray[np.float64],
    p2: NDArray[np.float64],
    p3: NDArray[np.float64],
    p4: NDArray[np.float64],
) -> float:
    """Calculate dihedral angle between four points in degrees."""
    b1 = p2 - p1
    b2 = p3 - p2
    b3 = p4 - p3

    n1 = normalize(np.cross(b1, b2))
    n2 = normalize(np.cross(b2, b3))

    m1 = np.cross(n1, normalize(b2))

    x = np.dot(n1, n2)
    y = np.dot(m1, n2)

    return float(np.degrees(np.arctan2(y, x)))


def rotation_matrix(
    axis: NDArray[np.float64],
    angle: float,
) -> NDArray[np.float64]:
    """Create rotation matrix for rotation around axis by angle (radians)."""
    axis = normalize(axis)
    c = np.cos(angle)
    s = np.sin(angle)
    t = 1 - c

    x, y, z = axis

    return np.array([
        [t*x*x + c,   t*x*y - s*z, t*x*z + s*y],
        [t*x*y + s*z, t*y*y + c,   t*y*z - s*x],
        [t*x*z - s*y, t*y*z + s*x, t*z*z + c],
    ])


def centroid(coords: NDArray[np.float64]) -> NDArray[np.float64]:
    """Calculate centroid of coordinates."""
    return coords.mean(axis=0)


def rmsd(
    coords1: NDArray[np.float64],
    coords2: NDArray[np.float64],
) -> float:
    """Calculate RMSD between two coordinate sets."""
    diff = coords1 - coords2
    return float(np.sqrt(np.mean(np.sum(diff ** 2, axis=1))))
```

---

## Summary

This implementation plan provides a complete blueprint for building a Python-based InterfaceAnalyzer with:

1. **Modern Python package structure** installable via pip/uv
2. **CLI interface** using Typer with main pipeline (`run`) and individual subtools
3. **Python API** for programmatic access to all functionality
4. **Comprehensive test suite** with unit, integration, and e2e tests

**Key external dependencies:**
- FreeSASA (SASA calculations)
- OpenMM (binding energy with AMBER ff14SB)
- MSMS (molecular surface generation for shape complementarity)

**Implementation follows 6 phases**, starting with core infrastructure and progressing to advanced metrics and the final CLI/API layer.

Total estimated files: ~35 Python modules + tests
