# Epitype

A Python toolkit for computing protein-protein interface metrics.

## Overview

Epitype provides both a Python API and a CLI for analyzing binding interfaces between protein chains. It computes metrics including:

- **Buried Surface Area (dSASA)**: Surface area buried upon complex formation
- **Shape Complementarity (Sc)**: How well the two surfaces fit together
- **Hydrogen Bonds**: Cross-interface hydrogen bond count
- **PackStat**: Packing quality at the interface
- **Binding Energy**: Estimated binding energy (requires OpenMM)

### Feature Availability by Platform

| Feature | Dependencies | Linux | macOS | Windows |
|---------|-------------|:-----:|:-----:|:-------:|
| SASA, H-bonds, PackStat | Core | ✅ | ✅ | ✅ |
| Energy calculations | Core | ✅ | ✅ | ✅ |
| Shape complementarity | NanoShaper | ✅ | ❌ | ❌ |

## Installation

```bash
pip install epitype
```

For development:

```bash
pip install -e ".[dev]"
```

### Optional Dependencies

All core features (SASA, H-bonds, PackStat) work out-of-the-box with the standard `pip` installation method. However, energy calculations and shape complementarity require additional dependencies that cannot be installed with pip or may need to be separately configured to work optimally with your system.

#### Energy Calculations

Energy calculations use **OpenMM** and **pdbfixer**, which are included in the core dependencies and install automatically. For CUDA-accelerated calculations or specific CUDA version matching, see the [OpenMM User Guide](https://docs.openmm.org/latest/userguide/index.html).

#### Shape Complementarity (NanoShaper)

Shape complementarity calculations require **NanoShaper**, which must be built from source.

> [!NOTE]
> Shape complementarity is only available on Linux x86-64 platforms, as NanoShaper does not support macOS or Windows.

**Prerequisites (Ubuntu/Debian):**
```bash
sudo apt-get install libboost-all-dev libgmp-dev libmpfr-dev
```

**Prerequisites (Fedora/RHEL):**
```bash
sudo dnf install boost-devel gmp-devel mpfr-devel
```

**Build NanoShaper:**
```bash
git clone https://gitlab.iit.it/SDecherchi/nanoshaper.git
cd nanoshaper
python setup.py
# Follow prompts to build standalone executable
# Add the resulting binary to your PATH
```

**Verify installation:**
```bash
NanoShaper  # Should print usage information
```

For more information, see the [NanoShaper documentation](https://gitlab.iit.it/SDecherchi/nanoshaper).

## Quick Start

### CLI

```bash
# Full interface analysis
epitype run complex.pdb --chains HL_A

# Save results to JSON or CSV
epitype run complex.pdb --chains HL_A -o results.json
epitype run complex.pdb --chains HL_A -o results.csv

# Fast analysis (skip energy and shape complementarity)
epitype run complex.pdb --chains HL_A --no-energy --no-shape

# Individual metrics
epitype sasa complex.pdb --chains HL_A
epitype hbonds complex.pdb --chains HL_A
epitype packstat complex.pdb --chains HL_A

# Show version
epitype version
```

### Python API

```python
from epitype import analyze, load_structure, detect_interface

# Full analysis pipeline
metrics = analyze("complex.pdb", chains="HL_A")
print(f"Buried surface area: {metrics.dSASA_int:.1f} A^2")
print(f"Interface H-bonds: {metrics.hbonds_int}")
print(f"PackStat score: {metrics.packstat:.3f}")

# Working with Structure objects
structure = load_structure("complex.pdb")
interface = detect_interface(structure, ["H", "L"], ["A"])
print(f"Interface residues: {interface.num_interface_residues}")

# Individual metrics
from epitype.metrics import sasa, hbonds, packstat

sasa_result = sasa.calculate_sasa(structure)
print(f"Total SASA: {sasa_result.total:.1f} A^2")

hbond_result = hbonds.analyze_hbonds(structure, interface)
print(f"Cross-interface H-bonds: {hbond_result.interface_hbonds}")
```

### Configuration Options

```python
from epitype import analyze, AnalysisConfig

# Custom configuration
config = AnalysisConfig(
    chain_spec="AB_C",
    interface_cutoff=10.0,      # Distance cutoff for interface detection
    probe_radius=1.4,           # SASA probe radius
    compute_energy=False,       # Skip energy (requires OpenMM)
    compute_shape=False,        # Skip shape complementarity (requires NanoShaper)
    compute_packstat=True,      # Include PackStat
)

metrics = analyze("complex.pdb", chains="AB_C", **config.__dict__)
```

### Output Format

```python
# Convert to dictionary
result_dict = metrics.to_dict()

# Available fields
{
    "dG_separated": 0.0,           # Binding energy (REU)
    "dSASA_int": 1850.5,           # Buried surface area (A^2)
    "dG_per_dSASA": 0.0,           # Energy per buried surface
    "sc_value": 0.68,              # Shape complementarity (0-1)
    "packstat": 0.75,              # Packing quality (0-1)
    "hbonds_int": 12,              # Interface H-bonds
    "delta_unsat_hbonds": 2,       # Buried unsatisfied H-bond donors/acceptors
    "interface_residues_1": [...], # Group 1 interface residue IDs
    "interface_residues_2": [...], # Group 2 interface residue IDs
    "sasa_complex": 45000.0,       # Complex SASA (A^2)
    "sasa_separated": 46850.5,     # Separated SASA (A^2)
}
```

## Metrics Reference

| Metric | Description | Range/Units | Interpretation |
|--------|-------------|-------------|----------------|
| `dG_separated` | Binding energy | REU | Lower (more negative) = stronger binding |
| `dSASA_int` | Buried surface area | Angstrom^2 | Larger = more extensive interface |
| `sc_value` | Shape complementarity | 0-1 | > 0.65 indicates good fit |
| `packstat` | Packing quality | 0-1 | > 0.65 indicates tight packing |
| `hbonds_int` | Cross-interface H-bonds | count | More = better electrostatic complementarity |
| `delta_unsat_hbonds` | Unsatisfied buried H-bond partners | count | Lower = better (fewer buried polar groups without partners) |

## Chain Specification

Chains are specified using underscore-separated groups:

- `A_B` - Chain A vs Chain B
- `HL_A` - Chains H+L vs Chain A (antibody heavy/light vs antigen)
- `AB_CD` - Chains A+B vs Chains C+D

## API Reference

### Top-level Functions

```python
# Load a structure
structure = load_structure("file.pdb")  # or .cif

# Detect interface between chain groups
interface = detect_interface(structure, group1=["H", "L"], group2=["A"], cutoff=8.0)

# Parse chain specification
group1, group2 = parse_chain_groups("HL_A")  # Returns (["H", "L"], ["A"])

# Full analysis
metrics = analyze("file.pdb", chains="HL_A")

# Analysis with pre-loaded structure
metrics = analyze_structure(structure, chain_spec="HL_A")
```

### Metrics Modules

```python
from epitype.metrics import sasa, hbonds, packstat, shape, energy

# SASA calculations
sasa_result = sasa.calculate_sasa(structure)
dsasa = sasa.calculate_dsasa(complex_struct, separated_struct, interface)

# Hydrogen bonds
hbond_result = hbonds.analyze_hbonds(structure, interface)

# PackStat (packing quality)
packstat_result = packstat.calculate_packstat(structure, interface)

# Shape complementarity (requires NanoShaper)
sc_result = shape.calculate_shape_complementarity(structure, interface)

# Binding energy (requires OpenMM)
dG = energy.calculate_binding_energy(complex_struct, separated_struct, interface)
```

## Development

### Running Tests

```bash
# All tests
pytest

# With coverage
pytest --cov=epitype --cov-report=term-missing

# Specific test file
pytest tests/unit/test_sasa.py -v
```

### Type Checking

```bash
mypy src/epitype
```

### Linting

```bash
ruff check src/epitype tests
```

## License

MIT License - see LICENSE file for details.

## Citation

If you use this software, please cite:

```
Epitype: A Python toolkit for protein-protein interface analysis
https://github.com/briney/epitype
```

## Acknowledgments

The algorithms are based on:

- Lawrence & Colman (1993) for shape complementarity
- Shrake & Rupley (1973) for SASA calculation
