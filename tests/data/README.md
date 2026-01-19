# Test Data

This directory contains PDB structures used for testing the interface analyzer.

## Structures

### 1YY9.pdb
- **Description**: Structure of the extracellular domain of the epidermal growth factor receptor (EGFR) in complex with the Fab fragment of Cetuximab (IMC-C225)
- **Resolution**: 2.8 Å
- **Chains**:
  - H, L: Antibody heavy and light chains (Fab)
  - A: EGFR extracellular domain (antigen)
- **Interface**: H,L vs A (antibody-antigen interface)
- **Expected metrics**:
  - dSASA_int: ~1800-2200 Å²
  - H-bonds: ~8-20
  - Shape complementarity: >0.60

### 4FQI.pdb
- **Description**: Antibody-antigen complex structure
- **Chains**: Multiple chains forming antibody-antigen interface
- **Interface**: Antibody-antigen binding interface
- **Expected metrics**:
  - dSASA_int: ~1500-2000 Å²
  - H-bonds: ~5-15

## Usage

These structures are automatically downloaded from the RCSB PDB database during test setup.
They are used by the pytest fixtures defined in `tests/conftest.py`.

## Sources

All structures are obtained from the RCSB Protein Data Bank:
- https://www.rcsb.org/structure/1YY9
- https://www.rcsb.org/structure/4FQI
