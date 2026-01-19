# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-01-18

### Added

- Initial release of Epitype
- Core structure representation (`Structure`, `Chain`, `Residue`, `Atom`)
- Interface detection using KD-tree spatial indexing
- Chain separation for binding energy calculations
- SASA (Solvent Accessible Surface Area) calculation using FreeSASA
  - Per-atom and per-residue SASA
  - Polar/apolar decomposition
  - dSASA (buried surface area) calculation
- Hydrogen bond detection and analysis
  - Cross-interface H-bond counting
  - Unsatisfied buried H-bond detection
- PackStat packing quality metric
  - Multi-probe SASA calculation
  - Interface-specific packing scores
- Shape complementarity calculation (requires MSMS)
  - Lawrence & Colman (1993) algorithm implementation
- Binding energy calculation (requires OpenMM)
  - AMBER ff14SB force field
  - Optional energy minimization
- Analysis pipeline for complete interface analysis
- Python API with high-level convenience functions
- Comprehensive test suite
  - Unit tests for all modules
  - Integration tests for API
  - End-to-end tests for workflows

### Dependencies

- Required: numpy, scipy, biopython, freesasa
- Optional: openmm, pdbfixer (for energy calculations)
- Optional: MSMS binary (for shape complementarity)

### Notes

- CLI implementation planned for future release
- Documentation available in README.md
