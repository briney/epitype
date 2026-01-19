# Epitype: Implementation Task List

> **Reference**: See `WORKPLAN.md` for detailed implementation specifications and code.

---

## Progress Summary

| Phase | Description | Tasks | Completed |
|-------|-------------|-------|-----------|
| 1 | Project Setup & Core Infrastructure | 15 | 15 |
| 2 | Interface Detection & Separation | 12 | 12 |
| 3 | Core Metrics (+ PackStat Optimization) | 21 | 21 |
| 4 | Advanced Metrics | 12 | 12 |
| 5 | Pipeline & CLI | 18 | 18 |
| 6 | API, Tests & Documentation | 16 | 16 |
| **Total** | | **94** | **94** |

---

## Phase 1: Project Setup & Core Infrastructure

### 1.1 Directory Structure
> Reference: WORKPLAN.md Section 1

- [x] Create root project directory structure
- [x] Create `src/epitype/` package directory
- [x] Create `src/epitype/core/` subdirectory
- [x] Create `src/epitype/io/` subdirectory
- [x] Create `src/epitype/metrics/` subdirectory
- [x] Create `src/epitype/surface/` subdirectory
- [x] Create `src/epitype/analysis/` subdirectory
- [x] Create `src/epitype/cli/` subdirectory
- [x] Create `src/epitype/utils/` subdirectory
- [x] Create `tests/` directory structure (unit/, integration/, e2e/, data/)

### 1.2 Configuration Files
> Reference: WORKPLAN.md Section 2

- [x] Create `pyproject.toml` with dependencies and build configuration
- [x] Create `py.typed` marker file for PEP 561
- [x] Create `.gitignore` file
- [x] Create `LICENSE` file (MIT)
- [x] Create placeholder `README.md`

### 1.3 Core Types
> Reference: WORKPLAN.md Section 3.1

- [x] Implement `src/epitype/core/types.py`
  - [x] Define `Coords`, `AtomIndex`, `ResidueIndex`, `ChainID` type aliases
  - [x] Define physical constants (`PROBE_RADIUS`, `DEFAULT_INTERFACE_CUTOFF`, `SEPARATION_DISTANCE`)
  - [x] Implement `AtomType` enum
  - [x] Implement `VdWRadii` dataclass
  - [x] Implement `InterfaceMetrics` dataclass with `to_dict()` method

---

## Phase 2: Interface Detection & Separation

### 2.1 Structure Representation
> Reference: WORKPLAN.md Section 3.2

- [x] Implement `src/epitype/core/structure.py`
  - [x] Implement `Atom` dataclass
    - [x] Properties: `vdw_radius`, `is_backbone`, `is_hydrogen`
    - [x] Method: `distance_to()`
  - [x] Implement `Residue` dataclass
    - [x] Properties: `ca`, `cb`, `center_of_mass`, `full_id`
    - [x] Method: `get_atom()`
  - [x] Implement `Chain` dataclass
    - [x] Properties: `atoms` (iterator), `coords`
  - [x] Implement `Structure` dataclass
    - [x] Properties: `atoms`, `residues`, `coords`, `num_atoms`, `num_residues`
    - [x] Methods: `get_chain()`, `get_chains()`, `subset()`, `copy()`, `translate()`

- [x] Write `tests/unit/test_structure.py`
  - [x] Test `Atom` creation and properties
  - [x] Test `Residue` CA/CB access and glycine fallback
  - [x] Test `Chain` iteration and coordinate extraction
  - [x] Test `Structure` subsetting, copying, and translation

### 2.2 Interface Detection
> Reference: WORKPLAN.md Section 3.3

- [x] Implement `src/epitype/core/interface.py`
  - [x] Implement `InterfaceRegion` dataclass
  - [x] Implement `detect_interface()` function using KD-tree
  - [x] Implement `parse_chain_groups()` function

- [x] Write `tests/unit/test_interface.py`
  - [x] Test chain group parsing (simple, multi-chain, invalid)
  - [x] Test interface detection with various cutoffs
  - [x] Test contact pair detection

### 2.3 Chain Separation
> Reference: WORKPLAN.md Section 3.4

- [x] Implement `src/epitype/core/separation.py`
  - [x] Implement `compute_interface_normal()` function
  - [x] Implement `separate_chains()` function
  - [x] Implement `compute_separation_distance()` function

- [x] Write `tests/unit/test_separation.py`
  - [x] Test interface normal computation
  - [x] Test separation increases distance
  - [x] Test original structure unchanged after separation

### 2.4 Geometry Utilities
> Reference: WORKPLAN.md Section 14.1

- [x] Implement `src/epitype/utils/geometry.py`
  - [x] Implement `normalize()` function
  - [x] Implement `angle_between()` function
  - [x] Implement `dihedral_angle()` function
  - [x] Implement `rotation_matrix()` function
  - [x] Implement `centroid()` function
  - [x] Implement `rmsd()` function

- [x] Write `tests/unit/test_geometry.py`
  - [x] Test vector normalization
  - [x] Test angle calculations (parallel, perpendicular, antiparallel)

---

## Phase 3: Core Metrics

### 3.1 I/O Parsers
> Reference: WORKPLAN.md Section 4.1

- [x] Implement `src/epitype/io/parsers.py`
  - [x] Implement `parse_pdb()` function using BioPython
  - [x] Implement `parse_mmcif()` function
  - [x] Implement `parse_structure()` auto-detect function
  - [x] Implement `_convert_biopython_structure()` helper

- [x] Write `tests/unit/test_parsers.py`
  - [x] Test PDB parsing
  - [x] Test chain ID preservation
  - [x] Test auto-format detection

### 3.2 Test Data Setup
> Reference: WORKPLAN.md Section 9.5

- [x] Create `tests/data/README.md`
- [x] Download `1YY9.pdb` test structure
- [x] Download `4FQI.pdb` test structure
- [x] Create `tests/conftest.py` with fixtures
  - [x] `test_data_dir` fixture
  - [x] `pdb_1yy9` and `pdb_4fqi` fixtures
  - [x] `temp_dir` fixture
  - [x] `minimal_structure` fixture
  - [x] `interface_ab` fixture

### 3.3 SASA Calculation
> Reference: WORKPLAN.md Section 5.1

- [x] Implement `src/epitype/metrics/sasa.py`
  - [x] Implement `SASAResult` dataclass
  - [x] Implement `calculate_sasa()` using FreeSASA
  - [x] Implement `calculate_dsasa()` function
  - [x] Implement `calculate_interface_sasa()` function

- [x] Write `tests/unit/test_sasa.py`
  - [x] Test SASA is positive
  - [x] Test polar + apolar = total
  - [x] Test per-atom and per-residue SASA
  - [x] Test dSASA calculation

### 3.4 Hydrogen Bond Detection
> Reference: WORKPLAN.md Section 5.3

- [x] Implement `src/epitype/metrics/hbonds.py`
  - [x] Define `DONOR_ATOMS` and `ACCEPTOR_ATOMS` sets
  - [x] Implement `HydrogenBond` dataclass
  - [x] Implement `HBondResult` dataclass
  - [x] Implement `detect_hydrogen_bonds()` function
  - [x] Implement `count_interface_hbonds()` function
  - [x] Implement `count_unsatisfied_hbonds()` function
  - [x] Implement `analyze_hbonds()` function

- [x] Write `tests/unit/test_hbonds.py`
  - [x] Test H-bond detection
  - [x] Test distance cutoff enforcement
  - [x] Test cross-interface counting
  - [x] Test no same-residue H-bonds

### 3.5 PackStat
> Reference: WORKPLAN.md Section 5.5

- [x] Implement `src/epitype/metrics/packstat.py`
  - [x] Define `PROBE_RADII` array (30 values, 0.1-3.0 Å)
  - [x] Implement `PackStatResult` dataclass
  - [x] Implement `calculate_packstat()` function using FreeSASA
  - [x] Implement `_freesasa_total()` helper for fast SASA calculation
  - [x] Implement `_reference_sasa()` helper
  - [x] Implement `_interface_packstat()` helper
  - [x] Implement `_estimate_cavity_volume()` helper

- [x] Write `tests/unit/test_packstat.py`
  - [x] Test packstat returns 0-1 range
  - [x] Test per-probe scores
  - [x] Test interface-specific packstat
  - [x] Test performance (< 30 seconds for real structure)

### 3.6 PackStat Optimization
> Optimized PackStat implementation to use FreeSASA for ~100x speedup

- [x] Replace custom `_shrake_rupley_sasa()` with FreeSASA C library
- [x] Remove `_fibonacci_sphere()` helper (no longer needed)
- [x] Remove `n_sphere_points` parameter from function signatures
- [x] Update tests to remove `@pytest.mark.slow` markers
- [x] Update WORKPLAN.md Section 5.5 with new implementation
- [x] Verify all tests pass with improved performance (~4s vs 12+ min)

---

## Phase 4: Advanced Metrics (External Dependencies)

### 4.1 MSMS Surface Wrapper
> Reference: WORKPLAN.md Section 5.6

- [x] Implement `src/epitype/surface/types.py`
  - [x] Define `SurfacePoint` dataclass

- [x] Implement `src/epitype/surface/msms.py`
  - [x] Implement `generate_surface()` function
  - [x] Implement `_write_xyzr()` helper
  - [x] Implement `_parse_msms_output()` helper
  - [x] Implement `check_msms_available()` function

- [x] Write `tests/unit/test_msms.py` (skip if MSMS unavailable)
  - [x] Test MSMS availability check
  - [x] Test surface generation

### 4.2 Shape Complementarity
> Reference: WORKPLAN.md Section 5.4

- [x] Implement `src/epitype/metrics/shape.py`
  - [x] Implement `ShapeComplementarityResult` dataclass
  - [x] Implement `calculate_shape_complementarity()` function
  - [x] Implement `_trim_surface_to_interface()` helper
  - [x] Implement `_calculate_sc_one_direction()` helper

- [x] Write `tests/unit/test_shape.py`
  - [x] Test Sc returns 0-1 range
  - [x] Test Sc calculation with known geometry

### 4.3 Binding Energy (OpenMM)
> Reference: WORKPLAN.md Section 5.2

- [x] Implement `src/epitype/metrics/energy.py`
  - [x] Define `KJ_TO_REU` conversion constant
  - [x] Implement `prepare_structure_for_openmm()` function
  - [x] Implement `calculate_potential_energy()` function
  - [x] Implement `calculate_binding_energy()` function
  - [x] Implement `_write_temp_pdb()` helper

- [x] Write `tests/unit/test_energy.py` (skip if OpenMM unavailable)
  - [x] Test energy calculation returns float
  - [x] Test binding energy is negative for favorable binding

---

## Phase 5: Pipeline & CLI

### 5.1 Output Writers
> Reference: WORKPLAN.md Section 4.2

- [x] Implement `src/epitype/io/writers.py`
  - [x] Implement `write_json()` function
  - [x] Implement `write_csv()` function
  - [x] Implement `format_metrics_table()` function

- [x] Write `tests/unit/test_writers.py`
  - [x] Test JSON output format
  - [x] Test CSV output format
  - [x] Test table formatting

### 5.2 Analysis Pipeline
> Reference: WORKPLAN.md Section 6.1

- [x] Implement `src/epitype/analysis/pipeline.py`
  - [x] Implement `AnalysisConfig` dataclass
  - [x] Define `ProgressCallback` type alias
  - [x] Implement `analyze_interface()` function
  - [x] Implement `analyze_structure()` convenience function

- [x] Write `tests/integration/test_pipeline.py`
  - [x] Test full analysis pipeline
  - [x] Test analysis with selective metrics
  - [x] Test invalid chain spec error handling
  - [x] Test progress callback

### 5.3 CLI Main Entry Point
> Reference: WORKPLAN.md Section 7.1

- [x] Implement `src/epitype/cli/main.py`
  - [x] Create Typer app
  - [x] Add subcommand typers
  - [x] Implement `version` command
  - [x] Implement main callback

### 5.4 CLI Run Command
> Reference: WORKPLAN.md Section 7.2

- [x] Implement `src/epitype/cli/run.py`
  - [x] Implement `run` command with all options
  - [x] Handle JSON/CSV output
  - [x] Handle progress display

### 5.5 CLI Subcommands
> Reference: WORKPLAN.md Sections 7.3, 13.1-13.4

- [x] Implement `src/epitype/cli/sasa.py`
- [x] Implement `src/epitype/cli/energy.py`
- [x] Implement `src/epitype/cli/hbonds.py`
- [x] Implement `src/epitype/cli/shape.py`
- [x] Implement `src/epitype/cli/packstat.py`

### 5.6 CLI Common Utilities
> Reference: WORKPLAN.md Section 7.4

- [x] Implement `src/epitype/cli/common.py`
  - [x] Define common option types
  - [x] Implement `validate_chain_spec()` function

### 5.7 CLI Tests
> Reference: WORKPLAN.md Section 9.3.2

- [x] Write `tests/integration/test_cli.py`
  - [x] Test `run` command basic execution
  - [x] Test JSON output
  - [x] Test CSV output
  - [x] Test `sasa` command
  - [x] Test `version` command

---

## Phase 6: API, Tests & Documentation

### 6.1 Package Exports
> Reference: WORKPLAN.md Sections 8.1, 8.2

- [x] Implement `src/epitype/__init__.py`
  - [x] Define `__version__`
  - [x] Export core data structures
  - [x] Export high-level API functions
  - [x] Define `__all__`

- [x] Implement `src/epitype/metrics/__init__.py`
  - [x] Export metric submodules

- [x] Create all other `__init__.py` files
  - [x] `src/epitype/core/__init__.py`
  - [x] `src/epitype/io/__init__.py`
  - [x] `src/epitype/surface/__init__.py`
  - [x] `src/epitype/analysis/__init__.py`
  - [x] `src/epitype/cli/__init__.py`
  - [x] `src/epitype/utils/__init__.py`

### 6.2 API Integration Tests
> Reference: WORKPLAN.md Section 9.3.3

- [x] Write `tests/integration/test_api.py`
  - [x] Test version accessibility
  - [x] Test `analyze()` function
  - [x] Test `load_structure()` function
  - [x] Test `detect_interface()` function
  - [x] Test metrics submodule access

### 6.3 End-to-End Tests
> Reference: WORKPLAN.md Section 9.4

- [x] Write `tests/e2e/test_full_analysis.py`
  - [x] Test complete workflow (file → analysis → output)
  - [x] Test known structure metric ranges

- [x] Create expected values files
  - [x] `tests/data/expected/1YY9_metrics.json`
  - [x] `tests/data/expected/4FQI_metrics.json`

### 6.4 Documentation
> Reference: WORKPLAN.md Section 1 (docs/)

- [x] Update `README.md` with full documentation
  - [x] Installation instructions
  - [x] Quick start guide
  - [x] CLI usage examples
  - [x] Python API examples

- [x] Create `CHANGELOG.md`

- [ ] Create `docs/` directory (optional, for MkDocs)
  - [ ] `docs/index.md`
  - [ ] `docs/installation.md`
  - [ ] `docs/cli.md`
  - [ ] `docs/api.md`
  - [ ] `docs/metrics.md`

### 6.5 Quality Assurance

- [x] Run full test suite: `pytest --cov=epitype`
- [x] Run type checker: `mypy src/epitype`
- [x] Run linter: `ruff check src/epitype tests`
- [ ] Verify CLI installation: `pip install -e . && epitype version`

### 6.6 Validation Against Known Structures

- [x] Validate 1YY9 metrics
  - [x] dSASA in range 2500-3500 Ų (actual ~3100)
  - [x] H-bonds in range 10-25 (actual 16)
  - [ ] Sc > 0.6 (if MSMS available)

- [x] Validate 4FQI metrics
  - [x] dSASA in range 4500-6500 Ų (actual ~5689)
  - [x] H-bonds in range 2-15 (actual 4)

---

## External Dependencies Checklist

### Required for Core Functionality
- [ ] Python 3.10+ installed
- [ ] Install core dependencies: `pip install numpy scipy biopython freesasa typer rich pydantic`

### Required for Energy Calculations
- [ ] Install OpenMM: `conda install -c conda-forge openmm` or `pip install openmm`
- [ ] Install PDBFixer: `pip install pdbfixer`
- [ ] Verify: `python -c "import openmm; print(openmm.__version__)"`

### Required for Shape Complementarity
- [ ] Download MSMS from https://ccsb.scripps.edu/msms/downloads/
- [ ] Install MSMS to PATH
- [ ] Verify: `msms -h`

### Development Dependencies
- [ ] Install dev dependencies: `pip install -e ".[dev]"`

---

## Notes

### Marking Tasks Complete

Update task status by changing `- [ ]` to `- [x]`:
```markdown
- [x] Completed task
- [ ] Pending task
```

### Updating Progress Summary

After completing tasks, update the "Completed" column in the Progress Summary table.

### Skipping Optional Tasks

Some tasks may be skipped if dependencies are unavailable:
- Shape complementarity tests (requires MSMS)
- Energy calculation tests (requires OpenMM)
- Documentation tasks (can be done post-MVP)

### Task Dependencies

- Phase 2 depends on Phase 1 completion
- Phase 3 depends on Phase 2 completion (structure parsing needed)
- Phase 4 depends on Phase 3 completion (core metrics first)
- Phase 5 depends on Phases 3-4 completion (all metrics needed for pipeline)
- Phase 6 depends on Phase 5 completion (full functionality needed for tests)

### Running Tests Incrementally

```bash
# After Phase 2
pytest tests/unit/test_structure.py tests/unit/test_interface.py -v

# After Phase 3
pytest tests/unit/ -v

# After Phase 5
pytest tests/unit/ tests/integration/ -v

# Full suite
pytest --cov=epitype --cov-report=term-missing
```
