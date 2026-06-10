# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

Fault Diagnosis Toolbox (`faultdiagnosistoolbox`) is a Python package for analysis and design of model-based fault
diagnosis systems for dynamic systems described by differential-algebraic equations. Core capabilities: structural
analysis of large-scale models, fault isolability/diagnosability analysis, sensor placement, and C/C++/Python code
generation for residual generators. It is a Python port of a Matlab toolbox.

## Build / install

The package is a mixed Python + C++ project built with `scikit-build-core` + `nanobind` + CMake. The compiled piece
is the `dmpermlib` extension module (Dulmage-Mendelsohn permutation / structural analysis library), built from
`src/src_extension/*.cc` plus the bundled `CSparse` library in `src/extern/CSparse`.

```bash
# Editable/dev install (rebuilds the C++ extension)
pip install -e .

# Build a wheel
python -m build --wheel        # or: make build
```

`CMakeLists.txt` at the repo root drives the CMake build of `dmpermlib`, linking against `csparse` and requiring
`numpy` headers. Requires Python >= 3.10. Prebuilt `.so` files for several Python versions/platforms are checked
into `src/faultdiagnosistoolbox/` (e.g. `dmpermlib.cpython-313-darwin.so`).

## Tests

```bash
make test            # python -m pytest -v
make test-cov        # with coverage for faultdiagnosistoolbox
make test-cov-report # coverage with HTML report

# Run a single test file
python -m pytest -v tests/watertank_test.py

# Run a single test function
python -m pytest -v tests/watertank_test.py::test_name
```

Tests live in `tests/` and exercise the toolbox against example models defined alongside the tests
(`watertank_model.py`, `electricmotor_model.py`).

## Documentation

```bash
make doc_build   # sphinx-build -M html docs doc_build
```

## Architecture

The Python package lives in `src/faultdiagnosistoolbox/`:

- **`DiagnosisModel.py`** - the central `DiagnosisModel` class. Models are constructed from a dict (`modelDef`) of
  type `"MatrixStruc"` (pure structure: `X`/`F`/`Z` incidence matrices) or `"Symbolic"` (sympy expressions over
  variables `x` (unknowns), `f` (faults), `z` (known signals), plus `parameters` and `rels`). Differential
  constraints between variables are declared with `DiffConstraint(dvar, ivar)`. This module provides diagnosability
  analysis, plotting hooks, and orchestrates code generation and sensor placement.
- **`dmperm.py`** - structural analysis built on Dulmage-Mendelsohn decomposition (`GetDMParts`, `PSODecomposition`,
  `MSO`, `Mplus`, observability/index checks, MTES computation). Wraps the compiled `dmpermlib` extension
  (`dmpermlib.cpython-*.so`) for the underlying sparse-matrix permutation algorithms implemented in
  `src/src_extension/` (`StructuralAnalysisModel`, `MSOAlg`, `SparseMatrix`).
- **`Matching.py`** - bipartite matching of equations to variables (`Matching`, `HallComponentMatching`,
  `MixedCausalityMatching`) used to derive computational causality for residual generators.
- **`CodeGeneration.py`** - generates executable residual generator code (Python and C) from a matched/causalized
  model: builds Hall components (algebraic, integral, mixed causality), writes per-equation code, and assembles full
  residual generator modules/files including setup scripts for compiling generated C extensions.
- **`SensorPlacement.py`** - algorithms for determining minimal sensor sets to achieve required
  detectability/isolability specifications.
- **`TestSelection.py` / `testselection_rf.py`** - test (residual) selection, including a random-forest-based
  approach.
- **`MHS.py`** - minimal hitting set computation, used in isolability/diagnosability analysis.
- **`StructurePlotting.py`** - plotting of model structure, Dulmage-Mendelsohn decomposition, and isolability
  matrices (matplotlib).
- **`diag_util.py`** - diagnosis utilities: isolability matrices, confusion matrices, plotting helpers.
- **`models/`** - bundled example models.

`src/Engine/` and `src/lab/` and `src/code_examples/` contain example/experimental scripts (e.g. residual generator
code-generation outputs `ResGen_*.cc`/`.so`, an electric-motor case study) — useful as reference for what generated
code looks like, but not part of the installed package.

## Typical workflow in this codebase

1. Define a model (`DiagnosisModel`) either structurally (incidence matrices `X`/`F`/`Z`) or symbolically (sympy
   relations + `DiffConstraint`).
2. Run structural analysis (`dmperm`) to find overdetermined parts / MSOs, check observability/index.
3. Select tests/residuals (`TestSelection`, `MHS`, `SensorPlacement`).
4. Generate residual generator code (`CodeGeneration`) for the chosen MSOs/causality.
5. Evaluate diagnosis performance via `diag_util` (isolability/confusion matrices).
