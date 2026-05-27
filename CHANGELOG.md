# Changelog

## [0.2.0] — 2026-05-27

### Added

- `tests/test_format_equivalence.py`: cross-loader equivalence on 1FYT —
  asserts `from_pdb(.pdb)`, `from_pdb(.cif)`, `from_structure(parsed)`,
  `from_biotite(biotite CIF parse)`, and `from_boltzgen_structure(numpy
  layout)` all produce the same SC and atom selection.
- `_validate_boltzgen_layout()` runtime check: `from_boltzgen_structure`
  now raises a clear `TypeError` if the passed object is missing one of
  the required `.chains` / `.residues` / `.atoms` attributes or a required
  field within them.
- `validation/rosetta_interface_sc.py` + `tests/test_rosetta_validation.py`:
  optional comparison against Rosetta `InterfaceAnalyzerMover` (skipped
  unless PyRosetta is installed).
- `tests/test_batch.py`: dedicated unit tests for `score_many` covering
  both default (raise) and `on_error="record"` modes.
- "Intake-format equivalence" section in README documenting the cross-
  loader test matrix.
- "A note on the SC range" section explaining the math and conventions.

### Changed

- **Breaking:** `score_many` raises by default on any failed file; use
  `on_error="record"` for the previous per-row error-reporting behavior.
- **Breaking:** PDB/mmCIF parsing now fails loudly for non-finite
  coordinates, missing chains, overlapping chain groups, ambiguous
  alternate locations, and empty post-filter atom groups.
- **Breaking:** alternate-location handling is now explicit via
  `altloc_policy` (`"fail"` / `"highest_occupancy"` / `"sc_rs"`); default
  is `"fail"` for safety. `"sc_rs"` matches the original CLI behavior.
- `from_boltzgen_refold` docstring now explains the `label_asym_id` vs
  `auth_asym_id` distinction and why this wrapper exists.
- `from_boltzgen_structure` is **no longer exported** from the public API
  pending end-to-end validation against a Structure produced by an actual
  BoltzGen inference run. The implementation remains in
  `shape_complementarity.io` for power users and tests.

### Fixed

- `score_many` no longer silently swallows errors by default.

## [0.1.0] — 2026-05-25

### Added

- `compute_sc()`: low-level PyO3 binding to `sc-rs`'s `ScCalculator`, accepting
  pre-parsed atom coordinate arrays.
- `ScResult`: read-only Python class exposing `sc`, `median_distance`,
  `trimmed_area`, `atoms_a`, `atoms_b`.
- `from_pdb()` / `from_structure()`: biopython-based PDB/mmCIF parsing that
  mirrors the filtering logic of the `sc-rs` CLI exactly (ATOM-only, altloc A/'
  ' , heavy atoms only by default).
- `score_many()`: multiprocessing batch scorer returning a `pd.DataFrame`;
  per-file exceptions are caught and reported in a `status`/`error` column.
- Full test suite: unit tests for the Rust core, IO layer, and parity tests
  comparing `shape_complementarity` output against the `sc-rs` CLI binary.
- Benchmark script (`benchmark/speed.py`).
- Pinned to `sc-rs` v1.0.0 via git tag.

### Changed

- Scientific scoring now fails loudly for non-finite coordinates, missing
  chains, overlapping chain groups, ambiguous alternate locations, and empty
  post-filter atom groups.
- `score_many()` raises by default on failed files; use `on_error="record"` for
  the previous per-row error reporting behavior.
- PDB/mmCIF alternate-location handling is explicit via `altloc_policy`.

### Validation

- Added a `.venv`-based development setup path and optional PyRosetta
  `InterfaceAnalyzerMover` validation script.
