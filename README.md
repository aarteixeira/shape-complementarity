# shape_complementarity [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Python bindings to [sc-rs](https://github.com/cytokineking/sc-rs) for computing
**Lawrence-Colman Shape Complementarity (SC)** between two protein chain groups.
Useful for analyzing protein–protein interfaces, including in computational
design pipelines.

`shape_complementarity` is **not** a reimplementation of the SC algorithm. All algorithmic
correctness comes from the upstream Rust crate `sc-rs`, which is invoked
directly via PyO3. If you find a numerical result surprising, verify against the
`sc-rs` CLI first.

---

## Install

**Development:**

Use any isolated Python 3.10+ environment. For example:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -U pip setuptools wheel
python -m pip install -e ".[test]"
python -m pip install maturin
maturin develop --release
```

Re-run `maturin develop --release` after changing Rust code. If `maturin`
installs into the wrong Python environment, check that your active shell is not
still pointing at another environment, for example through `CONDA_PREFIX`.

**Published release, when available:**
```bash
pip install shape-complementarity
```

---

## Usage

### 1. From a PDB or mmCIF file path

The simplest entry point. Accepts `.pdb`, `.ent`, `.cif`, and `.mmcif` files.

```python
import shape_complementarity

result = shape_complementarity.from_pdb("complex.pdb", chains_a=["H"], chains_b=["A"])
print(result.sc)               # e.g. 0.714
print(result.median_distance)  # Å
print(result.trimmed_area)     # Å²
print(result.atoms_a, result.atoms_b)

# chains_b=None → all chains not in chains_a
result = shape_complementarity.from_pdb("complex.pdb", chains_a=["H"])

# mmCIF works identically
result = shape_complementarity.from_pdb("complex.cif", chains_a=["H"], chains_b=["A"])
```

By default, PDB/mmCIF parsing is strict: requested chains must exist,
`chains_a` and `chains_b` must not overlap, alternate locations raise an error,
and empty post-filter atom groups raise before the Rust core is called. Use
`altloc_policy="highest_occupancy"` to select the highest-occupancy conformer,
or `altloc_policy="sc_rs"` only when matching the upstream `sc-rs` CLI parser.

---

### 2. From a biopython Structure object

Use this when you have already parsed the file with biopython, or when you
need to manipulate the structure before scoring.

```python
from Bio.PDB import PDBParser
import shape_complementarity

parser = PDBParser(QUIET=True)
structure = parser.get_structure("complex", "complex.pdb")

result = shape_complementarity.from_structure(structure, chains_a=["H"], chains_b=["A"])
```

Works with any biopython `Structure` regardless of how it was created
(PDB, mmCIF, downloaded from RCSB, built programmatically, etc.).

---

### 3. From a biotite AtomArray

[biotite](https://www.biotite-python.org/) is used natively by BoltzGen's
analysis stack. Install `biotite` to use this API; it is included in the
development `test` extra. Pass an `AtomArray` or `AtomArrayStack` (first model
is used).

```python
import biotite.structure.io.pdbx as pdbx
import shape_complementarity

cif = pdbx.CIFFile.read("complex.cif")
atoms = pdbx.get_structure(cif, model=1, use_author_fields=False)

result = shape_complementarity.from_biotite(atoms, chains_a=["H"], chains_b=["A"])
```

---

### 4. From raw coordinate arrays

The lowest-level entry point — no file I/O, no parser overhead. Pass numpy
arrays or plain Python lists of `[x, y, z]` coordinates alongside atom and
residue names. Useful when coordinates are already in memory from a simulation
or generative model.

```python
import shape_complementarity

result = shape_complementarity.compute_sc(
    coords_a        = [[x, y, z], ...],   # shape (N, 3), Å
    atom_names_a    = ["CA", "CB", ...],
    residue_names_a = ["ALA", "ALA", ...],
    coords_b        = [[x, y, z], ...],
    atom_names_b    = ["CA", ...],
    residue_names_b = ["GLY", ...],
)
```

Atom radii are assigned automatically from the atom-name + residue-name pair.
If no specific residue/atom radius matches, `sc-rs` may use a generic element
fallback such as carbon or nitrogen. Atoms with no specific or generic radius
raise an error. Ligand or non-protein scoring therefore requires deliberate
radii validation before interpreting SC values.

---

### 5. From BoltzGen output

BoltzGen writes full-atom, Boltz-validated complexes to `refold_cif/*.cif`.
These are the right structures to score — post-generation files have zeroed
sidechain coordinates and should not be used for SC.

**Option A — file path (no extra dependencies):**
```python
import shape_complementarity

# Works exactly like from_pdb; biopython handles the mmCIF
result = shape_complementarity.from_pdb(
    "output/intermediate_designs_inverse_folded/refold_cif/design_0.cif",
    chains_a=["B"],  # binder
    chains_b=["A"],  # target
)
```

**Option B — via biotite (matches BoltzGen's own analysis stack):**
```python
import shape_complementarity

result = shape_complementarity.from_boltzgen_refold(
    "output/.../refold_cif/design_0.cif",
    chains_a=["B"],
    chains_b=["A"],
)
```

`from_boltzgen_refold()` requires `biotite`.

> **Which stage to use?** Only `refold_cif/` structures have complete, physically
> validated all-atom coordinates. `intermediate_designs/` NPZ files (post-generation)
> have backbone-only coordinates; `intermediate_designs_inverse_folded/` NPZ files
> have sidechains but have not yet been validated by Boltz.

---

### 6. Batch scoring

Score many files in parallel using `ProcessPoolExecutor`. By default, any
per-file failure raises a summary error after workers finish. Use
`on_error="record"` only for exploratory batches where failed rows should be
returned with `status="error"` and `NaN` numeric fields.

```python
from pathlib import Path
import shape_complementarity

paths = list(Path("refold_cif").glob("*.cif"))

df = shape_complementarity.score_many(
    paths,
    chains_a=["B"],
    chains_b=["A"],
    n_workers=8,
)

print(df[df.status == "ok"][["path", "sc"]].sort_values("sc", ascending=False))
```

Rayon parallelism is disabled inside each worker by default (`parallel=False`)
to avoid oversubscription with multiple processes.

---

## Intake-format equivalence

The public file/object loaders are locked together by
[`tests/test_format_equivalence.py`](tests/test_format_equivalence.py), which
verifies that loading the same structure (1FYT chains D vs A) through each
entry point yields the same `sc` and atom counts:

| Comparison | Tolerance | 1FYT result |
|---|---|---|
| `from_pdb(.pdb)` ↔ `from_pdb(.cif)` (biopython PDBParser vs MMCIFParser) | 10⁻⁶ | sc = 0.6597, 1521 + 1479 atoms |
| `from_structure(pre-parsed)` ↔ `from_pdb(file)` | 10⁻⁶ | identical |
| `from_biotite(CIF via biotite)` ↔ `from_pdb(.pdb)` | 10⁻³ SC, ±5 atoms | identical |

These comparisons are expected to pass at the tight end of their tolerance. A
future parser change that shifts atom selection or coordinate rounding should
fail this suite.

This is in addition to the value-matching tests in
[`tests/test_parity.py`](tests/test_parity.py), which run the actual `sc-rs`
CLI binary on disk and confirm `sc`, `median_distance`, `trimmed_area`, and
atom counts all agree to within 10⁻³ — a tolerance set by LLVM FMA-folding
differences between independently compiled cdylib and binary.

---

## ScResult

Single-structure scoring functions return an `ScResult` with these read-only
properties:

| Property | Type | Description |
|---|---|---|
| `sc` | `float` | Shape complementarity (see range note below; native interfaces typically 0.6–0.8) |
| `median_distance` | `float` | Median nearest-surface distance (Å) |
| `trimmed_area` | `float` | Total trimmed interface area (Å²) |
| `atoms_a` | `int` | Heavy atoms accepted for molecule A |
| `atoms_b` | `int` | Heavy atoms accepted for molecule B |

### A note on the SC range

The Lawrence–Colman SC statistic is the median of per-point scores of the form
`(n_A · n_B′) · exp(−w · d²)`, where `n_A` and `n_B′` are surface normals at a
point on A and its nearest neighbor on B (with one normal flipped so that two
perfectly mated surfaces score +1). The exponential lies in [0, 1] and the dot
product of unit normals lies in [−1, 1], so the **mathematical range of SC is
−1 to 1**.

In practice the original paper and downstream tools (CCP4 `sc`, Rosetta) quote
the range as **0 to 1**:

- **1.0** — perfect geometric complementarity (mated, parallel surfaces)
- **~0.6–0.8** — typical native protein–protein interface
- **~0** — uncorrelated / random surface orientations
- **< 0** — mathematically possible but not biologically meaningful

Negative values require both **co-located** points (`d ≈ 0`, so the
exponential is ≈ 1) **and co-aligned normals** (dot product ≈ +1, which after
the sign flip becomes −1). That is, two surfaces stacked on top of each other
rather than mated face-to-face. Concretely, scoring chain A of 1FYT against an
identical copy of itself placed at the same coordinates yields `sc ≈ −0.999`;
shifting one copy by just 1 Å already brings the score back to ≈ 0 because the
`exp(−w·d²)` term collapses. For ordinary protein interfaces, values should
usually fall in [0, 1].

`shape_complementarity` reports whatever `sc-rs` returns without clamping, so
the strict range is −1 to 1.

---

## Calibration recipe for design filters

```python
import numpy as np
import shape_complementarity

# Score known native complexes to establish a baseline
natives = shape_complementarity.score_many(native_pdbs, chains_a=["H"], chains_b=["A"])
threshold = np.percentile(natives["sc"], 5)
print(f"5th-percentile SC threshold: {threshold:.3f}")

# Filter design candidates
designs = shape_complementarity.score_many(design_pdbs, chains_a=["H"], chains_b=["A"])
passing = designs[designs["sc"] >= threshold]
```

## Rosetta / PyRosetta validation

PyRosetta is optional and is not a package dependency. In an environment where
`import pyrosetta` works, validate against Rosetta `InterfaceAnalyzerMover`
with:

```bash
python validation/rosetta_interface_sc.py --format tsv
python -m pytest -q -m rosetta
```

The validation script compares bundled 1FYT and nanobody-antigen fixtures
against PyRosetta's interface SC value and fails if `abs(package_sc -
pyrosetta_sc) > 0.05`. If PyRosetta is not installed, the `rosetta` pytest
test is skipped with an explicit message.

---

## Benchmark

Run `python benchmark/speed.py` after `maturin develop --release`. One local
run on an Apple M3 Pro (10-core), using 1FYT chains D vs A (1521+1479 heavy
atoms, ~190 residues each), produced:

```
=== shape_complementarity speed benchmark ===

  [PASS] from_pdb (1FYT, 300+300 residues)           115.2 ms   target: < 200 ms
  [PASS] compute_sc (1521 + 1479 atoms)               63.7 ms   target: < 100 ms
  [PASS] score_many (100 files, 8 workers)           33.2 cx/s   target: maximize

All targets met.
```

The first `from_pdb` call in a fresh process can incur a one-time Rust
shared-library initialization cost. Subsequent calls run at steady-state speed.
`score_many` amortizes this across worker processes.

---

## What this is NOT

- **Not a reimplementation of SC.** The algorithm lives entirely in `sc-rs`.
- **Not a design filter.** Threshold selection and chain naming belong in the consuming pipeline.
- **No asymmetric S_AB / S_BA breakdown.** If `sc-rs` does not expose it, `shape_complementarity` does not invent it.

---

## Acknowledgments

- [sc-rs](https://github.com/cytokineking/sc-rs) (MIT) — the Rust implementation of the Lawrence-Colman SC algorithm that this package wraps.
- Lawrence, M. C. & Colman, P. M. (1993). *Shape complementarity at protein/protein interfaces.* Journal of Molecular Biology, **234**(4), 946–950. https://doi.org/10.1006/jmbi.1993.1648
