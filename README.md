# shape_complementarity [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Python bindings to [sc-rs](https://github.com/cytokineking/sc-rs) for computing
**Lawrence-Colman Shape Complementarity (SC)** between two protein chains.
Useful for analysing protein–protein interfaces, including in computational
design pipelines.

`shape_complementarity` is **not** a reimplementation of the SC algorithm. All algorithmic
correctness comes from the upstream Rust crate `sc-rs`, which is invoked
directly via PyO3. If you find a numerical result surprising, verify against the
`sc-rs` CLI first.

---

## Install

**Development (editable, re-run after Rust changes):**
```bash
pip install maturin
maturin develop --release
```

**Release / production:**
```bash
pip install shape-complementarity
```

---

## Usage

### 1. From a PDB or mmCIF file path

The simplest entry point. Accepts `.pdb`, `.ent`, and `.cif` files.

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
analysis stack. Pass an `AtomArray` or `AtomArrayStack` (first model is used).

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
Atoms whose combination is not in the sc-rs radius table are silently dropped
(same behavior as the sc-rs CLI).

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

**Option C — from an in-memory BoltzGen `Structure` object:**

No boltzgen import required; the function duck-types against the numpy
structured-array layout of `boltzgen.data.data.Structure`.

```python
import shape_complementarity

# structure is a boltzgen.data.data.Structure loaded elsewhere in the pipeline
result = shape_complementarity.from_boltzgen_structure(
    structure,
    chains_a=["B"],
    chains_b=["A"],
)
```

> **Which stage to use?** Only `refold_cif/` structures have complete, physically
> validated all-atom coordinates. `intermediate_designs/` NPZ files (post-generation)
> have backbone-only coordinates; `intermediate_designs_inverse_folded/` NPZ files
> have sidechains but have not yet been validated by Boltz.

---

### 6. Batch scoring

Score many files in parallel using `ProcessPoolExecutor`. Returns a
`pd.DataFrame` with one row per file; exceptions are caught per-file and
reported in the `status` / `error` columns rather than crashing the batch.

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

## ScResult

All functions return an `ScResult` with these read-only properties:

| Property | Type | Description |
|---|---|---|
| `sc` | `float` | Shape complementarity (−1 to 1; native interfaces typically 0.6–0.8) |
| `median_distance` | `float` | Median nearest-surface distance (Å) |
| `trimmed_area` | `float` | Total trimmed interface area (Å²) |
| `atoms_a` | `int` | Heavy atoms accepted for molecule A |
| `atoms_b` | `int` | Heavy atoms accepted for molecule B |

---

## Calibration recipe for design filters

```python
import numpy as np
import shape_complementarity

# Score known native complexes to establish a baseline
natives = shape_complementarity.score_many(native_pdbs, chains_a=["H"], chains_b=["A"])
threshold = np.percentile(natives[natives.status == "ok"]["sc"], 5)
print(f"5th-percentile SC threshold: {threshold:.3f}")

# Filter design candidates
designs = shape_complementarity.score_many(design_pdbs, chains_a=["H"], chains_b=["A"])
passing = designs[designs["sc"] >= threshold]
```

---

## Benchmark

Run `python benchmark/speed.py` after `maturin develop --release`. Example output
on an Apple M3 Pro (10-core) with 1FYT chains D vs A (1521+1479 heavy atoms,
~190 residues each):

```
=== shape_complementarity speed benchmark ===

  [PASS] from_pdb (1FYT, 300+300 residues)           115.2 ms   target: < 200 ms
  [PASS] compute_sc (1521 + 1479 atoms)               63.7 ms   target: < 100 ms
  [PASS] score_many (100 files, 8 workers)           33.2 cx/s   target: maximize

All targets met.
```

The first `from_pdb` call in a fresh process incurs a one-time ~400 ms cost for
Rust shared-library initialisation. Subsequent calls run at steady-state speed
as shown above. `score_many` amortises this across worker processes.

---

## What this is NOT

- **Not a reimplementation of SC.** The algorithm lives entirely in `sc-rs`.
- **Not a design filter.** Threshold selection and chain naming belong in the consuming pipeline.
- **No asymmetric S_AB / S_BA breakdown.** If `sc-rs` does not expose it, `shape_complementarity` does not invent it.

---

## Acknowledgments

- [sc-rs](https://github.com/cytokineking/sc-rs) (MIT) — the Rust implementation of the Lawrence-Colman SC algorithm that this package wraps.
- Lawrence, M. C. & Colman, P. M. (1993). *Shape complementarity at protein/protein interfaces.* Journal of Molecular Biology, **234**(4), 946–950. https://doi.org/10.1006/jmbi.1993.1648
