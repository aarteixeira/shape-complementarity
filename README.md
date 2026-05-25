# pysc [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Python bindings to [sc-rs](https://github.com/cytokineking/sc-rs) for computing
**Lawrence-Colman Shape Complementarity (SC)** on protein interfaces. Intended
for filtering candidate nanobody:antigen complexes during computational design
pipelines.

`pysc` is **not** a reimplementation of the SC algorithm. All algorithmic
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
pip install .
```

`pip install .` invokes maturin automatically and produces a wheel with the
compiled Rust extension bundled inside.

---

## Quickstart

```python
import pysc

result = pysc.from_pdb("complex.pdb", chains_a=["H"], chains_b=["A"])
print(result.sc)              # e.g. 0.631
print(result.median_distance) # Å
print(result.trimmed_area)    # Å²
```

---

## API Reference

### `from_pdb(pdb_path, chains_a, chains_b=None, *, model=0, include_hetatm=False, include_hydrogens=False, parallel=True) → ScResult`

Compute SC directly from a PDB or mmCIF file.

| Parameter | Default | Description |
|---|---|---|
| `pdb_path` | — | Path to `.pdb`, `.ent`, or `.cif` file |
| `chains_a` | — | List of chain IDs for molecule A |
| `chains_b` | `None` | List of chain IDs for molecule B; `None` = all chains not in `chains_a` |
| `model` | `0` | Model index (0-based); defaults to first model |
| `include_hetatm` | `False` | Include HETATM residues (ligands, ions, water) |
| `include_hydrogens` | `False` | Include hydrogen atoms |
| `parallel` | `True` | Enable Rayon parallelism inside sc-rs |

Raises `ValueError` if a group produces no atoms or the PDB cannot be parsed.

**Note on alternate locations:** Only atoms with altloc `' '` (blank) or `'A'`
are used, matching the sc-rs CLI behavior.

---

### `compute_sc(coords_a, atom_names_a, residue_names_a, coords_b, atom_names_b, residue_names_b, parallel=True) → ScResult`

Low-level binding. Accepts pre-parsed atom arrays directly, bypassing biopython.
Use this when you already hold atom coordinates in memory (e.g., from a design
pipeline) and want to skip file I/O overhead.

```python
from pysc import compute_sc

result = compute_sc(
    coords_a=[[x, y, z], ...],
    atom_names_a=["CA", "CB", ...],
    residue_names_a=["ALA", "ALA", ...],
    coords_b=[[x, y, z], ...],
    atom_names_b=["CA", ...],
    residue_names_b=["GLY", ...],
)
```

Each `coords_*` element is `[x, y, z]` in Ångströms. All three lists for each
group must have the same length. Atoms whose residue+atom-name combination is
not in the sc-rs radius table are silently dropped (same behavior as the CLI).

---

### `score_many(pdb_paths, chains_a, chains_b=None, n_workers=8, parallel=False, **kwargs) → pd.DataFrame`

Parallel batch scoring using `ProcessPoolExecutor`. Returns a DataFrame with
columns: `path`, `sc`, `median_distance`, `trimmed_area`, `atoms_a`, `atoms_b`,
`status` (`'ok'` or `'error'`), `error`.

```python
from pysc import score_many
import numpy as np

df = score_many(pdb_paths, chains_a=["H"], chains_b=["A"], n_workers=8)
print(df[df.status == "ok"][["path", "sc"]].head())
```

`parallel=False` (default) disables Rayon inside each worker to avoid
oversubscription with many worker processes.

---

### `ScResult`

Returned by all compute functions. Read-only properties:

| Property | Type | Description |
|---|---|---|
| `sc` | `float` | Shape complementarity score (−1 to 1; typical protein interfaces: 0.6–0.8) |
| `median_distance` | `float` | Median nearest-surface distance (Å) |
| `trimmed_area` | `float` | Total trimmed interface area (Å²) |
| `atoms_a` | `int` | Atoms accepted for molecule A |
| `atoms_b` | `int` | Atoms accepted for molecule B |

---

## Calibration recipe for design filters

```python
import numpy as np
from pysc import score_many

# Score a set of known native complexes
natives = score_many(native_pdbs, chains_a=["H"], chains_b=["A"])
threshold = np.percentile(natives[natives.status == "ok"]["sc"], 5)
print(f"5th-percentile SC threshold: {threshold:.3f}")

# Apply to design candidates
designs = score_many(design_pdbs, chains_a=["H"], chains_b=["A"])
passing = designs[designs["sc"] >= threshold]
```

---

## Benchmark

Run `python benchmark/speed.py` after `maturin develop --release`. Example output
on an Apple M3 Pro (10-core) with 1FYT chains D vs A (1521+1479 heavy atoms,
~190 residues each) as the test complex:

```
=== pysc speed benchmark ===

  [PASS] from_pdb (1FYT, 300+300 residues)           115.2 ms   target: < 200 ms
  [PASS] compute_sc (1521 + 1479 atoms)               63.7 ms   target: < 100 ms
  [PASS] score_many (100 files, 8 workers)           33.2 cx/s   target: maximize

All targets met.
```

Note: the first `from_pdb` call in a fresh process incurs a one-time ~400 ms
cost for Rust shared-library initialization. Subsequent calls run at steady-state
speed as shown above. `score_many` amortizes this cost across worker processes.

---

## What this is NOT

- **Not a reimplementation of SC.** The algorithm lives entirely in `sc-rs`.
  `pysc` only wraps it.
- **Not a design filter.** Threshold selection, chain naming conventions, and
  filtering logic belong in the consuming pipeline.
- **No asymmetric S_AB / S_BA breakdown.** If `sc-rs` does not expose it,
  `pysc` does not invent it.

---

## Acknowledgments

- [sc-rs](https://github.com/cytokineking/sc-rs) (MIT) — the Rust implementation
  of the Lawrence-Colman SC algorithm that this package binds.
- Lawrence, M. C. & Colman, P. M. (1993). *Shape complementarity at
  protein/protein interfaces.* Journal of Molecular Biology, **234**(4),
  946–950. https://doi.org/10.1006/jmbi.1993.1648
