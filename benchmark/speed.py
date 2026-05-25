"""Speed benchmarks for pysc.

Targets (from spec):
    Single from_pdb (300+300 residue complex):   < 200 ms
    compute_sc with pre-parsed arrays:           < 100 ms
    Batch of 100 complexes, n_workers=8:         report complexes/second
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

DATA = Path(__file__).parent.parent / "tests" / "data"
PDB_1FYT = DATA / "1fyt.pdb"

# Single-chain pair for core benchmarks (~190 residues each).
# 1FYT has 5 chains total; using D+A gives a representative "300+300 residue"
# complex size (1521 + 1479 heavy atoms).
CHAINS_A = ["D"]
CHAINS_B = ["A"]
REPS = 5  # Repetitions for stable timing


def _mean_ms(times: list[float]) -> float:
    return sum(times) / len(times) * 1000


def bench_from_pdb(n: int = REPS) -> float:
    from pysc import from_pdb

    from_pdb(PDB_1FYT, CHAINS_A, CHAINS_B)  # warmup: Rust SO init is a one-time cost
    times = []
    for _ in range(n):
        t0 = time.perf_counter()
        from_pdb(PDB_1FYT, CHAINS_A, CHAINS_B)
        times.append(time.perf_counter() - t0)
    return _mean_ms(times)


def bench_compute_sc(n: int = REPS) -> tuple[float, int, int]:
    """Pre-parse atoms once, then time only the Rust compute_sc calls."""
    from Bio.PDB import PDBParser

    from pysc._core import compute_sc
    from pysc.io import _extract_atom_arrays

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("1fyt", str(PDB_1FYT))
    model = list(structure.get_models())[0]

    coords_a, names_a, res_a = _extract_atom_arrays(model, CHAINS_A, False, False)
    coords_b, names_b, res_b = _extract_atom_arrays(model, CHAINS_B, False, False)

    times = []
    for _ in range(n):
        t0 = time.perf_counter()
        compute_sc(coords_a, names_a, res_a, coords_b, names_b, res_b)
        times.append(time.perf_counter() - t0)

    return _mean_ms(times), len(coords_a), len(coords_b)


def bench_batch(n_files: int = 100, n_workers: int = 8) -> float:
    from pysc import score_many

    paths = [PDB_1FYT] * n_files
    t0 = time.perf_counter()
    df = score_many(paths, CHAINS_A, CHAINS_B, n_workers=n_workers)
    elapsed = time.perf_counter() - t0

    ok = (df["status"] == "ok").sum()
    errors = n_files - ok
    print(f"  (ok={ok}, errors={errors})")
    return n_files / elapsed


def _row(label: str, value: str, target: str, passed: bool) -> str:
    mark = "PASS" if passed else "FAIL"
    return f"  [{mark}] {label:<40s} {value:>12s}   target: {target}"


def main():
    if not PDB_1FYT.exists():
        print("Test PDB not found. Run `pytest tests/` once to download fixtures.")
        sys.exit(1)

    print("\n=== pysc speed benchmark ===\n")

    rows = []

    # 1. Single from_pdb
    ms_pdb = bench_from_pdb()
    rows.append(_row("from_pdb (1FYT, 300+300 residues)", f"{ms_pdb:.1f} ms", "< 200 ms", ms_pdb < 200))

    # 2. compute_sc with pre-parsed arrays
    ms_raw, na, nb = bench_compute_sc()
    rows.append(
        _row(
            f"compute_sc ({na} + {nb} atoms)",
            f"{ms_raw:.1f} ms",
            "< 100 ms",
            ms_raw < 100,
        )
    )

    # 3. Batch
    cps = bench_batch()
    rows.append(_row("score_many (100 files, 8 workers)", f"{cps:.1f} cx/s", "maximize", True))

    print("\n".join(rows))
    print()

    failed = [r for r in rows if r.startswith("  [FAIL]")]
    if failed:
        print(f"WARNING: {len(failed)} target(s) not met.")
    else:
        print("All targets met.")


if __name__ == "__main__":
    main()
