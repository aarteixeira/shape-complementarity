"""Validate package SC values against PyRosetta InterfaceAnalyzerMover.

PyRosetta is intentionally not a package dependency. Run this script from an
environment where ``import pyrosetta`` succeeds.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path


DEFAULT_CASES = [
    ("tests/data/1fyt.pdb", "D", "A"),
    ("tests/data/nb_ag_test.pdb", "A", "L"),
]

_PYROSETTA = None
_ROSETTA = None


def _init_pyrosetta():
    global _PYROSETTA, _ROSETTA
    if _PYROSETTA is not None and _ROSETTA is not None:
        return _PYROSETTA, _ROSETTA
    try:
        import pyrosetta
        from pyrosetta import rosetta
    except ImportError as exc:
        raise ImportError(
            "PyRosetta is required for Rosetta validation. Install it in the "
            "active environment, then rerun this script."
        ) from exc

    pyrosetta.init("-mute all -ignore_unrecognized_res true")
    _PYROSETTA = pyrosetta
    _ROSETTA = rosetta
    return _PYROSETTA, _ROSETTA


def _require_method(obj, name: str):
    method = getattr(obj, name, None)
    if method is None:
        raise RuntimeError(f"PyRosetta object {type(obj).__name__} lacks {name}()")
    return method


def _pyrosetta_interface_sc(pdb_path: Path, chains_a: str, chains_b: str) -> float:
    pyrosetta, rosetta = _init_pyrosetta()
    pose = pyrosetta.pose_from_pdb(str(pdb_path))
    scorefxn = pyrosetta.get_fa_scorefxn()
    mover = rosetta.protocols.analysis.InterfaceAnalyzerMover()

    _require_method(mover, "set_interface")(f"{chains_a}_{chains_b}")
    _require_method(mover, "set_scorefunction")(scorefxn)
    _require_method(mover, "set_compute_interface_sc")(True)
    _require_method(mover, "set_pack_input")(False)
    _require_method(mover, "set_pack_separated")(False)
    _require_method(mover, "set_compute_packstat")(False)
    mover.apply(pose)

    data = _require_method(mover, "get_all_data")()
    if not hasattr(data, "sc_value"):
        raise RuntimeError("InterfaceAnalyzerMover data lacks sc_value")
    return float(data.sc_value)


def score_case(pdb_path: Path, chains_a: str, chains_b: str) -> dict:
    import shape_complementarity

    package = shape_complementarity.from_pdb(
        pdb_path,
        chains_a=list(chains_a),
        chains_b=list(chains_b),
    )
    rosetta_sc = _pyrosetta_interface_sc(pdb_path, chains_a, chains_b)
    return {
        "path": str(pdb_path),
        "chains_a": chains_a,
        "chains_b": chains_b,
        "package_sc": package.sc,
        "pyrosetta_sc": rosetta_sc,
        "abs_delta": abs(package.sc - rosetta_sc),
        "median_distance": package.median_distance,
        "trimmed_area": package.trimmed_area,
        "atoms_a": package.atoms_a,
        "atoms_b": package.atoms_b,
    }


def _parse_case(value: str) -> tuple[Path, str, str]:
    parts = value.split(":")
    if len(parts) != 3:
        raise argparse.ArgumentTypeError(
            "cases must have format /path/to/file.pdb:CHAINS_A:CHAINS_B"
        )
    return Path(parts[0]), parts[1], parts[2]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case",
        action="append",
        type=_parse_case,
        help=(
            "Case as /path/to/file.pdb:CHAINS_A:CHAINS_B. "
            "May be supplied more than once. Defaults to bundled fixtures."
        ),
    )
    parser.add_argument("--format", choices=("tsv", "json"), default="tsv")
    parser.add_argument("--max-delta", type=float, default=0.05)
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    cases = args.case or [
        (repo_root / pdb, chains_a, chains_b)
        for pdb, chains_a, chains_b in DEFAULT_CASES
    ]

    rows = [score_case(Path(pdb), chains_a, chains_b) for pdb, chains_a, chains_b in cases]
    if args.format == "json":
        print(json.dumps(rows, indent=2, sort_keys=True))
    else:
        columns = list(rows[0])
        print("\t".join(columns))
        for row in rows:
            print("\t".join(str(row[col]) for col in columns))

    failures = [row for row in rows if row["abs_delta"] > args.max_delta]
    if failures:
        raise SystemExit(
            f"{len(failures)} case(s) exceeded max delta {args.max_delta}: "
            + "; ".join(f"{row['path']} delta={row['abs_delta']:.4f}" for row in failures)
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
