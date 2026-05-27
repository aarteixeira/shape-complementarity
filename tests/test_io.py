"""Tests for pysc.io (biopython-based PDB parsing)."""
from __future__ import annotations

from pathlib import Path

import pytest

from shape_complementarity import from_pdb, from_structure
from shape_complementarity._core import ScResult

DATA = Path(__file__).parent / "data"
FYT = DATA / "1fyt.pdb"

# 1FYT: TCR (chains D, E) vs pMHC (chains A, B, C)
CHAINS_A = ["D", "E"]
CHAINS_B = ["A", "B", "C"]


def test_from_pdb_1fyt():
    """Canonical CCP4/Lawrence-Colman test case. SC should be in [0.45, 0.65]."""
    result = from_pdb(FYT, CHAINS_A, CHAINS_B)
    assert isinstance(result, ScResult)
    assert 0.45 <= result.sc <= 0.65, f"SC={result.sc:.4f} outside expected [0.45, 0.65]"
    assert result.atoms_a > 0
    assert result.atoms_b > 0


def test_from_pdb_default_chains_b():
    """chains_b=None should resolve to all chains not in chains_a."""
    r_explicit = from_pdb(FYT, CHAINS_A, CHAINS_B)
    r_implicit = from_pdb(FYT, CHAINS_A, chains_b=None)
    assert r_explicit.sc == pytest.approx(r_implicit.sc, abs=1e-9)
    assert r_explicit.atoms_b == r_implicit.atoms_b


def test_from_pdb_missing_chain_fails_loudly():
    with pytest.raises(ValueError, match="chains_a.*not found"):
        from_pdb(FYT, ["Z"], CHAINS_B)


def test_from_pdb_overlapping_chains_fail_loudly():
    with pytest.raises(ValueError, match="overlapping chain"):
        from_pdb(FYT, ["D"], ["D"])


def test_from_structure():
    """from_structure() with an already-parsed biopython Structure."""
    from Bio.PDB import PDBParser

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("1fyt", str(FYT))
    result = from_structure(structure, CHAINS_A, CHAINS_B)
    assert 0.45 <= result.sc <= 0.65


def test_from_structure_missing_chain_fails_loudly():
    from Bio.PDB import PDBParser

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("1fyt", str(FYT))
    with pytest.raises(ValueError, match="chains_b.*not found"):
        from_structure(structure, ["D"], ["Z"])


def _pdb_atom_line(serial, atom_name, altloc, res_name, chain, resseq, x, y, z, occ, element):
    return (
        f"ATOM  {serial:5d} {atom_name:<4}{altloc}{res_name:>3} {chain}"
        f"{resseq:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{20.0:6.2f}"
        f"          {element:>2}  \n"
    )


def _write_altloc_fixture(path: Path):
    serial = 1
    lines = []
    atom_defs = [
        ("N", "N", -1.2, 0.0, 0.0),
        ("CA", "C", 0.0, 0.0, 0.0),
        ("C", "C", 1.1, 0.0, 0.0),
        ("O", "O", 1.6, 1.0, 0.0),
        ("CB", "C", 0.0, -1.5, 0.0),
    ]
    for chain, y_offset in [("A", 0.0), ("B", 5.0)]:
        for resseq in range(1, 11):
            x_base = (resseq - 1) * 3.8
            for atom_name, element, dx, dy, dz in atom_defs:
                if chain == "A" and resseq == 3 and atom_name == "CB":
                    lines.append(_pdb_atom_line(
                        serial, atom_name, "A", "ALA", chain, resseq,
                        x_base + dx, y_offset + dy, dz, 0.40, element,
                    ))
                    serial += 1
                    lines.append(_pdb_atom_line(
                        serial, atom_name, "B", "ALA", chain, resseq,
                        x_base + dx, y_offset + dy + 0.2, dz, 0.60, element,
                    ))
                    serial += 1
                    continue
                lines.append(_pdb_atom_line(
                    serial, atom_name, " ", "ALA", chain, resseq,
                    x_base + dx, y_offset + dy, dz, 1.00, element,
                ))
                serial += 1
    lines.append("END\n")
    path.write_text("".join(lines))


def test_altloc_default_fails(tmp_path):
    pdb = tmp_path / "altloc.pdb"
    _write_altloc_fixture(pdb)
    with pytest.raises(ValueError, match="alternate locations"):
        from_pdb(pdb, ["A"], ["B"])


def test_altloc_highest_occupancy_and_sc_rs_are_explicit(tmp_path):
    pdb = tmp_path / "altloc.pdb"
    _write_altloc_fixture(pdb)
    highest = from_pdb(pdb, ["A"], ["B"], altloc_policy="highest_occupancy")
    sc_rs = from_pdb(pdb, ["A"], ["B"], altloc_policy="sc_rs")
    assert highest.atoms_a == sc_rs.atoms_a == 50
    assert highest.atoms_b == sc_rs.atoms_b == 50


def test_hetatm_excluded_by_default():
    """Adding HETATM atoms to the PDB must not change the SC result."""
    import tempfile

    # Read original PDB
    lines = FYT.read_text().splitlines(keepends=True)

    # Append a fake HETATM line (SO4 sulfate ion far from the interface)
    fake_hetatm = (
        "HETATM 9999  S   SO4 A 999      99.000  99.000  99.000  1.00  0.00           S  \n"
    )
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".pdb", delete=False, dir=DATA
    ) as tmp:
        tmp.writelines(lines)
        tmp.write(fake_hetatm)
        tmp_path = tmp.name

    try:
        r_orig = from_pdb(FYT, CHAINS_A, CHAINS_B)
        r_mod = from_pdb(tmp_path, CHAINS_A, CHAINS_B)
        # HETATM excluded by default → same result
        assert r_orig.sc == pytest.approx(r_mod.sc, abs=1e-9)
    finally:
        Path(tmp_path).unlink(missing_ok=True)


def test_model_index():
    """Requesting model=0 (default) and explicitly model=0 should agree."""
    r1 = from_pdb(FYT, CHAINS_A, CHAINS_B)
    r2 = from_pdb(FYT, CHAINS_A, CHAINS_B, model=0)
    assert r1.sc == pytest.approx(r2.sc, abs=1e-9)


def test_invalid_model_index():
    """Out-of-range model index should raise ValueError."""
    with pytest.raises(ValueError, match="model index"):
        from_pdb(FYT, CHAINS_A, CHAINS_B, model=999)
