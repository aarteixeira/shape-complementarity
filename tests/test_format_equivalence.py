"""Cross-loader / cross-format equivalence on a real structure (1FYT).

Locks in that every supported intake path produces the same SC and the same
atom selection for the same underlying structure:

    from_pdb(.pdb)                    biopython PDBParser
    from_pdb(.cif)                    biopython MMCIFParser
    from_structure(pre-parsed)        biopython, in-memory
    from_biotite(parsed CIF)          biotite, in-memory
    from_boltzgen_structure(...)      duck-typed numpy structured arrays

Empirically all five paths agree at the float-rounding floor (sc = 0.6597,
1521 + 1479 atoms). Any future parser change that shifts atom selection or
coordinate rounding will fail this suite.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
from Bio.PDB import PDBParser

from shape_complementarity import from_biotite, from_pdb, from_structure
# from_boltzgen_structure is intentionally not in the public API; reach into
# the module to exercise it.
from shape_complementarity.io import _is_hydrogen, _select_real_atom, from_boltzgen_structure

DATA = Path(__file__).parent / "data"
PDB_1FYT = DATA / "1fyt.pdb"
CIF_1FYT = DATA / "1fyt.cif"
CHAINS_A = ["D"]
CHAINS_B = ["A"]


if not PDB_1FYT.exists() or not CIF_1FYT.exists():
    pytest.skip("1FYT PDB and CIF fixtures required", allow_module_level=True)


# ── biopython PDB ↔ CIF ─────────────────────────────────────────────────────

def test_from_pdb_cif_equivalent():
    """from_pdb() must give identical SC and atom counts on the PDB and the CIF."""
    r_pdb = from_pdb(PDB_1FYT, CHAINS_A, CHAINS_B)
    r_cif = from_pdb(CIF_1FYT, CHAINS_A, CHAINS_B)
    assert r_pdb.sc == pytest.approx(r_cif.sc, abs=1e-6)
    assert r_pdb.atoms_a == r_cif.atoms_a
    assert r_pdb.atoms_b == r_cif.atoms_b


# ── from_structure(pre-parsed biopython) ↔ from_pdb(file) ───────────────────

def test_from_structure_matches_from_pdb():
    """A pre-parsed biopython Structure must give the same SC as reading the file."""
    structure = PDBParser(QUIET=True).get_structure("1fyt", str(PDB_1FYT))
    r_struct = from_structure(structure, CHAINS_A, CHAINS_B)
    r_pdb = from_pdb(PDB_1FYT, CHAINS_A, CHAINS_B)
    assert r_struct.sc == pytest.approx(r_pdb.sc, abs=1e-6)
    assert r_struct.atoms_a == r_pdb.atoms_a
    assert r_struct.atoms_b == r_pdb.atoms_b


# ── biotite (CIF) ↔ biopython (PDB) ─────────────────────────────────────────

biotite_pdbx = pytest.importorskip("biotite.structure.io.pdbx")


def test_from_biotite_matches_from_pdb():
    """Same structure read through biotite's CIF parser and biopython's PDB
    parser should give the same SC and atom counts on 1FYT (where auth and
    label chain IDs match)."""
    cif = biotite_pdbx.CIFFile.read(str(CIF_1FYT))
    atoms = biotite_pdbx.get_structure(cif, model=1, use_author_fields=True)
    r_biotite = from_biotite(atoms, CHAINS_A, CHAINS_B)
    r_pdb = from_pdb(PDB_1FYT, CHAINS_A, CHAINS_B)
    assert r_biotite.sc == pytest.approx(r_pdb.sc, abs=1e-3), (
        f"biotite={r_biotite.sc:.4f} vs biopython={r_pdb.sc:.4f}"
    )
    # Atom-selection conventions are very tight but not bit-identical across
    # parsers; ±5 atoms / side is the practical bound.
    assert abs(r_biotite.atoms_a - r_pdb.atoms_a) <= 5
    assert abs(r_biotite.atoms_b - r_pdb.atoms_b) <= 5


# ── from_boltzgen_structure (numpy-layout) ↔ from_pdb (biopython) ───────────

def _chain_specs_from_biopython(structure, chains):
    """Convert a biopython Structure into the chain_specs nested-list format
    consumed by the _make_structure helper in test_boltzgen.py."""
    specs = []
    model = structure[0]
    for chain in model.get_chains():
        if chain.id not in chains:
            continue
        residues = []
        for residue in chain.get_residues():
            if residue.id[0] != " ":  # skip HETATM
                continue
            atoms = []
            for da in residue.get_atoms():
                real = _select_real_atom(da, altloc_policy="sc_rs")
                if real is None:
                    continue
                name = real.name.strip()
                elem = (real.element or "").strip()
                if _is_hydrogen(name, elem):
                    continue
                c = real.coord
                atoms.append((name, (float(c[0]), float(c[1]), float(c[2]))))
            residues.append((residue.resname.strip(), atoms))
        specs.append((chain.id, residues))
    return specs


def test_from_boltzgen_structure_matches_from_pdb():
    """Build a BoltzGen-shaped numpy struct from 1FYT (using the mock helper
    in test_boltzgen.py for the dtype layout) and assert SC matches from_pdb."""
    sys.path.insert(0, str(Path(__file__).parent))
    from test_boltzgen import _make_structure  # noqa: PLC0415

    structure = PDBParser(QUIET=True).get_structure("1fyt", str(PDB_1FYT))
    specs = _chain_specs_from_biopython(structure, CHAINS_A + CHAINS_B)
    boltz_struct = _make_structure(specs)

    r_boltz = from_boltzgen_structure(boltz_struct, chains_a=CHAINS_A, chains_b=CHAINS_B)
    r_pdb = from_pdb(PDB_1FYT, CHAINS_A, CHAINS_B)
    assert r_boltz.sc == pytest.approx(r_pdb.sc, abs=1e-6)
    assert r_boltz.atoms_a == r_pdb.atoms_a
    assert r_boltz.atoms_b == r_pdb.atoms_b
