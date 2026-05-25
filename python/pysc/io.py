"""PDB/CIF parsing for pysc.

Parsing logic mirrors src/bin/sc.rs in sc-rs exactly so that parity tests pass:
- ATOM records only (HETATM skipped unless include_hetatm=True)
- Alternate locations: keep ' ' and 'A', skip all others
- Hydrogens: excluded by default using the same heuristic as the CLI
"""
from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from Bio.PDB.Structure import Structure as BioStructure

from pysc._core import ScResult, compute_sc


def _is_hydrogen(atom_name: str, element: str) -> bool:
    """Mirror the hydrogen-detection logic in sc-rs bin/sc.rs."""
    elem = element.strip().upper()
    name = atom_name.strip()
    if elem == "H":
        return True
    if name.startswith("H"):
        return True
    if name.endswith("H"):
        return True
    # Catch names like "1H", "2HB" (digit-prefixed hydrogen names)
    if "H" in name and name[:1].isdigit():
        return True
    return False


def _select_real_atom(atom):
    """Return a concrete (non-disordered) Atom for altloc ' ' or 'A'.

    Returns None if no acceptable altloc exists.
    """
    if atom.is_disordered():
        child_dict = atom.child_dict
        for altloc in ("A", " "):
            if altloc in child_dict:
                return child_dict[altloc]
        return None
    altloc = atom.altloc
    if altloc not in (" ", "A"):
        return None
    return atom


def _extract_atom_arrays(
    model,
    chains: list[str],
    include_hetatm: bool,
    include_hydrogens: bool,
) -> tuple[list[list[float]], list[str], list[str]]:
    coords: list[list[float]] = []
    atom_names: list[str] = []
    res_names: list[str] = []

    for chain in model.get_chains():
        if chain.id not in chains:
            continue
        for residue in chain.get_residues():
            # residue.id[0] == ' ' for standard ATOM records
            het = residue.id[0]
            if not include_hetatm and het != " ":
                continue
            for disordered_or_atom in residue.get_atoms():
                real = _select_real_atom(disordered_or_atom)
                if real is None:
                    continue
                atom_name = real.name.strip()
                element = (real.element or "").strip()
                if not include_hydrogens and _is_hydrogen(atom_name, element):
                    continue
                c = real.coord
                coords.append([float(c[0]), float(c[1]), float(c[2])])
                atom_names.append(atom_name)
                res_names.append(residue.resname.strip())

    return coords, atom_names, res_names


def _load_structure(path: str | Path):
    from Bio.PDB import MMCIFParser, PDBParser

    path = Path(path)
    suffix = path.suffix.lower()
    if suffix in (".cif", ".mmcif"):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    return parser.get_structure(path.stem, str(path))


def from_structure(
    structure: "BioStructure",
    chains_a: list[str],
    chains_b: list[str] | None = None,
    model: int = 0,
    include_hetatm: bool = False,
    include_hydrogens: bool = False,
    parallel: bool = True,
) -> ScResult:
    """Compute SC from a biopython Structure object.

    Args:
        structure:          biopython Structure (any source)
        chains_a:           chain IDs for molecule A
        chains_b:           chain IDs for molecule B; None = all chains not in chains_a
        model:              model index (0-based)
        include_hetatm:     include HETATM residues (default False, matches sc-rs)
        include_hydrogens:  include hydrogen atoms (default False, matches sc-rs)
        parallel:           enable Rayon parallelism inside sc-rs
    """
    models = list(structure.get_models())
    if model >= len(models):
        raise ValueError(
            f"model index {model} out of range (structure has {len(models)} model(s))"
        )
    m = models[model]

    if chains_b is None:
        all_chain_ids = {ch.id for ch in m.get_chains()}
        chains_b = sorted(all_chain_ids - set(chains_a))

    coords_a, names_a, res_a = _extract_atom_arrays(m, chains_a, include_hetatm, include_hydrogens)
    coords_b, names_b, res_b = _extract_atom_arrays(m, chains_b, include_hetatm, include_hydrogens)

    return compute_sc(coords_a, names_a, res_a, coords_b, names_b, res_b, parallel)


def from_pdb(
    pdb_path: str | Path,
    chains_a: list[str],
    chains_b: list[str] | None = None,
    model: int = 0,
    include_hetatm: bool = False,
    include_hydrogens: bool = False,
    parallel: bool = True,
) -> ScResult:
    """Compute SC from a PDB or mmCIF file.

    Args:
        pdb_path:           path to .pdb, .ent, or .cif file
        chains_a:           chain IDs for molecule A
        chains_b:           chain IDs for molecule B; None = all chains not in chains_a
        model:              model index (0-based, default first model)
        include_hetatm:     include HETATM residues (default False)
        include_hydrogens:  include hydrogen atoms (default False)
        parallel:           enable Rayon parallelism inside sc-rs
    """
    structure = _load_structure(pdb_path)
    return from_structure(
        structure,
        chains_a,
        chains_b,
        model=model,
        include_hetatm=include_hetatm,
        include_hydrogens=include_hydrogens,
        parallel=parallel,
    )
