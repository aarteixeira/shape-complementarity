from pathlib import Path

import pytest


DATA = Path(__file__).parent / "data"


@pytest.mark.rosetta
def test_pyrosetta_interface_analyzer_validation():
    pytest.importorskip(
        "pyrosetta",
        reason="PyRosetta is not installed in this environment",
    )
    from validation.rosetta_interface_sc import score_case

    cases = [
        (DATA / "1fyt.pdb", "D", "A"),
        (DATA / "nb_ag_test.pdb", "A", "L"),
    ]
    rows = [score_case(path, chains_a, chains_b) for path, chains_a, chains_b in cases]
    failures = [row for row in rows if row["abs_delta"] > 0.05]
    assert not failures, failures
