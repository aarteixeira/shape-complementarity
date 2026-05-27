from pathlib import Path

import pytest

from shape_complementarity import score_many


DATA = Path(__file__).parent / "data"
FYT = DATA / "1fyt.pdb"


def test_score_many_raises_by_default_on_failed_file(tmp_path):
    missing = tmp_path / "missing.pdb"
    with pytest.raises(ValueError, match="batch failed.*missing.pdb"):
        score_many([FYT, missing], ["D"], ["A"], n_workers=1)


def test_score_many_records_errors_when_requested(tmp_path):
    missing = tmp_path / "missing.pdb"
    df = score_many([FYT, missing], ["D"], ["A"], n_workers=1, on_error="record")
    assert list(df["status"]) == ["ok", "error"]
    assert df.loc[1, "path"].endswith("missing.pdb")
    assert df.loc[1, "error"]
