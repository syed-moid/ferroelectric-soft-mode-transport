"""fekappa.provenance: the three-line header of the derived tables."""
import gzip

import pandas as pd
import pytest

from fekappa import __version__
from fekappa.provenance import external_dir, stamp


def test_stamp_prepends_and_replaces_header(tmp_path):
    p = tmp_path / "t.csv"
    p.write_text("# units: W/mK\na,b\n1,2\n")
    stamp(p, __file__, ["in1.csv", "zenodo:x.tar.gz"])
    stamp(p, __file__, ["in1.csv", "zenodo:x.tar.gz"])  # restamping replaces, does not stack
    lines = p.read_text().splitlines()
    assert lines[0] == "# script: tests/test_provenance.py"
    assert lines[1] == "# inputs: in1.csv; zenodo:x.tar.gz"
    assert lines[2].startswith("# generated: ") and f"fekappa {__version__}" in lines[2]
    assert lines[3:] == ["# units: W/mK", "a,b", "1,2"]
    assert pd.read_csv(p, comment="#").to_dict("list") == {"a": [1], "b": [2]}


def test_stamp_gzip(tmp_path):
    p = tmp_path / "t.csv.gz"
    p.write_bytes(gzip.compress(b"a\n1\n"))
    stamp(p, __file__, ["x"])
    assert gzip.decompress(p.read_bytes()).decode().splitlines()[3:] == ["a", "1"]


def test_external_dir(monkeypatch, tmp_path):
    monkeypatch.delenv("FEKAPPA_TEST_DIR", raising=False)
    with pytest.raises(SystemExit):
        external_dir("FEKAPPA_TEST_DIR", "a test directory")
    monkeypatch.setenv("FEKAPPA_TEST_DIR", str(tmp_path))
    assert external_dir("FEKAPPA_TEST_DIR", "a test directory") == tmp_path
