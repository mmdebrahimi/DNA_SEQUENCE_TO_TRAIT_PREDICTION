"""Offline tests for `dna-decode morphology` (dna_decode/pigment/morphology_cli.py).

No genotype data / no network — synthetic dosages only. Pins: size polygenic call + ear call + abstention
list + JSON shape + dispatch from the unified CLI. Runnable via pytest OR standalone.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pytest  # noqa: E402

import dna_decode.cli as uni  # noqa: E402
from dna_decode.pigment.morphology_cli import _parse_dosages, main  # noqa: E402


def test_parse_dosages():
    assert _parse_dosages("IGF1=2,HMGA2=1,EAR=0") == {"IGF1": 2, "HMGA2": 1, "EAR": 0}
    with pytest.raises(ValueError):
        _parse_dosages("IGF1")            # no '='
    with pytest.raises(ValueError):
        _parse_dosages("IGF1=x")          # non-int


def test_full_panel_large_erect(capsys):
    rc = main(["--dosages", "IGF1=2,HMGA2=2,STC2=2,GHR=2,EAR=2", "--json"])
    assert rc == 0
    d = json.loads(capsys.readouterr().out)
    assert d["trait"] == "morphology"
    assert d["height"]["size_rank"] == "large/giant" and d["height"]["polygenic_score"] == 8
    assert d["ear"]["ear_type"].startswith("erect") and d["ear"]["confidence"] == "medium"
    assert len(d["abstains_on"]) == 3


def test_all_zero_small_drop(capsys):
    rc = main(["--dosages", "IGF1=0,HMGA2=0,STC2=0,GHR=0,EAR=0", "--json"])
    assert rc == 0
    d = json.loads(capsys.readouterr().out)
    assert d["height"]["size_rank"] == "toy/small"
    assert d["ear"]["ear_type"].startswith("drop")


def test_partial_size_panel_no_ear(capsys):
    rc = main(["--dosages", "HMGA2=1,IGF1=1", "--json"])
    assert rc == 0
    d = json.loads(capsys.readouterr().out)
    assert d["height"]["n_loci_scored"] == 2 and d["ear"] is None


def test_ear_only(capsys):
    rc = main(["--dosages", "EAR=1", "--json"])
    assert rc == 0
    d = json.loads(capsys.readouterr().out)
    assert d["height"] is None and d["ear"]["ear_type"].startswith("semi")


def test_bad_dose_and_unknown_locus_return_2():
    assert main(["--dosages", "IGF1=3"]) == 2         # dose out of range
    assert main(["--dosages", "NOSUCH=1"]) == 2       # unknown locus
    assert main(["--dosages", ""]) == 2               # nothing scorable


def test_dispatch_from_unified_cli():
    captured = {}

    def fake_main(argv):
        captured["argv"] = argv
        return 0

    import dna_decode.pigment.morphology_cli as mcli
    orig = mcli.main
    mcli.main = fake_main
    try:
        rc = uni.main(["morphology", "--dosages", "HMGA2=2,EAR=1"])
    finally:
        mcli.main = orig
    assert rc == 0 and captured["argv"] == ["--dosages", "HMGA2=2,EAR=1"]


if __name__ == "__main__":
    print("run via pytest (uses capsys fixture)")
