"""Offline tests for the FULL horse coat-colour CELL (dna_decode/pigment/horse_coat.py) + `dna-decode
horsecolor`. Distinct from tests/test_horse_coat.py (which pins the base E x A rule in data/horse_coat.py
that this cell REUSES).

Synthetic allele calls only — no data / no network. Pins the OMIA-curated epistatic rule, the TWO epistasis
anchors a naive rule mis-calls (e/e hides Agouti; grey is epistatic), cream dose-dependence, dun, abstention,
and CLI dispatch. Runnable via pytest OR standalone.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pytest  # noqa: E402

import dna_decode.cli as uni  # noqa: E402
from dna_decode.pigment.horse_coat import (  # noqa: E402
    HorseInputError,
    call_horse_coat,
    reference_integrity_ok,
)
from dna_decode.pigment.horse_coat_cli import main as horse_main  # noqa: E402


def test_reference_integrity():
    assert reference_integrity_ok() is True


def test_base_colors():
    assert call_horse_coat({"E": "E/e", "A": "A/a"}).coat_color == "bay"
    assert call_horse_coat({"E": "E/E", "A": "a/a"}).coat_color == "black"
    assert call_horse_coat({"E": "e/e", "A": "a/a"}).coat_color == "chestnut"


def test_anchor_ee_hides_agouti():
    c = call_horse_coat({"E": "e/e", "A": "A/A"})
    assert c.coat_color == "chestnut" and c.base_color == "chestnut"
    assert any("MIS-CALL" in n for n in c.notes)


def test_anchor_grey_epistatic():
    c = call_horse_coat({"E": "E/e", "A": "A/a", "G": "G/n"})
    assert c.greying is True and c.coat_color.startswith("grey (born bay") and c.base_color == "bay"


def test_cream_dose_dependent():
    assert call_horse_coat({"E": "e/e", "A": "a/a", "CR": "Cr/N"}).coat_color == "palomino"
    assert call_horse_coat({"E": "E/e", "A": "A/a", "CR": "Cr/N"}).coat_color == "buckskin"
    assert call_horse_coat({"E": "E/E", "A": "a/a", "CR": "Cr/N"}).coat_color == "smoky black"
    assert call_horse_coat({"E": "e/e", "A": "a/a", "CR": "Cr/Cr"}).coat_color.startswith("cremello")
    assert call_horse_coat({"E": "E/e", "A": "A/a", "CR": "Cr/Cr"}).coat_color.startswith("perlino")


def test_dun():
    assert call_horse_coat({"E": "E/E", "A": "a/a", "D": "D/nd1"}).coat_color.startswith("grullo")
    assert call_horse_coat({"E": "e/e", "A": "a/a", "D": "D/D"}).coat_color == "red dun"
    assert "dun" not in call_horse_coat({"E": "E/e", "A": "A/a", "D": "nd1/nd2"}).dilutions


def test_e_required():
    with pytest.raises(HorseInputError):
        call_horse_coat({"A": "A/a"})


def test_bay_black_unknown_without_agouti():
    c = call_horse_coat({"E": "E/e"})
    assert c.base_color == "undetermined" and c.confidence == "medium"
    assert any("A:" in a for a in c.abstains_on)


def test_reuses_deployed_base_rule():
    # the base E x A must AGREE with the deployed data/horse_coat rule it reuses
    from dna_decode.data.horse_coat import call_horse_base_colour
    for e, a, expect_base in [("Ee", "Aa", "bay"), ("EE", "aa", "black"), ("ee", "AA", "chestnut")]:
        assert call_horse_base_colour(e, a) == expect_base
        got = call_horse_coat({"E": f"{e[0]}/{e[1]}", "A": f"{a[0]}/{a[1]}"})
        assert (got.base_color if got.base_color != "chestnut" else "chestnut") == expect_base


def test_unmodelled_locus_abstains():
    c = call_horse_coat({"E": "E/e", "A": "A/a"}, present_loci=["Z", "TO"])
    assert any("Z:" in a for a in c.abstains_on) and any("TO:" in a for a in c.abstains_on)
    assert c.confidence == "medium"


def test_unmodelled_locus_in_loci_raises():
    with pytest.raises(HorseInputError):
        call_horse_coat({"E": "E/e", "Z": "Z/n"})


def test_cli_json_and_dispatch(capsys):
    rc = horse_main(["--loci", "E=E/e,A=A/a,CR=Cr/N", "--json"])
    assert rc == 0
    d = json.loads(capsys.readouterr().out)
    assert d["organism"] == "Equus_caballus" and d["coat_color"] == "buckskin"
    assert "knowledge_baseline" in d["evidence_tier"]
    captured = {}

    def fake(argv):
        captured["argv"] = argv
        return 0

    import dna_decode.pigment.horse_coat_cli as hc
    orig = hc.main
    hc.main = fake
    try:
        assert uni.main(["horsecolor", "--loci", "E=E/e,A=A/a"]) == 0
    finally:
        hc.main = orig
    assert captured["argv"] == ["--loci", "E=E/e,A=A/a"]


def test_cli_bad_input_returns_2():
    assert horse_main(["--loci", "A=A/a"]) == 2
    assert horse_main(["--loci", "E=X/x"]) == 2


if __name__ == "__main__":
    print("run via pytest (uses capsys fixture)")
