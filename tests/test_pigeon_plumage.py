"""Offline tests for the pigeon plumage cell (dna_decode/pigment/pigeon_plumage.py) + `dna-decode pigeoncolor`.

Synthetic allele calls only — no data / no network. Pins the Shapiro-lab-sourced epistatic rule, the anchors
a naive rule mis-calls (SOX10 e/e recessive-red epistatic over TYRP1; Z-linked reversed-hemizygous B locus),
the base series, dilute, wing pattern, abstention, and CLI dispatch. Runnable via pytest OR standalone.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pytest  # noqa: E402

import dna_decode.cli as uni  # noqa: E402
from dna_decode.pigment.pigeon_plumage import (  # noqa: E402
    PigeonInputError,
    call_pigeon_plumage,
    reference_integrity_ok,
)
from dna_decode.pigment.pigeon_plumage_cli import main as pigeon_main  # noqa: E402


def test_reference_integrity():
    assert reference_integrity_ok() is True


def test_base_series():
    assert call_pigeon_plumage({"B": "BA/B+"}).base_color == "ash-red"      # dominant
    assert call_pigeon_plumage({"B": "B+/b"}).base_color == "blue/black"    # brown recessive hidden
    assert call_pigeon_plumage({"B": "b/b"}).base_color == "brown"


def test_anchor_recessive_red_epistatic():
    # SOX10 e/e -> red REGARDLESS of a blue TYRP1 base
    c = call_pigeon_plumage({"B": "B+/B+", "E": "e/e"})
    assert c.is_recessive_red and c.base_color == "red (recessive)"
    assert any("mis-calls" in n for n in c.notes)


def test_anchor_zlinked_female_hemizygous():
    hen = call_pigeon_plumage({"B": "BA"})            # 1 allele -> female (ZW)
    cock = call_pigeon_plumage({"B": "BA/B+"})        # 2 alleles -> male (ZZ)
    assert "female" in hen.sex_basis and hen.base_color == "ash-red"
    assert "male" in cock.sex_basis


def test_dilute():
    assert call_pigeon_plumage({"B": "B+/B+", "D": "d/d"}).plumage.startswith("dun")       # dilute blue
    assert call_pigeon_plumage({"B": "BA/B+", "D": "d/d"}).plumage.startswith("ash-yellow")  # dilute ash-red
    assert call_pigeon_plumage({"B": "b/b", "D": "d/d"}).plumage.startswith("khaki")        # dilute brown


def test_wing_pattern():
    assert call_pigeon_plumage({"B": "B+/B+", "C": "CT/+"}).wing_pattern == "T-check"
    assert call_pigeon_plumage({"B": "B+/B+", "C": "C/+"}).wing_pattern == "checker"
    assert call_pigeon_plumage({"B": "B+/B+", "C": "+/+"}).wing_pattern == "bar"
    barless = call_pigeon_plumage({"B": "B+/B+", "C": "c/c"})
    assert barless.wing_pattern == "barless" and any("vision" in n for n in barless.notes)


def test_unmodelled_locus_abstains_and_raises():
    c = call_pigeon_plumage({"B": "B+/B+"}, present_loci=["S", "G"])
    assert any("S:" in a for a in c.abstains_on) and any("G:" in a for a in c.abstains_on)
    with pytest.raises(PigeonInputError):
        call_pigeon_plumage({"B": "B+/B+", "S": "S/s"})


def test_cli_json_and_dispatch(capsys):
    rc = pigeon_main(["--loci", "B=B+/B+,E=e/e", "--json"])
    assert rc == 0
    d = json.loads(capsys.readouterr().out)
    assert d["organism"] == "Columba_livia" and d["is_recessive_red"] is True
    captured = {}

    def fake(argv):
        captured["argv"] = argv
        return 0

    import dna_decode.pigment.pigeon_plumage_cli as pc
    orig = pc.main
    pc.main = fake
    try:
        assert uni.main(["pigeoncolor", "--loci", "B=B+/B+"]) == 0
    finally:
        pc.main = orig
    assert captured["argv"] == ["--loci", "B=B+/B+"]


def test_cli_bad_input_returns_2():
    assert pigeon_main(["--loci", "B=Z/Z"]) == 2               # unknown allele
    assert pigeon_main(["--loci", "E=E+/e/e"]) == 2            # 3 alleles at an autosomal locus


if __name__ == "__main__":
    print("run via pytest (uses capsys fixture)")
