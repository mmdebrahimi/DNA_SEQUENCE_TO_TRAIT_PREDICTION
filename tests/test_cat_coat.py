"""Offline tests for the cat coat-colour cell (dna_decode/pigment/cat_coat.py) + `dna-decode catcolor`.

Synthetic allele calls only — no data / no network. Pins the OMIA-curated epistatic rule, the THREE epistasis
anchors a naive rule mis-calls (dominant-white; X-linked tortoiseshell; orange-over-brown), the sex-dependent
X-linked O locus, dilute, colorpoint, calico, abstention, and CLI dispatch. Runnable via pytest OR standalone.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pytest  # noqa: E402

import dna_decode.cli as uni  # noqa: E402
from dna_decode.pigment.cat_coat import (  # noqa: E402
    CatInputError,
    call_cat_coat,
    reference_integrity_ok,
)
from dna_decode.pigment.cat_coat_cli import main as cat_main  # noqa: E402


def test_reference_integrity():
    assert reference_integrity_ok() is True


def test_anchor_dominant_white_epistatic():
    c = call_cat_coat({"W": "W/w", "O": "o/o", "B": "b/b"})
    assert c.is_epistatic_white and c.coat_color.startswith("white (dominant")


def test_anchor_tortoiseshell_xlinked():
    c = call_cat_coat({"O": "O/o", "B": "B/B", "D": "D/D"})
    assert c.is_tortoiseshell and c.base_color.startswith("tortoiseshell")
    assert "female" in c.sex_basis and any("X-inactivation" in n for n in c.notes)


def test_anchor_calico_is_tortie_plus_white():
    c = call_cat_coat({"O": "O/o", "W": "ws/w", "B": "B/B"})
    assert c.white_pattern == "calico" and c.coat_color.startswith("calico")


def test_anchor_orange_epistatic_over_brown():
    # a b/b (chocolate-genotype) ORANGE male is RED, not chocolate
    c = call_cat_coat({"O": "O", "B": "b/b"})
    assert c.base_color == "red" and "male" in c.sex_basis


def test_base_eumelanin_colors():
    assert call_cat_coat({"O": "o", "A": "a/a", "B": "B/B", "D": "D/D"}).base_color == "black"
    assert call_cat_coat({"O": "o", "A": "a/a", "B": "B/B", "D": "d/d"}).base_color == "blue"     # dilute
    assert call_cat_coat({"O": "o", "A": "a/a", "B": "b/b", "D": "D/D"}).base_color == "chocolate"
    assert call_cat_coat({"O": "o", "A": "a/a", "B": "bl/bl", "D": "D/D"}).base_color == "cinnamon"


def test_sex_inference_from_o_zygosity():
    assert "male" in call_cat_coat({"O": "O"}).sex_basis          # 1 allele -> male
    assert "female" in call_cat_coat({"O": "O/O"}).sex_basis      # 2 alleles -> female


def test_explicit_sex_flag():
    c = call_cat_coat({"O": "O/o"}, sex="female")
    assert "female" in c.sex_basis and c.is_tortoiseshell


def test_tabby_and_colorpoint():
    assert call_cat_coat({"O": "o", "A": "A/a", "B": "B/B"}).tabby is True
    assert call_cat_coat({"O": "o", "A": "a/a", "B": "B/B"}).tabby is False
    assert call_cat_coat({"O": "o", "C": "cs/cs", "B": "B/B"}).colorpoint == "siamese_points"
    assert call_cat_coat({"O": "o", "C": "cs/cb", "B": "B/B"}).colorpoint == "mink"


def test_albino_masks():
    c = call_cat_coat({"O": "o", "C": "c/c"})
    assert c.is_epistatic_white and "albino" in c.coat_color


def test_unmodelled_locus_abstains_and_raises():
    c = call_cat_coat({"O": "o", "B": "B/B"}, present_loci=["TA", "I"])
    assert any("TA:" in a for a in c.abstains_on) and any("I:" in a for a in c.abstains_on)
    with pytest.raises(CatInputError):
        call_cat_coat({"O": "o", "TA": "x/y"})


def test_cli_json_and_dispatch(capsys):
    rc = cat_main(["--loci", "O=O/o,W=ws/w,B=B/B", "--json"])
    assert rc == 0
    d = json.loads(capsys.readouterr().out)
    assert d["organism"] == "Felis_catus" and d["white_pattern"] == "calico"
    captured = {}

    def fake(argv):
        captured["argv"] = argv
        return 0

    import dna_decode.pigment.cat_coat_cli as cc
    orig = cc.main
    cc.main = fake
    try:
        assert uni.main(["catcolor", "--loci", "O=o,B=B/B"]) == 0
    finally:
        cc.main = orig
    assert captured["argv"] == ["--loci", "O=o,B=B/B"]


def test_cli_bad_input_returns_2():
    assert cat_main(["--loci", "O=X"]) == 2                 # unknown allele
    assert cat_main(["--loci", "B=B/B/b"]) == 2             # 3 alleles at an autosomal locus


if __name__ == "__main__":
    print("run via pytest (uses capsys fixture)")
