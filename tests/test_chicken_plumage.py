"""Offline tests for the chicken plumage cell (dna_decode/pigment/chicken_plumage.py) + `dna-decode plumage`.

Synthetic allele calls only — no data / no network. Pins the OMIA-curated epistatic rule, the epistasis
anchors a naive rule mis-calls (Extension canvas; Z-linked REVERSED-hemizygous barring; white masks), the
sex-dependent Z-linked loci, dilutions, silver, abstention, and CLI dispatch. Runnable via pytest OR standalone.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pytest  # noqa: E402

import dna_decode.cli as uni  # noqa: E402
from dna_decode.pigment.chicken_plumage import (  # noqa: E402
    ChickenInputError,
    call_chicken_plumage,
    reference_integrity_ok,
)
from dna_decode.pigment.chicken_plumage_cli import main as plumage_main  # noqa: E402


def test_reference_integrity():
    assert reference_integrity_ok() is True


def test_extension_canvas():
    assert call_chicken_plumage({"E": "E/E"}).eumelanin_canvas == "extended_black"
    assert call_chicken_plumage({"E": "EWh/EWh"}).eumelanin_canvas == "wheaten"
    assert call_chicken_plumage({"E": "e+/e+"}).eumelanin_canvas == "wild_partridge"


def test_anchor_extension_is_canvas():
    # barring on a wheaten canvas barely shows -> note fires (Extension is the canvas)
    c = call_chicken_plumage({"E": "EWh/EWh", "B": "B/b+"})
    assert any("naive" in n for n in c.notes)


def test_anchor_zlinked_reversed_hemizygosity():
    # a FEMALE is ZW hemizygous -> 1 allele; a MALE ZZ -> 2 alleles (reversed from mammals)
    hen = call_chicken_plumage({"E": "E/E", "B": "B"})
    cock = call_chicken_plumage({"E": "E/E", "B": "B/B"})
    assert "female" in hen.sex_basis and hen.barred
    assert "male" in cock.sex_basis and cock.barred


def test_anchor_white_masks():
    assert call_chicken_plumage({"E": "E/E", "I": "I/i+"}).is_white_masked
    assert call_chicken_plumage({"E": "E/E", "C": "c/c"}).white_type == "recessive_white"


def test_dilutions_and_silver():
    assert call_chicken_plumage({"E": "E/E", "BL": "Bl/bl+"}).dilution == "blue"
    assert call_chicken_plumage({"E": "E/E", "BL": "Bl/Bl"}).dilution == "splash"
    assert call_chicken_plumage({"E": "E/E", "LAV": "lav/lav"}).dilution == "lavender"
    assert call_chicken_plumage({"E": "ER/ER", "S": "S/S"}).silver is True


def test_explicit_sex_flag():
    c = call_chicken_plumage({"E": "E/E", "B": "B/b+"}, sex="male")
    assert "male" in c.sex_basis and c.barred


def test_unmodelled_locus_abstains_and_raises():
    c = call_chicken_plumage({"E": "E/E"}, present_loci=["CO", "MO"])
    assert any("CO:" in a for a in c.abstains_on) and any("MO:" in a for a in c.abstains_on)
    with pytest.raises(ChickenInputError):
        call_chicken_plumage({"E": "E/E", "CO": "x/y"})


def test_cli_json_and_dispatch(capsys):
    rc = plumage_main(["--loci", "E=E/E,B=B", "--sex", "female", "--json"])
    assert rc == 0
    d = json.loads(capsys.readouterr().out)
    assert d["organism"] == "Gallus_gallus" and d["barred"] is True and "female" in d["sex_basis"]
    captured = {}

    def fake(argv):
        captured["argv"] = argv
        return 0

    import dna_decode.pigment.chicken_plumage_cli as pc
    orig = pc.main
    pc.main = fake
    try:
        assert uni.main(["plumage", "--loci", "E=E/E"]) == 0
    finally:
        pc.main = orig
    assert captured["argv"] == ["--loci", "E=E/E"]


def test_cli_bad_input_returns_2():
    assert plumage_main(["--loci", "E=Q/Q"]) == 2              # unknown allele
    assert plumage_main(["--loci", "I=I/i+/i+"]) == 2          # 3 alleles at an autosomal locus


if __name__ == "__main__":
    print("run via pytest (uses capsys fixture)")
