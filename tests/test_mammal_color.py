"""Offline tests for the shared mammalian-colour engine (dna_decode/pigment/mammal_color.py) + the 5
`dna-decode {rabbit,mouse,cattle,pig,sheep}color` CLIs.

Synthetic allele calls only — no data / no network. Pins each organism's OMIA catalog reference-integrity
(the pinned epistasis anchors), canonical colours, and CLI dispatch. Runnable via pytest OR standalone.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pytest  # noqa: E402

import dna_decode.cli as uni  # noqa: E402
from dna_decode.pigment.mammal_color import (  # noqa: E402
    MAMMAL_CATALOGS,
    MammalInputError,
    call_mammal_color,
    reference_integrity_ok,
)
import dna_decode.pigment.mammal_color_cli as mcc  # noqa: E402


def test_all_catalogs_reference_integrity():
    for name, cat in MAMMAL_CATALOGS.items():
        assert reference_integrity_ok(cat) is True, name


def test_organisms_present():
    assert set(MAMMAL_CATALOGS) == {"rabbit", "mouse", "cattle", "pig", "sheep", "goat", "alpaca"}


def test_goat_asip():
    assert call_mammal_color(MAMMAL_CATALOGS["goat"], {"A": "AWt/a"}).pattern == "white/tan"
    assert call_mammal_color(MAMMAL_CATALOGS["goat"], {"A": "a/a", "B": "b/b"}).base_eumelanin == "brown/chocolate"


def test_alpaca_recessive_white():
    # the camelid twist: e/e -> white REGARDLESS of ASIP
    assert call_mammal_color(MAMMAL_CATALOGS["alpaca"], {"E": "e/e", "A": "A/A"}).is_white_masked
    assert call_mammal_color(MAMMAL_CATALOGS["alpaca"], {"E": "E/E", "A": "a/a"}).coat_color == "black"


def test_rabbit_canonical():
    c = call_mammal_color(MAMMAL_CATALOGS["rabbit"], {"A": "a/a", "B": "b/b", "D": "d/d", "E": "E/E"})
    assert c.base_eumelanin == "lilac"                     # chocolate + dilute
    assert call_mammal_color(MAMMAL_CATALOGS["rabbit"], {"C": "c/c", "E": "E/E"}).is_white_masked


def test_mouse_pink_eye():
    c = call_mammal_color(MAMMAL_CATALOGS["mouse"], {"A": "a/a", "B": "b/b", "P": "p/p", "E": "E/E"})
    assert "pink-eyed dilution" in c.dilutions and "brown/chocolate" in c.coat_color


def test_cattle_extension_and_dilution():
    assert call_mammal_color(MAMMAL_CATALOGS["cattle"], {"E": "ED/e"}).coat_color == "black"
    assert call_mammal_color(MAMMAL_CATALOGS["cattle"], {"E": "e/e"}).pigment_type == "phaeomelanin"
    dun = call_mammal_color(MAMMAL_CATALOGS["cattle"], {"E": "E+/E+", "PMEL": "Dh/Dh"})
    assert any("dilution" in d for d in dun.dilutions)


def test_pig_kit_epistatic():
    assert call_mammal_color(MAMMAL_CATALOGS["pig"], {"KIT": "I/i+", "E": "ED/ED"}).is_white_masked
    assert call_mammal_color(MAMMAL_CATALOGS["pig"], {"KIT": "i+/i+", "E": "e/e"}).pigment_type == "phaeomelanin"


def test_sheep_asip_extension_interplay():
    # ED dominant black OVERRIDES ASIP dominant-white
    assert call_mammal_color(MAMMAL_CATALOGS["sheep"], {"A": "AWt/a", "E": "ED/E+"}).coat_color == "black"
    # A^Wt white when not ED
    assert call_mammal_color(MAMMAL_CATALOGS["sheep"], {"A": "AWt/a", "E": "E+/E+"}).pattern == "white/tan"
    # a/a recessive black
    assert call_mammal_color(MAMMAL_CATALOGS["sheep"], {"A": "a/a", "E": "E+/E+"}).coat_color == "black"


def test_bad_input_raises():
    with pytest.raises(MammalInputError):
        call_mammal_color(MAMMAL_CATALOGS["rabbit"], {"E": "Z/Z"})          # unknown allele
    with pytest.raises(MammalInputError):
        call_mammal_color(MAMMAL_CATALOGS["rabbit"], {"NOSUCH": "A/a"})     # unknown locus


def test_cli_json_and_dispatch(capsys):
    rc = mcc.rabbit_main(["--loci", "A=A/A,C=C/C,E=E/E", "--json"])
    assert rc == 0
    d = json.loads(capsys.readouterr().out)
    assert d["organism"] == "Oryctolagus_cuniculus" and d["pattern"] == "agouti"
    # dispatch through the unified CLI for each organism
    for trait, fn_name in [("rabbitcolor", "rabbit_main"), ("mousecolor", "mouse_main"),
                           ("cattlecolor", "cattle_main"), ("pigcolor", "pig_main"), ("sheepcolor", "sheep_main"),
                           ("goatcolor", "goat_main"), ("alpacacolor", "alpaca_main")]:
        captured = {}
        orig = getattr(mcc, fn_name)
        setattr(mcc, fn_name, lambda argv, _c=captured: (_c.update(argv=argv) or 0))
        try:
            assert uni.main([trait, "--loci", "E=E/E"]) == 0
        finally:
            setattr(mcc, fn_name, orig)
        assert captured["argv"] == ["--loci", "E=E/E"]


def test_cli_bad_input_returns_2():
    assert mcc.cattle_main(["--loci", "E=X/X"]) == 2


if __name__ == "__main__":
    print("run via pytest (uses capsys fixture)")
