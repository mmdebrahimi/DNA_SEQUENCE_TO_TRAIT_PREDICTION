"""Offline tests for the dog coat-colour decoder (dna_decode/pigment/dog_coat.py) — pure epistasis logic.

Pins: the reference-integrity biology anchors (known breed genotypes -> colours), the E-locus EPISTASIS
anchor a naive rule mis-calls, eumelanin colour B x D, distribution K x A, pattern-locus abstention, and
input-validation errors. No network / no Docker / no data. Runnable via pytest OR standalone.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.pigment.dog_coat import (  # noqa: E402
    CoatInputError,
    call_coat_color,
    parse_genotype,
    reference_integrity_ok,
)


def test_reference_integrity_biology_contract():
    assert reference_integrity_ok() is True


def test_yellow_lab_recessive_red():
    r = call_coat_color({"E": "e/e", "B": "B/B", "D": "D/D"})
    assert r.coat_color == "red/yellow"
    assert r.pigment_type == "phaeomelanin"
    assert r.eumelanin_color is None


def test_e_locus_epistasis_anchor():
    # e/e is red/yellow EVEN with dominant-black K^B and brown b/b — the case a naive rule gets wrong.
    r = call_coat_color({"E": "e/e", "K": "KB/KB", "B": "b/b", "D": "D/D"})
    assert r.coat_color == "red/yellow"
    assert r.pigment_type == "phaeomelanin"
    assert any("EPISTATIC" in n or "MIS-CALL" in n for n in r.notes)
    # b/b still lightens the nose even on a red coat — surfaced as a note, not a coat call
    assert any("nose" in n.lower() for n in r.notes)


def test_eumelanin_color_bxd():
    assert call_coat_color({"E": "E/E", "K": "KB/KB", "B": "B/B", "D": "D/D"}).eumelanin_color == "black"
    assert call_coat_color({"E": "E/E", "K": "KB/KB", "B": "b/b", "D": "D/D"}).eumelanin_color == "brown/liver"
    assert call_coat_color({"E": "E/E", "K": "KB/KB", "B": "B/B", "D": "d/d"}).eumelanin_color == "blue/grey"
    assert call_coat_color({"E": "E/E", "K": "KB/KB", "B": "b/b", "D": "d/d"}).eumelanin_color == "isabella/lilac"


def test_b_family_alleles_collapse_recessive():
    # bs (p.Gln331Ter) is a b-family allele: bs/bs -> brown same as b/b
    assert call_coat_color({"E": "E/E", "K": "KB/KB", "B": "bs/bs", "D": "D/D"}).eumelanin_color == "brown/liver"
    # heterozygous B/bs -> still black (B dominant)
    assert call_coat_color({"E": "E/E", "K": "KB/KB", "B": "B/bs", "D": "D/D"}).eumelanin_color == "black"


def test_distribution_kxa():
    # ky/ky expresses A: tan-points, sable, agouti, recessive-black
    assert call_coat_color({"E": "E/E", "K": "ky/ky", "A": "at/at", "B": "B/B", "D": "D/D"}).distribution == "tan_points"
    assert call_coat_color({"E": "E/E", "K": "ky/ky", "A": "Ay/at", "B": "B/B", "D": "D/D"}).distribution == "sable"
    assert call_coat_color({"E": "E/E", "K": "ky/ky", "A": "aw/aw", "B": "B/B", "D": "D/D"}).distribution == "agouti"
    assert call_coat_color({"E": "E/E", "K": "ky/ky", "A": "a/a", "B": "B/B", "D": "D/D"}).distribution == "recessive_black"
    # K^B masks A -> solid regardless of A
    assert call_coat_color({"E": "E/E", "K": "KB/ky", "A": "at/at", "B": "B/B", "D": "D/D"}).distribution == "solid"


def test_agouti_unknown_when_A_absent_under_kyky():
    r = call_coat_color({"E": "E/E", "K": "ky/ky", "B": "B/B", "D": "D/D"})
    assert r.distribution == "agouti_unknown"
    assert any(ax.startswith("A:") for ax in r.abstains_on)
    assert r.confidence == "medium"


def test_pattern_locus_abstains():
    r = call_coat_color({"E": "E/E", "K": "KB/KB", "B": "B/B", "D": "D/D"}, present_loci=["M", "S"])
    axes = " ".join(r.abstains_on)
    assert "M:" in axes and "S:" in axes
    assert r.coat_color == "solid black"   # still a colour call; the pattern axes are withheld separately


def test_pattern_locus_in_loci_dict_raises_with_guidance():
    try:
        call_coat_color({"E": "E/E", "M": "M/m"})
    except CoatInputError as e:
        assert "present_loci" in str(e)
        return
    raise AssertionError("expected CoatInputError steering merle to present_loci")


def test_missing_E_locus_raises():
    try:
        call_coat_color({"B": "b/b", "D": "d/d"})
    except CoatInputError as e:
        assert "E" in str(e) and "pigment-type" in str(e)
        return
    raise AssertionError("expected CoatInputError when E locus absent")


def test_missing_color_locus_caps_confidence():
    # no B/D given -> wild-type assumed, confidence capped medium + a note
    r = call_coat_color({"E": "E/E", "K": "KB/KB"})
    assert r.eumelanin_color == "black"
    assert r.confidence == "medium"
    assert any("wild-type" in n for n in r.notes)


def test_parse_and_unknown_allele():
    assert parse_genotype("E", "e/E") == ("e", "E")
    try:
        parse_genotype("B", "B/z")
    except CoatInputError:
        return
    raise AssertionError("expected CoatInputError on unknown allele")


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn(); print(f"PASS {fn.__name__}")
    print(f"\n{len(fns)} passed")
