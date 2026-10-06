"""Guards for the phrasing-floor subtraction — a FLOOR, provenance-gated, and never clamped.

Context that sets the bar: `wiki/essentiality_missed_vocabulary_2026-10-05` measured that vocabulary
overlap CANNOT adjudicate the middle of the phrasing-vs-host-specific question (three defensible bars,
two verdict bands; the loosest bar produced a confident verdict that was retracted pre-publication). So
what ships is only the part that survives every bar, and a broad synonym vocabulary deliberately does not
exist.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from dna_decode.essentiality import phrasing_floor as pf
from dna_decode.essentiality.core_decoder import _CORE, score_gene

ROOT = Path(__file__).resolve().parents[1]
CACHE = Path("D:/dna_decode_cache/essentiality")
_HAVE = all((CACHE / f).exists() for f in ("CEGv2.txt", "NEGv1.txt", "Homo_sapiens.gene_info.gz"))


# --------------------------------------------------------------------------------------------------
# provenance gate
# --------------------------------------------------------------------------------------------------
def test_every_floor_word_is_in_the_decoders_own_patterns():
    txt = " ".join(p.pattern.lower() for _w, p in _CORE)
    for w in pf.ROBUST_II:
        assert w in txt, f"{w!r} absent from _CORE -- that would be asserted biology"


def test_the_floor_is_not_empty_so_the_gate_cannot_pass_vacuously():
    """A gate looping over an empty set passes forever. This is the shape I shipped once already."""
    assert len(pf.ROBUST_II) >= 4


def test_a_planted_non_pattern_word_RAISES_the_provenance_gate():
    """NON-VACUITY of the import-time gate, exercised without mutating module state permanently."""
    original = pf.ROBUST_II
    try:
        pf.ROBUST_II = frozenset(set(original) | {"zzzznotapattern"})
        with pytest.raises(pf.PhrasingFloorProvenanceError, match="zzzznotapattern"):
            pf._enforce_provenance()
    finally:
        pf.ROBUST_II = original
    pf._enforce_provenance()          # control: restored state still passes


def test_stop_words_contain_no_biology():
    """Boilerplate + English function words only. A biological stopword would silently drop a function."""
    for w in ("polymerase", "helicase", "ribosomal", "proteasome", "replication", "division"):
        assert w not in pf.STOP, w


# --------------------------------------------------------------------------------------------------
# floor selection
# --------------------------------------------------------------------------------------------------
def test_is_phrasing_floor_fires_on_function_words_only():
    assert pf.is_phrasing_floor("DNA polymerase delta subunit 1")
    assert pf.is_phrasing_floor("ATP-dependent RNA helicase")
    assert not pf.is_phrasing_floor("proteasome 26S subunit")
    assert not pf.is_phrasing_floor("small nuclear ribonucleoprotein")
    # the generic modifiers that inflated the retracted bar must NOT qualify
    for junk in ("binding partner", "small subunit", "beta chain", "alpha unit", "repeat region"):
        assert not pf.is_phrasing_floor(junk), junk


def test_phrasing_floor_genes_is_scoped_to_MISSED_genes():
    missed = [("POLD1", "DNA polymerase delta 1"), ("PSMA1", "proteasome subunit alpha 1")]
    assert pf.phrasing_floor_genes(missed) == {"POLD1"}


# --------------------------------------------------------------------------------------------------
# adjusted lift — reported honestly, never clamped
# --------------------------------------------------------------------------------------------------
def _zero(_g, _t):
    return 0.0


def test_the_adjusted_lift_RISES_and_is_not_clamped():
    """Floor genes are by construction MISSES, so removing them raises coverage(essential) and the lift.
    The plan explicitly refuses to assert adjusted <= raw: clamping would hide the very quantity this
    exists to expose."""
    ess = [("POLD1", "DNA polymerase delta 1"),        # floor, missed
           ("PSMA1", "proteasome subunit alpha"),      # not floor, missed
           ("RPL3", "ribosomal protein L3")]           # caught by the real decoder
    non = [("NEG1", "olfactory receptor"), ("NEG2", "keratin associated")]

    def score(g, t):
        return score_gene(g, t).core_score

    raw, adj, detail = pf.coverage_lift_adjusted(ess, non, score)
    assert detail["n_phrasing_floor"] == 1
    assert adj > raw, (raw, adj)          # rises, as the mechanism predicts
    # 1/3 -> 1/2 on the essential side, non-essential unchanged
    assert detail["coverage_essential_raw"] == pytest.approx(1 / 3, abs=1e-4)
    assert detail["coverage_essential_adjusted"] == pytest.approx(0.5, abs=1e-4)


def test_it_refuses_when_the_floor_would_remove_the_entire_essential_class():
    ess = [("POLD1", "DNA polymerase delta"), ("MCM2", "DNA helicase")]
    non = [("NEG1", "olfactory receptor")]
    with pytest.raises(ValueError, match="ENTIRE essential class"):
        pf.coverage_lift_adjusted(ess, non, _zero)


def test_it_refuses_an_empty_class():
    with pytest.raises(ValueError, match="non-empty"):
        pf.coverage_lift_adjusted([], [("a", "b")], _zero)


def test_the_detail_states_that_this_is_a_floor_not_an_estimate():
    ess = [("POLD1", "DNA polymerase delta"), ("RPL3", "ribosomal protein L3")]
    non = [("NEG1", "olfactory receptor")]
    _r, _a, detail = pf.coverage_lift_adjusted(ess, non, lambda g, t: score_gene(g, t).core_score)
    assert "FLOOR" in detail["floor_is_a_FLOOR_not_an_estimate"].upper()
    assert "0.25" in detail["floor_is_a_FLOOR_not_an_estimate"]


# --------------------------------------------------------------------------------------------------
# the promotion is real, and the published number reproduces
# --------------------------------------------------------------------------------------------------
def test_the_script_imports_the_vocabulary_from_the_package():
    src = (ROOT / "scripts" / "essentiality_missed_vocabulary.py").read_text(encoding="utf-8")
    assert "from dna_decode.essentiality.phrasing_floor import" in src
    assert 'ROBUST_II = {"polymerase"' not in src, "the local literal must be gone, not shadowing"
    assert "def content_words(text):" not in src, "content_words must come from the package"


@pytest.mark.skipif(not _HAVE, reason="D: essentiality cache absent")
def test_the_floor_reproduces_57_genes_on_the_real_human_missed_set():
    """The published figure, pinned against the real data rather than restated."""
    import sys

    sys.path.insert(0, str(ROOT / "scripts"))
    from essentiality_missed_vocabulary import human_rows

    ceg, _neg = human_rows()
    missed = [(g, d) for g, d in ceg if score_gene(g, d).core_score == 0]
    assert len(missed) == 566, len(missed)
    assert len(pf.phrasing_floor_genes(missed)) == 57
