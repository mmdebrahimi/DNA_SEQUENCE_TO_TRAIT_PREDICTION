"""Guards for the missed-vocabulary probe — the measurement that RETRACTED its own verdict.

The probe's first version used one bar ("any content word shared with a caught E. coli gene"), reported
frac_ii 0.5548 and a verdict of PHRASING_GAP_DOMINATES_LADDER_IS_WRONG_INSTRUMENT. That bar is inflated
>2x by generic modifiers (binding / small / beta / alpha / cell / repeat) -- sharing "alpha" with E. coli
says nothing about whether the catalogue could see a gene -- and the principled bar lands in a different
band. These tests pin the two ways that failure could come back.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "essentiality_missed_vocabulary.py"
SRC = SCRIPT.read_text(encoding="utf-8")


def test_a_single_bar_cannot_produce_the_verdict():
    """THE retracted failure. Three bars must be computed and the verdict must be able to refuse."""
    for bar in ("loose_any_shared_word", "strict_regex_vocabulary", "strict_minus_regex_fragments"):
        assert bar in SRC, bar
    assert "INDETERMINATE_BAR_SENSITIVE_VOCABULARY_CANNOT_ADJUDICATE" in SRC
    # and the refusal must be driven by band DISAGREEMENT, not hardcoded
    assert "len(set(bands.values())) > 1" in SRC


def test_no_curated_biology_is_written_into_the_probe():
    """The fabrication rail. Mechanism (ii) is adjudicated by corpus comparison, so the ONLY biological
    word list in the probe is the robust-(ii) floor -- and every word in it must already appear in the
    decoder's OWN regex patterns rather than being invented here.

    Parsed from SOURCE, not imported: importing the script would need sys.path surgery, and an import
    that fails would make this test error rather than check. The non-emptiness assert matters as much as
    the membership one -- an empty set would satisfy a bare `for` loop and pass vacuously forever.
    """
    from dna_decode.essentiality.core_decoder import _CORE
    from dna_decode.essentiality.phrasing_floor import ROBUST_II

    # PROMOTED 2026-10-06: the floor moved to dna_decode/essentiality/phrasing_floor.py, where an
    # IMPORT-TIME gate raises if any word is absent from _CORE. This test still checks the property
    # directly rather than trusting that gate -- a gate and its check should not be the same code.
    words = set(ROBUST_II)
    assert len(words) >= 4, f"floor suspiciously small ({words}) -- would pass vacuously"
    assert "ROBUST_II" in SRC, "the script must still reference the floor it reports"

    pattern_text = " ".join(p.pattern.lower() for _w, p in _CORE)
    for w in words:
        assert w in pattern_text, f"{w!r} is not in the decoder's own patterns -- that is asserted biology"


def test_robust_floors_do_not_depend_on_the_reach_definition():
    """An ABSENT word is absent under any reach definition, and the (ii) floor words are unambiguous.
    If a future edit made the floor depend on `eco_caught_words`, the floor stops being robust.

    REPOINTED 2026-10-06: the floor moved to dna_decode/essentiality/phrasing_floor.py (transfer-ladder
    Step 4). The old version sliced THIS script's source from the first `ROBUST_II` occurrence, which is
    now the import line at the top, so the slice swallowed the whole file including the loose-reach
    computation. Asserting against the module that actually computes the floor is both the correct scope
    and a stronger check than a source slice.
    """
    floor_src = (ROOT / "dna_decode" / "essentiality" / "phrasing_floor.py").read_text(encoding="utf-8")
    assert "eco_caught_words" not in floor_src, \
        "the mechanism-(ii) floor must not be computed from the loose reach set"
    assert "ROBUST_II" in floor_src and "def phrasing_floor_genes" in floor_src


def test_the_artifact_records_the_self_correction_and_is_parseable():
    """The retraction has to survive in the artifact, not just in a commit message."""
    hits = sorted(ROOT.glob("wiki/essentiality_missed_vocabulary_*.json"))
    if not hits:
        pytest.skip("probe artifact not emitted on this host")
    art = json.loads(hits[-1].read_text(encoding="utf-8"))   # parseable, per the FBA-artifact lesson
    assert art["record"] == "essentiality-missed-vocabulary-v1"
    assert "self_correction" in art and "0.5548" in art["self_correction"]
    assert set(art["bar_sensitivity"]) == {
        "loose_any_shared_word", "strict_regex_vocabulary", "strict_minus_regex_fragments"}
    # the loose bar must still be visibly the largest -- that IS the retracted inflation, kept on record
    b = art["bar_sensitivity"]
    assert b["loose_any_shared_word"] > b["strict_regex_vocabulary"]
    assert art["honest_limits"], "limits must not be emptied"
    # a word is not a function -- the caveat that bounds every number here
    assert any("WORD is not a function" in x for x in art["honest_limits"])


def test_probe_is_read_only_without_emit():
    """--emit is the ONLY write path; a bare run must not touch wiki/."""
    assert SRC.count("write_text") == 1, "exactly one write path expected"
    i = SRC.index("write_text")
    assert "if a.emit" in SRC[max(0, i - 400):i], "the single write must be behind --emit"
