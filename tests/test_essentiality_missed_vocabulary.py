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
    import re

    from dna_decode.essentiality.core_decoder import _CORE

    m = re.search(r"ROBUST_II\s*=\s*\{([^}]*)\}", SRC)
    assert m, "ROBUST_II literal not found -- the floor must stay auditable in source"
    words = {w.strip().strip("\"'") for w in m.group(1).split(",") if w.strip()}
    assert len(words) >= 4, f"floor suspiciously small ({words}) -- would pass vacuously"

    pattern_text = " ".join(p.pattern.lower() for _w, p in _CORE)
    for w in words:
        assert w in pattern_text, f"{w!r} is not in the decoder's own patterns -- that is asserted biology"


def test_robust_floors_do_not_depend_on_the_reach_definition():
    """An ABSENT word is absent under any reach definition, and the (ii) floor words are unambiguous.
    If a future edit made a floor depend on `eco_caught_words`, the floor stops being robust."""
    seg = SRC[SRC.index("ROBUST_II"):SRC.index("n_missed = len(missed)") + 400] \
        if "n_missed = len(missed)" in SRC[SRC.index("ROBUST_II"):] else SRC[SRC.index("ROBUST_II"):]
    assert "eco_caught_words" not in seg.split("robust_iii")[0], \
        "the mechanism-(ii) floor must not be computed from the loose reach set"


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
