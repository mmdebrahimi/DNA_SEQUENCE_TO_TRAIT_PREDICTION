"""The 0.911 figure must never read as CROSS-ORGANISM transfer evidence.

`wiki/essentiality_e3_human_2026-07-28.json` records 0.9107 as `learned_gbm_auroc` under
`"cv": "5-fold stratified"` on `"organism": "Homo sapiens"` -- a WITHIN-human CV whose own
baseline on that same split is `conserved_core_auroc` 0.5718. The CROSS-ORGANISM transfer
number is a different quantity: 0.5805, the conserved-core decoder applied UNCHANGED to human.

Three shipped surfaces used to say the learned complement "lifts it: E. coli 0.795 / human
0.911" in a sentence whose immediately preceding clause was the 0.580 TRANSFER number -- so a
reader took 0.911 as the transfer figure improved. `dna_decode/eval/regime.py` already warns
in prose that misreading it that way is "a live trap", and the trap was sitting in `--help`,
the most-read surface there is.

THE COINCIDENCE IS WHY IT READ PLAUSIBLY, and it is the reason a value-only check cannot catch
this class: the transfer number (0.5805) and the within-human baseline (0.5718) are two
DIFFERENT quantities whose values nearly coincide. That is the documented
`contract_number_audit` limitation -- a cited number can match the right VALUE in the wrong
QUANTITY -- so this guard checks the PROSE RELATIONSHIP, not the numerals.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

# The three shipped surfaces that carried the over-claim.
SURFACES = (
    "dna_decode/cli.py",
    "dna_decode/essentiality/cli.py",
    "dna_decode/data/cell_registry.py",
)

# Verbatim pre-fix strings, kept ONLY to prove this guard is non-vacuous. If the guard passed
# on these it would be testing nothing.
PRE_FIX_STRINGS = (
    "Cross-organism transfer to human (BAGEL CEG2/NEG) AUROC 0.580. The learned E3 complement "
    "(aa-composition+length+core, 5-fold CV) lifts it: E. coli 0.795 / human 0.911. NOT clinical.",
    "tail — the learned E3 complement lifts that: E. coli 0.795 / human 0.911). NOT a clinical tool.",
    "human (BAGEL CEG2/NEG) AUROC 0.580 (wiki/essentiality_report_card.json). The learned E3 "
    "complement lifts it (E. coli 0.795 wiki/essentiality_e3_learned_2026-07-28.json / "
    "human 0.911 wiki/essentiality_e3_human_2026-07-28.json)",
)

_WITHIN_MARKER = re.compile(r"within-organism|within-human", re.I)
_BARE_LIFT = re.compile(r"\blifts?\s+(it|that)\b", re.I)


def _segments_mentioning(text: str, needle: str) -> list[str]:
    """Paragraph-ish segments of `text` that mention `needle`.

    Splits on blank lines so a docstring paragraph or a contract string block is one segment --
    the unit a human actually reads, which is the unit the claim is made in.
    """
    return [seg for seg in re.split(r"\n\s*\n", text) if needle in seg]


def _violations(text: str) -> list[str]:
    """Every way a segment can present 0.911 as cross-organism transfer evidence."""
    bad: list[str] = []
    for seg in _segments_mentioning(text, "0.911"):
        if not _WITHIN_MARKER.search(seg):
            bad.append("0.911 appears with no within-organism marker")
        if _BARE_LIFT.search(seg) and "not lift" not in seg.lower():
            bad.append("a bare 'lifts it/that' attaches the learned CV to the transfer number")
    return bad


@pytest.mark.parametrize("rel", SURFACES)
def test_surface_does_not_present_0911_as_transfer_evidence(rel):
    text = (ROOT / rel).read_text(encoding="utf-8")
    if "0.911" not in text:
        pytest.skip(f"{rel} no longer cites 0.911")
    assert _violations(text) == [], f"{rel}: {_violations(text)}"


@pytest.mark.parametrize("pre", PRE_FIX_STRINGS, ids=("unified_cli", "essentiality_cli", "registry"))
def test_the_guard_is_non_vacuous_it_rejects_the_pre_fix_strings(pre):
    """Proven non-vacuous: every pre-fix string MUST be flagged.

    Without this, a guard that silently matched nothing would pass forever and the over-claim
    could be reintroduced unnoticed.
    """
    assert _violations(pre), f"guard failed to flag a known over-claim: {pre[:70]!r}"


def test_artifact_still_records_0911_as_within_human_cv():
    """The premise of the fix, pinned against the artifact rather than asserted from memory.

    If a future artifact rewrite made 0.911 a genuine transfer number, this guard's whole
    rationale would be void -- so the rationale is checked, not trusted.
    """
    import json

    a = json.loads((ROOT / "wiki/essentiality_e3_human_2026-07-28.json").read_text(encoding="utf-8"))
    assert a["organism"] == "Homo sapiens"
    assert "fold" in a["cv"], a["cv"]                      # a CV, not a cross-organism holdout
    assert round(a["learned_gbm_auroc"], 3) == 0.911
    # its OWN baseline on that same within-human split -- NOT the 0.5805 transfer number
    assert round(a["conserved_core_auroc"], 3) == 0.572


def test_regime_module_still_carries_the_live_trap_warning():
    """`regime.py` is where the trap is documented; if that warning is ever deleted the
    surfaces lose their cross-reference and this guard becomes the only record."""
    t = (ROOT / "dna_decode/eval/regime.py").read_text(encoding="utf-8")
    assert "0.911" in t and "trap" in t.lower()
