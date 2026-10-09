"""CLAUDE.md is auto-loaded EVERY session and EVERY compaction. Its size is a running cost, not a detail.

WHY THIS EXISTS (measured 2026-09-29). CLAUDE.md grew 41k -> 116k -> 277k chars since June (6.7x) and now
costs **~68,665 tokens on every single session load** -- ~80k once the parent CLAUDE.md and the
cross-project MEMORY.md are counted, which is ~40% of a 200k window gone before any work starts. That is
self-defeating: the document written to PREVENT rework is the main consumer of the window whose exhaustion
CAUSES the rework.

WHAT THIS TEST IS, AND IS NOT. It is a **growth ceiling**, not a diet. It does NOT claim the file has been
shrunk -- it pins today's size as a CEILING so the next addition is a DELIBERATE choice (raise the number,
on the record) instead of silent drift. `scripts/claude_md_weight.py` already reports that 39 bullets /
26,718 words (72% of the file) are PROVABLY stored elsewhere in wiki/, so the diet is available; this
guard just stops the trend while it is pending.

Pinned in WORDS, not tokens: words are deterministic, a token estimate depends on a tokenizer.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
CLAUDE_MD = ROOT / "CLAUDE.md"

# Measured 2026-09-29: 34,758 words after the first proof-of-format compression (was 36,861 / ~68,665 tok).
# LOWERED with that gain -- leaving a ceiling sized for the pre-diet file would silently stop guarding.
# Headroom is deliberately THIN (~2%). A thick allowance is how 41k became 277k.
#
# RAISED 2026-10-08: 35_400 -> 35_460, on the record exactly as this guard's docstring prescribes.
# Cost: a 60-word pointer recording that GLM families G-A and G-C closed as valid negatives, which exists
# to stop a future session re-proposing a conv encoder or re-running self-distillation (the paragraph it
# attaches to says self-distillation is "forbidden until a learned representation lands", and leaving
# that unqualified now reads as an open invitation). Derivations are NOT in CLAUDE.md -- they are in two
# wiki memos the pointer names.
#
# THE SIGNAL WORTH ACTING ON: the file arrived at this edit with **12 words of headroom**, so the ceiling
# is now binding on ordinary findings rather than on drift, and every future addition will face the same
# choice. `scripts/claude_md_weight.py` reports **43 bullets / 24,892 words (70% of the file) PROVABLY
# stored elsewhere in wiki/**, so the diet is available and is a USER AUTHORITY call (NEXT.md fork 5),
# not something to take unilaterally. Prefer compressing a candidate bullet over raising this number
# again.
WORD_CEILING = 35_460
MEASURED_AT_WRITE = 34_758


def _weight_module():
    spec = importlib.util.spec_from_file_location("cmw", ROOT / "scripts" / "claude_md_weight.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _words(text: str) -> int:
    return len(text.split())


def test_claude_md_stays_under_the_growth_ceiling():
    """The ceiling converts silent growth into a deliberate edit of THIS number."""
    words = _words(CLAUDE_MD.read_text(encoding="utf-8"))
    assert words <= WORD_CEILING, (
        f"CLAUDE.md is {words:,} words, over the {WORD_CEILING:,} ceiling. This file loads on EVERY "
        f"session and EVERY compaction, so growth is a recurring context cost.\n"
        f"Before raising the ceiling, run `PYTHONIOENCODING=utf-8 uv run python "
        f"scripts/claude_md_weight.py` -- it lists the bullets that are PROVABLY stored elsewhere in "
        f"wiki/ and can be compressed to a pointer instead. Raising this number is a real decision; "
        f"make it on purpose.")


def test_the_ceiling_is_not_vacuously_slack():
    """A ceiling far above the real size cannot detect anything. Pin that the headroom stays thin.

    If a legitimate diet lands, LOWER both constants -- do not leave a ceiling sized for the old file,
    or the guard silently stops guarding.
    """
    words = _words(CLAUDE_MD.read_text(encoding="utf-8"))
    slack = WORD_CEILING - words
    assert slack <= 4_000, (
        f"ceiling has {slack:,} words of slack over the live file ({words:,}). That is too loose to "
        f"detect drift -- lower WORD_CEILING to roughly the current size.")


def test_every_cited_artifact_in_claude_md_actually_resolves():
    """A pointer into wiki/ is the compression mechanism. A pointer that does not resolve is worse than
    the prose it replaced -- it reads as provenance and delivers nothing."""
    mod = _weight_module()
    text = CLAUDE_MD.read_text(encoding="utf-8")
    missing: list[str] = []
    for bullet in mod.bullets(text):
        info = mod.analyse(bullet)
        for cited in info.get("cited", ()) or ():
            if not mod.resolves(cited):
                missing.append(cited)
    assert not missing, f"{len(missing)} cited path(s) in CLAUDE.md do not resolve: {sorted(set(missing))[:8]}"


def test_the_bullets_that_are_the_ONLY_copy_are_identified_and_not_silently_compressible():
    """The safety rail on any future diet: some bullets cite nothing, so CLAUDE.md is the only place that
    knowledge exists. Those must never be compressed to a pointer. This test asserts the weight tool can
    still tell them apart -- if it ever reports zero, the diet has lost its guard rail."""
    mod = _weight_module()
    text = CLAUDE_MD.read_text(encoding="utf-8")
    long_no_store = [
        b for b in mod.bullets(text)
        if (info := mod.analyse(b))["words"] >= mod.LONG_BULLET_WORDS and not info.get("stored")
    ]
    assert long_no_store, (
        "the weight tool reports NO long-uncited bullets. Either the file changed shape or `analyse` "
        "broke; either way a diet driven by it would treat single-copy knowledge as compressible.")


@pytest.mark.parametrize("path", ["NEXT.md"])
def test_the_transient_store_stays_transient(path: str):
    """NEXT.md is documented as transient ('prune it, don't grow it'). It is not auto-loaded, so this is a
    soft bound -- but the same growth mechanism applies, and it is the file most likely to accrete."""
    words = _words((ROOT / path).read_text(encoding="utf-8"))
    assert words <= 9_000, (
        f"{path} is {words:,} words. It is declared transient; prune resolved items rather than growing it.")
