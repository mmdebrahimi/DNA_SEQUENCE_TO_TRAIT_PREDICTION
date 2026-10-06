"""Subtract the catalogue's PHRASING gap from the ladder's decay, so a bug is not read as phylogeny.

THE MEASUREMENT THIS IMPLEMENTS, and the one it refuses to implement.

`wiki/essentiality_missed_vocabulary_2026-10-05.{md,json}` asked whether the decoder's 83%-miss on human
essentials is really a catalogue PHRASING gap -- a conserved function under a bacteria-specific NAME --
which would make a transfer ladder the wrong instrument. The verdict was
`INDETERMINATE_BAR_SENSITIVE_VOCABULARY_CANNOT_ADJUDICATE`: three defensible bars spanning two verdict
bands (loose 0.5548 / strict-regex 0.2580 / strict-minus-fragments 0.2509). The first, loosest bar had
produced a confident `PHRASING_GAP_DOMINATES` that was retracted before publication -- it was inflated
>2x by generic modifiers (`binding`, `small`, `beta`, `alpha`, `cell`, `repeat`).

So TWO things follow and both are load-bearing:

  * A broad per-rung functional-class synonym vocabulary is NOT built here. It was the v1 plan's Step 4,
    was flagged in review as a second adjudication system, and was then MEASURED unable to adjudicate.
    Building it anyway would be a post-hoc rescue table.
  * What DOES ship is the FLOOR: the words that survive every bar. Host-specific absence (141/566 =
    0.2491) is ~2.5x the phrasing floor (57/566 = 0.1007), so the ladder is not invalidated -- but ~10%
    of the human miss has nothing to do with phylogeny and must be subtracted before any rung is scored,
    or it is silently attributed to distance or to the domain boundary.

PROVENANCE GATE, enforced at import rather than documented: every word in `ROBUST_II` must already
appear in the decoder's OWN `_CORE` regexes. That is what makes this a subtraction of the catalogue's
own vocabulary rather than an unsourced biological claim about which human genes "really" count --
exactly the hazard `dna_decode/data/ecoff_catalog.py` exists to refuse.
"""
from __future__ import annotations

import re

from dna_decode.essentiality.core_decoder import _CORE

# ---------------------------------------------------------------------------------------------------
# Vocabulary, promoted out of scripts/essentiality_missed_vocabulary.py so the ladder does not import a
# scripts/ module. The script now imports these from here (single definition, pinned by test).
# ---------------------------------------------------------------------------------------------------
WORD = re.compile(r"[a-z][a-z-]{2,}")

# ENGLISH function words + annotation boilerplate. Deliberately contains NO biology: dropping "protein"
# or "domain" changes which words surface, so they are listed as boilerplate, not as a judgement about
# what is or is not a function.
STOP = {
    "a", "an", "the", "of", "and", "or", "in", "to", "for", "with", "by", "on", "at", "from", "as",
    "is", "are", "be", "that", "this", "it", "its", "has", "have", "not", "no",
    "protein", "proteins", "domain", "containing", "family", "member", "like", "type", "subunit",
    "putative", "probable", "uncharacterized", "hypothetical", "homolog", "associated", "related",
    "component", "complex", "factor", "chain", "unit", "group", "part", "form",
}

# THE FLOOR. Each word is (a) in the decoder's own regex vocabulary and (b) unambiguously a FUNCTION
# rather than a modifier, so it is "in catalogue reach" under ANY definition of reach -- which is what
# makes the 57-gene figure survive all three bars. Not an estimate: a floor.
ROBUST_II = frozenset({"polymerase", "helicase", "replication", "division", "topoisomerase", "primase"})


class PhrasingFloorProvenanceError(ValueError):
    """A floor word that is not in the decoder's own patterns — i.e. asserted biology."""


def _pattern_text() -> str:
    return " ".join(p.pattern.lower() for _w, p in _CORE)


def _enforce_provenance() -> None:
    txt = _pattern_text()
    missing = sorted(w for w in ROBUST_II if w not in txt)
    if missing:
        raise PhrasingFloorProvenanceError(
            "floor words absent from the decoder's own _CORE patterns: %s. A floor word that the "
            "catalogue does not itself name is an unsourced biological claim about which genes 'really' "
            "count -- add it to _CORE deliberately (which ends the transfer experiment) or drop it here."
            % missing)


_enforce_provenance()          # at IMPORT, not in a function nobody calls


def content_words(text: str) -> set[str]:
    return {w for w in WORD.findall((text or "").lower()) if w not in STOP}


def is_phrasing_floor(text: str) -> bool:
    """Does this gene's description name a function the catalogue demonstrably reaches?"""
    return bool(content_words(text) & ROBUST_II)


def phrasing_floor_genes(missed_rows) -> set[str]:
    """`missed_rows` = (gene, description) for genes the decoder scored ZERO. -> the floor's gene set.

    Scoped to MISSED genes on purpose: a gene the decoder already catches is not part of any gap.
    """
    return {g for g, d in missed_rows if is_phrasing_floor("%s %s" % (g or "", d or ""))}


def coverage_lift_adjusted(essential_rows, nonessential_rows, score_fn):
    """coverage_lift with the phrasing floor removed from the ESSENTIAL class.

    Returns (raw, adjusted, detail). Floor genes are by construction ones the decoder MISSED, so
    removing them RAISES coverage(essential) and therefore the lift. That direction is reported, NEVER
    clamped: the adjusted number answers "what is the decay net of the catalogue's own phrasing bug",
    and clamping it to <= raw would hide exactly the quantity it exists to expose.

    `score_fn(gene, text) -> float` is injected so this stays pure and testable without the real decoder.
    """
    ess = [(g, t, score_fn(g, t)) for g, t in essential_rows]
    non = [(g, t, score_fn(g, t)) for g, t in nonessential_rows]
    if not ess or not non:
        raise ValueError("coverage_lift_adjusted needs a non-empty essential AND non-essential class")

    cov_non = sum(1 for _g, _t, s in non if s != 0) / len(non)
    cov_ess = sum(1 for _g, _t, s in ess if s != 0) / len(ess)
    raw = cov_ess - cov_non

    missed = [(g, t) for g, t, s in ess if s == 0]
    floor = phrasing_floor_genes(missed)
    kept = [(g, t, s) for g, t, s in ess if g not in floor]
    if not kept:
        raise ValueError("the phrasing floor removed the ENTIRE essential class -- refusing to report "
                         "an adjusted lift computed on nothing")
    cov_ess_adj = sum(1 for _g, _t, s in kept if s != 0) / len(kept)
    adjusted = cov_ess_adj - cov_non

    detail = {
        "n_essential": len(ess), "n_essential_missed": len(missed),
        "n_phrasing_floor": len(floor),
        "phrasing_floor_fraction_of_missed": (round(len(floor) / len(missed), 4) if missed else None),
        "coverage_essential_raw": round(cov_ess, 4),
        "coverage_essential_adjusted": round(cov_ess_adj, 4),
        "coverage_nonessential": round(cov_non, 4),
        "floor_words": sorted(ROBUST_II),
        "floor_is_a_FLOOR_not_an_estimate": (
            "this is a FLOOR, not an estimate: 57 genes survived an adversarial bar on the human rung, "
            "the true mechanism-(ii) share sits between ~0.10 and ~0.25, and vocabulary overlap cannot "
            "narrow it -- so the adjustment deliberately UNDER-corrects"),
    }
    return round(raw, 4), round(adjusted, 4), detail
