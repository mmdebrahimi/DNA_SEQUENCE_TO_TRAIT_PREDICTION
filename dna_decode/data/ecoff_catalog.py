"""EUCAST epidemiological cut-off values (ECOFFs), with provenance required per entry.

WHY THIS MODULE REFUSES RATHER THAN DEFAULTS. An ECOFF is a biological reference value. Writing one
from memory is the standing fabrication hazard in this repo, and it is a particularly dangerous one
here: a wrong cut-off produces a fully self-consistent evaluation -- every isolate classifies, every
rate computes, every test passes -- that is entirely wrong, with nothing downstream able to catch it.
So a value may only enter this catalog with a `source_url` AND a `verbatim_quote`, and `ecoff_for`
raises on anything unsourced instead of returning a plausible default.

NOT FROZEN. This module is deliberately outside the prospective-lock surface
(`dna_decode/eval/prospective_lock.py` pins `mic_tiers.py` and four others). It exists so an
ECOFF-anchored arm can be evaluated WITHOUT editing the frozen clinical-breakpoint catalog; adopting
ECOFF anchoring in the shipped decoder is a separate, user-authority decision.

STATUS 2026-09-11: EVERY ENTRY IS UNSOURCED. The values were not obtainable from this environment --
EUCAST publishes no downloadable ECOFF table and no stable per-result URL, and the numbers live only
behind an interactive JavaScript database (`https://mic.eucast.org/search/`). The EUCAST rationale
documents do not carry the values either; the ciprofloxacin one (v2.0, 2021-01-01) says in full:

    "2. MIC distributions and epidemiological cut-off (ECOFF) values
     MIC distributions and ECOFFs can be found at https://mic.eucast.org/..."

That is a redirect, not a value. See `wiki/ecoff_sources_2026-09-11.json` for the attempt log.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

# The drugs the evaluation would cover -- the four E. coli agents in the frozen clinical catalog that
# also have numeric MICs in the Oxford deposit.
SOURCED = "sourced"
UNSOURCED = "unsourced"


class UnsourcedEcoffError(RuntimeError):
    """Raised when an ECOFF is requested that carries no verifiable source.

    Deliberately an ERROR and not a warning-plus-default: a caller that proceeds on a defaulted
    cut-off would produce numbers indistinguishable from real ones.
    """


class UnknownDrugError(KeyError):
    """Raised for a drug absent from the catalog entirely."""


@dataclass(frozen=True)
class EcoffEntry:
    """One ECOFF, or one honest hole where an ECOFF would go."""

    drug: str
    organism: str
    value: Optional[float] = None          # mg/L; None while unsourced
    source_url: Optional[str] = None
    eucast_version: Optional[str] = None
    retrieved: Optional[str] = None
    verbatim_quote: Optional[str] = None
    status: str = UNSOURCED
    note: str = ""

    def is_usable(self) -> bool:
        """A usable entry has a value AND the provenance that makes the value checkable."""
        return (self.status == SOURCED
                and self.value is not None
                and bool(self.source_url)
                and bool(self.verbatim_quote))


# ---------------------------------------------------------------------------
# The catalog. Add a value ONLY together with the URL it came from and the text quoted verbatim from
# that page. An entry with a value but no quote is rejected by `validate_catalog`, which is what stops
# a remembered number from being pasted in later under the appearance of provenance.
# ---------------------------------------------------------------------------
# SOURCED 2026-09-11. The route that worked, recorded because the obvious ones do not: the EUCAST
# search IS addressable by query parameter -- the earlier failure was a wrong species id, not a JS wall.
# `search[species]=261` is Escherichia coli (taken from EUCAST's own VetCAST page, then CONFIRMED by
# fetching it and seeing "Escherichia coli" returned). Each agent row links to /search/diagram/<id>,
# which serves a PNG: the ECOFF is RENDERED INTO THE GRAPH IMAGE, in its bottom-left corner. That is why
# no amount of HTML scraping finds it, and it is why these values were read from the images themselves
# rather than from text. Read from the primary source; NOT recalled, NOT taken from a search summary.
_DB = "International MIC distribution - Reference database 2026-09-12, based on aggregated distributions"

ECOFFS: dict[str, EcoffEntry] = {
    "ciprofloxacin": EcoffEntry(
        "ciprofloxacin", "Escherichia coli", value=0.06, status=SOURCED,
        source_url="https://mic.eucast.org/search/diagram/763",
        eucast_version=_DB, retrieved="2026-09-11",
        verbatim_quote="Epidemiological cut-off (ECOFF): 0.06 mg/L / Wildtype (WT) organisms: "
                       "<= 0.06 mg/L / Confidence interval: 0.03 - 0.06 / 15667 observations "
                       "(53 data sources)",
        note="Established ECOFF (not parenthesised). Well powered."),
    "gentamicin": EcoffEntry(
        "gentamicin", "Escherichia coli", value=2.0, status=SOURCED,
        source_url="https://mic.eucast.org/search/diagram/652",
        eucast_version=_DB, retrieved="2026-09-11",
        verbatim_quote="Epidemiological cut-off (ECOFF): 2 mg/L / Wildtype (WT) organisms: <= 2 mg/L / "
                       "Confidence interval: 1 - 2 / 78136 observations (82 data sources)",
        note="Established ECOFF (not parenthesised). The best-powered of the four."),
    "ceftriaxone": EcoffEntry(
        "ceftriaxone", "Escherichia coli", value=0.125, status=SOURCED,
        source_url="https://mic.eucast.org/search/diagram/4605",
        eucast_version=_DB, retrieved="2026-09-11",
        verbatim_quote="Epidemiological cut-off (ECOFF): (0.125) mg/L / Wildtype (WT) organisms: "
                       "<= 0.125 mg/L / Confidence interval: 0.03 - 0.5 / 908 observations "
                       "(4 data sources)",
        note="TENTATIVE -- printed in PARENTHESES, which is EUCAST's own marker for a TECOFF set on "
             "3-4 distributions rather than the >=5 an ECOFF needs. 908 observations from 4 sources, "
             "against 78,136 from 82 for gentamicin, and a confidence interval spanning four "
             "doublings (0.03-0.5). Treat any ceftriaxone result as provisional on the ANCHOR, "
             "separately from whatever the cohort shows."),
    "tetracycline": EcoffEntry(
        "tetracycline", "Escherichia coli", value=8.0, status=SOURCED,
        source_url="https://mic.eucast.org/search/diagram/3351",
        eucast_version=_DB, retrieved="2026-09-11",
        verbatim_quote="Epidemiological cut-off (ECOFF): 8 mg/L / Wildtype (WT) organisms: <= 8 mg/L / "
                       "Confidence interval: 2 - 4 / 18917 observations (64 data sources)",
        note="Established ECOFF. The printed confidence interval (2-4) does NOT bracket the ECOFF (8); "
             "that is EUCAST's own rendering and was verified against the graph image rather than "
             "assumed to be a transcription error -- the interval describes the fitted wild-type "
             "distribution, not an interval on the cut-off. Do not 'correct' it."),
}


def ecoff_for(drug: str) -> float:
    """The ECOFF in mg/L, or RAISE. Never returns a default.

    `UnknownDrugError` for a drug not in the catalog; `UnsourcedEcoffError` for one present but
    without a verifiable value.
    """
    key = (drug or "").strip().lower()
    if key not in ECOFFS:
        raise UnknownDrugError(
            f"{drug!r} is not in the ECOFF catalog. Known: {sorted(ECOFFS)}")
    entry = ECOFFS[key]
    if not entry.is_usable():
        raise UnsourcedEcoffError(
            f"No sourced ECOFF for {key!r}. {entry.note or ''} "
            "Supply `value` together with `source_url` and `verbatim_quote` before scoring.")
    return float(entry.value)


def entry_for(drug: str) -> EcoffEntry:
    """The full entry, sourced or not -- for reporting what is missing and why."""
    key = (drug or "").strip().lower()
    if key not in ECOFFS:
        raise UnknownDrugError(f"{drug!r} is not in the ECOFF catalog. Known: {sorted(ECOFFS)}")
    return ECOFFS[key]


def sourced_drugs() -> list[str]:
    """Drugs whose ECOFF is usable. Empty is a legitimate answer, and today it IS the answer."""
    return sorted(k for k, v in ECOFFS.items() if v.is_usable())


def unsourced_drugs() -> list[str]:
    return sorted(k for k, v in ECOFFS.items() if not v.is_usable())


def validate_catalog() -> list[str]:
    """Return a list of provenance violations; empty means the catalog is internally honest.

    The violation that matters: a value present WITHOUT a source_url or verbatim_quote. That is the
    shape a remembered number takes when someone fills the catalog in later, so it is checked rather
    than trusted.
    """
    problems: list[str] = []
    for key, e in ECOFFS.items():
        if e.value is not None and not e.source_url:
            problems.append(f"{key}: has a value ({e.value}) but no source_url")
        if e.value is not None and not e.verbatim_quote:
            problems.append(f"{key}: has a value ({e.value}) but no verbatim_quote")
        if e.status == SOURCED and e.value is None:
            problems.append(f"{key}: marked sourced but carries no value")
        if e.value is not None and e.value <= 0:
            problems.append(f"{key}: non-positive ECOFF ({e.value})")
    return problems
