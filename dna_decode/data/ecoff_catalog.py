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
_UNREACHABLE = ("EUCAST publishes no downloadable ECOFF table and no stable per-result URL; the value "
                "is only rendered by the interactive JS database at https://mic.eucast.org/search/, "
                "which this environment cannot drive. Not recalled from memory on purpose.")

ECOFFS: dict[str, EcoffEntry] = {
    "ciprofloxacin": EcoffEntry("ciprofloxacin", "Escherichia coli", note=_UNREACHABLE),
    "gentamicin": EcoffEntry("gentamicin", "Escherichia coli", note=_UNREACHABLE),
    "ceftriaxone": EcoffEntry("ceftriaxone", "Escherichia coli", note=_UNREACHABLE),
    "tetracycline": EcoffEntry("tetracycline", "Escherichia coli", note=_UNREACHABLE),
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
