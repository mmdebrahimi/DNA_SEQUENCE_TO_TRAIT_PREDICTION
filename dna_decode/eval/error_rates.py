"""VME / ME — the same confusion matrix, expressed in the direction that can hurt someone.

WHAT THIS ADDS, STATED HONESTLY: no new information. `VME = 1 - sens` and `ME = 1 - spec` EXACTLY, from
counts the report card already carries. This module exists for FRAMING, and framing is the whole point of
the clinical-microbiology convention it borrows: an AST report is read by someone choosing a drug, so the
headline number is the rate at which the method would tell them a resistant organism is susceptible.

Why that reframing is not cosmetic here, on this project's own committed numbers:

    Klebsiella x meropenem -- sens 0.467  reads as "moderate, needs work"
                              VME  0.533  reads as "more than half of carbapenem-resistant isolates
                                          would be reported susceptible"

Same two numbers. The second is the one a reader can act on, and it is the one sens/spec buries.
Sourced from the clinical-microbiology convention recorded in `wiki/prior_art_decoder_landscape_2026-09-03.md`
(ESCMID/AMRrules use VME/ME as the headline metric rather than sens/spec).

VOCABULARY, and the direction matters:
  * VME (very major error)  = R called S  = FN / (TP + FN). The DANGEROUS direction: under-treatment.
  * ME  (major error)       = S called R  = FP / (FP + TN). Costly (unnecessary broad-spectrum therapy),
                              but it does not leave an infection untreated.

DELIBERATELY NO ACCEPTANCE BAR. Regulatory frameworks do publish numeric VME/ME ceilings for device
clearance, and this module does NOT assert one: the exact thresholds and the denominators they apply to
are not verifiable from anything in this repo, and writing a remembered number beside a real measurement
is the fabrication hazard this project guards against. Rates are reported; judging them against a
standard is a sourcing job nobody has done here yet.

Pure. No I/O, no network. Wilson intervals reuse `clonality.wilson_ci` rather than reimplementing it.
"""
from __future__ import annotations

from dataclasses import dataclass

from .clonality import wilson_ci


@dataclass(frozen=True)
class ErrorRates:
    """VME/ME with their own denominators. A rate whose denominator is 0 is None, never 0.0."""
    vme: float | None
    vme_n: int                      # the resistant isolates -- VME's denominator
    vme_ci: tuple[float, float] | None
    me: float | None
    me_n: int                       # the susceptible isolates -- ME's denominator
    me_ci: tuple[float, float] | None

    def as_dict(self) -> dict:
        return {
            "vme": self.vme, "vme_n_resistant": self.vme_n,
            "vme_ci": list(self.vme_ci) if self.vme_ci else None,
            "me": self.me, "me_n_susceptible": self.me_n,
            "me_ci": list(self.me_ci) if self.me_ci else None,
            "definition": ("VME = R called S = FN/(TP+FN), the dangerous direction; "
                           "ME = S called R = FP/(FP+TN). VME = 1-sens and ME = 1-spec exactly -- "
                           "this is a reframing of the committed counts, not a new measurement."),
        }


def vme_me(tp: int | None, fp: int | None, tn: int | None, fn: int | None) -> ErrorRates | None:
    """Error rates from a confusion matrix, or None when the counts are absent.

    Each rate carries its OWN denominator. Pooling them into a single "error rate" would hide exactly
    the asymmetry the convention exists to expose -- a cell can be near-perfect on one direction and
    dangerous on the other, which is the shape of this project's real cells.
    """
    if any(v is None for v in (tp, fp, tn, fn)):
        return None
    tp, fp, tn, fn = int(tp), int(fp), int(tn), int(fn)
    n_r, n_s = tp + fn, fp + tn
    return ErrorRates(
        vme=(fn / n_r) if n_r else None, vme_n=n_r,
        vme_ci=wilson_ci(fn, n_r) if n_r else None,
        me=(fp / n_s) if n_s else None, me_n=n_s,
        me_ci=wilson_ci(fp, n_s) if n_s else None,
    )


def worst_vme_first(cells: list[dict]) -> list[dict]:
    """Cells ordered by VME descending -- the dangerous direction at the top.

    Ordering by the error rather than by name is the other half of the reframing: a table sorted
    alphabetically makes a 0.533 VME sit wherever the alphabet puts it.
    """
    return sorted((c for c in cells if c.get("vme") is not None),
                  key=lambda c: (-c["vme"], c.get("organism", ""), c.get("drug", "")))
