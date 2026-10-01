"""Kill-count reporting — an inert filter must declare itself.

WHY. This repo has twice shipped a filter that removed nothing and presented it as a control, and both
had to be retracted:

* the HIV censoring guard filtered operator-prefixed folds (`<` / `>`) when there are **zero** such
  values in any column of any Stanford dataset;
* the ResFinder comparison excluded `POINT` rows as a control when there are **zero** POINT rows across
  all 1,818 AMRFinder runs — nothing was excluded, and calling it an active control over-claimed it.

A constraint layer without kill-count reporting reproduces that failure at scale, so every constraint
reports what it actually removed. Three states are DERIVED, never set:

* `active` — refuted at least one case. It is doing work.
* `inert` — evaluated cases and refuted none. **This is not protection.** It must be acknowledged with a
  reason, or `assert_no_silent_inert` raises.
* `not_evaluated` — never ran. Distinct from `inert` for the same reason `inapplicable` is distinct from
  `satisfied`: "did not check" and "checked and found nothing" are different claims.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from dna_decode.constraints.registry import INAPPLICABLE, REFUTED, SATISFIED, ConstraintVerdict

ACTIVE = "active"
INERT = "inert"
NOT_EVALUATED = "not_evaluated"


class SilentInertConstraintError(AssertionError):
    """Raised when a constraint removed nothing and nobody acknowledged it."""


@dataclass
class ConstraintCounts:
    n_evaluated: int = 0
    n_refuted: int = 0
    n_satisfied: int = 0
    n_inapplicable: int = 0

    @property
    def status(self) -> str:
        if self.n_evaluated == 0:
            return NOT_EVALUATED
        return ACTIVE if self.n_refuted > 0 else INERT


@dataclass
class ConstraintReport:
    """Per-constraint counts plus the acknowledgement ledger for inert constraints."""

    counts: dict[str, ConstraintCounts] = field(default_factory=dict)
    acknowledged: dict[str, str] = field(default_factory=dict)

    def record(self, key: str, verdict: ConstraintVerdict) -> None:
        c = self.counts.setdefault(key, ConstraintCounts())
        if verdict.verdict == INAPPLICABLE:
            c.n_inapplicable += 1
            return                     # an inapplicable case was not evaluated against the rule
        c.n_evaluated += 1
        if verdict.verdict == REFUTED:
            c.n_refuted += 1
        elif verdict.verdict == SATISFIED:
            c.n_satisfied += 1

    def note_not_evaluated(self, key: str) -> None:
        """Register a constraint that never ran, so it appears in the artifact rather than vanishing."""
        self.counts.setdefault(key, ConstraintCounts())

    def acknowledge_inert(self, key: str, reason: str) -> None:
        if not reason:
            raise ValueError("an inert acknowledgement needs a reason")
        self.acknowledged[key] = reason

    def status(self, key: str) -> str:
        return self.counts.get(key, ConstraintCounts()).status

    def inert_keys(self) -> list[str]:
        return sorted(k for k, c in self.counts.items() if c.status == INERT)

    def active_keys(self) -> list[str]:
        return sorted(k for k, c in self.counts.items() if c.status == ACTIVE)

    def assert_no_silent_inert(self) -> None:
        """Make inertness a DECISION rather than an oversight."""
        unacked = [k for k in self.inert_keys() if k not in self.acknowledged]
        if unacked:
            raise SilentInertConstraintError(
                "these constraints evaluated cases and refuted none, and were not acknowledged: "
                + ", ".join(unacked)
                + ". A filter that removes nothing is not a control -- either acknowledge it with a "
                  "reason via acknowledge_inert(), or remove it.")

    def as_dict(self) -> dict:
        """Every constraint's counts and status, so an inert one is VISIBLE in the artifact."""
        return {
            "constraints": {
                k: {"n_evaluated": c.n_evaluated, "n_refuted": c.n_refuted,
                    "n_satisfied": c.n_satisfied, "n_inapplicable": c.n_inapplicable,
                    "status": c.status,
                    "acknowledged_reason": self.acknowledged.get(k)}
                for k, c in sorted(self.counts.items())
            },
            "summary": {"n_active": len(self.active_keys()), "n_inert": len(self.inert_keys()),
                        "n_not_evaluated": sum(1 for c in self.counts.values()
                                               if c.status == NOT_EVALUATED)},
            "reading": ("`inert` means the constraint ran and removed nothing -- it is NOT evidence of "
                        "protection. `not_evaluated` means it never ran. The two are reported "
                        "separately on purpose."),
        }
