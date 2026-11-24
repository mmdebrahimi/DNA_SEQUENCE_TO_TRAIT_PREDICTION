"""Is an eval set CAPABLE of carrying a model-vs-null verdict at all? Ask before scoring.

THE POINT. A model-minus-null delta only means something if at least one of them beats a random
guess. Measured 2026-11-22 on the orthology-free (accessory-gene) set, the diffusion model scored
1.4098 nats and the strong order-3 null 1.4291, both WORSE than uniform ln(4)=1.3863: a model
confident about the wrong base loses to one that commits to nothing. Ranking two sub-chance
predictors is the same error class as a broken detector reading as a clean result, so the only
honest move is to REFUSE a verdict.

WHY THIS IS A MODULE AND NOT A COMMENT. The logic already existed, inline, at
`scripts/glm_beyond_kmer.py:228-245` -- but only AFTER a checkpoint had been loaded and every
window scored. Extracted here it is answerable from the NULL ALONE (`model_ce_nats=None`), which
makes a candidate eval set screenable for free, before any GPU work. That is the entire delta:
this is an EXTRACTION, not an invention, and the committed refusal is pinned by test.

WHAT THIS DOES NOT DO. It says a set is ADMISSIBLE, never that a set is UNCONFOUNDED. A set can
clear uniform and still leak near-orthologs from training; admissibility is necessary, not
sufficient. No verdict about any model belongs in this module.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

# Uniform cross-entropy for a 4-letter alphabet, in nats. The comparison constant, not a threshold:
# it is what a predictor that commits to nothing achieves, so it is the floor any informative
# predictor must clear. Kept a keyword argument on the predicate so a different alphabet size is
# expressible without editing this module.
UNIFORM_CE_NATS_DNA = math.log(4)

ADMISSIBLE = "ADMISSIBLE"
REFUSED_NULL_BELOW_UNIFORM = "REFUSED_NULL_BELOW_UNIFORM"
REFUSED_MODEL_BELOW_UNIFORM = "REFUSED_MODEL_BELOW_UNIFORM"
REFUSED_BOTH_BELOW_UNIFORM = "REFUSED_BOTH_BELOW_UNIFORM"

#: Every state this predicate can emit. A caller rendering states should read this rather than
#: restating the list, which is how a fourth state goes unhandled.
STATES = (ADMISSIBLE, REFUSED_NULL_BELOW_UNIFORM, REFUSED_MODEL_BELOW_UNIFORM,
          REFUSED_BOTH_BELOW_UNIFORM)


@dataclass(frozen=True)
class Admissibility:
    """The verdict plus the margins it was decided on, so a caller never re-derives them."""

    state: str
    uniform_ce_nats: float
    null_ce_nats: float
    null_margin_nats: float
    model_ce_nats: float | None
    model_margin_nats: float | None
    model_assessed: bool

    @property
    def admissible(self) -> bool:
        return self.state == ADMISSIBLE

    def to_dict(self) -> dict:
        return {
            "state": self.state,
            "admissible": self.admissible,
            "uniform_ce_nats": round(self.uniform_ce_nats, 4),
            "null_ce_nats": self.null_ce_nats,
            "null_margin_nats": round(self.null_margin_nats, 4),
            "model_ce_nats": self.model_ce_nats,
            "model_margin_nats": (None if self.model_margin_nats is None
                                  else round(self.model_margin_nats, 4)),
            "model_assessed": self.model_assessed,
        }


def admissibility(
    null_ce_nats: float,
    model_ce_nats: float | None = None,
    *,
    model_uniform_delta_nats: float | None = None,
    uniform: float = UNIFORM_CE_NATS_DNA,
) -> Admissibility:
    """Can this eval set carry a verdict? Decided on cross-entropy against uniform.

    `null_ce_nats` is the STRONG null's mean cross-entropy in nats. `model_ce_nats` is OPTIONAL and
    defaults to None -- the null-only call is the pre-scoring screen, and it asks only whether the
    set is predictable at all by something that is not the model.

    `model_uniform_delta_nats` is the model's PAIRED delta against the uniform null (positive =
    model better). Supply it instead of `model_ce_nats` to reproduce `glm_beyond_kmer.py`'s original
    test bit-for-bit: that scorer tested `uni_cell["paired_delta_nats"] > 0`, a value rounded to 4
    decimals BEFORE comparison, whereas `uniform - round(model_ce, 4)` rounds after. The two agree
    except within ~5e-5 of the boundary, and passing the delta removes even that difference. It
    takes precedence over `model_ce_nats` when both are given.

    Margins are `uniform - x`, so POSITIVE means "beats a random guess". The comparison is STRICT:
    a predictor sitting exactly at uniform has learned nothing and does not clear the bar.
    """
    null_margin = uniform - float(null_ce_nats)
    null_ok = null_margin > 0

    if model_uniform_delta_nats is not None:
        model_margin: float | None = float(model_uniform_delta_nats)
        model_ce_out = None if model_ce_nats is None else float(model_ce_nats)
        model_assessed = True
    elif model_ce_nats is not None:
        model_ce_out = float(model_ce_nats)
        model_margin = uniform - model_ce_out
        model_assessed = True
    else:
        model_ce_out = None
        model_margin = None
        model_assessed = False

    # An unassessed model cannot disqualify the set. The null-only screen is deliberately WEAKER
    # than the full gate: it can rule a set out before any scoring, never rule one in for good.
    model_ok = True if model_margin is None else model_margin > 0

    if not null_ok and not model_ok:
        state = REFUSED_BOTH_BELOW_UNIFORM
    elif not null_ok:
        state = REFUSED_NULL_BELOW_UNIFORM
    elif not model_ok:
        state = REFUSED_MODEL_BELOW_UNIFORM
    else:
        state = ADMISSIBLE

    return Admissibility(
        state=state,
        uniform_ce_nats=uniform,
        null_ce_nats=float(null_ce_nats),
        null_margin_nats=null_margin,
        model_ce_nats=model_ce_out,
        model_margin_nats=model_margin,
        model_assessed=model_assessed,
    )
