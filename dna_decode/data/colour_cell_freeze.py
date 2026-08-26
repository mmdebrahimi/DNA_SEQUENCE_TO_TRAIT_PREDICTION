"""The animal colour/plumage family is FROZEN at 19 cells (2026-08-26).

WHY. The family reached 19 CLI cells / 65 loci, every one shipping as KNOWLEDGE_BASELINE, before anyone
asked whether they COULD be validated. `dna_decode/pigment/substrate_screen.py` derives the answer from
the committed catalogs, and the family is blocked twice over:

  * 40 of 65 loci (62%) record NO causal variant at all -- only allele symbols and dominance order, as in
    "rabbit C (TYR): C full > chinchilla > Himalayan > c albino". SEVEN of 19 cells record none for ANY
    locus (alpaca, cattle, mouse, pig, pigeon, rabbit, sheep). Those cells are UNVALIDATABLE AS WRITTEN:
    a locus whose causal variant is unspecified cannot be scored against any genotype file, so no cohort
    helps. The blocker is CURATION, not data.
  * Of the 25 loci that DO record one, 14 (56%) are indel/structural -- off-panel for any SNP array or
    imputed biallelic-SNV panel. That is empirically what sank the one measured cell: the dog cell scored
    black 160/161 = 0.994 on Darwin's Ark and left every other base colour unscorable.

So adding colour cell #20 adds a rule that cannot be validated on any substrate. This module makes that
cost VISIBLE at the moment someone tries, instead of after nineteen more.

HONESTY -- THIS IS AN ATTENTION/SCOPE FREEZE, NOT ENFORCEMENT. It makes adding a cell FAIL A TEST
(`tests/test_colour_cell_freeze.py`), which a determined edit can still ratify by updating
`FROZEN_COLOUR_ROUTES` deliberately. That is the intended friction: a scope decision should be explicit,
not silent. Do NOT describe this as a hard invariant -- the same rail as the self-init population cap.

Gate references are G9/G10 in `wiki/negative_results_map_2026-06-13.md`: G1-G8 gate whether a usable
LABEL exists; G9/G10 gate whether the decoder's own rule is SCOREABLE against a genotype at all.
"""
from __future__ import annotations

from dna_decode.pigment.substrate_screen import collect, summarise, trait_for_species

FREEZE_DATE = "2026-08-26"

FREEZE_RATIONALE = (
    "40 of 65 colour-cell loci (62%) record no causal variant at all, leaving 7 of 19 cells unvalidatable "
    "as written; of the 25 loci that do record one, 14 (56%) are indel/structural and off-panel for any "
    "SNP array. Adding a 20th cell adds a rule that cannot be validated on any substrate. "
    "See wiki/colour_cell_substrate_screen_2026-08-26.md and gates G9/G10 in "
    "wiki/negative_results_map_2026-06-13.md."
)

# The 19 CLI trait names, as of the freeze. This literal IS the deliberate friction point -- but a
# hand-listed roster drifts silently (third instance of that bug class in this repo), so
# `tests/test_colour_cell_freeze.py` cross-checks it against the LIVE catalogs and the live CLI registry.
# The literal is what a human must consciously edit; the derived check is what stops it rotting.
FROZEN_COLOUR_ROUTES: frozenset[str] = frozenset({
    "alpacacolor", "buffalocolor", "camelcolor", "catcolor", "cattlecolor", "coatcolor",
    "donkeycolor", "foxcolor", "goatcolor", "guineapigcolor", "horsecolor", "minkcolor",
    "mousecolor", "pigcolor", "pigeoncolor", "plumage", "rabbitcolor", "roedeercolor",
    "sheepcolor",
})

# Headline counts, pinned in ONE place so the memo and the code cannot diverge silently.
EXPECTED_TOTALS = {"n_cells": 19, "n_loci": 65, "n_unrecorded": 40, "n_snv_panel_blocked": 14}


def gates_for(n_unrecorded: int, n_snv_panel_blocked: int) -> tuple[str, ...]:
    """Which rejection gates apply to a cell, from its locus counts. PURE.

    DERIVED FROM COUNTS, not from the verdict string. The plan sketched a `gates_for_verdict(verdict)`
    signature, but its own spec for PARTIALLY_SNV_TRACTABLE ("G9+G10 when it has BOTH unrecorded and
    off-panel loci, else the applicable one") cannot be satisfied from a verdict alone -- the verdict
    does not carry the two counts. Deriving from counts is the honest shape.

      G9  = causal variant unrecorded  -> the rule is not scoreable against any genotype
      G10 = variant class off-panel    -> recorded, but a SNP array cannot represent it
    """
    gates = []
    if n_unrecorded > 0:
        gates.append("G9")
    if n_snv_panel_blocked > 0:
        gates.append("G10")
    return tuple(gates)


def screen_summaries() -> dict[str, dict]:
    """CLI trait name -> that cell's full screen summary (counts + verdict). Derived, never re-typed."""
    return {trait_for_species(sp): summarise(rows) for sp, rows in collect().items()}


def freeze_status(route: str) -> dict:
    """Is this trait frozen, and what did the screen say about it?

    `route` accepts either the bare trait ("rabbitcolor") or the CLI route ("dna-rabbitcolor").
    An unknown trait returns frozen=False with no verdict -- it is not in the family, so the freeze has
    nothing to say about it.
    """
    trait = route.removeprefix("dna-")
    summaries = screen_summaries()
    s = summaries.get(trait)
    if s is None:
        return {"trait": trait, "frozen": False, "screen_verdict": None, "gates": ()}
    return {
        "trait": trait,
        "frozen": trait in FROZEN_COLOUR_ROUTES,
        "screen_verdict": s["verdict"],
        "gates": gates_for(s["n_unrecorded"], s["n_snv_panel_blocked"]),
    }
