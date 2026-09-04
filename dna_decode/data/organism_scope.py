"""Where a deployed rule is MEASURED to over-call outside the organism it was validated on.

THE GAP THIS EXISTS FOR, found 2026-09-03. The gentamicin v2 `symbol_rescue` was validated on E. coli
(N=131, sens 0.523 -> 0.892). But `calibrated_rule_for` has NO gentamicin entry for any organism, so every
organism falls through to the default `DRUG_RULE` and gets the rescue -- the rule is ORGANISM-AGNOSTIC in
code while its evidence is organism-specific. BV-BRC then measured what that costs: in Klebsiella,
`rmt` carriers are 58R/64S, a PPV of 0.475, so the rescue over-calls on more than half of them. In E. coli
it is 12/12 here and 146/146 on NCBI-PD.

WHY THIS IS L2 AND NOT AN L1 FIX. Restricting the rescue to E. coli would edit the frozen surface and
invalidate the v2 prospective lock -- a user-authority call, exactly as deploying the rescue was. L2
qualifies a call without competing with it, so an over-call warning can ship immediately while the L1
question is surfaced rather than decided. This module NEVER changes a call, a tier, or a metric.

EVERY FIELD IS TRACEABLE. Nothing here is asserted from mechanism or memory: the counts come from the
committed BV-BRC artifacts, the labels survived a pre-registered `aac(3)` control, and the carrier calls
were confirmed by a SECOND tool (AMRFinder agreed with CARD on 66/66).

SCOPE QUALIFICATION ADDED 2026-09-03, and it narrows what this warning may claim. Running this project's
OWN source-diversity bar (`source_diverse_validate.MAX_SOURCE_SHARE`, 0.60) over the same artifact, the
Klebsiella carrier set FAILS it: largest-source share 0.664, and 98.4% of the SUSCEPTIBLE carriers --
the ones carrying the entire finding -- come from a single study. Excluding that source, Klebsiella
carriers are 40R/1S (PPV 0.976), indistinguishable from the validated E. coli scope. The warning STAYS
(an over-cautious over-call warning costs far less than a missing one, and inside that study the signal
looks real: it splits its own carriers 18R/63S rather than calling them all susceptible, and the
pre-registered control passed). What changed is the CLAIM: the over-call is established in ONE
population, NOT for Klebsiella generally. See `wiki/rmt_source_concentration_2026-09-03.json`.
"""
from __future__ import annotations

# (drug, organism-key-fragment) -> the measured over-call. Organism matching is substring-insensitive
# because callers pass wildly different organism strings ("Klebsiella", "Klebsiella_pneumoniae",
# "Klebsiella pneumoniae subsp. pneumoniae").
OVERCALL_SCOPE: tuple[dict, ...] = (
    {
        "drug": "gentamicin",
        "organism_match": "klebsiella",
        "rule_component": "symbol_rescue ^(rmt[A-H]\\d*|npmA\\d*)$",
        "validated_on": "Escherichia coli / Shigella",
        "measured_ppv_here": 0.475,
        "measured_counts": {"R": 58, "S": 64},
        "validated_scope_ppv": 1.000,
        "validated_scope_counts": {"R": 12, "S": 0},
        "archive": "BV-BRC (measured MICs; independent of NCBI-PD)",
        "artifact": "wiki/gentamicin_rmt_bvbrc_2026-09-03.md",
        "controls_passed": [
            "carrier call is not the artefact -- susceptible carriers have BETTER CARD hits "
            "(identity >=99, coverage 100, 0 partial) than resistant ones (17 partial)",
            "labels are not the artefact -- the dominant study calls aac(3) carriers R 99% vs 83% "
            "elsewhere and no-determinant isolates R 4% (verdict SPECIFIC_TO_RMT)",
            "a SECOND caller agrees -- AMRFinder independently calls an RMTase in 66/66",
        ],
        "not_settled": "the mechanism: a full-length rmtB at gentamicin MIC <=1 is unexplained "
                       "(silencing / expression / plasmid context untested)",
        # The measured over-call does NOT clear this project's own source-diversity bar. Carried as a
        # first-class field so the warning can never be read as an archive-level claim.
        "source_concentration": {
            "n_sources": 6,
            "largest_source_share": 0.664,
            "susceptible_largest_source_share": 0.984,
            "own_bar": 0.60,
            "passes_own_bar": False,
            "ppv_excluding_largest_source": 0.976,
            "counts_excluding_largest_source": {"R": 40, "S": 1},
            "meaning": "the over-call is established in ONE study's population, NOT for Klebsiella "
                       "generally; excluding that source the PPV is indistinguishable from the "
                       "validated E. coli scope",
            "artifact": "wiki/rmt_source_concentration_2026-09-03.json",
        },
    },
)


def _norm(s: str | None) -> str:
    return (s or "").strip().lower().replace("_", " ")


def overcall_for(drug: str, organism: str | None) -> dict | None:
    """The measured over-call disclosure for this (drug, organism), or None.

    None means "nothing measured here" -- NEVER "measured and clean". A caller that needs that
    distinction should read the artifact; this returns a warning or nothing.
    """
    d, o = _norm(drug), _norm(organism)
    if not o:
        return None
    for row in OVERCALL_SCOPE:
        if _norm(row["drug"]) != d:
            continue
        if row["organism_match"] in o:
            return {
                "status": "measured_overcall_outside_validated_scope",
                "rule_component": row["rule_component"],
                "validated_on": row["validated_on"],
                "ppv_in_this_organism": row["measured_ppv_here"],
                "counts_in_this_organism": dict(row["measured_counts"]),
                "ppv_in_validated_scope": row["validated_scope_ppv"],
                "counts_in_validated_scope": dict(row["validated_scope_counts"]),
                "archive": row["archive"],
                "controls_passed": list(row["controls_passed"]),
                "not_settled": row["not_settled"],
                "artifact": row["artifact"],
                "source_concentration": dict(row["source_concentration"]),
                "note": ("DISCLOSURE ONLY -- this never changes the call, the tier, or any metric. The "
                         "rule still fires; this says its positive predictive value was MEASURED far "
                         "lower in this organism than in the one it was validated on."),
            }
    return None


def one_line(block: dict | None) -> str | None:
    """Human-readable warning, or None when there is nothing measured to say."""
    if not block:
        return None
    sc = block.get("source_concentration") or {}
    # The concentration caveat is NOT optional prose -- without it the line reads as an archive-level
    # claim, which is exactly the over-reach this field was added to prevent.
    caveat = ""
    if sc and not sc.get("passes_own_bar", True):
        caveat = (f" SINGLE-SOURCE-DOMINATED: {sc['susceptible_largest_source_share']:.0%} of the "
                  f"susceptible carriers come from one study and excluding it the PPV is "
                  f"{sc['ppv_excluding_largest_source']:.3f}, so this is established in ONE population, "
                  f"not for this organism generally.")
    return (f"ORGANISM-SCOPE WARNING: this rule's {block['rule_component']} was validated on "
            f"{block['validated_on']} (PPV {block['ppv_in_validated_scope']:.3f}) but MEASURED at PPV "
            f"{block['ppv_in_this_organism']:.3f} here "
            f"({block['counts_in_this_organism']['R']}R/{block['counts_in_this_organism']['S']}S) -- "
            f"a resistant call on this determinant is the least trustworthy kind.{caveat} "
            f"See {block['artifact']}")
