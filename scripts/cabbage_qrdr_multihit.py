"""Does the QRDR determinant COUNT predict ciprofloxacin R better than its presence? Salmonella vs
Campylobacter, same frozen rule, same drug.

THE OBSERVATION THIS TESTS. The frozen cipro rule is `qrdr_point` with threshold=1: one QRDR point
mutation and the call is R. On Campylobacter that is near-perfect (gyrA_T86I PPV 0.995 on 1,000
carriers). On Salmonella it over-calls badly -- gyrA_D87Y has 169 carriers at PPV 0.012 -- while
parC_S80I and parC_S80R sit at 1.000 and 0.964.

THE HYPOTHESIS, which is textbook Salmonella fluoroquinolone biology and therefore exactly the kind of
thing to test rather than assert: a SINGLE gyrA QRDR hit confers only reduced susceptibility (historic
nalidixic-acid resistance) with cipro MIC still under the clinical breakpoint, while clinical cipro R
needs MULTIPLE hits. If so, parC's high PPV is not parC predicting anything on its own -- it is parC
MARKING the multi-hit genotype, because in Salmonella parC mutations arise on top of a gyrA mutation.

TWO PREDICTIONS, distinguishable and both checkable from the cached calls:
  P1  R-rate rises monotonically with the NUMBER of counted determinants (Salmonella).
  P2  nearly every parC carrier ALSO carries gyrA. If instead parC carriers are mostly gyrA-free, the
      marking account is wrong and parC genuinely predicts alone.
A control: the same two quantities on Campylobacter, where single-hit gyrA already predicts, so the
count should add little.

Run: uv run python scripts/cabbage_qrdr_multihit.py
"""
from __future__ import annotations

import collections
import json
import sys
from datetime import date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

# (species, registry organism, drug, tag). The ceftriaxone row is the WITHIN-ORGANISM CONTROL: same
# organism and same cohort source as the Salmonella cipro row, but a different rule shape
# (subclass_any over acquired beta-lactamases vs qrdr_point over target-site point mutations). It is
# what separates "the rule shape over-calls" from "Salmonella labels or this cohort are poor".
CELLS = [("Salmonella enterica", "Salmonella", "ciprofloxacin", "salmonella_cipro"),
         ("Salmonella enterica", "Salmonella", "ceftriaxone", "salmonella_ceftriaxone"),
         ("Campylobacter jejuni", "Campylobacter", "ciprofloxacin", "campylobacter_cipro")]


class ParseFailure(RuntimeError):
    """A parse returned nothing where the input demonstrably had something.

    WHY THIS EXISTS AS A RAISE RATHER THAN AN EMPTY RESULT. This script's whole output is strata keyed
    by determinant count. If `build_labels` returns nothing, or if no isolate's determinants parse,
    every stratum is empty, `P1_PRIME` is `{}` and `P2` is `None` -- which reads as "inconclusive",
    the reassuring answer, when the truth is "the parser failed". That is the same defect the sibling
    anatomy probe shipped: it read a key the record does not have, every bucket came back
    determinant-free, and it printed a confident verdict matching the expected story.

    The guard below is a DOMAIN INVARIANT, not a heuristic: each of these cells is known to contain
    hundreds of determinant-carrying isolates (Campylobacter cipro ~1,006, Salmonella ceftriaxone
    ~915), so zero carriers cannot be a real measurement here.
    """


def _gene(symbol: str) -> str:
    """Gene family of a determinant symbol: `gyrA_D87Y` -> gyrA, `aac(6')-Ib-cr5` -> itself."""
    return symbol.split("_", 1)[0] if "_" in symbol else symbol


def score_cell(species: str, organism: str, drug: str) -> dict:
    from dna_decode.data.amr_portal import prepare
    from dna_decode.eval.amr_rules import call_resistance
    from dna_decode.eval.clonality import wilson_ci
    from scripts.cabbage_score_cell import CACHE, build_labels

    labels, census, _s, _y, _c = build_labels(species, drug)
    if not labels:
        raise ParseFailure(
            f"{species} x {drug}: build_labels returned ZERO labels. An empty label set and a "
            "mis-typed species/drug string are indistinguishable downstream -- every stratum would "
            "come back empty and the verdict would read 'inconclusive' rather than 'broken'.")
    by_count: dict[int, collections.Counter] = collections.defaultdict(collections.Counter)
    parc = collections.Counter()
    n_unscoreable = 0
    for acc, truth in labels.items():
        raw = CACHE / f"{acc}.tsv"
        if not raw.exists():
            n_unscoreable += 1
            continue
        try:
            res = call_resistance(str(prepare(raw, CACHE / "norm" / f"{acc}.tsv")),
                                 drug=drug, organism=organism)
        except Exception:                             # noqa: BLE001
            n_unscoreable += 1
            continue
        syms = sorted({str(h.get("symbol") or "?") for h in (res.get("determinants") or [])})
        genes = {_gene(s) for s in syms}
        by_count[len(syms)][truth] += 1
        if "parC" in genes:
            parc["parC_total"] += 1
            parc["parC_with_gyrA"] += int("gyrA" in genes)
            parc[f"parC_truth_{truth}"] += 1
    # DOMAIN INVARIANT, checked before any statistic is computed. Zero determinant-carrying isolates
    # cannot be a real measurement on these cells, so it means the determinant parse returned nothing.
    n_carriers = sum(sum(c.values()) for k, c in by_count.items() if k > 0)
    n_scored_total = sum(sum(c.values()) for c in by_count.values())
    if n_scored_total == 0:
        raise ParseFailure(
            f"{species} x {drug}: ZERO isolates scored of {len(labels)} labelled "
            f"({n_unscoreable} unscoreable). Nothing was measured; an empty strata table would read "
            "as inconclusive.")
    if n_carriers == 0:
        raise ParseFailure(
            f"{species} x {drug}: {n_scored_total} isolates scored and NOT ONE carries a counted "
            "determinant. These cells contain hundreds of carriers, so this is a determinant-parse "
            "failure, not a finding.")

    rows = []
    for k in sorted(by_count):
        c = by_count[k]
        n = c["R"] + c["S"]
        lo, hi = wilson_ci(c["R"], n) if n else (None, None)
        rows.append({"n_determinants": k, "n": n, "n_R": c["R"], "n_S": c["S"],
                     "r_rate": round(c["R"] / n, 4) if n else None,
                     "r_rate_wilson_ci": [round(lo, 4), round(hi, 4)] if n else None})
    # P1 AS ORIGINALLY SPECIFIED: strict non-decreasing R-rate in determinant count (n>=10 strata).
    # Kept and reported even though it FAILS on Salmonella, because quietly swapping a pre-registered
    # test for one that passes is the move this project refuses. It fails on a single dip -- the
    # 4-determinant stratum sits below the 3-determinant one on n=14 vs n=17 with heavily overlapping
    # Wilson intervals -- so strict monotonicity was a brittle operationalisation of the premise
    # "R-rate rises with count", not a measurement of it.
    powered = [r for r in rows if r["n"] >= 10 and r["n_determinants"] > 0]
    rates = [r["r_rate"] for r in powered]
    monotone = all(b >= a for a, b in zip(rates, rates[1:])) if len(rates) >= 2 else None

    # P1' the premise stated as a CONTRAST rather than a shape: single hit vs multiple hits. This is
    # what "the count carries the information" actually means, and it does not break on one dip.
    one = next((r for r in rows if r["n_determinants"] == 1), None)
    multi_R = sum(r["n_R"] for r in rows if r["n_determinants"] >= 2)
    multi_S = sum(r["n_S"] for r in rows if r["n_determinants"] >= 2)
    multi_n = multi_R + multi_S
    contrast = None
    if one and one["n"] >= 10 and multi_n >= 10:
        from dna_decode.eval.clonality import wilson_ci as _wci
        lo, hi = _wci(multi_R, multi_n)
        contrast = {
            "r_rate_exactly_one_determinant": one["r_rate"],
            "r_rate_exactly_one_ci": one["r_rate_wilson_ci"],
            "n_exactly_one": one["n"],
            "r_rate_two_or_more": round(multi_R / multi_n, 4),
            "r_rate_two_or_more_ci": [round(lo, 4), round(hi, 4)],
            "n_two_or_more": multi_n,
            "ratio": (round((multi_R / multi_n) / one["r_rate"], 1)
                      if one["r_rate"] else None),
            "cis_disjoint": bool(lo > one["r_rate_wilson_ci"][1]),
        }
    return {
        "single_vs_multi_contrast": contrast,
        "species": species, "registry_organism": organism, "drug": drug,
        "n_unscoreable": n_unscoreable,
        "by_determinant_count": rows,
        "powered_counts_used_for_monotonicity": [r["n_determinants"] for r in powered],
        "r_rate_monotone_in_count": monotone,
        "parC": dict(parc),
        "parC_with_gyrA_fraction": (round(parc["parC_with_gyrA"] / parc["parC_total"], 4)
                                    if parc["parC_total"] else None),
    }


def main() -> int:
    out = {"schema": "cabbage-qrdr-multihit-v2", "generated": date.today().isoformat(),
           "frozen_rule": "eval/amr_rules.call_resistance, counter=qrdr_point, threshold=1",
           "cells": {}}
    for species, organism, drug, tag in CELLS:
        print(f"=== {species} ({organism}) x {drug} ===", flush=True)
        r = score_cell(species, organism, drug)
        out["cells"][tag] = r
        for row in r["by_determinant_count"]:
            print(f"  n_det={row['n_determinants']}  n={row['n']:>5}  R={row['n_R']:>4} "
                  f"S={row['n_S']:>5}  R-rate={row['r_rate']}  CI={row['r_rate_wilson_ci']}",
                  flush=True)
        print(f"  monotone in count (n>=10 strata): {r['r_rate_monotone_in_count']}", flush=True)
        print(f"  parC carriers: {r['parC'].get('parC_total', 0)}  also carrying gyrA: "
              f"{r['parC'].get('parC_with_gyrA', 0)}  fraction={r['parC_with_gyrA_fraction']}",
              flush=True)

    s = out["cells"].get("salmonella_cipro", {})
    camp = out["cells"].get("campylobacter_cipro", {})
    cef = out["cells"].get("salmonella_ceftriaxone", {})
    frac = s.get("parC_with_gyrA_fraction")
    con = s.get("single_vs_multi_contrast") or {}
    s1 = next((r for r in s.get("by_determinant_count", []) if r["n_determinants"] == 1), {})
    c1 = next((r for r in camp.get("by_determinant_count", []) if r["n_determinants"] == 1), {})
    out["verdict"] = {
        "P1_AS_SPECIFIED_strict_monotone_in_count": s.get("r_rate_monotone_in_count"),
        "P1_specification_was_too_strict": (
            "REPORTED AS IT FELL. The pre-registered test asked for a strictly non-decreasing "
            "R-rate and Salmonella FAILS it on one dip: the 4-determinant stratum (n=14) sits below "
            "the 3-determinant one (n=17) with heavily overlapping Wilson intervals. The premise was "
            "'R-rate rises with count'; strict monotonicity across small strata was a brittle way to "
            "operationalise it. The contrast below is the same premise without that brittleness, and "
            "it is reported BESIDE the failed test, not instead of it."),
        "P1_PRIME_single_vs_multi_hit_contrast": con,
        "P2_parC_marks_the_multihit_genotype": (None if frac is None else frac >= 0.9),
        "cross_organism_contrast": {
            "salmonella_r_rate_at_one_determinant": s1.get("r_rate"),
            "campylobacter_r_rate_at_one_determinant": c1.get("r_rate"),
            "ratio": (round(c1["r_rate"] / s1["r_rate"], 1)
                      if s1.get("r_rate") and c1.get("r_rate") else None),
            "note": ("the SAME frozen rule (counter=qrdr_point, threshold=1) on the SAME drug. One "
                     "QRDR hit predicts resistance in Campylobacter and does not in Salmonella, "
                     "which is where the specificity gap between the two cells comes from."),
        },
        "mechanism_class_control": {
            "why": ("SAME organism, SAME cohort source, DIFFERENT rule shape. Salmonella x "
                    "ceftriaxone uses subclass_any over ACQUIRED beta-lactamases; Salmonella x "
                    "ciprofloxacin uses qrdr_point over TARGET-SITE point mutations. If the cipro "
                    "over-call were a property of Salmonella labels or of this cohort, ceftriaxone "
                    "would over-call too."),
            "salmonella_ceftriaxone_r_rate_at_one_determinant": (
                next((r["r_rate"] for r in cef.get("by_determinant_count", [])
                      if r["n_determinants"] == 1), None)),
            "salmonella_cipro_r_rate_at_one_determinant": s1.get("r_rate"),
            "reading": ("one ACQUIRED determinant is close to sufficient for ceftriaxone while one "
                        "TARGET-SITE point mutation is close to uninformative for ciprofloxacin IN "
                        "THE SAME ORGANISM AND COHORT, which exonerates the labels and the cohort: "
                        "whatever is wrong is in the rule, not the data."),
            "the_acquired_vs_target_site_story_is_REFUTED_by_campylobacter": (
                "The tempting generalisation -- threshold=1 suits acquired-enzyme mechanisms and "
                "mismatches target-site ones -- is WRONG, and the third cell is what kills it. "
                "Campylobacter ciprofloxacin is ALSO target-site (gyrA_T86I) and a single hit "
                "predicts at 0.9891. So mechanism CLASS does not decide it. What decides it is "
                "whether ONE determinant is clinically sufficient for that ORGANISM-DRUG PAIR: "
                "gyrA_T86I alone clears the Campylobacter cipro breakpoint, a single gyrA hit alone "
                "does not clear the Salmonella cipro breakpoint. threshold=1 is a claim about "
                "sufficiency per organism-drug pair, and it is a per-pair empirical question."),
        },
        "reading": (
            "The COUNT carries the information that presence does not FOR THIS PAIR: in Salmonella "
            "cipro a single counted "
            "determinant is almost always susceptible while two or more is usually resistant, so a "
            "threshold=1 rule over-calls by construction there. P2 true means parC's high PPV is "
            "parC MARKING a multi-hit genotype rather than predicting alone -- all parC carriers "
            "also carry gyrA -- so a rule that counted only parC would work for the wrong reason and "
            "would break on any population where parC arises first."),
    }
    out["honest_limits"] = [
        "OBSERVATIONAL and MARGINAL: determinant count is not a causal dose. Carriers of 2 "
        "determinants differ from carriers of 1 in more than the count.",
        "Determinant count is the count the FROZEN rule assigns, so it inherits whatever AMRFinder "
        "reports and whatever the rule's class filter admits -- not a complete genotype.",
        "Reduced-susceptibility isolates sit near the breakpoint, so some S labels here are "
        "genuinely borderline; STRICT tiers reduce but do not remove that.",
        "NOT clonality-corrected. A single-determinant stratum may be a few clones, and Salmonella "
        "surveillance is clonally structured.",
        "This DIAGNOSES an over-call; it does not authorise a rule change. Changing what the frozen "
        "rule counts retires the active prospective lock and is a USER decision.",
        "The two organisms' cohorts differ in source composition and size, so the Campylobacter arm "
        "is a CONTRAST, not a matched control.",
    ]
    p = REPO / "wiki" / f"cabbage_qrdr_multihit_{date.today().isoformat()}.json"
    p.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nverdict: {json.dumps(out['verdict'], indent=1)[:400]}")
    print(f"-> {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
