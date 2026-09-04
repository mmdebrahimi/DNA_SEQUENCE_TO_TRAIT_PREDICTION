"""Apply our OWN source-diversity bar to the Klebsiella over-call number we published.

THE PROBLEM THIS FOUND. On 2026-09-03 we published a measured over-call for the gentamicin `rmt` rescue
in Klebsiella: PPV 0.475 (58R/64S) on BV-BRC, and shipped it as an `organism_scope` warning. The memo
noted that 94% of the susceptible carriers came from one study and argued the `aac(3)` control cleared
it. What the memo did NOT do is apply the bar this project already enforces elsewhere:
`scripts/source_diverse_validate.py` REFUSES to report a cell whose largest source holds more than
`MAX_SOURCE_SHARE` of it, precisely because a concentrated cohort cannot speak for a population.

Run that bar over the same artifact and the Klebsiella cell FAILS it: largest-source share 0.664 against
a 0.60 bar, and 98.4% of the SUSCEPTIBLE carriers -- the ones carrying the entire finding -- come from a
single study. Excluding that source, Klebsiella rmt carriers are 40R/1S, PPV 0.976, which is
indistinguishable from the E. coli scope the rule was validated on.

WHAT THIS DOES AND DOES NOT OVERTURN. It does NOT retract the warning. An over-call warning that turns
out to be over-cautious costs far less than a missing one, and two facts still argue the signal is real
inside that study's population: the study splits its own rmt carriers 18R/63S rather than calling them
all susceptible, and the pre-registered `aac(3)` control passed (99% R inside vs 83% elsewhere). What it
DOES overturn is the SCOPE of the claim. "In Klebsiella the rescue over-calls on more than half of
carriers" is an archive-level claim that 6 sources at 66% concentration cannot support. The honest
claim is narrower: the over-call is established in ONE study's population and is NOT established for
Klebsiella generally.

Offline: reads the committed BV-BRC hunt artifact. Writes wiki/rmt_source_concentration_<date>.json.
"""
from __future__ import annotations

import argparse
import collections
import json
import sys
from datetime import date as _date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

# Import the bar rather than restating it -- a restated threshold drifts from the one that enforces.
from source_diverse_validate import MAX_SOURCE_SHARE  # noqa: E402

ORGANISM_KEYS = ("klebsiella", "escherichia", "salmonella", "enterobacter",
                 "citrobacter", "pseudomonas", "acinetobacter")


def source_of(hit: dict) -> str:
    """BV-BRC's `pmid` is the publication that curated the phenotype -- the source unit here.

    It arrives as a scalar OR a list depending on the record, and a list is unhashable: counting it
    naively raises rather than miscounts, which is the safe direction, but it must be normalised.
    """
    p = hit.get("pmid")
    if isinstance(p, list):
        p = ",".join(str(x) for x in p) if p else ""
    return str(p) if p else "NO_PMID"


def organism_of(hit: dict) -> str:
    name = (hit.get("genome_name") or "").lower()
    return next((k for k in ORGANISM_KEYS if k in name), "other")


def profile(hits: list[dict]) -> dict:
    """Concentration profile for one organism's carrier set, plus the leave-largest-source-out check."""
    scored = [h for h in hits if h.get("phenotype") in ("Resistant", "Susceptible")]
    if not scored:
        return {"n": 0}
    R = sum(1 for h in scored if h["phenotype"] == "Resistant")
    S = len(scored) - R
    src = collections.Counter(source_of(h) for h in scored)
    top, top_n = src.most_common(1)[0]
    sus = [h for h in scored if h["phenotype"] == "Susceptible"]
    sus_src = collections.Counter(source_of(h) for h in sus)

    rest = [h for h in scored if source_of(h) != top]
    rest_R = sum(1 for h in rest if h["phenotype"] == "Resistant")

    share = top_n / len(scored)
    return {
        "n": len(scored), "R": R, "S": S, "ppv": R / len(scored),
        "n_sources": len(src),
        "largest_source": top, "largest_source_share": share,
        "susceptible_largest_source": (sus_src.most_common(1)[0][0] if sus else None),
        "susceptible_largest_source_share": (sus_src.most_common(1)[0][1] / len(sus)) if sus else None,
        "excluding_largest_source": (
            {"n": len(rest), "R": rest_R, "S": len(rest) - rest_R,
             "ppv": rest_R / len(rest)} if rest else None),
        # The SAME predicate source_diverse_validate applies before it will report a cell at all.
        "passes_own_diversity_bar": share <= MAX_SOURCE_SHARE,
        "bar": MAX_SOURCE_SHARE,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--hunt", type=Path, default=ROOT / "wiki" / "gentamicin_rmt_bvbrc_hunt.json")
    ap.add_argument("--out", type=Path,
                    default=ROOT / "wiki" / f"rmt_source_concentration_{_date.today().isoformat()}.json")
    a = ap.parse_args()

    hits = json.loads(a.hunt.read_text(encoding="utf-8"))["all_hits"]
    by_org: dict[str, list[dict]] = collections.defaultdict(list)
    for h in hits:
        by_org[organism_of(h)].append(h)

    cells = {o: profile(hs) for o, hs in by_org.items()}
    cells = {o: p for o, p in cells.items() if p.get("n", 0) > 0}

    print(f"applying the project's own bar (largest-source share <= {MAX_SOURCE_SHARE:.0%})\n")
    for o, p in sorted(cells.items(), key=lambda kv: -kv[1]["n"]):
        flag = "PASSES" if p["passes_own_diversity_bar"] else "FAILS "
        print(f"  {flag} {o:<14} n={p['n']:<4} R={p['R']:<4} S={p['S']:<4} PPV={p['ppv']:.3f}  "
              f"sources={p['n_sources']} largest={p['largest_source_share']:.3f}")
        if p["S"] and p["susceptible_largest_source_share"] is not None:
            print(f"           susceptible carriers are {p['susceptible_largest_source_share']:.1%} "
                  f"from {p['susceptible_largest_source']}")
        if p["excluding_largest_source"]:
            e = p["excluding_largest_source"]
            print(f"           excluding it: n={e['n']} R={e['R']} S={e['S']} PPV={e['ppv']:.3f}")

    kleb = cells.get("klebsiella", {})
    if not kleb:
        print("\nno Klebsiella carriers in the artifact -- nothing to qualify", file=sys.stderr)
        return 2

    if kleb["passes_own_diversity_bar"]:
        verdict = "OVERCALL_SOURCE_DIVERSE"
        why = "the Klebsiella over-call rests on a source-diverse carrier set and stands as published."
    else:
        e = kleb["excluding_largest_source"]
        verdict = "OVERCALL_SINGLE_SOURCE_DOMINATED"
        why = (f"the Klebsiella over-call FAILS this project's own source-diversity bar "
               f"({kleb['largest_source_share']:.1%} from {kleb['largest_source']}, bar "
               f"{MAX_SOURCE_SHARE:.0%}); {kleb['susceptible_largest_source_share']:.1%} of the "
               f"SUSCEPTIBLE carriers that carry the entire finding come from that one study, and "
               f"excluding it the PPV is {e['ppv']:.3f} ({e['R']}R/{e['S']}S) -- indistinguishable from "
               "the validated E. coli scope. The over-call is established in ONE population, NOT for "
               "Klebsiella generally.")
    print(f"\nVERDICT: {verdict}\n  {why}")

    out = {
        "schema": "rmt-source-concentration-v1",
        "question": "does the published Klebsiella rmt over-call clear the source-diversity bar this "
                    "project enforces on its own validation cells?",
        "bar_imported_from": "scripts/source_diverse_validate.py::MAX_SOURCE_SHARE",
        "bar": MAX_SOURCE_SHARE,
        "cells": cells, "verdict": verdict, "why": why,
        "does_not_overturn": [
            "The organism_scope WARNING stays. An over-call warning that proves over-cautious costs far "
            "less than a missing one, and the safe direction is to keep it.",
            "The signal looks real INSIDE that study's population: it splits its own rmt carriers "
            "18R/63S rather than calling them all susceptible, and the pre-registered aac(3) control "
            "passed (99% R inside vs 83% elsewhere), so a global label inversion is excluded.",
            "The E. coli scope is untouched -- 12/12 here, 146/146 on NCBI-PD.",
        ],
        "what_it_does_overturn": (
            "the SCOPE of the claim. 'In Klebsiella the rescue over-calls on more than half of carriers' "
            "is an archive-level statement that a 66%-single-source carrier set cannot support."
        ),
        "honest_limits": [
            "Source unit is the curating publication (BV-BRC `pmid`), not BioProject -- comparable in "
            "spirit to the BioProject-based bar but not the identical unit.",
            "A single source cannot distinguish a real regional/population effect from something "
            "specific to that study's isolates. Both remain open.",
            "Excluding the largest source leaves n=41, which is a small comparator, not a refutation.",
        ],
    }
    a.out.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
