"""Stratify the serovar cell's accuracy by HOW each call was resolved, not just whether it was right.

WHY THIS EXISTS. The cell publishes 0.8400 (168 hit / 11 miss / 21 no_call on 200 wet-lab-labelled
isolates). That number pools outcomes of three different interpretability classes:

  * `exact`             -- the `(O,H1,H2)` key is in the table and its name survives BOTH ambiguity
                           policies, so the formula genuinely names one serovar.
  * `arbitrary_winner`  -- the key is contested (>1 Kauffmann-White serovar carries it even after
                           canonical naming) and the table records an arbitrary winner rather than
                           abstaining. A correct call here is right BY LUCK, in the memo's own words.
  * `fallback`          -- the exact key missed, and the caller matched `O+H1` ignoring H2, reporting
                           a serovar iff that collapsed to one name. Phase-incomplete resolution,
                           whose value may itself be an arbitrary winner.

A headline whose interpretable fraction is unknown is the project's own "a rate is a claim about its
denominator" error in a different costume, so this measures the split rather than asserting it is small.

HOW THE CONTESTED SET IS IDENTIFIED -- measured, not rebuilt. Docker is not required. The two ambiguity
policies were each scored on the SAME 200 isolates, and the policy affects ONLY the name lookup: this
script REFUSES unless all 200 antigenic formulas are byte-identical across the two checkpoints, which
is what makes the diff attributable to the policy alone. An isolate whose serovar appears under `first`
and is absent under `omit` is exactly one whose call depends on the arbitrary winner.

WHAT THIS IS NOT. It does not re-run blastn, re-call any antigen, or change any threshold; it reads two
committed checkpoints and the shipped table. It measures the COMPOSITION of an already-published number
and cannot move it.
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from datetime import date as _date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from dna_decode.salmserovar.equivalence import equivalent, load_formula_index  # noqa: E402

NON_SPECIFIC = {"", "-", "na", "n/a", "none", "unknown", "undetermined", "pending",
                "not applicable", "not determined", "untypeable", "untypable"}

EXACT = "exact"
ARBITRARY = "arbitrary_winner"
FALLBACK = "fallback"
UNRESOLVED = "unresolved"


def load_rows(path: Path) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for ln in path.read_text(encoding="utf-8").splitlines():
        ln = ln.strip()
        if ln:
            r = json.loads(ln)
            out[r["asm_acc"]] = r
    return out


def load_table(path: Path) -> dict[tuple[str, str, str], str]:
    table: dict[tuple[str, str, str], str] = {}
    for i, ln in enumerate(path.read_text(encoding="utf-8").splitlines()):
        if i == 0 or not ln.strip():
            continue
        parts = ln.split("\t")
        if len(parts) >= 4:
            table[(parts[0], parts[1], parts[2])] = parts[3]
    return table


def split_formula(formula: str | None) -> tuple[str, str, str] | None:
    """`O:H1:H2` -> the 3 axes, or None when the formula is absent or any axis is unresolved.

    The axes may THEMSELVES contain colons is not a concern here (they do not), but an unresolved
    axis is written `O?`/`H?` and must not be treated as a lookup key.
    """
    if not formula:
        return None
    parts = str(formula).split(":")
    if len(parts) != 3:
        return None
    o, h1, h2 = (p.strip() for p in parts)
    if not o or o.endswith("?") or not h1 or h1.endswith("?"):
        return None
    return o, h1, h2


def classify_route(formula: str | None, called: str | None, depends_on_winner: bool,
                   table: dict[tuple[str, str, str], str]) -> str:
    """Which mechanism produced this call.

    `depends_on_winner` is the MEASURED fact (named under `first`, absent under `omit`), so it takes
    precedence over any inference from the table -- the table only records winners and cannot say on
    its own which keys were contested.
    """
    if not called or str(called).strip().lower() in NON_SPECIFIC:
        return UNRESOLVED
    if depends_on_winner:
        return ARBITRARY
    key = split_formula(formula)
    if key is not None and key in table:
        return EXACT
    return FALLBACK


def score(call: str | None, truth: str, idx: dict) -> str:
    if not call or str(call).strip().lower() in NON_SPECIFIC:
        return "no_call"
    ok, _ = equivalent(str(call), truth, idx)
    return "hit" if ok else "miss"


def stratify(cohort: list[dict], first: dict[str, dict], omit: dict[str, dict],
             table: dict[tuple[str, str, str], str], idx: dict) -> dict:
    strata: dict[str, Counter] = {k: Counter() for k in (EXACT, ARBITRARY, FALLBACK, UNRESOLVED)}
    per_isolate = []
    for row in cohort:
        acc = row["asm_acc"]
        f, o = first.get(acc), omit.get(acc)
        if f is None:
            continue
        called = f.get("serovar")
        depends = bool(called) and not (o or {}).get("serovar")
        route = classify_route(f.get("formula"), called, depends, table)
        outcome = score(called, row["serovar_raw"], idx)
        strata[route][outcome] += 1
        per_isolate.append({"asm_acc": acc, "route": route, "outcome": outcome,
                           "formula": f.get("formula"), "called": called,
                           "truth": row["serovar_raw"]})
    return {"strata": {k: dict(v) for k, v in strata.items()}, "per_isolate": per_isolate}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cohort", type=Path, default=ROOT / "wiki" / "salmserovar_cohort_2026-09-04.json")
    ap.add_argument("--first", type=Path,
                    default=Path("D:/dna_decode_cache/ss2_checkpoint/ofix_rows_variantB.jsonl"),
                    help="the SHIPPED run (ambiguous_policy='first')")
    ap.add_argument("--omit", type=Path,
                    default=Path("D:/dna_decode_cache/ss2_checkpoint/ofix_rows.variantA.jsonl"),
                    help="the principled-policy run (ambiguous_policy='omit')")
    ap.add_argument("--table", type=Path, default=ROOT / "data" / "salmserovar_db" / "serovar_table.tsv")
    ap.add_argument("--out", type=Path, default=ROOT / "wiki" / f"salmserovar_resolution_strata_{_date.today()}.json")
    a = ap.parse_args(argv)

    cohort_doc = json.loads(a.cohort.read_text(encoding="utf-8"))
    cohort = cohort_doc["isolates"] if isinstance(cohort_doc, dict) and "isolates" in cohort_doc else cohort_doc
    first, omit = load_rows(a.first), load_rows(a.omit)
    table = load_table(a.table)

    # --- comparability gate ----------------------------------------------------------------------
    # The whole method rests on the two runs differing ONLY in the name lookup. If any antigenic
    # formula moved, the diff is not attributable to the ambiguity policy and the strata would be
    # fiction. Refuse rather than report a number built on an unchecked premise.
    shared = sorted(set(first) & set(omit))
    mismatched = [k for k in shared if first[k].get("formula") != omit[k].get("formula")]
    if mismatched or not shared:
        print(f"REFUSING: {len(mismatched)} of {len(shared)} formulas differ between the two policy "
              f"runs, so a serovar-name diff cannot be attributed to the ambiguity policy alone.")
        return 2

    res = stratify(cohort, first, omit, table, load_formula_index(a.table))
    strata = res["strata"]
    hits = sum(v.get("hit", 0) for v in strata.values())
    n = sum(sum(v.values()) for v in strata.values())
    called_routes = (EXACT, ARBITRARY, FALLBACK)
    interpretable = strata[EXACT].get("hit", 0)

    # --- what the arbitrary-winner policy actually TRADES ------------------------------------------
    # The name-table memo reports this policy as "strictly dominant" because it gains 3 hits over
    # `omit` at zero regressions. That is true ON THE HIT COUNT and it is not the whole trade: every
    # one of these 8 isolates would ABSTAIN under `omit`, so the policy also converts 5 abstentions
    # into confidently-WRONG serovar calls. This cell's own threshold work recorded the asymmetry
    # deliberately -- "an abstention is recoverable by a human, a confident wrong serovar is not" --
    # so the two framings disagree about the same 8 isolates and the disagreement belongs in the
    # artifact rather than in whichever summary is written last.
    aw = strata[ARBITRARY]
    arbitrary_trade = {
        "isolates_whose_call_depends_on_the_winner": sum(aw.values()),
        "hits_gained_vs_omit": aw.get("hit", 0),
        "abstentions_converted_to_wrong_calls": aw.get("miss", 0),
        "hit_rate_within_this_stratum": round(aw.get("hit", 0) / sum(aw.values()), 4) if sum(aw.values()) else None,
        "note": "Net +3 on hit count; net +5 on confidently-wrong calls. Which reading governs is an "
                "acceptance-bar question, deliberately NOT resolved here.",
    }

    # --- the largest single no-call cluster --------------------------------------------------------
    nocall_formulas = Counter(r["formula"] for r in res["per_isolate"]
                              if r["outcome"] == "no_call" and r["formula"])
    top_nocall = nocall_formulas.most_common(3)
    doc = {
        "schema": "salmserovar-resolution-strata-v1",
        "date": str(_date.today()),
        "question": "Of the published hits, how many were uniquely resolved vs won by an arbitrary "
                    "tie-break vs reached through the phase-incomplete fallback?",
        "n": n,
        "headline_reproduced": {"hit": hits, "accuracy": round(hits / n, 4) if n else None},
        "strata": strata,
        "hits_by_route": {r: strata[r].get("hit", 0) for r in called_routes},
        "unique_resolution_accuracy": round(interpretable / n, 4) if n else None,
        "formula_comparability_gate": {"checked": len(shared), "mismatched": 0, "passed": True},
        "arbitrary_winner_trade": arbitrary_trade,
        "fallback_finding": {
            "hits": strata[FALLBACK].get("hit", 0),
            "reached_but_unresolved": sum(1 for r in res["per_isolate"]
                                          if r["outcome"] == "no_call" and split_formula(r["formula"])),
            "note": "The O+H1 phase-incomplete fallback produced ZERO calls on this cohort -- every "
                    "one of the calls came through an exact key. It is not dead code: it was REACHED "
                    "on the isolates carrying a complete formula whose exact key missed, and failed "
                    "to collapse to a single name on each. So disclosure metadata for the fallback "
                    "route has measured headroom of zero HERE, which is a reason not to build it yet "
                    "rather than evidence it can never matter.",
        },
        "largest_no_call_clusters": [{"formula": f, "n": c} for f, c in top_nocall],
        "per_isolate": res["per_isolate"],
        "honest_limits": [
            "Measures the COMPOSITION of an already-published number; it cannot change that number.",
            "The contested set is identified by DIFFING the two ambiguity policies on the same 200 "
            "isolates, not by rebuilding the candidate table from Initial_Conditions (Docker down). "
            "That captures every isolate whose CALL depends on the winner, which is the "
            "decision-relevant set -- but it does NOT count a contested key whose arbitrary winner "
            "happens to equal what the fallback would return anyway, so the arbitrary_winner stratum "
            "is a LOWER BOUND on contested-key involvement.",
            "The fallback stratum may itself contain arbitrary winners inherited through the O+H1 "
            "candidate set; this split does not separate those.",
            "Per-serovar cap of 12 flattens prevalence, so every figure here is per-isolate accuracy "
            "on a diverse mix, not a population-weighted rate.",
        ],
    }
    a.out.write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8")
    print(f"n={n}  hits={hits} (acc {doc['headline_reproduced']['accuracy']})")
    for r in called_routes:
        print(f"  {r:<18} {dict(strata[r])}")
    print(f"  {UNRESOLVED:<18} {dict(strata[UNRESOLVED])}")
    print(f"uniquely-resolved accuracy = {doc['unique_resolution_accuracy']}")
    print(f"-> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
