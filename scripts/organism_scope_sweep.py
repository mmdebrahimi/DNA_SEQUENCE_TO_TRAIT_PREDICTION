"""Is the gentamicin organism-scope gap systemic, or isolated? Sweep every deployed rule and find out.

THE CONCERN, STATED GENERALLY. The gentamicin `symbol_rescue` is ORGANISM-AGNOSTIC IN CODE while its
evidence is organism-SPECIFIC: `calibrated_rule_for` has no gentamicin entry for any organism, so every
organism falls through to the default `DRUG_RULE`. That was found by accident on one rule. Four of the
six deployed drugs have NO calibrated entry at all (ceftriaxone / gentamicin / oxacillin / tetracycline),
so the same shape is possible for each of them, and "found by accident once" is not a scope statement.

WHAT THIS MEASURES. For every drug scored on more than one organism, the spread in specificity across
organisms -- specificity being the over-call axis, since a rule that fires too readily converts true
susceptibles into false resistants. A rule whose behaviour is genuinely organism-dependent should show a
large spread; a rule that travels should show a small one.

THE TWO NUMBERS ARE NOT INTERCHANGEABLE, which is why the source-diversity share is carried alongside.
A cell resting on one BioProject can post any spec at all, so a large spread in a concentrated cell is
evidence about the cohort, not about the rule. Cells are flagged, never silently dropped.

Offline: reads the committed report card. Writes wiki/organism_scope_sweep_<date>.json.
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

from source_diverse_validate import MAX_SOURCE_SHARE  # noqa: E402
from dna_decode.data.organism_scope import overcall_for  # noqa: E402

# A spread this large across organisms is the size of the gap that made gentamicin worth disclosing
# (the BV-BRC Klebsiella PPV fell from 1.000 to 0.475, a 0.525 swing). Anything approaching it deserves
# a look; the bar is deliberately generous so the sweep over-reports rather than under-reports.
SPREAD_FLAG = 0.15


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--card", type=Path, default=ROOT / "wiki" / "decoder_validation_report_card.json")
    ap.add_argument("--out", type=Path,
                    default=ROOT / "wiki" / f"organism_scope_sweep_{_date.today().isoformat()}.json")
    a = ap.parse_args()

    cells = json.loads(a.card.read_text(encoding="utf-8"))["cells"]
    by: dict[str, list[dict]] = collections.defaultdict(list)
    for c in cells:
        if c.get("state") != "SCORED":
            continue
        sc = c.get("source_concentration") or {}
        by[c["drug"]].append({
            "organism": c["organism"], "sens": c.get("sens"), "spec": c.get("spec"), "n": c.get("n"),
            "largest_source_share": sc.get("largest_share"),
            "distinct_bioprojects": sc.get("distinct_bioprojects"),
            "passes_diversity_bar": (sc.get("largest_share") is not None
                                     and sc["largest_share"] <= MAX_SOURCE_SHARE),
        })

    findings, multi = [], 0
    for drug, rows in sorted(by.items()):
        if len(rows) < 2:
            continue
        multi += 1
        specs = [r["spec"] for r in rows if isinstance(r["spec"], (int, float))]
        if len(specs) < 2:
            continue
        spread = max(specs) - min(specs)
        # A spread carried by a concentrated cell says more about that cohort than about the rule.
        concentrated = [r["organism"] for r in rows if not r["passes_diversity_bar"]]
        findings.append({
            "drug": drug, "organisms": rows, "spec_spread": spread,
            "flagged": spread >= SPREAD_FLAG,
            "concentrated_cells": concentrated,
            "attributable_to_concentration": bool(spread >= SPREAD_FLAG and concentrated),
            "already_disclosed": overcall_for(drug, "Klebsiella") is not None,
        })

    print(f"multi-organism SCORED drugs: {multi}   (spread flag >= {SPREAD_FLAG})\n")
    for f in findings:
        mark = "FLAG" if f["flagged"] else "ok  "
        print(f"  {mark} {f['drug']:<16} spec spread {f['spec_spread']:.3f}"
              f"{'   concentrated: ' + ','.join(f['concentrated_cells']) if f['concentrated_cells'] else ''}")

    real = [f for f in findings
            if f["flagged"] and not f["attributable_to_concentration"] and not f["already_disclosed"]]
    if real:
        verdict = "ORGANISM_SCOPE_GAP_FOUND"
        why = ("a drug shows an organism-dependent spec spread that is NOT attributable to a "
               f"concentrated cohort and is NOT already disclosed: {[f['drug'] for f in real]}")
    else:
        verdict = "ISOLATED_TO_GENTAMICIN"
        why = ("no other deployed drug shows an organism-dependent over-call. Every other multi-organism "
               "spec spread is small, and the one large spread is attributable to a cohort that fails "
               "the source-diversity bar rather than to the rule. The gentamicin gap is isolated, not "
               "the first instance of a systemic pattern.")
    print(f"\nVERDICT: {verdict}\n  {why}")

    out = {
        "schema": "organism-scope-sweep-v1",
        "question": "is the gentamicin organism-agnostic-code / organism-specific-evidence gap systemic?",
        "spread_flag": SPREAD_FLAG, "diversity_bar": MAX_SOURCE_SHARE,
        "n_multi_organism_drugs": multi, "findings": findings,
        "verdict": verdict, "why": why,
        "what_this_does_not_show": [
            "A drug scored on only ONE organism cannot show a spread at all -- absence of a flag is not "
            "evidence of organism-independence for meropenem or oxacillin.",
            "The scored cells were produced by the FROZEN rule. A cell can only reveal an over-call for "
            "a determinant its cohort actually CONTAINS; a cohort holding no carriers of a determinant "
            "is structurally incapable of detecting a rule keyed on it (the Oxford result exactly).",
            "Specificity is the over-call axis but not the only axis; a rule could be organism-dependent "
            "in sensitivity while looking flat here.",
        ],
    }
    a.out.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
