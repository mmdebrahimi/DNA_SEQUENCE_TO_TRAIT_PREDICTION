"""Flowering cell vs MEASURED AraPheno flowering time — an honest SPOT-CHECK (not a score).

The flowering cell (dna_decode/organism_rules/arabidopsis_flowering.py) is faithful-to-literature: it applies
published FRI/FLC functional-allele assignments. This confronts it with REAL MEASURED data — the free AraPheno
DTF1 phenotype (days to first flowering buds, 1001 Genomes accessions, CSV/free).

**WHY THIS IS A SPOT-CHECK, NOT A SCORED CELL (the load-bearing honesty rail):**
A real SCORE needs a per-accession FRI/FLC allele mapping for a powered cohort. That mapping EXISTS
(Zhang 2020 Plant J Table S1 = 1016 accessions x FRI allele; Shindo 2005 = 176 accessions with FRI haplotype +
flowering time + FLC level; Werner 2005 = 145 accessions genotyped for the 2 common FRI lesions) but EVERY
free route to it is walled: PMC -> CAPTCHA, OUP -> abstract-only, bioRxiv supplement -> figures only (no
Table S1), Wiley -> HTTP 402 Payment Required (a MONEY gate — not paid), and no CSV/repo mirror exists.
So this joins ONLY the handful of accessions whose FRI/FLC status the cell's own curated catalog already
carries from the literature => N is tiny and NOT a cohort. It is an ANECDOTE-level reality check that can
FALSIFY, not a validation that can confirm.

Threshold discipline (R2 — derive, don't assert): the early/late split is the MEDIAN of the measured DTF1
distribution over all accessions in the file, not a number I picked.

Run: uv run python scripts/flowering_arapheno_spotcheck.py [--dtf1 <cached csv>]
Free, offline once the CSV is cached. Frozen decoder surface untouched.
"""
from __future__ import annotations

import argparse
import csv
import json
import statistics
import sys
from datetime import date as _date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

DTF1_URL = "https://arapheno.1001genomes.org/rest/phenotype/703/values.csv"

# Accessions whose FRI/FLC status the cell's curated catalog carries, mapped to their AraPheno name.
# (fri_allele, flc_allele) as named in dna_decode.organism_rules.arabidopsis_flowering.
CATALOG_ACCESSIONS: dict[str, tuple[str, str]] = {
    "Col-0": ("Col", "Col"),        # FRI-LoF (16-bp del), strong FLC -> predict early
    "Van-0": ("Sf-2", "Van-0"),     # functional FRI, FLC nonsense -> predict early (the FLC route)
    "Bur-0": ("Sf-2", "Bur-0"),     # functional FRI, FLC aberrant-splice(null-behaving) -> predict early
    "Sf-2": ("Sf-2", "Col"),        # functional FRI + strong FLC -> predict late (winter annual)
    "Wil-2": ("Wil-2", "Col"),      # FRI-LoF (substitution), strong FLC -> predict early
}


def load_dtf1(path: str | None) -> list[dict]:
    if path and Path(path).exists():
        text = Path(path).read_text(encoding="utf-8", errors="replace")
    else:
        import urllib.request
        with urllib.request.urlopen(DTF1_URL, timeout=60) as r:   # free, public, read-only
            text = r.read().decode("utf-8", errors="replace")
    return list(csv.DictReader(text.splitlines()))


def main() -> int:
    ap = argparse.ArgumentParser(description="Flowering-cell spot-check vs measured AraPheno DTF1.")
    ap.add_argument("--dtf1", help="cached AraPheno DTF1 values.csv (else fetched from the free API)")
    args = ap.parse_args()

    from dna_decode.organism_rules.arabidopsis_flowering import call_flowering_habit

    rows = load_dtf1(args.dtf1)
    measured: dict[str, float] = {}
    for r in rows:
        try:
            measured[r["accession_name"].strip()] = float(r["phenotype_value"])
        except (KeyError, TypeError, ValueError):
            continue
    if not measured:
        print("error: no DTF1 values parsed", file=sys.stderr)
        return 2

    # R2: the early/late split is DERIVED from the measured distribution, not asserted.
    med = statistics.median(measured.values())

    checks = []
    for acc, (fri, flc) in CATALOG_ACCESSIONS.items():
        if acc not in measured:
            checks.append({"accession": acc, "status": "NOT_IN_DTF1"})
            continue
        call = call_flowering_habit(fri, flc)
        days = measured[acc]
        obs = "late" if days > med else "early"
        pred = "late" if call.habit == "winter_annual_late" else (
            "early" if call.habit == "summer_annual_early" else "abstain")
        checks.append({
            "accession": acc, "fri_allele": fri, "flc_allele": flc,
            "predicted_habit": call.habit, "predicted": pred, "confidence": call.confidence,
            "measured_days_to_flower": days, "observed": obs,
            "agree": (pred == obs) if pred != "abstain" else None,
        })

    scored = [c for c in checks if c.get("agree") is not None]
    n_agree = sum(1 for c in scored if c["agree"])
    misses = [c for c in scored if not c["agree"]]

    # NULL-BASELINE GATE (the threshold-vs-null-baseline discipline — compute the trivial baseline BEFORE
    # believing any agreement count). A binary call cannot be scored on a set with only ONE observed class:
    # a constant predictor attains the majority rate by construction.
    obs_classes = {c["observed"] for c in scored}
    maj = max(("early", "late"), key=lambda k: sum(1 for c in scored if c["observed"] == k))
    n_null = sum(1 for c in scored if c["observed"] == maj)
    degenerate = len(obs_classes) < 2
    verdict = "DEGENERATE_NO_LABEL_VARIATION" if degenerate else (
        "SPOT_CHECK_BEATS_NULL" if n_agree > n_null else "SPOT_CHECK_NOT_ABOVE_NULL")

    result = {
        "cell": "arabidopsis_flowering_arapheno_spotcheck",
        "status": "SPOT_CHECK_NOT_SCORED",
        "verdict": verdict,
        "null_baseline": {
            "constant_predictor": maj, "n_null_agree": n_null, "n_cell_agree": n_agree, "n": len(scored),
            "observed_classes_present": sorted(obs_classes),
            "reading": ("the joined set carries only ONE observed class -> a constant predictor attains "
                        f"{n_null}/{len(scored)} by construction and the cell's {n_agree}/{len(scored)} is "
                        "NOT evidence of anything. Reporting an agreement RATE here would be meaningless."
                        if degenerate else
                        f"constant-'{maj}' null attains {n_null}/{len(scored)}; cell attains {n_agree}"),
        },
        "why_not_scored": (
            "a real score needs a per-accession FRI/FLC allele mapping for a powered cohort. It EXISTS "
            "(Zhang 2020 Plant J Table S1 = 1016 accessions; Shindo 2005 = 176; Werner 2005 = 145) but every "
            "FREE route is walled: PMC=CAPTCHA, OUP=abstract-only, bioRxiv-supp=figures-only, "
            "Wiley=HTTP 402 Payment Required (MONEY gate, not paid), no CSV/repo mirror. So only the "
            "catalog's own literature-known accessions could be joined => N is an anecdote, not a cohort."),
        "measured_source": {"db": "AraPheno (1001 Genomes)", "phenotype": "DTF1 (days to visible flowering "
                            "buds from sowing)", "url": DTF1_URL, "n_accessions_in_file": len(measured)},
        "threshold_discipline": {"rule": "early/late split = MEDIAN of the measured DTF1 distribution "
                                 "(derived, not asserted)", "median_days": round(med, 3)},
        "n_joined": len(scored), "n_agree": n_agree,
        "checks": checks,
        "misses": misses,
    }
    out = REPO / "wiki" / f"arabidopsis_flowering_spotcheck_{_date.today().isoformat()}.json"
    out.write_text(json.dumps(result, indent=2), encoding="utf-8")

    print(f"[spotcheck] AraPheno DTF1: {len(measured)} accessions | median split {med:.2f} days")
    print(f"[spotcheck] joined {len(scored)} catalog accessions (N is an ANECDOTE, not a cohort)")
    print(f"{'accession':10}{'FRI':8}{'FLC':10}{'pred':7}{'obs':7}{'days':>8}  agree")
    for c in checks:
        if c.get("status") == "NOT_IN_DTF1":
            print(f"{c['accession']:10}{'-':8}{'-':10}{'-':7}{'-':7}{'-':>8}  not-in-DTF1")
        else:
            print(f"{c['accession']:10}{c['fri_allele']:8}{c['flc_allele']:10}{c['predicted']:7}"
                  f"{c['observed']:7}{c['measured_days_to_flower']:>8.2f}  {c['agree']}")
    nb = result["null_baseline"]
    print(f"[spotcheck] cell agree {n_agree}/{len(scored)} | NULL constant-'{maj}' agree "
          f"{n_null}/{len(scored)} | observed classes present: {sorted(obs_classes)}")
    for m in misses:
        print(f"[spotcheck] MISS: {m['accession']} predicted {m['predicted']} but measured "
              f"{m['measured_days_to_flower']} days ({m['observed']})")
    print(f"[spotcheck] VERDICT: {verdict}")
    print(f"[spotcheck] {nb['reading']}")
    print(f"artifact -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
