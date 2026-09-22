"""The Stanford PhenoSense fold-change tables are RIGHT-CENSORED at 100, and nothing here knew it.

WHAT THIS IS. A measurement of assay censoring in the project's one free, independent, isolate-level,
CONTINUOUS wet-lab label -- and of what that censoring can do to three PUBLISHED headline numbers.

WHY IT WAS INVISIBLE. The repo ships `assay_degeneracy()` (scripts/forward_inverse_roundtrip.py) which
REJECTS a continuous assay whose mode-share exceeds 25% -- the guard that excluded CcdB from the
forward/inverse sweep after it posted the sweep's best number purely because 79.3% of its variants sat
at the assay ceiling. That guard lives in the forward/inverse arm and has never been pointed at the HIV
arm's fold-change tables, which live in a different part of the tree.

THE EXISTING GUARD IS VACUOUS, WHICH IS ITS OWN FINDING. `hiv_nnrti_mutant_catalog.py` and
`hiv_epistasis.py` both state that censored folds ("<" / ">" prefixed) are kept at their numeric bound.
There are ZERO operator-prefixed values in any drug column of any of the four datasets. The guard
protects against a form of censoring that does not occur here, while the censoring that DOES occur --
a bare `100` with nothing above it -- passes straight through. Same shape as the self-corrected vacuous
ResFinder POINT-row filter: a filter that removes nothing is not a control.

WHAT IT DOES **NOT** CLAIM, and this scoping is load-bearing:

  * It is NOT a decoder change. No v0.1 mutant catalog is deployed -- `dna_decode/data/hiv_amr.py`
    routes NRTI/PI/INSTI through POSITION-based classes (`NRTI_MAJOR_POSITIONS`, `PI_CLASS`,
    `INSTI_CLASS`). This is a methodological audit of published numbers.
  * It does NOT reopen the declined NNRTI curation. That verdict was scored at
    `ILLUSTRATIVE_FOLD_CUTOFF = 3.0`, and the ceiling sits at 100 >> 3, so every censored observation
    is R under either reading: the R/S LABELS, and therefore the sens/spec that drove the verdict, are
    untouched by censoring. Censoring bites only in the COEFFICIENT-ESTIMATION step.
  * It does NOT re-explain the NNRTI OLS dropping Y181C. `hiv_amr.py` SHIPS Y181C/I/A, so that is a
    fact about a derived catalog that was never deployed; and the metric that drove the verdict is
    blind-spot recovery, defined over isolates carrying NO catalogued DRM, which Y181C carriers are
    not. Censoring and the co-occurrence account are not rivals -- if two DRMs both pin fold at the
    ceiling, censoring REMOVES the residual variance OLS needs to separate them, which AMPLIFIES the
    co-occurrence story rather than replacing it.

WHAT IT DOES CLAIM. The NRTI / PI / INSTI v0.1 catalogs were each selected by thresholding a
multivariate-OLS log10-fold coefficient at >= log10(1.5), fit on a response right-censored at 100 on
drugs where up to 45% of observations sit at the ceiling. Right-censoring attenuates coefficients
toward zero for the STRONGEST effects specifically -- exactly the majors such a catalog exists to
find -- so a coefficient threshold can systematically exclude true majors. Those three gain figures
are published headlines. Whether they move under a censoring-aware (Tobit) refit is a separate,
decidable question this audit does not answer.

Offline, read-only, seconds. Frozen AMR surface untouched.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from datetime import date as _date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

# The bars are IMPORTED from the shipped screen, never restated -- a restated constant drifts from the
# thing it is supposed to mirror, and then the audit silently measures a different rule.
from scripts.forward_inverse_roundtrip import (  # noqa: E402
    MAX_MODE_SHARE, MIN_DISTINCT_VALUES, assay_degeneracy,
)

DATASETS = ("NNRTI", "NRTI", "PI", "INI")
MIN_N = 30                    # below this a mode-share is noise, not a distribution
MIN_PLAUSIBLE_MAX = 1.5       # a fold column must reach at least this; excludes index/flag columns
ID_UNIQUENESS = 0.9           # all-integer AND this fraction distinct => an identifier, not a measurement


def fold_columns(rows: list[dict]) -> dict[str, list[float]]:
    """Every column that parses as a fold-change measurement, with its values.

    Derived from the data, NOT hand-listed: a hand-enumerated drug list is the documented
    `hardcoded_exclusion_list_undercovers` failure, and it would silently skip a drug added upstream.
    A column qualifies iff it is fully numeric over >= MIN_N non-missing values, reaches
    MIN_PLAUSIBLE_MAX, and is not IDENTIFIER-shaped.

    The identifier rule matters and was found by running this: `SeqID` is all-integer with 2,272
    distinct values in 2,272 rows and a maximum of 789,388. Admitted as a "fold column" it became the
    dataset's maximum, so the shared-upper-bound statistic was computed against an accession number
    and reported `None`. An id column is all-integer AND near-unique; a fold column is neither.
    """
    out: dict[str, list[float]] = {}
    for col in (rows[0] if rows else {}):
        vals = [r[col] for r in rows if r.get(col) not in ("", "NA", "-", None)]
        if len(vals) < MIN_N:
            continue
        try:
            nums = [float(v) for v in vals]
        except ValueError:
            continue
        if max(nums) < MIN_PLAUSIBLE_MAX:
            continue
        if all(float(x).is_integer() for x in nums) and len(set(nums)) / len(nums) > ID_UNIQUENESS:
            continue                       # identifier-shaped: all-integer and near-unique
        out[col] = nums
    return out


def censoring_profile(values: list[float]) -> dict:
    """Is this column right-censored, and where?

    A CEILING is not merely 'the maximum value'. It is a maximum that a large share of observations
    share while NOTHING exceeds it -- a pile-up at the top with an empty tail above. Reporting the max
    alone would call every column censored; reporting mode-share alone would flag a genuine biological
    mode. Both conditions must hold, and `n_above_max` is emitted so a reader can check the second.
    """
    n = len(values)
    counts = Counter(values)
    mx = max(values)
    at_max = counts[mx]
    return {
        "n": n,
        "max": mx,
        "n_at_max": at_max,
        "share_at_max": round(at_max / n, 4),
        "n_strictly_above_max": 0,          # true by construction; stated so the claim is legible
        "n_distinct": len(counts),
        "looks_right_censored": bool(at_max / n >= 0.05 and len(counts) >= MIN_DISTINCT_VALUES),
    }


def audit_dataset(path: Path) -> dict:
    with open(path, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    cols = fold_columns(rows)
    per_col = {}
    for col, nums in cols.items():
        deg = assay_degeneracy(nums)
        cen = censoring_profile(nums)
        per_col[col] = {
            "n": deg["n"],
            "n_distinct_values": deg["n_distinct_values"],
            "mode_value": deg["mode_value"],
            "mode_share": deg["mode_share"],
            "fails_shipped_degeneracy_bar": deg["degenerate"],
            "censoring": cen,
            "mode_is_the_ceiling": bool(deg["mode_value"] == cen["max"]),
        }
    maxima = {c["censoring"]["max"] for c in per_col.values()}
    # The right statistic is a shared UPPER BOUND, not a shared maximum. Some drugs (TDF 51, BIC 66)
    # simply never reach the cap in this cohort, so "every column maxes at the same value" is FALSE
    # while "no value anywhere exceeds the cap, and every column that reaches it piles up there" is
    # true and is the actual claim. Demanding equality reported `None` and read as a refutation.
    cap = max(maxima) if maxima else None
    reach = {c: v for c, v in per_col.items() if v["censoring"]["max"] == cap}
    return {
        "n_rows": len(rows),
        "n_fold_columns": len(per_col),
        "upper_bound": cap,
        "n_columns_reaching_the_bound": len(reach),
        "n_columns_below_the_bound": len(per_col) - len(reach),
        "columns_below_the_bound": {c: v["censoring"]["max"]
                                    for c, v in per_col.items() if c not in reach},
        "no_column_exceeds_the_bound": True,      # by construction of `cap`; restated for legibility
        "min_pileup_share_among_reaching_columns": (
            round(min(v["censoring"]["share_at_max"] for v in reach.values()), 4) if reach else None),
        "columns_failing_the_bar": sorted(
            c for c, v in per_col.items() if v["fails_shipped_degeneracy_bar"]),
        "per_column": per_col,
    }


def operator_prefixed_count(path: Path) -> int:
    """How many values carry a '<' / '>' operator -- the form the EXISTING guard handles.

    This is the non-vacuity check on that guard. If this is 0 everywhere, the guard cannot fire and
    the censoring that is actually present is unhandled.
    """
    with open(path, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return sum(1 for r in rows for v in r.values()
               if isinstance(v, str) and v.strip()[:1] in ("<", ">"))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", type=Path, default=ROOT / "data" / "raw" / "hiv")
    ap.add_argument("--out", type=Path,
                    default=ROOT / "wiki" / f"hiv_fold_censoring_audit_{_date.today()}.json")
    a = ap.parse_args(argv)

    per_ds, missing, op_counts = {}, [], {}
    for ds in DATASETS:
        p = a.data_dir / f"{ds}_DataSet.txt"
        if not p.exists():
            missing.append(ds)
            continue
        per_ds[ds] = audit_dataset(p)
        op_counts[ds] = operator_prefixed_count(p)

    if not per_ds:
        print("no HIV datasets found (they are gitignored); nothing to audit")
        return 2

    bounds = {v["upper_bound"] for v in per_ds.values()}
    failing = sorted({f"{ds}:{c}" for ds, v in per_ds.items() for c in v["columns_failing_the_bar"]})

    doc = {
        "schema": "hiv-fold-censoring-audit-v1",
        "date": str(_date.today()),
        "question": ("is the project's one free continuous wet-lab label right-censored, and does it "
                     "fail the repo's OWN shipped assay-degeneracy bar?"),
        "bars_imported_not_restated": {
            "max_mode_share": MAX_MODE_SHARE, "min_distinct_values": MIN_DISTINCT_VALUES,
            "source": "scripts/forward_inverse_roundtrip.py",
        },
        "datasets_missing_gitignored": missing,
        "results": per_ds,
        "headline": {
            "one_upper_bound_across_all_four_datasets": len(bounds) == 1 and None not in bounds,
            "upper_bound": (sorted(bounds)[0] if len(bounds) == 1 else None),
            "no_value_anywhere_exceeds_it": all(v["no_column_exceeds_the_bound"]
                                                for v in per_ds.values()),
            "columns_reaching_it": sum(v["n_columns_reaching_the_bound"] for v in per_ds.values()),
            "columns_below_it": sum(v["n_columns_below_the_bound"] for v in per_ds.values()),
            "columns_failing_the_shipped_bar": failing,
            "reading": ("a shared UPPER BOUND that NOTHING exceeds anywhere. Two drugs (TDF 51, BIC "
                        "66) never reach it in this cohort -- not evidence against the cap, just "
                        "drugs whose observed folds stay below it. The pile-up AT the bound varies "
                        "widely among the columns that do reach it (from 0.1% to 45%), so 'every "
                        "reaching column piles up' would be FALSE; what is shared is the bound, and "
                        "the heavy pile-ups are where the degeneracy bar fails."),
        },
        "existing_guard_is_vacuous": {
            "guard": ("hiv_nnrti_mutant_catalog.py / hiv_epistasis.py state that censored folds "
                      "('<'/'>') are kept at their numeric bound"),
            "operator_prefixed_values_found": op_counts,
            "verdict": ("VACUOUS -- the guard handles a censoring form that does not occur, while the "
                        "censoring that does occur (a bare ceiling value) passes through unflagged"
                        if not any(op_counts.values()) else "guard is live; re-read before citing"),
        },
        "what_this_does_not_claim": [
            "NOT a decoder change: no v0.1 mutant catalog is deployed; hiv_amr.py routes NRTI/PI/INSTI "
            "through POSITION-based classes.",
            "Does NOT reopen the declined NNRTI curation: that was scored at a fold cutoff of 3.0 and "
            "the ceiling is 100, so every censored observation is R either way -- the R/S labels and "
            "the sens/spec behind the verdict are unaffected.",
            "Does NOT re-explain the OLS dropping Y181C: hiv_amr.py ships Y181C/I/A, and censoring "
            "AMPLIFIES rather than replaces the co-occurrence account (a shared ceiling removes the "
            "residual variance needed to separate two co-occurring DRMs).",
            "Does NOT show that any published v0.1 gain figure is wrong. It shows the response those "
            "figures were fit on is censored; whether the numbers MOVE under a Tobit refit is a "
            "separate decidable question this audit does not answer.",
        ],
        "why_it_matters": (
            "the NRTI (+0.06..+0.14 on 5/6), PI (+0.056 mean, 8/8) and INSTI (+0.087 mean, 5/5) v0.1 "
            "catalogs were each selected by thresholding a multivariate-OLS log10-fold coefficient at "
            ">= log10(1.5). Right-censoring attenuates coefficients toward zero for the STRONGEST "
            "effects specifically -- the majors such a catalog exists to find -- so a coefficient "
            "threshold fit on a censored response can systematically exclude true majors."),
        "honest_limits": [
            "The ceiling is INFERRED from a pile-up with an empty tail above it, not read off assay "
            "documentation. It is consistent with a reporting cap; the audit does not prove the assay's "
            "dynamic range ends there.",
            "Attenuation is a property of fitting OLS to a censored response; this audit does NOT "
            "measure how much any particular coefficient moved.",
            "The HIV datasets are gitignored, so this cannot run in a clean checkout -- it degrades to "
            "exit 2 rather than reporting a vacuous clean result.",
        ],
    }
    a.out.write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8")

    h = doc["headline"]
    print(f"upper bound across all 4 datasets: {h['upper_bound']} "
          f"| nothing exceeds it: {h['no_value_anywhere_exceeds_it']} "
          f"| {h['columns_reaching_it']} cols reach it, {h['columns_below_it']} stay below")
    for ds, v in per_ds.items():
        worst = max(v["per_column"].items(), key=lambda kv: kv[1]["mode_share"])
        print(f"  {ds:<6} {v['n_fold_columns']:>2} fold cols | worst {worst[0]:<5} "
              f"mode-share {worst[1]['mode_share']:.3f} at {worst[1]['mode_value']} "
              f"| fails bar: {v['columns_failing_the_bar'] or 'none'}")
    print(f"\noperator-prefixed values (what the EXISTING guard handles): {op_counts}")
    print(f"-> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
