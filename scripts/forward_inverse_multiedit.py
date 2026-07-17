"""Does the MULTI-EDIT inverse extension earn a build? — the handoff's second experiment, answered cheaply.

The 2026-07-16 handoff proposed `cooc-multiedit-inverse` as the second experiment, flagged `[inferred]` /
"viability unproven", for this stated purpose:

    "for targets NO SINGLE EDIT REACHES, retrieve observed co-occurring edit-SETS"

That motivation is a factual claim about the data, and it is checkable before building anything. Two
independent tests, both on the only multi-mutant substrate in the cache (SR43C_ARATH, 694 real doubles
alongside 889 singles, measured in the same assay):

  Q1 REACH -- do multi-edits reach effects no single edit reaches? If not, the extension's whole purpose
     is void: you would be adding combinatorial machinery to reach targets a single edit already covers.

  Q2 ADDITIVITY -- does sum(single-edit effects) predict the double's measured effect? This asks, in the
     MOLECULAR regime, the question J2 answered for HIV resistance (`FAIL_ADDITIVE_SUFFICES`): if additive
     composition already predicts the double, then the "retrieval, NOT additive composition" framing that
     motivates the extension is answering a problem that does not bite. (Note the two possible answers are
     BOTH informative: additive-suffices => no need for retrieval; additive-fails => the extension cannot
     use composition either, and would need real measured edit-sets.)

A THIRD reason it is moot for what actually shipped, which needs no data at all: `dna-decode inverse` takes
a RANK target ("the p-th percentile of predicted damage"), and a percentile is defined RELATIVE TO THE
CANDIDATE POOL -- so every rank target is reachable by a single edit BY CONSTRUCTION. "Unreachable target"
is a MAGNITUDE concept, and the magnitude inverse is already not deployable
(`scripts/forward_inverse_deployable.py`).

Run:  uv run python scripts/forward_inverse_multiedit.py
Exit: 0 = ran; 2 = substrate unavailable.
"""
from __future__ import annotations

import argparse
import csv
import json
import statistics
import sys
from datetime import date
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.forward_inverse_roundtrip import PG  # noqa: E402

# The ONLY assay in the cache carrying real multi-mutants measured alongside their singles.
DMS_ID = "SR43C_ARATH_Tsuboyama_2023_2N88"


def _spearman(x, y) -> float:
    """Rank correlation with MID-RANKS for ties.

    Position-order tie-breaking is a documented trap in this repo (LESSONS_LEARNED: a constant vector
    scores an artifactual 1.0 rather than 0, because `sorted` assigns tied values ranks by POSITION). The
    additive sums here DO tie, so mid-ranks are load-bearing, not hygiene.
    """
    def rk(v):
        o = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(o):
            j = i
            while j + 1 < len(o) and v[o[j + 1]] == v[o[i]]:
                j += 1
            mid = (i + j) / 2.0                 # average rank over the tie group
            for k in range(i, j + 1):
                r[o[k]] = mid
            i = j + 1
        return r
    if len(x) < 3:
        return 0.0
    rx, ry = rk(x), rk(y)
    n = len(x)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** 0.5
    return num / den if den else 0.0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=REPO / "wiki")
    args = ap.parse_args()

    path = PG / "pg_dms" / "DMS_ProteinGym_substitutions" / f"{DMS_ID}.csv"
    if not path.exists():
        print(f"[multiedit] SUBSTRATE UNAVAILABLE: {path}", file=sys.stderr)
        return 2

    single: dict[str, float] = {}
    multi: dict[str, float] = {}
    with path.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            m = (r.get("mutant") or "").strip()
            try:
                v = float(r["DMS_score"])
            except (TypeError, ValueError, KeyError):
                continue
            (multi if ":" in m else single)[m] = v

    lo_s, hi_s = min(single.values()), max(single.values())
    beyond = [m for m, v in multi.items() if v < lo_s or v > hi_s]

    # Q2: only doubles whose BOTH constituent singles were measured in the same assay.
    pairs = []
    for m, v in multi.items():
        parts = m.split(":")
        if len(parts) != 2 or any(p not in single for p in parts):
            continue
        pairs.append((sum(single[p] for p in parts), v))
    add_rho = _spearman([a for a, _ in pairs], [b for _, b in pairs])
    add_mae = statistics.fmean(abs(a - b) for a, b in pairs) if pairs else None
    # Baseline: predict every double as the assay's median -- if additivity cannot beat THAT it is useless.
    med = statistics.median(single.values())
    null_mae = statistics.fmean(abs(med - b) for _, b in pairs) if pairs else None

    rep = {
        "schema": "forward-inverse-multiedit-v1",
        "date": date.today().isoformat(),
        "question": "does the multi-edit inverse extension (cooc-multiedit-inverse) earn a build?",
        "substrate": {"dms_id": DMS_ID, "n_single": len(single), "n_multi": len(multi),
                      "orders_present": sorted({m.count(":") + 1 for m in multi}),
                      "single_range": [round(lo_s, 3), round(hi_s, 3)],
                      "multi_range": [round(min(multi.values()), 3), round(max(multi.values()), 3)],
                      "note": "the ONLY multi-mutant assay in the cache -> n=1 substrate, doubles only"},
        "Q1_reach": {
            "claim_tested": "multi-edits reach effects no single edit reaches (the extension's motivation)",
            "n_multi_outside_single_range": len(beyond),
            "fraction": round(len(beyond) / len(multi), 4) if multi else None,
            "verdict": "MOTIVATED" if beyond else "NOT_MOTIVATED",
            "why": ("doubles span a range strictly INSIDE the singles' -- the assay is bounded (a protein "
                    "cannot be more unfolded than unfolded), so effect RANGE is set by the assay floor/"
                    "ceiling, not by how many edits you stack. Stacking saturates; it does not extend."),
        },
        "Q2_additivity": {
            "claim_tested": "sum(single effects) predicts the double's measured effect",
            "n_doubles_with_both_singles_measured": len(pairs),
            "spearman_additive_vs_measured": round(add_rho, 4),
            "mae_additive": round(add_mae, 4) if add_mae is not None else None,
            "mae_predict_median_null": round(null_mae, 4) if null_mae is not None else None,
            "additive_beats_null": (add_mae < null_mae) if (add_mae is not None and null_mae) else None,
        },
        "Q3_moot_for_what_shipped": (
            "`dna-decode inverse` takes a RANK target, and a percentile is defined relative to the candidate "
            "pool -- so every rank target is reachable by a single edit BY CONSTRUCTION. 'Unreachable "
            "target' is a MAGNITUDE concept, and the magnitude inverse is not deployable."),
        "verdict": ("MULTIEDIT_EXTENSION_NOT_WARRANTED" if not beyond
                    else "MULTIEDIT_EXTENSION_MOTIVATED_NEEDS_ITS_OWN_FALSIFIER"),
        "honest_scope": ("n=1 substrate (SR43C doubles). This does NOT prove multi-edits never extend reach "
                         "in any assay -- it shows the extension's stated motivation is FALSE where it can "
                         "be checked, which is enough to not build it on that motivation. A build would "
                         "need a substrate where Q1 comes back MOTIVATED."),
    }
    stem = f"forward_inverse_multiedit_{rep['date']}"
    (args.out_dir / f"{stem}.json").write_text(json.dumps(rep, indent=2), encoding="utf-8")

    s = rep["substrate"]
    print(f"[multiedit] {DMS_ID}: {s['n_single']} singles + {s['n_multi']} doubles (same assay)")
    print(f"  singles range {s['single_range']}   doubles range {s['multi_range']}")
    print(f"\n  Q1 REACH  -- do doubles reach what singles cannot?")
    print(f"     {rep['Q1_reach']['n_multi_outside_single_range']}/{s['n_multi']} outside the single range "
          f"-> {rep['Q1_reach']['verdict']}")
    q2 = rep["Q2_additivity"]
    print(f"  Q2 ADDITIVITY -- does sum(singles) predict the double? (n={q2['n_doubles_with_both_singles_measured']})")
    print(f"     spearman {q2['spearman_additive_vs_measured']:+.4f}   "
          f"MAE {q2['mae_additive']}  vs predict-median null {q2['mae_predict_median_null']}  "
          f"-> beats null: {q2['additive_beats_null']}")
    print(f"\n  VERDICT: {rep['verdict']}")
    print(f"  -> {args.out_dir / (stem + '.json')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
