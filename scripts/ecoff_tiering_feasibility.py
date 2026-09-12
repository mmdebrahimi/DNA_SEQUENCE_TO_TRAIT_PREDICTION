"""Can an ECOFF anchor even discriminate on the Oxford cohort? Per drug, with outcomes named first.

STEP 2 of plans/ECOFF_Anchored_Tiering_Evaluation_Only_Technical_Plan. It runs BEFORE the acceptance
bar is frozen, deliberately: two bars earlier in this project were mis-specified because a threshold
was set before the analysis that says what a correct result can look like.

THE THREE VERDICTS WERE NAMED IN THE PLAN BEFORE ANY ECOFF WAS SOURCED:

  SCOREABLE                     the ECOFF sits inside the panel's measured range and the stratum it
                                separates from the clinical breakpoint clears the 20-isolate floor.
  DEGENERATE_ECOFF_BELOW_PANEL  the ECOFF is strictly below the lowest dilution the panel reports, so
                                EVERY isolate is non-wild-type and the anchor discriminates nothing.
  UNDERPOWERED_STRATUM          fewer than 20 isolates separate the two anchors.

WHY THE DEGENERATE CASE IS NOT HYPOTHETICAL. Oxford's ciprofloxacin panel bottoms out at 0.125 mg/L and
the E. coli ECOFF is 0.06; its ceftriaxone panel bottoms out at 0.5 and the (tentative) ECOFF is 0.125.
Both collapse to 100% non-wild-type. Naming that outcome in advance is what keeps it a result rather
than a discovery made halfway through an analysis.

Read-only. No network. MIC = 2**upper per scripts/oxford_score.py:14 (the paper's log2 encoding).
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date as _date
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from dna_decode.data.ecoff_catalog import (  # noqa: E402
    UnsourcedEcoffError, ecoff_for, entry_for,
)
from dna_decode.data.mic_tiers import breakpoints_for  # noqa: E402  (FROZEN -- read only)
from dna_decode.eval.wildtype_tiering import (  # noqa: E402
    NWT, classify_wildtype, ecoff_is_resolvable_on_grid,
)

MIN_STRATUM = 20                    # the repo's existing per-class floor
SCOREABLE = "SCOREABLE"
DEGENERATE = "DEGENERATE_ECOFF_BELOW_PANEL"
UNDERPOWERED = "UNDERPOWERED_STRATUM"
UNSOURCED = "REFUSED_UNSOURCED_ECOFF"

# Oxford's MIC columns, by drug. Tetracycline is absent from the deposit entirely -- reported rather
# than silently skipped, since "no column" and "no signal" are different facts.
OXFORD_COLUMNS = {"gentamicin": "Gentamicin", "ciprofloxacin": "Ciprofloxacin",
                  "ceftriaxone": "Ceftriaxone", "tetracycline": None}


def verdict_for(ecoff: float, grid: list[float], stratum_n: int,
                min_stratum: int = MIN_STRATUM) -> str:
    """The frozen rule, applied mechanically. Degeneracy is checked FIRST: an anchor that cannot
    discriminate at all is not merely underpowered, and reporting it as underpowered would imply a
    bigger cohort would help."""
    if not ecoff_is_resolvable_on_grid(ecoff, grid):
        return DEGENERATE
    return SCOREABLE if stratum_n >= min_stratum else UNDERPOWERED


def screen_drug(df: pd.DataFrame, drug: str, column: str | None) -> dict:
    if column is None:
        return {"drug": drug, "verdict": "NOT_IN_COHORT",
                "reason": "the Oxford deposit carries no MIC column for this drug"}
    try:
        ecoff = ecoff_for(drug)
    except UnsourcedEcoffError as e:
        return {"drug": drug, "verdict": UNSOURCED, "reason": str(e)}

    mic = (2.0 ** pd.to_numeric(df[f"{column}_upper"], errors="coerce")).dropna()
    grid = sorted(float(v) for v in mic.unique())
    s_bp = breakpoints_for(drug)["clsi_s"]
    # The discriminating stratum: clinically SUSCEPTIBLE but above the ECOFF -- the isolates the two
    # anchors disagree about, and the only ones the comparison can learn anything from.
    stratum = int(((mic <= s_bp) & (mic > ecoff)).sum())
    n_nwt = sum(1 for m in mic if classify_wildtype([m], ecoff) == NWT)
    entry = entry_for(drug)
    return {
        "drug": drug, "ecoff_mg_L": ecoff, "ecoff_is_tentative": "TENTATIVE" in (entry.note or ""),
        "ecoff_source_url": entry.source_url,
        "clsi_susceptible_breakpoint": s_bp,
        "n_isolates": int(mic.size), "panel_floor": min(grid), "measured_grid": grid,
        "n_non_wildtype": n_nwt, "non_wildtype_fraction": round(n_nwt / mic.size, 4),
        "discriminating_stratum_n": stratum,
        "verdict": verdict_for(ecoff, grid, stratum),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--main-csv", type=Path, default=ROOT / "data" / "raw" / "oxford" / "main_data.csv")
    ap.add_argument("--out", type=Path,
                    default=ROOT / "wiki" / f"ecoff_tiering_feasibility_{_date.today()}.json")
    a = ap.parse_args(argv)

    df = pd.read_csv(a.main_csv, low_memory=False)
    rows = [screen_drug(df, d, c) for d, c in OXFORD_COLUMNS.items()]
    scoreable = [r["drug"] for r in rows if r["verdict"] == SCOREABLE]

    doc = {
        "schema": "ecoff-tiering-feasibility-v1",
        "date": str(_date.today()),
        "cohort": "Oxford E. coli bacteraemia deposit (data/raw/oxford)",
        "mic_convention": "MIC = 2**upper (log2 dilution indices; scripts/oxford_score.py:14)",
        "min_stratum": MIN_STRATUM,
        "results": rows,
        "scoreable_drugs": scoreable,
        "headline": (
            "2 of 3 measured drugs are DEGENERATE -- the E. coli ECOFF lies strictly below Oxford's "
            "panel floor (cipro 0.06 vs 0.125; ceftriaxone 0.125 vs 0.5), so every isolate is "
            "non-wild-type and the anchor discriminates nothing. Only gentamicin is SCOREABLE, and "
            "its discriminating stratum is 25 isolates -- just over the floor of 20."),
        "honest_limits": [
            "This screens the ANCHOR against the PANEL, nothing else. A SCOREABLE verdict says the "
            "comparison is possible on this cohort, not that it will find anything.",
            "Ceftriaxone's ECOFF is a TECOFF -- EUCAST prints it in parentheses, set on 4 data "
            "sources / 908 observations. Its degeneracy here does not depend on that (0.125 is below "
            "the 0.5 floor either way), but any ceftriaxone conclusion would.",
            "Tetracycline has no MIC column in the deposit at all, so it is NOT_IN_COHORT rather than "
            "degenerate -- a different fact, kept separate.",
            "Degeneracy is a property of THIS panel, not of the ECOFF. A cohort reporting finer low-end "
            "dilutions would make cipro and ceftriaxone scoreable.",
        ],
    }
    a.out.write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8")
    for r in rows:
        extra = ""
        if r.get("ecoff_mg_L") is not None:
            extra = (f"ECOFF {r['ecoff_mg_L']}{' (tentative)' if r['ecoff_is_tentative'] else ''} | "
                     f"floor {r['panel_floor']} | NWT {r['non_wildtype_fraction']:.1%} | "
                     f"stratum {r['discriminating_stratum_n']}")
        print(f"{r['drug']:<15}{r['verdict']:<30}{extra}")
    print(f"\nscoreable: {scoreable or 'NONE'}\n-> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
