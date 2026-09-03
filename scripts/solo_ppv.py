"""SOLO-PPV: grade a determinant on the occurrences where it appears ALONE (the WHO TB catalogue method).

WHY THIS, AND WHY IT MATTERS HERE. A pooled PPV over every carrier of a determinant silently credits the
determinant with resistance that a CO-OCCURRING determinant may have caused. The WHO *M. tuberculosis*
mutation catalogue solves this by grading only "solo" occurrences -- isolates where the variant appears
without another known resistance determinant -- then scoring PPV with a confidence interval. It is the
published, adopted answer to the confound, and it is a strictly stronger claim than a pooled rate.

The project met this confound twice from the other side and did not name it:

  - The NNRTI curation (2026-09-01) used a MULTIVARIATE OLS to deconfound co-occurrence, and its
    best-scoring variant DELETED canonical Y181C -- because Y181C co-occurs heavily with K103N, which
    absorbed the coefficient. SOLO handles the same confound by EXCLUSION rather than by regression, so
    it cannot silently reassign a determinant's effect to its co-traveller.
  - The gentamicin `rmt` hunt (2026-09-02) reported a POOLED PPV of 146/206. That pooled number credits
    `rmt` for 47 carriers that ALSO carry `aac(3)`, the classic gentamicin-modifying enzyme.

WHAT IT IS NOT. Solo-grading is not universally better: a Nature Communications 2025 analysis found
multivariable penalised regression improves sensitivity over the WHO SOLO catalogue and recovers
compensatory and hypersusceptibility variants that solo-only grading misses. Solo is the CONSERVATIVE
estimator -- it discards data to buy an uncontaminated denominator. Report both when both are available.

Pure + offline. Reads committed artifacts; writes wiki/solo_ppv_<date>.json.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# WHO catalogue grade-1 ("associated with resistance") bar, applied here to a determinant rather than a
# TB variant: at least 5 solo occurrences AND a PPV whose 95% CI lower bound clears 0.25.
WHO_MIN_SOLO = 5
WHO_MIN_PPV_CI_LOWER = 0.25


def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """Wilson score interval -- correct at the extremes, where a Wald interval degenerates.

    This matters directly: a perfect 99/99 gets a Wald interval of [1.0, 1.0], which asserts certainty
    the data cannot support. Wilson gives [0.963, 1.0].
    """
    if n <= 0:
        return (0.0, 1.0)
    p = k / n
    den = 1.0 + z * z / n
    centre = p + z * z / (2 * n)
    margin = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (max(0.0, (centre - margin) / den), min(1.0, (centre + margin) / den))


def score(r: int, s: int) -> dict:
    n = r + s
    lo, hi = wilson(r, n)
    return {"n": n, "resistant": r, "susceptible": s,
            "ppv": (r / n) if n else None, "ci95_lower": lo, "ci95_upper": hi,
            "meets_who_grade1_bar": bool(n >= WHO_MIN_SOLO and lo >= WHO_MIN_PPV_CI_LOWER)}


def gentamicin_rmt(exclude_project: str | None) -> dict:
    """SOLO-grade the deployed gentamicin `rmt` rescue against its own committed evidence."""
    art = ROOT / "wiki" / "gentamicin_rmt_specificity_hunt.json"
    if not art.is_file():
        return {"error": f"missing {art}"}
    d = json.loads(art.read_text(encoding="utf-8"))
    R, S = d["rmt_R_records"], d["rmt_S_records"]

    def split(recs, solo_only, drop_project):
        out = []
        for x in recs:
            if solo_only and x.get("has_gent_aac3"):
                continue                     # co-carries the classic gentamicin enzyme -> not solo
            if drop_project and x.get("bioproject_acc") == drop_project:
                continue
            out.append(x)
        return out

    strata = {
        "pooled_all_carriers": score(len(R), len(S)),
        "solo_no_aac3": score(len(split(R, True, None)), len(split(S, True, None))),
        "co_carriage_with_aac3": score(len([x for x in R if x.get("has_gent_aac3")]),
                                       len([x for x in S if x.get("has_gent_aac3")])),
    }
    if exclude_project:
        strata["solo_no_aac3_artifact_project_excluded"] = score(
            len(split(R, True, exclude_project)), len(split(S, True, exclude_project)))
    return {"determinant": "rmt* / npmA (gentamicin v2 symbol_rescue)",
            "excluded_project": exclude_project, "strata": strata}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--exclude-project", default="PRJNA1322038",
                    help="the LABEL_ARTIFACT submission (see wiki/gentamicin_rmt_project_control.json)")
    ap.add_argument("--out", type=Path, default=ROOT / "wiki" / "solo_ppv_2026-09-03.json")
    a = ap.parse_args()

    res = gentamicin_rmt(a.exclude_project)
    if "error" in res:
        print(res["error"])
        return 2

    print(f"SOLO-PPV grading of {res['determinant']}\n")
    for name, st in res["strata"].items():
        if st["n"] == 0:
            print(f"  {name:44} (empty stratum)")
            continue
        print(f"  {name:44} R={st['resistant']:4d} S={st['susceptible']:3d}  "
              f"PPV={st['ppv']:.4f}  95%CI [{st['ci95_lower']:.3f}, {st['ci95_upper']:.3f}]"
              f"{'  [clears WHO grade-1 bar]' if st['meets_who_grade1_bar'] else ''}")

    out = {"schema": "solo-ppv-v1",
           "method": ("WHO M. tuberculosis mutation catalogue SOLO grading, applied to a determinant "
                      "family: score only occurrences where the determinant appears WITHOUT another "
                      "known determinant for the same drug. Grade-1 bar: >=5 solo occurrences and a "
                      "PPV 95% CI lower bound >= 0.25."),
           "who_bar": {"min_solo": WHO_MIN_SOLO, "min_ppv_ci_lower": WHO_MIN_PPV_CI_LOWER},
           "result": res,
           "honest_limits": [
               "SOLO is the CONSERVATIVE estimator: it discards co-carriage data to buy an "
               "uncontaminated denominator. Multivariable penalised regression has been shown to beat "
               "SOLO grading on sensitivity for TB, recovering compensatory variants solo misses.",
               "'Solo' here means no co-carried aac(3), the classic gentamicin-modifying enzyme -- it "
               "is NOT a full co-determinant screen, because the committed hunt records only that one "
               "co-carriage flag per isolate. A wider screen would lower the solo count further.",
               "The carrier call is NCBI's own AMRFinder (PD AMR_genotypes), a tool-derived feature; "
               "only the phenotype is measured.",
               "Excluding the artifact project is justified by a separate pre-registered control "
               "(wiki/gentamicin_rmt_project_control.json), not by the numbers being inconvenient.",
           ]}
    a.out.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
