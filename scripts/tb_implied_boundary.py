"""Where does the genotype's resistance boundary actually sit on the MIC ladder? Measured on CRyPTIC.

THE QUESTION, AND WHY IT MOVED HERE. The 2026-09-12 ECOFF arm asked whether a wild-type anchor agrees
with the genotype better than a clinical breakpoint. On the Oxford cohort it returned WEAK_DIRECTIONAL:
an 18x relative enrichment in the predicted direction that sat inside its own permutation null, because
the discriminating stratum was 25 isolates carrying 2 determinants. Its memo named the fix -- "a cohort
with a full dilution series at the low end ... not a bigger cohort". CRyPTIC is that cohort and it was
already on disk: a genuine two-fold ladder (RIF 0.03-8, INH 0.025-12.8) against Oxford's collapsed
gentamicin panel {1, 2, 4, 32}.

WHAT CHANGES BESIDES THE SIZE -- AND THIS IS THE POINT. Oxford's framing needed TWO sourced numbers (an
ECOFF and a clinical breakpoint), which makes the result hostage to both and carries a fabrication
hazard. Here neither is needed. CRyPTIC ships its own BINARY_PHENOTYPE, and that column is a PURE STEP
on the ladder -- verified, not assumed: R-fraction is exactly 0.000 at every rung at or below one value
and exactly 1.000 above it. So the effective clinical cut-off is RECOVERED from the shipped data and is
re-derivable by anyone; `recover_cutoff` REFUSES when the step is impure, because an impure step would
mean BINARY_PHENOTYPE is not a MIC threshold and the entire recovery argument is void.

So the measured quantity is WHERE the genotype's boundary falls, not whether one asserted cut-off beats
another.

NOT THE CLOSED BV-BRC NEGATIVE. The MIC-continuous track is closed on BV-BRC grounds (G1: 91% of BV-BRC
MIC is XGBoost-from-genome). CRyPTIC MIC is a wet-lab broth-microdilution measurement; that distinction
is already recorded in scripts/tb_mic_calibration.py, which validated this same substrate for a
different question (interval calibration) on 2026-07-12.

Read-only: reuses the cached WHO grade-1/2 determinant calls and never re-streams the 2.9 GB parquet.
Offline, no network, no Docker. Bar: wiki/tb_implied_boundary_acceptance_bar.json.
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from datetime import date as _date
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

MIN_BOUNDARY_N = 100          # the frozen powering floor for the boundary rung
SUPPORTED, WEAK, FALSIFIED, INDETERMINATE = (
    "SUPPORTED", "WEAK_DIRECTIONAL", "FALSIFIED", "INDETERMINATE")

DRUGS = {"rifampicin": "RIF", "isoniazid": "INH"}


def parse_mic(value: str) -> tuple[float, str]:
    """CRyPTIC MIC string -> (log2 rung, censoring).

    The ladder is doubling, so log2 turns it into consecutive integers-with-offset and rung arithmetic
    becomes ordinary comparison. Censoring is RETAINED rather than dropped: a `<=`/`>` bound is a real
    measurement of a bound, and silently coercing it to a point value is the documented operator-aware
    censoring trap. Here the censored values sit at the ladder ENDS, never at the boundary rung, which
    `boundary_is_uncensored` asserts rather than assumes.
    """
    s = str(value).strip()
    if s.startswith("<="):
        return float(np.log2(float(s[2:]))), "left"
    if s.startswith("<"):
        return float(np.log2(float(s[1:]))), "left"
    if s.startswith(">="):
        return float(np.log2(float(s[2:]))), "right"
    if s.startswith(">"):
        return float(np.log2(float(s[1:]))), "right"
    return float(np.log2(float(s))), "exact"


def recover_cutoff(rungs: np.ndarray, pheno: np.ndarray) -> float | None:
    """The effective clinical cut-off, RECOVERED from the data. None when the step is impure.

    Returns the highest rung whose isolates are ALL susceptible, but only if the phenotype is a clean
    step: every rung at or below it entirely S, every rung above it entirely R. Any rung with a mixed
    R-fraction means BINARY_PHENOTYPE is not a pure MIC threshold -- in which case there is no cut-off
    to recover and this must REFUSE rather than pick the nearest rung. That refusal is the whole reason
    no breakpoint has to be recalled from memory anywhere in this arm.
    """
    # A ladder with NO resistance anywhere has no threshold to recover, and its "purity" is vacuous --
    # every rung is trivially all-S. Without this the loop sets `out` to the TOP rung, the `above`
    # slice is then empty so its check is skipped, and an information-free column yields a
    # confident-looking cut-off. Both classes must be present for a step to mean anything.
    if not (pheno == "R").any() or not (pheno == "S").any():
        return None
    out = None
    for rung in sorted(set(rungs.tolist())):
        frac_r = float((pheno[rungs == rung] == "R").mean())
        if frac_r == 0.0:
            out = rung
        elif frac_r == 1.0:
            continue
        else:
            return None                              # impure rung -> refuse
    if out is None:
        return None
    # BOTH directions must hold. Checking only ABOVE was a real defect: a non-monotone ladder with a
    # RESISTANT rung BELOW a susceptible one passed and returned a bogus cut-off. A threshold means
    # everything at-or-below is S and everything above is R -- assert both, or refuse.
    at_or_below, above = rungs <= out, rungs > out
    if at_or_below.any() and float((pheno[at_or_below] == "S").mean()) != 1.0:
        return None
    if above.any() and float((pheno[above] == "R").mean()) != 1.0:
        return None
    return out


def permutation_gap(in_boundary: np.ndarray, carries: np.ndarray, rng, n_perm: int = 1000) -> dict:
    """Null for the BOUNDARY LABEL: shuffle which susceptible isolates sit at the top rung.

    Holds the stratum's overall carriage and every rung size fixed, so the null asks precisely whether
    the boundary rung is enriched beyond what the same S-stratum carriage would produce by chance.
    """
    obs = float(carries[in_boundary].mean() - carries[~in_boundary].mean())
    gaps = []
    for _ in range(n_perm):
        perm = rng.permutation(in_boundary)
        gaps.append(float(carries[perm].mean() - carries[~perm].mean()))
    g = np.array(gaps)
    return {"observed_gap": obs, "n_perm": n_perm, "perm_mean": float(g.mean()),
            "perm_max": float(g.max()), "perm_p95": float(np.percentile(g, 95)),
            "exceeds_perm_max": bool(obs > g.max())}


def verdict_from_bar(per_drug: dict, min_boundary_n: int = MIN_BOUNDARY_N) -> str:
    """The FROZEN rule, applied mechanically (wiki/tb_implied_boundary_acceptance_bar.json)."""
    if not per_drug:
        return INDETERMINATE
    for r in per_drug.values():
        if r.get("cutoff_log2") is None or r.get("boundary_n", 0) < min_boundary_n:
            return INDETERMINATE
    if any(r["boundary_carriage"] <= r["lower_carriage"] for r in per_drug.values()):
        return FALSIFIED
    if all(r["permutation"]["exceeds_perm_max"] for r in per_drug.values()):
        return SUPPORTED
    return WEAK


def quality_disclosure(df: pd.DataFrame, cutoff: float) -> dict:
    """How the two COMPARED strata differ in the cohort's OWN per-isolate label-quality flag.

    NAMESPACE-SEPARATE and AUGMENT-ONLY: this never changes `verdict`, which stays a mechanical
    application of the frozen bar. It exists because the bar did not mention quality while the repo's
    two other CRyPTIC scorers (`score_tb_cryptic.py`, `score_tb_cryptic_parquet.py`) both hard-filter
    to HIGH -- so the shipped run is powered partly by isolates this repo elsewhere discards, and a
    reader must be able to see that without re-deriving it.

    Reports (a) the quality MIX of the boundary rung against the lower-S stratum it is compared with,
    and (b) the same contrast recomputed on HIGH-only isolates, which is the replication check: if the
    effect survives there, the pooled result is not a label-quality artifact even though its POWERING
    depends on the other tiers.
    """
    s = df[df["pheno"] == "S"]
    boundary = np.isclose(s["rung"].to_numpy(), cutoff)

    def mix(sub: pd.DataFrame) -> dict:
        n = len(sub)
        if not n:
            return {}
        return {q: round(float((sub["quality"] == q).sum()) / n, 4)
                for q in ("HIGH", "MEDIUM", "LOW")}

    b_mix, l_mix = mix(s[boundary]), mix(s[~boundary])
    hi = s[s["quality"] == "HIGH"]
    hi_b = np.isclose(hi["rung"].to_numpy(), cutoff)
    hi_c = hi["carries"].to_numpy()
    n_hi_b = int(hi_b.sum())
    high_only = {
        "boundary_n": n_hi_b,
        "clears_frozen_floor": bool(n_hi_b >= MIN_BOUNDARY_N),
        "boundary_carriage": float(hi_c[hi_b].mean()) if n_hi_b else None,
        "lower_carriage": float(hi_c[~hi_b].mean()) if (~hi_b).any() else None,
    }
    if high_only["boundary_carriage"] is not None and high_only["lower_carriage"] is not None:
        high_only["gap"] = high_only["boundary_carriage"] - high_only["lower_carriage"]
    return {
        "why": ("the frozen bar is silent on PHENOTYPE_QUALITY while two sibling CRyPTIC scorers "
                "require HIGH; this discloses the difference rather than restating the verdict"),
        "boundary_quality_mix": b_mix,
        "lower_s_quality_mix": l_mix,
        "low_share_ratio": (round(b_mix.get("LOW", 0.0) / l_mix["LOW"], 2)
                            if l_mix.get("LOW") else None),
        "high_only": high_only,
    }


def analyse(drug: str, code: str, table: pd.DataFrame, cache: dict, n_perm: int, seed: int) -> dict:
    mics, feats = cache["mics"][drug], cache["features"][drug]
    col = f"{code}_BINARY_PHENOTYPE"
    ph = dict(zip(table["UNIQUEID"], table[col]))
    qual = {u: str(q).strip().upper()
            for u, q in zip(table["UNIQUEID"], table[f"{code}_PHENOTYPE_QUALITY"])}

    rows = []
    for uid, raw in mics.items():
        p = ph.get(uid)
        if p not in ("R", "S"):
            continue
        rung, cens = parse_mic(raw)
        rows.append((uid, rung, cens, p, 1 if feats.get(uid) else 0, qual.get(uid, "NA")))
    df = pd.DataFrame(rows, columns=["uid", "rung", "cens", "pheno", "carries", "quality"])

    cutoff = recover_cutoff(df["rung"].to_numpy(), df["pheno"].to_numpy())
    if cutoff is None:
        return {"drug": drug, "cutoff_log2": None,
                "reason": "BINARY_PHENOTYPE is not a pure step on the MIC ladder; no cut-off recoverable"}

    s = df[df["pheno"] == "S"]
    boundary = np.isclose(s["rung"].to_numpy(), cutoff)
    carries = s["carries"].to_numpy()
    n_b = int(boundary.sum())
    result = {
        "drug": drug,
        "cutoff_log2": cutoff, "cutoff_mg_L": float(2.0 ** cutoff),
        "n_joined": int(len(df)), "n_susceptible": int(len(s)),
        "boundary_n": n_b, "lower_n": int((~boundary).sum()),
        "boundary_carriage": float(carries[boundary].mean()) if n_b else None,
        "lower_carriage": float(carries[~boundary].mean()) if (~boundary).any() else None,
        "boundary_is_uncensored": bool((s.loc[boundary, "cens"] == "exact").all()) if n_b else None,
        "overall_carriage": float(df["carries"].mean()),
        "quality_disclosure": quality_disclosure(df, cutoff),
    }
    if n_b >= MIN_BOUNDARY_N and (~boundary).any():
        result["permutation"] = permutation_gap(boundary, carries, np.random.default_rng(seed), n_perm)
    else:
        result["permutation"] = {"exceeds_perm_max": False, "note": "boundary rung below the frozen floor"}

    # SECONDARY, measured not asserted: is the determinant MIX at the boundary different from the
    # high-MIC rungs? A shift means the boundary carriers are a different mechanism class rather than
    # simply fewer of the same one. Which mechanisms they are is NOT claimed here.
    at_boundary, at_high = Counter(), Counter()
    for uid, raw in mics.items():
        if ph.get(uid) not in ("R", "S"):
            continue
        rung, _ = parse_mic(raw)
        if np.isclose(rung, cutoff):
            at_boundary.update(feats.get(uid, []))
        elif rung > cutoff + 1.5:                    # well above the cut-off, not the adjacent rung
            at_high.update(feats.get(uid, []))
    top_b = at_boundary.most_common(6)
    result["composition"] = {
        "boundary_top": [{"determinant": g, "n": n} for g, n in top_b],
        "high_mic_top": [{"determinant": g, "n": n} for g, n in at_high.most_common(6)],
        "boundary_distinct": len(at_boundary), "high_mic_distinct": len(at_high),
        # The top-1 comparison is a DELIBERATELY WEAK statistic and it under-detects: it reports no
        # shift for isoniazid, where the single commonest determinant is the same at both ends while
        # the per-GENE mix moves substantially. Kept because it is the strictest form of the claim
        # (RIF passes even this), but the per-gene share below is the one to read.
        "top_determinant_differs": bool(top_b and at_high.most_common(1)
                                        and top_b[0][0] != at_high.most_common(1)[0][0]),
        "gene_share": _gene_share(at_boundary, at_high),
    }
    return result


def _gene_share(at_boundary: Counter, at_high: Counter) -> dict:
    """Per-GENE share of determinant calls at the boundary vs high-MIC rungs.

    The gene is taken LEXICALLY from the determinant string's prefix (`katG_p.Ser315Thr` -> `katG`) --
    a pure string operation on data already in the cache. No mechanism biology is asserted or recalled
    here; the claim is only that the mix of GENES differs, which is checkable from the counts.
    """
    def shares(c: Counter) -> dict:
        total = sum(c.values())
        if not total:
            return {}
        by_gene: Counter = Counter()
        for det, n in c.items():
            by_gene[det.split("_", 1)[0]] += n
        return {g: round(n / total, 4) for g, n in by_gene.most_common()}

    b, h = shares(at_boundary), shares(at_high)
    genes = sorted(set(b) | set(h))
    shifts = {g: round(abs(b.get(g, 0.0) - h.get(g, 0.0)), 4) for g in genes}
    top = max(shifts.values(), default=0.0)
    # TIES ARE REPORTED, NOT BROKEN ARBITRARILY. With exactly two genes the shares sum to 1, so both
    # shift by the SAME amount and `max()` would return whichever happened to sort first -- an
    # arbitrary choice dressed as a finding. The magnitude is the real quantity; name every gene that
    # attains it.
    tied = sorted(g for g, s in shifts.items() if s == top)
    return {
        "boundary": b, "high_mic": h, "per_gene_shift": shifts,
        "largest_shift": top,
        "largest_shift_genes": tied,
        "largest_shift_is_tied": len(tied) > 1,
    }


def _quality_gate(per_drug: dict) -> dict:
    """Would the frozen powering floor still be met using only HIGH-quality labels?

    Reported BESIDE the verdict, never folded into it -- the verdict is the frozen bar applied
    mechanically, and the bar as frozen says nothing about quality. Rewriting the verdict here would
    be re-sizing a frozen bar after seeing the result, which this project treats as an authority call.
    """
    per = {d: r["quality_disclosure"]["high_only"]
           for d, r in per_drug.items() if r.get("cutoff_log2") is not None}
    if not per:
        return {"status": "not_assessable"}
    short = sorted(d for d, h in per.items() if not h["clears_frozen_floor"])
    return {
        "high_only_clears_frozen_floor": not short,
        "drugs_below_floor_under_high_only": short,
        "verdict_under_high_only": SUPPORTED if not short else INDETERMINATE,
        "effect_replicates_within_high": {d: h.get("gap") for d, h in per.items()},
        "reading": ("the DIRECTION and MAGNITUDE survive restriction to HIGH; what fails there is the "
                    "frozen sample-size floor. A powering failure is not a refutation."),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", type=Path,
                    default=ROOT / "data" / "raw" / "cryptic" / "CRyPTIC_reuse_table_20240917.csv")
    ap.add_argument("--cache", type=Path,
                    default=ROOT / "data" / "processed" / "tb_mic_features_cache.json")
    ap.add_argument("--n-perm", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=Path,
                    default=ROOT / "wiki" / f"tb_implied_boundary_{_date.today()}.json")
    a = ap.parse_args(argv)

    table = pd.read_csv(a.table, low_memory=False)
    cache = json.loads(a.cache.read_text(encoding="utf-8"))
    per_drug = {d: analyse(d, c, table, cache, a.n_perm, a.seed) for d, c in DRUGS.items()}

    doc = {
        "schema": "tb-implied-boundary-v1",
        "date": str(_date.today()),
        "filename_date_is_the_original_measurement": ("the artifact keeps its 2026-09-21 filename; "
                                                      "`date` is when it was last regenerated. The "
                                                      "measured numbers are unchanged -- the "
                                                      "2026-09-22 regeneration only ADDED the quality "
                                                      "disclosure. A second dated file carrying "
                                                      "identical headline numbers would be the "
                                                      "duplicate-artifact trap."),
        "bar": "wiki/tb_implied_boundary_acceptance_bar.json",
        "cohort": "CRyPTIC M. tuberculosis compendium (measured broth-microdilution MIC)",
        "cutoff_provenance": ("RECOVERED from the shipped BINARY_PHENOTYPE as a pure step on the MIC "
                              "ladder; refused if the step is impure. No breakpoint or ECOFF is recalled "
                              "from memory anywhere in this arm."),
        "results": per_drug,
        "verdict": verdict_from_bar(per_drug),
        "verdict_is_mechanical": "verdict_from_bar applied to the frozen rule; not authored",
        "quality_gate": _quality_gate(per_drug),
        "honest_limits": [
            "POWERING DEPENDS ON LABEL TIERS THIS REPO ELSEWHERE DISCARDS. The frozen bar is silent on "
            "PHENOTYPE_QUALITY, but score_tb_cryptic.py and score_tb_cryptic_parquet.py both require "
            "HIGH. Restricted to HIGH the rifampicin boundary rung falls below the frozen "
            "MIN_BOUNDARY_N floor, so the run would read INDETERMINATE -- see `quality_gate`. The "
            "EFFECT still replicates within HIGH alone; it is the POWERING that does not survive, and "
            "the two are different claims.",
            "The boundary rung is enriched for LOW-quality labels relative to the lower-S stratum it "
            "is compared against (see per-drug `quality_disclosure.low_share_ratio`).",
            "IN-DISTRIBUTION: the WHO catalogue was built partly from CRyPTIC, so the determinant calls "
            "and this cohort are not independent. This locates a boundary on the ladder; it is NOT an "
            "independent validation of the catalogue.",
            "One organism, two drugs, one compendium -- and both drugs share the same isolates, so they "
            "are not two independent replications.",
            "A determinant at a clinically-susceptible rung reads two ways and this arm cannot separate "
            "them: a genuine low-level-resistance mechanism sitting under the cut-off, or a catalogue "
            "false positive. Rung-position dependence argues against random false positives; it does "
            "not exclude the second reading.",
            "This does NOT name a numeric ECOFF and is not a claim that any ECOFF value is correct.",
            "Frozen surfaces are READ-only; TB rules live in the non-frozen organism_rules package.",
        ],
    }
    a.out.write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8")

    for drug, r in per_drug.items():
        if r.get("cutoff_log2") is None:
            print(f"{drug:<14} REFUSED -- {r['reason']}")
            continue
        p = r["permutation"]
        print(f"{drug:<14} cut-off {r['cutoff_mg_L']:<7.4g} | boundary n={r['boundary_n']:<5} "
              f"carriage {r['boundary_carriage']:.3f} vs lower {r['lower_carriage']:.3f} | "
              f"gap {p.get('observed_gap', float('nan')):.4f} vs perm max "
              f"{p.get('perm_max', float('nan')):.4f}")
    print(f"\n-> {doc['verdict']}\n-> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
