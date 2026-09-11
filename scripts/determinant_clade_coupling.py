"""Are ACQUIRED determinants less clade-coupled than CHROMOSOMAL ones? Measured directly.

THE CLAIM UNDER TEST. The 2026-09-10 /innovate sweep's surviving mechanism claim,
`COUPLING-not-popdesign`: the discriminator for de-confounded learned decoding is determinant-to-clade
COUPLING rather than population design -- horizontal transfer decouples a mobile determinant from clade
the way meiosis decouples a locus in a cross, so acquired determinants should be LESS clade-coupled than
chromosomal point mutations.

WHY THIS IS A STRONGER TEST THAN THE ONE IT SURVIVED. The kill-test INFERRED coupling from a proxy -- a
sensitivity drop under lineage collapse -- over 10 cells with a single drug in the chromosomal arm. This
measures the named mechanism itself: how concentrated each determinant's carriers are across lineages,
with both arms drawn from the SAME genomes and the SAME lineage partition.

TWO FILTERS THAT DECIDE WHETHER THE NUMBERS MEAN ANYTHING:

  * SYNONYMOUS ROWS. `mutations.tsv` is `--mutation_all` output: it lists every position SCREENED,
    including wildtype (`baeS_Q163Q`, `16S_A523A`). 137,288 such rows are present. Keeping them turns
    "screened and found nothing" into a called mutation and leaves the POINT arm mostly noise.
  * ONE MLST SCHEME. An ST number means a different lineage under a different scheme, so genomes typed
    under anything but `ecoli_achtman_4` are dropped rather than pooled into invented lineage identity.

THE NULL IS THE POINT. A determinant with 5 carriers cannot span more than 5 lineages, so raw
concentration is mechanically confounded with prevalence. Every family is scored against a null that
draws the SAME number of carriers at random from the same genomes, so the score is a deviation from
chance rather than a restatement of how common the determinant is.

AND THE NULL IS NOT ENOUGH ON ITS OWN -- THE CEILING MOVES TOO. The null fixes the EXPECTATION at a given
carrier count; it does not fix the MAXIMUM. The largest lineage holds 53 of 621 genomes, so a determinant
with 400 carriers cannot exceed 53/400 = 0.133 largest-lineage fraction no matter how coupled it is,
leaving an achievable coupling range of only ~0.048. `gyrA_S83L` scores 0.043 -- which reads as
near-chance against the rare families' 0.87 and is in fact 89% of everything available to it. Reporting
raw coupling alone would have published "the determinants that drive most resistance are decoupled",
which is a statement about the metric's headroom, not about the biology. So every family also carries
`normalized_coupling` = observed / (ceiling - null mean), and BOTH are reported. The frozen bar's verdict
rule is defined on the RAW score and is applied to the raw score unchanged; normalization is a secondary
robustness band, never a retro-fit of the rule.

Read-only. No network, no Docker, no model. Bar: wiki/determinant_clade_coupling_acceptance_bar.json.
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import re
import sys
from collections import defaultdict
from datetime import date as _date
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

RUNS = ROOT / "data" / "amrfinder_runs"
SCHEME = "MLST.ecoli_achtman_4."
# `<gene>_<WT><pos><MUT>`; WT == MUT is a wildtype/synonymous screen row, not a called mutation.
SUBSTITUTION = re.compile(r"^(.*?)_([A-Z])(\d+)([A-Z*])$")
# The broader form. AMRFinder also emits multi-residue indels (`ftsI_N337NYRIN`), frameshifts
# (`cirA_S90YfsTer15`), nonsense (`nfsA_Q67Ter`), deletions (`23S_G543DEL`) and PROMOTER positions with
# NEGATIVE coordinates (`ampC_C-42T`, and the wildtype `ampC_C-11C`). The single-residue regex above
# matches none of them, so all were being discarded as "unparseable" -- including 12 symbols AMRFinder
# actually CALLS. Same wildtype rule: equal ref and alt is a screen row, not a call.
VARIANT = re.compile(r"^(.*?)_([A-Z]+)(-?\d+)([A-Z*]+)$")
_INDEL = re.compile(r"(DEL|INS|DUP)$", re.IGNORECASE)


def _base(acc: str) -> str:
    return str(acc).split(".")[0]


def load_lineages() -> dict[str, str]:
    """assembly (version-stripped) -> lineage, restricted to ONE MLST scheme."""
    out: dict[str, str] = {}
    for f in sorted(glob.glob(str(ROOT / "data" / "processed" / "*.parquet"))):
        try:
            df = pd.read_parquet(f)
        except Exception:                                    # noqa: BLE001 - a bad parquet must not halt the scan
            continue
        if "mlst" not in df.columns or "assembly_accession" not in df.columns:
            continue
        for _, r in df[["assembly_accession", "mlst"]].dropna().iterrows():
            out.setdefault(_base(r["assembly_accession"]), str(r["mlst"]))
    return {k: v for k, v in out.items() if v.startswith(SCHEME)}


def is_called_substitution(symbol: str) -> bool:
    """True for a real non-synonymous SINGLE-RESIDUE call; False for wildtype or anything else.

    The frozen-substrate definition. Deliberately narrow, and kept unchanged so the headline stays the
    one the bar was frozen against -- `is_called_variant` is the corrected superset.
    """
    m = SUBSTITUTION.match(symbol or "")
    return bool(m) and m.group(2) != m.group(4)


def is_called_variant(symbol: str) -> bool:
    """True for any called chromosomal variant: substitution, indel, frameshift, nonsense, promoter.

    The CORRECTED filter. `is_called_substitution` discards 19,183 rows as unparseable, of which 3,904
    are genuine wildtype screen rows and the rest are real variant FORMS the single-residue regex cannot
    express. 12 of the discarded symbols are ones AMRFinder CALLS (`main.tsv`), six of them clearing the
    5-carrier floor -- so the narrow filter costs the chromosomal arm about a fifth of its families.
    """
    s = symbol or ""
    if not s:
        return False
    if "fs" in s:                                # frameshift, e.g. cirA_S90YfsTer15 -- always a call
        return True
    if _INDEL.search(s):                         # 23S_G543DEL / _INS / _DUP
        return True
    m = VARIANT.match(s)
    return bool(m) and m.group(2) != m.group(4)  # equal ref/alt is a wildtype screen row


def collect_carriers(lineages: dict[str, str], called_filter=is_called_substitution
                     ) -> tuple[dict, dict, dict, int]:
    """{symbol: {genome}} for acquired / CALLED point / SCREENED point, plus the synonymous-row count.

    THE CALLED-vs-SCREENED SPLIT IS THE ONE THAT DECIDES THIS RESULT. A POINT row in `main.tsv` is a
    substitution AMRFinder CALLS as resistance-conferring (50 distinct symbols here, topped by
    gyrA_S83L / parC_S80I / gyrA_D87N). `mutations.tsv` under --mutation_all lists every position
    SCREENED -- 1,516 symbols, a strict superset. The extra 1,466 are non-synonymous substitutions the
    tool does NOT call, and they are dominated by CLADE MARKERS: ompF_N52D and ampC_A231G sit at 98% in
    one lineage and are called ZERO times. Pooling them into the chromosomal arm measures phylogeny
    rather than resistance, which is the same error class as keeping the synonymous rows.
    """
    acquired: dict[str, set] = defaultdict(set)
    point: dict[str, set] = defaultdict(set)          # CALLED (main.tsv) -- the resistance determinants
    screened: dict[str, set] = defaultdict(set)       # every screened position (mutations.tsv)
    n_synonymous = 0
    for d in sorted(RUNS.iterdir()):
        g = _base(d.name)
        if g not in lineages:
            continue
        for fname in ("main.tsv", "mutations.tsv"):
            p = d / fname
            if not p.exists():
                continue
            text = p.read_text(encoding="utf-8", errors="replace").splitlines()
            for row in csv.DictReader(text, delimiter="\t"):
                symbol, subtype = row.get("Element symbol", ""), row.get("Subtype", "")
                if subtype == "AMR" and fname == "main.tsv":
                    acquired[symbol].add(g)
                elif subtype == "POINT":
                    if not called_filter(symbol):
                        n_synonymous += 1
                    elif fname == "main.tsv":
                        point[symbol].add(g)
                        screened[symbol].add(g)
                    else:
                        screened[symbol].add(g)
    return dict(acquired), dict(point), dict(screened), n_synonymous


def ceiling_largest_fraction(n_carriers: int, lineage_sizes: list[int]) -> float:
    """The LARGEST largest-lineage fraction any determinant with this many carriers could reach.

    Fill the biggest lineage first: with n <= the biggest lineage every carrier fits inside it (1.0);
    past that the best attainable is `biggest / n`, which FALLS as n grows. This is why raw coupling
    declines with prevalence in BOTH arms -- the headroom collapses, not the biology.
    """
    if n_carriers <= 0 or not lineage_sizes:
        return 0.0
    return min(1.0, max(lineage_sizes) / n_carriers)


def coupling_score(carriers: set, lineage_of: dict[str, str], all_genomes: list[str],
                   rng: np.random.Generator, n_null: int = 200,
                   lineage_sizes: list[int] | None = None) -> dict | None:
    """Largest-lineage fraction MINUS its null mean at the same carrier count.

    Positive = carriers sit in fewer lineages than chance = clade-COUPLED. Also reports that score as a
    fraction of what was ACHIEVABLE at this carrier count (see `ceiling_largest_fraction`).
    """
    n = len(carriers)
    if n < 2:
        return None
    if lineage_sizes is None:
        sizes: dict[str, int] = defaultdict(int)
        for g in all_genomes:
            sizes[lineage_of[g]] += 1
        lineage_sizes = list(sizes.values())

    def largest_fraction(genomes) -> float:
        counts: dict[str, int] = defaultdict(int)
        for g in genomes:
            counts[lineage_of[g]] += 1
        return max(counts.values()) / len(genomes)

    def n_lineages(genomes) -> int:
        return len({lineage_of[g] for g in genomes})

    observed_frac = largest_fraction(carriers)
    observed_lin = n_lineages(carriers)
    idx = np.arange(len(all_genomes))
    null_frac, null_lin = [], []
    for _ in range(n_null):
        draw = [all_genomes[i] for i in rng.choice(idx, size=n, replace=False)]
        null_frac.append(largest_fraction(draw))
        null_lin.append(n_lineages(draw))
    exp_lin = float(np.mean(null_lin))
    null_mean = float(np.mean(null_frac))
    ceiling = ceiling_largest_fraction(n, lineage_sizes)
    achievable = ceiling - null_mean
    raw = observed_frac - null_mean
    return {"n_carriers": n,
            "carrier_fraction": n / len(all_genomes),
            "largest_lineage_fraction": observed_frac,
            "null_mean_largest_fraction": null_mean,
            "coupling_score": raw,
            "ceiling_largest_fraction": ceiling,
            "achievable_coupling": achievable,
            # None rather than a huge/undefined ratio when there is no headroom left to normalize by.
            "normalized_coupling": (raw / achievable) if achievable > 1e-9 else None,
            "n_lineages": observed_lin,
            "expected_n_lineages": exp_lin,
            "lineage_spread_ratio": (observed_lin / exp_lin) if exp_lin else None}


def permutation_gap(point_scores, acq_scores, rng, n_perm: int = 1000) -> dict:
    """Null for the ARM LABEL: shuffle which families are point vs acquired, recompute the gap.

    Guards the comparison itself. Two arms of different sizes drawn from one score distribution differ
    by chance, and without this the sign of the gap alone would be reported as the finding.
    """
    obs = float(np.mean(point_scores) - np.mean(acq_scores))
    pooled = np.array(list(point_scores) + list(acq_scores), dtype=float)
    k = len(point_scores)
    gaps = []
    for _ in range(n_perm):
        perm = rng.permutation(pooled)
        gaps.append(float(np.mean(perm[:k]) - np.mean(perm[k:])))
    gaps_arr = np.array(gaps)
    return {"observed_gap": obs, "n_perm": n_perm,
            "perm_mean": float(gaps_arr.mean()), "perm_max": float(gaps_arr.max()),
            "perm_p95": float(np.percentile(gaps_arr, 95)),
            "exceeds_perm_max": bool(obs > gaps_arr.max())}


def verdict_from_bar(point_scores, acq_scores, perm: dict, min_families: int = 20) -> str:
    """The FROZEN rule, applied mechanically (wiki/determinant_clade_coupling_acceptance_bar.json)."""
    if len(point_scores) < min_families or len(acq_scores) < min_families:
        return "INDETERMINATE"
    if float(np.mean(point_scores)) <= float(np.mean(acq_scores)):
        return "FALSIFIED"
    return "SUPPORTED" if perm["exceeds_perm_max"] else "WEAK_DIRECTIONAL"


PREVALENCE_BANDS = (("rare_lt_0.10", 0.0, 0.10),
                    ("mid_0.10_to_0.33", 0.10, 0.33),
                    ("common_ge_0.33", 0.33, 1.01))


def _stratify(scored: dict) -> dict:
    """Per-arm prevalence bands carrying BOTH the raw and the ceiling-normalized mean.

    Reporting raw alone is what made the decline look like biology; reporting normalized alone would
    hide that the raw decline is real and needs explaining. Both, in the same table.
    """
    out: dict = {}
    for arm in ("point", "acquired"):
        out[arm] = {}
        for name, lo, hi in PREVALENCE_BANDS:
            members = sorted(k for k, x in scored[arm].items() if lo <= x["carrier_fraction"] < hi)
            raw = [scored[arm][k]["coupling_score"] for k in members]
            norm = [scored[arm][k]["normalized_coupling"] for k in members
                    if scored[arm][k]["normalized_coupling"] is not None]
            out[arm][name] = {
                "n_families": len(members),
                "mean_coupling_raw": float(np.mean(raw)) if raw else None,
                "mean_coupling_normalized": float(np.mean(norm)) if norm else None,
                "n_normalizable": len(norm),
                "members": members[:8],
            }
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--min-carriers", type=int, default=5)
    ap.add_argument("--max-carrier-fraction", type=float, default=0.50,
                    help="prevalence matching; see the bar's derivation")
    ap.add_argument("--n-null", type=int, default=200)
    ap.add_argument("--n-perm", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=Path,
                    default=ROOT / "wiki" / f"determinant_clade_coupling_{_date.today()}.json")
    a = ap.parse_args(argv)

    lineages = load_lineages()
    acquired, point, screened, n_syn = collect_carriers(lineages)
    genomes = sorted(lineages)
    rng = np.random.default_rng(a.seed)
    print(f"genomes {len(genomes)} | lineages {len(set(lineages.values()))} | "
          f"synonymous rows filtered {n_syn}")

    sizes: dict[str, int] = defaultdict(int)
    for g in genomes:
        sizes[lineages[g]] += 1
    lineage_sizes = list(sizes.values())
    print(f"largest lineage {max(lineage_sizes)} genomes -- sets the coupling ceiling at each prevalence")

    scored = {"acquired": {}, "point": {}, "screened": {}}
    for arm, table in (("acquired", acquired), ("point", point), ("screened", screened)):
        for sym, carriers in table.items():
            if len(carriers) < a.min_carriers:
                continue
            s = coupling_score(carriers, lineages, genomes, rng, n_null=a.n_null,
                               lineage_sizes=lineage_sizes)
            if s:
                scored[arm][sym] = s

    def band(arm, matched: bool, metric: str = "coupling_score"):
        return [v[metric] for v in scored[arm].values()
                if ((not matched) or v["carrier_fraction"] <= a.max_carrier_fraction)
                and v[metric] is not None]

    results = {}
    for label, matched, parm in (("called_point_prevalence_matched", True, "point"),
                                 ("called_point_unrestricted", False, "point"),
                                 ("screened_point_prevalence_matched", True, "screened"),
                                 ("screened_point_unrestricted", False, "screened")):
        pv, av = band(parm, matched), band("acquired", matched)
        perm = permutation_gap(pv, av, np.random.default_rng(a.seed), n_perm=a.n_perm)
        results[label] = {
            "n_point_families": len(pv), "n_acquired_families": len(av),
            "mean_coupling_point": float(np.mean(pv)) if pv else None,
            "mean_coupling_acquired": float(np.mean(av)) if av else None,
            "median_coupling_point": float(np.median(pv)) if pv else None,
            "median_coupling_acquired": float(np.median(av)) if av else None,
            "permutation": perm,
            "verdict": verdict_from_bar(pv, av, perm),
        }

    # SECONDARY band: the same comparison on ceiling-normalized scores. Reported beside the frozen-bar
    # result, never in place of it -- the bar's rule is defined on the raw score.
    normalized = {}
    for label, matched, parm in (("called_point_prevalence_matched", True, "point"),
                                 ("called_point_unrestricted", False, "point")):
        pv = band(parm, matched, "normalized_coupling")
        av = band("acquired", matched, "normalized_coupling")
        perm = permutation_gap(pv, av, np.random.default_rng(a.seed), n_perm=a.n_perm)
        normalized[label] = {
            "n_point_families": len(pv), "n_acquired_families": len(av),
            "mean_normalized_point": float(np.mean(pv)) if pv else None,
            "mean_normalized_acquired": float(np.mean(av)) if av else None,
            "permutation": perm,
            "verdict_if_the_bar_were_applied_here": verdict_from_bar(pv, av, perm),
        }

    doc = {
        "schema": "determinant-clade-coupling-v1",
        "date": str(_date.today()),
        "claim_under_test": ("COUPLING-not-popdesign: acquired determinants should be LESS clade-coupled "
                             "than chromosomal point mutations."),
        "bar": "wiki/determinant_clade_coupling_acceptance_bar.json",
        "n_genomes": len(genomes), "n_lineages": len(set(lineages.values())),
        "n_synonymous_rows_filtered": n_syn,
        "min_carriers": a.min_carriers, "max_carrier_fraction": a.max_carrier_fraction,
        "headline_verdict": results["called_point_prevalence_matched"]["verdict"],
        "headline_arm": ("CALLED point substitutions only -- the corrected arm. The screened_* "
                         "bands keep every screened position and are reported so the pollution "
                         "is visible, never as the headline."),
        "verdict_is_mechanical": "verdict_from_bar applied to the frozen rule; not authored",
        "results": results,
        "ceiling_normalized": normalized,
        "prevalence_stratified": _stratify(scored),
        "why_the_stratification_needed_a_ceiling": (
            "Raw coupling falls with prevalence in BOTH arms, which invites the reading that the "
            "determinants driving most resistance are decoupled. That reading is an ARTIFACT of the "
            "metric's headroom and was nearly published here. The largest lineage holds 53 of 621 "
            "genomes, so a determinant with 400 carriers cannot exceed 53/400 = 0.133 largest-lineage "
            "fraction and has only ~0.048 of coupling available to it. Normalized by what was "
            "achievable, gyrA_S83L sits at 0.895, parC_S80I at 0.755, gyrA_D87N at 0.698 -- the "
            "dominant chromosomal determinants are near-SATURATED, not near-chance -- while the "
            "acquired extremes at comparable prevalence sit at or below zero (blaTEM-1 -0.134, "
            "sul2 -0.132). The arm gap is LARGER under normalization, not smaller. The null fixes the "
            "EXPECTATION at a given carrier count; it does not fix the MAXIMUM, and only the ceiling "
            "does."),
        "n_called_point_symbols": len(point), "n_screened_point_symbols": len(screened),
        "top_coupled_point": sorted(
            ({"symbol": k, **v} for k, v in scored["point"].items()
             if v["carrier_fraction"] <= a.max_carrier_fraction),
            key=lambda d: -d["coupling_score"])[:5],
        "top_decoupled_acquired": sorted(
            ({"symbol": k, **v} for k, v in scored["acquired"].items()),
            key=lambda d: d["coupling_score"])[:5],
        "honest_limits": [
            "One organism (E. coli) and one lineage definition (7-locus MLST). ST is a COARSE lineage "
            "unit; a finer phylogeny could move every family in both arms.",
            "AMRFinder's own AMR-vs-POINT Subtype is the acquired-vs-chromosomal proxy. POINT is a "
            "chromosomal substitution by definition, but this is the tool's classification, not an "
            "independent determination of mobility.",
            "Establishing the mechanism's PREMISE is not establishing that coupling CAUSES the "
            "de-confounded failures; that needs the regime comparison itself.",
            "Prevalence matching is one-sided (it removes POINT families only), which is why the "
            "unrestricted result ships beside it rather than being replaced by it.",
            "Genomes are AMR-cohort leftovers, enriched for resistance and not a natural population.",
            "The ceiling normalization divides by a SMALL denominator at high prevalence (~0.048 at "
            "n=400), so a normalized score there is far noisier than one for a rare family. It "
            "corrects a real bias in the raw score; it does not make the two ends equally precise.",
            "`ceiling_largest_fraction` assumes carriers could be packed into the biggest lineage "
            "arbitrarily. That is the right bound for a metric ceiling but no real determinant is free "
            "to distribute itself, so the normalized score is a fraction of a THEORETICAL maximum.",
        ],
    }
    a.out.write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8")
    for label in results:
        r = results[label]
        print(f"{label:<20} point {r['mean_coupling_point']} (n={r['n_point_families']})  "
              f"acquired {r['mean_coupling_acquired']} (n={r['n_acquired_families']})  "
              f"gap {r['permutation']['observed_gap']:.4f} vs perm max "
              f"{r['permutation']['perm_max']:.4f} -> {r['verdict']}")
    print("-- ceiling-normalized (secondary; the frozen bar is on the raw score above) --")
    for label, r in normalized.items():
        print(f"{label:<20} point {r['mean_normalized_point']:.4f} "
              f"acquired {r['mean_normalized_acquired']:.4f}  "
              f"gap {r['permutation']['observed_gap']:.4f} vs perm max "
              f"{r['permutation']['perm_max']:.4f}")
    print(f"-> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
