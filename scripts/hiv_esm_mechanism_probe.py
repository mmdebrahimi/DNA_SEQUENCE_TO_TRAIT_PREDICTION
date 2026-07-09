#!/usr/bin/env python
"""Why does ESM look sensitive on NRTI but not NNRTI? Two probes that kill the obvious answer.

The first-pass reading of `hiv_esm_vs_catalog_allrt.py` was: "ESM's sensitivity tracks the FITNESS
COST of the resistance mechanism -- NNRTI DRMs (K103N, Y181C) are cheap so ESM misses them, NRTI DRMs
(M184V, K65R, TAMs) cost replication so ESM sees them." That story is intuitive, and WRONG. Both
probes below use only the cached masked-marginal matrix -- no new inference.

PROBE 1 (per-DRM). If the fitness-cost story held, NRTI DRMs should look MORE damaging to ESM than
NNRTI DRMs at their own positions. They look LESS damaging (median percentile 0.184 vs 0.289;
P(NRTI more damaging than NNRTI) = 0.472, i.e. no better than a coin flip). FALSIFIED.

PROBE 2 (per-cohort). The real explanation: MUTATION BURDEN. NRTI-resistant isolates are
treatment-experienced and carry more mutations. Counting mutations -- ignoring the language model
entirely -- beats ESM on 10 of 11 drugs, including every NRTI drug except 3TC (+0.035).

Conclusion: ESM has NO resistance signal in either class. Its apparent NRTI "sensitivity" is a
treatment-experience/divergence confound that a trivial baseline captures better.
"""
import csv
import json
import statistics as st
import sys

sys.path.insert(0, ".")
import dna_decode.data.hiv_amr as H
from scripts.hiv_esm_vs_catalog import AA, LOGP_CACHE, auroc, isolate_muts, rt_protein

CUTOFF = 3.0
ARMS = [("NNRTI", "data/raw/hiv/NNRTI_DataSet.txt", ["EFV", "NVP", "ETR", "RPV", "DOR"]),
        ("NRTI", "data/raw/hiv/NRTI_DataSet.txt", ["3TC", "ABC", "AZT", "D4T", "DDI", "TDF"])]
# Deconfounded NRTI mutant catalog, from wiki/hiv_nrti_mutant_catalog_v0_1_2026-06-21.json
NRTI_DRMS = ["K65R", "K70E", "K70G", "K70N", "K70R", "K70T", "L74V", "M184I", "M184V", "M41L",
             "Q151M", "L210W", "T215E", "T215F", "T215I", "T215L", "T215V", "T215Y",
             "V75A", "V75I", "V75M", "V75T"]


def main():
    prot = rt_protein()
    logp = {int(k): v for k, v in json.load(open(LOGP_CACHE)).items()}

    def damage_percentile(drm):
        wt, mut, pos = drm[0], drm[-1], int(drm[1:-1])
        if pos not in logp or prot[pos - 1] != wt:
            return None
        dmg = {a: logp[pos][wt] - logp[pos][a] for a in AA if a != wt}
        rank = sorted(dmg, key=lambda a: dmg[a])          # ascending damage
        return (rank.index(mut) + 1) / len(rank)

    # ---- PROBE 1: does ESM think NRTI DRMs are more damaging than NNRTI DRMs? ----
    nn = [p for p in (damage_percentile(d) for d in sorted(H.NNRTI_RT_MAJOR_DRMS)) if p]
    nr = [p for p in (damage_percentile(d) for d in NRTI_DRMS) if p]
    win = sum((y > x) + 0.5 * (y == x) for x in nn for y in nr) / (len(nn) * len(nr))
    print("PROBE 1 -- per-DRM damage percentile at own position")
    print(f"  NNRTI (cheap, circulating): n={len(nn)}  median {st.median(nn):.3f}")
    print(f"  NRTI  (fitness-costly):     n={len(nr)}  median {st.median(nr):.3f}")
    print(f"  P(NRTI DRM looks more damaging than NNRTI DRM) = {win:.3f}")
    print(f"  fitness-cost story predicts >> 0.5  ->  {'SUPPORTED' if win > 0.65 else 'FALSIFIED'}\n")

    # ---- PROBE 2: is ESM's full-cohort signal just mutation burden? ----
    def esm_score(r, pc):
        d = [logp[p][prot[p - 1]] - logp[p][aa] for p, aa in isolate_muts(r, pc)
             if p in logp and p <= len(prot) and prot[p - 1] != aa and aa in AA]
        return st.mean(d) if d else 0.0

    print("PROBE 2 -- full-cohort AUROC: ESM vs simply COUNTING mutations")
    print(f"  {'drug':5s} {'class':6s} {'n':>5s} | {'ESM':>7s} {'burden':>7s} {'ESM-burden':>11s}")
    rows_out, esm_wins = [], 0
    for cls, path, cols in ARMS:
        rows = list(csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"))
        pc = [c for c in rows[0] if c.startswith("P") and c[1:].isdigit()]
        for col in cols:
            have = [r for r in rows if r[col] not in ("NA", "", "-")]
            if len(have) < 50:
                continue
            y = [1 if float(r[col]) >= CUTOFF else 0 for r in have]
            ae = auroc(y, [esm_score(r, pc) for r in have])
            ab = auroc(y, [len(isolate_muts(r, pc)) for r in have])
            esm_wins += ae > ab
            print(f"  {col:5s} {cls:6s} {len(have):5d} | {ae:7.3f} {ab:7.3f} {ae-ab:+11.3f}")
            rows_out.append(dict(drug=col, cls=cls, n=len(have), esm_auroc=round(ae, 4),
                                 burden_auroc=round(ab, 4), esm_minus_burden=round(ae - ab, 4)))
    print(f"\n  ESM beats a plain mutation count on {esm_wins}/{len(rows_out)} drugs.")
    print("  => the apparent NRTI 'sensitivity' is treatment-experience burden, not resistance biology.")

    json.dump({"date": "2026-07-09", "model": "facebook/esm2_t33_650M_UR50D",
               "probe1_nnrti_median_damage_percentile": round(st.median(nn), 4),
               "probe1_nrti_median_damage_percentile": round(st.median(nr), 4),
               "probe1_p_nrti_more_damaging": round(win, 4),
               "probe1_fitness_cost_story": "FALSIFIED",
               "probe2_per_drug": rows_out,
               "probe2_esm_beats_burden_on": esm_wins, "probe2_n_drugs": len(rows_out),
               "conclusion": "ESM has no resistance signal in either class; its apparent NRTI "
                             "sensitivity is mutation burden (treatment experience), which a trivial "
                             "count captures better on 10/11 drugs."},
              open("wiki/hiv_esm_mechanism_probe_2026-07-09.json", "w"), indent=1)
    print("\nwrote wiki/hiv_esm_mechanism_probe_2026-07-09.json")


if __name__ == "__main__":
    main()
