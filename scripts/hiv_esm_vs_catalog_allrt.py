#!/usr/bin/env python
"""Generalize the ESM-vs-catalog negative across every RT-targeting HIV drug (NNRTI + NRTI).

Both classes target the SAME protein (reverse transcriptase), so the ESM2-650M masked-marginal
matrix cached by `scripts/hiv_esm_vs_catalog.py` scores all of them for free. This turns the
single-drug EFV result into an 11-drug, 2-class result at near-zero compute.

Comparator is the DEPLOYED rule (`hiv_amr.call_hiv_observed`) -- verified to agree with a
hand-reimplementation on 2168/2168 EFV isolates before this script was written.

Cutoff: a UNIFORM ILLUSTRATIVE fold >= 3x for every drug. No per-drug clinical cutoff exists
in-repo for these classes, and the project's PI v0.1 work set the precedent of scoring all arms
at one illustrative cutoff and reporting the COMPARISON, not an absolute-calibration claim.

Reports per drug:
  catalog AUROC (deployed rule, binary)      -- the incumbent
  ESM AUROC on the FULL cohort               -- POSITIVE CONTROL: is the LM sensitive at all?
  ESM AUROC on the CATALOG-NEGATIVE subset   -- the question: does it help where the catalog is blind?
  mutation-burden AUROC on that subset       -- the confounder
"""
import json
import statistics as st
import sys

sys.path.insert(0, ".")
import dna_decode.data.hiv_amr as H
from scripts.hiv_esm_vs_catalog import (AA, LOGP_CACHE, auroc, isolate_muts, rt_protein)

import csv

DRIFT = {122, 214, 272, 277}          # HXB2 != consensus B (derived in hiv_esm_vs_catalog.py)
CUTOFF = 3.0
ARMS = [
    ("NNRTI", "data/raw/hiv/NNRTI_DataSet.txt",
     [("EFV", "efavirenz"), ("NVP", "nevirapine"), ("ETR", "etravirine"),
      ("RPV", "rilpivirine"), ("DOR", "doravirine")]),
    ("NRTI", "data/raw/hiv/NRTI_DataSet.txt",
     [("3TC", "lamivudine"), ("ABC", "abacavir"), ("AZT", "zidovudine"),
      ("D4T", "stavudine"), ("DDI", "didanosine"), ("TDF", "tenofovir")]),
]


def main():
    prot = rt_protein()
    logp = {int(k): v for k, v in json.load(open(LOGP_CACHE)).items()}
    print(f"cached masked-marginals: {len(logp)} RT positions\n")

    def subs_of(r, pcols):
        return {f"{prot[p-1]}{p}{aa}" for p, aa in isolate_muts(r, pcols)
                if p <= len(prot) and p not in DRIFT and prot[p - 1] != aa}

    def esm_score(r, pcols):
        d = [logp[p][prot[p - 1]] - logp[p][aa] for p, aa in isolate_muts(r, pcols)
             if p in logp and p <= len(prot) and prot[p - 1] != aa and aa in AA]
        return st.mean(d) if d else 0.0

    print(f"{'drug':6s} {'class':6s} {'n':>5s} {'R':>5s} | {'catalog':>8s} {'ESMfull':>8s} | "
          f"{'sub_n':>6s} {'sub_R':>5s} {'ESMsub':>7s} {'burden':>7s}")
    print("-" * 82)
    out = []
    for cls, path, drugs in ARMS:
        rows = list(csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"))
        pcols = [c for c in rows[0] if c.startswith("P") and c[1:].isdigit()]
        for col, drug in drugs:
            if col not in rows[0]:
                print(f"{col:6s} {cls:6s}  -- column absent")
                continue
            have = [r for r in rows if r[col] not in ("NA", "", "-")]
            if len(have) < 50:
                print(f"{col:6s} {cls:6s} n={len(have)} too few")
                continue
            y = [1 if float(r[col]) >= CUTOFF else 0 for r in have]
            cat = [1 if H.call_hiv_observed(drug, {"RT": subs_of(r, pcols)}).prediction == "R" else 0
                   for r in have]
            esm = [esm_score(r, pcols) for r in have]
            a_cat, a_esm_full = auroc(y, cat), auroc(y, esm)

            sub = [r for r, c in zip(have, cat) if not c]
            ys = [1 if float(r[col]) >= CUTOFF else 0 for r in sub]
            if sum(ys) < 10 or len(ys) - sum(ys) < 10:
                print(f"{col:6s} {cls:6s} {len(have):5d} {sum(y):5d} | {a_cat:8.3f} {a_esm_full:8.3f} | "
                      f"{len(sub):6d} {sum(ys):5d}   UNDERPOWERED")
                continue
            a_esm_sub = auroc(ys, [esm_score(r, pcols) for r in sub])
            a_burden = auroc(ys, [len(isolate_muts(r, pcols)) for r in sub])
            print(f"{col:6s} {cls:6s} {len(have):5d} {sum(y):5d} | {a_cat:8.3f} {a_esm_full:8.3f} | "
                  f"{len(sub):6d} {sum(ys):5d} {a_esm_sub:7.3f} {a_burden:7.3f}")
            out.append(dict(drug=col, drug_full=drug, cls=cls, n=len(have), n_r=sum(y),
                            catalog_auroc=round(a_cat, 4), esm_full_auroc=round(a_esm_full, 4),
                            subset_n=len(sub), subset_r=sum(ys),
                            esm_subset_auroc=round(a_esm_sub, 4), burden_auroc=round(a_burden, 4)))

    full = [r["esm_full_auroc"] for r in out]
    subs = [r["esm_subset_auroc"] for r in out]
    cats = [r["catalog_auroc"] for r in out]
    print("-" * 82)
    print(f"median over {len(out)} drugs:  catalog {st.median(cats):.3f}   "
          f"ESM-full {st.median(full):.3f}   ESM-subset {st.median(subs):.3f}")
    print(f"drugs where ESM-full > 0.60 (any sensitivity to resistance): "
          f"{sum(1 for x in full if x > 0.60)}/{len(full)}")
    print(f"drugs where ESM-subset >= 0.65 (the PASS bar):              "
          f"{sum(1 for x in subs if x >= 0.65)}/{len(subs)}")

    # A raw ">=0.65" count is a trap. Require the lift to be real: beat the burden confounder by a
    # margin, on a subset that is actually powered, against a catalog that is itself informative.
    def genuine(r):
        return (r["esm_subset_auroc"] >= 0.65
                and r["esm_subset_auroc"] - r["burden_auroc"] >= 0.05
                and r["subset_r"] >= 20 and r["subset_n"] - r["subset_r"] >= 20
                and r["catalog_auroc"] >= 0.65)
    real = [r["drug"] for r in out if genuine(r)]
    nominal = [r["drug"] for r in out if r["esm_subset_auroc"] >= 0.65]
    print(f"  of those, GENUINE (beats burden by >=0.05, >=20/class, catalog informative): "
          f"{len(real)}/{len(subs)}  {real}")
    for r in out:
        if r["drug"] in nominal and r["drug"] not in real:
            print(f"  !! {r['drug']} clears 0.65 NOMINALLY only: esm={r['esm_subset_auroc']:.3f} vs "
                  f"burden={r['burden_auroc']:.3f} (lift {r['esm_subset_auroc']-r['burden_auroc']:+.3f}), "
                  f"n={r['subset_n']} R={r['subset_r']}, catalog={r['catalog_auroc']:.3f}")

    json.dump({"date": "2026-07-09", "model": "facebook/esm2_t33_650M_UR50D",
               "cutoff_fold": CUTOFF, "cutoff_status": "UNIFORM_ILLUSTRATIVE (no per-drug clinical "
               "cutoff in-repo for these classes; comparison, not absolute calibration)",
               "comparator": "deployed hiv_amr.call_hiv_observed",
               "per_drug": out,
               "median_catalog_auroc": round(st.median(cats), 4),
               "median_esm_full_auroc": round(st.median(full), 4),
               "median_esm_subset_auroc": round(st.median(subs), 4),
               "n_drugs_esm_sensitive": sum(1 for x in full if x > 0.60),
               "n_drugs_passing_bar_nominal": sum(1 for x in subs if x >= 0.65),
               "n_drugs_passing_bar_genuine": len(real),
               "genuine_pass_drugs": real,
               "genuine_criteria": "esm_subset>=0.65 AND esm_subset-burden>=0.05 AND >=20 per class "
                                  "AND catalog_auroc>=0.65",
               "nominal_pass_caveat": "a raw >=0.65 count is a trap: DOR clears it on n=37 with a "
                                     "burden baseline of 0.783 and a near-chance catalog (0.549)."},
              open("wiki/hiv_esm_vs_catalog_allrt_2026-07-09.json", "w"), indent=1)
    print("\nwrote wiki/hiv_esm_vs_catalog_allrt_2026-07-09.json")


if __name__ == "__main__":
    main()
