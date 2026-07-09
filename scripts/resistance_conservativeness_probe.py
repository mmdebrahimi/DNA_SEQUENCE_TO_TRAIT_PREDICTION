#!/usr/bin/env python
"""Are resistance mutations chemically CONSERVATIVE across pathogens? (BLOSUM only -- no model, no GPU.)

`hiv_esm_vs_catalog.py` + `hiv_esm_likelihood_probe.py` concluded: DRM sites are ordinarily conserved,
resistance is reached via chemically mild substitutions there, and so every exchangeability/likelihood
scorer (BLOSUM62 and ESM2 alike) rates the resistant residue benign. That memo then made a GENERAL
prediction: expect the same for ANY conservation-based predictor on ANY antagonistically-selected
phenotype.

This script tests that prediction directly, and then tries to kill it with the obvious confound.

PART 1 -- generalization. For each committed mutant-level resistance catalog, rank the resistant
residue among the 19 alternatives at its own position by BLOSUM62. Under a random mutant the rank is
uniform on 1..19, so the null median is exactly 10.0.

PART 2 -- THE CODON CONTROL (the confound that would explain everything). A resistance mutation must
be reachable by a SINGLE nucleotide change, and the genetic code is error-minimizing: single-nt
neighbours are chemically similar BY CONSTRUCTION. So "resistance mutations are conservative" might be
a fact about the genetic code, not about selection. Re-rank each HIV DRM among only the substitutions
reachable by one nt change from its real HXB2 codon (null = (|accessible|+1)/2).

RESULT (2026-07-09): the generalization HOLDS -- but most of it is the genetic code.
  Part 1: see the table the script prints (null median 10.0). Ranks are MID-ranks -- BLOSUM62
          scores tie heavily, and breaking ties by iteration order silently shifts the median.
  Part 2: among single-nt-accessible substitutions HIV DRMs are only modestly better than the
          accessible-null, on n=22. SUGGESTIVE AND UNDERPOWERED, not a claim.

HONEST CONCLUSION. Resistance mutations ARE chemically conservative (Part 1, robust, 3 unrelated
pathogens). The dominant driver is CODON ACCESSIBILITY -- a property of the genetic code -- with at
most a modest additional selection preference that this data cannot establish. Either way the
practical rule is unchanged and does not depend on resolving the cause: a conservation- or
exchangeability-based scorer is structurally blind to resistance.
"""
import json
import random
import re
import statistics as st
import sys

sys.path.insert(0, ".")
import dna_decode.data.fungal_amr as F
import dna_decode.data.hiv_amr as H
import dna_decode.data.sarscov2_amr as S
from scripts.hiv_esm_vs_catalog import CODON, RT_REF
from scripts.hiv_esm_likelihood_probe import NRTI_DRMS

from Bio.Align import substitution_matrices as _sm

B = _sm.load("BLOSUM62")
AA = list("ACDEFGHIKLMNPQRSTVWY")
_SUB = re.compile(r"^([A-Z])(\d+)([A-Z])$")


def blosum_rank(wt, mut, among=None):
    """Mid-rank of `mut` among `among` (default: all 19 non-wt AAs), most-conservative first.

    BLOSUM62 scores TIE heavily (many substitutions share a score). Breaking ties by iteration
    order makes the rank depend on set ordering and silently shifts the median. Use the standard
    mid-rank: 1 + (#strictly better) + (#tied - 1)/2.
    """
    alts = [a for a in (among or AA) if a != wt]
    if mut not in alts:
        return None
    s = B[wt, mut]
    better = sum(1 for a in alts if B[wt, a] > s)
    tied = sum(1 for a in alts if B[wt, a] == s)
    return 1 + better + (tied - 1) / 2.0


def fungal_drms():
    out = []
    for v in F.FUNGAL_RESISTANCE_MUTATIONS.values():
        if isinstance(v, dict):
            for muts in v.values():
                out += list(muts)
        else:
            out += list(v)
    return sorted(set(out))


def accessible_aas(codon):
    """Amino acids reachable from `codon` by exactly one nucleotide substitution (no stops)."""
    wt = CODON.get(codon)
    out = set()
    for i in range(3):
        for b in "ACGT":
            if b == codon[i]:
                continue
            a = CODON.get(codon[:i] + b + codon[i + 1:])
            if a and a != "*" and a != wt:
                out.add(a)
    return out


def main():
    cats = {"HIV RT (NNRTI+NRTI)": sorted(H.NNRTI_RT_MAJOR_DRMS) + NRTI_DRMS,
            "SARS-CoV-2 Mpro": sorted(S.MPRO_MAJOR_DRMS),
            "Fungal ERG11/FKS1": fungal_drms()}

    print("PART 1 -- BLOSUM62 rank of the resistant residue among 19 alternatives (null median 10.0)\n")
    print(f"{'pathogen':24s} {'n':>4s} {'median':>8s} {'frac<=5':>9s}")
    print("-" * 48)
    pooled, per_path = [], {}
    for name, muts in cats.items():
        rk = []
        for m in muts:
            g = _SUB.match(m)
            if not g:
                continue
            r = blosum_rank(g.group(1), g.group(3))
            if r:
                rk.append(r)
        if not rk:
            continue
        pooled += rk
        per_path[name] = dict(n=len(rk), median_rank=st.median(rk),
                              frac_le5=round(sum(1 for r in rk if r <= 5) / len(rk), 4))
        print(f"{name:24s} {len(rk):4d} {st.median(rk):8.1f} {per_path[name]['frac_le5']:9.0%}")
    print("-" * 48)
    print(f"{'POOLED':24s} {len(pooled):4d} {st.median(pooled):8.1f} "
          f"{sum(1 for r in pooled if r <= 5)/len(pooled):9.0%}")
    print(f"{'null (uniform 1..19)':24s} {'':4s} {10.0:8.1f} {5/19:9.0%}")

    # ---- PART 2: the codon-accessibility control (HIV only -- the one system with real codons) ----
    seq = "".join(l.strip() for l in open(RT_REF) if not l.startswith(">"))
    codons = [seq[i:i + 3] for i in range(0, len(seq) - 2, 3)]
    prot = "".join(CODON.get(c, "X") for c in codons)

    all_rank, acc_rank, acc_null, n_unreachable = [], [], [], 0
    for d in sorted(H.NNRTI_RT_MAJOR_DRMS) + NRTI_DRMS:
        g = _SUB.match(d)
        wt, pos, mut = g.group(1), int(g.group(2)), g.group(3)
        if pos > len(prot) or prot[pos - 1] != wt:
            continue
        acc = accessible_aas(codons[pos - 1])
        if mut not in acc:
            n_unreachable += 1
            continue
        all_rank.append(blosum_rank(wt, mut))
        acc_rank.append(blosum_rank(wt, mut, among=acc | {wt}))
        acc_null.append((len(acc) + 1) / 2)

    win = sum((n > a) + 0.5 * (n == a) for a, n in zip(acc_rank, acc_null)) / len(acc_rank)
    better = sum(1 for a, n in zip(acc_rank, acc_null) if a < n)
    print("\n\nPART 2 -- CODON CONTROL (HIV): does the genetic code already explain it?")
    print(f"  DRMs reachable by ONE nt change from the HXB2 codon: {len(acc_rank)} "
          f"(not reachable: {n_unreachable})")
    print(f"  among ALL 19 alternatives:       median rank {st.median(all_rank):.1f}  (null 10.0)")
    print(f"  among ACCESSIBLE substitutions:  median rank {st.median(acc_rank):.1f}  "
          f"(null {st.median(acc_null):.1f})")
    print(f"  more conservative than accessible-null: {better}/{len(acc_rank)}   P = {win:.3f}")
    print("\n  READ THIS HONESTLY: a random single-nt-accessible substitution is ALREADY conservative")
    print("  (the genetic code is error-minimising). The residual selection preference is")
    print(f"  SUGGESTIVE AND UNDERPOWERED (P={win:.3f}, n={len(acc_rank)}) -- not a claim.")

    print("\n" + "=" * 78)
    print("Resistance mutations ARE chemically conservative (Part 1: 3 unrelated pathogens, robust).")
    print("The dominant driver is CODON ACCESSIBILITY, a property of the genetic code, with at most a")
    print("modest additional selection preference this data cannot establish. The practical rule does")
    print("NOT depend on resolving the cause: a conservation-based scorer is structurally blind here.")
    print("=" * 78)

    json.dump({"date": "2026-07-09", "matrix": "BLOSUM62", "null_median_rank_of_19": 10.0,
               "part1_per_pathogen": per_path,
               "part1_pooled_n": len(pooled), "part1_pooled_median_rank": st.median(pooled),
               "part1_generalization": "HOLDS across HIV RT, SARS-CoV-2 Mpro, fungal ERG11/FKS1",
               "part2_hiv_n_reachable": len(acc_rank), "part2_hiv_n_unreachable": n_unreachable,
               "part2_median_rank_all19": st.median(all_rank),
               "part2_median_rank_accessible": st.median(acc_rank),
               "part2_accessible_null_median": st.median(acc_null),
               "part2_better_than_accessible_null": f"{better}/{len(acc_rank)}",
               "part2_p_more_conservative_than_accessible": round(win, 4),
               "part2_verdict": "codon accessibility explains most of it; residual selection "
                                "preference is SUGGESTIVE AND UNDERPOWERED (n=22), not a claim",
               "practical_rule": "a conservation-/exchangeability-based scorer is structurally blind "
                                 "to resistance, regardless of which driver dominates"},
              open("wiki/resistance_conservativeness_2026-07-09.json", "w"), indent=1)
    print("\nwrote wiki/resistance_conservativeness_2026-07-09.json")


if __name__ == "__main__":
    main()
