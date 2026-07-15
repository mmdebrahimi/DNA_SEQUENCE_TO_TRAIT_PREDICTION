"""Pre-resistance PREDICTIVE backtest — the real test that separates a codon-table lookup from a forecaster.

The /idea-validation-council's decisive question (Contrarian + First-Principles + Competitor): does the
"one-nt-away" flag have PREDICTIVE content, or is it a static restatement of the genetic code? The only clean
way to answer is a temporal holdout: compute the flags on catalogue v(N) and ask whether the DRMs ADDED in
v(N+1) are enriched at the flagged codons ABOVE the base rate of all adjacent codons.

Substrate: the WHO M. tuberculosis mutation catalogue, which shipped v1 (2021) -> v2 (2023) with the diff
ENCODED in the committed v2 master file (`CHANGES vs ver1` column: New AwR / New AwRI / UP from Uncertain to
AwR|AwRI = mutations newly graded resistance in v2; the v1-R set = the seed the forecaster "knew" at v1).
Deterministic, offline, no network, no GPU. Frozen decoder surface untouched.

The forecaster's structural ceiling (stated honestly): a "one nt from a known resistance codon" flag can ONLY
fire at a position that ALREADY had a v1 resistance allele — it can NEVER predict a brand-new mechanism at a
new position. So we report BOTH:
  (1) CEILING  — what fraction of v2-added-R mutations even occur at a v1-R position (the max the forecaster
                 could possibly catch); the rest are unreachable by construction.
  (2) ENRICHMENT (the real test) — a position-matched Fisher exact: among protein substitutions AT v1-R
                 positions, are the v2-added-R ones MORE one-nt-adjacent to the v1-R allele than the background
                 of non-added substitutions at those SAME positions? This controls for the "new DRMs cluster
                 at known hotspots" confound — hotspot membership alone is NOT credited; only codon-adjacency.

Aa-level adjacency (census method): a substitution is "one nt from the v1-R codon" iff SOME codon of the new
residue is one nt from SOME codon of a v1-R residue at that position. This is a GENEROUS upper bound (isolate
exact codon unknown at aa resolution) — so NO enrichment here is a STRONG fail (even the generous test finds
nothing); enrichment would flag codon-exact confirmation (H37Rv reference codons) as the follow-up.
"""
from __future__ import annotations

import csv
import json
import re
import sys
from datetime import date as _date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
MASTER = REPO / "data" / "raw" / "who_tb_catalogue" / "WHO-UCN-TB-2023.6-eng_catalogue_master_file.txt"

CODON = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
AAS = sorted(set(CODON.values()) - {"*"})
CODONS_FOR = {aa: [c for c, a in CODON.items() if a == aa] for aa in set(CODON.values())}
THREE_TO_ONE = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G",
    "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
    "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
}
PMUT = re.compile(r"^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$")


def min_codon_dist(aa_a: str, aa_b: str) -> int:
    return min(sum(x != y for x, y in zip(ca, cb)) for ca in CODONS_FOR[aa_a] for cb in CODONS_FOR[aa_b])


def one_nt_adjacent(aa_a: str, aa_b: str) -> bool:
    return min_codon_dist(aa_a, aa_b) == 1


R_GRADES = ("1) Assoc w R", "2) Assoc w R - Interim")
# CHANGES-vs-ver1 values meaning the mutation was ALREADY R in v1 (the forecaster's seed).
V1_R_CHANGES = {"UP from AwRI to AwR", "DOWN from AwR to AwRI", "DOWN from AwRI to Uncertain",
                "DOWN from AwR to Uncertain", "SWITCH from AwRI to NotAwRI", "DOWN from AwRI to NotAwRI"}
# CHANGES values meaning newly-R in v2 (NOT R in v1) — the prediction target.
V2_ADDED_CHANGES = {"New AwR", "New AwRI", "UP from Uncertain to AwR", "UP from Uncertain to AwRI",
                    "UP from NotAwR to AwR", "UP from NotAwR to AwRI", "UP from NotAwRI to AwR",
                    "UP from NotAwRI to AwRI"}


def load_rows():
    rows = list(csv.reader(open(MASTER, encoding="utf-8", errors="replace"), delimiter="\t"))
    h = rows[0]
    ci = {n: next(i for i, c in enumerate(h) if c.strip() == n)
          for n in ("drug", "gene", "mutation", "FINAL CONFIDENCE GRADING", "CHANGES vs ver1")}
    out = []
    for r in rows[1:]:
        if max(ci.values()) >= len(r):
            continue
        m = PMUT.match(r[ci["mutation"]].strip())
        if not m:
            continue                       # protein single-substitution only (codon concept defined)
        wt3, pos, mut3 = m.group(1), int(m.group(2)), m.group(3)
        if wt3 not in THREE_TO_ONE or mut3 not in THREE_TO_ONE:
            continue
        out.append({
            "gene": r[ci["gene"]].strip(), "pos": pos,
            "wt": THREE_TO_ONE[wt3], "mut": THREE_TO_ONE[mut3],
            "final": r[ci["FINAL CONFIDENCE GRADING"]].strip(),
            "changes": r[ci["CHANGES vs ver1"]].strip(),
        })
    return out


def main() -> int:
    from math import comb

    rows = load_rows()
    for r in rows:
        r["is_R"] = r["final"] in R_GRADES
        r["was_R_v1"] = (r["is_R"] and r["changes"] in ("No change", "no change")) or r["changes"] in V1_R_CHANGES
        r["v2_added"] = r["is_R"] and (r["changes"] in V2_ADDED_CHANGES) and not r["was_R_v1"]

    # v1 resistance alleles per (gene,pos)
    v1R = {}
    for r in rows:
        if r["was_R_v1"]:
            v1R.setdefault((r["gene"], r["pos"]), {"wt": r["wt"], "alleles": set()})["alleles"].add(r["mut"])

    v2_added = [r for r in rows if r["v2_added"]]
    n_added_total = len(v2_added)
    added_at_v1R = [r for r in v2_added if (r["gene"], r["pos"]) in v1R]
    ceiling = len(added_at_v1R) / n_added_total if n_added_total else float("nan")

    # Position-matched Fisher: at v1-R positions, flagged=one-nt-adjacent to a v1-R allele.
    def flagged(gene, pos, mut):
        alleles = v1R[(gene, pos)]["alleles"]
        return mut not in alleles and any(one_nt_adjacent(mut, a) for a in alleles)

    added_keys = {(r["gene"], r["pos"], r["mut"]) for r in v2_added}
    A_flag = sum(1 for r in added_at_v1R if flagged(r["gene"], r["pos"], r["mut"]))
    A_tot = len(added_at_v1R)
    # Background pool B: every possible substitution at a v1-R position that is NOT WT, NOT a v1-R allele,
    # and NOT a v2-added-R mutation — the "random new substitution at a known resistance codon" null.
    B_flag = B_tot = 0
    for (gene, pos), info in v1R.items():
        wt, alleles = info["wt"], info["alleles"]
        for aa in AAS:
            if aa == wt or aa in alleles or (gene, pos, aa) in added_keys:
                continue
            B_tot += 1
            if any(one_nt_adjacent(aa, a) for a in alleles):
                B_flag += 1
    base_rate = B_flag / B_tot if B_tot else float("nan")
    test_rate = A_flag / A_tot if A_tot else float("nan")
    enrichment = (test_rate / base_rate) if (base_rate and A_tot) else float("nan")

    # Fisher exact (one-sided, enrichment) on [[A_flag, A_tot-A_flag],[B_flag, B_tot-B_flag]].
    def fisher_greater(a, b, c, d):
        n = a + b + c + d
        r1, c1 = a + b, a + c
        def hyp(k):
            return comb(r1, k) * comb(n - r1, c1 - k) / comb(n, c1)
        kmax = min(r1, c1)
        return sum(hyp(k) for k in range(a, kmax + 1))

    p_value = (fisher_greater(A_flag, A_tot - A_flag, B_flag, B_tot - B_flag)
               if (A_tot and B_tot) else float("nan"))

    powered = A_tot >= 10
    material = (enrichment == enrichment) and enrichment >= 2.0 and (p_value == p_value) and p_value < 0.05
    verdict = "PASS" if (powered and material) else ("UNDERPOWERED" if not powered else "FAIL")

    result = {
        "cell": "pre_resistance_predictive_backtest",
        "council_question": "Do v(N) one-nt-away flags predict v(N+1)-added DRMs above the base rate?",
        "substrate": "WHO M. tuberculosis catalogue v1(2021) -> v2(2023), protein substitutions, CHANGES-vs-ver1 diff",
        "scope": "protein single-substitution mutations only (codon concept defined); non-coding/rRNA excluded",
        "counts": {
            "protein_substitutions_parsed": len(rows),
            "v1_R_positions": len(v1R),
            "v2_added_R_total": n_added_total,
            "v2_added_R_at_v1R_position": A_tot,
        },
        "ceiling": {
            "fraction_v2_added_at_a_known_v1R_position": round(ceiling, 4) if ceiling == ceiling else None,
            "reading": ("the MAX a one-nt-from-known-codon forecaster could catch; the rest are new "
                        "positions/mechanisms unreachable by construction"),
        },
        "enrichment_test_position_matched_fisher": {
            "v2added_flagged": A_flag, "v2added_total_at_v1R": A_tot, "test_rate": round(test_rate, 4) if A_tot else None,
            "background_flagged": B_flag, "background_total": B_tot, "base_rate": round(base_rate, 4) if B_tot else None,
            "enrichment_ratio": round(enrichment, 3) if enrichment == enrichment else None,
            "fisher_p_one_sided": round(p_value, 5) if p_value == p_value else None,
            "powered_min10": powered, "material_ge2x_and_p05": material,
        },
        "verdict": verdict,
        "interpretation": (
            f"Of {n_added_total} v2-added resistance protein-substitutions, {A_tot} "
            f"({ceiling:.0%}) fall at a position that already had a v1 resistance allele (the forecaster's "
            f"structural ceiling). Among substitutions at those known resistance codons, v2-added ones are "
            f"one-nt-adjacent to the v1-R allele at rate {test_rate if A_tot else float('nan'):.3f} vs a "
            f"position-matched base rate {base_rate if B_tot else float('nan'):.3f} "
            f"(enrichment {enrichment:.2f}x, Fisher p={p_value:.4g}). Verdict {verdict}."
        ) if (A_tot and B_tot) else f"Insufficient v2-added protein substitutions at v1-R positions (A_tot={A_tot}).",
    }
    out = REPO / "wiki" / f"pre_resistance_predictive_backtest_{_date.today().isoformat()}.json"
    out.write_text(json.dumps(result, indent=2), encoding="utf-8")

    c = result["counts"]; e = result["enrichment_test_position_matched_fisher"]
    print(f"[backtest] parsed {c['protein_substitutions_parsed']} protein substitutions; "
          f"v1-R positions {c['v1_R_positions']}; v2-added-R {c['v2_added_R_total']}")
    print(f"[backtest] CEILING: {A_tot}/{n_added_total} v2-additions at a known v1-R position "
          f"({(ceiling*100 if ceiling==ceiling else 0):.0f}%) — rest unreachable by construction")
    print(f"[backtest] ENRICHMENT (position-matched): test {A_flag}/{A_tot}="
          f"{(test_rate if A_tot else 0):.3f} vs base {B_flag}/{B_tot}={(base_rate if B_tot else 0):.3f} "
          f"-> {e['enrichment_ratio']}x  Fisher p={e['fisher_p_one_sided']}")
    print(f"[backtest] powered(>=10)={powered}  material(>=2x,p<.05)={material}")
    print(f"[backtest] VERDICT: {verdict}")
    print(f"artifact -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
