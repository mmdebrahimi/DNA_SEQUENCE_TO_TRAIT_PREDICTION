#!/usr/bin/env python
"""GO/NO-GO pre-check for the ΔΔG bet: are the EFV catalog-blind-spot resistant isolates resistant
via NNRTI-POCKET (binding-site) mutations, or via non-pocket mechanisms a binding-ΔΔG scorer can't see?

Zero-tool, D-free. No FoldX/Rosetta, no structures, no embeddings — pure position analysis over
data/raw/hiv/NNRTI_DataSet.txt (Stanford HIVDB genotype+PhenoSense fold-change).

WHY the naive test fails: the blind-spot subset is DEFINED as isolates with no major DRM at the 8
primary NNRTI-pocket positions {100,101,103,106,181,188,190,230}. The user's pocket windows
(100-110 / 179-190 / 225-236) cover 35/318 = 11% of RT, so with ~12 mutations/isolate, "≥1 window
mutation" is TRUE by chance for almost every isolate — uninformative. The honest question needs BOTH:
  (a) KNOWN-FUNCTIONAL secondary-pocket NNRTI mutations (curated), not just any window position, and
  (b) an R-vs-S ENRICHMENT null: are the blind-spot R enriched for these vs blind-spot S? If R and S
      carry pocket mutations at the same rate, the pocket is NOT what makes the R isolates resistant.

VERDICT LOGIC:
  GO   = blind-spot R are ENRICHED for known-functional secondary-pocket NNRTI mutations vs S
         (uncatalogued binding-site DRMs exist → ddG has real targets on the blind spot).
  NO-GO= no enrichment (blind-spot R resistance is NOT pocket-mediated → binding-ΔΔG is blind to it,
         so the HIV blind spot is the WRONG test set for a competitive-binding physics scorer).
"""
import csv
import json
import statistics as st
import sys

sys.path.insert(0, ".")
import dna_decode.data.hiv_amr as H

DATASET = "data/raw/hiv/NNRTI_DataSet.txt"
DRUG = "EFV"
CUTOFF = 3.0
DRIFT = {122, 214, 272, 277}  # HXB2 != consensus B (from hiv_esm_vs_catalog.py)
AA = set("ACDEFGHIKLMNPQRSTVWY")

# Curated KNOWN-FUNCTIONAL secondary/accessory NNRTI-pocket mutations (Stanford HIVDB), WITHIN the
# user's windows (100-110 / 179-190 / 225-236) but NOT in the 8 major-DRM positions.
KNOWN_SECONDARY_IN_WINDOW = {
    108: set("I"),                 # V108I
    179: set("DEFTL"),             # V179D/E/F/T/L
    225: set("H"),                 # P225H
    227: set("CL"),                # F227C/L
    234: set("I"),                 # L234I
    236: set("L"),                 # P236L
}
# Functional NNRTI-associated but OUTSIDE the windows (reported as an extended set + caveat).
KNOWN_FUNCTIONAL_OUT_OF_WINDOW = {
    98: set("G"),                  # A98G
    138: set("KAGQR"),             # E138K/A/G/Q/R (ETR/RPV)
    221: set("Y"),                 # H221Y
    238: set("TN"),                # K238T/N
    318: set("F"),                 # Y318F
}
# Consensus-B wild-type at the curated positions (for readable mutation labels; _RT_WT covers only
# the 8 major positions, so these are supplied here).
CURATED_WT = {98: "A", 108: "V", 138: "E", 179: "V", 221: "H", 225: "P",
              227: "F", 234: "L", 236: "P", 238: "K", 318: "Y"}


def load():
    rows = list(csv.DictReader(open(DATASET, encoding="utf-8"), delimiter="\t"))
    pcols = [c for c in rows[0] if c.startswith("P") and c[1:].isdigit()]
    return rows, pcols


def isolate_muts(r, pcols):
    out = []
    for c in pcols:
        v = r[c]
        if v in ("-", ""):
            continue
        for aa in v:
            if aa in AA:
                out.append((int(c[1:]), aa))
    return out


def main():
    rows, pcols = load()
    majors = H.NNRTI_RT_MAJOR_DRMS
    win_positions = set(range(100, 111)) | set(range(179, 191)) | set(range(225, 237))
    maxpos = max(int(c[1:]) for c in pcols)

    def catalog_negative(r):
        return not any(f"{H._RT_WT.get(p,'?')}{p}{aa}" in majors
                       for p, aa in isolate_muts(r, pcols))

    have = [r for r in rows if r[DRUG] not in ("NA", "", "-")]
    blind = [r for r in have if catalog_negative(r)]
    R = [r for r in blind if float(r[DRUG]) >= CUTOFF]
    S = [r for r in blind if float(r[DRUG]) < CUTOFF]

    def has_known_secondary(r):  # curated in-window functional
        return any(p in KNOWN_SECONDARY_IN_WINDOW and aa in KNOWN_SECONDARY_IN_WINDOW[p]
                   for p, aa in isolate_muts(r, pcols) if p not in DRIFT)

    def has_known_functional_ext(r):  # in-window + out-of-window functional
        for p, aa in isolate_muts(r, pcols):
            if p in DRIFT:
                continue
            if p in KNOWN_SECONDARY_IN_WINDOW and aa in KNOWN_SECONDARY_IN_WINDOW[p]:
                return True
            if p in KNOWN_FUNCTIONAL_OUT_OF_WINDOW and aa in KNOWN_FUNCTIONAL_OUT_OF_WINDOW[p]:
                return True
        return False

    def has_any_window(r):  # the naive/uninformative metric, for contrast
        return any(p in win_positions and p not in DRIFT for p, aa in isolate_muts(r, pcols))

    def frac(rs, f):
        return sum(1 for r in rs if f(r)) / len(rs) if rs else float("nan")

    def n(rs, f):
        return sum(1 for r in rs if f(r))

    # Which specific known-functional secondary mutations appear in R (with counts + R-vs-S)?
    from collections import Counter
    def known_muts(rs):
        c = Counter()
        for r in rs:
            for p, aa in isolate_muts(r, pcols):
                if p in DRIFT:
                    continue
                if (p in KNOWN_SECONDARY_IN_WINDOW and aa in KNOWN_SECONDARY_IN_WINDOW[p]) or \
                   (p in KNOWN_FUNCTIONAL_OUT_OF_WINDOW and aa in KNOWN_FUNCTIONAL_OUT_OF_WINDOW[p]):
                    wt = H._RT_WT.get(p) or CURATED_WT.get(p, "?")
                    c[f"{wt}{p}{aa}"] += 1
        return c
    Rk, Sk = known_muts(R), known_muts(S)

    r_sec, s_sec = frac(R, has_known_secondary), frac(S, has_known_secondary)
    r_ext, s_ext = frac(R, has_known_functional_ext), frac(S, has_known_functional_ext)
    r_win, s_win = frac(R, has_any_window), frac(S, has_any_window)
    enrich = (r_ext / s_ext) if s_ext and s_ext == s_ext else float("nan")
    # Burden-adjusted: the naive any-window R/S ratio captures the mutation-burden effect (R isolates
    # carry more mutations). Dividing the functional enrichment by it isolates the FUNCTIONAL-specific
    # signal beyond burden.
    naive_ratio = (r_win / s_win) if s_win else float("nan")
    burden_adj = (enrich / naive_ratio) if naive_ratio and naive_ratio == naive_ratio else float("nan")

    # burden context
    r_burden = st.median([len(isolate_muts(r, pcols)) for r in R]) if R else float("nan")

    print(f"EFV catalog-blind-spot subset: R={len(R)} (fold>={CUTOFF}x)  S={len(S)}  "
          f"(median #mut/R-isolate {r_burden})\n")
    print(f"{'metric':52s} {'R':>10s} {'S':>10s}")
    print("-" * 74)
    print(f"{'>=1 known secondary-pocket NNRTI mut (in-window curated)':52s} "
          f"{r_sec:>9.1%} {s_sec:>9.1%}")
    print(f"{'>=1 known FUNCTIONAL NNRTI mut (+138/221/238/318)':52s} "
          f"{r_ext:>9.1%} {s_ext:>9.1%}")
    print(f"{'>=1 ANY window-position mut (NAIVE, ~11% chance/mut)':52s} "
          f"{r_win:>9.1%} {s_win:>9.1%}")
    print(f"\nR-vs-S enrichment (known-functional): {enrich:.2f}x  "
          f"(R {n(R, has_known_functional_ext)}/{len(R)} vs S {n(S, has_known_functional_ext)}/{len(S)})")
    print(f"naive any-window R/S ratio (BURDEN control): {naive_ratio:.2f}x  ->  "
          f"burden-ADJUSTED functional enrichment = {burden_adj:.2f}x")
    print(f"\nKnown-functional NNRTI mutations in blind-spot R (count R | count S):")
    for m in sorted(set(Rk) | set(Sk), key=lambda x: -Rk[x]):
        flag = "in-window" if any(m[1:-1] and int("".join(filter(str.isdigit, m))) in KNOWN_SECONDARY_IN_WINDOW for _ in [0]) else "ext"
        print(f"  {m:10s} R={Rk[m]:3d}  S={Sk[m]:3d}  [{flag}]")

    # VERDICT: enrichment >= 2x AND >=15% of R carry a known-functional mutation => GO; else NO-GO.
    go = (enrich >= 2.0) and (r_ext >= 0.15) and (n(R, has_known_functional_ext) >= 5)
    verdict = "GO" if go else "NO-GO"
    print("\n" + "=" * 74)
    if go:
        print(f"VERDICT: GO — blind-spot R are {enrich:.1f}x enriched for known-functional NNRTI-pocket")
        print("  mutations vs S, and {:.0%} carry one. Uncatalogued binding-site DRMs EXIST on the".format(r_ext))
        print("  blind spot → a binding-ΔΔG scorer has real targets. A FoldX/Rosetta pilot is warranted.")
    else:
        print("VERDICT: NO-GO — blind-spot R are NOT meaningfully enriched for known-functional")
        print(f"  NNRTI-pocket mutations vs S (enrichment {enrich:.2f}x, only {r_ext:.0%} of R carry one).")
        print("  Their resistance is NOT pocket-mediated → binding-ΔΔG is structurally blind to it.")
        print("  The HIV EFV blind spot is the WRONG test set for a competitive-binding physics scorer;")
        print("  do NOT install FoldX/Rosetta for this test.")
    print("=" * 74)

    json.dump({
        "date": "2026-07-09", "drug": DRUG, "cutoff_fold": CUTOFF,
        "subset": "EFV catalog-blind-spot (no major DRM at 8 pocket positions)",
        "n_R": len(R), "n_S": len(S), "median_muts_per_R": r_burden,
        "windows": "100-110 / 179-190 / 225-236 (35/318 = 11% of RT)",
        "frac_R_known_secondary_in_window": round(r_sec, 4),
        "frac_S_known_secondary_in_window": round(s_sec, 4),
        "frac_R_known_functional_ext": round(r_ext, 4),
        "frac_S_known_functional_ext": round(s_ext, 4),
        "frac_R_any_window_naive": round(r_win, 4),
        "frac_S_any_window_naive": round(s_win, 4),
        "enrichment_R_over_S_known_functional": round(enrich, 3),
        "naive_window_R_over_S_ratio_burden_control": round(naive_ratio, 3),
        "burden_adjusted_functional_enrichment": round(burden_adj, 3),
        "n_R_with_known_functional": n(R, has_known_functional_ext),
        "known_functional_muts_in_R": dict(Rk), "known_functional_muts_in_S": dict(Sk),
        "verdict": verdict,
        "verdict_rule": "GO iff enrichment>=2.0 AND frac_R>=0.15 AND n_R>=5; else NO-GO",
        "caveat": "window-curated + extended (138/221/238/318) known-functional NNRTI mutations only; "
                  "a novel uncatalogued pocket mutation not in either set would be missed (conservative "
                  "toward GO-undercount). Naive any-window metric shown only to demonstrate it is "
                  "chance-dominated and must not be used.",
    }, open("wiki/hiv_blindspot_pocket_localization_2026-07-09.json", "w"), indent=1)
    print("\nwrote wiki/hiv_blindspot_pocket_localization_2026-07-09.json")


if __name__ == "__main__":
    main()
