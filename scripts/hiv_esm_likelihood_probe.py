#!/usr/bin/env python
"""WHY is ESM blind to HIV drug resistance? Decompose the effect, and test the likelihood story.

`hiv_esm_vs_catalog.py` established the negative and Control B showed the 16 NNRTI DRMs sit at
median damage percentile 0.29 at their own positions (i.e. ESM thinks they are among the *less*
damaging substitutions there). `hiv_esm_mechanism_probe.py` then killed the "fitness cost" story.

The memo/CLAUDE.md now assert a SECOND causal claim: "masked-marginals score evolutionary
likelihood; DRMs are viable, circulating, drug-selected variants abundant in UniRef, so the LM
calls the resistant variant likely, not damaging." That claim is still UNTESTED, and a rival
explanation predicts the same 0.29:

  RIVAL (positional tolerance): DRM positions are simply variable/unconserved. Any substitution
  there looks benign, and nothing about the DRM residue itself is special.

They are distinguishable, and the distinction changes the design rule:
  - positional tolerance  -> "don't trust an LM at tolerant sites"  (a site-level caveat)
  - mutant-specific likelihood -> "the LM has SEEN the resistant variant" (a label/leakage-shaped
    caveat that no amount of model scale or site-filtering fixes)

Three probes, all on the cached masked-marginal matrix. No new inference, no network.

PROBE A (positional tolerance). Shannon entropy of ESM's 20-AA distribution at DRM positions vs
all RT positions. If DRM sites are unusually tolerant, their entropy percentile >> 0.5.

PROBE B (mutant specificity). Rank of the DRM residue among the 19 non-WT alternatives, by ESM
log-prob (rank 1 = ESM's single most likely substitution). Null median rank = 10.

PROBE C (the likelihood story, direct). Across every (position, substitution) ESM scores, does the
model's log-prob track the substitution's EMPIRICAL FREQUENCY in the 2272-isolate cohort?

PROBE D (the BLOSUM control -- the one that matters). ESM favours chemically conservative
substitutions, and conservative substitutions are common in ANY protein. So probes B and C might
measure generic amino-acid exchangeability rather than anything ESM learned about HIV. Re-run both
against BLOSUM62:
  - B: BLOSUM62 ranks the same DRMs at median 4.0/19 vs ESM's 4.5, and ESM is more favourable on
    only 12/38. **Probe B is fully explained by exchangeability.** ESM is NOT specifically
    up-weighting the resistant residue; the DRMs are just chemically mild substitutions.
  - C: partial rho(ESM, count | BLOSUM) = +0.218 (vs partial rho(BLOSUM, count | ESM) = +0.161), so
    ESM does retain an independent, modest correlation with circulating frequency.

CONCLUSION (what survives all four): DRM sites are ordinarily conserved (A), and resistance is
achieved through chemically CONSERVATIVE substitutions there (B+D). Any exchangeability- or
likelihood-based scorer -- BLOSUM62 and ESM alike -- therefore rates the resistant residue benign.
That is a structural blindness, not a model-capacity problem, and it predicts the same failure for
any conservation-based predictor on any antagonistically-selected phenotype.
"""
import csv
import json
import math
import statistics as st
import sys
from collections import Counter

sys.path.insert(0, ".")
import dna_decode.data.hiv_amr as H
from scripts.hiv_esm_vs_catalog import AA, LOGP_CACHE, isolate_muts, load_rows, rt_protein

DRIFT = {122, 214, 272, 277}
NRTI_DRMS = ["K65R", "K70E", "K70G", "K70N", "K70R", "K70T", "L74V", "M184I", "M184V", "M41L",
             "Q151M", "L210W", "T215E", "T215F", "T215I", "T215L", "T215V", "T215Y",
             "V75A", "V75I", "V75M", "V75T"]


def spearman(xs, ys):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(v):
            j = i
            while j < len(v) and v[order[j]] == v[order[i]]:
                j += 1
            for k in range(i, j):
                r[order[k]] = (i + j - 1) / 2.0
            i = j
        return r
    rx, ry = rank(xs), rank(ys)
    n = len(xs)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = sum((a - mx) ** 2 for a in rx) ** 0.5
    dy = sum((b - my) ** 2 for b in ry) ** 0.5
    return num / (dx * dy) if dx and dy else float("nan")


def main():
    prot = rt_protein()
    logp = {int(k): v for k, v in json.load(open(LOGP_CACHE)).items()}
    rows, pcols = load_rows()

    def entropy(pos):
        lp = [logp[pos][a] for a in AA]
        m = max(lp)
        w = [math.exp(x - m) for x in lp]
        z = sum(w)
        p = [x / z for x in w]
        return -sum(x * math.log(x) for x in p if x > 0)

    positions = sorted(logp)
    ent = {p: entropy(p) for p in positions}
    all_ent = sorted(ent.values())

    def ent_pct(p):
        return sum(1 for e in all_ent if e < ent[p]) / len(all_ent)

    drm_pos = sorted({int(d[1:-1]) for d in H.NNRTI_RT_MAJOR_DRMS} |
                     {int(d[1:-1]) for d in NRTI_DRMS})
    drm_pos = [p for p in drm_pos if p in logp]

    print("PROBE A -- positional tolerance (ESM entropy at DRM sites vs all RT sites)")
    print(f"  all {len(positions)} positions: median entropy {st.median(all_ent):.3f}")
    dp = [ent[p] for p in drm_pos]
    pcts = [ent_pct(p) for p in drm_pos]
    print(f"  {len(drm_pos)} DRM positions:  median entropy {st.median(dp):.3f}  "
          f"median entropy-percentile {st.median(pcts):.3f}")
    tolerant = st.median(pcts) > 0.65
    print(f"  DRM sites unusually tolerant? {'YES' if tolerant else 'NO'} "
          f"(null percentile 0.5)\n")

    def mut_rank(drm):
        wt, mut, pos = drm[0], drm[-1], int(drm[1:-1])
        if pos not in logp or prot[pos - 1] != wt:
            return None
        alts = [a for a in AA if a != wt]
        order = sorted(alts, key=lambda a: -logp[pos][a])   # most likely first
        return order.index(mut) + 1

    print("PROBE B -- mutant specificity (rank of the DRM residue among 19 alternatives)")
    for name, lst in (("NNRTI", sorted(H.NNRTI_RT_MAJOR_DRMS)), ("NRTI", NRTI_DRMS)):
        rk = [r for r in (mut_rank(d) for d in lst) if r]
        print(f"  {name}: n={len(rk)}  median rank {st.median(rk):.1f} / 19   (null median 10)")
    allrk = [r for r in (mut_rank(d) for d in sorted(H.NNRTI_RT_MAJOR_DRMS)) if r] + \
            [r for r in (mut_rank(d) for d in NRTI_DRMS) if r]
    med_rank = st.median(allrk)
    top5 = sum(1 for r in allrk if r <= 5) / len(allrk)
    print(f"  ALL DRMs: median rank {med_rank:.1f}/19, {top5:.0%} are in ESM's top-5 most likely")
    specific = med_rank < 8
    print(f"  ESM specifically favours the resistant residue? {'YES' if specific else 'NO'}\n")

    # PROBE C: does ESM log-prob track empirical frequency of the substitution in the cohort?
    freq = Counter()
    for r in rows:
        for p, aa in isolate_muts(r, pcols):
            if p in logp and p <= len(prot) and p not in DRIFT and prot[p - 1] != aa:
                freq[(p, aa)] += 1
    obs = [(p, aa) for (p, aa) in freq if freq[(p, aa)] >= 5]
    x = [logp[p][aa] for p, aa in obs]
    y = [freq[(p, aa)] for p, aa in obs]
    rho_obs = spearman(x, y)

    # Contrast: over ALL possible substitutions, does ESM rank the OBSERVED ones above the unseen?
    seen, unseen = [], []
    for p in positions:
        if p in DRIFT or p > len(prot):
            continue
        for aa in AA:
            if aa == prot[p - 1]:
                continue
            (seen if freq.get((p, aa), 0) >= 5 else unseen).append(logp[p][aa])
    print("PROBE C -- does ESM's log-prob track what is OBSERVED circulating?")
    print(f"  substitutions seen >=5x in the cohort: {len(seen)};  rarely/never seen: {len(unseen)}")
    print(f"  mean ESM log-prob  seen {st.mean(seen):+.3f}   unseen {st.mean(unseen):+.3f}   "
          f"delta {st.mean(seen)-st.mean(unseen):+.3f}")
    print(f"  Spearman(ESM log-prob, empirical count) over {len(obs)} observed substitutions "
          f"= {rho_obs:+.3f}")
    tracks = st.mean(seen) > st.mean(unseen) and rho_obs > 0.2
    print(f"  ESM reproduces circulating variation? {'YES' if tracks else 'NO'}")

    # ---- PROBE D: the BLOSUM control. Is any of this ESM-specific, or just exchangeability? ----
    from Bio.Align import substitution_matrices as _sm
    B = _sm.load("BLOSUM62")

    def bl_rank(drm):
        wt, mut = drm[0], drm[-1]
        alts = [a for a in AA if a != wt]
        return sorted(alts, key=lambda a: -B[wt, a]).index(mut) + 1

    drms = sorted(H.NNRTI_RT_MAJOR_DRMS) + NRTI_DRMS
    pairs = [(mut_rank(d), bl_rank(d)) for d in drms if mut_rank(d)]
    med_bl = st.median([b for _, b in pairs])
    esm_better = sum(1 for e, b in pairs if e < b)

    bl_obs = [B[prot[p - 1], aa] for p, aa in obs]
    r_ec, r_eb, r_bc = rho_obs, spearman(x, bl_obs), spearman(bl_obs, y)

    def partial(rxy, rxz, ryz):
        d = ((1 - rxz ** 2) * (1 - ryz ** 2)) ** 0.5
        return (rxy - rxz * ryz) / d if d else float("nan")
    p_esm = partial(r_ec, r_eb, r_bc)
    p_bl = partial(r_bc, r_eb, r_ec)

    print("\nPROBE D -- BLOSUM62 control (is ANY of this ESM-specific?)")
    print(f"  DRM rank: ESM median {med_rank:.1f}/19   BLOSUM median {med_bl:.1f}/19   "
          f"ESM better on {esm_better}/{len(pairs)}")
    print(f"  rho(BLOSUM, empirical count) = {r_bc:+.3f}  (ESM was {r_ec:+.3f})")
    print(f"  PARTIAL rho(ESM, count | BLOSUM) = {p_esm:+.3f}   "
          f"PARTIAL rho(BLOSUM, count | ESM) = {p_bl:+.3f}")
    b_is_esm_specific = med_rank < med_bl - 0.5 and esm_better > len(pairs) * 0.6
    c_survives = p_esm > 0.15
    print(f"  Probe B ESM-specific? {'YES' if b_is_esm_specific else 'NO -- explained by exchangeability'}")
    print(f"  Probe C survives the control? {'YES' if c_survives else 'NO'}")

    print("\n" + "=" * 74)
    if tolerant and not specific:
        verdict = ("POSITIONAL-TOLERANCE story: DRM sites are unconserved. "
                   "The likelihood claim must be corrected.")
    elif specific and not b_is_esm_specific:
        verdict = ("CONSERVATIVE-SUBSTITUTION story. DRM sites are ordinarily conserved (A), and "
                   "the resistant residues are chemically MILD substitutions there -- BLOSUM62 "
                   f"ranks them {med_bl:.1f}/19, no worse than ESM's {med_rank:.1f}/19 (D). So the "
                   "apparent 'ESM up-weights the resistant residue' is generic exchangeability, "
                   "NOT evidence the model memorized circulating DRMs. ESM does retain a modest "
                   f"independent tie to circulating frequency (partial rho {p_esm:+.3f}), but that "
                   "is not what drives the DRM ranking. Any exchangeability- or likelihood-based "
                   "scorer is structurally blind here -- a property of the phenotype, not of model "
                   "capacity.")
    elif specific and b_is_esm_specific and tracks:
        verdict = ("MEMORIZATION story SUPPORTED: ESM favours the resistant residue beyond what "
                   "BLOSUM explains, and tracks circulating variation.")
    else:
        verdict = "MIXED -- report the components, claim nothing more."
    print(verdict)
    print("=" * 74)

    json.dump({"date": "2026-07-09", "model": "facebook/esm2_t33_650M_UR50D",
               "probeA_all_positions_median_entropy": round(st.median(all_ent), 4),
               "probeA_drm_positions_median_entropy": round(st.median(dp), 4),
               "probeA_drm_median_entropy_percentile": round(st.median(pcts), 4),
               "probeA_drm_sites_unusually_tolerant": bool(tolerant),
               "probeB_drm_median_rank_of_19": med_rank,
               "probeB_frac_in_top5": round(top5, 4),
               "probeB_esm_favours_resistant_residue": bool(specific),
               "probeC_mean_logp_seen": round(st.mean(seen), 4),
               "probeC_mean_logp_unseen": round(st.mean(unseen), 4),
               "probeC_spearman_logp_vs_empirical_count": round(rho_obs, 4),
               "probeC_esm_reproduces_circulating_variation": bool(tracks),
               "probeD_blosum_drm_median_rank_of_19": med_bl,
               "probeD_esm_drm_median_rank_of_19": med_rank,
               "probeD_esm_better_than_blosum_on": f"{esm_better}/{len(pairs)}",
               "probeD_rho_blosum_vs_empirical_count": round(r_bc, 4),
               "probeD_partial_rho_esm_given_blosum": round(p_esm, 4),
               "probeD_partial_rho_blosum_given_esm": round(p_bl, 4),
               "probeD_probeB_is_esm_specific": bool(b_is_esm_specific),
               "probeD_probeC_survives_control": bool(c_survives),
               "verdict": verdict},
              open("wiki/hiv_esm_likelihood_probe_2026-07-09.json", "w"), indent=1)
    print("\nwrote wiki/hiv_esm_likelihood_probe_2026-07-09.json")


if __name__ == "__main__":
    main()
