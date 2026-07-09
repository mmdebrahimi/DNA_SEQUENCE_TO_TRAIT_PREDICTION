#!/usr/bin/env python
"""Does a protein language model add signal exactly where the curated DRM catalog is blind?

PRE-REGISTERED, 2026-07-09. The project's one free, independent, isolate-level, CONTINUOUS
wet-lab label is the Stanford PhenoSense fold-change (data/raw/hiv/*_DataSet.txt). The
deterministic NNRTI major-DRM catalog scores EFV at sensitivity 0.947 / specificity 0.904 at a
3x fold-change cutoff -- but it structurally CANNOT see 53 resistant isolates that carry no
catalog DRM (19 of them >=10x fold-change).

QUESTION: on the CATALOG-NEGATIVE subset (n=1111; 53 R vs 1058 S), does ESM2-650M zero-shot
masked-marginal scoring separate resistant from susceptible?

  PASS  = ESM AUROC >= 0.65 AND > mutation-burden baseline AND shuffled-null < 0.55
  FAIL  = otherwise. Report the number honestly whatever it is.

Baselines are DERIVED, not asserted (R2):
  - mutation burden (count of mutations) AUROC = 0.453  -- measured, BELOW chance, so a positive
    ESM result cannot be explained by resistant isolates simply carrying more mutations.
  - shuffled-label null.

HONEST PRIOR: ESM masked-marginals score EVOLUTIONARY likelihood, not drug binding. Resistance
mutations are often fitness-costly, so ESM should flag them as "damaging" -- but "damaging" is not
"resistance-conferring". Expect ESM to detect consequential mutations without distinguishing the
drug-resistant ones. A FAIL here is a real, publishable-internally negative that closes the
"should the decoder add a learned variant scorer?" question for this regime.

Scoring is LENGTH-NORMALIZED (mean over mutations, not sum) so the score cannot grow with
mutation count by construction.

CPU-only. The masked-marginal matrix is computed ONCE over the RT reference, then every isolate
scores by table lookup -- the whole run is minutes, not hours.
"""
import csv
import json
import random
import statistics as st
import sys

sys.path.insert(0, ".")
import dna_decode.data.hiv_amr as H

MODEL = "facebook/esm2_t33_650M_UR50D"   # the PEAK ESM2 checkpoint (650M > 3B > 15B on ProteinGym)
RT_REF = "data/hiv_ref/HIV1_RT_HXB2_cds.fna"
DATASET = "data/raw/hiv/NNRTI_DataSet.txt"
LOGP_CACHE = "data/processed/hiv_rt_esm650m_masked_marginals.json"
DRUG = "EFV"
CUTOFF = 3.0
BATCH = 8
AA = list("ACDEFGHIKLMNPQRSTVWY")
CODON = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "TAT":"Y","TAC":"Y","CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","AAT":"N","AAC":"N",
    "AAA":"K","AAG":"K","GAT":"D","GAC":"D","GAA":"E","GAG":"E","TGT":"C","TGC":"C",
    "TGG":"W","CGT":"R","CGC":"R","CGA":"R","CGG":"R","AGT":"S","AGC":"S","AGA":"R",
    "AGG":"R","GGT":"G","GGC":"G","GGA":"G","GGG":"G","TAA":"*","TAG":"*","TGA":"*",
}


def auroc(y, s):
    pr = sorted(zip(s, y))
    n = len(pr)
    rk = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j < n and pr[j][0] == pr[i][0]:
            j += 1
        for k in range(i, j):
            rk[k] = (i + 1 + j) / 2.0
        i = j
    P = sum(y)
    N = n - P
    if not P or not N:
        return float("nan")
    sr = sum(r for r, (_, l) in zip(rk, pr) if l == 1)
    return (sr - P * (P + 1) / 2) / (P * N)


def rt_protein():
    seq = "".join(l.strip() for l in open(RT_REF) if not l.startswith(">"))
    return "".join(CODON.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))


def load_rows():
    rows = list(csv.DictReader(open(DATASET, encoding="utf-8"), delimiter="\t"))
    pcols = [c for c in rows[0] if c.startswith("P") and c[1:].isdigit()]
    return rows, pcols


def isolate_muts(r, pcols):
    out = []
    for c in pcols:
        v = r[c]
        if v in ("-", ""):
            continue
        for aa in v:                        # mixtures are listed as several letters
            if aa.isalpha() and aa in AA:
                out.append((int(c[1:]), aa))
    return out


def main():
    prot = rt_protein()
    print(f"RT reference: {len(prot)} aa")

    # INTEGRITY GATE: HXB2 must match the catalog WT at every catalog DRM position.
    for pos, wt in H._RT_WT.items():
        assert prot[pos - 1] == wt, f"reference/catalog WT mismatch at {pos}: {prot[pos-1]} != {wt}"
    print(f"integrity gate OK: HXB2 translation == catalog WT at all {len(H._RT_WT)} DRM positions")

    rows, pcols = load_rows()

    # The dataset's '-' means "consensus B". Where HXB2 != consensus B, isolates carrying HXB2's
    # residue must LIST it as a mutation. Those positions have the wrong WT here -> exclude.
    #
    # Detect on PURE single-letter calls only. A mixture entry ("KR") contains the consensus
    # residue by construction, so a substring test flags nearly every position (it excluded 285/318
    # and silently emptied the score). Empirically the split is unambiguous: 314/318 positions have
    # ZERO pure HXB2 calls; the 4 that drift show up in 11-40% of isolates. Threshold 1%.
    drifted = set()
    for c in pcols:
        p = int(c[1:])
        if p > len(prot):
            continue
        pure = sum(1 for r in rows if (r[c] or "").strip() == prot[p - 1])
        if pure > 0.01 * len(rows):
            drifted.add(p)
    print(f"positions where HXB2 != consensus B (excluded): {len(drifted)} -> {sorted(drifted)}")
    assert len(drifted) < 20, f"drift detector over-fired ({len(drifted)} positions) — score would be empty"

    have = [r for r in rows if r[DRUG] not in ("NA", "", "-")]
    majors = H.NNRTI_RT_MAJOR_DRMS

    def catalog_call(r):
        return any(f"{H._RT_WT.get(p,'?')}{p}{aa}" in majors for p, aa in isolate_muts(r, pcols))

    sub = [r for r in have if not catalog_call(r)]
    y = [1 if float(r[DRUG]) >= CUTOFF else 0 for r in sub]
    print(f"\ncatalog-negative subset: n={len(sub)}  R={sum(y)}  S={len(y)-sum(y)}")

    # ---- masked-marginal matrix over the RT reference (computed ONCE, then cached) ----
    maxpos = min(max(int(c[1:]) for c in pcols), len(prot))
    positions = [p for p in range(1, maxpos + 1) if p not in drifted]
    aa_ids = {a: a for a in AA}                       # only membership is used downstream
    import os
    if os.path.exists(LOGP_CACHE):
        logp = {int(k): v for k, v in json.load(open(LOGP_CACHE)).items()}
        print(f"loaded cached masked-marginals for {len(logp)} positions ({LOGP_CACHE})")
    else:
        import torch
        from transformers import AutoModelForMaskedLM, AutoTokenizer
        tok = AutoTokenizer.from_pretrained(MODEL)
        model = AutoModelForMaskedLM.from_pretrained(MODEL).eval()
        print(f"model loaded: {MODEL} (cpu)")
        tok_ids = {a: tok.convert_tokens_to_ids(a) for a in AA}
        ids = tok(prot, return_tensors="pt")["input_ids"]
        logp = {}
        with torch.no_grad():
            for s in range(0, len(positions), BATCH):
                chunk = positions[s:s + BATCH]
                batch = ids.repeat(len(chunk), 1)
                for r, p in enumerate(chunk):
                    batch[r, p] = tok.mask_token_id   # token idx == pos (CLS at 0)
                out = model(batch).logits.float().log_softmax(-1)
                for r, p in enumerate(chunk):
                    logp[p] = {a: out[r, p, i].item() for a, i in tok_ids.items()}
                if (s // BATCH) % 10 == 0:
                    print(f"  masked {s + len(chunk)}/{len(positions)} positions", flush=True)
        json.dump(logp, open(LOGP_CACHE, "w"))
        print(f"cached masked-marginals -> {LOGP_CACHE}")

    # ---- score isolates: LENGTH-NORMALIZED mean damage ----
    def esm_score(r):
        d = []
        for p, aa in isolate_muts(r, pcols):
            if p in logp and p <= len(prot):
                wt = prot[p - 1]
                if wt in aa_ids and aa in aa_ids and wt != aa:
                    d.append(logp[p][wt] - logp[p][aa])    # higher = more damaging
        return st.mean(d) if d else 0.0

    def n_scored(r):
        return sum(1 for p, aa in isolate_muts(r, pcols)
                   if p in logp and p <= len(prot) and prot[p - 1] != aa and aa in aa_ids)

    cov = [n_scored(r) for r in sub]
    zero = sum(1 for c in cov if c == 0)
    print(f"\ncoverage: mutations actually scored per isolate — median {st.median(cov)}, "
          f"mean {st.mean(cov):.1f}, isolates with ZERO scored mutations: {zero}/{len(sub)}")
    # A near-empty score would produce a meaningless AUROC ~0.5. Refuse to report one.
    assert st.median(cov) >= 2, f"degenerate score: median scored mutations = {st.median(cov)}"
    assert zero < 0.25 * len(sub), f"{zero} isolates have no scored mutation"

    esm = [esm_score(r) for r in sub]
    burden = [len(isolate_muts(r, pcols)) for r in sub]
    esm_sum = [esm_score(r) * max(1, len(isolate_muts(r, pcols))) for r in sub]

    a_esm = auroc(y, esm)
    a_burden = auroc(y, burden)
    a_sum = auroc(y, esm_sum)
    rng = random.Random(0)
    nulls = []
    for _ in range(200):
        ys = y[:]
        rng.shuffle(ys)
        nulls.append(auroc(ys, esm))
    a_null = st.median(nulls)

    # ---- POSITIVE CONTROLS: is the instrument sensitive at all? ----
    # Without these a null on the catalog-negative subset is UNINTERPRETABLE: it cannot
    # distinguish "no signal beyond the catalog" from "ESM cannot see drug resistance anywhere".

    # Control A: FULL cohort (catalog DRMs included). ESM should do well here if it sees resistance.
    y_full = [1 if float(r[DRUG]) >= CUTOFF else 0 for r in have]
    esm_full = [esm_score(r) for r in have]
    cat_full = [1 if catalog_call(r) else 0 for r in have]
    a_esm_full = auroc(y_full, esm_full)
    a_cat_full = auroc(y_full, cat_full)

    # Control B: do the 16 catalog major DRMs rank as damaging among all 19 substitutions
    # at their own position? Reports the median percentile of the true DRM's damage score.
    pct = []
    for drm in sorted(majors):
        wt, mut = drm[0], drm[-1]
        pos = int(drm[1:-1])
        if pos not in logp:
            continue
        dmg = {a: logp[pos][wt] - logp[pos][a] for a in AA if a != wt}
        rank = sorted(dmg, key=lambda a: dmg[a])           # ascending damage
        pct.append((drm, (rank.index(mut) + 1) / len(rank)))
    med_pct = st.median([p for _, p in pct]) if pct else float("nan")

    print("\n--- POSITIVE CONTROLS ---")
    print(f"A. FULL cohort (n={len(have)}, R={sum(y_full)}): "
          f"ESM AUROC={a_esm_full:.3f}   catalog AUROC={a_cat_full:.3f}")
    print(f"B. catalog DRMs' damage percentile at own position (median over {len(pct)}): {med_pct:.2f}")
    for drm, p in sorted(pct, key=lambda x: -x[1])[:5]:
        print(f"     {drm:7s} percentile {p:.2f}")
    instrument_ok = a_esm_full >= 0.60 or med_pct >= 0.70
    print(f"   instrument sensitive to resistance? {'YES' if instrument_ok else 'NO'}")
    if not instrument_ok:
        print("   => a null on the catalog-negative subset means the LM does not see drug")
        print("      resistance AT ALL (evolutionary likelihood != drug binding), NOT that the")
        print("      catalog's blind spot is irreducible. Report it as such.")

    print("\n" + "=" * 62)
    print(f"CATALOG-NEGATIVE SUBSET  n={len(sub)}  R={sum(y)}  S={len(y)-sum(y)}  ({DRUG} >= {CUTOFF}x)")
    print(f"  ESM2-650M mean-damage AUROC      = {a_esm:.3f}   <- the headline")
    print(f"  ESM2-650M sum-damage  AUROC      = {a_sum:.3f}   (burden-confounded variant)")
    print(f"  mutation-burden baseline AUROC   = {a_burden:.3f}")
    print(f"  shuffled-label null (median/200) = {a_null:.3f}")
    passed = a_esm >= 0.65 and a_esm > a_burden and a_null < 0.55
    print(f"\n  PASS bar: ESM >= 0.65 AND > burden AND null < 0.55")
    print(f"  VERDICT: {'PASS - the LM adds signal where the catalog is blind' if passed else 'FAIL - no usable signal beyond the catalog'}")
    print("=" * 62)

    json.dump({
        "date": "2026-07-09", "model": MODEL, "drug": DRUG, "cutoff_fold": CUTOFF,
        "subset": "catalog-negative (no NNRTI major DRM)",
        "n": len(sub), "n_resistant": sum(y), "n_susceptible": len(y) - sum(y),
        "esm_mean_damage_auroc": round(a_esm, 4),
        "esm_sum_damage_auroc": round(a_sum, 4),
        "mutation_burden_auroc": round(a_burden, 4),
        "shuffled_null_auroc": round(a_null, 4),
        "pass_bar": "esm>=0.65 AND esm>burden AND null<0.55",
        "passed": bool(passed),
        "control_a_full_cohort_esm_auroc": round(a_esm_full, 4),
        "control_a_full_cohort_catalog_auroc": round(a_cat_full, 4),
        "control_a_full_cohort_n": len(have), "control_a_full_cohort_n_resistant": sum(y_full),
        "control_b_drm_damage_percentile_median": round(med_pct, 4),
        "instrument_sensitive_to_resistance": bool(instrument_ok),
        "excluded_positions_hxb2_ne_consensus": sorted(drifted),
        "honest_prior": "ESM scores evolutionary likelihood, not drug binding; resistance mutations "
                        "are often fitness-costly, so 'damaging' != 'resistance-conferring'.",
    }, open("wiki/hiv_esm_vs_catalog_2026-07-09.json", "w"), indent=1)
    print("wrote wiki/hiv_esm_vs_catalog_2026-07-09.json")


if __name__ == "__main__":
    main()
