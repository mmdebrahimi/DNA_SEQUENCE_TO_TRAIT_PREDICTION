"""Can a model predict the effect of a mutation token it has NEVER SEEN? The fair test of an alphabet.

WHY THIS AND NOT THE EARLIER TEST. The 2026-09-29 CTX-M-14 run tested a context-free substitution table on
one protein and it lost to BLOSUM62. That was a strawman for a Genome LM: a tokenizer is never meant to
carry meaning alone. And the deeper objection to a GLM here is not that test, it is HEADROOM -- on HIV the
additive bag-of-tokens model already reaches Spearman 0.911 / R2 0.844 (n=2,168, EFV) and
`scripts/hiv_epistasis.py` showed pairwise interactions add nothing (2/24 drug-cells CI-positive,
`wiki/hiv_epistasis_result_2026-07-11.md`). A transformer cannot add much where a linear model explains 88%
of the variance.

BUT THE HEADROOM ARGUMENT HAS ONE PRECISE BLIND SPOT, AND IT IS EXACTLY THE ALPHABET IDEA. Every model in
that comparison is keyed on token IDENTITY (one-hot over `<pos><AA>`). A mutation absent from training is an
all-zero row and scores EXACTLY ZERO by construction -- it cannot be wrong, it simply cannot answer. That is
why the catalog's only measured failure mode is COMPLETENESS, and why both confirmed gaps (rmtE1, V179F)
were found by accident when new labels arrived.

A model keyed on token FEATURES -- position, chemistry, pretrained conservation -- still has a feature row
for a token it has never seen, so it CAN answer. Whether that answer is any good is the open question, and
it is the one thing the 0.911 ceiling does not settle.

THE TEST. Two stages, both on the free independent Stanford HIVDB PhenoSense labels:
  stage 1  fit the additive identity model (the incumbent) -> a fitted effect coefficient per mutation
           token. These coefficients are the target: "how much does this mutation move log10 fold-change".
  stage 2  LEAVE-ONE-TOKEN-OUT: predict a held-out token's coefficient from its FEATURES alone, having
           never seen that token's data. The identity model structurally cannot do this; that is the point.

PRE-REGISTERED BARS (frozen before any number was seen):
  B1  the feature model must beat predicting the MEAN coefficient (the best an identity-keyed model can do
      on an unseen token). If this fails, alphabets carry no generalisable signal at all and the idea is
      genuinely closed rather than merely untested.
  B2  the feature model must beat the CATALOG-POSITION baseline -- "is this position in the deployed
      catalog?" -- which is the real incumbent and is nearly free.
  B3  reported but NOT a bar: which feature family carries it (chemistry vs position vs pretrained ESM2).

HONESTY RAILS:
  * The target is itself a FITTED coefficient, so it carries estimation noise. Tokens are kept only at
    >= MIN_CARRIERS, and every metric is also reported carrier-count-weighted.
  * Leave-one-token-out on a coefficient vector from ONE joint fit leaves residual coupling between
    targets (collinear tokens share signal). Named, not hidden; a fully clean version refits per fold.
  * ESM2 sees no fold-change data, so it is a legitimate pretrained feature and not leakage.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

DATA = ROOT / "data" / "raw" / "hiv" / "NNRTI_DataSet.txt"
RT_REF = ROOT / "data" / "hiv_ref" / "HIV1_RT_HXB2_cds.fna"
ESM_CACHE = Path("D:/dna_decode_cache/hiv/rt_esm2_650M_marginals.json")
OUT = ROOT / "wiki" / "hiv_alphabet_generalization_2026-09-29.json"

AAS = "ACDEFGHIKLMNPQRSTVWY"
MIN_CARRIERS = 10          # the DRMcv.R convention already used by hiv_nnrti_baseline.py
RIDGE_LAMBDA = 1.0

PREREG = {
    "frozen_before_any_number_was_seen": True,
    "B1_beats_mean_coefficient": "feature model Spearman/R2 on held-out token coefficients > predicting "
                                 "the mean (what an identity-keyed model must do on an unseen token)",
    "B2_beats_catalog_position": "feature model > 'is this position in the deployed catalog' baseline",
    "B3_reported_not_a_bar": "which feature family carries the signal",
    "why_not_the_0.911_ceiling": "that ceiling is for SEEN tokens; an unseen token scores exactly 0 under "
                                 "an identity-keyed model, so the ceiling does not bound this task",
}


def load_dataset(drug: str):
    lines = DATA.read_text(encoding="utf-8", errors="replace").splitlines()
    head = lines[0].split("\t")
    di = head.index(drug)
    pos_cols = [(i, int(c[1:])) for i, c in enumerate(head) if c.startswith("P") and c[1:].isdigit()]
    rows = []
    for ln in lines[1:]:
        f = ln.split("\t")
        if len(f) != len(head):
            continue
        raw = f[di].strip()
        if raw in ("", "NA"):
            continue
        try:
            fold = float(raw.lstrip("<>="))
        except ValueError:
            continue
        if fold <= 0:
            continue
        muts = set()
        for i, pos in pos_cols:
            aa = f[i].strip()
            # '-' = wild type. Mixtures are multi-letter; take each letter as present.
            if aa in ("", "-", ".", "X", "*"):
                continue
            for a in aa:
                if a in AAS:
                    muts.add(f"{pos}{a}")
        rows.append({"y": math.log10(fold), "muts": muts})
    return rows


def ridge_fit(X, y, lam=RIDGE_LAMBDA):
    Xa = np.column_stack([X, np.ones(len(X))])
    A = Xa.T @ Xa + lam * np.eye(Xa.shape[1])
    A[-1, -1] -= lam
    return np.linalg.solve(A, Xa.T @ y)


def rt_protein() -> str:
    """Translate the COMMITTED HXB2 RT reference. Verified against the shipped catalog wild-types below."""
    sys.path.insert(0, str(ROOT / "scripts"))
    from pear_genotype_alphabet import CODON
    seq = "".join(l.strip() for l in RT_REF.read_text().splitlines() if not l.startswith(">"))
    prot = "".join(CODON.get(seq[i:i + 3], "X")
                   for i in range(0, len(seq) - len(seq) % 3, 3)).rstrip("*")
    # REFERENCE INTEGRITY: the repo already pins that this translation matches the catalog WT at every
    # DRM position. Re-assert it here rather than trusting the file -- a frame error must fail loudly.
    from dna_decode.data.hiv_amr import _RT_WT
    bad = [(pos, wt, prot[pos - 1]) for pos, wt in _RT_WT.items()
           if pos - 1 < len(prot) and prot[pos - 1] != wt]
    if bad:
        raise SystemExit(f"REFUSING: RT translation disagrees with catalog WT at {bad[:5]} -- frame error.")
    print(f"  RT reference verified against {len(_RT_WT)} catalog wild-types; length {len(prot)}")
    return prot


def esm2_marginals(protein: str) -> dict:
    if ESM_CACHE.exists():
        return json.loads(ESM_CACHE.read_text(encoding="utf-8"))
    import torch
    from transformers import AutoModelForMaskedLM, AutoTokenizer
    name = "facebook/esm2_t33_650M_UR50D"
    tok = AutoTokenizer.from_pretrained(name)
    model = AutoModelForMaskedLM.from_pretrained(name)
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    model = model.to(dev).eval()
    ids = tok(protein, return_tensors="pt")["input_ids"].to(dev)
    out: dict[str, dict] = {}
    with torch.no_grad():
        for i in range(len(protein)):
            m = ids.clone()
            m[0, i + 1] = tok.mask_token_id
            lp = torch.log_softmax(model(input_ids=m).logits[0, i + 1].float(), dim=-1)
            out[str(i + 1)] = {aa: float(lp[tok.convert_tokens_to_ids(aa)]) for aa in AAS}
            if (i + 1) % 100 == 0:
                print(f"    esm2 {i+1}/{len(protein)}", flush=True)
    ESM_CACHE.parent.mkdir(parents=True, exist_ok=True)
    ESM_CACHE.write_text(json.dumps(out), encoding="utf-8")
    return out


def spearman(a, b):
    from scipy.stats import spearmanr
    return float(spearmanr(a, b).statistic)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--drug", default="EFV")
    ap.add_argument("--no-esm", action="store_true", help="skip the pretrained feature (ablation)")
    args = ap.parse_args(argv)

    rows = load_dataset(args.drug)
    counts: dict[str, int] = {}
    for r in rows:
        for m in r["muts"]:
            counts[m] = counts.get(m, 0) + 1
    tokens = sorted([m for m, c in counts.items() if c >= MIN_CARRIERS],
                    key=lambda m: (int(m[:-1]), m[-1]))
    print(f"{args.drug}: {len(rows)} isolates | {len(tokens)} tokens with >= {MIN_CARRIERS} carriers")

    # ---------- stage 1: the incumbent additive identity model -> per-token effect coefficients ----------
    idx = {m: i for i, m in enumerate(tokens)}
    X = np.zeros((len(rows), len(tokens)))
    y = np.array([r["y"] for r in rows])
    for i, r in enumerate(rows):
        for m in r["muts"]:
            if m in idx:
                X[i, idx[m]] = 1.0
    w = ridge_fit(X, y)
    pred = np.column_stack([X, np.ones(len(X))]) @ w
    coef = w[:-1]
    print(f"  stage 1 additive identity model: in-sample Spearman {spearman(pred, y):+.4f} "
          f"(published OOF ~0.911 -- this is a fit, not a held-out claim)")

    # ---------- features per token ----------
    protein = None if args.no_esm else rt_protein()
    esm = None if args.no_esm else esm2_marginals(protein)
    from dna_decode.forward.variant_effect import blosum62_score
    # Deployed catalog positions for THE DRUG'S OWN CLASS. Read from the shipped module, never
    # hand-typed -- but CLASS-MATCHED, which the first version got wrong and is worth recording.
    #
    # THE DEFECT (found 2026-09-29): the first version scanned hiv_amr for ANY collection of ints and
    # unioned them, giving 17 positions -- the 8 NNRTI sites (_RT_WT) PLUS the 9 NRTI sites
    # (NRTI_MAJOR_POSITIONS). Both are in RT, which is exactly why it looked plausible, but NRTI sites
    # (M184V, the 41/210/215 TAMs) are COMMON in treatment-experienced isolates and NNRTI-NEUTRAL. So
    # against an EFV label the contaminated baseline scored Spearman -0.6167 -- strongly NEGATIVE, which
    # is the signature that gave it away. A non-vacuity guard (">= 5 positions") passed the whole time:
    # it checked the baseline was not EMPTY, never that it was the RIGHT one.
    import dna_decode.data.hiv_amr as ha
    CLASS_CATALOG = {
        "NNRTI": ("_RT_WT", ("EFV", "NVP", "ETR", "RPV", "DOR")),
        "NRTI": ("NRTI_MAJOR_POSITIONS", ("3TC", "ABC", "AZT", "D4T", "DDI", "TDF")),
    }
    klass = next((k for k, (_, drugs) in CLASS_CATALOG.items() if args.drug in drugs), None)
    if klass is None:
        raise SystemExit(f"REFUSING: no class-matched catalog for drug {args.drug}; a cross-class "
                         "baseline is not a control (see the defect note above).")
    attr = CLASS_CATALOG[klass][0]
    cat_obj = getattr(ha, attr)
    cat_pos = set(cat_obj)
    other = set()
    for k, (a2, _) in CLASS_CATALOG.items():
        if k != klass:
            other |= set(getattr(ha, a2))
    if cat_pos & other:
        raise SystemExit(f"REFUSING: {attr} overlaps another class's positions {sorted(cat_pos & other)} "
                         "-- the class-match guard cannot separate them.")
    if len(cat_pos) < 5:
        raise SystemExit(f"REFUSING: only {len(cat_pos)} positions in {attr} -- baseline would be vacuous.")
    print(f"  catalog baseline = {klass} via hiv_amr.{attr}: {len(cat_pos)} positions {sorted(cat_pos)}")
    print(f"  EXCLUDED as cross-class (would contaminate): {sorted(other)}")

    feats, tgt, wts, keep = [], [], [], []
    for m in tokens:
        pos, mut = int(m[:-1]), m[-1]
        wt = protein[pos - 1] if protein and pos - 1 < len(protein) else None
        row = [pos / 560.0, math.log10(counts[m])]
        row += [1.0 if a == mut else 0.0 for a in AAS]
        row += [1.0 if (wt == a) else 0.0 for a in AAS]
        row.append(float(blosum62_score(wt, mut)) if wt and wt in AAS and mut in AAS else 0.0)
        if esm is not None and wt:
            e = esm.get(str(pos))
            row.append((e[mut] - e[wt]) if e and mut in e and wt in e else 0.0)
            row.append(-float(np.sum([math.exp(v) * v for v in e.values()])) if e else 0.0)  # entropy-ish
        feats.append(row)
        tgt.append(coef[idx[m]])
        wts.append(counts[m])
        keep.append(m)
    F = np.asarray(feats)
    T = np.asarray(tgt)
    W = np.asarray(wts, dtype=float)
    print(f"  feature matrix {F.shape}")

    # ---------- stage 2: LEAVE-ONE-TOKEN-OUT ----------
    pred_feat, pred_mean, pred_cat = [], [], []
    for i in range(len(keep)):
        tr = np.arange(len(keep)) != i
        wf = ridge_fit(F[tr], T[tr])
        pred_feat.append(float(np.append(F[i], 1.0) @ wf))
        pred_mean.append(float(T[tr].mean()))
        p = int(keep[i][:-1])
        inc = np.array([1.0 if int(k[:-1]) in cat_pos else 0.0 for k in np.array(keep)[tr]])
        wc = ridge_fit(inc.reshape(-1, 1), T[tr])
        pred_cat.append(float(np.array([1.0 if p in cat_pos else 0.0, 1.0]) @ wc))
    pf, pm, pc = map(np.asarray, (pred_feat, pred_mean, pred_cat))

    def r2(p):
        return float(1.0 - np.sum((T - p) ** 2) / np.sum((T - T.mean()) ** 2))

    def wr2(p):
        return float(1.0 - np.sum(W * (T - p) ** 2) / np.sum(W * (T - np.average(T, weights=W)) ** 2))

    # SPEARMAN IS ARTIFACT-PRONE FOR THE GROUP-MEAN BASELINES UNDER LEAVE-ONE-OUT, and this was measured,
    # not assumed. The catalog baseline predicts a 2-group mean. Under LOO, dropping a HIGH-coefficient
    # catalog token lowers the remaining catalog mean, so that token receives a LOWER prediction -- an
    # anti-correlation between a held-out value and its own prediction, induced purely by excluding the
    # point from the mean it is scored against. With only ~23 catalog tokens it dominates: LOO Spearman
    # reads -0.743 while the SAME quantity computed from a single full fit is +0.301 (catalog group mean
    # +0.593 vs elsewhere +0.008, i.e. biologically correct). R2 is near-immune because it is driven by
    # between-group separation. So R2 is the metric of record here and the LOO Spearman is reported
    # alongside its full-fit counterpart, never alone.
    inc_all = np.array([1.0 if int(m[:-1]) in cat_pos else 0.0 for m in keep])
    w_full = ridge_fit(inc_all.reshape(-1, 1), T)
    pc_full = inc_all * w_full[0] + w_full[1]
    res = {
        "feature_model": {"spearman_loo": spearman(pf, T), "r2": r2(pf), "weighted_r2": wr2(pf)},
        "mean_baseline": {"spearman_loo": 0.0, "r2": r2(pm), "weighted_r2": wr2(pm)},
        "catalog_position_baseline": {
            "spearman_loo": spearman(pc, T), "r2": r2(pc), "weighted_r2": wr2(pc),
            "spearman_full_fit": spearman(pc_full, T),
            "loo_spearman_is_an_artifact": True,
            "group_mean_at_catalog": float(T[inc_all == 1].mean()),
            "group_mean_elsewhere": float(T[inc_all == 0].mean()),
            "n_tokens_at_catalog": int(inc_all.sum()),
        },
        "metric_of_record": "r2 (and weighted_r2); the LOO Spearman on a group-mean predictor is "
                            "anti-correlated by construction -- see the comment in the source",
    }
    verdict = {
        "B1_beats_mean_coefficient": res["feature_model"]["r2"] > res["mean_baseline"]["r2"],
        "B2_beats_catalog_position": res["feature_model"]["r2"] > res["catalog_position_baseline"]["r2"],
    }
    print("\n=== LEAVE-ONE-TOKEN-OUT: predicting an UNSEEN mutation's effect ===")
    for k, v in res.items():
        if not isinstance(v, dict):
            continue
        extra = (f"   spearman_full_fit {v['spearman_full_fit']:+.4f}" if "spearman_full_fit" in v else "")
        print(f"  {k:28s} R2 {v['r2']:+.4f}   wR2 {v['weighted_r2']:+.4f}   "
              f"spearman_loo {v['spearman_loo']:+.4f}{extra}")
    print("\n=== PRE-REGISTERED VERDICTS ===")
    for k, v in verdict.items():
        print(f"  {k:32s} {v}")

    OUT.write_text(json.dumps({
        "schema": "hiv-alphabet-generalization-v1", "analysis_date": "2026-09-29",
        "drug": args.drug, "n_isolates": len(rows), "n_tokens": len(keep),
        "min_carriers": MIN_CARRIERS, "esm_used": esm is not None,
        "stage1_identity_model_in_sample_spearman": spearman(pred, y),
        "results": res, "prereg": PREREG, "verdicts": verdict,
        "honest_limits": [
            "The target is a FITTED coefficient from one joint ridge fit, so it carries estimation noise "
            "and residual collinearity between tokens; leave-one-token-out does not remove that coupling.",
            "An identity-keyed model scores an unseen token at exactly 0 by construction, so the "
            "mean-coefficient baseline is the generous stand-in for it.",
            "One drug, one gene, one class. ~96% subtype B, which removes the ancestry confound but also "
            "means nothing here speaks to cross-subtype transfer.",
            "ESM2 sees no fold-change data -- a legitimate pretrained feature, not leakage.",
        ],
    }, indent=1), encoding="utf-8")
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
