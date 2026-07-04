"""Phase 2 of the hybrid plan — the DECISIVE test: does a supervised head on ESM embeddings BEAT both the
deterministic catalog AND zero-shot ESM on HELD-OUT positions (incl. recovering the NNRTI pocket signal
zero-shot gets backwards)?

Zero-shot ESM = evolutionary FITNESS (blind to pocket-evasion resistance; fair test AUC 0.24 on NNRTI). A
head trained on fold LABELS learns RESISTANCE directly and could recover it. The honest test is
LEAVE-ONE-POSITION-OUT CV (hold out ALL variants at a position; predict them from a model trained on OTHER
positions) — proves generalization to UNSEEN resistance sites, not memorization.

Data: HIV single-mutant isolates (attributable per-variant fold), 4 classes, label R=fold>=3.
Features per variant: [zero-shot LLR, BLOSUM62(wt,mut), ESM WT per-position embedding (1280d)].
Model: StandardScaler -> PCA(30) -> LogisticRegression(balanced). Baselines on the SAME held-out variants:
zero-shot -LLR (AUC) + the deterministic catalog (call_hiv_observed -> R/S, balacc). Weights -> D:.

Pre-committed verdict (per plans/Hybrid_Learned_Deterministic_Decoder_Plan.md Phase 2):
  WIN     : head balacc >= catalog balacc + 0.03  AND  head NNRTI AUC clearly > zero-shot NNRTI AUC (0.24)
  PARTIAL : head pooled AUC > zero-shot pooled AUC  but  head balacc < catalog balacc (novel-variant fallback only)
  FAIL    : head does not beat the catalog on held-out positions (learned scoring adds no deployable value)
"""
from __future__ import annotations

import json
import sys
from datetime import date as _date
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from Bio.Align import substitution_matrices  # noqa: E402
from dna_decode.data.hiv_amr import call_hiv_observed  # noqa: E402
from scripts.esm_catalog_extension_test import _auc  # noqa: E402
from scripts.esm_hiv_resistance_test import CLASSES, _single_mutant_variants, _wt_protein  # noqa: E402
from scripts.hiv_nnrti_validate import load_rows  # noqa: E402

AAS = "ACDEFGHIKLMNPQRSTVWY"
_BLOSUM = substitution_matrices.load("BLOSUM62")


def _blosum(wt: str, mut: str) -> float:
    try:
        return float(_BLOSUM[wt, mut])
    except (KeyError, IndexError):
        return 0.0


def _esm_wt_features(wt: str, model_name: str):
    """One ESM forward pass over WT -> ({pos: {aa: logprob}}, {pos: embedding-vector}). BOS at idx0."""
    import torch
    torch.hub.set_dir("D:/dna_decode_cache/torch")
    import esm
    model, alphabet = getattr(esm.pretrained, model_name)()
    model.eval()
    _, _, toks = alphabet.get_batch_converter()([("wt", wt)])
    with torch.no_grad():
        out = model(toks, repr_layers=[model.num_layers])
    logits = out["logits"][0]
    reps = out["representations"][model.num_layers][0].numpy()
    lp = torch.log_softmax(logits, dim=-1)
    idx = {aa: alphabet.get_idx(aa) for aa in AAS}
    logprobs = {p: {aa: float(lp[p, idx[aa]]) for aa in AAS} for p in range(1, len(wt) + 1)}
    embs = {p: reps[p] for p in range(1, len(wt) + 1)}
    return logprobs, embs


def build_dataset(ref_dir: Path, data_dir: Path, model_name: str) -> list[dict]:
    cache: dict[str, tuple] = {}
    rows = []
    for cls, fname, gene, plen, drug in CLASSES:
        path = data_dir / fname
        if not path.exists():
            continue
        wt = _wt_protein(ref_dir, gene)
        if gene not in cache:
            cache[gene] = _esm_wt_features(wt, model_name)
        logprobs, embs = cache[gene]
        for v in _single_mutant_variants(load_rows(path), wt, plen, drug):
            pos, w, m = v["pos"], v["wt"], v["mut"]
            v.update({"cls": cls, "gene": gene, "drug": drug,
                      "llr": logprobs[pos][m] - logprobs[pos][w],
                      "blosum": _blosum(w, m), "emb": embs[pos],
                      "group": f"{gene}:{pos}"})
            rows.append(v)
    return rows


# fold-column code -> drug NAME that call_hiv_observed accepts (the code alone -> INDETERMINATE -> false all-S)
_DRUG_NAME = {
    "EFV": "efavirenz", "NVP": "nevirapine", "ETR": "etravirine", "RPV": "rilpivirine", "DOR": "doravirine",
    "3TC": "lamivudine", "ABC": "abacavir", "AZT": "zidovudine", "D4T": "stavudine", "DDI": "didanosine",
    "TDF": "tenofovir", "NFV": "nelfinavir", "FPV": "fosamprenavir", "ATV": "atazanavir", "IDV": "indinavir",
    "LPV": "lopinavir", "SQV": "saquinavir", "TPV": "tipranavir", "DRV": "darunavir", "RAL": "raltegravir",
    "EVG": "elvitegravir", "DTG": "dolutegravir", "BIC": "bictegravir", "CAB": "cabotegravir",
}


def _catalog_call(v) -> str:
    obs = {v["gene"]: {f"{v['wt']}{v['pos']}{v['mut']}"}}
    return call_hiv_observed(_DRUG_NAME.get(v["drug"], v["drug"]), obs).prediction


def _balacc(y_true, y_pred) -> float:
    y_true, y_pred = np.asarray(y_true, bool), np.asarray(y_pred, bool)
    tp = int((y_true & y_pred).sum()); fn = int((y_true & ~y_pred).sum())
    tn = int((~y_true & ~y_pred).sum()); fp = int((~y_true & y_pred).sum())
    sens = tp / (tp + fn) if (tp + fn) else 0.0
    spec = tn / (tn + fp) if (tn + fp) else 0.0
    return round((sens + spec) / 2, 3)


def run(ref_dir: Path, data_dir: Path, model_name: str) -> dict:
    from sklearn.decomposition import PCA
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import LeaveOneGroupOut
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler

    data = build_dataset(ref_dir, data_dir, model_name)
    X = np.array([[v["llr"], v["blosum"], *v["emb"]] for v in data])
    y = np.array([v["label_R"] for v in data], dtype=bool)
    groups = np.array([v["group"] for v in data])
    logo = LeaveOneGroupOut()
    oof = np.full(len(data), np.nan)
    for tr, te in logo.split(X, y, groups):
        if len(np.unique(y[tr])) < 2:                       # a fold whose training set is single-class -> skip head
            continue
        pipe = Pipeline([("sc", StandardScaler()),
                         ("pca", PCA(n_components=min(30, len(tr) - 1, X.shape[1]))),
                         ("lr", LogisticRegression(class_weight="balanced", max_iter=2000, C=0.5))])
        pipe.fit(X[tr], y[tr])
        oof[te] = pipe.predict_proba(X[te])[:, 1]
    scored = ~np.isnan(oof)
    for i, v in enumerate(data):
        v["head_proba"] = round(float(oof[i]), 3) if scored[i] else None
        v["catalog_R"] = _catalog_call(v) == "R"

    def metrics(mask):
        idx = np.where(mask & scored)[0]
        if len(idx) < 8 or len(set(y[idx])) < 2:
            return {"n": int(len(idx)), "note": "under-powered"}
        yy = y[idx]
        head_auc = _auc(list(oof[idx][yy]), list(oof[idx][~yy]))
        zs_auc = _auc([-data[i]["llr"] for i in idx if y[i]], [-data[i]["llr"] for i in idx if not y[i]])
        head_bal = _balacc(yy, oof[idx] >= 0.5)
        cat_bal = _balacc(yy, np.array([data[i]["catalog_R"] for i in idx]))
        zs_bal = _balacc(yy, np.array([-data[i]["llr"] for i in idx]) > 0)
        return {"n": int(len(idx)), "n_R": int(yy.sum()),
                "head_auc": round(head_auc, 3) if head_auc is not None else None,
                "zeroshot_auc": round(zs_auc, 3) if zs_auc is not None else None,
                "head_balacc": head_bal, "catalog_balacc": cat_bal, "zeroshot_balacc": zs_bal}

    per_class = {c: metrics(np.array([v["cls"] == c for v in data])) for c in {v["cls"] for v in data}}
    pooled = metrics(np.ones(len(data), bool))
    nn = per_class.get("NNRTI", {})
    # verdict (pre-committed)
    verdict = "FAIL"
    if pooled.get("head_balacc") is not None:
        beats_catalog = pooled["head_balacc"] >= pooled["catalog_balacc"] + 0.03
        recovers_nnrti = (nn.get("head_auc") is not None and nn.get("zeroshot_auc") is not None
                          and nn["head_auc"] > nn["zeroshot_auc"] + 0.15)
        beats_zeroshot = (pooled.get("head_auc") is not None and pooled.get("zeroshot_auc") is not None
                          and pooled["head_auc"] > pooled["zeroshot_auc"])
        if beats_catalog and recovers_nnrti:
            verdict = "WIN"
        elif beats_zeroshot:
            verdict = "PARTIAL"
    return {
        "artifact": "esm_supervised_resistance_head", "schema": "esm-supervised-head-v1",
        "date": _date.today().isoformat(), "model": model_name,
        "design": "leave-one-POSITION-out CV; head=StandardScaler+PCA(30)+LogReg(balanced) on "
                  "[zero-shot LLR, BLOSUM62, ESM WT per-position embedding]; baselines scored on the SAME "
                  "held-out variants (zero-shot -LLR AUC + deterministic catalog balacc). Label = fold>=3.",
        "n_total": len(data), "n_scored_oof": int(scored.sum()), "n_R": int(y.sum()),
        "verdict": verdict, "pooled": pooled, "per_class": per_class,
        "interpretation": {
            "WIN": "the supervised head recovers resistance the catalog + zero-shot miss (incl. NNRTI pocket) "
                   "on UNSEEN positions -> the learned-scoring hybrid branch (V1) is greenlit for Phase 4.",
            "PARTIAL": "the head beats zero-shot but not the curated catalog -> it is a novel-variant FALLBACK "
                       "(score catalog-silent variants), NOT a replacement; integrate cautiously in Phase 4.",
            "FAIL": "the head does not beat the deterministic catalog on held-out positions -> learned scoring "
                    "adds no deployable value here; the curated catalog remains the sole scorer (honest close).",
        }[verdict],
        "honest_caveats": [
            "leave-one-position-out is the HARD, honest CV (generalize to unseen resistance sites); random "
            "split would inflate via memorization",
            "N=293 single-mutant variants (the attributable-label price); in-distribution HIVDB fold",
            "the deterministic catalog is a strong baseline = essentially a known-DRM lookup; the head's value "
            "is generalizing to catalog-SILENT variants",
            "single-class LOPO folds (a position with only R or only S) are skipped for the head (reported n)",
        ],
    }


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ref-dir", type=Path, default=REPO / "data/hiv_ref")
    ap.add_argument("--data-dir", type=Path, default=REPO / "data/raw/hiv")
    ap.add_argument("--model", default="esm2_t33_650M_UR50D")
    ap.add_argument("--out", type=Path, default=REPO / "wiki" / f"esm_supervised_head_result_{_date.today().isoformat()}.json")
    a = ap.parse_args(argv)
    res = run(a.ref_dir, a.data_dir, a.model)
    a.out.write_text(json.dumps(res, indent=2, default=float), encoding="utf-8")
    print(json.dumps({k: res[k] for k in ("verdict", "pooled", "per_class", "n_total", "n_R")}, indent=2, default=float))
    print(f"\n[wrote {a.out}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
