"""E3 learned complement — does a supervised classifier on cheap sequence features beat the conserved-core
decoder on essentiality? Cheap-first (no GPU): aa-composition + protein length + the conserved-core score,
5-fold stratified CV on the Goodall (E. coli) gold-standard. Compares 3 models to the deterministic baseline.
"""
from __future__ import annotations
import gzip, json
import numpy as np
import openpyxl
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from dna_decode.essentiality.core_decoder import score_gene

TGT = "D:/dna_decode_cache/essentiality"
AA = "ACDEFGHIKLMNPQRSTVWY"


def load_labels():
    wb = openpyxl.load_workbook(f"{TGT}/goodall_TableS1_essential.xlsx", read_only=True)
    rows = list(wb.active.iter_rows(values_only=True))
    lab = {}
    for r in rows[2:]:
        g = r[0]
        if not g: continue
        ess = str(r[3]).strip().lower() == "true"; non = str(r[4]).strip().lower() == "true"
        unc = str(r[5]).strip().lower() == "true"
        if unc or (ess == non): continue
        lab[g.strip()] = 1 if ess else 0
    return lab


def main():
    lab = load_labels()
    # NP_accession -> (symbol, product, product_length) from the feature table
    acc2 = {}
    with gzip.open(f"{TGT}/ecoli_k12_feature_table.txt.gz", "rt") as f:
        h = f.readline().rstrip("\n").split("\t")
        ip = h.index("product_accession"); isym = h.index("symbol"); ina = h.index("name"); ipl = h.index("product_length")
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) > ipl and p[0] == "CDS" and p[ip] and p[isym]:
                try: pl = int(p[ipl])
                except ValueError: pl = 0
                acc2[p[ip]] = (p[isym], p[ina], pl)
    # protein sequences by NP accession
    seq = {}; cur = None; buf = []
    with gzip.open(f"{TGT}/ecoli_k12_protein.faa.gz", "rt") as f:
        for line in f:
            if line.startswith(">"):
                if cur: seq[cur] = "".join(buf)
                cur = line[1:].split()[0]; buf = []
            else: buf.append(line.strip())
        if cur: seq[cur] = "".join(buf)
    # build feature matrix for genes with a label + sequence
    X_aa, X_len, X_core, y, genes = [], [], [], [], []
    for acc, (sym, prod, pl) in acc2.items():
        if sym not in lab or acc not in seq: continue
        s = seq[acc]; n = len(s)
        if n < 10: continue
        comp = [s.count(a)/n for a in AA]
        X_aa.append(comp); X_len.append(np.log1p(n))
        X_core.append(score_gene(sym, prod).core_score); y.append(lab[sym]); genes.append(sym)
    X_aa = np.array(X_aa); X_len = np.array(X_len).reshape(-1,1); X_core = np.array(X_core).reshape(-1,1)
    y = np.array(y)
    print(f"training set: {len(y)} genes ({int(y.sum())} essential / {int((1-y).sum())} non); base rate {y.mean():.3f}")

    def cv_auroc(X, model_fn):
        skf = StratifiedKFold(5, shuffle=True, random_state=0)
        oof = np.zeros(len(y))
        for tr, te in skf.split(X, y):
            m = model_fn(); m.fit(X[tr], y[tr]); oof[te] = m.predict_proba(X[te])[:,1]
        return roc_auc_score(y, oof), oof
    def recall_at_spec(scores, target_spec=0.98):
        # threshold where specificity >= target; report recall there
        neg = scores[y==0]; thr = np.quantile(neg, target_spec)
        return float((scores[y==1] >= thr).mean())

    LR = lambda: LogisticRegression(max_iter=2000, class_weight="balanced")
    GB = lambda: HistGradientBoostingClassifier(max_iter=200, class_weight="balanced", random_state=0)
    results = {}
    # baseline: conserved-core score alone (deterministic, no training) -> AUROC
    core_auroc = roc_auc_score(y, X_core.ravel())
    results["conserved_core_baseline"] = {"auroc": round(float(core_auroc),4),
        "recall_at_spec98": round(recall_at_spec(X_core.ravel()),4), "trained": False}
    # learned models
    for name, X in [("aa_composition_LR", X_aa), ("aa+len+core_LR", np.hstack([X_aa,X_len,X_core])),
                    ("aa+len+core_GBM", np.hstack([X_aa,X_len,X_core]))]:
        fn = GB if name.endswith("GBM") else LR
        au, oof = cv_auroc(X, fn)
        results[name] = {"auroc": round(float(au),4), "recall_at_spec98": round(recall_at_spec(oof),4), "trained": True, "cv": "5fold"}
    print("\n=== E3 learned-complement vs conserved-core (E. coli, Goodall gold-standard) ===")
    for k,v in results.items():
        print(f"  {k:<24} AUROC={v['auroc']:.3f}  recall@spec98={v['recall_at_spec98']:.3f}  {'(trained)' if v['trained'] else '(deterministic baseline)'}")
    best=max((v['auroc'],k) for k,v in results.items() if results[k]['trained'])
    lift=best[0]-core_auroc
    verdict = ("LEARNED_COMPLEMENT_EARNS_KEEP" if lift>0.02 else "MARGINAL" if lift>0 else "NO_LIFT")
    out={"organism":"E. coli K-12","reference":"Goodall 2018 mBio Table S1 (TraDIS)","n":len(y),
         "base_rate":round(float(y.mean()),4),"results":results,
         "best_learned":best[1],"auroc_lift_vs_core":round(float(lift),4),"verdict":verdict}
    json.dump(out, open("wiki/essentiality_e3_learned_2026-07-28.json","w"), indent=2)
    print(f"\nbest learned = {best[1]} (AUROC {best[0]:.3f}); lift vs conserved-core {lift:+.3f} -> {verdict}")
    print("wrote wiki/essentiality_e3_learned_2026-07-28.json")


if __name__ == "__main__":
    main()
