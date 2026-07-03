"""C. elegans ben-1 -> benzimidazole (albendazole) resistance decoder — the deterministic decoder in a NEW
kingdom (multicellular animal / nematode).

This is the same shape as the project's validated resistance cells (fungal ERG11, HIV RT, bacterial AMR): a
DETERMINISTIC determinant scan of a curated causal gene, validated against a measured drug-response phenotype.
Here the gene is ben-1 (beta-tubulin, the canonical benzimidazole target) and the phenotype is the Andersen-lab
albendazole high-throughput-assay (HTA) resistance call across wild C. elegans isolates.

Data (free, MIT / CC): AndersenLab/ce_cb_ct_betatubulin `isotype_variant_summary.tsv` — a fully joined per-wild-
isolate table with ben-1 variant columns (impactful SNV incl. missense F200Y/Q131L/..., structural variants,
low-expression) AND the measured `abz_hta_*_res` resistance calls. On D: (gitignored).

Determinant rule (transparent LoF/impactful scan): a strain is ben-1-determinant-positive if it has an impactful
ben-1 SNV OR a ben-1 structural variant OR low ben-1 expression -> predict RESISTANT. Validated against the HTA
call. The headline is the HONEST confusion: ben-1 variation is a highly SPECIFIC, high-PPV positive predictor
of resistance, but LOW-SENSITIVITY because most HTA-"resistant" wild isolates carry wild-type ben-1 (the
threshold admits marginal/polygenic resistance) — the project's high-spec/low-sens "scrutinize the label"
pattern, quantified by the strong-vs-marginal diagnostic (determinant+ strains are much more strongly resistant).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

PHENO_ASSAYS = ["abz_hta_norm_res", "abz_hta_2018_res", "abz_hta_2024_res"]
PHENO_CONTINUOUS = {"abz_hta_norm_res": "abz_hta_norm_pheno",
                    "abz_hta_2018_res": "abz_hta_2018_pheno", "abz_hta_2024_res": "abz_hta_2024_pheno"}


def load(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def _notna(s: pd.Series) -> pd.Series:
    return s.notna() & (s.astype(str).str.strip().str.upper() != "NA")


def ben1_determinant(v: pd.DataFrame) -> pd.Series:
    """Transparent ben-1 determinant scan: impactful SNV (incl. missense) OR structural variant OR low
    expression. This is the deployable positive predictor of benzimidazole resistance."""
    hi = _notna(v["ben-1_high_impact_SNV"])
    sv = _notna(v["ben-1_SVs"])
    lo = v["low_ben1_exp"] == True                                  # noqa: E712 (pandas Series compare)
    return (hi | sv | lo).rename("ben1_determinant")


def _binarize(y: pd.Series) -> pd.Series:
    return y.map({True: 1, "TRUE": 1, "True": 1, False: 0, "FALSE": 0, "False": 0})


def score(v: pd.DataFrame, det: pd.Series, pheno_col: str) -> dict:
    yb = _binarize(v[pheno_col])
    m = yb.notna()
    yy = yb[m].astype(int).values
    dd = det[m].astype(int).values
    tp = int(((dd == 1) & (yy == 1)).sum()); fp = int(((dd == 1) & (yy == 0)).sum())
    fn = int(((dd == 0) & (yy == 1)).sum()); tn = int(((dd == 0) & (yy == 0)).sum())
    sens = tp / (tp + fn) if (tp + fn) else float("nan")
    spec = tn / (tn + fp) if (tn + fp) else float("nan")
    ppv = tp / (tp + fp) if (tp + fp) else float("nan")
    acc = (tp + tn) / len(yy) if len(yy) else float("nan")
    rec = {"phenotype": pheno_col, "n_scored": int(m.sum()), "n_resistant": int((yy == 1).sum()),
           "n_sensitive": int((yy == 0).sum()), "n_determinant_pos": int(dd.sum()),
           "TP": tp, "FP": fp, "FN": fn, "TN": tn,
           "sensitivity": round(sens, 3), "specificity": round(spec, 3),
           "ppv": round(ppv, 3), "accuracy": round(acc, 3)}
    # strong-vs-marginal diagnostic: are determinant+ resistant strains MORE strongly resistant than
    # determinant- resistant strains? (higher HTA pheno = more resistant). Confirms ben-1 explains STRONG R.
    cont = PHENO_CONTINUOUS.get(pheno_col)
    if cont and cont in v.columns:
        p = pd.to_numeric(v.loc[m, cont], errors="coerce").values
        rpos = p[(dd == 1) & (yy == 1)]; rneg = p[(dd == 0) & (yy == 1)]; s = p[yy == 0]
        rec["strong_vs_marginal"] = {
            "mean_pheno_resistant_determinant_pos": round(float(np.nanmean(rpos)), 2) if len(rpos) else None,
            "mean_pheno_resistant_determinant_neg": round(float(np.nanmean(rneg)), 2) if len(rneg) else None,
            "mean_pheno_sensitive": round(float(np.nanmean(s)), 2) if len(s) else None}
    return rec


def determinant_class_breakdown(v: pd.DataFrame, pheno_col: str) -> dict:
    """How ben-1_clean_call variant classes distribute across R/S (the curated lab call, for cross-check)."""
    yb = _binarize(v[pheno_col]); m = yb.notna()
    sub = v[m].copy(); sub["_r"] = yb[m].astype(int).values
    out = {}
    for cls, g in sub.groupby(v["ben-1_clean_call"]):
        out[str(cls)] = {"n": len(g), "n_resistant": int((g["_r"] == 1).sum())}
    return out


def run(path: Path) -> dict:
    v = load(path)
    det = ben1_determinant(v)
    out = {"n_strains": len(v), "n_ben1_determinant_pos": int(det.sum()), "assays": []}
    for pheno in PHENO_ASSAYS:
        if pheno in v.columns:
            rec = score(v, det, pheno)
            rec["class_breakdown"] = determinant_class_breakdown(v, pheno)
            out["assays"].append(rec)
    return out


def main(argv=None) -> int:
    import argparse
    import json
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data", type=Path, default=Path("D:/dna_decode_cache/celegans_ben1/variants.tsv"))
    ap.add_argument("--out", type=Path, default=REPO / "wiki" / "celegans_ben1_scores.json")
    a = ap.parse_args(argv)
    res = run(a.data)
    a.out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(f"strains={res['n_strains']} ben-1 determinant+={res['n_ben1_determinant_pos']}")
    for r in res["assays"]:
        print(f"\n== {r['phenotype']} (n={r['n_scored']}, {r['n_resistant']}R/{r['n_sensitive']}S) ==")
        print(f"   sens={r['sensitivity']} spec={r['specificity']} PPV={r['ppv']} acc={r['accuracy']} "
              f"(TP={r['TP']} FP={r['FP']} FN={r['FN']} TN={r['TN']})")
        sm = r.get("strong_vs_marginal", {})
        print(f"   strong-vs-marginal pheno: det+ R={sm.get('mean_pheno_resistant_determinant_pos')} | "
              f"det- R={sm.get('mean_pheno_resistant_determinant_neg')} | S={sm.get('mean_pheno_sensitive')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
