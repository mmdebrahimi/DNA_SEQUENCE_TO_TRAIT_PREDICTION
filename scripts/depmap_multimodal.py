"""DepMap MULTIMODAL attribution — mutation vs copy-number vs expression per biomarker (the feature-match test).

The mutation-only DepMap decoder recovered BRAF/TP53 (point-mutation mechanisms) but MISSED the EGFR-TKIs
(erlotinib/gefitinib) because EGFR-inhibitor sensitivity is driven by EGFR AMPLIFICATION + EXPRESSION, not
point mutation. This adds CCLE copy-number + expression and tests, per (drug, biomarker), WHICH feature type
carries the de-confounded signal -- the direct multi-modal proof of "feature type must match mechanism type".

Continuous clade-centered Spearman (residualize feature + LFC by lineage) + within-lineage permutation feel is
inherited from the yeast CNV test. Data on D: (gitignored-class); reuses depmap_decoder for PRISM + lineage.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
import scripts.depmap_decoder as DD  # noqa: E402

# (drug, biomarker gene, CCLE column, expected_direction_of_feature_vs_LFC): sensitivity biomarker -> feature
# UP means MORE sensitive -> LOWER LFC -> NEGATIVE rho; resistance biomarker -> POSITIVE rho.
CASES = [
    ("vemurafenib", "BRAF", "BRAF (673)", "sensitize"),
    ("dabrafenib", "BRAF", "BRAF (673)", "sensitize"),
    ("erlotinib", "EGFR", "EGFR (1956)", "sensitize"),
    ("gefitinib", "EGFR", "EGFR (1956)", "sensitize"),
    ("lapatinib", "ERBB2", "ERBB2 (2064)", "sensitize"),
    ("nutlin-3", "TP53", "TP53 (7157)", "resist"),
]
_COLS = sorted({c[2] for c in CASES})


def _load_modalities(data: Path):
    lfc, treat, cells, mut = DD.load(data)
    r2d = cells.set_index("row_name")["depmap_id"].to_dict()
    tis = cells.set_index("depmap_id")["primary_tissue"].to_dict()
    lfc = lfc.rename(index=r2d); lfc = lfc[~lfc.index.duplicated()]
    G = DD.gene_matrix(mut, set(lfc.index))                                    # mutation (binary), cell x gene
    cn = _read_target_cols(data / "ccle_gene_cn.csv")
    ex = _read_target_cols(data / "ccle_expression.csv")
    return lfc, treat, cells, tis, G, cn, ex


def _read_target_cols(path: Path) -> pd.DataFrame:
    """Read only col 0 (cell-line id -> index) + the target gene columns, by INTEGER position (robust to the
    unnamed-first-column quirk that a name/callable usecols mishandles)."""
    with open(path, encoding="utf-8") as fh:
        cols = fh.readline().rstrip("\n").split(",")            # gene names have no commas -> safe split
    keep = set(_COLS)
    pos = [0] + [i for i, c in enumerate(cols) if c in keep]
    return pd.read_csv(path, usecols=pos, index_col=0)


# De-confounding primitive PROMOTED to the installable package (2026-07-02); this script is a thin runner.
from dna_decode.deconfound import group_centered_association as _clade_centered_rho  # noqa: E402


def run(data: Path) -> dict:
    lfc, treat, cells, tis, G, cn, ex = _load_modalities(data)
    out = {"cases": []}
    for drug, gene, col, kind in CASES:
        cols = DD._drug_cols(treat, drug)
        if not cols:
            continue
        y_all = lfc[cols].mean(axis=1)
        rec = {"drug": drug, "gene": gene, "kind": kind, "modalities": {}}
        # align each modality to the drug-response cells
        for mod, src in (("mutation", G[gene] if gene in G.columns else None),
                         ("copy_number", cn[col] if col in cn.columns else None),
                         ("expression", ex[col] if col in ex.columns else None)):
            if src is None:
                rec["modalities"][mod] = {"available": False}
                continue
            common = y_all.index.intersection(src.index)
            y = y_all.loc[common].astype(float).values
            x = pd.to_numeric(src.loc[common], errors="coerce").astype(float).values
            lineage = np.array([tis.get(c, "unknown") for c in common])
            g, w, n = _clade_centered_rho(y, x, lineage)
            # "correct" = de-confounded rho present + sign matches biology (sensitize->neg, resist->pos)
            want_neg = kind == "sensitize"
            correct = bool(not np.isnan(w) and abs(w) > 0.1 and ((w < 0) == want_neg))
            rec["modalities"][mod] = {"available": True, "n": n, "global_rho": round(g, 3),
                                      "within_lineage_rho": round(w, 3) if not np.isnan(w) else None,
                                      "de_confounded_correct": correct}
        # which modality best carries the de-confounded signal (by |within_lineage_rho|, correct-signed only)
        best = max((m for m in rec["modalities"] if rec["modalities"][m].get("de_confounded_correct")),
                   key=lambda m: abs(rec["modalities"][m]["within_lineage_rho"]), default=None)
        rec["best_modality"] = best
        out["cases"].append(rec)
    return out


def main(argv=None) -> int:
    import argparse
    import json
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=REPO / "wiki" / "depmap_multimodal_scores.json")
    a = ap.parse_args(argv)
    res = run(a.data)
    a.out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    for r in res["cases"]:
        print(f"\n{r['drug']} x {r['gene']} ({r['kind']}): BEST de-confounded modality = {r['best_modality']}")
        for mod, m in r["modalities"].items():
            if m.get("available"):
                mark = "OK" if m["de_confounded_correct"] else "  "
                print(f"   [{mark}] {mod:12} global_rho={m['global_rho']:+.3f} within-lineage_rho={m['within_lineage_rho']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
