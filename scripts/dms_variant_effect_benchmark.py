"""Multi-gene DMS variant-effect benchmark — the deterministic substitution-severity baseline across the human
molecular-phenotype landscape.

Extends the single PTEN cell (`proteingym_variant_effect.py`) to ~24 human Deep Mutational Scanning assays
(MAVEDB) spanning THREE phenotype modalities — abundance / function / binding — across 15 genes (incl. the PGx
bridges CYP2C19 / CYP2C9 / TPMT / NUDT15 / GCK). For each assay it correlates a DETERMINISTIC substitution
score (BLOSUM62, authoritative from Biopython; BLOSUM45/80 for robustness) with the measured effect, and asks:
does the simple substitution-severity rule capture protein-variant effect, and does that vary by modality?

METHODOLOGY (the load-bearing fix): DMS assays have OPPOSITE polarities — high score = functional in some
(PTEN VAMP-seq), = damaging in others (MSH2 LOF). Each assay's polarity is anchored by its NONSENSE variants
(maximally damaging): whichever extreme the nonsense mean sits at defines "damaging". The reported
`spearman_polarity_corrected` is signed so that POSITIVE always means "conservative substitution -> preserved
function" — making all assays comparable + aggregatable. Assays without a nonsense anchor report `|rho|` +
polarity=unknown. Data on D: (gitignored).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from scripts.proteingym_variant_effect import parse_missense  # noqa: E402  (reuse the HGVS parser)


def _matrix(name: str):
    from Bio.Align import substitution_matrices
    return substitution_matrices.load(name)


def _score(mat, a: str, b: str) -> float:
    try:
        return float(mat[a, b])
    except (KeyError, IndexError):
        return float(mat[b, a])


def assay_result(csv: Path, gene: str = "", modality: str = "", min_n: int = 30) -> dict | None:
    df = pd.read_csv(csv)
    if "hgvs_pro" not in df.columns or "score" not in df.columns:
        return None
    df["s"] = pd.to_numeric(df["score"], errors="coerce")
    mats = {n: _matrix(n) for n in ("BLOSUM62", "BLOSUM45", "BLOSUM80")}
    feats = {n: [] for n in mats}
    ys = []
    for hg, s in zip(df["hgvs_pro"].astype(str), df["s"]):
        if np.isnan(s):
            continue
        mm = parse_missense(hg)
        if mm:
            for n, mat in mats.items():
                feats[n].append(_score(mat, mm[0], mm[1]))
            ys.append(float(s))
    if len(ys) < min_n:
        return None
    ys = np.array(ys)
    # polarity anchor: nonsense variants are maximally damaging
    ter = pd.to_numeric(df[df["hgvs_pro"].astype(str).str.contains("Ter|\\*")]["score"], errors="coerce").dropna()
    if len(ter) >= 5:
        # if nonsense mean < missense mean -> high score = functional -> conservative(high BLOSUM) => positive corr
        polarity = "functional_high" if ter.mean() < ys.mean() else "damaging_high"
        sign = 1.0 if polarity == "functional_high" else -1.0
    else:
        polarity, sign = "unknown", 1.0
    out = {"gene": gene, "modality": modality, "n_missense": len(ys), "n_nonsense": int(len(ter)),
           "polarity": polarity}
    for n in mats:
        rho = float(spearmanr(np.array(feats[n]), ys)[0])
        out[f"spearman_{n.lower()}_raw"] = round(rho, 4)
    out["spearman_polarity_corrected"] = round(sign * out["spearman_blosum62_raw"], 4)  # headline (BLOSUM62)
    return out


def run(manifest_path: Path, min_n: int = 30) -> dict:
    manifest = json.loads(Path(manifest_path).read_text(encoding="utf-8"))
    rows = []
    for m in manifest:
        r = assay_result(Path(m["file"]), m.get("gene", ""), m.get("modality", ""), min_n=min_n)
        if r:
            r["urn"] = m["urn"]
            rows.append(r)
    corr = [r["spearman_polarity_corrected"] for r in rows if r["polarity"] != "unknown"]
    by_mod = {}
    for mod in ("abundance", "function", "binding"):
        vals = [r["spearman_polarity_corrected"] for r in rows
                if r["modality"] == mod and r["polarity"] != "unknown"]
        if vals:
            by_mod[mod] = {"n_assays": len(vals), "median_spearman": round(float(np.median(vals)), 4),
                           "mean_spearman": round(float(np.mean(vals)), 4)}
    return {
        "n_assays": len(rows),
        "n_genes": len({r["gene"].split()[0] for r in rows}),
        "deterministic_feature": "BLOSUM62 substitution score (polarity-corrected via nonsense anchor)",
        "overall_median_spearman": round(float(np.median(corr)), 4) if corr else None,
        "overall_mean_spearman": round(float(np.mean(corr)), 4) if corr else None,
        "by_modality": by_mod,
        "learned_contrast": "published ProteinGym zero-shot models (ESM-1v/EVE) reach ~0.4-0.5 Spearman on "
                            "such assays (CITED, not run -- GPU). BLOSUM is the deterministic floor.",
        "assays": sorted(rows, key=lambda r: -r["spearman_polarity_corrected"]),
    }


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", type=Path,
                    default=Path("D:/dna_decode_cache/proteingym/benchmark_manifest.json"))
    ap.add_argument("--out", type=Path, default=REPO / "wiki" / "dms_variant_effect_benchmark_scores.json")
    a = ap.parse_args(argv)
    res = run(a.manifest)
    a.out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(f"assays={res['n_assays']} genes={res['n_genes']} "
          f"overall median Spearman={res['overall_median_spearman']}")
    print("by modality:", json.dumps(res["by_modality"]))
    print("\nper-assay (polarity-corrected BLOSUM62 Spearman):")
    for r in res["assays"]:
        print(f"  {r['spearman_polarity_corrected']:+.3f} | {r['gene']:9} {r['modality']:9} "
              f"n={r['n_missense']:5} pol={r['polarity']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
