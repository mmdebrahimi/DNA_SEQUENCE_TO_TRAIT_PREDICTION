"""AlphaMissense vs the deterministic BLOSUM floor, on the same human DMS benchmark — quantifying the
deterministic -> strong-precomputed-predictor gap for the molecular-phenotype modality.

The multi-gene DMS benchmark (`dms_variant_effect_benchmark.py`) established the DETERMINISTIC floor
(BLOSUM62 ~0.22 median) and CITED published learned models at ~0.4-0.5 — not run (GPU). This turns that cited
contrast into a MEASURED one using AlphaMissense (Cheng et al. 2023, Science): a free PRECOMPUTED variant-effect
predictor (no GPU) covering all human missense variants. AM is a LEARNED (AlphaFold-distilled) predictor — it is
NOT the project's deterministic ethos; it is benchmarked here as the "strong usable predictor" to size the gap.

Join (the load-bearing correctness details):
  * per-assay UniProt id + OFFSET from MAVEDB metadata; UniProt position = DMS position + offset (construct/
    isoform numbering correction). A wrong offset -> near-zero AM match -> flagged by match_rate.
  * AM key = (uniprot_id, "{from1}{uniprot_pos}{to1}") e.g. ("P60484","V133A").
  * AM `am_pathogenicity` is a DAMAGING score (higher = more pathogenic). To make it comparable to BLOSUM's
    "positive = predicts functionality" axis, `am_predictive = -sign * spearman(AM, dms)` where `sign` maps the
    assay's DMS score onto functionality (from the nonsense anchor). So a POSITIVE am_predictive means AM
    correctly ranks damaging variants as low-function -- directly comparable to blosum_corrected.
Data on D: (gitignored).
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
from scripts.proteingym_variant_effect import parse_missense  # noqa: E402


def load_am(am_tsv: Path) -> dict:
    """am_filtered.tsv (no header; cols: uniprot, protein_variant, am_pathogenicity, am_class) -> lookup."""
    lut = {}
    with open(am_tsv, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            try:
                lut[(p[0], p[1])] = float(p[2])
            except ValueError:
                continue
    return lut


def assay_am(csv: Path, uniprot: str, offset: int, am: dict, min_n: int = 30) -> dict | None:
    df = pd.read_csv(csv)
    if "hgvs_pro" not in df.columns or "score" not in df.columns or not uniprot:
        return None
    df["s"] = pd.to_numeric(df["score"], errors="coerce")
    ys, ams = [], []
    n_missense = 0
    for hg, s in zip(df["hgvs_pro"].astype(str), df["s"]):
        if np.isnan(s):
            continue
        mm = parse_missense(hg)
        if not mm:
            continue
        n_missense += 1
        import re
        m = re.match(r"^p\.[A-Za-z]{3}(\d+)[A-Za-z]{3}$", hg)
        if not m:
            continue
        upos = int(m.group(1)) + int(offset or 0)
        val = am.get((uniprot, f"{mm[0]}{upos}{mm[1]}"))
        if val is not None:
            ys.append(float(s))
            ams.append(val)
    if len(ys) < min_n:
        return None
    ys = np.array(ys)
    ter = pd.to_numeric(df[df["hgvs_pro"].astype(str).str.contains("Ter|\\*")]["score"], errors="coerce").dropna()
    if len(ter) >= 5:
        sign = 1.0 if ter.mean() < ys.mean() else -1.0
        polarity = "functional_high" if sign > 0 else "damaging_high"
    else:
        sign, polarity = 1.0, "unknown"
    rho = float(spearmanr(np.array(ams), ys)[0])
    return {"uniprot": uniprot, "n_am_matched": len(ys), "n_missense": n_missense,
            "match_rate": round(len(ys) / max(n_missense, 1), 3), "polarity": polarity,
            "am_predictive_spearman": round(-sign * rho, 4)}


def run(manifest_uniprot: Path, am_tsv: Path, blosum_scores: Path | None = None) -> dict:
    man = json.loads(Path(manifest_uniprot).read_text(encoding="utf-8"))
    am = load_am(am_tsv)
    blos = {}
    if blosum_scores and Path(blosum_scores).exists():
        for r in json.loads(Path(blosum_scores).read_text(encoding="utf-8")).get("assays", []):
            blos[r["urn"]] = r["spearman_polarity_corrected"]
    rows = []
    for m in man:
        r = assay_am(Path(m["file"]), m.get("uniprot"), m.get("offset", 0), am)
        if r and r["match_rate"] >= 0.5:                 # require a real join (offset sanity)
            r.update({"gene": m["gene"], "modality": m["modality"], "urn": m["urn"],
                      "blosum_corrected_spearman": blos.get(m["urn"])})
            rows.append(r)
    am_vals = [r["am_predictive_spearman"] for r in rows if r["polarity"] != "unknown"]
    paired = [(r["blosum_corrected_spearman"], r["am_predictive_spearman"]) for r in rows
              if r["polarity"] != "unknown" and r["blosum_corrected_spearman"] is not None]
    return {
        "n_assays_am_joined": len(rows),
        "n_genes": len({r["gene"].split()[0] for r in rows}),
        "predictor": "AlphaMissense am_pathogenicity (precomputed; free; NOT deterministic -- learned/AF-distilled)",
        "am_median_spearman": round(float(np.median(am_vals)), 4) if am_vals else None,
        "am_mean_spearman": round(float(np.mean(am_vals)), 4) if am_vals else None,
        "blosum_median_spearman_paired": round(float(np.median([b for b, _ in paired])), 4) if paired else None,
        "am_median_spearman_paired": round(float(np.median([a for _, a in paired])), 4) if paired else None,
        "median_gain_am_over_blosum": round(float(np.median([a - b for b, a in paired])), 4) if paired else None,
        "assays": sorted(rows, key=lambda r: -(r["am_predictive_spearman"])),
    }


def main(argv=None) -> int:
    import argparse
    DG = Path("D:/dna_decode_cache/proteingym")
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", type=Path, default=DG / "benchmark_manifest_uniprot.json")
    ap.add_argument("--am", type=Path, default=DG / "am_filtered.tsv")
    ap.add_argument("--blosum", type=Path, default=REPO / "wiki" / "dms_variant_effect_benchmark_scores.json")
    ap.add_argument("--out", type=Path, default=REPO / "wiki" / "dms_alphamissense_benchmark_scores.json")
    a = ap.parse_args(argv)
    res = run(a.manifest, a.am, a.blosum)
    a.out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(f"AM-joined assays={res['n_assays_am_joined']} genes={res['n_genes']}")
    print(f"AM median Spearman={res['am_median_spearman']} | BLOSUM(paired)={res['blosum_median_spearman_paired']}"
          f" | median gain AM over BLOSUM={res['median_gain_am_over_blosum']}")
    print("\nper-assay (AM vs BLOSUM, polarity-corrected):")
    for r in res["assays"]:
        print(f"  AM={r['am_predictive_spearman']:+.3f} BLOSUM={r['blosum_corrected_spearman']} | "
              f"{r['gene']:9} {r['modality']:9} match={r['match_rate']} n={r['n_am_matched']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
