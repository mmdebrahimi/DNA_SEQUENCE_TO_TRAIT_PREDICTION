"""AlphaMissense on the held-out MaveDB DMS set — the deployable decoder's molecular-fitness generalization (move 3).

The full ESM2 prospective holdout got median |Spearman| 0.478 over 2383 held-out MaveDB assays (0.492 human).
This places AlphaMissense (the deployable, no-GPU decoder) on the SAME held-out task: for each held-out human
MaveDB assay (gene NOT in the ProteinGym benchmark), correlate AM pathogenicity vs the wet-lab DMS score, with
the MaveDB-shipped offset applied so AM's UniProt-numbered variants align with the assay.

LEAKAGE FRAMING (honest, load-bearing):
  - For ESM2 (zero-shot MLM), "held-out" = the assay's gene is not in ProteinGym (what ESM/ProSST were tuned
    against); ESM never saw any DMS label.
  - AlphaMissense is NOT sequence-held-out — it was trained proteome-wide, so it saw these protein sequences.
    BUT the DMS-FITNESS labels are independent of AM's training (AM trained on population/clinical pathogenicity
    weak-labels, NOT on DMS scores). So AM-vs-DMS-fitness is a fair LABEL-independent generalization test, and
    directly comparable to ESM2's number as "which deployable predictor ranks molecular fitness better."
  - The full leakage-free ESM2+ProSST HYBRID at this scale is GPU-bound -> the named Kaggle follow-up.

  uv run python scripts/mavedb_am_holdout.py --limit 60        # AM |Spearman| over N held-out human assays

Frozen AMR surface byte-unchanged (READ-only).
"""
from __future__ import annotations

import argparse
import json
import statistics as st
import sys
from datetime import date as _date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.clinical_gene_landscape_census import enumerate_human_mavedb, fetch_dms_offset  # noqa: E402
from scripts.clinical_am_hybrid_auroc import load_am, build_am_filter, AM_FILTERED  # noqa: E402
from scripts.mavedb_prospective_holdout import proteingym_gene_symbols  # noqa: E402


def _spearman(a: list[float], b: list[float]) -> float:
    from scipy.stats import spearmanr
    return float(spearmanr(a, b)[0])


def am_holdout_spearman(gene: str, meta: dict) -> dict | None:
    """|Spearman| between AM pathogenicity and -DMS (both oriented higher=damaging) on the offset-aligned join."""
    up = meta.get("uniprot")
    if not up:
        return None
    dms = fetch_dms_offset(meta["urn"], meta.get("offset", 0))
    am = load_am(up)
    shared = sorted(set(dms) & set(am))
    if len(shared) < 20:
        return None
    dvals = [dms[k] for k in shared]
    amvals = [am[k] for k in shared]  # higher = more pathogenic
    rho = _spearman(amvals, dvals)     # sign is assay-direction-dependent -> report |rho|
    return {"gene": gene, "uniprot": up, "urn": meta["urn"], "n": len(shared),
            "abs_spearman": round(abs(rho), 4)}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--limit", type=int, default=60, help="max held-out human genes to score")
    a = ap.parse_args()

    pg = {s.upper() for s in proteingym_gene_symbols()}
    landscape = enumerate_human_mavedb(use_cache=True)
    # held-out human genes with a UniProt, ranked by DMS size (most-informative first), NOT in ProteinGym
    held = [(g, m) for g, m in landscape.items()
            if m.get("uniprot") and g.upper() not in pg and (m.get("n_variants") or 0) >= 500]
    held.sort(key=lambda gm: -(gm[1].get("n_variants") or 0))
    held = held[:a.limit]
    print(f"held-out human MaveDB genes (not in ProteinGym, >=500 var): scoring up to {len(held)}", flush=True)

    ups = {m["uniprot"] for _, m in held}
    have = {ln.split("\t", 1)[0] for ln in AM_FILTERED.open(encoding="utf-8")} if AM_FILTERED.exists() else set()
    if ups - have:
        print(f"extending AM filter for {len(ups)} held-out UniProts (one AM-gz stream) ...", flush=True)
        build_am_filter(ups | have)

    scored = []
    for g, m in held:
        try:
            r = am_holdout_spearman(g, m)
        except Exception as e:  # noqa: BLE001
            r = None
            print(f"  {g}: err {str(e)[:60]}")
        if r:
            scored.append(r)
            print(f"  {g:12s} n={r['n']:5d} |rho|={r['abs_spearman']}", flush=True)

    vals = [r["abs_spearman"] for r in scored]
    median = round(st.median(vals), 4) if vals else None
    art = {"_schema": "mavedb-am-holdout-v1", "date": _date.today().isoformat(),
           "task": "AlphaMissense pathogenicity vs held-out MaveDB DMS fitness (|Spearman|), offset-aligned",
           "comparator_esm2": {"median_abs_spearman": 0.478, "median_human": 0.492,
                               "source": "wiki/mavedb_full_esm2_2026-07-22 (2383 assays / 978 human)"},
           "leakage_note": "AM is NOT sequence-held-out (proteome-wide training); the DMS-FITNESS labels ARE "
                           "independent of AM training. ESM2 is zero-shot. Full leakage-free ESM2+ProSST hybrid "
                           "at scale = Kaggle follow-up.",
           "n_scored": len(scored), "am_median_abs_spearman": median, "results": scored,
           "frozen_surface_changed": False}
    out = Path(f"wiki/mavedb_am_holdout_{_date.today().isoformat()}.json")
    out.write_text(json.dumps(art, indent=2), encoding="utf-8")
    print(f"\nAM median |Spearman| over {len(scored)} held-out human assays: {median}  "
          f"(ESM2 comparator 0.492 human)")
    print(f"artifact: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
