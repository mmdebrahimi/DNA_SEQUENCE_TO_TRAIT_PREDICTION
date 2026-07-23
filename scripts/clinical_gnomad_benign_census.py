"""gnomAD-benign census — break the benign-sparsity wall with a frequency-benign class (move 2 of the data hunt).

The ClinVar census left 18/28 clinical genes SINGLE_CLASS (pathogenic-dominated ClinVar missense). This supplies
the missing BENIGN class from gnomAD common variants (AF >= 1e-4 = the ACMG frequency-benign principle) and
re-scores: positives = ClinVar-PATHOGENIC missense, negatives = gnomAD-frequency-benign missense, both
intersected with the gene's DMS (offset-applied) so every decoder scores the same variants.

**CIRCULARITY (load-bearing — see scripts/gnomad_benign.py):** AlphaMissense was TRAINED with population-common
variants as its benign weak-label, so AM's AUROC on a gnomAD-benign set is CIRCULAR (reported, flagged, NOT
headlined). The FAIR (non-circular) decoders on this benign source are:
  - **DMS-itself ceiling** (wet-lab; never saw gnomAD) — does the molecular assay separate clinical-pathogenic
    from population-common? = the fitness-alignment ceiling, now measurable on constrained genes.
  - **BLOSUM62 floor** (never saw gnomAD).
ESM2 is also fair (self-supervised MLM) but GPU-slow at census scale — deferred (compute per gene on demand).

So the honest deliverable: gnomAD unlocks the DMS-CEILING (fitness-alignment) validation on the constrained
genes ClinVar alone could not score — confirming (or not) that R2 applies broadly — WITHOUT validating the
deployable AM decoder (which stays ClinVar-only, the 4-gene result).

  uv run python scripts/clinical_gnomad_benign_census.py                    # clinical census list
  uv run python scripts/clinical_gnomad_benign_census.py --genes PTEN,KRAS,GCK,TP53,MSH2

Frozen AMR surface byte-unchanged (READ-only).
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date as _date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.forward.variant_effect import blosum62_score  # noqa: E402
from scripts.clinical_variant_effect_validate import fetch_clinvar_missense, auroc, _spearman_sign, MIN_PER_CLASS  # noqa: E402
from scripts.clinical_am_hybrid_auroc import load_am, build_am_filter, AM_FILTERED  # noqa: E402
from scripts.clinical_gene_landscape_census import (  # noqa: E402
    enumerate_human_mavedb, fetch_dms_offset, CLINICAL_GENES, CURATED_UNIPROT,
)
from scripts.gnomad_benign import fetch_gnomad_benign, DEFAULT_AF_MIN  # noqa: E402


def census_gene_gnomad(gene: str, meta: dict, af_min: float) -> dict:
    up = (meta or {}).get("uniprot") or CURATED_UNIPROT.get(gene)
    rec = {"gene": gene, "uniprot": up, "offset": (meta or {}).get("offset", 0), "urn": (meta or {}).get("urn")}
    if not up or not meta:
        rec["state"] = "NO_UNIPROT_OR_DMS"
        return rec
    dms = fetch_dms_offset(meta["urn"], meta.get("offset", 0))
    clin = fetch_clinvar_missense(gene, use_cache=True)
    gnb = fetch_gnomad_benign(gene, af_min=af_min)
    am = load_am(up)

    path = {k for k, v in clin.items() if v == "PATH"} & set(dms)
    benign = (set(gnb) & set(dms)) - path          # exclude contradictory (common AND ClinVar-path) variants
    both = sorted(path | benign)
    labels = [k in path for k in both]              # True = pathogenic (positive)
    n_path, n_benign = sum(labels), len(both) - sum(labels)
    rec.update({"benign_source": "gnomAD", "af_min": af_min, "n_path_clinvar": n_path,
                "n_benign_gnomad": n_benign, "n_scored": len(both)})
    if n_path < MIN_PER_CLASS or n_benign < MIN_PER_CLASS:
        rec["state"] = "STILL_UNDERPOWERED"
        return rec

    # DMS-ceiling (non-circular): orient label-free vs BLOSUM, predict pathogenic = -preserved
    dvals = [dms[k] for k in both]
    blos = [blosum62_score(k[0], k[2]) for k in both]
    dsign = _spearman_sign(dvals, blos) or 1.0
    rec["dms_ceiling_auroc"] = round(auroc(labels, [-(dsign * v) for v in dvals]), 4)
    rec["blosum_floor_auroc"] = round(auroc(labels, [-v for v in blos]), 4)
    # AM (CIRCULAR on gnomAD-benign — reported, flagged, not headlined)
    am_keys = [k for k in both if k in am]
    if len(am_keys) == len(both):
        rec["am_auroc_CIRCULAR"] = round(auroc(labels, [am[k] for k in both]), 4)
    else:
        rec["am_note"] = f"AM covered {len(am_keys)}/{len(both)} — skipped"
    rec["state"] = "SCORED"
    return rec


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genes")
    ap.add_argument("--af-min", type=float, default=DEFAULT_AF_MIN)
    a = ap.parse_args()

    landscape = enumerate_human_mavedb(use_cache=True)
    genes = a.genes.split(",") if a.genes else CLINICAL_GENES
    ups = set()
    for g in genes:
        up = (landscape.get(g, {}) or {}).get("uniprot") or CURATED_UNIPROT.get(g)
        if up:
            ups.add(up)
    have = {ln.split("\t", 1)[0] for ln in AM_FILTERED.open(encoding="utf-8")} if AM_FILTERED.exists() else set()
    if ups - have:
        print(f"extending AM filter for {len(ups)} UniProts ...", flush=True)
        build_am_filter(ups | have)

    results = []
    for g in genes:
        rec = census_gene_gnomad(g, landscape.get(g), a.af_min)
        results.append(rec)
        if rec["state"] == "SCORED":
            circ = rec.get("am_auroc_CIRCULAR")
            print(f"  {g:8s} SCORED n={rec['n_scored']} ({rec['n_path_clinvar']}P/{rec['n_benign_gnomad']}B) "
                  f"DMS-ceiling={rec['dms_ceiling_auroc']} floor={rec['blosum_floor_auroc']} "
                  f"AM(circular)={circ}", flush=True)
        else:
            print(f"  {g:8s} {rec['state']} ({rec.get('n_path_clinvar')}P/{rec.get('n_benign_gnomad')}B)", flush=True)

    scored = [r for r in results if r["state"] == "SCORED"]
    art = {"_schema": "clinical-gnomad-benign-census-v1", "date": _date.today().isoformat(),
           "benign_source": f"gnomAD r4 frequency-benign (AF >= {a.af_min}); positives = ClinVar pathogenic",
           "circularity": "AlphaMissense trained on population-common-as-benign -> AM AUROC on gnomAD-benign is "
                          "CIRCULAR (reported as am_auroc_CIRCULAR, NOT headlined). Fair non-circular decoders: "
                          "DMS-itself ceiling (wet-lab) + BLOSUM62 floor. ESM2 fair but GPU-deferred.",
           "n_genes": len(results), "n_scored": len(scored),
           "scored_genes": [r["gene"] for r in scored], "results": results, "frozen_surface_changed": False}
    out = Path(f"wiki/clinical_gnomad_benign_census_{_date.today().isoformat()}.json")
    out.write_text(json.dumps(art, indent=2), encoding="utf-8")
    print(f"\nSCORED with gnomAD benign: {len(scored)}/{len(results)} -> {[r['gene'] for r in scored]}")
    print(f"artifact: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
