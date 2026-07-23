"""Clinical-gene landscape census — how many actionable human genes can the R2 decoder be validated on?

The data hunt (2026-07-22): MaveDB has 2,798 score sets → 662 distinct human protein-coding genes → 559 with
a UniProt id. The join to ClinVar was previously blocked by a numbering OFFSET (MaveDB assays are often
domain/construct-numbered), but MaveDB SHIPS the offset in the score-set metadata
(`targetGenes[].externalIdentifiers[] {dbName: UniProt, offset: N}`) — so `UniProt_pos = mavedb_pos + offset`
is AUTO-derivable, no manual curation. This censuses the clinically-actionable genes: for each, pick the
best (most-variant) protein_coding assay, apply the offset, join ClinVar path/benign + AlphaMissense, and
report which genes are AUROC-VIABLE (both classes ≥ MIN_PER_CLASS in the DMS-covered region) = the expandable
R2 clinical-decoder substrate beyond TP53/MSH2.

HONEST: the offset unlocks the JOIN; AUROC-viability still depends on per-gene class balance in the
DMS-covered region (a domain-only assay covering a benign-sparse region stays single-class — MLH1's
C-terminal-domain assay is 27R/2B even after the offset). This census reports the honest per-gene state.

  uv run python scripts/clinical_gene_landscape_census.py                 # real MaveDB + ClinVar + AM (D:)
  uv run python scripts/clinical_gene_landscape_census.py --genes TP53,MSH2,MLH1,KRAS

Frozen AMR surface byte-unchanged (READ-only).
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import sys
import urllib.request
from collections import defaultdict
from datetime import date as _date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.forward.variant_effect import blosum62_score  # noqa: E402
from scripts.mavedb_prospective_holdout import parse_hgvs_pro, _fetch_score_sets  # noqa: E402
from scripts.clinical_variant_effect_validate import fetch_clinvar_missense, auroc, MIN_PER_CLASS  # noqa: E402
from scripts.clinical_am_hybrid_auroc import load_am, build_am_filter, AM_FILTERED  # noqa: E402

API = "https://api.mavedb.org/api/v1"
LANDSCAPE_CACHE = Path("D:/dna_decode_cache/mavedb/human_landscape.json")

# Clinically-actionable human genes to census (ACMG-SF-adjacent + well-known DMS genes). The MaveDB scan
# fills in the best assay + offset; a gene absent from MaveDB-human is reported as no-DMS.
CLINICAL_GENES = [
    "TP53", "MSH2", "MLH1", "MSH6", "PMS2", "PTEN", "BRCA1", "BRCA2", "PALB2", "CHEK2", "VHL", "STK11",
    "TSC2", "KRAS", "LDLR", "MTHFR", "CBS", "GCK", "KCNH2", "SCN5A", "CALM1", "G6PD", "F9", "TPMT",
    "NUDT15", "CYP2C9", "RAD51C", "HMGCR",
]

# Curated gene->UniProt fallback: MaveDB's per-assay UniProt cross-ref is inconsistent (a gene's highest-variant
# assay may omit it — e.g. TP53's 00001234 has no UniProt id), so a known clinical UniProt is never lost. Used
# ONLY to fill a missing id; the OFFSET still comes from the chosen assay's metadata (0 when the assay omits it,
# which is correct for a canonical-numbered assay — the join's low-overlap fails safe if that assumption is wrong).
CURATED_UNIPROT = {
    "TP53": "P04637", "MSH2": "P43246", "MLH1": "P40692", "MSH6": "P52701", "PMS2": "P54278",
    "PTEN": "P60484", "BRCA1": "P38398", "BRCA2": "P51587", "PALB2": "Q86YC2", "CHEK2": "O96017",
    "VHL": "P40337", "STK11": "Q15831", "TSC2": "P49815", "KRAS": "P01116", "LDLR": "P01130",
    "MTHFR": "P42898", "CBS": "P35520", "GCK": "P35557", "KCNH2": "Q12809", "SCN5A": "Q14524",
    "CALM1": "P0DP23", "G6PD": "P11413", "F9": "P00740", "TPMT": "P51580", "NUDT15": "Q9NV35",
    "CYP2C9": "P11712", "RAD51C": "O43502", "HMGCR": "P04035",
}


def enumerate_human_mavedb(use_cache: bool = True) -> dict:
    """{GENE_UPPER: {'uniprot','offset','urn','n_variants'}} — best (most-variant) human protein_coding assay
    per gene, with the MaveDB-shipped UniProt offset. Cached to D:."""
    if use_cache and LANDSCAPE_CACHE.exists():
        return json.loads(LANDSCAPE_CACHE.read_text(encoding="utf-8"))
    ss = _fetch_score_sets(100, all_pages=True)
    by_gene: dict[str, list[dict]] = defaultdict(list)
    for s in ss:
        for tg in (s.get("targetGenes") or []):
            if tg.get("category") != "protein_coding":
                continue
            tax = ((tg.get("targetSequence", {}) or {}).get("taxonomy", {}) or {})
            if tax.get("organismName") != "Homo sapiens":
                continue
            up, off = None, 0
            for ei in tg.get("externalIdentifiers", []):
                if (ei.get("identifier", {}) or {}).get("dbName") == "UniProt":
                    up = ei["identifier"]["identifier"]
                    off = int(ei.get("offset", 0) or 0)
            name = (tg.get("name") or "").upper()
            by_gene[name].append({"uniprot": up, "offset": off, "urn": s.get("urn"),
                                  "n_variants": s.get("numVariants") or 0})
    best: dict[str, dict] = {}
    for name, lst in by_gene.items():
        b = max(lst, key=lambda r: r["n_variants"])
        # normalize gene name to the symbol (first token; MaveDB names like "BRCA1 RING domain")
        sym = name.split()[0] if name else name
        if sym not in best or b["n_variants"] > best[sym]["n_variants"]:
            best[sym] = b
    LANDSCAPE_CACHE.parent.mkdir(parents=True, exist_ok=True)
    LANDSCAPE_CACHE.write_text(json.dumps(best), encoding="utf-8")
    return best


def fetch_dms_offset(urn: str, offset: int) -> dict[tuple[str, int, str], float]:
    """{(wt, UniProt_pos, alt): score} — single-missense DMS with the MaveDB offset applied (UniProt_pos = mavedb_pos + offset)."""
    with urllib.request.urlopen(f"{API}/score-sets/{urn}/scores", timeout=120) as r:
        text = r.read().decode("utf-8", "replace")
    out: dict[tuple[str, int, str], float] = {}
    for row in csv.DictReader(io.StringIO(text)):
        p = parse_hgvs_pro(row.get("hgvs_pro", ""))
        raw = row.get("score", "")
        if p is None or raw in ("", "NA", None):
            continue
        try:
            out[(p[0], p[1] + offset, p[2])] = float(raw)
        except ValueError:
            continue
    return out


def census_gene(gene: str, meta: dict) -> dict:
    up = meta.get("uniprot") or CURATED_UNIPROT.get(gene)
    uniprot_source = "mavedb" if meta.get("uniprot") else ("curated" if up else None)
    rec = {"gene": gene, "uniprot": up, "uniprot_source": uniprot_source,
           "offset": meta.get("offset", 0), "urn": meta.get("urn"), "n_variants": meta.get("n_variants")}
    if not up:
        rec["state"] = "NO_UNIPROT_ID"  # MaveDB assay lacks a UniProt cross-ref AND gene not in curated map
        return rec
    dms = fetch_dms_offset(meta["urn"], meta.get("offset", 0))
    clin = fetch_clinvar_missense(gene, use_cache=True)
    am = load_am(up)
    shared = sorted(set(dms) & set(clin) & set(am))
    n_path = sum(1 for k in shared if clin[k] == "PATH")
    n_benign = len(shared) - n_path
    rec.update({"n_clinvar": len(clin), "n_am": len(am), "n_joined": len(shared),
                "n_path": n_path, "n_benign": n_benign})
    if not am:
        rec["state"] = "NO_AM_COVERAGE"
    elif n_path >= MIN_PER_CLASS and n_benign >= MIN_PER_CLASS:
        labels = [clin[k] == "PATH" for k in shared]
        rec["am_auroc"] = round(auroc(labels, [am[k] for k in shared]), 4)
        rec["blosum_floor_auroc"] = round(auroc(labels, [-blosum62_score(k[0], k[2]) for k in shared]), 4)
        rec["state"] = "AUROC_VIABLE"
    elif len(shared) == 0:
        rec["state"] = "NO_JOIN"
    else:
        rec["state"] = "SINGLE_CLASS"
    return rec


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genes", help="comma-separated gene subset (default: the clinical census list)")
    ap.add_argument("--no-cache", action="store_true", help="force fresh MaveDB landscape enumeration")
    a = ap.parse_args()

    print("enumerating MaveDB human protein_coding landscape ...", flush=True)
    landscape = enumerate_human_mavedb(use_cache=not a.no_cache)
    print(f"  human genes with a DMS: {len(landscape)} ({sum(1 for v in landscape.values() if v.get('uniprot'))} with UniProt)")

    genes = a.genes.split(",") if a.genes else CLINICAL_GENES
    # ensure AM filter covers all the UniProts we will census (mavedb-derived OR curated fallback; one re-stream)
    ups = set()
    for g in genes:
        up = (landscape.get(g, {}) or {}).get("uniprot") or CURATED_UNIPROT.get(g)
        if up:
            ups.add(up)
    have = set()
    if AM_FILTERED.exists():
        have = {ln.split("\t", 1)[0] for ln in AM_FILTERED.open(encoding="utf-8")}
    missing = ups - have
    if missing:
        print(f"  extending AM filter for {len(ups)} UniProts (streaming the AM gz once) ...", flush=True)
        build_am_filter(ups | have)

    results = []
    for g in genes:
        if g not in landscape:
            results.append({"gene": g, "state": "NO_MAVEDB_HUMAN_DMS"})
            print(f"  {g:8s} NO_MAVEDB_HUMAN_DMS")
            continue
        rec = census_gene(g, landscape[g])
        results.append(rec)
        extra = ""
        if rec["state"] == "AUROC_VIABLE":
            extra = f" AM={rec['am_auroc']} floor={rec['blosum_floor_auroc']} (n={rec['n_joined']}, {rec['n_path']}P/{rec['n_benign']}B)"
        elif rec.get("n_joined") is not None:
            extra = f" (n={rec['n_joined']}, {rec.get('n_path')}P/{rec.get('n_benign')}B, off={rec['offset']})"
        print(f"  {g:8s} {rec['state']}{extra}", flush=True)

    viable = [r for r in results if r["state"] == "AUROC_VIABLE"]
    art = {"_schema": "clinical-gene-landscape-census-v1", "date": _date.today().isoformat(),
           "mavedb_human_genes_with_dms": len(landscape),
           "mavedb_human_genes_with_uniprot": sum(1 for v in landscape.values() if v.get("uniprot")),
           "offset_source": "MaveDB targetGenes[].externalIdentifiers[] UniProt offset (auto-derived; UniProt_pos = mavedb_pos + offset)",
           "n_censused": len(genes), "n_auroc_viable": len(viable),
           "auroc_viable_genes": [r["gene"] for r in viable], "results": results,
           "frozen_surface_changed": False}
    out = Path(f"wiki/clinical_gene_landscape_census_{_date.today().isoformat()}.json")
    out.write_text(json.dumps(art, indent=2), encoding="utf-8")
    print(f"\nAUROC-viable clinical genes: {len(viable)} -> {[r['gene'] for r in viable]}")
    print(f"artifact: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
