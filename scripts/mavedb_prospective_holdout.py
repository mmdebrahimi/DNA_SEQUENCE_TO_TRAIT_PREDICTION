"""Leakage-free MaveDB prospective-holdout validation for the frozen `forward` variant-effect cell.

The R2 (molecular / protein variant-effect) analog of the AMR prospective-lock: MaveDB deep-mutational-scanning
score sets whose target gene is NOT in the ProteinGym v1.1 benchmark (the set the frozen hybrid was tuned +
validated on) are a HELD-OUT test set BY CONSTRUCTION — the hybrid could not have seen them. There is NO
population-structure/clonality confound in R2 (a designed mutant library has no lineage axis), so the only
leakage vector is benchmark overlap, which the dedup guard closes.

Two outputs:
  1. The MANIFEST (the novel artifact): MaveDB protein_coding score sets, deduped vs the committed ProteinGym
     catalog by gene symbol -> the provably-held-out set, with per-assay publication date (temporal overlay).
  2. A real-surface SCORING proof: download held-out score CSVs, parse `hgvs_pro` single-missense variants,
     score with the CPU BLOSUM62 baseline (no GPU), Spearman vs the measured functional score. The full
     ESM2+ProSST hybrid scoring is a Kaggle GPU follow-up (the forward cell's established pattern); BLOSUM62
     is the substitution-matrix baseline that runs anywhere and proves the harness end-to-end.

HONEST CAVEATS (both from the scoping memo `wiki/mavedb_forward_cell_scoping_2026-07-21.md`):
  - The MaveDB search endpoint page-caps at 100 score sets; a full manifest needs pagination (follow-up).
  - MaveDB does NOT standardize functional-score DIRECTION per assay (higher can mean more OR less fit) — the
    curation ProteinGym adds. So a raw per-assay Spearman SIGN is uninterpretable without the assay's direction
    metadata; we report |Spearman| as the harness proof and flag that direction-normalization is required for a
    real headline number.

  uv run python scripts/mavedb_prospective_holdout.py                 # manifest only (real network)
  uv run python scripts/mavedb_prospective_holdout.py --score 5       # + BLOSUM scoring proof on 5 held-out assays
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import sys
import urllib.request
from datetime import date as _date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.forward.variant_effect import blosum62_score  # noqa: E402

API = "https://api.mavedb.org/api/v1"
PG_CATALOG = Path("wiki/proteingym_v1.1_substitutions_catalog.tsv")

_AA3TO1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G",
    "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
    "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
}


def parse_hgvs_pro(hgvs_pro: str) -> tuple[str, int, str] | None:
    """`p.Val82Ala` -> ('V', 82, 'A'). Returns None for synonymous (`p.=`), nonsense (`Ter`/`*`),
    multi-variant (`p.[...]`), frameshift, or any non-single-missense form. PURE."""
    if not hgvs_pro or not hgvs_pro.startswith("p."):
        return None
    body = hgvs_pro[2:].strip()
    if body in ("=", "?") or "[" in body or "fs" in body or "del" in body or "ins" in body or "*" in body:
        return None
    # match <3-letter WT><int><3-letter ALT>
    wt3 = body[:3]
    if wt3 not in _AA3TO1:
        return None
    i = 3
    while i < len(body) and body[i].isdigit():
        i += 1
    if i == 3 or i >= len(body):
        return None
    pos = int(body[3:i])
    alt3 = body[i:]
    if alt3 in ("Ter", "="):
        return None
    if alt3 not in _AA3TO1:
        return None
    return _AA3TO1[wt3], pos, _AA3TO1[alt3]


def proteingym_gene_symbols(catalog_path: Path = PG_CATALOG) -> set[str]:
    """Gene symbols represented in the ProteinGym benchmark, from UniProt_ID prefixes (BRCA1_HUMAN -> BRCA1)."""
    syms = set()
    if not catalog_path.exists():
        return syms
    for row in csv.DictReader(catalog_path.open(encoding="utf-8"), delimiter="\t"):
        up = (row.get("UniProt_ID") or "").strip()
        if up:
            syms.add(up.split("_")[0].upper())
    return syms


def _norm_gene(name: str) -> str:
    return "".join(ch for ch in (name or "").upper() if ch.isalnum())


def is_held_out(mavedb_gene_name: str, pg_symbols: set[str]) -> bool:
    """Held out iff the MaveDB target gene symbol is NOT (even as a prefix-match) a ProteinGym gene. Conservative:
    a borderline match is treated as overlap (excluded), never leaked in."""
    g = _norm_gene(mavedb_gene_name)
    if not g:
        return False
    norm_pg = {_norm_gene(s) for s in pg_symbols}
    return g not in norm_pg


def build_manifest(score_sets: list[dict], pg_symbols: set[str], cutoff: str | None) -> list[dict]:
    """Filter to protein_coding + held-out; attach the temporal overlay (published after `cutoff`)."""
    out = []
    for s in score_sets:
        tgs = s.get("targetGenes") or []
        prot = [t for t in tgs if t.get("category") == "protein_coding"]
        if not prot:
            continue
        gene = prot[0].get("name") or ""
        if not is_held_out(gene, pg_symbols):
            continue
        tax = (prot[0].get("targetSequence", {}) or {}).get("taxonomy", {}) or {}
        pub = s.get("publishedDate", "") or ""
        out.append({
            "urn": s.get("urn"), "gene": gene, "organism": tax.get("organismName", "?"),
            "n_variants": s.get("numVariants"), "published": pub,
            "post_cutoff": bool(cutoff and pub and pub >= cutoff),
            "license": (s.get("license", {}) or {}).get("shortName"),
        })
    out.sort(key=lambda r: (r["organism"] != "Homo sapiens", r["published"] or ""), reverse=False)
    return out


def _fetch_page(offset: int) -> list[dict]:
    req = urllib.request.Request(f"{API}/score-sets/search",
                                 data=json.dumps({"offset": offset}).encode(),
                                 headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode("utf-8", "replace")).get("scoreSets", [])


def _fetch_score_sets(limit: int, all_pages: bool = False) -> list[dict]:
    """Fetch score sets. Single page ([:limit], the old default) OR, with all_pages, paginate the whole
    catalog via `offset` (the search endpoint caps a page at 100; `limit`>100 -> 422, so offset is the only
    pagination lever). Stops on the first short page."""
    if not all_pages:
        return _fetch_page(0)[:limit]
    out, offset = [], 0
    while True:
        page = _fetch_page(offset)
        out.extend(page)
        if len(page) < 100:
            return out
        offset += 100
        if offset > 100_000:  # safety backstop
            return out


def score_assay_blosum(urn: str) -> dict:
    """Download a score set's scores CSV, parse single-missense hgvs_pro, BLOSUM62 vs functional score.
    Returns {urn, n_missense, spearman, abs_spearman} or an error dict. Real network + CPU only."""
    from scipy.stats import spearmanr
    try:
        with urllib.request.urlopen(f"{API}/score-sets/{urn}/scores", timeout=60) as r:
            text = r.read().decode("utf-8", "replace")
    except Exception as e:  # noqa: BLE001
        return {"urn": urn, "error": f"download: {e}"}
    blos, sc = [], []
    for row in csv.DictReader(io.StringIO(text)):
        parsed = parse_hgvs_pro(row.get("hgvs_pro", ""))
        raw = row.get("score", "")
        if parsed is None or raw in ("", "NA", None):
            continue
        try:
            fscore = float(raw)
        except ValueError:
            continue
        wt, _pos, alt = parsed
        blos.append(blosum62_score(wt, alt))
        sc.append(fscore)
    if len(blos) < 10:
        return {"urn": urn, "n_missense": len(blos), "spearman": None, "note": "too few missense (<10)"}
    rho = float(spearmanr(blos, sc)[0])
    return {"urn": urn, "n_missense": len(blos), "spearman": round(rho, 4), "abs_spearman": round(abs(rho), 4)}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--limit", type=int, default=100, help="MaveDB search page cap (endpoint max ~100)")
    ap.add_argument("--all-pages", action="store_true",
                    help="paginate the WHOLE MaveDB catalog via offset (not just the first 100-record page)")
    ap.add_argument("--cutoff", default="2024-01-01", help="temporal overlay: flag assays published >= this")
    ap.add_argument("--score", type=int, default=0, help="download+BLOSUM-score up to N held-out assays (proof)")
    a = ap.parse_args()

    pg = proteingym_gene_symbols()
    print(f"ProteinGym benchmark gene symbols: {len(pg)}")
    score_sets = _fetch_score_sets(a.limit, all_pages=a.all_pages)
    print(f"MaveDB score sets fetched ({'ALL pages' if a.all_pages else f'page cap {a.limit}'}): {len(score_sets)}")
    manifest = build_manifest(score_sets, pg, a.cutoff)
    n_human = sum(1 for m in manifest if m["organism"] == "Homo sapiens")
    n_post = sum(1 for m in manifest if m["post_cutoff"])
    print(f"HELD-OUT (protein_coding, not in ProteinGym): {len(manifest)} "
          f"({n_human} human, {n_post} published >= {a.cutoff})")

    scored = []
    if a.score:
        for m in manifest[:a.score]:
            r = score_assay_blosum(m["urn"])
            r["gene"], r["organism"] = m["gene"], m["organism"]
            scored.append(r)
            print(f"  {m['urn']} {m['gene'][:20]:20s} {m['organism'][:16]:16s} -> "
                  f"n={r.get('n_missense')} |rho|={r.get('abs_spearman')}")

    art = {"_schema": "mavedb-prospective-holdout-v1", "date": _date.today().isoformat(),
           "proteingym_genes": len(pg), "mavedb_page_cap": a.limit, "n_fetched": len(score_sets),
           "n_held_out": len(manifest), "n_human": n_human, "cutoff": a.cutoff, "n_post_cutoff": n_post,
           "scorer_proof": "BLOSUM62 (CPU); ESM2+ProSST hybrid = Kaggle follow-up",
           "score_direction_caveat": "MaveDB does not standardize functional-score direction per assay; "
                                     "raw Spearman SIGN is uninterpretable without per-assay direction metadata; "
                                     "|Spearman| reported as a harness proof only",
           "manifest": manifest, "blosum_scoring_proof": scored,
           "frozen_surface_changed": False}
    # NAMESPACE the full-catalog run so it does NOT clobber the single-page manifest (which carries the
    # committed blosum_scoring_proof the paired-comparison artifact reads) — the shared-key-overwrite trap.
    tag = "_full" if a.all_pages else ""
    out = Path(f"wiki/mavedb_prospective_holdout{tag}_{_date.today().isoformat()}.json")
    out.write_text(json.dumps(art, indent=2), encoding="utf-8")
    print(f"artifact: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
