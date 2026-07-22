"""Clinical-significance validation of the R2 `forward` variant-effect cell on actionable human genes.

The R2 (molecular / protein variant-effect) cell is validated on DMS *fitness* (Spearman vs the wet-lab
functional score). This asks the CLINICAL question the epoch is about: does the decoder's variant-effect
score separate ClinVar **pathogenic vs benign** missense variants on the genes clinicians actually interpret
(TP53 / MSH2 / PTEN / BRCA1)? That is the VUS-resolution payoff of R2-for-humans.

THREE numbers per gene, each honest about what it is:
  1. **DMS-functional AUROC** — the FITNESS-ALIGNMENT CEILING. Does the wet-lab molecular assay itself
     separate ClinVar path/benign? If yes, the phenotype is fitness-aligned and R2 is the right regime; if
     not, no molecular decoder can help (the regime question, answered per gene).
  2. **BLOSUM62 AUROC** — the deterministic, no-GPU, no-network decoder FLOOR (the same scorer the forward
     cell ships as its offline baseline).
  3. [AlphaMissense / ESM2+ProSST hybrid AUROC — the deployable learned decoder] — a follow-up: the committed
     am_filtered.tsv covers only ProteinGym-overlap variants, and the full-proteome AM / GPU hybrid run is the
     forward cell's established Kaggle pattern. Named, not faked.

HONEST RAILS:
  - LEAKAGE: these are canonical ProteinGym proteins, so the *hybrid* is NOT leakage-free on their DMS. BUT
    ClinVar path/benign labels are INDEPENDENT of the DMS-fitness tuning — the decoder was never fit to
    ClinVar — so the ClinVar-AUROC is a partly-independent clinical readout. Tier = in-distribution clinical,
    NOT held-out. (The leakage-free tier is the gene-agnostic MaveDB prospective holdout.)
  - PARTIAL CIRCULARITY: ClinGen uses DMS as PS3/BS3 functional evidence, so a few recent ClinVar calls may
    be DMS-informed. Most path/benign calls rest on segregation/population/clinical evidence, so the
    DMS-AUROC is largely (not perfectly) independent. Flagged, not hidden.
  - CLASS BALANCE is gene-specific and load-bearing: TP53/MSH2 have both classes; BRCA1 pathogenicity is
    truncating-dominated (missense-path rare) and PTEN missense-path dominate — single-class genes are
    reported as AUROC-INAPPLICABLE (a finding), never a fake number.
  - ORIENTATION is LABEL-FREE: DMS sign is not standardized across MaveDB assays, so we orient each DMS by
    its rank-correlation with BLOSUM62 (both -> higher=preserved) — an orientation that never consults the
    clinical labels, so it cannot inflate the AUROC toward them.

  uv run python scripts/clinical_variant_effect_validate.py            # real network (MaveDB + ClinVar E-utils)
  uv run python scripts/clinical_variant_effect_validate.py --gene TP53

Frozen AMR surface byte-unchanged (READ-only; imports blosum62_score + parse_hgvs_pro).
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from datetime import date as _date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.forward.variant_effect import blosum62_score  # noqa: E402
from scripts.mavedb_prospective_holdout import parse_hgvs_pro, proteingym_gene_symbols  # noqa: E402

MAVEDB = "https://api.mavedb.org/api/v1"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
CLINVAR_CACHE = Path("D:/dna_decode_cache/clinvar/eutils")

_AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G",
        "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
        "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"}
_CV_TITLE = re.compile(r"\(p\.([A-Za-z]{3})(\d+)([A-Za-z]{3})\)")

# Actionable human genes with a free MaveDB functional DMS. The chosen score set is the gene's canonical
# loss-of-function / activity / abundance assay (fitness-aligned with tumor-suppressor pathogenicity).
CLINICAL_GENES = {
    "TP53":  {"uniprot": "P04637", "urn": "urn:mavedb:00001234-a-1", "assay": "TP53 functional score (MaveDB 00001234)"},
    "MSH2":  {"uniprot": "P43246", "urn": "urn:mavedb:00000050-a-1", "assay": "MSH2 LOF (Jia 2020, MaveDB 00000050)"},
    "PTEN":  {"uniprot": "P60484", "urn": "urn:mavedb:00000013-a-1", "assay": "PTEN lipid-phosphatase (Mighell 2018, MaveDB 00000013)"},
    "BRCA1": {"uniprot": "P38398", "urn": "urn:mavedb:00000003-a-1", "assay": "BRCA1 RING SGE (Findlay 2018, MaveDB 00000003)"},
}
MIN_PER_CLASS = 15  # AUROC needs both classes; below this a gene is single-class-dominated (a finding)


# ------------------------------------------------------------------ pure helpers (offline-testable) --------

def auroc(labels_pos: list[bool], scores: list[float]) -> float:
    """AUROC that `scores` (higher = more likely POSITIVE) separates the positive class. Mann-Whitney U with
    mid-ranks (tie-correct). Pure; raises on a degenerate single-class set."""
    n = len(labels_pos)
    if n != len(scores):
        raise ValueError("labels/scores length mismatch")
    npos = sum(1 for x in labels_pos if x)
    nneg = n - npos
    if npos == 0 or nneg == 0:
        raise ValueError("AUROC undefined: only one class present")
    order = sorted(range(n), key=lambda i: scores[i])
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and scores[order[j + 1]] == scores[order[i]]:
            j += 1
        mid = (i + j) / 2.0 + 1.0  # 1-based mid-rank
        for k in range(i, j + 1):
            ranks[order[k]] = mid
        i = j + 1
    sum_pos = sum(ranks[i] for i in range(n) if labels_pos[i])
    u = sum_pos - npos * (npos + 1) / 2.0
    return u / (npos * nneg)


def _spearman_sign(a: list[float], b: list[float]) -> float:
    """Sign of the rank-correlation between a and b (+1/-1/0). Label-free DMS orientation anchor."""
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(v):
            j = i
            while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]:
                j += 1
            mid = (i + j) / 2.0
            for k in range(i, j + 1):
                r[order[k]] = mid
            i = j + 1
        return r
    ra, rb = rank(a), rank(b)
    ma, mb = sum(ra) / len(ra), sum(rb) / len(rb)
    cov = sum((x - ma) * (y - mb) for x, y in zip(ra, rb))
    va = sum((x - ma) ** 2 for x in ra) ** 0.5
    vb = sum((y - mb) ** 2 for y in rb) ** 0.5
    if va == 0 or vb == 0:
        return 0.0
    r = cov / (va * vb)
    return 1.0 if r > 0 else (-1.0 if r < 0 else 0.0)


# ------------------------------------------------------------------ real-surface fetchers ------------------

def fetch_mavedb_dms(urn: str) -> dict[tuple[str, int, str], float]:
    """{(wt,pos,alt): functional_score} for single-missense variants of a MaveDB score set (UniProt numbering)."""
    with urllib.request.urlopen(f"{MAVEDB}/score-sets/{urn}/scores", timeout=120) as r:
        text = r.read().decode("utf-8", "replace")
    out: dict[tuple[str, int, str], float] = {}
    for row in csv.DictReader(io.StringIO(text)):
        parsed = parse_hgvs_pro(row.get("hgvs_pro", ""))
        raw = row.get("score", "")
        if parsed is None or raw in ("", "NA", None):
            continue
        try:
            out[parsed] = float(raw)
        except ValueError:
            continue
    return out


def fetch_clinvar_missense(gene: str, use_cache: bool = True) -> dict[tuple[str, int, str], str]:
    """{(wt,pos,alt): 'PATH'|'BENIGN'} germline missense classifications for one gene (ClinVar E-utilities).
    Conflicting-classification records are EXCLUDED. Cached per-gene to D: so reruns are offline."""
    cache = CLINVAR_CACHE / f"{gene}.json"
    if use_cache and cache.exists():
        raw = json.loads(cache.read_text(encoding="utf-8"))
        return {tuple(json.loads(k)): v for k, v in raw.items()}

    def esearch(term: str) -> list[str]:
        q = urllib.parse.urlencode({"db": "clinvar", "term": term, "retmax": 500, "retmode": "json"})
        with urllib.request.urlopen(f"{EUTILS}/esearch.fcgi?{q}", timeout=60) as r:
            return json.loads(r.read().decode())["esearchresult"].get("idlist", [])

    def esummary(ids: list[str]) -> dict:
        q = urllib.parse.urlencode({"db": "clinvar", "id": ",".join(ids), "retmode": "json"})
        with urllib.request.urlopen(f"{EUTILS}/esummary.fcgi?{q}", timeout=60) as r:
            return json.loads(r.read().decode())["result"]

    labels: dict[tuple[str, int, str], str] = {}
    for sig, tag in (("pathogenic", "PATH"), ("benign", "BENIGN")):
        term = (f'{gene}[gene] AND "missense variant"[molecular consequence] '
                f'AND "{sig}"[clinical significance]')
        ids = esearch(term)
        for i in range(0, min(len(ids), 500), 100):
            for uid, rec in esummary(ids[i:i + 100]).items():
                if uid == "uids":
                    continue
                m = _CV_TITLE.search(rec.get("title", ""))
                desc = ((rec.get("germline_classification", {}) or {}).get("description") or "").lower()
                if not m or sig not in desc or "conflict" in desc:
                    continue
                wt3, pos, alt3 = m.group(1), int(m.group(2)), m.group(3)
                if wt3 in _AA3 and alt3 in _AA3:
                    labels[(_AA3[wt3], pos, _AA3[alt3])] = tag
            time.sleep(0.34)  # NCBI 3 req/s courtesy limit
    if use_cache:
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_text(json.dumps({json.dumps(list(k)): v for k, v in labels.items()}), encoding="utf-8")
    return labels


# ------------------------------------------------------------------ per-gene validation --------------------

def score_from_tables(gene: str, meta: dict, dms: dict[tuple[str, int, str], float],
                      clin: dict[tuple[str, int, str], str], pg_syms: set[str]) -> dict:
    """PURE (no network): join a DMS table + a ClinVar label table -> the per-gene validation record.
    Extracted from validate_gene so the scoring/orientation/gating logic is offline-testable."""
    shared = sorted(set(dms) & set(clin))
    n_path = sum(1 for k in shared if clin[k] == "PATH")
    n_benign = len(shared) - n_path
    rec = {
        "gene": gene, "uniprot": meta["uniprot"], "urn": meta["urn"], "assay": meta["assay"],
        "n_dms_missense": len(dms), "n_clinvar_labels": len(clin), "n_joined": len(shared),
        "n_path": n_path, "n_benign": n_benign,
        "in_proteingym": gene.upper() in {s.upper() for s in pg_syms},
        "auroc_applicable": n_path >= MIN_PER_CLASS and n_benign >= MIN_PER_CLASS,
    }
    if not rec["auroc_applicable"]:
        rec["note"] = (f"single-class-dominated ({n_path} path / {n_benign} benign) — AUROC inapplicable; "
                       f"this gene's clinical-missense signal is one-sided (a finding, not a failure)")
        return rec
    blos = [blosum62_score(k[0], k[2]) for k in shared]          # higher = conservative = preserved/benign
    dvals = [dms[k] for k in shared]
    labels_path = [clin[k] == "PATH" for k in shared]
    # LABEL-FREE orientation: make DMS agree with BLOSUM (both higher=preserved). Never consults labels.
    dsign = _spearman_sign(dvals, blos) or 1.0
    dms_pres = [dsign * v for v in dvals]
    # pathogenic predictor = -preserved (higher preserved -> less pathogenic)
    rec["dms_auroc"] = round(auroc(labels_path, [-v for v in dms_pres]), 4)
    rec["blosum_auroc"] = round(auroc(labels_path, [-v for v in blos]), 4)
    rec["dms_orientation"] = "as-is" if dsign > 0 else "flipped-to-agree-with-blosum"
    return rec


def validate_gene(gene: str, meta: dict, pg_syms: set[str], use_cache: bool = True) -> dict:
    dms = fetch_mavedb_dms(meta["urn"])
    clin = fetch_clinvar_missense(gene, use_cache=use_cache)
    return score_from_tables(gene, meta, dms, clin, pg_syms)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gene", help="validate a single gene (default: all in CLINICAL_GENES)")
    ap.add_argument("--no-cache", action="store_true", help="force fresh ClinVar E-utilities fetch")
    a = ap.parse_args()

    pg = proteingym_gene_symbols()
    genes = {a.gene: CLINICAL_GENES[a.gene]} if a.gene else CLINICAL_GENES
    results = []
    for g, meta in genes.items():
        print(f"[{g}] fetching MaveDB {meta['urn']} + ClinVar ...", flush=True)
        rec = validate_gene(g, meta, pg, use_cache=not a.no_cache)
        results.append(rec)
        if rec["auroc_applicable"]:
            print(f"  joined={rec['n_joined']} (path={rec['n_path']}/benign={rec['n_benign']})  "
                  f"DMS-AUROC={rec['dms_auroc']}  BLOSUM-AUROC={rec['blosum_auroc']}  "
                  f"{'[in ProteinGym]' if rec['in_proteingym'] else ''}")
        else:
            print(f"  {rec['note']}")

    scored = [r for r in results if r.get("auroc_applicable")]
    art = {
        "_schema": "clinical-variant-effect-validation-v1", "date": _date.today().isoformat(),
        "question": "Does the R2 forward variant-effect decoder separate ClinVar pathogenic vs benign "
                    "missense on actionable human genes?",
        "tier": "in_distribution_clinical (proteins in ProteinGym; ClinVar labels independent of DMS tuning)",
        "ceiling_metric": "DMS-functional AUROC = fitness-alignment ceiling (can the molecular assay itself "
                          "separate clinical labels?)",
        "floor_metric": "BLOSUM62 AUROC = deterministic no-GPU decoder floor",
        "deployable_learned_followup": "AlphaMissense (full-proteome) + ESM2+ProSST hybrid AUROC — Kaggle/full-AM "
                                       "follow-up; committed am_filtered.tsv covers only ProteinGym-overlap variants",
        "circularity_caveat": "ClinGen may use DMS as PS3/BS3 evidence for a minority of recent calls; most "
                              "path/benign rest on segregation/population evidence -> DMS-AUROC largely independent",
        "orientation": "DMS oriented label-free via rank-agreement with BLOSUM62 (never consults clinical labels)",
        "min_per_class": MIN_PER_CLASS, "n_genes": len(results), "n_auroc_scored": len(scored),
        "results": results, "frozen_surface_changed": False,
    }
    out = Path(f"wiki/clinical_variant_effect_validation_{_date.today().isoformat()}.json")
    out.write_text(json.dumps(art, indent=2), encoding="utf-8")
    print(f"\nartifact: {out}  ({len(scored)}/{len(results)} genes AUROC-scored)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
