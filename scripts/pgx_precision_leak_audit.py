"""PGx precision-leak audit (READ-ONLY) — quantify the silent non-core -> *1 mis-call exposure.

Turns the "no sentinel layer (v0)" caveat into a NUMBER DNA-11 can prioritize against. For each PGx gene
with a COMMITTED GeT-RM truth set, count the real 1000G-overlap samples whose GeT-RM consensus diplotype
carries a NON-CORE star-allele the v0 core-SNP proxy cannot resolve -> today silently mis-called *1 (unless
a sentinel guards it). This is a property of the TRUTH set + the caller's core allele set, so it is computed
OFFLINE from committed data (NO 1000G VCF / NO network / NO Docker) and touches NO pgx code.

Faithful to the caller's own core definition: the per-gene `core` sets are IMPORTED from
`scripts.pgx_getrm_concordance.GENES` (single source of truth); the tier rule mirrors that script
(non-core = a star whose primary token is outside `core`; CYP2D6 structural/ambiguous truth excluded).

GUARDED genes (CYP2C19/CYP2C9 populate a SENTINELS list) -> their non-core samples are correctly WITHHELD
(a control: the sentinel doing its job). LEAK genes (SENTINELS=[]) -> their non-core count IS the exposure.

Read-only; writes only wiki/pgx_precision_leak_audit_<date>.{md,json} (Soraya-owned artifacts).
"""
from __future__ import annotations

import csv
import json
import re
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.pgx_getrm_concordance import GENES, TRUTH  # core sets + truth-file locations (single source)

# A gene is GUARDED iff its catalog ships a non-empty SENTINELS list (derived live, not hardcoded, so the
# audit tracks population as it lands -- e.g. TPMT moved leak->guarded once its 10 sentinels were populated).
def _guarded(gene: str) -> bool:
    import importlib
    try:
        return bool(getattr(importlib.import_module(f"dna_decode.pgx.{gene}_catalog"), "SENTINELS", []))
    except Exception:
        return False


GUARDED = {g for g in ("cyp2c19", "cyp2c9", "cyp2c8", "cyp3a5", "tpmt", "cyp2b6", "cyp2d6") if _guarded(g)}
# CYP2D6 structural stars (deletion/dup/hybrid) — NOT SNP-decodable; excluded from the SNP exposure count.
_CYP2D6_STRUCTURAL = {"*5", "*13", "*36", "*61", "*63", "*68"}


def _stars(diplo: str) -> list[str]:
    """Split a diplotype on / | + into raw star tokens (keeps suballele letters, e.g. '*3A', '*2A')."""
    out = []
    for tok in re.split(r"[/|+()]", diplo):
        m = re.match(r"\s*(\*\w+)", tok)
        if m:
            out.append(m.group(1))
    return out


def _primary(star: str) -> str:
    """'*3A' -> '*3'; '*17' -> '*17'."""
    d = "".join(ch for ch in star.lstrip("*") if ch.isdigit())
    return f"*{d}" if d else star


def _in_core(star: str, core: set[str]) -> bool:
    """Faithful membership: exact token in core (TPMT *3A/*3B/*3C) OR its primary in core (*2A -> *2)."""
    return star in core or _primary(star) in core


def audit_gene(key: str, cfg: dict) -> dict:
    core = cfg["core"]
    equiv = cfg.get("ref_equiv", {})   # reference-equivalent aliases (e.g. CYP2C19 *38 == *1) -> treat as core
    truth_file = cfg.get("truth_file", TRUTH)
    col = cfg["truth_col"]
    is_2d6 = key == "cyp2d6"
    n_core = n_noncore = n_structural = n_ambiguous = 0
    noncore_alleles: Counter = Counter()
    if not Path(truth_file).exists():
        return {"gene": key, "status": "no_committed_truth", "core_set": sorted(core)}
    with open(truth_file, encoding="utf-8") as fh:
        for rec in csv.DictReader(fh, delimiter="\t"):
            raw = (rec.get(col, "") or "").strip()
            if not raw or raw in ("NA", ".", ""):
                continue
            stars = _stars(raw)
            if not stars:
                continue
            if is_2d6 and (re.search(r"[xX]\d", raw) or any(s in _CYP2D6_STRUCTURAL for s in stars)):
                n_structural += 1
                continue
            if "(" in raw:            # parenthetical alternative -> genuinely ambiguous truth
                n_ambiguous += 1
                continue
            stars = [equiv.get(s, s) for s in stars]   # map reference-equivalent aliases to core (*38 -> *1)
            non = [s for s in stars if not _in_core(s, core)]
            if non:
                n_noncore += 1
                noncore_alleles.update(_primary(s) for s in non)
            else:
                n_core += 1
    scored = n_core + n_noncore
    return {
        "gene": key, "status": "scored", "guarded": key in GUARDED, "core_set": sorted(core),
        "n_scored_samples": scored, "n_core": n_core, "n_noncore": n_noncore,
        "n_structural_excluded": n_structural, "n_ambiguous_excluded": n_ambiguous,
        "noncore_rate": round(n_noncore / scored, 4) if scored else None,
        "noncore_alleles": dict(noncore_alleles.most_common()),
        "interpretation": ("non-core samples correctly WITHHELD by the sentinel layer (control)"
                           if key in GUARDED else
                           "non-core samples SILENTLY MIS-CALLED *1 today (the leak this exposure quantifies)"),
    }


def main() -> int:
    results = [audit_gene(k, cfg) for k, cfg in GENES.items()]
    scored = [r for r in results if r["status"] == "scored"]
    leak = [r for r in scored if not r["guarded"] and r["gene"] != "cyp2d6"]
    guarded = [r for r in scored if r["guarded"]]
    d2d6 = next((r for r in scored if r["gene"] == "cyp2d6"), None)

    leak_exposure = sum(r["n_noncore"] for r in leak)
    leak_scored = sum(r["n_scored_samples"] for r in leak)
    guarded_withheld = sum(r["n_noncore"] for r in guarded)

    summary = {
        "audit": "pgx_precision_leak", "date": "2026-07-28", "method": "offline_from_committed_getrm_truth",
        "headline_leak_exposure_samples": leak_exposure,
        "headline_leak_scored_samples": leak_scored,
        "headline_leak_noncore_rate": round(leak_exposure / leak_scored, 4) if leak_scored else None,
        "guarded_noncore_correctly_withheld": guarded_withheld,
        "leak_genes": [r["gene"] for r in leak],
        "genes_no_committed_truth": [r["gene"] for r in results if r["status"] == "no_committed_truth"],
        "per_gene": results,
        "caveats": [
            "exposure = GeT-RM-consensus samples carrying a non-core star the v0 core-SNP proxy silently "
            "calls *1; it is a LOWER BOUND on real-world leak (the ~87-sample GeT-RM 1000G overlap is small "
            "+ non-core alleles are rarer in other populations).",
            "NUDT15/UGT1A1/DPYD have NO committed GeT-RM truth (their 'external wall') -> not quantifiable "
            "offline here; their leak is real but must be sized after a v0.1 truth fetch.",
            "CYP2D6 non-core SNP alleles are reported separately; its STRUCTURAL alleles (*5/*13/... = "
            "n_structural_excluded) are a different surface handled by the CRAM structural stack.",
        ],
    }
    out_json = Path("wiki/pgx_precision_leak_audit_2026-07-28.json")
    out_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps({k: summary[k] for k in (
        "headline_leak_exposure_samples", "headline_leak_scored_samples", "headline_leak_noncore_rate",
        "guarded_noncore_correctly_withheld", "leak_genes", "genes_no_committed_truth")}, indent=2))
    print("\nper-gene:")
    for r in scored:
        tag = "GUARDED" if r["guarded"] else ("2D6" if r["gene"] == "cyp2d6" else "LEAK")
        print(f"  [{tag:7}] {r['gene']:8} non-core {r['n_noncore']}/{r['n_scored_samples']} "
              f"(rate {r['noncore_rate']}) alleles={r['noncore_alleles']}")
    print(f"\nwrote {out_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
