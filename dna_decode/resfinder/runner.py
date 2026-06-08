"""Acquired-AMR-gene caller — ResFinder allele DB via the shared blastn engine.

ResFinder allele headers: '<gene>_<allele#>_<accession>', e.g. 'blaNDM-19_1_MF370080', 'aac(6')-Ib_2_M23634'
(gene names carry hyphens/parens/primes). The gene is captured lazily up to the first _<digits>_ separator.
A gene is CALLED when its best allele clears thresholds (ResFinder defaults 90% identity / 60% coverage).
Per-class grouping comes from which DB file the allele lived in (caller passes class labels alongside).

This is an INDEPENDENT acquired-gene caller (different curated DB than AMRFinder) — use the per-gene calls as
a cross-tool concordance check against dna-amr, not as a redundant copy. Offline-safe via the engine.
"""
from __future__ import annotations

import re
from pathlib import Path

from dna_decode.typing.blast_caller import call_alleles

RES_IDENTITY_THRESHOLD = 90.0   # ResFinder default min percent identity
RES_COVERAGE_THRESHOLD = 60.0   # ResFinder default min coverage

# gene name = lazily up to the first _<digits>_ (allele-number separator); names keep hyphens/parens/primes.
_GENE_RE = re.compile(r"^(?P<gene>.+?)_\d+_")


def gene_of(allele_id: str) -> str:
    """'blaNDM-19_1_MF370080' -> 'blaNDM-19'; \"aac(6')-Ib_2_M23634\" -> \"aac(6')-Ib\"."""
    m = _GENE_RE.match(allele_id)
    return m.group("gene") if m else allele_id.split("_", 1)[0]


def call_resistance_genes(fasta: str | Path, db: str | Path, *, drug_class: str | None = None,
                          identity_threshold: float = RES_IDENTITY_THRESHOLD,
                          coverage_threshold: float = RES_COVERAGE_THRESHOLD,
                          blastn_bin: str | None = None, timeout: int = 600) -> dict:
    """blastn a ResFinder class allele DB vs `fasta`; return the called acquired-AMR genes (best allele each).

    `drug_class` (e.g. 'beta-lactam') labels every gene from this DB file; pass it when running per-class DBs.
    """
    res = call_alleles(fasta, db, identity_threshold=identity_threshold,
                       coverage_threshold=coverage_threshold, blastn_bin=blastn_bin, timeout=timeout)
    if res["status"] != "ok":
        return {"status": "unavailable", "tool": res.get("tool"), "genes": [], "reason": res.get("reason")}

    gene_best: dict[str, dict] = {}
    for allele_id, hit in res["per_allele"].items():
        if not hit["called"]:
            continue
        g = gene_of(allele_id)
        cov = hit["percent_coverage"]
        cur = gene_best.get(g)
        if cur is None or cov > cur["percent_coverage"]:
            gene_best[g] = {"gene": g, "drug_class": drug_class, "best_allele": allele_id,
                            "percent_identity": hit["percent_identity"], "percent_coverage": cov}
    genes = sorted(gene_best.values(), key=lambda r: (-r["percent_coverage"], r["gene"]))
    return {
        "status": "ok", "tool": "blastn", "method": "resfinder_blastn_v0",
        "parameters": {"identity_threshold": identity_threshold, "coverage_threshold": coverage_threshold,
                       "drug_class": drug_class},
        "genes": genes,
    }
