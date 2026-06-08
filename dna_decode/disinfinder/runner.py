"""Biocide/disinfectant resistance gene caller — DisinFinder allele DB via the shared blastn engine.

DisinFinder allele headers: '<gene>_<allele#>_<accession>' (e.g. 'qacA_1_AB566410', 'formA_1_X73835') —
the same CGE shape as ResFinder, so it reuses resfinder's `gene_of` parser. A gene is CALLED when its best
allele clears thresholds (90% identity / 60% coverage). Offline-safe via the engine.
"""
from __future__ import annotations

from pathlib import Path

from dna_decode.resfinder.runner import gene_of  # identical CGE <gene>_<allele#>_<acc> header parser
from dna_decode.typing.blast_caller import call_alleles

DISIN_IDENTITY_THRESHOLD = 90.0
DISIN_COVERAGE_THRESHOLD = 60.0


def call_disinfectant_genes(fasta: str | Path, db: str | Path, *,
                            identity_threshold: float = DISIN_IDENTITY_THRESHOLD,
                            coverage_threshold: float = DISIN_COVERAGE_THRESHOLD,
                            blastn_bin: str | None = None, timeout: int = 600) -> dict:
    """blastn the DisinFinder allele DB vs `fasta`; return the called biocide-resistance genes."""
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
            gene_best[g] = {"gene": g, "best_allele": allele_id,
                            "percent_identity": hit["percent_identity"], "percent_coverage": cov}
    genes = sorted(gene_best.values(), key=lambda r: (-r["percent_coverage"], r["gene"]))
    return {"status": "ok", "tool": "blastn", "method": "disinfinder_blastn_v0", "genes": genes}
