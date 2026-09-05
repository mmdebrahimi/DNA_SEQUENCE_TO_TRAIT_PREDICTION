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


LOCUS_OVERLAP_FRACTION = 0.50   # two hits are the same locus at >=50% overlap of the SHORTER interval


def _interval(hit: dict) -> tuple[str, int, int] | None:
    """(contig, start, end) with start<=end. A minus-strand HSP arrives as sstart>send."""
    if hit.get("contig") is None or hit.get("sstart") is None or hit.get("send") is None:
        return None
    s, e = int(hit["sstart"]), int(hit["send"])
    return (str(hit["contig"]), min(s, e), max(s, e))


def _overlap_fraction(a: tuple[str, int, int], b: tuple[str, int, int]) -> float:
    """Overlap as a fraction of the SHORTER interval; 0.0 across different contigs."""
    if a[0] != b[0]:
        return 0.0
    inter = min(a[2], b[2]) - max(a[1], b[1]) + 1
    if inter <= 0:
        return 0.0
    return inter / min(a[2] - a[1] + 1, b[2] - b[1] + 1)


def cluster_alleles_by_locus(called: list[tuple[str, dict]]) -> list[list[tuple[str, dict]]]:
    """Group called alleles that hit the SAME genomic locus. GREEDY-REPRESENTATIVE, not single-linkage.

    WHY THIS FUNCTION EXISTS. beta-lactamase variants differ by one to three point mutations, so a single
    blaTEM locus matches ~180 catalog TEM alleles above the 90% identity bar. Keying the output on the
    ALLELE name reports every one of them as a separately present gene -- including ESBLs like blaTEM-52
    in a genome that carries only the narrow-spectrum blaTEM-1.

    WHY NOT GROUP BY GENE NAME. blaOXA-1 and blaOXA-48 share a name prefix and are functionally unrelated
    (narrow-spectrum vs carbapenemase), and a genome can genuinely carry both. Position is the honest
    criterion: same locus -> one gene; different loci -> genuinely two.

    WHY GREEDY-REPRESENTATIVE. Single-linkage CHAINS -- an adjacent tandem array (this project has seen a
    real 7-copy blaTEM array) would merge into one call through a chain of pairwise overlaps. Each hit is
    instead compared against the cluster's REPRESENTATIVE, which is the highest-identity member, so a
    tandem copy that does not overlap the representative starts its own locus. Same method as
    `eval/clonality.greedy_representative_clusters_from_matrix`, for the same reason.

    An allele with NO position (a caller that did not return coordinates) becomes its own singleton
    cluster: what cannot be placed must not be silently merged.
    """
    ordered = sorted(called, key=lambda kv: (-kv[1]["percent_identity"], -kv[1]["percent_coverage"],
                                             kv[0]))
    clusters: list[list[tuple[str, dict]]] = []
    reps: list[tuple[str, int, int] | None] = []
    for allele_id, hit in ordered:
        iv = _interval(hit)
        placed = False
        if iv is not None:
            for i, rep in enumerate(reps):
                if rep is not None and _overlap_fraction(iv, rep) >= LOCUS_OVERLAP_FRACTION:
                    clusters[i].append((allele_id, hit))
                    placed = True
                    break
        if not placed:
            clusters.append([(allele_id, hit)])
            reps.append(iv)
    return clusters


def call_resistance_genes(fasta: str | Path, db: str | Path, *, drug_class: str | None = None,
                          identity_threshold: float = RES_IDENTITY_THRESHOLD,
                          coverage_threshold: float = RES_COVERAGE_THRESHOLD,
                          blastn_bin: str | None = None, timeout: int = 600) -> dict:
    """blastn a ResFinder class allele DB vs `fasta`; return the called acquired-AMR genes, one per LOCUS.

    `drug_class` (e.g. 'beta-lactam') labels every gene from this DB file; pass it when running per-class DBs.

    One entry per genomic locus, naming the best-matching allele there. Selection is IDENTITY-PRIMARY:
    within a locus the variants are all at ~100% coverage, so a coverage-first tiebreak is decided by
    dict order rather than by sequence -- the same defect fixed in the serotype and salmserovar callers.
    `n_alleles_at_locus` reports how many catalog alleles matched, so the collapse stays auditable.
    """
    res = call_alleles(fasta, db, identity_threshold=identity_threshold,
                       coverage_threshold=coverage_threshold, blastn_bin=blastn_bin, timeout=timeout,
                       with_positions=True)
    if res["status"] != "ok":
        return {"status": "unavailable", "tool": res.get("tool"), "genes": [], "reason": res.get("reason")}

    called = [(aid, h) for aid, h in res["per_allele"].items() if h["called"]]
    genes = []
    for cluster in cluster_alleles_by_locus(called):
        allele_id, hit = max(cluster, key=lambda kv: (kv[1]["percent_identity"],
                                                      kv[1]["percent_coverage"]))
        genes.append({"gene": gene_of(allele_id), "drug_class": drug_class, "best_allele": allele_id,
                      "percent_identity": hit["percent_identity"],
                      "percent_coverage": hit["percent_coverage"],
                      "n_alleles_at_locus": len(cluster),
                      "contig": hit.get("contig"), "sstart": hit.get("sstart"), "send": hit.get("send")})
    genes.sort(key=lambda r: (-r["percent_identity"], -r["percent_coverage"], r["gene"]))
    return {
        "status": "ok", "tool": "blastn", "method": "resfinder_blastn_v1_locus_collapsed",
        "parameters": {"identity_threshold": identity_threshold, "coverage_threshold": coverage_threshold,
                       "drug_class": drug_class, "locus_overlap_fraction": LOCUS_OVERLAP_FRACTION},
        "n_alleles_called": len(called),
        "genes": genes,
    }
