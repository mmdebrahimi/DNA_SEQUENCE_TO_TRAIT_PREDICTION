"""E. coli serotype caller — SerotypeFinder allele DB via the shared blastn engine.

SerotypeFinder allele headers: '<gene>_<allele#>_<accession>_<antigen>', e.g. 'wzx_1_GU299791_O1',
'fliC_307_AY249994_H9'. The serotype antigen is the LAST '_'-token (O# / H#). The O antigen is
determined by any of wzx/wzy/wzm/wzt; the H antigen by fliC (+ a few others). The call = best-coverage
called O antigen + best-coverage called H antigen -> "O1:H9" (either may be missing -> reported as O?/H?).

PlasmidFinder/SerotypeFinder default thresholds: identity 85 / coverage 60 (serotype alleles are more
divergent than plasmid replicons; SerotypeFinder uses 85% identity). Offline-safe via the engine.
"""
from __future__ import annotations

from pathlib import Path

from dna_decode.typing.blast_caller import call_alleles

SERO_IDENTITY_THRESHOLD = 85.0    # SerotypeFinder default min percent identity
SERO_COVERAGE_THRESHOLD = 60.0    # SerotypeFinder default min coverage


def antigen_of(allele_id: str) -> str | None:
    """'wzx_1_GU299791_O1' -> 'O1'; 'fliC_307_AY249994_H9' -> 'H9'. None if no O/H suffix."""
    tok = allele_id.rsplit("_", 1)[-1]
    return tok if (tok[:1] in ("O", "H") and tok[1:2].isdigit()) else None


def gene_of(allele_id: str) -> str:
    """'wzx_1_GU299791_O1' -> 'wzx'."""
    return allele_id.split("_", 1)[0]


def call_serotype(fasta: str | Path, db: str | Path, *,
                  identity_threshold: float = SERO_IDENTITY_THRESHOLD,
                  coverage_threshold: float = SERO_COVERAGE_THRESHOLD,
                  blastn_bin: str | None = None, timeout: int = 600) -> dict:
    """blastn the SerotypeFinder allele DB vs `fasta`; return the O:H serotype + supporting antigen hits."""
    res = call_alleles(fasta, db, identity_threshold=identity_threshold,
                       coverage_threshold=coverage_threshold, blastn_bin=blastn_bin, timeout=timeout)
    if res["status"] != "ok":
        return {"status": "unavailable", "tool": res.get("tool"), "o_antigen": None, "h_antigen": None,
                "serotype": None, "antigens": [], "reason": res.get("reason")}

    # Best called allele per antigen, ranked by IDENTITY then coverage.
    #
    # Identity-PRIMARY is load-bearing for the H antigens, and this is the same defect the sibling
    # Salmonella caller already fixed: flagellin (fliC) alleles CROSS-HYBRIDIZE at near-full COVERAGE
    # across different antigen types, so a coverage-only tiebreak picks the WRONG antigen while looking
    # perfectly confident. The true antigen is the highest-IDENTITY hit.
    #
    # Measured here before the change (N=150 wet-lab-labelled isolates, 2026-09-04): the H axis
    # resolved 98.7% of the time but was only 77.0% ACCURATE, and its errors were concentrated rather
    # than diffuse -- H21->H8 alone was 9 of 34 misses (26%), the top four pairs 65%. That is the
    # signature of systematic allele confusion, not noise. The O axis, which uses wzx/wzy rather than
    # flagellin and cross-hybridizes far less, was 93.1% accurate.
    ag_best: dict[str, dict] = {}
    for allele_id, hit in res["per_allele"].items():
        if not hit["called"]:
            continue
        ag = antigen_of(allele_id)
        if ag is None:
            continue
        key = (hit["percent_identity"], hit["percent_coverage"])
        cur = ag_best.get(ag)
        if cur is None or key > (cur["percent_identity"], cur["percent_coverage"]):
            ag_best[ag] = {"antigen": ag, "gene": gene_of(allele_id), "best_allele": allele_id,
                           "percent_identity": hit["percent_identity"],
                           "percent_coverage": hit["percent_coverage"]}

    def _top(prefix):
        cands = [v for k, v in ag_best.items() if k.startswith(prefix)]
        return max(cands, key=lambda v: (v["percent_identity"], v["percent_coverage"])) if cands else None

    o, h = _top("O"), _top("H")
    o_ag = o["antigen"] if o else None
    h_ag = h["antigen"] if h else None
    serotype = f"{o_ag or 'O?'}:{h_ag or 'H?'}" if (o_ag or h_ag) else None
    return {
        "status": "ok", "tool": "blastn", "method": "serotypefinder_blastn_v0",
        "parameters": {"identity_threshold": identity_threshold, "coverage_threshold": coverage_threshold},
        "o_antigen": o_ag, "h_antigen": h_ag, "serotype": serotype,
        "antigens": sorted(ag_best.values(), key=lambda v: (v["antigen"][:1], -v["percent_coverage"])),
    }
