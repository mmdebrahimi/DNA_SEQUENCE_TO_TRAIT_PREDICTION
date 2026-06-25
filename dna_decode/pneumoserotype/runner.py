"""Streptococcus pneumoniae capsular serotype caller — cps-reference DB via the shared blastn engine.

The pneumococcal capsule (cps) locus determines serotype; PneumoCaT / SeroBA / the GPS pipeline type a
genome by matching its cps locus against a curated per-serotype reference set (~90+ serotypes). This v0
mirrors that as the "smallest credible slice": blastn a curated cps-reference DB (one representative
sequence per serotype) vs the assembly and call the serotype of the best-matching reference (max coverage
among identity-passing hits). Same shape as dna-ktype (best allele -> type) and dna-serotype.

DB layout: a single FASTA `cps_references.fasta`, headers `serotype__<SEROTYPE>__<id>`, e.g.
`serotype__19F__01`, `serotype__6B__01`.

HONESTY (load-bearing): faithful to the cps-reference typing method (blastn best-match); NOT an independent
baseline. A single-best-reference v0 resolves SEROGROUP reliably but within-serogroup pairs that differ by
a single locus/SNP (e.g. 6A vs 6B at wciP, 19A vs 19F) need allele-level logic the full tools (PneumoCaT
CTVdb / SeroBA) add -> reported as the best-match serotype with that documented ceiling (the GPS pipeline
itself reports ~89% in-silico-vs-Quellung concordance). S. pneumoniae only. Offline-safe via the shared
engine (missing blastn/DB -> status 'unavailable').
"""
from __future__ import annotations

from pathlib import Path

from dna_decode.typing.blast_caller import call_alleles

# cps references are long, near-exact within a serotype; require solid coverage of the reference.
PNEUMO_IDENTITY_THRESHOLD = 90.0
PNEUMO_COVERAGE_THRESHOLD = 70.0


def serotype_of(ref_id: str) -> str | None:
    """`serotype__19F__01` -> '19F'. None if no serotype token."""
    parts = ref_id.split("__")
    return parts[1] if (len(parts) >= 2 and parts[0] == "serotype" and parts[1]) else None


def call_pneumo_serotype(fasta: str | Path, db_dir: str | Path, *,
                         identity_threshold: float = PNEUMO_IDENTITY_THRESHOLD,
                         coverage_threshold: float = PNEUMO_COVERAGE_THRESHOLD,
                         blastn_bin: str | None = None, timeout: int = 600) -> dict:
    """blastn the cps-reference DB vs `fasta`; call the serotype of the best-matching reference."""
    db_dir = Path(db_dir)
    cps_fasta = db_dir / "cps_references.fasta"
    if not cps_fasta.exists():
        return {"status": "unavailable", "serotype": None,
                "reason": f"cps reference DB not found in {db_dir} (need cps_references.fasta)"}
    res = call_alleles(fasta, cps_fasta, identity_threshold=identity_threshold,
                       coverage_threshold=coverage_threshold, blastn_bin=blastn_bin, timeout=timeout)
    if res["status"] != "ok":
        return {"status": "unavailable", "tool": res.get("tool"), "serotype": None,
                "reason": res.get("reason")}

    best = None  # (ref_id, (coverage, identity), hit) among CALLED refs
    for ref_id, hit in res["per_allele"].items():
        if not hit["called"]:
            continue
        key = (hit["percent_coverage"], hit["percent_identity"])
        if best is None or key > best[1]:
            best = (ref_id, key, hit)

    base = {"status": "ok", "tool": "blastn", "method": "cps_reference_blastn_v0",
            "parameters": {"identity_threshold": identity_threshold, "coverage_threshold": coverage_threshold}}
    if best is None:
        return {**base, "serotype": None, "best_reference": None,
                "note": "no cps reference matched at threshold (novel/non-typeable/partial cps locus)"}
    ref_id, _, hit = best
    st = serotype_of(ref_id)
    return {**base, "serotype": st, "best_reference": ref_id,
            "percent_identity": hit["percent_identity"], "percent_coverage": hit["percent_coverage"]}
