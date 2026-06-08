"""Plasmid replicon caller — real blastn of the PlasmidFinder allele DB vs an assembly.

Mirrors `pathotype/vf_runner.run_canonical_vf` (and reuses its `find_blastn` / `_find_makeblastdb`
resolvers — DRY): query = replicon alleles, subject DB = the assembly, so coverage is of the reference
allele. A replicon is CALLED when its best allele clears identity + coverage thresholds (PlasmidFinder
defaults 95% / 60%). Aggregates per-allele hits up to the replicon name and a coarse Inc family.

Offline-safe: if blastn/makeblastdb are absent OR the invocation fails, returns
`{status: "unavailable", reason}` (never raises), so the CLI + tests stay green without BLAST+.
"""
from __future__ import annotations

import re
from pathlib import Path

# shared curated-DB blastn engine (DRY across plasmid / serotype / resfinder typing decoders)
from dna_decode.typing.blast_caller import call_alleles

PLASMID_IDENTITY_THRESHOLD = 95.0   # PlasmidFinder default min percent identity
PLASMID_COVERAGE_THRESHOLD = 60.0   # PlasmidFinder default min reference-allele coverage

# PlasmidFinder allele headers: "<replicon>_<allele#>__<acc>" or "<replicon>_<allele#>_<note>_<acc>".
# The replicon name itself can contain digits/parens/slashes/hyphens (IncHI1B(R27), IncB/O/K/Z, IncI1-I(Alpha)),
# so capture lazily up to the FIRST _<digits>_ (the allele-number separator).
_ALLELE_RE = re.compile(r"^(?P<rep>.+?)_\d+(?:__|_)")


def replicon_family(allele_id: str) -> str:
    """'IncFIA_1__AP001918' -> 'IncFIA'; 'IncB/O/K/Z_2__GU256641' -> 'IncB/O/K/Z'. Fallback: pre-'_' token."""
    m = _ALLELE_RE.match(allele_id)
    if m:
        return m.group("rep")
    return allele_id.split("_", 1)[0]


def call_replicons(fasta: str | Path, db: str | Path, *,
                   identity_threshold: float = PLASMID_IDENTITY_THRESHOLD,
                   coverage_threshold: float = PLASMID_COVERAGE_THRESHOLD,
                   blastn_bin: str | None = None, timeout: int = 600) -> dict:
    """blastn the PlasmidFinder allele DB vs `fasta`; report called replicons + per-allele best hits."""
    res = call_alleles(fasta, db, identity_threshold=identity_threshold,
                       coverage_threshold=coverage_threshold, blastn_bin=blastn_bin, timeout=timeout)
    if res["status"] != "ok":
        return {"status": "unavailable", "tool": res.get("tool"), "replicons": [], "per_allele": {},
                "reason": res.get("reason")}

    per_allele: dict[str, dict] = {}
    rep_best: dict[str, dict] = {}
    for allele_id, hit in res["per_allele"].items():
        rep = replicon_family(allele_id)
        per_allele[allele_id] = {"replicon": rep, **hit}
        if not hit["called"]:
            continue
        cov = hit["percent_coverage"]
        cur = rep_best.get(rep)
        if cur is None or cov > cur["percent_coverage"]:
            rep_best[rep] = {"replicon": rep, "best_allele": allele_id,
                             "percent_identity": hit["percent_identity"], "percent_coverage": cov}
    replicons = sorted(rep_best.values(), key=lambda r: (-r["percent_coverage"], r["replicon"]))
    return {
        "status": "ok", "tool": "blastn", "method": "plasmidfinder_blastn_v0",
        "parameters": {"identity_threshold": identity_threshold, "coverage_threshold": coverage_threshold},
        "replicons": replicons,
        "per_allele": per_allele,
    }
