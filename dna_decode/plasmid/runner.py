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
import subprocess
import tempfile
from pathlib import Path

# reuse the pathotype BLAST resolvers — same install, no duplication
from dna_decode.pathotype.vf_runner import _find_makeblastdb, find_blastn

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
    blastn = blastn_bin or find_blastn()
    if not blastn:
        return {"status": "unavailable", "tool": "blastn", "replicons": [], "per_allele": {},
                "reason": "blastn not found (set $BLASTN_BIN or install NCBI BLAST+)"}
    makeblastdb = _find_makeblastdb(blastn)
    if not makeblastdb:
        return {"status": "unavailable", "tool": "makeblastdb", "replicons": [], "per_allele": {},
                "reason": "makeblastdb not found alongside blastn"}
    if not Path(db).exists():
        return {"status": "unavailable", "tool": "db", "replicons": [], "per_allele": {},
                "reason": f"PlasmidFinder DB not found at {db}"}

    with tempfile.TemporaryDirectory(prefix="plasmid_") as td:
        asm_db = str(Path(td) / "asm")
        try:
            subprocess.run([makeblastdb, "-in", str(fasta), "-dbtype", "nucl", "-out", asm_db],
                           check=True, capture_output=True, text=True, timeout=timeout)
            proc = subprocess.run(
                [blastn, "-query", str(db), "-db", asm_db,
                 "-outfmt", "6 qseqid pident length qlen", "-perc_identity", str(identity_threshold),
                 "-max_target_seqs", "5", "-evalue", "1e-20"],
                check=True, capture_output=True, text=True, timeout=timeout)
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired, OSError) as e:
            return {"status": "unavailable", "tool": "blastn", "replicons": [], "per_allele": {},
                    "reason": f"blastn invocation failed: {type(e).__name__}"}

    best: dict[str, tuple[float, float]] = {}          # allele_id -> (pident, coverage)
    for line in proc.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        qid, pident_s, length_s, qlen_s = parts[:4]
        try:
            pident = float(pident_s)
            cov = (float(length_s) / float(qlen_s)) * 100.0 if float(qlen_s) else 0.0
        except ValueError:
            continue
        prev = best.get(qid)
        if prev is None or cov > prev[1]:
            best[qid] = (pident, cov)

    per_allele: dict[str, dict] = {}
    rep_best: dict[str, dict] = {}
    for allele_id, (pident, cov) in best.items():
        called = pident >= identity_threshold and cov >= coverage_threshold
        rep = replicon_family(allele_id)
        per_allele[allele_id] = {"replicon": rep, "percent_identity": round(pident, 1),
                                 "percent_coverage": round(cov, 1), "called": called}
        if not called:
            continue
        cur = rep_best.get(rep)
        if cur is None or cov > cur["percent_coverage"]:
            rep_best[rep] = {"replicon": rep, "best_allele": allele_id,
                             "percent_identity": round(pident, 1), "percent_coverage": round(cov, 1)}
    replicons = sorted(rep_best.values(), key=lambda r: (-r["percent_coverage"], r["replicon"]))
    return {
        "status": "ok", "tool": "blastn", "method": "plasmidfinder_blastn_v0",
        "parameters": {"identity_threshold": identity_threshold, "coverage_threshold": coverage_threshold},
        "replicons": replicons,
        "per_allele": per_allele,
    }
