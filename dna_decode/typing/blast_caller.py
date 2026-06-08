"""Generic curated-DB allele caller — the shared core of the CGE-family typing decoders.

`call_alleles` BLASTs an allele DB (query) against an assembly (subject DB), returns the best HSP per
allele (max coverage among identity-passing hits) + a `called` flag vs the thresholds. Every CGE-style
decoder (plasmid replicon / serotype / resfinder gene presence) builds its aggregation on this — the only
per-decoder differences are the DB, the header→entity parser, and the report shape.

Reuses `pathotype/vf_runner`'s blastn resolvers (find_blastn / _find_makeblastdb) — DRY, one install path.
Offline-safe: missing blastn/makeblastdb/DB OR a failed invocation -> {status: "unavailable", reason}, never
raises (so CLIs + tests stay green without BLAST+).
"""
from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

from dna_decode.pathotype.vf_runner import _find_makeblastdb, find_blastn


def call_alleles(fasta: str | Path, db: str | Path, *, identity_threshold: float,
                 coverage_threshold: float, blastn_bin: str | None = None, timeout: int = 600) -> dict:
    """blastn the allele `db` vs `fasta`; return best hit per allele + `called` vs thresholds.

    Returns {status: "ok", per_allele: {allele_id: {percent_identity, percent_coverage, called}}}
    or {status: "unavailable", reason, per_allele: {}}.
    """
    blastn = blastn_bin or find_blastn()
    if not blastn:
        return {"status": "unavailable", "tool": "blastn", "per_allele": {},
                "reason": "blastn not found (set $BLASTN_BIN or install NCBI BLAST+)"}
    makeblastdb = _find_makeblastdb(blastn)
    if not makeblastdb:
        return {"status": "unavailable", "tool": "makeblastdb", "per_allele": {},
                "reason": "makeblastdb not found alongside blastn"}
    if not Path(db).exists():
        return {"status": "unavailable", "tool": "db", "per_allele": {},
                "reason": f"allele DB not found at {db}"}

    with tempfile.TemporaryDirectory(prefix="typing_") as td:
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
            return {"status": "unavailable", "tool": "blastn", "per_allele": {},
                    "reason": f"blastn invocation failed: {type(e).__name__}"}

    best: dict[str, tuple[float, float]] = {}
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

    per_allele = {aid: {"percent_identity": round(p, 1), "percent_coverage": round(c, 1),
                        "called": p >= identity_threshold and c >= coverage_threshold}
                  for aid, (p, c) in best.items()}
    return {"status": "ok", "tool": "blastn", "per_allele": per_allele}


def called_alleles(result: dict) -> dict[str, dict]:
    """Filter a call_alleles result to only the alleles that cleared both thresholds."""
    return {aid: h for aid, h in result.get("per_allele", {}).items() if h["called"]}
