"""Expression-context signal: IS-element-upstream-of-target-gene detection from an assembly.

The deterministic AMR decoder predicts R/S from PRESENCE of curated determinants. That is blind to
EXPRESSION-driven resistance — where an intrinsic/weak gene is OVEREXPRESSED by an upstream insertion
sequence (IS) providing a hybrid promoter. The canonical case is *Acinetobacter baumannii* x meropenem:
ISAba1 inserted immediately upstream of the intrinsic blaOXA-51 gene drives carbapenem resistance, yet
every isolate CARRIES OXA-51, so gene-presence cannot tell R from S and the decoder ABSTAINs.

This module reads the SAME assembly the decoder already has and detects the IS-upstream-of-target junction
deterministically (blastn the IS ref + the target-gene ref vs the assembly; an IS hit within `upstream_bp`
5-prime of a target hit ON THE SAME CONTIG = a junction). It is the consumer-side signal that can convert a
known ABSTAIN into a correct R for the EXPRESSION_FLOOR organism x drug — gated, validated independently
first (see plans/Expression_Context_Acinetobacter_Meropenem_Plan/).

PRIMARY rule = the EXACT frozen falsifier rule (`soraya_runs/2026-06-10-7qjq-expression-frontier-isaba1/`):
same-contig + strand-aware UPSTREAM-of-target proximity, NO IS-orientation constraint. That rule produced
junction-positive R 1/15, S 0/15 on the N=30 cohort; porting it verbatim keeps the validated signal
identical to the falsified one. `is_orientation=True` is an OPTIONAL refinement that adds an IS-orientation
constraint — it MUST be re-falsified on the original 30 (must still rescue GCA_000692095.1) before use, so
it is OFF by default and never enabled in the primary path.

Reuses the shared blastn resolvers (`vf_runner.find_blastn`/`_find_makeblastdb`); offline-safe — missing
blastn/makeblastdb/ref returns `{status: "unavailable", signal: False}` and never raises (so callers + tests
stay green without BLAST+). Unlike `typing/blast_caller.call_alleles`, this enumerates ALL hits per query
(NO `-max_target_seqs` truncation) because multi-copy IS elements matter (the rescued strain has 22 ISAba1
copies — a best-hit-only caller would miss the upstream copy).
"""
from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

from dna_decode.pathotype.vf_runner import _find_makeblastdb, find_blastn

# Frozen thresholds (from the 2026-06-10 falsifier — do NOT re-tune the primary path).
DEFAULT_UPSTREAM_BP = 400
DEFAULT_MIN_IDENTITY = 85.0
DEFAULT_MIN_LEN = 120
# No truncation: multi-copy IS elements (ISAba1 can be >20 copies) must all be enumerated.
_NO_TRUNCATION = "100000"


def _blast_all_hits(query: str, asm_db: str, blastn: str, *, min_identity: float, min_len: int,
                    timeout: int) -> list[dict]:
    """blastn `query` vs the assembly DB; return ALL identity/length-passing hits with contig + position.

    outfmt includes `sseqid` (the contig — REQUIRED for the same-contig test) + subject start/end (strand
    inferred from their order). No `-max_target_seqs` cap so multi-copy hits are never truncated.
    """
    proc = subprocess.run(
        [blastn, "-query", query, "-db", asm_db,
         "-outfmt", "6 qseqid sseqid sstart send pident length",
         "-perc_identity", str(min_identity), "-max_target_seqs", _NO_TRUNCATION, "-evalue", "1e-20"],
        check=True, capture_output=True, text=True, timeout=timeout)
    hits = []
    for line in proc.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 6:
            continue
        try:
            sstart, send, pident, length = int(f[2]), int(f[3]), float(f[4]), int(f[5])
        except ValueError:
            continue
        if pident < min_identity or length < min_len:
            continue
        hits.append({"contig": f[1], "lo": min(sstart, send), "hi": max(sstart, send),
                     "strand": "+" if sstart < send else "-", "pident": pident, "length": length})
    return hits


def detect_is_upstream_junction(genome_fasta: str | Path, *, is_ref: str | Path, target_ref: str | Path,
                                upstream_bp: int = DEFAULT_UPSTREAM_BP,
                                min_identity: float = DEFAULT_MIN_IDENTITY,
                                min_len: int = DEFAULT_MIN_LEN, is_orientation: bool = False,
                                blastn_bin: str | None = None, timeout: int = 600) -> dict:
    """Detect an IS-element inserted upstream of a target gene in `genome_fasta`.

    PRIMARY rule (`is_orientation=False`, the validated default): signal=True iff an IS hit lies within
    `upstream_bp` of the 5-prime (upstream) boundary of a target hit ON THE SAME CONTIG, strand-aware on the
    target. NO IS-orientation constraint (matches the frozen falsifier).

    `is_orientation=True` (OPTIONAL refinement, NOT validated): additionally require the IS element to be on
    the same strand as the target (a coarse promoter-direction proxy). Must be re-falsified before use.

    Returns {status: "ok"|"unavailable", signal: bool, evidence: {...}} — never raises.
    """
    blastn = blastn_bin or find_blastn()
    if not blastn:
        return {"status": "unavailable", "signal": False,
                "reason": "blastn not found (set $BLASTN_BIN or install NCBI BLAST+)", "evidence": {}}
    makeblastdb = _find_makeblastdb(blastn)
    if not makeblastdb:
        return {"status": "unavailable", "signal": False,
                "reason": "makeblastdb not found alongside blastn", "evidence": {}}
    for ref, label in ((is_ref, "is_ref"), (target_ref, "target_ref"), (genome_fasta, "genome_fasta")):
        if not Path(ref).exists():
            return {"status": "unavailable", "signal": False,
                    "reason": f"{label} not found at {ref}", "evidence": {}}

    with tempfile.TemporaryDirectory(prefix="exprctx_") as td:
        asm_db = str(Path(td) / "asm")
        try:
            subprocess.run([makeblastdb, "-in", str(genome_fasta), "-dbtype", "nucl", "-out", asm_db],
                           check=True, capture_output=True, text=True, timeout=timeout)
            is_hits = _blast_all_hits(str(is_ref), asm_db, blastn, min_identity=min_identity,
                                      min_len=min_len, timeout=timeout)
            tgt_hits = _blast_all_hits(str(target_ref), asm_db, blastn, min_identity=min_identity,
                                       min_len=min_len, timeout=timeout)
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired, OSError) as e:
            return {"status": "unavailable", "signal": False,
                    "reason": f"blastn invocation failed: {type(e).__name__}", "evidence": {}}

    junction = None
    for t in tgt_hits:
        # 5-prime (upstream) window of the target gene, strand-aware.
        if t["strand"] == "+":
            win_lo, win_hi = t["lo"] - upstream_bp, t["lo"]
        else:
            win_lo, win_hi = t["hi"], t["hi"] + upstream_bp
        for i in is_hits:
            if i["contig"] != t["contig"]:
                continue
            if i["hi"] < win_lo or i["lo"] > win_hi:        # IS hit overlaps the upstream window?
                continue
            if is_orientation and i["strand"] != t["strand"]:   # refinement (default off)
                continue
            distance = win_hi - i["hi"] if t["strand"] == "+" else i["lo"] - win_lo
            junction = {"contig": t["contig"], "distance_bp": max(0, distance),
                        "is_strand": i["strand"], "target_strand": t["strand"]}
            break
        if junction:
            break

    return {
        "status": "ok",
        "signal": junction is not None,
        "evidence": {
            "n_is_hits": len(is_hits),
            "n_target_hits": len(tgt_hits),
            "junction": junction,
            "upstream_bp": upstream_bp,
            "is_orientation_enforced": is_orientation,
            "raw_hits": {"is": is_hits, "target": tgt_hits},   # retained for reproducibility / off-target audit
        },
    }
