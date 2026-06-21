"""Canonical VirulenceFinder gene-caller (real blastn) + the v0.1 side-by-side diff.

This is the v0.1 add the ledger promised: a *canonical* gene-caller to compare the
fast k-mer-seed resolver (detect.py) against. It runs real `blastn` of the
VirulenceFinder E. coli allele DB against the input assembly — the same method
canonical VirulenceFinder uses (blastn over the allele DB with identity + coverage
thresholds), so it is faithful to "run via real blastn" without the heavier CGE
`virulencefinder` package + KMA plumbing.

HONESTY (ledger Decision T1; ratified interrogation Q2, 2026-06-04): BOTH callers
match the SAME VF DB (`virulence_ecoli.fsa`). So this diff is
lightweight-k-mer-seed-VF vs canonical-blastn-VF over ONE DB — a method-vs-method
AUDIT of the fast caller, NOT an independent baseline. `caller_is_independent_baseline`
stays False and the diff section carries the same caveat; we never emit a bare
headline "agreement %".

Offline-safe: when `blastn` is not installed the canonical call + the diff degrade to
`status: "unavailable"` with a reason (never silently dropped), so the CLI + its test
run on a host without the binary while still asserting the contract.
"""
from __future__ import annotations

import hashlib
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

from dna_decode.pathotype.markers import CLUSTER_MARKERS

# VirulenceFinder default thresholds (pinned into provenance).
VF_IDENTITY_THRESHOLD = 90.0   # min percent identity for a called allele hit
VF_COVERAGE_THRESHOLD = 60.0   # min percent of the reference allele covered by the hit
# resolver (k-mer) confident bar, for the per-cluster comparison framing.
RESOLVER_CONFIDENT_COV = 80.0

# Candidate locations searched when blastn is not on PATH (this host's native install
# lives here; harmless on other hosts — a missing dir just yields unavailable).
_KNOWN_BLAST_DIRS = [
    os.path.expanduser("~/ncbi-blast/bin"),
    "C:/Users/Farshad/ncbi-blast/bin",
]


def find_blastn() -> str | None:
    """Resolve a blastn executable: $BLASTN_BIN, then PATH, then known install dirs.
    Returns the path string, or None if unavailable (offline-safe trigger)."""
    env = os.environ.get("BLASTN_BIN")
    if env and Path(env).exists():
        return env
    onpath = shutil.which("blastn")
    if onpath:
        return onpath
    for d in _KNOWN_BLAST_DIRS:
        cand = Path(d) / ("blastn.exe" if os.name == "nt" else "blastn")
        if cand.exists():
            return str(cand)
    return None


def _find_makeblastdb(blastn_bin: str) -> str | None:
    """makeblastdb sits next to blastn in the same install; else fall back to PATH."""
    sib = Path(blastn_bin).with_name("makeblastdb.exe" if os.name == "nt" else "makeblastdb")
    if sib.exists():
        return str(sib)
    return shutil.which("makeblastdb")


def _cluster_for_allele(allele_id: str) -> str | None:
    """Map a VF allele header id to its marker cluster via the shared prefix catalog.
    Mirrors detect.build_vf_index so canonical + resolver use ONE cluster mapping."""
    gl = allele_id.split(":")[0].lower()
    for cluster, prefixes in CLUSTER_MARKERS.items():
        if any(gl.startswith(p) for p in prefixes):
            return cluster
    return None


def _db_sha256(db: str | Path) -> str | None:
    """sha256 (16-char prefix) of the VF DB file, or None if unreadable.

    Carried into every result so a downstream virulence tier / pathotype call can
    stamp exactly which curated DB produced its determinants (C2 provenance)."""
    try:
        return hashlib.sha256(Path(db).read_bytes()).hexdigest()[:16]
    except OSError:
        return None


def _interval_dedup(hsps: list[dict]) -> list[dict]:
    """Collapse overlapping HSPs of ONE copy; keep distinct copies at distinct coords (C3).

    `hsps` are normalized (start<=stop) HSPs for a single (allele_id, sseqid). Sorted by
    start; overlapping intervals merge (retaining the best-coverage HSP's identity /
    coverage / strand); disjoint intervals stay separate — a real tandem / multi-copy
    allele (e.g. a blaTEM array) yields one entry per copy rather than a single best-hit.
    """
    out: list[dict] = []
    for h in sorted(hsps, key=lambda x: (x["start"], x["stop"])):
        if out and h["start"] <= out[-1]["stop"]:
            last = out[-1]
            last["stop"] = max(last["stop"], h["stop"])
            if h["percent_coverage"] > last["percent_coverage"]:
                last["percent_identity"] = h["percent_identity"]
                last["percent_coverage"] = h["percent_coverage"]
                last["strand"] = h["strand"]
        else:
            out.append(dict(h))
    return out


def parse_blastn_outfmt6(stdout: str, identity_threshold: float,
                         coverage_threshold: float) -> dict:
    """Parse blastn `-outfmt 6 qseqid sseqid sstart send pident length qlen` into calls.

    PURE (no subprocess) so the genotype/coord/strand handling is unit-testable on
    synthetic stdout. Returns {per_gene, per_cluster, per_hit}:
      - per_gene / per_cluster: best-hit (max coverage) per cluster — byte-identical to
        the pre-change output so `build_vf_diff` is unaffected.
      - per_hit: every CALLED HSP (clustered OR not — C2), coords normalized
        (minus-strand sstart>send → start<=stop + strand), interval-deduped per copy.
    """
    best: dict[str, tuple[float, float]] = {}   # allele_id -> (pident, coverage)
    called_by_copy: dict[tuple[str, str], list[dict]] = {}
    for line in stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 7:
            continue
        qid, sseqid, sstart_s, send_s, pident_s, length_s, qlen_s = parts[:7]
        try:
            pident = float(pident_s)
            qlen = float(qlen_s)
            cov = (float(length_s) / qlen) * 100.0 if qlen else 0.0
            sstart, send = int(sstart_s), int(send_s)
        except ValueError:
            continue
        prev = best.get(qid)
        if prev is None or cov > prev[1]:
            best[qid] = (pident, cov)
        if pident >= identity_threshold and cov >= coverage_threshold:
            called_by_copy.setdefault((qid, sseqid), []).append({
                "start": min(sstart, send),
                "stop": max(sstart, send),
                "strand": "+" if sstart <= send else "-",
                "percent_identity": round(pident, 1),
                "percent_coverage": round(cov, 1),
            })

    per_gene: dict[str, dict] = {}
    per_cluster: dict[str, dict] = {}
    for allele_id, (pident, cov) in best.items():
        cluster = _cluster_for_allele(allele_id)
        if cluster is None:
            continue
        called = pident >= identity_threshold and cov >= coverage_threshold
        per_gene[allele_id] = {
            "cluster": cluster,
            "percent_identity": round(pident, 1),
            "percent_coverage": round(cov, 1),
            "called": called,
        }
        cur = per_cluster.get(cluster)
        # cluster best = the called allele with highest coverage (or, if none called, best seen).
        score = (called, cov)
        if cur is None or score > (cur["called"], cur["percent_coverage"]):
            per_cluster[cluster] = {
                "called": called, "best_gene": allele_id,
                "percent_identity": round(pident, 1), "percent_coverage": round(cov, 1),
            }
    # clusters with zero allele hits are absent in canonical too.
    for cluster in CLUSTER_MARKERS:
        per_cluster.setdefault(cluster, {"called": False, "best_gene": None,
                                         "percent_identity": None, "percent_coverage": None})

    # per_hit: every called allele (clustered OR not — C2), interval-deduped per copy.
    per_hit: list[dict] = []
    for (allele_id, sseqid), hsps in called_by_copy.items():
        cluster = _cluster_for_allele(allele_id)   # None when unclustered (kept, not dropped)
        vf_gene = allele_id.split(":")[0]
        for copy in _interval_dedup(hsps):
            per_hit.append({
                "allele_id": allele_id,
                "vf_gene": vf_gene,
                "cluster": cluster,
                "sseqid": sseqid,
                "start": copy["start"],
                "stop": copy["stop"],
                "strand": copy["strand"],
                "percent_identity": copy["percent_identity"],
                "percent_coverage": copy["percent_coverage"],
                "called": True,
            })
    return {"per_gene": per_gene, "per_cluster": per_cluster, "per_hit": per_hit}


def run_canonical_vf(fasta: str | Path, db: str | Path, *,
                     identity_threshold: float = VF_IDENTITY_THRESHOLD,
                     coverage_threshold: float = VF_COVERAGE_THRESHOLD,
                     blastn_bin: str | None = None,
                     all_hits: bool = False,
                     timeout: int = 600) -> dict:
    """Run canonical (blastn) VirulenceFinder gene-calling of `db` alleles vs `fasta`.

    Returns a dict with `status` "ok" | "unavailable", and on ok:
      - per-gene best hits (gene -> {cluster, percent_identity, percent_coverage, called})
      - per-cluster calls (cluster -> {called, best_gene, percent_identity, percent_coverage})
      - `per_hit`: ALL called HSPs (NOT cluster-filtered — C2: the DB has ~4942 alleles but
        only ~23 map to a resolver cluster; the genome-map virulence tier surfaces every
        called allele), each with retained coords {allele_id, vf_gene, cluster (None if
        unclustered), sseqid, start, stop, strand, percent_identity, percent_coverage,
        called}, interval-deduped per copy (tandem copies retained).
      - `db_sha`: the VF DB sha256 prefix (provenance).

    `all_hits=True` raises `-max_target_seqs` so tandem / multi-copy alleles are retained
    (the genome-map tier needs every copy); `per_gene` / `per_cluster` stay best-hit and
    byte-identical to the pre-change output (so `build_vf_diff` is UNCHANGED)."""
    db_sha = _db_sha256(db)
    blastn = blastn_bin or find_blastn()
    if not blastn:
        return {"status": "unavailable",
                "reason": "blastn not found (set $BLASTN_BIN or install NCBI BLAST+)",
                "tool": "blastn", "per_gene": {}, "per_cluster": {},
                "per_hit": [], "db_sha": db_sha}
    makeblastdb = _find_makeblastdb(blastn)
    if not makeblastdb:
        return {"status": "unavailable",
                "reason": "makeblastdb not found alongside blastn",
                "tool": "makeblastdb", "per_gene": {}, "per_cluster": {},
                "per_hit": [], "db_sha": db_sha}

    # `-max_target_seqs` bounds HSPs PER QUERY allele; raise it for the all-hits mode so a
    # tandem array of one allele is not truncated before _interval_dedup can split copies.
    max_target = "10000" if all_hits else "5"
    with tempfile.TemporaryDirectory(prefix="vf_canonical_") as td:
        asm_db = str(Path(td) / "asm")
        try:
            # Do NOT pass `-parse_seqids`: verified empirically 2026-06-19 that WITHOUT it
            # blastn `sseqid` is the EXACT FASTA header first-token (e.g. `CP021689.1`) —
            # the same token AMRFinder reports — so the shared genome_map
            # `build_contig_name_map` reconciles both overlays. WITH `-parse_seqids` the id
            # is mangled (`gb|CP021689.1|`) and the contig-name map silently breaks.
            subprocess.run([makeblastdb, "-in", str(fasta), "-dbtype", "nucl",
                            "-out", asm_db], check=True, capture_output=True,
                           text=True, timeout=timeout)
            # query = VF alleles, subject DB = the assembly → coverage is of the reference allele.
            proc = subprocess.run(
                [blastn, "-query", str(db), "-db", asm_db,
                 "-outfmt", "6 qseqid sseqid sstart send pident length qlen",
                 "-perc_identity", str(identity_threshold),
                 "-max_target_seqs", max_target, "-evalue", "1e-20"],
                check=True, capture_output=True, text=True, timeout=timeout)
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired, OSError) as e:
            return {"status": "unavailable",
                    "reason": f"blastn invocation failed: {type(e).__name__}",
                    "tool": "blastn", "per_gene": {}, "per_cluster": {},
                    "per_hit": [], "db_sha": db_sha}

        calls = parse_blastn_outfmt6(proc.stdout, identity_threshold, coverage_threshold)

    return {
        "status": "ok",
        "tool": "blastn",
        "method": "canonical_virulencefinder_blastn_v0_1",
        "all_hits": all_hits,
        "parameters": {"identity_threshold": identity_threshold,
                       "coverage_threshold": coverage_threshold},
        "per_gene": calls["per_gene"],
        "per_cluster": calls["per_cluster"],
        "per_hit": calls["per_hit"],
        "db_sha": db_sha,
    }


# the immutable honesty caveat carried into every diff section (ratified Q2, 2026-06-04).
NON_INDEPENDENCE_CAVEAT = (
    "Both callers match the SAME VirulenceFinder DB; high concordance is EXPECTED. "
    "This is an audit of the fast k-mer-seed caller against canonical blastn over one "
    "DB, NOT independent validation."
)


def build_vf_diff(resolver_profile: dict, canonical: dict, *,
                  resolver_marker_hits: list | None = None) -> dict:
    """Side-by-side diff: resolver k-mer cluster calls vs canonical blastn cluster calls.

    Per-cluster concordance (both_present / both_absent / disagree) at each caller's
    confident bar, a concordance summary, and the non-independence honesty flag + caveat.
    Degrades to status="unavailable" (section retained) when canonical is unavailable."""
    cov_by_cluster = {}
    for h in (resolver_marker_hits or []):
        if h.get("hit_status") == "CONFIDENT":
            cov_by_cluster[h["cluster"]] = h.get("percent_coverage")

    if canonical.get("status") != "ok":
        return {
            "status": "unavailable",
            "reason": canonical.get("reason", "canonical VirulenceFinder unavailable"),
            "caller_is_independent_baseline": False,
            "caveat": NON_INDEPENDENCE_CAVEAT,
            "resolver_confident_coverage": RESOLVER_CONFIDENT_COV,
            "per_gene": [],
            "per_cluster": [],
            "concordance": None,
        }

    canon_cluster = canonical.get("per_cluster", {})
    per_cluster = []
    both_present = both_absent = disagree = 0
    for cluster in sorted(CLUSTER_MARKERS):
        r_call = bool(resolver_profile.get(cluster, False))
        c_call = bool(canon_cluster.get(cluster, {}).get("called", False))
        if r_call and c_call:
            agreement, both_present = "both_present", both_present + 1
        elif (not r_call) and (not c_call):
            agreement, both_absent = "both_absent", both_absent + 1
        else:
            agreement, disagree = "disagree", disagree + 1
        per_cluster.append({
            "cluster": cluster,
            "resolver_called": r_call,
            "resolver_percent_coverage": cov_by_cluster.get(cluster),
            "canonical_called": c_call,
            "canonical_percent_identity": canon_cluster.get(cluster, {}).get("percent_identity"),
            "canonical_percent_coverage": canon_cluster.get(cluster, {}).get("percent_coverage"),
            "canonical_best_gene": canon_cluster.get(cluster, {}).get("best_gene"),
            "agreement": agreement,
        })

    per_gene = [
        {"gene": gid, "cluster": info["cluster"],
         "canonical_called": info["called"],
         "canonical_percent_identity": info["percent_identity"],
         "canonical_percent_coverage": info["percent_coverage"]}
        for gid, info in sorted(canonical.get("per_gene", {}).items())
    ]

    total = both_present + both_absent + disagree
    concordant = both_present + both_absent
    return {
        "status": "ok",
        "tool": canonical.get("tool", "blastn"),
        "method": canonical.get("method"),
        "parameters": canonical.get("parameters"),
        "caller_is_independent_baseline": False,
        "caveat": NON_INDEPENDENCE_CAVEAT,
        "resolver_confident_coverage": RESOLVER_CONFIDENT_COV,
        "per_gene": per_gene,
        "per_cluster": per_cluster,
        "concordance": {
            "n_clusters": total,
            "both_present": both_present,
            "both_absent": both_absent,
            "disagree": disagree,
            "per_cluster_concordance": round(concordant / total, 3) if total else None,
        },
    }
