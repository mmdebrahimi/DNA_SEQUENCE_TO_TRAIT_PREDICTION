"""v0 marker detection engine: genome FASTA -> per-cluster presence + provenance.

Pure-Python k-mer seed presence (k=15, both strands) against the VirulenceFinder
E. coli allele DB. Deliberately dependency-light (NO BLAST/KMA) so the CLI runs on
a laptop with only numpy-free stdlib. Detection is the ONLY heavy step; the resolver
(resolve.py) is pure logic on the profile this produces.

Honesty: coverage here is k-mer-seed coverage, a presence proxy — NOT BLAST percent
identity. We report `percent_coverage` and leave `percent_identity` null with
`method="kmer_seed_coverage_v0"`. A true VirulenceFinder/BLAST side-by-side diff
(ledger requirement) is a v0.1 add that needs the VF software installed.
"""
from __future__ import annotations

from pathlib import Path

from dna_decode.pathotype.markers import CLUSTER_MARKERS

K = 15
CONFIDENT_COV = 0.80
PARTIAL_COV = 0.65
_COMP = str.maketrans("ACGTacgt", "TGCAtgca")

# clusters whose presence requires a specific anchor gene (not just any member gene).
ANCHORS = {"LEE": "eae", "BFP_EAF": "bfp", "EAEC_REG": "aggr",
           "EAEC_TRANSPORT": "aata", "EAEC_T6SS": "aaic"}


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def parse_fasta(path: str | Path) -> list[tuple[str, str]]:
    out, name, buf = [], None, []
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(">"):
            if name is not None:
                out.append((name, "".join(buf)))
            name = line[1:].split()[0] if line[1:].split() else line[1:]
            buf = []
        else:
            buf.append(line.strip())
    if name is not None:
        out.append((name, "".join(buf)))
    return out


def build_vf_index(db_path: str | Path) -> dict[str, list[tuple[str, str]]]:
    """cluster -> list[(gene_name, allele_seq)] for the alleles that match the cluster."""
    alleles: list[tuple[str, str]] = []
    name, buf = None, []
    for line in Path(db_path).read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(">"):
            if name is not None:
                alleles.append((name, "".join(buf)))
            name = line[1:].split(":")[0]
            buf = []
        else:
            buf.append(line.strip())
    if name is not None:
        alleles.append((name, "".join(buf)))

    index: dict[str, list[tuple[str, str]]] = {c: [] for c in CLUSTER_MARKERS}
    for gene, seq in alleles:
        gl = gene.lower()
        for cluster, prefixes in CLUSTER_MARKERS.items():
            if any(gl.startswith(p) for p in prefixes):
                index[cluster].append((gene, seq))
    return index


def _genome_kmers(contigs: list[tuple[str, str]], k: int) -> set:
    out = set()
    for _, seq in contigs:
        s = seq.upper()
        for i in range(len(s) - k + 1):
            km = s[i:i + k]
            if "N" not in km:
                out.add(km)
    return out


def _coverage(seq: str, gset: set, k: int) -> float:
    a = seq.upper()
    if len(a) < k:
        return 0.0
    def cov(x):
        kms = [x[i:i + k] for i in range(len(x) - k + 1)]
        kms = [z for z in kms if "N" not in z]
        return sum(1 for z in kms if z in gset) / len(kms) if kms else 0.0
    return max(cov(a), cov(revcomp(a)))


def _locate(allele: str, contigs: list[tuple[str, str]]) -> tuple[str | None, int | None, str]:
    """Find a representative seed of allele in the contigs -> (contig, start, strand)."""
    a = allele.upper()
    seed_len = min(31, len(a))
    mid = max(0, len(a) // 2 - seed_len // 2)
    seed = a[mid:mid + seed_len]
    for cid, seq in contigs:
        s = seq.upper()
        pos = s.find(seed)
        if pos >= 0:
            return cid, pos, "+"
        rpos = s.find(revcomp(seed))
        if rpos >= 0:
            return cid, rpos, "-"
    return None, None, "+"


def assembly_qc(contigs: list[tuple[str, str]]) -> dict:
    lengths = sorted((len(s) for _, s in contigs), reverse=True)
    total = sum(lengths)
    half, run, n50 = total / 2, 0, 0
    for L in lengths:
        run += L
        if run >= half:
            n50 = L
            break
    if 3_500_000 <= total <= 7_000_000 and len(contigs) <= 2000:
        verdict = "PASS"
    elif 2_500_000 <= total <= 8_000_000:
        verdict = "WARN"
    else:
        verdict = "FAIL"
    return {"n_contigs": len(contigs), "total_bp": total, "n50": n50, "qc_verdict": verdict}


def detect(contigs: list[tuple[str, str]], vf_index: dict) -> dict:
    """Return {cluster_profile, partial_clusters, marker_hits}."""
    gset = _genome_kmers(contigs, K)
    profile: dict[str, bool] = {}
    partial: set[str] = set()
    hits: list[dict] = []

    for cluster, alleles in vf_index.items():
        if not alleles:
            profile[cluster] = False
            continue
        anchor = ANCHORS.get(cluster)
        # best-covered allele, optionally restricted to the anchor gene
        pool = [(g, s) for g, s in alleles if (anchor is None or g.lower().startswith(anchor))]
        if not pool:
            profile[cluster] = False
            continue
        best_gene, best_cov, best_seq = None, 0.0, None
        for gene, seq in pool:
            c = _coverage(seq, gset, K)
            if c > best_cov:
                best_gene, best_cov, best_seq = gene, c, seq
        present = best_cov >= CONFIDENT_COV
        profile[cluster] = present
        if not present and best_cov >= PARTIAL_COV:
            partial.add(cluster)
        if best_cov >= PARTIAL_COV:
            cid, start, strand = _locate(best_seq, contigs)
            hits.append({
                "cluster": cluster, "gene": best_gene, "allele": None,
                "source_db": "virulencefinder_ecoli",
                "percent_identity": None, "percent_coverage": round(best_cov * 100, 1),
                "method": "kmer_seed_coverage_v0",
                "contig": cid, "start": start,
                "end": (start + len(best_seq)) if start is not None else None,
                "strand": strand,
                "hit_status": "CONFIDENT" if present else "PARTIAL",
            })
    return {"cluster_profile": profile, "partial_clusters": frozenset(partial), "marker_hits": hits}
