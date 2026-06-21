"""Virulence-determinant overlay (Step 3) — the genome-map's 5th tier.

Surfaces curated VirulenceFinder (VF) allele PRESENCE behind the SAME coordinate-join
integrity gate as the AMR `determinant-phenotype` tier (presence of a curated
determinant, NEVER a learned pathogenicity claim), plus a SEPARATE genome-level
pathotype call computed via the FULL deployed `resolve_call` honesty contract.

Four pieces:
  - `parse_virulence_hits` : a canonical-VF `per_hit` list -> DeterminantHits (so the
    AMR coord-join machinery — `join_hits` / `build_contig_name_map` — is reused verbatim).
  - `join_virulence`       : thin reuse of `phenotype_overlay.join_hits`; emits VF-namespaced
    join-quality counts (M2 metric isolation) + an ambiguous-contig flag (M1).
  - `genome_pathotype_call`: the genome-level overlay (the virulence analog of the AMR
    R/S call), through `resolve_call` with qc_pass + cross-axis ExPEC support (C1 — a
    low-QC genome can never yield a confident commensal call).
  - `virulence_organism_in_scope` : E. coli / Shigella gate (the committed v1 DB is E. coli).

HONESTY: the VF caller matches the SAME curated DB the pathotype resolver uses — this is
located-determinant PRESENCE, not independent validation; every surfaced feature carries
`NON_INDEPENDENCE_CAVEAT` + the DB SHA. This module is READ-ONLY w.r.t. the frozen AMR
surface and never touches `build_vf_diff`'s cluster-audit semantics.
"""
from __future__ import annotations

from collections import Counter

from dna_decode.genome_map.phenotype_overlay import DeterminantHit, join_hits
from dna_decode.pathotype.detect import assembly_qc
from dna_decode.pathotype.expec_score import (
    EXPEC_SUPPORT_GENE_PREFIXES,
    meets_cross_axis_support,
    support_gene_count,
)
from dna_decode.pathotype.resolve import resolve_call
from dna_decode.pathotype.vf_runner import NON_INDEPENDENCE_CAVEAT

VIRULENCE_CLASS = "VIRULENCE"

# Which pathotype context each marker cluster contributes to (clustered hits only — used
# for the per-feature virulence wall; an unclustered VF allele has no pathotype context).
CLUSTER_PATHOTYPE_CONTEXT: dict[str, list[str]] = {
    "STX1": ["STEC/EHEC"], "STX2": ["STEC/EHEC"],
    "LEE": ["EPEC/EHEC"], "BFP_EAF": ["tEPEC"],
    "LT": ["ETEC"], "ST": ["ETEC"],
    "EAEC_REG": ["EAEC"], "EAEC_TRANSPORT": ["EAEC"], "EAEC_T6SS": ["EAEC"],
    "AAF_I": ["EAEC"], "AAF_II": ["EAEC"], "AAF_III": ["EAEC"],
    "AAF_IV": ["EAEC"], "AAF_V": ["EAEC"],
    "P_FIMBRIAE": ["ExPEC/UPEC"], "S_FIMBRIAE": ["ExPEC/UPEC"],
    "AFA_DRA": ["ExPEC/DAEC"], "HEMOLYSIN": ["ExPEC"], "CNF1": ["ExPEC"],
    "SIDEROPHORES": ["ExPEC-support"], "CAPSULE_SERUM": ["ExPEC-support"],
    "EIEC_FLAG": ["EIEC"], "DAEC_FLAG": ["DAEC"],
}


def cluster_pathotype_context(cluster: str | None) -> list[str]:
    """The pathotype(s) a marker cluster contributes to ([] for unclustered alleles)."""
    if not cluster:
        return []
    return list(CLUSTER_PATHOTYPE_CONTEXT.get(cluster, []))


def parse_virulence_hits(canonical_vf_result: dict) -> list[DeterminantHit]:
    """Map a canonical-VF result's `per_hit` list to DeterminantHits (called hits only).

    `protein_id` is always None (VF has no protein-id column) so the join falls to the
    coord path; `contig` = the FASTA first-token `sseqid` (the shared contig-name-map
    contract). `subclass` carries the marker cluster ("" when unclustered)."""
    out: list[DeterminantHit] = []
    for h in (canonical_vf_result or {}).get("per_hit", []):
        if not h.get("called"):
            continue
        out.append(DeterminantHit(
            symbol=str(h.get("vf_gene") or ""),
            name=str(h.get("allele_id") or ""),
            cls=VIRULENCE_CLASS,
            subclass=str(h.get("cluster") or ""),
            method="blastn",
            protein_id=None,
            contig=h.get("sseqid"),
            start=h.get("start"),
            stop=h.get("stop"),
            raw=dict(h),
        ))
    return out


def _count_ambiguous_contig_hits(hits: list[DeterminantHit], contig_names) -> int:
    """M1: # hits whose contig first-token is NON-UNIQUE among the genome's contigs.

    A duplicate first-token makes the coordinate join ambiguous (the contig-name map
    can't disambiguate two contigs with the same header token) — surfaced as degraded
    join quality, never silently coord-joined to the wrong contig."""
    if not contig_names:
        return 0
    dup = {name for name, c in Counter(contig_names).items() if c > 1}
    return sum(1 for h in hits if h.contig in dup)


def all_virulence_symbol_fallback(counts: dict) -> bool:
    """True iff there are VF rows AND none cleared a high-confidence join (the trap guard,
    VF-namespaced so it never feeds the AMR `all_joins_symbol_fallback`)."""
    return counts.get("n_vf_rows", 0) > 0 and counts.get("n_high_confidence_join", 0) == 0


def join_virulence(features, hits: list[DeterminantHit], *,
                   contig_name_map: dict[str, str] | None = None,
                   contig_names=None) -> tuple[list, dict]:
    """Join VF DeterminantHits to features (reusing the AMR coord-join machinery).

    Returns (joined_hits, counts) with VF-namespaced keys so no virulence number can
    feed the AMR join metrics or the AMR spike gate (M2):
      {n_vf_rows, n_high_confidence_join, n_symbol_fallback, n_unjoined,
       n_ambiguous_contig, all_virulence_joins_symbol_fallback}.
    `contig_names` is the multiset of FASTA contig first-tokens (for the M1 ambiguity flag).
    """
    joined, raw = join_hits(features, hits, contig_name_map=contig_name_map)
    counts = {
        "n_vf_rows": raw["n_main_rows"],
        "n_high_confidence_join": raw["n_high_confidence_join"],
        "n_symbol_fallback": raw["n_symbol_fallback"],
        "n_unjoined": raw["n_unjoined"],
        "n_ambiguous_contig": _count_ambiguous_contig_hits(hits, contig_names),
    }
    counts["all_virulence_joins_symbol_fallback"] = all_virulence_symbol_fallback(counts)
    return joined, counts


def _pergene_cov_from_per_gene(per_gene: dict) -> dict[str, float]:
    """Build the ExPEC per-gene coverage dict `{support-prefix: max coverage in [0,1]}`.

    Keys are the EXPEC_SUPPORT_GENE_PREFIXES (lowercase) the resolver's support scoring
    expects; coverage is the canonical-VF percent_coverage/100 (max across alleles that
    share a support prefix)."""
    out: dict[str, float] = {}
    for allele_id, info in (per_gene or {}).items():
        gl = str(allele_id).split(":")[0].lower()
        cov = float(info.get("percent_coverage") or 0.0) / 100.0
        for prefix in EXPEC_SUPPORT_GENE_PREFIXES:
            if gl.startswith(prefix):
                out[prefix] = max(out.get(prefix, 0.0), cov)
    return out


def genome_pathotype_call(canonical_vf_result: dict | None, contigs) -> dict:
    """Genome-level pathotype overlay via the FULL deployed resolve_call contract (C1).

    Computes the cluster profile from the canonical per-cluster calls (clustered only),
    the per-gene ExPEC support from per_gene, the assembly QC from `contigs`, and feeds
    ALL of qc_pass + support_gene_count + cross_axis_support into `resolve_call` — so a
    low-QC genome yields AMBIGUOUS_LOW_QC (never a confident commensal call), exactly as
    the deployed CLI does.

    `status="insufficient_context"` ONLY when the VF result is unavailable/errored OR
    `contigs` is missing — NEVER merely because support is absent (that is a real,
    confident COMMENSAL_LOW_MARKER_BURDEN call).
    """
    caveat = NON_INDEPENDENCE_CAVEAT
    db_sha = (canonical_vf_result or {}).get("db_sha")
    if not canonical_vf_result or canonical_vf_result.get("status") != "ok" or not contigs:
        return {
            "status": "insufficient_context",
            "reason": "VF caller unavailable/errored or no contigs to QC",
            "derived_call": None,
            "assembly_qc": None,
            "caveat": caveat,
            "db_sha": db_sha,
        }

    per_cluster = canonical_vf_result.get("per_cluster", {})
    profile = {c: bool(info.get("called")) for c, info in per_cluster.items()}
    pergene_cov = _pergene_cov_from_per_gene(canonical_vf_result.get("per_gene", {}))

    qc = assembly_qc(contigs)
    qc_pass = qc["qc_verdict"] != "FAIL"
    derived = resolve_call(
        profile,
        qc_pass=qc_pass,
        support_gene_count=support_gene_count(pergene_cov),
        cross_axis_support=meets_cross_axis_support(pergene_cov),
    )
    return {
        "status": "ok",
        "derived_call": derived,
        "assembly_qc": qc,
        "caveat": caveat,
        "db_sha": db_sha,
    }


def virulence_organism_in_scope(organism: str | None) -> bool:
    """True iff the AMRFinder `-O` organism is E. coli / Shigella (the committed v1 VF DB).

    Genus-token based, case/underscore-insensitive — so `Escherichia`,
    `Escherichia_coli`, `Escherichia coli`, `Escherichia_coli_Shigella`, and any
    `Shigella*` resolve True; Klebsiella / Salmonella / etc. resolve False."""
    if not organism:
        return False
    tokens = str(organism).replace("_", " ").strip().lower().split()
    if not tokens:
        return False
    return tokens[0] in {"escherichia", "shigella"}
