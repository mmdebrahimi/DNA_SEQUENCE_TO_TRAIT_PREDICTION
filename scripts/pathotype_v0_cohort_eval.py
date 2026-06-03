"""v0 resolver full-cohort evaluation (ledger Goal 5 / H4 abstention quality).

Runs detection + resolution over the 24 cached ExPEC/EPEC genomes and tabulates
calls vs labels: per-label call distribution, supported-class recall, confident-call
precision, abstention rate + enrichment. This is the data the ExPEC-sensitivity
calibration decision (H4) needs.

SPEED: the only expensive step is the per-genome k=15 k-mer set + per-cluster
coverage. We cache `cluster -> [best_gene, best_cov]` to data/pathotype_cov_cache/
per genome, so re-running after a resolver-rule change is instant (no re-detection).
Detection clusters are fixed by markers.py; if CLUSTER_MARKERS changes, bump CACHE_VER.

External-validity discipline (ledger v5): ExPEC labels are isolation-site-derived
(independent) -> ExPEC recall is a genotype->phenotype signal. EPEC labels are
DECA-curated. Report per-class, never blended into one accuracy number.
"""
from __future__ import annotations
import json, sys
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from scripts.pathotype_kmer_bakeoff import load_strains
from dna_decode.pathotype.detect import (
    _genome_kmers, _coverage, build_vf_index, ANCHORS, CONFIDENT_COV, PARTIAL_COV, K,
)
from dna_decode.pathotype.markers import CLUSTER_MARKERS
from dna_decode.pathotype.resolve import resolve_call

CACHE_VER = "v1"
CACHE = REPO / "data/pathotype_cov_cache"
DB = REPO / "data/virulencefinder_db/virulence_ecoli.fsa"
# expected supported call per label (0=ExPEC/Salipante isolation-site, 1=EPEC/Hazen)
EXPEC_OK = {"UPEC_COMPATIBLE", "ExPEC_COMPATIBLE"}
EPEC_OK = {"tEPEC_COMPATIBLE", "aEPEC_COMPATIBLE"}
ABSTAIN = {"AMBIGUOUS", "AMBIGUOUS_LOW_QC", "UNCLASSIFIED"}


def cluster_coverage(genome_id: str, seq: str, vf_index: dict) -> dict[str, list]:
    CACHE.mkdir(parents=True, exist_ok=True)
    cf = CACHE / f"{genome_id}_{CACHE_VER}.json"
    if cf.exists():
        return json.loads(cf.read_text())
    gset = _genome_kmers([("g", seq)], K)
    out = {}
    for cluster, alleles in vf_index.items():
        anchor = ANCHORS.get(cluster)
        pool = [(g, s) for g, s in alleles if (anchor is None or g.lower().startswith(anchor))]
        best_gene, best_cov = None, 0.0
        for gene, s in pool:
            c = _coverage(s, gset, K)
            if c > best_cov:
                best_gene, best_cov = gene, c
        out[cluster] = [best_gene, round(best_cov, 4)]
    cf.write_text(json.dumps(out))
    return out


def profile_from_coverage(cov: dict[str, list]):
    profile, partial = {}, set()
    for cluster, (_, c) in cov.items():
        if c >= CONFIDENT_COV:
            profile[cluster] = True
        elif c >= PARTIAL_COV:
            partial.add(cluster)
    return profile, frozenset(partial)


def main() -> int:
    seqs, labels, ids = load_strains()
    vf_index = build_vf_index(DB)
    print(f"[eval] {len(ids)} genomes; detecting (cached per-genome coverage)...", flush=True)

    rows = []
    for sid in ids:
        gid = sid.split("|")[0]
        cov = cluster_coverage(gid, seqs[sid], vf_index)
        profile, partial = profile_from_coverage(cov)
        call = resolve_call(profile, partial_clusters=partial, qc_pass=True)
        rows.append({"id": gid, "y": labels[sid], "call": call["primary"],
                     "tier": call["confidence_tier"], "validity": call["external_validity"],
                     "clusters": sorted(profile)})
        print(f"[eval] {gid} (y={labels[sid]}): {call['primary']} [{call['confidence_tier']}]", flush=True)

    expec = [r for r in rows if r["y"] == 0]
    epec = [r for r in rows if r["y"] == 1]

    def rate(rs, ok):
        return sum(1 for r in rs if r["call"] in ok) / len(rs) if rs else 0.0

    metrics = {
        "n_expec": len(expec), "n_epec": len(epec),
        "expec_recall": rate(expec, EXPEC_OK),
        "epec_recall": rate(epec, EPEC_OK),
        "expec_abstain_rate": rate(expec, ABSTAIN),
        "epec_abstain_rate": rate(epec, ABSTAIN),
        "expec_call_dist": dict(Counter(r["call"] for r in expec)),
        "epec_call_dist": dict(Counter(r["call"] for r in epec)),
    }
    # confident-call precision: among CONFIDENT supported calls, fraction matching label
    conf = [r for r in rows if r["tier"] == "CONFIDENT" and r["call"] in (EXPEC_OK | EPEC_OK)]
    correct = sum(1 for r in conf
                  if (r["y"] == 0 and r["call"] in EXPEC_OK) or (r["y"] == 1 and r["call"] in EPEC_OK))
    metrics["confident_supported_precision"] = (correct / len(conf)) if conf else None
    metrics["n_confident_supported"] = len(conf)

    print("\n[eval] === H4 metrics ===")
    print(f"[eval] EPEC recall (->tEPEC/aEPEC):  {metrics['epec_recall']:.2f}  dist={metrics['epec_call_dist']}")
    print(f"[eval] ExPEC recall (->UPEC/ExPEC):  {metrics['expec_recall']:.2f}  dist={metrics['expec_call_dist']}")
    print(f"[eval] abstain rate  ExPEC={metrics['expec_abstain_rate']:.2f}  EPEC={metrics['epec_abstain_rate']:.2f}")
    print(f"[eval] confident supported-call precision: {metrics['confident_supported_precision']} "
          f"(n={metrics['n_confident_supported']})")

    res = {"contrast": "ExPEC(Salipante) vs EPEC(Hazen)", "resolver_version": "see markers.RULES_VERSION",
           "metrics": metrics, "per_genome": rows}
    out = REPO / "research_outputs/pathotype_v0_cohort_eval_2026-06-03.json"
    out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(f"[eval] wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
