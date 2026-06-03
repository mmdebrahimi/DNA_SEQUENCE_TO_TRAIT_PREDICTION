"""Build a COMMITTED per-gene support-marker coverage cache for the 24-genome H4 cohort.

Mission ep4-v01-expec-recall: the existing cov cache (data/pathotype_cov_cache/) stores only the
single best-covered allele per CLUSTER, so multi-gene extraintestinal burden in SIDEROPHORES /
CAPSULE_SERUM is invisible. This computes per-GENE coverage (max over a gene's alleles) for the two
ExPEC support clusters across all 24 cached genomes and persists it to data/pathotype_pergene_cache/
so the resolver test runs offline (no ENA fetch, no re-detect).

Per-gene group = the member gene prefixes of SIDEROPHORES + CAPSULE_SERUM in markers.CLUSTER_MARKERS.
A gene is PRESENT iff its best allele coverage >= detect.CONFIDENT_COV (0.80) — same bar as cluster
presence, so the per-gene count is a strict refinement of the cluster boolean.
"""
from __future__ import annotations
import csv, json, re, sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from dna_decode.pathotype.detect import _genome_kmers, _coverage, build_vf_index, K
from dna_decode.pathotype.markers import CLUSTER_MARKERS

CACHE_VER = "v1"
CACHE = REPO / "data/pathotype_pergene_cache"
DB = REPO / "data/virulencefinder_db/virulence_ecoli.fsa"
SUPPORT_CLUSTERS = ["SIDEROPHORES", "CAPSULE_SERUM"]

# --- offline strain loader (replicates scripts.pathotype_kmer_bakeoff without its sklearn import) ---
F1 = REPO / "data/external/horesh2021_F1_genome_metadata.csv"
ENA_CACHE = REPO / "data/ena_wgs"
N_PER_CLASS = 12
WGS_MASTER = re.compile(r"^[A-Z]{4}\d{8}(\.fa)?$", re.I)


def _wgs_set_prefix(assembly_name: str) -> str:
    base = assembly_name.replace(".fa", "").replace(".fasta", "").strip()
    return base[:4] + "01"


def _concat_contigs(fasta_text: str) -> str:
    return "".join(l.strip() for l in fasta_text.splitlines() if not l.startswith(">"))


def load_strains():
    rows = list(csv.DictReader(open(F1, encoding="utf-8")))
    clean = [r for r in rows if "(predicted)" not in r["Pathotype"]
             and r["Pathotype"] not in ("Not determined", "")]
    ok = lambda r: bool(WGS_MASTER.match(r["Assembly_name"].strip()))
    expec = [r for r in clean if r["Pathotype"].strip().startswith("ExPEC")
             and r["Source"].startswith("Salipante") and ok(r)][:N_PER_CLASS]
    epec = [r for r in clean if r["Pathotype"].strip().startswith("EPEC")
            and r["Source"].startswith("Hazen") and ok(r)][:N_PER_CLASS]
    roster = [(r["ID"], r["Assembly_name"], 0) for r in expec] + \
             [(r["ID"], r["Assembly_name"], 1) for r in epec]
    seqs, labels, ids = {}, {}, []
    for sid, aname, y in roster:
        setp = _wgs_set_prefix(aname)
        f = ENA_CACHE / f"{setp}.fna"
        if not (f.exists() and f.stat().st_size > 1000):
            print(f"  skip {sid}: no cached assembly {setp}", flush=True)
            continue
        seq = _concat_contigs(f.read_text(encoding="utf-8"))
        if len(seq) < 1_000_000:
            print(f"  skip {sid}: too short ({len(seq)} bp)", flush=True)
            continue
        key = f"{sid}|{setp}"
        seqs[key], labels[key] = seq, y
        ids.append(key)
    return seqs, labels, ids


def support_gene_prefixes() -> dict[str, list[str]]:
    """cluster -> list of per-gene prefixes (each prefix == one 'gene' for counting)."""
    return {c: list(CLUSTER_MARKERS[c]) for c in SUPPORT_CLUSTERS}


def pergene_coverage(genome_id: str, seq: str, vf_index: dict, gset=None) -> dict[str, float]:
    """gene_prefix -> best allele coverage among that gene's alleles in the support clusters."""
    if gset is None:
        gset = _genome_kmers([("g", seq)], K)
    out: dict[str, float] = {}
    for cluster in SUPPORT_CLUSTERS:
        for prefix in CLUSTER_MARKERS[cluster]:
            best = 0.0
            for gene, s in vf_index[cluster]:
                if gene.lower().startswith(prefix):
                    c = _coverage(s, gset, K)
                    if c > best:
                        best = c
            out[prefix] = round(best, 4)
    return out


def main() -> int:
    seqs, labels, ids = load_strains()
    vf_index = build_vf_index(DB)
    CACHE.mkdir(parents=True, exist_ok=True)
    print(f"[pergene] {len(ids)} genomes; computing per-gene support coverage...", flush=True)
    for sid in ids:
        gid = sid.split("|")[0]
        cf = CACHE / f"{gid}_{CACHE_VER}.json"
        if cf.exists():
            print(f"[pergene] {gid}: cached", flush=True)
            continue
        cov = pergene_coverage(gid, seqs[sid], vf_index)
        cf.write_text(json.dumps(cov), encoding="utf-8")
        present = [g for g, c in cov.items() if c >= 0.80]
        print(f"[pergene] {gid} (y={labels[sid]}): present={present}", flush=True)
    print(f"[pergene] wrote {CACHE}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
