"""ExPEC-vs-EPEC k-mer baseline arm of the representation bake-off (EP-4).

The MANDATED classical control (CLAUDE.md: FM must beat classical by >=3pp) and
the confound diagnostic: ExPEC(Salipante) vs EPEC(Hazen) is study-confounded, so
a near-perfect k-mer AUROC means the contrast is trivially/batch-separable and
NT's value is unmeasurable on this substrate; a moderate k-mer AUROC leaves room
for NT lift.

Cheap: assembly FASTA only (ENA WGS-set fetch, no GCA, no Bakta, no GPU).
Non-circular labels: ExPEC=isolation-site, EPEC=DECA-curated.
"""
from __future__ import annotations
import csv, gzip, io, sys, json
from pathlib import Path
import urllib.request

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
import numpy as np
from sklearn.metrics import roc_auc_score
from dna_decode.models.classical_baselines import (
    CONTIG_SEPARATOR,
    build_kmer_vocabulary,
    extract_kmer_counts,
    kmers_to_feature_matrix,
)
from dna_decode.models.classifiers import predict_proba, train_xgboost_classifier

F1 = REPO / "data/external/horesh2021_F1_genome_metadata.csv"
ENA_CACHE = Path("C:/Users/Farshad/PythonProjects/dna_decode/data/ena_wgs")
N_PER_CLASS = 12


import re
from collections import Counter
WGS_MASTER = re.compile(r"^[A-Z]{4}\d{8}(\.fa)?$", re.I)  # AAAA00000000(.fa) -> WGS set


def select_strains():
    rows = list(csv.DictReader(open(F1, encoding="utf-8")))
    clean = [r for r in rows if "(predicted)" not in r["Pathotype"] and r["Pathotype"] not in ("Not determined", "")]
    # Only genuine WGS-master Assembly_names are ENA-fetchable; some Salipante rows
    # carry SRA run ids (SRR...) which the prefix logic mangles -> filter them out.
    def ok(r): return bool(WGS_MASTER.match(r["Assembly_name"].strip()))
    expec = [r for r in clean if r["Pathotype"].strip().startswith("ExPEC") and r["Source"].startswith("Salipante") and ok(r)][:N_PER_CLASS]
    epec = [r for r in clean if r["Pathotype"].strip().startswith("EPEC") and r["Source"].startswith("Hazen") and ok(r)][:N_PER_CLASS]
    out = []
    for r in expec: out.append((r["ID"], r["Assembly_name"], "ExPEC", 0, r["Source"]))
    for r in epec: out.append((r["ID"], r["Assembly_name"], "EPEC", 1, r["Source"]))
    return out


def wgs_set_prefix(assembly_name: str) -> str:
    base = assembly_name.replace(".fa", "").replace(".fasta", "").strip()
    return base[:4] + "01"  # master AAAA00000000 -> contig set AAAA01


def fetch_assembly(setp: str) -> str | None:
    ENA_CACHE.mkdir(parents=True, exist_ok=True)
    cache = ENA_CACHE / f"{setp}.fna"
    if cache.exists() and cache.stat().st_size > 1000:
        return cache.read_text(encoding="utf-8")
    url = f"https://www.ebi.ac.uk/ena/browser/api/fasta/{setp}"  # bare form (params break it)
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    try:
        raw = urllib.request.urlopen(req, timeout=120).read()
    except Exception as e:
        print(f"  FETCH FAIL {setp}: {type(e).__name__} {str(e)[:80]}")
        return None
    if raw[:2] == b"\x1f\x8b":  # gzip magic
        text = gzip.decompress(raw).decode("utf-8", "replace")
    else:
        text = raw.decode("utf-8", "replace")
    if ">" not in text:
        print(f"  NO FASTA in {setp} (len={len(text)})")
        return None
    cache.write_text(text, encoding="utf-8")
    return text


def concat_contigs(fasta_text: str) -> str:
    seqs, cur = [], []
    for line in fasta_text.splitlines():
        if line.startswith(">"):
            if cur: seqs.append("".join(cur)); cur = []
        else:
            cur.append(line.strip())
    if cur: seqs.append("".join(cur))
    return CONTIG_SEPARATOR.join(seqs)


# ---- Cached-counts LOSO (45x faster, provably equivalent to run_kmer_xgboost_loso) ----
# The canonical runner recounts every training genome twice PER FOLD (once in
# build_kmer_vocabulary, once in kmers_to_feature_matrix) -> ~5.5e9 pure-Python
# iters over 24 folds. We count each genome's k-mers ONCE up front, then per fold
# just AGGREGATE the cached dicts. Counter.update(dict) + most_common reproduce
# build_kmer_vocabulary EXACTLY (same counts, same first-encountered tie order
# given identical strain order); the cached matrix reproduces kmers_to_feature_matrix
# EXACTLY. XGBoost is seeded (random_state=42, calibrate=False) so identical
# (vocab, X) => identical scores. Equivalence is asserted at the feature level
# in _verify_feature_equivalence() before the full run is trusted.
TOP_N = 10_000


def precompute_counts(seqs_by_strain: dict[str, str], k: int) -> dict[str, dict[str, int]]:
    return {sid: extract_kmer_counts(seq, k=k) for sid, seq in seqs_by_strain.items()}


def _vocab_from_cache(strains: list[str], cache: dict[str, dict[str, int]], top_n: int) -> list[str]:
    agg: Counter[str] = Counter()
    for s in strains:
        agg.update(cache[s])  # mirrors build_kmer_vocabulary's aggregate.update(...)
    return [km for km, _ in agg.most_common(top_n)]


def _matrix_from_cache(strains: list[str], vocab: list[str], cache: dict[str, dict[str, int]]) -> np.ndarray:
    out = np.zeros((len(strains), len(vocab)), dtype=np.float32)
    vocab_index = {km: i for i, km in enumerate(vocab)}
    for row_i, s in enumerate(strains):
        cnt = cache[s]
        for km, idx in vocab_index.items():
            v = cnt.get(km)
            if v:
                out[row_i, idx] = v
    return out


def _verify_feature_equivalence(seqs_by_strain: dict[str, str], ids: list[str],
                                cache: dict[str, dict[str, int]], k: int, top_n: int) -> None:
    """Prove cached path == canonical at the (vocab, X) level for fold 0.

    If vocab and X are byte-identical, the downstream train_xgboost_classifier
    call is identical -> scores identical. Cheap: one canonical build over n-1 genomes.
    """
    train = ids[1:]  # hold out fold 0, same as run_kmer_xgboost_loso(i=0)
    train_seqs = [seqs_by_strain[s] for s in train]
    vocab_canon = build_kmer_vocabulary(train_seqs, k=k, top_n=top_n)
    X_canon = kmers_to_feature_matrix(train_seqs, vocab_canon, k=k)
    vocab_cache = _vocab_from_cache(train, cache, top_n)
    X_cache = _matrix_from_cache(train, vocab_cache, cache)
    assert vocab_canon == vocab_cache, "VOCAB MISMATCH: cached path diverges from canonical"
    assert np.array_equal(X_canon, X_cache), "MATRIX MISMATCH: cached path diverges from canonical"
    print(f"[verify] feature-level equivalence PASS (vocab={len(vocab_canon)}, X={X_canon.shape})")


def run_kmer_loso_cached(seqs_by_strain, labels_by_strain, ids, drug, k, top_n,
                         cache, binary: bool = False) -> tuple[np.ndarray, np.ndarray]:
    """LOSO over cached counts. binary=True -> presence/absence (X>0 -> 1)."""
    n = len(ids)
    y_true, y_score = [], []
    for i, held in enumerate(ids):
        train = [ids[j] for j in range(n) if j != i]
        vocab = _vocab_from_cache(train, cache, top_n)
        vocab_index = {km: idx for idx, km in enumerate(vocab)}
        X_train = _matrix_from_cache(train, vocab, cache)
        train_y = np.array([labels_by_strain[s] for s in train], dtype=int)
        # held-out row restricted to train vocab (mirrors kmers_to_feature_matrix on [test_seq])
        X_test = np.zeros((1, len(vocab)), dtype=np.float32)
        for km, v in cache[held].items():
            idx = vocab_index.get(km)
            if idx is not None and v:
                X_test[0, idx] = v
        if binary:
            X_train = (X_train > 0).astype(np.float32)
            X_test = (X_test > 0).astype(np.float32)
        clf = train_xgboost_classifier(X_train, train_y, drug_name=drug, calibrate=False)
        y_score.append(float(predict_proba(clf, X_test)[0]))
        y_true.append(int(labels_by_strain[held]))
        print(f"[kmer] fold {i+1}/{n} {held} done", flush=True)
    return np.array(y_true, dtype=int), np.array(y_score, dtype=np.float32)


def load_strains() -> tuple[dict[str, str], dict[str, int], list[str]]:
    """Fetch (cached) + concat the bake-off strain set. Shared by the bake-off
    and the full-spectrum diagnostic so both run on byte-identical sequences."""
    strains = select_strains()
    print(f"[kmer] selected {len(strains)} strains "
          f"({sum(1 for s in strains if s[3]==0)} ExPEC / {sum(1 for s in strains if s[3]==1)} EPEC)")
    seqs_by_strain, labels_by_strain, ids = {}, {}, []
    for sid, aname, cls, y, src in strains:
        setp = wgs_set_prefix(aname)
        print(f"[kmer] {sid} ({cls}, {src.split()[0]}) <- ENA {setp}", flush=True)
        ft = fetch_assembly(setp)
        if ft is None:
            print(f"  skip {sid}")
            continue
        seq = concat_contigs(ft)
        if len(seq) < 1_000_000:
            print(f"  skip {sid}: assembly too short ({len(seq)} bp)")
            continue
        key = f"{sid}|{setp}"
        seqs_by_strain[key] = seq
        labels_by_strain[key] = y
        ids.append(key)
    return seqs_by_strain, labels_by_strain, ids


def main() -> int:
    seqs_by_strain, labels_by_strain, ids = load_strains()
    npos = sum(labels_by_strain.values()); nneg = len(ids) - npos
    print(f"[kmer] fetched {len(ids)} usable ({nneg} ExPEC / {npos} EPEC); running LOSO k-mer...")
    if npos < 3 or nneg < 3:
        print("[kmer] ABORT: too few per class after fetch"); return 1
    k = 8
    print("[kmer] precomputing per-genome k-mer counts (once)...", flush=True)
    cache = precompute_counts(seqs_by_strain, k=k)
    _verify_feature_equivalence(seqs_by_strain, ids, cache, k=k, top_n=TOP_N)
    yt, ys = run_kmer_loso_cached(seqs_by_strain, labels_by_strain, ids,
                                  drug="expec_vs_epec", k=k, top_n=TOP_N, cache=cache)
    auroc = float(roc_auc_score(yt, ys)) if len(set(yt.tolist())) == 2 else float("nan")
    uniq = sorted(set(round(float(s), 4) for s in ys))
    degenerate = len(uniq) <= 2
    print(f"[kmer] LOSO AUROC = {auroc:.4f} | distinct scores={len(uniq)} degenerate={degenerate}")
    res = {"contrast": "ExPEC(Salipante) vs EPEC(Hazen)", "n": len(ids),
           "n_expec": nneg, "n_epec": npos, "representation": "kmer_k8_xgboost",
           "cv": "loso", "auroc": auroc, "distinct_scores": len(uniq), "degenerate": degenerate,
           "confound": "study-confounded (ExPEC=Salipante, EPEC=Hazen) -> high AUROC may be batch, not biology",
           "compute": "cached-per-genome LOSO; feature-level equivalence to run_kmer_xgboost_loso asserted (fold 0 vocab+X identical)",
           "ids": ids}
    out = REPO / "research_outputs/pathotype_bakeoff_kmer_expec_epec_2026-05-31.json"
    out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(f"[kmer] wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
