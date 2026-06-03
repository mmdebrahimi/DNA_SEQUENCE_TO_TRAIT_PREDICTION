"""ExPEC-vs-EPEC windowed-NT arm of the representation bake-off (EP-4).

Whole-genome WINDOWED NT mean-pool (NO Bakta annotation -> avoids the ~10min/genome
annotation bottleneck). Same 24-strain cohort + same LOSO as the k-mer arm, so the
AUROCs are directly comparable. The project's mandated bar (CLAUDE.md:211): the
foundation model must beat the classical (k-mer) baseline by >=3pp.

k-mer baseline on this contrast = 0.514 (chance, non-degenerate) -> beatable, max headroom.
Reuses cached ENA assemblies in data/ena_wgs/ (shared with the k-mer arm).
"""
from __future__ import annotations
import csv, re, sys, json
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
import numpy as np
from sklearn.metrics import roc_auc_score
from dna_decode.models.foundation import model_factory

F1 = REPO / "data/external/horesh2021_F1_genome_metadata.csv"
ENA_CACHE = REPO / "data/ena_wgs"
N_PER_CLASS = 12
WINDOW = 6000        # ~half NT v2 max_context (12288 nt); one embed call per window
N_WINDOWS = 200      # subsample windows/strain (matches smoke budget; keeps 860M fast)
WGS_MASTER = re.compile(r"^[A-Z]{4}\d{8}(\.fa)?$", re.I)


def select_strains():
    rows = list(csv.DictReader(open(F1, encoding="utf-8")))
    clean = [r for r in rows if "(predicted)" not in r["Pathotype"] and r["Pathotype"] not in ("Not determined", "")]
    ok = lambda r: bool(WGS_MASTER.match(r["Assembly_name"].strip()))
    expec = [r for r in clean if r["Pathotype"].strip().startswith("ExPEC") and r["Source"].startswith("Salipante") and ok(r)][:N_PER_CLASS]
    epec = [r for r in clean if r["Pathotype"].strip().startswith("EPEC") and r["Source"].startswith("Hazen") and ok(r)][:N_PER_CLASS]
    out = []
    for r in expec: out.append((r["ID"], r["Assembly_name"][:4] + "01", 0))
    for r in epec: out.append((r["ID"], r["Assembly_name"][:4] + "01", 1))
    return out


def load_concat(setp: str) -> str | None:
    p = ENA_CACHE / f"{setp}.fna"
    if not p.exists():
        return None
    seqs, cur = [], []
    for line in p.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if cur: seqs.append("".join(cur)); cur = []
        else:
            cur.append(line.strip())
    if cur: seqs.append("".join(cur))
    return "".join(seqs)  # contigs concatenated; windows tile the whole genome


def windowed_embedding(model, genome: str, rng) -> np.ndarray | None:
    n = len(genome)
    if n < WINDOW:
        return None
    starts = list(range(0, n - WINDOW, WINDOW))
    if not starts:
        return None
    if len(starts) > N_WINDOWS:
        starts = sorted(rng.choice(np.array(starts), size=N_WINDOWS, replace=False).tolist())
    vecs = []
    for s in starts:
        w = genome[s:s + WINDOW]
        if w.count("N") > 0.5 * len(w):   # skip mostly-N windows
            continue
        v = model.embed(w)
        v = np.asarray(v)
        if v.ndim == 2:
            v = v.mean(axis=0)
        vecs.append(v.astype(np.float32))
    return np.mean(np.stack(vecs), axis=0) if vecs else None


def loso_auroc(X, y):
    from xgboost import XGBClassifier
    n = len(y); proba = np.zeros(n)
    for i in range(n):
        tr = [j for j in range(n) if j != i]
        clf = XGBClassifier(n_estimators=120, max_depth=3, learning_rate=0.1,
                            subsample=0.9, eval_metric="logloss", n_jobs=4)
        clf.fit(X[tr], y[tr])
        proba[i] = clf.predict_proba(X[i:i+1])[0, 1]
    return float(roc_auc_score(y, proba)), proba


def main() -> int:
    rng = np.random.default_rng(0)
    strains = select_strains()
    print(f"[nt] {len(strains)} strains; loading NT v2 100M (cuda)...", flush=True)
    model = model_factory("nucleotide_transformer",
                          config_path=str(REPO / "config" / "datasources.yaml"), device="cuda")
    X, y, ids = [], [], []
    for sid, setp, label in strains:
        g = load_concat(setp)
        if g is None:
            print(f"[nt] skip {sid}: no cached {setp}"); continue
        print(f"[nt] {sid} ({setp}) embed {min(N_WINDOWS, len(g)//WINDOW)} windows...", flush=True)
        v = windowed_embedding(model, g, rng)
        if v is None:
            print(f"[nt] skip {sid}: no usable windows"); continue
        X.append(v); y.append(label); ids.append(f"{sid}|{setp}")
    X = np.stack(X); y = np.array(y)
    npos = int(y.sum()); nneg = int((y == 0).sum())
    print(f"[nt] embedded {len(y)} ({nneg} ExPEC / {npos} EPEC); LOSO...", flush=True)
    auroc, proba = loso_auroc(X, y)
    uniq = len(set(round(float(p), 4) for p in proba))
    kmer = 0.5139
    lift = auroc - kmer
    verdict = ("NT beats k-mer by >=3pp (real lift)" if lift >= 0.03
               else "NT does NOT beat k-mer by 3pp -> classical wins / no FM advantage")
    print(f"[nt] NT AUROC={auroc:.4f} | k-mer=0.514 | lift={lift:+.4f} | {verdict}")
    res = {"contrast": "ExPEC(Salipante) vs EPEC(Hazen)", "representation": "windowed_NT_v2_100m_meanpool",
           "window_bp": WINDOW, "n_windows": N_WINDOWS, "n": len(y), "n_expec": nneg, "n_epec": npos,
           "cv": "loso", "nt_auroc": auroc, "kmer_auroc_baseline": kmer, "lift_vs_kmer": lift,
           "distinct_scores": uniq, "verdict": verdict, "bar": ">=3pp over classical (CLAUDE.md:211)",
           "confound_note": "study-confounded contrast; k-mer at chance (0.514) so confound not k-mer-separable; "
                            "if NT lifts, check it's pathotype signal not batch via lineage-aware split next",
           "ids": ids}
    out = REPO / "research_outputs/pathotype_bakeoff_nt_windowed_expec_epec_2026-05-31.json"
    out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(f"[nt] wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
