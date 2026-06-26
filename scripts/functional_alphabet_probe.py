"""Non-neural functional-alphabet probe — does a function-level token alphabet beat k-mer WITHIN-LINEAGE?

The /hypothesise + /probe deliverable: the genomic-embedding bet is a closed 0-for-4 de-confounded
negative (learned lineage, not mechanism). Before paying for any neural/contrastive training, test the
ALPHABET half cheaply, NON-NEURALLY, CPU-only: a FUNCTIONAL-UNIT alphabet (codon/allele/mechanism tokens,
`functional_tokens.py`) vs the base-level K-MER alphabet, both scored on the project's de-confounded
WITHIN-LINEAGE metric under leave-one-MLST-out CV.

KILL CRITERION (pre-committed, frozen in `compute_probe_verdict`): the functional alphabet must EXCEED
the k-mer alphabet on within-lineage concordance, with the paired in-MLST label-permutation null giving
p < 0.05 for the gap. Any of BEATS_KMER / TIES / FAILS / UNDERPOWERED is a valid, honest completion.

sklearn-only (no xgboost / no GPU). Reuses point_baseline (functional tokens) + classical_baselines
(k-mer helpers) + within_lineage_diagnostic (the metric). Exit 0 always -- a research packet, not a gate.
"""
from __future__ import annotations

import argparse
import datetime
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from dna_decode.eval.functional_tokens import build_feature_matrix  # noqa: E402
from dna_decode.models.classical_baselines import (  # noqa: E402
    CONTIG_SEPARATOR,
    build_kmer_vocabulary,
    kmers_to_feature_matrix,
)
from scripts.within_lineage_diagnostic import within_lineage_concordance  # noqa: E402

MIN_PAIRS = 10  # below this the within-lineage metric is too underpowered to call


def _read_genome(refseq_cache: Path, accession: str) -> str | None:
    p = refseq_cache / accession / "genome.fna"
    if not p.exists():
        hits = list(refseq_cache.glob(f"{accession}*/genome.fna")) or list(refseq_cache.glob(f"{accession}*/*.fna"))
        if not hits:
            return None
        p = hits[0]
    contigs = []
    cur = []
    for line in p.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if cur:
                contigs.append("".join(cur))
                cur = []
        else:
            cur.append(line.strip())
    if cur:
        contigs.append("".join(cur))
    return CONTIG_SEPARATOR.join(contigs)


def _group_out_scores(X: np.ndarray, y: np.ndarray, groups: list[str]) -> list[float]:
    """Leave-one-MLST-out held-out scores via sklearn LogisticRegression (group = lineage)."""
    from sklearn.linear_model import LogisticRegression
    scores: list[float | None] = [None] * len(y)
    for g in sorted(set(groups)):
        test = [i for i, gg in enumerate(groups) if gg == g]
        train = [i for i in range(len(y)) if groups[i] != g]
        ytr = y[train]
        if len(set(ytr.tolist())) < 2:  # single-class train fold -> constant prior score
            for i in test:
                scores[i] = float(ytr.mean())
            continue
        clf = LogisticRegression(max_iter=1000, class_weight="balanced")
        clf.fit(X[train], ytr)
        df = np.atleast_1d(clf.decision_function(X[test]))
        for k, i in enumerate(test):
            scores[i] = float(df[k])
    return [s if s is not None else 0.0 for s in scores]


def _paired_within_lineage(func_scores, kmer_scores, y, mlst, *, n_perm=2000, seed=42) -> dict:
    """Within-lineage concordance for each arm + a paired in-MLST label-permutation null on the gap."""
    f_obs, n_pairs, n_shared = within_lineage_concordance(func_scores, y, mlst)
    k_obs, _, _ = within_lineage_concordance(kmer_scores, y, mlst)
    gap = f_obs - k_obs
    out = {"func_wl": f_obs, "kmer_wl": k_obs, "gap": gap, "n_pairs": n_pairs, "n_shared_lineages": n_shared}
    if not np.isfinite(gap):
        out.update({"gap_p": float("nan"), "gap_null_ci": [float("nan"), float("nan")]})
        return out
    rng = np.random.default_rng(seed)
    by = defaultdict(list)
    for i, m in enumerate(mlst):
        if m not in (None, "None", ""):
            by[m].append(i)
    shared = [idxs for idxs in by.values()
              if any(y[i] == 1 for i in idxs) and any(y[i] == 0 for i in idxs)]
    null = []
    for _ in range(n_perm):
        yp = list(y)
        for idxs in shared:
            lab = [y[i] for i in idxs]
            rng.shuffle(lab)
            for i, lv in zip(idxs, lab):
                yp[i] = lv
        fc, _, _ = within_lineage_concordance(func_scores, yp, mlst)
        kc, _, _ = within_lineage_concordance(kmer_scores, yp, mlst)
        null.append(fc - kc)
    null = np.asarray(null)
    out["gap_p"] = float((null >= gap).mean())  # one-sided: functional beats k-mer beyond the lineage null
    out["gap_null_ci"] = [float(np.percentile(null, 2.5)), float(np.percentile(null, 97.5))]
    return out


def compute_probe_verdict(func_wl, kmer_wl, gap, gap_p, n_pairs, *, min_pairs=MIN_PAIRS) -> dict:
    """FROZEN kill-criterion (pure). UNDERPOWERED if too few pairs / non-finite; else sign+significance."""
    finite = all(np.isfinite(v) for v in (func_wl, kmer_wl, gap))
    if (n_pairs or 0) < min_pairs or not finite:
        verdict = "UNDERPOWERED"
    elif gap > 0 and np.isfinite(gap_p) and gap_p < 0.05:
        verdict = "BEATS_KMER"
    elif gap > 0:
        verdict = "TIES"          # functional higher but not beyond the lineage-conditioned null
    else:
        verdict = "FAILS"         # functional <= k-mer within-lineage
    return {"verdict": verdict, "func_wl": func_wl, "kmer_wl": kmer_wl, "gap": gap,
            "gap_p": gap_p, "n_pairs": n_pairs}


def run_probe(cohort_path: Path, runs_root: Path, refseq_cache: Path, drug: str) -> dict:
    df = pd.read_parquet(cohort_path)
    label_col = f"ast_{drug}"
    df = df[df[label_col].notna()].copy()
    df[label_col] = df[label_col].astype(int)

    # FUNCTIONAL arm matrix (drops strains with no AMRFinder cache)
    Xf, vocab_f, func_ids, dropped_amr = build_feature_matrix(df, runs_root, drug)
    # K-MER arm: load genomes for the functional-admitted strains; intersect on genome availability
    acc_by_strain = dict(zip(df["strain_id"].astype(str), df["assembly_accession"].astype(str)))
    seqs: dict[str, str] = {}
    dropped_fasta = []
    for s in func_ids:
        g = _read_genome(refseq_cache, acc_by_strain[s])
        if g is None:
            dropped_fasta.append(s)
        else:
            seqs[s] = g
    admitted = [s for s in func_ids if s in seqs]           # SAME set, SAME order for both arms
    keep = [i for i, s in enumerate(func_ids) if s in seqs]
    Xf = Xf[keep]

    seq_list = [seqs[s] for s in admitted]
    kvocab = build_kmer_vocabulary(seq_list)
    Xk = kmers_to_feature_matrix(seq_list, kvocab)

    lab = dict(zip(df["strain_id"].astype(str), df[label_col]))
    mlst = dict(zip(df["strain_id"].astype(str), df["mlst"].astype(str)))
    y = np.asarray([int(lab[s]) for s in admitted])
    groups = [mlst[s] for s in admitted]

    assert Xf.shape[0] == Xk.shape[0] == len(admitted) == len(y), "functional/k-mer arms misaligned"

    func_scores = _group_out_scores(Xf, y, groups)
    kmer_scores = _group_out_scores(Xk, y, groups)
    metric = _paired_within_lineage(func_scores, kmer_scores, y, mlst=groups)
    verdict = compute_probe_verdict(metric["func_wl"], metric["kmer_wl"], metric["gap"],
                                    metric["gap_p"], metric["n_pairs"])

    return {
        "schema": "functional-alphabet-probe-v0",
        "cohort": cohort_path.name, "drug": drug,
        "n_admitted": len(admitted), "n_R": int((y == 1).sum()), "n_S": int((y == 0).sum()),
        "dropped_no_amrfinder": dropped_amr, "dropped_no_fasta": dropped_fasta,
        "functional_vocab_size": len(vocab_f), "kmer_vocab_size": len(kvocab),
        "metric": metric, "verdict": verdict,
        "kill_criterion": ("functional within-lineage concordance EXCEEDS k-mer with paired in-MLST "
                           "label-permutation gap p < 0.05; UNDERPOWERED when n_pairs < "
                           f"{MIN_PAIRS}. Any bucket is a valid honest completion."),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Non-neural functional-alphabet probe (within-lineage)")
    ap.add_argument("--cohort", type=Path, default=ROOT / "data/processed/stage2_n150_cipro_cohort.parquet")
    ap.add_argument("--smoke", action="store_true", help="use the N=40 gate_b mini cohort")
    ap.add_argument("--runs-root", type=Path, default=ROOT / "data/amrfinder_runs")
    ap.add_argument("--refseq-cache", type=Path, default=ROOT / "data/refseq_cache")
    ap.add_argument("--drug", default="ciprofloxacin")
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args(argv)

    cohort = (ROOT / "data/processed/gate_b_n40_cipro_cohort.parquet") if args.smoke else args.cohort
    rep = run_probe(cohort, args.runs_root, args.refseq_cache, args.drug)
    rep["analysis_date"] = datetime.date.today().isoformat()
    tag = "smoke" if args.smoke else "n147"
    base = args.out or (ROOT / "wiki" / f"functional_alphabet_probe_{tag}_{rep['analysis_date']}")
    Path(f"{base}.json").write_text(json.dumps(rep, indent=2), encoding="utf-8")

    m, v = rep["metric"], rep["verdict"]
    L = [f"# Functional-alphabet probe ({rep['analysis_date']}) -- {rep['cohort']}", "",
         f"**Verdict: {v['verdict']}**", "",
         f"- functional within-lineage concordance: {m['func_wl']}",
         f"- k-mer within-lineage concordance:      {m['kmer_wl']}",
         f"- gap (functional - k-mer):              {m['gap']}",
         f"- paired in-MLST permutation gap p:      {m.get('gap_p')}",
         f"- gap null 95% band:                     {m.get('gap_null_ci')}",
         f"- powering: n_shared_lineages={m['n_shared_lineages']}, n_pairs={m['n_pairs']}",
         f"- admitted N={rep['n_admitted']} (R={rep['n_R']}/S={rep['n_S']}); "
         f"dropped: {len(rep['dropped_no_amrfinder'])} no-AMRFinder, {len(rep['dropped_no_fasta'])} no-FASTA",
         f"- vocab: functional={rep['functional_vocab_size']} tokens, k-mer={rep['kmer_vocab_size']}", "",
         f"_{rep['kill_criterion']}_", "",
         "_Non-neural, CPU-only. A within-lineage win means 'beats k-mer within-lineage', NOT 'mechanism proven'._"]
    Path(f"{base}.md").write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L))
    print(f"[probe -> {base}.{{md,json}}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
