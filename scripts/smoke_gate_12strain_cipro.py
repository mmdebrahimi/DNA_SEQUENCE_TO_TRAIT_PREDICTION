"""Smoke gate runner for the 12-strain cipro mini cohort.

3 variants (NT-XGBoost + k-mer + gene-presence); AMRFinder deferred to Stage 1
since the cohort's persisted `plasmid_resistance_genes` / `chromosome_resistance_genes`
fields are empty for the mini cohort and no per-strain AMRFinderPlus CLI
infrastructure exists yet.

This is an ENGINEERING SMOKE / FALSIFICATION gate, NOT a powered classifier
comparison. At N=12 with balanced LOSO the per-strain noise floor is ±8.3%
and 95% AUROC CI width is ~±0.19. The 15-percentage-point gap acceptance bar
is an engineering heuristic to catch "NT obviously broken on real data" —
NOT statistically powered ranking. Real decision gate is Stage 1 N=50 → Stage 2 N=150.
"""
from __future__ import annotations

import argparse
import sys
from datetime import date
from pathlib import Path

import numpy as np

from dna_decode.data.annotations import parse_gff3
from dna_decode.data.cohort import load_cohort
from dna_decode.data.refseq import fasta_path, gff_path
from dna_decode.eval.cv import leave_one_strain_out_cv
from dna_decode.eval.metrics import compute_metrics
from dna_decode.models.cache import EmbeddingCache
from dna_decode.models.classical_baselines import (
    build_gene_presence_matrix,
    build_kmer_vocabulary,
    kmers_to_feature_matrix,
)
from dna_decode.models.classifiers import (
    ClassifierTrainingError,
    aggregate_strain_features,
    predict_proba,
    train_xgboost_classifier,
)


def _load_fasta_contigs(fasta_p: Path) -> list[str]:
    """Load all contigs as uppercase strings."""
    from Bio import SeqIO

    return [str(record.seq).upper() for record in SeqIO.parse(str(fasta_p), "fasta")]


def _extract_gene_ids(annotations_table) -> set[str]:
    """Extract unique gene IDs from a parsed GFF3 table (CDS rows only)."""
    cds = annotations_table[annotations_table["type"] == "CDS"]
    gene_ids: set[str] = set()
    for _, row in cds.iterrows():
        gid = row.get("gene_id") or row.get("locus_tag")
        if gid:
            gene_ids.add(str(gid))
    return gene_ids


def run_nt_xgboost(cohort, cache_path: Path, drug: str) -> dict:
    """NT-XGBoost LOSO using the existing populated cache."""
    cache = EmbeddingCache(
        cache_path,
        model_name="nucleotide_transformer",
        model_version="InstaDeepAI/nucleotide-transformer-v2-100m-multi-species",
        embedding_dim=512,
    )
    drug_lower = drug.lower()
    drug_strain_ids = cohort.per_drug_strain_ids[drug_lower]
    X_rows = []
    labels = []
    strain_order = []
    for sid in drug_strain_ids:
        s = cohort.strain_by_id(sid)
        if s is None or drug_lower not in s.ast_labels:
            continue
        gene_ids = cache.list_genes(sid)
        if not gene_ids:
            continue
        gene_matrix = cache.bulk_get([(sid, g) for g in gene_ids])
        X_rows.append(aggregate_strain_features(gene_matrix, "mean"))
        labels.append(s.ast_labels[drug_lower])
        strain_order.append(sid)
    X = np.stack(X_rows)
    y = np.array(labels, dtype=int)

    def _train(X_tr, y_tr):
        # calibrate=False for smoke at N=12. At N=11 training, CalibratedClassifierCV's
        # isotonic regression overcorrects and inverts predictions (AUROC = 0.000 vs
        # 0.750 without calibration on this cohort). Verified 2026-05-14.
        # Re-enable calibration at Stage 1 N=50+ where the calibration CV has enough
        # samples per fold to be stable.
        return train_xgboost_classifier(X_tr, y_tr, drug_name=drug, calibrate=False)

    def _predict(model, X_te):
        return predict_proba(model, X_te)

    cv = leave_one_strain_out_cv(X, y, strain_order, _train, _predict, drug=drug)
    m = compute_metrics(cv.all_y_true, cv.all_y_score)
    return {
        "variant": "NT-XGBoost (nucleotide_transformer)",
        "auroc": float(m.auroc),
        "auprc": float(m.auprc),
        "n_samples": int(m.n_samples),
        "label_balance": f"{int((y == 1).sum())}R / {int((y == 0).sum())}S",
    }


def run_kmer_xgboost(cohort, refseq_root: Path, drug: str, k: int = 8, top_n: int = 10000) -> dict:
    """k-mer + XGBoost LOSO with within-fold vocabulary rebuild."""
    drug_lower = drug.lower()
    strain_contigs: dict[str, list[str]] = {}
    labels_dict: dict[str, int] = {}
    for s in cohort.strains:
        if drug_lower not in s.ast_labels:
            continue
        fp = fasta_path(s.assembly_accession, refseq_root)
        strain_contigs[s.strain_id] = _load_fasta_contigs(fp)
        labels_dict[s.strain_id] = s.ast_labels[drug_lower]

    strain_order = sorted(strain_contigs.keys())
    sep = "N" * 100
    concatenated = [sep.join(strain_contigs[sid]) for sid in strain_order]
    y = np.array([labels_dict[sid] for sid in strain_order], dtype=int)

    n = len(strain_order)
    all_y_true: list[int] = []
    all_y_score: list[float] = []
    for i in range(n):
        train_seqs = [concatenated[j] for j in range(n) if j != i]
        train_y = np.array([y[j] for j in range(n) if j != i], dtype=int)
        test_seq = concatenated[i]
        test_y = int(y[i])

        vocab = build_kmer_vocabulary(train_seqs, k=k, top_n=top_n)
        X_train = kmers_to_feature_matrix(train_seqs, vocab, k=k)
        X_test = kmers_to_feature_matrix([test_seq], vocab, k=k)

        try:
            clf = train_xgboost_classifier(X_train, train_y, drug_name=drug, calibrate=False)
            y_score = float(predict_proba(clf, X_test)[0])
        except ClassifierTrainingError:
            y_score = float(train_y.mean())

        all_y_true.append(test_y)
        all_y_score.append(y_score)

    m = compute_metrics(np.array(all_y_true), np.array(all_y_score))
    return {
        "variant": f"k-mer (k={k}) + XGBoost",
        "auroc": float(m.auroc),
        "auprc": float(m.auprc),
        "n_samples": int(m.n_samples),
        "label_balance": f"{int((y == 1).sum())}R / {int((y == 0).sum())}S",
    }


def run_gene_presence_xgboost(cohort, refseq_root: Path, drug: str) -> dict:
    """Gene-presence + XGBoost LOSO with within-fold vocabulary rebuild."""
    drug_lower = drug.lower()
    strain_genes: dict[str, set[str]] = {}
    labels_dict: dict[str, int] = {}
    for s in cohort.strains:
        if drug_lower not in s.ast_labels:
            continue
        gp = gff_path(s.assembly_accession, refseq_root)
        ann = parse_gff3(gp)
        strain_genes[s.strain_id] = _extract_gene_ids(ann)
        labels_dict[s.strain_id] = s.ast_labels[drug_lower]

    strain_order = sorted(strain_genes.keys())
    y = np.array([labels_dict[sid] for sid in strain_order], dtype=int)

    n = len(strain_order)
    all_y_true: list[int] = []
    all_y_score: list[float] = []
    for i in range(n):
        train_gene_sets = [strain_genes[strain_order[j]] for j in range(n) if j != i]
        train_y = np.array([y[j] for j in range(n) if j != i], dtype=int)
        test_gene_set = strain_genes[strain_order[i]]
        test_y = int(y[i])

        X_train, vocab = build_gene_presence_matrix(train_gene_sets, gene_vocabulary=None)
        X_test, _ = build_gene_presence_matrix([test_gene_set], gene_vocabulary=vocab)

        try:
            clf = train_xgboost_classifier(X_train, train_y, drug_name=drug, calibrate=False)
            y_score = float(predict_proba(clf, X_test)[0])
        except ClassifierTrainingError:
            y_score = float(train_y.mean())

        all_y_true.append(test_y)
        all_y_score.append(y_score)

    m = compute_metrics(np.array(all_y_true), np.array(all_y_score))
    return {
        "variant": "Gene-presence + XGBoost",
        "auroc": float(m.auroc),
        "auprc": float(m.auprc),
        "n_samples": int(m.n_samples),
        "label_balance": f"{int((y == 1).sum())}R / {int((y == 0).sum())}S",
    }


def write_packet(
    results: list[dict],
    output_path: Path,
    cohort_path: Path,
    drug: str,
    gap_threshold_pp: float = 15.0,
) -> dict:
    """Write the smoke result packet + compute acceptance verdict."""
    nt = next((r for r in results if r["variant"].startswith("NT-")), None)
    classical = [r for r in results if not r["variant"].startswith("NT-")]
    best_classical = max(classical, key=lambda r: r["auroc"]) if classical else None

    if nt and best_classical:
        gap_pp = (best_classical["auroc"] - nt["auroc"]) * 100
        obviously_worse = gap_pp >= gap_threshold_pp
        verdict = "FAIL (NT obviously worse)" if obviously_worse else "PASS (NT not obviously worse)"
    else:
        gap_pp = None
        verdict = "INDETERMINATE (missing variant)"

    today = date.today().isoformat()
    lines = [
        f"# Smoke Gate — 12-strain cipro cohort ({today})",
        "",
        "> **This is an engineering smoke / falsification gate, NOT a powered classifier comparison.**",
        "> At N=12 with balanced LOSO the per-strain noise floor is ±8.3% and 95% AUROC CI width is ~±0.19.",
        "> The 15-percentage-point gap acceptance bar is an engineering heuristic to catch \"NT obviously broken\" —",
        "> NOT statistically powered ranking. Multiple-comparison statistical power at N=12 forbids classifier ranking.",
        "> Real decision gate scheduled at Stage 1 N=50 (local engineering screen) → Stage 2 N=150 (Databricks).",
        "",
        f"**Cohort:** `{cohort_path}` (12 strains, 6R/6S cipro)",
        f"**Drug:** {drug}",
        f"**Gap threshold:** ≥{gap_threshold_pp:.0f} pp (NT-XGBoost AUROC ≥ best-classical AUROC − {gap_threshold_pp:.0f} pp)",
        f"**Verdict:** {verdict}",
        "",
        "## Per-variant LOSO results",
        "",
        "| Variant | AUROC | AUPRC | N | Label balance |",
        "|---|---:|---:|---:|---|",
    ]
    for r in results:
        lines.append(
            f"| {r['variant']} | {r['auroc']:.3f} | {r['auprc']:.3f} | {r['n_samples']} | {r['label_balance']} |"
        )

    lines.extend(["", "## Gap analysis", ""])
    if nt and best_classical:
        lines.append(f"- NT-XGBoost AUROC: **{nt['auroc']:.3f}**")
        lines.append(f"- Best classical baseline ({best_classical['variant']}): **{best_classical['auroc']:.3f}**")
        lines.append(f"- Gap: **{gap_pp:+.1f} pp** (best_classical − NT)")
        lines.append(
            f"- Acceptance bar: gap < {gap_threshold_pp:.0f} pp → "
            f"{'PASS' if not obviously_worse else 'FAIL'}"
        )
    else:
        lines.append("- Missing variants; gap analysis skipped.")

    lines.extend([
        "",
        "## Notes",
        "",
        "- 3 variants run (NT-XGBoost + k-mer + gene-presence). AMRFinder deferred to Stage 1 prep — cohort's persisted `plasmid_resistance_genes` / `chromosome_resistance_genes` fields are empty for the mini cohort, and no per-strain AMRFinderPlus CLI infrastructure exists yet.",
        "- k-mer and gene-presence both use within-fold vocabulary rebuild (training-set-only) to prevent held-out leakage.",
        "- Top-K attribution genes (Tier 1-5 classification; gyrA/parC/parE presence check) NOT included in this smoke. Run as a separate `pipeline.py attribute` step if smoke passes.",
        "- Locked decisions reflected: B-B (smoke = 4→3 variants; clade-only dropped), `--per-class 20` for Stage 1 cohort (N=40), deterministic `zlib.crc32` for any future MLST hashing.",
        "",
        "## Next action",
        "",
        "- **PASS** → proceed to Stage 1 N=50 local engineering screen (~4 hours of GTX 860M time).",
        "- **FAIL** → NT is broken on real data. Demote NT track. Classical baselines become project spine.",
        "- **INDETERMINATE** → fix missing variant or rerun smoke.",
    ])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines), encoding="utf-8")
    return {"verdict": verdict, "gap_pp": gap_pp, "nt": nt, "best_classical": best_classical}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="12-strain cipro smoke gate runner (3-variant)")
    parser.add_argument("--cohort", type=Path, default=Path("data/processed/gate_b_mini_cohort.parquet"))
    parser.add_argument("--nt-cache", type=Path, default=Path("data/processed/mini_cipro_nt_cache.h5"))
    parser.add_argument("--refseq-cache", type=Path, default=Path("F:/dna_decode_cache/refseq"))
    parser.add_argument("--drug", default="ciprofloxacin")
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args(argv)

    if args.output is None:
        args.output = Path(f"wiki/smoke_gate_12strain_cipro_{date.today().isoformat()}.md")

    cohort = load_cohort(args.cohort)
    print(f"[smoke] cohort: {len(cohort.strains)} strains; drug={args.drug}")

    results: list[dict] = []

    print("[smoke] Running NT-XGBoost LOSO...")
    try:
        r = run_nt_xgboost(cohort, args.nt_cache, args.drug)
        results.append(r)
        print(f"  AUROC={r['auroc']:.3f} AUPRC={r['auprc']:.3f}")
    except Exception as e:
        print(f"  FAILED: {type(e).__name__}: {e}", file=sys.stderr)

    print("[smoke] Running k-mer + XGBoost LOSO...")
    try:
        r = run_kmer_xgboost(cohort, args.refseq_cache, args.drug)
        results.append(r)
        print(f"  AUROC={r['auroc']:.3f} AUPRC={r['auprc']:.3f}")
    except Exception as e:
        print(f"  FAILED: {type(e).__name__}: {e}", file=sys.stderr)

    print("[smoke] Running gene-presence + XGBoost LOSO...")
    try:
        r = run_gene_presence_xgboost(cohort, args.refseq_cache, args.drug)
        results.append(r)
        print(f"  AUROC={r['auroc']:.3f} AUPRC={r['auprc']:.3f}")
    except Exception as e:
        print(f"  FAILED: {type(e).__name__}: {e}", file=sys.stderr)

    print(f"[smoke] Writing result packet to {args.output}")
    summary = write_packet(results, args.output, args.cohort, args.drug)
    print(f"[smoke] {summary['verdict']}")
    if summary["gap_pp"] is not None:
        print(
            f"[smoke] gap: {summary['gap_pp']:+.1f} pp "
            f"(best classical: {summary['best_classical']['variant']})"
        )
    return 0 if summary["verdict"].startswith("PASS") else 1


if __name__ == "__main__":
    sys.exit(main())
