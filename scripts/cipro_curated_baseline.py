"""Cipro curated AMR baseline — Experiment 3 (PIVOT TRIGGER condition 4).

Builds a classical AMR baseline from AMRFinderPlus mechanism-audit outputs +
k-mer features + MLST one-hot. Trains LR (scaled) and XGB heads under LOSO
on N=38. Compares against the Stage 1b NT verdicts:

  Stage 1b NT-LR  AUROC = 0.673
  Stage 1b NT-XGB AUROC = 0.615
  Stage 1b k-mer-XGB    = 0.648

PIVOT TRIGGER condition 4 (declares frozen NT whole-genome pooling falsified):
  curated AMR baseline >= 0.80 AUROC at N=38  OR  beats NT by >= 10 pp.

Features (4 blocks):
  - POINT mutations:      multi-hot over union of QUINOLONE-class POINT mutation
                          symbols seen across cohort (gyrA_S83L, parC_S80I, etc.)
  - acquired AMR genes:   multi-hot over union of acquired gene symbols
                          (qnrS1, aac(6')-Ib-cr, blaCMY-2, sul1, ...)
  - k-mer top-N:          k=8 top-10000 multi-hot per strain (matching Stage 1b
                          k-mer-XGB feature space; rebuilt PER FOLD to prevent
                          test-leak into vocabulary selection)
  - MLST one-hot:         from cohort metadata

NOTE: Bakta gene presence is in the PIVOT TRIGGER spec but Bakta annotation
has not yet been run cohort-wide (toolchain validated but only on smoke
substrate). Bakta features deferred to a separate `--include-bakta` flag.

Run AFTER cipro_mechanism_audit.py lands (needs the JSON for AMRFinder features).
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import date as _date
from pathlib import Path
from typing import Literal

import numpy as np

from dna_decode.data.cohort import load_cohort
from dna_decode.data.refseq import fasta_path
from dna_decode.eval.cv import CVResult, FoldResult
from dna_decode.eval.metrics import compute_metrics
from dna_decode.models.classical_baselines import (
    CONTIG_SEPARATOR,
    _train_baseline_logreg,
    build_kmer_vocabulary,
    kmers_to_feature_matrix,
)
from dna_decode.models.classifiers import (
    ClassifierTrainingError,
    predict_proba,
    train_xgboost_classifier,
)


BASELINE_ABSOLUTE_GATE_AUROC = 0.80
BASELINE_COMPARATIVE_GATE_PP = 10.0

# Stage 1b baselines (frozen reference points for verdict gate)
STAGE1B_NT_LR_AUROC = 0.673
STAGE1B_NT_XGB_AUROC = 0.615
STAGE1B_KMER_XGB_AUROC = 0.648

# Amended PIVOT TRIGGER condition 4 gates (post-Experiment-2 NOISY verdict +
# Codex round-2 critique 2026-05-17): the all-feature gate is structurally
# circular because POINT mutations are essentially labels-in-genome-form;
# the load-bearing gate is no-POINT (which isolates non-textbook signal)
# OR mechanism-only at the absolute threshold (which proves AMRFinder's
# curated mechanism panel alone matches Stage 2 burst threshold).
AMENDED_NO_POINT_GATE_AUROC = max(0.75, STAGE1B_NT_LR_AUROC + 0.10)  # >= 0.773
AMENDED_MECHANISM_ONLY_GATE_AUROC = BASELINE_ABSOLUTE_GATE_AUROC  # 0.80


# Multi-block feature-set names. Each value is the tuple of blocks passed
# to run_loso. Used by the ablation panel.
ABLATION_FEATURE_SETS: dict[str, tuple[str, ...]] = {
    "all": ("point", "acquired", "kmer", "mlst"),
    "no_POINT": ("acquired", "kmer", "mlst"),
    "mechanism_only": ("point", "acquired"),
    "POINT_only": ("point",),
    "kmer_only": ("kmer",),
    "MLST_only": ("mlst",),
    "kmer_MLST_only": ("kmer", "mlst"),
}


@dataclass
class FeatureBundle:
    """Per-strain feature components + meta. Vocabularies are derived from the
    full cohort for AMR/MLST blocks (AMR features are intrinsic per-strain;
    leakage risk only applies to k-mer where train-fold vocab is rebuilt)."""

    point_mut_vocab: list[str]
    acquired_gene_vocab: list[str]
    mlst_vocab: list[str]
    point_mut: dict[str, np.ndarray]
    acquired: dict[str, np.ndarray]
    mlst: dict[str, np.ndarray]
    seqs: dict[str, str]
    labels: dict[str, int]
    strain_ids: list[str]


def _load_amrfinder_features(audit_json: Path) -> tuple[dict[str, set[str]], dict[str, set[str]]]:
    """Return (per-strain POINT mutation symbols, per-strain acquired gene symbols).

    Filters to QUINOLONE-class POINTs (cipro-relevant) + ALL acquired-gene hits
    (acquired-gene class is broader — qnr* + aac variants live there, plus
    co-resistance markers that may correlate with cipro phenotype).
    """
    data = json.loads(audit_json.read_text(encoding="utf-8"))
    point_by_strain: dict[str, set[str]] = {}
    acquired_by_strain: dict[str, set[str]] = {}
    for entry in data["per_strain"]:
        if entry.get("status") != "OK":
            continue
        sid = entry["strain_id"]
        points: set[str] = set()
        acquired: set[str] = set()
        for h in entry.get("hits", []):
            if h["kind"] == "mutation" and h.get("class", "").upper() in {"QUINOLONE", "FLUOROQUINOLONE"}:
                points.add(h["symbol"])
            elif h["kind"] == "acquired":
                acquired.add(h["symbol"])
        point_by_strain[sid] = points
        acquired_by_strain[sid] = acquired
    return point_by_strain, acquired_by_strain


def _load_fasta_contigs(fasta_p: Path) -> list[str]:
    from Bio import SeqIO
    return [str(record.seq).upper() for record in SeqIO.parse(str(fasta_p), "fasta")]


def build_features(
    cohort,
    audit_json: Path,
    refseq_root: Path,
    drug: str,
    decisive_subset: set[str] | None = None,
) -> FeatureBundle:
    point_by_strain, acquired_by_strain = _load_amrfinder_features(audit_json)
    drug_lower = drug.lower()

    seqs: dict[str, str] = {}
    labels: dict[str, int] = {}
    mlst_by_strain: dict[str, str] = {}
    strain_ids: list[str] = []
    for s in cohort.strains:
        if drug_lower not in s.ast_labels:
            continue
        if s.strain_id not in point_by_strain:
            continue
        if decisive_subset is not None and s.strain_id not in decisive_subset:
            continue
        fp = fasta_path(s.assembly_accession, refseq_root)
        if not fp.exists():
            continue
        seqs[s.strain_id] = CONTIG_SEPARATOR.join(_load_fasta_contigs(fp))
        labels[s.strain_id] = int(s.ast_labels[drug_lower])
        mlst_by_strain[s.strain_id] = s.mlst or "unknown"
        strain_ids.append(s.strain_id)

    # Cohort-wide vocabularies
    point_vocab = sorted({m for sid in strain_ids for m in point_by_strain[sid]})
    acquired_vocab = sorted({g for sid in strain_ids for g in acquired_by_strain[sid]})
    mlst_vocab = sorted(set(mlst_by_strain.values()))

    # Multi-hot per-strain matrices
    point_mat = {sid: np.array([1 if m in point_by_strain[sid] else 0 for m in point_vocab], dtype=np.float32) for sid in strain_ids}
    acquired_mat = {sid: np.array([1 if g in acquired_by_strain[sid] else 0 for g in acquired_vocab], dtype=np.float32) for sid in strain_ids}
    mlst_mat = {sid: np.array([1 if m == mlst_by_strain[sid] else 0 for m in mlst_vocab], dtype=np.float32) for sid in strain_ids}

    return FeatureBundle(
        point_mut_vocab=point_vocab,
        acquired_gene_vocab=acquired_vocab,
        mlst_vocab=mlst_vocab,
        point_mut=point_mat,
        acquired=acquired_mat,
        mlst=mlst_mat,
        seqs=seqs,
        labels=labels,
        strain_ids=strain_ids,
    )


def run_loso(
    bundle: FeatureBundle,
    head: Literal["logreg", "xgb"],
    drug: str,
    feature_blocks: tuple[str, ...] = ("point", "acquired", "kmer", "mlst"),
    k: int = 8,
    top_n: int = 10_000,
) -> CVResult:
    """LOSO with within-fold k-mer vocab rebuild + (optional) feature-block ablation."""
    result = CVResult(strategy="loso", drug=drug)
    n = len(bundle.strain_ids)
    for i, held in enumerate(bundle.strain_ids):
        train_ids = [bundle.strain_ids[j] for j in range(n) if j != i]
        train_y = np.array([bundle.labels[s] for s in train_ids], dtype=int)
        test_y = int(bundle.labels[held])

        # Stack feature blocks
        train_blocks: list[np.ndarray] = []
        test_blocks: list[np.ndarray] = []

        if "point" in feature_blocks:
            train_blocks.append(np.stack([bundle.point_mut[s] for s in train_ids]))
            test_blocks.append(bundle.point_mut[held][None, :])
        if "acquired" in feature_blocks:
            train_blocks.append(np.stack([bundle.acquired[s] for s in train_ids]))
            test_blocks.append(bundle.acquired[held][None, :])
        if "kmer" in feature_blocks:
            train_seqs = [bundle.seqs[s] for s in train_ids]
            vocab = build_kmer_vocabulary(train_seqs, k=k, top_n=top_n)
            X_km_train = kmers_to_feature_matrix(train_seqs, vocab, k=k).astype(np.float32)
            X_km_test = kmers_to_feature_matrix([bundle.seqs[held]], vocab, k=k).astype(np.float32)
            train_blocks.append(X_km_train)
            test_blocks.append(X_km_test)
        if "mlst" in feature_blocks:
            train_blocks.append(np.stack([bundle.mlst[s] for s in train_ids]))
            test_blocks.append(bundle.mlst[held][None, :])

        X_train = np.concatenate(train_blocks, axis=1)
        X_test = np.concatenate(test_blocks, axis=1)

        if head == "logreg":
            # scaled LR — feature blocks have very different scales (binary multi-hot
            # vs k-mer count); StandardScaler equalizes
            clf = _train_baseline_logreg(
                X_train, train_y, drug_name=f"{drug}_curated_lr",
                calibrate=False, scale_features=True,
            )
            score = float(clf.model.predict_proba(X_test)[:, 1][0])
        elif head == "xgb":
            clf = train_xgboost_classifier(
                X_train, train_y, drug_name=f"{drug}_curated_xgb", calibrate=False,
            )
            score = float(predict_proba(clf, X_test)[0])
        else:
            raise ValueError(f"unknown head: {head}")

        result.folds.append(
            FoldResult(
                held_out_id=held,
                held_out_indices=[i],
                train_indices=[j for j in range(n) if j != i],
                y_true=np.array([test_y]),
                y_score=np.array([score]),
                n_train=n - 1,
                n_test=1,
            )
        )
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort", type=Path, default=Path("data/processed/gate_b_n40_cipro_cohort.parquet"))
    parser.add_argument("--mech-audit", type=Path, default=Path("wiki/cipro_mechanism_audit_2026-05-17.json"))
    parser.add_argument("--mic-audit", type=Path, default=Path("wiki/cipro_mic_audit_2026-05-17.json"))
    parser.add_argument("--refseq-cache", type=Path, default=Path("D:/dna_decode_cache/refseq"))
    parser.add_argument("--drug", default="ciprofloxacin")
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--restrict-to-decisive", action="store_true",
                        help="Restrict LOSO to HIGH_R + HIGH_S strains from MIC audit")
    # Ablation panel always runs (load-bearing for amended condition 4); no flag needed.
    args = parser.parse_args(argv)

    if args.output is None:
        args.output = Path(f"wiki/cipro_curated_baseline_{_date.today().isoformat()}.md")

    if not args.mech_audit.exists():
        print(f"[curated] FAIL mech audit JSON not found: {args.mech_audit}")
        return 2

    cohort = load_cohort(args.cohort)

    decisive: set[str] | None = None
    if args.restrict_to_decisive:
        if not args.mic_audit.exists():
            print(f"[curated] FAIL mic audit JSON not found: {args.mic_audit}")
            return 2
        mic = json.loads(args.mic_audit.read_text(encoding="utf-8"))
        decisive = set(mic["decisive_R_ids"]) | set(mic["decisive_S_ids"])
        print(f"[curated] restricting to decisive subset: N={len(decisive)}")

    bundle = build_features(cohort, args.mech_audit, args.refseq_cache, args.drug, decisive_subset=decisive)
    n_r = sum(1 for s in bundle.strain_ids if bundle.labels[s] == 1)
    n_s = sum(1 for s in bundle.strain_ids if bundle.labels[s] == 0)
    print(f"[curated] N={len(bundle.strain_ids)} ({n_r}R/{n_s}S)")
    print(f"[curated] features: POINT={len(bundle.point_mut_vocab)} acquired={len(bundle.acquired_gene_vocab)} mlst={len(bundle.mlst_vocab)} kmer=top-10000")

    if n_r < 2 or n_s < 2:
        print(f"[curated] FAIL insufficient class balance for LOSO")
        return 3

    # Run full ablation panel (named multi-block feature sets) under XGB.
    # Always run — single-block-only ablation is insufficient to gate the
    # amended PIVOT TRIGGER condition 4 (Codex round-2 critique 2026-05-17).
    print(f"[curated] running ablation panel (XGB) over {len(ABLATION_FEATURE_SETS)} feature sets...")
    ablation_results: dict[str, float] = {}
    ablation_cvs: dict[str, CVResult] = {}
    for name, blocks in ABLATION_FEATURE_SETS.items():
        print(f"  {name:18s} blocks={list(blocks)}")
        cv = run_loso(bundle, head="xgb", drug=args.drug, feature_blocks=blocks)
        m = compute_metrics(cv.all_y_true, cv.all_y_score)
        ablation_results[name] = m.auroc
        ablation_cvs[name] = cv
        print(f"  {name:18s} XGB AUROC={m.auroc:.4f}")

    # LR head over the full feature stack (sanity reference)
    print(f"[curated] running LR head (all features)...")
    cv_lr_full = run_loso(bundle, head="logreg", drug=args.drug)
    m_lr = compute_metrics(cv_lr_full.all_y_true, cv_lr_full.all_y_score)
    cv_xgb_full = ablation_cvs["all"]
    m_xgb_auroc_full = ablation_results["all"]
    print(f"[curated] LR  (all) AUROC={m_lr.auroc:.4f} AUPRC={m_lr.auprc:.4f}")
    print(f"[curated] XGB (all) AUROC={m_xgb_auroc_full:.4f}")

    # Two-layer verdict:
    #   original_condition_4: frozen rule from return_decision_tree_2026-05-16.md
    #     (>= 0.80 absolute OR beats NT by >= 10pp) using all-feature best AUROC.
    #     Structurally circular when POINT mutations dominate; reported for
    #     audit-trail discipline only.
    #   amended_condition_4: load-bearing rule per Codex round-2 critique.
    #     Requires no_POINT >= 0.773 OR mechanism_only >= 0.80. Isolates
    #     non-textbook genomic signal vs label-shaped tautology.
    best_curated_full = max(m_lr.auroc, m_xgb_auroc_full)
    best_nt = max(STAGE1B_NT_LR_AUROC, STAGE1B_NT_XGB_AUROC)
    gap_pp = (best_curated_full - best_nt) * 100
    no_point_auroc = ablation_results["no_POINT"]
    mech_only_auroc = ablation_results["mechanism_only"]

    # Original (frozen) verdict
    if best_curated_full >= BASELINE_ABSOLUTE_GATE_AUROC:
        original_verdict = "ABSOLUTE_PASS"
        original_rationale = f"best all-feature AUROC {best_curated_full:.3f} >= {BASELINE_ABSOLUTE_GATE_AUROC}"
    elif gap_pp >= BASELINE_COMPARATIVE_GATE_PP:
        original_verdict = "COMPARATIVE_PASS"
        original_rationale = f"best all-feature beats best NT by {gap_pp:.1f}pp (>= {BASELINE_COMPARATIVE_GATE_PP}pp)"
    else:
        original_verdict = "FAIL"
        original_rationale = f"best all-feature {best_curated_full:.3f} vs best NT {best_nt:.3f}; gap {gap_pp:+.1f}pp"

    # Amended (load-bearing) verdict
    if no_point_auroc >= AMENDED_NO_POINT_GATE_AUROC:
        amended_verdict = "NO_POINT_PASS"
        amended_rationale = f"no_POINT AUROC {no_point_auroc:.3f} >= {AMENDED_NO_POINT_GATE_AUROC:.3f} — non-textbook signal exists"
    elif mech_only_auroc >= AMENDED_MECHANISM_ONLY_GATE_AUROC:
        amended_verdict = "MECHANISM_ONLY_PASS"
        amended_rationale = f"mechanism_only AUROC {mech_only_auroc:.3f} >= {AMENDED_MECHANISM_ONLY_GATE_AUROC} — AMRFinder panel matches absolute target"
    else:
        amended_verdict = "FAIL"
        amended_rationale = (
            f"no_POINT {no_point_auroc:.3f} (gate {AMENDED_NO_POINT_GATE_AUROC:.3f}) and "
            f"mechanism_only {mech_only_auroc:.3f} (gate {AMENDED_MECHANISM_ONLY_GATE_AUROC}) both below load-bearing thresholds"
        )

    print(f"\nOriginal (frozen) verdict: {original_verdict}")
    print(f"  {original_rationale}")
    print(f"Amended (load-bearing) verdict: {amended_verdict}")
    print(f"  {amended_rationale}")
    print(f"  PIVOT TRIGGER condition 4 (amended, load-bearing): {'MET' if amended_verdict != 'FAIL' else 'NOT MET'}")

    # JSON sidecar
    json_path = args.output.with_suffix(".json")
    payload = {
        "drug": args.drug,
        "cohort": str(args.cohort),
        "mech_audit": str(args.mech_audit),
        "decisive_subset": bool(decisive),
        "n_total": len(bundle.strain_ids),
        "n_R": n_r,
        "n_S": n_s,
        "feature_dims": {
            "point_mutations": len(bundle.point_mut_vocab),
            "acquired_genes": len(bundle.acquired_gene_vocab),
            "mlst": len(bundle.mlst_vocab),
            "kmer_top_n": 10000,
        },
        "point_mutation_vocab": bundle.point_mut_vocab,
        "acquired_gene_vocab": bundle.acquired_gene_vocab,
        "mlst_vocab": bundle.mlst_vocab,
        "lr_auroc": m_lr.auroc, "lr_auprc": m_lr.auprc,
        "xgb_all_auroc": m_xgb_auroc_full,
        "best_curated_all_auroc": best_curated_full,
        "stage1b_baselines": {
            "nt_lr_auroc": STAGE1B_NT_LR_AUROC,
            "nt_xgb_auroc": STAGE1B_NT_XGB_AUROC,
            "kmer_xgb_auroc": STAGE1B_KMER_XGB_AUROC,
        },
        "ablations_xgb": ablation_results,
        "amended_gates": {
            "no_point_threshold": AMENDED_NO_POINT_GATE_AUROC,
            "mechanism_only_threshold": AMENDED_MECHANISM_ONLY_GATE_AUROC,
        },
        "original_condition_4_verdict": original_verdict,
        "original_condition_4_rationale": original_rationale,
        "amended_condition_4_verdict": amended_verdict,
        "amended_condition_4_rationale": amended_rationale,
        "pivot_trigger_condition_4_load_bearing": amended_verdict != "FAIL",
        "per_strain_lr": [
            {"strain_id": f.held_out_id, "y_true": int(f.y_true[0]), "y_score": float(f.y_score[0])}
            for f in cv_lr_full.folds
        ],
        "per_strain_xgb_all": [
            {"strain_id": f.held_out_id, "y_true": int(f.y_true[0]), "y_score": float(f.y_score[0])}
            for f in cv_xgb_full.folds
        ],
    }
    json_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[curated] wrote JSON: {json_path}")

    # Markdown packet
    subset_tag = " (DECISIVE SUBSET)" if decisive else ""
    lines = [
        f"# Cipro curated AMR baseline{subset_tag} — N={len(bundle.strain_ids)} ({_date.today().isoformat()})",
        "",
        f"**Purpose:** PIVOT TRIGGER condition 4. Tests whether classical AMR features (AMRFinder POINT mutations + acquired genes + k-mer + MLST) match or beat the Stage 1b NT-frozen-embedding head.",
        f"**Cohort:** `{args.cohort}` ({n_r}R/{n_s}S{', decisive subset' if decisive else ''}).",
        f"**Inputs:** AMRFinder mechanism audit JSON (`{args.mech_audit}`).",
        f"**CV:** LOSO over {len(bundle.strain_ids)} strains.",
        "",
        f"## Verdicts (2-layer)",
        "",
        f"### Original (frozen pre-Experiment-2 rule): **{original_verdict}**",
        f"- Rationale: {original_rationale}",
        "",
        f"### Amended (load-bearing post-Codex-round-2): **{amended_verdict}**",
        f"- Rationale: {amended_rationale}",
        f"- PIVOT TRIGGER condition 4 (load-bearing): **{'MET' if amended_verdict != 'FAIL' else 'NOT MET'}**",
        "",
        "Original verdict uses the all-feature AUROC and is structurally circular when POINT mutations dominate (POINT mutations like gyrA-S83L are essentially labels-in-genome-form). Amended verdict gates on no-POINT (acquired+kmer+MLST) or mechanism-only (POINT+acquired) — isolating non-textbook genomic signal vs textbook tautology.",
        "",
        "## Head comparison",
        "",
        "| head | AUROC | AUPRC | reference |",
        "|---|---:|---:|---|",
        f"| **curated LR (all)** | **{m_lr.auroc:.3f}** | {m_lr.auprc:.3f} | scaled LR over all 4 blocks |",
        f"| **curated XGB (all)** | **{m_xgb_auroc_full:.3f}** | - | XGBoost over all 4 blocks |",
        f"| Stage 1b NT-LR | {STAGE1B_NT_LR_AUROC:.3f} | - | mean+max pool + scaled LR |",
        f"| Stage 1b NT-XGB | {STAGE1B_NT_XGB_AUROC:.3f} | - | mean+max pool + XGBoost |",
        f"| Stage 1b k-mer-XGB | {STAGE1B_KMER_XGB_AUROC:.3f} | - | k=8 top-10000 |",
        "",
        f"Gap (best_curated_all - best_NT): **{gap_pp:+.1f}pp**.",
        "",
        "## Feature dimensions",
        "",
        f"- POINT mutations (cipro-relevant QUINOLONE class): **{len(bundle.point_mut_vocab)}** features",
        f"- Acquired AMR genes (all classes; includes co-resistance markers): **{len(bundle.acquired_gene_vocab)}** features",
        f"- MLST one-hot: **{len(bundle.mlst_vocab)}** features",
        f"- k-mer top-N: **10000** features (k=8, within-fold vocab rebuild)",
        "",
        f"## Ablation panel (XGB; load-bearing for amended verdict)",
        "",
        f"Each row is XGB trained on the named feature set under LOSO. Amended-verdict gates: **no_POINT >= {AMENDED_NO_POINT_GATE_AUROC:.3f}** OR **mechanism_only >= {AMENDED_MECHANISM_ONLY_GATE_AUROC}**.",
        "",
        "| feature_set | blocks | AUROC | gate? |",
        "|---|---|---:|---|",
    ]
    for name, blocks in ABLATION_FEATURE_SETS.items():
        auroc = ablation_results[name]
        gate = ""
        if name == "no_POINT":
            gate = f"PASS (>= {AMENDED_NO_POINT_GATE_AUROC:.3f})" if auroc >= AMENDED_NO_POINT_GATE_AUROC else f"FAIL (< {AMENDED_NO_POINT_GATE_AUROC:.3f})"
        elif name == "mechanism_only":
            gate = f"PASS (>= {AMENDED_MECHANISM_ONLY_GATE_AUROC})" if auroc >= AMENDED_MECHANISM_ONLY_GATE_AUROC else f"FAIL (< {AMENDED_MECHANISM_ONLY_GATE_AUROC})"
        lines.append(f"| {name} | {','.join(blocks)} | {auroc:.3f} | {gate} |")
    lines.extend([
        "",
        "## Per-strain predictions (all-feature LR)",
        "",
        "| strain_id | y_true | y_score (LR-all) | y_score (XGB-all) |",
        "|---|---:|---:|---:|",
    ])
    xgb_by_strain = {f.held_out_id: float(f.y_score[0]) for f in cv_xgb_full.folds}
    for f in cv_lr_full.folds:
        sid = f.held_out_id
        lines.append(f"| {sid} | {int(f.y_true[0])} | {float(f.y_score[0]):.3f} | {xgb_by_strain.get(sid, float('nan')):.3f} |")
    lines.extend([
        "",
        "## How to interpret",
        "",
        "- **Amended NO_POINT_PASS:** classical features WITHOUT POINT mutations reach the load-bearing threshold. Non-textbook genomic signal exists; frozen-NT pooling is falsified. Next: invest in curated feature engineering + per-gene NT windows.",
        "- **Amended MECHANISM_ONLY_PASS:** AMRFinder POINT + acquired alone reaches the absolute target. NT pooling adds nothing over the curated mechanism panel. Same falsification; the right Phase 1 product is curated baselines.",
        "- **Amended FAIL:** neither non-textbook ablation passes. Two possibilities: (a) the N=38 label noise (Experiment 2 NOISY verdict) is the bottleneck — try cohort expansion to N=150 with strict MIC filters; or (b) the right features aren't in this stack — try Bakta annotation + per-gene NT windows on the CLEAN subset.",
        "- **Original verdict** is retained for audit-trail discipline only. PASS on original with FAIL on amended = signal is circular (POINT-dominated).",
        "",
        f"_JSON sidecar: `{json_path}`_",
    ])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines), encoding="utf-8")
    print(f"[curated] wrote packet: {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
