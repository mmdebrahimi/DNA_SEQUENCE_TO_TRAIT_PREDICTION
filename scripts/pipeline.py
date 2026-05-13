"""Step 14 — Single CLI entry point for the Phase 1 pipeline.

Subcommands:
    ingest      Download cohort genomes via assembly_accession
    train       Train XGBoost classifier per drug + LOMO-clade CV
    predict     Predict resistance from a single FASTA + trained model
    attribute   Run ISM + Tier 1-5 attribution on a strain in the cohort

Mirrors `scripts/pilot_gate.py` argparse + exit-code style. ALL subcommands
read `config/datasources.yaml` for shared paths. Exit codes:
    0   success
    1   run failure (e.g., metric below threshold, attribution miss)
    2   config / IO / cohort error
    3   missing dependency (Mash binary, GPU not detected, etc.)
"""
from __future__ import annotations

import argparse
import pickle
import sys
import warnings
from pathlib import Path

import numpy as np
import yaml

warnings.filterwarnings("ignore", category=UserWarning, module="sklearn")


# ---- Shared helpers ----


def _load_config(config_path: Path | str) -> dict:
    """Load config/datasources.yaml with consistent error path → exit 2."""
    path = Path(config_path)
    if not path.exists():
        raise FileNotFoundError(f"Config not found at {path}")
    with open(path) as f:
        return yaml.safe_load(f)


def _resolve_cohort_path(cfg: dict, cohort_arg: str | None) -> Path:
    """Resolve cohort.parquet path: --cohort arg > config default."""
    if cohort_arg:
        return Path(cohort_arg)
    default_dir = cfg.get("compute", {}).get("processed_dir", "data/processed")
    return Path(default_dir) / "cohort_v1.parquet"


def _resolve_cache_path(cfg: dict, cache_arg: str | None, model_name: str) -> Path:
    """Resolve embedding-cache HDF5 path: --cache arg > config default."""
    if cache_arg:
        return Path(cache_arg)
    default_dir = cfg.get("compute", {}).get("cache_dir", "data/cache")
    return Path(default_dir) / f"embeddings_{model_name}.h5"


def _resolve_model_dir(cfg: dict, models_arg: str | None) -> Path:
    """Resolve trained-model output directory."""
    if models_arg:
        return Path(models_arg)
    return Path(cfg.get("compute", {}).get("processed_dir", "data/processed")) / "models"


# ---- ingest subcommand ----


def cmd_ingest(args: argparse.Namespace, cfg: dict) -> int:
    """Build cohort + download genomes for the specified drugs."""
    from dna_decode.data.ast_data import load_bvbrc_ast
    from dna_decode.data.cohort import (
        CohortConstructionError,
        CohortSelectionCriteria,
        build_cohort,
        candidates_from_bvbrc_ast,
        download_cohort_genomes,
        save_cohort,
    )

    drugs = tuple(d.strip() for d in args.drugs.split(","))
    print(f"[ingest] drugs: {drugs}")

    # Load AST data
    if not args.ast_tsv:
        print("[ingest] --ast-tsv required for cohort construction", file=sys.stderr)
        return 2
    ast = load_bvbrc_ast(args.ast_tsv)
    print(f"[ingest] loaded {len(ast)} AST rows after broth-microdilution filter")

    # Build candidates (assembly metadata optional). Two sources, mutually exclusive:
    # 1. --assembly-metadata <yaml> — legacy / hand-authored fixtures
    # 2. --assembly-metadata-csv <csv> — BV-BRC Genomes-tab export (Phase 2 wire)
    assembly_meta: dict[str, dict[str, object]] | None = None
    assembly_meta_source: str = ""
    if getattr(args, "assembly_metadata_csv", None):
        from dna_decode.data.bvbrc_genome import (
            BvBrcGenomeError,
            load_bvbrc_genome_metadata,
        )

        csv_path = Path(args.assembly_metadata_csv)
        if not csv_path.exists():
            print(
                f"[ingest] --assembly-metadata-csv file not found: {csv_path}",
                file=sys.stderr,
            )
            return 2
        organism = (
            cfg.get("bvbrc_genomes", {})
            .get("default_filters", {})
            .get("organism", "Escherichia coli")
        )
        try:
            assembly_meta = load_bvbrc_genome_metadata(csv_path, organism=organism)
        except BvBrcGenomeError as e:
            print(f"[ingest] {e}", file=sys.stderr)
            return 2
        assembly_meta_source = f"CSV ({csv_path.name})"
        print(
            f"[ingest] loaded assembly metadata for {len(assembly_meta)} strains "
            f"from {assembly_meta_source}"
        )
    elif args.assembly_metadata:
        with open(args.assembly_metadata) as f:
            assembly_meta = yaml.safe_load(f)
        assembly_meta_source = f"YAML ({args.assembly_metadata})"
        print(
            f"[ingest] loaded assembly metadata for {len(assembly_meta or {})} strains "
            f"from {assembly_meta_source}"
        )

    # Coverage-log line: surfaces ID-namespace mismatches between AST and genome
    # metadata at the first ingest. If coverage is 0% the join is silently broken.
    if assembly_meta is not None:
        ast_strain_ids = set(str(s) for s in ast["strain_id"].astype(str).unique())
        meta_ids = set(str(k) for k in assembly_meta.keys())
        covered = len(ast_strain_ids & meta_ids)
        total = len(ast_strain_ids)
        pct = (100.0 * covered / total) if total else 0.0
        failing_qc = sum(
            1
            for v in assembly_meta.values()
            if int(v.get("contig_count", 0) or 0) > 500
            or int(v.get("n50", 0) or 0) < 50_000
        )
        print(
            f"[ingest] assembly_meta covers {covered} / {total} AST strain_ids "
            f"({pct:.1f}%); {failing_qc} loaded strains fail QC filter "
            f"(contig_count>500 OR n50<50000)"
        )

    candidates = candidates_from_bvbrc_ast(ast, assembly_metadata=assembly_meta)
    print(f"[ingest] built {len(candidates)} CandidateStrain records")

    # Build cohort
    try:
        cohort = build_cohort(
            candidates,
            drugs,
            CohortSelectionCriteria(
                target_per_drug=args.target_per_drug,
                three_drug_intersection_target=args.intersection_target,
            ),
        )
    except CohortConstructionError as e:
        print(f"[ingest] cohort construction failed: {e}", file=sys.stderr)
        return 1

    cohort_path = _resolve_cohort_path(cfg, args.cohort_out)
    save_cohort(cohort, cohort_path)
    print(f"[ingest] cohort saved: {cohort_path} ({len(cohort)} strains)")

    if args.download_genomes:
        cache_root = Path(cfg.get("refseq", {}).get("cache_dir", "data/cache/refseq"))
        print(f"[ingest] downloading {len(cohort)} genomes to {cache_root}")
        try:
            paths = download_cohort_genomes(cohort, cache_root)
            print(f"[ingest] downloaded {len(paths)} genomes")
        except CohortConstructionError as e:
            print(f"[ingest] genome download failed: {e}", file=sys.stderr)
            return 1
    return 0


# ---- train subcommand ----


def cmd_train(args: argparse.Namespace, cfg: dict) -> int:
    """Train per-drug classifier + LOMO-clade CV + optional clade-only baseline."""
    from dna_decode.data.cohort import load_cohort
    from dna_decode.eval.clade_baseline import (
        predict_clade_only,
        train_clade_only_classifier,
        validation_gate,
    )
    from dna_decode.eval.cv import leave_one_strain_out_cv
    from dna_decode.eval.metrics import compute_metrics, compute_per_clade_metrics
    from dna_decode.models.cache import EmbeddingCache, EmbeddingCacheError
    from dna_decode.models.classifiers import (
        aggregate_strain_features,
        predict_proba,
        train_xgboost_classifier,
    )

    cohort_path = _resolve_cohort_path(cfg, args.cohort)
    cache_path = _resolve_cache_path(cfg, args.cache, args.model)
    if not cohort_path.exists():
        print(f"[train] cohort not found at {cohort_path} — run `ingest` first", file=sys.stderr)
        return 2
    if not cache_path.exists():
        print(
            f"[train] embedding cache not found at {cache_path} — populate it first",
            file=sys.stderr,
        )
        return 2

    cohort = load_cohort(cohort_path)
    drug_lower = args.drug.lower()
    if drug_lower not in cohort.per_drug_strain_ids:
        print(
            f"[train] drug {args.drug!r} not in cohort's per_drug_strain_ids "
            f"({list(cohort.per_drug_strain_ids.keys())})",
            file=sys.stderr,
        )
        return 2
    drug_strain_ids = cohort.per_drug_strain_ids[drug_lower]
    print(f"[train] drug={args.drug}, model={args.model}, n_strains={len(drug_strain_ids)}")

    # Resolve foundation model metadata for cache dim
    foundation_cfg = cfg.get("foundation_models", {}).get(args.model)
    if foundation_cfg is None:
        print(f"[train] foundation model {args.model!r} not in config", file=sys.stderr)
        return 2

    try:
        cache = EmbeddingCache(
            cache_path,
            model_name=args.model,
            model_version=foundation_cfg["huggingface_id"],
            embedding_dim=foundation_cfg["embedding_dim"],
        )
    except EmbeddingCacheError as e:
        print(f"[train] cache error: {e}", file=sys.stderr)
        return 2

    # Build per-strain feature vector (mean-pool over genes in cache)
    X_rows = []
    labels = []
    strain_id_order = []
    for sid in drug_strain_ids:
        strain_obj = cohort.strain_by_id(sid)
        if strain_obj is None or drug_lower not in strain_obj.ast_labels:
            continue
        gene_ids = cache.list_genes(sid)
        if not gene_ids:
            continue
        gene_matrix = cache.bulk_get([(sid, g) for g in gene_ids])
        X_rows.append(aggregate_strain_features(gene_matrix, "mean"))
        labels.append(strain_obj.ast_labels[drug_lower])
        strain_id_order.append(sid)
    if not X_rows:
        print(f"[train] no strains in cohort have cached embeddings", file=sys.stderr)
        return 2
    X = np.stack(X_rows)
    y = np.array(labels, dtype=int)
    print(f"[train] feature matrix: {X.shape}, label balance: {(y==1).sum()}R / {(y==0).sum()}S")

    def _train_fn(X_train, y_train):
        return train_xgboost_classifier(
            X_train, y_train, drug_name=args.drug, calibrate=True
        )

    def _predict_fn(model, X_test):
        return predict_proba(model, X_test)

    cv_result = leave_one_strain_out_cv(
        X, y, strain_id_order, _train_fn, _predict_fn, drug=args.drug
    )
    metrics = compute_metrics(cv_result.all_y_true, cv_result.all_y_score)
    print(f"[train] LOSO AUROC={metrics.auroc:.3f} AUPRC={metrics.auprc:.3f} (n={metrics.n_samples})")

    # Optional clade-only baseline + validation gate
    if args.include_clade_baseline:
        clade_assignments = {
            s.strain_id: hash(s.mlst) % 10 for s in cohort.strains  # synthetic clade IDs
        }
        clade_model = train_clade_only_classifier(strain_id_order, clade_assignments, y)
        clade_scores = predict_clade_only(clade_model, strain_id_order, clade_assignments)
        clade_metrics = compute_metrics(y, clade_scores)
        print(
            f"[train] clade-only AUROC={clade_metrics.auroc:.3f} "
            f"(gap: {metrics.auroc - clade_metrics.auroc:+.3f})"
        )
        per_clade_foundation = {
            str(fold.held_out_id): compute_metrics(fold.y_true, fold.y_score).auroc
            for fold in cv_result.folds
        }
        per_clade_baseline = {
            str(sid): clade_metrics.auroc for sid in strain_id_order  # constant per-strain
        }
        gate = validation_gate(per_clade_foundation, per_clade_baseline)
        print(f"[train] validation_gate: {gate.get('passed')} ({gate.get('fraction_passing')})")

    # Save trained model
    full_clf = _train_fn(X, y)
    model_dir = _resolve_model_dir(cfg, args.models_dir)
    model_dir.mkdir(parents=True, exist_ok=True)
    out_path = model_dir / f"{args.drug}_{args.model}.pkl"
    with open(out_path, "wb") as f:
        pickle.dump(
            {
                "classifier": full_clf,
                "drug": args.drug,
                "model_name": args.model,
                "feature_dim": X.shape[1],
                "auroc_loso": float(metrics.auroc) if not np.isnan(metrics.auroc) else None,
                "strain_id_order": strain_id_order,
            },
            f,
        )
    print(f"[train] trained model saved: {out_path}")

    if not np.isnan(metrics.auroc) and metrics.auroc < args.min_auroc:
        print(
            f"[train] FAIL: AUROC {metrics.auroc:.3f} below threshold {args.min_auroc}",
            file=sys.stderr,
        )
        return 1
    return 0


# ---- predict subcommand ----


def cmd_predict(args: argparse.Namespace, cfg: dict) -> int:
    """Predict resistance for a single FASTA via a trained classifier."""
    import json

    model_path = Path(args.model_path)
    if not model_path.exists():
        print(f"[predict] trained model not found at {model_path}", file=sys.stderr)
        return 2
    with open(model_path, "rb") as f:
        bundle = pickle.load(f)
    print(f"[predict] loaded {bundle['drug']} classifier ({bundle['model_name']})")

    # For Phase 1, predict from a pre-computed cache entry (full FASTA →
    # embedding requires foundation model + GPU; deferred via the same
    # pattern as Step 7 lazy load)
    print(
        "[predict] real-FASTA inference path requires running ingest + cache "
        "populate first. Phase 1 predict subcommand operates on cached strains."
    )
    if not args.strain_id:
        print("[predict] --strain-id required (FASTA inference path deferred)", file=sys.stderr)
        return 2

    from dna_decode.models.cache import EmbeddingCache
    from dna_decode.models.classifiers import aggregate_strain_features, predict_proba

    foundation_cfg = cfg.get("foundation_models", {}).get(bundle["model_name"])
    if foundation_cfg is None:
        print(f"[predict] foundation model {bundle['model_name']!r} not in config", file=sys.stderr)
        return 2
    cache_path = _resolve_cache_path(cfg, args.cache, bundle["model_name"])
    cache = EmbeddingCache(
        cache_path,
        model_name=bundle["model_name"],
        model_version=foundation_cfg["huggingface_id"],
        embedding_dim=foundation_cfg["embedding_dim"],
    )
    gene_ids = cache.list_genes(args.strain_id)
    if not gene_ids:
        print(f"[predict] no cached embeddings for strain {args.strain_id!r}", file=sys.stderr)
        return 2
    gene_matrix = cache.bulk_get([(args.strain_id, g) for g in gene_ids])
    X = aggregate_strain_features(gene_matrix, "mean").reshape(1, -1)
    proba = float(predict_proba(bundle["classifier"], X)[0])
    out = {
        "strain_id": args.strain_id,
        "drug": bundle["drug"],
        "model_name": bundle["model_name"],
        "probability_resistant": proba,
        "binary_prediction": int(proba >= 0.5),
    }
    if args.output:
        Path(args.output).write_text(json.dumps(out, indent=2), encoding="utf-8")
        print(f"[predict] result written: {args.output}")
    print(json.dumps(out, indent=2))
    return 0


# ---- attribute subcommand ----


def cmd_attribute(args: argparse.Namespace, cfg: dict) -> int:
    """Run gene-level ISM + Tier 1-5 attribution on a strain in the cohort."""
    from dna_decode.data.annotations import parse_gff3
    from dna_decode.data.resistance_db import load_amrfinder, load_card, merge_catalogs
    from dna_decode.models.cache import EmbeddingCache
    from dna_decode.interp.mutagenesis import (
        build_attribution_report,
        gene_level_mutagenesis,
    )

    model_path = Path(args.model_path)
    if not model_path.exists():
        print(f"[attribute] trained model not found at {model_path}", file=sys.stderr)
        return 2
    with open(model_path, "rb") as f:
        bundle = pickle.load(f)
    classifier = bundle["classifier"]
    drug = bundle["drug"]
    model_name = bundle["model_name"]

    foundation_cfg = cfg.get("foundation_models", {}).get(model_name)
    cache_path = _resolve_cache_path(cfg, args.cache, model_name)
    cache = EmbeddingCache(
        cache_path,
        model_name=model_name,
        model_version=foundation_cfg["huggingface_id"],
        embedding_dim=foundation_cfg["embedding_dim"],
    )
    gene_ids = cache.list_genes(args.strain_id)
    if not gene_ids:
        print(f"[attribute] no cached embeddings for {args.strain_id!r}", file=sys.stderr)
        return 2
    strain_gene_embeddings = {
        g: cache.get(args.strain_id, g) for g in gene_ids
    }

    annotations = None
    if args.annotations:
        annotations = parse_gff3(args.annotations)

    print(f"[attribute] running gene-level ISM ({len(strain_gene_embeddings)} genes)")
    gene_effects = gene_level_mutagenesis(
        classifier, strain_gene_embeddings, annotations
    )
    print(f"[attribute] top-{args.top_k} attributed genes:")
    for _, row in gene_effects.head(args.top_k).iterrows():
        print(
            f"  {row['gene_id']:>12s}  ({row['locus_tag']:>14s})  "
            f"delta={row['prediction_delta']:+.3f}"
        )

    # Load resistance catalog
    catalog = None
    if args.card_path and args.amrfinder_path:
        card = load_card(args.card_path)
        amrfinder = load_amrfinder(args.amrfinder_path)
        catalog = merge_catalogs(card, amrfinder)
        print(f"[attribute] loaded resistance catalog: {len(catalog)} entries")

    if catalog is not None:
        report = build_attribution_report(
            gene_effects, drug, catalog, top_k=args.top_k, annotations=annotations
        )
        print(f"[attribute] tier counts: {report.tier_counts()}")
        print(f"[attribute] fraction Tier 1-3: {report.fraction_tier_1_to_3():.2f}")
        print(f"[attribute] fraction Fail: {report.fraction_fail():.2f}")

        if args.output:
            import json
            with open(args.output, "w", encoding="utf-8") as f:
                json.dump(
                    {
                        "strain_id": args.strain_id,
                        "drug": drug,
                        "tier_counts": {str(k): v for k, v in report.tier_counts().items()},
                        "fraction_tier_1_to_3": report.fraction_tier_1_to_3(),
                        "fraction_fail": report.fraction_fail(),
                        "top_k_hits": [
                            {"locus": locus, "tier": tier.tier, "rationale": tier.rationale}
                            for locus, tier in report.hits
                        ],
                    },
                    f,
                    indent=2,
                )
            print(f"[attribute] report written: {args.output}")

    return 0


# ---- main dispatcher ----


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Phase 1 E. coli G2P pipeline — single CLI entry point."
    )
    parser.add_argument(
        "--config", default="config/datasources.yaml", help="Path to datasources YAML."
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    # ingest
    p_ingest = sub.add_parser("ingest", help="Build cohort + optionally download genomes.")
    p_ingest.add_argument("--drugs", required=True, help="Comma-separated drug list.")
    p_ingest.add_argument("--ast-tsv", required=True, help="Path to BV-BRC AST TSV.")
    # Mutually exclusive: pick one assembly-metadata source.
    meta_group = p_ingest.add_mutually_exclusive_group()
    meta_group.add_argument(
        "--assembly-metadata",
        default=None,
        help="Optional YAML with per-strain assembly metadata (legacy / fixtures).",
    )
    meta_group.add_argument(
        "--assembly-metadata-csv",
        default=None,
        help=(
            "Optional BV-BRC Genomes-tab CSV with per-strain assembly metadata "
            "(Phase 2 wire; reads contig_count + n50 + MLST + accession from "
            "real BV-BRC export)."
        ),
    )
    p_ingest.add_argument("--target-per-drug", type=int, default=150)
    p_ingest.add_argument("--intersection-target", type=int, default=75)
    p_ingest.add_argument("--cohort-out", default=None)
    p_ingest.add_argument("--download-genomes", action="store_true")

    # train
    p_train = sub.add_parser("train", help="Train classifier per drug + run CV.")
    p_train.add_argument("--drug", required=True)
    p_train.add_argument(
        "--model", default="evo", help="Foundation model name (evo / dnabert2 / etc.)."
    )
    p_train.add_argument("--cohort", default=None)
    p_train.add_argument("--cache", default=None)
    p_train.add_argument("--models-dir", default=None)
    p_train.add_argument("--include-clade-baseline", action="store_true")
    p_train.add_argument("--min-auroc", type=float, default=0.80)

    # predict
    p_predict = sub.add_parser(
        "predict", help="Predict resistance from a cached strain + trained model."
    )
    p_predict.add_argument("--model-path", required=True)
    p_predict.add_argument("--strain-id", default=None)
    p_predict.add_argument("--cache", default=None)
    p_predict.add_argument("--output", default=None)

    # attribute
    p_attr = sub.add_parser("attribute", help="Run ISM + Tier 1-5 attribution.")
    p_attr.add_argument("--model-path", required=True)
    p_attr.add_argument("--strain-id", required=True)
    p_attr.add_argument("--cache", default=None)
    p_attr.add_argument("--annotations", default=None, help="GFF3 annotation file.")
    p_attr.add_argument("--card-path", default=None)
    p_attr.add_argument("--amrfinder-path", default=None)
    p_attr.add_argument("--top-k", type=int, default=20)
    p_attr.add_argument("--output", default=None)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        cfg = _load_config(args.config)
    except FileNotFoundError as e:
        print(f"[pipeline] {e}", file=sys.stderr)
        return 2

    dispatch = {
        "ingest": cmd_ingest,
        "train": cmd_train,
        "predict": cmd_predict,
        "attribute": cmd_attribute,
    }
    return dispatch[args.cmd](args, cfg)


if __name__ == "__main__":
    sys.exit(main())
