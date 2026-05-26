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
    from dna_decode.data.cohort import find_duplicate_accessions, load_cohort
    from dna_decode.eval.clade_baseline import (
        predict_clade_only,
        train_clade_only_classifier,
        validation_gate,
    )
    from dna_decode.eval.cv import leave_one_accession_out_cv, leave_one_strain_out_cv
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

    duplicate_accessions = find_duplicate_accessions(
        cohort, restrict_to_strain_ids=drug_strain_ids
    )
    cv_grouping = args.cv_grouping
    if cv_grouping == "auto":
        cv_grouping = "assembly_accession" if duplicate_accessions else "strain_id"
    if cv_grouping == "strain_id" and duplicate_accessions:
        sample_pairs = ", ".join(
            f"{accession}: {strain_ids}"
            for accession, strain_ids in sorted(duplicate_accessions.items())[:3]
        )
        print(
            "[train] duplicate assembly_accession values detected inside this drug pool. "
            "LOSO by strain_id would leak the same assembly across train and held-out folds. "
            "Use --cv-grouping assembly_accession or auto. "
            f"Sample duplicates: {sample_pairs}",
            file=sys.stderr,
        )
        return 2

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
    accession_assignments: dict[str, str] = {}
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
        accession = (strain_obj.assembly_accession or "").strip()
        accession_assignments[sid] = accession if accession else sid
    if not X_rows:
        print(f"[train] no strains in cohort have cached embeddings", file=sys.stderr)
        return 2
    X = np.stack(X_rows)
    y = np.array(labels, dtype=int)
    print(f"[train] feature matrix: {X.shape}, label balance: {(y==1).sum()}R / {(y==0).sum()}S")
    print(f"[train] CV grouping: {cv_grouping}")

    def _train_fn(X_train, y_train):
        return train_xgboost_classifier(
            X_train, y_train, drug_name=args.drug, calibrate=True
        )

    def _predict_fn(model, X_test):
        return predict_proba(model, X_test)

    if cv_grouping == "assembly_accession":
        cv_result = leave_one_accession_out_cv(
            X,
            y,
            strain_id_order,
            accession_assignments,
            _train_fn,
            _predict_fn,
            drug=args.drug,
        )
    else:
        cv_result = leave_one_strain_out_cv(
            X, y, strain_id_order, _train_fn, _predict_fn, drug=args.drug
        )
    metrics = compute_metrics(cv_result.all_y_true, cv_result.all_y_score)
    print(
        f"[train] {cv_result.strategy} AUROC={metrics.auroc:.3f} "
        f"AUPRC={metrics.auprc:.3f} (n={metrics.n_samples})"
    )

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

    # Save trained model + provenance fields for v0 predict output
    from datetime import date as _train_date
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
                "cv_strategy": cv_result.strategy,
                "cv_grouping": cv_grouping,
                "auroc_loso": float(metrics.auroc) if not np.isnan(metrics.auroc) else None,
                "auroc_cv_primary": float(metrics.auroc) if not np.isnan(metrics.auroc) else None,
                "auroc_lomo_clade_out": None,  # v1+ — LOMO-clade-out gated on more cohorts
                "strain_id_order": strain_id_order,
                "n_strains": len(strain_id_order),
                "training_cohort": cohort_path.stem,
                "trained_on": _train_date.today().isoformat(),
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


def _confidence_tier(proba: float) -> str:
    """Map calibrated probability to HIGH / MEDIUM / LOW confidence tier.

    HIGH:   proba >= 0.9 or <= 0.1 (clear call, large margin from threshold)
    MEDIUM: proba in [0.7, 0.9) or (0.1, 0.3] (clear-ish but smaller margin)
    LOW:    proba in (0.3, 0.7) (near decision threshold)

    Compares against proba directly (not abs(proba - 0.5)) to avoid the
    floating-point precision bug at 0.7 / 0.3 boundaries.
    """
    if proba >= 0.9 or proba <= 0.1:
        return "HIGH"
    if proba >= 0.7 or proba <= 0.3:
        return "MEDIUM"
    return "LOW"


def _load_audit_verdict(
    merge_json_path: Path | None,
    strain_id: str,
    fallback_to_cohort_gate: bool = False,
) -> dict | None:
    """Read the merge-gate JSON sidecar; extract this strain's audit verdict.

    Returns None if path not provided, file missing, or strain not in the audit.
    Wires the v0 success criterion #4 ("audit-aware honest output"): if the
    training cohort had `SUSPEND_CONDITION_4`, the prediction inherits that
    framing instead of silently overclaiming.
    """
    if not merge_json_path:
        return None
    p = Path(merge_json_path)
    if not p.exists():
        return None
    import json as _json
    with open(p, encoding="utf-8") as f:
        merge = _json.load(f)
    # Cohort-level verdict + per-strain row
    cohort_gate = (
        merge.get("gate_verdict")
        or merge.get("pre_curated_gate")
        or merge.get("recommended_next_step")
    )
    for row in merge.get("per_strain", []):
        if str(row.get("strain_id")) == str(strain_id):
            verdict: dict = {
                "noise_class": row.get("noise_class"),
                "mechanism_opacity_flag": row.get("mechanism_opacity_flag"),
                "mic_tier": row.get("mic_tier"),
                "primary_mechanisms": row.get("primary_mechanisms", []),
                "co_resistance_modifiers": row.get("co_resistance_modifiers", []),
                "cohort_gate_verdict": cohort_gate,
            }
            # Surface SUSPEND warning explicitly
            if cohort_gate and "SUSPEND" in str(cohort_gate).upper():
                verdict["suspend_gate_fired"] = True
                verdict["verdict_explanation"] = (
                    f"Training cohort fired {cohort_gate}; categorical labels carry noise. "
                    "Prediction is informational only; do not deploy as clinical decision support."
                )
            else:
                verdict["suspend_gate_fired"] = False
            return verdict
    if not fallback_to_cohort_gate or cohort_gate is None:
        return None  # strain not found in audit (e.g., held-out new strain)

    verdict = {
        "noise_class": None,
        "mechanism_opacity_flag": None,
        "mic_tier": None,
        "primary_mechanisms": [],
        "co_resistance_modifiers": [],
        "cohort_gate_verdict": cohort_gate,
    }
    if "SUSPEND" in str(cohort_gate).upper():
        verdict["suspend_gate_fired"] = True
        verdict["verdict_explanation"] = (
            f"Training cohort fired {cohort_gate}; this external genome is not present in the "
            "audit cohort, so only cohort-level noise framing is available. Prediction is "
            "informational only; do not deploy as clinical decision support."
        )
    else:
        verdict["suspend_gate_fired"] = False
    return verdict


def _resolve_predict_sample_id(
    strain_id: str | None,
    sample_id: str | None,
    genome_fasta: str | None,
) -> str | None:
    """Resolve the output label for predict."""
    if strain_id:
        return strain_id
    if sample_id:
        return sample_id
    if genome_fasta:
        return Path(genome_fasta).stem
    return None


def _embed_live_gene_sequences(
    model_name: str,
    config_path: Path,
    gene_sequences: dict[str, str],
) -> dict[str, np.ndarray]:
    """Embed CDS sequences in memory using the trained foundation model."""
    from dna_decode.models.foundation import model_factory

    gene_ids = list(gene_sequences.keys())
    if not gene_ids:
        return {}
    model = model_factory(model_name, config_path=config_path)
    batch_size = 4
    out: dict[str, np.ndarray] = {}
    for i in range(0, len(gene_ids), batch_size):
        chunk_gene_ids = gene_ids[i : i + batch_size]
        chunk_embeddings = model.embed_batch(
            [gene_sequences[gene_id] for gene_id in chunk_gene_ids]
        )
        for j, gene_id in enumerate(chunk_gene_ids):
            out[gene_id] = chunk_embeddings[j].astype(np.float32)
    return out


def _extract_top_k_attribution(
    classifier,
    strain_gene_embeddings: dict,
    annotations,
    drug: str,
    catalog,
    top_k: int,
) -> list[dict]:
    """Run gene-level ISM + tier the top-K via the resistance catalog.

    Falls back to ISM-only (no tier labels) when catalog is None. Wires the
    v0 success criterion #3 ("top-K attribution includes >=1 known-mechanism gene").
    """
    from dna_decode.interp.mutagenesis import (
        build_attribution_report,
        gene_level_mutagenesis,
    )

    gene_effects = gene_level_mutagenesis(classifier, strain_gene_embeddings, annotations)
    top_rows = gene_effects.head(top_k)

    if catalog is None:
        return [
            {
                "gene_id": str(row.get("gene_id", "")),
                "locus_tag": str(row.get("locus_tag", "")),
                "score": float(row.get("prediction_delta", 0.0)),
                "tier": "Tier ? (no resistance catalog provided)",
            }
            for _, row in top_rows.iterrows()
        ]

    report = build_attribution_report(gene_effects, drug, catalog, top_k=top_k, annotations=annotations)
    tier_by_locus = {locus: tier for locus, tier in report.hits}
    out: list[dict] = []
    for _, row in top_rows.iterrows():
        locus = str(row.get("locus_tag", "") or row.get("gene_id", ""))
        tier_obj = tier_by_locus.get(locus)
        out.append({
            "gene_id": str(row.get("gene_id", "")),
            "locus_tag": str(row.get("locus_tag", "")),
            "score": float(row.get("prediction_delta", 0.0)),
            "tier": f"Tier {tier_obj.tier} ({tier_obj.rationale})" if tier_obj else "Tier 5 (no catalog match)",
        })
    return out


def _render_predict_markdown(result: dict) -> str:
    """Format the v0 prediction as a human-readable markdown report."""
    lines = [
        f"# Decode result — strain `{result['strain_id']}` / drug `{result['drug']}`",
        "",
        f"**Prediction:** {result['prediction']}",
        f"**Calibrated probability (R):** {result['calibrated_probability']:.3f}",
        f"**Confidence tier:** {result['confidence_tier']}",
        "",
        "## Top-K gene attribution",
        "",
    ]
    if not result.get("top_k_attribution"):
        lines.extend(["(no attribution run — pass `--annotations` to enable)", ""])
    else:
        lines.extend([
            "| Rank | Gene | Locus | Score | Tier |",
            "|---:|---|---|---:|---|",
        ])
        for i, hit in enumerate(result["top_k_attribution"], start=1):
            lines.append(
                f"| {i} | {hit['gene_id']} | {hit['locus_tag']} | "
                f"{hit['score']:+.4f} | {hit['tier']} |"
            )
        lines.append("")

    av = result.get("audit_verdict")
    lines.extend(["## Audit verdict", ""])
    if av is None:
        lines.extend([
            "(no audit data — pass `--audit-merge-json` to surface training-cohort noise framing)",
            "",
        ])
    else:
        if av.get("suspend_gate_fired"):
            lines.extend([
                f"**SUSPEND gate fired on training cohort.** {av.get('verdict_explanation', '')}",
                "",
            ])
        lines.extend([
            f"- Cohort gate verdict: `{av.get('cohort_gate_verdict')}`",
            f"- Strain noise class: `{av.get('noise_class')}`",
            f"- MIC tier: `{av.get('mic_tier')}`",
            f"- Mechanism opacity flag: `{av.get('mechanism_opacity_flag')}`",
            f"- Primary mechanisms: {av.get('primary_mechanisms') or 'none'}",
            f"- Co-resistance modifiers: {av.get('co_resistance_modifiers') or 'none'}",
            "",
        ])

    if av is None:
        lines[-2] = (
            "**Non-canonical internal/debug run.** No audit data was supplied, so "
            "training-cohort noise framing is missing."
        )

    prov = result.get("provenance", {})
    lines.extend([
        "## Provenance",
        "",
        f"- Model: {prov.get('model', 'unknown')}",
        f"- Training cohort: {prov.get('training_cohort', 'unknown')}",
        f"- Reporting mode: {prov.get('reporting_mode', 'unknown')}",
        f"- CV strategy: {prov.get('cv_strategy', 'loso')}",
        f"- Primary CV AUROC: {prov.get('cv_auroc', prov.get('loso_auroc'))}",
        f"- Trained on: {prov.get('trained_on', 'unknown')}",
        "",
        "---",
        "",
        "Not a clinical decision support tool. Audit verdict + provenance block must accompany any downstream interpretation.",
    ])
    return "\n".join(lines)


def cmd_predict(args: argparse.Namespace, cfg: dict) -> int:
    """Predict resistance for a cached strain — v0 schema (JSON + markdown).

    Inputs:
      --model-path        trained classifier pickle (output of `train`)
      --strain-id         BV-BRC strain ID for cached-strain v0 mode
      --genome-fasta      genome FASTA for v0.1 genome-input mode
      --sample-id         optional output label for genome-input mode
      --cache             HDF5 embedding cache path (cached-strain mode only)
      --annotations       GFF3 / GenBank file (required for genome-input mode; optional but recommended for cached attribution)
      --card-path         CARD JSON for tier labels (optional)
      --amrfinder-path    AMRFinder TSV for tier labels (optional)
      --audit-merge-json  merge-gate JSON sidecar (required for canonical reporting)
      --allow-missing-audit  permit non-canonical internal/debug output without audit
      --top-k             top-K genes to attribute (default 10)
      --output            JSON output path (markdown sidecar auto-generated alongside)
      --output-md         explicit markdown output path (defaults to <output>.md)
      --no-attribution    skip ISM attribution (faster prediction)
    """
    import json

    if bool(args.strain_id) == bool(args.genome_fasta):
        print(
            "[predict] supply exactly one input mode: --strain-id for cached prediction "
            "or --genome-fasta for v0.1 genome-input prediction.",
            file=sys.stderr,
        )
        return 2
    if args.genome_fasta and not args.annotations:
        print(
            "[predict] --annotations required with --genome-fasta "
            "(supported: GFF3 or GenBank).",
            file=sys.stderr,
        )
        return 2
    if not args.audit_merge_json and not args.allow_missing_audit:
        print(
            "[predict] --audit-merge-json required for canonical v0 reporting. "
            "Use --allow-missing-audit only for non-canonical internal/debug runs.",
            file=sys.stderr,
        )
        return 2

    model_path = Path(args.model_path)
    if not model_path.exists():
        print(f"[predict] trained model not found at {model_path}", file=sys.stderr)
        return 2
    with open(model_path, "rb") as f:
        bundle = pickle.load(f)
    drug = bundle["drug"]
    model_name = bundle["model_name"]
    print(f"[predict] loaded {drug} classifier ({model_name})")

    from dna_decode.data.annotations import extract_cds_sequences, load_annotation_table
    from dna_decode.models.classifiers import aggregate_strain_features, predict_proba

    sample_id = _resolve_predict_sample_id(args.strain_id, args.sample_id, args.genome_fasta)
    if sample_id is None:
        print("[predict] unable to resolve sample identifier", file=sys.stderr)
        return 2

    if args.strain_id:
        foundation_cfg = cfg.get("foundation_models", {}).get(model_name)
        if foundation_cfg is None:
            print(f"[predict] foundation model {model_name!r} not in config", file=sys.stderr)
            return 2
        from dna_decode.models.cache import EmbeddingCache

        cache_path = _resolve_cache_path(cfg, args.cache, model_name)
        cache = EmbeddingCache(
            cache_path,
            model_name=model_name,
            model_version=foundation_cfg["huggingface_id"],
            embedding_dim=foundation_cfg["embedding_dim"],
        )
        gene_ids = cache.list_genes(args.strain_id)
        if not gene_ids:
            print(f"[predict] no cached embeddings for strain {args.strain_id!r}", file=sys.stderr)
            return 2
        gene_matrix = cache.bulk_get([(args.strain_id, g) for g in gene_ids])
        strain_gene_embeddings = {g: cache.get(args.strain_id, g) for g in gene_ids}
        annotations = load_annotation_table(args.annotations) if args.annotations else None
        print(f"[predict] mode=cached-strain sample={sample_id} genes={len(gene_ids)}")
    else:
        annotation_table = load_annotation_table(args.annotations)
        gene_sequences = extract_cds_sequences(args.genome_fasta, annotation_table)
        if not gene_sequences:
            print(
                "[predict] no CDS sequences extracted from genome input; cannot score sample",
                file=sys.stderr,
            )
            return 2
        strain_gene_embeddings = _embed_live_gene_sequences(
            model_name=model_name,
            config_path=Path(args.config),
            gene_sequences=gene_sequences,
        )
        gene_ids = list(strain_gene_embeddings.keys())
        gene_matrix = np.stack([strain_gene_embeddings[g] for g in gene_ids])
        annotations = annotation_table
        print(f"[predict] mode=genome-input sample={sample_id} genes={len(gene_ids)}")

    X = aggregate_strain_features(gene_matrix, "mean").reshape(1, -1)
    if X.shape[1] != bundle["classifier"].feature_dim:
        print(
            f"[predict] feature_dim mismatch: live path produced {X.shape[1]}, "
            f"trained classifier expects {bundle['classifier'].feature_dim}",
            file=sys.stderr,
        )
        return 2
    proba = float(predict_proba(bundle["classifier"], X)[0])
    prediction = "R" if proba >= 0.5 else "S"
    confidence = _confidence_tier(proba)
    print(f"[predict] sample={sample_id}: {prediction} (p={proba:.3f}, conf={confidence})")

    # Optional attribution
    top_k_attr: list[dict] = []
    if args.no_attribution:
        print("[predict] attribution skipped (--no-attribution)")
    else:
        catalog = None
        if args.card_path and args.amrfinder_path:
            from dna_decode.data.resistance_db import load_amrfinder, load_card, merge_catalogs
            catalog = merge_catalogs(load_card(args.card_path), load_amrfinder(args.amrfinder_path))
        top_k_attr = _extract_top_k_attribution(
            bundle["classifier"], strain_gene_embeddings, annotations, drug, catalog, args.top_k
        )
        print(f"[predict] top-{args.top_k} attribution computed ({len(top_k_attr)} hits)")

    # Audit verdict from merge-gate JSON if provided
    audit = _load_audit_verdict(
        Path(args.audit_merge_json) if args.audit_merge_json else None,
        sample_id,
        fallback_to_cohort_gate=bool(args.audit_merge_json),
    )

    # Provenance block from the training-pickle bundle
    provenance = {
        "model": f"{model_name} + XGBoost (frozen)",
        "training_cohort": bundle.get("training_cohort", "unknown"),
        "cv_strategy": bundle.get("cv_strategy", "loso"),
        "cv_auroc": bundle.get("auroc_cv_primary", bundle.get("auroc_loso")),
        "loso_auroc": bundle.get("auroc_loso"),
        "lomo_clade_out_auroc": bundle.get("auroc_lomo_clade_out"),
        "trained_on": bundle.get("trained_on", "unknown"),
        "input_mode": "genome_input" if args.genome_fasta else "cached_strain",
        "reporting_mode": (
            "canonical_audit_aware" if audit is not None else "non_canonical_missing_audit"
        ),
    }

    result = {
        "strain_id": sample_id,
        "drug": drug,
        "prediction": prediction,
        "calibrated_probability": proba,
        "confidence_tier": confidence,
        "top_k_attribution": top_k_attr,
        "audit_verdict": audit,
        "provenance": provenance,
    }

    # Write JSON + markdown sidecar
    if args.output:
        out_json_path = Path(args.output)
        out_json_path.parent.mkdir(parents=True, exist_ok=True)
        out_json_path.write_text(json.dumps(result, indent=2), encoding="utf-8")
        print(f"[predict] JSON written: {out_json_path}")

        out_md_path = Path(args.output_md) if args.output_md else out_json_path.with_suffix(".md")
        out_md_path.write_text(_render_predict_markdown(result), encoding="utf-8")
        print(f"[predict] markdown written: {out_md_path}")
    else:
        # Console-only mode
        print(json.dumps(result, indent=2))
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
    p_train.add_argument(
        "--cv-grouping",
        choices=("auto", "strain_id", "assembly_accession"),
        default="auto",
        help=(
            "Primary CV grouping. auto=use assembly_accession when duplicates exist, "
            "otherwise strain_id."
        ),
    )
    p_train.add_argument("--include-clade-baseline", action="store_true")
    p_train.add_argument("--min-auroc", type=float, default=0.80)

    # predict (v0 — emits the schema in wiki/decoder_v0_ux_and_success_criterion.md)
    p_predict = sub.add_parser(
        "predict", help="Predict resistance from a cached strain or genome input."
    )
    p_predict.add_argument("--model-path", required=True)
    p_predict.add_argument("--strain-id", default=None)
    p_predict.add_argument("--genome-fasta", default=None, help="Genome FASTA for v0.1 genome-input mode.")
    p_predict.add_argument("--sample-id", default=None, help="Optional output label for genome-input mode.")
    p_predict.add_argument("--cache", default=None)
    p_predict.add_argument(
        "--annotations",
        default=None,
        help="GFF3 or GenBank file. Required with --genome-fasta; enables attribution in cached mode.",
    )
    p_predict.add_argument("--card-path", default=None, help="CARD JSON for tier labels.")
    p_predict.add_argument("--amrfinder-path", default=None, help="AMRFinder TSV for tier labels.")
    p_predict.add_argument(
        "--audit-merge-json",
        default=None,
        help=(
            "Merge-gate JSON sidecar from `cipro_mechanism_phenotype_merge.py` "
            "(required for canonical reporting; propagates SUSPEND verdict + per-strain noise class)."
        ),
    )
    p_predict.add_argument(
        "--allow-missing-audit",
        action="store_true",
        help=(
            "Permit non-canonical internal/debug prediction without --audit-merge-json. "
            "Result will carry audit_verdict=null."
        ),
    )
    p_predict.add_argument("--top-k", type=int, default=10)
    p_predict.add_argument("--output", default=None, help="JSON output path. Markdown sidecar auto-generated as <output>.md.")
    p_predict.add_argument("--output-md", default=None, help="Explicit markdown output path (overrides default).")
    p_predict.add_argument(
        "--no-attribution", action="store_true", help="Skip ISM attribution (faster).",
    )

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
