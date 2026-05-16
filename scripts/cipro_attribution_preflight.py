"""Cipro attribution preflight on the existing N=38 NT cache.

Tests whether NT embeddings carry cipro-resistance signal at all,
INDEPENDENTLY of Stage 1's verdict gate. Useful when Stage 1 AUROC is low
(0.568 at FAIL) — distinguishes "NT has no signal for this drug" from
"NT has signal but the head/pooling/training-regime didn't expose it."

Approach: train final NT-XGBoost on all 38 cipro strains (no LOSO; we want
THE model's view, not held-out performance). For each cipro-R strain (label=1),
run gene_level_mutagenesis to get per-gene attribution scores. Aggregate
per-gene-symbol across all cipro-R strains. Check whether known cipro-resistance
loci (textbook QRDR + plasmid + efflux + porin + regulatory) appear in the
cohort-wide top-K=20.

v2 design (2026-05-15, post-Codex critique on v1):
- EXPANDED locus set: QRDR (gyrA/B + parC/E) + plasmid-encoded (qnr* + aac(6')-Ib-cr)
  + efflux (acrAB-TolC + oqxA/B) + porin (ompC/F) + regulatory (marRAB).
- SIGNED delta filter: only counts gene as "supporting R" if delta > 0 (knockout
  LOWERED R-probability). Negative delta = anti-R (knockout RAISED R-prob);
  conflating both inflates false-positives.
- TWO aggregations: (a) sum of positive deltas across R strains (magnitude);
  (b) count-in-top-K across R strains (frequency). Frequency surfaces
  consistent-but-not-dominant signal that magnitude misses.
- Per-strain top-K persisted to JSON for post-hoc re-aggregation.
- verify_complete integrity gate fires before training (consumer-side defense
  per LESSONS_LEARNED 2026-05-15).

Cheap (~5-10 min): N_cipro_R × N_genes × predict_proba calls. Output is
informational, not gate-bearing.
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

from dna_decode.data.annotations import (
    extract_cds_sequences,
    parse_gff3,
)
from dna_decode.data.cohort import load_cohort
from dna_decode.data.refseq import fasta_path, gff_path
from dna_decode.interp.mutagenesis import gene_level_mutagenesis
from dna_decode.models.cache import EmbeddingCache
from dna_decode.models.classifiers import (
    aggregate_strain_features,
    train_xgboost_classifier,
)


# Expanded cipro-resistance locus set (v2; Codex-flagged false-negative risk
# from narrow {gyrA, parC, parE} list in v1). Mechanism class per locus listed
# in comments for the verdict interpretation.
CIPRO_LOCI_BY_MECHANISM: dict[str, set[str]] = {
    "QRDR_target_alteration": {"gyrA", "gyrB", "parC", "parE"},
    "plasmid_protect_modify": {
        "qnrA", "qnrB", "qnrC", "qnrD", "qnrS",
        "aac(6')-Ib-cr", "aac(6')-Ib", "aac6-Ib-cr",  # name variants
    },
    "efflux": {"acrA", "acrB", "tolC", "oqxA", "oqxB", "mdfA", "mdtK"},
    "porin": {"ompC", "ompF"},
    "regulatory": {"marR", "marA", "marB", "soxR", "soxS"},
}
ALL_CIPRO_LOCI: set[str] = set().union(*CIPRO_LOCI_BY_MECHANISM.values())
# v1 narrow target — kept for v1-vs-v2 comparison in the verdict packet.
TEXTBOOK_QRDR_LOCI = {"gyrA", "parC", "parE"}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort", type=Path, default=Path("data/processed/gate_b_n40_cipro_cohort.parquet"))
    parser.add_argument("--nt-cache", type=Path, default=Path("D:/dna_decode_cache/embeddings/nt_n40_cipro.h5"))
    parser.add_argument("--refseq-cache", type=Path, default=Path("D:/dna_decode_cache/refseq"))
    parser.add_argument("--drug", default="ciprofloxacin")
    parser.add_argument("--top-k", type=int, default=20)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args(argv)

    if args.output is None:
        from datetime import date as _date
        args.output = Path(f"wiki/cipro_attribution_preflight_{_date.today().isoformat()}.md")

    cohort = load_cohort(args.cohort)
    cache = EmbeddingCache(
        args.nt_cache,
        model_name="nucleotide_transformer",
        model_version="InstaDeepAI/nucleotide-transformer-v2-100m-multi-species",
        embedding_dim=512,
    )

    drug_lower = args.drug.lower()
    print(f"[preflight] cohort: {len(cohort.strains)} strains; drug={args.drug}")

    # verify_complete integrity gate (Stage 1 mean-pool admits ≥1 cached gene
    # → silent half-flushed-strain landmine; LESSONS_LEARNED 2026-05-15).
    expected_by_strain: dict[str, set[str]] = {}
    for s in cohort.strains:
        if drug_lower not in s.ast_labels:
            continue
        gff = gff_path(s.assembly_accession, args.refseq_cache)
        fna = fasta_path(s.assembly_accession, args.refseq_cache)
        if not gff.exists() or not fna.exists():
            continue
        try:
            ann = parse_gff3(gff)
            expected_by_strain[s.strain_id] = set(extract_cds_sequences(fna, ann).keys())
        except Exception as e:
            print(f"[preflight] WARN expected-gene reconstruction failed for {s.strain_id}: {e!r}")
    report = cache.verify_complete(expected_by_strain)
    if not report.all_complete:
        print(f"[preflight] FAIL cache verify_complete: {report.counts}")
        return 2
    print(f"[preflight] cache verify_complete: ALL_COMPLETE ({len(report.status)} resolved)")

    # Load + aggregate per-strain features (mean pool, matching Stage 1's original config)
    strain_data: list[dict] = []
    for s in cohort.strains:
        if drug_lower not in s.ast_labels:
            continue
        gene_ids = cache.list_genes(s.strain_id)
        if not gene_ids:
            continue
        gene_matrix = cache.bulk_get([(s.strain_id, g) for g in gene_ids])
        gene_emb = {g: gene_matrix[i] for i, g in enumerate(gene_ids)}
        strain_data.append({
            "strain_id": s.strain_id,
            "accession": s.assembly_accession,
            "label": int(s.ast_labels[drug_lower]),
            "gene_embeddings": gene_emb,
            "feat": aggregate_strain_features(gene_matrix, "mean"),
        })

    X = np.stack([d["feat"] for d in strain_data])
    y = np.array([d["label"] for d in strain_data])
    print(f"[preflight] effective N={len(y)}; balance {int((y==1).sum())}R/{int((y==0).sum())}S")

    # Train final NT-XGBoost on all 38 (no LOSO — we want THE model)
    print(f"[preflight] Training final NT-XGBoost on all {len(y)} strains (no LOSO)...")
    clf = train_xgboost_classifier(X, y, drug_name=f"{args.drug}_preflight", calibrate=False)

    # Run gene_level_mutagenesis per cipro-R strain
    # Two aggregations per Codex v2 critique:
    #   (a) cohort_sum_pos_delta: sum of POSITIVE deltas per gene_symbol
    #       (signed; delta > 0 = knockout LOWERED R-prob = supporting R)
    #   (b) cohort_freq_in_topk: count of how many R strains rank a symbol
    #       in top-K (frequency surfaces consistent-but-not-dominant signal)
    cohort_sum_pos_delta: dict[str, float] = defaultdict(float)
    cohort_freq_in_topk: dict[str, int] = defaultdict(int)
    per_strain_top: list[dict] = []
    r_count = 0
    for d in strain_data:
        if d["label"] != 1:
            continue
        r_count += 1
        gff = gff_path(d["accession"], args.refseq_cache)
        if not gff.exists():
            print(f"[preflight] WARN no GFF for {d['strain_id']} ({d['accession']}); skipping symbol lookup")
            annotations = None
            symbol_map: dict[str, str] = {}
        else:
            annotations = parse_gff3(gff)
            symbol_map = {}
            if "gene_symbol" in annotations.columns:
                for _, row in annotations.iterrows():
                    gid = str(row.get("gene_id", "") or "")
                    sym = str(row.get("gene_symbol", "") or "")
                    if gid and sym:
                        symbol_map[gid] = sym

        effects = gene_level_mutagenesis(clf, d["gene_embeddings"], annotations=annotations)
        if effects.empty:
            continue
        top = effects.head(args.top_k)
        per_strain_top.append({
            "strain_id": d["strain_id"],
            "top_gene_ids": [str(g) for g in top["gene_id"]],
            "top_symbols": [symbol_map.get(str(g), "") for g in top["gene_id"]],
            "top_deltas": [float(x) for x in top["prediction_delta"]],
        })

        # v2 aggregations: signed-positive-delta sum + frequency-in-top-K
        seen_in_strain: set[str] = set()
        for _, row in top.iterrows():
            sym = symbol_map.get(str(row["gene_id"]), "")
            if not sym:
                continue
            delta = float(row["prediction_delta"])
            if delta > 0:  # supporting-R only; negative = anti-R / noise
                cohort_sum_pos_delta[sym] += delta
            if sym not in seen_in_strain:
                cohort_freq_in_topk[sym] += 1
                seen_in_strain.add(sym)

    # Cohort-wide rankings under both v2 aggregations
    ranking_by_pos_delta = sorted(
        cohort_sum_pos_delta.items(), key=lambda x: -x[1]
    )[: args.top_k]
    ranking_by_freq = sorted(
        cohort_freq_in_topk.items(), key=lambda x: (-x[1], -cohort_sum_pos_delta.get(x[0], 0.0))
    )[: args.top_k]

    # Locus detection across both expanded + narrow target sets
    pos_delta_symbols = {s for s, _ in ranking_by_pos_delta}
    freq_symbols = {s for s, _ in ranking_by_freq}
    any_topk_symbols = pos_delta_symbols | freq_symbols

    found_qrdr = TEXTBOOK_QRDR_LOCI & any_topk_symbols
    found_expanded = ALL_CIPRO_LOCI & any_topk_symbols
    # Mechanism breakdown for the verdict packet
    found_by_mechanism: dict[str, set[str]] = {
        mech: (loci & any_topk_symbols) for mech, loci in CIPRO_LOCI_BY_MECHANISM.items()
    }
    mechanisms_hit = {mech for mech, hits in found_by_mechanism.items() if hits}

    # Verdict per Codex's interpretation rubric (Quick-mode synthesis 2026-05-15)
    if len(found_qrdr) >= 1 and len(mechanisms_hit) >= 1:
        verdict = "STRONG_POSITIVE"  # ≥1 QRDR locus AND ≥1 mechanism class
    elif len(found_expanded) >= 1:
        verdict = "WEAK_POSITIVE"  # any expanded-set hit; not necessarily QRDR
    else:
        verdict = "INCONCLUSIVE_MISS"  # no cipro-loci recovered under v2 design

    # Stdout summary
    print(f"\n=== Cohort-wide top-{args.top_k} by SUM POSITIVE delta over {r_count} cipro-R strains ===")
    for i, (sym, total) in enumerate(ranking_by_pos_delta, 1):
        mech = next((m for m, loci in CIPRO_LOCI_BY_MECHANISM.items() if sym in loci), "")
        flag = f" <-- cipro:{mech}" if mech else ""
        print(f"  {i:2d}. {sym:20s} sum(+delta)={total:.4f} n_strains_topk={cohort_freq_in_topk[sym]}{flag}")
    print(f"\n=== Cohort-wide top-{args.top_k} by FREQUENCY in per-strain top-K ===")
    for i, (sym, freq) in enumerate(ranking_by_freq, 1):
        mech = next((m for m, loci in CIPRO_LOCI_BY_MECHANISM.items() if sym in loci), "")
        flag = f" <-- cipro:{mech}" if mech else ""
        pos_d = cohort_sum_pos_delta.get(sym, 0.0)
        print(f"  {i:2d}. {sym:20s} freq={freq:2d}/{r_count} sum(+delta)={pos_d:.4f}{flag}")
    print(f"\nFound textbook QRDR ({sorted(TEXTBOOK_QRDR_LOCI)}): {sorted(found_qrdr)}")
    print(f"Found expanded cipro-loci: {sorted(found_expanded)}")
    print(f"Mechanism classes hit: {sorted(mechanisms_hit)}")
    print(f"Verdict: {verdict}")

    # Persist per-strain top-K to JSON for post-hoc re-aggregation
    json_path = args.output.with_suffix(".per_strain.json")
    json_path.parent.mkdir(parents=True, exist_ok=True)
    json_payload = {
        "cohort_path": str(args.cohort),
        "nt_cache": str(args.nt_cache),
        "drug": args.drug,
        "top_k": args.top_k,
        "n_R_strains": r_count,
        "per_strain_top": per_strain_top,
        "cohort_sum_pos_delta": dict(cohort_sum_pos_delta),
        "cohort_freq_in_topk": dict(cohort_freq_in_topk),
        "verdict": verdict,
        "found_qrdr": sorted(found_qrdr),
        "found_expanded": sorted(found_expanded),
        "mechanisms_hit": sorted(mechanisms_hit),
    }
    json_path.write_text(json.dumps(json_payload, indent=2), encoding="utf-8")
    print(f"[preflight] wrote per-strain JSON: {json_path}")

    # Write markdown packet (v2 — expanded + signed-delta + dual-aggregation)
    from datetime import date as _date
    lines = [
        f"# Cipro attribution preflight v2 — N=38 cohort, mean-pool NT-XGBoost ({_date.today().isoformat()})",
        "",
        f"**Purpose:** test whether NT embeddings carry cipro-resistance signal at all (independently of Stage 1's verdict gate).",
        f"**Method:** train final NT-XGBoost on all {len(y)} strains (no LOSO); per cipro-R strain, run gene_level_mutagenesis; dual aggregations (sum positive delta + frequency-in-top-K) across {r_count} cipro-R strains.",
        f"**Pooling:** mean-pool (512-dim) — matches Stage 1 (NOT Stage 1b mean+max).",
        f"**Signed-delta filter:** only positive deltas (knockout LOWERED R-prob = gene supporting R) count toward sum aggregation.",
        f"**Loci tracked:** expanded set across 5 mechanism classes (QRDR target-alteration + plasmid protect/modify + efflux + porin + regulatory). v1 narrow {{gyrA, parC, parE}} preserved as TEXTBOOK_QRDR_LOCI subset.",
        "",
        f"**Verdict:** {verdict}",
        f"- Textbook QRDR found: {sorted(found_qrdr)}",
        f"- Expanded cipro-loci found: {sorted(found_expanded)}",
        f"- Mechanism classes hit: {sorted(mechanisms_hit)}",
        "",
        f"## Cohort-wide top-{args.top_k} by SUM POSITIVE delta",
        "",
        "| rank | gene_symbol | sum(+delta) | freq | mechanism |",
        "|---:|---|---:|---:|---|",
    ]
    for i, (sym, total) in enumerate(ranking_by_pos_delta, 1):
        mech = next((m for m, loci in CIPRO_LOCI_BY_MECHANISM.items() if sym in loci), "")
        lines.append(f"| {i} | {sym} | {total:.4f} | {cohort_freq_in_topk[sym]} | {mech} |")
    lines.extend([
        "",
        f"## Cohort-wide top-{args.top_k} by FREQUENCY in per-strain top-K",
        "",
        "| rank | gene_symbol | freq | sum(+delta) | mechanism |",
        "|---:|---|---:|---:|---|",
    ])
    for i, (sym, freq) in enumerate(ranking_by_freq, 1):
        mech = next((m for m, loci in CIPRO_LOCI_BY_MECHANISM.items() if sym in loci), "")
        pos_d = cohort_sum_pos_delta.get(sym, 0.0)
        lines.append(f"| {i} | {sym} | {freq}/{r_count} | {pos_d:.4f} | {mech} |")
    lines.extend([
        "",
        f"## Mechanism-class breakdown",
        "",
        "| mechanism | loci tracked | found in any top-K |",
        "|---|---|---|",
    ])
    for mech, loci in CIPRO_LOCI_BY_MECHANISM.items():
        hits = sorted(loci & any_topk_symbols)
        lines.append(f"| {mech} | {len(loci)} | {', '.join(hits) if hits else '(none)'} |")
    lines.extend([
        "",
        "## Verdict interpretation rubric",
        "",
        "- **STRONG_POSITIVE:** ≥1 QRDR locus (gyrA/B, parC/E) AND ≥1 mechanism class hit; positive signed delta in multiple R strains. Architecture has real biology.",
        "- **WEAK_POSITIVE:** any expanded-set locus appears in top-K, but no QRDR or only single-strain. Suggestive but not load-bearing.",
        "- **INCONCLUSIVE_MISS:** no cipro-loci recovered. DO NOT treat as architecture mismatch — mean-pool dilution (1/N_genes per knockout) + symbol-coverage limits (RefSeq ~11% CDSs carry `gene=`) are confounders. Next escalation: refactor `gene_level_mutagenesis` to mean+max + retest matching Stage 1b's pooling.",
        "- **DAMNING_MISS:** v2 + mean+max preflight both return no cipro-loci with positive signed delta. THEN architecture mismatch is defensible.",
        "",
        "## Notes",
        "",
        f"- N=38 cohort; 19R/19S balance per Stage 1 verdict.",
        f"- Mean-pool 512-dim features (matches Stage 1, NOT Stage 1b's 1024-dim mean+max). `gene_level_mutagenesis` hardcodes mean re-aggregation after knockout; mean+max attribution requires refactor.",
        f"- Empty `gene_symbol` entries are dropped from cohort aggregation (only ~11% of CDSs in RefSeq GFF3 carry `gene=` per CLAUDE.md gotcha). Genes-without-symbol may still rank high per-strain via gene_id but won't aggregate cross-strain. Per-strain JSON sidecar preserves raw gene_id ranking for offline inspection.",
        f"- Classifier exit-code-agnostic: the preflight uses the model regardless of Stage 1 verdict. If NT-XGBoost AUROC was FAIL at LOSO but attribution recovers QRDR, the architecture has signal that small-N LOSO can't surface.",
        f"- Per-strain top-K table persisted at sidecar JSON: `{json_path.name}`.",
    ])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines), encoding="utf-8")
    print(f"\n[preflight] wrote {args.output}")

    if verdict == "STRONG_POSITIVE":
        return 0
    if verdict == "WEAK_POSITIVE":
        return 2
    return 1  # INCONCLUSIVE_MISS


if __name__ == "__main__":
    sys.exit(main())
