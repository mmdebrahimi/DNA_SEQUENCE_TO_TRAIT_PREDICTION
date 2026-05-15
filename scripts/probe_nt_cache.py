"""Probe an NT embedding cache for completeness + integrity.

Thin CLI wrapper over `EmbeddingCache.verify_complete`. Reports per-strain
status: complete / partial / absent / corrupt / unresolved. Exit 0 only if
all resolved strains complete.

Motivation: cache.populate() skips at gene-dataset level; Stage 1 admits
any strain with >=1 cached gene + mean-pools whatever's there. A
crash-truncated strain becomes a silent landmine. Run this before any
Stage 1 invocation on a populate that may have been interrupted.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from dna_decode.data.annotations import extract_cds_sequences, parse_gff3
from dna_decode.data.cohort import load_cohort
from dna_decode.data.refseq import fasta_path, gff_path
from dna_decode.models.cache import (
    EmbeddingCache,
    EmbeddingCacheError,
    EmbeddingCacheVersionMismatch,
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort", required=True, type=Path)
    parser.add_argument("--cache", required=True, type=Path)
    parser.add_argument("--refseq-cache", required=True, type=Path)
    parser.add_argument("--embedding-dim", type=int, default=512)
    parser.add_argument(
        "--model-name",
        default="nucleotide_transformer",
    )
    parser.add_argument(
        "--model-version",
        default="InstaDeepAI/nucleotide-transformer-v2-100m-multi-species",
    )
    args = parser.parse_args(argv)

    if not args.cache.exists():
        print(f"FAIL: cache not found at {args.cache}", file=sys.stderr)
        return 2

    size = args.cache.stat().st_size
    print(f"[probe] cache file size: {size:,} bytes")

    try:
        cache = EmbeddingCache(
            args.cache,
            model_name=args.model_name,
            model_version=args.model_version,
            embedding_dim=args.embedding_dim,
        )
    except EmbeddingCacheVersionMismatch as e:
        print(f"FAIL: cache metadata mismatch: {e}", file=sys.stderr)
        return 2
    except EmbeddingCacheError as e:
        print(f"FAIL: cache open failed: {e}", file=sys.stderr)
        return 2
    except OSError as e:
        print(
            f"FAIL: HDF5 open failed (possible EOA-vs-size corruption or lock): {e}",
            file=sys.stderr,
        )
        return 2

    md = cache.metadata()
    print(
        f"[probe] root attrs: model_name={md.model_name!r} "
        f"version={md.model_version!r} dim={md.embedding_dim} "
        f"pooling={md.pooling_strategy!r} created_at={md.created_at!r}"
    )

    cohort = load_cohort(args.cohort)
    print(f"[probe] cohort has {len(cohort)} strains")

    expected_genes_by_strain: dict[str, set[str]] = {}
    unresolved: dict[str, str] = {}
    for strain in cohort.strains:
        sid = strain.strain_id
        acc = getattr(strain, "assembly_accession", None)
        if not acc:
            unresolved[sid] = "no assembly_accession"
            continue
        fna = fasta_path(acc, args.refseq_cache)
        gff = gff_path(acc, args.refseq_cache)
        if not fna.exists() or not gff.exists():
            unresolved[sid] = "missing fna/gff"
            continue
        try:
            ann = parse_gff3(gff)
            expected_genes_by_strain[sid] = set(extract_cds_sequences(fna, ann).keys())
        except Exception as e:
            unresolved[sid] = f"parse failed: {e!r}"

    report = cache.verify_complete(expected_genes_by_strain)

    print("\n=== Per-strain status ===")
    for sid in sorted(set(expected_genes_by_strain) | set(unresolved)):
        if sid in unresolved:
            print(f"  {sid:40s} unresolved: {unresolved[sid]}")
        else:
            s = report.status[sid]
            line = (
                f"  {sid:40s} {s:9s} "
                f"expected={report.expected_n[sid]} "
                f"cached={report.cached_n[sid]} "
                f"missing={len(report.missing_genes[sid])} "
                f"corrupt={len(report.corrupt_details[sid])}"
            )
            print(line)

    print("\n=== Summary ===")
    counts = report.counts
    counts["unresolved"] = len(unresolved)
    for status in ("complete", "partial", "absent", "corrupt", "unresolved"):
        if counts.get(status, 0) > 0:
            print(f"  {status}: {counts[status]}")

    if report.all_complete:
        print(
            f"\nVERDICT: ALL_COMPLETE ({len(report.status)} resolved strains complete)"
        )
        return 0
    print(
        f"\nVERDICT: INCOMPLETE — resume populate (partial/absent) or wipe+restart (corrupt)"
    )
    return 1


if __name__ == "__main__":
    sys.exit(main())
