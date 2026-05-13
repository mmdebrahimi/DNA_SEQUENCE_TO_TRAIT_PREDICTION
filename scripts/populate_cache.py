"""Phase 2 — standalone embedding-cache populate driver.

Bridges the strain_id ↔ assembly_accession key mismatch between
`pipeline.py ingest --download-genomes` (downloads by accession) and
`EmbeddingCache.populate` (keys cache by BV-BRC strain_id).

Reads a cohort parquet, resolves each strain's `genome.fna` + parsed
GFF3 annotations from the RefSeq cache via the strain's
`assembly_accession`, constructs a foundation model (mock or real), and
calls `cache.populate(model, strain_genomes, annotations)`.

Mock mode (`--model mock --allow-mock`) bypasses HuggingFace foundation
model loading and produces deterministic-hash embeddings — useful for
plumbing-only smoke runs (Gate A / Gate B infrastructure dry-run) but
biologically meaningless. The `--allow-mock` gate prevents accidental
no-signal runs.

Deferred until contracts settle: wrapping this into `pipeline.py embed`
as a fifth subcommand.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from dna_decode.data.annotations import parse_gff3
from dna_decode.data.cohort import load_cohort
from dna_decode.data.refseq import fasta_path, gff_path
from dna_decode.models.cache import EmbeddingCache
from dna_decode.models.foundation import MockFoundationModel, model_factory


DEFAULT_CONFIG_PATH = Path("config/datasources.yaml")


def resolve_strain_assets(
    cohort,
    refseq_cache_root: Path,
) -> tuple[dict, dict, list[tuple[str, str]]]:
    """Build strain_id → fasta_path + strain_id → parsed annotation tables.

    Returns (strain_genomes, strain_annotations, skipped) where `skipped` is a
    list of (strain_id, reason) tuples for strains that could not be resolved.
    Reasons: missing assembly_accession, no genome.fna, no annotations.gff3.
    """
    strain_genomes: dict[str, Path] = {}
    strain_annotations: dict = {}
    skipped: list[tuple[str, str]] = []

    for strain in cohort.strains:
        acc = getattr(strain, "assembly_accession", None)
        if not acc:
            skipped.append((strain.strain_id, "missing assembly_accession"))
            continue
        fna = fasta_path(acc, refseq_cache_root)
        gff = gff_path(acc, refseq_cache_root)
        if not fna.exists():
            skipped.append((strain.strain_id, f"no genome.fna at {fna}"))
            continue
        if not gff.exists():
            skipped.append((strain.strain_id, f"no annotations.gff3 at {gff}"))
            continue
        strain_genomes[strain.strain_id] = fna
        strain_annotations[strain.strain_id] = parse_gff3(gff)

    return strain_genomes, strain_annotations, skipped


def build_model(name: str, config_path: Path, device: str):
    """Construct a foundation-model wrapper for `name`.

    Mock short-circuits past the config (MockFoundationModel takes no
    HuggingFace deps). Real models go through `model_factory` which reads
    metadata from `config/datasources.yaml`.
    """
    if name == "mock":
        return MockFoundationModel(device=device)
    return model_factory(name, config_path=config_path, device=device)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Populate the embedding cache for a cohort. See module docstring."
    )
    parser.add_argument("--cohort", required=True, type=Path, help="Path to cohort parquet")
    parser.add_argument(
        "--model",
        default="dnabert2",
        help="Foundation model name (mock/evo/dnabert2/nucleotide_transformer/gena_lm)",
    )
    parser.add_argument(
        "--cache", required=True, type=Path, help="Path to HDF5 embedding cache file"
    )
    parser.add_argument(
        "--refseq-cache",
        required=True,
        type=Path,
        help="Path to RefSeq genome cache root (parent of per-accession dirs)",
    )
    parser.add_argument(
        "--config",
        default=DEFAULT_CONFIG_PATH,
        type=Path,
        help="datasources.yaml (for real foundation-model metadata)",
    )
    parser.add_argument(
        "--device", default="cpu", help="cpu or cuda (default cpu)"
    )
    parser.add_argument(
        "--allow-mock",
        action="store_true",
        help="Permit --model mock. Required to avoid accidental no-biology runs.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.model == "mock" and not args.allow_mock:
        print(
            "[populate_cache] --model mock requires --allow-mock (guards against accidental no-signal runs)",
            file=sys.stderr,
        )
        return 2

    if not args.cohort.exists():
        print(f"[populate_cache] cohort parquet not found: {args.cohort}", file=sys.stderr)
        return 2

    if not args.refseq_cache.exists():
        print(
            f"[populate_cache] refseq cache root not found: {args.refseq_cache} — did you run `pipeline ingest --download-genomes`?",
            file=sys.stderr,
        )
        return 2

    cohort = load_cohort(args.cohort)
    print(f"[populate_cache] loaded cohort: {len(cohort)} strains")

    strain_genomes, strain_annotations, skipped = resolve_strain_assets(
        cohort, args.refseq_cache
    )
    print(
        f"[populate_cache] resolved {len(strain_genomes)} strain(s); skipped {len(skipped)}"
    )
    for sid, reason in skipped[:10]:
        print(f"[populate_cache]   skip {sid}: {reason}")
    if len(skipped) > 10:
        print(f"[populate_cache]   ... and {len(skipped) - 10} more")
    if not strain_genomes:
        print(
            "[populate_cache] no strains resolved — refseq cache empty or assembly_accessions missing on cohort",
            file=sys.stderr,
        )
        return 1

    try:
        model = build_model(args.model, args.config, args.device)
    except Exception as e:  # FoundationModelError or upstream import failure
        print(f"[populate_cache] model construction failed: {type(e).__name__}: {e}", file=sys.stderr)
        return 1

    print(
        f"[populate_cache] model={args.model} (embedding_dim={model.metadata.embedding_dim}, "
        f"device={args.device})"
    )

    cache = EmbeddingCache(
        args.cache,
        model_name=args.model,
        model_version=model.metadata.huggingface_id,
        embedding_dim=model.metadata.embedding_dim,
    )

    def progress(strain_id: str, n_written: int, n_total: int) -> None:
        print(f"[populate_cache]   {strain_id}: wrote {n_written} / {n_total} embeddings")

    try:
        written_per_strain = cache.populate(
            model=model,
            strain_genomes=strain_genomes,
            annotations=strain_annotations,
            skip_existing=True,
            progress_callback=progress,
        )
    except Exception as e:
        print(f"[populate_cache] populate failed: {type(e).__name__}: {e}", file=sys.stderr)
        return 1

    total = sum(written_per_strain.values())
    print(
        f"[populate_cache] DONE: {total} new embeddings across "
        f"{len(written_per_strain)} strains → {args.cache}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
