"""Step 2 — tool-surface manifest (feasibility smoke + field/vocab/header inventory).

Runs Bakta (db-light) + AMRFinderPlus on ONE genome (via the Step-1 package) and
INVENTORIES the REAL tool surfaces so the downstream tier classifier (Step 3) and
determinant overlay (Step 4) are seeded from ground truth, never guessed
(brainstorm catch C5):

  - Bakta: which AnnotationTable fields are populated, the product VOCABULARY
    (db-light wording — drives the tier-boundary patterns), feature-type counts,
    and how many features are `hypothetical protein`/empty.
  - AMRFinder: the `main.tsv` HEADERS, and crucially whether protein-id / contig /
    coordinate columns are present — that determines whether Step-4 high-confidence
    joins are even possible (if absent, joins fall to symbol-fallback -> likely an
    honest NO-GO).

A tool wedge (the documented Docker-mount-corruption history) emits a BLOCKED
status, never a fabricated manifest. Output: wiki/genome_map_tool_surface_<date>.json.

This is the feasibility GATE: it must run before Steps 3 + 4 consume it.

Usage:
    uv run python -m scripts.genome_map_tool_surface \
        --fasta D:/dna_decode_cache/refseq/GCA_xxx/genome.fna \
        --accession GCA_xxx --organism Escherichia
    # or resolve the FASTA from an accession via NCBI Datasets:
    uv run python -m scripts.genome_map_tool_surface \
        --accession GCA_002180195.1 --organism Escherichia \
        --refseq-cache D:/dna_decode_cache/refseq

Exit: 0 = manifest written (status OK), 1 = BLOCKED (written, honest), 2 = usage/IO error.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date as _date
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.genome_map import ingest

# Low-confidence product wording (homology-only-hypothesis tier signal). Mirrored
# in tier_vocab (Step 3); surfaced here so the manifest can report how many of the
# real db-light products would land in each tier under these patterns.
_LOW_CONFIDENCE_TOKENS = (
    "putative", "probable", "by similarity", "domain-containing",
    "uncharacterized", "duf", "predicted",
)
_UNKNOWN_TOKENS = ("hypothetical protein",)

# AMRFinder column-name candidates (versions differ): the manifest reports which
# are actually present so Step 4 knows the real join keys.
_PROTEIN_ID_COLS = ("Protein identifier", "Protein id")
_CONTIG_COLS = ("Contig id", "Contig")
_START_COLS = ("Start",)
_STOP_COLS = ("Stop", "End")
_SYMBOL_COLS = ("Element symbol", "Gene symbol")


def _present(headers: list[str], candidates: tuple[str, ...]) -> str | None:
    """Return the first candidate column name present in headers, else None."""
    hset = set(headers)
    for c in candidates:
        if c in hset:
            return c
    return None


def build_bakta_inventory(table: pd.DataFrame, vocab_sample_n: int = 40) -> dict:
    """Inventory a parsed Bakta AnnotationTable (the output of ingest.load_genome_gff).

    Pure — no Docker. Reports populated fields, feature-type counts, a sample of
    the distinct product vocabulary, and per-tier-signal counts (how many products
    look low-confidence vs hypothetical under the v1 patterns).
    """
    n = len(table)
    populated = [c for c in table.columns if (table[c].astype(str).str.len() > 0).any()]
    type_counts = (
        table["type"].value_counts().to_dict() if "type" in table.columns else {}
    )
    products = (
        table["product"].astype(str).str.strip() if "product" in table.columns else pd.Series([], dtype=str)
    )
    distinct = sorted({p for p in products if p})
    low_conf = [p for p in distinct if any(t in p.lower() for t in _LOW_CONFIDENCE_TOKENS)]
    hypothetical = int(sum(1 for p in products if p.lower() in _UNKNOWN_TOKENS or p == ""))
    return {
        "n_features": int(n),
        "feature_type_counts": {str(k): int(v) for k, v in type_counts.items()},
        "populated_fields": populated,
        "n_distinct_products": len(distinct),
        "product_vocabulary_sample": distinct[:vocab_sample_n],
        "low_confidence_product_examples": low_conf[:vocab_sample_n],
        "hypothetical_or_empty_count": hypothetical,
    }


def build_amrfinder_inventory(main_tsv: Path | str) -> dict:
    """Inventory an AMRFinder main.tsv: headers + coordinate/protein-id presence.

    Pure — reads the TSV header + counts rows. The has_coords / has_protein_id /
    symbol_column fields are the load-bearing signal for Step 4's join feasibility.
    A missing/empty file reports has_* = False with n_rows = 0 (not a crash).
    """
    p = Path(main_tsv)
    if not p.exists() or p.stat().st_size == 0:
        return {
            "present": False, "headers": [], "n_rows": 0,
            "has_protein_id": False, "has_contig": False, "has_coords": False,
            "protein_id_column": None, "contig_column": None,
            "start_column": None, "stop_column": None, "symbol_column": None,
        }
    df = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
    headers = list(df.columns)
    protein_col = _present(headers, _PROTEIN_ID_COLS)
    contig_col = _present(headers, _CONTIG_COLS)
    start_col = _present(headers, _START_COLS)
    stop_col = _present(headers, _STOP_COLS)
    symbol_col = _present(headers, _SYMBOL_COLS)
    return {
        "present": True,
        "headers": headers,
        "n_rows": int(len(df)),
        "has_protein_id": protein_col is not None,
        "has_contig": contig_col is not None,
        "has_coords": (start_col is not None and stop_col is not None),
        "protein_id_column": protein_col,
        "contig_column": contig_col,
        "start_column": start_col,
        "stop_column": stop_col,
        "symbol_column": symbol_col,
    }


def build_manifest(
    accession: str,
    organism: str | None,
    bakta_table: pd.DataFrame | None,
    amrfinder_main_tsv: Path | str | None,
    *,
    bakta_status: str = "OK",
    amrfinder_status: str = "OK",
    generated: str | None = None,
) -> dict:
    """Assemble the full tool-surface manifest dict.

    status is OK iff BOTH tools succeeded; a wedged tool yields a BLOCKED status
    (BAKTA_ANNOTATION_BLOCKED / AMRFINDER_BLOCKED) and the corresponding section
    is null. The downstream steps consume this artifact and BLOCKED-gate honestly.
    """
    bakta_section = (
        build_bakta_inventory(bakta_table) if (bakta_status == "OK" and bakta_table is not None) else None
    )
    amr_section = (
        build_amrfinder_inventory(amrfinder_main_tsv)
        if (amrfinder_status == "OK" and amrfinder_main_tsv is not None)
        else None
    )
    if bakta_status != "OK":
        status = "BAKTA_ANNOTATION_BLOCKED"
    elif amrfinder_status != "OK":
        status = "AMRFINDER_BLOCKED"
    else:
        status = "OK"
    return {
        "artifact": "genome_map_tool_surface",
        "schema_version": "genome-map-tool-surface-v1",
        "generated": generated or _date.today().isoformat(),
        "genome_accession": accession,
        "amrfinder_organism": organism,
        "status": status,
        "bakta_status": bakta_status,
        "amrfinder_status": amrfinder_status,
        "bakta": bakta_section,
        "amrfinder": amr_section,
    }


def _resolve_fasta(args) -> Path:
    """Locate the genome FASTA: explicit --fasta, else download by --accession."""
    if args.fasta:
        return Path(args.fasta)
    if not args.refseq_cache:
        raise SystemExit("ERROR: provide --fasta OR (--accession AND --refseq-cache)")
    from dna_decode.data.refseq import download_genome

    cache_dir = download_genome(args.accession, args.refseq_cache)
    return Path(cache_dir) / "genome.fna"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="genome_map_tool_surface", description=__doc__)
    ap.add_argument("--accession", required=True, help="genome accession (labels the manifest)")
    ap.add_argument("--fasta", type=Path, default=None, help="genome FASTA path (else download by accession)")
    ap.add_argument("--refseq-cache", default=None, help="NCBI Datasets cache root (for download-by-accession)")
    ap.add_argument("--organism", default="Escherichia",
                    help="AMRFinder -O organism; pass 'none' for a generic (no -O) scan")
    ap.add_argument("--bakta-out", default=None, help="Bakta output dir (default: <refseq>/bakta/<acc>)")
    ap.add_argument("--amrfinder-out", default=None, help="AMRFinder output dir (default: <refseq>/amrfinder/<acc>)")
    ap.add_argument("--out", type=Path, default=None,
                    help="manifest JSON path (default: wiki/genome_map_tool_surface_<date>.json)")
    args = ap.parse_args(argv)

    from dna_decode.genome_map import annotate, amrfinder
    from tools.docker_runner import DockerRunnerError

    organism = None if str(args.organism).lower() == "none" else args.organism
    out_path = args.out or (REPO / "wiki" / f"genome_map_tool_surface_{_date.today().isoformat()}.json")

    try:
        fasta = _resolve_fasta(args)
    except Exception as e:  # noqa: BLE001 — usage/IO failure is fatal before any tool runs
        print(f"ERROR resolving FASTA: {e}", file=sys.stderr)
        return 2

    base = Path(args.fasta).parent if args.fasta else (Path(args.refseq_cache) / args.accession)
    bakta_out = Path(args.bakta_out) if args.bakta_out else base / "bakta"
    amr_out = Path(args.amrfinder_out) if args.amrfinder_out else base / "amrfinder"

    # --- Bakta ---
    bakta_table = None
    bakta_status = "OK"
    try:
        gff = annotate.run_bakta(fasta, bakta_out, prefix=args.accession)
        bakta_table = ingest.load_genome_gff(gff)
    except (DockerRunnerError, RuntimeError, OSError) as e:
        bakta_status = "BAKTA_ANNOTATION_BLOCKED"
        print(f"WARN Bakta wedged -> {bakta_status}: {e}", file=sys.stderr)

    # --- AMRFinder ---
    amr_main = None
    amr_status = "OK"
    try:
        amr_main, _ = amrfinder.run_amrfinder(fasta, amr_out, organism=organism)
    except (DockerRunnerError, RuntimeError, OSError) as e:
        amr_status = "AMRFINDER_BLOCKED"
        print(f"WARN AMRFinder wedged -> {amr_status}: {e}", file=sys.stderr)

    manifest = build_manifest(
        args.accession, organism, bakta_table, amr_main,
        bakta_status=bakta_status, amrfinder_status=amr_status,
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Wrote {out_path} (status={manifest['status']})")
    return 0 if manifest["status"] == "OK" else 1


if __name__ == "__main__":
    raise SystemExit(main())
