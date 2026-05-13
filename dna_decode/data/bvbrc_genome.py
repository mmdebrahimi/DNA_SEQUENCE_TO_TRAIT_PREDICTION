"""BV-BRC genome-metadata adapter — loads `BVBRC_genome.csv` (Genomes-tab export)
into the dict-of-dicts shape `cohort.candidates_from_bvbrc_ast` accepts via its
`assembly_metadata` parameter.

This bypasses the wrong-contract `pilot.fetch_ncbi_assembly_quality` scaffold
(which would return only `{contig_count, n50}` — 2 keys); the CSV supplies
7+ richer fields (`assembly_accession`, `mlst`, `country`, `year` + the 2),
all of which the cohort builder reads.

The cohort path's downstream consumer is `cohort.candidates_from_bvbrc_ast(
    ast_table, assembly_metadata=...)` at `dna_decode/data/cohort.py:324`.

Plasmid/chromosome resistance-gene tuples (which `candidates_from_bvbrc_ast`
also accepts) are AMRFinder-derived; the BV-BRC genome CSV cannot supply them.
They remain empty for strains loaded by this adapter.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

DEFAULT_ORGANISM = "Escherichia coli"

GENOME_METADATA_KEYS = (
    "assembly_accession",
    "mlst",
    "country",
    "year",
    "contig_count",
    "n50",
)

# Column-name mapping from BV-BRC Genomes export header → metadata-dict keys.
# Real BV-BRC export uses Title Case + spaces; the loader normalizes headers to
# lowercase+underscore before this lookup runs (same pattern as ast_data.py).
GENOME_COLUMN_MAP = {
    "genome_id": "strain_id",  # join key with AMR file's Genome ID
    "assembly_accession": "assembly_accession",
    "mlst": "mlst",
    "isolation_country": "country",
    "country": "country",
    "collection_year": "year",
    "contigs": "contig_count",
    "contig_count": "contig_count",
    "contig_n50": "n50",
    "n50": "n50",
}

# Columns used for organism filtering (preferred order). The first column present
# in the input is used; if neither is present, no organism filter is applied.
ORGANISM_FILTER_COLUMNS = ("species", "genome_name")


class BvBrcGenomeError(Exception):
    """Raised when the BV-BRC genome CSV is missing required columns or malformed."""


def _safe_int(raw: str) -> int:
    """Parse an int from a possibly-blank/N-A BV-BRC numeric string. Returns 0 on failure."""
    if not raw:
        return 0
    s = str(raw).strip()
    if not s or s.upper() in {"N/A", "NA", "NULL", "NONE", "-"}:
        return 0
    try:
        return int(float(s))  # tolerate "200000.0" floats
    except ValueError:
        return 0


def load_bvbrc_genome_metadata(
    csv_path: Path | str,
    organism: str = DEFAULT_ORGANISM,
) -> dict[str, dict[str, object]]:
    """Load BV-BRC `BVBRC_genome.csv` (Genomes tab export) → assembly_metadata dict.

    Args:
        csv_path: path to the CSV (or TSV — separator auto-detected).
        organism: substring-matched against the `Species` column (preferred) or
            `Genome Name` column (fallback). Case-insensitive. If neither column
            is present, no organism filter is applied (caller's responsibility).

    Returns:
        Dict keyed by string `strain_id` (BV-BRC `Genome ID`, e.g., "562.5691").
        Each value is a dict with keys: `assembly_accession`, `mlst`, `country`,
        `year`, `contig_count`, `n50`. Missing-column values default to "" / 0.

        The shape matches what `cohort.candidates_from_bvbrc_ast(...,
        assembly_metadata=...)` reads at `cohort.py:376-381`.

    Raises:
        FileNotFoundError: csv_path does not exist.
        BvBrcGenomeError: required key column missing after normalization.
    """
    path = Path(csv_path)
    if not path.exists():
        raise FileNotFoundError(f"BV-BRC genome CSV not found at {path}")

    raw = pd.read_csv(path, sep=None, engine="python", dtype=str, keep_default_na=False)

    # Normalize headers to lowercase+underscore so the column_map's literal-match
    # rename works against real BV-BRC exports ("Genome ID", "Contig N50",
    # "Isolation Country"). Same pattern as ast_data.py.
    raw.columns = [str(c).strip().lower().replace(" ", "_") for c in raw.columns]

    if "genome_id" not in raw.columns:
        raise BvBrcGenomeError(
            f"BV-BRC genome CSV missing required column 'Genome ID' "
            f"(normalized 'genome_id'). Found columns: {list(raw.columns)[:20]}"
        )

    # Apply organism filter using whichever column is present first.
    for col in ORGANISM_FILTER_COLUMNS:
        if col in raw.columns:
            raw = raw[raw[col].str.contains(organism, case=False, na=False)]
            break
    # If neither column present, no filter applied — caller's responsibility.

    # Rename only columns that exist. Duplicate targets are avoided in the map
    # design (e.g., country only comes from one source per row).
    rename_dict = {k: v for k, v in GENOME_COLUMN_MAP.items() if k in raw.columns}
    renamed = raw.rename(columns=rename_dict)

    out: dict[str, dict[str, object]] = {}
    duplicate_count = 0
    for _, row in renamed.iterrows():
        sid = str(row.get("strain_id", "")).strip()
        if not sid:
            continue
        if sid in out:
            duplicate_count += 1
            # Last-write-wins; consistent with `dict` semantics.
        out[sid] = {
            "assembly_accession": str(row.get("assembly_accession", "") or ""),
            "mlst": str(row.get("mlst", "") or ""),
            "country": str(row.get("country", "") or ""),
            "year": _safe_int(row.get("year", "") or ""),
            "contig_count": _safe_int(row.get("contig_count", "") or ""),
            "n50": _safe_int(row.get("n50", "") or ""),
        }

    if duplicate_count:
        import warnings as _warnings
        _warnings.warn(
            f"BV-BRC genome CSV had {duplicate_count} duplicate Genome ID(s); "
            f"last-write-wins applied (output keys: {len(out)}).",
            RuntimeWarning,
            stacklevel=2,
        )

    return out
