"""Step 5 — BV-BRC AST (antimicrobial susceptibility test) phenotype loader.

Loads strain × antibiotic × susceptibility-label rows from a BV-BRC AST TSV.
Phase 1 filters to broth-microdilution (per failure-mode #4 from
post-tech-plan brainstorm).
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

DEFAULT_ORGANISM = "Escherichia coli"
DEFAULT_METHOD_FILTER = ("broth_microdilution",)

AST_COLUMNS = (
    "strain_id",
    "antibiotic",
    "susceptibility_label",
    "mic_value",
    "mic_units",
    "measurement_method",
    "source",
)

# Phase 1 binary mapping: S/I → 0 (susceptible / intermediate); R → 1 (resistant)
# Intermediate treated as susceptible for v1 binary task. MIC regression (Phase 2)
# will use the raw MIC value separately. Accepts both single-letter and word-form
# labels (BV-BRC uses 'Resistant' / 'Susceptible' / 'Intermediate').
SUSCEPTIBILITY_BINARY_MAP = {
    "S": 0,
    "SUSCEPTIBLE": 0,
    "I": 0,
    "INTERMEDIATE": 0,
    "R": 1,
    "RESISTANT": 1,
}

ASTTable = pd.DataFrame


def binarize_susceptibility(label: str) -> int:
    """Map S/I/R label to binary 0/1. Unknown labels raise ValueError."""
    normalized = (label or "").strip().upper()
    if normalized in SUSCEPTIBILITY_BINARY_MAP:
        return SUSCEPTIBILITY_BINARY_MAP[normalized]
    raise ValueError(f"Unrecognized susceptibility label: {label!r}")


def _normalize_method(raw: str) -> str:
    """Normalize measurement-method strings to a canonical form."""
    if not raw:
        return ""
    lowered = raw.strip().lower()
    if "broth" in lowered or "mic" in lowered:
        return "broth_microdilution"
    if "disk" in lowered or "diffusion" in lowered:
        return "disk_diffusion"
    if "etest" in lowered or "e-test" in lowered:
        return "etest"
    return lowered.replace(" ", "_")


def _parse_mic_value(raw: str) -> float:
    """Parse an MIC string into a numeric value. Returns NaN on failure.

    Examples: '0.5' → 0.5; '<=2' → 2; '>32' → 32 (clamped to threshold);
    '' or 'NA' → NaN.
    """
    if not raw or raw.strip().upper() in {"NA", "N/A", "NULL", "NONE", "-"}:
        return math.nan
    s = raw.strip().lstrip("<>=").strip()
    try:
        return float(s)
    except ValueError:
        return math.nan


def load_bvbrc_ast(
    tsv_path: Path | str,
    organism: str = DEFAULT_ORGANISM,
    method_filter: tuple[str, ...] = DEFAULT_METHOD_FILTER,
) -> ASTTable:
    """Load BV-BRC AST TSV → typed DataFrame.

    Filters to the requested organism + measurement methods. Binarizes
    susceptibility labels into 0/1. Parses MIC value strings into floats
    (NaN for missing / non-numeric).
    """
    path = Path(tsv_path)
    if not path.exists():
        raise FileNotFoundError(f"BV-BRC AST TSV not found at {path}")

    # Auto-detect separator: BV-BRC exports CSV by default; older docs reference TSV.
    # `sep=None` + `engine="python"` sniffs comma/tab from the first line.
    raw = pd.read_csv(path, sep=None, engine="python", dtype=str, keep_default_na=False)

    # Normalize headers to lowercase+underscore so the column_map's literal-match rename
    # works against real BV-BRC exports (which use Title Case + spaces:
    # "Genome ID", "Resistant Phenotype", "Laboratory Typing Method"). Idempotent on
    # already-normalized headers — keeps existing TSV tests passing.
    raw.columns = [str(c).strip().lower().replace(" ", "_") for c in raw.columns]

    # Tolerant column mapping — BV-BRC field names vary by release.
    # Note: `measurement` and `measurement_value` MUST NOT both map to `mic_raw`
    # in the same rename (duplicate target → pd.rename produces collision and
    # `row.get(target)` returns a Series). Real BV-BRC exports have BOTH columns;
    # prefer `measurement_value` (numeric-only) when present, fall back to
    # `measurement` (operator-prefixed, e.g., "<=0.12") otherwise.
    column_map = {
        "genome_name": "organism",
        "organism_name": "organism",
        "strain": "strain_id",
        "genome_id": "strain_id",
        "antibiotic": "antibiotic",
        "drug": "antibiotic",
        "resistant_phenotype": "susceptibility_label",
        "phenotype": "susceptibility_label",
        "measurement_unit": "mic_units",
        "laboratory_typing_method": "measurement_method",
        "testing_standard": "source",
    }
    # Pick the MIC source column: prefer `measurement_value` (numeric-only),
    # fall back to `measurement` (operator-prefixed). Only one is renamed to `mic_raw`.
    if "measurement_value" in raw.columns:
        column_map["measurement_value"] = "mic_raw"
    elif "measurement" in raw.columns:
        column_map["measurement"] = "mic_raw"

    renamed = raw.rename(columns={k: v for k, v in column_map.items() if k in raw.columns})

    # Filter by organism if column present
    if "organism" in renamed.columns:
        renamed = renamed[renamed["organism"].str.contains(organism, case=False, na=False)]

    # Build the canonical AST DataFrame
    out_rows: list[dict[str, object]] = []
    for _, row in renamed.iterrows():
        method = _normalize_method(row.get("measurement_method", "") or "")
        if method_filter and method not in method_filter:
            continue
        label = (row.get("susceptibility_label", "") or "").strip().upper()
        if label not in SUSCEPTIBILITY_BINARY_MAP:
            continue  # silently drop unrecognized labels
        out_rows.append(
            {
                "strain_id": row.get("strain_id", ""),
                "antibiotic": (row.get("antibiotic", "") or "").strip().lower(),
                "susceptibility_label": label,
                "mic_value": _parse_mic_value(row.get("mic_raw", "") or ""),
                "mic_units": row.get("mic_units", ""),
                "measurement_method": method,
                "source": row.get("source", ""),
            }
        )

    return pd.DataFrame(out_rows, columns=list(AST_COLUMNS))


def get_drug_list(ast: ASTTable, min_strains: int = 50) -> list[str]:
    """Return antibiotics with >= min_strains labeled rows."""
    if len(ast) == 0:
        return []
    counts = ast.groupby("antibiotic").size()
    return sorted(counts[counts >= min_strains].index.tolist())


def get_binary_labels(ast: ASTTable, drug: str) -> dict[str, int]:
    """Return strain_id → binary R/S label for a single drug."""
    drug_rows = ast[ast["antibiotic"] == drug.lower()]
    return {
        row["strain_id"]: SUSCEPTIBILITY_BINARY_MAP[row["susceptibility_label"]]
        for _, row in drug_rows.iterrows()
    }
