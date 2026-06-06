"""BacDive carbon-source utilization phenotype loader (EP-6 substrate).

The next embedding-frontier substrate after AMR/pathotype both failed the
"beat domain knowledge on a de-confoundable, sampling-independent label" bar.
Carbon-source utilization (Biolog / API / BacDive growth-or-no-growth assays)
is the candidate because:

  - Labels are LAB MEASUREMENTS (does strain X grow on carbon source Y?), NOT a
    sampling-context category → clears the study==class confound that killed
    pathotype (see `feedback_sampling_defined_phenotype_intrinsic_confound`).
  - No single AMRFinder-style curated mechanism catalog to lose to (the niche
    where a learned decoder can plausibly earn its keep).
  - Large public paired genome+phenotype cohort (BacDive: thousands of sequenced
    strains × dozens of carbon sources; Li et al. 2023 PMC10729968 used 4397×58).

CAVEAT (the open question this loader EXISTS to test): Li et al. found carbon-
utilization prediction is largely PHYLOGENETIC — gene-content RF nails it
in-clade but fails out-of-clade until scale rescues it. That is the SAME
lineage-vs-mechanism crux we built `cohort_deconfound.py` for. So this loader
feeds the de-confound gate FIRST; a carbon source only becomes a substrate if a
within-lineage utilizer/non-utilizer cohort exists.

This module is the carbon-util analogue of `ast_data.py`: it parses a BacDive
export into a strain × carbon-source × utilization-label table, binarizes
(positive growth → 1, negative → 0), filters to E. coli, and exposes the same
`get_*_list` / `get_binary_labels` accessor shape. Pure parsing + pandas; no
network (BacDive API integration is deferred — a downloaded export is the input,
same posture as the BV-BRC AST TSV).
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

DEFAULT_ORGANISM = "Escherichia coli"

CARBON_COLUMNS = (
    "strain_id",
    "carbon_source",
    "utilization",          # binarized 0/1
    "utilization_label",    # canonical source token (positive / negative / ...)
    "organism",
    "assembly_accession",
)

# BacDive utilization values are reported with several vocabularies depending on
# the assay (Biolog +/-, API growth, literature curation). Map all to binary:
# utilizes → 1, does-not-utilize → 0. Ambiguous / borderline values are dropped
# (not forced) — same discipline as ast_data dropping unrecognized S/I/R labels.
UTILIZATION_BINARY_MAP = {
    "+": 1,
    "POSITIVE": 1,
    "POS": 1,
    "YES": 1,
    "Y": 1,
    "TRUE": 1,
    "GROWTH": 1,
    "UTILIZED": 1,
    "UTILIZES": 1,
    "1": 1,
    "-": 0,
    "NEGATIVE": 0,
    "NEG": 0,
    "NO": 0,
    "N": 0,
    "FALSE": 0,
    "NO_GROWTH": 0,
    "NOT_UTILIZED": 0,
    "0": 0,
}

# Values that explicitly mean "we don't know" → drop, never coerce to a class.
_AMBIGUOUS_TOKENS = {"", "NA", "N/A", "NULL", "NONE", "?", "W", "WEAK", "V",
                     "VARIABLE", "BORDERLINE", "ND", "NT"}

CarbonTable = pd.DataFrame


def binarize_utilization(label: str) -> int:
    """Map a BacDive utilization token to binary 0/1. Unknown → ValueError.

    Callers that want silent-drop semantics should check membership against
    UTILIZATION_BINARY_MAP first (load_bacdive_carbon does this).
    """
    normalized = (label or "").strip().upper().replace(" ", "_")
    if normalized in UTILIZATION_BINARY_MAP:
        return UTILIZATION_BINARY_MAP[normalized]
    raise ValueError(f"Unrecognized utilization label: {label!r}")


def _is_ambiguous(label: str) -> bool:
    return (label or "").strip().upper().replace(" ", "_") in _AMBIGUOUS_TOKENS


def load_bacdive_carbon(
    path: Path | str,
    organism: str = DEFAULT_ORGANISM,
) -> CarbonTable:
    """Load a BacDive carbon-utilization export → typed binarized DataFrame.

    Accepts a LONG-format export (one row per strain × carbon-source) with
    tolerant column names. Auto-detects CSV/TSV. Filters to `organism`,
    binarizes utilization, drops ambiguous/unrecognized values.

    Expected (tolerant) columns:
      strain_id   ← strain_id / bacdive_id / genome_id / id
      carbon_source ← carbon_source / substrate / metabolite / c_source
      utilization ← utilization / value / result / growth / phenotype
      organism    ← organism / species / genome_name        (optional)
      assembly_accession ← assembly_accession / ncbi_accession (optional)
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"BacDive carbon export not found at {path}")

    raw = pd.read_csv(path, sep=None, engine="python", dtype=str, keep_default_na=False)
    raw.columns = [str(c).strip().lower().replace(" ", "_") for c in raw.columns]

    column_map = {
        "bacdive_id": "strain_id",
        "genome_id": "strain_id",
        "id": "strain_id",
        "strain": "strain_id",
        "substrate": "carbon_source",
        "metabolite": "carbon_source",
        "c_source": "carbon_source",
        "value": "utilization_raw",
        "result": "utilization_raw",
        "growth": "utilization_raw",
        "phenotype": "utilization_raw",
        "utilization": "utilization_raw",
        "species": "organism",
        "genome_name": "organism",
        "organism_name": "organism",
        "ncbi_accession": "assembly_accession",
        "assembly": "assembly_accession",
    }
    renamed = raw.rename(columns={k: v for k, v in column_map.items() if k in raw.columns})

    if "organism" in renamed.columns:
        renamed = renamed[renamed["organism"].str.contains(organism, case=False, na=False)]

    out_rows: list[dict[str, object]] = []
    for _, row in renamed.iterrows():
        raw_label = (row.get("utilization_raw", "") or "").strip()
        if _is_ambiguous(raw_label):
            continue
        norm = raw_label.upper().replace(" ", "_")
        if norm not in UTILIZATION_BINARY_MAP:
            continue  # silently drop unrecognized (matches ast_data discipline)
        cs = (row.get("carbon_source", "") or "").strip().lower()
        sid = str(row.get("strain_id", "")).strip()
        if not cs or not sid:
            continue
        out_rows.append(
            {
                "strain_id": sid,
                "carbon_source": cs,
                "utilization": UTILIZATION_BINARY_MAP[norm],
                "utilization_label": norm,
                "organism": row.get("organism", ""),
                "assembly_accession": str(row.get("assembly_accession", "")).strip(),
            }
        )

    return pd.DataFrame(out_rows, columns=list(CARBON_COLUMNS))


def get_carbon_source_list(table: CarbonTable, min_strains: int = 50) -> list[str]:
    """Carbon sources with >= min_strains distinct labeled strains."""
    if len(table) == 0:
        return []
    counts = table.groupby("carbon_source")["strain_id"].nunique()
    return sorted(counts[counts >= min_strains].index.tolist())


def get_binary_labels(table: CarbonTable, carbon_source: str) -> dict[str, int]:
    """strain_id → binary utilization label for one carbon source.

    If a strain has conflicting rows for the same source, the MAJORITY value
    wins (ties → utilizes=1, the assay-positive-dominant convention); this is
    rare in BacDive but exports can carry duplicate assay records.
    """
    rows = table[table["carbon_source"] == carbon_source.lower()]
    by_strain: dict[str, list[int]] = {}
    for _, row in rows.iterrows():
        by_strain.setdefault(row["strain_id"], []).append(int(row["utilization"]))
    out: dict[str, int] = {}
    for sid, vals in by_strain.items():
        out[sid] = 1 if sum(vals) * 2 >= len(vals) else 0
    return out


@dataclass(frozen=True)
class CarbonSourceCensus:
    """Per-carbon-source feasibility row (consumed by the feasibility script)."""
    carbon_source: str
    n_strains: int
    n_positive: int
    n_negative: int
    n_with_accession: int
    minority_fraction: float   # min(pos,neg)/n — class balance proxy


def census_carbon_sources(table: CarbonTable, min_strains: int = 1) -> list[CarbonSourceCensus]:
    """Per-carbon-source counts: total / positive / negative / with-accession.

    The first feasibility filter BEFORE the de-confound gate: a carbon source
    needs both classes present and enough strains to bother running the gate.
    Sorted by n_strains descending.
    """
    out: list[CarbonSourceCensus] = []
    for cs in sorted(table["carbon_source"].unique()):
        labels = get_binary_labels(table, cs)
        if len(labels) < min_strains:
            continue
        pos = sum(1 for v in labels.values() if v == 1)
        neg = sum(1 for v in labels.values() if v == 0)
        n = len(labels)
        sub = table[table["carbon_source"] == cs]
        with_acc = sub[sub["assembly_accession"].str.len() > 0]["strain_id"].nunique()
        out.append(
            CarbonSourceCensus(
                carbon_source=cs,
                n_strains=n,
                n_positive=pos,
                n_negative=neg,
                n_with_accession=with_acc,
                minority_fraction=(min(pos, neg) / n) if n else 0.0,
            )
        )
    out.sort(key=lambda c: c.n_strains, reverse=True)
    return out
