"""Step 4 — Resistance database loaders (CARD + NCBI AMRFinderPlus).

Unifies two heterogeneous sources into a single ResistanceCatalog. Tests
inject fixture CARD JSON + AMRFinder TSV; production code refreshes via HTTP
cadence-gated (default: 30 days).
"""
from __future__ import annotations

import csv
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class ResistanceEntry:
    """Single resistance-gene record from CARD or AMRFinder."""

    gene_symbol: str  # normalized to title-case-ish; e.g., 'gyrA', 'CTX-M-15'
    gene_family: str
    drug_class: str
    resistance_mechanism: str
    source_db: str  # 'CARD' or 'AMRFinder'
    source_id: str


@dataclass
class ResistanceCatalog:
    """Collection of ResistanceEntry rows from one or both sources.

    Methods support case-insensitive gene-symbol lookup with simple alias
    mapping (CARD's lowercase first letter vs AMRFinder's uppercase).
    """

    entries: list[ResistanceEntry] = field(default_factory=list)
    _aliases: dict[str, str] = field(default_factory=dict)  # lower(symbol) -> canonical

    def __post_init__(self) -> None:
        # Build the alias index from inserted entries
        for e in self.entries:
            self._aliases.setdefault(e.gene_symbol.lower(), e.gene_symbol)

    def __len__(self) -> int:
        return len(self.entries)

    def add(self, entry: ResistanceEntry) -> None:
        self.entries.append(entry)
        self._aliases.setdefault(entry.gene_symbol.lower(), entry.gene_symbol)

    def extend(self, entries: Iterable[ResistanceEntry]) -> None:
        for e in entries:
            self.add(e)

    def map_gene_to_resistance(self, gene_symbol: str) -> list[ResistanceEntry]:
        """Case-insensitive lookup; returns all entries matching the symbol."""
        target = gene_symbol.lower()
        return [e for e in self.entries if e.gene_symbol.lower() == target]

    def all_gene_symbols(self) -> set[str]:
        return {e.gene_symbol for e in self.entries}

    def filter_by_drug_class(self, drug_class: str) -> list[ResistanceEntry]:
        """Case-insensitive substring match on drug_class."""
        target = drug_class.lower()
        return [e for e in self.entries if target in e.drug_class.lower()]


def load_card(card_json_path: Path | str) -> ResistanceCatalog:
    """Parse the CARD `card.json` model into a ResistanceCatalog.

    CARD's JSON is a deep nested structure keyed by ARO ID; each entry has
    `model_name` (gene symbol), `ARO_category` for drug class + resistance
    mechanism, and `ARO_accession` as source_id.
    """
    path = Path(card_json_path)
    if not path.exists():
        raise FileNotFoundError(f"CARD JSON not found at {path}")
    with open(path, encoding="utf-8") as f:
        data = json.load(f)

    catalog = ResistanceCatalog()
    for aro_id, record in data.items():
        if not isinstance(record, dict):
            continue
        gene_symbol = record.get("model_name", "")
        if not gene_symbol:
            continue
        categories = record.get("ARO_category", {})
        gene_family = ""
        drug_class = "unknown"
        mechanism = ""

        for cat in categories.values() if isinstance(categories, dict) else []:
            cat_class = cat.get("category_aro_class_name", "")
            cat_name = cat.get("category_aro_name", "")
            if cat_class == "AMR Gene Family":
                gene_family = cat_name
            elif cat_class == "Drug Class":
                drug_class = cat_name
            elif cat_class == "Resistance Mechanism":
                mechanism = cat_name

        catalog.add(
            ResistanceEntry(
                gene_symbol=gene_symbol,
                gene_family=gene_family,
                drug_class=drug_class or "unknown",
                resistance_mechanism=mechanism,
                source_db="CARD",
                source_id=record.get("ARO_accession", aro_id),
            )
        )

    return catalog


def load_amrfinder(amrfinder_tsv_path: Path | str) -> ResistanceCatalog:
    """Parse the NCBI AMRFinderPlus reference TSV into a ResistanceCatalog.

    AMRFinder TSV columns include: gene_symbol, gene_family, class
    (drug class), subclass, resistance_mechanism. Schema varies slightly
    by AMRFinder release; we read tolerantly via DictReader.
    """
    path = Path(amrfinder_tsv_path)
    if not path.exists():
        raise FileNotFoundError(f"AMRFinder TSV not found at {path}")

    catalog = ResistanceCatalog()
    with open(path, encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gene_symbol = (row.get("gene_symbol") or row.get("Gene symbol") or "").strip()
            if not gene_symbol:
                continue
            catalog.add(
                ResistanceEntry(
                    gene_symbol=gene_symbol,
                    gene_family=(row.get("gene_family") or row.get("Family") or "").strip(),
                    drug_class=(row.get("class") or row.get("Class") or "unknown").strip(),
                    resistance_mechanism=(
                        row.get("resistance_mechanism")
                        or row.get("Mechanism")
                        or ""
                    ).strip(),
                    source_db="AMRFinder",
                    source_id=(row.get("accession") or row.get("refseq_protein_accession") or gene_symbol).strip(),
                )
            )

    return catalog


def merge_catalogs(*catalogs: ResistanceCatalog) -> ResistanceCatalog:
    """Concatenate multiple ResistanceCatalog into one."""
    merged = ResistanceCatalog()
    for c in catalogs:
        merged.extend(c.entries)
    return merged
