"""Step 3 — Genome annotation parser (GFF3 + GenBank).

Produces a typed AnnotationTable (pandas DataFrame with stable schema) plus
CDS + intergenic sequence extractors. No biopython for GFF3 parsing (avoid the
bcbio-gff sidecar dependency); biopython only for GenBank + FASTA.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

import pandas as pd

# Stable column schema for AnnotationTable. Treat as a contract for downstream
# consumers (cohort builder, cohort filter, classical baselines).
ANNOTATION_COLUMNS = (
    "seqid",
    "source",
    "type",
    "start",
    "end",
    "strand",
    "gene_id",
    "locus_tag",
    "product",
)


class AnnotationParseError(Exception):
    """Malformed annotation file."""


# AnnotationTable is just an aliased pandas DataFrame for type hints; we don't
# subclass because pandas subclassing is fragile and downstream consumers
# benefit from native pandas API.
AnnotationTable = pd.DataFrame


def _parse_gff3_attrs(attr_string: str) -> dict[str, str]:
    """Parse a GFF3 9th-column attribute string (key=value; key=value)."""
    out: dict[str, str] = {}
    for kv in attr_string.split(";"):
        kv = kv.strip()
        if not kv or "=" not in kv:
            continue
        key, _, val = kv.partition("=")
        out[key.strip()] = val.strip()
    return out


def _gff3_lines(path: Path) -> Iterator[str]:
    """Yield non-comment, non-empty lines from a GFF3 file."""
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            yield line


def parse_gff3(path: Path | str) -> AnnotationTable:
    """Parse a GFF3 file into the stable AnnotationTable schema.

    Raises AnnotationParseError on malformed lines (with line-number context).
    """
    rows: list[dict[str, object]] = []
    path = Path(path)
    for lineno, line in enumerate(_gff3_lines(path), start=1):
        fields = line.split("\t")
        if len(fields) != 9:
            raise AnnotationParseError(
                f"GFF3 line {lineno} in {path.name} has {len(fields)} fields (expected 9): {line[:80]}"
            )
        seqid, source, type_, start, end, _score, strand, _phase, attrs = fields
        try:
            start_i = int(start)
            end_i = int(end)
        except ValueError as e:
            raise AnnotationParseError(
                f"GFF3 line {lineno}: non-integer start/end ({start}/{end})"
            ) from e
        attr_map = _parse_gff3_attrs(attrs)
        rows.append(
            {
                "seqid": seqid,
                "source": source,
                "type": type_,
                "start": start_i,
                "end": end_i,
                "strand": strand,
                "gene_id": attr_map.get("ID", attr_map.get("Name", "")),
                "locus_tag": attr_map.get("locus_tag", ""),
                "product": attr_map.get("product", ""),
            }
        )

    return pd.DataFrame(rows, columns=list(ANNOTATION_COLUMNS))


def parse_genbank(path: Path | str) -> AnnotationTable:
    """Parse a GenBank flat file into the stable AnnotationTable schema.

    Uses Biopython's SeqIO. Each feature with a non-empty location becomes a row.
    """
    from Bio import SeqIO  # local import — keep top-level import light

    rows: list[dict[str, object]] = []
    for record in SeqIO.parse(path, "genbank"):
        seqid = record.id
        for feat in record.features:
            quals = feat.qualifiers
            start = int(feat.location.start) + 1  # GenBank locations are 0-based half-open; emit 1-based inclusive
            end = int(feat.location.end)
            strand = "+" if feat.location.strand == 1 else "-" if feat.location.strand == -1 else "."
            rows.append(
                {
                    "seqid": seqid,
                    "source": "GenBank",
                    "type": feat.type,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "gene_id": quals.get("gene", [""])[0],
                    "locus_tag": quals.get("locus_tag", [""])[0],
                    "product": quals.get("product", [""])[0],
                }
            )

    return pd.DataFrame(rows, columns=list(ANNOTATION_COLUMNS))


def _revcomp(seq: str) -> str:
    """Reverse complement a DNA sequence (uppercase or mixed)."""
    table = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(table)[::-1]


def _load_genome_dict(fasta_path: Path | str) -> dict[str, str]:
    """Load a FASTA into a dict of seqid -> sequence string (uppercase)."""
    from Bio import SeqIO

    out: dict[str, str] = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        out[record.id] = str(record.seq).upper()
    return out


def extract_cds_sequences(
    genome_fasta: Path | str, annotations: AnnotationTable
) -> dict[str, str]:
    """Return gene_id (or locus_tag fallback) -> CDS nucleotide-sequence mapping.

    Reverse-complements `-` strand features. Only rows where type == 'CDS'.
    """
    seqs = _load_genome_dict(genome_fasta)
    out: dict[str, str] = {}

    cds_rows = annotations[annotations["type"] == "CDS"]
    for _, row in cds_rows.iterrows():
        seqid = row["seqid"]
        if seqid not in seqs:
            continue
        # GFF3 is 1-based inclusive on start, inclusive on end. Python slice is 0-based half-open.
        seq = seqs[seqid][int(row["start"]) - 1 : int(row["end"])]
        if row["strand"] == "-":
            seq = _revcomp(seq)
        key = row["gene_id"] or row["locus_tag"] or f"{seqid}:{row['start']}-{row['end']}"
        out[key] = seq

    return out


def extract_intergenic_regions(
    genome_fasta: Path | str,
    annotations: AnnotationTable,
    min_length: int = 30,
) -> dict[str, str]:
    """Return gene-pair-id -> intergenic-region nucleotide-sequence mapping.

    For each adjacent CDS pair on the same seqid (sorted by start), emits the
    sequence between gene-end and next-gene-start. Skips regions shorter than
    `min_length`.
    """
    seqs = _load_genome_dict(genome_fasta)
    out: dict[str, str] = {}

    cds = annotations[annotations["type"] == "CDS"].copy()
    cds.sort_values(["seqid", "start"], inplace=True)

    for seqid, group in cds.groupby("seqid", sort=False):
        if seqid not in seqs:
            continue
        rows = group.reset_index(drop=True)
        for i in range(len(rows) - 1):
            curr_end = int(rows.iloc[i]["end"])
            next_start = int(rows.iloc[i + 1]["start"])
            gap_len = next_start - curr_end - 1
            if gap_len < min_length:
                continue
            gap_seq = seqs[seqid][curr_end:next_start - 1]
            id_curr = rows.iloc[i]["gene_id"] or rows.iloc[i]["locus_tag"] or f"{seqid}:{curr_end}"
            id_next = (
                rows.iloc[i + 1]["gene_id"]
                or rows.iloc[i + 1]["locus_tag"]
                or f"{seqid}:{next_start}"
            )
            out[f"{id_curr}__{id_next}"] = gap_seq

    return out
