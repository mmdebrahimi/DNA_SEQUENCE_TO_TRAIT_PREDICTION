"""Shared ##FASTA-safe GFF loader (Step 1, brainstorm catch C2).

Bakta appends an embedded ``##FASTA`` block (the genome sequence) to the end
of its GFF3. The repo parser ``dna_decode.data.annotations.parse_gff3`` (and
``load_annotation_table``) split each line on TAB and REQUIRE exactly 9 fields
— the FASTA header + sequence lines have the wrong shape and raise
``AnnotationParseError`` (or silently corrupt the table). ``scripts/
pathotype_laptop_pipeline.parse_bakta_gff3`` already strips ``##FASTA`` before
parsing; this module is the SHARED, package-level form of that fix so BOTH the
Bakta-run path AND the offline provided-GFF path go through one loader.

v1 accepts Bakta-compatible GFF only (brainstorm catch M2) — arbitrary
GFF/GBK normalization is out of scope.
"""
from __future__ import annotations

from pathlib import Path

from dna_decode.data.annotations import AnnotationTable, parse_gff3

# The Bakta GFF embedded-sequence delimiter. Everything from this line onward is
# the FASTA dump, not feature rows.
FASTA_DELIMITER = "##FASTA"

# Suffix used for the stripped temp file written next to the source GFF.
_STRIPPED_SUFFIX = ".nofasta.gff3"


def strip_fasta_block(gff_text: str) -> str:
    """Return the GFF text up to (and excluding) the embedded ``##FASTA`` block.

    Idempotent: a GFF with no ``##FASTA`` block is returned unchanged. Splits on
    the first occurrence only (the block is always last in a Bakta GFF).
    """
    return gff_text.split(FASTA_DELIMITER)[0]


def load_genome_gff(path: Path | str) -> AnnotationTable:
    """Load a (possibly Bakta-embedded-FASTA) GFF3 into the stable AnnotationTable.

    Strips the ``##FASTA`` block first (cf. ``parse_bakta_gff3``), writes a
    sibling ``<name>.nofasta.gff3`` temp file, then delegates to the repo's
    ``parse_gff3`` (which carries the two-pass parent->CDS gene_symbol
    propagation Bakta needs). A plain GFF (no ``##FASTA``) is parsed unchanged.

    This is the SINGLE entry point for both the Bakta-run path and the offline
    provided-GFF path — never call ``parse_gff3`` / ``load_annotation_table``
    on a raw Bakta GFF directly.
    """
    path = Path(path)
    raw = path.read_text(encoding="utf-8")
    if FASTA_DELIMITER not in raw:
        # Plain GFF — no temp file needed; parse in place.
        return parse_gff3(path)
    stripped = strip_fasta_block(raw)
    tmp = path.with_suffix(_STRIPPED_SUFFIX)
    tmp.write_text(stripped, encoding="utf-8")
    return parse_gff3(tmp)
