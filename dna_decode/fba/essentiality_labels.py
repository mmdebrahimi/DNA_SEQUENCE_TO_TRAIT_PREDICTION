"""Per-organism experimental essential-gene gold standards (free, fetchable, model-gene-keyed).

The FBA engine generalizes across organisms, but a *validation number* needs a per-organism essential-gene
gold standard whose keys join to that GEM's gene ids. This module is the registry of those sources + their
pure parsers.

VALIDATED (label fetchable + keys join to the BiGG model genes):
- **escherichia_coli** (iML1515): Keio RB-TnSeq mutant fitness (Bernstein 2023); see `keio.py`. STRONG.
- **saccharomyces_cerevisiae / yeast** (iMM904): SGD `phenotype_data.tab` inviable-null ORFs (Giaever/SGD).
  Keyed by systematic ORF name (YXX###W) == iMM904 gene ids. WEAK discrimination (see the validation artifact).

LABEL-WALLED (engine runs, but no clean fetchable+keyed gold standard yet — EXTERNAL wall, documented):
- **staphylococcus_aureus** (iYS1720): the model's gene *ids* are STM#### (Salmonella-style locus tags) while
  its gene *names* are real symbols (ArgD, ...) -> a gold standard would need a name/symbol CROSSWALK.
- **pseudomonas_aeruginosa** (iJN1463): needs a fetchable PAO1 Tn-seq essential set keyed by PA-number.
"""
from __future__ import annotations

import re

_ORF = re.compile(r"^Y[A-P][LR]\d{3}[WC]$")  # S. cerevisiae systematic ORF name

# organism alias -> (kind, url). kind selects the parser below.
ESSENTIALITY_LABEL_SOURCES: dict[str, tuple[str, str]] = {
    "saccharomyces_cerevisiae": ("sgd", "http://sgd-archive.yeastgenome.org/curation/literature/phenotype_data.tab"),
    "yeast": ("sgd", "http://sgd-archive.yeastgenome.org/curation/literature/phenotype_data.tab"),
    "scerevisiae": ("sgd", "http://sgd-archive.yeastgenome.org/curation/literature/phenotype_data.tab"),
}

# organisms the engine runs on but that have no clean fetchable+keyed gold standard yet (honest walls)
LABEL_WALLED = {
    "staphylococcus_aureus": "iYS1720 gene ids are STM#### (Salmonella-style); needs a gene-name crosswalk",
    "saureus": "iYS1720 gene ids are STM#### (Salmonella-style); needs a gene-name crosswalk",
    "pseudomonas_aeruginosa": "needs a fetchable PAO1 Tn-seq essential set keyed by PA-number",
    "paeruginosa": "needs a fetchable PAO1 Tn-seq essential set keyed by PA-number",
}


def parse_sgd_essential(text: str) -> set[str]:
    """PURE: SGD phenotype_data.tab -> the set of essential systematic-ORF names.

    Essential = a NULL mutant with an 'inviable' phenotype. Columns (SGD phenotype_data.tab, tab-sep):
    [0]=systematic ORF, [6]=mutant_type, [9]=phenotype. We accept a row when col0 is an ORF name, the
    phenotype column contains 'inviable', and the row's mutant type is 'null' (robust to column drift by
    also checking the whole row for the 'null' token).
    """
    ess: set[str] = set()
    for ln in text.splitlines():
        c = ln.split("\t")
        if len(c) <= 9:
            continue
        orf = c[0].strip()
        if not _ORF.match(orf):
            continue
        phenotype = c[9].lower()
        if "inviable" not in phenotype:
            continue
        mutant = c[6].lower() if len(c) > 6 else ""
        if "null" in mutant or "\tnull\t" in ln.lower():
            ess.add(orf)
    return ess


def parse_essential(kind: str, text: str) -> set[str]:
    """Dispatch to the right parser by source kind."""
    if kind == "sgd":
        return parse_sgd_essential(text)
    raise ValueError(f"unknown essentiality-label kind: {kind}")
