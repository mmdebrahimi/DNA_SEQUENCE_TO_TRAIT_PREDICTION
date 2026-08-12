"""Per-organism experimental essential-gene gold standards (free, fetchable, model-gene-keyed).

The FBA engine generalizes across organisms, but a *validation number* needs a per-organism essential-gene
gold standard whose keys join to that GEM's gene ids. This module is the registry of those sources + their
pure parsers.

VALIDATED (label fetchable + keys join to the BiGG model genes):
- **escherichia_coli** (iML1515): Keio RB-TnSeq mutant fitness (Bernstein 2023); see `keio.py`. STRONG.
- **saccharomyces_cerevisiae / yeast** (iMM904): SGD `phenotype_data.tab` inviable-null ORFs (Giaever/SGD).
  Keyed by systematic ORF name (YXX###W) == iMM904 gene ids. WEAK discrimination (see the validation artifact).

LABEL-WALLED (engine runs, but no clean fetchable+keyed gold standard yet — EXTERNAL wall, documented):
- **staphylococcus_aureus** (iYS854, *S. aureus* USA300_TCH1516): gene ids are `USA300HOU_####`. The
  NTML / Nebraska transposon library (the natural gold standard) is USA300 **JE2/FPR3757**
  (`SAUSA300_####`) -> near-1:1 orthologs but still a CROSSWALK, not a free join.
- **salmonella** (iYS1720, Salmonella pan-reactome): gene ids are `STM####` (S. Typhimurium LT2).

NO MODEL (refused, not substituted):
- **pseudomonas_aeruginosa**: BiGG has NO P. aeruginosa reconstruction (checked 2026-08-07), so the
  alias RAISES rather than silently loading another species' model. The gold standard DOES exist and
  is fetchable (PLOS Comput Biol 2026 `pcbi.1013945.s011`, sheets GOLD_84 / GOLD_115) but is keyed to
  **PA14** locus tags -- so the blocker is the MODEL, not the label.

CORRECTED 2026-08-07 (see `wiki/fba_wrong_organism_model_bug_2026-08-07.md`): v0.11.0-v0.12.0 mapped
`saureus` -> iYS1720 (a *Salmonella* pan-reactome) and `paeruginosa` -> iJN1463 (*P. putida*). The
"needs a crosswalk / needs a PAO1 Tn-seq set" framing here was a MISDIAGNOSIS -- no label could ever
have joined, because the model was the wrong organism.
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

# The GROWTH CONDITION each gold standard was measured in. Load-bearing: essentiality is medium-dependent,
# so scoring a minimal-medium model against rich-medium labels reports biology as model error. SGD's
# inviable-null set comes from the deletion collection on YPD (RICH); measured effect of honouring this on
# yeast/iMM904 is MCC 0.2524 -> 0.3773 with false positives falling 67 -> 13. See `fba/medium.py`.
ESSENTIALITY_LABEL_CONDITION: dict[str, str] = {
    "saccharomyces_cerevisiae": "rich",
    "yeast": "rich",
    "scerevisiae": "rich",
}

# organisms the engine runs on but that have no clean fetchable+keyed gold standard yet (honest walls)
LABEL_WALLED = {
    "staphylococcus_aureus": "iYS854 ids are USA300HOU_####; NTML is USA300 JE2 (SAUSA300_####) -> crosswalk",
    "saureus": "iYS854 ids are USA300HOU_####; NTML is USA300 JE2 (SAUSA300_####) -> crosswalk",
    "salmonella": "iYS1720 ids are STM#### (S. Typhimurium LT2); no fetchable keyed gold standard wired yet",
    "pputida": "iJN1463 (P. putida KT2440); no fetchable keyed essentiality gold standard wired yet",
    "pseudomonas_putida": "iJN1463 (P. putida KT2440); no fetchable keyed essentiality gold standard wired yet",
}

# Organisms with NO genome-scale model at all -- distinct from LABEL_WALLED (which has a model but no
# joinable label). The blocker here is the MODEL. See `model._NO_BIGG_MODEL`.
MODEL_WALLED = {
    "pseudomonas_aeruginosa": "no P. aeruginosa GEM in BiGG; gold standard exists (PLOS pcbi.1013945 "
                              "GOLD_84/GOLD_115) but is PA14-keyed and has nothing to join to",
    "paeruginosa": "no P. aeruginosa GEM in BiGG; gold standard exists but has no model to join to",
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
