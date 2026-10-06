"""Decoder-side essentiality label registry — every parser GATED on a pinned W0 schema record.

NAMESPACE-SEPARATE from `dna_decode/fba/essentiality_labels.py` ON PURPOSE. That module is the
MODEL-GENE-KEYED registry: its `LABEL_WALLED` entry for S. aureus ("iYS854 ids are USA300HOU_####; NTML
is USA300 JE2 (SAUSA300_####) -> crosswalk") and `MODEL_WALLED` entry for P. aeruginosa ("no GEM in
BiGG ... the gold standard DOES exist and is fetchable but is keyed to PA14 locus tags -- so the blocker
is the MODEL, not the label") are about joining to a genome-scale model's gene ids. This decoder joins to
a genome ANNOTATION and needs no model, so those walls do not apply here -- but the wall MOVES to "label
ids must join the exact annotation table", which `annotation_join` gates.

NO PARSER RUNS AGAINST AN UNVERIFIED SCHEMA. `parse_or_refuse` reads the W0 record written by
`scripts/essentiality_ladder_w0_probe.py` and RAISES `SchemaUnverified` when it is missing or
`w0_verified` is false. This exists because two prose-written parsers failed silently on this plan and
both produced a NUMBER rather than an error.

A DISCREPANCY FOUND WHILE BUILDING THIS, and it shapes the contract: the W0 probe's `key_column_chosen`
is a LEXICAL CANDIDATE, not the join authority. For BAGEL it picks index 0 (`GENE`) because that column
is named gene-ish, but the real join key is index 2 (`ENTREZ_ID`) -- `Homo_sapiens.gene_info.gz` is keyed
by GeneID, which is what `build_essentiality_report_card.human_transfer_auroc` already joins on. So each
source DECLARES its own `key_column_index`, and the gate checks that the W0 record LISTS that index among
its candidates (a consistency check) rather than blindly adopting the probe's pick.
"""
from __future__ import annotations

import json
import sys
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WIKI = ROOT / "wiki"


class SchemaUnverified(RuntimeError):
    """No pinned W0 schema record for this source, or the record says it is unverified."""


class LabelParseError(RuntimeError):
    """A parser produced zero ids. An empty label set yields a clean-looking 0.0 coverage on no data."""


@dataclass(frozen=True)
class LabelSource:
    source_id: str
    organism: str
    strain: str
    key_space: str
    technology: str
    key_column_index: int        # DECLARED here; the W0 record corroborates, it does not decide
    url: str | None
    recorded_in: str


LABEL_SOURCES: dict[str, LabelSource] = {
    "goodall_tradis": LabelSource(
        "goodall_tradis", "Escherichia coli", "K-12 MG1655", "gene_symbol",
        "transposon_insertion", 0, None,
        "scripts/essentiality_e3_learned.py::load_labels (reused by identity)"),
    "bagel_ceg": LabelSource(
        "bagel_ceg", "Homo sapiens", "n/a", "entrez_gene_id", "crispr_ko", 2, None,
        "scripts/build_essentiality_report_card.py::human_transfer_auroc joins on ENTREZ_ID"),
    "bagel_neg": LabelSource(
        "bagel_neg", "Homo sapiens", "n/a", "entrez_gene_id", "crispr_ko", 2, None,
        "scripts/build_essentiality_report_card.py::human_transfer_auroc joins on ENTREZ_ID"),
    "sgd_phenotype": LabelSource(
        "sgd_phenotype", "Saccharomyces cerevisiae", "S288C", "systematic_orf",
        "deletion_collection", 0,
        "http://sgd-archive.yeastgenome.org/curation/literature/phenotype_data.tab",
        "dna_decode/fba/essentiality_labels.py::parse_sgd_essential (reused by identity)"),
    "ntml_nebraska": LabelSource(
        "ntml_nebraska", "Staphylococcus aureus", "USA300 JE2 / FPR3757", "sausa300_locus_tag",
        "transposon_insertion", 0, None,
        "dna_decode/fba/essentiality_labels.py LABEL_WALLED -- a MODEL-join wall only; no resolvable "
        "URL is recorded anywhere in-repo, so reachability is a HYPOTHESIS"),
    "plos_gold_pa14": LabelSource(
        "plos_gold_pa14", "Pseudomonas aeruginosa", "PA14", "pa14_locus_tag",
        "transposon_insertion", 0, None,
        "dna_decode/fba/essentiality_labels.py MODEL_WALLED -- names PLOS pcbi.1013945.s011 sheets "
        "GOLD_84/GOLD_115 but gives no resolvable URL, so reachability is a HYPOTHESIS"),
}


# ---------------------------------------------------------------------------------------------------
# Reused by IDENTITY -- no copies. Asserted by test so these cannot silently fork.
# ---------------------------------------------------------------------------------------------------
def _goodall_loader():
    """The CORRECT Goodall xlsx parser: header row 1, skips Unclear and ambiguous rows.

    HONEST LIMIT: it takes NO path argument -- it hardcodes the D: cache location -- so it cannot be
    pointed at a fixture without monkeypatching its module TGT. Reused anyway, because a second
    implementation of a parser whose header row is a known trap is strictly worse than this coupling.
    """
    sys.path.insert(0, str(ROOT / "scripts"))
    from essentiality_e3_learned import load_labels
    return load_labels


def _sgd_parser():
    from dna_decode.fba.essentiality_labels import parse_sgd_essential
    return parse_sgd_essential


# ---------------------------------------------------------------------------------------------------
# New PURE parsers, written against the W0 record's recorded header index and key column -- never prose.
# ---------------------------------------------------------------------------------------------------
def parse_bagel(text: str, key_column_index: int = 2) -> set[str]:
    """BAGEL CEGv2 / NEGv1 -> set of ENTREZ ids. Header is row 0; key is ENTREZ_ID at index 2."""
    out: set[str] = set()
    for i, line in enumerate(text.splitlines()):
        if i == 0:
            continue                                      # header, W0-confirmed at row 0
        parts = line.split("\t")
        if len(parts) > key_column_index:
            v = parts[key_column_index].strip()
            if v:
                out.add(v)
    return out


def parse_ntml(text: str, key_column_index: int = 0) -> set[str]:
    """NTML / Nebraska transposon library -> set of SAUSA300 locus tags.

    SHAPE IS UNPINNED: no NTML file has ever been on this host, so the header row, the label column and
    the call vocabulary are all unknown. This parser is therefore deliberately MINIMAL -- it extracts the
    key column of every non-header row -- and `parse_or_refuse` will not let it run until a W0 record
    exists. Finish it against the record, not against this guess.
    """
    out: set[str] = set()
    for i, line in enumerate(text.splitlines()):
        if i == 0 or not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) > key_column_index:
            v = parts[key_column_index].strip()
            if v:
                out.add(v)
    return out


def parse_plos_gold(rows, key_column_index: int = 0) -> set[str]:
    """PLOS pcbi.1013945.s011 GOLD_84 / GOLD_115 -> set of PA14 locus tags.

    SHAPE IS UNPINNED for the same reason as NTML: the workbook has never been on this host, so which
    SHEET and which header row carry the gold set is unknown. Gated identically.
    """
    out: set[str] = set()
    for i, row in enumerate(rows):
        if i == 0 or row is None:
            continue
        if len(row) > key_column_index and row[key_column_index] is not None:
            v = str(row[key_column_index]).strip()
            if v:
                out.add(v)
    return out


PARSER_FOR = {
    "bagel_ceg": parse_bagel,
    "bagel_neg": parse_bagel,
    "ntml_nebraska": parse_ntml,
    "plos_gold_pa14": parse_plos_gold,
}


# ---------------------------------------------------------------------------------------------------
# HOW EACH RUNG GETS ITS TWO CLASSES -- and they genuinely differ, which a single `label_source_id`
# cannot express. Surfaced by a Step-1/Step-5 mismatch: the ladder's human rung named one source while
# the human classes come from TWO files.
#
#   one_file_two_columns : one table carries both calls (Goodall's Essential / Non-essential columns)
#   two_files            : separate essential and non-essential sets (BAGEL CEGv2 / NEGv1)
#   essential_plus_complement : the source lists only essentials; everything else in the annotation is
#                        the negative class. NOTE this makes the negative class ANNOTATION-DEFINED, so
#                        its size depends on the annotation, not on a screen -- a real asymmetry the
#                        memo must disclose rather than average over.
# ---------------------------------------------------------------------------------------------------
RUNG_LABEL_PLAN: dict[str, dict] = {
    "ecoli": {"mode": "one_file_two_columns", "essential": "goodall_tradis",
              "nonessential": "goodall_tradis"},
    "paeruginosa": {"mode": "essential_plus_complement", "essential": "plos_gold_pa14",
                    "nonessential": None},
    "saureus": {"mode": "essential_plus_complement", "essential": "ntml_nebraska",
                "nonessential": None},
    "scerevisiae": {"mode": "essential_plus_complement", "essential": "sgd_phenotype",
                    "nonessential": None},
    "human": {"mode": "two_files", "essential": "bagel_ceg", "nonessential": "bagel_neg"},
}


def sources_for_rung(rung_key: str) -> tuple[str, ...]:
    """Every label source a rung consumes. One rung can need two files."""
    plan = RUNG_LABEL_PLAN[rung_key]
    return tuple(s for s in (plan["essential"], plan["nonessential"]) if s)


# ---------------------------------------------------------------------------------------------------
# The gate
# ---------------------------------------------------------------------------------------------------
def latest_w0_record(wiki_dir: Path | None = None) -> dict | None:
    wiki_dir = wiki_dir or WIKI
    hits = sorted(wiki_dir.glob("essentiality_ladder_w0_probe_*.json"))
    if not hits:
        return None
    return json.loads(hits[-1].read_text(encoding="utf-8"))


def w0_status(source_id: str, record: dict | None) -> tuple[bool, str]:
    if record is None:
        return False, ("no W0 probe record in wiki/ -- run "
                       "scripts/essentiality_ladder_w0_probe.py before any parser")
    for s in record.get("sources", []):
        if s.get("source_id") == source_id:
            if not s.get("w0_verified"):
                return False, "W0 record says unverified: %s" % s.get("reason", "(no reason)")
            return True, s.get("reason", "verified")
    return False, "source_id %r absent from the W0 record" % source_id


def parse_or_refuse(source_id: str, payload, *, record: dict | None = None) -> set[str]:
    """Parse a label payload ONLY against a pinned, verified schema. Fail-closed on every other path.

    `payload` is text for the tsv sources and a row iterable for xlsx ones. Returns the id set.
    """
    if source_id not in LABEL_SOURCES:
        raise SchemaUnverified("unknown label source %r" % source_id)
    src = LABEL_SOURCES[source_id]
    rec = record if record is not None else latest_w0_record()
    ok, why = w0_status(source_id, rec)
    if not ok:
        raise SchemaUnverified("%s: %s" % (source_id, why))

    # CONSISTENCY, not adoption: the probe's pick is lexical, so require only that our declared index is
    # among its candidates. For BAGEL the probe picks 0 (GENE) while the join key is 2 (ENTREZ_ID).
    entry = next(s for s in rec["sources"] if s["source_id"] == source_id)
    cands = (entry.get("schema") or {}).get("key_column_candidates")
    if cands is not None and src.key_column_index not in cands:
        raise SchemaUnverified(
            "%s: declared key_column_index %d is not among the W0 record's candidates %s -- the schema "
            "and the parser disagree about where the key lives" % (source_id, src.key_column_index, cands))

    if source_id == "goodall_tradis":
        ids = set(_goodall_loader()())
    elif source_id == "sgd_phenotype":
        ids = _sgd_parser()(payload)
    else:
        ids = PARSER_FOR[source_id](payload, src.key_column_index)

    if not ids:
        raise LabelParseError(
            "%s parsed ZERO ids. An empty label set produces a clean-looking 0.0 coverage on no data, "
            "so this refuses rather than returns an empty set." % source_id)
    return ids
